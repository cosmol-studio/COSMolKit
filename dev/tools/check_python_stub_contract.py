"""Check generated Python declarations against the compiled Rust registry.

The generator supplies BINDING_CONTRACT from its own linked cosmolkit build.
This module neither parses registry source nor maintains a second API list.
"""

from __future__ import annotations

import ast
import json
from typing import Literal, NotRequired, TypedDict, cast


class ContractEntry(TypedDict):
    semantic_id: str
    python_name: str
    feature: str
    item: Literal["callable", "type"]
    owner: Literal["module", "molecule", "type"]
    python_property: NotRequired[Literal["getter", "setter"] | None]


def python_document(document):
    """Select explicitly registered Python configuration schemas, when present."""
    if isinstance(document, list):
        return document
    document = dict(document)
    document["entries"] = [
        dict(row, fields=row["python_fields"], python_field_types=True)
        if row.get("python_fields") is not None else row
        for row in document["entries"]
    ]
    return document


def configuration_schema(stub, document):
    """Inspect compiled derive metadata for registration; never write a stub."""
    classes = {node.name: node for node in ast.parse(stub).body if isinstance(node, ast.ClassDef)}
    schemas = {}
    for row in document["entries"]:
        if row.get("role") != "parameter" or row["python_name"] not in classes:
            continue
        cls = classes[row["python_name"]]
        constructors = _functions(cls.body, "__new__") or _functions(cls.body, "__init__")
        if len(constructors) != 1:
            raise ValueError(f"{row['python_name']}: expected one actual constructor")
        schemas[row["semantic_id"]] = [
            {"name": name, "type": ast.unparse(arg.annotation),
             "default": None if default is None else ast.unparse(default)}
            for name, (arg, default, _) in _arguments(constructors[0]).items()
        ]
    return json.dumps(schemas)


def _native_projection_valid(projection, types):
    if projection in {"builtins.int", "builtins.str", "builtins.float", "builtins.bool", "builtins.bytes"}:
        return True
    if projection.startswith("list[") and projection.endswith("]"):
        name = projection[5:-1]
        return name.isidentifier() and name in types
    names = [name.strip() for name in projection.split("|")]
    return len(names) > 1 and all(name.isidentifier() and name in types for name in names)


def dynamic_type_declarations(module, stub, document):
    """Describe actual dynamic exceptions/enums without derive inventory.

    Never create a type or promise a missing export. Ordinary value classes
    still require their compiled derive metadata.
    """
    import enum
    declared = {node.name for node in ast.parse(stub).body if isinstance(node, ast.ClassDef)}
    declarations = []
    exports = []
    for row in document["entries"]:
        name = row["python_name"]
        if row["item"] != "type" or name in declared:
            continue
        cls = getattr(module, name, None)
        if not isinstance(cls, type) or cls.__name__ != name:
            continue
        if issubclass(cls, enum.Enum):
            base = "enum.IntEnum" if issubclass(cls, enum.IntEnum) else "builtins.str, enum.Enum" if issubclass(cls, str) else "enum.Enum"
            members = "".join(f"    {key} = {value.value!r}\n" for key, value in cls.__members__.items())
            declarations.append(f"\nclass {name}({base}):\n" + (members or "    ...\n"))
        elif issubclass(cls, BaseException):
            base = next(base for base in cls.__mro__[1:] if base.__module__ == "builtins")
            declarations.append(f"\nclass {name}(builtins.{base.__name__}): ...\n")
        else:
            continue
        declared.add(name)
        exports.append(name)
    if exports:
        declarations.append(f"\n__all__ += {exports!r}\n")
    return "".join(declarations)


def check_contract(stub: str, contract_json: str) -> list[str]:
    """Check declarations, configuration fields and same-name call forms.

    The complete document comes from the linked Rust registry. A list remains
    useful for focused checker regressions, not as a production API inventory.
    """
    document = python_document(json.loads(contract_json))
    entries = cast(list[ContractEntry], document if isinstance(document, list) else document["entries"])
    module = ast.parse(stub, filename="cosmolkit.pyi")
    functions: set[str] = set()
    classes: dict[str, set[str]] = {}
    properties: dict[str, dict[str, set[str]]] = {}
    for node in module.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            functions.add(node.name)
        elif isinstance(node, ast.ClassDef):
            accessors: dict[str, set[str]] = {}
            for member in node.body:
                if isinstance(member, ast.AnnAssign) and isinstance(member.target, ast.Name):
                    # stub-gen emits writable descriptors as annotated attributes.
                    annotation = ast.unparse(member.annotation)
                    if "ClassVar[" not in annotation:
                        accessors[member.target.id] = {"getter"} if "Final[" in annotation else {"getter", "setter"}
                elif isinstance(member, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    for decorator in member.decorator_list:
                        if isinstance(decorator, ast.Name) and decorator.id == "property":
                            accessors.setdefault(member.name, set()).add("getter")
                        elif isinstance(decorator, ast.Attribute) and decorator.attr == "setter" and isinstance(decorator.value, ast.Name):
                            accessors.setdefault(decorator.value.id, set()).add("setter")
            properties[node.name] = accessors
            classes[node.name] = {
                member.name
                for member in node.body
                if isinstance(member, (ast.FunctionDef, ast.AsyncFunctionDef))
                and not any(
                    (isinstance(decorator, ast.Name) and decorator.id == "property")
                    or (isinstance(decorator, ast.Attribute) and decorator.attr == "setter")
                    for decorator in member.decorator_list
                )
            }

    type_names = {
        entry["semantic_id"].rsplit(".", 1)[-1]: entry["python_name"]
        for entry in entries
        if entry["item"] == "type"
    }
    missing: list[str] = []
    for entry in entries:
        if entry["item"] != "callable":
            if entry.get("python_native") is not None:
                # Explicit canonical native projection; no ck wrapper is promised.
                if not _native_projection_valid(entry["python_native"], classes):
                    missing.append(f"{entry['semantic_id']}: invalid native value projection")
                continue
            if entry["python_name"] not in classes:
                missing.append(f"{entry['semantic_id']} -> cosmolkit.{entry['python_name']} [feature={entry['feature']}]")
            for field in entry.get("properties", []):
                if "getter" not in properties.get(entry["python_name"], {}).get(field["name"], set()):
                    missing.append(f"{entry['python_name']}.{field['name']}: registered property getter missing")
            continue
        name = entry["python_name"]
        if entry["owner"] == "module":
            path = name
            found = name in functions
        else:
            owner = (
                "Molecule"
                if entry["owner"] == "molecule"
                else entry["semantic_id"].rsplit(".", 1)[0]
            )
            owner = type_names.get(owner, owner)
            path = f"{owner}.{name}"
            access = entry.get("python_property")
            if access is not None and access not in ("getter", "setter"):
                raise ValueError(f"invalid Python property access: {access}")
            found = (
                access in properties.get(owner, {}).get(name, set())
                if access is not None
                else name in classes.get(owner, set())
            )
        if not found:
            missing.append(
                f"{entry['semantic_id']} -> cosmolkit.{path} "
                + f"[feature={entry['feature']}]"
            )
    if isinstance(document, dict):
        missing.extend(check_configuration(module, document))
        missing.extend(check_path_declarations(module, document))
        for adapter in document["python_adapters"]:
            owner = next(row for row in entries if row["semantic_id"] == adapter["type_semantic_id"])
            if adapter["name"] not in classes.get(owner["python_name"], set()):
                missing.append(f"{owner['python_name']}.{adapter['name']}: registered object adapter missing")
        missing.extend(check_alias_declarations(module, document))
        missing.extend(check_collection_declarations(module, document))
    return missing


def check_collection_declarations(tree, document):
    """Enforce registered Python collection projections on every producer."""
    classes = {node.name: node for node in tree.body if isinstance(node, ast.ClassDef)}
    entries = document["entries"]
    types = {row["semantic_id"].removeprefix("types."): row["python_name"]
             for row in entries if row["item"] == "type"}
    errors = []
    for collection in document.get("python_collections", []):
        name = collection["name"]
        cls = classes.get(name)
        if cls is None:
            errors.append(f"{name}: registered Rust-backed collection missing")
            continue
        for method in ("__new__", "__len__", "__getitem__", "__iter__", "__repr__", "to_numpy"):
            if not _functions(cls.body, method):
                errors.append(f"{name}.{method}: collection method missing")
        for method in _functions(cls.body, "to_numpy"):
            if _annotation(method.returns) != "numpy.NDArray[numpy.uint8]":
                errors.append(f"{name}.to_numpy: expected uint8 NumPy matrix annotation")
        indexed = {_annotation(method.returns) for method in _functions(cls.body, "__getitem__")}
        if indexed != {"Fingerprint|None", name}:
            errors.append(f"{name}.__getitem__: scalar/slice overloads must preserve their types")
        for row in entries:
            if row["item"] != "callable" or compact(row.get("output") or "") != compact(collection["rust_output"]):
                continue
            owner = row["semantic_id"].rsplit(".", 1)[0]
            owner = classes.get(types.get(owner, owner))
            methods = [] if owner is None else _functions(owner.body, row["python_name"])
            if not methods or any(_annotation(method.returns) != name for method in methods):
                errors.append(f"{row['semantic_id']}: expected registered {name} return type")
    return errors


def check_collection_runtime(module, document):
    errors = []
    for collection in document.get("python_collections", []):
        name = collection["name"]
        cls = getattr(module, name, None)
        if not isinstance(cls, type) or not callable(getattr(cls, "to_numpy", None)):
            errors.append(f"{name}: actual Rust-backed NumPy collection missing")
            continue
        try:
            import numpy as np
            fp = module.Fingerprint.from_on_bits(8, [1, 4])
            expected = np.array([0, 1, 0, 0, 1, 0, 0, 0], dtype=np.uint8)
            np.testing.assert_array_equal(fp.to_numpy(), expected)
            batch = cls([fp])
            array = batch.to_numpy()
            np.testing.assert_array_equal(array, expected.reshape(1, 8))
            if array.dtype != np.uint8 or not array.flags.c_contiguous or not array.flags.writeable:
                errors.append(f"{name}.to_numpy: expected independent C-contiguous uint8 storage")
            array[:] = 0
            np.testing.assert_array_equal(batch.to_numpy(), expected.reshape(1, 8))
            if batch[0] is not fp or not isinstance(batch[:], cls) or list(batch) != [fp]:
                errors.append(f"{name}: index, slice or iteration changed collection semantics")
        except Exception as error:
            errors.append(f"{name}: actual NumPy projection failed: {error}")
    return errors


def check_alias_declarations(tree, document):
    """Aliases inherit the original overloads; separate signatures are forbidden."""
    entries = document["entries"]
    types = {row["semantic_id"].removeprefix("types."): row["python_name"]
             for row in entries if row["item"] == "type"}
    registered = set()
    for row in entries:
        if row["item"] == "callable" and row["owner"] != "module" and row.get("receiver") is None and row.get("python_property") is None:
            owner = "Molecule" if row["owner"] == "molecule" else row["semantic_id"].rsplit(".", 1)[0]
            registered.add(f"{types.get(owner, owner)}.{row['python_name']}")
    for adapter in document.get("python_adapters", []):
        owner = next(row for row in entries if row["semantic_id"] == adapter["type_semantic_id"])
        registered.add(f"{owner['python_name']}.{adapter['name']}")
    classes = {node.name: node for node in tree.body if isinstance(node, ast.ClassDef)}
    errors = []
    seen = set()
    for alias in document.get("python_aliases", []):
        name, target = alias["name"], alias["target"]
        if not name.isidentifier() or name.startswith("_") or name in seen:
            errors.append(f"{name}: invalid or duplicate public alias")
        seen.add(name)
        owner, _, method = target.partition(".")
        if target not in registered or owner not in classes or not _functions(classes[owner].body, method):
            errors.append(f"{name}: alias target {target} is not a registered class callable")
        bindings = [node for node in tree.body if
                    isinstance(node, (ast.FunctionDef, ast.ClassDef)) and node.name == name
                    or isinstance(node, ast.Assign) and any(isinstance(value, ast.Name) and value.id == name for value in node.targets)
                    or isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name) and node.target.id == name]
        if len(bindings) != 1 or not isinstance(bindings[0], ast.Assign) or len(bindings[0].targets) != 1 or ast.unparse(bindings[0].value) != target:
            errors.append(f"{name}: expected direct assignment to {target}, not a separate signature or wrapper")
    return errors


def check_alias_runtime(module, document):
    errors = []
    for alias in document.get("python_aliases", []):
        owner, method = alias["target"].split(".")
        target = getattr(getattr(module, owner, None), method, None)
        actual = getattr(module, alias["name"], None)
        # Native classmethods may produce a fresh bound-method object on each
        # lookup; their equality compares the same function and bound class.
        if not callable(target) or not callable(actual) or actual != target:
            errors.append(f"{alias['name']}: actual export must be the existing {alias['target']} callable")
        if hasattr(module, "__all__") and alias["name"] not in module.__all__:
            errors.append(f"{alias['name']}: package export list must include the constructor alias")
    return errors


def sdf_path_calls(document):
    """Select SDF filesystem arguments from enabled canonical registrations."""
    for row in document["entries"]:
        if row["item"] != "callable":
            continue
        owner, _, name = row["semantic_id"].rpartition(".")
        if owner in ("SdfDataset", "SdfReader", "BatchExportReport") or owner == "MoleculeBatch" and name.startswith(("read_sdf", "write_sdf")):
            fields = [field for field in row.get("parameters", []) or [] if field["name"] in ("path", "directory", "report_path")]
            if fields:
                yield row, owner, fields


def _path_types(node):
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.BitOr):
        return _path_types(node.left) | _path_types(node.right)
    if isinstance(node, ast.Subscript) and _annotation(node.value) in ("Union", "Optional"):
        values = node.slice.elts if isinstance(node.slice, ast.Tuple) else [node.slice]
        result = set().union(*(_path_types(value) for value in values))
        return result | {"None"} if _annotation(node.value) == "Optional" else result
    return {_annotation(node)}


def check_path_declarations(stub, document):
    classes = {node.name: node for node in stub.body if isinstance(node, ast.ClassDef)}
    errors = []
    for row, owner, fields in sdf_path_calls(document):
        cls = classes.get(owner)
        if cls is None:
            continue  # Missing exports are reported by check_contract.
        for function in _functions(cls.body, row["python_name"]):
            arguments = _arguments(function)
            for name, (argument, default, _) in arguments.items():
                if name not in ("path", "directory", "out_dir", "report_path"):
                    continue
                # Rust's short method calls this `path`; the configured method
                # calls it `directory`, and Python's short form uses `out_dir`.
                registered = next((field for field in fields if field["name"] == name or name != "report_path" and field["name"] != "report_path"), None)
                optional = registered is not None and compact(registered["type"]).startswith("Option<") or isinstance(default, ast.Constant) and default.value is None
                required = {"str", "os.PathLike[str]"} | ({"None"} if optional else set())
                actual = _path_types(argument.annotation)
                # pathlib.Path is a redundant but valid PathLike[str] subtype.
                if not required <= actual or actual - required - {"pathlib.Path"}:
                    errors.append(f"{owner}.{function.name}({name}): expected str | os.PathLike[str]" + (" | None" if optional else ""))
            for field in fields:
                if field["name"] not in arguments and not (field["name"] != "report_path" and arguments.keys() & {"path", "directory", "out_dir"}):
                    errors.append(f"{owner}.{function.name}: missing registered path argument {field['name']}")
    return errors


def check_path_input(call, arguments, name, text="__ck_path_probe__"):
    """Verify fixed input forms, then probe rejection and error propagation."""
    from pathlib import Path
    class ObservedPath(Exception):
        pass
    class TextPath:
        def __fspath__(self):
            raise ObservedPath
    class BytesPath:
        def __fspath__(self):
            return b"__ck_path_probe__"
    class InvalidPath:
        def __fspath__(self):
            return 42
    class ValidPath:
        def __fspath__(self):
            return text
    errors = []
    for value in (text, Path(text), ValidPath()):
        try:
            call(**dict(arguments, **{name: value}))
        except Exception as error:
            errors.append(f"{name}: {type(value).__name__} text path failed ({type(error).__name__}: {error})")
    try:
        call(**dict(arguments, **{name: TextPath()}))
    except ObservedPath:
        pass
    except Exception as error:
        errors.append(f"{name}: PathLike protocol was not propagated ({type(error).__name__})")
    else:
        errors.append(f"{name}: PathLike protocol was not invoked")
    for value in (b"__ck_path_probe__", BytesPath(), InvalidPath(), object()):
        try:
            call(**dict(arguments, **{name: value}))
        except TypeError:
            pass
        except Exception as error:
            errors.append(f"{name}: invalid path must raise TypeError, got {type(error).__name__}")
        else:
            errors.append(f"{name}: invalid path was accepted")
    return errors


def check_path_runtime(module, document):
    import inspect
    from pathlib import Path
    from tempfile import TemporaryDirectory
    errors = []
    calls = list(sdf_path_calls(document))
    if not calls:
        return errors
    with TemporaryDirectory(prefix="ck-stub-paths-") as folder:
        batch = module.MoleculeBatch.from_smiles_list([])
        path = str(Path(folder) / "empty.sdf")
        try:
            report = batch.write_sdf(path)
        except Exception as error:
            return [f"SDF path smoke setup failed: {error}"]
        for row, owner, fields in calls:
            try:
                receiver = getattr(module, owner)
                if row.get("receiver") is not None:
                    if owner == "MoleculeBatch":
                        receiver = batch
                    elif owner == "BatchExportReport":
                        receiver = report
                call = getattr(receiver, row["python_name"])
                arguments = {}
                for name, parameter in inspect.signature(call).parameters.items():
                    if name in ("path", "directory", "out_dir"):
                        arguments[name] = str(Path(folder) / "records") if "files" in row["python_name"] else path
                    elif name == "report_path":
                        arguments[name] = None
                    elif parameter.default is inspect.Parameter.empty:
                        field = next(field for field in row["parameters"] if field["name"] == name)
                        arguments[name] = None if compact(field["type"]).startswith("Option<") else getattr(module, compact(field["type"]).split("::")[-1])()
                for name in arguments.keys() & {"path", "directory", "out_dir", "report_path"}:
                    text = str(Path(folder) / "counts.json") if owner == "BatchExportReport" or name == "report_path" else arguments[name]
                    errors.extend(f"{row['semantic_id']}: {error}" for error in check_path_input(call, arguments, name, text))
            except Exception as error:
                errors.append(f"{row['semantic_id']}: path check setup failed: {error}")
    return errors


def configuration_calls(entries, keywords=()):
    """Resolve canonical configured/default pairs, never a handwritten list."""
    by_id = {row["semantic_id"]: row for row in entries}
    explicit_bases = {row["target"]: row["semantic_id"] for row in keywords}
    types = {compact(row["rust_path"]): row for row in entries if row.get("role") == "parameter"}
    for row in entries:
        if row["item"] != "callable":
            continue
        identifier = row["semantic_id"]
        base = identifier
        for suffix in ("_with_params_", "_with_options_", "_with_params", "_with_options"):
            if identifier.endswith(suffix):
                base = identifier[:-len(suffix)] + ("_" if suffix.endswith("_") else "")
                break
        if base == identifier:
            base = explicit_bases.get(identifier)
        if base not in by_id:
            continue
        configurations = [(field, types[compact(field["type"])]) for field in row.get("parameters") or [] if compact(field["type"]) in types]
        if configurations:
            yield by_id[base], row, configurations


def compact(value):
    import re
    return re.sub(r"&\s*(?:'\w+\s+)?(?:mut\s+)?", "", value).replace(" ", "")


def _functions(body, name):
    return [node for node in body if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name]


def _arguments(function):
    positional = function.args.posonlyargs + function.args.args
    defaults = [None] * (len(positional) - len(function.args.defaults)) + list(function.args.defaults)
    values = {arg.arg: (arg, default, False) for arg, default in zip(positional, defaults) if arg.arg not in ("self", "cls")}
    values.update({arg.arg: (arg, default, True) for arg, default in zip(function.args.kwonlyargs, function.args.kw_defaults)})
    return values


def _annotation(node):
    if node is None:
        return None
    value = ast.unparse(node).replace("builtins.", "").replace("typing.", "").replace(" ", "")
    return value


def _read_type_matches(output, input_type, classes=None):
    """A constructor may accept None/Sequence/Mapping and expose owned values.

    This is a representation conversion, not permission to accept a different
    scalar type or promise a narrower setter than the constructor.
    """
    # Enum inputs accept strings, but the read boundary must stay typed.
    # Only recognize this widening for a real declared enum vocabulary.
    if isinstance(input_type, ast.BinOp) and isinstance(input_type.op, ast.BitOr):
        def alternatives(node):
            if isinstance(node, ast.BinOp) and isinstance(node.op, ast.BitOr):
                return alternatives(node.left) + alternatives(node.right)
            return [node]
        options = alternatives(input_type)
        cls = (classes or {}).get(_annotation(output))
        native_enum = cls is not None and any(
            isinstance(member, ast.AnnAssign) and "ClassVar[" in (_annotation(member.annotation) or "")
            for member in cls.body)
        integer_enum = cls is not None and any(_annotation(base) in ("enum.IntEnum", "IntEnum") for base in cls.bases)
        names = {_annotation(option) for option in options}
        if (native_enum or integer_enum) and names == {_annotation(output), "str"}:
            return True
        if integer_enum and names == {_annotation(output), "str", "int"}:
            return True
    output = _annotation(output)
    input_type = _annotation(input_type)
    if input_type is None or output is None:
        return False
    if output == input_type:
        return True
    output_optional = output.startswith("Optional[") or output.endswith("|None")
    input_optional = input_type.startswith("Optional[") or input_type.endswith("|None")
    if output_optional and not input_optional:
        return False
    if output.startswith("Optional["):
        output = output[9:-1]
    elif output.endswith("|None"):
        output = output[:-5]
    if input_type.startswith("Optional["):
        input_type = input_type[9:-1]
    elif input_type.endswith("|None"):
        input_type = input_type[:-5]
    def owned_type(value):
        return value.replace("list[", "Sequence[").replace("List[", "Sequence[").replace("dict[", "Mapping[").replace("Dict[", "Mapping[")
    if owned_type(output) == owned_type(input_type):
        return True
    # A real IntEnum getter is a valid int projection. Do not treat PyO3
    # eq_int vocabulary classes (which are not int subclasses) the same way.
    cls = (classes or {}).get(output)
    return input_type == "int" and cls is not None and any(
        _annotation(base) in ("enum.IntEnum", "IntEnum") for base in cls.bases
    )


def _default(value):
    if value is None:
        return ("required",)
    try:
        return ("value", ast.literal_eval(value))
    except (ValueError, TypeError):
        return ("expression", _annotation(value))


def _registered_default(field):
    """Normalize the registry's existing literal/default declaration syntax."""
    value = field["default"]
    if value is None:
        return ("required",)
    value = value.strip()
    if value == "none":
        return ("value", None)
    for wrapper in ("integer", "boolean", "string"):
        if value.startswith(wrapper + "(") and value.endswith(")"):
            value = value[len(wrapper) + 1:-1]
            break
    value = value.replace("true", "True").replace("false", "False")
    try:
        result = ast.literal_eval(value)
        if isinstance(result, str) and compact(field["type"]) in ("bool", "usize", "u32", "i32", "f64", "u64"):
            result = ast.literal_eval(result)
        return ("value", result)
    except (SyntaxError, ValueError, TypeError):
        return ("expression", compact(value))


def check_configuration(module, document):
    """A Parameter's constructor defines its configurable input fields.

    Output/computed getters do not become writable just because they exist.
    Constructor/field/type/overload consistency is checked without inventing
    stubs or executing a chemical algorithm.
    """
    entries = document["entries"]
    classes = {node.name: node for node in module.body if isinstance(node, ast.ClassDef)}
    types = {row["semantic_id"].rsplit(".", 1)[-1]: row["python_name"] for row in entries if row["item"] == "type"}
    errors = []
    fields_by_id = {}
    for row in entries:
        if row.get("role") != "parameter":
            continue
        name = row["python_name"]
        if row.get("fields") is None:
            errors.append(f"{name}: Parameter requires a registered configuration constructor schema")
            continue
        cls = classes.get(name)
        if cls is None:
            continue  # The surface check already reports this missing type.
        if requires_configuration_repr(row):
            methods = _functions(cls.body, "__repr__")
            if (len(methods) != 1
                or _annotation(methods[0].returns) != "str"
                or [arg.arg for arg in methods[0].args.posonlyargs + methods[0].args.args] != ["self"]
                or methods[0].args.kwonlyargs or methods[0].args.vararg
                or methods[0].args.kwarg or methods[0].args.defaults):
                errors.append(f"{name}: batch configuration requires __repr__(self) -> str")
        constructors = _functions(cls.body, "__new__") or _functions(cls.body, "__init__")
        if len(constructors) != 1:
            errors.append(f"{name}: expected one canonical constructor declaration")
            continue
        arguments = _arguments(constructors[0])
        fields_by_id[row["semantic_id"]] = arguments
        expected = {field["name"] for field in row["fields"]}
        if set(arguments) != expected:
            errors.append(f"{name}: constructor fields differ from registry: expected {sorted(expected)}, got {sorted(arguments)}")
        for field in row["fields"]:
            key = field["name"]
            argument = arguments.get(key)
            if argument is None:
                continue
            annotation, default, _ = argument
            if row.get("python_field_types") and _annotation(annotation.annotation) != _annotation(ast.parse(field["type"], mode="eval").body):
                errors.append(f"{name}.{key}: constructor type differs from registry")
            # Rust literal defaults have an unambiguous language projection.
            source_default = field["default"]
            expected_default = (
                _default(ast.parse(source_default, mode="eval").body)
                if row.get("python_field_types") and source_default is not None
                else _registered_default(field)
            )
            if source_default is None and default is not None:
                errors.append(f"{name}.{key}: required registry input has a default")
            elif source_default is not None:
                # Symbolic defaults are not guessed. They must agree between
                # the constructor and every field-based call below; actual
                # getters are checked separately against literal defaults.
                if (row.get("python_field_types") or expected_default[0] == "value") and _default(default) != expected_default:
                    errors.append(f"{name}.{key}: default differs from registry ({source_default})")
            attributes = [node for node in cls.body if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name) and node.target.id == key]
            getters = [fn for fn in _functions(cls.body, key) if any(isinstance(decorator, ast.Name) and decorator.id == "property" for decorator in fn.decorator_list)]
            setters = [fn for fn in _functions(cls.body, key) if any(isinstance(decorator, ast.Attribute) and decorator.attr == "setter" and isinstance(decorator.value, ast.Name) and decorator.value.id == key for decorator in fn.decorator_list)]
            writable = [attr for attr in attributes if not any(tag in (_annotation(attr.annotation) or "") for tag in ("Final[", "ClassVar["))]
            if not writable and not (getters and setters):
                errors.append(f"{name}.{key}: configuration input must be readable and writable")
            write_types = [attr.annotation for attr in writable] + [arg.annotation for fn in setters for arg in fn.args.args if arg.arg != "self"]
            read_types = [attr.annotation for attr in writable] + [fn.returns for fn in getters]
            if not read_types or any(not _read_type_matches(ty, annotation.annotation, classes) for ty in read_types) or any(_annotation(ty) != _annotation(annotation.annotation) for ty in write_types):
                errors.append(f"{name}.{key}: constructor/getter/setter types differ")
            if "Callable[" in field["type"] and not field.get("callback"):
                errors.append(f"{name}.{key}: callable input requires an executable callback contract")
            for alias in field.get("aliases", []):
                aliases = [node.annotation for node in cls.body if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name) and node.target.id == alias]
                if len(aliases) != 1 or _annotation(aliases[0]) != _annotation(annotation.annotation):
                    errors.append(f"{name}.{alias}: registered configuration alias declaration missing or mistyped")

    for base, target, configurations in configuration_calls(entries, document.get("keywords", ())):
        owner = "Molecule" if base["owner"] == "molecule" else types.get(base["semantic_id"].rsplit(".", 1)[0], base["semantic_id"].rsplit(".", 1)[0])
        body = module.body if base["owner"] == "module" else classes.get(owner, ast.ClassDef(name=owner, bases=[], keywords=[], body=[], decorator_list=[])).body
        forms = _functions(body, base["python_name"])
        path = f"{owner}.{base['python_name']}"
        if not forms:
            continue
        for parameter, configuration in configurations:
            constructor_fields = fields_by_id.get(configuration["semantic_id"])
            if constructor_fields is None:
                continue
            parameter_forms = []
            keyword_forms = []
            for form in forms:
                arguments = _arguments(form)
                if any(_annotation(arg.annotation) in (configuration["python_name"], f"Optional[{configuration['python_name']}]", f"{configuration['python_name']}|None", f"None|{configuration['python_name']}") for arg, _, _ in arguments.values()):
                    parameter_forms.append(form)
                    if set(constructor_fields) & set(arguments):
                        errors.append(f"{path}: parameter instance and field configuration must be mutually exclusive overloads")
                if set(constructor_fields) <= set(arguments) and all(arguments[field][2] for field in constructor_fields):
                    keyword_forms.append(form)
                    for field, (arg, default, _) in constructor_fields.items():
                        candidate, actual_default, _ = arguments[field]
                        if _annotation(candidate.annotation) != _annotation(arg.annotation) or _default(actual_default) != _default(default):
                            errors.append(f"{path}({field}=...): type/default differs from {configuration['python_name']}")
                if form.args.kwarg is not None or form.args.vararg is not None:
                    errors.append(f"{path}: untyped variadic arguments cannot replace explicit configuration call forms")
            if not parameter_forms:
                errors.append(f"{path}: missing {configuration['python_name']} instance call form")
            if not keyword_forms:
                errors.append(f"{path}: missing keyword-only {configuration['python_name']} field call form")
            if len(forms) < 2 or any(not any(isinstance(d, ast.Name) and d.id == "overload" or isinstance(d, ast.Attribute) and d.attr == "overload" for d in form.decorator_list) for form in forms):
                errors.append(f"{path}: configuration call forms require explicit overload declarations")
            for field in configuration["fields"]:
                for alias in field.get("aliases", []):
                    if not any(alias in _arguments(form) and _annotation(_arguments(form)[alias][0].annotation) == _annotation(ast.parse(field["type"], mode="eval").body) for form in forms):
                        errors.append(f"{path}: registered configuration keyword alias {alias!r} missing or mistyped")
    return errors


def _assignment_probes(previous, field):
    """Candidate values are validated by the real constructor before assignment.

    No domain-specific configuration names or second registry live here.
    """
    if isinstance(previous, bool):
        return [not previous]
    if isinstance(previous, int):
        return [previous + 1, 1, 0]
    if isinstance(previous, float):
        return [previous + 0.5, 0.5, 1.0]
    if isinstance(previous, str):
        return [previous + "_ck_probe", ""]
    if isinstance(previous, list):
        return [[], list(previous) + list(previous[:1])]
    if isinstance(previous, dict):
        return [{}]
    if previous is None and "Sequence[" in field["type"]:
        return [[]]
    return []


def check_effective_assignment(cls, row, field, schemas):
    errors = []
    for probe in _assignment_probes(getattr(cls(), field["name"]), field):
        value = cls()
        original = _configuration_snapshot(value, schemas)
        arguments = {item["name"]: getattr(value, item["name"]) for item in row["fields"]}
        try:
            candidate = cls(**dict(arguments, **{field["name"]: probe}))
        except (TypeError, ValueError, OverflowError):
            continue  # A source-defined constructor constraint, not a setter defect.
        expected = _configuration_snapshot(candidate, schemas)
        if expected == original:
            continue
        try:
            setattr(value, field["name"], probe)
            if _configuration_snapshot(value, schemas) != expected:
                errors.append(f"{row['python_name']}.{field['name']}: effective assignment differs from constructor (ignored input or changed sibling fields)")
        except Exception as error:
            errors.append(f"{row['python_name']}.{field['name']}: valid constructor input rejected by setter: {error}")
        break
    return errors


def check_callback_runtime(module, document):
    """Exercise registered callbacks through every matching call form.

    Two atoms are sufficient to detect missing dispatch; this is a binding
    contract probe, not a reference/corpus test or a new chemistry oracle.
    """
    errors = []
    for base, target, configurations in configuration_calls(document["entries"], document.get("keywords", ())):
        if base["owner"] != "molecule":
            continue
        for parameter, row in configurations:
            for field in row["fields"]:
                if not field.get("callback"):
                    continue
                path = f"{base['semantic_id']}({field['name']})"
                try:
                    molecule = module.Molecule.from_smiles("CC")
                    query = module.parse_smarts("C~C")
                    baseline = getattr(molecule, base["python_name"])(query)
                    if not baseline:
                        raise AssertionError("callback probe has no baseline match")
                    before = molecule.to_smiles()
                    for form, name in [("constructor", field["name"]), ("assignment", field["name"]), ("parameter", field["name"]), ("keyword", field["name"])] + [(form, alias) for alias in field.get("aliases", []) for form in ("assignment", "keyword")]:
                        for behavior in ("accept", "reject", "raise"):
                            seen = []
                            original_error = RuntimeError("binding callback probe")

                            def callback(*args):
                                kind = field["callback"]
                                if len(args) != 2:
                                    raise AssertionError("callback arity differs from its contract")
                                if kind == "final_match":
                                    if not isinstance(args[0], module.Molecule) or len(args[1]) != 2 or any(type(index) is not int or not 0 <= index < 2 for index in args[1]):
                                        raise AssertionError("invalid final-match callback arguments")
                                elif kind == "atom_match":
                                    if not isinstance(args[0], module.QueryAtom) or not isinstance(args[1], module.Atom):
                                        raise AssertionError("invalid atom-match callback arguments")
                                elif kind == "bond_match":
                                    if not all(isinstance(arg, module.Bond) for arg in args):
                                        raise AssertionError("invalid bond-match callback arguments")
                                else:
                                    raise AssertionError("unknown callback contract")
                                seen.append(args)
                                if behavior == "raise":
                                    raise original_error
                                return behavior == "accept"

                            cls = getattr(module, row["python_name"])
                            if form == "keyword":
                                call = lambda: getattr(molecule, base["python_name"])(query, **{name: callback})
                            else:
                                params = cls(**{name: callback}) if form in ("constructor", "parameter") else cls()
                                if form == "assignment":
                                    setattr(params, name, callback)
                                if getattr(params, name) is not callback:
                                    raise AssertionError(f"{name}: callback assignment was discarded")
                                entry = base if form == "parameter" else target
                                call = lambda: getattr(molecule, entry["python_name"])(query, params)
                            try:
                                result = call()
                            except Exception as error:
                                if behavior != "raise" or error is not original_error:
                                    raise AssertionError(f"{form}/{behavior}: original callback exception not propagated") from error
                            else:
                                def snapshot(value):
                                    if isinstance(value, list):
                                        return [snapshot(item) for item in value]
                                    if value is None or type(value) is bool:
                                        return value
                                    return tuple(value.atom_mapping())
                                expected = snapshot(baseline) if behavior == "accept" else [] if isinstance(baseline, list) else False if type(baseline) is bool else None
                                if behavior == "raise" or snapshot(result) != expected:
                                    raise AssertionError(f"{form}/{behavior}: callback result ignored")
                            if not seen or behavior == "raise" and len(seen) != 1:
                                raise AssertionError(f"{form}/{behavior}: callback not invoked or invoked again after exception")
                            if molecule.to_smiles() != before:
                                raise AssertionError("read-only matching changed the receiver")
                except Exception as error:
                    errors.append(f"{path}: callback contract failed: {error}")
    return errors


def _configuration_snapshot(value, schemas):
    """Compare registered values, not identities of fresh owned getter copies."""
    cls = type(value)
    if cls in schemas:
        return (cls, tuple((field["name"], _configuration_snapshot(getattr(value, field["name"]), schemas)) for field in schemas[cls]))
    if isinstance(value, float):
        import struct
        return (cls, struct.pack("!d", value))
    if isinstance(value, (list, tuple)):
        return (cls, tuple(_configuration_snapshot(item, schemas) for item in value))
    if isinstance(value, dict):
        return (cls, {key: _configuration_snapshot(item, schemas) for key, item in value.items()})
    return (cls, value)


def requires_configuration_repr(row):
    return row.get("role") == "parameter" and (row["feature"] == "cap-batch" or "cap-batch" in row.get("required_capabilities", []))


def check_configuration_repr(value, row, schemas):
    """All registered fields must reflect live values without changing state."""
    import inspect
    name = row["python_name"]
    if inspect.getattr_static(type(value), "__repr__") is object.__repr__:
        return [f"{name}: actual batch configuration __repr__ missing"]
    try:
        before = _configuration_snapshot(value, schemas)
        text = repr(value)
        expected = f"{name}(" + ", ".join(f"{field['name']}={getattr(value, field['name'])!r}" for field in row["fields"]) + ")"
        errors = []
        if text != expected:
            errors.append(f"{name}: repr must display all registered field values")
        if _configuration_snapshot(value, schemas) != before:
            errors.append(f"{name}: repr changed configuration state")
        return errors
    except Exception as error:
        return [f"{name}: configuration repr failed: {error}"]


def check_runtime(module, document):
    """Check actual extension exports and setters, not a promised stub.

    Default-constructible configuration values and an empty batch are used.
    Filesystem smoke checks use disposable empty files; invalid/protocol probes
    stop before IO. Callback checks use a two-atom in-memory molecule; no
    external fixture or reference installation is used.
    """
    import inspect
    document = python_document(document)
    entries = document["entries"]
    types = {row["semantic_id"].rsplit(".", 1)[-1]: row["python_name"] for row in entries if row["item"] == "type"}
    schemas = {
        getattr(module, row["python_name"]): row["fields"]
        for row in entries
        if row.get("role") == "parameter" and row.get("fields") is not None
        and isinstance(getattr(module, row["python_name"], None), type)
    }
    errors = []
    for row in entries:
        if row.get("python_native") is not None:
            import builtins
            native = row["python_native"]
            actual_types = {name: value for name, value in vars(module).items() if isinstance(value, type)}
            if not _native_projection_valid(native, actual_types) or native.startswith("builtins.") and not isinstance(getattr(builtins, native.split(".")[-1], None), type):
                errors.append(f"{row['semantic_id']}: invalid actual native value projection")
            continue
        owner = module if row["owner"] == "module" or row["item"] == "type" else getattr(module, "Molecule" if row["owner"] == "molecule" else types.get(row["semantic_id"].rsplit(".", 1)[0], ""), None)
        # Enum.name is an instance property whose class lookup intentionally
        # raises AttributeError. Inspect its real descriptor without invoking it.
        present = owner is not None and (
            inspect.getattr_static(owner, row["python_name"], None) is not None
            if row.get("python_property") is not None
            else hasattr(owner, row["python_name"])
        )
        if not present:
            errors.append(f"{row['semantic_id']}: missing actual Python export")
        if row.get("role") != "parameter" or row.get("fields") is None:
            continue
        cls = getattr(module, row["python_name"], None)
        if cls is None:
            continue
        if requires_configuration_repr(row) and inspect.getattr_static(cls, "__repr__") is object.__repr__:
            errors.append(f"{row['python_name']}: actual batch configuration __repr__ missing")
        for field in row["fields"]:
            descriptor = inspect.getattr_static(cls, field["name"], None)
            if descriptor is None or not hasattr(descriptor, "__set__") or isinstance(descriptor, property) and descriptor.fset is None:
                errors.append(f"{row['python_name']}.{field['name']}: actual setter missing")
        if any(field["default"] is None for field in row["fields"]):
            continue
        try:
            value = cls()
        except Exception as error:
            errors.append(f"{row['python_name']}: default construction failed: {error}")
            continue
        try:
            setattr(value, "ck_unknown_configuration_field", object())
        except AttributeError:
            pass
        except Exception as error:
            errors.append(f"{row['python_name']}: unknown configuration field raised the wrong error: {error}")
        else:
            errors.append(f"{row['python_name']}: unknown configuration field was silently accepted")
        if requires_configuration_repr(row):
            errors.extend(check_configuration_repr(value, row, schemas))
        for field in row["fields"]:
            try:
                previous = getattr(value, field["name"])
                before = _configuration_snapshot(value, schemas)
                expected = _registered_default(field)
                if expected[0] == "value" and expected[1] is not None and previous != expected[1]:
                    errors.append(f"{row['python_name']}.{field['name']}: actual default differs from registry")
                setattr(value, field["name"], previous)
                if _configuration_snapshot(value, schemas) != before:
                    errors.append(f"{row['python_name']}.{field['name']}: assignment did not preserve the value")
                errors.extend(check_effective_assignment(cls, row, field, schemas))
                vocabulary = getattr(type(previous), "_enum_string_values", None)
                if vocabulary is None and hasattr(type(previous), "__members__"):
                    vocabulary = [(name.lower(), member) for name, member in type(previous).__members__.items()]
                if vocabulary is not None:
                    arguments = {item["name"]: getattr(value, item["name"]) for item in row["fields"]}
                    for text, member in vocabulary:
                        candidate = cls(**dict(arguments, **{field["name"]: text}))
                        if type(getattr(candidate, field["name"])) is not type(member) or getattr(candidate, field["name"]) != member:
                            errors.append(f"{row['python_name']}.{field['name']}: string constructor differs from enum")
                        setattr(value, field["name"], text)
                        if type(getattr(value, field["name"])) is not type(member) or getattr(value, field["name"]) != member:
                            errors.append(f"{row['python_name']}.{field['name']}: string assignment differs from enum")
                    setattr(value, field["name"], previous)
                    for invalid, category in [("__invalid_enum_value__", ValueError), (object(), TypeError)]:
                        try:
                            setattr(value, field["name"], invalid)
                        except category:
                            if _configuration_snapshot(value, schemas) != before:
                                errors.append(f"{row['python_name']}.{field['name']}: rejected enum assignment changed state")
                        else:
                            errors.append(f"{row['python_name']}.{field['name']}: invalid enum input was not rejected")
                if compact(field["type"]) in ("bool", "usize", "u32", "i32", "f64", "String", "&str", "builtins.bool", "builtins.int", "builtins.float", "builtins.str", "builtins.bytes"):
                    try:
                        setattr(value, field["name"], object())
                    except Exception:
                        if _configuration_snapshot(value, schemas) != before:
                            errors.append(f"{row['python_name']}.{field['name']}: failed assignment changed the value")
                    else:
                        errors.append(f"{row['python_name']}.{field['name']}: setter accepted an invalid scalar type")
            except Exception as error:
                errors.append(f"{row['python_name']}.{field['name']}: actual assignment failed: {error}")
    errors.extend(check_path_runtime(module, document))
    errors.extend(check_callback_runtime(module, document))
    errors.extend(check_alias_runtime(module, document))
    errors.extend(check_collection_runtime(module, document))
    return errors
