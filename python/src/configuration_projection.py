"""Canonical language call forms, derived from the linked binding registry.

Only argument normalization and configuration display live here. Native
constructors validate options; the existing native configured method performs
the operation. The registry is
the sole API inventory. This module is embedded in the extension, not imported
from a developer checkout at runtime.
"""
import ast
import copy
import functools
import inspect
import re
import reprlib
import weakref


def _detach(value):
    """A replaced field's old view becomes an ordinary independent value."""
    if isinstance(value, (_ConfigurationList, _ConfigurationDict)):
        value._read = value._write = None
        for child in value.values() if isinstance(value, dict) else value:
            _detach(child)
    elif hasattr(type(value), "_configuration_descriptors"):
        value.__dict__.pop("_configuration_link", None)


def _refresh(value, fresh):
    if type(value) is type(fresh) and hasattr(type(value), "_configuration_descriptors"):
        value._configuration_replace(fresh)
        _refresh_fields(value)
        return value
    if isinstance(value, _ConfigurationList) and isinstance(fresh, list):
        value._refresh(fresh)
        return value
    if isinstance(value, _ConfigurationDict) and isinstance(fresh, dict):
        value._refresh(fresh)
        return value
    _detach(value)
    return fresh


def _refresh_fields(owner):
    cached = owner.__dict__.get("_configuration_views", {})
    for name, view in list(cached.items()):
        descriptor = type(owner)._configuration_descriptors[name]
        fresh = descriptor.__get__(owner, type(owner))
        updated = _refresh(view, fresh)
        if updated is not view:
            del cached[name]


def _linked(value, read, write):
    if hasattr(type(value), "_configuration_descriptors"):
        value.__dict__["_configuration_link"] = (read, write)
    elif isinstance(value, list):
        value = _ConfigurationList(value, read, write)
    elif isinstance(value, dict):
        value = _ConfigurationDict(value, read, write)
    return value


class _ConfigurationList(list):
    """A Python list whose edits commit through the native field setter."""
    def __init__(self, values, read, write):
        self._read, self._write = read, write
        list.__init__(self)
        self._refresh(values)

    def _refresh(self, values):
        old = list(self)
        children = []
        parent = weakref.ref(self)
        for index, value in enumerate(values):
            def read(index=index):
                view = parent()
                return view._read()[index] if view._read else list.__getitem__(view, index)

            def write(value, index=index):
                view = parent()
                view[index] = value

            if index < len(old):
                value = _refresh(old[index], value)
            children.append(_linked(value, read, write))
        for value in old[len(values):]:
            _detach(value)
        list.__setitem__(self, slice(None), children)

    def _mutate(self, method, *args, **kwargs):
        # Work on a snapshot first: conversion/validation failures change neither
        # the parent Rust record nor this list (including extend/slice updates).
        candidate = self._read() if self._read else list(self)
        result = getattr(list, method)(candidate, *args, **kwargs)
        if self._write:
            self._write(candidate)
        else:
            self._refresh(candidate)
        return self if method in ("__iadd__", "__imul__") else result


class _ConfigurationDict(dict):
    def __init__(self, values, read, write):
        self._read, self._write = read, write
        dict.__init__(self)
        self._refresh(values)

    def _refresh(self, values):
        parent = weakref.ref(self)
        children = {}
        for key, value in values.items():
            def read(key=key):
                view = parent()
                return view._read()[key] if view._read else dict.__getitem__(view, key)

            def write(value, key=key):
                view = parent()
                view[key] = value

            if key in self:
                value = _refresh(dict.__getitem__(self, key), value)
            children[key] = _linked(value, read, write)
        for key in self.keys() - values.keys():
            _detach(dict.__getitem__(self, key))
        dict.clear(self)
        dict.update(self, children)

    def _mutate(self, method, *args, **kwargs):
        candidate = self._read() if self._read else dict(self)
        result = getattr(dict, method)(candidate, *args, **kwargs)
        if self._write:
            self._write(candidate)
        else:
            self._refresh(candidate)
        return self if method == "__ior__" else result


def _container_mutator(method):
    def mutate(self, *args, **kwargs):
        return self._mutate(method, *args, **kwargs)
    return mutate


for _method in ("append", "extend", "insert", "pop", "remove", "clear", "reverse", "sort", "__setitem__", "__delitem__", "__iadd__", "__imul__"):
    setattr(_ConfigurationList, _method, _container_mutator(_method))
for _method in ("update", "setdefault", "pop", "popitem", "clear", "__setitem__", "__delitem__", "__ior__"):
    setattr(_ConfigurationDict, _method, _container_mutator(_method))


def _live_property(name, descriptor):
    def get(owner):
        value = descriptor.__get__(owner, type(owner))
        if not isinstance(value, (list, dict)) and not hasattr(type(value), "_configuration_descriptors"):
            return value
        cached = owner.__dict__.setdefault("_configuration_views", {})
        if name not in cached:
            child = None

            def released(_):
                view = child() if child is not None else None
                if view is not None:
                    _detach(view)

            parent = weakref.ref(owner, released)

            def read():
                current = parent()
                if current is None:
                    raise ReferenceError("configuration owner no longer exists")
                return descriptor.__get__(current, type(current))

            def write(value):
                current = parent()
                if current is None:
                    raise ReferenceError("configuration owner no longer exists")
                setattr(current, name, value)

            cached[name] = _linked(value, read, write)
            child = weakref.ref(cached[name])
        return cached[name]

    def set_value(owner, value):
        link = owner.__dict__.get("_configuration_link")
        if link:
            read, write = link
            candidate = read()
            descriptor.__set__(candidate, value)
            write(candidate)
        else:
            descriptor.__set__(owner, value)
            _refresh_fields(owner)

    return property(get, set_value, doc=descriptor.__doc__)


def _compact(value):
    return re.sub(r"&\s*(?:'\w+\s+)?(?:mut\s+)?", "", value).replace(" ", "")


def pairs(entries, keywords=()):
    by_id = {row["semantic_id"]: row for row in entries}
    explicit_bases = {row["target"]: row["semantic_id"] for row in keywords}
    types = {_compact(row["rust_path"]): row for row in entries if row.get("role") == "parameter"}
    for target in entries:
        if target["item"] != "callable":
            continue
        identifier = target["semantic_id"]
        for suffix in ("_with_params_", "_with_options_", "_with_params", "_with_options"):
            if identifier.endswith(suffix):
                base = identifier[:-len(suffix)] + ("_" if suffix.endswith("_") else "")
                break
        else:
            # Some canonical methods already take configuration, with no
            # separate short/_with_params pair. Their registry opts them in.
            base = explicit_bases.get(identifier)
            if base is None:
                continue
        if base not in by_id:
            continue
        configurations = [(field, types[_compact(field["type"])]) for field in target["parameters"] or [] if _compact(field["type"]) in types]
        if configurations:
            yield by_id[base], target, configurations


def _owner(module, entry, names):
    if entry["owner"] == "module":
        return module
    key = "Molecule" if entry["owner"] == "molecule" else entry["semantic_id"].rsplit(".", 1)[0]
    return getattr(module, names.get(key, key))


def _parameter_type(field, records):
    annotation = field["type"].replace(" ", "")
    if annotation.startswith("typing.Optional["):
        annotation = annotation[len("typing.Optional["):-1]
    return records.get(annotation)


def _keyword_options(configurations, entries):
    """Expose only unambiguous registered fields, including nested records."""
    records = {row["python_name"]: row for row in entries if row.get("role") == "parameter"}

    def walk(row, path=(), seen=()):
        for field in row["python_fields"]:
            yield path + (field["name"],), field
            for alias in field.get("aliases", []):
                yield path + (field["name"],), dict(field, name=alias)
            nested = _parameter_type(field, records)
            if nested is not None and nested["python_name"] not in seen:
                yield from walk(nested, path + (field["name"],), seen + (row["python_name"],))

    choices = [(owner["name"], path, field) for owner, row in configurations for path, field in walk(row)]
    direct = {field["name"] for _, path, field in choices if len(path) == 1}
    counts = {}
    for _, path, field in choices:
        counts[field["name"]] = counts.get(field["name"], 0) + 1
    return {owner["name"]: {
        field["name"]: (path, field) for name, path, field in choices if name == owner["name"]
        and (len(path) == 1 or field["name"] not in direct and counts[field["name"]] == 1)
    } for owner, _ in configurations}, records


def _configuration_value(row, supplied, module, records):
    paths = [path for path, _ in supplied]
    if len(paths) != len(set(paths)):
        raise TypeError("configuration field and its alias are mutually exclusive")
    direct = {path[0]: value for path, value in supplied if len(path) == 1}
    nested = {}
    for path, value in supplied:
        if len(path) > 1:
            nested.setdefault(path[0], []).append((path[1:], value))
    fields = {field["name"]: field for field in row["python_fields"]}
    for name, values in nested.items():
        if direct.get(name) is not None:
            raise TypeError(f"{name} and its configuration keywords are mutually exclusive")
        direct[name] = _configuration_value(_parameter_type(fields[name], records), values, module, records)
    return getattr(module, row["python_name"])(**direct)


def _wrapper(original, configured, target, configurations, module, entries):
    signature = inspect.signature(configured)
    classes = {field["name"]: getattr(module, row["python_name"]) for field, row in configurations}
    options, records = _keyword_options(configurations, entries)
    fields = {name: set(values) for name, values in options.items()}
    rows = {field["name"]: row for field, row in configurations}
    all_fields = set().union(*fields.values())
    if sum(map(len, fields.values())) != len(all_fields):
        raise TypeError(f"{target['semantic_id']}: ambiguous configuration fields")
    source_signature = inspect.signature(original)
    original_data = [p.name for p in source_signature.parameters.values() if p.name not in ("self", "cls") and p.name not in classes and p.name not in all_fields]
    target_data = [p.name for p in signature.parameters.values() if p.name not in ("self", "cls") and p.name not in classes]
    # Reordered data parameters retain their own names. Only genuinely renamed
    # inputs need an alias; positional zipping would swap filenames/report_path.
    renamed_source = [name for name in original_data if name not in target_data]
    renamed_target = [name for name in target_data if name not in original_data]
    aliases = dict(zip(renamed_source, renamed_target)) if len(renamed_source) == len(renamed_target) else {}
    optional = {field["name"] for field in target["parameters"] if _compact(field["type"]).startswith("Option<")}
    data_defaults = {aliases.get(name, name): value.default
                     for name, value in source_signature.parameters.items()
                     if name not in classes and name not in all_fields
                     and value.default is not inspect.Parameter.empty}
    configured_only = signature.parameters.keys() - source_signature.parameters.keys()

    @functools.wraps(original)
    def call(*args, **kwargs):
        configured_call = bool(set(kwargs) & (all_fields | classes.keys() | configured_only)) or any(isinstance(value, tuple(classes.values())) for value in args)
        if not configured_call:
            return original(*args, **kwargs)
        values = {key: value for key, value in kwargs.items() if key not in all_fields}
        for old, new in aliases.items():
            if old != new and old in values:
                if new in values:
                    raise TypeError(f"multiple values for {new!r}")
                values[new] = values.pop(old)
        bound = signature.bind_partial(*args, **values)
        for name, cls in classes.items():
            keywords = {key: kwargs[key] for key in fields[name] if key in kwargs}
            value = bound.arguments.get(name)
            if value is not None:
                if not isinstance(value, cls):
                    raise TypeError(f"{name} must be {cls.__name__}")
                if keywords:
                    raise TypeError(f"{name} and its configuration keywords are mutually exclusive")
            else:
                value = _configuration_value(rows[name], [(options[name][key][0], value) for key, value in keywords.items()], module, records)
            bound.arguments[name] = value
        for name in optional:
            bound.arguments.setdefault(name, None)
        for name, value in data_defaults.items():
            if name in signature.parameters:
                bound.arguments.setdefault(name, value)
        # Complete binding rejects missing required data before native dispatch.
        complete = signature.bind(*bound.args, **bound.kwargs)
        return configured(*complete.args, **complete.kwargs)

    call._configuration_contract = (target["semantic_id"], tuple((field["name"], row["semantic_id"]) for field, row in configurations))
    return call


def install(module, document):
    entries = document["entries"]
    names = {row["semantic_id"].removeprefix("types."): row["python_name"] for row in entries if row["item"] == "type"}
    # Preserve native classes and their extraction/validation. Only configuration
    # getters become live views; molecule/result snapshots remain independent.
    for row in entries:
        if row.get("role") != "parameter":
            continue
        cls = getattr(module, row["python_name"])
        if not hasattr(cls, "_configuration_replace"):
            raise TypeError(f"{row['python_name']}: native configuration view support missing")
        descriptors = {field["name"]: inspect.getattr_static(cls, field["name"]) for field in row["python_fields"]}
        cls._configuration_descriptors = descriptors
        allowed = set(descriptors)
        for name, descriptor in descriptors.items():
            setattr(cls, name, _live_property(name, descriptor))
        for field in row["python_fields"]:
            for alias in field.get("aliases", []):
                allowed.add(alias)
                def get(owner, name=field["name"]):
                    return getattr(owner, name)
                def set_value(owner, value, name=field["name"]):
                    setattr(owner, name, value)
                setattr(cls, alias, property(get, set_value))
        # __dict__ remains an implementation detail for live views/weakrefs;
        # arbitrary public attributes must never masquerade as native options.
        def set_attribute(owner, name, value, allowed=frozenset(allowed)):
            if name not in allowed:
                raise AttributeError(f"{type(owner).__name__}: unknown configuration field {name!r}")
            object.__setattr__(owner, name, value)
        cls.__setattr__ = set_attribute
    for row in entries:
        if row.get("role") == "parameter" and (row["feature"] == "cap-batch" or "cap-batch" in row.get("required_capabilities", [])):
            fields = tuple(field["name"] for field in row["python_fields"])
            setattr(getattr(module, row["python_name"]), "_configuration_repr", _configuration_repr(fields))
    for base, target, configurations in pairs(entries, document.get("keywords", ())):
        owner = _owner(module, base, names)
        original = getattr(owner, base["python_name"])
        configured = getattr(owner, target["python_name"])
        call = _wrapper(original, configured, target, configurations, module, entries)
        if owner is not module and "self" not in inspect.signature(original).parameters:
            call = staticmethod(call)
        setattr(owner, base["python_name"], call)
    # Install only after all configuration forms have been normalized. Each
    # export is the actual existing callable, never another dispatch wrapper.
    for alias in document.get("python_aliases", []):
        owner, method = alias["target"].split(".")
        setattr(module, alias["name"], getattr(getattr(module, owner), method))
        if hasattr(module, "__all__") and alias["name"] not in module.__all__:
            module.__all__ = [*module.__all__, alias["name"]]


def _configuration_repr(fields):
    @reprlib.recursive_repr(fillvalue="...")
    def describe(self) -> str:
        values = ", ".join(f"{field}={getattr(self, field)!r}" for field in fields)
        return f"{type(self).__name__}({values})"
    return describe


def declarations(module, stub, document):
    """Emit overload metadata only for installed, real normalization wrappers."""
    stub = enum_input_declarations(module, stub, document)
    tree = ast.parse(stub)
    classes = {node.name: node for node in tree.body if isinstance(node, ast.ClassDef)}
    alias_declarations = []
    for row in document["entries"]:
        if row.get("role") != "parameter":
            continue
        cls = classes[row["python_name"]]
        for field in row["python_fields"]:
            for alias in field.get("aliases", []):
                descriptor = inspect.getattr_static(getattr(module, row["python_name"]), alias)
                if not isinstance(descriptor, property) or descriptor.fset is None:
                    raise TypeError(f"{row['python_name']}.{alias}: real alias property missing")
                alias_declarations.append((cls.end_lineno, f"    {alias}: {field['type']}\n"))
    names = {row["semantic_id"].removeprefix("types."): row["python_name"] for row in document["entries"] if row["item"] == "type"}
    replacements = []
    for base, target, configurations in pairs(document["entries"], document.get("keywords", ())):
        options, _ = _keyword_options(configurations, document["entries"])
        owner = _owner(module, base, names)
        function = getattr(owner, base["python_name"])
        expected = (target["semantic_id"], tuple((field["name"], row["semantic_id"]) for field, row in configurations))
        if getattr(function, "_configuration_contract", None) != expected:
            raise TypeError(f"{base['semantic_id']}: actual configuration wrapper missing")
        scope = tree.body if owner is module else classes[owner.__name__].body
        originals = [node for node in scope if isinstance(node, ast.FunctionDef) and node.name == base["python_name"]]
        targets = [node for node in scope if isinstance(node, ast.FunctionDef) and node.name == target["python_name"]]
        if not originals or len(targets) != 1:
            raise TypeError(f"{base['semantic_id']}: native declaration unavailable")
        source = originals[0]
        explicit = copy.deepcopy(targets[0])
        explicit.name = source.name
        keywords = copy.deepcopy(targets[0])
        keywords.name = source.name
        config_names = {field["name"] for field, row in configurations}
        config_fields = {key for values in options.values() for key in values}
        default_configs = {
            field["name"]: row["python_name"]
            for field, row in configurations
            if all(value["default"] is not None for value in row["python_fields"])
        }
        # Carry the real short-form defaults of non-configuration inputs into
        # the configured forms (batch error mode, progress selection, etc.).
        source_args = source.args.posonlyargs + source.args.args
        source_defaults = dict(zip(
            [arg.arg for arg in source_args[len(source_args) - len(source.args.defaults):]],
            source.args.defaults,
        ))
        source_defaults.update({arg.arg: default for arg, default in zip(source.args.kwonlyargs, source.args.kw_defaults) if default is not None})
        optional_inputs = {field["name"] for field in target["parameters"] if _compact(field["type"]).startswith("Option<")}
        for form in (explicit, keywords):
            args = form.args.posonlyargs + form.args.args
            defaults = [None] * (len(args) - len(form.args.defaults)) + form.args.defaults
            defaults = [source_defaults.get(arg.arg, default) if arg.arg not in config_names | config_fields else default for arg, default in zip(args, defaults)]
            # The runtime wrapper supplies None for omitted optional inputs,
            # including Morgan's additional_output; the overload must agree.
            defaults = [ast.Constant(value=None) if default is None and arg.arg in optional_inputs else default for arg, default in zip(args, defaults)]
            if form is explicit:
                # The installed wrapper constructs omitted defaultable records,
                # including when only the execution params object is supplied.
                for index, arg in enumerate(args):
                    if arg.arg in default_configs:
                        arg.annotation = ast.parse(f"typing.Optional[{default_configs[arg.arg]}]", mode="eval").body
                        defaults[index] = ast.Constant(value=None)
            # A required config following optional data is keyword-only,
            # never accidentally assigned the preceding data's default.
            seen_default = False
            boundary = len(args)
            for index, default in enumerate(defaults):
                if default is None and seen_default:
                    boundary = index
                    break
                seen_default |= default is not None
            if boundary < len(args):
                positional_only = len(form.args.posonlyargs)
                form.args.posonlyargs = args[:min(boundary, positional_only)]
                form.args.args = args[positional_only:boundary]
                form.args.kwonlyargs = args[boundary:] + form.args.kwonlyargs
                form.args.kw_defaults = defaults[boundary:] + form.args.kw_defaults
            form.args.defaults = [default for default in defaults[:boundary] if default is not None]

        def remove_arguments(node, removed):
            positional = node.args.posonlyargs + node.args.args
            defaults = [None] * (len(positional) - len(node.args.defaults)) + node.args.defaults
            kept = [(arg, default) for arg, default in zip(positional, defaults) if arg.arg not in removed]
            posonly = {arg.arg for arg in node.args.posonlyargs}
            node.args.posonlyargs = [arg for arg, default in kept if arg.arg in posonly]
            node.args.args = [arg for arg, default in kept if arg.arg not in posonly]
            node.args.defaults = [default for arg, default in kept if default is not None]
            kept_keywords = [(arg, default) for arg, default in zip(node.args.kwonlyargs, node.args.kw_defaults) if arg.arg not in removed]
            node.args.kwonlyargs = [arg for arg, default in kept_keywords]
            node.args.kw_defaults = [default for arg, default in kept_keywords]

        default = copy.deepcopy(source)
        # Keep required positional inputs of the existing default API (e.g.
        # fragment atoms), but not conflicting convenience configuration defaults.
        remove_arguments(default, config_names | {arg.arg for arg in default.args.kwonlyargs if arg.arg in config_fields})
        remove_arguments(keywords, config_names)
        for values in options.values():
            for _, field in values.values():
                keywords.args.kwonlyargs.append(ast.arg(arg=field["name"], annotation=ast.parse(field["type"], mode="eval").body))
                keywords.args.kw_defaults.append(None if field["default"] is None else ast.parse(field["default"], mode="eval").body)
        forms = [default, explicit, keywords]
        # The configured native factory may be a staticmethod while the short
        # form is a classmethod. Its overload still needs the class receiver.
        if any(isinstance(node, ast.Name) and node.id == "classmethod" for node in source.decorator_list):
            receiver = source.args.posonlyargs[:1] or source.args.args[:1]
            for form in (explicit, keywords):
                if receiver and not any(arg.arg == receiver[0].arg for arg in form.args.posonlyargs + form.args.args):
                    form.args.args.insert(0, copy.deepcopy(receiver[0]))
        for form in forms:
            form.decorator_list = [node for node in source.decorator_list if not (isinstance(node, ast.Name) and node.id == "overload" or isinstance(node, ast.Attribute) and node.attr == "overload")]
            form.decorator_list.insert(0, ast.parse("typing.overload", mode="eval").body)
        indentation = "" if owner is module else "    "
        text = "\n".join(indentation + line for form in forms for line in ast.unparse(ast.fix_missing_locations(form)).splitlines()) + "\n"
        start = min(min([node.lineno] + [d.lineno for d in node.decorator_list]) for node in originals) - 1
        end = max(node.end_lineno for node in originals)
        replacements.append((start, end, text))
    lines = stub.splitlines(keepends=True)
    replacements.extend((line, line, text) for line, text in alias_declarations)
    for start, end, text in sorted(replacements, reverse=True):
        lines[start:end] = [text]
    stub = "".join(lines)
    for alias in document.get("python_aliases", []):
        stub += f"\n{alias['name']} = {alias['target']}\n"
    if document.get("python_aliases"):
        stub += f"\n__all__ += {[alias['name'] for alias in document['python_aliases']]!r}\n"
    return stub


def enum_input_declarations(module, stub, document):
    """Describe native enum extraction, without changing any output type."""
    enums = {
        name for name, value in vars(module).items()
        if isinstance(value, type) and hasattr(value, "_enum_string_values")
    }
    tree = ast.parse(stub)
    configurations = {
        row["python_name"]: {field["name"] for field in row["python_fields"]}
        for row in document["entries"]
        if row.get("role") == "parameter" and row.get("python_fields") is not None
    }

    class InputType(ast.NodeTransformer):
        def visit_Name(self, node):
            if node.id in enums:
                return ast.BinOp(left=node, op=ast.BitOr(), right=ast.parse("builtins.str", mode="eval").body)
            return node

    def input_type(annotation):
        return InputType().visit(copy.deepcopy(annotation))

    for scope in [tree] + [node for node in tree.body if isinstance(node, ast.ClassDef)]:
        body = []
        fields = configurations.get(getattr(scope, "name", None), set())
        for node in scope.body:
            if isinstance(node, ast.FunctionDef):
                # Return annotations, including enum-valued property getters,
                # intentionally remain untouched.
                for argument in node.args.posonlyargs + node.args.args + node.args.kwonlyargs:
                    if argument.annotation is not None:
                        argument.annotation = input_type(argument.annotation)
            elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name) and node.target.id in fields:
                expanded = input_type(node.annotation)
                if ast.dump(expanded) != ast.dump(node.annotation):
                    # A writable attribute cannot express a wider setter than
                    # getter; real native descriptors can, so describe them.
                    name = node.target.id
                    getter = ast.parse(f"@property\ndef {name}(self) -> {ast.unparse(node.annotation)}: ...").body[0]
                    setter = ast.parse(f"@{name}.setter\ndef {name}(self, value: {ast.unparse(expanded)}) -> None: ...").body[0]
                    body.extend((getter, setter))
                    continue
            body.append(node)
        scope.body = body
    return ast.unparse(ast.fix_missing_locations(tree)) + "\n"
