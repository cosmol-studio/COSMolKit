"""Check all compiled registry projections; never build or select an API subset."""
from __future__ import annotations

import argparse
import ast
from collections import Counter
import importlib
import inspect
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent


class RegistryMismatch(ValueError):
    def __init__(self, report):
        self.report = report
        super().__init__("\n".join(item["message"] for item in report["errors"]))


def expectations(snapshot):
    """Canonical records and their declared projections are the only inventory."""
    if snapshot.get("version") != 1 or not snapshot.get("records"):
        raise ValueError("Missing/empty canonical registry snapshot")
    paths, errors, native, ids = {}, [], [], set()
    for group in ("records", "conveniences", "adapters", "properties", "collections", "direct_errors"):
        for record in snapshot[group]:
            identity = record.get("id")
            if identity:
                if identity in ids:
                    errors.append(f"Duplicate canonical ID: {identity}")
                ids.add(identity)
            if record.get("native"):
                native.append({"id": identity, "projection": record["native"]})
                continue  # Explicit native projections do not promise wrapper classes.
            path = record.get("path")
            if record.get("error") or not path:
                errors.append(record.get("error") or f"Missing Python path: {identity}")
                continue
            previous = paths.get(path)
            if previous and previous["kind"] != record["kind"]:
                errors.append(f"Conflicting Python projection kinds: {path}")
                continue
            merged = {**(previous or {}), **record}
            merged["writable"] = record.get("writable", False) or (previous or {}).get("writable", False)
            paths[path] = merged
    return {"paths": paths, "metadata_errors": errors, "native_projections": native,
            "canonical_records": len(snapshot["records"])}


def stub_surface(source):
    nodes = {}
    for node in ast.parse(source).body:
        if isinstance(node, ast.ClassDef):
            nodes[node.name] = node
            for member in node.body:
                if isinstance(member, ast.FunctionDef):
                    nodes.setdefault(f"{node.name}.{member.name}", []).append(member)
        elif isinstance(node, ast.FunctionDef):
            nodes.setdefault(node.name, []).append(node)
    return nodes


def lookup(module, path):
    value = module
    for name in path.split("."):
        value = getattr(value, name)
    return value


def has_decorator(node, name):
    return any((isinstance(tag, ast.Name) and tag.id == name)
               or (isinstance(tag, ast.Attribute) and tag.attr == name)
               for tag in node.decorator_list)


def metadata(document):
    document = (document or "").strip()
    return json.loads(document) if document.startswith("{") else {}


def concrete_annotation(annotation):
    return annotation is not None and not any(
        (isinstance(part, ast.Name) and part.id in {"Any", "PyAny", "Incomplete"})
        or (isinstance(part, ast.Attribute) and part.attr in {"Any", "PyAny", "Incomplete"})
        for part in ast.walk(annotation))


def public_annotations(body, owner=""):
    """Read declarations, including native nested classes, not API names."""
    for node in body:
        if not isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AnnAssign)):
            continue
        if isinstance(node, ast.AnnAssign):
            yield owner + ast.unparse(node.target), node.annotation
        else:
            path = owner + node.name
            if isinstance(node, ast.ClassDef):
                yield from public_annotations(node.body, path + ".")
            else:
                arguments = node.args.posonlyargs + node.args.args + node.args.kwonlyargs
                arguments += [arg for arg in (node.args.vararg, node.args.kwarg) if arg]
                for arg in arguments:
                    if arg.annotation is not None:
                        yield f"{path}({arg.arg})", arg.annotation
                if node.returns is not None:
                    yield f"{path} return", node.returns


def companion_classes(annotation, *, payload=False):
    """Concrete classes, unions and explicitly typed native payload containers."""
    if isinstance(annotation, ast.Name):
        if annotation.id in {"Any", "PyAny", "BaseException", "Exception", "object"}:
            raise ValueError("annotation must name concrete classes")
        # A real native payload can be an enum in one Rust variant and a
        # primitive in another. Primitive union branches need no companion;
        # the payload caller still requires a concrete enum/error branch.
        if payload and annotation.id in {"bool", "int", "float", "str", "bytes"}:
            return set()
        return {annotation.id}
    if isinstance(annotation, ast.Constant) and annotation.value is None:
        return set()
    if isinstance(annotation, ast.BinOp) and isinstance(annotation.op, ast.BitOr):
        return companion_classes(annotation.left, payload=payload) | companion_classes(annotation.right, payload=payload)
    if payload and isinstance(annotation, ast.Subscript):
        if ast.unparse(annotation.value) not in {"list", "tuple", "dict", "set", "frozenset"}:
            raise ValueError("payload container must use a concrete built-in type")
        arguments = annotation.slice.elts if isinstance(annotation.slice, ast.Tuple) else [annotation.slice]
        classes = set()
        for argument in arguments:
            if isinstance(argument, ast.Constant) and argument.value is Ellipsis:
                if ast.unparse(annotation.value) != "tuple":
                    raise ValueError("only a tuple payload may use an ellipsis")
                continue
            classes.update(companion_classes(argument, payload=True))
        return classes
    raise ValueError("annotation must use concrete class names or their union")


def enum_payload(node, runtime):
    """Official native enum constants, not an arbitrary new business class."""
    constants = [member for member in node.body if isinstance(member, ast.AnnAssign)]
    if not constants:
        return False
    for member in constants:
        annotation = member.annotation
        if (not isinstance(member.target, ast.Name)
                or not isinstance(annotation, ast.Subscript)
                or not ((isinstance(annotation.value, ast.Name) and annotation.value.id == "Final")
                        or (isinstance(annotation.value, ast.Attribute) and annotation.value.attr == "Final"))
                or not isinstance(annotation.slice, ast.Name) or annotation.slice.id != node.name
                or not isinstance(getattr(runtime, member.target.id, None), runtime)):
            return False
    return all(not isinstance(member, ast.FunctionDef)
               or (member.name.startswith("__") and member.name.endswith("__")
                   and member.name not in {"__new__", "__init__"})
               for member in node.body)


def check_assignment_effect(cls, field):
    """Use construction as the oracle for a different scalar field value.

    Candidates are test inputs, not domain defaults or constraints. Report
    unexercised fields rather than inventing values for enums, callbacks or
    other structured configuration. Never infer writability from fset.
    """
    previous = getattr(cls(), field)
    if type(previous) is bool:
        candidates = [not previous]
    elif type(previous) is int:
        candidates = [previous + 1, previous - 1]
    elif type(previous) is float:
        candidates = [previous + 0.5, 0.5, 1.0]
    elif type(previous) is str:
        candidates = [previous + "_ck_probe", ""]
    else:
        return False
    for candidate in candidates:
        try:
            expected = getattr(cls(**{field: candidate}), field)
        except (TypeError, ValueError, OverflowError):
            continue
        # A constructor which normalizes the probe to its default provides
        # no evidence that assigning that probe must change the field.
        if type(expected) is not type(previous) or expected == previous:
            continue
        instance = cls()
        setattr(instance, field, candidate)
        actual = getattr(instance, field)
        if type(actual) is not type(expected) or actual != expected:
            raise ValueError(
                f"{cls.__name__}.{field}: setter readback {actual!r} differs "
                f"from constructor {expected!r} for input {candidate!r}"
            )
        return True
    return False


def validate(contract, module, stub):
    expected, nodes = dict(contract["paths"]), stub_surface(stub)
    errors = []
    def fail(category, message):
        errors.append({"category": category, "message": message})
    for message in contract["metadata_errors"]:
        fail("registry_metadata", message)
    for path, annotation in public_annotations(ast.parse(stub).body):
        if not concrete_annotation(annotation):
            fail("type_projection", f"Public annotation must be concrete: {path}: {ast.unparse(annotation)}")
    # An error's typed discriminant is a language representation of the
    # registered error, not an independently maintained API list. Its actual
    # official getter annotation identifies the one allowed companion class.
    canonical_paths = set(expected)
    error_queue = [path for path, record in expected.items() if record.get("error_type")]
    visited_errors = set()

    def include_payload_types(annotation, payload_path, *, require_companion):
        try:
            payload_classes = companion_classes(annotation, payload=True)
        except ValueError as error:
            fail("error_projection", f"Invalid payload annotation {payload_path}: {error}")
            return
        if require_companion and not payload_classes:
            fail("error_projection", f"Payload must reference a concrete value/error class: {payload_path}")
        for payload in sorted(payload_classes):
            payload_class = nodes.get(payload)
            payload_runtime = getattr(module, payload, None)
            if not isinstance(payload_class, ast.ClassDef) or not inspect.isclass(payload_runtime):
                fail("error_projection", f"Concrete payload class missing: {payload_path} -> {payload}")
                continue
            previous = expected.get(payload)
            if previous and previous["kind"] != "class":
                fail("error_projection", f"Payload class conflicts with registry: {payload}")
                continue
            is_error = issubclass(payload_runtime, BaseException)
            is_value = metadata(ast.get_docstring(payload_class)).get("role") == "error_value"
            if not previous and not is_error and not is_value and not enum_payload(payload_class, payload_runtime):
                fail("error_projection", f"Unregistered payload is not a concrete readonly value/enum/error: {payload_path} -> {payload}")
                continue
            expected[payload] = {**(previous or {}), "kind": "class"}
            if is_error or is_value:
                expected[payload]["error_type" if is_error else "error_value"] = True
                error_queue.append(payload)

    while error_queue:
        path = error_queue.pop()
        if path in visited_errors:
            continue
        visited_errors.add(path)
        record = expected[path]
        if not isinstance(nodes.get(path), ast.ClassDef):
            continue
        cls = nodes[path]
        try:
            runtime = lookup(module, path)
        except AttributeError:
            fail("error_projection", f"Registered exception class missing at runtime: {path}")
            continue
        if record.get("error_value"):
            if (metadata(runtime.__doc__).get("role") != "error_value"
                    or metadata(ast.get_docstring(cls)).get("role") != "error_value"
                    or issubclass(runtime, BaseException)):
                fail("error_projection", f"Rich payload must be an explicitly marked native value: {path}")
            getters = []
            for member in cls.body:
                if isinstance(member, (ast.Assign, ast.AnnAssign)):
                    fail("error_projection", f"Rich payload must expose readonly getters, not class fields: {path}")
                if not isinstance(member, ast.FunctionDef):
                    continue
                if member.name.startswith("__") and member.name.endswith("__"):
                    if member.name in {"__new__", "__init__"}:
                        fail("error_projection", f"Rich payload must not introduce a constructor: {path}")
                    continue
                getters.append(member)
                field_path = f"{path}.{member.name}"
                arguments = member.args.posonlyargs + member.args.args
                descriptor = inspect.getattr_static(runtime, member.name, None)
                if (len(nodes.get(field_path, [])) != 1 or len(arguments) != 1
                        or member.args.kwonlyargs or member.args.vararg or member.args.kwarg
                        or not has_decorator(member, "property") or has_decorator(member, "setter")
                        or not inspect.isgetsetdescriptor(descriptor)):
                    fail("error_projection", f"Rich payload field must be a native readonly getter: {field_path}")
                    continue
                include_payload_types(member.returns, field_path, require_companion=False)
            if not getters:
                fail("error_projection", f"Rich payload requires actual typed readonly fields: {path}")
            continue
        if not inspect.isclass(runtime) or not issubclass(runtime, BaseException):
            fail("error_projection", f"Registered error is not an exception class: {path}")
        domains = nodes.get(f"{path}.domain", [])
        domain_descriptor = inspect.getattr_static(runtime, "domain", None)
        if (len(domains) != 1 or not has_decorator(domains[0], "property")
                or not isinstance(domains[0].returns, ast.Name) or domains[0].returns.id != "str"
                or domain_descriptor is None or not hasattr(domain_descriptor, "__get__")
                or getattr(domain_descriptor, "fset", None) is not None):
            fail("error_projection", f"Exception requires a typed readonly domain: {path}")
        source_getters = [member for member in cls.body
                          if isinstance(member, ast.FunctionDef)
                          and metadata(ast.get_docstring(member)).get("role") == "error_source"]
        concrete_sources = set()
        if source_getters:
            getter = source_getters[0]
            source_path = f"{path}.{getter.name}"
            descriptor = inspect.getattr_static(runtime, getter.name, None)
            if (len(source_getters) != 1 or not has_decorator(getter, "property")
                    or len(nodes.get(source_path, [])) != 1
                    or has_decorator(getter, "setter") or descriptor is None
                    or not hasattr(descriptor, "__get__")
                    or getattr(descriptor, "fset", None) is not None
                    or metadata(getattr(descriptor, "__doc__", None)).get("role") != "error_source"):
                fail("error_projection", f"Source must be a native typed readonly getter: {source_path}")
            else:
                try:
                    concrete_sources = companion_classes(getter.returns)
                except ValueError as error:
                    fail("error_projection", f"Invalid source annotation {source_path}: {error}")
            for source in sorted(concrete_sources):
                source_class = nodes.get(source)
                source_runtime = getattr(module, source, None)
                if (not isinstance(source_class, ast.ClassDef) or not inspect.isclass(source_runtime)
                        or not issubclass(source_runtime, BaseException)):
                    fail("error_projection", f"Concrete source exception missing: {source_path} -> {source}")
                    continue
                previous = expected.get(source)
                if previous and previous["kind"] != "class":
                    fail("error_projection", f"Source exception conflicts with registry: {source}")
                    continue
                expected[source] = {**(previous or {}), "kind": "class", "error_type": True}
                error_queue.append(source)
        payload_getters = [member for member in cls.body
                           if isinstance(member, ast.FunctionDef)
                           and metadata(ast.get_docstring(member)).get("role") == "error_payload"]
        for getter in payload_getters:
            payload_path = f"{path}.{getter.name}"
            descriptor = inspect.getattr_static(runtime, getter.name, None)
            arguments = getter.args.posonlyargs + getter.args.args
            if (len(nodes.get(payload_path, [])) != 1 or len(arguments) != 1
                    or getter.args.kwonlyargs or getter.args.vararg or getter.args.kwarg
                    or not has_decorator(getter, "property") or has_decorator(getter, "setter")
                    or not inspect.isgetsetdescriptor(descriptor)
                    or metadata(descriptor.__doc__).get("role") != "error_payload"):
                fail("error_projection", f"Payload must be a native typed readonly getter: {payload_path}")
                continue
            include_payload_types(getter.returns, payload_path, require_companion=True)
        role = metadata(ast.get_docstring(cls)).get("role")
        if role == "error_leaf":
            # A concrete source-free opaque Rust leaf has no public payload or
            # discriminator to invent. Only a real marked error edge may add it.
            if path in canonical_paths:
                fail("error_projection", f"Leaf error must be reached through a concrete source/payload: {path}")
            if metadata(runtime.__doc__).get("role") != "error_leaf":
                fail("error_projection", f"Leaf error metadata differs at runtime: {path}")
            methods = [member for member in cls.body if isinstance(member, ast.FunctionDef)]
            display = [member for member in methods if member.name == "__str__"]
            if (source_getters or payload_getters
                    or any(member.name not in {"domain", "__str__"} for member in methods)
                    or any(isinstance(member, (ast.AnnAssign, ast.Assign)) for member in cls.body)
                    or len(display) != 1 or not isinstance(display[0].returns, ast.Name)
                    or display[0].returns.id != "str"):
                fail("error_projection", f"Leaf error may expose only native domain and Display: {path}")
            continue
        if role == "error_record":
            if metadata(runtime.__doc__).get("role") != "error_record":
                fail("error_projection", f"Record error metadata differs at runtime: {path}")
            if any(isinstance(member, ast.FunctionDef) and member.name == "variant" for member in cls.body):
                fail("error_projection", f"Record error must not invent a discriminant: {path}")
            payloads = {name: entry for name, entry in expected.items()
                        if name.rpartition(".")[0] == path and entry.get("read_accessor")}
            native_fields = [member for member in cls.body
                             if isinstance(member, ast.FunctionDef)
                             and has_decorator(member, "property")
                             and member.name not in {"domain", "source_message"}
                             and metadata(ast.get_docstring(member)).get("role") != "error_source"]
            if not payloads and not native_fields and not concrete_sources:
                fail("error_projection", f"Record error requires typed readonly payloads or a concrete source: {path}")
            for member in native_fields:
                name = f"{path}.{member.name}"
                declarations = nodes.get(name, [])
                descriptor = inspect.getattr_static(runtime, member.name, None)
                arguments = member.args.posonlyargs + member.args.args
                if (len(declarations) != 1 or len(arguments) != 1 or member.args.kwonlyargs
                        or member.args.vararg or member.args.kwarg
                        or not concrete_annotation(member.returns)
                        or descriptor is None or not hasattr(descriptor, "__get__")
                        or getattr(descriptor, "fset", None) is not None):
                    fail("error_projection", f"Native record payload must be typed and readonly: {name}")
            for name in payloads:
                declarations = nodes.get(name, [])
                accessor = inspect.getattr_static(runtime, name.rpartition(".")[2], None)
                if len(declarations) != 1 or not callable(accessor):
                    fail("error_projection", f"Record error payload accessor missing: {name}")
                    continue
                member = declarations[0]
                arguments = member.args.posonlyargs + member.args.args
                if (len(arguments) != 1 or member.args.kwonlyargs or member.args.vararg
                        or member.args.kwarg or not concrete_annotation(member.returns)
                        or has_decorator(member, "setter")):
                    fail("error_projection", f"Record error payload must be typed and readonly: {name}")
            continue
        getters = [member for member in nodes[path].body
                   if isinstance(member, ast.FunctionDef) and member.name == "variant"
                   and has_decorator(member, "property")]
        if len(getters) != 1 or not isinstance(getters[0].returns, ast.Name):
            fail("error_projection", f"Registered error requires a typed variant getter: {path}")
            continue
        companion = getters[0].returns.id
        if companion == path or not isinstance(nodes.get(companion), ast.ClassDef):
            fail("error_projection", f"Missing discriminant class for {path}.variant: {companion}")
            continue
        existing = expected.get(companion)
        if existing and existing["kind"] != "class":
            fail("error_projection", f"Discriminant projection conflicts with registry: {companion}")
            continue
        expected.setdefault(companion, {"kind": "class"})
    matched_runtime = matched_stub = 0
    configurations = fields_checked = 0
    required_configurations = []
    assignment_effect_checked = 0
    assignment_effect_unexercised = []
    for path, record in expected.items():
        kind = record["kind"]
        value = None
        try:
            # Do not execute a property getter on an unconstructed object.
            owner_path, _, name = path.rpartition(".")
            owner = lookup(module, owner_path) if owner_path else module
            value = inspect.getattr_static(owner, name) if kind == "property" else getattr(owner, name)
        except AttributeError:
            fail("missing_runtime", f"Missing runtime {kind}: {path}")
        else:
            valid = inspect.isclass(value) if kind == "class" else (
                hasattr(value, "__get__") if kind == "property" else callable(value))
            if valid:
                matched_runtime += 1
            else:
                fail("runtime_kind", f"Wrong runtime kind: {path} (expected {kind})")
            # __set__ is necessary but not proof of effective writability;
            # native readonly descriptors can expose it too. Check actual
            # configuration writes and their effects below.
            if kind == "property" and record.get("writable") and not hasattr(value, "__set__"):
                fail("configuration", f"Runtime property not writable: {path}")
        node = nodes.get(path)
        valid_stub = isinstance(node, ast.ClassDef) if kind == "class" else isinstance(node, list)
        if kind == "property" and valid_stub:
            valid_stub = any(has_decorator(member, "property") for member in node)
            if record.get("writable") and not any(has_decorator(member, "setter") for member in node):
                fail("configuration", f"Stub property not writable: {path}")
        if valid_stub:
            matched_stub += 1
        else:
            fail("missing_stub", f"Missing stub {kind}: {path}")
        fields = record.get("configuration")
        if fields is None:
            continue
        configurations += 1
        fields_checked += len(fields)
        if inspect.isclass(value):
            try:
                parameters = inspect.signature(value).parameters
                schema = record.get("configuration_schema")
                if schema is None or {field["name"] for field in schema} != set(fields):
                    fail("registry_metadata", f"Missing configuration schema: {path}")
                    continue
                aliases = {alias for field in schema for alias in field["aliases"]}
                if aliases & set(fields) or len(aliases) != sum(len(field["aliases"]) for field in schema):
                    fail("registry_metadata", f"Conflicting configuration aliases: {path}")
                if set(parameters) != set(fields) | aliases:
                    fail("configuration", f"Constructor fields/aliases differ from registry: {path}")
                required = {field["name"] for field in schema if field["required"]}
                actual_required = {name for name, parameter in parameters.items()
                                   if parameter.default is inspect.Parameter.empty}
                if actual_required != required:
                    fail("configuration", f"Required constructor inputs differ from registry: {path}")
                if required:
                    # Required domain inputs have no generic valid default.
                    # Check their declarations without inventing probe values.
                    required_configurations.append(path)
                    assignment_effect_unexercised.extend(f"{path}.{field}" for field in fields)
                    assignment_effect_unexercised.extend(f"{path}.{alias}" for alias in sorted(aliases))
                else:
                    default = value()
                    for field in [*fields, *sorted(aliases)]:
                        setattr(default, field, getattr(default, field))
                        if check_assignment_effect(value, field):
                            assignment_effect_checked += 1
                        else:
                            assignment_effect_unexercised.append(f"{path}.{field}")
            except (AttributeError, TypeError, ValueError, OverflowError) as error:
                fail("configuration", f"Configuration field unavailable/not writable: {path}: {error}")
        if isinstance(node, ast.ClassDef):
            setters = {member.name for member in node.body if isinstance(member, ast.FunctionDef) and has_decorator(member, "setter")}
            aliases = {alias for field in record.get("configuration_schema", []) for alias in field["aliases"]}
            if setters != set(fields) | aliases:
                fail("configuration", f"Stub writable fields differ from registry: {path}")
    roots = {path.split(".")[0] for path in expected}
    for name in sorted(set(getattr(module, "__all__", [])) - roots):
        fail("unregistered_export", f"Unregistered public export: {name}")
    report = {"canonical_records": contract["canonical_records"],
              "native_projections": contract["native_projections"],
              "expected_paths": len(expected), "matched_runtime": matched_runtime,
              "matched_stub": matched_stub, "configuration_records": configurations,
              "configuration_fields": fields_checked,
              "configuration_assignment_effect_checked": assignment_effect_checked,
              "configuration_assignment_effect_not_exercised": assignment_effect_unexercised,
              "required_configuration_assignment_not_exercised": required_configurations,
              "counts": dict(Counter(item["category"] for item in errors)), "errors": errors}
    if errors:
        raise RegistryMismatch(report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package", type=Path, default=HERE / "artifacts/smoke/package",
                        help="directory containing the already generated cosmolkit package")
    parser.add_argument("--snapshot-bin", type=Path, default=HERE.parent / "target/typed_facade_projection/debug/registry_snapshot")
    parser.add_argument("--report", type=Path, default=HERE / "artifacts/registry.json")
    options = parser.parse_args()
    if not options.snapshot_bin.is_file():
        raise ValueError("Registry reader is not built; run build.py first. This checker never compiles.")
    result = subprocess.run([str(options.snapshot_bin)], text=True, capture_output=True, check=True)
    contract = expectations(json.loads(result.stdout))
    sys.path.insert(0, str(options.package.resolve()))
    module = importlib.import_module("cosmolkit")
    if not Path(module.__file__).resolve().is_relative_to(options.package.resolve()):
        raise ValueError("Imported a different package instead of the requested generated package")
    stub = options.package / "cosmolkit/__init__.pyi"
    try:
        report = validate(contract, module, stub.read_text())
    except RegistryMismatch as error:
        report = error.report
    options.report.parent.mkdir(parents=True, exist_ok=True)
    options.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: value for key, value in report.items() if key not in {"errors", "native_projections"}}))
    for error in report["errors"][:10]:
        print(error["message"], file=sys.stderr)
    print(f"Full report: {options.report}")
    return int(bool(report["errors"]))


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, ImportError, OSError, subprocess.CalledProcessError) as error:
        print(f"registry check failed: {error}", file=sys.stderr)
        raise SystemExit(1)
