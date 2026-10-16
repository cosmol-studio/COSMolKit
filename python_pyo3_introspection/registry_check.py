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
    for group in ("records", "conveniences", "adapters", "properties", "collections"):
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


def validate(contract, module, stub):
    expected, nodes = contract["paths"], stub_surface(stub)
    errors = []
    def fail(category, message):
        errors.append({"category": category, "message": message})
    for message in contract["metadata_errors"]:
        fail("registry_metadata", message)
    matched_runtime = matched_stub = 0
    configurations = fields_checked = 0
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
            if kind == "property" and record.get("writable") and getattr(value, "fset", None) is None:
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
                if set(inspect.signature(value).parameters) != set(fields):
                    fail("configuration", f"Constructor fields differ from registry: {path}")
                default = value()
                for field in fields:
                    setattr(default, field, getattr(default, field))
            except (AttributeError, TypeError, ValueError) as error:
                fail("configuration", f"Configuration field unavailable/not writable: {path}: {error}")
        if isinstance(node, ast.ClassDef):
            setters = {member.name for member in node.body if isinstance(member, ast.FunctionDef) and has_decorator(member, "setter")}
            if setters != set(fields):
                fail("configuration", f"Stub writable fields differ from registry: {path}")
    roots = {path.split(".")[0] for path in expected}
    for name in sorted(set(getattr(module, "__all__", [])) - roots):
        fail("unregistered_export", f"Unregistered public export: {name}")
    report = {"canonical_records": contract["canonical_records"],
              "native_projections": contract["native_projections"],
              "expected_paths": len(expected), "matched_runtime": matched_runtime,
              "matched_stub": matched_stub, "configuration_records": configurations,
              "configuration_fields": fields_checked,
              "counts": dict(Counter(item["category"] for item in errors)), "errors": errors}
    if errors:
        raise RegistryMismatch(report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package", type=Path, default=HERE / "artifacts/smoke/package",
                        help="directory containing the already generated fixed_projection package")
    parser.add_argument("--snapshot-bin", type=Path, default=HERE.parent / "target/typed_facade_projection/debug/registry_snapshot")
    parser.add_argument("--report", type=Path, default=HERE / "artifacts/registry.json")
    options = parser.parse_args()
    if not options.snapshot_bin.is_file():
        raise ValueError("Registry reader is not built; run check.py first. This checker never compiles.")
    result = subprocess.run([str(options.snapshot_bin)], text=True, capture_output=True, check=True)
    contract = expectations(json.loads(result.stdout))
    sys.path.insert(0, str(options.package.resolve()))
    module = importlib.import_module("fixed_projection")
    if not Path(module.__file__).resolve().is_relative_to(options.package.resolve()):
        raise ValueError("Imported a different package instead of the requested generated package")
    stub = options.package / "fixed_projection/__init__.pyi"
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
