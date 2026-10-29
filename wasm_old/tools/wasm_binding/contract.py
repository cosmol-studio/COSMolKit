"""Check distribution projections against the actual compiled Rust registry."""

from __future__ import annotations

import importlib.util
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
_spec = importlib.util.spec_from_file_location("_python_projection_contract", ROOT / "dev/tools/check_python_stub_contract.py")
assert _spec is not None and _spec.loader is not None
_shared = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_shared)


def compiled_contract(api_root: Path, active: set[str], destination: Path, profile: str = "dev-test") -> dict:
    command = [
        "cargo", "run", "--quiet", "--locked", "--profile", profile,
        "-p", "cosmolkit", "--example", "binding_contract_manifest",
        "--no-default-features", "--features", ",".join(sorted(active)),
    ]
    result = subprocess.run(command, cwd=api_root, check=True, text=True, stdout=subprocess.PIPE)
    document = json.loads(result.stdout)
    # The only platform exclusion comes from the canonical export, not names
    # or a distribution-local list of missing methods.
    document["entries"] = [row for row in document["entries"] if row["platform"] != "native"]
    document["javascript_aliases"] = [row for row in document["javascript_aliases"] if row["platform"] != "native"]
    configuration_fields(document)
    destination.write_text(json.dumps(document), encoding="utf-8")
    return document


def configuration_fields(document: dict) -> None:
    """Use the registered configuration schema; Rust need not provide new().

    The language constructor schema is already declared in the canonical
    registry for public-field/Default Rust records. This only supplies the
    fields to check: it does not synthesize a JavaScript constructor, property
    or declaration. The actual JS/TS gates must still verify every field.
    """
    for row in document["entries"]:
        if row.get("role") == "parameter" and row.get("python_fields") is not None:
            row["fields"] = row["python_fields"]


def camel(name: str) -> str:
    return re.sub(r"_([a-z0-9])", lambda match: match[1].upper(), name)


def generate_conveniences(package: Path, library_name: str, document: dict) -> None:
    """Export the actual static factories; retain their exact TypeScript types."""
    javascript, declarations = [], []
    for row in document["javascript_aliases"]:
        name, target = camel(row["name"]), row["target"]
        comment = f"/** Alias of `{target}`. */\n"
        javascript.append(f"{comment}export const {name} = {target};\n")
        declarations.append(f"{comment}export const {name}: typeof {target};\n")
    for suffix, lines in [("js", javascript), ("d.ts", declarations)]:
        with (package / f"{library_name}.{suffix}").open("a", encoding="utf-8") as output:
            output.write("\n" + "".join(lines))


def convenience_typescript_probe(document: dict) -> str:
    lines = []
    for index, row in enumerate(document.get("javascript_aliases", [])):
        name, target = camel(row["name"]), row["target"]
        # Check both directions: no lost overload, parameter, default or result.
        lines.append(f"const alias{index}: typeof ck.{target} = ck.{name};")
        lines.append(f"const original{index}: typeof ck.{name} = ck.{target};")
    return "\n".join(lines)


def typescript_probe(document: dict) -> str:
    """Let TypeScript check its own declarations instead of regex parsing TS."""
    entries = document["entries"]
    types = {row["semantic_id"].rsplit(".", 1)[-1]: row["javascript_name"] for row in entries if row["item"] == "type"}
    lines = ['import * as ck from "cosmolkit-generated";']
    for index, row in enumerate(entries):
        name = row["javascript_name"]
        if row["item"] == "type" or row["owner"] == "module":
            lines.append(f"void ck.{name}; // {row['semantic_id']}")
        else:
            owner = "Molecule" if row["owner"] == "molecule" else types.get(row["semantic_id"].rsplit(".", 1)[0], row["semantic_id"].rsplit(".", 1)[0])
            if name != "new":
                holder = f"ck.{owner}.prototype" if row.get("receiver") else f"ck.{owner}"
                lines.append(f"void {holder}.{name}; // {row['semantic_id']}")
        if row["item"] == "type":
            for field in row.get("properties", []):
                lines.append(f"void ck.{name}.prototype.{camel(field['name'])};")
        if row.get("role") != "parameter":
            continue
        if row.get("fields") is None:
            raise ValueError(f"{row['semantic_id']}: Parameter needs a registered canonical constructor")
        lines.append(f"declare const p{index}: ck.{name};")
        for field in row["fields"]:
            key = camel(field["name"])
            # An assignment to itself proves declared writability without
            # inventing values or relying on TS's assignability of readonly
            # object types. The runtime gate checks the actual descriptor.
            lines.append(f"p{index}.{key} = p{index}.{key};")
    for index, (base, target, configurations) in enumerate(_shared.configuration_calls(entries, document.get("keywords", ()))):
        owner = "Molecule" if base["owner"] == "molecule" else types.get(base["semantic_id"].rsplit(".", 1)[0], base["semantic_id"].rsplit(".", 1)[0])
        if base["owner"] == "module":
            call = f"ck.{base['javascript_name']}"
        elif not base.get("receiver"):
            call = f"ck.{owner}.{base['javascript_name']}"
        else:
            lines.append(f"declare const receiver{index}: ck.{owner};")
            call = f"receiver{index}.{base['javascript_name']}"
        configuration_names = {field["name"] for field, _ in configurations}
        data_names = {field["name"] for field in base.get("parameters") or []} - configuration_names
        arguments = [f"null as unknown as Parameters<typeof {call}>[{position}]" for position, field in enumerate(base.get("parameters") or []) if field["name"] not in configuration_names]
        prefix = ", ".join(arguments)
        prefix += ", " if prefix else ""
        instances, objects = [], []
        for position, (parameter, configuration) in enumerate(configurations):
            variable = f"configuration{index}_{position}"
            lines.append(f"declare const {variable}: ck.{configuration['javascript_name']};")
            instances.append(variable)
            options = ", ".join(f"{camel(field['name'])}: {variable}.{camel(field['name'])}" for field in configuration["fields"])
            objects.append("{" + options + "}")
        lines.append(f"{call}({prefix}{', '.join(instances)});")
        lines.append(f"{call}({prefix}{', '.join(objects)});")
        for position, (parameter, configuration) in enumerate(configurations):
            for required in configuration["fields"]:
                # A convenience call may take this field as its explicit data
                # argument (e.g. fragment atoms), then fill the configuration
                # itself. That call is not missing the required input.
                if required["default"] is not None or required["name"] in data_names:
                    continue
                variable = f"configuration{index}_{position}"
                incomplete = ", ".join(f"{camel(field['name'])}: {variable}.{camel(field['name'])}" for field in configuration["fields"] if field["name"] != required["name"])
                rejected = list(objects)
                rejected[position] = "{" + incomplete + "}"
                lines.append(f"// @ts-expect-error required configuration field {camel(required['name'])}")
                lines.append(f"{call}({prefix}{', '.join(rejected)});")
    lines.append(convenience_typescript_probe(document))
    return "\n".join(lines) + "\n"
