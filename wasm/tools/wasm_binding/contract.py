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


def compiled_contract(api_root: Path, active: set[str], destination: Path) -> dict:
    command = [
        "cargo", "run", "--quiet", "--locked", "--profile", "dev-test",
        "-p", "cosmolkit", "--example", "binding_contract_manifest",
        "--no-default-features", "--features", ",".join(sorted(active)),
    ]
    result = subprocess.run(command, cwd=api_root, check=True, text=True, stdout=subprocess.PIPE)
    document = json.loads(result.stdout)
    # The only platform exclusion comes from the canonical export, not names
    # or a distribution-local list of missing methods.
    document["entries"] = [row for row in document["entries"] if row["platform"] != "native"]
    destination.write_text(json.dumps(document), encoding="utf-8")
    return document


def camel(name: str) -> str:
    return re.sub(r"_([a-z])", lambda match: match[1].upper(), name)


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
    for index, (base, target, configurations) in enumerate(_shared.configuration_calls(entries)):
        owner = "Molecule" if base["owner"] == "molecule" else types.get(base["semantic_id"].rsplit(".", 1)[0], base["semantic_id"].rsplit(".", 1)[0])
        if len(configurations) != 1:
            raise ValueError(f"{base['semantic_id']}: multiple configuration objects need an explicit canonical language projection")
        parameter, configuration = configurations[0]
        parameter_type = configuration["javascript_name"]
        if configuration["fields"] is None:
            continue  # The type check above already rejects this incomplete schema.
        if base["owner"] == "module":
            call = f"ck.{base['javascript_name']}"
        else:
            lines.append(f"declare const receiver{index}: ck.{owner};")
            call = f"receiver{index}.{base['javascript_name']}"
        arguments = [f"null as unknown as Parameters<typeof {call}>[{position}]" for position, _ in enumerate(base.get("parameters") or [])]
        prefix = ", ".join(arguments)
        prefix += ", " if prefix else ""
        lines.append(f"declare const configuration{index}: ck.{parameter_type};")
        lines.append(f"{call}({prefix}configuration{index});")
        options = ", ".join(f"{camel(field['name'])}: configuration{index}.{camel(field['name'])}" for field in configuration["fields"])
        lines.append(f"{call}({prefix}{{{options}}});")
    return "\n".join(lines) + "\n"
