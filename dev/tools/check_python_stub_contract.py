"""Check generated Python declarations against the compiled Rust registry.

The generator supplies BINDING_CONTRACT from its own linked cosmolkit build.
This module neither parses registry source nor maintains a second API list.
"""

from __future__ import annotations

import ast
import json
from typing import Literal, TypedDict, cast


class ContractEntry(TypedDict):
    semantic_id: str
    python_name: str
    feature: str
    item: Literal["callable", "type"]
    owner: Literal["module", "molecule", "type"]


def check_contract(stub: str, contract_json: str) -> list[str]:
    """Return every missing callable; malformed input fails rather than skips."""
    entries = cast(list[ContractEntry], json.loads(contract_json))
    module = ast.parse(stub, filename="cosmolkit.pyi")
    functions: set[str] = set()
    classes: dict[str, set[str]] = {}
    for node in module.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            functions.add(node.name)
        elif isinstance(node, ast.ClassDef):
            classes[node.name] = {
                member.name
                for member in node.body
                if isinstance(member, (ast.FunctionDef, ast.AsyncFunctionDef))
                and not any(
                    isinstance(decorator, ast.Name) and decorator.id == "property"
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
            found = name in classes.get(owner, set())
        if not found:
            missing.append(
                f"{entry['semantic_id']} -> cosmolkit.{path} "
                + f"[feature={entry['feature']}]"
            )
    return missing
