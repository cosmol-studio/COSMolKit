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


def check_contract(stub: str, contract_json: str) -> list[str]:
    """Return every missing callable; malformed input fails rather than skips."""
    entries = cast(list[ContractEntry], json.loads(contract_json))
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
                    accessors[member.target.id] = {"getter", "setter"}
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
    return missing
