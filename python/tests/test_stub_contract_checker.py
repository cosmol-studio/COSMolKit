"""Fixed regressions for the registry-driven stub-generation gate."""

from __future__ import annotations

import importlib.util
import json
from collections.abc import Callable
from pathlib import Path
from typing import cast

import pytest


_path = Path(__file__).resolve().parents[2] / "dev/tools/check_python_stub_contract.py"
_spec = importlib.util.spec_from_file_location("check_python_stub_contract", _path)
assert _spec is not None and _spec.loader is not None
_checker = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_checker)
check_contract = cast(Callable[[str, str], list[str]], _checker.check_contract)


def _entry(semantic_id: str, name: str, owner: str = "molecule", item: str = "callable"):
    return {
        "semantic_id": semantic_id,
        "python_name": name,
        "owner": owner,
        "item": item,
        "feature": "cap-valence",
        "status": "experimental",
    }


def test_module_class_factory_and_inplace_declarations_are_required() -> None:
    entries = [
        _entry("search.parse_smarts", "parse_smarts", "module"),
        _entry("types.QueryGraph", "QueryGraph", "type", "type"),
        _entry("QueryGraph.from_smarts", "from_smarts", "type"),
        _entry("Molecule.assign_valence_", "assign_valence_"),
    ]
    stub = """
def parse_smarts(text): ...
class QueryGraph:
    @staticmethod
    def from_smarts(text): ...
class Molecule:
    def assign_valence_(self): ...
"""
    assert check_contract(stub, json.dumps(entries)) == []


def test_runtime_and_stub_both_omitting_inplace_methods_cannot_pass() -> None:
    entries = [
        _entry("Molecule.assign_valence_", "assign_valence_"),
        _entry("Molecule.assign_valence_with_params_", "assign_valence_with_params_"),
    ]
    assert check_contract(
        "class Molecule:\n    def with_assigned_valence(self): ...\n", json.dumps(entries)
    ) == [
        "Molecule.assign_valence_ -> cosmolkit.Molecule.assign_valence_ [feature=cap-valence]",
        "Molecule.assign_valence_with_params_ -> cosmolkit.Molecule.assign_valence_with_params_ [feature=cap-valence]",
    ]


def test_comments_strings_attributes_and_wrong_owners_do_not_count() -> None:
    entries = [_entry("Molecule.assign_valence_", "assign_valence_")]
    for stub in [
        "# def assign_valence_(self): ...\n",
        'text = "def assign_valence_(self): ..."\n',
        "class Molecule:\n    assign_valence_: object\n",
        "def assign_valence_(self): ...\n",
        "class Other:\n    def assign_valence_(self): ...\n",
        "class Molecule:\n    @property\n    def assign_valence_(self): ...\n",
        "class Molecule:\n    @property\n    def assign_valence_(self): ...\n    @assign_valence_.setter\n    def assign_valence_(self, value): ...\n",
    ]:
        assert len(check_contract(stub, json.dumps(entries))) == 1


def test_only_compiled_rows_are_required_without_status_exemptions() -> None:
    # The generator supplies only cfg-enabled rows. No second feature list here.
    entries = [_entry("Molecule.assign_valence_", "assign_valence_")]
    assert check_contract(
        "class Molecule:\n    def assign_valence_(self): ...\n", json.dumps(entries)
    ) == []
    assert len(check_contract("", json.dumps(entries))) == 1


def test_registered_type_projection_and_exact_constructor_name() -> None:
    entries = [
        _entry("types.RustParams", "PythonParams", "type", "type"),
        _entry("RustParams.new", "new", "type"),
    ]
    assert check_contract(
        "class PythonParams:\n    @staticmethod\n    def new(): ...\n", json.dumps(entries)
    ) == []
    # Do not invent an unregistered new -> __new__ compatibility mapping.
    assert len(check_contract(
        "class PythonParams:\n    def __new__(cls): ...\n", json.dumps(entries)
    )) == 1


def test_invalid_stub_and_registry_fail_explicitly() -> None:
    with pytest.raises(SyntaxError):
        _ = check_contract("class:", "[]")
    with pytest.raises(json.JSONDecodeError):
        _ = check_contract("", "not json")


def test_explicit_constructor_projection() -> None:
    entries = [_entry("Params.new", "__new__", "type")]
    assert check_contract("class Params:\n    def __new__(cls): ...\n", json.dumps(entries)) == []
    assert len(check_contract("class Params:\n    def new(cls): ...\n", json.dumps(entries))) == 1


def test_explicit_property_access_requires_descriptor_not_method() -> None:
    getter = {**_entry("Params.limit", "limit", "type"), "python_property": "getter"}
    setter = {**_entry("Params.set_limit", "limit", "type"), "python_property": "setter"}
    readonly = "class Params:\n    @property\n    def limit(self) -> int: ...\n"
    writable = readonly + "    @limit.setter\n    def limit(self, value: int) -> None: ...\n"
    assert check_contract(readonly, json.dumps([getter])) == []
    assert len(check_contract(readonly, json.dumps([setter]))) == 1
    assert check_contract(writable, json.dumps([getter, setter])) == []
    assert check_contract("class Params:\n    limit: int\n", json.dumps([getter, setter])) == []
    assert len(check_contract("class Params:\n    def limit(self): ...\n", json.dumps([getter, setter]))) == 2
