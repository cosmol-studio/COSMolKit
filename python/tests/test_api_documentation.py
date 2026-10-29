"""Help text is part of the actual compiled Python API, not a separate inventory."""

import inspect
import json
from collections.abc import Mapping
from typing import cast

import cosmolkit as ck
import pytest


def has_description(documentation: str | None) -> bool:
    text = (documentation or "").strip()
    if text.startswith("{"):
        try:
            metadata, end = json.JSONDecoder().raw_decode(text)
        except ValueError:
            pass
        else:
            if isinstance(metadata, dict) and "role" in metadata:
                text = text[end:].strip()
    return bool(text)


def test_public_api_has_native_descriptions() -> None:
    missing: list[str] = []
    checked = 0
    visited: set[type[object]] = set()

    def check(name: str, value: object) -> None:
        nonlocal checked
        if inspect.isclass(value) and (value.__module__ == "builtins" or value.__module__.split(".")[0] == "cosmolkit"):
            if value in visited:
                return
            visited.add(value)
            checked += 1
            # Inherited str/ValueError help is not a description of a CK type.
            if not has_description(value.__doc__):
                missing.append(name)
            for member_name, member in cast(Mapping[str, object], vars(value)).items():
                if member_name.startswith("_"):
                    continue
                if isinstance(member, (staticmethod, classmethod)):
                    member = member.__func__
                check(f"{name}.{member_name}", member)
        elif callable(value) or inspect.isdatadescriptor(value):
            checked += 1
            if not has_description(inspect.getdoc(value)):
                missing.append(name)

    for name, value in cast(Mapping[str, object], vars(ck)).items():
        if not name.startswith("_"):
            check(name, value)
    assert checked > 0
    assert not missing, "Missing public API descriptions: " + ", ".join(missing)


@pytest.mark.parametrize(
    ("owner", "name", "details"),
    [
        (ck.Molecule, "coordinates_2d", ("float64", "(num_atoms, 2)", "None")),
        (ck.Conformer3D, "coordinates", ("float64", "(num_atoms, 3)", "independent")),
        (ck.Fingerprint, "to_numpy", ("uint8", "bit")),
        (ck.FingerprintBatch, "to_numpy", ("uint8", "None", "unequal widths")),
        (ck.Molecule, "atom", ("zero-based", "None")),
        (ck.Molecule, "with_hydrogens", ("explicit hydrogen", "source molecule is unchanged")),
        (ck.Molecule, "add_hydrogens_", ("In place", "Copy-on-write")),
        (ck.SubstructMatchParams, "final_match", ("target, atom_indices", "exceptions propagate")),
        (ck.SubstructMatchParams, "bond_match", ("query_bond, target_bond", "Bond values")),
        (ck.BatchValidationError, "record_errors", ("BatchError", "copy")),
        (ck.BioNcsOperator, "id", ("Identifier", "symmetry operator")),
    ],
)
def test_help_describes_observable_contract(
    owner: type[object], name: str, details: tuple[str, ...]
) -> None:
    documentation = inspect.getdoc(cast(object, getattr(owner, name)))
    assert documentation is not None
    for detail in details:
        assert detail in documentation


@pytest.mark.parametrize(
    ("field", "alias"),
    [
        ("final_match", "extra_final_check"),
        ("atom_match", "extra_atom_check"),
        ("bond_match", "extra_bond_check"),
    ],
)
def test_callback_alias_keeps_the_native_field_description(field: str, alias: str) -> None:
    assert inspect.getdoc(cast(object, getattr(ck.SubstructMatchParams, alias))) == inspect.getdoc(
        cast(object, getattr(ck.SubstructMatchParams, field))
    )
