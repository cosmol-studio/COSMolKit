#!/usr/bin/env python3
"""Generate complete modern RDKit CIPLabeler observable-state records."""

from __future__ import annotations

import argparse
import json
import multiprocessing
from importlib.metadata import version
from pathlib import Path
from typing import Any, Iterable

from rdkit import Chem, RDLogger
from rdkit.Chem import rdCIPLabeler


EXPECTED_RDKIT_VERSION = "2026.3.1"
UINT32_MAX = (1 << 32) - 1
REPO_ROOT = Path(__file__).resolve().parents[3]


def fixture_reference(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return resolved.as_posix()


def assert_rdkit_version() -> None:
    actual = version("rdkit")
    if actual != EXPECTED_RDKIT_VERSION:
        raise SystemExit(
            f"RDKit version mismatch: expected {EXPECTED_RDKIT_VERSION}, got {actual}"
        )


def property_names(owner: Any, include_computed: bool) -> set[str]:
    return set(
        owner.GetPropNames(
            includePrivate=True,
            includeComputed=include_computed,
        )
    )


def property_state(owner: Any, name: str, value_kind: str) -> dict[str, Any]:
    persistent = property_names(owner, include_computed=False)
    complete = property_names(owner, include_computed=True)
    if name not in complete:
        return {
            "present": False,
            "value": None,
            "value_type": None,
            "computed": False,
        }

    if value_kind == "boolean":
        value: Any = owner.GetBoolProp(name)
    elif value_kind == "unsigned":
        value = owner.GetUnsignedProp(name)
    elif value_kind == "unsigned_vector":
        value = [int(item) & UINT32_MAX for item in json.loads(owner.GetProp(name))]
    elif value_kind == "string":
        value = owner.GetProp(name)
    else:
        raise AssertionError(f"unsupported property kind: {value_kind}")
    return {
        "present": True,
        "value": value,
        "value_type": value_kind,
        "computed": name not in persistent,
    }


def snapshot(mol: Chem.Mol) -> dict[str, Any]:
    atoms = []
    for atom in mol.GetAtoms():
        atoms.append(
            {
                "index": atom.GetIdx(),
                "chiral_tag": str(atom.GetChiralTag()),
                "cip_code": property_state(atom, "_CIPCode", "string"),
                "cip_neighbor_order": property_state(
                    atom, "_CIPNeighborOrder", "unsigned_vector"
                ),
                "cip_rank": property_state(atom, "_CIPRank", "unsigned"),
            }
        )
    bonds = []
    for bond in mol.GetBonds():
        bonds.append(
            {
                "index": bond.GetIdx(),
                "begin": bond.GetBeginAtomIdx(),
                "end": bond.GetEndAtomIdx(),
                "stereo": str(bond.GetStereo()),
                "stereo_atoms_u32": [
                    int(item) & UINT32_MAX for item in bond.GetStereoAtoms()
                ],
                "cip_code": property_state(bond, "_CIPCode", "string"),
                "cip_neighbor_order": property_state(
                    bond, "_CIPNeighborOrder", "unsigned_vector"
                ),
            }
        )
    return {
        "molecule": {
            "cip_computed": property_state(
                mol, "_CIPComputed", "boolean"
            )
        },
        "atoms": atoms,
        "bonds": bonds,
    }


def set_property(owner: Any, name: str, value: Any, computed: bool = False) -> None:
    if isinstance(value, bool):
        owner.SetBoolProp(name, value, computed)
    elif isinstance(value, int):
        owner.SetUnsignedProp(name, value, computed)
    elif isinstance(value, list):
        owner.SetProp(name, json.dumps(value, separators=(",", ":")), computed)
    else:
        owner.SetProp(name, str(value), computed)


def apply_initial_state(mol: Chem.Mol, config: dict[str, Any]) -> None:
    for name in config.get("clear_molecule_props", []):
        if mol.HasProp(name):
            mol.ClearProp(name)
    for item in config.get("molecule_props", []):
        set_property(mol, item["name"], item["value"], item.get("computed", False))
    for item in config.get("atom_props", []):
        atom = mol.GetAtomWithIdx(int(item["index"]))
        if item.get("clear", False):
            if atom.HasProp(item["name"]):
                atom.ClearProp(item["name"])
        else:
            set_property(
                atom, item["name"], item["value"], item.get("computed", False)
            )
    for item in config.get("bond_props", []):
        bond = mol.GetBondWithIdx(int(item["index"]))
        if item.get("clear", False):
            if bond.HasProp(item["name"]):
                bond.ClearProp(item["name"])
        else:
            set_property(
                bond, item["name"], item["value"], item.get("computed", False)
            )


def parse_case(case: dict[str, Any]) -> Chem.Mol | None:
    input_kind = case.get("input_kind", "smiles")
    source = case["input"]
    sanitize = bool(case.get("sanitize", True))
    if input_kind == "smiles":
        params = Chem.SmilesParserParams()
        params.sanitize = sanitize
        params.removeHs = bool(case.get("remove_hs", True))
        params.allowCXSMILES = True
        params.parseName = False
        return Chem.MolFromSmiles(source, params)
    if input_kind == "molblock":
        return Chem.MolFromMolBlock(
            source,
            sanitize=sanitize,
            removeHs=bool(case.get("remove_hs", True)),
            strictParsing=bool(case.get("strict_parsing", True)),
        )
    raise ValueError(f"unsupported input_kind: {input_kind!r}")


def default_calls(mol: Chem.Mol) -> list[dict[str, Any]]:
    atom_selection = next(
        (
            [atom.GetIdx()]
            for atom in mol.GetAtoms()
            if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
        ),
        [],
    )
    bond_selection = next(
        (
            [bond.GetIdx()]
            for bond in mol.GetBonds()
            if bond.GetStereo() != Chem.BondStereo.STEREONONE
        ),
        [],
    )
    return [
        {
            "atoms_to_label": None,
            "bonds_to_label": None,
            "max_recursive_iterations": 0,
        },
        {
            "atoms_to_label": atom_selection,
            "bonds_to_label": None,
            "max_recursive_iterations": 0,
        },
        {
            "atoms_to_label": None,
            "bonds_to_label": bond_selection,
            "max_recursive_iterations": 0,
        },
        {
            "atoms_to_label": [],
            "bonds_to_label": [],
            "max_recursive_iterations": 1,
        },
    ]


def execute_call(mol: Chem.Mol, call: dict[str, Any]) -> dict[str, Any]:
    atoms = call.get("atoms_to_label")
    bonds = call.get("bonds_to_label")
    limit = int(call.get("max_recursive_iterations", 0))
    try:
        rdCIPLabeler.AssignCIPLabels(mol, atoms, bonds, limit)
        result = {"status": "ok", "error_type": None, "error_message": None}
    except Exception as error:  # noqa: BLE001 - exception identity is oracle output
        result = {
            "status": "error",
            "error_type": type(error).__name__,
            "error_message": str(error),
        }
    return {
        "atoms_to_label": atoms,
        "bonds_to_label": bonds,
        "max_recursive_iterations": limit,
        "result": result,
        "state": snapshot(mol),
    }


def record_for_case(case: dict[str, Any]) -> dict[str, Any]:
    source = {
        "fixture": case["fixture"],
        "input_kind": case.get("input_kind", "smiles"),
        "input": case["input"],
        "sanitize": bool(case.get("sanitize", True)),
        "remove_hs": bool(case.get("remove_hs", True)),
        "strict_parsing": bool(case.get("strict_parsing", True)),
        "update_property_cache": bool(case.get("update_property_cache", False)),
        "strict_property_cache": bool(case.get("strict_property_cache", False)),
    }
    try:
        mol = parse_case(case)
        if mol is None:
            return {
                "schema_version": 1,
                "case_id": case["case_id"],
                "source": source,
                "parse_status": "none",
                "parse_error": None,
                "initial_state": None,
                "calls": [],
            }
        if case.get("update_property_cache", False):
            mol.UpdatePropertyCache(strict=bool(case.get("strict_property_cache", False)))
        apply_initial_state(mol, case.get("initial_state", {}))
        initial = snapshot(mol)
        calls = [
            execute_call(mol, call)
            for call in case.get("calls", default_calls(mol))
        ]
        return {
            "schema_version": 1,
            "case_id": case["case_id"],
            "source": source,
            "parse_status": "ok",
            "parse_error": None,
            "initial_state": initial,
            "calls": calls,
        }
    except Exception as error:  # noqa: BLE001 - parse exception is oracle output
        return {
            "schema_version": 1,
            "case_id": case["case_id"],
            "source": source,
            "parse_status": "error",
            "parse_error": f"{type(error).__name__}: {error}",
            "initial_state": None,
            "calls": [],
        }


def corpus_cases(path: Path) -> Iterable[dict[str, Any]]:
    fixture = fixture_reference(path)
    for row, raw in enumerate(path.read_text(encoding="utf-8").splitlines()):
        text = raw.strip()
        if text and not text.startswith("#"):
            yield {
                "case_id": f"corpus-{row}",
                "fixture": f"{fixture}#{row}",
                "input_kind": "smiles",
                "input": text,
                "sanitize": True,
            }


def load_cases(path: Path) -> list[dict[str, Any]]:
    if path.suffix.lower() != ".json":
        return list(corpus_cases(path))
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema_version") != 1 or not isinstance(
        document.get("cases"), list
    ):
        raise SystemExit("focused CIPLabeler fixture must contain schema_version 1 cases")
    fixture = fixture_reference(path)
    cases = []
    for case in document["cases"]:
        copied = dict(case)
        copied.setdefault("fixture", fixture)
        input_file = copied.pop("input_file", None)
        if input_file is not None:
            source_path = REPO_ROOT / input_file
            copied["input"] = source_path.read_text(encoding="utf-8")
            copied["fixture"] = input_file
        cases.append(copied)
    return cases


def record_shard(
    indexed_cases: list[tuple[int, dict[str, Any]]],
) -> list[tuple[int, dict[str, Any]]]:
    return [(index, record_for_case(case)) for index, case in indexed_cases]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--shards",
        type=int,
        default=16,
        help="deterministic row-index-modulo shard count (default: 16)",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="parallel shard workers; output remains in input order (default: 1)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.shards < 1:
        raise SystemExit("--shards must be positive")
    if args.jobs < 1:
        raise SystemExit("--jobs must be positive")
    assert_rdkit_version()
    RDLogger.DisableLog("rdApp.*")
    cases = load_cases(args.input)
    shards: list[list[tuple[int, dict[str, Any]]]] = [
        [] for _ in range(args.shards)
    ]
    for index, case in enumerate(cases):
        shards[index % args.shards].append((index, case))
    nonempty_shards = [shard for shard in shards if shard]
    if args.jobs == 1 or len(nonempty_shards) <= 1:
        indexed_records = [
            item for shard in nonempty_shards for item in record_shard(shard)
        ]
    else:
        with multiprocessing.Pool(
            processes=min(args.jobs, len(nonempty_shards))
        ) as pool:
            indexed_records = [
                item
                for completed_shard in pool.imap_unordered(
                    record_shard, nonempty_shards, chunksize=1
                )
                for item in completed_shard
            ]
    indexed_records.sort(key=lambda item: item[0])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for _, record in indexed_records:
            handle.write(json.dumps(record, sort_keys=True, separators=(",", ":")))
            handle.write("\n")


if __name__ == "__main__":
    main()
