#!/usr/bin/env python3
"""Generate exact RDKit Layered fingerprint expected data for a SMILES corpus."""

from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path
from typing import Any, Iterable

from rdkit import Chem, DataStructs, RDLogger


EXPECTED_RDKIT_VERSION = "2026.3.1"
EXPECTED_ALGORITHM_VERSION = "0.7.0"
PROFILE_PATH = Path(__file__).with_name("layered_fingerprint_profile.json")


def assert_rdkit_version() -> None:
    actual = version("rdkit")
    if actual != EXPECTED_RDKIT_VERSION:
        raise SystemExit(
            f"RDKit version mismatch: expected {EXPECTED_RDKIT_VERSION}, got {actual}"
        )


def iter_smiles(path: Path) -> Iterable[str]:
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line and not line.startswith("#"):
            yield line


def load_profile() -> dict[str, Any]:
    profile = json.loads(PROFILE_PATH.read_text(encoding="utf-8"))
    if profile.get("schema_version") != 1:
        raise SystemExit("Layered fingerprint profile schema_version must be 1")
    if profile.get("rdkit_version") != EXPECTED_RDKIT_VERSION:
        raise SystemExit("Layered fingerprint profile has the wrong RDKit version")
    if profile.get("algorithm_version") != EXPECTED_ALGORITHM_VERSION:
        raise SystemExit("Layered fingerprint profile has the wrong algorithm version")
    source_paths = profile.get("source_paths")
    if not isinstance(source_paths, list) or not source_paths:
        raise SystemExit("Layered fingerprint profile must identify its RDKit sources")
    comparison_fields = profile.get("comparison_fields")
    if comparison_fields != ["num_bits", "on_bits", "atom_counts"]:
        raise SystemExit("Layered fingerprint comparison fields are incomplete")
    branches = profile.get("branches")
    if not isinstance(branches, list) or not branches:
        raise SystemExit("Layered fingerprint profile must contain branches")
    names = [branch.get("name") for branch in branches]
    if any(not isinstance(name, str) or not name for name in names):
        raise SystemExit("Layered fingerprint branch names must be non-empty strings")
    if len(set(names)) != len(names):
        raise SystemExit("Layered fingerprint branch names must be unique")
    return profile


def resolved_from_atoms(mol: Chem.Mol, mode: str | None) -> list[int] | None:
    if mode is None:
        return None
    atom_count = mol.GetNumAtoms()
    if mode == "first":
        return [0] if atom_count else []
    if mode == "terminal_pair":
        if atom_count == 0:
            return []
        if atom_count == 1:
            return [0]
        return [0, atom_count - 1]
    raise ValueError(f"unknown fromAtoms profile {mode!r}")


def resolved_atom_counts(atom_count: int, mode: str | None) -> list[int] | None:
    if mode is None:
        return None
    if mode == "zeros":
        return [0] * atom_count
    if mode == "index_plus_10":
        return [index + 10 for index in range(atom_count)]
    raise ValueError(f"unknown atomCounts profile {mode!r}")


def resolved_mask(fp_size: int, mode: str | None) -> DataStructs.ExplicitBitVect | None:
    if mode is None:
        return None
    mask = DataStructs.ExplicitBitVect(fp_size)
    if mode == "even":
        indexes = range(0, fp_size, 2)
    elif mode == "mod_three":
        indexes = range(0, fp_size, 3)
    else:
        raise ValueError(f"unknown setOnlyBits profile {mode!r}")
    for index in indexes:
        mask.SetBit(index)
    return mask


def call_fingerprint(mol: Chem.Mol, branch: dict[str, Any]) -> dict[str, Any]:
    fp_size = int(branch["fpSize"])
    layer_flags = int(branch["layerFlags"], 0)
    from_atoms = resolved_from_atoms(mol, branch.get("fromAtoms"))
    atom_counts = resolved_atom_counts(mol.GetNumAtoms(), branch.get("atomCounts"))
    atom_count_seed = None if atom_counts is None else list(atom_counts)
    set_only_bits = resolved_mask(fp_size, branch.get("setOnlyBits"))
    kwargs: dict[str, Any] = {
        "layerFlags": layer_flags,
        "minPath": int(branch["minPath"]),
        "maxPath": int(branch["maxPath"]),
        "fpSize": fp_size,
        "branchedPaths": bool(branch["branchedPaths"]),
    }
    if from_atoms is not None:
        kwargs["fromAtoms"] = from_atoms
    if atom_counts is not None:
        kwargs["atomCounts"] = atom_counts
    if set_only_bits is not None:
        kwargs["setOnlyBits"] = set_only_bits

    fingerprint = Chem.LayeredFingerprint(mol, **kwargs)
    return {
        "parameters": branch,
        "resolved_arguments": {
            "layerFlags": layer_flags,
            "fromAtoms": from_atoms,
            "atomCounts": atom_count_seed,
            "setOnlyBits": None
            if set_only_bits is None
            else list(set_only_bits.GetOnBits()),
        },
        "ok": True,
        "num_bits": int(fingerprint.GetNumBits()),
        "on_bits": list(fingerprint.GetOnBits()),
        "atom_counts": atom_counts,
        "error": None,
    }


def build_records(
    smiles_values: list[str], profile: dict[str, Any]
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for row, smiles in enumerate(smiles_values):
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        if mol is None:
            records.append(
                {
                    "row": row,
                    "smiles": smiles,
                    "rdkit_ok": False,
                    "branches": {},
                    "error": "MolFromSmiles returned None",
                }
            )
            continue

        branch_records: dict[str, Any] = {}
        for branch in profile["branches"]:
            name = branch["name"]
            try:
                branch_records[name] = call_fingerprint(mol, branch)
            except Exception as exc:  # noqa: BLE001
                branch_records[name] = {
                    "parameters": branch,
                    "resolved_arguments": None,
                    "ok": False,
                    "num_bits": None,
                    "on_bits": None,
                    "atom_counts": None,
                    "error": f"{type(exc).__name__}: {exc}",
                }
        records.append(
            {
                "row": row,
                "smiles": smiles,
                "rdkit_ok": True,
                "branches": branch_records,
                "error": None,
            }
        )
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    assert_rdkit_version()
    RDLogger.DisableLog("rdApp.*")
    profile = load_profile()
    records = build_records(list(iter_smiles(args.input)), profile)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, ensure_ascii=True, sort_keys=True))
            handle.write("\n")
    print(f"Wrote {len(records)} records to {args.output}")


if __name__ == "__main__":
    main()
