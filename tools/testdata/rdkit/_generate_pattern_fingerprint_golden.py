#!/usr/bin/env python3
"""Generate exact pinned-RDKit Pattern fingerprint expected data."""

from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path
from typing import Any, Iterable

from rdkit import Chem, DataStructs, RDLogger


EXPECTED_RDKIT_VERSION = "2026.3.1"
EXPECTED_FINGERPRINT_VERSION = "1.0.0"
PROFILE_PATH = Path(__file__).with_name("pattern_fingerprint_profile.json")


def assert_rdkit_version() -> None:
    actual = version("rdkit")
    if actual != EXPECTED_RDKIT_VERSION:
        raise SystemExit(
            f"RDKit version mismatch: expected {EXPECTED_RDKIT_VERSION}, got {actual}"
        )


def iter_inputs(path: Path) -> Iterable[tuple[str, str]]:
    for raw in path.read_text(encoding="utf-8").splitlines():
        value = raw.strip()
        if not value or value.startswith("#"):
            continue
        fields = value.split("\t", 1)
        if len(fields) == 2 and fields[0] in {"smiles", "smarts"}:
            yield fields[0], fields[1]
        else:
            yield "smiles", value


def load_profile() -> dict[str, Any]:
    profile = json.loads(PROFILE_PATH.read_text(encoding="utf-8"))
    if profile.get("schema_version") != 1:
        raise SystemExit("Pattern fingerprint profile schema_version must be 1")
    if profile.get("rdkit_version") != EXPECTED_RDKIT_VERSION:
        raise SystemExit("Pattern fingerprint profile has the wrong RDKit version")
    if profile.get("fingerprint_version") != EXPECTED_FINGERPRINT_VERSION:
        raise SystemExit("Pattern fingerprint source version must be 1.0.0")
    source_paths = profile.get("source_paths")
    if not isinstance(source_paths, list) or not source_paths:
        raise SystemExit("Pattern fingerprint profile must identify its RDKit sources")
    fields = profile.get("comparison_fields")
    if not isinstance(fields, list) or not fields:
        raise SystemExit("Pattern fingerprint profile must define comparison fields")
    branches = profile.get("branches")
    if not isinstance(branches, list) or not branches:
        raise SystemExit("Pattern fingerprint profile must contain branches")
    names = [branch.get("name") for branch in branches]
    if any(not isinstance(name, str) or not name for name in names):
        raise SystemExit("Pattern fingerprint branch names must be non-empty strings")
    if len(set(names)) != len(names):
        raise SystemExit("Pattern fingerprint branch names must be unique")
    required = {
        "name",
        "fpSize",
        "tautomericFingerprint",
        "atomCounts",
        "setOnlyBits",
    }
    if any(not required <= set(branch) for branch in branches):
        raise SystemExit("Pattern fingerprint branches are missing required fields")
    return profile


def make_atom_counts(molecule: Chem.Mol, state: str) -> list[int] | None:
    if state == "none":
        return None
    if state == "zeroed":
        return [0] * molecule.GetNumAtoms()
    if state == "index_plus_11":
        return [index + 11 for index in range(molecule.GetNumAtoms())]
    raise ValueError(f"unknown atomCounts state: {state}")


def make_set_only_bits(n_bits: int, state: str):
    if state == "none":
        return None
    width = n_bits + 1 if state == "wrong_width" else n_bits
    mask = DataStructs.ExplicitBitVect(width)
    if state == "all":
        for bit in range(width):
            mask.SetBit(bit)
    elif state == "sparse":
        for bit in (0, width // 2, width - 1):
            mask.SetBit(bit)
    elif state not in {"zero", "wrong_width"}:
        raise ValueError(f"unknown setOnlyBits state: {state}")
    return mask


def fingerprint_branch(
    molecule: Chem.Mol, branch: dict[str, Any]
) -> dict[str, Any]:
    n_bits = int(branch["fpSize"])
    atom_counts = make_atom_counts(molecule, branch["atomCounts"])
    atom_counts_before = None if atom_counts is None else list(atom_counts)
    mask = make_set_only_bits(n_bits, branch["setOnlyBits"])
    mask_bits = None if mask is None else list(mask.GetOnBits())
    try:
        fingerprint = Chem.PatternFingerprint(
            molecule,
            n_bits,
            [] if atom_counts is None else atom_counts,
            mask,
            bool(branch["tautomericFingerprint"]),
        )
        return {
            "ok": True,
            "n_bits": fingerprint.GetNumBits(),
            "on_bits": list(fingerprint.GetOnBits()),
            "atom_counts_before": atom_counts_before,
            "atom_counts_after": None if atom_counts is None else atom_counts,
            "set_only_bits": mask_bits,
            "error": None,
        }
    except Exception as error:  # noqa: BLE001 - the oracle records source errors
        return {
            "ok": False,
            "n_bits": None,
            "on_bits": None,
            "atom_counts_before": atom_counts_before,
            "atom_counts_after": None if atom_counts is None else atom_counts,
            "set_only_bits": mask_bits,
            "error": f"{type(error).__name__}: {error}",
        }


def build_records(
    inputs: list[tuple[str, str]], profile: dict[str, Any]
) -> list[dict[str, Any]]:
    records = []
    for row, (input_kind, smiles) in enumerate(inputs):
        molecule = (
            Chem.MolFromSmarts(smiles)
            if input_kind == "smarts"
            else Chem.MolFromSmiles(smiles, sanitize=True)
        )
        if molecule is None:
            records.append(
                {
                    "row": row,
                    "smiles": smiles,
                    "input_kind": input_kind,
                    "rdkit_ok": False,
                    "branches": {},
                    "error": "MolFromSmiles returned None",
                }
            )
            continue
        branches = {
            branch["name"]: {
                "parameters": branch,
                **fingerprint_branch(molecule, branch),
            }
            for branch in profile["branches"]
        }
        records.append(
            {
                "row": row,
                "smiles": smiles,
                "input_kind": input_kind,
                "rdkit_ok": True,
                "branches": branches,
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
    records = build_records(list(iter_inputs(args.input)), profile)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, ensure_ascii=True, sort_keys=True))
            handle.write("\n")
    print(f"Wrote {len(records)} records to {args.output}")


if __name__ == "__main__":
    main()
