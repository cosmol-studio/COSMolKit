#!/usr/bin/env python3
"""Generate exact RDKit RDKFingerprint expected data for a corpus."""

from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path
from typing import Iterable

from rdkit import Chem, RDLogger

EXPECTED_RDKIT_VERSION = "2026.3.1"
PROFILE_PATH = Path(__file__).with_name("rdkit_topological_fingerprint_profile.json")


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


def load_profile() -> list[dict[str, object]]:
    profile = json.loads(PROFILE_PATH.read_text(encoding="utf-8"))
    if profile.get("rdkit_version") != EXPECTED_RDKIT_VERSION:
        raise SystemExit("topological fingerprint profile is not pinned to RDKit 2026.3.1")
    branches = profile.get("branches")
    if not isinstance(branches, list) or not branches:
        raise SystemExit("topological fingerprint profile must contain branches")
    names = [branch.get("name") for branch in branches]
    if any(not isinstance(name, str) for name in names) or len(set(names)) != len(names):
        raise SystemExit("topological fingerprint profile branch names must be unique")
    return branches


def call_fingerprint(mol: Chem.Mol, branch: dict[str, object]):
    kwargs = {
        "minPath": branch["minPath"],
        "maxPath": branch["maxPath"],
        "fpSize": branch["fpSize"],
        "nBitsPerHash": branch["nBitsPerHash"],
        "useHs": branch["useHs"],
        "tgtDensity": branch["tgtDensity"],
        "minSize": branch["minSize"],
        "branchedPaths": branch["branchedPaths"],
        "useBondOrder": branch["useBondOrder"],
    }
    if branch.get("fromAtoms") == "first" and mol.GetNumAtoms() > 0:
        kwargs["fromAtoms"] = [0]
    if branch.get("atomInvariants") == "index_plus_one":
        kwargs["atomInvariants"] = [index + 1 for index in range(mol.GetNumAtoms())]
    atom_bits: list[list[int]] = []
    bit_info: dict[int, list[list[int]]] = {}
    if branch.get("provenance"):
        kwargs["atomBits"] = atom_bits
        kwargs["bitInfo"] = bit_info
    fp = Chem.RDKFingerprint(mol, **kwargs)
    record: dict[str, object] = {
        "ok": True,
        "on_bits": list(fp.GetOnBits()),
        "num_bits": fp.GetNumBits(),
        "num_on_bits": fp.GetNumOnBits(),
        "error": None,
    }
    if branch.get("provenance"):
        record["atom_bits"] = atom_bits
        record["bit_info"] = {
            str(bit): [list(path) for path in paths]
            for bit, paths in sorted(bit_info.items())
        }
    return record


def build_records(smiles_values: list[str], branches: list[dict[str, object]]):
    records = []
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
        branch_records = {}
        for branch in branches:
            name = branch["name"]
            try:
                result = call_fingerprint(mol, branch)
            except Exception as exc:  # noqa: BLE001
                result = {
                    "ok": False,
                    "on_bits": None,
                    "num_bits": None,
                    "num_on_bits": None,
                    "error": f"{type(exc).__name__}: {exc}",
                }
            branch_records[name] = {
                "parameters": branch,
                **result,
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
    branches = load_profile()
    records = build_records(list(iter_smiles(args.input)), branches)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, ensure_ascii=True, sort_keys=True))
            handle.write("\n")
    print(f"Wrote {len(records)} records to {args.output}")


if __name__ == "__main__":
    main()
