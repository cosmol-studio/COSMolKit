#!/usr/bin/env python3
"""Generate lightweight RDKit UFF/MMFF parameter-coverage golden data."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from _generate_forcefield_params_golden import (
    assert_rdkit_version,
    forcefield_result,
    iter_smiles,
    mmff_result,
)


def failed_mmff_result(error: str) -> dict[str, object]:
    return {
        "ok": False,
        "has_all": None,
        "atom_types": None,
        "formal_charges": None,
        "partial_charges": None,
        "error": error,
    }


def build_record(smiles: str) -> dict[str, object]:
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        error = "MolFromSmiles returned None"
        return {
            "smiles": smiles,
            "rdkit_ok": False,
            "uff": {"ok": False, "has_all": None, "error": error},
            "mmff": failed_mmff_result(error),
            "uff_explicit_h": {"ok": False, "has_all": None, "error": error},
            "mmff_explicit_h": failed_mmff_result(error),
            "error": error,
        }

    explicit_h_mol = Chem.AddHs(mol)
    return {
        "smiles": smiles,
        "rdkit_ok": True,
        "uff": forcefield_result(AllChem.UFFHasAllMoleculeParams, mol),
        "mmff": mmff_result(mol),
        "uff_explicit_h": forcefield_result(
            AllChem.UFFHasAllMoleculeParams, explicit_h_mol
        ),
        "mmff_explicit_h": mmff_result(explicit_h_mol),
        "error": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    assert_rdkit_version()
    RDLogger.DisableLog("rdApp.*")
    args.output.parent.mkdir(parents=True, exist_ok=True)

    count = 0
    with args.output.open("w", encoding="utf-8") as handle:
        for smiles in iter_smiles(args.input):
            handle.write(json.dumps(build_record(smiles), ensure_ascii=True))
            handle.write("\n")
            count += 1

    print(f"Wrote {count} records to {args.output}")


if __name__ == "__main__":
    main()
