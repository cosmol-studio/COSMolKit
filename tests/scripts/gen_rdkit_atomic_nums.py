#!/usr/bin/env python3
"""Generate RDKit graph-feature golden data (direct + AddHs) from tests/smiles.txt."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

from rdkit import Chem


def atom_features(mol: Chem.Mol) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for atom in mol.GetAtoms():
        rows.append(
            {
                "atomic_num": atom.GetAtomicNum(),
                "chirality": str(atom.GetChiralTag()),
                "degree": atom.GetTotalDegree(),
                "formal_charge": atom.GetFormalCharge(),
                "num_hs": atom.GetTotalNumHs(),
                "num_radical_electrons": atom.GetNumRadicalElectrons(),
                "hybridization": str(atom.GetHybridization()),
                "is_aromatic": atom.GetIsAromatic(),
                "is_in_ring": atom.IsInRing(),
                "explicit_valence": atom.GetExplicitValence(),
                "implicit_hs": atom.GetNumImplicitHs(),
                "total_valence": atom.GetTotalValence(),
            }
        )
    return rows


def bond_features(mol: Chem.Mol) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for bond in mol.GetBonds():
        rows.append(
            {
                "begin_atom": bond.GetBeginAtomIdx(),
                "end_atom": bond.GetEndAtomIdx(),
                "bond_type": str(bond.GetBondType()),
                "stereo": str(bond.GetStereo()),
                "is_conjugated": bond.GetIsConjugated(),
            }
        )
    return rows


def iter_smiles(path: Path) -> Iterable[str]:
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        yield line


def build_record(smiles: str) -> dict[str, object]:
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        return {
            "smiles": smiles,
            "rdkit_ok": False,
            "direct": None,
            "with_hs": None,
            "error": "MolFromSmiles returned None",
        }

    mol_h = Chem.AddHs(Chem.Mol(mol))
    return {
        "smiles": smiles,
        "rdkit_ok": True,
        "direct": {
            "atom_features": atom_features(mol),
            "bond_features": bond_features(mol),
        },
        "with_hs": {
            "atom_features": atom_features(mol_h),
            "bond_features": bond_features(mol_h),
        },
        "error": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("tests/smiles.txt"),
        help="input SMILES file (default: tests/smiles.txt)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("tests/golden/atomic_nums.jsonl"),
        help="output JSONL path (default: tests/golden/atomic_nums.jsonl)",
    )
    args = parser.parse_args()

    smiles_rows = list(iter_smiles(args.input))
    args.output.parent.mkdir(parents=True, exist_ok=True)

    with args.output.open("w", encoding="utf-8") as f:
        for smiles in smiles_rows:
            record = build_record(smiles)
            f.write(json.dumps(record, ensure_ascii=True))
            f.write("\n")

    print(f"Wrote {len(smiles_rows)} records to {args.output}")


if __name__ == "__main__":
    main()
