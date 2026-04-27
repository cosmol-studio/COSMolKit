#!/usr/bin/env python3
"""Generate RDKit graph-feature golden data (direct + AddHs) from tests/smiles.smi."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

from rdkit import Chem


def explicit_valence(atom: Chem.Atom) -> int:
    if hasattr(Chem, "ValenceType"):
        return int(atom.GetValence(Chem.ValenceType.EXPLICIT))
    return int(atom.GetExplicitValence())


def implicit_valence(atom: Chem.Atom) -> int:
    if hasattr(Chem, "ValenceType"):
        return int(atom.GetValence(Chem.ValenceType.IMPLICIT))
    return int(atom.GetImplicitValence())


def atom_features(mol: Chem.Mol) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for atom in mol.GetAtoms():
        cip_rank = atom.GetProp("_CIPRank") if atom.HasProp("_CIPRank") else None
        explicit = explicit_valence(atom)
        implicit = implicit_valence(atom)
        rows.append(
            {
                "atomic_num": atom.GetAtomicNum(),
                "isotope": atom.GetIsotope() or None,
                "chirality": str(atom.GetChiralTag()),
                "cip_code": atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None,
                "cip_rank": int(cip_rank) if cip_rank is not None else None,
                "chirality_possible": bool(atom.HasProp("_ChiralityPossible")),
                "degree": atom.GetTotalDegree(),
                "formal_charge": atom.GetFormalCharge(),
                "num_hs": atom.GetTotalNumHs(),
                "num_radical_electrons": atom.GetNumRadicalElectrons(),
                "hybridization": str(atom.GetHybridization()),
                "is_aromatic": atom.GetIsAromatic(),
                "is_in_ring": atom.IsInRing(),
                "explicit_valence": explicit,
                "implicit_hs": atom.GetNumImplicitHs(),
                "total_valence": explicit + implicit,
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

    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    mol_h = Chem.AddHs(Chem.Mol(mol))
    Chem.AssignStereochemistry(mol_h, cleanIt=True, force=True)
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
        default=Path("tests/smiles.smi"),
        help="input SMILES file (default: tests/smiles.smi)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("tests/golden/graph_features.jsonl"),
        help="output JSONL path (default: tests/golden/graph_features.jsonl)",
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
