#!/usr/bin/env python3
"""Generate RDKit DeleteSubstructs onlyFrags/chirality golden data."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


CASES: tuple[dict[str, object], ...] = (
    {
        "case": "carbonyl_delete_all_matches",
        "smiles": "CCC(=O).C=O",
        "smarts": "C=O",
        "only_frags": False,
        "use_chirality": False,
    },
    {
        "case": "chloride_fragment_delete_only_whole_fragment",
        "smiles": "CCO.Cl",
        "smarts": "[Cl;H1&X1,-]",
        "only_frags": True,
        "use_chirality": False,
    },
    {
        "case": "sodium_fragment_no_match_is_copy",
        "smiles": "CCO",
        "smarts": "[Na+]",
        "only_frags": True,
        "use_chirality": False,
    },
    {
        "case": "sodium_fragment_delete_from_acetate_salt",
        "smiles": "CC(=O)[O-].[Na+]",
        "smarts": "[Na+]",
        "only_frags": True,
        "use_chirality": False,
    },
    {
        "case": "acid_query_only_frags_preserves_larger_acid",
        "smiles": "CCC(=O)O.O=CO",
        "smarts": "C(=O)O",
        "only_frags": True,
        "use_chirality": False,
    },
    {
        "case": "acid_query_non_only_frags_deletes_all_matches",
        "smiles": "CCC(=O)O.O=CO",
        "smarts": "C(=O)O",
        "only_frags": False,
        "use_chirality": False,
    },
    {
        "case": "acid_query_no_match_is_copy",
        "smiles": "CCCO",
        "smarts": "C(=O)O",
        "only_frags": False,
        "use_chirality": False,
    },
    {
        "case": "partial_substructure_delete",
        "smiles": "CCOC",
        "smarts": "OC",
        "only_frags": False,
        "use_chirality": False,
    },
    {
        "case": "partial_substructure_only_frags_is_copy",
        "smiles": "CCOC",
        "smarts": "OC",
        "only_frags": True,
        "use_chirality": False,
    },
    {
        "case": "chlorine_only_frags_deletes_separate_chlorine",
        "smiles": "CCOCCl.Cl",
        "smarts": "Cl",
        "only_frags": True,
        "use_chirality": False,
    },
    {
        "case": "chlorine_non_only_frags_deletes_all_chlorine_atoms",
        "smiles": "CCOCCl.Cl",
        "smarts": "Cl",
        "only_frags": False,
        "use_chirality": False,
    },
    # Source-backed chiral cases from
    # third_party/rdkit/Code/GraphMol/ChemTransforms/testChemTransforms.cpp
    # lines 102-128.
    {
        "case": "rdkit_chiral_query_matching_center_use_chirality_true",
        "smiles": "CCO[C@H](N)(P)",
        "smarts": "O[C@H](N)(P)",
        "only_frags": False,
        "use_chirality": True,
    },
    {
        "case": "rdkit_chiral_query_matching_center_use_chirality_false",
        "smiles": "CCO[C@H](N)(P)",
        "smarts": "O[C@H](N)(P)",
        "only_frags": False,
        "use_chirality": False,
    },
    {
        "case": "rdkit_chiral_query_opposite_center_use_chirality_true",
        "smiles": "CCO[C@H](N)(P)",
        "smarts": "O[C@@H](N)(P)",
        "only_frags": False,
        "use_chirality": True,
    },
    {
        "case": "rdkit_chiral_query_opposite_center_use_chirality_false",
        "smiles": "CCO[C@H](N)(P)",
        "smarts": "O[C@@H](N)(P)",
        "only_frags": False,
        "use_chirality": False,
    },
)


def build_record(case: dict[str, object]) -> dict[str, object]:
    smiles = str(case["smiles"])
    smarts = str(case["smarts"])
    only_frags = bool(case["only_frags"])
    use_chirality = bool(case["use_chirality"])
    mol = Chem.MolFromSmiles(smiles)
    query = Chem.MolFromSmarts(smarts)
    if mol is None or query is None:
        return {
            **case,
            "rdkit_ok": False,
            "result_smiles": None,
            "num_atoms": None,
            "num_bonds": None,
            "error": "MolFromSmiles or MolFromSmarts returned None",
        }
    try:
        result = Chem.DeleteSubstructs(
            mol, query, onlyFrags=only_frags, useChirality=use_chirality
        )
        return {
            **case,
            "rdkit_ok": True,
            "result_smiles": Chem.MolToSmiles(result, isomericSmiles=True, canonical=True),
            "num_atoms": result.GetNumAtoms(),
            "num_bonds": result.GetNumBonds(),
            "error": None,
        }
    except Exception as exc:
        return {
            **case,
            "rdkit_ok": False,
            "result_smiles": None,
            "num_atoms": None,
            "num_bonds": None,
            "error": type(exc).__name__ + ": " + str(exc),
        }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("tests/smiles.smi"),
        help="accepted for unified generator compatibility; not used",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("tests/golden/delete_substructs_onlyfrags_chirality.jsonl"),
        help="output JSONL path",
    )
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as f:
        for case in CASES:
            f.write(json.dumps(build_record(case), ensure_ascii=True))
            f.write("\n")
    print(f"Wrote {len(CASES)} DeleteSubstructs records to {args.output}")


if __name__ == "__main__":
    main()
