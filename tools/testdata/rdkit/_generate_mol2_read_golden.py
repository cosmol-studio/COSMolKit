#!/usr/bin/env python3
"""Generate RDKit MOL2 read-parity golden data from RDKit MOL2 fixtures."""

from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path

from rdkit import Chem, RDLogger

EXPECTED_RDKIT_VERSION = "2026.3.1"


def assert_rdkit_version() -> None:
    actual = version("rdkit")
    assert actual == EXPECTED_RDKIT_VERSION, (
        f"RDKit version mismatch: expected {EXPECTED_RDKIT_VERSION}, got {actual}"
    )


def atom_features(mol: Chem.Mol) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for atom in mol.GetAtoms():
        isotope = atom.GetIsotope()
        atom_map = atom.GetAtomMapNum()
        rows.append(
            {
                "atomic_num": atom.GetAtomicNum(),
                "isotope": isotope if isotope else None,
                "formal_charge": atom.GetFormalCharge(),
                "is_aromatic": atom.GetIsAromatic(),
                "atom_map_num": atom_map if atom_map else None,
                "chiral_tag": str(atom.GetChiralTag()),
                "tripos_atom_name": atom.GetProp("_TriposAtomName")
                if atom.HasProp("_TriposAtomName")
                else None,
                "tripos_atom_type": atom.GetProp("_TriposAtomType")
                if atom.HasProp("_TriposAtomType")
                else None,
                "tripos_partial_charge": atom.GetProp("_TriposPartialCharge")
                if atom.HasProp("_TriposPartialCharge")
                else None,
            }
        )
    return rows


def bond_features(mol: Chem.Mol) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for bond in mol.GetBonds():
        rows.append(
            {
                "begin": bond.GetBeginAtomIdx(),
                "end": bond.GetEndAtomIdx(),
                "bond_type": str(bond.GetBondType()),
                "is_aromatic": bond.GetIsAromatic(),
                "direction": str(bond.GetBondDir()),
                "stereo": str(bond.GetStereo()),
                "stereo_atoms": list(bond.GetStereoAtoms()),
                "tripos_bond_type": bond.GetProp("_TriposBondType")
                if bond.HasProp("_TriposBondType")
                else None,
            }
        )
    return rows


def positions(mol: Chem.Mol) -> list[list[float]] | None:
    if mol.GetNumConformers() == 0:
        return None
    return mol.GetConformer().GetPositions().tolist()


def smiles_pair(mol: Chem.Mol) -> dict[str, str]:
    return {
        "canonical": Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True),
        "noncanonical": Chem.MolToSmiles(mol, isomericSmiles=True, canonical=False),
    }


def parse_mol2(
    block: str,
    *,
    sanitize: bool,
    remove_hs: bool,
    cleanup_substructures: bool,
) -> tuple[Chem.Mol | None, str | None]:
    try:
        mol = Chem.MolFromMol2Block(
            block,
            sanitize=sanitize,
            removeHs=remove_hs,
            cleanupSubstructures=cleanup_substructures,
        )
    except Exception as exc:
        return None, str(exc)
    if mol is None:
        return None, "MolFromMol2Block returned None"
    return mol, None


def build_case(
    fixture_path: Path,
    root: Path,
    *,
    sanitize: bool,
    remove_hs: bool,
    cleanup_substructures: bool,
) -> dict[str, object]:
    block = fixture_path.read_text(encoding="utf-8")
    case_id = (
        f"sanitize_{str(sanitize).lower()}_"
        f"remove_hs_{str(remove_hs).lower()}_"
        f"cleanup_{str(cleanup_substructures).lower()}"
    )
    rel_path = fixture_path.relative_to(root).as_posix()
    base = {
        "fixture": rel_path,
        "case_id": case_id,
        "sanitize": sanitize,
        "remove_hs": remove_hs,
        "cleanup_substructures": cleanup_substructures,
        "variant": "CORINA",
        "mol2": block,
    }
    mol, error = parse_mol2(
        block,
        sanitize=sanitize,
        remove_hs=remove_hs,
        cleanup_substructures=cleanup_substructures,
    )
    if mol is None:
        return {
            **base,
            "rdkit_ok": False,
            "rdkit_null": False,
            "atoms": None,
            "bonds": None,
            "positions": None,
            "chiral_tags": None,
            "smiles_out": None,
            "error": error,
        }
    return {
        **base,
        "rdkit_ok": True,
        "rdkit_null": False,
        "atoms": atom_features(mol),
        "bonds": bond_features(mol),
        "positions": positions(mol),
        "chiral_tags": [str(atom.GetChiralTag()) for atom in mol.GetAtoms()],
        "smiles_out": smiles_pair(mol),
        "error": None,
    }


def iter_fixtures(input_dir: Path) -> list[Path]:
    return sorted(path for path in input_dir.glob("*.mol2") if path.is_file())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("testdata/mol2/fixtures/rdkit"),
        help=(
            "input MOL2 fixture directory (default: testdata/mol2/fixtures/rdkit). "
            "When the unified golden runner passes testdata/smiles/corpus/smiles_small.smi, this script "
            "falls back to the default fixture directory."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("testdata/mol2/expected/rdkit/smiles_small/mol2_read.jsonl"),
        help="output JSONL path (default: testdata/mol2/expected/rdkit/smiles_small/mol2_read.jsonl)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="optional maximum fixture count for local debugging",
    )
    args = parser.parse_args()

    assert_rdkit_version()
    RDLogger.DisableLog("rdApp.*")

    input_dir = args.input
    if input_dir.is_file():
        input_dir = Path("testdata/mol2/fixtures/rdkit")

    fixtures = iter_fixtures(input_dir)
    if args.limit is not None:
        fixtures = fixtures[: args.limit]
    rows: list[dict[str, object]] = []
    for fixture in fixtures:
        for sanitize, remove_hs, cleanup_substructures in [
            (True, True, True),
            (True, False, True),
            (False, False, True),
            (False, False, False),
        ]:
            rows.append(
                build_case(
                    fixture,
                    input_dir,
                    sanitize=sanitize,
                    remove_hs=remove_hs,
                    cleanup_substructures=cleanup_substructures,
                )
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for record in rows:
            handle.write(json.dumps(record, ensure_ascii=True))
            handle.write("\n")

    print(f"Wrote {len(rows)} records to {args.output}")


if __name__ == "__main__":
    main()
