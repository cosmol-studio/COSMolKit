#!/usr/bin/env python3
"""Generate pinned-RDKit MolToInchi golden data for a SMILES corpus."""

from __future__ import annotations

import argparse
import hashlib
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


def load_smiles(path: Path) -> list[str]:
    return [
        line
        for raw in path.read_text(encoding="utf-8").splitlines()
        if (line := raw.strip()) and not line.startswith("#")
    ]


def explicit_valence(atom: Chem.Atom) -> int:
    return int(atom.GetValence(Chem.ValenceType.EXPLICIT))


def implicit_valence(atom: Chem.Atom) -> int:
    return int(atom.GetValence(Chem.ValenceType.IMPLICIT))


def mol_from_inchi_digest(molecule: Chem.Mol) -> str:
    atoms = []
    for atom in molecule.GetAtoms():
        explicit = explicit_valence(atom)
        implicit = implicit_valence(atom)
        atoms.append(
            [
                atom.GetAtomicNum(),
                atom.GetFormalCharge(),
                str(atom.GetChiralTag()),
                atom.GetIsotope() or None,
                atom.GetIsAromatic(),
                atom.GetNumExplicitHs(),
                atom.GetDegree(),
                explicit,
                implicit,
                atom.GetTotalNumHs(),
                explicit + implicit,
                atom.GetNumRadicalElectrons(),
                atom.GetNoImplicit(),
            ]
        )
    bonds = [
        [
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            str(bond.GetBondType()),
            str(bond.GetBondDir()),
            str(bond.GetStereo()).removeprefix("STEREO"),
            list(bond.GetStereoAtoms()),
            bond.GetIsAromatic(),
        ]
        for bond in molecule.GetBonds()
    ]
    payload = [Chem.MolToSmiles(Chem.Mol(molecule)), atoms, bonds]
    encoded = json.dumps(payload, ensure_ascii=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def mol_from_inchi_branches(inchi: str) -> dict[str, dict[str, object]]:
    branches = {}
    for sanitize in (False, True):
        for remove_hs in (False, True):
            name = f"sanitize{int(sanitize)}_remove_hs{int(remove_hs)}"
            try:
                molecule = Chem.MolFromInchi(
                    inchi, sanitize=sanitize, removeHs=remove_hs
                )
            except Exception as error:
                branches[name] = {
                    "status": "error",
                    "digest": None,
                    "error_kind": type(error).__name__,
                }
                continue
            if molecule is None:
                branches[name] = {
                    "status": "none",
                    "digest": None,
                    "error_kind": None,
                }
                continue
            branches[name] = {
                "status": "ok",
                "digest": mol_from_inchi_digest(molecule),
                "error_kind": None,
            }
    return branches


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("testdata/smiles/corpus/smiles_small.smi"),
        help="input SMILES file (default: testdata/smiles/corpus/smiles_small.smi)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("testdata/inchi/expected/rdkit/smiles_small/inchi.jsonl"),
        help="output JSONL path",
    )
    args = parser.parse_args()

    assert_rdkit_version()
    RDLogger.DisableLog("rdApp.*")
    smiles_rows = load_smiles(args.input)
    records: list[dict[str, object]] = []

    for row, smiles in enumerate(smiles_rows, start=1):
        molecule = Chem.MolFromSmiles(smiles, sanitize=True)
        if molecule is None:
            records.append(
                {
                    "row": row,
                    "smiles": smiles,
                    "rdkit_ok": False,
                    "inchi": None,
                    "error": "MolFromSmiles returned None",
                    "mol_from_inchi_branches": {},
                }
            )
            continue

        try:
            inchi = Chem.MolToInchi(molecule)
        except Exception as error:
            records.append(
                {
                    "row": row,
                    "smiles": smiles,
                    "rdkit_ok": False,
                    "inchi": None,
                    "error": f"{type(error).__name__}: {error}",
                    "mol_from_inchi_branches": {},
                }
            )
            continue

        if not inchi:
            records.append(
                {
                    "row": row,
                    "smiles": smiles,
                    "rdkit_ok": False,
                    "inchi": None,
                    "error": "MolToInchi returned an empty string",
                    "mol_from_inchi_branches": {},
                }
            )
            continue

        records.append(
            {
                "row": row,
                "smiles": smiles,
                "rdkit_ok": True,
                "inchi": inchi,
                "error": None,
                "mol_from_inchi_branches": mol_from_inchi_branches(inchi),
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as output:
        for record in records:
            output.write(json.dumps(record, ensure_ascii=True))
            output.write("\n")

    print(f"Wrote {len(records)} records to {args.output}")


if __name__ == "__main__":
    main()
