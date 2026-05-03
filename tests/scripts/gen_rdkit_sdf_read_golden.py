#!/usr/bin/env python3
"""Generate RDKit SDF read-parity golden data.

Each input SMILES is expanded across:
- 2D and 3D coordinate sources
- V2000 and V3000 molfile encodings
- original stereo markers and marker-stripped coordinate-only records

The marker-stripped 3D cases exercise RDKit's coordinate-inferred chirality
path. The marker-stripped 2D cases document that coordinates alone do not
carry tetrahedral chirality.
"""

from __future__ import annotations

import argparse
import json
import re
from importlib.metadata import version
from pathlib import Path
from typing import Iterable

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

EXPECTED_RDKIT_VERSION = "2026.3.1"


def assert_rdkit_version() -> None:
    actual = version("rdkit")
    assert actual == EXPECTED_RDKIT_VERSION, (
        f"RDKit version mismatch: expected {EXPECTED_RDKIT_VERSION}, got {actual}"
    )


def iter_smiles(path: Path) -> Iterable[str]:
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.replace("\u200b", "").strip()
        if not line or line.startswith("#"):
            continue
        yield line


def strip_v2000_stereo_markers(block: str) -> str:
    lines = block.splitlines()
    if len(lines) < 4:
        return block
    counts = lines[3]
    n_atoms = int(counts[0:3])
    n_bonds = int(counts[3:6])
    atom_start = 4
    bond_start = atom_start + n_atoms
    stripped = list(lines)
    for idx in range(atom_start, bond_start):
        line = stripped[idx]
        if len(line) >= 42:
            stripped[idx] = f"{line[:39]}  0{line[42:]}"
    for idx in range(bond_start, bond_start + n_bonds):
        line = stripped[idx]
        if len(line) >= 12:
            stripped[idx] = f"{line[:9]}  0{line[12:]}"
    return "\n".join(stripped) + "\n"


def strip_v3000_stereo_markers(block: str) -> str:
    stripped: list[str] = []
    for line in block.splitlines():
        if line.startswith("M  V30 "):
            line = re.sub(r" CFG=-?\d+", "", line)
        stripped.append(line)
    return "\n".join(stripped) + "\n"


def strip_stereo_markers(block: str, fmt: str) -> str:
    if fmt == "V2000":
        return strip_v2000_stereo_markers(block)
    return strip_v3000_stereo_markers(block)


def chiral_tags(mol: Chem.Mol) -> list[str]:
    return [str(atom.GetChiralTag()) for atom in mol.GetAtoms()]


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


def make_2d_source(mol: Chem.Mol) -> tuple[Chem.Mol | None, str | None]:
    source = Chem.Mol(mol)
    try:
        AllChem.Compute2DCoords(source)
        if source.GetNumConformers():
            Chem.WedgeMolBonds(source, source.GetConformer())
        return source, None
    except Exception as exc:
        return None, str(exc)


def make_3d_source(mol: Chem.Mol) -> tuple[Chem.Mol | None, str | None]:
    source = Chem.Mol(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    try:
        status = AllChem.EmbedMolecule(source, params)
    except Exception as exc:
        return None, str(exc)
    if status != 0:
        return None, f"EmbedMolecule failed with status {status}"
    return source, None


def render_molblock(mol: Chem.Mol, fmt: str) -> tuple[str | None, str | None]:
    try:
        return Chem.MolToMolBlock(mol, forceV3000=(fmt == "V3000")), None
    except Exception as exc:
        return None, str(exc)


def parse_molblock(block: str) -> tuple[Chem.Mol | None, str | None]:
    try:
        parsed = Chem.MolFromMolBlock(block, sanitize=True, removeHs=False)
    except Exception as exc:
        return None, str(exc)
    if parsed is None:
        return None, "MolFromMolBlock returned None"
    return parsed, None


def build_case(
    smiles: str,
    source: Chem.Mol | None,
    source_error: str | None,
    *,
    dimension: str,
    fmt: str,
    markers: str,
) -> dict[str, object]:
    case_id = f"{dimension.lower()}_{fmt.lower()}_{markers}"
    base = {
        "smiles": smiles,
        "case_id": case_id,
        "dimension": dimension,
        "format": fmt,
        "stereo_markers": markers,
    }
    if source is None:
        return {
            **base,
            "rdkit_ok": False,
            "sdf": None,
            "error": source_error or "source molecule generation failed",
        }

    block, render_error = render_molblock(source, fmt)
    if block is None:
        return {
            **base,
            "rdkit_ok": False,
            "sdf": None,
            "error": render_error,
        }
    if markers == "coords_only":
        block = strip_stereo_markers(block, fmt)

    parsed, parse_error = parse_molblock(block)
    if parsed is None:
        return {
            **base,
            "rdkit_ok": False,
            "sdf": block + "$$$$\n",
            "error": parse_error,
        }

    return {
        **base,
        "rdkit_ok": True,
        "sdf": block + "$$$$\n",
        "atoms": atom_features(parsed),
        "bonds": bond_features(parsed),
        "positions": positions(parsed),
        "chiral_tags": chiral_tags(parsed),
        "smiles_out": smiles_pair(parsed),
        "error": None,
    }


def build_records(smiles: str) -> list[dict[str, object]]:
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        return [
            {
                "smiles": smiles,
                "case_id": case_id,
                "dimension": dimension,
                "format": fmt,
                "stereo_markers": markers,
                "rdkit_ok": False,
                "sdf": None,
                "error": "MolFromSmiles returned None",
            }
            for dimension in ["2D", "3D"]
            for fmt in ["V2000", "V3000"]
            for markers in ["with_markers", "coords_only"]
            for case_id in [f"{dimension.lower()}_{fmt.lower()}_{markers}"]
        ]

    source_2d, error_2d = make_2d_source(mol)
    source_3d, error_3d = make_3d_source(mol)
    records: list[dict[str, object]] = []
    for dimension, source, source_error in [
        ("2D", source_2d, error_2d),
        ("3D", source_3d, error_3d),
    ]:
        for fmt in ["V2000", "V3000"]:
            for markers in ["with_markers", "coords_only"]:
                records.append(
                    build_case(
                        smiles,
                        source,
                        source_error,
                        dimension=dimension,
                        fmt=fmt,
                        markers=markers,
                    )
                )
    return records


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
        default=Path("tests/golden/sdf_read.jsonl"),
        help="output JSONL path (default: tests/golden/sdf_read.jsonl)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="optional maximum input SMILES rows for local debugging",
    )
    args = parser.parse_args()

    assert_rdkit_version()
    RDLogger.DisableLog("rdApp.*")
    rows = []
    for idx, smiles in enumerate(iter_smiles(args.input)):
        if args.limit is not None and idx >= args.limit:
            break
        rows.extend(build_records(smiles))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for record in rows:
            handle.write(json.dumps(record, ensure_ascii=True))
            handle.write("\n")

    print(f"Wrote {len(rows)} records to {args.output}")


if __name__ == "__main__":
    main()
