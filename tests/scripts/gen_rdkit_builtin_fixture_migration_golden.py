#!/usr/bin/env python3
"""Generate RDKit golden data for migrated built-in fixture corpora."""

from __future__ import annotations

import argparse
import gzip
import json
from importlib.metadata import version
from pathlib import Path

from rdkit import Chem, RDLogger

EXPECTED_RDKIT_VERSION = "2026.3.1"
REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FIXTURE_ROOT = REPO_ROOT / "tests" / "fixtures" / "rdkit_builtin"


def assert_rdkit_version() -> None:
    actual = version("rdkit")
    assert actual == EXPECTED_RDKIT_VERSION, (
        f"RDKit version mismatch: expected {EXPECTED_RDKIT_VERSION}, got {actual}"
    )


def mol_summary(mol: Chem.Mol | None) -> dict[str, object]:
    if mol is None:
        return {"rdkit_ok": False, "error": "RDKit returned None"}
    return {
        "rdkit_ok": True,
        "atom_count": mol.GetNumAtoms(),
        "bond_count": mol.GetNumBonds(),
        "atomic_numbers": [atom.GetAtomicNum() for atom in mol.GetAtoms()],
        "formal_charges": [atom.GetFormalCharge() for atom in mol.GetAtoms()],
        "bond_types": [str(bond.GetBondType()) for bond in mol.GetBonds()],
        "canonical_smiles": Chem.MolToSmiles(
            mol, isomericSmiles=True, canonical=True
        ),
        "error": None,
    }


def safe_mol_from_mol_file(path: Path) -> dict[str, object]:
    try:
        return mol_summary(Chem.MolFromMolFile(str(path), sanitize=True, removeHs=False))
    except Exception as exc:
        return {"rdkit_ok": False, "error": str(exc)}


def safe_mol_from_mol2_file(path: Path) -> dict[str, object]:
    try:
        return mol_summary(
            Chem.MolFromMol2File(
                str(path), sanitize=True, removeHs=False, cleanupSubstructures=True
            )
        )
    except Exception as exc:
        return {"rdkit_ok": False, "error": str(exc)}


def safe_mol_from_xyz_file(path: Path) -> dict[str, object]:
    try:
        return mol_summary(Chem.MolFromXYZFile(str(path)))
    except Exception as exc:
        return {"rdkit_ok": False, "error": str(exc)}


def split_sdf_records(text: str) -> list[str]:
    records: list[str] = []
    chunks = text.split("$$$$")
    for chunk in chunks:
        if chunk.strip():
            records.append(chunk.rstrip("\n") + "\n$$$$\n")
    return records


def sdf_record_summary(record_text: str) -> dict[str, object]:
    molblock = record_text.split("$$$$", 1)[0]
    try:
        return mol_summary(Chem.MolFromMolBlock(molblock, sanitize=True, removeHs=False))
    except Exception as exc:
        return {"rdkit_ok": False, "error": str(exc)}


def iter_sdf_records(path: Path) -> list[str]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            return split_sdf_records(handle.read())
    return split_sdf_records(path.read_text(encoding="utf-8", errors="replace"))


def smiles_from_line(raw: str) -> str | None:
    line = raw.replace("\u200b", "").strip()
    if not line or line.startswith("#"):
        return None
    return line.split()[0]


def safe_mol_from_smiles(smiles: str) -> dict[str, object]:
    try:
        return mol_summary(Chem.MolFromSmiles(smiles, sanitize=True))
    except Exception as exc:
        return {"rdkit_ok": False, "error": str(exc)}


def inventory_record(fixture_root: Path, path: Path) -> dict[str, object]:
    raw = path.read_bytes()
    return {
        "kind": "inventory",
        "fixture": path.relative_to(fixture_root).as_posix(),
        "byte_len": len(raw),
        "nonempty": len(raw) > 0,
    }


def build_records(fixture_root: Path) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for path in sorted(fixture_root.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(fixture_root).as_posix()
        lower = path.name.lower()
        if lower.endswith(".mol"):
            records.append(
                {
                    "kind": "mol",
                    "fixture": rel,
                    "record_index": 0,
                    **safe_mol_from_mol_file(path),
                }
            )
        elif lower.endswith(".mol2"):
            records.append(
                {
                    "kind": "mol2",
                    "fixture": rel,
                    "record_index": 0,
                    **safe_mol_from_mol2_file(path),
                }
            )
        elif lower.endswith(".xyz"):
            records.append(
                {
                    "kind": "xyz",
                    "fixture": rel,
                    "record_index": 0,
                    **safe_mol_from_xyz_file(path),
                }
            )
        elif lower.endswith(".sdf") or lower.endswith(".sdf.gz"):
            try:
                sdf_records = iter_sdf_records(path)
            except Exception as exc:
                records.append(
                    {
                        "kind": "sdf",
                        "fixture": rel,
                        "record_index": 0,
                        "rdkit_ok": False,
                        "error": str(exc),
                    }
                )
                continue
            if not sdf_records:
                records.append(
                    {
                        "kind": "sdf",
                        "fixture": rel,
                        "record_index": 0,
                        "rdkit_ok": False,
                        "error": "no SDF records",
                    }
                )
            for record_index, record_text in enumerate(sdf_records):
                records.append(
                    {
                        "kind": "sdf",
                        "fixture": rel,
                        "record_index": record_index,
                        **sdf_record_summary(record_text),
                    }
                )
        elif lower.endswith(".smi"):
            for line_index, raw in enumerate(
                path.read_text(encoding="utf-8", errors="replace").splitlines()
            ):
                smiles = smiles_from_line(raw)
                if smiles is None:
                    continue
                records.append(
                    {
                        "kind": "smi",
                        "fixture": rel,
                        "record_index": line_index,
                        "smiles": smiles,
                        **safe_mol_from_smiles(smiles),
                    }
                )
        else:
            records.append(inventory_record(fixture_root, path))
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_FIXTURE_ROOT,
        help=(
            "fixture root. When invoked by gen_all_rdkit_goldens.py this may be "
            f"the shared SMILES input file; in that case {DEFAULT_FIXTURE_ROOT} is used."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "tests" / "golden" / "rdkit_builtin_fixture_migration.jsonl",
        help="output JSONL path",
    )
    args = parser.parse_args()

    assert_rdkit_version()
    RDLogger.DisableLog("rdApp.*")
    fixture_root = args.input
    if fixture_root.is_file():
        fixture_root = DEFAULT_FIXTURE_ROOT
    records = build_records(fixture_root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, ensure_ascii=True, sort_keys=True))
            handle.write("\n")
    print(f"Wrote {len(records)} records to {args.output}")


if __name__ == "__main__":
    main()
