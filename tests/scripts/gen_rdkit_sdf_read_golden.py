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
import io
import json
import re
import tempfile
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


def supplier_summary(mol: Chem.Mol) -> dict[str, object]:
    return {
        "atoms": atom_features(mol),
        "bonds": bond_features(mol),
        "positions": positions(mol),
        "chiral_tags": chiral_tags(mol),
        "smiles_out": smiles_pair(mol),
        "properties": {
            name: mol.GetProp(name)
            for name in mol.GetPropNames(includePrivate=True, includeComputed=True)
        },
        "atom_properties": [
            {
                name: atom.GetProp(name)
                for name in atom.GetPropNames(includePrivate=True, includeComputed=True)
            }
            for atom in mol.GetAtoms()
        ],
        "bond_properties": [
            {
                name: bond.GetProp(name)
                for name in bond.GetPropNames(includePrivate=True, includeComputed=True)
            }
            for bond in mol.GetBonds()
        ],
    }


def error_row(base: dict[str, object], sdf: str | None, error: str | None) -> dict[str, object]:
    return {
        **base,
        "rdkit_ok": False,
        "sdf": sdf,
        "error": error or "RDKit operation failed",
    }


def success_row(base: dict[str, object], sdf: str, mol: Chem.Mol) -> dict[str, object]:
    return {
        **base,
        "rdkit_ok": True,
        "sdf": sdf,
        **supplier_summary(mol),
        "error": None,
    }


def read_forward_supplier(
    sdf: str,
    *,
    sanitize: bool,
    remove_hs: bool,
    strict_parsing: bool,
    process_property_lists: bool,
) -> tuple[Chem.Mol | None, str | None]:
    try:
        supplier = Chem.ForwardSDMolSupplier(
            io.BytesIO(sdf.encode("utf-8")),
            sanitize=sanitize,
            removeHs=remove_hs,
            strictParsing=strict_parsing,
        )
        supplier.SetProcessPropertyLists(process_property_lists)
        parsed = next(supplier)
    except StopIteration:
        return None, "ForwardSDMolSupplier returned no records"
    except Exception as exc:
        return None, str(exc)
    if parsed is None:
        return None, "ForwardSDMolSupplier returned None"
    return parsed, None


def read_indexed_supplier(
    sdf: str,
    *,
    sanitize: bool,
    remove_hs: bool,
    strict_parsing: bool,
    process_property_lists: bool,
) -> tuple[Chem.Mol | None, str | None]:
    with tempfile.NamedTemporaryFile("w", suffix=".sdf", encoding="utf-8") as handle:
        handle.write(sdf)
        handle.flush()
        try:
            supplier = Chem.SDMolSupplier(
                handle.name,
                sanitize=sanitize,
                removeHs=remove_hs,
                strictParsing=strict_parsing,
            )
            supplier.SetProcessPropertyLists(process_property_lists)
            parsed = supplier[0]
        except Exception as exc:
            return None, str(exc)
    if parsed is None:
        return None, "SDMolSupplier returned None"
    return parsed, None


def append_data_fields(block: str) -> str:
    return (
        block
        + ">  <ID>\nrecord-1\n\n"
        + ">  <BLANK_VALUE>\n\n"
        + ">  <REPEATED>\nfirst\n\n"
        + ">  <REPEATED>\nsecond\n\n"
        + ">  <atom.prop.note>\n[missing] alpha beta\n\n"
        + ">  <atom.iprop.rank>\n10 20 30 40 50 60 70 80 90 100\n\n"
        + ">  <bond.prop.flag>\none two three four five six seven eight nine\n\n"
        + "$$$$\n"
    )


def malformed_sdf_records() -> list[tuple[str, str]]:
    valid_minimal = (
        "malformed\n"
        "  COSMolKit      2D\n"
        "\n"
        "  1  0  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "M  END\n"
    )
    return [
        ("empty", ""),
        ("missing_delimiter", valid_minimal),
        ("invalid_data_header", valid_minimal + "> no label\nvalue\n\n$$$$\n"),
        ("spurious_data", valid_minimal + "spurious\n\n$$$$\n"),
        ("missing_m_end", valid_minimal.replace("M  END\n", "") + "$$$$\n"),
    ]


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
        "api": "ForwardSDMolSupplier",
        "operation": "read",
        "sanitize": True,
        "remove_hs": False,
        "strict_parsing": True,
        "process_property_lists": True,
    }
    if source is None:
        return error_row(base, None, source_error or "source molecule generation failed")

    block, render_error = render_molblock(source, fmt)
    if block is None:
        return error_row(base, None, render_error)
    if markers == "coords_only":
        block = strip_stereo_markers(block, fmt)

    parsed, parse_error = parse_molblock(block)
    if parsed is None:
        return error_row(base, block + "$$$$\n", parse_error)

    return success_row(base, block + "$$$$\n", parsed)


def build_supplier_cases(
    smiles: str,
    sdf: str,
    *,
    dimension: str,
    fmt: str,
    markers: str,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    sdf_variants = [
        ("plain", sdf),
        ("data_fields", append_data_fields(sdf.removesuffix("$$$$\n"))),
    ]
    for variant, sdf_text in sdf_variants:
        for api_name, parser_fn in [
            ("ForwardSDMolSupplier", read_forward_supplier),
            ("SDMolSupplier", read_indexed_supplier),
        ]:
            for sanitize in [True, False]:
                for remove_hs in [True, False]:
                    for strict_parsing in [True, False]:
                        for process_property_lists in [True, False]:
                            case_id = (
                                f"{dimension.lower()}_{fmt.lower()}_{markers}_{variant}_"
                                f"{api_name.lower()}_sanitize_{int(sanitize)}_"
                                f"removehs_{int(remove_hs)}_strict_{int(strict_parsing)}_"
                                f"proplists_{int(process_property_lists)}"
                            )
                            base = {
                                "smiles": smiles,
                                "case_id": case_id,
                                "dimension": dimension,
                                "format": fmt,
                                "stereo_markers": markers,
                                "api": api_name,
                                "operation": "read",
                                "sanitize": sanitize,
                                "remove_hs": remove_hs,
                                "strict_parsing": strict_parsing,
                                "process_property_lists": process_property_lists,
                            }
                            parsed, parse_error = parser_fn(
                                sdf_text,
                                sanitize=sanitize,
                                remove_hs=remove_hs,
                                strict_parsing=strict_parsing,
                                process_property_lists=process_property_lists,
                            )
                            rows.append(
                                error_row(base, sdf_text, parse_error)
                                if parsed is None
                                else success_row(base, sdf_text, parsed)
                            )
    return rows


def build_malformed_cases() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for label, sdf_text in malformed_sdf_records():
        for api_name, parser_fn in [
            ("ForwardSDMolSupplier", read_forward_supplier),
            ("SDMolSupplier", read_indexed_supplier),
        ]:
            for strict_parsing in [True, False]:
                base = {
                    "smiles": None,
                    "case_id": f"malformed_{label}_{api_name.lower()}_strict_{int(strict_parsing)}",
                    "dimension": "invalid",
                    "format": "invalid",
                    "stereo_markers": "invalid",
                    "api": api_name,
                    "operation": "malformed",
                    "sanitize": True,
                    "remove_hs": True,
                    "strict_parsing": strict_parsing,
                    "process_property_lists": True,
                }
                parsed, parse_error = parser_fn(
                    sdf_text,
                    sanitize=True,
                    remove_hs=True,
                    strict_parsing=strict_parsing,
                    process_property_lists=True,
                )
                rows.append(
                    error_row(base, sdf_text, parse_error)
                    if parsed is None
                    else success_row(base, sdf_text, parsed)
                )
    return rows


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
                sdf = records[-1].get("sdf")
                if isinstance(sdf, str):
                    records.extend(
                        build_supplier_cases(
                            smiles,
                            sdf,
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
    rows.extend(build_malformed_cases())

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for record in rows:
            handle.write(json.dumps(record, ensure_ascii=True))
            handle.write("\n")

    print(f"Wrote {len(rows)} records to {args.output}")


if __name__ == "__main__":
    main()
