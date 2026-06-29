#!/usr/bin/env python3
"""Generate RDKit MACCS fingerprint golden data."""

from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path
from typing import Iterable

from rdkit import Chem, RDLogger
from rdkit.Chem import MACCSkeys

EXPECTED_RDKIT_VERSION = "2026.3.1"
RAW_BITS = 167
PUBLIC_BITS = RAW_BITS - 1

FIXTURES: list[tuple[str, str]] = [
    ("empty", ""),
    ("methane", "C"),
    ("key_002_rf", "[Rf]"),
    ("key_003_ge", "[Ge]"),
    ("key_004_ac", "[Ac]"),
    ("key_005_sc", "[Sc]"),
    ("key_006_la", "[La]"),
    ("key_007_v", "[V]"),
    ("key_009_fe", "[Fe]"),
    ("key_010_be", "[Be]"),
    ("key_012_cu", "[Cu]"),
    ("key_018_boron", "B"),
    ("key_020_silicon", "[SiH4]"),
    ("key_027_iodine", "I"),
    ("key_029_phosphorus", "P"),
    ("key_035_lithium", "[Li]"),
    ("fluorine_atom", "F"),
    ("sulfur_atom", "S"),
    ("chlorine_atom", "Cl"),
    ("salt_fragments", "CCO.Cl"),
    ("benzene", "c1ccccc1"),
    ("biphenyl", "c1ccccc1c2ccccc2"),
    ("pyridine", "c1ncccc1"),
    ("morpholine", "O1CCNCC1"),
    ("ammonium", "[NH4+]"),
    ("acetate", "CC(=O)[O-]"),
    ("isotopic_methane", "[13CH4]"),
    ("nitro", "N=O"),
    ("cyclopropanol", "C1CC1O"),
    ("fragment_methanes", "C.C"),
    ("all_key_low_mix", "NCCO"),
    ("all_key_high_mix", "OCOCOCO"),
    ("key_008_hetero_four_ring", "O1CCC1"),
    ("key_011_four_ring", "C1CCC1"),
    ("key_013_o_n_c_c", "ON(C)C"),
    ("key_014_disulfide", "CSSC"),
    ("key_015_o_c_o_o", "O=C(O)O"),
    ("key_016_hetero_three_ring", "O1CC1"),
    ("key_017_alkyne", "C#C"),
    ("key_019_seven_ring", "C1CCCCCC1"),
    ("key_021_alkene_dihetero", "C=C(O)O"),
    ("key_022_three_ring", "C1CC1"),
    ("key_023_n_c_o_o", "NC(=O)O"),
    ("key_025_n_c_n_n", "NC(N)N"),
    ("key_026_cyclic_alkene", "C1=C2CCCC2C1"),
    ("key_030_c_hetero_c_c_any", "C[S](C)(C)C"),
    ("key_031_hetero_halogen", "N[Pt](Cl)(Cl)N"),
    ("key_034_ch2_double", "C=C"),
    ("key_036_s_ring", "S1CC1"),
    ("key_037_n_c_o_n", "NC(=O)N"),
    ("key_038_n_c_c_n", "NC(C)N"),
    ("key_039_o_s_o_o", "COS(=O)(=O)O"),
    ("key_041_c_n_triple", "C#N"),
    ("key_044_exotic_element", "[SeH2]"),
    ("key_050_substituted_alkene", "CC(C)=C"),
    ("key_052_n_n", "NNO"),
    ("key_053_hetero_bridge_3", "NCCCO"),
    ("key_054_hetero_bridge_2", "NCCN"),
    ("key_056_o_n_o_c", "ON(O)C"),
    ("key_062_ring_nonring_ring", "C1CC1C1CC1"),
    ("key_063_n_o_double", "N=O"),
    ("key_066_quaternary_carbon", "CC(C)(C)C"),
    ("key_072_o_x_x_o", "OCCO"),
    ("key_075_exocyclic_n_ring", "CN1CC1"),
    ("key_080_n_x_x_x_n", "NCCCN"),
    ("key_083_hetero_five_ring", "O1CCCC1"),
    ("key_084_nh2", "N"),
    ("key_087_halogen_ring_chain", "FC1CC1"),
    ("key_089_o_bridge_3", "OCCCO"),
    ("key_092_o_c_n_c", "OC(N)C"),
    ("key_098_hetero_six_ring", "O1CCCCC1"),
    ("key_101_large_ring", "C1CCCCCCC1"),
    ("key_105_ring_branch_ring", "C1CC(C1)C1CC1"),
    ("key_106_hetero_three_neighbors", "N(O)(O)O"),
    ("key_107_halogen_three_neighbors", "FC(F)(F)F"),
    ("key_108_methyl_chain_ch2", "CCCC"),
    ("key_113_o_aromatic_chain", "Oc1ccccc1"),
    ("key_115_methyl_any_ch2_any", "CC(C)C"),
    ("key_121_n_ring", "N1CCCC1"),
    ("key_122_n_three_neighbors", "CN(C)C"),
    ("key_123_o_c_o", "COC"),
    ("key_130_two_hetero_pairs", "NOON"),
    ("key_132_o_ch2_chain", "OCCC"),
    ("key_133_ring_nonring_n", "C1CC1N"),
    ("key_135_n_aromatic_chain", "Nc1ccccc1"),
    ("key_136_two_o_double", "O=CC=O"),
    ("key_140_four_oxygens", "OCOCOCO"),
    ("key_141_three_methyl", "CC(C)(C)C"),
    ("key_142_two_nitrogens", "NN"),
    ("key_145_two_six_rings", "C1CCCCC1C1CCCCC1"),
    ("key_146_three_oxygens", "OCOCO"),
    ("key_150_nonring_ring_path", "C1CC1CC1CC1"),
    ("key_154_carbonyl", "C=O"),
]


def assert_rdkit_version() -> None:
    actual = version("rdkit")
    assert actual == EXPECTED_RDKIT_VERSION, (
        f"RDKit version mismatch: expected {EXPECTED_RDKIT_VERSION}, got {actual}"
    )


def iter_smiles(path: Path) -> Iterable[str]:
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        yield line


def maccs_record(record_type: str, label: str | None, smiles: str) -> dict[str, object]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "record_type": record_type,
            "label": label,
            "smiles": smiles,
            "rdkit_ok": False,
            "raw_n_bits": RAW_BITS,
            "public_n_bits": PUBLIC_BITS,
            "raw_on_bits": None,
            "public_on_bits": None,
            "error": "RDKit MolFromSmiles returned None",
        }

    fp = MACCSkeys.GenMACCSKeys(mol)
    bits = [idx for idx, char in enumerate(fp.ToBitString()) if char == "1"]
    public_bits = [bit - 1 for bit in bits if 1 <= bit < RAW_BITS]
    return {
        "record_type": record_type,
        "label": label,
        "smiles": smiles,
        "rdkit_ok": True,
        "raw_n_bits": RAW_BITS,
        "public_n_bits": PUBLIC_BITS,
        "raw_on_bits": bits,
        "public_on_bits": public_bits,
        "error": None,
    }


def build_records(smiles_rows: list[str]) -> list[dict[str, object]]:
    records = [
        maccs_record("fixture", label, smiles)
        for label, smiles in FIXTURES
    ]
    records.extend(
        maccs_record("corpus", None, smiles)
        for smiles in smiles_rows
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
        default=Path("tests/golden/maccs_fingerprint.jsonl"),
        help="output JSONL path (default: tests/golden/maccs_fingerprint.jsonl)",
    )
    args = parser.parse_args()

    assert_rdkit_version()
    RDLogger.DisableLog("rdApp.*")
    records = build_records(list(iter_smiles(args.input)))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as f:
        for record in records:
            f.write(json.dumps(record, ensure_ascii=True))
            f.write("\n")
    print(f"Wrote {len(records)} records to {args.output}")


if __name__ == "__main__":
    main()
