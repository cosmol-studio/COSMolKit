#!/usr/bin/env python3
"""Generate every RDKit golden file used by the Rust parity tests."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = REPO_ROOT / "tests" / "smiles.smi"
DEFAULT_GOLDEN_DIR = REPO_ROOT / "tests" / "golden"

GENERATORS: list[tuple[str, str]] = [
    ("gen_rdkit_v2000_minimal_golden.py", "molblock_v2000_minimal.jsonl"),
    ("gen_rdkit_kekulize_molblock_golden.py", "molblock_v2000_kekulized.jsonl"),
    ("gen_rdkit_kekulize_flags_golden.py", "kekulize_clear_flags_false.jsonl"),
    ("gen_rdkit_graph_features.py", "graph_features.jsonl"),
    ("gen_rdkit_tetrahedral_stereo_geometry.py", "tetrahedral_stereo_geometry.jsonl"),
    ("gen_rdkit_dg_bounds_golden.py", "dg_bounds_matrix.jsonl"),
    ("gen_rdkit_smiles_writer_golden.py", "smiles_writer.jsonl"),
    ("gen_rdkit_isomeric_smiles_golden.py", "isomeric_smiles.jsonl"),
    ("gen_rdkit_svg_golden.py", "svg_drawer.jsonl"),
    ("gen_rdkit_prepared_draw_golden.py", "prepared_draw_molecule.jsonl"),
    ("gen_rdkit_morgan_fingerprint_golden.py", "morgan_fingerprint.jsonl"),
    ("gen_rdkit_sdf_write_golden.py", "sdf_write.jsonl"),
    ("gen_rdkit_sdf_read_golden.py", "sdf_read.jsonl"),
    ("gen_rdkit_molfile_read_golden.py", "molfile_read.jsonl"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--python",
        default=sys.executable,
        help="Python interpreter used to run generator scripts (default: current interpreter)",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"input SMILES file (default: {DEFAULT_INPUT})",
    )
    parser.add_argument(
        "--golden-dir",
        type=Path,
        default=DEFAULT_GOLDEN_DIR,
        help=f"directory for generated golden JSONL files (default: {DEFAULT_GOLDEN_DIR})",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="delete existing target golden files before regenerating them",
    )
    return parser.parse_args()


def run_generator(python: str, script_name: str, input_path: Path, output_path: Path) -> None:
    script_path = REPO_ROOT / "tests" / "scripts" / script_name
    cmd = [
        python,
        str(script_path),
        "--input",
        str(input_path),
        "--output",
        str(output_path),
    ]
    print(f"[golden] {' '.join(cmd)}")
    subprocess.run(cmd, check=True, cwd=REPO_ROOT)


def main() -> None:
    args = parse_args()
    args.golden_dir.mkdir(parents=True, exist_ok=True)
    if args.clean:
        for _, output_name in GENERATORS:
            output_path = args.golden_dir / output_name
            if output_path.exists():
                output_path.unlink()

    for script_name, output_name in GENERATORS:
        run_generator(args.python, script_name, args.input, args.golden_dir / output_name)

    print(f"Generated {len(GENERATORS)} RDKit golden files in {args.golden_dir}")


if __name__ == "__main__":
    main()
