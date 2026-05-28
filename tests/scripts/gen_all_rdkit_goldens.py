#!/usr/bin/env python3
"""Generate every RDKit golden file used by the Rust parity tests."""

from __future__ import annotations

import argparse
import concurrent.futures
import os
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
    ("gen_rdkit_forcefield_params_golden.py", "forcefield_params.jsonl"),
    ("gen_rdkit_mmff_builtin_golden.py", "mmff_builtin.jsonl"),
    ("gen_rdkit_smiles_writer_golden.py", "smiles_writer.jsonl"),
    ("gen_rdkit_isomeric_smiles_golden.py", "isomeric_smiles.jsonl"),
    ("gen_rdkit_svg_golden.py", "svg_drawer.jsonl"),
    ("gen_rdkit_prepared_draw_golden.py", "prepared_draw_molecule.jsonl"),
    ("gen_rdkit_morgan_fingerprint_golden.py", "morgan_fingerprint.jsonl"),
    ("gen_rdkit_sdf_write_golden.py", "sdf_write.jsonl"),
    ("gen_rdkit_sdf_read_golden.py", "sdf_read.jsonl"),
    ("gen_rdkit_mol2_read_golden.py", "mol2_read.jsonl"),
    ("gen_rdkit_molfile_read_golden.py", "molfile_read.jsonl"),
    ("gen_rdkit_xyz_read_golden.py", "xyz_read.jsonl"),
    (
        "gen_rdkit_builtin_fixture_migration_golden.py",
        "rdkit_builtin_fixture_migration.jsonl",
    ),
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
    parser.add_argument(
        "--jobs",
        type=int,
        default=min(4, os.cpu_count() or 1, len(GENERATORS)),
        help="number of generator subprocesses to run in parallel (default: min(4, cpu_count, generator_count))",
    )
    return parser.parse_args()


def run_generator(
    python: str, script_name: str, input_path: Path, output_path: Path
) -> subprocess.CompletedProcess[str]:
    script_path = REPO_ROOT / "tests" / "scripts" / script_name
    cmd = [
        python,
        str(script_path),
        "--input",
        str(input_path),
        "--output",
        str(output_path),
    ]
    return subprocess.run(
        cmd,
        cwd=REPO_ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def print_completed(
    script_name: str, output_name: str, result: subprocess.CompletedProcess[str]
) -> None:
    print(f"[golden] {script_name} -> {output_name}")
    if result.stdout:
        print(result.stdout, end="" if result.stdout.endswith("\n") else "\n")
    if result.stderr:
        print(result.stderr, end="" if result.stderr.endswith("\n") else "\n", file=sys.stderr)
    if result.returncode != 0:
        raise subprocess.CalledProcessError(
            result.returncode,
            result.args,
            output=result.stdout,
            stderr=result.stderr,
        )


def main() -> None:
    args = parse_args()
    if args.jobs < 1:
        raise SystemExit("--jobs must be >= 1")

    args.golden_dir.mkdir(parents=True, exist_ok=True)
    if args.clean:
        for _, output_name in GENERATORS:
            output_path = args.golden_dir / output_name
            if output_path.exists():
                output_path.unlink()

    if args.jobs == 1:
        for script_name, output_name in GENERATORS:
            print(f"[golden] starting {script_name} -> {output_name}")
            result = run_generator(
                args.python, script_name, args.input, args.golden_dir / output_name
            )
            print_completed(script_name, output_name, result)
    else:
        print(f"[golden] running {len(GENERATORS)} generators with --jobs {args.jobs}")
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
            futures = {
                executor.submit(
                    run_generator,
                    args.python,
                    script_name,
                    args.input,
                    args.golden_dir / output_name,
                ): (script_name, output_name)
                for script_name, output_name in GENERATORS
            }
            for future in concurrent.futures.as_completed(futures):
                script_name, output_name = futures[future]
                result = future.result()
                print_completed(script_name, output_name, result)

    print(f"Generated {len(GENERATORS)} RDKit golden files in {args.golden_dir}")


if __name__ == "__main__":
    main()
