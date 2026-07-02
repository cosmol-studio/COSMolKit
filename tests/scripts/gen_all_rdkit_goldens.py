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
PROFILE_INPUTS = {
    "smiles_small": REPO_ROOT / "tests" / "smiles.smi",
    "smiles_5000": REPO_ROOT / "tests" / "smiles_5000.smi",
}

PROFILE_GOLDEN_DIRS = {
    "smiles_small": REPO_ROOT / "tests" / "golden" / "smiles_small",
    "smiles_5000": REPO_ROOT / "tests" / "golden" / "smiles_5000",
}

GENERATOR_SPECS: list[dict[str, object]] = [
    {
        "script": "gen_rdkit_v2000_minimal_golden.py",
        "output": "molblock_v2000_minimal.jsonl",
        "suites": {"default", "strict-corpus", "all"},
    },
    {
        "script": "gen_rdkit_kekulize_molblock_golden.py",
        "output": "molblock_v2000_kekulized.jsonl",
        "suites": {"default", "strict-corpus", "all"},
    },
    {
        "script": "gen_rdkit_kekulize_flags_golden.py",
        "output": "kekulize_clear_flags_false.jsonl",
        "suites": {"default", "strict-corpus", "all"},
    },
    {
        "script": "gen_rdkit_graph_features.py",
        "output": "graph_features.jsonl",
        "suites": {"default", "strict-corpus", "all"},
    },
    {
        "script": "gen_rdkit_molecular_descriptors_golden.py",
        "output": "molecular_descriptors.jsonl",
        "suites": {"default", "strict-corpus", "all"},
    },
    {
        "script": "gen_rdkit_delete_substructs_golden.py",
        "output": "delete_substructs_onlyfrags_chirality.jsonl",
        "suites": {"default", "strict-corpus", "delete-substructs", "all"},
    },
    {
        "script": "gen_rdkit_tetrahedral_stereo_geometry.py",
        "output": "tetrahedral_stereo_geometry.jsonl",
        "suites": {"iterative", "all"},
    },
    {
        "script": "gen_rdkit_dg_bounds_golden.py",
        "output": "dg_bounds_matrix.jsonl",
        "suites": {"default", "strict-corpus", "all"},
    },
    {
        "script": "gen_rdkit_conformer_generation_golden.py",
        "output": "conformer_generation.jsonl",
        "suites": {"iterative", "all"},
    },
    {
        "script": "gen_rdkit_conformer_generation_library_golden.py",
        "output": "conformer_generation_library.jsonl",
        "suites": {"iterative", "all"},
    },
    {
        "script": "gen_rdkit_confseq_embed_golden.py",
        "output": "confseq_embed_template.jsonl",
        "suites": {"iterative", "all"},
    },
    {
        "script": "gen_rdkit_forcefield_params_golden.py",
        "output": "forcefield_params.jsonl",
        "suites": {"iterative", "all"},
    },
    {
        "script": "gen_rdkit_mmff_builtin_golden.py",
        "output": "mmff_builtin.jsonl",
        "suites": {"fixtures", "all"},
    },
    {
        "script": "gen_rdkit_smiles_writer_golden.py",
        "output": "smiles_writer.jsonl",
        "suites": {"default", "strict-corpus", "all"},
    },
    {
        "script": "gen_rdkit_isomeric_smiles_golden.py",
        "output": "isomeric_smiles.jsonl",
        "suites": {"default", "strict-corpus", "all"},
    },
    {
        "script": "gen_rdkit_svg_golden.py",
        "output": "svg_drawer.jsonl",
        "suites": {"default", "strict-corpus", "all"},
    },
    {
        "script": "gen_rdkit_prepared_draw_golden.py",
        "output": "prepared_draw_molecule.jsonl",
        "suites": {"default", "strict-corpus", "all"},
    },
    {
        "script": "gen_rdkit_morgan_fingerprint_golden.py",
        "output": "morgan_fingerprint.jsonl",
        "suites": {"fingerprint", "all"},
    },
    {
        "script": "gen_rdkit_maccs_fingerprint_golden.py",
        "output": "maccs_fingerprint.jsonl",
        "suites": {"fingerprint", "all"},
    },
    {
        "script": "gen_rdkit_sdf_write_golden.py",
        "output": "sdf_write.jsonl",
        "suites": {"default", "strict-corpus", "all"},
    },
    {
        "script": "gen_rdkit_sdf_read_golden.py",
        "output": "sdf_read.jsonl",
        "suites": {"default", "strict-corpus", "all"},
    },
    {
        "script": "gen_rdkit_mol2_read_golden.py",
        "output": "mol2_read.jsonl",
        "suites": {"fixtures", "all"},
    },
    {
        "script": "gen_rdkit_molfile_read_golden.py",
        "output": "molfile_read.jsonl",
        "suites": {"default", "strict-corpus", "all"},
    },
    {
        "script": "gen_rdkit_xyz_read_golden.py",
        "output": "xyz_read.jsonl",
        "suites": {"default", "strict-corpus", "all"},
    },
    {
        "script": "gen_rdkit_builtin_fixture_migration_golden.py",
        "output": "rdkit_builtin_fixture_migration.jsonl",
        "suites": {"fixtures", "all"},
    },
]

SUITE_ALIASES = {
    "small": "default",
    "daily": "default",
    "strict": "strict-corpus",
    "5000": "strict-corpus",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--python",
        default=sys.executable,
        help="Python interpreter used to run generator scripts (default: current interpreter)",
    )
    parser.add_argument(
        "--profile",
        choices=sorted(PROFILE_INPUTS),
        default="smiles_small",
        help="SMILES corpus/golden profile (default: smiles_small)",
    )
    parser.add_argument(
        "--suite",
        default="default",
        choices=[
            "default",
            "strict-corpus",
            "fixtures",
            "iterative",
            "fingerprint",
            "delete-substructs",
            "all",
            *sorted(SUITE_ALIASES),
        ],
        help="generator suite to run (default: default)",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="override the profile input SMILES file",
    )
    parser.add_argument(
        "--golden-dir",
        type=Path,
        default=None,
        help="override the profile golden directory",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="delete existing target golden files before regenerating them",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=None,
        help="number of generator subprocesses to run in parallel (default: min(4, cpu_count, generator_count))",
    )
    return parser.parse_args()


def selected_generators(suite: str) -> list[tuple[str, str]]:
    suite = SUITE_ALIASES.get(suite, suite)
    return [
        (spec["script"], spec["output"])
        for spec in GENERATOR_SPECS
        if suite in spec["suites"]
    ]


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
    generators = selected_generators(args.suite)
    input_path = args.input or PROFILE_INPUTS[args.profile]
    golden_dir = args.golden_dir or PROFILE_GOLDEN_DIRS[args.profile]
    jobs = args.jobs or min(4, os.cpu_count() or 1, len(generators))

    if not generators:
        raise SystemExit(f"--suite {args.suite!r} did not select any generators")
    if jobs < 1:
        raise SystemExit("--jobs must be >= 1")

    golden_dir.mkdir(parents=True, exist_ok=True)
    if args.clean:
        for _, output_name in generators:
            output_path = golden_dir / output_name
            if output_path.exists():
                output_path.unlink()

    if jobs == 1:
        for script_name, output_name in generators:
            print(f"[golden] starting {script_name} -> {output_name}")
            result = run_generator(
                args.python, script_name, input_path, golden_dir / output_name
            )
            print_completed(script_name, output_name, result)
    else:
        print(f"[golden] running {len(generators)} generators with --jobs {jobs}")
        with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
            futures = {
                executor.submit(
                    run_generator,
                    args.python,
                    script_name,
                    input_path,
                    golden_dir / output_name,
                ): (script_name, output_name)
                for script_name, output_name in generators
            }
            for future in concurrent.futures.as_completed(futures):
                script_name, output_name = futures[future]
                result = future.result()
                print_completed(script_name, output_name, result)

    print(
        f"Generated {len(generators)} RDKit golden files for profile {args.profile} "
        f"suite {SUITE_ALIASES.get(args.suite, args.suite)} in {golden_dir}"
    )


if __name__ == "__main__":
    main()
