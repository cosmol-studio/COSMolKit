#!/usr/bin/env python3
"""Benchmark completed high-feasibility descriptors against pinned RDKit.

Each engine/case/descriptor tuple runs in a fresh process. Exact output is
verified before timing. Cold measurements use molecules parsed before the
timer but never passed to the selected descriptor; warm measurements reuse one
descriptor-warmed molecule. Machine-specific JSON belongs under ``tmp/``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import statistics
import struct
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Callable


REPO_ROOT = Path(__file__).resolve().parents[2]
EXPECTED_RDKIT_VERSION = "2026.03.1"
EXPECTED_RDKIT_SOURCE_REVISION = "351f8f378f8ad6bbd517980c38896e66bf907af8"

CASES = (
    {"name": "path_16", "smiles": "C" * 16},
    {"name": "path_64", "smiles": "C" * 64},
    {"name": "fused_rings", "smiles": "c1ccc2c(c1)ccc1c3ccccc3ccc21"},
    {"name": "bridged_spiro", "smiles": "C1CC2CCC1C2.C1CCC2(CC1)CCCC2"},
    {"name": "druglike", "smiles": "CC(O)c1ccncc1C(=O)NCCCl"},
)

BENCHMARKS = (
    "chi_order_2",
    "chi_order_6",
    "ring_complexity",
    "mqn",
    "labute_asa",
    "slogp_vsa",
    "smr_vsa",
)


def f64_bits(value: float) -> str:
    return struct.pack(">d", float(value)).hex()


def value_bits(value: Any) -> Any:
    if isinstance(value, float):
        return f64_bits(value)
    if isinstance(value, (tuple, list)):
        return [value_bits(item) for item in value]
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_engine(
    engine: str, benchmark: str, smiles: str
) -> tuple[Callable[[], Any], Callable[[Any], Any], str]:
    if engine == "rdkit":
        from rdkit import Chem, rdBase
        from rdkit.Chem import rdMolDescriptors

        if rdBase.rdkitVersion != EXPECTED_RDKIT_VERSION:
            raise RuntimeError(
                "RDKit version mismatch: expected "
                f"{EXPECTED_RDKIT_VERSION}, got {rdBase.rdkitVersion}"
            )

        def parse() -> Any:
            molecule = Chem.MolFromSmiles(smiles)
            if molecule is None:
                raise RuntimeError(f"RDKit rejected benchmark SMILES {smiles!r}")
            return molecule

        calls: dict[str, Callable[[Any], Any]] = {
            "chi_order_2": lambda mol: rdMolDescriptors.CalcChiNv(mol, 2, False),
            "chi_order_6": lambda mol: rdMolDescriptors.CalcChiNv(mol, 6, False),
            "ring_complexity": lambda mol: (
                rdMolDescriptors.CalcNumRings(mol),
                rdMolDescriptors.CalcNumHeterocycles(mol),
                rdMolDescriptors.CalcNumSpiroAtoms(mol),
                rdMolDescriptors.CalcNumBridgeheadAtoms(mol),
            ),
            "mqn": lambda mol: list(rdMolDescriptors.MQNs_(mol, False)),
            "labute_asa": lambda mol: rdMolDescriptors.CalcLabuteASA(
                mol, True, False
            ),
            "slogp_vsa": lambda mol: list(
                rdMolDescriptors.SlogP_VSA_(mol, [], False)
            ),
            "smr_vsa": lambda mol: list(
                rdMolDescriptors.SMR_VSA_(mol, [], False)
            ),
        }
        return parse, calls[benchmark], rdBase.rdkitVersion

    import cosmolkit

    def parse() -> Any:
        return cosmolkit.Molecule.from_smiles(smiles)

    calls = {
        "chi_order_2": lambda mol: cosmolkit.calc_chi_nv(mol, 2, force=False),
        "chi_order_6": lambda mol: cosmolkit.calc_chi_nv(mol, 6, force=False),
        "ring_complexity": lambda mol: (
            cosmolkit.calc_num_rings(mol),
            cosmolkit.calc_num_heterocycles(mol),
            cosmolkit.calc_num_spiro_atoms(mol),
            cosmolkit.calc_num_bridgehead_atoms(mol),
        ),
        "mqn": lambda mol: cosmolkit.calc_mqns(mol),
        "labute_asa": lambda mol: cosmolkit.calc_labute_asa(
            mol, include_hydrogens=True, force=False
        ),
        "slogp_vsa": lambda mol: cosmolkit.calc_slogp_vsa(mol, force=False),
        "smr_vsa": lambda mol: cosmolkit.calc_smr_vsa(mol, force=False),
    }
    return parse, calls[benchmark], cosmolkit.__version__


def worker(
    engine: str,
    case_name: str,
    benchmark: str,
    iterations: int,
    samples: int,
) -> dict[str, Any]:
    case = next(item for item in CASES if item["name"] == case_name)
    parse, call, version = load_engine(engine, benchmark, str(case["smiles"]))
    expected = value_bits(call(parse()))

    cold_samples: list[float] = []
    for _ in range(samples):
        molecules = [parse() for _ in range(iterations)]
        started = time.perf_counter_ns()
        observed = [value_bits(call(molecule)) for molecule in molecules]
        elapsed = time.perf_counter_ns() - started
        if any(value != expected for value in observed):
            raise RuntimeError("cold descriptor output changed during benchmark")
        cold_samples.append(elapsed / iterations)

    warm_molecule = parse()
    if value_bits(call(warm_molecule)) != expected:
        raise RuntimeError("warmup descriptor output changed during benchmark")
    warm_samples: list[float] = []
    for _ in range(samples):
        started = time.perf_counter_ns()
        observed = [value_bits(call(warm_molecule)) for _ in range(iterations)]
        elapsed = time.perf_counter_ns() - started
        if any(value != expected for value in observed):
            raise RuntimeError("warm descriptor output changed during benchmark")
        warm_samples.append(elapsed / iterations)

    return {
        "engine": engine,
        "engine_version": version,
        "case": case_name,
        "benchmark": benchmark,
        "output": expected,
        "iterations_per_sample": iterations,
        "samples": samples,
        "cold_ns_per_call": cold_samples,
        "cold_median_ns": statistics.median(cold_samples),
        "warm_ns_per_call": warm_samples,
        "warm_median_ns": statistics.median(warm_samples),
    }


def run_worker(
    engine: str,
    case_name: str,
    benchmark: str,
    iterations: int,
    samples: int,
) -> dict[str, Any]:
    environment = os.environ.copy()
    environment.update(
        {
            "PYTHONHASHSEED": "0",
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
            "RAYON_NUM_THREADS": "1",
        }
    )
    completed = subprocess.run(
        [
            sys.executable,
            str(Path(__file__).resolve()),
            "--worker",
            engine,
            "--case",
            case_name,
            "--benchmark",
            benchmark,
            "--iterations",
            str(iterations),
            "--samples",
            str(samples),
        ],
        cwd=REPO_ROOT,
        env=environment,
        check=True,
        capture_output=True,
        text=True,
    )
    return json.loads(completed.stdout)


def case_complexity(case: dict[str, str]) -> dict[str, Any]:
    from rdkit import Chem

    molecule = Chem.MolFromSmiles(case["smiles"])
    if molecule is None:
        raise RuntimeError(f"RDKit rejected benchmark case {case['name']!r}")
    ring_info = molecule.GetRingInfo()
    return {
        "atoms": molecule.GetNumAtoms(),
        "bonds": molecule.GetNumBonds(),
        "atom_paths_by_length": {
            str(length): len(
                Chem.FindAllPathsOfLengthN(
                    molecule, length, useBonds=False, useHs=True
                )
            )
            for length in range(1, 8)
        },
        "atom_rings": len(ring_info.AtomRings()),
        "bond_rings": len(ring_info.BondRings()),
    }


def controller(output: Path, iterations: int, samples: int) -> dict[str, Any]:
    comparisons: list[dict[str, Any]] = []
    for case in CASES:
        for benchmark in BENCHMARKS:
            rdkit_result = run_worker(
                "rdkit", case["name"], benchmark, iterations, samples
            )
            cosmolkit_result = run_worker(
                "cosmolkit", case["name"], benchmark, iterations, samples
            )
            if rdkit_result["output"] != cosmolkit_result["output"]:
                raise RuntimeError(
                    "benchmark parity preflight failed for "
                    f"{case['name']}/{benchmark}"
                )
            comparisons.append(
                {
                    "case": case["name"],
                    "benchmark": benchmark,
                    "exact_output": True,
                    "rdkit": rdkit_result,
                    "cosmolkit": cosmolkit_result,
                    "cosmolkit_to_rdkit_cold_ratio": (
                        cosmolkit_result["cold_median_ns"]
                        / rdkit_result["cold_median_ns"]
                    ),
                    "cosmolkit_to_rdkit_warm_ratio": (
                        cosmolkit_result["warm_median_ns"]
                        / rdkit_result["warm_median_ns"]
                    ),
                }
            )

    script = Path(__file__).resolve()
    result = {
        "schema": "cosmolkit-high-feasibility-descriptor-benchmark-v1",
        "created_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "identity": {
            "python_executable": sys.executable,
            "python_version": platform.python_version(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "cpu_count": os.cpu_count(),
            "rdkit_version": EXPECTED_RDKIT_VERSION,
            "rdkit_source_revision": EXPECTED_RDKIT_SOURCE_REVISION,
            "script": str(script.relative_to(REPO_ROOT)),
            "script_sha256": sha256(script),
        },
        "definition": {
            "iterations_per_sample": iterations,
            "samples": samples,
            "cases": [
                {**case, "complexity": case_complexity(case)} for case in CASES
            ],
            "benchmarks": list(BENCHMARKS),
            "cold_definition": (
                "molecules are parsed before timing and have not been passed "
                "to the selected descriptor"
            ),
            "warm_definition": (
                "one molecule is descriptor-warmed before repeated timed calls"
            ),
        },
        "comparisons": comparisons,
        "passed": all(item["exact_output"] for item in comparisons),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    temporary.replace(output)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--iterations", type=int, default=50)
    parser.add_argument("--samples", type=int, default=7)
    parser.add_argument("--worker", choices=("rdkit", "cosmolkit"))
    parser.add_argument("--case", choices=tuple(item["name"] for item in CASES))
    parser.add_argument("--benchmark", choices=BENCHMARKS)
    args = parser.parse_args()
    if args.iterations < 1 or args.samples < 1:
        raise SystemExit("iterations and samples must be positive")
    if args.worker:
        if args.case is None or args.benchmark is None:
            raise SystemExit("--worker requires --case and --benchmark")
        print(
            json.dumps(
                worker(
                    args.worker,
                    args.case,
                    args.benchmark,
                    args.iterations,
                    args.samples,
                ),
                sort_keys=True,
            )
        )
        return
    if args.output is None:
        raise SystemExit("--output is required for the complete benchmark")
    result = controller(args.output, args.iterations, args.samples)
    print(json.dumps({"output": str(args.output), "passed": result["passed"]}))


if __name__ == "__main__":
    main()
