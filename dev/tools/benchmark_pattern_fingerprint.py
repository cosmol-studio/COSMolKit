#!/usr/bin/env python3
"""Deterministic COSMolKit/RDKit Pattern fingerprint benchmark.

The controller starts one fresh process per engine and case.  This keeps the
first call cold with respect to each implementation's compiled Pattern table,
then measures repeated calls in that same process to observe cache reuse.
Machine-specific JSON belongs under ``tmp/`` and is intentionally untracked.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
import platform
import resource
import statistics
import subprocess
import sys
import time
import tracemalloc
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
EXPECTED_RDKIT_VERSION = "2026.03.1"
CASES = (
    {
        "name": "path_rich",
        "kind": "smiles",
        "text": "C" * 64,
    },
    {
        "name": "ring_rich",
        "kind": "smiles",
        "text": "c1ccc2c(c1)ccc1c3ccccc3ccc21",
    },
    {
        "name": "query_bearing",
        "kind": "smarts",
        "text": "[$([#6]),$([#7])]~[!#1]",
    },
    {
        "name": "disconnected",
        "kind": "smiles",
        "text": "CCO.c1ccccc1.[Na+]",
    },
    {
        "name": "large",
        "kind": "smiles",
        "text": "C" * 160,
    },
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resident_bytes() -> int | None:
    status = Path("/proc/self/status")
    if not status.exists():
        return None
    for line in status.read_text(encoding="utf-8").splitlines():
        if line.startswith("VmRSS:"):
            return int(line.split()[1]) * 1024
    return None


def maximum_resident_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return int(value)
    return int(value) * 1024


def load_engine(engine: str, kind: str, text: str) -> tuple[Any, Any, str]:
    if engine == "cosmolkit":
        import cosmolkit

        molecule = (
            cosmolkit.parse_smarts(text)
            if kind == "smarts"
            else cosmolkit.Molecule.from_smiles(text)
        )

        def fingerprint() -> list[int]:
            return molecule.pattern_fingerprint(
                n_bits=2048, tautomeric=False
            ).on_bits()

        return molecule, fingerprint, cosmolkit.__version__

    from rdkit import Chem, rdBase

    if rdBase.rdkitVersion != EXPECTED_RDKIT_VERSION:
        raise RuntimeError(
            "RDKit version mismatch: "
            f"expected {EXPECTED_RDKIT_VERSION}, got {rdBase.rdkitVersion}"
        )
    molecule = (
        Chem.MolFromSmarts(text) if kind == "smarts" else Chem.MolFromSmiles(text)
    )
    if molecule is None:
        raise RuntimeError(f"RDKit rejected benchmark case {text!r}")

    def fingerprint() -> list[int]:
        return list(Chem.PatternFingerprint(molecule, fpSize=2048).GetOnBits())

    return molecule, fingerprint, rdBase.rdkitVersion


def worker(engine: str, case_name: str, iterations: int, rounds: int) -> dict[str, Any]:
    case = next((case for case in CASES if case["name"] == case_name), None)
    if case is None:
        raise ValueError(f"unknown benchmark case: {case_name}")

    molecule, fingerprint, version = load_engine(
        engine, str(case["kind"]), str(case["text"])
    )
    atom_count = int(molecule.num_atoms() if engine == "cosmolkit" else molecule.GetNumAtoms())
    bond_count = int(molecule.num_bonds() if engine == "cosmolkit" else molecule.GetNumBonds())

    gc.collect()
    rss_before = resident_bytes()
    tracemalloc.start()
    start = time.perf_counter_ns()
    expected = fingerprint()
    cold_ns = time.perf_counter_ns() - start
    traced_after_cold, traced_peak_after_cold = tracemalloc.get_traced_memory()
    rss_after_cold = resident_bytes()

    # The first call initializes the fixed compiled Pattern table.  These
    # untimed calls stabilize interpreter dispatch and allocator freelists;
    # they do not alter the immutable molecule or the compiled table.
    for _ in range(3):
        if fingerprint() != expected:
            raise RuntimeError("Pattern fingerprint changed during warmup")

    round_ns: list[int] = []
    for _ in range(rounds):
        start = time.perf_counter_ns()
        for _ in range(iterations):
            if fingerprint() != expected:
                raise RuntimeError("Pattern fingerprint changed during benchmark")
        round_ns.append(time.perf_counter_ns() - start)

    traced_current, traced_peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    rss_after = resident_bytes()
    warm_call_ns = [elapsed / iterations for elapsed in round_ns]

    return {
        "engine": engine,
        "engine_version": version,
        "case": case_name,
        "kind": case["kind"],
        "input_sha256": hashlib.sha256(str(case["text"]).encode()).hexdigest(),
        "atom_count": atom_count,
        "bond_count": bond_count,
        "n_bits": 2048,
        "on_bits": expected,
        "on_bit_count": len(expected),
        "iterations_per_round": iterations,
        "rounds": rounds,
        "cold_call_ns": cold_ns,
        "warm_call_ns": warm_call_ns,
        "warm_median_ns": statistics.median(warm_call_ns),
        "cold_to_warm_ratio": cold_ns / statistics.median(warm_call_ns),
        "allocation_observation": {
            "python_traced_current_after_cold_bytes": traced_after_cold,
            "python_traced_peak_after_cold_bytes": traced_peak_after_cold,
            "python_traced_current_final_bytes": traced_current,
            "python_traced_peak_final_bytes": traced_peak,
            "rss_before_bytes": rss_before,
            "rss_after_cold_bytes": rss_after_cold,
            "rss_final_bytes": rss_after,
            "maximum_rss_bytes": maximum_resident_bytes(),
        },
    }


def run_worker(engine: str, case: str, iterations: int, rounds: int) -> dict[str, Any]:
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--worker",
        "--engine",
        engine,
        "--case",
        case,
        "--iterations",
        str(iterations),
        "--rounds",
        str(rounds),
    ]
    completed = subprocess.run(
        command,
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return json.loads(completed.stdout)


def controller(output: Path, iterations: int, rounds: int) -> dict[str, Any]:
    measurements: list[dict[str, Any]] = []
    for case in CASES:
        for engine in ("cosmolkit", "rdkit"):
            measurements.append(
                run_worker(engine, str(case["name"]), iterations, rounds)
            )

    comparisons: list[dict[str, Any]] = []
    for case in CASES:
        name = str(case["name"])
        cosmolkit_result = next(
            row
            for row in measurements
            if row["case"] == name and row["engine"] == "cosmolkit"
        )
        rdkit_result = next(
            row
            for row in measurements
            if row["case"] == name and row["engine"] == "rdkit"
        )
        if cosmolkit_result["on_bits"] != rdkit_result["on_bits"]:
            raise RuntimeError(f"Pattern output mismatch in benchmark case {name}")
        comparisons.append(
            {
                "case": name,
                "exact_on_bits": True,
                "cosmolkit_to_rdkit_warm_ratio": cosmolkit_result["warm_median_ns"]
                / rdkit_result["warm_median_ns"],
                "cosmolkit_to_rdkit_cold_ratio": cosmolkit_result["cold_call_ns"]
                / rdkit_result["cold_call_ns"],
            }
        )

    script = Path(__file__).resolve()
    installed_extension = next(
        (
            path
            for path in (REPO_ROOT / ".venv").glob(
                "lib/python*/site-packages/cosmolkit/cosmolkit*.so"
            )
        ),
        None,
    )
    result = {
        "schema": "cosmolkit-pattern-fingerprint-benchmark-v1",
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "identity": {
            "python_executable": sys.executable,
            "python_version": platform.python_version(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "processor": platform.processor(),
            "cpu_count": os.cpu_count(),
            "script": str(script.relative_to(REPO_ROOT)),
            "script_sha256": sha256(script),
            "pattern_core_sha256": sha256(
                REPO_ROOT
                / "crates/cosmolkit-core/src/properties/fingerprint/pattern.rs"
            ),
            "installed_extension": (
                str(installed_extension) if installed_extension is not None else None
            ),
            "installed_extension_sha256": (
                sha256(installed_extension) if installed_extension is not None else None
            ),
        },
        "definition": {
            "iterations_per_round": iterations,
            "rounds": rounds,
            "cases": list(CASES),
            "memory_note": (
                "tracemalloc observes Python-managed allocations; RSS fields "
                "also include native Rust/C++ allocations and process noise"
            ),
        },
        "measurements": measurements,
        "comparisons": comparisons,
        "passed": all(row["exact_on_bits"] for row in comparisons),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    temporary.replace(output)
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "tmp/pattern-fingerprint-bench/benchmark.json",
    )
    parser.add_argument("--iterations", type=int, default=25)
    parser.add_argument("--rounds", type=int, default=7)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument(
        "--engine", choices=("cosmolkit", "rdkit"), help=argparse.SUPPRESS
    )
    parser.add_argument(
        "--case", choices=tuple(str(case["name"]) for case in CASES), help=argparse.SUPPRESS
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.iterations <= 0 or args.rounds <= 0:
        raise ValueError("iterations and rounds must be greater than zero")
    if args.worker:
        if args.engine is None or args.case is None:
            raise ValueError("worker mode requires --engine and --case")
        print(json.dumps(worker(args.engine, args.case, args.iterations, args.rounds)))
        return 0
    result = controller(args.output.resolve(), args.iterations, args.rounds)
    print(json.dumps({"output": str(args.output), "passed": result["passed"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
