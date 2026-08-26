#!/usr/bin/env python3
"""Benchmark source-backed Layered fingerprints against pinned RDKit.

The benchmark validates exact output before timing. Native allocation pressure
is reported as isolated-worker resident-set and high-water-mark measurements;
it is deliberately not described as an allocator call count.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import statistics
import subprocess
import sys
import time
from importlib.metadata import version as distribution_version
from pathlib import Path
from typing import Any

import cosmolkit
from rdkit import Chem, rdBase


EXPECTED_RDKIT_VERSION = "2026.03.1"
EXPECTED_RDKIT_SOURCE_REVISION = "351f8f378f8ad6bbd517980c38896e66bf907af8"
PROFILE_PATH = (
    Path(__file__).resolve().parents[2]
    / "tools/testdata/rdkit/layered_fingerprint_profile.json"
)
CASES = (
    {
        "name": "path_rich",
        "kind": "smiles",
        "text": "CCCCCCCCCCCCCCCCCCCCCCCC",
    },
    {
        "name": "ring_rich",
        "kind": "smiles",
        "text": "c1ccc2c(c1)ccc1c3ccccc3ccc21",
    },
    {
        "name": "query_bearing",
        "kind": "smarts",
        "text": "[#6,#7]-[O,N;H1]-c1ccccc1",
    },
    {
        "name": "disconnected",
        "kind": "smiles",
        "text": "CC.CCCC.C1CCCCC1",
    },
    {
        "name": "large",
        "kind": "smiles",
        "text": "C" * 80,
    },
)


def load_profile() -> dict[str, Any]:
    profile = json.loads(PROFILE_PATH.read_text(encoding="utf-8"))
    if profile.get("rdkit_version") != "2026.3.1":
        raise ValueError("Layered benchmark profile RDKit version mismatch")
    if profile.get("rdkit_source_revision") != EXPECTED_RDKIT_SOURCE_REVISION:
        raise ValueError("Layered benchmark profile source revision mismatch")
    if profile.get("algorithm_version") != "0.7.0":
        raise ValueError("Layered benchmark profile algorithm version mismatch")
    return profile


def parse_case(case: dict[str, str]) -> tuple[Any, Any]:
    if case["kind"] == "smarts":
        rd_mol = Chem.MolFromSmarts(case["text"])
        ck_mol = cosmolkit.parse_smarts(case["text"])
    else:
        rd_mol = Chem.MolFromSmiles(case["text"])
        ck_mol = cosmolkit.Molecule.from_smiles(case["text"])
    if rd_mol is None:
        raise ValueError(f"RDKit rejected benchmark case {case['name']!r}")
    return rd_mol, ck_mol


def exact_observations(case: dict[str, str]) -> dict[str, Any]:
    rd_mol, ck_mol = parse_case(case)
    rd_counts = [0] * rd_mol.GetNumAtoms()
    rd_fp = Chem.LayeredFingerprint(rd_mol, atomCounts=rd_counts)
    ck_result = ck_mol.fingerprint_layered_with_output(
        atom_counts=[0] * rd_mol.GetNumAtoms()
    )
    ck_fp = ck_result.fingerprint()
    rd_value = {
        "n_bits": rd_fp.GetNumBits(),
        "on_bits": list(rd_fp.GetOnBits()),
        "atom_counts": rd_counts,
    }
    ck_value = {
        "n_bits": ck_fp.n_bits(),
        "on_bits": ck_fp.on_bits(),
        "atom_counts": ck_result.atom_counts(),
    }
    if rd_value != ck_value:
        raise ValueError(
            f"Layered benchmark parity preflight failed for {case['name']!r}"
        )
    return {
        "atoms": rd_mol.GetNumAtoms(),
        "bonds": rd_mol.GetNumBonds(),
        **rd_value,
    }


def enumeration_measurements(case: dict[str, str]) -> dict[str, Any]:
    rd_mol, _ = parse_case(case)
    linear: dict[str, int] = {}
    branched: dict[str, int] = {}
    for length in range(1, 8):
        linear[str(length)] = len(
            Chem.FindAllPathsOfLengthN(
                rd_mol, length, useBonds=True, useHs=True
            )
        )
        branched[str(length)] = len(
            Chem.FindAllSubgraphsOfLengthN(rd_mol, length, useHs=True)
        )
    return {
        "source": "pinned RDKit FindAllPathsOfLengthN/FindAllSubgraphsOfLengthN",
        "linear_by_length": linear,
        "branched_by_length": branched,
        "linear_total": sum(linear.values()),
        "branched_total": sum(branched.values()),
    }


def proc_memory_kib() -> tuple[int, int]:
    fields: dict[str, int] = {}
    for line in Path("/proc/self/status").read_text(encoding="utf-8").splitlines():
        if line.startswith(("VmRSS:", "VmHWM:")):
            name, value, unit = line.split()
            if unit != "kB":
                raise ValueError(f"unexpected /proc memory unit: {unit}")
            fields[name.rstrip(":")] = int(value)
    return fields["VmRSS"], fields["VmHWM"]


def benchmark_worker(
    engine: str,
    case_name: str,
    iterations: int,
    samples: int,
    warmup: int,
) -> dict[str, Any]:
    case = next(case for case in CASES if case["name"] == case_name)
    rd_mol, ck_mol = parse_case(case)

    def call() -> int:
        if engine == "rdkit":
            return Chem.LayeredFingerprint(rd_mol).GetNumBits()
        return ck_mol.fingerprint_layered().n_bits()

    checksum = 0
    for _ in range(warmup):
        checksum ^= call()
    rss_before, hwm_before = proc_memory_kib()
    elapsed_samples_ns: list[int] = []
    for _ in range(samples):
        started = time.perf_counter_ns()
        for _ in range(iterations):
            checksum ^= call()
        elapsed_samples_ns.append(time.perf_counter_ns() - started)
    rss_after, hwm_after = proc_memory_kib()
    per_call = [elapsed / iterations for elapsed in elapsed_samples_ns]
    return {
        "engine": engine,
        "case": case_name,
        "iterations_per_sample": iterations,
        "samples": samples,
        "warmup_iterations": warmup,
        "elapsed_samples_ns": elapsed_samples_ns,
        "median_ns_per_call": statistics.median(per_call),
        "minimum_ns_per_call": min(per_call),
        "rss_before_kib": rss_before,
        "rss_after_kib": rss_after,
        "hwm_before_kib": hwm_before,
        "hwm_after_kib": hwm_after,
        "peak_working_set_delta_kib": max(0, hwm_after - max(hwm_before, rss_before)),
        "checksum": checksum,
    }


def run_worker(
    engine: str,
    case_name: str,
    iterations: int,
    samples: int,
    warmup: int,
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
            "--iterations",
            str(iterations),
            "--samples",
            str(samples),
            "--warmup",
            str(warmup),
        ],
        check=True,
        capture_output=True,
        text=True,
        env=environment,
    )
    return json.loads(completed.stdout)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--iterations", type=int, default=100)
    parser.add_argument("--samples", type=int, default=7)
    parser.add_argument("--warmup", type=int, default=20)
    parser.add_argument("--worker", choices=("rdkit", "cosmolkit"))
    parser.add_argument("--case", choices=tuple(case["name"] for case in CASES))
    args = parser.parse_args()
    if min(args.iterations, args.samples, args.warmup) < 1:
        raise SystemExit("iterations, samples, and warmup must be positive")
    if rdBase.rdkitVersion != EXPECTED_RDKIT_VERSION:
        raise SystemExit(
            f"RDKit version mismatch: expected {EXPECTED_RDKIT_VERSION}, "
            f"got {rdBase.rdkitVersion}"
        )
    load_profile()
    if args.worker:
        if args.case is None:
            raise SystemExit("--worker requires --case")
        print(
            json.dumps(
                benchmark_worker(
                    args.worker,
                    args.case,
                    args.iterations,
                    args.samples,
                    args.warmup,
                ),
                sort_keys=True,
            )
        )
        return
    if args.output is None:
        raise SystemExit("--output is required for the complete benchmark")

    case_results: list[dict[str, Any]] = []
    for case in CASES:
        exact = exact_observations(case)
        enumeration = enumeration_measurements(case)
        rdkit_result = run_worker(
            "rdkit", case["name"], args.iterations, args.samples, args.warmup
        )
        cosmolkit_result = run_worker(
            "cosmolkit", case["name"], args.iterations, args.samples, args.warmup
        )
        case_results.append(
            {
                "case": case,
                "exact_observation": exact,
                "enumeration": enumeration,
                "rdkit": rdkit_result,
                "cosmolkit": cosmolkit_result,
                "cosmolkit_to_rdkit_median_ratio": (
                    cosmolkit_result["median_ns_per_call"]
                    / rdkit_result["median_ns_per_call"]
                ),
            }
        )

    output = {
        "schema": "cosmolkit-layered-fingerprint-benchmark-v1",
        "created_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "rdkit_version": rdBase.rdkitVersion,
        "rdkit_distribution_version": distribution_version("rdkit"),
        "rdkit_source_revision": EXPECTED_RDKIT_SOURCE_REVISION,
        "cosmolkit_version": cosmolkit.__version__,
        "python": sys.version,
        "platform": platform.platform(),
        "logical_cpus": os.cpu_count(),
        "profile_path": str(PROFILE_PATH),
        "profile_sha256": hashlib.sha256(PROFILE_PATH.read_bytes()).hexdigest(),
        "memory_metric": (
            "isolated worker /proc VmRSS and VmHWM; peak_working_set_delta_kib "
            "is a native allocation-pressure proxy, not an allocator-call count"
        ),
        "timing_boundary": (
            "pre-parsed molecule to explicit Layered fingerprint with default "
            "source parameters; one width read consumes each result"
        ),
        "cases": case_results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + f".tmp-{os.getpid()}")
    temporary.write_text(
        json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(args.output)
    print(json.dumps({"output": str(args.output), "cases": len(case_results)}))


if __name__ == "__main__":
    main()
