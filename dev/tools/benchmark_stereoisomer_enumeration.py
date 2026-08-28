#!/usr/bin/env python3
"""Deterministic COSMolKit/RDKit stereoisomer-enumeration benchmark.

The benchmark uses public APIs and starts a fresh process for each engine.  It
records observable stage boundaries instead of pretending that public API
timings can isolate private implementation instructions: candidate discovery,
configuration/count setup, one-configuration finalization, lazy-prefix,
exhaustive, bounded-random, embedding, and parallel end-to-end costs.
Every timed workload has an exact cross-engine output preflight.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Callable


REPO_ROOT = Path(__file__).resolve().parents[2]
EXPECTED_RDKIT_VERSION = "2026.03.1"
CASES = {
    "candidate": "CC=CC(C=CC)C(C)C(C=CC)C=CC",
    "single_finalize": "FC(Cl)Br",
    "full_exhaustive": "BrC=CC1OC(C2)(F)C2(Cl)C1",
    "bounded_random": "Br" + "[CH](Cl)" * 20 + "F",
    "embedding": "BrC=CC1OC(C2)(F)C2(Cl)C1",
}
PROFILES = (
    "candidate_discovery",
    "configuration_count",
    "single_configuration_finalization",
    "lazy_prefix",
    "full_exhaustive",
    "bounded_random",
    "embedding",
    "parallel",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def digest_json(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


def normalize_rd_stereo_info(info: Any) -> tuple[Any, ...]:
    stereo_type = str(info.type).split(".")[-1].lower()
    specified = str(info.specified).split(".")[-1].lower()
    descriptor = str(info.descriptor).split(".")[-1].lower()
    if descriptor == "novalue":
        descriptor = "none"
    centered_on = int(info.centeredOn)
    controlling = tuple(
        None if int(atom) == 4_294_967_295 else int(atom)
        for atom in info.controllingAtoms
    )
    return stereo_type, specified, centered_on, descriptor, int(info.permutation), controlling


def normalize_ck_stereo_info(info: Any) -> tuple[Any, ...]:
    return (
        info.stereo_type,
        info.specified,
        int(info.center_index),
        info.descriptor,
        int(info.permutation),
        tuple(info.controlling_atoms),
    )


def load_engine(engine: str) -> tuple[str, dict[str, Callable[[], Any]]]:
    if engine == "cosmolkit":
        import cosmolkit

        def molecule(name: str) -> Any:
            return cosmolkit.Molecule.from_smiles(CASES[name])

        candidate = molecule("candidate")
        single = molecule("single_finalize")
        full = molecule("full_exhaustive")
        random_case = molecule("bounded_random")
        embedding = molecule("embedding")

        def outputs(source: Any, **kwargs: Any) -> list[str]:
            options = cosmolkit.StereoisomerOptions(**kwargs)
            return [isomer.to_smiles() for isomer in source.stereoisomers(options)]

        functions: dict[str, Callable[[], Any]] = {
            "candidate_discovery": lambda: [
                normalize_ck_stereo_info(info)
                for info in candidate.analyze_potential_stereo().stereo_info
            ],
            "configuration_count": lambda: int(random_case.stereoisomer_count()),
            "single_configuration_finalization": lambda: outputs(
                single, max_isomers=1, unique=False
            ),
            "lazy_prefix": lambda: outputs(
                random_case, max_isomers=5, rand=0xF00D, unique=False
            ),
            "full_exhaustive": lambda: outputs(
                full, max_isomers=0, unique=False
            ),
            "bounded_random": lambda: outputs(
                random_case, max_isomers=16, rand=0xF00D, unique=False
            ),
            "embedding": lambda: outputs(
                embedding, try_embedding=True, max_isomers=16, unique=True
            ),
        }
        functions["parallel"] = lambda: parallel_workload(
            lambda text: outputs(
                cosmolkit.Molecule.from_smiles(text),
                max_isomers=4,
                rand=0xF00D,
                unique=False,
            )
        )
        return cosmolkit.__version__, functions

    from rdkit import Chem, rdBase
    from rdkit.Chem.EnumerateStereoisomers import (
        EnumerateStereoisomers,
        GetStereoisomerCount,
        StereoEnumerationOptions,
    )

    if rdBase.rdkitVersion != EXPECTED_RDKIT_VERSION:
        raise RuntimeError(
            f"RDKit version mismatch: expected {EXPECTED_RDKIT_VERSION}, "
            f"got {rdBase.rdkitVersion}"
        )

    def molecule(name: str) -> Any:
        result = Chem.MolFromSmiles(CASES[name])
        if result is None:
            raise RuntimeError(f"RDKit rejected benchmark case {name}")
        return result

    candidate = molecule("candidate")
    single = molecule("single_finalize")
    full = molecule("full_exhaustive")
    random_case = molecule("bounded_random")
    embedding = molecule("embedding")

    def outputs(source: Any, **kwargs: Any) -> list[str]:
        options = StereoEnumerationOptions(
            tryEmbedding=kwargs.get("try_embedding", False),
            onlyUnassigned=kwargs.get("only_unassigned", True),
            maxIsomers=kwargs.get("max_isomers", 1024),
            rand=kwargs.get("rand"),
            unique=kwargs.get("unique", True),
            onlyStereoGroups=kwargs.get("only_stereo_groups", False),
        )
        return [
            Chem.MolToSmiles(isomer, canonical=True, isomericSmiles=True)
            for isomer in EnumerateStereoisomers(source, options=options)
        ]

    functions = {
        "candidate_discovery": lambda: [
            normalize_rd_stereo_info(info)
            for info in Chem.FindPotentialStereo(Chem.Mol(candidate))
        ],
        "configuration_count": lambda: int(GetStereoisomerCount(random_case)),
        "single_configuration_finalization": lambda: outputs(
            single, max_isomers=1, unique=False
        ),
        "lazy_prefix": lambda: outputs(
            random_case, max_isomers=5, rand=0xF00D, unique=False
        ),
        "full_exhaustive": lambda: outputs(full, max_isomers=0, unique=False),
        "bounded_random": lambda: outputs(
            random_case, max_isomers=16, rand=0xF00D, unique=False
        ),
        "embedding": lambda: outputs(
            embedding, try_embedding=True, max_isomers=16, unique=True
        ),
    }
    functions["parallel"] = lambda: parallel_workload(
        lambda text: outputs(
            Chem.MolFromSmiles(text), max_isomers=4, rand=0xF00D, unique=False
        )
    )
    return rdBase.rdkitVersion, functions


def parallel_workload(enumerate_one: Callable[[str], list[str]]) -> list[list[str]]:
    inputs = (
        CASES["single_finalize"],
        "FC=CF",
        "CC(F)C(Cl)Br",
        "CC(F)C(C)C(C)F",
    ) * 8
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        return list(executor.map(enumerate_one, inputs))


def benchmark_profile(function: Callable[[], Any], rounds: int) -> dict[str, Any]:
    expected = function()
    for _ in range(2):
        if function() != expected:
            raise RuntimeError("benchmark workload is not deterministic")
    timings = []
    for _ in range(rounds):
        start = time.perf_counter_ns()
        actual = function()
        timings.append(time.perf_counter_ns() - start)
        if actual != expected:
            raise RuntimeError("benchmark output changed during timing")
    return {
        "output_digest": digest_json(expected),
        "output_items": len(expected) if isinstance(expected, list) else 1,
        "round_ns": timings,
        "median_ns": statistics.median(timings),
        "minimum_ns": min(timings),
    }


def worker(engine: str, rounds: int) -> dict[str, Any]:
    version, functions = load_engine(engine)
    measurements = {
        profile: benchmark_profile(functions[profile], rounds)
        for profile in PROFILES
    }
    return {
        "engine": engine,
        "engine_version": version,
        "rounds": rounds,
        "measurements": measurements,
    }


def run_worker(engine: str, rounds: int) -> dict[str, Any]:
    completed = subprocess.run(
        [
            sys.executable,
            str(Path(__file__).resolve()),
            "--worker",
            "--engine",
            engine,
            "--rounds",
            str(rounds),
        ],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return json.loads(completed.stdout)


def controller(output: Path, rounds: int) -> dict[str, Any]:
    measurements = [run_worker(engine, rounds) for engine in ("cosmolkit", "rdkit")]
    by_engine = {row["engine"]: row for row in measurements}
    comparisons = []
    for profile in PROFILES:
        cosmolkit_result = by_engine["cosmolkit"]["measurements"][profile]
        rdkit_result = by_engine["rdkit"]["measurements"][profile]
        exact = cosmolkit_result["output_digest"] == rdkit_result["output_digest"]
        if not exact:
            raise RuntimeError(f"benchmark parity preflight failed for {profile}")
        comparisons.append(
            {
                "profile": profile,
                "exact_output": True,
                "cosmolkit_to_rdkit_median_ratio": (
                    cosmolkit_result["median_ns"] / rdkit_result["median_ns"]
                ),
            }
        )

    script = Path(__file__).resolve()
    result = {
        "schema": "cosmolkit-stereoisomer-enumeration-benchmark-v1",
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "identity": {
            "python_executable": sys.executable,
            "python_version": platform.python_version(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "cpu_count": os.cpu_count(),
            "script": str(script.relative_to(REPO_ROOT)),
            "script_sha256": sha256(script),
            "stereo_core_sha256": sha256(
                REPO_ROOT
                / "crates/cosmolkit-core/src/chemistry/stereo_enumerate.rs"
            ),
            "potential_stereo_core_sha256": sha256(
                REPO_ROOT / "crates/cosmolkit-core/src/chemistry/potential_stereo.rs"
            ),
        },
        "definition": {
            "rounds": rounds,
            "profiles": list(PROFILES),
            "cases": CASES,
            "stage_boundary_note": (
                "Public API boundaries are measured intact. Configuration count includes "
                "candidate/flipper setup; one-configuration finalization includes setup, "
                "one source configuration, and source-defined finalization."
            ),
        },
        "measurements": measurements,
        "comparisons": comparisons,
        "passed": True,
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
        default=REPO_ROOT / "tmp/stereoisomer-enumeration-bench/benchmark.json",
    )
    parser.add_argument("--rounds", type=int, default=7)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument(
        "--engine", choices=("cosmolkit", "rdkit"), help=argparse.SUPPRESS
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.rounds <= 0:
        raise ValueError("--rounds must be greater than zero")
    if args.worker:
        if args.engine is None:
            raise ValueError("worker mode requires --engine")
        print(json.dumps(worker(args.engine, args.rounds)))
        return 0
    result = controller(args.output.resolve(), args.rounds)
    print(json.dumps({"output": str(args.output), "passed": result["passed"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
