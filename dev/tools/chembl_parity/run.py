#!/usr/bin/env python3
"""Run a checksummed, resumable ChEMBL parity profile."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
import tomllib
from collections import Counter
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from datetime import datetime, timezone
from importlib.metadata import version as distribution_version
from pathlib import Path
from typing import Any

import cosmolkit
from rdkit import rdBase


RUN_SCHEMA = "cosmolkit-chembl-parity-run-v1"
PROFILE_SCHEMA = "cosmolkit-chembl-parity-profile-v1"
CORPUS_SCHEMA = "cosmolkit-chembl-parity-corpus-v1"
KNOWN_SCRIPTS = {
    "audit_core.py",
    "audit_surfaces.py",
    "audit_combinations.py",
    "audit_fingerprints.py",
}
SCRIPT_MODULES = {
    script: f"dev.tools.chembl_parity.{Path(script).stem}"
    for script in KNOWN_SCRIPTS
}
SCRIPT_MODES = {
    "audit_core.py": {"core", "descriptors", "fingerprints"},
    "audit_surfaces.py": {
        "binary_roundtrip",
        "bounds",
        "conformer",
        "coordinates_2d",
        "explicit_h",
        "forcefield",
        "fragments",
        "inchi_options",
        "inchi_parse",
        "io",
        "smiles",
        "svg",
        "topology_operations",
    },
    "audit_combinations.py": {"batch", "concurrent", "scalar"},
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, value: Any) -> None:
    temporary = path.with_suffix(path.suffix + f".tmp-{os.getpid()}")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected a JSON object: {path}")
    return value


def workspace_version(root: Path) -> str:
    cargo = tomllib.loads((root / "Cargo.toml").read_text(encoding="utf-8"))
    return str(cargo["workspace"]["package"]["version"])


def extension_artifacts() -> list[Path]:
    module_path = Path(cosmolkit.__file__).resolve()
    candidates = sorted(
        path
        for pattern in ("cosmolkit*.so", "cosmolkit*.pyd", "cosmolkit*.dylib")
        for path in module_path.parent.glob(pattern)
        if path.is_file()
    )
    return candidates or [module_path]


def validate_profile(profile: dict[str, Any]) -> None:
    if profile.get("schema") != PROFILE_SCHEMA:
        raise ValueError(f"profile schema must be {PROFILE_SCHEMA!r}")
    if not isinstance(profile.get("name"), str) or not profile["name"]:
        raise ValueError("profile name must be a nonempty string")
    if not isinstance(profile.get("rdkit_version"), str):
        raise ValueError("profile rdkit_version must be a string")
    if not isinstance(profile.get("python_version"), str):
        raise ValueError("profile python_version must be a string")
    if not isinstance(profile.get("phases"), list) or not profile["phases"]:
        raise ValueError("profile phases must be a nonempty list")
    informational = profile.get("informational_mismatch_metrics", [])
    if not isinstance(informational, list) or not all(
        isinstance(value, str) and value.startswith("mismatch.")
        for value in informational
    ):
        raise ValueError("informational mismatch metrics must be exact mismatch keys")
    names: set[str] = set()
    for phase in profile["phases"]:
        if not isinstance(phase, dict):
            raise ValueError("every phase must be an object")
        name = phase.get("name")
        script = phase.get("script")
        if not isinstance(name, str) or not name or name in names:
            raise ValueError(f"phase name is empty or duplicated: {name!r}")
        names.add(name)
        if script not in KNOWN_SCRIPTS:
            raise ValueError(f"phase {name!r} uses unknown script {script!r}")
        for key in (
            "max_atoms",
            "records_per_shard",
            "worker_divisor",
            "batch_size",
            "expected_processed",
        ):
            value = phase.get(key)
            if value is not None and (not isinstance(value, int) or value < 1):
                raise ValueError(f"phase {name!r} has invalid {key}: {value!r}")
        if script != "audit_fingerprints.py":
            mode = phase.get("mode")
            if mode not in SCRIPT_MODES[script]:
                raise ValueError(
                    f"phase {name!r} has invalid mode {mode!r} for {script}"
                )
        if not isinstance(phase.get("expected_processed"), int):
            raise ValueError(f"phase {name!r} requires expected_processed")


def validate_corpus(
    corpus_dir: Path, profile: dict[str, Any], *, verify_files: bool
) -> tuple[dict[str, Any], str]:
    manifest_path = corpus_dir / "manifest.json"
    manifest = load_json(manifest_path)
    if manifest.get("schema") != CORPUS_SCHEMA:
        raise ValueError(f"corpus schema must be {CORPUS_SCHEMA!r}")
    expected_sha = profile.get("corpus_source_sha256")
    if expected_sha and manifest.get("source_sha256") != expected_sha:
        raise ValueError("corpus source checksum does not match the profile")
    expected_records = profile.get("corpus_records")
    if expected_records is not None and manifest.get("records") != expected_records:
        raise ValueError("corpus record count does not match the profile")
    files = manifest.get("files")
    if not isinstance(files, list) or not files:
        raise ValueError("corpus manifest has no shard files")
    shard_ids: list[int] = []
    records = 0
    for entry in files:
        if not isinstance(entry, dict):
            raise ValueError("invalid corpus shard entry")
        shard = entry.get("shard")
        file_name = entry.get("file")
        if not isinstance(shard, int) or not isinstance(file_name, str):
            raise ValueError("corpus shard and file fields have invalid types")
        if Path(file_name).name != file_name:
            raise ValueError(f"corpus shard path must be a file name: {file_name!r}")
        shard_ids.append(shard)
        records += int(entry.get("records", -1))
        if verify_files:
            path = corpus_dir / file_name
            if not path.is_file():
                raise ValueError(f"missing corpus shard: {path}")
            if path.stat().st_size != int(entry.get("bytes", -1)):
                raise ValueError(f"corpus shard size mismatch: {path}")
            if sha256_file(path) != entry.get("sha256"):
                raise ValueError(f"corpus shard checksum mismatch: {path}")
    if sorted(shard_ids) != list(range(len(files))):
        raise ValueError("corpus shard identifiers must be contiguous from zero")
    if records != manifest.get("records"):
        raise ValueError("corpus shard record counts do not match the manifest total")
    return manifest, sha256_file(manifest_path)


def tracked_diff(root: Path) -> bytes:
    return subprocess.check_output(["git", "diff", "--binary", "HEAD"], cwd=root)


def repository_status(root: Path) -> bytes:
    return subprocess.check_output(
        ["git", "status", "--porcelain=v1", "--untracked-files=all"], cwd=root
    )


def require_external_or_ignored(root: Path, path: Path, label: str) -> None:
    try:
        path.relative_to(root)
    except ValueError:
        return
    ignored = subprocess.run(
        ["git", "check-ignore", "--quiet", "--no-index", str(path)],
        cwd=root,
        input=b"",
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    if ignored.returncode != 0 and not path.exists():
        ignored = subprocess.run(
            ["git", "check-ignore", "--quiet", "--no-index", str(path.parent)],
            cwd=root,
            input=b"",
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
    if ignored.returncode != 0:
        raise ValueError(
            f"{label} is inside the repository but is not ignored by Git: {path}"
        )


def git_head(root: Path) -> str:
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=root, text=True
    ).strip()


def build_identity(
    root: Path,
    tool_dir: Path,
    profile_path: Path,
    profile: dict[str, Any],
    corpus_manifest: dict[str, Any],
    corpus_manifest_sha256: str,
    diff: bytes,
    status: bytes,
) -> dict[str, Any]:
    artifacts = extension_artifacts()
    scripts = sorted(KNOWN_SCRIPTS | {"run.py"})
    effective_profile = json.dumps(
        profile, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return {
        "git_head": git_head(root),
        "tracked_diff_sha256": sha256_bytes(diff),
        "repository_status_sha256": sha256_bytes(status),
        "profile": str(profile_path),
        "profile_sha256": sha256_file(profile_path),
        "effective_profile_sha256": sha256_bytes(effective_profile),
        "selected_phases": [str(phase["name"]) for phase in profile["phases"]],
        "corpus_manifest_sha256": corpus_manifest_sha256,
        "corpus_source_sha256": corpus_manifest["source_sha256"],
        "cosmolkit_version": cosmolkit.__version__,
        "cosmolkit_artifacts": [
            {"path": str(path), "sha256": sha256_file(path)} for path in artifacts
        ],
        "rdkit_version": rdBase.rdkitVersion,
        "python_version": sys.version,
        "python_implementation": sys.implementation.name,
        "python_cache_tag": sys.implementation.cache_tag,
        "python_executable": sys.executable,
        "numpy_version": distribution_version("numpy"),
        "scripts": {
            script: sha256_file(tool_dir / script) for script in scripts
        },
        "reference_profiles": {
            str(path.relative_to(root)): sha256_file(path)
            for path in (
                root
                / "tools/testdata/rdkit/rdkit_topological_fingerprint_profile.json",
                root / "tools/testdata/rdkit/avalon_fingerprint_profile.json",
            )
        },
    }


def write_repository_state(
    run_dir: Path,
    identity: dict[str, Any],
    diff: bytes,
    status: bytes,
    *,
    resume: bool,
) -> None:
    artifacts = {
        "tracked-diff.patch": (
            diff,
            str(identity["tracked_diff_sha256"]),
        ),
        "repository-status.txt": (
            status,
            str(identity["repository_status_sha256"]),
        ),
    }
    for name, (content, expected_sha) in artifacts.items():
        path = run_dir / name
        if resume:
            if not path.is_file() or sha256_file(path) != expected_sha:
                raise ValueError(f"cannot resume: repository-state artifact changed: {name}")
        else:
            path.write_bytes(content)


def validate_runtime(root: Path, profile: dict[str, Any]) -> str:
    expected_workspace_version = workspace_version(root)
    if cosmolkit.__version__ != expected_workspace_version:
        raise ValueError(
            "installed COSMolKit does not match the workspace: "
            f"installed {cosmolkit.__version__}, workspace {expected_workspace_version}"
        )
    if rdBase.rdkitVersion != profile["rdkit_version"]:
        raise ValueError(
            "RDKit version mismatch: expected "
            f"{profile['rdkit_version']}, got {rdBase.rdkitVersion}"
        )
    if sys.implementation.name != "cpython":
        raise ValueError("the complete parity profile requires CPython")
    actual_python = platform.python_version()
    if actual_python != profile["python_version"]:
        raise ValueError(
            "Python version mismatch: expected "
            f"{profile['python_version']}, got {actual_python}"
        )
    return expected_workspace_version


def command_for(
    root: Path,
    phase: dict[str, Any],
    shard_path: Path,
    output_path: Path,
    shard_records: int,
    max_examples: int,
) -> list[str]:
    script = str(phase["script"])
    command = [
        sys.executable,
        "-m",
        SCRIPT_MODULES[script],
        "--input",
        str(shard_path),
        "--output",
        str(output_path),
    ]
    if script == "audit_fingerprints.py":
        command.extend(
            (
                "--topological-profile",
                str(root / "tools/testdata/rdkit/rdkit_topological_fingerprint_profile.json"),
                "--avalon-profile",
                str(root / "tools/testdata/rdkit/avalon_fingerprint_profile.json"),
            )
        )
    else:
        command.extend(("--mode", str(phase["mode"])))
    limit = phase.get("records_per_shard")
    if limit is not None:
        limit = min(shard_records, int(limit))
    elif script in {"audit_surfaces.py"}:
        limit = shard_records
    if limit is not None:
        command.extend(("--limit", str(limit)))
    if phase.get("max_atoms") is not None:
        command.extend(("--max-atoms", str(phase["max_atoms"])))
    if phase.get("batch_size") is not None:
        command.extend(("--batch-size", str(phase["batch_size"])))
    command.extend(("--max-examples", str(max_examples)))
    return command


def run_task(
    root: Path,
    command: list[str],
    output_path: Path,
    log_path: Path,
    environment: dict[str, str],
    shard: int,
) -> dict[str, Any]:
    started = time.monotonic()
    started_at = utc_now()
    with log_path.open("w", encoding="utf-8") as log:
        log.write(json.dumps({"command": command, "started": started_at}) + "\n")
        log.flush()
        completed = subprocess.run(
            command,
            cwd=root,
            env=environment,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
            text=True,
        )
    result: dict[str, Any] = {
        "shard": shard,
        "command": command,
        "started": started_at,
        "finished": utc_now(),
        "elapsed_seconds": time.monotonic() - started,
        "returncode": completed.returncode,
        "output": str(output_path),
        "log": str(log_path),
    }
    if output_path.is_file():
        result["output_sha256"] = sha256_file(output_path)
    findings = output_path.with_suffix(".findings.jsonl")
    if findings.is_file():
        result["findings_sha256"] = sha256_file(findings)
    return result


def completed_task_is_valid(result: Any, command: list[str]) -> bool:
    if not isinstance(result, dict) or result.get("returncode") != 0:
        return False
    if result.get("command") != command:
        return False
    output = Path(str(result.get("output", "")))
    if not output.is_file() or sha256_file(output) != result.get("output_sha256"):
        return False
    findings = output.with_suffix(".findings.jsonl")
    if not findings.is_file() or sha256_file(findings) != result.get(
        "findings_sha256"
    ):
        return False
    try:
        summary = load_json(output)
    except (OSError, ValueError, json.JSONDecodeError):
        return False
    return isinstance(summary.get("processed"), int) and isinstance(
        summary.get("counts"), dict
    )


def aggregate_phase(
    phase_dir: Path,
    task_results: dict[str, Any],
    shard_count: int,
    informational_metrics: set[str],
    expected_processed: int,
) -> dict[str, Any]:
    counts: Counter[str] = Counter()
    processed = 0
    summaries = 0
    successful_tasks = 0
    retained_examples = 0
    for shard in range(shard_count):
        result = task_results.get(str(shard))
        if not isinstance(result, dict) or result.get("returncode") != 0:
            continue
        output = Path(str(result["output"]))
        if not output.is_file():
            continue
        summary = load_json(output)
        summaries += 1
        successful_tasks += 1
        processed += int(summary["processed"])
        counts.update({key: int(value) for key, value in summary["counts"].items()})
        findings = output.with_suffix(".findings.jsonl")
        if findings.is_file():
            with findings.open(encoding="utf-8") as source:
                retained_examples += sum(1 for _ in source)
    mismatch_counts = {
        key: value
        for key, value in sorted(counts.items())
        if key.startswith("mismatch.") and value
    }
    blocking = {
        key: value
        for key, value in mismatch_counts.items()
        if key not in informational_metrics
    }
    processed_as_expected = processed == expected_processed
    complete = (
        successful_tasks == shard_count
        and summaries == shard_count
        and processed_as_expected
    )
    aggregate = {
        "expected_tasks": shard_count,
        "successful_tasks": successful_tasks,
        "failed_or_missing_tasks": shard_count - successful_tasks,
        "summaries": summaries,
        "processed": processed,
        "expected_processed": expected_processed,
        "processed_as_expected": processed_as_expected,
        "counts": dict(sorted(counts.items())),
        "mismatch_counts": mismatch_counts,
        "informational_mismatch_counts": {
            key: value
            for key, value in mismatch_counts.items()
            if key in informational_metrics
        },
        "blocking_mismatch_counts": blocking,
        "retained_finding_examples": retained_examples,
        "complete": complete,
        "passed": complete and not blocking,
    }
    atomic_json(phase_dir / "aggregate.json", aggregate)
    return aggregate


def initialize_manifest(
    run_dir: Path,
    identity: dict[str, Any],
    profile: dict[str, Any],
    workers: int,
    resume: bool,
) -> dict[str, Any]:
    manifest_path = run_dir / "manifest.json"
    if run_dir.exists() and not resume:
        raise ValueError(f"run directory already exists; use --resume: {run_dir}")
    if resume:
        manifest = load_json(manifest_path)
        if manifest.get("schema") != RUN_SCHEMA:
            raise ValueError("run manifest schema mismatch")
        if manifest.get("identity") != identity:
            raise ValueError("cannot resume because the run identity changed")
        manifest["resumed"] = utc_now()
        manifest["workers"] = workers
        return manifest
    run_dir.mkdir(parents=True)
    return {
        "schema": RUN_SCHEMA,
        "profile_name": profile["name"],
        "started": utc_now(),
        "identity": identity,
        "platform": platform.platform(),
        "hostname": platform.node(),
        "logical_cpus": os.cpu_count(),
        "workers": workers,
        "phases": {},
    }


def main() -> None:
    default_root = Path(__file__).resolve().parents[3]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=default_root)
    parser.add_argument("--corpus-dir", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument(
        "--profile",
        type=Path,
        default=Path(__file__).resolve().parent / "profiles/complete.json",
    )
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    parser.add_argument("--max-examples", type=int, default=12)
    parser.add_argument("--phase", action="append", dest="selected_phases")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if args.workers < 1 or args.max_examples < 0:
        raise SystemExit("--workers must be positive and --max-examples nonnegative")

    root = args.root.resolve()
    tool_dir = Path(__file__).resolve().parent
    corpus_dir = args.corpus_dir.resolve()
    run_dir = args.run_dir.resolve()
    require_external_or_ignored(root, corpus_dir, "corpus directory")
    require_external_or_ignored(root, run_dir, "run directory")
    profile_path = args.profile.resolve()
    profile = load_json(profile_path)
    validate_profile(profile)
    if args.selected_phases:
        requested = set(args.selected_phases)
        known = {str(phase["name"]) for phase in profile["phases"]}
        unknown = requested - known
        if unknown:
            raise SystemExit(f"unknown phases: {', '.join(sorted(unknown))}")
        profile = dict(profile)
        profile["phases"] = [
            phase for phase in profile["phases"] if phase["name"] in requested
        ]
    installed_version = validate_runtime(root, profile)
    corpus_manifest, corpus_manifest_sha = validate_corpus(
        corpus_dir, profile, verify_files=True
    )
    diff = tracked_diff(root)
    status = repository_status(root)
    identity = build_identity(
        root,
        tool_dir,
        profile_path,
        profile,
        corpus_manifest,
        corpus_manifest_sha,
        diff,
        status,
    )
    manifest = initialize_manifest(run_dir, identity, profile, args.workers, args.resume)
    write_repository_state(
        run_dir, identity, diff, status, resume=args.resume
    )
    manifest_path = run_dir / "manifest.json"
    atomic_json(manifest_path, manifest)

    environment = os.environ.copy()
    environment.update(
        {
            "COSMOLKIT_EXPECTED_VERSION": installed_version,
            "PYTHONHASHSEED": "0",
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
            "RAYON_NUM_THREADS": "1",
        }
    )
    informational = set(profile.get("informational_mismatch_metrics", []))
    shards = corpus_manifest["files"]
    started = time.monotonic()
    for phase in profile["phases"]:
        phase_name = str(phase["name"])
        phase_dir = run_dir / phase_name
        phase_dir.mkdir(exist_ok=True)
        phase_state = manifest["phases"].setdefault(phase_name, {"tasks": {}})
        task_results = phase_state.setdefault("tasks", {})
        worker_count = min(
            len(shards), max(1, args.workers // int(phase.get("worker_divisor", 1)))
        )
        pending: list[tuple[dict[str, Any], list[str], Path, Path]] = []
        for shard_entry in shards:
            shard = int(shard_entry["shard"])
            output = phase_dir / f"shard-{shard:03d}.json"
            log = phase_dir / f"shard-{shard:03d}.log"
            command = command_for(
                root,
                phase,
                corpus_dir / shard_entry["file"],
                output,
                int(shard_entry["records"]),
                args.max_examples,
            )
            if not completed_task_is_valid(task_results.get(str(shard)), command):
                pending.append((shard_entry, command, output, log))
        print(
            json.dumps(
                {
                    "event": "phase_start",
                    "phase": phase_name,
                    "workers": worker_count,
                    "pending_tasks": len(pending),
                    "total_tasks": len(shards),
                },
                sort_keys=True,
            ),
            flush=True,
        )
        futures: dict[Any, int] = {}
        next_task = 0
        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            while next_task < len(pending) or futures:
                while next_task < len(pending) and len(futures) < worker_count:
                    entry, command, output, log = pending[next_task]
                    shard = int(entry["shard"])
                    future = executor.submit(
                        run_task,
                        root,
                        command,
                        output,
                        log,
                        environment,
                        shard,
                    )
                    futures[future] = shard
                    next_task += 1
                done, _ = wait(futures, return_when=FIRST_COMPLETED)
                for future in done:
                    shard = futures.pop(future)
                    result = future.result()
                    task_results[str(shard)] = result
                    phase_state["last_update"] = utc_now()
                    atomic_json(manifest_path, manifest)
                    print(
                        json.dumps(
                            {
                                "event": "task_finish",
                                "phase": phase_name,
                                "shard": shard,
                                "returncode": result["returncode"],
                                "elapsed_seconds": result["elapsed_seconds"],
                            },
                            sort_keys=True,
                        ),
                        flush=True,
                    )
        aggregate = aggregate_phase(
            phase_dir,
            task_results,
            len(shards),
            informational,
            int(phase["expected_processed"]),
        )
        phase_state["definition"] = phase
        phase_state["workers"] = worker_count
        phase_state["aggregate"] = aggregate
        phase_state["finished"] = utc_now()
        atomic_json(manifest_path, manifest)
        print(
            json.dumps(
                {"event": "phase_finish", "phase": phase_name, **aggregate},
                sort_keys=True,
            ),
            flush=True,
        )

    aggregates = [
        manifest["phases"][str(phase["name"])].get("aggregate", {})
        for phase in profile["phases"]
    ]
    manifest["finished"] = utc_now()
    manifest["elapsed_seconds"] = time.monotonic() - started
    manifest["complete"] = all(value.get("complete") is True for value in aggregates)
    manifest["passed"] = all(value.get("passed") is True for value in aggregates)
    atomic_json(manifest_path, manifest)
    print(
        json.dumps(
            {
                "event": "run_finish",
                "complete": manifest["complete"],
                "passed": manifest["passed"],
                "manifest": str(manifest_path),
            },
            sort_keys=True,
        ),
        flush=True,
    )
    if not manifest["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
