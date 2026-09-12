"""Run one manifest-selected migration test target and record its evidence."""

from __future__ import annotations

import argparse
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence


REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_MANIFEST = (
    REPOSITORY_ROOT / "dev/gap_reports/crate_migration/unit_manifest.json"
)
DEFAULT_EVIDENCE_DIRECTORY = (
    REPOSITORY_ROOT / "dev/gap_reports/crate_migration/validation"
)
TEST_RESULT_RE = re.compile(
    r"test result: (?:ok|FAILED)\. (?P<passed>\d+) passed; "
    r"(?P<failed>\d+) failed; (?P<ignored>\d+) ignored;"
)


class ConfigurationError(ValueError):
    """Raised when a manifest entry cannot safely select one test command."""


@dataclass(frozen=True)
class SelectedValidation:
    unit_id: str
    layer: str
    package: str
    target: str
    features: tuple[str, ...]
    profile: str
    require_nonzero_tests: bool
    require_runtime_strict: bool


def load_manifest(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except OSError as error:
        raise ConfigurationError(f"cannot read manifest {path}: {error}") from error
    except json.JSONDecodeError as error:
        raise ConfigurationError(f"invalid manifest JSON in {path}: {error}") from error
    if not isinstance(value, dict):
        raise ConfigurationError("manifest root must be an object")
    return value


def select_validation(
    manifest: dict[str, Any], unit_id: str, layer: str
) -> SelectedValidation:
    if layer not in {"detached", "public"}:
        raise ConfigurationError(f"unknown validation layer {layer!r}")
    matches = [
        unit
        for unit in manifest.get("units", [])
        if isinstance(unit, dict) and unit.get("id") == unit_id
    ]
    if len(matches) != 1:
        raise ConfigurationError(
            f"unit {unit_id!r} must resolve to exactly one manifest entry; found {len(matches)}"
        )
    unit = matches[0]
    tests = unit.get("tests")
    config = tests.get(layer) if isinstance(tests, dict) else None
    if not isinstance(config, dict):
        raise ConfigurationError(f"unit {unit_id} has no {layer} test configuration")
    if config.get("status") != "planned":
        raise ConfigurationError(
            f"unit {unit_id} layer {layer} is not planned: {config.get('status')!r}"
        )
    if layer == "public" and config.get("feature_audit_state") != "frozen":
        raise ConfigurationError(
            f"unit {unit_id} public features are not frozen by its audit"
        )
    package = config.get("package")
    target = config.get("target")
    features = config.get("features")
    profile = config.get("profile")
    if not isinstance(package, str) or not package:
        raise ConfigurationError(f"unit {unit_id} layer {layer} has no package")
    if not isinstance(target, str) or not target:
        raise ConfigurationError(f"unit {unit_id} layer {layer} has no test target")
    if not isinstance(features, list) or any(
        not isinstance(feature, str) or not feature for feature in features
    ):
        raise ConfigurationError(f"unit {unit_id} layer {layer} has invalid features")
    if profile not in {"debug", "release"}:
        raise ConfigurationError(
            f"unit {unit_id} layer {layer} has invalid profile {profile!r}"
        )
    require_nonzero = config.get("require_nonzero_tests")
    require_runtime_strict = config.get("require_runtime_strict")
    if require_nonzero is not True:
        raise ConfigurationError(
            f"unit {unit_id} layer {layer} does not require a nonzero test count"
        )
    if not isinstance(require_runtime_strict, bool):
        raise ConfigurationError(
            f"unit {unit_id} layer {layer} lacks an explicit runtime-strict policy"
        )
    if layer == "public" and not require_runtime_strict:
        raise ConfigurationError(
            f"unit {unit_id} public validation must require runtime strict checks"
        )
    return SelectedValidation(
        unit_id=unit_id,
        layer=layer,
        package=package,
        target=target,
        features=tuple(features),
        profile=profile,
        require_nonzero_tests=require_nonzero,
        require_runtime_strict=require_runtime_strict,
    )


def cargo_metadata_command() -> list[str]:
    return ["cargo", "metadata", "--locked", "--no-deps", "--format-version", "1"]


def cargo_test_command(selected: SelectedValidation) -> list[str]:
    command = ["cargo", "test", "-p", selected.package]
    if selected.profile == "release":
        command.append("--release")
    if selected.features:
        command.extend(["--features", ",".join(selected.features)])
    command.extend(["--test", selected.target])
    return command


def package_features(metadata: dict[str, Any], package: str) -> set[str]:
    packages = [
        item
        for item in metadata.get("packages", [])
        if isinstance(item, dict) and item.get("name") == package
    ]
    if len(packages) != 1:
        raise ConfigurationError(
            f"package {package!r} must resolve exactly once in Cargo metadata; found {len(packages)}"
        )
    features = packages[0].get("features")
    if not isinstance(features, dict):
        raise ConfigurationError(f"package {package!r} has invalid Cargo feature metadata")
    return set(features)


def validate_features(
    selected: SelectedValidation, available_features: set[str]
) -> None:
    unknown = sorted(set(selected.features) - available_features)
    if unknown:
        raise ConfigurationError(
            f"unit {selected.unit_id} selects unknown {selected.package} features: "
            + ", ".join(unknown)
        )
    if selected.require_runtime_strict:
        if selected.package != "cosmolkit":
            raise ConfigurationError(
                "runtime-strict public validation must target package cosmolkit"
            )
        if "op-contracts-strict" not in selected.features:
            raise ConfigurationError(
                f"unit {selected.unit_id} requires runtime strict checks but does not select "
                "cosmolkit feature op-contracts-strict"
            )
        if "op-contracts-strict" not in available_features:
            raise ConfigurationError(
                "cosmolkit does not define the required runtime feature op-contracts-strict"
            )


def parse_executed_test_count(output: str) -> int:
    return sum(
        int(match.group("passed")) + int(match.group("failed"))
        for match in TEST_RESULT_RE.finditer(output)
    )


def run_process(command: Sequence[str], repository_root: Path) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    environment["CARGO_TERM_COLOR"] = "never"
    return subprocess.run(
        list(command),
        cwd=repository_root,
        env=environment,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )


def _command_text(command: Sequence[str]) -> str:
    return shlex.join(command)


def write_evidence(path: Path, evidence: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(evidence, indent=2, sort_keys=True) + "\n"
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            delete=False,
        ) as temporary:
            temporary.write(payload)
            temporary_name = temporary.name
        Path(temporary_name).replace(path)
    finally:
        if temporary_name is not None:
            temporary_path = Path(temporary_name)
            if temporary_path.exists():
                temporary_path.unlink()


def evidence_path_for(
    evidence_directory: Path, selected: SelectedValidation
) -> Path:
    return evidence_directory / f"{selected.unit_id}-{selected.layer}.json"


def execute_validation(
    selected: SelectedValidation,
    repository_root: Path,
    evidence_path: Path,
    *,
    dry_run: bool = False,
) -> int:
    started_at = datetime.now(timezone.utc).isoformat()
    metadata_command = cargo_metadata_command()
    test_command = cargo_test_command(selected)
    evidence: dict[str, Any] = {
        "schema_version": 1,
        "unit": selected.unit_id,
        "layer": selected.layer,
        "package": selected.package,
        "target": selected.target,
        "features": list(selected.features),
        "require_runtime_strict": selected.require_runtime_strict,
        "started_at": started_at,
        "commands": [_command_text(metadata_command), _command_text(test_command)],
        "status": "started",
        "exit_code": None,
        "test_count": None,
    }
    if dry_run:
        evidence.update(status="dry_run", exit_code=0)
        write_evidence(evidence_path, evidence)
        print(_command_text(test_command))
        return 0

    metadata_result = run_process(metadata_command, repository_root)
    if metadata_result.returncode != 0:
        evidence.update(
            status="cargo_metadata_failed",
            exit_code=metadata_result.returncode,
            output=metadata_result.stdout,
            finished_at=datetime.now(timezone.utc).isoformat(),
        )
        write_evidence(evidence_path, evidence)
        sys.stdout.write(metadata_result.stdout)
        return metadata_result.returncode or 1
    try:
        metadata = json.loads(metadata_result.stdout)
        validate_features(selected, package_features(metadata, selected.package))
    except (json.JSONDecodeError, ConfigurationError) as error:
        evidence.update(
            status="configuration_failed",
            exit_code=2,
            error=str(error),
            finished_at=datetime.now(timezone.utc).isoformat(),
        )
        write_evidence(evidence_path, evidence)
        print(f"validation configuration error: {error}", file=sys.stderr)
        return 2

    print(_command_text(test_command))
    test_result = run_process(test_command, repository_root)
    sys.stdout.write(test_result.stdout)
    test_count = parse_executed_test_count(test_result.stdout)
    status = "passed"
    exit_code = test_result.returncode
    error: str | None = None
    if test_result.returncode != 0:
        status = "test_command_failed"
    elif selected.require_nonzero_tests and test_count == 0:
        status = "zero_tests_rejected"
        exit_code = 2
        error = "the selected target executed zero tests"
        print(f"validation error: {error}", file=sys.stderr)
    evidence.update(
        status=status,
        exit_code=exit_code,
        test_count=test_count,
        output=test_result.stdout,
        finished_at=datetime.now(timezone.utc).isoformat(),
    )
    if error is not None:
        evidence["error"] = error
    write_evidence(evidence_path, evidence)
    return exit_code


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run one manifest-selected migration validation target."
    )
    parser.add_argument("--unit", required=True)
    parser.add_argument("--layer", required=True, choices=("detached", "public"))
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--repository-root", type=Path, default=REPOSITORY_ROOT)
    parser.add_argument("--evidence-log", type=Path)
    parser.add_argument("--dry-run", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        manifest = load_manifest(args.manifest)
        selected = select_validation(manifest, args.unit, args.layer)
    except ConfigurationError as error:
        print(f"validation configuration error: {error}", file=sys.stderr)
        return 2
    evidence_path = args.evidence_log or evidence_path_for(
        DEFAULT_EVIDENCE_DIRECTORY, selected
    )
    return execute_validation(
        selected,
        args.repository_root,
        evidence_path,
        dry_run=args.dry_run,
    )


if __name__ == "__main__":
    raise SystemExit(main())
