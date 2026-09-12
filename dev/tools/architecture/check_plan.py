"""Validate the sole split-crate migration plan and its unit manifest.

The checker is intentionally independent of Git.  It validates structural
execution rules and inventory evidence from the checked-out files only.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_PLAN = REPOSITORY_ROOT / "dev/plans/crate_architecture_completion_plan.md"
DEFAULT_MANIFEST = (
    REPOSITORY_ROOT / "dev/gap_reports/crate_migration/unit_manifest.json"
)

REQUIRED_READ_PATHS = (
    "dev/crate_architecture.md",
    "dev/agent_plan_standard.md",
    "dev/policy_invariants.md",
    "dev/source_reproduction_protocol.md",
    "dev/public_api_design.md",
    "dev/README.md",
)
ALLOWED_DISPOSITIONS = {"scheduled", "internal", "source_unsupported"}
STEP_RE = re.compile(r"^Step (?P<number>[1-9][0-9]*) \[(?P<checked>[ xX])\]: (?P<action>.+)$")
UNIT_RE = re.compile(r"^#### (?P<id>[A-Z][A-Za-z0-9_-]*) — (?P<scope>.+)$")
STAGE_RE = re.compile(r"^### (?P<id>[A-Z][A-Z0-9_]*) — ")
HISTORICAL_ROW_RE = re.compile(r"^\| `(?P<path>[^`]+)` \| .+\|$")


@dataclass(frozen=True)
class PlanStep:
    number: int
    checked: bool
    action: str
    line_number: int


@dataclass(frozen=True)
class ValidationIssue:
    code: str
    message: str

    def render(self) -> str:
        return f"{self.code}: {self.message}"


def _read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as error:
        raise ValueError(f"cannot read {path}: {error}") from error


def load_manifest(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(_read_text(path))
    except json.JSONDecodeError as error:
        raise ValueError(f"invalid JSON in {path}: {error}") from error
    if not isinstance(value, dict):
        raise ValueError(f"manifest root in {path} must be an object")
    return value


def parse_steps(plan_text: str) -> tuple[list[PlanStep], list[ValidationIssue]]:
    steps: list[PlanStep] = []
    issues: list[ValidationIssue] = []
    for line_number, line in enumerate(plan_text.splitlines(), start=1):
        if not line.startswith("Step "):
            continue
        match = STEP_RE.fullmatch(line)
        if match is None:
            issues.append(
                ValidationIssue(
                    "PLAN_STEP_FORMAT",
                    f"line {line_number} is not one complete canonical Step line",
                )
            )
            continue
        steps.append(
            PlanStep(
                number=int(match.group("number")),
                checked=match.group("checked").lower() == "x",
                action=match.group("action"),
                line_number=line_number,
            )
        )
    if not steps:
        issues.append(ValidationIssue("PLAN_NO_STEPS", "the plan has no Step lines"))
    return steps, issues


def parse_unit_headings(plan_text: str) -> list[str]:
    return [
        match.group("id")
        for line in plan_text.splitlines()
        if (match := UNIT_RE.fullmatch(line)) is not None
    ]


def parse_stage_headings(plan_text: str) -> list[str]:
    return [
        match.group("id")
        for line in plan_text.splitlines()
        if (match := STAGE_RE.match(line)) is not None
        and match.group("id") not in {"CLOSE", "BIND"}
    ]


def parse_historical_paths(plan_text: str) -> list[str]:
    marker = "## 附录：历史文件覆盖与最终 owner（不是第二执行队列）"
    if marker not in plan_text:
        return []
    appendix = plan_text.split(marker, maxsplit=1)[1]
    return [
        match.group("path")
        for line in appendix.splitlines()
        if (match := HISTORICAL_ROW_RE.fullmatch(line)) is not None
    ]


def validate_step_sequence(steps: list[PlanStep]) -> list[ValidationIssue]:
    issues: list[ValidationIssue] = []
    for expected, step in enumerate(steps, start=1):
        if step.number != expected:
            issues.append(
                ValidationIssue(
                    "PLAN_STEP_SEQUENCE",
                    f"line {step.line_number} has Step {step.number}; expected Step {expected}",
                )
            )
    return issues


def _is_read_step(step: PlanStep) -> bool:
    return step.action.startswith("Read ")


def _is_add_tests_step(step: PlanStep) -> bool:
    lowered = step.action.lower()
    return step.action.startswith("Add ") and ("test" in lowered or "测试" in step.action)


def _is_run_step(step: PlanStep) -> bool:
    return step.action.startswith("Run ")


def validate_read_prerequisites(steps: list[PlanStep]) -> list[ValidationIssue]:
    issues: list[ValidationIssue] = []
    for index, step in enumerate(steps):
        if _is_read_step(step):
            missing = [path for path in REQUIRED_READ_PATHS if f"`{path}`" not in step.action]
            if missing:
                issues.append(
                    ValidationIssue(
                        "PLAN_READ_CONTENT",
                        f"Step {step.number} omits required read paths: {', '.join(missing)}",
                    )
                )
            continue
        if index == 0 or not _is_read_step(steps[index - 1]):
            issues.append(
                ValidationIssue(
                    "PLAN_READ_PREDECESSOR",
                    f"real task Step {step.number} is not immediately preceded by a Read step",
                )
            )
    return issues


def validate_test_followups(steps: list[PlanStep]) -> list[ValidationIssue]:
    issues: list[ValidationIssue] = []
    for index, step in enumerate(steps):
        if not _is_add_tests_step(step):
            continue
        next_real = next(
            (candidate for candidate in steps[index + 1 :] if not _is_read_step(candidate)),
            None,
        )
        if next_real is None or not _is_run_step(next_real):
            suffix = "does not have a following real task" if next_real is None else (
                f"is followed by real task Step {next_real.number}, which is not Run"
            )
            issues.append(
                ValidationIssue("PLAN_TEST_RUN_ORDER", f"Add-tests Step {step.number} {suffix}")
            )
    return issues


def _require_keys(
    value: Any,
    keys: Iterable[str],
    location: str,
    issues: list[ValidationIssue],
) -> bool:
    if not isinstance(value, dict):
        issues.append(ValidationIssue("MANIFEST_TYPE", f"{location} must be an object"))
        return False
    missing = [key for key in keys if key not in value]
    if missing:
        issues.append(
            ValidationIssue(
                "MANIFEST_FIELDS",
                f"{location} is missing required fields: {', '.join(missing)}",
            )
        )
        return False
    return True


def validate_manifest_shape(manifest: dict[str, Any]) -> list[ValidationIssue]:
    issues: list[ValidationIssue] = []
    contract = manifest.get("contract")
    units = manifest.get("units")
    stages = manifest.get("stages")
    if not isinstance(contract, dict):
        issues.append(ValidationIssue("MANIFEST_CONTRACT", "contract must be an object"))
        return issues
    if not isinstance(stages, dict) or not stages:
        issues.append(ValidationIssue("MANIFEST_STAGES", "stages must be a non-empty object"))
    if not isinstance(units, list) or not units:
        issues.append(ValidationIssue("MANIFEST_UNITS", "units must be a non-empty array"))
        return issues

    required_units = contract.get("required_unit_fields", [])
    required_sources = contract.get("required_source_item_fields", [])
    seen: set[str] = set()
    for index, unit in enumerate(units):
        location = f"units[{index}]"
        if not _require_keys(unit, required_units, location, issues):
            continue
        unit_id = unit.get("id")
        if not isinstance(unit_id, str) or not unit_id:
            issues.append(ValidationIssue("MANIFEST_UNIT_ID", f"{location}.id is empty"))
            continue
        if unit_id in seen:
            issues.append(ValidationIssue("MANIFEST_UNIT_DUPLICATE", f"duplicate unit {unit_id}"))
        seen.add(unit_id)
        owner = unit.get("owner")
        if not isinstance(owner, str) or not owner.strip():
            issues.append(ValidationIssue("MANIFEST_OWNER", f"unit {unit_id} has no owner"))
        if unit.get("disposition") not in ALLOWED_DISPOSITIONS:
            issues.append(
                ValidationIssue(
                    "MANIFEST_DISPOSITION",
                    f"unit {unit_id} has invalid disposition {unit.get('disposition')!r}",
                )
            )
        source_items = unit.get("source_items")
        if not isinstance(source_items, list) or not source_items:
            issues.append(
                ValidationIssue("MANIFEST_SOURCE_ITEMS", f"unit {unit_id} has zero source coverage")
            )
        else:
            for source_index, source in enumerate(source_items):
                source_location = f"unit {unit_id} source_items[{source_index}]"
                if not _require_keys(source, required_sources, source_location, issues):
                    continue
                source_owner = source.get("owner")
                if not isinstance(source_owner, str) or not source_owner.strip():
                    issues.append(
                        ValidationIssue("MANIFEST_OWNER", f"{source_location} has no owner")
                    )
                if source.get("disposition") not in ALLOWED_DISPOSITIONS:
                    issues.append(
                        ValidationIssue(
                            "MANIFEST_DISPOSITION",
                            f"{source_location} has invalid disposition",
                        )
                    )
                if source.get("audit_state") == "frozen" and (
                    not source.get("source_path") or not source.get("source_symbol")
                ):
                    issues.append(
                        ValidationIssue(
                            "MANIFEST_FROZEN_SOURCE",
                            f"{source_location} is frozen without an exact path and symbol",
                        )
                    )
        _validate_test_layers(unit_id, unit.get("tests"), issues)
        _validate_evidence(unit_id, unit.get("evidence"), issues)
    return issues


def _validate_test_layers(
    unit_id: str, tests: Any, issues: list[ValidationIssue]
) -> None:
    if not isinstance(tests, dict):
        issues.append(ValidationIssue("MANIFEST_TESTS", f"unit {unit_id} tests must be an object"))
        return
    for layer in ("detached", "public"):
        config = tests.get(layer)
        if not isinstance(config, dict):
            issues.append(
                ValidationIssue("MANIFEST_TESTS", f"unit {unit_id} lacks {layer} test config")
            )
            continue
        if config.get("status") == "planned":
            if not config.get("package") or not config.get("target"):
                issues.append(
                    ValidationIssue(
                        "MANIFEST_TEST_TARGET",
                        f"unit {unit_id} planned {layer} tests lack package or target",
                    )
                )
            if config.get("require_nonzero_tests") is not True:
                issues.append(
                    ValidationIssue(
                        "MANIFEST_TEST_COUNT",
                        f"unit {unit_id} planned {layer} tests do not require a nonzero count",
                    )
                )
            if not isinstance(config.get("features"), list):
                issues.append(
                    ValidationIssue(
                        "MANIFEST_TEST_FEATURES",
                        f"unit {unit_id} planned {layer} features must be an array",
                    )
                )


def _validate_evidence(
    unit_id: str, evidence: Any, issues: list[ValidationIssue]
) -> None:
    required = ("unit_report", "source_inventory", "commands", "test_count", "completion")
    if not _require_keys(evidence, required, f"unit {unit_id} evidence", issues):
        return
    if evidence.get("completion") == "verified":
        if evidence.get("source_inventory") != "frozen":
            issues.append(
                ValidationIssue(
                    "MANIFEST_COMPLETION_EVIDENCE",
                    f"unit {unit_id} is verified without frozen source inventory",
                )
            )
        if not evidence.get("commands") or not isinstance(evidence.get("test_count"), int) or evidence.get("test_count", 0) <= 0:
            issues.append(
                ValidationIssue(
                    "MANIFEST_COMPLETION_EVIDENCE",
                    f"unit {unit_id} is verified without commands and a positive test count",
                )
            )


def _detect_cycle(graph: dict[str, list[str]]) -> list[str] | None:
    visiting: set[str] = set()
    visited: set[str] = set()
    stack: list[str] = []

    def visit(node: str) -> list[str] | None:
        if node in visited:
            return None
        if node in visiting:
            start = stack.index(node)
            return stack[start:] + [node]
        visiting.add(node)
        stack.append(node)
        for dependency in graph.get(node, []):
            cycle = visit(dependency)
            if cycle is not None:
                return cycle
        stack.pop()
        visiting.remove(node)
        visited.add(node)
        return None

    for node in graph:
        cycle = visit(node)
        if cycle is not None:
            return cycle
    return None


def validate_dependencies(manifest: dict[str, Any]) -> list[ValidationIssue]:
    issues: list[ValidationIssue] = []
    stages = manifest.get("stages")
    units = manifest.get("units")
    if not isinstance(stages, dict) or not isinstance(units, list):
        return issues

    stage_order = {stage: index for index, stage in enumerate(stages)}
    stage_graph: dict[str, list[str]] = {}
    for stage, config in stages.items():
        dependencies = config.get("depends_on_stages", []) if isinstance(config, dict) else []
        if not isinstance(dependencies, list):
            issues.append(
                ValidationIssue("MANIFEST_STAGE_DEPENDENCY", f"stage {stage} dependencies must be an array")
            )
            continue
        stage_graph[stage] = dependencies
        for dependency in dependencies:
            if dependency not in stages:
                issues.append(
                    ValidationIssue(
                        "MANIFEST_STAGE_DEPENDENCY",
                        f"stage {stage} references unknown dependency {dependency}",
                    )
                )
            elif stage_order[dependency] >= stage_order[stage]:
                issues.append(
                    ValidationIssue(
                        "MANIFEST_REVERSE_DEPENDENCY",
                        f"stage {stage} depends on non-prior stage {dependency}",
                    )
                )
    if cycle := _detect_cycle(stage_graph):
        issues.append(ValidationIssue("MANIFEST_DEPENDENCY_CYCLE", " -> ".join(cycle)))

    unit_order = {
        unit.get("id"): index
        for index, unit in enumerate(units)
        if isinstance(unit, dict) and isinstance(unit.get("id"), str)
    }
    unit_graph: dict[str, list[str]] = {}
    for unit in units:
        if not isinstance(unit, dict) or unit.get("id") not in unit_order:
            continue
        unit_id = unit["id"]
        dependencies = unit.get("dependencies", {})
        dependency_units = dependencies.get("units", []) if isinstance(dependencies, dict) else []
        unit_graph[unit_id] = dependency_units if isinstance(dependency_units, list) else []
        for dependency in unit_graph[unit_id]:
            if dependency not in unit_order:
                issues.append(
                    ValidationIssue(
                        "MANIFEST_UNIT_DEPENDENCY",
                        f"unit {unit_id} references unknown dependency {dependency}",
                    )
                )
            elif unit_order[dependency] >= unit_order[unit_id]:
                issues.append(
                    ValidationIssue(
                        "MANIFEST_REVERSE_DEPENDENCY",
                        f"unit {unit_id} depends on non-prior unit {dependency}",
                    )
                )
    if cycle := _detect_cycle(unit_graph):
        issues.append(ValidationIssue("MANIFEST_DEPENDENCY_CYCLE", " -> ".join(cycle)))
    return issues


def validate_plan_manifest_alignment(
    plan_text: str, manifest: dict[str, Any]
) -> list[ValidationIssue]:
    issues: list[ValidationIssue] = []
    plan_units = parse_unit_headings(plan_text)
    manifest_units = [
        unit.get("id") for unit in manifest.get("units", []) if isinstance(unit, dict)
    ]
    if plan_units != manifest_units:
        missing = [unit for unit in plan_units if unit not in manifest_units]
        extra = [unit for unit in manifest_units if unit not in plan_units]
        issues.append(
            ValidationIssue(
                "PLAN_MANIFEST_UNITS",
                f"unit order/coverage differs; missing={missing}, extra={extra}",
            )
        )
    plan_stages = parse_stage_headings(plan_text)
    manifest_stages = list(manifest.get("stages", {}))
    if plan_stages != manifest_stages:
        issues.append(
            ValidationIssue(
                "PLAN_MANIFEST_STAGES",
                f"stage order differs; plan={plan_stages}, manifest={manifest_stages}",
            )
        )
    return issues


def validate_historical_coverage(
    plan_text: str, manifest: dict[str, Any]
) -> list[ValidationIssue]:
    coverage = manifest.get("historical_coverage", {})
    expected_paths = coverage.get("expected_unique_paths") if isinstance(coverage, dict) else None
    expected_tests = coverage.get("expected_test_paths") if isinstance(coverage, dict) else None
    paths = parse_historical_paths(plan_text)
    unique_paths = set(paths)
    issues: list[ValidationIssue] = []
    if not paths:
        return [ValidationIssue("PLAN_HISTORICAL_COVERAGE", "historical appendix has zero paths")]
    if not isinstance(expected_paths, int) or len(unique_paths) != expected_paths:
        issues.append(
            ValidationIssue(
                "PLAN_HISTORICAL_COVERAGE",
                f"historical appendix has {len(unique_paths)} unique paths; expected {expected_paths}",
            )
        )
    test_paths = {path for path in unique_paths if path.startswith("tests/")}
    if not isinstance(expected_tests, int) or len(test_paths) != expected_tests:
        issues.append(
            ValidationIssue(
                "PLAN_HISTORICAL_TEST_COVERAGE",
                f"historical appendix has {len(test_paths)} unique test paths; expected {expected_tests}",
            )
        )
    return issues


def check_plan(plan_path: Path, manifest_path: Path) -> list[ValidationIssue]:
    plan_text = _read_text(plan_path)
    manifest = load_manifest(manifest_path)
    steps, issues = parse_steps(plan_text)
    issues.extend(validate_step_sequence(steps))
    issues.extend(validate_read_prerequisites(steps))
    issues.extend(validate_test_followups(steps))
    issues.extend(validate_manifest_shape(manifest))
    issues.extend(validate_dependencies(manifest))
    issues.extend(validate_plan_manifest_alignment(plan_text, manifest))
    issues.extend(validate_historical_coverage(plan_text, manifest))
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate the sequential split-crate plan and unit manifest."
    )
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--json", action="store_true", help="emit issues as JSON")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        issues = check_plan(args.plan, args.manifest)
    except ValueError as error:
        issues = [ValidationIssue("INPUT_ERROR", str(error))]
    if args.json:
        print(json.dumps([issue.__dict__ for issue in issues], indent=2))
    elif issues:
        for issue in issues:
            print(issue.render(), file=sys.stderr)
    else:
        print("Plan validation passed.")
    return 1 if issues else 0


if __name__ == "__main__":
    raise SystemExit(main())
