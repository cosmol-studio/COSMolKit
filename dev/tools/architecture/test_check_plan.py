from __future__ import annotations

import copy
import json
import tempfile
import unittest
from pathlib import Path

from check_plan import check_plan, validate_dependencies


READ_ACTION = (
    "Read `dev/crate_architecture.md`, `dev/agent_plan_standard.md`, "
    "`dev/policy_invariants.md`, `dev/source_reproduction_protocol.md`, "
    "`dev/public_api_design.md` and `dev/README.md` to reload the rules."
)


def valid_plan() -> str:
    return f"""# Test Plan

### VAL — cosmolkit-types

#### VAL-vocab — vocabulary

Step 1 [ ]: {READ_ACTION}
Step 2 [ ]: Audit VAL and write its report.

## 附录：历史文件覆盖与最终 owner（不是第二执行队列）

| 历史文件 | 最终 owner / 阶段 |
|---|---|
| `src/a.rs` | types / VAL |
| `tests/a.rs` | types / VAL |
"""


def valid_manifest() -> dict[str, object]:
    return {
        "schema_version": 1,
        "contract": {
            "required_unit_fields": [
                "id",
                "stage",
                "sequence",
                "owner",
                "target_root",
                "disposition",
                "source_items",
                "dependencies",
                "naming",
                "tests",
                "evidence",
            ],
            "required_source_item_fields": [
                "source_kind",
                "source_path",
                "source_symbol",
                "scope",
                "disposition",
                "owner",
                "audit_state",
            ],
        },
        "historical_coverage": {
            "expected_unique_paths": 2,
            "expected_test_paths": 1,
        },
        "stages": {
            "VAL": {
                "owner": "cosmolkit-types",
                "target_root": "crates/cosmolkit-types/src",
                "depends_on_stages": [],
            }
        },
        "units": [
            {
                "id": "VAL-vocab",
                "stage": "VAL",
                "sequence": 1,
                "owner": "cosmolkit-types",
                "target_root": "crates/cosmolkit-types/src",
                "disposition": "scheduled",
                "source_items": [
                    {
                        "source_kind": "historical_and_pinned_upstream",
                        "source_path": None,
                        "source_symbol": None,
                        "scope": "vocabulary",
                        "disposition": "scheduled",
                        "owner": "cosmolkit-types",
                        "audit_state": "audit_required",
                    }
                ],
                "dependencies": {"stages": [], "units": []},
                "naming": {},
                "tests": {
                    "detached": {
                        "status": "planned",
                        "package": "cosmolkit-types",
                        "target": "migration_val_vocab",
                        "features": [],
                        "require_nonzero_tests": True,
                    },
                    "public": {
                        "status": "not_scheduled_for_this_unit",
                        "package": None,
                        "target": None,
                        "features": [],
                        "require_nonzero_tests": False,
                    },
                },
                "evidence": {
                    "unit_report": "report.md",
                    "source_inventory": "audit_required",
                    "commands": [],
                    "test_count": None,
                    "completion": "not_verified",
                },
            }
        ],
    }


class CheckPlanNegativeTests(unittest.TestCase):
    def run_check(
        self, plan: str, manifest: dict[str, object] | None = None
    ) -> set[str]:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            plan_path = root / "plan.md"
            manifest_path = root / "manifest.json"
            plan_path.write_text(plan, encoding="utf-8")
            manifest_path.write_text(
                json.dumps(manifest or valid_manifest()), encoding="utf-8"
            )
            return {issue.code for issue in check_plan(plan_path, manifest_path)}

    def test_rejects_missing_read_predecessor(self) -> None:
        plan = valid_plan().replace(f"Step 1 [ ]: {READ_ACTION}", "Step 1 [ ]: Audit first.")
        self.assertIn("PLAN_READ_PREDECESSOR", self.run_check(plan))

    def test_rejects_broken_step_numbering(self) -> None:
        plan = valid_plan().replace("Step 2 [ ]:", "Step 3 [ ]:")
        self.assertIn("PLAN_STEP_SEQUENCE", self.run_check(plan))

    def test_rejects_reverse_stage_dependency(self) -> None:
        manifest = valid_manifest()
        stages = manifest["stages"]
        assert isinstance(stages, dict)
        stages["VAL"]["depends_on_stages"] = ["FUTURE"]
        stages["FUTURE"] = {
            "owner": "cosmolkit-model",
            "target_root": "crates/cosmolkit-model/src",
            "depends_on_stages": [],
        }
        codes = {issue.code for issue in validate_dependencies(manifest)}
        self.assertIn("MANIFEST_REVERSE_DEPENDENCY", codes)

    def test_rejects_missing_historical_file(self) -> None:
        manifest = copy.deepcopy(valid_manifest())
        coverage = manifest["historical_coverage"]
        assert isinstance(coverage, dict)
        coverage["expected_unique_paths"] = 3
        self.assertIn("PLAN_HISTORICAL_COVERAGE", self.run_check(valid_plan(), manifest))

    def test_rejects_zero_historical_coverage(self) -> None:
        plan = valid_plan().replace(
            "| `src/a.rs` | types / VAL |\n| `tests/a.rs` | types / VAL |\n", ""
        )
        self.assertIn("PLAN_HISTORICAL_COVERAGE", self.run_check(plan))


if __name__ == "__main__":
    unittest.main()
