from __future__ import annotations

import contextlib
import io
import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import run_unit_validation as runner


def selected(
    *,
    package: str = "cosmolkit-types",
    features: tuple[str, ...] = (),
    runtime_strict: bool = False,
) -> runner.SelectedValidation:
    return runner.SelectedValidation(
        unit_id="VAL-vocab",
        layer="public" if runtime_strict else "detached",
        package=package,
        target="migration_val_vocab",
        features=features,
        profile="release",
        require_nonzero_tests=True,
        require_runtime_strict=runtime_strict,
    )


def metadata_result(package: str, features: tuple[str, ...] = ()) -> subprocess.CompletedProcess[str]:
    metadata = {
        "packages": [
            {
                "name": package,
                "features": {feature: [] for feature in features},
            }
        ]
    }
    return subprocess.CompletedProcess(
        args=runner.cargo_metadata_command(),
        returncode=0,
        stdout=json.dumps(metadata),
    )


class RunUnitValidationRegressionTests(unittest.TestCase):
    def execute_with_results(
        self,
        validation: runner.SelectedValidation,
        *results: subprocess.CompletedProcess[str],
    ) -> tuple[int, dict[str, object]]:
        with tempfile.TemporaryDirectory() as directory:
            evidence_path = Path(directory) / "evidence.json"
            with patch.object(runner, "run_process", side_effect=results):
                with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(
                    io.StringIO()
                ):
                    exit_code = runner.execute_validation(
                        validation,
                        Path(directory),
                        evidence_path,
                    )
            evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
        return exit_code, evidence

    def test_rejects_successful_target_that_executes_zero_tests(self) -> None:
        validation = selected()
        cargo_test = subprocess.CompletedProcess(
            args=runner.cargo_test_command(validation),
            returncode=0,
            stdout=(
                "test result: ok. 0 passed; 0 failed; 0 ignored; "
                "0 measured; 0 filtered out; finished in 0.00s\n"
            ),
        )
        exit_code, evidence = self.execute_with_results(
            validation, metadata_result(validation.package), cargo_test
        )
        self.assertEqual(exit_code, 2)
        self.assertEqual(evidence["status"], "zero_tests_rejected")
        self.assertEqual(evidence["test_count"], 0)

    def test_rejects_public_validation_without_runtime_strict_feature(self) -> None:
        validation = selected(package="cosmolkit", runtime_strict=True)
        with self.assertRaisesRegex(runner.ConfigurationError, "does not select"):
            runner.validate_features(validation, {"op-contracts-strict"})

    def test_rejects_unknown_manifest_feature(self) -> None:
        validation = selected(features=("not-a-real-feature",))
        with self.assertRaisesRegex(runner.ConfigurationError, "unknown"):
            runner.validate_features(validation, set())

    def test_dependency_selectors_are_scoped_to_declared_dependencies(self) -> None:
        metadata = {
            "packages": [
                {
                    "name": "cosmolkit",
                    "features": {"valence": []},
                    "dependencies": [
                        {"name": "cosmolkit-core", "rename": None},
                        {"name": "renamed-package", "rename": "renamed"},
                    ],
                },
                {
                    "name": "cosmolkit-core",
                    "features": {"op-contracts-strict": []},
                    "dependencies": [],
                },
                {
                    "name": "renamed-package",
                    "features": {"strict": []},
                    "dependencies": [],
                },
                {
                    "name": "unrelated",
                    "features": {"must-not-leak": []},
                    "dependencies": [],
                },
            ]
        }
        self.assertEqual(
            runner.dependency_feature_selectors(metadata, "cosmolkit"),
            {"cosmolkit-core/op-contracts-strict", "renamed/strict"},
        )

    def test_propagates_failing_test_subprocess_and_records_evidence(self) -> None:
        validation = selected()
        cargo_test = subprocess.CompletedProcess(
            args=runner.cargo_test_command(validation),
            returncode=7,
            stdout="compilation or test failure\n",
        )
        exit_code, evidence = self.execute_with_results(
            validation, metadata_result(validation.package), cargo_test
        )
        self.assertEqual(exit_code, 7)
        self.assertEqual(evidence["exit_code"], 7)
        self.assertEqual(evidence["status"], "test_command_failed")
        self.assertIn("compilation or test failure", evidence["output"])


if __name__ == "__main__":
    unittest.main()
