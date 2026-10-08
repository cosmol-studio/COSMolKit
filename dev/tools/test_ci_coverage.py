"""Check coverage-stage shell syntax and the shared instrumentation boundary."""

from pathlib import Path
import subprocess
import tomllib
import unittest


ROOT = Path(__file__).resolve().parents[2]


def shell_steps():
    # This workflow's literal run blocks have fixed YAML indentation. Do not
    # execute expressions or provision dependencies to inspect shell syntax.
    steps = {}
    name = None
    script = None
    for line in (ROOT / ".github/workflows/coverage.yml").read_text().splitlines():
        if script is not None and line and not line.startswith("                  "):
            steps[name] = "\n".join(script) + "\n"
            script = None
        if line.startswith("            - name: "):
            name = line.removeprefix("            - name: ")
        elif line == "              run: |":
            script = []
        elif script is not None:
            script.append(line[18:] if line else "")
    if script is not None:
        steps[name] = "\n".join(script) + "\n"
    return steps


class CoverageWorkflowTests(unittest.TestCase):
    def test_core_has_no_runtime_features_or_forwarders(self):
        core = tomllib.loads((ROOT / "crates/cosmolkit-core/Cargo.toml").read_text())
        self.assertEqual(core.get("features", {}), {})
        facade = tomllib.loads((ROOT / "crates/cosmolkit/Cargo.toml").read_text())
        self.assertEqual(
            facade["features"]["op-contracts-strict"],
            ["runtime-invariants", "op-contracts"],
        )
        for feature_values in facade["features"].values():
            self.assertFalse(any("cosmolkit-core" in value and "/op-contracts" in value for value in feature_values))

    def test_feature_matrix_uses_declared_leaf_capabilities_and_user_bundles(self):
        workflow = (ROOT / ".github/workflows/features.yml").read_text()
        self.assertIn("cargo test -p cosmolkit --no-default-features --profile dev-test --features op-contracts-strict", workflow)
        self.assertIn("cargo test -p cosmolkit --no-default-features --profile dev-test --features core,op-contracts-strict --test feature_selection", workflow)
        loop = workflow.split("for feature in \\\n", 1)[1].split("                  do", 1)[0]
        selected = loop.replace("\\", "").split()
        features = tomllib.loads((ROOT / "crates/cosmolkit/Cargo.toml").read_text())["features"]
        self.assertEqual(len(selected), len(set(selected)))
        self.assertTrue(set(selected) <= features.keys(), set(selected) - features.keys())
        self.assertTrue({"core", "bio", "conformer", "forcefields", "search", "inchi", "fingerprints", "depict", "batch"} <= set(selected))
        self.assertTrue({"cap-io", "cap-serialization", "cap-descriptors", "cap-stereoisomers", "cap-confseq", "cap-hashing"} <= set(selected))

    def test_instrumented_build_and_profile_directories_agree_before_show_env(self):
        text = (ROOT / ".github/workflows/coverage.yml").read_text()
        env = text.split("        env:\n", 1)[1].split("        steps:\n", 1)[0]
        values = {}
        for line in env.splitlines():
            if line.strip() and not line.strip().startswith("#"):
                key, value = line.strip().split(":", 1)
                values[key] = value.strip()
        self.assertEqual(values["CARGO_TARGET_DIR"], values["CARGO_LLVM_COV_TARGET_DIR"])
        self.assertEqual(values["CARGO_TARGET_DIR"], "target/coverage-build")

    def test_every_multiline_shell_step_parses(self):
        steps = shell_steps()
        self.assertGreaterEqual(len(steps), 7)
        for name, script in steps.items():
            with self.subTest(step=name):
                result = subprocess.run(["bash", "-n"], input=script, text=True, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stderr)

    def test_build_regressions_preparation_and_corpus_share_instrumentation(self):
        steps = shell_steps()
        stages = {
            "Build libraries before fetching third-party test sources": "cargo build",
            "Run all default crate regression suites with coverage": "cargo test",
            "Prepare and validate all reference values": "cargo build",
            "Run all parity integration targets with coverage": "cargo test",
        }
        for name, command in stages.items():
            with self.subTest(step=name):
                script = steps[name]
                self.assertEqual(script.count("cargo llvm-cov show-env --sh"), 1)
                self.assertIn('eval "$coverage_env"', script)
                self.assertLess(script.index('export CARGO_TARGET_DIR="$CARGO_LLVM_COV_TARGET_DIR"'), script.index(command))
                self.assertIn("--profile dev-test", script)
                self.assertNotIn("--release", script)
                self.assertIn("cosmolkit/op-contracts-strict", script)
                self.assertNotIn("--no-clean", script)
                self.assertNotIn("--no-report", script)
        self.assertNotIn("--test reference_parity", steps["Run all default crate regression suites with coverage"])
        parity = steps["Run all parity integration targets with coverage"]
        self.assertIn("cargo test -p cosmolkit-parity-tests", parity)
        self.assertIn("--test '*' --no-fail-fast", parity)
        prepare = steps["Prepare and validate all reference values"]
        self.assertIn('$CARGO_TARGET_DIR/dev-test/cosmolkit-parity-tests', prepare)
        self.assertNotIn('$CARGO_TARGET_DIR/release/', prepare)

    def test_report_is_separate_and_test_failures_remain_failures(self):
        steps = shell_steps()
        self.assertIn("cargo llvm-cov report", steps["Generate coverage reports"])
        self.assertIn("--profile dev-test", steps["Generate coverage reports"])
        self.assertIn('exit "$status"', steps["Run all default crate regression suites with coverage"])
        self.assertIn("set -o pipefail", steps["Run all parity integration targets with coverage"])
        self.assertIn("run: exit 1", (ROOT / ".github/workflows/coverage.yml").read_text())

    def test_development_tests_and_distribution_have_separate_profiles(self):
        profiles = tomllib.loads((ROOT / "Cargo.toml").read_text())["profile"]
        self.assertEqual(profiles["release"]["opt-level"], 3)
        self.assertEqual(profiles["release"]["lto"], "fat")
        self.assertEqual(profiles["release"]["codegen-units"], 1)
        self.assertEqual(profiles["dev-test"]["inherits"], "release")
        self.assertEqual(profiles["dev-test"]["lto"], "off")
        self.assertEqual(profiles["dev-test"]["codegen-units"], 16)
        workflow = (ROOT / ".github/workflows/publish.yml").read_text()
        self.assertEqual(workflow.count("--release"), 5)
        self.assertNotIn("--profile dist", workflow)
        self.assertNotIn("--profile dev-test", workflow)

    def test_each_wheel_build_installs_the_inherited_sccache_wrapper(self):
        workflow = (ROOT / ".github/workflows/publish.yml").read_text()
        wheel_steps = workflow.split("    build-wheels:\n", 1)[1].split("    build-sdist:\n", 1)[0]
        builds = wheel_steps.split("uses: PyO3/maturin-action@v1")[1:]
        self.assertEqual(len(builds), 5)
        for build in builds:
            with self.subTest(step=build.split("args:", 1)[0]):
                self.assertIn("sccache: true", build.split("args:", 1)[0])


if __name__ == "__main__":
    unittest.main()
