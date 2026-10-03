"""Release bump regressions; no repository files are written."""

import importlib.util
from pathlib import Path
import subprocess
import sys
import tomllib
import unittest


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("bump", ROOT / "dev/tools/bump_rust_version.py")
bump = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bump)


class ReleaseBumpTests(unittest.TestCase):
    def setUp(self):
        files = ["Cargo.toml", "python/pyproject.toml", bump.README]
        files += bump.PACKAGE_VERSIONS
        files += [file for file, _, _ in bump.DEPENDENCIES]
        self.original = {file: (ROOT / file).read_text() for file in files}

    def test_python_package_dependency_and_distribution_move_together(self):
        updated = bump.prepare_updates(self.original, "0.5.0-rc.9")
        manifest = tomllib.loads(updated["python/Cargo.toml"])
        self.assertEqual(manifest["package"]["version"], "0.5.0-rc.9")
        self.assertEqual(manifest["dependencies"]["cosmolkit"]["version"], "0.5.0-rc.9")
        self.assertEqual(tomllib.loads(updated["python/pyproject.toml"])["project"]["version"], "0.5.0rc9")

    def test_every_existing_internal_path_dependency_is_covered(self):
        updated = bump.prepare_updates(self.original, "0.5.1-rc.10")
        workspace = tomllib.loads(updated["Cargo.toml"])
        self.assertEqual(workspace["workspace"]["package"]["version"], "0.5.1-rc.10")
        checked = 0
        for member in workspace["workspace"]["members"]:
            file = member + "/Cargo.toml"
            manifest = tomllib.loads(updated.get(file, (ROOT / file).read_text()))
            for section in ("dependencies", "dev-dependencies", "build-dependencies"):
                for name, dependency in manifest.get(section, {}).items():
                    if name.startswith("cosmolkit") and isinstance(dependency, dict) and "path" in dependency and "version" in dependency:
                        self.assertEqual(dependency["version"], "0.5.1-rc.10", (file, name))
                        checked += 1
        self.assertGreater(checked, 40)

    def test_stable_version_and_idempotence(self):
        updated = bump.prepare_updates(self.original, "0.5.0")
        self.assertEqual(tomllib.loads(updated["python/pyproject.toml"])["project"]["version"], "0.5.0")
        self.assertEqual(bump.prepare_updates(updated, "0.5.0"), updated)

    def test_only_marked_installation_example_changes(self):
        updated = bump.prepare_updates(self.original, "0.5.1-rc.10")
        before, body = self.original[bump.README].split(bump.BEGIN)
        _, after = body.split(bump.END)
        new_before, new_body = updated[bump.README].split(bump.BEGIN)
        _, new_after = new_body.split(bump.END)
        self.assertEqual((before, after), (new_before, new_after))

    def test_missing_package_field_fails_before_any_write(self):
        original = self.original.copy()
        version = tomllib.loads(original["python/Cargo.toml"])["package"]["version"]
        original["python/Cargo.toml"] = original["python/Cargo.toml"].replace('version = "' + version + '"\n', "", 1)
        with self.assertRaisesRegex(ValueError, "package.version"):
            bump.prepare_updates(original, "0.5.0-rc.9")
        self.assertEqual(self.original, {file: (ROOT / file).read_text() for file in self.original})

    def test_field_edit_preserves_other_values_and_comments(self):
        text = '[package]\nversion = "0.1.0" # release\nname = "example"\n[dependencies]\nx = "0.1.0"\n'
        self.assertEqual(bump.update_field(text, "package", "version", "0.5.0-rc.9"), text.replace('version = "0.1.0"', 'version = "0.5.0-rc.9"'))

    def test_lock_refresh_does_not_update_all_third_party_dependencies(self):
        self.assertEqual(bump.LOCK_REFRESH, ["cargo", "update", "--workspace"])

    def test_cli_dry_run_and_invalid_version_never_write(self):
        path = ROOT / "dev/tools/bump_rust_version.py"
        result = subprocess.run([sys.executable, str(path), "0.5.1-rc.10", "--dry-run"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--- python/Cargo.toml", result.stdout)
        self.assertIn("Would then run cargo update --workspace.", result.stdout)
        for version in ("0.5", "0.5.0-rc.", "0.5.0-rc.09", "0.5.0-beta.1"):
            with self.subTest(version=version):
                result = subprocess.run([sys.executable, str(path), version], capture_output=True, text=True)
                self.assertEqual(result.returncode, 2)
        self.assertEqual(self.original, {file: (ROOT / file).read_text() for file in self.original})


if __name__ == "__main__":
    unittest.main()
