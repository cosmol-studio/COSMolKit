"""Build-mode regressions; also validate real extension and host executables."""

from pathlib import Path
import tomllib
import unittest


ROOT = Path(__file__).resolve().parents[2]


class PythonBuildModeTests(unittest.TestCase):
    def setUp(self):
        self.manifest = tomllib.loads((ROOT / "python/Cargo.toml").read_text())

    def enabled_features(self, profile):
        enabled = set()
        pending = [profile]
        while pending:
            feature = pending.pop()
            if feature not in enabled:
                enabled.add(feature)
                pending.extend(self.manifest["features"].get(feature, []))
        return enabled

    def test_extension_profiles_do_not_embed_python(self):
        for profile in ("extension", "release-abi3-py39", "release-cpython"):
            with self.subTest(profile=profile):
                features = self.enabled_features(profile)
                self.assertIn("pyo3/extension-module", features)
        self.assertIn("pyo3/abi3-py39", self.enabled_features("release-abi3-py39"))
        self.assertNotIn("pyo3/abi3-py39", self.enabled_features("release-cpython"))
        self.assertNotIn(
            "extension-module", self.manifest["dependencies"]["pyo3"]["features"]
        )

    def test_host_profiles_remain_separate_and_linkable(self):
        self.assertNotIn("pyo3/extension-module", self.enabled_features("default"))
        self.assertNotIn("abi3-py39", self.manifest["dependencies"]["pyo3"]["features"])
        self.assertIn("experimental-inspect", self.manifest["dependencies"]["pyo3"]["features"])

    def test_host_loader_uses_the_selected_python_configuration(self):
        build = (ROOT / "python/build.rs").read_text()
        self.assertIn("pyo3-build-config", self.manifest["build-dependencies"])
        self.assertIn("pyo3_build_config::add_libpython_rpath_link_args()", build)
        self.assertIn('CARGO_FEATURE_EXTENSION', build)
        for path in ("/home/", "/usr/lib/", "libpython3."):
            self.assertNotIn(path, build)

    def test_maturin_supports_command_local_extension_mode(self):
        distribution = tomllib.loads((ROOT / "python/pyproject.toml").read_text())
        environment = tomllib.loads((ROOT / "pyproject.toml").read_text())
        lock = tomllib.loads((ROOT / "uv.lock").read_text())
        self.assertIn("maturin==1.15.0", distribution["build-system"]["requires"])
        self.assertIn("maturin>=1.9.4,<2.0", environment["dependency-groups"]["dev"])
        maturin = next(package for package in lock["package"] if package["name"] == "maturin")
        version = tuple(int(part) for part in maturin["version"].split("."))
        self.assertGreaterEqual(version, (1, 9, 4))
        self.assertLess(version, (2, 0, 0))
        self.assertEqual(distribution["tool"]["maturin"]["module-name"], "cosmolkit._native")


if __name__ == "__main__":
    unittest.main()
