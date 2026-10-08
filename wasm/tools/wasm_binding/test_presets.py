import json
from pathlib import Path
import tempfile
import tomllib
import unittest
from presets import PRESETS, MODULE_GROUPS, active_features, npm_release, selected_names
from run import CONFIG, TEST_DIR, prepare_package, selected_tests


class PresetTests(unittest.TestCase):
    def test_every_module_and_suite_is_assigned(self):
        config = tomllib.loads(CONFIG.read_text())["crates"][0]["wasm"]
        self.assertEqual(set(config["custom_rust_modules"]), {n for names in MODULE_GROUPS.values() for n in names.split()})
        for suffix in (".mjs", ".ts"):
            self.assertEqual(set(selected_tests(suffix, active_features("full"))), set(TEST_DIR.glob("*" + suffix)))

    def test_prerequisites_and_absent_domains(self):
        for preset in PRESETS:
            active = active_features(preset)
            self.assertTrue({"core", "search", "smiles", "depict"} <= active)
            if preset != "full":
                self.assertFalse({"full", "batch", "serialization"} & active)
        self.assertTrue({"forcefields", "alignment"} <= active_features("core-3d"))
        self.assertIn("inchi", active_features("core-inchi"))
        self.assertNotIn("inchi", active_features("core"))
        self.assertNotIn("bio_readers", selected_names(MODULE_GROUPS, active_features("core")))

    def test_versions_and_default_tag(self):
        for release in ("0.5.0", "0.5.0-rc.15"):
            self.assertEqual(len({npm_release(release, p)[0] for p in PRESETS}), len(PRESETS))
        self.assertEqual(npm_release("0.5.0", "full"), ("0.5.0", "latest"))
        self.assertEqual(npm_release("0.5.0-rc.15", "full"), ("0.5.0-rc.15", "rc"))
        self.assertEqual(npm_release("0.5.0", "core-bio"), ("0.5.0-core-bio.0", "core-bio"))

    def test_package_records_exact_preset(self):
        with tempfile.TemporaryDirectory() as directory:
            package = Path(directory)
            prepare_package(package, "wasm_wasm", {
                "version": "0.5.0-rc.15", "license": "MIT",
                "repository": "https://github.com/cosmol-studio/COSMolKit",
            }, "core-bio")
            metadata = json.loads((package / "package.json").read_text())
            self.assertEqual(metadata["version"], "0.5.0-rc.15.core-bio.0")
            self.assertEqual(metadata["cosmolkitVersion"], "0.5.0-rc.15")
            self.assertEqual(metadata["cosmolkitPreset"], "core-bio")
            self.assertEqual(metadata["publishConfig"]["tag"], "core-bio")
