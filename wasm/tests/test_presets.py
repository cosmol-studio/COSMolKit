"""Package selection and release labels, independent of chemistry APIs."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import features
from presets import PRESETS, npm_release


class PresetTests(unittest.TestCase):
    def test_facade_feature_forwarding(self):
        features.sync()

    def test_versions_and_default_tag(self):
        for release in ("0.5.1", "0.5.1-rc.5", "0.6.0-dev.1", "0.6.0-dev.4"):
            self.assertEqual(len({npm_release(release, p)[0] for p in PRESETS}), len(PRESETS))
            for preset in PRESETS:
                self.assertEqual(npm_release(release, preset), (
                    "@cosmol-studio/cosmolkit" + ("" if preset == "full" else f"-{preset}"),
                    release, "dev" if "-dev." in release else "rc" if "-rc." in release else "latest"))
        for version, preset in (("invalid", "full"), ("0.5.1", "unknown")):
            with self.assertRaises(ValueError):
                npm_release(version, preset)
