"""Small build-policy regressions, without rebuilding WASM."""

import tomllib
import unittest

from run import workspace_manifest


class ReleaseBuildTests(unittest.TestCase):
    def test_isolated_workspace_has_explicit_release_optimization(self):
        metadata = {"version": "0.5.0-rc.17", "edition": "2024", "license": "MIT"}
        manifest = tomllib.loads(workspace_manifest(metadata))
        self.assertEqual(manifest["workspace"]["package"], metadata)
        self.assertEqual(manifest["profile"]["release"], {
            "opt-level": 3, "lto": "fat", "codegen-units": 1,
        })


if __name__ == "__main__":
    unittest.main()
