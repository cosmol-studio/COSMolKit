"""Version catalog publication and the real static browser module regressions."""

import copy
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from version_catalog import CATALOG_PATH, check_version_assets, load_catalog, validate_catalog, write_version_assets


class VersionTests(unittest.TestCase):
    def test_manual_catalog_and_exact_publication(self):
        catalog = load_catalog()
        self.assertEqual(catalog["versions"][0], {"version": "latest", "url": "https://kit.cosmol.org/"})
        self.assertIn({"version": "0.3.0", "url": "https://c6862989.cosmolkit-docs-web.pages.dev/"}, catalog["versions"])
        with tempfile.TemporaryDirectory() as temporary:
            public = Path(temporary)
            write_version_assets(public)
            self.assertEqual((public / "versions.json").read_bytes(), CATALOG_PATH.read_bytes())
            check_version_assets(public)
            headers = public / "_headers"
            headers.write_text(headers.read_text().replace("Access-Control-Allow-Origin: *", ""))
            with self.assertRaisesRegex(ValueError, "CORS"):
                check_version_assets(public)

    def test_rejects_bad_catalog(self):
        catalog = load_catalog()
        invalid = [None, {}, {"schema_version": 1, "versions": []}]
        for url in ("javascript:alert(1)", "http://example.com/", "https://user@example.com/", "https://example.com/topic", "https://example.com/?x=1", "//example.com/"):
            modified = copy.deepcopy(catalog)
            modified["versions"].append({"version": "unsafe", "url": url})
            invalid.append(modified)
        duplicate = copy.deepcopy(catalog)
        duplicate["versions"].append(duplicate["versions"][1])
        invalid.append(duplicate)
        for value in invalid:
            with self.subTest(value=value), self.assertRaises(ValueError):
                validate_catalog(value)

    @unittest.skipUnless(shutil.which("node"), "Node required for DOM boundary tests")
    def test_actual_browser_module(self):
        subprocess.run(["node", str(Path(__file__).with_name("test_version_switch.cjs"))], check=True)


if __name__ == "__main__":
    unittest.main()
