"""Contract identity, migration aliases, and binding pairing regressions."""

import copy
from pathlib import Path
import tempfile
import unittest

from route_contract import MANIFEST, PAGES, counterpart, load_contract, redirect_rules


class RouteContractTests(unittest.TestCase):
    def test_all_legacy_forms_redirect_directly_to_python(self):
        rules = redirect_rules(PAGES)
        for page in PAGES:
            for legacy in page.get("legacy_paths", []):
                for suffix in ("", "/", ".html"):
                    self.assertEqual(rules[legacy + suffix], (page["path"], "301"))
                    self.assertNotIn(page["path"], rules)
        self.assertEqual(rules["/api"], ("/python/api", "301"))

    def test_missing_counterpart_and_common_page_landings(self):
        api = next(p for p in PAGES if p["path"] == "/python/api")
        self.assertIsNone(counterpart(PAGES, api, "javascript"))
        home = next(p for p in PAGES if p["path"] == "/")
        self.assertEqual(counterpart(PAGES, home, "javascript")["path"], "/javascript")
        self.assertEqual(counterpart(PAGES, home, "python")["path"], "/python")

    def test_published_counterparts_pair_by_topic_not_slug(self):
        pages = copy.deepcopy(PAGES)
        api = next(p for p in pages if p["path"] == "/python/api")
        js = dict(api, component="JavaScriptApi", path="/javascript/reference", binding="javascript")
        pages.append(js)
        self.assertEqual(counterpart(pages, api, "javascript"), js)
        self.assertEqual(counterpart(pages, js, "python"), api)
        js["status"] = "placeholder"
        self.assertIsNone(counterpart(pages, api, "javascript"))

    def test_invalid_contracts_fail_before_generation(self):
        original = MANIFEST.read_text(encoding="utf-8")
        changes = (
            ('path = "/python/api"', 'path = "/python/molecule"'),
            ('component = "Api"', 'component = "Molecule"'),
            ('topic = "api"', 'topic = "molecule"'),
            ('docname = "api"', 'docname = "molecule"'),
            ('path = "/python/api"', 'path = "/api"'),
            ('legacy_paths = ["/api"]', 'legacy_paths = ["/python"]'),
            ('legacy_paths = ["/api"]', 'legacy_paths = ["/molecule"]'),
            ('path = "/python/api"', 'path = "/python/../api"'),
            ('status = "placeholder"\nindexable = false', 'status = "placeholder"\nindexable = true'),
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "routes.toml"
            for old, new in changes:
                with self.subTest(new=new):
                    self.assertIn(old, original)
                    path.write_text(original.replace(old, new), encoding="utf-8")
                    with self.assertRaises(ValueError):
                        load_contract(path)
