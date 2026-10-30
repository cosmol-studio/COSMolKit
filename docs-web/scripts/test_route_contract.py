"""Contract identity, migration aliases, and binding pairing regressions."""

import copy
from pathlib import Path
import tempfile
import unittest

from route_contract import MANIFEST, PAGES, counterpart, load_contract, redirect_rules


class RouteContractTests(unittest.TestCase):
    def test_python_navigation_prioritizes_guides_and_places_coming_soon_last(self):
        navigation = sorted(
            (p for p in PAGES if p["binding"] == "python" and "order" in p
             and p["status"] == "ready"),
            key=lambda p: p["order"],
        )
        names = [p["docname"].removeprefix("python/") for p in navigation]
        self.assertEqual(names[:4], ["quickstart", "molecule", "io", "batch"])
        self.assertNotIn("installation", names)
        self.assertNotIn("py-modindex", names)
        self.assertEqual(names[-1], "confseq")
        self.assertIn("Coming soon", navigation[-1]["label"])
        api = names.index("api")
        for guide in ("batch", "fingerprints", "descriptors", "mcs", "reaction",
                      "forcefields", "protein"):
            self.assertLess(names.index(guide), api)
        self.assertEqual(len({p["order"] for p in navigation}), len(navigation))

        index = MANIFEST.parent.parent / "docs/source/python/index.rst"
        guide_tree = index.read_text(encoding="utf-8").split(".. toctree::")[1]
        sphinx_order = [line.strip() for line in guide_tree.splitlines()
                        if line.startswith("   ") and not line.strip().startswith(":")]
        sidebar_guides = [p["docname"].removeprefix("python/") for p in navigation if p.get("summary")]
        self.assertEqual(sphinx_order, sidebar_guides)

    def test_installation_is_first_quickstart_section(self):
        source = MANIFEST.parent.parent / "docs/source/python/quickstart.rst"
        text = source.read_text(encoding="utf-8")
        self.assertLess(text.index("Installation\n"), text.index("Create and inspect a molecule\n"))
        self.assertIn("pip install cosmolkit", text)
        rules = redirect_rules(PAGES)
        self.assertEqual(rules["/python/installation"], ("/python/quickstart", "301"))

    def test_generated_api_definitions_only_live_in_api_pages(self):
        source = MANIFEST.parent.parent / "docs/source"
        for page in source.rglob("*.rst"):
            if page.relative_to(source).as_posix() in {"python/api.rst", "javascript/api.rst"}:
                continue
            with self.subTest(page=page.name):
                self.assertNotRegex(page.read_text(encoding="utf-8"),
                                    r"(?m)^\s*\.\.\s+(?:auto(?:module|class|function|method|attribute|data)|js:auto\w+)::")

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
        self.assertEqual(counterpart(PAGES, api, "javascript")["path"], "/javascript/api")
        guide = next(p for p in PAGES if p["path"] == "/python/molecule")
        self.assertIsNone(counterpart(PAGES, guide, "javascript"))
        home = next(p for p in PAGES if p["path"] == "/")
        self.assertEqual(counterpart(PAGES, home, "javascript")["path"], "/javascript")
        self.assertEqual(counterpart(PAGES, home, "python")["path"], "/python")

    def test_published_counterparts_pair_by_topic_not_slug(self):
        pages = copy.deepcopy(PAGES)
        pages = [p for p in pages if p["path"] != "/javascript/api"]
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
            ('docname = "python/api"', 'docname = "python/molecule"'),
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

    def test_language_sources_are_separate_with_unchanged_public_destinations(self):
        source = MANIFEST.parent.parent / "docs/source"
        self.assertFalse((MANIFEST.parent.parent / "python/docs/source").exists())
        self.assertEqual(next(p for p in PAGES if p["path"] == "/")["docname"], "index")
        for page in PAGES:
            docname = page.get("docname")
            if not docname or docname in {"search", "genindex", "index"}:
                continue
            with self.subTest(page=page["path"]):
                self.assertTrue(docname.startswith(page["binding"] + "/"))
                self.assertTrue((source / (docname + ".rst")).is_file())
                public = "/" + docname.removesuffix("/index")
                self.assertEqual(page["path"], public)
        for language in ("python", "javascript"):
            index = (source / language / "index.rst").read_text(encoding="utf-8")
            self.assertIn(".. toctree::", index)
            self.assertIn("   api", index)
