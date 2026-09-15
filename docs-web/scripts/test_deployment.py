"""Offline deployment regression tests; no credentials or network required."""

from pathlib import Path
import shutil
import tempfile
import unittest
from xml.etree import ElementTree

from check_ssg_output import check_output
from flatten_html_routes import ROUTES, flatten_html_routes
from generate_sitemap import BASE_URL, SITEMAP_NAMESPACE, write_sitemap
from strip_client_runtime import strip_client_runtime


class DeploymentTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.public = Path(self.temporary.name)

    def page(self, relative, route, runtime=False):
        path = self.public / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        html = f'<link rel="canonical" href="{BASE_URL}{route}">'
        html += '<link rel="stylesheet" href="/style.css">' * 3
        html += '<main>Documentation</main>'
        if runtime:
            html += '<script>window.hydrate_queue=[];</script>'
            html += '<script>window.initial_dioxus_hydration_data=[];</script>'
            html += '<link rel="preload" as="script" href="/assets/cosmolkit-docs-web.js">'
            html += '<script type="module" src="/assets/cosmolkit-docs-web.js"></script>'
        path.write_text(html, encoding="utf-8")
        return path

    def prepare(self):
        self.page("index.html", "", runtime=True)
        for route in ROUTES:
            self.page(f"{route}/index.html", route, runtime=True)
        self.assertEqual(flatten_html_routes(self.public), 17)
        self.assertEqual(strip_client_runtime(self.public), 18)
        shutil.copyfile(
            Path(__file__).resolve().parents[1] / "deployment/public/_redirects",
            self.public / "_redirects",
        )
        self.assertEqual(write_sitemap(self.public), 15)

    def test_full_preparation_and_idempotence(self):
        self.prepare()
        self.assertEqual(check_output(self.public), 18)
        self.assertEqual(flatten_html_routes(self.public), 0)
        self.assertEqual(strip_client_runtime(self.public), 0)
        for route in ROUTES:
            self.assertFalse((self.public / route).exists())
            self.assertTrue((self.public / f"{route}.html").is_file())

    def test_flatten_preserves_unrelated_directories(self):
        module = self.page("_modules/index.html", "_modules/")
        asset = self.page("assets/index.html", "assets/")
        self.assertEqual(flatten_html_routes(self.public), 0)
        self.assertTrue(module.is_file())
        self.assertTrue(asset.is_file())

    def test_conflict_preflight_preserves_all_inputs(self):
        first = self.page("installation/index.html", "installation")
        second = self.page("validation/index.html", "validation")
        destination = self.page("validation.html", "validation")
        with self.assertRaises(FileExistsError):
            flatten_html_routes(self.public)
        self.assertTrue(first.is_file())
        self.assertTrue(second.is_file())
        self.assertTrue(destination.is_file())

    def test_unexpected_route_contents_are_not_removed(self):
        index = self.page("installation/index.html", "installation")
        extra = index.parent / "extra.txt"
        extra.write_text("keep", encoding="utf-8")
        with self.assertRaises(ValueError):
            flatten_html_routes(self.public)
        self.assertTrue(index.is_file())
        self.assertEqual(extra.read_text(encoding="utf-8"), "keep")

    def test_sitemap_excludes_flat_and_directory_utility_pages(self):
        self.page("index.html", "")
        self.page("installation/index.html", "installation")
        for route in ("search", "genindex", "py-modindex"):
            self.page(f"{route}.html", route)
            self.page(f"{route}/index.html", route)
        self.page("_modules/index.html", "_modules/")
        self.assertEqual(write_sitemap(self.public), 2)
        root = ElementTree.parse(self.public / "sitemap.xml")
        urls = [node.text for node in root.findall(f".//{{{SITEMAP_NAMESPACE}}}loc")]
        self.assertEqual(urls, [BASE_URL, BASE_URL + "installation"])

    def test_output_rejects_legacy_canonical(self):
        self.prepare()
        self.page("installation.html", "installation.html")
        with self.assertRaisesRegex(ValueError, "clean canonical"):
            check_output(self.public)

    def test_output_rejects_directory_slash_loop(self):
        self.prepare()
        self.page("installation/index.html", "installation")
        with self.assertRaisesRegex(ValueError, "trailing slash"):
            check_output(self.public)

    def test_output_requires_all_routes(self):
        self.prepare()
        (self.public / "javascript.html").unlink()
        with self.assertRaises(FileNotFoundError):
            check_output(self.public)

    def test_output_rejects_runtime(self):
        self.prepare()
        self.page("search.html", "search", runtime=True)
        with self.assertRaisesRegex(ValueError, "runtime"):
            check_output(self.public)

    def test_output_rejects_invalid_shell(self):
        self.prepare()
        page = self.public / "python.html"
        page.write_text(page.read_text(encoding="utf-8") + "<main></main>", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "SSG shell"):
            check_output(self.public)

    def test_output_rejects_nonpermanent_or_query_fragment_overrides(self):
        self.prepare()
        redirects = self.public / "_redirects"
        original = redirects.read_text(encoding="utf-8")
        for replacement in ("/search 302", "/search?fixed=1 301", "/search#fixed 301"):
            with self.subTest(replacement=replacement):
                redirects.write_text(original.replace("/search 301", replacement), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, "redirects must"):
                    check_output(self.public)

    def test_output_rejects_wrong_sitemap(self):
        self.prepare()
        sitemap = self.public / "sitemap.xml"
        sitemap.write_text(
            sitemap.read_text(encoding="utf-8").replace("/installation", "/search"),
            encoding="utf-8",
        )
        with self.assertRaisesRegex(ValueError, "sitemap"):
            check_output(self.public)


if __name__ == "__main__":
    unittest.main()
