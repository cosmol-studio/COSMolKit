"""Offline deployment regression tests; no credentials or network required."""

from pathlib import Path
import json
import shutil
import struct
import tempfile
import unittest
from unittest.mock import MagicMock, patch
from urllib.error import HTTPError
from email.message import Message
from xml.etree import ElementTree

from check_ssg_output import INDEXNOW_KEY, PROJECT_LINKS, SEARCH_BUNDLE_PREFIX, WEBSITE_JSON_LD, SOCIAL_IMAGE_URL, check_output
from flatten_html_routes import ROUTES, flatten_html_routes
from generate_sitemap import BASE_URL, EXCLUDED_ROUTES, SITEMAP_NAMESPACE, write_sitemap
from strip_client_runtime import strip_client_runtime
from prepare_deployment import SOCIAL_IMAGE_SOURCE, download_social_image, write_route_assets


SEARCH_SCRIPTS = (SEARCH_BUNDLE_PREFIX + "-fixture.js", SEARCH_BUNDLE_PREFIX + "_bg-fixture.wasm")
SEARCH_LOADER = f'<script id="docs-search-loader" type="module" data-search-bindings="{SEARCH_SCRIPTS[0]}" data-search-wasm="{SEARCH_SCRIPTS[1]}"></script>'
PNG_HEADER = b"\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR" + struct.pack(">II", 1200, 630)


class DeploymentTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.public = Path(self.temporary.name)

    def page(self, relative, route, runtime=False):
        path = self.public / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        title = f"COSMolKit {route or 'Documentation'}"
        description = f"Read the COSMolKit {route or 'documentation'} guide and reference."
        robots = "noindex, follow" if route in EXCLUDED_ROUTES else "index, follow"
        meta = {
            "description": description, "robots": robots,
            "og:type": "website", "og:site_name": "COSMolKit Documentation",
            "og:title": title, "og:description": description, "og:url": BASE_URL + route,
            "og:image": SOCIAL_IMAGE_URL, "og:image:type": "image/png",
            "og:image:width": "1200", "og:image:height": "630", "og:image:alt": "Documentation",
            "twitter:card": "summary_large_image", "twitter:title": title,
            "twitter:description": description, "twitter:image": SOCIAL_IMAGE_URL,
            "twitter:image:alt": "Documentation",
        }
        html = f'<!doctype html><html><head><title>{title}</title><link rel="canonical" href="{BASE_URL}{route}">'
        for key, value in meta.items():
            attribute = "property" if key.startswith("og:") else "name"
            html += f'<meta {attribute}="{key}" content="{value}">'
        html += '<link rel="stylesheet" href="/style.css">' * 3
        if route == "":
            html += f'<script type="application/ld+json">{json.dumps(WEBSITE_JSON_LD)}</script>'
        html += '</head><body><main role="main"><h1>Documentation</h1></main>'
        html += '<nav class="cosmolkit-project-links">'
        html += ''.join(f'<a href="{link}">Project</a>' for link in sorted(PROJECT_LINKS))
        html += '</nav>'
        if route == "python/search":
            html += '<form id="docs-search-form" action="/python/search" method="get"><input name="q" type="search"></form><p id="docs-search-status"></p><div id="search-results"></div>'
            html += SEARCH_LOADER
        if runtime:
            html += '<script>window.hydrate_queue=[];</script>'
            html += '<script>window.initial_dioxus_hydration_data=[];</script>'
            html += '<link rel="preload" as="script" href="/assets/cosmolkit-docs-web.js">'
            html += '<script type="module" src="/assets/cosmolkit-docs-web.js"></script>'
        html += '</body></html>'
        path.write_text(html, encoding="utf-8")
        return path

    def prepare(self):
        self.page("index.html", "", runtime=True)
        for route in ROUTES:
            self.page(f"{route}/index.html", route, runtime=True)
        self.assertEqual(flatten_html_routes(self.public), 17)
        self.assertEqual(strip_client_runtime(self.public), 18)
        docs_web = Path(__file__).resolve().parents[1]
        (self.public / "social-card.png").write_bytes(PNG_HEADER)
        for name in ("404.html",):
            shutil.copyfile(docs_web / "deployment/public" / name, self.public / name)
        shutil.copyfile(
            docs_web.parent / "python/docs/source/robots.txt", self.public / "robots.txt"
        )
        shutil.copyfile(
            docs_web.parent / f"python/docs/source/{INDEXNOW_KEY}.txt", self.public / f"{INDEXNOW_KEY}.txt"
        )
        write_route_assets(self.public)
        for src in SEARCH_SCRIPTS:
            path = self.public / src.lstrip("/")
            if not path.exists():
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b"\x00asm\x01\x00\x00\x00" if path.suffix == ".wasm" else b"// Generated binding fixture")
        self.assertEqual(write_sitemap(self.public), 13)

    def test_full_preparation_and_idempotence(self):
        self.prepare()
        self.assertEqual(check_output(self.public), 18)
        self.assertEqual(flatten_html_routes(self.public), 0)
        self.assertEqual(strip_client_runtime(self.public), 0)
        for route in ROUTES:
            self.assertFalse((self.public / route / "index.html").exists())
            self.assertTrue((self.public / f"{route}.html").is_file())

    def test_flatten_preserves_unrelated_directories(self):
        module = self.page("python/_modules/index.html", "python/_modules/")
        asset = self.page("assets/index.html", "assets/")
        self.assertEqual(flatten_html_routes(self.public), 0)
        self.assertTrue(module.is_file())
        self.assertTrue(asset.is_file())

    def test_conflict_preflight_preserves_all_inputs(self):
        first = self.page("python/installation/index.html", "python/installation")
        second = self.page("validation/index.html", "validation")
        destination = self.page("validation.html", "validation")
        with self.assertRaises(FileExistsError):
            flatten_html_routes(self.public)
        self.assertTrue(first.is_file())
        self.assertTrue(second.is_file())
        self.assertTrue(destination.is_file())

    def test_unexpected_route_contents_are_not_removed(self):
        index = self.page("python/installation/index.html", "python/installation")
        extra = index.parent / "extra.txt"
        extra.write_text("keep", encoding="utf-8")
        with self.assertRaises(ValueError):
            flatten_html_routes(self.public)
        self.assertTrue(index.is_file())
        self.assertEqual(extra.read_text(encoding="utf-8"), "keep")

    def test_sitemap_excludes_flat_and_directory_utility_pages(self):
        self.page("index.html", "")
        self.page("python/installation/index.html", "python/installation")
        for route in ("python/search", "python/genindex", "python/py-modindex", "javascript", "benchmarks"):
            self.page(f"{route}.html", route)
            self.page(f"{route}/index.html", route)
        self.page("python/_modules/index.html", "python/_modules/")
        self.assertEqual(write_sitemap(self.public), 2)
        root = ElementTree.parse(self.public / "sitemap.xml")
        urls = [node.text for node in root.findall(f".//{{{SITEMAP_NAMESPACE}}}loc")]
        self.assertEqual(urls, [BASE_URL, BASE_URL + "python/installation"])

    def test_sitemap_excludes_not_found_page(self):
        self.page("index.html", "")
        # Error pages do not need a canonical link and must never enter the sitemap.
        (self.public / "404.html").write_text("<h1>Page not found</h1>", encoding="utf-8")
        self.assertEqual(write_sitemap(self.public), 1)
        root = ElementTree.parse(self.public / "sitemap.xml")
        self.assertEqual(
            [node.text for node in root.findall(f".//{{{SITEMAP_NAMESPACE}}}loc")],
            [BASE_URL],
        )

    def test_output_requires_crawler_assets(self):
        self.prepare()
        for name in ("404.html", "robots.txt", f"{INDEXNOW_KEY}.txt"):
            with self.subTest(name=name):
                path = self.public / name
                original = path.read_bytes()
                path.unlink()
                with self.assertRaises(FileNotFoundError):
                    check_output(self.public)
                path.write_bytes(original)

    def test_output_rejects_indexable_or_dead_end_not_found_page(self):
        self.prepare()
        page = self.public / "404.html"
        original = page.read_text(encoding="utf-8")
        for old, new in (("noindex", "index"), ('href="/"', 'href="/missing"')):
            with self.subTest(old=old):
                page.write_text(original.replace(old, new), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, "404.html"):
                    check_output(self.public)

    def test_output_rejects_wrong_robots_sitemap_or_blocked_crawlers(self):
        self.prepare()
        robots = self.public / "robots.txt"
        original = robots.read_text(encoding="utf-8")
        for old, new in (
            (BASE_URL + "sitemap.xml", "https://tools.cosmol.org/sitemap.xml"),
            ("Allow: /", "Disallow: /"),
            ("User-agent: *", "User-agent: unknownbot"),
        ):
            with self.subTest(new=new):
                robots.write_text(original.replace(old, new), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, "robots.txt"):
                    check_output(self.public)

    def test_output_rejects_legacy_canonical(self):
        self.prepare()
        self.page("python/installation.html", "python/installation.html")
        with self.assertRaisesRegex(ValueError, "clean canonical"):
            check_output(self.public)

    def test_output_rejects_directory_slash_loop(self):
        self.prepare()
        self.page("python/installation/index.html", "python/installation")
        with self.assertRaisesRegex(ValueError, "trailing slash"):
            check_output(self.public)

    def test_output_requires_all_routes(self):
        self.prepare()
        (self.public / "javascript.html").unlink()
        with self.assertRaises(FileNotFoundError):
            check_output(self.public)

    def test_output_rejects_runtime(self):
        self.prepare()
        self.page("python/search.html", "python/search", runtime=True)
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
        for replacement in ("/python/search 302", "/python/search?fixed=1 301", "/search#fixed 301"):
            with self.subTest(replacement=replacement):
                redirects.write_text(original.replace("/python/search 301", replacement), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, "redirects must"):
                    check_output(self.public)

    def test_output_rejects_wrong_sitemap(self):
        self.prepare()
        sitemap = self.public / "sitemap.xml"
        sitemap.write_text(
            sitemap.read_text(encoding="utf-8").replace("/python/installation", "/python/search"),
            encoding="utf-8",
        )
        with self.assertRaisesRegex(ValueError, "sitemap"):
            check_output(self.public)

    def test_output_validates_final_head_metadata(self):
        self.prepare()
        path = self.public / "python/installation.html"
        original = path.read_text(encoding="utf-8")
        for old, new, message in (
            ('name="description"', 'name="removed-description"', "description"),
            ('<title>', '<title></title><title>', "title"),
            ('name="robots"', 'name="removed-robots"', "robots"),
            ('index, follow', 'noindex, follow', "robots"),
            ('property="og:type"', 'property="removed-og:type"', "og:type"),
            ('name="twitter:card"', 'name="removed-twitter:card"', "twitter:card"),
            ('<h1>', '<h1></h1><h1>', "h1"),
            ('role="main"', 'role="article"', "SSG shell"),
        ):
            with self.subTest(old=old):
                path.write_text(original.replace(old, new), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, message):
                    check_output(self.public)

    def test_output_rejects_duplicate_titles_and_descriptions(self):
        self.prepare()
        path = self.public / "python/installation.html"
        original = path.read_text(encoding="utf-8")
        for old, new, message in (
            ("COSMolKit python/installation", "COSMolKit python/api", "duplicates.*title"),
            ("Read the COSMolKit python/installation guide and reference.", "Read the COSMolKit python/api guide and reference.", "duplicates.*description"),
        ):
            with self.subTest(message=message):
                path.write_text(original.replace(old, new), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, message):
                    check_output(self.public)

    def test_output_requires_noindex_on_utility_and_placeholder_pages(self):
        self.prepare()
        for route in ("python/search", "python/genindex", "python/py-modindex", "javascript", "benchmarks"):
            with self.subTest(route=route):
                path = self.public / f"{route}.html"
                original = path.read_text(encoding="utf-8")
                path.write_text(original.replace("noindex, follow", "index, follow"), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, "robots"):
                    check_output(self.public)
                path.write_text(original, encoding="utf-8")

    def test_sitemap_honors_rendered_noindex(self):
        self.page("index.html", "")
        path = self.page("python/installation.html", "python/installation")
        path.write_text(path.read_text().replace("index, follow", "noindex, follow"), encoding="utf-8")
        self.assertEqual(write_sitemap(self.public), 1)

    def test_output_requires_homepage_structured_data(self):
        self.prepare()
        path = self.public / "index.html"
        path.write_text(path.read_text().replace('"@type": "WebSite"', '"@type": "Article"'), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "JSON-LD"):
            check_output(self.public)

    def test_output_requires_project_links_on_copied_module_pages(self):
        self.prepare()
        module = self.page("python/_modules/index.html", "python/_modules/")
        self.assertEqual(check_output(self.public), 18)
        module.write_text(module.read_text().replace("cosmolkit-project-links", "other-links"), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Project and Rust crates"):
            check_output(self.public)

    def test_output_requires_deployed_image_and_correct_key(self):
        self.prepare()
        image = self.public / "social-card.png"
        image.unlink()
        with self.assertRaisesRegex(ValueError, "missing deployed resource"):
            check_output(self.public)
        image.write_bytes(b'not an image')
        with self.assertRaisesRegex(ValueError, "1200 x 630"):
            check_output(self.public)
        image.write_bytes(PNG_HEADER)
        (self.public / f"{INDEXNOW_KEY}.txt").write_text("wrong-key", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "IndexNow"):
            check_output(self.public)

    def test_social_image_download_validates_and_saves_entire_response(self):
        png = PNG_HEADER + b'fixture payload beyond the 24 byte header'
        image = self.public / 'social-card.png'
        with patch("prepare_deployment.urlopen") as fetch:
            response = MagicMock()
            fetch.return_value.__enter__.return_value = response
            response.status = 200
            response.url = SOCIAL_IMAGE_SOURCE
            response.headers = Message()
            response.headers['Content-Type'] = 'image/png'
            response.read.return_value = png
            download_social_image(self.public)
            self.assertEqual(image.read_bytes(), png)
            response.read.assert_called_once_with()
            request = fetch.call_args.args[0]
            self.assertEqual(request.full_url, SOCIAL_IMAGE_SOURCE)
            self.assertEqual(request.get_header('User-agent'), 'COSMolKit-Docs-CI/1.0')
            self.assertEqual(fetch.call_args.kwargs, {'timeout': 30})
            for invalid in (b'not a PNG', png[:16] + struct.pack('>II', 600, 315)):
                response.read.return_value = invalid
                with self.assertRaisesRegex(ValueError, '1200 x 630'):
                    download_social_image(self.public)
                self.assertEqual(image.read_bytes(), png)
            response.read.return_value = png
            response.headers.replace_header('Content-Type', 'text/html')
            with self.assertRaisesRegex(ValueError, 'image/png'):
                download_social_image(self.public)
            fetch.side_effect = HTTPError(SOCIAL_IMAGE_SOURCE, 403, 'Forbidden', {}, None)
            with self.assertRaises(HTTPError):
                download_social_image(self.public)
            self.assertEqual(image.read_bytes(), png)

    def test_output_rejects_missing_or_unordered_search_dependencies(self):
        self.prepare()
        path = self.public / "python/search.html"
        original = path.read_text(encoding="utf-8")
        first = SEARCH_LOADER
        second = '<script defer src="/unexpected.js"></script>'
        for changed in (
            original.replace(first, ""),
            original.replace(first, second + first),
            original.replace('type="module"', 'type="text/javascript"'),
            original.replace('id="search-results"', 'id="missing-results"'),
            original.replace('action="/python/search"', 'action="/missing"'),
        ):
            with self.subTest(changed=changed[:60]):
                path.write_text(changed, encoding="utf-8")
                with self.assertRaisesRegex(ValueError, "search.html"):
                    check_output(self.public)
        path.write_text(original, encoding="utf-8")
        (self.public / SEARCH_SCRIPTS[0].lstrip("/")).unlink()
        with self.assertRaisesRegex(ValueError, "missing deployed resource"):
            check_output(self.public)

    def test_rejects_undeclared_noindex_html(self):
        from check_ssg_output import check_output
        self.prepare()
        extra = self.page("unexpected.html", "unexpected")
        extra.write_text(extra.read_text().replace("index, follow", "noindex, follow"), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "undeclared HTML"):
            check_output(self.public)

    def test_internal_links_must_resolve_without_legacy_redirects(self):
        self.prepare()
        path = self.public / "python/molecule.html"
        original = path.read_text(encoding="utf-8")
        for href, message in (("/api", "legacy route"), ("/python/missing", "broken internal link")):
            with self.subTest(href=href):
                path.write_text(original + f'<a href="{href}">API</a>', encoding="utf-8")
                with self.assertRaisesRegex(ValueError, message):
                    check_output(self.public)

    def test_python_landing_and_children_coexist_after_flatten(self):
        self.prepare()
        self.assertTrue((self.public / "python.html").is_file())
        self.assertTrue((self.public / "python/api.html").is_file())
        self.assertTrue((self.public / SEARCH_SCRIPTS[0].lstrip("/")).is_file())
        self.assertFalse((self.public / "python/index.html").exists())
        self.assertEqual(flatten_html_routes(self.public), 0)

    def test_search_wasm_must_be_valid_and_search_only(self):
        self.prepare()
        wasm = self.public / SEARCH_SCRIPTS[1].lstrip("/")
        original = wasm.read_bytes()
        wasm.write_bytes(b"not a wasm module")
        with self.assertRaisesRegex(ValueError, "WebAssembly"):
            check_output(self.public)
        wasm.write_bytes(original)
        home = self.public / "index.html"
        home.write_text(home.read_text(encoding="utf-8") + SEARCH_LOADER, encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "only on the search page"):
            check_output(self.public)


if __name__ == "__main__":
    unittest.main()
