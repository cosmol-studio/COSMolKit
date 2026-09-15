"""Validate the prepared Cloudflare Pages artifact, not legacy route bodies."""

from pathlib import Path
import sys
from xml.etree import ElementTree

from flatten_html_routes import ROUTES
from generate_sitemap import BASE_URL, EXCLUDED_ROUTES, SITEMAP_NAMESPACE, canonical_for


def check_output(public: Path) -> int:
    expected_urls = set()
    for route in ("", *ROUTES):
        relative = f"{route}.html" if route else "index.html"
        page = public / relative
        if route and (public / route).exists():
            raise ValueError(f"unflattened route would cause a trailing slash: {route}")
        html = page.read_text(encoding="utf-8")
        main_count = html.count("<main")
        stylesheet_count = html.count('rel="stylesheet"') + html.count("rel='stylesheet'")
        if main_count != 1 or stylesheet_count != 3:
            raise ValueError(
                f"invalid SSG shell for {relative}: "
                f"main={main_count}, stylesheets={stylesheet_count}"
            )
        if any(marker in html for marker in (
            "cosmolkit-docs-web", "hydrate_queue", "initial_dioxus_hydration_data",
        )):
            raise ValueError(f"client runtime was not stripped from {relative}")
        expected = BASE_URL + route
        if canonical_for(page) != expected:
            raise ValueError(f"expected clean canonical {expected} in {relative}")
        if route not in EXCLUDED_ROUTES:
            expected_urls.add(expected)

    rules = {}
    for line in (public / "_redirects").read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        source, destination, status = line.split()
        if source in rules:
            raise ValueError(f"duplicate redirect: {source}")
        rules[source] = (destination, status)
    expected_rules = {"/index.html": ("/", "301")}
    for route in ROUTES:
        for suffix in (".html", "/"):
            expected_rules[f"/{route}{suffix}"] = (f"/{route}", "301")
    if rules != expected_rules:
        raise ValueError("redirects must map only legacy/slash URLs to clean URLs with 301")

    root = ElementTree.parse(public / "sitemap.xml").getroot()
    urls = [node.text for node in root.findall(f"{{{SITEMAP_NAMESPACE}}}url/{{{SITEMAP_NAMESPACE}}}loc")]
    if set(urls) != expected_urls or len(urls) != len(expected_urls):
        raise ValueError("sitemap must contain exactly the indexable clean canonical URLs")
    return len(ROUTES) + 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: check_ssg_output.py PUBLIC_DIR")
    print(f"Validated {check_output(Path(sys.argv[1]))} clean SSG pages")
