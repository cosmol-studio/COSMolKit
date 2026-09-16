"""Validate the prepared Cloudflare Pages artifact, not legacy route bodies."""

import argparse
import json
from pathlib import Path
from urllib.parse import unquote, urljoin, urlsplit
from urllib.robotparser import RobotFileParser
from xml.etree import ElementTree

from route_contract import ROUTES, PAGES, SEARCH_PATH, SOCIAL_IMAGE_URL, redirect_rules
from generate_sitemap import BASE_URL, EXCLUDED_ROUTES, SITEMAP_NAMESPACE, canonical_for
from html_metadata import read_metadata
from prepare_deployment import validate_social_image

INDEXNOW_KEY = "283f8e40-25ad-4213-b83b-e5c6d3b3c5e0"
PROJECT_LINKS = {
    "https://github.com/cosmol-studio/COSMolKit",
    "https://crates.io/crates/cosmolkit",
    "https://crates.io/crates/cosmolkit-core",
    "https://crates.io/crates/cosmolkit-inchi",
    "https://crates.io/crates/cosmolkit-ringdecomposer",
    "https://docs.rs/cosmolkit/latest/cosmolkit/",
    "https://pypi.org/project/cosmolkit/",
    BASE_URL,
    "https://tools.cosmol.org/",
}
WEBSITE_JSON_LD = {
    "@context": "https://schema.org",
    "@type": "WebSite",
    "name": "COSMolKit",
    "alternateName": "COSMolKit Documentation",
    "url": BASE_URL,
}
SEARCH_BUNDLE_PREFIX = "/assets/docs_search"


def require_one(values: list[str], label: str) -> str:
    if len(values) != 1 or not values[0].strip():
        raise ValueError(f"{label}: expected one nonempty value, found {values!r}")
    return values[0].strip()


def local_resource(public: Path, url: str) -> Path:
    parsed = urlsplit(url)
    if parsed.netloc and (parsed.scheme != "https" or parsed.netloc != "kit.cosmol.org"):
        raise ValueError(f"resource must be deployed on {BASE_URL}: {url}")
    path = (public / parsed.path.lstrip("/")).resolve()
    if not path.is_relative_to(public.resolve()) or not path.is_file():
        raise ValueError(f"missing deployed resource: {url}")
    return path


def check_page_metadata(page, path: Path, route: str) -> tuple[str, str]:
    title = require_one(page.titles, f"{path.name} title")
    description = require_one(page.descriptions, f"{path.name} description")
    if "COSMolKit" not in title:
        raise ValueError(f"{path.name}: title does not identify COSMolKit")
    if page.h1_count != 1:
        raise ValueError(f"{path.name}: expected one h1, found {page.h1_count}")
    require_one(page.metas.get("robots", []), f"{path.name} robots")
    expected_robots = {"noindex" if route in EXCLUDED_ROUTES else "index", "follow"}
    if page.robots != expected_robots:
        raise ValueError(f"{path.name}: expected robots {sorted(expected_robots)}")

    if route == "":
        if len(page.json_ld) != 1 or json.loads(page.json_ld[0]) != WEBSITE_JSON_LD:
            raise ValueError("homepage must contain one correct WebSite JSON-LD block")
    elif page.json_ld:
        raise ValueError(f"{path.name}: WebSite JSON-LD is homepage-only")

    expected_meta = {
        "og:type": "website",
        "og:site_name": "COSMolKit Documentation",
        "og:title": title,
        "og:description": description,
        "og:url": BASE_URL + route,
        "og:image": SOCIAL_IMAGE_URL,
        "og:image:type": "image/png",
        "og:image:width": "1200",
        "og:image:height": "630",
        "twitter:card": "summary_large_image",
        "twitter:title": title,
        "twitter:description": description,
        "twitter:image": SOCIAL_IMAGE_URL,
    }
    for key, value in expected_meta.items():
        if page.metas.get(key) != [value]:
            raise ValueError(f"{path.name}: incorrect or missing {key}")
    alt = require_one(page.metas.get("og:image:alt", []), f"{path.name} og:image:alt")
    if page.metas.get("twitter:image:alt") != [alt]:
        raise ValueError(f"{path.name}: twitter:image:alt differs from og:image:alt")
    return title, description


def check_project_links(public: Path) -> None:
    # Preserve the former Sphinx check on both wrapped pages and copied source pages.
    for path in sorted(public.rglob("*.html")):
        if path.relative_to(public).as_posix() == "404.html":
            continue
        page = read_metadata(path)
        if page.project_link_sections != 1 or len(page.project_links) != len(PROJECT_LINKS):
            raise ValueError(f"{path}: expected one complete Project and Rust crates section")
        if set(page.project_links) != PROJECT_LINKS:
            raise ValueError(f"{path}: incorrect Project and Rust crates links")


def check_search(public: Path) -> None:
    page = read_metadata(public / (SEARCH_PATH.lstrip("/") + ".html"))
    scripts = [script for script in page.scripts if script.get("id") == "docs-search-loader"]
    if len(scripts) != 1 or any(script.get("src") for script in page.scripts):
        raise ValueError("search.html: expected one WASM module loader")
    loader = scripts[0]
    if loader.get("type") != "module" or "async" in loader:
        raise ValueError("search.html: search loader must be a deferred module")
    for key, suffix in (("data-search-bindings", ".js"), ("data-search-wasm", ".wasm")):
        resource = loader.get(key, "")
        if not resource.startswith(SEARCH_BUNDLE_PREFIX) or not resource.endswith(suffix):
            raise ValueError(f"search.html: invalid {key}")
        path = local_resource(public, resource)
        if not path.stat().st_size:
            raise ValueError(f"empty search resource: {resource}")
        if suffix == ".wasm" and path.read_bytes()[:8] != b"\x00asm\x01\x00\x00\x00":
            raise ValueError("search.html: invalid WebAssembly binary")
    for route in ("", *ROUTES):
        if "/" + route == SEARCH_PATH:
            continue
        other = read_metadata(public / (f"{route}.html" if route else "index.html"))
        if any(script.get("id") == "docs-search-loader" for script in other.scripts):
            raise ValueError(f"search WASM must load only on the search page: {route}")
    if not {"search-results", "docs-search-form", "docs-search-status"}.issubset(page.ids):
        raise ValueError("search.html: missing search results container")
    if not any(form.get("action") == SEARCH_PATH and form.get("method", "get").lower() == "get" for form in page.forms):
        raise ValueError("search.html: missing GET search form")
    if not any(field.get("name") == "q" and field.get("type") == "search" for field in page.inputs):
        raise ValueError("search.html: missing query input")


def check_crawler_assets(public: Path) -> None:
    not_found = read_metadata(public / "404.html")
    if "noindex" not in not_found.robots or "/" not in not_found.links:
        raise ValueError("404.html must contain a robots noindex meta tag and a home link")

    robots = RobotFileParser()
    lines = (public / "robots.txt").read_text(encoding="utf-8").splitlines()
    if not {"user-agent: *", "allow: /"}.issubset({line.strip().lower() for line in lines}):
        raise ValueError("robots.txt must explicitly allow general crawlers")
    robots.parse(lines)
    if robots.site_maps() != [BASE_URL + "sitemap.xml"]:
        raise ValueError("robots.txt must advertise the canonical sitemap URL")
    for agent in ("*", "Googlebot", "bingbot"):
        for route in ("", *ROUTES, "404.html"):
            if not robots.can_fetch(agent, BASE_URL + route):
                raise ValueError(f"robots.txt blocks {agent} from crawling /{route}")

    if (public / f"{INDEXNOW_KEY}.txt").read_text(encoding="utf-8").splitlines() != [INDEXNOW_KEY]:
        raise ValueError("IndexNow verification file must contain exactly its public key")


def check_output(public: Path) -> int:
    validate_social_image(local_resource(public, SOCIAL_IMAGE_URL).read_bytes())
    check_crawler_assets(public)
    expected_urls = set()
    titles = {}
    descriptions = {}
    expected_rules = redirect_rules(PAGES)
    for route in ("", *ROUTES):
        relative = f"{route}.html" if route else "index.html"
        page = public / relative
        if route and (public / route / "index.html").exists():
            raise ValueError(f"unflattened route would cause a trailing slash: {route}")
        html = page.read_text(encoding="utf-8")
        metadata = read_metadata(page)
        main_count = metadata.main_count
        stylesheet_count = len(metadata.stylesheets)
        if main_count != 1 or metadata.main_role_count != 1 or stylesheet_count != 3:
            raise ValueError(
                f"invalid SSG shell for {relative}: "
                f"main={main_count}, main roles={metadata.main_role_count}, stylesheets={stylesheet_count}"
            )
        if any(marker in html for marker in (
            "cosmolkit-docs-web", "hydrate_queue", "initial_dioxus_hydration_data",
        )):
            raise ValueError(f"client runtime was not stripped from {relative}")
        expected = BASE_URL + route
        if canonical_for(page) != expected:
            raise ValueError(f"expected clean canonical {expected} in {relative}")
        for href in metadata.links:
            target = urlsplit(urljoin(expected, href))
            if target.scheme not in ("http", "https") or target.netloc != urlsplit(BASE_URL).netloc:
                continue
            if target.path in expected_rules:
                raise ValueError(f"{relative}: internal link uses a legacy route: {href}")
            target_path = (public / unquote(target.path).lstrip("/")).resolve()
            if not target_path.is_relative_to(public.resolve()):
                raise ValueError(f"{relative}: link escapes the deployment directory: {href}")
            if target.path == "/":
                target_path = public / "index.html"
            elif not target_path.suffix:
                target_path = target_path.with_suffix(".html")
            if not target_path.is_file():
                raise ValueError(f"{relative}: broken internal link: {href}")
        title, description = check_page_metadata(metadata, page, route)
        if title in titles:
            raise ValueError(f"{relative}: duplicates {titles[title]} title")
        if description in descriptions:
            raise ValueError(f"{relative}: duplicates {descriptions[description]} description")
        titles[title] = relative
        descriptions[description] = relative
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
    if rules != expected_rules:
        raise ValueError("redirects must map only legacy/slash URLs to clean URLs with 301")

    root = ElementTree.parse(public / "sitemap.xml").getroot()
    urls = [node.text for node in root.findall(f"{{{SITEMAP_NAMESPACE}}}url/{{{SITEMAP_NAMESPACE}}}loc")]
    if set(urls) != expected_urls or len(urls) != len(expected_urls):
        raise ValueError("sitemap must contain exactly the indexable clean canonical URLs")
    allowed_html = {"index.html", "404.html", *(f"{route}.html" for route in ROUTES)}
    for path in public.rglob("*.html"):
        relative = path.relative_to(public).as_posix()
        if relative not in allowed_html and not relative.startswith("python/_modules/"):
            raise ValueError(f"undeclared HTML route in final artifact: {relative}")
    check_project_links(public)
    check_search(public)
    return len(ROUTES) + 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("public", type=Path)
    args = parser.parse_args()
    print(f"Validated SEO, links, crawler assets, and search dependencies for {check_output(args.public)} final SSG pages")
