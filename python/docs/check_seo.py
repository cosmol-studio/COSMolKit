#!/usr/bin/env python3
"""Validate search-facing artifacts produced by the Sphinx HTML build."""

from __future__ import annotations

import argparse
import json
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urljoin
from xml.etree import ElementTree

BASE_URL = "https://kit.cosmol.org/"
EXCLUDED_HTML = {"genindex.html", "py-modindex.html", "search.html"}
PROJECT_LINKS = {
    "https://github.com/cosmol-studio/COSMolKit",
    "https://crates.io/crates/cosmolkit",
    "https://crates.io/crates/cosmolkit-core",
    "https://crates.io/crates/cosmolkit-inchi",
    "https://crates.io/crates/cosmolkit-ringdecomposer",
    "https://docs.rs/cosmolkit/latest/cosmolkit/",
    "https://pypi.org/project/cosmolkit/",
    "https://kit.cosmol.org/",
    "https://tools.cosmol.org/",
}


class HeadParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.in_head = False
        self.canonicals: list[str] = []
        self.descriptions: list[str] = []
        self.json_ld: list[str] = []
        self.titles: list[str] = []
        self.project_link_sections = 0
        self.project_links: list[str] = []
        self._json_chunks: list[str] | None = None
        self._project_link_depth: int | None = None
        self._title_chunks: list[str] | None = None

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        values = dict(attrs)
        if self._project_link_depth is not None:
            self._project_link_depth += 1
            if tag == "a" and values.get("href"):
                self.project_links.append(str(values["href"]))
        elif tag == "nav" and "cosmolkit-project-links" in str(
            values.get("class", "")
        ).split():
            self.project_link_sections += 1
            self._project_link_depth = 1

        if tag == "head":
            self.in_head = True
        elif self.in_head and tag == "link" and values.get("rel") == "canonical":
            if values.get("href"):
                self.canonicals.append(str(values["href"]))
        elif self.in_head and tag == "meta" and values.get("name") == "description":
            if values.get("content"):
                self.descriptions.append(str(values["content"]))
        elif self.in_head and tag == "title":
            self._title_chunks = []
        elif (
            self.in_head
            and tag == "script"
            and values.get("type") == "application/ld+json"
        ):
            self._json_chunks = []

    def handle_endtag(self, tag: str) -> None:
        if self._project_link_depth is not None:
            self._project_link_depth -= 1
            if self._project_link_depth == 0:
                self._project_link_depth = None

        if tag == "head":
            self.in_head = False
        elif tag == "script" and self._json_chunks is not None:
            self.json_ld.append("".join(self._json_chunks))
            self._json_chunks = None
        elif tag == "title" and self._title_chunks is not None:
            self.titles.append("".join(self._title_chunks))
            self._title_chunks = None

    def handle_data(self, data: str) -> None:
        if self._json_chunks is not None:
            self._json_chunks.append(data)
        if self._title_chunks is not None:
            self._title_chunks.append(data)


def public_html_pages(build_dir: Path) -> list[Path]:
    return sorted(
        path
        for path in build_dir.glob("*.html")
        if path.name not in EXCLUDED_HTML
    )


def expected_url(path: Path) -> str:
    return BASE_URL if path.name == "index.html" else urljoin(BASE_URL, path.name)


def validate_html(build_dir: Path) -> set[str]:
    expected_urls: set[str] = set()
    descriptions: dict[str, str] = {}
    titles: dict[str, str] = {}
    for path in public_html_pages(build_dir):
        parser = HeadParser()
        parser.feed(path.read_text(encoding="utf-8"))
        expected = expected_url(path)
        expected_urls.add(expected)
        assert parser.canonicals == [expected], (
            f"{path.name}: expected one canonical {expected!r}, "
            f"found {parser.canonicals!r}"
        )
        assert len(parser.descriptions) == 1, (
            f"{path.name}: expected one meta description, "
            f"found {len(parser.descriptions)}"
        )
        description = parser.descriptions[0].strip()
        assert description, f"{path.name}: meta description is empty"
        assert description not in descriptions, (
            f"{path.name}: duplicates {descriptions[description]} meta description"
        )
        descriptions[description] = path.name

        assert len(parser.titles) == 1, (
            f"{path.name}: expected one title, found {len(parser.titles)}"
        )
        title = parser.titles[0].strip()
        assert title, f"{path.name}: title is empty"
        assert "COSMolKit" in title, f"{path.name}: title does not identify COSMolKit"
        assert title not in titles, f"{path.name}: duplicates {titles[title]} title"
        titles[title] = path.name

        if path.name == "index.html":
            assert len(parser.json_ld) == 1, "homepage must contain one JSON-LD block"
            website = json.loads(parser.json_ld[0])
            assert website == {
                "@context": "https://schema.org",
                "@type": "WebSite",
                "name": "COSMolKit",
                "alternateName": "COSMolKit Documentation",
                "url": BASE_URL,
            }
        else:
            assert not parser.json_ld, f"{path.name}: WebSite JSON-LD is homepage-only"
    return expected_urls


def validate_sitemap(build_dir: Path, expected_urls: set[str]) -> None:
    sitemap_path = build_dir / "sitemap.xml"
    assert sitemap_path.is_file(), "sitemap.xml was not generated"
    root = ElementTree.parse(sitemap_path).getroot()
    namespace = {"s": "http://www.sitemaps.org/schemas/sitemap/0.9"}
    urls = [
        node.text
        for node in root.findall("s:url/s:loc", namespace)
        if node.text is not None
    ]
    assert len(urls) == len(set(urls)), "sitemap.xml contains duplicate URLs"
    assert set(urls) == expected_urls, (
        f"sitemap URL set differs: missing={sorted(expected_urls - set(urls))!r}, "
        f"extra={sorted(set(urls) - expected_urls)!r}"
    )


def validate_project_links(build_dir: Path) -> int:
    pages = sorted(build_dir.rglob("*.html"))
    assert pages, "documentation build contains no HTML pages"
    for path in pages:
        parser = HeadParser()
        parser.feed(path.read_text(encoding="utf-8"))
        relative_path = path.relative_to(build_dir)
        assert parser.project_link_sections == 1, (
            f"{relative_path}: expected one Project and Rust crates section, "
            f"found {parser.project_link_sections}"
        )
        assert len(parser.project_links) == len(PROJECT_LINKS), (
            f"{relative_path}: expected {len(PROJECT_LINKS)} project links, "
            f"found {parser.project_links!r}"
        )
        assert set(parser.project_links) == PROJECT_LINKS, (
            f"{relative_path}: Project and Rust crates links differ: "
            f"missing={sorted(PROJECT_LINKS - set(parser.project_links))!r}, "
            f"extra={sorted(set(parser.project_links) - PROJECT_LINKS)!r}"
        )
    return len(pages)


def validate_robots(build_dir: Path) -> None:
    robots_path = build_dir / "robots.txt"
    assert robots_path.is_file(), "robots.txt was not copied"
    lines = {
        line.strip()
        for line in robots_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
    assert "User-agent: *" in lines
    assert "Allow: /" in lines
    assert f"Sitemap: {BASE_URL}sitemap.xml" in lines


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("build_dir", type=Path)
    args = parser.parse_args()
    build_dir = args.build_dir.resolve()
    expected_urls = validate_html(build_dir)
    validate_sitemap(build_dir, expected_urls)
    project_link_page_count = validate_project_links(build_dir)
    validate_robots(build_dir)
    print(
        f"validated SEO artifacts for {len(expected_urls)} public pages and "
        f"project links for {project_link_page_count} generated pages"
    )


if __name__ == "__main__":
    main()
