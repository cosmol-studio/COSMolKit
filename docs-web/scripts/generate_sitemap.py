"""Generate a sitemap from the canonical URLs in the static Dioxus site."""

from __future__ import annotations

import sys
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urlparse
from xml.etree import ElementTree

BASE_URL = "https://kit.cosmol.org/"
SITEMAP_NAMESPACE = "http://www.sitemaps.org/schemas/sitemap/0.9"
EXCLUDED_FILES = {"search.html", "genindex.html", "py-modindex.html"}


class CanonicalParser(HTMLParser):
    """Read the canonical link without depending on the page renderer."""

    def __init__(self) -> None:
        super().__init__()
        self.canonicals: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag != "link":
            return
        attributes = dict(attrs)
        if attributes.get("rel") == "canonical" and attributes.get("href"):
            self.canonicals.append(attributes["href"])


def iter_public_pages(public_dir: Path) -> list[Path]:
    pages = []
    for path in sorted(public_dir.rglob("*.html")):
        relative = path.relative_to(public_dir)
        if relative.parts[0] == "_modules":
            continue
        if len(relative.parts) == 1 and relative.name in EXCLUDED_FILES:
            continue
        pages.append(path)
    return pages


def canonical_for(path: Path) -> str:
    parser = CanonicalParser()
    parser.feed(path.read_text(encoding="utf-8"))
    if len(parser.canonicals) != 1:
        raise SystemExit(
            f"{path}: expected one canonical link, found {len(parser.canonicals)}"
        )
    canonical = parser.canonicals[0]
    parsed = urlparse(canonical)
    if parsed.scheme != "https" or parsed.netloc != "kit.cosmol.org":
        raise SystemExit(f"{path}: canonical URL is outside {BASE_URL}: {canonical}")
    return canonical


def write_sitemap(public_dir: Path) -> int:
    canonicals = sorted({canonical_for(path) for path in iter_public_pages(public_dir)})
    if not canonicals:
        raise SystemExit("no canonical pages found in the static public directory")

    ElementTree.register_namespace("", SITEMAP_NAMESPACE)
    root = ElementTree.Element(f"{{{SITEMAP_NAMESPACE}}}urlset")
    for canonical in canonicals:
        url = ElementTree.SubElement(root, f"{{{SITEMAP_NAMESPACE}}}url")
        ElementTree.SubElement(url, f"{{{SITEMAP_NAMESPACE}}}loc").text = canonical

    ElementTree.indent(root, space="  ")
    ElementTree.ElementTree(root).write(
        public_dir / "sitemap.xml",
        encoding="utf-8",
        xml_declaration=True,
    )
    return len(canonicals)


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit(f"usage: {Path(sys.argv[0]).name} PUBLIC_DIR")
    public_dir = Path(sys.argv[1]).resolve()
    if not public_dir.is_dir():
        raise SystemExit(f"public directory does not exist: {public_dir}")
    print(f"Generated sitemap.xml with {write_sitemap(public_dir)} canonical pages")


if __name__ == "__main__":
    main()
