"""Generate a sitemap from the canonical URLs in the static Dioxus site."""

from __future__ import annotations

import sys
from pathlib import Path
from urllib.parse import urlparse
from xml.etree import ElementTree

from html_metadata import read_metadata

from route_contract import BASE_URL, EXCLUDED_ROUTES
SITEMAP_NAMESPACE = "http://www.sitemaps.org/schemas/sitemap/0.9"
EXCLUDED_FILES = {f"{route}.html" for route in EXCLUDED_ROUTES} | {
    f"{route}/index.html" for route in EXCLUDED_ROUTES
} | {"404.html"}


def iter_public_pages(public_dir: Path) -> list[Path]:
    pages = []
    for path in sorted(public_dir.rglob("*.html")):
        relative = path.relative_to(public_dir)
        if "_modules" in relative.parts:
            continue
        if relative.as_posix() in EXCLUDED_FILES:
            continue
        if read_metadata(path).robots & {"noindex", "none"}:
            continue
        pages.append(path)
    return pages


def canonical_for(path: Path) -> str:
    parser = read_metadata(path)
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
