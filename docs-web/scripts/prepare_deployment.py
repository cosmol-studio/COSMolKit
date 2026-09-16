"""Prepare the final static artifact using the documentation route contract."""

from pathlib import Path
import shutil
import struct
import sys
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

from flatten_html_routes import flatten_html_routes
from strip_client_runtime import strip_client_runtime
from generate_sitemap import write_sitemap
from route_contract import PAGES, SOCIAL_IMAGE_SOURCE, redirect_rules


def validate_social_image(image):
    if len(image) < 24 or image[:16] != b"\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR" or struct.unpack(
        ">II", image[16:24]
    ) != (1200, 630):
        raise ValueError("social-card.png must be a 1200 x 630 PNG")


def download_social_image(public):
    request = Request(SOCIAL_IMAGE_SOURCE, headers={"User-Agent": "COSMolKit-Docs-CI/1.0"})
    with urlopen(request, timeout=30) as response:
        if response.status != 200 or response.headers.get_content_type() != "image/png":
            raise ValueError("social image must return HTTP 200 and image/png")
        if urlsplit(response.url).scheme != "https":
            raise ValueError("social image must remain HTTPS after redirects")
        image = response.read()
    validate_social_image(image)
    (public / "social-card.png").write_bytes(image)


def write_route_assets(public):
    rules = redirect_rules(PAGES)
    (public / "_redirects").write_text(
        "# Generated from routes.toml; do not edit.\n" + "".join(
            f"{source} {destination} {status}\n" for source, (destination, status) in sorted(rules.items())
        ), encoding="utf-8")


def prepare(public, sphinx):
    docs_web = Path(__file__).resolve().parents[1]
    # Fetch before changing the SSG layout; failures stop artifact preparation.
    download_social_image(public)
    print(f"Flattened {flatten_html_routes(public)} routes")
    print(f"Stripped runtime from {strip_client_runtime(public)} pages")
    shutil.copytree(docs_web / "deployment/public", public, dirs_exist_ok=True,
                    ignore=shutil.ignore_patterns("social-card.png"))
    namespace = public / "python"
    namespace.mkdir(exist_ok=True)
    for name in ("_static", "_sources", "_modules"):
        shutil.copytree(sphinx / name, namespace / name, dirs_exist_ok=True)
    for name in ("_images", "_downloads"):
        if (sphinx / name).is_dir():
            shutil.copytree(sphinx / name, namespace / name, dirs_exist_ok=True)
    shutil.copyfile(sphinx / "objects.inv", namespace / "objects.inv")
    for name in ("robots.txt", "283f8e40-25ad-4213-b83b-e5c6d3b3c5e0.txt"):
        shutil.copyfile(sphinx / name, public / name)
    write_route_assets(public)
    print(f"Generated sitemap with {write_sitemap(public)} pages")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: prepare_deployment.py PUBLIC_DIR SPHINX_HTML_DIR")
    prepare(Path(sys.argv[1]), Path(sys.argv[2]))
