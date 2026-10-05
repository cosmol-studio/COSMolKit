"""Manually maintained documentation versions; no release auto-discovery."""

import json
from pathlib import Path
from urllib.parse import urlsplit

DOCS_WEB = Path(__file__).resolve().parents[1]
CATALOG_PATH = DOCS_WEB / "versions.json"
CATALOG_URL = "https://kit.cosmol.org/versions.json"


def validate_catalog(catalog):
    if not isinstance(catalog, dict) or catalog.get("schema_version") != 1:
        raise ValueError("invalid docs version catalog schema")
    entries = catalog.get("versions")
    if not isinstance(entries, list) or not entries:
        raise ValueError("docs versions must be a nonempty list")
    names, urls = set(), set()
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError("invalid docs version entry")
        name, url = entry.get("version"), entry.get("url")
        if not isinstance(name, str) or not name.strip() or not isinstance(url, str):
            raise ValueError("docs version requires name and URL")
        parsed = urlsplit(url)
        if (parsed.scheme != "https" or not parsed.netloc or parsed.username
                or parsed.password or parsed.path != "/" or parsed.query or parsed.fragment
                or name in names or url in urls):
            raise ValueError("unsafe or duplicate docs version URL")
        names.add(name)
        urls.add(url)
    if entries[0] != {"version": "latest", "url": "https://kit.cosmol.org/"}:
        raise ValueError("docs catalog must start with the canonical latest site")
    return catalog


def load_catalog(path=CATALOG_PATH):
    return validate_catalog(json.loads(path.read_text(encoding="utf-8")))


def write_version_assets(public):
    load_catalog()  # Fail artifact preparation instead of publishing bad links.
    (public / "versions.json").write_bytes(CATALOG_PATH.read_bytes())
    (public / "_headers").write_bytes((DOCS_WEB / "deployment/public/_headers").read_bytes())


def check_version_assets(public):
    if load_catalog(public / "versions.json") != load_catalog():
        raise ValueError("deployed docs version catalog differs from maintained versions.json")
    headers = (public / "_headers").read_text(encoding="utf-8")
    section = headers.split("/versions.json\n", 1)
    if len(section) != 2:
        raise ValueError("versions.json requires cross-origin catalog headers")
    lines = []
    for line in section[1].splitlines():
        if line and not line[0].isspace():
            break
        lines.append(line.strip())
    if not {"Access-Control-Allow-Origin: *", "Cache-Control: no-store"} <= set(lines):
        raise ValueError("versions.json requires CORS and fresh-cache headers")
