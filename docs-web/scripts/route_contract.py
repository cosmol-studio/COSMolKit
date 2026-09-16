"""Load and validate the documentation route contract shared by every consumer."""

from pathlib import Path
import re
import tomllib
from urllib.parse import urlsplit

MANIFEST = Path(__file__).resolve().parents[1] / "routes.toml"


def load_contract(path=MANIFEST):
    contract = tomllib.loads(path.read_text(encoding="utf-8"))
    base = urlsplit(contract["base_url"])
    if set(contract) != {"base_url", "social_image_source", "pages"} or base.scheme != "https" or not base.netloc or base.path != "/" or base.query or base.fragment:
        raise ValueError("invalid contract header or canonical origin")
    image = urlsplit(contract["social_image_source"])
    if image.scheme != "https" or not image.netloc or image.username or image.password or image.query or image.fragment:
        raise ValueError("social image must have a public HTTPS URL without credentials or query parameters")
    pages = contract["pages"]
    seen = {key: set() for key in ("path", "component", "identity", "source")}
    for page in pages:
        allowed = {"component", "path", "binding", "topic", "status", "indexable", "docname", "source", "title", "legacy_paths", "label", "order", "summary", "query"}
        if set(page) - allowed:
            raise ValueError(f"unknown page fields: {set(page) - allowed}")
        if "query" in page and (not isinstance(page["query"], str) or not re.fullmatch(r"[a-z][a-z0-9_]*", page["query"])):
            raise ValueError("query must name a query parameter")
        route = page["path"]
        if route != "/" and not re.fullmatch(r"/[a-z0-9-]+(?:/[a-z0-9-]+)*", route):
            raise ValueError(f"invalid canonical route: {route}")
        if not re.fullmatch(r"[A-Z][A-Za-z0-9]*", page["component"]):
            raise ValueError("invalid Rust component")
        if page["binding"] not in ("common", "python", "javascript"):
            raise ValueError("invalid binding")
        if page["binding"] != "common" and not (route == "/" + page["binding"] or route.startswith("/" + page["binding"] + "/")):
            raise ValueError(f"route outside binding namespace: {route}")
        if page["status"] not in ("ready", "placeholder") or type(page["indexable"]) is not bool:
            raise ValueError("invalid publication policy")
        if page["status"] == "placeholder" and page["indexable"]:
            raise ValueError("placeholder cannot be indexable")
        for key, value in (("path", route), ("component", page["component"]), ("identity", (page["binding"], page["topic"]))):
            if value in seen[key]:
                raise ValueError(f"duplicate {key}: {value}")
            seen[key].add(value)
        if "docname" in page:
            docname = page["docname"]
            if not re.fullmatch(r"[a-z0-9-]+(?:/[a-z0-9-]+)*", docname):
                raise ValueError("invalid Sphinx docname")
            identity = (page["binding"], docname)
            if identity in seen["source"]:
                raise ValueError("duplicate Sphinx source")
            seen["source"].add(identity)
        if page.get("source", "native") not in ("native", "sphinx"):
            raise ValueError("unknown page source")
        if page.get("source") == "sphinx" and not all(key in page for key in ("docname", "title", "label", "order")):
            raise ValueError("incomplete Sphinx page")
        if page.get("source") == "sphinx" and page["binding"] != "python":
            raise ValueError("the current Sphinx source belongs to Python")
        if "order" in page and (type(page["order"]) is not int or page["order"] < 0 or not page.get("label")):
            raise ValueError("invalid navigation entry")
    if "/" not in seen["path"]:
        raise ValueError("missing documentation home")
    redirect_rules(pages)
    return contract


def redirect_rules(pages):
    canonical = {p["path"] for p in pages}
    rules = {"/index.html": ("/", "301")}
    for page in pages:
        route = page["path"]
        aliases = [route + suffix for suffix in (".html", "/")] if route != "/" else []
        for legacy in page.get("legacy_paths", []):
            if not re.fullmatch(r"/[a-z0-9-]+(?:/[a-z0-9-]+)*", legacy):
                raise ValueError(f"invalid legacy path: {legacy}")
            aliases.extend(legacy + suffix for suffix in ("", ".html", "/"))
        for alias in aliases:
            if alias in canonical or alias in rules:
                raise ValueError(f"redirect collides with a route or alias: {alias}")
            rules[alias] = (route, "301")
    return rules


def counterpart(pages, current, binding):
    topic = "landing" if current["binding"] == "common" else current["topic"]
    return next((p for p in pages if p["binding"] == binding and p["topic"] == topic
                 and (p["status"] == "ready" or current["binding"] == "common")), None)


CONTRACT = load_contract()
PAGES = CONTRACT["pages"]
BASE_URL = CONTRACT["base_url"]
SOCIAL_IMAGE_SOURCE = CONTRACT["social_image_source"]
SOCIAL_IMAGE_URL = BASE_URL + "social-card.png"
ROUTES = tuple(p["path"].lstrip("/") for p in PAGES if p["path"] != "/")
EXCLUDED_ROUTES = {p["path"].lstrip("/") for p in PAGES if not p["indexable"]}
SEARCH_PATH = next(p["path"] for p in PAGES if (p["binding"], p["topic"]) == ("python", "search"))
