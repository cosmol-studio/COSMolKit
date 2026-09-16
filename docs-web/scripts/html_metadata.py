"""Read document metadata and deployment dependencies from rendered HTML."""

from html.parser import HTMLParser
from pathlib import Path


class PageMetadata(HTMLParser):
    VOID_ELEMENTS = {"area", "base", "br", "col", "embed", "hr", "img", "input", "link", "meta", "param", "source", "track", "wbr"}

    def __init__(self):
        super().__init__()
        self.in_head = False
        self.titles = []
        self.metas = {}
        self.canonicals = []
        self.json_ld = []
        self.stylesheets = []
        self.scripts = []
        self.links = []
        self.forms = []
        self.inputs = []
        self.ids = set()
        self.main_count = 0
        self.main_role_count = 0
        self.h1_count = 0
        self.project_link_sections = 0
        self.project_links = []
        self._project_depth = None
        self._title_chunks = None
        self._json_chunks = None

    @property
    def descriptions(self):
        return self.metas.get("description", [])

    @property
    def robots(self):
        return {
            directive
            for value in self.metas.get("robots", [])
            for directive in value.lower().replace(",", " ").split()
        }

    def handle_starttag(self, tag, attrs):
        values = {name: value or "" for name, value in attrs}
        if self._project_depth is not None:
            if tag not in self.VOID_ELEMENTS:
                self._project_depth += 1
            if tag == "a" and values.get("href"):
                self.project_links.append(values["href"])
        elif tag == "nav" and "cosmolkit-project-links" in values.get("class", "").split():
            self.project_link_sections += 1
            self._project_depth = 1

        if tag == "head":
            self.in_head = True
        elif self.in_head:
            if tag == "title":
                self._title_chunks = []
            elif tag == "meta":
                key = values.get("name", values.get("property", "")).lower()
                self.metas.setdefault(key, []).append(values.get("content", ""))
            elif tag == "link" and values.get("rel") == "canonical":
                self.canonicals.append(values.get("href", ""))
            elif tag == "script" and values.get("type") == "application/ld+json":
                self._json_chunks = []

        if tag == "link" and values.get("rel") == "stylesheet":
            self.stylesheets.append(values.get("href", ""))
        if tag == "script":
            self.scripts.append(values)
        if tag == "a" and values.get("href"):
            self.links.append(values["href"])
        if tag == "form":
            self.forms.append(values)
        if tag == "input":
            self.inputs.append(values)
        if values.get("id"):
            self.ids.add(values["id"])
        self.main_count += tag == "main"
        self.main_role_count += values.get("role") == "main"
        self.h1_count += tag == "h1"

    def handle_endtag(self, tag):
        if self._project_depth is not None and tag not in self.VOID_ELEMENTS:
            self._project_depth -= 1
            if self._project_depth == 0:
                self._project_depth = None
        if tag == "head":
            self.in_head = False
        elif tag == "title" and self._title_chunks is not None:
            self.titles.append("".join(self._title_chunks))
            self._title_chunks = None
        elif tag == "script" and self._json_chunks is not None:
            self.json_ld.append("".join(self._json_chunks))
            self._json_chunks = None

    def handle_startendtag(self, tag, attrs):
        self.handle_starttag(tag, attrs)
        if tag not in self.VOID_ELEMENTS:
            self.handle_endtag(tag)

    def handle_data(self, data):
        if self._title_chunks is not None:
            self._title_chunks.append(data)
        if self._json_chunks is not None:
            self._json_chunks.append(data)


def read_metadata(path: Path) -> PageMetadata:
    result = PageMetadata()
    result.feed(path.read_text(encoding="utf-8"))
    return result
