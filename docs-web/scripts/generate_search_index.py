"""Extract search records from Sphinx content and the documentation route contract."""

import json
from html.parser import HTMLParser
from pathlib import Path
import sys

from html_metadata import read_metadata
from route_contract import PAGES


class ArticleText(HTMLParser):
    def __init__(self):
        super().__init__()
        self.active = False
        self.chunks = []

    def handle_starttag(self, tag, attrs):
        if tag == "article" and dict(attrs).get("id") == "furo-main-content":
            self.active = True

    def handle_endtag(self, tag):
        if tag == "article":
            self.active = False

    def handle_data(self, text):
        if self.active:
            self.chunks.append(text)


def generate_index(sphinx, destination):
    raw = (sphinx / "searchindex.js").read_text(encoding="utf-8").strip().removesuffix(";")
    if not raw.startswith("Search.setIndex(") or not raw.endswith(")"):
        raise ValueError("unsupported Sphinx search index format")
    index = json.loads(raw[len("Search.setIndex("):-1])
    pages = {p["docname"]: p for p in PAGES if p["binding"] == "python" and "docname" in p}
    records = {}
    for docname in index["docnames"]:
        if docname not in pages:
            raise ValueError(f"search document is absent from routes.toml: {docname}")
        page = pages[docname]
        html = sphinx / f"{docname}.html"
        parser = ArticleText()
        parser.feed(html.read_text(encoding="utf-8"))
        metadata = read_metadata(html)
        summary = metadata.descriptions[0] if metadata.descriptions else ""
        title = index["titles"][index["docnames"].index(docname)]
        records[page["path"]] = dict(title=title, url=page["path"], summary=summary, text=" ".join(parser.chunks))

    def add(doc_id, title, anchor, summary):
        docname = index["docnames"][doc_id]
        route = pages[docname]["path"]
        # The Sphinx introduction is represented by a native landing page.
        if pages[docname].get("source") != "sphinx":
            return
        url = route + ("#" + anchor if anchor else "")
        if url not in records:
            records[url] = dict(title=title, url=url, summary=summary, text=title + " " + summary)

    for prefix, objects in index["objects"].items():
        for doc_id, type_id, _priority, anchor, name in objects:
            title = prefix + "." + name if prefix else name
            kind = index["objnames"][str(type_id)]
            if anchor == "":
                anchor = title
            elif anchor == "-":
                anchor = kind[1] + "-" + title
            add(doc_id, title, anchor, kind[2] + " in " + index["titles"][doc_id])
    for title, entries in index["alltitles"].items():
        for doc_id, anchor in entries:
            add(doc_id, title, anchor, index["titles"][doc_id])
    destination.write_text(json.dumps(list(records.values()), ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    return len(records)


if __name__ == "__main__":
    print(f"Generated {generate_index(Path(sys.argv[1]), Path(sys.argv[2]))} search records")
