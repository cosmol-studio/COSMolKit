import os
import sys

sys.path.insert(0, os.path.abspath("../.."))

project = "COSMolKit"
author = "COSMolKit Contributors"
html_title = "COSMolKit — Rust-native cheminformatics toolkit"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_sitemap",
]

autosummary_generate = True
autodoc_typehints = "signature"
autoclass_content = "both"

templates_path = ["_templates"]
exclude_patterns = []

html_theme = "furo"
html_baseurl = "https://kit.cosmol.org/"
html_extra_path = [
    "robots.txt",
    "283f8e40-25ad-4213-b83b-e5c6d3b3c5e0.txt",
]
html_static_path = ["_static"]
html_css_files = ["cosmolkit.css"]

sitemap_url_scheme = "{link}"
sitemap_locales = [None]
sitemap_excludes = [
    "search.html",
    "genindex.html",
    "py-modindex.html",
    "_modules/*",
]
sitemap_indent = 2
# The documentation workflow uses a shallow checkout, so Git-derived lastmod
# timestamps would not be complete or reliable.

html_theme_options = {
    "source_repository": "https://github.com/cosmol-studio/COSMolKit/",
    "source_branch": "main",
    "source_directory": "python/docs/source/",
    "top_of_page_buttons": ["view"],
}


def canonicalize_homepage(app, pagename, templatename, context, doctree):
    """Keep the deployed homepage canonical at the subdomain root."""
    del templatename, doctree
    if pagename == "index":
        context["pageurl"] = app.config.html_baseurl


def canonicalize_sitemap_homepage(app, exception):
    """Make the sitemap use the same root URL as the homepage canonical."""
    if exception is not None:
        return

    from pathlib import Path
    from xml.etree import ElementTree

    sitemap_path = Path(app.outdir, app.config.sitemap_filename)
    if not sitemap_path.is_file():
        return

    namespace = "http://www.sitemaps.org/schemas/sitemap/0.9"
    ElementTree.register_namespace("", namespace)
    tree = ElementTree.parse(sitemap_path)
    index_url = f"{app.config.html_baseurl}index.html"
    for location in tree.findall(f".//{{{namespace}}}loc"):
        if location.text == index_url:
            location.text = app.config.html_baseurl
    ElementTree.indent(tree, space="  ")
    tree.write(sitemap_path, encoding="utf-8", xml_declaration=True)


def setup(app):
    app.connect("html-page-context", canonicalize_homepage)
    app.connect("build-finished", canonicalize_sitemap_homepage)
