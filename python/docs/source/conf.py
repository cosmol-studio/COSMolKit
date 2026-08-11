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
]

autosummary_generate = True
autodoc_typehints = "signature"
autoclass_content = "both"

templates_path = ["_templates"]
exclude_patterns = []

html_theme = "furo"

html_theme_options = {
    "source_repository": "https://github.com/cosmol-studio/COSMolKit/",
    "source_branch": "main",
    "source_directory": "python/docs/source/",
    "top_of_page_buttons": ["view"],
}
