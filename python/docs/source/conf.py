import os
import sys

sys.path.insert(0, os.path.abspath("../.."))

project = "cosmolkit"
author = "COSMolKit Contributors"

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
