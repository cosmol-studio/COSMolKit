"""Carry Sphinx descriptions into the Dioxus build without copying HTML heads."""

from pathlib import Path
import sys

from html_metadata import read_metadata


# Sphinx-generated utility pages have no author-written description.
UTILITY_DESCRIPTIONS = {
    "search": "Search COSMolKit Python guides and API reference for molecules, descriptors, fingerprints, and molecular workflows.",
    "genindex": "Browse the alphabetical index of COSMolKit Python classes, methods, and documentation topics.",
    "py-modindex": "Find COSMolKit Python modules and navigate to their API reference documentation.",
}


def rust_literal(value: str) -> str:
    marker = "#"
    while '"' + marker in value:
        marker += "#"
    return f'r{marker}"{value}"{marker}'


def generate_metadata(directory: Path, pages: list[str]) -> str:
    lines = [
        "pub fn sphinx_metadata(page: &str) -> &'static str {",
        "    match page {",
    ]
    for page in pages:
        metadata = read_metadata(directory / f"{page}.html")
        if not metadata.descriptions and page in UTILITY_DESCRIPTIONS:
            description = UTILITY_DESCRIPTIONS[page]
        elif len(metadata.descriptions) == 1 and metadata.descriptions[0].strip():
            description = metadata.descriptions[0]
        else:
            raise ValueError(f"{page}: expected one nonempty Sphinx description")
        lines.append(f"        {rust_literal(page)} => {rust_literal(description)},")
    lines.extend(['        _ => panic!("unknown Sphinx page: {page}"),', "    }", "}"])
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise SystemExit("usage: extract_sphinx_metadata.py HTML_DIR PAGE [PAGE ...]")
    sys.stdout.buffer.write(generate_metadata(Path(sys.argv[1]), sys.argv[2:]).encode("utf-8"))
