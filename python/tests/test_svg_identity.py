"""Public SVG exports must identify COSMolKit, not the reference renderer."""

from io import StringIO
from pathlib import Path
import xml.etree.ElementTree as ET

import cosmolkit


def assert_cosmolkit_svg_identity(svg: str) -> None:
    namespaces = list(ET.iterparse(StringIO(svg), events=("start-ns",)))
    declarations = [namespace for _, namespace in namespaces]
    assert ("ck", "https://kit.cosmol.org/") in declarations, declarations
    assert all(prefix != "rdkit" for prefix, _ in declarations), declarations
    assert "http://www.rdkit.org/xml" not in svg
    assert "https://www.cosmol.org" not in svg
    assert ET.fromstring(svg).tag == "{http://www.w3.org/2000/svg}svg"


def test_to_svg_uses_cosmolkit_website_namespace():
    molecule = cosmolkit.Molecule.from_smiles("CCO")
    assert_cosmolkit_svg_identity(molecule.to_svg(300, 200))


def test_write_svg_uses_cosmolkit_website_namespace(tmp_path: Path):
    molecule = cosmolkit.Molecule.from_smiles("CCO")
    path = tmp_path / "ethanol.svg"
    molecule.write_svg(str(path), 300, 200)
    assert_cosmolkit_svg_identity(path.read_text(encoding="utf-8"))
