"""Representative FFI smoke checks, not a repeated chemistry parity corpus."""

import struct
import xml.etree.ElementTree as ET

import cosmolkit
import pytest


def test_current_selected_extension_is_imported():
    assert cosmolkit._binding_profile == "drawing-bindings"
    assert cosmolkit.__version__ == "0.5.0-rc.9"
    assert cosmolkit.Molecule.__module__ == "cosmolkit"


def test_svg_png_are_real_and_preserve_the_receiver():
    molecule = cosmolkit.Molecule.from_smiles("CCO")
    before = (molecule.to_smiles(), molecule.num_atoms(), molecule.num_bonds())
    assert before == ("CCO", 3, 2)
    assert molecule.coordinates_2d() is None
    for _ in range(2):
        svg = molecule.to_svg(120, 80)
        assert isinstance(svg, str)
        root = ET.fromstring(svg)
        assert root.tag == "{http://www.w3.org/2000/svg}svg"
        assert root.attrib["width"] == "120px"
        assert root.attrib["height"] == "80px"
        # Exact current Rust emitter (draw.rs::init_drawing); binding must not
        # rewrite XML branding. The older Python branding test is separate.
        assert "xmlns:rdkit='http://www.rdkit.org/xml'" in svg
        png = molecule.to_png(120, 80)
        assert type(png) is bytes
        assert png[:8] == b"\x89PNG\r\n\x1a\n"
        assert png[12:16] == b"IHDR"
        assert struct.unpack(">II", png[16:24]) == (120, 80)
        assert (molecule.to_smiles(), molecule.num_atoms(), molecule.num_bonds()) == before
        assert molecule.coordinates_2d() is None


def test_layout_value_transform_preserves_original():
    original = cosmolkit.Molecule.from_smiles("CCO")
    positioned = original.with_2d_coordinates()
    assert original.coordinates_2d() is None
    before = positioned.coordinates_2d()
    assert before is not None and len(before) == 3
    assert len(positioned.to_svg(120, 80)) > 0
    assert positioned.to_png(120, 80).startswith(b"\x89PNG\r\n\x1a\n")
    assert positioned.coordinates_2d() == before
    assert original.coordinates_2d() is None


@pytest.mark.parametrize("method", ["to_svg", "to_png"])
def test_dimension_errors_and_python_integer_conversion(method):
    molecule = cosmolkit.Molecule.from_smiles("CCO")
    operation = getattr(molecule, method)
    with pytest.raises(cosmolkit.DrawingError, match="invalid drawing dimensions 0x80"):
        operation(0, 80)
    with pytest.raises(OverflowError):
        operation(-1, 80)
    with pytest.raises(OverflowError):
        operation(2**32, 80)
    with pytest.raises(TypeError):
        operation("120", 80)
    with pytest.raises(TypeError):
        operation()
    assert molecule.coordinates_2d() is None
