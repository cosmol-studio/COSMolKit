"""Selected FFI conversions and ownership; no chemistry corpus reproduction."""

import struct
import errno
import xml.etree.ElementTree as ET
from pathlib import Path
from collections.abc import Callable, Sequence
from typing import Literal, cast

import cosmolkit
import pytest


def coordinate_bits(rows: Sequence[Sequence[float]] | None) -> tuple[bytes, ...] | None:
    return None if rows is None else tuple(struct.pack(">d", v) for row in rows for v in row)


def state(molecule: cosmolkit.Molecule) -> tuple[str, int, int, tuple[bytes, ...] | None]:
    return (molecule.to_smiles(), molecule.num_atoms(), molecule.num_bonds(),
            coordinate_bits(molecule.coordinates_2d()))


def invalid_ffi_call(operation: object) -> Callable[..., object]:
    # Only malformed-input tests erase a native callable's signature. The
    # frozen AST/runtime controls verify that signature; no value is coerced.
    assert callable(operation)
    return operation


def test_all_parameter_defaults_and_immutable_fields():
    params = cosmolkit.Coordinate2DParams()
    empty_map: dict[int, list[float]] = {}
    expected: dict[str, object] = dict(coordinate_map=empty_map, canonical_orientation=False,
                    clear_existing_2d=True, flips_per_sample=0, samples=0,
                    sample_seed=0, permute_degree_four=False, force_rdkit=False,
                    use_ring_templates=False)
    assert {key: cast(object, getattr(params, key)) for key in expected} == expected
    for key in expected:
        setattr(params, key, expected[key])
        assert getattr(params, key) == expected[key]


def test_explicit_parameter_values_and_map_copy_bits():
    mapping = {0: [-0.0, 1.5]}
    params = cosmolkit.Coordinate2DParams(
        mapping, canonical_orientation=True, clear_existing_2d=False,
        flips_per_sample=7, samples=9, sample_seed=-13, permute_degree_four=True,
        force_rdkit=True, use_ring_templates=True)
    mapping[0][0] = 42.0
    returned = params.coordinate_map
    assert coordinate_bits([returned[0]]) == coordinate_bits([[-0.0, 1.5]])
    returned[0][0] = 23.0
    returned[1] = [4.0, 5.0]
    assert params.coordinate_map.keys() == {0}
    assert coordinate_bits([params.coordinate_map[0]]) == coordinate_bits([[-0.0, 1.5]])
    assert (params.canonical_orientation, params.clear_existing_2d,
            params.flips_per_sample, params.samples, params.sample_seed,
            params.permute_degree_four, params.force_rdkit,
            params.use_ring_templates) == (True, False, 7, 9, -13, True, True, True)


@pytest.mark.parametrize("mapping", [{0: [1.0]}, {0: [1.0, 2.0, 3.0]}])
def test_coordinate_map_requires_exact_two_element_rows(mapping: dict[int, list[float]]) -> None:
    with pytest.raises(ValueError):
        _ = cosmolkit.Coordinate2DParams(mapping)


@pytest.mark.parametrize("name,value", [
    ("flips_per_sample", -1), ("flips_per_sample", 2**32),
    ("samples", -1), ("samples", 2**32),
    ("sample_seed", -(2**31)-1), ("sample_seed", 2**31),
])
def test_parameter_integer_extraction(name: Literal["flips_per_sample", "samples", "sample_seed"], value: int) -> None:
    with pytest.raises(OverflowError):
        _ = invalid_ffi_call(cosmolkit.Coordinate2DParams)(**{name: value})


@pytest.mark.parametrize("key", [-1, 2 ** (8 * struct.calcsize("P"))])
def test_map_key_usize_extraction(key: int) -> None:
    with pytest.raises(OverflowError):
        _ = cosmolkit.Coordinate2DParams({key: [0.0, 0.0]})


def test_default_and_configured_layout_are_value_transforms():
    short_source = cosmolkit.Molecule.from_smiles("CCO")
    configured_source = cosmolkit.Molecule.from_smiles("CCO")
    before = ("CCO", 3, 2, None)
    assert state(short_source) == state(configured_source) == before
    short = short_source.with_2d_coordinates()
    assert state(short_source) == before
    configured = configured_source.with_2d_coordinates_with_params(cosmolkit.Coordinate2DParams())
    assert state(configured_source) == before
    assert coordinate_bits(short.coordinates_2d()) == coordinate_bits(configured.coordinates_2d())
    rows = short.coordinates_2d()
    assert rows is not None
    assert len(rows) == 3


def test_explicit_two_atom_map_forwarding_preserves_original():
    source = cosmolkit.Molecule.from_smiles("CCO")
    before = state(source)
    params = cosmolkit.Coordinate2DParams({0: [-0.0, 0.0], 1: [1.5, 0.0]}, force_rdkit=True)
    positioned = source.with_2d_coordinates_with_params(params)
    assert state(source) == before == ("CCO", 3, 2, None)
    assert positioned.num_atoms() == 3 and positioned.num_bonds() == 2
    rows = positioned.coordinates_2d()
    assert rows is not None
    assert len(rows) == 3
    assert coordinate_bits([params.coordinate_map[0]]) == coordinate_bits([[-0.0, 0.0]])


@pytest.mark.parametrize("method", ["to_svg", "to_png"])
@pytest.mark.parametrize("dimensions", [(0, 80), (120, 0)])
def test_drawing_error_typed_payload_and_absent_cause(method: Literal["to_svg", "to_png"], dimensions: tuple[int, int]) -> None:
    source = cosmolkit.Molecule.from_smiles("CCO").with_2d_coordinates()
    before = state(source)
    with pytest.raises(cosmolkit.DrawingError) as caught:
        _ = (source.to_svg if method == "to_svg" else source.to_png)(*dimensions)
    assert state(source) == before
    error = caught.value
    assert isinstance(error, ValueError)
    assert error.domain == "drawing" and error.kind == "InvalidDimensions"
    assert (error.width, error.height) == dimensions
    assert str(error) == f"invalid drawing dimensions {dimensions[0]}x{dimensions[1]}"
    assert error.__cause__ is None


def test_operation_error_keeps_actual_source_chain():
    source = cosmolkit.Molecule.from_smiles("CCO")
    before = state(source)
    params = cosmolkit.Coordinate2DParams({0: [0.0, 0.0], 99: [1.5, 0.0]}, force_rdkit=True)
    with pytest.raises(cosmolkit.OperationError) as caught:
        _ = source.with_2d_coordinates_with_params(params)
    assert state(source) == before
    error = caught.value
    assert isinstance(error, ValueError)
    assert error.domain == "operation" and error.kind == "Coordinate2D"
    assert str(error) == "2D coordinate generation failed: fragment layout failed: atom 99 is out of range for 3 atoms"
    assert isinstance(error.__cause__, ValueError)
    assert str(error.__cause__) == "fragment layout failed: atom 99 is out of range for 3 atoms"
    assert isinstance(error.__cause__.__cause__, ValueError)
    assert str(error.__cause__.__cause__) == "atom 99 is out of range for 3 atoms"
    assert error.__cause__.__cause__.__cause__ is None


@pytest.mark.parametrize("suffix", ["svg", "png"])
@pytest.mark.parametrize("dimensions", [(120, 80), (300, 200)])
def test_file_methods_return_none_and_exact_bytes(tmp_path: Path, suffix: Literal["svg", "png"], dimensions: tuple[int, int]) -> None:
    source = cosmolkit.Molecule.from_smiles("CCO").with_2d_coordinates()
    before = state(source)
    rendered = source.to_svg(*dimensions).encode("utf-8") if suffix == "svg" else source.to_png(*dimensions)
    assert state(source) == before
    path = tmp_path / f"ethanol.{suffix}"
    result = (source.write_svg if suffix == "svg" else source.write_png)(str(path), *dimensions)
    assert state(source) == before
    assert result is None
    assert path.read_bytes() == rendered


@pytest.mark.parametrize("suffix", ["svg", "png"])
@pytest.mark.parametrize("destination,expected_errno", [("missing/file", errno.ENOENT), ("directory", errno.EISDIR)])
def test_file_oserror_retains_errno_and_path(tmp_path: Path, suffix: Literal["svg", "png"], destination: Literal["missing/file", "directory"], expected_errno: int) -> None:
    source = cosmolkit.Molecule.from_smiles("CCO")
    before = state(source)
    path = tmp_path / destination
    if destination == "directory":
        path.mkdir()
    with pytest.raises(OSError) as caught:
        _ = (source.write_svg if suffix == "svg" else source.write_png)(str(path), 120, 80)
    assert state(source) == before
    assert caught.value.errno == expected_errno
    assert cast(object, caught.value.filename) == str(path)


@pytest.mark.parametrize("suffix", ["svg", "png"])
def test_render_error_never_opens_file(tmp_path: Path, suffix: Literal["svg", "png"]) -> None:
    source = cosmolkit.Molecule.from_smiles("CCO")
    before = state(source)
    path = tmp_path / f"invalid.{suffix}"
    with pytest.raises(cosmolkit.DrawingError):
        _ = (source.write_svg if suffix == "svg" else source.write_png)(str(path), 0, 80)
    assert state(source) == before
    assert not path.exists()
    _ = path.write_bytes(b"existing sentinel")
    with pytest.raises(cosmolkit.DrawingError):
        _ = (source.write_svg if suffix == "svg" else source.write_png)(str(path), 120, 0)
    assert state(source) == before
    assert path.read_bytes() == b"existing sentinel"


@pytest.mark.parametrize("suffix", ["svg", "png"])
def test_file_tilde_expansion(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, suffix: Literal["svg", "png"]) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))
    source = cosmolkit.Molecule.from_smiles("CCO")
    before = state(source)
    rendered = source.to_svg(120, 80).encode() if suffix == "svg" else source.to_png(120, 80)
    assert state(source) == before
    assert (source.write_svg if suffix == "svg" else source.write_png)(f"~/ethanol.{suffix}", 120, 80) is None
    assert state(source) == before
    assert (tmp_path / f"ethanol.{suffix}").read_bytes() == rendered


@pytest.mark.parametrize("smiles", ["CCO", "c1ccccc1"])
@pytest.mark.parametrize("dimensions", [(120, 80), (300, 200)])
@pytest.mark.parametrize("positioned", [False, True])
def test_svg_png_exact_repeats_and_observable_state(smiles: str, dimensions: tuple[int, int], positioned: bool) -> None:
    source = cosmolkit.Molecule.from_smiles(smiles)
    if positioned:
        source = source.with_2d_coordinates()
    before = state(source)
    svg_outputs: list[str] = []
    png_outputs: list[bytes] = []
    for _ in range(2):
        assert state(source) == before
        svg = source.to_svg(*dimensions)
        assert state(source) == before
        assert type(svg) is str
        assert svg.encode()[:6] == b"<?xml "
        root = ET.fromstring(svg)
        assert root.tag == "{http://www.w3.org/2000/svg}svg"
        assert (root.attrib["width"], root.attrib["height"]) == tuple(f"{v}px" for v in dimensions)
        assert "xmlns:ck='https://kit.cosmol.org/'" in svg
        assert "xmlns:rdkit" not in svg
        assert "www.rdkit.org" not in svg
        svg_outputs.append(svg)
        assert state(source) == before
        png = source.to_png(*dimensions)
        assert state(source) == before
        assert type(png) is bytes
        assert png[:8] == b"\x89PNG\r\n\x1a\n" and png[12:16] == b"IHDR"
        assert struct.unpack(">II", png[16:24]) == dimensions
        png_outputs.append(png)
    assert svg_outputs[0] == svg_outputs[1]
    assert png_outputs[0] == png_outputs[1]
    if positioned:
        returned = source.coordinates_2d()
        assert returned is not None
        returned[0][0] = 999.0
        assert state(source) == before


@pytest.mark.parametrize("method", ["to_svg", "to_png", "write_svg", "write_png"])
@pytest.mark.parametrize("arguments,category", [
    ((-1, 80), OverflowError), ((2**32, 80), OverflowError),
    (("120", 80), TypeError), ((), TypeError), ((120,), TypeError),
])
def test_required_dimensions_and_native_extraction(tmp_path: Path, method: Literal["to_svg", "to_png", "write_svg", "write_png"], arguments: tuple[object, ...], category: type[OverflowError] | type[TypeError]) -> None:
    source = cosmolkit.Molecule.from_smiles("CCO").with_2d_coordinates()
    before = state(source)
    path = tmp_path / "unopened"
    args = (str(path), *arguments) if method.startswith("write_") else arguments
    with pytest.raises(category):
        _ = invalid_ffi_call(cast(object, getattr(source, method)))(*args)
    assert state(source) == before
    assert not path.exists()


def test_constructor_parse_error_remains_value_error():
    with pytest.raises(ValueError) as caught:
        _ = cosmolkit.Molecule.from_smiles("(")
    assert not isinstance(caught.value, (cosmolkit.DrawingError, cosmolkit.OperationError))
