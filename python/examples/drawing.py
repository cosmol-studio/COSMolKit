"""Generate 2D coordinates and export SVG/PNG images."""

from pathlib import Path
from tempfile import TemporaryDirectory

import cosmolkit as ck


def main() -> None:
    original = ck.mol_from_smiles("CCO")
    laid_out = original.with_2d_coordinates()
    params = ck.Coordinate2DParams(coordinate_map={0: [-0.0, 1.5]})
    configured = original.with_2d_coordinates_with_params(params)
    assert original.coordinates_2d() is None
    assert laid_out.coordinates_2d() is not None
    assert configured.coordinates_2d() is not None

    svg = laid_out.to_svg(300, 200)
    png = laid_out.to_png(300, 200)
    with TemporaryDirectory(prefix="cosmolkit-drawing-") as directory:
        svg_path = Path(directory) / "ethanol.svg"
        png_path = Path(directory) / "ethanol.png"
        laid_out.write_svg(str(svg_path), 300, 200)
        laid_out.write_png(str(png_path), 300, 200)
        assert svg_path.read_bytes() == svg.encode("utf-8")
        assert png_path.read_bytes() == png
        print(f"SVG: {len(svg.encode('utf-8'))} bytes; PNG: {len(png)} bytes")

    try:
        _ = original.to_svg(0, 200)
    except ck.DrawingError as error:
        assert error.domain == "drawing"
        assert error.kind == "InvalidDimensions"
        assert (error.width, error.height) == (0, 200)
        print(f"{error.kind}: {error}")
    else:
        raise AssertionError("zero width must fail")


if __name__ == "__main__":
    main()
