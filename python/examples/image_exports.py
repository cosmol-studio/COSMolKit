"""COSMolKit usage: SVG string export + SVG/PNG file export."""

from pathlib import Path

from cosmolkit import Molecule

output_dir = Path(__file__).resolve().parent / "output"
output_dir.mkdir(parents=True, exist_ok=True)

mol = Molecule.from_smiles("c1ccccc1O", sanitize=True).compute_2d_coords()

svg = mol.to_svg(width=400, height=300)
print("SVG length:", len(svg))

svg_path = output_dir / "phenol.svg"
png_path = output_dir / "phenol.png"

mol.write_svg(str(svg_path), width=400, height=300)
mol.write_png(str(png_path), width=400, height=300)

print("SVG saved:", svg_path)
print("PNG saved:", png_path)
