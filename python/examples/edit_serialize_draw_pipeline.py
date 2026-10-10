"""Edit, serialize, inspect, and draw a molecule value.

Usage:
    .venv/bin/python python/examples/edit_serialize_draw_pipeline.py

This example shows an explicit edit boundary followed by serialization,
coordinate generation, drawing parity inspection, fingerprint generation, and
file exports.
"""

from __future__ import annotations

from pathlib import Path

import cosmolkit as ck


OUTPUT_DIR = Path(__file__).resolve().parent / "output" / "edit_pipeline"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

base = ck.mol_from_smiles("c1ccccc1", sanitize=True)
builder = base.to_builder()
oxygen = builder.add_atom(ck.AtomSpec(ck.Element.O))
_ = builder.add_bond(ck.BondSpec(0, oxygen, ck.BondOrder.SINGLE))
phenol = builder.build().sanitize()

payload = phenol.to_binary()
restored = ck.mol_from_binary(payload)
prepared = restored.with_2d_coordinates()

print("base smiles:", base.to_smiles())
print("edited smiles:", phenol.to_smiles())
print("restored smiles:", restored.to_smiles())
print("binary bytes:", len(payload))
coordinates = prepared.coordinates_2d()
assert coordinates is not None
print("2d coords shape:", coordinates.shape)

atoms = prepared.atoms()
bonds = prepared.bonds()
print("atom table:")
for atom in atoms:
    print(
        "  atom",
        atom.id(),
        "Z=",
        atom.atomic_number(),
        "charge=",
        atom.formal_charge(),
        "aromatic=",
        atom.is_aromatic(),
    )

print("bond table:")
for bond in bonds:
    print(
        "  bond",
        bond.id(),
        bond.begin(),
        bond.end(),
        bond.order_name(),
        "aromatic=",
        bond.is_aromatic(),
    )

additional = ck.FingerprintAdditionalOutput()
additional.allocate_atom_counts()
additional.allocate_bit_info_map()
fingerprint = prepared.fingerprint_morgan(radius=2, fp_size=512, additional_output=additional)
print("fingerprint bits:", fingerprint.on_bits()[:12])
print("atom counts:", additional.atom_counts())
bit_info = additional.bit_info_map()
assert bit_info is not None
print("bit info size:", len(bit_info))

svg_path = OUTPUT_DIR / "phenol.svg"
png_path = OUTPUT_DIR / "phenol.png"
sdf_path = OUTPUT_DIR / "phenol.sdf"

prepared.write_svg(str(svg_path), width=420, height=300)
prepared.write_png(str(png_path), width=420, height=300)
prepared.write_sdf(str(sdf_path), format="v2000")

print("wrote:", svg_path)
print("wrote:", png_path)
print("wrote:", sdf_path)
