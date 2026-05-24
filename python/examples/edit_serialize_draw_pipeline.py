"""Edit, serialize, inspect, and draw a molecule value.

Usage:
    .venv/bin/python python/examples/edit_serialize_draw_pipeline.py

This example shows an explicit edit boundary followed by serialization,
coordinate generation, drawing parity inspection, fingerprint generation, and
file exports.
"""

from __future__ import annotations

from pathlib import Path

from cosmolkit import Molecule


OUTPUT_DIR = Path(__file__).resolve().parent / "output" / "edit_pipeline"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

base = Molecule.from_smiles("c1ccccc1", sanitize=True)
editor = base.edit()
oxygen = editor.add_atom("O")
editor.add_bond(0, oxygen, order="single")
phenol = editor.commit(sanitize=True)

payload = phenol.mol_to_binary()
restored = Molecule.mol_from_binary(payload)
prepared = restored.with_2d_coords()

print("base smiles:", base.to_smiles())
print("edited smiles:", phenol.to_smiles())
print("restored smiles:", restored.to_smiles())
print("binary bytes:", len(payload))
print("2d coords shape:", prepared.coords_2d().shape)

atoms = prepared.atoms()
bonds = prepared.bonds()
print("atom table:")
for atom in atoms:
    print(
        "  atom",
        atom.idx(),
        "Z=",
        atom.atomic_num(),
        "charge=",
        atom.formal_charge(),
        "aromatic=",
        atom.is_aromatic(),
    )

print("bond table:")
for bond in bonds:
    print(
        "  bond",
        bond.idx(),
        bond.begin_atom_idx(),
        bond.end_atom_idx(),
        bond.bond_type_name(),
        "aromatic=",
        bond.is_aromatic(),
    )

fp_result = prepared.fingerprint_morgan_with_output(radius=2, n_bits=512)
additional = fp_result.additional_output()
print("fingerprint bits:", fp_result.fingerprint().on_bits()[:12])
print("atom counts:", additional.atom_counts())
print("bit info size:", len(additional.bit_info_map()))

svg_path = OUTPUT_DIR / "phenol.svg"
png_path = OUTPUT_DIR / "phenol.png"
sdf_path = OUTPUT_DIR / "phenol.sdf"

prepared.write_svg(str(svg_path), width=420, height=300)
prepared.write_png(str(png_path), width=420, height=300)
prepared.write_sdf(str(sdf_path), format="v2000")

print("wrote:", svg_path)
print("wrote:", png_path)
print("wrote:", sdf_path)
