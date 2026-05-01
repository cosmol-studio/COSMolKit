"""COSMolKit usage: IO + basic properties."""

from cosmolkit import Molecule

mol = Molecule.from_smiles("CCO", sanitize=True).with_2d_coords()
sdf_text = mol.to_sdf_string(format="v2000")
print("SDF length:", len(sdf_text))

saved_path = mol.write_sdf_to_directory(
    "python/examples", file_name="ethanol.sdf", format="v2000"
)
print("Saved:", saved_path)

lig = Molecule.read_sdf(saved_path, sanitize=True)
print("Loaded:", lig)
