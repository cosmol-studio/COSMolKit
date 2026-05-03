"""COSMolKit usage: parallel batch workflows."""

from pathlib import Path

from cosmolkit import MoleculeBatch

output_dir = Path(__file__).resolve().parent / "output" / "batch"
output_dir.mkdir(parents=True, exist_ok=True)

smiles = [
    "CCO",
    "c1ccccc1O",
    "[13CH3:7][C@H](F)Cl",
    "not-a-smiles",
    "CC(=O)O",
]

batch = MoleculeBatch.from_smiles_list(smiles, errors="keep", n_jobs=4)
print("records:", len(batch))
print("valid mask:", batch.valid_mask())
print("errors:", [error.as_dict() for error in batch.errors()])

prepared = batch.add_hydrogens(errors="keep", n_jobs=4).compute_2d_coords(
    errors="keep",
    n_jobs=4,
)

canonical = prepared.to_smiles_list(n_jobs=4)
noncanonical = prepared.to_smiles_list(canonical=False, n_jobs=4)
without_maps = prepared.to_smiles_list(ignore_atom_map_numbers=True, n_jobs=4)

print("canonical:", canonical)
print("noncanonical:", noncanonical)
print("without atom maps:", without_maps)

svgs = prepared.to_svg_list(width=320, height=240, n_jobs=4)
print("svg count:", sum(svg is not None for svg in svgs))

image_report = prepared.to_images(
    str(output_dir / "images"),
    format="svg",
    size=(320, 240),
    errors="skip",
    n_jobs=4,
    report_path=str(output_dir / "image_errors.json"),
)
print("image report:", image_report)

sdf_report = prepared.to_sdf(
    str(output_dir / "molecules.sdf"),
    format="v2000",
    errors="skip",
    n_jobs=4,
    report_path=str(output_dir / "sdf_errors.json"),
)
print("sdf report:", sdf_report)
