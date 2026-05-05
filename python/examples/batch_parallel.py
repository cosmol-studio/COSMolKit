"""COSMolKit usage: parallel batch workflows."""

from pathlib import Path

from cosmolkit import BatchErrorMode, BatchErrorType, MoleculeBatch

output_dir = Path(__file__).resolve().parent / "output" / "batch"
output_dir.mkdir(parents=True, exist_ok=True)

smiles = [
    "CCO",
    "c1ccccc1O",
    "[13CH3:7][C@H](F)Cl",
    "not-a-smiles",
    "CC(=O)O",
]

batch = MoleculeBatch.from_smiles_list(
    smiles,
    errors=BatchErrorMode.KEEP,
).with_parallel_jobs(4)
print("records:", len(batch))
print("parallel jobs:", batch.parallel_jobs())
print("valid mask:", batch.valid_mask())
for error in batch.errors():
    if error.error_type() == BatchErrorType.SMILES_PARSE:
        print("parse error:", error.index(), error.error_type().name, error.message())

prepared = batch.add_hydrogens(errors=BatchErrorMode.KEEP).compute_2d_coords(
    errors=BatchErrorMode.KEEP
)

canonical = prepared.to_smiles_list()
noncanonical = prepared.to_smiles_list(canonical=False)
without_maps = prepared.to_smiles_list(ignore_atom_map_numbers=True)

print("canonical:", canonical)
print("noncanonical:", noncanonical)
print("without atom maps:", without_maps)

for idx, molecule in enumerate(prepared):
    if molecule is not None:
        print("iterated molecule:", idx, molecule.to_smiles())

valid_molecules = [molecule for molecule in prepared.to_list() if molecule is not None]
print("valid molecule count from to_list:", len(valid_molecules))
print("tail smiles:", prepared[2:].to_smiles_list())
print("valid smiles:", prepared[prepared.valid_mask()].to_smiles_list())

svgs = prepared.to_svg_list(width=320, height=240)
print("svg count:", sum(svg is not None for svg in svgs))

image_report = prepared.to_images(
    str(output_dir / "images"),
    format="svg",
    size=(320, 240),
    errors=BatchErrorMode.SKIP,
    filenames=["ethanol", "phenol.svg", "chiral", "invalid.svg", "acetate"],
    report_path=str(output_dir / "image_errors.json"),
)
print("image report:", image_report)

sdf_report = prepared.to_sdf(
    str(output_dir / "molecules.sdf"),
    format="v2000",
    errors=BatchErrorMode.SKIP,
    report_path=str(output_dir / "sdf_errors.json"),
)
print("sdf report:", sdf_report)

sdf_file_report = prepared.to_sdf_files(
    str(output_dir / "sdf_records"),
    format="v2000",
    errors=BatchErrorMode.SKIP,
    filenames=["ethanol", "phenol.sdf", "chiral", "invalid.sdf", "acetate"],
    report_path=str(output_dir / "sdf_file_errors.json"),
)
print("sdf file report:", sdf_file_report)
