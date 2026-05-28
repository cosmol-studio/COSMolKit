"""Indexed SDF processing with chunked batch transforms.

Usage:
    .venv/bin/python python/examples/sdf_indexed_processing_pipeline.py

This example writes a small local SDF, reopens it as an indexed dataset, reads
selected records and chunks, computes fingerprints and coordinates, and exports
both image and SDF outputs from the filtered batch.
"""

from __future__ import annotations

from pathlib import Path

from cosmolkit import BatchErrorMode, Molecule, MoleculeBatch, SdfDataset


OUTPUT_DIR = Path(__file__).resolve().parent / "output" / "indexed_sdf"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
SDF_PATH = OUTPUT_DIR / "library.sdf"

SEED_SMILES = [
    "CCO",
    "CC(=O)O",
    "c1ccccc1",
    "c1ccccc1O",
    "N#Cc1ccccc1",
]


def build_seed_sdf(path: Path) -> None:
    records: list[str] = []
    for idx, smiles in enumerate(SEED_SMILES):
        mol = Molecule.from_smiles(smiles, sanitize=True).with_2d_coordinates()
        record = mol.to_2d_sdf_string(format="v2000")
        record = record.replace(
            "$$$$\n",
            f">  <source_index>\n{idx}\n\n>  <input_smiles>\n{smiles}\n\n$$$$\n",
        )
        records.append(record)
    path.write_text("".join(records), encoding="utf-8")


build_seed_sdf(SDF_PATH)

dataset = SdfDataset.open(str(SDF_PATH), coordinate_dim="2d")
print("dataset path:", dataset.path())
print("records:", len(dataset))

for index in range(len(dataset)):
    metadata = dataset.metadata(index)
    record = dataset[index]
    print(
        "record:",
        index,
        "bytes=",
        metadata.byte_len(),
        "source=",
        record.data_field("source_index"),
        "smiles=",
        record.molecule().to_smiles(),
    )

selected = dataset[[0, 2, 4]]
selected = selected.with_2d_coordinates(
    errors=BatchErrorMode.KEEP,
    n_jobs=2,
)
print("selected valid mask:", selected.valid_mask())
print("selected smiles:", selected.to_smiles_list())

fps = selected.fingerprint_morgan_list(radius=2, n_bits=512, n_jobs=2)
for index, fingerprint in enumerate(fps):
    if fingerprint is not None:
        print("selected fingerprint:", index, "bits=", len(fingerprint.on_bits()))

all_chunks: list[MoleculeBatch] = []
for chunk in dataset.batches(size=2, errors=BatchErrorMode.KEEP, n_jobs=2):
    prepared = chunk.sanitize(errors=BatchErrorMode.KEEP).with_2d_coordinates(
        errors=BatchErrorMode.KEEP,
    )
    all_chunks.append(prepared)
    print("chunk:", len(prepared), prepared.to_smiles_list())

report = selected.to_images(
    str(OUTPUT_DIR / "selected_images"),
    format="png",
    size=(320, 240),
    filenames=["ethanol", "benzene", "benzonitrile"],
    errors=BatchErrorMode.KEEP,
)
print("selected image export:", report)

sdf_report = selected.to_sdf_files(
    str(OUTPUT_DIR / "selected_sdf"),
    format="v2000",
    filenames=["ethanol", "benzene", "benzonitrile"],
    errors=BatchErrorMode.KEEP,
)
print("selected sdf export:", sdf_report)
