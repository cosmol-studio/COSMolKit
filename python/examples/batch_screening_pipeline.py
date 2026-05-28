"""End-to-end batch screening workflow.

Usage:
    .venv/bin/python python/examples/batch_screening_pipeline.py

This example combines parsing, value-style preparation, fingerprints,
similarity ranking, coordinate generation, SVG export, and structured per-record
errors. Invalid records stay aligned with their input positions.
"""

from __future__ import annotations

from pathlib import Path

from cosmolkit import BatchErrorMode, Molecule, MoleculeBatch


OUTPUT_DIR = Path(__file__).resolve().parent / "output" / "screening"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

SMILES = [
    "CCO",
    "c1ccccc1O",
    "CC(=O)O",
    "N#Cc1ccccc1",
    "not-a-smiles",
    "CCN(CC)CC",
    "O=C(O)c1ccccc1",
]


def first_valid_fingerprint(batch: MoleculeBatch):
    for fingerprint in batch.fingerprint_morgan_list(radius=2, n_bits=1024):
        if fingerprint is not None:
            return fingerprint
    raise RuntimeError("batch does not contain any valid molecules")


query = Molecule.from_smiles("c1ccccc1O", sanitize=True)
query_fp = query.fingerprint_morgan(radius=2, n_bits=1024)

batch = (
    MoleculeBatch.from_smiles_list(
        SMILES,
        sanitize=True,
        errors=BatchErrorMode.KEEP,
        n_jobs=4,
    )
    .with_parallel_jobs(4)
    .with_progress_bar(False)
)

prepared = (
    batch.sanitize(errors=BatchErrorMode.KEEP)
    .with_2d_coordinates(errors=BatchErrorMode.KEEP)
)

fingerprints = prepared.fingerprint_morgan_list(radius=2, n_bits=1024)
ranked: list[tuple[float, int, str]] = []
for index, (fingerprint, smiles) in enumerate(
    zip(fingerprints, prepared.to_smiles_list(), strict=True)
):
    if fingerprint is None or smiles is None:
        continue
    ranked.append((fingerprint.tanimoto(query_fp), index, smiles))

ranked.sort(reverse=True)

print("input records:", len(prepared))
print("valid records:", prepared.valid_count())
print("invalid records:", prepared.invalid_count())
for error in prepared.errors():
    print("error:", error.index(), error.operation(), error.message())

print("top hits:")
for score, index, smiles in ranked[:5]:
    print(f"  #{index}: {score:.3f} {smiles}")

hit_indices = [index for _, index, _ in ranked[:3]]
hits = prepared[hit_indices]
print("hit smiles:", hits.to_smiles_list())

image_report = hits.to_images(
    str(OUTPUT_DIR / "hits"),
    format="svg",
    size=(360, 260),
    filenames=[f"hit_{index}" for index in hit_indices],
    errors=BatchErrorMode.KEEP,
)
print("image export:", image_report)

sdf_report = hits.to_sdf(
    str(OUTPUT_DIR / "hits.sdf"),
    format="v2000",
    errors=BatchErrorMode.KEEP,
)
print("sdf export:", sdf_report)

reference_fp = first_valid_fingerprint(hits)
print("first hit self similarity:", reference_fp.tanimoto(reference_fp))
