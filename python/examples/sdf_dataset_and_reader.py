"""COSMolKit usage: SDF dataset indexing and batch readers."""

from __future__ import annotations

import argparse
from pathlib import Path

from cosmolkit import MoleculeBatch, SdfDataset, SdfReader


SMALL_SDF = """carbon
     COSMolKit      2D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
>  <supplier_id>
example-001

$$$$
oxygen
     COSMolKit      2D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
nitrogen
     COSMolKit      2D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"""


def default_sdf_path() -> Path:
    output_dir = Path(__file__).resolve().parent / "output" / "sdf"
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "small.sdf"
    path.write_text(SMALL_SDF)
    return path


def summarize_batches(batches) -> tuple[int, int, int]:
    batch_count = 0
    record_count = 0
    valid_count = 0
    for batch in batches:
        batch_count += 1
        record_count += len(batch)
        valid_count += sum(batch.valid_mask())
    return batch_count, record_count, valid_count


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("sdf", nargs="?", type=Path)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--progress-bar", action="store_true")
    args = parser.parse_args()

    path = args.sdf if args.sdf is not None else default_sdf_path()
    batch_size = args.batch_size

    print("SDF:", path)
    print("file size bytes:", path.stat().st_size)

    dataset = SdfDataset.open(str(path), coordinate_dim="auto")
    print("indexed records:", len(dataset))

    first_meta = dataset.metadata(0)
    print("first title:", first_meta.title())
    print("first byte range:", first_meta.byte_range())
    print("first line range:", first_meta.line_range())

    first = dataset[0]
    last = dataset[-1]
    print("first record title:", first.title())
    print("last record title:", last.title())
    print("first supplier_id:", first.data_field("supplier_id"))

    head = dataset[: min(2, len(dataset))]
    print("slice batch records:", len(head))
    print("slice valid mask:", head.valid_mask())

    dataset_batches = dataset.batches(
        size=batch_size,
        errors="keep",
        progress_bar=args.progress_bar,
    )
    batch_count, record_count, valid_count = summarize_batches(dataset_batches)
    print("dataset batch count:", batch_count)
    print("dataset records read:", record_count)
    print("dataset valid records:", valid_count)

    reader_batches = SdfReader.open(str(path), coordinate_dim="auto").batches(
        size=batch_size,
        errors="keep",
    )
    batch_count, record_count, valid_count = summarize_batches(reader_batches)
    print("reader batch count:", batch_count)
    print("reader records read:", record_count)
    print("reader valid records:", valid_count)

    if path.stat().st_size <= 50 * 1024 * 1024:
        all_records = MoleculeBatch.read_sdf(
            str(path),
            errors="keep",
            progress_bar=args.progress_bar,
        )
        print("whole-file batch records:", len(all_records))
        print("whole-file valid records:", sum(all_records.valid_mask()))
    else:
        print("whole-file batch skipped for large file")


if __name__ == "__main__":
    main()
