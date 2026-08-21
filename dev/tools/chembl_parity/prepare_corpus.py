#!/usr/bin/env python3
"""Materialize deterministic ChEMBL 37 JSONL shards with a checksummed manifest."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any


CHEMBL_37_SOURCE_SHA256 = (
    "ea6181ce8dc7af41974e35b92e1febb0c9dcbe2c62f7ccc4a5d983ac19f696e7"
)
CHEMBL_37_SOURCE_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/"
    "chembl_37/chembl_37_chemreps.txt.gz"
)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stable_shard(value: str, shards: int) -> int:
    digest = hashlib.blake2b(value.encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, "little") % shards


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shards", type=int, default=128)
    parser.add_argument(
        "--allow-source-mismatch",
        action="store_true",
        help="permit a different input only for an explicitly separate audit",
    )
    args = parser.parse_args()
    if args.shards < 1:
        raise SystemExit("--shards must be positive")
    if args.output_dir.exists():
        raise SystemExit(f"output already exists: {args.output_dir}")
    source_sha256 = file_sha256(args.input)
    if source_sha256 != CHEMBL_37_SOURCE_SHA256 and not args.allow_source_mismatch:
        raise SystemExit(
            "input is not the pinned ChEMBL 37 source; "
            "use --allow-source-mismatch only for a separately documented corpus"
        )

    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    temp_dir = Path(
        tempfile.mkdtemp(prefix=args.output_dir.name + ".prepare-", dir=args.output_dir.parent)
    )
    handles = []
    counts = [0] * args.shards
    total_rows = 0
    missing_smiles = 0
    try:
        for shard in range(args.shards):
            handles.append(
                (temp_dir / f"shard-{shard:03d}.jsonl").open("w", encoding="utf-8")
            )
        with gzip.open(args.input, "rt", encoding="utf-8", newline="") as source:
            reader = csv.DictReader(source, delimiter="\t")
            required = {
                "chembl_id",
                "canonical_smiles",
                "standard_inchi",
                "standard_inchi_key",
            }
            if not required.issubset(reader.fieldnames or ()):
                missing = sorted(required - set(reader.fieldnames or ()))
                raise SystemExit(f"source is missing required columns: {missing}")
            for row in reader:
                total_rows += 1
                chembl_id = row["chembl_id"]
                if not chembl_id:
                    raise SystemExit(f"source row {total_rows} has an empty chembl_id")
                smiles = row["canonical_smiles"]
                if not smiles:
                    missing_smiles += 1
                    continue
                shard = stable_shard(chembl_id, args.shards)
                record: dict[str, Any] = {
                    "row": total_rows,
                    "chembl_id": chembl_id,
                    "smiles": smiles,
                    "standard_inchi": row["standard_inchi"],
                    "standard_inchi_key": row["standard_inchi_key"],
                }
                handles[shard].write(json.dumps(record, ensure_ascii=True) + "\n")
                counts[shard] += 1
        for handle in handles:
            handle.close()
        handles.clear()

        shard_records = []
        for shard, count in enumerate(counts):
            path = temp_dir / f"shard-{shard:03d}.jsonl"
            shard_records.append(
                {
                    "shard": shard,
                    "records": count,
                    "bytes": path.stat().st_size,
                    "sha256": file_sha256(path),
                    "file": path.name,
                }
            )
        manifest = {
            "schema": "cosmolkit-chembl-parity-corpus-v1",
            "release": "ChEMBL 37",
            "source_url": CHEMBL_37_SOURCE_URL,
            "source_file": args.input.name,
            "source_bytes": args.input.stat().st_size,
            "source_sha256": source_sha256,
            "total_source_rows": total_rows,
            "records": sum(counts),
            "missing_smiles": missing_smiles,
            "shards": args.shards,
            "selection": "all rows with nonempty canonical_smiles",
            "assignment": "blake2b-64(chembl_id) modulo shards",
            "files": shard_records,
        }
        (temp_dir / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        os.rename(temp_dir, args.output_dir)
        print(json.dumps(manifest, sort_keys=True))
    except BaseException:
        for handle in handles:
            handle.close()
        shutil.rmtree(temp_dir, ignore_errors=True)
        raise


if __name__ == "__main__":
    main()
