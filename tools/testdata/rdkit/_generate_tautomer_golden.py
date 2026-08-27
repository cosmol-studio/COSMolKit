#!/usr/bin/env python3
"""Generate exact pinned-RDKit tautomer evidence for a supported corpus."""

from __future__ import annotations

import argparse
import json
import multiprocessing
import shutil
import tempfile
from pathlib import Path
from typing import Any

from _tautomer_oracle import (
    assert_rdkit_version,
    build_record,
    iter_input_records,
    load_profile,
)


def branch_names(input_path: Path, profile: dict[str, object]) -> list[str]:
    if input_path.suffix == ".json":
        return [str(branch["name"]) for branch in profile["branches"]]  # type: ignore[index]
    return [str(name) for name in profile["corpus_branches"]]  # type: ignore[index]


def write_shard(
    arguments: tuple[
        int,
        list[dict[str, Any]],
        dict[str, Any],
        list[str],
        Path,
    ],
) -> tuple[int, Path]:
    shard_index, records, profile, selected, output = arguments
    with output.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(
                json.dumps(
                    build_record(record, profile, selected),
                    ensure_ascii=True,
                    sort_keys=True,
                    separators=(",", ":"),
                )
                + "\n"
            )
    return shard_index, output


def generate_records(
    records: list[dict[str, Any]],
    profile: dict[str, Any],
    selected: list[str],
    output: Path,
    *,
    shards: int,
    jobs: int,
) -> None:
    if shards < 1:
        raise ValueError("shards must be positive")
    if jobs < 1:
        raise ValueError("jobs must be positive")
    output.parent.mkdir(parents=True, exist_ok=True)
    shard_count = min(shards, max(1, len(records)))
    rows_per_shard = (len(records) + shard_count - 1) // shard_count
    with tempfile.TemporaryDirectory(
        prefix=f".{output.name}.shards-", dir=output.parent
    ) as temporary:
        temporary_dir = Path(temporary)
        arguments = []
        for shard_index in range(shard_count):
            start = shard_index * rows_per_shard
            stop = min(len(records), start + rows_per_shard)
            if start >= stop:
                break
            arguments.append(
                (
                    shard_index,
                    records[start:stop],
                    profile,
                    selected,
                    temporary_dir / f"{shard_index:04}.jsonl",
                )
            )
        if jobs == 1 or len(arguments) <= 1:
            completed = [write_shard(argument) for argument in arguments]
        else:
            with multiprocessing.Pool(processes=min(jobs, len(arguments))) as pool:
                completed = list(pool.imap_unordered(write_shard, arguments, chunksize=1))
        with output.open("wb") as handle:
            for _, shard_path in sorted(completed):
                with shard_path.open("rb") as shard:
                    shutil.copyfileobj(shard, handle, length=1024 * 1024)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--shards", type=int, default=16)
    parser.add_argument("--jobs", type=int, default=1)
    args = parser.parse_args()
    assert_rdkit_version()
    profile = load_profile()
    selected = branch_names(args.input, profile)
    generate_records(
        list(iter_input_records(args.input)),
        profile,
        selected,
        args.output,
        shards=args.shards,
        jobs=args.jobs,
    )


if __name__ == "__main__":
    main()
