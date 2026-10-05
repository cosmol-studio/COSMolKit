#!/usr/bin/env python3
"""Prepare pinned native Canonicalize observations independently of Enumerate."""
from __future__ import annotations
import argparse
import json
import multiprocessing
from pathlib import Path
from _tautomer_oracle import assert_rdkit_version, canonicalize_branch, error_record, load_profile, parse_molecule


def build_record(arguments):
    row, smiles, branches = arguments
    base = {"schema_version": 1, "row": row, "case_id": f"smiles_5000:{row}", "smiles": smiles, "sanitize": True, "remove_hs": True}
    try:
        molecule = parse_molecule(base)
        if molecule is None:
            return {**base, "parse": {"ok": False, "error": {"type": "NullMolecule", "message": "MolFromSmiles returned None"}}, "branches": {}}
    except Exception as error:
        return {**base, "parse": {"ok": False, "error": error_record(error)}, "branches": {}}
    return {**base, "parse": {"ok": True, "error": None}, "branches": {branch["name"]: {"parameters": branch, **canonicalize_branch(molecule, branch)} for branch in branches}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=1)
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("jobs must be positive")
    assert_rdkit_version()
    profile = load_profile()
    branches = [branch for branch in profile["branches"] if branch["name"] in profile["corpus_branches"]]
    lines = args.input.read_text().splitlines()
    # Preserve each original record, including whitespace, duplicate or invalid text.
    arguments = [(row, text, branches) for row, text in enumerate(lines)]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with multiprocessing.Pool(args.jobs) as pool, args.output.open("w") as output:
        for record in pool.imap(build_record, arguments, chunksize=1):
            output.write(json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    main()
