#!/usr/bin/env python3
"""Native Gemmi reference adapter for BIO PDB coordinate output parity.

Calls the pinned native Gemmi oracle binary (already built) to produce
reference outputs for the bio_pdb_output_pdb and bio_pdb_output_cif tasks.
Deterministic: processes cases sequentially in input order.
"""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

REFERENCE_LIBRARY = "gemmi"
REFERENCE_VERSION = "0.7.5"
REFERENCE_COMMIT = "5cc1c23c6007e0e6cbd69289c6f7c0bff50e943e"

_ORACLE = Path(__file__).parent.parent.parent / "target" / "gemmi-mmcif-writer-golden" / "pdb_coordinate_oracle_v2"


def _run_oracle(mode: str, input_text: str) -> str:
    """Call the native oracle for one case, return the raw output."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".in", delete=False) as fin:
        fin.write(input_text + "\n")
        fin.flush()
        with tempfile.NamedTemporaryFile(mode="r", suffix=".out", delete=False) as fout:
            result = subprocess.run(
                [str(_ORACLE), mode, fin.name, fout.name],
                capture_output=True,
                text=True,
                timeout=30,
            )
            if result.returncode != 0:
                raise RuntimeError(f"oracle {mode} failed: {result.stderr}")
            return Path(fout.name).read_text().strip()


def generate_bio_pdb_output_pdb(corpus, parameters, threads=8):
    """Generate reference outputs for bio_pdb_output_pdb."""
    del threads
    records = []
    for case in corpus:
        for profile in _all_profiles():
            profile_str = " ".join(str(int(b)) for b in profile)
            output = _run_oracle("pdb", case["text"] + " " + profile_str)
            records.append({
                "case_id": case["id"],
                "parameters": {"profile": profile_str},
                "output": {"text": output},
            })
    return records


def generate_bio_pdb_output_cif(corpus, parameters, threads=8):
    """Generate reference outputs for bio_pdb_output_cif."""
    del threads
    records = []
    for case in corpus:
        for profile in _all_profiles():
            profile_str = " ".join(str(int(b)) for b in profile)
            output = _run_oracle("cif", case["text"] + " " + profile_str)
            records.append({
                "case_id": case["id"],
                "parameters": {"profile": profile_str},
                "output": {"text": output},
            })
    return records


def _all_profiles():
    """All 32 five-bool profiles (ter, numbered, ignores, preserve, end)."""
    return [
        (ter, num, ign, pres, end)
        for ter in (False, True)
        for num in (False, True)
        for ign in (False, True)
        for pres in (False, True)
        for end in (False, True)
    ]


if __name__ == "__main__":
    import sys
    corpus = json.loads(sys.stdin.read()) if len(sys.argv) > 1 else []
    print(json.dumps(generate_bio_pdb_output_pdb(corpus, {}), indent=2))
