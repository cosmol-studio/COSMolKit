#!/usr/bin/env python3
"""Preparation-only transport to the pinned native Gemmi coordinate oracle.

Inputs are the Rust registry's ordered recipes. The native C++ oracle owns
parsing, writing and the declared seven-record projection. This adapter only
frames input paths and decodes its byte-escaped transport.
"""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

REFERENCE_LIBRARY = "gemmi"
REFERENCE_VERSION = "0.7.5"
REFERENCE_COMMIT = "5cc1c23c6007e0e6cbd69289c6f7c0bff50e943e"
PROFILE_FIELDS = ("ter_records", "numbered_ter", "ter_ignores_type",
                  "preserve_serial", "end_record")


def decode_output(transport: bytes) -> str:
    """Decode escape_bytes from pdb_coordinate_oracle.cpp without trimming PDB."""
    if not transport.endswith(b"\n") or transport.count(b"\n") != 1:
        raise ValueError("native oracle must return exactly one transport row")
    size, escaped = transport[:-1].split(b"\t", 1)
    result = bytearray()
    index = 0
    while index < len(escaped):
        value = escaped[index]
        index += 1
        if value != ord("\\"):
            result.append(value)
            continue
        escape = escaped[index]
        index += 1
        if escape == ord("x"):
            result.append(int(escaped[index:index + 2], 16))
            index += 2
        elif escape in (ord("n"), ord("t"), ord("\\")):
            result.append({ord("n"): 10, ord("t"): 9, ord("\\"): 92}[escape])
        else:
            raise ValueError("unknown native transport escape")
    if len(result) != int(size):
        raise ValueError("native output byte count mismatch")
    return result.decode("utf-8")


def generate(inputs: list[dict], oracle: Path) -> list[dict]:
    records = []
    with tempfile.TemporaryDirectory(prefix="gemmi-bio-reference-") as directory:
        directory = Path(directory)
        source = directory / "structure"
        request = directory / "request.in"
        output = directory / "reference.out"
        for original in inputs:
            row = original["BioPdbOutput"]
            case, profile = row["case"], row["profile"]
            mode = case["format"]
            if mode not in ("pdb", "cif"):
                raise ValueError(f"unknown explicit BIO input format: {mode}")
            source.write_text(case["text"], encoding="utf-8")
            flags = " ".join(str(int(profile[field])) for field in PROFILE_FIELDS)
            request.write_text(f"{source} {flags}\n", encoding="utf-8")
            result = subprocess.run([str(oracle), mode, str(request), str(output)],
                                    capture_output=True, timeout=30)
            if result.returncode != 0:
                # Preparation fails visibly; no default text or invented error row.
                raise RuntimeError(f"native Gemmi {case['id']} failed: "
                                   f"{result.stderr.decode('utf-8', errors='replace')}")
            text = decode_output(output.read_bytes())
            records.append({"input": original,
                            "output": {"BioPdbOutput": {"text": text, "error": None}}})
    return records


if __name__ == "__main__":
    json.dump(generate(json.load(sys.stdin), Path(sys.argv[1]).resolve()), sys.stdout)
