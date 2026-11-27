"""Receive binary64 results from the revision-specific C++ descriptor executable."""

import hashlib
import json
import struct
import subprocess
import tempfile
from collections.abc import Callable, Iterator, Sequence
from pathlib import Path


def verify_build(executable: Path, identity: dict) -> None:
    build = identity["native_build"]
    with executable.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    if digest != build["binary_sha256"]:
        raise ValueError("descriptor reference executable checksum mismatch")

    observed = json.loads(subprocess.check_output([executable, "--build-info"]))
    if observed != build["build_info"] or observed["revision"] != identity["rdkit"]:
        raise ValueError("descriptor reference executable/build identity mismatch")


def _encode_bits(value: object) -> bytes:
    if type(value) is int and 0 <= value < (1 << 64):
        return b"\x04" + struct.pack("<Q", value)
    if isinstance(value, list):
        return b"\x07" + struct.pack("<I", len(value)) + b"".join(map(_encode_bits, value))
    raise ValueError("expected unsigned binary64 bits or an ordered array of bits")


def smiles_input(line: str) -> str:
    match json.loads(line):
        case {"smiles": str(text)} as record if len(record) == 1:
            return text
        case _:
            raise ValueError('descriptor JSONL input requires one "smiles" string')


def iter_outcomes(executable: Path, task: str, profiles: Sequence[dict], smiles: Sequence[str],
                  encode_error: Callable[[str], bytes], max_line: int) -> Iterator[list[bytes]]:
    # Spool this worker's bounded input chunk; drain stdout as each molecule finishes.
    with tempfile.TemporaryFile("w+b") as requests:
        for i, text in enumerate(smiles):
            for j, profile in enumerate(profiles):
                record = {"id": [i, j], "smiles": text, "task": task, "profile": profile}
                requests.write(json.dumps(record, allow_nan=False).encode() + b"\n")
        requests.seek(0)

        with subprocess.Popen([executable], stdin=requests, stdout=subprocess.PIPE) as process:
            stream = process.stdout
            assert stream is not None
            try:
                for i in range(len(smiles)):
                    row = []
                    for j in range(len(profiles)):
                        line = stream.readline(max_line + 1)
                        if len(line) > max_line or not line.endswith(b"\n"):
                            raise ValueError(f"incomplete or oversized {task} reference at {i}/{j}")
                        result = json.loads(line)
                        if result.pop("id") != [i, j]:
                            raise ValueError(f"out-of-order {task} reference at {i}/{j}")
                        match result:
                            case {"bits": bits} if len(result) == 1:
                                row.append(b"\x00" + _encode_bits(bits))
                            case {"rejected": "parse_rejected"} if len(result) == 1:
                                row.append(b"\x01" + encode_error("parse_rejected"))
                            case {"error": str(message)} if len(result) == 1:
                                row.append(b"\x05" + encode_error(message))
                            case _:
                                raise ValueError(f"invalid {task} reference outcome at {i}/{j}")
                    yield row
                if stream.read(1) or process.wait() != 0:
                    raise ValueError(f"unexpected trailing output or unsuccessful {task} reference")
            finally:
                if process.poll() is None:
                    process.kill()
