#!/usr/bin/env python3
"""Prepare pinned-Gemmi structural-writer expected data."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
import tempfile
from pathlib import Path


EXPECTED_GEMMI_COMMIT = "5cc1c23c6007e0e6cbd69289c6f7c0bff50e943e"
EXPECTED_GEMMI_VERSION = "0.7.5"
PROFILE_NAME = "bio_mmcif_writer"
SUITE_NAME = "mmcif_writer"
PROFILE_PATH = Path("testdata/bio/gemmi_mmcif_writer_profile.json")
OUTPUT_DIR = Path("testdata/bio/expected/gemmi/bio_mmcif_writer")

# The generator intentionally does not invoke Git. These hashes bind the
# executable source path to the pinned commit while respecting the port plan's
# no-Git execution rule.
GEMMI_SOURCE_HASHES = {
    "third_party/gemmi/CMakeLists.txt": "0b8552d75f795ed6af740a42ce43740f6cb6315d865016158c5ea0cf3dbe9aae",
    "third_party/gemmi/include/gemmi/cifdoc.hpp": "2622215dbbd7e093285877fec7e36a1986326e5f1da3a33f2116a3fed2cd9b34",
    "third_party/gemmi/include/gemmi/to_cif.hpp": "8d27f60b91295c5703ae949f4c0d09df9f539cb90ed85b3e2336fa32b8e17866",
    "third_party/gemmi/include/gemmi/to_mmcif.hpp": "d2388a87f3b05ff6477955fc28e13c27dc191a92039c2ffb62fa2741480b7e7b",
    "third_party/gemmi/include/gemmi/version.hpp": "0aad13cb253dccaada4819b5593b5cc88e49590d1775dfbbf62fb2b4ac284c1f",
    "third_party/gemmi/src/mmcif.cpp": "676d3b6fbed822bbe6084cfa255aba828f1ad8b54e4a8019cd9b67c24dc30a5a",
    "third_party/gemmi/src/pdb.cpp": "ab64447896fc7575537e0fe22f3d984b9c9d7f1b2b20da10a3607133fb2446e6",
    "third_party/gemmi/src/to_mmcif.cpp": "d747726769c5941a026575e929ca690aaa9049187b534c6681a5bc5fa969aade",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(64 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_identity(root: Path, path: Path) -> dict[str, str]:
    absolute = path if path.is_absolute() else root / path
    relative = absolute.resolve().relative_to(root.resolve())
    return {"path": relative.as_posix(), "sha256": sha256(absolute)}


def data_line_count(path: Path) -> int:
    return sum(
        1
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    )


def require_pinned_sources(root: Path) -> None:
    for relative, expected in GEMMI_SOURCE_HASHES.items():
        path = root / relative
        actual = sha256(path)
        if actual != expected:
            raise SystemExit(
                f"Gemmi source identity mismatch for {relative}: "
                f"expected {expected}, got {actual}; expected commit "
                f"{EXPECTED_GEMMI_COMMIT}"
            )


def build_oracle(root: Path, cmake: str, cxx: str, jobs: int) -> Path:
    source = root / "third_party/gemmi"
    build = root / "target/gemmi-mmcif-writer-golden"
    subprocess.run(
        [
            cmake,
            "-S",
            str(source),
            "-B",
            str(build),
            "-DONLY_PROGRAM=ON",
            "-DBUILD_GEMMI_PROGRAM=OFF",
            "-DCMAKE_BUILD_TYPE=Release",
        ],
        check=True,
    )
    subprocess.run(
        [cmake, "--build", str(build), "--target", "gemmi_cpp", "-j", str(jobs)],
        check=True,
    )
    binary = build / ("mmcif_writer_oracle.exe" if os.name == "nt" else "mmcif_writer_oracle")
    subprocess.run(
        [
            cxx,
            "-std=c++14",
            "-O2",
            f"-I{source / 'include'}",
            str(root / "tools/testdata/gemmi/mmcif_writer_oracle.cpp"),
            str(build / "libgemmi_cpp.a"),
            "-o",
            str(binary),
        ],
        check=True,
    )
    if not binary.is_file():
        raise SystemExit(f"Gemmi converter was not built at {binary}")
    return binary


def load_profile(root: Path) -> dict[str, object]:
    profile = json.loads((root / PROFILE_PATH).read_text(encoding="utf-8"))
    if profile.get("schema_version") != 1:
        raise SystemExit("unsupported Gemmi writer profile schema")
    if profile.get("profile") != PROFILE_NAME:
        raise SystemExit("Gemmi writer profile name mismatch")
    if profile.get("gemmi_commit") != EXPECTED_GEMMI_COMMIT:
        raise SystemExit("Gemmi writer profile commit mismatch")
    return profile


def prepare(root: Path, binary: Path, profile: dict[str, object]) -> None:
    output_dir = root / OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)
    profile_identity = file_identity(root, PROFILE_PATH)
    generator_path = Path("tools/testdata/gemmi/generate_all.py")
    generator_identities = [
        file_identity(root, generator_path),
        file_identity(root, Path("tools/testdata/gemmi/mmcif_writer_oracle.cpp")),
    ] + [
        file_identity(root, Path(relative)) for relative in GEMMI_SOURCE_HASHES
    ]
    outputs: list[dict[str, object]] = []

    with tempfile.TemporaryDirectory(
        prefix="gemmi-mmcif-writer-", dir=root / "target"
    ) as temporary:
        temporary_dir = Path(temporary)
        for raw_case in profile.get("cases", []):
            if not isinstance(raw_case, dict):
                raise SystemExit("Gemmi writer profile case must be an object")
            input_path = Path(str(raw_case["input"]))
            output_name = str(raw_case["output"])
            input_format = str(raw_case["input_format"])
            arguments = [str(value) for value in raw_case.get("arguments", [])]
            source = root / input_path
            staged = temporary_dir / output_name
            if input_format not in {"pdb", "mmcif"} or arguments:
                raise SystemExit(f"unsupported Gemmi writer case: {raw_case}")
            command = [str(binary), str(source), str(staged)]
            subprocess.run(command, check=True)
            destination = output_dir / output_name
            os.replace(staged, destination)
            outputs.append(
                {
                    "path": output_name,
                    "output_schema_version": 1,
                    "generator": generator_identities,
                    "inputs": [profile_identity, file_identity(root, source)],
                    "options": {
                        "profile": PROFILE_NAME,
                        "arguments": [
                            "--input",
                            PROFILE_PATH.as_posix(),
                            "--output",
                            output_name,
                        ],
                    },
                    "platform": {
                        "system": platform.system(),
                        "machine": platform.machine(),
                    },
                    "sha256": sha256(destination),
                    "records": data_line_count(destination),
                }
            )

    manifest = {
        "schema_version": 1,
        "family": "gemmi",
        "domain": "bio",
        "profile": PROFILE_NAME,
        "reference_implementation": {
            "name": "gemmi",
            "version": EXPECTED_GEMMI_VERSION,
        },
        "reference_runtime": {
            "implementation": "gemmi-cpp",
            "version": EXPECTED_GEMMI_VERSION,
        },
        "input": profile_identity,
        "outputs": outputs,
    }
    staged_manifest = output_dir / "manifest.json.tmp"
    staged_manifest.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(staged_manifest, output_dir / "manifest.json")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", default=PROFILE_NAME, choices=[PROFILE_NAME])
    parser.add_argument("--suite", default=SUITE_NAME, choices=[SUITE_NAME])
    parser.add_argument("--jobs", type=int, default=max(1, os.cpu_count() or 1))
    parser.add_argument("--cmake", default=shutil.which("cmake") or "cmake")
    parser.add_argument("--cxx", default=shutil.which("c++") or "c++")
    args = parser.parse_args()
    if args.jobs < 1:
        raise SystemExit("--jobs must be at least 1")

    root = Path(__file__).resolve().parents[3]
    require_pinned_sources(root)
    profile = load_profile(root)
    binary = build_oracle(root, args.cmake, args.cxx, args.jobs)
    prepare(root, binary, profile)


if __name__ == "__main__":
    main()
