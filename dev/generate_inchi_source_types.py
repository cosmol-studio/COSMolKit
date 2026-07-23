#!/usr/bin/env python3
"""Regenerate owned Rust source types from the configured libinchi target."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
BASE = ROOT / "third_party" / "InChI" / "INCHI-1-SRC" / "INCHI_BASE" / "src"
API = (
    ROOT
    / "third_party"
    / "InChI"
    / "INCHI-1-SRC"
    / "INCHI_API"
    / "libinchi"
    / "src"
)
WRAPPER = ROOT / "dev" / "inchi_type_inventory_wrapper.h"
GENERATOR = ROOT / "dev" / "inchi-source-type-generator" / "Cargo.toml"
OUTPUT = ROOT / "crates" / "cosmolkit-inchi" / "src" / "source_types" / "generated.rs"


def bindgen_command(input_path: Path, allowlist: str, recursive: bool) -> list[str]:
    bindgen = os.environ.get("BINDGEN") or shutil.which("bindgen")
    if bindgen is None:
        raise RuntimeError("set BINDGEN to the bindgen-cli executable")
    command = [
        bindgen,
        str(input_path),
        "--allowlist-file",
        allowlist,
        "--no-layout-tests",
        "--no-doc-comments",
        "--no-derive-debug",
        "--no-derive-copy",
        "--no-derive-default",
        "--ignore-functions",
        "--with-derive-partialeq",
    ]
    if not recursive:
        command.append("--no-recursive-allowlist")
    command.extend(
        [
            "--",
            "-std=c11",
            "-DCOMPILE_ANSI_ONLY",
            "-DTARGET_API_LIB",
            f"-I{BASE}",
            f"-I{API}",
            f"-I{API / 'ixa'}",
        ]
    )
    return command


def generate_bindings(
    input_path: Path,
    output_path: Path,
    allowlist: str,
    recursive: bool,
    environment: dict[str, str],
) -> None:
    result = subprocess.run(
        bindgen_command(input_path, allowlist, recursive),
        cwd=ROOT,
        env=environment,
        check=True,
        capture_output=True,
        text=True,
    )
    output_path.write_text(result.stdout, encoding="utf-8")


def main() -> None:
    c_files = sorted(BASE.glob("*.c")) + sorted(API.glob("*.c")) + sorted(
        (API / "ixa").glob("*.c")
    )
    if len(c_files) != 60:
        raise RuntimeError(f"expected 60 production C files, found {len(c_files)}")

    environment = os.environ.copy()
    environment.setdefault("LIBCLANG_PATH", "/usr/lib/llvm-18/lib")
    bindgen = os.environ.get("BINDGEN") or shutil.which("bindgen")
    if bindgen is None:
        raise RuntimeError("set BINDGEN to the bindgen-cli executable")
    bindgen_version = subprocess.run(
        [bindgen, "--version"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    if bindgen_version != "bindgen 0.72.1":
        raise RuntimeError(f"expected bindgen 0.72.1, found {bindgen_version}")
    with tempfile.TemporaryDirectory(prefix="cosmolkit-inchi-types-") as temporary:
        temporary_root = Path(temporary)
        generated_inputs: list[Path] = []

        header_output = temporary_root / "000_headers.rs"
        generate_bindings(WRAPPER, header_output, ".*InChI.*", True, environment)
        generated_inputs.append(header_output)

        for index, c_file in enumerate(c_files, start=1):
            local_output = temporary_root / f"{index:03d}_{c_file.stem}.rs"
            source_pattern = rf".*{re.escape(c_file.name)}"
            generate_bindings(
                c_file,
                local_output,
                source_pattern,
                False,
                environment,
            )
            generated_inputs.append(local_output)

        subprocess.run(
            [
                "cargo",
                "run",
                "--release",
                "--manifest-path",
                str(GENERATOR),
                "--",
                *(str(path) for path in generated_inputs),
                str(OUTPUT),
            ],
            cwd=ROOT,
            env=environment,
            check=True,
        )
        subprocess.run(
            ["rustfmt", "--edition", "2024", str(OUTPUT)],
            cwd=ROOT,
            env=environment,
            check=True,
        )


if __name__ == "__main__":
    main()
