"""Validate release versions; manual publication is restricted to matching RCs."""

import os
from pathlib import Path
import re
import tomllib


ROOT = Path(__file__).resolve().parents[2]


def check_versions(rust: str, python: str, wasm: str, event: str, ref: str, requested: bool) -> bool:
    rc = re.fullmatch(r"(\d+\.\d+\.\d+)-rc\.(\d+)", rust)
    stable = re.fullmatch(r"\d+\.\d+\.\d+", rust)
    if rc is None and stable is None:
        raise ValueError("Publishing requires a stable or RC Rust version")
    expected_python = f"{rc[1]}rc{rc[2]}" if rc else rust
    if python != expected_python or wasm != rust:
        raise ValueError("Rust, Python and WASM release versions must match")
    if event == "push" and ref != f"refs/tags/v{rust}":
        raise ValueError("Release tag must equal v followed by the Rust version")
    if event == "workflow_dispatch" and requested and rc is None:
        raise ValueError("Manual publication is RC-only; publish stable releases by pushing a tag")
    return rc is not None


def main() -> None:
    def manifest(path: str) -> dict:
        with (ROOT / path).open("rb") as source:
            return tomllib.load(source)

    try:
        allowed = check_versions(
            manifest("Cargo.toml")["workspace"]["package"]["version"],
            manifest("python/pyproject.toml")["project"]["version"],
            manifest("wasm/Cargo.toml")["package"]["version"],
            os.environ.get("GITHUB_EVENT_NAME", "workflow_dispatch"),
            os.environ.get("GITHUB_REF", ""),
            os.environ.get("PUBLISH_REQUESTED") == "true",
        )
    except ValueError as error:
        raise SystemExit(f"::error::{error}") from error
    if output := os.environ.get("GITHUB_OUTPUT"):
        with open(output, "a", encoding="utf-8") as destination:
            destination.write(f"manual_publish_allowed={str(allowed).lower()}\n")
    print(f"Release versions match; manual RC publication allowed: {allowed}")


if __name__ == "__main__":
    main()
