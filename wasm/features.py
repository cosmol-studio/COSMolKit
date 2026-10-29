"""Forward the facade's feature graph without maintaining a second graph."""

import argparse
import json
from pathlib import Path
import re
import tomllib

ROOT = Path(__file__).resolve().parent


def sync(write=False):
    source = tomllib.loads((ROOT.parent / "crates/cosmolkit/Cargo.toml").read_text())["features"]
    lines = ["[features]", "# Generated from the Rust facade by features.py --write."]
    for name, edges in source.items():
        forwarded = ([] if name == "default" else [f"cosmolkit/{name}"])
        forwarded += [edge for edge in edges if edge in source]
        lines.append(f"{name} = {json.dumps(forwarded)}")
    lines.append("registry-check = []")
    lines.append("protocol-tests = []")
    manifest = ROOT / "Cargo.toml"
    original = manifest.read_text()
    updated, count = re.subn(r"^\[features\]\n.*?(?=^\[)", lambda _: "\n".join(lines) + "\n\n", original, flags=re.M | re.S)
    if count != 1:
        raise RuntimeError("Expected one binding feature table")
    if write:
        manifest.write_text(updated)
    elif updated != original:
        raise RuntimeError("Run python wasm/features.py --write to synchronize facade features")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    sync(parser.parse_args().write)
