"""Project the Rust facade's feature declarations; never define another tree."""

from pathlib import Path
import argparse
import json
import re
import tomllib

ROOT = Path(__file__).resolve().parents[3]
MANIFEST = ROOT / "wasm/Cargo.toml"
# The sole approved platform exclusion. Archive code/dependencies also have
# target gates in their owning Rust crates, independent of this projection.
EXCLUDED = {"cap-serialization"}


def facade_features():
    return tomllib.loads((ROOT / "crates/cosmolkit/Cargo.toml").read_text())["features"]


def projected_features():
    source = facade_features()
    return {
        name: ([f"cosmolkit/{name}"] if name != "default" else [])
        + [edge for edge in edges if edge in source and edge not in EXCLUDED]
        for name, edges in source.items()
        if name not in EXCLUDED
    }


def resolve(selected):
    source = projected_features()
    active, pending = set(), list(selected)
    while pending:
        name = pending.pop()
        if name not in active:
            active.add(name)
            pending.extend(edge for edge in source[name] if edge in source)
    return active


def sync_manifest(write=False):
    original = MANIFEST.read_text()
    table = "[features]\n# Generated from crates/cosmolkit/Cargo.toml by tools/wasm_binding/features.py.\n"
    table += "\n".join(f"{name} = {json.dumps(edges)}" for name, edges in projected_features().items()) + "\n\n"
    updated, count = re.subn(r"^\[features\]\n.*?(?=^\[)", lambda _: table, original, flags=re.M | re.S)
    if count != 1:
        raise RuntimeError("Expected one WASM feature table")
    if write:
        MANIFEST.write_text(updated)
    elif original != updated:
        raise RuntimeError("WASM feature projection is stale; run python3 wasm/tools/wasm_binding/features.py --write")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    sync_manifest(parser.parse_args().write)
