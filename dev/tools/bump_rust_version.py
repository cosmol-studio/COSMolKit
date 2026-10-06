#!/usr/bin/env python3
"""Update explicitly listed release versions, then refresh workspace lock entries."""

import argparse
import difflib
from pathlib import Path
import subprocess
import tomllib


# File, TOML section, dependency names. No README or directory scanning.
DEPENDENCIES = [
    ("crates/cosmolkit/Cargo.toml", "dependencies",
     "cosmolkit-model cosmolkit-macros cosmolkit-core cosmolkit-descriptors "
     "cosmolkit-alignment cosmolkit-batch cosmolkit-bio cosmolkit-conformer "
     "cosmolkit-depict cosmolkit-fingerprints cosmolkit-forcefields cosmolkit-inchi "
     "cosmolkit-io cosmolkit-search cosmolkit-smiles cosmolkit-stereo cosmolkit-tautomer"),
    ("crates/cosmolkit-alignment/Cargo.toml", "dependencies", "cosmolkit-core cosmolkit-model cosmolkit-types cosmolkit-search"),
    ("crates/cosmolkit-batch/Cargo.toml", "dependencies", "cosmolkit-model"),
    ("crates/cosmolkit-bio/Cargo.toml", "dependencies", "cosmolkit-model cosmolkit-types"),
    ("crates/cosmolkit-conformer/Cargo.toml", "dependencies", "cosmolkit-model cosmolkit-core cosmolkit-forcefields cosmolkit-alignment cosmolkit-search"),
    ("crates/cosmolkit-conformer/Cargo.toml", "dev-dependencies", "cosmolkit-smiles cosmolkit-io"),
    ("crates/cosmolkit-core/Cargo.toml", "dependencies", "cosmolkit-model cosmolkit-ringdecomposer cosmolkit-types"),
    ("crates/cosmolkit-depict/Cargo.toml", "dependencies", "cosmolkit-core cosmolkit-model cosmolkit-search"),
    ("crates/cosmolkit-descriptors/Cargo.toml", "dependencies", "cosmolkit-core cosmolkit-model cosmolkit-search"),
    ("crates/cosmolkit-descriptors/Cargo.toml", "dev-dependencies", "cosmolkit-smiles"),
    ("crates/cosmolkit-fingerprints/Cargo.toml", "dependencies", "cosmolkit-core cosmolkit-model cosmolkit-search cosmolkit-stereo"),
    ("crates/cosmolkit-fingerprints/Cargo.toml", "dev-dependencies", "cosmolkit-smiles"),
    ("crates/cosmolkit-forcefields/Cargo.toml", "dependencies", "cosmolkit-core cosmolkit-model cosmolkit-search"),
    ("crates/cosmolkit-forcefields/Cargo.toml", "dev-dependencies", "cosmolkit-smiles cosmolkit-io"),
    ("crates/cosmolkit-io/Cargo.toml", "dependencies", "cosmolkit-bio cosmolkit-core cosmolkit-model cosmolkit-search cosmolkit-types"),
    ("crates/cosmolkit-model/Cargo.toml", "dependencies", "cosmolkit-types"),
    ("crates/cosmolkit-search/Cargo.toml", "dependencies", "cosmolkit-cx cosmolkit-core cosmolkit-model cosmolkit-types"),
    ("crates/cosmolkit-smiles/Cargo.toml", "dependencies", "cosmolkit-core cosmolkit-cx cosmolkit-model cosmolkit-types"),
    ("crates/cosmolkit-stereo/Cargo.toml", "dependencies", "cosmolkit-core cosmolkit-model cosmolkit-types"),
    ("crates/cosmolkit-tautomer/Cargo.toml", "dependencies", "cosmolkit-core cosmolkit-model cosmolkit-search cosmolkit-smiles cosmolkit-types"),
    ("python/Cargo.toml", "dependencies", "cosmolkit"),
    ("wasm/Cargo.toml", "dependencies", "cosmolkit"),
]
README = "crates/cosmolkit/README.md"
# Binding crate package version is separate from its cosmolkit dependency.
PACKAGE_VERSIONS = ["python/Cargo.toml"]
BEGIN = "<!-- rust-install-version:start -->"
END = "<!-- rust-install-version:end -->"
LOCK_REFRESH = ["cargo", "update", "--workspace"]


def update_field(text, section, key, version):
    """Change one named field, preserving comments and all unrelated text."""
    header = "[" + section + "]\n"
    if text.count(header) != 1:
        raise ValueError(f"expected one section: {section}")
    before, body = text.split(header)
    lines = body.splitlines(keepends=True)
    matches = []
    for index, line in enumerate(lines):
        if line.startswith("["):
            break
        if line.startswith(key + " = "):
            matches.append(index)
    if len(matches) != 1:
        raise ValueError(f"expected one field: {section}.{key}")
    index = matches[0]
    prefix, old, suffix = lines[index].split('"', 2)
    if key != "version" and prefix != key + ' = { version = ':
        raise ValueError(f"unexpected dependency layout: {section}.{key}")
    lines[index] = prefix + '"' + version + '"' + suffix
    return before + header + "".join(lines)


def prepare_updates(original, version):
    """Prepare the complete release edit set before any file is changed."""
    base, separator, rc = version.partition("-rc.")
    updated = original.copy()
    # All current publishable crate versions inherit this single field.
    updated["Cargo.toml"] = update_field(updated["Cargo.toml"], "workspace.package", "version", version)
    for file in PACKAGE_VERSIONS:
        updated[file] = update_field(updated[file], "package", "version", version)
    # Python release candidates use rcN; stable versions stay X.Y.Z.
    python_version = base + "rc" + rc if separator else base
    updated["python/pyproject.toml"] = update_field(
        updated["python/pyproject.toml"], "project", "version", python_version
    )
    for file, section, names in DEPENDENCIES:
        for name in names.split():
            updated[file] = update_field(updated[file], section, name, version)

    # Only this explicitly marked installation example; historical prose is untouched.
    text = updated[README]
    if text.count(BEGIN) != 1 or text.count(END) != 1:
        raise ValueError("README installation markers must each occur once")
    before, example = text.split(BEGIN)
    example, after = example.split(END)
    old = tomllib.loads(example.strip().removeprefix("```toml").removesuffix("```").strip())["cosmolkit"]["version"]
    token = 'version = "' + old + '"'
    if example.count(token) != 1:
        raise ValueError("expected exactly one installation version")
    updated[README] = before + BEGIN + example.replace(token, 'version = "' + version + '"') + END + after
    return updated


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("version", help="for example: 0.5.0-rc.9")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    base, separator, rc = args.version.partition("-rc.")
    numbers = base.split(".") + ([rc] if separator else [])
    if len(base.split(".")) != 3 or not all(
        n.isascii() and n.isdecimal() and (n == "0" or not n.startswith("0"))
        for n in numbers
    ):
        parser.error("expected X.Y.Z or X.Y.Z-rc.N")

    root = Path(__file__).resolve().parents[2]
    files = ["Cargo.toml", "python/pyproject.toml", README] + PACKAGE_VERSIONS + [file for file, _, _ in DEPENDENCIES]
    original = {file: (root / file).read_text() for file in files}
    try:
        updated = prepare_updates(original, args.version)
    except ValueError as error:
        parser.error(str(error))

    # Prepare all edits before writing any file.
    for file in original:
        if original[file] == updated[file]:
            continue
        if args.dry_run:
            print("".join(difflib.unified_diff(original[file].splitlines(True), updated[file].splitlines(True), file, file)), end="")
        else:
            (root / file).write_text(updated[file])
            print(f"Updated {file}")
    if args.dry_run:
        print("Would then run cargo update --workspace.")
    elif subprocess.run(LOCK_REFRESH, cwd=root).returncode:
        parser.exit(1, "cargo update --workspace failed; version edits remain.\n")


if __name__ == "__main__":
    main()
