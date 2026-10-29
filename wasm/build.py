"""Build the WASM package and declarations; optionally check binding protocols."""
from pathlib import Path
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import tomllib

import features
from presets import PRESETS, npm_release

ROOT = Path(__file__).resolve().parent
LIBRARY = "cosmolkit_wasm"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preset", choices=PRESETS, default="full")
    parser.add_argument("--features", help="explicit comma-separated Rust facade features")
    parser.add_argument("--check", action="store_true", help="check registry and fixed framework samples")
    parser.add_argument("--registry", action="store_true", help="include registry inspection metadata")
    parser.add_argument("--release", action="store_true")
    parser.add_argument("--target", choices=["nodejs", "web", "bundler"], default="nodejs")
    parser.add_argument("--out-dir", type=Path, default=ROOT / "pkg")
    args = parser.parse_args()
    if (args.check or args.registry) and args.target != "nodejs":
        parser.error("inspection builds require --target nodejs")
    features.sync()
    selected = args.features.split(",") if args.features is not None else PRESETS[args.preset].copy()
    full = "full" in selected
    if args.check or args.registry:
        selected.append("registry-check")
    if args.check and full:
        selected.append("protocol-tests")
    target = Path(os.environ.get("CARGO_TARGET_DIR", ROOT / "target")).resolve()
    out = args.out_dir.resolve()
    subprocess.run(["cargo", "build", "--locked", "--manifest-path", str(ROOT / "Cargo.toml"),
                    "--target", "wasm32-unknown-unknown", "--no-default-features",
                    "--target-dir", str(target),
                    *(["--features", ",".join(selected)] if selected else []),
                    *(["--release"] if args.release else [])], check=True)
    profile = "release" if args.release else "debug"
    # Package only this generation's files, even when an output directory is reused.
    with tempfile.TemporaryDirectory(prefix="cosmolkit-wasm-") as directory:
        staging = Path(directory)
        subprocess.run(["wasm-bindgen", str(target / f"wasm32-unknown-unknown/{profile}/{LIBRARY}.wasm"),
                        "--target", args.target, "--out-dir", str(staging)], check=True)
        generated = sorted(str(path.relative_to(staging)) for path in staging.rglob("*") if path.is_file())
        shutil.copytree(staging, out, dirs_exist_ok=True)
    version = tomllib.loads((ROOT / "Cargo.toml").read_text())["package"]["version"]
    name, version, tag = npm_release(version, args.preset)
    package = {
        "name": name, "version": version,
        "description": "Rust-native cheminformatics for JavaScript and TypeScript",
        "license": "MIT", "repository": "https://github.com/cosmol-studio/COSMolKit",
        "type": "commonjs" if args.target == "nodejs" else "module",
        "main": f"{LIBRARY}.js", "types": f"{LIBRARY}.d.ts",
        "files": generated + ["snippets/package.json", "README.md", "LICENSE"],
        "publishConfig": {"access": "public", "tag": tag},
    }
    if args.check or args.registry or args.features is not None:
        package["private"] = True
    (out / "package.json").write_text(json.dumps(package, indent=2) + "\n")
    # wasm-bindgen's imported JS snippets remain ES modules for every target.
    if (out / "snippets").is_dir():
        (out / "snippets/package.json").write_text('{"type":"module"}\n')
    shutil.copyfile(ROOT / "README.md", out / "README.md")
    shutil.copyfile(ROOT.parent / "LICENSE", out / "LICENSE")
    if args.check:
        subprocess.run([sys.executable, str(ROOT / "registry_check.py"),
                        "--package", str(out / f"{LIBRARY}.js")], check=True)
        if full:
            subprocess.run(["node", str(ROOT / "tests/runtime.cjs")], check=True,
                           env={**os.environ, "COSMOLKIT_WASM_PACKAGE": str(out / f"{LIBRARY}.js")})
        config = out / "tsconfig.check.json"
        options = {
            "extends": str(ROOT / "tests/tsconfig.json"),
            "compilerOptions": {"paths": {"cosmolkit": [str(out / f"{LIBRARY}.d.ts")]}}
        }
        if not full:
            options["files"] = [str(out / f"{LIBRARY}.d.ts")]
        config.write_text(json.dumps(options))
        subprocess.run(["npx", "--yes", "--package", "typescript@5.8.3",
                        "tsc", "-p", str(config)], check=True)


if __name__ == "__main__":
    main()
