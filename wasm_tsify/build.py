"""Build once, then let wasm-bindgen generate the Node package and declarations."""
from pathlib import Path
import argparse
import subprocess
import features

ROOT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--registry", action="store_true", help="include canonical registry inspection metadata")
parser.add_argument("--features", default="full", help="comma-separated Rust facade features, e.g. core,fingerprints")
parser.add_argument("--release", action="store_true")
parser.add_argument("--target", choices=["nodejs", "web", "bundler"], default="nodejs")
parser.add_argument("--out-dir", type=Path, default=ROOT / "pkg")
args = parser.parse_args()
features.sync()
selected = args.features.split(",") if args.features else []
if args.registry:
    selected.append("registry-check")
subprocess.run(["cargo", "build", "--manifest-path", str(ROOT / "Cargo.toml"),
                "--target", "wasm32-unknown-unknown", "--no-default-features",
                "--target-dir", str(ROOT / "target"),
                *(["--features", ",".join(selected)] if selected else []),
                *(["--release"] if args.release else [])], check=True)
profile = "release" if args.release else "debug"
subprocess.run(["wasm-bindgen", str(ROOT / f"target/wasm32-unknown-unknown/{profile}/cosmolkit_tsify_probe.wasm"),
                "--target", args.target, "--out-dir", str(args.out_dir)], check=True)
