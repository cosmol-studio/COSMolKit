"""Build once, then let wasm-bindgen generate the Node package and declarations."""
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parent
subprocess.run(["cargo", "build", "--manifest-path", str(ROOT / "Cargo.toml"),
                "--target", "wasm32-unknown-unknown"], check=True)
subprocess.run(["wasm-bindgen", str(ROOT / "target/wasm32-unknown-unknown/debug/cosmolkit_tsify_probe.wasm"),
                "--target", "nodejs", "--out-dir", str(ROOT / "pkg")], check=True)
