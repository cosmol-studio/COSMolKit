"""Build the optional search library target and generate wasm-bindgen assets."""

import os
from pathlib import Path
import shutil
import subprocess
import sys
import tomllib

from extract_sphinx_metadata import rust_literal


def binding_version(root):
    manifest = tomllib.loads((root / "Cargo.toml").read_text(encoding="utf-8"))
    return manifest["dependencies"]["wasm-bindgen"]["version"].removeprefix("=")


def binding_tool(root):
    version = binding_version(root)
    suffix = ".exe" if os.name == "nt" else ""
    candidates = [os.environ.get("COSMOLKIT_WASM_BINDGEN"),
                  str(root / "target/search-tools/bin" / ("wasm-bindgen" + suffix)),
                  shutil.which("wasm-bindgen")]
    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            result = subprocess.run([candidate, "--version"], capture_output=True, text=True, check=True)
            if result.stdout.strip() == "wasm-bindgen " + version:
                return candidate
    raise RuntimeError(f"wasm-bindgen-cli {version} is required; install with cargo binstall wasm-bindgen-cli --version {version} --root {root / 'target/search-tools'} --no-confirm --disable-strategies compile")


def build_search_bundle(out: Path):
    root = Path(__file__).resolve().parents[1]
    tool = binding_tool(root)
    # A separate target directory avoids contending with the enclosing Cargo
    # build. --lib + search-engine takes build.rs's index-only branch.
    target = root / "target/search-engine"
    environment = {key: value for key, value in os.environ.items()
                   if not key.startswith(("CARGO_FEATURE_", "CARGO_CFG_"))
                   and key not in ("OUT_DIR", "CARGO_ENCODED_RUSTFLAGS", "RUSTFLAGS")}
    command = [os.environ.get("CARGO", "cargo"), "build", "--manifest-path", str(root / "Cargo.toml"),
               "--lib", "--no-default-features", "--features", "search-engine",
               "--target", "wasm32-unknown-unknown", "--target-dir", str(target), "--release", "--locked",
               "--config", 'profile.release.opt-level="z"', "--config", "profile.release.lto=true",
               "--config", "profile.release.codegen-units=1"]
    subprocess.run(command, cwd=root, env=environment, check=True, stdout=subprocess.PIPE)
    subprocess.run([tool, str(target / "wasm32-unknown-unknown/release/cosmolkit_docs_search.wasm"),
                    "--target", "web", "--out-dir", str(out), "--out-name", "docs_search", "--no-typescript"],
                   check=True, stdout=subprocess.PIPE)
    # asset!("/...") resolves from the crate root, not the filesystem root.
    # Filesystem-absolute Unix paths are otherwise prefixed with the crate twice.
    # build.rs can invoke this script through a Windows extended-length path.
    asset_root = Path(str(root).removeprefix("\\\\?\\"))
    asset_out = Path(str(out.resolve()).removeprefix("\\\\?\\"))
    asset_dir = asset_out.relative_to(asset_root).as_posix()
    bindings = rust_literal(f"/{asset_dir}/docs_search.js")
    wasm = rust_literal(f"/{asset_dir}/docs_search_bg.wasm")
    (out / "search_asset.rs").write_text(
        f"const SEARCH_BINDINGS: Asset = asset!({bindings});\nconst SEARCH_WASM: Asset = asset!({wasm});\n",
        encoding="utf-8")


if __name__ == "__main__":
    if sys.argv[1:] != ["--install-tools"]:
        raise SystemExit("usage: build_search_bundle.py --install-tools")
    root = Path(__file__).resolve().parents[1]
    try:
        tool = binding_tool(root)
    except RuntimeError:
        subprocess.run(["cargo", "binstall", "wasm-bindgen-cli", "--version", binding_version(root),
                        "--root", str(root / "target/search-tools"), "--no-confirm", "--disable-strategies", "compile"], check=True)
        tool = binding_tool(root)
    print(f"Search binding generator: {tool}")
