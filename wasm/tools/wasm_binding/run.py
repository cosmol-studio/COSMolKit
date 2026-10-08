"""Build and validate the generated COSMolKit WASM JavaScript/TypeScript API."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import tempfile
import tomllib
from pathlib import Path

from presets import PRESETS, MODULE_GROUPS, TEST_GROUPS, active_features, npm_release, selected_names, source_modules
from features import sync_manifest


ROOT = Path(__file__).resolve().parents[3]
TEST_DIR = ROOT / "wasm" / "tests"
CONFIG = ROOT / "wasm" / "tools" / "wasm_binding" / "alef.toml"


def command(name: str, environment_name: str) -> str:
    configured = os.environ.get(environment_name)
    if configured:
        return configured
    resolved = shutil.which(name)
    if resolved:
        return resolved
    raise SystemExit(
        f"{name} is required; install it or set {environment_name} to its executable path"
    )


def run(*args: str, cwd: Path, env: dict[str, str] | None = None) -> None:
    subprocess.run(args, cwd=cwd, env=env, check=True)


def molecule_methods(sources) -> set[str]:
    return {
        name
        for source in sources
        for body in re.findall(
            r"^impl (?:crate::)?Molecule \{\n(.*?)^\}",
            source.read_text(encoding="utf-8"),
            re.MULTILINE | re.DOTALL,
        )
        for name in re.findall(r"pub fn (\w+)\s*\(", body)
    }


def check_generated_surface(source: Path, active: set[str] | None = None) -> None:
    # Check every current binding-facing Molecule method, including custom
    # modules. Historical names and future APIs are not this projection's ABI.
    expected = molecule_methods(source_modules(active or active_features("full")))
    if not expected:
        raise SystemExit("No binding-facing Molecule methods found")
    generated_methods = molecule_methods(source.parent.glob("*.rs"))
    missing = sorted(expected - generated_methods)
    if missing:
        raise SystemExit(
            "Alef omitted required ABI-safe Molecule methods: "
            + ", ".join(missing)
        )


def prepare_package(package: Path, library_name: str, metadata: dict, preset: str = "full") -> None:
    module = f"{library_name}.js"
    declaration = f"{library_name}.d.ts"
    version, tag = npm_release(metadata["version"], preset)
    (package / "package.json").write_text(
        json.dumps({
            "name": "@cosmol-studio/cosmolkit",
            "version": version,
            "cosmolkitPreset": preset,
            "cosmolkitVersion": metadata["version"],
            "publishConfig": {"access": "public", "tag": tag},
            "description": "WebAssembly bindings for COSMolKit",
            "type": "module",
            "main": module,
            "module": module,
            "types": declaration,
            "exports": {
                ".": {"types": f"./{declaration}", "default": f"./{module}"},
                "./*": "./*",
            },
            "files": ["*.js", "*.d.ts", "*.wasm", "snippets", "README.md", "LICENSE"],
            "license": metadata["license"],
            "repository": {"type": "git", "url": metadata["repository"]},
        }, indent=2) + "\n",
        encoding="utf-8",
    )
    shutil.copy2(ROOT / "wasm" / "README.md", package / "README.md")
    shutil.copy2(ROOT / "LICENSE", package / "LICENSE")


def check_typescript(workspace: Path, declaration: Path, active: set[str] | None = None) -> None:
    shim = workspace / "wasm-generated.d.ts"
    # Resolve the generated public declaration without depending on baseUrl.
    shim.write_text(
        f'export * from {json.dumps(str(declaration.with_suffix("")))};\n',
        encoding="utf-8",
    )
    type_config = workspace / "tsconfig.json"
    type_config.write_text(
        json.dumps({
            "compilerOptions": {
                "strict": True,
                "target": "ES2022",
                "module": "NodeNext",
                "moduleResolution": "NodeNext",
                "lib": ["ES2022", "DOM", "ESNext.Disposable"],
                "noEmit": True,
                "paths": {"cosmolkit-generated": [str(shim)]},
            },
            "files": list(map(str, selected_tests(".ts", active or active_features("full")))),
        }),
        encoding="utf-8",
    )
    # Runner images may ship a different tsc. Use the same pinned compiler
    # locally and in CI unless the caller explicitly requests another one.
    tsc = os.environ.get("TSC_BIN")
    if tsc:
        run(tsc, "--project", str(type_config), cwd=ROOT)
    else:
        run("npx", "--yes", "--package", "typescript@5.8.3", "tsc", "--project", str(type_config), cwd=ROOT)
    print("TypeScript declaration checks passed", flush=True)


def selected_tests(suffix: str, active: set[str]) -> list[Path]:
    known = {name for names in TEST_GROUPS.values() for name in names.split()}
    paths = sorted(TEST_DIR.glob("*" + suffix))
    unknown = {path.stem for path in paths} - known
    if unknown:
        raise SystemExit(f"Assign new WASM tests to a feature group: {sorted(unknown)}")
    selected = selected_names(TEST_GROUPS, active)
    return [path for path in paths if path.stem in selected]


def main() -> None:
    sync_manifest()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, help="Export the tested npm package to a new directory")
    parser.add_argument("--preset", choices=PRESETS, default="full", help="Fixed distribution preset (default: full)")
    args = parser.parse_args()
    active = active_features(args.preset)
    print(f"WASM preset: {args.preset}; features: {', '.join(sorted(active))}", flush=True)
    alef = command("alef", "ALEF_BIN")
    wasm_bindgen = command("wasm-bindgen", "WASM_BINDGEN_BIN")
    bun = os.environ.get("BUN_BIN") or shutil.which("bun")
    node = None if bun else command("node", "NODE_BIN")

    with tempfile.TemporaryDirectory(prefix="cosmolkit-wasm-") as temporary:
        workspace = Path(temporary)
        api_root = Path(os.environ.get("COSMOLKIT_API_ROOT", ROOT)).resolve()
        # Alef synchronizes versions. Only disposable input copies may be
        # writable through its workspace; source manifests must never be links.
        shutil.copytree(ROOT / "wasm" / "src", workspace / "wasm" / "src")
        manifest_text = (ROOT / "wasm" / "Cargo.toml").read_text(encoding="utf-8")
        manifest_text = manifest_text.replace(
            'path = "../crates/cosmolkit"',
            "path = " + json.dumps(str(api_root / "crates" / "cosmolkit")),
        )
        (workspace / "wasm" / "Cargo.toml").write_text(manifest_text, encoding="utf-8")
        with (api_root / "Cargo.toml").open("rb") as source_manifest:
            package_metadata = tomllib.load(source_manifest)["workspace"]["package"]
        (workspace / "Cargo.toml").write_text(
            '[workspace]\nmembers = ["wasm", "tmp/alef/wasm"]\nresolver = "2"\n'
            "[workspace.package]\n"
            + "\n".join(f"{key} = {json.dumps(package_metadata[key])}" for key in ("version", "edition", "license"))
            + "\n",
            encoding="utf-8",
        )
        config_text = CONFIG.read_text(encoding="utf-8")
        custom_modules = tomllib.loads(config_text)["crates"][0]["wasm"].get("custom_rust_modules", [])
        known = {name for names in MODULE_GROUPS.values() for name in names.split()}
        if unknown := set(custom_modules) - known:
            raise SystemExit(f"Assign new custom WASM modules to a feature group: {sorted(unknown)}")
        selected = selected_names(MODULE_GROUPS, active)
        custom_modules = [name for name in custom_modules if name in selected]
        config_text = re.sub(r"^features = .*", "features = " + json.dumps(sorted(active)), config_text, flags=re.MULTILINE)
        config_text = re.sub(r"^custom_rust_modules = .*", "custom_rust_modules = " + json.dumps(custom_modules), config_text, flags=re.MULTILINE)
        (workspace / "alef.toml").write_text(config_text, encoding="utf-8")

        private_bin = workspace / "private-bin"
        private_bin.mkdir()
        git_guard = private_bin / "git"
        git_guard.write_text("#!/bin/sh\nexit 127\n", encoding="utf-8")
        git_guard.chmod(0o755)
        build_env = os.environ.copy()
        build_env["PATH"] = str(private_bin) + os.pathsep + build_env.get("PATH", "")
        target = Path(build_env.get("CARGO_TARGET_DIR", workspace / "target")).resolve()
        build_env["CARGO_TARGET_DIR"] = str(target)

        generated = workspace / "tmp" / "alef" / "wasm"
        run(
            alef,
            "--config",
            str(workspace / "alef.toml"),
            "generate",
            "--crate",
            "cosmolkit-wasm",
            "--lang",
            "wasm",
            "--clean",
            cwd=workspace,
            env=build_env,
        )
        for module_name in custom_modules:
            shutil.copy2(
                ROOT / "wasm" / "src" / "js" / f"{module_name}.rs",
                generated / "src" / f"{module_name}.rs",
            )
        check_generated_surface(generated / "src" / "lib.rs", active)
        manifest = generated / "Cargo.toml"
        # cfg in hand-written projections uses the same resolved source features.
        # Disable dependency defaults so an omitted domain cannot leak in as full.
        manifest_text = manifest.read_text(encoding="utf-8")
        dependency = 'cosmolkit-wasm = { path = ' + json.dumps(str(workspace / "wasm")) + ', default-features = false, features = ' + json.dumps(PRESETS[args.preset]) + ' }'
        manifest_text, replacements = re.subn(r"^cosmolkit-wasm\s*=.*$", lambda _: dependency, manifest_text, flags=re.MULTILINE)
        if replacements != 1:
            raise SystemExit("Expected exactly one generated cosmolkit-wasm dependency")
        feature_table = "[features]\ndefault = " + json.dumps(sorted(active)) + "\n"
        feature_table += "\n".join(f"{feature} = []" for feature in tomllib.loads((ROOT / "wasm/Cargo.toml").read_text())["features"] if feature != "default") + "\n"
        manifest_text, tables = re.subn(r"^\[features\]\n.*?(?=^\[|\Z)", lambda _: feature_table + "\n", manifest_text, flags=re.MULTILINE | re.DOTALL)
        if tables == 0:
            manifest_text += "\n" + feature_table
        elif tables != 1:
            raise SystemExit("Expected at most one generated features table")
        manifest.write_text(manifest_text, encoding="utf-8")
        run(
            "cargo",
            "build",
            "--manifest-path",
            str(manifest),
            "--target",
            "wasm32-unknown-unknown",
            "--release",
            cwd=workspace,
            env=build_env,
        )

        with manifest.open("rb") as generated_manifest:
            manifest_data = tomllib.load(generated_manifest)
        library_name = manifest_data.get("lib", {}).get(
            "name", manifest_data["package"]["name"]
        ).replace("-", "_")
        wasm_binary = (
            target
            / "wasm32-unknown-unknown"
            / "release"
            / f"{library_name}.wasm"
        )
        package = generated / "pkg"
        package.mkdir(parents=True, exist_ok=True)
        run(
            wasm_bindgen,
            str(wasm_binary),
            "--target",
            "web",
            "--out-dir",
            str(package),
            cwd=workspace,
        )

        module = package / f"{library_name}.js"
        declaration = package / f"{library_name}.d.ts"
        background_binary = package / f"{library_name}_bg.wasm"
        for path in (module, declaration, background_binary):
            if not path.is_file():
                raise SystemExit(f"wasm-bindgen did not produce expected file: {path}")
        prepare_package(package, library_name, package_metadata, args.preset)

        runtime_env = build_env.copy()
        runtime_env.update(
            {
                "COSMOLKIT_WASM_MODULE": str(module),
                "COSMOLKIT_WASM_BINARY": str(background_binary),
                "COSMOLKIT_WASM_PRESET": args.preset,
            }
        )
        runtime_failure = None
        try:
            if bun:
                run(bun, "test", *map(str, selected_tests(".mjs", active)), cwd=ROOT, env=runtime_env)
            else:
                run(node, "--test", *map(str, selected_tests(".mjs", active)), cwd=ROOT, env=runtime_env)
        except subprocess.CalledProcessError as error:
            # Still check TypeScript when runtime tests fail; retain failure.
            runtime_failure = error

        check_typescript(workspace, declaration, active)
        if runtime_failure is not None:
            raise runtime_failure
        if args.out_dir is not None:
            shutil.copytree(package, args.out_dir.resolve())
            print(f"Tested npm package exported to {args.out_dir}", flush=True)
        print(f"WASM binary bytes ({args.preset}): {background_binary.stat().st_size}", flush=True)

    print("WASM JavaScript runtime and TypeScript declaration checks passed")


if __name__ == "__main__":
    main()
