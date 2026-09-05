"""Build and validate the generated COSMolKit WASM JavaScript/TypeScript API."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
import tomllib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
TEST_DIR = ROOT / "wasm" / "tests"
CONFIG = ROOT / "wasm" / "tools" / "wasm_binding" / "alef.toml"

# These are the stable, ABI-safe operations promised by the binding-facing
# projection.  The Rust facade remains a complete re-export; this list covers
# the methods that must survive the Alef -> wasm-bindgen projection.  Optional
# source results have explicit `*_or_error` counterparts and are checked
# through those counterparts instead of being silently dropped by a generator.
EXPECTED_METHODS = {
    "new",
    "from_smiles",
    "from_smiles_with_sanitize",
    "from_mol_block",
    "from_xyz_block",
    "from_pdb_block",
    "from_mmcif_block",
    "from_protein_sequence",
    "from_nucleic_sequence",
    "from_mol2_or_error",
    "from_inchi_or_error",
    "from_binary",
    "to_smiles",
    "to_smiles_with_isomeric",
    "to_mol_block_v2000",
    "to_mol_block_v3000",
    "to_sdf_record_v2000",
    "to_sdf_record_v3000",
    "to_smarts",
    "to_cx_smarts",
    "to_inchi",
    "to_inchi_key",
    "num_atoms",
    "num_bonds",
    "name_or_empty",
    "with_name",
    "with_property",
    "with_sdf_data_field",
    "property_or_empty",
    "property_keys",
    "sdf_data_field_names",
    "sdf_data_field_or_empty",
    "source_coordinate_dimension_or_empty",
    "atomic_numbers",
    "coordinates_2d",
    "num_conformers_3d",
    "coordinates_3d",
    "to_svg",
    "to_png",
    "molecular_weight",
    "exact_molecular_weight",
    "crippen_log_p",
    "crippen_molar_refractivity",
    "tpsa",
    "h_bond_acceptors",
    "h_bond_donors",
    "formula",
    "fraction_csp3",
    "num_heavy_atoms",
    "num_rings",
    "num_aromatic_rings",
    "num_rotatable_bonds",
    "qed",
    "hall_kier_alpha",
    "kappa_1",
    "kappa_2",
    "kappa_3",
    "chi",
    "chi_0",
    "chi_1",
    "phi",
    "num_spiro_atoms",
    "num_bridgehead_atoms",
    "labute_asa",
    "num_atom_stereo_centers",
    "num_unspecified_atom_stereo_centers",
    "pattern_fingerprint",
    "morgan_fingerprint",
    "atom_pair_fingerprint",
    "layered_fingerprint",
    "topological_fingerprint",
    "maccs_fingerprint",
    "avalon_fingerprint",
    "topological_torsion_fingerprint",
    "with_hydrogens",
    "without_hydrogens",
    "with_kekulized_bonds",
    "sanitize",
    "with_assigned_valence",
    "with_assigned_rings",
    "with_assigned_ring_families",
    "with_assigned_aromaticity",
    "with_assigned_radicals",
    "with_2d_coordinates",
    "with_3d_conformer",
    "perceive_stereochemistry",
    "stereoisomer_count",
    "stereoisomer_count_with_options",
    "hash",
    "to_binary",
    "to_pdb_block",
    "dg_bounds_matrix",
    "fragments",
    "largest_fragment",
    "murcko_scaffold",
    "net_scaffold",
    "has_substruct_match",
}


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


def check_generated_surface(source: Path) -> None:
    generated_methods = {
        line.strip().split("(", 1)[0].removeprefix("pub fn ")
        for line in source.read_text(encoding="utf-8").splitlines()
        if line.strip().startswith("pub fn ")
    }
    missing = sorted(EXPECTED_METHODS - generated_methods)
    if missing:
        raise SystemExit(
            "Alef omitted required ABI-safe Molecule methods: "
            + ", ".join(missing)
        )


def main() -> None:
    alef = command("alef", "ALEF_BIN")
    wasm_bindgen = command("wasm-bindgen", "WASM_BINDGEN_BIN")
    bun = os.environ.get("BUN_BIN") or shutil.which("bun")
    node = None if bun else command("node", "NODE_BIN")

    with tempfile.TemporaryDirectory(prefix="cosmolkit-wasm-") as temporary:
        workspace = Path(temporary)
        (workspace / "crates").symlink_to(ROOT / "crates", target_is_directory=True)
        (workspace / "wasm").symlink_to(ROOT / "wasm", target_is_directory=True)
        (workspace / "Cargo.toml").symlink_to(ROOT / "Cargo.toml")
        (workspace / "alef.toml").write_text(CONFIG.read_text(encoding="utf-8"), encoding="utf-8")

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
        )
        check_generated_surface(generated / "src" / "lib.rs")
        manifest = generated / "Cargo.toml"
        run(
            "cargo",
            "build",
            "--manifest-path",
            str(manifest),
            "--target",
            "wasm32-unknown-unknown",
            "--release",
            cwd=workspace,
        )

        with manifest.open("rb") as generated_manifest:
            manifest_data = tomllib.load(generated_manifest)
        library_name = manifest_data.get("lib", {}).get(
            "name", manifest_data["package"]["name"]
        ).replace("-", "_")
        wasm_binary = (
            generated
            / "target"
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

        runtime_env = os.environ.copy()
        runtime_env.update(
            {
                "COSMOLKIT_WASM_MODULE": str(module),
                "COSMOLKIT_WASM_BINARY": str(background_binary),
            }
        )
        if bun:
            run(bun, "test", str(TEST_DIR / "binding_surface.mjs"), cwd=ROOT, env=runtime_env)
        else:
            run(node, "--test", str(TEST_DIR / "binding_surface.mjs"), cwd=ROOT, env=runtime_env)

        shim = workspace / "wasm-generated.d.ts"
        # Import the generated module without its `.d.ts` suffix so TypeScript
        # resolves the declaration as the module's public type surface.
        shim.write_text(
            f'export * from {json.dumps(str(declaration.with_suffix("")))};\n',
            encoding="utf-8",
        )
        type_config = workspace / "tsconfig.json"
        type_config.write_text(
            json.dumps(
                {
                    "compilerOptions": {
                        "strict": True,
                        "target": "ES2022",
                        "module": "NodeNext",
                        "moduleResolution": "NodeNext",
                        "lib": ["ES2022", "DOM", "ESNext.Disposable"],
                        "noEmit": True,
                        "baseUrl": str(TEST_DIR),
                        "paths": {"cosmolkit-generated": [str(shim)]},
                    },
                    "files": [str(TEST_DIR / "binding_surface.ts")],
                }
            ),
            encoding="utf-8",
        )
        tsc = shutil.which("tsc")
        if tsc:
            run(tsc, "--project", str(type_config), cwd=ROOT)
        else:
            run("npx", "--yes", "--package", "typescript@5.8.3", "tsc", "--project", str(type_config), cwd=ROOT)

    print("WASM JavaScript runtime and TypeScript declaration checks passed")


if __name__ == "__main__":
    main()
