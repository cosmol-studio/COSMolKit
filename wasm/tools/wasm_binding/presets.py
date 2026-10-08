"""Fixed npm distributions and feature-scoped binding/test selection."""

from pathlib import Path
import re
import tomllib

ROOT = Path(__file__).resolve().parents[3]
PRESETS = {
    "core": ["core"],
    "core-fingerprints": ["core", "fingerprints"],
    "core-analysis": ["core", "fingerprints", "descriptors"],
    "core-reaction": ["core", "reaction"],
    "core-3d": ["core", "conformer"],
    "core-bio": ["core", "bio"],
    "core-inchi": ["core", "inchi"],
    "full": ["full"],
}

# Custom modules are language projections, not another chemistry API registry.
# Each module is assigned to its existing binding-facing Cargo feature.
MODULE_GROUPS = {
    "core": "valence_operations element_boundary aromaticity_boundary smiles_boundary host_values smiles_parameters transform_parameters transform_errors drawing_errors image_errors query_construction property_values group_values io_errors sdf_reading property_strings xyz_mol2 xyz_mol2_errors mol_sdf_writing mol_write_errors sdf_datasets kekulize matrices matrix_errors bond_order radicals rings sanitize valence_errors search search_errors hydrogens alignment_values",
    "alignment": "alignment_parameters alignment_operations",
    "batch": "batch_boundary batch_queries batch_images sdf_batches",
    "fingerprints": "fingerprint_values fingerprint_errors fingerprint_source_errors atom_pair_parameters layered_pattern_values layered_pattern_errors morgan_parameters path_codes fingerprint_additional atom_pair_molecule morgan_molecule topological_torsion_parameters topological_torsion_molecule maccs layered_pattern_molecule topological",
    "batch+fingerprints": "batch_fingerprint_values batch_atom_pair batch_layered_pattern batch_morgan",
    "bio": "bio_residue bio_readers bio_read_errors bio_hierarchy bio_hierarchy_errors bio_parts bio_metadata bio_selection bio_writers bio_write_errors bio_protein bio_conversion",
    "conformer": "embed_parameters conformer_operations conformer_errors",
    "depict": "depict_operations depict_errors",
    "descriptors": "descriptor_values descriptor_errors descriptors",
    "forcefields": "forcefield_properties uff mmff",
    "hashing": "hashing",
    "serialization": "serialization pickle_errors",
    "reaction": "reaction reaction_parameters reaction_errors",
    "inchi": "inchi",
}

TEST_GROUPS = {
    "core": "element_projection aromaticity_projection property_strings query_construction depict_projection hydrogens sdf_reading xyz_mol2 mol_sdf_writing sdf_datasets kekulize matrices radicals sanitize search preset_surface",
    "core+fingerprints": "valence",
    "core+descriptors": "rings",
    "alignment": "alignment_parameters alignment_operations",
    "batch": "batch_foundation batch_queries batch_transforms batch_images sdf_batches",
    "fingerprints": "fingerprint_values atom_pair_molecule morgan_molecule layered_pattern maccs path_codes topological torsion hashing",
    "batch+fingerprints": "batch_atom_pair batch_layered_pattern batch_morgan",
    "bio": "bio_residue bio_readers bio_metadata bio_hierarchy bio_parts bio_selection bio_writers bio_protein bio_conversion",
    "conformer": "conformer_projection",
    "descriptors": "descriptors_projection",
    "forcefields": "forcefield_properties uff mmff",
    "serialization": "serialization",
    "reaction": "reaction",
    "inchi": "inchi",
    # These original mixed-domain suites remain intact and run in full.
    "full": "binding_surface chemical_boundary_regressions",
}


def active_features(preset: str) -> set[str]:
    with (ROOT / "wasm/Cargo.toml").open("rb") as source:
        features = tomllib.load(source)["features"]
    pending = list(PRESETS[preset])
    active = set()
    while pending:
        feature = pending.pop()
        if feature not in active:
            active.add(feature)
            pending.extend(value for value in features[feature] if "/" not in value and not value.startswith("dep:"))
    return active


def selected_names(groups: dict[str, str], active: set[str]) -> set[str]:
    return {name for group, names in groups.items() if set(group.split("+")) <= active for name in names.split()}


def npm_release(version: str, preset: str) -> tuple[str, str]:
    if not re.fullmatch(r"\d+\.\d+\.\d+(?:-rc\.\d+)?", version):
        raise ValueError(f"Unsupported release version: {version}")
    if preset not in PRESETS:
        raise ValueError(f"Unknown WASM preset: {preset}")
    if preset == "full":
        return version, "rc" if "-" in version else "latest"
    separator = "." if "-" in version else "-"
    return f"{version}{separator}{preset}.0", preset


def source_modules(active: set[str]) -> list[Path]:
    """Use the source module declarations for the existing completeness gate."""
    root = ROOT / "wasm/src/lib.rs"
    modules = [root]
    for attrs, name in re.findall(r"((?:#\[[^\n]+\]\n)*)mod (\w+);", root.read_text()):
        if "test" in attrs:
            continue
        required = set(re.findall(r'feature = "([^"]+)"', attrs))
        if required <= active:
            path = re.search(r'#\[path = "([^"]+)"\]', attrs)
            modules.append(root.parent / (path[1] if path else name + ".rs"))
    return modules
