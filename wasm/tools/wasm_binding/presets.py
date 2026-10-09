"""Fixed npm distributions and feature-scoped binding/test selection."""

from pathlib import Path
import re
import tomllib
from features import resolve

ROOT = Path(__file__).resolve().parents[3]
PRESETS = {
    "core": ["core"],
    "core-search": ["core", "search"],
    "core-fingerprints": ["core", "fingerprints"],
    "core-analysis": ["core", "search", "fingerprints", "descriptors"],
    "core-reaction": ["core", "reaction"],
    "core-depict": ["core", "depict"],
    "core-3d": ["core", "conformer"],
    "core-bio": ["core", "bio"],
    "core-inchi": ["core", "inchi"],
    "full": ["full"],
}

# Custom modules are language projections, not another chemistry API registry.
# Each module is assigned to its existing binding-facing Cargo feature.
MODULE_GROUPS = {
    "cap-transforms": "fragments",
    "cap-fingerprints+core": "avalon",
    "cap-stereo": "stereo_queries",
    "core": "query_values valence_operations element_boundary aromaticity_boundary smiles_boundary host_values smiles_parameters transform_parameters transform_errors property_values group_values io_errors sdf_reading property_strings xyz_mol2 xyz_mol2_errors mol_sdf_writing mol_write_errors sdf_datasets kekulize matrices matrix_errors bond_order radicals rings sanitize valence_errors hydrogens",
    "cap-search": "query_construction search search_errors mcs",
    "cap-alignment": "alignment_parameters alignment_operations",
    "cap-io": "alignment_values",
    "cap-batch": "batch_boundary batch_queries sdf_batches batch_images",
    "cap-fingerprints": "fingerprint_values fingerprint_errors fingerprint_source_errors atom_pair_parameters layered_pattern_values layered_pattern_errors morgan_parameters path_codes fingerprint_additional atom_pair_molecule morgan_molecule topological_torsion_parameters topological_torsion_molecule maccs layered_pattern_molecule topological",
    "cap-batch+cap-fingerprints": "batch_fingerprint_values batch_atom_pair batch_layered_pattern batch_morgan",
    "bio": "bio_residue bio_readers bio_read_errors bio_hierarchy bio_hierarchy_errors bio_parts bio_metadata bio_selection bio_writers bio_write_errors bio_protein bio_conversion",
    "conformer": "embed_parameters conformer_operations conformer_errors",
    "cap-depict": "depict_operations depict_errors drawing_errors image_errors",
    "descriptors": "descriptor_values descriptor_errors descriptors",
    "cap-forcefields": "forcefield_properties uff mmff",
    "cap-hashing": "hashing",
    "reaction": "reaction reaction_parameters reaction_errors",
    "inchi": "inchi",
}

TEST_GROUPS = {
    "cap-fingerprints+core": "avalon",
    "cap-transforms": "fragments",
    "cap-stereo": "stereo_queries",
    "core": "element_projection aromaticity_projection property_strings hydrogens xyz_mol2 kekulize matrices radicals sanitize preset_surface",
    "cap-search": "query_construction search mcs",
    "cap-search+cap-depict": "sdf_reading sdf_datasets mol_sdf_writing",
    "cap-depict": "depict_projection",
    "core+fingerprints": "valence",
    "core+descriptors": "rings",
    "cap-alignment": "alignment_parameters alignment_operations",
    "cap-batch": "batch_foundation",
    "cap-batch+cap-depict": "batch_transforms batch_images",
    "cap-batch+cap-conformer+cap-depict": "batch_queries",
    "cap-batch+cap-search+cap-depict": "sdf_batches",
    "fingerprints": "fingerprint_values atom_pair_molecule maccs path_codes torsion hashing fingerprint_presets",
    "cap-fingerprints+cap-search": "morgan_molecule layered_pattern topological",
    "cap-batch+cap-fingerprints+cap-search": "batch_morgan_query",
    "cap-batch+cap-fingerprints": "batch_atom_pair batch_layered_pattern batch_morgan",
    "bio": "bio_residue bio_readers bio_metadata bio_hierarchy bio_parts bio_selection bio_writers bio_protein bio_conversion",
    "conformer": "conformer_projection",
    "descriptors": "descriptors_projection",
    "cap-forcefields": "forcefield_properties uff mmff",
    "reaction": "reaction",
    "inchi": "inchi",
    # These original mixed-domain suites remain intact and run in full.
    "full": "binding_surface chemical_boundary_regressions",
}


def active_features(preset: str) -> set[str]:
    return resolve(PRESETS[preset])


def selected_names(groups: dict[str, str], active: set[str]) -> set[str]:
    return {name for group, names in groups.items() if set(group.split("+")) <= active for name in names.split()}


def npm_release(version: str, preset: str) -> tuple[str, str, str]:
    if not re.fullmatch(r"\d+\.\d+\.\d+(?:-rc\.\d+)?", version):
        raise ValueError(f"Unsupported release version: {version}")
    if preset not in PRESETS:
        raise ValueError(f"Unknown WASM preset: {preset}")
    name = "@cosmol-studio/cosmolkit"
    if preset != "full":
        name += f"-{preset}"
    return name, version, "rc" if "-" in version else "latest"


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
