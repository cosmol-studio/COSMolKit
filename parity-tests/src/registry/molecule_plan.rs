//! Frozen molecular matrices and future task scopes.
//!
//! One task is one observable chemistry operation, not a historical audit phase.
//! Corpus size, binding language and value/in-place projection are NOT chemistry
//! parameter axes. See ../../TASKS.md for preparation and comparison contracts.
//! The parent registry explicitly selects executable rows from these matrices.

/// Stable behavior categories; readiness never changes category.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Category {
    Notation,
    MolecularIo,
    Chemistry,
    Stereo,
    Descriptors,
    FingerprintGeneration,
    FingerprintValues,
    Query,
    Depiction,
    Conformers,
    ForceFields,
    Alignment,
    Tautomers,
    Identifiers,
    Bio,
    Composition,
    Native,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Reference {
    Rdkit,
    Gemmi,
    RustScalar,
    Roundtrip,
    NoExternalReference,
}

/// A future task's scope, NOT an executable or frozen parameter matrix.
/// Exact values, combinations, errors, reference pin and public boundary must
/// be frozen before promotion to typed Profile variants. Unknown counts are
/// None, never zero or an invented Cartesian product of every option.
#[derive(Debug, Clone, Copy)]
pub struct FutureTask {
    pub id: &'static str,
    pub category: Category,
    pub reference: Reference,
    pub parameter_axes: &'static str,
    pub comparison: &'static str,
}

pub const FUTURE_TASKS: &[FutureTask] = &[
    FutureTask {
        id: "smiles_write",
        category: Category::Notation,
        reference: Reference::Rdkit,
        parameter_axes: "canonical, isomeric, kekule, explicit bonds/H, root and random traversal",
        comparison: "Text and atom/bond output order",
    },
    FutureTask {
        id: "cxsmiles_read",
        category: Category::Notation,
        reference: Reference::Rdkit,
        parameter_axes: "CX fields, strictness, names, coordinate records",
        comparison: "Typed topology, query predicates, properties and coordinates",
    },
    FutureTask {
        id: "cxsmiles_write",
        category: Category::Notation,
        reference: Reference::Rdkit,
        parameter_axes: "CX field mask, coordinate selector, canonical/isomeric options",
        comparison: "Text and output mappings; approved coordinate-selection difference explicit",
    },
    FutureTask {
        id: "smarts_read",
        category: Category::Notation,
        reference: Reference::Rdkit,
        parameter_axes: "merge H, replacements, CX/name policy",
        comparison: "Ordered query graph and predicate trees",
    },
    FutureTask {
        id: "smarts_write",
        category: Category::Notation,
        reference: Reference::Rdkit,
        parameter_axes: "root, isomeric and CX fields",
        comparison: "Exact text and query semantics",
    },
    FutureTask {
        id: "molblock_read",
        category: Category::MolecularIo,
        reference: Reference::Rdkit,
        parameter_axes: "V2000/V3000, sanitize, remove H, strict parsing, coordinate mode",
        comparison: "Concrete/query result, graph, stereo, SGroups, coordinates and errors",
    },
    FutureTask {
        id: "molblock_write",
        category: Category::MolecularIo,
        reference: Reference::Rdkit,
        parameter_axes: "V2000/V3000, stereo, kekule, coordinate selector",
        comparison: "Text and preserved graph fields",
    },
    FutureTask {
        id: "sdf_read",
        category: Category::MolecularIo,
        reference: Reference::Rdkit,
        parameter_axes: "record boundaries, invalid records, duplicate data fields, reader options",
        comparison: "Ordered records, errors, data fields and graph category",
    },
    FutureTask {
        id: "sdf_write",
        category: Category::MolecularIo,
        reference: Reference::Rdkit,
        parameter_axes: "record/data-field order, molblock options",
        comparison: "Serialized records and fields",
    },
    FutureTask {
        id: "mol2_read",
        category: Category::MolecularIo,
        reference: Reference::Rdkit,
        parameter_axes: "sanitize, remove H, substructure cleanup",
        comparison: "Topology, atom types, charge and outcome",
    },
    FutureTask {
        id: "binary_roundtrip",
        category: Category::MolecularIo,
        reference: Reference::Roundtrip,
        parameter_axes: "serialization version, property/coordinate selection",
        comparison: "Before/after CK public state; no RDKit byte-compatibility claim",
    },
    FutureTask {
        id: "sanitize_stages",
        category: Category::Chemistry,
        reference: Reference::Rdkit,
        parameter_axes: "named single stages, prerequisite sequences and approved stage combinations",
        comparison: "Stage result and structured failure; separate from initial ALL profile",
    },
    FutureTask {
        id: "aromaticity",
        category: Category::Chemistry,
        reference: Reference::Rdkit,
        parameter_axes: "each supported model, prepared kekule input",
        comparison: "Atom/bond aromatic flags and bond orders",
    },
    FutureTask {
        id: "radicals",
        category: Category::Chemistry,
        reference: Reference::Rdkit,
        parameter_axes: "raw versus prepared input",
        comparison: "Per-atom radical counts and errors",
    },
    FutureTask {
        id: "conjugation",
        category: Category::Chemistry,
        reference: Reference::Rdkit,
        parameter_axes: "bond types and prepared chemistry state",
        comparison: "Per-bond conjugation",
    },
    FutureTask {
        id: "hybridization",
        category: Category::Chemistry,
        reference: Reference::Rdkit,
        parameter_axes: "prepared valence/aromaticity state",
        comparison: "Per-atom hybridization",
    },
    FutureTask {
        id: "cleanup",
        category: Category::Chemistry,
        reference: Reference::Rdkit,
        parameter_axes: "ordinary, organometallic and atropisomer cleanup as separate profiles",
        comparison: "Topology changes and outcome",
    },
    FutureTask {
        id: "rings_fast",
        category: Category::Chemistry,
        reference: Reference::Rdkit,
        parameter_axes: "uninitialized and existing ring state",
        comparison: "Ring rows, membership and find type",
    },
    FutureTask {
        id: "rings_sssr",
        category: Category::Chemistry,
        reference: Reference::Rdkit,
        parameter_axes: "dative/hydrogen bond inclusion",
        comparison: "Ordered rings and membership",
    },
    FutureTask {
        id: "rings_symmetrized",
        category: Category::Chemistry,
        reference: Reference::Rdkit,
        parameter_axes: "dative/hydrogen bond inclusion",
        comparison: "Symmetrized rings and membership",
    },
    FutureTask {
        id: "ring_families",
        category: Category::Chemistry,
        reference: Reference::Rdkit,
        parameter_axes: "bond inclusion, fused/spiro/bridged inputs",
        comparison: "Family rows, counts and relevant cycles",
    },
    FutureTask {
        id: "fragments",
        category: Category::Chemistry,
        reference: Reference::Rdkit,
        parameter_axes: "index groups versus molecule outputs, sanitize fragments",
        comparison: "Components, source mappings and output graphs",
    },
    FutureTask {
        id: "distance_matrix_3d",
        category: Category::Chemistry,
        reference: Reference::Rdkit,
        parameter_axes: "supplied 3D coordinates, conformer selection and atom weights",
        comparison: "Dimensions and every matrix entry",
    },
    FutureTask {
        id: "hydrogen_options",
        category: Category::Chemistry,
        reference: Reference::Rdkit,
        parameter_axes: "each remaining AddHs/RemoveHs option and source-required interactions",
        comparison: "Topology, isotopes, mapping, coordinates and stereo",
    },
    FutureTask {
        id: "structure_stereo",
        category: Category::Stereo,
        reference: Reference::Rdkit,
        parameter_axes: "2D/3D selection, replace existing tags",
        comparison: "Atom tags and bond stereo",
    },
    FutureTask {
        id: "stereo_enumeration",
        category: Category::Stereo,
        reference: Reference::Rdkit,
        parameter_axes: "unassigned-only, unique, enhanced groups, embedding, limits and fixed seed",
        comparison: "Enumerated structures, ordering and completion state",
    },
    FutureTask {
        id: "connectivity",
        category: Category::Descriptors,
        reference: Reference::Rdkit,
        parameter_axes: "Chi families/orders, Hall-Kier, Kappa, Phi",
        comparison: "Scalar descriptor values",
    },
    FutureTask {
        id: "lipinski_counts",
        category: Category::Descriptors,
        reference: Reference::Rdkit,
        parameter_axes: "HBD/HBA, rotatable-bond definitions, atom/ring classifications",
        comparison: "Exact integer counts",
    },
    FutureTask {
        id: "mqn",
        category: Category::Descriptors,
        reference: Reference::Rdkit,
        parameter_axes: "complete component vector",
        comparison: "Every vector element",
    },
    FutureTask {
        id: "crippen",
        category: Category::Descriptors,
        reference: Reference::Rdkit,
        parameter_axes: "hydrogen handling, per-atom contributions",
        comparison: "LogP/MR totals and contribution arrays",
    },
    FutureTask {
        id: "surface_area",
        category: Category::Descriptors,
        reference: Reference::Rdkit,
        parameter_axes: "Labute hydrogen policy; VSA family and default/custom bins",
        comparison: "Scalar area, bins and contributions",
    },
    FutureTask {
        id: "morgan",
        category: Category::FingerprintGeneration,
        reference: Reference::Rdkit,
        parameter_axes: "radius, bit/count/sparse output, chirality, features, invariants, roots, provenance",
        comparison: "Vector and complete additional output",
    },
    FutureTask {
        id: "maccs",
        category: Category::FingerprintGeneration,
        reference: Reference::Rdkit,
        parameter_axes: "standard key definition",
        comparison: "Exact bit vector",
    },
    FutureTask {
        id: "rdk_fingerprint",
        category: Category::FingerprintGeneration,
        reference: Reference::Rdkit,
        parameter_axes: "path limits, branched paths, H, bond order, size, roots, provenance",
        comparison: "Vector and atom/bond provenance",
    },
    FutureTask {
        id: "avalon",
        category: Category::FingerprintGeneration,
        reference: Reference::Rdkit,
        parameter_axes: "size, query mode and supported flag profiles",
        comparison: "Exact bit vector",
    },
    FutureTask {
        id: "atom_pair",
        category: Category::FingerprintGeneration,
        reference: Reference::Rdkit,
        parameter_axes: "distance range, 2D/3D, chirality, roots/ignored atoms, bit/count/sparse output",
        comparison: "Vector and additional output",
    },
    FutureTask {
        id: "topological_torsion",
        category: Category::FingerprintGeneration,
        reference: Reference::Rdkit,
        parameter_axes: "torsion length, chirality, roots, count simulation, output representation",
        comparison: "Vector and additional output",
    },
    FutureTask {
        id: "pattern",
        category: Category::FingerprintGeneration,
        reference: Reference::Rdkit,
        parameter_axes: "size, tautomer mode, atom counts and set-only mask",
        comparison: "Vector and updated counts",
    },
    FutureTask {
        id: "layered",
        category: Category::FingerprintGeneration,
        reference: Reference::Rdkit,
        parameter_axes: "layer mask, paths, roots, count seeds, set-only mask",
        comparison: "Vector and updated counts",
    },
    FutureTask {
        id: "fingerprint_arithmetic",
        category: Category::FingerprintValues,
        reference: Reference::Rdkit,
        parameter_axes: "representation/index width, binary/scalar operation, signed/zero counts",
        comparison: "Length, exact entries and error",
    },
    FutureTask {
        id: "fingerprint_similarity",
        category: Category::FingerprintValues,
        reference: Reference::Rdkit,
        parameter_axes: "named metric, count/bit semantics, empty vectors, metric parameters",
        comparison: "Scalar result and error",
    },
    FutureTask {
        id: "fingerprint_serialization",
        category: Category::FingerprintValues,
        reference: Reference::Rdkit,
        parameter_axes: "supported text/binary form and index width",
        comparison: "Bytes/text and reconstructed entries",
    },
    FutureTask {
        id: "substructure_match",
        category: Category::Query,
        reference: Reference::Rdkit,
        parameter_axes: "chirality, query-query, recursion, uniqueness, max matches, properties, callbacks",
        comparison: "Boolean, ordered mappings and errors",
    },
    FutureTask {
        id: "compiled_query",
        category: Category::Query,
        reference: Reference::Rdkit,
        parameter_axes: "same matching profiles as ordinary matching",
        comparison: "Same result as reference and ordinary CK path",
    },
    FutureTask {
        id: "generic_groups",
        category: Category::Query,
        reference: Reference::Rdkit,
        parameter_axes: "each supported generic label and generic-matcher option",
        comparison: "Match results",
    },
    FutureTask {
        id: "mcs",
        category: Category::Query,
        reference: Reference::Rdkit,
        parameter_axes: "atom/bond comparators, ring/stereo constraints, objective, seed, threshold, limits, ties",
        comparison: "Counts, query, SMARTS, degenerate results and completion",
    },
    FutureTask {
        id: "depiction_options",
        category: Category::Depiction,
        reference: Reference::Rdkit,
        parameter_axes: "orientation, templates, constrained maps, sampling/seed, mimic distances",
        comparison: "Coordinates and source topology",
    },
    FutureTask {
        id: "depiction_normalize",
        category: Category::Depiction,
        reference: Reference::Rdkit,
        parameter_axes: "normalization, straightening, conformer selection and scaling",
        comparison: "Coordinates and returned transform/scale",
    },
    FutureTask {
        id: "svg",
        category: Category::Depiction,
        reference: Reference::Rdkit,
        parameter_axes: "labels, highlights, stereo, SGroups, annotations, dimensions",
        comparison: "Declared SVG structure/text boundary; approved CK branding explicit",
    },
    FutureTask {
        id: "png",
        category: Category::Depiction,
        reference: Reference::Rdkit,
        parameter_axes: "same drawing profiles, image size and renderer",
        comparison: "Declared decoded pixel/rendering boundary, not unexplained byte equality",
    },
    FutureTask {
        id: "bounds",
        category: Category::Conformers,
        reference: Reference::Rdkit,
        parameter_axes: "bond policy, smoothing, macrocycle and fragment settings",
        comparison: "Every bound entry and outcome",
    },
    FutureTask {
        id: "embedding",
        category: Category::Conformers,
        reference: Reference::Rdkit,
        parameter_axes: "ETKDG version, count, fixed seed, chirality, random coordinates, constraints",
        comparison: "Outcome, conformer count, coordinates and diagnostics",
    },
    FutureTask {
        id: "conformer_pruning",
        category: Category::Conformers,
        reference: Reference::Rdkit,
        parameter_axes: "RMS threshold, symmetry, heavy atoms and terminal groups",
        comparison: "Retained conformers and ordering",
    },
    FutureTask {
        id: "uff_parameters",
        category: Category::ForceFields,
        reference: Reference::Rdkit,
        parameter_axes: "element/bond environments and interfragment options",
        comparison: "Parameter availability and values",
    },
    FutureTask {
        id: "uff_energy_gradient",
        category: Category::ForceFields,
        reference: Reference::Rdkit,
        parameter_axes: "fixed input geometry and constraints",
        comparison: "Energy and gradient vector",
    },
    FutureTask {
        id: "uff_optimize",
        category: Category::ForceFields,
        reference: Reference::Rdkit,
        parameter_axes: "iteration limit, tolerances, constraints and conformer selection",
        comparison: "Status, final energy and coordinates",
    },
    FutureTask {
        id: "mmff_parameters",
        category: Category::ForceFields,
        reference: Reference::Rdkit,
        parameter_axes: "MMFF94/MMFF94s and parameter families",
        comparison: "Types, charges, parameter availability and values",
    },
    FutureTask {
        id: "mmff_energy_gradient",
        category: Category::ForceFields,
        reference: Reference::Rdkit,
        parameter_axes: "variant, fixed geometry and constraints",
        comparison: "Energy and gradient vector",
    },
    FutureTask {
        id: "mmff_optimize",
        category: Category::ForceFields,
        reference: Reference::Rdkit,
        parameter_axes: "variant, iteration limit, tolerances and constraints",
        comparison: "Status, final energy and coordinates",
    },
    FutureTask {
        id: "alignment",
        category: Category::Alignment,
        reference: Reference::Rdkit,
        parameter_axes: "mapping, weights, reflection, conformer IDs and iteration limit",
        comparison: "RMSD, transform and coordinates",
    },
    FutureTask {
        id: "best_rmsd",
        category: Category::Alignment,
        reference: Reference::Rdkit,
        parameter_axes: "symmetry, terminal groups, mapping limits and H policy",
        comparison: "Best mapping and RMSD",
    },
    FutureTask {
        id: "tautomer_enumeration",
        category: Category::Tautomers,
        reference: Reference::Rdkit,
        parameter_axes: "default/v1 transforms, enumeration limits, stereo/isotope retention and reassignment",
        comparison: "Ordered results, modified atoms/bonds and completion",
    },
    FutureTask {
        id: "tautomer_canonicalization",
        category: Category::Tautomers,
        reference: Reference::Rdkit,
        parameter_axes: "scoring profiles and enumeration policy",
        comparison: "Scores and canonical selected structure",
    },
    FutureTask {
        id: "inchi_generate",
        category: Category::Identifiers,
        reference: Reference::Rdkit,
        parameter_axes: "standard/nonstandard, stereo/isotope/fixed-H options",
        comparison: "InChI, AuxInfo, return status and messages",
    },
    FutureTask {
        id: "inchi_read",
        category: Category::Identifiers,
        reference: Reference::Rdkit,
        parameter_axes: "sanitize, remove H and supported parse options",
        comparison: "Molecule state and errors",
    },
    FutureTask {
        id: "inchi_key",
        category: Category::Identifiers,
        reference: Reference::Rdkit,
        parameter_axes: "valid/invalid InChI and supported key modes",
        comparison: "Exact key and outcome",
    },
    FutureTask {
        id: "pdb_read",
        category: Category::Bio,
        reference: Reference::Gemmi,
        parameter_axes: "models, altloc, connectivity, element inference and metadata",
        comparison: "Typed hierarchy, source IDs, coordinates and metadata",
    },
    FutureTask {
        id: "mmcif_read",
        category: Category::Bio,
        reference: Reference::Gemmi,
        parameter_axes: "models, altloc, author/label IDs, entities, assemblies and missing values",
        comparison: "Typed hierarchy, identifiers, coordinates and metadata",
    },
    FutureTask {
        id: "bio_write",
        category: Category::Bio,
        reference: Reference::Gemmi,
        parameter_axes: "PDB/mmCIF, model/altloc selection and formatting policy",
        comparison: "Source-backed preserved fields and roundtrip",
    },
    FutureTask {
        id: "bio_selection",
        category: Category::Bio,
        reference: Reference::Gemmi,
        parameter_axes: "protein/DNA/RNA/water/ligand, chains, residues, models and altloc",
        comparison: "Selected original IDs and topology",
    },
    FutureTask {
        id: "bio_sequence",
        category: Category::Bio,
        reference: Reference::Gemmi,
        parameter_axes: "protein/DNA/RNA, modified residues, missing/unknown residues",
        comparison: "Sequence, residue classification and source mapping",
    },
    FutureTask {
        id: "bio_geometry",
        category: Category::Bio,
        reference: Reference::Gemmi,
        parameter_axes: "distances, neighbors, crystal/assembly transforms and selection",
        comparison: "Coordinates, contacts and mappings",
    },
    FutureTask {
        id: "batch",
        category: Category::Composition,
        reference: Reference::RustScalar,
        parameter_axes: "same chemistry profiles, invalid-record policy, ordering and worker count",
        comparison: "Every indexed result/error versus scalar Rust",
    },
    FutureTask {
        id: "operation_sequences",
        category: Category::Composition,
        reference: Reference::Rdkit,
        parameter_axes: "explicit source-backed operation sequences",
        comparison: "State after each stage, not only the last result",
    },
    FutureTask {
        id: "concurrent_reads",
        category: Category::Composition,
        reference: Reference::RustScalar,
        parameter_axes: "declared read operations and fixed schedules",
        comparison: "Same outputs as scalar Rust, unchanged input",
    },
    FutureTask {
        id: "confseq",
        category: Category::Native,
        reference: Reference::NoExternalReference,
        parameter_axes: "native sequence/geometry parameters",
        comparison: "Local correctness and binding consistency; not RDKit parity",
    },
];

impl FutureTask {
    pub const fn profile_count(&self) -> Option<usize> {
        None
    }
}

impl TaskId {
    pub const fn name(self) -> &'static str {
        match self {
            Self::SmilesRead => "smiles_read",
            Self::Sanitize => "sanitize",
            Self::Kekulize => "kekulize",
            Self::MolecularWeight => "molecular_weight",
            Self::ExactMolecularWeight => "exact_molecular_weight",
            Self::MolecularFormula => "molecular_formula",
            Self::AddHydrogens => "add_hydrogens",
            Self::RemoveHydrogens => "remove_hydrogens",
            Self::CipLabels => "cip_labels",
            Self::PotentialStereo => "potential_stereo",
            Self::Coordinates2d => "coordinates_2d",
            Self::Valence => "valence",
            Self::DistanceMatrix => "distance_matrix",
            Self::NumHeavyAtoms => "num_heavy_atoms",
            Self::TotalAtomCount => "total_atom_count",
            Self::LipinskiHBA => "lipinski_hba",
            Self::LipinskiHBD => "lipinski_hbd",
            Self::FractionCSP3 => "fraction_csp3",
            Self::NumRings => "num_rings",
            Self::NumHeterocycles => "num_heterocycles",
            Self::NumAromaticRings => "num_aromatic_rings",
            Self::NumSaturatedRings => "num_saturated_rings",
            Self::NumAliphaticRings => "num_aliphatic_rings",
            Self::NumAromaticHeterocycles => "num_aromatic_heterocycles",
            Self::NumAromaticCarbocycles => "num_aromatic_carbocycles",
            Self::NumAliphaticHeterocycles => "num_aliphatic_heterocycles",
            Self::NumAliphaticCarbocycles => "num_aliphatic_carbocycles",
            Self::NumSaturatedHeterocycles => "num_saturated_heterocycles",
            Self::NumSaturatedCarbocycles => "num_saturated_carbocycles",
        }
    }
    pub const fn category(self) -> Category {
        match self {
            Self::SmilesRead => Category::Notation,
            Self::Sanitize
            | Self::Kekulize
            | Self::AddHydrogens
            | Self::RemoveHydrogens
            | Self::Valence => Category::Chemistry,
            Self::MolecularWeight | Self::ExactMolecularWeight | Self::MolecularFormula => {
                Category::Descriptors
            }
            Self::NumHeavyAtoms => Category::Descriptors,
            Self::TotalAtomCount => Category::Descriptors,
            Self::LipinskiHBA | Self::LipinskiHBD | Self::FractionCSP3 => Category::Descriptors,
            Self::NumRings
            | Self::NumHeterocycles
            | Self::NumAromaticRings
            | Self::NumSaturatedRings
            | Self::NumAliphaticRings
            | Self::NumAromaticHeterocycles
            | Self::NumAromaticCarbocycles
            | Self::NumAliphaticHeterocycles
            | Self::NumAliphaticCarbocycles
            | Self::NumSaturatedHeterocycles
            | Self::NumSaturatedCarbocycles => Category::Descriptors,
            Self::CipLabels | Self::PotentialStereo => Category::Stereo,
            Self::Coordinates2d => Category::Depiction,
            Self::DistanceMatrix => Category::Chemistry,
        }
    }
}

#[cfg(test)]
mod catalog_tests {
    use super::*;
    #[test]
    fn future_catalog_has_unique_ids_and_no_fabricated_counts() {
        for (i, row) in FUTURE_TASKS.iter().enumerate() {
            assert!(
                !row.id.is_empty() && !row.parameter_axes.is_empty() && !row.comparison.is_empty()
            );
            assert!(
                !FUTURE_TASKS[..i]
                    .iter()
                    .any(|previous| previous.id == row.id)
            );
            assert_eq!(row.profile_count(), None);
        }
    }
    #[test]
    fn native_and_consistency_tasks_do_not_claim_rdkit_parity() {
        for row in FUTURE_TASKS {
            if row.category == Category::Native {
                assert_eq!(row.reference, Reference::NoExternalReference);
            }
            if row.id == "batch" || row.id == "concurrent_reads" {
                assert_eq!(row.reference, Reference::RustScalar);
            }
            if row.id == "binary_roundtrip" {
                assert_eq!(row.reference, Reference::Roundtrip);
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum TaskId {
    DistanceMatrix,
    SmilesRead,
    Sanitize,
    Kekulize,
    MolecularWeight,
    ExactMolecularWeight,
    MolecularFormula,
    AddHydrogens,
    RemoveHydrogens,
    CipLabels,
    PotentialStereo,
    Coordinates2d,
    Valence,
    NumHeavyAtoms,
    TotalAtomCount,
    LipinskiHBA,
    LipinskiHBD,
    FractionCSP3,
    NumRings,
    NumHeterocycles,
    NumAromaticRings,
    NumSaturatedRings,
    NumAliphaticRings,
    NumAromaticHeterocycles,
    NumAromaticCarbocycles,
    NumAliphaticHeterocycles,
    NumAliphaticCarbocycles,
    NumSaturatedHeterocycles,
    NumSaturatedCarbocycles,
}

/// `None` and an explicitly empty selection must never be conflated.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum CipSelection {
    All,
    FirstTaggedAtom,
    FirstStereoBond,
    Empty,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Profile {
    DistanceMatrix {
        use_bond_order: bool,
        use_atom_weights: bool,
    },
    SmilesRead {
        sanitize: bool,
        remove_hydrogens: bool,
    },
    SanitizeAll,
    Kekulize {
        clear_aromatic_flags: bool,
    },
    MolecularWeight {
        only_heavy: bool,
    },
    ExactMolecularWeight {
        only_heavy: bool,
    },
    MolecularFormula {
        separate_isotopes: bool,
        abbreviate_h_isotopes: bool,
    },
    AddHydrogens {
        explicit_only: bool,
    },
    RemoveHydrogens {
        sanitize: bool,
    },
    CipLabels {
        selection: CipSelection,
        max_recursive_iterations: u32,
    },
    PotentialStereo {
        clean: bool,
        flag_possible: bool,
        allow_nontetrahedral: bool,
    },
    Coordinates2dDefault,
    Valence {
        strict: bool,
    },
    NumHeavyAtoms {
        remove_hydrogens: bool,
    },
    TotalAtomCount {
        remove_hydrogens: bool,
    },
    LipinskiHBA {
        remove_hydrogens: bool,
    },
    LipinskiHBD {
        remove_hydrogens: bool,
    },
    FractionCSP3 {
        remove_hydrogens: bool,
    },
    NumRings {
        remove_hydrogens: bool,
    },
    NumHeterocycles {
        remove_hydrogens: bool,
    },
    NumAromaticRings {
        remove_hydrogens: bool,
    },
    NumSaturatedRings {
        remove_hydrogens: bool,
    },
    NumAliphaticRings {
        remove_hydrogens: bool,
    },
    NumAromaticHeterocycles {
        remove_hydrogens: bool,
    },
    NumAromaticCarbocycles {
        remove_hydrogens: bool,
    },
    NumAliphaticHeterocycles {
        remove_hydrogens: bool,
    },
    NumAliphaticCarbocycles {
        remove_hydrogens: bool,
    },
    NumSaturatedHeterocycles {
        remove_hydrogens: bool,
    },
    NumSaturatedCarbocycles {
        remove_hydrogens: bool,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputState {
    SmilesText,
    UnsanitizedHydrogensRetained,
    SanitizedHydrogensRemoved,
    SanitizedThenAddAllHydrogens,
    SanitizedHydrogensPerProfile,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Comparison {
    MatrixBits,
    TopologyAndOutcome,
    Float64Bits,
    ExactText,
    CipLabelsAndOutcome,
    StereoInfoAndCleanedTopology,
    CoordinatesAndTopology,
    ValenceRowsAndOutcome,
    Unsigned,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Prerequisite {
    MolecularPipeline,
    PublicValenceReadoutAndMolecularPipeline,
}

#[derive(Debug, Clone, Copy)]
pub struct Task {
    pub id: TaskId,
    pub input: InputState,
    pub comparison: Comparison,
    pub prerequisite: Prerequisite,
}

use Comparison::*;
use InputState::*;
use Prerequisite::*;
use TaskId::*;

/// Frozen matrices. Executable registration lives in the parent module;
/// unregistered matrices remain planned, not passing evidence.
pub const TASKS: &[Task] = &[
    Task {
        id: DistanceMatrix,
        input: SanitizedHydrogensRemoved,
        comparison: MatrixBits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SmilesRead,
        input: SmilesText,
        comparison: TopologyAndOutcome,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Sanitize,
        input: UnsanitizedHydrogensRetained,
        comparison: TopologyAndOutcome,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Kekulize,
        input: SanitizedHydrogensRemoved,
        comparison: TopologyAndOutcome,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: MolecularWeight,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: ExactMolecularWeight,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: MolecularFormula,
        input: SanitizedHydrogensRemoved,
        comparison: ExactText,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumHeavyAtoms,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: TotalAtomCount,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: LipinskiHBA,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: LipinskiHBD,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: FractionCSP3,
        input: SanitizedHydrogensPerProfile,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumRings,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumHeterocycles,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumAromaticRings,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumSaturatedRings,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumAliphaticRings,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumAromaticHeterocycles,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumAromaticCarbocycles,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumAliphaticHeterocycles,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumAliphaticCarbocycles,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumSaturatedHeterocycles,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumSaturatedCarbocycles,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: AddHydrogens,
        input: SanitizedHydrogensRemoved,
        comparison: TopologyAndOutcome,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: RemoveHydrogens,
        input: SanitizedThenAddAllHydrogens,
        comparison: TopologyAndOutcome,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: CipLabels,
        input: SanitizedHydrogensRemoved,
        comparison: CipLabelsAndOutcome,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: PotentialStereo,
        input: SanitizedHydrogensRemoved,
        comparison: StereoInfoAndCleanedTopology,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Coordinates2d,
        input: SanitizedHydrogensRemoved,
        comparison: CoordinatesAndTopology,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Valence,
        input: UnsanitizedHydrogensRetained,
        comparison: ValenceRowsAndOutcome,
        prerequisite: PublicValenceReadoutAndMolecularPipeline,
    },
];

impl TaskId {
    /// All combinations below are deliberate; there is no implicit sampling.
    /// Unlisted API options are pinned in TASKS.md, not silently expanded.
    pub fn profiles(self) -> Vec<Profile> {
        let booleans = [false, true];
        match self {
            DistanceMatrix => booleans
                .into_iter()
                .flat_map(|use_bond_order| {
                    booleans.map(move |use_atom_weights| Profile::DistanceMatrix {
                        use_bond_order,
                        use_atom_weights,
                    })
                })
                .collect(),
            SmilesRead => booleans
                .into_iter()
                .flat_map(|sanitize| {
                    booleans.map(move |remove_hydrogens| Profile::SmilesRead {
                        sanitize,
                        remove_hydrogens,
                    })
                })
                .collect(),
            Sanitize => vec![Profile::SanitizeAll],
            Kekulize => booleans
                .map(|clear_aromatic_flags| Profile::Kekulize {
                    clear_aromatic_flags,
                })
                .into(),
            MolecularWeight => booleans
                .map(|only_heavy| Profile::MolecularWeight { only_heavy })
                .into(),
            ExactMolecularWeight => booleans
                .map(|only_heavy| Profile::ExactMolecularWeight { only_heavy })
                .into(),
            MolecularFormula => booleans
                .into_iter()
                .flat_map(|separate_isotopes| {
                    booleans.map(move |abbreviate_h_isotopes| Profile::MolecularFormula {
                        separate_isotopes,
                        abbreviate_h_isotopes,
                    })
                })
                .collect(),
            AddHydrogens => booleans
                .map(|explicit_only| Profile::AddHydrogens { explicit_only })
                .into(),
            RemoveHydrogens => booleans
                .map(|sanitize| Profile::RemoveHydrogens { sanitize })
                .into(),
            CipLabels => vec![
                Profile::CipLabels {
                    selection: CipSelection::All,
                    max_recursive_iterations: 0,
                },
                Profile::CipLabels {
                    selection: CipSelection::FirstTaggedAtom,
                    max_recursive_iterations: 0,
                },
                Profile::CipLabels {
                    selection: CipSelection::FirstStereoBond,
                    max_recursive_iterations: 0,
                },
                Profile::CipLabels {
                    selection: CipSelection::Empty,
                    max_recursive_iterations: 1,
                },
            ],
            PotentialStereo => booleans
                .into_iter()
                .flat_map(|clean| {
                    booleans.into_iter().flat_map(move |flag_possible| {
                        booleans.map(move |allow_nontetrahedral| Profile::PotentialStereo {
                            clean,
                            flag_possible,
                            allow_nontetrahedral,
                        })
                    })
                })
                .collect(),
            Coordinates2d => vec![Profile::Coordinates2dDefault],
            Valence => booleans.map(|strict| Profile::Valence { strict }).into(),
            NumHeavyAtoms => booleans
                .map(|remove_hydrogens| Profile::NumHeavyAtoms { remove_hydrogens })
                .into(),
            TotalAtomCount => booleans
                .map(|remove_hydrogens| Profile::TotalAtomCount { remove_hydrogens })
                .into(),
            LipinskiHBA => booleans
                .map(|remove_hydrogens| Profile::LipinskiHBA { remove_hydrogens })
                .into(),
            LipinskiHBD => booleans
                .map(|remove_hydrogens| Profile::LipinskiHBD { remove_hydrogens })
                .into(),
            FractionCSP3 => booleans
                .map(|remove_hydrogens| Profile::FractionCSP3 { remove_hydrogens })
                .into(),
            NumRings => booleans
                .map(|remove_hydrogens| Profile::NumRings { remove_hydrogens })
                .into(),
            NumHeterocycles => booleans
                .map(|remove_hydrogens| Profile::NumHeterocycles { remove_hydrogens })
                .into(),
            NumAromaticRings => booleans
                .map(|remove_hydrogens| Profile::NumAromaticRings { remove_hydrogens })
                .into(),
            NumSaturatedRings => booleans
                .map(|remove_hydrogens| Profile::NumSaturatedRings { remove_hydrogens })
                .into(),
            NumAliphaticRings => booleans
                .map(|remove_hydrogens| Profile::NumAliphaticRings { remove_hydrogens })
                .into(),
            NumAromaticHeterocycles => booleans
                .map(|remove_hydrogens| Profile::NumAromaticHeterocycles { remove_hydrogens })
                .into(),
            NumAromaticCarbocycles => booleans
                .map(|remove_hydrogens| Profile::NumAromaticCarbocycles { remove_hydrogens })
                .into(),
            NumAliphaticHeterocycles => booleans
                .map(|remove_hydrogens| Profile::NumAliphaticHeterocycles { remove_hydrogens })
                .into(),
            NumAliphaticCarbocycles => booleans
                .map(|remove_hydrogens| Profile::NumAliphaticCarbocycles { remove_hydrogens })
                .into(),
            NumSaturatedHeterocycles => booleans
                .map(|remove_hydrogens| Profile::NumSaturatedHeterocycles { remove_hydrogens })
                .into(),
            NumSaturatedCarbocycles => booleans
                .map(|remove_hydrogens| Profile::NumSaturatedCarbocycles { remove_hydrogens })
                .into(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn molecular_plan_has_unique_tasks_and_profiles() {
        for (i, task) in TASKS.iter().enumerate() {
            assert!(!TASKS[..i].iter().any(|other| other.id == task.id));
            let profiles = task.id.profiles();
            assert!(!profiles.is_empty());
            for (j, profile) in profiles.iter().enumerate() {
                assert!(!profiles[..j].contains(profile));
            }
        }
    }

    #[test]
    fn molecular_plan_counts_are_explicit() {
        assert_eq!(
            TASKS
                .iter()
                .map(|t| t.id.profiles().len())
                .collect::<Vec<_>>(),
            // RING-LIVE-PUBLIC T1: eleven ring-state tasks, each with the
            // two explicit remove_hydrogens profiles.
            [
                4, 4, 1, 2, 2, 2, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 8, 1,
                2
            ]
        );
        assert_eq!(
            TASKS
                .iter()
                .filter(|t| t.prerequisite == PublicValenceReadoutAndMolecularPipeline)
                .map(|t| t.id)
                .collect::<Vec<_>>(),
            [Valence]
        );
    }

    #[test]
    fn molecular_plan_does_not_silently_register_unimplemented_runners() {
        let executable = super::super::select(None).unwrap();
        assert_eq!(executable.len(), 28);
        assert_eq!(executable[0].operation, super::super::Operation::FuzzyAnd);
        assert_eq!(executable[1].operation, super::super::Operation::FuzzyOr);
    }
}
