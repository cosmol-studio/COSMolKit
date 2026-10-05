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
    // Ordinary matching is registered with explicit typed profiles in search.rs.
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
        id: "svg_options",
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
        id: "uff_optimize_options",
        category: Category::ForceFields,
        reference: Reference::Rdkit,
        parameter_axes: "additional tolerances, constraints and multi-conformer/thread dispatch beyond the executable single-conformer matrix",
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
            Self::TautomerEnumeration => "tautomer_enumeration",
            Self::TautomerCanonicalization => "tautomer_canonicalization",
            Self::Chi0 => "chi_0",
            Self::Chi1 => "chi_1",
            Self::HallKierAlpha => "hall_kier_alpha",
            Self::HallKierAlphaWithContributions => "hall_kier_alpha_with_contributions",
            Self::Kappa1 => "kappa_1",
            Self::Kappa2 => "kappa_2",
            Self::Kappa3 => "kappa_3",
            Self::Phi => "phi",
            Self::Mqns => "mqns",
            Self::Chi0V => "chi_0_v",
            Self::Chi1V => "chi_1_v",
            Self::Chi2V => "chi_2_v",
            Self::Chi3V => "chi_3_v",
            Self::Chi4V => "chi_4_v",
            Self::Chi0N => "chi_0_n",
            Self::Chi1N => "chi_1_n",
            Self::Chi2N => "chi_2_n",
            Self::Chi3N => "chi_3_n",
            Self::Chi4N => "chi_4_n",
            Self::ChiNV => "chi_n_v",
            Self::ChiNN => "chi_n_n",
            Self::NumAmideBonds => "num_amide_bonds",
            Self::NumSpiroAtoms => "num_spiro_atoms",
            Self::NumBridgeheadAtoms => "num_bridgehead_atoms",
            Self::NumAtomStereoCenters => "num_atom_stereo_centers",
            Self::NumUnspecifiedAtomStereoCenters => "num_unspecified_atom_stereo_centers",
            Self::NumRotatableBonds => "num_rotatable_bonds",
            Self::CrippenDescriptors => "crippen_descriptors",
            Self::LabuteAsa => "labute_asa",
            Self::LabuteAsaContributions => "labute_asa_contributions",
            Self::Tpsa => "tpsa",
            Self::SlogpVsa => "slogp_vsa",
            Self::SmrVsa => "smr_vsa",
            Self::SlogpVsa1 => "slogp_vsa_1",
            Self::SlogpVsa2 => "slogp_vsa_2",
            Self::SlogpVsa3 => "slogp_vsa_3",
            Self::SlogpVsa4 => "slogp_vsa_4",
            Self::SlogpVsa5 => "slogp_vsa_5",
            Self::SlogpVsa6 => "slogp_vsa_6",
            Self::SlogpVsa7 => "slogp_vsa_7",
            Self::SlogpVsa8 => "slogp_vsa_8",
            Self::SlogpVsa9 => "slogp_vsa_9",
            Self::SlogpVsa10 => "slogp_vsa_10",
            Self::SlogpVsa11 => "slogp_vsa_11",
            Self::SlogpVsa12 => "slogp_vsa_12",
            Self::SmrVsa1 => "smr_vsa_1",
            Self::SmrVsa2 => "smr_vsa_2",
            Self::SmrVsa3 => "smr_vsa_3",
            Self::SmrVsa4 => "smr_vsa_4",
            Self::SmrVsa5 => "smr_vsa_5",
            Self::SmrVsa6 => "smr_vsa_6",
            Self::SmrVsa7 => "smr_vsa_7",
            Self::SmrVsa8 => "smr_vsa_8",
            Self::SmrVsa9 => "smr_vsa_9",
            Self::SmrVsa10 => "smr_vsa_10",
            Self::Qed => "qed",
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
            Self::Svg => "svg",
            Self::Valence => "valence",
            Self::DistanceMatrix => "distance_matrix",
            Self::NumHeavyAtoms => "num_heavy_atoms",
            Self::TotalAtomCount => "total_atom_count",
            Self::LipinskiHBA => "lipinski_hba",
            Self::LipinskiHBD => "lipinski_hbd",
            Self::FractionCSP3 => "fraction_csp3",
            Self::NumHeteroatoms => "num_heteroatoms",
            Self::NumHba => "num_hba",
            Self::NumHbd => "num_hbd",
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
            Self::MorganFingerprint => "morgan_fingerprint",
            Self::MorganSparseFingerprint => "morgan_sparse_fingerprint",
            Self::MorganCountFingerprint => "morgan_count_fingerprint",
            Self::MorganSparseCountFingerprint => "morgan_sparse_count_fingerprint",
        }
    }
    pub const fn category(self) -> Category {
        match self {
            Self::TautomerEnumeration | Self::TautomerCanonicalization => Category::Tautomers,
            Self::Chi0
            | Self::Chi1
            | Self::HallKierAlpha
            | Self::HallKierAlphaWithContributions
            | Self::Kappa1
            | Self::Kappa2
            | Self::Kappa3
            | Self::Phi
            | Self::Mqns
            | Self::Chi0V
            | Self::Chi1V
            | Self::Chi2V
            | Self::Chi3V
            | Self::Chi4V
            | Self::Chi0N
            | Self::Chi1N
            | Self::Chi2N
            | Self::Chi3N
            | Self::Chi4N
            | Self::ChiNV
            | Self::ChiNN => Category::Descriptors,
            Self::NumAmideBonds
            | Self::NumSpiroAtoms
            | Self::NumBridgeheadAtoms
            | Self::NumAtomStereoCenters
            | Self::NumUnspecifiedAtomStereoCenters
            | Self::NumRotatableBonds
            | Self::CrippenDescriptors
            | Self::LabuteAsa
            | Self::LabuteAsaContributions
            | Self::Tpsa
            | Self::SlogpVsa
            | Self::SmrVsa
            | Self::SlogpVsa1
            | Self::SlogpVsa2
            | Self::SlogpVsa3
            | Self::SlogpVsa4
            | Self::SlogpVsa5
            | Self::SlogpVsa6
            | Self::SlogpVsa7
            | Self::SlogpVsa8
            | Self::SlogpVsa9
            | Self::SlogpVsa10
            | Self::SlogpVsa11
            | Self::SlogpVsa12
            | Self::SmrVsa1
            | Self::SmrVsa2
            | Self::SmrVsa3
            | Self::SmrVsa4
            | Self::SmrVsa5
            | Self::SmrVsa6
            | Self::SmrVsa7
            | Self::SmrVsa8
            | Self::SmrVsa9
            | Self::SmrVsa10
            | Self::Qed => Category::Descriptors,
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
            Self::NumHeteroatoms => Category::Descriptors,
            Self::NumHba => Category::Descriptors,
            Self::NumHbd => Category::Descriptors,
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
            Self::Coordinates2d | Self::Svg => Category::Depiction,
            Self::DistanceMatrix => Category::Chemistry,
            Self::MorganFingerprint
            | Self::MorganSparseFingerprint
            | Self::MorganCountFingerprint
            | Self::MorganSparseCountFingerprint => Category::FingerprintGeneration,
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
                super::super::TASKS
                    .iter()
                    .all(|task| task.operation.name() != row.id),
                "executable task still listed as future: {}",
                row.id
            );
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
    TautomerEnumeration,
    TautomerCanonicalization,
    NumAmideBonds,
    NumSpiroAtoms,
    NumBridgeheadAtoms,
    NumAtomStereoCenters,
    NumUnspecifiedAtomStereoCenters,
    NumRotatableBonds,
    CrippenDescriptors,
    LabuteAsa,
    LabuteAsaContributions,
    Tpsa,
    SlogpVsa,
    SmrVsa,
    SlogpVsa1,
    SlogpVsa2,
    SlogpVsa3,
    SlogpVsa4,
    SlogpVsa5,
    SlogpVsa6,
    SlogpVsa7,
    SlogpVsa8,
    SlogpVsa9,
    SlogpVsa10,
    SlogpVsa11,
    SlogpVsa12,
    SmrVsa1,
    SmrVsa2,
    SmrVsa3,
    SmrVsa4,
    SmrVsa5,
    SmrVsa6,
    SmrVsa7,
    SmrVsa8,
    SmrVsa9,
    SmrVsa10,
    Qed,

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
    Svg,
    Valence,
    NumHeavyAtoms,
    TotalAtomCount,
    LipinskiHBA,
    LipinskiHBD,
    FractionCSP3,
    NumHeteroatoms,
    NumHba,
    NumHbd,
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
    MorganFingerprint,
    MorganSparseFingerprint,
    MorganCountFingerprint,
    MorganSparseCountFingerprint,
    Chi0,
    Chi1,
    HallKierAlpha,
    HallKierAlphaWithContributions,
    Kappa1,
    Kappa2,
    Kappa3,
    Phi,
    Mqns,
    Chi0V,
    Chi1V,
    Chi2V,
    Chi3V,
    Chi4V,
    Chi0N,
    Chi1N,
    Chi2N,
    Chi3N,
    Chi4N,
    ChiNV,
    ChiNN,
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
pub enum MorganOutputKind {
    DenseBits,
    SparseBits,
    HashedCounts,
    SparseCounts,
}

impl MorganOutputKind {
    pub const fn task_name(self) -> &'static str {
        match self {
            Self::DenseBits => "morgan_fingerprint",
            Self::SparseBits => "morgan_sparse_fingerprint",
            Self::HashedCounts => "morgan_count_fingerprint",
            Self::SparseCounts => "morgan_sparse_count_fingerprint",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum MorganInvariantKind {
    Connectivity,
    Features,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum TautomerCatalog {
    Current,
    V1,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TautomerProfile {
    pub catalog: TautomerCatalog,
    pub max_tautomers: u32,
    pub max_transforms: u32,
    pub remove_sp3_stereo: bool,
    pub remove_bond_stereo: bool,
    pub remove_isotopic_hydrogens: bool,
    pub reassign_stereo: bool,
}
impl TautomerProfile {
    pub const fn branch(self) -> &'static str {
        match self.catalog {
            TautomerCatalog::Current => "default",
            TautomerCatalog::V1 => "v1",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum RotatableBondMode {
    Default,
    NonStrict,
    Strict,
    StrictLinkages,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum VsaBins {
    Default,
    CustomDuplicates,
}
impl VsaBins {
    pub const fn values(self) -> Option<&'static [f64]> {
        match self {
            Self::Default => None,
            Self::CustomDuplicates => Some(&[-0.2, 0.0, 0.25, 0.25, 0.8]),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Profile {
    TautomerEnumeration {
        parameters: TautomerProfile,
    },
    TautomerCanonicalization {
        parameters: TautomerProfile,
    },
    NumAmideBonds,
    NumSpiroAtoms,
    NumBridgeheadAtoms,
    NumAtomStereoCenters,
    NumUnspecifiedAtomStereoCenters,
    NumRotatableBonds {
        mode: Option<RotatableBondMode>,
    },
    CrippenDescriptors {
        include_hydrogens: Option<bool>,
        force: bool,
    },
    LabuteAsa {
        include_hydrogens: Option<bool>,
        force: bool,
    },
    LabuteAsaContributions {
        include_hydrogens: Option<bool>,
        force: bool,
    },
    Tpsa {
        include_sulfur_phosphorus: Option<bool>,
        force: bool,
    },
    SlogpVsa {
        bins: VsaBins,
        force: Option<bool>,
    },
    SmrVsa {
        bins: VsaBins,
        force: Option<bool>,
    },
    SlogpVsa1,
    SlogpVsa2,
    SlogpVsa3,
    SlogpVsa4,
    SlogpVsa5,
    SlogpVsa6,
    SlogpVsa7,
    SlogpVsa8,
    SlogpVsa9,
    SlogpVsa10,
    SlogpVsa11,
    SlogpVsa12,
    SmrVsa1,
    SmrVsa2,
    SmrVsa3,
    SmrVsa4,
    SmrVsa5,
    SmrVsa6,
    SmrVsa7,
    SmrVsa8,
    SmrVsa9,
    SmrVsa10,
    Qed,
    Chi0VWithParams {
        force: bool,
    },
    Chi1VWithParams {
        force: bool,
    },
    Chi2VWithParams {
        force: bool,
    },
    Chi3VWithParams {
        force: bool,
    },
    Chi4VWithParams {
        force: bool,
    },
    Chi0NWithParams {
        force: bool,
    },
    Chi1NWithParams {
        force: bool,
    },
    Chi2NWithParams {
        force: bool,
    },
    Chi3NWithParams {
        force: bool,
    },
    Chi4NWithParams {
        force: bool,
    },
    ChiNVWithParams {
        order: u32,
        force: bool,
    },
    ChiNNWithParams {
        order: u32,
        force: bool,
    },
    LabuteAsaCacheSequence,
    SlogpVsaCacheSequence,
    SmrVsaCacheSequence,
    ChiNVCacheSequence,
    ChiNNCacheSequence,

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
    Chi0,
    Chi1,
    HallKierAlpha,
    HallKierAlphaWithContributions,
    Kappa1,
    Kappa2,
    Kappa3,
    Phi,
    Mqns {
        force: bool,
    },
    Chi0V,
    Chi1V,
    Chi2V,
    Chi3V,
    Chi4V,
    Chi0N,
    Chi1N,
    Chi2N,
    Chi3N,
    Chi4N,
    ChiNV {
        order: u32,
    },
    ChiNN {
        order: u32,
    },
    Coordinates2dDefault,
    /// The only drawing profile: 300x300, default source preparation.
    SvgDefault,
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
    NumHeteroatoms {
        remove_hydrogens: bool,
    },
    NumHba {
        remove_hydrogens: bool,
    },
    NumHbd {
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
    Morgan {
        output: MorganOutputKind,
        radius: u32,
        include_chirality: bool,
        invariants: MorganInvariantKind,
        count_simulation: bool,
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
    TautomerFullEnumeration,
    TautomerFullCanonicalization,
    MatrixBits,
    TopologyAndOutcome,
    Float64Bits,
    Float64PairBits,
    Float64VectorBits,
    LabuteAsaContributionsBits,
    Float64ContributionsBits,
    UnsignedVector,
    ExactText,
    SvgText,
    CipLabelsAndOutcome,
    StereoInfoAndCleanedTopology,
    CoordinatesAndTopology,
    ValenceRowsAndOutcome,
    Unsigned,
    MorganFingerprintAndAdditionalOutput,
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
        id: NumHeteroatoms,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumHba,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumHbd,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
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
        id: Svg,
        input: SanitizedHydrogensRemoved,
        comparison: SvgText,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Valence,
        input: UnsanitizedHydrogensRetained,
        comparison: ValenceRowsAndOutcome,
        prerequisite: PublicValenceReadoutAndMolecularPipeline,
    },
    Task {
        id: MorganFingerprint,
        input: SanitizedHydrogensRemoved,
        comparison: MorganFingerprintAndAdditionalOutput,
        prerequisite: PublicValenceReadoutAndMolecularPipeline,
    },
    Task {
        id: MorganSparseFingerprint,
        input: SanitizedHydrogensRemoved,
        comparison: MorganFingerprintAndAdditionalOutput,
        prerequisite: PublicValenceReadoutAndMolecularPipeline,
    },
    Task {
        id: MorganCountFingerprint,
        input: SanitizedHydrogensRemoved,
        comparison: MorganFingerprintAndAdditionalOutput,
        prerequisite: PublicValenceReadoutAndMolecularPipeline,
    },
    Task {
        id: MorganSparseCountFingerprint,
        input: SanitizedHydrogensRemoved,
        comparison: MorganFingerprintAndAdditionalOutput,
        prerequisite: PublicValenceReadoutAndMolecularPipeline,
    },
    Task {
        id: Chi0,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Chi1,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: HallKierAlpha,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: HallKierAlphaWithContributions,
        input: SanitizedHydrogensRemoved,
        comparison: Float64ContributionsBits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Kappa1,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Kappa2,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Kappa3,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Phi,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Mqns,
        input: SanitizedHydrogensRemoved,
        comparison: UnsignedVector,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Chi0V,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Chi1V,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Chi2V,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Chi3V,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Chi4V,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Chi0N,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Chi1N,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Chi2N,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Chi3N,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Chi4N,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: ChiNV,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: ChiNN,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: TautomerEnumeration,
        input: SanitizedHydrogensRemoved,
        comparison: TautomerFullEnumeration,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: TautomerCanonicalization,
        input: SanitizedHydrogensRemoved,
        comparison: TautomerFullCanonicalization,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumAmideBonds,
        input: SanitizedHydrogensRemoved,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumSpiroAtoms,
        input: SanitizedHydrogensRemoved,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumBridgeheadAtoms,
        input: SanitizedHydrogensRemoved,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumAtomStereoCenters,
        input: SanitizedHydrogensRemoved,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumUnspecifiedAtomStereoCenters,
        input: SanitizedHydrogensRemoved,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumRotatableBonds,
        input: SanitizedHydrogensRemoved,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: CrippenDescriptors,
        input: SanitizedHydrogensRemoved,
        comparison: Float64PairBits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: LabuteAsa,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: LabuteAsaContributions,
        input: SanitizedHydrogensRemoved,
        comparison: LabuteAsaContributionsBits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Tpsa,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SlogpVsa,
        input: SanitizedHydrogensRemoved,
        comparison: Float64VectorBits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SmrVsa,
        input: SanitizedHydrogensRemoved,
        comparison: Float64VectorBits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SlogpVsa1,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SlogpVsa2,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SlogpVsa3,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SlogpVsa4,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SlogpVsa5,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SlogpVsa6,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SlogpVsa7,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SlogpVsa8,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SlogpVsa9,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SlogpVsa10,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SlogpVsa11,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SlogpVsa12,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SmrVsa1,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SmrVsa2,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SmrVsa3,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SmrVsa4,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SmrVsa5,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SmrVsa6,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SmrVsa7,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SmrVsa8,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SmrVsa9,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SmrVsa10,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Qed,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
];

impl TaskId {
    /// All combinations below are deliberate; there is no implicit sampling.
    /// Unlisted API options are pinned in TASKS.md, not silently expanded.
    pub fn profiles(self) -> Vec<Profile> {
        let booleans = [false, true];
        match self {
            TautomerEnumeration | TautomerCanonicalization => {
                [TautomerCatalog::Current, TautomerCatalog::V1]
                    .map(|catalog| {
                        let parameters = TautomerProfile {
                            catalog,
                            max_tautomers: 1000,
                            max_transforms: 1000,
                            remove_sp3_stereo: true,
                            remove_bond_stereo: true,
                            remove_isotopic_hydrogens: true,
                            reassign_stereo: true,
                        };
                        if self == TautomerEnumeration {
                            Profile::TautomerEnumeration { parameters }
                        } else {
                            Profile::TautomerCanonicalization { parameters }
                        }
                    })
                    .into()
            }
            NumAmideBonds => vec![Profile::NumAmideBonds],
            NumSpiroAtoms => vec![Profile::NumSpiroAtoms],
            NumBridgeheadAtoms => vec![Profile::NumBridgeheadAtoms],
            NumAtomStereoCenters => vec![Profile::NumAtomStereoCenters],
            NumUnspecifiedAtomStereoCenters => vec![Profile::NumUnspecifiedAtomStereoCenters],
            NumRotatableBonds => vec![
                Profile::NumRotatableBonds { mode: None },
                Profile::NumRotatableBonds {
                    mode: Some(RotatableBondMode::Default),
                },
                Profile::NumRotatableBonds {
                    mode: Some(RotatableBondMode::NonStrict),
                },
                Profile::NumRotatableBonds {
                    mode: Some(RotatableBondMode::Strict),
                },
                Profile::NumRotatableBonds {
                    mode: Some(RotatableBondMode::StrictLinkages),
                },
            ],
            CrippenDescriptors => [Profile::CrippenDescriptors {
                include_hydrogens: None,
                force: false,
            }]
            .into_iter()
            .chain(booleans.into_iter().flat_map(|flag| {
                booleans.map(move |force| Profile::CrippenDescriptors {
                    include_hydrogens: Some(flag),
                    force,
                })
            }))
            .collect(),
            LabuteAsa => [Profile::LabuteAsa {
                include_hydrogens: None,
                force: false,
            }]
            .into_iter()
            .chain(booleans.into_iter().flat_map(|flag| {
                booleans.map(move |force| Profile::LabuteAsa {
                    include_hydrogens: Some(flag),
                    force,
                })
            }))
            .chain([Profile::LabuteAsaCacheSequence])
            .collect(),
            LabuteAsaContributions => [Profile::LabuteAsaContributions {
                include_hydrogens: None,
                force: false,
            }]
            .into_iter()
            .chain(booleans.into_iter().flat_map(|flag| {
                booleans.map(move |force| Profile::LabuteAsaContributions {
                    include_hydrogens: Some(flag),
                    force,
                })
            }))
            .collect(),
            Tpsa => [Profile::Tpsa {
                include_sulfur_phosphorus: None,
                force: false,
            }]
            .into_iter()
            .chain(booleans.into_iter().flat_map(|flag| {
                booleans.map(move |force| Profile::Tpsa {
                    include_sulfur_phosphorus: Some(flag),
                    force,
                })
            }))
            .collect(),
            SlogpVsa => [Profile::SlogpVsa {
                bins: VsaBins::Default,
                force: None,
            }]
            .into_iter()
            .chain(
                [VsaBins::Default, VsaBins::CustomDuplicates]
                    .into_iter()
                    .flat_map(|bins| {
                        booleans.map(move |force| Profile::SlogpVsa {
                            bins,
                            force: Some(force),
                        })
                    }),
            )
            .chain([Profile::SlogpVsaCacheSequence])
            .collect(),
            SmrVsa => [Profile::SmrVsa {
                bins: VsaBins::Default,
                force: None,
            }]
            .into_iter()
            .chain(
                [VsaBins::Default, VsaBins::CustomDuplicates]
                    .into_iter()
                    .flat_map(|bins| {
                        booleans.map(move |force| Profile::SmrVsa {
                            bins,
                            force: Some(force),
                        })
                    }),
            )
            .chain([Profile::SmrVsaCacheSequence])
            .collect(),
            SlogpVsa1 => vec![Profile::SlogpVsa1],
            SlogpVsa2 => vec![Profile::SlogpVsa2],
            SlogpVsa3 => vec![Profile::SlogpVsa3],
            SlogpVsa4 => vec![Profile::SlogpVsa4],
            SlogpVsa5 => vec![Profile::SlogpVsa5],
            SlogpVsa6 => vec![Profile::SlogpVsa6],
            SlogpVsa7 => vec![Profile::SlogpVsa7],
            SlogpVsa8 => vec![Profile::SlogpVsa8],
            SlogpVsa9 => vec![Profile::SlogpVsa9],
            SlogpVsa10 => vec![Profile::SlogpVsa10],
            SlogpVsa11 => vec![Profile::SlogpVsa11],
            SlogpVsa12 => vec![Profile::SlogpVsa12],
            SmrVsa1 => vec![Profile::SmrVsa1],
            SmrVsa2 => vec![Profile::SmrVsa2],
            SmrVsa3 => vec![Profile::SmrVsa3],
            SmrVsa4 => vec![Profile::SmrVsa4],
            SmrVsa5 => vec![Profile::SmrVsa5],
            SmrVsa6 => vec![Profile::SmrVsa6],
            SmrVsa7 => vec![Profile::SmrVsa7],
            SmrVsa8 => vec![Profile::SmrVsa8],
            SmrVsa9 => vec![Profile::SmrVsa9],
            SmrVsa10 => vec![Profile::SmrVsa10],
            Qed => vec![Profile::Qed],
            Chi0 => vec![Profile::Chi0],
            Chi1 => vec![Profile::Chi1],
            HallKierAlpha => vec![Profile::HallKierAlpha],
            HallKierAlphaWithContributions => vec![Profile::HallKierAlphaWithContributions],
            Kappa1 => vec![Profile::Kappa1],
            Kappa2 => vec![Profile::Kappa2],
            Kappa3 => vec![Profile::Kappa3],
            Phi => vec![Profile::Phi],
            Mqns => booleans.map(|force| Profile::Mqns { force }).into(),
            Chi0V => vec![
                Profile::Chi0V,
                Profile::Chi0VWithParams { force: false },
                Profile::Chi0VWithParams { force: true },
            ],
            Chi1V => vec![
                Profile::Chi1V,
                Profile::Chi1VWithParams { force: false },
                Profile::Chi1VWithParams { force: true },
            ],
            Chi2V => vec![
                Profile::Chi2V,
                Profile::Chi2VWithParams { force: false },
                Profile::Chi2VWithParams { force: true },
            ],
            Chi3V => vec![
                Profile::Chi3V,
                Profile::Chi3VWithParams { force: false },
                Profile::Chi3VWithParams { force: true },
            ],
            Chi4V => vec![
                Profile::Chi4V,
                Profile::Chi4VWithParams { force: false },
                Profile::Chi4VWithParams { force: true },
            ],
            Chi0N => vec![
                Profile::Chi0N,
                Profile::Chi0NWithParams { force: false },
                Profile::Chi0NWithParams { force: true },
            ],
            Chi1N => vec![
                Profile::Chi1N,
                Profile::Chi1NWithParams { force: false },
                Profile::Chi1NWithParams { force: true },
            ],
            Chi2N => vec![
                Profile::Chi2N,
                Profile::Chi2NWithParams { force: false },
                Profile::Chi2NWithParams { force: true },
            ],
            Chi3N => vec![
                Profile::Chi3N,
                Profile::Chi3NWithParams { force: false },
                Profile::Chi3NWithParams { force: true },
            ],
            Chi4N => vec![
                Profile::Chi4N,
                Profile::Chi4NWithParams { force: false },
                Profile::Chi4NWithParams { force: true },
            ],
            ChiNV => (0..=6)
                .map(|order| Profile::ChiNV { order })
                .chain((0..=6).flat_map(|order| {
                    booleans.map(move |force| Profile::ChiNVWithParams { order, force })
                }))
                .chain([Profile::ChiNVCacheSequence])
                .collect(),
            ChiNN => (0..=6)
                .map(|order| Profile::ChiNN { order })
                .chain((0..=6).flat_map(|order| {
                    booleans.map(move |force| Profile::ChiNNWithParams { order, force })
                }))
                .chain([Profile::ChiNNCacheSequence])
                .collect(),
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
            Svg => vec![Profile::SvgDefault],
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
            NumHeteroatoms => booleans
                .map(|remove_hydrogens| Profile::NumHeteroatoms { remove_hydrogens })
                .into(),
            NumHba => booleans
                .map(|remove_hydrogens| Profile::NumHba { remove_hydrogens })
                .into(),
            NumHbd => booleans
                .map(|remove_hydrogens| Profile::NumHbd { remove_hydrogens })
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
            MorganFingerprint
            | MorganSparseFingerprint
            | MorganCountFingerprint
            | MorganSparseCountFingerprint => {
                let output = match self {
                    MorganFingerprint => MorganOutputKind::DenseBits,
                    MorganSparseFingerprint => MorganOutputKind::SparseBits,
                    MorganCountFingerprint => MorganOutputKind::HashedCounts,
                    MorganSparseCountFingerprint => MorganOutputKind::SparseCounts,
                    _ => unreachable!("the enclosing task id is a Morgan output family"),
                };
                [2_u32, 3]
                    .into_iter()
                    .flat_map(|radius| {
                        [false, true]
                            .into_iter()
                            .flat_map(move |include_chirality| {
                                [
                                    MorganInvariantKind::Connectivity,
                                    MorganInvariantKind::Features,
                                ]
                                .into_iter()
                                .flat_map(move |invariants| {
                                    [false, true].map(move |count_simulation| Profile::Morgan {
                                        output,
                                        radius,
                                        include_chirality,
                                        invariants,
                                        count_simulation,
                                    })
                                })
                            })
                    })
                    .collect()
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parity_morgan_registry_molecular_plan_has_unique_tasks_and_profiles() {
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
    fn parity_morgan_registry_molecular_plan_counts_are_explicit() {
        assert_eq!(
            TASKS
                .iter()
                .map(|t| t.id.profiles().len())
                .collect::<Vec<_>>(),
            [
                4, 4, 1, 2, 2, 2, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                4, 8, 1, 1, 2, 16, 16, 16, 16, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 3, 3, 3, 3, 3, 3,
                3, 3, 22, 22, 2, 2, 1, 1, 1, 1, 1, 5, 5, 6, 5, 5, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
            ]
        );
        assert_eq!(
            TASKS
                .iter()
                .filter(|t| t.prerequisite == PublicValenceReadoutAndMolecularPipeline)
                .map(|t| t.id)
                .collect::<Vec<_>>(),
            [
                Valence,
                MorganFingerprint,
                MorganSparseFingerprint,
                MorganCountFingerprint,
                MorganSparseCountFingerprint
            ]
        );
    }

    #[test]
    fn parity_morgan_registry_molecular_plan_tracks_executable_registration() {
        let executable = super::super::select(None).unwrap();
        assert_eq!(executable.len(), 100);
        assert_eq!(
            executable[62].operation,
            super::super::Operation::SubstructureMatch
        );
        assert_eq!(
            executable[0].operation,
            super::super::Operation::BioPdbOutput
        );
        assert_eq!(
            executable[1].operation,
            super::super::Operation::BioPdbOutput
        );
        assert_eq!(executable[2].operation, super::super::Operation::FuzzyAnd);
        assert_eq!(executable[3].operation, super::super::Operation::FuzzyOr);
        for id in [
            MorganFingerprint,
            MorganSparseFingerprint,
            MorganCountFingerprint,
            MorganSparseCountFingerprint,
        ] {
            assert!(
                executable
                    .iter()
                    .any(|task| task.operation == super::super::Operation::Molecular(id))
            );
        }
    }
}
