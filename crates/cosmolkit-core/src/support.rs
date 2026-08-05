#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SupportStatus {
    Supported,
    SupportedWithRdkitParity { rdkit_version: &'static str },
    PreservedOnly,
    Experimental,
    Unsupported { reason: &'static str },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FeatureCategory {
    Core,
    TopologyOperation,
    Io,
    Fingerprint,
    Drawing,
    Stereo,
    Valence,
    Batch,
    BioHierarchy,
    BioCoordinate,
    BioSelection,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FeatureSpec {
    pub name: &'static str,
    pub category: FeatureCategory,
    pub status: SupportStatus,
    pub parity_sensitive: bool,
    pub docs: &'static str,
}

impl FeatureSpec {
    #[must_use]
    pub const fn unsupported(
        name: &'static str,
        category: FeatureCategory,
        parity_sensitive: bool,
        reason: &'static str,
        docs: &'static str,
    ) -> Self {
        Self {
            name,
            category,
            status: SupportStatus::Unsupported { reason },
            parity_sensitive,
            docs,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[error("unsupported feature {feature}: {reason}")]
pub struct UnsupportedFeatureError {
    pub feature: &'static str,
    pub reason: &'static str,
}

impl UnsupportedFeatureError {
    #[must_use]
    pub const fn from_spec(feature: &'static FeatureSpec) -> Self {
        let reason = match feature.status {
            SupportStatus::Unsupported { reason } => reason,
            _ => "feature is not available in this build",
        };
        Self {
            feature: feature.name,
            reason,
        }
    }
}

pub const SMILES_PARSE_FEATURE: FeatureSpec = FeatureSpec {
    name: "smiles.parse",
    category: FeatureCategory::Io,
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "Parse SMILES into Molecule with registered sanitize, valence, aromaticity, kekulize, ring, hydrogen, stereo, and supported CX postprocessing. The documented parse scope is checked against pinned RDKit graph-state golden data; unmodeled branches return explicit errors instead of guessed chemistry.",
};

pub const SMILES_WRITE_FEATURE: FeatureSpec = FeatureSpec {
    name: "smiles.write",
    category: FeatureCategory::Io,
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "Write canonical, noncanonical, rooted, isomeric, kekule, explicit-bond, explicit-hydrogen, dative, atom-map, and supported CX SMILES branches. The documented branch matrix is checked exactly against pinned RDKit output; unsupported branches return explicit errors.",
};

pub const MOLBLOCK_IO_FEATURE: FeatureSpec = FeatureSpec {
    name: "molblock.io",
    category: FeatureCategory::Io,
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "Read and write the documented V2000/V3000 MolBlock and SDF branches with sanitize, remove-H, strict-parsing, coordinate-dimension, stereo, SGroup, RGroup, alias, value-line, and aromatic-bond handling. Covered fields and output branches are checked against pinned RDKit; unsupported extensions fail explicitly.",
};

pub const MOL2_READ_FEATURE: FeatureSpec = FeatureSpec {
    name: "mol2.read",
    category: FeatureCategory::Io,
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "Read the exposed Tripos MOL2 profile with sanitize, remove-H, CORINA-variant, and cleanup-substructure controls. Topology, atom and bond fields, coordinates, chirality, and SMILES output are checked against pinned RDKit.",
};

pub const HYDROGENS_FEATURE: FeatureSpec = FeatureSpec {
    name: "molecule.hydrogens",
    category: FeatureCategory::TopologyOperation,
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "Value-style and explicit in-place AddHs/RemoveHs operations with parameterized removal, atom and coordinate remapping, isotope tracking, stereo maintenance, and operation-contract checks. Modeled branches are covered by graph-state parity and focused source tests; unsupported states fail explicitly.",
};

pub const COORDINATE_2D_FEATURE: FeatureSpec = FeatureSpec {
    name: "coordinates.2d",
    category: FeatureCategory::TopologyOperation,
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "Generate 2D coordinates through registered value-style and in-place operations, including the local RDKit coordinate path, templates, constrained depiction, normalization, and straightening used by drawing and MolBlock workflows. Prepared coordinates and final depiction outputs are checked against pinned RDKit; unavailable CoordGen runtime branches fail explicitly.",
};

pub const COORDINATE_EDIT_FEATURE: FeatureSpec = FeatureSpec {
    name: "coordinates.edit",
    category: FeatureCategory::TopologyOperation,
    status: SupportStatus::Supported,
    parity_sensitive: true,
    docs: "Edit 2D and 3D coordinate blocks through registered value-style and explicit in-place operations, including replacement, clearing, and appending conformers with operation-contract validation.",
};

pub const CONFORMER_GENERATION_FEATURE: FeatureSpec = FeatureSpec {
    name: "coordinates.3d.conformer_generation",
    category: FeatureCategory::TopologyOperation,
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "Generate one or multiple 3D conformers through the exposed DG/KDG/ETDG/ETKDG presets and parameter controls, including explicit seeds, pruning, coordMap, CPCI, custom bounds matrices, stereo checks, macrocycles, and small-ring torsions. Deterministic covered branches are compared with pinned RDKit coordinates and status fields.",
};

pub const SANITIZE_FEATURE: FeatureSpec = FeatureSpec {
    name: "molecule.sanitize",
    category: FeatureCategory::TopologyOperation,
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "Run modeled RDKit sanitize stages as registered weak topology-state operations, including cleanup, properties, rings, kekulization, radicals, aromaticity, conjugation, hybridization, chirality cleanup, and hydrogen adjustment. End-to-end graph parity and focused stage tests cover the supported branches; unsupported requested states fail explicitly.",
};

pub const KEKULIZE_FEATURE: FeatureSpec = FeatureSpec {
    name: "molecule.with_kekulized_bonds",
    category: FeatureCategory::TopologyOperation,
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "Rewrite modeled aromatic systems through registered value-style and in-place kekulization, including fragment filtering, fused-ring candidate selection, backtracking, dummy-query permutations, clear-aromatic-flags behavior, and restoration on failure. Covered output branches are checked against pinned RDKit.",
};

pub const FINGERPRINT_FEATURE: FeatureSpec = FeatureSpec {
    name: "fingerprint.morgan",
    category: FeatureCategory::Fingerprint,
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "Provide the enumerated Morgan sparse-count, sparse-bit, hashed-count, explicit-bit, AdditionalOutput, and MACCS raw/public projection branches with exact pinned-RDKit parity. RDKFingerprint/topological, Avalon, and unmodeled preparation branches return structured unsupported errors rather than approximate vectors.",
};

pub const DESCRIPTORS_FEATURE: FeatureSpec = FeatureSpec {
    name: "descriptors.molecular",
    category: FeatureCategory::Core,
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "RDKit-aligned molecular descriptors exposed through the Rust API and Python bindings. The supported surface covers average and exact molecular weight (including only-heavy mode), formula options, H-bond donor/acceptor counts, fraction Csp3, every exposed Crippen include-H/force combination, every exposed TPSA force/S/P combination, aromatic-ring count, rotatable-bond modes, and QED. Supported corpus rows are compared field-by-field and bit-for-bit against pinned RDKit golden data, while inputs rejected by RDKit are required to fail in COSMolKit as well.",
};

pub const SUBSTRUCTURE_FEATURE: FeatureSpec = FeatureSpec {
    name: "substructure.match",
    category: FeatureCategory::Core,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Substructure matching is unfinished until the exposed molecule-query and future SMARTS-query surfaces pass strict RDKit parity tests. The current VF2 implementation must not be presented as RDKit-compatible while any atom/bond compatibility or query-matching branch remains approximate or marker-open.",
};

pub const DRAWING_FEATURE: FeatureSpec = FeatureSpec {
    name: "drawing.depiction",
    category: FeatureCategory::Drawing,
    status: SupportStatus::Supported,
    parity_sensitive: true,
    docs: "Render the documented molecule depiction surface as SVG and PNG. Covered SVG output is compared exactly with pinned RDKit after normalizing only tool identifiers; PNG is deterministic local rasterization of that SVG through usvg and resvg, not a Cairo or Qt byte-parity claim. Unsupported annotation branches fail explicitly.",
};

pub const STEREO_FEATURE: FeatureSpec = FeatureSpec {
    name: "stereo.perception",
    category: FeatureCategory::Stereo,
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "Inspect modeled atom and bond stereo state and assign atom ChiralTag values from a selected 3D conformer through the registered with_chiral_tags_from_structure operation. This stable 3D assignment surface has exact full-state parity with pinned RDKit 2026.03.1 assignChiralTypesFrom3D across all 77 fixed oracle records, covering default or explicit conformer selection, replacement control, tetrahedral C/S/Se centers, environment-enabled square-planar/trigonal-bipyramidal/octahedral centers, property updates, no-op paths, and source-defined errors. It preserves topology and coordinates and commits no partial state on error. The broader assignStereochemistryFrom3D workflow, 3D double-bond direction/E-Z assignment, CIP orchestration, and distinct-substituent validation are separate capabilities and are not claimed by this feature.",
};

pub const VALENCE_FEATURE: FeatureSpec = FeatureSpec {
    name: "valence.assignment",
    category: FeatureCategory::Valence,
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "Assign modeled RDKit valence, implicit hydrogens, property-cache state, and radicals through registered operations used by parsing and sanitization. Covered graph fields and focused source branches match pinned RDKit; unmodeled query or dative states fail explicitly.",
};

pub const RINGS_FEATURE: FeatureSpec = FeatureSpec {
    name: "rings.symm_sssr",
    category: FeatureCategory::TopologyOperation,
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "Provide modeled SSSR, symmetrized SSSR, fast ring traversal, and URF-backed ring-family and relevant-cycle perception. Ring membership is covered by pinned-RDKit graph parity and focused source tests cover the supported SSSR and ring-family graph states.",
};

pub const AROMATICITY_FEATURE: FeatureSpec = FeatureSpec {
    name: "aromaticity.assignment",
    category: FeatureCategory::TopologyOperation,
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "Assign modeled RDKit aromatic atom and bond state through registered operations. Covered graph states match pinned RDKit and unsupported aromaticity models or source states fail explicitly.",
};

pub const INCHI_FEATURE: FeatureSpec = FeatureSpec {
    name: "inchi.scalar",
    category: FeatureCategory::Io,
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "The four audited scalar APIs (MolToInchi, MolToInchiKey, InchiToInchiKey, and MolFromInchi) reproduce pinned RDKit 2026.03.1 and official InChI v1.07.5 behavior exactly where the official C source defines behavior. The NormalizeAndCompare initial-buffer allocation-failure path executes undefined C behavior, so COSMolKit returns a deterministic structured allocation error. Other official InChI API families are outside this public feature.",
};

pub const BATCH_FEATURE: FeatureSpec = FeatureSpec {
    name: "batch.operations",
    category: FeatureCategory::Batch,
    status: SupportStatus::Supported,
    parity_sensitive: false,
    docs: "Stable batch construction from SMILES lists, ordered transformations via registered molecule operations, \
           error modes (Strict/KeepErrors), valid mask, filter valid, SMILES export with params, \
           PNG image export, configurable progress reporting, and deterministic ordered parallel \
           execution through Rayon with batch-level or per-call worker-count selection.",
};

pub const BIO_STRUCTURE_FEATURE: FeatureSpec = FeatureSpec {
    name: "bio.structure",
    category: FeatureCategory::BioHierarchy,
    status: SupportStatus::Supported,
    parity_sensitive: false,
    docs: "Flat-row BioStructure hierarchy and coordinate storage for protein, PDB, mmCIF, mmJSON, and chem-comp workflows. Public access is read-only; mutation goes through crate-internal builders or registered BioStructure operations.",
};

pub const BIO_PDB_COORDINATE_SUBSET_READ_FEATURE: FeatureSpec = FeatureSpec {
    name: "bio.pdb.coordinate_subset.read",
    category: FeatureCategory::Io,
    status: SupportStatus::Supported,
    parity_sensitive: true,
    docs: "Read the documented Gemmi-aligned PDB structural record scope into BioStructure, including coordinates, models, anisotropic data, hierarchy, entities, connections, selected metadata, crystallography, and NCS transforms. Unmodeled records fail explicitly in strict mode.",
};

pub const BIO_MMCIF_ATOM_SITE_SUBSET_READ_FEATURE: FeatureSpec = FeatureSpec {
    name: "bio.mmcif.atom_site_subset.read",
    category: FeatureCategory::Io,
    status: SupportStatus::Supported,
    parity_sensitive: true,
    docs: "Read the documented Gemmi-aligned mmCIF, mmJSON, and chem-comp structural scope into BioStructure, including atom sites, entities and sequences, references, connectivity, assemblies, SIFTS mappings, crystallography, and NCS transforms. Unmodeled structural categories fail explicitly where required.",
};

pub const PUBLIC_FEATURES: &[&FeatureSpec] = &[
    &SMILES_PARSE_FEATURE,
    &SMILES_WRITE_FEATURE,
    &MOLBLOCK_IO_FEATURE,
    &MOL2_READ_FEATURE,
    &HYDROGENS_FEATURE,
    &COORDINATE_2D_FEATURE,
    &COORDINATE_EDIT_FEATURE,
    &CONFORMER_GENERATION_FEATURE,
    &SANITIZE_FEATURE,
    &KEKULIZE_FEATURE,
    &FINGERPRINT_FEATURE,
    &DESCRIPTORS_FEATURE,
    &SUBSTRUCTURE_FEATURE,
    &DRAWING_FEATURE,
    &STEREO_FEATURE,
    &VALENCE_FEATURE,
    &RINGS_FEATURE,
    &AROMATICITY_FEATURE,
    &INCHI_FEATURE,
    &BATCH_FEATURE,
    &DG_BOUNDS_FEATURE,
    &BIO_STRUCTURE_FEATURE,
    &BIO_PDB_COORDINATE_SUBSET_READ_FEATURE,
    &BIO_MMCIF_ATOM_SITE_SUBSET_READ_FEATURE,
    &BIO_SELECTION_FEATURE,
];

pub const DG_BOUNDS_FEATURE: FeatureSpec = FeatureSpec {
    name: "distgeom.bounds_matrix",
    category: FeatureCategory::Core,
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "Generate the documented RDKit-style distance-geometry bounds matrix, including triangle smoothing, topological 1-2 through 1-5 bounds, ring and macrocycle handling, and VDW lower bounds. Matrix shape and every covered lower and upper entry are compared with pinned RDKit.",
};

pub const BIO_SELECTION_FEATURE: FeatureSpec = FeatureSpec {
    name: "bio.selection",
    category: FeatureCategory::BioSelection,
    status: SupportStatus::Experimental,
    parity_sensitive: false,
    docs: "Experimental BioStructure selection and filtering operations (e.g. remove_waters).",
};
