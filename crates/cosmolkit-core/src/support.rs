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
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Parse SMILES into Molecule with sanitize integration through registered operations (kekulize, valence, aromaticity). CX extensions (coords, labels, values, props, radicals, stereo, SGroups) are parsed. Remove-H with updateExplicitCount is implemented.",
};

pub const SMILES_WRITE_FEATURE: FeatureSpec = FeatureSpec {
    name: "smiles.write",
    category: FeatureCategory::Io,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Plain SMILES output (canonical and noncanonical) is implemented. Aromatic atoms (lowercase) and BondOrder::Aromatic bonds are supported. CX extensions are implemented: atom labels, molfile values, 2D coordinates, radicals, atom properties, enhanced stereo groups, SGroups, coordinate/hydrogen/zero bonds. Fragment API, random SMILES, and stereochemical perception branch remain unported.",
};

pub const MOLBLOCK_IO_FEATURE: FeatureSpec = FeatureSpec {
    name: "molblock.io",
    category: FeatureCategory::Io,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Experimental V2000/V3000 MolBlock/SDF writer with parity flag, bond-stereo, SGroup, RGroup, alias, value lines, and aromatic-bond bookkeeping. Reader has partial V2000 parsing. Unsupported branches (complex SMARTS queries, atropisomer wedge-bonds) fail closed.",
};

pub const HYDROGENS_FEATURE: FeatureSpec = FeatureSpec {
    name: "molecule.hydrogens",
    category: FeatureCategory::TopologyOperation,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Experimental value-style explicit hydrogen operations. Remove-H is being ported through the operation-contract path; unsupported source branches fail closed.",
};

pub const COORDINATE_2D_FEATURE: FeatureSpec = FeatureSpec {
    name: "coordinates.2d",
    category: FeatureCategory::TopologyOperation,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Generate 2D coordinates while preserving molecule value semantics. Invalid inputs and unported depiction branches fail explicitly.",
};

pub const SANITIZE_FEATURE: FeatureSpec = FeatureSpec {
    name: "molecule.sanitize",
    category: FeatureCategory::TopologyOperation,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Run supported RDKit-aligned sanitization steps as a weak topology-state operation; unported requested steps fail closed.",
};

pub const KEKULIZE_FEATURE: FeatureSpec = FeatureSpec {
    name: "molecule.with_kekulized_bonds",
    category: FeatureCategory::TopologyOperation,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Experimental operation-pipeline skeleton for kekulized bond rewriting.",
};

pub const FINGERPRINT_FEATURE: FeatureSpec = FeatureSpec {
    name: "fingerprint.morgan",
    category: FeatureCategory::Fingerprint,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Compute Morgan-style fingerprints with connectivity invariants (RDKit component-vector hash style). Environment propagation uses RDKit's seed=layer + sorted neighbor-pair hashing. Chirality support, feature invariants (element/property classification), custom atom/bond invariants, count-simulation with configurable bounds. Hash-value alignment is structurally compatible but not bit-identical (uses own hash_combine instead of gboost::hash).",
};

pub const DRAWING_FEATURE: FeatureSpec = FeatureSpec {
    name: "drawing.depiction",
    category: FeatureCategory::Drawing,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "SVG/PNG molecule renderer ported from RDKit MolDraw2D. \
           Includes atom labels (isotope/charge/H/map), bond geometry \
           (single/double/triple/wedge/aromatic), radical dots, clash \
           detection, scale calculation. SVG output via native XML; \
           PNG via usvg+resvg rasterization.",
};

pub const STEREO_FEATURE: FeatureSpec = FeatureSpec {
    name: "stereo.perception",
    category: FeatureCategory::Stereo,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Tetrahedral stereo detection from typed state (ChiralTag + chiral_permutation). Double-bond E/Z potential detection. Full geometric perception from 3D coordinates and CIP labeler not yet ported.",
};

pub const VALENCE_FEATURE: FeatureSpec = FeatureSpec {
    name: "valence.assignment",
    category: FeatureCategory::Valence,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Experimental RDKit-aligned valence and implicit hydrogen assignment.",
};

pub const RINGS_FEATURE: FeatureSpec = FeatureSpec {
    name: "rings.symm_sssr",
    category: FeatureCategory::TopologyOperation,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Experimental RDKit-aligned SSSR and symmetrized SSSR ring perception.",
};

pub const AROMATICITY_FEATURE: FeatureSpec = FeatureSpec {
    name: "aromaticity.assignment",
    category: FeatureCategory::TopologyOperation,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Experimental RDKit-aligned aromaticity assignment scaffold with fail-closed unsupported branches.",
};

pub const BATCH_FEATURE: FeatureSpec = FeatureSpec {
    name: "batch.operations",
    category: FeatureCategory::Batch,
    status: SupportStatus::Experimental,
    parity_sensitive: false,
    docs: "Batch construction from SMILES list, ordered transformations via registered molecule operations, \
           error modes (Strict/KeepErrors), valid mask, filter valid, SMILES export with params, \
           and PNG image export. Batch scheduling and parallel execution are not yet implemented.",
};

pub const BIO_STRUCTURE_FEATURE: FeatureSpec = FeatureSpec {
    name: "bio.structure",
    category: FeatureCategory::BioHierarchy,
    status: SupportStatus::Experimental,
    parity_sensitive: false,
    docs: "Experimental flat-row BioStructure hierarchy and coordinate storage. Public access is read-only; mutation must go through crate-internal builders or registered BioStructure operations.",
};

pub const BIO_PDB_COORDINATE_SUBSET_READ_FEATURE: FeatureSpec = FeatureSpec {
    name: "bio.pdb.coordinate_subset.read",
    category: FeatureCategory::Io,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Experimental Gemmi-aligned PDB coordinate subset reader. It models ATOM/HETATM, MODEL/ENDMDL, ANISOU, residue identity, atom identity, element, charge, occupancy, B factor, and coordinates. PDB metadata, connectivity, sequence, TER/entity semantics, secondary-structure, and other unported records are not complete PDB support.",
};

pub const BIO_MMCIF_ATOM_SITE_SUBSET_READ_FEATURE: FeatureSpec = FeatureSpec {
    name: "bio.mmcif.atom_site_subset.read",
    category: FeatureCategory::Io,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Experimental Gemmi-aligned mmCIF _atom_site subset reader. It is not complete mmCIF structure semantics; unsupported source branches remain visible in io::bio.",
};

pub const PUBLIC_FEATURES: &[&FeatureSpec] = &[
    &SMILES_PARSE_FEATURE,
    &SMILES_WRITE_FEATURE,
    &MOLBLOCK_IO_FEATURE,
    &HYDROGENS_FEATURE,
    &COORDINATE_2D_FEATURE,
    &SANITIZE_FEATURE,
    &KEKULIZE_FEATURE,
    &FINGERPRINT_FEATURE,
    &DRAWING_FEATURE,
    &STEREO_FEATURE,
    &VALENCE_FEATURE,
    &RINGS_FEATURE,
    &AROMATICITY_FEATURE,
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
    status: SupportStatus::Experimental,
    parity_sensitive: false,
    docs: "Distance geometry bounds matrix generation. Uses topology-based 1-2/1-3/1-4 bounds \
           with hybridization angle estimates, VDW lower bounds, and triangle inequality smoothing. \
           Does not include UFF bond-length parameters or full RDKit triangle smoother. \
           The bounds matrix is an n×n Vec<Vec<f64>> of upper bounds.",
};

pub const BIO_SELECTION_FEATURE: FeatureSpec = FeatureSpec {
    name: "bio.selection",
    category: FeatureCategory::BioSelection,
    status: SupportStatus::Experimental,
    parity_sensitive: false,
    docs: "Experimental BioStructure selection and filtering operations (e.g. remove_waters).",
};
