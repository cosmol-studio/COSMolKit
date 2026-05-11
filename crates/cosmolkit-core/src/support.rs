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

pub const SMILES_PARSE_FEATURE: FeatureSpec = FeatureSpec::unsupported(
    "smiles.parse",
    FeatureCategory::Io,
    true,
    "SMILES parser has not been ported to the redesigned core",
    "Parse SMILES into Molecule.",
);

pub const SMILES_WRITE_FEATURE: FeatureSpec = FeatureSpec::unsupported(
    "smiles.write",
    FeatureCategory::Io,
    true,
    "SMILES writer has not been ported to the redesigned core",
    "Serialize Molecule as SMILES.",
);

pub const MOLBLOCK_IO_FEATURE: FeatureSpec = FeatureSpec::unsupported(
    "molblock.io",
    FeatureCategory::Io,
    true,
    "MolBlock/SDF I/O has not been ported to the redesigned core",
    "Read and write MolBlock/SDF records.",
);

pub const HYDROGENS_FEATURE: FeatureSpec = FeatureSpec::unsupported(
    "molecule.hydrogens",
    FeatureCategory::TopologyOperation,
    true,
    "hydrogen addition/removal has not been ported to the redesigned core",
    "Add or remove explicit hydrogens as value-style topology operations.",
);

pub const COORDINATE_2D_FEATURE: FeatureSpec = FeatureSpec::unsupported(
    "coordinates.2d",
    FeatureCategory::TopologyOperation,
    true,
    "2D coordinate generation has not been ported to the redesigned core",
    "Generate 2D coordinates while preserving molecule value semantics.",
);

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

pub const FINGERPRINT_FEATURE: FeatureSpec = FeatureSpec::unsupported(
    "fingerprint.morgan",
    FeatureCategory::Fingerprint,
    true,
    "Morgan fingerprint generation has not been ported to the redesigned core",
    "Compute Morgan fingerprints and optional additional output.",
);

pub const DRAWING_FEATURE: FeatureSpec = FeatureSpec::unsupported(
    "drawing.depiction",
    FeatureCategory::Drawing,
    true,
    "drawing and depiction have not been ported to the redesigned core",
    "Prepare and render molecular depictions.",
);

pub const STEREO_FEATURE: FeatureSpec = FeatureSpec::unsupported(
    "stereo.perception",
    FeatureCategory::Stereo,
    true,
    "stereochemistry perception has not been ported to the redesigned core",
    "Perceive and validate stereochemistry.",
);

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

pub const BATCH_FEATURE: FeatureSpec = FeatureSpec::unsupported(
    "batch.operations",
    FeatureCategory::Batch,
    false,
    "batch operation execution has not been ported to the redesigned core",
    "Run ordered, traceable batch operations.",
);

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
    &BIO_STRUCTURE_FEATURE,
    &BIO_PDB_COORDINATE_SUBSET_READ_FEATURE,
    &BIO_MMCIF_ATOM_SITE_SUBSET_READ_FEATURE,
    &BIO_SELECTION_FEATURE,
];

pub const BIO_SELECTION_FEATURE: FeatureSpec = FeatureSpec {
    name: "bio.selection",
    category: FeatureCategory::BioSelection,
    status: SupportStatus::Experimental,
    parity_sensitive: false,
    docs: "Experimental BioStructure selection and filtering operations (e.g. remove_waters).",
};
