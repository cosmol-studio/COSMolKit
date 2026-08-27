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
    docs: "Parse SMILES into Molecule with registered sanitize, valence, aromaticity, kekulize, ring, hydrogen, stereo, and supported CX postprocessing. The documented parse scope is checked against pinned RDKit graph-state golden data; separately scoped upstream extensions are not advertised as part of this boundary.",
};

pub const SMILES_WRITE_FEATURE: FeatureSpec = FeatureSpec {
    name: "smiles.write",
    category: FeatureCategory::Io,
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "Write canonical, noncanonical, rooted, isomeric, kekule, explicit-bond, explicit-hydrogen, dative, atom-map, and supported CX SMILES branches. The documented branch matrix is checked exactly against pinned RDKit output; other upstream writer extensions are separate capability boundaries.",
};

pub const MOLBLOCK_IO_FEATURE: FeatureSpec = FeatureSpec {
    name: "molblock.io",
    category: FeatureCategory::Io,
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "Read and write the documented V2000/V3000 MolBlock and SDF branches with sanitize, remove-H, strict-parsing, coordinate-dimension, stereo, SGroup, RGroup, alias, value-line, and aromatic-bond handling. Covered fields and output branches are checked against pinned RDKit; unrelated CTAB extensions are outside this named boundary.",
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
    docs: "Value-style and explicit in-place AddHs/RemoveHs operations with parameterized removal, atom and coordinate remapping, isotope tracking, stereo maintenance, and operation-contract checks. The documented state model is covered by graph-state parity and focused source tests.",
};

pub const COORDINATE_2D_FEATURE: FeatureSpec = FeatureSpec {
    name: "coordinates.2d",
    category: FeatureCategory::TopologyOperation,
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "Generate 2D coordinates through registered value-style and in-place operations, including the local RDKit coordinate path, templates, constrained depiction, normalization, and straightening used by drawing and MolBlock workflows. Prepared coordinates and final depiction outputs are checked against pinned RDKit; the separate CoordGen runtime is outside this local-coordinate boundary.",
};

pub const COORDINATE_EDIT_FEATURE: FeatureSpec = FeatureSpec {
    name: "coordinates.edit",
    category: FeatureCategory::TopologyOperation,
    status: SupportStatus::Supported,
    parity_sensitive: true,
    docs: "Edit 2D and 3D coordinate blocks through registered value-style and explicit in-place operations, including replacement, clearing, and appending conformers with operation-contract validation.",
};

pub const MOLALIGN_FEATURE: FeatureSpec = FeatureSpec {
    name: "coordinates.molalign",
    category: FeatureCategory::TopologyOperation,
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "Measure and align 3D conformers through the source-backed ordinary RDKit MolAlign boundary, including explicit and automatic atom maps, weighted and reflected alignment, best-map selection, coordinate-frame RMSD, and all-conformer ordering. Measurement APIs are read-only; coordinate changes use registered value-style or explicitly named in-place operations. O3A is a separate capability.",
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
    docs: "Run modeled RDKit sanitize stages as registered weak topology-state operations, including cleanup, properties, rings, kekulization, radicals, aromaticity, conjugation, hybridization, chirality cleanup, and hydrogen adjustment. End-to-end graph parity and focused stage tests cover the documented state and option boundary.",
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
    docs: "Provide the enumerated Morgan sparse-count, sparse-bit, hashed-count, explicit-bit, AdditionalOutput, and MACCS raw/public projection branches with exact pinned-RDKit parity. Other fingerprint families have separate feature specifications and do not silently fall back to Morgan.",
};

pub const ATOM_PAIR_FINGERPRINT_FEATURE: FeatureSpec = FeatureSpec {
    name: "fingerprint.atom_pair",
    category: FeatureCategory::Fingerprint,
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "Provide the source-backed AtomPair sparse-count, folded-count, sparse-bit, explicit-bit, AdditionalOutput, 2D/3D distance, chirality, selector, custom-invariant, count-simulation, metadata/JSON, and ordered batch branches through one generator core. Exact validation covers every mutually parseable ChEMBL 37 record across ten profiles and 118,809,964 comparisons against pinned RDKit.",
};

pub const PATTERN_FINGERPRINT_FEATURE: FeatureSpec = FeatureSpec {
    name: "fingerprint.pattern",
    category: FeatureCategory::Fingerprint,
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "Provide the source-backed ordinary-molecule Pattern fingerprint with the fixed 13-query table, ordinary and tautomer-aware hashing, query-bearing molecule behavior, variable nonzero widths, and ordered batch execution through one compile-once core. Focused, small, maintained 5,000-row, and complete ChEMBL 37 validation covers ten full-corpus profiles and 28,978,040 exact comparisons against pinned RDKit. RDKit identifies Pattern version 1.0.0 as experimental; the separate MolBundle intersection overload remains outside COSMolKit's molecule and ordered-batch model.",
};

pub const TOPOLOGICAL_TORSION_FINGERPRINT_FEATURE: FeatureSpec = FeatureSpec {
    name: "fingerprint.topological_torsion",
    category: FeatureCategory::Fingerprint,
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "Provide the modern Topological Torsion sparse-count, sparse-bit, folded-count, explicit-bit, AdditionalOutput, JSON, scalar, ordered bulk, and Rust batch branches plus the three legacy compatibility adapters. Focused and maintained 5,000-row validation covers modern and legacy options, selections, chirality, path behavior, provenance, collisions, and errors. The complete ChEMBL 37 audit matches every mutually parseable record across nine profiles and 127,503,376 exact vector/provenance comparisons against pinned RDKit. RDKFingerprint, Atom Pair fingerprints, and unrelated fingerprint families remain separate capabilities.",
};

pub const AVALON_FINGERPRINT_FEATURE: FeatureSpec = FeatureSpec {
    name: "fingerprint.avalon",
    category: FeatureCategory::Fingerprint,
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "Provide the source-backed Avalon/REACCS explicit-bit path for n_bits, is_query, and bit_flags, including source byte rounding, aromaticity passes, query branches, and non-SSS feature families. The pinned RDKit adapter matches exactly across the maintained 5,000-row, 23-profile validation matrix; count/string overloads and unrelated Avalon APIs remain out of scope.",
};

pub const LAYERED_FINGERPRINT_FEATURE: FeatureSpec = FeatureSpec {
    name: "fingerprint.layered",
    category: FeatureCategory::Fingerprint,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Provide the single source-backed RDKit Layered fingerprint 0.7.0 core with all six active layers, rooted linear and branched path selection, masks, seeded atom counts, scalar calls, and ordered batch execution. All 2,897,804 mutually parseable ChEMBL 37 records match pinned RDKit 2026.03.1 across 18 profiles and 52,160,472 exact comparisons. RDKit labels this family experimental, so COSMolKit preserves that status. COSMolKit intentionally returns the documented bond-path result for unrooted linear calls instead of reproducing the pinned source's input-dependent process crash.",
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
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "Match ordinary-molecule query-bearing Molecule values through the canonical VF2 matcher, including SMARTS parser/writer round trips, recursive query execution, atom/bond predicates, stereo, callbacks, match ordering, counts, errors, and parameters covered by the pinned RDKit SMARTS and substructure corpora. Reaction SMARTS and database/container SMARTS remain separate unsupported capability boundaries.",
};

pub const TAUTOMER_ENUMERATION_FEATURE: FeatureSpec = FeatureSpec {
    name: "tautomer.enumeration",
    category: FeatureCategory::Core,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Enumerate source-ordered molecular tautomers through the registered multi-molecule operation path. The implementation is under source-level RDKit parity validation and is not promoted to a parity-supported status until the complete tautomer plan finishes.",
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

pub const CIP_LABELER_FEATURE: FeatureSpec = FeatureSpec {
    name: "stereo.cip_labeler",
    category: FeatureCategory::Stereo,
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "Assign modern molecular-context CIP descriptors through the source-backed RDKit 2026.03.1 CIPLabeler path. Full-molecule and selected atom/bond assignment share one implementation and expose typed R/S/r/s, E/Z, and M/P/m/p descriptors. Focused, maintained 5,000-row, and complete ChEMBL 37 full/selected state comparisons pass exactly. Non-tetrahedral configurations outside the pinned findConfigs dispatcher and source-equivalent cancellation remain separate unsupported boundaries.",
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

pub const BIO_PDB_READ_FEATURE: FeatureSpec = FeatureSpec {
    name: "bio.pdb.read",
    category: FeatureCategory::Io,
    status: SupportStatus::Supported,
    parity_sensitive: true,
    docs: "Read the documented Gemmi-aligned PDB structural record scope into BioStructure, including decimal and hybrid-36 identifiers, coordinates, models, anisotropic data, hierarchy, entities, connections, selected metadata, crystallography, and NCS transforms.",
};

pub const BIO_MMCIF_READ_FEATURE: FeatureSpec = FeatureSpec {
    name: "bio.mmcif.read",
    category: FeatureCategory::Io,
    status: SupportStatus::Supported,
    parity_sensitive: true,
    docs: "Read the documented Gemmi-aligned mmCIF, mmJSON, and chem-comp structural scope into BioStructure, including atom sites, entities and sequences, references, connectivity, assemblies, SIFTS mappings, crystallography, and NCS transforms. Unmodeled structural categories fail explicitly where required.",
};

pub const BIO_MMCIF_WRITE_FEATURE: FeatureSpec = FeatureSpec {
    name: "bio.mmcif.write",
    category: FeatureCategory::Io,
    status: SupportStatus::Supported,
    parity_sensitive: true,
    docs: "Serialize BioStructure as canonical Gemmi-aligned mmCIF with source-defined category groups and CIF formatting options. The writer emits state represented by BioStructure and belongs only to BioStructure, not Protein or Molecule.",
};

#[deprecated(note = "use BIO_PDB_READ_FEATURE")]
pub const BIO_PDB_COORDINATE_SUBSET_READ_FEATURE: FeatureSpec = BIO_PDB_READ_FEATURE;

#[deprecated(note = "use BIO_MMCIF_READ_FEATURE")]
pub const BIO_MMCIF_ATOM_SITE_SUBSET_READ_FEATURE: FeatureSpec = BIO_MMCIF_READ_FEATURE;

pub const PUBLIC_FEATURES: &[&FeatureSpec] = &[
    &SMILES_PARSE_FEATURE,
    &SMILES_WRITE_FEATURE,
    &MOLBLOCK_IO_FEATURE,
    &MOL2_READ_FEATURE,
    &HYDROGENS_FEATURE,
    &COORDINATE_2D_FEATURE,
    &COORDINATE_EDIT_FEATURE,
    &MOLALIGN_FEATURE,
    &CONFORMER_GENERATION_FEATURE,
    &SANITIZE_FEATURE,
    &KEKULIZE_FEATURE,
    &FINGERPRINT_FEATURE,
    &ATOM_PAIR_FINGERPRINT_FEATURE,
    &PATTERN_FINGERPRINT_FEATURE,
    &TOPOLOGICAL_TORSION_FINGERPRINT_FEATURE,
    &AVALON_FINGERPRINT_FEATURE,
    &LAYERED_FINGERPRINT_FEATURE,
    &DESCRIPTORS_FEATURE,
    &SUBSTRUCTURE_FEATURE,
    &TAUTOMER_ENUMERATION_FEATURE,
    &DRAWING_FEATURE,
    &STEREO_FEATURE,
    &CIP_LABELER_FEATURE,
    &VALENCE_FEATURE,
    &RINGS_FEATURE,
    &AROMATICITY_FEATURE,
    &INCHI_FEATURE,
    &BATCH_FEATURE,
    &DG_BOUNDS_FEATURE,
    &BIO_STRUCTURE_FEATURE,
    &BIO_PDB_READ_FEATURE,
    &BIO_MMCIF_READ_FEATURE,
    &BIO_MMCIF_WRITE_FEATURE,
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

#[cfg(test)]
mod tests {
    use super::{
        ATOM_PAIR_FINGERPRINT_FEATURE, PATTERN_FINGERPRINT_FEATURE, PUBLIC_FEATURES, SupportStatus,
    };

    #[test]
    fn atom_pair_support_metadata_records_the_validated_rdkit_boundary() {
        assert_eq!(
            ATOM_PAIR_FINGERPRINT_FEATURE.status,
            SupportStatus::SupportedWithRdkitParity {
                rdkit_version: "2026.03.1",
            }
        );
        assert!(
            PUBLIC_FEATURES
                .iter()
                .any(|feature| **feature == ATOM_PAIR_FINGERPRINT_FEATURE)
        );
    }

    #[test]
    fn pattern_support_metadata_records_the_validated_rdkit_boundary() {
        assert_eq!(
            PATTERN_FINGERPRINT_FEATURE.status,
            SupportStatus::SupportedWithRdkitParity {
                rdkit_version: "2026.03.1",
            }
        );
        assert!(
            PUBLIC_FEATURES
                .iter()
                .any(|feature| **feature == PATTERN_FINGERPRINT_FEATURE)
        );
    }
}
