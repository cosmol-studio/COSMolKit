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
    docs: "Parse SMILES into Molecule with sanitize integration through registered operations (kekulize, valence, aromaticity, rings). RDKit-aligned postprocessing includes first-2D/first-3D conformer selection, wedged/3D stereo assignment (including non-tetrahedral branches), atropisomer chirality mutation paths, CX wiggly-bond direction cleanup, and _NeedsQueryScan ring/non-ring query completion. CX extensions (coords, labels, values, props, radicals, stereo, SGroups, hierarchy, polymer, linknodes) are parsed. Remove-H isotope tracking and the targeted fixture-backed reader parity gaps from the current checklist are closed, but the reader is not marker-complete: `notation/smiles.rs` still contains 1 `RDKit❌❌`, 2 `RDKit❗❗`, 14 `RDKit✔️❌`, and 713 `RDKit❗✔️` copied-source lines across the remaining parser/helper blockers tracked by the gap report. Remaining unported or unresolved branches fail closed or remain explicitly tracked by gap reports.",
};

pub const SMILES_WRITE_FEATURE: FeatureSpec = FeatureSpec {
    name: "smiles.write",
    category: FeatureCategory::Io,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Plain SMILES output (canonical and noncanonical) is implemented, including the checklist-closed parity cases for noncanonical/rooted/connected/ring/fused/CIP-tie double-bond direction output and non-tetrahedral class emission/permutation recomputation. Aromatic atoms (lowercase) and BondOrder::Aromatic bonds are supported. CX writer blocks are implemented for bond wedge/dash config, ring-bond cis/trans config, linknodes, polymer SGroups, SGroup hierarchy, atropisomer bonds, atom labels, molfile values, 2D coordinates, radicals, atom properties, enhanced stereo groups, and coordinate/hydrogen/zero bonds. Writer behavior depends on the chemistry-core sanitize/valence/kekulize/ring state pipeline, and writer-internal unsupported stage guards were replaced by concrete invariant/validation errors where reachable. The frozen writer file is marker-closed for the current checklist scope, but the feature remains experimental and depends on broader parser/chemistry parity surfaces that are still open elsewhere.",
};

pub const MOLBLOCK_IO_FEATURE: FeatureSpec = FeatureSpec {
    name: "molblock.io",
    category: FeatureCategory::Io,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Experimental V2000/V3000 MolBlock/SDF writer with parity flag, bond-stereo, SGroup, RGroup, alias, value lines, and aromatic-bond bookkeeping. Reader has partial V2000 parsing. The writer and reader remain dependent on explicit valence/kekulize/ring state management. Unsupported branches (complex SMARTS queries, atropisomer wedge-bonds) fail closed.",
};

pub const HYDROGENS_FEATURE: FeatureSpec = FeatureSpec {
    name: "molecule.hydrogens",
    category: FeatureCategory::TopologyOperation,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Experimental value-style explicit hydrogen operations. Remove-H is being ported through the operation-contract path and depends on valence/kekulize/ring state being available; unsupported source branches fail closed.",
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
    docs: "Run supported RDKit-aligned sanitization steps as a weak topology-state operation, sequencing the explicit valence/kekulize/ring handoff used by the SMILES reader and other operations. Full RDKit flag/error/cleanup closure is still pending in the broader operation-orchestration surface: `operations/ops.rs` still contains 216 `RDKit✔️❌` copied-source lines across the remaining sanitize/property/cleanup orchestration blocks and helper routines tracked by the gap report. Unported requested steps fail closed.",
};

pub const KEKULIZE_FEATURE: FeatureSpec = FeatureSpec {
    name: "molecule.with_kekulized_bonds",
    category: FeatureCategory::TopologyOperation,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Experimental operation-pipeline for kekulized bond rewriting. This is the dependency used by fused aromatic assignment and KekulizeIfPossible restoration. Fragment filtering, fused aromatic candidate selection, worker ordering/backtracking, dummy-question permutation, and value-style `KekulizeIfPossible` restoration have focused regression coverage, but broader operation-state interaction closure is still pending and `chemistry/kekulize.rs` still contains 397 `RDKit✔️❌` copied-source lines in the current frozen-scope audit; unsupported branches fail closed.",
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
           (single/double/triple/wedge/aromatic/dative), radical dots, clash \
           detection, scale calculation, and smoothed bond joins. \
           Annotations: CIP codes (R/S, E/Z), atom notes, bond notes, \
           SGroup data, brackets, variable bonds, link nodes, close-contact \
           markers, and highlights. SVG metadata, data-tag attributes, and \
           CSS class output for atoms/bonds. \
           SVG output via native XML; PNG via usvg+resvg rasterization.",
};

pub const STEREO_FEATURE: FeatureSpec = FeatureSpec {
    name: "stereo.perception",
    category: FeatureCategory::Stereo,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Tetrahedral stereo detection from typed state (ChiralTag + chiral_permutation). \
           CIP ranking system (assignAtomCIPRanks with iterative neighbor-rank refinement) ported. \
           R/S label assignment (assignAtomChiralCodes) from ChiralTag + permutation. \
           Double-bond E/Z potential detection. Pseudo-3D wedge-based chiral tag detection \
           (atomChiralTypeFromBondDirPseudo3D). Full non-tetrahedral stereo infrastructure \
           (SquarePlanar, TrigonalBipyramidal, Octahedral swap tables and across-atom lookup). \
           Ring stereochemistry special-case detection. Full CIP-based bond stereo codes \
           and assignLegacyCIPLabels dispatcher ported. assignAtomChiralTagsFromStructure \
           (full 3D coordinate-based ChiralTag assignment) remains blocked on Conformer \
           infrastructure completeness.",
};

pub const VALENCE_FEATURE: FeatureSpec = FeatureSpec {
    name: "valence.assignment",
    category: FeatureCategory::Valence,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Experimental RDKit-aligned valence and implicit hydrogen assignment. This is a shared dependency for sanitize, kekulize, and SMILES postprocessing. `chemistry/valence.rs` now only retains 4 `RDKit✔️❌` copied-source lines in `ValenceContext::new`, and remaining work is concentrated in property-cache maintenance, radicals, dative/query edge cases, and broader entrypoint/orchestration logic in `operations/ops.rs`. Unsupported branches fail closed.",
};

pub const RINGS_FEATURE: FeatureSpec = FeatureSpec {
    name: "rings.symm_sssr",
    category: FeatureCategory::TopologyOperation,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Experimental RDKit-aligned SSSR, symmetrized SSSR, fast ring traversal, and URF-enabled ring-family/relevant-cycle perception via `cosmolkit_ringdecomposer`. SSSR active-bond filtering, D2 duplicate-candidate handling, D3/extra-ring discovery, symmetrized K4 storage, fastFindRings DFS traversal, and the URF-enabled ring-family/relevant-cycle path have focused regression coverage. The frozen ring-perception file is marker-closed for the current checklist scope, but the feature remains experimental and is not a blanket claim of complete RDKit ring parity outside that audited scope.",
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
    docs: "Experimental flat-row BioStructure hierarchy and coordinate storage. This is COSMolKit's single public structural model for protein/PDB/mmCIF work. Public access is read-only; mutation must go through crate-internal builders or registered BioStructure operations.",
};

pub const BIO_PDB_COORDINATE_SUBSET_READ_FEATURE: FeatureSpec = FeatureSpec {
    name: "bio.pdb.coordinate_subset.read",
    category: FeatureCategory::Io,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Experimental Gemmi-aligned PDB coordinate subset reader into BioStructure. This is the structural IO path and the required front end for future RDKit-compatible molecule input. The reader models ATOM/HETATM, MODEL/ENDMDL, ANISOU, residue identity, atom identity, element, charge, occupancy, B factor, coordinates, SEQRES entities, selected header metadata, AUTHOR, and CRYST1. Connectivity, TER/entity splitting semantics, secondary-structure, and other unported records are not complete PDB support.",
};

pub const BIO_MMCIF_ATOM_SITE_SUBSET_READ_FEATURE: FeatureSpec = FeatureSpec {
    name: "bio.mmcif.atom_site_subset.read",
    category: FeatureCategory::Io,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Experimental Gemmi-aligned mmCIF _atom_site subset reader into BioStructure. This is the structural IO path and the required front end for any future molecule compatibility input. The reader also models _entity, _entity_poly, _entity_poly_seq, and _struct_asym enough to link chains to entities and preserve polymer sequence. RDKit-derived mmCIF parser work remains deferred unless a Molecule compatibility need is approved. It is not complete mmCIF structure semantics; unsupported source branches remain visible in io::bio.",
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
    parity_sensitive: true,
    docs: "Experimental distance-geometry bounds matrix generation. The current Rust DG bounds surface is \
           source-backed across the selected RDKit baseline: raw BoundsMatrix upper/lower triangle storage, \
           triangle smoothing, 1-2/1-3/1-4/1-5 bound setting, VDW lower bounds, collectBondsAndAngles, \
           both setTopolBounds overloads, and GetMoleculeBoundsMatrix-style wrapper defaults are implemented \
           with focused strict tests. This is a port-closure statement for the audited DG bounds scope, not \
           a blanket RDKit parity guarantee for every possible molecule/input outside that baseline.",
};

pub const BIO_SELECTION_FEATURE: FeatureSpec = FeatureSpec {
    name: "bio.selection",
    category: FeatureCategory::BioSelection,
    status: SupportStatus::Experimental,
    parity_sensitive: false,
    docs: "Experimental BioStructure selection and filtering operations (e.g. remove_waters).",
};
