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

pub const MOL2_READ_FEATURE: FeatureSpec = FeatureSpec {
    name: "mol2.read",
    category: FeatureCategory::Io,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "RDKit-compatible Tripos MOL2 reading is source-ported from `Mol2FileParser.cpp` for the exposed `Mol2FileToMol`/`Mol2BlockToMol` profile, including `Mol2ParserParams` controls for sanitize, removeHs, CORINA variant, and cleanupSubstructures. The feature remains experimental while broader fixture parity and marker audit work continues.",
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
    docs: "Experimental RDKit-aligned 2D depiction surface with value semantics. The active Rust path includes parameterized compute2DCoords entrypoints, preferCoordGen/forceRDKit routing, ring-template registry loading, mimic-distance embedding, constrained 2D/3D depiction matching, normalize/straighten helpers, and registered with_2d_coordinates exposure used by batch, MolBlock, and drawing callers. CoordGen-backed runtime branches are not available in this build and fail explicitly instead of silently diverging; final whole-surface audit/validation remains tracked separately.",
};

pub const COORDINATE_EDIT_FEATURE: FeatureSpec = FeatureSpec {
    name: "coordinates.edit",
    category: FeatureCategory::TopologyOperation,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Experimental coordinate block editing via registered molecule operations. This surface covers replacement of the 2D coordinate block, replacement of an existing 3D conformer, and appending an additional 3D conformer. The goal is to keep coordinate mutation inside the operation-contract path so tests can surface state/contract mismatches instead of bypassing them with ad hoc storage edits.",
};

pub const CONFORMER_GENERATION_FEATURE: FeatureSpec = FeatureSpec {
    name: "coordinates.3d.conformer_generation",
    category: FeatureCategory::TopologyOperation,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Experimental RDKit-aligned distance-geometry conformer generation. The exposed surface uses the source-ported EmbedParameters presets, DG/KDG/ETDG/ETKDG entry points, source-backed seeded and unseeded RNG setup, deterministic explicit-seed single-conformer path, deterministic batch seed policy for multi-conformer generation, pruning, terminal-group symmetrization during symmetry-aware pruning, coordMap, CPCI, custom bounds-matrix size validation, stereo/chiral checks, macrocycle and small-ring torsion paths. Final marker audit: no first-axis `RDKit❌❌` block remains in the audited conformer-generation path; residual `RDKit✔️❌`, `RDKit✔️❗`, and `RDKit❗✔️` markers remain in the bounds-builder helper surface and should not be overstated as blanket parity closure.",
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
    docs: "Mixed fingerprint surface. The enumerated Morgan branches are source-backed and tested for exact RDKit equality across sparse-count, sparse-bit, hashed-count, explicit-bit, AdditionalOutput, and the committed branch matrix. MACCS is tested as the exact RDKit raw 167-bit vector and COSMolKit 166-bit projection on targeted fixtures plus the small and strict SMILES profiles. These are validation boundaries, not a claim that every RDKit fingerprint API or every unmodeled input-state preparation branch is ported. Morgan branches that require unported preparation or supplied-pattern behavior return FingerprintError. RDKFingerprintMol/topological and Avalon are not implemented: their former local path-hash approximations were removed and their public calls now return structured unsupported errors. No approximate, correlation-compatible, or 99.9%-matching vector is accepted.",
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
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Experimental SVG/PNG molecule renderer with strict RDKit MolDraw2D final-SVG golden coverage for the shared corpus. \
           The passing SVG parity boundary must not be generalized to marker-open branches. \
           Atom labels, common bond geometry, radical dots, scale calculation, selected annotations, SGroup data, brackets, variable bonds, close-contact markers, highlights, SVG metadata, data-tag attributes, and CSS class output are modeled for the covered surface. \
           Link-node drawing, StereoGroup masking, atomRegions, and other marker-open branches remain unfinished until their RDKit source helpers are ported and covered by exact final-output parity tests. \
           PNG output is a rasterization of the local SVG via usvg+resvg, not a RDKit Cairo/Qt bit-parity target.",
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

pub const INCHI_FEATURE: FeatureSpec = FeatureSpec {
    name: "inchi.scalar",
    category: FeatureCategory::Io,
    status: SupportStatus::SupportedWithRdkitParity {
        rdkit_version: "2026.03.1",
    },
    parity_sensitive: true,
    docs: "The four audited scalar APIs (MolToInchi, MolToInchiKey, InchiToInchiKey, and MolFromInchi) reproduce pinned RDKit 2026.03.1 and official InChI v1.07.5 behavior exactly where the official C source defines behavior. The active NormalizeAndCompare initial-buffer allocation-failure path executes undefined C behavior; COSMolKit intentionally returns a deterministic structured allocation error there. MolBlock, SDF/V3000, IXA, AuxInfo conversion, incremental INCHIGEN, version-query, and extended-polymer InChI APIs are outside this public feature and remain frozen or unsupported.",
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
    docs: "Experimental Gemmi-aligned PDB structural reader into BioStructure. This is the structural IO path and the required front end for future RDKit-compatible molecule input. The public feature name keeps the historical subset label for API stability, but the current reader surface covers ATOM/HETATM, MODEL/ENDMDL, ANISOU, residue and chain identity, TER semantics, SEQRES entities, DBREF, SSBOND/LINK/CISPEP, MODRES, selected header metadata, AUTHOR, CRYST1, SCALE, ORIGX, and MTRIX/NCS records. Remaining unsupported Gemmi branches fail explicitly and stay marked in io::bio.",
};

pub const BIO_MMCIF_ATOM_SITE_SUBSET_READ_FEATURE: FeatureSpec = FeatureSpec {
    name: "bio.mmcif.atom_site_subset.read",
    category: FeatureCategory::Io,
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Experimental Gemmi-aligned mmCIF/mmJSON structural reader into BioStructure. This is the structural IO path and the required front end for any future molecule compatibility input. The public feature name keeps the historical atom-site subset label for API stability, but the current reader surface also covers mmJSON dispatch, _entity, _entity_poly, _entity_poly_seq, _struct_ref/_struct_ref_seq, _struct_asym, _struct_conn, _struct_mon_prot_cis, _pdbx_struct_mod_residue, _pdbx_struct_assembly*, _pdbx_sifts_xref_db, _struct_ncs_oper, crystallographic transforms, and chem-comp CIF handoff through the same dispatch path. RDKit-derived macromolecular parser work remains deferred unless a Molecule compatibility need is approved. Remaining unsupported Gemmi branches fail explicitly and stay marked in io::bio.",
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
    status: SupportStatus::Experimental,
    parity_sensitive: true,
    docs: "Experimental distance-geometry bounds matrix generation. The current Rust DG bounds surface is \
           source-backed across the selected RDKit baseline: raw BoundsMatrix upper/lower triangle storage, \
           triangle smoothing, 1-2/1-3/1-4/1-5 bound setting, VDW lower bounds, collectBondsAndAngles, \
           both setTopolBounds overloads, and GetMoleculeBoundsMatrix-style wrapper defaults are implemented \
           with focused strict tests. The final DG bounds audit found no remaining first-axis `RDKit❌*` \
           gap in the audited call chain, but deliberate `RDKit✔️❌`, `RDKit✔️❗`, and `RDKit❗✔️` markers \
           remain visible for performance and helper-abstraction caveats. This is a port-closure statement \
           for the audited DG bounds scope, not a blanket RDKit parity guarantee for every possible \
           molecule/input outside that baseline.",
};

pub const BIO_SELECTION_FEATURE: FeatureSpec = FeatureSpec {
    name: "bio.selection",
    category: FeatureCategory::BioSelection,
    status: SupportStatus::Experimental,
    parity_sensitive: false,
    docs: "Experimental BioStructure selection and filtering operations (e.g. remove_waters).",
};
