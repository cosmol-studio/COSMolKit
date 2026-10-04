//! Private selection-copying child of `selection.rs` (BIO-COPY packet).
//!
//! Source pins: third_party/gemmi/include/gemmi/model.hpp
//! (`Atom::empty_copy`, `Residue::empty_copy`, `Chain::empty_copy`,
//! `Model::empty_copy`, `Structure::empty_copy`) and
//! third_party/gemmi/include/gemmi/select.hpp
//! (`add_matching_children`, `Selection::copy_selection`). This module is
//! a private semantic child of the canonical CID selection module: it
//! reuses the parent's canonical `Selection` and borrowed row predicates
//! without widening their visibility, and it introduces no second row
//! model, parser or matcher. The narrow cross-crate bridge returns only
//! the copier's replacement blocks and structured errors, not a live value.

use crate::hierarchy::{
    BioAtomRow, BioChainId, BioChainRow, BioCoordinateFormat, BioModelId, BioModelRow,
    BioResidueId, BioResidueRow, BioRowSpan, BioStructureError,
};
use crate::structure_metadata::BioStructureSourceState;
use std::collections::BTreeMap;
use std::error::Error;

/// The single temporary copy output buffer (BIO-COPY packet): exactly the
/// five copied tables — models, chains, residues, atoms and coordinates —
/// and nothing else. It is not a second retained BIO model: it carries no
/// metadata, entities, assemblies, connections or validation state, and it
/// exists only to be assembled into the final detached
/// `BioStructureParts` by the A08 composition.
#[derive(Default)]
pub(crate) struct CopyBuffer {
    pub(crate) models: Vec<BioModelRow>,
    pub(crate) chains: Vec<BioChainRow>,
    pub(crate) residues: Vec<BioResidueRow>,
    pub(crate) atoms: Vec<BioAtomRow>,
    pub(crate) coordinates: Vec<[f64; 3]>,
}

/// Selection-copy cause wrapping the existing typed structure and
/// traversal errors (BIO-COPY A04/A06). No strings: the source chain keeps
/// the underlying typed cause reachable for callers.
#[derive(Debug)]
pub enum SelectionCopyCause {
    Structure(BioStructureError),
    Traverse(crate::selection::BioRowTraverseError),
}

#[derive(Debug)]
pub struct SelectionCopyError {
    pub(crate) cause: SelectionCopyCause,
}

impl SelectionCopyError {
    /// Typed underlying copier failure, without string conversion.
    pub fn cause(&self) -> &SelectionCopyCause {
        &self.cause
    }

    pub(crate) fn from_structure(cause: BioStructureError) -> Self {
        Self {
            cause: SelectionCopyCause::Structure(cause),
        }
    }

    pub(crate) fn from_traverse_chain(cause: crate::selection::BioRowChainError) -> Self {
        Self {
            cause: SelectionCopyCause::Traverse(crate::selection::BioRowTraverseError::Chain(
                cause,
            )),
        }
    }

    pub(crate) fn from_traverse_model(cause: crate::selection::BioRowModelError) -> Self {
        Self {
            cause: SelectionCopyCause::Traverse(crate::selection::BioRowTraverseError::Model(
                cause,
            )),
        }
    }
}

impl std::fmt::Display for SelectionCopyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "selection copy failed: {}", self.cause)
    }
}

impl std::fmt::Display for SelectionCopyCause {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Structure(e) => write!(f, "{e}"),
            Self::Traverse(e) => write!(f, "{e}"),
        }
    }
}

impl Error for SelectionCopyError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        Some(match &self.cause {
            SelectionCopyCause::Structure(e) => e as &(dyn Error + 'static),
            SelectionCopyCause::Traverse(e) => e as &(dyn Error + 'static),
        })
    }
}

impl Error for SelectionCopyCause {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        Some(match self {
            Self::Structure(cause) => cause,
            Self::Traverse(cause) => cause,
        })
    }
}

/// A03 output: the reconstructed row plus the atom's exact source
/// position for the copy buffer's coordinate transport. Gemmi's `Atom`
/// carries `Position pos` inline; CK stores positions in the shared
/// `BioCoordinateBlock`, so the detached copy pairs them explicitly.
#[must_use]
#[derive(Clone, Debug)]
pub(crate) struct CopiedAtom {
    pub(crate) row: BioAtomRow,
    pub(crate) position: [f64; 3],
}

/// Port of Gemmi `Atom::empty_copy` plus coordinate transport (BIO-COPY
/// A03).
///
/// Gemmi's `Atom` stores `Position pos` inline, so the implicit
/// memberwise copy carries the position with the row. CK keeps positions
/// in the shared `BioCoordinateBlock`, so the detached copy returns the
/// reconstructed row together with the atom's exact source position for
/// the copy buffer to append.
#[must_use]
pub(crate) fn copy_selected_atom(
    atom: &BioAtomRow,
    new_residue_id: BioResidueId,
    position: [f64; 3],
) -> CopiedAtom {
    // Gemmi✔️✔️: Atom empty_copy() const { return Atom(*this); }
    // Behavior review: the implicit C++ copy constructor copies every
    // member (`name`, `altloc`, `charge`, `element`, `calc_flag`, `flag`,
    // `tls_group_id`, `serial`, `fraction`, `pos`, `occ`, `b_iso`,
    // `aniso`) memberwise. CK reconstructs the canonical `BioAtomRow`
    // through its value constructor with the accepted source row's exact
    // represented fields — name, element, isotope mass number, altloc,
    // formal charge, calc flag, occupancy, b_iso, anisou, tls group id,
    // fraction and the `AtomSourceIds` source identifiers — so every
    // stored f64 bit (including NaN/negative-zero occupancy/anisou/fraction)
    // and the raw PDB serial survive unchanged. The only difference from
    // the source copy is the reassigned local parent `residue_id`
    // argument, which Gemmi obtains structurally by placing the copy
    // inside the new residue; the source serial is NOT reinterpreted as
    // the new local id. The `flag` member has no CK representation and
    // has no stored counterpart to preserve.
    // Complexity review: O(1) field moves plus one copied `AtomName`,
    // `Element` and `AtomSourceIds` (a `Copy`-sized `Option<PdbAtomSerial>`)
    // — exactly the implicit memberwise copy's cost; no allocation.
    CopiedAtom {
        row: BioAtomRow::new(
            new_residue_id,
            atom.name(),
            atom.element(),
            atom.isotope_mass_number(),
            atom.altloc(),
            atom.formal_charge(),
            atom.calc_flag(),
            atom.occupancy(),
            atom.b_iso(),
            *atom.anisou(),
            atom.tls_group_id(),
            atom.fraction(),
            atom.source().clone(),
        ),
        position,
    }
}

/// Port of Gemmi `Residue::empty_copy` plus the residue-to-atom
/// `add_matching_children` recursion (BIO-COPY A04).
///
/// The residue has already been accepted by the parent chain-level gate;
/// this helper evaluates the canonical borrowed atom predicate exactly
/// once per visited source atom, appends the A03 copies in source order,
/// and constructs the copied residue with remapped parent/span.
#[must_use]
pub(crate) fn append_selected_residue(
    buffer: &mut CopyBuffer,
    selection: &super::Selection,
    residue: &BioResidueRow,
    new_chain_id: BioChainId,
    atoms: &[BioAtomRow],
    positions: &[[f64; 3]],
    input_format: BioCoordinateFormat,
) -> Result<BioResidueId, SelectionCopyError> {
    // Gemmi✔️✔️: Residue empty_copy() const {
    // Gemmi✔️✔️:   Residue res((ResidueId&)*this);
    // Gemmi✔️✔️:   res.subchain = subchain;
    // Gemmi✔️✔️:   res.entity_id = entity_id;
    // Gemmi✔️✔️:   res.label_seq = label_seq;
    // Gemmi✔️✔️:   res.entity_type = entity_type;
    // Gemmi✔️✔️:   res.het_flag = het_flag;
    // Gemmi✔️✔️:   res.flag = flag;
    // Gemmi✔️✔️:   res.sifts_unp = sifts_unp;
    // Gemmi✔️✔️:   return res;
    // Gemmi✔️✔️: }
    // Behavior review: the source empty copy transfers the residue
    // identity (ResidueId: auth seq + icode + segment, carried by CK's
    // `ResidueSourceIds::seq_id`/`segment_id`), `subchain` and
    // `label_seq` (`label_seq_id`/`subchain_id`), `entity_id`/`entity_type`
    // (`entity_id`/`entity_kind`, with the source label entity string in
    // `label_entity_id`), `het_flag`, the custom `flag` (no CK
    // representation, nothing stored to preserve) and `sifts_unp`
    // verbatim. The parent `chain_id` and the local atom span are CK's
    // structural reassignments, which Gemmi obtains by placing the copy
    // inside the target chain; the source seqid/segment names are NOT
    // reinterpreted as the new local ids. The copied row stores the
    // caller-provided NEW chain id. Accepted residues with zero
    // matching atoms are kept, exactly like the source `push_back` before
    // the recursion.
    // Complexity review: O(residue atoms) predicate evaluations plus one
    // amortized push per copied atom and one residue construction. The
    // borrowed per-atom predicate work is allocation-free; the residue
    // construction necessarily clones `ResidueSourceIds`, whose two owned
    // strings (`subchain_id`, `label_entity_id`) allocate — the same
    // per-copied-residue string copies the source `empty_copy` performs —
    // and vector growth is amortized. No per-atom allocation and no
    // repeated whole-table scans: the atoms are the caller's borrowed
    // contiguous span.
    //
    // Gemmi✔️✔️:   for (const auto& orig_child : orig.children())
    // Gemmi✔️✔️:     if (matches(orig_child)) {
    // Gemmi✔️✔️:       target.children().push_back(orig_child.empty_copy());
    // Gemmi✔️✔️:       add_matching_children(orig_child, target.children().back());
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:   void add_matching_children(const Atom&, Atom&) const {}
    // Behavior review: `matches(atom)` is the canonical borrowed R08 atom
    // adapter evaluated exactly once per visited atom. The atom-level
    // recursion is the empty overload — `copy_selected_atom` (A03) is the
    // `empty_copy()` call in the same loop. The flag byte handed to the
    // adapter is the already-proven inert zero input (empty CID flag
    // patterns accept every value; C01/C24 contract) — no invented flags.
    // Complexity review: one pass in source order; each accepted atom
    // appends one row and one position.
    let new_residue_id = BioResidueId::new(buffer.residues.len() as u32);
    let atom_start = buffer.atoms.len() as u32;
    for (atom, position) in atoms.iter().zip(positions) {
        if super::bio_atom_row_matches(selection, atom, input_format, 0) {
            let copied = copy_selected_atom(atom, new_residue_id, *position);
            buffer.atoms.push(copied.row);
            buffer.coordinates.push(copied.position);
        }
    }
    let atom_span = BioRowSpan::new(atom_start, buffer.atoms.len() as u32 - atom_start)
        .map_err(SelectionCopyError::from_structure)?;
    buffer.residues.push(BioResidueRow::new(
        new_chain_id,
        atom_span,
        residue.name(),
        residue.residue_info_kind(),
        residue.entity_kind(),
        residue.entity_id(),
        residue.het_flag(),
        residue.source().clone(),
        residue.sifts_unp(),
    ));
    Ok(new_residue_id)
}

/// Port of Gemmi `Chain::empty_copy` plus the chain-to-residue
/// `add_matching_children` recursion (BIO-COPY A05).
///
/// The chain has already been accepted by the parent model-level gate;
/// this helper scans the chain's original residue span once, reuses the
/// canonical borrowed residue predicate once per visited residue, appends
/// only accepted residues through A04, and keeps an accepted chain even
/// when it ends up with zero accepted residues.
pub(crate) fn append_selected_chain(
    buffer: &mut CopyBuffer,
    selection: &super::Selection,
    chain: &BioChainRow,
    new_model_id: BioModelId,
    residues: &[BioResidueRow],
    atoms: &[BioAtomRow],
    positions: &[[f64; 3]],
    input_format: BioCoordinateFormat,
) -> Result<BioChainId, SelectionCopyError> {
    // Gemmi✔️✔️: Chain empty_copy() const { return Chain(name); }
    // Behavior review: the source empty copy keeps ONLY the chain name —
    // the copied chain starts with no residues. CK's `BioChainRow` carries
    // the name as the `ChainSourceIds` auth/label pair, which is cloned
    // verbatim with no normalization, plus the entity linkage and kind;
    // the local `model_id` parent and residue span are the structural
    // reassignments Gemmi obtains by placing the copy inside the target
    // model. Zero accepted residues still produce the chain row, exactly
    // like the source `push_back` before its recursion.
    // Complexity review: one chain-row construction (cloning the source id
    // pair, the only allocation the source string copy also performs) plus
    // one pass over the chain's residue span.
    //
    // Gemmi✔️✔️:   for (const auto& orig_child : orig.children())
    // Gemmi✔️✔️:     if (matches(orig_child)) {
    // Gemmi✔️✔️:       target.children().push_back(orig_child.empty_copy());
    // Gemmi✔️✔️:       add_matching_children(orig_child, target.children().back());
    // Gemmi✔️✔️:     }
    // Behavior review: `matches(residue)` is the canonical borrowed R06
    // residue adapter evaluated exactly once per visited residue; the
    // residue flag byte is the proven inert zero input (empty CID flag
    // patterns accept every value) — no invented flags, no name
    // normalization. Accepted residues recurse through A04
    // (`append_selected_residue`) with the NEW local chain id and each
    // residue's own contiguous atom span — no selected_atom_ids
    // flatten-and-rebuild and no repeated whole-table scans.
    // Complexity review: single pass over the chain's residue span in
    // source order; each accepted residue costs its own atom-span walk.
    let new_chain_id = BioChainId::new(buffer.chains.len() as u32);
    let residue_start = buffer.residues.len() as u32;
    let chain_span = chain.residue_span();
    for residue in
        &residues[chain_span.start() as usize..(chain_span.start() + chain_span.len()) as usize]
    {
        if super::bio_residue_row_matches(selection, residue, input_format, 0) {
            let atom_span = residue.atom_span();
            let atom_start = atom_span.start() as usize;
            let atom_end = atom_start + atom_span.len() as usize;
            append_selected_residue(
                buffer,
                selection,
                residue,
                new_chain_id,
                &atoms[atom_start..atom_end],
                &positions[atom_start..atom_end],
                input_format,
            )?;
        }
    }
    let residue_span = BioRowSpan::new(residue_start, buffer.residues.len() as u32 - residue_start)
        .map_err(SelectionCopyError::from_structure)?;
    buffer.chains.push(BioChainRow::new(
        new_model_id,
        chain.entity_id(),
        residue_span,
        chain.kind(),
        chain.source().clone(),
    ));
    Ok(new_chain_id)
}

/// Port of Gemmi `Model::empty_copy` plus the model-to-chain
/// `add_matching_children` recursion (BIO-COPY A06).
///
/// The model has already been accepted by the structure-level gate; this
/// helper scans the model's chain span once, reuses the canonical borrowed
/// chain predicate once per visited chain, appends only accepted chains
/// through A05, preserves the actual stored model number, and emits a model
/// row with its remapped chain span even when no chain matches.
///
/// # Errors
///
/// Returns the existing typed structure error (wrapped in
/// [`SelectionCopyError`]) from A04/A05 span construction, and propagates
/// the canonical chain gate's typed missing canonical chain-name error
/// unchanged — no guessed replacement name, no silent skip.
pub(crate) fn append_selected_model(
    buffer: &mut CopyBuffer,
    selection: &super::Selection,
    model: &BioModelRow,
    chains: &[BioChainRow],
    residues: &[BioResidueRow],
    atoms: &[BioAtomRow],
    positions: &[[f64; 3]],
    input_format: BioCoordinateFormat,
) -> Result<(), SelectionCopyError> {
    // Gemmi✔️✔️: Model empty_copy() const { return Model(num); }
    // Behavior review: the source empty copy keeps ONLY the model number;
    // the copied model starts with no chains. CK's `BioModelRow` stores
    // the reader-assigned `source_model_number` verbatim (including
    // `None`, which the A08 model gate rejects with the typed
    // missing-representation error BEFORE this helper runs, exactly as
    // the source predicate always sees a reader-assigned number); the
    // chain span is the structural reassignment Gemmi obtains by placing
    // the copy inside the target structure. Zero accepted chains still
    // emit the model row, exactly like the source `push_back` before its
    // recursion. No guessed replacement number is ever substituted.
    // Complexity review: one model-row construction plus one pass over
    // the model's chain span; no allocation of its own beyond the
    // amortized vector growth in the appended child rows.
    //
    // Gemmi✔️✔️:   for (const auto& orig_child : orig.children())
    // Gemmi✔️✔️:     if (matches(orig_child)) {
    // Gemmi✔️✔️:       target.children().push_back(orig_child.empty_copy());
    // Gemmi✔️✔️:       add_matching_children(orig_child, target.children().back());
    // Gemmi✔️✔️:     }
    // Behavior review: `matches(chain)` is the canonical borrowed R04
    // chain adapter evaluated exactly once per visited chain. This is
    // the REAL copier error propagation point: a chain row with neither
    // auth nor label identity returns the typed
    // `BioRowChainError::MissingCanonicalChainName` from inside the copy
    // recursion and is propagated unchanged (source order preserved —
    // the error surfaces at the first offending chain, no later chain
    // is visited, matching the source's exception out of the loop).
    // Accepted chains recurse through A05 (`append_selected_chain`).
    // Complexity review: single pass over the model's chain span in
    // source order; each accepted chain costs its own residue-span walk.
    let model_id = BioModelId::new(buffer.models.len() as u32);
    let chain_start = buffer.chains.len() as u32;
    let model_span = model.chain_span();
    for chain in
        &chains[model_span.start() as usize..(model_span.start() + model_span.len()) as usize]
    {
        if super::bio_chain_row_matches(selection, chain, input_format)
            .map_err(SelectionCopyError::from_traverse_chain)?
        {
            append_selected_chain(
                buffer,
                selection,
                chain,
                model_id,
                residues,
                atoms,
                positions,
                input_format,
            )?;
        }
    }
    let chain_span = BioRowSpan::new(chain_start, buffer.chains.len() as u32 - chain_start)
        .map_err(SelectionCopyError::from_structure)?;
    buffer
        .models
        .push(BioModelRow::new(chain_span, model.source_model_number()));
    Ok(())
}

/// Port of Gemmi `Structure::empty_copy`'s source-state assignments
/// (BIO-COPY A07).
///
/// The pinned `empty_copy` assigns `name`, `has_origx`, `origx`, `info`,
/// `raw_remarks` and `resolution` (among its other structure members) and
/// does NOT assign `conect_map`, `has_d_fraction`, `non_ascii_line` or
/// `ter_status`, which therefore keep the fresh `Structure st;` member
/// state.
#[must_use]
pub(crate) fn empty_copy_source_state(source: &BioStructureSourceState) -> BioStructureSourceState {
    // Gemmi✔️✔️: Structure st;
    // Gemmi✔️✔️: st.name = name;
    // Gemmi✔️✔️: st.has_origx = has_origx;
    // Gemmi✔️✔️: st.origx = origx;
    // Gemmi✔️✔️: st.info = info;
    // Gemmi✔️✔️: st.raw_remarks = raw_remarks;
    // Gemmi✔️✔️: st.resolution = resolution;
    // Gemmi✔️✔️: return st;
    // Behavior review: the assigned fields copy memberwise — name,
    // has_origx, origx (a Copy transform), info, raw_remarks (source order
    // and duplicates preserved) and resolution (bit-exact f64). The four
    // unassigned members (conect_map, has_d_fraction, non_ascii_line,
    // ter_status) keep the fresh-structure member state — empty map, false,
    // 0, NUL — exactly because the pinned empty_copy does not assign them.
    // This is Gemmi empty_copy behavior, not data-loss cleanup, and it is
    // NOT the Protein projection's independently approved preservation
    // contract: no unrelated metadata is reset.
    // Complexity review: O(len) copies of the assigned name/info/remarks
    // storage (the same std::string/map/vector copies the source performs);
    // the defaulted members are constant-time.
    BioStructureSourceState {
        name: source.name.clone(),
        resolution: source.resolution,
        conect_map: BTreeMap::new(),
        has_d_fraction: false,
        non_ascii_line: 0,
        ter_status: 0,
        has_origx: source.has_origx,
        origx: source.origx,
        info: source.info.clone(),
        raw_remarks: source.raw_remarks.clone(),
    }
}

/// Port of Gemmi `Selection::copy_selection` over a whole structure
/// (BIO-COPY A08): the one canonical detached copier.
///
/// Gates each source model once with the canonical borrowed model adapter,
/// recurses through A06/A05/A04/A03 into the five-vector construction
/// buffer, replaces the source state through A07, shares every unchanged
/// metadata block immutably, and validates the resulting CK row hierarchy
/// before returning the detached `BioStructureData`.
///
/// # Errors
///
/// Propagates the typed model/chain gate errors and span-structure errors
/// unchanged — the same REAL returned error chain proven at A06.
pub(crate) fn copy_bio_selection<'a>(
    source: impl Into<crate::hierarchy::BioStructureCopySource<'a>>,
    selection: &crate::selection::BioSelectionData,
) -> Result<crate::hierarchy::BioStructureData, SelectionCopyError> {
    let source = source.into();
    // Gemmi✔️✔️:   // copy all but models (in general, empty_copy copies all but children)
    // Gemmi✔️✔️:   Structure empty_copy() const {
    // Gemmi✔️✔️:     Structure st;
    // Gemmi✔️✔️:     st.name = name;
    // Gemmi✔️✔️:     st.cell = cell;
    // Gemmi✔️✔️:     st.spacegroup_hm = spacegroup_hm;
    // Gemmi✔️✔️:     st.ncs = ncs;
    // Gemmi✔️✔️:     st.entities = entities;
    // Gemmi✔️✔️:     st.connections = connections;
    // Gemmi✔️✔️:     st.cispeps = cispeps;
    // Gemmi✔️✔️:     st.mod_residues = mod_residues;
    // Gemmi✔️✔️:     st.helices = helices;
    // Gemmi✔️✔️:     st.sheets = sheets;
    // Gemmi✔️✔️:     st.assemblies = assemblies;
    // Gemmi✔️✔️:     st.meta = meta;
    // Gemmi✔️✔️:     st.input_format = input_format;
    // Gemmi✔️✔️:     st.has_origx = has_origx;
    // Gemmi✔️✔️:     st.origx = origx;
    // Gemmi✔️✔️:     st.info = info;
    // Gemmi✔️✔️:     st.raw_remarks = raw_remarks;
    // Gemmi✔️✔️:     st.resolution = resolution;
    // Gemmi✔️✔️:     return st;
    // Gemmi✔️✔️:   }
    // Behavior review (Structure::empty_copy assignments): the source
    // memberwise-copies every non-child member. CK materializes the
    // structure-level assignments as the shared immutable blocks below —
    // entities/connections/cispeps/mod_residues/helices/sheets/assemblies
    // (st.entities..st.assemblies), crystal+ncs_operators (st.cell /
    // st.spacegroup_hm / st.ncs), metadata (st.meta) and input_format —
    // while the `Structure st;` member-default members that the source
    // does NOT re-assign in the source_state subset are handled by A07
    // (`st.has_origx/origx/info/raw_remarks/resolution` copies plus the
    // defaulted conect_map/has_d_fraction/non_ascii_line/ter_status). The
    // `st.name` copy lands in A07's name. Every source metadata name is
    // retained even for an empty output; no assembly expansion/resolution.
    // Complexity review (Structure::empty_copy assignments): the source
    // performs one memberwise copy per container; CK Arc-shares each
    // immutable block (refcount bump) instead — an intentional allocation
    // improvement that preserves semantics because detached blocks are
    // immutable, so no copied bytes can diverge from the source's.
    //
    // Gemmi✔️✔️:   void add_matching_children(const T& orig, T& target) const {
    // Gemmi✔️✔️:     for (const auto& orig_child : orig.children())
    // Gemmi✔️✔️:       if (matches(orig_child)) {
    // Gemmi✔️✔️:         target.children().push_back(orig_child.empty_copy());
    // Gemmi✔️✔️:         add_matching_children(orig_child, target.children().back());
    // Gemmi✔️✔️:       }
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   void add_matching_children(const Atom&, Atom&) const {}
    // Behavior review (structure-to-model recursion): the generic template
    // instantiated at T=Structure recurses Structure->Model->Chain->Residue
    // and terminates at the Atom no-op. CK runs the canonical borrowed model
    // gate once per visited source model (typed MissingModelNumber cause
    // propagates — no guessed number) and recurses accepted models through
    // A06, which repeats the same loop at the next level down to A03 (the
    // Atom empty_copy inside the residue loop). Source order and
    // selected-but-empty parents are preserved by construction.
    // Complexity review (recursion): one pass per visited row, matching the
    // source recursion; no repeated whole-table scans or flatten/rebuild.
    //
    // Typed CK input validation (BIO-COPY-A08-INPUT): the pinned Gemmi
    // recursion indexes live in-memory children and cannot observe a
    // malformed hierarchy, but CK rows are sliced by stored spans BEFORE
    // any traversal — so the source hierarchy is validated up front with
    // the existing typed `validate` (BioStructureView), failing closed with
    // `BioStructureError` (e.g. RowSpanOutOfBounds) instead of panicking in
    // a slice or silently disappearing behind a non-matching selection.
    // This is a CK input-safety boundary on top of the source behavior,
    // not a Gemmi behavior claim; the final `copied.validate()` below
    // independently re-validates the produced hierarchy.
    // Complexity review (input validation): the up-front source.validate()
    // is a FULL O(total rows) pass — table-length checks, per-level
    // span/parent/coverage walks over every model/chain/residue row and
    // per-reference entity lookups whose subchain membership scans the
    // referenced entity's subchain list. It runs even when the model
    // predicate accepts nothing, so whole-operation cost is NOT universally
    // dominated by the selection recursion; the added input pass is the CK
    // input-safety price over the source behavior.
    // Complexity review (recursion): the gated traversal itself visits only
    // VISITED rows — one pass per visited row matching the source recursion,
    // with per-visited-row predicate list scans — so it scales with the
    // selection, not the whole table.
    source
        .validate()
        .map_err(SelectionCopyError::from_structure)?;
    // Gemmi✔️✔️:   template<typename T>
    // Gemmi✔️✔️:   T copy_selection(const T& orig) const {
    // Gemmi✔️✔️:     T copied = orig.empty_copy();
    // Gemmi✔️✔️:     add_matching_children(orig, copied);
    // Gemmi✔️✔️:     return copied;
    // Gemmi✔️✔️:   }
    // Behavior review: composed entry — the validated source's empty_copy
    // state plus the gated recursion below produce the one canonical
    // detached output; the final validate below rejects any structural
    // defect in the produced rows before the value is returned.
    // Complexity review (whole operation, separate parts): (1) full input
    // validation as accounted above — a fixed whole-table cost that also
    // runs when nothing matches; (2) visited-row recursion as accounted
    // above; (3) SELECTED-OUTPUT validation — the final copied.validate()
    // is a second full pass, but over the OUTPUT rows only (O(selected
    // rows), with the same entity subchain-list membership scans); (4) the
    // A07 source-state replacement's owned string/container copies (name,
    // info map, raw remarks) allocate per copied container, exactly the
    // source empty_copy's memberwise copies; (5) the five copied tables'
    // amortized vector growth. The Arc-sharing note above covers ONLY the
    // unchanged metadata blocks; it does not by itself establish
    // whole-operation performance equivalence with the source — the two
    // validation passes and the string copies are additive CK/source costs.
    let inner = &selection.selection;
    let mut buffer = CopyBuffer::default();
    for model in source.models.iter() {
        if super::bio_model_row_matches(inner, model)
            .map_err(SelectionCopyError::from_traverse_model)?
        {
            append_selected_model(
                &mut buffer,
                inner,
                model,
                &source.chains,
                &source.residues,
                &source.atoms,
                source.coordinates.positions(),
                source.input_format,
            )?;
        }
    }
    let copied = crate::hierarchy::BioStructureData {
        input_format: source.input_format,
        models: std::sync::Arc::new(std::mem::take(&mut buffer.models)),
        chains: std::sync::Arc::new(std::mem::take(&mut buffer.chains)),
        residues: std::sync::Arc::new(std::mem::take(&mut buffer.residues)),
        atoms: std::sync::Arc::new(std::mem::take(&mut buffer.atoms)),
        entities: std::sync::Arc::clone(&source.entities),
        connections: std::sync::Arc::clone(&source.connections),
        cispeps: std::sync::Arc::clone(&source.cispeps),
        mod_residues: std::sync::Arc::clone(&source.mod_residues),
        helices: std::sync::Arc::clone(&source.helices),
        sheets: std::sync::Arc::clone(&source.sheets),
        metadata: std::sync::Arc::clone(&source.metadata),
        source_state: std::sync::Arc::new(empty_copy_source_state(&source.source_state)),
        coordinates: std::sync::Arc::new(crate::hierarchy::BioCoordinateBlock::new(
            std::mem::take(&mut buffer.coordinates),
        )),
        crystal: std::sync::Arc::clone(&source.crystal),
        ncs_operators: std::sync::Arc::clone(&source.ncs_operators),
        assemblies: std::sync::Arc::clone(&source.assemblies),
    };
    copied
        .validate()
        .map_err(SelectionCopyError::from_structure)?;
    Ok(copied)
}

/// Complete replacement output of the detached selection copier.
/// Unchanged metadata stays with the caller; only the declared six blocks are returned.
#[derive(Debug)]
pub struct BioSelectionBlocks {
    pub models: std::sync::Arc<Vec<crate::BioModelRow>>,
    pub chains: std::sync::Arc<Vec<crate::BioChainRow>>,
    pub residues: std::sync::Arc<Vec<crate::BioResidueRow>>,
    pub atoms: std::sync::Arc<Vec<crate::BioAtomRow>>,
    pub source_state: std::sync::Arc<crate::BioStructureSourceState>,
    pub coordinates: std::sync::Arc<crate::BioCoordinateBlock>,
}

/// Delegate to the sole validated Gemmi copier and move its replacement blocks.
/// This is a domain boundary, not a second implementation or a live-value constructor.
pub fn copy_selection_blocks(
    source: crate::hierarchy::BioStructureCopySource<'_>,
    selection: &crate::selection::BioSelectionData,
) -> Result<BioSelectionBlocks, SelectionCopyError> {
    let copied = copy_bio_selection(source, selection)?;
    Ok(BioSelectionBlocks {
        models: copied.models,
        chains: copied.chains,
        residues: copied.residues,
        atoms: copied.atoms,
        source_state: copied.source_state,
        coordinates: copied.coordinates,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn bio_copy_public_block_boundary_moves_complete_output_and_keeps_source() {
        use crate::hierarchy::{BioStructureCopySource, BioStructureData};
        use std::sync::Arc;
        let source =
            BioStructureData::from_parts(a09_parts(Some(b"A"), Some("S1"), Some(b"B"), Some("S2")))
                .unwrap();
        let selection = crate::selection::BioSelectionData {
            selection: super::super::Selection::default(),
        };
        let peer = source.clone();
        let positions = source
            .coordinates
            .positions()
            .iter()
            .map(|p| p.map(f64::to_bits))
            .collect::<Vec<_>>();
        let blocks =
            copy_selection_blocks(BioStructureCopySource::from(&source), &selection).unwrap();
        assert_eq!(blocks.models.len(), 2);
        assert_eq!(blocks.chains.len(), 2);
        assert_eq!(blocks.residues.len(), 2);
        assert_eq!(blocks.atoms.len(), 2);
        assert_eq!(
            blocks
                .coordinates
                .positions()
                .iter()
                .map(|p| p.map(f64::to_bits))
                .collect::<Vec<_>>(),
            positions
        );
        assert!(!Arc::ptr_eq(&blocks.atoms, &source.atoms));
        assert!(!Arc::ptr_eq(&blocks.coordinates, &source.coordinates));
        assert_eq!(blocks.source_state.name, source.source_state.name);
        assert!(blocks.source_state.conect_map.is_empty());
        assert!(Arc::ptr_eq(&source.atoms, &peer.atoms));
        assert!(Arc::ptr_eq(&source.coordinates, &peer.coordinates));
        assert!(Arc::ptr_eq(&source.assemblies, &peer.assemblies));
        assert_eq!(
            source
                .coordinates
                .positions()
                .iter()
                .map(|p| p.map(f64::to_bits))
                .collect::<Vec<_>>(),
            positions
        );
        let mut invalid = source.clone();
        Arc::make_mut(&mut invalid.models)[0] =
            crate::BioModelRow::new(crate::BioRowSpan::new(0, 5).unwrap(), Some(1));
        let error =
            copy_selection_blocks(BioStructureCopySource::from(&invalid), &selection).unwrap_err();
        assert!(matches!(
            error.cause(),
            SelectionCopyCause::Structure(BioStructureError::RowSpanOutOfBounds {
                start: 0,
                len: 5,
                table_len: 2
            })
        ));
        assert!(
            error
                .source()
                .unwrap()
                .downcast_ref::<BioStructureError>()
                .is_some()
        );
    }

    use crate::hierarchy::BioCalcFlag;
    use crate::source_ids::{AltLocLabel, AtomName, AtomSourceIds, PdbAtomSerial};
    use cosmolkit_types::Element;

    fn row(
        residue_id: BioResidueId,
        name: &[u8],
        element: Element,
        isotope: Option<u16>,
        altloc: Option<AltLocLabel>,
        charge: i8,
        calc_flag: BioCalcFlag,
        occupancy: f64,
        b_iso: f64,
        anisou: [f64; 6],
        tls: i16,
        fraction: f64,
        serial: Option<i32>,
    ) -> BioAtomRow {
        BioAtomRow::new(
            residue_id,
            AtomName::from_ascii(name).unwrap(),
            element,
            isotope,
            altloc,
            charge,
            calc_flag,
            occupancy,
            b_iso,
            anisou,
            tls,
            fraction,
            AtomSourceIds::new(serial.map(PdbAtomSerial::new)),
        )
    }

    #[test]
    fn bio_copy_a03_atom_fields_source_ids_and_parent_reassignment() {
        // H with explicit isotope mass number Some(2) (deuterium), present
        // altloc, NaN occupancy, -0.0 b_iso, NaN/signed-zero anisou bits,
        // a raw PDB serial, and a signed-zero/NaN position — every
        // represented field compared getter-by-getter with bit equality,
        // plus the new local parent id.
        let source = row(
            BioResidueId::new(3),
            b" H  ",
            Element::from_atomic_number(1).unwrap(),
            Some(2),
            Some(AltLocLabel::new(b'A')),
            -2,
            BioCalcFlag::Calculated,
            f64::NAN,
            -0.0,
            [f64::NAN, -0.0, 1.25, -1.25, f64::INFINITY, -f64::INFINITY],
            7,
            0.25,
            Some(4242),
        );
        let position = [-0.0, 1.5, f64::NAN];
        let copied = copy_selected_atom(&source, BioResidueId::new(9), position);

        assert_eq!(copied.row.residue_id(), BioResidueId::new(9));
        assert_ne!(copied.row.residue_id(), source.residue_id());
        assert_eq!(copied.row.name().as_bytes(), b" H  ");
        assert_eq!(copied.row.element(), source.element());
        assert_eq!(copied.row.isotope_mass_number(), Some(2));
        assert_eq!(copied.row.altloc(), source.altloc());
        assert_eq!(copied.row.formal_charge(), -2);
        assert_eq!(copied.row.calc_flag(), BioCalcFlag::Calculated);
        assert_eq!(copied.row.occupancy().to_bits(), f64::NAN.to_bits());
        assert_eq!(copied.row.b_iso().to_bits(), (-0.0_f64).to_bits());
        for (copied_component, source_component) in copied.row.anisou().iter().zip(source.anisou())
        {
            assert_eq!(copied_component.to_bits(), source_component.to_bits());
        }
        assert_eq!(copied.row.tls_group_id(), 7);
        assert_eq!(copied.row.fraction().to_bits(), 0.25_f64.to_bits());
        assert_eq!(copied.row.source().serial(), source.source().serial());
        for (copied_component, source_component) in copied.position.iter().zip(position) {
            assert_eq!(copied_component.to_bits(), source_component.to_bits());
        }
    }

    #[test]
    fn bio_copy_a03_missing_altloc_and_isotope_copy_unchanged() {
        // Missing altloc (None) and unspecified isotope (None): the source
        // representation gaps survive verbatim; a finite position copies
        // component-wise with exact bits.
        let source = row(
            BioResidueId::new(0),
            b" CA ",
            Element::from_atomic_number(6).unwrap(),
            None,
            None,
            0,
            BioCalcFlag::default(),
            1.0,
            20.0,
            [0.0; 6],
            -1,
            0.0,
            None,
        );
        let position = [8.5, -9.25, 0.0];
        let copied = copy_selected_atom(&source, BioResidueId::new(1), position);

        assert_eq!(copied.row.residue_id(), BioResidueId::new(1));
        assert_eq!(copied.row.name().as_bytes(), b" CA ");
        assert_eq!(copied.row.altloc(), None);
        assert_eq!(copied.row.isotope_mass_number(), None);
        assert_eq!(
            copied.row.element(),
            Element::from_atomic_number(6).unwrap()
        );
        assert_eq!(copied.row.formal_charge(), 0);
        assert_eq!(copied.row.calc_flag(), BioCalcFlag::NotSet);
        assert_eq!(copied.row.occupancy().to_bits(), 1.0_f64.to_bits());
        assert_eq!(copied.row.b_iso().to_bits(), 20.0_f64.to_bits());
        assert_eq!(copied.row.anisou(), &[0.0; 6]);
        assert_eq!(copied.row.tls_group_id(), -1);
        assert_eq!(copied.row.fraction().to_bits(), 0.0_f64.to_bits());
        assert_eq!(copied.row.source().serial(), None);
        for (copied_component, source_component) in copied.position.iter().zip(position) {
            assert_eq!(copied_component.to_bits(), source_component.to_bits());
        }
    }

    #[test]
    fn bio_copy_a04_atom_masks_order_fields_spans_parents_and_bits() {
        use crate::hierarchy::{
            BioChainId, BioChainRow, BioCoordinateFormat, BioEntityId, BioModelId, BioResidueRow,
            BioRowSpan, BioSiftsUnpResidue, ChainKind, EntityKind,
        };
        use crate::residue::ResidueInfoKind;
        use crate::source_ids::{
            ChainSourceIds, PdbChainId, PdbSeqId, ResidueName, ResidueSourceIds,
        };

        let source_residue = BioResidueRow::new(
            BioChainId::new(0),
            BioRowSpan::new(0, 2).unwrap(),
            ResidueName::from_ascii(b"ALA").unwrap(),
            ResidueInfoKind::Unknown,
            EntityKind::Polymer,
            Some(BioEntityId::new(5)),
            Some(b'H'),
            ResidueSourceIds::new(
                Some(PdbSeqId::new(101, Some(b'A'))),
                Some(7),
                Some(*b"SEG1"),
                Some("SUB1".to_owned()),
                Some("E9".to_owned()),
            )
            .unwrap(),
            BioSiftsUnpResidue::default(),
        );
        let atoms = [
            row(
                BioResidueId::new(0),
                b" CA ",
                Element::from_atomic_number(6).unwrap(),
                None,
                None,
                0,
                BioCalcFlag::default(),
                1.0,
                20.0,
                [0.0; 6],
                -1,
                0.0,
                Some(11),
            ),
            row(
                BioResidueId::new(0),
                b" N  ",
                Element::from_atomic_number(7).unwrap(),
                None,
                None,
                0,
                BioCalcFlag::default(),
                0.5,
                19.5,
                [0.0; 6],
                -1,
                0.0,
                Some(12),
            ),
        ];
        let positions = [[-0.0, 1.0, 2.0], [3.0, -0.0, f64::NAN]];
        let pdb = BioCoordinateFormat::Pdb;

        let mask = |names: &str| super::super::Selection {
            atom_names: super::super::SelectionList {
                all: false,
                inverted: false,
                list: names.to_owned(),
            },
            ..super::super::Selection::default()
        };
        for (names, expected) in [("O", 0usize), ("CA", 1), ("N", 1), ("", 2)] {
            // The empty-string mask is `all=false, list=""` which matches
            // nothing — for the BOTH case use the default all-accepting
            // selection instead.
            let selection = if names.is_empty() {
                super::super::Selection::default()
            } else {
                mask(names)
            };
            let mut buffer = CopyBuffer::default();
            let new_id = append_selected_residue(
                &mut buffer,
                &selection,
                &source_residue,
                BioChainId::new(2),
                &atoms,
                &positions,
                pdb,
            )
            .unwrap();
            assert_eq!(new_id, BioResidueId::new(0));
            assert_eq!(buffer.residues.len(), 1);
            assert_eq!(buffer.atoms.len(), expected);
            assert_eq!(buffer.coordinates.len(), expected);
            let copied = &buffer.residues[0];
            assert_eq!(copied.chain_id(), BioChainId::new(2));
            assert_eq!(copied.atom_span().start(), 0);
            assert_eq!(copied.atom_span().len() as usize, expected);
            assert_eq!(copied.name().as_bytes(), b"ALA");
            assert_eq!(copied.residue_info_kind(), ResidueInfoKind::Unknown);
            assert_eq!(copied.entity_kind(), EntityKind::Polymer);
            assert_eq!(copied.entity_id(), Some(BioEntityId::new(5)));
            assert_eq!(copied.het_flag(), Some(b'H'));
            assert_eq!(copied.source().seq_id(), source_residue.source().seq_id());
            assert_eq!(copied.source().label_seq_id(), Some(7));
            assert_eq!(copied.source().segment_id(), Some(b"SEG1"));
            assert_eq!(copied.source().subchain_id(), Some("SUB1"));
            assert_eq!(copied.source().label_entity_id().as_deref(), Some("E9"));
            assert_eq!(copied.sifts_unp(), source_residue.sifts_unp());
            for (index, atom) in buffer.atoms.iter().enumerate() {
                assert_eq!(atom.residue_id(), new_id);
                let source_atom = &atoms[if expected == 1 {
                    if names == "N" { 1 } else { 0 }
                } else {
                    index
                }];
                assert_eq!(atom.name().as_bytes(), source_atom.name().as_bytes());
                assert_eq!(atom.source().serial(), source_atom.source().serial());
                for (copied_component, source_component) in buffer.coordinates[index].iter().zip(
                    positions[if expected == 1 {
                        if names == "N" { 1 } else { 0 }
                    } else {
                        index
                    }],
                ) {
                    assert_eq!(copied_component.to_bits(), source_component.to_bits());
                }
            }
        }

        // Existing-nonempty buffer: spans and residue ids continue past the
        // already-copied rows, source order preserved across appends.
        let mut buffer = CopyBuffer::default();
        append_selected_residue(
            &mut buffer,
            &mask("CA"),
            &source_residue,
            BioChainId::new(2),
            &atoms,
            &positions,
            pdb,
        )
        .unwrap();
        let second = append_selected_residue(
            &mut buffer,
            &mask("N"),
            &source_residue,
            BioChainId::new(2),
            &atoms,
            &positions,
            pdb,
        )
        .unwrap();
        assert_eq!(second, BioResidueId::new(1));
        assert_eq!(buffer.residues[1].atom_span().start(), 1);
        assert_eq!(buffer.residues[1].atom_span().len(), 1);
        assert_eq!(buffer.atoms[1].name().as_bytes(), b" N  ");
        assert_eq!(buffer.atoms[1].residue_id(), second);

        // Supplementary matcher-boundary evidence (BIO-COPY-A04-CLOSE):
        // this exercises the canonical chain matcher directly and wraps
        // its typed error manually — it proves the gate's typed cause and
        // source chain, NOT propagation through selection copying (A04
        // itself has no missing-name path). Real copier error propagation
        // is proven at the A06 chain-gate caller and the A08 composition.
        let bad_chain = BioChainRow::new(
            BioModelId::new(0),
            None,
            BioRowSpan::new(0, 1).unwrap(),
            ChainKind::Protein,
            ChainSourceIds::new(None, None),
        );
        let error = super::super::bio_chain_row_matches(
            &super::super::Selection::default(),
            &bad_chain,
            BioCoordinateFormat::Mmcif,
        )
        .unwrap_err();
        let wrapped = crate::selection::BioSelectionMatchError::from_traverse_error(
            crate::selection::BioRowTraverseError::Chain(error),
        );
        let chain_source = std::error::Error::source(&wrapped).unwrap();
        assert!(chain_source.is::<crate::selection::BioRowTraverseError>());
        let leaf = chain_source.source().unwrap();
        assert!(leaf.is::<crate::selection::BioRowChainError>());
    }

    #[test]
    fn bio_copy_a04_source_count_product_eight_calls() {
        use crate::hierarchy::{
            BioChainId, BioCoordinateFormat, BioEntityId, BioResidueRow, BioRowSpan,
            BioSiftsUnpResidue, EntityKind,
        };
        use crate::residue::ResidueInfoKind;
        use crate::source_ids::{PdbSeqId, ResidueName, ResidueSourceIds};

        // BIO-COPY-A04-CLOSE: the exact four-mask x SOURCE-residue-atom-count
        // (0/2) product = 8 actual append_selected_residue calls. This value
        // is the field TEMPLATE only; each iteration below reconstructs the
        // ACTUAL source row (BIO-COPY-A04-SPAN) with the real span for that
        // source count, so the empty-source case is a genuine span-0 row and
        // every mask — including both — yields an accepted empty residue.
        let source_residue = BioResidueRow::new(
            BioChainId::new(7),
            BioRowSpan::new(0, 2).unwrap(),
            ResidueName::from_ascii(b"ARG").unwrap(),
            ResidueInfoKind::Unknown,
            EntityKind::Polymer,
            Some(BioEntityId::new(9)),
            Some(b'W'),
            ResidueSourceIds::new(
                Some(PdbSeqId::new(203, Some(b'Z'))),
                Some(17),
                Some(*b"SEG9"),
                Some("SUB9".to_owned()),
                Some("E17".to_owned()),
            )
            .unwrap(),
            BioSiftsUnpResidue::new(Some(b'W'), 3, 42),
        );
        let two_atoms = [
            row(
                BioResidueId::new(7),
                b" CA ",
                Element::from_atomic_number(6).unwrap(),
                None,
                None,
                1,
                BioCalcFlag::default(),
                f64::NAN,
                -0.0,
                [1.5, -0.0, f64::NAN, 0.0, 2.5, -2.5],
                11,
                0.75,
                Some(71),
            ),
            row(
                BioResidueId::new(7),
                b" N  ",
                Element::from_atomic_number(7).unwrap(),
                None,
                None,
                -1,
                BioCalcFlag::default(),
                0.25,
                31.5,
                [0.0; 6],
                -2,
                -0.0,
                Some(72),
            ),
        ];
        let two_positions = [[-0.0, f64::NAN, 4.5], [5.5, 6.5, -0.0]];
        let empty_atoms: [BioAtomRow; 0] = [];
        let empty_positions: [[f64; 3]; 0] = [];
        let pdb = BioCoordinateFormat::Pdb;

        let mask = |names: &str| super::super::Selection {
            atom_names: super::super::SelectionList {
                all: false,
                inverted: false,
                list: names.to_owned(),
            },
            ..super::super::Selection::default()
        };
        // Four literal mask rows: expected selected SOURCE atom indices.
        let rows: [(&str, Option<&str>, &[usize]); 4] = [
            ("none", Some("O"), &[]),
            ("first", Some("CA"), &[0]),
            ("second", Some("N"), &[1]),
            ("both", None, &[0, 1]),
        ];
        let mut calls = 0usize;
        for source_count in [0usize, 2] {
            let (atoms, positions) = if source_count == 0 {
                (&empty_atoms[..], &empty_positions[..])
            } else {
                (&two_atoms[..], &two_positions[..])
            };
            // BIO-COPY-A04-SPAN: each iteration reconstructs the ACTUAL
            // source row from the template's getters with the real span —
            // start 0, length equal to both borrowed slice lengths — so the
            // empty-source case is a genuine span-0 row, not merely empty
            // slices beside a stale span(0,2) template.
            let source_residue = BioResidueRow::new(
                source_residue.chain_id(),
                BioRowSpan::new(0, source_count as u32).unwrap(),
                source_residue.name(),
                source_residue.residue_info_kind(),
                source_residue.entity_kind(),
                source_residue.entity_id(),
                source_residue.het_flag(),
                source_residue.source().clone(),
                source_residue.sifts_unp(),
            );
            assert_eq!(source_residue.atom_span().start(), 0);
            assert_eq!(source_residue.atom_span().len() as usize, atoms.len());
            assert_eq!(source_residue.atom_span().len() as usize, positions.len());
            assert_eq!(
                source_residue.sifts_unp(),
                BioSiftsUnpResidue::new(Some(b'W'), 3, 42)
            );
            for (_label, names, nominal_indices) in rows {
                // For an empty SOURCE residue every mask's selected set is
                // empty: the nominal indices index into a source that has
                // no atoms.
                let expected_indices: &[usize] = if source_count == 0 {
                    &[]
                } else {
                    nominal_indices
                };
                let selection = match names {
                    Some(names) => mask(names),
                    None => super::super::Selection::default(),
                };
                let mut buffer = CopyBuffer::default();
                let new_id = append_selected_residue(
                    &mut buffer,
                    &selection,
                    &source_residue,
                    BioChainId::new(4),
                    atoms,
                    positions,
                    pdb,
                )
                .unwrap();
                calls += 1;

                // Count proof and empty-source invariants.
                assert_eq!(buffer.residues.len(), 1);
                assert_eq!(buffer.atoms.len(), expected_indices.len());
                assert_eq!(buffer.coordinates.len(), expected_indices.len());
                let copied = &buffer.residues[0];
                assert_eq!(copied.atom_span().len() as usize, expected_indices.len());
                if source_count == 0 {
                    assert!(expected_indices.is_empty());
                    assert_eq!(copied.atom_span().start(), 0);
                }

                // Every residue getter/source field, including nondefault
                // SIFTS, in every row.
                assert_eq!(copied.chain_id(), BioChainId::new(4));
                assert_eq!(copied.name().as_bytes(), b"ARG");
                assert_eq!(copied.residue_info_kind(), ResidueInfoKind::Unknown);
                assert_eq!(copied.entity_kind(), EntityKind::Polymer);
                assert_eq!(copied.entity_id(), Some(BioEntityId::new(9)));
                assert_eq!(copied.het_flag(), Some(b'W'));
                assert_eq!(copied.source().seq_id(), source_residue.source().seq_id());
                assert_eq!(copied.source().label_seq_id(), Some(17));
                assert_eq!(copied.source().segment_id(), Some(b"SEG9"));
                assert_eq!(copied.source().subchain_id(), Some("SUB9"));
                assert_eq!(copied.source().label_entity_id(), Some("E17"));
                assert_eq!(
                    copied.sifts_unp(),
                    BioSiftsUnpResidue::new(Some(b'W'), 3, 42)
                );
                assert_ne!(copied.sifts_unp(), BioSiftsUnpResidue::default());

                // Ordered selected source atom ids and coordinate bits.
                for (position, expected_source_index) in expected_indices.iter().enumerate() {
                    let atom = &buffer.atoms[position];
                    assert_eq!(atom.residue_id(), new_id);
                    assert_eq!(
                        atom.source().serial(),
                        two_atoms[*expected_source_index].source().serial()
                    );
                    assert_eq!(
                        atom.name().as_bytes(),
                        two_atoms[*expected_source_index].name().as_bytes()
                    );
                    for (copied_component, source_component) in buffer.coordinates[position]
                        .iter()
                        .zip(two_positions[*expected_source_index])
                    {
                        assert_eq!(copied_component.to_bits(), source_component.to_bits());
                    }
                }
            }
        }
        assert_eq!(calls, 8);
    }

    #[test]
    fn bio_copy_a05_residue_masks_order_spans_and_retention() {
        use crate::hierarchy::{
            BioChainId, BioChainRow, BioCoordinateFormat, BioEntityId, BioModelId, BioResidueRow,
            BioRowSpan, BioSiftsUnpResidue, ChainKind, EntityKind,
        };
        use crate::residue::ResidueInfoKind;
        use crate::source_ids::{
            ChainSourceIds, PdbChainId, PdbSeqId, ResidueName, ResidueSourceIds,
        };

        // Two residues (ALA then GLY) with one atom each; the chain carries
        // both auth and label source names plus an entity linkage.
        let residues = [
            BioResidueRow::new(
                BioChainId::new(0),
                BioRowSpan::new(0, 1).unwrap(),
                ResidueName::from_ascii(b"ALA").unwrap(),
                ResidueInfoKind::Unknown,
                EntityKind::Polymer,
                Some(BioEntityId::new(3)),
                None,
                ResidueSourceIds::new(Some(PdbSeqId::new(101, None)), None, None, None, None)
                    .unwrap(),
                BioSiftsUnpResidue::default(),
            ),
            BioResidueRow::new(
                BioChainId::new(0),
                BioRowSpan::new(1, 1).unwrap(),
                ResidueName::from_ascii(b"GLY").unwrap(),
                ResidueInfoKind::Unknown,
                EntityKind::Polymer,
                Some(BioEntityId::new(3)),
                None,
                ResidueSourceIds::new(Some(PdbSeqId::new(102, None)), None, None, None, None)
                    .unwrap(),
                BioSiftsUnpResidue::default(),
            ),
        ];
        let atoms = [
            row(
                BioResidueId::new(0),
                b" CA ",
                Element::from_atomic_number(6).unwrap(),
                None,
                None,
                0,
                BioCalcFlag::default(),
                1.0,
                20.0,
                [0.0; 6],
                -1,
                0.0,
                Some(21),
            ),
            row(
                BioResidueId::new(1),
                b" N  ",
                Element::from_atomic_number(7).unwrap(),
                None,
                None,
                0,
                BioCalcFlag::default(),
                1.0,
                20.0,
                [0.0; 6],
                -1,
                0.0,
                Some(22),
            ),
        ];
        let positions = [[-0.0, 1.0, f64::NAN], [2.0, -0.0, 3.0]];
        let chain = BioChainRow::new(
            BioModelId::new(0),
            Some(BioEntityId::new(3)),
            BioRowSpan::new(0, 2).unwrap(),
            ChainKind::Protein,
            ChainSourceIds::new(
                Some(PdbChainId::from_ascii(b"A").unwrap()),
                Some("labelA".to_owned()),
            ),
        );
        let pdb = BioCoordinateFormat::Pdb;

        let residue_mask = |names: &str| super::super::Selection {
            residue_names: super::super::SelectionList {
                all: false,
                inverted: false,
                list: names.to_owned(),
            },
            ..super::super::Selection::default()
        };
        let atom_reject = |names: &str| super::super::Selection {
            atom_names: super::super::SelectionList {
                all: false,
                inverted: false,
                list: names.to_owned(),
            },
            ..super::super::Selection::default()
        };

        // Masks over accepted-residue counts (none/first/second/both),
        // each on a fresh empty buffer.
        for (names, expected_residues, expected_atoms) in [
            ("XXX", 0usize, 0usize),
            ("ALA", 1, 1),
            ("GLY", 1, 1),
            ("", 2, 2),
        ] {
            let selection = if names.is_empty() {
                super::super::Selection::default()
            } else {
                residue_mask(names)
            };
            let mut buffer = CopyBuffer::default();
            let new_id = append_selected_chain(
                &mut buffer,
                &selection,
                &chain,
                BioModelId::new(4),
                &residues,
                &atoms,
                &positions,
                pdb,
            )
            .unwrap();
            assert_eq!(new_id, BioChainId::new(0));
            assert_eq!(buffer.chains.len(), 1);
            assert_eq!(buffer.residues.len(), expected_residues);
            assert_eq!(buffer.atoms.len(), expected_atoms);
            let copied = &buffer.chains[0];
            assert_eq!(copied.model_id(), BioModelId::new(4));
            assert_eq!(copied.entity_id(), Some(BioEntityId::new(3)));
            assert_eq!(copied.kind(), ChainKind::Protein);
            assert_eq!(
                copied.source().auth_chain_id(),
                chain.source().auth_chain_id()
            );
            assert_eq!(
                copied.source().label_asym_id(),
                chain.source().label_asym_id()
            );
            assert_eq!(copied.residue_span().start(), 0);
            assert_eq!(copied.residue_span().len() as usize, expected_residues);
            for (index, residue) in buffer.residues.iter().enumerate() {
                assert_eq!(residue.chain_id(), new_id);
                let expected_index = match (names, expected_residues) {
                    ("ALA", 1) => 0,
                    ("GLY", 1) => 1,
                    (_, 2) => index,
                    _ => unreachable!(),
                };
                assert_eq!(
                    residue.name().as_bytes(),
                    residues[expected_index].name().as_bytes()
                );
                assert_eq!(
                    residue.source().seq_id(),
                    residues[expected_index].source().seq_id()
                );
                assert_eq!(residue.entity_id(), Some(BioEntityId::new(3)));
                assert_eq!(residue.atom_span().len(), 1);
                assert_eq!(
                    buffer.atoms[index].residue_id(),
                    BioResidueId::new(index as u32)
                );
            }
        }

        // Atom-level rejection inside an accepted residue leaves the
        // accepted EMPTY residue and the accepted chain with correct spans.
        let mut buffer = CopyBuffer::default();
        append_selected_chain(
            &mut buffer,
            &atom_reject("O"),
            &chain,
            BioModelId::new(4),
            &residues,
            &atoms,
            &positions,
            pdb,
        )
        .unwrap();
        assert_eq!(buffer.chains.len(), 1);
        assert_eq!(buffer.residues.len(), 2);
        assert!(buffer.atoms.is_empty());
        assert_eq!(buffer.residues[0].atom_span().len(), 0);
        assert_eq!(buffer.residues[1].atom_span().start(), 0);
        assert_eq!(buffer.residues[1].atom_span().len(), 0);
        assert_eq!(buffer.chains[0].residue_span().len(), 2);

        // Nonempty continuation: a second chain appends after the first,
        // with continuing chain/residue ids and preserved source order.
        let mut buffer = CopyBuffer::default();
        append_selected_chain(
            &mut buffer,
            &super::super::Selection::default(),
            &chain,
            BioModelId::new(4),
            &residues,
            &atoms,
            &positions,
            pdb,
        )
        .unwrap();
        let second = append_selected_chain(
            &mut buffer,
            &residue_mask("GLY"),
            &chain,
            BioModelId::new(5),
            &residues,
            &atoms,
            &positions,
            pdb,
        )
        .unwrap();
        assert_eq!(second, BioChainId::new(1));
        assert_eq!(buffer.chains[1].residue_span().start(), 2);
        assert_eq!(buffer.chains[1].residue_span().len(), 1);
        assert_eq!(buffer.chains[1].model_id(), BioModelId::new(5));
        assert_eq!(buffer.residues[2].name().as_bytes(), b"GLY");
        assert_eq!(buffer.residues[2].chain_id(), second);
        assert_eq!(buffer.residues[2].atom_span().start(), 2);
        assert_eq!(
            buffer.atoms[2].source().serial(),
            Some(PdbAtomSerial::new(22))
        );
        assert_eq!(buffer.coordinates[2][2].to_bits(), 3.0_f64.to_bits());
    }

    #[test]
    fn bio_copy_a06_chain_masks_empty_model_retention_and_error_chain() {
        use crate::hierarchy::{
            BioChainId, BioChainRow, BioCoordinateFormat, BioModelId, BioModelRow, BioRowSpan,
            ChainKind,
        };
        use crate::source_ids::{ChainSourceIds, PdbChainId};

        // Two distinct source chains (auth A then B) under one model with
        // source number 3; each chain carries one residue with one atom via
        // the shared two-residue/two-atom tables from the A05 test shape.
        let residues = [
            BioResidueRow::new(
                BioChainId::new(0),
                BioRowSpan::new(0, 1).unwrap(),
                crate::source_ids::ResidueName::from_ascii(b"ALA").unwrap(),
                crate::residue::ResidueInfoKind::Unknown,
                crate::hierarchy::EntityKind::Polymer,
                None,
                None,
                crate::source_ids::ResidueSourceIds::new(
                    Some(crate::source_ids::PdbSeqId::new(101, None)),
                    None,
                    None,
                    None,
                    None,
                )
                .unwrap(),
                crate::hierarchy::BioSiftsUnpResidue::default(),
            ),
            BioResidueRow::new(
                BioChainId::new(1),
                BioRowSpan::new(1, 1).unwrap(),
                crate::source_ids::ResidueName::from_ascii(b"GLY").unwrap(),
                crate::residue::ResidueInfoKind::Unknown,
                crate::hierarchy::EntityKind::Polymer,
                None,
                None,
                crate::source_ids::ResidueSourceIds::new(
                    Some(crate::source_ids::PdbSeqId::new(102, None)),
                    None,
                    None,
                    None,
                    None,
                )
                .unwrap(),
                crate::hierarchy::BioSiftsUnpResidue::default(),
            ),
        ];
        let atoms = [
            row(
                BioResidueId::new(0),
                b" CA ",
                Element::from_atomic_number(6).unwrap(),
                None,
                None,
                0,
                BioCalcFlag::default(),
                1.0,
                20.0,
                [0.0; 6],
                -1,
                0.0,
                Some(31),
            ),
            row(
                BioResidueId::new(1),
                b" N  ",
                Element::from_atomic_number(7).unwrap(),
                None,
                None,
                0,
                BioCalcFlag::default(),
                1.0,
                20.0,
                [0.0; 6],
                -1,
                0.0,
                Some(32),
            ),
        ];
        let positions = [[-0.0, 1.0, f64::NAN], [2.5, -0.0, 3.5]];
        let chains = [
            BioChainRow::new(
                BioModelId::new(0),
                None,
                BioRowSpan::new(0, 1).unwrap(),
                ChainKind::Protein,
                ChainSourceIds::new(Some(PdbChainId::from_ascii(b"A").unwrap()), None),
            ),
            BioChainRow::new(
                BioModelId::new(0),
                None,
                BioRowSpan::new(1, 1).unwrap(),
                ChainKind::WaterOnly,
                ChainSourceIds::new(Some(PdbChainId::from_ascii(b"B").unwrap()), None),
            ),
        ];
        let model = BioModelRow::new(BioRowSpan::new(0, 2).unwrap(), Some(3));
        let pdb = BioCoordinateFormat::Pdb;
        let chain_mask = |auth: &str| super::super::Selection {
            chain_ids: super::super::SelectionList {
                all: false,
                inverted: false,
                list: auth.to_owned(),
            },
            ..super::super::Selection::default()
        };

        // Chain acceptance masks over the two distinct source chains.
        for (list, expected_chains, expected_residues) in
            [("Z", 0usize, 0usize), ("A", 1, 1), ("B", 1, 1), ("", 2, 2)]
        {
            let selection = if list.is_empty() {
                super::super::Selection::default()
            } else {
                chain_mask(list)
            };
            let mut buffer = CopyBuffer::default();
            append_selected_model(
                &mut buffer,
                &selection,
                &model,
                &chains,
                &residues,
                &atoms,
                &positions,
                pdb,
            )
            .unwrap();
            assert_eq!(buffer.models.len(), 1);
            assert_eq!(buffer.models[0].source_model_number(), Some(3));
            assert_eq!(buffer.models[0].chain_span().start(), 0);
            assert_eq!(
                buffer.models[0].chain_span().len() as usize,
                expected_chains
            );
            assert_eq!(buffer.chains.len(), expected_chains);
            assert_eq!(buffer.residues.len(), expected_residues);
            assert_eq!(buffer.atoms.len(), expected_residues);
            for (index, chain) in buffer.chains.iter().enumerate() {
                assert_eq!(chain.model_id(), BioModelId::new(0));
                let expected = match (list, expected_chains) {
                    ("A", 1) => &chains[0],
                    ("B", 1) => &chains[1],
                    (_, 2) => &chains[index],
                    _ => unreachable!(),
                };
                assert_eq!(
                    chain.source().auth_chain_id(),
                    expected.source().auth_chain_id()
                );
                assert_eq!(chain.kind(), expected.kind());
                assert_eq!(chain.residue_span().len(), 1);
            }
        }

        // Accepted originally EMPTY model (span 0, zero chains): emitted
        // with remapped empty span and preserved source number.
        let empty_model = BioModelRow::new(BioRowSpan::new(0, 0).unwrap(), Some(9));
        let mut buffer = CopyBuffer::default();
        append_selected_model(
            &mut buffer,
            &super::super::Selection::default(),
            &empty_model,
            &chains,
            &residues,
            &atoms,
            &positions,
            pdb,
        )
        .unwrap();
        assert_eq!(buffer.models.len(), 1);
        assert_eq!(buffer.models[0].source_model_number(), Some(9));
        assert_eq!(buffer.models[0].chain_span().start(), 0);
        assert_eq!(buffer.models[0].chain_span().len(), 0);
        assert!(buffer.chains.is_empty());

        // REAL copier error propagation: a chain missing both canonical
        // names fails INSIDE append_selected_model with the typed cause and
        // the standard source chain (no manual wrapping). The model row
        // carries the single-chain span matching this one-row table.
        let bad_chains = [BioChainRow::new(
            BioModelId::new(0),
            None,
            BioRowSpan::new(0, 1).unwrap(),
            ChainKind::Protein,
            ChainSourceIds::new(None, None),
        )];
        let bad_model = BioModelRow::new(BioRowSpan::new(0, 1).unwrap(), Some(3));
        let error = append_selected_model(
            &mut CopyBuffer::default(),
            &super::super::Selection::default(),
            &bad_model,
            &bad_chains,
            &residues,
            &atoms,
            &positions,
            BioCoordinateFormat::Mmcif,
        )
        .unwrap_err();
        assert!(matches!(
            error.cause,
            SelectionCopyCause::Traverse(crate::selection::BioRowTraverseError::Chain(
                crate::selection::BioRowChainError::MissingCanonicalChainName
            ))
        ));
        let source = std::error::Error::source(&error).unwrap();
        assert!(source.is::<crate::selection::BioRowTraverseError>());
        let leaf = source.source().unwrap();
        assert!(leaf.is::<crate::selection::BioRowChainError>());
    }

    #[test]
    fn bio_copy_a07_all16_defaulted_field_combinations() {
        use crate::structure_metadata::BioStructureSourceState;
        use std::collections::BTreeMap;

        // BIO-COPY A07: all 16 combinations of the four source-defaulted
        // fields (conect_map nonempty/empty x has_d_fraction true/false x
        // non_ascii_line nonzero/0 x ter_status nonzero/NUL). Every
        // source-copied field stays nondefault in every combination; the
        // four unassigned members keep fresh-structure defaults — Gemmi
        // empty_copy behavior, not data-loss cleanup.
        let base = || {
            let mut info = BTreeMap::new();
            info.insert("_key".to_owned(), "value".to_owned());
            BioStructureSourceState {
                name: "source".to_owned(),
                resolution: f64::NAN,
                conect_map: BTreeMap::new(),
                has_d_fraction: false,
                non_ascii_line: 0,
                ter_status: 0,
                has_origx: true,
                origx: crate::hierarchy::BioTransform::new(
                    [
                        [-0.0, 1.5, f64::NAN],
                        [2.5, -0.0, 3.5],
                        [f64::INFINITY, -1.0, -0.0],
                    ],
                    [-0.0, f64::NAN, 4.5],
                ),
                info,
                raw_remarks: vec!["r1".to_owned(), "r1".to_owned(), "r2".to_owned()],
            }
        };
        let mut cases = 0usize;
        for conect_nonempty in [false, true] {
            for d_fraction in [false, true] {
                for non_ascii in [0, 42] {
                    for ter in [0u8, b'T'] {
                        let mut source = base();
                        if conect_nonempty {
                            source.conect_map.insert(7, vec![8, 8]);
                        }
                        source.has_d_fraction = d_fraction;
                        source.non_ascii_line = non_ascii;
                        source.ter_status = ter;
                        let copied = empty_copy_source_state(&source);
                        cases += 1;

                        // Copied fields: exact values, order, duplicates, bits.
                        assert_eq!(copied.name, "source");
                        assert_eq!(copied.resolution.to_bits(), f64::NAN.to_bits());
                        assert!(copied.has_origx);
                        for (copied_row, source_row) in
                            copied.origx.matrix().iter().zip(source.origx.matrix())
                        {
                            for (copied_component, source_component) in
                                copied_row.iter().zip(source_row)
                            {
                                assert_eq!(copied_component.to_bits(), source_component.to_bits());
                            }
                        }
                        for (copied_component, source_component) in copied
                            .origx
                            .translation()
                            .iter()
                            .zip(source.origx.translation())
                        {
                            assert_eq!(copied_component.to_bits(), source_component.to_bits());
                        }
                        assert_eq!(copied.info.get("_key").map(String::as_str), Some("value"));
                        assert_eq!(copied.raw_remarks, source.raw_remarks);

                        // Defaulted fields: fresh-structure member state.
                        assert!(copied.conect_map.is_empty());
                        assert!(!copied.has_d_fraction);
                        assert_eq!(copied.non_ascii_line, 0);
                        assert_eq!(copied.ter_status, 0);

                        // Source immutability.
                        assert_eq!(source.conect_map.len(), usize::from(conect_nonempty));
                        assert_eq!(source.has_d_fraction, d_fraction);
                        assert_eq!(source.non_ascii_line, non_ascii);
                        assert_eq!(source.ter_status, ter);
                        assert_eq!(source.raw_remarks.len(), 3);
                    }
                }
            }
        }
        assert_eq!(cases, 16);
    }

    #[test]
    fn bio_copy_a08_composed_selections_metadata_and_error_chain() {
        use crate::hierarchy::{
            BioAssembly, BioAssemblyGenerator, BioAssemblySpecialKind, BioChainRow,
            BioCoordinateBlock, BioCoordinateFormat, BioEntityRow, BioModelRow, BioResidueRow,
            BioRowSpan, BioSiftsUnpResidue, BioStructureData, BioStructureParts, ChainKind,
            EntityKind, PolymerKind,
        };
        use crate::residue::ResidueInfoKind;
        use crate::source_ids::{
            ChainSourceIds, EntitySourceIds, PdbChainId, PdbSeqId, ResidueName, ResidueSourceIds,
        };
        use crate::structure_metadata::BioStructureSourceState;

        // Two-model fixture: model 1 = chain A/ALA/ CA (serial 41),
        // model 2 = chain B/GLY/ N  (serial 42), one entity covering both
        // subchains, one assembly with source names, nondefault source_state.
        let mut source_state = BioStructureSourceState::default();
        source_state.name = "src".to_owned();
        source_state.resolution = f64::NAN;
        source_state.conect_map.insert(41, vec![42]);
        source_state.has_d_fraction = true;
        source_state.non_ascii_line = 7;
        source_state.ter_status = b'X';
        source_state.has_origx = true;
        source_state.raw_remarks = vec!["r1".to_owned()];
        let empty_parts = || BioStructureParts {
            input_format: BioCoordinateFormat::Pdb,
            models: Vec::new(),
            chains: Vec::new(),
            residues: Vec::new(),
            atoms: Vec::new(),
            entities: Vec::new(),
            connections: Vec::new(),
            cispeps: Vec::new(),
            mod_residues: Vec::new(),
            helices: Vec::new(),
            sheets: Vec::new(),
            metadata: Default::default(),
            source_state: BioStructureSourceState::default(),
            coordinates: BioCoordinateBlock::new(Vec::new()),
            crystal: None,
            ncs_operators: Vec::new(),
            assemblies: Vec::new(),
        };
        let parts = BioStructureParts {
            models: vec![
                BioModelRow::new(BioRowSpan::new(0, 1).unwrap(), Some(1)),
                BioModelRow::new(BioRowSpan::new(1, 1).unwrap(), Some(2)),
            ],
            chains: vec![
                BioChainRow::new(
                    crate::hierarchy::BioModelId::new(0),
                    None,
                    BioRowSpan::new(0, 1).unwrap(),
                    ChainKind::Protein,
                    ChainSourceIds::new(Some(PdbChainId::from_ascii(b"A").unwrap()), None),
                ),
                BioChainRow::new(
                    crate::hierarchy::BioModelId::new(1),
                    None,
                    BioRowSpan::new(1, 1).unwrap(),
                    ChainKind::Protein,
                    ChainSourceIds::new(Some(PdbChainId::from_ascii(b"B").unwrap()), None),
                ),
            ],
            residues: vec![
                BioResidueRow::new(
                    crate::hierarchy::BioChainId::new(0),
                    BioRowSpan::new(0, 1).unwrap(),
                    ResidueName::from_ascii(b"ALA").unwrap(),
                    ResidueInfoKind::Unknown,
                    EntityKind::Polymer,
                    None,
                    None,
                    ResidueSourceIds::new(
                        Some(PdbSeqId::new(101, None)),
                        None,
                        None,
                        Some("S1".to_owned()),
                        None,
                    )
                    .unwrap(),
                    BioSiftsUnpResidue::default(),
                ),
                BioResidueRow::new(
                    crate::hierarchy::BioChainId::new(1),
                    BioRowSpan::new(1, 1).unwrap(),
                    ResidueName::from_ascii(b"GLY").unwrap(),
                    ResidueInfoKind::Unknown,
                    EntityKind::Polymer,
                    None,
                    None,
                    ResidueSourceIds::new(
                        Some(PdbSeqId::new(202, None)),
                        None,
                        None,
                        Some("S2".to_owned()),
                        None,
                    )
                    .unwrap(),
                    BioSiftsUnpResidue::default(),
                ),
            ],
            atoms: vec![
                row(
                    BioResidueId::new(0),
                    b" CA ",
                    Element::from_atomic_number(6).unwrap(),
                    None,
                    None,
                    0,
                    BioCalcFlag::default(),
                    1.0,
                    20.0,
                    [0.0; 6],
                    -1,
                    0.0,
                    Some(41),
                ),
                row(
                    BioResidueId::new(1),
                    b" N  ",
                    Element::from_atomic_number(7).unwrap(),
                    None,
                    None,
                    0,
                    BioCalcFlag::default(),
                    0.5,
                    21.0,
                    [0.0; 6],
                    -1,
                    0.0,
                    Some(42),
                ),
            ],
            entities: vec![BioEntityRow::new(
                EntityKind::Polymer,
                PolymerKind::PeptideL,
                false,
                Vec::new(),
                Vec::new(),
                Vec::new(),
                vec!["S1".to_owned(), "S2".to_owned()],
                EntitySourceIds::new("1".to_owned()),
            )],
            coordinates: BioCoordinateBlock::new(vec![[-0.0, 1.0, f64::NAN], [2.0, -0.0, 3.0]]),
            assemblies: vec![BioAssembly::new(
                "1".to_owned(),
                true,
                false,
                BioAssemblySpecialKind::NotApplicable,
                1,
                String::new(),
                String::new(),
                0.0,
                0.0,
                0.0,
                vec![BioAssemblyGenerator::new(
                    vec!["A".to_owned(), "B".to_owned()],
                    vec!["S1".to_owned()],
                    Vec::new(),
                )],
            )],
            source_state,
            ..empty_parts()
        };
        let source = BioStructureData::from_parts(parts.clone()).unwrap();
        let source_serials: Vec<_> = source.atoms.iter().map(|a| a.source().serial()).collect();

        let with = |modify: &dyn Fn(super::super::Selection) -> super::super::Selection| {
            crate::selection::BioSelectionData {
                selection: modify(super::super::Selection::default()),
            }
        };
        let by_chain = |v: &str| {
            with(&|s| super::super::Selection {
                chain_ids: super::super::SelectionList {
                    all: false,
                    inverted: false,
                    list: v.to_owned(),
                },
                ..s
            })
        };
        let by_residue = |v: &str| {
            with(&|s| super::super::Selection {
                residue_names: super::super::SelectionList {
                    all: false,
                    inverted: false,
                    list: v.to_owned(),
                },
                ..s
            })
        };
        let by_atom = |v: &str| {
            with(&|s| super::super::Selection {
                atom_names: super::super::SelectionList {
                    all: false,
                    inverted: false,
                    list: v.to_owned(),
                },
                ..s
            })
        };
        let by_model = |v: i32| with(&|s| super::super::Selection { mdl: v, ..s });

        // Seven selections with literal expected shapes:
        // (models, chains, residues, atoms, copied-position -> original atom index).
        let cases: Vec<(
            &str,
            crate::selection::BioSelectionData,
            [usize; 4],
            Vec<usize>,
        )> = vec![
            ("default-all", with(&|s| s), [2, 2, 2, 2], vec![0, 1]),
            ("no-model", by_model(99), [0, 0, 0, 0], vec![]),
            ("no-chain", by_chain("Z"), [2, 0, 0, 0], vec![]),
            ("no-residue", by_residue("XXX"), [2, 2, 0, 0], vec![]),
            ("no-atom", by_atom("O"), [2, 2, 2, 0], vec![]),
            ("one-chain", by_chain("A"), [2, 1, 1, 1], vec![0]),
            ("one-residue", by_residue("GLY"), [2, 2, 1, 1], vec![1]),
        ];
        for (label, selection, expected, expected_source_atoms) in cases {
            let copied = copy_bio_selection(&source, &selection).unwrap();
            assert_eq!(copied.models.len(), expected[0], "{label}");
            assert_eq!(copied.chains.len(), expected[1], "{label}");
            assert_eq!(copied.residues.len(), expected[2], "{label}");
            assert_eq!(copied.atoms.len(), expected[3], "{label}");
            assert_eq!(copied.coordinates.positions().len(), expected[3], "{label}");
            // Exact remapped spans/parents/source numbers per level.
            for (index, model) in copied.models.iter().enumerate() {
                let expected_span_len = if expected[1] == 2 {
                    1
                } else if expected[1] == 1 {
                    usize::from(index == 0)
                } else {
                    0
                };
                assert_eq!(
                    model.chain_span().len() as usize,
                    expected_span_len,
                    "{label} model {index}"
                );
                assert_eq!(model.source_model_number(), Some(index as i32 + 1));
            }
            for (index, chain) in copied.chains.iter().enumerate() {
                assert_eq!(
                    chain.model_id(),
                    crate::hierarchy::BioModelId::new(index as u32),
                    "{label}"
                );
                let expected_chain_residues = match label {
                    "no-residue" => 0,
                    "one-residue" => usize::from(index == 1),
                    _ => 1,
                };
                assert_eq!(
                    chain.residue_span().len() as usize,
                    expected_chain_residues,
                    "{label} chain {index}"
                );
            }
            for (index, residue) in copied.residues.iter().enumerate() {
                assert_eq!(
                    residue.atom_span().len() as usize,
                    if expected[3] == 0 { 0 } else { 1 },
                    "{label} residue {index}"
                );
                let expected_chain_id = if label == "one-residue" { 1 } else { index };
                assert_eq!(
                    residue.chain_id(),
                    crate::hierarchy::BioChainId::new(expected_chain_id as u32),
                    "{label}"
                );
            }
            // Positions map by original atom table index (bit-exact).
            for (copied_position, original_index) in copied
                .coordinates
                .positions()
                .iter()
                .zip(&expected_source_atoms)
            {
                for (copied_component, source_component) in copied_position
                    .iter()
                    .zip(source.coordinates.positions()[*original_index])
                {
                    assert_eq!(copied_component.to_bits(), source_component.to_bits());
                }
            }
            if expected[3] == 1 {
                let original = expected_source_atoms[0];
                assert_eq!(
                    copied.atoms[0].source().serial(),
                    source.atoms[original].source().serial(),
                    "{label}"
                );
                assert_eq!(
                    copied.atoms[0].residue_id(),
                    BioResidueId::new(0),
                    "{label}"
                );
            }
            // All retained metadata categories stay exact and SHARED.
            assert!(std::sync::Arc::ptr_eq(&copied.entities, &source.entities));
            assert!(std::sync::Arc::ptr_eq(
                &copied.connections,
                &source.connections
            ));
            assert!(std::sync::Arc::ptr_eq(&copied.cispeps, &source.cispeps));
            assert!(std::sync::Arc::ptr_eq(
                &copied.mod_residues,
                &source.mod_residues
            ));
            assert!(std::sync::Arc::ptr_eq(&copied.helices, &source.helices));
            assert!(std::sync::Arc::ptr_eq(&copied.sheets, &source.sheets));
            assert!(std::sync::Arc::ptr_eq(&copied.metadata, &source.metadata));
            assert!(std::sync::Arc::ptr_eq(&copied.crystal, &source.crystal));
            assert!(std::sync::Arc::ptr_eq(
                &copied.ncs_operators,
                &source.ncs_operators
            ));
            assert!(std::sync::Arc::ptr_eq(
                &copied.assemblies,
                &source.assemblies
            ));
            assert_eq!(copied.input_format, BioCoordinateFormat::Pdb);
            assert_eq!(copied.assemblies[0].generators[0].chains, vec!["A", "B"]);
            // A07 source-state semantics inside the composed result.
            assert_eq!(copied.source_state.name, "src");
            assert_eq!(copied.source_state.resolution.to_bits(), f64::NAN.to_bits());
            assert!(copied.source_state.has_origx);
            assert_eq!(copied.source_state.raw_remarks, vec!["r1"]);
            assert!(copied.source_state.conect_map.is_empty());
            assert!(!copied.source_state.has_d_fraction);
            assert_eq!(copied.source_state.non_ascii_line, 0);
            assert_eq!(copied.source_state.ter_status, 0);
            // The composed output is a validated CK hierarchy.
            copied.validate().unwrap();
        }

        // Original-source stability: nothing above mutated the source.
        assert_eq!(
            source_serials,
            vec![Some(PdbAtomSerial::new(41)), Some(PdbAtomSerial::new(42))]
        );
        assert_eq!(source.models.len(), 2);
        assert_eq!(source.source_state.conect_map.len(), 1);
        assert!(source.source_state.has_d_fraction);

        // REAL composed error chain, model leaf: a validated source whose
        // stored model number is missing fails inside copy_bio_selection.
        let none_parts = BioStructureParts {
            models: vec![BioModelRow::new(BioRowSpan::new(0, 0).unwrap(), None)],
            ..empty_parts()
        };
        let none_source =
            BioStructureData::from_parts(none_parts).expect("none-model fixture must construct");
        {
            let error = copy_bio_selection(&none_source, &with(&|s| s)).unwrap_err();
            assert!(matches!(
                error.cause,
                SelectionCopyCause::Traverse(crate::selection::BioRowTraverseError::Model(
                    crate::selection::BioRowModelError::MissingModelNumber
                ))
            ));
            let source = std::error::Error::source(&error).unwrap();
            assert!(source.is::<crate::selection::BioRowTraverseError>());
            let leaf = source.source().unwrap();
            assert!(leaf.is::<crate::selection::BioRowModelError>());
        }
        // REAL composed error chain, chain leaf (same chain as A06).
        let bad_parts = BioStructureParts {
            input_format: BioCoordinateFormat::Mmcif,
            models: vec![BioModelRow::new(BioRowSpan::new(0, 1).unwrap(), Some(3))],
            chains: vec![BioChainRow::new(
                crate::hierarchy::BioModelId::new(0),
                None,
                BioRowSpan::new(0, 0).unwrap(),
                ChainKind::Protein,
                ChainSourceIds::new(None, None),
            )],
            ..empty_parts()
        };
        let bad_source =
            BioStructureData::from_parts(bad_parts).expect("bad-chain fixture must construct");
        {
            let error = copy_bio_selection(&bad_source, &with(&|s| s)).unwrap_err();
            assert!(matches!(
                error.cause,
                SelectionCopyCause::Traverse(crate::selection::BioRowTraverseError::Chain(
                    crate::selection::BioRowChainError::MissingCanonicalChainName
                ))
            ));
            let source = std::error::Error::source(&error).unwrap();
            assert!(source.is::<crate::selection::BioRowTraverseError>());
            let leaf = source.source().unwrap();
            assert!(leaf.is::<crate::selection::BioRowChainError>());
        }

        // BIO-COPY-A08-INPUT / A08-PROOF: the ACTUAL valid composition
        // fixture, cloned, with ONLY its models block's first model chain
        // span changed to (0,5) beyond the 2-row chain table — the second
        // model and every other source block (chains, entities, assembly,
        // metadata, source_state, coordinates) are the original fixture's.
        // The malformed BioStructureData is constructed directly because
        // from_parts rejects it. The REAL copy_bio_selection must return the
        // exact typed RowSpanOutOfBounds payload under BOTH the default-all
        // and no-model-match selectors — the input check precedes traversal,
        // so a non-matching selection cannot mask it. No panic catching, no
        // manual wrapping, no default output.
        let mut corrupted_parts = parts;
        corrupted_parts.models[0] = BioModelRow::new(BioRowSpan::new(0, 5).unwrap(), Some(1));
        let corrupted = BioStructureData {
            input_format: corrupted_parts.input_format,
            models: std::sync::Arc::new(corrupted_parts.models),
            chains: std::sync::Arc::new(corrupted_parts.chains),
            residues: std::sync::Arc::new(corrupted_parts.residues),
            atoms: std::sync::Arc::new(corrupted_parts.atoms),
            entities: std::sync::Arc::new(corrupted_parts.entities),
            connections: std::sync::Arc::new(corrupted_parts.connections),
            cispeps: std::sync::Arc::new(corrupted_parts.cispeps),
            mod_residues: std::sync::Arc::new(corrupted_parts.mod_residues),
            helices: std::sync::Arc::new(corrupted_parts.helices),
            sheets: std::sync::Arc::new(corrupted_parts.sheets),
            metadata: std::sync::Arc::new(corrupted_parts.metadata),
            source_state: std::sync::Arc::new(corrupted_parts.source_state),
            coordinates: std::sync::Arc::new(corrupted_parts.coordinates),
            crystal: std::sync::Arc::new(corrupted_parts.crystal),
            ncs_operators: std::sync::Arc::new(corrupted_parts.ncs_operators),
            assemblies: std::sync::Arc::new(corrupted_parts.assemblies),
        };
        for (label, selection) in [("default-all", with(&|s| s)), ("no-model", by_model(99))] {
            let error = copy_bio_selection(&corrupted, &selection).unwrap_err();
            assert!(
                matches!(
                    &error.cause,
                    SelectionCopyCause::Structure(BioStructureError::RowSpanOutOfBounds {
                        start: 0,
                        len: 5,
                        table_len: 2,
                    })
                ),
                "{label}"
            );
            let source = std::error::Error::source(&error).unwrap();
            assert!(source.is::<BioStructureError>(), "{label}");
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn a09_parts(
        auth: Option<&[u8]>,
        label: Option<&str>,
        auth2: Option<&[u8]>,
        label2: Option<&str>,
    ) -> crate::hierarchy::BioStructureParts {
        use crate::hierarchy::{
            BioAssembly, BioAssemblyGenerator, BioAssemblySpecialKind, BioChainRow,
            BioCoordinateBlock, BioCoordinateFormat, BioEntityRow, BioModelRow, BioResidueRow,
            BioRowSpan, BioSiftsUnpResidue, BioStructureParts, ChainKind, EntityKind, PolymerKind,
        };
        use crate::source_ids::{
            ChainSourceIds, EntitySourceIds, PdbChainId, PdbSeqId, ResidueName, ResidueSourceIds,
        };
        use crate::structure_metadata::BioStructureSourceState;
        let mut source_state = BioStructureSourceState::default();
        source_state.name = "meta-src".to_owned();
        source_state.resolution = f64::NAN;
        source_state.conect_map.insert(41, vec![42]);
        source_state.has_d_fraction = true;
        source_state.non_ascii_line = 9;
        source_state.ter_status = b'Z';
        source_state.has_origx = true;
        source_state.raw_remarks = vec!["rem".to_owned()];
        let chain_ids = |auth: Option<&[u8]>, label: Option<&str>| {
            ChainSourceIds::new(
                auth.map(|a| PdbChainId::from_ascii(a).unwrap()),
                label.map(|l| l.to_owned()),
            )
        };
        BioStructureParts {
            input_format: BioCoordinateFormat::Pdb,
            models: vec![
                BioModelRow::new(BioRowSpan::new(0, 1).unwrap(), Some(1)),
                BioModelRow::new(BioRowSpan::new(1, 1).unwrap(), Some(2)),
            ],
            chains: vec![
                BioChainRow::new(
                    crate::hierarchy::BioModelId::new(0),
                    None,
                    BioRowSpan::new(0, 1).unwrap(),
                    ChainKind::Protein,
                    chain_ids(auth, label),
                ),
                BioChainRow::new(
                    crate::hierarchy::BioModelId::new(1),
                    None,
                    BioRowSpan::new(1, 1).unwrap(),
                    ChainKind::WaterOnly,
                    chain_ids(auth2, label2),
                ),
            ],
            residues: vec![
                BioResidueRow::new(
                    crate::hierarchy::BioChainId::new(0),
                    BioRowSpan::new(0, 1).unwrap(),
                    ResidueName::from_ascii(b"ALA").unwrap(),
                    crate::residue::ResidueInfoKind::Unknown,
                    EntityKind::Polymer,
                    None,
                    None,
                    ResidueSourceIds::new(
                        Some(PdbSeqId::new(101, None)),
                        None,
                        None,
                        Some("S1".to_owned()),
                        None,
                    )
                    .unwrap(),
                    BioSiftsUnpResidue::default(),
                ),
                BioResidueRow::new(
                    crate::hierarchy::BioChainId::new(1),
                    BioRowSpan::new(1, 1).unwrap(),
                    ResidueName::from_ascii(b"GLY").unwrap(),
                    crate::residue::ResidueInfoKind::Unknown,
                    EntityKind::Polymer,
                    None,
                    None,
                    ResidueSourceIds::new(
                        Some(PdbSeqId::new(202, None)),
                        None,
                        None,
                        Some("S2".to_owned()),
                        None,
                    )
                    .unwrap(),
                    BioSiftsUnpResidue::default(),
                ),
            ],
            atoms: vec![
                row(
                    BioResidueId::new(0),
                    b" CA ",
                    Element::from_atomic_number(6).unwrap(),
                    None,
                    None,
                    0,
                    BioCalcFlag::default(),
                    f64::NAN,
                    -0.0,
                    [1.5, -0.0, f64::NAN, 0.0, 2.5, -2.5],
                    4,
                    0.5,
                    Some(41),
                ),
                row(
                    BioResidueId::new(1),
                    b" N  ",
                    Element::from_atomic_number(7).unwrap(),
                    None,
                    None,
                    0,
                    BioCalcFlag::default(),
                    1.0,
                    30.0,
                    [0.0; 6],
                    -1,
                    0.0,
                    Some(42),
                ),
            ],
            // The entity references an UNSELECTED subchain alongside the
            // two used ones — preserved verbatim, never pruned.
            entities: vec![BioEntityRow::new(
                EntityKind::Polymer,
                PolymerKind::PeptideL,
                false,
                Vec::new(),
                Vec::new(),
                Vec::new(),
                vec![
                    "S1".to_owned(),
                    "S2".to_owned(),
                    "NEVER_SELECTED".to_owned(),
                ],
                EntitySourceIds::new("1".to_owned()),
            )],
            assemblies: vec![BioAssembly::new(
                "asm".to_owned(),
                true,
                false,
                BioAssemblySpecialKind::NotApplicable,
                2,
                "od".to_owned(),
                "sw".to_owned(),
                -0.0,
                f64::NAN,
                f64::MIN_POSITIVE,
                vec![BioAssemblyGenerator::new(
                    vec!["CH1".to_owned(), "CH2".to_owned()],
                    vec!["S1".to_owned()],
                    Vec::new(),
                )],
            )],
            source_state,
            coordinates: BioCoordinateBlock::new(vec![[-0.0, 1.0, f64::NAN], [2.0, -0.0, 3.0]]),
            ..{
                let empty = BioStructureParts {
                    input_format: BioCoordinateFormat::Pdb,
                    models: Vec::new(),
                    chains: Vec::new(),
                    residues: Vec::new(),
                    atoms: Vec::new(),
                    entities: Vec::new(),
                    connections: Vec::new(),
                    cispeps: Vec::new(),
                    mod_residues: Vec::new(),
                    helices: Vec::new(),
                    sheets: Vec::new(),
                    metadata: Default::default(),
                    source_state: BioStructureSourceState::default(),
                    coordinates: BioCoordinateBlock::new(Vec::new()),
                    crystal: None,
                    ncs_operators: Vec::new(),
                    assemblies: Vec::new(),
                };
                empty
            }
        }
    }

    #[test]
    fn bio_copy_a09_selector_reference_shapes_28() {
        use crate::hierarchy::{
            BioAssembly, BioAssemblyGenerator, BioAssemblySpecialKind, BioChainRow,
            BioCoordinateBlock, BioCoordinateFormat, BioEntityRow, BioModelRow, BioResidueRow,
            BioRowSpan, BioSiftsUnpResidue, BioStructureData, BioStructureParts, ChainKind,
            EntityKind, PolymerKind,
        };
        use crate::source_ids::{
            ChainSourceIds, EntitySourceIds, PdbChainId, PdbSeqId, ResidueName, ResidueSourceIds,
        };
        use crate::structure_metadata::BioStructureSourceState;

        // Four source-reference states: (chain-1 identity, chain-2
        // identity, the name a chain-targeting selection resolves to).
        // resolved-auth: auth "A" selected by "A"; unresolved-auth:
        // auth "A" while the selection references "Z"; resolved-subchain:
        // no auth, label "LB" selected by "LB"; unresolved-subchain: no
        // auth, label "LB" while the selection references "ZZ".
        type State<'a> = (
            Option<&'a [u8]>,
            Option<&'a str>,
            Option<&'a [u8]>,
            Option<&'a str>,
            &'a str,
            bool,
        );
        let states: [State; 4] = [
            (Some(b"A"), None, Some(b"A2"), None, "A", true),
            (Some(b"A"), None, Some(b"A2"), None, "Z", false),
            (None, Some("LB"), None, Some("LB2"), "LB", true),
            (None, Some("LB"), None, Some("LB2"), "ZZ", false),
        ];
        let build = a09_parts;
        let with = |modify: &dyn Fn(super::super::Selection) -> super::super::Selection| {
            crate::selection::BioSelectionData {
                selection: modify(super::super::Selection::default()),
            }
        };
        let mut calls = 0usize;
        for (auth, label, auth2, label2, chain_ref, resolved) in states {
            let source = BioStructureData::from_parts(build(auth, label, auth2, label2))
                .expect("state fixture");
            let by_chain = |v: &str| {
                with(&|s| super::super::Selection {
                    chain_ids: super::super::SelectionList {
                        all: false,
                        inverted: false,
                        list: v.to_owned(),
                    },
                    ..s
                })
            };
            let by_residue = |v: &str| {
                with(&|s| super::super::Selection {
                    residue_names: super::super::SelectionList {
                        all: false,
                        inverted: false,
                        list: v.to_owned(),
                    },
                    ..s
                })
            };
            let by_atom = |v: &str| {
                with(&|s| super::super::Selection {
                    atom_names: super::super::SelectionList {
                        all: false,
                        inverted: false,
                        list: v.to_owned(),
                    },
                    ..s
                })
            };
            // The seven A08 selection profiles. `one-chain` targets the
            // state's chain reference, so unresolved states select nothing
            // there while every other profile keeps its A08 meaning.
            let one_chain_expected: [usize; 4] = if resolved { [2, 1, 1, 1] } else { [2, 0, 0, 0] };
            let profiles: Vec<(
                &str,
                crate::selection::BioSelectionData,
                [usize; 4],
                Vec<usize>,
            )> = vec![
                ("default-all", with(&|s| s), [2, 2, 2, 2], vec![0, 1]),
                (
                    "no-model",
                    with(&|s| super::super::Selection { mdl: 99, ..s }),
                    [0, 0, 0, 0],
                    vec![],
                ),
                ("no-chain", by_chain("NOPE"), [2, 0, 0, 0], vec![]),
                ("no-residue", by_residue("XXX"), [2, 2, 0, 0], vec![]),
                ("no-atom", by_atom("O"), [2, 2, 2, 0], vec![]),
                ("one-chain", by_chain(chain_ref), one_chain_expected, {
                    if resolved { vec![0] } else { vec![] }
                }),
                ("one-residue", by_residue("GLY"), [2, 2, 1, 1], vec![1]),
            ];
            for (label_p, selection, expected, expected_source_atoms) in profiles {
                let copied = copy_bio_selection(&source, &selection).unwrap();
                calls += 1;
                assert_eq!(copied.models.len(), expected[0], "{label_p}");
                assert_eq!(copied.chains.len(), expected[1], "{label_p}");
                assert_eq!(copied.residues.len(), expected[2], "{label_p}");
                assert_eq!(copied.atoms.len(), expected[3], "{label_p}");
                assert_eq!(copied.coordinates.positions().len(), expected[3]);
                for (index, model) in copied.models.iter().enumerate() {
                    let expected_span = if expected[1] == 2 {
                        1
                    } else if expected[1] == 1 {
                        usize::from(index == 0)
                    } else {
                        0
                    };
                    assert_eq!(model.chain_span().len() as usize, expected_span);
                    assert_eq!(model.source_model_number(), Some(index as i32 + 1));
                }
                for (index, chain) in copied.chains.iter().enumerate() {
                    assert_eq!(
                        chain.model_id(),
                        crate::hierarchy::BioModelId::new(index as u32)
                    );
                    let expected_chain_residues = match label_p {
                        "no-residue" => 0,
                        "one-residue" => usize::from(index == 1),
                        _ => 1,
                    };
                    assert_eq!(chain.residue_span().len() as usize, expected_chain_residues);
                    // Raw source identity preserved verbatim (no rewrite).
                    assert_eq!(
                        chain.source().auth_chain_id(),
                        source.chains[index].source().auth_chain_id()
                    );
                    assert_eq!(
                        chain.source().label_asym_id(),
                        source.chains[index].source().label_asym_id()
                    );
                }
                for (index, residue) in copied.residues.iter().enumerate() {
                    assert_eq!(
                        residue.atom_span().len() as usize,
                        if expected[3] == 0 { 0 } else { 1 }
                    );
                    let expected_chain_id = if label_p == "one-residue" { 1 } else { index };
                    assert_eq!(
                        residue.chain_id(),
                        crate::hierarchy::BioChainId::new(expected_chain_id as u32)
                    );
                }
                for (copied_position, original_index) in copied
                    .coordinates
                    .positions()
                    .iter()
                    .zip(&expected_source_atoms)
                {
                    for (copied_component, source_component) in copied_position
                        .iter()
                        .zip(source.coordinates.positions()[*original_index])
                    {
                        assert_eq!(copied_component.to_bits(), source_component.to_bits());
                    }
                }
                if expected[3] == 1 {
                    let original = expected_source_atoms[0];
                    assert_eq!(
                        copied.atoms[0].source().serial(),
                        source.atoms[original].source().serial()
                    );
                    if original == 0 {
                        assert_eq!(copied.atoms[0].occupancy().to_bits(), f64::NAN.to_bits());
                        assert_eq!(copied.atoms[0].b_iso().to_bits(), (-0.0_f64).to_bits());
                    }
                }
                // Exact preserved raw metadata/float bits and sharing.
                assert!(std::sync::Arc::ptr_eq(&copied.entities, &source.entities));
                assert!(std::sync::Arc::ptr_eq(
                    &copied.assemblies,
                    &source.assemblies
                ));
                assert!(std::sync::Arc::ptr_eq(&copied.metadata, &source.metadata));
                let assembly = &copied.assemblies[0];
                assert_eq!(assembly.name, "asm");
                assert_eq!(assembly.generators[0].chains, vec!["CH1", "CH2"]);
                assert_eq!(assembly.oligomeric_details, "od");
                assert_eq!(assembly.software_name, "sw");
                assert_eq!(assembly.buried_surface_area.to_bits(), (-0.0_f64).to_bits());
                assert!(assembly.surface_area.is_nan());
                assert_eq!(
                    assembly.solvent_free_energy_change.to_bits(),
                    f64::MIN_POSITIVE.to_bits()
                );
                assert_eq!(
                    copied.entities[0].subchains(),
                    &[
                        "S1".to_owned(),
                        "S2".to_owned(),
                        "NEVER_SELECTED".to_owned()
                    ]
                );
                // A07 source-state semantics inside the composed result.
                assert_eq!(copied.source_state.name, "meta-src");
                assert_eq!(copied.source_state.resolution.to_bits(), f64::NAN.to_bits());
                assert!(copied.source_state.has_origx);
                assert_eq!(copied.source_state.raw_remarks, vec!["rem"]);
                assert!(copied.source_state.conect_map.is_empty());
                assert!(!copied.source_state.has_d_fraction);
                assert_eq!(copied.source_state.non_ascii_line, 0);
                assert_eq!(copied.source_state.ter_status, 0);
                copied.validate().unwrap();
            }
            // Original sharing/storage stability after all seven profiles.
            assert_eq!(source.models.len(), 2);
            assert_eq!(source.chains.len(), 2);
            assert_eq!(source.atoms.len(), 2);
            assert_eq!(source.source_state.conect_map.len(), 1);
            assert!(source.source_state.has_d_fraction);
            assert_eq!(source.source_state.ter_status, b'Z');
            assert_eq!(source.entities[0].subchains().len(), 3);
        }
        assert_eq!(calls, 28);
    }

    #[test]
    fn bio_copy_a09_source_metadata_reference_states_28() {
        use crate::hierarchy::{
            BioAssemblyGenerator, BioAssemblyOperator, BioAssemblySpecialKind, BioChainId,
            BioCoordinateFormat, BioCrystalCell, BioCrystalInfo, BioModelId, BioNcsOperator,
            BioStructureData, BioTransform,
        };
        use crate::metadata::BioMetadata;
        use crate::relationships::{BioCisPep, BioConnection, BioModRes};
        use crate::secondary_structure::{BioHelix, BioSheet};

        // Frozen hierarchy: auth A/B chains with label names S1/S2, models
        // 1/2, residues ALA/GLY, the two original atoms/coordinates. The
        // four reference states change ONLY the one assembly generator's
        // source lists; every other block is identical across states.
        let state_lists: [(&str, Vec<String>, Vec<String>); 4] = [
            ("auth-resolved", vec!["A".into(), "B".into()], vec![]),
            ("auth-unresolved", vec!["Z".into(), "Q".into()], vec![]),
            ("subchain-resolved", vec![], vec!["S1".into(), "S2".into()]),
            (
                "subchain-unresolved",
                vec![],
                vec!["Z1".into(), "Z2".into()],
            ),
        ];
        let identity = BioTransform::identity();
        let build_state = |chains: Vec<String>, subchains: Vec<String>| {
            let mut parts = a09_parts(Some(b"A"), Some("S1"), Some(b"B"), Some("S2"));
            // ONLY the generator's source lists change per state.
            parts.assemblies[0].generators[0] = BioAssemblyGenerator::new(
                chains,
                subchains,
                vec![BioAssemblyOperator::new(
                    Some("op".to_owned()),
                    Some("type".to_owned()),
                    identity,
                )],
            );
            // Nonempty retained categories (this fixture only).
            parts.connections = vec![BioConnection {
                name: "conn".to_owned(),
                reported_distance: -0.0,
                ..Default::default()
            }];
            parts.cispeps = vec![BioCisPep {
                model_num: 2,
                reported_angle: f64::NAN,
                ..Default::default()
            }];
            parts.mod_residues = vec![BioModRes {
                parent_comp_id: "ALA".to_owned(),
                mod_id: "mod".to_owned(),
                details: "detail".to_owned(),
                ..Default::default()
            }];
            parts.helices = vec![BioHelix {
                length: 7,
                ..Default::default()
            }];
            parts.sheets = vec![BioSheet::new("sheet")];
            let mut metadata = BioMetadata::default();
            metadata.authors = vec!["author".to_owned(), "author".to_owned()];
            metadata.solved_by = "solver".to_owned();
            parts.metadata = metadata;
            parts.crystal = Some(BioCrystalInfo::new(
                BioCrystalCell::default(),
                Some("P 1".to_owned()),
                Some("2".to_owned()),
                identity,
                identity,
                false,
                0,
                Vec::new(),
            ));
            parts.ncs_operators = vec![BioNcsOperator::new("ncs".to_owned(), true, identity)];
            parts
                .source_state
                .info
                .insert("key".to_owned(), "value".to_owned());
            parts.source_state.raw_remarks = vec!["rem".to_owned(), "rem".to_owned()];
            parts
        };
        let with = |modify: &dyn Fn(super::super::Selection) -> super::super::Selection| {
            crate::selection::BioSelectionData {
                selection: modify(super::super::Selection::default()),
            }
        };
        let sel = |field: &str, v: &str| {
            with(&|s| {
                let mut next = s;
                match field {
                    "chain" => {
                        next.chain_ids = super::super::SelectionList {
                            all: false,
                            inverted: false,
                            list: v.to_owned(),
                        }
                    }
                    "residue" => {
                        next.residue_names = super::super::SelectionList {
                            all: false,
                            inverted: false,
                            list: v.to_owned(),
                        }
                    }
                    "atom" => {
                        next.atom_names = super::super::SelectionList {
                            all: false,
                            inverted: false,
                            list: v.to_owned(),
                        }
                    }
                    _ => {}
                }
                next
            })
        };
        // Seven FIXED selectors; one-chain is always "A".
        let profiles: Vec<(
            &str,
            crate::selection::BioSelectionData,
            Vec<(u32, u32)>,
            Vec<(u32, u32)>,
            Vec<u32>,
            Vec<(u32, u32)>,
            Vec<u32>,
            Vec<u32>,
            Vec<usize>,
        )> = vec![
            (
                "all",
                with(&|s| s),
                vec![(0, 1), (1, 1)],
                vec![(0, 1), (1, 1)],
                vec![0, 1],
                vec![(0, 1), (1, 1)],
                vec![0, 1],
                vec![0, 1],
                vec![0, 1],
            ),
            (
                "no-model",
                with(&|s| super::super::Selection { mdl: 99, ..s }),
                vec![],
                vec![],
                vec![],
                vec![],
                vec![],
                vec![],
                vec![],
            ),
            (
                "no-chain",
                sel("chain", "NOPE"),
                vec![(0, 0), (0, 0)],
                vec![],
                vec![],
                vec![],
                vec![],
                vec![],
                vec![],
            ),
            (
                "no-residue",
                sel("residue", "XXX"),
                vec![(0, 1), (1, 1)],
                vec![(0, 0), (0, 0)],
                vec![0, 1],
                vec![],
                vec![],
                vec![],
                vec![],
            ),
            (
                "no-atom",
                sel("atom", "O"),
                vec![(0, 1), (1, 1)],
                vec![(0, 1), (1, 1)],
                vec![0, 1],
                vec![(0, 0), (0, 0)],
                vec![0, 1],
                vec![],
                vec![],
            ),
            (
                "one-chain",
                sel("chain", "A"),
                vec![(0, 1), (1, 0)],
                vec![(0, 1)],
                vec![0],
                vec![(0, 1)],
                vec![0],
                vec![0],
                vec![0],
            ),
            (
                "one-residue",
                sel("residue", "GLY"),
                vec![(0, 1), (1, 1)],
                vec![(0, 0), (0, 1)],
                vec![0, 1],
                vec![(0, 1)],
                vec![1],
                vec![0],
                vec![1],
            ),
        ];
        let mut calls = 0usize;
        for (state_label, chains, subchains) in state_lists {
            let parts = build_state(chains, subchains);
            assert!(!parts.connections.is_empty());
            assert!(!parts.cispeps.is_empty());
            assert!(!parts.mod_residues.is_empty());
            assert!(!parts.helices.is_empty());
            assert!(!parts.sheets.is_empty());
            assert_eq!(parts.metadata.authors.len(), 2);
            assert!(parts.crystal.is_some());
            assert_eq!(parts.ncs_operators.len(), 1);
            let source = BioStructureData::from_parts(parts).expect("state fixture");
            // Step 10 baselines for this state's seven calls — captured once,
            // never refreshed, checked after EVERY copier call below.
            let baseline_models_ptr = std::sync::Arc::as_ptr(&source.models);
            let baseline_atoms_ptr = std::sync::Arc::as_ptr(&source.atoms);
            let baseline_assemblies_ptr = std::sync::Arc::as_ptr(&source.assemblies);
            let baseline_source_state_ptr = std::sync::Arc::as_ptr(&source.source_state);
            let baseline_chains_ptr = std::sync::Arc::as_ptr(&source.chains);
            let baseline_residues_ptr = std::sync::Arc::as_ptr(&source.residues);
            let baseline_coordinates_ptr = std::sync::Arc::as_ptr(&source.coordinates);
            // Fixed two-row source coordinate-bit snapshot (never refreshed).
            assert_eq!(source.coordinates.positions().len(), 2);
            let baseline_position_bits: [[u64; 3]; 2] = [
                source.coordinates.positions()[0].map(f64::to_bits),
                source.coordinates.positions()[1].map(f64::to_bits),
            ];
            let baseline_atom0_bits = source.atoms[0].occupancy().to_bits();
            let baseline_resolution_bits = source.source_state.resolution.to_bits();
            for (
                profile,
                selection,
                model_spans,
                chain_spans,
                chain_model_ids,
                residue_spans,
                residue_chain_ids,
                atom_residue_ids,
                source_atom_indices,
            ) in &profiles
            {
                let copied = copy_bio_selection(&source, selection).unwrap();
                calls += 1;
                assert_eq!(
                    copied.models.len(),
                    model_spans.len(),
                    "{state_label}/{profile}"
                );
                for (row, (start, len)) in copied.models.iter().zip(model_spans) {
                    assert_eq!(
                        (row.chain_span().start(), row.chain_span().len()),
                        (*start, *len)
                    );
                }
                // Model numbers 1/2 retained even for empty accepted parents.
                for (index, row) in copied.models.iter().enumerate() {
                    assert_eq!(row.source_model_number(), Some(index as i32 + 1));
                }
                assert_eq!(copied.chains.len(), chain_spans.len());
                for ((row, (start, len)), model_id) in
                    copied.chains.iter().zip(chain_spans).zip(chain_model_ids)
                {
                    assert_eq!(
                        (row.residue_span().start(), row.residue_span().len()),
                        (*start, *len)
                    );
                    assert_eq!(row.model_id(), BioModelId::new(*model_id));
                }
                assert_eq!(copied.residues.len(), residue_spans.len());
                for ((row, (start, len)), chain_id) in copied
                    .residues
                    .iter()
                    .zip(residue_spans)
                    .zip(residue_chain_ids)
                {
                    assert_eq!(
                        (row.atom_span().start(), row.atom_span().len()),
                        (*start, *len)
                    );
                    assert_eq!(row.chain_id(), BioChainId::new(*chain_id));
                }
                assert_eq!(copied.atoms.len(), atom_residue_ids.len());
                for (row, residue_id) in copied.atoms.iter().zip(atom_residue_ids) {
                    assert_eq!(row.residue_id(), BioResidueId::new(*residue_id));
                }
                for (copied_position, original) in copied
                    .coordinates
                    .positions()
                    .iter()
                    .zip(source_atom_indices)
                {
                    for (a, b) in copied_position
                        .iter()
                        .zip(source.coordinates.positions()[*original])
                    {
                        assert_eq!(a.to_bits(), b.to_bits());
                    }
                }
                // All 10 unchanged metadata Arc identities in EVERY call.
                assert!(std::sync::Arc::ptr_eq(&copied.entities, &source.entities));
                assert!(std::sync::Arc::ptr_eq(
                    &copied.connections,
                    &source.connections
                ));
                assert!(std::sync::Arc::ptr_eq(&copied.cispeps, &source.cispeps));
                assert!(std::sync::Arc::ptr_eq(
                    &copied.mod_residues,
                    &source.mod_residues
                ));
                assert!(std::sync::Arc::ptr_eq(&copied.helices, &source.helices));
                assert!(std::sync::Arc::ptr_eq(&copied.sheets, &source.sheets));
                assert!(std::sync::Arc::ptr_eq(&copied.metadata, &source.metadata));
                assert!(std::sync::Arc::ptr_eq(&copied.crystal, &source.crystal));
                assert!(std::sync::Arc::ptr_eq(
                    &copied.ncs_operators,
                    &source.ncs_operators
                ));
                assert!(std::sync::Arc::ptr_eq(
                    &copied.assemblies,
                    &source.assemblies
                ));
                assert_eq!(copied.input_format, BioCoordinateFormat::Pdb);
                copied.validate().unwrap();

                // ---- Step 10 exact bit/field proof (new matrix only) ----
                // Retained generator source lists, operator strings and
                // transform components by to_bits.
                let generator = &copied.assemblies[0].generators[0];
                assert_eq!(generator.chains, source.assemblies[0].generators[0].chains);
                assert_eq!(
                    generator.subchains,
                    source.assemblies[0].generators[0].subchains
                );
                let operator = &generator.operators[0];
                assert_eq!(operator.name.as_deref(), Some("op"));
                assert_eq!(operator.operator_type.as_deref(), Some("type"));
                for (copied_row, source_row) in operator.transform.matrix().iter().zip(
                    source.assemblies[0].generators[0].operators[0]
                        .transform
                        .matrix(),
                ) {
                    for (a, b) in copied_row.iter().zip(source_row) {
                        assert_eq!(a.to_bits(), b.to_bits());
                    }
                }
                for (a, b) in operator.transform.translation().iter().zip(
                    source.assemblies[0].generators[0].operators[0]
                        .transform
                        .translation(),
                ) {
                    assert_eq!(a.to_bits(), b.to_bits());
                }
                let assembly = &copied.assemblies[0];
                assert_eq!(assembly.name, "asm");
                assert!(assembly.author_determined);
                assert!(!assembly.software_determined);
                assert_eq!(assembly.special_kind, BioAssemblySpecialKind::NotApplicable);
                assert_eq!(assembly.oligomeric_count, 2);
                assert_eq!(assembly.oligomeric_details, "od");
                assert_eq!(assembly.software_name, "sw");
                assert_eq!(
                    assembly.buried_surface_area.to_bits(),
                    source.assemblies[0].buried_surface_area.to_bits()
                );
                assert!(assembly.surface_area.is_nan());
                assert_eq!(
                    assembly.surface_area.to_bits(),
                    source.assemblies[0].surface_area.to_bits()
                );
                assert_eq!(
                    assembly.solvent_free_energy_change.to_bits(),
                    source.assemblies[0].solvent_free_energy_change.to_bits()
                );
                // Fixed named category values and float carriers.
                assert_eq!(copied.connections[0].name, "conn");
                assert_eq!(
                    copied.connections[0].reported_distance.to_bits(),
                    source.connections[0].reported_distance.to_bits()
                );
                assert_eq!(
                    copied.connections[0].reported_distance.to_bits(),
                    (-0.0_f64).to_bits()
                );
                assert_eq!(copied.cispeps[0].model_num, 2);
                assert_eq!(
                    copied.cispeps[0].reported_angle.to_bits(),
                    source.cispeps[0].reported_angle.to_bits()
                );
                assert_eq!(copied.mod_residues[0].parent_comp_id, "ALA");
                assert_eq!(copied.mod_residues[0].mod_id, "mod");
                assert_eq!(copied.mod_residues[0].details, "detail");
                assert_eq!(copied.helices[0].length, 7);
                assert_eq!(copied.sheets[0].name, "sheet");
                assert_eq!(copied.metadata.authors, vec!["author", "author"]);
                assert_eq!(copied.metadata.solved_by, "solver");
                let crystal = (*copied.crystal).as_ref().unwrap();
                assert_eq!(crystal.space_group_hm().as_deref(), Some("P 1"));
                assert_eq!(crystal.z_pdb().as_deref(), Some("2"));
                assert_eq!(copied.ncs_operators[0].id, "ncs");
                assert!(copied.ncs_operators[0].given);
                // source_state: retained fields, info map, ordered duplicate
                // remarks, origx component bits and the four resets.
                assert_eq!(copied.source_state.name, "meta-src");
                assert_eq!(
                    copied.source_state.resolution.to_bits(),
                    source.source_state.resolution.to_bits()
                );
                assert!(copied.source_state.has_origx);
                for (a, b) in copied
                    .source_state
                    .origx
                    .matrix()
                    .iter()
                    .zip(source.source_state.origx.matrix())
                {
                    for (x, y) in a.iter().zip(b) {
                        assert_eq!(x.to_bits(), y.to_bits());
                    }
                }
                for (a, b) in copied
                    .source_state
                    .origx
                    .translation()
                    .iter()
                    .zip(source.source_state.origx.translation())
                {
                    assert_eq!(a.to_bits(), b.to_bits());
                }
                assert_eq!(
                    copied.source_state.info.get("key").map(String::as_str),
                    Some("value")
                );
                assert_eq!(copied.source_state.raw_remarks, vec!["rem", "rem"]);
                assert!(copied.source_state.conect_map.is_empty());
                assert!(!copied.source_state.has_d_fraction);
                assert_eq!(copied.source_state.non_ascii_line, 0);
                assert_eq!(copied.source_state.ter_status, 0);
                // Every selected atom: full field bits + position components
                // against its frozen original source index, original
                // serial/name and the remapped parent (asserted above).
                for (atom, original) in copied.atoms.iter().zip(source_atom_indices) {
                    let source_atom = &source.atoms[*original];
                    assert_eq!(
                        atom.occupancy().to_bits(),
                        source_atom.occupancy().to_bits()
                    );
                    assert_eq!(atom.b_iso().to_bits(), source_atom.b_iso().to_bits());
                    assert_eq!(atom.fraction().to_bits(), source_atom.fraction().to_bits());
                    for (a, b) in atom.anisou().iter().zip(source_atom.anisou()) {
                        assert_eq!(a.to_bits(), b.to_bits());
                    }
                    assert_eq!(atom.source().serial(), source_atom.source().serial());
                    assert_eq!(atom.name().as_bytes(), source_atom.name().as_bytes());
                }
                // Per-call unchanged-source checkpoint: storage pointers and
                // scalar sentinels captured before this state's seven calls
                // (28 across four states) remain unchanged after every
                // copier call.
                assert_eq!(std::sync::Arc::as_ptr(&source.chains), baseline_chains_ptr);
                assert_eq!(
                    std::sync::Arc::as_ptr(&source.residues),
                    baseline_residues_ptr
                );
                assert_eq!(
                    std::sync::Arc::as_ptr(&source.coordinates),
                    baseline_coordinates_ptr
                );
                assert_eq!(source.coordinates.positions().len(), 2);
                for (row, baseline_row) in source
                    .coordinates
                    .positions()
                    .iter()
                    .zip(baseline_position_bits)
                {
                    for (component, baseline_component) in row.iter().zip(baseline_row) {
                        assert_eq!(component.to_bits(), baseline_component);
                    }
                }
                assert_eq!(std::sync::Arc::as_ptr(&source.models), baseline_models_ptr);
                assert_eq!(std::sync::Arc::as_ptr(&source.atoms), baseline_atoms_ptr);
                assert_eq!(
                    std::sync::Arc::as_ptr(&source.assemblies),
                    baseline_assemblies_ptr
                );
                assert_eq!(
                    std::sync::Arc::as_ptr(&source.source_state),
                    baseline_source_state_ptr
                );
                assert_eq!(source.atoms[0].occupancy().to_bits(), baseline_atom0_bits);
                assert_eq!(
                    source.source_state.resolution.to_bits(),
                    baseline_resolution_bits
                );
            }
        }
        assert_eq!(calls, 28);
    }
}
