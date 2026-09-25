//! Detached attachment-point chemistry shared with the Molfile finalizer.
//! No live molecule or runtime commit authority is available in this owner.

use cosmolkit_model::{
    Atom, AtomId, AtomQueryPredicate, AtomSpec, Bond, BondId, BondOrder, BondQueryPredicate,
    BondSpec, CoordinateBlock, CoordinateValidationError, Element, QueryAtom, QueryBond, QueryNode,
    QueryStateError, QueryStateRef, TopologyBlock, TopologyMapping, TopologyValidationError,
    remap_query_rows_with_appended,
};

use crate::{
    hydrogens::{HydrogenError, grow_coordinate_rows, place_terminal_attachment_coordinates},
    valence::{ValenceError, assign_valence_state_for_atom_from_parts},
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AttachmentWarning {
    pub atom: AtomId,
    pub value: i32,
}

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum AttachmentExpansionError {
    #[error("invalid source topology: {0}")]
    Topology(#[from] TopologyValidationError),
    #[error("invalid source coordinates: {0}")]
    Coordinates(#[from] CoordinateValidationError),
    #[error("attachment value rows have length {actual}, expected {expected}")]
    ValueCount { actual: usize, expected: usize },
    #[error("attachment query transport failed: {0}")]
    Query(#[from] QueryStateError),
    #[error("attachment atom property failed: {0}")]
    Property(#[from] cosmolkit_model::AtomPropertyError),
    #[error("attachment topology edit failed: {0}")]
    Edit(#[from] cosmolkit_model::TopologyEditError),
    #[error("attachment property-cache calculation failed: {0}")]
    Valence(#[from] ValenceError),
    #[error("attachment coordinate placement failed: {0}")]
    Placement(#[from] HydrogenError),
}

#[derive(Debug, Clone)]
pub struct AttachmentExpansionResult {
    pub topology: TopologyBlock,
    pub coordinates: CoordinateBlock,
    pub mapping: TopologyMapping,
    pub query_rows: Option<(Vec<QueryAtom>, Vec<QueryBond>)>,
    pub appended_query_atoms: Vec<QueryAtom>,
    pub appended_query_bonds: Vec<QueryBond>,
    pub appended_valence: Vec<(AtomId, i32, i32)>,
    pub warnings: Vec<AttachmentWarning>,
}

/// Expand the parser's already converted integer attachment values in one
/// detached transaction. The IO owner alone interprets text fields; this
/// source-algorithm owner receives typed values aligned to the original rows.
pub fn expand_attachment_points(
    mut topology: TopologyBlock,
    mut coordinates: CoordinateBlock,
    query_state: Option<QueryStateRef<'_>>,
    attachment_values: &[Option<i32>],
    add_as_queries: bool,
    add_coords: bool,
) -> Result<AttachmentExpansionResult, AttachmentExpansionError> {
    // BEGIN RDKIT CPP FUNCTION details::addExplicitAttachmentPoint (append owner)
    // RDKit❗❌: unsigned int addExplicitAttachmentPoint(RWMol &mol, unsigned int atomIdx,
    // RDKit❗❌:                                         unsigned int val, bool addAsQuery,
    // RDKit❗❌:                                         bool addCoords) {
    // The query-construction branch of this function is anchored in
    // attachment_query_rows, its actual typed-predicate implementing helper.
    // Behavior review: this detached append follows the source's property,
    // bond, non-strict valence and optional-coordinate order. Degree-one
    // mixed-conformer direction isolation is intentional CK-COORD-001.
    // Complexity review: rebuilding the validated detached graph per append
    // is materially costlier than mutating RDKit's RWMol in place.
    // BEGIN RDKIT CPP FUNCTION MolOps::expandAttachmentPoints
    // RDKit❗❌: void expandAttachmentPoints(RWMol &mol, bool addAsQueries, bool addCoords) {
    // RDKit❗❌:   for (auto atom : mol.atoms()) {
    // RDKit❗❌:     int value;
    // RDKit❗❌:     if (atom->getPropIfPresent(common_properties::molAttachPoint, value)) {
    // RDKit❗❌:       std::vector<int> tgtVals;
    // RDKit❗❌:       if (value == 1 || value == -1) {
    // RDKit❗❌:         tgtVals.push_back(1);
    // RDKit❗❌:       }
    // RDKit❗❌:       if (value == 2 || value == -1) {
    // RDKit❗❌:         tgtVals.push_back(2);
    // RDKit❗❌:       }
    // RDKit❗❌:       if (tgtVals.empty()) {
    // RDKit❗❌:         BOOST_LOG(rdWarningLog)
    // RDKit❗❌:             << "Invalid value for molAttachPoint: " << value << " on atom "
    // RDKit❗❌:             << atom->getIdx() << ". Not expanding this atttachment point."
    // RDKit❗❌:             << std::endl;
    // RDKit❗❌:         continue;
    // RDKit❗❌:       }
    // RDKit❗❌:       for (auto tval : tgtVals) {
    // RDKit❗❌:         atom->clearProp(common_properties::molAttachPoint);
    // RDKit❗❌:         details::addExplicitAttachmentPoint(mol, atom->getIdx(), tval,
    // RDKit❗❌:                                             addAsQueries, addCoords);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION MolOps::expandAttachmentPoints
    // Behavior review: original-row iteration, 1-then-2 dispatch, property
    // clearing and typed warnings follow the source. Coordinate case 1 has
    // the explicitly approved CK-COORD-001 per-conformer isolation.
    // Complexity review: detached topology reconstruction after each append
    // clones/rebuilds graph state, materially costlier than RWMol mutation;
    // no claim of equivalent performance is made.
    topology.validate()?;
    coordinates.validate_for_atom_count(topology.atoms.len())?;
    let old_atoms = topology.atoms.len();
    let old_bonds = topology.bonds.len();
    if attachment_values.len() != old_atoms {
        return Err(AttachmentExpansionError::ValueCount {
            actual: attachment_values.len(),
            expected: old_atoms,
        });
    }
    if let Some(state) = query_state {
        state.validate_for_topology(&topology)?;
    }
    let mut appended_query_atoms = Vec::new();
    let mut appended_query_bonds = Vec::new();
    let mut appended_valence = Vec::new();
    let mut warnings = Vec::new();
    for (parent_index, value) in attachment_values.iter().copied().enumerate() {
        let Some(value) = value else { continue };
        let labels: &[i32] = match value {
            1 => &[1],
            2 => &[2],
            -1 => &[1, 2],
            _ => {
                warnings.push(AttachmentWarning {
                    atom: AtomId::new(parent_index),
                    value,
                });
                continue;
            }
        };
        let parent = AtomId::new(parent_index);
        for &label in labels {
            topology.atoms[parent_index].clear_prop("molAttachPoint");
            let atom_id = AtomId::new(topology.atoms.len());
            let bond_id = BondId::new(topology.bonds.len());
            let mut atom = Atom::from_spec(atom_id, AtomSpec::new(Element::DUMMY));
            // RDKit❗❌:   newAtom->setProp(common_properties::_fromAttachPoint, val);
            atom.set_prop("_fromAttchpt", label.to_string())?;
            let bond = Bond::from_spec(bond_id, BondSpec::new(parent, atom_id, BondOrder::Single));
            // RDKit❗❌:   bool updateLabel = false;
            // RDKit❗❌:   bool takeOwnership = true;
            // RDKit❗❌:   auto idx = mol.addAtom(newAtom, updateLabel, takeOwnership);
            // RDKit❗❌:   mol.addBond(atomIdx, idx, Bond::SINGLE);
            topology.atoms.push(atom);
            topology.bonds.push(bond);
            topology = TopologyBlock::try_from_parts(
                topology.atoms,
                topology.bonds,
                topology.substance_groups,
                topology.stereo_groups,
            )?;
            coordinates = grow_coordinate_rows(coordinates, 1);
            // The source updates only the new atom's cache, non-strictly,
            // immediately after its bond is added. The detached model has no
            // mutable cache; retain the calculated facts for the IO finalizer.
            // RDKit❗❌:   mol.getAtomWithIdx(idx)->updatePropertyCache(false);
            let (explicit, implicit) = assign_valence_state_for_atom_from_parts(
                &topology.atoms,
                &topology.bonds,
                &topology.adjacency,
                atom_id,
                false,
            )?;
            appended_valence.push((atom_id, explicit, implicit));
            // RDKit❗❌:   if (addCoords) {
            // RDKit❗❌:     setTerminalAtomCoords(mol, idx, atomIdx);
            // RDKit❗❌:   }
            if add_coords {
                coordinates = place_terminal_attachment_coordinates(
                    &topology,
                    coordinates,
                    atom_id,
                    parent,
                    bond_id,
                )?;
            }
            let (query_atom, query_bond) = attachment_query_rows(
                &topology.atoms[atom_id.index()],
                &topology.bonds[bond_id.index()],
                add_as_queries,
            );
            appended_query_atoms.push(query_atom);
            appended_query_bonds.push(query_bond);
            // RDKit❗❌:   return idx;
            // RDKit❗❌: }
            // END RDKIT CPP FUNCTION details::addExplicitAttachmentPoint (append owner)
        }
    }
    let mapping = TopologyMapping::with_appended(
        old_atoms,
        old_bonds,
        appended_query_atoms.len(),
        appended_query_bonds.len(),
    );
    let query_rows = query_state
        .map(|state| {
            remap_query_rows_with_appended(
                state,
                &topology,
                &mapping,
                &appended_query_atoms,
                &appended_query_bonds,
            )
        })
        .transpose()?;
    coordinates.validate_for_atom_count(topology.atoms.len())?;
    Ok(AttachmentExpansionResult {
        topology,
        coordinates,
        mapping,
        query_rows,
        appended_query_atoms,
        appended_query_bonds,
        appended_valence,
        warnings,
    })
}

/// Construct the typed query overlay for one source-appended attachment row.
/// The ordinary branch remains carrier-derived when the surrounding Molfile
/// record already needs uniform QueryGraph storage; the query branch is an
/// explicit null-query atom. In both branches the new single bond is ordinary.
pub fn attachment_query_rows(
    atom: &Atom,
    bond: &Bond,
    add_as_query: bool,
) -> (QueryAtom, QueryBond) {
    // BEGIN RDKIT CPP FUNCTION details::addExplicitAttachmentPoint (query branch)
    // RDKit❗❌: unsigned int addExplicitAttachmentPoint(RWMol &mol, unsigned int atomIdx,
    // RDKit❗❌:                                         unsigned int val, bool addAsQuery,
    // RDKit❗❌:                                         bool addCoords) {
    // RDKit❗❌:   Atom *newAtom = nullptr;
    // RDKit❗❌:   if (addAsQuery) {
    // RDKit❗❌:     newAtom = new QueryAtom(0);
    // RDKit❗❌:     newAtom->setQuery(RDKit::makeAtomNullQuery());
    // RDKit❗❌:   } else {
    // RDKit❗❌:     newAtom = new Atom(0);
    // RDKit❗❌:   }
    // END RDKIT CPP FUNCTION details::addExplicitAttachmentPoint (query branch)
    // Behavior review: source query identity belongs to the new atom only;
    // the new bond is an ordinary single bond even in a QueryGraph record.
    // Existing query rows are transported separately, never rebuilt by shape.
    // Complexity review: two bounded leaf nodes and two row clones per append;
    // no graph scans or predicate-tree clones occur in this constructor.
    let query_atom = if add_as_query {
        QueryAtom::from_parts(atom.clone(), QueryNode::predicate(AtomQueryPredicate::Any))
    } else {
        QueryAtom::from_carrier_parts(
            atom.clone(),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(atom.atomic_number())),
        )
    };
    let query_bond = QueryBond::from_carrier_parts(
        bond.clone(),
        QueryNode::predicate(BondQueryPredicate::Order(bond.order())),
    );
    (query_atom, query_bond)
}
