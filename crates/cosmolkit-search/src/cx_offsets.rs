//! Reaction-wide source-index windows for canonical query CX lowering.

use crate::{CxQueryLoweringError, apply_cx_to_query_graph};
use cosmolkit_cx::{CxRecord, ParsedCxExtensions};
use cosmolkit_model::QueryGraph;

/// Lower a complete reaction-global CX block onto one template. The existing
/// zero-offset lowerer remains the sole implementation of graph effects.
#[doc(hidden)]
pub fn apply_cx_to_query_graph_with_offsets(
    graph: &mut QueryGraph,
    parsed: &ParsedCxExtensions,
    start_atom: usize,
    start_bond: usize,
) -> Result<(), CxQueryLoweringError> {
    // RDKit❗❌: #define VALID_ATIDX(_atidx_) \
    // RDKit❗❌:   ((_atidx_) >= startAtomIdx && (_atidx_) < startAtomIdx + mol.getNumAtoms())
    // RDKit❗❌: #define VALID_BNDIDX(_bidx_) \
    // RDKit❗❌:   ((_bidx_) >= startBondIdx && (_bidx_) < startBondIdx + mol.getNumBonds())
    // RDKit❗❌: if (VALID_ATIDX(aidx) && VALID_BNDIDX(bidx)) {
    // RDKit❗❌:   auto bnd = get_bond_with_smiles_idx(mol, bidx - startBondIdx);
    // Complexity review: this explicit typed projection adds an O(records +
    // items) buffer to the already-decoded grammar records. There is no raw
    // text rewriting, sentinel index or clone of live molecule state.
    if start_atom == 0 && start_bond == 0 {
        return apply_cx_to_query_graph(graph, parsed);
    }
    let atom_count = graph.num_atoms();
    let bond_count = graph.num_bonds();
    let local_atom = |index: usize| index.checked_sub(start_atom).filter(|&i| i < atom_count);
    let local_bond = |index: usize| index.checked_sub(start_bond).filter(|&i| i < bond_count);
    let source_sub = |index: usize| -> Result<usize, CxQueryLoweringError> {
        let index = u32::try_from(index).map_err(|_| {
            CxQueryLoweringError::InvalidGraph("CX index exceeds source unsigned range".into())
        })?;
        let start = u32::try_from(start_atom).map_err(|_| {
            CxQueryLoweringError::InvalidGraph(
                "CX startAtomIdx exceeds source unsigned range".into(),
            )
        })?;
        Ok(index.wrapping_sub(start) as usize)
    };
    let mut records = Vec::with_capacity(parsed.records().len());
    for record in parsed.records() {
        let mut record = record.clone();
        match &mut record {
            CxRecord::Coordinates(values) => {
                values.values = values
                    .values
                    .iter()
                    .skip(start_atom)
                    .take(atom_count)
                    .copied()
                    .collect();
            }
            CxRecord::AtomLabels(values) => {
                *values = values
                    .iter()
                    .skip(start_atom)
                    .take(atom_count)
                    .cloned()
                    .collect();
            }
            CxRecord::AtomValues(values) => {
                // RDKit❗✔️: if (tkn != "" && VALID_ATIDX(atIdx)) {
                // RDKit❗✔️:   mol.getAtomWithIdx(atIdx)->setProp(RDKit::common_properties::molFileValue,
                // RDKit❗✔️:                                      tkn);
                // Unlike labels, this pinned helper uses the GLOBAL row after
                // window validation. Keep its in-range writes and structural
                // failure instead of correcting the source subtraction bug.
                for (index, value) in values.iter_mut().enumerate() {
                    if local_atom(index).is_none() {
                        *value = None;
                    } else if index >= atom_count && value.as_ref().is_some_and(|v| !v.is_empty()) {
                        return Err(CxQueryLoweringError::AtomIndex { index });
                    }
                }
            }
            CxRecord::AtomProperties(values) => values.retain_mut(|value| {
                if let Some(index) = local_atom(value.atom) {
                    value.atom = index;
                    true
                } else {
                    false
                }
            }),
            CxRecord::CoordinateBonds(values) => values.bonds.retain_mut(|value| {
                if let (Some(atom), Some(bond)) = (local_atom(value.atom), local_bond(value.bond)) {
                    value.atom = atom;
                    value.bond = bond;
                    true
                } else {
                    false
                }
            }),
            CxRecord::ZeroBonds(values)
            | CxRecord::DoubleBondStereo(cosmolkit_cx::CxDoubleBondStereo {
                bonds: values, ..
            }) => {
                *values = values
                    .iter()
                    .filter_map(|&value| local_bond(value))
                    .collect();
            }
            CxRecord::Unsaturation(values) => {
                *values = values
                    .iter()
                    .filter_map(|&value| local_atom(value))
                    .collect();
            }
            CxRecord::EnhancedStereo(values) => {
                values.atoms = values
                    .atoms
                    .iter()
                    .filter_map(|&value| local_atom(value))
                    .collect();
            }
            CxRecord::RingBonds(values) => values.retain_mut(|value| {
                if let Some(atom) = local_atom(value.atom) {
                    value.atom = atom;
                    true
                } else {
                    false
                }
            }),
            CxRecord::Substitution(values) => values.retain_mut(|value| {
                if let Some(atom) = local_atom(value.atom) {
                    value.atom = atom;
                    true
                } else {
                    false
                }
            }),
            CxRecord::WedgedBonds(values) => values.retain_mut(|value| {
                if let (Some(atom), Some(bond)) = (local_atom(value.atom), local_bond(value.bond)) {
                    value.atom = atom;
                    value.bond = bond;
                    true
                } else {
                    false
                }
            }),
            CxRecord::Radicals(values) => values.retain_mut(|value| {
                if let Some(atom) = local_atom(value.atom) {
                    value.atom = atom;
                    true
                } else {
                    false
                }
            }),
            CxRecord::LinkNodes(values) => {
                let mut nodes = Vec::with_capacity(values.len());
                for mut node in std::mem::take(values) {
                    let Some(atom) = local_atom(node.atom) else {
                        continue;
                    };
                    node.atom = atom;
                    // RDKit❗✔️: idx1 = *nbrs.first;
                    // RDKit❗✔️: nbrs.first++;
                    // RDKit❗✔️: idx2 = *nbrs.first;
                    // RDKit❗✔️: (idx1 - startAtomIdx + 1)
                    // Implicit neighbors are already local in the source and
                    // nevertheless undergo unsigned subtraction in its output.
                    let outer = if let Some(outer) = node.outer_atoms {
                        outer
                    } else {
                        let neighbors = &graph.adjacency()[atom];
                        if neighbors.len() != 2 {
                            return Err(CxQueryLoweringError::InvalidGraph(format!(
                                "CX link-node atom {atom} does not have degree two"
                            )));
                        }
                        [neighbors[0].0, neighbors[1].0]
                    };
                    node.outer_atoms = Some([source_sub(outer[0])?, source_sub(outer[1])?]);
                    nodes.push(node);
                }
                *values = nodes;
            }
            CxRecord::DataSGroup(values) => {
                values.atoms = values
                    .atoms
                    .iter()
                    .filter_map(|&value| local_atom(value))
                    .collect();
            }
            CxRecord::PolymerSGroup(values) => {
                values.atoms = values
                    .atoms
                    .iter()
                    .filter_map(|&value| local_atom(value))
                    .collect();
                // RDKit❗✔️: if (VALID_ATIDX(cidx)) {
                // RDKit❗✔️:   cidx -= startAtomIdx;
                // RDKit❗✔️: } else {
                // RDKit❗✔️:   keepSGroup = false;
                // These source crossing fields deliberately use ATOM bounds,
                // even though downstream finalization interprets bond rows.
                if values
                    .head_crossings
                    .iter()
                    .chain(&values.tail_crossings)
                    .any(|&v| local_atom(v).is_none())
                {
                    values.atoms.clear();
                } else {
                    values.head_crossings = values
                        .head_crossings
                        .iter()
                        .filter_map(|&v| local_atom(v))
                        .collect();
                    values.tail_crossings = values
                        .tail_crossings
                        .iter()
                        .filter_map(|&v| local_atom(v))
                        .collect();
                }
            }
            CxRecord::VariableAttachments(values) => values.retain_mut(|value| {
                if let Some(atom) = local_atom(value.atom) {
                    value.atom = atom;
                    value.endpoints = value
                        .endpoints
                        .iter()
                        .filter_map(|&v| local_atom(v))
                        .collect();
                    true
                } else {
                    false
                }
            }),
            // SGroup sequence IDs are local to each independent parse_it call;
            // dropping empty SGroup records would corrupt these ordinals.
            CxRecord::SGroupHierarchy(_) | CxRecord::Unknown(_) => {}
        }
        records.push(record);
    }
    apply_cx_to_query_graph(graph, &ParsedCxExtensions::new(records, parsed.consumed()))
}
