//! Search behaviour for the query graph value owned by `cosmolkit-model`.

pub use cosmolkit_model::{
    AtomQueryPredicate, BondQueryPredicate, QueryAtom, QueryBond, QueryGraph, QueryGraphError,
    QueryNode,
};

use crate::Molecule;

pub use cosmolkit_model::{
    AtomRangeBounds, AtomRangeDataFunction, AtomRangeQuery, RecursiveStructureQuery,
};

/// Operations that interpret a query graph.
pub struct QueryGraphOperator<'a> {
    inner: &'a QueryGraph,
}

impl std::fmt::Debug for QueryGraphOperator<'_> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("QueryGraphOperator")
            .field("inner", &self.inner)
            .finish()
    }
}

impl<'a> QueryGraphOperator<'a> {
    #[must_use]
    pub const fn new(inner: &'a QueryGraph) -> Self {
        Self { inner }
    }

    #[must_use]
    pub const fn inner(&self) -> &'a QueryGraph {
        self.inner
    }

    pub fn to_smarts(
        &self,
        params: &crate::SmilesWriteParams,
    ) -> Result<String, crate::SmartsWriteError> {
        crate::search::smarts_write::query_graph_to_smarts(self.inner, params)
    }

    pub fn to_cx_smarts(
        &self,
        params: &crate::SmilesWriteParams,
    ) -> Result<String, crate::SmartsWriteError> {
        crate::search::smarts_write::query_graph_to_cx_smarts(self.inner, params)
    }

    pub fn fragment_to_smarts(
        &self,
        params: &crate::SmilesWriteParams,
        atom_ids: &[crate::AtomId],
        bond_ids: Option<&[crate::BondId]>,
    ) -> Result<String, crate::SmartsWriteError> {
        crate::search::smarts_write::query_graph_fragment_to_smarts(
            self.inner, params, atom_ids, bond_ids,
        )
    }

    pub fn fragment_to_cx_smarts(
        &self,
        params: &crate::SmilesWriteParams,
        atom_ids: &[crate::AtomId],
        bond_ids: Option<&[crate::BondId]>,
    ) -> Result<String, crate::SmartsWriteError> {
        crate::search::smarts_write::query_graph_fragment_to_cx_smarts(
            self.inner, params, atom_ids, bond_ids,
        )
    }

    /// Serialize one query atom while preserving its query predicate tree.
    pub fn atom_to_smarts(
        &self,
        atom_id: crate::AtomId,
        params: &crate::SmilesWriteParams,
    ) -> Result<String, crate::SmartsWriteError> {
        let atom = self.inner.atom(atom_id.index()).ok_or(
            crate::SmartsWriteError::FragmentAtomOutOfRange {
                atom: atom_id.index(),
            },
        )?;
        crate::search::smarts_write::query_atom_to_smarts(atom, params)
    }

    /// Serialize one query bond while preserving its query predicate tree.
    pub fn bond_to_smarts(
        &self,
        bond_id: crate::BondId,
        params: &crate::SmilesWriteParams,
        atom_to_left: Option<crate::AtomId>,
    ) -> Result<String, crate::SmartsWriteError> {
        let bond = self.inner.bond(bond_id.index()).ok_or(
            crate::SmartsWriteError::FragmentBondOutOfRange {
                bond: bond_id.index(),
            },
        )?;
        let atom_to_left = atom_to_left
            .map(|atom| {
                if atom != bond.begin() && atom != bond.end() {
                    Err(crate::SmartsWriteError::BondAtomNotEndpoint {
                        bond: bond_id.index(),
                        atom: atom.index(),
                    })
                } else {
                    Ok(atom.index())
                }
            })
            .transpose()?;
        crate::search::smarts_write::query_bond_to_smarts(bond, params, atom_to_left)
    }

    #[must_use]
    pub fn compile(&self) -> CompiledQuery {
        CompiledQuery::compile(self.inner.clone())
    }

    pub fn matches(&self, molecule: &Molecule) -> Vec<MatchResult> {
        crate::search::substruct::get_substruct_matches(molecule, self.inner)
    }
}

/// Build an exact-attribute query graph from a concrete molecule.
pub(crate) fn query_graph_from_concrete_molecule(
    molecule: &Molecule,
) -> Result<QueryGraph, QueryGraphError> {
    let atoms = molecule
        .atoms()
        .iter()
        .map(|atom| {
            let predicate =
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(atom.atomic_number()));
            QueryAtom::from_parts(atom.clone(), predicate)
        })
        .collect::<Vec<_>>();
    let bonds = molecule
        .bonds()
        .iter()
        .map(|bond| {
            let predicate = if bond.order() == crate::BondOrder::Unspecified {
                QueryNode::predicate(BondQueryPredicate::Any)
            } else {
                QueryNode::predicate(BondQueryPredicate::Order(bond.order()))
            };
            QueryBond::from_parts(bond.clone(), predicate)
        })
        .collect::<Vec<_>>();
    let mut props = molecule.properties().props().clone();
    if let Some(name) = molecule.properties().name() {
        props.insert("_Name".to_owned(), name.to_owned());
    }
    QueryGraph::from_parts(
        atoms,
        bonds,
        props,
        molecule.coordinates_2d().map_or_else(Vec::new, |coords| {
            vec![cosmolkit_model::Conformer2D::new(0, coords.to_vec())]
        }),
        molecule.conformers_3d().to_vec(),
        molecule.stereo_groups().to_vec(),
    )
}

/// Remove parser-only properties after SMARTS parsing.
pub(crate) fn cleanup_query_graph_parser_state(graph: &mut QueryGraph) {
    for atom in graph.atoms_mut() {
        atom.atom_mut().clear_prop("_RingClosures");
        atom.atom_mut()
            .clear_prop(crate::notation::smiles::SMILES_START_PROP);
    }
    for bond in graph.bonds_mut() {
        bond.bond_mut()
            .clear_prop(crate::notation::smiles::UNSPECIFIED_ORDER_PROP);
        bond.bond_mut()
            .clear_prop(crate::notation::smiles::CXSMILES_BOND_IDX_PROP);
    }
}

pub(crate) fn query_graph_has_directional_bonds(graph: &QueryGraph) -> bool {
    graph
        .bonds()
        .iter()
        .any(|bond| bond.bond().direction() != crate::BondDirection::None)
}

fn query_graph_neighboring_directed_bond(
    graph: &QueryGraph,
    atom: crate::AtomId,
) -> Option<crate::BondId> {
    // BEGIN RDKIT CPP FUNCTION getNeighboringDirectedBond
    // RDKit✔️✔️: const Bond *getNeighboringDirectedBond(const ROMol &mol, const Atom *atom) {
    // RDKit✔️✔️:   PRECONDITION(atom, "no atom");
    // RDKit✔️✔️:   for (const auto &bondIdx :
    // RDKit✔️✔️:        boost::make_iterator_range(mol.getAtomBonds(atom))) {
    // RDKit✔️✔️:     const Bond *bond = mol[bondIdx];
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (bond->getBondType() != Bond::BondType::DOUBLE &&
    // RDKit✔️✔️:         hasStereoBondDir(bond)) {
    // RDKit✔️✔️:       return bond;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return nullptr;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION getNeighboringDirectedBond
    for (_, bond_index) in graph.adjacency().get(atom.index())?.iter().copied() {
        let bond = graph.bonds().get(bond_index)?;
        if bond.bond().order() != crate::BondOrder::Double
            && matches!(
                bond.bond().direction(),
                crate::BondDirection::EndDownRight | crate::BondDirection::EndUpRight
            )
        {
            return Some(bond.id());
        }
    }
    None
}

fn query_graph_opposite_direction(direction: crate::BondDirection) -> crate::BondDirection {
    match direction {
        crate::BondDirection::EndUpRight => crate::BondDirection::EndDownRight,
        crate::BondDirection::EndDownRight => crate::BondDirection::EndUpRight,
        other => other,
    }
}

/// Apply RDKit's directional-bond double-bond stereo inference directly to a
/// detached query graph.
pub(crate) fn set_query_graph_bond_stereo_from_directions(
    graph: &mut QueryGraph,
) -> Result<(), crate::StereoError> {
    // BEGIN RDKIT CPP FUNCTION setBondStereoFromDirections
    // RDKit✔️✔️: void setBondStereoFromDirections(ROMol &mol) {
    // RDKit✔️✔️:   mol.clearProp("_needsDetectBondStereo");
    // RDKit✔️✔️:   for (Bond *bond : mol.bonds()) {
    // RDKit✔️✔️:     if (bond->getBondType() == Bond::DOUBLE &&
    // RDKit✔️✔️:         bond->getStereo() != Bond::STEREOANY) {
    // RDKit✔️✔️:       const Atom *stereoBondBeginAtom = bond->getBeginAtom();
    // RDKit✔️✔️:       const Atom *stereoBondEndAtom = bond->getEndAtom();
    // RDKit✔️✔️:       const Bond *directedBondAtBegin =
    // RDKit✔️✔️:           Chirality::getNeighboringDirectedBond(mol, stereoBondBeginAtom);
    // RDKit✔️✔️:       const Bond *directedBondAtEnd =
    // RDKit✔️✔️:           Chirality::getNeighboringDirectedBond(mol, stereoBondEndAtom);
    // RDKit✔️✔️:       if (directedBondAtBegin != nullptr && directedBondAtEnd != nullptr) {
    // RDKit✔️✔️:         unsigned beginSideStereoAtom =
    // RDKit✔️✔️:             directedBondAtBegin->getOtherAtomIdx(stereoBondBeginAtom->getIdx());
    // RDKit✔️✔️:         unsigned endSideStereoAtom =
    // RDKit✔️✔️:             directedBondAtEnd->getOtherAtomIdx(stereoBondEndAtom->getIdx());
    // RDKit✔️✔️:         bond->setStereoAtoms(beginSideStereoAtom, endSideStereoAtom);
    // RDKit✔️✔️:         auto beginSideBondDirection = directedBondAtBegin->getBondDir();
    // RDKit✔️✔️:         if (directedBondAtBegin->getBeginAtom() == stereoBondBeginAtom) {
    // RDKit✔️✔️:           beginSideBondDirection = getOppositeBondDir(beginSideBondDirection);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         auto endSideBondDirection = directedBondAtEnd->getBondDir();
    // RDKit✔️✔️:         if (directedBondAtEnd->getEndAtom() == stereoBondEndAtom) {
    // RDKit✔️✔️:           endSideBondDirection = getOppositeBondDir(endSideBondDirection);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         if (beginSideBondDirection == endSideBondDirection) {
    // RDKit✔️✔️:           bond->setStereo(Bond::STEREOTRANS);
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           bond->setStereo(Bond::STEREOCIS);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION setBondStereoFromDirections
    graph.clear_prop("_needsDetectBondStereo");
    let mut updates = Vec::new();
    for bond in graph.bonds() {
        if bond.bond().order() != crate::BondOrder::Double
            || bond.bond().stereo() == crate::BondStereo::Any
        {
            continue;
        }
        let begin = bond.begin();
        let end = bond.end();
        let Some(begin_dir_id) = query_graph_neighboring_directed_bond(graph, begin) else {
            continue;
        };
        let Some(end_dir_id) = query_graph_neighboring_directed_bond(graph, end) else {
            continue;
        };
        let Some(begin_dir) = graph.bonds().get(begin_dir_id.index()) else {
            continue;
        };
        let Some(end_dir) = graph.bonds().get(end_dir_id.index()) else {
            continue;
        };
        let begin_side_atom = if begin_dir.begin() == begin {
            begin_dir.end()
        } else {
            begin_dir.begin()
        };
        let end_side_atom = if end_dir.begin() == end {
            end_dir.end()
        } else {
            end_dir.begin()
        };
        let begin_direction = if begin_dir.begin() == begin {
            query_graph_opposite_direction(begin_dir.bond().direction())
        } else {
            begin_dir.bond().direction()
        };
        let end_direction = if end_dir.end() == end {
            query_graph_opposite_direction(end_dir.bond().direction())
        } else {
            end_dir.bond().direction()
        };
        let stereo = if begin_direction == end_direction {
            crate::BondStereo::Trans
        } else {
            crate::BondStereo::Cis
        };
        updates.push((bond.id(), [begin_side_atom, end_side_atom], stereo));
    }
    for (bond_id, stereo_atoms, stereo) in updates {
        let Some(bond) = graph.bonds_mut().get_mut(bond_id.index()) else {
            continue;
        };
        bond.bond_mut().set_stereo_atoms(Some(stereo_atoms));
        bond.bond_mut().set_stereo(stereo);
    }
    Ok(())
}

/// A compiled query ready for repeated matching.
#[derive(Debug, Clone, PartialEq)]
pub struct CompiledQuery {
    graph: QueryGraph,
    /// Pre-built topology view consumed by the VF2 planner.  Keeping this
    /// detached and immutable makes compilation useful across repeated target
    /// matches without exposing matcher internals in the model graph.
    vf2_graph: crate::search::substruct::CompiledQueryGraph,
    atom_order: Vec<usize>,
}

impl CompiledQuery {
    #[must_use]
    pub fn compile(graph: QueryGraph) -> Self {
        let vf2_graph = crate::search::substruct::compile_query_graph(&graph);
        let atom_order = crate::search::substruct::compile_query_order_from_graph(&vf2_graph);
        Self {
            graph,
            vf2_graph,
            atom_order,
        }
    }

    #[must_use]
    pub fn graph(&self) -> &QueryGraph {
        &self.graph
    }

    #[must_use]
    pub fn num_atoms(&self) -> usize {
        self.graph.num_atoms()
    }

    #[must_use]
    pub fn num_bonds(&self) -> usize {
        self.graph.num_bonds()
    }

    #[must_use]
    pub fn atom_order(&self) -> &[usize] {
        &self.atom_order
    }

    pub fn matches(&self, molecule: &Molecule) -> Vec<MatchResult> {
        crate::search::substruct::get_substruct_matches_with_compiled_query(
            molecule,
            &self.graph,
            &crate::search::substruct::SubstructMatchParams::default(),
            &self.vf2_graph,
            &self.atom_order,
        )
        .unwrap_or_default()
    }
}

/// A match mapping from query graph indices to target molecule indices.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MatchResult {
    pub atom_mapping: Vec<usize>,
    pub bond_mapping: Vec<usize>,
}

impl MatchResult {
    #[must_use]
    pub fn new(atom_mapping: Vec<usize>, bond_mapping: Vec<usize>) -> Self {
        Self {
            atom_mapping,
            bond_mapping,
        }
    }
}

/// Result metadata for a maximum-common-substructure search.
#[derive(Debug, Clone, PartialEq)]
pub struct McsResult {
    pub query: QueryGraph,
    pub atom_count: usize,
    pub bond_count: usize,
    pub completed: bool,
}

impl McsResult {
    #[must_use]
    pub fn new(query: QueryGraph, atom_count: usize, bond_count: usize, completed: bool) -> Self {
        Self {
            query,
            atom_count,
            bond_count,
            completed,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compiled_query_reuses_vf2_plan_without_changing_matches() {
        let query = crate::search::smarts_parse::parse_smarts(
            "CC",
            &crate::search::smarts_parse::SmartsParseParams::default(),
        )
        .expect("query parses");
        let molecule = Molecule::from_smiles("CCCO").expect("target parses");
        let compiled = CompiledQuery::compile(query.clone());

        assert_eq!(compiled.vf2_graph.num_atoms(), query.num_atoms());
        assert_eq!(compiled.vf2_graph.num_bonds(), query.num_bonds());
        assert_eq!(compiled.atom_order.len(), query.num_atoms());
        assert_eq!(
            crate::search::substruct::get_substruct_matches_with_params(
                &molecule,
                &query,
                &crate::search::substruct::SubstructMatchParams::default(),
            ),
            crate::search::substruct::get_substruct_matches_with_compiled_query(
                &molecule,
                &query,
                &crate::search::substruct::SubstructMatchParams::default(),
                &compiled.vf2_graph,
                &compiled.atom_order,
            )
            .expect("compiled query matches"),
        );
    }
}
