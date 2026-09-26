//! Explicit detached target data consumed by search algorithms.

use cosmolkit_core::{RingInfo, ValenceAssignment};
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomQueryPredicate, Bond, BondId, BondQueryPredicate, Conformer3D,
    CoordinateBlock, QueryNode, QueryStateError, QueryStateRef, StereoGroup, TopologyBlock,
};

/// Read-only target capability required by query evaluation.
///
/// This trait deliberately exposes model blocks and value assignments rather
/// than a live `Molecule`. Implementations cannot provide mutation or runtime
/// cache authority to the search crate.
pub trait SearchTargetAccess {
    fn topology_block(&self) -> &TopologyBlock;
    fn coordinate_block(&self) -> &CoordinateBlock;
    fn ring_info(&self) -> Option<&RingInfo>;
    fn valence(&self) -> Option<&ValenceAssignment>;

    #[doc(hidden)]
    fn atom_has_query(&self, _atom: AtomId) -> bool {
        false
    }

    #[doc(hidden)]
    fn atom_query_predicate(&self, _atom: AtomId) -> Option<&QueryNode<AtomQueryPredicate>> {
        None
    }

    #[doc(hidden)]
    fn bond_has_query(&self, _bond: BondId) -> bool {
        false
    }

    #[doc(hidden)]
    fn bond_query_predicate(&self, _bond: BondId) -> Option<&QueryNode<BondQueryPredicate>> {
        None
    }

    /// Atomic number visible to query evaluation. Ordinary targets use the
    /// validated model value; a detached depiction target may carry RDKit's
    /// temporary non-element sentinel without changing its topology.
    fn query_atomic_number(&self, atom: &Atom) -> u8 {
        atom.atomic_number()
    }

    fn stereo_groups(&self) -> &[StereoGroup] {
        &self.topology_block().stereo_groups
    }

    fn atoms(&self) -> &[Atom] {
        &self.topology_block().atoms
    }

    fn bonds(&self) -> &[Bond] {
        &self.topology_block().bonds
    }

    fn adjacency(&self) -> &AdjacencyList {
        &self.topology_block().adjacency
    }

    fn num_atoms(&self) -> usize {
        self.atoms().len()
    }

    fn num_bonds(&self) -> usize {
        self.bonds().len()
    }

    fn conformers_3d(&self) -> &[Conformer3D] {
        &self.coordinate_block().conformers_3d
    }
}

/// Borrowed detached search input assembled by the facade.
#[derive(Debug, Clone, Copy)]
pub struct SearchTarget<'a> {
    topology: &'a TopologyBlock,
    coordinates: &'a CoordinateBlock,
    stereo_groups: &'a [StereoGroup],
    ring_info: Option<&'a RingInfo>,
    valence: Option<&'a ValenceAssignment>,
    atomic_number_overrides: Option<&'a [Option<u8>]>,
    query_state: Option<QueryStateRef<'a>>,
}

impl<'a> SearchTarget<'a> {
    #[must_use]
    pub const fn new(
        topology: &'a TopologyBlock,
        coordinates: &'a CoordinateBlock,
        stereo_groups: &'a [StereoGroup],
        ring_info: Option<&'a RingInfo>,
        valence: Option<&'a ValenceAssignment>,
    ) -> Self {
        Self {
            topology,
            coordinates,
            stereo_groups,
            ring_info,
            valence,
            atomic_number_overrides: None,
            query_state: None,
        }
    }

    pub(crate) fn with_ring_info<'b>(&self, ring_info: &'b RingInfo) -> SearchTarget<'b>
    where
        'a: 'b,
    {
        SearchTarget {
            topology: self.topology,
            coordinates: self.coordinates,
            stereo_groups: self.stereo_groups,
            ring_info: Some(ring_info),
            valence: self.valence,
            atomic_number_overrides: self.atomic_number_overrides,
            query_state: self.query_state,
        }
    }

    /// Attach query predicates and source query origins aligned with this
    /// target's current topology. Carrier chemistry remains authoritative in
    /// `topology`; the validated overlay never exposes its stored carriers.
    #[doc(hidden)]
    pub fn try_with_query_state(
        mut self,
        query_state: QueryStateRef<'a>,
    ) -> Result<Self, QueryStateError> {
        query_state.validate_for_topology(self.topology)?;
        self.query_state = Some(query_state);
        Ok(self)
    }

    /// Attach temporary query-visible atomic numbers aligned with target atom
    /// indices. This does not modify canonical `Element` or the target block.
    pub fn with_atomic_number_overrides(mut self, overrides: &'a [Option<u8>]) -> Self {
        assert_eq!(overrides.len(), self.topology.atoms.len());
        self.atomic_number_overrides = Some(overrides);
        self
    }
}

impl SearchTargetAccess for SearchTarget<'_> {
    fn topology_block(&self) -> &TopologyBlock {
        self.topology
    }

    fn coordinate_block(&self) -> &CoordinateBlock {
        self.coordinates
    }

    fn ring_info(&self) -> Option<&RingInfo> {
        self.ring_info
    }

    fn valence(&self) -> Option<&ValenceAssignment> {
        self.valence
    }

    fn atom_has_query(&self, atom: AtomId) -> bool {
        self.query_state
            .is_some_and(|query_state| query_state.atom_has_query(atom))
    }

    fn atom_query_predicate(&self, atom: AtomId) -> Option<&QueryNode<AtomQueryPredicate>> {
        self.query_state
            .map(|query_state| query_state.atom_predicate(atom))
    }

    fn bond_has_query(&self, bond: BondId) -> bool {
        self.query_state
            .is_some_and(|query_state| query_state.bond_has_query(bond))
    }

    fn bond_query_predicate(&self, bond: BondId) -> Option<&QueryNode<BondQueryPredicate>> {
        self.query_state
            .map(|query_state| query_state.bond_predicate(bond))
    }

    fn query_atomic_number(&self, atom: &Atom) -> u8 {
        // RDKit❗✔️: constexpr int DUMMY_ATOMIC_NUM = 200;
        // RDKit❗✔️: for (auto &at : rs_mol.atoms()) {
        // RDKit❗✔️:   if (!rs_atoms.test(at->getIdx())) {
        // RDKit❗✔️:     at->setAtomicNum(DUMMY_ATOMIC_NUM);
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // Behavior: only the query-visible number changes; unlike the source
        // clone, the typed model atom remains a valid element. Other source
        // effects of setAtomicNum must be audited by the owning caller.
        // Complexity: one indexed optional read per query atom, with no clone
        // or allocation in this hot path, versus an O(V+E) source mol clone.
        self.atomic_number_overrides
            .and_then(|overrides| overrides[atom.id().index()])
            .unwrap_or_else(|| atom.atomic_number())
    }

    fn stereo_groups(&self) -> &[StereoGroup] {
        self.stereo_groups
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, BondSpec, QueryAtom, QueryAtomIdentity, QueryBond};
    use cosmolkit_types::{BondOrder, Element};

    fn topology(elements: &[Element], bond: Option<(usize, usize)>) -> TopologyBlock {
        let atoms = elements
            .iter()
            .enumerate()
            .map(|(index, element)| Atom::from_spec(AtomId::new(index), AtomSpec::new(*element)))
            .collect();
        let bonds = bond
            .map(|(begin, end)| {
                vec![Bond::from_spec(
                    BondId::new(0),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                )]
            })
            .unwrap_or_default();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
    }

    fn target<'a>(
        topology: &'a TopologyBlock,
        coordinates: &'a CoordinateBlock,
    ) -> SearchTarget<'a> {
        SearchTarget::new(topology, coordinates, &topology.stereo_groups, None, None)
    }

    #[test]
    fn q86_target_state_defaults_and_origins_are_independent_of_tree_shape() {
        let topology = topology(&[Element::C, Element::C], Some((0, 1)));
        let coordinates = CoordinateBlock::default();
        let plain = target(&topology, &coordinates);
        assert!(!plain.atom_has_query(AtomId::new(0)));
        assert!(!plain.bond_has_query(BondId::new(0)));
        assert_eq!(plain.atom_query_predicate(AtomId::new(0)), None);
        assert_eq!(plain.bond_query_predicate(BondId::new(0)), None);

        let same_atom_tree = QueryNode::predicate(AtomQueryPredicate::Any);
        let same_bond_tree = QueryNode::predicate(BondQueryPredicate::Any);
        let atoms = vec![
            QueryAtom::from_parts(topology.atoms[0].clone(), same_atom_tree.clone()),
            QueryAtom::from_carrier_parts(topology.atoms[1].clone(), same_atom_tree.clone()),
        ];
        let bonds = vec![QueryBond::from_parts(
            topology.bonds[0].clone(),
            same_bond_tree.clone(),
        )];
        let state = QueryStateRef::try_for_topology(&atoms, &bonds, &topology).unwrap();
        let attached = target(&topology, &coordinates)
            .try_with_query_state(state)
            .unwrap();

        assert!(attached.atom_has_query(AtomId::new(0)));
        assert!(!attached.atom_has_query(AtomId::new(1)));
        assert!(attached.bond_has_query(BondId::new(0)));
        assert_eq!(
            attached.atom_query_predicate(AtomId::new(0)),
            Some(&same_atom_tree)
        );
        assert_eq!(
            attached.bond_query_predicate(BondId::new(0)),
            Some(&same_bond_tree)
        );
    }

    #[test]
    fn q86_target_state_attachment_revalidates_counts_ids_and_endpoints() {
        let source = topology(&[Element::C, Element::O], Some((0, 1)));
        let coordinates = CoordinateBlock::default();
        let atoms: Vec<_> = source
            .atoms
            .iter()
            .cloned()
            .map(|atom| QueryAtom::from_parts(atom, QueryNode::predicate(AtomQueryPredicate::Any)))
            .collect();
        let bonds = vec![QueryBond::from_parts(
            source.bonds[0].clone(),
            QueryNode::predicate(BondQueryPredicate::Any),
        )];
        let atoms_before = atoms.clone();
        let bonds_before = bonds.clone();
        let source_before = source.clone();
        let state = QueryStateRef::try_for_topology(&atoms, &bonds, &source).unwrap();

        let short = topology(&[Element::C], None);
        assert_eq!(
            target(&short, &coordinates)
                .try_with_query_state(state)
                .unwrap_err(),
            QueryStateError::AtomCount {
                actual: 2,
                expected: 1,
            }
        );

        let mut wrong_atom_id = source.clone();
        wrong_atom_id.atoms[0] = Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C));
        assert_eq!(
            target(&wrong_atom_id, &coordinates)
                .try_with_query_state(state)
                .unwrap_err(),
            QueryStateError::AtomId {
                position: 0,
                actual: AtomId::new(0),
                expected: AtomId::new(1),
            }
        );

        let mut wrong_bond_id = source.clone();
        wrong_bond_id.bonds[0] = Bond::from_spec(
            BondId::new(1),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        );
        assert_eq!(
            target(&wrong_bond_id, &coordinates)
                .try_with_query_state(state)
                .unwrap_err(),
            QueryStateError::BondId {
                position: 0,
                actual: BondId::new(0),
                expected: BondId::new(1),
            }
        );

        let reversed = topology(&[Element::C, Element::O], Some((1, 0)));
        assert_eq!(
            target(&reversed, &coordinates)
                .try_with_query_state(state)
                .unwrap_err(),
            QueryStateError::BondEndpoints {
                position: 0,
                actual: (AtomId::new(0), AtomId::new(1)),
                expected: (AtomId::new(1), AtomId::new(0)),
            }
        );
        assert_eq!(atoms, atoms_before);
        assert_eq!(bonds, bonds_before);
        assert_eq!(source, source_before);
    }

    #[test]
    fn q86_target_state_keeps_current_carriers_and_non_element_boundary() {
        let topology = topology(&[Element::C], None);
        let coordinates = CoordinateBlock::default();
        let old_carrier = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::O).with_formal_charge(-1),
        );
        let atoms = vec![QueryAtom::from_parts(
            old_carrier,
            QueryNode::predicate(AtomQueryPredicate::Any),
        )];
        let state = QueryStateRef::try_for_topology(&atoms, &[], &topology).unwrap();
        let attached = target(&topology, &coordinates)
            .try_with_query_state(state)
            .unwrap();
        assert_eq!(attached.query_atomic_number(&attached.atoms()[0]), 6);
        assert_eq!(attached.atoms()[0].formal_charge(), 0);

        let non_element = vec![QueryAtom::from_identity_parts(
            AtomId::new(0),
            QueryAtomIdentity::AtomicNumber(200),
            QueryNode::predicate(AtomQueryPredicate::Any),
        )];
        assert_eq!(
            QueryStateRef::try_for_topology(&non_element, &[], &topology).unwrap_err(),
            QueryStateError::NonElementAtomIdentity {
                position: 0,
                atom: AtomId::new(0),
                atomic_number: 200,
            }
        );
    }
}
