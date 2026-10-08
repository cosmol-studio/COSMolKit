//! Shared graph-carrier finalization; no syntax parsing or runtime authority.
use cosmolkit_model::{Atom, AtomId, Bond, BondId, PropertyValue, QueryAtom, QueryBond};
use cosmolkit_types::{BondOrder, ChiralTag};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum ParserCarrierError {
    #[error("parser property string conversion failed: {0}")]
    PropertyString(#[from] crate::PropertyStringError),
    #[error("atom property operation failed: {0}")]
    AtomProperty(#[from] cosmolkit_model::AtomPropertyError),
    #[error("invalid detached model: {0}")]
    Model(String),
}

mod sealed {
    pub trait AtomCarrier {}
    impl AtomCarrier for cosmolkit_model::Atom {}
    impl AtomCarrier for cosmolkit_model::QueryAtom {}
    pub trait BondCarrier {}
    impl BondCarrier for cosmolkit_model::Bond {}
    impl BondCarrier for cosmolkit_model::QueryBond {}
}

pub trait ParserAtom: sealed::AtomCarrier {
    fn chiral_tag(&self) -> ChiralTag;
    fn explicit_hydrogens(&self) -> u8;
    fn chiral_permutation(&self) -> Option<u32>;
    fn atomic_number(&self) -> u8;
    fn prop(&self, key: &str) -> Option<&PropertyValue>;
    fn clear_prop(&mut self, key: &str) -> Result<(), cosmolkit_model::AtomPropertyError>;
    fn set_attachment_point(
        &mut self,
        value: i32,
    ) -> Result<(), cosmolkit_model::AtomPropertyError>;
    fn set_chiral_tag(&mut self, value: ChiralTag);
    fn set_chiral_permutation(&mut self, value: Option<u32>);
}
macro_rules! parser_atom_access {
    ($ty:ty) => {
        impl ParserAtom for $ty {
            fn chiral_tag(&self) -> ChiralTag {
                <$ty>::chiral_tag(self)
            }
            fn explicit_hydrogens(&self) -> u8 {
                <$ty>::explicit_hydrogens(self)
            }
            fn chiral_permutation(&self) -> Option<u32> {
                <$ty>::chiral_permutation(self)
            }
            fn atomic_number(&self) -> u8 {
                <$ty>::atomic_number(self)
            }
            fn prop(&self, key: &str) -> Option<&PropertyValue> {
                <$ty>::prop(self, key)
            }
            fn clear_prop(&mut self, key: &str) -> Result<(), cosmolkit_model::AtomPropertyError> {
                <$ty>::clear_prop(self, key)
            }
            fn set_attachment_point(
                &mut self,
                value: i32,
            ) -> Result<(), cosmolkit_model::AtomPropertyError> {
                <$ty>::set_prop(self, "_fromAttchpt", value)
            }
            fn set_chiral_tag(&mut self, value: ChiralTag) {
                <$ty>::set_chiral_tag(self, value);
            }
            fn set_chiral_permutation(&mut self, value: Option<u32>) {
                <$ty>::set_chiral_permutation(self, value);
            }
        }
    };
}
parser_atom_access!(Atom);
parser_atom_access!(QueryAtom);
pub trait ParserBond: sealed::BondCarrier {
    fn begin(&self) -> AtomId;
    fn end(&self) -> AtomId;
    fn order(&self) -> BondOrder;
}
impl ParserBond for Bond {
    fn begin(&self) -> AtomId {
        Bond::begin(self)
    }
    fn end(&self) -> AtomId {
        Bond::end(self)
    }
    fn order(&self) -> BondOrder {
        Bond::order(self)
    }
}
impl ParserBond for QueryBond {
    fn begin(&self) -> AtomId {
        self.bond().begin()
    }
    fn end(&self) -> AtomId {
        self.bond().end()
    }
    fn order(&self) -> BondOrder {
        self.bond().order()
    }
}

fn get_bond_ordering<B: ParserBond>(
    atom: AtomId,
    bonds: &[B],
    incident: impl Iterator<Item = (usize, BondId)>,
    ring_closures: &[BondId],
) -> Result<(Vec<BondId>, usize), ParserCarrierError> {
    // BEGIN RDKIT CPP FUNCTION GetBondOrdering
    // RDKit✔️✔️: unsigned int GetBondOrdering(INT_LIST &bondOrdering, const RDKit::RWMol *mol,
    // RDKit✔️✔️:                              const RDKit::Atom *atom) {
    // RDKit✔️✔️:   INT_VECT ringClosures;
    // RDKit✔️✔️:   atom->getPropIfPresent(common_properties::_RingClosures, ringClosures);
    // RDKit✔️✔️:   std::list<SIZET_PAIR> neighbors;
    // RDKit✔️✔️:   neighbors.emplace_back(atom->getIdx(), -1);
    // RDKit✔️✔️:   for (auto nbrIdx : boost::make_iterator_range(mol->getAtomNeighbors(atom))) {
    // RDKit✔️✔️:     const Bond *nbrBond = mol->getBondBetweenAtoms(atom->getIdx(), nbrIdx);
    // RDKit✔️✔️:     if (std::find(ringClosures.begin(), ringClosures.end(),
    // RDKit✔️✔️:                   static_cast<int>(nbrBond->getIdx())) == ringClosures.end()) {
    // RDKit✔️✔️:       neighbors.emplace_back(nbrIdx, nbrBond->getIdx());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   neighbors.sort();
    // RDKit✔️✔️:   auto selfPos = neighbors.begin();
    // RDKit✔️✔️:   if (selfPos->first != atom->getIdx()) {
    // RDKit✔️✔️:     ++selfPos;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   CHECK_INVARIANT(selfPos->first == atom->getIdx(), "weird atom ordering");
    // RDKit✔️✔️:   for (auto neighborIt = neighbors.begin(); neighborIt != neighbors.end();
    // RDKit✔️✔️:        ++neighborIt) {
    // RDKit✔️✔️:     if (neighborIt != selfPos) {
    // RDKit✔️✔️:       bondOrdering.push_back(rdcast<int>(neighborIt->second));
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       bondOrdering.insert(bondOrdering.end(), ringClosures.begin(),
    // RDKit✔️✔️:                           ringClosures.end());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return ringClosures.size();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION GetBondOrdering
    let mut neighbors = vec![(atom.index(), None)];
    for (neighbor, bond) in incident {
        if !ring_closures.contains(&bond) {
            neighbors.push((neighbor, Some(bond)));
        }
    }
    neighbors.sort_by_key(|(neighbor, _)| *neighbor);
    let self_position = neighbors
        .iter()
        .position(|(neighbor, bond)| *neighbor == atom.index() && bond.is_none())
        .ok_or_else(|| {
            ParserCarrierError::Model("SMILES atom is absent from bond ordering".into())
        })?;
    let mut ordering = Vec::with_capacity(neighbors.len().saturating_sub(1) + ring_closures.len());
    for (position, (_, bond)) in neighbors.into_iter().enumerate() {
        if position == self_position {
            ordering.extend(ring_closures.iter().copied());
        } else {
            ordering.push(bond.expect("only the self sentinel has no bond"));
        }
    }
    // Validate that every generated ID still refers to an incident bond. This
    // turns the source invariant into a structured detached-model error.
    if ordering.iter().any(|bond| {
        let bond = &bonds[bond.index()];
        bond.begin() != atom && bond.end() != atom
    }) {
        return Err(ParserCarrierError::Model(
            "SMILES bond ordering contains a non-incident bond".into(),
        ));
    }
    Ok((ordering, ring_closures.len()))
}

pub fn parser_chirality_assignments<
    A: ParserAtom,
    B: ParserBond,
    I: Iterator<Item = (usize, BondId)>,
>(
    atoms: &[A],
    bonds: &[B],
    neighbors: impl Fn(usize) -> I,
    ring_closures_by_atom: &[Vec<BondId>],
    smiles_start_atoms: &[bool],
) -> Result<Vec<(ChiralTag, Option<u32>)>, ParserCarrierError> {
    // BEGIN RDKIT CPP FUNCTION AdjustAtomChiralityFlags
    // RDKit✔️✔️: void AdjustAtomChiralityFlags(RWMol *mol) {
    // RDKit✔️✔️:   for (auto atom : mol->atoms()) {
    // RDKit✔️✔️:     Atom::ChiralType chiralType = atom->getChiralTag();
    // RDKit✔️✔️:     if (chiralType == Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit✔️✔️:         chiralType == Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit✔️✔️:       INT_LIST bondOrdering;
    // RDKit✔️✔️:       unsigned int numClosures = GetBondOrdering(bondOrdering, mol, atom);
    // RDKit✔️✔️:       int nSwaps = atom->getPerturbationOrder(bondOrdering);
    // RDKit✔️✔️:       if (Canon::chiralAtomNeedsTagInversion(
    // RDKit✔️✔️:               *mol, atom, atom->hasProp(common_properties::_SmilesStart),
    // RDKit✔️✔️:               numClosures)) {
    // RDKit✔️✔️:         ++nSwaps;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (nSwaps % 2) {
    // RDKit✔️✔️:         atom->invertChirality();
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (chiralType == Atom::CHI_SQUAREPLANAR ||
    // RDKit✔️✔️:                chiralType == Atom::CHI_TRIGONALBIPYRAMIDAL ||
    // RDKit✔️✔️:                chiralType == Atom::CHI_OCTAHEDRAL) {
    // RDKit✔️✔️:       INT_LIST bonds;
    // RDKit✔️✔️:       GetBondOrdering(bonds, mol, atom);
    // RDKit✔️✔️:
    // RDKit✔️✔️:       unsigned int ref_max = Chirality::getMaxNbors(chiralType);
    // RDKit✔️✔️:
    // RDKit✔️✔️:       // insert (-1) for hydrogens or missing ligands, where these are placed
    // RDKit✔️✔️:       // depends on if it is the first atom or not
    // RDKit✔️✔️:       if (bonds.size() < ref_max) {
    // RDKit✔️✔️:         if (atom->hasProp(common_properties::_SmilesStart)) {
    // RDKit✔️✔️:           bonds.insert(bonds.begin(), ref_max - bonds.size(), -1);
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           bonds.insert(++bonds.begin(), ref_max - bonds.size(), -1);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       atom->setProp(common_properties::_chiralPermutation,
    // RDKit✔️✔️:                     Chirality::getChiralPermutation(atom, bonds, true));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION AdjustAtomChiralityFlags
    let mut invert = vec![false; atoms.len()];
    let mut nontetrahedral_permutations = vec![None; atoms.len()];
    for atom_index in 0..atoms.len() {
        let atom = &atoms[atom_index];
        let atom_id = AtomId::new(atom_index);
        match atom.chiral_tag() {
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw => {
                let (ordering, num_closures) = get_bond_ordering(
                    atom_id,
                    bonds,
                    neighbors(atom_index),
                    &ring_closures_by_atom[atom_index],
                )?;
                let storage_order = neighbors(atom_index)
                    .map(|(_, bond)| bond)
                    .collect::<Vec<_>>();
                let mut swaps = crate::parser_stereo_order::count_swaps_to_interconvert(
                    &ordering,
                    &storage_order,
                )
                .map_err(|_| {
                    ParserCarrierError::Model(
                        "SMILES and storage bond orderings are not permutations".into(),
                    )
                })?;
                let unsaturated = neighbors(atom_index).any(|(_, bond)| {
                    crate::parser_stereo_order::bond_order_as_double(bonds[bond.index()].order())
                        > 1.0
                });
                if crate::parser_stereo_order::chiral_atom_needs_tag_inversion(
                    neighbors(atom_index).count(),
                    atom.explicit_hydrogens(),
                    smiles_start_atoms[atom_index],
                    crate::parser_stereo_order::atom_has_fourth_valence(
                        atom.explicit_hydrogens(),
                        false,
                    ),
                    num_closures,
                    unsaturated,
                ) {
                    swaps += 1;
                }
                invert[atom_index] = swaps % 2 == 1;
            }
            ChiralTag::SquarePlanar | ChiralTag::TrigonalBipyramidal | ChiralTag::Octahedral => {
                let (ordering, _) = get_bond_ordering(
                    atom_id,
                    bonds,
                    neighbors(atom_index),
                    &ring_closures_by_atom[atom_index],
                )?;
                let mut probe = ordering.into_iter().map(Some).collect::<Vec<_>>();
                crate::parser_stereo_order::insert_implicit_nontetrahedral_neighbors(
                    &mut probe,
                    atom.chiral_tag(),
                    smiles_start_atoms[atom_index],
                );
                let incident = neighbors(atom_index)
                    .map(|(_, bond)| bond)
                    .collect::<Vec<_>>();
                nontetrahedral_permutations[atom_index] = Some(
                    crate::parser_stereo_order::nontetrahedral_chiral_permutation(
                        atom.chiral_permutation().unwrap_or(0),
                        atom.chiral_tag(),
                        bonds.len(),
                        &incident,
                        &probe,
                        true,
                    )
                    .map_err(|error| ParserCarrierError::Model(error.to_string()))?,
                );
            }
            _ => {}
        }
    }
    Ok(atoms
        .iter()
        .zip(invert)
        .zip(nontetrahedral_permutations)
        .map(|((atom, invert), permutation)| {
            let tag = if invert {
                crate::parser_stereo_order::invert_tetrahedral_tag(atom.chiral_tag())
            } else {
                atom.chiral_tag()
            };
            (tag, permutation.or(atom.chiral_permutation()))
        })
        .collect())
}

pub fn cleanup_parser_atoms<A: ParserAtom>(atoms: &mut [A]) -> Result<(), ParserCarrierError> {
    // RDKit✔️✔️: void CleanupAfterParsing(RWMol *mol) {
    // RDKit✔️✔️:   PRECONDITION(mol, "no molecule");
    // RDKit✔️✔️:   for (auto atom : mol->atoms()) {
    // RDKit✔️✔️:     atom->clearProp(common_properties::_RingClosures);
    // RDKit✔️✔️:     atom->clearProp(common_properties::_SmilesStart);
    // RDKit✔️✔️:     std::string label;
    // RDKit✔️✔️:     if (atom->getAtomicNum() == 0 &&
    // RDKit✔️✔️:         atom->getPropIfPresent(common_properties::atomLabel, label)) {
    // RDKit✔️✔️:       // marvinsketch can output higher labels than _AP1 and _AP2, but they
    // RDKit✔️✔️:       // aren't part of the MOL file spec so we don't treat them as attachment
    // RDKit✔️✔️:       // points
    // RDKit✔️✔️:       if (label == "_AP1") {
    // RDKit✔️✔️:         atom->setProp(common_properties::_fromAttachPoint, 1);
    // RDKit✔️✔️:       } else if (label == "_AP2") {
    // RDKit✔️✔️:         atom->setProp(common_properties::_fromAttachPoint, 2);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (auto bond : mol->bonds()) {
    // RDKit✔️✔️:     bond->clearProp(common_properties::_unspecifiedOrder);
    // RDKit✔️✔️:     bond->clearProp("_cxsmilesBondIdx");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (auto sg : RDKit::getSubstanceGroups(*mol)) {
    // RDKit✔️✔️:     sg.clearProp("_cxsmilesindex");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!Chirality::getAllowNontetrahedralChirality()) {
    // RDKit✔️✔️:     bool needWarn = false;
    // RDKit✔️✔️:     for (auto atom : mol->atoms()) {
    // RDKit✔️✔️:       if (atom->hasProp(common_properties::_chiralPermutation)) {
    // RDKit✔️✔️:         needWarn = true;
    // RDKit✔️✔️:         atom->clearProp(common_properties::_chiralPermutation);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (atom->getChiralTag() > Atom::ChiralType::CHI_OTHER) {
    // RDKit✔️✔️:         needWarn = true;
    // RDKit✔️✔️:         atom->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (needWarn) {
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:           << "ignoring non-tetrahedral stereo specification since setAllowNontetrahedralChirality() is false."
    // RDKit✔️✔️:           << std::endl;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Linear carrier passes and constant-sized property keys preserve source cost.
    for atom in atoms.iter_mut() {
        atom.clear_prop("_RingClosures")?;
        atom.clear_prop("_SmilesStart")?;
        if atom.atomic_number() == 0 {
            let label = atom
                .prop("atomLabel")
                .map(crate::property_value_to_string)
                .transpose()?;
            match label.as_ref().map(|value| value.as_bytes()) {
                Some(b"_AP1") => atom.set_attachment_point(1)?,
                Some(b"_AP2") => atom.set_attachment_point(2)?,
                _ => {}
            }
        }
    }
    // This is the source first atom pass only. Both destination callers clean
    // their actual bonds and copied SGroups before reaching the shared final
    // non-tetrahedral pass; keeping it here reordered fallible property reads.
    Ok(())
}

/// Final source cleanup pass, reached only after atom, bond and copied-group cleanup.
#[doc(hidden)]
pub fn cleanup_parser_nontetrahedral_atoms<A: ParserAtom>(
    atoms: &mut [A],
) -> Result<(), ParserCarrierError> {
    // RDKit✔️✔️:   if (!Chirality::getAllowNontetrahedralChirality()) {
    // RDKit✔️✔️:     bool needWarn = false;
    // RDKit✔️✔️:     for (auto atom : mol->atoms()) {
    // RDKit✔️✔️:       if (atom->hasProp(common_properties::_chiralPermutation)) {
    // RDKit✔️✔️:         needWarn = true;
    // RDKit✔️✔️:         atom->clearProp(common_properties::_chiralPermutation);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (atom->getChiralTag() > Atom::ChiralType::CHI_OTHER) {
    // RDKit✔️✔️:         needWarn = true;
    // RDKit✔️✔️:         atom->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (needWarn) {
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:           << "ignoring non-tetrahedral stereo specification since setAllowNontetrahedralChirality() is false."
    // RDKit✔️✔️:           << std::endl;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // Raw dictionary presence includes wrong tags, because hasProp does not
    // cast. The existing typed permutation field is the parser's detached
    // projection of that same native property and also requires clearing.
    // Each reached clearProp retains canonical computed-entry errors before
    // changing the tag; the one warning follows the entire successful pass.
    // Cost: one borrowed atom pass, fixed property keys, O(1) temporary state;
    // no buffering or full-graph clone. The environment policy is read here,
    // after both caller bond and copied-SGroup loops, exactly source order.
    if !crate::nontetrahedral_enabled() {
        let mut need_warn = false;
        for atom in atoms {
            if atom.prop("_chiralPermutation").is_some() || atom.chiral_permutation().is_some() {
                need_warn = true;
                atom.clear_prop("_chiralPermutation")?;
                atom.set_chiral_permutation(None);
            }
            if atom.chiral_tag().rdkit_code() > ChiralTag::Other.rdkit_code() {
                need_warn = true;
                atom.set_chiral_tag(ChiralTag::Unspecified);
            }
        }
        if need_warn {
            eprintln!(
                "ignoring non-tetrahedral stereo specification since setAllowNontetrahedralChirality() is false."
            );
        }
    }
    Ok(())
}
pub fn cleanup_parser_substance_groups(
    groups: &[cosmolkit_model::SubstanceGroup],
) -> Result<(), cosmolkit_model::MoleculePropertyError> {
    // RDKit✔️❌:   for (auto sg : RDKit::getSubstanceGroups(*mol)) {
    // RDKit✔️❌:     sg.clearProp("_cxsmilesindex");
    // RDKit✔️❌:   }
    // C++ auto (without &) copies each group. Clearing the local copy cannot
    // change the original carrier. Clone one group at a time, matching the
    // source per-group copy and linear property cost. Canonical tree storage
    // plus source-order key copies add allocations versus the source Dict.
    for group in groups {
        let mut local = group.clone();
        local.clear_prop("_cxsmilesindex")?;
    }
    Ok(())
}
