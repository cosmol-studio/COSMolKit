use crate::{ReactionInput, ReactionProductError, ReactionRole, ReactionRowOrigin};
use cosmolkit_model::{
    Atom, AtomId, Bond, BondId, BondQueryPredicate, BondSpec, NeighborRef, QueryGraph, QueryNode,
    TopologyBlock,
};
use cosmolkit_types::{BondDirection, BondOrder, ChiralTag};
use std::collections::{BTreeMap, BTreeSet, VecDeque};

/// Canonical detached product rows under construction, with temporary ordered
/// neighbor indices. These indices are discarded when canonical adjacency is
/// finalized; no alternative query AST or live runtime authority exists here.
pub(crate) struct ProductBuilder {
    pub topology: TopologyBlock,
    pub neighbors: Vec<Vec<NeighborRef>>,
    pub bookmarks: BTreeMap<i32, Vec<usize>>,
    pub atom_origins: Vec<Option<ReactionRowOrigin<AtomId>>>,
    pub bond_origins: Vec<Option<ReactionRowOrigin<BondId>>>,
}

impl ProductBuilder {
    fn add_atom(&mut self, atom: Atom, origin: Option<ReactionRowOrigin<AtomId>>) -> usize {
        let row = self.topology.atoms.len();
        self.topology.atoms.push(atom.with_id(AtomId::new(row)));
        self.neighbors.push(Vec::new());
        self.atom_origins.push(origin);
        row
    }

    pub(crate) fn bond_between(&self, begin: usize, end: usize) -> Option<BondId> {
        self.neighbors[begin]
            .iter()
            .find(|n| n.atom_index == end)
            .map(|n| n.bond)
    }

    fn add_bond(
        &mut self,
        begin: usize,
        end: usize,
        order: BondOrder,
        origin: Option<ReactionRowOrigin<BondId>>,
    ) -> Result<BondId, ReactionProductError> {
        // RDKit❗✔️: unsigned int RWMol::addBond(unsigned int atomIdx1, unsigned int atomIdx2,
        // RDKit❗✔️:                             Bond::BondType bondType) {
        // RDKit❗✔️:   // if the atom indices are bad, the next two calls will catch that.
        // RDKit❗✔️:   auto beginAtom = getAtomWithIdx(atomIdx1);
        // RDKit❗✔️:   auto endAtom = getAtomWithIdx(atomIdx2);
        // RDKit❗✔️:   PRECONDITION(atomIdx1 != atomIdx2, "attempt to add self-bond");
        // RDKit❗✔️:   PRECONDITION(!(boost::edge(atomIdx1, atomIdx2, d_graph).second),
        // RDKit❗✔️:                "bond already exists");
        // RDKit❗✔️:
        // RDKit❗✔️:   auto *b = new Bond(bondType);
        // RDKit❗✔️:   b->setOwningMol(this);
        // RDKit❗✔️:   if (bondType == Bond::AROMATIC) {
        // RDKit❗✔️:     b->setIsAromatic(1);
        // RDKit❗✔️:     //
        // RDKit❗✔️:     // assume that aromatic bonds connect aromatic atoms
        // RDKit❗✔️:     //   This is relevant for file formats like MOL, where there
        // RDKit❗✔️:     //   is no such thing as an aromatic atom, but bonds can be
        // RDKit❗✔️:     //   marked aromatic.
        // RDKit❗✔️:     //
        // RDKit❗✔️:     beginAtom->setIsAromatic(1);
        // RDKit❗✔️:     endAtom->setIsAromatic(1);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   auto [which, ok] = boost::add_edge(atomIdx1, atomIdx2, d_graph);
        // RDKit❗✔️:   d_graph[which] = b;
        // RDKit❗✔️:   ++numBonds;
        // RDKit❗✔️:   b->setIdx(numBonds - 1);
        // RDKit❗✔️:   b->setBeginAtomIdx(atomIdx1);
        // RDKit❗✔️:   b->setEndAtomIdx(atomIdx2);
        // RDKit❗✔️:
        // RDKit❗✔️:   // the valence values on the begin and end atoms need to be updated:
        // RDKit❗✔️:   beginAtom->clearPropertyCache();
        // RDKit❗✔️:   endAtom->clearPropertyCache();
        // RDKit❗✔️:
        // RDKit❗✔️:   // we're in a batch edit, and at least one of the bond ends is scheduled
        // RDKit❗✔️:   // for deletion, so mark the new bond for deletion too:
        // RDKit❗✔️:   if (dp_delAtoms &&
        // RDKit❗✔️:       ((atomIdx1 < dp_delAtoms->size() && dp_delAtoms->test(atomIdx1)) ||
        // RDKit❗✔️:        (atomIdx2 < dp_delAtoms->size() && dp_delAtoms->test(atomIdx2)))) {
        // RDKit❗✔️:     if (dp_delBonds->size() < numBonds) {
        // RDKit❗✔️:       dp_delBonds->resize(numBonds);
        // RDKit❗✔️:     }
        // RDKit❗✔️:     dp_delBonds->set(numBonds - 1);
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   return numBonds;
        // RDKit❗✔️: }
        // Ordered O(degree) edge check and amortized row appends match the
        // source adjacency-list graph; no full topology/adjacency rebuild.
        if begin >= self.topology.atoms.len()
            || end >= self.topology.atoms.len()
            || begin == end
            || self.bond_between(begin, end).is_some()
        {
            return Err(invariant(
                "addBond",
                "invalid endpoints or duplicate bond",
                None,
                Some(begin),
                None,
            ));
        }
        let id = BondId::new(self.topology.bonds.len());
        let aromatic = order == BondOrder::Aromatic;
        if aromatic {
            self.topology.atoms[begin].set_aromatic(true);
            self.topology.atoms[end].set_aromatic(true);
        }
        let bond = Bond::from_spec(
            id,
            BondSpec::new(AtomId::new(begin), AtomId::new(end), order).with_aromatic(aromatic),
        );
        self.topology.bonds.push(bond);
        self.neighbors[begin].push(NeighborRef {
            atom_index: end,
            bond: id,
        });
        self.neighbors[end].push(NeighborRef {
            atom_index: begin,
            bond: id,
        });
        self.bond_origins.push(origin);
        Ok(id)
    }
}

pub(crate) struct ReactantProductMapping {
    pub mapped: Vec<bool>,
    pub skipped: Vec<bool>,
    pub reactant_to_product: BTreeMap<usize, Vec<usize>>,
    pub product_to_reactant: BTreeMap<usize, usize>,
    pub product_atom_bond: BTreeMap<usize, usize>,
    pub template_bonds: BTreeSet<(u32, u32)>,
}

pub(crate) fn invariant(
    stage: &'static str,
    detail: &'static str,
    reactant_atom: Option<usize>,
    product_atom: Option<usize>,
    bond: Option<BondId>,
) -> ReactionProductError {
    ReactionProductError::Invariant {
        stage,
        detail,
        reactant_atom,
        product_atom,
        bond,
    }
}

pub(crate) trait ReactionAtomPropertyRead {
    fn property_id(&self) -> AtomId;
    fn property(&self, key: &str) -> Option<&cosmolkit_model::PropertyValue>;
    fn typed_inversion(&self) -> Option<i32>;
}
impl ReactionAtomPropertyRead for Atom {
    fn property_id(&self) -> AtomId {
        self.id()
    }
    fn property(&self, key: &str) -> Option<&cosmolkit_model::PropertyValue> {
        self.prop(key)
    }
    fn typed_inversion(&self) -> Option<i32> {
        self.mol_inversion_flag()
    }
}
impl ReactionAtomPropertyRead for cosmolkit_model::QueryAtom {
    fn property_id(&self) -> AtomId {
        self.id()
    }
    fn property(&self, key: &str) -> Option<&cosmolkit_model::PropertyValue> {
        self.prop(key)
    }
    fn typed_inversion(&self) -> Option<i32> {
        self.mol_inversion_flag()
    }
}

pub(crate) fn int_prop(
    atom: &impl ReactionAtomPropertyRead,
    key: &'static str,
) -> Result<Option<i32>, ReactionProductError> {
    atom.property(key)
        .map(|value| {
            cosmolkit_core::property_value_to_int(value).map_err(|source| {
                ReactionProductError::PropertyInt {
                    atom: atom.property_id(),
                    key,
                    source,
                }
            })
        })
        .transpose()
}

pub(crate) fn uint_prop(
    atom: &impl ReactionAtomPropertyRead,
    key: &'static str,
) -> Result<Option<u32>, ReactionProductError> {
    atom.property(key)
        .map(|value| {
            cosmolkit_core::property_value_to_uint(value).map_err(|source| {
                ReactionProductError::PropertyUInt {
                    atom: atom.property_id(),
                    key,
                    source,
                }
            })
        })
        .transpose()
}

pub(crate) fn inversion_flag(
    atom: &impl ReactionAtomPropertyRead,
) -> Result<Option<i32>, ReactionProductError> {
    match atom.typed_inversion() {
        Some(value) => Ok(Some(value)),
        None => int_prop(atom, "molInversionFlag"),
    }
}

pub(crate) fn source_u32(kind: &'static str, index: usize) -> Result<u32, ReactionProductError> {
    u32::try_from(index).map_err(|_| ReactionProductError::RowOverflow { kind, index })
}

pub(crate) fn update_implicit_atom_properties(product: &mut Atom, reactant: &Atom) {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: updateImplicitAtomProperties
    // RDKit❗✔️: void updateImplicitAtomProperties(Atom *prodAtom, const Atom *reactAtom) {
    // RDKit❗✔️:   PRECONDITION(prodAtom, "no product atom");
    // RDKit❗✔️:   PRECONDITION(reactAtom, "no reactant atom");
    // RDKit❗✔️:   if (prodAtom->getAtomicNum() != reactAtom->getAtomicNum()) {
    // RDKit❗✔️:     // if we changed atom identity all bets are off, just
    // RDKit❗✔️:     // return
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (!prodAtom->hasProp(common_properties::_QueryFormalCharge)) {
    // RDKit❗✔️:     prodAtom->setFormalCharge(reactAtom->getFormalCharge());
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (!prodAtom->hasProp(common_properties::_QueryIsotope)) {
    // RDKit❗✔️:     prodAtom->setIsotope(reactAtom->getIsotope());
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (!prodAtom->hasProp(common_properties::_ReactionDegreeChanged)) {
    // RDKit❗✔️:     if (!prodAtom->hasProp(common_properties::_QueryHCount)) {
    // RDKit❗✔️:       prodAtom->setNumExplicitHs(reactAtom->getNumExplicitHs());
    // RDKit❗✔️:       prodAtom->setNoImplicit(reactAtom->getNoImplicit());
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    if product.atomic_number() != reactant.atomic_number() {
        return;
    }
    if product.prop("_QueryFormalCharge").is_none() {
        product.set_formal_charge(reactant.formal_charge());
    }
    if product.prop("_QueryIsotope").is_none() {
        product.set_isotope(reactant.isotope());
    }
    if product.prop("_ReactionDegreeChanged").is_none() && product.prop("_QueryHCount").is_none() {
        product.set_explicit_hydrogens(reactant.explicit_hydrogens());
        product.set_no_implicit(reactant.no_implicit());
    }
}

pub(crate) fn update_from_template(
    template: &impl ReactionAtomPropertyRead,
    atom: &mut Atom,
) -> Result<bool, ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: updatePropsFromImplicitProps
    // RDKit❗✔️: bool updatePropsFromImplicitProps(Atom *templateAtom, Atom *atom) {
    // RDKit❗✔️:   PRECONDITION(templateAtom, "no atom");
    // RDKit❗✔️:   PRECONDITION(atom, "no atom");
    // RDKit❗✔️:   bool res = false;
    // RDKit❗✔️:   int val;
    // RDKit❗✔️:   if (templateAtom->getPropIfPresent(common_properties::_QueryFormalCharge,
    // RDKit❗✔️:                                      val) &&
    // RDKit❗✔️:       val != atom->getFormalCharge()) {
    // RDKit❗✔️:     atom->setFormalCharge(val);
    // RDKit❗✔️:     res = true;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   unsigned int uval;
    // RDKit❗✔️:   if (templateAtom->getPropIfPresent(common_properties::_QueryHCount, uval)) {
    // RDKit❗✔️:     if (!atom->getNoImplicit() || atom->getNumExplicitHs() != uval) {
    // RDKit❗✔️:       atom->setNumExplicitHs(uval);
    // RDKit❗✔️:       atom->setNoImplicit(true);  // this was github #1544
    // RDKit❗✔️:       res = true;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (templateAtom->getPropIfPresent(common_properties::_QueryMass, uval)) {
    // RDKit❗✔️:     // FIX: technically should do something with this
    // RDKit❗✔️:     // atom->setMass(val);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (templateAtom->getPropIfPresent(common_properties::_QueryIsotope, uval) &&
    // RDKit❗✔️:       uval != atom->getIsotope()) {
    // RDKit❗✔️:     atom->setIsotope(uval);
    // RDKit❗✔️:     res = true;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    let mut changed = false;
    if let Some(value) = int_prop(template, "_QueryFormalCharge")?
        && value != i32::from(atom.formal_charge())
    {
        atom.set_formal_charge(value as i8);
        changed = true;
    }
    if let Some(value) = uint_prop(template, "_QueryHCount")?
        && (!atom.no_implicit() || u32::from(atom.explicit_hydrogens()) != value)
    {
        atom.set_explicit_hydrogens(value as u8);
        atom.set_no_implicit(true);
        changed = true;
    }
    // Source reads/converts this value even though it does not set mass.
    let _mass = uint_prop(template, "_QueryMass")?;
    if let Some(value) = uint_prop(template, "_QueryIsotope")?
        && u32::from(atom.isotope().unwrap_or(0)) != value
    {
        atom.set_isotope(Some(value as u16));
        changed = true;
    }
    Ok(changed)
}

pub(crate) fn convert_template(
    template: &QueryGraph,
    template_index: usize,
) -> Result<ProductBuilder, ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: convertTemplateToMol
    // RDKit❗✔️: RWMOL_SPTR convertTemplateToMol(const ROMOL_SPTR prodTemplateSptr) {
    // RDKit❗✔️:   const ROMol *prodTemplate = prodTemplateSptr.get();
    // RDKit❗✔️:   auto *res = new RWMol();
    // RDKit❗✔️:
    // RDKit❗✔️:   // --------- --------- --------- --------- --------- ---------
    // RDKit❗✔️:   // Initialize by making a copy of the product template as a normal molecule.
    // RDKit❗✔️:   // NOTE that we can't just use a normal copy because we do not want to end up
    // RDKit❗✔️:   // with query atoms or bonds in the product.
    // RDKit❗✔️:
    // RDKit❗✔️:   // copy in the atoms:
    // RDKit❗✔️:   ROMol::ATOM_ITER_PAIR atItP = prodTemplate->getVertices();
    // RDKit❗✔️:   while (atItP.first != atItP.second) {
    // RDKit❗✔️:     const Atom *oAtom = (*prodTemplate)[*(atItP.first++)];
    // RDKit❗✔️:     auto *newAtom = new Atom(*oAtom);
    // RDKit❗✔️:     res->addAtom(newAtom, false, true);
    // RDKit❗✔️:     int mapNum;
    // RDKit❗✔️:     if (newAtom->getPropIfPresent(common_properties::molAtomMapNumber,
    // RDKit❗✔️:                                   mapNum)) {
    // RDKit❗✔️:       // set bookmarks for the mapped atoms:
    // RDKit❗✔️:       res->setAtomBookmark(newAtom, mapNum);
    // RDKit❗✔️:       // now clear the molAtomMapNumber property so that it doesn't
    // RDKit❗✔️:       // end up in the products (this was bug 3140490):
    // RDKit❗✔️:       newAtom->clearProp(common_properties::molAtomMapNumber);
    // RDKit❗✔️:       newAtom->setProp<int>(common_properties::reactionMapNum, mapNum);
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     newAtom->setChiralTag(Atom::CHI_UNSPECIFIED);
    // RDKit❗✔️:     // if the product-template atom has the inversion flag set
    // RDKit❗✔️:     // to 4 (=SET), then bring its stereochem over, otherwise we'll
    // RDKit❗✔️:     // ignore it:
    // RDKit❗✔️:     int iFlag;
    // RDKit❗✔️:     if (oAtom->getPropIfPresent(common_properties::molInversionFlag, iFlag)) {
    // RDKit❗✔️:       if (iFlag == 4) {
    // RDKit❗✔️:         newAtom->setChiralTag(oAtom->getChiralTag());
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     // check for properties we need to set:
    // RDKit❗✔️:     updatePropsFromImplicitProps(newAtom, newAtom);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // and the bonds:
    // RDKit❗✔️:   ROMol::BOND_ITER_PAIR bondItP = prodTemplate->getEdges();
    // RDKit❗✔️:   while (bondItP.first != bondItP.second) {
    // RDKit❗✔️:     const Bond *oldB = (*prodTemplate)[*(bondItP.first++)];
    // RDKit❗✔️:     unsigned int bondIdx;
    // RDKit❗✔️:     bondIdx = res->addBond(oldB->getBeginAtomIdx(), oldB->getEndAtomIdx(),
    // RDKit❗✔️:                            oldB->getBondType()) -
    // RDKit❗✔️:               1;
    // RDKit❗✔️:     // make sure we don't lose the bond dir information:
    // RDKit❗✔️:     Bond *newB = res->getBondWithIdx(bondIdx);
    // RDKit❗✔️:     newB->setBondDir(oldB->getBondDir());
    // RDKit❗✔️:     // Special case/hack:
    // RDKit❗✔️:     //  The product has been processed by the SMARTS parser.
    // RDKit❗✔️:     //  The SMARTS parser tags unspecified bonds as single, but then adds
    // RDKit❗✔️:     //  a query so that they match single or double
    // RDKit❗✔️:     //  This caused Issue 1748846
    // RDKit❗✔️:     //   http://sourceforge.net/tracker/index.php?func=detail&aid=1748846&group_id=160139&atid=814650
    // RDKit❗✔️:     //  We need to fix that little problem now:
    // RDKit❗✔️:     if (oldB->hasQuery()) {
    // RDKit❗✔️:       //  remember that the product has been processed by the SMARTS parser.
    // RDKit❗✔️:       std::string queryDescription = oldB->getQuery()->getDescription();
    // RDKit❗✔️:       if (queryDescription == "BondOr" && oldB->getBondType() == Bond::SINGLE) {
    // RDKit❗✔️:         //  We need to fix that little problem now:
    // RDKit❗✔️:         if (newB->getBeginAtom()->getIsAromatic() &&
    // RDKit❗✔️:             newB->getEndAtom()->getIsAromatic()) {
    // RDKit❗✔️:           newB->setBondType(Bond::AROMATIC);
    // RDKit❗✔️:           newB->setIsAromatic(true);
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           newB->setBondType(Bond::SINGLE);
    // RDKit❗✔️:           newB->setIsAromatic(false);
    // RDKit❗✔️:         }
    // RDKit❗✔️:       } else if (queryDescription == "BondNull") {
    // RDKit❗✔️:         newB->setProp(common_properties::NullBond, 1);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     // Double bond stereo: if a double bond has at least one bond on each side,
    // RDKit❗✔️:     // and none of those has a direction, then mark it as unknown stereo to have
    // RDKit❗✔️:     // it reset later on. This has to be done before the reactant atoms are
    // RDKit❗✔️:     // added,
    // RDKit❗✔️:     if (oldB->getBondType() == Bond::BondType::DOUBLE) {
    // RDKit❗✔️:       const Atom *startAtom = oldB->getBeginAtom();
    // RDKit❗✔️:       const Atom *endAtom = oldB->getEndAtom();
    // RDKit❗✔️:
    // RDKit❗✔️:       if (startAtom->getDegree() > 1 && endAtom->getDegree() > 1 &&
    // RDKit❗✔️:           (Chirality::getNeighboringDirectedBond(*prodTemplate, startAtom) ==
    // RDKit❗✔️:                nullptr ||
    // RDKit❗✔️:            Chirality::getNeighboringDirectedBond(*prodTemplate, endAtom) ==
    // RDKit❗✔️:                nullptr)) {
    // RDKit❗✔️:         newB->setProp(_UnknownStereoRxnBond, 1);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     // copy properties over:
    // RDKit❗✔️:     bool preserveExisting = true;
    // RDKit❗✔️:     newB->updateProps(*static_cast<const RDProps *>(oldB), preserveExisting);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return RWMOL_SPTR(res);
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    // Product source metadata starts empty; source template coordinates and
    // molecule properties are not implicitly copied into a new product.
    let mut result = ProductBuilder {
        topology: TopologyBlock::default(),
        neighbors: Vec::with_capacity(template.num_atoms()),
        bookmarks: BTreeMap::new(),
        atom_origins: Vec::with_capacity(template.num_atoms()),
        bond_origins: Vec::with_capacity(template.num_bonds()),
    };
    for (row, query_atom) in template.atoms().iter().enumerate() {
        let mut atom = query_atom.try_to_atom()?;
        if let Some(map) =
            crate::validation::atom_map(query_atom, ReactionRole::Product, template_index)?
        {
            result.bookmarks.entry(map).or_default().push(row);
            atom.set_atom_map(None);
            atom.clear_prop("molAtomMapNumber");
            atom.set_prop("old_mapno", map)?;
        }
        atom.set_chiral_tag(ChiralTag::Unspecified);
        if inversion_flag(query_atom)? == Some(4) {
            atom.set_chiral_tag(query_atom.chiral_tag());
        }
        // Source uses one object as both template and output. Reads precede
        // writes and keys are disjoint; the borrowed original has identical
        // property values, avoiding a second property-value clone.
        update_from_template(query_atom, &mut atom)?;
        result.add_atom(atom, None);
    }
    for query_bond in template.bonds() {
        let old = query_bond.bond();
        let id = result.add_bond(old.begin().index(), old.end().index(), old.order(), None)?;
        let mut root = query_bond.predicate();
        while let QueryNode::Not(child) = root {
            root = child;
        }
        {
            let bond = &mut result.topology.bonds[id.index()];
            bond.set_direction(old.direction());
            if !query_bond.predicate_is_carrier_derived() {
                if matches!(root, QueryNode::Or(_)) && old.order() == BondOrder::Single {
                    let aromatic = result.topology.atoms[old.begin().index()].is_aromatic()
                        && result.topology.atoms[old.end().index()].is_aromatic();
                    bond.set_order(if aromatic {
                        BondOrder::Aromatic
                    } else {
                        BondOrder::Single
                    });
                    bond.set_aromatic(aromatic);
                } else if matches!(root, QueryNode::Predicate(BondQueryPredicate::Any)) {
                    bond.set_prop("NullBond", 1)?;
                }
            }
            if old.order() == BondOrder::Double {
                let a = old.begin().index();
                let b = old.end().index();
                let directed = |at: usize| {
                    template.adjacency()[at].iter().any(|&(_, bond_index)| {
                        let bond = template.bonds()[bond_index].bond();
                        bond.order() != BondOrder::Double
                            && matches!(
                                bond.direction(),
                                BondDirection::EndDownRight | BondDirection::EndUpRight
                            )
                    })
                };
                if template.adjacency()[a].len() > 1
                    && template.adjacency()[b].len() > 1
                    && (!directed(a) || !directed(b))
                {
                    bond.set_prop("_UnknownStereoRxnBond", 1)?;
                }
            }
            bond.update_properties_from(old, true);
        }
    }
    Ok(result)
}

pub(crate) fn atom_mappings(
    matched: &[usize],
    template: &QueryGraph,
    template_index: usize,
    product: &ProductBuilder,
    reactant_atoms: usize,
) -> Result<ReactantProductMapping, ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: getAtomMappingsReactantProduct
    // RDKit❗✔️: ReactantProductAtomMapping *getAtomMappingsReactantProduct(
    // RDKit❗✔️:     const MatchVectType &match, const ROMol &reactantTemplate,
    // RDKit❗✔️:     RWMOL_SPTR product, unsigned numReactAtoms) {
    // RDKit❗✔️:   auto *mapping = new ReactantProductAtomMapping(numReactAtoms);
    // RDKit❗✔️:
    // RDKit❗✔️:   // keep track of which mapped atoms in the reactant template are bonded to
    // RDKit❗✔️:   // each other.
    // RDKit❗✔️:   // This is part of the fix for #1387
    // RDKit❗✔️:   {
    // RDKit❗✔️:     ROMol::EDGE_ITER firstB, lastB;
    // RDKit❗✔️:     boost::tie(firstB, lastB) = reactantTemplate.getEdges();
    // RDKit❗✔️:     while (firstB != lastB) {
    // RDKit❗✔️:       const Bond *bond = reactantTemplate[*firstB];
    // RDKit❗✔️:       // this will put in pairs with 0s for things that aren't mapped, but we
    // RDKit❗✔️:       // don't care about that
    // RDKit❗✔️:       int a1mapidx = bond->getBeginAtom()->getAtomMapNum();
    // RDKit❗✔️:       int a2mapidx = bond->getEndAtom()->getAtomMapNum();
    // RDKit❗✔️:       if (a1mapidx > a2mapidx) {
    // RDKit❗✔️:         std::swap(a1mapidx, a2mapidx);
    // RDKit❗✔️:       }
    // RDKit❗✔️:       mapping->reactantTemplateAtomBonds[std::make_pair(a1mapidx, a2mapidx)] =
    // RDKit❗✔️:           1;
    // RDKit❗✔️:       ++firstB;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   for (const auto &i : match) {
    // RDKit❗✔️:     const Atom *templateAtom = reactantTemplate.getAtomWithIdx(i.first);
    // RDKit❗✔️:     int molAtomMapNumber;
    // RDKit❗✔️:     if (templateAtom->getPropIfPresent(common_properties::molAtomMapNumber,
    // RDKit❗✔️:                                        molAtomMapNumber)) {
    // RDKit❗✔️:       if (product->hasAtomBookmark(molAtomMapNumber)) {
    // RDKit❗✔️:         RWMol::ATOM_PTR_LIST atomIdxs =
    // RDKit❗✔️:             product->getAllAtomsWithBookmark(molAtomMapNumber);
    // RDKit❗✔️:         for (auto a : atomIdxs) {
    // RDKit❗✔️:           unsigned int pIdx = a->getIdx();
    // RDKit❗✔️:           mapping->reactProdAtomMap[i.second].push_back(pIdx);
    // RDKit❗✔️:           mapping->mappedAtoms[i.second] = 1;
    // RDKit❗✔️:           CHECK_INVARIANT(pIdx < product->getNumAtoms(), "yikes!");
    // RDKit❗✔️:           mapping->prodReactAtomMap[pIdx] = i.second;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         // this skippedAtom has an atomMapNumber, but it's not in this product
    // RDKit❗✔️:         // (it's either in another product or it's not mapped at all).
    // RDKit❗✔️:         mapping->skippedAtoms[i.second] = 1;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       // This skippedAtom appears in the match, but not in a product:
    // RDKit❗✔️:       mapping->skippedAtoms[i.second] = 1;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return mapping;
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    let mut mapping = ReactantProductMapping {
        mapped: vec![false; reactant_atoms],
        skipped: vec![false; reactant_atoms],
        reactant_to_product: BTreeMap::new(),
        product_to_reactant: BTreeMap::new(),
        product_atom_bond: BTreeMap::new(),
        template_bonds: BTreeSet::new(),
    };
    for bond in template.bonds() {
        let mut pair = [
            crate::validation::atom_map(
                &template.atoms()[bond.begin().index()],
                ReactionRole::Reactant,
                template_index,
            )?
            .unwrap_or(0),
            crate::validation::atom_map(
                &template.atoms()[bond.end().index()],
                ReactionRole::Reactant,
                template_index,
            )?
            .unwrap_or(0),
        ];
        pair.sort();
        mapping
            .template_bonds
            .insert((pair[0] as u32, pair[1] as u32));
    }
    for (query_row, &target_row) in matched.iter().enumerate() {
        if query_row >= template.num_atoms() || target_row >= reactant_atoms {
            return Err(invariant(
                "getAtomMappingsReactantProduct",
                "matched row out of range",
                Some(target_row),
                None,
                None,
            ));
        }
        let map = crate::validation::atom_map(
            &template.atoms()[query_row],
            ReactionRole::Reactant,
            template_index,
        )?;
        if let Some(products) = map.and_then(|map| product.bookmarks.get(&map)) {
            for &product_row in products {
                mapping
                    .reactant_to_product
                    .entry(target_row)
                    .or_default()
                    .push(product_row);
                mapping.mapped[target_row] = true;
                if product_row >= product.topology.atoms.len() {
                    return Err(invariant(
                        "getAtomMappingsReactantProduct",
                        "product row out of range",
                        Some(target_row),
                        Some(product_row),
                        None,
                    ));
                }
                mapping.product_to_reactant.insert(product_row, target_row);
            }
        } else {
            mapping.skipped[target_row] = true;
        }
    }
    Ok(mapping)
}

pub(crate) fn transfer_bond_properties(
    product: &mut ProductBuilder,
    input: ReactionInput<'_>,
    input_index: usize,
    mapping: &ReactantProductMapping,
) -> Result<(), ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: setReactantBondPropertiesToProduct
    // RDKit❗✔️: void setReactantBondPropertiesToProduct(RWMOL_SPTR product,
    // RDKit❗✔️:                                         const ROMol &reactant,
    // RDKit❗✔️:                                         ReactantProductAtomMapping *mapping) {
    // RDKit❗✔️:   for (unsigned int bidx = 0; bidx < product->getNumBonds(); ++bidx) {
    // RDKit❗✔️:     auto pBond = product->getBondWithIdx(bidx);
    // RDKit❗✔️:     auto rBondBegin = mapping->prodReactAtomMap.find(pBond->getBeginAtomIdx());
    // RDKit❗✔️:     auto rBondEnd = mapping->prodReactAtomMap.find(pBond->getEndAtomIdx());
    // RDKit❗✔️:
    // RDKit❗✔️:     if (rBondBegin == mapping->prodReactAtomMap.end() ||
    // RDKit❗✔️:         rBondEnd == mapping->prodReactAtomMap.end()) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     // the bond is between two mapped atoms from this reactant:
    // RDKit❗✔️:     const Bond *rBond =
    // RDKit❗✔️:         reactant.getBondBetweenAtoms(rBondBegin->second, rBondEnd->second);
    // RDKit❗✔️:     if (!rBond) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (!pBond->hasProp(common_properties::NullBond) &&
    // RDKit❗✔️:         !pBond->hasProp(common_properties::_MolFileBondQuery) &&
    // RDKit❗✔️:         !rBond->hasQuery()) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     if (!rBond->hasQuery()) {
    // RDKit❗✔️:       pBond->setBondType(rBond->getBondType());
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       QueryBond qBond(rBond->getBondType());
    // RDKit❗✔️:       qBond.setQuery(rBond->getQuery()->copy());
    // RDKit❗✔️:       // replaceBond copies, so we are safe passing a pointer
    // RDKit❗✔️:       // to a local:
    // RDKit❗✔️:       product->replaceBond(bidx, &qBond);
    // RDKit❗✔️:       pBond = product->getBondWithIdx(bidx);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (rBond->getBondType() == Bond::DOUBLE &&
    // RDKit❗✔️:         rBond->getBondDir() == Bond::EITHERDOUBLE) {
    // RDKit❗✔️:       pBond->setBondDir(Bond::EITHERDOUBLE);
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     pBond->setIsAromatic(rBond->getIsAromatic());
    // RDKit❗✔️:
    // RDKit❗✔️:     pBond->updateProps(*rBond);
    // RDKit❗✔️:     if (pBond->hasProp(common_properties::NullBond)) {
    // RDKit❗✔️:       pBond->clearProp(common_properties::NullBond);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    for bond in &mut product.topology.bonds {
        let (Some(&begin), Some(&end)) = (
            mapping.product_to_reactant.get(&bond.begin().index()),
            mapping.product_to_reactant.get(&bond.end().index()),
        ) else {
            continue;
        };
        let Some(neighbor) = input
            .topology
            .adjacency
            .neighbors_of(begin)
            .iter()
            .find(|n| n.atom_index == end)
        else {
            continue;
        };
        let original = &input.topology.bonds[neighbor.bond.index()];
        if original.query().is_some() {
            return Err(ReactionProductError::QueryReactantBond {
                bond: original.id(),
            });
        }
        if bond.prop("NullBond").is_none() && bond.prop("_MolFileBondQuery").is_none() {
            continue;
        }
        bond.set_order(original.order());
        if original.order() == BondOrder::Double
            && original.direction() == BondDirection::EitherDouble
        {
            bond.set_direction(BondDirection::EitherDouble);
        }
        bond.set_aromatic(original.is_aromatic());
        bond.update_properties_from(original, false);
        if bond.prop("NullBond").is_some() {
            bond.clear_prop("NullBond");
        }
        product.bond_origins[bond.id().index()] = Some(ReactionRowOrigin {
            input: input_index,
            row: original.id(),
        });
    }
    Ok(())
}

pub(crate) fn check_product_chirality(
    reactant_tag: ChiralTag,
    atom: &mut Atom,
) -> Result<(), ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: checkProductChirality
    // RDKit❗✔️: void checkProductChirality(Atom::ChiralType reactantChirality,
    // RDKit❗✔️:                            Atom *productAtom) {
    // RDKit❗✔️:   int flagVal;
    // RDKit❗✔️:   productAtom->getProp(common_properties::molInversionFlag, flagVal);
    // RDKit❗✔️:   switch (flagVal) {
    // RDKit❗✔️:     case 0:
    // RDKit❗✔️:       // reaction doesn't have anything to say about the chirality
    // RDKit❗✔️:       // FIX: should we clear the chirality or leave it alone? for now we leave
    // RDKit❗✔️:       // it alone
    // RDKit❗✔️:       productAtom->setChiralTag(reactantChirality);
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case 1:
    // RDKit❗✔️:       // reaction inverts chirality
    // RDKit❗✔️:       if (reactantChirality != Atom::CHI_TETRAHEDRAL_CW &&
    // RDKit❗✔️:           reactantChirality != Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit❗✔️:         BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:             << "unsupported chiral type on reactant atom ignored\n";
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         productAtom->setChiralTag(reactantChirality);
    // RDKit❗✔️:         productAtom->invertChirality();
    // RDKit❗✔️:       }
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case 2:
    // RDKit❗✔️:       // reaction retains chirality:
    // RDKit❗✔️:       // retention: just set to the reactant
    // RDKit❗✔️:       productAtom->setChiralTag(reactantChirality);
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case 3:
    // RDKit❗✔️:       // reaction destroys chirality:
    // RDKit❗✔️:       // remove stereo
    // RDKit❗✔️:       productAtom->setChiralTag(Atom::CHI_UNSPECIFIED);
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case 4:
    // RDKit❗✔️:       // reaction creates chirality.
    // RDKit❗✔️:       // set stereo, so leave it the way it was in the product template
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     default:
    // RDKit❗✔️:       BOOST_LOG(rdWarningLog) << "unrecognized chiral inversion/retention flag "
    // RDKit❗✔️:                                  "on product atom ignored\n";
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    let flag = inversion_flag(atom)?.ok_or(ReactionProductError::MissingProperty {
        atom: atom.id(),
        key: "molInversionFlag",
    })?;
    match flag {
        0 | 2 => atom.set_chiral_tag(reactant_tag),
        1 => match reactant_tag {
            ChiralTag::TetrahedralCw => atom.set_chiral_tag(ChiralTag::TetrahedralCcw),
            ChiralTag::TetrahedralCcw => atom.set_chiral_tag(ChiralTag::TetrahedralCw),
            _ => eprintln!("unsupported chiral type on reactant atom ignored"),
        },
        3 => atom.set_chiral_tag(ChiralTag::Unspecified),
        4 => {}
        _ => eprintln!("unrecognized chiral inversion/retention flag on product atom ignored"),
    }
    Ok(())
}

pub(crate) fn transfer_atom_properties(
    product: &mut Atom,
    reactant: &Atom,
    implicit: bool,
    input_index: usize,
) -> Result<(), ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: setReactantAtomPropertiesToProduct
    // RDKit❗✔️: void setReactantAtomPropertiesToProduct(Atom *productAtom,
    // RDKit❗✔️:                                         const Atom &reactantAtom,
    // RDKit❗✔️:                                         bool setImplicitProperties,
    // RDKit❗✔️: 					unsigned int reactantId) {
    // RDKit❗✔️:   // which properties need to be set from the reactant?
    // RDKit❗✔️:   if (productAtom->getAtomicNum() <= 0 ||
    // RDKit❗✔️:       productAtom->hasProp(common_properties::_MolFileAtomQuery)) {
    // RDKit❗✔️:     productAtom->setAtomicNum(reactantAtom.getAtomicNum());
    // RDKit❗✔️:     productAtom->setIsAromatic(reactantAtom.getIsAromatic());
    // RDKit❗✔️:     // don't copy isotope information over from dummy atoms
    // RDKit❗✔️:     // (part of github #243) unless we're setting implicit properties,
    // RDKit❗✔️:     // in which case we do need to copy them in (github #1269)
    // RDKit❗✔️:     if (!setImplicitProperties) {
    // RDKit❗✔️:       productAtom->setIsotope(reactantAtom.getIsotope());
    // RDKit❗✔️:     }
    // RDKit❗✔️:     // remove dummy labels (if present)
    // RDKit❗✔️:     if (productAtom->hasProp(common_properties::dummyLabel)) {
    // RDKit❗✔️:       productAtom->clearProp(common_properties::dummyLabel);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (productAtom->hasProp(common_properties::_MolFileRLabel)) {
    // RDKit❗✔️:       productAtom->clearProp(common_properties::_MolFileRLabel);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     productAtom->setProp(WAS_DUMMY, true);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     // remove bookkeeping labels (if present)
    // RDKit❗✔️:     if (productAtom->hasProp(WAS_DUMMY)) {
    // RDKit❗✔️:       productAtom->clearProp(WAS_DUMMY);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   productAtom->setProp<unsigned int>(common_properties::reactantAtomIdx,
    // RDKit❗✔️:                                      reactantAtom.getIdx());
    // RDKit❗✔️:   productAtom->setProp<unsigned int>(common_properties::reactantIdx,
    // RDKit❗✔️: 				     reactantId);
    // RDKit❗✔️:   if (setImplicitProperties) {
    // RDKit❗✔️:     updateImplicitAtomProperties(productAtom, &reactantAtom);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // One might be tempted to copy over the reactant atom's chirality into the
    // RDKit❗✔️:   // product atom if chirality is not specified on the product. This would be a
    // RDKit❗✔️:   // very bad idea because the order of bonds will almost certainly change on
    // RDKit❗✔️:   // the atom and the chirality is referenced to bond order.
    // RDKit❗✔️:
    // RDKit❗✔️:   // --------- --------- --------- --------- --------- ---------
    // RDKit❗✔️:   // While we're here, set the stereochemistry
    // RDKit❗✔️:   // FIX: this should be free-standing, not in this function.
    // RDKit❗✔️:   if (reactantAtom.getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit❗✔️:       reactantAtom.getChiralTag() != Atom::CHI_OTHER &&
    // RDKit❗✔️:       productAtom->hasProp(common_properties::molInversionFlag)) {
    // RDKit❗✔️:     checkProductChirality(reactantAtom.getChiralTag(), productAtom);
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // copy over residue information if it's there. This was github #1632
    // RDKit❗✔️:   if (reactantAtom.getMonomerInfo()) {
    // RDKit❗✔️:     productAtom->setMonomerInfo(reactantAtom.getMonomerInfo()->copy());
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    if product.atomic_number() == 0 || product.prop("_MolFileAtomQuery").is_some() {
        product.set_element(reactant.element());
        product.set_aromatic(reactant.is_aromatic());
        if !implicit {
            product.set_isotope(reactant.isotope());
        }
        product.clear_prop("dummyLabel");
        product.clear_prop("_MolFileRLabel");
        product.set_prop("was_dummy", true)?;
    } else {
        product.clear_prop("was_dummy");
    }
    product.set_prop(
        "react_atom_idx",
        source_u32("reactant atom", reactant.id().index())?,
    )?;
    product.set_prop("react_idx", source_u32("reactant", input_index)?)?;
    if implicit {
        update_implicit_atom_properties(product, reactant);
    }
    if !matches!(
        reactant.chiral_tag(),
        ChiralTag::Unspecified | ChiralTag::Other
    ) && inversion_flag(product)?.is_some()
    {
        check_product_chirality(reactant.chiral_tag(), product)?;
    }
    if let Some(info) = reactant.pdb_residue_info() {
        product.set_pdb_residue_info(Some(info.clone()));
    }
    Ok(())
}

fn add_source_bond(
    original: &Bond,
    product: &mut ProductBuilder,
    begin: usize,
    end: usize,
    input_index: usize,
) -> Result<BondId, ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: addBondToProduct
    // RDKit❗✔️: Bond *addBondToProduct(const Bond &origB, RWMol &product,
    // RDKit❗✔️:                        unsigned int begAtomIdx, unsigned int endAtomIdx) {
    // RDKit❗✔️:   if (!origB.hasQuery()) {
    // RDKit❗✔️:     auto idx = product.addBond(begAtomIdx, endAtomIdx, origB.getBondType());
    // RDKit❗✔️:     return product.getBondWithIdx(idx - 1);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     QueryBond *qbond = new QueryBond(origB.getBondType());
    // RDKit❗✔️:     qbond->setBeginAtomIdx(begAtomIdx);
    // RDKit❗✔️:     qbond->setEndAtomIdx(endAtomIdx);
    // RDKit❗✔️:     qbond->setQuery(origB.getQuery()->copy());
    // RDKit❗✔️:     bool takeOwnership = true;
    // RDKit❗✔️:     product.addBond(qbond, takeOwnership);
    // RDKit❗✔️:     return qbond;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    if original.query().is_some() {
        return Err(ReactionProductError::QueryReactantBond {
            bond: original.id(),
        });
    }
    product.add_bond(
        begin,
        end,
        original.order(),
        Some(ReactionRowOrigin {
            input: input_index,
            row: original.id(),
        }),
    )
}

fn add_missing_bonds(
    original: &Bond,
    product: &mut ProductBuilder,
    mapping: &ReactantProductMapping,
    input_index: usize,
) -> Result<(), ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: addMissingProductBonds
    // RDKit❗✔️: void addMissingProductBonds(const Bond &origB, RWMOL_SPTR product,
    // RDKit❗✔️:                             ReactantProductAtomMapping *mapping) {
    // RDKit❗✔️:   unsigned int begIdx = origB.getBeginAtomIdx();
    // RDKit❗✔️:   unsigned int endIdx = origB.getEndAtomIdx();
    // RDKit❗✔️:
    // RDKit❗✔️:   std::vector<unsigned> prodBeginIdxs = mapping->reactProdAtomMap[begIdx];
    // RDKit❗✔️:   std::vector<unsigned> prodEndIdxs = mapping->reactProdAtomMap[endIdx];
    // RDKit❗✔️:   CHECK_INVARIANT(prodBeginIdxs.size() == prodEndIdxs.size(),
    // RDKit❗✔️:                   "Different number of start-end points for product bonds.");
    // RDKit❗✔️:   for (unsigned i = 0; i < prodBeginIdxs.size(); i++) {
    // RDKit❗✔️:     addBondToProduct(origB, *product, prodBeginIdxs.at(i), prodEndIdxs.at(i));
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    let begins = mapping
        .reactant_to_product
        .get(&original.begin().index())
        .map_or(&[][..], Vec::as_slice);
    let ends = mapping
        .reactant_to_product
        .get(&original.end().index())
        .map_or(&[][..], Vec::as_slice);
    if begins.len() != ends.len() {
        return Err(invariant(
            "addMissingProductBonds",
            "different number of start-end points",
            Some(original.begin().index()),
            None,
            Some(original.id()),
        ));
    }
    for (&begin, &end) in begins.iter().zip(ends) {
        add_source_bond(original, product, begin, end, input_index)?;
    }
    Ok(())
}

fn add_missing_atom(
    reactant_atom: &Atom,
    reactant_neighbor: usize,
    product_neighbor: usize,
    product: &mut ProductBuilder,
    input: ReactionInput<'_>,
    mapping: &mut ReactantProductMapping,
    input_index: usize,
) -> Result<(), ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: addMissingProductAtom
    // RDKit❗✔️: void addMissingProductAtom(const Atom &reactAtom, unsigned reactNeighborIdx,
    // RDKit❗✔️:                            unsigned prodNeighborIdx, RWMOL_SPTR product,
    // RDKit❗✔️:                            const ROMol &reactant,
    // RDKit❗✔️:                            ReactantProductAtomMapping *mapping, unsigned int reactantId) {
    // RDKit❗✔️:   Atom *newAtom = nullptr;
    // RDKit❗✔️:   if (!reactAtom.hasQuery()) {
    // RDKit❗✔️:     newAtom = new Atom(reactAtom);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     newAtom = new QueryAtom(dynamic_cast<const QueryAtom &>(reactAtom));
    // RDKit❗✔️:   }
    // RDKit❗✔️:   unsigned reactAtomIdx = reactAtom.getIdx();
    // RDKit❗✔️:   newAtom->setProp<unsigned int>(common_properties::reactantAtomIdx,
    // RDKit❗✔️:                                  reactAtomIdx);
    // RDKit❗✔️:   newAtom->setProp<unsigned int>(common_properties::reactantIdx,
    // RDKit❗✔️:                                  reactantId);
    // RDKit❗✔️:   unsigned productIdx = product->addAtom(newAtom, false, true);
    // RDKit❗✔️:   mapping->reactProdAtomMap[reactAtomIdx].push_back(productIdx);
    // RDKit❗✔️:   mapping->prodReactAtomMap[productIdx] = reactAtomIdx;
    // RDKit❗✔️:   // add the bonds
    // RDKit❗✔️:   const Bond *origB =
    // RDKit❗✔️:       reactant.getBondBetweenAtoms(reactNeighborIdx, reactAtomIdx);
    // RDKit❗✔️:   unsigned int begIdx = productIdx;
    // RDKit❗✔️:   unsigned int endIdx = prodNeighborIdx;
    // RDKit❗✔️:   if (origB->getBeginAtomIdx() == reactNeighborIdx) {
    // RDKit❗✔️:     std::swap(begIdx, endIdx);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   Bond *prodB = addBondToProduct(*origB, *product, begIdx, endIdx);
    // RDKit❗✔️:   if (origB->getBondType() == Bond::DOUBLE &&
    // RDKit❗✔️:       origB->getBondDir() == Bond::EITHERDOUBLE) {
    // RDKit❗✔️:     prodB->setBondDir(Bond::EITHERDOUBLE);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   bool preserveExisting = true;
    // RDKit❗✔️:   prodB->updateProps(*origB, preserveExisting);
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    let mut atom = reactant_atom.clone();
    let target = atom.id().index();
    atom.set_prop("react_atom_idx", source_u32("reactant atom", target)?)?;
    atom.set_prop("react_idx", source_u32("reactant", input_index)?)?;
    let row = product.add_atom(
        atom,
        Some(ReactionRowOrigin {
            input: input_index,
            row: reactant_atom.id(),
        }),
    );
    mapping
        .reactant_to_product
        .entry(target)
        .or_default()
        .push(row);
    mapping.product_to_reactant.insert(row, target);
    let neighbor = input
        .topology
        .adjacency
        .neighbors_of(reactant_neighbor)
        .iter()
        .find(|n| n.atom_index == target)
        .ok_or_else(|| {
            invariant(
                "addMissingProductAtom",
                "missing reactant neighbor bond",
                Some(target),
                Some(row),
                None,
            )
        })?;
    let original = &input.topology.bonds[neighbor.bond.index()];
    let (begin, end) = if original.begin().index() == reactant_neighbor {
        (product_neighbor, row)
    } else {
        (row, product_neighbor)
    };
    let id = add_source_bond(original, product, begin, end, input_index)?;
    let bond = &mut product.topology.bonds[id.index()];
    if original.order() == BondOrder::Double && original.direction() == BondDirection::EitherDouble
    {
        bond.set_direction(BondDirection::EitherDouble);
    }
    bond.update_properties_from(original, true);
    Ok(())
}

pub(crate) fn add_neighbors(
    input: ReactionInput<'_>,
    start: usize,
    product: &mut ProductBuilder,
    visited: &mut [bool],
    chiral_to_check: &mut Vec<usize>,
    mapping: &mut ReactantProductMapping,
    input_index: usize,
) -> Result<(), ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: addReactantNeighborsToProduct
    // RDKit❗✔️: void addReactantNeighborsToProduct(
    // RDKit❗✔️:     const ROMol &reactant, const Atom &reactantAtom, RWMOL_SPTR product,
    // RDKit❗✔️:     boost::dynamic_bitset<> &visitedAtoms,
    // RDKit❗✔️:     std::vector<const Atom *> &chiralAtomsToCheck,
    // RDKit❗✔️:     ReactantProductAtomMapping *mapping, unsigned int reactantId) {
    // RDKit❗✔️:   std::list<const Atom *> atomStack;
    // RDKit❗✔️:   atomStack.push_back(&reactantAtom);
    // RDKit❗✔️:
    // RDKit❗✔️:   // std::cerr << "-------------------" << std::endl;
    // RDKit❗✔️:   // std::cerr << "  add reactant neighbors from: " << reactantAtom.getIdx()
    // RDKit❗✔️:   //           << std::endl;
    // RDKit❗✔️:   // #if 1
    // RDKit❗✔️:   //   product->updatePropertyCache(false);
    // RDKit❗✔️:   //   product->debugMol(std::cerr);
    // RDKit❗✔️:   //   std::cerr << "-------------------" << std::endl;
    // RDKit❗✔️:   // #endif
    // RDKit❗✔️:
    // RDKit❗✔️:   while (!atomStack.empty()) {
    // RDKit❗✔️:     const Atom *lReactantAtom = atomStack.front();
    // RDKit❗✔️:     // std::cerr << "    front: " << lReactantAtom->getIdx() << std::endl;
    // RDKit❗✔️:     atomStack.pop_front();
    // RDKit❗✔️:
    // RDKit❗✔️:     // each atom in the stack is guaranteed to already be in the product:
    // RDKit❗✔️:     CHECK_INVARIANT(mapping->reactProdAtomMap.find(lReactantAtom->getIdx()) !=
    // RDKit❗✔️:                         mapping->reactProdAtomMap.end(),
    // RDKit❗✔️:                     "reactant atom on traversal stack not present in product.");
    // RDKit❗✔️:
    // RDKit❗✔️:     std::vector<unsigned> lReactantAtomProductIndex =
    // RDKit❗✔️:         mapping->reactProdAtomMap[lReactantAtom->getIdx()];
    // RDKit❗✔️:     unsigned lreactIdx = lReactantAtom->getIdx();
    // RDKit❗✔️:     visitedAtoms[lreactIdx] = 1;
    // RDKit❗✔️:     // Check our neighbors:
    // RDKit❗✔️:     ROMol::ADJ_ITER nbrIdx, endNbrs;
    // RDKit❗✔️:     boost::tie(nbrIdx, endNbrs) = reactant.getAtomNeighbors(lReactantAtom);
    // RDKit❗✔️:     while (nbrIdx != endNbrs) {
    // RDKit❗✔️:       // Four possibilities here. The neighbor:
    // RDKit❗✔️:       //  0) has been visited already: do nothing
    // RDKit❗✔️:       //  1) is part of the match (thus already in the product): set a bond to
    // RDKit❗✔️:       //  it
    // RDKit❗✔️:       //  2) has been added: set a bond to it
    // RDKit❗✔️:       //  3) has not yet been added: add it, set a bond to it, and push it
    // RDKit❗✔️:       //     onto the stack
    // RDKit❗✔️:       // std::cerr << "       nbr: " << *nbrIdx << std::endl;
    // RDKit❗✔️:       // std::cerr << "              visited: " << visitedAtoms[*nbrIdx]
    // RDKit❗✔️:       //           << "  skipped: " << mapping->skippedAtoms[*nbrIdx]
    // RDKit❗✔️:       //           << " mapped: " << mapping->mappedAtoms[*nbrIdx]
    // RDKit❗✔️:       //           << " mappedO: " << mapping->mappedAtoms[lreactIdx] <<
    // RDKit❗✔️:       //           std::endl;
    // RDKit❗✔️:       if (!visitedAtoms[*nbrIdx] && !mapping->skippedAtoms[*nbrIdx]) {
    // RDKit❗✔️:         if (mapping->mappedAtoms[*nbrIdx]) {
    // RDKit❗✔️:           // this is case 1 (neighbor in match); set a bond to the neighbor if
    // RDKit❗✔️:           // this atom
    // RDKit❗✔️:           // is not also in the match (match-match bonds were set when the
    // RDKit❗✔️:           // product template was
    // RDKit❗✔️:           // copied in to start things off).;
    // RDKit❗✔️:           if (!mapping->mappedAtoms[lreactIdx]) {
    // RDKit❗✔️:             CHECK_INVARIANT(mapping->reactProdAtomMap.find(*nbrIdx) !=
    // RDKit❗✔️:                                 mapping->reactProdAtomMap.end(),
    // RDKit❗✔️:                             "reactant atom not present in product.");
    // RDKit❗✔️:             const Bond *origB =
    // RDKit❗✔️:                 reactant.getBondBetweenAtoms(lreactIdx, *nbrIdx);
    // RDKit❗✔️:             addMissingProductBonds(*origB, product, mapping);
    // RDKit❗✔️:           } else {
    // RDKit❗✔️:             // both mapped atoms are in the match.
    // RDKit❗✔️:             // they are bonded in the reactant (otherwise we wouldn't be here),
    // RDKit❗✔️:             //
    // RDKit❗✔️:             // If they do not have already have a bond in the product and did
    // RDKit❗✔️:             // not have one in the reactant template then set one here
    // RDKit❗✔️:             // If they do have a bond in the reactant template, then we
    // RDKit❗✔️:             // assume that this is an intentional bond break, so we don't do
    // RDKit❗✔️:             // anything
    // RDKit❗✔️:             //
    // RDKit❗✔️:             // this was github #1387
    // RDKit❗✔️:             unsigned prodBeginIdx = mapping->reactProdAtomMap[lreactIdx][0];
    // RDKit❗✔️:             unsigned prodEndIdx = mapping->reactProdAtomMap[*nbrIdx][0];
    // RDKit❗✔️:             if (!product->getBondBetweenAtoms(prodBeginIdx, prodEndIdx)) {
    // RDKit❗✔️:               // They must be mapped
    // RDKit❗✔️:               CHECK_INVARIANT(
    // RDKit❗✔️:                   product->getAtomWithIdx(prodBeginIdx)
    // RDKit❗✔️:                           ->hasProp(common_properties::reactionMapNum) &&
    // RDKit❗✔️:                       product->getAtomWithIdx(prodEndIdx)
    // RDKit❗✔️:                           ->hasProp(common_properties::reactionMapNum),
    // RDKit❗✔️:                   "atoms should be mapped in product");
    // RDKit❗✔️:               int a1mapidx =
    // RDKit❗✔️:                   product->getAtomWithIdx(prodBeginIdx)
    // RDKit❗✔️:                       ->getProp<int>(common_properties::reactionMapNum);
    // RDKit❗✔️:               int a2mapidx =
    // RDKit❗✔️:                   product->getAtomWithIdx(prodEndIdx)
    // RDKit❗✔️:                       ->getProp<int>(common_properties::reactionMapNum);
    // RDKit❗✔️:               if (a1mapidx > a2mapidx) {
    // RDKit❗✔️:                 std::swap(a1mapidx, a2mapidx);
    // RDKit❗✔️:               }
    // RDKit❗✔️:               if (mapping->reactantTemplateAtomBonds.find(
    // RDKit❗✔️:                       std::make_pair(a1mapidx, a2mapidx)) ==
    // RDKit❗✔️:                   mapping->reactantTemplateAtomBonds.end()) {
    // RDKit❗✔️:                 const Bond *origB =
    // RDKit❗✔️:                     reactant.getBondBetweenAtoms(lreactIdx, *nbrIdx);
    // RDKit❗✔️:                 addMissingProductBonds(*origB, product, mapping);
    // RDKit❗✔️:               }
    // RDKit❗✔️:             }
    // RDKit❗✔️:           }
    // RDKit❗✔️:         } else if (mapping->reactProdAtomMap.find(*nbrIdx) !=
    // RDKit❗✔️:                    mapping->reactProdAtomMap.end()) {
    // RDKit❗✔️:           // case 2, the neighbor has been added and we just need to set a bond
    // RDKit❗✔️:           // to it:
    // RDKit❗✔️:           const Bond *origB = reactant.getBondBetweenAtoms(lreactIdx, *nbrIdx);
    // RDKit❗✔️:           addMissingProductBonds(*origB, product, mapping);
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           // case 3, add the atom, a bond to it, and push the atom onto the
    // RDKit❗✔️:           // stack
    // RDKit❗✔️:           const Atom *neighbor = reactant.getAtomWithIdx(*nbrIdx);
    // RDKit❗✔️:           for (unsigned int i : lReactantAtomProductIndex) {
    // RDKit❗✔️:             addMissingProductAtom(*neighbor, lreactIdx, i, product, reactant,
    // RDKit❗✔️:                                   mapping, reactantId);
    // RDKit❗✔️:           }
    // RDKit❗✔️:           // update the stack:
    // RDKit❗✔️:           atomStack.push_back(neighbor);
    // RDKit❗✔️:           // if the atom is chiral, we need to check its bond ordering later:
    // RDKit❗✔️:           if (neighbor->getChiralTag() != Atom::CHI_UNSPECIFIED) {
    // RDKit❗✔️:             chiralAtomsToCheck.push_back(neighbor);
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:       nbrIdx++;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }  // end of atomStack traversal
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    let mut stack = VecDeque::from([start]);
    while let Some(current) = stack.pop_front() {
        let current_products = mapping
            .reactant_to_product
            .get(&current)
            .ok_or_else(|| {
                invariant(
                    "addReactantNeighborsToProduct",
                    "traversal atom not present in product",
                    Some(current),
                    None,
                    None,
                )
            })?
            .clone();
        visited[current] = true;
        for neighbor in input.topology.adjacency.neighbors_of(current) {
            let target = neighbor.atom_index;
            if visited[target] || mapping.skipped[target] {
                continue;
            }
            let original = &input.topology.bonds[neighbor.bond.index()];
            if mapping.mapped[target] {
                if !mapping.mapped[current] {
                    if !mapping.reactant_to_product.contains_key(&target) {
                        return Err(invariant(
                            "addReactantNeighborsToProduct",
                            "reactant atom not present in product",
                            Some(target),
                            None,
                            Some(original.id()),
                        ));
                    }
                    add_missing_bonds(original, product, mapping, input_index)?;
                } else {
                    let begin = mapping
                        .reactant_to_product
                        .get(&current)
                        .and_then(|rows| rows.first())
                        .copied()
                        .ok_or_else(|| {
                            invariant(
                                "addReactantNeighborsToProduct",
                                "missing mapped begin row",
                                Some(current),
                                None,
                                Some(original.id()),
                            )
                        })?;
                    let end = mapping
                        .reactant_to_product
                        .get(&target)
                        .and_then(|rows| rows.first())
                        .copied()
                        .ok_or_else(|| {
                            invariant(
                                "addReactantNeighborsToProduct",
                                "missing mapped end row",
                                Some(target),
                                None,
                                Some(original.id()),
                            )
                        })?;
                    if product.bond_between(begin, end).is_none() {
                        let a = &product.topology.atoms[begin];
                        let b = &product.topology.atoms[end];
                        let mut maps = [
                            int_prop(a, "old_mapno")?.ok_or(
                                ReactionProductError::MissingProperty {
                                    atom: a.id(),
                                    key: "old_mapno",
                                },
                            )?,
                            int_prop(b, "old_mapno")?.ok_or(
                                ReactionProductError::MissingProperty {
                                    atom: b.id(),
                                    key: "old_mapno",
                                },
                            )?,
                        ];
                        maps.sort();
                        if !mapping
                            .template_bonds
                            .contains(&(maps[0] as u32, maps[1] as u32))
                        {
                            add_missing_bonds(original, product, mapping, input_index)?;
                        }
                    }
                }
            } else if mapping.reactant_to_product.contains_key(&target) {
                add_missing_bonds(original, product, mapping, input_index)?;
            } else {
                let atom = &input.topology.atoms[target];
                for &product_neighbor in &current_products {
                    add_missing_atom(
                        atom,
                        current,
                        product_neighbor,
                        product,
                        input,
                        mapping,
                        input_index,
                    )?;
                }
                stack.push_back(target);
                if atom.chiral_tag() != ChiralTag::Unspecified {
                    chiral_to_check.push(target);
                }
            }
        }
    }
    Ok(())
}
