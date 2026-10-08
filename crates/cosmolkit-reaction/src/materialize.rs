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
        let count = cosmolkit_model::add_source_bond_order(
            &mut self.topology.atoms,
            &mut self.topology.bonds,
            cosmolkit_model::SourceBondNeighbors {
                original: None,
                appended: &mut self.neighbors,
            },
            None,
            AtomId::new(begin),
            AtomId::new(end),
            order,
        )?;
        // Product provenance belongs to this consumer, after native append.
        let id = BondId::new(count - 1);
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

fn new_reactant_product_mapping(length: u32) -> ReactantProductMapping {
    // RDKit❗❌:   ReactantProductAtomMapping(unsigned lenghtBitSet) {
    // RDKit❗❌:     mappedAtoms.resize(lenghtBitSet);
    // RDKit❗❌:     skippedAtoms.resize(lenghtBitSet);
    // RDKit❗❌:   }
    // Boost❗❌: template <typename Block, typename Allocator>
    // Boost❗❌: void dynamic_bitset<Block, Allocator>::
    // Boost❗❌: resize(size_type num_bits, bool value) // strong guarantee
    // Boost❗❌: {
    // Boost❗❌:
    // Boost❗❌:   const size_type old_num_blocks = num_blocks();
    // Boost❗❌:   const size_type required_blocks = calc_num_blocks(num_bits);
    // Boost❗❌:
    // Boost❗❌:   const block_type v = value? detail::dynamic_bitset_impl::max_limit<Block>::value : Block(0);
    // Boost❗❌:
    // Boost❗❌:   if (required_blocks != old_num_blocks) {
    // Boost❗❌:     m_bits.resize(required_blocks, v); // s.g. (copy)
    // Boost❗❌:   }
    // Boost❗❌:
    // Boost❗❌:
    // Boost❗❌:   // At this point:
    // Boost❗❌:   //
    // Boost❗❌:   //  - if the buffer was shrunk, we have nothing more to do,
    // Boost❗❌:   //    except a call to m_zero_unused_bits()
    // Boost❗❌:   //
    // Boost❗❌:   //  - if it was enlarged, all the (used) bits in the new blocks have
    // Boost❗❌:   //    the correct value, but we have not yet touched those bits, if
    // Boost❗❌:   //    any, that were 'unused bits' before enlarging: if value == true,
    // Boost❗❌:   //    they must be set.
    // Boost❗❌:
    // Boost❗❌:   if (value && (num_bits > m_num_bits)) {
    // Boost❗❌:
    // Boost❗❌:     const block_width_type extra_bits = count_extra_bits();
    // Boost❗❌:     if (extra_bits) {
    // Boost❗❌:         assert(old_num_blocks >= 1 && old_num_blocks <= m_bits.size());
    // Boost❗❌:
    // Boost❗❌:         // Set them.
    // Boost❗❌:         m_bits[old_num_blocks - 1] |= (v << extra_bits);
    // Boost❗❌:     }
    // Boost❗❌:
    // Boost❗❌:   }
    // Boost❗❌:
    // Boost❗❌:   m_num_bits = num_bits;
    // Boost❗❌:   m_zero_unused_bits();
    // Boost❗❌:
    // Boost❗❌: }
    // Only constructor-reached empty-bitset resize(false) is projected here;
    // all new bits are false and the default-constructed maps are empty.
    // Vec<bool> stores bytes rather than packed Boost blocks: materially more
    // storage/initialization traffic, explicitly cost ❌ until final comparison.
    ReactantProductMapping {
        mapped: vec![false; length as usize],
        skipped: vec![false; length as usize],
        reactant_to_product: BTreeMap::new(),
        product_to_reactant: BTreeMap::new(),
        product_atom_bond: BTreeMap::new(),
        template_bonds: BTreeSet::new(),
    }
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
    fn typed_map(&self) -> Option<u32>;
    fn raw_property_exists(&self, key: &[u8]) -> bool;
    fn typed_permutation(&self) -> Option<u32>;
    fn typed_parity(&self) -> Option<i32>;
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
    fn typed_map(&self) -> Option<u32> {
        self.atom_map()
    }
    fn raw_property_exists(&self, key: &[u8]) -> bool {
        self.prop(key).is_some()
    }
    fn typed_permutation(&self) -> Option<u32> {
        self.chiral_permutation()
    }
    fn typed_parity(&self) -> Option<i32> {
        self.mol_parity()
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
    fn typed_map(&self) -> Option<u32> {
        self.atom_map()
    }
    fn raw_property_exists(&self, key: &[u8]) -> bool {
        self.prop(key).is_some()
    }
    fn typed_permutation(&self) -> Option<u32> {
        self.chiral_permutation()
    }
    fn typed_parity(&self) -> Option<i32> {
        self.mol_parity()
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

fn update_from_template_source<T: ReactionAtomPropertyRead>(
    template: Option<&T>,
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
    // None explicitly represents the native same-object template/output alias;
    // each property read completes before the reached scalar write, with no
    // cloned atom, delayed mutation or property snapshot.
    let mut changed = false;
    if let Some(value) = (match template {
        Some(template) => int_prop(template, "_QueryFormalCharge"),
        None => int_prop(atom, "_QueryFormalCharge"),
    })? && value != i32::from(atom.formal_charge())
    {
        atom.set_formal_charge(value as i8);
        changed = true;
    }
    if let Some(value) = (match template {
        Some(template) => uint_prop(template, "_QueryHCount"),
        None => uint_prop(atom, "_QueryHCount"),
    })? && (!atom.no_implicit() || u32::from(atom.explicit_hydrogens()) != value)
    {
        atom.set_explicit_hydrogens(value as u8);
        atom.set_no_implicit(true);
        changed = true;
    }
    // Source reads/converts this value even though it does not set mass.
    let _mass = (match template {
        Some(template) => uint_prop(template, "_QueryMass"),
        None => uint_prop(atom, "_QueryMass"),
    })?;
    if let Some(value) = (match template {
        Some(template) => uint_prop(template, "_QueryIsotope"),
        None => uint_prop(atom, "_QueryIsotope"),
    })? && u32::from(atom.isotope().unwrap_or(0)) != value
    {
        atom.set_isotope(Some(value as u16));
        changed = true;
    }
    Ok(changed)
}

pub(crate) fn update_from_template(
    template: &impl ReactionAtomPropertyRead,
    atom: &mut Atom,
) -> Result<bool, ReactionProductError> {
    update_from_template_source(Some(template), atom)
}

fn template_neighboring_directed_bond_source(
    template: &QueryGraph,
    atom: AtomId,
) -> Result<bool, ReactionProductError> {
    let neighbors = template.adjacency().get(atom.index()).ok_or_else(|| {
        invariant(
            "getNeighboringDirectedBond",
            "atom row out of range",
            None,
            Some(atom.index()),
            None,
        )
    })?;
    let bonds = neighbors.iter().map(|&(_, row)| {
        template
            .bonds()
            .get(row)
            .map(|bond| bond.bond())
            .ok_or_else(|| {
                invariant(
                    "getNeighboringDirectedBond",
                    "incident bond row out of range",
                    None,
                    Some(atom.index()),
                    Some(BondId::new(row)),
                )
            })
    });
    // Reuse the sole CORE source traversal with borrowed original incident
    // order; no whole QueryGraph-to-TopologyBlock copy or validation pass.
    Ok(cosmolkit_core::neighboring_directed_bond_from_incident(bonds)?.is_some())
}

pub(crate) fn convert_template(
    template: &QueryGraph,
    template_index: usize,
) -> Result<ProductBuilder, ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: convertTemplateToMol
    // RDKit❗❌: RWMOL_SPTR convertTemplateToMol(const ROMOL_SPTR prodTemplateSptr) {
    // RDKit❗❌:   const ROMol *prodTemplate = prodTemplateSptr.get();
    // RDKit❗❌:   auto *res = new RWMol();
    // RDKit❗❌:
    // RDKit❗❌:   // --------- --------- --------- --------- --------- ---------
    // RDKit❗❌:   // Initialize by making a copy of the product template as a normal molecule.
    // RDKit❗❌:   // NOTE that we can't just use a normal copy because we do not want to end up
    // RDKit❗❌:   // with query atoms or bonds in the product.
    // RDKit❗❌:
    // RDKit❗❌:   // copy in the atoms:
    // RDKit❗❌:   ROMol::ATOM_ITER_PAIR atItP = prodTemplate->getVertices();
    // RDKit❗❌:   while (atItP.first != atItP.second) {
    // RDKit❗❌:     const Atom *oAtom = (*prodTemplate)[*(atItP.first++)];
    // RDKit❗❌:     auto *newAtom = new Atom(*oAtom);
    // RDKit❗❌:     res->addAtom(newAtom, false, true);
    // RDKit❗❌:     int mapNum;
    // RDKit❗❌:     if (newAtom->getPropIfPresent(common_properties::molAtomMapNumber,
    // RDKit❗❌:                                   mapNum)) {
    // RDKit❗❌:       // set bookmarks for the mapped atoms:
    // RDKit❗❌:       res->setAtomBookmark(newAtom, mapNum);
    // RDKit❗❌:       // now clear the molAtomMapNumber property so that it doesn't
    // RDKit❗❌:       // end up in the products (this was bug 3140490):
    // RDKit❗❌:       newAtom->clearProp(common_properties::molAtomMapNumber);
    // RDKit❗❌:       newAtom->setProp<int>(common_properties::reactionMapNum, mapNum);
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     newAtom->setChiralTag(Atom::CHI_UNSPECIFIED);
    // RDKit❗❌:     // if the product-template atom has the inversion flag set
    // RDKit❗❌:     // to 4 (=SET), then bring its stereochem over, otherwise we'll
    // RDKit❗❌:     // ignore it:
    // RDKit❗❌:     int iFlag;
    // RDKit❗❌:     if (oAtom->getPropIfPresent(common_properties::molInversionFlag, iFlag)) {
    // RDKit❗❌:       if (iFlag == 4) {
    // RDKit❗❌:         newAtom->setChiralTag(oAtom->getChiralTag());
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     // check for properties we need to set:
    // RDKit❗❌:     updatePropsFromImplicitProps(newAtom, newAtom);
    // RDKit❗❌:   }
    // RDKit❗❌:   // and the bonds:
    // RDKit❗❌:   ROMol::BOND_ITER_PAIR bondItP = prodTemplate->getEdges();
    // RDKit❗❌:   while (bondItP.first != bondItP.second) {
    // RDKit❗❌:     const Bond *oldB = (*prodTemplate)[*(bondItP.first++)];
    // RDKit❗❌:     unsigned int bondIdx;
    // RDKit❗❌:     bondIdx = res->addBond(oldB->getBeginAtomIdx(), oldB->getEndAtomIdx(),
    // RDKit❗❌:                            oldB->getBondType()) -
    // RDKit❗❌:               1;
    // RDKit❗❌:     // make sure we don't lose the bond dir information:
    // RDKit❗❌:     Bond *newB = res->getBondWithIdx(bondIdx);
    // RDKit❗❌:     newB->setBondDir(oldB->getBondDir());
    // RDKit❗❌:     // Special case/hack:
    // RDKit❗❌:     //  The product has been processed by the SMARTS parser.
    // RDKit❗❌:     //  The SMARTS parser tags unspecified bonds as single, but then adds
    // RDKit❗❌:     //  a query so that they match single or double
    // RDKit❗❌:     //  This caused Issue 1748846
    // RDKit❗❌:     //   http://sourceforge.net/tracker/index.php?func=detail&aid=1748846&group_id=160139&atid=814650
    // RDKit❗❌:     //  We need to fix that little problem now:
    // RDKit❗❌:     if (oldB->hasQuery()) {
    // RDKit❗❌:       //  remember that the product has been processed by the SMARTS parser.
    // RDKit❗❌:       std::string queryDescription = oldB->getQuery()->getDescription();
    // RDKit❗❌:       if (queryDescription == "BondOr" && oldB->getBondType() == Bond::SINGLE) {
    // RDKit❗❌:         //  We need to fix that little problem now:
    // RDKit❗❌:         if (newB->getBeginAtom()->getIsAromatic() &&
    // RDKit❗❌:             newB->getEndAtom()->getIsAromatic()) {
    // RDKit❗❌:           newB->setBondType(Bond::AROMATIC);
    // RDKit❗❌:           newB->setIsAromatic(true);
    // RDKit❗❌:         } else {
    // RDKit❗❌:           newB->setBondType(Bond::SINGLE);
    // RDKit❗❌:           newB->setIsAromatic(false);
    // RDKit❗❌:         }
    // RDKit❗❌:       } else if (queryDescription == "BondNull") {
    // RDKit❗❌:         newB->setProp(common_properties::NullBond, 1);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     // Double bond stereo: if a double bond has at least one bond on each side,
    // RDKit❗❌:     // and none of those has a direction, then mark it as unknown stereo to have
    // RDKit❗❌:     // it reset later on. This has to be done before the reactant atoms are
    // RDKit❗❌:     // added,
    // RDKit❗❌:     if (oldB->getBondType() == Bond::BondType::DOUBLE) {
    // RDKit❗❌:       const Atom *startAtom = oldB->getBeginAtom();
    // RDKit❗❌:       const Atom *endAtom = oldB->getEndAtom();
    // RDKit❗❌:
    // RDKit❗❌:       if (startAtom->getDegree() > 1 && endAtom->getDegree() > 1 &&
    // RDKit❗❌:           (Chirality::getNeighboringDirectedBond(*prodTemplate, startAtom) ==
    // RDKit❗❌:                nullptr ||
    // RDKit❗❌:            Chirality::getNeighboringDirectedBond(*prodTemplate, endAtom) ==
    // RDKit❗❌:                nullptr)) {
    // RDKit❗❌:         newB->setProp(_UnknownStereoRxnBond, 1);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     // copy properties over:
    // RDKit❗❌:     bool preserveExisting = true;
    // RDKit❗❌:     newB->updateProps(*static_cast<const RDProps *>(oldB), preserveExisting);
    // RDKit❗❌:   }
    // RDKit❗❌:   return RWMOL_SPTR(res);
    // RDKit❗❌: }
    // END RDKIT COMPLETE CPP FUNCTION
    // Product source metadata starts empty; source template coordinates and
    // molecule properties are not implicitly copied into a new product.
    // Extra detached provenance vectors and copied dictionary key/order storage
    // add material allocation costs versus native conversion; cost stays ❌.
    let mut result = ProductBuilder {
        topology: TopologyBlock::default(),
        neighbors: Vec::with_capacity(template.num_atoms()),
        bookmarks: BTreeMap::new(),
        atom_origins: Vec::with_capacity(template.num_atoms()),
        bond_origins: Vec::with_capacity(template.num_bonds()),
    };
    for query_atom in template.atoms() {
        let row = result.add_atom(query_atom.try_to_atom()?, None);
        if let Some(map) = crate::validation::atom_map(
            &result.topology.atoms[row],
            ReactionRole::Product,
            template_index,
        )? {
            result.bookmarks.entry(map).or_default().push(row);
            let atom = &mut result.topology.atoms[row];
            atom.clear_prop("molAtomMapNumber")?;
            atom.set_atom_map(None);
            atom.set_prop("old_mapno", map)?;
        }
        let atom = &mut result.topology.atoms[row];
        atom.set_chiral_tag(ChiralTag::Unspecified);
        if inversion_flag(query_atom)? == Some(4) {
            atom.set_chiral_tag(query_atom.chiral_tag());
        }
        update_from_template_source::<Atom>(None, atom)?;
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
                let degree = |row: usize| {
                    template.try_atom_degree(AtomId::new(row)).ok_or_else(|| {
                        invariant(
                            "convertTemplateToMol",
                            "degree atom row out of range",
                            None,
                            Some(row),
                            None,
                        )
                    })
                };
                if degree(a)? > 1
                    && degree(b)? > 1
                    && (!template_neighboring_directed_bond_source(template, old.begin())?
                        || !template_neighboring_directed_bond_source(template, old.end())?)
                {
                    bond.set_prop("_UnknownStereoRxnBond", 1)?;
                }
            }
            bond.update_properties_from(old, true);
        }
    }
    Ok(result)
}

pub(crate) fn atom_mappings_source(
    matched: impl IntoIterator<Item = Result<(i32, i32), ReactionProductError>>,
    template: &QueryGraph,
    template_index: usize,
    product: &ProductBuilder,
    reactant_atoms: u32,
) -> Result<ReactantProductMapping, ReactionProductError> {
    // RDKit❗❌: ReactantProductAtomMapping *getAtomMappingsReactantProduct(
    // RDKit❗❌:     const MatchVectType &match, const ROMol &reactantTemplate,
    // RDKit❗❌:     RWMOL_SPTR product, unsigned numReactAtoms) {
    // RDKit❗❌:   auto *mapping = new ReactantProductAtomMapping(numReactAtoms);
    // RDKit❗❌:
    // RDKit❗❌:   // keep track of which mapped atoms in the reactant template are bonded to
    // RDKit❗❌:   // each other.
    // RDKit❗❌:   // This is part of the fix for #1387
    // RDKit❗❌:   {
    // RDKit❗❌:     ROMol::EDGE_ITER firstB, lastB;
    // RDKit❗❌:     boost::tie(firstB, lastB) = reactantTemplate.getEdges();
    // RDKit❗❌:     while (firstB != lastB) {
    // RDKit❗❌:       const Bond *bond = reactantTemplate[*firstB];
    // RDKit❗❌:       // this will put in pairs with 0s for things that aren't mapped, but we
    // RDKit❗❌:       // don't care about that
    // RDKit❗❌:       int a1mapidx = bond->getBeginAtom()->getAtomMapNum();
    // RDKit❗❌:       int a2mapidx = bond->getEndAtom()->getAtomMapNum();
    // RDKit❗❌:       if (a1mapidx > a2mapidx) {
    // RDKit❗❌:         std::swap(a1mapidx, a2mapidx);
    // RDKit❗❌:       }
    // RDKit❗❌:       mapping->reactantTemplateAtomBonds[std::make_pair(a1mapidx, a2mapidx)] =
    // RDKit❗❌:           1;
    // RDKit❗❌:       ++firstB;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   for (const auto &i : match) {
    // RDKit❗❌:     const Atom *templateAtom = reactantTemplate.getAtomWithIdx(i.first);
    // RDKit❗❌:     int molAtomMapNumber;
    // RDKit❗❌:     if (templateAtom->getPropIfPresent(common_properties::molAtomMapNumber,
    // RDKit❗❌:                                        molAtomMapNumber)) {
    // RDKit❗❌:       if (product->hasAtomBookmark(molAtomMapNumber)) {
    // RDKit❗❌:         RWMol::ATOM_PTR_LIST atomIdxs =
    // RDKit❗❌:             product->getAllAtomsWithBookmark(molAtomMapNumber);
    // RDKit❗❌:         for (auto a : atomIdxs) {
    // RDKit❗❌:           unsigned int pIdx = a->getIdx();
    // RDKit❗❌:           mapping->reactProdAtomMap[i.second].push_back(pIdx);
    // RDKit❗❌:           mapping->mappedAtoms[i.second] = 1;
    // RDKit❗❌:           CHECK_INVARIANT(pIdx < product->getNumAtoms(), "yikes!");
    // RDKit❗❌:           mapping->prodReactAtomMap[pIdx] = i.second;
    // RDKit❗❌:         }
    // RDKit❗❌:       } else {
    // RDKit❗❌:         // this skippedAtom has an atomMapNumber, but it's not in this product
    // RDKit❗❌:         // (it's either in another product or it's not mapped at all).
    // RDKit❗❌:         mapping->skippedAtoms[i.second] = 1;
    // RDKit❗❌:       }
    // RDKit❗❌:     } else {
    // RDKit❗❌:       // This skippedAtom appears in the match, but not in a product:
    // RDKit❗❌:       mapping->skippedAtoms[i.second] = 1;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return mapping;
    // RDKit❗❌: }
    // Source construction precedes all template-bond/property/match reads.
    let mut mapping = new_reactant_product_mapping(reactant_atoms);
    for bond in template.bonds() {
        // RDKit❗✔️:   int getAtomMapNum() const {
        // RDKit❗✔️:     int mapno = 0;
        // RDKit❗✔️:     getPropIfPresent(common_properties::molAtomMapNumber, mapno);
        // RDKit❗✔️:     return mapno;
        // RDKit❗✔️:   }
        let begin = template.atoms().get(bond.begin().index()).ok_or_else(|| {
            invariant(
                "getAtomMappingsReactantProduct",
                "template bond begin atom out of range",
                None,
                None,
                Some(bond.id()),
            )
        })?;
        let mut first = crate::validation::atom_map(begin, ReactionRole::Reactant, template_index)?
            .unwrap_or(0);
        let end = template.atoms().get(bond.end().index()).ok_or_else(|| {
            invariant(
                "getAtomMappingsReactantProduct",
                "template bond end atom out of range",
                None,
                None,
                Some(bond.id()),
            )
        })?;
        let mut second =
            crate::validation::atom_map(end, ReactionRole::Reactant, template_index)?.unwrap_or(0);
        if first > second {
            std::mem::swap(&mut first, &mut second);
        }
        // Native sorts signed map numbers before conversion to unsigned pair.
        mapping.template_bonds.insert((first as u32, second as u32));
    }
    for matched in matched {
        let (query_signed, target_signed) = matched?;
        // Native getAtomWithIdx converts the signed query index to uint32.
        // Its source range precondition precedes property reads and bit access.
        // RDKit❗✔️: const Atom *ROMol::getAtomWithIdx(unsigned int idx) const {
        // RDKit❗✔️:   URANGE_CHECK(idx, getNumAtoms());
        // RDKit❗✔️:
        // RDKit❗✔️:   auto vd = boost::vertex(idx, d_graph);
        // RDKit❗✔️:   const auto res = d_graph[vd];
        // RDKit❗✔️:
        // RDKit❗✔️:   POSTCONDITION(res, "");
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        let query_row = (query_signed as u32) as usize;
        let template_atom = template.atoms().get(query_row).ok_or_else(|| {
            invariant(
                "getAtomMappingsReactantProduct",
                "matched template atom out of range",
                Some((target_signed as u32) as usize),
                None,
                None,
            )
        })?;
        let map =
            crate::validation::atom_map(template_atom, ReactionRole::Reactant, template_index)?;
        let target_key = (target_signed as u32) as usize;
        let target_bit = target_signed as usize;
        if let Some(products) = map.and_then(|map| product.bookmarks.get(&map)) {
            for &product_row in products {
                mapping
                    .reactant_to_product
                    .entry(target_key)
                    .or_default()
                    .push(product_row);
                // Boost❗❌:     reference operator[](size_type pos) {
                // Boost❗❌:         return reference(m_bits[block_index(pos)], bit_index(pos));
                // Boost❗❌:     }
                // Native unchecked bad bit indices are undefined, not a
                // source-defined fallback. Keep a structural error at this
                // exact reached access; an empty bookmark never reaches it.
                *mapping.mapped.get_mut(target_bit).ok_or_else(|| {
                    invariant(
                        "getAtomMappingsReactantProduct",
                        "mapped bit index out of range",
                        Some(target_bit),
                        Some(product_row),
                        None,
                    )
                })? = true;
                if product_row >= product.topology.atoms.len() {
                    return Err(invariant(
                        "getAtomMappingsReactantProduct",
                        "product row out of range",
                        Some(target_key),
                        Some(product_row),
                        None,
                    ));
                }
                mapping.product_to_reactant.insert(product_row, target_key);
            }
        } else {
            *mapping.skipped.get_mut(target_bit).ok_or_else(|| {
                invariant(
                    "getAtomMappingsReactantProduct",
                    "skipped bit index out of range",
                    Some(target_bit),
                    None,
                    None,
                )
            })? = true;
        }
    }
    Ok(mapping)
}

pub(crate) fn atom_mappings(
    matched: &[usize],
    template: &QueryGraph,
    template_index: usize,
    product: &ProductBuilder,
    reactant_atoms: usize,
) -> Result<ReactantProductMapping, ReactionProductError> {
    let matched = matched.iter().enumerate().map(|(query, &target)| {
        Ok((
            source_u32("matched query row", query)? as i32,
            source_u32("matched reactant row", target)? as i32,
        ))
    });
    atom_mappings_source(
        matched,
        template,
        template_index,
        product,
        source_u32("reactant atom count", reactant_atoms)?,
    )
}

pub(crate) fn transfer_bond_properties(
    product: &mut ProductBuilder,
    input: ReactionInput<'_>,
    input_index: usize,
    mapping: &ReactantProductMapping,
) -> Result<(), ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: setReactantBondPropertiesToProduct
    // RDKit❗❌: void setReactantBondPropertiesToProduct(RWMOL_SPTR product,
    // RDKit❗❌:                                         const ROMol &reactant,
    // RDKit❗❌:                                         ReactantProductAtomMapping *mapping) {
    // RDKit❗❌:   for (unsigned int bidx = 0; bidx < product->getNumBonds(); ++bidx) {
    // RDKit❗❌:     auto pBond = product->getBondWithIdx(bidx);
    // RDKit❗❌:     auto rBondBegin = mapping->prodReactAtomMap.find(pBond->getBeginAtomIdx());
    // RDKit❗❌:     auto rBondEnd = mapping->prodReactAtomMap.find(pBond->getEndAtomIdx());
    // RDKit❗❌:
    // RDKit❗❌:     if (rBondBegin == mapping->prodReactAtomMap.end() ||
    // RDKit❗❌:         rBondEnd == mapping->prodReactAtomMap.end()) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     // the bond is between two mapped atoms from this reactant:
    // RDKit❗❌:     const Bond *rBond =
    // RDKit❗❌:         reactant.getBondBetweenAtoms(rBondBegin->second, rBondEnd->second);
    // RDKit❗❌:     if (!rBond) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:     if (!pBond->hasProp(common_properties::NullBond) &&
    // RDKit❗❌:         !pBond->hasProp(common_properties::_MolFileBondQuery) &&
    // RDKit❗❌:         !rBond->hasQuery()) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     if (!rBond->hasQuery()) {
    // RDKit❗❌:       pBond->setBondType(rBond->getBondType());
    // RDKit❗❌:     } else {
    // RDKit❗❌:       QueryBond qBond(rBond->getBondType());
    // RDKit❗❌:       qBond.setQuery(rBond->getQuery()->copy());
    // RDKit❗❌:       // replaceBond copies, so we are safe passing a pointer
    // RDKit❗❌:       // to a local:
    // RDKit❗❌:       product->replaceBond(bidx, &qBond);
    // RDKit❗❌:       pBond = product->getBondWithIdx(bidx);
    // RDKit❗❌:     }
    // RDKit❗❌:     if (rBond->getBondType() == Bond::DOUBLE &&
    // RDKit❗❌:         rBond->getBondDir() == Bond::EITHERDOUBLE) {
    // RDKit❗❌:       pBond->setBondDir(Bond::EITHERDOUBLE);
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     pBond->setIsAromatic(rBond->getIsAromatic());
    // RDKit❗❌:
    // RDKit❗❌:     pBond->updateProps(*rBond);
    // RDKit❗❌:     if (pBond->hasProp(common_properties::NullBond)) {
    // RDKit❗❌:       pBond->clearProp(common_properties::NullBond);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END RDKIT COMPLETE CPP FUNCTION
    // Cost: canonical property-tree/query cloning and extra provenance writes
    // exceed the native vector dictionary/pointer bookkeeping; marker cost ❌.
    // Native replaceBond preserves endpoints and stable bond index; the
    // detached MODEL kernel owns that algorithm for every consumer.
    for row in 0..product.topology.bonds.len() {
        let bond = &product.topology.bonds[row];
        let (Some(&begin), Some(&end)) = (
            mapping.product_to_reactant.get(&bond.begin().index()),
            mapping.product_to_reactant.get(&bond.end().index()),
        ) else {
            continue;
        };
        let Some(original_id) = cosmolkit_model::source_bond_between_atoms(
            input.topology.atoms.len(),
            AtomId::new(begin),
            AtomId::new(end),
            || input.topology.adjacency.neighbors_of(begin).iter(),
        )?
        else {
            continue;
        };
        let original = input
            .topology
            .bonds
            .get(original_id.index())
            .ok_or_else(|| {
                invariant(
                    "setReactantBondPropertiesToProduct",
                    "reactant edge bond row out of range",
                    Some(begin),
                    None,
                    Some(original_id),
                )
            })?;
        if bond.prop("NullBond").is_none()
            && bond.prop("_MolFileBondQuery").is_none()
            && original.query().is_none()
        {
            continue;
        }
        if let Some(query) = original.query() {
            // QueryBond(BondType) starts from ordinary Bond defaults; setQuery
            // replaces its temporary order query with a deep source query copy.
            // No reagent stereo/direction/conjugation is copied into qBond.
            // RDKit❗❌: QueryBond::QueryBond(BondType bT) : Bond(bT) {
            // RDKit❗❌:   if (bT != Bond::UNSPECIFIED) {
            // RDKit❗❌:     dp_query = makeBondOrderEqualsQuery(bT);
            // RDKit❗❌:   } else {
            // RDKit❗❌:     dp_query = makeBondNullQuery();
            // RDKit❗❌:   }
            // RDKit❗❌: };
            // RDKit❗✔️: Bond::Bond(BondType bT) : RDProps() {
            // RDKit❗✔️:   initBond();
            // RDKit❗✔️:   d_bondType = bT;
            // RDKit❗✔️: };
            // RDKit❗✔️: void Bond::initBond() {
            // RDKit❗✔️:   d_bondType = UNSPECIFIED;
            // RDKit❗✔️:   d_dirTag = NONE;
            // RDKit❗✔️:   d_stereo = STEREONONE;
            // RDKit❗✔️:   dp_mol = nullptr;
            // RDKit❗✔️:   d_beginAtomIdx = 0;
            // RDKit❗✔️:   d_endAtomIdx = 0;
            // RDKit❗✔️:   df_isAromatic = 0;
            // RDKit❗✔️:   d_index = 0;
            // RDKit❗✔️:   df_isConjugated = 0;
            // RDKit❗✔️:   dp_stereoAtoms = nullptr;
            // RDKit❗✔️: };
            // RDKit❗✔️:   void setQuery(QUERYBOND_QUERY *what) override {
            // RDKit❗✔️:     // free up any existing query (Issue255):
            // RDKit❗✔️:     delete dp_query;
            // RDKit❗✔️:     dp_query = what;
            // RDKit❗✔️:   }
            let replacement = Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(0), original.order())
                    .with_query(query.clone()),
            );
            cosmolkit_model::replace_source_bond(
                &mut product.topology,
                BondId::new(row),
                &replacement,
                false,
                true,
                |order| {
                    cosmolkit_core::bond_type_as_double(order).map_err(ReactionProductError::from)
                },
                |value| {
                    cosmolkit_core::property_value_to_uint(value)
                        .map_err(ReactionProductError::from)
                },
            )?;
        } else {
            product.topology.bonds[row].set_order(original.order());
        }
        let bond = &mut product.topology.bonds[row];
        if original.order() == BondOrder::Double
            && original.direction() == BondDirection::EitherDouble
        {
            bond.set_direction(BondDirection::EitherDouble);
        }
        bond.set_aromatic(original.is_aromatic());
        bond.update_properties_from(original, false);
        if bond.prop("NullBond").is_some() {
            bond.clear_prop("NullBond")?;
        }
        product.bond_origins[row] = Some(ReactionRowOrigin {
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
        1 => {
            if matches!(
                reactant_tag,
                ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
            ) {
                atom.set_chiral_tag(reactant_tag);
                cosmolkit_core::invert_atom_chirality(atom)?;
            } else {
                eprintln!("unsupported chiral type on reactant atom ignored");
            }
        }
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
    // PDB residue info is the modeled independent monomer capability. Generic
    // non-PDB monomer subclasses have no detached input representation.
    // Presence-only checks precede reached clear/read operations; failed clears
    // must retain prior scalar writes and prevent later bookkeeping writes.
    if product.atomic_number() == 0 || product.prop("_MolFileAtomQuery").is_some() {
        product.set_element(reactant.element());
        product.set_aromatic(reactant.is_aromatic());
        if !implicit {
            product.set_isotope(reactant.isotope());
        }
        if product.prop("dummyLabel").is_some() {
            product.clear_prop("dummyLabel")?;
        }
        if product.prop("_MolFileRLabel").is_some() {
            product.clear_prop("_MolFileRLabel")?;
        }
        product.set_prop("was_dummy", true)?;
    } else if product.prop("was_dummy").is_some() {
        product.clear_prop("was_dummy")?;
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
    ) && crate::management::atom_has_property_source(product, b"molInversionFlag")
    {
        check_product_chirality(reactant.chiral_tag(), product)?;
    }
    if let Some(info) = reactant.pdb_residue_info() {
        // RDKit❗✔️:   AtomMonomerInfo *copy() const override {
        // RDKit❗✔️:     return static_cast<AtomMonomerInfo *>(new AtomPDBResidueInfo(*this));
        // RDKit❗✔️:   }
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
    let origin = Some(ReactionRowOrigin {
        input: input_index,
        row: original.id(),
    });
    if let Some(query) = original.query() {
        // QueryBond(type) defaults plus source setQuery(copy). In particular an
        // aromatic query type must not set aromatic endpoint flags or clear
        // source valence facts: ROMol's pointer overload has no such effects.
        // RDKit❗✔️: QueryBond::QueryBond(BondType bT) : Bond(bT) {
        // RDKit❗✔️:   if (bT != Bond::UNSPECIFIED) {
        // RDKit❗✔️:     dp_query = makeBondOrderEqualsQuery(bT);
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     dp_query = makeBondNullQuery();
        // RDKit❗✔️:   }
        // RDKit❗✔️: };
        // RDKit❗✔️: Bond::Bond(BondType bT) : RDProps() {
        // RDKit❗✔️:   initBond();
        // RDKit❗✔️:   d_bondType = bT;
        // RDKit❗✔️: };
        // RDKit❗✔️: void Bond::initBond() {
        // RDKit❗✔️:   d_bondType = UNSPECIFIED;
        // RDKit❗✔️:   d_dirTag = NONE;
        // RDKit❗✔️:   d_stereo = STEREONONE;
        // RDKit❗✔️:   dp_mol = nullptr;
        // RDKit❗✔️:   d_beginAtomIdx = 0;
        // RDKit❗✔️:   d_endAtomIdx = 0;
        // RDKit❗✔️:   df_isAromatic = 0;
        // RDKit❗✔️:   d_index = 0;
        // RDKit❗✔️:   df_isConjugated = 0;
        // RDKit❗✔️:   dp_stereoAtoms = nullptr;
        // RDKit❗✔️: };
        // RDKit❗✔️:   void setQuery(QUERYBOND_QUERY *what) override {
        // RDKit❗✔️:     // free up any existing query (Issue255):
        // RDKit❗✔️:     delete dp_query;
        // RDKit❗✔️:     dp_query = what;
        // RDKit❗✔️:   }
        let qbond = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(begin), AtomId::new(end), original.order())
                .with_query(query.clone()),
        );
        let count = cosmolkit_model::add_source_bond_value(
            product.topology.atoms.len(),
            &mut product.topology.bonds,
            cosmolkit_model::SourceBondNeighbors {
                original: None,
                appended: &mut product.neighbors,
            },
            std::borrow::Cow::Owned(qbond),
        )?;
        product.bond_origins.push(origin);
        Ok(BondId::new(count - 1))
    } else {
        product.add_bond(begin, end, original.order(), origin)
    }
}

fn add_missing_bonds(
    original: &Bond,
    product: &mut ProductBuilder,
    mapping: &mut ReactantProductMapping,
    input_index: usize,
) -> Result<(), ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: addMissingProductBonds
    // RDKit❗🔝: void addMissingProductBonds(const Bond &origB, RWMOL_SPTR product,
    // RDKit❗🔝:                             ReactantProductAtomMapping *mapping) {
    // RDKit❗🔝:   unsigned int begIdx = origB.getBeginAtomIdx();
    // RDKit❗🔝:   unsigned int endIdx = origB.getEndAtomIdx();
    // RDKit❗🔝:
    // RDKit❗🔝:   std::vector<unsigned> prodBeginIdxs = mapping->reactProdAtomMap[begIdx];
    // RDKit❗🔝:   std::vector<unsigned> prodEndIdxs = mapping->reactProdAtomMap[endIdx];
    // RDKit❗🔝:   CHECK_INVARIANT(prodBeginIdxs.size() == prodEndIdxs.size(),
    // RDKit❗🔝:                   "Different number of start-end points for product bonds.");
    // RDKit❗🔝:   for (unsigned i = 0; i < prodBeginIdxs.size(); i++) {
    // RDKit❗🔝:     addBondToProduct(origB, *product, prodBeginIdxs.at(i), prodEndIdxs.at(i));
    // RDKit❗🔝:   }
    // RDKit❗🔝: }
    // END RDKIT COMPLETE CPP FUNCTION
    // Native operator[] inserts both missing keys before the size invariant.
    mapping
        .reactant_to_product
        .entry(original.begin().index())
        .or_default();
    mapping
        .reactant_to_product
        .entry(original.end().index())
        .or_default();
    let begins = &mapping.reactant_to_product[&original.begin().index()];
    let ends = &mapping.reactant_to_product[&original.end().index()];
    // Cost improvement: add_source_bond cannot mutate mapping, so borrowed
    // lists preserve every source encounter/error effect without native's two
    // temporary vector copies or O(N) additional storage.
    if begins.len() != ends.len() {
        return Err(invariant(
            "addMissingProductBonds",
            "different number of start-end points",
            Some(original.begin().index()),
            None,
            Some(original.id()),
        ));
    }
    let mut row = 0u32;
    while (row as usize) < begins.len() {
        add_source_bond(
            original,
            product,
            begins[row as usize],
            ends[row as usize],
            input_index,
        )?;
        row = row.wrapping_add(1);
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
    // RDKit❗❌: void addMissingProductAtom(const Atom &reactAtom, unsigned reactNeighborIdx,
    // RDKit❗❌:                            unsigned prodNeighborIdx, RWMOL_SPTR product,
    // RDKit❗❌:                            const ROMol &reactant,
    // RDKit❗❌:                            ReactantProductAtomMapping *mapping, unsigned int reactantId) {
    // RDKit❗❌:   Atom *newAtom = nullptr;
    // RDKit❗❌:   if (!reactAtom.hasQuery()) {
    // RDKit❗❌:     newAtom = new Atom(reactAtom);
    // RDKit❌❌:   } else {
    // RDKit❌❌:     newAtom = new QueryAtom(dynamic_cast<const QueryAtom &>(reactAtom));
    // RDKit❗❌:   }
    // RDKit❗❌:   unsigned reactAtomIdx = reactAtom.getIdx();
    // RDKit❗❌:   newAtom->setProp<unsigned int>(common_properties::reactantAtomIdx,
    // RDKit❗❌:                                  reactAtomIdx);
    // RDKit❗❌:   newAtom->setProp<unsigned int>(common_properties::reactantIdx,
    // RDKit❗❌:                                  reactantId);
    // RDKit❗❌:   unsigned productIdx = product->addAtom(newAtom, false, true);
    // RDKit❗❌:   mapping->reactProdAtomMap[reactAtomIdx].push_back(productIdx);
    // RDKit❗❌:   mapping->prodReactAtomMap[productIdx] = reactAtomIdx;
    // RDKit❗❌:   // add the bonds
    // RDKit❗❌:   const Bond *origB =
    // RDKit❗❌:       reactant.getBondBetweenAtoms(reactNeighborIdx, reactAtomIdx);
    // RDKit❗❌:   unsigned int begIdx = productIdx;
    // RDKit❗❌:   unsigned int endIdx = prodNeighborIdx;
    // RDKit❗❌:   if (origB->getBeginAtomIdx() == reactNeighborIdx) {
    // RDKit❗❌:     std::swap(begIdx, endIdx);
    // RDKit❗❌:   }
    // RDKit❗❌:   Bond *prodB = addBondToProduct(*origB, *product, begIdx, endIdx);
    // RDKit❗❌:   if (origB->getBondType() == Bond::DOUBLE &&
    // RDKit❗❌:       origB->getBondDir() == Bond::EITHERDOUBLE) {
    // RDKit❗❌:     prodB->setBondDir(Bond::EITHERDOUBLE);
    // RDKit❗❌:   }
    // RDKit❗❌:   bool preserveExisting = true;
    // RDKit❗❌:   prodB->updateProps(*origB, preserveExisting);
    // RDKit❗❌: }
    // Independent G7 query-bearing Atom reagents are not representable by
    // ReactionInput's canonical concrete Atom rows; that branch stays ❌❌,
    // rather than silently lowering a query to an ordinary atom. Ordinary
    // clone uses MODEL's one canonical atom/property/PDB copy implementation.
    // Cost ❌: detached dictionary-copy and neighbor/provenance storage add
    // work to native graph/pointer construction; no whole product clone.

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
    let original_id = cosmolkit_model::source_bond_between_atoms(
        input.topology.atoms.len(),
        AtomId::new(reactant_neighbor),
        AtomId::new(target),
        || {
            input
                .topology
                .adjacency
                .neighbors_of(reactant_neighbor)
                .iter()
        },
    )?
    .ok_or_else(|| {
        invariant(
            "addMissingProductAtom",
            "missing reactant neighbor bond",
            Some(target),
            Some(row),
            None,
        )
    })?;
    let original = input
        .topology
        .bonds
        .get(original_id.index())
        .ok_or_else(|| {
            invariant(
                "addMissingProductAtom",
                "reactant edge bond row out of range",
                Some(target),
                Some(row),
                Some(original_id),
            )
        })?;
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
    // RDKit❗❌: void addReactantNeighborsToProduct(
    // RDKit❗❌:     const ROMol &reactant, const Atom &reactantAtom, RWMOL_SPTR product,
    // RDKit❗❌:     boost::dynamic_bitset<> &visitedAtoms,
    // RDKit❗❌:     std::vector<const Atom *> &chiralAtomsToCheck,
    // RDKit❗❌:     ReactantProductAtomMapping *mapping, unsigned int reactantId) {
    // RDKit❗❌:   std::list<const Atom *> atomStack;
    // RDKit❗❌:   atomStack.push_back(&reactantAtom);
    // RDKit❗❌:
    // RDKit❗❌:   // std::cerr << "-------------------" << std::endl;
    // RDKit❗❌:   // std::cerr << "  add reactant neighbors from: " << reactantAtom.getIdx()
    // RDKit❗❌:   //           << std::endl;
    // RDKit❗❌:   // #if 1
    // RDKit❗❌:   //   product->updatePropertyCache(false);
    // RDKit❗❌:   //   product->debugMol(std::cerr);
    // RDKit❗❌:   //   std::cerr << "-------------------" << std::endl;
    // RDKit❗❌:   // #endif
    // RDKit❗❌:
    // RDKit❗❌:   while (!atomStack.empty()) {
    // RDKit❗❌:     const Atom *lReactantAtom = atomStack.front();
    // RDKit❗❌:     // std::cerr << "    front: " << lReactantAtom->getIdx() << std::endl;
    // RDKit❗❌:     atomStack.pop_front();
    // RDKit❗❌:
    // RDKit❗❌:     // each atom in the stack is guaranteed to already be in the product:
    // RDKit❗❌:     CHECK_INVARIANT(mapping->reactProdAtomMap.find(lReactantAtom->getIdx()) !=
    // RDKit❗❌:                         mapping->reactProdAtomMap.end(),
    // RDKit❗❌:                     "reactant atom on traversal stack not present in product.");
    // RDKit❗❌:
    // RDKit❗❌:     std::vector<unsigned> lReactantAtomProductIndex =
    // RDKit❗❌:         mapping->reactProdAtomMap[lReactantAtom->getIdx()];
    // RDKit❗❌:     unsigned lreactIdx = lReactantAtom->getIdx();
    // RDKit❗❌:     visitedAtoms[lreactIdx] = 1;
    // RDKit❗❌:     // Check our neighbors:
    // RDKit❗❌:     ROMol::ADJ_ITER nbrIdx, endNbrs;
    // RDKit❗❌:     boost::tie(nbrIdx, endNbrs) = reactant.getAtomNeighbors(lReactantAtom);
    // RDKit❗❌:     while (nbrIdx != endNbrs) {
    // RDKit❗❌:       // Four possibilities here. The neighbor:
    // RDKit❗❌:       //  0) has been visited already: do nothing
    // RDKit❗❌:       //  1) is part of the match (thus already in the product): set a bond to
    // RDKit❗❌:       //  it
    // RDKit❗❌:       //  2) has been added: set a bond to it
    // RDKit❗❌:       //  3) has not yet been added: add it, set a bond to it, and push it
    // RDKit❗❌:       //     onto the stack
    // RDKit❗❌:       // std::cerr << "       nbr: " << *nbrIdx << std::endl;
    // RDKit❗❌:       // std::cerr << "              visited: " << visitedAtoms[*nbrIdx]
    // RDKit❗❌:       //           << "  skipped: " << mapping->skippedAtoms[*nbrIdx]
    // RDKit❗❌:       //           << " mapped: " << mapping->mappedAtoms[*nbrIdx]
    // RDKit❗❌:       //           << " mappedO: " << mapping->mappedAtoms[lreactIdx] <<
    // RDKit❗❌:       //           std::endl;
    // RDKit❗❌:       if (!visitedAtoms[*nbrIdx] && !mapping->skippedAtoms[*nbrIdx]) {
    // RDKit❗❌:         if (mapping->mappedAtoms[*nbrIdx]) {
    // RDKit❗❌:           // this is case 1 (neighbor in match); set a bond to the neighbor if
    // RDKit❗❌:           // this atom
    // RDKit❗❌:           // is not also in the match (match-match bonds were set when the
    // RDKit❗❌:           // product template was
    // RDKit❗❌:           // copied in to start things off).;
    // RDKit❗❌:           if (!mapping->mappedAtoms[lreactIdx]) {
    // RDKit❗❌:             CHECK_INVARIANT(mapping->reactProdAtomMap.find(*nbrIdx) !=
    // RDKit❗❌:                                 mapping->reactProdAtomMap.end(),
    // RDKit❗❌:                             "reactant atom not present in product.");
    // RDKit❗❌:             const Bond *origB =
    // RDKit❗❌:                 reactant.getBondBetweenAtoms(lreactIdx, *nbrIdx);
    // RDKit❗❌:             addMissingProductBonds(*origB, product, mapping);
    // RDKit❗❌:           } else {
    // RDKit❗❌:             // both mapped atoms are in the match.
    // RDKit❗❌:             // they are bonded in the reactant (otherwise we wouldn't be here),
    // RDKit❗❌:             //
    // RDKit❗❌:             // If they do not have already have a bond in the product and did
    // RDKit❗❌:             // not have one in the reactant template then set one here
    // RDKit❗❌:             // If they do have a bond in the reactant template, then we
    // RDKit❗❌:             // assume that this is an intentional bond break, so we don't do
    // RDKit❗❌:             // anything
    // RDKit❗❌:             //
    // RDKit❗❌:             // this was github #1387
    // RDKit❗❌:             unsigned prodBeginIdx = mapping->reactProdAtomMap[lreactIdx][0];
    // RDKit❗❌:             unsigned prodEndIdx = mapping->reactProdAtomMap[*nbrIdx][0];
    // RDKit❗❌:             if (!product->getBondBetweenAtoms(prodBeginIdx, prodEndIdx)) {
    // RDKit❗❌:               // They must be mapped
    // RDKit❗❌:               CHECK_INVARIANT(
    // RDKit❗❌:                   product->getAtomWithIdx(prodBeginIdx)
    // RDKit❗❌:                           ->hasProp(common_properties::reactionMapNum) &&
    // RDKit❗❌:                       product->getAtomWithIdx(prodEndIdx)
    // RDKit❗❌:                           ->hasProp(common_properties::reactionMapNum),
    // RDKit❗❌:                   "atoms should be mapped in product");
    // RDKit❗❌:               int a1mapidx =
    // RDKit❗❌:                   product->getAtomWithIdx(prodBeginIdx)
    // RDKit❗❌:                       ->getProp<int>(common_properties::reactionMapNum);
    // RDKit❗❌:               int a2mapidx =
    // RDKit❗❌:                   product->getAtomWithIdx(prodEndIdx)
    // RDKit❗❌:                       ->getProp<int>(common_properties::reactionMapNum);
    // RDKit❗❌:               if (a1mapidx > a2mapidx) {
    // RDKit❗❌:                 std::swap(a1mapidx, a2mapidx);
    // RDKit❗❌:               }
    // RDKit❗❌:               if (mapping->reactantTemplateAtomBonds.find(
    // RDKit❗❌:                       std::make_pair(a1mapidx, a2mapidx)) ==
    // RDKit❗❌:                   mapping->reactantTemplateAtomBonds.end()) {
    // RDKit❗❌:                 const Bond *origB =
    // RDKit❗❌:                     reactant.getBondBetweenAtoms(lreactIdx, *nbrIdx);
    // RDKit❗❌:                 addMissingProductBonds(*origB, product, mapping);
    // RDKit❗❌:               }
    // RDKit❗❌:             }
    // RDKit❗❌:           }
    // RDKit❗❌:         } else if (mapping->reactProdAtomMap.find(*nbrIdx) !=
    // RDKit❗❌:                    mapping->reactProdAtomMap.end()) {
    // RDKit❗❌:           // case 2, the neighbor has been added and we just need to set a bond
    // RDKit❗❌:           // to it:
    // RDKit❗❌:           const Bond *origB = reactant.getBondBetweenAtoms(lreactIdx, *nbrIdx);
    // RDKit❗❌:           addMissingProductBonds(*origB, product, mapping);
    // RDKit❗❌:         } else {
    // RDKit❗❌:           // case 3, add the atom, a bond to it, and push the atom onto the
    // RDKit❗❌:           // stack
    // RDKit❗❌:           const Atom *neighbor = reactant.getAtomWithIdx(*nbrIdx);
    // RDKit❗❌:           for (unsigned int i : lReactantAtomProductIndex) {
    // RDKit❗❌:             addMissingProductAtom(*neighbor, lreactIdx, i, product, reactant,
    // RDKit❗❌:                                   mapping, reactantId);
    // RDKit❗❌:           }
    // RDKit❗❌:           // update the stack:
    // RDKit❗❌:           atomStack.push_back(neighbor);
    // RDKit❗❌:           // if the atom is chiral, we need to check its bond ordering later:
    // RDKit❗❌:           if (neighbor->getChiralTag() != Atom::CHI_UNSPECIFIED) {
    // RDKit❗❌:             chiralAtomsToCheck.push_back(neighbor);
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       nbrIdx++;
    // RDKit❗❌:     }
    // RDKit❗❌:   }  // end of atomStack traversal
    // RDKit❗❌: }
    // END RDKIT COMPLETE CPP FUNCTION
    // Native keeps pointers on a FIFO list; borrow the canonical Atom rows so
    // each popped pointer's actual getIdx(), rather than a projected CSR row,
    // determines all subsequent mapping/bit/neighbor accesses.
    // Cost ❌: byte bitsets and detached property/provenance/neighbor transport
    // retain the known memory/work gaps despite VecDeque avoiding list nodes.
    let stage = "addReactantNeighborsToProduct";
    let bit = |bits: &[bool], row: usize, detail: &'static str| {
        bits.get(row)
            .copied()
            .ok_or_else(|| invariant(stage, detail, Some(row), None, None))
    };
    let original_bond = |begin: usize, end: usize| -> Result<&Bond, ReactionProductError> {
        let id = cosmolkit_model::source_bond_between_atoms(
            input.topology.atoms.len(),
            AtomId::new(begin),
            AtomId::new(end),
            || input.topology.adjacency.neighbors_of(begin).iter(),
        )?
        .ok_or_else(|| {
            invariant(
                stage,
                "missing reactant neighbor bond",
                Some(end),
                None,
                None,
            )
        })?;
        input.topology.bonds.get(id.index()).ok_or_else(|| {
            invariant(
                stage,
                "reactant edge bond row out of range",
                Some(end),
                None,
                Some(id),
            )
        })
    };
    let start = input.topology.atoms.get(start).ok_or_else(|| {
        invariant(
            stage,
            "start atom row out of range",
            Some(start),
            None,
            None,
        )
    })?;
    let mut stack = VecDeque::from([start]);
    while let Some(current_atom) = stack.pop_front() {
        let current = current_atom.id().index();
        // Snapshot precedes visited writes and must survive mapping changes
        // made while the original neighbor range is being processed.
        let current_products = mapping
            .reactant_to_product
            .get(&current)
            .ok_or_else(|| {
                invariant(
                    stage,
                    "traversal atom not present in product",
                    Some(current),
                    None,
                    None,
                )
            })?
            .clone();
        *visited.get_mut(current).ok_or_else(|| {
            invariant(
                stage,
                "visited bit index out of range",
                Some(current),
                None,
                None,
            )
        })? = true;
        let neighbors = input
            .topology
            .adjacency
            .try_neighbors_of(current)
            .ok_or_else(|| {
                invariant(
                    stage,
                    "incident neighbor row out of range",
                    Some(current),
                    None,
                    None,
                )
            })?;
        for neighbor in neighbors {
            let target = neighbor.atom_index;
            // Native && short-circuit: do not inspect skipped/mapped/edge/atom
            // state when an earlier condition has excluded this neighbor.
            if bit(visited, target, "visited bit index out of range")? {
                continue;
            }
            if bit(&mapping.skipped, target, "skipped bit index out of range")? {
                continue;
            }
            if bit(&mapping.mapped, target, "mapped bit index out of range")? {
                if !bit(&mapping.mapped, current, "mapped bit index out of range")? {
                    if !mapping.reactant_to_product.contains_key(&target) {
                        return Err(invariant(
                            stage,
                            "reactant atom not present in product",
                            Some(target),
                            None,
                            None,
                        ));
                    }
                    add_missing_bonds(
                        original_bond(current, target)?,
                        product,
                        mapping,
                        input_index,
                    )?;
                } else {
                    // Current entry is known to exist from the stack invariant.
                    // Native operator[] creates the target entry before its
                    // unchecked [0]; keep that mutation even on structural error.
                    let begin = mapping.reactant_to_product[&current]
                        .first()
                        .copied()
                        .ok_or_else(|| {
                            invariant(stage, "missing mapped begin row", Some(current), None, None)
                        })?;
                    let end = mapping
                        .reactant_to_product
                        .entry(target)
                        .or_default()
                        .first()
                        .copied()
                        .ok_or_else(|| {
                            invariant(stage, "missing mapped end row", Some(target), None, None)
                        })?;
                    // ProductBuilder construction keeps neighbor rows aligned;
                    // the canonical lookup checks both atom bounds first.
                    let existing = cosmolkit_model::source_bond_between_atoms(
                        product.topology.atoms.len(),
                        AtomId::new(begin),
                        AtomId::new(end),
                        || product.neighbors[begin].iter(),
                    )?;
                    if existing.is_none() {
                        let a = &product.topology.atoms[begin];
                        if a.prop("old_mapno").is_none() {
                            return Err(invariant(
                                stage,
                                "atoms should be mapped in product",
                                Some(current),
                                Some(begin),
                                None,
                            ));
                        }
                        let b = &product.topology.atoms[end];
                        if b.prop("old_mapno").is_none() {
                            return Err(invariant(
                                stage,
                                "atoms should be mapped in product",
                                Some(target),
                                Some(end),
                                None,
                            ));
                        }
                        // Both presence checks complete before either numeric
                        // conversion, exactly as the source CHECK_INVARIANT.
                        let mut first = int_prop(a, "old_mapno")?.ok_or(
                            ReactionProductError::MissingProperty {
                                atom: a.id(),
                                key: "old_mapno",
                            },
                        )?;
                        let mut second = int_prop(b, "old_mapno")?.ok_or(
                            ReactionProductError::MissingProperty {
                                atom: b.id(),
                                key: "old_mapno",
                            },
                        )?;
                        if first > second {
                            std::mem::swap(&mut first, &mut second);
                        }
                        if !mapping
                            .template_bonds
                            .contains(&(first as u32, second as u32))
                        {
                            add_missing_bonds(
                                original_bond(current, target)?,
                                product,
                                mapping,
                                input_index,
                            )?;
                        }
                    }
                }
            } else if mapping.reactant_to_product.contains_key(&target) {
                add_missing_bonds(
                    original_bond(current, target)?,
                    product,
                    mapping,
                    input_index,
                )?;
            } else {
                let atom = input.topology.atoms.get(target).ok_or_else(|| {
                    invariant(
                        stage,
                        "neighbor atom row out of range",
                        Some(target),
                        None,
                        None,
                    )
                })?;
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
                stack.push_back(atom);
                if atom.chiral_tag() != ChiralTag::Unspecified {
                    chiral_to_check.push(atom.id().index());
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod source_add_bond_consumer_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, Element, SourceAtomValenceFacts, TopologyEditError};
    #[test]
    fn source_add_bond_consumer_resets_native_cached_valences_and_retains_provenance() {
        let mut builder = ProductBuilder {
            topology: TopologyBlock::default(),
            neighbors: vec![],
            bookmarks: BTreeMap::new(),
            atom_origins: vec![],
            bond_origins: vec![],
        };
        for i in 0..3 {
            let mut atom = Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C));
            atom.set_source_valence_facts(SourceAtomValenceFacts {
                explicit_valence: 3,
                implicit_valence: 1,
            });
            builder.add_atom(atom, None);
        }
        let id = builder.add_bond(0, 1, BondOrder::Aromatic, None).unwrap();
        assert_eq!(id, BondId::new(0));
        assert_eq!(builder.bond_origins, [None]);
        assert_eq!(builder.bond_between(1, 0), Some(id));
        assert_eq!(
            builder.topology.atoms[0].source_valence_facts(),
            SourceAtomValenceFacts::UNINITIALIZED
        );
        assert_eq!(
            builder.topology.atoms[1].source_valence_facts(),
            SourceAtomValenceFacts::UNINITIALIZED
        );
        assert_eq!(
            builder.topology.atoms[2]
                .source_valence_facts()
                .explicit_valence,
            3
        );
        assert!(builder.topology.atoms[0].is_aromatic() && builder.topology.atoms[1].is_aromatic());
        let error = builder.add_bond(9, 8, BondOrder::Single, None).unwrap_err();
        assert!(
            matches!(error, ReactionProductError::TopologyEdit(TopologyEditError::AtomOutOfRange { atom, .. }) if atom == AtomId::new(9))
        );
        assert_eq!(builder.topology.bonds.len(), 1);
        assert_eq!(builder.bond_origins, [None]);
    }
}

#[cfg(test)]
mod complete_reactant_product_mapping_constructor_source_tests {
    use super::*;

    #[test]
    fn zero_width_has_no_bits_and_all_default_maps_are_empty() {
        let mapping = new_reactant_product_mapping(0);
        assert!(mapping.mapped.is_empty());
        assert!(mapping.skipped.is_empty());
        assert!(mapping.reactant_to_product.is_empty());
        assert!(mapping.product_to_reactant.is_empty());
        assert!(mapping.product_atom_bond.is_empty());
        assert!(mapping.template_bonds.is_empty());
    }

    #[test]
    fn source_resize_initializes_all_bits_false_across_block_boundaries() {
        for width in [1, 7, 8, 31, 32, 63, 64, 65, 129] {
            let mapping = new_reactant_product_mapping(width);
            assert_eq!(mapping.mapped.len(), width as usize);
            assert_eq!(mapping.skipped.len(), width as usize);
            assert!(mapping.mapped.iter().all(|&value| !value));
            assert!(mapping.skipped.iter().all(|&value| !value));
        }
    }

    #[test]
    fn mapped_and_skipped_bits_and_separate_instances_have_independent_storage() {
        let mut first = new_reactant_product_mapping(65);
        let second = new_reactant_product_mapping(65);
        first.mapped[64] = true;
        first.skipped[0] = true;
        first.reactant_to_product.insert(64, vec![3]);
        assert!(!first.skipped[64]);
        assert!(!first.mapped[0]);
        assert!(!second.mapped[64]);
        assert!(!second.skipped[0]);
        assert!(second.reactant_to_product.is_empty());
        assert_eq!(first.mapped.iter().filter(|&&v| v).count(), 1);
        assert_eq!(first.skipped.iter().filter(|&&v| v).count(), 1);
    }
}

#[cfg(test)]
mod complete_update_implicit_atom_properties_source_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, PropertyValue};
    use cosmolkit_types::Element;

    fn atom(
        element: Element,
        charge: i8,
        isotope: Option<u16>,
        hydrogens: u8,
        no_implicit: bool,
    ) -> Atom {
        let mut atom = Atom::from_spec(AtomId::new(0), AtomSpec::new(element));
        atom.set_formal_charge(charge);
        atom.set_isotope(isotope);
        atom.set_explicit_hydrogens(hydrogens);
        atom.set_no_implicit(no_implicit);
        atom
    }
    fn fields(atom: &Atom) -> (i8, Option<u16>, u8, bool) {
        (
            atom.formal_charge(),
            atom.isotope(),
            atom.explicit_hydrogens(),
            atom.no_implicit(),
        )
    }

    #[test]
    fn identity_change_returns_before_all_inheritance_and_retains_product_annotations() {
        let mut product = atom(Element::N, 1, Some(15), 2, true);
        let reactant = atom(Element::C, -2, Some(13), 5, false);
        product
            .set_prop("unrelated", PropertyValue::Bool(false))
            .unwrap();
        update_implicit_atom_properties(&mut product, &reactant);
        assert_eq!(fields(&product), (1, Some(15), 2, true));
        assert_eq!(product.prop("unrelated"), Some(&PropertyValue::Bool(false)));
    }

    #[test]
    fn every_presence_guard_combination_preserves_native_independent_copy_rules_without_value_casts()
     {
        let keys = [
            "_QueryFormalCharge",
            "_QueryIsotope",
            "_ReactionDegreeChanged",
            "_QueryHCount",
        ];
        for mask in 0..16u32 {
            let mut product = atom(Element::C, 1, Some(12), 2, true);
            let reactant = atom(Element::C, -2, Some(13), 5, false);
            for (bit, key) in keys.iter().enumerate() {
                if mask & (1 << bit) != 0 {
                    product
                        .set_prop(*key, PropertyValue::String("wrong-numeric-tag".into()))
                        .unwrap();
                }
            }
            update_implicit_atom_properties(&mut product, &reactant);
            let h_blocked = mask & 12 != 0;
            assert_eq!(
                fields(&product),
                (
                    if mask & 1 != 0 { 1 } else { -2 },
                    if mask & 2 != 0 { Some(12) } else { Some(13) },
                    if h_blocked { 2 } else { 5 },
                    h_blocked,
                ),
                "presence mask {mask}"
            );
            for (bit, key) in keys.iter().enumerate() {
                assert_eq!(product.prop(*key).is_some(), mask & (1 << bit) != 0);
            }
        }
    }

    #[test]
    fn source_zero_values_clear_previous_fields_including_no_implicit_and_isotope() {
        let mut product = atom(Element::C, 2, Some(13), 7, true);
        let reactant = atom(Element::C, 0, None, 0, false);
        update_implicit_atom_properties(&mut product, &reactant);
        assert_eq!(fields(&product), (0, None, 0, false));
    }

    #[test]
    fn reactant_query_annotations_do_not_guard_product_or_copy_generic_properties() {
        let mut product = atom(Element::C, 1, Some(12), 2, true);
        let mut reactant = atom(Element::C, -2, Some(13), 5, false);
        for key in [
            "_QueryFormalCharge",
            "_QueryIsotope",
            "_ReactionDegreeChanged",
            "_QueryHCount",
            "source-only",
        ] {
            reactant.set_prop(key, PropertyValue::Int(0)).unwrap();
        }
        product
            .set_prop("product-only", PropertyValue::Int(9))
            .unwrap();
        update_implicit_atom_properties(&mut product, &reactant);
        assert_eq!(fields(&product), (-2, Some(13), 5, false));
        assert_eq!(product.prop("source-only"), None);
        assert_eq!(product.prop("_QueryHCount"), None);
        assert_eq!(product.prop("product-only"), Some(&PropertyValue::Int(9)));
        assert_eq!(fields(&reactant), (-2, Some(13), 5, false));
    }
}

#[cfg(test)]
mod complete_update_from_template_source_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, PropertyValue, QueryAtom};
    use cosmolkit_types::Element;

    fn atom() -> Atom {
        Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C))
    }
    fn fields(atom: &Atom) -> (i8, Option<u16>, u8, bool) {
        (
            atom.formal_charge(),
            atom.isotope(),
            atom.explicit_hydrogens(),
            atom.no_implicit(),
        )
    }

    #[test]
    fn absent_properties_and_valid_mass_only_leave_fields_and_changed_false() {
        let mut product = atom();
        let mut template = atom();
        let initial = fields(&product);
        assert!(!update_from_template(&template, &mut product).unwrap());
        template
            .set_prop("_QueryMass", PropertyValue::UInt(12))
            .unwrap();
        assert!(!update_from_template(&template, &mut product).unwrap());
        assert_eq!(fields(&product), initial);
        assert_eq!(product.prop("_QueryMass"), None);
    }

    #[test]
    fn native_compare_precedes_narrow_scalar_assignment_and_repeated_call_can_still_report_changed()
    {
        let mut template = atom();
        template
            .set_prop("_QueryFormalCharge", PropertyValue::Int(256))
            .unwrap();
        template
            .set_prop("_QueryHCount", PropertyValue::UInt(256))
            .unwrap();
        template
            .set_prop("_QueryIsotope", PropertyValue::UInt(65536))
            .unwrap();
        let mut product = atom();
        for _ in 0..2 {
            assert!(update_from_template(&template, &mut product).unwrap());
            assert_eq!(fields(&product), (0, None, 0, true));
        }
    }

    #[test]
    fn equal_hydrogen_count_still_sets_no_implicit_then_second_call_is_unchanged() {
        let mut template = atom();
        template
            .set_prop("_QueryFormalCharge", PropertyValue::Int(0))
            .unwrap();
        template
            .set_prop("_QueryHCount", PropertyValue::UInt(0))
            .unwrap();
        template
            .set_prop("_QueryIsotope", PropertyValue::UInt(0))
            .unwrap();
        let mut product = atom();
        assert!(!product.no_implicit());
        assert!(update_from_template(&template, &mut product).unwrap());
        assert!(product.no_implicit());
        assert!(!update_from_template(&template, &mut product).unwrap());
    }

    #[test]
    fn conversion_failure_order_retains_exact_prior_mutation_prefix_including_unused_mass_read() {
        let keys = [
            "_QueryFormalCharge",
            "_QueryHCount",
            "_QueryMass",
            "_QueryIsotope",
        ];
        for stage in 0..4 {
            let mut template = atom();
            template.set_prop(keys[0], PropertyValue::Int(2)).unwrap();
            template.set_prop(keys[1], PropertyValue::UInt(4)).unwrap();
            template.set_prop(keys[2], PropertyValue::UInt(12)).unwrap();
            template.set_prop(keys[3], PropertyValue::UInt(13)).unwrap();
            template
                .set_prop(
                    keys[stage],
                    PropertyValue::String("bad-numeric-value".into()),
                )
                .unwrap();
            let mut product = atom();
            let error = update_from_template(&template, &mut product).unwrap_err();
            match error {
                ReactionProductError::PropertyInt { key, .. } if stage == 0 => {
                    assert_eq!(key, keys[stage])
                }
                ReactionProductError::PropertyUInt { key, .. } if stage != 0 => {
                    assert_eq!(key, keys[stage])
                }
                other => panic!("unexpected source conversion error: {other:?}"),
            }
            assert_eq!(
                fields(&product),
                (
                    if stage == 0 { 0 } else { 2 },
                    None,
                    if stage >= 2 { 4 } else { 0 },
                    stage >= 2
                ),
                "failure stage {stage}"
            );
            assert_eq!(product.prop(keys[stage]), None);
        }
    }

    #[test]
    fn query_atom_and_atom_templates_share_canonical_typed_property_conversion() {
        let mut query = QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C));
        let mut template = atom();
        for (key, value) in [
            ("_QueryFormalCharge", PropertyValue::String("-2".into())),
            ("_QueryHCount", PropertyValue::String("3".into())),
            ("_QueryMass", PropertyValue::UInt(12)),
            ("_QueryIsotope", PropertyValue::String("13".into())),
        ] {
            query.set_prop(key, value.clone()).unwrap();
            template.set_prop(key, value).unwrap();
        }
        let mut from_query = atom();
        let mut from_atom = atom();
        assert!(update_from_template(&query, &mut from_query).unwrap());
        assert!(update_from_template(&template, &mut from_atom).unwrap());
        assert_eq!(fields(&from_query), (-2, Some(13), 3, true));
        assert_eq!(fields(&from_query), fields(&from_atom));
        assert_eq!(from_query.prop("_QueryHCount"), None);
    }
}

#[cfg(test)]
mod complete_convert_template_source_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, BondSpec, PropertyValue, QueryAtom, QueryBond};
    use cosmolkit_types::{BondStereo, Element};

    fn graph(count: usize, edges: &[(usize, usize, BondOrder)]) -> QueryGraph {
        QueryGraph::from_parts(
            (0..count)
                .map(|i| QueryAtom::new(AtomId::new(i), AtomSpec::new(Element::C)))
                .collect(),
            edges
                .iter()
                .enumerate()
                .map(|(i, &(a, b, order))| {
                    QueryBond::new(
                        BondId::new(i),
                        BondSpec::new(AtomId::new(a), AtomId::new(b), order),
                    )
                })
                .collect(),
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn chain() -> QueryGraph {
        graph(
            4,
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Double),
                (2, 3, BondOrder::Single),
            ],
        )
    }
    fn unknown(graph: &QueryGraph) -> Option<PropertyValue> {
        convert_template(graph, 3).unwrap().topology.bonds[1]
            .prop("_UnknownStereoRxnBond")
            .cloned()
    }

    #[test]
    fn zero_negative_and_duplicate_maps_create_ordered_bookmarks_then_clear_maps_and_write_old_mapno()
     {
        let mut template = graph(3, &[]);
        for (row, map) in [0, -2, -2].into_iter().enumerate() {
            template.atoms_mut()[row]
                .set_prop("molAtomMapNumber", PropertyValue::Int(map))
                .unwrap();
        }
        let product = convert_template(&template, 3).unwrap();
        assert_eq!(product.bookmarks.get(&0), Some(&vec![0]));
        assert_eq!(product.bookmarks.get(&-2), Some(&vec![1, 2]));
        for (row, expected) in [0, -2, -2].into_iter().enumerate() {
            let atom = &product.topology.atoms[row];
            assert_eq!(atom.id(), AtomId::new(row));
            assert_eq!(atom.atom_map(), None);
            assert_eq!(atom.prop("molAtomMapNumber"), None);
            assert_eq!(atom.prop("old_mapno"), Some(&PropertyValue::Int(expected)));
            assert!(template.atoms()[row].prop("molAtomMapNumber").is_some());
        }
    }

    #[test]
    fn map_clear_error_propagates_before_reaction_map_write_instead_of_being_discarded() {
        let mut template = graph(1, &[]);
        template.atoms_mut()[0]
            .set_prop("molAtomMapNumber", PropertyValue::Int(1))
            .unwrap();
        template.atoms_mut()[0]
            .set_prop("__computedProps", PropertyValue::Bool(false))
            .unwrap();
        assert!(matches!(
            convert_template(&template, 3),
            Err(ReactionProductError::AtomProperty(_))
        ));
        assert!(template.atoms()[0].prop("molAtomMapNumber").is_some());
        assert_eq!(template.atoms()[0].prop("old_mapno"), None);
    }

    #[test]
    fn inversion_set_four_alone_keeps_template_chirality_and_other_flags_clear_it() {
        for flag in [None, Some(0), Some(1), Some(2), Some(3), Some(4), Some(5)] {
            let mut template = graph(1, &[]);
            template.atoms_mut()[0].set_chiral_tag(ChiralTag::TetrahedralCw);
            template.atoms_mut()[0].set_mol_inversion_flag(flag);
            let product = convert_template(&template, 3).unwrap();
            assert_eq!(
                product.topology.atoms[0].chiral_tag(),
                if flag == Some(4) {
                    ChiralTag::TetrahedralCw
                } else {
                    ChiralTag::Unspecified
                }
            );
        }
    }

    #[test]
    fn self_alias_annotation_reads_apply_to_the_copied_atom_in_source_order() {
        let mut template = graph(1, &[]);
        for (key, value) in [
            ("_QueryFormalCharge", PropertyValue::Int(-2)),
            ("_QueryHCount", PropertyValue::UInt(3)),
            ("_QueryMass", PropertyValue::UInt(12)),
            ("_QueryIsotope", PropertyValue::UInt(13)),
        ] {
            template.atoms_mut()[0].set_prop(key, value).unwrap();
        }
        let product = convert_template(&template, 3).unwrap();
        let atom = &product.topology.atoms[0];
        assert_eq!(
            (
                atom.formal_charge(),
                atom.explicit_hydrogens(),
                atom.no_implicit(),
                atom.isotope()
            ),
            (-2, 3, true, Some(13))
        );
        assert_eq!(template.atoms()[0].formal_charge(), 0);
        let mut copied = template.atoms()[0].try_to_atom().unwrap();
        copied
            .set_prop("_QueryMass", PropertyValue::String("bad".into()))
            .unwrap();
        assert!(matches!(
            update_from_template_source::<Atom>(None, &mut copied),
            Err(ReactionProductError::PropertyUInt {
                key: "_QueryMass",
                ..
            })
        ));
        assert_eq!(
            (
                copied.formal_charge(),
                copied.explicit_hydrogens(),
                copied.no_implicit(),
                copied.isotope()
            ),
            (-2, 3, true, None)
        );
    }

    #[test]
    fn bond_or_uses_destination_aromaticity_and_negated_description_while_carrier_query_stays_absent()
     {
        for aromatic in [false, true] {
            let mut template = graph(2, &[(0, 1, BondOrder::Single)]);
            for atom in template.atoms_mut() {
                atom.set_aromatic(aromatic);
            }
            template.bonds_mut()[0].set_predicate(QueryNode::Not(Box::new(QueryNode::Or(vec![
                QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Single)),
                QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Aromatic)),
            ]))));
            let product = convert_template(&template, 3).unwrap();
            assert_eq!(
                product.topology.bonds[0].order(),
                if aromatic {
                    BondOrder::Aromatic
                } else {
                    BondOrder::Single
                }
            );
            assert_eq!(product.topology.bonds[0].is_aromatic(), aromatic);
        }
        let mut template = graph(2, &[(0, 1, BondOrder::Single)]);
        for atom in template.atoms_mut() {
            atom.set_aromatic(true);
        }
        let old = template.bonds()[0].bond().clone();
        template.bonds_mut()[0] = QueryBond::from_carrier_parts(
            old,
            QueryNode::Or(vec![
                QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Single)),
                QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Aromatic)),
            ]),
        );
        assert!(template.bonds()[0].predicate_is_carrier_derived());
        let product = convert_template(&template, 3).unwrap();
        assert_eq!(product.topology.bonds[0].order(), BondOrder::Single);
        assert!(!product.topology.bonds[0].is_aromatic());
    }

    #[test]
    fn native_property_update_overwrites_new_null_and_unknown_flags_with_existing_source_values() {
        let mut template = graph(2, &[(0, 1, BondOrder::Single)]);
        template.bonds_mut()[0].set_predicate(QueryNode::Predicate(BondQueryPredicate::Any));
        template.bonds_mut()[0]
            .bond_mut()
            .set_prop("NullBond", PropertyValue::Int(7))
            .unwrap();
        assert_eq!(
            convert_template(&template, 3).unwrap().topology.bonds[0].prop("NullBond"),
            Some(&PropertyValue::Int(7))
        );
        let mut template = chain();
        template.bonds_mut()[1]
            .bond_mut()
            .set_prop("_UnknownStereoRxnBond", PropertyValue::Int(0))
            .unwrap();
        assert_eq!(unknown(&template), Some(PropertyValue::Int(0)));
    }

    #[test]
    fn unknown_double_stereo_requires_both_degrees_and_missing_direction_on_either_side() {
        for (left, right, marked) in [
            (BondDirection::None, BondDirection::None, true),
            (BondDirection::EndUpRight, BondDirection::None, true),
            (BondDirection::None, BondDirection::EndDownRight, true),
            (
                BondDirection::EndUpRight,
                BondDirection::EndDownRight,
                false,
            ),
            (BondDirection::Unknown, BondDirection::EndDownRight, true),
        ] {
            let mut template = chain();
            template.bonds_mut()[0].bond_mut().set_direction(left);
            template.bonds_mut()[2].bond_mut().set_direction(right);
            assert_eq!(unknown(&template).is_some(), marked);
        }
        let short = graph(3, &[(0, 1, BondOrder::Double), (1, 2, BondOrder::Single)]);
        assert_eq!(
            convert_template(&short, 3).unwrap().topology.bonds[0].prop("_UnknownStereoRxnBond"),
            None
        );
    }

    #[test]
    fn ordinary_bond_projection_copies_direction_and_dictionary_but_not_template_stereo_or_conjugation()
     {
        let mut template = chain();
        let old = template.bonds_mut()[1].bond_mut();
        old.set_direction(BondDirection::EitherDouble);
        old.set_stereo(BondStereo::E).unwrap();
        old.set_conjugated(true);
        old.set_prop("raw", PropertyValue::String("preserved".into()))
            .unwrap();
        let product = convert_template(&template, 3).unwrap();
        let bond = &product.topology.bonds[1];
        assert_eq!(bond.direction(), BondDirection::EitherDouble);
        assert_eq!(bond.stereo(), BondStereo::None);
        assert!(!bond.is_conjugated());
        assert_eq!(
            bond.prop("raw"),
            Some(&PropertyValue::String("preserved".into()))
        );
        assert!(product.atom_origins.iter().all(Option::is_none));
        assert!(product.bond_origins.iter().all(Option::is_none));
    }
    mod complete_atom_mappings_source_tests {
        use super::*;

        fn mapped_template(maps: &[Option<i32>]) -> QueryGraph {
            let mut template = graph(maps.len(), &[]);
            for (row, map) in maps.iter().enumerate() {
                if let Some(map) = map {
                    template.atoms_mut()[row]
                        .set_prop("molAtomMapNumber", PropertyValue::Int(*map))
                        .unwrap();
                }
            }
            template
        }
        fn product(count: usize, bookmarks: &[(i32, Vec<usize>)]) -> ProductBuilder {
            let mut product = convert_template(&graph(count, &[]), 0).unwrap();
            product.bookmarks = bookmarks.iter().cloned().collect();
            product
        }
        fn run(
            pairs: &[(i32, i32)],
            template: &QueryGraph,
            product: &ProductBuilder,
            count: u32,
        ) -> Result<ReactantProductMapping, ReactionProductError> {
            atom_mappings_source(pairs.iter().copied().map(Ok), template, 4, product, count)
        }

        #[test]
        fn sparse_and_reordered_pair_indices_keep_one_to_many_bookmark_order() {
            let template = mapped_template(&[None, None, Some(1)]);
            let mapping = run(&[(2, 1)], &template, &product(3, &[(1, vec![2, 0])]), 3).unwrap();
            assert_eq!(mapping.mapped, [false, true, false]);
            assert_eq!(mapping.skipped, [false, false, false]);
            assert_eq!(mapping.reactant_to_product.get(&1), Some(&vec![2, 0]));
            assert_eq!(
                mapping.product_to_reactant,
                BTreeMap::from([(0, 1), (2, 1)])
            );
            assert!(mapping.product_atom_bond.is_empty());
        }

        #[test]
        fn duplicated_matches_preserve_all_rows_last_reverse_mapping_and_independent_skip_bits() {
            let template = mapped_template(&[Some(1), Some(2), None]);
            let mapping = run(
                &[(1, 1), (0, 1), (2, 1), (0, 0)],
                &template,
                &product(2, &[(1, vec![0]), (2, vec![1])]),
                2,
            )
            .unwrap();
            assert_eq!(mapping.reactant_to_product.get(&1), Some(&vec![1, 0]));
            assert_eq!(mapping.reactant_to_product.get(&0), Some(&vec![0]));
            assert_eq!(mapping.product_to_reactant.get(&0), Some(&0));
            assert_eq!(mapping.mapped, [true, true]);
            assert_eq!(mapping.skipped, [false, true]);
        }

        #[test]
        fn existing_empty_bookmark_does_not_access_even_negative_target_bit_or_mark_skipped() {
            let template = mapped_template(&[Some(1)]);
            let mapping = run(&[(0, -1)], &template, &product(0, &[(1, vec![])]), 0).unwrap();
            assert!(mapping.mapped.is_empty());
            assert!(mapping.skipped.is_empty());
            assert!(mapping.reactant_to_product.is_empty());
            assert!(matches!(
                run(&[(0, -1)], &template, &product(0, &[]), 0),
                Err(ReactionProductError::Invariant {
                    detail: "skipped bit index out of range",
                    ..
                })
            ));
        }

        #[test]
        fn template_property_error_precedes_unreached_invalid_target_access() {
            let mut template = mapped_template(&[None]);
            template.atoms_mut()[0]
                .set_prop("molAtomMapNumber", PropertyValue::String("bad-map".into()))
                .unwrap();
            assert!(matches!(
                run(&[(0, 99)], &template, &product(0, &[]), 0),
                Err(ReactionProductError::TemplateProperty(_))
            ));
        }

        #[test]
        fn reached_bit_access_precedes_product_check_and_query_range_precedes_both() {
            let template = mapped_template(&[Some(1)]);
            let product = product(1, &[(1, vec![99])]);
            assert!(matches!(
                run(&[(0, 99)], &template, &product, 1),
                Err(ReactionProductError::Invariant {
                    detail: "mapped bit index out of range",
                    ..
                })
            ));
            assert!(matches!(
                run(&[(0, 0)], &template, &product, 1),
                Err(ReactionProductError::Invariant {
                    detail: "product row out of range",
                    ..
                })
            ));
            assert!(matches!(
                run(&[(-1, 99)], &template, &product, 1),
                Err(ReactionProductError::Invariant {
                    detail: "matched template atom out of range",
                    ..
                })
            ));
        }

        #[test]
        fn template_bonds_are_built_before_empty_match_with_signed_sort_then_unsigned_key_conversion()
         {
            let mut template = graph(3, &[(0, 1, BondOrder::Single), (1, 2, BondOrder::Single)]);
            template.atoms_mut()[0]
                .set_prop("molAtomMapNumber", PropertyValue::Int(-2))
                .unwrap();
            template.atoms_mut()[1]
                .set_prop("molAtomMapNumber", PropertyValue::Int(1))
                .unwrap();
            let mapping = run(&[], &template, &product(0, &[]), 0).unwrap();
            assert_eq!(
                mapping.template_bonds,
                BTreeSet::from([((-2i32) as u32, 1), (0, 1)])
            );
            template.atoms_mut()[0]
                .set_prop("molAtomMapNumber", PropertyValue::String("bad-map".into()))
                .unwrap();
            assert!(matches!(
                run(&[], &template, &product(0, &[]), 0),
                Err(ReactionProductError::TemplateProperty(_))
            ));
        }

        #[test]
        fn template_bond_begin_property_is_read_before_end_atom_lookup() {
            let mut template = graph(2, &[(0, 1, BondOrder::Single)]);
            let old = template.bonds()[0].bond().clone().remapped(
                BondId::new(0),
                AtomId::new(0),
                AtomId::new(99),
                None,
            );
            template.bonds_mut()[0] =
                QueryBond::from_parts(old, QueryNode::Predicate(BondQueryPredicate::Any));
            template.atoms_mut()[0]
                .set_prop("molAtomMapNumber", PropertyValue::String("bad-map".into()))
                .unwrap();
            assert!(matches!(
                run(&[], &template, &product(0, &[]), 0),
                Err(ReactionProductError::TemplateProperty(_))
            ));
            template.atoms_mut()[0]
                .set_prop("molAtomMapNumber", PropertyValue::Int(1))
                .unwrap();
            assert!(matches!(
                run(&[], &template, &product(0, &[]), 0),
                Err(ReactionProductError::Invariant {
                    detail: "template bond end atom out of range",
                    ..
                })
            ));
        }

        #[test]
        fn dense_match_projection_delegates_source_pair_kernel_without_a_second_mapping_algorithm()
        {
            let template = mapped_template(&[Some(1), Some(2)]);
            let product = product(2, &[(1, vec![0]), (2, vec![1])]);
            let native = run(&[(0, 1), (1, 0)], &template, &product, 2).unwrap();
            let projected = atom_mappings(&[1, 0], &template, 4, &product, 2).unwrap();
            assert_eq!(projected.mapped, native.mapped);
            assert_eq!(projected.skipped, native.skipped);
            assert_eq!(projected.reactant_to_product, native.reactant_to_product);
            assert_eq!(projected.product_to_reactant, native.product_to_reactant);
            assert_eq!(projected.template_bonds, native.template_bonds);
        }
    }
    mod complete_transfer_bond_properties_source_tests {
        use super::*;
        use cosmolkit_model::{
            AdjacencyList, CoordinateBlock, MoleculeProperties, TopologyEditError,
        };

        fn product(order: BondOrder) -> ProductBuilder {
            let mut product = convert_template(&graph(2, &[(0, 1, order)]), 0).unwrap();
            product.topology.bonds[0].clear_prop("NullBond").unwrap();
            product
        }
        fn input(order: BondOrder, query: bool) -> TopologyBlock {
            let mut topology = product(order).topology;
            if query {
                topology.bonds[0] = Bond::from_spec(
                    BondId::new(0),
                    BondSpec::new(AtomId::new(0), AtomId::new(1), order).with_query(
                        QueryNode::Not(Box::new(QueryNode::Predicate(BondQueryPredicate::Order(
                            order,
                        )))),
                    ),
                );
            }
            topology.adjacency = AdjacencyList::from_topology(2, &topology.bonds);
            topology
        }
        fn mapping() -> ReactantProductMapping {
            let mut mapping = new_reactant_product_mapping(2);
            mapping.product_to_reactant = BTreeMap::from([(0, 0), (1, 1)]);
            mapping
        }
        fn run(
            product: &mut ProductBuilder,
            topology: &TopologyBlock,
            mapping: &ReactantProductMapping,
        ) -> Result<(), ReactionProductError> {
            transfer_bond_properties(
                product,
                ReactionInput {
                    topology,
                    coordinates: &CoordinateBlock::default(),
                    properties: &MoleculeProperties::default(),
                    rings: None,
                    valence: None,
                },
                7,
                mapping,
            )
        }
        #[test]
        fn unmapped_endpoints_and_missing_reactant_bond_skip_all_mutations() {
            let mut product = product(BondOrder::Single);
            product.topology.bonds[0]
                .set_prop("NullBond", PropertyValue::Int(0))
                .unwrap();
            let before = product.topology.clone();
            let mut mapping = mapping();
            mapping.product_to_reactant.remove(&1);
            run(&mut product, &input(BondOrder::Double, true), &mapping).unwrap();
            assert_eq!(product.topology, before);
            mapping.product_to_reactant.insert(1, 1);
            let mut reactant = input(BondOrder::Double, true);
            reactant.bonds.clear();
            reactant.adjacency = AdjacencyList::from_topology(2, &[]);
            run(&mut product, &reactant, &mapping).unwrap();
            assert_eq!(product.topology, before);
            assert_eq!(product.bond_origins, [None]);
        }
        #[test]
        fn ordinary_unmarked_bond_skips_transfer_even_when_reactant_order_and_properties_differ() {
            let mut product = product(BondOrder::Single);
            let before = product.topology.clone();
            let mut reactant = input(BondOrder::Double, false);
            reactant.bonds[0]
                .set_prop("reagent", PropertyValue::Int(8))
                .unwrap();
            run(&mut product, &reactant, &mapping()).unwrap();
            assert_eq!(product.topology, before);
            assert_eq!(product.bond_origins, [None]);
        }
        #[test]
        fn presence_only_markers_trigger_order_and_dictionary_transfer_without_hydrogen_adjustment()
        {
            for key in ["NullBond", "_MolFileBondQuery"] {
                let mut product = product(BondOrder::Single);
                product.topology.atoms[0].set_explicit_hydrogens(3);
                product.topology.bonds[0]
                    .set_prop(key, PropertyValue::String("not-an-int".into()))
                    .unwrap();
                product.topology.bonds[0]
                    .set_prop("product-only", PropertyValue::Int(1))
                    .unwrap();
                product.topology.bonds[0].set_direction(BondDirection::BeginWedge);
                product.topology.bonds[0].set_conjugated(true);
                let mut reactant = input(BondOrder::Triple, false);
                reactant.bonds[0].set_aromatic(true);
                reactant.bonds[0]
                    .set_prop("reactant-only", PropertyValue::Int(2))
                    .unwrap();
                run(&mut product, &reactant, &mapping()).unwrap();
                let bond = &product.topology.bonds[0];
                assert_eq!(bond.order(), BondOrder::Triple);
                assert!(bond.is_aromatic());
                assert!(bond.is_conjugated());
                assert_eq!(bond.direction(), BondDirection::BeginWedge);
                assert_eq!(bond.prop("product-only"), None);
                assert_eq!(bond.prop(key), None);
                assert_eq!(bond.prop("reactant-only"), Some(&PropertyValue::Int(2)));
                assert_eq!(product.topology.atoms[0].explicit_hydrogens(), 3);
                assert_eq!(
                    product.bond_origins,
                    [Some(ReactionRowOrigin {
                        input: 7,
                        row: BondId::new(0)
                    })]
                );
            }
        }
        #[test]
        fn query_replacement_copies_query_resets_old_metadata_and_adjusts_fractional_hydrogens() {
            let mut product = product(BondOrder::Single);
            for atom in &mut product.topology.atoms {
                atom.set_explicit_hydrogens(2);
            }
            product.topology.bonds[0].set_direction(BondDirection::BeginWedge);
            product.topology.bonds[0].set_conjugated(true);
            product.topology.bonds[0]
                .set_stereo(BondStereo::Any)
                .unwrap();
            product.topology.bonds[0]
                .set_prop("product-only", PropertyValue::Int(1))
                .unwrap();
            let mut reactant = input(BondOrder::OneAndHalf, true);
            reactant.bonds[0].set_direction(BondDirection::EndDownRight);
            reactant.bonds[0].set_conjugated(true);
            reactant.bonds[0].set_stereo(BondStereo::Any).unwrap();
            reactant.bonds[0]
                .set_prop("reactant-only", PropertyValue::Int(2))
                .unwrap();
            run(&mut product, &reactant, &mapping()).unwrap();
            let bond = &product.topology.bonds[0];
            assert_eq!(bond.query(), reactant.bonds[0].query());
            assert!(!std::ptr::eq(
                bond.query().unwrap(),
                reactant.bonds[0].query().unwrap()
            ));
            assert_eq!(
                (bond.id(), bond.begin(), bond.end()),
                (BondId::new(0), AtomId::new(0), AtomId::new(1))
            );
            assert_eq!(bond.direction(), BondDirection::None);
            assert_eq!(bond.stereo(), BondStereo::None);
            assert!(!bond.is_conjugated());
            assert_eq!(bond.prop("product-only"), None);
            assert_eq!(bond.prop("reactant-only"), Some(&PropertyValue::Int(2)));
            assert_eq!(
                product
                    .topology
                    .atoms
                    .iter()
                    .map(Atom::explicit_hydrogens)
                    .collect::<Vec<_>>(),
                [1, 1]
            );
            assert_eq!(product.neighbors[0][0].bond, BondId::new(0));
        }
        #[test]
        fn only_double_either_direction_is_inherited_for_ordinary_and_query_branches() {
            for query in [false, true] {
                for (order, direction, expected) in [
                    (
                        BondOrder::Double,
                        BondDirection::EitherDouble,
                        Some(BondDirection::EitherDouble),
                    ),
                    (BondOrder::Single, BondDirection::EitherDouble, None),
                    (BondOrder::Double, BondDirection::BeginWedge, None),
                ] {
                    let mut product = product(BondOrder::Single);
                    product.topology.bonds[0]
                        .set_prop("NullBond", PropertyValue::Int(0))
                        .unwrap();
                    product.topology.bonds[0].set_direction(BondDirection::EndUpRight);
                    let mut reactant = input(order, query);
                    reactant.bonds[0].set_direction(direction);
                    run(&mut product, &reactant, &mapping()).unwrap();
                    assert_eq!(
                        product.topology.bonds[0].direction(),
                        expected.unwrap_or(if query {
                            BondDirection::None
                        } else {
                            BondDirection::EndUpRight
                        })
                    );
                }
            }
        }
        #[test]
        fn malformed_null_clear_propagates_after_native_prior_writes_and_before_provenance() {
            let mut product = product(BondOrder::Single);
            product.topology.bonds[0]
                .set_prop("NullBond", PropertyValue::Int(0))
                .unwrap();
            let mut reactant = input(BondOrder::Double, false);
            reactant.bonds[0].set_aromatic(true);
            reactant.bonds[0]
                .set_prop("NullBond", PropertyValue::Int(9))
                .unwrap();
            reactant.bonds[0]
                .set_prop("__computedProps", PropertyValue::String("bad-list".into()))
                .unwrap();
            assert!(matches!(
                run(&mut product, &reactant, &mapping()),
                Err(ReactionProductError::BondValue(
                    cosmolkit_model::BondValueError::ComputedListKind(_)
                ))
            ));
            assert_eq!(product.topology.bonds[0].order(), BondOrder::Double);
            assert!(product.topology.bonds[0].is_aromatic());
            assert_eq!(
                product.topology.bonds[0].prop("NullBond"),
                Some(&PropertyValue::Int(9))
            );
            assert_eq!(product.bond_origins, [None]);
        }
        #[test]
        fn source_endpoint_range_checks_precede_edge_lookup_and_transfer_branch_selection() {
            let mut product = product(BondOrder::Single);
            let before = product.topology.clone();
            let mut mapping = mapping();
            mapping.product_to_reactant.insert(0, 99);
            mapping.product_to_reactant.insert(1, 98);
            assert!(
                matches!(run(&mut product,&input(BondOrder::Double,false),&mapping),Err(ReactionProductError::TopologyEdit(TopologyEditError::AtomOutOfRange {atom,..})) if atom==AtomId::new(99))
            );
            mapping.product_to_reactant.insert(0, 0);
            assert!(
                matches!(run(&mut product,&input(BondOrder::Double,false),&mapping),Err(ReactionProductError::TopologyEdit(TopologyEditError::AtomOutOfRange {atom,..})) if atom==AtomId::new(98))
            );
            assert_eq!(product.topology, before);
        }
        #[test]
        fn query_replace_order_read_failure_precedes_replacement_and_keeps_product_state() {
            let mut product = product(BondOrder::Single);
            let before = product.topology.clone();
            assert!(matches!(
                run(&mut product, &input(BondOrder::Other, true), &mapping()),
                Err(ReactionProductError::Valence(
                    cosmolkit_core::ValenceError::BadBondType { .. }
                ))
            ));
            assert_eq!(product.topology, before);
            assert_eq!(product.bond_origins, [None]);
        }
    }
}

#[cfg(test)]
mod complete_check_product_chirality_source_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, PropertyValue};
    use cosmolkit_types::Element;
    const TAGS: [ChiralTag; 9] = [
        ChiralTag::Unspecified,
        ChiralTag::TetrahedralCw,
        ChiralTag::TetrahedralCcw,
        ChiralTag::Other,
        ChiralTag::Tetrahedral,
        ChiralTag::Allene,
        ChiralTag::SquarePlanar,
        ChiralTag::TrigonalBipyramidal,
        ChiralTag::Octahedral,
    ];
    fn atom(flag: i32, tag: ChiralTag) -> Atom {
        let mut atom = Atom::from_spec(AtomId::new(3), AtomSpec::new(Element::C));
        atom.set_chiral_tag(tag);
        atom.set_prop("molInversionFlag", PropertyValue::Int(flag))
            .unwrap();
        atom
    }
    #[test]
    fn zero_and_two_copy_every_source_chiral_tag_including_non_tetrahedral_tags() {
        for flag in [0, 2] {
            for reactant in TAGS {
                for product in TAGS {
                    let mut atom = atom(flag, product);
                    check_product_chirality(reactant, &mut atom).unwrap();
                    assert_eq!(atom.chiral_tag(), reactant);
                }
            }
        }
    }
    #[test]
    fn inversion_only_accepts_cw_ccw_and_other_reactant_tags_leave_existing_product_tag() {
        for reactant in TAGS {
            for product in TAGS {
                let mut atom = atom(1, product);
                check_product_chirality(reactant, &mut atom).unwrap();
                assert_eq!(
                    atom.chiral_tag(),
                    match reactant {
                        ChiralTag::TetrahedralCw => ChiralTag::TetrahedralCcw,
                        ChiralTag::TetrahedralCcw => ChiralTag::TetrahedralCw,
                        _ => product,
                    }
                );
            }
        }
    }
    #[test]
    fn destruction_creation_and_unrecognized_flags_have_literal_source_effects() {
        for flag in [3, 4, -1, 5, i32::MIN, i32::MAX] {
            for reactant in TAGS {
                for product in TAGS {
                    let mut atom = atom(flag, product);
                    check_product_chirality(reactant, &mut atom).unwrap();
                    assert_eq!(
                        atom.chiral_tag(),
                        if flag == 3 {
                            ChiralTag::Unspecified
                        } else {
                            product
                        }
                    );
                }
            }
        }
    }
    #[test]
    fn required_flag_missing_or_wrong_kind_errors_before_any_chiral_write() {
        let mut atom = Atom::from_spec(AtomId::new(3), AtomSpec::new(Element::C));
        atom.set_chiral_tag(ChiralTag::Allene);
        assert!(
            matches!(check_product_chirality(ChiralTag::TetrahedralCw,&mut atom),Err(ReactionProductError::MissingProperty {atom: id,key:"molInversionFlag"}) if id==AtomId::new(3))
        );
        assert_eq!(atom.chiral_tag(), ChiralTag::Allene);
        atom.set_prop("molInversionFlag", PropertyValue::String("bad".into()))
            .unwrap();
        assert!(matches!(
            check_product_chirality(ChiralTag::TetrahedralCw, &mut atom),
            Err(ReactionProductError::PropertyInt {
                key: "molInversionFlag",
                ..
            })
        ));
        assert_eq!(atom.chiral_tag(), ChiralTag::Allene);
    }
    #[test]
    fn existing_typed_projection_and_raw_numeric_conversion_preserve_other_atom_state() {
        for raw in [false, true] {
            let mut atom = Atom::from_spec(AtomId::new(3), AtomSpec::new(Element::N));
            atom.set_chiral_tag(ChiralTag::Allene);
            atom.set_chiral_permutation(Some(17));
            atom.set_prop("unrelated", PropertyValue::Int(9)).unwrap();
            if raw {
                atom.set_prop("molInversionFlag", PropertyValue::UInt(2))
                    .unwrap();
            } else {
                atom.set_mol_inversion_flag(Some(2));
            }
            let mut expected = atom.clone();
            expected.set_chiral_tag(ChiralTag::Octahedral);
            check_product_chirality(ChiralTag::Octahedral, &mut atom).unwrap();
            assert_eq!(atom, expected);
        }
    }
}

#[cfg(test)]
mod complete_transfer_atom_properties_source_tests {
    use super::*;
    use cosmolkit_model::{AtomPdbResidueInfo, AtomSpec, PropertyValue};
    use cosmolkit_types::Element;
    fn atom(element: Element) -> Atom {
        Atom::from_spec(AtomId::new(5), AtomSpec::new(element))
    }
    fn reactant() -> Atom {
        let mut atom = atom(Element::N);
        atom.set_aromatic(true);
        atom.set_isotope(Some(15));
        atom.set_formal_charge(-1);
        atom.set_explicit_hydrogens(2);
        atom.set_no_implicit(true);
        atom
    }
    fn info(name: &str) -> AtomPdbResidueInfo {
        AtomPdbResidueInfo::new(name, 17, "GLY", 4, "A", true).with_occupancy(0.75)
    }
    #[test]
    fn dummy_and_presence_only_query_marker_replace_identity_and_apply_implicit_source_rules() {
        for dummy in [false, true] {
            for implicit in [false, true] {
                let mut product = atom(if dummy {
                    Element::from_atomic_number(0).unwrap()
                } else {
                    Element::O
                });
                if !dummy {
                    product
                        .set_prop("_MolFileAtomQuery", PropertyValue::String("bad-tag".into()))
                        .unwrap();
                }
                product.set_isotope(Some(18));
                product.set_prop("dummyLabel", "X").unwrap();
                product
                    .set_prop("_MolFileRLabel", PropertyValue::Int(0))
                    .unwrap();
                transfer_atom_properties(&mut product, &reactant(), implicit, 7).unwrap();
                assert_eq!(product.element(), Element::N);
                assert!(product.is_aromatic());
                assert_eq!(product.isotope(), Some(15));
                assert_eq!(product.formal_charge(), if implicit { -1 } else { 0 });
                assert_eq!(product.explicit_hydrogens(), if implicit { 2 } else { 0 });
                assert_eq!(product.no_implicit(), implicit);
                assert_eq!(product.prop("dummyLabel"), None);
                assert_eq!(product.prop("_MolFileRLabel"), None);
                assert_eq!(product.prop("was_dummy"), Some(&PropertyValue::Bool(true)));
                assert_eq!(
                    product.prop("react_atom_idx"),
                    Some(&PropertyValue::UInt(5))
                );
                assert_eq!(product.prop("react_idx"), Some(&PropertyValue::UInt(7)));
            }
        }
    }
    #[test]
    fn concrete_identity_clears_only_present_bookkeeping_and_identity_change_prevents_implicit_inheritance()
     {
        for implicit in [false, true] {
            let mut product = atom(Element::O);
            product.set_isotope(Some(18));
            product
                .set_prop("was_dummy", PropertyValue::Bool(false))
                .unwrap();
            product.set_prop("dummyLabel", "retained").unwrap();
            transfer_atom_properties(&mut product, &reactant(), implicit, 0).unwrap();
            assert_eq!(product.element(), Element::O);
            assert_eq!(product.isotope(), Some(18));
            assert_eq!(product.formal_charge(), 0);
            assert!(!product.is_aromatic());
            assert_eq!(product.prop("was_dummy"), None);
            assert!(product.prop("dummyLabel").is_some());
        }
    }
    #[test]
    fn absent_clear_keys_do_not_read_malformed_computed_membership() {
        for dummy in [false, true] {
            let mut product = atom(if dummy {
                Element::from_atomic_number(0).unwrap()
            } else {
                Element::C
            });
            product
                .set_prop("__computedProps", PropertyValue::Bool(false))
                .unwrap();
            transfer_atom_properties(&mut product, &reactant(), false, 0).unwrap();
            assert_eq!(
                product.prop("__computedProps"),
                Some(&PropertyValue::Bool(false))
            );
            assert_eq!(product.prop("react_idx"), Some(&PropertyValue::UInt(0)));
        }
    }
    #[test]
    fn reached_dummy_label_clear_errors_preserve_identity_isotope_prefix_and_prevent_bookkeeping() {
        for label in ["dummyLabel", "_MolFileRLabel"] {
            let mut product = atom(Element::from_atomic_number(0).unwrap());
            product.set_prop(label, "X").unwrap();
            product
                .set_prop("__computedProps", PropertyValue::Bool(false))
                .unwrap();
            assert!(matches!(
                transfer_atom_properties(&mut product, &reactant(), false, 0),
                Err(ReactionProductError::AtomProperty(_))
            ));
            assert_eq!(product.element(), Element::N);
            assert_eq!(product.isotope(), Some(15));
            assert!(product.is_aromatic());
            assert!(product.prop(label).is_some());
            assert_eq!(product.prop("was_dummy"), None);
            assert_eq!(product.prop("react_atom_idx"), None);
        }
    }
    #[test]
    fn reached_concrete_bookkeeping_clear_failure_precedes_origin_and_implicit_writes() {
        let mut product = atom(Element::N);
        product
            .set_prop("was_dummy", PropertyValue::Bool(false))
            .unwrap();
        product
            .set_prop("__computedProps", PropertyValue::Bool(false))
            .unwrap();
        assert!(matches!(
            transfer_atom_properties(&mut product, &reactant(), true, 0),
            Err(ReactionProductError::AtomProperty(_))
        ));
        assert_eq!(product.formal_charge(), 0);
        assert_eq!(product.prop("react_atom_idx"), None);
        assert!(product.prop("was_dummy").is_some());
    }
    #[test]
    fn chirality_guards_skip_bad_flags_for_unspecified_other_and_error_before_residue_copy_when_reached()
     {
        for tag in [
            ChiralTag::Unspecified,
            ChiralTag::Other,
            ChiralTag::SquarePlanar,
        ] {
            let mut product = atom(Element::N);
            product.set_chiral_tag(ChiralTag::Allene);
            product
                .set_prop("molInversionFlag", PropertyValue::String("bad".into()))
                .unwrap();
            product.set_pdb_residue_info(Some(info("OLD")));
            let mut reactant = reactant();
            reactant.set_chiral_tag(tag);
            reactant.set_pdb_residue_info(Some(info("NEW")));
            let result = transfer_atom_properties(&mut product, &reactant, true, 7);
            assert_eq!(product.formal_charge(), -1);
            assert_eq!(product.prop("react_idx"), Some(&PropertyValue::UInt(7)));
            assert_eq!(product.chiral_tag(), ChiralTag::Allene);
            if tag == ChiralTag::SquarePlanar {
                assert!(matches!(
                    result,
                    Err(ReactionProductError::PropertyInt {
                        key: "molInversionFlag",
                        ..
                    })
                ));
                assert_eq!(product.pdb_residue_info(), Some(&info("OLD")));
            } else {
                result.unwrap();
                assert_eq!(product.pdb_residue_info(), Some(&info("NEW")));
            }
        }
    }
    #[test]
    fn typed_flag_presence_reuses_canonical_checker_and_residue_copy_occurs_only_when_present() {
        let mut product = atom(Element::N);
        product.set_mol_inversion_flag(Some(1));
        product.set_pdb_residue_info(Some(info("OLD")));
        let mut reactant = reactant();
        reactant.set_chiral_tag(ChiralTag::TetrahedralCw);
        transfer_atom_properties(&mut product, &reactant, false, 0).unwrap();
        assert_eq!(product.chiral_tag(), ChiralTag::TetrahedralCcw);
        assert_eq!(product.pdb_residue_info(), Some(&info("OLD")));
        reactant.set_pdb_residue_info(Some(info("NEW")));
        transfer_atom_properties(&mut product, &reactant, false, 0).unwrap();
        assert_eq!(product.pdb_residue_info(), reactant.pdb_residue_info());
        assert!(!std::ptr::eq(
            product.pdb_residue_info().unwrap(),
            reactant.pdb_residue_info().unwrap()
        ));
    }
    #[test]
    fn source_origin_indices_are_unsigned_and_overflow_errors_retain_reached_prefix() {
        let mut product = atom(Element::from_atomic_number(0).unwrap());
        transfer_atom_properties(&mut product, &reactant(), false, u32::MAX as usize).unwrap();
        assert_eq!(
            product.prop("react_idx"),
            Some(&PropertyValue::UInt(u32::MAX))
        );
        if usize::BITS > 32 {
            let mut product = atom(Element::from_atomic_number(0).unwrap());
            assert!(matches!(
                transfer_atom_properties(&mut product, &reactant(), false, u32::MAX as usize + 1),
                Err(ReactionProductError::RowOverflow {
                    kind: "reactant",
                    ..
                })
            ));
            assert_eq!(product.element(), Element::N);
            assert_eq!(
                product.prop("react_atom_idx"),
                Some(&PropertyValue::UInt(5))
            );
            assert_eq!(product.prop("react_idx"), None);
        }
    }
}

#[cfg(test)]
mod complete_add_source_bond_source_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, PropertyValue, SourceAtomValenceFacts};
    use cosmolkit_types::{BondStereo, Element};
    fn product() -> ProductBuilder {
        ProductBuilder {
            topology: TopologyBlock {
                atoms: (0..3)
                    .map(|i| {
                        let mut atom = Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C));
                        atom.set_source_valence_facts(SourceAtomValenceFacts {
                            explicit_valence: 3,
                            implicit_valence: 1,
                        });
                        atom
                    })
                    .collect(),
                ..TopologyBlock::default()
            },
            neighbors: vec![vec![]; 3],
            bookmarks: BTreeMap::new(),
            atom_origins: vec![None; 3],
            bond_origins: vec![],
        }
    }
    fn original(order: BondOrder, query: bool) -> Bond {
        let spec = BondSpec::new(AtomId::new(8), AtomId::new(9), order)
            .with_aromatic(true)
            .with_conjugated(true)
            .with_direction(BondDirection::BeginWedge)
            .with_stereo(BondStereo::Any)
            .with_prop("original-only", PropertyValue::Int(8))
            .unwrap();
        Bond::from_spec(
            BondId::new(11),
            if query {
                spec.with_query(QueryNode::Not(Box::new(QueryNode::Predicate(
                    BondQueryPredicate::Any,
                ))))
            } else {
                spec
            },
        )
    }
    #[test]
    fn ordinary_bond_only_inherits_type_and_uses_source_order_overload_cache_reset() {
        let mut product = product();
        let id =
            add_source_bond(&original(BondOrder::Double, false), &mut product, 2, 0, 7).unwrap();
        let bond = &product.topology.bonds[0];
        assert_eq!(id, BondId::new(0));
        assert_eq!(bond.order(), BondOrder::Double);
        assert!(bond.query().is_none());
        assert_eq!((bond.begin(), bond.end()), (AtomId::new(2), AtomId::new(0)));
        assert!(!bond.is_aromatic());
        assert!(!bond.is_conjugated());
        assert_eq!(bond.stereo(), BondStereo::None);
        assert_eq!(bond.direction(), BondDirection::None);
        assert_eq!(bond.prop("original-only"), None);
        for row in [0, 2] {
            assert_eq!(
                product.topology.atoms[row]
                    .source_valence_facts()
                    .explicit_valence,
                -1
            );
        }
        assert_eq!(
            product.topology.atoms[1]
                .source_valence_facts()
                .explicit_valence,
            3
        );
        assert_eq!(
            product.bond_origins,
            [Some(ReactionRowOrigin {
                input: 7,
                row: BondId::new(11)
            })]
        );
    }
    #[test]
    fn aromatic_query_uses_pointer_overload_without_atom_flags_or_cache_effects() {
        let mut product = product();
        let atoms = product.topology.atoms.clone();
        let original = original(BondOrder::Aromatic, true);
        add_source_bond(&original, &mut product, 1, 2, 7).unwrap();
        let bond = &product.topology.bonds[0];
        assert_eq!(product.topology.atoms, atoms);
        assert_eq!(bond.order(), BondOrder::Aromatic);
        assert!(!bond.is_aromatic());
        assert_eq!(bond.query(), original.query());
        assert!(!std::ptr::eq(
            bond.query().unwrap(),
            original.query().unwrap()
        ));
        assert_eq!(bond.direction(), BondDirection::None);
        assert_eq!(bond.stereo(), BondStereo::None);
        assert!(!bond.is_conjugated());
        assert_eq!(bond.prop("original-only"), None);
    }
    #[test]
    fn ordinary_aromatic_type_sets_only_source_endpoint_aromatic_flags() {
        let mut product = product();
        add_source_bond(&original(BondOrder::Aromatic, false), &mut product, 0, 2, 0).unwrap();
        assert!(product.topology.bonds[0].is_aromatic());
        assert!(product.topology.atoms[0].is_aromatic());
        assert!(product.topology.atoms[2].is_aromatic());
        assert!(!product.topology.atoms[1].is_aromatic());
    }
    #[test]
    fn appends_keep_stable_rows_origin_alignment_and_query_other_order_without_valence_read() {
        let mut product = product();
        add_source_bond(&original(BondOrder::Single, false), &mut product, 0, 1, 3).unwrap();
        assert_eq!(
            add_source_bond(&original(BondOrder::Other, true), &mut product, 1, 2, 4).unwrap(),
            BondId::new(1)
        );
        assert_eq!(product.topology.bonds[1].order(), BondOrder::Other);
        assert!(product.topology.bonds[1].query().is_some());
        assert_eq!(
            product
                .bond_origins
                .iter()
                .map(|x| x.unwrap().input)
                .collect::<Vec<_>>(),
            [3, 4]
        );
        assert_eq!(
            product.neighbors[1]
                .iter()
                .map(|x| x.bond)
                .collect::<Vec<_>>(),
            [BondId::new(0), BondId::new(1)]
        );
    }
    #[test]
    fn both_overloads_preserve_rows_neighbors_and_origins_on_range_self_or_duplicate_error() {
        for query in [false, true] {
            let mut product = product();
            let before = product.topology.clone();
            for (begin, end) in [(9, 8), (0, 8), (1, 1)] {
                assert!(
                    add_source_bond(
                        &original(BondOrder::Aromatic, query),
                        &mut product,
                        begin,
                        end,
                        7
                    )
                    .is_err()
                );
                assert_eq!(product.topology, before);
                assert!(product.bond_origins.is_empty());
                assert!(product.neighbors.iter().all(Vec::is_empty));
            }
            add_source_bond(&original(BondOrder::Single, query), &mut product, 0, 1, 7).unwrap();
            let before = product.topology.clone();
            let neighbors = product.neighbors.clone();
            let origins = product.bond_origins.clone();
            assert!(matches!(
                add_source_bond(&original(BondOrder::Aromatic, query), &mut product, 1, 0, 8),
                Err(ReactionProductError::TopologyEdit(
                    cosmolkit_model::TopologyEditError::DuplicateBond { .. }
                ))
            ));
            assert_eq!(product.topology, before);
            assert_eq!(product.neighbors, neighbors);
            assert_eq!(product.bond_origins, origins);
        }
    }
}

#[cfg(test)]
mod complete_add_missing_bonds_source_tests {
    use super::*;
    use cosmolkit_model::AtomSpec;
    use cosmolkit_types::Element;
    fn product() -> ProductBuilder {
        ProductBuilder {
            topology: TopologyBlock {
                atoms: (0..4)
                    .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
                    .collect(),
                ..TopologyBlock::default()
            },
            neighbors: vec![vec![]; 4],
            bookmarks: BTreeMap::new(),
            atom_origins: vec![None; 4],
            bond_origins: vec![],
        }
    }
    fn original() -> Bond {
        Bond::from_spec(
            BondId::new(9),
            BondSpec::new(AtomId::new(4), AtomId::new(5), BondOrder::Single),
        )
    }
    fn mapping(begins: Option<Vec<usize>>, ends: Option<Vec<usize>>) -> ReactantProductMapping {
        let mut mapping = new_reactant_product_mapping(6);
        if let Some(x) = begins {
            mapping.reactant_to_product.insert(4, x);
        }
        if let Some(x) = ends {
            mapping.reactant_to_product.insert(5, x);
        }
        mapping
    }
    #[test]
    fn absent_keys_are_inserted_even_when_both_lists_empty_and_same_key_is_inserted_once() {
        let mut product = product();
        let before = product.topology.clone();
        let mut mapping = mapping(None, None);
        add_missing_bonds(&original(), &mut product, &mut mapping, 7).unwrap();
        assert_eq!(
            mapping.reactant_to_product,
            BTreeMap::from([(4, vec![]), (5, vec![])])
        );
        assert_eq!(product.topology, before);
        let original = Bond::from_spec(
            BondId::new(9),
            BondSpec::new(AtomId::new(4), AtomId::new(4), BondOrder::Single),
        );
        let mut same = new_reactant_product_mapping(5);
        add_missing_bonds(&original, &mut product, &mut same, 0).unwrap();
        assert_eq!(same.reactant_to_product, BTreeMap::from([(4, vec![])]));
    }
    #[test]
    fn size_invariant_runs_after_both_insertions_and_before_product_or_bond_access() {
        for (begins, ends) in [(Some(vec![99]), None), (None, Some(vec![99]))] {
            let mut product = product();
            product.neighbors.clear();
            let before = product.topology.clone();
            let mut mapping = mapping(begins, ends);
            assert!(matches!(
                add_missing_bonds(&original(), &mut product, &mut mapping, 0),
                Err(ReactionProductError::Invariant {
                    stage: "addMissingProductBonds",
                    detail: "different number of start-end points",
                    ..
                })
            ));
            assert!(mapping.reactant_to_product.contains_key(&4));
            assert!(mapping.reactant_to_product.contains_key(&5));
            assert_eq!(product.topology, before);
            assert!(product.bond_origins.is_empty());
        }
    }
    #[test]
    fn aligned_lists_pair_by_physical_position_without_cartesian_product_or_sorting() {
        let mut product = product();
        let mut mapping = mapping(Some(vec![2, 0]), Some(vec![3, 1]));
        let before = mapping.reactant_to_product.clone();
        add_missing_bonds(&original(), &mut product, &mut mapping, 7).unwrap();
        assert_eq!(
            product
                .topology
                .bonds
                .iter()
                .map(|b| (b.begin().index(), b.end().index()))
                .collect::<Vec<_>>(),
            [(2, 3), (0, 1)]
        );
        assert_eq!(
            product.bond_origins,
            [Some(ReactionRowOrigin {
                input: 7,
                row: BondId::new(9)
            }); 2]
        );
        assert_eq!(mapping.reactant_to_product, before);
    }
    #[test]
    fn later_add_failure_keeps_prior_bond_neighbor_and_origin_prefix_without_visiting_later_pair() {
        let mut product = product();
        let mut mapping = mapping(Some(vec![0, 0, 2]), Some(vec![1, 1, 3]));
        assert!(matches!(
            add_missing_bonds(&original(), &mut product, &mut mapping, 7),
            Err(ReactionProductError::TopologyEdit(
                cosmolkit_model::TopologyEditError::DuplicateBond { .. }
            ))
        ));
        assert_eq!(product.topology.bonds.len(), 1);
        assert_eq!(product.neighbors[0].len(), 1);
        assert!(product.neighbors[2].is_empty());
        assert_eq!(product.bond_origins.len(), 1);
        assert_eq!(mapping.reactant_to_product.get(&4), Some(&vec![0, 0, 2]));
    }
    #[test]
    fn query_bond_pairs_use_the_shared_pointer_overload_and_retain_construction_defaults() {
        let mut product = product();
        let atoms = product.topology.atoms.clone();
        let mut mapping = mapping(Some(vec![0, 2]), Some(vec![1, 3]));
        let original = Bond::from_spec(
            BondId::new(9),
            BondSpec::new(AtomId::new(4), AtomId::new(5), BondOrder::Aromatic)
                .with_query(QueryNode::Predicate(BondQueryPredicate::Any)),
        );
        add_missing_bonds(&original, &mut product, &mut mapping, 7).unwrap();
        assert_eq!(product.topology.atoms, atoms);
        assert!(
            product
                .topology
                .bonds
                .iter()
                .all(|b| b.query() == original.query() && !b.is_aromatic())
        );
    }
}

#[cfg(test)]
mod complete_add_missing_atom_source_tests {
    use super::*;
    use cosmolkit_model::{
        AdjacencyList, AtomPdbResidueInfo, AtomSpec, CoordinateBlock, MoleculeProperties,
        PropertyValue,
    };
    use cosmolkit_types::{BondStereo, Element};
    fn input(
        reversed: bool,
        query: bool,
        order: BondOrder,
        direction: BondDirection,
    ) -> TopologyBlock {
        let mut atom = Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::N));
        atom.set_isotope(Some(15));
        atom.set_chiral_tag(ChiralTag::Allene);
        atom.set_prop("source-only", PropertyValue::Int(7)).unwrap();
        atom.set_pdb_residue_info(Some(AtomPdbResidueInfo::new(
            " N  ", 2, "GLY", 3, "A", false,
        )));
        let (begin, end) = if reversed { (1, 0) } else { (0, 1) };
        let spec = BondSpec::new(AtomId::new(begin), AtomId::new(end), order)
            .with_direction(direction)
            .with_stereo(BondStereo::Any)
            .with_aromatic(true)
            .with_conjugated(true)
            .with_prop("source-bond", PropertyValue::Int(8))
            .unwrap()
            .with_prop("NullBond", PropertyValue::Int(0))
            .unwrap();
        let bond = Bond::from_spec(
            BondId::new(0),
            if query {
                spec.with_query(QueryNode::Predicate(BondQueryPredicate::Any))
            } else {
                spec
            },
        );
        TopologyBlock {
            atoms: vec![
                Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
                atom,
            ],
            bonds: vec![bond.clone()],
            adjacency: AdjacencyList::from_topology(2, &[bond]),
            ..TopologyBlock::default()
        }
    }
    fn product() -> ProductBuilder {
        ProductBuilder {
            topology: TopologyBlock {
                atoms: vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C))],
                ..TopologyBlock::default()
            },
            neighbors: vec![vec![]],
            bookmarks: BTreeMap::new(),
            atom_origins: vec![None],
            bond_origins: vec![],
        }
    }
    fn run(
        reactant: &TopologyBlock,
        atom: &Atom,
        neighbor: usize,
        product_neighbor: usize,
        product: &mut ProductBuilder,
        mapping: &mut ReactantProductMapping,
    ) -> Result<(), ReactionProductError> {
        add_missing_atom(
            atom,
            neighbor,
            product_neighbor,
            product,
            ReactionInput {
                topology: reactant,
                coordinates: &CoordinateBlock::default(),
                properties: &MoleculeProperties::default(),
                rings: None,
                valence: None,
            },
            mapping,
            7,
        )
    }
    #[test]
    fn ordinary_atom_copy_origin_maps_and_source_bond_orientation_are_preserved() {
        for reversed in [false, true] {
            let input = input(
                reversed,
                false,
                BondOrder::Single,
                BondDirection::BeginWedge,
            );
            let before = input.clone();
            let mut product = product();
            let mut mapping = new_reactant_product_mapping(2);
            run(&input, &input.atoms[1], 0, 0, &mut product, &mut mapping).unwrap();
            let atom = &product.topology.atoms[1];
            assert_eq!(atom.element(), Element::N);
            assert_eq!(atom.isotope(), Some(15));
            assert_eq!(atom.chiral_tag(), ChiralTag::Allene);
            assert_eq!(atom.pdb_residue_info(), input.atoms[1].pdb_residue_info());
            assert_eq!(atom.prop("source-only"), Some(&PropertyValue::Int(7)));
            assert_eq!(atom.prop("react_atom_idx"), Some(&PropertyValue::UInt(1)));
            assert_eq!(atom.prop("react_idx"), Some(&PropertyValue::UInt(7)));
            assert_eq!(
                product.atom_origins,
                [
                    None,
                    Some(ReactionRowOrigin {
                        input: 7,
                        row: AtomId::new(1)
                    })
                ]
            );
            assert_eq!(mapping.reactant_to_product, BTreeMap::from([(1, vec![1])]));
            assert_eq!(mapping.product_to_reactant, BTreeMap::from([(1, 1)]));
            let bond = &product.topology.bonds[0];
            assert_eq!(
                (bond.begin().index(), bond.end().index()),
                if reversed { (1, 0) } else { (0, 1) }
            );
            assert_eq!(bond.direction(), BondDirection::None);
            assert_eq!(bond.stereo(), BondStereo::None);
            assert!(!bond.is_aromatic());
            assert!(!bond.is_conjugated());
            assert_eq!(bond.prop("source-bond"), Some(&PropertyValue::Int(8)));
            assert_eq!(bond.prop("NullBond"), Some(&PropertyValue::Int(0)));
            assert_eq!(input, before);
        }
    }
    #[test]
    fn source_double_either_direction_and_query_copy_use_shared_bond_kernel() {
        for query in [false, true] {
            for (order, direction, expected) in [
                (
                    BondOrder::Double,
                    BondDirection::EitherDouble,
                    BondDirection::EitherDouble,
                ),
                (
                    BondOrder::Single,
                    BondDirection::EitherDouble,
                    BondDirection::None,
                ),
                (
                    BondOrder::Double,
                    BondDirection::EndUpRight,
                    BondDirection::None,
                ),
            ] {
                let input = input(false, query, order, direction);
                let mut product = product();
                let mut mapping = new_reactant_product_mapping(2);
                run(&input, &input.atoms[1], 0, 0, &mut product, &mut mapping).unwrap();
                let bond = &product.topology.bonds[0];
                assert_eq!(bond.direction(), expected);
                assert_eq!(bond.query(), input.bonds[0].query());
                assert_eq!(bond.prop("NullBond"), Some(&PropertyValue::Int(0)));
            }
        }
    }
    #[test]
    fn missing_original_bond_retains_added_atom_and_bidirectional_mapping_before_structural_null_error()
     {
        let mut input = input(false, false, BondOrder::Single, BondDirection::None);
        input.bonds.clear();
        input.adjacency = AdjacencyList::from_topology(2, &[]);
        let mut product = product();
        let mut mapping = new_reactant_product_mapping(2);
        assert!(matches!(
            run(&input, &input.atoms[1], 0, 0, &mut product, &mut mapping),
            Err(ReactionProductError::Invariant {
                detail: "missing reactant neighbor bond",
                ..
            })
        ));
        assert_eq!(product.topology.atoms.len(), 2);
        assert_eq!(product.atom_origins.len(), 2);
        assert_eq!(mapping.reactant_to_product.get(&1), Some(&vec![1]));
        assert_eq!(mapping.product_to_reactant.get(&1), Some(&1));
        assert!(product.topology.bonds.is_empty());
    }
    #[test]
    fn original_endpoint_range_failure_occurs_after_atom_mapping_writes_with_begin_checked_first() {
        let input = input(false, false, BondOrder::Single, BondDirection::None);
        let foreign = input.atoms[1].clone().with_id(AtomId::new(9));
        let mut product = product();
        let mut mapping = new_reactant_product_mapping(2);
        assert!(
            matches!(run(&input,&foreign,99,0,&mut product,&mut mapping),Err(ReactionProductError::TopologyEdit(cosmolkit_model::TopologyEditError::AtomOutOfRange {atom,..})) if atom==AtomId::new(99))
        );
        assert_eq!(product.topology.atoms.len(), 2);
        assert_eq!(mapping.reactant_to_product.get(&9), Some(&vec![1]));
        assert_eq!(mapping.product_to_reactant.get(&1), Some(&9));
    }
    #[test]
    fn product_bond_failure_occurs_after_atom_append_without_undoing_origin_or_mapping_prefix() {
        let input = input(false, false, BondOrder::Single, BondDirection::None);
        let mut product = product();
        let mut mapping = new_reactant_product_mapping(2);
        assert!(
            matches!(run(&input,&input.atoms[1],0,99,&mut product,&mut mapping),Err(ReactionProductError::TopologyEdit(cosmolkit_model::TopologyEditError::AtomOutOfRange {atom,..})) if atom==AtomId::new(99))
        );
        assert_eq!(product.topology.atoms.len(), 2);
        assert_eq!(product.neighbors.len(), 2);
        assert_eq!(product.atom_origins.len(), 2);
        assert_eq!(mapping.reactant_to_product.get(&1), Some(&vec![1]));
        assert!(product.bond_origins.is_empty());
    }
    #[test]
    fn repeated_atom_add_appends_mapping_rows_and_never_changes_mapped_or_skipped_bits() {
        let input = input(false, false, BondOrder::Single, BondDirection::None);
        let mut product = product();
        let mut mapping = new_reactant_product_mapping(2);
        mapping.mapped[0] = true;
        mapping.skipped[1] = true;
        mapping.reactant_to_product.insert(1, vec![99]);
        for _ in 0..2 {
            run(&input, &input.atoms[1], 0, 0, &mut product, &mut mapping).unwrap();
        }
        assert_eq!(mapping.reactant_to_product.get(&1), Some(&vec![99, 1, 2]));
        assert_eq!(
            mapping.product_to_reactant,
            BTreeMap::from([(1, 1), (2, 1)])
        );
        assert_eq!(mapping.mapped, [true, false]);
        assert_eq!(mapping.skipped, [false, true]);
        assert_eq!(product.topology.bonds.len(), 2);
    }
}

#[cfg(test)]
mod complete_add_neighbors_source_tests {
    use super::*;
    use cosmolkit_model::{
        AdjacencyList, AtomSpec, CoordinateBlock, MoleculeProperties, PropertyValue,
    };
    use cosmolkit_types::Element;
    fn input(count: usize, edges: &[(usize, usize)]) -> TopologyBlock {
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(row, &(begin, end))| {
                Bond::from_spec(
                    BondId::new(row),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                )
            })
            .collect::<Vec<_>>();
        TopologyBlock {
            atoms: (0..count)
                .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
                .collect(),
            adjacency: AdjacencyList::from_topology(count, &bonds),
            bonds,
            ..TopologyBlock::default()
        }
    }
    fn product(count: usize) -> ProductBuilder {
        ProductBuilder {
            topology: TopologyBlock {
                atoms: (0..count)
                    .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
                    .collect(),
                ..TopologyBlock::default()
            },
            neighbors: vec![vec![]; count],
            bookmarks: BTreeMap::new(),
            atom_origins: vec![None; count],
            bond_origins: vec![],
        }
    }
    fn mapping(count: usize, rows: &[(usize, Vec<usize>)]) -> ReactantProductMapping {
        let mut mapping = new_reactant_product_mapping(count as u32);
        mapping.reactant_to_product = rows.iter().cloned().collect();
        mapping
    }
    fn run(
        input: &TopologyBlock,
        product: &mut ProductBuilder,
        visited: &mut [bool],
        chiral: &mut Vec<usize>,
        mapping: &mut ReactantProductMapping,
    ) -> Result<(), ReactionProductError> {
        add_neighbors(
            ReactionInput {
                topology: input,
                coordinates: &CoordinateBlock::default(),
                properties: &MoleculeProperties::default(),
                rings: None,
                valence: None,
            },
            0,
            product,
            visited,
            chiral,
            mapping,
            7,
        )
    }
    fn malformed_edge() -> TopologyBlock {
        let mut input = input(2, &[(0, 1)]);
        input.bonds[0] =
            input.bonds[0]
                .clone()
                .remapped(BondId::new(99), AtomId::new(0), AtomId::new(1), None);
        input.adjacency = AdjacencyList::from_topology(2, &input.bonds);
        input
    }
    fn both_mapped() -> ReactantProductMapping {
        let mut mapping = mapping(2, &[(0, vec![0]), (1, vec![1])]);
        mapping.mapped.fill(true);
        mapping
    }
    #[test]
    fn fifo_neighbor_encounter_order_controls_product_rows_and_chiral_check_order() {
        let mut input = input(4, &[(0, 1), (0, 2), (1, 3)]);
        for atom in &mut input.atoms[1..] {
            atom.set_chiral_tag(ChiralTag::Other);
        }
        let mut product = product(1);
        let mut mapping = mapping(4, &[(0, vec![0])]);
        mapping.mapped[0] = true;
        let mut visited = vec![false; 4];
        let mut chiral = vec![];
        run(
            &input,
            &mut product,
            &mut visited,
            &mut chiral,
            &mut mapping,
        )
        .unwrap();
        assert_eq!(chiral, [1, 2, 3]);
        assert_eq!(visited, [true; 4]);
        assert_eq!(
            product.atom_origins[1..]
                .iter()
                .map(|x| x.unwrap().row.index())
                .collect::<Vec<_>>(),
            [1, 2, 3]
        );
        assert_eq!(
            product
                .topology
                .bonds
                .iter()
                .map(|b| (b.begin().index(), b.end().index()))
                .collect::<Vec<_>>(),
            [(0, 1), (0, 2), (1, 3)]
        );
    }
    #[test]
    fn visited_then_skipped_short_circuit_never_reads_later_bitsets_or_invalid_edge_row() {
        for already_visited in [false, true] {
            let input = malformed_edge();
            let mut product = product(1);
            let before = product.topology.clone();
            let mut mapping = mapping(2, &[(0, vec![0])]);
            mapping.mapped.clear();
            let mut visited = vec![false, already_visited];
            if already_visited {
                mapping.skipped.clear();
            } else {
                mapping.skipped[1] = true;
            }
            run(
                &input,
                &mut product,
                &mut visited,
                &mut vec![],
                &mut mapping,
            )
            .unwrap();
            assert_eq!(product.topology, before);
            assert!(visited[0]);
            assert_eq!(visited[1], already_visited);
        }
    }
    #[test]
    fn existing_product_bond_skips_original_edge_lookup_and_all_map_property_reads() {
        let input = malformed_edge();
        let mut product = product(2);
        product.add_bond(0, 1, BondOrder::Single, None).unwrap();
        product.topology.atoms[0]
            .set_prop("old_mapno", PropertyValue::String("bad".into()))
            .unwrap();
        let before = product.topology.clone();
        run(
            &input,
            &mut product,
            &mut [false; 2],
            &mut vec![],
            &mut both_mapped(),
        )
        .unwrap();
        assert_eq!(product.topology, before);
    }
    #[test]
    fn intentional_template_bond_break_uses_signed_sort_then_unsigned_key_and_defers_original_read()
    {
        let input = malformed_edge();
        let mut product = product(2);
        product.topology.atoms[0]
            .set_prop("old_mapno", PropertyValue::Int(2))
            .unwrap();
        product.topology.atoms[1]
            .set_prop("old_mapno", PropertyValue::Int(-2))
            .unwrap();
        let mut mapping = both_mapped();
        mapping.template_bonds.insert(((-2i32) as u32, 2));
        run(
            &input,
            &mut product,
            &mut [false; 2],
            &mut vec![],
            &mut mapping,
        )
        .unwrap();
        assert!(product.topology.bonds.is_empty());
    }
    #[test]
    fn mapped_neighbors_not_bonded_in_template_get_aligned_missing_bond() {
        let input = input(2, &[(0, 1)]);
        let mut product = product(2);
        for (row, map) in [2, -2].into_iter().enumerate() {
            product.topology.atoms[row]
                .set_prop("old_mapno", PropertyValue::Int(map))
                .unwrap();
        }
        run(
            &input,
            &mut product,
            &mut [false; 2],
            &mut vec![],
            &mut both_mapped(),
        )
        .unwrap();
        assert_eq!(product.topology.bonds.len(), 1);
        assert_eq!(
            (
                product.topology.bonds[0].begin().index(),
                product.topology.bonds[0].end().index()
            ),
            (0, 1)
        );
    }
    #[test]
    fn both_property_presence_checks_precede_either_numeric_conversion() {
        let input = input(2, &[(0, 1)]);
        let mut product = product(2);
        product.topology.atoms[0]
            .set_prop("old_mapno", PropertyValue::String("bad".into()))
            .unwrap();
        assert!(matches!(
            run(
                &input,
                &mut product,
                &mut [false; 2],
                &mut vec![],
                &mut both_mapped()
            ),
            Err(ReactionProductError::Invariant {
                detail: "atoms should be mapped in product",
                ..
            })
        ));
        product.topology.atoms[1]
            .set_prop("old_mapno", PropertyValue::Int(2))
            .unwrap();
        assert!(matches!(
            run(
                &input,
                &mut product,
                &mut [false; 2],
                &mut vec![],
                &mut both_mapped()
            ),
            Err(ReactionProductError::PropertyInt {
                key: "old_mapno",
                ..
            })
        ));
        assert!(product.topology.bonds.is_empty());
    }
    #[test]
    fn target_operator_index_inserts_missing_key_after_begin_first_element_access() {
        let input = input(2, &[(0, 1)]);
        for empty_begin in [false, true] {
            let mut product = product(2);
            let mut mapping = mapping(2, &[(0, if empty_begin { vec![] } else { vec![0] })]);
            mapping.mapped.fill(true);
            let mut visited = [false; 2];
            let error = run(
                &input,
                &mut product,
                &mut visited,
                &mut vec![],
                &mut mapping,
            )
            .unwrap_err();
            assert!(
                matches!(error,ReactionProductError::Invariant {detail,..} if detail==if empty_begin {"missing mapped begin row"}else{"missing mapped end row"})
            );
            assert_eq!(mapping.reactant_to_product.contains_key(&1), !empty_begin);
            assert_eq!(visited, [true, false]);
        }
    }
    #[test]
    fn mapped_neighbor_invariant_precedes_original_edge_access_when_current_is_unmapped() {
        let input = malformed_edge();
        let mut product = product(1);
        let mut mapping = mapping(2, &[(0, vec![0])]);
        mapping.mapped[1] = true;
        assert!(matches!(
            run(
                &input,
                &mut product,
                &mut [false; 2],
                &mut vec![],
                &mut mapping
            ),
            Err(ReactionProductError::Invariant {
                detail: "reactant atom not present in product",
                ..
            })
        ));
        assert!(!mapping.reactant_to_product.contains_key(&1));
        assert!(product.topology.bonds.is_empty());
    }
    #[test]
    fn already_added_unmapped_neighbor_adds_bond_without_atom_or_queue_chiral_insertions() {
        let mut input = input(2, &[(0, 1)]);
        input.atoms[1].set_chiral_tag(ChiralTag::Allene);
        let mut product = product(2);
        let mut mapping = mapping(2, &[(0, vec![0]), (1, vec![1])]);
        let mut visited = [false; 2];
        let mut chiral = vec![];
        run(
            &input,
            &mut product,
            &mut visited,
            &mut chiral,
            &mut mapping,
        )
        .unwrap();
        assert_eq!(product.topology.atoms.len(), 2);
        assert_eq!(product.topology.bonds.len(), 1);
        assert_eq!(visited, [true, false]);
        assert!(chiral.is_empty());
    }
    #[test]
    fn empty_product_snapshot_still_enqueues_neighbor_and_chiral_check_then_hits_stack_invariant() {
        let mut input = input(2, &[(0, 1)]);
        input.atoms[1].set_chiral_tag(ChiralTag::Allene);
        let mut product = product(0);
        let mut mapping = mapping(2, &[(0, vec![])]);
        let mut visited = [false; 2];
        let mut chiral = vec![];
        assert!(matches!(
            run(
                &input,
                &mut product,
                &mut visited,
                &mut chiral,
                &mut mapping
            ),
            Err(ReactionProductError::Invariant {
                detail: "traversal atom not present in product",
                reactant_atom: Some(1),
                ..
            })
        ));
        assert!(product.topology.atoms.is_empty());
        assert_eq!(chiral, [1]);
        assert_eq!(visited, [true, false]);
        assert!(!mapping.reactant_to_product.contains_key(&1));
    }
    #[test]
    fn one_to_many_current_snapshot_preserves_physical_copy_order_across_new_mapping_updates() {
        let input = input(3, &[(0, 1), (1, 2)]);
        let mut product = product(2);
        let mut mapping = mapping(3, &[(0, vec![1, 0])]);
        mapping.mapped[0] = true;
        run(
            &input,
            &mut product,
            &mut [false; 3],
            &mut vec![],
            &mut mapping,
        )
        .unwrap();
        assert_eq!(mapping.reactant_to_product.get(&0), Some(&vec![1, 0]));
        assert_eq!(mapping.reactant_to_product.get(&1), Some(&vec![2, 3]));
        assert_eq!(mapping.reactant_to_product.get(&2), Some(&vec![4, 5]));
        assert_eq!(
            product
                .topology
                .bonds
                .iter()
                .map(|b| (b.begin().index(), b.end().index()))
                .collect::<Vec<_>>(),
            [(1, 2), (0, 3), (2, 4), (3, 5)]
        );
    }
    #[test]
    fn cycle_adds_cross_bond_once_using_already_added_unvisited_neighbor_case() {
        let input = input(3, &[(0, 1), (0, 2), (1, 2)]);
        let mut product = product(1);
        let mut mapping = mapping(3, &[(0, vec![0])]);
        mapping.mapped[0] = true;
        let mut visited = [false; 3];
        run(
            &input,
            &mut product,
            &mut visited,
            &mut vec![],
            &mut mapping,
        )
        .unwrap();
        assert_eq!(visited, [true; 3]);
        assert_eq!(product.topology.atoms.len(), 3);
        assert_eq!(product.topology.bonds.len(), 3);
    }
    #[test]
    fn later_copy_failure_keeps_prior_atom_mapping_bond_prefix_but_does_not_enqueue_neighbor() {
        let mut input = input(2, &[(0, 1)]);
        input.atoms[1].set_chiral_tag(ChiralTag::Allene);
        let mut product = product(2);
        let mut mapping = mapping(2, &[(0, vec![0, 99])]);
        mapping.mapped[0] = true;
        let mut visited = [false; 2];
        let mut chiral = vec![];
        assert!(
            matches!(run(&input,&mut product,&mut visited,&mut chiral,&mut mapping),Err(ReactionProductError::TopologyEdit(cosmolkit_model::TopologyEditError::AtomOutOfRange {atom,..})) if atom==AtomId::new(99))
        );
        assert_eq!(product.topology.atoms.len(), 4);
        assert_eq!(product.topology.bonds.len(), 1);
        assert_eq!(mapping.reactant_to_product.get(&1), Some(&vec![2, 3]));
        assert_eq!(visited, [true, false]);
        assert!(chiral.is_empty());
    }
}
