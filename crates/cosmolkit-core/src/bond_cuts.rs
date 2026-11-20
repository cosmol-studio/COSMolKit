//! Selected-bond fragmentation against the original bond index space.

use std::collections::HashSet;

use cosmolkit_model::{
    AtomId, AtomSpec, Bond, BondDirection, BondId, BondSpec, ChiralTag, CoordinateBlock, Element,
    MoleculeProperties, PropertyValue, PropertyValueKind, TopologyBatchEdit, TopologyBlock,
    TopologyMapping,
};

/// Labels identify the atom to which each dummy will attach.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BondCut {
    bond: BondId,
    labels: CutLabels,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum CutLabels {
    OppositeIndices,
    AttachmentSides([usize; 2]),
    MatchedEnds([(AtomId, usize); 2]),
}

impl BondCut {
    pub const fn at_bond(
        bond: BondId,
        begin_label: Option<usize>,
        end_label: Option<usize>,
    ) -> Result<Self, BondCutError> {
        let labels = match (begin_label, end_label) {
            (None, None) => CutLabels::OppositeIndices,
            (Some(begin), Some(end)) => CutLabels::AttachmentSides([begin, end]),
            _ => return Err(BondCutError::IncompleteLabels { bond }),
        };
        Ok(Self { bond, labels })
    }

    /// Keep each match label paired with its source atom, regardless of bond direction.
    #[doc(hidden)]
    pub const fn from_matched_ends(bond: BondId, ends: [(AtomId, usize); 2]) -> Self {
        Self {
            bond,
            labels: CutLabels::MatchedEnds(ends),
        }
    }
}

#[derive(Clone, Debug, PartialEq, thiserror::Error)]
pub enum BondCutError {
    #[error("provide both attachment labels for bond {bond}, or omit both")]
    IncompleteLabels { bond: BondId },
    #[error("select at least one bond to cut")]
    EmptySelection,
    #[error("bond {bond} occurs more than once in the selection")]
    Duplicate { bond: BondId },
    #[error("selected bond {bond} is outside the input's {count} bonds")]
    BondOutOfRange { bond: BondId, count: usize },
    #[error("selected bond {bond} joins {actual:?}; expected {expected:?}")]
    EndpointsChanged {
        bond: BondId,
        expected: [AtomId; 2],
        actual: [AtomId; 2],
    },
    #[error("dummy isotope {label} on the atom {attachment} side of bond {bond} exceeds 65535")]
    IsotopeOutOfRange {
        bond: BondId,
        attachment: AtomId,
        label: usize,
    },
    #[error(
        "atom {atom} has a pre-existing _newBondOrder of type {kind:?}, expected an incident-bond list"
    )]
    IncidentOrderProperty {
        atom: AtomId,
        kind: PropertyValueKind,
    },
    #[error(transparent)]
    Edit(#[from] cosmolkit_model::TopologyEditError),
    #[error(transparent)]
    Coordinates(#[from] cosmolkit_model::CoordinateValidationError),
    #[error(transparent)]
    Stereo(#[from] crate::StereoOrderError),
    #[error(transparent)]
    BondProperty(#[from] cosmolkit_model::BondValueError),
    #[error(transparent)]
    AtomProperty(#[from] cosmolkit_model::AtomPropertyError),
    #[error(transparent)]
    NumericProperty(#[from] crate::PropertyUIntReadError),
    #[error(transparent)]
    TextProperty(#[from] crate::PropertyStringError),
    #[error(transparent)]
    Valence(#[from] cosmolkit_model::ValenceError),
}

struct ResolvedBondCut<'a> {
    bond: &'a Bond,
    attachment_labels: [u16; 2],
}

fn resolve_cuts<'a>(
    bonds: &'a [Bond],
    cuts: &[BondCut],
) -> Result<Vec<ResolvedBondCut<'a>>, BondCutError> {
    // RDKit✔️✔️:   if (!bondIndices.get()) {
    // RDKit✔️✔️:     throw_value_error("empty bond indices");
    // RDKit✔️✔️:   }
    if cuts.is_empty() {
        return Err(BondCutError::EmptySelection);
    }

    // RDKit✔️✔️:   if (std::unordered_set<unsigned int>(bondIndices.begin(), bondIndices.end())
    // RDKit✔️✔️:           .size() < bondIndices.size()) {
    // RDKit✔️✔️:     throw ValueErrorException("bondIndices contains duplicates");
    // RDKit✔️✔️:   }
    let mut seen = HashSet::with_capacity(cuts.len());
    for cut in cuts {
        if !seen.insert(cut.bond) {
            return Err(BondCutError::Duplicate { bond: cut.bond });
        }
    }

    cuts.iter()
        .map(|cut| {
            // RDKit✔️✔️:     bondsToRemove.push_back(res->getBondWithIdx(bondIdx));
            let bond = bonds
                .get(cut.bond.index())
                .ok_or(BondCutError::BondOutOfRange {
                    bond: cut.bond,
                    count: bonds.len(),
                })?;
            let actual = [bond.begin(), bond.end()];
            let labels = match cut.labels {
                CutLabels::OppositeIndices => [bond.end().index(), bond.begin().index()],
                CutLabels::AttachmentSides(pair) => pair,
                CutLabels::MatchedEnds([(a, x), (b, y)]) if [a, b] == actual => [x, y],
                CutLabels::MatchedEnds([(a, x), (b, y)]) if [b, a] == actual => [y, x],
                CutLabels::MatchedEnds([(a, _), (b, _)]) => {
                    return Err(BondCutError::EndpointsChanged {
                        bond: cut.bond,
                        expected: [a, b],
                        actual,
                    });
                }
            };

            // RDKit❗✔️: void Atom::setIsotope(unsigned int what) { d_isotope = what; }
            // RDKit❗✔️:   std::uint16_t d_isotope;
            // COSMolKit rejects unrepresentable labels instead of wrapping to a different isotope.
            let isotope = |side: usize| {
                u16::try_from(labels[side]).map_err(|_| BondCutError::IsotopeOutOfRange {
                    bond: cut.bond,
                    attachment: actual[side],
                    label: labels[side],
                })
            };
            Ok(ResolvedBondCut {
                bond,
                attachment_labels: [isotope(0)?, isotope(1)?],
            })
        })
        .collect()
}

pub struct BondCutResult {
    pub topology: TopologyBlock,
    pub coordinates: CoordinateBlock,
    pub properties: MoleculeProperties,
    pub mapping: TopologyMapping,
}

fn adjust_tetrahedral_order(
    source: &TopologyBlock,
    edit: &mut TopologyBatchEdit,
    atom: AtomId,
    cut: BondId,
) -> Result<(), BondCutError> {
    let original = &source.atoms[atom.index()];
    let tag = original.chiral_tag();
    if !matches!(tag, ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw) {
        return Ok(());
    }

    // RDKit✔️✔️:   const int check_bond_index = static_cast<int>(bond->getIdx());
    // RDKit✔️✔️:   if (nAt->getPropIfPresent(newBondOrder, incomingOrder)) {
    // RDKit✔️✔️:     for (int bidx : incomingOrder) {
    // RDKit✔️✔️:       if (bidx != check_bond_index) {
    // RDKit✔️✔️:         newOrder.push_back(bidx);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     for (auto obond : mol.atomBonds(oAt)) {
    // RDKit✔️✔️:       if (obond != bond) {
    // RDKit✔️✔️:         newOrder.push_back(obond->getIdx());
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   newOrder.push_back(bond->getIdx());
    // RDKit❗✔️:   nAt->setProp(newBondOrder, newOrder, true);
    // The private native std::list<int> survives successive cuts, then is cleared at commit.
    // A caller-supplied vector is not the same RDProps type as a list.
    let incident = source.adjacency.neighbors_of(atom.index());
    let signed = |id: BondId| {
        u32::try_from(id.index()).map(|x| x as i32).map_err(|_| {
            crate::StereoOrderError::BondIndexSourceWidth {
                bond_index: id.index(),
            }
        })
    };
    if let Some(value) = original.prop("_newBondOrder") {
        return Err(BondCutError::IncidentOrderProperty {
            atom,
            kind: value.kind(),
        });
    }
    let target = edit.atom_mut(atom)?;
    let mut order = match target.prop("_newBondOrder") {
        Some(PropertyValue::IntVector(previous)) => previous.clone(),
        Some(other) => {
            return Err(BondCutError::IncidentOrderProperty {
                atom,
                kind: other.kind(),
            });
        }
        None => incident
            .iter()
            .map(|n| signed(n.bond))
            .collect::<Result<Vec<_>, _>>()?,
    };
    let id = signed(cut)?;
    order.retain(|&i| i != id);
    order.push(id);

    // RDKit✔️✔️:   unsigned int nSwaps = oAt->getPerturbationOrder(newOrder);
    // RDKit✔️✔️:   nAt->setChiralTag(oAt->getChiralTag());
    // RDKit✔️✔️:   if (nSwaps % 2) {
    // RDKit✔️✔️:     nAt->invertChirality();
    // RDKit✔️✔️:   }
    let swaps = crate::atom_perturbation_order(&order, incident.iter().map(|n| n.bond.index()))?;
    target.set_computed_prop("_newBondOrder", PropertyValue::IntVector(order))?;
    target.set_chiral_tag(tag);
    if swaps % 2 != 0 {
        crate::invert_atom_chirality(target)?;
    }
    Ok(())
}

/// Produce one detached molecular graph; a ring cut need not split it into components.
pub fn cut_selected_bonds(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    properties: &MoleculeProperties,
    cuts: &[BondCut],
) -> Result<BondCutResult, BondCutError> {
    let cuts = resolve_cuts(&topology.bonds, cuts)?;
    for conf in &coordinates.conformers_2d {
        conf.validate_for_atom_count(topology.atoms.len())?;
    }
    for conf in &coordinates.conformers_3d {
        conf.validate_for_atom_count(topology.atoms.len())?;
    }

    // RDKit✔️✔️:   std::unique_ptr<RWMol> res(new RWMol(mol));
    // RDKit✔️✔️:   res->beginBatchEdit();
    let mut edit = topology.begin_batch_edit()?;
    let mut coordinates = coordinates.clone();
    let mut properties = properties.clone();
    let mut touched = Vec::new();
    let mut seen = HashSet::new();

    for cut in cuts {
        let bond = cut.bond;
        // RDKit✔️✔️:     unsigned int bidx = bond->getBeginAtomIdx();
        // RDKit✔️✔️:     unsigned int eidx = bond->getEndAtomIdx();
        // RDKit✔️✔️:     Bond::BondType bT = bond->getBondType();
        // RDKit✔️✔️:     Bond::BondDir bD = bond->getBondDir();
        // RDKit✔️✔️:     res->removeBond(bidx, eidx);
        let [begin, end] = [bond.begin(), bond.end()];
        edit.remove_bond(bond.id())?;

        // RDKit✔️✔️:       std::unique_ptr<Atom> at1(new Atom(0));
        // RDKit✔️✔️:       std::unique_ptr<Atom> at2(new Atom(0));
        // RDKit❗✔️:       if (dummyLabels) {
        // RDKit❗✔️:         at1->setIsotope((*dummyLabels)[i].first);
        // RDKit❗✔️:         at2->setIsotope((*dummyLabels)[i].second);
        // RDKit❗✔️:       } else {
        // RDKit❗✔️:         at1->setIsotope(bidx);
        // RDKit❗✔️:         at2->setIsotope(eidx);
        // RDKit❗✔️:       }
        // RDKit✔️✔️:       unsigned int idx1 = res->addAtom(at1.release(), false, true);
        // RDKit✔️✔️:       unsigned int idx2 = res->addAtom(at2.release(), false, true);
        let [on_begin, on_end] = cut.attachment_labels;
        let first = edit.add_atom(AtomSpec::new(Element::DUMMY).with_isotope(on_end));
        let second = edit.add_atom(AtomSpec::new(Element::DUMMY).with_isotope(on_begin));

        // RDKit✔️✔️:       bondidx = res->addBond(idx1, eidx, bT) - 1;
        // RDKit✔️✔️:       if (bD == Bond::ENDDOWNRIGHT || bD == Bond::ENDUPRIGHT) {
        // RDKit✔️✔️:         res->getBondWithIdx(bondidx)->setBondDir(bD);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       bondidx = res->addBond(bidx, idx2, bT) - 1;
        // RDKit✔️✔️:       res->getBondWithIdx(bondidx)->setBondDir(bD);
        let direction = bond.direction();
        let first_dir = match direction {
            BondDirection::EndDownRight | BondDirection::EndUpRight => direction,
            _ => BondDirection::None,
        };
        for (a, b, d) in [(first, end, first_dir), (begin, second, direction)] {
            let id = BondId::new(edit.add_bond_by_order(a, b, bond.order())? - 1);
            let new = edit.bond_mut(id)?;
            new.set_direction(d);
            // RDKit✔️✔️:       if (outQuery) {
            // RDKit✔️✔️:         auto qb = std::make_unique<QueryBond>(*res->getBondWithIdx(bondidx));
            // RDKit✔️✔️:         qb->setQuery(outQuery->copy());
            // RDKit✔️✔️:         res->replaceBond(bondidx, qb.get());
            // RDKit✔️✔️:       }
            if let Some(query) = bond.query() {
                let spec = BondSpec::new(a, b, bond.order())
                    .with_aromatic(new.is_aromatic())
                    .with_direction(d)
                    .with_query(query.clone());
                *new = Bond::from_spec(id, spec);
            }
        }

        // RDKit✔️✔️:   for (const auto *atom : {bgn, end}) {
        // RDKit✔️✔️:     for (auto obnd : mol.atomBonds(atom)) {
        // RDKit✔️✔️:       if (obnd->getIdx() != bnd->getIdx() && !obnd->getStereoAtoms().empty()) {
        // RDKit✔️✔️:         obnd->setProp(molfragSaveStereo, obnd->getStereo());
        // RDKit✔️✔️:         res.emplace_back(obnd, obnd->getStereoAtoms());
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       for (auto &stereo_atoms : nbr_bond_stereo) {
        // RDKit✔️✔️:         std::replace(stereo_atoms.second.begin(), stereo_atoms.second.end(),
        // RDKit✔️✔️:                      bidx, idx1);
        // RDKit✔️✔️:         std::replace(stereo_atoms.second.begin(), stereo_atoms.second.end(),
        // RDKit✔️✔️:                      eidx, idx2);
        // RDKit✔️✔️:         stereo_atoms.first->getStereoAtoms().swap(stereo_atoms.second);
        // RDKit✔️✔️:         stereo_atoms.first->clearProp(molfragSaveStereo);
        // Batch removal only sets a mask until finish. New dummy bonds have no stereo references.
        for atom in [begin, end] {
            for neighbor in topology.adjacency.neighbors_of(atom.index()) {
                if neighbor.bond == bond.id() {
                    continue;
                }
                let adjacent = edit.bond_mut(neighbor.bond)?;
                if adjacent.stereo_atom_references().is_empty() {
                    continue;
                }
                let refs = adjacent
                    .stereo_atom_references()
                    .iter()
                    .map(|&x| {
                        if x == begin {
                            first
                        } else if x == end {
                            second
                        } else {
                            x
                        }
                    })
                    .collect();
                adjacent.set_source_stereo_atom_references(refs);
                adjacent.clear_prop("_molfragSaveStereo")?;
            }
            adjust_tetrahedral_order(topology, &mut edit, atom, bond.id())?;
        }

        // RDKit✔️✔️:         conf->setAtomPos(idx1, conf->getAtomPos(bidx));
        // RDKit✔️✔️:         conf->setAtomPos(idx2, conf->getAtomPos(eidx));
        for (dummy, original) in [(first, begin), (second, end)] {
            for conf in &mut coordinates.conformers_2d {
                conf.source_set_atom_position(dummy.index(), conf.coordinates()[original.index()])?;
            }
            for conf in &mut coordinates.conformers_3d {
                conf.source_set_atom_position(dummy.index(), conf.coordinates()[original.index()])?;
            }
        }
        for id in [first, second, begin, end] {
            if seen.insert(id) {
                touched.push(id);
            }
        }
    }

    // RDKit❗❌:   res->commitBatchEdit();
    // RDKit❗✔️:   res->clearComputedProps();
    // The shared dense-row batch kernel compacts bonds and clears computed properties.
    // It can scan/shift surviving rows for each removed bond, unlike native edge deletion.
    let (mut topology, mapping) = edit.finish_source::<BondCutError>(
        &mut coordinates,
        &mut properties,
        &mut |v| Ok(crate::property_value_to_uint(v)?),
        &mut |v| Ok(crate::property_value_to_string(v)?),
    )?;
    properties.remap_topology(mapping.atoms().new_to_old(), mapping.bonds().new_to_old());

    // RDKit✔️✔️:   for (auto [idx, count] : atomsToUpdate) {
    // RDKit✔️✔️:     Atom *atom = res->getAtomWithIdx(idx);
    // RDKit✔️✔️:     atom->updatePropertyCache(false);
    for atom in touched {
        crate::valence::update_source_atom_cache(&mut topology, atom, false)?;
    }
    Ok(BondCutResult {
        topology,
        coordinates,
        properties,
        mapping,
    })
}
