//! Source-indexed transformations of detached molecule blocks.
use std::{
    borrow::Cow,
    collections::{BTreeMap, BTreeSet},
};

use crate::{property_value_to_string, property_value_to_uint};
use cosmolkit_model::{
    Atom, AtomId, AtomMapping, AtomSpec, Bond, BondId, BondMapping, BondSpec, Conformer2D,
    Conformer3D, CoordinateBlock, MoleculeProperties, PropertyValue, SubstanceGroupId,
    TopologyBlock, TopologyEditError, TopologyMapping, TopologyValidationError,
};

pub type IndexedEdit = (
    TopologyBlock,
    CoordinateBlock,
    MoleculeProperties,
    TopologyMapping,
);

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum IndexedGraphError {
    #[error("the atom permutation must not be empty")]
    EmptyAtomOrder,
    #[error("atom {atom} has an invalid indexed property {property}")]
    InvalidIndexProperty {
        atom: AtomId,
        property: &'static str,
    },
    #[error("bond {bond} has an invalid stereo neighbor {atom}")]
    InvalidStereoNeighbor { bond: BondId, atom: AtomId },
    #[error("cannot match {dimension}D coordinate set {id}: {reason}")]
    CoordinateKey {
        dimension: u8,
        id: usize,
        reason: &'static str,
    },
    #[error("{0}")]
    TopologyEdit(#[from] TopologyEditError),
    #[error("{0}")]
    Topology(#[from] TopologyValidationError),
    #[error("{0}")]
    Subset(#[from] crate::AtomPathSubsetCopyError),
    #[error("{0}")]
    NumericProperty(#[from] crate::PropertyUIntReadError),
    #[error("{0}")]
    TextProperty(#[from] crate::PropertyStringError),
    #[error("{0}")]
    Valence(#[from] cosmolkit_model::ValenceError),
}

fn remap_ring_stereo(
    atom: &mut Atom,
    old_to_new: &[Option<AtomId>],
) -> Result<(), IndexedGraphError> {
    // BEGIN RDKIT CPP FUNCTION Renumber.cpp :: renumberAtoms
    // RDKit❗✔️:       for (auto &val : nAtoms) {
    // RDKit❗✔️:         if (val < 0) {
    // RDKit❗✔️:           val = -1 * (revOrder[(-val - 1)] + 1);
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           val = revOrder[val - 1] + 1;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:       nAtom->setProp(common_properties::_ringStereoAtoms, nAtoms, true);
    // RDKit❗✔️:       otherAtom = revOrder[otherAtom];
    // RDKit❗✔️:       nAtom->setProp(common_properties::_ringStereoOtherAtom, otherAtom, true);
    // Reject unrepresentable or stale references instead of indexing outside the source graph.
    let id = atom.id();
    let bad = |property| IndexedGraphError::InvalidIndexProperty { atom: id, property };
    let mapped = |i: usize| old_to_new.get(i).copied().flatten().map(AtomId::index);
    if let Some(value) = atom.prop("_ringStereoAtoms") {
        let PropertyValue::IntVector(values) = value else {
            return Err(bad("_ringStereoAtoms"));
        };
        let values = values
            .iter()
            .map(|v| {
                let old = (v.unsigned_abs() as usize).checked_sub(1);
                let new = old
                    .and_then(mapped)
                    .and_then(|n| i32::try_from(n + 1).ok())
                    .ok_or_else(|| bad("_ringStereoAtoms"))?;
                Ok(if *v < 0 { -new } else { new })
            })
            .collect::<Result<Vec<_>, IndexedGraphError>>()?;
        atom.set_computed_prop("_ringStereoAtoms", PropertyValue::IntVector(values))
            .map_err(TopologyEditError::from)?;
    }
    if let Some(value) = atom.prop("_ringStereoOtherAtom") {
        let new = mapped(property_value_to_uint(value)? as usize)
            .and_then(|n| u32::try_from(n).ok())
            .ok_or_else(|| bad("_ringStereoOtherAtom"))?;
        atom.set_computed_prop("_ringStereoOtherAtom", new)
            .map_err(TopologyEditError::from)?;
    }
    // END RDKIT CPP FUNCTION renumberAtoms
    Ok(())
}

pub fn renumber_atoms(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    order: &[AtomId],
) -> Result<IndexedEdit, IndexedGraphError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::renumberAtoms (Release_2026_09_1)
    // RDKit❗✔️: ROMol *renumberAtoms(const ROMol &mol,
    // RDKit❗✔️:                      const std::vector<unsigned int> &newOrder) {
    // RDKit❗✔️:   unsigned int nAts = mol.getNumAtoms();
    // RDKit❗✔️:   PRECONDITION(newOrder.size() == nAts, "bad newOrder size");
    //
    // RDKit❗✔️:   std::vector<unsigned int> revOrder(nAts, nAts);
    // RDKit❗✔️:   for (unsigned int nIdx = 0; nIdx < nAts; ++nIdx) {
    // RDKit❗✔️:     unsigned int oIdx = newOrder[nIdx];
    // RDKit❗✔️:     if (oIdx >= nAts) {
    // RDKit❗✔️:       throw ValueErrorException("idx value exceeds numAtoms");
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (revOrder[oIdx] != nAts) {
    // RDKit❗✔️:       throw ValueErrorException("newOrder contains duplicate indices");
    // RDKit❗✔️:     }
    // RDKit❗✔️:     revOrder[oIdx] = nIdx;
    // RDKit❗✔️:   }
    //
    // RDKit❗✔️:   // ------
    // RDKit❗✔️:   // newOrder[i] : which atom should be in position i of the new mol
    // RDKit❗✔️:   // revOrder[i] : where atom i of the original mol landed in the new mol
    // RDKit❗✔️:   auto *res = new RWMol();
    //
    // RDKit❗✔️:   // copy over the atoms:
    // RDKit❗✔️:   for (unsigned int nIdx = 0; nIdx < nAts; ++nIdx) {
    // RDKit❗✔️:     unsigned int oIdx = newOrder[nIdx];
    // RDKit❗✔️:     const Atom *oAtom = mol.getAtomWithIdx(oIdx);
    // RDKit❗✔️:     Atom *nAtom = oAtom->copy();
    // RDKit❗✔️:     res->addAtom(nAtom, false, true);
    // reordered_atoms validates the inverse map and copies typed topology. Empty input
    // follows the public Python precondition. Ring and valence caches are invalidated.
    if order.is_empty() {
        return Err(IndexedGraphError::EmptyAtomOrder);
    }
    let (mut t, map) = topology.reordered_atoms(order)?;
    for atom in &mut t.atoms {
        remap_ring_stereo(atom, &map.atoms.old_to_new)?;
    }
    // RDKit❗✔️:   // now the bonds:
    // RDKit❗✔️:   for (const auto oBond : mol.bonds()) {
    // RDKit❗✔️:     Bond *nBond = oBond->copy();
    // RDKit❗✔️:     nBond->setBeginAtomIdx(revOrder[oBond->getBeginAtomIdx()]);
    // RDKit❗✔️:     nBond->setEndAtomIdx(revOrder[oBond->getEndAtomIdx()]);
    // RDKit❗✔️:     res->addBond(nBond, true);
    // RDKit❗✔️:     // take care of atom-numbering-dependent properties:
    // RDKit❗✔️:     for (auto &idx : nBond->getStereoAtoms()) {
    // RDKit❗✔️:       idx = revOrder[idx];
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    //
    // RDKit❗✔️:   // Conformers:
    // RDKit❗✔️:   for (auto oConf = mol.beginConformers(); oConf != mol.endConformers();
    // RDKit❗✔️:        ++oConf) {
    // RDKit❗✔️:     auto *nConf = new Conformer(nAts);
    // RDKit❗✔️:     for (unsigned int i = 0; i < nAts; ++i) {
    // RDKit❗✔️:       nConf->setAtomPos(i, (*oConf)->getAtomPos(newOrder[i]));
    // RDKit❗✔️:     }
    // RDKit❗✔️:     nConf->setId((*oConf)->getId());
    // RDKit❗✔️:     nConf->set3D((*oConf)->is3D());
    // RDKit❗✔️:     res->addConformer(nConf);
    // RDKit❗✔️:   }
    for (old, new) in topology.bonds.iter().zip(&mut t.bonds) {
        new.set_source_stereo_atom_references(
            old.stereo_atom_references()
                .iter()
                .map(|a| map.atoms.old_to_new[a.index()].expect("complete permutation"))
                .collect(),
        );
    }
    // RDKit❗✔️:   if (mol.getStereoGroups().size()) {
    // RDKit❗✔️:     std::vector<StereoGroup> nsgs;
    // RDKit❗✔️:     nsgs.reserve(mol.getStereoGroups().size());
    // RDKit❗✔️:     for (const auto &osg : mol.getStereoGroups()) {
    // RDKit❗✔️:       std::vector<Atom *> ats;
    // RDKit❗✔️:       std::vector<Bond *> bds;
    // RDKit❗✔️:       ats.reserve(osg.getAtoms().size());
    // RDKit❗✔️:       bds.reserve(osg.getBonds().size());
    // RDKit❗✔️:       for (const auto aptr : osg.getAtoms()) {
    // RDKit❗✔️:         ats.push_back(res->getAtomWithIdx(revOrder[aptr->getIdx()]));
    // RDKit❗✔️:       }
    // RDKit❗✔️:       for (const auto bptr : osg.getBonds()) {
    // RDKit❗✔️:         bds.push_back(
    // RDKit❗✔️:             res->getBondWithIdx(bptr->getIdx()));  // bonds do not change order
    // RDKit❗✔️:       }
    //
    // RDKit❗✔️:       nsgs.emplace_back(osg.getGroupType(), ats, bds, osg.getReadId());
    // RDKit❗✔️:     }
    // RDKit❗✔️:     res->setStereoGroups(std::move(nsgs));
    // RDKit❗✔️:   }
    //
    // RDKit❗✔️:   return dynamic_cast<ROMol *>(res);
    // RDKit creates fresh molecule/conformer records and enhanced groups with read IDs only.
    // END RDKIT CPP FUNCTION MolOps::renumberAtoms
    t.substance_groups.clear();
    for group in &mut t.stereo_groups {
        cosmolkit_model::set_stereo_group_write_id(group, 0);
    }
    let rows = map.retained_atom_indices();
    let coords = CoordinateBlock {
        conformers_2d: coordinates
            .conformers_2d
            .iter()
            .map(|c| Conformer2D::new(c.id(), rows.iter().map(|&i| c.coordinates()[i]).collect()))
            .collect(),
        conformers_3d: coordinates
            .conformers_3d
            .iter()
            .map(|c| {
                Conformer3D::new(
                    c.id(),
                    rows.iter().map(|&i| c.coordinates()[i]).collect(),
                    c.is_3d(),
                )
            })
            .collect(),
        source_coordinate_dim: coordinates.source_coordinate_dim,
        source_conformer_order: coordinates.source_conformer_order.clone(),
    };
    Ok((t, coords, MoleculeProperties::default(), map))
}

fn compact_rows<I: Copy>(
    keep: impl Iterator<Item = bool>,
    id: impl Fn(usize) -> I,
) -> (Vec<Option<I>>, Vec<Option<I>>) {
    let mut retained = Vec::new();
    let forward = keep
        .enumerate()
        .map(|(old, keep)| {
            keep.then(|| {
                let new = id(retained.len());
                retained.push(Some(id(old)));
                new
            })
        })
        .collect();
    (forward, retained)
}

pub fn remove_atom(
    t: &TopologyBlock,
    c: &CoordinateBlock,
    p: &MoleculeProperties,
    atom: AtomId,
) -> Result<IndexedEdit, IndexedGraphError> {
    if atom.index() >= t.atoms.len() {
        return Err(TopologyEditError::AtomOutOfRange {
            atom,
            atom_count: t.atoms.len(),
        }
        .into());
    }
    let (old_to_new, new_to_old) =
        compact_rows(t.atoms.iter().map(|a| a.id() != atom), AtomId::new);
    let atoms = AtomMapping {
        old_to_new,
        new_to_old,
    };
    let (old_to_new, new_to_old) = compact_rows(
        t.bonds.iter().map(|b| b.begin() != atom && b.end() != atom),
        BondId::new,
    );
    let map = TopologyMapping {
        atoms,
        bonds: BondMapping {
            old_to_new,
            new_to_old,
        },
    };
    let (mut t, mut c, mut p) = (t.clone(), c.clone(), p.clone());
    cosmolkit_model::remove_indexed_atom::<IndexedGraphError>(
        &mut t,
        atom,
        &mut c,
        &mut p,
        |v| Ok(property_value_to_uint(v)?),
        |v| Ok(property_value_to_string(v)?),
    )?;
    p.remap_topology(&map.atoms.new_to_old, &map.bonds.new_to_old);
    Ok((t, c, p, map))
}

pub fn prepare_atom_replacement(
    t: &TopologyBlock,
    target: AtomId,
    spec: &AtomSpec,
    preserve_props: bool,
) -> Result<Atom, IndexedGraphError> {
    let old = t
        .atoms
        .get(target.index())
        .ok_or(TopologyEditError::AtomOutOfRange {
            atom: target,
            atom_count: t.atoms.len(),
        })?;
    // BEGIN RDKIT CPP FUNCTION RWMol::replaceAtom
    // RDKit❗✔️:   if (preserveProps) {
    // RDKit❗✔️:     const bool replaceExistingData = false;
    // RDKit❗✔️:     atom_p->updateProps(*d_graph[vd], replaceExistingData);
    // RDKit❗✔️:   }
    let mut atom = Atom::from_spec(target, spec.clone());
    if preserve_props {
        atom.replace_source_properties_from(old);
    }
    // These are detached local-ID inputs, not native pointers. Validate the
    // replacement's references before a live runtime can apply the prepared row.
    if let Some(order) = atom.template_attachment_order() {
        order
            .validate_for_atom_count(t.atoms.len())
            .map_err(|source| TopologyValidationError::TemplateAttachmentOrder {
                atom: target,
                source,
            })?;
    }
    // END RDKIT CPP FUNCTION RWMol::replaceAtom
    Ok(atom)
}

pub fn replace_atom(topology: &mut Cow<'_, TopologyBlock>, atom: Atom) {
    // RDKit❗✔️:   d_graph[vd] = atom_p;
    // Row identity and endpoints are unchanged; no coordinate/property mapping
    // or clone of an unrelated block is required.
    let index = atom.id().index();
    topology.to_mut().atoms[index] = atom;
}

pub fn prepare_bond_replacement(
    t: &TopologyBlock,
    target: BondId,
    spec: &BondSpec,
    preserve_props: bool,
    keep_sgroups: bool,
) -> Result<PreparedBondReplacement, IndexedGraphError> {
    let old = t
        .bonds
        .get(target.index())
        .ok_or(TopologyEditError::BondOutOfRange {
            bond: target,
            bond_count: t.bonds.len(),
        })?;
    let replacement = Bond::from_spec(target, spec.clone());
    if matches!(
        replacement.stereo(),
        cosmolkit_model::BondStereo::Cis | cosmolkit_model::BondStereo::Trans
    ) && replacement.stereo_atoms().is_none()
    {
        return Err(TopologyValidationError::StereoAtomsRequired {
            bond: target,
            stereo: replacement.stereo(),
        }
        .into());
    }
    if let Some(refs) = replacement.stereo_atoms() {
        for ((side, opposite), atom) in [(old.begin(), old.end()), (old.end(), old.begin())]
            .into_iter()
            .zip(refs)
        {
            if atom == opposite
                || !t
                    .adjacency
                    .neighbors_of(side.index())
                    .iter()
                    .any(|n| n.atom_index == atom.index())
            {
                return Err(IndexedGraphError::InvalidStereoNeighbor { bond: target, atom });
            }
        }
    }
    // RDKit❗✔️:   auto orderDifference =
    // RDKit❗✔️:       bond_p->getBondTypeAsDouble() - obond->getBondTypeAsDouble();
    // RDKit❗✔️:   if (!keepSGroups) {
    // RDKit❗✔️:     removeSubstanceGroupsReferencingBond(*this, idx);
    // RDKit❗✔️:   }
    for &atom in replacement.stereo_atom_references() {
        if atom.index() >= t.atoms.len() {
            return Err(TopologyValidationError::StereoReferenceOutOfRange {
                bond: target,
                atom,
                atom_count: t.atoms.len(),
            }
            .into());
        }
    }
    // Perform fallible preparation on the affected groups only. A bad raw
    // PARENT value must not leave live topology partly edited.
    crate::bond_type_as_double(replacement.order())?;
    crate::bond_type_as_double(old.order())?;
    let groups = if !keep_sgroups && !t.substance_groups.is_empty() {
        let mut groups = t.substance_groups.clone();
        cosmolkit_model::remove_source_groups_referencing_bond::<IndexedGraphError>(
            &mut groups,
            target,
            &mut |v| Ok(property_value_to_uint(v)?),
        )?;
        Some(groups)
    } else {
        None
    };
    Ok(PreparedBondReplacement {
        target,
        replacement,
        preserve_props,
        groups,
    })
}

pub struct PreparedBondReplacement {
    target: BondId,
    replacement: Bond,
    preserve_props: bool,
    groups: Option<Vec<cosmolkit_model::SubstanceGroup>>,
}

pub fn replace_bond(
    topology: &mut Cow<'_, TopologyBlock>,
    prepared: PreparedBondReplacement,
) -> Result<(), IndexedGraphError> {
    let topology = topology.to_mut();
    // The same model kernel owns bond replacement. Group conversion and bond
    // order reads have succeeded before mutation; no second group pass is needed.
    cosmolkit_model::replace_source_bond::<IndexedGraphError>(
        topology,
        prepared.target,
        &prepared.replacement,
        prepared.preserve_props,
        true,
        |order| Ok(crate::bond_type_as_double(order)?),
        |v| Ok(property_value_to_uint(v)?),
    )?;
    if let Some(groups) = prepared.groups {
        topology.substance_groups = groups;
    }
    Ok(())
}

fn match_coordinates<'a, C>(
    left: &[C],
    right: &'a [C],
    dim: u8,
    id: impl Fn(&C) -> usize,
) -> Result<Vec<&'a C>, IndexedGraphError> {
    let err = |id, reason| IndexedGraphError::CoordinateKey {
        dimension: dim,
        id,
        reason,
    };
    let mut lookup = BTreeMap::new();
    for c in right {
        if lookup.insert(id(c), c).is_some() {
            return Err(err(id(c), "duplicate right ID"));
        }
    }
    let mut seen = BTreeSet::new();
    let mut matched = Vec::with_capacity(left.len());
    for c in left {
        let key = id(c);
        if !seen.insert(key) {
            return Err(err(key, "duplicate left ID"));
        }
        matched.push(
            lookup
                .remove(&key)
                .ok_or_else(|| err(key, "missing from right"))?,
        );
    }
    if let Some((&key, _)) = lookup.first_key_value() {
        return Err(err(key, "missing from left"));
    }
    Ok(matched)
}

enum AppendedCoordinates<'a> {
    Unchanged,
    Replace(&'a CoordinateBlock),
    Extend {
        two_d: Vec<&'a Conformer2D>,
        three_d: Vec<&'a Conformer3D>,
    },
}

pub struct PreparedAppend<'a> {
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    groups: Vec<cosmolkit_model::SubstanceGroup>,
    stereo: Option<Vec<cosmolkit_model::StereoGroup>>,
    coordinates: AppendedCoordinates<'a>,
    mapping: TopologyMapping,
}

pub fn prepare_molecule_append<'a>(
    left: &TopologyBlock,
    lc: &CoordinateBlock,
    right: &TopologyBlock,
    rc: &'a CoordinateBlock,
) -> Result<PreparedAppend<'a>, IndexedGraphError> {
    // BEGIN RDKIT CPP FUNCTION RWMol::insertMol (Release_2026_09_1)
    // RDKit❗❌: void RWMol::insertMol(const ROMol &other) {
    // RDKit❗❌:   auto origNumAtoms = getNumAtoms();
    // RDKit❗❌:   auto origNumBonds = getNumBonds();
    // RDKit❗❌:   for (const auto oatom : other.atoms()) {
    // RDKit❗❌:     Atom *newAt = oatom->copy();
    // RDKit❗❌:     const bool updateLabel = false;
    // RDKit❗❌:     const bool takeOwnership = true;
    // RDKit❗❌:     addAtom(newAt, updateLabel, takeOwnership);
    // RDKit❗❌:     // take care of atom-numbering-dependent properties:
    // RDKit❗❌:     if (INT_VECT nAtoms;
    // RDKit❗❌:         newAt->getPropIfPresent(common_properties::_ringStereoAtoms, nAtoms)) {
    // RDKit❗❌:       for (auto &val : nAtoms) {
    // RDKit❗❌:         if (val < 0) {
    // RDKit❗❌:           val -= rdcast<int>(origNumAtoms);
    // RDKit❗❌:         } else {
    // RDKit❗❌:           val += origNumAtoms;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       newAt->setProp(common_properties::_ringStereoAtoms, nAtoms, true);
    // RDKit❗❌:     }
    // RDKit❗❌:     if (unsigned int val;
    // RDKit❗❌:         oatom->getPropIfPresent(common_properties::_ringStereoOtherAtom, val)) {
    // RDKit❗❌:       newAt->setProp(common_properties::_ringStereoOtherAtom,
    // RDKit❗❌:                      val + origNumAtoms, true);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // Prepare only incoming rows and fallible remapping. Runtime supplies the
    // owned-or-borrowed left blocks; applying this plan appends in place.
    let (na, nb) = (left.atoms.len(), left.bonds.len());
    let coordinates = if na == 0 && lc.conformers_2d.is_empty() && lc.conformers_3d.is_empty() {
        if rc.conformers_2d.is_empty() && rc.conformers_3d.is_empty() && lc == rc {
            AppendedCoordinates::Unchanged
        } else {
            AppendedCoordinates::Replace(rc)
        }
    } else if right.atoms.is_empty() && rc.conformers_2d.is_empty() && rc.conformers_3d.is_empty() {
        AppendedCoordinates::Unchanged
    } else {
        let two_d = match_coordinates(&lc.conformers_2d, &rc.conformers_2d, 2, Conformer2D::id)?;
        let three_d = match_coordinates(&lc.conformers_3d, &rc.conformers_3d, 3, Conformer3D::id)?;
        if two_d.is_empty() && three_d.is_empty() {
            AppendedCoordinates::Unchanged
        } else {
            AppendedCoordinates::Extend { two_d, three_d }
        }
    };
    let atom_id = |a: AtomId| AtomId::new(na + a.index());
    let bond_id = |b: BondId| BondId::new(nb + b.index());
    let right_map = right
        .atoms
        .iter()
        .map(|a| Some(atom_id(a.id())))
        .collect::<Vec<_>>();
    let mut atoms = Vec::with_capacity(right.atoms.len());
    for old in &right.atoms {
        let mut a = old.clone().with_id(atom_id(old.id()));
        a.remap_template_attachment_order(&right_map)
            .map_err(|source| TopologyEditError::TemplateAttachmentRemap {
                carrier: old.id(),
                source,
            })?;
        remap_ring_stereo(&mut a, &right_map)?;
        atoms.push(a);
    }
    let mut bonds = Vec::with_capacity(right.bonds.len());
    // RDKit❗❌:   for (const auto obond : other.bonds()) {
    // RDKit❗❌:     Bond *bond_p = obond->copy();
    // RDKit❗❌:     unsigned int idx1, idx2;
    // RDKit❗❌:     idx1 = bond_p->getBeginAtomIdx() + origNumAtoms;
    // RDKit❗❌:     idx2 = bond_p->getEndAtomIdx() + origNumAtoms;
    // RDKit❗❌:     bond_p->setOwningMol(this);
    // RDKit❗❌:     bond_p->setBeginAtomIdx(idx1);
    // RDKit❗❌:     bond_p->setEndAtomIdx(idx2);
    // RDKit❗❌:     for (auto &v : bond_p->getStereoAtoms()) {
    // RDKit❗❌:       v += origNumAtoms;
    // RDKit❗❌:     }
    // RDKit❗❌:     const bool takeOwnership = true;
    // RDKit❗❌:     addBond(bond_p, takeOwnership);
    // RDKit❗❌:   }
    for old in &right.bonds {
        let mut b = old.clone().remapped(
            bond_id(old.id()),
            atom_id(old.begin()),
            atom_id(old.end()),
            None,
        );
        b.set_source_stereo_atom_references(
            old.stereo_atom_references()
                .iter()
                .copied()
                .map(atom_id)
                .collect(),
        );
        bonds.push(b);
    }
    let offset = left.substance_groups.len();
    //
    // RDKit❗✔️:   // add stereo groups
    // RDKit❗✔️:   insertStereoGroups(*this, other, origNumAtoms, origNumBonds);
    // RDKit❗✔️:   // add substance groups
    // RDKit❗✔️:   insertSubstanceGroups(*this, other, origNumAtoms, origNumBonds);
    // Substance group crossing bonds and parent IDs are typed indices in this model.
    // END RDKIT CPP FUNCTION RWMol::insertMol
    let mut groups = Vec::with_capacity(right.substance_groups.len());
    for (i, old) in right.substance_groups.iter().enumerate() {
        let mut g = old
            .clone()
            .with_inserted_offsets(SubstanceGroupId::new(i + offset), na, nb)
            .with_head_crossing_bonds(
                old.head_crossing_bonds()
                    .iter()
                    .copied()
                    .map(bond_id)
                    .collect(),
            )
            .with_crossing_bond_correspondence(
                old.crossing_bond_correspondence()
                    .iter()
                    .copied()
                    .map(bond_id)
                    .collect(),
            );
        if let Some(parent) = old.parent() {
            g = g.with_parent(SubstanceGroupId::new(parent.index() + offset));
        }
        groups.push(g);
    }
    let stereo = if right.stereo_groups.is_empty() {
        None
    } else {
        Some(
            cosmolkit_model::insert_stereo_groups(
                &left.stereo_groups,
                &right.stereo_groups,
                na,
                nb,
            )
            .map_err(TopologyEditError::from)?,
        )
    };
    let mapping = TopologyMapping::with_appended(na, nb, right.atoms.len(), right.bonds.len());
    Ok(PreparedAppend {
        atoms,
        bonds,
        groups,
        stereo,
        coordinates,
        mapping,
    })
}

pub fn append_molecule(
    topology: &mut Cow<'_, TopologyBlock>,
    coordinates: &mut Cow<'_, CoordinateBlock>,
    properties: &mut Cow<'_, MoleculeProperties>,
    prepared: PreparedAppend<'_>,
) -> TopologyMapping {
    // RDKit❗❌:     addAtom(newAt, updateLabel, takeOwnership);
    // RDKit❗❌:     addBond(bond_p, takeOwnership);
    // Source ownership is respected: append new rows to unique left storage;
    // shared left storage detaches once. CSR adjacency and explicit mapping
    // still require O(V+E) rebuilding, unlike the upstream mutable graph.
    if !prepared.atoms.is_empty()
        || !prepared.bonds.is_empty()
        || !prepared.groups.is_empty()
        || prepared.stereo.is_some()
    {
        let topology = topology.to_mut();
        topology.atoms.extend(prepared.atoms);
        topology.bonds.extend(prepared.bonds);
        topology.substance_groups.extend(prepared.groups);
        if let Some(stereo) = prepared.stereo {
            topology.stereo_groups = stereo;
        }
        topology.adjacency =
            cosmolkit_model::AdjacencyList::from_topology(topology.atoms.len(), &topology.bonds);
    }
    match prepared.coordinates {
        AppendedCoordinates::Unchanged => {}
        AppendedCoordinates::Replace(right) => *coordinates = Cow::Owned(right.clone()),
        AppendedCoordinates::Extend { two_d, three_d } => {
            let left = coordinates.to_mut();
            for (a, b) in left.conformers_2d.iter_mut().zip(two_d) {
                a.append_atom_positions(b.coordinates());
            }
            for (a, b) in left.conformers_3d.iter_mut().zip(three_d) {
                a.append_atom_positions(b.coordinates());
            }
        }
    }
    if !properties.sdf_property_lists().is_empty() {
        properties.to_mut().remap_topology(
            &prepared.mapping.atoms.new_to_old,
            &prepared.mapping.bonds.new_to_old,
        );
    }
    prepared.mapping
}
