//! Lipinski descriptor implementation boundary.
//!
//! This module exclusively owns direct Lipinski counts, descriptor SMARTS,
//! ring-family classification, spiro and bridgehead counts, and stereo-center
//! projections. Ring perception, SMARTS matching, and stereo assignment remain
//! shared chemistry infrastructure.

use super::{DescriptorError, DescriptorResult, rdkit_count_query_matches, rdkit_count_smarts_matches};
use crate::{AtomId, BondOrder, ChiralTag, Molecule};
use std::borrow::Cow;

const RDKIT_NUM_HETEROATOMS_SMARTS: &str = "[!#6;!#1]";
const RDKIT_NUM_AMIDE_BONDS_SMARTS: &str = "C(=[O;!R])N";

fn has_stereo_assigned(molecule: &Molecule) -> bool {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::hasStereoAssigned
    // RDKit✔️✔️: bool hasStereoAssigned(const ROMol &mol) {
    // RDKit✔️✔️:   return mol.hasProp(common_properties::_StereochemDone);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::hasStereoAssigned
    molecule.prop("_StereochemDone").is_some()
}

fn stereo_assigned_for_descriptor(molecule: &Molecule) -> DescriptorResult<Cow<'_, Molecule>> {
    if has_stereo_assigned(molecule) {
        return Ok(Cow::Borrowed(molecule));
    }

    let mut working = molecule.clone();
    crate::smiles_write::assign_stereochemistry_on_working_copy(&mut working, true).map_err(|_| {
        DescriptorError::Unsupported {
            function: "stereo_center_descriptors",
            rdkit_function: "MolOps::assignStereochemistry",
        }
    })?;
    Ok(Cow::Owned(working))
}

pub(super) fn num_atom_stereo_centers(molecule: &Molecule) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::numAtomStereoCenters
    // RDKit✔️✔️: unsigned int numAtomStereoCenters(const ROMol &mol) {
    // RDKit✔️✔️:   std::unique_ptr<ROMol> tmol;
    // RDKit✔️✔️:   const ROMol *mptr = &mol;
    // RDKit✔️✔️:   if (!hasStereoAssigned(mol)) {
    // RDKit✔️✔️:     tmol.reset(new ROMol(mol));
    // RDKit✔️✔️:     constexpr bool cleanIt = true;
    // RDKit✔️✔️:     constexpr bool force = true;
    // RDKit✔️✔️:     constexpr bool flagPossible = true;
    // RDKit✔️✔️:     MolOps::assignStereochemistry(*tmol, cleanIt, force, flagPossible);
    // RDKit✔️✔️:     mptr = tmol.get();
    // RDKit✔️✔️:   }
    let counted = stereo_assigned_for_descriptor(molecule)?;

    // RDKit✔️✔️:   unsigned int res = 0;
    let mut result = 0_u32;
    // RDKit✔️✔️:   for (const auto &atom : mptr->atoms()) {
    for atom in counted.atoms() {
        // RDKit✔️✔️:     if (atom->hasProp(common_properties::_ChiralityPossible)) {
        // RDKit✔️✔️:       ++res;
        // RDKit✔️✔️:     }
        if atom.prop("_ChiralityPossible").is_some() {
            result += 1;
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::numAtomStereoCenters
    Ok(result)
}

pub(super) fn num_unspecified_atom_stereo_centers(molecule: &Molecule) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::numUnspecifiedAtomStereoCenters
    // RDKit✔️✔️: unsigned int numUnspecifiedAtomStereoCenters(const ROMol &mol) {
    // RDKit✔️✔️:   std::unique_ptr<ROMol> tmol;
    // RDKit✔️✔️:   const ROMol *mptr = &mol;
    // RDKit✔️✔️:   if (!hasStereoAssigned(mol)) {
    // RDKit✔️✔️:     tmol.reset(new ROMol(mol));
    // RDKit✔️✔️:     constexpr bool cleanIt = true;
    // RDKit✔️✔️:     constexpr bool force = true;
    // RDKit✔️✔️:     constexpr bool flagPossible = true;
    // RDKit✔️✔️:     MolOps::assignStereochemistry(*tmol, cleanIt, force, flagPossible);
    // RDKit✔️✔️:     mptr = tmol.get();
    // RDKit✔️✔️:   }
    let counted = stereo_assigned_for_descriptor(molecule)?;

    // RDKit✔️✔️:   unsigned int res = 0;
    let mut result = 0_u32;
    // RDKit✔️✔️:   for (const auto &atom : mptr->atoms()) {
    for atom in counted.atoms() {
        // RDKit✔️✔️:     if (atom->hasProp(common_properties::_ChiralityPossible) &&
        // RDKit✔️✔️:         atom->getChiralTag() == Atom::CHI_UNSPECIFIED) {
        // RDKit✔️✔️:       ++res;
        // RDKit✔️✔️:     }
        if atom.prop("_ChiralityPossible").is_some() && atom.chiral_tag() == ChiralTag::Unspecified {
            result += 1;
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::numUnspecifiedAtomStereoCenters
    Ok(result)
}

pub(super) fn calc_lipinski_hba(molecule: &Molecule) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcLipinskiHBA
    // RDKit✔️✔️: unsigned int calcLipinskiHBA(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    let mut result = 0_u32;
    // RDKit✔️✔️:   for (ROMol::ConstAtomIterator iter = mol.beginAtoms(); iter != mol.endAtoms();
    // RDKit✔️✔️:        ++iter) {
    for atom in molecule.atoms() {
        // RDKit✔️✔️:     if ((*iter)->getAtomicNum() == 7 || (*iter)->getAtomicNum() == 8) {
        // RDKit✔️✔️:       ++res;
        // RDKit✔️✔️:     }
        if atom.atomic_number() == 7 || atom.atomic_number() == 8 {
            result += 1;
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcLipinskiHBA
    Ok(result)
}

pub(super) fn calc_lipinski_hbd(molecule: &Molecule) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcLipinskiHBD
    // RDKit✔️✔️: unsigned int calcLipinskiHBD(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    let mut result = 0_u32;
    // RDKit✔️✔️:   for (ROMol::ConstAtomIterator iter = mol.beginAtoms(); iter != mol.endAtoms();
    // RDKit✔️✔️:        ++iter) {
    for atom in molecule.atoms() {
        // RDKit✔️✔️:     if (((*iter)->getAtomicNum() == 7 || (*iter)->getAtomicNum() == 8)) {
        if atom.atomic_number() == 7 || atom.atomic_number() == 8 {
            // RDKit✔️✔️:       res += (*iter)->getTotalNumHs(true);
            result += crate::chemistry::valence::total_num_hydrogens(molecule, atom.id(), true).map_err(|_| {
                DescriptorError::Unsupported {
                    function: "calc_lipinski_hbd",
                    rdkit_function: "Atom::getTotalNumHs(true)",
                }
            })?;
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcLipinskiHBD
    Ok(result)
}

pub(super) fn calc_num_heteroatoms(molecule: &Molecule) -> DescriptorResult<u32> {
    // RDKit✔️✔️: #define SMARTSCOUNTFUNC(nm, pattern, vers)         \
    // RDKit✔️✔️:   const std::string nm##Version = vers;            \
    // RDKit✔️✔️:   unsigned int calc##nm(const RDKit::ROMol &mol) { \
    // RDKit✔️✔️:     pattern_flyweight m(pattern);                  \
    // RDKit✔️✔️:     return m.get().countMatches(mol);              \
    // RDKit✔️✔️:   }                                                \
    // RDKit✔️✔️:   extern int no_such_variable
    // RDKit✔️✔️: SMARTSCOUNTFUNC(NumHeteroatoms, "[!#6;!#1]", "1.0.1");
    rdkit_count_smarts_matches("calc_num_heteroatoms", RDKIT_NUM_HETEROATOMS_SMARTS)
        .and_then(|query| rdkit_count_query_matches(molecule, &query, "calc_num_heteroatoms"))
}

pub(super) fn calc_num_amide_bonds(molecule: &Molecule) -> DescriptorResult<u32> {
    // RDKit✔️✔️: #define SMARTSCOUNTFUNC(nm, pattern, vers)         \
    // RDKit✔️✔️:   const std::string nm##Version = vers;            \
    // RDKit✔️✔️:   unsigned int calc##nm(const RDKit::ROMol &mol) { \
    // RDKit✔️✔️:     pattern_flyweight m(pattern);                  \
    // RDKit✔️✔️:     return m.get().countMatches(mol);              \
    // RDKit✔️✔️:   }                                                \
    // RDKit✔️✔️:   extern int no_such_variable
    // RDKit✔️✔️: SMARTSCOUNTFUNC(NumAmideBonds, "C(=[O;!R])N", "1.0.0");
    rdkit_count_smarts_matches("calc_num_amide_bonds", RDKIT_NUM_AMIDE_BONDS_SMARTS)
        .and_then(|query| rdkit_count_query_matches(molecule, &query, "calc_num_amide_bonds"))
}

pub(super) fn calc_num_heavy_atoms(molecule: &Molecule) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumHeavyAtoms
    // RDKit✔️✔️: unsigned int calcNumHeavyAtoms(const ROMol &mol) {
    // RDKit✔️✔️:   return mol.getNumHeavyAtoms();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumHeavyAtoms
    // BEGIN RDKIT CPP FUNCTION: RDKit::ROMol::getNumHeavyAtoms
    // RDKit✔️✔️: unsigned int ROMol::getNumHeavyAtoms() const {
    // RDKit✔️✔️:   unsigned int res = 0;
    let mut result = 0_u32;
    // RDKit✔️✔️:   for (const auto atom : atoms()) {
    for atom in molecule.atoms() {
        // RDKit✔️✔️:     if (atom->getAtomicNum() > 1) {
        // RDKit✔️✔️:       ++res;
        // RDKit✔️✔️:     }
        if atom.atomic_number() > 1 {
            result += 1;
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION: RDKit::ROMol::getNumHeavyAtoms
    Ok(result)
}

pub(super) fn calc_num_atoms(molecule: &Molecule) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumAtoms
    // RDKit✔️✔️: unsigned int calcNumAtoms(const ROMol &mol) {
    // RDKit✔️✔️:   bool onlyExplicit = false;
    // RDKit✔️✔️:   return mol.getNumAtoms(onlyExplicit);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumAtoms
    // BEGIN RDKIT CPP FUNCTION: RDKit::ROMol::getNumAtoms
    // RDKit✔️✔️: unsigned int ROMol::getNumAtoms(bool onlyExplicit) const {
    // RDKit✔️✔️:   int res = rdcast<int>(boost::num_vertices(d_graph));
    let mut result = u32::try_from(molecule.num_atoms()).map_err(|_| DescriptorError::Unsupported {
        function: "calc_num_atoms",
        rdkit_function: "ROMol::getNumAtoms explicit vertex count",
    })?;
    // RDKit✔️✔️:   if (!onlyExplicit) {
    // RDKit✔️✔️:     // if we are interested in hydrogens as well add them up from
    // RDKit✔️✔️:     // each
    // RDKit✔️✔️:     for (const auto atom : atoms()) {
    for atom in molecule.atoms() {
        // RDKit✔️✔️:       res += atom->getTotalNumHs();
        result = result
            .checked_add(
                crate::chemistry::valence::total_num_hydrogens(molecule, atom.id(), false).map_err(|_| {
                    DescriptorError::Unsupported {
                        function: "calc_num_atoms",
                        rdkit_function: "Atom::getTotalNumHs(false)",
                    }
                })?,
            )
            .ok_or(DescriptorError::Unsupported {
                function: "calc_num_atoms",
                rdkit_function: "ROMol::getNumAtoms unsigned result",
            })?;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION: RDKit::ROMol::getNumAtoms
    Ok(result)
}

pub(super) fn calc_num_rings(molecule: &Molecule) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumRings
    // RDKit✔️✔️: unsigned int calcNumRings(const ROMol &mol) {
    // RDKit✔️✔️:   return mol.getRingInfo()->numRings();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumRings
    let rings = super::descriptor_ring_info(molecule, false).map_err(|_| DescriptorError::Unsupported {
        function: "calc_num_rings",
        rdkit_function: "RingInfo::numRings",
    })?;
    u32::try_from(rings.num_rings()).map_err(|_| DescriptorError::Unsupported {
        function: "calc_num_rings",
        rdkit_function: "RingInfo::numRings unsigned result",
    })
}

pub(super) fn calc_num_heterocycles(molecule: &Molecule) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumHeterocycles
    // RDKit✔️✔️: unsigned int calcNumHeterocycles(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    let mut result = 0_u32;
    let rings = super::descriptor_ring_info(molecule, false).map_err(|_| DescriptorError::Unsupported {
        function: "calc_num_heterocycles",
        rdkit_function: "RingInfo::atomRings",
    })?;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->atomRings()) {
    for ring in rings.atom_rings() {
        // RDKit✔️✔️:     for (auto i : iv) {
        for atom_id in ring {
            // RDKit✔️✔️:       if (mol.getAtomWithIdx(i)->getAtomicNum() != 6) {
            // RDKit✔️✔️:         ++res;
            // RDKit✔️✔️:         break;
            // RDKit✔️✔️:       }
            if molecule.atoms()[atom_id.index()].atomic_number() != 6 {
                result += 1;
                break;
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumHeterocycles
    Ok(result)
}

pub(super) fn calc_num_saturated_rings(molecule: &Molecule) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumSaturatedRings
    // RDKit✔️✔️: unsigned int calcNumSaturatedRings(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    let mut result = 0_u32;
    let rings = super::descriptor_ring_info(molecule, false).map_err(|_| DescriptorError::Unsupported {
        function: "calc_num_saturated_rings",
        rdkit_function: "RingInfo::bondRings",
    })?;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->bondRings()) {
    for ring in rings.bond_rings() {
        // RDKit✔️✔️:     ++res;
        result += 1;
        // RDKit✔️✔️:     for (int i : iv) {
        for bond_id in ring {
            let bond = &molecule.bonds()[bond_id.index()];
            // RDKit✔️✔️:       if (mol.getBondWithIdx(i)->getBondType() != Bond::SINGLE ||
            // RDKit✔️✔️:           mol.getBondWithIdx(i)->getIsAromatic()) {
            // RDKit✔️✔️:         --res;
            // RDKit✔️✔️:         break;
            // RDKit✔️✔️:       }
            if bond.order() != BondOrder::Single || bond.is_aromatic() {
                result -= 1;
                break;
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumSaturatedRings
    Ok(result)
}

pub(super) fn calc_num_aliphatic_rings(molecule: &Molecule) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumAliphaticRings
    // RDKit✔️✔️: unsigned int calcNumAliphaticRings(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    let mut result = 0_u32;
    let rings = super::descriptor_ring_info(molecule, false).map_err(|_| DescriptorError::Unsupported {
        function: "calc_num_aliphatic_rings",
        rdkit_function: "RingInfo::bondRings",
    })?;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->bondRings()) {
    for ring in rings.bond_rings() {
        // RDKit✔️✔️:     for (auto i : iv) {
        for bond_id in ring {
            // RDKit✔️✔️:       if (!mol.getBondWithIdx(i)->getIsAromatic()) {
            // RDKit✔️✔️:         ++res;
            // RDKit✔️✔️:         break;
            // RDKit✔️✔️:       }
            if !molecule.bonds()[bond_id.index()].is_aromatic() {
                result += 1;
                break;
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumAliphaticRings
    Ok(result)
}

pub(super) fn calc_num_aromatic_heterocycles(molecule: &Molecule) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumAromaticHeterocycles
    // RDKit✔️✔️: unsigned int calcNumAromaticHeterocycles(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    let mut result = 0_u32;
    let rings = super::descriptor_ring_info(molecule, false).map_err(|_| DescriptorError::Unsupported {
        function: "calc_num_aromatic_heterocycles",
        rdkit_function: "RingInfo::bondRings",
    })?;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->bondRings()) {
    for ring in rings.bond_rings() {
        // RDKit✔️✔️:     bool countIt = false;
        let mut count_it = false;
        // RDKit✔️✔️:     for (auto i : iv) {
        for bond_id in ring {
            let bond = &molecule.bonds()[bond_id.index()];
            // RDKit✔️✔️:       if (!mol.getBondWithIdx(i)->getIsAromatic()) {
            // RDKit✔️✔️:         countIt = false;
            // RDKit✔️✔️:         break;
            // RDKit✔️✔️:       }
            if !bond.is_aromatic() {
                count_it = false;
                break;
            }
            // RDKit✔️✔️:       // we're checking each atom twice, which is kind of doofy, but this
            // RDKit✔️✔️:       // function is hopefully not going to be a big time sink.
            // RDKit✔️✔️:       if (!countIt &&
            // RDKit✔️✔️:           (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
            // RDKit✔️✔️:            mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6)) {
            // RDKit✔️✔️:         countIt = true;
            // RDKit✔️✔️:       }
            if !count_it
                && (molecule.atoms()[bond.begin().index()].atomic_number() != 6
                    || molecule.atoms()[bond.end().index()].atomic_number() != 6)
            {
                count_it = true;
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:     if (countIt) {
        // RDKit✔️✔️:       ++res;
        // RDKit✔️✔️:     }
        if count_it {
            result += 1;
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumAromaticHeterocycles
    Ok(result)
}

pub(super) fn calc_num_aromatic_carbocycles(molecule: &Molecule) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumAromaticCarbocycles
    // RDKit✔️✔️: unsigned int calcNumAromaticCarbocycles(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    let mut result = 0_u32;
    let rings = super::descriptor_ring_info(molecule, false).map_err(|_| DescriptorError::Unsupported {
        function: "calc_num_aromatic_carbocycles",
        rdkit_function: "RingInfo::bondRings",
    })?;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->bondRings()) {
    for ring in rings.bond_rings() {
        // RDKit✔️✔️:     bool countIt = true;
        let mut count_it = true;
        // RDKit✔️✔️:     for (auto i : iv) {
        for bond_id in ring {
            let bond = &molecule.bonds()[bond_id.index()];
            // RDKit✔️✔️:       if (!mol.getBondWithIdx(i)->getIsAromatic()) {
            // RDKit✔️✔️:         countIt = false;
            // RDKit✔️✔️:         break;
            // RDKit✔️✔️:       }
            if !bond.is_aromatic() {
                count_it = false;
                break;
            }
            // RDKit✔️✔️:       // we're checking each atom twice, which is kind of doofy, but this
            // RDKit✔️✔️:       // function is hopefully not going to be a big time sync.
            // RDKit✔️✔️:       if (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
            // RDKit✔️✔️:           mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6) {
            // RDKit✔️✔️:         countIt = false;
            // RDKit✔️✔️:         break;
            // RDKit✔️✔️:       }
            if molecule.atoms()[bond.begin().index()].atomic_number() != 6
                || molecule.atoms()[bond.end().index()].atomic_number() != 6
            {
                count_it = false;
                break;
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:     if (countIt) {
        // RDKit✔️✔️:       ++res;
        // RDKit✔️✔️:     }
        if count_it {
            result += 1;
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumAromaticCarbocycles
    Ok(result)
}

pub(super) fn calc_num_aliphatic_heterocycles(molecule: &Molecule) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumAliphaticHeterocycles
    // RDKit✔️✔️: unsigned int calcNumAliphaticHeterocycles(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    let mut result = 0_u32;
    let rings = super::descriptor_ring_info(molecule, false).map_err(|_| DescriptorError::Unsupported {
        function: "calc_num_aliphatic_heterocycles",
        rdkit_function: "RingInfo::bondRings",
    })?;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->bondRings()) {
    for ring in rings.bond_rings() {
        // RDKit✔️✔️:     bool hasAliph = false;
        // RDKit✔️✔️:     bool hasHetero = false;
        let mut has_aliphatic = false;
        let mut has_heteroatom = false;
        // RDKit✔️✔️:     for (auto i : iv) {
        for bond_id in ring {
            let bond = &molecule.bonds()[bond_id.index()];
            // RDKit✔️✔️:       if (!mol.getBondWithIdx(i)->getIsAromatic()) {
            // RDKit✔️✔️:         hasAliph = true;
            // RDKit✔️✔️:       }
            if !bond.is_aromatic() {
                has_aliphatic = true;
            }
            // RDKit✔️✔️:       // we're checking each atom twice, which is kind of doofy, but this
            // RDKit✔️✔️:       // function is hopefully not going to be a big time sink.
            // RDKit✔️✔️:       if (!hasHetero &&
            // RDKit✔️✔️:           (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
            // RDKit✔️✔️:            mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6)) {
            // RDKit✔️✔️:         hasHetero = true;
            // RDKit✔️✔️:       }
            if !has_heteroatom
                && (molecule.atoms()[bond.begin().index()].atomic_number() != 6
                    || molecule.atoms()[bond.end().index()].atomic_number() != 6)
            {
                has_heteroatom = true;
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:     if (hasHetero && hasAliph) {
        // RDKit✔️✔️:       ++res;
        // RDKit✔️✔️:     }
        if has_heteroatom && has_aliphatic {
            result += 1;
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumAliphaticHeterocycles
    Ok(result)
}

pub(super) fn calc_num_aliphatic_carbocycles(molecule: &Molecule) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumAliphaticCarbocycles
    // RDKit✔️✔️: unsigned int calcNumAliphaticCarbocycles(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    let mut result = 0_u32;
    let rings = super::descriptor_ring_info(molecule, false).map_err(|_| DescriptorError::Unsupported {
        function: "calc_num_aliphatic_carbocycles",
        rdkit_function: "RingInfo::bondRings",
    })?;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->bondRings()) {
    for ring in rings.bond_rings() {
        // RDKit✔️✔️:     bool hasAliph = false;
        // RDKit✔️✔️:     bool hasHetero = false;
        let mut has_aliphatic = false;
        let mut has_heteroatom = false;
        // RDKit✔️✔️:     for (auto i : iv) {
        for bond_id in ring {
            let bond = &molecule.bonds()[bond_id.index()];
            // RDKit✔️✔️:       if (!mol.getBondWithIdx(i)->getIsAromatic()) {
            // RDKit✔️✔️:         hasAliph = true;
            // RDKit✔️✔️:       }
            if !bond.is_aromatic() {
                has_aliphatic = true;
            }
            // RDKit✔️✔️:       // we're checking each atom twice, which is kind of doofy, but this
            // RDKit✔️✔️:       // function is hopefully not going to be a big time sync.
            // RDKit✔️✔️:       if (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
            // RDKit✔️✔️:           mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6) {
            // RDKit✔️✔️:         hasHetero = true;
            // RDKit✔️✔️:         break;
            // RDKit✔️✔️:       }
            if molecule.atoms()[bond.begin().index()].atomic_number() != 6
                || molecule.atoms()[bond.end().index()].atomic_number() != 6
            {
                has_heteroatom = true;
                break;
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:     if (hasAliph && !hasHetero) {
        // RDKit✔️✔️:       ++res;
        // RDKit✔️✔️:     }
        if has_aliphatic && !has_heteroatom {
            result += 1;
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumAliphaticCarbocycles
    Ok(result)
}

pub(super) fn calc_num_saturated_heterocycles(molecule: &Molecule) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumSaturatedHeterocycles
    // RDKit✔️✔️: unsigned int calcNumSaturatedHeterocycles(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    let mut result = 0_u32;
    let rings = super::descriptor_ring_info(molecule, false).map_err(|_| DescriptorError::Unsupported {
        function: "calc_num_saturated_heterocycles",
        rdkit_function: "RingInfo::bondRings",
    })?;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->bondRings()) {
    for ring in rings.bond_rings() {
        // RDKit✔️✔️:     bool countIt = false;
        let mut count_it = false;
        // RDKit✔️✔️:     for (auto i : iv) {
        for bond_id in ring {
            let bond = &molecule.bonds()[bond_id.index()];
            // RDKit✔️✔️:       if (mol.getBondWithIdx(i)->getBondType() != Bond::SINGLE ||
            // RDKit✔️✔️:           mol.getBondWithIdx(i)->getIsAromatic()) {
            // RDKit✔️✔️:         countIt = false;
            // RDKit✔️✔️:         break;
            // RDKit✔️✔️:       }
            if bond.order() != BondOrder::Single || bond.is_aromatic() {
                count_it = false;
                break;
            }
            // RDKit✔️✔️:       // we're checking each atom twice, which is kind of doofy, but this
            // RDKit✔️✔️:       // function is hopefully not going to be a big time sync.
            // RDKit✔️✔️:       if (!countIt &&
            // RDKit✔️✔️:           (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
            // RDKit✔️✔️:            mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6)) {
            // RDKit✔️✔️:         countIt = true;
            // RDKit✔️✔️:       }
            if !count_it
                && (molecule.atoms()[bond.begin().index()].atomic_number() != 6
                    || molecule.atoms()[bond.end().index()].atomic_number() != 6)
            {
                count_it = true;
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:     if (countIt) {
        // RDKit✔️✔️:       ++res;
        // RDKit✔️✔️:     }
        if count_it {
            result += 1;
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumSaturatedHeterocycles
    Ok(result)
}

pub(super) fn calc_num_saturated_carbocycles(molecule: &Molecule) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumSaturatedCarbocycles
    // RDKit✔️✔️: unsigned int calcNumSaturatedCarbocycles(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    let mut result = 0_u32;
    let rings = super::descriptor_ring_info(molecule, false).map_err(|_| DescriptorError::Unsupported {
        function: "calc_num_saturated_carbocycles",
        rdkit_function: "RingInfo::bondRings",
    })?;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->bondRings()) {
    for ring in rings.bond_rings() {
        // RDKit✔️✔️:     bool countIt = true;
        let mut count_it = true;
        // RDKit✔️✔️:     for (auto i : iv) {
        for bond_id in ring {
            let bond = &molecule.bonds()[bond_id.index()];
            // RDKit✔️✔️:       if (mol.getBondWithIdx(i)->getBondType() != Bond::SINGLE ||
            // RDKit✔️✔️:           mol.getBondWithIdx(i)->getIsAromatic()) {
            // RDKit✔️✔️:         countIt = false;
            // RDKit✔️✔️:         break;
            // RDKit✔️✔️:       }
            if bond.order() != BondOrder::Single || bond.is_aromatic() {
                count_it = false;
                break;
            }
            // RDKit✔️✔️:       // we're checking each atom twice, which is kind of doofy, but this
            // RDKit✔️✔️:       // function is hopefully not going to be a big time sync.
            // RDKit✔️✔️:       if (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
            // RDKit✔️✔️:           mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6) {
            // RDKit✔️✔️:         countIt = false;
            // RDKit✔️✔️:         break;
            // RDKit✔️✔️:       }
            if molecule.atoms()[bond.begin().index()].atomic_number() != 6
                || molecule.atoms()[bond.end().index()].atomic_number() != 6
            {
                count_it = false;
                break;
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:     if (countIt) {
        // RDKit✔️✔️:       ++res;
        // RDKit✔️✔️:     }
        if count_it {
            result += 1;
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumSaturatedCarbocycles
    Ok(result)
}

pub(super) fn calc_num_spiro_atoms(molecule: &Molecule, atoms: Option<&mut Vec<AtomId>>) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumSpiroAtoms
    // RDKit✔️✔️: unsigned int calcNumSpiroAtoms(const ROMol &mol,
    // RDKit✔️✔️:                                std::vector<unsigned int> *atoms) {
    // RDKit✔️✔️:   if (!mol.getRingInfo() || !mol.getRingInfo()->isSssrOrBetter()) {
    // RDKit✔️✔️:     MolOps::findSSSR(mol);
    // RDKit✔️✔️:   }
    let rings = super::descriptor_ring_info(molecule, true).map_err(|_| DescriptorError::Unsupported {
        function: "calc_num_spiro_atoms",
        rdkit_function: "MolOps::findSSSR/RingInfo::atomRings",
    })?;
    // RDKit✔️✔️:   const RingInfo *rInfo = mol.getRingInfo();
    // RDKit✔️✔️:   std::vector<unsigned int> lAtoms;
    let mut local_atoms = Vec::new();
    // RDKit✔️✔️:   if (!atoms) {
    // RDKit✔️✔️:     atoms = &lAtoms;
    // RDKit✔️✔️:   }
    let atoms = atoms.unwrap_or(&mut local_atoms);

    // RDKit✔️✔️:   for (unsigned int i = 0; i < rInfo->atomRings().size(); ++i) {
    for left_index in 0..rings.atom_rings().len() {
        // RDKit✔️✔️:     const INT_VECT &ri = rInfo->atomRings()[i];
        let left_ring = &rings.atom_rings()[left_index];
        // RDKit✔️✔️:     for (unsigned int j = i + 1; j < rInfo->atomRings().size(); ++j) {
        for right_index in left_index + 1..rings.atom_rings().len() {
            // RDKit✔️✔️:       const INT_VECT &rj = rInfo->atomRings()[j];
            let right_ring = &rings.atom_rings()[right_index];
            // RDKit✔️✔️:       // EFF: using intersect here does more work and memory allocation than is
            // RDKit✔️✔️:       // required
            // RDKit✔️✔️:       INT_VECT inter;
            // RDKit✔️✔️:       Intersect(ri, rj, inter);
            let intersection = left_ring
                .iter()
                .copied()
                .filter(|atom| right_ring.contains(atom))
                .collect::<Vec<_>>();
            // RDKit✔️✔️:       if (inter.size() == 1) {
            if intersection.len() == 1 {
                // RDKit✔️✔️:         if (std::find(atoms->begin(), atoms->end(), inter[0]) == atoms->end()) {
                // RDKit✔️✔️:           atoms->push_back(inter[0]);
                // RDKit✔️✔️:         }
                if !atoms.contains(&intersection[0]) {
                    atoms.push(intersection[0]);
                }
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return atoms->size();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumSpiroAtoms
    u32::try_from(atoms.len()).map_err(|_| DescriptorError::Unsupported {
        function: "calc_num_spiro_atoms",
        rdkit_function: "calcNumSpiroAtoms unsigned result",
    })
}

pub(super) fn calc_num_bridgehead_atoms(molecule: &Molecule, atoms: Option<&mut Vec<AtomId>>) -> DescriptorResult<u32> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumBridgeheadAtoms
    // RDKit✔️✔️: unsigned int calcNumBridgeheadAtoms(const ROMol &mol,
    // RDKit✔️✔️:                                     std::vector<unsigned int> *atoms) {
    // RDKit✔️✔️:   if (!mol.getRingInfo() || !mol.getRingInfo()->isSssrOrBetter()) {
    // RDKit✔️✔️:     MolOps::findSSSR(mol);
    // RDKit✔️✔️:   }
    let rings = super::descriptor_ring_info(molecule, true).map_err(|_| DescriptorError::Unsupported {
        function: "calc_num_bridgehead_atoms",
        rdkit_function: "MolOps::findSSSR/RingInfo::bondRings",
    })?;
    // RDKit✔️✔️:   const RingInfo *rInfo = mol.getRingInfo();
    // RDKit✔️✔️:   std::vector<unsigned int> lAtoms;
    let mut local_atoms = Vec::new();
    // RDKit✔️✔️:   if (!atoms) {
    // RDKit✔️✔️:     atoms = &lAtoms;
    // RDKit✔️✔️:   }
    let atoms = atoms.unwrap_or(&mut local_atoms);

    // RDKit✔️✔️:   for (unsigned int i = 0; i < rInfo->bondRings().size(); ++i) {
    for left_index in 0..rings.bond_rings().len() {
        // RDKit✔️✔️:     const INT_VECT &ri = rInfo->bondRings()[i];
        let left_ring = &rings.bond_rings()[left_index];
        // RDKit✔️✔️:     for (unsigned int j = i + 1; j < rInfo->bondRings().size(); ++j) {
        for right_index in left_index + 1..rings.bond_rings().len() {
            // RDKit✔️✔️:       const INT_VECT &rj = rInfo->bondRings()[j];
            let right_ring = &rings.bond_rings()[right_index];
            // RDKit✔️✔️:       // EFF: using intersect here does more work and memory allocation than is
            // RDKit✔️✔️:       // required
            // RDKit✔️✔️:       INT_VECT inter;
            // RDKit✔️✔️:       Intersect(ri, rj, inter);
            let intersection = left_ring
                .iter()
                .copied()
                .filter(|bond| right_ring.contains(bond))
                .collect::<Vec<_>>();
            // RDKit✔️✔️:       if (inter.size() > 1) {
            if intersection.len() > 1 {
                // RDKit✔️✔️:         INT_VECT atomCounts(mol.getNumAtoms(), 0);
                let mut atom_counts = vec![0_u32; molecule.atoms().len()];
                // RDKit✔️✔️:         for (auto ii : inter) {
                for bond_id in intersection {
                    let bond = &molecule.bonds()[bond_id.index()];
                    // RDKit✔️✔️:           atomCounts[mol.getBondWithIdx(ii)->getBeginAtomIdx()] += 1;
                    // RDKit✔️✔️:           atomCounts[mol.getBondWithIdx(ii)->getEndAtomIdx()] += 1;
                    atom_counts[bond.begin().index()] += 1;
                    atom_counts[bond.end().index()] += 1;
                    // RDKit✔️✔️:         }
                }
                // RDKit✔️✔️:         for (unsigned int ti = 0; ti < atomCounts.size(); ++ti) {
                for (atom_index, count) in atom_counts.into_iter().enumerate() {
                    // RDKit✔️✔️:           if (atomCounts[ti] == 1) {
                    if count == 1 {
                        let atom_id = AtomId::new(atom_index);
                        // RDKit✔️✔️:             if (std::find(atoms->begin(), atoms->end(), ti) == atoms->end()) {
                        // RDKit✔️✔️:               atoms->push_back(ti);
                        // RDKit✔️✔️:             }
                        if !atoms.contains(&atom_id) {
                            atoms.push(atom_id);
                        }
                        // RDKit✔️✔️:           }
                    }
                    // RDKit✔️✔️:         }
                }
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return atoms->size();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcNumBridgeheadAtoms
    u32::try_from(atoms.len()).map_err(|_| DescriptorError::Unsupported {
        function: "calc_num_bridgehead_atoms",
        rdkit_function: "calcNumBridgeheadAtoms unsigned result",
    })
}

#[cfg(test)]
mod tests {
    use super::{
        calc_lipinski_hba, calc_lipinski_hbd, calc_num_aliphatic_carbocycles, calc_num_aliphatic_heterocycles,
        calc_num_aliphatic_rings, calc_num_amide_bonds, calc_num_aromatic_carbocycles, calc_num_aromatic_heterocycles,
        calc_num_atoms, calc_num_bridgehead_atoms, calc_num_heavy_atoms, calc_num_heteroatoms, calc_num_heterocycles,
        calc_num_rings, calc_num_saturated_carbocycles, calc_num_saturated_heterocycles, calc_num_saturated_rings,
        calc_num_spiro_atoms, has_stereo_assigned, num_atom_stereo_centers, num_unspecified_atom_stereo_centers,
    };
    use crate::{AtomId, AtomSpec, BondOrder, BondSpec, Element, Molecule, MoleculeBuilder};

    fn explicit_methylamine() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let nitrogen = builder.add_atom(AtomSpec::new(Element::N));
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(nitrogen, carbon, BondOrder::Single))
            .expect("N-C bond");
        for _ in 0..2 {
            let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
            builder
                .add_bond(BondSpec::new(nitrogen, hydrogen, BondOrder::Single))
                .expect("N-H bond");
        }
        builder
            .build()
            .expect("methylamine graph")
            .sanitize()
            .expect("sanitized methylamine")
    }

    #[test]
    fn direct_and_standard_lipinski_counts_match_pinned_rdkit_distinct_definitions() {
        const CASES: [(&str, [u32; 8]); 7] = [
            ("", [0, 0, 0, 0, 0, 0, 0, 0]),
            ("CCO", [1, 1, 1, 1, 1, 0, 3, 9]),
            ("NC(=O)C", [2, 2, 1, 1, 2, 1, 4, 9]),
            ("NC(=O)N", [3, 4, 1, 2, 3, 2, 4, 8]),
            ("NCC(=O)O", [3, 3, 2, 2, 3, 0, 5, 10]),
            ("c1ncc[nH]1", [2, 1, 1, 1, 2, 0, 5, 9]),
            ("[Na+].[Cl-]", [0, 0, 0, 0, 2, 0, 2, 2]),
        ];

        for (smiles, expected) in CASES {
            let molecule = Molecule::from_smiles(smiles).expect("Lipinski fixture");
            let actual = [
                calc_lipinski_hba(&molecule).unwrap(),
                calc_lipinski_hbd(&molecule).unwrap(),
                crate::properties::descriptors::calc_num_hba(&molecule).unwrap(),
                crate::properties::descriptors::calc_num_hbd(&molecule).unwrap(),
                calc_num_heteroatoms(&molecule).unwrap(),
                calc_num_amide_bonds(&molecule).unwrap(),
                calc_num_heavy_atoms(&molecule).unwrap(),
                calc_num_atoms(&molecule).unwrap(),
            ];
            assert_eq!(actual, expected, "{smiles:?} Lipinski count vector");
        }
    }

    #[test]
    fn lipinski_hydrogen_counts_include_explicit_neighbors_without_double_counting_atoms() {
        let molecule = explicit_methylamine();
        assert_eq!(calc_lipinski_hba(&molecule), Ok(1));
        assert_eq!(calc_lipinski_hbd(&molecule), Ok(2));
        assert_eq!(crate::properties::descriptors::calc_num_hba(&molecule), Ok(1));
        assert_eq!(crate::properties::descriptors::calc_num_hbd(&molecule), Ok(1));
        assert_eq!(calc_num_heteroatoms(&molecule), Ok(1));
        assert_eq!(calc_num_amide_bonds(&molecule), Ok(0));
        assert_eq!(calc_num_heavy_atoms(&molecule), Ok(2));
        assert_eq!(calc_num_atoms(&molecule), Ok(7));
    }

    #[test]
    fn ring_classification_matches_pinned_rdkit_across_every_topology_family() {
        const CASES: [(&str, [u32; 11]); 11] = [
            ("CC", [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
            ("c1ccccc1", [1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0]),
            ("c1ncccc1", [1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0]),
            ("C1CCCCC1", [1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1]),
            ("C1=CCCCC1", [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0]),
            ("O1CCCCC1", [1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0]),
            ("C1CCC2CCCCC2C1", [2, 0, 0, 2, 2, 0, 0, 0, 2, 0, 2]),
            ("C1CC2CCC1C2", [2, 0, 0, 2, 2, 0, 0, 0, 2, 0, 2]),
            ("C1CCC2(CC1)CCCCC2", [2, 0, 0, 2, 2, 0, 0, 0, 2, 0, 2]),
            ("C1CCCCCCCCCCC1", [1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1]),
            ("c1ccc2ncccc2c1", [2, 1, 2, 0, 0, 1, 1, 0, 0, 0, 0]),
        ];

        for (smiles, expected) in CASES {
            let molecule = Molecule::from_smiles(smiles).expect("ring fixture");
            let actual = [
                calc_num_rings(&molecule).unwrap(),
                calc_num_heterocycles(&molecule).unwrap(),
                crate::properties::descriptors::calc_num_aromatic_rings(&molecule).unwrap(),
                calc_num_saturated_rings(&molecule).unwrap(),
                calc_num_aliphatic_rings(&molecule).unwrap(),
                calc_num_aromatic_heterocycles(&molecule).unwrap(),
                calc_num_aromatic_carbocycles(&molecule).unwrap(),
                calc_num_aliphatic_heterocycles(&molecule).unwrap(),
                calc_num_aliphatic_carbocycles(&molecule).unwrap(),
                calc_num_saturated_heterocycles(&molecule).unwrap(),
                calc_num_saturated_carbocycles(&molecule).unwrap(),
            ];
            assert_eq!(actual, expected, "{smiles:?} ring descriptor vector");
        }
    }

    #[test]
    fn spiro_and_bridgehead_atoms_match_pinned_rdkit_order_and_deduplication() {
        const CASES: [(&str, &[usize], &[usize]); 5] = [
            ("C1CC2CCC1CC2", &[], &[2, 5]),
            ("C1CCC2(C1)CC1CCC2CC1", &[3], &[6, 9]),
            ("C1CCC2(CC1)CCC1(CC2)CCCC1", &[3, 8], &[]),
            ("C1C2CC3CC1CC(C2)C3", &[], &[1, 5, 3, 7]),
            ("C1CCC2CCCCC2C1", &[], &[]),
        ];

        for (smiles, expected_spiro, expected_bridgehead) in CASES {
            let molecule = Molecule::from_smiles(smiles).expect("spiro/bridgehead fixture");
            let mut spiro_atoms = Vec::new();
            let mut bridgehead_atoms = Vec::new();

            assert_eq!(
                calc_num_spiro_atoms(&molecule, Some(&mut spiro_atoms)).unwrap(),
                u32::try_from(expected_spiro.len()).unwrap(),
                "{smiles:?} spiro count"
            );
            assert_eq!(
                spiro_atoms,
                expected_spiro.iter().copied().map(AtomId::new).collect::<Vec<_>>(),
                "{smiles:?} spiro atom order"
            );
            assert_eq!(
                calc_num_bridgehead_atoms(&molecule, Some(&mut bridgehead_atoms)).unwrap(),
                u32::try_from(expected_bridgehead.len()).unwrap(),
                "{smiles:?} bridgehead count"
            );
            assert_eq!(
                bridgehead_atoms,
                expected_bridgehead.iter().copied().map(AtomId::new).collect::<Vec<_>>(),
                "{smiles:?} bridgehead atom order"
            );
        }

        let molecule = Molecule::from_smiles("C1CCC2(CC1)CCC1(CC2)CCCC1").expect("multi-spiro fixture");
        let mut atoms = vec![AtomId::new(3)];
        assert_eq!(calc_num_spiro_atoms(&molecule, Some(&mut atoms)), Ok(2));
        assert_eq!(atoms, vec![AtomId::new(3), AtomId::new(8)]);
    }

    #[test]
    fn stereo_center_counts_match_pinned_rdkit_without_mutating_the_caller() {
        const CASES: [(&str, u32, u32); 7] = [
            ("", 0, 0),
            ("CC(F)Cl", 1, 1),
            ("C[C@H](F)Cl", 1, 0),
            ("CC(F)C", 0, 0),
            ("CC(F)Cl.C[C@@H](Br)I", 2, 1),
            ("FC(Cl)(Br)I", 1, 1),
            ("F[C@](Cl)(Br)I", 1, 0),
        ];

        for (smiles, expected_all, expected_unspecified) in CASES {
            let molecule = Molecule::from_smiles(smiles).expect("stereo-center fixture");
            let before = molecule.clone();
            assert_eq!(
                num_atom_stereo_centers(&molecule),
                Ok(expected_all),
                "{smiles:?} all stereo centers"
            );
            assert_eq!(
                num_unspecified_atom_stereo_centers(&molecule),
                Ok(expected_unspecified),
                "{smiles:?} unspecified stereo centers"
            );
            assert_eq!(
                num_atom_stereo_centers(&molecule),
                Ok(expected_all),
                "{smiles:?} repeated all stereo centers"
            );
            assert_eq!(molecule, before, "{smiles:?} descriptor caller state");
        }

        let mut preassigned = Molecule::from_smiles("C[C@H](F)Cl").expect("preassigned stereo-center fixture");
        preassigned.properties_mut().set_computed_prop("_StereochemDone", "1");
        assert!(has_stereo_assigned(&preassigned));
        let before = preassigned.clone();
        assert_eq!(num_atom_stereo_centers(&preassigned), Ok(1));
        assert_eq!(num_unspecified_atom_stereo_centers(&preassigned), Ok(0));
        assert_eq!(preassigned, before);
    }
}
