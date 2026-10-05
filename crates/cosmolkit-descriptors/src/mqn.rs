//! Molecular Quantum Numbers over borrowed, prepared detached input.
//! Pinned RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8 (BSD), MQN.cpp.

use cosmolkit_core::total_hydrogen_count_from_validated;
use cosmolkit_model::BondOrder;

use crate::{DescriptorError, DescriptorInput, DescriptorResult, RotatableBondsOptions};

/// MQN version exported by the pinned source header.
pub const MQN_VERSION: &str = "1.0.0";

/// Computes all forty-two source MQNs from supplied FINAL valence/ring rows.
///
/// `force` is ignored by the pinned implementation. Neither value triggers
/// chemistry preparation, computed-state writes or persistent memoization.
/// The result preserves unsigned source arithmetic, component order and all
/// explicit atom/bond rows, including wildcard and explicit hydrogen rows.
/// Invalid prepared input and dependency errors retain structural causes.
pub fn mqns(input: &DescriptorInput<'_>, force: bool) -> DescriptorResult<Vec<u32>> {
    // Complete verbatim source function, MQN.cpp:20 onward:
    // RDKit✔️❌: std::vector<unsigned int> calcMQNs(const ROMol &mol, bool) {
    // RDKit✔️❌:   // FIX: use force value to enable caching
    // RDKit✔️❌:   std::vector<unsigned int> res(42, 0);
    // RDKit✔️❌:
    // RDKit✔️❌:   // ---------------------------------------------------
    // RDKit✔️❌:   // atom-centered things
    // RDKit✔️❌:   // Note: We're not doing exactly the same thing
    // RDKit✔️❌:   //       as the original paper on polarity counts
    // RDKit✔️❌:   //       since we're using different donor and acceptor
    // RDKit✔️❌:   //       definitions.
    // RDKit✔️❌:   ROMol::VERTEX_ITER atBegin, atEnd;
    // RDKit✔️❌:   boost::tie(atBegin, atEnd) = mol.getVertices();
    // RDKit✔️❌:   while (atBegin != atEnd) {
    // RDKit✔️❌:     const Atom *at = mol[*atBegin];
    // RDKit✔️❌:     ++atBegin;
    // RDKit✔️❌:     unsigned int nHs = at->getTotalNumHs();
    // RDKit✔️❌:     unsigned int nRings = mol.getRingInfo()->numAtomRings(at->getIdx());
    // RDKit✔️❌:     switch (at->getAtomicNum()) {
    // RDKit✔️❌:       case 0:
    // RDKit✔️❌:       case 1:
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       case 6:
    // RDKit✔️❌:         res[0]++;
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       case 9:
    // RDKit✔️❌:         res[1]++;
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       case 17:
    // RDKit✔️❌:         res[2]++;
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       case 35:
    // RDKit✔️❌:         res[3]++;
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       case 53:
    // RDKit✔️❌:         res[4]++;
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       case 16:
    // RDKit✔️❌:         res[5]++;
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       case 15:
    // RDKit✔️❌:         res[6]++;
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       case 7:
    // RDKit✔️❌:         if (!nRings) {
    // RDKit✔️❌:           res[7]++;
    // RDKit✔️❌:         } else {
    // RDKit✔️❌:           res[8]++;
    // RDKit✔️❌:         }
    // RDKit✔️❌:         if (at->getDegree() != 4) {
    // RDKit✔️❌:           res[19]++;  // number of acceptor sites
    // RDKit✔️❌:           res[20]++;  // number of acceptor atoms
    // RDKit✔️❌:         }
    // RDKit✔️❌:         if (nHs) {
    // RDKit✔️❌:           res[21] += nHs;  // number of donor sites
    // RDKit✔️❌:           res[22]++;       // number of donor atoms
    // RDKit✔️❌:         }
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       case 8:
    // RDKit✔️❌:         if (!nRings) {
    // RDKit✔️❌:           res[9]++;
    // RDKit✔️❌:         } else {
    // RDKit✔️❌:           res[10]++;
    // RDKit✔️❌:         }
    // RDKit✔️❌:         res[20]++;  // number of acceptor atoms
    // RDKit✔️❌:         if (at->getFormalCharge() != -1) {
    // RDKit✔️❌:           res[19] += 2;  // number of acceptor sites
    // RDKit✔️❌:         } else {
    // RDKit✔️❌:           res[19] += 3;  // number of acceptor sites
    // RDKit✔️❌:         }
    // RDKit✔️❌:         if (nHs) {
    // RDKit✔️❌:           res[21] += nHs;  // number of donor sites
    // RDKit✔️❌:           res[22]++;       // number of donor atoms
    // RDKit✔️❌:         }
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       default:
    // RDKit✔️❌:         break;
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     if (at->getFormalCharge() > 0) {
    // RDKit✔️❌:       res[24]++;  // positive charges
    // RDKit✔️❌:     } else if (at->getFormalCharge() < 0) {
    // RDKit✔️❌:       res[23]++;  // negative charges
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     if (at->getAtomicNum() != 1) {
    // RDKit✔️❌:       switch (at->getDegree()) {
    // RDKit✔️❌:         case 1:
    // RDKit✔️❌:           res[25]++;
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         case 2:
    // RDKit✔️❌:           if (!nRings) {
    // RDKit✔️❌:             res[26]++;
    // RDKit✔️❌:           } else {
    // RDKit✔️❌:             res[29]++;
    // RDKit✔️❌:           }
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         case 3:
    // RDKit✔️❌:           if (!nRings) {
    // RDKit✔️❌:             res[27]++;
    // RDKit✔️❌:           } else {
    // RDKit✔️❌:             res[30]++;
    // RDKit✔️❌:           }
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         case 4:
    // RDKit✔️❌:           if (!nRings) {
    // RDKit✔️❌:             res[28]++;
    // RDKit✔️❌:           } else {
    // RDKit✔️❌:             res[31]++;
    // RDKit✔️❌:           }
    // RDKit✔️❌:           break;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       if (nRings >= 2) {
    // RDKit✔️❌:         res[40]++;
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   res[11] = mol.getNumHeavyAtoms();
    // RDKit✔️❌:
    // RDKit✔️❌:   // ---------------------------------------------------
    // RDKit✔️❌:   // bond counts:
    // RDKit✔️❌:   unsigned int nAromatic = 0;
    // RDKit✔️❌:   ROMol::EDGE_ITER firstB, lastB;
    // RDKit✔️❌:   boost::tie(firstB, lastB) = mol.getEdges();
    // RDKit✔️❌:   while (firstB != lastB) {
    // RDKit✔️❌:     const Bond *bond = mol[*firstB];
    // RDKit✔️❌:     if (bond->getIsAromatic()) {
    // RDKit✔️❌:       ++nAromatic;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     unsigned int nRings = mol.getRingInfo()->numBondRings(bond->getIdx());
    // RDKit✔️❌:     switch (bond->getBondType()) {
    // RDKit✔️❌:       case Bond::SINGLE:
    // RDKit✔️❌:         if (!nRings) {
    // RDKit✔️❌:           res[12]++;
    // RDKit✔️❌:         } else {
    // RDKit✔️❌:           res[15]++;
    // RDKit✔️❌:         }
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       case Bond::DOUBLE:
    // RDKit✔️❌:         if (!nRings) {
    // RDKit✔️❌:           res[13]++;
    // RDKit✔️❌:         } else {
    // RDKit✔️❌:           res[16]++;
    // RDKit✔️❌:         }
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       case Bond::TRIPLE:
    // RDKit✔️❌:         if (!nRings) {
    // RDKit✔️❌:           res[14]++;
    // RDKit✔️❌:         } else {
    // RDKit✔️❌:           res[17]++;
    // RDKit✔️❌:         }
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       default:
    // RDKit✔️❌:         break;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (nRings >= 2) {
    // RDKit✔️❌:       res[41]++;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     ++firstB;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   // rather than do the work to kekulize the molecule, we cheat
    // RDKit✔️❌:   // by just dividing the number of aromatic bonds evenly among the
    // RDKit✔️❌:   // cyclic single bond and cyclic double bond bins and give any
    // RDKit✔️❌:   // remainder to the single bonds
    // RDKit✔️❌:   res[15] += nAromatic / 2;
    // RDKit✔️❌:   res[16] += nAromatic / 2;
    // RDKit✔️❌:   if (nAromatic % 2) {
    // RDKit✔️❌:     res[15]++;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   res[18] = calcNumRotatableBonds(mol);
    // RDKit✔️❌:
    // RDKit✔️❌:   // ---------------------------------------------------
    // RDKit✔️❌:   //  ring size counts
    // RDKit✔️❌:   for (const auto &iv : mol.getRingInfo()->atomRings()) {
    // RDKit✔️❌:     if (iv.size() < 10) {
    // RDKit✔️❌:       res[iv.size() + 29]++;
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       res[39]++;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   return res;
    // RDKit✔️❌: }
    // Behavior: directly reproduce source predicates and wrapping u32 work;
    // reuse the sole hydrogen, heavy-count and Default rotor owners. Each
    // supplied ring row contributes once, even duplicates or small rows.
    // Performance: the same fixed Vec(42) and atom/bond/ring passes, but
    // preflight, heavy-count and rotor boundaries validate graph/rows again.
    // Their adjacency allocations/scans are an explicit extra boundary cost.
    // Context and valence borrow FINAL rows; no chemistry or input clones.
    let _ = force;
    let topology = input.topology();
    let valence = crate::prepared_valence(topology, Some(input.valence()), "mqns")?;
    let _validated = crate::patterns::prepared_context(input, "mqns")?;
    let rings = input.ring_info();
    let mut res = vec![0u32; 42];

    for atom in &topology.atoms {
        let n_hs = total_hydrogen_count_from_validated(topology, &valence, atom.id(), false)
            .map_err(|source| DescriptorError::Valence {
                function: "mqns",
                source,
            })?;
        // Source numAtomRings returns unsigned int, not size_t.
        let n_rings = rings.num_atom_rings(atom.id()) as u32;
        // Atom::getDegree is also unsigned int in the source.
        let degree = topology.adjacency.neighbors_of(atom.id().index()).len() as u32;
        match atom.atomic_number() {
            6 => res[0] = res[0].wrapping_add(1),
            9 => res[1] = res[1].wrapping_add(1),
            17 => res[2] = res[2].wrapping_add(1),
            35 => res[3] = res[3].wrapping_add(1),
            53 => res[4] = res[4].wrapping_add(1),
            16 => res[5] = res[5].wrapping_add(1),
            15 => res[6] = res[6].wrapping_add(1),
            7 => {
                let bin = if n_rings == 0 { 7 } else { 8 };
                res[bin] = res[bin].wrapping_add(1);
                if degree != 4 {
                    res[19] = res[19].wrapping_add(1);
                    res[20] = res[20].wrapping_add(1);
                }
                if n_hs != 0 {
                    res[21] = res[21].wrapping_add(n_hs);
                    res[22] = res[22].wrapping_add(1);
                }
            }
            8 => {
                let bin = if n_rings == 0 { 9 } else { 10 };
                res[bin] = res[bin].wrapping_add(1);
                res[20] = res[20].wrapping_add(1);
                res[19] = res[19].wrapping_add(if atom.formal_charge() != -1 { 2 } else { 3 });
                if n_hs != 0 {
                    res[21] = res[21].wrapping_add(n_hs);
                    res[22] = res[22].wrapping_add(1);
                }
            }
            _ => {}
        }
        if atom.formal_charge() > 0 {
            res[24] = res[24].wrapping_add(1);
        } else if atom.formal_charge() < 0 {
            res[23] = res[23].wrapping_add(1);
        }
        if atom.atomic_number() != 1 {
            let bin = match degree {
                1 => Some(25),
                2 if n_rings == 0 => Some(26),
                2 => Some(29),
                3 if n_rings == 0 => Some(27),
                3 => Some(30),
                4 if n_rings == 0 => Some(28),
                4 => Some(31),
                _ => None,
            };
            if let Some(bin) = bin {
                res[bin] = res[bin].wrapping_add(1);
            }
            if n_rings >= 2 {
                res[40] = res[40].wrapping_add(1);
            }
        }
    }
    // Preflight already validated this immutable topology. The count kernel's
    // remaining known flattened failure is its checked source-u32 overflow.
    // Preserve other dependency errors, never invent Unsupported for MQNs.
    res[11] = crate::num_heavy_atoms_prepared(input).map_err(|error| match error {
        DescriptorError::Unsupported {
            function: "num_heavy_atoms",
            ..
        } => DescriptorError::CountOverflow {
            function: "mqns",
            field: "heavy_atoms",
        },
        other => other,
    })?;

    let mut n_aromatic = 0u32;
    for bond in &topology.bonds {
        if bond.is_aromatic() {
            n_aromatic = n_aromatic.wrapping_add(1);
        }
        let n_rings = rings.num_bond_rings(bond.id()) as u32;
        let bin = match (bond.order(), n_rings == 0) {
            (BondOrder::Single, true) => Some(12),
            (BondOrder::Single, false) => Some(15),
            (BondOrder::Double, true) => Some(13),
            (BondOrder::Double, false) => Some(16),
            (BondOrder::Triple, true) => Some(14),
            (BondOrder::Triple, false) => Some(17),
            _ => None,
        };
        if let Some(bin) = bin {
            res[bin] = res[bin].wrapping_add(1);
        }
        if n_rings >= 2 {
            res[41] = res[41].wrapping_add(1);
        }
    }
    res[15] = res[15].wrapping_add(n_aromatic / 2);
    res[16] = res[16].wrapping_add(n_aromatic / 2);
    if n_aromatic % 2 != 0 {
        res[15] = res[15].wrapping_add(1);
    }
    res[18] = crate::num_rotatable_bonds_prepared(input, RotatableBondsOptions::Default)?;

    for atoms in rings.atom_rings() {
        // Unlike the membership-count getter, source iv.size() is size_t.
        let bin = if atoms.len() < 10 {
            atoms.len() + 29
        } else {
            39
        };
        res[bin] = res[bin].wrapping_add(1);
    }
    Ok(res)
}
