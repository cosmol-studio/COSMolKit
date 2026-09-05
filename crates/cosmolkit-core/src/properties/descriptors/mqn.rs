//! Molecular quantum number descriptor implementation boundary.
//!
//! This module exclusively owns the single fixed-order 42-component MQN core.
//! It reuses shared ring, rotatable-bond, hydrogen, and graph state while
//! retaining the source's MQN-specific N/O polarity counts in this one core.

use super::{DescriptorError, DescriptorResult, NumRotatableBondsOptions};
use crate::{BondOrder, Molecule};

pub(super) fn calc_mqns_core(molecule: &Molecule) -> DescriptorResult<[u32; 42]> {
    // BEGIN RDKIT CPP FUNCTION: RDKit::Descriptors::calcMQNs
    // RDKit✔️✔️: std::vector<unsigned int> calcMQNs(const ROMol &mol, bool) {
    // RDKit✔️✔️:   // FIX: use force value to enable caching
    // RDKit✔️✔️:   std::vector<unsigned int> res(42, 0);
    let mut result = [0_u32; 42];
    let rings =
        super::descriptor_ring_info(molecule, false).map_err(|_| DescriptorError::Unsupported {
            function: "calc_mqns",
            rdkit_function: "MolOps::symmetrizeSSSR/RingInfo",
        })?;

    // RDKit✔️✔️:   // ---------------------------------------------------
    // RDKit✔️✔️:   // atom-centered things
    // RDKit✔️✔️:   // Note: We're not doing exactly the same thing
    // RDKit✔️✔️:   //       as the original paper on polarity counts
    // RDKit✔️✔️:   //       since we're using different donor and acceptor
    // RDKit✔️✔️:   //       definitions.
    // RDKit✔️✔️:   ROMol::VERTEX_ITER atBegin, atEnd;
    // RDKit✔️✔️:   boost::tie(atBegin, atEnd) = mol.getVertices();
    // RDKit✔️✔️:   while (atBegin != atEnd) {
    for atom in molecule.atoms() {
        // RDKit✔️✔️:     const Atom *at = mol[*atBegin];
        // RDKit✔️✔️:     ++atBegin;
        // RDKit✔️✔️:     unsigned int nHs = at->getTotalNumHs();
        let hydrogen_count = crate::valence::total_num_hydrogens(molecule, atom.id(), false)
            .map_err(|_| DescriptorError::Unsupported {
                function: "calc_mqns",
                rdkit_function: "Atom::getTotalNumHs",
            })?;
        // RDKit✔️✔️:     unsigned int nRings = mol.getRingInfo()->numAtomRings(at->getIdx());
        let ring_count = rings.num_atom_rings(atom.id());
        let degree = molecule.adjacency().neighbors_of(atom.id().index()).len();
        // RDKit✔️✔️:     switch (at->getAtomicNum()) {
        match atom.atomic_number() {
            // RDKit✔️✔️:       case 0:
            // RDKit✔️✔️:       case 1:
            // RDKit✔️✔️:         break;
            0 | 1 => {}
            // RDKit✔️✔️:       case 6:
            // RDKit✔️✔️:         res[0]++;
            // RDKit✔️✔️:         break;
            6 => result[0] += 1,
            // RDKit✔️✔️:       case 9:
            // RDKit✔️✔️:         res[1]++;
            // RDKit✔️✔️:         break;
            9 => result[1] += 1,
            // RDKit✔️✔️:       case 17:
            // RDKit✔️✔️:         res[2]++;
            // RDKit✔️✔️:         break;
            17 => result[2] += 1,
            // RDKit✔️✔️:       case 35:
            // RDKit✔️✔️:         res[3]++;
            // RDKit✔️✔️:         break;
            35 => result[3] += 1,
            // RDKit✔️✔️:       case 53:
            // RDKit✔️✔️:         res[4]++;
            // RDKit✔️✔️:         break;
            53 => result[4] += 1,
            // RDKit✔️✔️:       case 16:
            // RDKit✔️✔️:         res[5]++;
            // RDKit✔️✔️:         break;
            16 => result[5] += 1,
            // RDKit✔️✔️:       case 15:
            // RDKit✔️✔️:         res[6]++;
            // RDKit✔️✔️:         break;
            15 => result[6] += 1,
            // RDKit✔️✔️:       case 7:
            7 => {
                // RDKit✔️✔️:         if (!nRings) {
                // RDKit✔️✔️:           res[7]++;
                // RDKit✔️✔️:         } else {
                // RDKit✔️✔️:           res[8]++;
                // RDKit✔️✔️:         }
                if ring_count == 0 {
                    result[7] += 1;
                } else {
                    result[8] += 1;
                }
                // RDKit✔️✔️:         if (at->getDegree() != 4) {
                // RDKit✔️✔️:           res[19]++;  // number of acceptor sites
                // RDKit✔️✔️:           res[20]++;  // number of acceptor atoms
                // RDKit✔️✔️:         }
                if degree != 4 {
                    result[19] += 1;
                    result[20] += 1;
                }
                // RDKit✔️✔️:         if (nHs) {
                // RDKit✔️✔️:           res[21] += nHs;  // number of donor sites
                // RDKit✔️✔️:           res[22]++;       // number of donor atoms
                // RDKit✔️✔️:         }
                if hydrogen_count != 0 {
                    result[21] += hydrogen_count;
                    result[22] += 1;
                }
                // RDKit✔️✔️:         break;
            }
            // RDKit✔️✔️:       case 8:
            8 => {
                // RDKit✔️✔️:         if (!nRings) {
                // RDKit✔️✔️:           res[9]++;
                // RDKit✔️✔️:         } else {
                // RDKit✔️✔️:           res[10]++;
                // RDKit✔️✔️:         }
                if ring_count == 0 {
                    result[9] += 1;
                } else {
                    result[10] += 1;
                }
                // RDKit✔️✔️:         res[20]++;  // number of acceptor atoms
                result[20] += 1;
                // RDKit✔️✔️:         if (at->getFormalCharge() != -1) {
                // RDKit✔️✔️:           res[19] += 2;  // number of acceptor sites
                // RDKit✔️✔️:         } else {
                // RDKit✔️✔️:           res[19] += 3;  // number of acceptor sites
                // RDKit✔️✔️:         }
                if atom.formal_charge() != -1 {
                    result[19] += 2;
                } else {
                    result[19] += 3;
                }
                // RDKit✔️✔️:         if (nHs) {
                // RDKit✔️✔️:           res[21] += nHs;  // number of donor sites
                // RDKit✔️✔️:           res[22]++;       // number of donor atoms
                // RDKit✔️✔️:         }
                if hydrogen_count != 0 {
                    result[21] += hydrogen_count;
                    result[22] += 1;
                }
                // RDKit✔️✔️:         break;
            }
            // RDKit✔️✔️:       default:
            // RDKit✔️✔️:         break;
            _ => {} // RDKit✔️✔️:     }
        }

        // RDKit✔️✔️:     if (at->getFormalCharge() > 0) {
        // RDKit✔️✔️:       res[24]++;  // positive charges
        // RDKit✔️✔️:     } else if (at->getFormalCharge() < 0) {
        // RDKit✔️✔️:       res[23]++;  // negative charges
        // RDKit✔️✔️:     }
        if atom.formal_charge() > 0 {
            result[24] += 1;
        } else if atom.formal_charge() < 0 {
            result[23] += 1;
        }

        // RDKit✔️✔️:     if (at->getAtomicNum() != 1) {
        if atom.atomic_number() != 1 {
            // RDKit✔️✔️:       switch (at->getDegree()) {
            match degree {
                // RDKit✔️✔️:         case 1:
                // RDKit✔️✔️:           res[25]++;
                // RDKit✔️✔️:           break;
                1 => result[25] += 1,
                // RDKit✔️✔️:         case 2:
                2 => {
                    // RDKit✔️✔️:           if (!nRings) {
                    // RDKit✔️✔️:             res[26]++;
                    // RDKit✔️✔️:           } else {
                    // RDKit✔️✔️:             res[29]++;
                    // RDKit✔️✔️:           }
                    if ring_count == 0 {
                        result[26] += 1;
                    } else {
                        result[29] += 1;
                    }
                    // RDKit✔️✔️:           break;
                }
                // RDKit✔️✔️:         case 3:
                3 => {
                    // RDKit✔️✔️:           if (!nRings) {
                    // RDKit✔️✔️:             res[27]++;
                    // RDKit✔️✔️:           } else {
                    // RDKit✔️✔️:             res[30]++;
                    // RDKit✔️✔️:           }
                    if ring_count == 0 {
                        result[27] += 1;
                    } else {
                        result[30] += 1;
                    }
                    // RDKit✔️✔️:           break;
                }
                // RDKit✔️✔️:         case 4:
                4 => {
                    // RDKit✔️✔️:           if (!nRings) {
                    // RDKit✔️✔️:             res[28]++;
                    // RDKit✔️✔️:           } else {
                    // RDKit✔️✔️:             res[31]++;
                    // RDKit✔️✔️:           }
                    if ring_count == 0 {
                        result[28] += 1;
                    } else {
                        result[31] += 1;
                    }
                    // RDKit✔️✔️:           break;
                }
                // RDKit✔️✔️:       }
                _ => {}
            }
            // RDKit✔️✔️:       if (nRings >= 2) {
            // RDKit✔️✔️:         res[40]++;
            // RDKit✔️✔️:       }
            if ring_count >= 2 {
                result[40] += 1;
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   res[11] = mol.getNumHeavyAtoms();
    result[11] = super::lipinski::calc_num_heavy_atoms(molecule)?;

    // RDKit✔️✔️:   // ---------------------------------------------------
    // RDKit✔️✔️:   // bond counts:
    // RDKit✔️✔️:   unsigned int nAromatic = 0;
    let mut aromatic_count = 0_u32;
    // RDKit✔️✔️:   ROMol::EDGE_ITER firstB, lastB;
    // RDKit✔️✔️:   boost::tie(firstB, lastB) = mol.getEdges();
    // RDKit✔️✔️:   while (firstB != lastB) {
    for bond in molecule.bonds() {
        // RDKit✔️✔️:     const Bond *bond = mol[*firstB];
        // RDKit✔️✔️:     if (bond->getIsAromatic()) {
        // RDKit✔️✔️:       ++nAromatic;
        // RDKit✔️✔️:     }
        if bond.is_aromatic() {
            aromatic_count += 1;
        }
        // RDKit✔️✔️:     unsigned int nRings = mol.getRingInfo()->numBondRings(bond->getIdx());
        let ring_count = rings.num_bond_rings(bond.id());
        // RDKit✔️✔️:     switch (bond->getBondType()) {
        match bond.order() {
            // RDKit✔️✔️:       case Bond::SINGLE:
            BondOrder::Single => {
                // RDKit✔️✔️:         if (!nRings) {
                // RDKit✔️✔️:           res[12]++;
                // RDKit✔️✔️:         } else {
                // RDKit✔️✔️:           res[15]++;
                // RDKit✔️✔️:         }
                if ring_count == 0 {
                    result[12] += 1;
                } else {
                    result[15] += 1;
                }
                // RDKit✔️✔️:         break;
            }
            // RDKit✔️✔️:       case Bond::DOUBLE:
            BondOrder::Double => {
                // RDKit✔️✔️:         if (!nRings) {
                // RDKit✔️✔️:           res[13]++;
                // RDKit✔️✔️:         } else {
                // RDKit✔️✔️:           res[16]++;
                // RDKit✔️✔️:         }
                if ring_count == 0 {
                    result[13] += 1;
                } else {
                    result[16] += 1;
                }
                // RDKit✔️✔️:         break;
            }
            // RDKit✔️✔️:       case Bond::TRIPLE:
            BondOrder::Triple => {
                // RDKit✔️✔️:         if (!nRings) {
                // RDKit✔️✔️:           res[14]++;
                // RDKit✔️✔️:         } else {
                // RDKit✔️✔️:           res[17]++;
                // RDKit✔️✔️:         }
                if ring_count == 0 {
                    result[14] += 1;
                } else {
                    result[17] += 1;
                }
                // RDKit✔️✔️:         break;
            }
            // RDKit✔️✔️:       default:
            // RDKit✔️✔️:         break;
            _ => {} // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:     if (nRings >= 2) {
        // RDKit✔️✔️:       res[41]++;
        // RDKit✔️✔️:     }
        if ring_count >= 2 {
            result[41] += 1;
        }
        // RDKit✔️✔️:     ++firstB;
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   // rather than do the work to kekulize the molecule, we cheat
    // RDKit✔️✔️:   // by just dividing the number of aromatic bonds evenly among the
    // RDKit✔️✔️:   // cyclic single bond and cyclic double bond bins and give any
    // RDKit✔️✔️:   // remainder to the single bonds
    // RDKit✔️✔️:   res[15] += nAromatic / 2;
    result[15] += aromatic_count / 2;
    // RDKit✔️✔️:   res[16] += nAromatic / 2;
    result[16] += aromatic_count / 2;
    // RDKit✔️✔️:   if (nAromatic % 2) {
    // RDKit✔️✔️:     res[15]++;
    // RDKit✔️✔️:   }
    if aromatic_count % 2 != 0 {
        result[15] += 1;
    }
    // RDKit✔️✔️:   res[18] = calcNumRotatableBonds(mol);
    result[18] = super::calc_num_rotatable_bonds(molecule, NumRotatableBondsOptions::Default)?;

    // RDKit✔️✔️:   // ---------------------------------------------------
    // RDKit✔️✔️:   //  ring size counts
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->atomRings()) {
    for ring in rings.atom_rings() {
        // RDKit✔️✔️:     if (iv.size() < 10) {
        // RDKit✔️✔️:       res[iv.size() + 29]++;
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       res[39]++;
        // RDKit✔️✔️:     }
        if ring.len() < 10 {
            result[ring.len() + 29] += 1;
        } else {
            result[39] += 1;
        }
        // RDKit✔️✔️:   }
    }

    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: RDKit::Descriptors::calcMQNs
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::calc_mqns_core;
    use crate::Molecule;

    #[test]
    fn mqn_matches_pinned_rdkit_complete_vectors_for_focused_branches() {
        const CASES: [(&str, [u32; 42]); 5] = [
            ("", [0; 42]),
            (
                "CCO",
                [
                    2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 3, 2, 0, 0, 0, 0, 0, 0, 2, 1, 1, 1, 0, 0, 2,
                    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                ],
            ),
            (
                "c1ccccc1",
                [
                    6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 6, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                ],
            ),
            (
                "[Na+].[Cl-]",
                [
                    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                ],
            ),
            (
                "C1CC2CCC1C2",
                [
                    7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 5, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 3, 2,
                ],
            ),
        ];

        for (smiles, expected) in CASES {
            let molecule = Molecule::from_smiles(smiles).expect("focused MQN fixture");
            assert_eq!(
                calc_mqns_core(&molecule),
                Ok(expected),
                "{smiles:?} MQN vector"
            );
        }
    }

    #[test]
    fn mqn_complete_branch_fixture_exercises_and_matches_every_index() {
        const SMILES: &str = concat!(
            "CC(F)(Cl)Br.CPI.N1CCOCC1.[O-]C(=O)[NH3+].c1ncccc1.C#CC=C.",
            "C1#CC1.C1CCC2(CC1)CCCCC2.C1C2CC3CC1CC(C2)C3.C1CC1.C1CCC1.",
            "C1CCCC1.C1CCCCCC1.C1CCCCCCC1.C1CCCCCCCC1.C1CCCCCCCCCC1.CCCC.CS"
        );
        const EXPECTED: [u32; 42] = [
            93, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 105, 13, 2, 1, 82, 3, 1, 1, 10, 6, 4, 2, 1, 1, 15, 5,
            1, 1, 78, 4, 1, 2, 1, 1, 8, 1, 1, 1, 1, 11, 12,
        ];

        assert!(EXPECTED.iter().all(|&value| value != 0));
        let molecule = Molecule::from_smiles(SMILES).expect("complete MQN branch fixture");
        assert_eq!(calc_mqns_core(&molecule), Ok(EXPECTED));
    }
}
