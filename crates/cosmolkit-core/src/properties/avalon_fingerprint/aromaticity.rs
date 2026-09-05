//! Source-shaped Avalon aromatic-bond perception used by fingerprint preprocessing.

use crate::FingerprintError;

use super::reaccs::MoleculeState;
use super::rings::{combine_rings, proper_ring_pairs, ring_list};

const SINGLE: i32 = 1;
const DOUBLE: i32 = 2;
const TRIPLE: i32 = 3;
const AROMATIC: i32 = 4;

pub(super) fn perceive_aromatic_bonds(molecule: &mut MoleculeState) -> Result<(), FingerprintError> {
    perceive_aromatic_bonds_with(molecule, |_| {})
}

fn perceive_aromatic_bonds_with(
    molecule: &mut MoleculeState,
    mut record_mutation: impl FnMut(usize),
) -> Result<(), FingerprintError> {
    // Avalon❗✔️: void PerceiveAromaticBonds(struct reaccs_molecule_t *mp)
    // Avalon❗✔️: /*
    // Avalon❗✔️:  * Sets bond types of the bonds of *mp to "AROMATIC" if in
    // Avalon❗✔️:  * six-ring of sp2 atoms.
    // Avalon❗✔️:  */
    // Avalon❗✔️: {
    // Avalon❗✔️:    unsigned i;
    // Avalon❗✔️:    int ii;
    // Avalon❗✔️:    int changed;
    // Avalon❗✔️:    int ndouble, nsingle, naromatic;
    // Avalon❗✔️:    bond_set_node *ring_list, *ring_pairs, *plist, *plisth;
    // Avalon❗✔️:    pair_t *graph;
    // Avalon❗✔️:    int *bond_is_in_ring;
    // Avalon❗✔️:    int *sp_count;
    // Avalon❗✔️:    int is_cumulene;
    // Avalon❗✔️:    int ring_is_aromatic;
    // Avalon❗✔️:
    // Avalon❗✔️:    graph = TypeAlloc(mp->n_bonds+1, pair_t);
    let mut graph = Vec::with_capacity(molecule.bonds.len());
    // Avalon❗✔️:    sp_count = TypeAlloc(mp->n_atoms+1, int);
    let mut sp_count = vec![0_i32; molecule.atoms.len() + 1];
    // Avalon❗✔️:    bond_is_in_ring = TypeAlloc(mp->n_bonds, int);
    let mut bond_is_in_ring = vec![false; molecule.bonds.len()];
    // Avalon❗✔️:
    // Avalon❗✔️:    for (i=0; i<mp->n_bonds; i++)
    // Avalon❗✔️:    {
    for (bond_index, bond) in molecule.bonds.iter().enumerate() {
        // Avalon❗✔️:       bond_is_in_ring[i] = FALSE;
        // The zero-initialized Rust vector already holds FALSE.
        // Avalon❗✔️:       graph[i][0] = mp->bond_array[i].atoms[0];
        let first = validated_source_atom_number(bond.atoms[0], molecule.atoms.len(), bond_index)?;
        // Avalon❗✔️:       graph[i][1] = mp->bond_array[i].atoms[1];
        let second = validated_source_atom_number(bond.atoms[1], molecule.atoms.len(), bond_index)?;
        graph.push([first, second]);
        // Avalon❗✔️:    }
    }
    // Avalon❗✔️:    ring_list = RingList(graph,mp->n_bonds);
    let mut rings = ring_list(&graph);
    // Avalon❗✔️:    ring_list = CombineRings(ring_list);
    combine_rings(&mut rings);
    // Avalon❗✔️:
    // Avalon❗✔️:    for (plist=ring_list; plist!=(bond_set_node *)NULL; plist=plist->next)
    // Avalon❗✔️:       for (ii=NextMember(plist->bond_set, 0); ii>=0; ii=NextMember(plist->bond_set, (unsigned)ii+1))
    // Avalon❗✔️:          bond_is_in_ring[ii] = TRUE;
    for ring in &rings {
        for (bond_index, in_ring) in bond_is_in_ring.iter_mut().enumerate() {
            if ring.bond_set.contains(bond_index) {
                *in_ring = true;
            }
        }
    }
    // Avalon❗✔️:
    // Avalon❗✔️:    /* add the additional rings for the analysis */
    // Avalon❗✔️:    ring_pairs = ProperRingPairs(ring_list, mp->n_atoms, graph);
    let mut ring_pairs = proper_ring_pairs(&rings, molecule.atoms.len(), &graph);
    // Avalon❗✔️:    while (ring_pairs!=(bond_set_node *)NULL)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       plisth = ring_pairs->next;
    // Avalon❗✔️:       ring_pairs->next = ring_list;
    // Avalon❗✔️:       ring_list=ring_pairs;
    // Avalon❗✔️:       ring_pairs = plisth;
    // Avalon❗✔️:    }
    // ProperRingPairs returns prepend order. Moving each head to the front of
    // the basis list reverses it again, yielding pair-enumeration order first.
    ring_pairs.reverse();
    ring_pairs.append(&mut rings);
    rings = ring_pairs;
    // Avalon❗✔️:
    // Avalon❗✔️:    do
    // Avalon❗✔️:    {
    loop {
        // Avalon❗✔️:       changed = FALSE;
        let mut changed = false;
        // Avalon❗✔️:       for (plist=ring_list; plist!=(bond_set_node *)NULL; plist=plist->next)
        // Avalon❗✔️:       {
        for ring in &rings {
            // Avalon❗✔️:          for (i=0; i<=mp->n_atoms; i++) sp_count[i]=0;
            sp_count.fill(0);
            // Avalon❗✔️:          ndouble = nsingle = naromatic = 0;
            let mut ndouble = 0_i32;
            let mut nsingle = 0_i32;
            let mut naromatic = 0_i32;
            // Avalon❗✔️:          for (ii=NextMember(plist->bond_set, 0); ii>=0; ii=NextMember(plist->bond_set, (unsigned)ii+1))
            // Avalon❗✔️:          {
            for bond_index in 0..=ring.bond_set.max_member() {
                if !ring.bond_set.contains(bond_index) {
                    continue;
                }
                // Avalon❗✔️:             i = (unsigned)ii;
                // Avalon❗✔️:                if (mp->bond_array[i].bond_type == SINGLE)        nsingle++;
                if molecule.bonds[bond_index].bond_type == SINGLE {
                    nsingle += 1;
                // Avalon❗✔️:                // else if (mp->bond_array[i].bond_type == DOUBLE)   ndouble++;
                // Avalon❗✔️:                else if (mp->bond_array[i].bond_type == AROMATIC) naromatic++;
                } else if molecule.bonds[bond_index].bond_type == AROMATIC {
                    naromatic += 1;
                }
                // Avalon❗✔️:          }
            }
            // Avalon❗✔️:
            // Avalon❗✔️:          /* count number of double bonds per ring atom and mark good ones as sp2 */
            // Avalon❗✔️:          is_cumulene = FALSE;
            let mut is_cumulene = false;
            // Avalon❗✔️:          for (ii=NextMember(plist->bond_set, 0); ii>=0; ii=NextMember(plist->bond_set, (unsigned)ii+1))
            // Avalon❗✔️:          {
            for bond_index in 0..=ring.bond_set.max_member() {
                if !ring.bond_set.contains(bond_index) {
                    continue;
                }
                // Avalon❗✔️:             i = (unsigned)ii;
                // Avalon❗✔️:             if (mp->bond_array[i].bond_type == DOUBLE)
                // Avalon❗✔️:             {
                let bond = &molecule.bonds[bond_index];
                if bond.bond_type == DOUBLE {
                    // Avalon❗✔️:                ndouble++;
                    ndouble += 1;
                    // Avalon❗✔️:                sp_count[mp->bond_array[i].atoms[0]]++;
                    sp_count[bond.atoms[0] as usize] += 1;
                    // Avalon❗✔️:                if (sp_count[mp->bond_array[i].atoms[0]] > 1) is_cumulene = TRUE;
                    if sp_count[bond.atoms[0] as usize] > 1 {
                        is_cumulene = true;
                    }
                    // Avalon❗✔️:                sp_count[mp->bond_array[i].atoms[1]]++;
                    sp_count[bond.atoms[1] as usize] += 1;
                    // Avalon❗✔️:                if (sp_count[mp->bond_array[i].atoms[1]] > 1) is_cumulene = TRUE;
                    if sp_count[bond.atoms[1] as usize] > 1 {
                        is_cumulene = true;
                    }
                    // Avalon❗✔️:             }
                    // Avalon❗✔️:             else if (mp->bond_array[i].bond_type == TRIPLE)
                    // Avalon❗✔️:             {
                } else if bond.bond_type == TRIPLE {
                    // Avalon❗✔️:                is_cumulene = TRUE;
                    is_cumulene = true;
                    // Avalon❗✔️:             }
                }
                // Avalon❗✔️:          }
            }
            // Avalon❗✔️:          // mark proven aromatic atoms as sp2
            // Avalon❗✔️:          for (ii=NextMember(plist->bond_set, 0); ii>=0; ii=NextMember(plist->bond_set, (unsigned)ii+1))
            // Avalon❗✔️:          {
            for bond_index in 0..=ring.bond_set.max_member() {
                if !ring.bond_set.contains(bond_index) {
                    continue;
                }
                // Avalon❗✔️:             i = (unsigned)ii;
                // Avalon❗✔️:             if (mp->bond_array[i].bond_type == AROMATIC  &&  IsMember(plist->bond_set,i))
                // Avalon❗✔️:             {
                let bond = &molecule.bonds[bond_index];
                if bond.bond_type == AROMATIC && ring.bond_set.contains(bond_index) {
                    // Avalon❗✔️:                if (sp_count[mp->bond_array[i].atoms[0]] == 0) sp_count[mp->bond_array[i].atoms[0]]++;
                    if sp_count[bond.atoms[0] as usize] == 0 {
                        sp_count[bond.atoms[0] as usize] += 1;
                    }
                    // Avalon❗✔️:                if (sp_count[mp->bond_array[i].atoms[1]] == 0) sp_count[mp->bond_array[i].atoms[1]]++;
                    if sp_count[bond.atoms[1] as usize] == 0 {
                        sp_count[bond.atoms[1] as usize] += 1;
                    }
                    // Avalon❗✔️:             }
                }
                // Avalon❗✔️:          }
            }
            // Avalon❗✔️:          // check if all ring atoms are sp2 */
            // Avalon❗✔️:          ring_is_aromatic =
            // Avalon❗✔️:             !is_cumulene  &&                    // cumulenes are not good
            // Avalon❗✔️:             ((plist->cardinality-2)%4) == 0;    // Hueckle aromatic ring sizes
            let mut ring_is_aromatic = !is_cumulene && ring.cardinality.wrapping_sub(2) % 4 == 0;
            // Avalon❗✔️:          for (ii=NextMember(plist->bond_set, 0); ii>=0; ii=NextMember(plist->bond_set, (unsigned)ii+1))
            // Avalon❗✔️:          {
            for bond_index in 0..=ring.bond_set.max_member() {
                if !ring.bond_set.contains(bond_index) {
                    continue;
                }
                // Avalon❗✔️:             i = (unsigned)ii;
                let bond = &molecule.bonds[bond_index];
                // Avalon❗✔️:             /* all bonds must connect sp2 atoms */
                // Avalon❗✔️:             if (sp_count[mp->bond_array[i].atoms[0]] != 1) ring_is_aromatic = FALSE;
                if sp_count[bond.atoms[0] as usize] != 1 {
                    ring_is_aromatic = false;
                }
                // Avalon❗✔️:             if (sp_count[mp->bond_array[i].atoms[1]] != 1) ring_is_aromatic = FALSE;
                if sp_count[bond.atoms[1] as usize] != 1 {
                    ring_is_aromatic = false;
                }
                // Avalon❗✔️:          }
            }
            // Avalon❗✔️:
            // Avalon❗✔️:          if (ring_is_aromatic  &&  (ndouble > 0  ||  nsingle > 0))
            // Avalon❗✔️:          {
            if ring_is_aromatic && (ndouble > 0 || nsingle > 0) {
                // Avalon❗✔️:             for (ii=NextMember(plist->bond_set, 0); ii>=0; ii=NextMember(plist->bond_set, (unsigned)ii+1))
                // Avalon❗✔️:             {
                for (bond_index, in_basis_ring) in bond_is_in_ring.iter().copied().enumerate() {
                    if !ring.bond_set.contains(bond_index) {
                        continue;
                    }
                    // Avalon❗✔️:                i = (unsigned)ii;
                    // Avalon❗✔️:                if (bond_is_in_ring[i]  &&  mp->bond_array[i].bond_type != AROMATIC)
                    // Avalon❗✔️:                {
                    if in_basis_ring && molecule.bonds[bond_index].bond_type != AROMATIC {
                        // Avalon❗✔️:                   changed = TRUE;
                        changed = true;
                        // Avalon❗✔️:                   mp->bond_array[i].bond_type = AROMATIC;
                        molecule.bonds[bond_index].bond_type = AROMATIC;
                        record_mutation(bond_index);
                        // Avalon❗✔️:                }
                    }
                    // Avalon❗✔️:             }
                }
                // Avalon❗✔️:          }
            }
            // Avalon❗✔️:       }
        }
        // Avalon❗✔️:    } while (changed);
        if !changed {
            break;
        }
    }
    // Avalon❗✔️:
    // Avalon❗✔️:    DisposeBondSetList(ring_list);
    // Avalon❗✔️:
    // Avalon❗✔️:    MyFree((char *)bond_is_in_ring);
    // Avalon❗✔️:    MyFree((char *)sp_count);
    // Avalon❗✔️:    MyFree((char *)graph);
    // Rust releases the ring list and all three working vectors on scope exit.
    // Avalon❗✔️: }
    Ok(())
}

fn validated_source_atom_number(
    atom_number: i32,
    atom_count: usize,
    bond_index: usize,
) -> Result<usize, FingerprintError> {
    usize::try_from(atom_number)
        .ok()
        .filter(|&number| number > 0 && number <= atom_count)
        .ok_or_else(|| FingerprintError::AvalonConversion {
            reason: format!("Avalon bond {} references invalid atom {}", bond_index + 1, atom_number),
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::properties::avalon_fingerprint::reaccs::{Atom, Bond};

    fn molecule(atom_count: usize, endpoints: &[[i32; 2]], types: &[i32]) -> MoleculeState {
        MoleculeState {
            atoms: vec![Atom::default(); atom_count],
            bonds: endpoints
                .iter()
                .zip(types)
                .map(|(&atoms, &bond_type)| Bond {
                    atoms,
                    bond_type,
                    ..Bond::default()
                })
                .collect(),
            ..MoleculeState::default()
        }
    }

    fn bond_types(molecule: &MoleculeState) -> Vec<i32> {
        molecule.bonds.iter().map(|bond| bond.bond_type).collect()
    }

    #[test]
    fn alternating_six_ring_matches_native_mutation_order() {
        let endpoints = [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1]];
        let mut molecule = molecule(6, &endpoints, &[SINGLE, DOUBLE, SINGLE, DOUBLE, SINGLE, DOUBLE]);
        let mut mutations = Vec::new();

        perceive_aromatic_bonds_with(&mut molecule, |bond| mutations.push(bond)).unwrap();

        assert_eq!(bond_types(&molecule), vec![AROMATIC; 6]);
        assert_eq!(mutations, vec![0, 1, 2, 3, 4, 5]);
    }

    #[test]
    fn four_ring_cumulene_and_triple_cases_match_native_exclusions() {
        let four_endpoints = [[1, 2], [2, 3], [3, 4], [4, 1]];
        let four_types = [SINGLE, DOUBLE, SINGLE, DOUBLE];
        let mut four = molecule(4, &four_endpoints, &four_types);
        perceive_aromatic_bonds(&mut four).unwrap();
        assert_eq!(bond_types(&four), four_types);

        let six_endpoints = [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1]];
        let cumulene_types = [DOUBLE, DOUBLE, SINGLE, DOUBLE, SINGLE, SINGLE];
        let mut cumulene = molecule(6, &six_endpoints, &cumulene_types);
        perceive_aromatic_bonds(&mut cumulene).unwrap();
        assert_eq!(bond_types(&cumulene), cumulene_types);

        let triple_types = [TRIPLE, SINGLE, SINGLE, DOUBLE, SINGLE, DOUBLE];
        let mut triple = molecule(6, &six_endpoints, &triple_types);
        perceive_aromatic_bonds(&mut triple).unwrap();
        assert_eq!(bond_types(&triple), triple_types);
    }

    #[test]
    fn existing_aromatic_bonds_supply_source_sp2_marks() {
        let endpoints = [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1]];
        let mut molecule = molecule(6, &endpoints, &[AROMATIC, AROMATIC, SINGLE, DOUBLE, SINGLE, DOUBLE]);
        let mut mutations = Vec::new();

        perceive_aromatic_bonds_with(&mut molecule, |bond| mutations.push(bond)).unwrap();

        assert_eq!(bond_types(&molecule), vec![AROMATIC; 6]);
        assert_eq!(mutations, vec![2, 3, 4, 5]);
    }

    #[test]
    fn fused_proper_ring_is_processed_before_basis_rings() {
        let endpoints = [[1, 2], [2, 3], [3, 4], [4, 1], [3, 5], [5, 6], [6, 4]];
        let mut molecule = molecule(6, &endpoints, &[DOUBLE, SINGLE, SINGLE, SINGLE, DOUBLE, SINGLE, DOUBLE]);
        let mut mutations = Vec::new();

        perceive_aromatic_bonds_with(&mut molecule, |bond| mutations.push(bond)).unwrap();

        assert_eq!(
            bond_types(&molecule),
            vec![AROMATIC, AROMATIC, SINGLE, AROMATIC, AROMATIC, AROMATIC, AROMATIC]
        );
        assert_eq!(mutations, vec![0, 1, 3, 4, 5, 6]);
    }

    #[test]
    fn proper_ring_order_requires_and_preserves_repeated_propagation() {
        let endpoints = [
            [1, 2],
            [2, 3],
            [3, 4],
            [4, 1],
            [3, 5],
            [5, 6],
            [6, 4],
            [5, 7],
            [7, 8],
            [8, 6],
        ];
        let mut molecule = molecule(
            8,
            &endpoints,
            &[
                SINGLE, DOUBLE, SINGLE, DOUBLE, SINGLE, DOUBLE, SINGLE, SINGLE, DOUBLE, SINGLE,
            ],
        );
        let mut mutations = Vec::new();

        perceive_aromatic_bonds_with(&mut molecule, |bond| mutations.push(bond)).unwrap();

        assert_eq!(bond_types(&molecule), vec![AROMATIC; 10]);
        assert_eq!(mutations, vec![0, 1, 3, 4, 5, 6, 2, 7, 8, 9]);
    }

    #[test]
    fn invalid_one_based_bond_endpoint_is_structured_error() {
        for invalid in [0, -1, 3] {
            let mut molecule = molecule(2, &[[1, invalid]], &[SINGLE]);
            assert!(matches!(
                perceive_aromatic_bonds(&mut molecule),
                Err(FingerprintError::AvalonConversion { .. })
            ));
        }
    }
}
