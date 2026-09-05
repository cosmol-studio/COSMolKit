//! Source-shaped Daylight aromaticity used by Avalon fingerprint preprocessing.

use crate::FingerprintError;

use super::preprocess::Neighbourhood;
use super::reaccs::MoleculeState;
use super::rings::{BondSetNode, combine_rings, prepend_fused_ring_pairs, proper_ring_pairs, ring_list};
use super::symbols::{atom_symbol_match, get_symbol_list};

const SINGLE: i32 = 1;
const DOUBLE: i32 = 2;
const TRIPLE: i32 = 3;
const AROMATIC: i32 = 4;

struct DyAromaticCandidates {
    rings: Vec<BondSetNode>,
    bond_is_in_ring: Vec<bool>,
    atom_is_in_ring: Vec<bool>,
    aromaticity_candidate: Vec<bool>,
}

fn build_dy_aromatic_candidates(molecule: &MoleculeState) -> Result<Option<DyAromaticCandidates>, FingerprintError> {
    // Avalon❗✔️:    unsigned i, j;
    // Avalon❗✔️:    int ii;
    // Avalon❗✔️:    int changed;
    // Avalon❗✔️:    bond_set_node *ring_list, *ring_pairs, *plist, *plisth;
    // Avalon❗✔️:    bond_set_node *p, *aromatic_candidates;
    // Avalon❗✔️:    bit_set_t *set;
    // Avalon❗✔️:    pair_t *graph;
    // Avalon❗✔️: //    unsigned graph[MAXBONDS][2];
    // Avalon❗✔️:    int *bond_is_in_ring;
    // Avalon❗✔️: //   int bond_is_in_ring[MAXBONDS];
    // Avalon❗✔️:    int *atom_is_in_ring;
    // Avalon❗✔️: //   int atom_is_in_ring[MAXATOMS];
    // Avalon❗✔️:    int *aromaticity_candidate;
    // Avalon❗✔️: //   int aromaticity_candidate[MAXATOMS];
    // Remaining evaluation-local declarations are reproduced in the
    // evaluation unit that consumes this complete candidate state.
    // Avalon❗✔️:
    // Avalon❗✔️: //   if (mp->n_bonds > MAXBONDS) return;
    // Avalon❗✔️: //   if (mp->n_atoms > MAXATOMS) return;
    // Avalon❗✔️:
    // Avalon❗✔️:    bond_is_in_ring = TypeAlloc(mp->n_bonds, int);
    let mut bond_is_in_ring = vec![false; molecule.bonds.len()];
    // Avalon❗✔️:    graph = TypeAlloc(mp->n_bonds, pair_t);
    let mut graph = Vec::with_capacity(molecule.bonds.len());
    // Avalon❗✔️:    /* Process all rings */
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
    let initial_rings = ring_list(&graph);
    // Avalon❗✔️:    if (IsNULL(ring_list))
    // Avalon❗✔️:    {
    // Avalon❗✔️:       MyFree((char *)bond_is_in_ring);
    // Avalon❗✔️:       MyFree((char *)graph);
    // Avalon❗✔️:       return;
    // Avalon❗✔️:    }
    if initial_rings.is_empty() {
        return Ok(None);
    }
    // Avalon❗✔️:
    // Avalon❗✔️:    atom_is_in_ring = TypeAlloc(mp->n_atoms, int);
    let mut atom_is_in_ring = vec![false; molecule.atoms.len()];
    // Avalon❗✔️:    aromaticity_candidate = TypeAlloc(mp->n_atoms, int);
    let mut aromaticity_candidate = vec![false; molecule.atoms.len()];
    // Avalon❗✔️:    for (plist=ring_list; plist!=(bond_set_node *)NULL; plist=plist->next)
    // Avalon❗✔️:       for (ii=NextMember(plist->bond_set, 0); ii>=0; ii=NextMember(plist->bond_set, (unsigned)ii+1))
    // Avalon❗✔️:          bond_is_in_ring[ii] = TRUE;
    for ring in &initial_rings {
        for (bond_index, in_ring) in bond_is_in_ring.iter_mut().enumerate() {
            if ring.bond_set.contains(bond_index) {
                *in_ring = true;
            }
        }
    }
    // Avalon❗✔️:
    // Avalon❗✔️:    for (i=0; i<mp->n_atoms; i++)
    // Avalon❗✔️:    {
    for (atom_index, atom) in molecule.atoms.iter().enumerate() {
        // Avalon❗✔️:       atom_is_in_ring[i] = FALSE;
        // The zero-initialized Rust vector already holds FALSE.
        // Avalon❗✔️:       if (0 != strcmp("C", mp->atom_array[i].atom_symbol))
        // Avalon❗✔️:          aromaticity_candidate[i] = TRUE;
        // Avalon❗✔️:       else
        // Avalon❗✔️:          aromaticity_candidate[i] = FALSE;
        aromaticity_candidate[atom_index] = atom.atom_symbol != "C";
        // Avalon❗✔️:    }
    }
    // Avalon❗✔️:    for (i=0; i<mp->n_bonds; i++)
    // Avalon❗✔️:       if (mp->bond_array[i].bond_type >  SINGLE  &&
    // Avalon❗✔️:           mp->bond_array[i].bond_type != TRIPLE)
    // Avalon❗✔️:       {
    for bond in &molecule.bonds {
        if bond.bond_type > SINGLE && bond.bond_type != TRIPLE {
            // Avalon❗✔️:          aromaticity_candidate[mp->bond_array[i].atoms[0]-1] = TRUE;
            aromaticity_candidate[bond.atoms[0] as usize - 1] = true;
            // Avalon❗✔️:          aromaticity_candidate[mp->bond_array[i].atoms[1]-1] = TRUE;
            aromaticity_candidate[bond.atoms[1] as usize - 1] = true;
            // Avalon❗✔️:       }
        }
    }
    // Avalon❗✔️:    /* Process just aromaticity candidates */
    // Avalon❗✔️:    DisposeBondSetList(ring_list);
    // Rust releases the initial ring list after its last use.
    // Avalon❗✔️:    for (i=0; i<mp->n_bonds; i++)
    // Avalon❗✔️:    {
    for (bond_index, bond) in molecule.bonds.iter().enumerate() {
        // Avalon❗✔️:       graph[i][0] = mp->bond_array[i].atoms[0];
        graph[bond_index][0] = bond.atoms[0] as usize;
        // Avalon❗✔️:       graph[i][1] = mp->bond_array[i].atoms[1];
        graph[bond_index][1] = bond.atoms[1] as usize;
        // Avalon❗✔️:       if (!bond_is_in_ring[i])  /* only ring bonds can be aromatic */
        // Avalon❗✔️:       {
        if !bond_is_in_ring[bond_index] {
            // Avalon❗✔️:          graph[i][0] = 0;
            graph[bond_index][0] = 0;
            // Avalon❗✔️:          graph[i][1] = 0;
            graph[bond_index][1] = 0;
            // Avalon❗✔️:       }
            // Avalon❗✔️:       else
            // Avalon❗✔️:       {
        } else {
            // Avalon❗✔️:          atom_is_in_ring[graph[i][0]-1] = TRUE;
            atom_is_in_ring[graph[bond_index][0] - 1] = true;
            // Avalon❗✔️:          atom_is_in_ring[graph[i][1]-1] = TRUE;
            atom_is_in_ring[graph[bond_index][1] - 1] = true;
            // Avalon❗✔️:          /* make ring perception ignore bonds with non-candidates */
            // Avalon❗✔️:          if (!aromaticity_candidate[graph[i][0]-1]) graph[i][0] = 0;
            if !aromaticity_candidate[graph[bond_index][0] - 1] {
                graph[bond_index][0] = 0;
            }
            // Avalon❗✔️:          if (!aromaticity_candidate[graph[i][1]-1]) graph[i][1] = 0;
            if !aromaticity_candidate[graph[bond_index][1] - 1] {
                graph[bond_index][1] = 0;
            }
            // Avalon❗✔️:       }
        }
        // Avalon❗✔️:    }
    }
    // Avalon❗✔️:    ring_list = RingList(graph,mp->n_bonds);
    let mut rings = ring_list(&graph);
    // Avalon❗✔️:    if (IsNULL(ring_list))
    // Avalon❗✔️:    {
    // Avalon❗✔️:       MyFree((char *)bond_is_in_ring);
    // Avalon❗✔️:       MyFree((char *)atom_is_in_ring);
    // Avalon❗✔️:       MyFree((char *)aromaticity_candidate);
    // Avalon❗✔️:       MyFree((char *)graph);
    // Avalon❗✔️:       return;
    // Avalon❗✔️:    }
    if rings.is_empty() {
        return Ok(None);
    }
    // Avalon❗✔️:    ring_list = CombineRings(ring_list);
    combine_rings(&mut rings);
    // Avalon❗✔️:
    // Avalon❗✔️:    ring_pairs = ProperRingPairs(ring_list, mp->n_atoms, graph);
    let mut ring_pairs = proper_ring_pairs(&rings, molecule.atoms.len(), &graph);
    // Avalon❗✔️:
    // Avalon❗✔️:    /* add the additional rings for the analysis */
    // Avalon❗✔️:    while (ring_pairs!=(bond_set_node *)NULL)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       plisth = ring_pairs->next;
    // Avalon❗✔️:       ring_pairs->next = ring_list;
    // Avalon❗✔️:       ring_list=ring_pairs;
    // Avalon❗✔️:       ring_pairs = plisth;
    // Avalon❗✔️:    }
    ring_pairs.reverse();
    ring_pairs.append(&mut rings);
    rings = ring_pairs;
    // The complete fused-pair source block is reproduced in the ring helper.
    prepend_fused_ring_pairs(&mut rings);

    Ok(Some(DyAromaticCandidates {
        rings,
        bond_is_in_ring,
        atom_is_in_ring,
        aromaticity_candidate,
    }))
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DyMutation {
    QueryHydrogen(usize),
    RingBond(usize),
    ClosureBond(usize),
}

pub(super) fn perceive_daylight_aromaticity(
    molecule: &mut MoleculeState,
    neighbours: &[Neighbourhood],
) -> Result<(), FingerprintError> {
    perceive_daylight_aromaticity_with(molecule, neighbours, |_| {}).map(|_| ())
}

fn perceive_daylight_aromaticity_with(
    molecule: &mut MoleculeState,
    neighbours: &[Neighbourhood],
    mut record_mutation: impl FnMut(DyMutation),
) -> Result<usize, FingerprintError> {
    // Avalon❗✔️: void PerceiveDYAromaticity(struct reaccs_molecule_t *mp,
    // Avalon❗✔️: 			   neighbourhood_t nbp[])
    // Avalon❗✔️: /*
    // Avalon❗✔️:  * Converts bonds in rings that are perceived as aromatic by Daylight
    // Avalon❗✔️:  * programs into AROMATIC bonds. In addition, hydrogen counts of aromatic
    // Avalon❗✔️:  * non-carbon atoms are remembered in the query_H_count fields.
    // Avalon❗✔️:  */
    // Avalon❗✔️: {
    if neighbours.len() != molecule.atoms.len() {
        return Err(FingerprintError::AvalonConversion {
            reason: "Avalon neighbourhood array has the wrong length".to_string(),
        });
    }
    let Some(mut candidates) = build_dy_aromatic_candidates(molecule)? else {
        return Ok(0);
    };
    // Avalon❗✔️:    int in_ring_double;
    // Avalon❗✔️:    int in_ring_aromatic;
    // Avalon❗✔️:    int exo_pull;
    // Avalon❗✔️:    int conjugated;
    // Avalon❗✔️:    int is_in_ring;
    // Avalon❗✔️:    int local_pi, npi;
    // Avalon❗✔️:
    // Avalon❗✔️:    struct reaccs_atom_t *ap, *alp;
    // Avalon❗✔️:    struct reaccs_bond_t *blp;
    // Avalon❗✔️:    neighbourhood_t *nbph;
    // Rust uses source-order indexed slices for the same records.
    // Avalon❗✔️:
    // Avalon❗✔️:    do
    // Avalon❗✔️:    {
    let mut passes = 0;
    loop {
        passes += 1;
        // Avalon❗✔️:       changed = FALSE;
        let mut changed = false;
        // Avalon❗✔️:       for (plist=aromatic_candidates;
        // Avalon❗✔️:            plist!=(bond_set_node *)NULL;
        // Avalon❗✔️: 	   plist=plist->next)
        // Avalon❗✔️:       {
        for ring in &candidates.rings {
            // Avalon❗✔️: 	 npi = 0;
            let mut npi = 0_i32;
            // Avalon❗✔️: 	 conjugated = TRUE;
            let mut conjugated = true;
            // Avalon❗✔️:          /* Look for atoms that are members of this ring */
            // Avalon❗✔️: 	 for (i=0, ap=mp->atom_array, nbph=nbp;
            // Avalon❗✔️: 	      i<mp->n_atoms;
            // Avalon❗✔️: 	      i++, ap++, nbph++)
            // Avalon❗✔️: 	 {
            for atom_index in 0..molecule.atoms.len() {
                // Avalon❗✔️:             if (!atom_is_in_ring[i]) continue;
                if !candidates.atom_is_in_ring[atom_index] {
                    continue;
                }
                // Avalon❗✔️:             /* Counts # of double bonds in this ring attached to ap */
                // Avalon❗✔️: 	    in_ring_double = 0;
                let mut in_ring_double = 0_i32;
                // Avalon❗✔️:             /* Counts # of aromatic bonds in this ring attached to ap */
                // Avalon❗✔️: 	    in_ring_aromatic = 0;
                let mut in_ring_aromatic = 0_i32;
                // Avalon❗✔️:             /* becomes true if there is an exocyclic DB that pulls electrons */
                // Avalon❗✔️: 	    exo_pull = FALSE;
                let mut exo_pull = false;
                // Avalon❗✔️: 	    is_in_ring = FALSE;
                let mut is_in_ring = false;
                // Avalon❗✔️: 	    local_pi = 0;       /* count pi electrons of this atom */
                let mut local_pi = 0_i32;
                // Avalon❗✔️: 	    for (j=0; j<nbph->n_ligands; j++)
                // Avalon❗✔️: 	    {
                for (&ligand_atom, &ligand_bond) in neighbours[atom_index]
                    .atoms()
                    .iter()
                    .zip(neighbours[atom_index].bonds())
                {
                    // Avalon❗✔️: 	       alp = &mp->atom_array[nbph->atoms[j]];
                    let ligand = &molecule.atoms[ligand_atom];
                    // Avalon❗✔️: 	       blp = &mp->bond_array[nbph->bonds[j]];
                    let bond = &molecule.bonds[ligand_bond];
                    // Avalon❗✔️: 	       if (bond_is_in_ring[nbph->bonds[j]]  &&  IsMember(plist->bond_set, nbph->bonds[j]))
                    // Avalon❗✔️: 	       {
                    if candidates.bond_is_in_ring[ligand_bond] && ring.bond_set.contains(ligand_bond) {
                        // Avalon❗✔️: 		  is_in_ring = TRUE;
                        is_in_ring = true;
                        // Avalon❗✔️: 		  if (blp->bond_type == AROMATIC) in_ring_aromatic++;
                        if bond.bond_type == AROMATIC {
                            in_ring_aromatic += 1;
                        }
                        // Avalon❗✔️: 		  if (blp->bond_type == DOUBLE)   in_ring_double++;
                        if bond.bond_type == DOUBLE {
                            in_ring_double += 1;
                        }
                        // Avalon❗✔️: 	       }
                        // Avalon❗✔️: 	       else if (mp->bond_array[nbph->bonds[j]].bond_type == DOUBLE)
                        // Avalon❗✔️: 	       {
                    } else if bond.bond_type == DOUBLE {
                        // Avalon❗✔️: 		  if (ap->atom_symbol[0] == 'C'  &&  ap->atom_symbol[1] == '\0'  &&
                        // Avalon❗✔️:                       // AtomSymbolMatch(ap->atom_symbol, "C")  &&
                        // Avalon❗✔️:                       !bond_is_in_ring[nbph->bonds[j]]  && /* patches DY bug */
                        // Avalon❗✔️: 		      AtomSymbolMatch(alp->atom_symbol, "O,S,P,N,L"))
                        // Avalon❗✔️: 		     exo_pull = TRUE;
                        if molecule.atoms[atom_index].atom_symbol == "C"
                            && !candidates.bond_is_in_ring[ligand_bond]
                            && atom_symbol_match(&ligand.atom_symbol, "O,S,P,N,L")
                        {
                            exo_pull = true;
                        }
                        // Avalon❗✔️: 	       }
                    }
                    // Avalon❗✔️: 	    }
                }
                // Avalon❗✔️: 	    if (!is_in_ring) continue;  /* no ligands in ring found */
                if !is_in_ring {
                    continue;
                }
                let atom_symbol = molecule.atoms[atom_index].atom_symbol.as_str();
                let atom_charge = molecule.atoms[atom_index].charge;
                let mut set_query_h = false;
                // Avalon❗✔️: 	    if ((in_ring_aromatic >= 1  ||  in_ring_double  == 1)  &&
                // Avalon❗✔️: 	        (AtomSymbolMatch(ap->atom_symbol, "C,N,A,*") ||
                // Avalon❗✔️: 		 0 == strcmp(ap->atom_symbol, "L"))) /* lists may be aromatic */
                // Avalon❗✔️: 	    {
                if (in_ring_aromatic >= 1 || in_ring_double == 1)
                    && (atom_symbol_match(atom_symbol, "C,N,A,*") || atom_symbol == "L")
                {
                    // Avalon❗✔️:                /* Count this <b>atom</b> for one electron */
                    // Avalon❗✔️: 	       local_pi = 1;
                    local_pi = 1;
                    // Avalon❗✔️: 	    }
                    // Avalon❗✔️: 	    else if ((in_ring_aromatic == 0  &&  in_ring_double  == 0)  &&
                    // Avalon❗✔️: 		     ap->charge == 0                                    &&
                    // Avalon❗✔️: 		     (AtomSymbolMatch(ap->atom_symbol, "N,S,O") ||
                    // Avalon❗✔️:                      /* term that catches lists with O,N,S */
                    // Avalon❗✔️:                       AtomSymbolMatch(ap->atom_symbol, "L")  &&
                    // Avalon❗✔️:                       getList(mp,i) != NULL  &&
                    // Avalon❗✔️:                       getList(mp,i)->logic  &&  /* no not-lists */
                    // Avalon❗✔️:                       strchr(getList(mp, i)->string,'C') == NULL))
                    // Avalon❗✔️: 	    {
                } else {
                    let list_supplies_two = atom_symbol_match(atom_symbol, "L")
                        && get_symbol_list(molecule, atom_index)
                            .is_some_and(|list| list.inclusive && !list.symbols.contains('C'));
                    if in_ring_aromatic == 0
                        && in_ring_double == 0
                        && atom_charge == 0
                        && (atom_symbol_match(atom_symbol, "N,S,O") || list_supplies_two)
                    {
                        // Avalon❗✔️: 	       local_pi = 2;
                        local_pi = 2;
                        // Avalon❗✔️: 	       if (nbph->n_ligands  == 2  &&
                        // Avalon❗✔️: 	           in_ring_aromatic == 0  &&
                        // Avalon❗✔️: 	           in_ring_double   == 0  &&
                        // Avalon❗✔️: 	           0 == strcmp(ap->atom_symbol, "N"))
                        // Avalon❗✔️: 		  ap->query_H_count = 1+1;      /* offset 1! */
                        set_query_h = neighbours[atom_index].atoms().len() == 2
                            && in_ring_aromatic == 0
                            && in_ring_double == 0
                            && atom_symbol == "N";
                        // Avalon❗✔️: 	    }
                        // Avalon❗✔️: 	    else if ((in_ring_aromatic == 0  &&  in_ring_double  == 0)  &&
                        // Avalon❗✔️: 		     ap->charge == 0                                    &&
                        // Avalon❗✔️: 		     exo_pull						&&
                        // Avalon❗✔️: 		     0 == strcmp(ap->atom_symbol, "C"))
                        // Avalon❗✔️: 	    {
                    } else if in_ring_aromatic == 0
                        && in_ring_double == 0
                        && atom_charge == 0
                        && exo_pull
                        && atom_symbol == "C"
                    {
                        // Avalon❗✔️: 	       local_pi = 0;
                        local_pi = 0;
                        // Avalon❗✔️: 	    }
                        // Avalon❗✔️: 	    else
                        // Avalon❗✔️:             {
                    } else {
                        // Avalon❗✔️: 	       conjugated = FALSE;
                        conjugated = false;
                        // Avalon❗✔️:             }
                    }
                }
                // Avalon❗✔️:
                // Avalon❗✔️: 	    /* catch special cases that Daylight don't consider aromatic */
                // Avalon❗✔️: 	    if (ap->charge < 0  &&  AtomSymbolMatch(ap->atom_symbol, "C,N"))
                // Avalon❗✔️: 	       conjugated = FALSE;
                if atom_charge < 0 && atom_symbol_match(atom_symbol, "C,N") {
                    conjugated = false;
                }
                // Avalon❗✔️:
                // Avalon❗✔️: 	    npi += local_pi;
                npi += local_pi;
                if set_query_h {
                    molecule.atoms[atom_index].query_h_count = 2;
                    record_mutation(DyMutation::QueryHydrogen(atom_index));
                }
                // Avalon❗✔️: 	 }
            }
            // Avalon❗✔️: 	 if (conjugated == FALSE) continue;	/* system not conjugated */
            if !conjugated {
                continue;
            }
            // Avalon❗✔️: 	 if (npi%4 != 2)          continue;	/* not a Huckel system */
            if npi % 4 != 2 {
                continue;
            }
            // Avalon❗✔️:
            // Avalon❗✔️: 	 /* now we have an aromatic ring => modify bond_types */
            // Avalon❗✔️:          for (i=0; i<mp->n_bonds; i++)
            // Avalon❗✔️:             if (bond_is_in_ring[i]  &&  IsMember(plist->bond_set,i))
            // Avalon❗✔️: 	    {
            for (bond_index, in_ring) in candidates.bond_is_in_ring.iter().copied().enumerate() {
                if in_ring && ring.bond_set.contains(bond_index) {
                    // Avalon❗✔️: 	       if (mp->bond_array[i].bond_type != AROMATIC)
                    // Avalon❗✔️: 	       {
                    if molecule.bonds[bond_index].bond_type != AROMATIC {
                        // Avalon❗✔️: 		  mp->bond_array[i].bond_type = AROMATIC;
                        molecule.bonds[bond_index].bond_type = AROMATIC;
                        record_mutation(DyMutation::RingBond(bond_index));
                        // Avalon❗✔️: 		  changed = TRUE;
                        changed = true;
                        // Avalon❗✔️: 	       }
                    }
                    // Avalon❗✔️: 	    }
                }
            }
            // Avalon❗✔️:       }
        }
        // Avalon❗✔️:    } while (changed);
        if !changed {
            break;
        }
    }
    // Avalon❗✔️:
    // Avalon❗✔️:    DisposeBondSetList(aromatic_candidates);
    candidates.rings.clear();
    // Avalon❗✔️:    // add ring bonds between aromatic atoms as aromatic
    // Avalon❗✔️:    for (i=0; i<mp->n_atoms; i++)
    // Avalon❗✔️:       aromaticity_candidate[i] = FALSE;
    candidates.aromaticity_candidate.fill(false);
    // Avalon❗✔️:    for (i=0, blp=mp->bond_array; i<mp->n_bonds; i++, blp++)
    // Avalon❗✔️:    {
    for bond in &molecule.bonds {
        // Avalon❗✔️:       if (blp->bond_type == AROMATIC)
        // Avalon❗✔️:       {
        if bond.bond_type == AROMATIC {
            // Avalon❗✔️: 	 aromaticity_candidate[blp->atoms[0]-1] = TRUE;
            candidates.aromaticity_candidate[bond.atoms[0] as usize - 1] = true;
            // Avalon❗✔️: 	 aromaticity_candidate[blp->atoms[1]-1] = TRUE;
            candidates.aromaticity_candidate[bond.atoms[1] as usize - 1] = true;
            // Avalon❗✔️:       }
        }
        // Avalon❗✔️:    }
    }
    // Avalon❗✔️:    for (i=0, blp=mp->bond_array; i<mp->n_bonds; i++, blp++)
    // Avalon❗✔️:       if (aromaticity_candidate[blp->atoms[0]-1]  &&
    // Avalon❗✔️:           aromaticity_candidate[blp->atoms[1]-1]  &&
    // Avalon❗✔️: 	  bond_is_in_ring[i]                      &&
    // Avalon❗✔️: 	  blp->bond_type == SINGLE)
    // Avalon❗✔️: 	 blp->bond_type = AROMATIC;
    for (bond_index, bond) in molecule.bonds.iter_mut().enumerate() {
        if candidates.aromaticity_candidate[bond.atoms[0] as usize - 1]
            && candidates.aromaticity_candidate[bond.atoms[1] as usize - 1]
            && candidates.bond_is_in_ring[bond_index]
            && bond.bond_type == SINGLE
        {
            bond.bond_type = AROMATIC;
            record_mutation(DyMutation::ClosureBond(bond_index));
        }
    }
    // Avalon❗✔️:
    // Avalon❗✔️:    MyFree((char *)atom_is_in_ring);
    // Avalon❗✔️:    MyFree((char *)aromaticity_candidate);
    // Avalon❗✔️:    MyFree((char *)bond_is_in_ring);
    // Avalon❗✔️:    MyFree((char *)graph);
    // Rust releases all candidate state on scope exit; graph was released
    // after candidate construction because evaluation never reads it.
    // Avalon❗✔️: }
    Ok(passes)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::properties::avalon_fingerprint::preprocess::setup_neighbourhood;
    use crate::properties::avalon_fingerprint::reaccs::{Atom, Bond};

    fn molecule(symbols: &[&str], endpoints: &[[i32; 2]], types: &[i32]) -> MoleculeState {
        MoleculeState {
            atoms: symbols
                .iter()
                .map(|symbol| Atom {
                    atom_symbol: (*symbol).to_string(),
                    ..Atom::default()
                })
                .collect(),
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

    fn member_lists(rings: &[BondSetNode]) -> Vec<Vec<usize>> {
        rings
            .iter()
            .map(|ring| {
                (0..=ring.bond_set.max_member())
                    .filter(|&bond| ring.bond_set.contains(bond))
                    .collect()
            })
            .collect()
    }

    fn run_daylight(molecule: &mut MoleculeState) -> Vec<DyMutation> {
        let neighbours = setup_neighbourhood(molecule, molecule.atoms.len()).unwrap();
        let mut mutations = Vec::new();
        perceive_daylight_aromaticity_with(molecule, &neighbours, |mutation| mutations.push(mutation)).unwrap();
        mutations
    }

    fn bond_types(molecule: &MoleculeState) -> Vec<i32> {
        molecule.bonds.iter().map(|bond| bond.bond_type).collect()
    }

    #[test]
    fn benzene_candidate_state_matches_native() {
        let endpoints = [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1]];
        let molecule = molecule(&["C"; 6], &endpoints, &[SINGLE, 2, SINGLE, 2, SINGLE, 2]);

        let state = build_dy_aromatic_candidates(&molecule).unwrap().unwrap();

        assert_eq!(member_lists(&state.rings), vec![vec![0, 1, 2, 3, 4, 5]]);
        assert_eq!(state.bond_is_in_ring, vec![true; 6]);
        assert_eq!(state.atom_is_in_ring, vec![true; 6]);
        assert_eq!(state.aromaticity_candidate, vec![true; 6]);
    }

    #[test]
    fn saturated_carbon_ring_has_no_reconstructed_candidate() {
        let endpoints = [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1]];
        let molecule = molecule(&["C"; 6], &endpoints, &[SINGLE; 6]);

        assert!(build_dy_aromatic_candidates(&molecule).unwrap().is_none());
    }

    #[test]
    fn three_fused_candidate_duplicates_and_order_match_native() {
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
        let molecule = molecule(&["N"; 8], &endpoints, &[SINGLE; 10]);

        let state = build_dy_aromatic_candidates(&molecule).unwrap().unwrap();

        assert_eq!(
            member_lists(&state.rings),
            vec![
                vec![0, 1, 3, 4, 5, 6],
                vec![2, 4, 6, 7, 8, 9],
                vec![0, 1, 3, 4, 6, 7, 8, 9],
                vec![0, 1, 3, 4, 6, 7, 8, 9],
                vec![2, 4, 6, 7, 8, 9],
                vec![0, 1, 3, 4, 5, 6],
                vec![5, 7, 8, 9],
                vec![2, 4, 5, 6],
                vec![0, 1, 2, 3],
            ]
        );
        assert_eq!(state.bond_is_in_ring, vec![true; 10]);
        assert_eq!(state.atom_is_in_ring, vec![true; 8]);
        assert_eq!(state.aromaticity_candidate, vec![true; 8]);
    }

    #[test]
    fn candidate_builder_rejects_invalid_one_based_endpoints() {
        for invalid in [0, -1, 3] {
            let molecule = molecule(&["C"; 2], &[[1, invalid]], &[SINGLE]);
            assert!(matches!(
                build_dy_aromatic_candidates(&molecule),
                Err(FingerprintError::AvalonConversion { .. })
            ));
        }
    }

    #[test]
    fn benzene_and_pyridine_match_native_bond_and_mutation_order() {
        let endpoints = [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1]];
        for symbols in [["C", "C", "C", "C", "C", "C"], ["N", "C", "C", "C", "C", "C"]] {
            let mut molecule = molecule(&symbols, &endpoints, &[DOUBLE, SINGLE, DOUBLE, SINGLE, DOUBLE, SINGLE]);
            let mutations = run_daylight(&mut molecule);

            assert_eq!(bond_types(&molecule), vec![AROMATIC; 6]);
            assert_eq!(mutations, (0..6).map(DyMutation::RingBond).collect::<Vec<_>>());
            assert!(molecule.atoms.iter().all(|atom| atom.query_h_count == 0));
        }
    }

    #[test]
    fn pyrrole_furan_and_symbol_lists_match_native_local_pi_rules() {
        let endpoints = [[1, 2], [2, 3], [3, 4], [4, 5], [5, 1]];
        let types = [SINGLE, DOUBLE, SINGLE, DOUBLE, SINGLE];

        let mut pyrrole = molecule(&["N", "C", "C", "C", "C"], &endpoints, &types);
        let mutations = run_daylight(&mut pyrrole);
        assert_eq!(bond_types(&pyrrole), vec![AROMATIC; 5]);
        assert_eq!(pyrrole.atoms[0].query_h_count, 2);
        assert_eq!(mutations[0], DyMutation::QueryHydrogen(0));
        assert_eq!(&mutations[1..], &(0..5).map(DyMutation::RingBond).collect::<Vec<_>>());

        let mut furan = molecule(&["O", "C", "C", "C", "C"], &endpoints, &types);
        assert_eq!(
            run_daylight(&mut furan),
            (0..5).map(DyMutation::RingBond).collect::<Vec<_>>()
        );
        assert_eq!(bond_types(&furan), vec![AROMATIC; 5]);

        let mut list = molecule(&["L", "C", "C", "C", "C"], &endpoints, &types);
        list.symbol_lists.push(super::super::reaccs::SymbolList {
            atom: 1,
            inclusive: true,
            symbols: "N,O,S".to_string(),
        });
        run_daylight(&mut list);
        assert_eq!(bond_types(&list), vec![AROMATIC; 5]);

        let mut carbon_list = molecule(&["L", "C", "C", "C", "C"], &endpoints, &types);
        carbon_list.symbol_lists.push(super::super::reaccs::SymbolList {
            atom: 1,
            inclusive: true,
            symbols: "N,Cl".to_string(),
        });
        assert!(run_daylight(&mut carbon_list).is_empty());
        assert_eq!(bond_types(&carbon_list), types);
    }

    #[test]
    fn negative_carbon_or_nitrogen_disables_daylight_aromaticity() {
        let endpoints = [[1, 2], [2, 3], [3, 4], [4, 5], [5, 1]];
        let types = [SINGLE, DOUBLE, SINGLE, DOUBLE, SINGLE];
        let mut molecule = molecule(&["N", "C", "C", "C", "C"], &endpoints, &types);
        molecule.atoms[0].charge = -1;

        assert!(run_daylight(&mut molecule).is_empty());
        assert_eq!(bond_types(&molecule), types);
        assert_eq!(molecule.atoms[0].query_h_count, 0);
    }

    #[test]
    fn final_closure_marks_shared_single_bond_between_aromatic_atoms() {
        let endpoints = [[1, 2], [2, 3], [3, 4], [4, 1], [3, 5], [5, 6], [6, 4]];
        let mut molecule = molecule(
            &["C"; 6],
            &endpoints,
            &[DOUBLE, SINGLE, SINGLE, SINGLE, DOUBLE, SINGLE, DOUBLE],
        );

        let mutations = run_daylight(&mut molecule);

        assert_eq!(bond_types(&molecule), vec![AROMATIC; 7]);
        assert_eq!(mutations.last(), Some(&DyMutation::ClosureBond(2)));
        assert_eq!(
            &mutations[..mutations.len() - 1],
            &[0, 1, 3, 4, 5, 6]
                .into_iter()
                .map(DyMutation::RingBond)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn exocyclic_double_bond_pull_contributes_zero_pi_electrons() {
        let endpoints = [[1, 2], [2, 3], [3, 1], [1, 4]];
        let mut molecule = molecule(&["C", "C", "C", "O"], &endpoints, &[SINGLE, DOUBLE, SINGLE, DOUBLE]);

        let mutations = run_daylight(&mut molecule);

        assert_eq!(bond_types(&molecule), vec![AROMATIC, AROMATIC, AROMATIC, DOUBLE]);
        assert_eq!(
            mutations,
            vec![
                DyMutation::RingBond(0),
                DyMutation::RingBond(1),
                DyMutation::RingBond(2),
            ]
        );
    }

    #[test]
    fn query_h_is_written_before_a_later_nonconjugated_atom_rejects_the_ring() {
        let endpoints = [[1, 2], [2, 3], [3, 4], [4, 5], [5, 1]];
        let mut molecule = molecule(
            &["N", "Xe", "C", "C", "Xe"],
            &endpoints,
            &[SINGLE, SINGLE, DOUBLE, SINGLE, SINGLE],
        );

        let mutations = run_daylight(&mut molecule);

        assert_eq!(mutations, vec![DyMutation::QueryHydrogen(0)]);
        assert_eq!(bond_types(&molecule), vec![SINGLE, SINGLE, DOUBLE, SINGLE, SINGLE]);
        assert_eq!(molecule.atoms[0].query_h_count, 2);
    }

    #[test]
    fn daylight_requires_matching_neighbourhood_length() {
        let mut molecule = molecule(&["C"; 2], &[[1, 2]], &[SINGLE]);
        assert!(matches!(
            perceive_daylight_aromaticity(&mut molecule, &[]),
            Err(FingerprintError::AvalonConversion { reason })
                if reason == "Avalon neighbourhood array has the wrong length"
        ));
    }

    #[test]
    fn fused_rings_require_two_mutating_passes_and_one_stability_pass() {
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
            &["N", "N", "C", "C", "C", "C", "C", "C"],
            &endpoints,
            &[
                SINGLE, SINGLE, DOUBLE, SINGLE, SINGLE, DOUBLE, SINGLE, SINGLE, SINGLE, SINGLE,
            ],
        );
        let neighbours = setup_neighbourhood(&molecule, molecule.atoms.len()).unwrap();
        let mut mutations = Vec::new();

        let passes =
            perceive_daylight_aromaticity_with(&mut molecule, &neighbours, |mutation| mutations.push(mutation))
                .unwrap();

        assert_eq!(passes, 3);
        assert_eq!(
            mutations,
            vec![
                DyMutation::QueryHydrogen(0),
                DyMutation::QueryHydrogen(1),
                DyMutation::QueryHydrogen(0),
                DyMutation::QueryHydrogen(1),
                DyMutation::QueryHydrogen(0),
                DyMutation::QueryHydrogen(1),
                DyMutation::RingBond(0),
                DyMutation::RingBond(1),
                DyMutation::RingBond(2),
                DyMutation::RingBond(3),
                DyMutation::RingBond(4),
                DyMutation::RingBond(5),
                DyMutation::RingBond(6),
            ]
        );
        assert_eq!(
            bond_types(&molecule),
            vec![
                AROMATIC, AROMATIC, AROMATIC, AROMATIC, AROMATIC, AROMATIC, AROMATIC, SINGLE, SINGLE, SINGLE
            ]
        );
        assert_eq!(
            molecule.atoms.iter().map(|atom| atom.query_h_count).collect::<Vec<_>>(),
            vec![2, 2, 0, 0, 0, 0, 0, 0]
        );
    }
}
