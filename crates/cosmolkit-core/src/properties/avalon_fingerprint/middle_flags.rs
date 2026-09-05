//! Source-order implementation of Avalon fingerprint flags `0x000020..=0x000400`.

use crate::FingerprintError;

use super::AvalonFingerprintFlags;
use super::fingerprint_state::{FingerprintPreprocessingState, with_prepared_fingerprint_state};
use super::hash::next_hash;
use super::low_flags::{
    prepare_atom_class_state, prepare_bond_path_state, prepare_ring_path_state,
};
use super::reaccs::MoleculeState;
use super::traversal::{
    AUGMENTED_ATOM_SEED, AUGMENTED_BOND_SEED, BOND_PATH_SEED, HCOUNT_CLASS_PATH_SEED,
    HCOUNT_PAIR_SEED, HCOUNT_PATH_SEED, PathFlags, add_bit, set_path_bits_rec,
};

const SUB_AS_IS: i32 = -2;
const SUB_ONE: i32 = 1;
const SUB_MORE: i32 = 6;
const SPECIAL_RING: u32 = 0xfc & !(1 << 6);
const ACCUMULATE_BITS: i32 = 0x0002;

pub(super) fn count_middle_flag_families(
    molecule: &mut MoleculeState,
    counts: &mut [i32],
    which_bits: AvalonFingerprintFlags,
    as_query: bool,
    fpflags: i32,
    exclude_atom: i32,
) -> Result<i32, FingerprintError> {
    if counts.is_empty() {
        return Err(FingerprintError::InvalidArguments {
            reason: "Avalon fingerprint count array must not be empty",
        });
    }
    // Avalon❗✔️:    if (0 == (fpflags & ACCUMULATE_BITS))
    // Avalon❗✔️:       for (i=0; i<ncounts; i++) fp_counts[i] = 0;
    if fpflags & ACCUMULATE_BITS == 0 {
        counts.fill(0);
    }
    with_prepared_fingerprint_state(
        molecule,
        which_bits,
        as_query,
        fpflags,
        exclude_atom,
        |working, state| {
            Ok(count_middle_flag_families_prepared(
                working,
                state,
                counts,
                which_bits,
                as_query,
                exclude_atom,
            ))
        },
    )
}

pub(super) fn count_middle_flag_families_prepared(
    working: &mut MoleculeState,
    state: &FingerprintPreprocessingState,
    counts: &mut [i32],
    which_bits: AvalonFingerprintFlags,
    as_query: bool,
    exclude_atom: i32,
) -> i32 {
    let mut result = 0_i32;
    if which_bits.contains(AvalonFingerprintFlags::AUGMENTED_ATOM) {
        result += count_augmented_atoms(working, state, counts, as_query, exclude_atom);
    }
    if which_bits.contains(AvalonFingerprintFlags::AUGMENTED_BOND) {
        result += count_augmented_bonds(working, state, counts, exclude_atom);
    }
    if which_bits.contains(AvalonFingerprintFlags::HCOUNT_PAIR) {
        result += count_hydrogen_pairs(working, state, counts, exclude_atom);
    }
    if which_bits.contains(AvalonFingerprintFlags::HCOUNT_PATH) {
        result += count_hydrogen_paths(working, state, counts, exclude_atom);
    }
    prepare_ring_path_state(working, state, exclude_atom);
    prepare_bond_path_state(working, exclude_atom);
    if which_bits.contains(AvalonFingerprintFlags::BOND_PATH) {
        result += count_bond_paths(working, state, counts, exclude_atom);
    }
    prepare_atom_class_state(working, exclude_atom);
    if which_bits.contains(AvalonFingerprintFlags::HCOUNT_CLASS_PATH) {
        result += count_hydrogen_class_paths(working, state, counts, exclude_atom);
    }
    result
}

fn count_augmented_atoms(
    molecule: &MoleculeState,
    state: &FingerprintPreprocessingState,
    counts: &mut [i32],
    as_query: bool,
    exclude_atom: i32,
) -> i32 {
    // Avalon❗✔️: // fprintf(stderr, "USE_AUGMENTED_ATOM\n");
    // Avalon❗✔️:    if (which_bits & USE_AUGMENTED_ATOM)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       strcpy(prefix, "AA:");
    // Avalon❗✔️:       /* Set bits for all triples of atoms connected to a common atom */
    let mut result = 0_i32;
    // Avalon❗✔️:       for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    // Avalon❗✔️:       {
    for (atom_index, atom) in molecule.atoms.iter().enumerate() {
        // Avalon❗✔️:          if (ap->color <= 0) continue;
        // Avalon❗✔️: 	 if (i+1 == exclude_atom) continue;
        if atom.color <= 0 || atom_index as i32 + 1 == exclude_atom {
            continue;
        }

        // Avalon❗✔️: 	 /* Set bit for atoms with more than one double or a triple bond */
        // Avalon❗✔️: 	 if (nspecial[i] >= 2)
        // Avalon❗✔️: 	 {
        if state.special_neighbours[atom_index] >= 2 {
            // Avalon❗✔️: 	    ADD_BIT(fp_counts, ncounts, NEXT_SEED(ap->color*AUGMENTED_ATOM_SEED,101));
            add_bit(
                counts,
                next_hash((atom.color as u64).wrapping_mul(AUGMENTED_ATOM_SEED), 101),
            );
            // Avalon❗✔️: 	    // ADD_BIT(fp_counts, ncounts, NEXT_SEED(ap->color*AUGMENTED_ATOM_SEED,301));
            // Avalon❗✔️: 	    result += 1;
            result += 1;
            // Avalon❗✔️: 	 }
        }
        // Avalon❗✔️:          old_seed = seed;
        // Avalon❗✔️:          // Add some bits for hydrogen counted or hetero central atoms with
        // Avalon❗✔️:          // hetero neighbours
        // Avalon❗✔️:          if ((H_count[i+1] > 0 || ap->color != 6)  && degree[i] >= 2)
        if (state.hydrogen_counts[atom_index + 1] > 0 || atom.color != 6)
            && state.degree[atom_index] >= 2
        {
            // Avalon❗✔️:             for (i1=0; i1<nbp[i].n_ligands; i1++)
            // Avalon❗✔️: 	    {
            let neighbourhood = &state.neighbours[atom_index];
            for first_position in 0..neighbourhood.atoms().len() {
                let first_atom = neighbourhood.atoms()[first_position];
                let first_bond = neighbourhood.bonds()[first_position];
                // Avalon❗✔️: 	       if (mp->atom_array[nbp[i].atoms[i1]].color == 0) continue;
                // Avalon❗✔️: 	       if (nbp[i].atoms[i1]+1 == exclude_atom) continue;
                // Avalon❗✔️: 	       if (mp->bond_array[nbp[i].bonds[i1]].color == 0) continue;
                if molecule.atoms[first_atom].color == 0
                    || first_atom as i32 + 1 == exclude_atom
                    || molecule.bonds[first_bond].color == 0
                {
                    continue;
                }
                // Avalon❗✔️:                for (i2=i1+1; i2<nbp[i].n_ligands; i2++)
                // Avalon❗✔️:                {
                for second_position in first_position + 1..neighbourhood.atoms().len() {
                    let second_atom = neighbourhood.atoms()[second_position];
                    let second_bond = neighbourhood.bonds()[second_position];
                    // Avalon❗✔️:                   seed = NEXT_SEED(AUGMENTED_ATOM_SEED, 97);
                    let mut seed = next_hash(AUGMENTED_ATOM_SEED, 97);
                    // Avalon❗✔️:                   seed = NEXT_SEED(seed, ap->color);
                    seed = next_hash(seed, atom.color as u64);
                    // Avalon❗✔️:                   if (mp->atom_array[nbp[i].atoms[i2]].color == 0) continue;
                    // Avalon❗✔️: 		  if (nbp[i].atoms[i2]+1 == exclude_atom) continue;
                    // Avalon❗✔️:                   if (mp->bond_array[nbp[i].bonds[i2]].color == 0) continue;
                    if molecule.atoms[second_atom].color == 0
                        || second_atom as i32 + 1 == exclude_atom
                        || molecule.bonds[second_bond].color == 0
                    {
                        continue;
                    }
                    // Avalon❗✔️:                   if (mp->atom_array[nbp[i].atoms[i1]].color == 6  &&
                    // Avalon❗✔️:                       mp->atom_array[nbp[i].atoms[i2]].color == 6) continue;
                    if molecule.atoms[first_atom].color == 6
                        && molecule.atoms[second_atom].color == 6
                    {
                        continue;
                    }
                    // Avalon❗✔️:                   sum = 0;
                    // Avalon❗✔️:                   sum += mp->atom_array[nbp[i].atoms[i1]].color *
                    // Avalon❗✔️:                          mp->bond_array[nbp[i].bonds[i1]].color;
                    // Avalon❗✔️:                   sum += mp->atom_array[nbp[i].atoms[i2]].color *
                    // Avalon❗✔️:                          mp->bond_array[nbp[i].bonds[i2]].color;
                    let sum = molecule.atoms[first_atom].color * molecule.bonds[first_bond].color
                        + molecule.atoms[second_atom].color * molecule.bonds[second_bond].color;
                    // Avalon❗✔️:                   prod = 1;
                    // Avalon❗✔️:                   prod *= mp->atom_array[nbp[i].atoms[i1]].color; prod &= 0xFFF;
                    // Avalon❗✔️:                   prod *= mp->atom_array[nbp[i].atoms[i2]].color; prod &= 0xFFF;
                    let mut product = molecule.atoms[first_atom].color & 0x0fff;
                    product = product.wrapping_mul(molecule.atoms[second_atom].color) & 0x0fff;
                    // Avalon❗✔️:                   seed = NEXT_SEED(seed, sum);
                    seed = next_hash(seed, sum as u64);
                    // Avalon❗✔️:                   seed = NEXT_SEED(seed, prod);
                    seed = next_hash(seed, product as u64);
                    // Avalon❗✔️:                   ADD_BIT(fp_counts, ncounts, seed);
                    add_bit(counts, seed);
                    // Avalon❗✔️:                   // seed = NEXT_SEED(seed, (sum*prod)&0xFFF);
                    // Avalon❗✔️:                   // ADD_BIT(fp_counts, ncounts, seed);
                    // Avalon❗✔️: 		  result += 1;
                    result += 1;
                    // Avalon❗✔️: 		  nmulti = 0;
                    let mut multiple = 0_i32;
                    // Avalon❗✔️:                   if (mp->bond_array[nbp[i].bonds[i1]].color >= 2) nmulti++;
                    if molecule.bonds[first_bond].color >= 2 {
                        multiple += 1;
                    }
                    // Avalon❗✔️:                   if (mp->bond_array[nbp[i].bonds[i2]].color >= 2) nmulti++;
                    if molecule.bonds[second_bond].color >= 2 {
                        multiple += 1;
                    }
                    // Avalon❗✔️:                   // add bit for hetero-substituted unsaturated hetero atom if degree>1 is well-defined => catch e.g. nitroso vs. nitro
                    // Avalon❗✔️:                   if (ap->color != 6  &&  nmulti > 0  &&
                    // Avalon❗✔️:                       (!as_query                 ||
                    // Avalon❗✔️:                        ap->sub_desc == SUB_AS_IS ||
                    // Avalon❗✔️:                        (ap->sub_desc != NONE     &&
                    // Avalon❗✔️:                         ap->sub_desc != SUB_MORE &&
                    // Avalon❗✔️:                         ap->sub_desc == degree[i]+SUB_ONE-1)))
                    if atom.color != 6
                        && multiple > 0
                        && (!as_query
                            || atom.sub_desc == SUB_AS_IS
                            || (atom.sub_desc != 0
                                && atom.sub_desc != SUB_MORE
                                && atom.sub_desc == state.degree[atom_index] + SUB_ONE - 1))
                    {
                        // Avalon❗✔️:                   {
                        // Avalon❗✔️:                      ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, nmulti+173*degree[i]));
                        add_bit(
                            counts,
                            next_hash(seed, (multiple + 173 * state.degree[atom_index]) as u64),
                        );
                        // Avalon❗✔️:                      result += 1;
                        result += 1;
                        // Avalon❗✔️:                   }
                    }
                    // Avalon❗✔️:                }
                }
                // Avalon❗✔️:             }
            }
        }

        // Avalon❗✔️:          if (degree[i] <= 2) continue;
        if state.degree[atom_index] <= 2 {
            continue;
        }
        let neighbourhood = &state.neighbours[atom_index];
        // Avalon❗✔️:          for (i1=0; i1<nbp[i].n_ligands; i1++)
        // Avalon❗✔️:          {
        for first_position in 0..neighbourhood.atoms().len() {
            let first_atom = neighbourhood.atoms()[first_position];
            let first_bond = neighbourhood.bonds()[first_position];
            // Avalon❗✔️:             if (nbp[i].atoms[i1]+1 == exclude_atom) continue;
            if first_atom as i32 + 1 == exclude_atom {
                continue;
            }
            // Avalon❗✔️:             for (i2=i1+1; i2<nbp[i].n_ligands; i2++)
            // Avalon❗✔️: 	    {
            for second_position in first_position + 1..neighbourhood.atoms().len() {
                let second_atom = neighbourhood.atoms()[second_position];
                let second_bond = neighbourhood.bonds()[second_position];
                // Avalon❗✔️: 	       if (nbp[i].atoms[i2]+1 == exclude_atom) continue;
                if second_atom as i32 + 1 == exclude_atom {
                    continue;
                }
                // Avalon❗✔️:                for (i3=i2+1; i3<nbp[i].n_ligands; i3++)
                // Avalon❗✔️:                {
                for third_position in second_position + 1..neighbourhood.atoms().len() {
                    let third_atom = neighbourhood.atoms()[third_position];
                    let third_bond = neighbourhood.bonds()[third_position];
                    // Avalon❗✔️: 		  if (nbp[i].atoms[i3]+1 == exclude_atom) continue;
                    if third_atom as i32 + 1 == exclude_atom {
                        continue;
                    }
                    // Avalon❗✔️:                   seed = AUGMENTED_ATOM_SEED;
                    // Avalon❗✔️:                   seed = NEXT_SEED(seed, ap->color);
                    let mut seed = next_hash(AUGMENTED_ATOM_SEED, atom.color as u64);
                    // Avalon❗✔️:                   if (mp->atom_array[nbp[i].atoms[i1]].color == 0) continue;
                    // Avalon❗✔️:                   if (mp->atom_array[nbp[i].atoms[i2]].color == 0) continue;
                    // Avalon❗✔️:                   if (mp->atom_array[nbp[i].atoms[i3]].color == 0) continue;
                    // Avalon❗✔️:                   if (mp->bond_array[nbp[i].bonds[i1]].color == 0) continue;
                    // Avalon❗✔️:                   if (mp->bond_array[nbp[i].bonds[i2]].color == 0) continue;
                    // Avalon❗✔️:                   if (mp->bond_array[nbp[i].bonds[i3]].color == 0) continue;
                    if molecule.atoms[first_atom].color == 0
                        || molecule.atoms[second_atom].color == 0
                        || molecule.atoms[third_atom].color == 0
                        || molecule.bonds[first_bond].color == 0
                        || molecule.bonds[second_bond].color == 0
                        || molecule.bonds[third_bond].color == 0
                    {
                        continue;
                    }
                    // Avalon❗✔️: 		  /* count hetero neighbours */
                    // Avalon❗✔️: 		  nqtmp = 0;
                    let mut hetero = 0_i32;
                    // Avalon❗✔️:                   if (mp->atom_array[nbp[i].atoms[i1]].color != 6) nqtmp++;
                    // Avalon❗✔️:                   if (mp->atom_array[nbp[i].atoms[i2]].color != 6) nqtmp++;
                    // Avalon❗✔️:                   if (mp->atom_array[nbp[i].atoms[i3]].color != 6) nqtmp++;
                    for neighbour_atom in [first_atom, second_atom, third_atom] {
                        if molecule.atoms[neighbour_atom].color != 6 {
                            hetero += 1;
                        }
                    }

                    // Avalon❗✔️: 		  /* make sure to add some bits for really odd ones */
                    // Avalon❗✔️: 		  nmulti = 0;
                    let mut multiple = 0_i32;
                    // Avalon❗✔️:                   if (mp->bond_array[nbp[i].bonds[i1]].color == 2) nmulti++;
                    // Avalon❗✔️:                   if (mp->bond_array[nbp[i].bonds[i1]].color == 3) nmulti++;
                    // Avalon❗✔️:                   if (mp->bond_array[nbp[i].bonds[i2]].color == 2) nmulti++;
                    // Avalon❗✔️:                   if (mp->bond_array[nbp[i].bonds[i2]].color == 3) nmulti++;
                    // Avalon❗✔️:                   if (mp->bond_array[nbp[i].bonds[i3]].color == 2) nmulti++;
                    // Avalon❗✔️:                   if (mp->bond_array[nbp[i].bonds[i3]].color == 3) nmulti++;
                    for neighbour_bond in [first_bond, second_bond, third_bond] {
                        if matches!(molecule.bonds[neighbour_bond].color, 2 | 3) {
                            multiple += 1;
                        }
                    }

                    // Avalon❗✔️:                   sum = 0;
                    // Avalon❗✔️:                   sum += mp->atom_array[nbp[i].atoms[i1]].color *
                    // Avalon❗✔️:                          mp->bond_array[nbp[i].bonds[i1]].color;
                    // Avalon❗✔️:                   sum += mp->atom_array[nbp[i].atoms[i2]].color *
                    // Avalon❗✔️:                          mp->bond_array[nbp[i].bonds[i2]].color;
                    // Avalon❗✔️:                   sum += mp->atom_array[nbp[i].atoms[i3]].color *
                    // Avalon❗✔️:                          mp->bond_array[nbp[i].bonds[i3]].color;
                    let sum = molecule.atoms[first_atom].color * molecule.bonds[first_bond].color
                        + molecule.atoms[second_atom].color * molecule.bonds[second_bond].color
                        + molecule.atoms[third_atom].color * molecule.bonds[third_bond].color;
                    // Avalon❗✔️:                   prod = 1;
                    // Avalon❗✔️:                   prod *= mp->atom_array[nbp[i].atoms[i1]].color; prod &= 0xFFF;
                    // Avalon❗✔️:                   prod *= mp->atom_array[nbp[i].atoms[i2]].color; prod &= 0xFFF;
                    // Avalon❗✔️:                   prod *= mp->atom_array[nbp[i].atoms[i3]].color; prod &= 0xFFF;
                    let mut product = molecule.atoms[first_atom].color & 0x0fff;
                    product = product.wrapping_mul(molecule.atoms[second_atom].color) & 0x0fff;
                    product = product.wrapping_mul(molecule.atoms[third_atom].color) & 0x0fff;
                    // Avalon❗✔️:                   seed = NEXT_SEED(seed, sum);
                    seed = next_hash(seed, sum as u64);
                    // Avalon❗✔️:                   seed = NEXT_SEED(seed, prod);
                    seed = next_hash(seed, product as u64);
                    // Avalon❗✔️:                   if ((nqtmp > 2  ||  nmulti >= 2)  &&  ap->color == 6)
                    if (hetero > 2 || multiple >= 2) && atom.color == 6 {
                        // Avalon❗✔️:                   {
                        // Avalon❗✔️:                      ADD_BIT(fp_counts, ncounts, seed);
                        add_bit(counts, seed);
                        // Avalon❗✔️:                      ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, 73));
                        add_bit(counts, next_hash(seed, 73));
                        // Avalon❗✔️: 		     result += 2;
                        result += 2;
                        // Avalon❗✔️:                   }
                    }
                    // Avalon❗✔️:                   seed = NEXT_SEED(seed, (sum*prod)&0xFFF);
                    seed = next_hash(seed, sum.wrapping_mul(product) as u64 & 0x0fff);
                    // Avalon❗✔️:                   ADD_BIT(fp_counts, ncounts, seed);
                    add_bit(counts, seed);
                    // Avalon❗✔️: 		  result++;
                    result += 1;
                    // Avalon❗✔️: 		  if (nmulti >= 2  ||  ap->color > 6)   // make sure R-NO2 is covered
                    if multiple >= 2 || atom.color > 6 {
                        // Avalon❗✔️:                   {
                        // Avalon❗✔️:                       ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, 53)); result++;
                        add_bit(counts, next_hash(seed, 53));
                        result += 1;
                        // Avalon❗✔️:                   }
                    }
                    // Avalon❗✔️:                }
                }
                // Avalon❗✔️: 	    }
            }
            // Avalon❗✔️: 	 }
        }
        // Avalon❗✔️:       }
    }
    // Avalon❗✔️:    }
    result
}

fn count_augmented_bonds(
    molecule: &MoleculeState,
    state: &FingerprintPreprocessingState,
    counts: &mut [i32],
    exclude_atom: i32,
) -> i32 {
    // Avalon❗✔️: // fprintf(stderr, "USE_AUGMENTED_BOND\n");
    // Avalon❗✔️:    if (which_bits & USE_AUGMENTED_BOND)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       strcpy(prefix, "AB:");
    // Avalon❗✔️:       /* Set bits for all bonds with both end-degrees > 2 */
    let mut result = 0_i32;
    // Avalon❗✔️:       for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
    // Avalon❗✔️:       {
    for (bond_index, bond) in molecule.bonds.iter().enumerate() {
        // Avalon❗✔️:          if (bp->atoms[0] == exclude_atom) continue;
        // Avalon❗✔️:          if (bp->atoms[1] == exclude_atom) continue;
        if bond.atoms.contains(&exclude_atom) {
            continue;
        }
        // Avalon❗✔️:          ai1 = bp->atoms[0]-1;
        // Avalon❗✔️:          ai2 = bp->atoms[1]-1;
        let first_endpoint = bond.atoms[0] as usize - 1;
        let second_endpoint = bond.atoms[1] as usize - 1;
        // Avalon❗✔️:          if (degree[ai1] <= 2) continue;
        // Avalon❗✔️:          if (degree[ai2] <= 2) continue;
        if state.degree[first_endpoint] <= 2 || state.degree[second_endpoint] <= 2 {
            continue;
        }
        let first_neighbours = &state.neighbours[first_endpoint];
        let second_neighbours = &state.neighbours[second_endpoint];
        // Avalon❗✔️:          for (i1=0; i1<nbp[ai1].n_ligands; i1++)
        for first_position in 0..first_neighbours.atoms().len() {
            // Avalon❗✔️:             for (i2=i1+1; i2<nbp[ai1].n_ligands; i2++)
            // Avalon❗✔️:             {
            for second_position in first_position + 1..first_neighbours.atoms().len() {
                let first_atom = first_neighbours.atoms()[first_position];
                let second_atom = first_neighbours.atoms()[second_position];
                let first_bond = first_neighbours.bonds()[first_position];
                let second_bond = first_neighbours.bonds()[second_position];
                // Avalon❗✔️:                /* don't reuse current bond */
                // Avalon❗✔️:                if (nbp[ai1].bonds[i1] == i) continue;
                // Avalon❗✔️:                if (nbp[ai1].bonds[i2] == i) continue;
                if first_bond == bond_index || second_bond == bond_index {
                    continue;
                }
                // Avalon❗✔️:                if (nbp[ai1].atoms[i1]+1 == exclude_atom) continue;
                // Avalon❗✔️:                if (nbp[ai1].atoms[i2]+1 == exclude_atom) continue;
                if first_atom as i32 + 1 == exclude_atom || second_atom as i32 + 1 == exclude_atom {
                    continue;
                }
                // Avalon❗✔️:                /* don't use masked atoms/bonds */
                // Avalon❗✔️:                if (mp->atom_array[nbp[ai1].atoms[i1]].color==0) continue;
                // Avalon❗✔️:                if (mp->atom_array[nbp[ai1].atoms[i2]].color==0) continue;
                // Avalon❗✔️:                if (mp->bond_array[nbp[ai1].bonds[i1]].color==0) continue;
                // Avalon❗✔️:                if (mp->bond_array[nbp[ai1].bonds[i2]].color==0) continue;
                if molecule.atoms[first_atom].color == 0
                    || molecule.atoms[second_atom].color == 0
                    || molecule.bonds[first_bond].color == 0
                    || molecule.bonds[second_bond].color == 0
                {
                    continue;
                }
                // Avalon❗✔️:                sumi = 0;
                // Avalon❗✔️:                sumi += mp->atom_array[nbp[ai1].atoms[i1]].color*
                // Avalon❗✔️:                        mp->bond_array[nbp[ai1].bonds[i1]].color;
                // Avalon❗✔️:                sumi += mp->atom_array[nbp[ai1].atoms[i2]].color*
                // Avalon❗✔️:                        mp->bond_array[nbp[ai1].bonds[i2]].color;
                let mut first_sum = molecule.atoms[first_atom].color
                    * molecule.bonds[first_bond].color
                    + molecule.atoms[second_atom].color * molecule.bonds[second_bond].color;
                // Avalon❗✔️:                prodi = 1;
                // Avalon❗✔️:                prodi *= mp->atom_array[nbp[ai1].atoms[i1]].color+
                // Avalon❗✔️:                         mp->bond_array[nbp[ai1].bonds[i1]].color;
                // Avalon❗✔️:                prodi *= mp->atom_array[nbp[ai1].atoms[i2]].color+
                // Avalon❗✔️:                         mp->bond_array[nbp[ai1].bonds[i2]].color;
                let mut first_product = (molecule.atoms[first_atom].color
                    + molecule.bonds[first_bond].color)
                    .wrapping_mul(
                        molecule.atoms[second_atom].color + molecule.bonds[second_bond].color,
                    );
                // Avalon❗✔️:                sumi  &= 0x0FFF;
                // Avalon❗✔️:                prodi &= 0x0FFF;
                first_sum &= 0x0fff;
                first_product &= 0x0fff;
                // Avalon❗✔️:                for (j1=0; j1<nbp[ai2].n_ligands; j1++)
                for third_position in 0..second_neighbours.atoms().len() {
                    // Avalon❗✔️:                   for (j2=j1+1; j2<nbp[ai2].n_ligands; j2++)
                    // Avalon❗✔️:                   {
                    for fourth_position in third_position + 1..second_neighbours.atoms().len() {
                        let third_atom = second_neighbours.atoms()[third_position];
                        let fourth_atom = second_neighbours.atoms()[fourth_position];
                        let third_bond = second_neighbours.bonds()[third_position];
                        let fourth_bond = second_neighbours.bonds()[fourth_position];
                        // Avalon❗✔️:                      /* don't reuse current bond */
                        // Avalon❗✔️:                      if (nbp[ai2].bonds[j1] == i) continue;
                        // Avalon❗✔️:                      if (nbp[ai2].bonds[j2] == i) continue;
                        if third_bond == bond_index || fourth_bond == bond_index {
                            continue;
                        }
                        // Avalon❗✔️:                      if (nbp[ai2].atoms[j1]+1 == exclude_atom) continue;
                        // Avalon❗✔️:                      if (nbp[ai2].atoms[j2]+1 == exclude_atom) continue;
                        if third_atom as i32 + 1 == exclude_atom
                            || fourth_atom as i32 + 1 == exclude_atom
                        {
                            continue;
                        }
                        // Avalon❗✔️:                      /* don't use masked atoms/bonds */
                        // Avalon❗✔️:                      if (mp->atom_array[nbp[ai2].atoms[j1]].color==0) continue;
                        // Avalon❗✔️:                      if (mp->atom_array[nbp[ai2].atoms[j2]].color==0) continue;
                        // Avalon❗✔️:                      if (mp->bond_array[nbp[ai2].bonds[j1]].color==0) continue;
                        // Avalon❗✔️:                      if (mp->bond_array[nbp[ai2].bonds[j2]].color==0) continue;
                        if molecule.atoms[third_atom].color == 0
                            || molecule.atoms[fourth_atom].color == 0
                            || molecule.bonds[third_bond].color == 0
                            || molecule.bonds[fourth_bond].color == 0
                        {
                            continue;
                        }
                        // Avalon❗✔️:                      sumj = 0;
                        // Avalon❗✔️:                      sumj += mp->atom_array[nbp[ai2].atoms[j1]].color*
                        // Avalon❗✔️:                              mp->bond_array[nbp[ai2].bonds[j1]].color;
                        // Avalon❗✔️:                      sumj += mp->atom_array[nbp[ai2].atoms[j2]].color*
                        // Avalon❗✔️:                              mp->bond_array[nbp[ai2].bonds[j2]].color;
                        let mut second_sum = molecule.atoms[third_atom].color
                            * molecule.bonds[third_bond].color
                            + molecule.atoms[fourth_atom].color * molecule.bonds[fourth_bond].color;
                        // Avalon❗✔️:                      prodj = 1;
                        // Avalon❗✔️:                      prodj *= mp->atom_array[nbp[ai2].atoms[j1]].color+
                        // Avalon❗✔️:                               mp->bond_array[nbp[ai2].bonds[j1]].color;
                        // Avalon❗✔️:                      prodj *= mp->atom_array[nbp[ai2].atoms[j2]].color+
                        // Avalon❗✔️:                               mp->bond_array[nbp[ai2].bonds[j2]].color;
                        let mut second_product = (molecule.atoms[third_atom].color
                            + molecule.bonds[third_bond].color)
                            .wrapping_mul(
                                molecule.atoms[fourth_atom].color
                                    + molecule.bonds[fourth_bond].color,
                            );
                        // Avalon❗✔️:                      sumj  &= 0x0FFF;
                        // Avalon❗✔️:                      prodj &= 0x0FFF;
                        second_sum &= 0x0fff;
                        second_product &= 0x0fff;
                        // Avalon❗✔️:                      seed = AUGMENTED_BOND_SEED;
                        // Avalon❗✔️:                      seed = NEXT_SEED(seed, bp->color);
                        let mut seed = next_hash(AUGMENTED_BOND_SEED, bond.color as u64);
                        // Avalon❗✔️:                      seed = NEXT_SEED(seed, prodi+prodj);
                        seed = next_hash(seed, first_product.wrapping_add(second_product) as u64);
                        // Avalon❗✔️:                      seed = NEXT_SEED(seed, sumi*sumj);
                        seed = next_hash(seed, first_sum.wrapping_mul(second_sum) as u64);
                        // Avalon❗✔️:                      ADD_BIT(fp_counts, ncounts, seed);
                        add_bit(counts, seed);
                        // Avalon❗✔️: 		     result++;
                        result += 1;
                        // Avalon❗✔️:                   }
                    }
                }
                // Avalon❗✔️:             }
            }
        }
        // Avalon❗✔️:       }
    }
    // Avalon❗✔️:    }
    result
}

fn count_hydrogen_pairs(
    molecule: &MoleculeState,
    state: &FingerprintPreprocessingState,
    counts: &mut [i32],
    exclude_atom: i32,
) -> i32 {
    // Avalon❗✔️: // fprintf(stderr, "USE_HCOUNT_PAIR\n");
    // Avalon❗✔️:    if (1*which_bits & USE_HCOUNT_PAIR)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       strcpy(prefix, "HPR:");
    // Avalon❗✔️:       /* generate bits for hydrogen counted described bonds */
    let mut result = 0_i32;
    // Avalon❗✔️:       for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
    // Avalon❗✔️:       {
    for bond in &molecule.bonds {
        let first = bond.atoms[0] as usize - 1;
        let second = bond.atoms[1] as usize - 1;
        let first_hydrogens = state.hydrogen_counts[bond.atoms[0] as usize];
        let second_hydrogens = state.hydrogen_counts[bond.atoms[1] as usize];
        // Avalon❗✔️:          if (H_count[bp->atoms[0]] == 0  &&  H_count[bp->atoms[1]] == 0)
        // Avalon❗✔️:             continue;
        if first_hydrogens == 0 && second_hydrogens == 0 {
            continue;
        }
        // Avalon❗✔️:          if (bp->atoms[0] == exclude_atom) continue;
        // Avalon❗✔️:          if (bp->atoms[1] == exclude_atom) continue;
        if bond.atoms.contains(&exclude_atom) {
            continue;
        }
        // Avalon❗✔️:          /* Don't consider CC single bonds */
        // Avalon❗✔️:          if (mp->atom_array[bp->atoms[0]-1].color == 6  &&
        // Avalon❗✔️:              mp->atom_array[bp->atoms[1]-1].color == 6  &&
        // Avalon❗✔️:              bp->color == 1) continue;
        if molecule.atoms[first].color == 6 && molecule.atoms[second].color == 6 && bond.color == 1
        {
            continue;
        }
        // Avalon❗✔️:          /* Don't consider explicit AH bonds */
        // Avalon❗✔️:          if (mp->atom_array[bp->atoms[0]-1].color == 0) continue;
        // Avalon❗✔️:          if (mp->atom_array[bp->atoms[1]-1].color == 0) continue;
        if molecule.atoms[first].color == 0 || molecule.atoms[second].color == 0 {
            continue;
        }
        // Avalon❗✔️:          for (j1=0; j1<=H_count[bp->atoms[0]]; j1++)
        for first_count in 0..=first_hydrogens {
            // Avalon❗✔️:             for (j2=0; j2<=H_count[bp->atoms[1]]; j2++)
            // Avalon❗✔️:             {
            for second_count in 0..=second_hydrogens {
                // Avalon❗✔️:                if (j1+j2 == 0) continue;        /* at least one hydrogen */
                if first_count + second_count == 0 {
                    continue;
                }
                // Avalon❗✔️:                // unsaturation triggers bit like a hydrogen
                // Avalon❗✔️:                if (!unsaturated[bp->atoms[0]-1]  &&  !unsaturated[bp->atoms[1]-1]  &&
                // Avalon❗✔️:                    j1*j2 == 0  && bp->color <= 1) continue;
                if !state.unsaturated[first]
                    && !state.unsaturated[second]
                    && first_count * second_count == 0
                    && bond.color <= 1
                {
                    continue;
                }
                // Avalon❗✔️:                if (j1+j2 > 3) continue;         /* at most 3 hydrogens */
                if first_count + second_count > 3 {
                    continue;
                }
                // Avalon❗✔️:                // seed = HCOUNT_PAIR_SEED + 53*(j1+j2) + 7*(j1+1)*(j2+1);
                // Avalon❗✔️:                seed = NEXT_SEED(HCOUNT_PAIR_SEED, 7*(j1+1)*(j2+1));
                let mut seed = next_hash(
                    HCOUNT_PAIR_SEED,
                    (7 * (first_count + 1) * (second_count + 1)) as u64,
                );
                // Avalon❗✔️:                seed = NEXT_SEED(seed, 53*(j1+j2));
                seed = next_hash(seed, (53 * (first_count + second_count)) as u64);
                // Avalon❗✔️:                seed = NEXT_SEED(seed, bp->color);
                seed = next_hash(seed, bond.color as u64);
                // Avalon❗✔️:                seed = NEXT_SEED(seed, mp->atom_array[bp->atoms[0]-1].color +
                // Avalon❗✔️:                                       mp->atom_array[bp->atoms[1]-1].color);
                seed = next_hash(
                    seed,
                    (molecule.atoms[first].color + molecule.atoms[second].color) as u64,
                );
                // Avalon❗✔️:                seed = NEXT_SEED(seed, mp->atom_array[bp->atoms[0]-1].color *
                // Avalon❗✔️:                                       mp->atom_array[bp->atoms[1]-1].color);
                seed = next_hash(
                    seed,
                    molecule.atoms[first]
                        .color
                        .wrapping_mul(molecule.atoms[second].color) as u64,
                );
                // Avalon❗✔️:                ADD_BIT(fp_counts, ncounts, seed);
                add_bit(counts, seed);
                // Avalon❗✔️: 	       result++;
                result += 1;
                // Avalon❗✔️:                if (bp->color > 1)
                // Avalon❗✔️:                {
                if bond.color > 1 {
                    // Avalon❗✔️:                   seed = NEXT_SEED(seed, 83);
                    seed = next_hash(seed, 83);
                    // Avalon❗✔️:                   ADD_BIT(fp_counts, ncounts, seed);
                    add_bit(counts, seed);
                    // Avalon❗✔️: 		  result++;
                    result += 1;
                    // Avalon❗✔️:                   /* Add more bits for really special pairs */
                    // Avalon❗✔️:                   if (j1+j2 == 1  &&
                    // Avalon❗✔️:                       (mp->atom_array[bp->atoms[0]-1].color != 6  ||
                    // Avalon❗✔️:                        mp->atom_array[bp->atoms[1]-1].color != 6) &&
                    // Avalon❗✔️:                       bp->color <= 3)
                    if first_count + second_count == 1
                        && (molecule.atoms[first].color != 6 || molecule.atoms[second].color != 6)
                        && bond.color <= 3
                    {
                        // Avalon❗✔️:                   {
                        // Avalon❗✔️:                      seed = NEXT_SEED(seed, 91);
                        seed = next_hash(seed, 91);
                        // Avalon❗✔️:                      ADD_BIT(fp_counts, ncounts, seed);
                        add_bit(counts, seed);
                        // Avalon❗✔️:                      seed = NEXT_SEED(seed, 97);
                        seed = next_hash(seed, 97);
                        // Avalon❗✔️:                      ADD_BIT(fp_counts, ncounts, seed);
                        add_bit(counts, seed);
                        // Avalon❗✔️:                      seed = NEXT_SEED(seed, 103);
                        seed = next_hash(seed, 103);
                        // Avalon❗✔️:                      ADD_BIT(fp_counts, ncounts, seed);
                        add_bit(counts, seed);
                        // Avalon❗✔️: 		     result += 3;
                        result += 3;
                        // Avalon❗✔️:                   }
                    }
                    // Avalon❗✔️:                }
                }
                // Avalon❗✔️:             }
            }
        }
        // Avalon❗✔️:       }
    }
    // Avalon❗✔️:    }
    result
}

fn count_hydrogen_paths(
    molecule: &MoleculeState,
    state: &FingerprintPreprocessingState,
    counts: &mut [i32],
    exclude_atom: i32,
) -> i32 {
    // Avalon❗✔️: // fprintf(stderr, "USE_HCOUNT_PATH\n");
    // Avalon❗✔️:    if (which_bits & USE_HCOUNT_PATH)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       strcpy(prefix, "HPH:");
    // Avalon❗✔️:       /* generate a short path for each atom that has a hydrogen */
    // Avalon❗✔️:       seed = HCOUNT_PATH_SEED;
    let parent_seed = HCOUNT_PATH_SEED;
    let mut result = 0_i32;
    let mut touched = vec![0_i32; molecule.atoms.len()];
    // Avalon❗✔️:       for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    // Avalon❗✔️:       {
    for (atom_index, atom) in molecule.atoms.iter().enumerate() {
        // Avalon❗✔️:          if (ap->color <= 0) continue;
        // Avalon❗✔️:          if (i+1 == exclude_atom) continue;
        // Avalon❗✔️:          if (H_count[i+1] == 0) continue;
        if atom.color <= 0
            || atom_index as i32 + 1 == exclude_atom
            || state.hydrogen_counts[atom_index + 1] == 0
        {
            continue;
        }
        // Avalon❗✔️:          /* don't consider non methyl carbon atoms */
        // Avalon❗✔️:          if (ap->color == 6  &&  H_count[i+1] < 3  && degree[i] < 3) continue;
        if atom.color == 6
            && state.hydrogen_counts[atom_index + 1] < 3
            && state.degree[atom_index] < 3
        {
            continue;
        }

        // Avalon❗✔️:          touched_indices[i] = 1; /* updating */
        touched[atom_index] = 1;
        // Avalon❗✔️:          old_seed = seed;
        // Avalon❗✔️:          seed = NEXT_SEED(seed, ap->color);
        let seed = next_hash(parent_seed, atom.color as u64);
        // Avalon❗✔️:          if (ap->color != 6)
        if atom.color != 6 {
            // Avalon❗✔️:             result += SetPathBitsRec(mp, nbp,
            // Avalon❗✔️:                                      fp_counts, ncounts,
            // Avalon❗✔️:                                      seed, touched_indices,
            // Avalon❗✔️:                                      1,
            // Avalon❗✔️:                                      2, 5,  /* path length 2 to 4 */ // EVG
            // Avalon❗✔️:                                      i, 0, -1,
            // Avalon❗✔️:                                      IGNORE_PATH_SYMBOL |
            // Avalon❗✔️: // IGNORE_TERM_SYMBOL |
            // Avalon❗✔️:                                      PROCESS_CHAINS,
            // Avalon❗✔️:                                      exclude_atom,
            // Avalon❗✔️:                                      prefix);
            result += set_path_bits_rec(
                molecule,
                &state.neighbours,
                counts,
                seed,
                &mut touched,
                1,
                2,
                5,
                atom_index,
                0,
                -1,
                PathFlags::IGNORE_PATH_SYMBOL | PathFlags::PROCESS_CHAINS,
                exclude_atom,
            );
        // Avalon❗✔️:          else
        // Avalon❗✔️:          {
        } else {
            // Avalon❗✔️:             if (degree[i] > 2)  // tertiary hydrogen
            // Avalon❗✔️:             {
            if state.degree[atom_index] > 2 {
                // Avalon❗✔️:                if (0)          // Class disabled to save bit density
                // Avalon❗✔️:                result += SetPathBitsRec(mp, nbp,
                // Avalon❗✔️:                                         fp_counts, ncounts,
                // Avalon❗✔️:                                         NEXT_SEED(seed, 101), touched_indices,
                // Avalon❗✔️:                                         1,
                // Avalon❗✔️:                                         2, 2,  /* path length 2 to 3 */
                // Avalon❗✔️:                                         i, 0, -1,
                // Avalon❗✔️:                                         IGNORE_PATH_SYMBOL |
                // Avalon❗✔️:                                         IGNORE_TERM_SYMBOL |
                // Avalon❗✔️:                                         PROCESS_CHAINS,
                // Avalon❗✔️:                                         exclude_atom,
                // Avalon❗✔️:                                         prefix);
                // Disabled in the source by the literal `if (0)`.
                // Avalon❗✔️:                result += SetPathBitsRec(mp, nbp,
                // Avalon❗✔️:                                         fp_counts, ncounts,
                // Avalon❗✔️:                                         NEXT_SEED(seed, 103), touched_indices,
                // Avalon❗✔️:                                         1,
                // Avalon❗✔️:                                         2, 5,  /* path length 2 to 3 */
                // Avalon❗✔️:                                         i, 0, -1,
                // Avalon❗✔️:                                         IGNORE_PATH_SYMBOL |
                // Avalon❗✔️:                                         FORCED_HETERO_END |  /* i-Pr...Q */
                // Avalon❗✔️:                                         PROCESS_CHAINS,
                // Avalon❗✔️:                                         exclude_atom,
                // Avalon❗✔️:                                         prefix);
                result += set_path_bits_rec(
                    molecule,
                    &state.neighbours,
                    counts,
                    next_hash(seed, 103),
                    &mut touched,
                    1,
                    2,
                    5,
                    atom_index,
                    0,
                    -1,
                    PathFlags::IGNORE_PATH_SYMBOL
                        | PathFlags::FORCED_HETERO_END
                        | PathFlags::PROCESS_CHAINS,
                    exclude_atom,
                );
                // Avalon❗✔️:             }
            }
            // Avalon❗✔️:             if (H_count[i+1] >= 3)      // methyl
            // Avalon❗✔️:             {
            if state.hydrogen_counts[atom_index + 1] >= 3 {
                // Avalon❗✔️:                if (0)          // Class disabled to save bit density
                // Avalon❗✔️:                result += SetPathBitsRec(mp, nbp,
                // Avalon❗✔️:                                         fp_counts, ncounts,
                // Avalon❗✔️:                                         NEXT_SEED(seed,1103), touched_indices,
                // Avalon❗✔️:                                         1,
                // Avalon❗✔️:                                         2, 3,  /* path length 2 to 3 */
                // Avalon❗✔️:                                         i, 0, -1,
                // Avalon❗✔️:                                         IGNORE_PATH_SYMBOL |
                // Avalon❗✔️:                                         PROCESS_CHAINS,
                // Avalon❗✔️:                                         exclude_atom,
                // Avalon❗✔️:                                         prefix);
                // Disabled in the source by the literal `if (0)`.
                // Avalon❗✔️:                result += SetPathBitsRec(mp, nbp,
                // Avalon❗✔️:                                         fp_counts, ncounts,
                // Avalon❗✔️:                                         NEXT_SEED(seed,1105), touched_indices,
                // Avalon❗✔️:                                         1,
                // Avalon❗✔️:                                         3, 6,  /* path length 4 to 4 */
                // Avalon❗✔️:                                         i, 0, -1,
                // Avalon❗✔️:                                         IGNORE_PATH_SYMBOL |
                // Avalon❗✔️:                                         FORCED_HETERO_END |  /* Me...Q */
                // Avalon❗✔️:                                         PROCESS_CHAINS,
                // Avalon❗✔️:                                         exclude_atom,
                // Avalon❗✔️:                                         prefix);
                result += set_path_bits_rec(
                    molecule,
                    &state.neighbours,
                    counts,
                    next_hash(seed, 1105),
                    &mut touched,
                    1,
                    3,
                    6,
                    atom_index,
                    0,
                    -1,
                    PathFlags::IGNORE_PATH_SYMBOL
                        | PathFlags::FORCED_HETERO_END
                        | PathFlags::PROCESS_CHAINS,
                    exclude_atom,
                );
                // Avalon❗✔️:             }
            }
            // Avalon❗✔️:          }
        }
        // Avalon❗✔️:          touched_indices[i] = 0; /* down-dating */
        touched[atom_index] = 0;
        // Avalon❗✔️:          seed = old_seed;
        // Avalon❗✔️:          if (ap->color == 6) continue;
        if atom.color == 6 {
            continue;
        }

        // Avalon❗✔️:          if (H_count[i+1] > 1)  /* catch the difference between NH and NH2! */
        // Avalon❗✔️:          {
        if state.hydrogen_counts[atom_index + 1] > 1 {
            // Avalon❗✔️:             touched_indices[i] = 1; /* updating */
            touched[atom_index] = 1;
            // Avalon❗✔️:             old_seed = seed;
            // Avalon❗✔️:             seed = NEXT_SEED(seed, 113);
            // Avalon❗✔️:             seed = NEXT_SEED(seed, ap->color);
            let seed = next_hash(next_hash(parent_seed, 113), atom.color as u64);
            // Avalon❗✔️:             result += SetPathBitsRec(mp, nbp,
            // Avalon❗✔️:                                      fp_counts, ncounts,
            // Avalon❗✔️:                                      seed, touched_indices,
            // Avalon❗✔️:                                      1,
            // Avalon❗✔️:                                      1, 5,  /* path length 1 to 5 */
            // Avalon❗✔️:                                      i, 0, -1,
            // Avalon❗✔️: // FORCED_HETERO_END |
            // Avalon❗✔️:                                      IGNORE_PATH_SYMBOL |
            // Avalon❗✔️:                                      PROCESS_CHAINS,
            // Avalon❗✔️:                                      exclude_atom,
            // Avalon❗✔️:                                      prefix);
            result += set_path_bits_rec(
                molecule,
                &state.neighbours,
                counts,
                seed,
                &mut touched,
                1,
                1,
                5,
                atom_index,
                0,
                -1,
                PathFlags::IGNORE_PATH_SYMBOL | PathFlags::PROCESS_CHAINS,
                exclude_atom,
            );
            // Avalon❗✔️:             touched_indices[i] = 0; /* down-dating */
            touched[atom_index] = 0;
            // Avalon❗✔️:             seed = old_seed;
            // Avalon❗✔️:          }
        }

        let mut seed = parent_seed;
        // Avalon❗✔️:          for (j=1; j<H_count[i+1]; j++)
        // Avalon❗✔️:          {
        for hydrogen_index in 1..state.hydrogen_counts[atom_index + 1] {
            // Avalon❗✔️:             seed = NEXT_SEED(seed, 61*j);
            seed = next_hash(seed, (61 * hydrogen_index) as u64);
            // Avalon❗✔️:             seed = NEXT_SEED(seed, ap->color);
            seed = next_hash(seed, atom.color as u64);
            // Avalon❗✔️:             ADD_BIT(fp_counts, ncounts, seed);
            add_bit(counts, seed);
            // Avalon❗✔️: 	    result++;
            result += 1;
            // Avalon❗✔️:          }
        }
        // Avalon❗✔️:          seed = old_seed;
        // Avalon❗✔️:       }
    }
    // Avalon❗✔️:    }
    result
}

fn count_bond_paths(
    molecule: &MoleculeState,
    state: &FingerprintPreprocessingState,
    counts: &mut [i32],
    exclude_atom: i32,
) -> i32 {
    // Avalon❗✔️: // fprintf(stderr, "USE_BOND_PATH\n");
    // Avalon❗✔️:    // Here, we have all non-trivial atoms mapped to ANY_COLOR while the
    // Avalon❗✔️:    // bond type is retained.
    // Avalon❗✔️:    if (which_bits & USE_BOND_PATH)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       strcpy(prefix, "BP:");
    // Avalon❗✔️:       seed = BOND_PATH_SEED;
    let parent_seed = BOND_PATH_SEED;
    let mut result = 0_i32;
    let mut touched = vec![0_i32; molecule.atoms.len()];
    // Avalon❗✔️:       for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    // Avalon❗✔️:       {
    for (atom_index, atom) in molecule.atoms.iter().enumerate() {
        // Avalon❗✔️:          if (ap->color <= 0) continue;
        // Avalon❗✔️:          if (i+1 == exclude_atom) continue;
        if atom.color <= 0 || atom_index as i32 + 1 == exclude_atom {
            continue;
        }
        // Avalon❗✔️:          touched_indices[i] = 1; /* updating */
        touched[atom_index] = 1;
        // Avalon❗✔️:          old_seed = seed;
        // Avalon❗✔️:          // start at branch node on ring
        // Avalon❗✔️:          if (degree[i] > 2  &&  atom_status[i] > 1)
        // Avalon❗✔️:          {
        if state.degree[atom_index] > 2 && state.atom_status[atom_index] > 1 {
            // Avalon❗✔️:             result += SetPathBitsRec(mp, nbp,
            // Avalon❗✔️:                                      fp_counts, ncounts,
            // Avalon❗✔️:                                      seed, touched_indices,
            // Avalon❗✔️:                                      1,
            // Avalon❗✔️:                                      4, 4,         /* was 4 to 4 */
            // Avalon❗✔️:                                      i, 0, -1,
            // Avalon❗✔️:                                      PROCESS_CHAINS,
            // Avalon❗✔️:                                      exclude_atom,
            // Avalon❗✔️:                                      prefix);
            result += set_path_bits_rec(
                molecule,
                &state.neighbours,
                counts,
                parent_seed,
                &mut touched,
                1,
                4,
                4,
                atom_index,
                0,
                -1,
                PathFlags::PROCESS_CHAINS,
                exclude_atom,
            );
            // Avalon❗✔️:             if (ap->rsize_flags&SPECIAL_RING)
            if atom.rsize_flags & SPECIAL_RING != 0 {
                // Avalon❗✔️: 	           result += SetPathBitsRec(mp, nbp,
                // Avalon❗✔️:                                         fp_counts, ncounts,
                // Avalon❗✔️:                                         NEXT_SEED(seed,217), touched_indices,
                // Avalon❗✔️:                                         1,
                // Avalon❗✔️:                                         5, 5,
                // Avalon❗✔️:                                         i, 0, -1,
                // Avalon❗✔️:                                         IGNORE_PATH_SYMBOL |
                // Avalon❗✔️:                                         PROCESS_CHAINS,
                // Avalon❗✔️:                                         exclude_atom,
                // Avalon❗✔️:                                         prefix);
                result += set_path_bits_rec(
                    molecule,
                    &state.neighbours,
                    counts,
                    next_hash(parent_seed, 217),
                    &mut touched,
                    1,
                    5,
                    5,
                    atom_index,
                    0,
                    -1,
                    PathFlags::IGNORE_PATH_SYMBOL | PathFlags::PROCESS_CHAINS,
                    exclude_atom,
                );
            }
            // Avalon❗✔️:          }
        }
        // Avalon❗✔️:          /* Add other bits to catch poorly specified ring closures */
        // Avalon❗✔️:          seed = old_seed;
        // Avalon❗✔️:          seed = NEXT_SEED(seed, 11);
        let ring_seed = next_hash(parent_seed, 11);
        // Avalon❗✔️:          result += SetPathBitsRec(mp, nbp,
        // Avalon❗✔️:                                   fp_counts, ncounts,
        // Avalon❗✔️:                                   seed, touched_indices,
        // Avalon❗✔️:                                   1,
        // Avalon❗✔️:                                   4, 6,         /* was 5 to 6 */
        // Avalon❗✔️:                                   i, 0, -1,
        // Avalon❗✔️:                                   // DEBUG_PATH |
        // Avalon❗✔️:                                   FORCED_RING_PATH |
        // Avalon❗✔️:                                   IGNORE_PATH_SYMBOL |
        // Avalon❗✔️:                                   PROCESS_RING_CLOSURES,
        // Avalon❗✔️:                                   exclude_atom,
        // Avalon❗✔️:                                   prefix);
        result += set_path_bits_rec(
            molecule,
            &state.neighbours,
            counts,
            ring_seed,
            &mut touched,
            1,
            4,
            6,
            atom_index,
            0,
            -1,
            PathFlags::FORCED_RING_PATH
                | PathFlags::IGNORE_PATH_SYMBOL
                | PathFlags::PROCESS_RING_CLOSURES,
            exclude_atom,
        );
        // Avalon❗✔️:          /* Add bits for paths starting with rare bond orders */
        // Avalon❗✔️:          for (j=0; j<nbp[i].n_ligands; j++)
        // Avalon❗✔️:          {
        for (&neighbour_atom, &neighbour_bond) in state.neighbours[atom_index]
            .atoms()
            .iter()
            .zip(state.neighbours[atom_index].bonds())
        {
            // Avalon❗✔️:             bp = &mp->bond_array[nbp[i].bonds[j]];
            // Avalon❗✔️:             ai = nbp[i].atoms[j];
            let bond = &molecule.bonds[neighbour_bond];
            // Avalon❗✔️:             if (ai+1 == exclude_atom) continue;
            // Avalon❗✔️:             if (bp->color == 0) continue;
            if neighbour_atom as i32 + 1 == exclude_atom || bond.color == 0 {
                continue;
            }
            // Avalon❗✔️:             if (bp->bond_type != DOUBLE  &&
            // Avalon❗✔️:                 bp->bond_type != TRIPLE)
            // Avalon❗✔️:                continue;
            if !matches!(bond.bond_type, 2 | 3) {
                continue;
            }
            // Avalon❗✔️:             if (atom_status[i] <= 0  &&  bp->bond_type != TRIPLE)
            // Avalon❗✔️:                continue;
            if state.atom_status[atom_index] <= 0 && bond.bond_type != 3 {
                continue;
            }
            // Avalon❗✔️:             seed = old_seed;
            // Avalon❗✔️:             seed = NEXT_SEED(seed, bp->color*413);
            let seed = next_hash(parent_seed, (bond.color * 413) as u64);
            // Avalon❗✔️:             touched_indices[ai] = 1;    /* updating */
            touched[neighbour_atom] = 1;
            // Avalon❗✔️:             result += SetPathBitsRec(mp, nbp,
            // Avalon❗✔️:                                      fp_counts, ncounts,
            // Avalon❗✔️:                                      seed, touched_indices,
            // Avalon❗✔️:                                      2,
            // Avalon❗✔️:                                      4, 5,
            // Avalon❗✔️:                                      ai, 0, i,
            // Avalon❗✔️:                                      IGNORE_PATH_SYMBOL |
            // Avalon❗✔️:                                      PROCESS_RING_CLOSURES |
            // Avalon❗✔️:                                      PROCESS_CHAINS,
            // Avalon❗✔️:                                      exclude_atom,
            // Avalon❗✔️:                                      prefix);
            result += set_path_bits_rec(
                molecule,
                &state.neighbours,
                counts,
                seed,
                &mut touched,
                2,
                4,
                5,
                neighbour_atom,
                0,
                atom_index as i32,
                PathFlags::IGNORE_PATH_SYMBOL
                    | PathFlags::PROCESS_RING_CLOSURES
                    | PathFlags::PROCESS_CHAINS,
                exclude_atom,
            );
            // Avalon❗✔️:             touched_indices[ai] = 0;    /* down-dating */
            touched[neighbour_atom] = 0;
            // Avalon❗✔️:          }
        }
        // Avalon❗✔️:          seed = old_seed;
        // Avalon❗✔️:          touched_indices[i] = 0; /* down-dating */
        touched[atom_index] = 0;
        // Avalon❗✔️:       }
    }
    // Avalon❗✔️:    }
    result
}

fn count_hydrogen_class_paths(
    molecule: &MoleculeState,
    state: &FingerprintPreprocessingState,
    counts: &mut [i32],
    exclude_atom: i32,
) -> i32 {
    // Avalon❗✔️: // fprintf(stderr, "USE_HCOUNT_CLASS_PATH\n");
    // Avalon❗✔️:    // Here, we unify atom types to carbon/hetero distinction
    // Avalon❗✔️:    if (which_bits & USE_HCOUNT_CLASS_PATH)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       strcpy(prefix, "HCP:");
    // Avalon❗✔️:       /* generate a short path for each atom that has a hydrogen */
    // Avalon❗✔️:       seed = HCOUNT_CLASS_PATH_SEED;
    let parent_seed = HCOUNT_CLASS_PATH_SEED;
    let mut result = 0_i32;
    let mut touched = vec![0_i32; molecule.atoms.len()];
    // Avalon❗✔️:       for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    // Avalon❗✔️:       {
    for (atom_index, atom) in molecule.atoms.iter().enumerate() {
        // Avalon❗✔️:          if (ap->color <= 0) continue;
        // Avalon❗✔️:          if (i+1 == exclude_atom) continue;
        // Avalon❗✔️:          if (H_count[i+1] == 0) continue;
        // Avalon❗✔️:          if (ap->color == 6  &&  H_count[i+1] < 2) continue;
        if atom.color <= 0
            || atom_index as i32 + 1 == exclude_atom
            || state.hydrogen_counts[atom_index + 1] == 0
            || (atom.color == 6 && state.hydrogen_counts[atom_index + 1] < 2)
        {
            continue;
        }
        // Avalon❗✔️:          old_seed = seed;
        // Avalon❗✔️:          seed = NEXT_SEED(seed, ap->color);
        let seed = next_hash(parent_seed, atom.color as u64);
        // Avalon❗✔️:          touched_indices[i] = 1; /* updating */
        touched[atom_index] = 1;
        // Avalon❗✔️:          if (ap->color == 6)
        if atom.color == 6 {
            // Avalon❗✔️:             result += SetPathBitsRec(mp, nbp,
            // Avalon❗✔️:                                      fp_counts, ncounts,
            // Avalon❗✔️:                                      seed, touched_indices,
            // Avalon❗✔️:                                      1,
            // Avalon❗✔️:                                      2, 4,  /* path length 1 to 3 */
            // Avalon❗✔️:                                      i, 0, -1,
            // Avalon❗✔️:                                      FORCED_HETERO_END |
            // Avalon❗✔️:                                      PROCESS_CHAINS,
            // Avalon❗✔️: 				     exclude_atom,
            // Avalon❗✔️: 				     prefix);
            result += set_path_bits_rec(
                molecule,
                &state.neighbours,
                counts,
                seed,
                &mut touched,
                1,
                2,
                4,
                atom_index,
                0,
                -1,
                PathFlags::FORCED_HETERO_END | PathFlags::PROCESS_CHAINS,
                exclude_atom,
            );
        // Avalon❗✔️:          else
        } else {
            // Avalon❗✔️:             result += SetPathBitsRec(mp, nbp,
            // Avalon❗✔️:                                      fp_counts, ncounts,
            // Avalon❗✔️:                                      seed, touched_indices,
            // Avalon❗✔️:                                      1,
            // Avalon❗✔️:                                      2, 5,
            // Avalon❗✔️:                                      i, 0, -1,
            // Avalon❗✔️:                                      IGNORE_PATH_SYMBOL |
            // Avalon❗✔️:                                      FORCED_HETERO_END |
            // Avalon❗✔️:                                      PROCESS_CHAINS,
            // Avalon❗✔️: 				     exclude_atom,
            // Avalon❗✔️: 				     prefix);
            result += set_path_bits_rec(
                molecule,
                &state.neighbours,
                counts,
                seed,
                &mut touched,
                1,
                2,
                5,
                atom_index,
                0,
                -1,
                PathFlags::IGNORE_PATH_SYMBOL
                    | PathFlags::FORCED_HETERO_END
                    | PathFlags::PROCESS_CHAINS,
                exclude_atom,
            );
        }
        // Avalon❗✔️:          seed = old_seed;
        // Avalon❗✔️:          touched_indices[i] = 0; /* down-dating */
        touched[atom_index] = 0;
        // Avalon❗✔️:       }
    }
    // Avalon❗✔️:    }
    result
}

#[cfg(test)]
mod tests {
    use super::super::fingerprint_state::USE_DY_AROMATICITY;
    use super::super::reaccs::{Atom, Bond};
    use super::*;

    const MIDDLE_MASKS: &[u32] = &[
        0x020, 0x040, 0x080, 0x100, 0x200, 0x400, 0x060, 0x0a0, 0x120, 0x220, 0x420, 0x0c0, 0x140,
        0x240, 0x440, 0x180, 0x280, 0x480, 0x300, 0x500, 0x600,
    ];
    const QUERY_MASKS: &[u32] = &[0x020, 0x040, 0x080, 0x100, 0x200, 0x400, 0x7e0];
    const HYDROGEN_MASKS: &[u32] = &[0x040, 0x080, 0x100, 0x1c0];

    // Generated with AvalonToolkit_2.0.5-pre.3 `SetFingerprintCountsWithFocus`
    // and `SetFingerprintBits` from `tmp/parity-audit/avalon_ring_state.c`.
    const NATIVE_GOLDENS: &str = r#"middle_matrix=mask_020_q0_dy0 result=4 counts=14:1,50:1,54:1,59:1, bits=0040000000004408
middle_matrix=mask_040_q0_dy0 result=0 counts= bits=0000000000000000
middle_matrix=mask_080_q0_dy0 result=0 counts= bits=0000000000000000
middle_matrix=mask_100_q0_dy0 result=8 counts=26:1,29:2,30:1,37:2,40:1,59:1, bits=0000006420010008
middle_matrix=mask_200_q0_dy0 result=36 counts=5:4,6:4,9:4,11:4,16:16,19:4,31:12,49:4,55:4,59:4,63:4, bits=600a098000008288
middle_matrix=mask_400_q0_dy0 result=1 counts=33:1, bits=0000000002000000
middle_matrix=mask_060_q0_dy0 result=4 counts=14:1,50:1,54:1,59:1, bits=0040000000004408
middle_matrix=mask_0a0_q0_dy0 result=4 counts=14:1,50:1,54:1,59:1, bits=0040000000004408
middle_matrix=mask_120_q0_dy0 result=12 counts=14:1,26:1,29:2,30:1,37:2,40:1,50:1,54:1,59:2, bits=0040006420014408
middle_matrix=mask_220_q0_dy0 result=40 counts=5:4,6:4,9:4,11:4,14:1,16:16,19:4,31:12,49:4,50:1,54:1,55:4,59:5,63:4, bits=604a09800000c688
middle_matrix=mask_420_q0_dy0 result=5 counts=14:1,33:1,50:1,54:1,59:1, bits=0040000002004408
middle_matrix=mask_0c0_q0_dy0 result=0 counts= bits=0000000000000000
middle_matrix=mask_140_q0_dy0 result=8 counts=26:1,29:2,30:1,37:2,40:1,59:1, bits=0000006420010008
middle_matrix=mask_240_q0_dy0 result=36 counts=5:4,6:4,9:4,11:4,16:16,19:4,31:12,49:4,55:4,59:4,63:4, bits=600a098000008288
middle_matrix=mask_440_q0_dy0 result=1 counts=33:1, bits=0000000002000000
middle_matrix=mask_180_q0_dy0 result=8 counts=26:1,29:2,30:1,37:2,40:1,59:1, bits=0000006420010008
middle_matrix=mask_280_q0_dy0 result=36 counts=5:4,6:4,9:4,11:4,16:16,19:4,31:12,49:4,55:4,59:4,63:4, bits=600a098000008288
middle_matrix=mask_480_q0_dy0 result=1 counts=33:1, bits=0000000002000000
middle_matrix=mask_300_q0_dy0 result=44 counts=5:4,6:4,9:4,11:4,16:16,19:4,26:1,29:2,30:1,31:12,37:2,40:1,49:4,55:4,59:5,63:4, bits=600a09e420018288
middle_matrix=mask_500_q0_dy0 result=9 counts=26:1,29:2,30:1,33:1,37:2,40:1,59:1, bits=0000006422010008
middle_matrix=mask_600_q0_dy0 result=37 counts=5:4,6:4,9:4,11:4,16:16,19:4,31:12,33:1,49:4,55:4,59:4,63:4, bits=600a098002008288
middle_matrix=mask_020_q1_dy0 result=2 counts=14:1,59:1, bits=0040000000000008
middle_matrix=mask_040_q1_dy0 result=0 counts= bits=0000000000000000
middle_matrix=mask_080_q1_dy0 result=0 counts= bits=0000000000000000
middle_matrix=mask_100_q1_dy0 result=0 counts= bits=0000000000000000
middle_matrix=mask_200_q1_dy0 result=36 counts=5:4,6:4,9:4,11:4,16:16,19:4,31:12,49:4,55:4,59:4,63:4, bits=600a098000008288
middle_matrix=mask_400_q1_dy0 result=1 counts=33:1, bits=0000000002000000
middle_matrix=mask_060_q1_dy0 result=2 counts=14:1,59:1, bits=0040000000000008
middle_matrix=mask_0a0_q1_dy0 result=2 counts=14:1,59:1, bits=0040000000000008
middle_matrix=mask_120_q1_dy0 result=2 counts=14:1,59:1, bits=0040000000000008
middle_matrix=mask_220_q1_dy0 result=38 counts=5:4,6:4,9:4,11:4,14:1,16:16,19:4,31:12,49:4,55:4,59:5,63:4, bits=604a098000008288
middle_matrix=mask_420_q1_dy0 result=3 counts=14:1,33:1,59:1, bits=0040000002000008
middle_matrix=mask_0c0_q1_dy0 result=0 counts= bits=0000000000000000
middle_matrix=mask_140_q1_dy0 result=0 counts= bits=0000000000000000
middle_matrix=mask_240_q1_dy0 result=36 counts=5:4,6:4,9:4,11:4,16:16,19:4,31:12,49:4,55:4,59:4,63:4, bits=600a098000008288
middle_matrix=mask_440_q1_dy0 result=1 counts=33:1, bits=0000000002000000
middle_matrix=mask_180_q1_dy0 result=0 counts= bits=0000000000000000
middle_matrix=mask_280_q1_dy0 result=36 counts=5:4,6:4,9:4,11:4,16:16,19:4,31:12,49:4,55:4,59:4,63:4, bits=600a098000008288
middle_matrix=mask_480_q1_dy0 result=1 counts=33:1, bits=0000000002000000
middle_matrix=mask_300_q1_dy0 result=36 counts=5:4,6:4,9:4,11:4,16:16,19:4,31:12,49:4,55:4,59:4,63:4, bits=600a098000008288
middle_matrix=mask_500_q1_dy0 result=1 counts=33:1, bits=0000000002000000
middle_matrix=mask_600_q1_dy0 result=37 counts=5:4,6:4,9:4,11:4,16:16,19:4,31:12,33:1,49:4,55:4,59:4,63:4, bits=600a098002008288
middle_matrix=mask_020_q0_dy1 result=4 counts=13:1,22:1,24:1,33:1, bits=0020400102000000
middle_matrix=mask_040_q0_dy1 result=0 counts= bits=0000000000000000
middle_matrix=mask_080_q0_dy1 result=0 counts= bits=0000000000000000
middle_matrix=mask_100_q0_dy1 result=8 counts=2:1,9:2,26:1,39:1,44:1,53:2, bits=0402000480102000
middle_matrix=mask_200_q0_dy1 result=64 counts=0:2,1:2,2:2,4:4,7:2,8:6,9:2,10:2,12:2,16:2,17:2,18:2,20:4,24:4,26:2,30:4,32:6,33:2,36:6,38:4,40:4,41:6,44:6,45:2,48:2,49:6,53:6,58:6,59:2, bits=971717455333230c
middle_matrix=mask_400_q0_dy1 result=1 counts=58:1, bits=0000000000000004
middle_matrix=mask_060_q0_dy1 result=4 counts=13:1,22:1,24:1,33:1, bits=0020400102000000
middle_matrix=mask_0a0_q0_dy1 result=4 counts=13:1,22:1,24:1,33:1, bits=0020400102000000
middle_matrix=mask_120_q0_dy1 result=12 counts=2:1,9:2,13:1,22:1,24:1,26:1,33:1,39:1,44:1,53:2, bits=0422400582102000
middle_matrix=mask_220_q0_dy1 result=68 counts=0:2,1:2,2:2,4:4,7:2,8:6,9:2,10:2,12:2,13:1,16:2,17:2,18:2,20:4,22:1,24:5,26:2,30:4,32:6,33:3,36:6,38:4,40:4,41:6,44:6,45:2,48:2,49:6,53:6,58:6,59:2, bits=973757455333230c
middle_matrix=mask_420_q0_dy1 result=5 counts=13:1,22:1,24:1,33:1,58:1, bits=0020400102000004
middle_matrix=mask_0c0_q0_dy1 result=0 counts= bits=0000000000000000
middle_matrix=mask_140_q0_dy1 result=8 counts=2:1,9:2,26:1,39:1,44:1,53:2, bits=0402000480102000
middle_matrix=mask_240_q0_dy1 result=64 counts=0:2,1:2,2:2,4:4,7:2,8:6,9:2,10:2,12:2,16:2,17:2,18:2,20:4,24:4,26:2,30:4,32:6,33:2,36:6,38:4,40:4,41:6,44:6,45:2,48:2,49:6,53:6,58:6,59:2, bits=971717455333230c
middle_matrix=mask_440_q0_dy1 result=1 counts=58:1, bits=0000000000000004
middle_matrix=mask_180_q0_dy1 result=8 counts=2:1,9:2,26:1,39:1,44:1,53:2, bits=0402000480102000
middle_matrix=mask_280_q0_dy1 result=64 counts=0:2,1:2,2:2,4:4,7:2,8:6,9:2,10:2,12:2,16:2,17:2,18:2,20:4,24:4,26:2,30:4,32:6,33:2,36:6,38:4,40:4,41:6,44:6,45:2,48:2,49:6,53:6,58:6,59:2, bits=971717455333230c
middle_matrix=mask_480_q0_dy1 result=1 counts=58:1, bits=0000000000000004
middle_matrix=mask_300_q0_dy1 result=72 counts=0:2,1:2,2:3,4:4,7:2,8:6,9:4,10:2,12:2,16:2,17:2,18:2,20:4,24:4,26:3,30:4,32:6,33:2,36:6,38:4,39:1,40:4,41:6,44:7,45:2,48:2,49:6,53:8,58:6,59:2, bits=97171745d333230c
middle_matrix=mask_500_q0_dy1 result=9 counts=2:1,9:2,26:1,39:1,44:1,53:2,58:1, bits=0402000480102004
middle_matrix=mask_600_q0_dy1 result=65 counts=0:2,1:2,2:2,4:4,7:2,8:6,9:2,10:2,12:2,16:2,17:2,18:2,20:4,24:4,26:2,30:4,32:6,33:2,36:6,38:4,40:4,41:6,44:6,45:2,48:2,49:6,53:6,58:7,59:2, bits=971717455333230c
middle_matrix=mask_020_q1_dy1 result=2 counts=22:1,33:1, bits=0000400002000000
middle_matrix=mask_040_q1_dy1 result=0 counts= bits=0000000000000000
middle_matrix=mask_080_q1_dy1 result=0 counts= bits=0000000000000000
middle_matrix=mask_100_q1_dy1 result=0 counts= bits=0000000000000000
middle_matrix=mask_200_q1_dy1 result=64 counts=0:2,1:2,2:2,4:4,7:2,8:6,9:2,10:2,12:2,16:2,17:2,18:2,20:4,24:4,26:2,30:4,32:6,33:2,36:6,38:4,40:4,41:6,44:6,45:2,48:2,49:6,53:6,58:6,59:2, bits=971717455333230c
middle_matrix=mask_400_q1_dy1 result=1 counts=58:1, bits=0000000000000004
middle_matrix=mask_060_q1_dy1 result=2 counts=22:1,33:1, bits=0000400002000000
middle_matrix=mask_0a0_q1_dy1 result=2 counts=22:1,33:1, bits=0000400002000000
middle_matrix=mask_120_q1_dy1 result=2 counts=22:1,33:1, bits=0000400002000000
middle_matrix=mask_220_q1_dy1 result=66 counts=0:2,1:2,2:2,4:4,7:2,8:6,9:2,10:2,12:2,16:2,17:2,18:2,20:4,22:1,24:4,26:2,30:4,32:6,33:3,36:6,38:4,40:4,41:6,44:6,45:2,48:2,49:6,53:6,58:6,59:2, bits=971757455333230c
middle_matrix=mask_420_q1_dy1 result=3 counts=22:1,33:1,58:1, bits=0000400002000004
middle_matrix=mask_0c0_q1_dy1 result=0 counts= bits=0000000000000000
middle_matrix=mask_140_q1_dy1 result=0 counts= bits=0000000000000000
middle_matrix=mask_240_q1_dy1 result=64 counts=0:2,1:2,2:2,4:4,7:2,8:6,9:2,10:2,12:2,16:2,17:2,18:2,20:4,24:4,26:2,30:4,32:6,33:2,36:6,38:4,40:4,41:6,44:6,45:2,48:2,49:6,53:6,58:6,59:2, bits=971717455333230c
middle_matrix=mask_440_q1_dy1 result=1 counts=58:1, bits=0000000000000004
middle_matrix=mask_180_q1_dy1 result=0 counts= bits=0000000000000000
middle_matrix=mask_280_q1_dy1 result=64 counts=0:2,1:2,2:2,4:4,7:2,8:6,9:2,10:2,12:2,16:2,17:2,18:2,20:4,24:4,26:2,30:4,32:6,33:2,36:6,38:4,40:4,41:6,44:6,45:2,48:2,49:6,53:6,58:6,59:2, bits=971717455333230c
middle_matrix=mask_480_q1_dy1 result=1 counts=58:1, bits=0000000000000004
middle_matrix=mask_300_q1_dy1 result=64 counts=0:2,1:2,2:2,4:4,7:2,8:6,9:2,10:2,12:2,16:2,17:2,18:2,20:4,24:4,26:2,30:4,32:6,33:2,36:6,38:4,40:4,41:6,44:6,45:2,48:2,49:6,53:6,58:6,59:2, bits=971717455333230c
middle_matrix=mask_500_q1_dy1 result=1 counts=58:1, bits=0000000000000004
middle_matrix=mask_600_q1_dy1 result=65 counts=0:2,1:2,2:2,4:4,7:2,8:6,9:2,10:2,12:2,16:2,17:2,18:2,20:4,24:4,26:2,30:4,32:6,33:2,36:6,38:4,40:4,41:6,44:6,45:2,48:2,49:6,53:6,58:7,59:2, bits=971717455333230c
middle_query=mask_020_q0_dy0 result=1 counts=59:1, bits=0000000000000008
middle_query=mask_040_q0_dy0 result=4 counts=22:1,24:1,49:1,51:1, bits=0000400100000a00
middle_query=mask_080_q0_dy0 result=1 counts=33:1, bits=0000000002000000
middle_query=mask_100_q0_dy0 result=0 counts= bits=0000000000000000
middle_query=mask_200_q0_dy0 result=0 counts= bits=0000000000000000
middle_query=mask_400_q0_dy0 result=0 counts= bits=0000000000000000
middle_query=mask_7e0_q0_dy0 result=6 counts=22:1,24:1,33:1,49:1,51:1,59:1, bits=0000400102000a08
middle_query=mask_020_q1_dy0 result=1 counts=59:1, bits=0000000000000008
middle_query=mask_040_q1_dy0 result=4 counts=22:1,24:1,49:1,51:1, bits=0000400100000a00
middle_query=mask_080_q1_dy0 result=0 counts= bits=0000000000000000
middle_query=mask_100_q1_dy0 result=0 counts= bits=0000000000000000
middle_query=mask_200_q1_dy0 result=0 counts= bits=0000000000000000
middle_query=mask_400_q1_dy0 result=0 counts= bits=0000000000000000
middle_query=mask_7e0_q1_dy0 result=5 counts=22:1,24:1,49:1,51:1,59:1, bits=0000400100000a08
middle_query=mask_020_q0_dy1 result=1 counts=59:1, bits=0000000000000008
middle_query=mask_040_q0_dy1 result=4 counts=22:1,24:1,49:1,51:1, bits=0000400100000a00
middle_query=mask_080_q0_dy1 result=1 counts=33:1, bits=0000000002000000
middle_query=mask_100_q0_dy1 result=0 counts= bits=0000000000000000
middle_query=mask_200_q0_dy1 result=0 counts= bits=0000000000000000
middle_query=mask_400_q0_dy1 result=0 counts= bits=0000000000000000
middle_query=mask_7e0_q0_dy1 result=6 counts=22:1,24:1,33:1,49:1,51:1,59:1, bits=0000400102000a08
middle_query=mask_020_q1_dy1 result=1 counts=59:1, bits=0000000000000008
middle_query=mask_040_q1_dy1 result=4 counts=22:1,24:1,49:1,51:1, bits=0000400100000a00
middle_query=mask_080_q1_dy1 result=0 counts= bits=0000000000000000
middle_query=mask_100_q1_dy1 result=0 counts= bits=0000000000000000
middle_query=mask_200_q1_dy1 result=0 counts= bits=0000000000000000
middle_query=mask_400_q1_dy1 result=0 counts= bits=0000000000000000
middle_query=mask_7e0_q1_dy1 result=5 counts=22:1,24:1,49:1,51:1,59:1, bits=0000400100000a08
middle_h=implicit_mask_040_dy0 result=1 counts=6:1, bits=4000000000000000
middle_h=implicit_mask_080_dy0 result=1 counts=30:1, bits=0000004000000000
middle_h=implicit_mask_100_dy0 result=2 counts=42:2, bits=0000000000040000
middle_h=implicit_mask_1c0_dy0 result=4 counts=6:1,30:1,42:2, bits=4000004000040000
middle_h=explicit_mask_040_dy0 result=1 counts=6:1, bits=4000000000000000
middle_h=explicit_mask_080_dy0 result=1 counts=30:1, bits=0000004000000000
middle_h=explicit_mask_100_dy0 result=2 counts=42:2, bits=0000000000040000
middle_h=explicit_mask_1c0_dy0 result=4 counts=6:1,30:1,42:2, bits=4000004000040000
middle_h=implicit_mask_040_dy1 result=1 counts=6:1, bits=4000000000000000
middle_h=implicit_mask_080_dy1 result=1 counts=30:1, bits=0000004000000000
middle_h=implicit_mask_100_dy1 result=2 counts=42:2, bits=0000000000040000
middle_h=implicit_mask_1c0_dy1 result=4 counts=6:1,30:1,42:2, bits=4000004000040000
middle_h=explicit_mask_040_dy1 result=1 counts=6:1, bits=4000000000000000
middle_h=explicit_mask_080_dy1 result=1 counts=30:1, bits=0000004000000000
middle_h=explicit_mask_100_dy1 result=2 counts=42:2, bits=0000000000040000
middle_h=explicit_mask_1c0_dy1 result=4 counts=6:1,30:1,42:2, bits=4000004000040000"#;

    fn matrix_fixture() -> MoleculeState {
        let symbols = ["N", "C", "C", "C", "O", "C"];
        let endpoints = [[1, 2], [2, 3], [3, 4], [4, 1], [3, 5], [5, 6], [6, 4]];
        let bond_types = [2, 1, 1, 1, 2, 1, 2];
        MoleculeState {
            atoms: symbols
                .map(|symbol| Atom {
                    atom_symbol: symbol.to_string(),
                    ..Atom::default()
                })
                .to_vec(),
            bonds: endpoints
                .into_iter()
                .zip(bond_types)
                .map(|(atoms, bond_type)| Bond {
                    atoms,
                    bond_type,
                    ..Bond::default()
                })
                .collect(),
            ..MoleculeState::default()
        }
    }

    fn query_fixture() -> MoleculeState {
        let mut nitrogen = Atom {
            atom_symbol: "N".to_string(),
            ..Atom::default()
        };
        nitrogen.sub_desc = SUB_MORE;
        nitrogen.query_h_count = 3;
        let mut oxygen = Atom {
            atom_symbol: "O".to_string(),
            ..Atom::default()
        };
        oxygen.sub_desc = SUB_AS_IS;
        oxygen.query_h_count = 1;
        MoleculeState {
            atoms: vec![
                nitrogen,
                oxygen,
                Atom {
                    atom_symbol: "C".to_string(),
                    ..Atom::default()
                },
            ],
            bonds: vec![
                Bond {
                    atoms: [1, 2],
                    bond_type: 1,
                    ..Bond::default()
                },
                Bond {
                    atoms: [2, 3],
                    bond_type: 1,
                    ..Bond::default()
                },
            ],
            ..MoleculeState::default()
        }
    }

    fn hydrogen_fixture(explicit_hydrogen: bool) -> MoleculeState {
        let mut atoms = ["C", "C", "O"]
            .map(|symbol| Atom {
                atom_symbol: symbol.to_string(),
                ..Atom::default()
            })
            .to_vec();
        let mut bonds = vec![
            Bond {
                atoms: [1, 2],
                bond_type: 1,
                ..Bond::default()
            },
            Bond {
                atoms: [2, 3],
                bond_type: 1,
                ..Bond::default()
            },
        ];
        if explicit_hydrogen {
            atoms.push(Atom {
                atom_symbol: "H".to_string(),
                ..Atom::default()
            });
            bonds.push(Bond {
                atoms: [3, 4],
                bond_type: 1,
                ..Bond::default()
            });
        }
        MoleculeState {
            atoms,
            bonds,
            ..MoleculeState::default()
        }
    }

    fn sparse_text(counts: &[i32]) -> String {
        counts
            .iter()
            .copied()
            .enumerate()
            .filter(|&(_, count)| count != 0)
            .map(|(index, count)| format!("{index}:{count},"))
            .collect()
    }

    fn packed_bit_hex(counts: &[i32]) -> String {
        let mut bytes = vec![0_u8; counts.len().div_ceil(8)];
        for (index, &count) in counts.iter().enumerate() {
            if count > 0 {
                bytes[index / 8] |= 1 << (index % 8);
            }
        }
        bytes
            .into_iter()
            .map(|byte| format!("{byte:02x}"))
            .collect()
    }

    fn run_case(
        mut molecule: MoleculeState,
        mask: u32,
        as_query: bool,
        daylight: bool,
    ) -> (i32, String, String) {
        let mut counts = vec![0_i32; 64];
        let result = count_middle_flag_families(
            &mut molecule,
            &mut counts,
            AvalonFingerprintFlags::from_bits(mask).unwrap(),
            as_query,
            i32::from(daylight) * USE_DY_AROMATICITY,
            0,
        )
        .unwrap();
        (result, sparse_text(&counts), packed_bit_hex(&counts))
    }

    fn expected_lines(prefix: &str) -> impl Iterator<Item = &'static str> + '_ {
        NATIVE_GOLDENS
            .lines()
            .filter(move |line| line.starts_with(prefix))
    }

    #[test]
    fn middle_single_and_pairwise_matrix_matches_native_counts_and_bits() {
        assert_eq!(MIDDLE_MASKS.len(), 21);
        assert_eq!(
            MIDDLE_MASKS
                .iter()
                .filter(|mask| mask.count_ones() == 1)
                .count(),
            6
        );
        assert_eq!(
            MIDDLE_MASKS
                .iter()
                .filter(|mask| mask.count_ones() == 2)
                .count(),
            15
        );
        let mut expected = expected_lines("middle_matrix=");
        for daylight in [false, true] {
            for as_query in [false, true] {
                for &mask in MIDDLE_MASKS {
                    let (result, counts, bits) =
                        run_case(matrix_fixture(), mask, as_query, daylight);
                    let actual = format!(
                        "middle_matrix=mask_{mask:03x}_q{}_dy{} result={result} counts={counts} bits={bits}",
                        u8::from(as_query),
                        u8::from(daylight),
                    );
                    assert_eq!(Some(actual.as_str()), expected.next(), "{actual}");
                }
            }
        }
        assert_eq!(expected.next(), None);
    }

    #[test]
    fn middle_query_branches_match_native_counts_and_bits() {
        let mut expected = expected_lines("middle_query=");
        for daylight in [false, true] {
            for as_query in [false, true] {
                for &mask in QUERY_MASKS {
                    let (result, counts, bits) =
                        run_case(query_fixture(), mask, as_query, daylight);
                    let actual = format!(
                        "middle_query=mask_{mask:03x}_q{}_dy{} result={result} counts={counts} bits={bits}",
                        u8::from(as_query),
                        u8::from(daylight),
                    );
                    assert_eq!(Some(actual.as_str()), expected.next(), "{actual}");
                }
            }
        }
        assert_eq!(expected.next(), None);
    }

    #[test]
    fn explicit_and_implicit_hydrogen_paths_match_native_counts_and_bits() {
        let mut expected = expected_lines("middle_h=");
        for daylight in [false, true] {
            for explicit_hydrogen in [false, true] {
                for &mask in HYDROGEN_MASKS {
                    let (result, counts, bits) =
                        run_case(hydrogen_fixture(explicit_hydrogen), mask, false, daylight);
                    let fixture = if explicit_hydrogen {
                        "explicit"
                    } else {
                        "implicit"
                    };
                    let actual = format!(
                        "middle_h={fixture}_mask_{mask:03x}_dy{} result={result} counts={counts} bits={bits}",
                        u8::from(daylight),
                    );
                    assert_eq!(Some(actual.as_str()), expected.next(), "{actual}");
                }
            }
        }
        assert_eq!(expected.next(), None);
    }
}
