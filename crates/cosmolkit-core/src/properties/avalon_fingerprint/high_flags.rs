//! Source-order implementation of Avalon fingerprint flags `0x000800..=0x004000`.

use crate::FingerprintError;

use super::AvalonFingerprintFlags;
use super::fingerprint_state::{FingerprintPreprocessingState, avalon_atomic_number, with_prepared_fingerprint_state};
use super::hash::next_hash;
use super::reaccs::MoleculeState;
use super::symbols::atom_symbol_match;
use super::traversal::{
    CLASS_SPIDER_SEED, DEGREE_PATH_SEED, PathFlags, RING_SIZE_SEED, add_bit, add_bit_count, build_path_length_matrix,
    collect_special_neighbours, set_feature_bits, set_path_bits_rec,
};

const SUB_AS_IS: i32 = -2;
const SUB_ONE: i32 = 1;
const SUB_MORE: i32 = 6;
const ANY_BOND: i32 = 8;
const ACCUMULATE_BITS: i32 = 0x0002;
const MAX_SPIDER: usize = 7;
const CSP3: i32 = 19;
const HETERO: i32 = 23;
const GENERIC: i32 = -1;
const SPECIAL_RING: u32 = 0xfc & !(1 << 6);
const HETERO_FLAG: i32 = 0x0100;
const RING_SUBST_FLAG: i32 = 0x0200;
const QUART_FLAG: i32 = 0x0400;
const CSP3_FLAG: i32 = 0x0800;
const RS_SPECIAL_FLAG: i32 = 0x1000;
const C_FLAG: i32 = 0x0001;
const O_FLAG: i32 = 0x0002;
const N_FLAG: i32 = 0x0003;
const S_FLAG: i32 = 0x0004;
const P_FLAG: i32 = 0x0005;
const X_FLAG: i32 = 0x0006;

pub(super) fn count_high_flag_families(
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
            Ok(count_high_flag_families_prepared(
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

pub(super) fn count_high_flag_families_prepared(
    working: &mut MoleculeState,
    state: &FingerprintPreprocessingState,
    counts: &mut [i32],
    which_bits: AvalonFingerprintFlags,
    as_query: bool,
    exclude_atom: i32,
) -> i32 {
    let mut result = 0_i32;
    if which_bits.contains(AvalonFingerprintFlags::RING_SIZE_COUNTS) {
        result += count_ring_size_counts(working, counts, exclude_atom);
    }

    prepare_degree_path_state(working, state, as_query, exclude_atom);
    if which_bits.contains(AvalonFingerprintFlags::DEGREE_PATH) {
        result += count_degree_paths(working, state, counts, exclude_atom);
    }

    // Avalon❗✔️: // fprintf(stderr, "USE_CLASS_SPIDERS | USE_FEATURE_PAIRS | USE_NON_SSS_BITS\n");
    // Avalon❗✔️:    if (which_bits & (USE_CLASS_SPIDERS | USE_FEATURE_PAIRS | USE_NON_SSS_BITS))
    // Avalon❗✔️:    {
    // Avalon❗✔️:       /* Collect length_matrix */
    // Avalon❗✔️:       /* allocate storage length_matrix */
    // Avalon❗✔️:       length_tmp = TypeAlloc(mp->n_atoms*mp->n_atoms, int);
    // Avalon❗✔️:       /* allocat indices */
    // Avalon❗✔️:       length_matrix = TypeAlloc(mp->n_atoms, int*);
    // Avalon❗✔️:       for (i=0; i<mp->n_atoms; i++)
    // Avalon❗✔️:       {
    // Avalon❗✔️:          /* set relative pointers */
    // Avalon❗✔️:          length_matrix[i] = length_tmp+i*mp->n_atoms;
    // Avalon❗✔️:          touched_indices[i] = 0;
    // Avalon❗✔️:       }
    // Avalon❗✔️:       for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    // Avalon❗✔️:       {
    // Avalon❗✔️:           if (i+1 == exclude_atom) continue;
    // Avalon❗✔️:           touched_indices[i] = 1; /* updating */
    // Avalon❗✔️: // if (FALSE) fprintf(stderr, "starting path search at atom %d(%d)\n", i+1, ap->color);
    // Avalon❗✔️:           SetPathLengthFlags(mp, touched_indices,
    // Avalon❗✔️:                              i, 0, i,
    // Avalon❗✔️:                              12,      /* path perception distance <= 12 */
    // Avalon❗✔️:                              length_matrix, nbp,
    // Avalon❗✔️:                              exclude_atom,
    // Avalon❗✔️: 			     prefix);
    // Avalon❗✔️:           touched_indices[i] = 0; /* down-dating */
    // Avalon❗✔️:       }
    // Avalon❗✔️:    }
    let feature_mask = AvalonFingerprintFlags::CLASS_SPIDERS | AvalonFingerprintFlags::FEATURE_PAIRS;
    let feature_active = which_bits.bits() & feature_mask.bits() != 0;
    let length_matrix = feature_active.then(|| build_path_length_matrix(working, &state.neighbours, 12, exclude_atom));

    if feature_active {
        prepare_class_spider_state(working, state, exclude_atom);
        if which_bits.contains(AvalonFingerprintFlags::CLASS_SPIDERS) {
            result += count_class_spiders(working, state, counts, exclude_atom);
        }
        if which_bits.contains(AvalonFingerprintFlags::FEATURE_PAIRS) {
            prepare_feature_pair_state(working, state, exclude_atom);
            result += count_feature_pairs(
                working,
                counts,
                length_matrix
                    .as_deref()
                    .expect("feature flags allocate the source length matrix"),
                exclude_atom,
            );
        }
    }
    result
}
fn count_ring_size_counts(molecule: &MoleculeState, counts: &mut [i32], exclude_atom: i32) -> i32 {
    // Avalon❗✔️: // fprintf(stderr, "USE_RING_SIZE_COUNTS\n");
    // Avalon❗✔️:    if (which_bits & USE_RING_SIZE_COUNTS)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       strcpy(prefix, "RSC:");
    let mut result = 0_i32;
    // Avalon❗✔️:       for (j=3; j<10; j++)      /* loop through ring_sizes */
    // Avalon❗✔️:       {
    for ring_size in 3_i32..10 {
        // Avalon❗✔️:          nrbonds = 0;
        let mut ring_bonds = 0_i32;
        // Avalon❗✔️:          for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
        // Avalon❗✔️:             if (bp->atoms[0] != exclude_atom  &&
        // Avalon❗✔️:                 bp->atoms[1] != exclude_atom  &&
        // Avalon❗✔️:                 (bp->rsize_flags&(1<<j))) nrbonds++;
        for bond in &molecule.bonds {
            if bond.atoms[0] != exclude_atom
                && bond.atoms[1] != exclude_atom
                && bond.rsize_flags & (1_u32 << ring_size) != 0
            {
                ring_bonds += 1;
            }
        }
        // Avalon❗✔️:          seed = RING_SIZE_SEED;
        // Avalon❗✔️:          seed = NEXT_SEED(seed, j*13);
        let mut seed = next_hash(RING_SIZE_SEED, (ring_size * 13) as u64);
        // Avalon❗✔️:          for (i=1; i<100; i*=2)
        // Avalon❗✔️:          {
        let mut level = 1_i32;
        while level < 100 {
            // Avalon❗✔️:             if (nrbonds >= j*i)
            // Avalon❗✔️:             {
            if ring_bonds >= ring_size * level {
                // Avalon❗✔️:                seed = NEXT_SEED(seed, i);
                seed = next_hash(seed, level as u64);
                // Avalon❗✔️:                /* don't set a bit if just one 5- or 6-ring */
                // Avalon❗✔️:                if ((j != 6  &&  j != 5)  || nrbonds > j*i)
                // Avalon❗✔️: 	       {
                if (ring_size != 6 && ring_size != 5) || ring_bonds > ring_size * level {
                    // Avalon❗✔️:                   ADD_BIT_COUNT(fp_counts, ncounts, seed, nrbonds);
                    add_bit_count(counts, seed, ring_bonds);
                    // Avalon❗✔️: 		  result++;
                    result += 1;
                    // Avalon❗✔️: 	       }
                }
                // Avalon❗✔️:                if (j != 6  &&  j != 5)    /* set one more bit for odd rings */
                // Avalon❗✔️:                {
                if ring_size != 6 && ring_size != 5 {
                    // Avalon❗✔️:                   seed = NEXT_SEED(seed, 17);
                    seed = next_hash(seed, 17);
                    // Avalon❗✔️:                   ADD_BIT_COUNT(fp_counts, ncounts, seed, nrbonds);
                    add_bit_count(counts, seed, ring_bonds);
                    // Avalon❗✔️: 		  result++;
                    result += 1;
                    // Avalon❗✔️:                }
                }
                // Avalon❗✔️:             }
                level *= 2;
                // Avalon❗✔️:             else
                // Avalon❗✔️:                break;
            } else {
                break;
            }
            // Avalon❗✔️:          }
        }
        // Avalon❗✔️:       }
    }

    // Avalon❗✔️:       /* Set bits for different ring sizes connected by a bond */
    // Avalon❗✔️:       for (j=0; j<15; j++)      /* loop through ring_sizes */
    // Avalon❗✔️:          for (k=0; k<15; k++)      /* loop through ring_sizes */
    // Avalon❗✔️:             rscounts[j][k] = 0;
    let mut ring_size_counts = [[0_i32; 15]; 15];
    // Avalon❗✔️:       /* count connections */
    // Avalon❗✔️:       for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
    // Avalon❗✔️:       {
    for bond in &molecule.bonds {
        // Avalon❗✔️:          if (bp->atoms[0] == exclude_atom) continue;
        // Avalon❗✔️:          if (bp->atoms[1] == exclude_atom) continue;
        if bond.atoms.contains(&exclude_atom) {
            continue;
        }
        let first = bond.atoms[0] as usize - 1;
        let second = bond.atoms[1] as usize - 1;
        // Avalon❗✔️:          for (j=3; j<15; j++)      /* loop through ring_sizes */
        // Avalon❗✔️:             for (k=3; k<15; k++)      /* loop through ring_sizes */
        // Avalon❗✔️:             {
        for first_size in 3_usize..15 {
            for second_size in 3_usize..15 {
                // Avalon❗✔️:                if (j==k) continue;
                if first_size == second_size {
                    continue;
                }
                // Avalon❗✔️:                if ((mp->atom_array[bp->atoms[0]-1].rsize_flags&(1<<j))  &&
                // Avalon❗✔️:                    (mp->atom_array[bp->atoms[1]-1].rsize_flags&(1<<k)))
                // Avalon❗✔️:                {
                if molecule.atoms[first].rsize_flags & (1_u32 << first_size) != 0
                    && molecule.atoms[second].rsize_flags & (1_u32 << second_size) != 0
                {
                    // Avalon❗✔️:                   rscounts[j][k]++;
                    ring_size_counts[first_size][second_size] += 1;
                    // Avalon❗✔️:                   rscounts[k][j]++;
                    ring_size_counts[second_size][first_size] += 1;
                    // Avalon❗✔️:                }
                }
                // Avalon❗✔️:             }
            }
        }
        // Avalon❗✔️:       }
    }
    // Avalon❗✔️:       /* set bits */
    // Avalon❗✔️:       for (j=3; j<9; j++)      /* loop through not too large ring_sizes */
    // Avalon❗✔️:          for (k=j+1; k<9; k++)      /* loop through not too large ring_sizes */
    // Avalon❗✔️:          {
    for first_size in 3_usize..9 {
        for second_size in first_size + 1..9 {
            // Avalon❗✔️:             if (rscounts[j][k] == 0) continue;
            if ring_size_counts[first_size][second_size] == 0 {
                continue;
            }
            // Avalon❗✔️:             seed = 2*RING_SIZE_SEED;
            // Avalon❗✔️:             seed = NEXT_SEED(seed, (j+k)*11);
            // Avalon❗✔️:             seed = NEXT_SEED(seed, j*k*13);
            // Avalon❗✔️:             seed = NEXT_SEED(seed, 19);
            let mut seed = next_hash(2 * RING_SIZE_SEED, ((first_size + second_size) * 11) as u64);
            seed = next_hash(seed, (first_size * second_size * 13) as u64);
            seed = next_hash(seed, 19);
            // Avalon❗✔️:             ADD_BIT(fp_counts, ncounts, seed);
            add_bit(counts, seed);
            // Avalon❗✔️: 	    result++;
            result += 1;
            // Avalon❗✔️:             /* more special links => more bits */
            // Avalon❗✔️:             if (rscounts[j][k]==1  ||  j==6  ||  k==6) continue;
            if ring_size_counts[first_size][second_size] == 1 || first_size == 6 || second_size == 6 {
                continue;
            }
            // Avalon❗✔️:             seed = 2*RING_SIZE_SEED+23;
            // Avalon❗✔️:             seed = NEXT_SEED(seed, (j+k)*11);
            // Avalon❗✔️:             seed = NEXT_SEED(seed, j*k*13);
            // Avalon❗✔️:             seed = NEXT_SEED(seed, 19);
            seed = next_hash(2 * RING_SIZE_SEED + 23, ((first_size + second_size) * 11) as u64);
            seed = next_hash(seed, (first_size * second_size * 13) as u64);
            seed = next_hash(seed, 19);
            // Avalon❗✔️:             ADD_BIT(fp_counts, ncounts, seed);
            add_bit(counts, seed);
            // Avalon❗✔️: 	    result++;
            result += 1;
            // Avalon❗✔️:          }
        }
    }
    // Avalon❗✔️:    }
    result
}

fn prepare_degree_path_state(
    molecule: &mut MoleculeState,
    state: &FingerprintPreprocessingState,
    as_query: bool,
    exclude_atom: i32,
) {
    // Avalon❗✔️:    /* Set the color property to represent all different atom types */
    // Avalon❗✔️:    for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    // Avalon❗✔️:    {
    for (atom_index, atom) in molecule.atoms.iter_mut().enumerate() {
        // Avalon❗✔️:       ap->color = StringToInt(periodic_table, ap->atom_symbol);
        atom.color = avalon_atomic_number(&atom.atom_symbol);
        // Avalon❗✔️:       if (ap->color <= 1) ap->color = 0;        /* ignore hydrogens */
        if atom.color <= 1 {
            atom.color = 0;
        }
        // Avalon❗✔️:       /* mark special atom types */
        // Avalon❗✔️:       if (ap->color > 115) ap->color = -1;
        if atom.color > 115 {
            atom.color = -1;
        }
        // Avalon❗✔️:       if (0 == strcmp("A", ap->atom_symbol)) ap->color = -1;
        if atom.atom_symbol == "A" {
            atom.color = -1;
        }
        // Avalon❗✔️:       if (i+1 == exclude_atom) ap->color = 0;
        if atom_index as i32 + 1 == exclude_atom {
            atom.color = 0;
        }
        // Avalon❗✔️:       if (ap->color > 1  &&
        // Avalon❗✔️:           (!as_query                 ||
        // Avalon❗✔️:            ap->sub_desc == SUB_AS_IS ||
        // Avalon❗✔️:            (ap->sub_desc != NONE     &&
        // Avalon❗✔️:             ap->sub_desc != SUB_MORE &&
        // Avalon❗✔️:             ap->sub_desc == degree[i]+SUB_ONE-1)))
        // Avalon❗✔️:       {
        if atom.color > 1
            && (!as_query
                || atom.sub_desc == SUB_AS_IS
                || (atom.sub_desc != 0
                    && atom.sub_desc != SUB_MORE
                    && atom.sub_desc == state.degree[atom_index] + SUB_ONE - 1))
        {
            // Avalon❗✔️:          ap->color += 32*degree[i];
            atom.color += 32 * state.degree[atom_index];
            // Avalon❗✔️:       }
            // Avalon❗✔️:       else
            // Avalon❗✔️:          ap->color = 0;
        } else {
            atom.color = 0;
        }
        // Avalon❗✔️:    }
    }
    // Avalon❗✔️:    for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
    // Avalon❗✔️:    {
    for bond in &mut molecule.bonds {
        // Avalon❗✔️:       if (SINGLE <= bp->bond_type  &&  bp->bond_type <= ANY_BOND)
        // Avalon❗✔️:          bp->color = 5;
        // Avalon❗✔️:       else
        // Avalon❗✔️:          bp->color = 0;
        bond.color = i32::from((1..=ANY_BOND).contains(&bond.bond_type)) * 5;
        // Avalon❗✔️:       if (bp->atoms[0] == exclude_atom) bp->color = 0;
        // Avalon❗✔️:       if (bp->atoms[1] == exclude_atom) bp->color = 0;
        if bond.atoms.contains(&exclude_atom) {
            bond.color = 0;
        }
        // Avalon❗✔️:    }
    }
}

fn count_degree_paths(
    molecule: &MoleculeState,
    state: &FingerprintPreprocessingState,
    counts: &mut [i32],
    exclude_atom: i32,
) -> i32 {
    // Avalon❗✔️: // fprintf(stderr, "USE_DEGREE_PATH\n");
    // Avalon❗✔️:    /*
    // Avalon❗✔️:     * add special bits for paths starting at atoms with >= 4 neighbours or
    // Avalon❗✔️:     * methyl atoms
    // Avalon❗✔️:     *
    // Avalon❗✔️:     * atom color represent degree bond color identical for all bonds to non-H
    // Avalon❗✔️:     */
    // Avalon❗✔️:    if (which_bits & USE_DEGREE_PATH)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       strcpy(prefix, "DP:");
    // Avalon❗✔️:       seed = DEGREE_PATH_SEED;
    let mut result = 0_i32;
    let seed = DEGREE_PATH_SEED;
    let mut touched_indices = vec![0_i32; molecule.atoms.len()];

    // Avalon❗✔️:       // set bits for degree paths starting with special carbon atoms
    // Avalon❗✔️:       for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    // Avalon❗✔️:       {
    for (atom_index, atom) in molecule.atoms.iter().enumerate() {
        // Avalon❗✔️:          if (ap->color <= 0) continue;  // only process if degree defined
        // Avalon❗✔️:          if (i+1 == exclude_atom) continue;
        if atom.color <= 0 || atom_index as i32 + 1 == exclude_atom {
            continue;
        }
        // Avalon❗✔️:          /* don't start on usual atoms */
        // Avalon❗✔️:          if (degree[i] <= 3  &&		// include if high degree
        // Avalon❗✔️:              atom_status[i] <= 2  &&	// include if ring fusion
        // Avalon❗✔️:              degree[i] != 1)		// include if terminal atom
        // Avalon❗✔️:             continue;
        if state.degree[atom_index] <= 3 && state.atom_status[atom_index] <= 2 && state.degree[atom_index] != 1 {
            continue;
        }
        // Avalon❗✔️:          /* only keep terminals if methyl */
        // Avalon❗✔️:          if (degree[i] == 1  &&  0 != strcmp(ap->atom_symbol, "C"))
        // Avalon❗✔️:             continue;
        if state.degree[atom_index] == 1 && atom.atom_symbol != "C" {
            continue;
        }
        // Avalon❗✔️: // fprintf(stderr,"degree of %s atom %d is %d, status is %d\n",
        // Avalon❗✔️: // ap->atom_symbol, i+1, degree[i], atom_status[i]);
        // Avalon❗✔️:          touched_indices[i] = 1; /* updating */
        touched_indices[atom_index] = 1;
        // Avalon❗✔️:          old_seed = seed;
        // Avalon❗✔️:          seed = NEXT_SEED(seed, ap->color);
        let start_seed = next_hash(seed, atom.color as u64);
        // Avalon❗✔️: // fprintf(stderr, "%s(%d) has degree %d(%d)\n",
        // Avalon❗✔️: // ap->atom_symbol, i+1, degree[i], ap->sub_desc);
        // Avalon❗✔️:          result += SetPathBitsRec(mp, nbp,
        // Avalon❗✔️:                                   fp_counts, ncounts,
        // Avalon❗✔️:                                   seed, touched_indices,
        // Avalon❗✔️:                                   1,
        // Avalon❗✔️:                                   2, 4,  /* path length 1 to 3 */
        // Avalon❗✔️:                                   i, 0, -1,
        // Avalon❗✔️:                                   // IGNORE_TERM_SYMBOL |
        // Avalon❗✔️:                                   IGNORE_PATH_SYMBOL |
        // Avalon❗✔️:                                   PROCESS_CHAINS,
        // Avalon❗✔️:                                   exclude_atom,
        // Avalon❗✔️: 				   prefix);
        result += set_path_bits_rec(
            molecule,
            &state.neighbours,
            counts,
            start_seed,
            &mut touched_indices,
            1,
            2,
            4,
            atom_index,
            0,
            -1,
            PathFlags::IGNORE_PATH_SYMBOL | PathFlags::PROCESS_CHAINS,
            exclude_atom,
        );
        // Avalon❗✔️:          /* special CH fusion atoms */
        // Avalon❗✔️:          if (atom_status[i] > 2  &&  H_count[i+1] >= 1)
        // Avalon❗✔️:          {
        if state.atom_status[atom_index] > 2 && state.hydrogen_counts[atom_index + 1] >= 1 {
            // Avalon❗✔️:             seed = NEXT_SEED(seed, 219);
            // Avalon❗✔️:             result += SetPathBitsRec(mp, nbp,
            // Avalon❗✔️:                                      fp_counts, ncounts,
            // Avalon❗✔️:                                      seed, touched_indices,
            // Avalon❗✔️:                                      1,
            // Avalon❗✔️:                                      2, 5,
            // Avalon❗✔️:                                      i, 0, -1,
            // Avalon❗✔️:                                      IGNORE_PATH_SYMBOL |
            // Avalon❗✔️:                                      PROCESS_CHAINS,
            // Avalon❗✔️: 				     exclude_atom,
            // Avalon❗✔️: 				     prefix);
            result += set_path_bits_rec(
                molecule,
                &state.neighbours,
                counts,
                next_hash(start_seed, 219),
                &mut touched_indices,
                1,
                2,
                5,
                atom_index,
                0,
                -1,
                PathFlags::IGNORE_PATH_SYMBOL | PathFlags::PROCESS_CHAINS,
                exclude_atom,
            );
            // Avalon❗✔️:          }
        }
        // Avalon❗✔️:          seed = old_seed;
        // The immutable source seed remains unchanged.
        // Avalon❗✔️:          touched_indices[i] = 0; /* down-dating */
        touched_indices[atom_index] = 0;
        // Avalon❗✔️:       }
    }

    // Avalon❗✔️:       // set bits for degree paths starting with hetero atoms
    // Avalon❗✔️:       for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    // Avalon❗✔️:       {
    for (atom_index, atom) in molecule.atoms.iter().enumerate() {
        // Avalon❗✔️:          if (i+1 == exclude_atom) continue;
        if atom_index as i32 + 1 == exclude_atom {
            continue;
        }
        // Avalon❗✔️:          tmp = StringToInt(periodic_table, ap->atom_symbol);
        let atomic_number = avalon_atomic_number(&atom.atom_symbol);
        // Avalon❗✔️:          if (1 >= tmp  ||  tmp >= 115) continue;
        // Avalon❗✔️:          if (tmp == 6) continue;
        // Avalon❗✔️:          if (tmp < 10)	continue; // exclude common hetero atoms
        if atomic_number <= 1 || atomic_number >= 115 || atomic_number == 6 || atomic_number < 10 {
            continue;
        }
        // Avalon❗✔️:          touched_indices[i] = 1; /* updating */
        touched_indices[atom_index] = 1;
        // Avalon❗✔️:          old_seed = seed;
        // Avalon❗✔️:          seed = NEXT_SEED(seed, tmp);
        let start_seed = next_hash(seed, atomic_number as u64);
        // Avalon❗✔️: // fprintf(stderr, "%s(%d) has degree %d(%d)\n",
        // Avalon❗✔️: // ap->atom_symbol, i+1, degree[i], ap->sub_desc);
        // Avalon❗✔️:          result += SetPathBitsRec(mp, nbp,
        // Avalon❗✔️:                                   fp_counts, ncounts,
        // Avalon❗✔️:                                   seed, touched_indices,
        // Avalon❗✔️:                                   1,
        // Avalon❗✔️:                                   2, 2,
        // Avalon❗✔️:                                   i, 0, -1,
        // Avalon❗✔️:                                   PROCESS_CHAINS,
        // Avalon❗✔️:                                   exclude_atom,
        // Avalon❗✔️: 				   prefix);
        result += set_path_bits_rec(
            molecule,
            &state.neighbours,
            counts,
            start_seed,
            &mut touched_indices,
            1,
            2,
            2,
            atom_index,
            0,
            -1,
            PathFlags::PROCESS_CHAINS,
            exclude_atom,
        );
        // Avalon❗✔️:          if (0) // might overly populate complexes
        // Avalon❗✔️:          result += SetPathBitsRec(mp, nbp,
        // Avalon❗✔️:                                   fp_counts, ncounts,
        // Avalon❗✔️:                                   101+seed, touched_indices,
        // Avalon❗✔️:                                   1,
        // Avalon❗✔️:                                   2, 2,
        // Avalon❗✔️:                                   i, 0, -1,
        // Avalon❗✔️:                                   FORCED_RING_PATH |
        // Avalon❗✔️:                                   PROCESS_CHAINS,
        // Avalon❗✔️:                                   exclude_atom,
        // Avalon❗✔️: 				   prefix);
        // The source branch is compile-time disabled.
        // Avalon❗✔️:          seed = old_seed;
        // The immutable source seed remains unchanged.
        // Avalon❗✔️:          touched_indices[i] = 0; /* down-dating */
        touched_indices[atom_index] = 0;
        // Avalon❗✔️:       }
    }
    // Avalon❗✔️:    }
    result
}

fn prepare_class_spider_state(molecule: &mut MoleculeState, state: &FingerprintPreprocessingState, exclude_atom: i32) {
    // Avalon❗✔️:    /*
    // Avalon❗✔️:     * This screen class will catch non-linear fragments composed of rather
    // Avalon❗✔️:     * frequent linear sub-fragments
    // Avalon❗✔️:     */
    // Avalon❗✔️:    if (which_bits & (USE_CLASS_SPIDERS | USE_FEATURE_PAIRS))
    // Avalon❗✔️:    {
    // Avalon❗✔️:       for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    // Avalon❗✔️:       {
    for (atom_index, atom) in molecule.atoms.iter_mut().enumerate() {
        // Avalon❗✔️:          ap->color = StringToInt(periodic_table, ap->atom_symbol);
        atom.color = avalon_atomic_number(&atom.atom_symbol);
        // Avalon❗✔️:          if (0 == strcmp("H", ap->atom_symbol))
        // Avalon❗✔️:             ap->color = 0;      /* ignore hydrogens */
        // Avalon❗✔️:          else if (0 == strcmp("D", ap->atom_symbol))
        // Avalon❗✔️:             ap->color = 0;      /* ignore hydrogens */
        // Avalon❗✔️:          else if (0 == strcmp("T", ap->atom_symbol))
        // Avalon❗✔️:             ap->color = 0;      /* ignore hydrogens */
        if matches!(atom.atom_symbol.as_str(), "H" | "D" | "T") {
            atom.color = 0;
        // Avalon❗✔️:          else if (0 == strcmp("Q", ap->atom_symbol))
        // Avalon❗✔️:             ap->color = HETERO;
        } else if atom.atom_symbol == "Q" {
            atom.color = HETERO;
        // Avalon❗✔️:          else if (0 == strcmp("A", ap->atom_symbol))
        // Avalon❗✔️:             ap->color = GENERIC;
        // Avalon❗✔️:          else if (0 == strcmp("L", ap->atom_symbol))
        // Avalon❗✔️:             ap->color = GENERIC;
        } else if matches!(atom.atom_symbol.as_str(), "A" | "L") {
            atom.color = GENERIC;
        // Avalon❗✔️:          else if (0 == strcmp("C", ap->atom_symbol))
        // Avalon❗✔️:          {
        } else if atom.atom_symbol == "C" {
            // Avalon❗✔️:             ap->color = 6;      /* carbon second row elements are one class */
            atom.color = 6;
            // Avalon❗✔️:             if (cdegree[i] >= 3)
            // Avalon❗✔️:                ap->color = CSP3;
            if state.carbon_degree[atom_index] >= 3 {
                atom.color = CSP3;
            }
            // Avalon❗✔️:          }
            // Avalon❗✔️:          else if (ap->color > 1  &&  ap->color < 115)
            // Avalon❗✔️:                ap->color = HETERO;
        } else if atom.color > 1 && atom.color < 115 {
            atom.color = HETERO;
        // Avalon❗✔️:          else        /* This could be R atoms or other odd things */
        // Avalon❗✔️:             ap->color = 0;
        } else {
            atom.color = 0;
        }
        // Avalon❗✔️:          if (i+1 == exclude_atom) ap->color = 0;
        if atom_index as i32 + 1 == exclude_atom {
            atom.color = 0;
        }
        // Avalon❗✔️: // fprintf(stderr, "atom %d has color %d\n", i+1, ap->color);
        // Avalon❗✔️:       }
    }
    // Avalon❗✔️:       /*
    // Avalon❗✔️:        * Bond colors are already set OK, i.e. equal for A-H bonds
    // Avalon❗✔️:        */
    // Avalon❗✔️:       /* NOP */
}

fn count_class_spiders(
    molecule: &MoleculeState,
    state: &FingerprintPreprocessingState,
    counts: &mut [i32],
    exclude_atom: i32,
) -> i32 {
    let mut result = 0_i32;
    // Avalon❗✔️:       /* Now we start setting bits */
    // Avalon❗✔️:       for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    // Avalon❗✔️:       {
    for (atom_index, atom) in molecule.atoms.iter().enumerate() {
        // Avalon❗✔️:          if (i+1 == exclude_atom) continue;
        // Avalon❗✔️:          /* Spiders have at least three legs (;-) */
        // Avalon❗✔️:          if (degree[i] < 3) continue;
        // Avalon❗✔️:          /* Spider needs to be special atom or carbon */
        // Avalon❗✔️:          if (ap->color != CSP3  && ap->color != 6) continue;
        if atom_index as i32 + 1 == exclude_atom
            || state.degree[atom_index] < 3
            || (atom.color != CSP3 && atom.color != 6)
        {
            continue;
        }
        // Avalon❗✔️:          touched_indices[i] = 1; /* updating */
        // Avalon❗✔️:          for (j=0; j<=MAX_SPIDER; j++)
        // Avalon❗✔️:             hetero[j] = csp3[j] = 0;
        // Avalon❗✔️:          if (which_bits & USE_CLASS_SPIDERS)
        // Avalon❗✔️:          {
        // Avalon❗✔️:             strcpy(prefix, "CS:");
        // Avalon❗✔️:             SpecialNeighboursRec(mp, touched_indices,
        // Avalon❗✔️:                                  1, i, MAX_SPIDER,
        // Avalon❗✔️:                                  csp3, hetero, nbp,
        // Avalon❗✔️:                                  exclude_atom,
        // Avalon❗✔️: 				 prefix);
        // Avalon❗✔️:          }
        let (csp3, hetero) =
            collect_special_neighbours(molecule, &state.neighbours, atom_index, MAX_SPIDER, exclude_atom);
        // Avalon❗✔️:          touched_indices[i] = 0; /* down-dating */
        // Avalon❗✔️:          /* set bits for spiders with one CSP3 atom and two heteros */
        // Avalon❗✔️:          if (which_bits & USE_CLASS_SPIDERS)
        // Avalon❗✔️:          for (j=1; j<=MAX_SPIDER; j++)
        // Avalon❗✔️:          {
        for (csp3_distance, &csp3_count) in csp3.iter().enumerate().take(MAX_SPIDER + 1).skip(1) {
            // Avalon❗✔️:             strcpy(prefix, "CS:");
            // Avalon❗✔️:             if (csp3[j] == 0) continue;
            if csp3_count == 0 {
                continue;
            }
            // Avalon❗✔️:             seed = CLASS_SPIDER_SEED;
            let mut seed = CLASS_SPIDER_SEED;
            // Avalon❗✔️:             if (ap->color == HETERO)
            // Avalon❗✔️:             {
            // Avalon❗✔️: // seed = CLASS_SPIDER_SEED+HETERO*8+CSP3*11;
            // Avalon❗✔️:                seed = NEXT_SEED(seed, HETERO*8);
            // Avalon❗✔️:                seed = NEXT_SEED(seed, CSP3*11);
            // Avalon❗✔️:             }
            // Avalon❗✔️:             else
            // Avalon❗✔️:             {
            // Avalon❗✔️: // seed = CLASS_SPIDER_SEED+6*8+CSP3*11;
            // Avalon❗✔️:                seed = NEXT_SEED(seed, 6*8);
            // Avalon❗✔️:                seed = NEXT_SEED(seed, CSP3*11);
            // Avalon❗✔️:             }
            seed = next_hash(seed, (if atom.color == HETERO { HETERO } else { 6 } * 8) as u64);
            seed = next_hash(seed, (CSP3 * 11) as u64);
            // Avalon❗✔️:             for (j1=1; j1<=MAX_SPIDER; j1++)
            // Avalon❗✔️:             {
            for first_distance in 1..=MAX_SPIDER {
                // Avalon❗✔️:                tmp1 = hetero[j1];
                // Avalon❗✔️:                if (tmp1 <= 0) continue;
                if hetero[first_distance] <= 0 {
                    continue;
                }
                // Avalon❗✔️:                for (j2=j1; j2<=MAX_SPIDER; j2++)
                // Avalon❗✔️:                {
                for second_distance in first_distance..=MAX_SPIDER {
                    // Avalon❗✔️:                   tmp2 = hetero[j2];
                    let mut second_count = hetero[second_distance];
                    // Avalon❗✔️:                   if (j2 == j1) tmp2--; /* consumed in outer loop */
                    if second_distance == first_distance {
                        second_count -= 1;
                    }
                    // Avalon❗✔️:                   if (tmp2 <= 0) continue;
                    if second_count <= 0 {
                        continue;
                    }
                    // Avalon❗✔️:                   old_seed = seed;
                    // Avalon❗✔️:                   seed = NEXT_SEED(seed, j);
                    // Avalon❗✔️:                   seed = NEXT_SEED(seed, j1+j2);
                    let branch_seed = next_hash(
                        next_hash(seed, csp3_distance as u64),
                        (first_distance + second_distance) as u64,
                    );
                    // Avalon❗✔️:                   ADD_BIT(fp_counts, ncounts, seed);
                    add_bit(counts, branch_seed);
                    // Avalon❗✔️: 		  result++;
                    result += 1;
                    // Avalon❗✔️: // fprintf(stderr, " setting CSP3 bit for %d: %d|%d|%d\n", ap->color, j, j1, j2);
                    // Avalon❗✔️:                   seed = old_seed;
                    // The immutable parent seed remains unchanged.
                    // Avalon❗✔️:                }
                }
                // Avalon❗✔️:             }
            }
            // Avalon❗✔️:          }
        }

        // Avalon❗✔️:          /* don't hetero-spider normal carbons */
        // Avalon❗✔️:          if (ap->color != CSP3) continue;
        if atom.color != CSP3 {
            continue;
        }
        // Avalon❗✔️:          /* set bits for spiders with three defined HETERO atoms */
        // Avalon❗✔️:          if (which_bits & USE_CLASS_SPIDERS)
        // Avalon❗✔️:          for (j=1; j<=MAX_SPIDER; j++)
        // Avalon❗✔️:          {
        for first_distance in 1..=MAX_SPIDER {
            // Avalon❗✔️:             strcpy(prefix, "CS:");
            // Avalon❗✔️:             if (hetero[j] == 0) continue;
            if hetero[first_distance] == 0 {
                continue;
            }
            // Avalon❗✔️:             seed = CLASS_SPIDER_SEED;
            let mut seed = CLASS_SPIDER_SEED;
            // Avalon❗✔️:             if (ap->color == HETERO)
            // Avalon❗✔️:             {
            // Avalon❗✔️: // seed = CLASS_SPIDER_SEED+HETERO*8+HETERO*11;
            // Avalon❗✔️:                seed = NEXT_SEED(seed, HETERO*8);
            // Avalon❗✔️:                seed = NEXT_SEED(seed, HETERO*11);
            // Avalon❗✔️:             }
            // Avalon❗✔️:             else
            // Avalon❗✔️:             {
            // Avalon❗✔️: // seed = CLASS_SPIDER_SEED+6*8+HETERO*11;
            // Avalon❗✔️:                seed = NEXT_SEED(seed, 6*8);
            // Avalon❗✔️:                seed = NEXT_SEED(seed, HETERO*11);
            // Avalon❗✔️:             }
            seed = next_hash(seed, (if atom.color == HETERO { HETERO } else { 6 } * 8) as u64);
            seed = next_hash(seed, (HETERO * 11) as u64);
            // Avalon❗✔️:             for (j1=j; j1<=MAX_SPIDER; j1++)
            // Avalon❗✔️:             {
            for second_distance in first_distance..=MAX_SPIDER {
                // Avalon❗✔️:                tmp1 = hetero[j1];
                let mut second_count = hetero[second_distance];
                // Avalon❗✔️:                if (j1 == j) tmp1--; /* we've consumed this one in outer loop */
                if second_distance == first_distance {
                    second_count -= 1;
                }
                // Avalon❗✔️:                if (tmp1 <= 0) continue;
                if second_count <= 0 {
                    continue;
                }
                // Avalon❗✔️:                for (j2=j1; j2<=MAX_SPIDER; j2++)
                // Avalon❗✔️:                {
                for third_distance in second_distance..=MAX_SPIDER {
                    // Avalon❗✔️:                   tmp2 = hetero[j2];
                    let mut third_count = hetero[third_distance];
                    // Avalon❗✔️:                   if (j2 == j)  tmp2--; /* consumed in outer loop */
                    if third_distance == first_distance {
                        third_count -= 1;
                    }
                    // Avalon❗✔️:                   if (j2 == j1) tmp2--; /* consumed in outer loop */
                    if third_distance == second_distance {
                        third_count -= 1;
                    }
                    // Avalon❗✔️:                   if (tmp2 <= 0) continue;
                    if third_count <= 0 {
                        continue;
                    }
                    // Avalon❗✔️:                   old_seed = seed;
                    // Avalon❗✔️:                   seed = NEXT_SEED(seed, j*j1*j2);
                    let branch_seed = next_hash(seed, (first_distance * second_distance * third_distance) as u64);
                    // Avalon❗✔️:                   ADD_BIT(fp_counts, ncounts, seed);
                    add_bit(counts, branch_seed);
                    // Avalon❗✔️: 		  result++;
                    result += 1;
                    // Avalon❗✔️: // fprintf(stderr, " setting HETERO bit for %d: %d|%d|%d\n", ap->color, j, j1, j2);
                    // Avalon❗✔️:                   // Additional bits for quarternary centers
                    // Avalon❗✔️:                   if (degree[i] > 3) {ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, 501)); result++;}
                    if state.degree[atom_index] > 3 {
                        add_bit(counts, next_hash(branch_seed, 501));
                        result += 1;
                    }
                    // Avalon❗✔️:                   seed = old_seed;
                    // The immutable parent seed remains unchanged.
                    // Avalon❗✔️:                }
                }
                // Avalon❗✔️:             }
            }
            // Avalon❗✔️:          }
        }
        // Avalon❗✔️:       }
    }
    result
}

fn prepare_feature_pair_state(molecule: &mut MoleculeState, state: &FingerprintPreprocessingState, exclude_atom: i32) {
    // Avalon❗✔️:       /**
    // Avalon❗✔️:        * Collect bits that represent feature/path_length/feature triples.
    // Avalon❗✔️:        */
    // Avalon❗✔️:       if (which_bits & USE_FEATURE_PAIRS)
    // Avalon❗✔️:       {
    // Avalon❗✔️:          strcpy(prefix, "FP:");
    // Avalon❗✔️:          /* set feature flags in atom colors */
    // Avalon❗✔️:          for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    // Avalon❗✔️:          {
    for (atom_index, atom) in molecule.atoms.iter_mut().enumerate() {
        let was_hetero = atom.color == HETERO;
        // Avalon❗✔️:             if (0 == strcmp(ap->atom_symbol, "C"))      flags = C_FLAG;
        // Avalon❗✔️:             else if (0 == strcmp(ap->atom_symbol, "O")) flags = O_FLAG;
        // Avalon❗✔️:             else if (0 == strcmp(ap->atom_symbol, "N")) flags = N_FLAG;
        // Avalon❗✔️:             else if (0 == strcmp(ap->atom_symbol, "S")) flags = S_FLAG;
        // Avalon❗✔️:             else if (0 == strcmp(ap->atom_symbol, "P")) flags = P_FLAG;
        // Avalon❗✔️:             else if (AtomSymbolMatch(ap->atom_symbol,"F,Cl,Br,I,At"))
        // Avalon❗✔️:                flags = X_FLAG;
        // Avalon❗✔️:             else
        // Avalon❗✔️:                flags = 0;
        let mut flags = match atom.atom_symbol.as_str() {
            "C" => C_FLAG,
            "O" => O_FLAG,
            "N" => N_FLAG,
            "S" => S_FLAG,
            "P" => P_FLAG,
            symbol if atom_symbol_match(symbol, "F,Cl,Br,I,At") => X_FLAG,
            _ => 0,
        };
        // Avalon❗✔️:             if (ap->color == HETERO)
        // Avalon❗✔️: 	    {
        // Avalon❗✔️: 	       flags |= HETERO_FLAG;
        // Avalon❗✔️: 	    }
        if was_hetero {
            flags |= HETERO_FLAG;
        }
        // Avalon❗✔️: // fprintf(stderr, "HETERO_FLAG set for atom %d\n", i+1);
        // Avalon❗✔️:             if (cdegree[i] >= 3)     flags |= CSP3_FLAG;
        if state.carbon_degree[atom_index] >= 3 {
            flags |= CSP3_FLAG;
        }
        // Avalon❗✔️: // fprintf(stderr, "%s atom %d has degree %d\n", ap->atom_symbol, i+1, degree[i]);
        // Avalon❗✔️:             if (degree[i] >= 4)
        // Avalon❗✔️:             {
        // Avalon❗✔️:                flags |= QUART_FLAG;
        // Avalon❗✔️:             }
        if state.degree[atom_index] >= 4 {
            flags |= QUART_FLAG;
        }
        // Avalon❗✔️:             if (atom_status[i] > 0  &&  degree[i] >= 3)
        // Avalon❗✔️:             {
        if state.atom_status[atom_index] > 0 && state.degree[atom_index] >= 3 {
            // Avalon❗✔️:                flags |= RING_SUBST_FLAG;
            flags |= RING_SUBST_FLAG;
            // Avalon❗✔️:                if (0 != (ap->rsize_flags&SPECIAL_RING))
            // Avalon❗✔️:                {
            // Avalon❗✔️:                    flags |= RS_SPECIAL_FLAG;
            // Avalon❗✔️:                }
            if atom.rsize_flags & SPECIAL_RING != 0 {
                flags |= RS_SPECIAL_FLAG;
            }
            // Avalon❗✔️: // fprintf(stderr, "RS_SPECIAL_FLAG set for atom %d\n", i+1);
            // Avalon❗✔️:             }
        }
        // Avalon❗✔️:             if (i+1 == exclude_atom) flags = 0;
        if atom_index as i32 + 1 == exclude_atom {
            flags = 0;
        }
        // Avalon❗✔️:             ap->color = flags;
        atom.color = flags;
        // Avalon❗✔️:          }
    }
}

fn count_feature_pairs(
    molecule: &MoleculeState,
    counts: &mut [i32],
    length_matrix: &[Vec<i32>],
    exclude_atom: i32,
) -> i32 {
    let mut result = 0_i32;
    // Avalon❗✔️:          /* collect bits for selected feature pairs */
    // Avalon❗✔️:          if(0)         // Class disabled to save bit density
    // Avalon❗✔️:          result +=
    // Avalon❗✔️:             SetFeatureBits(mp, fp_counts, ncounts,
    // Avalon❗✔️:                                  CSP3_FLAG,
    // Avalon❗✔️:                                  HETERO_FLAG,         /* to hetero or ring subst */
    // Avalon❗✔️:                                  2, 3,               /* with path length 1 to 9 */
    // Avalon❗✔️:                                  FALSE,               /* don't use count */
    // Avalon❗✔️:                                  TRUE,                /* use atom type flags */
    // Avalon❗✔️:                                  length_matrix,
    // Avalon❗✔️:                                  1237,                /* seed = 1237 */
    // Avalon❗✔️: 			   exclude_atom,
    // Avalon❗✔️: 			   prefix);
    // Avalon❗✔️:          if(0)         // Class disabled to save bit density
    // Avalon❗✔️:          result +=
    // Avalon❗✔️:             SetFeatureBits(mp, fp_counts, ncounts,
    // Avalon❗✔️:                                  HETERO_FLAG,         /* from ring substitution */
    // Avalon❗✔️:                                  HETERO_FLAG,         /* to hetero or ring subst */
    // Avalon❗✔️:                                  1, 12,               /* with path length 1 to 10 */
    // Avalon❗✔️:                                  TRUE,                /* use count */
    // Avalon❗✔️:                                  TRUE,                /* use atom type flags */
    // Avalon❗✔️:                                  length_matrix,
    // Avalon❗✔️:                                  1237,                /* seed = 1237 */
    // Avalon❗✔️: 			   exclude_atom,
    // Avalon❗✔️: 			   prefix);
    // The two source branches above are compile-time disabled.
    // Avalon❗✔️:          if (1)
    // Avalon❗✔️:          result +=
    // Avalon❗✔️:             SetFeatureBits(mp, fp_counts, ncounts,
    // Avalon❗✔️:                                  RING_SUBST_FLAG,     /* from ring substitution */
    // Avalon❗✔️:                                  RING_SUBST_FLAG,     /* to ring substitution */
    // Avalon❗✔️:                                  5, 7,               /* with path length 1 to 12 */
    // Avalon❗✔️:                                  FALSE,               /* don't use count */
    // Avalon❗✔️:                                  TRUE,                /* use atom type flags */
    // Avalon❗✔️:                                  length_matrix,
    // Avalon❗✔️:                                  2237,                /* seed = 2237 */
    // Avalon❗✔️: 			   exclude_atom,
    // Avalon❗✔️: 			   prefix);
    result += set_feature_bits(
        molecule,
        counts,
        RING_SUBST_FLAG,
        RING_SUBST_FLAG,
        5,
        7,
        false,
        true,
        length_matrix,
        2237,
        exclude_atom,
    );
    // Avalon❗✔️:          if (0)         // Class disabled to save bit density
    // Avalon❗✔️:          result +=
    // Avalon❗✔️:             SetRingSizePairBits(mp, fp_counts, ncounts,
    // Avalon❗✔️: 			    RING_SUBST_FLAG,     /* from ring substitution */
    // Avalon❗✔️:                                  HETERO_FLAG,         /* to hetero atom */
    // Avalon❗✔️:                                  2, 5,               /* with path length 1 to 12 */
    // Avalon❗✔️:                                  length_matrix,
    // Avalon❗✔️:                                  3237,                /* seed = 2237 */
    // Avalon❗✔️: 			    exclude_atom,
    // Avalon❗✔️: 			    prefix);
    // Avalon❗✔️:          if (0)         // Class disabled to save bit density
    // Avalon❗✔️:          result +=
    // Avalon❗✔️:             SetFeatureBits(mp, fp_counts, ncounts,
    // Avalon❗✔️: 			   RS_SPECIAL_FLAG,/* from special ring substitution */
    // Avalon❗✔️:                                  HETERO_FLAG,    /* to hetero atom */
    // Avalon❗✔️:                                  2, 4,                /* with path length 1 to 5 */
    // Avalon❗✔️:                                  TRUE,               /* don't use count */
    // Avalon❗✔️:                                  FALSE,                /* use atom type flags */
    // Avalon❗✔️:                                  length_matrix,
    // Avalon❗✔️:                                  3237,                /* seed = 3237 */
    // Avalon❗✔️: 			   exclude_atom,
    // Avalon❗✔️: 			   prefix);
    // The two source branches above are compile-time disabled.
    // Avalon❗✔️:          if (1)
    // Avalon❗✔️:          result +=
    // Avalon❗✔️:             SetFeatureBits(mp, fp_counts, ncounts,
    // Avalon❗✔️:                                  QUART_FLAG,          /* from quartenary atom */
    // Avalon❗✔️:                                  HETERO_FLAG,         /* to hetero or ring subst */
    // Avalon❗✔️:                                  1, 8,                /* with path length 1 to 8 */
    // Avalon❗✔️:                                  FALSE,               /* don't use count */
    // Avalon❗✔️:                                  TRUE,                /* use atom type flags */
    // Avalon❗✔️:                                  length_matrix,
    // Avalon❗✔️:                                  4237,                /* seed = 4237 */
    // Avalon❗✔️: 			   exclude_atom,
    // Avalon❗✔️: 			   prefix);
    result += set_feature_bits(
        molecule,
        counts,
        QUART_FLAG,
        HETERO_FLAG,
        1,
        8,
        false,
        true,
        length_matrix,
        4237,
        exclude_atom,
    );
    // Avalon❗✔️:          if(1)
    // Avalon❗✔️:          result +=
    // Avalon❗✔️:             SetFeatureBits(mp, fp_counts, ncounts,
    // Avalon❗✔️:                                  QUART_FLAG,          /* from quartenary atom */
    // Avalon❗✔️:                                  RING_SUBST_FLAG,
    // Avalon❗✔️:                                  1, 6,                /* with path length 1 to 8 */
    // Avalon❗✔️:                                  FALSE,               /* don't use count */
    // Avalon❗✔️:                                  TRUE,                /* use atom type flags */
    // Avalon❗✔️:                                  length_matrix,
    // Avalon❗✔️:                                  5237,                /* seed = 4237 */
    // Avalon❗✔️: 			   exclude_atom,
    // Avalon❗✔️: 			   prefix);
    result += set_feature_bits(
        molecule,
        counts,
        QUART_FLAG,
        RING_SUBST_FLAG,
        1,
        6,
        false,
        true,
        length_matrix,
        5237,
        exclude_atom,
    );
    // Avalon❗✔️:          if(1)
    // Avalon❗✔️:          result +=
    // Avalon❗✔️:             SetFeatureBits(mp, fp_counts, ncounts,
    // Avalon❗✔️:                                  X_FLAG,          /* from halogen atom */
    // Avalon❗✔️:                                  CSP3_FLAG,
    // Avalon❗✔️:                                  1, 1,                /* with path length 1 to 1 */
    // Avalon❗✔️:                                  FALSE,               /* don't use count */
    // Avalon❗✔️:                                  TRUE,                /* use atom type flags */
    // Avalon❗✔️:                                  length_matrix,
    // Avalon❗✔️:                                  15237,                /* seed = 15237 */
    // Avalon❗✔️: 			   exclude_atom,
    // Avalon❗✔️: 			   prefix);
    result += set_feature_bits(
        molecule,
        counts,
        X_FLAG,
        CSP3_FLAG,
        1,
        1,
        false,
        true,
        length_matrix,
        15237,
        exclude_atom,
    );
    // Avalon❗✔️:          if (0)         // Class disabled to save bit density
    // Avalon❗✔️:          result +=
    // Avalon❗✔️:             SetFeatureBits(mp, fp_counts, ncounts,
    // Avalon❗✔️: 			   RS_SPECIAL_FLAG,/* from special ring substitution */
    // Avalon❗✔️:                                  RING_SUBST_FLAG,    /* to hetero atom */
    // Avalon❗✔️:                                  2, 4,                /* with path length 1 to 5 */
    // Avalon❗✔️:                                  TRUE,               /* don't use count */
    // Avalon❗✔️:                                  FALSE,                /* use atom type flags */
    // Avalon❗✔️:                                  length_matrix,
    // Avalon❗✔️:                                  6237,                /* seed = 6237 */
    // Avalon❗✔️: 			   exclude_atom,
    // Avalon❗✔️: 			   prefix);
    // Avalon❗✔️:          if(0)         // Class disabled to save bit density
    // Avalon❗✔️:          result +=
    // Avalon❗✔️:             SetFeatureBits(mp, fp_counts, ncounts,
    // Avalon❗✔️:                                  HETERO_FLAG,         /* from hetero */
    // Avalon❗✔️:                                  RING_SUBST_FLAG,     /* to ring substitution */
    // Avalon❗✔️:                                  1, 6,                /* with path length 1 to 8 */
    // Avalon❗✔️:                                  FALSE,               /* don't use count */
    // Avalon❗✔️:                                  TRUE,                /* use atom type flags */
    // Avalon❗✔️:                                  length_matrix,
    // Avalon❗✔️:                                  7237,                /* seed = 4237 */
    // Avalon❗✔️: 			   exclude_atom,
    // Avalon❗✔️: 			   prefix);
    // Avalon❗✔️:          /* Set bits for ring-subst/ring-subst/hetero triples */
    // Avalon❗✔️:          if (0) // too many spurious bits
    // Avalon❗✔️:          {
    // Avalon❗✔️:             for (i1=0, ap1=mp->atom_array; i1<mp->n_atoms; i1++, ap1++)
    // Avalon❗✔️:             {
    // Avalon❗✔️:                if (i1+1 == exclude_atom) continue;
    // Avalon❗✔️:                if (0 == (ap1->color&RING_SUBST_FLAG)) continue;
    // Avalon❗✔️:                /* first atom must be in a non-sixmembered ring */
    // Avalon❗✔️:                if (!(ap1->rsize_flags&SPECIAL_RING)) continue;
    // Avalon❗✔️:                for (i2=0, ap2=mp->atom_array; i2<mp->n_atoms; i2++, ap2++)
    // Avalon❗✔️:                {
    // Avalon❗✔️:                   if (i1 == i2) continue;
    // Avalon❗✔️:                   if (i2+1 == exclude_atom) continue;
    // Avalon❗✔️:                   if (0 == (ap2->color&RING_SUBST_FLAG)) continue;
    // Avalon❗✔️:                   for (i3=0, ap3=mp->atom_array; i3<mp->n_atoms; i3++, ap3++)
    // Avalon❗✔️:                   {
    // Avalon❗✔️:                      if (i1 == i3) continue;
    // Avalon❗✔️:                      if (i2 == i3) continue;
    // Avalon❗✔️:                      if (i3+1 == exclude_atom) continue;
    // Avalon❗✔️:                      if (0 == (ap3->color&HETERO_FLAG)) continue;
    // Avalon❗✔️:                      for (j=2; j<=6; j++)
    // Avalon❗✔️:                      {
    // Avalon❗✔️:                         if (0 == (length_matrix[i1][i2]&(1<<j))) continue;
    // Avalon❗✔️:                         for (j1=2; j1<=5; j1++)
    // Avalon❗✔️:                         {
    // Avalon❗✔️:                            if (0 == (length_matrix[i2][i3]&(1<<j1))) continue;
    // Avalon❗✔️:                            for (j2=2; j2<=5; j2++)
    // Avalon❗✔️:                            {
    // Avalon❗✔️:                               if (0 == (length_matrix[i3][i1]&(1<<j2)))
    // Avalon❗✔️:                                  continue;
    // Avalon❗✔️:                               // Make sure we have a real triangle
    // Avalon❗✔️:                               if (j+j1  == j2) continue;
    // Avalon❗✔️:                               if (j+j2  == j1) continue;
    // Avalon❗✔️:                               if (j1+j1 == j)  continue;
    // Avalon❗✔️:                               if (j+j1+j2 > 11) continue;  // spider too large
    // Avalon❗✔️:                               if (j+j1+j2 <  8) continue;  // spider too small
    // Avalon❗✔️:                               seed = CLASS_SPIDER_SEED;
    // Avalon❗✔️:                               seed = NEXT_SEED(seed, j+j1+j2);
    // Avalon❗✔️:                               seed = NEXT_SEED(seed, j*j1*j2);
    // Avalon❗✔️:                               // distinguish ring sizes of second ring-subst
    // Avalon❗✔️:                               // for (k=3; k<15; k++)
    // Avalon❗✔️:                               for (k=3; k<9; k++)
    // Avalon❗✔️:                               {
    // Avalon❗✔️:                                  if (!(ap2->rsize_flags & (1<<k))) continue;
    // Avalon❗✔️:                                  ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, k*213));
    // Avalon❗✔️: 				 result++;
    // Avalon❗✔️:                               }
    // Avalon❗✔️: // fprintf(stderr, "processing tripple (%d|%x)-%d-(%d|%x)-%d-(%d|%s)-%d\n",
    // Avalon❗✔️: // i1+1, ap1->rsize_flags, j,
    // Avalon❗✔️: // i2+1, ap2->rsize_flags, j1,
    // Avalon❗✔️: // i3+1, ap3->atom_symbol, j2);
    // Avalon❗✔️:                            }
    // Avalon❗✔️:                         }
    // Avalon❗✔️:                      }
    // Avalon❗✔️:                   }
    // Avalon❗✔️:                }
    // Avalon❗✔️:             }
    // The literal `if (0)` makes this complete source traversal unreachable.
    // Avalon❗✔️:          }
    // Avalon❗✔️:       }
    // Avalon❗✔️:    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::properties::avalon_fingerprint::fingerprint_state::USE_DY_AROMATICITY;
    use crate::properties::avalon_fingerprint::reaccs::{Atom, Bond};

    const HIGH_MASKS: &[u32] = &[
        0x0800, 0x1000, 0x2000, 0x4000, 0x1800, 0x2800, 0x4800, 0x3000, 0x5000, 0x6000,
    ];

    const NATIVE_GOLDENS: &str = r#"high_matrix=spiro_mask_0800_q0_dy0 result=4 counts=28:3,38:1,42:1,55:3, bits=0000001040048000
high_matrix=spiro_mask_1000_q0_dy0 result=37 counts=0:1,2:3,3:3,4:4,18:5,25:3,31:2,39:1,40:1,47:2,48:1,51:1,52:5,53:1,54:1,58:2,63:1, bits=1d00048280817984
high_matrix=spiro_mask_2000_q0_dy0 result=4 counts=7:1,9:1,23:1,35:1, bits=8002800008000000
high_matrix=spiro_mask_4000_q0_dy0 result=20 counts=11:1,14:1,30:1,31:1,32:1,41:1,53:1,55:1,60:1, bits=004800c00102a010
high_matrix=spiro_mask_1800_q0_dy0 result=41 counts=0:1,2:3,3:3,4:4,18:5,25:3,28:3,31:2,38:1,39:1,40:1,42:1,47:2,48:1,51:1,52:5,53:1,54:1,55:3,58:2,63:1, bits=1d000492c085f984
high_matrix=spiro_mask_2800_q0_dy0 result=8 counts=7:1,9:1,23:1,28:3,35:1,38:1,42:1,55:3, bits=8002801048048000
high_matrix=spiro_mask_4800_q0_dy0 result=24 counts=11:1,14:1,28:3,30:1,31:1,32:1,38:1,41:1,42:1,53:1,55:4,60:1, bits=004800d04106a010
high_matrix=spiro_mask_3000_q0_dy0 result=41 counts=0:1,2:3,3:3,4:4,7:1,9:1,18:5,23:1,25:3,31:2,35:1,39:1,40:1,47:2,48:1,51:1,52:5,53:1,54:1,58:2,63:1, bits=9d02848288817984
high_matrix=spiro_mask_5000_q0_dy0 result=57 counts=0:1,2:3,3:3,4:4,11:1,14:1,18:5,25:3,30:1,31:3,32:1,39:1,40:1,41:1,47:2,48:1,51:1,52:5,53:2,54:1,55:1,58:2,60:1,63:1, bits=1d4804c28183f994
high_matrix=spiro_mask_6000_q0_dy0 result=24 counts=7:1,9:1,11:1,14:1,23:1,30:1,31:1,32:1,35:1,41:1,53:1,55:1,60:1, bits=804a80c00902a010
high_matrix=spiro_mask_0800_q1_dy0 result=4 counts=28:3,38:1,42:1,55:3, bits=0000001040048000
high_matrix=spiro_mask_1000_q1_dy0 result=0 counts= bits=0000000000000000
high_matrix=spiro_mask_2000_q1_dy0 result=4 counts=7:1,9:1,23:1,35:1, bits=8002800008000000
high_matrix=spiro_mask_4000_q1_dy0 result=2 counts=41:1, bits=0000000000020000
high_matrix=spiro_mask_1800_q1_dy0 result=4 counts=28:3,38:1,42:1,55:3, bits=0000001040048000
high_matrix=spiro_mask_2800_q1_dy0 result=8 counts=7:1,9:1,23:1,28:3,35:1,38:1,42:1,55:3, bits=8002801048048000
high_matrix=spiro_mask_4800_q1_dy0 result=6 counts=28:3,38:1,41:1,42:1,55:3, bits=0000001040068000
high_matrix=spiro_mask_3000_q1_dy0 result=4 counts=7:1,9:1,23:1,35:1, bits=8002800008000000
high_matrix=spiro_mask_5000_q1_dy0 result=2 counts=41:1, bits=0000000000020000
high_matrix=spiro_mask_6000_q1_dy0 result=6 counts=7:1,9:1,23:1,35:1,41:1, bits=8002800008020000
high_matrix=spiro_mask_0800_q0_dy1 result=4 counts=28:3,38:1,42:1,55:3, bits=0000001040048000
high_matrix=spiro_mask_1000_q0_dy1 result=37 counts=0:1,2:3,3:3,4:4,18:5,25:3,31:2,39:1,40:1,47:2,48:1,51:1,52:5,53:1,54:1,58:2,63:1, bits=1d00048280817984
high_matrix=spiro_mask_2000_q0_dy1 result=4 counts=7:1,9:1,23:1,35:1, bits=8002800008000000
high_matrix=spiro_mask_4000_q0_dy1 result=20 counts=11:1,14:1,30:1,31:1,32:1,41:1,53:1,55:1,60:1, bits=004800c00102a010
high_matrix=spiro_mask_1800_q0_dy1 result=41 counts=0:1,2:3,3:3,4:4,18:5,25:3,28:3,31:2,38:1,39:1,40:1,42:1,47:2,48:1,51:1,52:5,53:1,54:1,55:3,58:2,63:1, bits=1d000492c085f984
high_matrix=spiro_mask_2800_q0_dy1 result=8 counts=7:1,9:1,23:1,28:3,35:1,38:1,42:1,55:3, bits=8002801048048000
high_matrix=spiro_mask_4800_q0_dy1 result=24 counts=11:1,14:1,28:3,30:1,31:1,32:1,38:1,41:1,42:1,53:1,55:4,60:1, bits=004800d04106a010
high_matrix=spiro_mask_3000_q0_dy1 result=41 counts=0:1,2:3,3:3,4:4,7:1,9:1,18:5,23:1,25:3,31:2,35:1,39:1,40:1,47:2,48:1,51:1,52:5,53:1,54:1,58:2,63:1, bits=9d02848288817984
high_matrix=spiro_mask_5000_q0_dy1 result=57 counts=0:1,2:3,3:3,4:4,11:1,14:1,18:5,25:3,30:1,31:3,32:1,39:1,40:1,41:1,47:2,48:1,51:1,52:5,53:2,54:1,55:1,58:2,60:1,63:1, bits=1d4804c28183f994
high_matrix=spiro_mask_6000_q0_dy1 result=24 counts=7:1,9:1,11:1,14:1,23:1,30:1,31:1,32:1,35:1,41:1,53:1,55:1,60:1, bits=804a80c00902a010
high_matrix=spiro_mask_0800_q1_dy1 result=4 counts=28:3,38:1,42:1,55:3, bits=0000001040048000
high_matrix=spiro_mask_1000_q1_dy1 result=0 counts= bits=0000000000000000
high_matrix=spiro_mask_2000_q1_dy1 result=4 counts=7:1,9:1,23:1,35:1, bits=8002800008000000
high_matrix=spiro_mask_4000_q1_dy1 result=2 counts=41:1, bits=0000000000020000
high_matrix=spiro_mask_1800_q1_dy1 result=4 counts=28:3,38:1,42:1,55:3, bits=0000001040048000
high_matrix=spiro_mask_2800_q1_dy1 result=8 counts=7:1,9:1,23:1,28:3,35:1,38:1,42:1,55:3, bits=8002801048048000
high_matrix=spiro_mask_4800_q1_dy1 result=6 counts=28:3,38:1,41:1,42:1,55:3, bits=0000001040068000
high_matrix=spiro_mask_3000_q1_dy1 result=4 counts=7:1,9:1,23:1,35:1, bits=8002800008000000
high_matrix=spiro_mask_5000_q1_dy1 result=2 counts=41:1, bits=0000000000020000
high_matrix=spiro_mask_6000_q1_dy1 result=6 counts=7:1,9:1,23:1,35:1,41:1, bits=8002800008020000
high_boundary=fused_mask_0800_q0_dy0 result=3 counts=32:7,41:1,59:7, bits=0000000001020008
high_boundary=spider6_7_mask_2000_q0_dy0 result=2 counts=6:1,60:1, bits=4000000000000010
high_boundary=feature8_mask_4000_q0_dy0 result=2 counts=56:1, bits=0000000000000001
high_boundary=ring9_10_mask_0800_q0_dy0 result=2 counts=10:9,31:9, bits=0004008000000000"#;

    fn molecule(symbols: &[&str], endpoints: &[[i32; 2]]) -> MoleculeState {
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
                .copied()
                .map(|atoms| Bond {
                    atoms,
                    bond_type: 1,
                    ..Bond::default()
                })
                .collect(),
            ..MoleculeState::default()
        }
    }

    fn high_matrix_fixture() -> MoleculeState {
        molecule(
            &[
                "C", "C", "C", "C", "Cl", "N", "O", "S", "C", "C", "C", "C", "C", "C", "C", "C", "C",
            ],
            &[
                [1, 2],
                [1, 3],
                [1, 4],
                [1, 5],
                [2, 6],
                [3, 7],
                [4, 8],
                [9, 10],
                [10, 11],
                [11, 9],
                [9, 12],
                [12, 13],
                [13, 14],
                [14, 15],
                [15, 9],
                [10, 16],
                [13, 17],
            ],
        )
    }

    fn fused_fixture() -> MoleculeState {
        molecule(
            &["C", "C", "C", "C", "C", "C"],
            &[[1, 2], [2, 3], [3, 4], [4, 1], [3, 5], [5, 6], [6, 4]],
        )
    }

    fn spider_boundary_fixture() -> MoleculeState {
        let mut symbols = vec!["C"];
        let mut endpoints = Vec::new();
        let mut atom = 1_i32;
        for branch in 0..4 {
            let mut previous = 1_i32;
            let branch_length = if branch == 3 { 7 } else { 6 };
            for depth in 1..=branch_length {
                atom += 1;
                symbols.push(if depth == branch_length {
                    match branch {
                        0 => "N",
                        1 => "O",
                        _ => "S",
                    }
                } else {
                    "C"
                });
                endpoints.push([previous, atom]);
                previous = atom;
            }
        }
        molecule(&symbols, &endpoints)
    }

    fn feature_boundary_fixture() -> MoleculeState {
        let mut symbols = vec!["C"; 12];
        symbols[8] = "N";
        let mut endpoints = (1_i32..=8).map(|atom| [atom, atom + 1]).collect::<Vec<_>>();
        endpoints.extend([[1, 10], [1, 11], [1, 12]]);
        molecule(&symbols, &endpoints)
    }

    fn ring_boundary_fixture() -> MoleculeState {
        let mut endpoints = (1_i32..9).map(|atom| [atom, atom + 1]).collect::<Vec<_>>();
        endpoints.push([9, 1]);
        endpoints.extend((10_i32..19).map(|atom| [atom, atom + 1]));
        endpoints.push([19, 10]);
        molecule(&vec!["C"; 19], &endpoints)
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
        bytes.into_iter().map(|byte| format!("{byte:02x}")).collect()
    }

    fn run_case(mut molecule: MoleculeState, mask: u32, as_query: bool, daylight: bool) -> (i32, String, String) {
        let mut counts = vec![0_i32; 64];
        let result = count_high_flag_families(
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
        NATIVE_GOLDENS.lines().filter(move |line| line.starts_with(prefix))
    }

    #[test]
    fn high_single_and_pairwise_matrix_matches_native_counts_and_bits() {
        assert_eq!(HIGH_MASKS.len(), 10);
        assert_eq!(HIGH_MASKS.iter().filter(|mask| mask.count_ones() == 1).count(), 4);
        assert_eq!(HIGH_MASKS.iter().filter(|mask| mask.count_ones() == 2).count(), 6);
        let mut expected = expected_lines("high_matrix=");
        for daylight in [false, true] {
            for as_query in [false, true] {
                for &mask in HIGH_MASKS {
                    let (result, counts, bits) = run_case(high_matrix_fixture(), mask, as_query, daylight);
                    let actual = format!(
                        "high_matrix=spiro_mask_{mask:04x}_q{}_dy{} result={result} counts={counts} bits={bits}",
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
    fn fused_and_source_size_boundaries_match_native_counts_and_bits() {
        let cases = [
            ("fused", 0x0800, fused_fixture()),
            ("spider6_7", 0x2000, spider_boundary_fixture()),
            ("feature8", 0x4000, feature_boundary_fixture()),
            ("ring9_10", 0x0800, ring_boundary_fixture()),
        ];
        let mut expected = expected_lines("high_boundary=");
        for (fixture, mask, molecule) in cases {
            let (result, counts, bits) = run_case(molecule, mask, false, false);
            let actual =
                format!("high_boundary={fixture}_mask_{mask:04x}_q0_dy0 result={result} counts={counts} bits={bits}");
            assert_eq!(Some(actual.as_str()), expected.next(), "{actual}");
        }
        assert_eq!(expected.next(), None);
    }
}
