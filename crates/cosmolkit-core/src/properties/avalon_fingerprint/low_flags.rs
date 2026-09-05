//! Source-order implementation of Avalon fingerprint flags `0x000001..=0x000010`.

use crate::FingerprintError;

use super::AvalonFingerprintFlags;
use super::fingerprint_state::{FingerprintPreprocessingState, avalon_atomic_number, with_prepared_fingerprint_state};
use super::hash::next_hash;
use super::reaccs::MoleculeState;
use super::traversal::{
    ATOM_CLASS_PATH_SEED, ATOM_COUNT_SEED, ATOM_SYMBOL_PATH_SEED, PathFlags, RING_PATH_SEED, RING_PATTERN_SEED,
    add_bit, add_bit_count, set_path_bits_rec,
};

const NCOUNT_HASH: usize = 128;
const NCOUNT_SEED_HASH: usize = 128 * 128;
const SUB_AS_IS: i32 = -2;
const SUB_ONE: i32 = 1;
const SUB_MORE: i32 = 6;
const ANY_BOND: i32 = 8;
const SPECIAL_RING: u32 = 0xfc & !(1 << 6);
const ACCUMULATE_BITS: i32 = 0x0002;

pub(super) fn count_low_flag_families(
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
            Ok(count_low_flag_families_prepared(
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

pub(super) fn count_low_flag_families_prepared(
    working: &mut MoleculeState,
    state: &FingerprintPreprocessingState,
    counts: &mut [i32],
    which_bits: AvalonFingerprintFlags,
    as_query: bool,
    exclude_atom: i32,
) -> i32 {
    let mut result = 0_i32;
    if which_bits.contains(AvalonFingerprintFlags::ATOM_COUNT) {
        result += count_atom_types(working, state, counts, exclude_atom);
    }
    if which_bits.contains(AvalonFingerprintFlags::ATOM_SYMBOL_PATH) {
        result += count_atom_symbol_paths(working, state, counts, as_query, exclude_atom);
    }
    prepare_ring_path_state(working, state, exclude_atom);
    if which_bits.contains(AvalonFingerprintFlags::RING_PATH) {
        result += count_ring_paths(working, state, counts, exclude_atom);
    }
    prepare_bond_path_state(working, exclude_atom);
    prepare_atom_class_state(working, exclude_atom);
    if which_bits.contains(AvalonFingerprintFlags::ATOM_CLASS_PATH) {
        result += count_atom_class_paths_first_pass(working, state, counts, exclude_atom);
    }
    unify_bond_classes(working, exclude_atom);
    if which_bits.contains(AvalonFingerprintFlags::ATOM_CLASS_PATH) {
        result += count_atom_class_paths_second_pass(working, state, counts, exclude_atom);
    }
    prepare_ring_pattern_state(working, state, exclude_atom);
    if which_bits.contains(AvalonFingerprintFlags::RING_PATTERN) {
        result += count_ring_patterns(working, state, counts, exclude_atom);
    }
    result
}

pub(super) fn prepare_ring_path_state(
    molecule: &mut MoleculeState,
    state: &FingerprintPreprocessingState,
    exclude_atom: i32,
) {
    // Avalon❗✔️:    /* Compute ring paths */
    // Avalon❗✔️:    for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    // Avalon❗✔️:    {
    for (atom_index, atom) in molecule.atoms.iter_mut().enumerate() {
        // Avalon❗✔️:       if (atom_status[i] <= 0) ap->color = 0;
        if state.atom_status[atom_index] <= 0 {
            atom.color = 0;
        }
        // Avalon❗✔️:       if (i+1 == exclude_atom) ap->color = 0;
        if atom_index as i32 + 1 == exclude_atom {
            atom.color = 0;
        }
        // Avalon❗✔️:       if (ap->color == 0) continue;
        // Avalon❗✔️:    }
    }
    // Avalon❗✔️:    /* remove all bonds with only non-ring atoms from consideration */
    // Avalon❗✔️:    for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
    // Avalon❗✔️:    {
    for (bond_index, bond) in molecule.bonds.iter_mut().enumerate() {
        let first = bond.atoms[0] as usize - 1;
        let second = bond.atoms[1] as usize - 1;
        // Avalon❗✔️:       if (bond_status[i] <= 0  &&
        // Avalon❗✔️:           atom_status[bp->atoms[0]-1] == 0  &&
        // Avalon❗✔️:           atom_status[bp->atoms[1]-1] == 0) bp->color = 0;
        if state.bond_status[bond_index] <= 0 && state.atom_status[first] == 0 && state.atom_status[second] == 0 {
            bond.color = 0;
        }
        // Avalon❗✔️:       if (bp->atoms[0] == exclude_atom) bp->color = 0;
        // Avalon❗✔️:       if (bp->atoms[1] == exclude_atom) bp->color = 0;
        if bond.atoms.contains(&exclude_atom) {
            bond.color = 0;
        }
        // Avalon❗✔️:    }
    }
}

fn count_ring_paths(
    molecule: &MoleculeState,
    state: &FingerprintPreprocessingState,
    counts: &mut [i32],
    exclude_atom: i32,
) -> i32 {
    // Avalon❗✔️: // fprintf(stderr, "USE_RING_PATH\n");
    // Avalon❗✔️:    if (which_bits & USE_RING_PATH)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       strcpy(prefix, "RP:");
    // Avalon❗✔️:       seed = RING_PATH_SEED;
    let parent_seed = RING_PATH_SEED;
    let mut result = 0_i32;
    let mut touched = vec![0_i32; molecule.atoms.len()];
    // Avalon❗✔️:       for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    // Avalon❗✔️:       {
    for (atom_index, atom) in molecule.atoms.iter().enumerate() {
        // Avalon❗✔️:          if (ap->color <= 0) continue;
        if atom.color <= 0 {
            continue;
        }
        // Avalon❗✔️:          touched_indices[i] = 1; /* updating */
        touched[atom_index] = 1;
        // Avalon❗✔️:          old_seed = seed;
        // Avalon❗✔️:          seed = NEXT_SEED(seed, ap->color);
        let mut seed = next_hash(parent_seed, atom.color as u64);
        // Avalon❗✔️:          if (0)         // Class disabled to save bit density
        // Avalon❗✔️:          result += SetPathBitsRec(mp, nbp,
        // Avalon❗✔️:                                   fp_counts, ncounts,
        // Avalon❗✔️:                                   seed, touched_indices,
        // Avalon❗✔️:                                   1,
        // Avalon❗✔️:                                   2, 3,
        // Avalon❗✔️:                                   i, 0, -1,
        // Avalon❗✔️:                                   PROCESS_CHAINS,
        // Avalon❗✔️:                                   exclude_atom,
        // Avalon❗✔️:                                   prefix);
        // Avalon❗✔️:          if (ap->color > 5  &&  ap->color < 10  &&  atom_status[i]>2)      // only start at common light atoms
        if atom.color > 5 && atom.color < 10 && state.atom_status[atom_index] > 2 {
            // Avalon❗✔️:          {
            // Avalon❗✔️:             seed = NEXT_SEED(seed, 61);
            seed = next_hash(seed, 61);
            // Avalon❗✔️:             result += SetPathBitsRec(mp, nbp,
            // Avalon❗✔️:                                      fp_counts, ncounts,
            // Avalon❗✔️:                                      seed, touched_indices,
            // Avalon❗✔️:                                      1,
            // Avalon❗✔️:                                      3, 8,  /* 3 to 8 */
            // Avalon❗✔️:                                      i, 0, -1,
            // Avalon❗✔️:                                      STOP_AT_HEAVY_ATOM |
            // Avalon❗✔️:                                      PROCESS_RING_CLOSURES,
            // Avalon❗✔️:                                      exclude_atom,
            // Avalon❗✔️:                                      prefix);
            result += set_path_bits_rec(
                molecule,
                &state.neighbours,
                counts,
                seed,
                &mut touched,
                1,
                3,
                8,
                atom_index,
                0,
                -1,
                PathFlags::STOP_AT_HEAVY_ATOM | PathFlags::PROCESS_RING_CLOSURES,
                exclude_atom,
            );
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

pub(super) fn prepare_bond_path_state(molecule: &mut MoleculeState, exclude_atom: i32) {
    // Avalon❗✔️:    /* Set the color property to represent all different atom types */
    // Avalon❗✔️:    for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    // Avalon❗✔️:    {
    for (atom_index, atom) in molecule.atoms.iter_mut().enumerate() {
        // Avalon❗✔️:       ap->color = StringToInt(periodic_table, ap->atom_symbol);
        atom.color = avalon_atomic_number(&atom.atom_symbol);
        // Avalon❗✔️:       if (ap->color == 1  ||  i+1 == exclude_atom)
        if atom.color == 1 || atom_index as i32 + 1 == exclude_atom {
            // Avalon❗✔️:          ap->color = 0;   /* ignore hydrogens */
            atom.color = 0;
            // Avalon❗✔️:       else
        } else {
            // Avalon❗✔️:          ap->color = ANY_COLOR; /* treat all other atoms alike */
            atom.color = 113;
        }
        // Avalon❗✔️:    }
    }

    // Avalon❗✔️:    /* Set the color property to represent the different bond type classes */
    // Avalon❗✔️:    for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
    // Avalon❗✔️:    {
    for bond in &mut molecule.bonds {
        // Avalon❗✔️:       if (bp->bond_type == SINGLE)        bp->color = 1;
        // Avalon❗✔️:       else if (bp->bond_type == DOUBLE)   bp->color = 2;
        // Avalon❗✔️:       else if (bp->bond_type == TRIPLE)   bp->color = 3;
        // Avalon❗✔️:       else if (bp->bond_type == AROMATIC) bp->color = 4;
        // Avalon❗✔️:       else                                bp->color = 0;
        bond.color = match bond.bond_type {
            1 => 1,
            2 => 2,
            3 => 3,
            4 => 4,
            _ => 0,
        };
        // Avalon❗✔️:       if (bp->atoms[0] == exclude_atom) bp->color = 0;
        // Avalon❗✔️:       if (bp->atoms[1] == exclude_atom) bp->color = 0;
        if bond.atoms.contains(&exclude_atom) {
            bond.color = 0;
        }
        // Avalon❗✔️:    }
    }
}

pub(super) fn prepare_atom_class_state(molecule: &mut MoleculeState, exclude_atom: i32) {
    // Avalon❗✔️:    /* Set the color property to represent the different atom type classes */
    // Avalon❗✔️:    for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    // Avalon❗✔️:    {
    for atom_index in 0..molecule.atoms.len() {
        // Avalon❗✔️:       if (i+1 == exclude_atom)
        // Avalon❗✔️:       {
        if atom_index as i32 + 1 == exclude_atom {
            // Avalon❗✔️:          ap->color = 0;
            molecule.atoms[atom_index].color = 0;
            // Avalon❗✔️:          continue;
            continue;
            // Avalon❗✔️:       }
        }
        let symbol = molecule.atoms[atom_index].atom_symbol.as_str();
        // Avalon❗✔️:       if (0 == strcmp(ap->atom_symbol,"L"))     // Check if the atom list
        // Avalon❗✔️:       {                                         // is completely in one class
        if symbol == "L" {
            // Avalon❗✔️:          ap->color = 0;
            let mut color = 0_i32;
            // Avalon❗✔️:          symbol_lists = mp->symbol_lists;
            // Avalon❗✔️:          while (!IsNULL(symbol_lists))
            // Avalon❗✔️:          {
            for symbol_list in &molecule.symbol_lists {
                // Avalon❗✔️:             if (symbol_lists->atom == i+1)      /* found a match */
                if symbol_list.atom != atom_index as i32 + 1 {
                    continue;
                }
                // Avalon❗✔️:             {
                // Avalon❗✔️:                ap->color = 0;
                color = 0;
                // Avalon❗✔️:                if (symbol_lists->logic == EXCLUSIVE) break;
                if !symbol_list.inclusive {
                    break;
                }
                // Avalon❗✔️:                /* Step through all tokens in symbol list */
                // Avalon❗✔️:                strcpy(pat_buf,symbol_lists->string);
                // Avalon❗✔️:                for (tokp = strtok(pat_buf,",");
                // Avalon❗✔️:                     !IsNULL(tokp);
                // Avalon❗✔️:                     tokp = strtok((char *)NULL,","))
                // Avalon❗✔️:                {
                for token in symbol_list.symbols.split(',').filter(|token| !token.is_empty()) {
                    // Avalon❗✔️: #define HYDROGEN_FOUND 1
                    // Avalon❗✔️: #define CARBON_FOUND   2
                    // Avalon❗✔️: #define ONS_FOUND      4
                    // Avalon❗✔️: #define OTHER_FOUND    8
                    // Avalon❗✔️:                   if      (0 == strcmp("H",tokp)) ap->color |= HYDROGEN_FOUND;
                    // Avalon❗✔️:                   else if (0 == strcmp("C",tokp)) ap->color |= CARBON_FOUND;
                    // Avalon❗✔️:                   else if (0 == strcmp("O",tokp)) ap->color |= ONS_FOUND;
                    // Avalon❗✔️:                   else if (0 == strcmp("N",tokp)) ap->color |= ONS_FOUND;
                    // Avalon❗✔️:                   else if (0 == strcmp("S",tokp)) ap->color |= ONS_FOUND;
                    // Avalon❗✔️:                   else                            ap->color |= OTHER_FOUND;
                    color |= match token {
                        "H" => 1,
                        "C" => 2,
                        "O" | "N" | "S" => 4,
                        _ => 8,
                    };
                    // Avalon❗✔️:                }
                }
                // Avalon❗✔️:                if      (ap->color == CARBON_FOUND) ap->color = 6;
                // Avalon❗✔️:                else if (ap->color == ONS_FOUND)    ap->color = 8;
                // Avalon❗✔️:                else if (ap->color == OTHER_FOUND)  ap->color = 8;
                // Avalon❗✔️:                else                                ap->color = 0;
                color = match color {
                    2 => 6,
                    4 | 8 => 8,
                    _ => 0,
                };
                // Avalon❗✔️:             }
                // Avalon❗✔️:             symbol_lists = symbol_lists->next;
                // Avalon❗✔️:          }
            }
            molecule.atoms[atom_index].color = color;
            // Avalon❗✔️:       }
            // Avalon❗✔️:       else
        } else {
            // Avalon❗✔️:          if (0 == strcmp("H", ap->atom_symbol))
            // Avalon❗✔️:             ap->color = 0;      /* ignore hydrogens */
            // Avalon❗✔️:          else if (0 == strcmp("D", ap->atom_symbol))
            // Avalon❗✔️:             ap->color = 0;      /* ignore hydrogens */
            // Avalon❗✔️:          else if (0 == strcmp("T", ap->atom_symbol))
            // Avalon❗✔️:             ap->color = 0;      /* ignore hydrogens */
            // Avalon❗✔️:          else if (0 == strcmp("C", ap->atom_symbol))
            // Avalon❗✔️:             ap->color = 6;      /* carbon second row elements are one class */
            // Avalon❗✔️:          else if (0 == strcmp("N", ap->atom_symbol))
            // Avalon❗✔️:             ap->color = 8;      /* nitrogen, oxigen, and sulfur are one class */
            // Avalon❗✔️:          else if (0 == strcmp("O", ap->atom_symbol))
            // Avalon❗✔️:             ap->color = 8;      /* nitrogen, oxigen, and sulfur are one class */
            // Avalon❗✔️:          else if (0 == strcmp("S", ap->atom_symbol))
            // Avalon❗✔️:             ap->color = 8;      /* nitrogen, oxigen, and sulfur are one class */
            // Avalon❗✔️:          else if (0 == strcmp("Q", ap->atom_symbol))
            // Avalon❗✔️:             ap->color = 8;      /* nitrogen, oxigen, and sulfur are one class */
            // Avalon❗✔️:          else if (0 == strcmp("A", ap->atom_symbol))
            // Avalon❗✔️:             ap->color = 0;
            // Avalon❗✔️:          else
            // Avalon❗✔️:          {
            let color = match symbol {
                "H" | "D" | "T" | "A" => 0,
                "C" => 6,
                "N" | "O" | "S" | "Q" => 8,
                _ => {
                    // Avalon❗✔️:             tmp = StringToInt(periodic_table, ap->atom_symbol);
                    let atomic_number = avalon_atomic_number(symbol);
                    // Avalon❗✔️:             if (1 < tmp  &&  tmp < 115)
                    // Avalon❗✔️:                ap->color = 8;     /* non-carbon is the only second class */
                    // Avalon❗✔️:             else        /* This could be R atoms or other odd things */
                    // Avalon❗✔️:                ap->color = 0;
                    i32::from(atomic_number > 1 && atomic_number < 115) * 8
                }
            };
            // Avalon❗✔️:          }
            molecule.atoms[atom_index].color = color;
        }
        // Avalon❗✔️:       if (ap->color > 115) ap->color = 0;       /* ignore special atom types */
        if molecule.atoms[atom_index].color > 115 {
            molecule.atoms[atom_index].color = 0;
        }
        // Avalon❗✔️:    }
    }
}

fn count_atom_class_paths_first_pass(
    molecule: &MoleculeState,
    state: &FingerprintPreprocessingState,
    counts: &mut [i32],
    exclude_atom: i32,
) -> i32 {
    // Avalon❗✔️: // fprintf(stderr, "USE_ATOM_CLASS_PATH\n");
    // Avalon❗✔️:    /* Do a small first pass with defined bond types */
    // Avalon❗✔️:    if (which_bits & USE_ATOM_CLASS_PATH)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       strcpy(prefix, "ACP:");
    // Avalon❗✔️:       // seed = ATOM_CLASS_PATH_SEED+117;
    // Avalon❗✔️:       seed = NEXT_SEED(ATOM_CLASS_PATH_SEED,117);
    let parent_seed = next_hash(ATOM_CLASS_PATH_SEED, 117);
    let mut touched = vec![0_i32; molecule.atoms.len()];
    let mut result = 0_i32;
    // Avalon❗✔️:       for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    // Avalon❗✔️:       {
    for (atom_index, atom) in molecule.atoms.iter().enumerate() {
        // Avalon❗✔️:          if (ap->color <= 0) continue;
        // Avalon❗✔️:          if (i+1 == exclude_atom) continue;
        // Avalon❗✔️:          if (ap->color == 6  &&  degree[i] < 3) continue;
        if atom.color <= 0 || atom_index as i32 + 1 == exclude_atom || (atom.color == 6 && state.degree[atom_index] < 3)
        {
            continue;
        }
        // Avalon❗✔️:          touched_indices[i] = 1; /* updating */
        touched[atom_index] = 1;
        // Avalon❗✔️:          old_seed = seed;
        // Avalon❗✔️:          seed = NEXT_SEED(seed, ap->color);
        let seed = next_hash(parent_seed, atom.color as u64);
        // Avalon❗✔️:          if (ap->color == 6)
        // Avalon❗✔️:          {
        // Avalon❗✔️:             if (0)         // Class disabled to save bit density
        // Avalon❗✔️:             result += SetPathBitsRec(mp, nbp,
        // Avalon❗✔️:                                      fp_counts, ncounts,
        // Avalon❗✔️:                                      seed, touched_indices,
        // Avalon❗✔️:                                      1,
        // Avalon❗✔️:                                      3, 4,
        // Avalon❗✔️:                                      i, 0, -1,
        // Avalon❗✔️:                                      FORCED_HETERO_END |
        // Avalon❗✔️:                                      IGNORE_PATH_SYMBOL |
        // Avalon❗✔️:                                      PROCESS_CHAINS,
        // Avalon❗✔️:                                      exclude_atom,
        // Avalon❗✔️:                                      prefix);
        // Avalon❗✔️:          }
        // Avalon❗✔️:          else
        // Avalon❗✔️:          {
        if atom.color != 6 {
            // Avalon❗✔️:             if (0)         // Class disabled to save bit density
            // Avalon❗✔️:             result += SetPathBitsRec(mp, nbp,
            // Avalon❗✔️:                                      fp_counts, ncounts,
            // Avalon❗✔️:                                      seed, touched_indices,
            // Avalon❗✔️:                                      1,
            // Avalon❗✔️:                                      2, 3,  /* path length 3 to 3 */
            // Avalon❗✔️:                                      i, 0, -1,
            // Avalon❗✔️:                                      // IGNORE_PATH_SYMBOL |
            // Avalon❗✔️:                                      PROCESS_CHAINS,
            // Avalon❗✔️:                                      exclude_atom,
            // Avalon❗✔️:                                      prefix);
            // Avalon❗✔️:             if (0)         // Class disabled to save bit density
            // Avalon❗✔️:             result += SetPathBitsRec(mp, nbp,
            // Avalon❗✔️:                                      fp_counts, ncounts,
            // Avalon❗✔️:                                      seed, touched_indices,
            // Avalon❗✔️:                                      1,
            // Avalon❗✔️:                                      4, 5,  /* path length 4 to 5 */
            // Avalon❗✔️:                                      i, 0, -1,
            // Avalon❗✔️:                                      IGNORE_PATH_SYMBOL |
            // Avalon❗✔️:                                      FORCED_HETERO_END |
            // Avalon❗✔️:                                      PROCESS_CHAINS,
            // Avalon❗✔️:                                      exclude_atom,
            // Avalon❗✔️:                                      prefix);
            // Avalon❗✔️:             result += SetPathBitsRec(mp, nbp,
            // Avalon❗✔️:                                      fp_counts, ncounts,
            // Avalon❗✔️:                                      seed, touched_indices,
            // Avalon❗✔️:                                      1,
            // Avalon❗✔️:                                      3, 4,  /* path length 4 to 7 */
            // Avalon❗✔️:                                      i, 0, -1,
            // Avalon❗✔️:                                      FORCED_RING_PATH |
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
                1,
                3,
                4,
                atom_index,
                0,
                -1,
                PathFlags::FORCED_RING_PATH | PathFlags::PROCESS_RING_CLOSURES | PathFlags::PROCESS_CHAINS,
                exclude_atom,
            );
            // Avalon❗✔️:          }
        }
        // Avalon❗✔️:          touched_indices[i] = 0;
        touched[atom_index] = 0;
        // Avalon❗✔️:          seed = old_seed;
        // Avalon❗✔️:       }
    }
    // Avalon❗✔️:    }
    result
}

fn unify_bond_classes(molecule: &mut MoleculeState, exclude_atom: i32) {
    // Avalon❗✔️:    /* Set the color property to only a single class */
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

fn count_atom_class_paths_second_pass(
    molecule: &MoleculeState,
    state: &FingerprintPreprocessingState,
    counts: &mut [i32],
    exclude_atom: i32,
) -> i32 {
    // Avalon❗✔️: // fprintf(stderr, "USE_ATOM_CLASS_PATH\n");
    // Avalon❗✔️:    // Here, we've unified atom types to carbon/hetero and made all bond types
    // Avalon❗✔️:    // identical
    // Avalon❗✔️:    if (which_bits & USE_ATOM_CLASS_PATH)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       strcpy(prefix, "ACP:");
    // Avalon❗✔️:       seed = ATOM_CLASS_PATH_SEED;
    let parent_seed = ATOM_CLASS_PATH_SEED;
    let mut touched = vec![0_i32; molecule.atoms.len()];
    let mut result = 0_i32;
    // Avalon❗✔️:       for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    // Avalon❗✔️:       {
    for (atom_index, atom) in molecule.atoms.iter().enumerate() {
        // Avalon❗✔️:          if (ap->color <= 0) continue;
        // Avalon❗✔️:          if (i+1 == exclude_atom) continue;
        // Avalon❗✔️:          if (ap->color == 6) continue;
        if atom.color <= 0 || atom_index as i32 + 1 == exclude_atom || atom.color == 6 {
            continue;
        }
        // Avalon❗✔️:          touched_indices[i] = 1; /* updating */
        touched[atom_index] = 1;
        // Avalon❗✔️:          old_seed = seed;
        // Avalon❗✔️:          seed = NEXT_SEED(seed, ap->color);
        let mut seed = next_hash(parent_seed, atom.color as u64);
        // Avalon❗✔️:          // if (ap->color == 6)
        // Avalon❗✔️:          {
        // Avalon❗✔️:             if (0*degree[i] > 2)         // Class disabled to save bit density
        // Avalon❗✔️:                result += SetPathBitsRec(mp, nbp,
        // Avalon❗✔️:                                         fp_counts, ncounts,
        // Avalon❗✔️:                                         seed, touched_indices,
        // Avalon❗✔️:                                         1,
        // Avalon❗✔️:                                         3, 3,  /* path length 3 to 3 */
        // Avalon❗✔️:                                         i, 0, -1,
        // Avalon❗✔️:                                         FORCED_HETERO_END |
        // Avalon❗✔️:                                         PROCESS_CHAINS,
        // Avalon❗✔️:                                         exclude_atom,
        // Avalon❗✔️:                                         prefix);
        // Avalon❗✔️:          }
        // Avalon❗✔️:          // else
        // Avalon❗✔️:          {
        // Avalon❗✔️:             if(0)         // Class disabled to save bit density
        // Avalon❗✔️:             result += SetPathBitsRec(mp, nbp,
        // Avalon❗✔️:                                      fp_counts, ncounts,
        // Avalon❗✔️:                                      seed, touched_indices,
        // Avalon❗✔️:                                      1,
        // Avalon❗✔️:                                      3, 4,  /* path length 3 to 4 */
        // Avalon❗✔️:                                      i, 0, -1,
        // Avalon❗✔️:                                      IGNORE_PATH_SYMBOL |
        // Avalon❗✔️:                                      PROCESS_CHAINS,
        // Avalon❗✔️:                                      exclude_atom,
        // Avalon❗✔️:                                      prefix);
        // Avalon❗✔️:             seed = NEXT_SEED(seed, 23+ap->color*19);
        seed = next_hash(seed, (23 + atom.color * 19) as u64);
        // Avalon❗✔️:             result += SetPathBitsRec(mp, nbp,
        // Avalon❗✔️:                                      fp_counts, ncounts,
        // Avalon❗✔️:                                      seed, touched_indices,
        // Avalon❗✔️:                                      1,
        // Avalon❗✔️:                                      3, 9,  /* path length 3 to 9 */
        // Avalon❗✔️:                                      i, 0, -1,
        // Avalon❗✔️:                                      IGNORE_PATH_SYMBOL |
        // Avalon❗✔️:                                      PROCESS_RING_CLOSURES,
        // Avalon❗✔️:                                      exclude_atom,
        // Avalon❗✔️:                                      prefix);
        result += set_path_bits_rec(
            molecule,
            &state.neighbours,
            counts,
            seed,
            &mut touched,
            1,
            3,
            9,
            atom_index,
            0,
            -1,
            PathFlags::IGNORE_PATH_SYMBOL | PathFlags::PROCESS_RING_CLOSURES,
            exclude_atom,
        );
        // Avalon❗✔️:          }
        // Avalon❗✔️:          seed = old_seed;
        // Avalon❗✔️:          touched_indices[i] = 0; /* down-dating */
        touched[atom_index] = 0;
        // Avalon❗✔️:       }
    }
    // Avalon❗✔️:       /* Q-Q and Q-C ring bond count */
    // Avalon❗✔️:       qq_count=0;
    let mut qq_count = 0_i32;
    // Avalon❗✔️:       qc_count=0;
    let mut qc_count = 0_i32;
    // Avalon❗✔️:       for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
    // Avalon❗✔️:       {
    for (bond_index, bond) in molecule.bonds.iter().enumerate() {
        // Avalon❗✔️:          if (bond_status[i] == 0) continue;
        // Avalon❗✔️:          if (bp->color == 0) continue;
        // Avalon❗✔️:          if (bp->atoms[0] == exclude_atom) continue;
        // Avalon❗✔️:          if (bp->atoms[1] == exclude_atom) continue;
        if state.bond_status[bond_index] == 0 || bond.color == 0 || bond.atoms.contains(&exclude_atom) {
            continue;
        }
        // Avalon❗✔️:          ai1 = mp->atom_array[bp->atoms[0]-1].color;
        let first_color = molecule.atoms[bond.atoms[0] as usize - 1].color;
        // Avalon❗✔️:          ai2 = mp->atom_array[bp->atoms[1]-1].color;
        let second_color = molecule.atoms[bond.atoms[1] as usize - 1].color;
        // Avalon❗✔️:          if (ai1 == 0  ||  ai2 == 0) continue;
        // Avalon❗✔️:          if (ai1 == 6  &&  ai2 == 6) continue;
        if first_color == 0 || second_color == 0 || (first_color == 6 && second_color == 6) {
            continue;
        }
        // Avalon❗✔️:          if (ai1 != 6  &&  ai2 != 6)
        if first_color != 6 && second_color != 6 {
            // Avalon❗✔️:          {
            // Avalon❗✔️:             qq_count++;
            qq_count += 1;
            // Avalon❗✔️:             for (j=3; j<9; j++)        /* set bits for not too large ring size */
            for ring_size in 3..9 {
                // Avalon❗✔️:                if (bp->rsize_flags&(1<<j))
                // Avalon❗✔️:                {
                if bond.rsize_flags & (1 << ring_size) != 0 {
                    // Avalon❗✔️:                   ADD_BIT(fp_counts, ncounts, NEXT_SEED(ATOM_CLASS_PATH_SEED*17, j*8));
                    add_bit(counts, next_hash(ATOM_CLASS_PATH_SEED * 17, (ring_size * 8) as u64));
                    // Avalon❗✔️:                   result++;
                    result += 1;
                    // Avalon❗✔️:                }
                }
            }
            // Avalon❗✔️:          }
            // Avalon❗✔️:          else
        } else {
            // Avalon❗✔️:          {
            // Avalon❗✔️:             qc_count++;
            qc_count += 1;
            // Avalon❗✔️:             for (j=3; j<9; j++)        /* set bits for not too large ring size */
            for ring_size in 3..9 {
                // Avalon❗✔️:                if (bp->rsize_flags&(1<<j))
                // Avalon❗✔️:                {
                if bond.rsize_flags & (1 << ring_size) != 0 {
                    // Avalon❗✔️:                   ADD_BIT(fp_counts, ncounts, NEXT_SEED(ATOM_CLASS_PATH_SEED*19, j*8));
                    add_bit(counts, next_hash(ATOM_CLASS_PATH_SEED * 19, (ring_size * 8) as u64));
                    // Avalon❗✔️:                   result++;
                    result += 1;
                    // Avalon❗✔️:                }
                }
            }
            // Avalon❗✔️:          }
        }
        // Avalon❗✔️:       }
    }
    // Avalon❗✔️:       seed = 2*ATOM_CLASS_PATH_SEED+3;
    let mut seed = 2 * ATOM_CLASS_PATH_SEED + 3;
    // Avalon❗✔️:       for (i=1; i<=qq_count; i=(int)(1+i*1.5))
    let mut count_level = 1_i32;
    while count_level <= qq_count {
        // Avalon❗✔️:          seed = NEXT_SEED(seed, i*153);
        seed = next_hash(seed, (count_level * 153) as u64);
        // Avalon❗✔️:          ADD_BIT(fp_counts, ncounts, seed);
        add_bit(counts, seed);
        // Avalon❗✔️:          result++;
        result += 1;
        // Avalon❗✔️:          seed = NEXT_SEED(seed, 53);
        seed = next_hash(seed, 53);
        // Avalon❗✔️:          if (i <= 1)
        if count_level <= 1 {
            // Avalon❗✔️:             ADD_BIT(fp_counts, ncounts, seed);
            add_bit(counts, seed);
            // Avalon❗✔️:             result++;
            result += 1;
        }
        count_level = 1 + (f64::from(count_level) * 1.5) as i32;
    }
    // Avalon❗✔️:       seed = 3*ATOM_CLASS_PATH_SEED+5;
    seed = 3 * ATOM_CLASS_PATH_SEED + 5;
    // Avalon❗✔️:       for (i=1; i<=MIN(qc_count,2); i++)
    for level in 1..=qc_count.min(2) {
        // Avalon❗✔️:          seed = NEXT_SEED(seed, i*157);
        seed = next_hash(seed, (level * 157) as u64);
        // Avalon❗✔️:          ADD_BIT(fp_counts, ncounts, seed);
        add_bit(counts, seed);
        // Avalon❗✔️:          result++;
        result += 1;
    }
    // Avalon❗✔️:       for (i=3; i<=qc_count; i=(int)(i*1.8))
    count_level = 3;
    while count_level <= qc_count {
        // Avalon❗✔️:          seed = NEXT_SEED(seed, i*157);
        seed = next_hash(seed, (count_level * 157) as u64);
        // Avalon❗✔️:          ADD_BIT(fp_counts, ncounts, seed);
        add_bit(counts, seed);
        // Avalon❗✔️:          result++;
        result += 1;
        count_level = (f64::from(count_level) * 1.8) as i32;
    }
    // Avalon❗✔️:    }
    result
}

fn prepare_ring_pattern_state(molecule: &mut MoleculeState, state: &FingerprintPreprocessingState, exclude_atom: i32) {
    // Avalon❗✔️:    /* Compute ring patters with at least one cycle */
    // Avalon❗✔️:    /* remove bonds from consideration that don't have at least one ring atom */
    // Avalon❗✔️:    for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
    // Avalon❗✔️:    {
    for bond in &mut molecule.bonds {
        // Avalon❗✔️:       if (bp->atoms[0] == exclude_atom) continue;
        // Avalon❗✔️:       if (bp->atoms[1] == exclude_atom) continue;
        if bond.atoms.contains(&exclude_atom) {
            continue;
        }
        // Avalon❗✔️:       bp->color = 5;
        bond.color = 5;
        let first = bond.atoms[0] as usize - 1;
        let second = bond.atoms[1] as usize - 1;
        // Avalon❗✔️:       /* ignore non-ring bonds */
        // Avalon❗✔️:       if (atom_status[bp->atoms[0]-1] == 0  &&
        // Avalon❗✔️:           atom_status[bp->atoms[1]-1] == 0)
        if state.atom_status[first] == 0 && state.atom_status[second] == 0 {
            // Avalon❗✔️:          bp->color = 0;
            bond.color = 0;
            // Avalon❗✔️:       else
        } else {
            // Avalon❗✔️:          // may be redundant since atom types have already been unified
            // Avalon❗✔️:          ap = &mp->atom_array[bp->atoms[0]-1];
            // Avalon❗✔️:          if (ap->color != 0)
            // Avalon❗✔️:          {
            if molecule.atoms[first].color != 0 {
                // Avalon❗✔️:             if (ap->color != 6) ap->color = 8;  /* non-carbons in one class */
                if molecule.atoms[first].color != 6 {
                    molecule.atoms[first].color = 8;
                }
                // Avalon❗✔️:             if (ap->rsize_flags == 0) ap->color = 0;
                if molecule.atoms[first].rsize_flags == 0 {
                    molecule.atoms[first].color = 0;
                }
                // Avalon❗✔️:          }
            }
            // Avalon❗✔️:          // may be redundant since atom types have already been unified
            // Avalon❗✔️:          ap = &mp->atom_array[bp->atoms[1]-1];
            // Avalon❗✔️:          if (ap->color != 0)
            // Avalon❗✔️:          {
            if molecule.atoms[second].color != 0 {
                // Avalon❗✔️:             if (ap->color != 6) ap->color = 8;  /* non-carbons in one class */
                if molecule.atoms[second].color != 6 {
                    molecule.atoms[second].color = 8;
                }
                // Avalon❗✔️:             if (ap->rsize_flags == 0) ap->color = 0;
                if molecule.atoms[second].rsize_flags == 0 {
                    molecule.atoms[second].color = 0;
                }
                // Avalon❗✔️:          }
            }
            // Avalon❗✔️:       }
        }
        // Avalon❗✔️:    }
    }
}

fn count_ring_patterns(
    molecule: &mut MoleculeState,
    state: &FingerprintPreprocessingState,
    counts: &mut [i32],
    exclude_atom: i32,
) -> i32 {
    // Avalon❗✔️: // fprintf(stderr, "USE_RING_PATTERN\n");
    // Avalon❗✔️:    // Here, we have bond type ignored and atom types mapped to C and Q
    // Avalon❗✔️:    if (which_bits & USE_RING_PATTERN)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       strcpy(prefix, "RPT:");
    // Avalon❗✔️:       /* first process ring bond paths with atom classes */
    // Avalon❗✔️:       seed = RING_PATTERN_SEED;
    let parent_seed = RING_PATTERN_SEED;
    let mut touched = vec![0_i32; molecule.atoms.len()];
    let mut result = 0_i32;
    // Avalon❗✔️:       for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    for atom_index in 0..molecule.atoms.len() {
        let atom = &molecule.atoms[atom_index];
        // Avalon❗✔️:          if (ap->color <= 0) continue;
        // Avalon❗✔️:          if (i+1 == exclude_atom) continue;
        // Avalon❗✔️:          /* Don't process fragments starting at carbon in 6-ring only */
        // Avalon❗✔️:          if (ap->color == 6  &&  0 == (ap->rsize_flags&SPECIAL_RING)) continue;
        if atom.color <= 0
            || atom_index as i32 + 1 == exclude_atom
            || (atom.color == 6 && atom.rsize_flags & SPECIAL_RING == 0)
        {
            continue;
        }
        touched[atom_index] = 1;
        let seed = next_hash(parent_seed, atom.color as u64);
        // Avalon❗✔️:          result += SetPathBitsRec(mp, nbp,
        // Avalon❗✔️:                                   fp_counts, ncounts,
        // Avalon❗✔️:                                   seed, touched_indices,
        // Avalon❗✔️:                                   1,
        // Avalon❗✔️:                                   3, 3,  /* ring bond path size 3 to 3 */
        // Avalon❗✔️:                                   i, 0, -1,
        // Avalon❗✔️:                                   PROCESS_CHAINS,
        // Avalon❗✔️:                                   exclude_atom,
        // Avalon❗✔️:                                   prefix);
        result += set_path_bits_rec(
            molecule,
            &state.neighbours,
            counts,
            seed,
            &mut touched,
            1,
            3,
            3,
            atom_index,
            0,
            -1,
            PathFlags::PROCESS_CHAINS,
            exclude_atom,
        );
        touched[atom_index] = 0;
    }

    // Avalon❗✔️:       /* Now, we only include complete rings but ignore atom-type */
    // Avalon❗✔️:       /* 'A' atoms are now included nodes */
    // Avalon❗✔️:       for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    // Avalon❗✔️:       {
    for (atom_index, atom) in molecule.atoms.iter_mut().enumerate() {
        // Avalon❗✔️:          /* add 'A' atom to standard class */
        // Avalon❗✔️:          if (0 == strcmp("A", ap->atom_symbol))
        // Avalon❗✔️:          {
        if atom.atom_symbol == "A" {
            // Avalon❗✔️:             ap->color = 9;
            atom.color = 9;
            // Avalon❗✔️:          }
        }
        // Avalon❗✔️:          if (i+1 == exclude_atom) ap->color = 0;
        if atom_index as i32 + 1 == exclude_atom {
            atom.color = 0;
        }
        // Avalon❗✔️:          if (ap->color == 0) continue;
        if atom.color == 0 {
            continue;
        }
        // Avalon❗✔️:          ap->color = 9;    /* all ring atoms in same class */
        atom.color = 9;
        // Avalon❗✔️:       }
    }
    // Avalon❗✔️:       for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
    // Avalon❗✔️:       {
    for (bond_index, bond) in molecule.bonds.iter_mut().enumerate() {
        // Avalon❗✔️:          if (bp->atoms[0] == exclude_atom) bp->color = 0;
        // Avalon❗✔️:          if (bp->atoms[1] == exclude_atom) bp->color = 0;
        if bond.atoms.contains(&exclude_atom) {
            bond.color = 0;
        }
        // Avalon❗✔️:          if (bond_status[i] <= 0) bp->color = 0;
        if state.bond_status[bond_index] <= 0 {
            bond.color = 0;
        }
        // Avalon❗✔️:       }
    }
    // Avalon❗✔️: // seed = RING_PATTERN_SEED+23;
    // Avalon❗✔️:       seed = NEXT_SEED(RING_PATTERN_SEED, 23);
    let mut seed = next_hash(RING_PATTERN_SEED, 23);
    // Avalon❗✔️:       for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    // Avalon❗✔️:       {
    for atom_index in 0..molecule.atoms.len() {
        let atom = &molecule.atoms[atom_index];
        // Avalon❗✔️:          if (ap->color == 0) continue;
        if atom.color == 0 {
            continue;
        }
        // Avalon❗✔️:          if (i+1 == exclude_atom) continue;
        if atom_index as i32 + 1 == exclude_atom {
            continue;
        }
        // Avalon❗✔️:          touched_indices[i] = 1; /* updating */
        touched[atom_index] = 1;
        // Avalon❗✔️:          old_seed = seed;
        let old_seed = seed;
        // Avalon❗✔️:          seed = NEXT_SEED(seed, ap->color);
        seed = next_hash(seed, atom.color as u64);
        // Avalon❗✔️:          if (0)         // Class disabled to save bit density
        // Avalon❗✔️:          result += SetPathBitsRec(mp, nbp,
        // Avalon❗✔️:                                   fp_counts, ncounts,
        // Avalon❗✔️:                                   seed, touched_indices,
        // Avalon❗✔️:                                   1,
        // Avalon❗✔️:                                   4, 17,  /* ring size 4 to 17 */
        // Avalon❗✔️:                                   i, 0, -1,
        // Avalon❗✔️:                                   IGNORE_PATH_SYMBOL |
        // Avalon❗✔️:                                   PROCESS_RING_CLOSURES,
        // Avalon❗✔️:                                   exclude_atom,
        // Avalon❗✔️:                                   prefix);

        // Avalon❗✔️:          seed = old_seed;
        seed = old_seed;
        // Avalon❗✔️:          if (0)         // Class disabled to save bit density
        // Avalon❗✔️:          if (atom_status[i] > 2)        // start at ring fusion
        // Avalon❗✔️:          {
        // Avalon❗✔️:             seed = NEXT_SEED(seed, 61);
        // Avalon❗✔️:             result += SetPathBitsRec(mp, nbp,
        // Avalon❗✔️:                                      fp_counts, ncounts,
        // Avalon❗✔️:                                      seed, touched_indices,
        // Avalon❗✔️:                                      1,
        // Avalon❗✔️:                                      6, 17,  /* ring path size 6 to 17 */
        // Avalon❗✔️:                                      i, 0, -1,
        // Avalon❗✔️:                                      IGNORE_PATH_SYMBOL |
        // Avalon❗✔️:                                      PROCESS_RING_CLOSURES,
        // Avalon❗✔️:                                      exclude_atom,
        // Avalon❗✔️:                                      prefix);
        // Avalon❗✔️:          }

        // Avalon❗✔️:          seed = old_seed;
        seed = old_seed;
        // Avalon❗✔️:          if (0)         // Class disabled to save bit density
        // Avalon❗✔️:          if (degree[i] > 2  && atom_status[i] >= 2)     // start at ring substituents
        // Avalon❗✔️:          {
        // Avalon❗✔️:             seed = NEXT_SEED(seed, 67);
        // Avalon❗✔️:             result += SetPathBitsRec(mp, nbp,
        // Avalon❗✔️:                                      fp_counts, ncounts,
        // Avalon❗✔️:                                      seed, touched_indices,
        // Avalon❗✔️:                                      1,
        // Avalon❗✔️:                                      6, 17,  /* ring path size 6 to 17 */
        // Avalon❗✔️:                                      i, 0, -1,
        // Avalon❗✔️:                                      IGNORE_PATH_SYMBOL |
        // Avalon❗✔️:                                      PROCESS_RING_CLOSURES,
        // Avalon❗✔️:                                      exclude_atom,
        // Avalon❗✔️:                                      prefix);
        // Avalon❗✔️:          }
        // Avalon❗✔️:          seed = old_seed;
        seed = old_seed;
        // Keep the final source assignment observable to the Rust compiler's
        // unused-assignment analysis without changing fingerprint behavior.
        let _ = seed;

        // Avalon❗✔️:          touched_indices[i] = 0; /* down-dating */
        touched[atom_index] = 0;
        // Avalon❗✔️:       }
    }
    // Avalon❗✔️:    }
    result
}

fn count_atom_symbol_paths(
    molecule: &MoleculeState,
    state: &FingerprintPreprocessingState,
    counts: &mut [i32],
    as_query: bool,
    exclude_atom: i32,
) -> i32 {
    // Avalon❗✔️: // fprintf(stderr, "USE_ATOM_SYMBOL_PATH\n");
    // Avalon❗✔️:    if (which_bits & USE_ATOM_SYMBOL_PATH)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       strcpy(prefix, "ASP:");
    // Avalon❗✔️:       seed = ATOM_SYMBOL_PATH_SEED;
    let parent_seed = ATOM_SYMBOL_PATH_SEED;
    let mut result = 0_i32;
    let mut touched = vec![0_i32; molecule.atoms.len()];
    // Avalon❗✔️:       for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    // Avalon❗✔️:       {
    for (atom_index, atom) in molecule.atoms.iter().enumerate() {
        // Avalon❗✔️:          if (i+1 == exclude_atom) continue;
        if atom_index as i32 + 1 == exclude_atom {
            continue;
        }
        // Avalon❗✔️:          if (ap->color <= 0) continue;
        if atom.color <= 0 {
            continue;
        }
        // Avalon❗✔️:          touched_indices[i] = 1; /* updating */
        touched[atom_index] = 1;
        // Avalon❗✔️:          old_seed = seed;
        // Avalon❗✔️:          seed = NEXT_SEED(seed, ap->color);
        let mut seed = next_hash(parent_seed, atom.color as u64);
        // Avalon❗✔️:          /* Ignore common atom types */
        // Avalon❗✔️:          if (ap->color != 6  &&  ap->color != 7  &&  ap->color != 8)
        // Avalon❗✔️:          {
        if !matches!(atom.color, 6 | 7 | 8) {
            // Avalon❗✔️:             ADD_BIT(fp_counts, ncounts, seed);
            add_bit(counts, seed);
            // Avalon❗✔️:             result++;
            result += 1;
            // Avalon❗✔️:          }
        }
        // Avalon❗✔️:          /* fingerprint two atom pairs if substitution count is defined */
        // Avalon❗✔️: if (1) // [TODO] 1
        // Avalon❗✔️:          if (ap->color !=6  &&
        // Avalon❗✔️:              (!as_query                 ||
        // Avalon❗✔️:               ap->sub_desc == SUB_AS_IS ||
        // Avalon❗✔️:               (ap->sub_desc != NONE     &&
        // Avalon❗✔️:                ap->sub_desc != SUB_MORE &&
        // Avalon❗✔️:                ap->sub_desc == degree[i]+SUB_ONE-1)))
        if atom.color != 6
            && (!as_query
                || atom.sub_desc == SUB_AS_IS
                || (atom.sub_desc != 0
                    && atom.sub_desc != SUB_MORE
                    && atom.sub_desc == state.degree[atom_index] + SUB_ONE - 1))
        {
            // Avalon❗✔️:          {
            // Avalon❗✔️:             result += SetPathBitsRec(mp, nbp,
            // Avalon❗✔️:                                      fp_counts, ncounts,
            // Avalon❗✔️:                                      seed+12347, touched_indices,
            // Avalon❗✔️:                                      1,
            // Avalon❗✔️:                                      1, 1,  /* path length 1 to 1 */
            // Avalon❗✔️:                                      i, 0, -1,
            // Avalon❗✔️:                                      FORCED_HETERO_END |
            // Avalon❗✔️:                                      PROCESS_CHAINS,
            // Avalon❗✔️:                                      exclude_atom,
            // Avalon❗✔️:                                      prefix);
            result += set_path_bits_rec(
                molecule,
                &state.neighbours,
                counts,
                seed.wrapping_add(12_347),
                &mut touched,
                1,
                1,
                1,
                atom_index,
                0,
                -1,
                PathFlags::FORCED_HETERO_END | PathFlags::PROCESS_CHAINS,
                exclude_atom,
            );
        }
        // Avalon❗✔️:          /* 2 bond paths for not very common starts */
        // Avalon❗✔️: if (1) // [TODO] 2
        // Avalon❗✔️:          if (ap->color >= 10  ||
        // Avalon❗✔️:              (ap->color == 7  &&  degree[i] > 0)  ||
        // Avalon❗✔️:              (ap->color == 8  &&  degree[i] > 1)  ||
        // Avalon❗✔️:              (ap->color == 6  &&  degree[i] > 2  &&  atom_status[i] > 0))
        if atom.color >= 10
            || (atom.color == 7 && state.degree[atom_index] > 0)
            || (atom.color == 8 && state.degree[atom_index] > 1)
            || (atom.color == 6 && state.degree[atom_index] > 2 && state.atom_status[atom_index] > 0)
        {
            // Avalon❗✔️:             result += SetPathBitsRec(mp, nbp,
            // Avalon❗✔️:                                      fp_counts, ncounts,
            // Avalon❗✔️:                                      seed, touched_indices,
            // Avalon❗✔️:                                      1,
            // Avalon❗✔️:                                      1, 2,  /* path length 1 to 2 */
            // Avalon❗✔️:                                      i, 0, -1,
            // Avalon❗✔️:                                      STOP_AT_HEAVY_ATOM |            // don't cross very heavy atoms
            // Avalon❗✔️:                                      FORCED_HETERO_END |
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
                2,
                atom_index,
                0,
                -1,
                PathFlags::STOP_AT_HEAVY_ATOM | PathFlags::FORCED_HETERO_END | PathFlags::PROCESS_CHAINS,
                exclude_atom,
            );
        }
        // Avalon❗✔️:          /* Add more paths starting at special atoms */
        // Avalon❗✔️: if (1) // [TODO] 3
        // Avalon❗✔️:          if (ap->color > 6)
        // Avalon❗✔️:          {
        if atom.color > 6 {
            // Avalon❗✔️:             result += SetPathBitsRec(mp, nbp,
            // Avalon❗✔️:                                      fp_counts, ncounts,
            // Avalon❗✔️:                                      217*ap->color+seed, touched_indices,
            // Avalon❗✔️:                                      1,
            // Avalon❗✔️:                                      3, 4,  /* path length 3 to 4 */
            // Avalon❗✔️:                                      i, 0, -1,
            // Avalon❗✔️:                                      IGNORE_PATH_SYMBOL |
            // Avalon❗✔️:                                      PROCESS_RING_CLOSURES |
            // Avalon❗✔️:                                      STOP_AT_HEAVY_ATOM |            // don't cross very heavy atoms
            // Avalon❗✔️:                                      PROCESS_CHAINS,
            // Avalon❗✔️:                                      exclude_atom,
            // Avalon❗✔️:                                      prefix);
            result += set_path_bits_rec(
                molecule,
                &state.neighbours,
                counts,
                seed.wrapping_add((217_i64 * i64::from(atom.color)) as u64),
                &mut touched,
                1,
                3,
                4,
                atom_index,
                0,
                -1,
                PathFlags::IGNORE_PATH_SYMBOL
                    | PathFlags::PROCESS_RING_CLOSURES
                    | PathFlags::STOP_AT_HEAVY_ATOM
                    | PathFlags::PROCESS_CHAINS,
                exclude_atom,
            );
            // Avalon❗✔️:             if (ap->color > 10  &&  ap->color <= 18)     // only third row of periodic table
            if atom.color > 10 && atom.color <= 18 {
                // Avalon❗✔️:                  result += SetPathBitsRec(mp, nbp,
                // Avalon❗✔️:                                           fp_counts, ncounts,
                // Avalon❗✔️:                                           17+seed, touched_indices,
                // Avalon❗✔️:                                           1,
                // Avalon❗✔️:                                           5, 7,
                // Avalon❗✔️:                                           i, 0, -1,
                // Avalon❗✔️:                                           FORCED_HETERO_END |
                // Avalon❗✔️:                                           IGNORE_PATH_SYMBOL |
                // Avalon❗✔️:                                           STOP_AT_HEAVY_ATOM |            // don't cross very heavy atoms
                // Avalon❗✔️:                                           // PROCESS_RING_CLOSURES |      // rings containing metal cause too many bits
                // Avalon❗✔️:                                           PROCESS_CHAINS,
                // Avalon❗✔️:                                           exclude_atom,
                // Avalon❗✔️:                                           prefix);
                result += set_path_bits_rec(
                    molecule,
                    &state.neighbours,
                    counts,
                    seed.wrapping_add(17),
                    &mut touched,
                    1,
                    5,
                    7,
                    atom_index,
                    0,
                    -1,
                    PathFlags::FORCED_HETERO_END
                        | PathFlags::IGNORE_PATH_SYMBOL
                        | PathFlags::STOP_AT_HEAVY_ATOM
                        | PathFlags::PROCESS_CHAINS,
                    exclude_atom,
                );
            }
            // Avalon❗✔️:          }
        }
        // Avalon❗✔️: if (1) // [TODO] 4
        // Avalon❗✔️:          if ((ap->color == 7 ||  ap->color == 8) &&      // only do this for common hetero elements
        // Avalon❗✔️:              cdegree[i]>2)
        // Avalon❗✔️:          {
        if matches!(atom.color, 7 | 8) && state.carbon_degree[atom_index] > 2 {
            // Avalon❗✔️:             seed = NEXT_SEED(seed, ap->color*23);
            seed = next_hash(seed, (atom.color * 23) as u64);
            // Avalon❗✔️:             result += SetPathBitsRec(mp, nbp,
            // Avalon❗✔️:                                      fp_counts, ncounts,
            // Avalon❗✔️:                                      seed, touched_indices,
            // Avalon❗✔️:                                      1,
            // Avalon❗✔️:                                      2, 6,  /* path length 1 to 4 */
            // Avalon❗✔️:                                      i, 0, -1,
            // Avalon❗✔️:                                      STOP_AT_HEAVY_ATOM |            // don't cross very heavy atoms
            // Avalon❗✔️:                                      IGNORE_TERM_SYMBOL |
            // Avalon❗✔️:                                      // IGNORE_PATH_SYMBOL |
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
                6,
                atom_index,
                0,
                -1,
                PathFlags::STOP_AT_HEAVY_ATOM | PathFlags::IGNORE_TERM_SYMBOL | PathFlags::PROCESS_CHAINS,
                exclude_atom,
            );
            // Avalon❗✔️:          }
        }
        // Avalon❗✔️:          seed = old_seed;
        // Avalon❗✔️:          /* Add another set of bits for paths starting at more crowded atoms */
        seed = parent_seed;
        // Avalon❗✔️: if (1) // [TODO] 5
        // Avalon❗✔️:          if (degree[i] >= 4  &&
        // Avalon❗✔️:              ap->color != 5  &&         // don't do it for boron
        // Avalon❗✔️:              ap->color < 18)            // don't do it for transition metals
        if state.degree[atom_index] >= 4 && atom.color != 5 && atom.color < 18 {
            // Avalon❗✔️:             seed = NEXT_SEED(seed, ap->color);
            seed = next_hash(seed, atom.color as u64);
            // Avalon❗✔️:             if (1)
            // Avalon❗✔️:             result += SetPathBitsRec(mp, nbp,
            // Avalon❗✔️:                                      fp_counts, ncounts,
            // Avalon❗✔️:                                      NEXT_SEED(seed,3*107), touched_indices,
            // Avalon❗✔️:                                      1,
            // Avalon❗✔️:                                      2, 4,  /* path length 2 to 4 */
            // Avalon❗✔️:                                      i, 0, -1,
            // Avalon❗✔️:                                      STOP_AT_HEAVY_ATOM |            // don't cross very heavy atoms
            // Avalon❗✔️:                                      IGNORE_TERM_SYMBOL |     /* gedeck.mol */
            // Avalon❗✔️:                                      PROCESS_CHAINS,
            // Avalon❗✔️:                                      exclude_atom,
            // Avalon❗✔️:                                      prefix);
            result += set_path_bits_rec(
                molecule,
                &state.neighbours,
                counts,
                next_hash(seed, 3 * 107),
                &mut touched,
                1,
                2,
                4,
                atom_index,
                0,
                -1,
                PathFlags::STOP_AT_HEAVY_ATOM | PathFlags::IGNORE_TERM_SYMBOL | PathFlags::PROCESS_CHAINS,
                exclude_atom,
            );
        }
        // Avalon❗✔️:          seed = old_seed;
        seed = parent_seed;
        // Avalon❗✔️:          /* Add another set of bits for paths starting at spiro atoms */
        // Avalon❗✔️: if (1) // [TODO] 6
        // Avalon❗✔️:          if (atom_status[i] >= 4 &&
        // Avalon❗✔️:              ap->color != 5      &&         // don't do it for boron
        // Avalon❗✔️:              ap->color < 18)                // don't do it for transition metals
        if state.atom_status[atom_index] >= 4 && atom.color != 5 && atom.color < 18 {
            // Avalon❗✔️: // fprintf(stderr, "spiro atom %s(%d) found\n", ap->atom_symbol, i+1);
            // Avalon❗✔️:             seed = NEXT_SEED(seed, ap->color+55);
            seed = next_hash(seed, (atom.color + 55) as u64);
            // Avalon❗✔️: // if (1)
            // Avalon❗✔️:             result += SetPathBitsRec(mp, nbp,
            // Avalon❗✔️:                                      fp_counts, ncounts,
            // Avalon❗✔️:                                      NEXT_SEED(seed,3*109), touched_indices,
            // Avalon❗✔️:                                      1,
            // Avalon❗✔️:                                      2, 4,  /* path length 1 to 4 */
            // Avalon❗✔️:                                      i, 0, -1,
            // Avalon❗✔️:                                      STOP_AT_HEAVY_ATOM |            // don't cross very heavy atoms
            // Avalon❗✔️:                                      IGNORE_TERM_SYMBOL |
            // Avalon❗✔️:                                      IGNORE_PATH_SYMBOL |
            // Avalon❗✔️:                                      PROCESS_RING_CLOSURES |
            // Avalon❗✔️:                                      PROCESS_CHAINS,
            // Avalon❗✔️:                                      exclude_atom,
            // Avalon❗✔️:                                      prefix);
            result += set_path_bits_rec(
                molecule,
                &state.neighbours,
                counts,
                next_hash(seed, 3 * 109),
                &mut touched,
                1,
                2,
                4,
                atom_index,
                0,
                -1,
                PathFlags::STOP_AT_HEAVY_ATOM
                    | PathFlags::IGNORE_TERM_SYMBOL
                    | PathFlags::IGNORE_PATH_SYMBOL
                    | PathFlags::PROCESS_RING_CLOSURES
                    | PathFlags::PROCESS_CHAINS,
                exclude_atom,
            );
        }
        // Avalon❗✔️:          seed = old_seed;
        // Avalon❗✔️:          /* Add bits for paths starting with rare bond orders */
        // Avalon❗✔️:          if (1)
        // Avalon❗✔️:          for (j=0; j<nbp[i].n_ligands; j++)
        // Avalon❗✔️:          {
        for (&neighbour_atom, &bond_index) in state.neighbours[atom_index]
            .atoms()
            .iter()
            .zip(state.neighbours[atom_index].bonds())
        {
            // Avalon❗✔️:             bp = &mp->bond_array[nbp[i].bonds[j]];
            // Avalon❗✔️:             ai = nbp[i].atoms[j];
            let bond = &molecule.bonds[bond_index];
            // Avalon❗✔️:             if (ai+1 == exclude_atom) continue;
            if neighbour_atom as i32 + 1 == exclude_atom {
                continue;
            }
            // Avalon❗✔️:             if (bp->color == 0) continue;
            // Avalon❗✔️:             if (bp->bond_type != DOUBLE  &&
            // Avalon❗✔️:                 bp->bond_type != TRIPLE)
            // Avalon❗✔️:                continue;
            if bond.color == 0 || !matches!(bond.bond_type, 2 | 3) {
                continue;
            }
            seed = next_hash(parent_seed, atom.color as u64);
            // Avalon❗✔️:             seed = NEXT_SEED(seed, bp->color*613);
            seed = next_hash(seed, (bond.color * 613) as u64);
            // Avalon❗✔️:             touched_indices[ai] = 1;    /* updating */
            touched[neighbour_atom] = 1;
            // Avalon❗✔️:             result += SetPathBitsRec(mp, nbp,
            // Avalon❗✔️:                                      fp_counts, ncounts,
            // Avalon❗✔️:                                      seed, touched_indices,
            // Avalon❗✔️:                                      2,
            // Avalon❗✔️:                                      5, 5,
            // Avalon❗✔️:                                      ai, 0, i,
            // Avalon❗✔️:                                      STOP_AT_HEAVY_ATOM |            // don't cross very heavy atoms
            // Avalon❗✔️:                                      IGNORE_PATH_SYMBOL |
            // Avalon❗✔️:                                      // IGNORE_TERM_SYMBOL |
            // Avalon❗✔️:                                      PROCESS_RING_CLOSURES |
            // Avalon❗✔️:                                      (0*PROCESS_CHAINS),
            // Avalon❗✔️:                                      exclude_atom,
            // Avalon❗✔️:                                      prefix);
            result += set_path_bits_rec(
                molecule,
                &state.neighbours,
                counts,
                seed,
                &mut touched,
                2,
                5,
                5,
                neighbour_atom,
                0,
                atom_index as i32,
                PathFlags::STOP_AT_HEAVY_ATOM | PathFlags::IGNORE_PATH_SYMBOL | PathFlags::PROCESS_RING_CLOSURES,
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

fn count_atom_types(
    molecule: &MoleculeState,
    state: &FingerprintPreprocessingState,
    counts: &mut [i32],
    exclude_atom: i32,
) -> i32 {
    // Avalon❗✔️: // fprintf(stderr, "USE_ATOM_COUNT\n");
    // Avalon❗✔️:    if (which_bits & USE_ATOM_COUNT)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       strcpy(prefix, "AC:");
    // Avalon❗✔️:       /* Collect hashed counts of atom types with hydrogen counts */
    // Avalon❗✔️:       for (i=0; i<NCOUNT_HASH; i++)
    // Avalon❗✔️:          atom_type_count_hash[i] = 0;
    let mut atom_type_count_hash = [0_i32; NCOUNT_HASH];
    // Avalon❗✔️:       for (i=0; i<NCOUNT_SEED_HASH; i++)
    // Avalon❗✔️:          atom_type_count_seed_hash[i] = 0;
    let mut atom_type_count_seed_hash = vec![0_i32; NCOUNT_SEED_HASH];
    // Avalon❗✔️:       nringch2=0; nfusionch=0; nspiro=0;
    let mut ring_ch2 = 0_i32;
    let mut fusion_ch = 0_i32;
    let mut spiro = 0_i32;
    // Avalon❗✔️:       for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    // Avalon❗✔️:       {
    for (atom_index, atom) in molecule.atoms.iter().enumerate() {
        // Avalon❗✔️:          if (i+1 == exclude_atom) continue;
        if atom_index as i32 + 1 == exclude_atom {
            continue;
        }
        // Avalon❗✔️:          if (ap->color == 0) continue;
        if atom.color == 0 {
            continue;
        }
        // Avalon❗✔️:          if (ap->color == 6)
        // Avalon❗✔️:          {
        if atom.color == 6 {
            // Avalon❗✔️:             /* count CH2 in ring */
            // Avalon❗✔️:             if (atom_status[i] > 0  &&  H_count[i+1] >= 2) nringch2++;
            if state.atom_status[atom_index] > 0 && state.hydrogen_counts[atom_index + 1] >= 2 {
                ring_ch2 += 1;
            }
            // Avalon❗✔️:             /* count ring fusion CH atoms */
            // Avalon❗✔️:             if (atom_status[i] > 2  &&  H_count[i+1] >= 1) nfusionch++;
            if state.atom_status[atom_index] > 2 && state.hydrogen_counts[atom_index + 1] >= 1 {
                fusion_ch += 1;
            }
            // Avalon❗✔️:             if (atom_status[i] > 3) nspiro++;
            if state.atom_status[atom_index] > 3 {
                spiro += 1;
            }
            // Avalon❗✔️:             continue;   // carbon has a special count processing
            continue;
            // Avalon❗✔️:          }
        }
        // Avalon❗✔️:          hash = 0;
        let mut hash = 0_i32;
        // Avalon❗✔️:          seed = NEXT_SEED(317*ATOM_COUNT_SEED,507);
        let mut seed = next_hash(317 * ATOM_COUNT_SEED, 507);
        // Avalon❗✔️:          seed = NEXT_SEED(seed,ap->color+17);
        seed = next_hash(seed, (atom.color + 17) as u64);
        // Avalon❗✔️:          atom_type_count_seed_hash[seed%NCOUNT_SEED_HASH]++;
        // `NCOUNT_SEED_HASH` expands from the unparenthesized C macro
        // `128*128`, so `%` and `*` associate left-to-right here.
        atom_type_count_seed_hash[(seed as usize % NCOUNT_HASH) * NCOUNT_HASH] += 1;
        // Avalon❗✔️:          /* normal hetero only with hydrogen */
        // Avalon❗✔️:          if ((ap->color == 7  ||  ap->color == 8)  &&  H_count[i+1] <= 0)
        // Avalon❗✔️:             continue;
        if matches!(atom.color, 7 | 8) && state.hydrogen_counts[atom_index + 1] <= 0 {
            continue;
        }
        // Avalon❗✔️:          hash = hash*7 + ap->color+13;
        hash = hash * 7 + atom.color + 13;
        // Avalon❗✔️:          atom_type_count_hash[hash%NCOUNT_HASH]++;
        atom_type_count_hash[hash as usize % NCOUNT_HASH] += 1;
        // Avalon❗✔️:          seed = NEXT_SEED(seed,ap->color+13);
        seed = next_hash(seed, (atom.color + 13) as u64);
        // Avalon❗✔️:          atom_type_count_seed_hash[seed%NCOUNT_SEED_HASH]++;
        atom_type_count_seed_hash[(seed as usize % NCOUNT_HASH) * NCOUNT_HASH] += 1;
        // Avalon❗✔️:          /* one more bit for rare types */
        // Avalon❗✔️:          if (ap->color != 7  &&  ap->color != 8)
        // Avalon❗✔️:          {
        if !matches!(atom.color, 7 | 8) {
            // Avalon❗✔️:             hash = hash*7 + ap->color+2*13;
            hash = hash * 7 + atom.color + 2 * 13;
            // Avalon❗✔️:             atom_type_count_hash[hash%NCOUNT_HASH]++;
            atom_type_count_hash[hash as usize % NCOUNT_HASH] += 1;
            // Avalon❗✔️:             seed = NEXT_SEED(seed,ap->color+2*13);
            seed = next_hash(seed, (atom.color + 2 * 13) as u64);
            // Avalon❗✔️:             atom_type_count_seed_hash[seed%NCOUNT_SEED_HASH]++;
            atom_type_count_seed_hash[(seed as usize % NCOUNT_HASH) * NCOUNT_HASH] += 1;
            // Avalon❗✔️:          }
        }
        // Avalon❗✔️:       }
    }
    // Avalon❗✔️:       nbits = 0;
    let mut bit_count = 0_i32;
    // Avalon❗✔️:       /* Now, we set the corresponding bits */
    // Avalon❗✔️:
    // Avalon❗✔️:       if (ncounts <= 2048)      // old atom count fingerprints
    if counts.len() <= 2048 {
        // Avalon❗✔️:           for (i=0; i<NCOUNT_HASH; i++)
        // Avalon❗✔️:           {
        for (index, &count) in atom_type_count_hash.iter().enumerate() {
            // Avalon❗✔️:              if (atom_type_count_hash[i] > 0)
            // Avalon❗✔️:              {
            if count > 0 {
                // Avalon❗✔️:                 ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*19,3), atom_type_count_hash[i]);
                add_bit_count(counts, next_hash((index * 19) as u64, 3), count);
                // Avalon❗✔️:                 nbits++;
                bit_count += 1;
                // Avalon❗✔️:              }
            }
            // Avalon❗✔️:              if (atom_type_count_hash[i] > 1)
            // Avalon❗✔️:              {
            if count > 1 {
                // Avalon❗✔️:                 ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*23,5), atom_type_count_hash[i]);
                add_bit_count(counts, next_hash((index * 23) as u64, 5), count);
                // Avalon❗✔️:                 nbits++;
                bit_count += 1;
                // Avalon❗✔️:              }
            }
            // Avalon❗✔️:              if (atom_type_count_hash[i] > 2)
            // Avalon❗✔️:              {
            if count > 2 {
                // Avalon❗✔️:                 ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*29,7), atom_type_count_hash[i]);
                add_bit_count(counts, next_hash((index * 29) as u64, 7), count);
                // Avalon❗✔️:                 nbits++;
                bit_count += 1;
                // Avalon❗✔️:              }
            }
            // Avalon❗✔️:              if (atom_type_count_hash[i] > 4)
            // Avalon❗✔️:              {
            if count > 4 {
                // Avalon❗✔️:                 ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*31,11), atom_type_count_hash[i]);
                add_bit_count(counts, next_hash((index * 31) as u64, 11), count);
                // Avalon❗✔️:                 nbits++;
                bit_count += 1;
                // Avalon❗✔️:              }
            }
            // Avalon❗✔️:              if (atom_type_count_hash[i] > 8)
            // Avalon❗✔️:              {
            if count > 8 {
                // Avalon❗✔️:                 ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*37,13), atom_type_count_hash[i]);
                add_bit_count(counts, next_hash((index * 37) as u64, 13), count);
                // Avalon❗✔️:                 nbits++;
                bit_count += 1;
                // Avalon❗✔️:              }
            }
            // Avalon❗✔️:              if (atom_type_count_hash[i] > 16)
            // Avalon❗✔️:              {
            if count > 16 {
                // Avalon❗✔️:                 ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*41,17), atom_type_count_hash[i]);
                add_bit_count(counts, next_hash((index * 41) as u64, 17), count);
                // Avalon❗✔️:                 nbits++;
                bit_count += 1;
                // Avalon❗✔️:              }
            }
        }
        // Avalon❗✔️:           }
        // Avalon❗✔️:       else      // new atom count fingerprints
    } else {
        // Avalon❗✔️:           for (i=0; i<NCOUNT_SEED_HASH; i++)
        // Avalon❗✔️:           {
        for (index, &count) in atom_type_count_seed_hash.iter().enumerate() {
            // Avalon❗✔️:              if (atom_type_count_seed_hash[i] > 0)
            // Avalon❗✔️:              {
            if count > 0 {
                // Avalon❗✔️:                 ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*19,3), atom_type_count_seed_hash[i]);
                add_bit_count(counts, next_hash((index * 19) as u64, 3), count);
                // Avalon❗✔️:                 nbits++;
                bit_count += 1;
                // Avalon❗✔️:              }
            }
            // Avalon❗✔️:              if (atom_type_count_seed_hash[i] > 1)
            // Avalon❗✔️:              {
            if count > 1 {
                // Avalon❗✔️:                 ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*23,5), atom_type_count_seed_hash[i]);
                add_bit_count(counts, next_hash((index * 23) as u64, 5), count);
                // Avalon❗✔️:                 nbits++;
                bit_count += 1;
                // Avalon❗✔️:              }
            }
            // Avalon❗✔️:              if (atom_type_count_seed_hash[i] > 2)
            // Avalon❗✔️:              {
            if count > 2 {
                // Avalon❗✔️:                 ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*29,7), atom_type_count_seed_hash[i]);
                add_bit_count(counts, next_hash((index * 29) as u64, 7), count);
                // Avalon❗✔️:                 nbits++;
                bit_count += 1;
                // Avalon❗✔️:              }
            }
            // Avalon❗✔️:              if (atom_type_count_seed_hash[i] > 3)
            // Avalon❗✔️:              {
            if count > 3 {
                // Avalon❗✔️:                 ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*29,71), atom_type_count_seed_hash[i]);
                add_bit_count(counts, next_hash((index * 29) as u64, 71), count);
                // Avalon❗✔️:                 nbits++;
                bit_count += 1;
                // Avalon❗✔️:              }
            }
            // Avalon❗✔️:              if (atom_type_count_seed_hash[i] > 4)
            // Avalon❗✔️:              {
            if count > 4 {
                // Avalon❗✔️:                 ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*31,11), atom_type_count_seed_hash[i]);
                add_bit_count(counts, next_hash((index * 31) as u64, 11), count);
                // Avalon❗✔️:                 nbits++;
                bit_count += 1;
                // Avalon❗✔️:              }
            }
            // Avalon❗✔️:              if (atom_type_count_seed_hash[i] > 6)
            // Avalon❗✔️:              {
            if count > 6 {
                // Avalon❗✔️:                 ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*31,113), atom_type_count_seed_hash[i]);
                add_bit_count(counts, next_hash((index * 31) as u64, 113), count);
                // Avalon❗✔️:                 nbits++;
                bit_count += 1;
                // Avalon❗✔️:              }
            }
            // Avalon❗✔️:              if (atom_type_count_seed_hash[i] > 8)
            // Avalon❗✔️:              {
            if count > 8 {
                // Avalon❗✔️:                 ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*37,13), atom_type_count_seed_hash[i]);
                add_bit_count(counts, next_hash((index * 37) as u64, 13), count);
                // Avalon❗✔️:                 nbits++;
                bit_count += 1;
                // Avalon❗✔️:              }
            }
            // Avalon❗✔️:              if (atom_type_count_seed_hash[i] > 12)
            // Avalon❗✔️:              {
            if count > 12 {
                // Avalon❗✔️:                 ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*37,133), atom_type_count_seed_hash[i]);
                add_bit_count(counts, next_hash((index * 37) as u64, 133), count);
                // Avalon❗✔️:                 nbits++;
                bit_count += 1;
                // Avalon❗✔️:              }
            }
            // Avalon❗✔️:              if (atom_type_count_seed_hash[i] > 16)
            // Avalon❗✔️:              {
            if count > 16 {
                // Avalon❗✔️:                 ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*41,17), atom_type_count_seed_hash[i]);
                add_bit_count(counts, next_hash((index * 41) as u64, 17), count);
                // Avalon❗✔️:                 nbits++;
                bit_count += 1;
                // Avalon❗✔️:              }
            }
        }
        // Avalon❗✔️:           }
    }

    // Avalon❗✔️:       if (naromatic >= 10) {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 341), naromatic); nbits++;}
    add_count_bit_if_at_least(counts, &mut bit_count, state.aromatic_bond_count, 10, 341);
    // Avalon❗✔️:       if (naromatic >= 14) {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 441), naromatic); nbits++;}
    add_count_bit_if_at_least(counts, &mut bit_count, state.aromatic_bond_count, 14, 441);
    // Avalon❗✔️:       if (naromatic >= 18) {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 541), naromatic); nbits++;}
    add_count_bit_if_at_least(counts, &mut bit_count, state.aromatic_bond_count, 18, 541);
    // Avalon❗✔️:       if (naromatic >= 22) {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 641), naromatic); nbits++;}
    add_count_bit_if_at_least(counts, &mut bit_count, state.aromatic_bond_count, 22, 641);
    // Avalon❗✔️:       if (naromatic >= 26) {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 741), naromatic); nbits++;}
    add_count_bit_if_at_least(counts, &mut bit_count, state.aromatic_bond_count, 26, 741);

    // Preserve the source's use of `naromatic`, rather than `ndouble`, as the
    // count added by these three double-bond threshold branches.
    // Avalon❗✔️:       if (ndouble >= 3) {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 371), naromatic); nbits++;}
    add_count_bit_if_condition(
        counts,
        &mut bit_count,
        state.double_bond_count >= 3,
        state.aromatic_bond_count,
        371,
    );
    // Avalon❗✔️:       if (ndouble >= 5) {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 471), naromatic); nbits++;}
    add_count_bit_if_condition(
        counts,
        &mut bit_count,
        state.double_bond_count >= 5,
        state.aromatic_bond_count,
        471,
    );
    // Avalon❗✔️:       if (ndouble >= 8) {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 571), naromatic); nbits++;}
    add_count_bit_if_condition(
        counts,
        &mut bit_count,
        state.double_bond_count >= 8,
        state.aromatic_bond_count,
        571,
    );

    // Avalon❗✔️:       if (nringch2 >= 6)  {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 411), nringch2); nbits++;}
    add_count_bit_if_at_least(counts, &mut bit_count, ring_ch2, 6, 411);
    // Avalon❗✔️:       if (nringch2 >= 12) {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 511), nringch2); nbits++;}
    add_count_bit_if_at_least(counts, &mut bit_count, ring_ch2, 12, 511);
    // Avalon❗✔️:       if (nringch2 >= 22) {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 611), nringch2); nbits++;}
    add_count_bit_if_at_least(counts, &mut bit_count, ring_ch2, 22, 611);

    // Avalon❗✔️:       if (nfusionch >= 2)  {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 421), nfusionch); nbits++;}
    add_count_bit_if_at_least(counts, &mut bit_count, fusion_ch, 2, 421);
    // Avalon❗✔️:       if (nfusionch >= 4)  {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 521), nfusionch); nbits++;}
    add_count_bit_if_at_least(counts, &mut bit_count, fusion_ch, 4, 521);

    // Avalon❗✔️:       // Note spiro atoms
    // Avalon❗✔️:       if (nspiro >= 1)  {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 322), nspiro); nbits++;}
    add_count_bit_if_at_least(counts, &mut bit_count, spiro, 1, 322);
    // Avalon❗✔️:       if (nspiro >= 1)  {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 422), nspiro); nbits++;}
    add_count_bit_if_at_least(counts, &mut bit_count, spiro, 1, 422);
    // Avalon❗✔️:       if (nspiro >= 2)  {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 522), nspiro); nbits++;}
    add_count_bit_if_at_least(counts, &mut bit_count, spiro, 2, 522);
    // Avalon❗✔️:       if (nspiro >= 2)  {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 622), nspiro); nbits++;}
    add_count_bit_if_at_least(counts, &mut bit_count, spiro, 2, 622);

    // Avalon❗✔️:       if (nfusionb >= 1)  {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 425), nfusionb); nbits++;}
    add_count_bit_if_at_least(counts, &mut bit_count, state.fusion_bond_count, 1, 425);
    // Avalon❗✔️:       if (nfusionb >= 2)  {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 525), nfusionb); nbits++;}
    add_count_bit_if_at_least(counts, &mut bit_count, state.fusion_bond_count, 2, 525);
    // Avalon❗✔️:       if (nfusionb >= 3)  {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 625), nfusionb); nbits++;}
    add_count_bit_if_at_least(counts, &mut bit_count, state.fusion_bond_count, 3, 625);
    // Avalon❗✔️:       if (nfusionb >= 5)  {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 825), nfusionb); nbits++;}
    add_count_bit_if_at_least(counts, &mut bit_count, state.fusion_bond_count, 5, 825);
    // Avalon❗✔️:       if (nrare_atoms > 0)
    // Avalon❗✔️:       {
    if state.rare_atom_count > 0 {
        // Avalon❗✔️:          seed = NEXT_SEED(317*ATOM_COUNT_SEED,5);
        let seed = next_hash(317 * ATOM_COUNT_SEED, 5);
        // Avalon❗✔️:          ADD_BIT_COUNT(fp_counts, ncounts, seed, nrare_atoms); nbits++;
        add_bit_count(counts, seed, state.rare_atom_count);
        bit_count += 1;
        // Avalon❗✔️:          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(seed,101), nrare_atoms); nbits++;
        add_bit_count(counts, next_hash(seed, 101), state.rare_atom_count);
        bit_count += 1;
        // Avalon❗✔️:          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(seed,103), nrare_atoms); nbits++;
        add_bit_count(counts, next_hash(seed, 103), state.rare_atom_count);
        bit_count += 1;
        // Avalon❗✔️:          if (nrare_atoms > 1)
        if state.rare_atom_count > 1 {
            // Avalon❗✔️:             ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(seed,203), nrare_atoms);
            add_bit_count(counts, next_hash(seed, 203), state.rare_atom_count);
            // Avalon❗✔️:             nbits++;
            bit_count += 1;
        }
        // Avalon❗✔️:          if (nrare_atoms > 3)
        if state.rare_atom_count > 3 {
            // Avalon❗✔️:             ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(seed,211), nrare_atoms);
            add_bit_count(counts, next_hash(seed, 211), state.rare_atom_count);
            // Avalon❗✔️:             nbits++;
            bit_count += 1;
        }
        // Avalon❗✔️:       }
    }
    // Avalon❗✔️:       result += nbits;
    // Avalon❗✔️:    }
    bit_count
}

#[inline]
fn add_count_bit_if_at_least(counts: &mut [i32], bit_count: &mut i32, value: i32, threshold: i32, increment: u64) {
    add_count_bit_if_condition(counts, bit_count, value >= threshold, value, increment);
}

#[inline]
fn add_count_bit_if_condition(counts: &mut [i32], bit_count: &mut i32, condition: bool, value: i32, increment: u64) {
    if condition {
        add_bit_count(counts, next_hash(ATOM_COUNT_SEED, increment), value);
        *bit_count += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::properties::avalon_fingerprint::reaccs::{Atom, Bond};

    fn low_flag_fixture() -> MoleculeState {
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

    fn sparse(counts: &[i32]) -> Vec<(usize, i32)> {
        counts
            .iter()
            .copied()
            .enumerate()
            .filter(|&(_, count)| count != 0)
            .collect()
    }

    fn sparse_text(counts: &[i32]) -> String {
        sparse(counts)
            .into_iter()
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

    struct LowFlagGolden {
        mask: u32,
        basic_result: i32,
        basic_counts: &'static str,
        basic_bits: &'static str,
        daylight_result: i32,
        daylight_counts: &'static str,
        daylight_bits: &'static str,
    }

    // Generated with AvalonToolkit_2.0.5-pre.3 `SetFingerprintCountsWithFocus`
    // and `SetFingerprintBits` from `tmp/parity-audit/avalon_ring_state.c`.
    const LOW_FLAG_GOLDENS: &[LowFlagGolden] = &[
        LowFlagGolden {
            mask: 0x01,
            basic_result: 28,
            basic_counts: "4:6,15:6,17:4,40:6,42:2,62:4,",
            basic_bits: "1080020000050040",
            daylight_result: 28,
            daylight_counts: "4:6,15:6,17:4,40:6,42:2,62:4,",
            daylight_bits: "1080020000050040",
        },
        LowFlagGolden {
            mask: 0x02,
            basic_result: 12,
            basic_counts: "2:1,7:1,17:1,25:1,27:1,28:1,29:1,30:1,31:1,32:1,33:1,37:1,38:1,39:1,41:1,44:1,46:1,47:2,50:1,53:1,54:1,57:2,",
            basic_bits: "840002fae3d26402",
            daylight_result: 12,
            daylight_counts: "1:1,5:2,9:1,10:1,11:1,14:3,17:1,18:1,22:1,26:2,30:1,32:1,33:1,42:1,43:1,44:1,49:1,52:1,58:1,59:1,",
            daylight_bits: "224e4644031c120c",
        },
        LowFlagGolden {
            mask: 0x04,
            basic_result: 28,
            basic_counts: "3:1,4:1,8:2,12:2,13:1,18:1,22:1,25:1,26:3,29:1,31:2,33:1,38:1,39:1,41:1,43:1,44:1,45:1,51:2,54:1,55:1,56:1,58:1,59:1,62:1,63:1,",
            basic_bits: "183144a6c23ac8cd",
            daylight_result: 38,
            daylight_counts: "1:1,2:1,3:1,4:3,8:1,12:1,13:3,14:2,15:3,16:2,18:3,21:1,22:3,23:2,24:1,25:1,27:1,33:1,34:1,35:2,36:1,38:1,39:1,40:1,41:1,47:3,51:1,52:1,54:1,59:5,61:1,63:1,",
            daylight_bits: "1ef1e50bde8358a8",
        },
        LowFlagGolden {
            mask: 0x08,
            basic_result: 41,
            basic_counts: "4:4,9:2,17:2,18:2,19:2,21:2,25:1,27:4,30:2,34:4,35:6,40:4,45:4,46:2,48:2,52:4,56:2,63:4,",
            basic_bits: "10022e4a0c611181",
            daylight_result: 41,
            daylight_counts: "2:1,4:5,8:1,9:1,10:1,13:3,17:1,20:1,23:1,24:2,25:1,32:1,35:4,36:2,38:1,40:4,45:4,48:3,49:1,52:5,57:2,60:2,62:1,63:5,",
            daylight_bits: "14279203592113d2",
        },
        LowFlagGolden {
            mask: 0x10,
            basic_result: 1,
            basic_counts: "35:1,",
            basic_bits: "0000000008000000",
            daylight_result: 2,
            daylight_counts: "35:1,",
            daylight_bits: "0000000008000000",
        },
        LowFlagGolden {
            mask: 0x03,
            basic_result: 40,
            basic_counts: "2:1,4:6,7:1,15:6,17:5,25:1,27:1,28:1,29:1,30:1,31:1,32:1,33:1,37:1,38:1,39:1,40:6,41:1,42:2,44:1,46:1,47:2,50:1,53:1,54:1,57:2,62:4,",
            basic_bits: "948002fae3d76442",
            daylight_result: 40,
            daylight_counts: "1:1,4:6,5:2,9:1,10:1,11:1,14:3,15:6,17:5,18:1,22:1,26:2,30:1,32:1,33:1,40:6,42:3,43:1,44:1,49:1,52:1,58:1,59:1,62:4,",
            daylight_bits: "32ce4644031d124c",
        },
        LowFlagGolden {
            mask: 0x05,
            basic_result: 56,
            basic_counts: "3:1,4:7,8:2,12:2,13:1,15:6,17:4,18:1,22:1,25:1,26:3,29:1,31:2,33:1,38:1,39:1,40:6,41:1,42:2,43:1,44:1,45:1,51:2,54:1,55:1,56:1,58:1,59:1,62:5,63:1,",
            basic_bits: "18b146a6c23fc8cd",
            daylight_result: 66,
            daylight_counts: "1:1,2:1,3:1,4:9,8:1,12:1,13:3,14:2,15:9,16:2,17:4,18:3,21:1,22:3,23:2,24:1,25:1,27:1,33:1,34:1,35:2,36:1,38:1,39:1,40:7,41:1,42:2,47:3,51:1,52:1,54:1,59:5,61:1,62:4,63:1,",
            daylight_bits: "1ef1e70bde8758e8",
        },
        LowFlagGolden {
            mask: 0x09,
            basic_result: 69,
            basic_counts: "4:10,9:2,15:6,17:6,18:2,19:2,21:2,25:1,27:4,30:2,34:4,35:6,40:10,42:2,45:4,46:2,48:2,52:4,56:2,62:4,63:4,",
            basic_bits: "10822e4a0c6511c1",
            daylight_result: 69,
            daylight_counts: "2:1,4:11,8:1,9:1,10:1,13:3,15:6,17:5,20:1,23:1,24:2,25:1,32:1,35:4,36:2,38:1,40:10,42:2,45:4,48:3,49:1,52:5,57:2,60:2,62:5,63:5,",
            daylight_bits: "14a79203592513d2",
        },
        LowFlagGolden {
            mask: 0x11,
            basic_result: 29,
            basic_counts: "4:6,15:6,17:4,35:1,40:6,42:2,62:4,",
            basic_bits: "1080020008050040",
            daylight_result: 30,
            daylight_counts: "4:6,15:6,17:4,35:1,40:6,42:2,62:4,",
            daylight_bits: "1080020008050040",
        },
        LowFlagGolden {
            mask: 0x06,
            basic_result: 40,
            basic_counts: "2:1,3:1,4:1,7:1,8:2,12:2,13:1,17:1,18:1,22:1,25:2,26:3,27:1,28:1,29:2,30:1,31:3,32:1,33:2,37:1,38:2,39:2,41:2,43:1,44:2,45:1,46:1,47:2,50:1,51:2,53:1,54:2,55:1,56:1,57:2,58:1,59:1,62:1,63:1,",
            basic_bits: "9c3146fee3faeccf",
            daylight_result: 50,
            daylight_counts: "1:2,2:1,3:1,4:3,5:2,8:1,9:1,10:1,11:1,12:1,13:3,14:5,15:3,16:2,17:1,18:4,21:1,22:4,23:2,24:1,25:1,26:2,27:1,30:1,32:1,33:2,34:1,35:2,36:1,38:1,39:1,40:1,41:1,42:1,43:1,44:1,47:3,49:1,51:1,52:2,54:1,58:1,59:6,61:1,63:1,",
            daylight_bits: "3effe74fdf9f5aac",
        },
        LowFlagGolden {
            mask: 0x0a,
            basic_result: 53,
            basic_counts: "2:1,4:4,7:1,9:2,17:3,18:2,19:2,21:2,25:2,27:5,28:1,29:1,30:3,31:1,32:1,33:1,34:4,35:6,37:1,38:1,39:1,40:4,41:1,44:1,45:4,46:3,47:2,48:2,50:1,52:4,53:1,54:1,56:2,57:2,63:4,",
            basic_bits: "94022efaeff37583",
            daylight_result: 53,
            daylight_counts: "1:1,2:1,4:5,5:2,8:1,9:2,10:2,11:1,13:3,14:3,17:2,18:1,20:1,22:1,23:1,24:2,25:1,26:2,30:1,32:2,33:1,35:4,36:2,38:1,40:4,42:1,43:1,44:1,45:4,48:3,49:2,52:6,57:2,58:1,59:1,60:2,62:1,63:5,",
            daylight_bits: "366fd6475b3d13de",
        },
        LowFlagGolden {
            mask: 0x12,
            basic_result: 13,
            basic_counts: "2:1,7:1,17:1,25:1,27:1,28:1,29:1,30:1,31:1,32:1,33:1,35:1,37:1,38:1,39:1,41:1,44:1,46:1,47:2,50:1,53:1,54:1,57:2,",
            basic_bits: "840002faebd26402",
            daylight_result: 14,
            daylight_counts: "1:1,5:2,9:1,10:1,11:1,14:3,17:1,18:1,22:1,26:2,30:1,32:1,33:1,35:1,42:1,43:1,44:1,49:1,52:1,58:1,59:1,",
            daylight_bits: "224e46440b1c120c",
        },
        LowFlagGolden {
            mask: 0x0c,
            basic_result: 69,
            basic_counts: "3:1,4:5,8:2,9:2,12:2,13:1,17:2,18:3,19:2,21:2,22:1,25:2,26:3,27:4,29:1,30:2,31:2,33:1,34:4,35:6,38:1,39:1,40:4,41:1,43:1,44:1,45:5,46:2,48:2,51:2,52:4,54:1,55:1,56:3,58:1,59:1,62:1,63:5,",
            basic_bits: "18336eeece7bd9cd",
            daylight_result: 79,
            daylight_counts: "1:1,2:2,3:1,4:8,8:2,9:1,10:1,12:1,13:6,14:2,15:3,16:2,17:1,18:3,20:1,21:1,22:3,23:3,24:3,25:2,27:1,32:1,33:1,34:1,35:6,36:3,38:2,39:1,40:5,41:1,45:4,47:3,48:3,49:1,51:1,52:6,54:1,57:2,59:5,60:2,61:1,62:1,63:6,",
            daylight_bits: "1ef7f70bdfa35bfa",
        },
        LowFlagGolden {
            mask: 0x14,
            basic_result: 29,
            basic_counts: "3:1,4:1,8:2,12:2,13:1,18:1,22:1,25:1,26:3,29:1,31:2,33:1,35:1,38:1,39:1,41:1,43:1,44:1,45:1,51:2,54:1,55:1,56:1,58:1,59:1,62:1,63:1,",
            basic_bits: "183144a6ca3ac8cd",
            daylight_result: 40,
            daylight_counts: "1:1,2:1,3:1,4:3,8:1,12:1,13:3,14:2,15:3,16:2,18:3,21:1,22:3,23:2,24:1,25:1,27:1,33:1,34:1,35:3,36:1,38:1,39:1,40:1,41:1,47:3,51:1,52:1,54:1,59:5,61:1,63:1,",
            daylight_bits: "1ef1e50bde8358a8",
        },
        LowFlagGolden {
            mask: 0x18,
            basic_result: 42,
            basic_counts: "4:4,9:2,17:2,18:2,19:2,21:2,25:1,27:4,30:2,34:4,35:7,40:4,45:4,46:2,48:2,52:4,56:2,63:4,",
            basic_bits: "10022e4a0c611181",
            daylight_result: 43,
            daylight_counts: "2:1,4:5,8:1,9:1,10:1,13:3,17:1,20:1,23:1,24:2,25:1,32:1,35:5,36:2,38:1,40:4,45:4,48:3,49:1,52:5,57:2,60:2,62:1,63:5,",
            daylight_bits: "14279203592113d2",
        },
    ];

    fn query_sensitive_fixture() -> MoleculeState {
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

    #[test]
    fn low_flag_single_and_pairwise_matrix_matches_native_counts_and_bits() {
        assert_eq!(LOW_FLAG_GOLDENS.len(), 15);
        assert_eq!(
            LOW_FLAG_GOLDENS
                .iter()
                .filter(|golden| golden.mask.count_ones() == 1)
                .count(),
            5
        );
        assert_eq!(
            LOW_FLAG_GOLDENS
                .iter()
                .filter(|golden| golden.mask.count_ones() == 2)
                .count(),
            10
        );

        for golden in LOW_FLAG_GOLDENS {
            let flags = AvalonFingerprintFlags::from_bits(golden.mask).unwrap();
            for as_query in [false, true] {
                for (fpflags, expected_result, expected_counts, expected_bits) in [
                    (0, golden.basic_result, golden.basic_counts, golden.basic_bits),
                    (
                        super::super::fingerprint_state::USE_DY_AROMATICITY,
                        golden.daylight_result,
                        golden.daylight_counts,
                        golden.daylight_bits,
                    ),
                ] {
                    let mut molecule = low_flag_fixture();
                    let mut counts = vec![0; 64];
                    let result =
                        count_low_flag_families(&mut molecule, &mut counts, flags, as_query, fpflags, 0).unwrap();
                    let case = format!(
                        "mask={:#04x}, as_query={as_query}, daylight={}",
                        golden.mask,
                        fpflags != 0
                    );

                    assert_eq!(result, expected_result, "{case}");
                    assert_eq!(sparse_text(&counts), expected_counts, "{case}");
                    assert_eq!(packed_bit_hex(&counts), expected_bits, "{case}");
                }
            }
        }
    }

    #[test]
    fn query_mode_low_flags_match_native_query_hydrogen_branches() {
        let cases = [
            (0x04, false, 4, "0:1,27:1,49:1,57:1,", "0100000800000202"),
            (0x04, true, 3, "0:1,27:1,57:1,", "0100000800000002"),
            (0x10, false, 1, "61:1,", "0000000000000020"),
            (0x10, true, 1, "61:1,", "0000000000000020"),
            (0x14, false, 5, "0:1,27:1,49:1,57:1,61:1,", "0100000800000222"),
            (0x14, true, 4, "0:1,27:1,57:1,61:1,", "0100000800000022"),
        ];

        for fpflags in [0, super::super::fingerprint_state::USE_DY_AROMATICITY] {
            for (mask, as_query, expected_result, expected_counts, expected_bits) in cases {
                let mut molecule = query_sensitive_fixture();
                let mut counts = vec![0; 64];
                let result = count_low_flag_families(
                    &mut molecule,
                    &mut counts,
                    AvalonFingerprintFlags::from_bits(mask).unwrap(),
                    as_query,
                    fpflags,
                    0,
                )
                .unwrap();
                let case = format!("mask={mask:#04x}, as_query={as_query}, daylight={}", fpflags != 0);

                assert_eq!(result, expected_result, "{case}");
                assert_eq!(sparse_text(&counts), expected_counts, "{case}");
                assert_eq!(packed_bit_hex(&counts), expected_bits, "{case}");
            }
        }
    }

    #[test]
    fn low_flag_families_match_native_exported_counts() {
        let fixtures = [
            (
                AvalonFingerprintFlags::RING_PATTERN,
                28,
                vec![(4, 6), (15, 6), (17, 4), (40, 6), (42, 2), (62, 4)],
            ),
            (
                AvalonFingerprintFlags::RING_PATH,
                12,
                vec![
                    (2, 1),
                    (7, 1),
                    (17, 1),
                    (25, 1),
                    (27, 1),
                    (28, 1),
                    (29, 1),
                    (30, 1),
                    (31, 1),
                    (32, 1),
                    (33, 1),
                    (37, 1),
                    (38, 1),
                    (39, 1),
                    (41, 1),
                    (44, 1),
                    (46, 1),
                    (47, 2),
                    (50, 1),
                    (53, 1),
                    (54, 1),
                    (57, 2),
                ],
            ),
            (
                AvalonFingerprintFlags::ATOM_SYMBOL_PATH,
                28,
                vec![
                    (3, 1),
                    (4, 1),
                    (8, 2),
                    (12, 2),
                    (13, 1),
                    (18, 1),
                    (22, 1),
                    (25, 1),
                    (26, 3),
                    (29, 1),
                    (31, 2),
                    (33, 1),
                    (38, 1),
                    (39, 1),
                    (41, 1),
                    (43, 1),
                    (44, 1),
                    (45, 1),
                    (51, 2),
                    (54, 1),
                    (55, 1),
                    (56, 1),
                    (58, 1),
                    (59, 1),
                    (62, 1),
                    (63, 1),
                ],
            ),
            (
                AvalonFingerprintFlags::ATOM_CLASS_PATH,
                41,
                vec![
                    (4, 4),
                    (9, 2),
                    (17, 2),
                    (18, 2),
                    (19, 2),
                    (21, 2),
                    (25, 1),
                    (27, 4),
                    (30, 2),
                    (34, 4),
                    (35, 6),
                    (40, 4),
                    (45, 4),
                    (46, 2),
                    (48, 2),
                    (52, 4),
                    (56, 2),
                    (63, 4),
                ],
            ),
            (AvalonFingerprintFlags::ATOM_COUNT, 1, vec![(35, 1)]),
        ];

        for (flag, expected_result, expected_counts) in fixtures {
            let mut molecule = low_flag_fixture();
            let mut counts = vec![0; 64];
            let result = count_low_flag_families(&mut molecule, &mut counts, flag, false, 0, 0).unwrap();
            assert_eq!(result, expected_result, "flag {:#x}", flag.bits());
            assert_eq!(sparse(&counts), expected_counts, "flag {:#x}", flag.bits());
        }
    }

    #[test]
    fn atom_count_large_slot_branch_matches_native_exported_counts() {
        let mut molecule = low_flag_fixture();
        let mut counts = vec![0; 4096];

        let result = count_low_flag_families(
            &mut molecule,
            &mut counts,
            AvalonFingerprintFlags::ATOM_COUNT,
            false,
            0,
            0,
        )
        .unwrap();

        assert_eq!(result, 3);
        assert_eq!(sparse(&counts), vec![(53, 1), (196, 1), (4067, 1)]);
    }
}
