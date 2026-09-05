//! Source-order implementation of Avalon non-SSS fingerprint flags `0x100000..=0x800000`.

use crate::FingerprintError;

use super::AvalonFingerprintFlags;
use super::fingerprint_state::{FingerprintPreprocessingState, with_prepared_fingerprint_state};
use super::hash::next_hash;
use super::reaccs::MoleculeState;
use super::symbols::atom_symbol_match;
use super::traversal::{NON_SSS_SEED, add_bit, build_path_length_matrix};

const ACCUMULATE_BITS: i32 = 0x0002;

pub(super) fn count_non_sss_flag_families(
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
            Ok(count_non_sss_flag_families_prepared(
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

pub(super) fn count_non_sss_flag_families_prepared(
    working: &MoleculeState,
    state: &FingerprintPreprocessingState,
    counts: &mut [i32],
    which_bits: AvalonFingerprintFlags,
    as_query: bool,
    exclude_atom: i32,
) -> i32 {
    let length_matrix = build_path_length_matrix(working, &state.neighbours, 12, exclude_atom);
    count_non_sss_bits(
        working,
        state,
        counts,
        which_bits,
        as_query,
        &length_matrix,
        exclude_atom,
    )
}

fn count_non_sss_bits(
    molecule: &MoleculeState,
    state: &FingerprintPreprocessingState,
    counts: &mut [i32],
    which_bits: AvalonFingerprintFlags,
    as_query: bool,
    length_matrix: &[Vec<i32>],
    exclude_atom: i32,
) -> i32 {
    // Avalon❗✔️:    /* Collect bits that describe scaffolds. Those bits cannot (yet?) be used for SSS screening */
    // Avalon❗✔️:    /* The method first collects a variant of the extended connectivity but derived only        */
    // Avalon❗✔️:    /* on the ring bond graph. The bits are then set by combining atom features with the        */
    // Avalon❗✔️:    /* corresponding scaffold extended connectivity of that atom. In addition, we add           */
    // Avalon❗✔️:    /* bits to code connectivity of scaffolds by certain numbers of bonds using the             */
    // Avalon❗✔️:    /* previously collected length_matrix. */
    // Avalon❗✔️:    if (!as_query  &&  (which_bits & (USE_NON_SSS_BITS)))
    // Avalon❗✔️:    {
    if as_query || which_bits.bits() & AvalonFingerprintFlags::NON_SSS_BITS.bits() == 0 {
        return 0;
    }
    // Avalon❗✔️:       strcpy(prefix, "NSB:");
    // Avalon❗✔️:       seed = NON_SSS_SEED;
    let seed = NON_SSS_SEED;
    // Avalon❗✔️:       extcon = TypeAlloc(mp->n_atoms, int);
    let mut extcon = vec![0_i32; molecule.atoms.len()];
    // Avalon❗✔️:       extcon2 = TypeAlloc(mp->n_atoms, int);
    let mut extcon2 = vec![0_i32; molecule.atoms.len()];
    // Avalon❗✔️:       /* initialized extended connectivity */
    // Avalon❗✔️:       for (j=0, ap=mp->atom_array; j<mp->n_atoms; j++, ap++)
    // Avalon❗✔️:       {
    for (atom_index, atom) in molecule.atoms.iter().enumerate() {
        // Avalon❗✔️:           if (atom_status[j] <= 0) continue;
        if state.atom_status[atom_index] <= 0 {
            continue;
        }
        // Avalon❗✔️:           extcon[j] = ap->rsize_flags;
        extcon[atom_index] = atom.rsize_flags as i32;
        // Avalon❗✔️:       }
    }
    // Avalon❗✔️:       /* propagate extended connectivity to neighbours for a few cycles */
    // Avalon❗✔️:       for (i=0; i<32; i++)
    // Avalon❗✔️:       {
    propagate_ring_connectivity(state, &mut extcon, &mut extcon2);

    // Avalon❗✔️:       /* propagate smallest hash to all members of ring system */
    propagate_smallest_ring_hash(state, &mut extcon);

    let mut result = 0_i32;
    // Avalon❗✔️:       /* Now, use extcon to set bits */
    // Avalon❗✔️:       for (j=0, ap=mp->atom_array; j<mp->n_atoms; j++, ap++)
    // Avalon❗✔️:       {
    for atom_index in 0..molecule.atoms.len() {
        // Avalon❗✔️:           if (j+1 == exclude_atom) continue;
        if atom_index as i32 + 1 == exclude_atom {
            continue;
        }
        // Avalon❗✔️:           /* skip non-ring atoms */
        // Avalon❗✔️:           if (atom_status[j] <= 0) continue;
        if state.atom_status[atom_index] <= 0 {
            continue;
        }
        // Avalon❗✔️:           if (which_bits & USE_SCAFFOLD_IDS)
        // Avalon❗✔️:           {
        if which_bits.contains(AvalonFingerprintFlags::SCAFFOLD_IDS) {
            // Avalon❗✔️:                strcpy(prefix, "SI:");
            // Avalon❗✔️: // ADD_BIT(fp_counts, ncounts, extcon[j]*1013 + seed);
            // Avalon❗✔️:               ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, extcon[j]*10013));
            add_bit(counts, next_hash(seed, extcon[atom_index].wrapping_mul(10013) as u64));
            // Avalon❗✔️: // ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, extcon[j]*40013));
            // Avalon❗✔️: 	      result += 1;
            result += 1;
            // Avalon❗✔️:           }
        }

        // Avalon❗✔️:           if (which_bits & USE_SCAFFOLD_LINKS)
        // Avalon❗✔️:           {
        if which_bits.contains(AvalonFingerprintFlags::SCAFFOLD_LINKS) {
            // Avalon❗✔️:                strcpy(prefix, "SL:");
            // Avalon❗✔️:               /* needs to be a substitution point */
            // Avalon❗✔️:               if (degree[j] <= atom_status[j]) continue;
            if state.degree[atom_index] <= state.atom_status[atom_index] {
                continue;
            }
            // Avalon❗✔️:               for (jj=0; jj<mp->n_atoms; jj++)
            // Avalon❗✔️:               {
            for other_index in 0..molecule.atoms.len() {
                // Avalon❗✔️:                   if (jj+1 == exclude_atom) continue;
                if other_index as i32 + 1 == exclude_atom {
                    continue;
                }
                // Avalon❗✔️:                   /* skip non-ring atoms */
                // Avalon❗✔️:                   if (atom_status[jj] <= 0) continue;
                if state.atom_status[other_index] <= 0 {
                    continue;
                }
                // Avalon❗✔️:                   if (j == jj) continue;
                if atom_index == other_index {
                    continue;
                }
                // Avalon❗✔️:                   /* needs to be a substitution point */
                // Avalon❗✔️:                   if (degree[jj] <= atom_status[jj]) continue;
                if state.degree[other_index] <= state.atom_status[other_index] {
                    continue;
                }
                // Avalon❗✔️:                   for (k = 0; k<12; k++)
                // Avalon❗✔️:                       if (length_matrix[j][jj] == (1<<k))
                // Avalon❗✔️:                           break;
                let mut path_length = 0_i32;
                while path_length < 12 && length_matrix[atom_index][other_index] != 1_i32 << path_length {
                    path_length += 1;
                }
                // Avalon❗✔️:                   if (k >= 12) continue;    // multiple paths => no chain connection
                if path_length >= 12 {
                    continue;
                }
                // Avalon❗✔️:                   if (k >= 3) continue;    // only short one counts
                if path_length >= 3 {
                    continue;
                }
                // Avalon❗✔️:                   // bit for ring system
                // Avalon❗✔️:                   ADD_BIT(fp_counts, ncounts, NEXT_SEED(NEXT_SEED(seed,3*k), extcon[j]*1013 + extcon[jj]*2003));
                let ring_pair = extcon[atom_index]
                    .wrapping_mul(1013)
                    .wrapping_add(extcon[other_index].wrapping_mul(2003));
                add_bit(
                    counts,
                    next_hash(next_hash(seed, path_length.wrapping_mul(3) as u64), ring_pair as u64),
                );
                // Avalon❗✔️:                   // bit for position
                // Avalon❗✔️:                   ADD_BIT(fp_counts, ncounts, NEXT_SEED(NEXT_SEED(seed,5*k), extcon2[j]*2013 + extcon[jj]*1003));
                let first_position = extcon2[atom_index]
                    .wrapping_mul(2013)
                    .wrapping_add(extcon[other_index].wrapping_mul(1003));
                add_bit(
                    counts,
                    next_hash(
                        next_hash(seed, path_length.wrapping_mul(5) as u64),
                        first_position as u64,
                    ),
                );
                // Avalon❗✔️:                   ADD_BIT(fp_counts, ncounts, NEXT_SEED(NEXT_SEED(seed,5*k), extcon[j]*2013 + extcon2[jj]*1003));
                let second_position = extcon[atom_index]
                    .wrapping_mul(2013)
                    .wrapping_add(extcon2[other_index].wrapping_mul(1003));
                add_bit(
                    counts,
                    next_hash(
                        next_hash(seed, path_length.wrapping_mul(5) as u64),
                        second_position as u64,
                    ),
                );
                // Avalon❗✔️:                   ADD_BIT(fp_counts, ncounts, NEXT_SEED(NEXT_SEED(seed,7*k), extcon2[j]*3013 + extcon2[jj]*3003));
                let positions = extcon2[atom_index]
                    .wrapping_mul(3013)
                    .wrapping_add(extcon2[other_index].wrapping_mul(3003));
                add_bit(
                    counts,
                    next_hash(next_hash(seed, path_length.wrapping_mul(7) as u64), positions as u64),
                );
                // Avalon❗✔️: 		  result += 4;
                result += 4;
                // Avalon❗✔️:               }
            }
            // Avalon❗✔️:           }
        }
        // Avalon❗✔️:       }
    }
    // Avalon❗✔️: // fprintf(stderr, "USE_SCAFFOLD_COLORS\n");
    // Avalon❗✔️:       if (which_bits & USE_SCAFFOLD_COLORS)
    // Avalon❗✔️:       {
    if which_bits.contains(AvalonFingerprintFlags::SCAFFOLD_COLORS) {
        // Avalon❗✔️:          strcpy(prefix, "SC:");
        // Avalon❗✔️:           for (j=0, ap=mp->atom_array; j<mp->n_atoms; j++, ap++)
        // Avalon❗✔️:           {
        for (atom_index, atom) in molecule.atoms.iter().enumerate() {
            // Avalon❗✔️:               if (j+1 == exclude_atom) continue;
            if atom_index as i32 + 1 == exclude_atom {
                continue;
            }
            // Avalon❗✔️:               if (atom_status[j] <= 0) continue;
            if state.atom_status[atom_index] <= 0 {
                continue;
            }
            // Avalon❗✔️:               tmp1 = 0;
            let mut offset = 0_i32;
            // Avalon❗✔️:               if (0 == strcmp(ap->atom_symbol, "C")) tmp1 = 101;
            if atom.atom_symbol == "C" {
                offset = 101;
            // Avalon❗✔️:               else if (0 == strcmp(ap->atom_symbol, "O")) tmp1 = 301;
            } else if atom.atom_symbol == "O" {
                offset = 301;
            // Avalon❗✔️:               else if (0 == strcmp(ap->atom_symbol, "N")) tmp1 = 401;
            } else if atom.atom_symbol == "N" {
                offset = 401;
            // Avalon❗✔️:               else if (0 == strcmp(ap->atom_symbol, "S")) tmp1 = 601;
            } else if atom.atom_symbol == "S" {
                offset = 601;
            // Avalon❗✔️:               else if (0 == strcmp(ap->atom_symbol, "P")) tmp1 = 701;
            } else if atom.atom_symbol == "P" {
                offset = 701;
            // Avalon❗✔️:               else if (AtomSymbolMatch(ap->atom_symbol,"F,Cl,Br,I,At"))
            } else if atom_symbol_match(&atom.atom_symbol, "F,Cl,Br,I,At") {
                // Avalon❗✔️:                   tmp1 = 901;
                offset = 901;
            }
            // Avalon❗✔️:               extcon[j] = ap->rsize_flags+tmp1;
            extcon[atom_index] = (atom.rsize_flags as i32).wrapping_add(offset);
            // Avalon❗✔️:           }
        }
        // Avalon❗✔️:           for (i=0; i<32; i++)
        // Avalon❗✔️:           {
        // Avalon❗✔️:               for (j=0; j<mp->n_atoms; j++) extcon2[j] = 0;
        // Avalon❗✔️:               for (j=0; j<mp->n_atoms; j++)
        // Avalon❗✔️:               {
        // Avalon❗✔️:                   /* skip non-ring atoms */
        // Avalon❗✔️:                   if (atom_status[j] <= 0) continue;
        // Avalon❗✔️:                   extcon2[j] = atom_status[j]*3+(extcon[j]*0xF);
        // Avalon❗✔️:                   sum = 0; prod = 0;
        // Avalon❗✔️:                   for (jj=0; jj<nbp[j].n_ligands; jj++)
        // Avalon❗✔️:                   {
        // Avalon❗✔️:                       // only propagate through ring bonds
        // Avalon❗✔️:                       if (bond_status[nbp[j].bonds[jj]] <= 0) continue;
        // Avalon❗✔️:                       sum += extcon[nbp[j].atoms[jj]];
        // Avalon❗✔️:                       prod *= 0xFF&(1+extcon[nbp[j].atoms[jj]]);
        // Avalon❗✔️:                   }
        // Avalon❗✔️:                   extcon2[j] = ((extcon2[j]+(sum*191)+prod)<<8) + ((extcon2[j]&0xFF0000)>>16);
        // Avalon❗✔️:                   extcon2[j] &= 0xFFFFFF;
        // Avalon❗✔️:               }
        // Avalon❗✔️:               for (j=0; j<mp->n_atoms; j++) extcon[j] = extcon2[j];
        // Avalon❗✔️:           }
        propagate_ring_connectivity(state, &mut extcon, &mut extcon2);
        // Avalon❗✔️:           /* propagate smallest hash to all members of ring system */
        // Avalon❗✔️:           for (;;)
        // Avalon❗✔️:           {
        // Avalon❗✔️:               changed = FALSE;
        // Avalon❗✔️:               for (j=0; j<mp->n_atoms; j++)
        // Avalon❗✔️:               {
        // Avalon❗✔️:                   /* skip non-ring atoms */
        // Avalon❗✔️:                   if (atom_status[j] <= 0) continue;
        // Avalon❗✔️:                   for (jj=0; jj<nbp[j].n_ligands; jj++)
        // Avalon❗✔️:                   {
        // Avalon❗✔️:                       // only propagate through ring bonds
        // Avalon❗✔️:                       if (bond_status[nbp[j].bonds[jj]] <= 0) continue;
        // Avalon❗✔️:                       if (extcon[j] < extcon[nbp[j].atoms[jj]])
        // Avalon❗✔️:                       {
        // Avalon❗✔️:                           changed = TRUE;
        // Avalon❗✔️:                           extcon[nbp[j].atoms[jj]] = extcon[j];
        // Avalon❗✔️:                       }
        // Avalon❗✔️:                   }
        // Avalon❗✔️:               }
        // Avalon❗✔️:               if (!changed) break;
        // Avalon❗✔️:           }
        propagate_smallest_ring_hash(state, &mut extcon);
        // Avalon❗✔️:           for (j=0; j<mp->n_atoms; j++)
        // Avalon❗✔️:           {
        for atom_index in 0..molecule.atoms.len() {
            // Avalon❗✔️:               if (j+1 == exclude_atom) continue;
            if atom_index as i32 + 1 == exclude_atom {
                continue;
            }
            // Avalon❗✔️:               ADD_BIT(fp_counts, ncounts, NEXT_SEED(NEXT_SEED(seed,4), extcon[j]*3013));
            add_bit(
                counts,
                next_hash(next_hash(seed, 4), extcon[atom_index].wrapping_mul(3013) as u64),
            );
            // Avalon❗✔️: // ADD_BIT(fp_counts, ncounts, NEXT_SEED(NEXT_SEED(seed,4), extcon[j]*7013));
            // Avalon❗✔️: // ADD_BIT(fp_counts, ncounts, extcon[j]*16013 + 4*seed);
            // Avalon❗✔️: 	      result += 1;
            result += 1;
            // Avalon❗✔️:           }
        }
        // Avalon❗✔️:       }
    }
    // Avalon❗✔️:        MyFree((char *)extcon);
    // Avalon❗✔️:        MyFree((char *)extcon2);
    // Rust drops both vectors here.
    // Avalon❗✔️:    }
    result
}

fn propagate_ring_connectivity(state: &FingerprintPreprocessingState, extcon: &mut [i32], extcon2: &mut [i32]) {
    // Avalon❗✔️:       for (i=0; i<32; i++)
    // Avalon❗✔️:       {
    for _ in 0..32 {
        // Avalon❗✔️:           for (j=0; j<mp->n_atoms; j++) extcon2[j] = 0;
        extcon2.fill(0);
        // Avalon❗✔️:           for (j=0; j<mp->n_atoms; j++)
        // Avalon❗✔️:           {
        for atom_index in 0..extcon.len() {
            // Avalon❗✔️:               /* skip non-ring atoms */
            // Avalon❗✔️:               if (atom_status[j] <= 0) continue;
            if state.atom_status[atom_index] <= 0 {
                continue;
            }
            // Avalon❗✔️:               extcon2[j] = atom_status[j]*3+(extcon[j]*0xF);
            extcon2[atom_index] = state.atom_status[atom_index]
                .wrapping_mul(3)
                .wrapping_add(extcon[atom_index].wrapping_mul(0x0f));
            // Avalon❗✔️:               sum = 0; prod = 0;
            let mut sum = 0_i32;
            let mut product = 0_i32;
            // Avalon❗✔️:               for (jj=0; jj<nbp[j].n_ligands; jj++)
            // Avalon❗✔️:               {
            for (&neighbour, &bond_index) in state.neighbours[atom_index]
                .atoms()
                .iter()
                .zip(state.neighbours[atom_index].bonds())
            {
                // Avalon❗✔️:                   // only propagate through ring bonds
                // Avalon❗✔️:                   if (bond_status[nbp[j].bonds[jj]] <= 0) continue;
                if state.bond_status[bond_index] <= 0 {
                    continue;
                }
                // Avalon❗✔️:                   sum += extcon[nbp[j].atoms[jj]];
                sum = sum.wrapping_add(extcon[neighbour]);
                // Avalon❗✔️:                   prod *= 0xFF&(1+extcon[nbp[j].atoms[jj]]);
                product = product.wrapping_mul(0xff & extcon[neighbour].wrapping_add(1));
                // Avalon❗✔️:               }
            }
            let base = extcon2[atom_index];
            // Avalon❗✔️:               extcon2[j] =
            // Avalon❗✔️:                   ((extcon2[j]+(sum*191)+prod)<<8) +
            // Avalon❗✔️:                   ((extcon2[j]&0xFF0000)>>16);
            extcon2[atom_index] = base
                .wrapping_add(sum.wrapping_mul(191))
                .wrapping_add(product)
                .wrapping_shl(8)
                .wrapping_add((base & 0x00ff_0000) >> 16);
            // Avalon❗✔️:               extcon2[j] &= 0xFFFFFF;
            extcon2[atom_index] &= 0x00ff_ffff;
            // Avalon❗✔️:           }
        }
        // Avalon❗✔️:           for (j=0; j<mp->n_atoms; j++) extcon[j] = extcon2[j];
        extcon.copy_from_slice(extcon2);
        // Avalon❗✔️:       }
    }
}

fn propagate_smallest_ring_hash(state: &FingerprintPreprocessingState, extcon: &mut [i32]) {
    // Avalon❗✔️:       /* propagate smallest hash to all members of ring system */
    // Avalon❗✔️:       for (;;)
    // Avalon❗✔️:       {
    loop {
        // Avalon❗✔️:           changed = FALSE;
        let mut changed = false;
        // Avalon❗✔️:           for (j=0; j<mp->n_atoms; j++)
        // Avalon❗✔️:           {
        for atom_index in 0..extcon.len() {
            // Avalon❗✔️:               /* skip non-ring atoms */
            // Avalon❗✔️:               if (atom_status[j] <= 0) continue;
            if state.atom_status[atom_index] <= 0 {
                continue;
            }
            // Avalon❗✔️:               for (jj=0; jj<nbp[j].n_ligands; jj++)
            // Avalon❗✔️:               {
            for (&neighbour, &bond_index) in state.neighbours[atom_index]
                .atoms()
                .iter()
                .zip(state.neighbours[atom_index].bonds())
            {
                // Avalon❗✔️:                   // only propagate through ring bonds
                // Avalon❗✔️:                   if (bond_status[nbp[j].bonds[jj]] <= 0) continue;
                if state.bond_status[bond_index] <= 0 {
                    continue;
                }
                // Avalon❗✔️:                   if (extcon[j] < extcon[nbp[j].atoms[jj]])
                // Avalon❗✔️:                   {
                if extcon[atom_index] < extcon[neighbour] {
                    // Avalon❗✔️:                       changed = TRUE;
                    changed = true;
                    // Avalon❗✔️:                       extcon[nbp[j].atoms[jj]] = extcon[j];
                    extcon[neighbour] = extcon[atom_index];
                    // Avalon❗✔️:                   }
                }
                // Avalon❗✔️:               }
            }
            // Avalon❗✔️:           }
        }
        // Avalon❗✔️:           if (!changed) break;
        if !changed {
            break;
        }
        // Avalon❗✔️:       }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::properties::avalon_fingerprint::fingerprint_state::USE_DY_AROMATICITY;
    use crate::properties::avalon_fingerprint::fingerprint_state::prepare_fingerprint_state;
    use crate::properties::avalon_fingerprint::hash::hash_string;
    use crate::properties::avalon_fingerprint::reaccs::{Atom, Bond};

    const NON_SSS_MASKS: &[u32] = &[0x100000, 0x200000, 0x400000, 0x800000, 0xf00000];

    // Generated by `tmp/parity-audit/avalon_non_sss_goldens.c` against the
    // checksum-pinned Avalon source and the pinned RDKit native library.
    const NATIVE_GOLDENS: &str = r#"non_sss=scaffold_mask_100000_q0_dy0 result=7 counts=7:3,21:4, bits=8000200000000000
non_sss=scaffold_mask_200000_q0_dy0 result=8 counts=3:1,21:4,22:3, bits=0800600000000000
non_sss=scaffold_mask_400000_q0_dy0 result=8 counts=10:1,17:1,25:2,26:2,39:1,44:1, bits=0004020680100000
non_sss=scaffold_mask_800000_q0_dy0 result=0 counts= bits=0000000000000000
non_sss=scaffold_mask_f00000_q0_dy0 result=23 counts=3:1,7:3,10:1,17:1,21:8,22:3,25:2,26:2,39:1,44:1, bits=8804620680100000
non_sss=scaffold_mask_100000_q1_dy0 result=0 counts= bits=0000000000000000
non_sss=scaffold_mask_200000_q1_dy0 result=0 counts= bits=0000000000000000
non_sss=scaffold_mask_400000_q1_dy0 result=0 counts= bits=0000000000000000
non_sss=scaffold_mask_800000_q1_dy0 result=0 counts= bits=0000000000000000
non_sss=scaffold_mask_f00000_q1_dy0 result=0 counts= bits=0000000000000000
non_sss=scaffold_mask_100000_q0_dy1 result=7 counts=7:3,21:4, bits=8000200000000000
non_sss=scaffold_mask_200000_q0_dy1 result=8 counts=3:1,21:4,22:3, bits=0800600000000000
non_sss=scaffold_mask_400000_q0_dy1 result=8 counts=10:1,17:1,25:2,26:2,39:1,44:1, bits=0004020680100000
non_sss=scaffold_mask_800000_q0_dy1 result=0 counts= bits=0000000000000000
non_sss=scaffold_mask_f00000_q0_dy1 result=23 counts=3:1,7:3,10:1,17:1,21:8,22:3,25:2,26:2,39:1,44:1, bits=8804620680100000
non_sss=scaffold_mask_100000_q1_dy1 result=0 counts= bits=0000000000000000
non_sss=scaffold_mask_200000_q1_dy1 result=0 counts= bits=0000000000000000
non_sss=scaffold_mask_400000_q1_dy1 result=0 counts= bits=0000000000000000
non_sss=scaffold_mask_800000_q1_dy1 result=0 counts= bits=0000000000000000
non_sss=scaffold_mask_f00000_q1_dy1 result=0 counts= bits=0000000000000000
non_sss_exclude=scaffold_atom_0_dy0 result=23 counts=3:1,7:3,10:1,17:1,21:8,22:3,25:2,26:2,39:1,44:1,
non_sss_exclude=scaffold_atom_1_dy0 result=13 counts=3:1,7:2,21:8,28:2,
non_sss_exclude=scaffold_atom_4_dy0 result=13 counts=3:1,7:3,15:3,21:3,22:3,
non_sss_exclude=scaffold_atom_8_dy0 result=22 counts=7:3,10:1,17:1,21:8,22:3,25:2,26:2,39:1,44:1,
non_sss_exclude=scaffold_atom_0_dy1 result=23 counts=3:1,7:3,10:1,17:1,21:8,22:3,25:2,26:2,39:1,44:1,
non_sss_exclude=scaffold_atom_1_dy1 result=13 counts=3:1,7:2,21:8,28:2,
non_sss_exclude=scaffold_atom_4_dy1 result=13 counts=3:1,7:3,15:3,21:3,22:3,
non_sss_exclude=scaffold_atom_8_dy1 result=22 counts=7:3,10:1,17:1,21:8,22:3,25:2,26:2,39:1,44:1,"#;

    fn scaffold_fixture() -> MoleculeState {
        let symbols = ["C", "N", "O", "C", "S", "P", "Cl", "R"];
        let endpoints = [[1, 2], [2, 3], [3, 1], [4, 5], [5, 6], [6, 7], [7, 4], [1, 4], [2, 8]];
        let mut molecule = MoleculeState {
            atoms: symbols
                .into_iter()
                .map(|symbol| Atom {
                    atom_symbol: symbol.to_string(),
                    ..Atom::default()
                })
                .collect(),
            bonds: endpoints
                .into_iter()
                .map(|atoms| Bond {
                    atoms,
                    bond_type: 1,
                    ..Bond::default()
                })
                .collect(),
            ..MoleculeState::default()
        };
        molecule.atoms[7].atext = "ALA".to_string();
        molecule
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

    fn run_case(mask: u32, as_query: bool, daylight: bool, exclude_atom: i32) -> (i32, String, String) {
        let mut molecule = scaffold_fixture();
        let mut counts = vec![0_i32; 64];
        let result = count_non_sss_flag_families(
            &mut molecule,
            &mut counts,
            AvalonFingerprintFlags::from_bits(mask).unwrap(),
            as_query,
            i32::from(daylight) * USE_DY_AROMATICITY,
            exclude_atom,
        )
        .unwrap();
        (result, sparse_text(&counts), packed_bit_hex(&counts))
    }

    #[test]
    fn non_sss_single_combined_query_and_aromaticity_matrix_matches_native() {
        let mut expected = NATIVE_GOLDENS.lines().filter(|line| line.starts_with("non_sss="));
        for daylight in [false, true] {
            for as_query in [false, true] {
                for &mask in NON_SSS_MASKS {
                    let (result, counts, bits) = run_case(mask, as_query, daylight, 0);
                    let actual = format!(
                        "non_sss=scaffold_mask_{mask:06x}_q{}_dy{} result={result} counts={counts} bits={bits}",
                        u8::from(as_query),
                        u8::from(daylight)
                    );
                    assert_eq!(Some(actual.as_str()), expected.next(), "mask {mask:06x}");
                }
            }
        }
        assert_eq!(None, expected.next());
    }

    #[test]
    fn non_sss_exclusion_paths_match_native_counts() {
        let mut expected = NATIVE_GOLDENS
            .lines()
            .filter(|line| line.starts_with("non_sss_exclude="));
        for daylight in [false, true] {
            for exclude_atom in [0, 1, 4, 8] {
                let (result, counts, _) = run_case(0xf00000, false, daylight, exclude_atom);
                let actual = format!(
                    "non_sss_exclude=scaffold_atom_{exclude_atom}_dy{} result={result} counts={counts}",
                    u8::from(daylight)
                );
                assert_eq!(Some(actual.as_str()), expected.next(), "exclude {exclude_atom}");
            }
        }
        assert_eq!(None, expected.next());
    }

    #[test]
    fn shortcut_label_flag_controls_the_source_hash_color() {
        let mut without_shortcut = scaffold_fixture();
        prepare_fingerprint_state(&mut without_shortcut, AvalonFingerprintFlags::SCAFFOLD_IDS, false, 0, 0).unwrap();
        assert_eq!(without_shortcut.atoms[7].color, 0);

        let mut with_shortcut = scaffold_fixture();
        prepare_fingerprint_state(&mut with_shortcut, AvalonFingerprintFlags::SHORTCUT_LABELS, false, 0, 0).unwrap();
        assert_eq!(
            with_shortcut.atoms[7].color,
            ((hash_string("ALA") & 0x00ff_ff00) | 119) as i32
        );
    }
}
