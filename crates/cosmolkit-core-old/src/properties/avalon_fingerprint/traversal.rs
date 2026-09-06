//! Source-shaped traversal and count hashing shared by Avalon feature families.

use std::ops::{BitOr, BitOrAssign};

use super::hash::{hash_position, next_hash};
use super::preprocess::Neighbourhood;
use super::reaccs::MoleculeState;

// Avalon❗✔️: #define RING_PATTERN_SEED      11
#[allow(dead_code)]
pub(super) const RING_PATTERN_SEED: u64 = 11;
// Avalon❗✔️: #define RING_PATH_SEED         13
#[allow(dead_code)]
pub(super) const RING_PATH_SEED: u64 = 13;
// Avalon❗✔️: #define ATOM_SYMBOL_PATH_SEED  17
pub(super) const ATOM_SYMBOL_PATH_SEED: u64 = 17;
// Avalon❗✔️: #define ATOM_CLASS_PATH_SEED   23
#[allow(dead_code)]
pub(super) const ATOM_CLASS_PATH_SEED: u64 = 23;
// Avalon❗✔️: #define ATOM_COUNT_SEED        31
#[allow(dead_code)]
pub(super) const ATOM_COUNT_SEED: u64 = 31;
// Avalon❗✔️: #define AUGMENTED_ATOM_SEED    37
#[allow(dead_code)]
pub(super) const AUGMENTED_ATOM_SEED: u64 = 37;
// Avalon❗✔️: #define HCOUNT_PATH_SEED       41
#[allow(dead_code)]
pub(super) const HCOUNT_PATH_SEED: u64 = 41;
// Avalon❗✔️: #define HCOUNT_CLASS_PATH_SEED 43
#[allow(dead_code)]
pub(super) const HCOUNT_CLASS_PATH_SEED: u64 = 43;
// Avalon❗✔️: #define HCOUNT_PAIR_SEED       47
#[allow(dead_code)]
pub(super) const HCOUNT_PAIR_SEED: u64 = 47;
// Avalon❗✔️: #define BOND_PATH_SEED         53
#[allow(dead_code)]
pub(super) const BOND_PATH_SEED: u64 = 53;
// Avalon❗✔️: #define AUGMENTED_BOND_SEED    61
#[allow(dead_code)]
pub(super) const AUGMENTED_BOND_SEED: u64 = 61;
// Avalon❗✔️: #define RING_SIZE_SEED         67
#[allow(dead_code)]
pub(super) const RING_SIZE_SEED: u64 = 67;
// Avalon❗✔️: #define DEGREE_PATH_SEED       71
#[allow(dead_code)]
pub(super) const DEGREE_PATH_SEED: u64 = 71;
// Avalon❗✔️: #define CLASS_SPIDER_SEED      79
#[allow(dead_code)]
pub(super) const CLASS_SPIDER_SEED: u64 = 79;
// Avalon❗✔️: #define RING_CLOSURE_SEED      101
const RING_CLOSURE_SEED: i32 = 101;
// Avalon❗✔️: #define NON_SSS_SEED          179
#[allow(dead_code)]
pub(super) const NON_SSS_SEED: u64 = 179;
const ANY_COLOR: i32 = 113;
const CSP3: i32 = 19;
const HETERO: i32 = 23;
const TYPE_MASK: i32 = 0x00ff;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct PathFlags(i32);

impl PathFlags {
    pub const PROCESS_RING_CLOSURES: Self = Self(0x0001);
    pub const PROCESS_CHAINS: Self = Self(0x0002);
    pub const FORCED_HETERO_END: Self = Self(0x0004);
    pub const IGNORE_PATH_SYMBOL: Self = Self(0x0008);
    pub const IGNORE_TERM_SYMBOL: Self = Self(0x0010);
    pub const FORCED_RING_PATH: Self = Self(0x0020);
    pub const STOP_AT_HEAVY_ATOM: Self = Self(0x0040);
    pub const DEBUG_PATH: Self = Self(0x0100);

    const fn contains(self, other: Self) -> bool {
        self.0 & other.0 != 0
    }
}

impl BitOr for PathFlags {
    type Output = Self;

    fn bitor(self, rhs: Self) -> Self::Output {
        Self(self.0 | rhs.0)
    }
}

impl BitOrAssign for PathFlags {
    fn bitor_assign(&mut self, rhs: Self) {
        self.0 |= rhs.0;
    }
}

pub(super) fn add_bit(counts: &mut [i32], seed: u64) {
    // Avalon❗✔️: #define ADD_BIT(counts, ncounts, seed) (counts[hash_position(seed,ncounts)]++)
    let index = hash_position(seed, counts.len());
    counts[index] += 1;
}

pub(super) fn add_bit_count(counts: &mut [i32], seed: u64, count: i32) {
    // Avalon❗✔️: #define ADD_BIT_COUNT(counts, ncounts, seed, count) (counts[hash_position(seed,ncounts)]+=count)
    let index = hash_position(seed, counts.len());
    counts[index] += count;
}

#[allow(clippy::too_many_arguments)]
fn set_path_length_flags(
    molecule: &MoleculeState,
    touched_indices: &mut [i32],
    start_index: usize,
    path_length: usize,
    current_index: usize,
    max_size: usize,
    length_matrix: &mut [Vec<i32>],
    neighbours: &[Neighbourhood],
    exclude_atom: i32,
) {
    // Avalon❗✔️: static
    // Avalon❗✔️: void SetPathLengthFlags(struct reaccs_molecule_t *mp,
    // Avalon❗✔️:                         int touched_indices[],
    // Avalon❗✔️:                         int start_index, int path_length, int current_index,
    // Avalon❗✔️:                         int max_size,
    // Avalon❗✔️:                         int **length_matrix,
    // Avalon❗✔️:                         neighbourhood_t nbp[],
    // Avalon❗✔️:                         int exclude_atom,
    // Avalon❗✔️:                         char *prefix)
    // Avalon❗✔️: {
    // Avalon❗✔️:    int i, ai;
    // Avalon❗✔️:
    // Avalon❗✔️:    for (i=0; i<nbp[current_index].n_ligands; i++)
    // Avalon❗✔️:    {
    for &atom_index in neighbours[current_index].atoms() {
        // Avalon❗✔️:       if (path_length+1 > max_size)             /* don't go too far */
        // Avalon❗✔️:          continue;
        if path_length + 1 > max_size {
            continue;
        }
        // Avalon❗✔️:       ai = nbp[current_index].atoms[i];
        // Avalon❗✔️:       if (ai+1 == exclude_atom) continue;
        if atom_index as i32 + 1 == exclude_atom {
            continue;
        }
        // Avalon❗✔️:       if (touched_indices[ai]) continue;   /* don't walk backwards */
        if touched_indices[atom_index] != 0 {
            continue;
        }
        // Avalon❗✔️:       if (mp->atom_array[ai].color == 0) continue;
        if molecule.atoms[atom_index].color == 0 {
            continue;
        }
        // Avalon❗✔️:       touched_indices[ai] = 1;  /* updating */
        touched_indices[atom_index] = 1;
        // Avalon❗✔️:       length_matrix[start_index][ai] |= 1<<(path_length+1);
        length_matrix[start_index][atom_index] |= 1_i32 << (path_length + 1);
        // Avalon❗✔️:       SetPathLengthFlags(mp, touched_indices,
        // Avalon❗✔️:                          start_index, path_length+1, ai,
        // Avalon❗✔️:                          max_size, length_matrix, nbp,
        // Avalon❗✔️:                          exclude_atom,
        // Avalon❗✔️:                          prefix);
        set_path_length_flags(
            molecule,
            touched_indices,
            start_index,
            path_length + 1,
            atom_index,
            max_size,
            length_matrix,
            neighbours,
            exclude_atom,
        );
        // Avalon❗✔️:       touched_indices[ai] = 0;  /* down-dating */
        touched_indices[atom_index] = 0;
        // Avalon❗✔️:    }
    }
    // Avalon❗✔️: }
}

pub(super) fn build_path_length_matrix(
    molecule: &MoleculeState,
    neighbours: &[Neighbourhood],
    max_size: usize,
    exclude_atom: i32,
) -> Vec<Vec<i32>> {
    // Avalon❗✔️:       length_tmp = TypeAlloc(mp->n_atoms*mp->n_atoms, int);
    // Avalon❗✔️:       length_matrix = TypeAlloc(mp->n_atoms, int*);
    let mut length_matrix = vec![vec![0_i32; molecule.atoms.len()]; molecule.atoms.len()];
    // Avalon❗✔️:       for (i=0; i<mp->n_atoms; i++)
    // Avalon❗✔️:       {
    // Avalon❗✔️:          length_matrix[i] = length_tmp+i*mp->n_atoms;
    // Rust rows already point at their zero-initialized owned storage.
    let mut touched_indices = vec![0_i32; molecule.atoms.len()];
    // Avalon❗✔️:          touched_indices[i] = 0;
    // Avalon❗✔️:       }
    // Avalon❗✔️:       for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    // Avalon❗✔️:       {
    for atom_index in 0..molecule.atoms.len() {
        // Avalon❗✔️:           if (i+1 == exclude_atom) continue;
        if atom_index as i32 + 1 == exclude_atom {
            continue;
        }
        // Avalon❗✔️:           touched_indices[i] = 1; /* updating */
        touched_indices[atom_index] = 1;
        // Avalon❗✔️:           SetPathLengthFlags(mp, touched_indices,
        // Avalon❗✔️:                              i, 0, i,
        // Avalon❗✔️:                              12,
        // Avalon❗✔️:                              length_matrix, nbp,
        // Avalon❗✔️:                              exclude_atom,
        // Avalon❗✔️:                              prefix);
        set_path_length_flags(
            molecule,
            &mut touched_indices,
            atom_index,
            0,
            atom_index,
            max_size,
            &mut length_matrix,
            neighbours,
            exclude_atom,
        );
        // Avalon❗✔️:           touched_indices[i] = 0; /* down-dating */
        touched_indices[atom_index] = 0;
        // Avalon❗✔️:       }
    }
    length_matrix
}

#[allow(clippy::too_many_arguments)]
fn special_neighbours_rec(
    molecule: &MoleculeState,
    touched_indices: &mut [i32],
    path_length: usize,
    current_index: usize,
    max_size: usize,
    csp3: &mut [i32],
    hetero: &mut [i32],
    neighbours: &[Neighbourhood],
    exclude_atom: i32,
) {
    // Avalon❗✔️: static
    // Avalon❗✔️: void SpecialNeighboursRec(struct reaccs_molecule_t *mp,
    // Avalon❗✔️:                           int touched_indices[],
    // Avalon❗✔️:                           int path_length,
    // Avalon❗✔️:                           int current_index,
    // Avalon❗✔️:                           int max_size,
    // Avalon❗✔️:                           int csp3[],
    // Avalon❗✔️:                           int hetero[],
    // Avalon❗✔️:                           neighbourhood_t nbp[],
    // Avalon❗✔️:                           int exclude_atom,
    // Avalon❗✔️:                           char *prefix)
    // Avalon❗✔️: {
    // Avalon❗✔️:    int i, ai;
    // Avalon❗✔️:
    // Avalon❗✔️:    for (i=0; i<nbp[current_index].n_ligands; i++)
    // Avalon❗✔️:    {
    for &atom_index in neighbours[current_index].atoms() {
        // Avalon❗✔️:       if (path_length+1 > max_size)             /* don't go too far */
        // Avalon❗✔️:          continue;
        if path_length + 1 > max_size {
            continue;
        }
        // Avalon❗✔️:       ai = nbp[current_index].atoms[i];
        // Avalon❗✔️:       if (ai+1 == exclude_atom) continue;
        if atom_index as i32 + 1 == exclude_atom {
            continue;
        }
        // Avalon❗✔️:       if (touched_indices[ai]) continue;   /* don't walk backwards */
        if touched_indices[atom_index] != 0 {
            continue;
        }
        // Avalon❗✔️:       if (mp->atom_array[ai].color <= 0) continue;
        let color = molecule.atoms[atom_index].color;
        if color <= 0 {
            continue;
        }
        // Avalon❗✔️:       if (mp->atom_array[ai].color == CSP3) csp3[path_length]++;
        if color == CSP3 {
            csp3[path_length] += 1;
        // Avalon❗✔️:       else if (mp->atom_array[ai].color == HETERO) hetero[path_length]++;
        } else if color == HETERO {
            hetero[path_length] += 1;
        }
        // Avalon❗✔️:       /* only carbon-connected SPIDERS (beta atoms are always included) */
        // Avalon❗✔️:       if (mp->atom_array[ai].color == HETERO  &&  path_length > 1)
        // Avalon❗✔️:          continue;
        if color == HETERO && path_length > 1 {
            continue;
        }
        // Avalon❗✔️:       touched_indices[ai] = 1;  /* updating */
        touched_indices[atom_index] = 1;
        // Avalon❗✔️:       SpecialNeighboursRec(mp, touched_indices,
        // Avalon❗✔️:                            path_length+1, ai,
        // Avalon❗✔️:                            max_size, csp3, hetero, nbp,
        // Avalon❗✔️:                            exclude_atom,
        // Avalon❗✔️:                            prefix);
        special_neighbours_rec(
            molecule,
            touched_indices,
            path_length + 1,
            atom_index,
            max_size,
            csp3,
            hetero,
            neighbours,
            exclude_atom,
        );
        // Avalon❗✔️:       touched_indices[ai] = 0;  /* down-dating */
        touched_indices[atom_index] = 0;
        // Avalon❗✔️:    }
    }
    // Avalon❗✔️: }
}

pub(super) fn collect_special_neighbours(
    molecule: &MoleculeState,
    neighbours: &[Neighbourhood],
    start_index: usize,
    max_size: usize,
    exclude_atom: i32,
) -> (Vec<i32>, Vec<i32>) {
    let mut csp3 = vec![0_i32; max_size + 1];
    let mut hetero = vec![0_i32; max_size + 1];
    let mut touched_indices = vec![0_i32; molecule.atoms.len()];
    touched_indices[start_index] = 1;
    special_neighbours_rec(
        molecule,
        &mut touched_indices,
        1,
        start_index,
        max_size,
        &mut csp3,
        &mut hetero,
        neighbours,
        exclude_atom,
    );
    (csp3, hetero)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn set_path_bits_rec(
    molecule: &MoleculeState,
    neighbours: &[Neighbourhood],
    counts: &mut [i32],
    seed: u64,
    touched_indices: &mut [i32],
    nbonds: i32,
    minbonds: i32,
    maxbonds: i32,
    sprout_index: usize,
    _first_index: i32,
    last_index: i32,
    flags: PathFlags,
    exclude_atom: i32,
) -> i32 {
    // Avalon❗✔️: int SetPathBitsRec(struct reaccs_molecule_t *mp,
    // Avalon❗✔️:                    neighbourhood_t *nbp,
    // Avalon❗✔️:                    int *fp_counts, int ncounts,
    // Avalon❗✔️:                    uint64_t seed,
    // Avalon❗✔️:                    int *touched_indices,
    // Avalon❗✔️:                    int nbonds,
    // Avalon❗✔️:                    int minbonds, int maxbonds,
    // Avalon❗✔️:                    int sprout_index,
    // Avalon❗✔️:                    int first_index,
    // Avalon❗✔️:                    int last_index,
    // Avalon❗✔️:                    int flags,
    // Avalon❗✔️:                    int exclude_atom,
    // Avalon❗✔️:                    char *prefix)
    // Avalon❗✔️: {
    // Avalon❗✔️:    int result;
    // Avalon❗✔️:    int i, ai, bi, acolor, bcolor;
    // Avalon❗✔️:    uint64_t old_seed;
    // Avalon❗✔️:    struct reaccs_atom_t *ap;
    // Avalon❗✔️:
    // Avalon❗✔️:    result = 0;
    let mut result = 0_i32;
    // Avalon❗✔️:    old_seed = seed;
    // Rust keeps each branch seed in an immutable local.
    // Avalon❗✔️:    if (nbonds > maxbonds) return (result);
    if nbonds > maxbonds {
        return result;
    }
    // Avalon❗✔️:    for (i=0; i<nbp[sprout_index].n_ligands; i++)
    // Avalon❗✔️:    {
    for (&atom_index, &bond_index) in neighbours[sprout_index]
        .atoms()
        .iter()
        .zip(neighbours[sprout_index].bonds())
    {
        // Avalon❗✔️:       ai = nbp[sprout_index].atoms[i];
        // Avalon❗✔️:       if (ai == last_index) continue;
        if atom_index as i32 == last_index {
            continue;
        }
        // Avalon❗✔️:       if (ai+1 == exclude_atom) continue;
        if atom_index as i32 + 1 == exclude_atom {
            continue;
        }
        // Avalon❗✔️:       bi = nbp[sprout_index].bonds[i];
        // Avalon❗✔️:       ap = mp->atom_array+ai;
        // Avalon❗✔️:       if ((flags&FORCED_RING_PATH) && mp->bond_array[bi].rsize_flags == 0)
        // Avalon❗✔️:          continue;
        if flags.contains(PathFlags::FORCED_RING_PATH)
            && molecule.bonds[bond_index].rsize_flags == 0
        {
            continue;
        }
        let atom_color = molecule.atoms[atom_index].color;
        // Avalon❗✔️:       if (ap->color >= 18  &&  ap->color != ANY_COLOR  &&  0 != (flags & STOP_AT_HEAVY_ATOM)) continue;
        if atom_color >= 18
            && atom_color != ANY_COLOR
            && flags.contains(PathFlags::STOP_AT_HEAVY_ATOM)
        {
            continue;
        }
        // Avalon❗✔️:       if (touched_indices[ai] > 0)      /* ring closure */
        // Avalon❗✔️:       {
        if touched_indices[atom_index] > 0 {
            // Avalon❗✔️:          if (0 == (flags & PROCESS_RING_CLOSURES)) continue;
            if !flags.contains(PathFlags::PROCESS_RING_CLOSURES) {
                continue;
            }
            // Avalon❗✔️:          if (touched_indices[ai] > 1) continue; // not just a plain ring
            if touched_indices[atom_index] > 1 {
                continue;
            }
            // Avalon❗✔️:          bcolor = mp->bond_array[bi].color;
            let bond_color = molecule.bonds[bond_index].color;
            // Avalon❗✔️:          if (bcolor == 0) continue;
            if bond_color == 0 {
                continue;
            }
            // Avalon❗✔️:          acolor = mp->atom_array[ai].color;
            // Avalon❗✔️:          if (acolor ==  0  &&  0 == (flags & IGNORE_PATH_SYMBOL)) continue;
            if atom_color == 0 && !flags.contains(PathFlags::IGNORE_PATH_SYMBOL) {
                continue;
            }
            // Avalon❗✔️:          old_seed = seed;
            // Avalon❗✔️:          seed = NEXT_SEED(seed, bcolor*16);
            let mut branch_seed = next_hash(seed, (bond_color * 16) as u64);
            // Avalon❗✔️:          if (0 == (flags & IGNORE_PATH_SYMBOL))
            // Avalon❗✔️:          {
            if !flags.contains(PathFlags::IGNORE_PATH_SYMBOL) {
                // Avalon❗✔️:             seed = NEXT_SEED(seed, acolor);
                branch_seed = next_hash(branch_seed, atom_color as u64);
                // Avalon❗✔️:          }
            }
            // Avalon❗✔️:          seed = NEXT_SEED(seed, touched_indices[ai]);
            branch_seed = next_hash(branch_seed, touched_indices[atom_index] as u64);
            // Avalon❗✔️:          if (flags & DEBUG_PATH) printTouched(mp, touched_indices, ai, seed, "ring");
            // Debug output has no fingerprint behavior.
            // Avalon❗✔️:          ADD_BIT(fp_counts, ncounts, seed);
            add_bit(counts, branch_seed);
            // Avalon❗✔️:          seed = NEXT_SEED(seed, RING_CLOSURE_SEED*(nbonds-touched_indices[ai]));
            branch_seed = next_hash(
                branch_seed,
                (RING_CLOSURE_SEED * (nbonds - touched_indices[atom_index])) as u64,
            );
            // Avalon❗✔️:          ADD_BIT(fp_counts, ncounts, seed);
            add_bit(counts, branch_seed);
            // Avalon❗✔️:          result++;
            result += 1;
            // Avalon❗✔️:          seed = old_seed;
            // The immutable parent seed remains unchanged.
            // Avalon❗✔️:       }
            // Avalon❗✔️:       else                              /* normal path */
            // Avalon❗✔️:       {
        } else {
            // Avalon❗✔️:          bcolor = mp->bond_array[bi].color;
            let bond_color = molecule.bonds[bond_index].color;
            // Avalon❗✔️:          if (bcolor == 0) continue;
            if bond_color == 0 {
                continue;
            }
            // Avalon❗✔️:          acolor = mp->atom_array[ai].color;
            // Avalon❗✔️:          if (acolor == 0  &&  0 == (flags & IGNORE_PATH_SYMBOL)) continue;
            if atom_color == 0 && !flags.contains(PathFlags::IGNORE_PATH_SYMBOL) {
                continue;
            }
            // Avalon❗✔️:          /* saving and updating */
            // Avalon❗✔️:          touched_indices[ai] = nbonds+1;
            touched_indices[atom_index] = nbonds + 1;
            // Avalon❗✔️:          old_seed = seed;
            // Avalon❗✔️:          seed = NEXT_SEED(seed, bcolor*16);
            let mut branch_seed = next_hash(seed, (bond_color * 16) as u64);
            // Avalon❗✔️:          if (0 == (flags & IGNORE_PATH_SYMBOL))
            // Avalon❗✔️:          {
            if !flags.contains(PathFlags::IGNORE_PATH_SYMBOL) {
                // Avalon❗✔️:             seed = NEXT_SEED(seed, acolor);
                branch_seed = next_hash(branch_seed, atom_color as u64);
                // Avalon❗✔️:          }
            }
            // Avalon❗✔️:          if (nbonds >= minbonds  &&
            // Avalon❗✔️:              (acolor > 0  ||  0 != (flags & IGNORE_TERM_SYMBOL))  &&
            // Avalon❗✔️:              0 != (flags & PROCESS_CHAINS))
            // Avalon❗✔️:          {
            if nbonds >= minbonds
                && (atom_color > 0 || flags.contains(PathFlags::IGNORE_TERM_SYMBOL))
                && flags.contains(PathFlags::PROCESS_CHAINS)
            {
                // Avalon❗✔️:             if (0 == (flags & FORCED_HETERO_END)  ||
                // Avalon❗✔️:                 0 != (flags & IGNORE_TERM_SYMBOL) ||
                // Avalon❗✔️:                 acolor != 6)
                // Avalon❗✔️:             {
                if !flags.contains(PathFlags::FORCED_HETERO_END)
                    || flags.contains(PathFlags::IGNORE_TERM_SYMBOL)
                    || atom_color != 6
                {
                    // Avalon❗✔️:                if (acolor > 1)         /* don't hit hydrogens! */
                    // Avalon❗✔️:                {
                    if atom_color > 1 {
                        // Avalon❗✔️:                   if (0 != (flags & IGNORE_TERM_SYMBOL))
                        // Avalon❗✔️:                   {
                        if flags.contains(PathFlags::IGNORE_TERM_SYMBOL) {
                            // Avalon❗✔️:                      ADD_BIT(fp_counts, ncounts, seed);
                            add_bit(counts, branch_seed);
                            // Avalon❗✔️:                   }
                            // Avalon❗✔️:                   else
                            // Avalon❗✔️:                   {
                        } else {
                            // Avalon❗✔️:                      ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed,17*acolor));
                            add_bit(counts, next_hash(branch_seed, (17 * atom_color) as u64));
                            // Avalon❗✔️:                   }
                        }
                        // Avalon❗✔️:                   result++;
                        result += 1;
                        // Avalon❗✔️:                }
                    }
                    // Avalon❗✔️:             }
                }
                // Avalon❗✔️:          }
            }
            // Avalon❗✔️:          if (nbonds+1 <= maxbonds)       /* continue recursion */
            // Avalon❗✔️:          {
            if nbonds + 1 <= maxbonds {
                // Avalon❗✔️:             result += SetPathBitsRec(mp, nbp,
                // Avalon❗✔️:                                      fp_counts, ncounts,
                // Avalon❗✔️:                                      seed, touched_indices,
                // Avalon❗✔️:                                      nbonds+1,
                // Avalon❗✔️:                                      minbonds, maxbonds,
                // Avalon❗✔️:                                      ai,
                // Avalon❗✔️:                                      first_index, sprout_index,
                // Avalon❗✔️:                                      flags, exclude_atom,
                // Avalon❗✔️:                                      prefix);
                result += set_path_bits_rec(
                    molecule,
                    neighbours,
                    counts,
                    branch_seed,
                    touched_indices,
                    nbonds + 1,
                    minbonds,
                    maxbonds,
                    atom_index,
                    _first_index,
                    sprout_index as i32,
                    flags,
                    exclude_atom,
                );
                // Avalon❗✔️:          }
            }
            // Avalon❗✔️:
            // Avalon❗✔️:          /* restoring and down dating */
            // Avalon❗✔️:          touched_indices[ai] = 0;
            touched_indices[atom_index] = 0;
            // Avalon❗✔️:          seed = old_seed;
            // The immutable parent seed remains unchanged.
            // Avalon❗✔️:       }
        }
        // Avalon❗✔️:    }
    }
    // Avalon❗✔️:
    // Avalon❗✔️:    return result;
    result
    // Avalon❗✔️: }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn set_feature_bits(
    molecule: &MoleculeState,
    counts: &mut [i32],
    start_flags: i32,
    end_flags: i32,
    path_min: usize,
    path_max: usize,
    use_counts: bool,
    use_atom_types: bool,
    length_matrix: &[Vec<i32>],
    start_seed: u64,
    exclude_atom: i32,
) -> i32 {
    // Avalon❗✔️: int SetFeatureBits(struct reaccs_molecule_t *mp,
    // Avalon❗✔️:                    int *fp_counts,
    // Avalon❗✔️:                    int ncounts,
    // Avalon❗✔️:                    int start_flags,
    // Avalon❗✔️:                    int end_flags,
    // Avalon❗✔️:                    int path_min,
    // Avalon❗✔️:                    int path_max,
    // Avalon❗✔️:                    int use_counts,
    // Avalon❗✔️:                    int use_atom_types,
    // Avalon❗✔️:                    int **length_matrix,
    // Avalon❗✔️:                    uint64_t start_seed,
    // Avalon❗✔️:                    int exclude_atom,
    // Avalon❗✔️:                    char *prefix)
    // Avalon❗✔️: {
    // Avalon❗✔️:    int result = 0;
    let mut result = 0_i32;
    // Avalon❗✔️:    int coli, colj;
    // Avalon❗✔️:    int i, j, k;
    // Avalon❗✔️:    uint64_t seed_i, seed;
    // Avalon❗✔️:    int *counts;
    // Avalon❗✔️:
    // Avalon❗✔️:    counts = TypeAlloc(ncounts*4, int); /* allocate tmp array for counts */
    let mut feature_counts = vec![0_i32; counts.len() * 4];
    // Avalon❗✔️:    for (i=0; i<mp->n_atoms; i++)
    // Avalon❗✔️:    {
    for atom_i in 0..molecule.atoms.len() {
        // Avalon❗✔️:       if (i+1 == exclude_atom) continue;
        if atom_i as i32 + 1 == exclude_atom {
            continue;
        }
        // Avalon❗✔️:       coli = mp->atom_array[i].color;
        let color_i = molecule.atoms[atom_i].color;
        // Avalon❗✔️:       if (0 == (coli&start_flags)) continue;
        if color_i & start_flags == 0 {
            continue;
        }
        // Avalon❗✔️:       if (use_atom_types)
        // Avalon❗✔️:       {
        let seed_i = if use_atom_types {
            // Avalon❗✔️:          if (0 == (coli&TYPE_MASK)) continue; // ignore generic atoms
            if color_i & TYPE_MASK == 0 {
                continue;
            }
            // Avalon❗✔️:          seed_i = NEXT_SEED(start_seed, coli&TYPE_MASK);
            next_hash(start_seed, (color_i & TYPE_MASK) as u64)
            // Avalon❗✔️:       }
            // Avalon❗✔️:       else                seed_i = start_seed;
        } else {
            start_seed
        };
        // Avalon❗✔️:       for (j=0; j<mp->n_atoms; j++)
        // Avalon❗✔️:       {
        for atom_j in 0..molecule.atoms.len() {
            // Avalon❗✔️:          if (j+1 == exclude_atom) continue;
            if atom_j as i32 + 1 == exclude_atom {
                continue;
            }
            // Avalon❗✔️:          colj = mp->atom_array[j].color;
            let color_j = molecule.atoms[atom_j].color;
            // Avalon❗✔️:          if (0 == (colj&end_flags)) continue;
            if color_j & end_flags == 0 {
                continue;
            }
            // Avalon❗✔️:          if (use_atom_types)
            // Avalon❗✔️:          {
            let seed = if use_atom_types {
                // Avalon❗✔️:             if (0 == (colj&TYPE_MASK)) continue; // ignore generic atoms
                if color_j & TYPE_MASK == 0 {
                    continue;
                }
                // Avalon❗✔️:             seed = NEXT_SEED(seed_i, colj&TYPE_MASK);
                next_hash(seed_i, (color_j & TYPE_MASK) as u64)
                // Avalon❗✔️:          }
                // Avalon❗✔️:          else                seed = seed_i;
            } else {
                seed_i
            };
            // Avalon❗✔️:          for (k=path_min; k<=path_max; k++)
            for path_length in path_min..=path_max {
                // Avalon❗✔️:             if ((1<<k)&length_matrix[i][j])
                // Avalon❗✔️:             {
                if (1_i32 << path_length) & length_matrix[atom_i][atom_j] != 0 {
                    // Avalon❗✔️:                /* count the features */
                    // Avalon❗✔️:                counts[(k*19+seed)%(ncounts*4)]++;
                    let index = ((path_length as u64 * 19).wrapping_add(seed)
                        % feature_counts.len() as u64) as usize;
                    feature_counts[index] += 1;
                    // Avalon❗✔️:                result++;
                    result += 1;
                    // Avalon❗✔️:             }
                }
            }
            // Avalon❗✔️:       }
        }
        // Avalon❗✔️:    }
    }
    // Avalon❗✔️:    /* Set the bits */
    // Avalon❗✔️:    for (i=0; i<ncounts*4; i++)
    // Avalon❗✔️:    {
    for (index, &count) in feature_counts.iter().enumerate() {
        // Avalon❗✔️:       if (counts[i] > 0)
        // Avalon❗✔️:       {
        if count > 0 {
            // Avalon❗✔️:          ADD_BIT(fp_counts, ncounts, i);
            add_bit(counts, index as u64);
            // Avalon❗✔️:          result++;
            result += 1;
            // Avalon❗✔️:       }
        }
        // Avalon❗✔️:       if (use_counts  &&  counts[i] > 1)
        // Avalon❗✔️:       {
        if use_counts && count > 1 {
            // Avalon❗✔️:          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i,19), counts[i]);
            add_bit_count(counts, next_hash(index as u64, 19), count);
            // Avalon❗✔️:          result++;
            result += 1;
            // Avalon❗✔️:       }
        }
        // Avalon❗✔️:       if (use_counts  &&  counts[i] > 2)
        // Avalon❗✔️:       {
        if use_counts && count > 2 {
            // Avalon❗✔️:          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i,29), counts[i]);
            add_bit_count(counts, next_hash(index as u64, 29), count);
            // Avalon❗✔️:          result++;
            result += 1;
            // Avalon❗✔️:       }
        }
        // Avalon❗✔️:       if (use_counts  &&  counts[i] > 4)
        // Avalon❗✔️:       {
        if use_counts && count > 4 {
            // Avalon❗✔️:          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i,59), counts[i]);
            add_bit_count(counts, next_hash(index as u64, 59), count);
            // Avalon❗✔️:          result++;
            result += 1;
            // Avalon❗✔️:       }
        }
        // Avalon❗✔️:       if (use_counts  &&  counts[i] > 8)
        // Avalon❗✔️:       {
        if use_counts && count > 8 {
            // Avalon❗✔️:          ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i,79), counts[i]);
            add_bit_count(counts, next_hash(index as u64, 79), count);
            // Avalon❗✔️:          result++;
            result += 1;
            // Avalon❗✔️:       }
        }
        // Avalon❗✔️:    }
    }
    // Avalon❗✔️:    MyFree((char *)counts);
    // Rust releases the temporary vector on scope exit.
    // Avalon❗✔️:    return (result);
    result
    // Avalon❗✔️: }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn set_ring_size_pair_bits(
    molecule: &MoleculeState,
    counts: &mut [i32],
    start_flags: i32,
    end_flags: i32,
    path_min: usize,
    path_max: usize,
    length_matrix: &[Vec<i32>],
    start_seed: u64,
    exclude_atom: i32,
) -> i32 {
    // Avalon❗✔️: int SetRingSizePairBits(struct reaccs_molecule_t *mp,
    // Avalon❗✔️:                         int *fp_counts,
    // Avalon❗✔️:                         int ncounts,
    // Avalon❗✔️:                         int start_flags,
    // Avalon❗✔️:                         int end_flags,
    // Avalon❗✔️:                         int path_min,
    // Avalon❗✔️:                         int path_max,
    // Avalon❗✔️:                         int **length_matrix,
    // Avalon❗✔️:                         uint64_t start_seed,
    // Avalon❗✔️:                         int exclude_atom,
    // Avalon❗✔️:                         char *prefix)
    // Avalon❗✔️: {
    // Avalon❗✔️:    int result = 0;
    let mut result = 0_i32;
    // Avalon❗✔️:    int coli, colj;
    // Avalon❗✔️:    int rsizei, rsizej;
    // Avalon❗✔️:    int i, j, k;
    // Avalon❗✔️:    int ri, rj;
    // Avalon❗✔️:    uint64_t seed_i, seed;
    // Avalon❗✔️:    for (i=0; i<mp->n_atoms; i++)
    // Avalon❗✔️:    {
    for atom_i in 0..molecule.atoms.len() {
        // Avalon❗✔️:       if (i+1 == exclude_atom) continue;
        if atom_i as i32 + 1 == exclude_atom {
            continue;
        }
        // Avalon❗✔️:       coli = mp->atom_array[i].color;
        // Avalon❗✔️:       if (0 == (coli&start_flags)) continue;
        if molecule.atoms[atom_i].color & start_flags == 0 {
            continue;
        }
        // Avalon❗✔️:       rsizei = mp->atom_array[i].rsize_flags;
        let ring_sizes_i = molecule.atoms[atom_i].rsize_flags;
        // Avalon❗✔️:       if (0 == rsizei) continue;
        if ring_sizes_i == 0 {
            continue;
        }
        // Avalon❗✔️:       seed_i = start_seed;
        let seed_i = start_seed;
        // Avalon❗✔️:       for (j=0; j<mp->n_atoms; j++)
        // Avalon❗✔️:       {
        for atom_j in 0..molecule.atoms.len() {
            // Avalon❗✔️:          if (j+1 == exclude_atom) continue;
            if atom_j as i32 + 1 == exclude_atom {
                continue;
            }
            // Avalon❗✔️:          colj = mp->atom_array[j].color;
            // Avalon❗✔️:          if (0 == (colj&end_flags)) continue;
            if molecule.atoms[atom_j].color & end_flags == 0 {
                continue;
            }
            // Avalon❗✔️:          rsizej = mp->atom_array[j].rsize_flags;
            let ring_sizes_j = molecule.atoms[atom_j].rsize_flags;
            // Avalon❗✔️:          if (0 == rsizej) continue;
            if ring_sizes_j == 0 {
                continue;
            }
            // Avalon❗✔️:          seed = seed_i;
            let seed = seed_i;
            // Avalon❗✔️:          for (k=path_min; k<=path_max; k++)
            for path_length in path_min..=path_max {
                // Avalon❗✔️:             if ((1<<k)&length_matrix[i][j])
                // Avalon❗✔️:             {
                if (1_i32 << path_length) & length_matrix[atom_i][atom_j] == 0 {
                    continue;
                }
                // Avalon❗✔️:                for (ri=3; ri<15; ri++)
                // Avalon❗✔️:                {
                for ring_i in 3_u32..15 {
                    // Avalon❗✔️:                   if (0 == (rsizei&(1<<ri))) continue;
                    if ring_sizes_i & (1_u32 << ring_i) == 0 {
                        continue;
                    }
                    // Avalon❗✔️:                   for (rj=3; rj<15; rj++)
                    // Avalon❗✔️:                   {
                    for ring_j in 3_u32..15 {
                        // Avalon❗✔️:                      if (0 == (rsizej&(1<<rj))) continue;
                        if ring_sizes_j & (1_u32 << ring_j) == 0 {
                            continue;
                        }
                        // Avalon❗✔️:                      ADD_BIT(fp_counts, ncounts, NEXT_SEED(NEXT_SEED(seed, (k/2)*119), 37*ri*rj));
                        let pair_seed = next_hash(
                            next_hash(seed, (path_length as u64 / 2) * 119),
                            37 * u64::from(ring_i) * u64::from(ring_j),
                        );
                        add_bit(counts, pair_seed);
                        // Avalon❗✔️:                      result++;
                        result += 1;
                        // Avalon❗✔️:                   }
                    }
                    // Avalon❗✔️:                }
                }
                // Avalon❗✔️:             }
            }
            // Avalon❗✔️:       }
        }
        // Avalon❗✔️:    }
    }
    // Avalon❗✔️:    return (result);
    result
    // Avalon❗✔️: }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::properties::avalon_fingerprint::preprocess::setup_neighbourhood;
    use crate::properties::avalon_fingerprint::reaccs::{Atom, Bond};

    fn molecule(atom_colors: &[i32], endpoints: &[[i32; 2]], bond_colors: &[i32]) -> MoleculeState {
        MoleculeState {
            atoms: atom_colors
                .iter()
                .map(|&color| Atom {
                    color,
                    ..Atom::default()
                })
                .collect(),
            bonds: endpoints
                .iter()
                .zip(bond_colors)
                .map(|(&atoms, &color)| Bond {
                    atoms,
                    color,
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

    #[test]
    fn chain_path_counts_match_native_exported_function() {
        let molecule = molecule(&[6, 7, 8, 6], &[[1, 2], [2, 3], [3, 4]], &[1, 2, 1]);
        let neighbours = setup_neighbourhood(&molecule, 4).unwrap();
        let mut touched = vec![0; 4];
        touched[0] = 1;
        let mut counts = vec![0; 64];

        let result = set_path_bits_rec(
            &molecule,
            &neighbours,
            &mut counts,
            ATOM_SYMBOL_PATH_SEED,
            &mut touched,
            1,
            1,
            3,
            0,
            0,
            -1,
            PathFlags::PROCESS_CHAINS,
            0,
        );

        assert_eq!(result, 3);
        assert_eq!(sparse(&counts), vec![(36, 1), (45, 1), (55, 1)]);
        assert_eq!(touched, vec![1, 0, 0, 0]);
    }

    #[test]
    fn triangle_ring_closure_counts_match_native_exported_function() {
        let molecule = molecule(&[6, 6, 6], &[[1, 2], [2, 3], [3, 1]], &[1, 1, 1]);
        let neighbours = setup_neighbourhood(&molecule, 3).unwrap();
        let mut touched = vec![0; 3];
        touched[0] = 1;
        let mut counts = vec![0; 64];

        let result = set_path_bits_rec(
            &molecule,
            &neighbours,
            &mut counts,
            ATOM_SYMBOL_PATH_SEED,
            &mut touched,
            1,
            1,
            3,
            0,
            0,
            -1,
            PathFlags::PROCESS_CHAINS | PathFlags::PROCESS_RING_CLOSURES,
            0,
        );

        assert_eq!(result, 6);
        assert_eq!(sparse(&counts), vec![(30, 2), (33, 2), (55, 2), (56, 2)]);
    }

    #[test]
    fn feature_pair_counts_match_native_exported_function() {
        let molecule = molecule(&[0x0401, 0x0103, 0x0102], &[], &[]);
        let length_matrix = vec![
            vec![0, 1 << 2, 1 << 3],
            vec![1 << 2, 0, 0],
            vec![1 << 3, 0, 0],
        ];
        let mut counts = vec![0; 64];

        let result = set_feature_bits(
            &molecule,
            &mut counts,
            0x0400,
            0x0100,
            1,
            4,
            true,
            true,
            &length_matrix,
            4237,
            0,
        );

        assert_eq!(result, 4);
        assert_eq!(sparse(&counts), vec![(9, 1), (55, 1)]);
    }

    #[test]
    fn feature_count_thresholds_match_native_exported_function() {
        let molecule = molecule(&[0x0101; 4], &[], &[]);
        let mut length_matrix = vec![vec![0; 4]; 4];
        for (atom_i, row) in length_matrix.iter_mut().enumerate() {
            for (atom_j, distance_flags) in row.iter_mut().enumerate() {
                if atom_i != atom_j {
                    *distance_flags = 1 << 1;
                }
            }
        }
        let mut counts = vec![0; 257];

        let result = set_feature_bits(
            &molecule,
            &mut counts,
            0x0100,
            0x0100,
            1,
            1,
            true,
            true,
            &length_matrix,
            4237,
            0,
        );

        assert_eq!(result, 17);
        assert_eq!(
            sparse(&counts),
            vec![(7, 12), (55, 12), (69, 12), (131, 1), (152, 12)]
        );
    }

    #[test]
    fn ring_size_pair_counts_match_native_exported_function() {
        let mut molecule = molecule(&[0x0201, 0x0103, 0x0102], &[], &[]);
        molecule.atoms[0].rsize_flags = (1 << 5) | (1 << 6);
        molecule.atoms[1].rsize_flags = 1 << 4;
        molecule.atoms[2].rsize_flags = 1 << 7;
        let length_matrix = vec![
            vec![0, 1 << 2, 1 << 3],
            vec![1 << 2, 0, 0],
            vec![1 << 3, 0, 0],
        ];
        let mut counts = vec![0; 64];

        let result = set_ring_size_pair_bits(
            &molecule,
            &mut counts,
            0x0200,
            0x0100,
            1,
            4,
            &length_matrix,
            3237,
            0,
        );

        assert_eq!(result, 4);
        assert_eq!(sparse(&counts), vec![(12, 1), (25, 1), (47, 1), (50, 1)]);
    }

    #[test]
    fn path_length_matrix_keeps_all_simple_path_lengths() {
        let molecule = molecule(&[6, 6, 6, 6], &[[1, 2], [2, 3], [3, 1], [3, 4]], &[1; 4]);
        let neighbours = setup_neighbourhood(&molecule, 4).unwrap();

        let lengths = build_path_length_matrix(&molecule, &neighbours, 12, 0);

        assert_eq!(lengths[0][1], (1 << 1) | (1 << 2));
        assert_eq!(lengths[0][2], (1 << 1) | (1 << 2));
        assert_eq!(lengths[0][3], (1 << 2) | (1 << 3));
    }

    #[test]
    fn special_neighbour_recursion_stops_past_beta_hetero_atoms() {
        let molecule = molecule(
            &[6, HETERO, CSP3, HETERO, CSP3],
            &[[1, 2], [2, 3], [3, 4], [4, 5]],
            &[1; 4],
        );
        let neighbours = setup_neighbourhood(&molecule, 5).unwrap();

        let (csp3, hetero) = collect_special_neighbours(&molecule, &neighbours, 0, 7, 0);

        assert_eq!(hetero[1], 1);
        assert_eq!(csp3[2], 1);
        assert_eq!(hetero[3], 1);
        assert_eq!(csp3[4], 0);
    }

    #[test]
    fn path_control_flags_match_native_exported_function() {
        for (colors, flags, expected_result, expected_counts) in [
            (
                vec![6, 26, 7],
                PathFlags::PROCESS_CHAINS | PathFlags::STOP_AT_HEAVY_ATOM,
                0,
                vec![],
            ),
            (
                vec![6, 0, 7],
                PathFlags::PROCESS_CHAINS | PathFlags::IGNORE_PATH_SYMBOL,
                1,
                vec![(17, 1)],
            ),
            (
                vec![6, 6, 7],
                PathFlags::PROCESS_CHAINS | PathFlags::FORCED_HETERO_END,
                1,
                vec![(37, 1)],
            ),
        ] {
            let molecule = molecule(&colors, &[[1, 2], [2, 3]], &[1, 1]);
            let neighbours = setup_neighbourhood(&molecule, 3).unwrap();
            let mut touched = vec![1, 0, 0];
            let mut counts = vec![0; 64];

            let result = set_path_bits_rec(
                &molecule,
                &neighbours,
                &mut counts,
                ATOM_SYMBOL_PATH_SEED,
                &mut touched,
                1,
                1,
                3,
                0,
                0,
                -1,
                flags,
                0,
            );

            assert_eq!(result, expected_result);
            assert_eq!(sparse(&counts), expected_counts);
        }
    }
}
