//! Deterministic input matrix, not expected results or a chemistry implementation.
use super::Pair;

pub const LENGTHS: [u64; 5] = [16, 256, 65_536, 1 << 31, u32::MAX as u64];
pub const COUNTS: [(&str, &[i32]); 5] = [
    ("positive", &[1, 2, 17]),
    ("negative", &[-1, -2, -17]),
    ("mixed", &[-17, 1, 23]),
    ("stored_zero", &[0, 0, -1, 1]),
    ("extreme", &[i32::MIN, i32::MAX, 0, -1, 1]),
];
// Masks select six sorted index slots; membership and count are independent.
pub const SHAPES: [(&str, u8, u8); 10] = [
    ("empty_both", 0, 0),
    ("left_empty", 0, 0b111111),
    ("right_empty", 0b111111, 0),
    ("same_keys", 0b111111, 0b111111),
    ("left_before_right", 0b000111, 0b111000),
    ("right_before_left", 0b111000, 0b000111),
    ("interleaved", 0b010101, 0b101010),
    ("partial_overlap", 0b001111, 0b111100),
    ("left_subset", 0b001100, 0b111111),
    ("right_subset", 0b111111, 0b001100),
];
pub const SAMPLES_PER_CELL: usize = 20;
pub const PAIRS: usize = SHAPES.len() * LENGTHS.len() * COUNTS.len() * SAMPLES_PER_CELL;

pub fn generate() -> Vec<Pair> {
    let mut pairs = Vec::with_capacity(PAIRS);
    for (shape, left_mask, right_mask) in SHAPES {
        for length in LENGTHS {
            for (counts, palette) in COUNTS {
                for sample in 0..SAMPLES_PER_CELL {
                    // Rotate six evenly spaced slots. No random state, platform
                    // hashing, SMILES parsing or fingerprint algorithm is involved.
                    let mut keys: Vec<_> = (0..6)
                        .map(|slot| (sample as u64 * 31 + slot * (length / 6)) % length)
                        .collect();
                    keys.sort_unstable();
                    let entries = |mask: u8, side: usize| {
                        keys.iter()
                            .enumerate()
                            .filter(|(slot, _)| mask & (1 << slot) != 0)
                            .map(|(slot, &key)| {
                                (key, palette[(sample + slot + side) % palette.len()])
                            })
                            .collect()
                    };
                    pairs.push(Pair {
                        id: format!("fp5000-v1/{shape}/{length}/{counts}/{sample}"),
                        length,
                        left: entries(left_mask, 0),
                        right: entries(right_mask, 1),
                    });
                }
            }
        }
    }
    pairs
}
