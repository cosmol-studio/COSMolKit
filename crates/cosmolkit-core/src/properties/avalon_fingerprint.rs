//! Source-backed Avalon explicit-bit fingerprints.
//!
//! This module follows RDKit's `GetAvalonFP(ROMol, ...)` adapter and the
//! pinned Avalon/REACCS engine for the exposed `n_bits`, `is_query`, and
//! `bit_flags` surface. The returned vector keeps the requested public size;
//! the adapter's four-byte internal rounding and second accumulation pass are
//! preserved in the implementation below.

use std::ops::{BitOr, BitOrAssign};

use crate::{Fingerprint, FingerprintError, Molecule};

mod aromaticity;
mod daylight_aromaticity;
mod fingerprint_state;
mod hash;
mod high_flags;
mod low_flags;
mod middle_flags;
mod non_sss_flags;
mod preprocess;
mod reaccs;
mod rings;
mod symbols;
mod traversal;

use self::fingerprint_state::with_prepared_fingerprint_state;
use self::high_flags::count_high_flag_families_prepared;
use self::low_flags::count_low_flag_families_prepared;
use self::middle_flags::count_middle_flag_families_prepared;
use self::non_sss_flags::count_non_sss_flag_families_prepared;
use self::reaccs::mol_to_reaccs;

/// Avalon feature flags from the source `ssmatch.h` definitions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct AvalonFingerprintFlags(u32);

impl AvalonFingerprintFlags {
    pub const RING_PATTERN: Self = Self(0x000001);
    pub const RING_PATH: Self = Self(0x000002);
    pub const ATOM_SYMBOL_PATH: Self = Self(0x000004);
    pub const ATOM_CLASS_PATH: Self = Self(0x000008);
    pub const ATOM_COUNT: Self = Self(0x000010);
    pub const AUGMENTED_ATOM: Self = Self(0x000020);
    pub const HCOUNT_PATH: Self = Self(0x000040);
    pub const HCOUNT_CLASS_PATH: Self = Self(0x000080);
    pub const HCOUNT_PAIR: Self = Self(0x000100);
    pub const BOND_PATH: Self = Self(0x000200);
    pub const AUGMENTED_BOND: Self = Self(0x000400);
    pub const RING_SIZE_COUNTS: Self = Self(0x000800);
    pub const DEGREE_PATH: Self = Self(0x001000);
    pub const CLASS_SPIDERS: Self = Self(0x002000);
    pub const FEATURE_PAIRS: Self = Self(0x004000);
    pub const SCAFFOLD_IDS: Self = Self(0x100000);
    pub const SCAFFOLD_COLORS: Self = Self(0x200000);
    pub const SCAFFOLD_LINKS: Self = Self(0x400000);
    pub const SHORTCUT_LABELS: Self = Self(0x800000);

    /// C++ `avalonSSSBits` default (`0x007fff`).
    pub const ALL_FEATURES: Self = Self(0x007fff);
    pub const NON_SSS_BITS: Self = Self(0xf00000);
    pub const SIMILARITY: Self = Self(0xf07fff);

    pub const fn bits(self) -> u32 {
        self.0
    }

    pub const fn from_bits(bits: u32) -> Option<Self> {
        if bits & !Self::SIMILARITY.0 == 0 {
            Some(Self(bits))
        } else {
            None
        }
    }

    pub const fn from_bits_retain(bits: u32) -> Self {
        Self(bits)
    }

    pub const fn contains(self, other: Self) -> bool {
        self.0 & other.0 == other.0
    }
}

impl BitOr for AvalonFingerprintFlags {
    type Output = Self;

    fn bitor(self, rhs: Self) -> Self::Output {
        Self(self.0 | rhs.0)
    }
}

impl BitOrAssign for AvalonFingerprintFlags {
    fn bitor_assign(&mut self, rhs: Self) {
        self.0 |= rhs.0;
    }
}

/// Parameters for the source-backed Avalon bit-vector API.
#[derive(Debug, Clone, PartialEq)]
pub struct AvalonFingerprintParams {
    pub n_bits: u32,
    pub is_query: bool,
    pub bit_flags: AvalonFingerprintFlags,
}

impl Default for AvalonFingerprintParams {
    fn default() -> Self {
        Self {
            n_bits: 512,
            is_query: false,
            bit_flags: AvalonFingerprintFlags::ALL_FEATURES,
        }
    }
}

impl AvalonFingerprintParams {
    pub fn validate(&self) -> Result<(), FingerprintError> {
        if self.n_bits < 8 {
            return Err(FingerprintError::InvalidArguments {
                reason: "Avalon n_bits must be at least 8",
            });
        }
        if AvalonFingerprintFlags::from_bits(self.bit_flags.bits()).is_none() {
            return Err(FingerprintError::InvalidArguments {
                reason: "Avalon bit_flags contains undefined source bits",
            });
        }
        Ok(())
    }
}

pub fn avalon_fingerprint(
    molecule: &Molecule,
    params: &AvalonFingerprintParams,
) -> Result<Fingerprint, FingerprintError> {
    params.validate()?;
    let mut molecule_state = mol_to_reaccs(molecule)?;
    let n_bytes = (params.n_bits / 8) as usize;
    let mut bytes = vec![0_u8; n_bytes];
    set_fingerprint_bits(
        &mut molecule_state,
        &mut bytes,
        params.bit_flags,
        params.is_query,
        false,
    )?;
    if !params.is_query {
        set_fingerprint_bits(
            &mut molecule_state,
            &mut bytes,
            params.bit_flags,
            false,
            true,
        )?;
    }
    Ok(Fingerprint::from_lsb_bytes(params.n_bits as usize, &bytes))
}

const ACCUMULATE_BITS: i32 = 0x0002;
const USE_DY_AROMATICITY: i32 = 0x0001;

fn count_fingerprint_patterns(
    molecule: &mut reaccs::MoleculeState,
    counts: &mut [i32],
    bit_flags: AvalonFingerprintFlags,
    is_query: bool,
    fpflags: i32,
) -> Result<i32, FingerprintError> {
    // Avalon✔️✔️: CountFingerprintPatterns allocates and prepares one shared
    // state, then executes every active family before restoring bond types.
    if counts.is_empty() {
        return Err(FingerprintError::InvalidArguments {
            reason: "Avalon fingerprint count array must not be empty",
        });
    }
    if fpflags & ACCUMULATE_BITS == 0 {
        counts.fill(0);
    }
    with_prepared_fingerprint_state(
        molecule,
        bit_flags,
        is_query,
        fpflags,
        0,
        |working, state| {
            let mut result = 0_i32;
            let saved_atom_colors = working
                .atoms
                .iter()
                .map(|atom| atom.color)
                .collect::<Vec<_>>();
            let saved_bond_colors = working
                .bonds
                .iter()
                .map(|bond| bond.color)
                .collect::<Vec<_>>();
            result +=
                count_low_flag_families_prepared(working, state, counts, bit_flags, is_query, 0);
            for (atom, color) in working.atoms.iter_mut().zip(saved_atom_colors) {
                atom.color = color;
            }
            for (bond, color) in working.bonds.iter_mut().zip(saved_bond_colors) {
                bond.color = color;
            }
            result +=
                count_middle_flag_families_prepared(working, state, counts, bit_flags, is_query, 0);
            result +=
                count_high_flag_families_prepared(working, state, counts, bit_flags, is_query, 0);
            if bit_flags.bits() & AvalonFingerprintFlags::NON_SSS_BITS.bits() != 0 {
                result += count_non_sss_flag_families_prepared(
                    working, state, counts, bit_flags, is_query, 0,
                );
            }
            Ok(result)
        },
    )
}

fn set_fingerprint_bits(
    molecule: &mut reaccs::MoleculeState,
    fingerprint: &mut [u8],
    bit_flags: AvalonFingerprintFlags,
    is_query: bool,
    accumulate: bool,
) -> Result<i32, FingerprintError> {
    // Avalon❗✔️: int* fp_counts = TypeAlloc(nbytes*8, int);
    // Avalon❗✔️: result = SetFingerprintCountsWithFocus(mp, fp_counts, nbytes*8,
    // Avalon❗✔️:                                    which_bits, as_query, fpflags, 0);
    let internal_nbytes = (fingerprint.len() + 3) & !3;
    let mut counts = vec![0_i32; internal_nbytes * 8];
    let fpflags = if accumulate {
        ACCUMULATE_BITS | USE_DY_AROMATICITY
    } else {
        0
    };
    let result = count_fingerprint_patterns(molecule, &mut counts, bit_flags, is_query, fpflags)?;
    // Avalon❗✔️: for (i=0; i<nbytes*8; i++)
    // Avalon❗✔️:    if (fp_counts[i] > 0) SET_BIT(fingerprint, nbytes, i);
    for (bit, &count) in counts.iter().take(fingerprint.len() * 8).enumerate() {
        if count > 0 {
            fingerprint[bit / 8] |= 1_u8 << (bit % 8);
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_molecule_reports_source_conversion_error() {
        let err = avalon_fingerprint(&Molecule::new(), &AvalonFingerprintParams::default())
            .expect_err("empty molecules cannot be serialized to Avalon REACCS");
        assert!(matches!(err, FingerprintError::AvalonConversion { .. }));
    }

    #[test]
    fn source_cpp_default_profile_is_typed() {
        let params = AvalonFingerprintParams::default();
        assert_eq!(params.n_bits, 512);
        assert!(!params.is_query);
        assert_eq!(params.bit_flags, AvalonFingerprintFlags::ALL_FEATURES);
    }

    #[test]
    fn source_python_profile_is_distinct_from_cpp_default() {
        assert_eq!(AvalonFingerprintFlags::SIMILARITY.bits(), 0xf07fff);
        assert_ne!(
            AvalonFingerprintParams::default().bit_flags,
            AvalonFingerprintFlags::SIMILARITY
        );
    }

    #[test]
    fn typed_flags_reject_unknown_bits() {
        assert!(AvalonFingerprintFlags::from_bits(0x80000000).is_none());
        let params = AvalonFingerprintParams {
            bit_flags: AvalonFingerprintFlags::from_bits_retain(0x80000000),
            ..Default::default()
        };
        assert!(matches!(
            avalon_fingerprint(&Molecule::new(), &params),
            Err(FingerprintError::InvalidArguments { reason })
                if reason == "Avalon bit_flags contains undefined source bits"
        ));
    }

    #[test]
    fn native_sub_byte_sizes_are_structured_errors() {
        for n_bits in [0, 1, 7] {
            let params = AvalonFingerprintParams {
                n_bits,
                ..Default::default()
            };
            assert!(matches!(
                avalon_fingerprint(&Molecule::new(), &params),
                Err(FingerprintError::InvalidArguments {
                    reason: "Avalon n_bits must be at least 8"
                })
            ));
        }
    }

    #[test]
    fn non_byte_aligned_sizes_remain_valid_source_arguments() {
        let molecule = Molecule::from_smiles("CC").expect("fixture");
        for n_bits in [9, 31, 32, 33, 511, 513] {
            let params = AvalonFingerprintParams {
                n_bits,
                ..Default::default()
            };
            let fingerprint = avalon_fingerprint(&molecule, &params).expect("source-valid size");
            assert_eq!(fingerprint.n_bits(), n_bits as usize);
        }
    }

    #[test]
    fn source_adapter_bits_match_native_profiles_and_second_pass() {
        let ethanol = Molecule::from_smiles("CCO").expect("fixture");
        let benzene = Molecule::from_smiles("c1ccccc1").expect("fixture");
        let default = AvalonFingerprintParams::default();
        let ethanol_fp = avalon_fingerprint(
            &ethanol,
            &AvalonFingerprintParams {
                n_bits: 64,
                ..default.clone()
            },
        )
        .expect("ethanol Avalon fingerprint");
        assert_eq!(ethanol_fp.on_bits(), vec![6, 14, 30, 31, 42]);
        let benzene_fp = avalon_fingerprint(
            &benzene,
            &AvalonFingerprintParams {
                n_bits: 64,
                ..default.clone()
            },
        )
        .expect("benzene Avalon fingerprint");
        assert_eq!(benzene_fp.on_bits(), vec![10, 16, 23, 29, 31, 37]);

        let non_sss = avalon_fingerprint(
            &ethanol,
            &AvalonFingerprintParams {
                n_bits: 64,
                bit_flags: AvalonFingerprintFlags::NON_SSS_BITS,
                ..default.clone()
            },
        )
        .expect("non-SSS fingerprint");
        assert_eq!(non_sss.on_bits(), vec![3]);

        let query = avalon_fingerprint(
            &Molecule::from_smiles("C[NH2+]C").unwrap(),
            &AvalonFingerprintParams {
                n_bits: 64,
                is_query: true,
                ..default.clone()
            },
        )
        .expect("query fingerprint");
        assert!(query.on_bits().is_empty());

        let narrow = avalon_fingerprint(
            &ethanol,
            &AvalonFingerprintParams {
                n_bits: 9,
                ..default
            },
        )
        .expect("rounded source size");
        assert_eq!(narrow.n_bits(), 9);
        // RDKit's adapter rounds its internal byte buffer to four bytes before
        // hashing, while the public vector retains the requested nine-bit size.
        assert_eq!(narrow.on_bits(), vec![6]);
        assert_eq!(ethanol.num_atoms(), 3);
    }

    #[test]
    fn adapter_reset_accumulate_repeat_and_size_boundaries_follow_source() {
        let molecule = Molecule::from_smiles("CCO").expect("fixture");
        for n_bits in [8, 31, 32, 33, 511, 512, 513] {
            let params = AvalonFingerprintParams {
                n_bits,
                ..Default::default()
            };
            let first = avalon_fingerprint(&molecule, &params).expect("valid source size");
            let second = avalon_fingerprint(&molecule, &params).expect("repeat is deterministic");
            assert_eq!(first, second);
            assert_eq!(first.n_bits(), n_bits as usize);
        }

        let mut state = super::reaccs::mol_to_reaccs(&molecule).expect("REACCS state");
        let mut bytes = vec![0_u8; 8];
        super::set_fingerprint_bits(
            &mut state,
            &mut bytes,
            AvalonFingerprintFlags::ALL_FEATURES,
            false,
            false,
        )
        .expect("first source pass");
        let first_pass = bytes.clone();
        super::set_fingerprint_bits(
            &mut state,
            &mut bytes,
            AvalonFingerprintFlags::ALL_FEATURES,
            false,
            true,
        )
        .expect("accumulated Daylight pass");
        for (before, after) in first_pass.iter().zip(bytes.iter()) {
            assert_eq!(after & before, *before);
        }

        let mut prefilled = vec![0xff_u8; 8];
        super::set_fingerprint_bits(
            &mut state,
            &mut prefilled,
            AvalonFingerprintFlags::ALL_FEATURES,
            false,
            false,
        )
        .expect("source SET_BIT does not clear caller storage");
        assert!(prefilled.iter().all(|byte| *byte == 0xff));
    }
}
