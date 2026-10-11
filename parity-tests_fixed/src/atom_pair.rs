//! Exact AtomPair observations through the existing public facade.
use crate::fingerprints::{Bits, CountVector};
use cosmolkit::{
    AtomPairFingerprintParams, AtomPairParams, AtomPairReadError, FingerprintAdditionalOutput,
    FingerprintError, Molecule,
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Form {
    SparseCount,
    Count,
    SparseBits,
    Bits,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Profile {
    pub form: Form,
    pub min_distance: u32,
    pub max_distance: u32,
    pub include_chirality: bool,
    pub use_2d: bool,
    pub count_simulation: bool,
    pub count_bounds: Vec<u32>,
    pub fp_size: u32,
    pub bits_per_feature: u32,
    pub from_first: bool,
    pub ignore_first: bool,
    pub index_plus_11: bool,
    pub additional_output: bool,
}

impl Profile {
    pub fn sparse_size(&self) -> u64 {
        // AtomPairGenerator.cpp getResultSize():
        // RDKit✔️✔️: OutputType result = 1;
        // RDKit✔️✔️: return (result << (numAtomPairFingerprintBits +
        // RDKit✔️✔️:                      2 * (this->dp_fingerprintArguments->df_includeChirality
        // RDKit✔️✔️:                               ? numChiralBits
        // RDKit✔️✔️:                               : 0)));
        // Representation size only; fingerprint chemistry stays in its owner.
        1u64 << (23 + if self.include_chirality { 4 } else { 0 })
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Additional {
    pub atom_counts: Vec<u32>,
    pub atom_to_bits: Vec<Vec<u64>>,
    pub bit_info_map: BTreeMap<u64, Vec<(u32, u32)>>,
    pub atoms_per_bit: BTreeMap<u64, Vec<Vec<i32>>>,
}

impl Additional {
    fn read(value: &FingerprintAdditionalOutput) -> Result<Self, String> {
        Ok(Self {
            atom_counts: value.atom_counts().ok_or("missing atom counts")?.to_vec(),
            atom_to_bits: value.atom_to_bits().ok_or("missing atom bits")?.to_vec(),
            bit_info_map: value.bit_info_map().ok_or("missing bit info")?.clone(),
            atoms_per_bit: value.atoms_per_bit().ok_or("missing bit atoms")?.clone(),
        })
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Observation {
    Counts {
        fingerprint: CountVector,
        additional_output: Option<Additional>,
    },
    Bits {
        fingerprint: Bits,
        additional_output: Option<Additional>,
    },
    /// Defined SparseIntVect bounds exception, not an upstream process crash.
    IndexError { index: u64, size: u64 },
}

pub fn run(mol: &Molecule, profile: &Profile) -> Result<Observation, String> {
    let n = mol.num_atoms();
    let first = || (0..n.min(1) as u32).collect();
    let params = AtomPairFingerprintParams {
        generator: AtomPairParams {
            min_distance: profile.min_distance,
            max_distance: profile.max_distance,
            include_chirality: profile.include_chirality,
            use_2d: profile.use_2d,
            count_simulation: profile.count_simulation,
            fp_size: profile.fp_size,
            count_bounds: profile.count_bounds.clone(),
            bits_per_feature: profile.bits_per_feature,
        },
        from_atoms: profile.from_first.then(first),
        ignore_atoms: profile.ignore_first.then(first),
        custom_atom_invariants: profile
            .index_plus_11
            .then(|| (0..n as u32).map(|i| i + 11).collect()),
        ..Default::default()
    };
    let mut additional = profile.additional_output.then(|| {
        let mut out = FingerprintAdditionalOutput::new();
        out.allocate_atom_counts();
        out.allocate_atom_to_bits();
        out.allocate_bit_info_map();
        out.allocate_atoms_per_bit();
        out
    });
    let result = match profile.form {
        Form::SparseCount => mol
            .fingerprint_atom_pair_sparse_count_with_params(&params, additional.as_mut())
            .map(|fp| Observation::Counts {
                fingerprint: CountVector {
                    length: fp.length(),
                    entries: fp
                        .nonzero_elements()
                        .iter()
                        .map(|(&i, &v)| (i, v))
                        .collect(),
                },
                additional_output: None,
            }),
        Form::Count => mol
            .fingerprint_atom_pair_count_with_params(&params, additional.as_mut())
            .map(|fp| Observation::Counts {
                fingerprint: CountVector {
                    length: u64::from(fp.length()),
                    entries: fp
                        .nonzero_elements()
                        .iter()
                        .map(|(&i, &v)| (u64::from(i), v))
                        .collect(),
                },
                additional_output: None,
            }),
        Form::SparseBits => mol
            .fingerprint_atom_pair_sparse_with_params(&params, additional.as_mut())
            .map(|fp| {
                let mut on_bits: Vec<u32> = fp.on_bits().into_iter().map(|i| i as u32).collect();
                on_bits.sort_unstable();
                Observation::Bits {
                    fingerprint: Bits {
                        length: fp.n_bits(),
                        on_bits,
                    },
                    additional_output: None,
                }
            }),
        Form::Bits => mol
            .fingerprint_atom_pair_with_params(&params, additional.as_mut())
            .map(|fp| Observation::Bits {
                fingerprint: fp.into(),
                additional_output: None,
            }),
    };
    let mut observed = match result {
        Ok(value) => value,
        Err(AtomPairReadError::Generator(error)) => {
            // The facade exposes the typed source chain, not the owner's enum.
            if let Some(FingerprintError::SparseIndexOutOfRange { index, size }) =
                std::error::Error::source(&error)
                    .and_then(|source| source.downcast_ref::<FingerprintError>())
                && profile.form == Form::SparseCount
                && profile.bits_per_feature > 1
                && *size == profile.sparse_size()
                && index >= size
            {
                return Ok(Observation::IndexError {
                    index: *index,
                    size: *size,
                });
            }
            return Err(error.to_string());
        }
        Err(error) => return Err(error.to_string()),
    };
    if let Some(value) = additional.as_ref() {
        let extra = Some(Additional::read(value)?);
        match &mut observed {
            Observation::Counts {
                additional_output, ..
            }
            | Observation::Bits {
                additional_output, ..
            } => *additional_output = extra,
            Observation::IndexError { .. } => unreachable!(),
        }
    }
    Ok(observed)
}

pub fn validate(profile: &Profile, value: &Observation) -> bool {
    let (length, additional) = match value {
        Observation::IndexError { index, size } => {
            return profile.form == Form::SparseCount
                && profile.bits_per_feature > 1
                && *size == profile.sparse_size()
                && index >= size;
        }
        Observation::Counts {
            fingerprint,
            additional_output,
        } => {
            if !matches!(profile.form, Form::Count | Form::SparseCount)
                || fingerprint
                    .entries
                    .iter()
                    .any(|&(i, v)| i >= fingerprint.length || v == 0)
                || !fingerprint.entries.windows(2).all(|w| w[0].0 < w[1].0)
            {
                return false;
            }
            (fingerprint.length, additional_output)
        }
        Observation::Bits {
            fingerprint,
            additional_output,
        } => {
            if !matches!(profile.form, Form::Bits | Form::SparseBits)
                || fingerprint.on_bits.iter().any(|&i| i >= fingerprint.length)
                || !fingerprint.on_bits.windows(2).all(|w| w[0] < w[1])
            {
                return false;
            }
            (u64::from(fingerprint.length), additional_output)
        }
    };
    let expected_length = if matches!(profile.form, Form::Count | Form::Bits) {
        u64::from(profile.fp_size)
    } else {
        profile.sparse_size()
    };
    length == expected_length
        && additional.is_some() == profile.additional_output
        && additional.as_ref().is_none_or(|out| {
            out.atom_counts.len() == out.atom_to_bits.len()
                && out.atom_to_bits.iter().flatten().all(|&i| i < length)
                && out
                    .bit_info_map
                    .keys()
                    .chain(out.atoms_per_bit.keys())
                    .all(|&i| i < length)
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn profiles() -> Vec<Profile> {
        let value: serde_json::Value = serde_json::from_str(include_str!(
            "../testdata/profiles/atom_pair_historical_0_3.json"
        ))
        .unwrap();
        value["profiles"]
            .as_array()
            .unwrap()
            .iter()
            .map(|row| serde_json::from_value(row["params"]["AtomPair"].clone()).unwrap())
            .collect()
    }

    #[test]
    fn atom_pair_every_historical_form_is_a_real_public_call() {
        let profiles = profiles();
        assert_eq!(profiles.len(), 40);
        for smiles in ["CCO", "C[C@H](N)C(=O)O", "c1ccccc1", ""] {
            let mol = Molecule::from_smiles(smiles).unwrap();
            for profile in &profiles {
                let observed = run(&mol, profile).unwrap();
                assert!(
                    validate(profile, &observed),
                    "{smiles} {profile:?} {observed:?}"
                );
                let extra = match &observed {
                    Observation::Counts {
                        additional_output, ..
                    }
                    | Observation::Bits {
                        additional_output, ..
                    } => additional_output,
                    Observation::IndexError { .. } => continue,
                };
                if let Some(extra) = extra {
                    assert_eq!(extra.atom_counts.len(), mol.num_atoms());
                    assert_eq!(extra.atom_to_bits.len(), mol.num_atoms());
                }
            }
        }
    }

    #[test]
    fn atom_pair_bounds_errors_keep_exact_index_and_size() {
        let profile = profiles()
            .into_iter()
            .find(|p| p.form == Form::SparseCount && p.bits_per_feature == 2)
            .unwrap();
        let actual = run(&Molecule::from_smiles("CCO").unwrap(), &profile).unwrap();
        let Observation::IndexError { index, size } = actual else {
            panic!("missing native bounds error")
        };
        assert!(index >= size);
        assert!(validate(&profile, &actual));
        assert!(!validate(
            &profile,
            &Observation::IndexError {
                index: size - 1,
                size
            }
        ));
        assert!(!validate(
            &profile,
            &Observation::IndexError {
                index,
                size: size + 1
            }
        ));
        let mut other = profile.clone();
        other.form = Form::Count;
        assert!(!validate(&other, &actual));
        other = profile;
        other.bits_per_feature = 1;
        assert!(!validate(&other, &actual));
    }
}
