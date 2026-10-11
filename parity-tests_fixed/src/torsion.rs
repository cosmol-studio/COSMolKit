//! Exact TopologicalTorsion observations through the existing public facade.
use crate::fingerprints::{Bits, CountVector};
use cosmolkit::{
    FingerprintAdditionalOutput, Molecule, TopologicalTorsionFingerprintParams,
    TopologicalTorsionParams,
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
    pub torsion_atom_count: u32,
    pub only_shortest_paths: bool,
    pub include_chirality: bool,
    pub count_simulation: bool,
    pub count_bounds: Vec<u32>,
    pub fp_size: u32,
    pub bits_per_feature: u32,
    pub from_first: bool,
    pub ignore_first: bool,
    pub index_plus_17: bool,
    pub additional_output: bool,
}

impl Profile {
    pub fn sparse_size(&self) -> u64 {
        // TopologicalTorsionGenerator.cpp getResultSize():
        // RDKit✔️✔️:   OutputType result = 1;
        // RDKit✔️✔️:   return (result << ((
        // RDKit✔️✔️:               dynamic_cast<const TopologicalTorsionArguments *>(
        // RDKit✔️✔️:                   this->dp_fingerprintArguments)
        // RDKit✔️✔️:                   ->d_torsionAtomCount *
        // RDKit✔️✔️:               (codeSize + (dynamic_cast<const TopologicalTorsionArguments *>(
        // RDKit✔️✔️:                                this->dp_fingerprintArguments)
        // RDKit✔️✔️:                                    ->df_includeChirality
        // RDKit✔️✔️:                                ? numChiralBits
        // RDKit✔️✔️:                                : 0)))));
        // Representation only: codeSize=9 and numChiralBits=2; chemistry
        // and allocation remain owned by the public fingerprint operation.
        1u64 << (self.torsion_atom_count * (9 + if self.include_chirality { 2 } else { 0 }))
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Additional {
    pub bit_paths: BTreeMap<u64, Vec<Vec<i32>>>,
    pub atom_counts: Vec<u32>,
    pub atom_to_bits: Vec<Vec<u64>>,
    pub bit_info_map: BTreeMap<u64, Vec<(u32, u32)>>,
    pub atoms_per_bit: BTreeMap<u64, Vec<Vec<i32>>>,
}

impl Additional {
    fn read(value: &FingerprintAdditionalOutput) -> Result<Self, String> {
        Ok(Self {
            bit_paths: value.bit_paths().ok_or("missing bit paths")?.clone(),
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
}

pub fn run(mol: &Molecule, profile: &Profile) -> Result<Observation, String> {
    let n = mol.num_atoms();
    let first = || (0..n.min(1) as u32).collect();
    let params = TopologicalTorsionFingerprintParams {
        generator: TopologicalTorsionParams {
            torsion_atom_count: profile.torsion_atom_count,
            only_shortest_paths: profile.only_shortest_paths,
            include_chirality: profile.include_chirality,
            count_simulation: profile.count_simulation,
            fp_size: profile.fp_size,
            count_bounds: profile.count_bounds.clone(),
            bits_per_feature: profile.bits_per_feature,
        },
        from_atoms: profile.from_first.then(first),
        ignore_atoms: profile.ignore_first.then(first),
        custom_atom_invariants: profile
            .index_plus_17
            .then(|| (0..n as u32).map(|i| i + 17).collect()),
        ..Default::default()
    };
    let mut additional = profile.additional_output.then(|| {
        let mut out = FingerprintAdditionalOutput::new();
        out.allocate_bit_paths();
        out.allocate_atom_counts();
        out.allocate_atom_to_bits();
        out.allocate_bit_info_map();
        out.allocate_atoms_per_bit();
        out
    });
    let result = match profile.form {
        Form::SparseCount => mol
            .fingerprint_topological_torsion_sparse_count_with_params(&params, additional.as_mut())
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
            .fingerprint_topological_torsion_count_with_params(&params, additional.as_mut())
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
            .fingerprint_topological_torsion_sparse_with_params(&params, additional.as_mut())
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
            .fingerprint_topological_torsion_with_params(&params, additional.as_mut())
            .map(|fp| Observation::Bits {
                fingerprint: fp.into(),
                additional_output: None,
            }),
    };
    let mut observed = result.map_err(|error| error.to_string())?;
    if let Some(value) = additional.as_ref() {
        let extra = Some(Additional::read(value)?);
        match &mut observed {
            Observation::Counts {
                additional_output, ..
            }
            | Observation::Bits {
                additional_output, ..
            } => *additional_output = extra,
        }
    }
    Ok(observed)
}

pub fn validate(profile: &Profile, value: &Observation) -> bool {
    let (length, additional) = match value {
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
        if profile.form == Form::SparseBits {
            profile.sparse_size().min(u64::from(u32::MAX))
        } else {
            profile.sparse_size()
        }
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
                    .chain(out.bit_paths.keys())
                    .all(|&i| i < length)
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn profiles() -> Vec<Profile> {
        let value: serde_json::Value = serde_json::from_str(include_str!(
            "../testdata/profiles/torsion_historical_0_3.json"
        ))
        .unwrap();
        value["profiles"]
            .as_array()
            .unwrap()
            .iter()
            .map(|row| serde_json::from_value(row["params"]["TopologicalTorsion"].clone()).unwrap())
            .collect()
    }

    #[test]
    fn torsion_every_historical_form_is_a_real_public_call() {
        let profiles = profiles();
        assert_eq!(profiles.len(), 36);
        for smiles in ["CCCC", "C[C@H](N)C(=O)O", "c1ccccc1", ""] {
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
                };
                if let Some(extra) = extra {
                    assert_eq!(extra.atom_counts.len(), mol.num_atoms());
                    assert_eq!(extra.atom_to_bits.len(), mol.num_atoms());
                }
            }
        }
    }
}
