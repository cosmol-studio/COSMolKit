//! Fingerprint profiles and explicit public calls, without domain algorithms.
use crate::{Result, value::Value};
use cosmolkit::*;
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::collections::BTreeMap;

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub struct MorganProfile {
    pub form: String,
    pub radius: u32,
    pub include_chirality: bool,
    pub use_bond_types: bool,
    pub include_ring_membership: bool,
    pub only_nonzero_invariants: bool,
    pub include_redundant_environments: bool,
    pub fp_size: u32,
    pub count_simulation: bool,
    pub count_bounds: Vec<u32>,
    pub bits_per_feature: u32,
    pub feature_invariants: bool,
    pub bond_chirality: Option<bool>,
    pub from_first: bool,
    pub ignore_first: bool,
    pub custom_atoms: bool,
    pub custom_bonds: bool,
    pub custom_atom_offset: u32,
    pub zero_every: usize,
    pub custom_bond_offset: u32,
    pub additional_output: bool,
}
impl Default for MorganProfile {
    fn default() -> Self {
        Self {
            form: "bits".into(),
            radius: 2,
            include_chirality: false,
            use_bond_types: true,
            include_ring_membership: true,
            only_nonzero_invariants: false,
            include_redundant_environments: false,
            fp_size: 2048,
            count_simulation: false,
            count_bounds: vec![1, 2, 4, 8],
            bits_per_feature: 1,
            feature_invariants: false,
            bond_chirality: None,
            from_first: false,
            ignore_first: false,
            custom_atoms: false,
            custom_bonds: false,
            custom_atom_offset: 1,
            zero_every: 2,
            custom_bond_offset: 7,
            additional_output: true,
        }
    }
}

pub fn morgan_profiles() -> Result<Vec<serde_json::Value>> {
    // Historical 15 generator/call recipes, each in all four output forms.
    // Explicit overrides make the parameter matrix readable and editable.
    let changes = [
        json!({}),
        json!({"include_chirality":true}),
        json!({"radius":3}),
        json!({"use_bond_types":false,"fp_size":1024}),
        json!({"count_simulation":true}),
        json!({"include_redundant_environments":true}),
        json!({"feature_invariants":true}),
        json!({"include_ring_membership":false}),
        json!({"bond_chirality":true}),
        json!({"bits_per_feature":2}),
        json!({"from_first":true}),
        json!({"ignore_first":true}),
        json!({"custom_atoms":true}),
        json!({"custom_atoms":true,"only_nonzero_invariants":true}),
        json!({"custom_bonds":true}),
    ];
    let mut profiles = Vec::new();
    for changes in changes {
        for form in ["bits", "sparse_bits", "count", "sparse_count"] {
            let mut profile = serde_json::to_value(MorganProfile::default())?;
            profile
                .as_object_mut()
                .unwrap()
                .extend(changes.as_object().unwrap().clone());
            profile["form"] = json!(form);
            profiles.push(profile);
        }
    }
    Ok(profiles)
}

impl MorganProfile {
    fn generator_params(&self) -> MorganGeneratorParams {
        MorganGeneratorParams {
            radius: self.radius,
            include_chirality: self.include_chirality,
            use_bond_types: self.use_bond_types,
            include_ring_membership: self.include_ring_membership,
            only_nonzero_invariants: self.only_nonzero_invariants,
            include_redundant_environments: self.include_redundant_environments,
            fp_size: self.fp_size,
            count_simulation: self.count_simulation,
            count_bounds: self.count_bounds.clone(),
            bits_per_feature: self.bits_per_feature,
        }
    }
    fn atom_generator(&self) -> MorganAtomInvariantsGenerator {
        if self.feature_invariants {
            MorganAtomInvariantsGenerator::features(None)
        } else {
            MorganAtomInvariantsGenerator::connectivity(self.include_ring_membership)
        }
    }
    fn bond_generator(&self) -> Option<MorganBondInvariantsGenerator> {
        self.bond_chirality
            .map(|chirality| MorganBondInvariantsGenerator::new(self.use_bond_types, chirality))
    }
    fn call(&self, atoms: usize, bonds: usize) -> MorganCallParams {
        MorganCallParams::new(
            self.from_first
                .then(|| if atoms == 0 { vec![] } else { vec![0] }),
            self.ignore_first
                .then(|| if atoms == 0 { vec![] } else { vec![0] }),
            self.custom_atoms.then(|| {
                (0..atoms)
                    .map(|index| {
                        if self.only_nonzero_invariants && index % self.zero_every == 0 {
                            0
                        } else {
                            index as u32 + self.custom_atom_offset
                        }
                    })
                    .collect()
            }),
            self.custom_bonds.then(|| {
                (0..bonds)
                    .map(|index| index as u32 + self.custom_bond_offset)
                    .collect()
            }),
            -1,
        )
    }
}

pub fn bits(fp: &Fingerprint) -> Value {
    Value::List(vec![
        Value::UInt(fp.n_bits().into()),
        Value::List(
            fp.on_bits()
                .into_iter()
                .map(|bit| Value::Int(bit.into()))
                .collect(),
        ),
    ])
}
pub fn sparse_bits(fp: &SparseBitFingerprint) -> Value {
    Value::List(vec![
        Value::UInt(fp.n_bits().into()),
        Value::List(
            fp.on_bits()
                .into_iter()
                .map(|bit| Value::Int(bit.into()))
                .collect(),
        ),
    ])
}
pub fn counts(length: u64, entries: impl Iterator<Item = (u64, i32)>) -> Value {
    Value::List(vec![
        Value::UInt(length),
        Value::List(
            entries
                .map(|(bit, count)| Value::List(vec![Value::UInt(bit), Value::Int(count.into())]))
                .collect(),
        ),
    ])
}
pub(crate) fn provenance(
    counts: Option<&[u32]>,
    atoms: Option<&[Vec<u64>]>,
    info: Option<&BTreeMap<u64, Vec<(u32, u32)>>>,
    per_bit: Option<&BTreeMap<u64, Vec<Vec<i32>>>>,
) -> Value {
    Value::List(vec![
        counts.map_or(Value::Null, |values| {
            Value::List(
                values
                    .iter()
                    .map(|value| Value::Int((*value).into()))
                    .collect(),
            )
        }),
        atoms.map_or(Value::Null, |rows| {
            Value::List(
                rows.iter()
                    .map(|row| Value::List(row.iter().copied().map(Value::UInt).collect()))
                    .collect(),
            )
        }),
        info.map_or(Value::Null, |entries| {
            Value::List(
                entries
                    .iter()
                    .map(|(key, rows)| {
                        Value::List(vec![
                            Value::UInt(*key),
                            Value::List(
                                rows.iter()
                                    .map(|(atom, radius)| {
                                        Value::List(vec![
                                            Value::Int((*atom).into()),
                                            Value::Int((*radius).into()),
                                        ])
                                    })
                                    .collect(),
                            ),
                        ])
                    })
                    .collect(),
            )
        }),
        per_bit.map_or(Value::Null, |entries| {
            Value::List(
                entries
                    .iter()
                    .map(|(key, rows)| {
                        Value::List(vec![
                            Value::UInt(*key),
                            Value::List(
                                rows.iter()
                                    .map(|row| {
                                        Value::List(
                                            row.iter()
                                                .map(|atom| Value::Int((*atom).into()))
                                                .collect(),
                                        )
                                    })
                                    .collect(),
                            ),
                        ])
                    })
                    .collect(),
            )
        }),
    ])
}

pub fn morgan(mol: &Molecule, profile: &serde_json::Value) -> Result<Value> {
    let profile: MorganProfile = serde_json::from_value(profile.clone())?;
    let generator = MorganFingerprintGenerator::new(
        Some(&profile.generator_params()),
        Some(&profile.atom_generator()),
        profile.bond_generator().as_ref(),
    )?;
    let call = profile.call(mol.num_atoms(), mol.num_bonds());
    let mut output = FingerprintAdditionalOutput::new();
    if profile.additional_output {
        output.allocate_atom_counts();
        output.allocate_atom_to_bits();
        output.allocate_bit_info_map();
        output.allocate_atoms_per_bit();
    }
    let extra = profile.additional_output.then_some(&mut output);
    let fp = match profile.form.as_str() {
        "bits" => bits(&mol.fingerprint_morgan_with_generator(&generator, Some(&call), extra)?),
        "sparse_bits" => sparse_bits(&mol.fingerprint_morgan_sparse_with_generator(
            &generator,
            Some(&call),
            extra,
        )?),
        "count" => {
            let fp = mol.fingerprint_morgan_count_with_generator(&generator, Some(&call), extra)?;
            counts(
                fp.length().into(),
                fp.nonzero_elements()
                    .iter()
                    .map(|(key, count)| ((*key).into(), *count)),
            )
        }
        "sparse_count" => {
            let fp =
                mol.fingerprint_morgan_sparse_count_with_generator(&generator, Some(&call), extra)?;
            counts(
                fp.length(),
                fp.nonzero_elements()
                    .iter()
                    .map(|(key, count)| (*key, *count)),
            )
        }
        _ => return Err("unknown Morgan form".into()),
    };
    Ok(Value::List(vec![
        fp,
        if profile.additional_output {
            provenance(
                output.atom_counts(),
                output.atom_to_bits(),
                output.bit_info_map(),
                output.atoms_per_bit(),
            )
        } else {
            Value::Null
        },
    ]))
}

pub fn morgan_has_batch(profile: &serde_json::Value) -> bool {
    profile
        .get("form")
        .and_then(serde_json::Value::as_str)
        .unwrap_or("bits")
        == "bits"
}

pub fn morgan_batch(
    batch: &MoleculeBatch,
    profile: &serde_json::Value,
    atoms: usize,
    bonds: usize,
    runtime: &BatchQueryParams,
) -> std::result::Result<Vec<Option<Value>>, BatchValidationError> {
    // Profiles are validated before execution, never silently defaulted here.
    let profile: MorganProfile =
        serde_json::from_value(profile.clone()).expect("validated Morgan profile");
    batch
        .fingerprint_morgan_with_output_with_generator_params(
            &profile.generator_params(),
            Some(&profile.atom_generator()),
            profile.bond_generator().as_ref(),
            &profile.call(atoms, bonds),
            profile.additional_output,
            runtime,
        )
        .map(|values| {
            values
                .into_iter()
                .map(|value| {
                    value.map(|value| {
                        let extra = value
                            .additional_output
                            .as_ref()
                            .map_or(Value::Null, |extra| {
                                provenance(
                                    extra.atom_counts(),
                                    extra.atom_to_bits(),
                                    extra.bit_info_map(),
                                    extra.atoms_per_bit(),
                                )
                            });
                        Value::List(vec![bits(&value.fingerprint), extra])
                    })
                })
                .collect()
        })
}

pub fn validate_morgan(profile: &serde_json::Value) -> Result<()> {
    let profile: MorganProfile = serde_json::from_value(profile.clone())?;
    if !["bits", "sparse_bits", "count", "sparse_count"].contains(&profile.form.as_str()) {
        return Err("unknown Morgan form".into());
    }
    if profile.zero_every == 0 {
        return Err("zero_every must be positive".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn morgan_matrix_has_four_forms_and_all_historical_recipes() {
        let profiles = morgan_profiles().unwrap();
        assert_eq!(profiles.len(), 60);
        for profile in profiles {
            validate_morgan(&profile).unwrap();
        }
        let mut p = MorganProfile::default();
        p.custom_atoms = true;
        p.custom_bonds = true;
        p.only_nonzero_invariants = true;
        let call = p.call(4, 3);
        assert_eq!(call.custom_atom_invariants, Some(vec![0, 2, 0, 4]));
        assert_eq!(call.custom_bond_invariants, Some(vec![7, 8, 9]));
        assert!(validate_morgan(&json!({"raduis":2})).is_err());
        assert!(validate_morgan(&json!({"form":"missing"})).is_err());
    }
}
