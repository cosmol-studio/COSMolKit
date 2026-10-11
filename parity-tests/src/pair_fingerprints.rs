//! AtomPair and torsion generator calls. Output forms are compared separately.
use crate::{
    Result,
    fingerprints::{bits, counts, sparse_bits},
    value::Value,
};
use cosmolkit::*;
use serde::{Deserialize, Serialize};
use serde_json::json;

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub struct PairProfile {
    pub form: String,
    pub min_distance: u32,
    pub max_distance: u32,
    pub include_chirality: bool,
    pub use_2d: bool,
    pub count_simulation: bool,
    pub fp_size: u32,
    pub count_bounds: Vec<u32>,
    pub bits_per_feature: u32,
    pub from_first: bool,
    pub ignore_first: bool,
    pub custom_atoms: bool,
    pub additional_output: bool,
}
impl Default for PairProfile {
    fn default() -> Self {
        Self {
            form: "bits".into(),
            min_distance: 1,
            max_distance: 30,
            include_chirality: false,
            use_2d: true,
            count_simulation: true,
            fp_size: 2048,
            count_bounds: vec![1, 2, 4, 8],
            bits_per_feature: 1,
            from_first: false,
            ignore_first: false,
            custom_atoms: false,
            additional_output: false,
        }
    }
}
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub struct TorsionProfile {
    pub form: String,
    pub torsion_atom_count: u32,
    pub only_shortest_paths: bool,
    pub include_chirality: bool,
    pub count_simulation: bool,
    pub fp_size: u32,
    pub count_bounds: Vec<u32>,
    pub bits_per_feature: u32,
    pub from_first: bool,
    pub ignore_first: bool,
    pub custom_atoms: bool,
    pub additional_output: bool,
}
impl Default for TorsionProfile {
    fn default() -> Self {
        Self {
            form: "bits".into(),
            torsion_atom_count: 4,
            only_shortest_paths: false,
            include_chirality: false,
            count_simulation: true,
            fp_size: 2048,
            count_bounds: vec![1, 2, 4, 8],
            bits_per_feature: 1,
            from_first: false,
            ignore_first: false,
            custom_atoms: false,
            additional_output: false,
        }
    }
}

pub fn profiles(task: &str) -> Result<Vec<serde_json::Value>> {
    let (base, changes) = if task == "fingerprint_atom_pair" {
        (
            serde_json::to_value(PairProfile::default())?,
            vec![
                json!({}),
                json!({"include_chirality":true}),
                json!({"count_simulation":false}),
                json!({"count_bounds":[1,3,5],"fp_size":257}),
                json!({"min_distance":2,"max_distance":4}),
                json!({"from_first":true}),
                json!({"ignore_first":true}),
                json!({"custom_atoms":true}),
                json!({"bits_per_feature":2}),
                json!({"include_chirality":true,"additional_output":true}),
            ],
        )
    } else {
        (
            serde_json::to_value(TorsionProfile::default())?,
            vec![
                json!({}),
                json!({"count_simulation":false}),
                json!({"torsion_atom_count":3,"count_bounds":[1,3,5],"fp_size":1536}),
                json!({"include_chirality":true,"fp_size":4096}),
                json!({"torsion_atom_count":5,"only_shortest_paths":true,"fp_size":4096}),
                json!({"count_simulation":false,"fp_size":257,"bits_per_feature":3}),
                json!({"custom_atoms":true}),
                json!({"count_simulation":false,"fp_size":1000,"additional_output":true}),
                json!({"additional_output":true}),
            ],
        )
    };
    let mut profiles = Vec::new();
    for change in changes {
        for form in ["bits", "sparse_bits", "count", "sparse_count"] {
            let mut profile = base.clone();
            profile
                .as_object_mut()
                .unwrap()
                .extend(change.as_object().unwrap().clone());
            profile["form"] = json!(form);
            profiles.push(profile);
        }
    }
    Ok(profiles)
}
pub fn normalize(task: &str, profile: serde_json::Value) -> Result<serde_json::Value> {
    let value = if task == "fingerprint_atom_pair" {
        serde_json::to_value(serde_json::from_value::<PairProfile>(profile)?)?
    } else {
        serde_json::to_value(serde_json::from_value::<TorsionProfile>(profile)?)?
    };
    if !["bits", "sparse_bits", "count", "sparse_count"].contains(&value["form"].as_str().unwrap())
    {
        return Err("unknown fingerprint form".into());
    }
    Ok(value)
}
impl PairProfile {
    fn params(&self, atoms: usize) -> AtomPairFingerprintParams {
        AtomPairFingerprintParams {
            generator: AtomPairParams {
                min_distance: self.min_distance,
                max_distance: self.max_distance,
                include_chirality: self.include_chirality,
                use_2d: self.use_2d,
                count_simulation: self.count_simulation,
                fp_size: self.fp_size,
                count_bounds: self.count_bounds.clone(),
                bits_per_feature: self.bits_per_feature,
            },
            from_atoms: self
                .from_first
                .then(|| if atoms == 0 { vec![] } else { vec![0] }),
            ignore_atoms: self
                .ignore_first
                .then(|| if atoms == 0 { vec![] } else { vec![0] }),
            custom_atom_invariants: self
                .custom_atoms
                .then(|| (0..atoms as u32).map(|i| i + 11).collect()),
            ..Default::default()
        }
    }
}
impl TorsionProfile {
    fn params(&self, atoms: usize) -> TopologicalTorsionFingerprintParams {
        TopologicalTorsionFingerprintParams {
            generator: TopologicalTorsionParams {
                torsion_atom_count: self.torsion_atom_count,
                only_shortest_paths: self.only_shortest_paths,
                include_chirality: self.include_chirality,
                count_simulation: self.count_simulation,
                fp_size: self.fp_size,
                count_bounds: self.count_bounds.clone(),
                bits_per_feature: self.bits_per_feature,
            },
            from_atoms: self
                .from_first
                .then(|| if atoms == 0 { vec![] } else { vec![0] }),
            ignore_atoms: self
                .ignore_first
                .then(|| if atoms == 0 { vec![] } else { vec![0] }),
            custom_atom_invariants: self
                .custom_atoms
                .then(|| (0..atoms as u32).map(|i| i + 17).collect()),
            ..Default::default()
        }
    }
}
fn extra_output(collect: bool) -> FingerprintAdditionalOutput {
    let mut output = FingerprintAdditionalOutput::new();
    if collect {
        output.allocate_atom_counts();
        output.allocate_atom_to_bits();
        output.allocate_bit_info_map();
        output.allocate_bit_paths();
        output.allocate_atoms_per_bit();
    }
    output
}
fn extra_value(output: &FingerprintAdditionalOutput) -> Value {
    let mut fields = match crate::fingerprints::provenance(
        output.atom_counts(),
        output.atom_to_bits(),
        output.bit_info_map(),
        output.atoms_per_bit(),
    ) {
        Value::List(fields) => fields,
        _ => unreachable!(),
    };
    fields.push(output.bit_paths().map_or(Value::Null, |entries| {
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
                                        row.iter().map(|i| Value::Int((*i).into())).collect(),
                                    )
                                })
                                .collect(),
                        ),
                    ])
                })
                .collect(),
        )
    }));
    Value::List(fields)
}
pub fn scalar(mol: &Molecule, task: &str, profile: &serde_json::Value) -> Result<Value> {
    let collect = profile["additional_output"].as_bool().unwrap();
    let mut output = extra_output(collect);
    let extra = collect.then_some(&mut output);
    let value = if task == "fingerprint_atom_pair" {
        let p: PairProfile = serde_json::from_value(profile.clone())?;
        let params = p.params(mol.num_atoms());
        match p.form.as_str() {
            "bits" => bits(&mol.fingerprint_atom_pair_with_params(&params, extra)?),
            "sparse_bits" => {
                sparse_bits(&mol.fingerprint_atom_pair_sparse_with_params(&params, extra)?)
            }
            "count" => {
                let fp = mol.fingerprint_atom_pair_count_with_params(&params, extra)?;
                counts(
                    fp.length().into(),
                    fp.nonzero_elements().iter().map(|(k, v)| ((*k).into(), *v)),
                )
            }
            "sparse_count" => {
                let fp = mol.fingerprint_atom_pair_sparse_count_with_params(&params, extra)?;
                counts(
                    fp.length(),
                    fp.nonzero_elements().iter().map(|(k, v)| (*k, *v)),
                )
            }
            _ => return Err("unknown AtomPair form".into()),
        }
    } else {
        let p: TorsionProfile = serde_json::from_value(profile.clone())?;
        let params = p.params(mol.num_atoms());
        match p.form.as_str() {
            "bits" => bits(&mol.fingerprint_topological_torsion_with_params(&params, extra)?),
            "sparse_bits" => sparse_bits(
                &mol.fingerprint_topological_torsion_sparse_with_params(&params, extra)?,
            ),
            "count" => {
                let fp = mol.fingerprint_topological_torsion_count_with_params(&params, extra)?;
                counts(
                    fp.length().into(),
                    fp.nonzero_elements().iter().map(|(k, v)| ((*k).into(), *v)),
                )
            }
            "sparse_count" => {
                let fp =
                    mol.fingerprint_topological_torsion_sparse_count_with_params(&params, extra)?;
                counts(
                    fp.length(),
                    fp.nonzero_elements().iter().map(|(k, v)| (*k, *v)),
                )
            }
            _ => return Err("unknown torsion form".into()),
        }
    };
    Ok(Value::List(vec![
        value,
        if collect {
            extra_value(&output)
        } else {
            Value::Null
        },
    ]))
}
pub fn has_batch(task: &str, profile: &serde_json::Value) -> bool {
    !profile["additional_output"].as_bool().unwrap()
        && (task == "fingerprint_atom_pair" || profile["form"] == "bits")
}
pub fn batch(
    batch: &MoleculeBatch,
    task: &str,
    profile: &serde_json::Value,
    atoms: usize,
    runtime: &BatchQueryParams,
) -> std::result::Result<Vec<Option<Value>>, BatchValidationError> {
    let values = if task == "fingerprint_atom_pair" {
        let p: PairProfile =
            serde_json::from_value(profile.clone()).expect("validated pair profile");
        let p2 = p.params(atoms);
        match p.form.as_str() {
            "bits" => batch
                .fingerprint_atom_pair_list_with_params(&p2, runtime)
                .map(|v| v.into_iter().map(|v| v.map(|v| bits(&v))).collect()),
            "sparse_bits" => batch
                .fingerprint_atom_pair_sparse_bits_list_with_params(&p2, runtime)
                .map(|v| v.into_iter().map(|v| v.map(|v| sparse_bits(&v))).collect()),
            "count" => batch
                .fingerprint_atom_pair_count_list_with_params(&p2, runtime)
                .map(|v| {
                    v.into_iter()
                        .map(|v| {
                            v.map(|v| {
                                counts(
                                    v.length().into(),
                                    v.nonzero_elements().iter().map(|(k, c)| ((*k).into(), *c)),
                                )
                            })
                        })
                        .collect()
                }),
            "sparse_count" => batch
                .fingerprint_atom_pair_sparse_count_list_with_params(&p2, runtime)
                .map(|v| {
                    v.into_iter()
                        .map(|v| {
                            v.map(|v| {
                                counts(
                                    v.length(),
                                    v.nonzero_elements().iter().map(|(k, c)| (*k, *c)),
                                )
                            })
                        })
                        .collect()
                }),
            _ => unreachable!(),
        }
    } else {
        let p: TorsionProfile =
            serde_json::from_value(profile.clone()).expect("validated torsion profile");
        batch
            .fingerprint_topological_torsion_list_with_params(&p.params(atoms), runtime)
            .map(|v| v.into_iter().map(|v| v.map(|v| bits(&v))).collect())
    };
    values.map(|values: Vec<Option<Value>>| {
        values
            .into_iter()
            .map(|v| v.map(|v| Value::List(vec![v, Value::Null])))
            .collect()
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn retained_pair_and_torsion_matrices() {
        for (task, count) in [("fingerprint_atom_pair", 40), ("fingerprint_torsion", 36)] {
            let profiles = profiles(task).unwrap();
            assert_eq!(profiles.len(), count);
            for profile in profiles {
                normalize(task, profile).unwrap();
            }
        }
    }
}
