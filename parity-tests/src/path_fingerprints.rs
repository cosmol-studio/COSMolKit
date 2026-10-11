//! Path, Pattern, Avalon and MACCS fingerprint calls and readable profile matrices.
use crate::{Result, fingerprints::bits, value::Value};
use cosmolkit::*;
use serde::{Deserialize, Serialize};
use serde_json::json;

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub struct TopologicalProfile {
    pub min_path: u32,
    pub max_path: u32,
    pub fp_size: u32,
    pub bits_per_feature: u32,
    pub use_hs: bool,
    pub target_density: f64,
    pub min_size: u32,
    pub branched_paths: bool,
    pub use_bond_order: bool,
    pub custom_atoms: bool,
    pub roots: String,
}
impl Default for TopologicalProfile {
    fn default() -> Self {
        Self {
            min_path: 1,
            max_path: 7,
            fp_size: 2048,
            bits_per_feature: 2,
            use_hs: true,
            target_density: 0.0,
            min_size: 128,
            branched_paths: true,
            use_bond_order: true,
            custom_atoms: false,
            roots: "all".into(),
        }
    }
}
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub struct LayeredProfile {
    pub layers: u32,
    pub min_path: u32,
    pub max_path: u32,
    pub fp_size: u32,
    pub branched_paths: bool,
    pub roots: String,
    pub counts: String,
    pub mask: String,
}
impl Default for LayeredProfile {
    fn default() -> Self {
        Self {
            layers: u32::MAX,
            min_path: 1,
            max_path: 7,
            fp_size: 2048,
            branched_paths: true,
            roots: "all".into(),
            counts: "absent".into(),
            mask: "absent".into(),
        }
    }
}
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub struct PatternProfile {
    pub fp_size: usize,
    pub tautomeric: bool,
    pub native_counts: String,
    pub native_mask: String,
}
impl Default for PatternProfile {
    fn default() -> Self {
        Self {
            fp_size: 2048,
            tautomeric: false,
            native_counts: "absent".into(),
            native_mask: "absent".into(),
        }
    }
}
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(default, deny_unknown_fields)]
pub struct AvalonProfile {
    pub fp_size: u32,
    pub is_query: bool,
    pub bit_flags: u32,
}

#[derive(Clone, Debug, Deserialize, Serialize, Default)]
#[serde(default, deny_unknown_fields)]
pub struct MaccsProfile {
    pub raw: bool,
}
impl Default for AvalonProfile {
    fn default() -> Self {
        Self {
            fp_size: 512,
            is_query: false,
            bit_flags: 32767,
        }
    }
}

fn matrix(base: serde_json::Value, changes: Vec<serde_json::Value>) -> Vec<serde_json::Value> {
    changes
        .into_iter()
        .map(|change| {
            let mut p = base.clone();
            p.as_object_mut()
                .unwrap()
                .extend(change.as_object().unwrap().clone());
            p
        })
        .collect()
}
pub fn profiles(task: &str) -> Result<Vec<serde_json::Value>> {
    Ok(match task {
        "fingerprint_topological" => matrix(
            serde_json::to_value(TopologicalProfile::default())?,
            vec![
                json!({}),
                json!({"fp_size":1024}),
                json!({"bits_per_feature":1}),
                json!({"min_path":2,"max_path":4}),
                json!({"use_hs":false}),
                json!({"branched_paths":false}),
                json!({"use_bond_order":false}),
                json!({"target_density":0.2}),
                json!({"roots":"first"}),
                json!({"custom_atoms":true}),
                json!({"max_path":3,"fp_size":1024,"bits_per_feature":1,"branched_paths":false,"use_bond_order":false}),
                json!({"min_path":2,"max_path":5,"fp_size":1024,"use_hs":false}),
                json!({"fp_size":4096,"target_density":0.35,"min_size":256}),
                json!({"branched_paths":false,"roots":"first"}),
            ],
        ),
        "fingerprint_layered" => matrix(
            serde_json::to_value(LayeredProfile::default())?,
            vec![
                json!({}),
                json!({"layers":63}),
                json!({"layers":1}),
                json!({"layers":2}),
                json!({"layers":4}),
                json!({"layers":8}),
                json!({"layers":16}),
                json!({"layers":32}),
                json!({"layers":7}),
                json!({"layers":63,"min_path":2,"max_path":4,"fp_size":257}),
                json!({"layers":63,"fp_size":64}),
                json!({"layers":63,"roots":"first"}),
                json!({"layers":63,"roots":"terminals"}),
                json!({"layers":63,"roots":"first","branched_paths":false}),
                json!({"layers":63,"fp_size":256,"counts":"zero","mask":"even"}),
                json!({"layers":63,"fp_size":257,"counts":"seeded","mask":"mod_three"}),
                json!({"layers":4294967232u32,"counts":"zero"}),
                json!({"layers":0,"counts":"zero"}),
            ],
        ),
        "fingerprint_pattern" => matrix(
            serde_json::to_value(PatternProfile::default())?,
            vec![
                json!({}),
                json!({"tautomeric":true}),
                json!({"fp_size":1}),
                json!({"fp_size":7,"tautomeric":true}),
                json!({"fp_size":127}),
                json!({"fp_size":257,"native_counts":"zero"}),
                json!({"fp_size":257,"tautomeric":true,"native_counts":"index_plus_11"}),
                json!({"fp_size":257,"native_mask":"empty"}),
                json!({"fp_size":257,"tautomeric":true,"native_mask":"all"}),
                json!({"fp_size":257,"native_counts":"index_plus_11","native_mask":"sparse"}),
            ],
        ),
        "fingerprint_avalon" => {
            let mut changes = vec![
                json!({}),
                json!({"bit_flags":15761407}),
                json!({"bit_flags":15728640}),
                json!({"is_query":true}),
            ];
            for bit in (0..15).chain(20..24) {
                changes.push(json!({"bit_flags":1u32<<bit}));
            }
            matrix(serde_json::to_value(AvalonProfile::default())?, changes)
        }
        "fingerprint_maccs" => vec![json!({"raw":false}), json!({"raw":true})],
        _ => return Err(format!("unknown path fingerprint {task}").into()),
    })
}
pub fn normalize(task: &str, profile: serde_json::Value) -> Result<serde_json::Value> {
    Ok(match task {
        "fingerprint_topological" => {
            let p: TopologicalProfile = serde_json::from_value(profile)?;
            roots(&p.roots, 1)?;
            serde_json::to_value(p)?
        }
        "fingerprint_layered" => {
            let p: LayeredProfile = serde_json::from_value(profile)?;
            roots(&p.roots, 1)?;
            if !["absent", "zero", "seeded"].contains(&p.counts.as_str())
                || !["absent", "even", "mod_three", "empty"].contains(&p.mask.as_str())
            {
                return Err("invalid Layered counts/mask".into());
            }
            serde_json::to_value(p)?
        }
        "fingerprint_pattern" => {
            let p: PatternProfile = serde_json::from_value(profile)?;
            if !["absent", "zero", "index_plus_11"].contains(&p.native_counts.as_str())
                || !["absent", "empty", "all", "sparse"].contains(&p.native_mask.as_str())
            {
                return Err("invalid native Pattern options".into());
            }
            serde_json::to_value(p)?
        }
        "fingerprint_avalon" => {
            serde_json::to_value(serde_json::from_value::<AvalonProfile>(profile)?)?
        }
        "fingerprint_maccs" => {
            serde_json::to_value(serde_json::from_value::<MaccsProfile>(profile)?)?
        }
        _ => return Err(format!("invalid {task} profile").into()),
    })
}
fn roots(mode: &str, atoms: usize) -> Result<Option<Vec<u32>>> {
    Ok(match mode {
        "all" => None,
        "first" => Some(if atoms == 0 { vec![] } else { vec![0] }),
        "terminals" => Some(match atoms {
            0 => vec![],
            1 => vec![0],
            _ => vec![0, atoms as u32 - 1],
        }),
        _ => return Err("roots must be all, first or terminals".into()),
    })
}
fn layered_params(p: &LayeredProfile, atoms: usize) -> Result<LayeredFingerprintParams> {
    Ok(LayeredFingerprintParams {
        layers: LayeredFingerprintLayers::from_bits_retain(p.layers),
        min_path: p.min_path,
        max_path: p.max_path,
        fp_size: p.fp_size,
        branched_paths: p.branched_paths,
        from_atoms: roots(&p.roots, atoms)?,
        atom_counts: match p.counts.as_str() {
            "absent" => None,
            "zero" => Some(vec![0; atoms]),
            "seeded" => Some((0..atoms as u32).map(|i| i + 10).collect()),
            _ => unreachable!(),
        },
        set_only_bits: if p.mask == "absent" {
            None
        } else {
            Some(Fingerprint::from_on_bits(
                p.fp_size,
                (0..p.fp_size).filter(|i| match p.mask.as_str() {
                    "even" => i % 2 == 0,
                    "mod_three" => i % 3 == 0,
                    _ => false,
                }),
            )?)
        },
    })
}
fn layered_value(result: LayeredFingerprintResult) -> Value {
    Value::List(vec![
        bits(&result.fingerprint),
        result.atom_counts.map_or(Value::Null, |v| {
            Value::List(v.into_iter().map(|i| Value::Int(i.into())).collect())
        }),
    ])
}
pub fn scalar(mol: &Molecule, task: &str, profile: &serde_json::Value) -> Result<Value> {
    Ok(match task {
        "fingerprint_topological" => {
            let p: TopologicalProfile = serde_json::from_value(profile.clone())?;
            let result = mol.fingerprint_topological_with_output_with_params(
                &TopologicalFingerprintParams {
                    min_path: p.min_path,
                    max_path: p.max_path,
                    fp_size: p.fp_size,
                    num_bits_per_feature: p.bits_per_feature,
                    use_hs: p.use_hs,
                    target_density: p.target_density,
                    min_size: p.min_size,
                    branched_paths: p.branched_paths,
                    use_bond_order: p.use_bond_order,
                    atom_invariants: p
                        .custom_atoms
                        .then(|| (1..=mol.num_atoms() as u32).collect()),
                    from_atoms: roots(&p.roots, mol.num_atoms())?,
                    ignore_atoms: None,
                },
                TopologicalFingerprintOutputRequest {
                    atom_bits: true,
                    bit_info: true,
                },
            )?;
            Value::List(vec![
                bits(&result.fingerprint),
                Value::List(
                    result
                        .atom_bits()?
                        .iter()
                        .map(|row| {
                            Value::List(row.iter().map(|i| Value::UInt((*i).into())).collect())
                        })
                        .collect(),
                ),
                Value::List(
                    result
                        .bit_info()?
                        .iter()
                        .map(|(key, rows)| {
                            Value::List(vec![
                                Value::UInt((*key).into()),
                                Value::List(
                                    rows.iter()
                                        .map(|row| {
                                            Value::List(
                                                row.iter()
                                                    .map(|i| Value::Int((*i).into()))
                                                    .collect(),
                                            )
                                        })
                                        .collect(),
                                ),
                            ])
                        })
                        .collect(),
                ),
            ])
        }
        "fingerprint_layered" => {
            let p: LayeredProfile = serde_json::from_value(profile.clone())?;
            layered_value(
                mol.fingerprint_layered_with_output_with_params(&layered_params(
                    &p,
                    mol.num_atoms(),
                )?)?,
            )
        }
        "fingerprint_pattern" => {
            let p: PatternProfile = serde_json::from_value(profile.clone())?;
            bits(
                &mol.fingerprint_pattern_with_params(&PatternFingerprintParams {
                    n_bits: p.fp_size,
                    tautomeric: p.tautomeric,
                })?,
            )
        }
        "fingerprint_avalon" => {
            let p: AvalonProfile = serde_json::from_value(profile.clone())?;
            bits(
                &mol.fingerprint_avalon_with_params(&AvalonFingerprintParams {
                    n_bits: p.fp_size,
                    is_query: p.is_query,
                    bit_flags: AvalonFingerprintFlags::from_bits_retain(p.bit_flags),
                })?,
            )
        }
        "fingerprint_maccs" => {
            let p: MaccsProfile = serde_json::from_value(profile.clone())?;
            bits(&if p.raw {
                mol.fingerprint_maccs_raw()?
            } else {
                mol.fingerprint_maccs()?
            })
        }
        _ => return Err(format!("unknown path fingerprint {task}").into()),
    })
}
pub fn has_batch(task: &str) -> bool {
    matches!(task, "fingerprint_layered" | "fingerprint_pattern")
}
pub fn batch(
    batch: &MoleculeBatch,
    task: &str,
    profile: &serde_json::Value,
    atoms: usize,
    runtime: &BatchQueryParams,
) -> std::result::Result<Vec<Option<Value>>, BatchValidationError> {
    if task == "fingerprint_pattern" {
        let p: PatternProfile =
            serde_json::from_value(profile.clone()).expect("validated Pattern profile");
        batch
            .fingerprint_pattern_list_with_params(
                &PatternFingerprintParams {
                    n_bits: p.fp_size,
                    tautomeric: p.tautomeric,
                },
                runtime,
            )
            .map(|v| v.into_iter().map(|v| v.map(|v| bits(&v))).collect())
    } else {
        let p: LayeredProfile =
            serde_json::from_value(profile.clone()).expect("validated Layered profile");
        // Structural validation of caller profiles occurs before any references are opened.
        let params = layered_params(&p, atoms).map_err(|error| {
            BatchValidationError::from_record_errors(vec![BatchError::with_source(
                0,
                "profile",
                std::io::Error::other(error),
            )])
        })?;
        batch
            .fingerprint_layered_with_output_list_with_params(&params, runtime)
            .map(|v| v.into_iter().map(|v| v.map(layered_value)).collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn historical_path_matrices_and_invalid_profiles() {
        for (task, count) in [
            ("fingerprint_topological", 14),
            ("fingerprint_layered", 18),
            ("fingerprint_pattern", 10),
            ("fingerprint_avalon", 23),
            ("fingerprint_maccs", 2),
        ] {
            let profiles = profiles(task).unwrap();
            assert_eq!(profiles.len(), count);
            for profile in profiles {
                normalize(task, profile).unwrap();
            }
        }
        assert!(normalize("fingerprint_layered", json!({"roots":"missing"})).is_err());
        assert!(normalize("fingerprint_pattern", json!({"fp_szie":1})).is_err());
    }
}
