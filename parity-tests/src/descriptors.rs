//! Public descriptor calls and their explicit parameter spaces.
use crate::{Result, value::Value};
use cosmolkit::{Molecule, RotatableBondsOptions};
use serde_json::{Value as Json, json};

pub const TASKS: &[(&str, &str)] = &[
    ("labute_asa_sequence", "descriptors"),
    ("slogp_vsa_sequence", "descriptors"),
    ("smr_vsa_sequence", "descriptors"),
    ("chi_n_v_sequence", "descriptors"),
    ("chi_n_n_sequence", "descriptors"),
    ("slogp_vsa_bins", "descriptors"),
    ("smr_vsa_bins", "descriptors"),
    ("molecular_weight", "descriptors"),
    ("exact_molecular_weight", "descriptors"),
    ("molecular_formula", "descriptors"),
    ("num_heavy_atoms", "descriptors"),
    ("total_atom_count", "descriptors"),
    ("lipinski_hba", "descriptors"),
    ("lipinski_hbd", "descriptors"),
    ("fraction_csp3", "descriptors"),
    ("num_heteroatoms", "descriptors"),
    ("num_hba", "descriptors"),
    ("num_hbd", "descriptors"),
    ("num_amide_bonds", "descriptors"),
    ("num_atom_stereo_centers", "descriptors"),
    ("num_unspecified_atom_stereo_centers", "descriptors"),
    ("num_rotatable_bonds", "descriptors"),
    ("chi_0", "descriptors"),
    ("chi_1", "descriptors"),
    ("chi_0_n", "descriptors"),
    ("chi_1_n", "descriptors"),
    ("chi_2_n", "descriptors"),
    ("chi_3_n", "descriptors"),
    ("chi_4_n", "descriptors"),
    ("chi_0_v", "descriptors"),
    ("chi_1_v", "descriptors"),
    ("chi_2_v", "descriptors"),
    ("chi_3_v", "descriptors"),
    ("chi_4_v", "descriptors"),
    ("chi_n_n", "descriptors"),
    ("chi_n_v", "descriptors"),
    ("hall_kier_alpha", "descriptors"),
    ("kappa_1", "descriptors"),
    ("kappa_2", "descriptors"),
    ("kappa_3", "descriptors"),
    ("phi", "descriptors"),
    ("mqns", "descriptors"),
    ("crippen", "descriptors"),
    ("labute_asa", "descriptors"),
    ("tpsa", "descriptors"),
    ("slogp_vsa", "descriptors"),
    ("smr_vsa", "descriptors"),
    ("qed", "descriptors"),
    ("num_rings", "rings"),
    ("num_heterocycles", "rings"),
    ("num_aromatic_rings", "rings"),
    ("num_saturated_rings", "rings"),
    ("num_aliphatic_rings", "rings"),
    ("num_aromatic_heterocycles", "rings"),
    ("num_aromatic_carbocycles", "rings"),
    ("num_aliphatic_heterocycles", "rings"),
    ("num_aliphatic_carbocycles", "rings"),
    ("num_saturated_heterocycles", "rings"),
    ("num_saturated_carbocycles", "rings"),
    ("num_spiro_atoms", "rings"),
    ("num_bridgehead_atoms", "rings"),
];

pub fn profiles(task: &str) -> Result<Vec<Json>> {
    let profiles: Vec<Json> = match task {
        "chi_0" | "chi_1" => vec![json!({"result_type":"f64"})],
        "molecular_weight" | "exact_molecular_weight" => {
            [false, true].map(|v| json!({"only_heavy":v})).into()
        }
        "molecular_formula" => [false, true]
            .into_iter()
            .flat_map(|a| {
                [false, true].map(move |b| json!({"separate_isotopes":a,"abbreviate_h_isotopes":b}))
            })
            .collect(),
        "crippen" | "labute_asa" => [false, true]
            .into_iter()
            .flat_map(|h| [false, true].map(move |f| json!({"include_hydrogens":h,"force":f})))
            .collect(),
        "tpsa" => [false, true]
            .into_iter()
            .flat_map(|sp| {
                [false, true].map(move |f| json!({"include_sulfur_phosphorus":sp,"force":f}))
            })
            .collect(),
        "slogp_vsa" | "smr_vsa" => [
            Json::Null,
            json!([-2.0, 0.0, 2.0, 4.0]),
            json!([0.0, 0.0, 1.0]),
            json!([-0.2, 0.0, 0.25, 0.25, 0.8]),
        ]
        .into_iter()
        .flat_map(|bins| [false, true].map(move |f| json!({"bins":bins,"force":f})))
        .collect(),
        "labute_asa_sequence" => {
            vec![json!({"calls":[[false,false],[true,false],[true,true],[false,false]]})]
        }
        "slogp_vsa_sequence" | "smr_vsa_sequence" => {
            vec![json!({"custom_bins":[-0.2,0.0,0.25,0.25,0.8],"calls":"cold_warm_forced_v1"})]
        }
        "chi_n_v_sequence" | "chi_n_n_sequence" => {
            vec![json!({"forces":[false,false,true],"orders":[0,1,2,3,4,5,6]})]
        }
        "num_rotatable_bonds" => ["default", "non_strict", "strict", "strict_linkages"]
            .map(|v| json!({"definition":v}))
            .into(),
        "chi_n_n" | "chi_n_v" => (0..=6)
            .flat_map(|order| [false, true].map(move |force| json!({"order":order,"force":force})))
            .collect(),
        task if task.starts_with("chi_") && (task.ends_with("_n") || task.ends_with("_v"))
            || task == "mqns" =>
        {
            [false, true].map(|force| json!({"force":force})).into()
        }
        task if TASKS.iter().any(|(name, _)| *name == task) => vec![json!({})],
        _ => return Err(format!("unknown descriptor {task}").into()),
    };
    if matches!(
        task,
        "num_heavy_atoms"
            | "total_atom_count"
            | "lipinski_hba"
            | "lipinski_hbd"
            | "fraction_csp3"
            | "num_heteroatoms"
            | "num_hba"
            | "num_hbd"
    ) || TASKS.iter().any(|(name, suite)| {
        *name == task
            && *suite == "rings"
            && !matches!(task, "num_spiro_atoms" | "num_bridgehead_atoms")
    }) {
        return Ok(profiles
            .into_iter()
            .flat_map(|p| {
                [false, true].map(move |remove_hs| {
                    let mut p = p.clone();
                    p["remove_hs"] = json!(remove_hs);
                    p
                })
            })
            .collect());
    }
    Ok(profiles)
}

pub fn normalize(task: &str, p: Json) -> Result<Json> {
    let defaults = profiles(task)?;
    if task.ends_with("_sequence") {
        if defaults.contains(&p) {
            return Ok(p);
        }
        return Err("descriptor sequence must retain its complete ordered calls".into());
    }
    let template = defaults[0]
        .as_object()
        .ok_or("invalid descriptor defaults")?;
    let fields = p
        .as_object()
        .ok_or("descriptor parameters must be an object")?;
    if fields.len() != template.len() || fields.keys().any(|key| !template.contains_key(key)) {
        return Err(format!("invalid parameter fields for {task}").into());
    }
    for (key, value) in fields {
        let valid = match key.as_str() {
            "bins" => {
                value.is_null()
                    || value.as_array().is_some_and(|rows| {
                        rows.iter().all(|v| v.as_f64().is_some_and(f64::is_finite))
                    })
            }
            "order" => value.as_u64().is_some_and(|n| n <= u32::MAX as u64),
            "definition" => matches!(
                value.as_str(),
                Some("default" | "non_strict" | "strict" | "strict_linkages")
            ),
            "result_type" => value == "f64",
            _ => value.is_boolean(),
        };
        if !valid {
            return Err(format!("invalid {task}.{key}").into());
        }
    }
    Ok(p)
}

fn flag(p: &Json, name: &str) -> Result<bool> {
    p[name]
        .as_bool()
        .ok_or_else(|| format!("missing boolean {name}").into())
}
fn floats(values: impl IntoIterator<Item = f64>) -> Value {
    Value::List(
        values
            .into_iter()
            .map(|v| Value::Float(v.to_bits()))
            .collect(),
    )
}

pub fn scalar(m: &Molecule, task: &str, p: &Json) -> Result<Value> {
    if task == "labute_asa_sequence" {
        return Ok(floats(
            p["calls"]
                .as_array()
                .ok_or("missing ASA calls")?
                .iter()
                .map(|call| {
                    m.labute_asa_with_params(call[0].as_bool().unwrap(), call[1].as_bool().unwrap())
                })
                .collect::<std::result::Result<Vec<_>, _>>()?,
        ));
    }
    if matches!(task, "slogp_vsa_sequence" | "smr_vsa_sequence") {
        let bins: Vec<f64> = serde_json::from_value(p["custom_bins"].clone())?;
        let cold = m.slogp_vsa()?;
        let warm = m.slogp_vsa()?;
        let forced = m.slogp_vsa_with_params(Some(&bins), true)?;
        return Ok(Value::List(if task == "slogp_vsa_sequence" {
            vec![floats(cold), floats(warm), floats(forced)]
        } else {
            vec![
                floats(m.smr_vsa()?),
                floats(m.smr_vsa_with_params(Some(&bins), true)?),
            ]
        }));
    }
    if matches!(task, "chi_n_v_sequence" | "chi_n_n_sequence") {
        return Ok(Value::List(
            p["forces"]
                .as_array()
                .ok_or("missing forces")?
                .iter()
                .map(|force| {
                    let results = p["orders"]
                        .as_array()
                        .ok_or("missing orders")?
                        .iter()
                        .map(|order| {
                            if task == "chi_n_v_sequence" {
                                m.chi_n_v_with_params(
                                    order.as_u64().unwrap() as u32,
                                    force.as_bool().unwrap(),
                                )
                            } else {
                                m.chi_n_n_with_params(
                                    order.as_u64().unwrap() as u32,
                                    force.as_bool().unwrap(),
                                )
                            }
                        })
                        .collect::<std::result::Result<Vec<_>, _>>()?;
                    Ok(floats(results))
                })
                .collect::<Result<_>>()?,
        ));
    }
    if task == "slogp_vsa_bins" {
        return Ok(floats([
            m.slogp_vsa_1()?,
            m.slogp_vsa_2()?,
            m.slogp_vsa_3()?,
            m.slogp_vsa_4()?,
            m.slogp_vsa_5()?,
            m.slogp_vsa_6()?,
            m.slogp_vsa_7()?,
            m.slogp_vsa_8()?,
            m.slogp_vsa_9()?,
            m.slogp_vsa_10()?,
            m.slogp_vsa_11()?,
            m.slogp_vsa_12()?,
        ]));
    }
    if task == "smr_vsa_bins" {
        return Ok(floats([
            m.smr_vsa_1()?,
            m.smr_vsa_2()?,
            m.smr_vsa_3()?,
            m.smr_vsa_4()?,
            m.smr_vsa_5()?,
            m.smr_vsa_6()?,
            m.smr_vsa_7()?,
            m.smr_vsa_8()?,
            m.smr_vsa_9()?,
            m.smr_vsa_10()?,
        ]));
    }
    let force = if p.get("force").is_some() {
        flag(p, "force")?
    } else {
        false
    };
    let number = match task {
        "molecular_weight" => m.molecular_weight_with_params(flag(p, "only_heavy")?)?,
        "exact_molecular_weight" => m.exact_molecular_weight_with_params(flag(p, "only_heavy")?)?,
        "molecular_formula" => {
            return Ok(Value::Text(m.molecular_formula_with_params(
                flag(p, "separate_isotopes")?,
                flag(p, "abbreviate_h_isotopes")?,
            )?));
        }
        "chi_0" => m.chi_0()?,
        "chi_1" => m.chi_1()?,
        "chi_0_n" => m.chi_0_n_with_params(force)?,
        "chi_1_n" => m.chi_1_n_with_params(force)?,
        "chi_2_n" => m.chi_2_n_with_params(force)?,
        "chi_3_n" => m.chi_3_n_with_params(force)?,
        "chi_4_n" => m.chi_4_n_with_params(force)?,
        "chi_0_v" => m.chi_0_v_with_params(force)?,
        "chi_1_v" => m.chi_1_v_with_params(force)?,
        "chi_2_v" => m.chi_2_v_with_params(force)?,
        "chi_3_v" => m.chi_3_v_with_params(force)?,
        "chi_4_v" => m.chi_4_v_with_params(force)?,
        "chi_n_n" => {
            m.chi_n_n_with_params(p["order"].as_u64().ok_or("missing order")? as u32, force)?
        }
        "chi_n_v" => {
            m.chi_n_v_with_params(p["order"].as_u64().ok_or("missing order")? as u32, force)?
        }
        "fraction_csp3" => m.fraction_csp3()?,
        "kappa_1" => m.kappa_1()?,
        "kappa_2" => m.kappa_2()?,
        "kappa_3" => m.kappa_3()?,
        "phi" => m.phi()?,
        "hall_kier_alpha" => {
            let (sum, rows) = m.hall_kier_alpha_with_contributions()?;
            return Ok(Value::List(vec![Value::Float(sum.to_bits()), floats(rows)]));
        }
        "crippen" => {
            let v = m.crippen_descriptors_with_params(flag(p, "include_hydrogens")?, force)?;
            return Ok(floats([v.logp, v.molar_refractivity]));
        }
        "labute_asa" => {
            let h = flag(p, "include_hydrogens")?;
            let total = m.labute_asa_with_params(h, force)?;
            let v = m.labute_asa_contributions_with_params(h, force)?;
            return Ok(Value::List(vec![
                Value::Float(total.to_bits()),
                Value::Float(v.asa.to_bits()),
                floats(v.atom_contributions),
                Value::Float(v.hydrogen_contribution.to_bits()),
            ]));
        }
        "tpsa" => m.tpsa_with_params(flag(p, "include_sulfur_phosphorus")?, force)?,
        "slogp_vsa" | "smr_vsa" => {
            let bins: Option<Vec<f64>> = serde_json::from_value(p["bins"].clone())?;
            return Ok(floats(if task == "slogp_vsa" {
                m.slogp_vsa_with_params(bins.as_deref(), force)?
            } else {
                m.smr_vsa_with_params(bins.as_deref(), force)?
            }));
        }
        "qed" => m.qed()?,
        "mqns" => {
            return Ok(Value::List(
                m.mqns(force)?
                    .into_iter()
                    .map(|v| Value::Int(v as i64))
                    .collect(),
            ));
        }
        _ => {
            return Ok(Value::Int(match task {
                "num_heavy_atoms" => m.num_heavy_atoms()?,
                "total_atom_count" => m.total_atom_count()?,
                "lipinski_hba" => m.lipinski_hba()?,
                "lipinski_hbd" => m.lipinski_hbd()?,
                "num_heteroatoms" => m.num_heteroatoms()?,
                "num_hba" => m.num_hba()?,
                "num_hbd" => m.num_hbd()?,
                "num_amide_bonds" => m.num_amide_bonds()?,
                "num_atom_stereo_centers" => m.num_atom_stereo_centers()?,
                "num_unspecified_atom_stereo_centers" => m.num_unspecified_atom_stereo_centers()?,
                "num_rotatable_bonds" => {
                    m.num_rotatable_bonds_with_params(&match p["definition"]
                        .as_str()
                        .ok_or("missing definition")?
                    {
                        "default" => RotatableBondsOptions::Default,
                        "non_strict" => RotatableBondsOptions::NonStrict,
                        "strict" => RotatableBondsOptions::Strict,
                        "strict_linkages" => RotatableBondsOptions::StrictLinkages,
                        _ => return Err("invalid rotatable definition".into()),
                    })?
                }
                "num_rings" => m.num_rings()?,
                "num_heterocycles" => m.num_heterocycles()?,
                "num_aromatic_rings" => m.num_aromatic_rings()?,
                "num_saturated_rings" => m.num_saturated_rings()?,
                "num_aliphatic_rings" => m.num_aliphatic_rings()?,
                "num_aromatic_heterocycles" => m.num_aromatic_heterocycles()?,
                "num_aromatic_carbocycles" => m.num_aromatic_carbocycles()?,
                "num_aliphatic_heterocycles" => m.num_aliphatic_heterocycles()?,
                "num_aliphatic_carbocycles" => m.num_aliphatic_carbocycles()?,
                "num_saturated_heterocycles" => m.num_saturated_heterocycles()?,
                "num_saturated_carbocycles" => m.num_saturated_carbocycles()?,
                "num_spiro_atoms" => m.num_spiro_atoms()?,
                "num_bridgehead_atoms" => m.num_bridgehead_atoms()?,
                _ => return Err(format!("unknown descriptor {task}").into()),
            } as i64));
        }
    };
    Ok(Value::Float(number.to_bits()))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn all_descriptor_profiles_are_explicit_and_checked() {
        for (task, _) in TASKS {
            for p in profiles(task).unwrap() {
                assert_eq!(normalize(task, p.clone()).unwrap(), p);
            }
        }
        assert!(normalize("chi_n_n", json!({"force":false,"order":-1})).is_err());
        assert!(normalize("qed", json!({"ignored":true})).is_err());
        assert!(normalize("slogp_vsa", json!({"force":false,"bins":["bad"]})).is_err());
    }
}
