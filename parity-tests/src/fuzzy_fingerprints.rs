//! Paired sparse-count operations. The second input is explicit in each profile.
use crate::{Result, fingerprints::counts, value::Value};
use cosmolkit::*;
use serde::{Deserialize, Serialize};
use serde_json::json;

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Profile {
    pub other_smiles: String,
    pub wide: bool,
    pub radius: u32,
    pub fp_size: u32,
}
pub fn profiles() -> Vec<serde_json::Value> {
    ["CCO", "c1ccccc1", "CC(=O)N"]
        .into_iter()
        .flat_map(|other| {
            [false, true].into_iter().map(
                move |wide| json!({"other_smiles":other,"wide":wide,"radius":2,"fp_size":2048}),
            )
        })
        .collect()
}
pub fn normalize(profile: serde_json::Value) -> Result<serde_json::Value> {
    let p: Profile = serde_json::from_value(profile)?;
    Molecule::from_smiles(&p.other_smiles)?;
    Ok(serde_json::to_value(p)?)
}
pub fn scalar(mol: &Molecule, profile: &serde_json::Value) -> Result<Value> {
    let p: Profile = serde_json::from_value(profile.clone())?;
    let other = Molecule::from_smiles(&p.other_smiles)?;
    let params = MorganFingerprintParams {
        generator: MorganParams {
            radius: p.radius,
            fp_size: p.fp_size,
            ..Default::default()
        },
        ..Default::default()
    };
    let left = mol.fingerprint_morgan_count_with_params(&params, None)?;
    let right = other.fingerprint_morgan_count_with_params(&params, None)?;
    if p.wide {
        let make = |source: &SparseCountFingerprint32| -> Result<SparseCountFingerprint> {
            let mut value = SparseCountFingerprint::new(1u64 << 40);
            for (key, count) in source.nonzero_elements() {
                value.set_value(u64::from(*key) + (1u64 << 32), *count)?;
            }
            Ok(value)
        };
        let left = make(&left)?;
        let right = make(&right)?;
        let and = left.fuzzy_and(&right)?;
        let or = left.fuzzy_or(&right)?;
        Ok(Value::List(vec![
            counts(
                and.length(),
                and.nonzero_elements().iter().map(|(k, v)| (*k, *v)),
            ),
            counts(
                or.length(),
                or.nonzero_elements().iter().map(|(k, v)| (*k, *v)),
            ),
        ]))
    } else {
        let and = left.fuzzy_and(&right)?;
        let or = left.fuzzy_or(&right)?;
        Ok(Value::List(vec![
            counts(
                and.length().into(),
                and.nonzero_elements()
                    .iter()
                    .map(|(k, v)| ((*k).into(), *v)),
            ),
            counts(
                or.length().into(),
                or.nonzero_elements().iter().map(|(k, v)| ((*k).into(), *v)),
            ),
        ]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn paired_inputs_and_both_index_widths_are_explicit() {
        let profiles = profiles();
        assert_eq!(profiles.len(), 6);
        for profile in profiles {
            normalize(profile).unwrap();
        }
    }
}
