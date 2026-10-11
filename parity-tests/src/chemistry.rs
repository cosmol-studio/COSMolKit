//! Read state, complete transforms and matrix queries through public APIs.
use crate::{Result, tasks, value::Value};
use cosmolkit::{BatchErrorMode, BatchParams, BatchRecord, Molecule, MoleculeBatch};
use rayon::prelude::*;
use serde_json::{Value as Json, json};

pub const TASKS: &[(&str, &str)] = &[
    ("smiles_read", "core"),
    ("sanitize", "core"),
    ("kekulize", "core"),
    ("add_hs", "core"),
    ("remove_hs", "core"),
    ("potential_stereo", "core"),
    ("distance_matrix", "core"),
    ("fragments", "scaffolds"),
    ("largest_fragment", "scaffolds"),
    ("murcko_scaffold", "scaffolds"),
    ("net_scaffold", "scaffolds"),
    ("murcko_decompose", "scaffolds"),
];

pub fn profiles(task: &str) -> Result<Vec<Json>> {
    Ok(match task {
        "smiles_read" => [false,true].into_iter().flat_map(|sanitize|[false,true].map(move|remove_hs|json!({"sanitize":sanitize,"remove_hs":remove_hs}))).collect(),
        "kekulize"=>[false,true].map(|clear|json!({"clear_aromatic_flags":clear})).into(),
        "add_hs"=>[false,true].into_iter().flat_map(|explicit|[false,true].map(move|first|json!({"explicit_only":explicit,"first_atom_only":first}))).collect(),
        "remove_hs"=>[false,true].into_iter().flat_map(|sanitize|[false,true].map(move|isotopes|json!({"sanitize":sanitize,"remove_isotopes":isotopes,"remove_mapped":true,"remove_hydrides":false}))).collect(),
        "potential_stereo"=>[false,true].into_iter().flat_map(|clean|[false,true].into_iter().flat_map(move|flag|[false,true].map(move|allow|json!({"clean":clean,"flag_possible":flag,"allow_nontetrahedral":allow})))).collect(),
        "distance_matrix"=>[false,true].into_iter().flat_map(|bo|[false,true].map(move|aw|json!({"use_bond_order":bo,"use_atom_weights":aw}))).collect(),
        task if TASKS.iter().any(|(name,_)|*name==task)=>vec![json!({})],
        _=>return Err(format!("unknown core task {task}").into()),
    })
}
pub fn normalize(task: &str, p: Json) -> Result<Json> {
    let defaults = profiles(task)?;
    let fields = p.as_object().ok_or("core parameters must be an object")?;
    let expected = defaults[0].as_object().ok_or("invalid core defaults")?;
    if fields.len() != expected.len()
        || fields
            .iter()
            .any(|(key, value)| !expected.contains_key(key) || !value.is_boolean())
    {
        return Err(format!("invalid parameters for {task}").into());
    }
    Ok(p)
}
fn flag(p: &Json, key: &str) -> Result<bool> {
    p[key]
        .as_bool()
        .ok_or_else(|| format!("missing boolean {key}").into())
}
fn add_params(p: &Json) -> Result<cosmolkit::AddHsParams> {
    Ok(cosmolkit::AddHsParams {
        explicit_only: flag(p, "explicit_only")?,
        only_on_atoms: if flag(p, "first_atom_only")? {
            Some(vec![cosmolkit::AtomId::new(0)])
        } else {
            None
        },
        ..Default::default()
    })
}
fn remove_params(p: &Json) -> Result<cosmolkit::RemoveHsParams> {
    Ok(cosmolkit::RemoveHsParams {
        sanitize: flag(p, "sanitize")?,
        remove_isotopes: flag(p, "remove_isotopes")?,
        remove_mapped: flag(p, "remove_mapped")?,
        remove_hydrides: flag(p, "remove_hydrides")?,
        ..Default::default()
    })
}
fn kekulize_params(p: &Json) -> Result<cosmolkit::KekulizeParams> {
    Ok(cosmolkit::KekulizeParams {
        mark_atoms_bonds: flag(p, "clear_aromatic_flags")?,
        ..Default::default()
    })
}

/// Stable source fields, in atom/bond order. No private cache access or derived guesses.
pub fn state(m: &Molecule) -> Value {
    Value::List(vec![
        Value::List(
            m.atoms()
                .iter()
                .map(|a| {
                    Value::List(vec![
                        Value::Int(a.atomic_number() as i64),
                        Value::Int(a.formal_charge() as i64),
                        Value::Int(a.isotope().unwrap_or(0) as i64),
                        Value::Int(a.explicit_hydrogens() as i64),
                        Value::Bool(a.no_implicit()),
                        Value::Bool(a.is_aromatic()),
                        Value::Int(a.chiral_tag().rdkit_code()),
                        Value::Int(a.hybridization().rdkit_code()),
                        Value::Int(a.radical_electrons() as i64),
                        Value::Int(a.atom_map().unwrap_or(0) as i64),
                        Value::Int(a.chiral_permutation().unwrap_or(0) as i64),
                    ])
                })
                .collect(),
        ),
        Value::List(
            m.bonds()
                .iter()
                .map(|b| {
                    Value::List(vec![
                        Value::Int(b.begin().index() as i64),
                        Value::Int(b.end().index() as i64),
                        Value::Int(b.order().rdkit_code()),
                        Value::Bool(b.is_aromatic()),
                        Value::Bool(b.is_conjugated()),
                        Value::Int(b.direction().rdkit_code()),
                        Value::Int(b.stereo().rdkit_code()),
                        Value::List(
                            b.stereo_atom_references()
                                .iter()
                                .map(|a| Value::Int(a.index() as i64))
                                .collect(),
                        ),
                    ])
                })
                .collect(),
        ),
    ])
}

pub fn scalar(m: &Molecule, task: &str, p: &Json) -> Result<Value> {
    Ok(match task {
        "smiles_read" => state(m),
        "sanitize" => state(&m.sanitize()?),
        "kekulize" => state(&m.with_kekulized_bonds_with_params(&kekulize_params(p)?)?),
        "add_hs" => state(&m.with_hydrogens_with_params(&add_params(p)?)?),
        "remove_hs" => state(&m.without_hydrogens_with_params(&remove_params(p)?)?),
        "distance_matrix" => {
            let matrix = m.distance_matrix_with_params(&cosmolkit::DistanceMatrixParams {
                use_bond_order: flag(p, "use_bond_order")?,
                use_atom_weights: flag(p, "use_atom_weights")?,
            })?;
            Value::List(vec![
                Value::Int(matrix.dimension() as i64),
                Value::List(matrix.values().iter().map(|v| Value::float(*v)).collect()),
            ])
        }
        "potential_stereo" => {
            use cosmolkit::{
                PotentialStereoCenter as C, PotentialStereoDescriptor as D,
                PotentialStereoSpecified as S, PotentialStereoType as T,
            };
            let output = m.potential_stereo_with_params(&cosmolkit::PotentialStereoParams {
                clean: flag(p, "clean")?,
                flag_possible: flag(p, "flag_possible")?,
                allow_nontetrahedral: flag(p, "allow_nontetrahedral")?,
            })?;
            let rows = output
                .stereo
                .into_iter()
                .map(|s| {
                    Value::List(vec![
                        Value::Int(match s.stereo_type {
                            T::AtomTetrahedral => 1,
                            T::AtomSquarePlanar => 2,
                            T::AtomTrigonalBipyramidal => 3,
                            T::AtomOctahedral => 4,
                            T::BondDouble => 5,
                            T::BondCumuleneEven => 6,
                            T::BondAtropisomer => 7,
                        }),
                        Value::Int(match s.specified {
                            S::Unspecified => 0,
                            S::Specified => 1,
                            S::Unknown => 2,
                        }),
                        Value::Int(match s.centered_on {
                            C::Atom(a) => a.index(),
                            C::Bond(b) => b.index(),
                        } as i64),
                        Value::Int(match s.descriptor {
                            D::None => 0,
                            D::TetrahedralClockwise => 1,
                            D::TetrahedralCounterclockwise => 2,
                            D::BondCis => 3,
                            D::BondTrans => 4,
                            D::BondAtropCw => 5,
                            D::BondAtropCcw => 6,
                        }),
                        Value::Int(s.permutation as i64),
                        Value::List(
                            s.controlling_atoms
                                .into_iter()
                                .map(|a| {
                                    Value::Int(a.map_or(u32::MAX as usize, |a| a.index()) as i64)
                                })
                                .collect(),
                        ),
                    ])
                })
                .collect();
            Value::List(vec![
                Value::List(rows),
                output.cleaned_molecule.as_ref().map_or(Value::Null, state),
            ])
        }
        "fragments" => Value::List(m.fragments()?.iter().map(state).collect()),
        "largest_fragment" => state(&m.largest_fragment()?),
        "murcko_scaffold" => Value::Bytes(m.murcko_scaffold()?.to_smiles()?.into_bytes()),
        "net_scaffold" => Value::Bytes(m.net_scaffold()?.to_smiles()?.into_bytes()),
        "murcko_decompose" => state(&m.murcko_decompose()?),
        _ => return Err(format!("unknown core task {task}").into()),
    })
}

pub fn evaluate(
    task: &str,
    profiles: &[Json],
    rows: &[String],
    first_index: usize,
    threads: usize,
    pool: &rayon::ThreadPool,
) -> Result<Vec<Vec<Vec<u8>>>> {
    let ordinary = 10usize.saturating_sub(first_index).min(rows.len());
    let has_batch = matches!(
        task,
        "smiles_read" | "sanitize" | "kekulize" | "add_hs" | "remove_hs"
    );
    let mut results = vec![vec![Vec::new(); profiles.len()]; rows.len()];
    for (pi, p) in profiles.iter().enumerate() {
        let params = cosmolkit::SmilesParseParams {
            sanitize: if task == "smiles_read" {
                flag(p, "sanitize")?
            } else {
                task != "sanitize"
            },
            remove_hs: if task == "smiles_read" {
                flag(p, "remove_hs")?
            } else {
                task != "sanitize" && task != "remove_hs"
            },
            ..Default::default()
        };
        let n = if has_batch { ordinary } else { rows.len() };
        let singles = pool.install(|| {
            rows[..n]
                .par_iter()
                .map(
                    |text| match Molecule::from_smiles_with_params(text, &params) {
                        Ok(m) => scalar(&m, task, p)
                            .map(tasks::success)
                            .unwrap_or_else(|e| tasks::failed_error(e.as_ref())),
                        Err(e) => tasks::parse_failure(&e),
                    },
                )
                .collect::<Vec<_>>()
        });
        for (i, v) in singles.into_iter().enumerate() {
            results[i][pi] = v;
        }
        if n == rows.len() {
            continue;
        }
        let runtime = BatchParams {
            errors: Some(BatchErrorMode::KeepErrors),
            n_jobs: Some(threads),
            progress_bar: Some(false),
        };
        let parsed = MoleculeBatch::from_smiles_list_with_params(&rows[n..], &params, &runtime)?;
        let transformed = match task {
            "smiles_read" => parsed.clone(),
            "sanitize" => parsed.sanitize_with_params(&Default::default(), &runtime)?,
            "kekulize" => {
                parsed.with_kekulized_bonds_with_params(&kekulize_params(p)?, &runtime)?
            }
            "add_hs" => parsed.with_hydrogens_with_params(&add_params(p)?, &runtime)?,
            "remove_hs" => parsed.without_hydrogens_with_params(&remove_params(p)?, &runtime)?,
            _ => unreachable!(),
        };
        if transformed.len() != rows.len() - n {
            return Err("batch transform row count mismatch".into());
        }
        for (i, record) in transformed.records().iter().enumerate() {
            results[n + i][pi] = match record {
                BatchRecord::Molecule(m) => tasks::success(state(m)),
                BatchRecord::Error(e) => {
                    use std::error::Error;
                    if matches!(parsed.get(i), Some(BatchRecord::Error(_))) {
                        e.source()
                            .and_then(|e| e.downcast_ref::<cosmolkit::SmilesError>())
                            .map(tasks::parse_failure)
                            .unwrap_or_else(|| tasks::failed_error(e))
                    } else {
                        tasks::failed_error(e)
                    }
                }
            };
        }
    }
    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn core_profiles_and_batch_split_are_explicit() {
        for (task, _) in TASKS {
            for p in profiles(task).unwrap() {
                assert_eq!(normalize(task, p.clone()).unwrap(), p);
            }
        }
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .unwrap();
        let rows = vec!["CCO".into(); 12];
        let results = evaluate("add_hs", &profiles("add_hs").unwrap(), &rows, 0, 2, &pool).unwrap();
        assert_eq!(results.len(), 12);
        assert!(results.iter().all(|r| r == &results[0]));
    }
}
