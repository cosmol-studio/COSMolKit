//! Explicit public calls and parameter profiles. No chemistry algorithms live here.
use crate::{Result, value::Value};
use cosmolkit::{
    AtomId, BatchErrorMode, BatchParams, BatchQueryParams, BatchRecord, Molecule, MoleculeBatch,
    SmilesWriteParams,
};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::collections::BTreeMap;

pub const RDKIT: &str = "2026.09.1";
pub const TASKS: &[(&str, &str)] = &[
    ("smiles_write", "smiles"),
    ("cip", "cip"),
    ("valence", "valence"),
    ("fingerprint_morgan", "fingerprints"),
    ("fingerprint_atom_pair", "fingerprints"),
    ("fingerprint_torsion", "fingerprints"),
    ("fingerprint_topological", "fingerprints"),
    ("fingerprint_layered", "fingerprints"),
    ("fingerprint_pattern", "fingerprints"),
    ("fingerprint_avalon", "fingerprints"),
    ("fingerprint_maccs", "fingerprints"),
    ("fingerprint_fuzzy", "fingerprints"),
];

pub fn available() -> Vec<(&'static str, &'static str)> {
    TASKS
        .iter()
        .chain(crate::descriptors::TASKS)
        .chain(crate::chemistry::TASKS)
        .chain(crate::search_tautomer::TASKS)
        .chain(crate::geometry::TASKS)
        .chain(crate::bio::TASKS)
        .chain(crate::special::TASKS)
        .copied()
        .collect()
}

pub fn fixed_input(task: &str) -> Option<&'static str> {
    match task {
        "bio_pdb_output" => Some("testdata/bio/small.jsonl"),
        _ => crate::special::input(task),
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WriterProfile {
    pub isomeric_smiles: bool,
    pub kekule: bool,
    pub canonical: bool,
    pub clean_stereo: bool,
    pub all_bonds_explicit: bool,
    pub all_hydrogens_explicit: bool,
    pub include_dative_bonds: bool,
    pub ignore_atom_map_numbers: bool,
    pub root: String,
}

impl WriterProfile {
    fn params(&self, atom_count: usize) -> Result<SmilesWriteParams> {
        Ok(SmilesWriteParams {
            isomeric_smiles: self.isomeric_smiles,
            kekule: self.kekule,
            canonical: self.canonical,
            clean_stereo: self.clean_stereo,
            all_bonds_explicit: self.all_bonds_explicit,
            all_hydrogens_explicit: self.all_hydrogens_explicit,
            include_dative_bonds: self.include_dative_bonds,
            ignore_atom_map_numbers: self.ignore_atom_map_numbers,
            rooted_at_atom: match self.root.as_str() {
                "default" => None,
                "first" => Some(AtomId::new(0)),
                "last" => Some(AtomId::new(atom_count.saturating_sub(1))),
                _ => return Err("root must be default, first or last".into()),
            },
        })
    }
}

pub fn profiles(task: &str) -> Result<Vec<serde_json::Value>> {
    if crate::special::input(task).is_some() {
        return Ok(crate::special::profiles(task));
    }
    if task == "bio_pdb_output" {
        return Ok(crate::bio::profiles());
    }
    if crate::geometry::TASKS.iter().any(|(name, _)| *name == task) {
        return crate::geometry::profiles(task);
    }
    Ok(match task {
        "smiles_write" => {
            let mut profiles = Vec::new();
            for mask in 0..256 {
                for root in ["default", "first", "last"] {
                    profiles.push(serde_json::to_value(WriterProfile {
                        isomeric_smiles: mask & 1 != 0,
                        kekule: mask & 2 != 0,
                        canonical: mask & 4 != 0,
                        clean_stereo: mask & 8 != 0,
                        all_bonds_explicit: mask & 16 != 0,
                        all_hydrogens_explicit: mask & 32 != 0,
                        include_dative_bonds: mask & 64 != 0,
                        ignore_atom_map_numbers: mask & 128 != 0,
                        root: root.into(),
                    })?);
                }
            }
            profiles
        }
        "cip" => ["all", "first_tagged_atom", "first_stereo_bond", "empty"]
            .map(|selection| json!({"selection":selection,"max_recursive_iterations":if selection=="empty" {1} else {0}})).into(),
        "valence" => vec![json!({"recalculate": true}), json!({"recalculate": false})],
        "fingerprint_morgan" => crate::fingerprints::morgan_profiles()?,
        "fingerprint_atom_pair" | "fingerprint_torsion" => {
            crate::pair_fingerprints::profiles(task)?
        }
        "fingerprint_fuzzy" => crate::fuzzy_fingerprints::profiles(),
        task if task.starts_with("fingerprint_") => crate::path_fingerprints::profiles(task)?,
        task if crate::chemistry::TASKS
            .iter()
            .any(|(name, _)| *name == task) =>
        {
            return crate::chemistry::profiles(task);
        }
        "substructure" | "tautomer" => return crate::search_tautomer::profiles(task),
        task => return crate::descriptors::profiles(task),
    })
}

pub fn normalize_profile(task: &str, profile: serde_json::Value) -> Result<serde_json::Value> {
    if crate::special::input(task).is_some() {
        return crate::special::normalize(task, profile);
    }
    if task == "bio_pdb_output" {
        return crate::bio::normalize(profile);
    }
    if crate::geometry::TASKS.iter().any(|(name, _)| *name == task) {
        return crate::geometry::normalize(task, profile);
    }
    match task {
        "fingerprint_atom_pair" | "fingerprint_torsion" => {
            crate::pair_fingerprints::normalize(task, profile)
        }
        "fingerprint_morgan" => {
            crate::fingerprints::validate_morgan(&profile)?;
            Ok(serde_json::to_value(serde_json::from_value::<
                crate::fingerprints::MorganProfile,
            >(profile)?)?)
        }
        "fingerprint_fuzzy" => crate::fuzzy_fingerprints::normalize(profile),
        task if task.starts_with("fingerprint_") => {
            crate::path_fingerprints::normalize(task, profile)
        }
        "smiles_write" => {
            let parsed: WriterProfile = serde_json::from_value(profile.clone())?;
            parsed.params(1)?;
            Ok(profile)
        }
        "cip"
            if profile.as_object().is_some_and(|p| p.len() == 2)
                && matches!(
                    profile["selection"].as_str(),
                    Some("all" | "first_tagged_atom" | "first_stereo_bond" | "empty")
                )
                && profile["max_recursive_iterations"]
                    .as_u64()
                    .is_some_and(|v| v <= u32::MAX as u64) =>
        {
            Ok(profile)
        }
        "valence"
            if profile.as_object().is_some_and(|value| {
                value.len() == 1
                    && value
                        .get("recalculate")
                        .is_some_and(serde_json::Value::is_boolean)
            }) =>
        {
            Ok(profile)
        }
        task if crate::chemistry::TASKS
            .iter()
            .any(|(name, _)| *name == task) =>
        {
            crate::chemistry::normalize(task, profile)
        }
        "substructure" | "tautomer" => crate::search_tautomer::normalize(task, profile),
        _ => crate::descriptors::normalize(task, profile),
    }
}

/// Status 0 is a calculated value; 1 is an explicit source-defined rejection.
/// 2/3/4 are upstream undefined/timeout/crash; 5 is an unclassified failure.
pub fn outcome(status: u8, value: Value) -> Vec<u8> {
    let mut bytes = vec![status];
    bytes.extend(value.encode());
    bytes
}
pub(crate) fn success(value: Value) -> Vec<u8> {
    outcome(0, value)
}
fn failure(error: impl std::fmt::Display) -> Vec<u8> {
    outcome(5, Value::Text(error.to_string()))
}
pub(crate) fn failed_error(error: &(dyn std::error::Error + 'static)) -> Vec<u8> {
    let mut cause = Some(error);
    let mut optimizer_domain = false;
    while let Some(inner) = cause {
        optimizer_domain |= inner.is::<cosmolkit::ConformerRunError>()
            || inner.is::<cosmolkit::UffOptimizationError>();
        // The kernel enum is deliberately private. Require the real public
        // optimizer domain AND its exact terminal source diagnostic, never a
        // substring of an arbitrary top-level exception message.
        if optimizer_domain
            && inner.source().is_none()
            && inner.to_string() == "OptimizerBadDirection"
        {
            return outcome(
                1,
                Value::List(vec![
                    Value::Text("BfgsBadDirection".into()),
                    Value::Text(error.to_string()),
                ]),
            );
        }
        if let Some(cosmolkit::CipLabelerError::TooManyNodes { limit }) =
            inner.downcast_ref::<cosmolkit::CipLabelerError>()
        {
            return outcome(
                1,
                Value::List(vec![
                    Value::Text("CipTooManyNodes".into()),
                    Value::UInt(*limit as u64),
                    Value::Text(error.to_string()),
                ]),
            );
        }
        if matches!(
            inner.downcast_ref::<cosmolkit::CipLabelerError>(),
            Some(cosmolkit::CipLabelerError::MaxIterationsExceeded)
        ) {
            return outcome(
                1,
                Value::List(vec![
                    Value::Text("CipMaxIterationsExceeded".into()),
                    Value::Text(error.to_string()),
                ]),
            );
        }
        if let Some(cosmolkit::ValenceError::InvalidValence { atom, message, .. }) =
            inner.downcast_ref::<cosmolkit::ValenceError>()
        {
            return outcome(
                1,
                Value::List(vec![
                    Value::Text("AtomValenceException".into()),
                    Value::Int(atom.index() as i64),
                    Value::Text(message.clone()),
                    Value::Text(error.to_string()),
                ]),
            );
        }
        if matches!(inner.downcast_ref::<cosmolkit::MolWriteError>(),Some(cosmolkit::MolWriteError::Value(message)) if message=="Atomic number not found")
        {
            return outcome(
                1,
                Value::List(vec![
                    Value::Text("PeriodicTableAtomicNumber".into()),
                    Value::Text(error.to_string()),
                ]),
            );
        }
        if let Some(cosmolkit::FingerprintError::SparseIndexOutOfRange { index, size }) =
            inner.downcast_ref::<cosmolkit::FingerprintError>()
        {
            return outcome(
                1,
                Value::List(vec![
                    Value::Text("IndexError".into()),
                    Value::UInt(*index),
                    Value::UInt(*size),
                    Value::Text(error.to_string()),
                ]),
            );
        }
        if let Some(cosmolkit::KekulizeError::NotKekulizable { problem_atoms }) =
            inner.downcast_ref::<cosmolkit::KekulizeError>()
        {
            return outcome(
                1,
                Value::List(vec![
                    Value::Text("KekulizeException".into()),
                    Value::List(
                        problem_atoms
                            .iter()
                            .map(|atom| Value::Int(atom.index() as i64))
                            .collect(),
                    ),
                    Value::Text(error.to_string()),
                ]),
            );
        }
        cause = inner.source();
    }
    failure(error)
}

/// Numeric outputs remain bit-exact. Exception spelling is not a numeric value.
/// Only this explicit upstream exception / typed CK cause pair is comparable.
/// Retain the native reference verbatim; no regeneration for message formatting.
pub fn matching_error(task: &str, expected: &[u8], actual: &[u8]) -> bool {
    if expected.first() != Some(&5) || actual.first() != Some(&1) {
        return false;
    }
    let Ok(Value::Text(native)) = crate::value::decode(&expected[1..]) else {
        return false;
    };
    let Ok(Value::List(ck)) = crate::value::decode(&actual[1..]) else {
        return false;
    };
    match task {
        "embedding" | "uff_optimize" | "uff_optimize_conformers" => {
            // BFGSOpt.h::minimize, after linearSearch:
            // RDKit✔️✔️:     CHECK_INVARIANT(status >= 0, "bad direction in linearSearch");
            // Only exception representation is paired; no numerical value,
            // timeout, unrelated invariant or untyped failure is normalized.
            native
                .starts_with("RuntimeError: Invariant Violation\n\tbad direction in linearSearch\n")
                && native
                    .lines()
                    .any(|line| line.ends_with("in file Code/Numerics/Optimizer/BFGSOpt.h"))
                && native
                    .lines()
                    .any(|line| line.trim() == "Failed Expression: status >= 0")
                && matches!(ck.as_slice(), [Value::Text(kind), Value::Text(_)] if kind == "BfgsBadDirection")
        }
        "cip" => {
            // CIPLabeler/Digraph.cpp::addNode, before insertion:
            // RDKit✔️✔️: if constexpr (MAX_NODE_COUNT > 0) {
            // RDKit✔️✔️:   if (d_nodes.size() >= MAX_NODE_COUNT) {
            // RDKit✔️✔️:     std::stringstream errmsg;
            // RDKit✔️✔️:     errmsg << "Digraph generation failed: more than " << MAX_NODE_COUNT
            // RDKit✔️✔️:            << " nodes found.";
            // RDKit✔️✔️:     throw TooManyNodesException(errmsg.str());
            // RDKit✔️✔️:   }
            // RDKit✔️✔️:   }
            // Match only the typed cause and identical limit, not arbitrary errors.
            // One bounded diagnostic comparison; calculated values remain untouched.
            (native == "RuntimeError: Max Iterations Exceeded in CIP label calculation"
                && matches!(ck.as_slice(),[Value::Text(kind),Value::Text(_)] if kind=="CipMaxIterationsExceeded"))
                || matches!(ck.as_slice(), [Value::Text(kind), Value::UInt(limit), Value::Text(_)]
                    if kind == "CipTooManyNodes"
                        && native == format!("RuntimeError: Digraph generation failed: more than {limit} nodes found."))
        }
        "smiles_write" | "sanitize" | "kekulize" => {
            if let [Value::Text(class), Value::List(atoms), Value::Text(_)] = ck.as_slice() {
                return class == "KekulizeException"
                    && native.starts_with("KekulizeException: ")
                    && native
                        .split_once("Unkekulized atoms:")
                        .and_then(|(_, atoms)| {
                            atoms
                                .split_whitespace()
                                .map(|v| v.parse::<i64>().map(Value::Int))
                                .collect::<std::result::Result<Vec<_>, _>>()
                                .ok()
                        })
                        .as_ref()
                        == Some(atoms);
            }
            matches!(ck.as_slice(),[Value::Text(class),Value::Int(atom),Value::Text(message),Value::Text(_)]
                if task=="sanitize" && class=="AtomValenceException"
                    && native.strip_prefix("AtomValenceException: ")==Some(message.as_str())
                    && message.starts_with(&format!("Explicit valence for atom # {atom} ")))
        }
        "fingerprint_avalon" => {
            // PeriodicTable.h (e.g. getRvdw) uses this exact source precondition:
            // RDKit✔️✔️: PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
            // Compare the known failure kind; retain both original diagnostics.
            native.starts_with("RuntimeError: Pre-condition Violation\n\tAtomic number not found\n")
                && native.contains("Failed Expression: atomicNumber < byanum.size()")
                && matches!(ck.first(),Some(Value::Text(kind)) if kind=="PeriodicTableAtomicNumber")
        }
        "fingerprint_atom_pair" | "fingerprint_morgan" | "fingerprint_torsion" => {
            // SparseIntVect.h::getVal/setVal preserves its signed Python payload:
            // RDKit✔️✔️: throw IndexErrorException(static_cast<int>(idx));
            // O(1) typed error projection; no fingerprint value is normalized.
            matches!(ck.as_slice(), [Value::Text(class),Value::UInt(index),Value::UInt(_),Value::Text(_)]
                if class == "IndexError" && native.strip_prefix("IndexError: ").and_then(|text|text.parse::<i32>().ok()) == Some(*index as i32))
        }
        _ => false,
    }
}
fn rejected() -> Vec<u8> {
    outcome(1, Value::Text("parse_rejected".into()))
}

pub(crate) fn parse_failure(error: &cosmolkit::SmilesError) -> Vec<u8> {
    match error {
        cosmolkit::SmilesError::Parse(_)
        | cosmolkit::SmilesError::Sanitize(_)
        | cosmolkit::SmilesError::Hydrogen(_)
        | cosmolkit::SmilesError::Stereo(_) => rejected(),
        cosmolkit::SmilesError::Construction(_) => failure(error),
    }
}

fn scalar(mol: &Molecule, task: &str, profile: &serde_json::Value) -> Result<Value> {
    Ok(match task {
        "fingerprint_morgan" => crate::fingerprints::morgan(mol, profile)?,
        "fingerprint_atom_pair" | "fingerprint_torsion" => {
            crate::pair_fingerprints::scalar(mol, task, profile)?
        }
        "fingerprint_fuzzy" => crate::fuzzy_fingerprints::scalar(mol, profile)?,
        task if task.starts_with("fingerprint_") => {
            crate::path_fingerprints::scalar(mol, task, profile)?
        }
        "smiles_write" => {
            let profile: WriterProfile = serde_json::from_value(profile.clone())?;
            Value::Bytes(
                mol.to_smiles_with_params(&profile.params(mol.num_atoms())?)?
                    .into_bytes(),
            )
        }
        "valence" => Value::List(
            mol.atom_metadata(
                profile["recalculate"]
                    .as_bool()
                    .ok_or("missing recalculate")?,
            )?
            .iter()
            .map(|row| {
                Value::List(vec![
                    Value::Int(row.degree as i64),
                    Value::Int(row.explicit_valence as i64),
                    Value::Int(row.implicit_hydrogens as i64),
                    Value::Int(row.total_hydrogens as i64),
                    Value::Int(row.total_valence as i64),
                ])
            })
            .collect(),
        ),
        "cip" => {
            use cosmolkit::{BondStereo, ChiralTag, CipLabelOptions};
            let mut options = CipLabelOptions::default().with_max_recursive_iterations(
                profile["max_recursive_iterations"]
                    .as_u64()
                    .ok_or("missing CIP recursion limit")? as u32,
            );
            match profile["selection"].as_str() {
                Some("all") => (),
                Some("first_tagged_atom") => {
                    options = options
                        .with_atoms(
                            mol.atoms()
                                .iter()
                                .find(|a| a.chiral_tag() != ChiralTag::Unspecified)
                                .map(|a| a.id()),
                        )
                        .with_bonds([]);
                }
                Some("first_stereo_bond") => {
                    options = options.with_atoms([]).with_bonds(
                        mol.bonds()
                            .iter()
                            .find(|b| !matches!(b.stereo(), BondStereo::None | BondStereo::Any))
                            .map(|b| b.id()),
                    );
                }
                Some("empty") => options = options.with_atoms([]).with_bonds([]),
                _ => return Err("unknown CIP selection".into()),
            }
            let labeled = mol.with_cip_labels_with_options(&options)?;
            let atoms = labeled
                .atoms()
                .iter()
                .map(|atom| {
                    atom.cip_descriptor().map(|label| {
                        label.map_or(Value::Null, |label| Value::Text(label.as_str().into()))
                    })
                })
                .collect::<std::result::Result<Vec<_>, _>>()?;
            let bonds = labeled
                .bonds()
                .iter()
                .map(|bond| {
                    bond.cip_descriptor().map(|label| {
                        label.map_or(Value::Null, |label| Value::Text(label.as_str().into()))
                    })
                })
                .collect::<std::result::Result<Vec<_>, _>>()?;
            Value::List(vec![Value::List(atoms), Value::List(bonds)])
        }
        _ => crate::descriptors::scalar(mol, task, profile)?,
    })
}

/// Results remain case-major. A public batch receives one homogeneous profile.
pub fn evaluate(
    task: &str,
    profiles: &[serde_json::Value],
    rows: &[String],
    first_index: usize,
    threads: usize,
    pool: &rayon::ThreadPool,
) -> Result<Vec<Vec<Vec<u8>>>> {
    if crate::special::input(task).is_some() {
        return crate::special::evaluate(task, rows, pool);
    }
    if task == "bio_pdb_output" {
        return crate::bio::evaluate(profiles, rows, pool);
    }
    if crate::geometry::TASKS.iter().any(|(name, _)| *name == task) {
        return crate::geometry::evaluate(task, profiles, rows, first_index, threads, pool);
    }
    if crate::chemistry::TASKS
        .iter()
        .any(|(name, _)| *name == task)
    {
        return crate::chemistry::evaluate(task, profiles, rows, first_index, threads, pool);
    }
    if task != "smiles_write" && !task.starts_with("fingerprint_") {
        return Ok(pool.install(|| {
            rows.par_iter()
                .map(|text| {
                    let removed = Molecule::from_smiles(text);
                    let retained = profiles.iter().any(|p| p["remove_hs"] == false).then(|| {
                        Molecule::from_smiles_with_params(
                            text,
                            &cosmolkit::SmilesParseParams {
                                remove_hs: false,
                                ..Default::default()
                            },
                        )
                    });
                    profiles
                        .iter()
                        .map(|profile| {
                            let molecule = if profile["remove_hs"] == false {
                                retained.as_ref().expect("retained profile parsed")
                            } else {
                                &removed
                            };
                            match molecule {
                                Ok(mol) => scalar(mol, task, profile)
                                    .map(success)
                                    .unwrap_or_else(|e| failed_error(e.as_ref())),
                                Err(e) => parse_failure(e),
                            }
                        })
                        .collect()
                })
                .collect::<Vec<_>>()
        }));
    }
    let ordinary = 10usize.saturating_sub(first_index).min(rows.len());
    let mut molecules: Vec<std::result::Result<Molecule, Vec<u8>>> = rows[..ordinary]
        .iter()
        .map(|text| Molecule::from_smiles(text).map_err(|error| parse_failure(&error)))
        .collect();
    if ordinary < rows.len() {
        let batch = MoleculeBatch::from_smiles_list_with_params(
            &rows[ordinary..],
            &Default::default(),
            &BatchParams {
                errors: Some(BatchErrorMode::KeepErrors),
                n_jobs: Some(threads),
                progress_bar: Some(false),
            },
        )?;
        for record in batch.records() {
            molecules.push(match record {
                BatchRecord::Molecule(mol) => Ok(mol.clone()),
                BatchRecord::Error(error) => {
                    use std::error::Error;
                    Err(error
                        .source()
                        .and_then(|error| error.downcast_ref::<cosmolkit::SmilesError>())
                        .map(parse_failure)
                        .unwrap_or_else(|| failure(error)))
                }
            });
        }
    }
    let mut results: Vec<Vec<Vec<u8>>> = molecules
        .iter()
        .map(|mol| match mol {
            Ok(_) => vec![Vec::new(); profiles.len()],
            Err(error) => vec![error.clone(); profiles.len()],
        })
        .collect();
    for (profile_index, profile) in profiles.iter().enumerate() {
        let writer = if task == "smiles_write" {
            Some(serde_json::from_value::<WriterProfile>(profile.clone())?)
        } else {
            None
        };
        let has_batch = match task {
            "smiles_write" => true,
            "fingerprint_morgan" => crate::fingerprints::morgan_has_batch(profile),
            "fingerprint_atom_pair" | "fingerprint_torsion" => {
                crate::pair_fingerprints::has_batch(task, profile)
            }
            _ => crate::path_fingerprints::has_batch(task),
        };
        if !has_batch {
            let values = pool.install(|| {
                molecules
                    .par_iter()
                    .map(|mol| {
                        mol.as_ref().ok().map(|mol| {
                            scalar(mol, task, profile)
                                .map(success)
                                .unwrap_or_else(|error| failed_error(error.as_ref()))
                        })
                    })
                    .collect::<Vec<_>>()
            });
            for (index, value) in values.into_iter().enumerate() {
                if let Some(value) = value {
                    results[index][profile_index] = value;
                }
            }
            continue;
        }
        let mut groups = BTreeMap::<(usize, usize), Vec<usize>>::new();
        for (index, mol) in molecules.iter().enumerate() {
            let Ok(mol) = mol else { continue };
            if index < ordinary {
                results[index][profile_index] = scalar(mol, task, profile)
                    .map(success)
                    .unwrap_or_else(|error| failed_error(error.as_ref()));
            } else {
                // "last" resolves to a different absolute atom index per size.
                groups
                    .entry(
                        if task == "fingerprint_layered"
                            || (task.starts_with("fingerprint_")
                                && (profile["custom_atoms"] == true
                                    || profile["custom_bonds"] == true))
                        {
                            (mol.num_atoms(), mol.num_bonds())
                        } else if writer.as_ref().is_some_and(|writer| writer.root == "last") {
                            (mol.num_atoms(), 0)
                        } else {
                            (0, 0)
                        },
                    )
                    .or_default()
                    .push(index);
            }
        }
        let grouped_results = pool.install(|| {
            groups
                .into_par_iter()
                .map(|(_, mut indices)| -> Result<Vec<(usize, Vec<u8>)>> {
                    let mut output = Vec::with_capacity(indices.len());
                    while !indices.is_empty() {
                        let atom_count = molecules[indices[0]].as_ref().unwrap().num_atoms();
                        let bond_count = molecules[indices[0]].as_ref().unwrap().num_bonds();
                        let batch = MoleculeBatch::from_records(
                            indices
                                .iter()
                                .map(|index| {
                                    BatchRecord::Molecule(
                                        molecules[*index].as_ref().unwrap().clone(),
                                    )
                                })
                                .collect(),
                            BatchErrorMode::KeepErrors,
                        )?;
                        let runtime = BatchQueryParams {
                            n_jobs: Some(threads),
                            progress_bar: Some(false),
                            ..Default::default()
                        };
                        let values = if let Some(writer) = &writer {
                            batch
                                .to_smiles_list_with_params(&writer.params(atom_count)?, &runtime)
                                .map(|values| {
                                    values
                                        .into_iter()
                                        .map(|value| {
                                            value.map(|value| Value::Bytes(value.into_bytes()))
                                        })
                                        .collect()
                                })
                        } else if task == "fingerprint_morgan" {
                            crate::fingerprints::morgan_batch(
                                &batch, profile, atom_count, bond_count, &runtime,
                            )
                        } else if matches!(task, "fingerprint_atom_pair" | "fingerprint_torsion") {
                            crate::pair_fingerprints::batch(
                                &batch, task, profile, atom_count, &runtime,
                            )
                        } else {
                            crate::path_fingerprints::batch(
                                &batch, task, profile, atom_count, &runtime,
                            )
                        };
                        match values {
                            Ok(values) => {
                                if values.len() != indices.len() {
                                    return Err(
                                        "batch result count differs from input count".into()
                                    );
                                }
                                for (index, value) in indices.iter().zip(values) {
                                    output.push((
                                        *index,
                                        success(
                                            value.ok_or("valid batch row unexpectedly missing")?,
                                        ),
                                    ));
                                }
                                break;
                            }
                            Err(error) => {
                                // The query API returns no successful payloads on error.
                                // Retain each failure, then rerun only its valid complement via batch.
                                let mut failed = Vec::new();
                                for error in error.record_errors {
                                    let index = *indices
                                        .get(error.index)
                                        .ok_or("batch error index outside input")?;
                                    output.push((index, failed_error(&error)));
                                    failed.push(index);
                                }
                                if failed.is_empty() {
                                    return Err("batch failed without indexed errors".into());
                                }
                                indices.retain(|index| !failed.contains(index));
                            }
                        }
                    }
                    Ok(output)
                })
                .collect::<Result<Vec<_>>>()
        })?;
        for group in grouped_results {
            for (index, value) in group {
                results[index][profile_index] = value;
            }
        }
    }
    Ok(results)
}

pub fn pack_row(outcomes: &[Vec<u8>]) -> Vec<u8> {
    let mut row = (outcomes.len() as u32).to_le_bytes().to_vec();
    for outcome in outcomes {
        row.extend((outcome.len() as u32).to_le_bytes());
        row.extend(outcome);
    }
    row
}

pub fn unpack_row(mut row: &[u8], profiles: usize) -> Result<Vec<&[u8]>> {
    use std::io::Read;
    let mut number = [0; 4];
    row.read_exact(&mut number)?;
    if u32::from_le_bytes(number) as usize != profiles {
        return Err("reference profile count mismatch".into());
    }
    let mut values = Vec::with_capacity(profiles);
    for _ in 0..profiles {
        row.read_exact(&mut number)?;
        let length = u32::from_le_bytes(number) as usize;
        if length == 0 || length > row.len() {
            return Err("invalid outcome length".into());
        }
        values.push(&row[..length]);
        row = &row[length..];
    }
    if !row.is_empty() {
        return Err("trailing outcome bytes".into());
    }
    Ok(values)
}

#[cfg(test)]
mod tests {
    #[test]
    fn bfgs_direction_errors_require_the_real_domain_and_exact_native_invariant() {
        use super::*;
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap();
        let profiles = crate::geometry::profiles("uff_optimize_conformers").unwrap();
        let actual = crate::geometry::evaluate(
            "uff_optimize_conformers",
            &profiles,
            &["II".into()],
            0,
            1,
            &pool,
        )
        .unwrap()
        .remove(0)
        .remove(0);
        let message = "RuntimeError: Invariant Violation\n\tbad direction in linearSearch\n\tViolation occurred on line 234 in file Code/Numerics/Optimizer/BFGSOpt.h\n\tFailed Expression: status >= 0\n\tRDKIT: 2026.09.1\n";
        let native = outcome(5, Value::Text(message.into()));
        assert_eq!(actual[0], 1);
        assert!(matching_error("uff_optimize_conformers", &native, &actual));
        assert!(!matching_error("valence", &native, &actual));
        assert!(!matching_error(
            "uff_optimize_conformers",
            &native,
            &failure("OptimizerBadDirection")
        ));
        assert!(!matching_error(
            "uff_optimize_conformers",
            &outcome(5, Value::Text(message.replace("status >= 0", "status > 0"))),
            &actual
        ));
        assert!(!matching_error(
            "uff_optimize_conformers",
            &outcome(5, Value::Text(message.replace("BFGSOpt.h", "Other.h"))),
            &actual
        ));
        let params = vec![
            serde_json::json!({"method":"etdg","count":2,"seed":20260306,
            "add_hs":true,"random_coords":false,"max_iterations":100,"prune_rms":-1.0,"timeout":10}),
        ];
        let actual = crate::geometry::evaluate(
            "embedding",
            &params,
            &["Cc1ccc([C@]2(C)O[C@H]2C)oc1=O".into()],
            0,
            1,
            &pool,
        )
        .unwrap()
        .remove(0)
        .remove(0);
        assert!(matching_error("embedding", &native, &actual));
    }

    #[test]
    fn cip_node_limit_matches_only_the_same_typed_rejection() {
        use super::*;
        let native = outcome(
            5,
            Value::Text(
                "RuntimeError: Digraph generation failed: more than 100000 nodes found.".into(),
            ),
        );
        let actual = failed_error(&cosmolkit::CipLabelerError::TooManyNodes { limit: 100_000 });
        assert!(matching_error("cip", &native, &actual));
        assert!(!matching_error("valence", &native, &actual));
        assert!(!matching_error(
            "cip",
            &native,
            &failed_error(&cosmolkit::CipLabelerError::TooManyNodes { limit: 99_999 })
        ));
        assert!(!matching_error(
            "cip",
            &native,
            &failed_error(&cosmolkit::CipLabelerError::MaxIterationsExceeded)
        ));
        assert!(!matching_error(
            "cip",
            &native,
            &failure("Digraph generation failed: more than 100000 nodes found.")
        ));
    }

    #[test]
    fn cip_selection_and_recursion_are_not_lost() {
        let profiles = super::profiles("cip").unwrap();
        assert_eq!(profiles.len(), 4);
        for p in profiles {
            super::normalize_profile("cip", p).unwrap();
        }
        assert!(
            super::normalize_profile(
                "cip",
                serde_json::json!({"selection":"empty","max_recursive_iterations":-1})
            )
            .is_err()
        );
        assert_eq!(
            crate::chemistry::profiles("potential_stereo")
                .unwrap()
                .len(),
            8
        );
        let expected = super::outcome(
            5,
            crate::value::Value::Text(
                "RuntimeError: Max Iterations Exceeded in CIP label calculation".into(),
            ),
        );
        let actual = super::failed_error(&cosmolkit::CipLabelerError::MaxIterationsExceeded);
        assert!(super::matching_error("cip", &expected, &actual));
        assert!(!super::matching_error("valence", &expected, &actual));
    }
    use super::*;

    #[test]
    fn known_kekulization_errors_require_the_same_problem_atoms() {
        let native = outcome(
            5,
            Value::Text("KekulizeException: Can't kekulize mol. Unkekulized atoms: 8 9 11".into()),
        );
        let error = cosmolkit::KekulizeError::NotKekulizable {
            problem_atoms: vec![AtomId::new(8), AtomId::new(9), AtomId::new(11)],
        };
        assert!(matching_error(
            "smiles_write",
            &native,
            &failed_error(&error)
        ));
        assert!(!matching_error("cip", &native, &failed_error(&error)));
        assert!(!matching_error(
            "smiles_write",
            &native,
            &failed_error(&cosmolkit::KekulizeError::NotKekulizable {
                problem_atoms: vec![AtomId::new(8)]
            })
        ));
        assert!(!matching_error(
            "smiles_write",
            &native,
            &failure("KekulizeException: unrelated error")
        ));
        let unknown = outcome(5, Value::Text("RuntimeError: unrelated error".into()));
        assert!(!matching_error(
            "smiles_write",
            &unknown,
            &failed_error(&error)
        ));
    }

    #[test]
    fn profile_row_counts_cannot_silently_truncate() {
        let outputs = vec![success(Value::Int(42)), rejected()];
        let row = pack_row(&outputs);
        assert_eq!(unpack_row(&row, 2).unwrap(), outputs);
        assert!(unpack_row(&row, 1).is_err());
        assert!(unpack_row(&row[..row.len() - 1], 2).is_err());
    }

    #[test]
    fn sparse_index_errors_require_the_same_native_index() {
        let error = cosmolkit::FingerprintError::SparseIndexOutOfRange {
            index: u64::from(u32::MAX),
            size: 8_388_608,
        };
        let native = outcome(5, Value::Text("IndexError: -1".into()));
        assert!(matching_error(
            "fingerprint_atom_pair",
            &native,
            &failed_error(&error)
        ));
        assert!(!matching_error(
            "fingerprint_atom_pair",
            &outcome(5, Value::Text("IndexError: -2".into())),
            &failed_error(&error)
        ));
        assert!(!matching_error("valence", &native, &failed_error(&error)));
    }

    #[test]
    fn avalon_preconditions_match_only_the_specific_periodic_table_failure() {
        let error = cosmolkit::MolWriteError::Value("Atomic number not found".into());
        let native=outcome(5,Value::Text("RuntimeError: Pre-condition Violation\n\tAtomic number not found\n\tFailed Expression: atomicNumber < byanum.size()".into()));
        assert!(matching_error(
            "fingerprint_avalon",
            &native,
            &failed_error(&error)
        ));
        let other = cosmolkit::MolWriteError::Value("other MOL failure".into());
        assert!(!matching_error(
            "fingerprint_avalon",
            &native,
            &failed_error(&other)
        ));
        assert!(!matching_error(
            "smiles_write",
            &native,
            &failed_error(&error)
        ));
    }
}
