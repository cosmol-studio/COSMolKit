//! Search and tautomer profiles, with direct public calls and complete results.
use crate::{Result, tasks, value::Value};
use cosmolkit::{Molecule, SubstructMatchParams, TautomerParams};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use serde_json::{Value as Json, json};

pub const TASKS: &[(&str, &str)] = &[("substructure", "search"), ("tautomer", "tautomer")];

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SearchProfile {
    pub smarts: String,
    pub max_matches: usize,
    pub uniquify: bool,
    pub use_chirality: bool,
    pub use_enhanced_stereo: bool,
    pub specified_stereo_query_matches_unspecified: bool,
    pub use_query_query_matches: bool,
    pub recursion_possible: bool,
    pub max_recursive_matches: usize,
    pub num_threads: i32,
    pub aromatic_matches_conjugated: bool,
    pub aromatic_matches_single_or_double: bool,
    pub atom_properties: Vec<String>,
    pub bond_properties: Vec<String>,
    pub use_generic_matchers: bool,
    pub compiled_default: bool,
}
impl SearchProfile {
    fn new(smarts: &str) -> Self {
        Self {
            smarts: smarts.into(),
            max_matches: 1000,
            uniquify: true,
            use_chirality: false,
            use_enhanced_stereo: false,
            specified_stereo_query_matches_unspecified: false,
            use_query_query_matches: false,
            recursion_possible: true,
            max_recursive_matches: 1000,
            num_threads: 1,
            aromatic_matches_conjugated: false,
            aromatic_matches_single_or_double: false,
            atom_properties: vec![],
            bond_properties: vec![],
            use_generic_matchers: false,
            compiled_default: false,
        }
    }
    fn params(&self) -> SubstructMatchParams {
        SubstructMatchParams {
            max_matches: self.max_matches,
            uniquify: self.uniquify,
            use_chirality: self.use_chirality,
            use_enhanced_stereo: self.use_enhanced_stereo,
            specified_stereo_query_matches_unspecified: self
                .specified_stereo_query_matches_unspecified,
            use_query_query_matches: self.use_query_query_matches,
            recursion_possible: self.recursion_possible,
            max_recursive_matches: self.max_recursive_matches,
            num_threads: self.num_threads,
            aromatic_matches_conjugated: self.aromatic_matches_conjugated,
            aromatic_matches_single_or_double: self.aromatic_matches_single_or_double,
            atom_properties: self.atom_properties.clone(),
            bond_properties: self.bond_properties.clone(),
            use_generic_matchers: self.use_generic_matchers,
            ..Default::default()
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TautomerProfile {
    pub catalog: String,
    pub max_tautomers: u32,
    pub max_transforms: u32,
    pub remove_sp3_stereo: bool,
    pub remove_bond_stereo: bool,
    pub remove_isotopic_hydrogens: bool,
    pub reassign_stereo: bool,
}
impl TautomerProfile {
    pub(crate) fn params(&self) -> Result<TautomerParams> {
        let p = match self.catalog.as_str() {
            "current" => TautomerParams::default(),
            "v1" => TautomerParams::v1()?,
            _ => return Err("catalog must be current or v1".into()),
        };
        Ok(p.with_max_tautomers(self.max_tautomers)
            .with_max_transforms(self.max_transforms)
            .with_remove_sp3_stereo(self.remove_sp3_stereo)
            .with_remove_bond_stereo(self.remove_bond_stereo)
            .with_remove_isotopic_hydrogens(self.remove_isotopic_hydrogens)
            .with_reassign_stereo(self.reassign_stereo))
    }
}
pub fn profiles(task: &str) -> Result<Vec<Json>> {
    match task {
        "substructure" => {
            let mut rows = Vec::new();
            let mut add = |smarts: &str, edit: fn(&mut SearchProfile)| {
                let mut p = SearchProfile::new(smarts);
                edit(&mut p);
                rows.push(serde_json::to_value(p).expect("serializable search profile"));
            };
            add("CO", |p| p.compiled_default = true);
            add("c1ccccc1", |p| p.compiled_default = true);
            add("CC", |p| p.uniquify = false);
            add("[C,N,O]", |p| p.max_matches = 3);
            add("[$(C-O)]", |_| {});
            add("[$(C-O)]", |p| p.recursion_possible = false);
            add("[$(C-O)]", |p| p.max_recursive_matches = 1);
            add("N[C@H](C)C(=O)O", |_| {});
            add("N[C@H](C)C(=O)O", |p| p.use_chirality = true);
            add("N[C@H](C)C(=O)O", |p| {
                p.use_chirality = true;
                p.specified_stereo_query_matches_unspecified = true;
            });
            add("N[C@H](C)C(=O)O", |p| {
                p.use_chirality = true;
                p.use_enhanced_stereo = true;
            });
            add("C* |atomProp:1._QueryAtomGenericLabel.ALK|", |p| {
                p.use_generic_matchers = true
            });
            add("C |atomProp:0.foo.one|", |p| {
                p.atom_properties = vec!["foo".into()]
            });
            add("CC", |p| p.bond_properties = vec!["foo".into()]);
            add("[#6]:[#6]", |p| p.aromatic_matches_conjugated = true);
            add("[#6]:[#6]", |p| p.aromatic_matches_single_or_double = true);
            add("[#6]", |p| p.use_query_query_matches = true);
            add("CO", |p| p.num_threads = 2);
            Ok(rows)
        }
        "tautomer" => {
            let base = json!({"catalog":"current","max_tautomers":1000,"max_transforms":1000,"remove_sp3_stereo":true,"remove_bond_stereo":true,"remove_isotopic_hydrogens":true,"reassign_stereo":true});
            let mut rows = vec![base.clone()];
            for (field, value) in [
                ("catalog", json!("v1")),
                ("max_tautomers", json!(1)),
                ("max_transforms", json!(1)),
                ("remove_sp3_stereo", json!(false)),
                ("remove_bond_stereo", json!(false)),
                ("remove_isotopic_hydrogens", json!(false)),
                ("reassign_stereo", json!(false)),
            ] {
                let mut p = base.clone();
                p[field] = value;
                rows.push(p);
            }
            Ok(rows)
        }
        _ => Err(format!("unknown search/tautomer task {task}").into()),
    }
}
pub fn normalize(task: &str, p: Json) -> Result<Json> {
    match task {
        "substructure" => {
            let p: SearchProfile = serde_json::from_value(p)?;
            cosmolkit::search::parse_smarts(&p.smarts)?;
            if p.compiled_default {
                let mut default = SearchProfile::new(&p.smarts);
                default.compiled_default = true;
                if serde_json::to_value(&p)? != serde_json::to_value(default)? {
                    return Err("compiled query uses default matching parameters only".into());
                }
            }
            Ok(serde_json::to_value(p)?)
        }
        "tautomer" => {
            let p: TautomerProfile = serde_json::from_value(p)?;
            if !matches!(p.catalog.as_str(), "current" | "v1") {
                return Err("unknown tautomer catalog".into());
            }
            Ok(serde_json::to_value(p)?)
        }
        _ => Err(format!("unknown search/tautomer task {task}").into()),
    }
}
fn indices(rows: impl IntoIterator<Item = usize>) -> Value {
    Value::List(rows.into_iter().map(|v| Value::Int(v as i64)).collect())
}
fn score(m: &Molecule) -> Result<Value> {
    let s = m.tautomer_score()?;
    Ok(Value::List(vec![
        Value::Int(s.ring() as i64),
        Value::Int(s.substructure() as i64),
        Value::Int(s.hetero_hydrogen() as i64),
        Value::Int(s.total() as i64),
    ]))
}
// Full graph observations are fixed-width here: a 1,000-tautomer peptide
// must not allocate a nested Value for every small integer. Widths below
// match the public detached fields; no graph/CIP/stereo field is omitted.
fn tautomer_state(m: &Molecule) -> Result<Value> {
    let mut bytes = Vec::new();
    bytes.extend((m.num_atoms() as u32).to_le_bytes());
    for a in m.atoms() {
        bytes.push(a.atomic_number());
        bytes.extend(i16::from(a.formal_charge()).to_le_bytes());
        bytes.extend(a.isotope().unwrap_or(0).to_le_bytes());
        bytes.extend([
            a.explicit_hydrogens(),
            u8::from(a.no_implicit()),
            u8::from(a.is_aromatic()),
            a.chiral_tag().rdkit_code() as u8,
            a.hybridization().rdkit_code() as u8,
            a.radical_electrons(),
        ]);
        bytes.extend(a.atom_map().unwrap_or(0).to_le_bytes());
        bytes.extend(a.chiral_permutation().unwrap_or(0).to_le_bytes());
        if let Some(cip) = a.prop("_CIPCode") {
            let text = cip.as_string()?;
            let raw = text.as_bytes();
            bytes.extend(u32::try_from(raw.len())?.to_le_bytes());
            bytes.extend(raw);
        } else {
            bytes.extend(u32::MAX.to_le_bytes());
        }
    }
    bytes.extend((m.num_bonds() as u32).to_le_bytes());
    for b in m.bonds() {
        bytes.extend((b.begin().index() as u32).to_le_bytes());
        bytes.extend((b.end().index() as u32).to_le_bytes());
        bytes.extend([
            b.order().rdkit_code() as u8,
            u8::from(b.is_aromatic()),
            u8::from(b.is_conjugated()),
            b.direction().rdkit_code() as u8,
            b.stereo().rdkit_code() as u8,
        ]);
        bytes.extend((b.stereo_atom_references().len() as u32).to_le_bytes());
        for a in b.stereo_atom_references() {
            bytes.extend((a.index() as u32).to_le_bytes());
        }
    }
    Ok(Value::Bytes(bytes))
}
pub(crate) fn tautomer(m: &Molecule, p: &TautomerParams) -> Result<Value> {
    let enumeration = m.enumerate_tautomers_with_params(p)?;
    let rows = enumeration
        .entries()
        .map(|(key, m)| {
            Ok(Value::List(vec![
                Value::Bytes(key.as_bytes().to_vec()),
                Value::Bytes(m.to_smiles()?.into_bytes()),
                tautomer_state(m)?,
                score(m)?,
            ]))
        })
        .collect::<Result<Vec<_>>>()?;
    use cosmolkit::TautomerEnumerationStatus as S;
    let status = match enumeration.status() {
        S::Completed => 0,
        S::MaxTautomersReached => 1,
        S::MaxTransformsReached => 2,
        S::Canceled => 3,
    };
    let chosen = enumeration.canonical_tautomer_with_params(p)?;
    let direct = m.canonical_tautomer_with_params(p)?;
    Ok(Value::List(vec![
        Value::Int(status),
        indices(enumeration.modified_atoms().iter().map(|a| a.index())),
        indices(enumeration.modified_bonds().iter().map(|b| b.index())),
        Value::List(rows),
        Value::List(vec![
            Value::Bytes(chosen.to_smiles()?.into_bytes()),
            tautomer_state(&chosen)?,
            score(&chosen)?,
        ]),
        Value::List(vec![
            Value::Bytes(direct.to_smiles()?.into_bytes()),
            tautomer_state(&direct)?,
            score(&direct)?,
        ]),
    ]))
}
fn search(
    m: &Molecule,
    p: &SearchProfile,
    q: &cosmolkit::QueryGraph,
    compiled: Option<&cosmolkit::CompiledQuery>,
) -> Result<Value> {
    let params = p.params();
    let found = m.substruct_matches_with_params(q, &params)?;
    let first = m
        .substruct_match_with_params(q, &params)?
        .map_or(Value::Null, |v| indices(v.atom_mapping));
    let has = m.has_substruct_match_with_params(q, &params)?;
    let compiled = compiled
        .map(|q| {
            m.substruct_matches_compiled(q).map(|rows| {
                Value::List(rows.into_iter().map(|v| indices(v.atom_mapping)).collect())
            })
        })
        .transpose()?
        .unwrap_or(Value::Null);
    Ok(Value::List(vec![
        Value::List(
            found
                .iter()
                .map(|v| indices(v.atom_mapping.clone()))
                .collect(),
        ),
        Value::List(
            found
                .iter()
                .map(|v| indices(v.bond_mapping.clone()))
                .collect(),
        ),
        first,
        Value::Bool(has),
        compiled,
    ]))
}
pub struct Prepared {
    search_profiles: Vec<SearchProfile>,
    queries: Vec<cosmolkit::QueryGraph>,
    compiled: Vec<Option<cosmolkit::CompiledQuery>>,
    tautomer_params: Vec<TautomerParams>,
    profile_count: usize,
}

pub fn prepare(task: &str, profiles: &[Json]) -> Result<Prepared> {
    let search_profiles = if task == "substructure" {
        profiles
            .iter()
            .cloned()
            .map(serde_json::from_value::<SearchProfile>)
            .collect::<std::result::Result<Vec<_>, _>>()?
    } else {
        vec![]
    };
    let queries = search_profiles
        .iter()
        .map(|p| cosmolkit::search::parse_smarts(&p.smarts))
        .collect::<std::result::Result<Vec<_>, _>>()?;
    let compiled = queries
        .iter()
        .zip(&search_profiles)
        .map(|(q, p)| {
            p.compiled_default
                .then(|| cosmolkit::search::compile_query(q))
                .transpose()
        })
        .collect::<std::result::Result<Vec<_>, _>>()?;
    let tautomer_params = if task == "tautomer" {
        profiles
            .iter()
            .cloned()
            .map(|p| serde_json::from_value::<TautomerProfile>(p)?.params())
            .collect::<Result<Vec<_>>>()?
    } else {
        vec![]
    };
    Ok(Prepared {
        search_profiles,
        queries,
        compiled,
        tautomer_params,
        profile_count: profiles.len(),
    })
}

impl Prepared {
    pub fn evaluate(&self, rows: &[String], pool: &rayon::ThreadPool) -> Vec<Vec<Vec<u8>>> {
        pool.install(|| {
            rows.par_iter()
                .map(|text| match Molecule::from_smiles(text) {
                    Err(e) => vec![tasks::parse_failure(&e); self.profile_count],
                    Ok(m) => (0..self.profile_count)
                        .map(|i| {
                            let result = if !self.tautomer_params.is_empty() {
                                tautomer(&m, &self.tautomer_params[i])
                            } else {
                                search(
                                    &m,
                                    &self.search_profiles[i],
                                    &self.queries[i],
                                    self.compiled[i].as_ref(),
                                )
                            };
                            result
                                .map(tasks::success)
                                .unwrap_or_else(|e| tasks::failed_error(e.as_ref()))
                        })
                        .collect(),
                })
                .collect()
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn retained_search_and_tautomer_parameter_matrices() {
        assert_eq!(profiles("substructure").unwrap().len(), 18);
        assert_eq!(profiles("tautomer").unwrap().len(), 8);
        for (task, _) in TASKS {
            for p in profiles(task).unwrap() {
                assert_eq!(normalize(task, p.clone()).unwrap(), p);
            }
        }
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .unwrap();
        for task in ["substructure", "tautomer"] {
            let prepared = prepare(task, &profiles(task).unwrap()).unwrap();
            let result = prepared.evaluate(&["CCO".into()], &pool);
            assert_eq!(result[0].len(), profiles(task).unwrap().len());
            assert!(result[0].iter().all(|value| value[0] == 0));
        }
    }
}
