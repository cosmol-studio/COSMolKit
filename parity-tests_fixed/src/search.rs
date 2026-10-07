//! Registered standalone matching parity through canonical public APIs.
use crate::registry::{Input, Record, SmilesCase, Value};
use cosmolkit::{Molecule, SubstructMatchParams, search};
use serde::{Deserialize, Serialize};

/// Every scalar matcher parameter is frozen in each recipe. Callbacks are absent;
/// their Python lifetime/error paths have separate focused binding regressions.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Options {
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
    pub extra_atom_check_overrides_default_check: bool,
    pub extra_bond_check_overrides_default_check: bool,
    pub use_generic_matchers: bool,
}
impl Options {
    fn source_defaults() -> Self {
        let p = SubstructMatchParams::default();
        Self {
            max_matches: p.max_matches,
            uniquify: p.uniquify,
            use_chirality: p.use_chirality,
            use_enhanced_stereo: p.use_enhanced_stereo,
            specified_stereo_query_matches_unspecified: p
                .specified_stereo_query_matches_unspecified,
            use_query_query_matches: p.use_query_query_matches,
            recursion_possible: p.recursion_possible,
            max_recursive_matches: p.max_recursive_matches,
            num_threads: p.num_threads,
            aromatic_matches_conjugated: p.aromatic_matches_conjugated,
            aromatic_matches_single_or_double: p.aromatic_matches_single_or_double,
            atom_properties: p.atom_properties,
            bond_properties: p.bond_properties,
            extra_atom_check_overrides_default_check: p.extra_atom_check_overrides_default_check,
            extra_bond_check_overrides_default_check: p.extra_bond_check_overrides_default_check,
            use_generic_matchers: p.use_generic_matchers,
        }
    }
    fn canonical(&self) -> SubstructMatchParams {
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
            extra_atom_check_overrides_default_check: self.extra_atom_check_overrides_default_check,
            extra_bond_check_overrides_default_check: self.extra_bond_check_overrides_default_check,
            use_generic_matchers: self.use_generic_matchers,
            ..Default::default()
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Profile {
    pub id: String,
    pub smarts: String,
    pub options: Options,
    pub default_entrypoints: bool,
}
pub fn profiles() -> Vec<Profile> {
    let mut rows = Vec::new();
    let mut add = |id: &str, smarts: &str, default_entrypoints: bool, edit: fn(&mut Options)| {
        let mut options = Options::source_defaults();
        edit(&mut options);
        rows.push(Profile {
            id: id.into(),
            smarts: smarts.into(),
            options,
            default_entrypoints,
        });
    };
    add("default_CO", "CO", true, |_| {});
    add("default_aromatic_ring", "c1ccccc1", true, |_| {});
    add("oriented_non_unique", "CC", false, |p| p.uniquify = false);
    add("match_limit_three", "[C,N,O]", false, |p| p.max_matches = 3);
    add("recursive_enabled", "[$(C-O)]", false, |_| {});
    add("recursive_disabled", "[$(C-O)]", false, |p| {
        p.recursion_possible = false
    });
    add("recursive_limit_one", "[$(C-O)]", false, |p| {
        p.max_recursive_matches = 1
    });
    add("stereo_ignored", "N[C@H](C)C(=O)O", false, |_| {});
    add("stereo_required", "N[C@H](C)C(=O)O", false, |p| {
        p.use_chirality = true
    });
    add(
        "stereo_unspecified_allowed",
        "N[C@H](C)C(=O)O",
        false,
        |p| {
            p.use_chirality = true;
            p.specified_stereo_query_matches_unspecified = true;
        },
    );
    add("enhanced_stereo_required", "N[C@H](C)C(=O)O", false, |p| {
        p.use_chirality = true;
        p.use_enhanced_stereo = true;
    });
    add(
        "generic_alkyl",
        "C* |atomProp:1._QueryAtomGenericLabel.ALK|",
        false,
        |p| p.use_generic_matchers = true,
    );
    add(
        "exact_atom_property",
        "C |atomProp:0.foo.one|",
        false,
        |p| p.atom_properties = vec!["foo".into()],
    );
    add("exact_bond_property", "CC", false, |p| {
        p.bond_properties = vec!["foo".into()]
    });
    add("aromatic_conjugated", "[#6]:[#6]", false, |p| {
        p.aromatic_matches_conjugated = true
    });
    add("aromatic_single_double", "[#6]:[#6]", false, |p| {
        p.aromatic_matches_single_or_double = true
    });
    add("query_query_option", "[#6]", false, |p| {
        p.use_query_query_matches = true
    });
    add("thread_parameter_two", "CO", false, |p| p.num_threads = 2);
    rows
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SearchInput {
    pub case: SmilesCase,
    pub profile: Profile,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Outcome {
    Matches {
        atom_mappings: Vec<Vec<usize>>,
        bond_mappings: Vec<Vec<usize>>,
        first_mapping: Option<Vec<usize>>,
        has_match: bool,
        compiled_atom_mappings: Option<Vec<Vec<usize>>>,
    },
    Error {
        stage: String,
        kind: String,
    },
}
pub fn validate_output(input: &SearchInput, out: &Outcome) -> Result<(), String> {
    if !profiles().contains(&input.profile) {
        return Err("unregistered matching profile".into());
    }
    if let Outcome::Matches {
        atom_mappings,
        bond_mappings,
        first_mapping,
        has_match,
        compiled_atom_mappings,
    } = out
    {
        if atom_mappings.len() != bond_mappings.len()
            || atom_mappings.len() > input.profile.options.max_matches
        {
            return Err("matching rows or limit mismatch".into());
        }
        if *has_match != !atom_mappings.is_empty()
            || first_mapping.as_ref() != atom_mappings.first()
        {
            return Err("first/boolean/full result inconsistency".into());
        }
        if compiled_atom_mappings.is_some() != input.profile.default_entrypoints {
            return Err("compiled default coverage mismatch".into());
        }
        if compiled_atom_mappings
            .as_ref()
            .is_some_and(|v| v != atom_mappings)
        {
            return Err("compiled matching differs from default".into());
        }
    }
    Ok(())
}
pub fn run(input: &SearchInput) -> Result<Record, String> {
    let output = match Molecule::from_smiles(&input.case.smiles) {
        Err(_) => Outcome::Error {
            stage: "Parse".into(),
            kind: "SmilesParse".into(),
        },
        Ok(mol) => {
            let q = search::parse_smarts(&input.profile.smarts).map_err(|e| e.to_string())?;
            let matches = if input.profile.default_entrypoints {
                mol.substruct_matches(&q)
            } else {
                mol.substruct_matches_with_params(&q, &input.profile.options.canonical())
            }
            .map_err(|e| e.to_string())?;
            let atom_mappings = matches
                .iter()
                .map(|m| m.atom_mapping.clone())
                .collect::<Vec<_>>();
            let bond_mappings = matches.iter().map(|m| m.bond_mapping.clone()).collect();
            let (first_mapping, has_match, compiled_atom_mappings) =
                if input.profile.default_entrypoints {
                    let first = mol
                        .substruct_match(&q)
                        .map_err(|e| e.to_string())?
                        .map(|m| m.atom_mapping);
                    let has = mol.has_substruct_match(&q).map_err(|e| e.to_string())?;
                    let compiled = search::compile_query(&q).map_err(|e| e.to_string())?;
                    let compiled_matches = mol
                        .substruct_matches_compiled(&compiled)
                        .map_err(|e| e.to_string())?;
                    (
                        first,
                        has,
                        Some(
                            compiled_matches
                                .into_iter()
                                .map(|m| m.atom_mapping)
                                .collect(),
                        ),
                    )
                } else {
                    (
                        atom_mappings.first().cloned(),
                        !atom_mappings.is_empty(),
                        None,
                    )
                };
            Outcome::Matches {
                atom_mappings,
                bond_mappings,
                first_mapping,
                has_match,
                compiled_atom_mappings,
            }
        }
    };
    Ok(Record {
        input: Input::Search(input.clone()),
        output: Value::Search(output),
    })
}
