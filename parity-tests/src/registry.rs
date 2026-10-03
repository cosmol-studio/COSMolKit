//! Executable tasks, typed input families and parameter matrices.
//! Future catalog rows are not executable registration or parity claims.
pub mod fingerprint;
pub mod fingerprint_corpus;
pub mod molecule_plan;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Operation {
    FuzzyAnd,
    FuzzyOr,
    Molecular(molecule_plan::TaskId),
}

impl Operation {
    pub fn name(self) -> &'static str {
        match self {
            Self::FuzzyAnd => "fuzzy_and",
            Self::FuzzyOr => "fuzzy_or",
            Self::Molecular(id) => id.name(),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Width {
    U32,
    U64,
}

/// Input format is part of test identity, never inferred from a filename.
/// Only registered formats have executable loaders and tests.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CorpusType {
    Smiles,
    FingerprintPairs,
    Pdb,
    Cif,
    Mmcif,
    Sdf,
}

impl CorpusType {
    pub fn name(self) -> &'static str {
        match self {
            Self::Smiles => "smiles",
            Self::FingerprintPairs => "fingerprint_pairs",
            Self::Pdb => "pdb",
            Self::Cif => "cif",
            Self::Mmcif => "mmcif",
            Self::Sdf => "sdf",
        }
    }
}

pub struct Task {
    pub operation: Operation,
    pub corpus_type: CorpusType,
    pub generator: &'static str,
}

pub const TASKS: &[Task] = &[
    Task {
        operation: Operation::FuzzyAnd,
        corpus_type: CorpusType::FingerprintPairs,
        generator: "generate_fuzzy_and",
    },
    Task {
        operation: Operation::FuzzyOr,
        corpus_type: CorpusType::FingerprintPairs,
        generator: "generate_fuzzy_or",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::SmilesRead),
        corpus_type: CorpusType::Smiles,
        generator: "generate_smiles_read",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::Sanitize),
        corpus_type: CorpusType::Smiles,
        generator: "generate_sanitize",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::Kekulize),
        corpus_type: CorpusType::Smiles,
        generator: "generate_kekulize",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::MolecularWeight),
        corpus_type: CorpusType::Smiles,
        generator: "generate_molecular_weight",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::ExactMolecularWeight),
        corpus_type: CorpusType::Smiles,
        generator: "generate_exact_molecular_weight",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::MolecularFormula),
        corpus_type: CorpusType::Smiles,
        generator: "generate_molecular_formula",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::NumHeavyAtoms),
        corpus_type: CorpusType::Smiles,
        generator: "generate_num_heavy_atoms",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::TotalAtomCount),
        corpus_type: CorpusType::Smiles,
        generator: "generate_total_atom_count",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::LipinskiHBA),
        corpus_type: CorpusType::Smiles,
        generator: "generate_lipinski_hba",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::LipinskiHBD),
        corpus_type: CorpusType::Smiles,
        generator: "generate_lipinski_hbd",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::FractionCSP3),
        corpus_type: CorpusType::Smiles,
        generator: "generate_fraction_csp3",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::AddHydrogens),
        corpus_type: CorpusType::Smiles,
        generator: "generate_add_hydrogens",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::RemoveHydrogens),
        corpus_type: CorpusType::Smiles,
        generator: "generate_remove_hydrogens",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::Coordinates2d),
        corpus_type: CorpusType::Smiles,
        generator: "generate_coordinates_2d",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::DistanceMatrix),
        corpus_type: CorpusType::Smiles,
        generator: "generate_distance_matrix",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::NumRings),
        corpus_type: CorpusType::Smiles,
        generator: "generate_num_rings",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::NumHeterocycles),
        corpus_type: CorpusType::Smiles,
        generator: "generate_num_heterocycles",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::NumAromaticRings),
        corpus_type: CorpusType::Smiles,
        generator: "generate_num_aromatic_rings",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::NumSaturatedRings),
        corpus_type: CorpusType::Smiles,
        generator: "generate_num_saturated_rings",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::NumAliphaticRings),
        corpus_type: CorpusType::Smiles,
        generator: "generate_num_aliphatic_rings",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::NumAromaticHeterocycles),
        corpus_type: CorpusType::Smiles,
        generator: "generate_num_aromatic_heterocycles",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::NumAromaticCarbocycles),
        corpus_type: CorpusType::Smiles,
        generator: "generate_num_aromatic_carbocycles",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::NumAliphaticHeterocycles),
        corpus_type: CorpusType::Smiles,
        generator: "generate_num_aliphatic_heterocycles",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::NumAliphaticCarbocycles),
        corpus_type: CorpusType::Smiles,
        generator: "generate_num_aliphatic_carbocycles",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::NumSaturatedHeterocycles),
        corpus_type: CorpusType::Smiles,
        generator: "generate_num_saturated_heterocycles",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::NumSaturatedCarbocycles),
        corpus_type: CorpusType::Smiles,
        generator: "generate_num_saturated_carbocycles",
    },
];

pub const RDKIT_VERSION: &str = "2026.03.1";

// Only typed integer indices/counts cross the oracle boundary.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Pair {
    pub id: String,
    pub length: u64,
    pub left: Vec<(u64, i32)>,
    pub right: Vec<(u64, i32)>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FingerprintInput {
    pub case: Pair,
    pub operation: Operation,
    pub width: Width,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FingerprintValue {
    pub length: u64,
    pub entries: Vec<(u64, i32)>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Record {
    pub input: Input,
    pub output: Value,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SmilesCase {
    pub id: String,
    pub smiles: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Input {
    Fingerprint(FingerprintInput),
    Molecular {
        case: SmilesCase,
        profile: molecule_plan::Profile,
    },
}

impl Input {
    pub fn task_name(&self) -> &'static str {
        match self {
            Self::Fingerprint(input) => input.operation.name(),
            Self::Molecular { profile, .. } => {
                use molecule_plan::Profile::*;
                match profile {
                    SmilesRead { .. } => "smiles_read",
                    SanitizeAll => "sanitize",
                    Kekulize { .. } => "kekulize",
                    MolecularWeight { .. } => "molecular_weight",
                    ExactMolecularWeight { .. } => "exact_molecular_weight",
                    MolecularFormula { .. } => "molecular_formula",
                    NumHeavyAtoms { .. } => "num_heavy_atoms",
                    TotalAtomCount { .. } => "total_atom_count",
                    LipinskiHBA { .. } => "lipinski_hba",
                    LipinskiHBD { .. } => "lipinski_hbd",
                    FractionCSP3 { .. } => "fraction_csp3",
                    NumRings { .. } => "num_rings",
                    NumHeterocycles { .. } => "num_heterocycles",
                    NumAromaticRings { .. } => "num_aromatic_rings",
                    NumSaturatedRings { .. } => "num_saturated_rings",
                    NumAliphaticRings { .. } => "num_aliphatic_rings",
                    NumAromaticHeterocycles { .. } => "num_aromatic_heterocycles",
                    NumAromaticCarbocycles { .. } => "num_aromatic_carbocycles",
                    NumAliphaticHeterocycles { .. } => "num_aliphatic_heterocycles",
                    NumAliphaticCarbocycles { .. } => "num_aliphatic_carbocycles",
                    NumSaturatedHeterocycles { .. } => "num_saturated_heterocycles",
                    NumSaturatedCarbocycles { .. } => "num_saturated_carbocycles",
                    AddHydrogens { .. } => "add_hydrogens",
                    RemoveHydrogens { .. } => "remove_hydrogens",
                    Coordinates2dDefault => "coordinates_2d",
                    CipLabels { .. } => "cip_labels",
                    PotentialStereo { .. } => "potential_stereo",
                    Valence { .. } => "valence",
                    DistanceMatrix { .. } => "distance_matrix",
                }
            }
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Value {
    Fingerprint(FingerprintValue),
    Molecular(crate::molecular::Outcome),
}

#[derive(Clone, Debug, Default)]
pub struct Corpus {
    pub fingerprints: Vec<Pair>,
    pub molecules: Vec<SmilesCase>,
}

impl Task {
    pub fn key(&self) -> String {
        format!("{}_{}", self.operation.name(), self.corpus_type.name())
    }
    pub fn validate_reference(&self, input: &Input, output: &Value) -> Result<(), String> {
        match (self.corpus_type, input, output) {
            (
                CorpusType::FingerprintPairs,
                Input::Fingerprint(input),
                Value::Fingerprint(output),
            ) => fingerprint::validate_output(input, output),
            (CorpusType::Smiles, Input::Molecular { profile, .. }, Value::Molecular(output)) => {
                crate::molecular::validate_output(profile, output)
            }
            _ => Err("reference input/output kind mismatch".into()),
        }
    }
    pub fn count(&self, cases: &Corpus) -> usize {
        match self.operation {
            Operation::Molecular(id) => cases.molecules.len() * id.profiles().len(),
            _ => cases.fingerprints.len() * fingerprint::WIDTHS.len(),
        }
    }
}

pub fn select(name: Option<&str>) -> Result<Vec<&'static Task>, String> {
    let selected: Vec<_> = TASKS
        .iter()
        .filter(|t| name.is_none_or(|n| n == t.key() || n == t.operation.name()))
        .collect();
    if selected.is_empty() {
        return Err(format!("unknown task: {name:?}"));
    }
    Ok(selected)
}

pub fn validate(corpus: &Corpus, tasks: &[&Task]) -> Result<(), String> {
    use std::collections::BTreeSet;
    if tasks.is_empty() {
        return Err("empty task selection".into());
    }
    for task in tasks {
        if task.count(corpus) == 0 {
            return Err(format!(
                "{}: corpus/profile selection is empty",
                task.operation.name()
            ));
        }
    }
    let mut molecule_ids = BTreeSet::new();
    for case in &corpus.molecules {
        if case.id.is_empty() || !molecule_ids.insert(&case.id) {
            return Err("empty/duplicate molecular case ID".into());
        }
    }
    if tasks
        .iter()
        .any(|t| t.corpus_type == CorpusType::FingerprintPairs)
    {
        fingerprint::validate(&corpus.fingerprints)?;
    }
    Ok(())
}

pub fn expand(cases: &Corpus, task: &Task) -> Vec<Input> {
    if let Operation::Molecular(id) = task.operation {
        return cases
            .molecules
            .iter()
            .flat_map(|case| {
                id.profiles()
                    .into_iter()
                    .map(move |profile| Input::Molecular {
                        case: case.clone(),
                        profile,
                    })
            })
            .collect();
    }
    fingerprint::expand(&cases.fingerprints, task.operation)
}

pub fn builtin() -> Vec<Pair> {
    // Each row names a source branch, not an arbitrary sampling size.
    [
        ("empty_both", vec![], vec![]),
        ("left_empty_right_tail", vec![], vec![(2, 3)]),
        ("right_empty_remove_left", vec![(2, -3)], vec![]),
        (
            "shared_signed_min_max",
            vec![(1, 5), (3, -2), (8, 4)],
            vec![(1, 3), (3, -4), (9, 7)],
        ),
        (
            "disjoint_interleaved",
            vec![(1, 2), (3, -1)],
            vec![(0, 4), (2, -3), (4, 5)],
        ),
        (
            "explicit_zero_shared_and_exclusive",
            vec![(1, 0), (3, 0)],
            vec![(1, -2), (4, 0)],
        ),
    ]
    .into_iter()
    .map(|(id, left, right)| Pair {
        id: id.into(),
        length: 16,
        left,
        right,
    })
    .collect()
}
