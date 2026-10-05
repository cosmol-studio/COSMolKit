//! Executable tasks, typed input families and parameter matrices.
//! Future catalog rows are not executable registration or parity claims.
pub mod fingerprint;
pub mod fingerprint_corpus;
pub mod molecule_plan;
use serde::{Deserialize, Serialize};

/// Fixed source-regression recipes, distinct from function/corpus tasks.
/// This is the same authoritative registry; no SMILES expansion is performed.
pub struct SpecialRegression {
    pub key: &'static str,
    pub fixture: &'static str,
    pub generator: &'static str,
    pub output: &'static str,
    pub rows: usize,
    pub schema: SpecialRegressionSchema,
    pub generator_dependencies: &'static [&'static str],
    pub takes_input: bool,
}

#[derive(Clone, Copy)]
pub enum SpecialRegressionSchema {
    StructureTags,
    TautomerBranches,
}

pub const SPECIAL_REGRESSIONS: &[SpecialRegression] = &[
    SpecialRegression {
        key: "assign_chiral_tags_from_structure",
        fixture: "testdata/stereo/fixtures/assign_atom_chiral_tags_from_structure_cases.json",
        generator: "tools/testdata/rdkit/_generate_tetrahedral_stereo_geometry.py",
        output: "assign_atom_chiral_tags_from_structure.jsonl",
        rows: 77,
        schema: SpecialRegressionSchema::StructureTags,
        generator_dependencies: &[],
        takes_input: false,
    },
    SpecialRegression {
        key: "tautomer_long_conjugated",
        fixture: "testdata/tautomer/fixtures/rdkit/long_conjugated_cases.json",
        generator: "tools/testdata/rdkit/_generate_tautomer_special_regression.py",
        output: "long_conjugated.jsonl",
        rows: 1,
        schema: SpecialRegressionSchema::TautomerBranches,
        generator_dependencies: &[
            "tools/testdata/rdkit/_tautomer_oracle.py",
            "tools/testdata/rdkit/tautomer_profile.json",
        ],
        takes_input: true,
    },
];

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Operation {
    SubstructureMatch,
    FuzzyAnd,
    FuzzyOr,
    Molecular(molecule_plan::TaskId),
    BioPdbOutput,
    UffCoverage,
    UffOptimization,
    UffConformerOptimization,
}

impl Operation {
    pub fn name(self) -> &'static str {
        match self {
            Self::SubstructureMatch => "substructure_match",
            Self::FuzzyAnd => "fuzzy_and",
            Self::FuzzyOr => "fuzzy_or",
            Self::Molecular(id) => id.name(),
            Self::BioPdbOutput => "bio_pdb_output",
            Self::UffCoverage => "uff_has_all_molecule_params",
            Self::UffOptimization => "uff_optimize",
            Self::UffConformerOptimization => "uff_optimize_conformers",
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
        operation: Operation::BioPdbOutput,
        corpus_type: CorpusType::Pdb,
        generator: "generate_bio_pdb_output_pdb",
    },
    Task {
        operation: Operation::BioPdbOutput,
        corpus_type: CorpusType::Cif,
        generator: "generate_bio_pdb_output_cif",
    },
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
        operation: Operation::Molecular(molecule_plan::TaskId::NumHeteroatoms),
        corpus_type: CorpusType::Smiles,
        generator: "generate_num_heteroatoms",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::NumHba),
        corpus_type: CorpusType::Smiles,
        generator: "generate_num_hba",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::NumHbd),
        corpus_type: CorpusType::Smiles,
        generator: "generate_num_hbd",
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
        operation: Operation::Molecular(molecule_plan::TaskId::Svg),
        corpus_type: CorpusType::Smiles,
        generator: "generate_svg",
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
    Task {
        operation: Operation::UffCoverage,
        corpus_type: CorpusType::Smiles,
        generator: "generate_uff_has_all_molecule_params",
    },
    Task {
        operation: Operation::UffOptimization,
        corpus_type: CorpusType::Smiles,
        generator: "generate_uff_optimize",
    },
    Task {
        operation: Operation::UffConformerOptimization,
        corpus_type: CorpusType::Smiles,
        generator: "generate_uff_optimize_conformers",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::MorganFingerprint),
        corpus_type: CorpusType::Smiles,
        generator: "generate_morgan_fingerprint",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::MorganSparseFingerprint),
        corpus_type: CorpusType::Smiles,
        generator: "generate_morgan_sparse_fingerprint",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::MorganCountFingerprint),
        corpus_type: CorpusType::Smiles,
        generator: "generate_morgan_count_fingerprint",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::MorganSparseCountFingerprint),
        corpus_type: CorpusType::Smiles,
        generator: "generate_morgan_sparse_count_fingerprint",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::Chi0),
        corpus_type: CorpusType::Smiles,
        generator: "generate_chi_0",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::Chi1),
        corpus_type: CorpusType::Smiles,
        generator: "generate_chi_1",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::HallKierAlpha),
        corpus_type: CorpusType::Smiles,
        generator: "generate_hall_kier_alpha",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::HallKierAlphaWithContributions),
        corpus_type: CorpusType::Smiles,
        generator: "generate_hall_kier_alpha_with_contributions",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::Kappa1),
        corpus_type: CorpusType::Smiles,
        generator: "generate_kappa_1",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::Kappa2),
        corpus_type: CorpusType::Smiles,
        generator: "generate_kappa_2",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::Kappa3),
        corpus_type: CorpusType::Smiles,
        generator: "generate_kappa_3",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::Phi),
        corpus_type: CorpusType::Smiles,
        generator: "generate_phi",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::Mqns),
        corpus_type: CorpusType::Smiles,
        generator: "generate_mqns",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::Chi0V),
        corpus_type: CorpusType::Smiles,
        generator: "generate_chi_0_v",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::Chi1V),
        corpus_type: CorpusType::Smiles,
        generator: "generate_chi_1_v",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::Chi2V),
        corpus_type: CorpusType::Smiles,
        generator: "generate_chi_2_v",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::Chi3V),
        corpus_type: CorpusType::Smiles,
        generator: "generate_chi_3_v",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::Chi4V),
        corpus_type: CorpusType::Smiles,
        generator: "generate_chi_4_v",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::Chi0N),
        corpus_type: CorpusType::Smiles,
        generator: "generate_chi_0_n",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::Chi1N),
        corpus_type: CorpusType::Smiles,
        generator: "generate_chi_1_n",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::Chi2N),
        corpus_type: CorpusType::Smiles,
        generator: "generate_chi_2_n",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::Chi3N),
        corpus_type: CorpusType::Smiles,
        generator: "generate_chi_3_n",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::Chi4N),
        corpus_type: CorpusType::Smiles,
        generator: "generate_chi_4_n",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::ChiNV),
        corpus_type: CorpusType::Smiles,
        generator: "generate_chi_n_v",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::ChiNN),
        corpus_type: CorpusType::Smiles,
        generator: "generate_chi_n_n",
    },
    Task {
        operation: Operation::SubstructureMatch,
        corpus_type: CorpusType::Smiles,
        generator: "generate_substructure_match",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::TautomerEnumeration),
        corpus_type: CorpusType::Smiles,
        generator: "generate_tautomer_enumeration",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::TautomerCanonicalization),
        corpus_type: CorpusType::Smiles,
        generator: "generate_tautomer_canonicalization",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::NumAmideBonds),
        corpus_type: CorpusType::Smiles,
        generator: "generate_num_amide_bonds",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::NumSpiroAtoms),
        corpus_type: CorpusType::Smiles,
        generator: "generate_num_spiro_atoms",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::NumBridgeheadAtoms),
        corpus_type: CorpusType::Smiles,
        generator: "generate_num_bridgehead_atoms",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::NumAtomStereoCenters),
        corpus_type: CorpusType::Smiles,
        generator: "generate_num_atom_stereo_centers",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::NumUnspecifiedAtomStereoCenters),
        corpus_type: CorpusType::Smiles,
        generator: "generate_num_unspecified_atom_stereo_centers",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::NumRotatableBonds),
        corpus_type: CorpusType::Smiles,
        generator: "generate_num_rotatable_bonds",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::CrippenDescriptors),
        corpus_type: CorpusType::Smiles,
        generator: "generate_crippen_descriptors",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::LabuteAsa),
        corpus_type: CorpusType::Smiles,
        generator: "generate_labute_asa",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::LabuteAsaContributions),
        corpus_type: CorpusType::Smiles,
        generator: "generate_labute_asa_contributions",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::Tpsa),
        corpus_type: CorpusType::Smiles,
        generator: "generate_tpsa",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::SlogpVsa),
        corpus_type: CorpusType::Smiles,
        generator: "generate_slogp_vsa",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::SmrVsa),
        corpus_type: CorpusType::Smiles,
        generator: "generate_smr_vsa",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::SlogpVsa1),
        corpus_type: CorpusType::Smiles,
        generator: "generate_slogp_vsa_1",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::SlogpVsa2),
        corpus_type: CorpusType::Smiles,
        generator: "generate_slogp_vsa_2",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::SlogpVsa3),
        corpus_type: CorpusType::Smiles,
        generator: "generate_slogp_vsa_3",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::SlogpVsa4),
        corpus_type: CorpusType::Smiles,
        generator: "generate_slogp_vsa_4",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::SlogpVsa5),
        corpus_type: CorpusType::Smiles,
        generator: "generate_slogp_vsa_5",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::SlogpVsa6),
        corpus_type: CorpusType::Smiles,
        generator: "generate_slogp_vsa_6",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::SlogpVsa7),
        corpus_type: CorpusType::Smiles,
        generator: "generate_slogp_vsa_7",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::SlogpVsa8),
        corpus_type: CorpusType::Smiles,
        generator: "generate_slogp_vsa_8",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::SlogpVsa9),
        corpus_type: CorpusType::Smiles,
        generator: "generate_slogp_vsa_9",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::SlogpVsa10),
        corpus_type: CorpusType::Smiles,
        generator: "generate_slogp_vsa_10",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::SlogpVsa11),
        corpus_type: CorpusType::Smiles,
        generator: "generate_slogp_vsa_11",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::SlogpVsa12),
        corpus_type: CorpusType::Smiles,
        generator: "generate_slogp_vsa_12",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::SmrVsa1),
        corpus_type: CorpusType::Smiles,
        generator: "generate_smr_vsa_1",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::SmrVsa2),
        corpus_type: CorpusType::Smiles,
        generator: "generate_smr_vsa_2",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::SmrVsa3),
        corpus_type: CorpusType::Smiles,
        generator: "generate_smr_vsa_3",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::SmrVsa4),
        corpus_type: CorpusType::Smiles,
        generator: "generate_smr_vsa_4",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::SmrVsa5),
        corpus_type: CorpusType::Smiles,
        generator: "generate_smr_vsa_5",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::SmrVsa6),
        corpus_type: CorpusType::Smiles,
        generator: "generate_smr_vsa_6",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::SmrVsa7),
        corpus_type: CorpusType::Smiles,
        generator: "generate_smr_vsa_7",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::SmrVsa8),
        corpus_type: CorpusType::Smiles,
        generator: "generate_smr_vsa_8",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::SmrVsa9),
        corpus_type: CorpusType::Smiles,
        generator: "generate_smr_vsa_9",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::SmrVsa10),
        corpus_type: CorpusType::Smiles,
        generator: "generate_smr_vsa_10",
    },
    Task {
        operation: Operation::Molecular(molecule_plan::TaskId::Qed),
        corpus_type: CorpusType::Smiles,
        generator: "generate_qed",
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
    Search(crate::search::SearchInput),
    Uff(crate::uff::UffInput),
    Fingerprint(FingerprintInput),
    Molecular {
        case: SmilesCase,
        profile: molecule_plan::Profile,
    },
    BioPdbOutput {
        case: BioPdbCase,
        profile: BioPdbOutputProfile,
    },
}

/// Typed BIO text case for PDB output parity (frozen §5).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct BioPdbCase {
    pub id: String,
    pub text: String,
    pub format: BioPdbCorpusFormat,
}

impl BioPdbCase {
    pub(crate) fn matches_corpus(&self, corpus: CorpusType) -> bool {
        matches!(
            (self.format, corpus),
            (BioPdbCorpusFormat::Pdb, CorpusType::Pdb) | (BioPdbCorpusFormat::Cif, CorpusType::Cif)
        )
    }
}

/// Corpus format for BIO PDB output tasks.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BioPdbCorpusFormat {
    Pdb,
    Cif,
}

/// The 32 source option profiles (frozen §5: 2^5 bool combinations).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct BioPdbOutputProfile {
    pub ter_records: bool,
    pub numbered_ter: bool,
    pub ter_ignores_type: bool,
    pub preserve_serial: bool,
    pub end_record: bool,
}

impl BioPdbOutputProfile {
    pub const ALL: [Self; 32] = [
        Self::f(false, false, false, false, false),
        Self::f(false, false, false, false, true),
        Self::f(false, false, false, true, false),
        Self::f(false, false, false, true, true),
        Self::f(false, false, true, false, false),
        Self::f(false, false, true, false, true),
        Self::f(false, false, true, true, false),
        Self::f(false, false, true, true, true),
        Self::f(false, true, false, false, false),
        Self::f(false, true, false, false, true),
        Self::f(false, true, false, true, false),
        Self::f(false, true, false, true, true),
        Self::f(false, true, true, false, false),
        Self::f(false, true, true, false, true),
        Self::f(false, true, true, true, false),
        Self::f(false, true, true, true, true),
        Self::f(true, false, false, false, false),
        Self::f(true, false, false, false, true),
        Self::f(true, false, false, true, false),
        Self::f(true, false, false, true, true),
        Self::f(true, false, true, false, false),
        Self::f(true, false, true, false, true),
        Self::f(true, false, true, true, false),
        Self::f(true, false, true, true, true),
        Self::f(true, true, false, false, false),
        Self::f(true, true, false, false, true),
        Self::f(true, true, false, true, false),
        Self::f(true, true, false, true, true),
        Self::f(true, true, true, false, false),
        Self::f(true, true, true, false, true),
        Self::f(true, true, true, true, false),
        Self::f(true, true, true, true, true),
    ];

    const fn f(ter: bool, num: bool, ign: bool, pres: bool, end: bool) -> Self {
        Self {
            ter_records: ter,
            numbered_ter: num,
            ter_ignores_type: ign,
            preserve_serial: pres,
            end_record: end,
        }
    }
}

impl Input {
    pub fn task_name(&self) -> &'static str {
        match self {
            Self::Search(_) => "substructure_match",
            Self::Uff(row) => match row.profile {
                crate::uff::Profile::Coverage { .. } => "uff_has_all_molecule_params",
                crate::uff::Profile::Optimization { .. } => "uff_optimize",
                crate::uff::Profile::ConformerOptimization { .. } => "uff_optimize_conformers",
            },
            Self::Fingerprint(input) => input.operation.name(),
            Self::BioPdbOutput { .. } => "bio_pdb_output",
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
                    NumHeteroatoms { .. } => "num_heteroatoms",
                    NumHba { .. } => "num_hba",
                    NumHbd { .. } => "num_hbd",
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
                    SvgDefault => "svg",
                    CipLabels { .. } => "cip_labels",
                    PotentialStereo { .. } => "potential_stereo",
                    Valence { .. } => "valence",
                    DistanceMatrix { .. } => "distance_matrix",
                    NumAmideBonds => "num_amide_bonds",
                    NumSpiroAtoms => "num_spiro_atoms",
                    NumBridgeheadAtoms => "num_bridgehead_atoms",
                    NumAtomStereoCenters => "num_atom_stereo_centers",
                    NumUnspecifiedAtomStereoCenters => "num_unspecified_atom_stereo_centers",
                    NumRotatableBonds { .. } => "num_rotatable_bonds",
                    CrippenDescriptors { .. } => "crippen_descriptors",
                    LabuteAsa { .. } => "labute_asa",
                    LabuteAsaContributions { .. } => "labute_asa_contributions",
                    Tpsa { .. } => "tpsa",
                    SlogpVsa { .. } => "slogp_vsa",
                    SmrVsa { .. } => "smr_vsa",
                    SlogpVsa1 => "slogp_vsa_1",
                    SlogpVsa2 => "slogp_vsa_2",
                    SlogpVsa3 => "slogp_vsa_3",
                    SlogpVsa4 => "slogp_vsa_4",
                    SlogpVsa5 => "slogp_vsa_5",
                    SlogpVsa6 => "slogp_vsa_6",
                    SlogpVsa7 => "slogp_vsa_7",
                    SlogpVsa8 => "slogp_vsa_8",
                    SlogpVsa9 => "slogp_vsa_9",
                    SlogpVsa10 => "slogp_vsa_10",
                    SlogpVsa11 => "slogp_vsa_11",
                    SlogpVsa12 => "slogp_vsa_12",
                    SmrVsa1 => "smr_vsa_1",
                    SmrVsa2 => "smr_vsa_2",
                    SmrVsa3 => "smr_vsa_3",
                    SmrVsa4 => "smr_vsa_4",
                    SmrVsa5 => "smr_vsa_5",
                    SmrVsa6 => "smr_vsa_6",
                    SmrVsa7 => "smr_vsa_7",
                    SmrVsa8 => "smr_vsa_8",
                    SmrVsa9 => "smr_vsa_9",
                    SmrVsa10 => "smr_vsa_10",
                    Qed => "qed",
                    Chi0VWithParams { .. } => "chi_0_v",
                    Chi1VWithParams { .. } => "chi_1_v",
                    Chi2VWithParams { .. } => "chi_2_v",
                    Chi3VWithParams { .. } => "chi_3_v",
                    Chi4VWithParams { .. } => "chi_4_v",
                    Chi0NWithParams { .. } => "chi_0_n",
                    Chi1NWithParams { .. } => "chi_1_n",
                    Chi2NWithParams { .. } => "chi_2_n",
                    Chi3NWithParams { .. } => "chi_3_n",
                    Chi4NWithParams { .. } => "chi_4_n",
                    ChiNVWithParams { .. } => "chi_n_v",
                    ChiNNWithParams { .. } => "chi_n_n",
                    LabuteAsaCacheSequence => "labute_asa",
                    SlogpVsaCacheSequence => "slogp_vsa",
                    SmrVsaCacheSequence => "smr_vsa",
                    ChiNVCacheSequence => "chi_n_v",
                    ChiNNCacheSequence => "chi_n_n",
                    Chi0 => "chi_0",
                    Chi1 => "chi_1",
                    HallKierAlpha => "hall_kier_alpha",
                    HallKierAlphaWithContributions => "hall_kier_alpha_with_contributions",
                    Kappa1 => "kappa_1",
                    Kappa2 => "kappa_2",
                    Kappa3 => "kappa_3",
                    Phi => "phi",
                    Mqns { .. } => "mqns",
                    Chi0V => "chi_0_v",
                    Chi1V => "chi_1_v",
                    Chi2V => "chi_2_v",
                    Chi3V => "chi_3_v",
                    Chi4V => "chi_4_v",
                    Chi0N => "chi_0_n",
                    Chi1N => "chi_1_n",
                    Chi2N => "chi_2_n",
                    Chi3N => "chi_3_n",
                    Chi4N => "chi_4_n",
                    ChiNV { .. } => "chi_n_v",
                    ChiNN { .. } => "chi_n_n",
                    TautomerEnumeration { .. } => "tautomer_enumeration",
                    TautomerCanonicalization { .. } => "tautomer_canonicalization",
                    Morgan { output, .. } => output.task_name(),
                }
            }
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Value {
    Search(crate::search::Outcome),
    Uff(crate::uff::Observation),
    Fingerprint(FingerprintValue),
    Molecular(crate::molecular::Outcome),
    BioPdbOutput(BioPdbOutputValue),
}

/// Text output from a BIO PDB write operation.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct BioPdbOutputValue {
    pub text: String,
    pub error: Option<BioPdbOutputError>,
}

/// Stage-typed error for BIO PDB output.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "stage", rename_all = "snake_case")]
pub enum BioPdbOutputError {
    Parse {
        format: BioPdbCorpusFormat,
        message: String,
    },
    Write {
        message: String,
    },
}

impl std::fmt::Display for BioPdbOutputError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Parse { format, message } => {
                write!(f, "parse {format:?}: {message}")
            }
            Self::Write { message } => write!(f, "write: {message}"),
        }
    }
}

/// Native reference identity for the BIO PDB output tasks.
pub const BIO_PDB_REFERENCE: ReferenceBackend = ReferenceBackend {
    library: "gemmi",
    version: "0.7.5",
    commit: "5cc1c23c6007e0e6cbd69289c6f7c0bff50e943e",
};

pub struct ReferenceBackend {
    pub library: &'static str,
    pub version: &'static str,
    pub commit: &'static str,
}

#[derive(Clone, Debug, Default)]
pub struct Corpus {
    pub fingerprints: Vec<Pair>,
    pub molecules: Vec<SmilesCase>,
    pub bio_cases: Vec<BioPdbCase>,
}

impl Task {
    pub fn key(&self) -> String {
        format!("{}_{}", self.operation.name(), self.corpus_type.name())
    }
    pub fn validate_reference(
        &self,
        recipe: &Input,
        prepared: &Input,
        output: &Value,
    ) -> Result<(), String> {
        if !matches!(recipe, Input::Uff(_)) && recipe != prepared {
            return Err("reference case/parameter mismatch".into());
        }
        if recipe.task_name() != self.operation.name()
            || prepared.task_name() != self.operation.name()
        {
            return Err("reference task/input mismatch".into());
        }
        match (recipe, prepared, output) {
            (Input::Search(recipe), Input::Search(prepared), Value::Search(output))
                if recipe == prepared =>
            {
                crate::search::validate_output(recipe, output)
            }
            (Input::Uff(recipe), Input::Uff(prepared), Value::Uff(output)) => {
                crate::uff::validate_reference(recipe, prepared, output)
            }
            (
                Input::Fingerprint(recipe),
                Input::Fingerprint(prepared),
                Value::Fingerprint(output),
            ) if recipe == prepared => fingerprint::validate_output(recipe, output),
            (
                Input::Molecular { .. },
                Input::Molecular { profile, .. },
                Value::Molecular(output),
            ) if recipe == prepared => crate::molecular::validate_output(profile, output),
            (
                Input::BioPdbOutput { .. },
                Input::BioPdbOutput { case, .. },
                Value::BioPdbOutput(_),
            ) if recipe == prepared => match self.corpus_type {
                CorpusType::Pdb if case.format == BioPdbCorpusFormat::Pdb => Ok(()),
                CorpusType::Cif if case.format == BioPdbCorpusFormat::Cif => Ok(()),
                CorpusType::Pdb => Err("pdb task requires pdb format case".into()),
                CorpusType::Cif => Err("cif task requires cif format case".into()),
                _ => Err("reference input/output kind mismatch".into()),
            },
            _ => Err("reference input/output kind mismatch".into()),
        }
    }
    pub fn count(&self, cases: &Corpus) -> usize {
        match self.operation {
            Operation::SubstructureMatch => cases.molecules.len() * crate::search::profiles().len(),
            Operation::UffCoverage
            | Operation::UffOptimization
            | Operation::UffConformerOptimization => {
                cases.molecules.len() * crate::uff::profiles(self.operation).len()
            }
            Operation::Molecular(id) => cases.molecules.len() * id.profiles().len(),
            Operation::BioPdbOutput => {
                cases
                    .bio_cases
                    .iter()
                    .filter(|case| case.matches_corpus(self.corpus_type))
                    .count()
                    * BioPdbOutputProfile::ALL.len()
            }
            _ => cases.fingerprints.len() * fingerprint::WIDTHS.len(),
        }
    }
}

pub fn select(name: Option<&str>) -> Result<Vec<&'static Task>, String> {
    let selected: Vec<_> = TASKS
        .iter()
        .filter(|t| name.is_none_or(|n| {
            n == t.key() || n == t.operation.name()
                || (n == "descriptors" && matches!(t.operation,
                    Operation::Molecular(id) if id.category() == molecule_plan::Category::Descriptors))
        }))
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
    let mut bio_ids = BTreeSet::new();
    for case in &corpus.bio_cases {
        if case.id.is_empty() || !bio_ids.insert((case.format as u8, &case.id)) {
            return Err("empty/duplicate BIO case ID within input family".into());
        }
    }
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
    if task.operation == Operation::SubstructureMatch {
        return cases
            .molecules
            .iter()
            .flat_map(|case| {
                crate::search::profiles().into_iter().map(move |profile| {
                    Input::Search(crate::search::SearchInput {
                        case: case.clone(),
                        profile,
                    })
                })
            })
            .collect();
    }
    if task.operation == Operation::BioPdbOutput {
        return cases
            .bio_cases
            .iter()
            .filter(|case| case.matches_corpus(task.corpus_type))
            .flat_map(|case| {
                BioPdbOutputProfile::ALL
                    .into_iter()
                    .map(move |profile| Input::BioPdbOutput {
                        case: case.clone(),
                        profile,
                    })
            })
            .collect();
    }
    if matches!(
        task.operation,
        Operation::UffCoverage | Operation::UffOptimization | Operation::UffConformerOptimization
    ) {
        return cases
            .molecules
            .iter()
            .flat_map(|case| {
                crate::uff::profiles(task.operation)
                    .into_iter()
                    .map(move |profile| {
                        Input::Uff(crate::uff::UffInput {
                            case: case.clone(),
                            profile,
                            preparation: None,
                        })
                    })
            })
            .collect();
    }
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
