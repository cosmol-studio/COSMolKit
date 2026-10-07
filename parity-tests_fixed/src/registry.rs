//! Executable tasks, typed input families and parameter matrices.
//! Declarations also generate Cargo test functions; no future-task catalog.
pub mod molecule_plan;
use serde::{Deserialize, Serialize};

pub struct SpecialRegression {
    pub key: &'static str,
    pub fixture: &'static str,
    pub rows: usize,
    pub schema: SpecialRegressionSchema,
}
#[derive(Clone, Copy)]
pub enum SpecialRegressionSchema {
    StructureTags,
    TautomerBranches,
}
pub const SPECIAL_REGRESSIONS: &[SpecialRegression] = &[
    SpecialRegression {
        key: "structure_tags",
        fixture: "special/structure_tags.json",
        rows: 77,
        schema: SpecialRegressionSchema::StructureTags,
    },
    SpecialRegression {
        key: "tautomer_long_conjugated",
        fixture: "special/tautomer_long_conjugated.json",
        rows: 1,
        schema: SpecialRegressionSchema::TautomerBranches,
    },
];

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Operation {
    Fingerprint(crate::fingerprints::Kind),
    SmilesWrite,
    SubstructureMatch,
    Molecular(molecule_plan::TaskId),
    BioPdbOutput,
    UffCoverage,
    UffOptimization,
    UffConformerOptimization,
    MmffCoverage,
    MmffOptimization,
    MmffConformerOptimization,
}

impl Operation {
    pub fn name(self) -> &'static str {
        match self {
            Self::Fingerprint(kind) => kind.name(),
            Self::SmilesWrite => "smiles_write",
            Self::SubstructureMatch => "substructure_match",
            Self::Molecular(id) => id.name(),
            Self::BioPdbOutput => "bio_pdb_output",
            Self::UffCoverage => "uff_has_all_molecule_params",
            Self::UffOptimization => "uff_optimize",
            Self::UffConformerOptimization => "uff_optimize_conformers",
            Self::MmffCoverage => "mmff_has_all_molecule_params",
            Self::MmffOptimization => "mmff_optimize",
            Self::MmffConformerOptimization => "mmff_optimize_conformers",
        }
    }
}

/// Input format is part of test identity, never inferred from a filename.
/// Only registered formats have executable loaders and tests.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CorpusType {
    Smiles,
    Pdb,
    Cif,
}

impl CorpusType {
    pub fn name(self) -> &'static str {
        match self {
            Self::Smiles => "smiles",
            Self::Pdb => "pdb",
            Self::Cif => "cif",
        }
    }
}

pub struct Task {
    pub key: &'static str,
    pub operation: Operation,
    pub corpus_type: CorpusType,
    pub generator: &'static str,
}

/// One declaration generates both adapters and Cargo test functions.
#[macro_export]
macro_rules! corpus_tasks {
    ($apply:ident) => {
        $apply! {
            (bio_pdb_output_pdb, Operation::BioPdbOutput, CorpusType::Pdb, "generate_bio_pdb_output_pdb"),
            (bio_pdb_output_cif, Operation::BioPdbOutput, CorpusType::Cif, "generate_bio_pdb_output_cif"),
            (smiles_read_smiles, Operation::Molecular(molecule_plan::TaskId::SmilesRead), CorpusType::Smiles, "generate_smiles_read"),
            (smiles_write_smiles, Operation::SmilesWrite, CorpusType::Smiles, "generate_smiles_write"),
            (sanitize_smiles, Operation::Molecular(molecule_plan::TaskId::Sanitize), CorpusType::Smiles, "generate_sanitize"),
            (kekulize_smiles, Operation::Molecular(molecule_plan::TaskId::Kekulize), CorpusType::Smiles, "generate_kekulize"),
            (molecular_weight_smiles, Operation::Molecular(molecule_plan::TaskId::MolecularWeight), CorpusType::Smiles, "generate_molecular_weight"),
            (exact_molecular_weight_smiles, Operation::Molecular(molecule_plan::TaskId::ExactMolecularWeight), CorpusType::Smiles, "generate_exact_molecular_weight"),
            (molecular_formula_smiles, Operation::Molecular(molecule_plan::TaskId::MolecularFormula), CorpusType::Smiles, "generate_molecular_formula"),
            (num_heavy_atoms_smiles, Operation::Molecular(molecule_plan::TaskId::NumHeavyAtoms), CorpusType::Smiles, "generate_num_heavy_atoms"),
            (total_atom_count_smiles, Operation::Molecular(molecule_plan::TaskId::TotalAtomCount), CorpusType::Smiles, "generate_total_atom_count"),
            (lipinski_hba_smiles, Operation::Molecular(molecule_plan::TaskId::LipinskiHBA), CorpusType::Smiles, "generate_lipinski_hba"),
            (lipinski_hbd_smiles, Operation::Molecular(molecule_plan::TaskId::LipinskiHBD), CorpusType::Smiles, "generate_lipinski_hbd"),
            (fraction_csp3_smiles, Operation::Molecular(molecule_plan::TaskId::FractionCSP3), CorpusType::Smiles, "generate_fraction_csp3"),
            (num_heteroatoms_smiles, Operation::Molecular(molecule_plan::TaskId::NumHeteroatoms), CorpusType::Smiles, "generate_num_heteroatoms"),
            (num_hba_smiles, Operation::Molecular(molecule_plan::TaskId::NumHba), CorpusType::Smiles, "generate_num_hba"),
            (num_hbd_smiles, Operation::Molecular(molecule_plan::TaskId::NumHbd), CorpusType::Smiles, "generate_num_hbd"),
            (add_hydrogens_smiles, Operation::Molecular(molecule_plan::TaskId::AddHydrogens), CorpusType::Smiles, "generate_add_hydrogens"),
            (remove_hydrogens_smiles, Operation::Molecular(molecule_plan::TaskId::RemoveHydrogens), CorpusType::Smiles, "generate_remove_hydrogens"),
            (coordinates_2d_smiles, Operation::Molecular(molecule_plan::TaskId::Coordinates2d), CorpusType::Smiles, "generate_coordinates_2d"),
            (svg_smiles, Operation::Molecular(molecule_plan::TaskId::Svg), CorpusType::Smiles, "generate_svg"),
            (distance_matrix_smiles, Operation::Molecular(molecule_plan::TaskId::DistanceMatrix), CorpusType::Smiles, "generate_distance_matrix"),
            (num_rings_smiles, Operation::Molecular(molecule_plan::TaskId::NumRings), CorpusType::Smiles, "generate_num_rings"),
            (num_heterocycles_smiles, Operation::Molecular(molecule_plan::TaskId::NumHeterocycles), CorpusType::Smiles, "generate_num_heterocycles"),
            (num_aromatic_rings_smiles, Operation::Molecular(molecule_plan::TaskId::NumAromaticRings), CorpusType::Smiles, "generate_num_aromatic_rings"),
            (num_saturated_rings_smiles, Operation::Molecular(molecule_plan::TaskId::NumSaturatedRings), CorpusType::Smiles, "generate_num_saturated_rings"),
            (num_aliphatic_rings_smiles, Operation::Molecular(molecule_plan::TaskId::NumAliphaticRings), CorpusType::Smiles, "generate_num_aliphatic_rings"),
            (num_aromatic_heterocycles_smiles, Operation::Molecular(molecule_plan::TaskId::NumAromaticHeterocycles), CorpusType::Smiles, "generate_num_aromatic_heterocycles"),
            (num_aromatic_carbocycles_smiles, Operation::Molecular(molecule_plan::TaskId::NumAromaticCarbocycles), CorpusType::Smiles, "generate_num_aromatic_carbocycles"),
            (num_aliphatic_heterocycles_smiles, Operation::Molecular(molecule_plan::TaskId::NumAliphaticHeterocycles), CorpusType::Smiles, "generate_num_aliphatic_heterocycles"),
            (num_aliphatic_carbocycles_smiles, Operation::Molecular(molecule_plan::TaskId::NumAliphaticCarbocycles), CorpusType::Smiles, "generate_num_aliphatic_carbocycles"),
            (num_saturated_heterocycles_smiles, Operation::Molecular(molecule_plan::TaskId::NumSaturatedHeterocycles), CorpusType::Smiles, "generate_num_saturated_heterocycles"),
            (num_saturated_carbocycles_smiles, Operation::Molecular(molecule_plan::TaskId::NumSaturatedCarbocycles), CorpusType::Smiles, "generate_num_saturated_carbocycles"),
            (uff_has_all_molecule_params_smiles, Operation::UffCoverage, CorpusType::Smiles, "generate_uff_has_all_molecule_params"),
            (uff_optimize_smiles, Operation::UffOptimization, CorpusType::Smiles, "generate_uff_optimize"),
            (uff_optimize_conformers_smiles, Operation::UffConformerOptimization, CorpusType::Smiles, "generate_uff_optimize_conformers"),
            (mmff_has_all_molecule_params_smiles, Operation::MmffCoverage, CorpusType::Smiles, "generate_mmff_has_all_molecule_params"),
            (mmff_optimize_smiles, Operation::MmffOptimization, CorpusType::Smiles, "generate_mmff_optimize"),
            (mmff_optimize_conformers_smiles, Operation::MmffConformerOptimization, CorpusType::Smiles, "generate_mmff_optimize_conformers"),
            (morgan_fingerprint_smiles, Operation::Molecular(molecule_plan::TaskId::MorganFingerprint), CorpusType::Smiles, "generate_morgan_fingerprint"),
            (maccs_fingerprint_smiles, Operation::Fingerprint($crate::fingerprints::Kind::Maccs), CorpusType::Smiles, "generate_fingerprint"),
            (topological_fingerprint_smiles, Operation::Fingerprint($crate::fingerprints::Kind::Topological), CorpusType::Smiles, "generate_fingerprint"),
            (layered_fingerprint_smiles, Operation::Fingerprint($crate::fingerprints::Kind::Layered), CorpusType::Smiles, "generate_fingerprint"),
            (pattern_fingerprint_smiles, Operation::Fingerprint($crate::fingerprints::Kind::Pattern), CorpusType::Smiles, "generate_fingerprint"),
            (fuzzy_and_smiles, Operation::Fingerprint($crate::fingerprints::Kind::FuzzyAnd), CorpusType::Smiles, "generate_fingerprint"),
            (fuzzy_or_smiles, Operation::Fingerprint($crate::fingerprints::Kind::FuzzyOr), CorpusType::Smiles, "generate_fingerprint"),
            (morgan_sparse_fingerprint_smiles, Operation::Molecular(molecule_plan::TaskId::MorganSparseFingerprint), CorpusType::Smiles, "generate_morgan_sparse_fingerprint"),
            (morgan_count_fingerprint_smiles, Operation::Molecular(molecule_plan::TaskId::MorganCountFingerprint), CorpusType::Smiles, "generate_morgan_count_fingerprint"),
            (morgan_sparse_count_fingerprint_smiles, Operation::Molecular(molecule_plan::TaskId::MorganSparseCountFingerprint), CorpusType::Smiles, "generate_morgan_sparse_count_fingerprint"),
            (chi_0_smiles, Operation::Molecular(molecule_plan::TaskId::Chi0), CorpusType::Smiles, "generate_chi_0"),
            (chi_1_smiles, Operation::Molecular(molecule_plan::TaskId::Chi1), CorpusType::Smiles, "generate_chi_1"),
            (hall_kier_alpha_smiles, Operation::Molecular(molecule_plan::TaskId::HallKierAlpha), CorpusType::Smiles, "generate_hall_kier_alpha"),
            (hall_kier_alpha_with_contributions_smiles, Operation::Molecular(molecule_plan::TaskId::HallKierAlphaWithContributions), CorpusType::Smiles, "generate_hall_kier_alpha_with_contributions"),
            (kappa_1_smiles, Operation::Molecular(molecule_plan::TaskId::Kappa1), CorpusType::Smiles, "generate_kappa_1"),
            (kappa_2_smiles, Operation::Molecular(molecule_plan::TaskId::Kappa2), CorpusType::Smiles, "generate_kappa_2"),
            (kappa_3_smiles, Operation::Molecular(molecule_plan::TaskId::Kappa3), CorpusType::Smiles, "generate_kappa_3"),
            (phi_smiles, Operation::Molecular(molecule_plan::TaskId::Phi), CorpusType::Smiles, "generate_phi"),
            (mqns_smiles, Operation::Molecular(molecule_plan::TaskId::Mqns), CorpusType::Smiles, "generate_mqns"),
            (chi_0_v_smiles, Operation::Molecular(molecule_plan::TaskId::Chi0V), CorpusType::Smiles, "generate_chi_0_v"),
            (chi_1_v_smiles, Operation::Molecular(molecule_plan::TaskId::Chi1V), CorpusType::Smiles, "generate_chi_1_v"),
            (chi_2_v_smiles, Operation::Molecular(molecule_plan::TaskId::Chi2V), CorpusType::Smiles, "generate_chi_2_v"),
            (chi_3_v_smiles, Operation::Molecular(molecule_plan::TaskId::Chi3V), CorpusType::Smiles, "generate_chi_3_v"),
            (chi_4_v_smiles, Operation::Molecular(molecule_plan::TaskId::Chi4V), CorpusType::Smiles, "generate_chi_4_v"),
            (chi_0_n_smiles, Operation::Molecular(molecule_plan::TaskId::Chi0N), CorpusType::Smiles, "generate_chi_0_n"),
            (chi_1_n_smiles, Operation::Molecular(molecule_plan::TaskId::Chi1N), CorpusType::Smiles, "generate_chi_1_n"),
            (chi_2_n_smiles, Operation::Molecular(molecule_plan::TaskId::Chi2N), CorpusType::Smiles, "generate_chi_2_n"),
            (chi_3_n_smiles, Operation::Molecular(molecule_plan::TaskId::Chi3N), CorpusType::Smiles, "generate_chi_3_n"),
            (chi_4_n_smiles, Operation::Molecular(molecule_plan::TaskId::Chi4N), CorpusType::Smiles, "generate_chi_4_n"),
            (chi_n_v_smiles, Operation::Molecular(molecule_plan::TaskId::ChiNV), CorpusType::Smiles, "generate_chi_n_v"),
            (chi_n_n_smiles, Operation::Molecular(molecule_plan::TaskId::ChiNN), CorpusType::Smiles, "generate_chi_n_n"),
            (substructure_match_smiles, Operation::SubstructureMatch, CorpusType::Smiles, "generate_substructure_match"),
            (tautomer_enumeration_smiles, Operation::Molecular(molecule_plan::TaskId::TautomerEnumeration), CorpusType::Smiles, "generate_tautomer_enumeration"),
            (tautomer_canonicalization_smiles, Operation::Molecular(molecule_plan::TaskId::TautomerCanonicalization), CorpusType::Smiles, "generate_tautomer_canonicalization"),
            (num_amide_bonds_smiles, Operation::Molecular(molecule_plan::TaskId::NumAmideBonds), CorpusType::Smiles, "generate_num_amide_bonds"),
            (num_spiro_atoms_smiles, Operation::Molecular(molecule_plan::TaskId::NumSpiroAtoms), CorpusType::Smiles, "generate_num_spiro_atoms"),
            (num_bridgehead_atoms_smiles, Operation::Molecular(molecule_plan::TaskId::NumBridgeheadAtoms), CorpusType::Smiles, "generate_num_bridgehead_atoms"),
            (num_atom_stereo_centers_smiles, Operation::Molecular(molecule_plan::TaskId::NumAtomStereoCenters), CorpusType::Smiles, "generate_num_atom_stereo_centers"),
            (num_unspecified_atom_stereo_centers_smiles, Operation::Molecular(molecule_plan::TaskId::NumUnspecifiedAtomStereoCenters), CorpusType::Smiles, "generate_num_unspecified_atom_stereo_centers"),
            (num_rotatable_bonds_smiles, Operation::Molecular(molecule_plan::TaskId::NumRotatableBonds), CorpusType::Smiles, "generate_num_rotatable_bonds"),
            (crippen_descriptors_smiles, Operation::Molecular(molecule_plan::TaskId::CrippenDescriptors), CorpusType::Smiles, "generate_crippen_descriptors"),
            (labute_asa_smiles, Operation::Molecular(molecule_plan::TaskId::LabuteAsa), CorpusType::Smiles, "generate_labute_asa"),
            (labute_asa_contributions_smiles, Operation::Molecular(molecule_plan::TaskId::LabuteAsaContributions), CorpusType::Smiles, "generate_labute_asa_contributions"),
            (tpsa_smiles, Operation::Molecular(molecule_plan::TaskId::Tpsa), CorpusType::Smiles, "generate_tpsa"),
            (slogp_vsa_smiles, Operation::Molecular(molecule_plan::TaskId::SlogpVsa), CorpusType::Smiles, "generate_slogp_vsa"),
            (smr_vsa_smiles, Operation::Molecular(molecule_plan::TaskId::SmrVsa), CorpusType::Smiles, "generate_smr_vsa"),
            (slogp_vsa_1_smiles, Operation::Molecular(molecule_plan::TaskId::SlogpVsa1), CorpusType::Smiles, "generate_slogp_vsa_1"),
            (slogp_vsa_2_smiles, Operation::Molecular(molecule_plan::TaskId::SlogpVsa2), CorpusType::Smiles, "generate_slogp_vsa_2"),
            (slogp_vsa_3_smiles, Operation::Molecular(molecule_plan::TaskId::SlogpVsa3), CorpusType::Smiles, "generate_slogp_vsa_3"),
            (slogp_vsa_4_smiles, Operation::Molecular(molecule_plan::TaskId::SlogpVsa4), CorpusType::Smiles, "generate_slogp_vsa_4"),
            (slogp_vsa_5_smiles, Operation::Molecular(molecule_plan::TaskId::SlogpVsa5), CorpusType::Smiles, "generate_slogp_vsa_5"),
            (slogp_vsa_6_smiles, Operation::Molecular(molecule_plan::TaskId::SlogpVsa6), CorpusType::Smiles, "generate_slogp_vsa_6"),
            (slogp_vsa_7_smiles, Operation::Molecular(molecule_plan::TaskId::SlogpVsa7), CorpusType::Smiles, "generate_slogp_vsa_7"),
            (slogp_vsa_8_smiles, Operation::Molecular(molecule_plan::TaskId::SlogpVsa8), CorpusType::Smiles, "generate_slogp_vsa_8"),
            (slogp_vsa_9_smiles, Operation::Molecular(molecule_plan::TaskId::SlogpVsa9), CorpusType::Smiles, "generate_slogp_vsa_9"),
            (slogp_vsa_10_smiles, Operation::Molecular(molecule_plan::TaskId::SlogpVsa10), CorpusType::Smiles, "generate_slogp_vsa_10"),
            (slogp_vsa_11_smiles, Operation::Molecular(molecule_plan::TaskId::SlogpVsa11), CorpusType::Smiles, "generate_slogp_vsa_11"),
            (slogp_vsa_12_smiles, Operation::Molecular(molecule_plan::TaskId::SlogpVsa12), CorpusType::Smiles, "generate_slogp_vsa_12"),
            (smr_vsa_1_smiles, Operation::Molecular(molecule_plan::TaskId::SmrVsa1), CorpusType::Smiles, "generate_smr_vsa_1"),
            (smr_vsa_2_smiles, Operation::Molecular(molecule_plan::TaskId::SmrVsa2), CorpusType::Smiles, "generate_smr_vsa_2"),
            (smr_vsa_3_smiles, Operation::Molecular(molecule_plan::TaskId::SmrVsa3), CorpusType::Smiles, "generate_smr_vsa_3"),
            (smr_vsa_4_smiles, Operation::Molecular(molecule_plan::TaskId::SmrVsa4), CorpusType::Smiles, "generate_smr_vsa_4"),
            (smr_vsa_5_smiles, Operation::Molecular(molecule_plan::TaskId::SmrVsa5), CorpusType::Smiles, "generate_smr_vsa_5"),
            (smr_vsa_6_smiles, Operation::Molecular(molecule_plan::TaskId::SmrVsa6), CorpusType::Smiles, "generate_smr_vsa_6"),
            (smr_vsa_7_smiles, Operation::Molecular(molecule_plan::TaskId::SmrVsa7), CorpusType::Smiles, "generate_smr_vsa_7"),
            (smr_vsa_8_smiles, Operation::Molecular(molecule_plan::TaskId::SmrVsa8), CorpusType::Smiles, "generate_smr_vsa_8"),
            (smr_vsa_9_smiles, Operation::Molecular(molecule_plan::TaskId::SmrVsa9), CorpusType::Smiles, "generate_smr_vsa_9"),
            (smr_vsa_10_smiles, Operation::Molecular(molecule_plan::TaskId::SmrVsa10), CorpusType::Smiles, "generate_smr_vsa_10"),
            (qed_smiles, Operation::Molecular(molecule_plan::TaskId::Qed), CorpusType::Smiles, "generate_qed"),
        }
        $apply! { composition: batch_smiles }
    };
}
macro_rules! define_tasks {
    (composition: $($key:ident),*) => {};
    ($(($key:ident, $operation:expr, $corpus:expr, $generator:literal)),* $(,)?) => {
        pub const TASKS: &[Task] = &[$(Task { key: stringify!($key), operation: $operation, corpus_type: $corpus, generator: $generator }),*];
    };
}
crate::corpus_tasks!(define_tasks);

pub const RDKIT_VERSION: &str = "2026.03.1";

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
    Fingerprint(crate::fingerprints::FingerprintInput),
    SmilesWrite(crate::smiles_write::WriteInput),
    Search(crate::search::SearchInput),
    Uff(crate::uff::UffInput),
    Mmff(crate::mmff::MmffInput),
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
            Self::Fingerprint(row) => row.task_name(),
            Self::SmilesWrite(_) => "smiles_write",
            Self::Mmff(row) => row.profile.task_name(),
            Self::Search(_) => "substructure_match",
            Self::Uff(row) => match row.profile {
                crate::uff::Profile::Coverage { .. } => "uff_has_all_molecule_params",
                crate::uff::Profile::Optimization { .. } => "uff_optimize",
                crate::uff::Profile::ConformerOptimization { .. } => "uff_optimize_conformers",
            },
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
    Fingerprint(crate::fingerprints::Observation),
    SmilesWrite(crate::smiles_write::Outcome),
    Search(crate::search::Outcome),
    Uff(crate::uff::Observation),
    Mmff(crate::mmff::Observation),
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
    pub molecules: Vec<SmilesCase>,
    pub bio_cases: Vec<BioPdbCase>,
}

impl Task {
    pub fn validate_reference(
        &self,
        recipe: &Input,
        prepared: &Input,
        output: &Value,
    ) -> Result<(), String> {
        if !matches!(recipe, Input::Uff(_) | Input::Mmff(_)) && recipe != prepared {
            return Err("reference case/parameter mismatch".into());
        }
        if recipe.task_name() != self.operation.name()
            || prepared.task_name() != self.operation.name()
        {
            return Err("reference task/input mismatch".into());
        }
        match (recipe, prepared, output) {
            (Input::SmilesWrite(recipe), Input::SmilesWrite(prepared), Value::SmilesWrite(_))
                if recipe == prepared =>
            {
                Ok(())
            }
            (
                Input::Fingerprint(recipe),
                Input::Fingerprint(prepared),
                Value::Fingerprint(output),
            ) if recipe == prepared => crate::fingerprints::validate_output(recipe, output),
            (Input::Search(recipe), Input::Search(prepared), Value::Search(output))
                if recipe == prepared =>
            {
                crate::search::validate_output(recipe, output)
            }
            (Input::Uff(recipe), Input::Uff(prepared), Value::Uff(output)) => {
                crate::uff::validate_reference(recipe, prepared, output)
            }
            (Input::Mmff(recipe), Input::Mmff(prepared), Value::Mmff(output)) => {
                crate::mmff::validate_reference(recipe, prepared, output)
            }
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
            Operation::Fingerprint(_) => cases.molecules.len(),
            Operation::SmilesWrite => cases.molecules.len() * crate::smiles_write::profiles().len(),
            Operation::SubstructureMatch => cases.molecules.len() * crate::search::profiles().len(),
            Operation::UffCoverage
            | Operation::UffOptimization
            | Operation::UffConformerOptimization => {
                cases.molecules.len() * crate::uff::profiles(self.operation).len()
            }
            Operation::Molecular(id) => cases.molecules.len() * id.profiles().len(),
            Operation::MmffCoverage
            | Operation::MmffOptimization
            | Operation::MmffConformerOptimization => {
                cases.molecules.len() * crate::mmff::profiles(self.operation).len()
            }
            Operation::BioPdbOutput => {
                cases
                    .bio_cases
                    .iter()
                    .filter(|case| case.matches_corpus(self.corpus_type))
                    .count()
                    * BioPdbOutputProfile::ALL.len()
            }
        }
    }
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
    Ok(())
}

pub fn expand(cases: &Corpus, task: &Task) -> Vec<Input> {
    if let Operation::Fingerprint(kind) = task.operation {
        return crate::fingerprints::inputs(&cases.molecules, kind);
    }
    if matches!(
        task.operation,
        Operation::MmffCoverage
            | Operation::MmffOptimization
            | Operation::MmffConformerOptimization
    ) {
        return cases
            .molecules
            .iter()
            .flat_map(|case| {
                crate::mmff::profiles(task.operation)
                    .into_iter()
                    .map(move |profile| {
                        Input::Mmff(crate::mmff::MmffInput {
                            case: case.clone(),
                            profile,
                            preparation: None,
                        })
                    })
            })
            .collect();
    }
    if task.operation == Operation::SmilesWrite {
        let profiles = crate::smiles_write::profiles();
        return cases
            .molecules
            .iter()
            .flat_map(|case| {
                profiles.iter().map(move |profile| {
                    Input::SmilesWrite(crate::smiles_write::WriteInput {
                        case: case.clone(),
                        profile: *profile,
                    })
                })
            })
            .collect();
    }
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
    let Operation::Molecular(id) = task.operation else {
        unreachable!("all non-molecular operations were expanded above");
    };
    cases
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
        .collect()
}
