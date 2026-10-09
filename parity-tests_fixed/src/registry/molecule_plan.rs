//! Molecular parameter matrices and comparison contracts.

/// Stable behavior categories; readiness never changes category.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Category {
    Notation,
    MolecularIo,
    Chemistry,
    Stereo,
    Descriptors,
    FingerprintGeneration,
    FingerprintValues,
    Query,
    Depiction,
    Conformers,
    ForceFields,
    Alignment,
    Tautomers,
    Identifiers,
    Bio,
    Composition,
    Native,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Reference {
    Rdkit,
    Gemmi,
    RustScalar,
    Roundtrip,
    NoExternalReference,
}

impl TaskId {
    pub const fn name(self) -> &'static str {
        match self {
            Self::TautomerEnumeration => "tautomer_enumeration",
            Self::TautomerCanonicalization => "tautomer_canonicalization",
            Self::Chi0 => "chi_0",
            Self::Chi1 => "chi_1",
            Self::HallKierAlpha => "hall_kier_alpha",
            Self::HallKierAlphaWithContributions => "hall_kier_alpha_with_contributions",
            Self::Kappa1 => "kappa_1",
            Self::Kappa2 => "kappa_2",
            Self::Kappa3 => "kappa_3",
            Self::Phi => "phi",
            Self::Mqns => "mqns",
            Self::Chi0V => "chi_0_v",
            Self::Chi1V => "chi_1_v",
            Self::Chi2V => "chi_2_v",
            Self::Chi3V => "chi_3_v",
            Self::Chi4V => "chi_4_v",
            Self::Chi0N => "chi_0_n",
            Self::Chi1N => "chi_1_n",
            Self::Chi2N => "chi_2_n",
            Self::Chi3N => "chi_3_n",
            Self::Chi4N => "chi_4_n",
            Self::ChiNV => "chi_n_v",
            Self::ChiNN => "chi_n_n",
            Self::NumAmideBonds => "num_amide_bonds",
            Self::NumSpiroAtoms => "num_spiro_atoms",
            Self::NumBridgeheadAtoms => "num_bridgehead_atoms",
            Self::NumAtomStereoCenters => "num_atom_stereo_centers",
            Self::NumUnspecifiedAtomStereoCenters => "num_unspecified_atom_stereo_centers",
            Self::NumRotatableBonds => "num_rotatable_bonds",
            Self::CrippenDescriptors => "crippen_descriptors",
            Self::LabuteAsa => "labute_asa",
            Self::LabuteAsaContributions => "labute_asa_contributions",
            Self::Tpsa => "tpsa",
            Self::SlogpVsa => "slogp_vsa",
            Self::SmrVsa => "smr_vsa",
            Self::SlogpVsa1 => "slogp_vsa_1",
            Self::SlogpVsa2 => "slogp_vsa_2",
            Self::SlogpVsa3 => "slogp_vsa_3",
            Self::SlogpVsa4 => "slogp_vsa_4",
            Self::SlogpVsa5 => "slogp_vsa_5",
            Self::SlogpVsa6 => "slogp_vsa_6",
            Self::SlogpVsa7 => "slogp_vsa_7",
            Self::SlogpVsa8 => "slogp_vsa_8",
            Self::SlogpVsa9 => "slogp_vsa_9",
            Self::SlogpVsa10 => "slogp_vsa_10",
            Self::SlogpVsa11 => "slogp_vsa_11",
            Self::SlogpVsa12 => "slogp_vsa_12",
            Self::SmrVsa1 => "smr_vsa_1",
            Self::SmrVsa2 => "smr_vsa_2",
            Self::SmrVsa3 => "smr_vsa_3",
            Self::SmrVsa4 => "smr_vsa_4",
            Self::SmrVsa5 => "smr_vsa_5",
            Self::SmrVsa6 => "smr_vsa_6",
            Self::SmrVsa7 => "smr_vsa_7",
            Self::SmrVsa8 => "smr_vsa_8",
            Self::SmrVsa9 => "smr_vsa_9",
            Self::SmrVsa10 => "smr_vsa_10",
            Self::Qed => "qed",
            Self::SmilesRead => "smiles_read",
            Self::Sanitize => "sanitize",
            Self::Kekulize => "kekulize",
            Self::MolecularWeight => "molecular_weight",
            Self::ExactMolecularWeight => "exact_molecular_weight",
            Self::MolecularFormula => "molecular_formula",
            Self::AddHydrogens => "add_hydrogens",
            Self::RemoveHydrogens => "remove_hs",
            Self::CipLabels => "cip_labels",
            Self::PotentialStereo => "potential_stereo",
            Self::Coordinates2d => "coordinates_2d",
            Self::Svg => "svg",
            Self::Valence => "valence",
            Self::DistanceMatrix => "distance_matrix",
            Self::NumHeavyAtoms => "num_heavy_atoms",
            Self::TotalAtomCount => "total_atom_count",
            Self::LipinskiHBA => "lipinski_hba",
            Self::LipinskiHBD => "lipinski_hbd",
            Self::FractionCSP3 => "fraction_csp3",
            Self::NumHeteroatoms => "num_heteroatoms",
            Self::NumHba => "num_hba",
            Self::NumHbd => "num_hbd",
            Self::NumRings => "num_rings",
            Self::NumHeterocycles => "num_heterocycles",
            Self::NumAromaticRings => "num_aromatic_rings",
            Self::NumSaturatedRings => "num_saturated_rings",
            Self::NumAliphaticRings => "num_aliphatic_rings",
            Self::NumAromaticHeterocycles => "num_aromatic_heterocycles",
            Self::NumAromaticCarbocycles => "num_aromatic_carbocycles",
            Self::NumAliphaticHeterocycles => "num_aliphatic_heterocycles",
            Self::NumAliphaticCarbocycles => "num_aliphatic_carbocycles",
            Self::NumSaturatedHeterocycles => "num_saturated_heterocycles",
            Self::NumSaturatedCarbocycles => "num_saturated_carbocycles",
            Self::MorganFingerprint => "fingerprint_morgan",
            Self::MorganSparseFingerprint => "fingerprint_morgan_sparse",
            Self::MorganCountFingerprint => "fingerprint_morgan_count",
            Self::MorganSparseCountFingerprint => "fingerprint_morgan_sparse_count",
        }
    }
    pub const fn category(self) -> Category {
        match self {
            Self::TautomerEnumeration | Self::TautomerCanonicalization => Category::Tautomers,
            Self::Chi0
            | Self::Chi1
            | Self::HallKierAlpha
            | Self::HallKierAlphaWithContributions
            | Self::Kappa1
            | Self::Kappa2
            | Self::Kappa3
            | Self::Phi
            | Self::Mqns
            | Self::Chi0V
            | Self::Chi1V
            | Self::Chi2V
            | Self::Chi3V
            | Self::Chi4V
            | Self::Chi0N
            | Self::Chi1N
            | Self::Chi2N
            | Self::Chi3N
            | Self::Chi4N
            | Self::ChiNV
            | Self::ChiNN => Category::Descriptors,
            Self::NumAmideBonds
            | Self::NumSpiroAtoms
            | Self::NumBridgeheadAtoms
            | Self::NumAtomStereoCenters
            | Self::NumUnspecifiedAtomStereoCenters
            | Self::NumRotatableBonds
            | Self::CrippenDescriptors
            | Self::LabuteAsa
            | Self::LabuteAsaContributions
            | Self::Tpsa
            | Self::SlogpVsa
            | Self::SmrVsa
            | Self::SlogpVsa1
            | Self::SlogpVsa2
            | Self::SlogpVsa3
            | Self::SlogpVsa4
            | Self::SlogpVsa5
            | Self::SlogpVsa6
            | Self::SlogpVsa7
            | Self::SlogpVsa8
            | Self::SlogpVsa9
            | Self::SlogpVsa10
            | Self::SlogpVsa11
            | Self::SlogpVsa12
            | Self::SmrVsa1
            | Self::SmrVsa2
            | Self::SmrVsa3
            | Self::SmrVsa4
            | Self::SmrVsa5
            | Self::SmrVsa6
            | Self::SmrVsa7
            | Self::SmrVsa8
            | Self::SmrVsa9
            | Self::SmrVsa10
            | Self::Qed => Category::Descriptors,
            Self::SmilesRead => Category::Notation,
            Self::Sanitize
            | Self::Kekulize
            | Self::AddHydrogens
            | Self::RemoveHydrogens
            | Self::Valence => Category::Chemistry,
            Self::MolecularWeight | Self::ExactMolecularWeight | Self::MolecularFormula => {
                Category::Descriptors
            }
            Self::NumHeavyAtoms => Category::Descriptors,
            Self::TotalAtomCount => Category::Descriptors,
            Self::LipinskiHBA | Self::LipinskiHBD | Self::FractionCSP3 => Category::Descriptors,
            Self::NumHeteroatoms => Category::Descriptors,
            Self::NumHba => Category::Descriptors,
            Self::NumHbd => Category::Descriptors,
            Self::NumRings
            | Self::NumHeterocycles
            | Self::NumAromaticRings
            | Self::NumSaturatedRings
            | Self::NumAliphaticRings
            | Self::NumAromaticHeterocycles
            | Self::NumAromaticCarbocycles
            | Self::NumAliphaticHeterocycles
            | Self::NumAliphaticCarbocycles
            | Self::NumSaturatedHeterocycles
            | Self::NumSaturatedCarbocycles => Category::Descriptors,
            Self::CipLabels | Self::PotentialStereo => Category::Stereo,
            Self::Coordinates2d | Self::Svg => Category::Depiction,
            Self::DistanceMatrix => Category::Chemistry,
            Self::MorganFingerprint
            | Self::MorganSparseFingerprint
            | Self::MorganCountFingerprint
            | Self::MorganSparseCountFingerprint => Category::FingerprintGeneration,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum TaskId {
    TautomerEnumeration,
    TautomerCanonicalization,
    NumAmideBonds,
    NumSpiroAtoms,
    NumBridgeheadAtoms,
    NumAtomStereoCenters,
    NumUnspecifiedAtomStereoCenters,
    NumRotatableBonds,
    CrippenDescriptors,
    LabuteAsa,
    LabuteAsaContributions,
    Tpsa,
    SlogpVsa,
    SmrVsa,
    SlogpVsa1,
    SlogpVsa2,
    SlogpVsa3,
    SlogpVsa4,
    SlogpVsa5,
    SlogpVsa6,
    SlogpVsa7,
    SlogpVsa8,
    SlogpVsa9,
    SlogpVsa10,
    SlogpVsa11,
    SlogpVsa12,
    SmrVsa1,
    SmrVsa2,
    SmrVsa3,
    SmrVsa4,
    SmrVsa5,
    SmrVsa6,
    SmrVsa7,
    SmrVsa8,
    SmrVsa9,
    SmrVsa10,
    Qed,

    DistanceMatrix,
    SmilesRead,
    Sanitize,
    Kekulize,
    MolecularWeight,
    ExactMolecularWeight,
    MolecularFormula,
    AddHydrogens,
    RemoveHydrogens,
    CipLabels,
    PotentialStereo,
    Coordinates2d,
    Svg,
    Valence,
    NumHeavyAtoms,
    TotalAtomCount,
    LipinskiHBA,
    LipinskiHBD,
    FractionCSP3,
    NumHeteroatoms,
    NumHba,
    NumHbd,
    NumRings,
    NumHeterocycles,
    NumAromaticRings,
    NumSaturatedRings,
    NumAliphaticRings,
    NumAromaticHeterocycles,
    NumAromaticCarbocycles,
    NumAliphaticHeterocycles,
    NumAliphaticCarbocycles,
    NumSaturatedHeterocycles,
    NumSaturatedCarbocycles,
    MorganFingerprint,
    MorganSparseFingerprint,
    MorganCountFingerprint,
    MorganSparseCountFingerprint,
    Chi0,
    Chi1,
    HallKierAlpha,
    HallKierAlphaWithContributions,
    Kappa1,
    Kappa2,
    Kappa3,
    Phi,
    Mqns,
    Chi0V,
    Chi1V,
    Chi2V,
    Chi3V,
    Chi4V,
    Chi0N,
    Chi1N,
    Chi2N,
    Chi3N,
    Chi4N,
    ChiNV,
    ChiNN,
}

/// `None` and an explicitly empty selection must never be conflated.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum CipSelection {
    All,
    FirstTaggedAtom,
    FirstStereoBond,
    Empty,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum MorganOutputKind {
    DenseBits,
    SparseBits,
    HashedCounts,
    SparseCounts,
}

impl MorganOutputKind {
    pub const fn task_name(self) -> &'static str {
        match self {
            Self::DenseBits => "fingerprint_morgan",
            Self::SparseBits => "fingerprint_morgan_sparse",
            Self::HashedCounts => "fingerprint_morgan_count",
            Self::SparseCounts => "fingerprint_morgan_sparse_count",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum MorganInvariantKind {
    Connectivity,
    Features,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum TautomerCatalog {
    Current,
    V1,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TautomerProfile {
    pub catalog: TautomerCatalog,
    pub max_tautomers: u32,
    pub max_transforms: u32,
    pub remove_sp3_stereo: bool,
    pub remove_bond_stereo: bool,
    pub remove_isotopic_hydrogens: bool,
    pub reassign_stereo: bool,
}
impl TautomerProfile {
    pub const fn branch(self) -> &'static str {
        match self.catalog {
            TautomerCatalog::Current => "default",
            TautomerCatalog::V1 => "v1",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum RotatableBondMode {
    Default,
    NonStrict,
    Strict,
    StrictLinkages,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum VsaBins {
    Default,
    CustomDuplicates,
}
impl VsaBins {
    pub const fn values(self) -> Option<&'static [f64]> {
        match self {
            Self::Default => None,
            Self::CustomDuplicates => Some(&[-0.2, 0.0, 0.25, 0.25, 0.8]),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Profile {
    TautomerEnumeration {
        parameters: TautomerProfile,
    },
    TautomerCanonicalization {
        parameters: TautomerProfile,
    },
    NumAmideBonds,
    NumSpiroAtoms,
    NumBridgeheadAtoms,
    NumAtomStereoCenters,
    NumUnspecifiedAtomStereoCenters,
    NumRotatableBonds {
        mode: Option<RotatableBondMode>,
    },
    CrippenDescriptors {
        include_hydrogens: Option<bool>,
        force: bool,
    },
    LabuteAsa {
        include_hydrogens: Option<bool>,
        force: bool,
    },
    LabuteAsaContributions {
        include_hydrogens: Option<bool>,
        force: bool,
    },
    Tpsa {
        include_sulfur_phosphorus: Option<bool>,
        force: bool,
    },
    SlogpVsa {
        bins: VsaBins,
        force: Option<bool>,
    },
    SmrVsa {
        bins: VsaBins,
        force: Option<bool>,
    },
    SlogpVsa1,
    SlogpVsa2,
    SlogpVsa3,
    SlogpVsa4,
    SlogpVsa5,
    SlogpVsa6,
    SlogpVsa7,
    SlogpVsa8,
    SlogpVsa9,
    SlogpVsa10,
    SlogpVsa11,
    SlogpVsa12,
    SmrVsa1,
    SmrVsa2,
    SmrVsa3,
    SmrVsa4,
    SmrVsa5,
    SmrVsa6,
    SmrVsa7,
    SmrVsa8,
    SmrVsa9,
    SmrVsa10,
    Qed,
    Chi0VWithParams {
        force: bool,
    },
    Chi1VWithParams {
        force: bool,
    },
    Chi2VWithParams {
        force: bool,
    },
    Chi3VWithParams {
        force: bool,
    },
    Chi4VWithParams {
        force: bool,
    },
    Chi0NWithParams {
        force: bool,
    },
    Chi1NWithParams {
        force: bool,
    },
    Chi2NWithParams {
        force: bool,
    },
    Chi3NWithParams {
        force: bool,
    },
    Chi4NWithParams {
        force: bool,
    },
    ChiNVWithParams {
        order: u32,
        force: bool,
    },
    ChiNNWithParams {
        order: u32,
        force: bool,
    },
    LabuteAsaCacheSequence,
    SlogpVsaCacheSequence,
    SmrVsaCacheSequence,
    ChiNVCacheSequence,
    ChiNNCacheSequence,

    DistanceMatrix {
        use_bond_order: bool,
        use_atom_weights: bool,
    },
    SmilesRead {
        sanitize: bool,
        remove_hs: bool,
    },
    SanitizeAll,
    Kekulize {
        clear_aromatic_flags: bool,
    },
    MolecularWeight {
        only_heavy: bool,
    },
    ExactMolecularWeight {
        only_heavy: bool,
    },
    MolecularFormula {
        separate_isotopes: bool,
        abbreviate_h_isotopes: bool,
    },
    AddHydrogens {
        explicit_only: bool,
    },
    RemoveHydrogens {
        sanitize: bool,
    },
    CipLabels {
        selection: CipSelection,
        max_recursive_iterations: u32,
    },
    PotentialStereo {
        clean: bool,
        flag_possible: bool,
        allow_nontetrahedral: bool,
    },
    Chi0,
    Chi1,
    HallKierAlpha,
    HallKierAlphaWithContributions,
    Kappa1,
    Kappa2,
    Kappa3,
    Phi,
    Mqns {
        force: bool,
    },
    Chi0V,
    Chi1V,
    Chi2V,
    Chi3V,
    Chi4V,
    Chi0N,
    Chi1N,
    Chi2N,
    Chi3N,
    Chi4N,
    ChiNV {
        order: u32,
    },
    ChiNN {
        order: u32,
    },
    Coordinates2dDefault,
    /// The only drawing profile: 300x300, default source preparation.
    SvgDefault,
    Valence {
        strict: bool,
    },
    NumHeavyAtoms {
        remove_hs: bool,
    },
    TotalAtomCount {
        remove_hs: bool,
    },
    LipinskiHBA {
        remove_hs: bool,
    },
    LipinskiHBD {
        remove_hs: bool,
    },
    FractionCSP3 {
        remove_hs: bool,
    },
    NumHeteroatoms {
        remove_hs: bool,
    },
    NumHba {
        remove_hs: bool,
    },
    NumHbd {
        remove_hs: bool,
    },
    NumRings {
        remove_hs: bool,
    },
    NumHeterocycles {
        remove_hs: bool,
    },
    NumAromaticRings {
        remove_hs: bool,
    },
    NumSaturatedRings {
        remove_hs: bool,
    },
    NumAliphaticRings {
        remove_hs: bool,
    },
    NumAromaticHeterocycles {
        remove_hs: bool,
    },
    NumAromaticCarbocycles {
        remove_hs: bool,
    },
    NumAliphaticHeterocycles {
        remove_hs: bool,
    },
    NumAliphaticCarbocycles {
        remove_hs: bool,
    },
    NumSaturatedHeterocycles {
        remove_hs: bool,
    },
    NumSaturatedCarbocycles {
        remove_hs: bool,
    },
    Morgan {
        output: MorganOutputKind,
        radius: u32,
        include_chirality: bool,
        invariants: MorganInvariantKind,
        count_simulation: bool,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputState {
    SmilesText,
    UnsanitizedHydrogensRetained,
    SanitizedHydrogensRemoved,
    SanitizedThenAddAllHydrogens,
    SanitizedHydrogensPerProfile,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Comparison {
    TautomerFullEnumeration,
    TautomerFullCanonicalization,
    MatrixBits,
    TopologyAndOutcome,
    Float64Bits,
    Float64PairBits,
    Float64VectorBits,
    LabuteAsaContributionsBits,
    Float64ContributionsBits,
    UnsignedVector,
    ExactText,
    SvgText,
    CipLabelsAndOutcome,
    StereoInfoAndCleanedTopology,
    CoordinatesAndTopology,
    ValenceRowsAndOutcome,
    Unsigned,
    MorganFingerprintAndAdditionalOutput,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Prerequisite {
    MolecularPipeline,
    PublicValenceReadoutAndMolecularPipeline,
}

#[derive(Debug, Clone, Copy)]
pub struct Task {
    pub id: TaskId,
    pub input: InputState,
    pub comparison: Comparison,
    pub prerequisite: Prerequisite,
}

use Comparison::*;
use InputState::*;
use Prerequisite::*;
use TaskId::*;

/// Frozen matrices. Executable registration lives in the parent module;
/// unregistered matrices remain planned, not passing evidence.
pub const TASKS: &[Task] = &[
    Task {
        id: DistanceMatrix,
        input: SanitizedHydrogensRemoved,
        comparison: MatrixBits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SmilesRead,
        input: SmilesText,
        comparison: TopologyAndOutcome,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Sanitize,
        input: UnsanitizedHydrogensRetained,
        comparison: TopologyAndOutcome,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Kekulize,
        input: SanitizedHydrogensRemoved,
        comparison: TopologyAndOutcome,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: MolecularWeight,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: ExactMolecularWeight,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: MolecularFormula,
        input: SanitizedHydrogensRemoved,
        comparison: ExactText,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumHeavyAtoms,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: TotalAtomCount,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: LipinskiHBA,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: LipinskiHBD,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: FractionCSP3,
        input: SanitizedHydrogensPerProfile,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumHeteroatoms,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumHba,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumHbd,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumRings,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumHeterocycles,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumAromaticRings,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumSaturatedRings,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumAliphaticRings,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumAromaticHeterocycles,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumAromaticCarbocycles,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumAliphaticHeterocycles,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumAliphaticCarbocycles,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumSaturatedHeterocycles,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumSaturatedCarbocycles,
        input: SanitizedHydrogensPerProfile,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: AddHydrogens,
        input: SanitizedHydrogensRemoved,
        comparison: TopologyAndOutcome,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: RemoveHydrogens,
        input: SanitizedThenAddAllHydrogens,
        comparison: TopologyAndOutcome,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: CipLabels,
        input: SanitizedHydrogensRemoved,
        comparison: CipLabelsAndOutcome,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: PotentialStereo,
        input: SanitizedHydrogensRemoved,
        comparison: StereoInfoAndCleanedTopology,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Coordinates2d,
        input: SanitizedHydrogensRemoved,
        comparison: CoordinatesAndTopology,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Svg,
        input: SanitizedHydrogensRemoved,
        comparison: SvgText,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Valence,
        input: UnsanitizedHydrogensRetained,
        comparison: ValenceRowsAndOutcome,
        prerequisite: PublicValenceReadoutAndMolecularPipeline,
    },
    Task {
        id: MorganFingerprint,
        input: SanitizedHydrogensRemoved,
        comparison: MorganFingerprintAndAdditionalOutput,
        prerequisite: PublicValenceReadoutAndMolecularPipeline,
    },
    Task {
        id: MorganSparseFingerprint,
        input: SanitizedHydrogensRemoved,
        comparison: MorganFingerprintAndAdditionalOutput,
        prerequisite: PublicValenceReadoutAndMolecularPipeline,
    },
    Task {
        id: MorganCountFingerprint,
        input: SanitizedHydrogensRemoved,
        comparison: MorganFingerprintAndAdditionalOutput,
        prerequisite: PublicValenceReadoutAndMolecularPipeline,
    },
    Task {
        id: MorganSparseCountFingerprint,
        input: SanitizedHydrogensRemoved,
        comparison: MorganFingerprintAndAdditionalOutput,
        prerequisite: PublicValenceReadoutAndMolecularPipeline,
    },
    Task {
        id: Chi0,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Chi1,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: HallKierAlpha,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: HallKierAlphaWithContributions,
        input: SanitizedHydrogensRemoved,
        comparison: Float64ContributionsBits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Kappa1,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Kappa2,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Kappa3,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Phi,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Mqns,
        input: SanitizedHydrogensRemoved,
        comparison: UnsignedVector,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Chi0V,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Chi1V,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Chi2V,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Chi3V,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Chi4V,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Chi0N,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Chi1N,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Chi2N,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Chi3N,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Chi4N,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: ChiNV,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: ChiNN,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: TautomerEnumeration,
        input: SanitizedHydrogensRemoved,
        comparison: TautomerFullEnumeration,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: TautomerCanonicalization,
        input: SanitizedHydrogensRemoved,
        comparison: TautomerFullCanonicalization,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumAmideBonds,
        input: SanitizedHydrogensRemoved,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumSpiroAtoms,
        input: SanitizedHydrogensRemoved,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumBridgeheadAtoms,
        input: SanitizedHydrogensRemoved,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumAtomStereoCenters,
        input: SanitizedHydrogensRemoved,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumUnspecifiedAtomStereoCenters,
        input: SanitizedHydrogensRemoved,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: NumRotatableBonds,
        input: SanitizedHydrogensRemoved,
        comparison: Unsigned,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: CrippenDescriptors,
        input: SanitizedHydrogensRemoved,
        comparison: Float64PairBits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: LabuteAsa,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: LabuteAsaContributions,
        input: SanitizedHydrogensRemoved,
        comparison: LabuteAsaContributionsBits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Tpsa,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SlogpVsa,
        input: SanitizedHydrogensRemoved,
        comparison: Float64VectorBits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SmrVsa,
        input: SanitizedHydrogensRemoved,
        comparison: Float64VectorBits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SlogpVsa1,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SlogpVsa2,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SlogpVsa3,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SlogpVsa4,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SlogpVsa5,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SlogpVsa6,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SlogpVsa7,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SlogpVsa8,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SlogpVsa9,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SlogpVsa10,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SlogpVsa11,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SlogpVsa12,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SmrVsa1,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SmrVsa2,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SmrVsa3,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SmrVsa4,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SmrVsa5,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SmrVsa6,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SmrVsa7,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SmrVsa8,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SmrVsa9,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: SmrVsa10,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
    Task {
        id: Qed,
        input: SanitizedHydrogensRemoved,
        comparison: Float64Bits,
        prerequisite: MolecularPipeline,
    },
];

impl TaskId {
    /// All combinations below are deliberate; there is no implicit sampling.
    /// Unlisted API options are pinned in TASKS.md, not silently expanded.
    pub fn profiles(self) -> Vec<Profile> {
        let booleans = [false, true];
        match self {
            TautomerEnumeration | TautomerCanonicalization => {
                [TautomerCatalog::Current, TautomerCatalog::V1]
                    .map(|catalog| {
                        let parameters = TautomerProfile {
                            catalog,
                            max_tautomers: 1000,
                            max_transforms: 1000,
                            remove_sp3_stereo: true,
                            remove_bond_stereo: true,
                            remove_isotopic_hydrogens: true,
                            reassign_stereo: true,
                        };
                        if self == TautomerEnumeration {
                            Profile::TautomerEnumeration { parameters }
                        } else {
                            Profile::TautomerCanonicalization { parameters }
                        }
                    })
                    .into()
            }
            NumAmideBonds => vec![Profile::NumAmideBonds],
            NumSpiroAtoms => vec![Profile::NumSpiroAtoms],
            NumBridgeheadAtoms => vec![Profile::NumBridgeheadAtoms],
            NumAtomStereoCenters => vec![Profile::NumAtomStereoCenters],
            NumUnspecifiedAtomStereoCenters => vec![Profile::NumUnspecifiedAtomStereoCenters],
            NumRotatableBonds => vec![
                Profile::NumRotatableBonds { mode: None },
                Profile::NumRotatableBonds {
                    mode: Some(RotatableBondMode::Default),
                },
                Profile::NumRotatableBonds {
                    mode: Some(RotatableBondMode::NonStrict),
                },
                Profile::NumRotatableBonds {
                    mode: Some(RotatableBondMode::Strict),
                },
                Profile::NumRotatableBonds {
                    mode: Some(RotatableBondMode::StrictLinkages),
                },
            ],
            CrippenDescriptors => [Profile::CrippenDescriptors {
                include_hydrogens: None,
                force: false,
            }]
            .into_iter()
            .chain(booleans.into_iter().flat_map(|flag| {
                booleans.map(move |force| Profile::CrippenDescriptors {
                    include_hydrogens: Some(flag),
                    force,
                })
            }))
            .collect(),
            LabuteAsa => [Profile::LabuteAsa {
                include_hydrogens: None,
                force: false,
            }]
            .into_iter()
            .chain(booleans.into_iter().flat_map(|flag| {
                booleans.map(move |force| Profile::LabuteAsa {
                    include_hydrogens: Some(flag),
                    force,
                })
            }))
            .chain([Profile::LabuteAsaCacheSequence])
            .collect(),
            LabuteAsaContributions => [Profile::LabuteAsaContributions {
                include_hydrogens: None,
                force: false,
            }]
            .into_iter()
            .chain(booleans.into_iter().flat_map(|flag| {
                booleans.map(move |force| Profile::LabuteAsaContributions {
                    include_hydrogens: Some(flag),
                    force,
                })
            }))
            .collect(),
            Tpsa => [Profile::Tpsa {
                include_sulfur_phosphorus: None,
                force: false,
            }]
            .into_iter()
            .chain(booleans.into_iter().flat_map(|flag| {
                booleans.map(move |force| Profile::Tpsa {
                    include_sulfur_phosphorus: Some(flag),
                    force,
                })
            }))
            .collect(),
            SlogpVsa => [Profile::SlogpVsa {
                bins: VsaBins::Default,
                force: None,
            }]
            .into_iter()
            .chain(
                [VsaBins::Default, VsaBins::CustomDuplicates]
                    .into_iter()
                    .flat_map(|bins| {
                        booleans.map(move |force| Profile::SlogpVsa {
                            bins,
                            force: Some(force),
                        })
                    }),
            )
            .chain([Profile::SlogpVsaCacheSequence])
            .collect(),
            SmrVsa => [Profile::SmrVsa {
                bins: VsaBins::Default,
                force: None,
            }]
            .into_iter()
            .chain(
                [VsaBins::Default, VsaBins::CustomDuplicates]
                    .into_iter()
                    .flat_map(|bins| {
                        booleans.map(move |force| Profile::SmrVsa {
                            bins,
                            force: Some(force),
                        })
                    }),
            )
            .chain([Profile::SmrVsaCacheSequence])
            .collect(),
            SlogpVsa1 => vec![Profile::SlogpVsa1],
            SlogpVsa2 => vec![Profile::SlogpVsa2],
            SlogpVsa3 => vec![Profile::SlogpVsa3],
            SlogpVsa4 => vec![Profile::SlogpVsa4],
            SlogpVsa5 => vec![Profile::SlogpVsa5],
            SlogpVsa6 => vec![Profile::SlogpVsa6],
            SlogpVsa7 => vec![Profile::SlogpVsa7],
            SlogpVsa8 => vec![Profile::SlogpVsa8],
            SlogpVsa9 => vec![Profile::SlogpVsa9],
            SlogpVsa10 => vec![Profile::SlogpVsa10],
            SlogpVsa11 => vec![Profile::SlogpVsa11],
            SlogpVsa12 => vec![Profile::SlogpVsa12],
            SmrVsa1 => vec![Profile::SmrVsa1],
            SmrVsa2 => vec![Profile::SmrVsa2],
            SmrVsa3 => vec![Profile::SmrVsa3],
            SmrVsa4 => vec![Profile::SmrVsa4],
            SmrVsa5 => vec![Profile::SmrVsa5],
            SmrVsa6 => vec![Profile::SmrVsa6],
            SmrVsa7 => vec![Profile::SmrVsa7],
            SmrVsa8 => vec![Profile::SmrVsa8],
            SmrVsa9 => vec![Profile::SmrVsa9],
            SmrVsa10 => vec![Profile::SmrVsa10],
            Qed => vec![Profile::Qed],
            Chi0 => vec![Profile::Chi0],
            Chi1 => vec![Profile::Chi1],
            HallKierAlpha => vec![Profile::HallKierAlpha],
            HallKierAlphaWithContributions => vec![Profile::HallKierAlphaWithContributions],
            Kappa1 => vec![Profile::Kappa1],
            Kappa2 => vec![Profile::Kappa2],
            Kappa3 => vec![Profile::Kappa3],
            Phi => vec![Profile::Phi],
            Mqns => booleans.map(|force| Profile::Mqns { force }).into(),
            Chi0V => vec![
                Profile::Chi0V,
                Profile::Chi0VWithParams { force: false },
                Profile::Chi0VWithParams { force: true },
            ],
            Chi1V => vec![
                Profile::Chi1V,
                Profile::Chi1VWithParams { force: false },
                Profile::Chi1VWithParams { force: true },
            ],
            Chi2V => vec![
                Profile::Chi2V,
                Profile::Chi2VWithParams { force: false },
                Profile::Chi2VWithParams { force: true },
            ],
            Chi3V => vec![
                Profile::Chi3V,
                Profile::Chi3VWithParams { force: false },
                Profile::Chi3VWithParams { force: true },
            ],
            Chi4V => vec![
                Profile::Chi4V,
                Profile::Chi4VWithParams { force: false },
                Profile::Chi4VWithParams { force: true },
            ],
            Chi0N => vec![
                Profile::Chi0N,
                Profile::Chi0NWithParams { force: false },
                Profile::Chi0NWithParams { force: true },
            ],
            Chi1N => vec![
                Profile::Chi1N,
                Profile::Chi1NWithParams { force: false },
                Profile::Chi1NWithParams { force: true },
            ],
            Chi2N => vec![
                Profile::Chi2N,
                Profile::Chi2NWithParams { force: false },
                Profile::Chi2NWithParams { force: true },
            ],
            Chi3N => vec![
                Profile::Chi3N,
                Profile::Chi3NWithParams { force: false },
                Profile::Chi3NWithParams { force: true },
            ],
            Chi4N => vec![
                Profile::Chi4N,
                Profile::Chi4NWithParams { force: false },
                Profile::Chi4NWithParams { force: true },
            ],
            ChiNV => (0..=6)
                .map(|order| Profile::ChiNV { order })
                .chain((0..=6).flat_map(|order| {
                    booleans.map(move |force| Profile::ChiNVWithParams { order, force })
                }))
                .chain([Profile::ChiNVCacheSequence])
                .collect(),
            ChiNN => (0..=6)
                .map(|order| Profile::ChiNN { order })
                .chain((0..=6).flat_map(|order| {
                    booleans.map(move |force| Profile::ChiNNWithParams { order, force })
                }))
                .chain([Profile::ChiNNCacheSequence])
                .collect(),
            DistanceMatrix => booleans
                .into_iter()
                .flat_map(|use_bond_order| {
                    booleans.map(move |use_atom_weights| Profile::DistanceMatrix {
                        use_bond_order,
                        use_atom_weights,
                    })
                })
                .collect(),
            SmilesRead => booleans
                .into_iter()
                .flat_map(|sanitize| {
                    booleans.map(move |remove_hs| Profile::SmilesRead {
                        sanitize,
                        remove_hs,
                    })
                })
                .collect(),
            Sanitize => vec![Profile::SanitizeAll],
            Kekulize => booleans
                .map(|clear_aromatic_flags| Profile::Kekulize {
                    clear_aromatic_flags,
                })
                .into(),
            MolecularWeight => booleans
                .map(|only_heavy| Profile::MolecularWeight { only_heavy })
                .into(),
            ExactMolecularWeight => booleans
                .map(|only_heavy| Profile::ExactMolecularWeight { only_heavy })
                .into(),
            MolecularFormula => booleans
                .into_iter()
                .flat_map(|separate_isotopes| {
                    booleans.map(move |abbreviate_h_isotopes| Profile::MolecularFormula {
                        separate_isotopes,
                        abbreviate_h_isotopes,
                    })
                })
                .collect(),
            AddHydrogens => booleans
                .map(|explicit_only| Profile::AddHydrogens { explicit_only })
                .into(),
            RemoveHydrogens => booleans
                .map(|sanitize| Profile::RemoveHydrogens { sanitize })
                .into(),
            CipLabels => vec![
                Profile::CipLabels {
                    selection: CipSelection::All,
                    max_recursive_iterations: 0,
                },
                Profile::CipLabels {
                    selection: CipSelection::FirstTaggedAtom,
                    max_recursive_iterations: 0,
                },
                Profile::CipLabels {
                    selection: CipSelection::FirstStereoBond,
                    max_recursive_iterations: 0,
                },
                Profile::CipLabels {
                    selection: CipSelection::Empty,
                    max_recursive_iterations: 1,
                },
            ],
            PotentialStereo => booleans
                .into_iter()
                .flat_map(|clean| {
                    booleans.into_iter().flat_map(move |flag_possible| {
                        booleans.map(move |allow_nontetrahedral| Profile::PotentialStereo {
                            clean,
                            flag_possible,
                            allow_nontetrahedral,
                        })
                    })
                })
                .collect(),
            Coordinates2d => vec![Profile::Coordinates2dDefault],
            Svg => vec![Profile::SvgDefault],
            Valence => booleans.map(|strict| Profile::Valence { strict }).into(),
            NumHeavyAtoms => booleans
                .map(|remove_hs| Profile::NumHeavyAtoms { remove_hs })
                .into(),
            TotalAtomCount => booleans
                .map(|remove_hs| Profile::TotalAtomCount { remove_hs })
                .into(),
            LipinskiHBA => booleans
                .map(|remove_hs| Profile::LipinskiHBA { remove_hs })
                .into(),
            LipinskiHBD => booleans
                .map(|remove_hs| Profile::LipinskiHBD { remove_hs })
                .into(),
            FractionCSP3 => booleans
                .map(|remove_hs| Profile::FractionCSP3 { remove_hs })
                .into(),
            NumHeteroatoms => booleans
                .map(|remove_hs| Profile::NumHeteroatoms { remove_hs })
                .into(),
            NumHba => booleans
                .map(|remove_hs| Profile::NumHba { remove_hs })
                .into(),
            NumHbd => booleans
                .map(|remove_hs| Profile::NumHbd { remove_hs })
                .into(),
            NumRings => booleans
                .map(|remove_hs| Profile::NumRings { remove_hs })
                .into(),
            NumHeterocycles => booleans
                .map(|remove_hs| Profile::NumHeterocycles { remove_hs })
                .into(),
            NumAromaticRings => booleans
                .map(|remove_hs| Profile::NumAromaticRings { remove_hs })
                .into(),
            NumSaturatedRings => booleans
                .map(|remove_hs| Profile::NumSaturatedRings { remove_hs })
                .into(),
            NumAliphaticRings => booleans
                .map(|remove_hs| Profile::NumAliphaticRings { remove_hs })
                .into(),
            NumAromaticHeterocycles => booleans
                .map(|remove_hs| Profile::NumAromaticHeterocycles { remove_hs })
                .into(),
            NumAromaticCarbocycles => booleans
                .map(|remove_hs| Profile::NumAromaticCarbocycles { remove_hs })
                .into(),
            NumAliphaticHeterocycles => booleans
                .map(|remove_hs| Profile::NumAliphaticHeterocycles { remove_hs })
                .into(),
            NumAliphaticCarbocycles => booleans
                .map(|remove_hs| Profile::NumAliphaticCarbocycles { remove_hs })
                .into(),
            NumSaturatedHeterocycles => booleans
                .map(|remove_hs| Profile::NumSaturatedHeterocycles { remove_hs })
                .into(),
            NumSaturatedCarbocycles => booleans
                .map(|remove_hs| Profile::NumSaturatedCarbocycles { remove_hs })
                .into(),
            MorganFingerprint
            | MorganSparseFingerprint
            | MorganCountFingerprint
            | MorganSparseCountFingerprint => {
                let output = match self {
                    MorganFingerprint => MorganOutputKind::DenseBits,
                    MorganSparseFingerprint => MorganOutputKind::SparseBits,
                    MorganCountFingerprint => MorganOutputKind::HashedCounts,
                    MorganSparseCountFingerprint => MorganOutputKind::SparseCounts,
                    _ => unreachable!("the enclosing task id is a Morgan output family"),
                };
                [2_u32, 3]
                    .into_iter()
                    .flat_map(|radius| {
                        [false, true]
                            .into_iter()
                            .flat_map(move |include_chirality| {
                                [
                                    MorganInvariantKind::Connectivity,
                                    MorganInvariantKind::Features,
                                ]
                                .into_iter()
                                .flat_map(move |invariants| {
                                    [false, true].map(move |count_simulation| Profile::Morgan {
                                        output,
                                        radius,
                                        include_chirality,
                                        invariants,
                                        count_simulation,
                                    })
                                })
                            })
                    })
                    .collect()
            }
        }
    }
}
