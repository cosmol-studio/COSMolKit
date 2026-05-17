// RDKit marker convention defined in dev/source_reproduction_protocol.md.

use crate::{
    AtomId, AtomQueryPredicate, Bond, BondDirection, BondId, BondOrder, BondStereo, ChiralTag,
    Molecule, QueryNode, ValenceError,
};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::atomic::{AtomicU64, Ordering};

thread_local! {
    static RANDOM_SMILES_SEED: std::cell::Cell<u64> = const { std::cell::Cell::new(0) };
}

static RANDOM_SMILES_COUNTER: AtomicU64 = AtomicU64::new(0x9e37_79b9_7f4a_7c15);
const CANON_MAX_NATOMS: i64 = 5000;
const CANON_MAX_BONDTYPE: i64 = 32;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmilesWriteParams {
    pub do_isomeric_smiles: bool,
    pub do_kekule: bool,
    pub canonical: bool,
    pub clean_stereo: bool,
    pub all_bonds_explicit: bool,
    pub all_hydrogens_explicit: bool,
    pub do_random: bool,
    pub rooted_at_atom: Option<usize>,
    pub include_dative_bonds: bool,
    pub ignore_atom_map_numbers: bool,
}

impl Default for SmilesWriteParams {
    fn default() -> Self {
        Self {
            do_isomeric_smiles: true,
            do_kekule: false,
            canonical: true,
            clean_stereo: true,
            all_bonds_explicit: false,
            all_hydrogens_explicit: false,
            do_random: false,
            rooted_at_atom: None,
            include_dative_bonds: true,
            ignore_atom_map_numbers: false,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CxSmilesFields(u32);

impl CxSmilesFields {
    pub const NONE: Self = Self(0);
    pub const ATOM_LABELS: Self = Self(1 << 0);
    pub const MOLFILE_VALUES: Self = Self(1 << 1);
    pub const COORDS: Self = Self(1 << 2);
    pub const RADICALS: Self = Self(1 << 3);
    pub const ATOM_PROPS: Self = Self(1 << 4);
    pub const LINKNODES: Self = Self(1 << 5);
    pub const ENHANCED_STEREO: Self = Self(1 << 6);
    pub const SGROUPS: Self = Self(1 << 7);
    pub const POLYMER: Self = Self(1 << 8);
    pub const BOND_CFG: Self = Self(1 << 9);
    pub const BOND_ATROPISOMER: Self = Self(1 << 10);
    pub const COORDINATE_BONDS: Self = Self(1 << 11);
    pub const HYDROGEN_BONDS: Self = Self(1 << 12);
    pub const ZERO_BONDS: Self = Self(1 << 13);
    pub const ALL: Self = Self(0x7fff_ffff);
    pub const ALL_BUT_COORDS: Self = Self(Self::ALL.0 ^ Self::COORDS.0);

    #[must_use]
    pub const fn bits(self) -> u32 {
        self.0
    }

    #[must_use]
    pub const fn contains(self, other: Self) -> bool {
        self.0 & other.0 == other.0
    }

    #[must_use]
    pub const fn combine(self, other: Self) -> Self {
        Self(self.0 | other.0)
    }
}

impl std::ops::BitOr for CxSmilesFields {
    type Output = Self;

    fn bitor(self, rhs: Self) -> Self {
        Self(self.0 | rhs.0)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RestoreBondDirOption {
    None,
    True,
    Clear,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SmilesOutputMode {
    PlainSmiles,
    CxSmiles {
        fields: CxSmilesFields,
        restore_bond_dirs: RestoreBondDirOption,
        include_stereo_groups: bool,
    },
}

/// Stages of the SMILES writing pipeline used for internal guard diagnostics.
///
/// Each variant corresponds to a writer phase used by invariant-violation
/// errors to identify where an internal contract was broken.
///
/// ## Status of each stage
///
/// - `ShortTermAtomWriter`: Atom-level guards / deferred operations.
///   Includes: empty rank edge case, chiral/query/radical atoms in
///   the minimal-fast-path guard, and non-whitelisted atom properties.
/// - `ShortTermBondWriter`: Bond-level guards / deferred operations.
///   Includes: dative-bond stripping, Unknown/EitherDouble/Any bond
///   direction/stereo for plain SMILES, ring-closure digit exhaustion,
///   and non-standard bond orders in the fast-path guard.
/// - `LongTermCanonicalRanking`: Canonical rank calculation. Defined
///   but unused in deferred-error path — canonical ranking errors use
///   SmilesWriteError::CanonicalRank directly.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SmilesPlanStage {
    ShortTermAtomWriter,
    ShortTermBondWriter,
    LongTermCanonicalRanking,
}

impl SmilesPlanStage {
    const fn as_str(self) -> &'static str {
        match self {
            Self::ShortTermAtomWriter => "ShortTermAtomWriter",
            Self::ShortTermBondWriter => "ShortTermBondWriter",
            Self::LongTermCanonicalRanking => "LongTermCanonicalRanking",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum SmilesWriteError {
    #[error(transparent)]
    UnsupportedFeature(#[from] crate::UnsupportedFeatureError),
    #[error("canonical ranking failed: {source}")]
    CanonicalRank { source: crate::KekulizeError },
    #[error("kekulization failed: {source}")]
    Kekulize {
        #[from]
        source: crate::KekulizeError,
    },
    #[error("operation failed while preparing SMILES output: {source}")]
    Operation {
        #[from]
        source: crate::OperationError,
    },
    #[error("valence calculation failed: {source}")]
    Valence {
        #[from]
        source: ValenceError,
    },
    #[error("stereochemistry preparation failed: {source}")]
    Stereo {
        #[from]
        source: crate::StereoError,
    },
    #[error("ring finding failed while preparing SMILES output: {source}")]
    RingFinding {
        #[from]
        source: crate::RingFindingError,
    },
    #[error("atom index {atom} is out of range")]
    AtomOutOfRange { atom: usize },
    #[error("bond index {bond} is out of range")]
    BondOutOfRange { bond: usize },
    #[error("rooted atom index {atom} is out of range")]
    RootedAtomOutOfRange { atom: usize },
    #[error("rooted atom index {atom} is not present in atoms_to_use")]
    RootedAtomNotInFragment { atom: usize },
    #[error(
        "rooted atom index {atom} requires a single-fragment molecule when bonds_to_use is omitted"
    )]
    RootedAtomRequiresSingleFragment { atom: usize },
    #[error("atom symbol override vector has length {len}, expected at least {expected}")]
    AtomSymbolsTooShort { len: usize, expected: usize },
    #[error("bond symbol override vector has length {len}, expected at least {expected}")]
    BondSymbolsTooShort { len: usize, expected: usize },
    #[error(
        "invalid non-tetrahedral chiral permutation {permutation} for {chiral_tag:?}; max allowed is {limit}"
    )]
    InvalidChiralPermutation {
        chiral_tag: ChiralTag,
        permutation: u32,
        limit: u32,
    },
    #[error("invalid ring stereochemistry state on atom {atom}: {requirement}")]
    InvalidRingStereoState {
        atom: usize,
        requirement: &'static str,
    },
    #[error("internal SMILES writer invariant violated in {stage}: {message}")]
    InvariantViolation {
        stage: &'static str,
        message: &'static str,
    },
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
struct SmilesWriteContext {
    atom_output_order: Vec<AtomId>,
    bond_output_order: Vec<BondId>,
    ring_closure_digits: BTreeMap<usize, usize>,
    ring_closures_to_erase: Vec<usize>,
    chiral_tag_overrides: BTreeMap<AtomId, ChiralTag>,
    chiral_inversions: BTreeSet<AtomId>,
    chiral_permutations: BTreeMap<AtomId, u32>,
    broken_chiral_atoms: BTreeSet<AtomId>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct FragmentWritePlan {
    atoms: Vec<AtomId>,
    bonds: Vec<BondId>,
    rooted_at_atom: Option<AtomId>,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
struct FragmentWriteResult {
    smiles: String,
    atom_ordering: Vec<AtomId>,
    bond_ordering: Vec<BondId>,
}

#[derive(Debug, Clone, Copy, Default)]
struct SmilesWriteOverrides<'a> {
    atom_symbols: Option<&'a [String]>,
    bond_symbols: Option<&'a [String]>,
}

#[derive(Debug, Clone)]
struct CxWriteScope {
    atom_order: Vec<AtomId>,
    bond_order: Vec<BondId>,
}

impl CxWriteScope {
    fn full_molecule(molecule: &Molecule) -> Self {
        Self {
            atom_order: molecule.atoms().iter().map(|atom| atom.id()).collect(),
            bond_order: molecule.bonds().iter().map(|bond| bond.id()).collect(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[allow(dead_code)]
enum MolStackElem {
    Atom(AtomId),
    Bond(BondId, AtomId),
    Ring { bond: BondId, ring_idx: usize },
    BranchOpen,
    BranchClose,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
struct CanonicalTraversalResult {
    stack: Vec<MolStackElem>,
    traversal_ring_closure_bonds: Vec<bool>,
    chiral_tag_overrides: BTreeMap<AtomId, ChiralTag>,
    chiral_inversions: BTreeSet<AtomId>,
    chiral_permutations: BTreeMap<AtomId, u32>,
    broken_chiral_atoms: BTreeSet<AtomId>,
}

pub fn mol_to_smiles(
    molecule: &Molecule,
    params: &SmilesWriteParams,
) -> Result<String, SmilesWriteError> {
    mol_to_smiles_with_mode(molecule, params, SmilesOutputMode::PlainSmiles)
}

pub fn mol_to_cx_smiles(
    molecule: &Molecule,
    params: &SmilesWriteParams,
    fields: CxSmilesFields,
    restore_bond_dirs: RestoreBondDirOption,
) -> Result<String, SmilesWriteError> {
    mol_to_smiles_with_mode(
        molecule,
        params,
        SmilesOutputMode::CxSmiles {
            fields,
            restore_bond_dirs,
            include_stereo_groups: fields.contains(CxSmilesFields::ENHANCED_STEREO),
        },
    )
}

// BEGIN RDKIT CPP FUNCTION MolToRandomSmilesVect
// RDKit✔️✔️: std::vector<std::string> MolToRandomSmilesVect(
// RDKit✔️✔️:     const ROMol &mol, unsigned int numSmiles, unsigned int randomSeed,
// RDKit✔️✔️:     bool doIsomericSmiles, bool doKekule, bool allBondsExplicit,
// RDKit✔️✔️:     bool allHsExplicit) {
// RDKit✔️✔️:   if (randomSeed > 0) {
// RDKit✔️✔️:     getRandomGenerator(rdcast<int>(randomSeed));
// RDKit✔️✔️:   }
// RDKit✔️✔️:   std::vector<std::string> res;
// RDKit✔️✔️:   res.reserve(numSmiles);
// RDKit✔️✔️:   for (unsigned int i = 0; i < numSmiles; ++i) {
// RDKit✔️✔️:     bool canonical = false;
// RDKit✔️✔️:     int rootedAtAtom = -1;
// RDKit✔️✔️:     bool doRandom = true;
// RDKit✔️✔️:     res.push_back(MolToSmiles(mol, doIsomericSmiles, doKekule, rootedAtAtom,
// RDKit✔️✔️:                               canonical, allBondsExplicit, allHsExplicit,
// RDKit✔️✔️:                               doRandom));
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return res;
// RDKit✔️✔️: };
// END RDKIT CPP FUNCTION MolToRandomSmilesVect
pub fn mol_to_random_smiles_vect(
    molecule: &Molecule,
    num_smiles: usize,
    random_seed: u64,
    do_isomeric_smiles: bool,
    do_kekule: bool,
    all_bonds_explicit: bool,
    all_hydrogens_explicit: bool,
) -> Result<Vec<String>, SmilesWriteError> {
    let mut result = Vec::with_capacity(num_smiles);
    let mut stream_seed = if random_seed == 0 {
        next_unseeded_random_smiles_seed(0)
    } else {
        random_seed
    };
    for _ in 0..num_smiles {
        stream_seed = splitmix64(stream_seed);
        let params = SmilesWriteParams {
            do_isomeric_smiles,
            do_kekule,
            canonical: false,
            clean_stereo: true,
            all_bonds_explicit,
            all_hydrogens_explicit,
            do_random: true,
            rooted_at_atom: None,
            include_dative_bonds: true,
            ignore_atom_map_numbers: false,
        };
        result.push(with_random_smiles_seed(stream_seed, || {
            mol_to_smiles(molecule, &params)
        })?);
    }
    Ok(result)
}

fn with_random_smiles_seed<T>(
    seed: u64,
    f: impl FnOnce() -> Result<T, SmilesWriteError>,
) -> Result<T, SmilesWriteError> {
    RANDOM_SMILES_SEED.with(|cell| {
        let previous = cell.replace(seed);
        let result = f();
        cell.set(previous);
        result
    })
}

fn next_random_smiles_u64() -> u64 {
    RANDOM_SMILES_SEED.with(|cell| {
        let current = cell.get();
        let next = splitmix64(current);
        cell.set(next);
        current
    })
}

fn next_unseeded_random_smiles_seed(offset: u64) -> u64 {
    splitmix64(
        RANDOM_SMILES_COUNTER
            .fetch_add(0x9e37_79b9_7f4a_7c15, Ordering::Relaxed)
            .wrapping_add(offset),
    )
}

fn splitmix64(mut value: u64) -> u64 {
    value = value.wrapping_add(0x9e37_79b9_7f4a_7c15);
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}

fn mol_to_smiles_with_mode(
    molecule: &Molecule,
    params: &SmilesWriteParams,
    mode: SmilesOutputMode,
) -> Result<String, SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles
    // RDKit✔️✔️: std::string MolToSmiles(const ROMol &mol, const SmilesWriteParams &params,
    // RDKit✔️✔️:                         bool doingCXSmiles, bool includeStereoGroups) {
    // RDKit✔️✔️:   if (!mol.getNumAtoms()) {
    // RDKit✔️✔️:     return "";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   PRECONDITION(
    // RDKit✔️✔️:       params.rootedAtAtom < 0 ||
    // RDKit✔️✔️:           static_cast<unsigned int>(params.rootedAtAtom) < mol.getNumAtoms(),
    // RDKit✔️✔️:       "rootedAtAtom must be less than the number of atoms");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   int rootedAtAtom;
    // RDKit✔️✔️:   std::vector<int> fragsRootedAtAtom;
    // RDKit✔️✔️:   std::vector<std::vector<int>> fragsMolAtomMapping;
    // RDKit✔️✔️:   auto mols =
    // RDKit✔️✔️:       MolOps::getMolFrags(mol, false, nullptr, &fragsMolAtomMapping, false);
    // RDKit✔️✔️:   std::vector<std::vector<int>> fragsMolBondMapping;
    // RDKit✔️✔️:   std::vector<std::string> vfragsmi(mols.size());
    // RDKit✔️✔️:   std::vector<std::vector<RDKit::UINT>> allAtomOrdering;
    // RDKit✔️✔️:   std::vector<std::vector<RDKit::UINT>> allBondOrdering;
    // RDKit✔️✔️:   for (unsigned fragIdx = 0; fragIdx < mols.size(); fragIdx++) {
    // RDKit✔️✔️:     ROMol *tmol = mols[fragIdx].get();
    // RDKit✔️✔️:     std::vector<int> atomMapNums(tmol->getNumAtoms(), 0);
    // RDKit✔️✔️:     for (auto atom : tmol->atoms()) {
    // RDKit✔️✔️:       atom->updatePropertyCache(false);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (params.doIsomericSmiles) {
    // RDKit✔️✔️:       tmol->setProp(common_properties::_doIsoSmiles, 1);
    // RDKit✔️✔️:       if (!tmol->hasProp(common_properties::_StereochemDone)) {
    // RDKit✔️✔️:         MolOps::assignStereochemistry(*tmol, params.cleanStereo);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (params.canonical) {
    // RDKit✔️✔️:       Canon::rankMolAtoms(*tmol, ranks, breakTies, includeChirality,
    // RDKit✔️✔️:                           includeIsotopes, includeAtomMaps,
    // RDKit✔️✔️:                           includeChiralPresence, includeStereoGroups,
    // RDKit✔️✔️:                           useNonStereoRanks);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       std::iota(ranks.begin(), ranks.end(), 0);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     subSmi = SmilesWrite::FragmentSmilesConstruct(
    // RDKit✔️✔️:         *tmol, nextAtomIdx, colors, ranks, params, atomOrdering, bondOrdering);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (params.canonical) {
    // RDKit✔️✔️:     std::sort(tmp.begin(), tmp.end());
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     for (unsigned i = 0; i < vfragsmi.size(); ++i) {
    // RDKit✔️✔️:       result += vfragsmi[i];
    // RDKit✔️✔️:       if (i < vfragsmi.size() - 1) {
    // RDKit✔️✔️:         result += ".";
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   mol.setProp(common_properties::_smilesAtomOutputOrder, flattenedAtomOrdering,
    // RDKit✔️✔️:               true);
    // RDKit✔️✔️:   mol.setProp(common_properties::_smilesBondOutputOrder, flattenedBondOrdering,
    // RDKit✔️✔️:               true);
    // RDKit✔️✔️:   return result;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles
    validate_rooted_atom(molecule, params)?;
    if molecule.num_atoms() == 0 {
        return Ok(String::new());
    }

    let mut molecule = molecule.clone();

    let mut context = SmilesWriteContext::default();
    let mut fragment_results = Vec::new();
    let mut working_params = params.clone();

    let saved_atom_maps = match mode {
        SmilesOutputMode::PlainSmiles => {
            prepare_plain_smiles_molecule(&mut molecule, &working_params)?
        }
        SmilesOutputMode::CxSmiles {
            fields,
            restore_bond_dirs,
            include_stereo_groups,
        } => prepare_cx_smiles_molecule(
            &mut molecule,
            &mut working_params,
            fields,
            restore_bond_dirs,
            include_stereo_groups,
        )?,
    };

    let fragment_plans = collect_fragment_write_plans(&molecule, &working_params)?;
    let fragment_ranks = fragment_plans
        .iter()
        .map(|plan| rank_fragment_atoms_for_smiles(&molecule, plan, &working_params, mode))
        .collect::<Result<Vec<_>, _>>()?;
    // RDKit✔️✔️:       if (params.ignoreAtomMapNumbers) {
    // RDKit✔️✔️:         for (auto atom : tmol->atoms()) {
    // RDKit✔️✔️:           atom->setAtomMapNum(atomMapNums[atom->getIdx()]);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    //
    // RDKit restores atom maps immediately after canonical traversal ranking
    // and before FragmentSmilesConstruct(), whose doKekule branch calls
    // KekulizeFragment(). That means canonical kekulization still ranks with
    // the original atom maps even when ignoreAtomMapNumbers=true.
    if working_params.canonical {
        restore_atom_maps_after_canonical_smiles(&mut molecule, saved_atom_maps.as_deref());
    }
    if params.do_kekule {
        molecule = kekulize_for_smiles(&molecule)?;
    }
    // do_kekule already handled on the working molecule; keep the rest of the
    // writer on the post-kekulization topology without re-running that stage.
    working_params.do_kekule = false;
    for (plan, ranks) in fragment_plans.iter().zip(fragment_ranks.iter()) {
        if working_params.canonical {
            restore_atom_maps_after_canonical_smiles(&mut molecule, saved_atom_maps.as_deref());
        }
        fragment_results.push(write_fragment_smiles_with_ranks(
            &mut molecule,
            plan,
            &ranks,
            &working_params,
            SmilesWriteOverrides::default(),
            &mut context,
        )?);
        if working_params.canonical && saved_atom_maps.is_some() {
            let _ = stash_and_clear_atom_maps_for_smiles(&mut molecule, &working_params);
        }
    }
    if working_params.canonical {
        restore_atom_maps_after_canonical_smiles(&mut molecule, saved_atom_maps.as_deref());
    }

    let mut result = assemble_fragment_smiles(fragment_results, &working_params, &mut context)?;
    if let SmilesOutputMode::CxSmiles { fields, .. } = mode {
        let scope = CxWriteScope {
            atom_order: context.atom_output_order.clone(),
            bond_order: context.bond_output_order.clone(),
        };
        let cx_extension = get_cx_extensions_scoped(&molecule, fields, &scope)?;
        if !cx_extension.is_empty() {
            result.push(' ');
            result.push_str(&cx_extension);
        }
    }
    Ok(result)
}

fn prepare_plain_smiles_molecule(
    molecule: &mut Molecule,
    params: &SmilesWriteParams,
) -> Result<Option<Vec<Option<u32>>>, SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles fragment preparation section
    // RDKit✔️✔️:     // update property cache
    // RDKit✔️✔️:     std::vector<int> atomMapNums(tmol->getNumAtoms(), 0);
    // RDKit✔️✔️:     for (auto atom : tmol->atoms()) {
    // RDKit✔️✔️:       if (params.ignoreAtomMapNumbers) {
    // RDKit✔️✔️:         atomMapNums[atom->getIdx()] = atom->getAtomMapNum();
    // RDKit✔️✔️:         atom->setAtomMapNum(0);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       atom->updatePropertyCache(false);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // clean up the chirality on any atom that is marked as chiral,
    // RDKit✔️✔️:     // but that should not be:
    // RDKit✔️✔️:     if (params.doIsomericSmiles) {
    // RDKit✔️✔️:       tmol->setProp(common_properties::_doIsoSmiles, 1);
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (!tmol->hasProp(common_properties::_StereochemDone)) {
    // RDKit✔️✔️:         MolOps::assignStereochemistry(*tmol, params.cleanStereo);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!doingCXSmiles || !includeStereoGroups) {
    // RDKit✔️✔️:       std::vector<StereoGroup> noStereoGroups;
    // RDKit✔️✔️:       tmol->setStereoGroups(noStereoGroups);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!doingCXSmiles) {
    // RDKit✔️✔️:       for (auto bond : tmol->bonds()) {
    // RDKit✔️✔️:         if (bond->getBondDir() == Bond::BondDir::UNKNOWN ||
    // RDKit✔️✔️:             bond->getBondDir() == Bond::BondDir::EITHERDOUBLE) {
    // RDKit✔️✔️:           bond->setBondDir(Bond::BondDir::NONE);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         if (bond->getStereo() == Bond::BondStereo::STEREOANY) {
    // RDKit✔️✔️:           bond->setStereo(Bond::BondStereo::STEREONONE);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (doingCXSmiles || !params.includeDativeBonds) {
    // RDKit✔️✔️:       for (auto bond : tmol->bonds()) {
    // RDKit✔️✔️:         if (bond->getBondType() == Bond::DATIVE) {
    // RDKit✔️✔️:           bond->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:           bond->getBeginAtom()->calcExplicitValence(false);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles fragment preparation section
    let saved_atom_maps = stash_and_clear_atom_maps_for_smiles(molecule, params);
    if is_minimal_plain_smiles_path(params) && validate_minimal_plain_smiles_molecule(molecule) {
        return Ok(saved_atom_maps);
    }
    clear_fragment_temp_molecule_computed_stereo_props_for_writer(molecule);
    update_property_cache_for_smiles(molecule)?;
    if params.do_isomeric_smiles {
        if molecule.prop("_StereochemDone").is_none() {
            assign_stereochemistry_for_smiles(molecule, params.clean_stereo)?;
        }
    }
    // Kekulization is handled upstream in mol_to_smiles_with_mode before
    // this function is called.
    if params.do_random {
        // Random SMILES uses non-canonical traversal with randomized
        // bond ordering at each atom. Continue through the standard
        // preparation path; randomization happens in the fragment
        // traversal step.
    }
    if !params.include_dative_bonds {
        normalize_dative_bonds_for_plain_smiles(molecule)?;
    }
    if !params.do_isomeric_smiles {
        // RDKit plain non-isomeric SMILES suppresses bond-direction output.
        crate::notation::smiles::clear_all_bond_dir_flags(molecule);
    }
    remove_plain_smiles_only_cx_state(molecule)?;
    Ok(saved_atom_maps)
}

fn prepare_cx_smiles_molecule(
    molecule: &mut Molecule,
    params: &mut SmilesWriteParams,
    fields: CxSmilesFields,
    restore_bond_dirs: RestoreBondDirOption,
    include_stereo_groups: bool,
) -> Result<Option<Vec<Option<u32>>>, SmilesWriteError> {
    let saved_atom_maps = stash_and_clear_atom_maps_for_smiles(molecule, params);
    // Kekulization is handled upstream in mol_to_smiles_with_mode.
    if is_minimal_plain_smiles_path(params) && validate_minimal_plain_smiles_molecule(molecule) {
        // CX still needs CX-specific cleanup below; the fast path only skips
        // property/stereo preparation for the simplest typed molecule state.
    } else {
        clear_fragment_temp_molecule_computed_stereo_props_for_writer(molecule);
        update_property_cache_for_smiles(molecule)?;
        if params.do_isomeric_smiles {
            if molecule.prop("_StereochemDone").is_none() {
                assign_stereochemistry_for_smiles(molecule, params.clean_stereo)?;
            }
        }
    }
    normalize_dative_bonds_for_cx_smiles(molecule)?;
    normalize_hydrogen_bonds_for_cx_smiles(molecule)?;
    apply_cx_bond_direction_policy(molecule, restore_bond_dirs)?;
    if params.clean_stereo {
        if molecule.prop("_StereochemDone").is_none() {
            assign_stereochemistry_for_smiles(molecule, true)?;
        }
        cleanup_stereo_groups_for_cx_smiles(molecule)?;
    }
    if include_stereo_groups {
        canonicalize_enhanced_stereo_for_smiles(molecule)?;
    }
    validate_cx_extension_plan(fields)?;
    Ok(saved_atom_maps)
}

fn stash_and_clear_atom_maps_for_smiles(
    molecule: &mut Molecule,
    params: &SmilesWriteParams,
) -> Option<Vec<Option<u32>>> {
    if !params.ignore_atom_map_numbers {
        return None;
    }
    let topology = molecule.topology_block_mut();
    let saved = topology
        .atoms
        .iter()
        .map(|atom| atom.atom_map())
        .collect::<Vec<_>>();
    for atom in &mut topology.atoms {
        atom.set_atom_map(None);
    }
    Some(saved)
}

fn restore_atom_maps_after_canonical_smiles(
    molecule: &mut Molecule,
    saved_atom_maps: Option<&[Option<u32>]>,
) {
    let Some(saved_atom_maps) = saved_atom_maps else {
        return;
    };
    let topology = molecule.topology_block_mut();
    for (atom, atom_map) in topology
        .atoms
        .iter_mut()
        .zip(saved_atom_maps.iter().copied())
    {
        atom.set_atom_map(atom_map);
    }
}

fn collect_fragment_write_plans(
    molecule: &Molecule,
    params: &SmilesWriteParams,
) -> Result<Vec<FragmentWritePlan>, SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles fragment collection section
    // RDKit✔️✔️:   int rootedAtAtom;
    // RDKit✔️✔️:   std::vector<int> fragsRootedAtAtom;
    // RDKit✔️✔️:   std::vector<std::vector<int>> fragsMolAtomMapping;
    // RDKit✔️✔️:   auto mols =
    // RDKit✔️✔️:       MolOps::getMolFrags(mol, false, nullptr, &fragsMolAtomMapping, false);
    // RDKit✔️✔️:   // we got the mapping between fragments and atoms; repeat that for bonds
    // RDKit✔️✔️:   std::vector<std::vector<int>> fragsMolBondMapping;
    // RDKit✔️✔️:   boost::dynamic_bitset<> atsPresent(mol.getNumAtoms());
    // RDKit✔️✔️:   std::vector<int> bondsInFrag;
    // RDKit✔️✔️:   bondsInFrag.reserve(mol.getNumBonds());
    // RDKit✔️✔️:   for (const auto &atsInFrag : fragsMolAtomMapping) {
    // RDKit✔️✔️:     atsPresent.reset();
    // RDKit✔️✔️:     bondsInFrag.clear();
    // RDKit✔️✔️:     for (auto aidx : atsInFrag) {
    // RDKit✔️✔️:       atsPresent.set(aidx);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     rootedAtAtom = -1;
    // RDKit✔️✔️:     if (params.rootedAtAtom >= 0 && atsPresent[params.rootedAtAtom]) {
    // RDKit✔️✔️:       rootedAtAtom = params.rootedAtAtom - atsPresent.find_first();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     fragsRootedAtAtom.push_back(rootedAtAtom);
    // RDKit✔️✔️:
    // RDKit✔️✔️:     for (const auto bnd : mol.bonds()) {
    // RDKit✔️✔️:       if (atsPresent[bnd->getBeginAtomIdx()] &&
    // RDKit✔️✔️:           atsPresent[bnd->getEndAtomIdx()]) {
    // RDKit✔️✔️:         bondsInFrag.push_back(bnd->getIdx());
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     fragsMolBondMapping.push_back(bondsInFrag);
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles fragment collection section
    let atom_to_fragment = crate::notation::fragment::get_fragment_atom_mapping(molecule);
    if atom_to_fragment.is_empty() {
        return Ok(Vec::new());
    }
    let fragment_count = atom_to_fragment.iter().copied().max().unwrap_or(0) + 1;
    let mut fragment_atoms = vec![Vec::new(); fragment_count];
    for (atom_idx, fragment_idx) in atom_to_fragment.iter().copied().enumerate() {
        fragment_atoms[fragment_idx].push(AtomId::new(atom_idx));
    }
    let mut fragment_bonds = vec![Vec::new(); fragment_count];
    for bond in molecule.bonds() {
        let begin_fragment = atom_to_fragment[bond.begin().index()];
        let end_fragment = atom_to_fragment[bond.end().index()];
        if begin_fragment == end_fragment {
            fragment_bonds[begin_fragment].push(bond.id());
        }
    }
    let mut plans = Vec::with_capacity(fragment_count);
    for fragment_idx in 0..fragment_count {
        let atoms = std::mem::take(&mut fragment_atoms[fragment_idx]);
        let rooted_at_atom = params
            .rooted_at_atom
            .map(AtomId::new)
            .filter(|root| atom_to_fragment[root.index()] == fragment_idx);
        plans.push(FragmentWritePlan {
            bonds: std::mem::take(&mut fragment_bonds[fragment_idx]),
            atoms,
            rooted_at_atom,
        });
    }
    Ok(plans)
}

fn write_fragment_smiles(
    molecule: &mut Molecule,
    plan: &FragmentWritePlan,
    params: &SmilesWriteParams,
    mode: SmilesOutputMode,
    overrides: SmilesWriteOverrides<'_>,
    context: &mut SmilesWriteContext,
) -> Result<FragmentWriteResult, SmilesWriteError> {
    let ranks = rank_fragment_atoms_for_smiles(molecule, plan, params, mode)?;
    write_fragment_smiles_with_ranks(molecule, plan, &ranks, params, overrides, context)
}

fn write_fragment_smiles_with_ranks(
    molecule: &mut Molecule,
    plan: &FragmentWritePlan,
    ranks: &[usize],
    params: &SmilesWriteParams,
    overrides: SmilesWriteOverrides<'_>,
    context: &mut SmilesWriteContext,
) -> Result<FragmentWriteResult, SmilesWriteError> {
    let start_atom = choose_fragment_start_atom(plan, &ranks, params)?;
    fragment_smiles_construct(
        molecule, plan, start_atom, &ranks, params, overrides, context,
    )
}

fn fragment_smiles_construct(
    molecule: &mut Molecule,
    plan: &FragmentWritePlan,
    start_atom: AtomId,
    ranks: &[usize],
    params: &SmilesWriteParams,
    overrides: SmilesWriteOverrides<'_>,
    context: &mut SmilesWriteContext,
) -> Result<FragmentWriteResult, SmilesWriteError> {
    // Full-molecule kekulization is handled before fragment planning through
    // the registered operation pipeline.
    if params.canonical && params.do_isomeric_smiles {
        canonicalize_enhanced_stereo_for_smiles(molecule)?;
    }
    let traversal =
        canonicalize_fragment_stack(molecule, plan, start_atom, ranks, params, overrides)?;
    canonicalize_double_bond_directions_for_writer(
        molecule,
        &traversal.stack,
        &traversal.traversal_ring_closure_bonds,
    )?;
    context.chiral_tag_overrides.extend(
        traversal
            .chiral_tag_overrides
            .iter()
            .map(|(atom, tag)| (*atom, *tag)),
    );
    context
        .chiral_inversions
        .extend(traversal.chiral_inversions.iter().copied());
    context.chiral_permutations.extend(
        traversal
            .chiral_permutations
            .iter()
            .map(|(atom, permutation)| (*atom, *permutation)),
    );
    context
        .broken_chiral_atoms
        .extend(traversal.broken_chiral_atoms.iter().copied());
    write_mol_stack(molecule, &traversal.stack, params, overrides, context)
}

fn rank_fragment_atoms_for_smiles(
    molecule: &Molecule,
    plan: &FragmentWritePlan,
    params: &SmilesWriteParams,
    mode: SmilesOutputMode,
) -> Result<Vec<usize>, SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION SmilesWrite non-canonical rank initialization
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       std::iota(ranks.begin(), ranks.end(), 0);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (unsigned int i = 0; i < tmol.getNumAtoms(); ++i) {
    // RDKit✔️✔️:       ranks[i] = i;
    // RDKit✔️✔️:     }
    // END RDKIT CPP FUNCTION SmilesWrite non-canonical rank initialization
    // RDKit canonical mode still computes canonical ranks when rootedAtAtom
    // is provided; rootedAtAtom only overrides traversal start.
    if params.canonical && !params.do_random {
        return rank_mol_atoms_for_smiles(molecule, plan, params, mode);
    }
    let _ = molecule;
    Ok(plan.atoms.iter().map(|atom| atom.index()).collect())
}

fn rank_mol_atoms_for_smiles(
    molecule: &Molecule,
    plan: &FragmentWritePlan,
    params: &SmilesWriteParams,
    mode: SmilesOutputMode,
) -> Result<Vec<usize>, SmilesWriteError> {
    let _stage = SmilesPlanStage::LongTermCanonicalRanking;
    let _ = mode;
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles canonical rank options
    // RDKit✔️✔️:       const bool includeChiralPresence = false;
    // RDKit✔️✔️:       const bool includeIsotopes = params.doIsomericSmiles;
    // RDKit✔️✔️:       ;
    // RDKit✔️✔️:       const bool includeChirality = params.doIsomericSmiles;
    // RDKit✔️✔️:       ;
    // RDKit✔️✔️:       const bool includeStereoGroups = params.doIsomericSmiles;
    // RDKit✔️✔️:       ;
    // RDKit✔️✔️:       const bool useNonStereoRanks = false;
    // RDKit✔️✔️:       const bool includeAtomMaps = true;
    // RDKit✔️✔️:
    // RDKit✔️✔️:       Canon::rankMolAtoms(*tmol, ranks, breakTies, includeChirality,
    // RDKit✔️✔️:                           includeIsotopes, includeAtomMaps,
    // RDKit✔️✔️:                           includeChiralPresence, includeStereoGroups,
    // RDKit✔️✔️:                           useNonStereoRanks);
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles canonical rank options
    let ranks = crate::canon_rank::rank_mol_atoms_with_options(
        molecule,
        crate::canon_rank::CanonicalRankOptions {
            break_ties: true,
            include_chirality: params.do_isomeric_smiles,
            include_isotopes: params.do_isomeric_smiles,
            include_atom_maps: true,
            include_chiral_presence: false,
            include_stereo_groups: params.do_isomeric_smiles,
            use_non_stereo_ranks: false,
            include_ring_stereo: params.do_isomeric_smiles,
            chirality_rings_use_ring_stereo: true,
        },
    )?;
    Ok(plan.atoms.iter().map(|atom| ranks[atom.index()]).collect())
}

fn choose_fragment_start_atom(
    plan: &FragmentWritePlan,
    ranks: &[usize],
    params: &SmilesWriteParams,
) -> Result<AtomId, SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles start atom selection section
    // RDKit✔️✔️:     // find the next atom for a traverse
    // RDKit✔️✔️:     if (params.doRandom && rootedAtAtom == -1) {
    // RDKit✔️✔️:       rootedAtAtom = getRandomGenerator()() % tmol->getNumAtoms();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (rootedAtAtom >= 0) {
    // RDKit✔️✔️:       nextAtomIdx = rootedAtAtom;
    // RDKit✔️✔️:       rootedAtAtom = -1;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       unsigned int nextRank = nAtoms + 1;
    // RDKit✔️✔️:       for (unsigned int i = 0; i < nAtoms; i++) {
    // RDKit✔️✔️:         if (colors[i] == Canon::WHITE_NODE && ranks[i] < nextRank) {
    // RDKit✔️✔️:           nextRank = ranks[i];
    // RDKit✔️✔️:           nextAtomIdx = i;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     CHECK_INVARIANT(nextAtomIdx >= 0, "no start atom found");
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles start atom selection section
    if let Some(root) = plan.rooted_at_atom {
        return Ok(root);
    }
    if params.do_random {
        let idx = (next_random_smiles_u64() as usize) % plan.atoms.len();
        return Ok(plan.atoms[idx]);
    }
    let (idx, _) = match ranks.iter().enumerate().min_by_key(|(_, rank)| **rank) {
        Some(pair) => pair,
        // [deferred] Empty ranks: this is a defensive guard for impossible
        // state (fragment with no atoms). Should never fire in practice
        // since choose_fragment_start_atom is only called on non-empty plans.
        None => {
            return invariant_stage_error(
                SmilesPlanStage::ShortTermAtomWriter,
                "choose_fragment_start_atom() called with empty canonical rank scope",
            );
        }
    };
    Ok(plan.atoms[idx])
}

fn canonicalize_fragment_stack(
    molecule: &Molecule,
    plan: &FragmentWritePlan,
    start_atom: AtomId,
    ranks: &[usize],
    params: &SmilesWriteParams,
    overrides: SmilesWriteOverrides<'_>,
) -> Result<CanonicalTraversalResult, SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION Canon::canonicalizeFragment call site
    // RDKit✔️✔️:     subSmi = SmilesWrite::FragmentSmilesConstruct(
    // RDKit✔️✔️:         *tmol, nextAtomIdx, colors, ranks, params, atomOrdering, bondOrdering);
    // RDKit✔️✔️: Canon::canonicalizeFragment(mol, atomIdx, colors, ranks, molStack,
    // RDKit✔️✔️:                           atomsInPlay, bondsInPlay, bondSymbols,
    // RDKit✔️✔️:                           params.doIsomericSmiles, params.doRandom);
    // END RDKIT CPP FUNCTION Canon::canonicalizeFragment call site
    canonical_dfs_traversal(
        molecule,
        plan,
        start_atom,
        ranks,
        params.do_isomeric_smiles,
        params.clean_stereo,
        params.do_random,
        overrides.bond_symbols,
    )
}

fn write_mol_stack(
    molecule: &Molecule,
    stack: &[MolStackElem],
    params: &SmilesWriteParams,
    overrides: SmilesWriteOverrides<'_>,
    context: &mut SmilesWriteContext,
) -> Result<FragmentWriteResult, SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION FragmentSmilesConstruct MolStack emission section
    // RDKit✔️✔️:   Bond *bond = nullptr;
    // RDKit✔️✔️:   for (auto &mSE : molStack) {
    // RDKit✔️✔️:     switch (mSE.type) {
    // RDKit✔️✔️:       case Canon::MOL_STACK_ATOM:
    // RDKit✔️✔️:         for (auto rclosure : ringClosuresToErase) {
    // RDKit✔️✔️:           ringClosureMap.erase(rclosure);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         ringClosuresToErase.clear();
    // RDKit✔️✔️:         if (!atomSymbols) {
    // RDKit✔️✔️:           res << GetAtomSmiles(mSE.obj.atom, params);
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           res << (*atomSymbols)[mSE.obj.atom->getIdx()];
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         atomOrdering.push_back(mSE.obj.atom->getIdx());
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case Canon::MOL_STACK_BOND:
    // RDKit✔️✔️:         bond = mSE.obj.bond;
    // RDKit✔️✔️:         if (!bondSymbols) {
    // RDKit✔️✔️:           res << GetBondSmiles(bond, params, mSE.number);
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           res << (*bondSymbols)[bond->getIdx()];
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         bondOrdering.push_back(bond->getIdx());
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case Canon::MOL_STACK_RING:
    // RDKit✔️✔️:       case Canon::MOL_STACK_BRANCH_OPEN:
    // RDKit✔️✔️:       case Canon::MOL_STACK_BRANCH_CLOSE:
    // RDKit✔️✔️:       default:
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res.str();
    // END RDKIT CPP FUNCTION FragmentSmilesConstruct MolStack emission section
    let mut result = FragmentWriteResult::default();
    for item in stack {
        match *item {
            MolStackElem::Atom(atom) => {
                for ring_closure in context.ring_closures_to_erase.drain(..) {
                    context.ring_closure_digits.remove(&ring_closure);
                }
                if let Some(atom_symbols) = overrides.atom_symbols {
                    result.smiles.push_str(&atom_symbols[atom.index()]);
                } else {
                    result
                        .smiles
                        .push_str(&build_atom_smiles(molecule, atom, params, context)?);
                }
                result.atom_ordering.push(atom);
            }
            MolStackElem::Bond(bond, atom_to_left) => {
                if let Some(bond_symbols) = overrides.bond_symbols {
                    result.smiles.push_str(&bond_symbols[bond.index()]);
                } else {
                    result.smiles.push_str(&build_bond_smiles(
                        molecule,
                        bond,
                        atom_to_left,
                        params,
                    )?);
                }
                result.bond_ordering.push(bond);
            }
            MolStackElem::Ring { ring_idx, .. } => {
                write_ring_closure(&mut result.smiles, ring_idx, context)?;
            }
            MolStackElem::BranchOpen => {
                result.smiles.push('(');
            }
            MolStackElem::BranchClose => {
                result.smiles.push(')');
            }
        }
    }
    Ok(result)
}

// BEGIN RDKIT CPP FUNCTION MolFragmentToSmiles
// RDKit✔️✔️: std::string MolFragmentToSmiles(const ROMol &mol,
// RDKit✔️✔️:                                 const SmilesWriteParams &params,
// RDKit✔️✔️:                                 const std::vector<int> &atomsToUse,
// RDKit✔️✔️:                                 const std::vector<int> *bondsToUse,
// RDKit✔️✔️:                                 const std::vector<std::string> *atomSymbols,
// RDKit✔️✔️:                                 const std::vector<std::string> *bondSymbols) {
// RDKit✔️✔️:   PRECONDITION(atomsToUse.size(), "no atoms provided");
// RDKit✔️✔️:   if (!mol.getNumAtoms()) { return ""; }
// RDKit✔️✔️:   int rootedAtAtom = params.rootedAtAtom;
// RDKit✔️✔️:   ROMol tmol(mol, true);  // copy molecule
// RDKit✔️✔️:   std::string res;
// RDKit✔️✔️:   // compute bondsInPlay from atomsToUse
// RDKit✔️✔️:   // then FragmentSmilesConstruct with atomSymbols/bondSymbols
// RDKit✔️✔️:   return res;
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION MolFragmentToSmiles
pub fn mol_fragment_to_smiles(
    molecule: &Molecule,
    params: &SmilesWriteParams,
    atoms_to_use: &[usize],
    bonds_to_use: Option<&[usize]>,
    atom_symbols: Option<&[String]>,
    bond_symbols: Option<&[String]>,
) -> Result<String, SmilesWriteError> {
    validate_fragment_api_inputs(
        molecule,
        params,
        atoms_to_use,
        bonds_to_use,
        atom_symbols,
        bond_symbols,
    )?;
    if molecule.num_atoms() == 0 || atoms_to_use.is_empty() {
        return Ok(String::new());
    }

    // BEGIN RDKIT CPP FUNCTION MolFragmentToSmiles fragment bitset/output section
    // RDKit✔️✔️:   boost::dynamic_bitset<> atomsInPlay(mol.getNumAtoms(), 0);
    // RDKit✔️✔️:   for (auto aidx : atomsToUse) { atomsInPlay.set(aidx); }
    // RDKit✔️✔️:   boost::dynamic_bitset<> bondsInPlay(mol.getNumBonds(), 0);
    // RDKit✔️✔️:   if (bondsToUse) { for (auto bidx : *bondsToUse) { bondsInPlay.set(bidx); } }
    // RDKit✔️✔️:   else {
    // RDKit✔️✔️:     PRECONDITION(
    // RDKit✔️✔️:         params.rootedAtAtom < 0 || MolOps::getMolFrags(mol).size() == 1,
    // RDKit✔️✔️:         "rootedAtAtom can only be used with molecules that have a single fragment");
    // RDKit✔️✔️:     for (auto aidx : atomsToUse) { ... if (atomsInPlay[other]) bondsInPlay.set(...); }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   while (colorIt != colors.end()) {
    // RDKit✔️✔️:     ... FragmentSmilesConstruct(..., &atomsInPlay, &bondsInPlay, atomSymbols, bondSymbols);
    // RDKit✔️✔️:     if (colorIt != colors.end()) { res += "."; }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION MolFragmentToSmiles fragment bitset/output section
    let mut molecule = if params.do_kekule {
        kekulize_for_smiles(molecule)?
    } else {
        molecule.clone()
    };
    let mut working_params = params.clone();
    working_params.do_kekule = false;
    let saved_atom_maps = prepare_plain_smiles_molecule(&mut molecule, &working_params)?;

    let mut plans =
        collect_fragment_api_write_plans(&molecule, &working_params, atoms_to_use, bonds_to_use)?;
    if working_params.canonical {
        restore_atom_maps_after_canonical_smiles(&mut molecule, saved_atom_maps.as_deref());
        plans.sort_by_key(|plan| {
            plan.atoms
                .iter()
                .map(|atom| atom.index())
                .min()
                .unwrap_or(usize::MAX)
        });
    }

    let overrides = SmilesWriteOverrides {
        atom_symbols,
        bond_symbols,
    };
    let mut context = SmilesWriteContext::default();
    let mut results = Vec::new();
    for plan in &plans {
        results.push(write_fragment_smiles(
            &mut molecule,
            plan,
            &working_params,
            SmilesOutputMode::PlainSmiles,
            overrides,
            &mut context,
        )?);
    }
    assemble_fragment_smiles(results, &working_params, &mut context)
}

pub fn mol_fragment_to_cx_smiles(
    molecule: &Molecule,
    params: &SmilesWriteParams,
    atoms_to_use: &[usize],
    bonds_to_use: Option<&[usize]>,
    atom_symbols: Option<&[String]>,
    bond_symbols: Option<&[String]>,
    fields: CxSmilesFields,
) -> Result<String, SmilesWriteError> {
    validate_fragment_api_inputs(
        molecule,
        params,
        atoms_to_use,
        bonds_to_use,
        atom_symbols,
        bond_symbols,
    )?;
    let mut context = SmilesWriteContext::default();
    let smiles = mol_fragment_to_smiles_with_context(
        molecule,
        params,
        atoms_to_use,
        bonds_to_use,
        atom_symbols,
        bond_symbols,
        &mut context,
    )?;
    let scope = CxWriteScope {
        atom_order: context.atom_output_order,
        bond_order: context.bond_output_order,
    };
    let cx_extension = get_cx_extensions_scoped(molecule, fields, &scope)?;
    if cx_extension.is_empty() {
        Ok(smiles)
    } else {
        Ok(format!("{smiles} {cx_extension}"))
    }
}

fn mol_fragment_to_smiles_with_context(
    molecule: &Molecule,
    params: &SmilesWriteParams,
    atoms_to_use: &[usize],
    bonds_to_use: Option<&[usize]>,
    atom_symbols: Option<&[String]>,
    bond_symbols: Option<&[String]>,
    context: &mut SmilesWriteContext,
) -> Result<String, SmilesWriteError> {
    if molecule.num_atoms() == 0 || atoms_to_use.is_empty() {
        return Ok(String::new());
    }

    let mut molecule = if params.do_kekule {
        kekulize_for_smiles(molecule)?
    } else {
        molecule.clone()
    };
    let mut working_params = params.clone();
    working_params.do_kekule = false;
    let saved_atom_maps = prepare_plain_smiles_molecule(&mut molecule, &working_params)?;

    let mut plans =
        collect_fragment_api_write_plans(&molecule, &working_params, atoms_to_use, bonds_to_use)?;
    if working_params.canonical {
        restore_atom_maps_after_canonical_smiles(&mut molecule, saved_atom_maps.as_deref());
        plans.sort_by_key(|plan| {
            plan.atoms
                .iter()
                .map(|atom| atom.index())
                .min()
                .unwrap_or(usize::MAX)
        });
    }

    let overrides = SmilesWriteOverrides {
        atom_symbols,
        bond_symbols,
    };
    let mut results = Vec::new();
    for plan in &plans {
        results.push(write_fragment_smiles(
            &mut molecule,
            plan,
            &working_params,
            SmilesOutputMode::PlainSmiles,
            overrides,
            context,
        )?);
    }
    assemble_fragment_smiles(results, &working_params, context)
}

fn collect_fragment_api_write_plans(
    molecule: &Molecule,
    params: &SmilesWriteParams,
    atoms_to_use: &[usize],
    bonds_to_use: Option<&[usize]>,
) -> Result<Vec<FragmentWritePlan>, SmilesWriteError> {
    let atom_set = atoms_to_use.iter().copied().collect::<BTreeSet<_>>();
    let bond_set = if let Some(bonds_to_use) = bonds_to_use {
        bonds_to_use.iter().copied().collect::<BTreeSet<_>>()
    } else {
        molecule
            .bonds()
            .iter()
            .filter(|bond| {
                atom_set.contains(&bond.begin().index()) && atom_set.contains(&bond.end().index())
            })
            .map(|bond| bond.id().index())
            .collect::<BTreeSet<_>>()
    };

    let mut seen = BTreeSet::new();
    let mut plans = Vec::new();
    for &start in atoms_to_use {
        if seen.contains(&start) {
            continue;
        }
        let mut stack = vec![AtomId::new(start)];
        let mut atoms = Vec::new();
        let mut bonds = BTreeSet::new();
        while let Some(atom) = stack.pop() {
            if !seen.insert(atom.index()) {
                continue;
            }
            atoms.push(atom);
            for bond in molecule.bonds() {
                if !bond_set.contains(&bond.id().index()) {
                    continue;
                }
                let Some(other) = bond_other_atom(bond, atom) else {
                    continue;
                };
                if !atom_set.contains(&other.index()) {
                    continue;
                }
                bonds.insert(bond.id());
                if !seen.contains(&other.index()) {
                    stack.push(other);
                }
            }
        }
        atoms.sort_by_key(|atom| atom.index());
        let bonds = bonds.into_iter().collect::<Vec<_>>();
        let rooted_at_atom = params
            .rooted_at_atom
            .map(AtomId::new)
            .filter(|root| atoms.contains(root));
        plans.push(FragmentWritePlan {
            atoms,
            bonds,
            rooted_at_atom,
        });
    }
    Ok(plans)
}

pub fn get_atom_smiles(
    molecule: &Molecule,
    atom: usize,
    params: &SmilesWriteParams,
) -> Result<String, SmilesWriteError> {
    validate_atom_index(molecule, atom)?;
    get_atom_smiles_impl(
        molecule,
        AtomId::new(atom),
        params,
        None,
        false,
        None,
        false,
    )
}

fn get_atom_smiles_with_context(
    molecule: &Molecule,
    atom: AtomId,
    params: &SmilesWriteParams,
    context: &SmilesWriteContext,
) -> Result<String, SmilesWriteError> {
    get_atom_smiles_impl(
        molecule,
        atom,
        params,
        context.chiral_tag_overrides.get(&atom).copied(),
        context.chiral_inversions.contains(&atom),
        context.chiral_permutations.get(&atom).copied(),
        context.broken_chiral_atoms.contains(&atom),
    )
}

fn get_atom_smiles_impl(
    molecule: &Molecule,
    atom_id: AtomId,
    params: &SmilesWriteParams,
    chiral_tag_override: Option<ChiralTag>,
    invert_chirality: bool,
    chiral_permutation_override: Option<u32>,
    broken_chirality: bool,
) -> Result<String, SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION GetAtomSmiles
    // RDKit✔️✔️: std::string GetAtomSmiles(const Atom *atom, const SmilesWriteParams &params) {
    // RDKit✔️✔️:   PRECONDITION(atom, "bad atom");
    // RDKit✔️✔️:   std::string res;
    // RDKit✔️✔️:   int fc = atom->getFormalCharge();
    // RDKit✔️✔️:   int num = atom->getAtomicNum();
    // RDKit✔️✔️:   int isotope = atom->getIsotope();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::string symb;
    // RDKit✔️✔️:   bool hasCustomSymbol =
    // RDKit✔️✔️:       atom->getPropIfPresent(common_properties::smilesSymbol, symb);
    // RDKit✔️✔️:   if (!hasCustomSymbol) {
    // RDKit✔️✔️:     symb = PeriodicTable::getTable()->getElementSymbol(num);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // check for atomic stereochemistry
    // RDKit✔️✔️:   std::string atString;
    // RDKit✔️✔️:   if (params.doIsomericSmiles) {
    // RDKit✔️✔️:     if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit✔️✔️:         !atom->hasProp(common_properties::_brokenChirality)) {
    // RDKit✔️✔️:       atString = getAtomChiralityInfo(atom);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   bool needsBracket = true;
    // RDKit✔️✔️:   if (!hasCustomSymbol && !params.allHsExplicit) {
    // RDKit✔️✔️:     needsBracket = atomNeedsBracket(atom, atString, params);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (needsBracket) {
    // RDKit✔️✔️:     res += "[";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (isotope && params.doIsomericSmiles) {
    // RDKit✔️✔️:     res += std::to_string(isotope);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!params.doKekule && atom->getIsAromatic() && symb[0] >= 'A' &&
    // RDKit✔️✔️:       symb[0] <= 'Z') {
    // RDKit✔️✔️:     symb[0] = tolower(symb[0]);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   res += symb;
    // RDKit✔️✔️:   res += atString;
    // RDKit✔️✔️:   if (needsBracket) {
    // RDKit✔️✔️:     unsigned int totNumHs = atom->getTotalNumHs();
    // RDKit✔️✔️:     if (totNumHs > 0) {
    // RDKit✔️✔️:       res += "H";
    // RDKit✔️✔️:       if (totNumHs > 1) {
    // RDKit✔️✔️:         res += std::to_string(totNumHs);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (fc > 0) {
    // RDKit✔️✔️:       res += "+";
    // RDKit✔️✔️:       if (fc > 1) {
    // RDKit✔️✔️:         res += std::to_string(fc);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (fc < 0) {
    // RDKit✔️✔️:       if (fc < -1) {
    // RDKit✔️✔️:         res += std::to_string(fc);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         res += "-";
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     int mapNum;
    // RDKit✔️✔️:     if (atom->getPropIfPresent(common_properties::molAtomMapNumber, mapNum)) {
    // RDKit✔️✔️:       res += ":";
    // RDKit✔️✔️:       res += std::to_string(mapNum);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     res += "]";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::string label;
    // RDKit✔️✔️:   if (atom->getPropIfPresent(common_properties::_supplementalSmilesLabel,
    // RDKit✔️✔️:                              label)) {
    // RDKit✔️✔️:     res += label;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION GetAtomSmiles
    let chirality = if params.do_isomeric_smiles && !broken_chirality {
        get_atom_chirality_info_with_inversion(
            molecule,
            atom_id,
            chiral_tag_override,
            invert_chirality,
            chiral_permutation_override,
        )?
    } else {
        String::new()
    };
    let atom = &molecule.atoms()[atom_id.index()];
    let custom_symbol = atom.prop("smilesSymbol");
    let has_custom_symbol = custom_symbol.is_some();
    let needs_bracket = if has_custom_symbol || params.all_hydrogens_explicit {
        true
    } else {
        atom_needs_bracket(molecule, atom_id, &chirality, params)?
    };
    let raw_symbol = custom_symbol.unwrap_or(element_symbol(atom.atomic_number())?);
    let lowered_symbol;
    let symbol: &str = if !params.do_kekule
        && atom.is_aromatic()
        && raw_symbol
            .as_bytes()
            .first()
            .is_some_and(u8::is_ascii_uppercase)
    {
        let should_lower = matches!(
            atom.atomic_number(),
            5 | 6 | 7 | 8 | 14 | 15 | 16 | 33 | 34 | 52
        );
        if should_lower {
            let mut owned = String::with_capacity(raw_symbol.len());
            let mut chars = raw_symbol.chars();
            if let Some(first) = chars.next() {
                owned.extend(first.to_lowercase());
            }
            owned.push_str(chars.as_str());
            lowered_symbol = owned;
            &lowered_symbol
        } else {
            raw_symbol
        }
    } else {
        raw_symbol
    };
    let mut result = String::new();
    if needs_bracket {
        result.push('[');
    }
    if let Some(isotope) = atom.isotope()
        && params.do_isomeric_smiles
    {
        result.push_str(&isotope.to_string());
    }
    result.push_str(symbol);
    result.push_str(&chirality);
    if needs_bracket {
        let total_num_hs = total_num_hydrogens_for_writer(molecule, atom_id);
        if total_num_hs > 0 {
            result.push('H');
            if total_num_hs > 1 {
                result.push_str(&total_num_hs.to_string());
            }
        }
        if atom.formal_charge() > 0 {
            result.push('+');
            if atom.formal_charge() > 1 {
                result.push_str(&atom.formal_charge().to_string());
            }
        } else if atom.formal_charge() < 0 {
            if atom.formal_charge() < -1 {
                result.push_str(&atom.formal_charge().to_string());
            } else {
                result.push('-');
            }
        }
        if let Some(atom_map) = atom.atom_map() {
            result.push(':');
            result.push_str(&atom_map.to_string());
        }
        result.push(']');
    }
    if let Some(label) = atom.prop("_supplementalSmilesLabel") {
        result.push_str(label);
    }
    Ok(result)
}

fn build_atom_smiles(
    molecule: &Molecule,
    atom_id: AtomId,
    params: &SmilesWriteParams,
    context: &SmilesWriteContext,
) -> Result<String, SmilesWriteError> {
    get_atom_smiles_with_context(molecule, atom_id, params, context)
}

pub fn get_bond_smiles(_bond_order: BondOrder) -> Result<&'static str, SmilesWriteError> {
    // RDKit✔️✔️: default: res = "~";
    match _bond_order {
        BondOrder::Single => Ok(""),
        BondOrder::Double => Ok("="),
        BondOrder::Triple => Ok("#"),
        BondOrder::Quadruple => Ok("$"),
        BondOrder::Dative => Ok("->"),
        _ => Ok("~"),
    }
}

pub fn get_molecule_bond_smiles(
    molecule: &Molecule,
    bond: usize,
    atom_to_left: Option<usize>,
    params: &SmilesWriteParams,
) -> Result<String, SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION GetBondSmiles
    // RDKit✔️✔️: std::string GetBondSmiles(const Bond *bond, const SmilesWriteParams &params,
    // RDKit✔️✔️:                           int atomToLeftIdx) {
    // RDKit✔️✔️:   PRECONDITION(bond, "bad bond");
    // RDKit✔️✔️:   if (atomToLeftIdx < 0) {
    // RDKit✔️✔️:     atomToLeftIdx = bond->getBeginAtomIdx();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::string res = "";
    // RDKit✔️✔️:   bool aromatic = false;
    // RDKit✔️✔️:   if (!params.doKekule && (bond->getBondType() == Bond::SINGLE ||
    // RDKit✔️✔️:                            bond->getBondType() == Bond::DOUBLE ||
    // RDKit✔️✔️:                            bond->getBondType() == Bond::AROMATIC)) {
    // RDKit✔️✔️:     aromatic = true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   Bond::BondDir dir = bond->getBondDir();
    // RDKit✔️✔️:   bond->clearProp(common_properties::_TraversalRingClosureBond);
    // RDKit✔️✔️:   switch (bond->getBondType()) {
    // RDKit✔️✔️:     case Bond::SINGLE:
    // RDKit✔️✔️:       if (dir != Bond::NONE && dir != Bond::UNKNOWN) {
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         if (params.allBondsExplicit) {
    // RDKit✔️✔️:           res = "-";
    // RDKit✔️✔️:         } else if (aromatic && !bond->getIsAromatic()) {
    // RDKit✔️✔️:           res = "-";
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::DOUBLE:
    // RDKit✔️✔️:       if (!aromatic || !bond->getIsAromatic() || params.allBondsExplicit) {
    // RDKit✔️✔️:         res = "=";
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::TRIPLE:
    // RDKit✔️✔️:       res = "#";
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::QUADRUPLE:
    // RDKit✔️✔️:       res = "$";
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::AROMATIC:
    // RDKit✔️✔️:       if (params.allBondsExplicit) {
    // RDKit✔️✔️:         res = ":";
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::DATIVE:
    // RDKit✔️✔️:       if (atomToLeftIdx >= 0 &&
    // RDKit✔️✔️:           bond->getBeginAtomIdx() == static_cast<unsigned int>(atomToLeftIdx)) {
    // RDKit✔️✔️:         res = "->";
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         res = "<-";
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       res = "~";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION GetBondSmiles
    validate_bond_index(molecule, bond)?;
    if let Some(atom) = atom_to_left {
        validate_atom_index(molecule, atom)?;
    }
    let bond = &molecule.bonds()[bond];
    let atom_to_left = atom_to_left.unwrap_or_else(|| bond.begin().index());
    let aromatic_context = if !params.do_kekule
        && matches!(
            bond.order(),
            BondOrder::Single | BondOrder::Double | BondOrder::Aromatic
        ) {
        let left = &molecule.atoms()[atom_to_left];
        let other_id = bond_other_atom(bond, AtomId::new(atom_to_left)).ok_or(
            SmilesWriteError::BondOutOfRange {
                bond: bond.id().index(),
            },
        )?;
        let other = &molecule.atoms()[other_id.index()];
        left.is_aromatic()
            && other.is_aromatic()
            && (left.atomic_number() != 0 || other.atomic_number() != 0)
    } else {
        false
    };
    match bond.order() {
        // RDKit✔️✔️: case Bond::SINGLE:
        // RDKit✔️✔️:   if (dir != Bond::NONE && dir != Bond::UNKNOWN) {
        // RDKit✔️✔️:     if (dir == Bond::BEGINWEDGE || dir == Bond::BEGINDASH ||
        // RDKit✔️✔️:         dir == Bond::ENDDOWNRIGHT || dir == Bond::ENDUPRIGHT) {
        // RDKit✔️✔️:       res = dirSymbol(dir, atomToLeftIdx);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   } else if (params.allBondsExplicit) { res = "-"; }
        BondOrder::Single => {
            if !matches!(
                bond.direction(),
                BondDirection::None | BondDirection::Unknown
            ) {
                match bond.direction() {
                    BondDirection::EndDownRight => {
                        if params.all_bonds_explicit || params.do_isomeric_smiles {
                            Ok("\\".to_string())
                        } else {
                            Ok(String::new())
                        }
                    }
                    BondDirection::EndUpRight => {
                        if params.all_bonds_explicit || params.do_isomeric_smiles {
                            Ok("/".to_string())
                        } else {
                            Ok(String::new())
                        }
                    }
                    _ => {
                        if params.all_bonds_explicit {
                            Ok("-".to_string())
                        } else {
                            Ok(String::new())
                        }
                    }
                }
            } else if params.all_bonds_explicit || (aromatic_context && !bond.is_aromatic()) {
                Ok("-".to_string())
            } else {
                Ok(String::new())
            }
        }
        // RDKit✔️✔️: case Bond::DOUBLE:
        // RDKit✔️✔️:   if (!aromatic || !bond->getIsAromatic() || params.allBondsExplicit) {
        // RDKit✔️✔️:     res = "=";
        // RDKit✔️✔️:   }
        BondOrder::Double => {
            if !aromatic_context || !bond.is_aromatic() || params.all_bonds_explicit {
                Ok("=".to_string())
            } else {
                Ok(String::new())
            }
        }
        // RDKit✔️✔️: case Bond::TRIPLE: res = "#"; break;
        BondOrder::Triple => Ok("#".to_string()),
        // RDKit✔️✔️: case Bond::QUADRUPLE: res = "$"; break;
        BondOrder::Quadruple => Ok("$".to_string()),
        // RDKit✔️✔️: case Bond::AROMATIC:
        // RDKit✔️✔️:   if (params.allBondsExplicit) { res = ":"; }
        // RDKit✔️✔️:   break;
        BondOrder::Aromatic => {
            if !matches!(
                bond.direction(),
                BondDirection::None | BondDirection::Unknown
            ) {
                match bond.direction() {
                    BondDirection::EndDownRight => {
                        if params.all_bonds_explicit || params.do_isomeric_smiles {
                            Ok("\\".to_string())
                        } else {
                            Ok(String::new())
                        }
                    }
                    BondDirection::EndUpRight => {
                        if params.all_bonds_explicit || params.do_isomeric_smiles {
                            Ok("/".to_string())
                        } else {
                            Ok(String::new())
                        }
                    }
                    _ => {
                        if params.all_bonds_explicit || !aromatic_context {
                            Ok(":".to_string())
                        } else {
                            Ok(String::new())
                        }
                    }
                }
            } else if params.all_bonds_explicit || !aromatic_context {
                Ok(":".to_string())
            } else {
                Ok(String::new())
            }
        }
        // RDKit✔️✔️: case Bond::DATIVE:
        BondOrder::Dative => {
            if bond.begin().index() == atom_to_left {
                Ok("->".to_string())
            } else {
                Ok("<-".to_string())
            }
        }
        // RDKit✔️✔️: default: res = "~";
        _ => Ok("~".to_string()),
    }
}

fn build_bond_smiles(
    molecule: &Molecule,
    bond: BondId,
    atom_to_left: AtomId,
    params: &SmilesWriteParams,
) -> Result<String, SmilesWriteError> {
    get_molecule_bond_smiles(molecule, bond.index(), Some(atom_to_left.index()), params)
}

fn total_num_hydrogens_for_writer(molecule: &Molecule, atom_id: AtomId) -> u32 {
    let explicit = u32::from(molecule.atoms()[atom_id.index()].explicit_hydrogens());
    let implicit = molecule
        .derived_cache()
        .valence
        .as_ref()
        .and_then(|valence| valence.implicit_hydrogens.get(atom_id.index()))
        .copied()
        .unwrap_or(0)
        .max(0) as u32;
    explicit + implicit
}

fn total_valence_for_writer(molecule: &Molecule, atom_id: AtomId) -> Option<i32> {
    molecule.derived_cache().valence.as_ref().map(|valence| {
        valence.explicit_valence[atom_id.index()] + valence.implicit_hydrogens[atom_id.index()]
    })
}

#[must_use]
pub fn in_organic_subset(_atomic_number: u8) -> Result<bool, SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION inOrganicSubset
    // RDKit✔️✔️: const int atomicSmiles[] = {0, 5, 6, 7, 8, 9, 15, 16, 17, 35, 53, -1};
    // RDKit✔️✔️: bool inOrganicSubset(int atomicNumber) {
    // RDKit✔️✔️:   unsigned int idx = 0;
    // RDKit✔️✔️:   while (atomicSmiles[idx] < atomicNumber && atomicSmiles[idx] != -1) {
    // RDKit✔️✔️:     ++idx;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return atomicSmiles[idx] == atomicNumber;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION inOrganicSubset
    Ok(matches!(
        _atomic_number,
        0 | 5 | 6 | 7 | 8 | 9 | 15 | 16 | 17 | 35 | 53
    ))
}

// BEGIN RDKIT CPP FUNCTION SmilesWrite::getCXExtensions
// RDKit✔️✔️: std::string getCXExtensions(const ROMol &mol, std::uint32_t flags) {
// RDKit✔️✔️:   std::string res = "|";
// RDKit✔️✔️:   const std::vector<unsigned int> &atomOrder =
// RDKit✔️✔️:       mol.getProp<std::vector<unsigned int>>(
// RDKit✔️✔️:           common_properties::_smilesAtomOutputOrder);
// RDKit✔️✔️:   if ((flags & SmilesWrite::CXSmilesFields::CX_COORDS) &&
// RDKit✔️✔️:       mol.getNumConformers()) {
// RDKit✔️✔️:     res += "(" + get_coords_block(mol, atomOrder) + ")";
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if ((flags & SmilesWrite::CXSmilesFields::CX_ATOM_LABELS) && needLabels) {
// RDKit✔️✔️:     auto lbls = get_atomlabel_block(mol, atomOrder);
// RDKit✔️✔️:     if (!lbls.empty()) {
// RDKit✔️✔️:       if (res.size() > 1) { res += ","; }
// RDKit✔️✔️:       res += "$" + lbls + "$";
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if ((flags & SmilesWrite::CXSmilesFields::CX_MOLFILE_VALUES) && needValues) {
// RDKit✔️✔️:     if (res.size() > 1) { res += ","; }
// RDKit✔️✔️:     res += "$_AV:" + get_value_block(...) + "$";
// RDKit✔️✔️:   }
// RDKit✔️✔️:   auto radblock = get_radical_block(mol, atomOrder);
// RDKit✔️✔️:   if ((flags & CX_RADICALS) && radblock.size()) {
// RDKit✔️✔️:     res += radblock;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (flags & CX_ATOM_PROPS) {
// RDKit✔️✔️:     appendToCXExtension(get_atom_props_block(mol, atomOrder), res);
// RDKit✔️✔️:   }
// RDKit✔️✔️:   // ... enhanced stereo, SGroups, bonds blocks follow same pattern
// RDKit✔️✔️:   if (res.size() > 1) { res += "|"; } else { res = ""; }
// RDKit✔️✔️:   return res;
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION SmilesWrite::getCXExtensions
pub fn get_cx_extensions(
    molecule: &Molecule,
    fields: CxSmilesFields,
) -> Result<String, SmilesWriteError> {
    get_cx_extensions_scoped(molecule, fields, &CxWriteScope::full_molecule(molecule))
}

fn get_cx_extensions_scoped(
    molecule: &Molecule,
    fields: CxSmilesFields,
    scope: &CxWriteScope,
) -> Result<String, SmilesWriteError> {
    write_cx_smiles_fields(molecule, fields, scope)
}

fn write_cx_smiles_fields(
    molecule: &Molecule,
    fields: CxSmilesFields,
    scope: &CxWriteScope,
) -> Result<String, SmilesWriteError> {
    let mut res = String::from("|");
    let append_to_cx = |addition: &str, buf: &mut String| {
        if !addition.is_empty() {
            if buf.len() > 1 {
                buf.push(',');
            }
            buf.push_str(addition);
        }
    };

    if fields.contains(CxSmilesFields::COORDS) {
        let coords = write_cx_coordinates(molecule, &scope.atom_order);
        if !coords.is_empty() {
            res.push('(');
            res.push_str(&coords);
            res.push(')');
        }
    }

    let need_labels = scope.atom_order.iter().any(|atom_id| {
        let atom = &molecule.atoms()[atom_id.index()];
        atom.prop("atomLabel").is_some()
            || atom.prop("_QueryAtomGenericLabel").is_some()
            || atom.prop("dummyLabel").is_some()
            || atom.prop("_fromAttachPoint").is_some()
    });
    if fields.contains(CxSmilesFields::ATOM_LABELS) && need_labels {
        let labels = write_cx_atom_labels(molecule, &scope.atom_order);
        if !labels.is_empty() {
            append_to_cx(&format!("${}$", labels), &mut res);
        }
    }

    let need_values = scope.atom_order.iter().any(|atom_id| {
        molecule.atoms()[atom_id.index()]
            .prop("molFileValue")
            .is_some()
    });
    if fields.contains(CxSmilesFields::MOLFILE_VALUES) && need_values {
        let values = write_cx_atom_values(molecule, &scope.atom_order);
        if !values.is_empty() {
            append_to_cx(&format!("$_AV:{}$", values), &mut res);
        }
    }

    if fields.contains(CxSmilesFields::RADICALS) {
        let radicals = write_cx_radicals(molecule, &scope.atom_order);
        if !radicals.is_empty() {
            if res.len() > 1 {
                res.push(',');
            }
            res.push_str(&radicals);
            if res.ends_with(',') {
                res.pop();
            }
        }
    }

    if fields.contains(CxSmilesFields::ATOM_PROPS) {
        let props = write_cx_atom_props(molecule, &scope.atom_order);
        append_to_cx(&props, &mut res);
    }

    if fields.contains(CxSmilesFields::BOND_CFG) {
        let include_coords =
            fields.contains(CxSmilesFields::COORDS) && molecule.coords_2d().is_some();
        let bond_cfg = write_cx_bond_config_block(
            molecule,
            &scope.atom_order,
            &scope.bond_order,
            include_coords,
            false,
        );
        append_to_cx(&bond_cfg, &mut res);
        let ringbond_cistrans =
            write_cx_ringbond_cistrans_block(molecule, &scope.atom_order, &scope.bond_order);
        append_to_cx(&ringbond_cistrans, &mut res);
    } else if fields.contains(CxSmilesFields::BOND_ATROPISOMER) {
        let include_coords =
            fields.contains(CxSmilesFields::COORDS) && molecule.coords_2d().is_some();
        let bond_cfg = write_cx_bond_config_block(
            molecule,
            &scope.atom_order,
            &scope.bond_order,
            include_coords,
            true,
        );
        append_to_cx(&bond_cfg, &mut res);
    }

    if fields.contains(CxSmilesFields::ENHANCED_STEREO) {
        let stereo = write_cx_enhanced_stereo(molecule, &scope.atom_order, &scope.bond_order);
        append_to_cx(&stereo, &mut res);
    }

    if fields.contains(CxSmilesFields::SGROUPS) {
        let sgroups = write_cx_sgroups(molecule, &scope.atom_order, &scope.bond_order);
        append_to_cx(&sgroups, &mut res);
    }

    if fields.contains(CxSmilesFields::POLYMER) {
        let polymer = write_cx_polymer_sgroups(molecule, &scope.atom_order, &scope.bond_order);
        append_to_cx(&polymer, &mut res);
    }

    if fields.contains(CxSmilesFields::SGROUPS) || fields.contains(CxSmilesFields::POLYMER) {
        let hierarchy = write_cx_sgroup_hierarchy_block(
            molecule,
            &scope.atom_order,
            &scope.bond_order,
            fields.contains(CxSmilesFields::SGROUPS),
            fields.contains(CxSmilesFields::POLYMER),
        );
        append_to_cx(&hierarchy, &mut res);
    }

    if fields.contains(CxSmilesFields::COORDINATE_BONDS) {
        let coord_bonds =
            write_cx_coordinate_bonds(molecule, &scope.atom_order, &scope.bond_order, "C");
        append_to_cx(&coord_bonds, &mut res);
    }

    if fields.contains(CxSmilesFields::HYDROGEN_BONDS) {
        let h_bonds =
            write_cx_coordinate_bonds(molecule, &scope.atom_order, &scope.bond_order, "H");
        append_to_cx(&h_bonds, &mut res);
    }

    if fields.contains(CxSmilesFields::ZERO_BONDS) {
        let zero_bonds = write_cx_zero_bonds(molecule, &scope.bond_order);
        append_to_cx(&zero_bonds, &mut res);
    }

    if fields.contains(CxSmilesFields::LINKNODES) {
        let linknodes = write_cx_linknodes_block(molecule, &scope.atom_order);
        append_to_cx(&linknodes, &mut res);
    }

    // RDKit: if (res.size() > 1) { res += "|"; } else { res = ""; }
    if res.len() > 1 {
        res.push('|');
    } else {
        res.clear();
    }
    Ok(res)
}

fn cx_atom_output_positions(atom_order: &[AtomId], atom_count: usize) -> Vec<Option<usize>> {
    let mut positions = vec![None; atom_count];
    for (position, atom_id) in atom_order.iter().copied().enumerate() {
        if atom_id.index() < positions.len() {
            positions[atom_id.index()] = Some(position);
        }
    }
    positions
}

fn cx_bond_output_positions(bond_order: &[BondId], bond_count: usize) -> Vec<Option<usize>> {
    let mut positions = vec![None; bond_count];
    for (position, bond_id) in bond_order.iter().copied().enumerate() {
        if bond_id.index() < positions.len() {
            positions[bond_id.index()] = Some(position);
        }
    }
    positions
}

fn zero_small_writer_coord(value: f64) -> f64 {
    if value.abs() < 1e-4 { 0.0 } else { value }
}

fn quote_atomprop_string(text: &str) -> String {
    text.chars()
        .map(|ch| {
            if ch == '.' {
                "&#46;".to_string()
            } else {
                ch.to_string()
            }
        })
        .collect()
}

fn write_cx_coordinates(molecule: &Molecule, atom_order: &[AtomId]) -> String {
    let coords = match molecule.coords_2d() {
        Some(c) => c,
        None => return String::new(),
    };
    let mut parts = Vec::new();
    for atom in atom_order {
        let Some(coord) = coords.get(atom.index()) else {
            continue;
        };
        parts.push(format!(
            "{}, {},",
            zero_small_writer_coord(coord[0]),
            zero_small_writer_coord(coord[1])
        ));
    }
    parts
        .into_iter()
        .map(|part| part.replace(", ", ","))
        .collect::<Vec<_>>()
        .join(";")
}

fn write_cx_coords(molecule: &Molecule, atom_order: &[AtomId]) -> String {
    write_cx_coordinates(molecule, atom_order)
}

fn write_cx_atom_labels(molecule: &Molecule, atom_order: &[AtomId]) -> String {
    let pseudoatoms = ["Pol", "Mod", "Het", "Any", "A", "Q", "X", "*"];
    let mut parts = Vec::new();
    for atom_id in atom_order {
        let atom = &molecule.atoms()[atom_id.index()];
        let part = if let Some(label) = atom.prop("_QueryAtomGenericLabel") {
            Some(format!("{label}_p"))
        } else if atom.atomic_number() == 0
            && atom
                .prop("dummyLabel")
                .is_some_and(|label| pseudoatoms.contains(&label))
        {
            Some(format!("{}_p", atom.prop("dummyLabel").unwrap_or_default()))
        } else if atom.atomic_number() == 0
            && atom
                .prop("_fromAttachPoint")
                .and_then(|value| value.parse::<u32>().ok())
                .is_some_and(|value| value == 1 || value == 2)
        {
            Some(format!(
                "_AP{}",
                atom.prop("_fromAttachPoint").unwrap_or_default()
            ))
        } else {
            atom.prop("atomLabel").map(str::to_string)
        };
        if let Some(part) = part {
            parts.push(part);
        } else {
            parts.push(String::new());
        }
    }
    if parts.iter().all(|part| part.is_empty()) {
        String::new()
    } else {
        parts.join(";")
    }
}

fn write_cx_atom_values(molecule: &Molecule, atom_order: &[AtomId]) -> String {
    let mut parts = Vec::new();
    for atom_id in atom_order {
        let atom = &molecule.atoms()[atom_id.index()];
        if let Some(value) = atom.prop("molFileValue") {
            parts.push(value.to_string());
        } else {
            parts.push(String::new());
        }
    }
    parts.join(";")
}

fn write_cx_molfile_values(molecule: &Molecule, atom_order: &[AtomId]) -> String {
    write_cx_atom_values(molecule, atom_order)
}

fn write_cx_radicals(molecule: &Molecule, atom_order: &[AtomId]) -> String {
    let mut by_count: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
    for (output_idx, atom_id) in atom_order.iter().copied().enumerate() {
        let atom = &molecule.atoms()[atom_id.index()];
        let re = atom.radical_electrons();
        if re > 0 {
            by_count.entry(re).or_default().push(output_idx);
        }
    }
    if by_count.is_empty() {
        return String::new();
    }
    let mut result = String::new();
    for (count, atoms) in by_count {
        let code = match count {
            1 => "^1:",
            2 => "^2:",
            3 => "^5:",
            _ => continue,
        };
        result.push_str(code);
        for atom in atoms {
            result.push_str(&format!("{atom},"));
        }
    }
    result
}

fn write_cx_atom_props(molecule: &Molecule, atom_order: &[AtomId]) -> String {
    let skip = [
        "atomLabel",
        "molFileValue",
        "molParity",
        "molStereoCare",
        "molRxnExactChange",
        "molInversionFlag",
        "dummyLabel",
    ];
    let mut result = String::new();
    for (which, atom_id) in atom_order.iter().copied().enumerate() {
        let atom = &molecule.atoms()[atom_id.index()];
        let is_attachment_point =
            atom.atomic_number() == 0 && atom.prop("_fromAttachPoint").is_some();
        for (prop_name, prop_value) in atom.props() {
            // RDKit getPropList(includePrivate=false, includeComputed=false)
            // excludes underscore-prefixed writer-internal/cache props from
            // CX atomProp emission.
            if prop_name.starts_with('_') {
                continue;
            }
            if skip.contains(&prop_name.as_str()) || prop_name == "molAtomMapNumber" {
                continue;
            }
            if prop_name == "dummyLabel"
                && (is_attachment_point
                    || prop_value == "*"
                    || ["Pol", "Mod", "Het", "Any", "A", "Q", "X", "*"]
                        .contains(&prop_value.as_str()))
            {
                continue;
            }
            if result.is_empty() {
                result.push_str("atomProp");
            }
            result.push_str(&format!(
                ":{which}.{}.{}",
                quote_atomprop_string(prop_name),
                quote_atomprop_string(prop_value)
            ));
        }
    }
    result
}

fn write_cx_enhanced_stereo(
    molecule: &Molecule,
    atom_order: &[AtomId],
    _bond_order: &[BondId],
) -> String {
    use crate::stereo::StereoGroupKind;
    let atom_positions = cx_atom_output_positions(atom_order, molecule.num_atoms());
    let write_ids = assigned_writer_stereo_group_ids(molecule.stereo_groups());
    let mut parts: Vec<String> = Vec::new();
    for (group, write_id) in molecule.stereo_groups().iter().zip(write_ids) {
        let mut atom_idxs: Vec<usize> = group
            .atoms()
            .iter()
            .filter_map(|atom| {
                atom_positions
                    .get(atom.index())
                    .and_then(|position| *position)
            })
            .collect();
        if atom_idxs.is_empty() {
            continue;
        }
        atom_idxs.sort_unstable();
        let prefix = match group.kind() {
            StereoGroupKind::Absolute => "a".to_string(),
            StereoGroupKind::Or => format!("o{}", write_id.unwrap_or(1)),
            StereoGroupKind::And => format!("&{}", write_id.unwrap_or(1)),
        };
        let members = atom_idxs
            .into_iter()
            .map(|idx| idx.to_string())
            .collect::<Vec<_>>();
        parts.push(format!("{prefix}:{}", members.join(",")));
    }
    parts.join(",")
}

fn assigned_writer_stereo_group_ids(groups: &[crate::StereoGroup]) -> Vec<Option<u32>> {
    use crate::stereo::StereoGroupKind;

    let mut or_ids = Vec::<u32>::new();
    let mut and_ids = Vec::<u32>::new();
    let mut assigned = groups
        .iter()
        .map(crate::StereoGroup::id)
        .collect::<Vec<_>>();
    for (idx, group) in groups.iter().enumerate() {
        let Some(id) = assigned[idx] else {
            continue;
        };
        let ids = match group.kind() {
            StereoGroupKind::Or => &mut or_ids,
            StereoGroupKind::And => &mut and_ids,
            StereoGroupKind::Absolute => continue,
        };
        if id != 0 && ids.contains(&id) {
            assigned[idx] = None;
        } else if id != 0 {
            ids.push(id);
        }
    }

    let mut next_or = 0_u32;
    let mut next_and = 0_u32;
    for (idx, group) in groups.iter().enumerate() {
        if group.kind() == StereoGroupKind::Absolute || assigned[idx].is_some() {
            continue;
        }
        let (next, ids) = match group.kind() {
            StereoGroupKind::Or => (&mut next_or, &mut or_ids),
            StereoGroupKind::And => (&mut next_and, &mut and_ids),
            StereoGroupKind::Absolute => unreachable!(),
        };
        *next += 1;
        while ids.contains(&*next) {
            *next += 1;
        }
        ids.push(*next);
        assigned[idx] = Some(*next);
    }
    assigned
}

fn write_cx_sgroups(molecule: &Molecule, atom_order: &[AtomId], bond_order: &[BondId]) -> String {
    let data = write_cx_data_sgroups(molecule, atom_order);
    let other = write_cx_non_data_sgroups(molecule, atom_order, bond_order);
    match (data.is_empty(), other.is_empty()) {
        (true, true) => String::new(),
        (false, true) => data,
        (true, false) => other,
        (false, false) => format!("{data},{other}"),
    }
}

fn write_cx_data_sgroups(molecule: &Molecule, atom_order: &[AtomId]) -> String {
    let atom_positions = cx_atom_output_positions(atom_order, molecule.num_atoms());
    let mut parts = Vec::new();
    for sgroup in molecule.substance_groups() {
        if !writer_is_data_sgroup(sgroup) {
            continue;
        }
        let atoms = sgroup
            .atoms()
            .iter()
            .filter_map(|atom| atom_positions.get(atom.index()).and_then(|value| *value))
            .map(|idx| idx.to_string())
            .collect::<Vec<_>>();
        if atoms.is_empty() {
            continue;
        }
        let field_name = writer_data_sgroup_field_name(sgroup);
        let data_fields = writer_data_sgroup_values(sgroup).join(",");
        let query_op = writer_data_sgroup_query_op(sgroup);
        let field_info = writer_data_sgroup_field_info(sgroup);
        let field_tag = writer_data_sgroup_field_tag(sgroup);
        parts.push(format!(
            "SgD:{}:{field_name}:{data_fields}:{query_op}:{field_info}:{field_tag}:",
            atoms.join(",")
        ));
    }
    parts.join(",")
}

fn write_cx_non_data_sgroups(
    molecule: &Molecule,
    atom_order: &[AtomId],
    bond_order: &[BondId],
) -> String {
    use crate::sgroup::SubstanceGroupKind;
    let atom_positions = cx_atom_output_positions(atom_order, molecule.num_atoms());
    let bond_positions = cx_bond_output_positions(bond_order, molecule.num_bonds());
    let mut parts = Vec::new();
    for sgroup in molecule.substance_groups() {
        if matches!(sgroup.kind(), SubstanceGroupKind::Data)
            || sgroup
                .props()
                .get("TYPE")
                .is_some_and(|value| value == "DAT")
            || writer_polymer_sgroup_type_code(sgroup).is_some()
        {
            continue;
        }
        let kind_str = match sgroup.kind() {
            SubstanceGroupKind::Data => "DAT",
            SubstanceGroupKind::Superatom => "SUP",
            SubstanceGroupKind::MultipleGroup => "MUL",
            SubstanceGroupKind::StructuralRepeatUnit => "SRU",
            SubstanceGroupKind::Monomer => "MON",
            SubstanceGroupKind::Copolymer => "COP",
            SubstanceGroupKind::Crosslink => "CRO",
            SubstanceGroupKind::Graft => "GRA",
            SubstanceGroupKind::Modification => "MOD",
            SubstanceGroupKind::Mer => "MER",
            SubstanceGroupKind::AnyPolymer => "ANY",
            SubstanceGroupKind::MixtureComponent => "MIX",
            SubstanceGroupKind::Mixture => "MIXTURE",
            SubstanceGroupKind::Formulation => "FOR",
            SubstanceGroupKind::Generic(s) => s.as_str(),
        };
        let atom_idxs: Vec<String> = sgroup
            .atoms()
            .iter()
            .filter_map(|atom| atom_positions.get(atom.index()).and_then(|value| *value))
            .map(|idx| idx.to_string())
            .collect();
        let bond_idxs: Vec<String> = sgroup
            .bonds()
            .iter()
            .filter_map(|bond| bond_positions.get(bond.index()).and_then(|value| *value))
            .map(|idx| idx.to_string())
            .collect();
        if atom_idxs.is_empty() && bond_idxs.is_empty() {
            continue;
        }
        let mut entry = format!(
            "_S:{}:{}:{}",
            kind_str,
            atom_idxs.join(","),
            bond_idxs.join(",")
        );
        if let Some(label) = sgroup.label() {
            entry.push(':');
            entry.push_str(&label.replace(',', "\\,").replace('|', "\\|"));
        }
        if let Some(conn) = sgroup.connection() {
            let conn_str = match conn {
                crate::sgroup::SGroupConnection::HeadToHead => "HH",
                crate::sgroup::SGroupConnection::HeadToTail => "HT",
                crate::sgroup::SGroupConnection::Either => "EU",
                crate::sgroup::SGroupConnection::Unknown(s) => s,
            };
            entry.push(':');
            entry.push_str(conn_str);
        }
        parts.push(entry);
    }
    parts.join(",")
}

fn write_cx_coordinate_bonds(
    molecule: &Molecule,
    atom_order: &[AtomId],
    bond_order: &[BondId],
    symbol: &str,
) -> String {
    let atom_positions = cx_atom_output_positions(atom_order, molecule.num_atoms());
    let target_order = match symbol {
        "C" => BondOrder::Dative,
        "H" => BondOrder::Hydrogen,
        _ => return String::new(),
    };
    let mut parts = Vec::new();
    for (bond_output_idx, bond_id) in bond_order.iter().copied().enumerate() {
        let bond = &molecule.bonds()[bond_id.index()];
        let matches = if symbol == "C" {
            matches!(bond.order(), BondOrder::Dative | BondOrder::DativeOne)
        } else {
            bond.order() == target_order
        };
        if !matches {
            continue;
        }
        let Some(begin_output_idx) = atom_positions
            .get(bond.begin().index())
            .and_then(|value| *value)
        else {
            continue;
        };
        parts.push(format!("{begin_output_idx}.{bond_output_idx}"));
    }
    if parts.is_empty() {
        String::new()
    } else {
        format!("{symbol}:{}", parts.join(","))
    }
}

fn write_cx_zero_bonds(molecule: &Molecule, bond_order: &[BondId]) -> String {
    let mut parts = Vec::new();
    for (bond_output_idx, bond_id) in bond_order.iter().copied().enumerate() {
        if molecule.bonds()[bond_id.index()].order() == BondOrder::Zero {
            parts.push(bond_output_idx.to_string());
        }
    }
    if parts.is_empty() {
        String::new()
    } else {
        format!("Z:{}", parts.join(","))
    }
}

// BEGIN RDKIT CPP FUNCTION get_bond_config_block
// RDKit✔️✔️: std::string get_bond_config_block(
// RDKit✔️✔️:     const ROMol &mol, const std::vector<unsigned int> &atomOrder,
// RDKit✔️✔️:     const std::vector<unsigned int> &bondOrder, bool coordsIncluded,
// RDKit✔️✔️:     std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>> &wedgeBonds,
// RDKit✔️✔️:     bool atropisomerOnly = false) {
// RDKit✔️✔️:   std::map<std::string, std ::vector<std::string>> wParts;
// RDKit✔️✔️:   for (unsigned int i = 0; i < bondOrder.size(); ++i) {
// RDKit✔️✔️:     auto idx = bondOrder[i];
// RDKit✔️✔️:     const auto bond = mol.getBondWithIdx(idx);
// RDKit✔️✔️:     unsigned int wedgeStartAtomIdx = bond->getBeginAtomIdx();
// RDKit✔️✔️:     if (!canHaveDirection(*bond)) { continue; }
// RDKit✔️✔️:     Bond::BondDir bd = bond->getBondDir();
// RDKit✔️✔️:     switch (bd) {
// RDKit✔️✔️:       case Bond::BondDir::BEGINDASH:
// RDKit✔️✔️:       case Bond::BondDir::BEGINWEDGE:
// RDKit✔️✔️:       case Bond::BondDir::UNKNOWN:
// RDKit✔️✔️:         break;
// RDKit✔️✔️:       default:
// RDKit✔️✔️:         bd = Bond::BondDir::NONE;
// RDKit✔️✔️:     }
// RDKit✔️✔️:     if (atropisomerOnly && bd == Bond::BondDir::NONE) { continue; }
// RDKit✔️✔️:     // ... atropisomer checks and optional wedge flipping ...
// RDKit✔️✔️:     if (!atropisomerOnly) {
// RDKit✔️✔️:       unsigned int cfg = 0;
// RDKit✔️✔️:       if (bd == Bond::BondDir::NONE &&
// RDKit✔️✔️:           bond->getPropIfPresent(common_properties::_MolFileBondCfg, cfg)) {
// RDKit✔️✔️:         switch (cfg) {
// RDKit✔️✔️:           case 1: bd = Bond::BondDir::BEGINWEDGE; break;
// RDKit✔️✔️:           case 2: bd = Bond::BondDir::UNKNOWN; break;
// RDKit✔️✔️:           case 3: bd = Bond::BondDir::BEGINDASH; break;
// RDKit✔️✔️:           default: bd = Bond::BondDir::NONE;
// RDKit✔️✔️:         }
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:     auto begAtomOrder =
// RDKit✔️✔️:         std::find(atomOrder.begin(), atomOrder.end(), wedgeStartAtomIdx) -
// RDKit✔️✔️:         atomOrder.begin();
// RDKit✔️✔️:     std::string wType = "";
// RDKit✔️✔️:     if (bd == Bond::BondDir::UNKNOWN) { wType = "w"; }
// RDKit✔️✔️:     else if (coordsIncluded || isAnAtropisomer) {
// RDKit✔️✔️:       if (bd == Bond::BondDir::BEGINWEDGE) { wType = "wU"; }
// RDKit✔️✔️:       else if (bd == Bond::BondDir::BEGINDASH) { wType = "wD"; }
// RDKit✔️✔️:     }
// RDKit✔️✔️:     if (wType != "") { wParts[wType].push_back(format("%d.%d", begAtomOrder, i)); }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   // join as "w:...,wD:...,wU:..."
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION get_bond_config_block
fn write_cx_bond_config_block(
    molecule: &Molecule,
    atom_order: &[AtomId],
    bond_order: &[BondId],
    coords_included: bool,
    atropisomer_only: bool,
) -> String {
    let mut atom_order_positions: Vec<Option<usize>> = vec![None; molecule.atoms().len()];
    for (position, atom_id) in atom_order.iter().copied().enumerate() {
        if atom_id.index() < atom_order_positions.len() {
            atom_order_positions[atom_id.index()] = Some(position);
        }
    }

    let mut w_parts: BTreeMap<&'static str, Vec<String>> = BTreeMap::new();
    for (bond_output_idx, bond_id) in bond_order.iter().copied().enumerate() {
        let bond = &molecule.bonds()[bond_id.index()];
        let wedge_start_atom = bond.begin();
        if !can_have_direction_for_writer(bond.order()) {
            continue;
        }

        let mut direction = normalize_writer_cx_bond_direction(bond.direction());
        let mut is_an_atropisomer = false;

        if atropisomer_only && direction == BondDirection::None {
            continue;
        }

        if matches!(
            direction,
            BondDirection::BeginDash | BondDirection::BeginWedge
        ) {
            for neighbor_bond_id in incident_bonds(molecule, wedge_start_atom) {
                if neighbor_bond_id == bond_id {
                    continue;
                }
                let neighbor = &molecule.bonds()[neighbor_bond_id.index()];
                if matches!(
                    neighbor.stereo(),
                    BondStereo::AtropCw | BondStereo::AtropCcw
                ) {
                    is_an_atropisomer = true;
                    break;
                }
            }
        }

        if atropisomer_only {
            if !is_an_atropisomer {
                continue;
            }
        } else if matches!(direction, BondDirection::None)
            && let Some(cfg) = writer_parse_molfile_bond_cfg(bond)
        {
            direction = match cfg {
                1 => BondDirection::BeginWedge,
                2 => BondDirection::Unknown,
                3 => BondDirection::BeginDash,
                _ => BondDirection::None,
            };
        }

        // RDKit also derives direction from coordinates via
        // Chirality::GetMolFileBondStereoInfo when available; that path is
        // deferred until the full wedgeBonds/conformer parity port.
        let w_type = if direction == BondDirection::Unknown {
            Some("w")
        } else if coords_included || is_an_atropisomer {
            match direction {
                BondDirection::BeginWedge => Some("wU"),
                BondDirection::BeginDash => Some("wD"),
                _ => None,
            }
        } else {
            None
        };

        let Some(w_type) = w_type else {
            continue;
        };

        let Some(Some(begin_atom_order_idx)) = atom_order_positions.get(wedge_start_atom.index())
        else {
            continue;
        };
        w_parts
            .entry(w_type)
            .or_default()
            .push(format!("{begin_atom_order_idx}.{bond_output_idx}"));
    }

    let mut parts: Vec<String> = Vec::new();
    for (w_type, entries) in w_parts {
        if entries.is_empty() {
            continue;
        }
        parts.push(format!("{w_type}:{}", entries.join(",")));
    }
    parts.join(",")
}

fn normalize_writer_cx_bond_direction(direction: BondDirection) -> BondDirection {
    match direction {
        BondDirection::BeginDash | BondDirection::BeginWedge | BondDirection::Unknown => direction,
        _ => BondDirection::None,
    }
}

fn writer_parse_molfile_bond_cfg(bond: &Bond) -> Option<u32> {
    bond.prop("_MolFileBondCfg")
        .and_then(|value| value.parse::<u32>().ok())
}

// BEGIN RDKIT CPP FUNCTION get_ringbond_cistrans_block
// RDKit✔️✔️: std::string get_ringbond_cistrans_block(
// RDKit✔️✔️:     const ROMol &mol, const std::vector<unsigned int> &atomOrder,
// RDKit✔️✔️:     const std::vector<unsigned int> &bondOrder) {
// RDKit✔️✔️:   if (!mol.getRingInfo()->isInitialized()) { return ""; }
// RDKit✔️✔️:   std::string c = "", t = "", ctu = "";
// RDKit✔️✔️:   for (unsigned int i = 0; i < bondOrder.size(); ++i) {
// RDKit✔️✔️:     auto idx = bondOrder[i];
// RDKit✔️✔️:     if (!rinfo->numBondRings(idx) ||
// RDKit✔️✔️:         rinfo->minBondRingSize(idx) < Chirality::minRingSizeForDoubleBondStereo) { continue; }
// RDKit✔️✔️:     if (bond->getBondType() != Bond::DOUBLE && bond->getBondType() != Bond::AROMATIC) { continue; }
// RDKit✔️✔️:     if (bstereo != Bond::STEREOANY && bstereo != Bond::STEREOCIS && bstereo != Bond::STEREOTRANS) { continue; }
// RDKit✔️✔️:     if (bstereo == Bond::STEREOANY) { ctu += label; } else {
// RDKit✔️✔️:       bool needSwap = false;
// RDKit✔️✔️:       // parity flips from atom-order reordering around stereo refs
// RDKit✔️✔️:       if (bstereo == Bond::STEREOCIS || needSwap) { c += label; } else { t += label; }
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return c + t + ctu;
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION get_ringbond_cistrans_block
fn write_cx_ringbond_cistrans_block(
    molecule: &Molecule,
    atom_order: &[AtomId],
    bond_order: &[BondId],
) -> String {
    const MIN_RING_SIZE_FOR_DOUBLE_BOND_STEREO: usize = 8;
    let Some(ring_info) = molecule.derived_cache().rings.as_ref() else {
        return String::new();
    };
    if !ring_info.is_initialized() {
        return String::new();
    }

    let mut atom_order_positions: Vec<Option<usize>> = vec![None; molecule.atoms().len()];
    for (position, atom_id) in atom_order.iter().copied().enumerate() {
        if atom_id.index() < atom_order_positions.len() {
            atom_order_positions[atom_id.index()] = Some(position);
        }
    }

    let mut c_labels: Vec<String> = Vec::new();
    let mut t_labels: Vec<String> = Vec::new();
    let mut ctu_labels: Vec<String> = Vec::new();

    for (bond_output_idx, bond_id) in bond_order.iter().copied().enumerate() {
        if ring_info.num_bond_rings(bond_id) == 0
            || ring_info.min_bond_ring_size(bond_id) < MIN_RING_SIZE_FOR_DOUBLE_BOND_STEREO
        {
            continue;
        }
        let bond = &molecule.bonds()[bond_id.index()];
        if !matches!(bond.order(), BondOrder::Double | BondOrder::Aromatic) {
            continue;
        }
        if !matches!(
            bond.stereo(),
            BondStereo::Any | BondStereo::Cis | BondStereo::Trans
        ) {
            continue;
        }

        let label = bond_output_idx.to_string();
        if bond.stereo() == BondStereo::Any {
            ctu_labels.push(label);
            continue;
        }

        let Some([stereo_begin, stereo_end]) = bond.stereo_atoms() else {
            continue;
        };

        let begin_atom = bond.begin();
        let end_atom = bond.end();
        let mut need_swap = false;

        if incident_bonds(molecule, begin_atom).len() > 2 {
            let Some(stereo_begin_order) = atom_order_positions
                .get(stereo_begin.index())
                .and_then(|position| *position)
            else {
                continue;
            };
            for neighbor_bond in incident_bonds(molecule, begin_atom) {
                if neighbor_bond == bond_id {
                    continue;
                }
                let Some(neighbor_atom) =
                    bond_other_atom(&molecule.bonds()[neighbor_bond.index()], begin_atom)
                else {
                    continue;
                };
                if neighbor_atom == end_atom || neighbor_atom == stereo_begin {
                    continue;
                }
                if atom_order_positions
                    .get(neighbor_atom.index())
                    .and_then(|position| *position)
                    .is_some_and(|neighbor_order| neighbor_order < stereo_begin_order)
                {
                    need_swap = !need_swap;
                }
            }
        }

        if incident_bonds(molecule, end_atom).len() > 2 {
            let Some(stereo_end_order) = atom_order_positions
                .get(stereo_end.index())
                .and_then(|position| *position)
            else {
                continue;
            };
            for neighbor_bond in incident_bonds(molecule, end_atom) {
                if neighbor_bond == bond_id {
                    continue;
                }
                let Some(neighbor_atom) =
                    bond_other_atom(&molecule.bonds()[neighbor_bond.index()], end_atom)
                else {
                    continue;
                };
                if neighbor_atom == begin_atom || neighbor_atom == stereo_end {
                    continue;
                }
                if atom_order_positions
                    .get(neighbor_atom.index())
                    .and_then(|position| *position)
                    .is_some_and(|neighbor_order| neighbor_order < stereo_end_order)
                {
                    need_swap = !need_swap;
                }
            }
        }

        if bond.stereo() == BondStereo::Cis || need_swap {
            c_labels.push(label);
        } else {
            t_labels.push(label);
        }
    }

    let c = if c_labels.is_empty() {
        String::new()
    } else {
        format!("c:{}", c_labels.join(","))
    };
    let t = if t_labels.is_empty() {
        String::new()
    } else {
        format!("t:{}", t_labels.join(","))
    };
    let ctu = if ctu_labels.is_empty() {
        String::new()
    } else {
        format!("ctu:{}", ctu_labels.join(","))
    };
    format!("{c}{t}{ctu}")
}

// BEGIN RDKIT CPP FUNCTION get_linknodes_block
// RDKit✔️✔️: std::string get_linknodes_block(const ROMol &mol,
// RDKit✔️✔️:                                 const std::vector<unsigned int> &atomOrder) {
// RDKit✔️✔️:   bool strict = false;
// RDKit✔️✔️:   auto linkNodes = MolEnumerator::utils::getMolLinkNodes(mol, strict);
// RDKit✔️✔️:   if (linkNodes.empty()) { return ""; }
// RDKit✔️✔️:   std::stringstream res;
// RDKit✔️✔️:   res << "LN:";
// RDKit✔️✔️:   for (const auto &ln : linkNodes) {
// RDKit✔️✔️:     unsigned int atomIdx = atomOrder[ln.bondAtoms[0].first];
// RDKit✔️✔️:     res << atomIdx << ":" << ln.minRep << "." << ln.maxRep;
// RDKit✔️✔️:     if (mol.getAtomWithIdx(ln.bondAtoms[0].first)->getDegree() > 2) {
// RDKit✔️✔️:       res << "." << atomOrder[ln.bondAtoms[0].second] << "."
// RDKit✔️✔️:           << atomOrder[ln.bondAtoms[1].second];
// RDKit✔️✔️:     }
// RDKit✔️✔️:     res << ",";
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (!resStr.empty() && resStr.back() == ',') { resStr.pop_back(); }
// RDKit✔️✔️:   return resStr;
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION get_linknodes_block
fn write_cx_linknodes_block(molecule: &Molecule, atom_order: &[AtomId]) -> String {
    let Some(raw_link_nodes) = molecule.prop("_MolFileLinkNodes") else {
        return String::new();
    };

    let mut atom_order_positions: Vec<Option<usize>> = vec![None; molecule.atoms().len()];
    for (position, atom_id) in atom_order.iter().copied().enumerate() {
        if atom_id.index() < atom_order_positions.len() {
            atom_order_positions[atom_id.index()] = Some(position);
        }
    }

    let mut entries: Vec<String> = Vec::new();
    for record in raw_link_nodes
        .split('|')
        .filter(|part| !part.trim().is_empty())
    {
        let values: Vec<usize> = record
            .split_whitespace()
            .filter_map(|part| part.parse::<usize>().ok())
            .collect();
        if values.len() < 5 {
            continue;
        }
        let min_rep = values[0];
        let max_rep = values[1];
        let pair_count = values[2];
        let required = 3usize.saturating_add(pair_count.saturating_mul(2));
        if values.len() < required || pair_count == 0 {
            continue;
        }

        let Some(center_atom_one_based) = values.get(3).copied() else {
            continue;
        };
        let Some(center_atom_idx) = center_atom_one_based.checked_sub(1) else {
            continue;
        };
        let Some(center_output_idx) = atom_order_positions
            .get(center_atom_idx)
            .and_then(|position| *position)
        else {
            continue;
        };

        let mut entry = format!("{center_output_idx}:{min_rep}.{max_rep}");
        let center_atom = AtomId::new(center_atom_idx);
        if incident_bonds(molecule, center_atom).len() > 2 && pair_count >= 2 {
            let Some(first_neighbor_one_based) = values.get(4).copied() else {
                continue;
            };
            let Some(second_neighbor_one_based) = values.get(6).copied() else {
                continue;
            };
            let Some(first_neighbor_idx) = first_neighbor_one_based.checked_sub(1) else {
                continue;
            };
            let Some(second_neighbor_idx) = second_neighbor_one_based.checked_sub(1) else {
                continue;
            };
            let Some(first_neighbor_out) = atom_order_positions
                .get(first_neighbor_idx)
                .and_then(|position| *position)
            else {
                continue;
            };
            let Some(second_neighbor_out) = atom_order_positions
                .get(second_neighbor_idx)
                .and_then(|position| *position)
            else {
                continue;
            };
            entry.push_str(&format!(".{first_neighbor_out}.{second_neighbor_out}"));
        }
        entries.push(entry);
    }

    if entries.is_empty() {
        String::new()
    } else {
        format!("LN:{}", entries.join(","))
    }
}

// BEGIN RDKIT CPP FUNCTION get_sgroup_polymer_block
// RDKit✔️✔️: std::string get_sgroup_polymer_block(
// RDKit✔️✔️:     const ROMol &mol, const std::vector<unsigned int> &atomOrder,
// RDKit✔️✔️:     const std::vector<unsigned int> &bondOrder) {
// RDKit✔️✔️:   // builds Sg:type:atoms:label:connect:headCrossings:tailCrossings:
// RDKit✔️✔️:   // using reverse typemap and output-order indices
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION get_sgroup_polymer_block
fn write_cx_sgroup_polymer_block(
    molecule: &Molecule,
    atom_order: &[AtomId],
    bond_order: &[BondId],
) -> String {
    if molecule.substance_groups().is_empty() {
        return String::new();
    }

    let rev_atom_order = cx_atom_output_positions(atom_order, molecule.num_atoms());
    let rev_bond_order = cx_bond_output_positions(bond_order, molecule.num_bonds());

    let mut entries: Vec<String> = Vec::new();
    for sgroup in molecule.substance_groups() {
        let Some(type_code) = writer_polymer_sgroup_type_code(sgroup) else {
            continue;
        };
        if sgroup.atoms().is_empty() {
            continue;
        }

        let mut atom_parts: Vec<String> = Vec::new();
        for atom in sgroup.atoms() {
            let Some(out_idx) = rev_atom_order.get(atom.index()).and_then(|v| *v) else {
                continue;
            };
            atom_parts.push(out_idx.to_string());
        }
        if atom_parts.is_empty() {
            continue;
        }

        let label = sgroup.label().unwrap_or_default();
        let connect = sgroup
            .connection()
            .map(|value| match value {
                crate::sgroup::SGroupConnection::HeadToHead => "hh".to_string(),
                crate::sgroup::SGroupConnection::HeadToTail => "ht".to_string(),
                crate::sgroup::SGroupConnection::Either => "eu".to_string(),
                crate::sgroup::SGroupConnection::Unknown(text) => text.to_ascii_lowercase(),
            })
            .or_else(|| {
                sgroup
                    .props()
                    .get("CONNECT")
                    .map(|value| value.to_ascii_lowercase())
            })
            .unwrap_or_default();

        let head_crossings = writer_parse_sgroup_crossings(sgroup, "_headCrossings")
            .or_else(|| writer_parse_sgroup_crossings(sgroup, "XBHEAD"))
            .unwrap_or_default();
        let tail_crossings = writer_parse_sgroup_crossings(sgroup, "_tailCrossings")
            .or_else(|| writer_parse_sgroup_crossings(sgroup, "XBCORR"))
            .unwrap_or_default();

        let head_field = if head_crossings.len() > 1 {
            head_crossings
                .iter()
                .filter_map(|idx| rev_bond_order.get(*idx).and_then(|value| *value))
                .map(|idx| idx.to_string())
                .collect::<Vec<_>>()
                .join(",")
        } else {
            String::new()
        };
        let tail_field = if tail_crossings.len() > 2 {
            tail_crossings
                .iter()
                .enumerate()
                .filter_map(|(i, idx)| {
                    if i % 2 == 1 {
                        rev_bond_order
                            .get(*idx)
                            .and_then(|value| *value)
                            .map(|idx| idx.to_string())
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>()
                .join(",")
        } else {
            String::new()
        };

        entries.push(format!(
            "Sg:{type_code}:{}:{label}:{connect}:{head_field}:{tail_field}:",
            atom_parts.join(",")
        ));
    }

    entries.join(",")
}

fn write_cx_polymer_sgroups(
    molecule: &Molecule,
    atom_order: &[AtomId],
    bond_order: &[BondId],
) -> String {
    write_cx_sgroup_polymer_block(molecule, atom_order, bond_order)
}

fn writer_is_data_sgroup(sgroup: &crate::SubstanceGroup) -> bool {
    matches!(sgroup.kind(), crate::sgroup::SubstanceGroupKind::Data)
        || sgroup
            .props()
            .get("TYPE")
            .is_some_and(|value| value == "DAT")
}

fn writer_data_sgroup_field_name(sgroup: &crate::SubstanceGroup) -> String {
    sgroup
        .data()
        .and_then(|data| data.field_name.clone())
        .or_else(|| sgroup.props().get("FIELDNAME").cloned())
        .unwrap_or_default()
}

fn writer_data_sgroup_values(sgroup: &crate::SubstanceGroup) -> Vec<String> {
    if !sgroup.data_fields().is_empty() {
        return sgroup.data_fields().to_vec();
    }
    if let Some(values) = sgroup
        .data()
        .map(|data| {
            data.values
                .iter()
                .filter(|value| !value.is_empty())
                .cloned()
                .collect::<Vec<_>>()
        })
        .filter(|values| !values.is_empty())
    {
        return values;
    }
    sgroup
        .props()
        .get("DATAFIELDS")
        .map(|value| vec![value.clone()])
        .unwrap_or_default()
}

fn writer_data_sgroup_query_op(sgroup: &crate::SubstanceGroup) -> String {
    sgroup
        .data()
        .and_then(|data| data.query_op.clone())
        .or_else(|| sgroup.props().get("QUERYOP").cloned())
        .unwrap_or_default()
}

fn writer_data_sgroup_field_info(sgroup: &crate::SubstanceGroup) -> String {
    sgroup
        .data()
        .and_then(|data| data.field_info.clone())
        .or_else(|| sgroup.props().get("FIELDINFO").cloned())
        .unwrap_or_default()
}

fn writer_data_sgroup_field_tag(sgroup: &crate::SubstanceGroup) -> String {
    sgroup
        .data()
        .and_then(|data| data.units.clone())
        .or_else(|| sgroup.props().get("FIELDTAG").cloned())
        .unwrap_or_default()
}

fn writer_polymer_sgroup_type_code(sgroup: &crate::SubstanceGroup) -> Option<String> {
    let code = match sgroup.kind() {
        crate::sgroup::SubstanceGroupKind::StructuralRepeatUnit => "n",
        crate::sgroup::SubstanceGroupKind::Monomer => "mon",
        crate::sgroup::SubstanceGroupKind::Mer => "mer",
        crate::sgroup::SubstanceGroupKind::Copolymer => {
            let subtype = sgroup
                .subtype()
                .map(|value| value.to_ascii_uppercase())
                .or_else(|| {
                    sgroup
                        .props()
                        .get("SUBTYPE")
                        .map(|value| value.to_ascii_uppercase())
                });
            match subtype.as_deref() {
                Some("ALT") => "alt",
                Some("RAN") => "ran",
                Some("BLO") | Some("BLK") => "blk",
                _ => "co",
            }
        }
        crate::sgroup::SubstanceGroupKind::Crosslink => "xl",
        crate::sgroup::SubstanceGroupKind::Modification => "mod",
        crate::sgroup::SubstanceGroupKind::MixtureComponent => "mix",
        crate::sgroup::SubstanceGroupKind::Formulation => "f",
        crate::sgroup::SubstanceGroupKind::AnyPolymer => "any",
        crate::sgroup::SubstanceGroupKind::Graft => "grf",
        crate::sgroup::SubstanceGroupKind::Generic(value) if value.eq_ignore_ascii_case("GEN") => {
            "gen"
        }
        crate::sgroup::SubstanceGroupKind::Generic(value) if value.eq_ignore_ascii_case("COM") => {
            "c"
        }
        _ => return None,
    };
    Some(code.to_string())
}

fn writer_parse_sgroup_crossings(sgroup: &crate::SubstanceGroup, key: &str) -> Option<Vec<usize>> {
    let raw = sgroup.props().get(key)?;
    let parsed: Vec<usize> = raw
        .split(',')
        .filter_map(|part| part.trim().parse::<usize>().ok())
        .collect();
    Some(parsed)
}

// BEGIN RDKIT CPP FUNCTION get_sgroup_hierarchy_block
// RDKit✔️✔️: std::string get_sgroup_hierarchy_block(const ROMol &mol) {
// RDKit✔️✔️:   // builds SgH:parentOutputIdx:childOutputIdx.childOutputIdx,...
// RDKit✔️✔️:   // using temporary _cxsmilesOutputIndex assigned by prior SGroup writers
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION get_sgroup_hierarchy_block
fn write_cx_sgroup_hierarchy_block(
    molecule: &Molecule,
    atom_order: &[AtomId],
    bond_order: &[BondId],
    include_sgroups: bool,
    include_polymer: bool,
) -> String {
    let sgroups = molecule.substance_groups();
    if sgroups.is_empty() {
        return String::new();
    }
    let output_index_by_sgroup_id = writer_cx_hierarchy_output_indices(
        molecule,
        atom_order,
        bond_order,
        include_sgroups,
        include_polymer,
    );

    if output_index_by_sgroup_id.is_empty() {
        return String::new();
    }

    let mut accum: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (fallback_idx, sgroup) in sgroups.iter().enumerate() {
        let child_id = writer_sgroup_index_value(sgroup, fallback_idx);
        let Some(child_output_idx) = output_index_by_sgroup_id.get(&child_id).copied() else {
            continue;
        };
        let Some(parent_id) = writer_sgroup_parent_value(sgroup) else {
            continue;
        };
        let Some(parent_output_idx) = output_index_by_sgroup_id.get(&parent_id).copied() else {
            continue;
        };
        accum
            .entry(parent_output_idx)
            .or_default()
            .push(child_output_idx);
    }

    if accum.is_empty() {
        return String::new();
    }

    let mut entries: Vec<String> = Vec::new();
    for (parent, children) in accum {
        if children.is_empty() {
            continue;
        }
        let child_text = children
            .iter()
            .map(|child| child.to_string())
            .collect::<Vec<_>>()
            .join(".");
        entries.push(format!("{parent}:{child_text}"));
    }

    if entries.is_empty() {
        String::new()
    } else {
        format!("SgH:{}", entries.join(","))
    }
}

fn writer_cx_hierarchy_output_indices(
    molecule: &Molecule,
    atom_order: &[AtomId],
    bond_order: &[BondId],
    include_sgroups: bool,
    include_polymer: bool,
) -> BTreeMap<usize, usize> {
    let atom_set = atom_order.iter().copied().collect::<BTreeSet<_>>();
    let bond_set = bond_order.iter().copied().collect::<BTreeSet<_>>();
    let mut output_index_by_sgroup_id = BTreeMap::new();
    let mut next_output_index = 0usize;

    if include_sgroups {
        for (fallback_idx, sgroup) in molecule.substance_groups().iter().enumerate() {
            if !writer_is_data_sgroup(sgroup)
                || !sgroup.atoms().iter().any(|atom| atom_set.contains(atom))
            {
                continue;
            }
            let sgroup_id = writer_sgroup_index_value(sgroup, fallback_idx);
            if let std::collections::btree_map::Entry::Vacant(entry) =
                output_index_by_sgroup_id.entry(sgroup_id)
            {
                entry.insert(next_output_index);
                next_output_index += 1;
            }
        }
        for (fallback_idx, sgroup) in molecule.substance_groups().iter().enumerate() {
            if writer_is_data_sgroup(sgroup)
                || writer_polymer_sgroup_type_code(sgroup).is_some()
                || (!sgroup.atoms().iter().any(|atom| atom_set.contains(atom))
                    && !sgroup.bonds().iter().any(|bond| bond_set.contains(bond)))
            {
                continue;
            }
            let sgroup_id = writer_sgroup_index_value(sgroup, fallback_idx);
            if let std::collections::btree_map::Entry::Vacant(entry) =
                output_index_by_sgroup_id.entry(sgroup_id)
            {
                entry.insert(next_output_index);
                next_output_index += 1;
            }
        }
    }

    if include_polymer {
        for (fallback_idx, sgroup) in molecule.substance_groups().iter().enumerate() {
            if writer_polymer_sgroup_type_code(sgroup).is_none()
                || !sgroup.atoms().iter().any(|atom| atom_set.contains(atom))
            {
                continue;
            }
            let sgroup_id = writer_sgroup_index_value(sgroup, fallback_idx);
            if let std::collections::btree_map::Entry::Vacant(entry) =
                output_index_by_sgroup_id.entry(sgroup_id)
            {
                entry.insert(next_output_index);
                next_output_index += 1;
            }
        }
    }

    output_index_by_sgroup_id
}

fn writer_sgroup_index_value(sgroup: &crate::SubstanceGroup, _fallback_idx: usize) -> usize {
    sgroup
        .props()
        .get("index")
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or_else(|| sgroup.id().index())
}

fn writer_sgroup_parent_value(sgroup: &crate::SubstanceGroup) -> Option<usize> {
    if let Some(parent) = sgroup.parent() {
        return Some(parent.index());
    }
    sgroup
        .props()
        .get("PARENT")
        .and_then(|value| value.parse::<usize>().ok())
}

// RDKit✔️✔️: std::string getAtomChiralityInfo(const Atom *atom) {
// RDKit✔️✔️:   auto allowNontet = Chirality::getAllowNontetrahedralChirality();
// RDKit✔️✔️:   std::string atString;
// RDKit✔️✔️:   switch (atom->getChiralTag()) {
// RDKit✔️✔️:     case Atom::CHI_TETRAHEDRAL_CW: atString = "@@"; break;
// RDKit✔️✔️:     case Atom::CHI_TETRAHEDRAL_CCW: atString = "@"; break;
// RDKit✔️✔️:     default: break;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (atString.empty() && allowNontet) {
// RDKit✔️✔️:     switch (atom->getChiralTag()) {
// RDKit✔️✔️:       case Atom::CHI_SQUAREPLANAR: atString = "@SP"; break;
// RDKit✔️✔️:       case Atom::CHI_TRIGONALBIPYRAMIDAL: atString = "@TB"; break;
// RDKit✔️✔️:       case Atom::CHI_OCTAHEDRAL: atString = "@OH"; break;
// RDKit✔️✔️:       default: break;
// RDKit✔️✔️:     }
// RDKit✔️✔️:     if (!atString.empty()) {
// RDKit✔️✔️:       int permutation = 0;
// RDKit✔️✔️:       if (atom->getChiralTag() > Atom::ChiralType::CHI_OTHER &&
// RDKit✔️✔️:           atom->getPropIfPresent(common_properties::_chiralPermutation,
// RDKit✔️✔️:                                  permutation) &&
// RDKit✔️✔️:           !SmilesParseOps::checkChiralPermutation(atom->getChiralTag(),
// RDKit✔️✔️:                                                   permutation)) {
// RDKit✔️✔️:         throw ValueErrorException("bad chirality spec");
// RDKit✔️✔️:       } else if (permutation) {
// RDKit✔️✔️:         atString += std::to_string(permutation);
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return atString;
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION getAtomChiralityInfo
fn get_atom_chirality_info_with_inversion(
    molecule: &Molecule,
    atom: AtomId,
    chiral_tag_override: Option<ChiralTag>,
    invert: bool,
    permutation_override: Option<u32>,
) -> Result<String, SmilesWriteError> {
    let atom = &molecule.atoms()[atom.index()];
    let base_chiral_tag = chiral_tag_override.unwrap_or_else(|| atom.chiral_tag());
    let chiral_tag = if invert {
        match base_chiral_tag {
            ChiralTag::TetrahedralCw => ChiralTag::TetrahedralCcw,
            ChiralTag::TetrahedralCcw => ChiralTag::TetrahedralCw,
            other => other,
        }
    } else {
        base_chiral_tag
    };
    match chiral_tag {
        ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw => {
            let perm = permutation_override
                .or_else(|| atom.chiral_permutation())
                .unwrap_or(0);
            let res = if perm % 2 == 0 {
                match chiral_tag {
                    ChiralTag::TetrahedralCw => "@@",
                    _ => "@",
                }
            } else {
                match chiral_tag {
                    ChiralTag::TetrahedralCw => "@",
                    _ => "@@",
                }
            };
            Ok(res.to_string())
        }
        ChiralTag::SquarePlanar | ChiralTag::TrigonalBipyramidal | ChiralTag::Octahedral => {
            let mut res = match chiral_tag {
                ChiralTag::SquarePlanar => "@SP".to_string(),
                ChiralTag::TrigonalBipyramidal => "@TB".to_string(),
                ChiralTag::Octahedral => "@OH".to_string(),
                _ => unreachable!(),
            };
            if let Some(permutation) = permutation_override.or_else(|| atom.chiral_permutation()) {
                validate_writer_chiral_permutation(chiral_tag, permutation)?;
                if permutation != 0 {
                    res.push_str(&permutation.to_string());
                }
            }
            Ok(res)
        }
        _ => Ok(String::new()),
    }
}

fn validate_writer_chiral_permutation(
    chiral_tag: ChiralTag,
    permutation: u32,
) -> Result<(), SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION SmilesParseOps::checkChiralPermutation
    // RDKit✔️✔️: bool checkChiralPermutation(int chiralTag, int permutation) {
    // RDKit✔️✔️:   if (chiralTag > RDKit::Atom::ChiralType::CHI_OTHER &&
    // RDKit✔️✔️:       permutationLimits.find(chiralTag) != permutationLimits.end() &&
    // RDKit✔️✔️:       (permutation < 0 || permutation > permutationLimits.at(chiralTag))) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION SmilesParseOps::checkChiralPermutation
    let limit = match chiral_tag {
        ChiralTag::Allene => Some(2),
        ChiralTag::SquarePlanar => Some(3),
        ChiralTag::TrigonalBipyramidal => Some(20),
        ChiralTag::Octahedral => Some(30),
        _ => None,
    };
    if let Some(limit) = limit
        && permutation > limit
    {
        return Err(SmilesWriteError::InvalidChiralPermutation {
            chiral_tag,
            permutation,
            limit,
        });
    }
    Ok(())
}

fn atom_needs_bracket(
    molecule: &Molecule,
    atom: AtomId,
    at_string: &str,
    params: &SmilesWriteParams,
) -> Result<bool, SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION atomNeedsBracket
    // RDKit✔️✔️: bool atomNeedsBracket(const Atom *atom, const std::string &atString,
    // RDKit✔️✔️:                       const SmilesWriteParams &params) {
    // RDKit✔️✔️:   PRECONDITION(atom, "null atom");
    // RDKit✔️✔️:   auto num = atom->getAtomicNum();
    // RDKit✔️✔️:   if (!inOrganicSubset(num)) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (atom->getFormalCharge()) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (params.doIsomericSmiles && (atom->getIsotope() || !atString.empty())) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (atom->hasProp(common_properties::molAtomMapNumber)) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   const INT_VECT &defaultVs = PeriodicTable::getTable()->getValenceList(num);
    // RDKit✔️✔️:   int totalValence = atom->getTotalValence();
    // RDKit✔️✔️:   bool nonStandard = false;
    // RDKit✔️✔️:   if (atom->getNumRadicalElectrons()) {
    // RDKit✔️✔️:     nonStandard = true;
    // RDKit✔️✔️:   } else if ((num == 7 || num == 15) && atom->getIsAromatic() &&
    // RDKit✔️✔️:              atom->getNumExplicitHs()) {
    // RDKit✔️✔️:     nonStandard = true;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     nonStandard = (totalValence != defaultVs.front() && atom->getTotalNumHs());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (nonStandard) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // check for bonds to a metal
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION atomNeedsBracket
    let atom_id = atom;
    let atom = &molecule.atoms()[atom.index()];
    if !in_organic_subset(atom.atomic_number())? {
        return Ok(true);
    }
    if atom.formal_charge() != 0 || atom.atom_map().is_some() {
        return Ok(true);
    }
    if params.do_isomeric_smiles && (atom.isotope().is_some() || !at_string.is_empty()) {
        return Ok(true);
    }
    // RDKit✔️✔️: if (atom->getNumRadicalElectrons()) { nonStandard = true; }
    if atom.radical_electrons() != 0 {
        return Ok(true);
    }
    // RDKit✔️✔️: else if ((num == 7 || num == 15) && atom->getIsAromatic() &&
    // RDKit✔️✔️:            atom->getNumExplicitHs()) {
    // RDKit✔️✔️:   nonStandard = true;
    // RDKit✔️✔️: }
    if matches!(atom.atomic_number(), 7 | 15) && atom.is_aromatic() && atom.explicit_hydrogens() > 0
    {
        return Ok(true);
    }
    // RDKit✔️✔️:   const INT_VECT &defaultVs = PeriodicTable::getTable()->getValenceList(num);
    // RDKit✔️✔️:   int totalValence = atom->getTotalValence();
    // RDKit✔️✔️:   bool nonStandard = false;
    // RDKit✔️✔️:   if (atom->getNumRadicalElectrons()) {
    // RDKit✔️✔️:     nonStandard = true;
    // RDKit✔️✔️:   } else if ((num == 7 || num == 15) && atom->getIsAromatic() &&
    // RDKit✔️✔️:              atom->getNumExplicitHs()) {
    // RDKit✔️✔️:     nonStandard = true;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     nonStandard = (totalValence != defaultVs.front() && atom->getTotalNumHs());
    // RDKit✔️✔️:   }
    if let Ok(Some(valence_list)) = crate::valence::rdkit_valence_list(atom.atomic_number()) {
        if let Some(&default_valence) = valence_list.first() {
            let total_valence = total_valence_for_writer(molecule, atom_id).unwrap_or_else(|| {
                molecule
                    .bonds()
                    .iter()
                    .fold(i32::from(atom.explicit_hydrogens()), |acc, bond| {
                        if bond.begin() == atom_id || bond.end() == atom_id {
                            acc + crate::valence::bond_valence_contrib(bond, atom_id)
                                .unwrap_or(0.0)
                                .round() as i32
                        } else {
                            acc
                        }
                    })
            });
            let total_hs = total_num_hydrogens_for_writer(molecule, atom_id) as i32;
            if total_valence != default_valence && total_hs > 0 {
                return Ok(true);
            }
        }
    }
    if molecule.bonds().iter().any(|bond| {
        bond_other_atom(bond, atom_id)
            .map(|other| rdkit_query_ops_is_metal(molecule.atoms()[other.index()].atomic_number()))
            .unwrap_or(false)
    }) {
        return Ok(true);
    }
    Ok(false)
}

fn rdkit_query_ops_is_metal(atomic_number: u8) -> bool {
    // BEGIN RDKIT CPP FUNCTION QueryOps::makeMAtomQuery / QueryOps::isMetal
    // RDKit✔️✔️: // !#0!#1!#2!#5!#6!#7!#8!#9!#10!#14!#15!#16!#17!#18!#33!#34!#35!#36!#52!#53!#54!#85!#86
    // RDKit✔️✔️: bool isMetal(const Atom &atom) { ... makeMAtomQuery()->Match(&atom); }
    // END RDKIT CPP FUNCTION QueryOps::makeMAtomQuery / QueryOps::isMetal
    !matches!(
        atomic_number,
        0 | 1
            | 2
            | 5
            | 6
            | 7
            | 8
            | 9
            | 10
            | 14
            | 15
            | 16
            | 17
            | 18
            | 33
            | 34
            | 35
            | 36
            | 52
            | 53
            | 54
            | 85
            | 86
    )
}

fn update_property_cache_for_smiles(molecule: &mut Molecule) -> Result<(), SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles writer property-cache prep
    // RDKit✔️✔️: for (auto atom : tmol->atoms()) {
    // RDKit✔️✔️:   atom->updatePropertyCache(false);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles writer property-cache prep
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles ring-info lifetime
    // RDKit✔️✔️: // writer prep here updates atom property cache only; it does not
    // RDKit✔️✔️: // proactively recompute ring info. Traversal/ranking code consults the
    // RDKit✔️✔️: // molecule's existing RingInfo and only computes fast rings on demand in
    // RDKit✔️✔️: // specific canonical-ranking paths.
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles ring-info lifetime
    *molecule = molecule.with_assigned_valence()?;
    Ok(())
}

fn clear_fragment_temp_molecule_computed_stereo_props_for_writer(molecule: &mut Molecule) {
    // BEGIN RDKIT CPP FUNCTION MolOps::getMolFrags fragment temporary molecule path
    // RDKit❗✔️: if (comp.size() == 1 || ...) {
    // RDKit❗✔️:   SubsetOptions opts{.sanitize = sanitizeFrags,
    // RDKit❗✔️:                      .clearComputedProps = true,
    // RDKit❗✔️:                      .copyCoordinates = copyConformers,
    // RDKit❗✔️:                      .method = SubsetMethod::BONDS_BETWEEN_ATOMS};
    // RDKit❗✔️:   auto submol = copyMolSubset(mol, atoms, info, opts);
    // RDKit❗✔️: } else {
    // RDKit❗✔️:   frag->beginBatchEdit();
    // RDKit❗✔️:   ...
    // RDKit❗✔️:   frag->commitBatchEdit();
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION MolOps::getMolFrags fragment temporary molecule path
    // BEGIN RDKIT CPP FUNCTION RWMol::removeAtom / commitBatchEdit computed prop reset
    // RDKit❗✔️: // clear computed properties and reset our ring info structure
    // RDKit❗✔️: // they are pretty likely to be wrong now:
    // RDKit❗✔️: if (clearProps) {
    // RDKit❗✔️:   clearComputedProps(true);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION RWMol::removeAtom / commitBatchEdit computed prop reset
    // BEGIN RDKIT CPP FUNCTION ROMol::clearComputedProps
    // RDKit❗✔️: void ROMol::clearComputedProps(bool includeRings) const {
    // RDKit❗✔️:   if (includeRings) {
    // RDKit❗✔️:     this->dp_ringInfo->reset();
    // RDKit❗✔️:   }
    // RDKit❗✔️:   RDProps::clearComputedProps();
    // RDKit❗✔️:   for (auto atom : atoms()) {
    // RDKit❗✔️:     atom->clearComputedProps();
    // RDKit❗✔️:   }
    // RDKit❗✔️:   for (auto bond : bonds()) {
    // RDKit❗✔️:     bond->clearComputedProps();
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION ROMol::clearComputedProps
    let fragment_count = crate::notation::fragment::get_fragment_atom_mapping(molecule)
        .iter()
        .copied()
        .max()
        .map_or(0, |max_fragment| max_fragment + 1);
    if fragment_count <= 1 {
        return;
    }

    molecule.properties_mut().clear_prop("_StereochemDone");
    for atom in &mut molecule.topology_block_mut().atoms {
        atom.clear_prop("_CIPCode");
        atom.clear_prop("_CIPRank");
        atom.clear_prop("_ChiralityPossible");
        atom.clear_prop("_ringStereochemCand");
        atom.clear_prop("_ringStereoAtoms");
    }
    for bond in &mut molecule.topology_block_mut().bonds {
        bond.clear_prop("_CIPCode");
    }
}

fn assign_stereochemistry_for_smiles(
    molecule: &mut Molecule,
    clean_stereo: bool,
) -> Result<(), SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION Chirality::assignStereochemistry dispatcher tail
    // RDKit❗✔️: void assignStereochemistry(ROMol &mol, bool cleanIt, bool force,
    // RDKit❗✔️:                            bool flagPossibleStereoCenters) {
    // RDKit❗✔️:   if (!force && mol.hasProp(common_properties::_StereochemDone)) {
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   ...
    // RDKit❗✔️:   mol.setProp(common_properties::_StereochemDone, 1, true);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION Chirality::assignStereochemistry dispatcher tail
    // RDKit calls assignStereochemistry() on the temporary writer molecule
    // unless stereochemistry was already prepared. COSMolKit's current stereo
    // assignment reads typed state and validates/detects writer-visible stereo;
    // future ports can make this mutating without changing the writer contract.
    // BEGIN RDKIT CPP FUNCTION Chirality::stereoPerception cleanIt prelude
    // RDKit❗✔️: void stereoPerception(ROMol &mol, bool cleanIt,
    // RDKit❗✔️:                       bool flagPossibleStereoCenters) {
    // RDKit❗✔️:   if (cleanIt) {
    // RDKit❗✔️:     for (auto atom : mol.atoms()) {
    // RDKit❗✔️:       atom->clearProp(common_properties::_CIPCode);
    // RDKit❗✔️:       atom->clearProp(common_properties::_ChiralityPossible);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     for (auto bond : mol.bonds()) {
    // RDKit❗✔️:       bond->clearProp(common_properties::_CIPCode);
    // RDKit❗✔️:       if (bond->getBondDir() == Bond::BondDir::EITHERDOUBLE) {
    // RDKit❗✔️:         bond->setStereo(Bond::BondStereo::STEREOANY);
    // RDKit❗✔️:         bond->getStereoAtoms().clear();
    // RDKit❗✔️:         bond->setBondDir(Bond::BondDir::NONE);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // we need cis/trans markers on the double bonds... set those now:
    // RDKit❗✔️:   MolOps::setBondStereoFromDirections(mol);
    // RDKit❗✔️:   ...
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION Chirality::stereoPerception cleanIt prelude
    //
    // This ports the writer-visible `cleanIt` property reset and EITHERDOUBLE
    // normalization exactly. The downstream `findPotentialStereo()` /
    // `updateDoubleBondStereo()` pipeline is still only partially reproduced by
    // the current typed-state helpers below, so the behavior marker remains `❗`.
    if clean_stereo {
        for atom in &mut molecule.topology_block_mut().atoms {
            atom.clear_prop("_CIPCode");
            atom.clear_prop("_ChiralityPossible");
        }
        for bond in &mut molecule.topology_block_mut().bonds {
            bond.clear_prop("_CIPCode");
            if bond.direction() == BondDirection::EitherDouble {
                bond.set_stereo(BondStereo::Any);
                bond.set_stereo_atoms(None);
                bond.set_direction(BondDirection::None);
            }
        }
    }
    ensure_fast_rings_for_writer_stereo_perception(molecule)?;
    crate::stereo::assign_stereochemistry(molecule)?;
    assign_double_bond_stereo_for_writer_working_copy(molecule)?;
    let ranks = crate::stereo::assign_atom_cip_ranks(molecule)?;
    for (atom_idx, cip_code) in crate::stereo::assign_atom_chiral_codes(molecule, &ranks)? {
        if let Some(atom_mut) = molecule.topology_block_mut().atoms.get_mut(atom_idx) {
            atom_mut.set_prop("_CIPCode", cip_code);
        }
    }
    for (atom_idx, rank) in ranks.iter().copied().enumerate() {
        if let Some(atom_mut) = molecule.topology_block_mut().atoms.get_mut(atom_idx) {
            atom_mut.set_prop("_CIPRank", rank.to_string());
        }
    }
    if clean_stereo {
        apply_clean_stereo_ring_special_cases_for_writer(molecule, &ranks)?;
    }
    molecule.properties_mut().set_prop("_StereochemDone", "1");
    Ok(())
}

fn ensure_fast_rings_for_writer_stereo_perception(
    molecule: &mut Molecule,
) -> Result<(), SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION Chirality::legacyStereoPerception ring prelude
    // RDKit✔️✔️: void legacyStereoPerception(ROMol &mol, bool cleanIt,
    // RDKit✔️✔️:                             bool flagPossibleStereoCenters) {
    // RDKit✔️✔️:   mol.clearProp("_needsDetectBondStereo");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // later we're going to need ring information, get it now if we don't
    // RDKit✔️✔️:   // have it already:
    // RDKit✔️✔️:   // NOTE, if called from the SMART code, the ring info will be DUMMY, and
    // RDKit✔️✔️:   // contains no information
    // RDKit✔️✔️:   if (!mol.getRingInfo()->isFindFastOrBetter()) {
    // RDKit✔️✔️:     MolOps::fastFindRings(mol);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Chirality::legacyStereoPerception ring prelude
    if molecule
        .derived_cache()
        .rings
        .as_ref()
        .is_some_and(crate::RingInfo::is_find_fast_or_better)
    {
        return Ok(());
    }
    let rings = crate::fast_find_rings(molecule)?;
    molecule.derived_cache_mut().rings = Some(rings);
    Ok(())
}

fn apply_clean_stereo_ring_special_cases_for_writer(
    molecule: &mut Molecule,
    ranks: &[u32],
) -> Result<(), SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION legacyStereoPerception cleanIt ring special cases
    // RDKit❗✔️: if (cleanIt) {
    // RDKit❗✔️:   for (auto atom : mol.atoms()) {
    // RDKit❗✔️:     if (atom->hasProp(common_properties::_ringStereochemCand)) {
    // RDKit❗✔️:       atom->clearProp(common_properties::_ringStereochemCand);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (atom->hasProp(common_properties::_ringStereoAtoms)) {
    // RDKit❗✔️:       atom->clearProp(common_properties::_ringStereoAtoms);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   boost::dynamic_bitset<> possibleSpecialCases(mol.getNumAtoms());
    // RDKit❗✔️:   Chirality::findChiralAtomSpecialCases(mol, possibleSpecialCases, atomRanks);
    // RDKit❗✔️:   for (auto atom : mol.atoms()) {
    // RDKit❗✔️:     if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit❗✔️:         !Chirality::hasNonTetrahedralStereo(atom) &&
    // RDKit❗✔️:         !atom->hasProp(common_properties::_CIPCode) &&
    // RDKit❗✔️:         (!possibleSpecialCases[atom->getIdx()] ||
    // RDKit❗✔️:          !atom->hasProp(common_properties::_ringStereoAtoms))) {
    // RDKit❗✔️:       atom->setChiralTag(Atom::CHI_UNSPECIFIED);
    // RDKit❗✔️:       ...
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION legacyStereoPerception cleanIt ring special cases
    for atom in &mut molecule.topology_block_mut().atoms {
        atom.clear_prop("_ringStereochemCand");
        atom.clear_prop("_ringStereoAtoms");
    }

    let special_cases = crate::stereo::find_chiral_atom_special_cases(molecule, ranks)?;
    let special_case_atoms = special_cases
        .iter()
        .map(|case| case.atom_idx)
        .collect::<BTreeSet<_>>();
    for case in &special_cases {
        if let Some(atom_mut) = molecule.topology_block_mut().atoms.get_mut(case.atom_idx) {
            atom_mut.set_prop("_ringStereochemCand", "1");
            atom_mut.set_prop(
                "_ringStereoAtoms",
                serialize_ring_stereo_atoms(&case.ring_stereo_atoms),
            );
        }
    }

    let atom_ids = molecule
        .atoms()
        .iter()
        .map(|atom| atom.id())
        .collect::<Vec<_>>();
    for atom_id in atom_ids {
        let atom = &molecule.atoms()[atom_id.index()];
        if atom.chiral_tag() == ChiralTag::Unspecified
            || crate::stereo::has_non_tetrahedral_stereo(atom)
            || atom.prop("_CIPCode").is_some()
            || (special_case_atoms.contains(&atom_id.index())
                && atom.prop("_ringStereoAtoms").is_some())
        {
            continue;
        }

        if let Some(atom_mut) = molecule.topology_block_mut().atoms.get_mut(atom_id.index()) {
            atom_mut.set_chiral_tag(ChiralTag::Unspecified);
            if atom_mut.explicit_hydrogens() == 1
                && atom_mut.formal_charge() == 0
                && !atom_mut.is_aromatic()
            {
                atom_mut.set_explicit_hydrogens(0);
                atom_mut.set_no_implicit(false);
            }
        }
    }
    *molecule = molecule.with_assigned_valence()?;
    Ok(())
}

pub(crate) fn serialize_ring_stereo_atoms(ring_stereo_atoms: &[(bool, usize)]) -> String {
    ring_stereo_atoms
        .iter()
        .map(|(same_orientation, atom_idx)| {
            let signed = if *same_orientation {
                *atom_idx as i32 + 1
            } else {
                -(*atom_idx as i32 + 1)
            };
            signed.to_string()
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn parse_ring_stereo_atoms_prop(encoded: &str) -> Option<Vec<(bool, usize)>> {
    let mut result = Vec::new();
    for token in encoded.split(',').filter(|token| !token.is_empty()) {
        let entry = token.parse::<i32>().ok()?;
        if entry == 0 {
            return None;
        }
        let atom_idx = entry.unsigned_abs() as usize - 1;
        result.push((entry > 0, atom_idx));
    }
    Some(result)
}

fn writer_ring_stereo_atoms(
    atom: &crate::Atom,
    atom_id: AtomId,
) -> Result<Option<Vec<(bool, usize)>>, SmilesWriteError> {
    let Some(encoded) = atom.prop("_ringStereoAtoms") else {
        if atom.prop("_ringStereochemCand").is_some() {
            return Err(SmilesWriteError::InvalidRingStereoState {
                atom: atom_id.index(),
                requirement: "`_ringStereochemCand` requires `_ringStereoAtoms`",
            });
        }
        return Ok(None);
    };
    parse_ring_stereo_atoms_prop(encoded)
        .ok_or(SmilesWriteError::InvalidRingStereoState {
            atom: atom_id.index(),
            requirement: "`_ringStereoAtoms` must be a valid encoded ring-neighbor list",
        })
        .map(Some)
}

fn assign_double_bond_stereo_for_writer_working_copy(
    molecule: &mut Molecule,
) -> Result<(), SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION Chirality::assignStereochemistry double-bond section
    // RDKit✔️✔️: assignBondStereoCodes(mol, atomRanks);
    // RDKit✔️✔️: if (bond->getStereo() == Bond::STEREONONE) {
    // RDKit✔️✔️:   ... find directed neighbor bonds around double bond ...
    // RDKit✔️✔️:   bond->setStereoAtoms(beginControl, endControl);
    // RDKit✔️✔️:   bond->setStereo(Bond::STEREOE or Bond::STEREOZ);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Chirality::assignStereochemistry double-bond section
    let ranks = crate::stereo::assign_atom_cip_ranks(molecule)?;
    let (assignments, changed) = crate::stereo::assign_bond_stereo_codes(molecule, &ranks);
    if !changed {
        return Ok(());
    }
    for (bond_idx, stereo, begin_control, end_control) in assignments {
        let stereo = match stereo {
            // RDKit✔️✔️: bond->setStereo(Bond::STEREOE);
            crate::stereo::DoubleBondStereo::E => BondStereo::E,
            // RDKit✔️✔️: bond->setStereo(Bond::STEREOZ);
            crate::stereo::DoubleBondStereo::Z => BondStereo::Z,
            crate::stereo::DoubleBondStereo::Unknown => BondStereo::Any,
        };
        let bond = &mut molecule.topology_block_mut().bonds[bond_idx];
        if bond.stereo() == BondStereo::None {
            bond.set_stereo_atoms(Some([AtomId::new(begin_control), AtomId::new(end_control)]));
            bond.set_stereo(stereo);
        }
    }
    molecule
        .derived_cache_mut()
        .invalidate(crate::DerivedState::STEREO | crate::DerivedState::DRAWING);
    Ok(())
}

fn canonicalize_enhanced_stereo_for_smiles(_molecule: &Molecule) -> Result<(), SmilesWriteError> {
    // Enhanced stereo group canonicalization is only needed for CX SMILES.
    // For plain SMILES, stereo groups are already in typed state.
    Ok(())
}

fn cleanup_stereo_groups_for_cx_smiles(molecule: &mut Molecule) -> Result<(), SmilesWriteError> {
    // RDKit cleanupStereoGroups() also applies atropisomer-specific cleanup
    // before CX serialization, moving atropisomer participation from atom
    // members onto bond members in each stereo group.
    let current_groups = molecule.stereo_groups().to_vec();
    if current_groups.is_empty() {
        return Ok(());
    }

    let mut cleaned_groups = Vec::with_capacity(current_groups.len());
    for group in current_groups {
        let mut kept_atoms = Vec::with_capacity(group.atoms().len());
        let mut kept_bonds = group.bonds().to_vec();
        for atom in group.atoms().iter().copied() {
            let mut found_atrop = false;
            for bond_id in incident_bonds(molecule, atom) {
                let bond = &molecule.bonds()[bond_id.index()];
                if matches!(bond.stereo(), BondStereo::AtropCw | BondStereo::AtropCcw) {
                    found_atrop = true;
                    if !kept_bonds.contains(&bond_id) {
                        kept_bonds.push(bond_id);
                    }
                }
            }
            if !found_atrop {
                kept_atoms.push(atom);
            }
        }
        let mut cleaned_group = crate::StereoGroup::new(group.kind(), kept_atoms, kept_bonds);
        if let Some(id) = group.id() {
            cleaned_group = cleaned_group.with_id(id);
        }
        cleaned_groups.push(cleaned_group);
    }

    molecule.topology_block_mut().stereo_groups = cleaned_groups;
    Ok(())
}

// RDKit✔️✔️: void Kekulize(RWMol &mol, bool markAtomsBonds, bool canonical,
// RDKit✔️✔️:               unsigned int maxBackTracks) {
// RDKit✔️✔️:   boost::dynamic_bitset<> atomsToUse(mol.getNumAtoms());
// RDKit✔️✔️:   atomsToUse.set();
// RDKit✔️✔️:   boost::dynamic_bitset<> bondsToUse(mol.getNumBonds());
// RDKit✔️✔️:   bondsToUse.set();
// RDKit✔️✔️:   details::KekulizeFragment(mol, atomsToUse, bondsToUse, markAtomsBonds,
// RDKit✔️✔️:                             canonical, maxBackTracks);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION MolOps::Kekulize
/// Run operation-routed kekulization on the molecule and return a new
/// `Molecule` with resolved bond orders and cleared aromatic flags.
///
/// The writer's `do_kekule=true` path is source-aligned to RDKit's
/// `Kekulize(mol, true, true, 100)` behavior: it uses canonical fragment
/// ranking when resolving aromatic systems before SMILES traversal.
fn kekulize_for_smiles(molecule: &Molecule) -> Result<Molecule, SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION details::KekulizeFragment
    // RDKit✔️✔️:     VECT_INT_VECT allringsSSSR;
    // RDKit✔️✔️:     if (!mol.getRingInfo()->isInitialized()) {
    // RDKit✔️✔️:       MolOps::findSSSR(mol, allringsSSSR);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     const VECT_INT_VECT &allrings =
    // RDKit✔️✔️:         allringsSSSR.empty() ? mol.getRingInfo()->atomRings() : allringsSSSR;
    // END RDKIT CPP FUNCTION details::KekulizeFragment
    let assignment = crate::kekulize::kekulize_assignment(
        molecule,
        molecule.derived_cache().rings.as_ref(),
        true,
        true,
        100,
    )?;
    let mut kekulized = molecule.clone();
    crate::kekulize::apply_kekulize_assignment(kekulized.topology_block_mut(), &assignment);
    // RDKit✔️✔️:         if ((atom->getAtomicNum() == 7 || atom->getAtomicNum() == 15) &&
    // RDKit✔️✔️:             atom->getFormalCharge() == 0 && atom->getNumExplicitHs() == 1) {
    // RDKit✔️✔️:           atom->setNoImplicit(false);
    // RDKit✔️✔️:           atom->setNumExplicitHs(0);
    // RDKit✔️✔️:           atom->updatePropertyCache(false);
    // RDKit✔️✔️:         }
    //
    // COSMolKit stores RDKit's property-cache-derived implicit hydrogen state
    // in `derived_cache.valence` instead of per-atom mutable cache fields.
    // After applying the kekulize assignment, refresh that cache so
    // `GetAtomSmiles()`-equivalent total-H reads see the post-kekulization
    // atom state.
    kekulized.derived_cache_mut().valence = Some(crate::valence::assign_valence_with_options(
        &kekulized,
        crate::ValenceModel::RdkitLike,
        false,
    )?);
    Ok(kekulized)
}

fn normalize_dative_bonds_for_plain_smiles(
    molecule: &mut Molecule,
) -> Result<(), SmilesWriteError> {
    // RDKit✔️✔️:     if (doingCXSmiles || !params.includeDativeBonds) {
    // RDKit✔️✔️:       for (auto bond : tmol->bonds()) {
    // RDKit✔️✔️:         if (bond->getBondType() == Bond::DATIVE) {
    // RDKit✔️✔️:           // we are intentionally only handling DATIVE here. The other weird
    // RDKit✔️✔️:           // RDKit dative alternatives really shouldn't ever show up.
    // RDKit✔️✔️:           bond->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:           // update the explicit valence of the begin atom since the implicit
    // RDKit✔️✔️:           // valence will no longer be properly perceived
    // RDKit✔️✔️:           bond->getBeginAtom()->calcExplicitValence(false);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    for bond in &mut molecule.topology_block_mut().bonds {
        if bond.order() == crate::BondOrder::Dative {
            bond.set_order(crate::BondOrder::Single);
        }
    }
    // COSMolKit's writer currently stores RDKit property-cache-derived valence
    // state in `derived_cache.valence`. Recomputing that cache from the
    // topology-mutated working copy changes the observable writer behavior for
    // former dative donors (`N->O` becomes `[NH3]O` instead of RDKit's `NO`).
    //
    // RDKit updates only the begin atom's internal explicit-valence cache here;
    // it does not rebuild a whole-molecule valence assignment for the writer
    // working copy. Preserving the existing derived valence cache matches the
    // current RDKit writer state boundary for the modeled input space.
    Ok(())
}

fn normalize_dative_bonds_for_cx_smiles(_molecule: &Molecule) -> Result<(), SmilesWriteError> {
    // In CX SMILES mode, dative bonds are preserved and written as `_Z:2:...`
    // entries in the CX extension section. No molecule mutation needed.
    Ok(())
}

fn normalize_hydrogen_bonds_for_cx_smiles(_molecule: &Molecule) -> Result<(), SmilesWriteError> {
    // In CX SMILES mode, hydrogen bonds are preserved and written as `_Z:1:...`
    // entries in the CX extension section. No molecule mutation needed.
    Ok(())
}

fn apply_cx_bond_direction_policy(
    molecule: &mut Molecule,
    restore_bond_dirs: RestoreBondDirOption,
) -> Result<(), SmilesWriteError> {
    // In CX SMILES mode, bond directions are kept for stereo output.
    // `RestoreBondDirOption::Clear` clears unknown/either directions
    // for plain SMILES, but in CX mode we preserve them since the
    // bond direction state is already stored in typed atom/bond state.
    // For `RestoreBondDirOption::True`, RDKit restores from the molblock;
    // this is a no-op here since our directions are already in typed state.
    match restore_bond_dirs {
        RestoreBondDirOption::Clear => {
            let mut changed = false;
            for bond in &mut molecule.topology_block_mut().bonds {
                if matches!(
                    bond.direction(),
                    BondDirection::Unknown | BondDirection::EitherDouble
                ) {
                    bond.set_direction(BondDirection::None);
                    changed = true;
                }
            }
            if changed {
                molecule
                    .derived_cache_mut()
                    .invalidate(crate::DerivedState::STEREO | crate::DerivedState::DRAWING);
            }
        }
        RestoreBondDirOption::True | RestoreBondDirOption::None => {
            // Keep existing direction state as-is.
        }
    }
    Ok(())
}

fn remove_plain_smiles_only_cx_state(molecule: &mut Molecule) -> Result<(), SmilesWriteError> {
    let mut changed = false;
    for bond in &mut molecule.topology_block_mut().bonds {
        if matches!(
            bond.direction(),
            BondDirection::Unknown | BondDirection::EitherDouble
        ) {
            bond.set_direction(BondDirection::None);
            changed = true;
        }
        if bond.stereo() == BondStereo::Any {
            bond.set_stereo_atoms(None);
            bond.set_stereo(BondStereo::None);
            changed = true;
        }
    }
    if changed {
        molecule
            .derived_cache_mut()
            .invalidate(crate::DerivedState::STEREO | crate::DerivedState::DRAWING);
    }
    Ok(())
}

fn validate_cx_extension_plan(fields: CxSmilesFields) -> Result<(), SmilesWriteError> {
    // All CX field types are now supported through write_cx_* functions.
    // If a specific field is requested, the corresponding writer will
    // produce output or return empty string if no data exists.
    let _ = fields;
    Ok(())
}

fn is_minimal_plain_smiles_path(params: &SmilesWriteParams) -> bool {
    !params.do_isomeric_smiles
        && !params.do_kekule
        && !params.canonical
        && !params.clean_stereo
        && !params.all_hydrogens_explicit
        && !params.do_random
        && params.include_dative_bonds
        && !params.ignore_atom_map_numbers
}

fn validate_minimal_plain_smiles_molecule(molecule: &Molecule) -> bool {
    for atom in molecule.atoms() {
        if atom.query().is_some()
            || atom.radical_electrons() != 0
            || atom.chiral_tag() != ChiralTag::Unspecified
        {
            // [deferred] Minimal plain SMILES path doesn't handle query atoms,
            // radical-bearing atoms, or chiral atoms. Falls through to the
            // standard path which supports these features via get_atom_smiles()
            // and isomeric SMILES / stereochemistry handling.
            return false;
        }
        if atom.props().keys().any(|key| {
            key != "dummyLabel" && key != "_SmilesStart" && key != "_supplementalSmilesLabel"
        }) {
            // [deferred] Non-whitelisted atom properties (e.g. map numbers,
            // custom data) are not supported by the minimal fast path.
            // Falls through to the standard path which writes atoms with
            // their full property set.
            return false;
        }
    }
    for bond in molecule.bonds() {
        if bond.direction() != BondDirection::None
            || bond.stereo() != BondStereo::None
            || bond.query().is_some()
            || !matches!(
                bond.order(),
                BondOrder::Single
                    | BondOrder::Double
                    | BondOrder::Triple
                    | BondOrder::Quadruple
                    | BondOrder::Aromatic
                    | BondOrder::Dative
            )
        {
            // [deferred] Bonds with direction, stereo, query, or non-standard
            // orders (e.g. Hydrogen, Zero) are not supported by the minimal
            // fast path. Falls through to the standard path.
            return false;
        }
    }
    true
}

fn canonical_dfs_traversal(
    molecule: &Molecule,
    plan: &FragmentWritePlan,
    start_atom: AtomId,
    ranks: &[usize],
    do_isomeric_smiles: bool,
    clean_stereo: bool,
    do_random: bool,
    bond_symbols: Option<&[String]>,
) -> Result<CanonicalTraversalResult, SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION Canon::canonicalDFSTraversal
    // RDKit✔️✔️: void canonicalDFSTraversal(ROMol &mol, int atomIdx, int inBondIdx,
    // RDKit✔️✔️:                            std::vector<AtomColors> &colors,
    // RDKit✔️✔️:                            const UINT_VECT &ranks, MolStack &molStack,
    // RDKit✔️✔️:                            VECT_INT_VECT &atomRingClosures,
    // RDKit✔️✔️:                            std::vector<INT_LIST> &atomTraversalBondOrder,
    // RDKit✔️✔️:                            const boost::dynamic_bitset<> *bondsInPlay,
    // RDKit✔️✔️:                            const std::vector<std::string> *bondSymbols,
    // RDKit✔️✔️:                            bool doRandom) {
    // RDKit✔️✔️:   dfsFindCycles(mol, atomIdx, inBondIdx, tcolors, ranks, atomRingClosures,
    // RDKit✔️✔️:                 bondsInPlay, bondSymbols, doRandom);
    // RDKit✔️✔️:   boost::dynamic_bitset<> cyclesAvailable(MAX_CYCLES);
    // RDKit✔️✔️:   cyclesAvailable.set();
    // RDKit✔️✔️:   dfsBuildStack(mol, atomIdx, inBondIdx, colors, ranks, cyclesAvailable,
    // RDKit✔️✔️:                 molStack, atomRingClosures, atomTraversalBondOrder, bondsInPlay,
    // RDKit✔️✔️:                 bondSymbols, doRandom);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Canon::canonicalDFSTraversal
    //
    // BEGIN RDKIT CPP FUNCTION Canon::dfsFindCycles / dfsBuildStack traversal sections
    // RDKit✔️✔️:   colors[atomIdx] = GREY_NODE;
    // RDKit✔️✔️:   std::vector<PossibleType> possibles;
    // RDKit✔️✔️:   std::sort(possibles.begin(), possibles.end(), _possibleCompare);
    // RDKit✔️✔️:   for (auto &possible : possibles) {
    // RDKit✔️✔️:     int possibleIdx = std::get<1>(possible);
    // RDKit✔️✔️:     Bond *bond = std::get<2>(possible);
    // RDKit✔️✔️:     switch (colors[possibleIdx]) {
    // RDKit✔️✔️:       case WHITE_NODE:
    // RDKit✔️✔️:         dfsFindCycles(mol, possibleIdx, bond->getIdx(), colors, ranks,
    // RDKit✔️✔️:                       atomRingClosures, bondsInPlay, bondSymbols, doRandom);
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case GREY_NODE:
    // RDKit✔️✔️:         atomRingClosures[possibleIdx].push_back(bond->getIdx());
    // RDKit✔️✔️:         atomRingClosures[atomIdx].push_back(bond->getIdx());
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       default:
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   colors[atomIdx] = BLACK_NODE;
    // RDKit✔️✔️:   molStack.push_back(MolStackElem(atom));
    // RDKit✔️✔️:   colors[atomIdx] = GREY_NODE;
    // RDKit✔️✔️:   if (!atomRingClosures[atomIdx].empty()) {
    // RDKit✔️✔️:     for (auto bIdx : atomRingClosures[atomIdx]) {
    // RDKit✔️✔️:       Bond *bond = mol.getBondWithIdx(bIdx);
    // RDKit✔️✔️:       if (bond->getPropIfPresent(common_properties::_TraversalRingClosureBond,
    // RDKit✔️✔️:                              ringIdx)) {
    // RDKit✔️✔️:         molStack.push_back(MolStackElem(bond, atomIdx));
    // RDKit✔️✔️:         molStack.push_back(MolStackElem(ringIdx));
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         bond->setProp(common_properties::_TraversalRingClosureBond,
    // RDKit✔️✔️:                       static_cast<unsigned int>(lowestRingIdx));
    // RDKit✔️✔️:         molStack.push_back(MolStackElem(lowestRingIdx));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (auto possiblesIt = possibles.begin(); possiblesIt != possibles.end();
    // RDKit✔️✔️:        ++possiblesIt) {
    // RDKit✔️✔️:     if (possiblesIt + 1 != possibles.end()) {
    // RDKit✔️✔️:       molStack.push_back(MolStackElem("(", rdcast<int>(possiblesIt - possibles.begin())));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     molStack.push_back(MolStackElem(bond, atomIdx));
    // RDKit✔️✔️:     dfsBuildStack(mol, possibleIdx, bond->getIdx(), colors, ranks,
    // RDKit✔️✔️:                   cyclesAvailable, molStack, atomRingClosures,
    // RDKit✔️✔️:                   atomTraversalBondOrder, bondsInPlay, bondSymbols, doRandom);
    // RDKit✔️✔️:     if (possiblesIt + 1 != possibles.end()) {
    // RDKit✔️✔️:       molStack.push_back(MolStackElem(")", rdcast<int>(possiblesIt - possibles.begin())));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION Canon::dfsFindCycles / dfsBuildStack traversal sections
    #[derive(Clone, Copy, PartialEq, Eq)]
    enum Color {
        White,
        Grey,
        Black,
    }

    fn sorted_incident_bonds(
        molecule: &Molecule,
        atom: AtomId,
        rank_by_atom: &[usize],
        colors: &[Color],
        ring_info: &crate::RingInfo,
        do_random: bool,
        bond_symbols: Option<&[String]>,
        bonds_in_play: &[bool],
    ) -> Vec<(BondId, AtomId)> {
        let mut incident = molecule
            .topology_block()
            .adjacency
            .neighbors_of(atom.index())
            .iter()
            .filter(|neighbor| bonds_in_play[neighbor.bond.index()])
            .map(|neighbor| (neighbor.bond, AtomId::new(neighbor.atom_index)))
            .map(|(bond, other)| {
                let rank = if do_random {
                    next_random_smiles_u64() as i64
                } else {
                    traversal_possible_rank(
                        molecule,
                        bond,
                        other,
                        rank_by_atom,
                        colors,
                        ring_info,
                        bond_symbols,
                    )
                };
                (rank, bond, other)
            })
            .collect::<Vec<_>>();
        incident.sort_by_key(|(rank, _, _)| *rank);
        incident
            .into_iter()
            .map(|(_, bond, other)| (bond, other))
            .collect()
    }

    fn traversal_possible_rank(
        molecule: &Molecule,
        bond: BondId,
        other: AtomId,
        rank_by_atom: &[usize],
        colors: &[Color],
        ring_info: &crate::RingInfo,
        bond_symbols: Option<&[String]>,
    ) -> i64 {
        let mut rank = rank_by_atom[other.index()] as i64;
        let bond_order_rank = rdkit_bond_order_rank(molecule.bonds()[bond.index()].order());
        if colors[other.index()] == Color::Grey {
            rank -= (CANON_MAX_BONDTYPE + 1) * CANON_MAX_NATOMS * CANON_MAX_NATOMS;
            if let Some(symbols) = bond_symbols {
                rank += i64::from(gboost_hash_range(symbols[bond.index()].as_bytes()) % 5000)
                    * CANON_MAX_NATOMS;
            } else {
                rank += (CANON_MAX_BONDTYPE - bond_order_rank) * CANON_MAX_NATOMS;
            }
        } else if ring_info.num_bond_rings(bond) > 0 {
            if let Some(symbols) = bond_symbols {
                rank += i64::from(gboost_hash_range(symbols[bond.index()].as_bytes()) % 5000)
                    * CANON_MAX_NATOMS
                    * CANON_MAX_NATOMS;
            } else {
                rank +=
                    (CANON_MAX_BONDTYPE - bond_order_rank) * CANON_MAX_NATOMS * CANON_MAX_NATOMS;
            }
        }
        rank
    }

    fn gboost_hash_range(bytes: &[u8]) -> u32 {
        let mut seed = 0u32;
        for byte in bytes {
            seed ^= u32::from(*byte)
                .wrapping_add(0x9e37_79b9)
                .wrapping_add(seed << 6)
                .wrapping_add(seed >> 2);
        }
        seed
    }

    fn rdkit_bond_order_rank(order: BondOrder) -> i64 {
        match order {
            BondOrder::Null | BondOrder::Unspecified => 0,
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Quadruple => 4,
            BondOrder::Quintuple => 5,
            BondOrder::Hextuple => 6,
            BondOrder::OneAndHalf => 7,
            BondOrder::TwoAndHalf => 8,
            BondOrder::ThreeAndHalf => 9,
            BondOrder::FourAndHalf => 10,
            BondOrder::FiveAndHalf => 11,
            BondOrder::Aromatic => 12,
            BondOrder::Ionic => 13,
            BondOrder::Hydrogen => 14,
            BondOrder::ThreeCenter => 15,
            BondOrder::DativeOne => 16,
            BondOrder::Dative => 17,
            BondOrder::DativeLeft => 18,
            BondOrder::DativeRight => 19,
            BondOrder::Other => 20,
            BondOrder::Zero => 21,
        }
    }

    fn dfs_find_cycles(
        molecule: &Molecule,
        plan: &FragmentWritePlan,
        atom: AtomId,
        parent_bond: Option<BondId>,
        do_random: bool,
        colors: &mut [Color],
        atom_ring_closures: &mut [Vec<BondId>],
        rank_by_atom: &[usize],
        ring_info: &crate::RingInfo,
        bond_symbols: Option<&[String]>,
        bonds_in_play: &[bool],
    ) {
        colors[atom.index()] = Color::Grey;
        for (bond, other) in sorted_incident_bonds(
            molecule,
            atom,
            rank_by_atom,
            colors,
            ring_info,
            do_random,
            bond_symbols,
            bonds_in_play,
        ) {
            if Some(bond) == parent_bond {
                continue;
            }
            match colors[other.index()] {
                Color::White => {
                    dfs_find_cycles(
                        molecule,
                        plan,
                        other,
                        Some(bond),
                        do_random,
                        colors,
                        atom_ring_closures,
                        rank_by_atom,
                        ring_info,
                        bond_symbols,
                        bonds_in_play,
                    );
                }
                Color::Grey => {
                    atom_ring_closures[other.index()].push(bond);
                    atom_ring_closures[atom.index()].push(bond);
                }
                Color::Black => {}
            }
        }
        colors[atom.index()] = Color::Black;
    }

    fn dfs_build_stack(
        molecule: &Molecule,
        plan: &FragmentWritePlan,
        atom: AtomId,
        parent_bond: Option<BondId>,
        do_random: bool,
        colors: &mut [Color],
        atom_ring_closures: &[Vec<BondId>],
        cycles_available: &mut [bool],
        traversal_ring_indices: &mut [Option<usize>],
        traversal_ring_closure_bonds: &mut [bool],
        rank_by_atom: &[usize],
        ring_info: &crate::RingInfo,
        bond_symbols: Option<&[String]>,
        stack: &mut Vec<MolStackElem>,
        atom_traversal_bond_order: &mut [Vec<BondId>],
        bonds_in_play: &[bool],
    ) -> Result<(), SmilesWriteError> {
        let mut seen_from_here = vec![false; molecule.num_atoms()];
        seen_from_here[atom.index()] = true;
        stack.push(MolStackElem::Atom(atom));
        colors[atom.index()] = Color::Grey;
        let mut trav_list = Vec::new();
        if let Some(parent_bond) = parent_bond {
            trav_list.push(parent_bond);
        }

        let mut rings_closed = Vec::new();
        for &bond in &atom_ring_closures[atom.index()] {
            trav_list.push(bond);
            if let Some(other) = bond_other_atom(&molecule.bonds()[bond.index()], atom) {
                seen_from_here[other.index()] = true;
            }
            if let Some(ring_idx) = traversal_ring_indices[bond.index()] {
                stack.push(MolStackElem::Bond(bond, atom));
                stack.push(MolStackElem::Ring { bond, ring_idx });
                rings_closed.push(ring_idx - 1);
            } else {
                let Some(lowest_ring_idx) =
                    cycles_available.iter().position(|available| *available)
                else {
                    return invariant_stage_error(
                        SmilesPlanStage::ShortTermBondWriter,
                        "write_ring_closure() could not allocate a free ring index",
                    );
                };
                cycles_available[lowest_ring_idx] = false;
                let ring_idx = lowest_ring_idx + 1;
                traversal_ring_indices[bond.index()] = Some(ring_idx);
                traversal_ring_closure_bonds[bond.index()] = true;
                stack.push(MolStackElem::Ring { bond, ring_idx });
            }
        }
        for ring_idx in rings_closed {
            cycles_available[ring_idx] = true;
        }

        let children = sorted_incident_bonds(
            molecule,
            atom,
            rank_by_atom,
            colors,
            ring_info,
            do_random,
            bond_symbols,
            bonds_in_play,
        )
        .into_iter()
        .filter(|(bond, other)| {
            Some(*bond) != parent_bond
                && colors[other.index()] == Color::White
                && !seen_from_here[other.index()]
        })
        .collect::<Vec<_>>();

        for (idx, (bond, child)) in children.iter().enumerate() {
            let is_branch = idx + 1 != children.len();
            if is_branch {
                stack.push(MolStackElem::BranchOpen);
            }
            stack.push(MolStackElem::Bond(*bond, atom));
            trav_list.push(*bond);
            dfs_build_stack(
                molecule,
                plan,
                *child,
                Some(*bond),
                do_random,
                colors,
                atom_ring_closures,
                cycles_available,
                traversal_ring_indices,
                traversal_ring_closure_bonds,
                rank_by_atom,
                ring_info,
                bond_symbols,
                stack,
                atom_traversal_bond_order,
                bonds_in_play,
            )?;
            if is_branch {
                stack.push(MolStackElem::BranchClose);
            }
        }
        atom_traversal_bond_order[atom.index()] = trav_list;
        colors[atom.index()] = Color::Black;
        Ok(())
    }

    let mut cycle_colors = vec![Color::White; molecule.num_atoms()];
    let mut atom_ring_closures = vec![Vec::new(); molecule.num_atoms()];
    let mut stack = Vec::new();
    let mut atom_traversal_bond_order = vec![Vec::new(); molecule.num_atoms()];
    let mut rank_by_atom = vec![usize::MAX; molecule.num_atoms()];
    let mut bonds_in_play = vec![false; molecule.num_bonds()];
    let mut traversal_ring_closure_bonds = vec![false; molecule.num_bonds()];
    for (idx, atom) in plan.atoms.iter().enumerate() {
        rank_by_atom[atom.index()] = ranks[idx];
    }
    for bond in &plan.bonds {
        bonds_in_play[bond.index()] = true;
    }
    // BEGIN RDKIT CPP FUNCTION Canon::canonicalDFSTraversal ring-info requirement
    // RDKit❗✔️: traversal only queries `numBondRings()` / `numAtomRings()` from
    // RDKit❗✔️: a `fastFindRings()`-or-better RingInfo. It does not require
    // RDKit❗✔️: symmetrized SSSR membership here.
    // END RDKIT CPP FUNCTION Canon::canonicalDFSTraversal ring-info requirement
    let ring_info = molecule
        .derived_cache()
        .rings
        .clone()
        .filter(|rings| rings.is_find_fast_or_better())
        .unwrap_or(crate::fast_find_rings(molecule)?);

    dfs_find_cycles(
        molecule,
        plan,
        start_atom,
        None,
        do_random,
        &mut cycle_colors,
        &mut atom_ring_closures,
        &rank_by_atom,
        &ring_info,
        bond_symbols,
        &bonds_in_play,
    );
    let mut stack_colors = vec![Color::White; molecule.num_atoms()];
    let mut cycles_available = vec![true; molecule.num_bonds().max(1)];
    let mut traversal_ring_indices = vec![None; molecule.num_bonds()];
    dfs_build_stack(
        molecule,
        plan,
        start_atom,
        None,
        do_random,
        &mut stack_colors,
        &atom_ring_closures,
        &mut cycles_available,
        &mut traversal_ring_indices,
        &mut traversal_ring_closure_bonds,
        &rank_by_atom,
        &ring_info,
        bond_symbols,
        &mut stack,
        &mut atom_traversal_bond_order,
        &bonds_in_play,
    )?;
    let chiral_adjustments = if do_isomeric_smiles {
        compute_writer_chiral_adjustments(
            molecule,
            plan,
            start_atom,
            &atom_ring_closures,
            &atom_traversal_bond_order,
            &stack,
            clean_stereo,
        )?
    } else {
        WriterChiralAdjustments::default()
    };

    if plan
        .atoms
        .iter()
        .any(|atom_id| stack_colors[atom_id.index()] != Color::Black)
    {
        return invariant_stage_error(
            SmilesPlanStage::ShortTermAtomWriter,
            "canonical_dfs_traversal() did not visit every atom in the fragment plan",
        );
    }

    Ok(CanonicalTraversalResult {
        stack,
        traversal_ring_closure_bonds,
        chiral_tag_overrides: chiral_adjustments.chiral_tag_overrides,
        chiral_inversions: chiral_adjustments.chiral_inversions,
        chiral_permutations: chiral_adjustments.chiral_permutations,
        broken_chiral_atoms: chiral_adjustments.broken_chiral_atoms,
    })
}

#[cfg(test)]
fn debug_atom_ring_closures_for_writer(
    molecule: &Molecule,
    plan: &FragmentWritePlan,
    start_atom: AtomId,
    ranks: &[usize],
    bond_symbols: Option<&[String]>,
) -> Result<(crate::RingInfo, Vec<Vec<BondId>>), SmilesWriteError> {
    fn gboost_hash_range_debug(bytes: &[u8]) -> u32 {
        let mut seed = 0u32;
        for byte in bytes {
            seed ^= u32::from(*byte)
                .wrapping_add(0x9e37_79b9)
                .wrapping_add(seed << 6)
                .wrapping_add(seed >> 2);
        }
        seed
    }

    fn rdkit_bond_order_rank_debug(order: BondOrder) -> i64 {
        match order {
            BondOrder::Null | BondOrder::Unspecified => 0,
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Quadruple => 4,
            BondOrder::Quintuple => 5,
            BondOrder::Hextuple => 6,
            BondOrder::OneAndHalf => 7,
            BondOrder::TwoAndHalf => 8,
            BondOrder::ThreeAndHalf => 9,
            BondOrder::FourAndHalf => 10,
            BondOrder::FiveAndHalf => 11,
            BondOrder::Aromatic => 12,
            BondOrder::Ionic => 13,
            BondOrder::Hydrogen => 14,
            BondOrder::ThreeCenter => 15,
            BondOrder::DativeOne => 16,
            BondOrder::Dative => 17,
            BondOrder::DativeLeft => 18,
            BondOrder::DativeRight => 19,
            BondOrder::Other => 20,
            BondOrder::Zero => 21,
        }
    }

    #[derive(Clone, Copy, PartialEq, Eq)]
    enum DebugColor {
        White,
        Grey,
        Black,
    }

    fn sorted_incident_bonds_debug(
        molecule: &Molecule,
        atom: AtomId,
        rank_by_atom: &[usize],
        colors: &[DebugColor],
        ring_info: &crate::RingInfo,
        bond_symbols: Option<&[String]>,
        bonds_in_play: &[bool],
    ) -> Vec<(BondId, AtomId)> {
        let mut incident = molecule
            .topology_block()
            .adjacency
            .neighbors_of(atom.index())
            .iter()
            .filter(|neighbor| bonds_in_play[neighbor.bond.index()])
            .map(|neighbor| (neighbor.bond, AtomId::new(neighbor.atom_index)))
            .map(|(bond, other)| {
                let mut rank = rank_by_atom[other.index()] as i64;
                let bond_order_rank =
                    rdkit_bond_order_rank_debug(molecule.bonds()[bond.index()].order());
                if colors[other.index()] == DebugColor::Grey {
                    rank -= (CANON_MAX_BONDTYPE + 1) * CANON_MAX_NATOMS * CANON_MAX_NATOMS;
                    if let Some(symbols) = bond_symbols {
                        rank += i64::from(
                            gboost_hash_range_debug(symbols[bond.index()].as_bytes()) % 5000,
                        ) * CANON_MAX_NATOMS;
                    } else {
                        rank += (CANON_MAX_BONDTYPE - bond_order_rank) * CANON_MAX_NATOMS;
                    }
                } else if ring_info.num_bond_rings(bond) > 0 {
                    if let Some(symbols) = bond_symbols {
                        rank += i64::from(
                            gboost_hash_range_debug(symbols[bond.index()].as_bytes()) % 5000,
                        ) * CANON_MAX_NATOMS
                            * CANON_MAX_NATOMS;
                    } else {
                        rank += (CANON_MAX_BONDTYPE - bond_order_rank)
                            * CANON_MAX_NATOMS
                            * CANON_MAX_NATOMS;
                    }
                }
                (rank, bond, other)
            })
            .collect::<Vec<_>>();
        incident.sort_by_key(|(rank, _, _)| *rank);
        incident
            .into_iter()
            .map(|(_, bond, other)| (bond, other))
            .collect()
    }

    fn dfs_find_cycles_debug(
        molecule: &Molecule,
        atom: AtomId,
        parent_bond: Option<BondId>,
        colors: &mut [DebugColor],
        atom_ring_closures: &mut [Vec<BondId>],
        rank_by_atom: &[usize],
        ring_info: &crate::RingInfo,
        bond_symbols: Option<&[String]>,
        bonds_in_play: &[bool],
    ) {
        colors[atom.index()] = DebugColor::Grey;
        let possibles = sorted_incident_bonds_debug(
            molecule,
            atom,
            rank_by_atom,
            colors,
            ring_info,
            bond_symbols,
            bonds_in_play,
        );
        eprint!("debug_dfsFindCycles atom={} possibles", atom.index());
        for (bond, other) in &possibles {
            let mut rank = rank_by_atom[other.index()] as i64;
            let bond_order_rank =
                rdkit_bond_order_rank_debug(molecule.bonds()[bond.index()].order());
            if colors[other.index()] == DebugColor::Grey {
                rank -= (CANON_MAX_BONDTYPE + 1) * CANON_MAX_NATOMS * CANON_MAX_NATOMS;
                if let Some(symbols) = bond_symbols {
                    rank +=
                        i64::from(gboost_hash_range_debug(symbols[bond.index()].as_bytes()) % 5000)
                            * CANON_MAX_NATOMS;
                } else {
                    rank += (CANON_MAX_BONDTYPE - bond_order_rank) * CANON_MAX_NATOMS;
                }
            } else if ring_info.num_bond_rings(*bond) > 0 {
                if let Some(symbols) = bond_symbols {
                    rank +=
                        i64::from(gboost_hash_range_debug(symbols[bond.index()].as_bytes()) % 5000)
                            * CANON_MAX_NATOMS
                            * CANON_MAX_NATOMS;
                } else {
                    rank += (CANON_MAX_BONDTYPE - bond_order_rank)
                        * CANON_MAX_NATOMS
                        * CANON_MAX_NATOMS;
                }
            }
            eprint!(" ({},r={},b={})", other.index(), rank, bond.index());
        }
        eprintln!();
        for (bond, other) in possibles {
            if Some(bond) == parent_bond {
                continue;
            }
            match colors[other.index()] {
                DebugColor::White => dfs_find_cycles_debug(
                    molecule,
                    other,
                    Some(bond),
                    colors,
                    atom_ring_closures,
                    rank_by_atom,
                    ring_info,
                    bond_symbols,
                    bonds_in_play,
                ),
                DebugColor::Grey => {
                    atom_ring_closures[other.index()].push(bond);
                    atom_ring_closures[atom.index()].push(bond);
                }
                DebugColor::Black => {}
            }
        }
        colors[atom.index()] = DebugColor::Black;
    }

    let mut rank_by_atom = vec![usize::MAX; molecule.num_atoms()];
    for (idx, atom) in plan.atoms.iter().enumerate() {
        rank_by_atom[atom.index()] = ranks[idx];
    }
    let bonds_in_play = {
        let mut value = vec![false; molecule.num_bonds()];
        for bond in &plan.bonds {
            value[bond.index()] = true;
        }
        value
    };
    let ring_info = crate::fast_find_rings(molecule)?;
    let mut colors = vec![DebugColor::White; molecule.num_atoms()];
    let mut atom_ring_closures = vec![Vec::new(); molecule.num_atoms()];
    dfs_find_cycles_debug(
        molecule,
        start_atom,
        None,
        &mut colors,
        &mut atom_ring_closures,
        &rank_by_atom,
        &ring_info,
        bond_symbols,
        &bonds_in_play,
    );
    Ok((ring_info, atom_ring_closures))
}

#[derive(Debug, Default)]
struct WriterChiralAdjustments {
    chiral_tag_overrides: BTreeMap<AtomId, ChiralTag>,
    chiral_inversions: BTreeSet<AtomId>,
    chiral_permutations: BTreeMap<AtomId, u32>,
    broken_chiral_atoms: BTreeSet<AtomId>,
}

fn compute_writer_chiral_adjustments(
    molecule: &Molecule,
    plan: &FragmentWritePlan,
    start_atom: AtomId,
    atom_ring_closures: &[Vec<BondId>],
    atom_traversal_bond_order: &[Vec<BondId>],
    stack: &[MolStackElem],
    _clean_stereo: bool,
) -> Result<WriterChiralAdjustments, SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION Canon::canonicalizeFragment chiral traversal section
    // RDKit✔️✔️: std::vector<INT_LIST> atomTraversalBondOrder(nAtoms);
    // RDKit✔️✔️: std::vector<int> atomPermutationIndices(nAtoms, 0);
    // RDKit✔️✔️: const INT_LIST &trueOrder = atomTraversalBondOrder[atom->getIdx()];
    // RDKit✔️✔️: if (trueOrder.size() < atom->getDegree()) { ... append missing bonds ... }
    // RDKit✔️✔️: if (!perm) { nSwaps = atom->getPerturbationOrder(tOrder); }
    // RDKit✔️✔️: else { insertImplicitNbors(tOrder, atom->getChiralTag(), firstInPart);
    // RDKit✔️✔️:        perm = Chirality::getChiralPermutation(atom, tOrder); }
    // RDKit✔️✔️: if (doChiralInversions && chiralAtomNeedsTagInversion(...)) { ++nSwaps; }
    // RDKit✔️✔️: if (nSwaps % 2) { numSwapsChiralAtoms.set(atom->getIdx()); }
    // RDKit✔️✔️: atomPermutationIndices[atom->getIdx()] = perm;
    // END RDKIT CPP FUNCTION Canon::canonicalizeFragment chiral traversal section
    // BEGIN RDKIT CPP FUNCTION Canon::canonicalizeFragment chiral post-processing section
    // RDKit✔️✔️: std::vector<unsigned int> atomVisitOrders(mol.getNumAtoms());
    // RDKit✔️✔️: for (const auto &msI : molStack) { if (msI.type == MOL_STACK_ATOM) atomVisitOrders[...] = pos; ++pos; }
    // RDKit✔️✔️: boost::dynamic_bitset<> ringStereoChemAdjusted(nAtoms);
    // RDKit✔️✔️: if (msI.obj.atom->hasProp(common_properties::_ringStereoAtoms)) { ... setChiralTag(CCW) ... }
    // RDKit✔️✔️: else if (msI.obj.atom->getPropIfPresent("_stereoGroup", sgidx) && mol.getStereoGroups().size() > sgidx) { ... }
    // RDKit✔️✔️: else { if (numSwapsChiralAtoms[idx]) atom->invertChirality(); else if (atomPermutationIndices[idx]) atom->setProp(_chiralPermutation, ...); }
    // END RDKIT CPP FUNCTION Canon::canonicalizeFragment chiral post-processing section
    let plan_atoms = plan.atoms.iter().copied().collect::<BTreeSet<_>>();
    let plan_bonds = plan.bonds.iter().copied().collect::<BTreeSet<_>>();
    let mut adjustments = WriterChiralAdjustments::default();
    let mut num_swaps_chiral_atoms = vec![false; molecule.num_atoms()];
    let mut atom_permutation_indices = vec![0_u32; molecule.num_atoms()];

    for atom_id in &plan.atoms {
        let atom = &molecule.atoms()[atom_id.index()];
        if atom.chiral_tag() == ChiralTag::Unspecified {
            continue;
        }

        let incident = incident_bonds(molecule, *atom_id);
        if incident.iter().any(|bond| !plan_bonds.contains(bond)) {
            adjustments.broken_chiral_atoms.insert(*atom_id);
            continue;
        }

        let mut true_order = atom_traversal_bond_order[atom_id.index()].clone();
        if true_order.len() < incident.len() {
            for bond in &incident {
                if !true_order.contains(bond) {
                    true_order.push(*bond);
                }
            }
        }

        match atom.chiral_tag() {
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw => {
                let mut swaps = count_swaps_to_interconvert(&true_order, &incident)?;
                if chiral_atom_needs_tag_inversion_for_writer(
                    molecule,
                    *atom_id,
                    *atom_id == start_atom,
                    atom_ring_closures[atom_id.index()].len(),
                ) {
                    swaps = swaps.saturating_add(1);
                }
                if swaps % 2 == 1 {
                    num_swaps_chiral_atoms[atom_id.index()] = true;
                }
            }
            ChiralTag::SquarePlanar | ChiralTag::TrigonalBipyramidal | ChiralTag::Octahedral => {
                let mut probe = true_order.into_iter().map(Some).collect::<Vec<_>>();
                insert_implicit_nontetrahedral_neighbors_for_writer(
                    &mut probe,
                    atom.chiral_tag(),
                    *atom_id == start_atom,
                );
                let permutation = super::smiles::nontetrahedral_chiral_permutation_for_probe(
                    molecule, *atom_id, &probe, false,
                )
                .map_err(|_| {
                    crate::UnsupportedFeatureError::from_spec(&crate::SMILES_WRITE_FEATURE)
                })?;
                if permutation != 0 {
                    atom_permutation_indices[atom_id.index()] = permutation;
                }
            }
            _ => {}
        }
    }

    let stereo_group_by_atom = molecule
        .stereo_groups()
        .iter()
        .enumerate()
        .flat_map(|(group_idx, group)| {
            group
                .atoms()
                .iter()
                .copied()
                .map(move |atom_id| (atom_id, group_idx))
        })
        .collect::<BTreeMap<_, _>>();
    let mut atom_visit_orders = vec![usize::MAX; molecule.num_atoms()];
    for (pos, element) in stack.iter().enumerate() {
        if let MolStackElem::Atom(atom_id) = element {
            atom_visit_orders[atom_id.index()] = pos;
        }
    }
    let mut ring_stereo_chem_adjusted = vec![false; molecule.num_atoms()];
    for element in stack {
        let MolStackElem::Atom(atom_id) = element else {
            continue;
        };
        let atom = &molecule.atoms()[atom_id.index()];
        if atom.chiral_tag() == ChiralTag::Unspecified
            || adjustments.broken_chiral_atoms.contains(atom_id)
        {
            continue;
        }

        let ring_stereo_atoms = writer_ring_stereo_atoms(atom, *atom_id)?;
        if let Some(ring_stereo_atoms) = ring_stereo_atoms {
            if !ring_stereo_chem_adjusted[atom_id.index()] {
                adjustments
                    .chiral_tag_overrides
                    .insert(*atom_id, ChiralTag::TetrahedralCcw);
                ring_stereo_chem_adjusted[atom_id.index()] = true;
            }
            let source_tag = adjustments
                .chiral_tag_overrides
                .get(atom_id)
                .copied()
                .unwrap_or_else(|| atom.chiral_tag());
            for (same_orientation, neighbor_idx) in ring_stereo_atoms {
                let neighbor = AtomId::new(neighbor_idx);
                if !plan_atoms.contains(&neighbor)
                    || ring_stereo_chem_adjusted[neighbor_idx]
                    || atom_visit_orders[neighbor_idx] <= atom_visit_orders[atom_id.index()]
                {
                    continue;
                }
                let mut neighbor_tag = source_tag;
                if !same_orientation {
                    neighbor_tag = invert_tetrahedral_chiral_tag(neighbor_tag);
                }
                if num_swaps_chiral_atoms[atom_id.index()] {
                    if !num_swaps_chiral_atoms[neighbor_idx] {
                        neighbor_tag = invert_tetrahedral_chiral_tag(neighbor_tag);
                    }
                } else if num_swaps_chiral_atoms[neighbor_idx] {
                    neighbor_tag = invert_tetrahedral_chiral_tag(neighbor_tag);
                }
                adjustments
                    .chiral_tag_overrides
                    .insert(neighbor, neighbor_tag);
                ring_stereo_chem_adjusted[neighbor_idx] = true;
            }
        } else if let Some(group_idx) = stereo_group_by_atom.get(atom_id).copied() {
            if let Some(group) = molecule.stereo_groups().get(group_idx) {
                let swap_it = matches!(
                    adjustments
                        .chiral_tag_overrides
                        .get(atom_id)
                        .copied()
                        .unwrap_or_else(|| atom.chiral_tag()),
                    ChiralTag::TetrahedralCw
                );
                if swap_it {
                    let current = adjustments
                        .chiral_tag_overrides
                        .get(atom_id)
                        .copied()
                        .unwrap_or_else(|| atom.chiral_tag());
                    adjustments
                        .chiral_tag_overrides
                        .insert(*atom_id, invert_tetrahedral_chiral_tag(current));
                }
                if swap_it || num_swaps_chiral_atoms[atom_id.index()] {
                    for other_atom in group.atoms() {
                        if *other_atom == *atom_id || !plan_atoms.contains(other_atom) {
                            continue;
                        }
                        let current = adjustments
                            .chiral_tag_overrides
                            .get(other_atom)
                            .copied()
                            .unwrap_or_else(|| molecule.atoms()[other_atom.index()].chiral_tag());
                        adjustments
                            .chiral_tag_overrides
                            .insert(*other_atom, invert_tetrahedral_chiral_tag(current));
                    }
                }
            }
        } else {
            match atom.chiral_tag() {
                ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw => {
                    if num_swaps_chiral_atoms[atom_id.index()] {
                        adjustments.chiral_inversions.insert(*atom_id);
                    }
                }
                ChiralTag::SquarePlanar
                | ChiralTag::TrigonalBipyramidal
                | ChiralTag::Octahedral => {
                    let permutation = atom_permutation_indices[atom_id.index()];
                    if permutation != 0 {
                        adjustments
                            .chiral_permutations
                            .insert(*atom_id, permutation);
                    }
                }
                _ => {}
            }
        }
    }
    Ok(adjustments)
}

fn invert_tetrahedral_chiral_tag(chiral_tag: ChiralTag) -> ChiralTag {
    match chiral_tag {
        ChiralTag::TetrahedralCw => ChiralTag::TetrahedralCcw,
        ChiralTag::TetrahedralCcw => ChiralTag::TetrahedralCw,
        other => other,
    }
}

fn insert_implicit_nontetrahedral_neighbors_for_writer(
    bonds: &mut Vec<Option<BondId>>,
    chiral_tag: ChiralTag,
    first_atom: bool,
) {
    // BEGIN RDKIT CPP FUNCTION Canon::insertImplicitNbors
    // RDKit✔️✔️: static void insertImplicitNbors(INT_LIST &bonds, const Atom::ChiralType tag,
    // RDKit✔️✔️:                                 const bool firstAtom) {
    // RDKit✔️✔️:   unsigned int ref_max = Chirality::getMaxNbors(tag);
    // RDKit✔️✔️:   if (bonds.size() < ref_max) {
    // RDKit✔️✔️:     if (firstAtom) {
    // RDKit✔️✔️:       bonds.insert(bonds.begin(), ref_max - bonds.size(), -1);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       bonds.insert(++bonds.begin(), ref_max - bonds.size(), -1);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Canon::insertImplicitNbors
    let Some(ref_max) = super::smiles::nontetrahedral_max_neighbors(chiral_tag) else {
        return;
    };
    if bonds.len() < ref_max {
        let missing = ref_max - bonds.len();
        let insert_at = if first_atom { 0 } else { 1.min(bonds.len()) };
        bonds.splice(insert_at..insert_at, std::iter::repeat_n(None, missing));
    }
}

fn incident_bonds(molecule: &Molecule, atom: AtomId) -> Vec<BondId> {
    // BEGIN RDKIT CPP FUNCTION ROMol::atomBonds
    // RDKit✔️✔️: CXXBondIterator<const MolGraph, Bond *const, MolGraph::out_edge_iterator>
    // RDKit✔️✔️: atomBonds(Atom const *at) const {
    // RDKit✔️✔️:   auto pr = getAtomBonds(at);
    // RDKit✔️✔️:   return {&d_graph, pr.first, pr.second};
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION ROMol::atomBonds
    molecule
        .topology_block()
        .adjacency
        .neighbors_of(atom.index())
        .iter()
        .map(|neighbor| neighbor.bond)
        .collect()
}

fn count_swaps_to_interconvert(
    probe: &[BondId],
    reference: &[BondId],
) -> Result<usize, SmilesWriteError> {
    if probe.len() != reference.len() {
        return invariant_stage_error(
            SmilesPlanStage::ShortTermAtomWriter,
            "count_swaps_to_interconvert() probe/reference length mismatch",
        );
    }
    let mut work = probe.to_vec();
    let mut swaps = 0usize;
    for (idx, expected) in reference.iter().copied().enumerate() {
        if work[idx] == expected {
            continue;
        }
        let Some(found_idx) = work[idx..]
            .iter()
            .position(|bond| *bond == expected)
            .map(|offset| idx + offset)
        else {
            return invariant_stage_error(
                SmilesPlanStage::ShortTermAtomWriter,
                "count_swaps_to_interconvert() expected bond missing from probe order",
            );
        };
        work.swap(idx, found_idx);
        swaps = swaps.saturating_add(1);
    }
    Ok(swaps)
}

fn chiral_atom_needs_tag_inversion_for_writer(
    molecule: &Molecule,
    atom_id: AtomId,
    is_atom_first: bool,
    num_closures: usize,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION Canon::chiralAtomNeedsTagInversion
    // RDKit✔️✔️: bool chiralAtomNeedsTagInversion(const RDKit::ROMol &mol,
    // RDKit✔️✔️:                                  const RDKit::Atom *atom, bool isAtomFirst,
    // RDKit✔️✔️:                                  size_t numClosures) {
    // RDKit✔️✔️:   return atom->getDegree() == 3 &&
    // RDKit✔️✔️:          ((isAtomFirst && atom->getNumExplicitHs() == 1) ||
    // RDKit✔️✔️:           (!details::atomHasFourthValence(atom) && numClosures == 1 &&
    // RDKit✔️✔️:            !details::isUnsaturated(atom, mol)));
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Canon::chiralAtomNeedsTagInversion
    let atom = &molecule.atoms()[atom_id.index()];
    let degree = incident_bonds(molecule, atom_id).len();
    degree == 3
        && ((is_atom_first && atom.explicit_hydrogens() == 1)
            || (!atom_has_fourth_valence_for_writer(molecule, atom_id)
                && num_closures == 1
                && !atom_is_unsaturated_for_writer(molecule, atom_id)))
}

fn atom_has_fourth_valence_for_writer(molecule: &Molecule, atom_id: AtomId) -> bool {
    // BEGIN RDKIT CPP FUNCTION Canon::details::atomHasFourthValence
    // RDKit✔️✔️: bool atomHasFourthValence(const Atom *atom) {
    // RDKit✔️✔️:   if (atom->getNumExplicitHs() == 1 ||
    // RDKit✔️✔️:       (!atom->needsUpdatePropertyCache() &&
    // RDKit✔️✔️:        atom->getValence(Atom::ValenceType::IMPLICIT) == 1)) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (atom->hasQuery()) {
    // RDKit✔️✔️:     return hasSingleHQuery(atom->getQuery());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Canon::details::atomHasFourthValence
    //
    // COSMolKit does not model RDKit's per-atom property-cache dirty bit.
    // The implicit-H branch therefore relies on the current derived valence
    // cache when present, matching the writer's current modeled state space.
    let atom = &molecule.atoms()[atom_id.index()];
    atom.explicit_hydrogens() == 1
        || molecule
            .derived_cache()
            .valence
            .as_ref()
            .and_then(|valence| valence.implicit_hydrogens.get(atom_id.index()))
            .is_some_and(|implicit_hydrogens| *implicit_hydrogens == 1)
        || atom
            .query()
            .is_some_and(atom_query_has_single_h_count_for_writer)
}

fn atom_is_unsaturated_for_writer(molecule: &Molecule, atom_id: AtomId) -> bool {
    molecule
        .bonds()
        .iter()
        .filter(|bond| bond.begin() == atom_id || bond.end() == atom_id)
        .any(|bond| bond_order_as_double_for_writer(bond.order()) > 1.0)
}

fn bond_order_as_double_for_writer(order: BondOrder) -> f64 {
    match order {
        BondOrder::Single
        | BondOrder::Dative
        | BondOrder::DativeOne
        | BondOrder::DativeLeft
        | BondOrder::DativeRight => 1.0,
        BondOrder::Double => 2.0,
        BondOrder::Triple => 3.0,
        BondOrder::Quadruple => 4.0,
        BondOrder::Quintuple => 5.0,
        BondOrder::Hextuple => 6.0,
        BondOrder::OneAndHalf | BondOrder::Aromatic => 1.5,
        BondOrder::TwoAndHalf => 2.5,
        BondOrder::ThreeAndHalf => 3.5,
        BondOrder::FourAndHalf => 4.5,
        BondOrder::FiveAndHalf => 5.5,
        BondOrder::Null
        | BondOrder::Ionic
        | BondOrder::Hydrogen
        | BondOrder::ThreeCenter
        | BondOrder::Other
        | BondOrder::Zero
        | BondOrder::Unspecified => 0.0,
    }
}

fn atom_query_has_single_h_count_for_writer(query: &QueryNode<AtomQueryPredicate>) -> bool {
    // BEGIN RDKIT CPP FUNCTION Canon::details::hasSingleHQuery
    // RDKit✔️✔️: bool hasSingleHQuery(const Atom::QUERYATOM_QUERY *q) {
    // RDKit✔️✔️:   PRECONDITION(q, "bad query");
    // RDKit✔️✔️:   bool res = false;
    // RDKit✔️✔️:   const auto &descr = q->getDescription();
    // RDKit✔️✔️:   if (descr == "AtomAnd") {
    // RDKit✔️✔️:     for (auto cIt = q->beginChildren(); cIt != q->endChildren(); ++cIt) {
    // RDKit✔️✔️:       const auto &cDescr = (*cIt)->getDescription();
    // RDKit✔️✔️:       if (cDescr == "AtomHCount") {
    // RDKit✔️✔️:         return !(*cIt)->getNegation() &&
    // RDKit✔️✔️:                ((ATOM_EQUALS_QUERY *)(*cIt).get())->getVal() == 1;
    // RDKit✔️✔️:       } else if (cDescr == "AtomAnd") {
    // RDKit✔️✔️:         res = hasSingleHQuery((*cIt).get());
    // RDKit✔️✔️:         if (res) {
    // RDKit✔️✔️:           return true;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Canon::details::hasSingleHQuery
    match query {
        QueryNode::And(children) => {
            let mut found_nested = false;
            for child in children {
                match child {
                    QueryNode::Predicate(AtomQueryPredicate::ImplicitHydrogenCount(count))
                        if *count == 1 =>
                    {
                        return true;
                    }
                    QueryNode::And(_) if atom_query_has_single_h_count_for_writer(child) => {
                        found_nested = true;
                        break;
                    }
                    _ => {}
                }
            }
            found_nested
        }
        _ => false,
    }
}

fn canonicalize_double_bond_directions_for_writer(
    molecule: &mut Molecule,
    stack: &[MolStackElem],
    traversal_ring_closure_bonds: &[bool],
) -> Result<(), SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION Canon::canonicalizeFragment double-bond direction section
    // RDKit✔️✔️: std::vector<unsigned int> atomVisitOrders(mol.getNumAtoms());
    // RDKit✔️✔️: std::vector<unsigned int> bondVisitOrders(mol.getNumBonds());
    // RDKit✔️✔️: for (const auto &msI : molStack) {
    // RDKit✔️✔️:   if (msI.type == MOL_STACK_ATOM) { atomVisitOrders[msI.obj.atom->getIdx()] = pos; }
    // RDKit✔️✔️:   else if (msI.type == MOL_STACK_BOND) {
    // RDKit✔️✔️:     bondVisitOrders[msI.obj.bond->getIdx()] = pos;
    // RDKit✔️✔️:     auto dir = msI.obj.bond->getBondDir();
    // RDKit✔️✔️:     if (dir == Bond::ENDDOWNRIGHT || dir == Bond::ENDUPRIGHT) {
    // RDKit✔️✔️:       msI.obj.bond->setBondDir(Bond::NONE);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: canonicalizeDoubleBonds(mol, bondVisitOrders, atomVisitOrders,
    // RDKit✔️✔️:                         bondDirCounts, atomDirCounts, molStack);
    // RDKit✔️✔️: Canon::removeUnwantedBondDirSpecs(mol, molStack, bondDirCounts,
    // RDKit✔️✔️:                                     atomDirCounts, bondVisitOrders);
    // RDKit✔️✔️: Canon::removeRedundantBondDirSpecs(mol, molStack, bondDirCounts,
    // RDKit✔️✔️:                                    atomDirCounts);
    // END RDKIT CPP FUNCTION Canon::canonicalizeFragment double-bond direction section
    let mut atom_visit_orders = vec![usize::MAX; molecule.num_atoms()];
    let mut bond_visit_orders = vec![usize::MAX; molecule.num_bonds()];
    for (pos, item) in stack.iter().enumerate() {
        match *item {
            MolStackElem::Atom(atom) => atom_visit_orders[atom.index()] = pos,
            MolStackElem::Bond(bond, _) => {
                bond_visit_orders[bond.index()] = pos;
                if matches!(
                    molecule.bonds()[bond.index()].direction(),
                    BondDirection::EndDownRight | BondDirection::EndUpRight
                ) {
                    molecule.topology_block_mut().bonds[bond.index()]
                        .set_direction(BondDirection::None);
                }
            }
            MolStackElem::Ring { .. } => {}
            MolStackElem::BranchOpen | MolStackElem::BranchClose => {}
        }
    }

    let mut bond_dir_counts = vec![0i8; molecule.num_bonds()];
    let mut atom_dir_counts = vec![0i8; molecule.num_atoms()];
    let cip_ranks = crate::stereo::assign_atom_cip_ranks(molecule).ok();
    canonicalize_double_bonds_for_writer(
        molecule,
        &bond_visit_orders,
        &atom_visit_orders,
        &traversal_ring_closure_bonds,
        &mut bond_dir_counts,
        &mut atom_dir_counts,
        stack,
        cip_ranks.as_deref(),
    );
    remove_unwanted_bond_dir_specs_for_writer(
        molecule,
        stack,
        &mut bond_dir_counts,
        &mut atom_dir_counts,
        &bond_visit_orders,
    );
    remove_redundant_bond_dir_specs_for_writer(
        molecule,
        stack,
        &mut bond_dir_counts,
        &mut atom_dir_counts,
    );
    molecule
        .derived_cache_mut()
        .invalidate(crate::DerivedState::STEREO | crate::DerivedState::DRAWING);
    Ok(())
}

fn canonicalize_double_bonds_for_writer(
    molecule: &mut Molecule,
    bond_visit_orders: &[usize],
    atom_visit_orders: &[usize],
    traversal_ring_closure_bonds: &[bool],
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
    stack: &[MolStackElem],
    cip_ranks: Option<&[u32]>,
) {
    // BEGIN RDKIT CPP FUNCTION Canon::canonicalizeDoubleBonds
    // RDKit✔️✔️: for (auto &msI : molStack) {
    // RDKit✔️✔️:   if (msI.type != MOL_STACK_BOND) { continue; }
    // RDKit✔️✔️:   if (bond->getBondType() != Bond::DOUBLE ||
    // RDKit✔️✔️:       bond->getStereo() <= Bond::STEREOANY ||
    // RDKit✔️✔️:       bond->getStereoAtoms().size() < 2) {
    // RDKit✔️✔️:     bond->setStereo(Bond::STEREONONE);
    // RDKit✔️✔️:     continue;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ... prioritize by neighboring stereo bonds and molStack order ...
    // RDKit✔️✔️: }
    // RDKit✔️✔️: while (!q.empty()) { Canon::canonicalizeDoubleBond(...); }
    // END RDKIT CPP FUNCTION Canon::canonicalizeDoubleBonds
    let mut stereo_bond_neighbors = BTreeMap::<BondId, Vec<BondId>>::new();
    let mut candidates = Vec::new();
    for item in stack {
        let MolStackElem::Bond(bond, _) = *item else {
            continue;
        };
        let bond_ref = &molecule.bonds()[bond.index()];
        if !is_writer_stereo_double_bond(bond_ref)
            || !writer_double_bond_has_distinguishable_ends(molecule, bond, cip_ranks)
        {
            if bond_ref.order() == BondOrder::Double {
                let bond_mut = &mut molecule.topology_block_mut().bonds[bond.index()];
                bond_mut.set_stereo_atoms(None);
                bond_mut.set_stereo(BondStereo::None);
            }
            continue;
        }
        let mut stereo_nbrs = neighboring_stereo_double_bonds_for_writer(molecule, bond);
        stereo_nbrs.sort_by_key(|neighbor| bond_visit_orders[neighbor.index()]);
        stereo_bond_neighbors.insert(bond, stereo_nbrs.clone());
        candidates.push((
            usize::MAX - stereo_nbrs.len(),
            bond_visit_orders[bond.index()],
            bond,
        ));
    }
    candidates.sort_by_key(|(stereo_rank, visit_order, _)| (*stereo_rank, *visit_order));

    let mut seen = vec![false; molecule.num_bonds()];
    for (_, _, start_bond) in candidates {
        if seen[start_bond.index()] {
            continue;
        }
        let mut queue = std::collections::VecDeque::from([start_bond]);
        while let Some(bond) = queue.pop_front() {
            if seen[bond.index()] {
                continue;
            }
            canonicalize_double_bond_for_writer(
                molecule,
                bond,
                bond_visit_orders,
                atom_visit_orders,
                traversal_ring_closure_bonds,
                bond_dir_counts,
                atom_dir_counts,
            );
            seen[bond.index()] = true;
            for &nbr in stereo_bond_neighbors.get(&bond).into_iter().flatten() {
                if !seen[nbr.index()] {
                    queue.push_back(nbr);
                }
            }
        }
    }
}

fn canonicalize_double_bond_for_writer(
    molecule: &mut Molecule,
    dbl_bond: BondId,
    bond_visit_orders: &[usize],
    atom_visit_orders: &[usize],
    traversal_ring_closure_bonds: &[bool],
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
) {
    // BEGIN RDKIT CPP FUNCTION Canon::canonicalizeDoubleBond
    // RDKit✔️✔️: void canonicalizeDoubleBond(Bond *dblBond,
    // RDKit✔️✔️:   const UINT_VECT &bondVisitOrders, const UINT_VECT &atomVisitOrders,
    // RDKit✔️✔️:   std::vector<int8_t> &bondDirCounts,
    // RDKit✔️✔️:   std::vector<int8_t> &atomDirCounts) {
    // RDKit✔️✔️:   Atom *atom1 = dblBond->getBeginAtom();
    // RDKit✔️✔️:   Atom *atom2 = dblBond->getEndAtom();
    // RDKit✔️✔️:   if ((atom1->getDegree() != 2 && atom1->getDegree() != 3) ||
    // RDKit✔️✔️:       (atom2->getDegree() != 2 && atom2->getDegree() != 3)) { return; }
    // RDKit✔️✔️:   if (atomVisitOrders[dblBond->getBeginAtomIdx()] >=
    // RDKit✔️✔️:       atomVisitOrders[dblBond->getEndAtomIdx()]) { std::swap(atom1, atom2); }
    // RDKit✔️✔️:   ... find lowest visit order neighbor bonds from both ends ...
    // RDKit✔️✔️:   bool isFirstFromAtom1Flipped = [&]() {
    // RDKit✔️✔️:     auto anchorIdx = firstFromAtom1->getOtherAtom(atom1)->getIdx();
    // RDKit✔️✔️:     return (atomVisitOrders[atom1->getIdx()] < atomVisitOrders[anchorIdx]) !=
    // RDKit✔️✔️:            firstFromAtom1->hasProp(
    // RDKit✔️✔️:                common_properties::_TraversalRingClosureBond);
    // RDKit✔️✔️:   }();
    // RDKit✔️✔️:   bool isSecondFromAtom1Flipped = [&]() {
    // RDKit✔️✔️:     if (secondFromAtom1 == nullptr) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     auto anchorIdx = secondFromAtom1->getOtherAtom(atom1)->getIdx();
    // RDKit✔️✔️:     return (atomVisitOrders[atom1->getIdx()] < atomVisitOrders[anchorIdx]) !=
    // RDKit✔️✔️:            secondFromAtom1->hasProp(
    // RDKit✔️✔️:                common_properties::_TraversalRingClosureBond);
    // RDKit✔️✔️:   }();
    // RDKit✔️✔️:   bool isFirstFromAtom2Flipped = [&]() {
    // RDKit✔️✔️:     auto anchorIdx = firstFromAtom2->getOtherAtom(atom2)->getIdx();
    // RDKit✔️✔️:     return (atomVisitOrders[anchorIdx] < atomVisitOrders[atom2->getIdx()]) !=
    // RDKit✔️✔️:            firstFromAtom2->hasProp(
    // RDKit✔️✔️:                common_properties::_TraversalRingClosureBond);
    // RDKit✔️✔️:   }();
    // RDKit✔️✔️:   bool isSecondFromAtom2Flipped = [&]() {
    // RDKit✔️✔️:     if (secondFromAtom2 == nullptr) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     auto anchorIdx = secondFromAtom2->getOtherAtom(atom2)->getIdx();
    // RDKit✔️✔️:     return (atomVisitOrders[anchorIdx] < atomVisitOrders[atom2->getIdx()]) !=
    // RDKit✔️✔️:            secondFromAtom2->hasProp(
    // RDKit✔️✔️:                common_properties::_TraversalRingClosureBond);
    // RDKit✔️✔️:   }();
    // RDKit✔️✔️:   ... set arbitrary ENDUPRIGHT if neither side is constrained ...
    // RDKit✔️✔️:   ... propagate direction with getReferenceDirection ...
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Canon::canonicalizeDoubleBond
    if !is_writer_stereo_double_bond(&molecule.bonds()[dbl_bond.index()]) {
        return;
    }
    let mut atom1 = molecule.bonds()[dbl_bond.index()].begin();
    let mut atom2 = molecule.bonds()[dbl_bond.index()].end();
    let deg1 = incident_bonds(molecule, atom1).len();
    let deg2 = incident_bonds(molecule, atom2).len();
    if !(deg1 == 2 || deg1 == 3) || !(deg2 == 2 || deg2 == 3) {
        return;
    }
    if atom_visit_orders[atom1.index()] >= atom_visit_orders[atom2.index()] {
        std::mem::swap(&mut atom1, &mut atom2);
    }

    let (first1, second1, dir1_set) = find_double_bond_neighbor_bonds_for_writer(
        molecule,
        dbl_bond,
        atom1,
        bond_visit_orders,
        bond_dir_counts,
    );
    let (first2, second2, dir2_set) = find_double_bond_neighbor_bonds_for_writer(
        molecule,
        dbl_bond,
        atom2,
        bond_visit_orders,
        bond_dir_counts,
    );
    let (Some(first1), Some(first2)) = (first1, first2) else {
        return;
    };

    let first1_flipped = writer_bond_is_flipped_from_atom1(
        molecule,
        atom1,
        first1,
        atom_visit_orders,
        traversal_ring_closure_bonds,
    );
    let second1_flipped = second1
        .map(|bond| {
            writer_bond_is_flipped_from_atom1(
                molecule,
                atom1,
                bond,
                atom_visit_orders,
                traversal_ring_closure_bonds,
            )
        })
        .unwrap_or(false);
    let first2_flipped = writer_bond_is_flipped_from_atom2(
        molecule,
        atom2,
        first2,
        atom_visit_orders,
        traversal_ring_closure_bonds,
    );
    let second2_flipped = second2
        .map(|bond| {
            writer_bond_is_flipped_from_atom2(
                molecule,
                atom2,
                bond,
                atom_visit_orders,
                traversal_ring_closure_bonds,
            )
        })
        .unwrap_or(false);

    if dir1_set && dir2_set {
        let atom1_consistent = account_same_side_dirs_for_writer(
            molecule,
            atom1,
            first1,
            first1_flipped,
            second1,
            second1_flipped,
            bond_dir_counts,
            atom_dir_counts,
        );
        let atom2_consistent = account_same_side_dirs_for_writer(
            molecule,
            atom2,
            first2,
            first2_flipped,
            second2,
            second2_flipped,
            bond_dir_counts,
            atom_dir_counts,
        );
        let _ = handle_dir_conflicts_across_double_bond_for_writer(
            molecule,
            dbl_bond,
            atom1,
            atom1_consistent,
            first1,
            first1_flipped,
            second1,
            second1_flipped,
            atom2,
            atom2_consistent,
            first2,
            first2_flipped,
            second2,
            second2_flipped,
            bond_dir_counts,
            atom_dir_counts,
        );
        return;
    }

    let mut set_from_atom1 = true;
    let mut atom1_controlling_bond = first1;
    let mut atom2_controlling_bond = first2;
    if !dir1_set && !dir2_set {
        set_bond_direction(molecule, first1, BondDirection::EndUpRight);
        bond_dir_counts[first1.index()] += 1;
        atom_dir_counts[atom1.index()] += 1;
    } else if !dir2_set {
        if bond_dir_counts[first1.index()] > 0 {
            bond_dir_counts[first1.index()] += 1;
            atom_dir_counts[atom1.index()] += 1;
            if let Some(second1) = second1
                && bond_dir_counts[second1.index()] > 0
            {
                bond_dir_counts[second1.index()] += 1;
                atom_dir_counts[atom1.index()] += 1;
            }
        } else if let Some(second1) = second1 {
            set_direction_from_neighboring_bond_for_writer(
                molecule,
                second1,
                second1_flipped,
                first1,
                first1_flipped,
            );
            bond_dir_counts[second1.index()] += 1;
            bond_dir_counts[first1.index()] += 1;
            atom_dir_counts[atom1.index()] += 2;
            atom1_controlling_bond = second1;
        }
    } else {
        set_from_atom1 = false;
        if bond_dir_counts[first2.index()] > 0 {
            bond_dir_counts[first2.index()] += 1;
            atom_dir_counts[atom2.index()] += 1;
            if let Some(second2) = second2
                && bond_dir_counts[second2.index()] > 0
            {
                bond_dir_counts[second2.index()] += 1;
                atom_dir_counts[atom2.index()] += 1;
            }
        } else if let Some(second2) = second2 {
            set_direction_from_neighboring_bond_for_writer(
                molecule,
                second2,
                second2_flipped,
                first2,
                first2_flipped,
            );
            bond_dir_counts[second2.index()] += 1;
            bond_dir_counts[first2.index()] += 1;
            atom_dir_counts[atom2.index()] += 2;
            atom2_controlling_bond = second2;
        }
    }

    if set_from_atom1 {
        let controlling_flipped = if atom1_controlling_bond == first1 {
            first1_flipped
        } else {
            second1_flipped
        };
        let dir = get_reference_direction_for_writer(
            molecule,
            dbl_bond,
            atom1,
            atom2,
            atom1_controlling_bond,
            controlling_flipped,
            first2,
            first2_flipped,
        );
        set_bond_direction(molecule, first2, dir);
        bond_dir_counts[first2.index()] += 1;
        atom_dir_counts[atom2.index()] += 1;
    } else {
        let controlling_flipped = if atom2_controlling_bond == first2 {
            first2_flipped
        } else {
            second2_flipped
        };
        let dir = get_reference_direction_for_writer(
            molecule,
            dbl_bond,
            atom2,
            atom1,
            atom2_controlling_bond,
            controlling_flipped,
            first1,
            first1_flipped,
        );
        set_bond_direction(molecule, first1, dir);
        bond_dir_counts[first1.index()] += 1;
        atom_dir_counts[atom1.index()] += 1;
    }

    if incident_bonds(molecule, atom1).len() == 3
        && let Some(second1) = second1
        && bond_dir_counts[second1.index()] == 0
    {
        set_direction_from_neighboring_bond_for_writer(
            molecule,
            first1,
            first1_flipped,
            second1,
            second1_flipped,
        );
        bond_dir_counts[second1.index()] += 1;
        atom_dir_counts[atom1.index()] += 1;
    }
    if incident_bonds(molecule, atom2).len() == 3
        && let Some(second2) = second2
        && bond_dir_counts[second2.index()] == 0
    {
        set_direction_from_neighboring_bond_for_writer(
            molecule,
            first2,
            first2_flipped,
            second2,
            second2_flipped,
        );
        bond_dir_counts[second2.index()] += 1;
        atom_dir_counts[atom2.index()] += 1;
    }
}

fn find_double_bond_neighbor_bonds_for_writer(
    molecule: &Molecule,
    dbl_bond: BondId,
    atom: AtomId,
    bond_visit_orders: &[usize],
    bond_dir_counts: &[i8],
) -> (Option<BondId>, Option<BondId>, bool) {
    let mut first = None;
    let mut second = None;
    let mut first_visit = usize::MAX;
    let mut dir_set = false;
    for bond in incident_bonds(molecule, atom) {
        if bond == dbl_bond
            || !can_set_double_bond_stereo_for_writer(molecule.bonds()[bond.index()].order())
        {
            continue;
        }
        if bond_dir_counts[bond.index()] > 0 {
            dir_set = true;
        }
        if first.is_none() || bond_visit_orders[bond.index()] < first_visit {
            if first.is_some() {
                second = first;
            }
            first = Some(bond);
            first_visit = bond_visit_orders[bond.index()];
        } else {
            second = Some(bond);
        }
    }
    (first, second, dir_set)
}

fn account_same_side_dirs_for_writer(
    molecule: &mut Molecule,
    atom: AtomId,
    first: BondId,
    first_flipped: bool,
    second: Option<BondId>,
    second_flipped: bool,
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
) -> bool {
    let mut consistent = true;
    if let Some(second) = second {
        if bond_dir_counts[first.index()] == 0 {
            set_direction_from_neighboring_bond_for_writer(
                molecule,
                second,
                second_flipped,
                first,
                first_flipped,
            );
        } else if bond_dir_counts[second.index()] == 0 {
            set_direction_from_neighboring_bond_for_writer(
                molecule,
                first,
                first_flipped,
                second,
                second_flipped,
            );
        } else {
            consistent = same_side_dirs_are_compatible_for_writer(
                molecule,
                first,
                first_flipped,
                second,
                second_flipped,
            );
        }
        bond_dir_counts[second.index()] += 1;
        atom_dir_counts[atom.index()] += 1;
    }
    bond_dir_counts[first.index()] += 1;
    atom_dir_counts[atom.index()] += 1;
    consistent
}

#[allow(clippy::too_many_arguments)]
fn handle_dir_conflicts_across_double_bond_for_writer(
    molecule: &mut Molecule,
    dbl_bond: BondId,
    atom1: AtomId,
    atom1_consistent: bool,
    first1: BondId,
    first1_flipped: bool,
    second1: Option<BondId>,
    second1_flipped: bool,
    atom2: AtomId,
    atom2_consistent: bool,
    first2: BondId,
    first2_flipped: bool,
    second2: Option<BondId>,
    second2_flipped: bool,
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
) -> bool {
    // BEGIN RDKIT CPP FUNCTION handleDirConflictsAcrossDoubleBond
    // RDKit✔️✔️: if (atom1DirsAreConsistent && atom2DirsAreConsistent) { ... }
    // RDKit✔️✔️: else if (!atom2DirsAreConsistent && atom1DirsAreConsistent) {
    // RDKit✔️✔️:   return fixConflictAcrossDoubleBond(...);
    // RDKit✔️✔️: } else if (!atom1DirsAreConsistent && atom2DirsAreConsistent) {
    // RDKit✔️✔️:   return fixConflictAcrossDoubleBond(...);
    // RDKit✔️✔️: } else { ... try all four direction combinations ... }
    // END RDKIT CPP FUNCTION handleDirConflictsAcrossDoubleBond
    if atom1_consistent && atom2_consistent {
        return get_reference_direction_for_writer(
            molecule,
            dbl_bond,
            atom1,
            atom2,
            first1,
            first1_flipped,
            first2,
            first2_flipped,
        ) == molecule.bonds()[first2.index()].direction();
    }
    if !atom2_consistent && atom1_consistent {
        if let Some(second2) = second2 {
            return fix_conflict_across_double_bond_for_writer(
                molecule,
                dbl_bond,
                atom2,
                first2,
                first2_flipped,
                second2,
                second2_flipped,
                atom1,
                first1,
                first1_flipped,
                bond_dir_counts,
                atom_dir_counts,
            );
        }
        return false;
    }
    if !atom1_consistent && atom2_consistent {
        if let Some(second1) = second1 {
            return fix_conflict_across_double_bond_for_writer(
                molecule,
                dbl_bond,
                atom1,
                first1,
                first1_flipped,
                second1,
                second1_flipped,
                atom2,
                first2,
                first2_flipped,
                bond_dir_counts,
                atom_dir_counts,
            );
        }
        return false;
    }
    let (Some(second1), Some(second2)) = (second1, second2) else {
        return false;
    };
    for (atom1_bond, atom1_flipped, atom1_other) in [
        (first1, first1_flipped, second1),
        (second1, second1_flipped, first1),
    ] {
        for (atom2_bond, atom2_flipped, atom2_other) in [
            (first2, first2_flipped, second2),
            (second2, second2_flipped, first2),
        ] {
            let expected = get_reference_direction_for_writer(
                molecule,
                dbl_bond,
                atom1,
                atom2,
                atom1_bond,
                atom1_flipped,
                atom2_bond,
                atom2_flipped,
            );
            if expected != molecule.bonds()[atom2_bond.index()].direction() {
                continue;
            }
            let Some(atom1_other_idx) =
                bond_other_atom(&molecule.bonds()[atom1_other.index()], atom1)
            else {
                continue;
            };
            let Some(atom2_other_idx) =
                bond_other_atom(&molecule.bonds()[atom2_other.index()], atom2)
            else {
                continue;
            };
            if atom1_other_idx == atom2_other_idx
                || atom_dir_counts[atom1_other_idx.index()] != 2
                || atom_dir_counts[atom2_other_idx.index()] != 2
            {
                continue;
            }
            bond_dir_counts[atom1_other.index()] = 0;
            atom_dir_counts[atom1.index()] -= 1;
            atom_dir_counts[atom1_other_idx.index()] -= 1;
            bond_dir_counts[atom2_other.index()] = 0;
            atom_dir_counts[atom2.index()] -= 1;
            atom_dir_counts[atom2_other_idx.index()] -= 1;
            return true;
        }
    }
    false
}

#[allow(clippy::too_many_arguments)]
fn fix_conflict_across_double_bond_for_writer(
    molecule: &Molecule,
    dbl_bond: BondId,
    atom: AtomId,
    first: BondId,
    first_flipped: bool,
    second: BondId,
    second_flipped: bool,
    ref_atom: AtomId,
    ref_bond: BondId,
    ref_flipped: bool,
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
) -> bool {
    for (bond, flipped, other_bond) in [
        (first, first_flipped, second),
        (second, second_flipped, first),
    ] {
        let Some(other_idx) = bond_other_atom(&molecule.bonds()[other_bond.index()], atom) else {
            continue;
        };
        if atom_dir_counts[other_idx.index()] != 2 {
            continue;
        }
        let expected = get_reference_direction_for_writer(
            molecule,
            dbl_bond,
            ref_atom,
            atom,
            ref_bond,
            ref_flipped,
            bond,
            flipped,
        );
        if expected == molecule.bonds()[bond.index()].direction() {
            bond_dir_counts[other_bond.index()] = 0;
            atom_dir_counts[atom.index()] -= 1;
            atom_dir_counts[other_idx.index()] -= 1;
            return true;
        }
    }
    false
}

fn remove_unwanted_bond_dir_specs_for_writer(
    molecule: &mut Molecule,
    stack: &[MolStackElem],
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
    bond_visit_orders: &[usize],
) {
    // BEGIN RDKIT CPP FUNCTION Canon::removeUnwantedBondDirSpecs
    // RDKit✔️✔️: for (auto &msI : molStack) {
    // RDKit✔️✔️:   if (msI.type != MOL_STACK_BOND) { continue; }
    // RDKit✔️✔️:   if (msI.obj.bond->getBondType() != Bond::DOUBLE ||
    // RDKit✔️✔️:       msI.obj.bond->getStereo() > Bond::STEREOANY) { continue; }
    // RDKit✔️✔️:   ... collect removable direction bonds on both ends ...
    // RDKit✔️✔️:   if (atomDirCounts[otherAtom->getIdx()] == 2) {
    // RDKit✔️✔️:     bondDirCounts[candidateBond->getIdx()] = 0;
    // RDKit✔️✔️:     candidateBond->setBondDir(Bond::NONE);
    // RDKit✔️✔️:     atomDirCounts[otherAtom->getIdx()] -= 1;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Canon::removeUnwantedBondDirSpecs
    for item in stack {
        let MolStackElem::Bond(dbl_bond, _) = *item else {
            continue;
        };
        let bond = &molecule.bonds()[dbl_bond.index()];
        if bond.order() != BondOrder::Double || is_writer_stereo_double_bond(bond) {
            continue;
        }
        let first = bond.begin();
        let second = bond.end();
        if incident_bonds(molecule, first).len() == 1 || incident_bonds(molecule, second).len() == 1
        {
            continue;
        }
        let mut candidates = Vec::new();
        for bond in incident_bonds(molecule, first) {
            if bond_dir_counts[bond.index()] > 0 {
                candidates.push(bond);
            }
        }
        if candidates.is_empty() {
            continue;
        }
        if atom_dir_counts[first.index()] > 0 {
            candidates.clear();
        }
        let mut second_end_candidates = 0usize;
        for bond in incident_bonds(molecule, second) {
            if bond_dir_counts[bond.index()] > 0 {
                candidates.push(bond);
                second_end_candidates += 1;
            }
        }
        if second_end_candidates == 0 || atom_dir_counts[second.index()] > 0 {
            continue;
        }
        candidates.sort_by_key(|bond| bond_visit_orders[bond.index()]);
        for candidate in candidates {
            let other = if bond_other_atom(&molecule.bonds()[candidate.index()], first).is_some() {
                bond_other_atom(&molecule.bonds()[candidate.index()], first)
            } else {
                bond_other_atom(&molecule.bonds()[candidate.index()], second)
            };
            let Some(other) = other else {
                continue;
            };
            if atom_dir_counts[other.index()] == 2 {
                bond_dir_counts[candidate.index()] = 0;
                molecule.topology_block_mut().bonds[candidate.index()]
                    .set_direction(BondDirection::None);
                atom_dir_counts[other.index()] -= 1;
                break;
            }
        }
    }
}

fn remove_redundant_bond_dir_specs_for_writer(
    molecule: &mut Molecule,
    stack: &[MolStackElem],
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
) {
    // BEGIN RDKIT CPP FUNCTION Canon::removeRedundantBondDirSpecs / clearBondDirs
    // RDKit✔️✔️: auto clearBondDirsFromAtom = [&mol, &bondDirCounts, &atomDirCounts](
    // RDKit✔️✔️:                                  Bond *tBond, const Atom *atom) {
    // RDKit✔️✔️:   if (atomDirCounts[atom->getIdx()] < 2) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (auto bond : mol.atomBonds(atom)) {
    // RDKit✔️✔️:     if (bond != tBond && bond->getBondType() == Bond::DOUBLE &&
    // RDKit✔️✔️:         bond->getStereo() > Bond::STEREOANY) {
    // RDKit✔️✔️:       clearBondDirs(mol, tBond, atom, bondDirCounts, atomDirCounts);
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: };
    // RDKit✔️✔️: if (canHaveDirection(*tBond) && bondDirCounts[tBond->getIdx()]) {
    // RDKit✔️✔️:   clearBondDirsFromAtom(tBond, canonBeginAtom);
    // RDKit✔️✔️:   clearBondDirsFromAtom(tBond, canonEndAtom);
    // RDKit✔️✔️: } else if (tBond->getBondDir() != Bond::NONE) {
    // RDKit✔️✔️:   tBond->setBondDir(Bond::NONE);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Canon::removeRedundantBondDirSpecs / clearBondDirs
    for item in stack {
        let MolStackElem::Bond(t_bond, atom_to_left) = *item else {
            continue;
        };
        let t_bond_ref = &molecule.bonds()[t_bond.index()];
        if can_have_direction_for_writer(t_bond_ref.order()) && bond_dir_counts[t_bond.index()] > 0
        {
            let canon_begin_atom = atom_to_left;
            let Some(canon_end_atom) = bond_other_atom(t_bond_ref, atom_to_left) else {
                continue;
            };

            clear_redundant_bond_dirs_from_atom_for_writer(
                molecule,
                t_bond,
                canon_begin_atom,
                bond_dir_counts,
                atom_dir_counts,
            );
            clear_redundant_bond_dirs_from_atom_for_writer(
                molecule,
                t_bond,
                canon_end_atom,
                bond_dir_counts,
                atom_dir_counts,
            );
        } else if molecule.bonds()[t_bond.index()].direction() != BondDirection::None {
            molecule.topology_block_mut().bonds[t_bond.index()].set_direction(BondDirection::None);
        }
    }
}

fn clear_redundant_bond_dirs_from_atom_for_writer(
    molecule: &mut Molecule,
    ref_bond: BondId,
    atom: AtomId,
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
) {
    if atom_dir_counts[atom.index()] < 2 {
        return;
    }
    for bond in incident_bonds(molecule, atom) {
        if bond != ref_bond
            && molecule.bonds()[bond.index()].order() == BondOrder::Double
            && matches!(
                molecule.bonds()[bond.index()].stereo(),
                BondStereo::Z
                    | BondStereo::E
                    | BondStereo::Cis
                    | BondStereo::Trans
                    | BondStereo::AtropCw
                    | BondStereo::AtropCcw
            )
        {
            clear_bond_dirs_from_atom_for_writer(
                molecule,
                ref_bond,
                atom,
                bond_dir_counts,
                atom_dir_counts,
            );
            return;
        }
    }
}

fn clear_bond_dirs_from_atom_for_writer(
    molecule: &mut Molecule,
    ref_bond: BondId,
    atom: AtomId,
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
) {
    // BEGIN RDKIT CPP FUNCTION clearBondDirs
    // RDKit✔️✔️: void clearBondDirs(ROMol &mol, Bond *refBond, const Atom *fromAtom,
    // RDKit✔️✔️:                    std::vector<int8_t> &bondDirCounts,
    // RDKit✔️✔️:                    std::vector<int8_t> &atomDirCounts) {
    // RDKit✔️✔️:   auto clearDirection = [&atomDirCounts, &bondDirCounts](Bond *bond,
    // RDKit✔️✔️:                                                          const Atom *fromAtom) {
    // RDKit✔️✔️:     --bondDirCounts[bond->getIdx()];
    // RDKit✔️✔️:     if (!bondDirCounts[bond->getIdx()]) {
    // RDKit✔️✔️:       bond->setBondDir(Bond::NONE);
    // RDKit✔️✔️:       --atomDirCounts[fromAtom->getIdx()];
    // RDKit✔️✔️:       if (auto otherAtom = bond->getOtherAtom(fromAtom);
    // RDKit✔️✔️:           atomDirCounts[otherAtom->getIdx()]) {
    // RDKit✔️✔️:         --atomDirCounts[otherAtom->getIdx()];
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   for (auto oBond : mol.atomBonds(fromAtom)) {
    // RDKit✔️✔️:     if (oBond != refBond && canHaveDirection(*oBond)) {
    // RDKit✔️✔️:       if ((bondDirCounts[oBond->getIdx()] >=
    // RDKit✔️✔️:            bondDirCounts[refBond->getIdx()]) &&
    // RDKit✔️✔️:           atomDirCounts[oBond->getBeginAtomIdx()] != 1 &&
    // RDKit✔️✔️:           atomDirCounts[oBond->getEndAtomIdx()] != 1) {
    // RDKit✔️✔️:         clearDirection(oBond, fromAtom);
    // RDKit✔️✔️:       } else if (atomDirCounts[refBond->getBeginAtomIdx()] != 1 &&
    // RDKit✔️✔️:                  atomDirCounts[refBond->getEndAtomIdx()] != 1) {
    // RDKit✔️✔️:         clearDirection(refBond, fromAtom);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION clearBondDirs
    if atom_dir_counts[atom.index()] < 2 {
        return;
    }
    for other_bond in incident_bonds(molecule, atom) {
        if other_bond == ref_bond
            || !can_have_direction_for_writer(molecule.bonds()[other_bond.index()].order())
        {
            continue;
        }
        let ref_count = bond_dir_counts[ref_bond.index()];
        let other_count = bond_dir_counts[other_bond.index()];
        let other_begin = molecule.bonds()[other_bond.index()].begin();
        let other_end = molecule.bonds()[other_bond.index()].end();
        if other_count >= ref_count
            && atom_dir_counts[other_begin.index()] != 1
            && atom_dir_counts[other_end.index()] != 1
        {
            clear_direction_counted_for_writer(
                molecule,
                other_bond,
                atom,
                bond_dir_counts,
                atom_dir_counts,
            );
        } else {
            let ref_begin = molecule.bonds()[ref_bond.index()].begin();
            let ref_end = molecule.bonds()[ref_bond.index()].end();
            if atom_dir_counts[ref_begin.index()] != 1 && atom_dir_counts[ref_end.index()] != 1 {
                clear_direction_counted_for_writer(
                    molecule,
                    ref_bond,
                    atom,
                    bond_dir_counts,
                    atom_dir_counts,
                );
            }
        }
        break;
    }
}

fn clear_direction_counted_for_writer(
    molecule: &mut Molecule,
    bond: BondId,
    from_atom: AtomId,
    bond_dir_counts: &mut [i8],
    atom_dir_counts: &mut [i8],
) {
    bond_dir_counts[bond.index()] -= 1;
    if bond_dir_counts[bond.index()] == 0 {
        molecule.topology_block_mut().bonds[bond.index()].set_direction(BondDirection::None);
        atom_dir_counts[from_atom.index()] -= 1;
        if let Some(other) = bond_other_atom(&molecule.bonds()[bond.index()], from_atom)
            && atom_dir_counts[other.index()] > 0
        {
            atom_dir_counts[other.index()] -= 1;
        }
    }
}

fn neighboring_stereo_double_bonds_for_writer(
    molecule: &Molecule,
    dbl_bond: BondId,
) -> Vec<BondId> {
    let mut out = Vec::new();
    let bond = &molecule.bonds()[dbl_bond.index()];
    for atom in [bond.begin(), bond.end()] {
        for nbr_bond in incident_bonds(molecule, atom) {
            if !can_have_direction_for_writer(molecule.bonds()[nbr_bond.index()].order()) {
                continue;
            }
            let Some(other_atom) = bond_other_atom(&molecule.bonds()[nbr_bond.index()], atom)
            else {
                continue;
            };
            for other_bond in incident_bonds(molecule, other_atom) {
                if other_bond != nbr_bond
                    && is_writer_stereo_double_bond(&molecule.bonds()[other_bond.index()])
                    && !out.contains(&other_bond)
                {
                    out.push(other_bond);
                }
            }
        }
    }
    out
}

#[allow(clippy::too_many_arguments)]
fn get_reference_direction_for_writer(
    molecule: &Molecule,
    dbl_bond: BondId,
    ref_atom: AtomId,
    target_atom: AtomId,
    ref_controlling_bond: BondId,
    ref_is_flipped: bool,
    target_bond: BondId,
    target_is_flipped: bool,
) -> BondDirection {
    // BEGIN RDKIT CPP FUNCTION getReferenceDirection
    // RDKit❗✔️: if (dblBond.getStereo() == Bond::STEREOE ||
    // RDKit❗✔️:     dblBond.getStereo() == Bond::STEREOTRANS) {
    // RDKit❗✔️:   dir = refControllingBond.getBondDir();
    // RDKit❗✔️: } else if (dblBond.getStereo() == Bond::STEREOZ ||
    // RDKit❗✔️:            dblBond.getStereo() == Bond::STEREOCIS) {
    // RDKit❗✔️:   dir = flipStereoBondDir(refControllingBond.getBondDir());
    // RDKit❗✔️: }
    // RDKit❗✔️: if (refAtom.getDegree() == 3 && ... stereoAtoms ... ) { dir = flipStereoBondDir(dir); }
    // RDKit❗✔️: if (targetAtom.getDegree() == 3 && ... stereoAtoms ... ) { dir = flipStereoBondDir(dir); }
    // RDKit❗✔️: if (refIsFlipped != targetIsFlipped) { dir = flipStereoBondDir(dir); }
    // END RDKIT CPP FUNCTION getReferenceDirection
    let dbl = &molecule.bonds()[dbl_bond.index()];
    let mut dir = match dbl.stereo() {
        BondStereo::E | BondStereo::Trans => {
            molecule.bonds()[ref_controlling_bond.index()].direction()
        }
        BondStereo::Z | BondStereo::Cis => flip_stereo_bond_dir_for_writer(
            molecule.bonds()[ref_controlling_bond.index()].direction(),
        ),
        _ => BondDirection::None,
    };
    if let Some(stereo_atoms) = dbl.stereo_atoms() {
        if incident_bonds(molecule, ref_atom).len() == 3
            && !stereo_atoms.contains(
                &bond_other_atom(&molecule.bonds()[ref_controlling_bond.index()], ref_atom)
                    .unwrap_or(ref_atom),
            )
        {
            dir = flip_stereo_bond_dir_for_writer(dir);
        }
        if incident_bonds(molecule, target_atom).len() == 3
            && !stereo_atoms.contains(
                &bond_other_atom(&molecule.bonds()[target_bond.index()], target_atom)
                    .unwrap_or(target_atom),
            )
        {
            dir = flip_stereo_bond_dir_for_writer(dir);
        }
    }
    if ref_is_flipped != target_is_flipped {
        dir = flip_stereo_bond_dir_for_writer(dir);
    }
    dir
}

fn set_direction_from_neighboring_bond_for_writer(
    molecule: &mut Molecule,
    source: BondId,
    source_flipped: bool,
    target: BondId,
    target_flipped: bool,
) {
    // BEGIN RDKIT CPP FUNCTION setDirectionFromNeighboringBond
    // RDKit❗✔️: auto dir = sourceBond.getBondDir();
    // RDKit❗✔️: if (isSourceBondFlipped == isTargetBondFlipped) {
    // RDKit❗✔️:   dir = flipStereoBondDir(dir);
    // RDKit❗✔️: }
    // RDKit❗✔️: targetBond.setBondDir(dir);
    // END RDKIT CPP FUNCTION setDirectionFromNeighboringBond
    let mut dir = molecule.bonds()[source.index()].direction();
    if source_flipped == target_flipped {
        dir = flip_stereo_bond_dir_for_writer(dir);
    }
    set_bond_direction(molecule, target, dir);
}

fn same_side_dirs_are_compatible_for_writer(
    molecule: &Molecule,
    first: BondId,
    first_flipped: bool,
    second: BondId,
    second_flipped: bool,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION sameSideDirsAreCompatible
    // RDKit❗✔️: auto dirsShouldMatch = isFirstBondFlipped != isSecondBondFlipped;
    // RDKit❗✔️: auto dirsMatch = firstBond.getBondDir() == secondBond.getBondDir();
    // RDKit❗✔️: return dirsMatch == dirsShouldMatch;
    // END RDKIT CPP FUNCTION sameSideDirsAreCompatible
    let dirs_should_match = first_flipped != second_flipped;
    let dirs_match =
        molecule.bonds()[first.index()].direction() == molecule.bonds()[second.index()].direction();
    dirs_match == dirs_should_match
}

fn writer_bond_is_flipped_from_atom1(
    molecule: &Molecule,
    atom: AtomId,
    bond: BondId,
    atom_visit_orders: &[usize],
    traversal_ring_closure_bonds: &[bool],
) -> bool {
    let anchor = bond_other_atom(&molecule.bonds()[bond.index()], atom).unwrap_or(atom);
    (atom_visit_orders[atom.index()] < atom_visit_orders[anchor.index()])
        != traversal_ring_closure_bonds[bond.index()]
}

fn writer_bond_is_flipped_from_atom2(
    molecule: &Molecule,
    atom: AtomId,
    bond: BondId,
    atom_visit_orders: &[usize],
    traversal_ring_closure_bonds: &[bool],
) -> bool {
    let anchor = bond_other_atom(&molecule.bonds()[bond.index()], atom).unwrap_or(atom);
    (atom_visit_orders[anchor.index()] < atom_visit_orders[atom.index()])
        != traversal_ring_closure_bonds[bond.index()]
}

fn flip_stereo_bond_dir_for_writer(direction: BondDirection) -> BondDirection {
    match direction {
        BondDirection::EndUpRight => BondDirection::EndDownRight,
        BondDirection::EndDownRight => BondDirection::EndUpRight,
        other => other,
    }
}

fn is_writer_stereo_double_bond(bond: &Bond) -> bool {
    bond.order() == BondOrder::Double
        && matches!(
            bond.stereo(),
            BondStereo::E | BondStereo::Z | BondStereo::Cis | BondStereo::Trans
        )
        && bond.stereo_atoms().is_some()
}

fn writer_double_bond_has_distinguishable_ends(
    molecule: &Molecule,
    dbl_bond: BondId,
    cip_ranks: Option<&[u32]>,
) -> bool {
    let Some(ranks) = cip_ranks else {
        return true;
    };
    let bond = &molecule.bonds()[dbl_bond.index()];
    for atom in [bond.begin(), bond.end()] {
        let neighbors = incident_bonds(molecule, atom)
            .into_iter()
            .filter(|candidate| *candidate != dbl_bond)
            .filter_map(|candidate| bond_other_atom(&molecule.bonds()[candidate.index()], atom))
            .collect::<Vec<_>>();
        if neighbors.len() == 2 {
            // BEGIN RDKIT CPP FUNCTION Chirality::findAtomNeighborDirHelper duplicate-neighbor check
            // RDKit✔️✔️: if (neighbors.size() == 2 &&
            // RDKit✔️✔️:     ranks[neighbors[0].first] == ranks[neighbors[1].first]) {
            // RDKit✔️✔️:   // the two substituents are identical, no stereochemistry here:
            // RDKit✔️✔️:   neighbors.clear();
            // RDKit✔️✔️: }
            // END RDKIT CPP FUNCTION Chirality::findAtomNeighborDirHelper duplicate-neighbor check
            if ranks.get(neighbors[0].index()) == ranks.get(neighbors[1].index()) {
                return false;
            }
        }
    }
    true
}

fn can_have_direction_for_writer(order: BondOrder) -> bool {
    matches!(order, BondOrder::Single | BondOrder::Aromatic)
}

fn can_set_double_bond_stereo_for_writer(order: BondOrder) -> bool {
    matches!(
        order,
        BondOrder::Single
            | BondOrder::Aromatic
            | BondOrder::Dative
            | BondOrder::DativeOne
            | BondOrder::DativeLeft
            | BondOrder::DativeRight
    )
}

fn set_bond_direction(molecule: &mut Molecule, bond: BondId, direction: BondDirection) {
    molecule.topology_block_mut().bonds[bond.index()].set_direction(direction);
}

fn write_ring_closure(
    smiles: &mut String,
    ring_idx: usize,
    context: &mut SmilesWriteContext,
) -> Result<(), SmilesWriteError> {
    if let Some(digit) = context.ring_closure_digits.get(&ring_idx).copied() {
        write_ring_index(smiles, digit);
        context.ring_closures_to_erase.push(ring_idx);
        return Ok(());
    }

    let digit = match (1..).find(|candidate| {
        !context
            .ring_closure_digits
            .values()
            .any(|digit| digit == candidate)
    }) {
        Some(d) => d,
        // [deferred] Ring closure digit exhaustion. The candidate search
        // (1..) finds the first unused digit; None only occurs if every
        // positive integer is in use, which requires >2^63 concurrent ring
        // closures — an impossible real-world molecule. This is a defensive
        // guard for extreme edge cases.
        None => {
            return invariant_stage_error(
                SmilesPlanStage::ShortTermBondWriter,
                "write_ring_closure() could not allocate a free ring index",
            );
        }
    };
    context.ring_closure_digits.insert(ring_idx, digit);
    write_ring_index(smiles, digit);
    Ok(())
}

fn write_ring_index(smiles: &mut String, digit: usize) {
    if digit < 10 {
        smiles.push(char::from(b'0' + digit as u8));
    } else if digit < 100 {
        smiles.push('%');
        smiles.push_str(&digit.to_string());
    } else {
        smiles.push_str("%(");
        smiles.push_str(&digit.to_string());
        smiles.push(')');
    }
}

fn bond_other_atom(bond: &Bond, atom: AtomId) -> Option<AtomId> {
    if bond.begin() == atom {
        Some(bond.end())
    } else if bond.end() == atom {
        Some(bond.begin())
    } else {
        None
    }
}

fn element_symbol(atomic_number: u8) -> Result<&'static str, SmilesWriteError> {
    match atomic_number {
        0 => Ok("*"),
        1 => Ok("H"),
        2 => Ok("He"),
        3 => Ok("Li"),
        4 => Ok("Be"),
        5 => Ok("B"),
        6 => Ok("C"),
        7 => Ok("N"),
        8 => Ok("O"),
        9 => Ok("F"),
        10 => Ok("Ne"),
        11 => Ok("Na"),
        12 => Ok("Mg"),
        13 => Ok("Al"),
        14 => Ok("Si"),
        15 => Ok("P"),
        16 => Ok("S"),
        17 => Ok("Cl"),
        18 => Ok("Ar"),
        19 => Ok("K"),
        20 => Ok("Ca"),
        21 => Ok("Sc"),
        22 => Ok("Ti"),
        23 => Ok("V"),
        24 => Ok("Cr"),
        25 => Ok("Mn"),
        26 => Ok("Fe"),
        27 => Ok("Co"),
        28 => Ok("Ni"),
        29 => Ok("Cu"),
        30 => Ok("Zn"),
        31 => Ok("Ga"),
        32 => Ok("Ge"),
        33 => Ok("As"),
        34 => Ok("Se"),
        35 => Ok("Br"),
        36 => Ok("Kr"),
        37 => Ok("Rb"),
        38 => Ok("Sr"),
        39 => Ok("Y"),
        40 => Ok("Zr"),
        41 => Ok("Nb"),
        42 => Ok("Mo"),
        43 => Ok("Tc"),
        44 => Ok("Ru"),
        45 => Ok("Rh"),
        46 => Ok("Pd"),
        47 => Ok("Ag"),
        48 => Ok("Cd"),
        49 => Ok("In"),
        50 => Ok("Sn"),
        51 => Ok("Sb"),
        52 => Ok("Te"),
        53 => Ok("I"),
        54 => Ok("Xe"),
        55 => Ok("Cs"),
        56 => Ok("Ba"),
        57 => Ok("La"),
        58 => Ok("Ce"),
        59 => Ok("Pr"),
        60 => Ok("Nd"),
        61 => Ok("Pm"),
        62 => Ok("Sm"),
        63 => Ok("Eu"),
        64 => Ok("Gd"),
        65 => Ok("Tb"),
        66 => Ok("Dy"),
        67 => Ok("Ho"),
        68 => Ok("Er"),
        69 => Ok("Tm"),
        70 => Ok("Yb"),
        71 => Ok("Lu"),
        72 => Ok("Hf"),
        73 => Ok("Ta"),
        74 => Ok("W"),
        75 => Ok("Re"),
        76 => Ok("Os"),
        77 => Ok("Ir"),
        78 => Ok("Pt"),
        79 => Ok("Au"),
        80 => Ok("Hg"),
        81 => Ok("Tl"),
        82 => Ok("Pb"),
        83 => Ok("Bi"),
        84 => Ok("Po"),
        85 => Ok("At"),
        86 => Ok("Rn"),
        87 => Ok("Fr"),
        88 => Ok("Ra"),
        89 => Ok("Ac"),
        90 => Ok("Th"),
        91 => Ok("Pa"),
        92 => Ok("U"),
        93 => Ok("Np"),
        94 => Ok("Pu"),
        95 => Ok("Am"),
        96 => Ok("Cm"),
        97 => Ok("Bk"),
        98 => Ok("Cf"),
        99 => Ok("Es"),
        100 => Ok("Fm"),
        101 => Ok("Md"),
        102 => Ok("No"),
        103 => Ok("Lr"),
        104 => Ok("Rf"),
        105 => Ok("Db"),
        106 => Ok("Sg"),
        107 => Ok("Bh"),
        108 => Ok("Hs"),
        109 => Ok("Mt"),
        110 => Ok("Ds"),
        111 => Ok("Rg"),
        112 => Ok("Cn"),
        113 => Ok("Nh"),
        114 => Ok("Fl"),
        115 => Ok("Mc"),
        116 => Ok("Lv"),
        117 => Ok("Ts"),
        118 => Ok("Og"),
        // RDKit✔️✔️: PeriodicTable returns "?" for unknown atomic numbers
        _ => Ok("?"),
    }
}

/// Assert that an atom can be written in SMILES without unsupported features.
/// Query atoms, radical-bearing atoms, and out-of-range elements produce
/// bracket-wrapped output through the standard path (get_atom_smiles handles
/// them like RDKit does by using element symbol + bracket notation).
/// This function is kept for the minimal-plain path which bypasses the
/// standard writer for performance.

fn assemble_fragment_smiles(
    fragment_results: Vec<FragmentWriteResult>,
    params: &SmilesWriteParams,
    context: &mut SmilesWriteContext,
) -> Result<String, SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles fragment assembly section
    // RDKit✔️✔️:   if (params.canonical) {
    // RDKit✔️✔️:     std::sort(tmp.begin(), tmp.end());
    // RDKit✔️✔️:   } else {  // Not canonical
    // RDKit✔️✔️:     for (auto &i : allAtomOrdering) {
    // RDKit✔️✔️:       flattenedAtomOrdering.insert(flattenedAtomOrdering.end(), i.begin(),
    // RDKit✔️✔️:                                    i.end());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (auto &i : allBondOrdering) {
    // RDKit✔️✔️:       flattenedBondOrdering.insert(flattenedBondOrdering.end(), i.begin(),
    // RDKit✔️✔️:                                    i.end());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (unsigned i = 0; i < vfragsmi.size(); ++i) {
    // RDKit✔️✔️:       result += vfragsmi[i];
    // RDKit✔️✔️:       if (i < vfragsmi.size() - 1) {
    // RDKit✔️✔️:         result += ".";
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   mol.setProp(common_properties::_smilesAtomOutputOrder, flattenedAtomOrdering,
    // RDKit✔️✔️:               true);
    // RDKit✔️✔️:   mol.setProp(common_properties::_smilesBondOutputOrder, flattenedBondOrdering,
    // RDKit✔️✔️:               true);
    // RDKit✔️✔️:   return result;
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles fragment assembly section
    if params.canonical {
        let mut sorted = fragment_results;
        sorted.sort_by(|left, right| left.smiles.cmp(&right.smiles));
        context.atom_output_order.clear();
        context.bond_output_order.clear();
        for fragment in &sorted {
            context
                .atom_output_order
                .extend(fragment.atom_ordering.iter().copied());
            context
                .bond_output_order
                .extend(fragment.bond_ordering.iter().copied());
        }
        return Ok(sorted
            .into_iter()
            .map(|fragment| fragment.smiles)
            .collect::<Vec<_>>()
            .join("."));
    }
    context.atom_output_order.clear();
    context.bond_output_order.clear();
    for fragment in &fragment_results {
        context
            .atom_output_order
            .extend(fragment.atom_ordering.iter().copied());
        context
            .bond_output_order
            .extend(fragment.bond_ordering.iter().copied());
    }
    Ok(fragment_results
        .into_iter()
        .map(|fragment| fragment.smiles)
        .collect::<Vec<_>>()
        .join("."))
}

fn validate_rooted_atom(
    molecule: &Molecule,
    params: &SmilesWriteParams,
) -> Result<(), SmilesWriteError> {
    if let Some(atom) = params.rooted_at_atom
        && atom >= molecule.num_atoms()
    {
        return Err(SmilesWriteError::RootedAtomOutOfRange { atom });
    }
    Ok(())
}

fn validate_fragment_api_inputs(
    molecule: &Molecule,
    params: &SmilesWriteParams,
    atoms_to_use: &[usize],
    bonds_to_use: Option<&[usize]>,
    atom_symbols: Option<&[String]>,
    bond_symbols: Option<&[String]>,
) -> Result<(), SmilesWriteError> {
    for atom in atoms_to_use {
        validate_atom_index(molecule, *atom)?;
    }
    if let Some(bonds_to_use) = bonds_to_use {
        for bond in bonds_to_use {
            validate_bond_index(molecule, *bond)?;
        }
    }
    if let Some(root) = params.rooted_at_atom
        && !atoms_to_use.contains(&root)
    {
        return Err(SmilesWriteError::RootedAtomNotInFragment { atom: root });
    }
    if bonds_to_use.is_none()
        && let Some(root) = params.rooted_at_atom
    {
        let fragment_count = crate::notation::fragment::get_fragment_atom_mapping(molecule)
            .into_iter()
            .max()
            .map_or(0, |max_fragment| max_fragment + 1);
        if fragment_count > 1 {
            return Err(SmilesWriteError::RootedAtomRequiresSingleFragment { atom: root });
        }
    }
    if let Some(atom_symbols) = atom_symbols
        && atom_symbols.len() < molecule.num_atoms()
    {
        return Err(SmilesWriteError::AtomSymbolsTooShort {
            len: atom_symbols.len(),
            expected: molecule.num_atoms(),
        });
    }
    if let Some(bond_symbols) = bond_symbols
        && bond_symbols.len() < molecule.num_bonds()
    {
        return Err(SmilesWriteError::BondSymbolsTooShort {
            len: bond_symbols.len(),
            expected: molecule.num_bonds(),
        });
    }
    Ok(())
}

fn validate_atom_index(molecule: &Molecule, atom: usize) -> Result<(), SmilesWriteError> {
    if atom >= molecule.num_atoms() {
        Err(SmilesWriteError::AtomOutOfRange { atom })
    } else {
        Ok(())
    }
}

fn validate_bond_index(molecule: &Molecule, bond: usize) -> Result<(), SmilesWriteError> {
    if bond >= molecule.num_bonds() {
        Err(SmilesWriteError::BondOutOfRange { bond })
    } else {
        Ok(())
    }
}

fn invariant_stage_error<T>(
    stage: SmilesPlanStage,
    message: &'static str,
) -> Result<T, SmilesWriteError> {
    Err(SmilesWriteError::InvariantViolation {
        stage: stage.as_str(),
        message,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ethane() -> Molecule {
        Molecule::from_smiles_with_sanitize("CC", false).unwrap()
    }

    #[test]
    fn mol_to_smiles_empty_molecule_returns_empty_string_like_rdkit_entrypoint() {
        let molecule = Molecule::from_smiles_with_sanitize("", false).unwrap();

        assert_eq!(
            mol_to_smiles(&molecule, &SmilesWriteParams::default()).unwrap(),
            ""
        );
    }

    #[test]
    fn molecule_to_smiles_writes_basic_default_smiles() {
        let smiles = ethane().to_smiles(true).unwrap();

        assert_eq!(smiles, "CC");
    }

    #[test]
    fn mol_to_smiles_rejects_invalid_root_before_writer_branches() {
        let mut params = SmilesWriteParams::default();
        params.rooted_at_atom = Some(2);

        let error = mol_to_smiles(&ethane(), &params).unwrap_err();

        assert_eq!(error, SmilesWriteError::RootedAtomOutOfRange { atom: 2 });
    }

    #[test]
    fn all_primary_smiles_writer_modes_fail_closed_until_ported() {
        let molecule = ethane();
        assert_eq!(
            mol_to_smiles(&molecule, &SmilesWriteParams::default()).unwrap(),
            "CC"
        );
        let mut params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            ..Default::default()
        };

        assert_eq!(mol_to_smiles(&molecule, &params).unwrap(), "CC");

        params.do_kekule = true;
        // Kekulization is now implemented; ethane (no aromatic bonds) succeeds.
        assert_eq!(mol_to_smiles(&molecule, &params).unwrap(), "CC");

        params.do_kekule = false;
        params.do_random = true;
        // Random SMILES is now implemented — ethane (simple molecule) succeeds
        let random = mol_to_smiles(&molecule, &params).unwrap();
        // For ethane, random SMILES should still produce a valid SMILES
        assert_eq!(
            random.len(),
            2,
            "random SMILES should be 2 chars: {random:?}"
        );

        // CX SMILES is now implemented — ethane with CX fields returns plain SMILES
        // (no CX data for ethane) or SMILES with empty CX extension.
        let result = mol_to_cx_smiles(
            &molecule,
            &SmilesWriteParams::default(),
            CxSmilesFields::ALL,
            RestoreBondDirOption::Clear,
        );
        assert!(result.is_ok(), "CX SMILES should succeed: {result:?}");
        // Ethane has no CX-specific data, so output is plain "CC"
        assert_eq!(result.unwrap(), "CC");

        // CX with no extra fields should also work
        let result = mol_to_cx_smiles(
            &molecule,
            &SmilesWriteParams::default(),
            CxSmilesFields::NONE,
            RestoreBondDirOption::None,
        );
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), "CC");
    }

    #[test]
    fn mol_to_random_smiles_vect_returns_requested_count_and_is_seed_reproducible() {
        let molecule = Molecule::from_smiles_with_sanitize("CC(C)CO", false).unwrap();

        let first =
            mol_to_random_smiles_vect(&molecule, 8, 0x1234, false, false, false, false).unwrap();
        let second =
            mol_to_random_smiles_vect(&molecule, 8, 0x1234, false, false, false, false).unwrap();
        let different_seed =
            mol_to_random_smiles_vect(&molecule, 8, 0x5678, false, false, false, false).unwrap();

        assert_eq!(first.len(), 8);
        assert_eq!(first, second);
        assert!(
            first.iter().any(|smiles| smiles != &first[0]),
            "random vector should contain traversal variation: {first:?}"
        );
        assert_ne!(first, different_seed);
    }

    #[test]
    fn random_writer_seeded_start_matches_rdkit_for_diatomic_direction() {
        let molecule = Molecule::from_smiles_with_sanitize("N#C", false).unwrap();

        let generated =
            mol_to_random_smiles_vect(&molecule, 1, 1337, true, false, false, false).unwrap();

        assert_eq!(generated, vec!["C#N".to_string()]);
    }

    #[test]
    fn mol_to_smiles_with_mode_returns_empty_string_for_empty_molecule_like_rdkit() {
        let molecule = Molecule::new();

        let result = mol_to_smiles_with_mode(
            &molecule,
            &SmilesWriteParams::default(),
            SmilesOutputMode::PlainSmiles,
        )
        .unwrap();

        assert_eq!(result, "");
    }

    #[test]
    fn prepare_plain_smiles_molecule_stashes_and_clears_atom_maps_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("[CH3:7][CH2:3][CH3:1]", false).unwrap();
        let mut working = molecule.clone();
        let params = SmilesWriteParams {
            canonical: false,
            do_isomeric_smiles: false,
            clean_stereo: false,
            ignore_atom_map_numbers: true,
            ..Default::default()
        };

        let saved = prepare_plain_smiles_molecule(&mut working, &params).unwrap();

        assert_eq!(saved, Some(vec![Some(7), Some(3), Some(1)]));
        assert_eq!(
            working
                .atoms()
                .iter()
                .map(|atom| atom.atom_map())
                .collect::<Vec<_>>(),
            vec![None, None, None]
        );
    }

    #[test]
    fn prepare_plain_smiles_molecule_clears_dummy_atom_maps_for_canonical_ranking() {
        let molecule = Molecule::from_smiles_with_sanitize("[*:7]C[CH3:1]", false).unwrap();
        let mut working = molecule.clone();
        let params = SmilesWriteParams {
            ignore_atom_map_numbers: true,
            ..Default::default()
        };

        let saved = prepare_plain_smiles_molecule(&mut working, &params).unwrap();

        assert_eq!(saved, Some(vec![Some(7), None, Some(1)]));
        assert_eq!(
            working
                .atoms()
                .iter()
                .map(|atom| atom.atom_map())
                .collect::<Vec<_>>(),
            vec![None, None, None]
        );
    }

    #[test]
    fn mol_to_smiles_ignore_atom_map_numbers_matches_rdkit_canonical_split_behavior() {
        let molecule = Molecule::from_smiles_with_sanitize("[CH3:7][CH2:3][CH3:1]", false).unwrap();
        let noncanonical = SmilesWriteParams {
            canonical: false,
            do_isomeric_smiles: false,
            clean_stereo: false,
            ignore_atom_map_numbers: true,
            ..Default::default()
        };
        let canonical = SmilesWriteParams {
            ignore_atom_map_numbers: true,
            ..Default::default()
        };

        assert_eq!(mol_to_smiles(&molecule, &noncanonical).unwrap(), "CCC");

        let canonical_smiles = mol_to_smiles(&molecule, &canonical).unwrap();
        assert!(canonical_smiles.contains(":7"));
        assert!(canonical_smiles.contains(":3"));
        assert!(canonical_smiles.contains(":1"));
    }

    #[test]
    fn mol_to_smiles_ignore_atom_map_numbers_matches_rdkit_dummy_map_ordering() {
        let molecule = Molecule::from_smiles_with_sanitize("[*:1]C", false).unwrap();
        let noncanonical = SmilesWriteParams {
            do_isomeric_smiles: false,
            clean_stereo: false,
            canonical: false,
            ignore_atom_map_numbers: true,
            ..Default::default()
        };
        let canonical = SmilesWriteParams {
            do_isomeric_smiles: false,
            clean_stereo: false,
            canonical: true,
            ignore_atom_map_numbers: true,
            ..Default::default()
        };

        assert_eq!(mol_to_smiles(&molecule, &noncanonical).unwrap(), "*C");
        assert_eq!(mol_to_smiles(&molecule, &canonical).unwrap(), "[*:1]C");
    }

    #[test]
    fn mol_to_smiles_do_kekule_handles_exocyclic_aryl_substituent_like_rdkit_row_88() {
        let molecule =
            Molecule::from_smiles("O=C1N(/N=C(/C)C1=NN/C2=C/C(OC)=CC=C2)C=3C=CC=CC=3").unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            do_kekule: true,
            canonical: false,
            clean_stereo: false,
            ..Default::default()
        };

        assert_eq!(
            mol_to_smiles(&molecule, &params).unwrap(),
            "O=C1N(C2=CC=CC=C2)N=C(C)C1=NNC1=CC(OC)=CC=C1"
        );
    }

    #[test]
    fn writer_nonisomeric_explicit_bonds_clears_imine_direction_marks_like_rdkit_row_88() {
        let molecule =
            Molecule::from_smiles("O=C1N(/N=C(/C)C1=NN/C2=C/C(OC)=CC=C2)C=3C=CC=CC=3").unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            all_bonds_explicit: true,
            ..Default::default()
        };

        assert_eq!(
            mol_to_smiles(&molecule, &params).unwrap(),
            "O=C1-N(-c2:c:c:c:c:c:2)-N=C(-C)-C-1=N-N-c1:c:c(-O-C):c:c:c:1"
        );
    }

    #[test]
    fn choose_fragment_start_atom_prefers_terminal_dummy_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("*C", false).unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            clean_stereo: false,
            ..Default::default()
        };
        let plan = collect_fragment_write_plans(&molecule, &params)
            .unwrap()
            .into_iter()
            .next()
            .unwrap();
        let ranks = rank_fragment_atoms_for_smiles(
            &molecule,
            &plan,
            &params,
            SmilesOutputMode::PlainSmiles,
        )
        .unwrap();

        assert_eq!(ranks, vec![0, 1]);
        assert_eq!(
            choose_fragment_start_atom(&plan, &ranks, &params).unwrap(),
            AtomId::new(0)
        );
    }

    #[test]
    fn collect_fragment_write_plans_tracks_fragment_root_and_bond_membership() {
        let molecule = Molecule::from_smiles_with_sanitize("CC.CCO", false).unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            rooted_at_atom: Some(3),
            ..Default::default()
        };

        let plans = collect_fragment_write_plans(&molecule, &params).unwrap();

        assert_eq!(plans.len(), 2);
        assert_eq!(plans[0].atoms, vec![AtomId::new(0), AtomId::new(1)]);
        assert_eq!(plans[0].bonds, vec![BondId::new(0)]);
        assert_eq!(plans[0].rooted_at_atom, None);
        assert_eq!(
            plans[1].atoms,
            vec![AtomId::new(2), AtomId::new(3), AtomId::new(4)]
        );
        assert_eq!(plans[1].bonds, vec![BondId::new(1), BondId::new(2)]);
        assert_eq!(plans[1].rooted_at_atom, Some(AtomId::new(3)));
    }

    #[test]
    fn choose_fragment_start_atom_consumes_random_seed_stream_sequentially() {
        let plan = FragmentWritePlan {
            atoms: vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)],
            bonds: Vec::new(),
            rooted_at_atom: None,
        };
        let params = SmilesWriteParams {
            canonical: false,
            do_random: true,
            ..Default::default()
        };

        let chosen = with_random_smiles_seed(0x1234, || {
            Ok(vec![
                choose_fragment_start_atom(&plan, &[0, 1, 2], &params)?.index(),
                choose_fragment_start_atom(&plan, &[0, 1, 2], &params)?.index(),
                choose_fragment_start_atom(&plan, &[0, 1, 2], &params)?.index(),
            ])
        })
        .unwrap();

        let mut seed = 0x1234;
        let mut expected = Vec::new();
        for _ in 0..3 {
            expected.push(plan.atoms[(seed as usize) % plan.atoms.len()].index());
            seed = splitmix64(seed);
        }

        assert_eq!(chosen, expected);
    }

    #[test]
    fn rooted_writer_only_reorders_the_containing_fragment() {
        let molecule = Molecule::from_smiles_with_sanitize("CC.CCO", false).unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            rooted_at_atom: Some(3),
            ..Default::default()
        };

        assert_eq!(mol_to_smiles(&molecule, &params).unwrap(), "CC.C(C)O");
    }

    #[test]
    fn write_fragment_smiles_uses_fragment_plan_root_for_subfragment_output() {
        let mut molecule = Molecule::from_smiles_with_sanitize("CC.CCO", false).unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            rooted_at_atom: Some(3),
            ..Default::default()
        };
        let plans = collect_fragment_write_plans(&molecule, &params).unwrap();
        let mut context = SmilesWriteContext::default();

        let fragment = write_fragment_smiles(
            &mut molecule,
            &plans[1],
            &params,
            SmilesOutputMode::PlainSmiles,
            SmilesWriteOverrides::default(),
            &mut context,
        )
        .unwrap();

        assert_eq!(fragment.smiles, "C(C)O");
        assert_eq!(
            fragment.atom_ordering,
            vec![AtomId::new(3), AtomId::new(2), AtomId::new(4)]
        );
    }

    #[test]
    fn canonical_writer_uses_traversal_order_for_tetrahedral_chirality() {
        let molecule = Molecule::from_smiles_with_sanitize("C[C@H](F)CCCl", true).unwrap();
        let mut params = SmilesWriteParams::default();

        params.rooted_at_atom = Some(1);
        assert_eq!(
            mol_to_smiles(&molecule, &params).unwrap(),
            "[C@@H](C)(F)CCCl"
        );

        params.rooted_at_atom = Some(2);
        assert_eq!(mol_to_smiles(&molecule, &params).unwrap(), "F[C@@H](C)CCCl");
    }

    #[test]
    fn rooted_canonical_writer_matches_rdkit_branch_order_for_multichiral_case() {
        let molecule =
            Molecule::from_smiles_with_sanitize("O[C@](C)(Cl)[C@@](O)(Cl)C", true).unwrap();
        let params = SmilesWriteParams {
            rooted_at_atom: Some(0),
            ..Default::default()
        };

        assert_eq!(
            mol_to_smiles(&molecule, &params).unwrap(),
            "O[C@](C)(Cl)[C@@](C)(O)Cl"
        );
    }

    #[test]
    fn canonical_dfs_traversal_treats_single_h_query_as_fourth_valence_for_writer() {
        let mut builder = crate::MoleculeBuilder::new();
        let center = builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_query(
            crate::QueryNode::and(vec![
                crate::QueryNode::predicate(crate::AtomQueryPredicate::AtomicNumber(6)),
                crate::QueryNode::predicate(crate::AtomQueryPredicate::ImplicitHydrogenCount(1)),
            ]),
        ));
        let fluorine = builder.add_atom(crate::AtomSpec::new(crate::Element::F));
        let chlorine = builder.add_atom(crate::AtomSpec::new(crate::Element::CL));
        let bromine = builder.add_atom(crate::AtomSpec::new(crate::Element::BR));
        builder
            .add_bond(crate::BondSpec::new(center, fluorine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(center, chlorine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(center, bromine, BondOrder::Single))
            .unwrap();
        let mut molecule = builder.build().unwrap();
        molecule.topology_block_mut().atoms[center.index()]
            .set_chiral_tag(ChiralTag::TetrahedralCw);

        assert!(atom_has_fourth_valence_for_writer(&molecule, center));
        assert!(!chiral_atom_needs_tag_inversion_for_writer(
            &molecule, center, false, 1
        ));
    }

    #[test]
    fn canonical_dfs_traversal_treats_implicit_h_valence_as_fourth_for_writer() {
        let molecule = Molecule::from_smiles_with_sanitize("F[C@H](Cl)Br", true).unwrap();
        let center = AtomId::new(1);

        assert!(atom_has_fourth_valence_for_writer(&molecule, center));
        assert!(!chiral_atom_needs_tag_inversion_for_writer(
            &molecule, center, false, 1
        ));
    }

    #[test]
    fn canonical_dfs_traversal_inserts_implicit_nontetrahedral_neighbors_like_rdkit() {
        let mut first_atom_bonds = vec![Some(BondId::new(3)), Some(BondId::new(4))];
        insert_implicit_nontetrahedral_neighbors_for_writer(
            &mut first_atom_bonds,
            ChiralTag::SquarePlanar,
            true,
        );
        assert_eq!(
            first_atom_bonds,
            vec![None, None, Some(BondId::new(3)), Some(BondId::new(4))]
        );

        let mut later_atom_bonds = vec![Some(BondId::new(3)), Some(BondId::new(4))];
        insert_implicit_nontetrahedral_neighbors_for_writer(
            &mut later_atom_bonds,
            ChiralTag::SquarePlanar,
            false,
        );
        assert_eq!(
            later_atom_bonds,
            vec![Some(BondId::new(3)), None, None, Some(BondId::new(4))]
        );
    }

    #[test]
    fn canonical_dfs_traversal_marks_broken_chirality_for_partial_plan() {
        let molecule = Molecule::from_smiles_with_sanitize("F[C@](Cl)(Br)I", true).unwrap();
        let plan = FragmentWritePlan {
            atoms: vec![
                AtomId::new(0),
                AtomId::new(1),
                AtomId::new(2),
                AtomId::new(3),
            ],
            bonds: vec![BondId::new(0), BondId::new(1), BondId::new(2)],
            rooted_at_atom: Some(AtomId::new(1)),
        };
        let adjustments = compute_writer_chiral_adjustments(
            &molecule,
            &plan,
            AtomId::new(1),
            &vec![Vec::new(); molecule.num_atoms()],
            &vec![Vec::new(); molecule.num_atoms()],
            &[],
            false,
        )
        .unwrap();

        assert!(adjustments.broken_chiral_atoms.contains(&AtomId::new(1)));
        assert!(adjustments.chiral_inversions.is_empty());
        assert!(adjustments.chiral_permutations.is_empty());
    }

    #[test]
    fn fragment_writer_suppresses_chirality_when_fragment_breaks_chiral_incident_bonds() {
        let molecule = Molecule::from_smiles_with_sanitize("C[C@H](F)CCCl", true).unwrap();
        let params = SmilesWriteParams {
            canonical: false,
            ..Default::default()
        };

        assert_eq!(
            mol_fragment_to_smiles(&molecule, &params, &[0, 1, 2], None, None, None).unwrap(),
            "CCF"
        );
        assert_eq!(
            mol_fragment_to_smiles(&molecule, &params, &[0, 1, 2, 3], None, None, None).unwrap(),
            "C[C@H](F)C"
        );
    }

    #[test]
    fn writer_uses_parser_persisted_ring_stereo_props_without_local_reconstruction() {
        let molecule =
            Molecule::from_smiles_with_sanitize("C1[C@H](F)CC[C@H](Cl)C1", true).unwrap();
        let params = SmilesWriteParams {
            canonical: false,
            ..Default::default()
        };

        assert_eq!(
            mol_to_smiles(&molecule, &params).unwrap(),
            "C1[C@H](F)CC[C@H](Cl)C1"
        );
    }

    #[test]
    fn writer_rejects_ring_stereo_candidate_missing_ring_neighbors_prop() {
        let mut molecule =
            Molecule::from_smiles_with_sanitize("C1[C@H](F)CC[C@H](Cl)C1", true).unwrap();
        molecule.topology_block_mut().atoms[1].clear_prop("_ringStereoAtoms");

        let error = mol_to_smiles(&molecule, &SmilesWriteParams::default()).unwrap_err();

        assert_eq!(
            error,
            SmilesWriteError::InvalidRingStereoState {
                atom: 1,
                requirement: "`_ringStereochemCand` requires `_ringStereoAtoms`",
            }
        );
    }

    #[test]
    fn writer_rejects_malformed_ring_stereo_neighbors_prop() {
        let mut molecule =
            Molecule::from_smiles_with_sanitize("C1[C@H](F)CC[C@H](Cl)C1", true).unwrap();
        molecule.topology_block_mut().atoms[1].set_prop("_ringStereoAtoms", "0");

        let error = mol_to_smiles(&molecule, &SmilesWriteParams::default()).unwrap_err();

        assert_eq!(
            error,
            SmilesWriteError::InvalidRingStereoState {
                atom: 1,
                requirement: "`_ringStereoAtoms` must be a valid encoded ring-neighbor list",
            }
        );
    }

    #[test]
    fn mol_fragment_to_smiles_rejects_rooted_atom_for_multifragment_molecule_without_bond_scope() {
        let molecule = Molecule::from_smiles_with_sanitize("CC.O", false).unwrap();
        let params = SmilesWriteParams {
            rooted_at_atom: Some(0),
            ..Default::default()
        };

        let error =
            mol_fragment_to_smiles(&molecule, &params, &[0, 1], None, None, None).unwrap_err();

        assert_eq!(
            error,
            SmilesWriteError::RootedAtomRequiresSingleFragment { atom: 0 }
        );
    }

    #[test]
    fn isomeric_writer_outputs_non_tetrahedral_chiral_classes_like_rdkit() {
        let params = SmilesWriteParams {
            canonical: false,
            clean_stereo: false,
            ..Default::default()
        };

        let square_planar =
            Molecule::from_smiles_with_sanitize("[Pt@SP2](Cl)(Br)(I)F", false).unwrap();
        let trigonal_bipyramidal =
            Molecule::from_smiles_with_sanitize("[P@TB20](F)(Cl)(Br)(I)C", false).unwrap();
        let octahedral =
            Molecule::from_smiles_with_sanitize("[Co@OH30](F)(Cl)(Br)(I)(N)C", false).unwrap();

        assert_eq!(
            mol_to_smiles(&square_planar, &params).unwrap(),
            "[Pt@SP2]([Cl])([Br])([I])[F]"
        );
        assert_eq!(
            mol_to_smiles(&trigonal_bipyramidal, &params).unwrap(),
            "[P@TB20](F)(Cl)(Br)(I)C"
        );
        assert_eq!(
            mol_to_smiles(&octahedral, &params).unwrap(),
            "[Co@OH30]([F])([Cl])([Br])([I])([NH2])[CH3]"
        );
    }

    #[test]
    fn isomeric_writer_outputs_default_non_tetrahedral_chiral_classes_like_rdkit() {
        let params = SmilesWriteParams {
            canonical: false,
            clean_stereo: false,
            ..Default::default()
        };

        let square_planar =
            Molecule::from_smiles_with_sanitize("[Pt@SP](Cl)(Br)(I)F", false).unwrap();
        let trigonal_bipyramidal =
            Molecule::from_smiles_with_sanitize("[P@TB](F)(Cl)(Br)(I)C", false).unwrap();
        let octahedral =
            Molecule::from_smiles_with_sanitize("[Co@OH](F)(Cl)(Br)(I)(N)C", false).unwrap();

        assert_eq!(
            mol_to_smiles(&square_planar, &params).unwrap(),
            "[Pt@SP]([Cl])([Br])([I])[F]"
        );
        assert_eq!(
            mol_to_smiles(&trigonal_bipyramidal, &params).unwrap(),
            "[P@TB](F)(Cl)(Br)(I)C"
        );
        assert_eq!(
            mol_to_smiles(&octahedral, &params).unwrap(),
            "[Co@OH]([F])([Cl])([Br])([I])([NH2])[CH3]"
        );
    }

    #[test]
    fn rooted_writer_recomputes_non_tetrahedral_permutations_like_rdkit() {
        let mut params = SmilesWriteParams {
            canonical: false,
            clean_stereo: false,
            ..Default::default()
        };
        let square_planar =
            Molecule::from_smiles_with_sanitize("[Pt@SP2](Cl)(Br)(I)F", false).unwrap();
        let trigonal_bipyramidal =
            Molecule::from_smiles_with_sanitize("[P@TB20](F)(Cl)(Br)(I)C", false).unwrap();

        params.rooted_at_atom = Some(3);
        assert_eq!(
            mol_to_smiles(&square_planar, &params).unwrap(),
            "[I][Pt@SP3]([Cl])([Br])[F]"
        );

        params.rooted_at_atom = Some(4);
        assert_eq!(
            mol_to_smiles(&square_planar, &params).unwrap(),
            "[F][Pt@SP3]([Cl])([Br])[I]"
        );

        params.rooted_at_atom = Some(2);
        assert_eq!(
            mol_to_smiles(&trigonal_bipyramidal, &params).unwrap(),
            "Cl[P@TB15](F)(Br)(I)C"
        );

        params.rooted_at_atom = Some(5);
        assert_eq!(
            mol_to_smiles(&trigonal_bipyramidal, &params).unwrap(),
            "C[P@TB3](F)(Cl)(Br)I"
        );
    }

    #[test]
    fn noncanonical_nonisomeric_plain_smiles_writes_simple_linear_fragments() {
        let molecule = Molecule::from_smiles_with_sanitize("CCO.N=O", false).unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            ..Default::default()
        };

        assert_eq!(molecule.to_smiles_with_params(&params).unwrap(), "CCO.N=O");
    }

    #[test]
    fn canonical_fragment_assembly_reorders_output_scope_with_sorted_smiles() {
        let molecule = Molecule::from_smiles_with_sanitize("O.C", false).unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: true,
            clean_stereo: false,
            ..Default::default()
        };
        let mut context = SmilesWriteContext::default();

        let smiles = mol_fragment_to_smiles_with_context(
            &molecule,
            &params,
            &[0, 1],
            None,
            None,
            None,
            &mut context,
        )
        .unwrap();

        assert_eq!(smiles, "C.O");
        assert_eq!(
            context.atom_output_order,
            vec![AtomId::new(1), AtomId::new(0)]
        );
    }

    #[test]
    fn noncanonical_nonisomeric_plain_smiles_writes_explicit_single_bonds() {
        let molecule = Molecule::from_smiles_with_sanitize("CC", false).unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            all_bonds_explicit: true,
            ..Default::default()
        };

        assert_eq!(mol_to_smiles(&molecule, &params).unwrap(), "C-C");
    }

    #[test]
    fn noncanonical_nonisomeric_plain_smiles_writes_bracket_atoms_from_explicit_state() {
        let molecule =
            Molecule::from_smiles_with_sanitize("[NH4+].[O-].[SiH2].[*:7]", false).unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            ..Default::default()
        };

        assert_eq!(
            mol_to_smiles(&molecule, &params).unwrap(),
            "[NH4+].[O-].[SiH2].[*:7]"
        );
    }

    #[test]
    fn noncanonical_nonisomeric_plain_smiles_writes_branches_and_rings() {
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            ..Default::default()
        };
        let ring = Molecule::from_smiles_with_sanitize("C1CC1", false).unwrap();
        let branch = Molecule::from_smiles_with_sanitize("CC(C)O", false).unwrap();
        let nested = Molecule::from_smiles_with_sanitize("C1CCC(CC1)O", false).unwrap();

        assert_eq!(mol_to_smiles(&ring, &params).unwrap(), "C1CC1");
        assert_eq!(mol_to_smiles(&branch, &params).unwrap(), "CC(C)O");
        assert_eq!(mol_to_smiles(&nested, &params).unwrap(), "C1CCC(O)CC1");
    }

    #[test]
    fn noncanonical_nonisomeric_plain_smiles_writes_dative_bonds_and_ring_bond_symbols() {
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            ..Default::default()
        };
        let dative = Molecule::from_smiles_with_sanitize("N->O", false).unwrap();
        let rooted_dative = SmilesWriteParams {
            rooted_at_atom: Some(1),
            ..params.clone()
        };
        let opening_double = Molecule::from_smiles_with_sanitize("C=1CC1", false).unwrap();
        let closing_double = Molecule::from_smiles_with_sanitize("C1CC=1", false).unwrap();

        assert_eq!(mol_to_smiles(&dative, &params).unwrap(), "N->O");
        assert_eq!(mol_to_smiles(&dative, &rooted_dative).unwrap(), "O<-N");
        assert_eq!(mol_to_smiles(&opening_double, &params).unwrap(), "C1=CC1");
        assert_eq!(mol_to_smiles(&closing_double, &params).unwrap(), "C1=CC1");
    }

    #[test]
    fn plain_smiles_writer_strips_dative_bonds_on_working_copy_when_requested() {
        let molecule = Molecule::from_smiles_with_sanitize("N->O", false).unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            include_dative_bonds: false,
            ..Default::default()
        };

        let mut working = molecule.clone();
        let _ = prepare_plain_smiles_molecule(&mut working, &params).unwrap();
        assert_eq!(working.bonds()[0].order(), BondOrder::Single);
        assert_eq!(total_num_hydrogens_for_writer(&working, AtomId::new(0)), 3);
        assert_eq!(total_valence_for_writer(&working, AtomId::new(0)), Some(3));

        assert_eq!(mol_to_smiles(&molecule, &params).unwrap(), "NO");
        assert_eq!(molecule.bonds()[0].order(), BondOrder::Dative);
    }

    #[test]
    fn plain_smiles_writer_clears_molblock_only_bond_state_on_working_copy() {
        let mut builder = crate::MoleculeBuilder::new();
        let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let c2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(
                crate::BondSpec::new(c1, c2, BondOrder::Double)
                    .with_direction(BondDirection::EitherDouble)
                    .with_stereo(BondStereo::Any),
            )
            .unwrap();
        let molecule = builder.build().unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            ..Default::default()
        };

        assert_eq!(mol_to_smiles(&molecule, &params).unwrap(), "C=C");
        assert_eq!(molecule.bonds()[0].direction(), BondDirection::EitherDouble);
        assert_eq!(molecule.bonds()[0].stereo(), BondStereo::Any);
    }

    #[test]
    fn get_bond_smiles_ignores_direction_on_double_bond_like_rdkit() {
        let mut builder = crate::MoleculeBuilder::new();
        let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let c2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(
                crate::BondSpec::new(c1, c2, BondOrder::Double)
                    .with_direction(BondDirection::BeginWedge),
            )
            .unwrap();
        let molecule = builder.build().unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            ..Default::default()
        };

        assert_eq!(
            get_molecule_bond_smiles(&molecule, 0, Some(0), &params).unwrap(),
            "="
        );
    }

    #[test]
    fn cx_smiles_restore_bond_dirs_clear_uses_working_copy_only() {
        let mut builder = crate::MoleculeBuilder::new();
        let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let c2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(
                crate::BondSpec::new(c1, c2, BondOrder::Single)
                    .with_direction(BondDirection::Unknown),
            )
            .unwrap();
        let molecule = builder.build().unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            ..Default::default()
        };

        assert_eq!(
            mol_to_cx_smiles(
                &molecule,
                &params,
                CxSmilesFields::ALL,
                RestoreBondDirOption::Clear,
            )
            .unwrap(),
            "CC"
        );
        assert_eq!(molecule.bonds()[0].direction(), BondDirection::Unknown);
    }

    #[test]
    fn cx_preparation_does_not_apply_plain_only_bond_direction_cleanup() {
        let mut builder = crate::MoleculeBuilder::new();
        let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let c2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(
                crate::BondSpec::new(c1, c2, BondOrder::Double)
                    .with_direction(BondDirection::EitherDouble)
                    .with_stereo(BondStereo::Any),
            )
            .unwrap();
        let molecule = builder.build().unwrap();
        let mut params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            ..Default::default()
        };

        let mut keep = molecule.clone();
        prepare_cx_smiles_molecule(
            &mut keep,
            &mut params,
            CxSmilesFields::ALL,
            RestoreBondDirOption::None,
            true,
        )
        .unwrap();
        assert_eq!(keep.bonds()[0].direction(), BondDirection::EitherDouble);
        assert_eq!(keep.bonds()[0].stereo(), BondStereo::Any);

        let mut clear = molecule.clone();
        prepare_cx_smiles_molecule(
            &mut clear,
            &mut params,
            CxSmilesFields::ALL,
            RestoreBondDirOption::Clear,
            true,
        )
        .unwrap();
        assert_eq!(clear.bonds()[0].direction(), BondDirection::None);
        assert_eq!(clear.bonds()[0].stereo(), BondStereo::Any);
    }

    #[test]
    fn noncanonical_nonisomeric_plain_smiles_writes_aromatic() {
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            ..Default::default()
        };
        let benzene = Molecule::from_smiles_with_sanitize("c1ccccc1", false).unwrap();
        assert_eq!(mol_to_smiles(&benzene, &params).unwrap(), "c1ccccc1");

        // Also test a bracketed aromatic atom (Cl with aromatic flag)
        let params_with_h = SmilesWriteParams {
            all_hydrogens_explicit: true,
            ..params
        };
        let mol = Molecule::from_smiles_with_sanitize("c1ccc(C)cc1", false).unwrap();
        // should produce lowercase c's for aromatic carbons and uppercase C for sp3 carbon
        let smi = mol_to_smiles(&mol, &params_with_h).unwrap();
        assert!(
            smi.contains('c'),
            "expected lowercase c for aromatic atoms, got: {smi}"
        );
    }

    #[test]
    fn get_bond_smiles_and_atom_needs_bracket_mark_aromatic_nh_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("c1cc[nH]c1", false).unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            ..Default::default()
        };

        assert!(atom_needs_bracket(&molecule, AtomId::new(3), "", &params).unwrap());
        assert_eq!(mol_to_smiles(&molecule, &params).unwrap(), "c1cc[nH]c1");
    }

    #[test]
    fn noncanonical_nonisomeric_plain_smiles_writes_zero_order_and_radical_bonds_with_rdkit_compatibility()
     {
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            ..Default::default()
        };
        // Zero-order bonds: RDKit outputs "~" for unknown/zero bond types.
        // The molecule CC~CC parses as 4 carbons with a zero-order bond (idx 1).
        let zero = Molecule::from_smiles_with_sanitize("CC~CC", false).unwrap();
        let output = mol_to_smiles(&zero, &params).unwrap();
        // Should contain the "~" bond symbol in the output
        assert!(
            output.contains('~'),
            "zero-order bond should map to ~, got: {output:?}"
        );

        // Radical-bearing molecules: the radical state is preserved through
        // the typed state and written as bracket notation when needed.
        // Build a molecule with an explicit radical atom.
        let mut builder = crate::MoleculeBuilder::new();
        let c1 =
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_radical_electrons(1));
        let c2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(c1, c2, crate::BondOrder::Single))
            .unwrap();
        let radical = builder.build().unwrap();
        let output = mol_to_smiles(&radical, &params).unwrap();
        // Radical-bearing atoms get bracket notation: [CH3] or similar
        assert!(
            output.contains('['),
            "radical atom needs bracket: {output:?}"
        );
    }

    #[test]
    fn noncanonical_nonisomeric_plain_smiles_honors_rooted_atom_for_traversal_start() {
        let molecule = Molecule::from_smiles_with_sanitize("CCO", false).unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            rooted_at_atom: Some(1),
            ..Default::default()
        };

        assert_eq!(mol_to_smiles(&molecule, &params).unwrap(), "C(C)O");
    }

    #[test]
    fn isomeric_smiles_handles_tetrahedral_chirality() {
        // (R)-Alaninol: [C@@H](C)(N)CO  →  canonical output may differ
        // from input but `@`/`@@` marks must be present.
        let mut params = SmilesWriteParams::default();
        params.canonical = false;
        params.clean_stereo = false;
        let mol = Molecule::from_smiles_with_sanitize("C[C@@H](N)CO", false).unwrap();
        let smi = mol_to_smiles(&mol, &params).unwrap();
        assert!(
            smi.contains("@") || smi.contains("@@"),
            "chiral atom should produce @ or @@ mark, got: {smi}"
        );
        assert_eq!(smi, "C[C@@H](N)CO");
    }

    #[test]
    fn isomeric_smiles_handles_bond_stereo_direction() {
        let mut params = SmilesWriteParams::default();
        params.canonical = false;
        params.clean_stereo = false;
        let mol = Molecule::from_smiles_with_sanitize("C/C=C/C", false).unwrap();
        let smi = mol_to_smiles(&mol, &params).unwrap();
        assert_eq!(smi, "C/C=C/C");
    }

    #[test]
    fn isomeric_smiles_writes_double_bond_stereo_with_direction() {
        let mut params = SmilesWriteParams::default();
        params.canonical = false;
        params.clean_stereo = false;
        let mol = Molecule::from_smiles_with_sanitize("C/C=C/C", true).unwrap();
        let smi = mol_to_smiles(&mol, &params).unwrap();
        assert_eq!(smi, "C/C=C/C");
    }

    #[test]
    fn writer_canonicalizes_rooted_double_bond_directions_like_rdkit() {
        let mol = Molecule::from_smiles_with_sanitize("C/C=C/C", true).unwrap();
        let mut params = SmilesWriteParams {
            canonical: false,
            clean_stereo: false,
            rooted_at_atom: Some(1),
            ..Default::default()
        };

        assert_eq!(mol_to_smiles(&mol, &params).unwrap(), "C(/C)=C\\C");

        let z_mol = Molecule::from_smiles_with_sanitize("C/C=C\\C", true).unwrap();
        params.rooted_at_atom = Some(1);
        assert_eq!(mol_to_smiles(&z_mol, &params).unwrap(), "C(/C)=C/C");
    }

    #[test]
    fn writer_rooted_double_bond_directions_match_rdkit_without_sanitize() {
        let mut params = SmilesWriteParams {
            canonical: false,
            clean_stereo: false,
            rooted_at_atom: Some(1),
            ..Default::default()
        };

        let e_mol = Molecule::from_smiles_with_sanitize("C/C=C/C", false).unwrap();
        assert_eq!(mol_to_smiles(&e_mol, &params).unwrap(), "C(/C)=C\\C");

        let z_mol = Molecule::from_smiles_with_sanitize("C/C=C\\C", false).unwrap();
        params.rooted_at_atom = Some(1);
        assert_eq!(mol_to_smiles(&z_mol, &params).unwrap(), "C(/C)=C/C");
    }

    #[test]
    fn writer_rooted_terminal_double_bond_directions_match_rdkit_without_isomeric_smiles() {
        let mol = Molecule::from_smiles_with_sanitize("F/C=C\\F", true).unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            all_bonds_explicit: true,
            rooted_at_atom: Some(3),
            ..Default::default()
        };

        assert_eq!(mol_to_smiles(&mol, &params).unwrap(), "F/C=C\\F");
    }

    #[test]
    fn writer_preserves_rdkit_fused_ring_closure_digit_order_in_canonical_kekule_mode() {
        let mol = Molecule::from_smiles("Clc1ccc2ccc3cccc4ccc1c2c34").unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            do_kekule: true,
            canonical: true,
            clean_stereo: false,
            ..Default::default()
        };

        assert_eq!(
            mol_to_smiles(&mol, &params).unwrap(),
            "ClC1=C2C=CC3=CC=CC4=CC=C(C=C1)C2=C43"
        );
    }

    #[test]
    fn writer_kekule_all_hydrogens_explicit_preserves_pyrrolic_hydrogen_like_rdkit() {
        let mol = Molecule::from_smiles("[nH]1cccc1").unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            do_kekule: true,
            canonical: false,
            clean_stereo: false,
            all_hydrogens_explicit: true,
            ..Default::default()
        };

        assert_eq!(
            mol_to_smiles(&mol, &params).unwrap(),
            "[NH]1[CH]=[CH][CH]=[CH]1"
        );
    }

    #[test]
    fn writer_does_not_bracket_standard_aromatic_carbons_like_rdkit() {
        let mol = Molecule::from_smiles("c1ccccc1").unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            ..Default::default()
        };

        assert!(!atom_needs_bracket(&mol, AtomId::new(0), "", &params).unwrap());
        assert_eq!(mol_to_smiles(&mol, &params).unwrap(), "c1ccccc1");
    }

    #[test]
    fn writer_preserves_unbracketed_aromatic_ring_segment_like_rdkit_row_142() {
        let mol = Molecule::from_smiles("C12(C(C)c3ccccc3)NCC(C1(C)C)CC2").unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            do_kekule: false,
            canonical: false,
            clean_stereo: false,
            include_dative_bonds: false,
            ignore_atom_map_numbers: false,
            ..Default::default()
        };

        let smiles = mol_to_smiles(&mol, &params).unwrap();
        assert!(smiles.contains("c3ccccc3"), "got: {smiles}");
        assert!(!smiles.contains("[cH]"), "got: {smiles}");
    }

    #[test]
    fn writer_restores_atom_maps_before_kekulize_like_rdkit_row_142() {
        let input = "[C:12]12([CH:62]([CH3:65])[c:61]3[cH:64][cH:67][cH:68][cH:66][cH:63]3)[CH:20]4[c:30]5[c:40]6[c:49]7[c:57]8[c:60]([c:59]9[c:55]([c:47]([c:44]([c:52]9[c:51]([c:43]%10[c:35]%11[c:25]%12[c:19]%13%14)[c:53]8[c:45]%11[c:39]6[c:29]4%13)[c:34]([c:24]%15[c:15]%16[c:7]%17[c:3]%18%19)[c:33]%10[c:23]%16[c:16]%12[c:8]%18[c:11]%14[c:5]1%20)[c:37]([c:36]%21[c:26]%22[c:18]%23[c:10]%24[c:13]%25[c:6]%26%27)[c:27]%15[c:17]%22[c:9]%17[c:4]%24[c:1]%19[c:2]%20%26)[c:54]([c:46]%21[c:38]%28[c:28]%23[c:21]%25%29)[c:56]%30[c:48]%28[c:41]%31[c:31]%29[c:22]%32[c:14]2%27)[c:58]%30[c:50]7[c:42]%31[c:32]5%32";
        let mut molecule = Molecule::from_smiles(input).unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            do_kekule: true,
            canonical: true,
            clean_stereo: false,
            include_dative_bonds: false,
            ignore_atom_map_numbers: true,
            rooted_at_atom: Some(molecule.num_atoms() - 1),
            ..Default::default()
        };

        let saved_atom_maps = prepare_plain_smiles_molecule(&mut molecule, &params).unwrap();
        let plan = collect_fragment_write_plans(&molecule, &params)
            .unwrap()
            .into_iter()
            .next()
            .unwrap();
        let _ranks = rank_fragment_atoms_for_smiles(
            &molecule,
            &plan,
            &params,
            SmilesOutputMode::PlainSmiles,
        )
        .unwrap();
        restore_atom_maps_after_canonical_smiles(&mut molecule, saved_atom_maps.as_deref());
        let kekulized = kekulize_for_smiles(&molecule).unwrap();

        assert_eq!(kekulized.bonds()[11].order(), BondOrder::Double);
        assert_eq!(kekulized.bonds()[72].order(), BondOrder::Single);
    }

    #[test]
    fn writer_kekule_ignoring_dative_bonds_matches_rdkit_metal_porphyrin_branch() {
        let mol = Molecule::from_smiles(
            "O=C(O[Na])CC1=C(C(C(O[Na])=O)=C(C)C2=CC3=[N]4C(C(C=O)=C3CC)=CC5=C(C=C)C(C)=C6[N-]75)[N-]2[Cu+2]47[N](C8=C6)=C1C(C8C)CCC(O[Na])=O",
        )
        .unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            do_kekule: true,
            canonical: false,
            clean_stereo: false,
            include_dative_bonds: false,
            ..Default::default()
        };

        assert_eq!(
            mol_to_smiles(&mol, &params).unwrap(),
            "O=C([O][Na])CC1=C2[N]3[Cu+2]45[N]6C(=CC7=C(C)C(C([O][Na])=O)=C1[N-]74)C(CC)=C(C=O)C6=CC1=C(C=C)C(C)=C([N-]15)C=C3C(C)C2CCC([O][Na])=O"
        );
    }

    fn writer_cleans_nonstereo_double_bond_direction_specs_like_rdkit() {
        let mol = Molecule::from_smiles_with_sanitize("C/C=C(/C)C", true).unwrap();
        let params = SmilesWriteParams {
            canonical: true,
            clean_stereo: false,
            ..Default::default()
        };

        assert_eq!(mol_to_smiles(&mol, &params).unwrap(), "CC=C(C)C");
    }

    #[test]
    fn writer_canonical_nitro_start_atom_matches_rdkit_default_branch() {
        let mol = Molecule::from_smiles_with_sanitize("[O-][N+](=O)O", true).unwrap();
        let params = SmilesWriteParams::default();

        assert_eq!(mol_to_smiles(&mol, &params).unwrap(), "O=[N+]([O-])O");
    }

    #[test]
    fn writer_uses_cip_rank_ties_for_duplicate_double_bond_ligands_like_rdkit() {
        let params = SmilesWriteParams {
            canonical: true,
            clean_stereo: false,
            ..Default::default()
        };

        let fixtures = [
            ("C/C=C(/C(F))C(Cl)", "C/C=C(/CF)CCl"),
            ("C/C=C(/C(F))C(F)", "CC=C(CF)CF"),
            ("C/C=C(/CO)CN", "C/C=C(\\CN)CO"),
        ];

        for (input, expected) in fixtures {
            for sanitize in [true, false] {
                let mol = Molecule::from_smiles_with_sanitize(input, sanitize).unwrap();
                assert_eq!(
                    mol_to_smiles(&mol, &params).unwrap(),
                    expected,
                    "RDKit 2026.03.1 canonical output for {input} with sanitize={sanitize}"
                );
            }
        }
    }

    #[test]
    fn writer_canonicalizes_connected_double_bond_direction_queue_like_rdkit() {
        let mol = Molecule::from_smiles_with_sanitize("C/C=C/C=C\\C", true).unwrap();
        let params = SmilesWriteParams {
            canonical: true,
            clean_stereo: false,
            ..Default::default()
        };

        assert_eq!(mol_to_smiles(&mol, &params).unwrap(), "C/C=C\\C=C\\C");
    }

    #[test]
    fn writer_canonicalizes_connected_double_bond_direction_queue_without_sanitize() {
        let mol = Molecule::from_smiles_with_sanitize("C/C=C/C=C\\C", false).unwrap();
        let params = SmilesWriteParams {
            canonical: true,
            clean_stereo: false,
            ..Default::default()
        };

        assert_eq!(mol_to_smiles(&mol, &params).unwrap(), "C/C=C\\C=C\\C");
    }

    #[test]
    fn writer_canonicalizes_ring_double_bond_directions_like_rdkit() {
        let params = SmilesWriteParams {
            canonical: true,
            clean_stereo: false,
            ..Default::default()
        };

        let fixtures = [
            ("C1C/C=C/CCCCCCCC1", "C1=C/CCCCCCCCCC/1"),
            ("C1/C=C/C=C/CCCCCCCCC1", "C1=C/CCCCCCCCCC/C=C/1"),
            ("C1/C=C/CCCCC1", "C1=C/CCCCCC/1"),
            ("C1/C=C\\CCCCC1", "C1=C\\CCCCCC/1"),
        ];

        for (input, expected) in fixtures {
            for sanitize in [true, false] {
                let mol = Molecule::from_smiles_with_sanitize(input, sanitize).unwrap();
                assert_eq!(
                    mol_to_smiles(&mol, &params).unwrap(),
                    expected,
                    "RDKit 2026.03.1 canonical output for {input} with sanitize={sanitize}"
                );
            }
        }
    }

    #[test]
    fn writer_canonicalizes_fused_ring_double_bond_directions_like_rdkit() {
        let params = SmilesWriteParams {
            canonical: true,
            clean_stereo: false,
            ..Default::default()
        };

        let fixtures = [
            ("C1/C=C/C2CCCCC2C1", "C1=CC2CCCCC2CC1", "C1=C/C2CCCCC2CC/1"),
            (
                "C1/C=C\\C2CCCCC2C1",
                "C1=CC2CCCCC2CC1",
                "C1=C\\C2CCCCC2CC/1",
            ),
            (
                "C1C/C=C/CC2CCCCC12",
                "C1=CCC2CCCCC2CC1",
                "C1=C/CC2CCCCC2CC/1",
            ),
            (
                "C1/C=C/C2=C/CCCC2C1",
                "C1=CC2=CCCCC2CC1",
                "C1=C2/C=C/CCC2CCC1",
            ),
            (
                "C1/C=C\\C2=C/CCCC2C1",
                "C1=CC2=CCCCC2CC1",
                "C1=C2/C=C\\CCC2CCC1",
            ),
            (
                "C1/C=C/C=C/CC2CCCCC12",
                "C1=C/CC2CCCCC2C/C=C/1",
                "C1=C/CC2CCCCC2C/C=C/1",
            ),
            (
                "C1/C=C\\C=C/CC2CCCCC12",
                "C1=C\\CC2CCCCC2C\\C=C/1",
                "C1=C\\CC2CCCCC2C\\C=C/1",
            ),
            (
                "C1/C=C/C=C\\CC2CCCCC12",
                "C1=C\\CC2CCCCC2C/C=C/1",
                "C1=C\\CC2CCCCC2C/C=C/1",
            ),
        ];

        for (input, expected_sanitized, expected_unsanitized) in fixtures {
            for (sanitize, expected) in [(true, expected_sanitized), (false, expected_unsanitized)]
            {
                let mol = Molecule::from_smiles_with_sanitize(input, sanitize).unwrap();
                assert_eq!(
                    mol_to_smiles(&mol, &params).unwrap(),
                    expected,
                    "RDKit 2026.03.1 canonical output for {input} with sanitize={sanitize}"
                );
            }
        }
    }

    #[test]
    fn prepare_plain_smiles_molecule_initializes_fast_ring_info_for_fused_ring_stereo_like_rdkit() {
        let input = "C1/C=C/C2=C/CCCC2C1";
        let mut molecule = Molecule::from_smiles_with_sanitize(input, false).unwrap();
        let params = SmilesWriteParams {
            canonical: true,
            clean_stereo: false,
            ..Default::default()
        };

        let _saved_atom_maps = prepare_plain_smiles_molecule(&mut molecule, &params).unwrap();
        assert!(
            molecule
                .derived_cache()
                .rings
                .as_ref()
                .is_some_and(crate::RingInfo::is_find_fast_or_better)
        );
    }

    #[test]
    #[ignore = "debug helper for upstream/rust checkpoint alignment"]
    fn debug_probe_rust_writer_fused_ring_chain() {
        let input = "C1/C=C/C2=C/CCCC2C1";
        let params = SmilesWriteParams {
            canonical: true,
            clean_stereo: false,
            ..Default::default()
        };

        let focus = [0usize, 1usize, 2usize, 3usize, 4usize];
        let print_state = |name: &str, mol: &Molecule| {
            eprintln!(
                "checkpoint={name} rings_initialized={} is_symm_sssr={} stereo_done={}",
                mol.derived_cache().rings.is_some(),
                mol.derived_cache()
                    .rings
                    .as_ref()
                    .is_some_and(crate::RingInfo::is_symm_sssr),
                mol.prop("_StereochemDone").is_some()
            );
            for bond_idx in focus {
                let bond = &mol.bonds()[bond_idx];
                let ring_count = mol
                    .derived_cache()
                    .rings
                    .as_ref()
                    .map(|ri| ri.num_bond_rings(BondId::new(bond_idx)))
                    .unwrap_or(0);
                let min_ring = mol
                    .derived_cache()
                    .rings
                    .as_ref()
                    .map(|ri| ri.min_bond_ring_size(BondId::new(bond_idx)))
                    .unwrap_or(0);
                eprintln!(
                    "bond {} {}-{} dir={:?} stereo={:?} stereo_atoms={:?} ring_count={} min_ring={}",
                    bond_idx,
                    bond.begin().index(),
                    bond.end().index(),
                    bond.direction(),
                    bond.stereo(),
                    bond.stereo_atoms(),
                    ring_count,
                    min_ring
                );
            }
        };

        let mut raw = Molecule::from_smiles_with_sanitize(input, false).unwrap();
        print_state("raw_parse", &raw);

        update_property_cache_for_smiles(&mut raw).unwrap();
        print_state("post_update_property_cache", &raw);

        let mut ringed = raw.clone();
        ensure_fast_rings_for_writer_stereo_perception(&mut ringed).unwrap();
        print_state("post_fast_find_rings", &ringed);

        let mut assigned = ringed.clone();
        assign_stereochemistry_for_smiles(&mut assigned, false).unwrap();
        print_state("post_assign_stereochemistry", &assigned);

        eprintln!(
            "checkpoint=final_output smiles={}",
            mol_to_smiles(
                &Molecule::from_smiles_with_sanitize(input, false).unwrap(),
                &params
            )
            .unwrap()
        );
    }

    #[test]
    #[ignore = "debug helper for upstream/rust checkpoint alignment"]
    fn debug_probe_rust_writer_row_91() {
        let input = "O=C1/C(CC[C@@H](C)C1)=C(C)/C";
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            do_kekule: false,
            canonical: false,
            clean_stereo: false,
            all_bonds_explicit: true,
            all_hydrogens_explicit: false,
            include_dative_bonds: false,
            ignore_atom_map_numbers: false,
            ..Default::default()
        };

        let focus = [0usize, 1usize, 2usize, 3usize, 4usize];
        let print_state = |name: &str, mol: &Molecule| {
            eprintln!(
                "checkpoint={name} rings_initialized={} is_symm_sssr={} stereo_done={}",
                mol.derived_cache().rings.is_some(),
                mol.derived_cache()
                    .rings
                    .as_ref()
                    .is_some_and(crate::RingInfo::is_symm_sssr),
                mol.prop("_StereochemDone").is_some()
            );
            for bond_idx in focus {
                let bond = &mol.bonds()[bond_idx];
                let ring_count = mol
                    .derived_cache()
                    .rings
                    .as_ref()
                    .map(|ri| ri.num_bond_rings(BondId::new(bond_idx)))
                    .unwrap_or(0);
                let min_ring = mol
                    .derived_cache()
                    .rings
                    .as_ref()
                    .map(|ri| ri.min_bond_ring_size(BondId::new(bond_idx)))
                    .unwrap_or(0);
                eprintln!(
                    "bond {} {}-{} order={:?} dir={:?} stereo={:?} stereo_atoms={:?} ring_count={} min_ring={}",
                    bond_idx,
                    bond.begin().index(),
                    bond.end().index(),
                    bond.order(),
                    bond.direction(),
                    bond.stereo(),
                    bond.stereo_atoms(),
                    ring_count,
                    min_ring
                );
            }
        };

        let mut raw = Molecule::from_smiles_with_sanitize(input, false).unwrap();
        print_state("raw_parse", &raw);

        update_property_cache_for_smiles(&mut raw).unwrap();
        print_state("post_update_property_cache", &raw);

        let mut ringed = raw.clone();
        ensure_fast_rings_for_writer_stereo_perception(&mut ringed).unwrap();
        print_state("post_fast_find_rings", &ringed);

        let mut assigned = ringed.clone();
        assign_stereochemistry_for_smiles(&mut assigned, false).unwrap();
        print_state("post_assign_stereochemistry", &assigned);

        eprintln!(
            "checkpoint=final_output smiles={}",
            mol_to_smiles(
                &Molecule::from_smiles_with_sanitize(input, false).unwrap(),
                &params
            )
            .unwrap()
        );
    }

    #[test]
    #[ignore = "debug helper for upstream/rust checkpoint alignment"]
    fn debug_probe_rust_writer_row_91_traversal_state() {
        let input = "O=C1/C(CC[C@@H](C)C1)=C(C)/C";
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            do_kekule: false,
            canonical: false,
            clean_stereo: false,
            all_bonds_explicit: true,
            all_hydrogens_explicit: false,
            include_dative_bonds: false,
            ignore_atom_map_numbers: false,
            ..Default::default()
        };

        let mut molecule = Molecule::from_smiles_with_sanitize(input, false).unwrap();
        let _ = prepare_plain_smiles_molecule(&mut molecule, &params).unwrap();
        let plans = collect_fragment_write_plans(&molecule, &params).unwrap();
        let plan = &plans[0];
        let ranks =
            rank_fragment_atoms_for_smiles(&molecule, plan, &params, SmilesOutputMode::PlainSmiles)
                .unwrap();
        let start_atom = choose_fragment_start_atom(plan, &ranks, &params).unwrap();
        let traversal = canonicalize_fragment_stack(
            &molecule,
            plan,
            start_atom,
            &ranks,
            &params,
            SmilesWriteOverrides::default(),
        )
        .unwrap();
        eprintln!(
            "checkpoint=post_traversal bond1_dir={:?} bond1_stereo={:?} stack_len={}",
            molecule.bonds()[1].direction(),
            molecule.bonds()[1].stereo(),
            traversal.stack.len()
        );
        for item in &traversal.stack {
            match *item {
                MolStackElem::Bond(bond, atom_to_left) => {
                    if bond.index() == 1 {
                        eprintln!(
                            "stack_bond1 atom_to_left={} dir={:?} stereo={:?} stereo_atoms={:?}",
                            atom_to_left.index(),
                            molecule.bonds()[1].direction(),
                            molecule.bonds()[1].stereo(),
                            molecule.bonds()[1].stereo_atoms()
                        );
                    }
                }
                _ => {}
            }
        }
        let result = write_mol_stack(
            &molecule,
            &traversal.stack,
            &params,
            SmilesWriteOverrides::default(),
            &mut SmilesWriteContext::default(),
        )
        .unwrap();
        eprintln!("checkpoint=post_write smiles={}", result.smiles);
    }

    #[test]
    #[ignore = "debug helper for upstream/rust checkpoint alignment"]
    fn debug_probe_rust_writer_row_91_direction_canonicalization() {
        let input = "O=C1/C(CC[C@@H](C)C1)=C(C)/C";
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            do_kekule: false,
            canonical: false,
            clean_stereo: false,
            all_bonds_explicit: true,
            all_hydrogens_explicit: false,
            include_dative_bonds: false,
            ignore_atom_map_numbers: false,
            ..Default::default()
        };

        let mut molecule = Molecule::from_smiles_with_sanitize(input, false).unwrap();
        let _ = prepare_plain_smiles_molecule(&mut molecule, &params).unwrap();
        let plan = collect_fragment_write_plans(&molecule, &params)
            .unwrap()
            .remove(0);
        let ranks = rank_fragment_atoms_for_smiles(
            &molecule,
            &plan,
            &params,
            SmilesOutputMode::PlainSmiles,
        )
        .unwrap();
        let start_atom = choose_fragment_start_atom(&plan, &ranks, &params).unwrap();
        let traversal = canonicalize_fragment_stack(
            &molecule,
            &plan,
            start_atom,
            &ranks,
            &params,
            SmilesWriteOverrides::default(),
        )
        .unwrap();
        let mut canonicalized = molecule.clone();
        canonicalize_double_bond_directions_for_writer(
            &mut canonicalized,
            &traversal.stack,
            &traversal.traversal_ring_closure_bonds,
        )
        .unwrap();
        eprintln!(
            "checkpoint=post_direction_canonicalization bond1_dir={:?} bond7_dir={:?} bond8_dir={:?}",
            canonicalized.bonds()[1].direction(),
            canonicalized.bonds()[7].direction(),
            canonicalized.bonds()[8].direction()
        );
        eprintln!(
            "checkpoint=post_direction_canonicalization_smiles={}",
            write_mol_stack(
                &canonicalized,
                &traversal.stack,
                &params,
                SmilesWriteOverrides::default(),
                &mut SmilesWriteContext::default(),
            )
            .unwrap()
            .smiles
        );
    }

    #[test]
    fn canonicalize_double_bond_clear_bond_dirs_does_not_require_neighboring_stereo_double_bond() {
        let mut builder = crate::MoleculeBuilder::new();
        let center = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let anchor = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let fluorine = builder.add_atom(crate::AtomSpec::new(crate::Element::F));
        let chlorine = builder.add_atom(crate::AtomSpec::new(crate::Element::CL));
        let ref_bond = builder
            .add_bond(
                crate::BondSpec::new(center, fluorine, BondOrder::Single)
                    .with_direction(BondDirection::EndUpRight),
            )
            .unwrap();
        builder
            .add_bond(
                crate::BondSpec::new(center, chlorine, BondOrder::Single)
                    .with_direction(BondDirection::EndDownRight),
            )
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(center, anchor, BondOrder::Double))
            .unwrap();
        let mut molecule = builder.build().unwrap();

        let mut bond_dir_counts = vec![1, 1, 0];
        let mut atom_dir_counts = vec![2, 0, 2, 2];

        clear_bond_dirs_from_atom_for_writer(
            &mut molecule,
            ref_bond,
            center,
            &mut bond_dir_counts,
            &mut atom_dir_counts,
        );

        assert_eq!(bond_dir_counts, vec![1, 0, 0]);
        assert_eq!(atom_dir_counts, vec![1, 0, 2, 1]);
        assert_eq!(molecule.bonds()[1].direction(), BondDirection::None);
    }

    #[test]
    fn remove_redundant_bond_dir_specs_requires_neighboring_stereo_double_bond_like_rdkit() {
        let mut builder = crate::MoleculeBuilder::new();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let ref_bond = builder
            .add_bond(
                crate::BondSpec::new(a0, a1, BondOrder::Single)
                    .with_direction(BondDirection::EndUpRight),
            )
            .unwrap();
        let _dbl_bond = builder
            .add_bond(crate::BondSpec::new(a1, a2, BondOrder::Single))
            .unwrap();
        let mut molecule = builder.build().unwrap();
        let stack = vec![
            MolStackElem::Atom(a0),
            MolStackElem::Bond(ref_bond, a0),
            MolStackElem::Atom(a1),
        ];
        let mut bond_dir_counts = vec![1i8, 0i8];
        let mut atom_dir_counts = vec![2i8, 2i8, 2i8];

        remove_redundant_bond_dir_specs_for_writer(
            &mut molecule,
            &stack,
            &mut bond_dir_counts,
            &mut atom_dir_counts,
        );

        assert_eq!(molecule.bonds()[0].direction(), BondDirection::EndUpRight);
        assert_eq!(bond_dir_counts, vec![1, 0]);
        assert_eq!(atom_dir_counts, vec![2, 2, 2]);
    }

    #[test]
    fn writer_canonical_fragment_scope_preserves_aromatic_fused_ring_form_like_rdkit_row_94() {
        let mol = Molecule::from_smiles_with_sanitize(
            "Cl.Cl.COc1ccc2nccc([C@@H](O)[C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)c2c1",
            false,
        )
        .unwrap();
        let params = SmilesWriteParams::default();

        assert_eq!(
            mol_to_smiles(&mol, &params).unwrap(),
            "C=C[C@H]1CN2CC[C@H]1C[C@H]2[C@H](O)c1ccnc2ccc(OC)cc12.Cl.Cl"
        );
    }

    #[test]
    fn writer_plain_nonisomeric_row_91_matches_rdkit_after_direction_cleanup() {
        let mol =
            Molecule::from_smiles_with_sanitize("O=C1/C(CC[C@@H](C)C1)=C(C)/C", false).unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            do_kekule: false,
            canonical: false,
            clean_stereo: false,
            all_bonds_explicit: true,
            all_hydrogens_explicit: false,
            include_dative_bonds: false,
            ignore_atom_map_numbers: false,
            ..Default::default()
        };

        assert_eq!(
            mol_to_smiles(&mol, &params).unwrap(),
            "O=C1-C(=C(-C)-C)-C-C-C(-C)-C-1"
        );
    }

    #[test]
    #[ignore = "debug helper for upstream/rust checkpoint alignment"]
    fn debug_probe_rust_writer_row_91_direction_cleanup_substeps() {
        let input = "O=C1/C(CC[C@@H](C)C1)=C(C)/C";
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            do_kekule: false,
            canonical: false,
            clean_stereo: false,
            all_bonds_explicit: true,
            all_hydrogens_explicit: false,
            include_dative_bonds: false,
            ignore_atom_map_numbers: false,
            ..Default::default()
        };

        let mut molecule = Molecule::from_smiles_with_sanitize(input, false).unwrap();
        let _ = prepare_plain_smiles_molecule(&mut molecule, &params).unwrap();
        let plan = collect_fragment_write_plans(&molecule, &params)
            .unwrap()
            .remove(0);
        let ranks = rank_fragment_atoms_for_smiles(
            &molecule,
            &plan,
            &params,
            SmilesOutputMode::PlainSmiles,
        )
        .unwrap();
        let start_atom = choose_fragment_start_atom(&plan, &ranks, &params).unwrap();
        let traversal = canonicalize_fragment_stack(
            &molecule,
            &plan,
            start_atom,
            &ranks,
            &params,
            SmilesWriteOverrides::default(),
        )
        .unwrap();
        let mut canonicalized = molecule.clone();
        let mut atom_visit_orders = vec![usize::MAX; canonicalized.num_atoms()];
        let mut bond_visit_orders = vec![usize::MAX; canonicalized.num_bonds()];
        for (pos, item) in traversal.stack.iter().enumerate() {
            match *item {
                MolStackElem::Atom(atom) => atom_visit_orders[atom.index()] = pos,
                MolStackElem::Bond(bond, _) => {
                    bond_visit_orders[bond.index()] = pos;
                    if matches!(
                        canonicalized.bonds()[bond.index()].direction(),
                        BondDirection::EndDownRight | BondDirection::EndUpRight
                    ) {
                        canonicalized.topology_block_mut().bonds[bond.index()]
                            .set_direction(BondDirection::None);
                    }
                }
                MolStackElem::Ring { .. }
                | MolStackElem::BranchOpen
                | MolStackElem::BranchClose => {}
            }
        }

        let mut bond_dir_counts = vec![0i8; canonicalized.num_bonds()];
        let mut atom_dir_counts = vec![0i8; canonicalized.num_atoms()];
        let cip_ranks = crate::stereo::assign_atom_cip_ranks(&canonicalized).ok();
        canonicalize_double_bonds_for_writer(
            &mut canonicalized,
            &bond_visit_orders,
            &atom_visit_orders,
            &traversal.traversal_ring_closure_bonds,
            &mut bond_dir_counts,
            &mut atom_dir_counts,
            &traversal.stack,
            cip_ranks.as_deref(),
        );
        eprintln!(
            "checkpoint=after_double_bonds bond1_dir={:?} bond7_dir={:?} bond8_dir={:?} bond_dir_counts={:?} atom_dir_counts={:?}",
            canonicalized.bonds()[1].direction(),
            canonicalized.bonds()[7].direction(),
            canonicalized.bonds()[8].direction(),
            bond_dir_counts,
            atom_dir_counts
        );
        remove_unwanted_bond_dir_specs_for_writer(
            &mut canonicalized,
            &traversal.stack,
            &mut bond_dir_counts,
            &mut atom_dir_counts,
            &bond_visit_orders,
        );
        eprintln!(
            "checkpoint=after_remove_unwanted bond1_dir={:?} bond7_dir={:?} bond8_dir={:?} bond_dir_counts={:?} atom_dir_counts={:?}",
            canonicalized.bonds()[1].direction(),
            canonicalized.bonds()[7].direction(),
            canonicalized.bonds()[8].direction(),
            bond_dir_counts,
            atom_dir_counts
        );
        let bond1_begin = canonicalized.bonds()[1].begin();
        let bond1_end = canonicalized.bonds()[1].end();
        let bond8_end = canonicalized.bonds()[8].end();
        eprintln!(
            "checkpoint=before_manual_clear atom1_neighbors={:?} atom4_neighbors={:?} atom8_neighbors={:?}",
            incident_bonds(&canonicalized, bond1_begin)
                .iter()
                .map(|bond| (
                    bond.index(),
                    canonicalized.bonds()[bond.index()].order(),
                    canonicalized.bonds()[bond.index()].direction()
                ))
                .collect::<Vec<_>>(),
            incident_bonds(&canonicalized, bond1_end)
                .iter()
                .map(|bond| (
                    bond.index(),
                    canonicalized.bonds()[bond.index()].order(),
                    canonicalized.bonds()[bond.index()].direction()
                ))
                .collect::<Vec<_>>(),
            incident_bonds(&canonicalized, bond8_end)
                .iter()
                .map(|bond| (
                    bond.index(),
                    canonicalized.bonds()[bond.index()].order(),
                    canonicalized.bonds()[bond.index()].direction()
                ))
                .collect::<Vec<_>>()
        );
        clear_bond_dirs_from_atom_for_writer(
            &mut canonicalized,
            BondId::new(1),
            bond1_begin,
            &mut bond_dir_counts,
            &mut atom_dir_counts,
        );
        eprintln!(
            "checkpoint=after_manual_clear_bond1_begin bond1_dir={:?} bond7_dir={:?} bond8_dir={:?} bond_dir_counts={:?} atom_dir_counts={:?}",
            canonicalized.bonds()[1].direction(),
            canonicalized.bonds()[7].direction(),
            canonicalized.bonds()[8].direction(),
            bond_dir_counts,
            atom_dir_counts
        );
        clear_bond_dirs_from_atom_for_writer(
            &mut canonicalized,
            BondId::new(1),
            bond1_end,
            &mut bond_dir_counts,
            &mut atom_dir_counts,
        );
        eprintln!(
            "checkpoint=after_manual_clear_bond1_end bond1_dir={:?} bond7_dir={:?} bond8_dir={:?} bond_dir_counts={:?} atom_dir_counts={:?}",
            canonicalized.bonds()[1].direction(),
            canonicalized.bonds()[7].direction(),
            canonicalized.bonds()[8].direction(),
            bond_dir_counts,
            atom_dir_counts
        );
        remove_redundant_bond_dir_specs_for_writer(
            &mut canonicalized,
            &traversal.stack,
            &mut bond_dir_counts,
            &mut atom_dir_counts,
        );
        eprintln!(
            "checkpoint=after_remove_redundant bond1_dir={:?} bond7_dir={:?} bond8_dir={:?} bond_dir_counts={:?} atom_dir_counts={:?}",
            canonicalized.bonds()[1].direction(),
            canonicalized.bonds()[7].direction(),
            canonicalized.bonds()[8].direction(),
            bond_dir_counts,
            atom_dir_counts
        );
        eprintln!(
            "checkpoint=post_write smiles={}",
            write_mol_stack(
                &canonicalized,
                &traversal.stack,
                &params,
                SmilesWriteOverrides::default(),
                &mut SmilesWriteContext::default(),
            )
            .unwrap()
            .smiles
        );
    }

    #[test]
    fn writer_rejects_invalid_nontetrahedral_permutation_with_structured_error() {
        let error = validate_writer_chiral_permutation(ChiralTag::SquarePlanar, 4).unwrap_err();

        assert_eq!(
            error,
            SmilesWriteError::InvalidChiralPermutation {
                chiral_tag: ChiralTag::SquarePlanar,
                permutation: 4,
                limit: 3,
            }
        );
    }

    #[test]
    fn writer_start_atom_guard_reports_invariant_for_empty_rank_scope() {
        let plan = FragmentWritePlan {
            atoms: Vec::new(),
            bonds: Vec::new(),
            rooted_at_atom: None,
        };
        let params = SmilesWriteParams::default();

        let error = choose_fragment_start_atom(&plan, &[], &params).unwrap_err();

        assert_eq!(
            error,
            SmilesWriteError::InvariantViolation {
                stage: "ShortTermAtomWriter",
                message: "choose_fragment_start_atom() called with empty canonical rank scope",
            }
        );
    }

    #[test]
    fn writer_swap_counter_reports_invariant_for_length_mismatch() {
        let error =
            count_swaps_to_interconvert(&[BondId::new(0)], &[BondId::new(0), BondId::new(1)])
                .unwrap_err();

        assert_eq!(
            error,
            SmilesWriteError::InvariantViolation {
                stage: "ShortTermAtomWriter",
                message: "count_swaps_to_interconvert() probe/reference length mismatch",
            }
        );
    }

    #[test]
    fn writer_swap_counter_reports_invariant_for_missing_expected_bond() {
        let error = count_swaps_to_interconvert(
            &[BondId::new(0), BondId::new(1)],
            &[BondId::new(0), BondId::new(2)],
        )
        .unwrap_err();

        assert_eq!(
            error,
            SmilesWriteError::InvariantViolation {
                stage: "ShortTermAtomWriter",
                message: "count_swaps_to_interconvert() expected bond missing from probe order",
            }
        );
    }

    #[test]
    fn writer_bond_guard_uses_invariant_error_shape() {
        let error = invariant_stage_error::<()>(
            SmilesPlanStage::ShortTermBondWriter,
            "write_ring_closure() could not allocate a free ring index",
        )
        .unwrap_err();

        assert_eq!(
            error,
            SmilesWriteError::InvariantViolation {
                stage: "ShortTermBondWriter",
                message: "write_ring_closure() could not allocate a free ring index",
            }
        );
    }

    #[test]
    fn atom_bond_and_fragment_writer_entries_fail_closed_until_ported() {
        let molecule = ethane();

        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            ..Default::default()
        };

        assert_eq!(get_atom_smiles(&molecule, 0, &params).unwrap(), "C");
        assert_eq!(
            get_molecule_bond_smiles(&molecule, 0, Some(0), &params).unwrap(),
            ""
        );
        assert_eq!(get_bond_smiles(BondOrder::Single).unwrap(), "");
        assert_eq!(get_bond_smiles(BondOrder::Dative).unwrap(), "->");
        // Fragment API is now implemented — ethane fragment produces SMILES
        let fragment = mol_fragment_to_smiles(
            &molecule,
            &SmilesWriteParams::default(),
            &[0, 1],
            None,
            None,
            None,
        )
        .unwrap();
        assert_eq!(fragment, "CC", "ethane fragment should produce CC");

        let fragment_cx = mol_fragment_to_cx_smiles(
            &molecule,
            &SmilesWriteParams::default(),
            &[0, 1],
            None,
            None,
            None,
            CxSmilesFields::ALL,
        )
        .unwrap();
        assert_eq!(fragment_cx, "CC", "ethane fragment CX should be plain CC");
    }

    #[test]
    fn atom_bond_and_fragment_writer_entries_validate_indices_before_unsupported() {
        let molecule = ethane();

        assert_eq!(
            get_atom_smiles(&molecule, 2, &SmilesWriteParams::default()).unwrap_err(),
            SmilesWriteError::AtomOutOfRange { atom: 2 }
        );
        assert_eq!(
            get_molecule_bond_smiles(&molecule, 1, None, &SmilesWriteParams::default())
                .unwrap_err(),
            SmilesWriteError::BondOutOfRange { bond: 1 }
        );
        assert_eq!(
            mol_fragment_to_smiles(
                &molecule,
                &SmilesWriteParams::default(),
                &[2],
                None,
                None,
                None
            )
            .unwrap_err(),
            SmilesWriteError::AtomOutOfRange { atom: 2 }
        );
    }

    #[test]
    fn mol_fragment_to_smiles_uses_original_atom_and_bond_symbol_overrides() {
        let molecule = Molecule::from_smiles_with_sanitize("CCO", false).unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            ..Default::default()
        };
        let atom_symbols = vec!["A".to_string(), "B".to_string(), "C".to_string()];
        let bond_symbols = vec!["~".to_string(), "!".to_string()];

        let fragment = mol_fragment_to_smiles(
            &molecule,
            &params,
            &[0, 1],
            Some(&[0]),
            Some(&atom_symbols),
            Some(&bond_symbols),
        )
        .unwrap();

        assert_eq!(fragment, "A~B");
    }

    #[test]
    fn mol_fragment_to_smiles_splits_disconnected_atoms_within_fragment_scope() {
        let molecule = Molecule::from_smiles_with_sanitize("CCO", false).unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            ..Default::default()
        };

        let fragment =
            mol_fragment_to_smiles(&molecule, &params, &[0, 2], None, None, None).unwrap();

        assert_eq!(fragment, "C.O");
    }

    #[test]
    fn mol_fragment_to_cx_smiles_filters_cx_blocks_to_fragment_output_scope() {
        let mut builder = crate::MoleculeBuilder::new();
        let c0 = builder.add_atom(
            crate::AtomSpec::new(crate::Element::C).with_prop("_supplementalSmilesLabel", "keep"),
        );
        let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let o2 = builder.add_atom(
            crate::AtomSpec::new(crate::Element::O).with_prop("_supplementalSmilesLabel", "drop"),
        );
        builder
            .add_bond(crate::BondSpec::new(c0, c1, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(c1, o2, BondOrder::Dative))
            .unwrap();
        let molecule = builder.build().unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            ..Default::default()
        };

        let fragment = mol_fragment_to_cx_smiles(
            &molecule,
            &params,
            &[0, 1],
            Some(&[0]),
            None,
            None,
            CxSmilesFields::ATOM_LABELS | CxSmilesFields::COORDINATE_BONDS,
        )
        .unwrap();

        assert!(
            fragment.contains("keep"),
            "fragment CX should keep atom 0 label: {fragment}"
        );
        assert!(
            !fragment.contains("drop"),
            "fragment CX must not leak atom 2 label: {fragment}"
        );
        assert!(
            !fragment.contains("_Z:2"),
            "fragment CX must not leak out-of-scope dative bond: {fragment}"
        );
    }

    // ── CX SMILES Extension Tests ──────────────────────────────────────────

    fn cx_ethanol() -> Molecule {
        // CCO with atom label on O
        let mut builder = crate::MoleculeBuilder::new();
        let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let c2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let o = builder
            .add_atom(crate::AtomSpec::new(crate::Element::O).with_prop("atomLabel", "Hydroxy"));
        builder
            .add_bond(crate::BondSpec::new(c1, c2, crate::BondOrder::Single))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(c2, o, crate::BondOrder::Single))
            .unwrap();
        builder.build().unwrap()
    }

    #[test]
    fn cx_individual_coords_writes_atom_order_coordinates() {
        let mut builder = crate::MoleculeBuilder::new();
        let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let c2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(c1, c2, crate::BondOrder::Single))
            .unwrap();
        builder
            .set_2d_coordinates(vec![[0.0, 0.0], [1.5, 0.0]])
            .unwrap();
        let molecule = builder.build().unwrap();
        let scope = CxWriteScope::full_molecule(&molecule);
        let coords_str = write_cx_coords(&molecule, &scope.atom_order);
        assert_eq!(coords_str, "0,0,;1.5,0,");
    }

    #[test]
    fn cx_empty_coords_when_no_2d_coords_present() {
        let molecule = ethane();
        let scope = CxWriteScope::full_molecule(&molecule);
        assert_eq!(write_cx_coords(&molecule, &scope.atom_order), "");
    }

    #[test]
    fn cx_individual_atom_labels_writes_label_entries() {
        let molecule = cx_ethanol();
        let scope = CxWriteScope::full_molecule(&molecule);
        let labels = write_cx_atom_labels(&molecule, &scope.atom_order);
        assert_eq!(labels, ";;Hydroxy");
    }

    #[test]
    fn cx_individual_radicals_writes_entries() {
        // Build a molecule with a radical
        let mut builder = crate::MoleculeBuilder::new();
        let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let c2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(c1, c2, crate::BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();
        // Test write_cx_radicals directly (bypasses the SMILES writer which
        // rejects radical atoms)
        let scope = CxWriteScope::full_molecule(&mol);
        let radicals = write_cx_radicals(&mol, &scope.atom_order);
        // No radicals on plain ethane
        assert_eq!(radicals, "");
    }

    #[test]
    fn cx_no_radicals_when_none_present() {
        let molecule = ethane();
        let scope = CxWriteScope::full_molecule(&molecule);
        assert_eq!(write_cx_radicals(&molecule, &scope.atom_order), "");
    }

    #[test]
    fn cx_atom_props_writes_entries_for_atoms_with_properties() {
        let molecule = ethane();
        let scope = CxWriteScope::full_molecule(&molecule);
        let props = write_cx_atom_props(&molecule, &scope.atom_order);
        assert_eq!(props, "");
    }

    #[test]
    fn cx_coordinate_bonds_writes_entries() {
        // Build a molecule with a dative bond
        let mut builder = crate::MoleculeBuilder::new();
        let n = builder.add_atom(crate::AtomSpec::new(crate::Element::N));
        let o = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        builder
            .add_bond(crate::BondSpec::new(n, o, crate::BondOrder::Dative))
            .unwrap();
        let mol = builder.build().unwrap();

        let scope = CxWriteScope::full_molecule(&mol);
        let coord_bonds =
            write_cx_coordinate_bonds(&mol, &scope.atom_order, &scope.bond_order, "C");
        assert_eq!(coord_bonds, "C:0.0");
    }

    #[test]
    fn cx_bond_cfg_writes_wedge_and_dash_when_coords_are_included() {
        let mut builder = crate::MoleculeBuilder::new();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(
                crate::BondSpec::new(a0, a1, crate::BondOrder::Single)
                    .with_direction(crate::BondDirection::BeginWedge),
            )
            .unwrap();
        builder
            .add_bond(
                crate::BondSpec::new(a1, a2, crate::BondOrder::Single)
                    .with_direction(crate::BondDirection::BeginDash),
            )
            .unwrap();
        builder
            .set_2d_coordinates(vec![[0.0, 0.0], [1.5, 0.0], [3.0, 0.0]])
            .unwrap();
        let molecule = builder.build().unwrap();

        let extension =
            get_cx_extensions(&molecule, CxSmilesFields::COORDS | CxSmilesFields::BOND_CFG)
                .unwrap();
        assert!(
            extension.contains("wU:0.0"),
            "wedge-up bond config should be emitted, got: {extension:?}"
        );
        assert!(
            extension.contains("wD:1.1"),
            "wedge-down bond config should be emitted, got: {extension:?}"
        );
    }

    #[test]
    fn cx_bond_cfg_uses_molfile_bond_cfg_fallback_for_wedge_and_dash() {
        let mut builder = crate::MoleculeBuilder::new();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(
                crate::BondSpec::new(a0, a1, crate::BondOrder::Single)
                    .with_prop("_MolFileBondCfg", "3"),
            )
            .unwrap();
        builder
            .set_2d_coordinates(vec![[0.0, 0.0], [1.5, 0.0]])
            .unwrap();
        let molecule = builder.build().unwrap();

        let extension =
            get_cx_extensions(&molecule, CxSmilesFields::COORDS | CxSmilesFields::BOND_CFG)
                .unwrap();
        assert!(
            extension.contains("wD:0.0"),
            "MolFile cfg=3 should emit dash wedge config, got: {extension:?}"
        );
    }

    #[test]
    fn cx_bond_cfg_emits_unknown_without_coords_from_molfile_cfg() {
        let mut builder = crate::MoleculeBuilder::new();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(
                crate::BondSpec::new(a0, a1, crate::BondOrder::Single)
                    .with_prop("_MolFileBondCfg", "2"),
            )
            .unwrap();
        let molecule = builder.build().unwrap();

        let extension = get_cx_extensions(&molecule, CxSmilesFields::BOND_CFG).unwrap();
        assert_eq!(extension, "|w:0.0|");
    }

    #[test]
    fn cx_bond_atropisomer_writes_wedge_for_atrop_neighbor() {
        let mut builder = crate::MoleculeBuilder::new();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(
                crate::BondSpec::new(a0, a1, crate::BondOrder::Single)
                    .with_direction(crate::BondDirection::BeginWedge),
            )
            .unwrap();
        builder
            .add_bond(
                crate::BondSpec::new(a0, a2, crate::BondOrder::Single)
                    .with_stereo(crate::BondStereo::AtropCw),
            )
            .unwrap();
        let molecule = builder.build().unwrap();

        let extension = get_cx_extensions(&molecule, CxSmilesFields::BOND_ATROPISOMER).unwrap();
        assert_eq!(extension, "|wU:0.0|");
    }

    #[test]
    fn cx_bond_atropisomer_does_not_use_molfile_cfg_fallback() {
        let mut builder = crate::MoleculeBuilder::new();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(
                crate::BondSpec::new(a0, a1, crate::BondOrder::Single)
                    .with_prop("_MolFileBondCfg", "1"),
            )
            .unwrap();
        builder
            .add_bond(
                crate::BondSpec::new(a0, a2, crate::BondOrder::Single)
                    .with_stereo(crate::BondStereo::AtropCw),
            )
            .unwrap();
        let molecule = builder.build().unwrap();

        let extension = get_cx_extensions(&molecule, CxSmilesFields::BOND_ATROPISOMER).unwrap();
        assert_eq!(extension, "");
    }

    fn cx_ring_double_bond_molecule(stereo: BondStereo) -> Molecule {
        let mut builder = crate::MoleculeBuilder::new();
        let atoms: Vec<_> = (0..8)
            .map(|_| builder.add_atom(crate::AtomSpec::new(crate::Element::C)))
            .collect();
        builder
            .add_bond(
                crate::BondSpec::new(atoms[0], atoms[1], crate::BondOrder::Double)
                    .with_stereo(stereo)
                    .with_stereo_atoms(atoms[7], atoms[2]),
            )
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(
                atoms[1],
                atoms[2],
                crate::BondOrder::Single,
            ))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(
                atoms[2],
                atoms[3],
                crate::BondOrder::Single,
            ))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(
                atoms[3],
                atoms[4],
                crate::BondOrder::Single,
            ))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(
                atoms[4],
                atoms[5],
                crate::BondOrder::Single,
            ))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(
                atoms[5],
                atoms[6],
                crate::BondOrder::Single,
            ))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(
                atoms[6],
                atoms[7],
                crate::BondOrder::Single,
            ))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(
                atoms[7],
                atoms[0],
                crate::BondOrder::Single,
            ))
            .unwrap();
        let mut molecule = builder.build().unwrap();
        let ring_info = crate::symmetrize_sssr(&molecule).unwrap();
        molecule.derived_cache_mut().rings = Some(ring_info);
        molecule
    }

    #[test]
    fn write_cx_ringbond_cistrans_block_writes_c_block_for_cis_ring_double_bond() {
        let molecule = cx_ring_double_bond_molecule(BondStereo::Cis);
        let extension = get_cx_extensions(&molecule, CxSmilesFields::BOND_CFG).unwrap();
        assert!(
            extension.contains("c:0"),
            "cis ring double bond should emit c block, got: {extension:?}"
        );
    }

    #[test]
    fn write_cx_ringbond_cistrans_block_writes_ctu_block_for_any_ring_double_bond() {
        let molecule = cx_ring_double_bond_molecule(BondStereo::Any);
        let extension = get_cx_extensions(&molecule, CxSmilesFields::BOND_CFG).unwrap();
        assert!(
            extension.contains("ctu:0"),
            "any ring double bond should emit ctu block, got: {extension:?}"
        );
    }

    #[test]
    fn cx_linknodes_block_writes_compact_form_for_degree_two_center() {
        let mut builder = crate::MoleculeBuilder::new();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(a0, a1, crate::BondOrder::Single))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(a1, a2, crate::BondOrder::Single))
            .unwrap();
        let molecule = builder
            .build()
            .unwrap()
            .with_prop("_MolFileLinkNodes", "1 3 2 2 1 2 3");

        let scope = CxWriteScope::full_molecule(&molecule);
        let block = write_cx_linknodes_block(&molecule, &scope.atom_order);
        assert_eq!(block, "LN:1:1.3");
    }

    #[test]
    fn cx_linknodes_block_writes_outer_atoms_for_degree_three_center() {
        let mut builder = crate::MoleculeBuilder::new();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a3 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(a0, a1, crate::BondOrder::Single))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(a0, a2, crate::BondOrder::Single))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(a0, a3, crate::BondOrder::Single))
            .unwrap();
        let molecule = builder
            .build()
            .unwrap()
            .with_prop("_MolFileLinkNodes", "1 4 2 1 2 1 3");

        let scope = CxWriteScope::full_molecule(&molecule);
        let block = write_cx_linknodes_block(&molecule, &scope.atom_order);
        assert_eq!(block, "LN:0:1.4.1.2");
    }

    #[test]
    fn cx_polymer_block_writes_sru_label_and_connect() {
        let mut builder = crate::MoleculeBuilder::new();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(a0, a1, crate::BondOrder::Single))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(a1, a2, crate::BondOrder::Single))
            .unwrap();
        builder
            .add_substance_group(
                crate::SubstanceGroup::new(
                    crate::SubstanceGroupId::new(0),
                    crate::SubstanceGroupKind::StructuralRepeatUnit,
                )
                .with_atoms(vec![a0, a1])
                .with_label("SRU")
                .with_connection(crate::SGroupConnection::HeadToHead),
            )
            .unwrap();
        let molecule = builder.build().unwrap();

        let extension = get_cx_extensions(&molecule, CxSmilesFields::POLYMER).unwrap();
        assert_eq!(extension, "|Sg:n:0,1:SRU:hh:::|");
    }

    #[test]
    fn cx_polymer_block_writes_copolymer_subtype_and_crossings() {
        let mut builder = crate::MoleculeBuilder::new();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(a0, a1, crate::BondOrder::Single))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(a1, a2, crate::BondOrder::Single))
            .unwrap();
        builder
            .add_substance_group(
                crate::SubstanceGroup::new(
                    crate::SubstanceGroupId::new(0),
                    crate::SubstanceGroupKind::Copolymer,
                )
                .with_atoms(vec![a0, a1])
                .with_subtype("ALT")
                .with_label("COP")
                .with_connection(crate::SGroupConnection::Either)
                .with_prop("XBHEAD", "1,0")
                .with_prop("XBCORR", "0,1,0,1"),
            )
            .unwrap();
        let molecule = builder.build().unwrap();

        let extension = get_cx_extensions(&molecule, CxSmilesFields::POLYMER).unwrap();
        assert_eq!(extension, "|Sg:alt:0,1:COP:eu:1,0:1,1:|");
    }

    #[test]
    fn cx_sgroup_hierarchy_writes_block_for_polymer_parent_child() {
        let mut builder = crate::MoleculeBuilder::new();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(a0, a1, crate::BondOrder::Single))
            .unwrap();
        builder
            .add_bond(crate::BondSpec::new(a1, a2, crate::BondOrder::Single))
            .unwrap();
        builder
            .add_substance_group(
                crate::SubstanceGroup::new(
                    crate::SubstanceGroupId::new(0),
                    crate::SubstanceGroupKind::StructuralRepeatUnit,
                )
                .with_atoms(vec![a0, a1]),
            )
            .unwrap();
        builder
            .add_substance_group(
                crate::SubstanceGroup::new(
                    crate::SubstanceGroupId::new(1),
                    crate::SubstanceGroupKind::StructuralRepeatUnit,
                )
                .with_atoms(vec![a1, a2])
                .with_parent(crate::SubstanceGroupId::new(0)),
            )
            .unwrap();
        let molecule = builder.build().unwrap();

        let extension = get_cx_extensions(&molecule, CxSmilesFields::POLYMER).unwrap();
        assert!(
            extension.contains("SgH:0:1"),
            "polymer hierarchy should be emitted, got: {extension:?}"
        );
    }

    #[test]
    fn cx_sgroup_hierarchy_writes_block_for_sgroup_parent_prop() {
        let mut builder = crate::MoleculeBuilder::new();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(a0, a1, crate::BondOrder::Single))
            .unwrap();
        builder
            .add_substance_group(
                crate::SubstanceGroup::new(
                    crate::SubstanceGroupId::new(0),
                    crate::SubstanceGroupKind::Data,
                )
                .with_atoms(vec![a0])
                .with_prop("index", "0"),
            )
            .unwrap();
        builder
            .add_substance_group(
                crate::SubstanceGroup::new(
                    crate::SubstanceGroupId::new(1),
                    crate::SubstanceGroupKind::Data,
                )
                .with_atoms(vec![a1])
                .with_prop("index", "1")
                .with_prop("PARENT", "0"),
            )
            .unwrap();
        let molecule = builder.build().unwrap();

        let extension = get_cx_extensions(&molecule, CxSmilesFields::SGROUPS).unwrap();
        assert!(
            extension.contains("SgH:0:1"),
            "SGroup hierarchy should be emitted for SGROUPS fields, got: {extension:?}"
        );
    }

    #[test]
    fn cx_full_molecule_with_all_fields_returns_smiles_with_cx_extension() {
        let molecule = cx_ethanol();
        let result = mol_to_cx_smiles(
            &molecule,
            &SmilesWriteParams::default(),
            CxSmilesFields::ATOM_LABELS,
            RestoreBondDirOption::None,
        )
        .unwrap();
        // Should contain the SMILES and CX extension separated by space
        assert!(!result.is_empty(), "should produce output");
        assert!(
            result.contains("Hydroxy"),
            "CX label should appear: {result:?}"
        );
    }

    #[test]
    fn cx_no_cx_data_returns_plain_smiles() {
        let molecule = ethane();
        let result = mol_to_cx_smiles(
            &molecule,
            &SmilesWriteParams::default(),
            CxSmilesFields::ALL,
            RestoreBondDirOption::Clear,
        )
        .unwrap();
        assert_eq!(result, "CC", "ethane with no CX data should be plain CC");
    }

    #[test]
    fn cx_get_cx_extensions_returns_union_of_requested_fields() {
        let molecule = cx_ethanol();
        // Request ATOM_LABELS
        let result = get_cx_extensions(&molecule, CxSmilesFields::ATOM_LABELS).unwrap();
        assert!(result.contains("Hydroxy"), "should have atom labels");
        assert!(
            !result.contains('('),
            "no coords (molecule has no 2D coords)"
        );
    }

    #[test]
    fn cx_radical_fields_only_produces_expected_output() {
        let molecule = cx_ethanol();
        let result = get_cx_extensions(&molecule, CxSmilesFields::RADICALS).unwrap();
        // cx_ethanol has no radicals, so result should be empty
        assert_eq!(result, "", "no radicals expected: {result:?}");
    }

    #[test]
    fn cx_individual_atom_props_escapes_property_name_and_value_dots() {
        let mut builder = crate::MoleculeBuilder::new();
        let c = builder
            .add_atom(crate::AtomSpec::new(crate::Element::C).with_prop("foo.bar", "baz.qux"));
        let c2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(c, c2, crate::BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();

        let scope = CxWriteScope::full_molecule(&mol);
        let props = write_cx_atom_props(&mol, &scope.atom_order);
        assert_eq!(props, "atomProp:0.foo&#46;bar.baz&#46;qux");
    }

    #[test]
    fn cx_molfile_values_when_present() {
        let mut builder = crate::MoleculeBuilder::new();
        let c = builder.add_atom(
            crate::AtomSpec::new(crate::Element::C).with_prop("molFileValue", "test_value"),
        );
        let c2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(c, c2, crate::BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();

        let scope = CxWriteScope::full_molecule(&mol);
        let values = write_cx_molfile_values(&mol, &scope.atom_order);
        assert_eq!(values, "test_value;");
    }

    #[test]
    fn cx_stereo_group_writes_appropriate_code() {
        let mut builder = crate::MoleculeBuilder::new();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(a0, a1, crate::BondOrder::Single))
            .unwrap();
        builder
            .add_stereo_group(crate::StereoGroup::new(
                crate::StereoGroupKind::Or,
                vec![a1],
                vec![],
            ))
            .unwrap();
        builder
            .add_stereo_group(crate::StereoGroup::new(
                crate::StereoGroupKind::And,
                vec![a0],
                vec![],
            ))
            .unwrap();
        let molecule = builder.build().unwrap();
        let scope = CxWriteScope::full_molecule(&molecule);
        let stereo = write_cx_enhanced_stereo(&molecule, &scope.atom_order, &scope.bond_order);
        assert_eq!(stereo, "o1:1,&1:0");
    }

    #[test]
    fn cleanup_stereo_groups_for_cx_smiles_moves_atrop_atoms_to_bond_membership() {
        let mut builder = crate::MoleculeBuilder::new();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(a0, a1, crate::BondOrder::Single))
            .unwrap();
        builder
            .add_bond(
                crate::BondSpec::new(a1, a2, crate::BondOrder::Single)
                    .with_stereo(crate::BondStereo::AtropCw),
            )
            .unwrap();
        builder
            .add_stereo_group(crate::StereoGroup::new(
                crate::StereoGroupKind::Or,
                vec![a1],
                Vec::new(),
            ))
            .unwrap();
        let mut molecule = builder.build().unwrap();

        cleanup_stereo_groups_for_cx_smiles(&mut molecule).unwrap();

        assert_eq!(molecule.stereo_groups().len(), 1);
        let group = &molecule.stereo_groups()[0];
        assert!(group.atoms().is_empty());
        assert_eq!(group.bonds(), &[BondId::new(1)]);
    }

    #[test]
    fn cx_sgroups_empty_when_no_sgroups() {
        let molecule = ethane();
        let scope = CxWriteScope::full_molecule(&molecule);
        let sgroups = write_cx_sgroups(&molecule, &scope.atom_order, &scope.bond_order);
        assert_eq!(sgroups, "");
    }

    #[test]
    fn cx_data_sgroup_prefers_typed_data_fields() {
        let mut builder = crate::MoleculeBuilder::new();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(a0, a1, crate::BondOrder::Single))
            .unwrap();
        builder
            .add_substance_group(
                crate::SubstanceGroup::new(
                    crate::SubstanceGroupId::new(0),
                    crate::SubstanceGroupKind::Data,
                )
                .with_atoms(vec![a1, a0])
                .with_data(crate::SGroupData {
                    field_name: Some("FIELD".to_string()),
                    field_info: Some("INFO".to_string()),
                    query_op: Some("OP".to_string()),
                    units: Some("TAG".to_string()),
                    values: vec!["VALUE".to_string()],
                    ..Default::default()
                }),
            )
            .unwrap();
        let molecule = builder.build().unwrap();
        let scope = CxWriteScope::full_molecule(&molecule);

        assert_eq!(
            write_cx_data_sgroups(&molecule, &scope.atom_order),
            "SgD:1,0:FIELD:VALUE:OP:INFO:TAG:"
        );
    }

    #[test]
    fn build_smiles_helpers_match_main_writer_helpers() {
        let mut builder = crate::MoleculeBuilder::new();
        let c = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let n = builder.add_atom(crate::AtomSpec::new(crate::Element::N));
        builder
            .add_bond(crate::BondSpec::new(c, n, crate::BondOrder::Triple))
            .unwrap();
        let molecule = builder.build().unwrap();
        let params = SmilesWriteParams::default();
        let context = SmilesWriteContext::default();

        assert_eq!(
            build_atom_smiles(&molecule, c, &params, &context).unwrap(),
            get_atom_smiles(&molecule, c.index(), &params).unwrap()
        );
        assert_eq!(
            build_bond_smiles(&molecule, BondId::new(0), c, &params).unwrap(),
            get_molecule_bond_smiles(&molecule, 0, Some(c.index()), &params).unwrap()
        );
    }

    #[test]
    fn writer_rooted_ring_stereo_centers_follow_rdkit_canonicalize_fragment_postprocessing() {
        let molecule = Molecule::from_smiles(
            "O=C(NC[C@]12C[C@H]3C[C@H](C[C@H](C3)C1)C2)[C@@H]1C[C@H]2c3ccccc3[C@@H]1c1ccccc12",
        )
        .unwrap();
        let params = SmilesWriteParams {
            canonical: false,
            clean_stereo: false,
            rooted_at_atom: Some(molecule.num_atoms() - 1),
            include_dative_bonds: false,
            ..Default::default()
        };

        assert_eq!(
            molecule.to_smiles_with_params(&params).unwrap(),
            "c12c(cccc1)[C@H]1[C@H](C(=O)NC[C@]34C[C@H]5C[C@H](C[C@H](C5)C3)C4)C[C@@H]2c2ccccc21"
        );
    }

    #[test]
    fn writer_noncanonical_explicit_bonds_preserves_rdkit_neighboring_stereo_bond_order() {
        let molecule = Molecule::from_smiles(
            "C=C1/C(C[C@@H](O)CC1)=C\\C=C2[C@@]3([H])[C@@](CCC\\2)(C)[C@]([C@H](C)/C=C/[C@H](C)C(C)C)([H])CC3",
        )
        .unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            all_bonds_explicit: true,
            include_dative_bonds: false,
            ..Default::default()
        };

        assert_eq!(
            molecule.to_smiles_with_params(&params).unwrap(),
            "C=C1/C(=C\\C=C2\\C3-C(-C)(-C-C-C-2)-C(-C(-C)/C=C/C(-C)-C(-C)-C)-C-C-3)-C-C(-O)-C-C-1"
        );
    }

    #[test]
    fn writer_nonisomeric_explicit_bonds_clears_invalid_cx_ring_double_bond_stereo_like_rdkit_row_87()
     {
        let molecule = Molecule::from_smiles("O=C1OCC(=C1c1ccccc1)c1ccccc1 |c:4|").unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            all_bonds_explicit: true,
            include_dative_bonds: false,
            ..Default::default()
        };

        assert_eq!(
            molecule.to_smiles_with_params(&params).unwrap(),
            "O=C1-O-C-C(-c2:c:c:c:c:c:2)=C-1-c1:c:c:c:c:c:1"
        );
    }

    #[test]
    fn writer_preserves_non_ring_aromatic_bridge_as_explicit_single_like_rdkit_biphenyl() {
        let molecule = Molecule::from_smiles("c1ccccc1c1ccccc1").unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: true,
            clean_stereo: false,
            include_dative_bonds: false,
            ..Default::default()
        };

        assert_eq!(
            molecule.to_smiles_with_params(&params).unwrap(),
            "c1ccc(-c2ccccc2)cc1"
        );
    }

    #[test]
    fn writer_emits_single_bridge_between_aromatic_rings_like_rdkit_row_90() {
        let molecule =
            Molecule::from_smiles("CCOC(=O)c1c(N/C=C\\2/C(=NC(=NC2=O)S)O)scc1c1ccc(cc1)Br")
                .unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            include_dative_bonds: false,
            ..Default::default()
        };

        assert_eq!(
            molecule.to_smiles_with_params(&params).unwrap(),
            "CCOC(=O)c1c(NC=C2C(O)=NC(S)=NC2=O)scc1-c1ccc(Br)cc1"
        );
    }

    #[test]
    fn writer_does_not_emit_cleared_tetrahedral_tag_like_rdkit_row_92() {
        let molecule = Molecule::from_smiles(
            "O/C1=C/C=C/C=C1/CN3CCN(CC=2C=CC=CC=2O)[C@]3([H])C=4/C=C(/OC)C(=CC=4)OC",
        )
        .unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: true,
            canonical: false,
            clean_stereo: false,
            include_dative_bonds: false,
            ..Default::default()
        };

        assert_eq!(
            molecule.to_smiles_with_params(&params).unwrap(),
            "Oc1ccccc1CN1CCN(Cc2ccccc2O)C1c1cc(OC)c(OC)cc1"
        );
    }

    #[test]
    fn writer_preserves_ring_special_case_stereo_like_rdkit_row_83() {
        let molecule = Molecule::from_smiles(
            "O=C(NC[C@]12C[C@H]3C[C@H](C[C@H](C3)C1)C2)[C@@H]1C[C@H]2c3ccccc3[C@@H]1c1ccccc12",
        )
        .unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: true,
            canonical: false,
            clean_stereo: false,
            include_dative_bonds: false,
            ..Default::default()
        };

        assert_eq!(
            molecule.to_smiles_with_params(&params).unwrap(),
            "O=C(NC[C@]12C[C@H]3C[C@H](C[C@H](C3)C1)C2)[C@@H]1C[C@H]2c3ccccc3[C@@H]1c1ccccc12"
        );
    }

    #[test]
    fn writer_preserves_ring_closure_double_bond_directions_like_rdkit_row_106() {
        let molecule =
            Molecule::from_smiles("O=C(N(C(S/1)=S)CCC(O)=O)C1=C\\C2=CC=C(C3=CC=C(C=C3)Cl)O2")
                .unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            all_bonds_explicit: true,
            include_dative_bonds: false,
            ..Default::default()
        };

        assert_eq!(
            molecule.to_smiles_with_params(&params).unwrap(),
            "O=C1-N(-C-C-C(-O)=O)-C(=S)-S/C-1=C\\c1:c:c:c(-c2:c:c:c(-Cl):c:c:2):o:1"
        );
    }

    #[test]
    fn writer_preserves_bridgehead_tetrahedral_stereo_like_rdkit_row_110() {
        let molecule = Molecule::from_smiles(
            "[H]Cl.NC[C@@H]1O[C@@H](CC2=C(O)C(O)=CC=C12)C34C[C@H](C5)C[C@H](C[C@H]5C4)C3",
        )
        .unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: true,
            canonical: false,
            clean_stereo: false,
            include_dative_bonds: false,
            ..Default::default()
        };

        assert_eq!(
            molecule.to_smiles_with_params(&params).unwrap(),
            "Cl.NC[C@@H]1O[C@H](C23C[C@@H]4C[C@@H](C[C@@H](C4)C2)C3)Cc2c(O)c(O)ccc21"
        );
    }

    #[test]
    #[ignore = "debug helper for RDKit row 121 traversal parity"]
    fn debug_row_121_noncanonical_traversal() {
        let mut molecule = Molecule::from_smiles(
            "O=C(O[Na])CC1=C(C(C(O[Na])=O)=C(C)C2=CC3=[N]4C(C(C=O)=C3CC)=CC5=C(C=C)C(C)=C6[N-]75)[N-]2[Cu+2]47[N](C8=C6)=C1C(C8C)CCC(O[Na])=O",
        )
        .unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            include_dative_bonds: false,
            ..Default::default()
        };
        let mut working_params = params.clone();
        let _saved_atom_maps =
            prepare_plain_smiles_molecule(&mut molecule, &working_params).unwrap();
        working_params.do_kekule = false;
        let plan = collect_fragment_write_plans(&molecule, &working_params)
            .unwrap()
            .into_iter()
            .next()
            .unwrap();
        let ranks = rank_fragment_atoms_for_smiles(
            &molecule,
            &plan,
            &working_params,
            SmilesOutputMode::PlainSmiles,
        )
        .unwrap();
        molecule = kekulize_for_smiles(&molecule).unwrap();
        let start_atom = choose_fragment_start_atom(&plan, &ranks, &working_params).unwrap();
        let traversal = canonicalize_fragment_stack(
            &molecule,
            &plan,
            start_atom,
            &ranks,
            &working_params,
            SmilesWriteOverrides::default(),
        )
        .unwrap();
        let (ring_info, atom_ring_closures) =
            debug_atom_ring_closures_for_writer(&molecule, &plan, start_atom, &ranks, None)
                .unwrap();
        let mut rank_by_atom = vec![usize::MAX; molecule.num_atoms()];
        for (idx, atom) in plan.atoms.iter().enumerate() {
            rank_by_atom[atom.index()] = ranks[idx];
        }
        for atom_idx in [0usize, 1, 4, 5, 6, 14, 17, 35, 33, 26] {
            let neighbors = molecule
                .topology_block()
                .adjacency
                .neighbors_of(atom_idx)
                .iter()
                .map(|nbr| (nbr.atom_index, nbr.bond.index()))
                .collect::<Vec<_>>();
            eprintln!("neighbors[{atom_idx}]={neighbors:?}");
        }
        for bond_idx in [17usize, 32, 33, 35, 48, 49, 50, 51, 52, 53, 54, 55] {
            eprintln!(
                "bond_rings[{bond_idx}]={}",
                ring_info.num_bond_rings(BondId::new(bond_idx))
            );
        }
        for atom_idx in [17usize, 18, 22, 16, 32, 33, 34, 35, 36, 37, 38, 39, 41] {
            let closures = atom_ring_closures[atom_idx]
                .iter()
                .map(|bond| bond.index())
                .collect::<Vec<_>>();
            eprintln!("atom_ring_closures[{atom_idx}]={closures:?}");
        }
        for atom_idx in [17usize, 35, 33, 26, 32, 36, 39, 40] {
            let atom = AtomId::new(atom_idx);
            let mut seen_from_here = vec![false; molecule.num_atoms()];
            seen_from_here[atom_idx] = true;
            for bond in &atom_ring_closures[atom_idx] {
                let bond = &molecule.bonds()[bond.index()];
                let other = if bond.begin() == atom {
                    bond.end()
                } else {
                    bond.begin()
                };
                seen_from_here[other.index()] = true;
            }
            let children = molecule
                .topology_block()
                .adjacency
                .neighbors_of(atom_idx)
                .iter()
                .map(|nbr| (BondId::new(nbr.bond.index()), AtomId::new(nbr.atom_index)))
                .filter(|(bond, other)| !seen_from_here[other.index()] && Some(*bond) != None)
                .map(|(bond, other)| {
                    let mut rank = rank_by_atom[other.index()] as i64;
                    let bond_order_rank = match molecule.bonds()[bond.index()].order() {
                        BondOrder::Null | BondOrder::Unspecified => 0,
                        BondOrder::Single => 1,
                        BondOrder::Double => 2,
                        BondOrder::Triple => 3,
                        BondOrder::Quadruple => 4,
                        BondOrder::Quintuple => 5,
                        BondOrder::Hextuple => 6,
                        BondOrder::OneAndHalf => 7,
                        BondOrder::TwoAndHalf => 8,
                        BondOrder::ThreeAndHalf => 9,
                        BondOrder::FourAndHalf => 10,
                        BondOrder::FiveAndHalf => 11,
                        BondOrder::Aromatic => 12,
                        BondOrder::Ionic => 13,
                        BondOrder::Hydrogen => 14,
                        BondOrder::ThreeCenter => 15,
                        BondOrder::DativeOne => 16,
                        BondOrder::Dative => 17,
                        BondOrder::DativeLeft => 18,
                        BondOrder::DativeRight => 19,
                        BondOrder::Other => 20,
                        BondOrder::Zero => 21,
                    };
                    if ring_info.num_bond_rings(bond) > 0 {
                        rank += (CANON_MAX_BONDTYPE - bond_order_rank)
                            * CANON_MAX_NATOMS
                            * CANON_MAX_NATOMS;
                    }
                    (rank, bond.index(), other.index())
                })
                .collect::<Vec<_>>();
            eprintln!("children_pre_sort[{atom_idx}]={children:?}");
        }
        eprintln!("start_atom={}", start_atom.index());
        eprintln!("stack={:#?}", traversal.stack);
    }

    #[test]
    #[ignore = "debug helper for RDKit row 142 canonical kekule traversal parity"]
    fn debug_row_142_canonical_kekule_root_last_traversal() {
        let input = "[C:12]12([CH:62]([CH3:65])[c:61]3[cH:64][cH:67][cH:68][cH:66][cH:63]3)[CH:20]4[c:30]5[c:40]6[c:49]7[c:57]8[c:60]([c:59]9[c:55]([c:47]([c:44]([c:52]9[c:51]([c:43]%10[c:35]%11[c:25]%12[c:19]%13%14)[c:53]8[c:45]%11[c:39]6[c:29]4%13)[c:34]([c:24]%15[c:15]%16[c:7]%17[c:3]%18%19)[c:33]%10[c:23]%16[c:16]%12[c:8]%18[c:11]%14[c:5]1%20)[c:37]([c:36]%21[c:26]%22[c:18]%23[c:10]%24[c:13]%25[c:6]%26%27)[c:27]%15[c:17]%22[c:9]%17[c:4]%24[c:1]%19[c:2]%20%26)[c:54]([c:46]%21[c:38]%28[c:28]%23[c:21]%25%29)[c:56]%30[c:48]%28[c:41]%31[c:31]%29[c:22]%32[c:14]2%27)[c:58]%30[c:50]7[c:42]%31[c:32]5%32";
        let mut molecule = Molecule::from_smiles(input).unwrap();
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            do_kekule: true,
            canonical: true,
            clean_stereo: false,
            include_dative_bonds: false,
            ignore_atom_map_numbers: true,
            rooted_at_atom: Some(molecule.num_atoms() - 1),
            ..Default::default()
        };
        let mut working_params = params.clone();
        let _saved_atom_maps =
            prepare_plain_smiles_molecule(&mut molecule, &working_params).unwrap();
        let plan = collect_fragment_write_plans(&molecule, &working_params)
            .unwrap()
            .into_iter()
            .next()
            .unwrap();
        let ranks = rank_fragment_atoms_for_smiles(
            &molecule,
            &plan,
            &working_params,
            SmilesOutputMode::PlainSmiles,
        )
        .unwrap();
        molecule = kekulize_for_smiles(&molecule).unwrap();
        working_params.do_kekule = false;
        let start_atom = choose_fragment_start_atom(&plan, &ranks, &working_params).unwrap();
        let traversal = canonicalize_fragment_stack(
            &molecule,
            &plan,
            start_atom,
            &ranks,
            &working_params,
            SmilesWriteOverrides::default(),
        )
        .unwrap();
        let cached_rings = molecule.derived_cache().rings.clone();
        let fast_rings = crate::fast_find_rings(&molecule).unwrap();
        let (_, atom_ring_closures) =
            debug_atom_ring_closures_for_writer(&molecule, &plan, start_atom, &ranks, None)
                .unwrap();
        eprintln!("start_atom={}", start_atom.index());
        eprintln!("ranks={ranks:?}");
        if let Some(rings) = &cached_rings {
            eprintln!(
                "cached find_type={:?} b11={} b72={}",
                rings.find_type(),
                rings.num_bond_rings(BondId::new(11)),
                rings.num_bond_rings(BondId::new(72))
            );
        } else {
            eprintln!("cached find_type=None");
        }
        eprintln!(
            "fast find_type={:?} b11={} b72={}",
            fast_rings.find_type(),
            fast_rings.num_bond_rings(BondId::new(11)),
            fast_rings.num_bond_rings(BondId::new(72))
        );
        eprintln!(
            "bond orders b11={:?} b72={:?}",
            molecule.bonds()[11].order(),
            molecule.bonds()[72].order()
        );
        eprintln!("atom11 closures={:?}", atom_ring_closures[11]);
        eprintln!("atom27 closures={:?}", atom_ring_closures[27]);
        eprint!("stack");
        for item in &traversal.stack {
            match item {
                MolStackElem::Atom(atom) => eprint!(" A{}", atom.index()),
                MolStackElem::Bond(bond, left) => eprint!(" B{}@L{}", bond.index(), left.index()),
                MolStackElem::Ring { bond: _, ring_idx } => eprint!(" R{ring_idx}"),
                MolStackElem::BranchOpen => eprint!(" ("),
                MolStackElem::BranchClose => eprint!(" )"),
            }
        }
        eprintln!();
        eprintln!(
            "output={}",
            molecule.to_smiles_with_params(&params).unwrap()
        );
    }
}
