// RDKit marker convention defined in dev/source_reproduction_protocol.md.

mod cx;
mod direction;
mod stereo;

pub(crate) use self::cx::get_cx_extensions;
pub(crate) use self::stereo::{assign_stereochemistry_on_working_copy, update_property_cache_for_smiles};
use self::{cx::*, direction::*, stereo::*};

use crate::{AtomId, Bond, BondDirection, BondId, BondOrder, BondStereo, ChiralTag, Molecule, ValenceError};
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
    #[error("molecule invariant failed while preparing SMILES output: {source}")]
    MoleculeInvariant {
        #[from]
        source: crate::InvariantError,
    },
    #[error(transparent)]
    UnsupportedFeature(#[from] crate::UnsupportedFeatureError),
    #[error(transparent)]
    Fragment(#[from] crate::fragment::FragmentError),
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
    #[error("rooted atom index {atom} requires a single-fragment molecule when bonds_to_use is omitted")]
    RootedAtomRequiresSingleFragment { atom: usize },
    #[error("atom symbol override vector has length {len}, expected at least {expected}")]
    AtomSymbolsTooShort { len: usize, expected: usize },
    #[error("bond symbol override vector has length {len}, expected at least {expected}")]
    BondSymbolsTooShort { len: usize, expected: usize },
    #[error("invalid non-tetrahedral chiral permutation {permutation} for {chiral_tag:?}; max allowed is {limit}")]
    InvalidChiralPermutation {
        chiral_tag: ChiralTag,
        permutation: u32,
        limit: u32,
    },
    #[error("invalid ring stereochemistry state on atom {atom}: {requirement}")]
    InvalidRingStereoState { atom: usize, requirement: &'static str },
    #[error("internal SMILES writer invariant violated in {stage}: {message}")]
    InvariantViolation { stage: &'static str, message: &'static str },
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
pub(crate) enum MolStackElem {
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

pub fn mol_to_smiles(molecule: &Molecule, params: &SmilesWriteParams) -> Result<String, SmilesWriteError> {
    mol_to_smiles_with_mode(molecule, params, SmilesOutputMode::PlainSmiles)
}

pub(crate) fn mol_to_smiles_with_atom_output_order(
    molecule: &Molecule,
    params: &SmilesWriteParams,
) -> Result<(String, Vec<AtomId>), SmilesWriteError> {
    let mut context = SmilesWriteContext::default();
    let smiles = mol_to_smiles_with_mode_and_context(molecule, params, SmilesOutputMode::PlainSmiles, &mut context)?;
    Ok((smiles, context.atom_output_order))
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
    let mut context = SmilesWriteContext::default();
    mol_to_smiles_with_mode_and_context(molecule, params, mode, &mut context)
}

fn mol_to_smiles_with_mode_and_context(
    molecule: &Molecule,
    params: &SmilesWriteParams,
    mode: SmilesOutputMode,
    context: &mut SmilesWriteContext,
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

    let mut fragment_results = Vec::new();
    let mut working_params = params.clone();

    let saved_atom_maps = match mode {
        SmilesOutputMode::PlainSmiles => prepare_plain_smiles_molecule(&mut molecule, &working_params)?,
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
            context,
        )?);
        if working_params.canonical && saved_atom_maps.is_some() {
            let _ = stash_and_clear_atom_maps_for_smiles(&mut molecule, &working_params);
        }
    }
    if working_params.canonical {
        restore_atom_maps_after_canonical_smiles(&mut molecule, saved_atom_maps.as_deref());
    }

    let mut result = assemble_fragment_smiles(fragment_results, &working_params, context)?;
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
    let saved = topology.atoms.iter().map(|atom| atom.atom_map()).collect::<Vec<_>>();
    for atom in &mut topology.atoms {
        atom.set_atom_map(None);
    }
    Some(saved)
}

fn restore_atom_maps_after_canonical_smiles(molecule: &mut Molecule, saved_atom_maps: Option<&[Option<u32>]>) {
    let Some(saved_atom_maps) = saved_atom_maps else {
        return;
    };
    let topology = molecule.topology_block_mut();
    for (atom, atom_map) in topology.atoms.iter_mut().zip(saved_atom_maps.iter().copied()) {
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
    fragment_smiles_construct(molecule, plan, start_atom, &ranks, params, overrides, context)
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
    let traversal = canonicalize_fragment_stack(molecule, plan, start_atom, ranks, params, overrides)?;
    canonicalize_double_bond_directions_for_writer(
        molecule,
        &traversal.stack,
        &traversal.traversal_ring_closure_bonds,
    )?;
    context
        .chiral_tag_overrides
        .extend(traversal.chiral_tag_overrides.iter().map(|(atom, tag)| (*atom, *tag)));
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
    // RDKit✔️✔️:   auto mols =
    // RDKit✔️✔️:       MolOps::getMolFrags(mol, false, nullptr, &fragsMolAtomMapping, false);
    // RDKit✔️✔️:   for (unsigned fragIdx = 0; fragIdx < mols.size(); fragIdx++) {
    // RDKit✔️✔️:     ROMol *tmol = mols[fragIdx].get();
    let fragment_atom_indices = plan.atoms.iter().map(|atom| atom.index()).collect::<Vec<_>>();
    let fragment = crate::fragment::build_fragment_molecule(molecule, &fragment_atom_indices, false)?;
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
        &fragment,
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
    Ok(ranks)
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

pub(crate) fn canonicalize_fragment_stack_for_smarts(
    molecule: &Molecule,
    atoms: &[AtomId],
    bonds: &[BondId],
    start_atom: AtomId,
    ranks: &[usize],
    params: &SmilesWriteParams,
) -> Result<Vec<MolStackElem>, SmilesWriteError> {
    let plan = FragmentWritePlan {
        atoms: atoms.to_vec(),
        bonds: bonds.to_vec(),
        rooted_at_atom: Some(start_atom),
    };
    let mut traversal_params = params.clone();
    traversal_params.do_random = false;
    Ok(canonicalize_fragment_stack(
        molecule,
        &plan,
        start_atom,
        ranks,
        &traversal_params,
        SmilesWriteOverrides::default(),
    )?
    .stack)
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
                    result
                        .smiles
                        .push_str(&build_bond_smiles(molecule, bond, atom_to_left, params)?);
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
    validate_fragment_api_inputs(molecule, params, atoms_to_use, bonds_to_use, atom_symbols, bond_symbols)?;
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

    let mut plans = collect_fragment_api_write_plans(&molecule, &working_params, atoms_to_use, bonds_to_use)?;
    if working_params.canonical {
        restore_atom_maps_after_canonical_smiles(&mut molecule, saved_atom_maps.as_deref());
        plans.sort_by_key(|plan| plan.atoms.iter().map(|atom| atom.index()).min().unwrap_or(usize::MAX));
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
    validate_fragment_api_inputs(molecule, params, atoms_to_use, bonds_to_use, atom_symbols, bond_symbols)?;
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

    let mut plans = collect_fragment_api_write_plans(&molecule, &working_params, atoms_to_use, bonds_to_use)?;
    if working_params.canonical {
        restore_atom_maps_after_canonical_smiles(&mut molecule, saved_atom_maps.as_deref());
        plans.sort_by_key(|plan| plan.atoms.iter().map(|atom| atom.index()).min().unwrap_or(usize::MAX));
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
            .filter(|bond| atom_set.contains(&bond.begin().index()) && atom_set.contains(&bond.end().index()))
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
    get_atom_smiles_impl(molecule, AtomId::new(atom), params, None, false, None, false)
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
    let symbol: &str =
        if !params.do_kekule && atom.is_aromatic() && raw_symbol.as_bytes().first().is_some_and(u8::is_ascii_uppercase)
        {
            let should_lower = matches!(atom.atomic_number(), 5 | 6 | 7 | 8 | 14 | 15 | 16 | 33 | 34 | 52);
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
        let other_id = bond_other_atom(bond, AtomId::new(atom_to_left)).ok_or(SmilesWriteError::BondOutOfRange {
            bond: bond.id().index(),
        })?;
        let other = &molecule.atoms()[other_id.index()];
        left.is_aromatic() && other.is_aromatic() && (left.atomic_number() != 0 || other.atomic_number() != 0)
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
            if !matches!(bond.direction(), BondDirection::None | BondDirection::Unknown) {
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
            if !matches!(bond.direction(), BondDirection::None | BondDirection::Unknown) {
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
    molecule
        .derived_cache()
        .valence
        .as_ref()
        .map(|valence| valence.explicit_valence[atom_id.index()] + valence.implicit_hydrogens[atom_id.index()])
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
    Ok(matches!(_atomic_number, 0 | 5 | 6 | 7 | 8 | 9 | 15 | 16 | 17 | 35 | 53))
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

    let digit = match (1..).find(|candidate| !context.ring_closure_digits.values().any(|digit| digit == candidate)) {
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
    // Preserve the writer's existing out-of-table representation while the
    // source-backed valid-symbol lookup remains canonical and shared.
    Ok(crate::rdkit_element_symbol(atomic_number).unwrap_or("?"))
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
            context.atom_output_order.extend(fragment.atom_ordering.iter().copied());
            context.bond_output_order.extend(fragment.bond_ordering.iter().copied());
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
        context.atom_output_order.extend(fragment.atom_ordering.iter().copied());
        context.bond_output_order.extend(fragment.bond_ordering.iter().copied());
    }
    Ok(fragment_results
        .into_iter()
        .map(|fragment| fragment.smiles)
        .collect::<Vec<_>>()
        .join("."))
}

fn validate_rooted_atom(molecule: &Molecule, params: &SmilesWriteParams) -> Result<(), SmilesWriteError> {
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

fn invariant_stage_error<T>(stage: SmilesPlanStage, message: &'static str) -> Result<T, SmilesWriteError> {
    Err(SmilesWriteError::InvariantViolation {
        stage: stage.as_str(),
        message,
    })
}

#[cfg(test)]
mod tests;
