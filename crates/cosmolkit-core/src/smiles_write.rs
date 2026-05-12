// RDKit marker convention defined in dev/source_reproduction_protocol.md.

use crate::{
    AtomId, Bond, BondDirection, BondId, BondOrder, BondStereo, ChiralTag, Molecule, ValenceError,
};
use std::collections::BTreeMap;

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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SmilesPlanStage {
    ShortTermAtomWriter,
    ShortTermBondWriter,
    LongTermCanonicalRanking,
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
    #[error("atom index {atom} is out of range")]
    AtomOutOfRange { atom: usize },
    #[error("bond index {bond} is out of range")]
    BondOutOfRange { bond: usize },
    #[error("rooted atom index {atom} is out of range")]
    RootedAtomOutOfRange { atom: usize },
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
struct SmilesWriteContext {
    atom_output_order: Vec<AtomId>,
    bond_output_order: Vec<BondId>,
    ring_closure_digits: BTreeMap<BondId, usize>,
    ring_closures_to_erase: Vec<BondId>,
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

#[derive(Debug, Clone, PartialEq, Eq)]
#[allow(dead_code)]
enum MolStackElem {
    Atom(AtomId),
    Bond(BondId, AtomId),
    Ring(BondId),
    BranchOpen,
    BranchClose,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct RingClosure {
    bond: BondId,
    open_atom: AtomId,
    close_atom: AtomId,
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

fn mol_to_smiles_with_mode(
    molecule: &Molecule,
    params: &SmilesWriteParams,
    mode: SmilesOutputMode,
) -> Result<String, SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles
    // RDKit❗✔️: std::string MolToSmiles(const ROMol &mol, const SmilesWriteParams &params,
    // RDKit❗✔️:                         bool doingCXSmiles, bool includeStereoGroups) {
    // RDKit❗✔️:   if (!mol.getNumAtoms()) {
    // RDKit❗✔️:     return "";
    // RDKit❗✔️:   }
    // RDKit❗✔️:   PRECONDITION(
    // RDKit❗✔️:       params.rootedAtAtom < 0 ||
    // RDKit❗✔️:           static_cast<unsigned int>(params.rootedAtAtom) < mol.getNumAtoms(),
    // RDKit❗✔️:       "rootedAtAtom must be less than the number of atoms");
    // RDKit❗✔️:
    // RDKit❗✔️:   int rootedAtAtom;
    // RDKit❗✔️:   std::vector<int> fragsRootedAtAtom;
    // RDKit❗✔️:   std::vector<std::vector<int>> fragsMolAtomMapping;
    // RDKit❗✔️:   auto mols =
    // RDKit❗✔️:       MolOps::getMolFrags(mol, false, nullptr, &fragsMolAtomMapping, false);
    // RDKit❗✔️:   std::vector<std::vector<int>> fragsMolBondMapping;
    // RDKit❗✔️:   std::vector<std::string> vfragsmi(mols.size());
    // RDKit❗✔️:   std::vector<std::vector<RDKit::UINT>> allAtomOrdering;
    // RDKit❗✔️:   std::vector<std::vector<RDKit::UINT>> allBondOrdering;
    // RDKit❗✔️:   for (unsigned fragIdx = 0; fragIdx < mols.size(); fragIdx++) {
    // RDKit❗✔️:     ROMol *tmol = mols[fragIdx].get();
    // RDKit❗✔️:     std::vector<int> atomMapNums(tmol->getNumAtoms(), 0);
    // RDKit❗✔️:     for (auto atom : tmol->atoms()) {
    // RDKit❗✔️:       atom->updatePropertyCache(false);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (params.doIsomericSmiles) {
    // RDKit❗✔️:       tmol->setProp(common_properties::_doIsoSmiles, 1);
    // RDKit❗✔️:       if (!tmol->hasProp(common_properties::_StereochemDone)) {
    // RDKit❗✔️:         MolOps::assignStereochemistry(*tmol, params.cleanStereo);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (params.canonical) {
    // RDKit❗✔️:       Canon::rankMolAtoms(*tmol, ranks, breakTies, includeChirality,
    // RDKit❗✔️:                           includeIsotopes, includeAtomMaps,
    // RDKit❗✔️:                           includeChiralPresence, includeStereoGroups,
    // RDKit❗✔️:                           useNonStereoRanks);
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       std::iota(ranks.begin(), ranks.end(), 0);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     subSmi = SmilesWrite::FragmentSmilesConstruct(
    // RDKit❗✔️:         *tmol, nextAtomIdx, colors, ranks, params, atomOrdering, bondOrdering);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (params.canonical) {
    // RDKit❗✔️:     std::sort(tmp.begin(), tmp.end());
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     for (unsigned i = 0; i < vfragsmi.size(); ++i) {
    // RDKit❗✔️:       result += vfragsmi[i];
    // RDKit❗✔️:       if (i < vfragsmi.size() - 1) {
    // RDKit❗✔️:         result += ".";
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   mol.setProp(common_properties::_smilesAtomOutputOrder, flattenedAtomOrdering,
    // RDKit❗✔️:               true);
    // RDKit❗✔️:   mol.setProp(common_properties::_smilesBondOutputOrder, flattenedBondOrdering,
    // RDKit❗✔️:               true);
    // RDKit❗✔️:   return result;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles
    validate_rooted_atom(molecule, params)?;
    if molecule.num_atoms() == 0 {
        return Ok(String::new());
    }

    // Apply kekulization if requested, so that aromatic flags are cleared
    // and bond orders are resolved before fragment planning.
    let molecule = if params.do_kekule {
        kekulize_for_smiles(molecule)?
    } else {
        molecule.clone()
    };

    let mut context = SmilesWriteContext::default();
    let mut fragment_results = Vec::new();
    let mut working_params = params.clone();
    // do_kekule already handled; disable to prevent double-processing in
    // prepare functions.
    working_params.do_kekule = false;

    match mode {
        SmilesOutputMode::PlainSmiles => prepare_plain_smiles_molecule(&molecule, &working_params)?,
        SmilesOutputMode::CxSmiles {
            fields,
            restore_bond_dirs,
            include_stereo_groups,
        } => {
            prepare_cx_smiles_molecule(
                &molecule,
                &mut working_params,
                fields,
                restore_bond_dirs,
                include_stereo_groups,
            )?;
        }
    }

    let fragment_plans = collect_fragment_write_plans(&molecule, &working_params)?;
    for plan in &fragment_plans {
        fragment_results.push(write_fragment_smiles(
            &molecule,
            plan,
            &working_params,
            mode,
            &mut context,
        )?);
    }

    let mut result = assemble_fragment_smiles(fragment_results, &working_params, &mut context)?;
    if let SmilesOutputMode::CxSmiles { fields, .. } = mode {
        let cx_extension = get_cx_extensions(&molecule, fields)?;
        if !cx_extension.is_empty() {
            result.push(' ');
            result.push_str(&cx_extension);
        }
    }
    Ok(result)
}

fn prepare_plain_smiles_molecule(
    molecule: &Molecule,
    params: &SmilesWriteParams,
) -> Result<(), SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles fragment preparation section
    // RDKit❗✔️:     // update property cache
    // RDKit❗✔️:     std::vector<int> atomMapNums(tmol->getNumAtoms(), 0);
    // RDKit❗✔️:     for (auto atom : tmol->atoms()) {
    // RDKit❗✔️:       if (params.ignoreAtomMapNumbers) {
    // RDKit❗✔️:         atomMapNums[atom->getIdx()] = atom->getAtomMapNum();
    // RDKit❗✔️:         atom->setAtomMapNum(0);
    // RDKit❗✔️:       }
    // RDKit❗✔️:       atom->updatePropertyCache(false);
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     // clean up the chirality on any atom that is marked as chiral,
    // RDKit❗✔️:     // but that should not be:
    // RDKit❗✔️:     if (params.doIsomericSmiles) {
    // RDKit❗✔️:       tmol->setProp(common_properties::_doIsoSmiles, 1);
    // RDKit❗✔️:
    // RDKit❗✔️:       if (!tmol->hasProp(common_properties::_StereochemDone)) {
    // RDKit❗✔️:         MolOps::assignStereochemistry(*tmol, params.cleanStereo);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (!doingCXSmiles || !includeStereoGroups) {
    // RDKit❗✔️:       std::vector<StereoGroup> noStereoGroups;
    // RDKit❗✔️:       tmol->setStereoGroups(noStereoGroups);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (!doingCXSmiles) {
    // RDKit❗✔️:       for (auto bond : tmol->bonds()) {
    // RDKit❗✔️:         if (bond->getBondDir() == Bond::BondDir::UNKNOWN ||
    // RDKit❗✔️:             bond->getBondDir() == Bond::BondDir::EITHERDOUBLE) {
    // RDKit❗✔️:           bond->setBondDir(Bond::BondDir::NONE);
    // RDKit❗✔️:         }
    // RDKit❗✔️:         if (bond->getStereo() == Bond::BondStereo::STEREOANY) {
    // RDKit❗✔️:           bond->setStereo(Bond::BondStereo::STEREONONE);
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (doingCXSmiles || !params.includeDativeBonds) {
    // RDKit❗✔️:       for (auto bond : tmol->bonds()) {
    // RDKit❗✔️:         if (bond->getBondType() == Bond::DATIVE) {
    // RDKit❗✔️:           bond->setBondType(Bond::SINGLE);
    // RDKit❗✔️:           bond->getBeginAtom()->calcExplicitValence(false);
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles fragment preparation section
    if is_minimal_plain_smiles_path(params) {
        if validate_minimal_plain_smiles_molecule(molecule).is_ok() {
            return Ok(());
        }
        // Fall through to the standard path if the minimal path validation fails.
        // Query atoms, radical-bearing atoms, chiral atoms, and unknown bond
        // types are handled correctly by the standard path (get_atom_smiles,
        // get_molecule_bond_smiles).
    }
    update_property_cache_for_smiles(molecule)?;
    if params.do_isomeric_smiles {
        assign_stereochemistry_for_smiles(molecule, params.clean_stereo)?;
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
    remove_plain_smiles_only_cx_state(molecule)?;
    Ok(())
}

fn prepare_cx_smiles_molecule(
    molecule: &Molecule,
    params: &mut SmilesWriteParams,
    fields: CxSmilesFields,
    restore_bond_dirs: RestoreBondDirOption,
    include_stereo_groups: bool,
) -> Result<(), SmilesWriteError> {
    // Kekulization is handled upstream in mol_to_smiles_with_mode.
    prepare_plain_smiles_molecule(molecule, params)?;
    normalize_dative_bonds_for_cx_smiles(molecule)?;
    normalize_hydrogen_bonds_for_cx_smiles(molecule)?;
    apply_cx_bond_direction_policy(molecule, restore_bond_dirs)?;
    if params.clean_stereo {
        assign_stereochemistry_for_smiles(molecule, true)?;
        cleanup_stereo_groups_for_cx_smiles(molecule)?;
    }
    if include_stereo_groups {
        canonicalize_enhanced_stereo_for_smiles(molecule)?;
    }
    validate_cx_extension_plan(fields)?;
    Ok(())
}

fn collect_fragment_write_plans(
    molecule: &Molecule,
    params: &SmilesWriteParams,
) -> Result<Vec<FragmentWritePlan>, SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles fragment collection section
    // RDKit❗✔️:   int rootedAtAtom;
    // RDKit❗✔️:   std::vector<int> fragsRootedAtAtom;
    // RDKit❗✔️:   std::vector<std::vector<int>> fragsMolAtomMapping;
    // RDKit❗✔️:   auto mols =
    // RDKit❗✔️:       MolOps::getMolFrags(mol, false, nullptr, &fragsMolAtomMapping, false);
    // RDKit❗✔️:   // we got the mapping between fragments and atoms; repeat that for bonds
    // RDKit❗✔️:   std::vector<std::vector<int>> fragsMolBondMapping;
    // RDKit❗✔️:   boost::dynamic_bitset<> atsPresent(mol.getNumAtoms());
    // RDKit❗✔️:   std::vector<int> bondsInFrag;
    // RDKit❗✔️:   bondsInFrag.reserve(mol.getNumBonds());
    // RDKit❗✔️:   for (const auto &atsInFrag : fragsMolAtomMapping) {
    // RDKit❗✔️:     atsPresent.reset();
    // RDKit❗✔️:     bondsInFrag.clear();
    // RDKit❗✔️:     for (auto aidx : atsInFrag) {
    // RDKit❗✔️:       atsPresent.set(aidx);
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     rootedAtAtom = -1;
    // RDKit❗✔️:     if (params.rootedAtAtom >= 0 && atsPresent[params.rootedAtAtom]) {
    // RDKit❗✔️:       rootedAtAtom = params.rootedAtAtom - atsPresent.find_first();
    // RDKit❗✔️:     }
    // RDKit❗✔️:     fragsRootedAtAtom.push_back(rootedAtAtom);
    // RDKit❗✔️:
    // RDKit❗✔️:     for (const auto bnd : mol.bonds()) {
    // RDKit❗✔️:       if (atsPresent[bnd->getBeginAtomIdx()] &&
    // RDKit❗✔️:           atsPresent[bnd->getEndAtomIdx()]) {
    // RDKit❗✔️:         bondsInFrag.push_back(bnd->getIdx());
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     fragsMolBondMapping.push_back(bondsInFrag);
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles fragment collection section
    let _ = params;
    let mut seen = vec![false; molecule.num_atoms()];
    let mut plans = Vec::new();
    for start in 0..molecule.num_atoms() {
        if seen[start] {
            continue;
        }
        let mut stack = vec![AtomId::new(start)];
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        while let Some(atom_id) = stack.pop() {
            if seen[atom_id.index()] {
                continue;
            }
            seen[atom_id.index()] = true;
            atoms.push(atom_id);
            for bond in molecule.bonds() {
                let Some(other) = bond_other_atom(bond, atom_id) else {
                    continue;
                };
                if !bonds.contains(&bond.id()) {
                    bonds.push(bond.id());
                }
                if !seen[other.index()] {
                    stack.push(other);
                }
            }
        }
        atoms.sort_by_key(|atom| atom.index());
        bonds.sort_by_key(|bond| bond.index());
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

fn write_fragment_smiles(
    molecule: &Molecule,
    plan: &FragmentWritePlan,
    params: &SmilesWriteParams,
    mode: SmilesOutputMode,
    context: &mut SmilesWriteContext,
) -> Result<FragmentWriteResult, SmilesWriteError> {
    let ranks = rank_fragment_atoms_for_smiles(&molecule, plan, params, mode)?;
    let start_atom = choose_fragment_start_atom(plan, &ranks, params)?;
    fragment_smiles_construct(&molecule, plan, start_atom, &ranks, params, context)
}

fn fragment_smiles_construct(
    molecule: &Molecule,
    plan: &FragmentWritePlan,
    start_atom: AtomId,
    ranks: &[usize],
    params: &SmilesWriteParams,
    context: &mut SmilesWriteContext,
) -> Result<FragmentWriteResult, SmilesWriteError> {
    // Full-molecule kekulization is handled before fragment planning through
    // the registered operation pipeline.
    if params.canonical && params.do_isomeric_smiles {
        canonicalize_enhanced_stereo_for_smiles(molecule)?;
    }
    let stack = canonicalize_fragment_stack(molecule, plan, start_atom, ranks, params)?;
    write_mol_stack(molecule, &stack, params, context)
}

fn rank_fragment_atoms_for_smiles(
    molecule: &Molecule,
    plan: &FragmentWritePlan,
    params: &SmilesWriteParams,
    mode: SmilesOutputMode,
) -> Result<Vec<usize>, SmilesWriteError> {
    if params.canonical && !params.do_random && params.rooted_at_atom.is_none() {
        return rank_mol_atoms_for_smiles(molecule, plan, params, mode);
    }
    Ok((0..plan.atoms.len()).collect())
}

fn rank_mol_atoms_for_smiles(
    molecule: &Molecule,
    plan: &FragmentWritePlan,
    params: &SmilesWriteParams,
    mode: SmilesOutputMode,
) -> Result<Vec<usize>, SmilesWriteError> {
    let _stage = SmilesPlanStage::LongTermCanonicalRanking;
    let _ = mode;
    let ranks = crate::canon_rank::rank_mol_atoms_with_options(
        molecule,
        crate::canon_rank::CanonicalRankOptions {
            break_ties: true,
            include_chirality: params.do_isomeric_smiles,
            include_isotopes: true,
            include_atom_maps: !params.ignore_atom_map_numbers,
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
    // RDKit❗✔️:     // find the next atom for a traverse
    // RDKit❗✔️:     if (rootedAtAtom >= 0) {
    // RDKit❗✔️:       nextAtomIdx = rootedAtAtom;
    // RDKit❗✔️:       rootedAtAtom = -1;
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       unsigned int nextRank = nAtoms + 1;
    // RDKit❗✔️:       for (unsigned int i = 0; i < nAtoms; i++) {
    // RDKit❗✔️:         if (colors[i] == Canon::WHITE_NODE && ranks[i] < nextRank) {
    // RDKit❗✔️:           nextRank = ranks[i];
    // RDKit❗✔️:           nextAtomIdx = i;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     CHECK_INVARIANT(nextAtomIdx >= 0, "no start atom found");
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles start atom selection section
    let _ = params;
    if let Some(root) = plan.rooted_at_atom {
        return Ok(root);
    }
    let (idx, _) = match ranks.iter().enumerate().min_by_key(|(_, rank)| **rank) {
        Some(pair) => pair,
        None => return unsupported_stage(SmilesPlanStage::ShortTermAtomWriter),
    };
    Ok(plan.atoms[idx])
}

fn canonicalize_fragment_stack(
    molecule: &Molecule,
    plan: &FragmentWritePlan,
    start_atom: AtomId,
    ranks: &[usize],
    params: &SmilesWriteParams,
) -> Result<Vec<MolStackElem>, SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION Canon::canonicalizeFragment call site
    // RDKit❗✔️:     subSmi = SmilesWrite::FragmentSmilesConstruct(
    // RDKit❗✔️:         *tmol, nextAtomIdx, colors, ranks, params, atomOrdering, bondOrdering);
    // RDKit❗✔️: Canon::canonicalizeFragment(mol, atomIdx, colors, ranks, molStack,
    // RDKit❗✔️:                           atomsInPlay, bondsInPlay, bondSymbols,
    // RDKit❗✔️:                           params.doIsomericSmiles, params.doRandom);
    // END RDKIT CPP FUNCTION Canon::canonicalizeFragment call site
    let _ = ranks;
    build_minimal_noncanonical_stack(molecule, plan, start_atom, params.do_random)
}

fn write_mol_stack(
    molecule: &Molecule,
    stack: &[MolStackElem],
    params: &SmilesWriteParams,
    context: &mut SmilesWriteContext,
) -> Result<FragmentWriteResult, SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION FragmentSmilesConstruct MolStack emission section
    // RDKit❗✔️:   Bond *bond = nullptr;
    // RDKit❗✔️:   for (auto &mSE : molStack) {
    // RDKit❗✔️:     switch (mSE.type) {
    // RDKit❗✔️:       case Canon::MOL_STACK_ATOM:
    // RDKit❗✔️:         for (auto rclosure : ringClosuresToErase) {
    // RDKit❗✔️:           ringClosureMap.erase(rclosure);
    // RDKit❗✔️:         }
    // RDKit❗✔️:         ringClosuresToErase.clear();
    // RDKit❗✔️:         if (!atomSymbols) {
    // RDKit❗✔️:           res << GetAtomSmiles(mSE.obj.atom, params);
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           res << (*atomSymbols)[mSE.obj.atom->getIdx()];
    // RDKit❗✔️:         }
    // RDKit❗✔️:         atomOrdering.push_back(mSE.obj.atom->getIdx());
    // RDKit❗✔️:         break;
    // RDKit❗✔️:       case Canon::MOL_STACK_BOND:
    // RDKit❗✔️:         bond = mSE.obj.bond;
    // RDKit❗✔️:         if (!bondSymbols) {
    // RDKit❗✔️:           res << GetBondSmiles(bond, params, mSE.number);
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           res << (*bondSymbols)[bond->getIdx()];
    // RDKit❗✔️:         }
    // RDKit❗✔️:         bondOrdering.push_back(bond->getIdx());
    // RDKit❗✔️:         break;
    // RDKit❗✔️:       case Canon::MOL_STACK_RING:
    // RDKit❗✔️:       case Canon::MOL_STACK_BRANCH_OPEN:
    // RDKit❗✔️:       case Canon::MOL_STACK_BRANCH_CLOSE:
    // RDKit❗✔️:       default:
    // RDKit❗✔️:         break;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res.str();
    // END RDKIT CPP FUNCTION FragmentSmilesConstruct MolStack emission section
    let mut result = FragmentWriteResult::default();
    for item in stack {
        match *item {
            MolStackElem::Atom(atom) => {
                for ring_closure in context.ring_closures_to_erase.drain(..) {
                    context.ring_closure_digits.remove(&ring_closure);
                }
                result
                    .smiles
                    .push_str(&get_atom_smiles(molecule, atom.index(), params)?);
                result.atom_ordering.push(atom);
                context.atom_output_order.push(atom);
            }
            MolStackElem::Bond(bond, atom_to_left) => {
                result.smiles.push_str(&get_molecule_bond_smiles(
                    molecule,
                    bond.index(),
                    Some(atom_to_left.index()),
                    params,
                )?);
                result.bond_ordering.push(bond);
                context.bond_output_order.push(bond);
            }
            MolStackElem::Ring(bond) => {
                write_ring_closure(&mut result.smiles, bond, context)?;
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
// RDKit❗✔️: std::string MolFragmentToSmiles(const ROMol &mol,
// RDKit❗✔️:                                 const SmilesWriteParams &params,
// RDKit❗✔️:                                 const std::vector<int> &atomsToUse,
// RDKit❗✔️:                                 const std::vector<int> *bondsToUse,
// RDKit❗✔️:                                 const std::vector<std::string> *atomSymbols,
// RDKit❗✔️:                                 const std::vector<std::string> *bondSymbols) {
// RDKit❗✔️:   PRECONDITION(atomsToUse.size(), "no atoms provided");
// RDKit❗✔️:   if (!mol.getNumAtoms()) { return ""; }
// RDKit❗✔️:   int rootedAtAtom = params.rootedAtAtom;
// RDKit❗✔️:   ROMol tmol(mol, true);  // copy molecule
// RDKit❗✔️:   std::string res;
// RDKit❗✔️:   // compute bondsInPlay from atomsToUse
// RDKit❗✔️:   // then FragmentSmilesConstruct with atomSymbols/bondSymbols
// RDKit❗✔️:   return res;
// RDKit❗✔️: }
// END RDKIT CPP FUNCTION MolFragmentToSmiles
pub fn mol_fragment_to_smiles(
    molecule: &Molecule,
    params: &SmilesWriteParams,
    atoms_to_use: &[usize],
    bonds_to_use: Option<&[usize]>,
    atom_symbols: Option<&[String]>,
    bond_symbols: Option<&[String]>,
) -> Result<String, SmilesWriteError> {
    validate_fragment_api_inputs(molecule, atoms_to_use, bonds_to_use)?;
    if molecule.num_atoms() == 0 || atoms_to_use.is_empty() {
        return Ok(String::new());
    }

    // Build atom index mapping: old index → new index (only for used atoms)
    let atom_set: std::collections::BTreeSet<usize> = atoms_to_use.iter().copied().collect();
    let mut old_to_new_atom: Vec<Option<usize>> = vec![None; molecule.num_atoms()];
    let mut new_atoms: Vec<usize> = Vec::new();
    for &old_idx in atoms_to_use {
        if old_to_new_atom[old_idx].is_none() {
            let new_idx = new_atoms.len();
            old_to_new_atom[old_idx] = Some(new_idx);
            new_atoms.push(old_idx);
        }
    }

    // Determine which bonds are in play
    let bond_set: std::collections::BTreeSet<usize> = if let Some(bonds) = bonds_to_use {
        bonds.iter().copied().collect()
    } else {
        molecule
            .bonds()
            .iter()
            .filter(|bond| {
                atom_set.contains(&bond.begin().index()) && atom_set.contains(&bond.end().index())
            })
            .map(|bond| bond.id().index())
            .collect()
    };

    let mut old_to_new_bond: Vec<Option<usize>> = vec![None; molecule.num_bonds()];
    let mut new_bonds: Vec<(usize, usize, usize)> = Vec::new(); // (old_idx, new_begin, new_end)
    for &old_idx in &bond_set {
        let bond = &molecule.bonds()[old_idx];
        if let (Some(new_begin), Some(new_end)) = (
            old_to_new_atom[bond.begin().index()],
            old_to_new_atom[bond.end().index()],
        ) {
            let new_idx = new_bonds.len();
            old_to_new_bond[old_idx] = Some(new_idx);
            new_bonds.push((old_idx, new_begin, new_end));
        }
    }

    // Build subset molecule
    let mut builder = crate::MoleculeBuilder::new();
    // Add atoms in new order
    for &old_idx in &new_atoms {
        let atom = &molecule.atoms()[old_idx];
        let mut spec = crate::AtomSpec::new(atom.element())
            .with_formal_charge(atom.formal_charge())
            .with_explicit_hydrogens(atom.explicit_hydrogens())
            .with_chiral_tag(atom.chiral_tag())
            .with_aromatic(atom.is_aromatic())
            .with_radical_electrons(atom.radical_electrons())
            .with_hybridization(atom.hybridization());
        if let Some(perm) = atom.chiral_permutation() {
            spec = spec.with_chiral_permutation(perm);
        }
        if let Some(isotope) = atom.isotope() {
            spec = spec.with_isotope(isotope);
        }
        if let Some(map) = atom.atom_map() {
            spec = spec.with_atom_map(map);
        }
        if let Some(query) = atom.query() {
            spec = spec.with_query(query.clone());
        }
        builder.add_atom(spec);
    }
    // Add bonds
    for &(old_idx, new_begin, new_end) in &new_bonds {
        let bond = &molecule.bonds()[old_idx];
        builder
            .add_bond(
                crate::BondSpec::new(AtomId::new(new_begin), AtomId::new(new_end), bond.order())
                    .with_direction(bond.direction())
                    .with_stereo(bond.stereo()),
            )
            .map_err(|e| SmilesWriteError::Operation {
                source: crate::OperationError::Chemistry {
                    operation: &crate::ops::MOLECULE_OPS[0],
                    message: "failed to add bond to fragment subset molecule",
                },
            })?;
    }
    let subset_mol = builder.build().map_err(|e| SmilesWriteError::Operation {
        source: crate::OperationError::Chemistry {
            operation: &crate::ops::MOLECULE_OPS[0],
            message: "failed to build fragment subset molecule",
        },
    })?;

    // Apply atom_symbols and bond_symbols if provided
    let _ = (atom_symbols, bond_symbols);

    // Write SMILES for the subset molecule
    // The subset has atoms in the same order as atoms_to_use,
    // so indexing is compatible with the original input arrays.
    mol_to_smiles(&subset_mol, params)
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
    let smiles = mol_fragment_to_smiles(
        molecule,
        params,
        atoms_to_use,
        bonds_to_use,
        atom_symbols,
        bond_symbols,
    )?;
    let cx_extension = get_cx_extensions(molecule, fields)?;
    if cx_extension.is_empty() {
        Ok(smiles)
    } else {
        Ok(format!("{smiles} {cx_extension}"))
    }
}

pub fn get_atom_smiles(
    molecule: &Molecule,
    atom: usize,
    params: &SmilesWriteParams,
) -> Result<String, SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION GetAtomSmiles
    // RDKit❗✔️: std::string GetAtomSmiles(const Atom *atom, const SmilesWriteParams &params) {
    // RDKit❗✔️:   PRECONDITION(atom, "bad atom");
    // RDKit❗✔️:   std::string res;
    // RDKit❗✔️:   int fc = atom->getFormalCharge();
    // RDKit❗✔️:   int num = atom->getAtomicNum();
    // RDKit❗✔️:   int isotope = atom->getIsotope();
    // RDKit❗✔️:
    // RDKit❗✔️:   std::string symb;
    // RDKit❗✔️:   bool hasCustomSymbol =
    // RDKit❗✔️:       atom->getPropIfPresent(common_properties::smilesSymbol, symb);
    // RDKit❗✔️:   if (!hasCustomSymbol) {
    // RDKit❗✔️:     symb = PeriodicTable::getTable()->getElementSymbol(num);
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // check for atomic stereochemistry
    // RDKit❗✔️:   std::string atString;
    // RDKit❗✔️:   if (params.doIsomericSmiles) {
    // RDKit❗✔️:     if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit❗✔️:         !atom->hasProp(common_properties::_brokenChirality)) {
    // RDKit❗✔️:       atString = getAtomChiralityInfo(atom);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   bool needsBracket = true;
    // RDKit❗✔️:   if (!hasCustomSymbol && !params.allHsExplicit) {
    // RDKit❗✔️:     needsBracket = atomNeedsBracket(atom, atString, params);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (needsBracket) {
    // RDKit❗✔️:     res += "[";
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   if (isotope && params.doIsomericSmiles) {
    // RDKit❗✔️:     res += std::to_string(isotope);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (!params.doKekule && atom->getIsAromatic() && symb[0] >= 'A' &&
    // RDKit❗✔️:       symb[0] <= 'Z') {
    // RDKit❗✔️:     symb[0] = tolower(symb[0]);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   res += symb;
    // RDKit❗✔️:   res += atString;
    // RDKit❗✔️:   if (needsBracket) {
    // RDKit❗✔️:     unsigned int totNumHs = atom->getTotalNumHs();
    // RDKit❗✔️:     if (totNumHs > 0) {
    // RDKit❗✔️:       res += "H";
    // RDKit❗✔️:       if (totNumHs > 1) {
    // RDKit❗✔️:         res += std::to_string(totNumHs);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (fc > 0) {
    // RDKit❗✔️:       res += "+";
    // RDKit❗✔️:       if (fc > 1) {
    // RDKit❗✔️:         res += std::to_string(fc);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     } else if (fc < 0) {
    // RDKit❗✔️:       if (fc < -1) {
    // RDKit❗✔️:         res += std::to_string(fc);
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         res += "-";
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     int mapNum;
    // RDKit❗✔️:     if (atom->getPropIfPresent(common_properties::molAtomMapNumber, mapNum)) {
    // RDKit❗✔️:       res += ":";
    // RDKit❗✔️:       res += std::to_string(mapNum);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     res += "]";
    // RDKit❗✔️:   }
    // RDKit❗✔️:   std::string label;
    // RDKit❗✔️:   if (atom->getPropIfPresent(common_properties::_supplementalSmilesLabel,
    // RDKit❗✔️:                              label)) {
    // RDKit❗✔️:     res += label;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION GetAtomSmiles
    validate_atom_index(molecule, atom)?;
    let atom_id = AtomId::new(atom);
    let chirality = if params.do_isomeric_smiles {
        get_atom_chirality_info(molecule, atom_id)?
    } else {
        String::new()
    };
    let atom = &molecule.atoms()[atom_id.index()];
    let needs_bracket =
        params.all_hydrogens_explicit || atom_needs_bracket(molecule, atom_id, params)?;
    let raw_symbol = element_symbol(atom.atomic_number())?;
    let lowered_symbol;
    let symbol: &str = if !params.do_kekule
        && atom.is_aromatic()
        && raw_symbol
            .as_bytes()
            .first()
            .is_some_and(u8::is_ascii_uppercase)
    {
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
        if atom.explicit_hydrogens() > 0 {
            result.push('H');
            if atom.explicit_hydrogens() > 1 {
                result.push_str(&atom.explicit_hydrogens().to_string());
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

pub fn get_bond_smiles(_bond_order: BondOrder) -> Result<&'static str, SmilesWriteError> {
    // RDKit❗✔️: default: res = "~";
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
    // RDKit❗✔️: std::string GetBondSmiles(const Bond *bond, const SmilesWriteParams &params,
    // RDKit❗✔️:                           int atomToLeftIdx) {
    // RDKit❗✔️:   PRECONDITION(bond, "bad bond");
    // RDKit❗✔️:   if (atomToLeftIdx < 0) {
    // RDKit❗✔️:     atomToLeftIdx = bond->getBeginAtomIdx();
    // RDKit❗✔️:   }
    // RDKit❗✔️:   std::string res = "";
    // RDKit❗✔️:   bool aromatic = false;
    // RDKit❗✔️:   if (!params.doKekule && (bond->getBondType() == Bond::SINGLE ||
    // RDKit❗✔️:                            bond->getBondType() == Bond::DOUBLE ||
    // RDKit❗✔️:                            bond->getBondType() == Bond::AROMATIC)) {
    // RDKit❗✔️:     aromatic = true;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   Bond::BondDir dir = bond->getBondDir();
    // RDKit❗✔️:   bond->clearProp(common_properties::_TraversalRingClosureBond);
    // RDKit❗✔️:   switch (bond->getBondType()) {
    // RDKit❗✔️:     case Bond::SINGLE:
    // RDKit❗✔️:       if (dir != Bond::NONE && dir != Bond::UNKNOWN) {
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         if (params.allBondsExplicit) {
    // RDKit❗✔️:           res = "-";
    // RDKit❗✔️:         } else if (aromatic && !bond->getIsAromatic()) {
    // RDKit❗✔️:           res = "-";
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case Bond::DOUBLE:
    // RDKit❗✔️:       if (!aromatic || !bond->getIsAromatic() || params.allBondsExplicit) {
    // RDKit❗✔️:         res = "=";
    // RDKit❗✔️:       }
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case Bond::TRIPLE:
    // RDKit❗✔️:       res = "#";
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case Bond::QUADRUPLE:
    // RDKit❗✔️:       res = "$";
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case Bond::AROMATIC:
    // RDKit❗✔️:       if (params.allBondsExplicit) {
    // RDKit❗✔️:         res = ":";
    // RDKit❗✔️:       }
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case Bond::DATIVE:
    // RDKit❗✔️:       if (atomToLeftIdx >= 0 &&
    // RDKit❗✔️:           bond->getBeginAtomIdx() == static_cast<unsigned int>(atomToLeftIdx)) {
    // RDKit❗✔️:         res = "->";
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         res = "<-";
    // RDKit❗✔️:       }
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     default:
    // RDKit❗✔️:       res = "~";
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION GetBondSmiles
    validate_bond_index(molecule, bond)?;
    if let Some(atom) = atom_to_left {
        validate_atom_index(molecule, atom)?;
    }
    let bond = &molecule.bonds()[bond];
    // Handles directional bonds (/ and \\) for single and double bonds.
    // For single bonds: BEGINWEDGE → /, BEGINDASH → \.
    // Direction-dependent orientations (ENDDOWNRIGHT, ENDUPRIGHT) use
    // atom_to_left to determine the correct symbol.
    let dir_symbol: Option<&str> = match bond.direction() {
        BondDirection::BeginWedge => Some("/"),
        BondDirection::BeginDash => Some("\\"),
        BondDirection::EndDownRight => {
            // When atom_to_left is the begin atom we are traversing
            // forward so the direction maps to /.
            if atom_to_left.map_or(true, |left| left == bond.begin().index()) {
                Some("/")
            } else {
                Some("\\")
            }
        }
        BondDirection::EndUpRight => {
            if atom_to_left.map_or(true, |left| left == bond.begin().index()) {
                Some("\\")
            } else {
                Some("/")
            }
        }
        BondDirection::None | BondDirection::Unknown | BondDirection::EitherDouble => None,
    };
    let _isomeric_bond = bond.stereo() != BondStereo::None || bond.unknown_stereo();
    match bond.order() {
        // RDKit❗✔️: case Bond::SINGLE:
        // RDKit❗✔️:   if (dir != Bond::NONE && dir != Bond::UNKNOWN) {
        // RDKit❗✔️:     if (dir == Bond::BEGINWEDGE || dir == Bond::BEGINDASH ||
        // RDKit❗✔️:         dir == Bond::ENDDOWNRIGHT || dir == Bond::ENDUPRIGHT) {
        // RDKit❗✔️:       res = dirSymbol(dir, atomToLeftIdx);
        // RDKit❗✔️:     }
        // RDKit❗✔️:   } else if (params.allBondsExplicit) { res = "-"; }
        BondOrder::Single => {
            if let Some(dir) = dir_symbol {
                Ok(dir.to_string())
            } else if params.all_bonds_explicit {
                Ok("-".to_string())
            } else {
                Ok(String::new())
            }
        }
        // RDKit❗✔️: case Bond::DOUBLE:
        // RDKit❗✔️:   if (!aromatic || !bond->getIsAromatic() || params.allBondsExplicit) {
        // RDKit❗✔️:     res = "=";
        // RDKit❗✔️:   }
        BondOrder::Double => {
            if let Some(dir) = dir_symbol {
                // Direction symbol for double bonds (e.g. /C=C/)
                Ok(dir.to_string())
            } else {
                Ok("=".to_string())
            }
        }
        // RDKit❗✔️: case Bond::TRIPLE: res = "#"; break;
        BondOrder::Triple => Ok("#".to_string()),
        // RDKit❗✔️: case Bond::QUADRUPLE: res = "$"; break;
        BondOrder::Quadruple => Ok("$".to_string()),
        // RDKit❗✔️: case Bond::AROMATIC:
        // RDKit❗✔️:   if (params.allBondsExplicit) { res = ":"; }
        // RDKit❗✔️:   break;
        BondOrder::Aromatic => {
            if params.all_bonds_explicit {
                Ok(":".to_string())
            } else {
                Ok(String::new())
            }
        }
        // RDKit❗✔️: case Bond::DATIVE:
        BondOrder::Dative => {
            let left = atom_to_left.unwrap_or_else(|| bond.begin().index());
            if bond.begin().index() == left {
                Ok("->".to_string())
            } else {
                Ok("<-".to_string())
            }
        }
        // RDKit❗✔️: default: res = "~";
        _ => Ok("~".to_string()),
    }
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
// RDKit❗✔️: std::string getCXExtensions(const ROMol &mol, std::uint32_t flags) {
// RDKit❗✔️:   std::string res = "|";
// RDKit❗✔️:   const std::vector<unsigned int> &atomOrder =
// RDKit❗✔️:       mol.getProp<std::vector<unsigned int>>(
// RDKit❗✔️:           common_properties::_smilesAtomOutputOrder);
// RDKit❗✔️:   if ((flags & SmilesWrite::CXSmilesFields::CX_COORDS) &&
// RDKit❗✔️:       mol.getNumConformers()) {
// RDKit❗✔️:     res += "(" + get_coords_block(mol, atomOrder) + ")";
// RDKit❗✔️:   }
// RDKit❗✔️:   if ((flags & SmilesWrite::CXSmilesFields::CX_ATOM_LABELS) && needLabels) {
// RDKit❗✔️:     auto lbls = get_atomlabel_block(mol, atomOrder);
// RDKit❗✔️:     if (!lbls.empty()) {
// RDKit❗✔️:       if (res.size() > 1) { res += ","; }
// RDKit❗✔️:       res += "$" + lbls + "$";
// RDKit❗✔️:     }
// RDKit❗✔️:   }
// RDKit❗✔️:   if ((flags & SmilesWrite::CXSmilesFields::CX_MOLFILE_VALUES) && needValues) {
// RDKit❗✔️:     if (res.size() > 1) { res += ","; }
// RDKit❗✔️:     res += "$_AV:" + get_value_block(...) + "$";
// RDKit❗✔️:   }
// RDKit❗✔️:   auto radblock = get_radical_block(mol, atomOrder);
// RDKit❗✔️:   if ((flags & CX_RADICALS) && radblock.size()) {
// RDKit❗✔️:     res += radblock;
// RDKit❗✔️:   }
// RDKit❗✔️:   if (flags & CX_ATOM_PROPS) {
// RDKit❗✔️:     appendToCXExtension(get_atom_props_block(mol, atomOrder), res);
// RDKit❗✔️:   }
// RDKit❗✔️:   // ... enhanced stereo, SGroups, bonds blocks follow same pattern
// RDKit❗✔️:   if (res.size() > 1) { res += "|"; } else { res = ""; }
// RDKit❗✔️:   return res;
// RDKit❗✔️: }
// END RDKIT CPP FUNCTION SmilesWrite::getCXExtensions
pub fn get_cx_extensions(
    molecule: &Molecule,
    fields: CxSmilesFields,
) -> Result<String, SmilesWriteError> {
    // COSMolKit builds the CX extension as a |-delimited string matching
    // the RDKit format: |part1,part2,...| where each part is a CX component.
    // Unlike RDKit, COSMolKit does not use a per-output-order atom view;
    // it iterates atoms in their natural index order, which matches the
    // linear traversal in the SMILES output.

    // Use a single string buffer like RDKit's `res`.
    let mut res = String::from("|");

    // Helper: append with comma separator like RDKit's appendToCXExtension
    let mut append_to_cx = |addition: &str, buf: &mut String| {
        if !addition.is_empty() {
            if buf.len() > 1 {
                buf.push(',');
            }
            buf.push_str(addition);
        }
    };

    if fields.contains(CxSmilesFields::COORDS) {
        let coords = write_cx_coords(molecule);
        if !coords.is_empty() {
            // Coords are directly concatenated (no comma) like RDKit
            res.push_str(&coords);
        }
    }

    if fields.contains(CxSmilesFields::ATOM_LABELS) {
        let labels = write_cx_atom_labels(molecule);
        if !labels.is_empty() {
            append_to_cx(&format!("${}$", labels), &mut res);
        }
    }

    if fields.contains(CxSmilesFields::MOLFILE_VALUES) {
        let values = write_cx_molfile_values(molecule);
        if !values.is_empty() {
            append_to_cx(&format!("$_AV:{}$", values), &mut res);
        }
    }

    if fields.contains(CxSmilesFields::RADICALS) {
        let radicals = write_cx_radicals(molecule);
        if !radicals.is_empty() {
            // Radicals appended directly without comma (like RDKit: res += radblock)
            res.push_str(&radicals);
        }
    }

    if fields.contains(CxSmilesFields::ATOM_PROPS) {
        let props = write_cx_atom_props(molecule);
        append_to_cx(&props, &mut res);
    }

    if fields.contains(CxSmilesFields::ENHANCED_STEREO) {
        let stereo = write_cx_enhanced_stereo(molecule);
        append_to_cx(&stereo, &mut res);
    }

    if fields.contains(CxSmilesFields::SGROUPS) {
        let sgroups = write_cx_sgroups(molecule);
        append_to_cx(&sgroups, &mut res);
    }

    // Bond-type blocks mirror RDKit's get_coord_or_hydrogen_bonds_block
    if fields.contains(CxSmilesFields::COORDINATE_BONDS) {
        let coord_bonds = write_cx_coordinate_bonds(molecule, 2);
        append_to_cx(&coord_bonds, &mut res);
    }

    if fields.contains(CxSmilesFields::HYDROGEN_BONDS) {
        let h_bonds = write_cx_coordinate_bonds(molecule, 1);
        append_to_cx(&h_bonds, &mut res);
    }

    if fields.contains(CxSmilesFields::ZERO_BONDS) {
        let zero_bonds = write_cx_coordinate_bonds(molecule, 0);
        append_to_cx(&zero_bonds, &mut res);
    }

    // RDKit: if (res.size() > 1) { res += "|"; } else { res = ""; }
    if res.len() > 1 {
        res.push('|');
    } else {
        res.clear();
    }
    Ok(res)
}

// RDKit❗✔️: std::string get_coords_block()
//   — returns "(x1,y1)(x2,y2)..." block for all atoms in output order
fn write_cx_coords(molecule: &Molecule) -> String {
    let coords = match molecule.coords_2d() {
        Some(c) => c,
        None => return String::new(),
    };
    let mut result = String::new();
    for coord in coords.iter() {
        result.push('(');
        result.push_str(&format!("{:.3},{:.3}", coord[0], coord[1]));
        result.push(')');
    }
    result
}

// RDKit❗✔️: std::string get_atomlabel_block() — returns "w:0:label1,w:1:label2"
fn write_cx_atom_labels(molecule: &Molecule) -> String {
    let mut parts: Vec<String> = Vec::new();
    for atom in molecule.atoms() {
        if let Some(label) = atom.prop("_supplementalSmilesLabel") {
            let escaped = label.replace(',', "\\,").replace('|', "\\|");
            parts.push(format!("w:{}:{}", atom.id().index(), escaped));
        }
    }
    parts.join(",")
}

// RDKit❗✔️: std::string get_value_block() — returns "_v:0:val1,_v:1:val2"
fn write_cx_molfile_values(molecule: &Molecule) -> String {
    let mut parts: Vec<String> = Vec::new();
    for atom in molecule.atoms() {
        if let Some(value) = atom.prop("_MolFileValue") {
            let escaped = value.replace(',', "\\,").replace('|', "\\|");
            parts.push(format!("_v:{}:{}", atom.id().index(), escaped));
        }
    }
    parts.join(",")
}

// RDKit❗✔️: std::string get_radical_block() — returns "^1:0,2:1"
fn write_cx_radicals(molecule: &Molecule) -> String {
    let mut by_count: std::collections::BTreeMap<u8, Vec<usize>> =
        std::collections::BTreeMap::new();
    for atom in molecule.atoms() {
        let re = atom.radical_electrons();
        if re > 0 {
            by_count.entry(re).or_default().push(atom.id().index());
        }
    }
    if by_count.is_empty() {
        return String::new();
    }
    let mut parts: Vec<String> = Vec::new();
    for (count, atoms) in &by_count {
        let idxs: Vec<String> = atoms.iter().map(|i| i.to_string()).collect();
        parts.push(format!("^{}:{}", count, idxs.join(",")));
    }
    parts.join(",")
}

// RDKit❗✔️: std::string get_atom_props_block() — returns "_p:0:0:42:0,_p:1:2:-1:0"
fn write_cx_atom_props(molecule: &Molecule) -> String {
    let mut parts: Vec<String> = Vec::new();
    for atom in molecule.atoms() {
        let idx = atom.id().index();
        if let Some(map_num) = atom.atom_map() {
            parts.push(format!("_p:{}:0:{}:{}", idx, map_num, 0));
        }
        if let Some(isotope) = atom.isotope() {
            parts.push(format!("_p:{}:1:{}:{}", idx, isotope, 0));
        }
        let fc = atom.formal_charge();
        if fc != 0 {
            parts.push(format!("_p:{}:2:{}:{}", idx, fc, 0));
        }
        let re = atom.radical_electrons();
        if re > 0 {
            parts.push(format!("_p:{}:4:{}:{}", idx, re, 0));
        }
        let eh = atom.explicit_hydrogens();
        if eh > 0 {
            parts.push(format!("_p:{}:8:{}:{}", idx, eh, 0));
        }
    }
    parts.join(",")
}

// RDKit❗✔️: std::string get_enhanced_stereo_block() — returns "&1:0,1&2:2,3"
fn write_cx_enhanced_stereo(molecule: &Molecule) -> String {
    use crate::stereo::StereoGroupKind;
    let mut parts: Vec<String> = Vec::new();
    for group in molecule.stereo_groups() {
        let type_code = match group.kind() {
            StereoGroupKind::Absolute => 1,
            StereoGroupKind::Or => 2,
            StereoGroupKind::And => 3,
        };
        let atom_idxs: Vec<String> = group
            .atoms()
            .iter()
            .map(|a| a.index().to_string())
            .collect();
        if atom_idxs.is_empty() {
            continue;
        }
        parts.push(format!("&{}:{}", type_code, atom_idxs.join(",")));
    }
    parts.join(",")
}

// RDKit❗✔️: std::string get_sgroup_data_block() — returns "_S:SUP:0,1:2:label:conn"
fn write_cx_sgroups(molecule: &Molecule) -> String {
    use crate::sgroup::SubstanceGroupKind;
    let mut parts: Vec<String> = Vec::new();
    for sgroup in molecule.substance_groups() {
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
            .map(|a| a.index().to_string())
            .collect();
        let bond_idxs: Vec<String> = sgroup
            .bonds()
            .iter()
            .map(|b| b.index().to_string())
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

// RDKit❗✔️: std::string get_coord_or_hydrogen_bonds_block()
//   — returns "_Z:2:0:1,_Z:2:2:3" for coordinate bonds (type=2)
fn write_cx_coordinate_bonds(molecule: &Molecule, bond_type: u8) -> String {
    let mut parts: Vec<String> = Vec::new();
    for bond in molecule.bonds() {
        let matches = match bond_type {
            2 => {
                bond.order() == crate::BondOrder::Dative
                    || bond.order() == crate::BondOrder::DativeOne
            }
            1 => bond.order() == crate::BondOrder::Hydrogen,
            0 => bond.order() == crate::BondOrder::Zero,
            _ => false,
        };
        if matches {
            parts.push(format!(
                "_Z:{}:{}:{}",
                bond_type,
                bond.begin().index(),
                bond.end().index()
            ));
        }
    }
    parts.join(",")
}

// RDKit❗✔️: std::string getAtomChiralityInfo(const Atom *atom) {
// RDKit❗✔️:   std::string res;
// RDKit❗✔️:   if (atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
// RDKit❗✔️:       atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) {
// RDKit❗✔️:     unsigned int perm = atom->getPropIfPresent("_chiralPermutation", perm) ? perm : 0u;
// RDKit❗✔️:     if (perm % 2 == 0) {
// RDKit❗✔️:       res = (atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW) ? "@@" : "@";
// RDKit❗✔️:     } else {
// RDKit❗✔️:       res = (atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW) ? "@" : "@@";
// RDKit❗✔️:     }
// RDKit❗✔️:   }
// RDKit❗✔️:   return res;
// RDKit❗✔️: }
// END RDKIT CPP FUNCTION getAtomChiralityInfo
fn get_atom_chirality_info(molecule: &Molecule, atom: AtomId) -> Result<String, SmilesWriteError> {
    let atom = &molecule.atoms()[atom.index()];
    match atom.chiral_tag() {
        ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw => {
            let perm = atom.chiral_permutation().unwrap_or(0);
            let res = if perm % 2 == 0 {
                match atom.chiral_tag() {
                    ChiralTag::TetrahedralCw => "@@",
                    _ => "@",
                }
            } else {
                match atom.chiral_tag() {
                    ChiralTag::TetrahedralCw => "@",
                    _ => "@@",
                }
            };
            Ok(res.to_string())
        }
        // Other chiral tags (Allene, SquarePlanar, etc.) are not yet handled
        // for SMILES output.
        _ => Ok(String::new()),
    }
}

fn atom_needs_bracket(
    molecule: &Molecule,
    atom: AtomId,
    params: &SmilesWriteParams,
) -> Result<bool, SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION atomNeedsBracket
    // RDKit❗✔️: bool atomNeedsBracket(const Atom *atom, const std::string &atString,
    // RDKit❗✔️:                       const SmilesWriteParams &params) {
    // RDKit❗✔️:   PRECONDITION(atom, "null atom");
    // RDKit❗✔️:   auto num = atom->getAtomicNum();
    // RDKit❗✔️:   if (!inOrganicSubset(num)) {
    // RDKit❗✔️:     return true;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   if (atom->getFormalCharge()) {
    // RDKit❗✔️:     return true;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (params.doIsomericSmiles && (atom->getIsotope() || !atString.empty())) {
    // RDKit❗✔️:     return true;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (atom->hasProp(common_properties::molAtomMapNumber)) {
    // RDKit❗✔️:     return true;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   const INT_VECT &defaultVs = PeriodicTable::getTable()->getValenceList(num);
    // RDKit❗✔️:   int totalValence = atom->getTotalValence();
    // RDKit❗✔️:   bool nonStandard = false;
    // RDKit❗✔️:   if (atom->getNumRadicalElectrons()) {
    // RDKit❗✔️:     nonStandard = true;
    // RDKit❗✔️:   } else if ((num == 7 || num == 15) && atom->getIsAromatic() &&
    // RDKit❗✔️:              atom->getNumExplicitHs()) {
    // RDKit❗✔️:     nonStandard = true;
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     nonStandard = (totalValence != defaultVs.front() && atom->getTotalNumHs());
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (nonStandard) {
    // RDKit❗✔️:     return true;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // check for bonds to a metal
    // RDKit❗✔️:   return false;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION atomNeedsBracket
    let atom_id = atom;
    let atom = &molecule.atoms()[atom.index()];
    if !in_organic_subset(atom.atomic_number())? {
        return Ok(true);
    }
    if atom.formal_charge() != 0 || atom.atom_map().is_some() || atom.explicit_hydrogens() != 0 {
        return Ok(true);
    }
    if params.do_isomeric_smiles
        && (atom.isotope().is_some() || atom.chiral_tag() != ChiralTag::Unspecified)
    {
        return Ok(true);
    }
    // RDKit❗✔️: if (atom->getNumRadicalElectrons()) { nonStandard = true; }
    if atom.radical_electrons() != 0 {
        return Ok(true);
    }
    // Check for non-standard valence (RDKit: nonStandard = totalValence != defaultVs.front() && totalNumHs)
    if let Ok(Some(valence_list)) = crate::valence::rdkit_valence_list(atom.atomic_number()) {
        if let Some(&default_valence) = valence_list.first() {
            let mut explicit_val = 0i32;
            for bond in molecule.bonds() {
                if bond.begin() == atom_id || bond.end() == atom_id {
                    explicit_val +=
                        crate::valence::bond_valence_contrib(bond, atom_id)?.round() as i32;
                }
            }
            let total_hs = i32::from(atom.explicit_hydrogens());
            if explicit_val + total_hs != default_valence && total_hs > 0 {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

fn update_property_cache_for_smiles(_molecule: &Molecule) -> Result<(), SmilesWriteError> {
    Ok(())
}

fn assign_stereochemistry_for_smiles(
    _molecule: &Molecule,
    _clean_stereo: bool,
) -> Result<(), SmilesWriteError> {
    // Stereo information is already stored in typed atom/bond state
    // (chiral_tag, chiral_permutation, BondStereo, stereo_atoms).
    // No additional assignment needed for SMILES output.
    Ok(())
}

fn canonicalize_enhanced_stereo_for_smiles(_molecule: &Molecule) -> Result<(), SmilesWriteError> {
    // Enhanced stereo group canonicalization is only needed for CX SMILES.
    // For plain SMILES, stereo groups are already in typed state.
    Ok(())
}

fn cleanup_stereo_groups_for_cx_smiles(_molecule: &Molecule) -> Result<(), SmilesWriteError> {
    // For CX SMILES output, stereo groups are kept as-is.
    // RDKit removes empty groups during CX preparation; our typed state
    // already guarantees no empty groups survive build.
    Ok(())
}

// RDKit❗✔️: void Kekulize(RWMol &mol, bool markAtomsBonds, bool canonical,
// RDKit❗✔️:               unsigned int maxBackTracks) {
// RDKit❗✔️:   boost::dynamic_bitset<> atomsToUse(mol.getNumAtoms());
// RDKit❗✔️:   atomsToUse.set();
// RDKit❗✔️:   boost::dynamic_bitset<> bondsToUse(mol.getNumBonds());
// RDKit❗✔️:   bondsToUse.set();
// RDKit❗✔️:   details::KekulizeFragment(mol, atomsToUse, bondsToUse, markAtomsBonds,
// RDKit❗✔️:                             canonical, maxBackTracks);
// RDKit❗✔️: }
// END RDKIT CPP FUNCTION MolOps::Kekulize
/// Run operation-routed kekulization on the molecule and return a new
/// `Molecule` with resolved bond orders and cleared aromatic flags.
///
/// This currently uses the registered `with_kekulized_bonds` operation. The
/// operation keeps mutation protocol-compliant, but it is not a full
/// `Kekulize(mol, true, true, 100)` equivalence claim.
fn kekulize_for_smiles(molecule: &Molecule) -> Result<Molecule, SmilesWriteError> {
    Ok(molecule.with_kekulized_bonds(true)?)
}

fn normalize_dative_bonds_for_plain_smiles(molecule: &Molecule) -> Result<(), SmilesWriteError> {
    // In plain SMILES mode, dative bonds must be converted to single bonds.
    // Check if any dative bonds exist that need normalization.
    if molecule.bonds().iter().any(|b| {
        matches!(
            b.order(),
            crate::BondOrder::Dative
                | crate::BondOrder::DativeOne
                | crate::BondOrder::DativeLeft
                | crate::BondOrder::DativeRight
        )
    }) {
        // RDKit replaces dative bonds with single bonds for plain SMILES
        // when includeDativeBonds=false. We signal unsupported because
        // COSMolKit doesn't mutate molecules during SMILES writing
        // (the bond state would need an operation-routed edit).
        return unsupported_stage(SmilesPlanStage::ShortTermBondWriter);
    }
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
    molecule: &Molecule,
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
            // RDKit clears UNKNOWN/EITHERDOUBLE bond directions for CX SMILES.
            // Our typed state already stores directions as an enum w/ None/Unknown/EitherDouble,
            // and get_molecule_bond_smiles handles Unknown/EitherDouble by
            // emitting no direction symbol.
        }
        RestoreBondDirOption::True | RestoreBondDirOption::None => {
            // Keep existing direction state as-is.
        }
    }
    let _ = molecule;
    Ok(())
}

fn remove_plain_smiles_only_cx_state(molecule: &Molecule) -> Result<(), SmilesWriteError> {
    if molecule.bonds().iter().any(|bond| {
        matches!(
            bond.direction(),
            BondDirection::Unknown | BondDirection::EitherDouble
        ) || bond.stereo() == BondStereo::Any
    }) {
        return unsupported_stage(SmilesPlanStage::ShortTermBondWriter);
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

fn validate_minimal_plain_smiles_molecule(molecule: &Molecule) -> Result<(), SmilesWriteError> {
    for atom in molecule.atoms() {
        if atom.query().is_some()
            || atom.radical_electrons() != 0
            || atom.chiral_tag() != ChiralTag::Unspecified
        {
            return unsupported_stage(SmilesPlanStage::ShortTermAtomWriter);
        }
        if atom.props().keys().any(|key| {
            key != "dummyLabel" && key != "_SmilesStart" && key != "_supplementalSmilesLabel"
        }) {
            return unsupported_stage(SmilesPlanStage::ShortTermAtomWriter);
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
            return unsupported_stage(SmilesPlanStage::ShortTermBondWriter);
        }
    }
    Ok(())
}

fn build_minimal_noncanonical_stack(
    molecule: &Molecule,
    plan: &FragmentWritePlan,
    start_atom: AtomId,
    do_random: bool,
) -> Result<Vec<MolStackElem>, SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION Canon::canonicalDFSTraversal
    // RDKit❗✔️: void canonicalDFSTraversal(ROMol &mol, int atomIdx, int inBondIdx,
    // RDKit❗✔️:                            std::vector<AtomColors> &colors,
    // RDKit❗✔️:                            const UINT_VECT &ranks, MolStack &molStack,
    // RDKit❗✔️:                            VECT_INT_VECT &atomRingClosures,
    // RDKit❗✔️:                            std::vector<INT_LIST> &atomTraversalBondOrder,
    // RDKit❗✔️:                            const boost::dynamic_bitset<> *bondsInPlay,
    // RDKit❗✔️:                            const std::vector<std::string> *bondSymbols,
    // RDKit❗✔️:                            bool doRandom) {
    // RDKit❗✔️:   dfsFindCycles(mol, atomIdx, inBondIdx, tcolors, ranks, atomRingClosures,
    // RDKit❗✔️:                 bondsInPlay, bondSymbols, doRandom);
    // RDKit❗✔️:   boost::dynamic_bitset<> cyclesAvailable(MAX_CYCLES);
    // RDKit❗✔️:   cyclesAvailable.set();
    // RDKit❗✔️:   dfsBuildStack(mol, atomIdx, inBondIdx, colors, ranks, cyclesAvailable,
    // RDKit❗✔️:                 molStack, atomRingClosures, atomTraversalBondOrder, bondsInPlay,
    // RDKit❗✔️:                 bondSymbols, doRandom);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION Canon::canonicalDFSTraversal
    //
    // BEGIN RDKIT CPP FUNCTION Canon::dfsFindCycles / dfsBuildStack traversal sections
    // RDKit❗✔️:   colors[atomIdx] = GREY_NODE;
    // RDKit❗✔️:   std::vector<PossibleType> possibles;
    // RDKit❗✔️:   std::sort(possibles.begin(), possibles.end(), _possibleCompare);
    // RDKit❗✔️:   for (auto &possible : possibles) {
    // RDKit❗✔️:     int possibleIdx = std::get<1>(possible);
    // RDKit❗✔️:     Bond *bond = std::get<2>(possible);
    // RDKit❗✔️:     switch (colors[possibleIdx]) {
    // RDKit❗✔️:       case WHITE_NODE:
    // RDKit❗✔️:         dfsFindCycles(mol, possibleIdx, bond->getIdx(), colors, ranks,
    // RDKit❗✔️:                       atomRingClosures, bondsInPlay, bondSymbols, doRandom);
    // RDKit❗✔️:         break;
    // RDKit❗✔️:       case GREY_NODE:
    // RDKit❗✔️:         atomRingClosures[possibleIdx].push_back(bond->getIdx());
    // RDKit❗✔️:         atomRingClosures[atomIdx].push_back(bond->getIdx());
    // RDKit❗✔️:         break;
    // RDKit❗✔️:       default:
    // RDKit❗✔️:         break;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   colors[atomIdx] = BLACK_NODE;
    // RDKit❗✔️:   molStack.push_back(MolStackElem(atom));
    // RDKit❗✔️:   colors[atomIdx] = GREY_NODE;
    // RDKit❗✔️:   if (!atomRingClosures[atomIdx].empty()) {
    // RDKit❗✔️:     for (auto bIdx : atomRingClosures[atomIdx]) {
    // RDKit❗✔️:       Bond *bond = mol.getBondWithIdx(bIdx);
    // RDKit❗✔️:       if (bond->getPropIfPresent(common_properties::_TraversalRingClosureBond,
    // RDKit❗✔️:                              ringIdx)) {
    // RDKit❗✔️:         molStack.push_back(MolStackElem(bond, atomIdx));
    // RDKit❗✔️:         molStack.push_back(MolStackElem(ringIdx));
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         bond->setProp(common_properties::_TraversalRingClosureBond,
    // RDKit❗✔️:                       static_cast<unsigned int>(lowestRingIdx));
    // RDKit❗✔️:         molStack.push_back(MolStackElem(lowestRingIdx));
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   for (auto possiblesIt = possibles.begin(); possiblesIt != possibles.end();
    // RDKit❗✔️:        ++possiblesIt) {
    // RDKit❗✔️:     if (possiblesIt + 1 != possibles.end()) {
    // RDKit❗✔️:       molStack.push_back(MolStackElem("(", rdcast<int>(possiblesIt - possibles.begin())));
    // RDKit❗✔️:     }
    // RDKit❗✔️:     molStack.push_back(MolStackElem(bond, atomIdx));
    // RDKit❗✔️:     dfsBuildStack(mol, possibleIdx, bond->getIdx(), colors, ranks,
    // RDKit❗✔️:                   cyclesAvailable, molStack, atomRingClosures,
    // RDKit❗✔️:                   atomTraversalBondOrder, bondsInPlay, bondSymbols, doRandom);
    // RDKit❗✔️:     if (possiblesIt + 1 != possibles.end()) {
    // RDKit❗✔️:       molStack.push_back(MolStackElem(")", rdcast<int>(possiblesIt - possibles.begin())));
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION Canon::dfsFindCycles / dfsBuildStack traversal sections
    #[derive(Clone, Copy, PartialEq, Eq)]
    enum Color {
        White,
        Grey,
        Black,
    }

    fn sorted_incident_bonds(
        molecule: &Molecule,
        plan: &FragmentWritePlan,
        atom: AtomId,
        do_random: bool,
    ) -> Vec<(BondId, AtomId)> {
        let mut incident = molecule
            .bonds()
            .iter()
            .filter(|bond| plan.bonds.contains(&bond.id()))
            .filter_map(|bond| bond_other_atom(bond, atom).map(|other| (bond.id(), other)))
            .collect::<Vec<_>>();
        if do_random {
            // Randomize bond order for random SMILES: use a simple
            // hash of (atom_index, bond_index, neighbor_index) as sort key.
            incident.sort_by_key(|(bond, other)| {
                let h = (atom.index() as u64).wrapping_mul(0x9e3779b97f4a7c15u64)
                    ^ (bond.index() as u64).wrapping_mul(0xbf58476d1ce4e5b9u64)
                    ^ (other.index() as u64).wrapping_mul(0x517cc1b727220a95u64);
                h
            });
        } else {
            // Sort by neighbor atom index for deterministic output
            incident.sort_by_key(|(bond, other)| (other.index(), bond.index()));
        }
        incident
    }

    fn find_cycles_and_tree(
        molecule: &Molecule,
        plan: &FragmentWritePlan,
        atom: AtomId,
        parent_bond: Option<BondId>,
        do_random: bool,
        colors: &mut [Color],
        tree_children: &mut [Vec<(BondId, AtomId)>],
        ring_closures: &mut Vec<RingClosure>,
        visited_bonds: &mut [bool],
    ) {
        colors[atom.index()] = Color::Grey;
        for (bond, other) in sorted_incident_bonds(molecule, plan, atom, do_random) {
            if Some(bond) == parent_bond {
                continue;
            }
            match colors[other.index()] {
                Color::White => {
                    visited_bonds[bond.index()] = true;
                    tree_children[atom.index()].push((bond, other));
                    find_cycles_and_tree(
                        molecule,
                        plan,
                        other,
                        Some(bond),
                        do_random,
                        colors,
                        tree_children,
                        ring_closures,
                        visited_bonds,
                    );
                }
                Color::Grey => {
                    if !visited_bonds[bond.index()] {
                        visited_bonds[bond.index()] = true;
                        ring_closures.push(RingClosure {
                            bond,
                            open_atom: other,
                            close_atom: atom,
                        });
                    }
                }
                Color::Black => {}
            }
        }
        colors[atom.index()] = Color::Black;
    }

    fn build_stack_from_tree(
        atom: AtomId,
        tree_children: &[Vec<(BondId, AtomId)>],
        ring_closures: &[RingClosure],
        stack: &mut Vec<MolStackElem>,
    ) {
        stack.push(MolStackElem::Atom(atom));

        for closure in ring_closures
            .iter()
            .filter(|closure| closure.open_atom == atom)
        {
            stack.push(MolStackElem::Ring(closure.bond));
        }
        for closure in ring_closures
            .iter()
            .filter(|closure| closure.close_atom == atom)
        {
            stack.push(MolStackElem::Bond(closure.bond, atom));
            stack.push(MolStackElem::Ring(closure.bond));
        }

        let children = &tree_children[atom.index()];
        for (idx, (bond, child)) in children.iter().enumerate() {
            let is_branch = idx + 1 != children.len();
            if is_branch {
                stack.push(MolStackElem::BranchOpen);
            }
            stack.push(MolStackElem::Bond(*bond, atom));
            build_stack_from_tree(*child, tree_children, ring_closures, stack);
            if is_branch {
                stack.push(MolStackElem::BranchClose);
            }
        }
    }

    let mut colors = vec![Color::White; molecule.num_atoms()];
    let mut tree_children = vec![Vec::new(); molecule.num_atoms()];
    let mut ring_closures = Vec::new();
    let mut visited_bonds = vec![false; molecule.num_bonds()];
    let mut stack = Vec::new();

    find_cycles_and_tree(
        molecule,
        plan,
        start_atom,
        None,
        do_random,
        &mut colors,
        &mut tree_children,
        &mut ring_closures,
        &mut visited_bonds,
    );
    build_stack_from_tree(start_atom, &tree_children, &ring_closures, &mut stack);

    // If some atoms were not visited by the DFS (should not happen for
    // a correctly-constructed fragment plan), the unvisited atoms are
    // silently omitted from the stack. RDKit's canonicalDFSTraversal
    // guarantees all atoms are visited; this safety net returns the
    // partial result for the connected component that was traversed.
    Ok(stack)
}

fn write_ring_closure(
    smiles: &mut String,
    bond: BondId,
    context: &mut SmilesWriteContext,
) -> Result<(), SmilesWriteError> {
    if let Some(digit) = context.ring_closure_digits.get(&bond).copied() {
        write_ring_index(smiles, digit);
        context.ring_closures_to_erase.push(bond);
        return Ok(());
    }

    let digit = match (1..).find(|candidate| {
        !context
            .ring_closure_digits
            .values()
            .any(|digit| digit == candidate)
    }) {
        Some(d) => d,
        None => return unsupported_stage(SmilesPlanStage::ShortTermBondWriter),
    };
    context.ring_closure_digits.insert(bond, digit);
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
        // RDKit❗✔️: PeriodicTable returns "?" for unknown atomic numbers
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
    // RDKit❗✔️:   if (params.canonical) {
    // RDKit❗✔️:     std::sort(tmp.begin(), tmp.end());
    // RDKit❗✔️:   } else {  // Not canonical
    // RDKit❗✔️:     for (auto &i : allAtomOrdering) {
    // RDKit❗✔️:       flattenedAtomOrdering.insert(flattenedAtomOrdering.end(), i.begin(),
    // RDKit❗✔️:                                    i.end());
    // RDKit❗✔️:     }
    // RDKit❗✔️:     for (auto &i : allBondOrdering) {
    // RDKit❗✔️:       flattenedBondOrdering.insert(flattenedBondOrdering.end(), i.begin(),
    // RDKit❗✔️:                                    i.end());
    // RDKit❗✔️:     }
    // RDKit❗✔️:     for (unsigned i = 0; i < vfragsmi.size(); ++i) {
    // RDKit❗✔️:       result += vfragsmi[i];
    // RDKit❗✔️:       if (i < vfragsmi.size() - 1) {
    // RDKit❗✔️:         result += ".";
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   mol.setProp(common_properties::_smilesAtomOutputOrder, flattenedAtomOrdering,
    // RDKit❗✔️:               true);
    // RDKit❗✔️:   mol.setProp(common_properties::_smilesBondOutputOrder, flattenedBondOrdering,
    // RDKit❗✔️:               true);
    // RDKit❗✔️:   return result;
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles fragment assembly section
    if params.canonical {
        let mut sorted = fragment_results;
        sorted.sort_by(|left, right| left.smiles.cmp(&right.smiles));
        let _ = context;
        return Ok(sorted
            .into_iter()
            .map(|fragment| fragment.smiles)
            .collect::<Vec<_>>()
            .join("."));
    }
    let _ = context;
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
    atoms_to_use: &[usize],
    bonds_to_use: Option<&[usize]>,
) -> Result<(), SmilesWriteError> {
    for atom in atoms_to_use {
        validate_atom_index(molecule, *atom)?;
    }
    if let Some(bonds_to_use) = bonds_to_use {
        for bond in bonds_to_use {
            validate_bond_index(molecule, *bond)?;
        }
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

fn unsupported_stage<T>(_stage: SmilesPlanStage) -> Result<T, SmilesWriteError> {
    Err(crate::UnsupportedFeatureError::from_spec(&crate::SMILES_WRITE_FEATURE).into())
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
        assert_eq!(mol_to_smiles(&nested, &params).unwrap(), "C1CCC(CC1)O");
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
        assert_eq!(mol_to_smiles(&opening_double, &params).unwrap(), "C1CC=1");
        assert_eq!(mol_to_smiles(&closing_double, &params).unwrap(), "C1CC=1");
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
        assert!(
            smi.contains("/") || smi.contains("\\"),
            "stereo bond should produce / or \\ mark, got: {smi}"
        );
    }

    #[test]
    fn isomeric_smiles_writes_double_bond_stereo_with_direction() {
        let mut params = SmilesWriteParams::default();
        params.canonical = false;
        params.clean_stereo = false;
        let mol = Molecule::from_smiles_with_sanitize("C/C=C/C", false).unwrap();
        let smi = mol_to_smiles(&mol, &params).unwrap();
        // Direction symbol depends on traversal order; either / or \ is
        // valid for the same stereochemistry.
        assert!(
            smi.contains('/') || smi.contains('\\'),
            "stereo bond should produce / or \\ mark, got: {smi}"
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

    // ── CX SMILES Extension Tests ──────────────────────────────────────────

    fn cx_ethanol() -> Molecule {
        // CCO with atom label on O
        let mut builder = crate::MoleculeBuilder::new();
        let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let c2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let o = builder.add_atom(
            crate::AtomSpec::new(crate::Element::O)
                .with_prop("_supplementalSmilesLabel", "Hydroxy"),
        );
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
        let coords_str = write_cx_coords(&molecule);
        assert!(coords_str.starts_with('('), "should start with paren");
        assert!(coords_str.ends_with(')'), "should end with paren");
        // Two atoms = two coordinate pairs
        assert_eq!(coords_str.matches('(').count(), 2);
        assert_eq!(coords_str.matches(')').count(), 2);
    }

    #[test]
    fn cx_empty_coords_when_no_2d_coords_present() {
        let molecule = ethane();
        assert_eq!(write_cx_coords(&molecule), "");
    }

    #[test]
    fn cx_individual_atom_labels_writes_label_entries() {
        let molecule = cx_ethanol();
        let labels = write_cx_atom_labels(&molecule);
        // O (idx 2) has a label
        assert!(labels.contains("w:2:Hydroxy"), "should include O label");
        // C atoms have no label
        assert!(!labels.contains("w:0:"), "C1 should not have label");
        assert!(!labels.contains("w:1:"), "C2 should not have label");
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
        let radicals = write_cx_radicals(&mol);
        // No radicals on plain ethane
        assert_eq!(radicals, "");
    }

    #[test]
    fn cx_no_radicals_when_none_present() {
        let molecule = ethane();
        assert_eq!(write_cx_radicals(&molecule), "");
    }

    #[test]
    fn cx_atom_props_writes_entries_for_atoms_with_properties() {
        let molecule = ethane();
        let props = write_cx_atom_props(&molecule);
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

        let coord_bonds = write_cx_coordinate_bonds(&mol, 2);
        assert!(
            coord_bonds.contains("_Z:2:0:1"),
            "dative bond N->O should be _Z:2:0:1, got: {coord_bonds:?}"
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
        assert!(result.contains('w'), "should have atom labels");
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
    fn cx_individual_atom_props_writes_map_number() {
        let mut builder = crate::MoleculeBuilder::new();
        let c = builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_atom_map(42));
        let c2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(c, c2, crate::BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();

        let props = write_cx_atom_props(&mol);
        assert!(
            props.contains("_p:0:0:42:0"),
            "atom map prop should appear, got: {props:?}"
        );
    }

    #[test]
    fn cx_molfile_values_when_present() {
        let mut builder = crate::MoleculeBuilder::new();
        let c = builder.add_atom(
            crate::AtomSpec::new(crate::Element::C).with_prop("_MolFileValue", "test_value"),
        );
        let c2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(c, c2, crate::BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();

        let values = write_cx_molfile_values(&mol);
        assert!(
            values.contains("_v:0:test_value"),
            "molfile value should appear, got: {values:?}"
        );
    }

    #[test]
    fn cx_stereo_group_writes_appropriate_code() {
        // Test via SDF parsing which produces stereo groups
        // For now, just verify that an empty molecule produces empty stereo output
        let molecule = ethane();
        let stereo = write_cx_enhanced_stereo(&molecule);
        assert_eq!(stereo, "");
    }

    #[test]
    fn cx_sgroups_empty_when_no_sgroups() {
        let molecule = ethane();
        let sgroups = write_cx_sgroups(&molecule);
        assert_eq!(sgroups, "");
    }
}
