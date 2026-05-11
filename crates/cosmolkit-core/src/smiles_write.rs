use crate::{AtomId, Bond, BondDirection, BondId, BondOrder, BondStereo, ChiralTag, Molecule};
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
    ShortTermFragmentTraversal,
    LongTermCanonicalRanking,
    LongTermIsomericStereo,
    LongTermKekule,
    LongTermCxExtensions,
    LongTermRandomSmiles,
    LongTermFragmentApi,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum SmilesWriteError {
    #[error("SMILES writing is not implemented")]
    NotImplemented,
    #[error(transparent)]
    UnsupportedFeature(#[from] crate::UnsupportedFeatureError),
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
    // RDKit❌❌:
    // RDKit❌❌:   int rootedAtAtom;
    // RDKit❌❌:   std::vector<int> fragsRootedAtAtom;
    // RDKit❌❌:   std::vector<std::vector<int>> fragsMolAtomMapping;
    // RDKit❌❌:   auto mols =
    // RDKit❌❌:       MolOps::getMolFrags(mol, false, nullptr, &fragsMolAtomMapping, false);
    // RDKit❌❌:   std::vector<std::vector<int>> fragsMolBondMapping;
    // RDKit❌❌:   std::vector<std::string> vfragsmi(mols.size());
    // RDKit❌❌:   std::vector<std::vector<RDKit::UINT>> allAtomOrdering;
    // RDKit❌❌:   std::vector<std::vector<RDKit::UINT>> allBondOrdering;
    // RDKit❌❌:   for (unsigned fragIdx = 0; fragIdx < mols.size(); fragIdx++) {
    // RDKit❌❌:     ROMol *tmol = mols[fragIdx].get();
    // RDKit❌❌:     std::vector<int> atomMapNums(tmol->getNumAtoms(), 0);
    // RDKit❌❌:     for (auto atom : tmol->atoms()) {
    // RDKit❌❌:       atom->updatePropertyCache(false);
    // RDKit❌❌:     }
    // RDKit❌❌:     if (params.doIsomericSmiles) {
    // RDKit❌❌:       tmol->setProp(common_properties::_doIsoSmiles, 1);
    // RDKit❌❌:       if (!tmol->hasProp(common_properties::_StereochemDone)) {
    // RDKit❌❌:         MolOps::assignStereochemistry(*tmol, params.cleanStereo);
    // RDKit❌❌:       }
    // RDKit❌❌:     }
    // RDKit❌❌:     if (params.canonical) {
    // RDKit❌❌:       Canon::rankMolAtoms(*tmol, ranks, breakTies, includeChirality,
    // RDKit❌❌:                           includeIsotopes, includeAtomMaps,
    // RDKit❌❌:                           includeChiralPresence, includeStereoGroups,
    // RDKit❌❌:                           useNonStereoRanks);
    // RDKit❌❌:     } else {
    // RDKit❗✔️:       std::iota(ranks.begin(), ranks.end(), 0);
    // RDKit❌❌:     }
    // RDKit❌❌:     subSmi = SmilesWrite::FragmentSmilesConstruct(
    // RDKit❌❌:         *tmol, nextAtomIdx, colors, ranks, params, atomOrdering, bondOrdering);
    // RDKit❌❌:   }
    // RDKit❌❌:   if (params.canonical) {
    // RDKit❌❌:     std::sort(tmp.begin(), tmp.end());
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     for (unsigned i = 0; i < vfragsmi.size(); ++i) {
    // RDKit❗✔️:       result += vfragsmi[i];
    // RDKit❗✔️:       if (i < vfragsmi.size() - 1) {
    // RDKit❗✔️:         result += ".";
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❌❌:   mol.setProp(common_properties::_smilesAtomOutputOrder, flattenedAtomOrdering,
    // RDKit❌❌:               true);
    // RDKit❌❌:   mol.setProp(common_properties::_smilesBondOutputOrder, flattenedBondOrdering,
    // RDKit❌❌:               true);
    // RDKit❗✔️:   return result;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles
    validate_rooted_atom(molecule, params)?;
    if molecule.num_atoms() == 0 {
        return Ok(String::new());
    }

    let mut context = SmilesWriteContext::default();
    let mut fragment_results = Vec::new();
    let mut working_params = params.clone();

    match mode {
        SmilesOutputMode::PlainSmiles => prepare_plain_smiles_molecule(molecule, &working_params)?,
        SmilesOutputMode::CxSmiles {
            fields,
            restore_bond_dirs,
            include_stereo_groups,
        } => {
            prepare_cx_smiles_molecule(
                molecule,
                &mut working_params,
                fields,
                restore_bond_dirs,
                include_stereo_groups,
            )?;
        }
    }

    let fragment_plans = collect_fragment_write_plans(molecule, &working_params)?;
    for plan in &fragment_plans {
        fragment_results.push(write_fragment_smiles(
            molecule,
            plan,
            &working_params,
            mode,
            &mut context,
        )?);
    }

    let mut result = assemble_fragment_smiles(fragment_results, &working_params, &mut context)?;
    if let SmilesOutputMode::CxSmiles { fields, .. } = mode {
        let cx_extension = get_cx_extensions(molecule, fields)?;
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
    // RDKit❌❌:     // update property cache
    // RDKit❌❌:     std::vector<int> atomMapNums(tmol->getNumAtoms(), 0);
    // RDKit❌❌:     for (auto atom : tmol->atoms()) {
    // RDKit❌❌:       if (params.ignoreAtomMapNumbers) {
    // RDKit❌❌:         atomMapNums[atom->getIdx()] = atom->getAtomMapNum();
    // RDKit❌❌:         atom->setAtomMapNum(0);
    // RDKit❌❌:       }
    // RDKit❌❌:       atom->updatePropertyCache(false);
    // RDKit❌❌:     }
    // RDKit❌❌:
    // RDKit❌❌:     // clean up the chirality on any atom that is marked as chiral,
    // RDKit❌❌:     // but that should not be:
    // RDKit❌❌:     if (params.doIsomericSmiles) {
    // RDKit❌❌:       tmol->setProp(common_properties::_doIsoSmiles, 1);
    // RDKit❌❌:
    // RDKit❌❌:       if (!tmol->hasProp(common_properties::_StereochemDone)) {
    // RDKit❌❌:         MolOps::assignStereochemistry(*tmol, params.cleanStereo);
    // RDKit❌❌:       }
    // RDKit❌❌:     }
    // RDKit❌❌:     if (!doingCXSmiles || !includeStereoGroups) {
    // RDKit❌❌:       std::vector<StereoGroup> noStereoGroups;
    // RDKit❌❌:       tmol->setStereoGroups(noStereoGroups);
    // RDKit❌❌:     }
    // RDKit❌❌:     if (!doingCXSmiles) {
    // RDKit❌❌:       for (auto bond : tmol->bonds()) {
    // RDKit❌❌:         if (bond->getBondDir() == Bond::BondDir::UNKNOWN ||
    // RDKit❌❌:             bond->getBondDir() == Bond::BondDir::EITHERDOUBLE) {
    // RDKit❌❌:           bond->setBondDir(Bond::BondDir::NONE);
    // RDKit❌❌:         }
    // RDKit❌❌:         if (bond->getStereo() == Bond::BondStereo::STEREOANY) {
    // RDKit❌❌:           bond->setStereo(Bond::BondStereo::STEREONONE);
    // RDKit❌❌:         }
    // RDKit❌❌:       }
    // RDKit❌❌:     }
    // RDKit❌❌:     if (doingCXSmiles || !params.includeDativeBonds) {
    // RDKit❌❌:       for (auto bond : tmol->bonds()) {
    // RDKit❌❌:         if (bond->getBondType() == Bond::DATIVE) {
    // RDKit❌❌:           bond->setBondType(Bond::SINGLE);
    // RDKit❌❌:           bond->getBeginAtom()->calcExplicitValence(false);
    // RDKit❌❌:         }
    // RDKit❌❌:       }
    // RDKit❌❌:     }
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles fragment preparation section
    if is_minimal_plain_smiles_path(params) {
        return validate_minimal_plain_smiles_molecule(molecule);
    }
    update_property_cache_for_smiles(molecule)?;
    if params.do_isomeric_smiles {
        assign_stereochemistry_for_smiles(molecule, params.clean_stereo)?;
    }
    if params.do_kekule {
        kekulize_for_smiles(molecule)?;
    }
    if params.do_random {
        return unsupported_stage(SmilesPlanStage::LongTermRandomSmiles);
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
    if params.do_kekule {
        kekulize_for_smiles(molecule)?;
        params.do_kekule = false;
    }
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
    let ranks = rank_fragment_atoms_for_smiles(molecule, plan, params, mode)?;
    let start_atom = choose_fragment_start_atom(plan, &ranks, params)?;
    fragment_smiles_construct(molecule, plan, start_atom, &ranks, params, context)
}

fn fragment_smiles_construct(
    molecule: &Molecule,
    plan: &FragmentWritePlan,
    start_atom: AtomId,
    ranks: &[usize],
    params: &SmilesWriteParams,
    context: &mut SmilesWriteContext,
) -> Result<FragmentWriteResult, SmilesWriteError> {
    if params.do_kekule {
        kekulize_fragment_for_smiles(molecule, plan)?;
    }
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
    let _ = (molecule, plan, params, mode);
    unsupported_stage(SmilesPlanStage::LongTermCanonicalRanking)
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
    let (idx, _) = ranks
        .iter()
        .enumerate()
        .min_by_key(|(_, rank)| **rank)
        .ok_or(SmilesWriteError::NotImplemented)?;
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
    // RDKit❌❌:     subSmi = SmilesWrite::FragmentSmilesConstruct(
    // RDKit❌❌:         *tmol, nextAtomIdx, colors, ranks, params, atomOrdering, bondOrdering);
    // RDKit❌❌: Canon::canonicalizeFragment(mol, atomIdx, colors, ranks, molStack,
    // RDKit❌❌:                           atomsInPlay, bondsInPlay, bondSymbols,
    // RDKit❌❌:                           params.doIsomericSmiles, params.doRandom);
    // END RDKIT CPP FUNCTION Canon::canonicalizeFragment call site
    let _ = (ranks, params);
    build_minimal_noncanonical_stack(molecule, plan, start_atom)
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
    // RDKit❌❌:           res << (*atomSymbols)[mSE.obj.atom->getIdx()];
    // RDKit❌❌:         }
    // RDKit❗✔️:         atomOrdering.push_back(mSE.obj.atom->getIdx());
    // RDKit❗✔️:         break;
    // RDKit❗✔️:       case Canon::MOL_STACK_BOND:
    // RDKit❗✔️:         bond = mSE.obj.bond;
    // RDKit❗✔️:         if (!bondSymbols) {
    // RDKit❗✔️:           res << GetBondSmiles(bond, params, mSE.number);
    // RDKit❗✔️:         } else {
    // RDKit❌❌:           res << (*bondSymbols)[bond->getIdx()];
    // RDKit❌❌:         }
    // RDKit❗✔️:         bondOrdering.push_back(bond->getIdx());
    // RDKit❗✔️:         break;
    // RDKit❌❌:       case Canon::MOL_STACK_RING:
    // RDKit❌❌:       case Canon::MOL_STACK_BRANCH_OPEN:
    // RDKit❌❌:       case Canon::MOL_STACK_BRANCH_CLOSE:
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

pub fn mol_fragment_to_smiles(
    molecule: &Molecule,
    params: &SmilesWriteParams,
    atoms_to_use: &[usize],
    bonds_to_use: Option<&[usize]>,
    atom_symbols: Option<&[String]>,
    bond_symbols: Option<&[String]>,
) -> Result<String, SmilesWriteError> {
    let _ = (
        params,
        atoms_to_use,
        bonds_to_use,
        atom_symbols,
        bond_symbols,
    );
    validate_fragment_api_inputs(molecule, atoms_to_use, bonds_to_use)?;
    unsupported_stage(SmilesPlanStage::LongTermFragmentApi)
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
    let _ = fields;
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
    // RDKit❌❌:
    // RDKit❗✔️:   std::string symb;
    // RDKit❌❌:   bool hasCustomSymbol =
    // RDKit❌❌:       atom->getPropIfPresent(common_properties::smilesSymbol, symb);
    // RDKit❗✔️:   if (!hasCustomSymbol) {
    // RDKit❗✔️:     symb = PeriodicTable::getTable()->getElementSymbol(num);
    // RDKit❗✔️:   }
    // RDKit❌❌:
    // RDKit❌❌:   // check for atomic stereochemistry
    // RDKit❌❌:   std::string atString;
    // RDKit❌❌:   if (params.doIsomericSmiles) {
    // RDKit❌❌:     if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit❌❌:         !atom->hasProp(common_properties::_brokenChirality)) {
    // RDKit❌❌:       atString = getAtomChiralityInfo(atom);
    // RDKit❌❌:     }
    // RDKit❌❌:   }
    // RDKit❗✔️:   bool needsBracket = true;
    // RDKit❗✔️:   if (!hasCustomSymbol && !params.allHsExplicit) {
    // RDKit❗✔️:     needsBracket = atomNeedsBracket(atom, atString, params);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (needsBracket) {
    // RDKit❗✔️:     res += "[";
    // RDKit❗✔️:   }
    // RDKit❌❌:
    // RDKit❗✔️:   if (isotope && params.doIsomericSmiles) {
    // RDKit❗✔️:     res += std::to_string(isotope);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (!params.doKekule && atom->getIsAromatic() && symb[0] >= 'A' &&
    // RDKit❗✔️:       symb[0] <= 'Z') {
    // RDKit❌❌:   }
    // RDKit❗✔️:   res += symb;
    // RDKit❌❌:   res += atString;
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
    // RDKit❌❌:   std::string label;
    // RDKit❌❌:   if (atom->getPropIfPresent(common_properties::_supplementalSmilesLabel,
    // RDKit❌❌:                              label)) {
    // RDKit❌❌:     res += label;
    // RDKit❌❌:   }
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
    if atom.query().is_some()
        || atom.radical_electrons() != 0
        || atom.chiral_tag() != ChiralTag::Unspecified
        || atom.is_aromatic()
    {
        return unsupported_stage(SmilesPlanStage::ShortTermAtomWriter);
    }
    let needs_bracket =
        params.all_hydrogens_explicit || atom_needs_bracket(molecule, atom_id, params)?;
    let symbol = element_symbol(atom.atomic_number())?;
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
    match _bond_order {
        BondOrder::Single => Ok(""),
        BondOrder::Double => Ok("="),
        BondOrder::Triple => Ok("#"),
        BondOrder::Quadruple => Ok("$"),
        BondOrder::Dative => Ok("->"),
        _ => unsupported_stage(SmilesPlanStage::ShortTermBondWriter),
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
    // RDKit❌❌:   std::string res = "";
    // RDKit❌❌:   bool aromatic = false;
    // RDKit❌❌:   if (!params.doKekule && (bond->getBondType() == Bond::SINGLE ||
    // RDKit❌❌:                            bond->getBondType() == Bond::DOUBLE ||
    // RDKit❌❌:                            bond->getBondType() == Bond::AROMATIC)) {
    // RDKit❌❌:   }
    // RDKit❌❌:   Bond::BondDir dir = bond->getBondDir();
    // RDKit❌❌:   bond->clearProp(common_properties::_TraversalRingClosureBond);
    // RDKit❗✔️:   switch (bond->getBondType()) {
    // RDKit❗✔️:     case Bond::SINGLE:
    // RDKit❗✔️:       if (dir != Bond::NONE && dir != Bond::UNKNOWN) {
    // RDKit❌❌:       } else {
    // RDKit❗✔️:         if (params.allBondsExplicit) {
    // RDKit❗✔️:           res = "-";
    // RDKit❌❌:         } else if (aromatic && !bond->getIsAromatic()) {
    // RDKit❌❌:           res = "-";
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
    // RDKit❌❌:     case Bond::AROMATIC:
    // RDKit❗✔️:     case Bond::DATIVE:
    // RDKit❗✔️:       if (atomToLeftIdx >= 0 &&
    // RDKit❗✔️:           bond->getBeginAtomIdx() == static_cast<unsigned int>(atomToLeftIdx)) {
    // RDKit❗✔️:         res = "->";
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         res = "<-";
    // RDKit❗✔️:       }
    // RDKit❗✔️:       break;
    // RDKit❌❌:     default:
    // RDKit❌❌:       res = "~";
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION GetBondSmiles
    validate_bond_index(molecule, bond)?;
    if let Some(atom) = atom_to_left {
        validate_atom_index(molecule, atom)?;
    }
    let bond = &molecule.bonds()[bond];
    if bond.direction() != BondDirection::None
        || bond.stereo() != BondStereo::None
        || bond.query().is_some()
    {
        return unsupported_stage(SmilesPlanStage::ShortTermBondWriter);
    }
    if bond.is_aromatic() || bond.order() == BondOrder::Aromatic {
        return unsupported_stage(SmilesPlanStage::ShortTermBondWriter);
    }
    if bond.order() == BondOrder::Dative {
        let left = atom_to_left.unwrap_or_else(|| bond.begin().index());
        return if bond.begin().index() == left {
            Ok("->".to_string())
        } else {
            Ok("<-".to_string())
        };
    }
    let symbol = get_bond_smiles(bond.order())?;
    if symbol.is_empty() && params.all_bonds_explicit {
        Ok("-".to_string())
    } else {
        Ok(symbol.to_string())
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

pub fn get_cx_extensions(
    molecule: &Molecule,
    fields: CxSmilesFields,
) -> Result<String, SmilesWriteError> {
    let _ = (molecule, fields);
    unsupported_stage(SmilesPlanStage::LongTermCxExtensions)
}

fn get_atom_chirality_info(molecule: &Molecule, atom: AtomId) -> Result<String, SmilesWriteError> {
    let _ = (molecule, atom);
    unsupported_stage(SmilesPlanStage::LongTermIsomericStereo)
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
    // RDKit❌❌:   const INT_VECT &defaultVs = PeriodicTable::getTable()->getValenceList(num);
    // RDKit❌❌:   int totalValence = atom->getTotalValence();
    // RDKit❌❌:   bool nonStandard = false;
    // RDKit❌❌:   if (atom->getNumRadicalElectrons()) {
    // RDKit❌❌:     nonStandard = true;
    // RDKit❌❌:   } else if ((num == 7 || num == 15) && atom->getIsAromatic() &&
    // RDKit❌❌:              atom->getNumExplicitHs()) {
    // RDKit❌❌:     nonStandard = true;
    // RDKit❌❌:   } else {
    // RDKit❌❌:     nonStandard = (totalValence != defaultVs.front() && atom->getTotalNumHs());
    // RDKit❌❌:   }
    // RDKit❌❌:   if (nonStandard) {
    // RDKit❌❌:     return true;
    // RDKit❌❌:   }
    // RDKit❌❌:
    // RDKit❌❌:   // check for bonds to a metal
    // RDKit❗✔️:   return false;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION atomNeedsBracket
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
    Ok(false)
}

fn update_property_cache_for_smiles(_molecule: &Molecule) -> Result<(), SmilesWriteError> {
    unsupported_stage(SmilesPlanStage::ShortTermAtomWriter)
}

fn assign_stereochemistry_for_smiles(
    _molecule: &Molecule,
    _clean_stereo: bool,
) -> Result<(), SmilesWriteError> {
    unsupported_stage(SmilesPlanStage::LongTermIsomericStereo)
}

fn canonicalize_enhanced_stereo_for_smiles(_molecule: &Molecule) -> Result<(), SmilesWriteError> {
    unsupported_stage(SmilesPlanStage::LongTermIsomericStereo)
}

fn cleanup_stereo_groups_for_cx_smiles(_molecule: &Molecule) -> Result<(), SmilesWriteError> {
    unsupported_stage(SmilesPlanStage::LongTermCxExtensions)
}

fn kekulize_for_smiles(_molecule: &Molecule) -> Result<(), SmilesWriteError> {
    unsupported_stage(SmilesPlanStage::LongTermKekule)
}

fn kekulize_fragment_for_smiles(
    _molecule: &Molecule,
    _plan: &FragmentWritePlan,
) -> Result<(), SmilesWriteError> {
    unsupported_stage(SmilesPlanStage::LongTermKekule)
}

fn normalize_dative_bonds_for_plain_smiles(_molecule: &Molecule) -> Result<(), SmilesWriteError> {
    unsupported_stage(SmilesPlanStage::ShortTermBondWriter)
}

fn normalize_dative_bonds_for_cx_smiles(_molecule: &Molecule) -> Result<(), SmilesWriteError> {
    unsupported_stage(SmilesPlanStage::LongTermCxExtensions)
}

fn normalize_hydrogen_bonds_for_cx_smiles(_molecule: &Molecule) -> Result<(), SmilesWriteError> {
    unsupported_stage(SmilesPlanStage::LongTermCxExtensions)
}

fn apply_cx_bond_direction_policy(
    _molecule: &Molecule,
    _restore_bond_dirs: RestoreBondDirOption,
) -> Result<(), SmilesWriteError> {
    unsupported_stage(SmilesPlanStage::LongTermCxExtensions)
}

fn remove_plain_smiles_only_cx_state(_molecule: &Molecule) -> Result<(), SmilesWriteError> {
    unsupported_stage(SmilesPlanStage::LongTermCxExtensions)
}

fn validate_cx_extension_plan(_fields: CxSmilesFields) -> Result<(), SmilesWriteError> {
    unsupported_stage(SmilesPlanStage::LongTermCxExtensions)
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
            || atom.is_aromatic()
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
            || bond.is_aromatic()
            || bond.query().is_some()
            || !matches!(
                bond.order(),
                BondOrder::Single
                    | BondOrder::Double
                    | BondOrder::Triple
                    | BondOrder::Quadruple
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
    ) -> Vec<(BondId, AtomId)> {
        let mut incident = molecule
            .bonds()
            .iter()
            .filter(|bond| plan.bonds.contains(&bond.id()))
            .filter_map(|bond| bond_other_atom(bond, atom).map(|other| (bond.id(), other)))
            .collect::<Vec<_>>();
        incident.sort_by_key(|(bond, other)| (other.index(), bond.index()));
        incident
    }

    fn find_cycles_and_tree(
        molecule: &Molecule,
        plan: &FragmentWritePlan,
        atom: AtomId,
        parent_bond: Option<BondId>,
        colors: &mut [Color],
        tree_children: &mut [Vec<(BondId, AtomId)>],
        ring_closures: &mut Vec<RingClosure>,
        visited_bonds: &mut [bool],
    ) {
        colors[atom.index()] = Color::Grey;
        for (bond, other) in sorted_incident_bonds(molecule, plan, atom) {
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
        &mut colors,
        &mut tree_children,
        &mut ring_closures,
        &mut visited_bonds,
    );
    build_stack_from_tree(start_atom, &tree_children, &ring_closures, &mut stack);

    if plan
        .atoms
        .iter()
        .any(|atom| colors[atom.index()] == Color::White)
        || plan.bonds.iter().any(|bond| !visited_bonds[bond.index()])
    {
        return unsupported_stage(SmilesPlanStage::ShortTermFragmentTraversal);
    }
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

    let digit = (1..)
        .find(|candidate| {
            !context
                .ring_closure_digits
                .values()
                .any(|digit| digit == candidate)
        })
        .ok_or(SmilesWriteError::NotImplemented)?;
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
        _ => unsupported_stage(SmilesPlanStage::ShortTermAtomWriter),
    }
}

fn assemble_fragment_smiles(
    fragment_results: Vec<FragmentWriteResult>,
    params: &SmilesWriteParams,
    context: &mut SmilesWriteContext,
) -> Result<String, SmilesWriteError> {
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles fragment assembly section
    // RDKit❌❌:   if (params.canonical) {
    // RDKit❌❌:     std::sort(tmp.begin(), tmp.end());
    // RDKit❌❌:   } else {  // Not canonical
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
    // RDKit❌❌:   mol.setProp(common_properties::_smilesAtomOutputOrder, flattenedAtomOrdering,
    // RDKit❌❌:               true);
    // RDKit❌❌:   mol.setProp(common_properties::_smilesBondOutputOrder, flattenedBondOrdering,
    // RDKit❌❌:               true);
    // RDKit❗✔️:   return result;
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles fragment assembly section
    if params.canonical {
        return unsupported_stage(SmilesPlanStage::LongTermCanonicalRanking);
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
    fn molecule_to_smiles_enters_smiles_writer_and_fails_closed() {
        let error = ethane().to_smiles(true).unwrap_err();

        assert!(matches!(error, SmilesWriteError::UnsupportedFeature(_)));
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
        let mut params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            ..Default::default()
        };

        assert_eq!(mol_to_smiles(&molecule, &params).unwrap(), "CC");

        params.do_kekule = true;
        assert!(matches!(
            mol_to_smiles(&molecule, &params),
            Err(SmilesWriteError::UnsupportedFeature(_))
        ));

        params.do_kekule = false;
        params.do_random = true;
        assert!(matches!(
            mol_to_smiles(&molecule, &params),
            Err(SmilesWriteError::UnsupportedFeature(_))
        ));

        assert!(matches!(
            mol_to_cx_smiles(
                &molecule,
                &SmilesWriteParams::default(),
                CxSmilesFields::ALL,
                RestoreBondDirOption::Clear,
            ),
            Err(SmilesWriteError::UnsupportedFeature(_))
        ));
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
    fn noncanonical_nonisomeric_plain_smiles_rejects_aromatic_query_zero_and_radical_bonds() {
        let params = SmilesWriteParams {
            do_isomeric_smiles: false,
            canonical: false,
            clean_stereo: false,
            ..Default::default()
        };
        let aromatic = Molecule::from_smiles_with_sanitize("c1ccccc1", false).unwrap();
        let zero = Molecule::from_smiles_with_sanitize("CC~CC |Z:1|", false).unwrap();
        let radical = Molecule::from_smiles_with_sanitize("CCC |^1:0|", false).unwrap();

        assert!(matches!(
            mol_to_smiles(&aromatic, &params),
            Err(SmilesWriteError::UnsupportedFeature(_))
        ));
        assert!(matches!(
            mol_to_smiles(&zero, &params),
            Err(SmilesWriteError::UnsupportedFeature(_))
        ));
        assert!(matches!(
            mol_to_smiles(&radical, &params),
            Err(SmilesWriteError::UnsupportedFeature(_))
        ));
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
        assert!(matches!(
            mol_fragment_to_smiles(
                &molecule,
                &SmilesWriteParams::default(),
                &[0, 1],
                None,
                None,
                None
            ),
            Err(SmilesWriteError::UnsupportedFeature(_))
        ));
        assert!(matches!(
            mol_fragment_to_cx_smiles(
                &molecule,
                &SmilesWriteParams::default(),
                &[0, 1],
                None,
                None,
                None,
                CxSmilesFields::ALL,
            ),
            Err(SmilesWriteError::UnsupportedFeature(_))
        ));
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
}
