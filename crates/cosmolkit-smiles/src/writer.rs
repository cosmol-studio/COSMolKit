use std::collections::BTreeMap;

use cosmolkit_core::{KekulizeParams, ValenceAssignment, ValenceModel, fast_find_rings_from_parts};
use cosmolkit_model::{
    Atom, AtomId, Bond, BondId, PropertyText, PropertyValue, StereoGroup, TopologyBlock,
    set_stereo_group_write_id, stereo_group_write_id,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo, ChiralTag};

use crate::fragment::{
    FragmentSelectionMasks, FragmentWriteInputError, PreparedFragmentStereo,
    build_fragment_selection_masks, prepare_fragment_stereo, rank_prepared_fragment,
    validate_fragment_write_inputs,
};
use crate::{
    CxSmilesWriteParams, SmilesParseError, SmilesRecord, canonical_rank,
    cx_writer::write_cx_extensions, stereo,
};

mod direction;

const MAX_NATOMS: i64 = 5000;
const MAX_BONDTYPE: i64 = 32;
const MAX_CYCLES: usize = 1024;

fn source_string_property(
    value: &PropertyValue,
    _name: &str,
) -> Result<PropertyText, SmilesParseError> {
    // RDKit❗✔️: bool getValIfPresent(const std::string_view what, std::string &res) const {
    // RDKit❗✔️:     for (const auto &i : _data) {
    // RDKit❗✔️:       if (i.key == what) {
    // RDKit❗✔️:         rdvalue_tostring(i.val, res);
    // RDKit❗✔️:         return true;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     return false;
    // RDKit❗✔️:   }
    cosmolkit_core::property_value_to_string(value).map_err(SmilesParseError::WriterProperty)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[doc(hidden)]
pub enum AtomColor {
    White,
    Grey,
    Black,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[doc(hidden)]
pub enum MolStackElem {
    Atom(usize),
    Bond { bond: BondId, atom_to_left: usize },
    Ring(i32),
    BranchOpen(i32),
    BranchClose(i32),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Possible {
    rank: i64,
    atom: usize,
    bond: BondId,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
struct ChiralAdjustment {
    chiral_tag_override: Option<ChiralTag>,
    invert_tetrahedral: bool,
    nontetrahedral_permutation: Option<u32>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmilesWriteOutput {
    pub text: PropertyText,
    pub atom_order: Vec<AtomId>,
    pub bond_order: Vec<BondId>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct FragmentWriteOutput {
    text: PropertyText,
    atom_order: Vec<AtomId>,
    bond_order: Vec<BondId>,
}

#[derive(Debug)]
struct WriterStereoFragment {
    topology: TopologyBlock,
    source_atoms: Vec<AtomId>,
    source_bonds: Vec<BondId>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SmilesWriteParams {
    pub isomeric_smiles: bool,
    pub kekule: bool,
    pub canonical: bool,
    pub clean_stereo: bool,
    pub rooted_at_atom: Option<AtomId>,
    pub all_bonds_explicit: bool,
    pub all_hydrogens_explicit: bool,
    pub include_dative_bonds: bool,
    pub ignore_atom_map_numbers: bool,
}

impl Default for SmilesWriteParams {
    fn default() -> Self {
        Self {
            isomeric_smiles: true,
            kekule: false,
            canonical: true,
            clean_stereo: true,
            rooted_at_atom: None,
            all_bonds_explicit: false,
            all_hydrogens_explicit: false,
            include_dative_bonds: true,
            ignore_atom_map_numbers: false,
        }
    }
}

/// Source options accepted by RDKit's random-SMILES vector writer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RandomSmilesWriteParams {
    pub isomeric_smiles: bool,
    pub kekule: bool,
    pub all_bonds_explicit: bool,
    pub all_hydrogens_explicit: bool,
}

impl Default for RandomSmilesWriteParams {
    fn default() -> Self {
        Self {
            isomeric_smiles: true,
            kekule: false,
            all_bonds_explicit: false,
            all_hydrogens_explicit: false,
        }
    }
}

/// Writes canonical SMILES from detached values using RDKit-compatible atom
/// ranking and traversal.
pub fn write_smiles<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
) -> Result<PropertyText, SmilesParseError> {
    let record = record.into();
    write_smiles_with_params(record, &SmilesWriteParams::default())
}

/// Writes SMILES from detached values with explicit canonicalization policy.
pub fn write_smiles_with_params<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
    params: &SmilesWriteParams,
) -> Result<PropertyText, SmilesParseError> {
    let record = record.into();
    write_smiles_output(record, params, false).map(|output| output.text)
}

/// Writes SMILES with source random-traversal behavior when `do_random` is
/// true, borrowing the shared process stream for the complete write.
pub fn write_smiles_with_random<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
    params: &SmilesWriteParams,
    do_random: bool,
) -> Result<PropertyText, SmilesParseError> {
    let record = record.into();
    // BEGIN RDKIT CPP FUNCTION MolToSmiles params overload
    // RDKit❗❌: std::string MolToSmiles(const ROMol &mol, const SmilesWriteParams &params) {
    // RDKit❗❌:   bool doingCXSmiles = false;
    // RDKit❗❌:   return SmilesWrite::detail::MolToSmiles(mol, params, doingCXSmiles);
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION MolToSmiles params overload
    // Behavior: the domain `do_random` argument selects the source random
    // traversal path for this whole-molecule call; one seed-zero core borrow
    // preserves the process stream for every root, cycle, and stack draw.
    // Complexity: the closure acquires one O(1) mutex per writer call, with no
    // lock or allocation added to individual random draws.
    if do_random {
        cosmolkit_core::with_rdkit_random_generator(0, |random_stream| {
            write_smiles_output_with_random_stream(record, params, false, Some(random_stream))
        })
        .map(|output| output.text)
    } else {
        write_smiles_with_params(record, params)
    }
}

/// Write an ordered vector of source-compatible random SMILES.
///
/// A positive `random_seed` reseeds the shared stream once before producing
/// any rows. Zero preserves its current state. Duplicate strings are retained.
pub fn write_random_smiles_vector<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
    num_smiles: u32,
    random_seed: u32,
    params: &RandomSmilesWriteParams,
) -> Result<Vec<PropertyText>, SmilesParseError> {
    let record = record.into();
    // BEGIN RDKIT CPP FUNCTION MolToRandomSmilesVect
    // RDKit❗❌: std::vector<std::string> MolToRandomSmilesVect(
    // RDKit❗❌:     const ROMol &mol, unsigned int numSmiles, unsigned int randomSeed,
    // RDKit❗❌:     bool doIsomericSmiles, bool doKekule, bool allBondsExplicit,
    // RDKit❗❌:     bool allHsExplicit) {
    // RDKit❗❌:   if (randomSeed > 0) {
    // RDKit❗❌:     getRandomGenerator(rdcast<int>(randomSeed));
    // RDKit❗❌:   }
    // RDKit❗❌:   std::vector<std::string> res;
    // RDKit❗❌:   res.reserve(numSmiles);
    // RDKit❗❌:   for (unsigned int i = 0; i < numSmiles; ++i) {
    // RDKit❗❌:     bool canonical = false;
    // RDKit❗❌:     int rootedAtAtom = -1;
    // RDKit❗❌:     bool doRandom = true;
    // RDKit❗❌:     res.push_back(MolToSmiles(mol, doIsomericSmiles, doKekule, rootedAtAtom,
    // RDKit❗❌:                               canonical, allBondsExplicit, allHsExplicit,
    // RDKit❗❌:                               doRandom));
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: };
    // END RDKIT CPP FUNCTION MolToRandomSmilesVect
    // BEGIN RDKIT CPP FUNCTION MolToSmiles boolean-overload parameter defaults
    // RDKit❗❌: inline std::string MolToSmiles(const ROMol &mol, bool doIsomericSmiles = true,
    // RDKit❗❌:                                bool doKekule = false, int rootedAtAtom = -1,
    // RDKit❗❌:                                bool canonical = true,
    // RDKit❗❌:                                bool allBondsExplicit = false,
    // RDKit❗❌:                                bool allHsExplicit = false,
    // RDKit❗❌:                                bool doRandom = false,
    // RDKit❗❌:                                bool ignoreAtomMapNumbers = false) {
    // RDKit❗❌:   SmilesWriteParams ps;
    // RDKit❗❌:   ps.doIsomericSmiles = doIsomericSmiles;
    // RDKit❗❌:   ps.doKekule = doKekule;
    // RDKit❗❌:   ps.rootedAtAtom = rootedAtAtom;
    // RDKit❗❌:   ps.canonical = canonical;
    // RDKit❗❌:   ps.allBondsExplicit = allBondsExplicit;
    // RDKit❗❌:   ps.allHsExplicit = allHsExplicit;
    // RDKit❗❌:   ps.doRandom = doRandom;
    // RDKit❗❌:   ps.ignoreAtomMapNumbers = ignoreAtomMapNumbers;
    // RDKit❗❌:   return MolToSmiles(mol, ps);
    // RDKit❗❌: };
    // END RDKIT CPP FUNCTION MolToSmiles boolean-overload parameter defaults
    // BEGIN RDKIT CPP FUNCTION RDGeneral/Invariant.h rdcast configuration
    // RDKit❗❌: #ifdef RDDEBUG
    // RDKit❗❌: #define rdcast boost::numeric_cast
    // RDKit❗❌: #else
    // RDKit❗❌: #define rdcast static_cast
    // RDKit❗❌: #endif
    // END RDKIT CPP FUNCTION RDGeneral/Invariant.h rdcast configuration
    // Behavior: Rust's u32-to-i32 cast preserves the low 32-bit pattern, which
    // matches the loaded release reference's static_cast behavior: seeds with
    // the high bit set become nonpositive and do not reseed. The alternate
    // RDDEBUG checked-cast configuration is not the pinned loaded oracle.
    // Complexity: one cast and one outer O(1) core lock; the shared mutex adds
    // synchronization versus RDKit's global stream, while each draw is lock-free.
    let seed = random_seed as i32;
    cosmolkit_core::with_rdkit_random_generator(seed, |random_stream| {
        let writer_params = SmilesWriteParams {
            isomeric_smiles: params.isomeric_smiles,
            kekule: params.kekule,
            canonical: false,
            clean_stereo: true,
            rooted_at_atom: None,
            all_bonds_explicit: params.all_bonds_explicit,
            all_hydrogens_explicit: params.all_hydrogens_explicit,
            include_dative_bonds: true,
            ignore_atom_map_numbers: false,
        };
        let mut results = Vec::with_capacity(num_smiles as usize);
        for _ in 0..num_smiles {
            let output = write_smiles_output_with_random_stream(
                record,
                &writer_params,
                false,
                Some(&mut *random_stream),
            )?;
            results.push(output.text);
        }
        Ok(results)
    })
}

pub(crate) fn write_smiles_for_cx<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
    params: &SmilesWriteParams,
) -> Result<SmilesWriteOutput, SmilesParseError> {
    let record = record.into();
    write_smiles_output(record, params, true)
}

fn canonicalize_enhanced_stereo(
    topology: &mut TopologyBlock,
    ranks: &[i64],
) -> Result<BTreeMap<usize, usize>, SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION Canon.cpp::canonicalizeEnhancedStereo
    // RDKit❗✔️: void canonicalizeEnhancedStereo(ROMol &mol,
    // RDKit❗✔️:                                 const std::vector<unsigned int> *atomRanks) {
    // RDKit❗✔️:   const auto &sgs = mol.getStereoGroups();
    // RDKit❗✔️:   if (sgs.empty()) {
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   std::vector<unsigned int> lranks;
    // RDKit❗✔️:   if (!atomRanks) {
    // RDKit❌❌:     bool breakTies = true;
    // RDKit❌❌:     rankMolAtoms(mol, lranks, breakTies);
    // RDKit❌❌:     atomRanks = &lranks;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // one thing that makes this all easier is that the stereogroups are
    // RDKit❗✔️:   // independent of each other
    // RDKit❗✔️:   std::vector<StereoGroup> newSgs;
    // RDKit❗✔️:   for (auto &sg : sgs) {
    // RDKit❗✔️:     // we don't do anything to ABS groups
    // RDKit❗✔️:     if (sg.getGroupType() == StereoGroupType::STEREO_ABSOLUTE) {
    // RDKit❗✔️:       newSgs.push_back(sg);
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     // sort the atoms by rank:
    // RDKit❗✔️:     auto getAtomRank = [&atomRanks](const Atom *at1, const Atom *at2) {
    // RDKit❗✔️:       return atomRanks->at(at1->getIdx()) < atomRanks->at(at2->getIdx());
    // RDKit❗✔️:     };
    // RDKit❗✔️:     auto sgAtoms = sg.getAtoms();
    // RDKit❗✔️:     std::sort(sgAtoms.begin(), sgAtoms.end(), getAtomRank);
    // RDKit❗✔️:
    // RDKit❗✔️:     // sort the bonds by atom rank:
    // RDKit❗✔️:     auto getBondRank = [&atomRanks](const Bond *bd1, const Bond *bd2) {
    // RDKit❗✔️:       unsigned int bd1at1 = atomRanks->at(bd1->getBeginAtomIdx());
    // RDKit❗✔️:       unsigned int bd1at2 = atomRanks->at(bd1->getEndAtomIdx());
    // RDKit❗✔️:       unsigned int bd2at1 = atomRanks->at(bd2->getBeginAtomIdx());
    // RDKit❗✔️:       unsigned int bd2at2 = atomRanks->at(bd2->getEndAtomIdx());
    // RDKit❗✔️:       if (bd1at1 < bd1at2) {
    // RDKit❗✔️:         std::swap(bd1at1, bd1at2);
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (bd2at1 < bd2at2) {
    // RDKit❗✔️:         std::swap(bd2at1, bd2at2);
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (bd1at1 != bd2at1) {
    // RDKit❗✔️:         return bd1at1 < bd2at1;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       return bd1at2 < bd2at2;
    // RDKit❗✔️:     };
    // RDKit❗✔️:     auto sgBonds = sg.getBonds();
    // RDKit❗✔️:     std::sort(sgBonds.begin(), sgBonds.end(), getBondRank);
    // RDKit❗✔️:
    // RDKit❗✔️:     // find the reference (lowest-ranked) atom (or lowest-ranked bond)
    // RDKit❗✔️:     Atom::ChiralType foundRefState = Atom::ChiralType::CHI_TETRAHEDRAL_CCW;
    // RDKit❗✔️:     if (sgAtoms.size() > 0) {
    // RDKit❗✔️:       foundRefState = sgAtoms.front()->getChiralTag();
    // RDKit❗✔️:     } else if (sgBonds.size() > 0) {
    // RDKit❗✔️:       if (sgBonds.front()->getStereo() == Bond::BondStereo::STEREOATROPCCW) {
    // RDKit❗✔️:         foundRefState =
    // RDKit❗✔️:             Atom::ChiralType::CHI_TETRAHEDRAL_CCW;  // convert atropisomer CCW
    // RDKit❗✔️:                                                     // to atom CCW
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         foundRefState =
    // RDKit❗✔️:             Atom::ChiralType::CHI_TETRAHEDRAL_CW;  // convert atropisomer CW
    // RDKit❗✔️:                                                     // to atom CW
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     // we will use CCW as the "canonical" state for chirality, so if the
    // RDKit❗✔️:     // referenceAtom is already CCW then we don't need to do anything more
    // RDKit❗✔️:     // with this stereogroup
    // RDKit❗✔️:     auto refState = Atom::ChiralType::CHI_TETRAHEDRAL_CCW;
    // RDKit❗✔️:     if (foundRefState != refState) {
    // RDKit❗✔️:       // we need to flip everyone... so loop over the other atoms and bonds
    // RDKit❗✔️:       // and flip them all:
    // RDKit❗✔️:       for (auto atom : sgAtoms) {
    // RDKit❗✔️:         atom->invertChirality();
    // RDKit❗✔️:       }
    // RDKit❗✔️:       for (auto bond : sgBonds) {
    // RDKit❗✔️:         bond->invertChirality();
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     newSgs.emplace_back(
    // RDKit❗✔️:         StereoGroup(sg.getGroupType(), std::move(sgAtoms), std::move(sgBonds)));
    // RDKit❗✔️:
    // RDKit❗✔️:     // note that we do not forward the Group Ids: this is intentional, so that
    // RDKit❗✔️:     // the Ids are reassigned based on the canonicalized order.
    // RDKit❗✔️:     if (sgAtoms.size() > 0) {
    // RDKit❗✔️:       sgAtoms.front()->setProp("_stereoGroup", newSgs.size() - 1, true);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   mol.setStereoGroups(newSgs);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION Canon.cpp::canonicalizeEnhancedStereo
    // Behavior review: source rank comparators, group order, ABS pass-through,
    // atom-first/bond-only reference selection, full-group inversion and default
    // non-ABS IDs are preserved. The source's transient size_t `_stereoGroup`
    // writes are represented by this sparse usize map and consumed by the later
    // source post-processing branch; no model property narrowing occurs.
    // Complexity review: member vectors are copied once and sorted in O(A log A
    // + B log B); unstable Rust sorts and C++ std::sort have matching asymptotic
    // cost without secondary tie keys. Sparse BTreeMap writes retain O(G) state,
    // corresponding to the source's sparse per-reference-atom properties.

    if topology.stereo_groups.is_empty() {
        return Ok(BTreeMap::new());
    }

    let mut canonical_groups = Vec::with_capacity(topology.stereo_groups.len());
    let mut atom_group_references = BTreeMap::new();
    for group in &topology.stereo_groups {
        if group.kind() == cosmolkit_model::StereoGroupKind::Absolute {
            canonical_groups.push(group.clone());
            continue;
        }

        let mut atoms = group.atoms().to_vec();
        atoms.sort_unstable_by(|first, second| ranks[first.index()].cmp(&ranks[second.index()]));

        let rank_pair = |bond_id: BondId| {
            let bond = &topology.bonds[bond_id.index()];
            let begin_rank = ranks[bond.begin().index()];
            let end_rank = ranks[bond.end().index()];
            if begin_rank < end_rank {
                (end_rank, begin_rank)
            } else {
                (begin_rank, end_rank)
            }
        };
        let mut bonds = group.bonds().to_vec();
        bonds.sort_unstable_by_key(|bond_id| rank_pair(*bond_id));

        let found_reference_state = if let Some(reference_atom) = atoms.first() {
            topology.atoms[reference_atom.index()].chiral_tag()
        } else if let Some(reference_bond) = bonds.first() {
            if topology.bonds[reference_bond.index()].stereo() == BondStereo::AtropCcw {
                ChiralTag::TetrahedralCcw
            } else {
                ChiralTag::TetrahedralCw
            }
        } else {
            ChiralTag::TetrahedralCcw
        };

        if found_reference_state != ChiralTag::TetrahedralCcw {
            for atom_id in &atoms {
                cosmolkit_core::invert_atom_chirality(&mut topology.atoms[atom_id.index()])
                    .map_err(SmilesParseError::WriterStereoOrder)?;
            }
            for bond_id in &bonds {
                cosmolkit_core::invert_bond_chirality(&mut topology.bonds[bond_id.index()])
                    .map_err(SmilesParseError::WriterStereoBond)?;
            }
        }

        let group_index = canonical_groups.len();
        if let Some(reference_atom) = atoms.first() {
            atom_group_references.insert(reference_atom.index(), group_index);
        }
        canonical_groups.push(StereoGroup::new(group.kind(), atoms, bonds)?);
    }
    topology.stereo_groups = canonical_groups;
    Ok(atom_group_references)
}

/// Writes a selected original-index fragment without constructing a renumbered
/// topology. `source_rings` and `existing_valence` are optional detached state
/// supplied only by an owner that can prove they describe this exact topology;
/// a bare `SmilesRecord` supplies neither.
pub fn write_fragment_smiles_output<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
    params: &SmilesWriteParams,
    atoms_to_use: &[AtomId],
    bonds_to_use: Option<&[BondId]>,
    atom_symbols: Option<&[String]>,
    bond_symbols: Option<&[String]>,
    source_rings: Option<&cosmolkit_core::RingInfo>,
    existing_valence: Option<&ValenceAssignment>,
) -> Result<SmilesWriteOutput, FragmentWriteInputError> {
    let record = record.into();
    // BEGIN COMPLETE RDKit .6 CHEM27 SmilesParse::MolFragmentToSmiles
    // RDKit❗❌: std::string MolFragmentToSmiles(const ROMol &mol,
    // RDKit❗❌:                                 const SmilesWriteParams &params,
    // RDKit❗❌:                                 const std::vector<int> &atomsToUse,
    // RDKit❗❌:                                 const std::vector<int> *bondsToUse,
    // RDKit❗❌:                                 const std::vector<std::string> *atomSymbols,
    // RDKit❗❌:                                 const std::vector<std::string> *bondSymbols) {
    // RDKit❗❌:   PRECONDITION(atomsToUse.size(), "no atoms provided");
    // RDKit❗❌:   PRECONDITION(
    // RDKit❗❌:       params.rootedAtAtom < 0 ||
    // RDKit❗❌:           static_cast<unsigned int>(params.rootedAtAtom) < mol.getNumAtoms(),
    // RDKit❗❌:       "rootedAtomAtom must be less than the number of atoms");
    // RDKit❗❌:   PRECONDITION(params.rootedAtAtom < 0 ||
    // RDKit❗❌:                    std::find(atomsToUse.begin(), atomsToUse.end(),
    // RDKit❗❌:                              params.rootedAtAtom) != atomsToUse.end(),
    // RDKit❗❌:                "rootedAtAtom not found in atomsToUse");
    // RDKit❗❌:   PRECONDITION(!atomSymbols || atomSymbols->size() >= mol.getNumAtoms(),
    // RDKit❗❌:                "bad atomSymbols vector");
    // RDKit❗❌:   PRECONDITION(!bondSymbols || bondSymbols->size() >= mol.getNumBonds(),
    // RDKit❗❌:                "bad bondSymbols vector");
    // RDKit❗❌:   if (!mol.getNumAtoms()) {
    // RDKit❗❌:     return "";
    // RDKit❗❌:   }
    // RDKit❗❌:   int rootedAtAtom = params.rootedAtAtom;
    // RDKit❗❌:
    // RDKit❗❌:   ROMol tmol(mol, true);
    // RDKit❗❌:   if (params.doIsomericSmiles) {
    // RDKit❗❌:     tmol.setProp(common_properties::_doIsoSmiles, 1);
    // RDKit❗❌:   }
    // RDKit❗❌:   std::string res;
    // RDKit❗❌:
    // RDKit❗❌:   boost::dynamic_bitset<> atomsInPlay(mol.getNumAtoms(), 0);
    // RDKit❗❌:   for (auto aidx : atomsToUse) {
    // RDKit❗❌:     atomsInPlay.set(aidx);
    // RDKit❗❌:   }
    // RDKit❗❌:   // figure out which bonds are actually in play:
    // RDKit❗❌:   boost::dynamic_bitset<> bondsInPlay(mol.getNumBonds(), 0);
    // RDKit❗❌:   if (bondsToUse) {
    // RDKit❗❌:     for (auto bidx : *bondsToUse) {
    // RDKit❗❌:       bondsInPlay.set(bidx);
    // RDKit❗❌:     }
    // RDKit❗❌:   } else {
    // RDKit❗❌:     PRECONDITION(
    // RDKit❗❌:         params.rootedAtAtom < 0 || MolOps::getMolFrags(mol).size() == 1,
    // RDKit❗❌:         "rootedAtAtom can only be used with molecules that have a single fragment");
    // RDKit❗❌:
    // RDKit❗❌:     for (auto aidx : atomsToUse) {
    // RDKit❗❌:       for (const auto &bndi : boost::make_iterator_range(
    // RDKit❗❌:                mol.getAtomBonds(mol.getAtomWithIdx(aidx)))) {
    // RDKit❗❌:         const Bond *bond = mol[bndi];
    // RDKit❗❌:         if (atomsInPlay[bond->getOtherAtomIdx(aidx)]) {
    // RDKit❗❌:           bondsInPlay.set(bond->getIdx());
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // copy over the rings that only involve atoms/bonds in this fragment:
    // RDKit❗❌:   if (mol.getRingInfo()->isInitialized()) {
    // RDKit❗❌:     tmol.getRingInfo()->reset();
    // RDKit❗❌:     tmol.getRingInfo()->initialize();
    // RDKit❗❌:     for (unsigned int ridx = 0; ridx < mol.getRingInfo()->numRings(); ++ridx) {
    // RDKit❗❌:       const INT_VECT &aring = mol.getRingInfo()->atomRings()[ridx];
    // RDKit❗❌:       bool keepIt = true;
    // RDKit❗❌:       for (auto aidx : aring) {
    // RDKit❗❌:         if (!atomsInPlay[aidx]) {
    // RDKit❗❌:           keepIt = false;
    // RDKit❗❌:           break;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       if (keepIt) {
    // RDKit❗❌:         const INT_VECT &bring = mol.getRingInfo()->bondRings()[ridx];
    // RDKit❗❌:         for (auto bidx : bring) {
    // RDKit❗❌:           if (!bondsInPlay[bidx]) {
    // RDKit❗❌:             keepIt = false;
    // RDKit❗❌:             break;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:         if (keepIt) {
    // RDKit❗❌:           tmol.getRingInfo()->addRing(aring, bring);
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   if (tmol.needsUpdatePropertyCache()) {
    // RDKit❗❌:     for (auto atom : tmol.atoms()) {
    // RDKit❗❌:       atom->updatePropertyCache(false);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   UINT_VECT ranks(tmol.getNumAtoms());
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<unsigned int> atomOrdering;
    // RDKit❗❌:   std::vector<unsigned int> bondOrdering;
    // RDKit❗❌:
    // RDKit❗❌:   // check stereochemistry:
    // RDKit❗❌:   if (params.doIsomericSmiles) {
    // RDKit❗❌:     if (!mol.hasProp(common_properties::_StereochemDone)) {
    // RDKit❗❌:       MolOps::assignStereochemistry(tmol, true);
    // RDKit❗❌:     } else {
    // RDKit❗❌:       tmol.setProp(common_properties::_StereochemDone, 1);
    // RDKit❗❌:       // we need the CIP codes:
    // RDKit❗❌:       for (auto aidx : atomsToUse) {
    // RDKit❗❌:         const Atom *oAt = mol.getAtomWithIdx(aidx);
    // RDKit❗❌:         std::string cipCode;
    // RDKit❗❌:         if (oAt->getPropIfPresent(common_properties::_CIPCode, cipCode)) {
    // RDKit❗❌:           tmol.getAtomWithIdx(aidx)->setProp(common_properties::_CIPCode,
    // RDKit❗❌:                                              cipCode);
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     // check for double bonds where the atoms defining stereo are not included
    // RDKit❗❌:     for (auto bnd : tmol.bonds()) {
    // RDKit❗❌:       if (bondsInPlay[bnd->getIdx()] && bnd->getBondType() == Bond::DOUBLE &&
    // RDKit❗❌:           bnd->getStereo() != Bond::BondStereo::STEREONONE) {
    // RDKit❗❌:         const auto &stereoAtoms = bnd->getStereoAtoms();
    // RDKit❗❌:         if (stereoAtoms.size() != 2) {
    // RDKit❗❌:           continue;
    // RDKit❗❌:         }
    // RDKit❗❌:         // check at both ends of the bond to see if the stereo atom is in play.
    // RDKit❗❌:         // If not and there's another neighbor atom there that *is* in play,
    // RDKit❗❌:         // then keep the stereochemistry and swap the stereo atoms. If not,
    // RDKit❗❌:         // remove the stereochemistry.
    // RDKit❗❌:         const std::vector<std::pair<int, const Atom *>>
    // RDKit❗❌:             stereoAtomsAndBondAtoms = {
    // RDKit❗❌:                 std::make_pair(stereoAtoms[0], bnd->getBeginAtom()),
    // RDKit❗❌:                 std::make_pair(stereoAtoms[1], bnd->getEndAtom())};
    // RDKit❗❌:         for (auto [stereoAtomIdx, bondAtom] : stereoAtomsAndBondAtoms) {
    // RDKit❗❌:           if (!atomsInPlay[stereoAtomIdx]) {
    // RDKit❗❌:             if (bondAtom->getDegree() > 2) {
    // RDKit❗❌:               bool updated = false;
    // RDKit❗❌:               for (auto nbrAt : tmol.atomNeighbors(bondAtom)) {
    // RDKit❗❌:                 if (nbrAt->getIdx() !=
    // RDKit❗❌:                         static_cast<unsigned int>(stereoAtomIdx) &&
    // RDKit❗❌:                     atomsInPlay[nbrAt->getIdx()]) {
    // RDKit❗❌:                   if (stereoAtomIdx == stereoAtoms[0]) {
    // RDKit❗❌:                     bnd->setStereoAtoms(nbrAt->getIdx(), stereoAtoms[1]);
    // RDKit❗❌:                   } else {
    // RDKit❗❌:                     bnd->setStereoAtoms(stereoAtoms[0], nbrAt->getIdx());
    // RDKit❗❌:                   }
    // RDKit❗❌:                   updated = true;
    // RDKit❗❌:                   if (bnd->getStereo() == Bond::BondStereo::STEREOZ ||
    // RDKit❗❌:                       bnd->getStereo() == Bond::BondStereo::STEREOCIS) {
    // RDKit❗❌:                     bnd->setStereo(Bond::BondStereo::STEREOTRANS);
    // RDKit❗❌:                   } else if (bnd->getStereo() == Bond::BondStereo::STEREOE ||
    // RDKit❗❌:                              bnd->getStereo() ==
    // RDKit❗❌:                                  Bond::BondStereo::STEREOTRANS) {
    // RDKit❗❌:                     bnd->setStereo(Bond::BondStereo::STEREOCIS);
    // RDKit❗❌:                   }
    // RDKit❗❌:                   break;
    // RDKit❗❌:                 }
    // RDKit❗❌:               }
    // RDKit❗❌:               if (!updated) {
    // RDKit❗❌:                 bnd->setStereo(Bond::BondStereo::STEREONONE);
    // RDKit❗❌:                 break;
    // RDKit❗❌:               }
    // RDKit❗❌:             } else {
    // RDKit❗❌:               bnd->setStereo(Bond::BondStereo::STEREONONE);
    // RDKit❗❌:               break;
    // RDKit❗❌:             }
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   if (params.canonical) {
    // RDKit❗❌:     bool breakTies = true;
    // RDKit❗❌:     bool includeChiralPresence = false;
    // RDKit❗❌:     bool includeRingStereo = true;
    // RDKit❗❌:     Canon::rankFragmentAtoms(
    // RDKit❗❌:         tmol, ranks, atomsInPlay, bondsInPlay, atomSymbols, bondSymbols,
    // RDKit❗❌:         breakTies, params.doIsomericSmiles, params.doIsomericSmiles,
    // RDKit❗❌:         !params.ignoreAtomMapNumbers, includeChiralPresence, includeRingStereo);
    // RDKit❗❌:     // std::cerr << "RANKS: ";
    // RDKit❗❌:     // std::copy(ranks.begin(), ranks.end(),
    // RDKit❗❌:     //           std::ostream_iterator<int>(std::cerr, " "));
    // RDKit❗❌:     // std::cerr << std::endl;
    // RDKit❗❌:     // MolOps::rankAtomsInFragment(tmol,ranks,atomsInPlay,bondsInPlay,atomSymbols,bondSymbols);
    // RDKit❗❌:   } else {
    // RDKit❗❌:     for (unsigned int i = 0; i < tmol.getNumAtoms(); ++i) {
    // RDKit❗❌:       ranks[i] = i;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: #ifdef VERBOSE_CANON
    // RDKit❗❌:   for (unsigned int tmpI = 0; tmpI < ranks.size(); tmpI++) {
    // RDKit❗❌:     std::cout << tmpI << " " << ranks[tmpI] << " "
    // RDKit❗❌:               << *(tmol.getAtomWithIdx(tmpI)) << std::endl;
    // RDKit❗❌:   }
    // RDKit❗❌: #endif
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<Canon::AtomColors> colors(tmol.getNumAtoms(), Canon::BLACK_NODE);
    // RDKit❗❌:   for (auto aidx : atomsToUse) {
    // RDKit❗❌:     colors[aidx] = Canon::WHITE_NODE;
    // RDKit❗❌:   }
    // RDKit❗❌:   std::vector<Canon::AtomColors>::iterator colorIt;
    // RDKit❗❌:   colorIt = colors.begin();
    // RDKit❗❌:   // loop to deal with the possibility that there might be disconnected
    // RDKit❗❌:   // fragments
    // RDKit❗❌:   while (colorIt != colors.end()) {
    // RDKit❗❌:     int nextAtomIdx = -1;
    // RDKit❗❌:
    // RDKit❗❌:     // find the next atom for a traverse
    // RDKit❗❌:     if (rootedAtAtom >= 0) {
    // RDKit❗❌:       nextAtomIdx = rootedAtAtom;
    // RDKit❗❌:       rootedAtAtom = -1;
    // RDKit❗❌:     } else {
    // RDKit❗❌:       unsigned int nextRank = rdcast<unsigned int>(tmol.getNumAtoms()) + 1;
    // RDKit❗❌:       for (auto i : atomsToUse) {
    // RDKit❗❌:         if (colors[i] == Canon::WHITE_NODE && ranks[i] < nextRank) {
    // RDKit❗❌:           nextRank = ranks[i];
    // RDKit❗❌:           nextAtomIdx = i;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     CHECK_INVARIANT(nextAtomIdx >= 0, "no start atom found");
    // RDKit❗❌:     auto subSmi = SmilesWrite::FragmentSmilesConstruct(
    // RDKit❗❌:         tmol, nextAtomIdx, colors, ranks, params, atomOrdering, bondOrdering,
    // RDKit❗❌:         &atomsInPlay, &bondsInPlay, atomSymbols, bondSymbols);
    // RDKit❗❌:
    // RDKit❗❌:     res += subSmi;
    // RDKit❗❌:     colorIt = std::find(colors.begin(), colors.end(), Canon::WHITE_NODE);
    // RDKit❗❌:     if (colorIt != colors.end()) {
    // RDKit❗❌:       res += ".";
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   mol.setProp(common_properties::_smilesAtomOutputOrder, atomOrdering, true);
    // RDKit❗❌:   mol.setProp(common_properties::_smilesBondOutputOrder, bondOrdering, true);
    // RDKit❗❌:
    // RDKit❗❌:   return res;
    // RDKit❗❌: }  // end of MolFragmentToSmiles()
    // END COMPLETE RDKit .6 CHEM27 SmilesParse::MolFragmentToSmiles
    // This existing detached entry owns the full source call. Its preparation
    // owner performs the new controller prepass before this single whole-mask
    // rank call, outside the component loop; no component extraction is added.
    // Existing random/typed-cache behavior remains qualified and unchanged.
    // Behavior review: preflight, full-index masks, S59 preparation and S60
    // ranking stay on one cloned topology; each component uses the existing
    // source-shaped masked traversal and emits original AtomId/BondId rows.
    // Complexity review: selection/preparation/ranking are O(V+E) plus the
    // selected rank owner; per-component DFS uses the same full-size state
    // vectors and adjacency scans as the existing writer. No subgraph remap is
    // built, at the cost of retaining the original topology for the call.
    record
        .topology
        .validate()
        .map_err(|error| FragmentWriteInputError::InvalidTopology(error.to_string()))?;
    validate_fragment_write_inputs(
        &record.topology,
        params,
        atoms_to_use,
        atom_symbols,
        bond_symbols,
    )?;
    let masks =
        build_fragment_selection_masks(&record.topology, params, atoms_to_use, bonds_to_use)?;
    let prepared = prepare_fragment_stereo(
        record,
        params,
        atoms_to_use,
        &masks,
        source_rings,
        existing_valence,
    )?;
    let original_ranks =
        rank_prepared_fragment(&prepared, &masks, atom_symbols, bond_symbols, params)?;
    let PreparedFragmentStereo {
        topology: mut topology,
        retained_rings,
        ranking_rings,
        valence,
        ..
    } = prepared;
    let prepared_rings = ranking_rings.as_ref().or(retained_rings.as_ref());
    let ranks = original_ranks
        .into_iter()
        .map(|rank| rank as i64)
        .collect::<Vec<_>>();

    let mut colors = vec![AtomColor::Black; topology.atoms.len()];
    for (atom_index, selected) in masks.atoms_in_play.iter().copied().enumerate() {
        if selected {
            colors[atom_index] = AtomColor::White;
        }
    }
    let ring_bonds = find_ring_bonds(&topology);
    let mut rooted_at_atom = params.rooted_at_atom.map(AtomId::index);
    let mut text = PropertyText::new();
    let mut atom_order = Vec::new();
    let mut bond_order = Vec::new();

    while masks
        .atoms_in_play
        .iter()
        .enumerate()
        .any(|(atom_index, selected)| *selected && colors[atom_index] == AtomColor::White)
    {
        let start = if let Some(root) = rooted_at_atom.take() {
            root
        } else {
            atoms_to_use
                .iter()
                .map(|atom| atom.index())
                .filter(|atom| colors[*atom] == AtomColor::White)
                .min_by_key(|atom| ranks[*atom])
                .ok_or_else(|| {
                    FragmentWriteInputError::Writer(SmilesParseError::Model(
                        "fragment traversal has white atoms but no selected start atom".into(),
                    ))
                })?
        };

        // BEGIN RDKIT CPP FUNCTION SmilesWrite::FragmentSmilesConstruct selected kekulization
        // RDKit❗✔️:   if (params.doKekule) {
        // RDKit❗✔️:     if (atomsInPlay && bondsInPlay) {
        // RDKit❗✔️:       MolOps::details::KekulizeFragment(static_cast<RWMol &>(mol),
        // RDKit❗✔️:                                         *atomsInPlay, *bondsInPlay);
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // END RDKIT CPP FUNCTION SmilesWrite::FragmentSmilesConstruct selected kekulization
        if params.kekule {
            topology = cosmolkit_core::kekulize_selected_fragment(
                &topology,
                &masks.atoms_in_play,
                &masks.bonds_in_play,
                &KekulizeParams::default(),
            )
            .map_err(SmilesParseError::WriterKekulize)?
            .topology;
        }
        // BEGIN RDKIT CPP FUNCTION SmilesWrite::FragmentSmilesConstruct enhanced-stereo dispatch
        // RDKit❗✔️:   if (params.canonical && params.doIsomericSmiles) {
        // RDKit❗✔️:     Canon::canonicalizeEnhancedStereo(mol, &ranks);
        // RDKit❗✔️:   }
        // END RDKIT CPP FUNCTION SmilesWrite::FragmentSmilesConstruct enhanced-stereo dispatch
        // Behavior review: this source stage runs on the full prepared topology
        // after optional selected-mask kekulization and before fragment traversal.
        // Its sparse reference map carries the exact first rank-sorted atom used
        // by canonicalizeFragment's later `_stereoGroup` branch.
        // Complexity review: the map stores one reference per nonempty atom
        // group; normalization copies/sorts only group members and does not
        // allocate a topology-sized marker array.
        let stereo_group_references = if params.canonical && params.isomeric_smiles {
            canonicalize_enhanced_stereo(&mut topology, &ranks)?
        } else {
            BTreeMap::new()
        };

        let mut cycle_colors = colors.clone();
        let mut ring_closures = vec![Vec::new(); topology.atoms.len()];
        dfs_find_cycles(
            &topology,
            start,
            None,
            &mut cycle_colors,
            &ranks,
            &ring_bonds,
            Some(&masks.bonds_in_play),
            bond_symbols,
            &mut ring_closures,
            None,
        )?;
        let mut stack = Vec::with_capacity(topology.atoms.len() + topology.bonds.len());
        let mut ring_ids = vec![None; topology.bonds.len()];
        let mut traversal_ring_closure_bonds = vec![false; topology.bonds.len()];
        let mut available_ring_ids = vec![true; MAX_CYCLES];
        let mut atom_traversal_bond_order = vec![Vec::new(); topology.atoms.len()];
        let source_ranks = ranks.iter().map(|rank| *rank as u32).collect::<Vec<_>>();
        dfs_build_stack(
            &mut topology,
            start,
            None,
            &mut colors,
            &source_ranks,
            &ring_bonds,
            &ring_closures,
            &mut ring_ids,
            &mut available_ring_ids,
            &mut stack,
            &mut atom_traversal_bond_order,
            &mut traversal_ring_closure_bonds,
            Some(&masks.bonds_in_play),
            bond_symbols.map(DfsBondSymbols::Utf8),
            None,
        )?;
        let computed_valence = if params.kekule {
            Some(
                cosmolkit_core::assign_valence_with_options_for_topology(
                    &topology,
                    ValenceModel::RdkitLike,
                    false,
                )
                .map_err(|error| SmilesParseError::WriterValence(error.to_string()))?,
            )
        } else {
            None
        };
        let emitted_valence = computed_valence
            .as_ref()
            .unwrap_or_else(|| valence.as_ref());
        let chiral_adjustments = compute_chiral_adjustments(
            &mut topology,
            emitted_valence,
            prepared_rings,
            params.isomeric_smiles,
            start,
            &ring_closures,
            &atom_traversal_bond_order,
            &stack,
            &stereo_group_references,
            Some(&masks.atoms_in_play),
            Some(&masks.bonds_in_play),
            &traversal_ring_closure_bonds,
        )?;
        text.extend_bytes(
            (&write_mol_stack(
                &topology,
                emitted_valence,
                &stack,
                &chiral_adjustments,
                params,
                atom_symbols,
                bond_symbols,
            )?)
                .as_ref(),
        );
        atom_order.extend(stack.iter().filter_map(|element| match element {
            MolStackElem::Atom(atom) => Some(AtomId::new(*atom)),
            _ => None,
        }));
        bond_order.extend(stack.iter().filter_map(|element| match element {
            MolStackElem::Bond { bond, .. } => Some(*bond),
            _ => None,
        }));
        if masks
            .atoms_in_play
            .iter()
            .enumerate()
            .any(|(atom_index, selected)| *selected && colors[atom_index] == AtomColor::White)
        {
            text.push_byte(b'.');
        }
    }

    Ok(SmilesWriteOutput {
        text,
        atom_order,
        bond_order,
    })
}

/// Writes detached fragment CXSMILES using the fragment traversal's emitted rows.
pub fn write_fragment_cx_smiles<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
    params: &CxSmilesWriteParams,
    atoms_to_use: &[AtomId],
    bonds_to_use: Option<&[BondId]>,
    atom_symbols: Option<&[String]>,
    bond_symbols: Option<&[String]>,
    source_rings: Option<&cosmolkit_core::RingInfo>,
    existing_valence: Option<&ValenceAssignment>,
) -> Result<PropertyText, FragmentWriteInputError> {
    let record = record.into();
    // BEGIN RDKIT CPP FUNCTION SmilesWrite.cpp::MolFragmentToCXSmiles
    // RDKit❗✔️:   auto res = MolFragmentToSmiles(mol, params, atomsToUse,
    // RDKit❗✔️:                                  bondsToUse, atomSymbols,
    // RDKit❗✔️:                                  bondSymbols);
    // RDKit❗✔️:   auto cxext = SmilesWrite::getCXExtensions(mol);
    // RDKit❗✔️:   if (!cxext.empty()) {
    // RDKit❗✔️:     res += " " + cxext;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // END RDKIT CPP FUNCTION SmilesWrite.cpp::MolFragmentToCXSmiles
    // Behavior review: the fragment traversal supplies full-topology output
    // rows; CX fields are appended from the original immutable record, matching
    // the pinned wrapper's base-write-then-extension call order. The CX writer
    // applies the requested field mask, and its SGroup reverse maps preserve
    // the source's zero-initialized rows for members outside the fragment.
    // Complexity review: the existing fragment writer performs the source
    // topology clone/traversal once and returns its maps; CX emission reuses
    // those maps without a second topology copy or map reconstruction.
    let output = write_fragment_smiles_output(
        record,
        &params.smiles,
        atoms_to_use,
        bonds_to_use,
        atom_symbols,
        bond_symbols,
        source_rings,
        existing_valence,
    )?;
    let extension = write_cx_extensions(
        record,
        params.fields,
        &output.atom_order,
        &output.bond_order,
        None,
        params.coordinate_selection,
    )?;
    if extension.is_empty() {
        Ok(output.text)
    } else {
        {
            let mut text = output.text;
            text.push_byte(b' ');
            text.extend_bytes(extension.as_bytes());
            Ok(text)
        }
    }
}

fn write_smiles_output<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
    params: &SmilesWriteParams,
    doing_cx_smiles: bool,
) -> Result<SmilesWriteOutput, SmilesParseError> {
    let record = record.into();
    write_smiles_output_with_random_stream(record, params, doing_cx_smiles, None)
}

fn write_smiles_output_with_random_stream<'record>(
    record: impl Into<crate::SmilesRecordView<'record>>,
    params: &SmilesWriteParams,
    doing_cx_smiles: bool,
    mut random_stream: Option<&mut cosmolkit_core::RdkitRandomGenerator<'_>>,
) -> Result<SmilesWriteOutput, SmilesParseError> {
    let record = record.into();
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles root validation
    // RDKit✔️❌:   if (!mol.getNumAtoms()) {
    // RDKit✔️❌:     return "";
    // RDKit✔️❌:   }
    // RDKit✔️❌:   PRECONDITION(
    // RDKit✔️❌:       params.rootedAtAtom < 0 ||
    // RDKit✔️❌:           static_cast<unsigned int>(params.rootedAtAtom) < mol.getNumAtoms(),
    // RDKit✔️❌:       "rootedAtAtom must be less than the number of atoms");
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles root validation
    // Detached records need structural validation before indexed access; this
    // adds an O(V + E) pass before the source's constant-time root precondition.
    record
        .topology
        .validate()
        .map_err(|error| SmilesParseError::Model(error.to_string()))?;
    if record.topology.atoms.is_empty() {
        return Ok(SmilesWriteOutput {
            text: PropertyText::new(),
            atom_order: Vec::new(),
            bond_order: Vec::new(),
        });
    }
    if let Some(root) = params.rooted_at_atom {
        if root.index() >= record.topology.atoms.len() {
            return Err(SmilesParseError::WriterRootAtomOutOfRange {
                atom_index: root.index(),
                atom_count: record.topology.atoms.len(),
            });
        }
    }
    reject_unmodeled_stereochemical_writing(&record.topology, doing_cx_smiles)?;

    // RDKit's public writer performs stereo preparation and direction
    // canonicalization on a private molecule copy. Keep the detached record
    // immutable and make those same temporary changes on this working block.
    let mut topology = record.topology.clone();
    let original_atom_maps = if params.ignore_atom_map_numbers {
        let original = params.canonical.then(|| {
            topology
                .atoms
                .iter()
                .map(|atom| atom.atom_map().filter(|map| *map != 0))
                .collect::<Vec<_>>()
        });
        // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles map setup
        // RDKit✔️✔️:       if (params.ignoreAtomMapNumbers) {
        // RDKit✔️✔️:         atomMapNums[atom->getIdx()] = atom->getAtomMapNum();
        // RDKit✔️✔️:         atom->setAtomMapNum(0);
        // RDKit✔️✔️:       }
        // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles map setup
        for atom in &mut topology.atoms {
            atom.set_atom_map(None);
        }
        original
    } else {
        None
    };
    let source_valence = if doing_cx_smiles || !params.include_dative_bonds {
        Some(
            cosmolkit_core::assign_valence_with_options_for_topology(
                &topology,
                ValenceModel::RdkitLike,
                false,
            )
            .map_err(|error| SmilesParseError::WriterValence(error.to_string()))?,
        )
    } else {
        None
    };
    let components = connected_components(&topology);
    let component_uses_subset_stereo_copy = components
        .iter()
        .map(|component| {
            writer_component_uses_subset_stereo_copy(&topology, component, components.len())
        })
        .collect::<Vec<_>>();
    let stereochem_done_marker_is_computed = record
        .properties
        .prop("_StereochemDone")
        .map(|_| record.properties.is_prop_computed("_StereochemDone"))
        .transpose()?;
    prepare_writer_stereochemistry(
        &mut topology,
        stereochem_done_marker_is_computed,
        record.rings,
        params,
        &components,
    )?;
    let mut dative_donors = Vec::new();
    let mut hydrogen_bond_atoms = Vec::new();
    // RDKit✔️✔️: if (!doingCXSmiles || !includeStereoGroups) {
    // RDKit✔️✔️:   std::vector<StereoGroup> noStereoGroups;
    // RDKit✔️✔️:   tmol->setStereoGroups(noStereoGroups);
    // RDKit✔️✔️: }
    if doing_cx_smiles {
        // RDKit✔️✔️: if (doingCXSmiles || !params.includeDativeBonds) {
        // RDKit✔️✔️:   for (auto bond : tmol->bonds()) {
        // RDKit✔️✔️:     if (bond->getBondType() == Bond::DATIVE) {
        // RDKit✔️✔️:       bond->setBondType(Bond::SINGLE);
        // RDKit✔️✔️:       bond->getBeginAtom()->calcExplicitValence(false);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // RDKit✔️✔️: if (doingCXSmiles) {
        // RDKit✔️✔️:   for (auto bond : tmol->bonds()) {
        // RDKit✔️✔️:     if (bond->getBondType() == Bond::HYDROGEN) {
        // RDKit✔️✔️:       bond->setBondType(Bond::SINGLE);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        for bond in &mut topology.bonds {
            if bond.order() == BondOrder::Dative {
                dative_donors.push(bond.begin());
                bond.set_order(BondOrder::Single);
            } else if bond.order() == BondOrder::Hydrogen {
                hydrogen_bond_atoms.push(bond.begin());
                hydrogen_bond_atoms.push(bond.end());
                bond.set_order(BondOrder::Single);
            }
            if !matches!(
                bond.direction(),
                BondDirection::EndDownRight | BondDirection::EndUpRight
            ) {
                bond.set_direction(BondDirection::None);
            }
            if bond.stereo() == BondStereo::Any {
                bond.set_stereo(BondStereo::None)
                    .map_err(SmilesParseError::WriterStereoBond)?;
            }
        }
        topology.adjacency =
            cosmolkit_model::AdjacencyList::from_topology(topology.atoms.len(), &topology.bonds);
    }
    // BEGIN RDKIT CPP FUNCTION detail::MolToSmiles ordinary bond cleanup
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
    // END RDKIT CPP FUNCTION detail::MolToSmiles ordinary bond cleanup
    // Undefined stereo changes canonical bond invariants. Clear it on the
    // writer copy before ranking, as in the source's one O(E) cleanup pass.
    if !doing_cx_smiles {
        for bond in &mut topology.bonds {
            if matches!(
                bond.direction(),
                BondDirection::Unknown | BondDirection::EitherDouble
            ) {
                bond.set_direction(BondDirection::None);
            }
            if bond.stereo() == BondStereo::Any {
                bond.set_stereo(BondStereo::None)
                    .map_err(SmilesParseError::WriterStereoBond)?;
            }
        }
    }
    if !doing_cx_smiles && !params.include_dative_bonds {
        // BEGIN RDKIT CPP FUNCTION detail::MolToSmiles includeDativeBonds conversion
        // RDKit❗✔️:     if (doingCXSmiles || !params.includeDativeBonds) {
        // RDKit❗✔️:       for (auto bond : tmol->bonds()) {
        // RDKit❗✔️:         if (bond->getBondType() == Bond::DATIVE) {
        // RDKit❗✔️:           bond->setBondType(Bond::SINGLE);
        // RDKit❗✔️:           bond->getBeginAtom()->calcExplicitValence(false);
        // RDKit❗✔️:         }
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // END RDKIT CPP FUNCTION detail::MolToSmiles includeDativeBonds conversion
        for bond in &mut topology.bonds {
            if bond.order() == BondOrder::Dative {
                dative_donors.push(bond.begin());
                bond.set_order(BondOrder::Single);
            }
        }
    }
    if !doing_cx_smiles {
        topology.stereo_groups.clear();
    }
    // RDKit✔️✔️:     for (auto atom : tmol->atoms()) {
    // RDKit✔️✔️:       atom->updatePropertyCache(false);
    // RDKit✔️✔️:     }
    let mut valence = cosmolkit_core::assign_valence_with_options_for_topology(
        &topology,
        ValenceModel::RdkitLike,
        false,
    )
    .map_err(|error| SmilesParseError::WriterValence(error.to_string()))?;
    if let Some(source_valence) = &source_valence {
        // RDKit✔️✔️: bond->setBondType(Bond::SINGLE);
        // RDKit✔️✔️: // update the explicit valence of the begin atom since the implicit
        // RDKit✔️✔️: // valence will no longer be properly perceived
        // RDKit✔️✔️: bond->getBeginAtom()->calcExplicitValence(false);
        // The source deliberately retains the donor's already perceived
        // implicit valence while recalculating only its explicit valence.
        for atom in &dative_donors {
            valence.implicit_hydrogens[atom.index()] =
                source_valence.implicit_hydrogens[atom.index()];
        }
        // RDKit converts hydrogen bonds to single bonds without recalculating
        // either endpoint's cached valence.
        for atom in &hydrogen_bond_atoms {
            valence.explicit_valence[atom.index()] = source_valence.explicit_valence[atom.index()];
            valence.implicit_hydrogens[atom.index()] =
                source_valence.implicit_hydrogens[atom.index()];
        }
    }
    let rings =
        fast_find_rings_from_parts(topology.atoms.len(), &topology.bonds, &topology.adjacency)
            .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
    // RDKit✔️✔️:       if (!tmol->hasProp(common_properties::_StereochemDone)) {
    // RDKit✔️✔️:         MolOps::assignStereochemistry(*tmol, params.cleanStereo);
    // RDKit✔️✔️:       }
    // The marker-aware preparation above is the source's sole stereo assignment.
    // A second direction-to-stereo pass would reinterpret residual slashes on
    // newly transformed bonds even when the source marker forbids reassignment.
    let mut ranking_topology = topology.clone();
    if doing_cx_smiles {
        for atom in dative_donors
            .iter()
            .chain(hydrogen_bond_atoms.iter())
            .copied()
        {
            let total_hydrogens =
                usize::from(ranking_topology.atoms[atom.index()].explicit_hydrogens())
                    + usize::try_from(valence.implicit_hydrogens[atom.index()].max(0))
                        .unwrap_or(usize::MAX);
            ranking_topology.atoms[atom.index()]
                .set_explicit_hydrogens(u8::try_from(total_hydrogens).unwrap_or(u8::MAX));
            ranking_topology.atoms[atom.index()].set_no_implicit(true);
        }
    }
    let mut component_ranks = Vec::with_capacity(components.len());
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles canonical rank policy
    // RDKit✔️✔️:       const bool breakTies = true;
    // RDKit✔️✔️:       const bool includeChiralPresence = false;
    // RDKit✔️✔️:       const bool includeIsotopes = params.doIsomericSmiles;
    // RDKit✔️✔️:       const bool includeChirality = params.doIsomericSmiles;
    // RDKit✔️✔️:       const bool includeStereoGroups = params.doIsomericSmiles;
    // RDKit✔️✔️:       const bool useNonStereoRanks = false;
    // RDKit✔️✔️:       const bool includeAtomMaps = true;
    // RDKit✔️✔️:       Canon::rankMolAtoms(*tmol, ranks, breakTies, includeChirality,
    // RDKit✔️✔️:                           includeIsotopes, includeAtomMaps,
    // RDKit✔️✔️:                           includeChiralPresence, includeStereoGroups,
    // RDKit✔️✔️:                           useNonStereoRanks);
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles canonical rank policy
    // Complexity review: policy propagation adds constant-time field assignments per
    // component. Fragment construction, rank refinement, allocations, loop nesting,
    // and component mapping are unchanged; no scan, clone, or lookup was added.
    let mut kekulize_fragments = if params.kekule {
        // BEGIN RDKIT CPP FUNCTION MolOps::getTheFrags fragment copies before ranking
        // RDKit❗❌:   if (nFrags == 1) {
        // RDKit❗❌:     res.emplace_back(new RWMol(mol));
        // RDKit❗❌:   } else {
        // RDKit❗❌:     if (comp.size() == 1 ||
        // RDKit❗❌:         (nFrags > 3 && !fragmentHasChallengingFeatures(comp, atomsInFrag))) {
        // RDKit❗❌:       auto submol = copyMolSubset(mol, atoms, info, opts);
        // RDKit❗❌:       res.push_back(std::move(submol));
        // RDKit❗❌:     } else {
        // RDKit❗❌:       res.emplace_back(new RWMol(mol));
        // RDKit❗❌:       auto &frag = res.back();
        // RDKit❗❌:       frag->beginBatchEdit();
        // RDKit❗❌:       for (unsigned int idx = 0; idx < mol.getNumAtoms(); ++idx) {
        // RDKit❗❌:         if (!atomsInFrag[idx]) { frag->removeAtom(idx); }
        // RDKit❗❌:       }
        // RDKit❗❌:       frag->commitBatchEdit();
        // RDKit❗❌:     }
        // RDKit❗❌:   }
        // END RDKIT CPP FUNCTION MolOps::getTheFrags fragment copies before ranking
        // Local adapters keep source maps while allocating a membership mask and
        // copying each component topology before any source fragment is ranked.
        components
            .iter()
            .map(|component| {
                if components.len() == 1 {
                    Ok(WriterStereoFragment {
                        topology: topology.clone(),
                        source_atoms: (0..topology.atoms.len()).map(AtomId::new).collect(),
                        source_bonds: (0..topology.bonds.len()).map(BondId::new).collect(),
                    })
                } else {
                    let mut inside = vec![false; topology.atoms.len()];
                    for source_atom in component {
                        inside[*source_atom] = true;
                    }
                    let clears_computed_props = component.len() == 1
                        || (components.len() > 3
                            && !fragment_has_challenging_features(&topology, component, &inside));
                    if clears_computed_props {
                        extract_writer_subset_fragment(&topology, component, &inside)
                    } else {
                        extract_writer_preserving_fragment(&topology, &inside)
                    }
                }
            })
            .collect::<Result<Vec<_>, _>>()?
    } else {
        Vec::new()
    };
    for (component_index, component) in components.iter().enumerate() {
        component_ranks.push(if params.canonical {
            let mut policy = canonical_rank::CanonicalRankPolicy::default();
            policy.include_chirality = params.isomeric_smiles;
            policy.include_isotopes = params.isomeric_smiles;
            policy.include_stereo_groups = params.isomeric_smiles;
            canonical_rank::rank_component_atoms_with_policy(&ranking_topology, component, policy)
                .map_err(|error| match error {
                    cosmolkit_core::CanonicalRankError::StereoGroup(cause) => {
                        SmilesParseError::StereoGroup(cause)
                    }
                    other => SmilesParseError::CanonicalRank(other.to_string()),
                })?
                .into_iter()
                .map(|rank| rank as i64)
                .collect::<Vec<_>>()
        } else {
            (0..topology.atoms.len())
                .map(|index| index as i64)
                .collect::<Vec<_>>()
        });
        // RDKit restores ignored map values after this fragment's rank and
        // before constructing that fragment's SMILES. Keep the source rows
        // mapped to their saved values while later component ranks continue to
        // use the unchanged, map-cleared `ranking_topology`.
        if let Some(original_atom_maps) = original_atom_maps.as_ref() {
            // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles map restoration
            // RDKit✔️✔️:       if (params.ignoreAtomMapNumbers) {
            // RDKit✔️✔️:         for (auto atom : tmol->atoms()) {
            // RDKit✔️✔️:           atom->setAtomMapNum(atomMapNums[atom->getIdx()]);
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles map restoration
            // Source map number zero clears the map property on restoration.
            for source_atom in component {
                topology.atoms[*source_atom].set_atom_map(original_atom_maps[*source_atom]);
            }
        }
        if params.kekule {
            // BEGIN RDKIT CPP FUNCTION FragmentSmilesConstruct Kekulize branch
            // RDKit❗❌:   if (params.doKekule) {
            // RDKit❗❌:     if (atomsInPlay && bondsInPlay) {
            // RDKit❗❌:       MolOps::details::KekulizeFragment(static_cast<RWMol &>(mol), *atomsInPlay, *bondsInPlay);
            // RDKit❗❌:     } else {
            // RDKit❗❌:       MolOps::Kekulize(static_cast<RWMol &>(mol));
            // RDKit❗❌:     }
            // RDKit❗❌:   }
            // END RDKIT CPP FUNCTION FragmentSmilesConstruct Kekulize branch
            // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/Kekulize.cpp :: KekulizeFragment atom error
            // RDKit✔️✔️:           throw AtomKekulizeException(msg, atom->getIdx());
            // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/Kekulize.cpp :: KekulizeFragment atom error
            // FragmentSmilesConstruct receives an extracted component from
            // getTheFrags, so its full Kekulize call reports component-local
            // atom rows. Reuse the writer's source-mapped extraction branches
            // and project core Kekulize's atom/bond fields and complete
            // valence cache back to the private full-topology copy.
            let fragment = &mut kekulize_fragments[component_index];
            if let Some(original_atom_maps) = original_atom_maps.as_ref() {
                for (local_index, source_atom) in fragment.source_atoms.iter().copied().enumerate()
                {
                    fragment.topology.atoms[local_index]
                        .set_atom_map(original_atom_maps[source_atom.index()]);
                }
            }
            let kekulized =
                kekulize_writer_fragment(&fragment.topology, &valence, &fragment.source_atoms)?;
            project_kekulized_writer_valence(
                &mut valence,
                kekulized.final_valence.as_ref(),
                &fragment.source_atoms,
            )?;
            for (local_index, source_atom) in fragment.source_atoms.iter().copied().enumerate() {
                let kekulized_atom = &kekulized.topology.atoms[local_index];
                let target = &mut topology.atoms[source_atom.index()];
                target.set_aromatic(kekulized_atom.is_aromatic());
                target.set_explicit_hydrogens(kekulized_atom.explicit_hydrogens());
                target.set_no_implicit(kekulized_atom.no_implicit());
            }
            for (local_index, source_bond) in fragment.source_bonds.iter().copied().enumerate() {
                let kekulized_bond = &kekulized.topology.bonds[local_index];
                let target = &mut topology.bonds[source_bond.index()];
                target.set_order(kekulized_bond.order());
                target.set_aromatic(kekulized_bond.is_aromatic());
            }
        }
    }
    if params.canonical && !params.isomeric_smiles {
        prepare_canonical_nonisomeric_stereo_fallback(
            &mut topology,
            stereochem_done_marker_is_computed,
            record.rings,
            &components,
        )?;
    }
    let ring_bonds = writer_ring_bonds(&topology, record.rings, components.len());
    let mut colors = vec![AtomColor::White; topology.atoms.len()];
    let mut fragments = Vec::new();

    // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles canonical and non-canonical fragments
    // RDKit✔️✔️:     // Sort the vfragsmi, but also sort the atom and bond order vectors into
    // RDKit✔️✔️:     // the same order
    // RDKit✔️✔️:     typedef std::tuple<std::string, std::vector<unsigned int>,
    // RDKit✔️✔️:                        std::vector<unsigned int>> tplType;
    // RDKit✔️✔️:     std::vector<tplType> tmp(vfragsmi.size());
    // RDKit✔️✔️:     for (unsigned int ti = 0; ti < vfragsmi.size(); ++ti) {
    // RDKit✔️✔️:       tmp[ti] = std::make_tuple(vfragsmi[ti], allAtomOrdering[ti],
    // RDKit✔️✔️:                                 allBondOrdering[ti]);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     std::sort(tmp.begin(), tmp.end());
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
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles non-canonical fragments
    for (component, ranks) in components.into_iter().zip(component_ranks) {
        let component_stereo_groups = writer_component_stereo_groups(
            &topology,
            &component,
            component_uses_subset_stereo_copy[fragments.len()],
        )?;
        let original_stereo_groups =
            std::mem::replace(&mut topology.stereo_groups, component_stereo_groups);
        // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles rooted component
        // RDKit✔️🔝:     rootedAtAtom = -1;
        // RDKit✔️🔝:     if (params.rootedAtAtom >= 0 && atsPresent[params.rootedAtAtom]) {
        // RDKit✔️🔝:       rootedAtAtom = params.rootedAtAtom - atsPresent.find_first();
        // RDKit✔️🔝:     }
        // RDKit✔️🔝:     if (rootedAtAtom >= 0) {
        // RDKit✔️🔝:       nextAtomIdx = rootedAtAtom;
        // RDKit✔️🔝:       rootedAtAtom = -1;
        // RDKit✔️🔝:     } else {
        // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles rooted component
        // Sorted component rows make root membership/remapping O(V) across the
        // writer; the source rebuilds an O(V)-bit set once per component.
        let explicit_root = params
            .rooted_at_atom
            .filter(|root| component.contains(&root.index()));
        // BEGIN RDKIT CPP FUNCTION SmilesWrite.cpp::detail::MolToSmiles random root
        // RDKit❗❌:     rootedAtAtom = fragsRootedAtAtom[fragIdx];
        // RDKit❗❌:     if (params.doRandom && rootedAtAtom == -1) {
        // RDKit❗❌:       rootedAtAtom = getRandomGenerator()() % tmol->getNumAtoms();
        // RDKit❗❌:     }
        // END RDKIT CPP FUNCTION SmilesWrite.cpp::detail::MolToSmiles random root
        // Behavior: the presence of the borrowed shared stream represents the
        // private random mode. A component-specific explicit root consumes no
        // draw; otherwise one source draw is reduced modulo that component's
        // atom count and mapped back to its original atom ID. The public random
        // entry remains deferred until cycle and stack draws are wired.
        // Complexity: one constant-time draw/modulo and indexed component
        // lookup per eligible component, with no extra collection; the shared
        // stream adds one known outer mutex acquisition versus RDKit's global.
        let start = if let Some(root) = explicit_root {
            // RDKit subtracts the fragment's first source atom index, then
            // uses that value in the compact fragment's atom-index domain.
            let rooted_fragment_index = root.index() - component[0];
            *component.get(rooted_fragment_index).ok_or(
                SmilesParseError::WriterRootAtomOutOfRange {
                    atom_index: rooted_fragment_index,
                    atom_count: component.len(),
                },
            )?
        } else if let Some(random_stream) = random_stream.as_mut() {
            let local_root = (random_stream.next_u32() as usize) % component.len();
            component[local_root]
        } else if params.canonical {
            component
                .iter()
                .copied()
                .min_by_key(|atom| ranks[*atom])
                .expect("connected component is nonempty")
        } else {
            component[0]
        };
        // FragmentSmilesConstruct canonicalizes enhanced stereo after its
        // optional fragment kekulization and before canonicalizeFragment.
        // The whole writer reaches the same shared source stage here; ordinary
        // SMILES has already removed stereo groups, so this is a no-op there.
        let stereo_group_references = if params.canonical && params.isomeric_smiles {
            canonicalize_enhanced_stereo(&mut topology, &ranks)?
        } else {
            BTreeMap::new()
        };
        let mut cycle_colors = colors.clone();
        let mut ring_closures = vec![Vec::new(); topology.atoms.len()];
        dfs_find_cycles(
            &topology,
            start,
            None,
            &mut cycle_colors,
            &ranks,
            &ring_bonds,
            None,
            None,
            &mut ring_closures,
            random_stream.as_deref_mut(),
        )?;
        let mut stack = Vec::with_capacity(topology.atoms.len() + topology.bonds.len());
        let mut ring_ids = vec![None; topology.bonds.len()];
        let mut traversal_ring_closure_bonds = vec![false; topology.bonds.len()];
        let mut available_ring_ids = vec![true; MAX_CYCLES];
        let mut atom_traversal_bond_order = vec![Vec::new(); topology.atoms.len()];
        let source_ranks = ranks.iter().map(|rank| *rank as u32).collect::<Vec<_>>();
        dfs_build_stack(
            &mut topology,
            start,
            None,
            &mut colors,
            &source_ranks,
            &ring_bonds,
            &ring_closures,
            &mut ring_ids,
            &mut available_ring_ids,
            &mut stack,
            &mut atom_traversal_bond_order,
            &mut traversal_ring_closure_bonds,
            None,
            None,
            random_stream.as_deref_mut(),
        )?;
        let chiral_adjustments = compute_chiral_adjustments(
            &mut topology,
            &valence,
            Some(&rings),
            params.isomeric_smiles,
            start,
            &ring_closures,
            &atom_traversal_bond_order,
            &stack,
            &stereo_group_references,
            None,
            None,
            &traversal_ring_closure_bonds,
        )?;
        let text = write_mol_stack(
            &topology,
            &valence,
            &stack,
            &chiral_adjustments,
            params,
            None,
            None,
        )?;
        let atom_order = stack
            .iter()
            .filter_map(|element| match element {
                MolStackElem::Atom(atom) => Some(AtomId::new(*atom)),
                _ => None,
            })
            .collect();
        let bond_order = stack
            .iter()
            .filter_map(|element| match element {
                MolStackElem::Bond { bond, .. } => Some(*bond),
                _ => None,
            })
            .collect();
        fragments.push(FragmentWriteOutput {
            text,
            atom_order,
            bond_order,
        });
        topology.stereo_groups = original_stereo_groups;
    }
    if params.canonical {
        fragments.sort_by(|left, right| {
            left.text
                .cmp(&right.text)
                .then_with(|| left.atom_order.cmp(&right.atom_order))
                .then_with(|| left.bond_order.cmp(&right.bond_order))
        });
    }
    let mut text = PropertyText::new();
    for (index, fragment) in fragments.iter().enumerate() {
        if index != 0 {
            text.push_byte(b'.');
        }
        text.extend_bytes(fragment.text.as_bytes());
    }
    let atom_order = fragments
        .iter()
        .flat_map(|fragment| fragment.atom_order.iter().copied())
        .collect();
    let bond_order = fragments
        .iter()
        .flat_map(|fragment| fragment.bond_order.iter().copied())
        .collect();
    Ok(SmilesWriteOutput {
        text,
        atom_order,
        bond_order,
    })
}

// Carry the writer's actual per-atom caches into the detached component.
// Source calcImplicitValence only rebuilds explicit valence when its cache is
// negative; a cold Kekulize call would discard retained donor/H-bond caches.
fn kekulize_writer_fragment(
    fragment_topology: &TopologyBlock,
    writer_valence: &ValenceAssignment,
    source_atoms: &[AtomId],
) -> Result<cosmolkit_core::KekulizeAssignment, SmilesParseError> {
    // BEGIN RDKIT CPP BLOCK Atom::calcImplicitValence retained explicit cache
    // RDKit✔️❌:   if (d_explicitValence == -1) {
    // RDKit✔️❌:     calcExplicitValence(strict);
    // RDKit✔️❌:   }
    // END RDKIT CPP BLOCK Atom::calcImplicitValence retained explicit cache
    if source_atoms.len() != fragment_topology.atoms.len()
        || writer_valence.explicit_valence.len() != writer_valence.implicit_hydrogens.len()
        || source_atoms
            .iter()
            .any(|atom| atom.index() >= writer_valence.explicit_valence.len())
    {
        return Err(SmilesParseError::WriterValence(
            "writer fragment valence dimensions or source mapping are invalid".into(),
        ));
    }
    // This detached transport allocates two component-sized Vecs. Core also
    // clones the supplied cache for its owned result; upstream mutates its
    // existing Atom caches. Preserve the behavioral contract explicitly,
    // without claiming allocation equivalence or using a cold-cache fallback.
    let component_valence = ValenceAssignment {
        explicit_valence: source_atoms
            .iter()
            .map(|atom| writer_valence.explicit_valence[atom.index()])
            .collect(),
        implicit_hydrogens: source_atoms
            .iter()
            .map(|atom| writer_valence.implicit_hydrogens[atom.index()])
            .collect(),
    };
    cosmolkit_core::kekulize_with_query_state_and_ring_info(
        fragment_topology,
        &KekulizeParams::default(),
        None,
        None,
        Some(&component_valence),
    )
    .map_err(SmilesParseError::WriterKekulize)
}

// RDKit updates the same fragment Atom objects consumed by GetAtomSmiles.
// Detached Rust atoms and their two cache rows must travel together through the
// component-local -> source-atom mapping. Copy every returned row, including
// rows that were not the neutral aromatic N/P normalization site.
fn project_kekulized_writer_valence(
    writer_valence: &mut ValenceAssignment,
    fragment_valence: Option<&ValenceAssignment>,
    source_atoms: &[AtomId],
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP BLOCK KekulizeFragment cached hydrogen transfer
    // RDKit✔️✔️:           atom->setNoImplicit(false);
    // RDKit✔️✔️:           atom->setNumExplicitHs(0);
    // RDKit✔️✔️:           atom->updatePropertyCache(false);
    // RDKit✔️✔️: unsigned int Atom::getTotalNumHs(bool includeNeighbors) const {
    // RDKit✔️✔️:   int res = getNumExplicitHs() + getNumImplicitHs();
    // END RDKIT CPP BLOCK KekulizeFragment cached hydrogen transfer
    // The core assignment already contains the actual post-Kekulize cache;
    // recomputing it here would change the source's retained-cache semantics.
    // This adds two scalar stores to the existing O(V) detached fragment
    // projection and allocates no cache, mask, or fallback molecule.
    let fragment_valence = fragment_valence.ok_or_else(|| {
        SmilesParseError::WriterValence(
            "nonempty kekulized writer fragment returned no valence cache".into(),
        )
    })?;
    if fragment_valence.explicit_valence.len() != source_atoms.len()
        || fragment_valence.implicit_hydrogens.len() != source_atoms.len()
        || writer_valence.explicit_valence.len() != writer_valence.implicit_hydrogens.len()
        || source_atoms
            .iter()
            .any(|atom| atom.index() >= writer_valence.explicit_valence.len())
    {
        return Err(SmilesParseError::WriterValence(
            "kekulized writer fragment valence dimensions or source mapping are invalid".into(),
        ));
    }
    for (local_index, source_atom) in source_atoms.iter().copied().enumerate() {
        writer_valence.explicit_valence[source_atom.index()] =
            fragment_valence.explicit_valence[local_index];
        writer_valence.implicit_hydrogens[source_atom.index()] =
            fragment_valence.implicit_hydrogens[local_index];
    }
    Ok(())
}

fn prepare_writer_stereochemistry(
    topology: &mut TopologyBlock,
    stereochem_done_marker_is_computed: Option<bool>,
    source_rings: Option<&cosmolkit_core::RingInfo>,
    params: &SmilesWriteParams,
    components: &[Vec<usize>],
) -> Result<(), SmilesParseError> {
    // Keep malformed caller-supplied ring-relative state observable even when
    // clean stereo would clear computed ring relations before serialization.
    for atom in &topology.atoms {
        if let Some(encoded) = atom.prop("_ringStereoAtoms") {
            parse_ring_stereo_atoms(encoded, topology.atoms.len())?;
        }
    }
    if !params.isomeric_smiles {
        return Ok(());
    }

    // BEGIN RDKIT CPP FUNCTION MolOps::assignStereochemistry defaults
    // RDKit❗✔️: RDKIT_GRAPHMOL_EXPORT void assignStereochemistry(
    // RDKit❗✔️:     ROMol &mol, bool cleanIt = false, bool force = false,
    // RDKit❗✔️:     bool flagPossibleStereoCenters = false);
    // END RDKIT CPP FUNCTION MolOps::assignStereochemistry defaults
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles stereo preparation
    // RDKit❗✔️:     if (params.doIsomericSmiles) {
    // RDKit❗✔️:       if (!tmol->hasProp(common_properties::_StereochemDone)) {
    // RDKit❗✔️:         MolOps::assignStereochemistry(*tmol, params.cleanStereo);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles stereo preparation
    // Behavior review: `clean_it=false` uses the existing false/false wrapper;
    // `clean_it=true` currently uses true/true and fails the preserved pinned
    // modern large-ring case. No wrapper is selected by molecule contents.
    prepare_writer_stereo_components(
        topology,
        stereochem_done_marker_is_computed,
        source_rings,
        params.clean_stereo,
        components,
    )
}

fn prepare_canonical_nonisomeric_stereo_fallback(
    topology: &mut TopologyBlock,
    stereochem_done_marker_is_computed: Option<bool>,
    source_rings: Option<&cosmolkit_core::RingInfo>,
    components: &[Vec<usize>],
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION Canon::canonicalizeFragment stereo fallback
    // RDKit❗✔️:   if (!mol.hasProp(common_properties::_StereochemDone)) {
    // RDKit❗✔️:     MolOps::assignStereochemistry(mol, false);
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION Canon::canonicalizeFragment stereo fallback
    // The caller reaches this only after the writer has ranked components and
    // performed FragmentSmilesConstruct's optional kekulization.
    // Complexity review: this reuses the writer's existing component mapping
    // and core assignment; the detached adapter still materializes fragment
    // topology and mapping vectors that RDKit already owns as mutable fragments.
    prepare_writer_stereo_components(
        topology,
        stereochem_done_marker_is_computed,
        source_rings,
        false,
        components,
    )
}

fn prepare_writer_stereo_components(
    topology: &mut TopologyBlock,
    stereochem_done_marker_is_computed: Option<bool>,
    source_rings: Option<&cosmolkit_core::RingInfo>,
    clean_it: bool,
    components: &[Vec<usize>],
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::getTheFrags computed-property copy selection
    // RDKit❗❌:       if (comp.size() == 1 ||
    // RDKit❗❌:           (nFrags > 3 && !fragmentHasChallengingFeatures(comp, atomsInFrag))) {
    // RDKit❗❌:         SubsetOptions opts{.sanitize = sanitizeFrags,
    // RDKit❗❌:                            .clearComputedProps = true,
    // RDKit❗❌:                            .copyCoordinates = copyConformers,
    // RDKit❗❌:                            .method = SubsetMethod::BONDS_BETWEEN_ATOMS};
    // RDKit❗❌:         auto submol = copyMolSubset(mol, atoms, info, opts);
    // RDKit❗❌:         res.push_back(std::move(submol));
    // RDKit❗❌:       } else {
    // RDKit❗❌:         res.emplace_back(new RWMol(mol));
    // RDKit❗❌:         auto &frag = res.back();
    // RDKit❗❌:         frag->beginBatchEdit();
    // RDKit❗❌:         for (unsigned int idx = 0; idx < mol.getNumAtoms(); ++idx) {
    // RDKit❗❌:           if (!atomsInFrag[idx]) {
    // RDKit❗❌:             frag->removeAtom(idx);
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:         frag->commitBatchEdit();
    // RDKit❗❌:       }
    // END RDKIT CPP FUNCTION MolOps::getTheFrags computed-property copy selection
    // BEGIN RDKIT CPP FUNCTION Subset::copyMolSubset molecule-property transport
    // RDKit❗❌:   auto extracted_mol = std::make_unique<RWMol>();
    // RDKit❗❌:   if (options.clearComputedProps) {
    // RDKit❗❌:     extracted_mol->clearComputedProps();
    // RDKit❗❌:   } else {
    // RDKit❗❌:     copyComputedProps(mol, *extracted_mol);
    // RDKit❗❌:   }
    // END RDKIT CPP FUNCTION Subset::copyMolSubset molecule-property transport
    // BEGIN RDKIT CPP FUNCTION RWMol::commitBatchEdit computed-property clearing
    // RDKit❗❌:   batchRemoveBonds();
    // RDKit❗❌:   batchRemoveAtoms();
    // RDKit❗❌:   dp_ringInfo->reset();
    // RDKit❗❌:   clearComputedProps(true);
    // RDKit❗❌:   dp_delBonds.reset();
    // RDKit❗❌:   dp_delAtoms.reset();
    // END RDKIT CPP FUNCTION RWMol::commitBatchEdit computed-property clearing
    // BEGIN RDKIT CPP FUNCTION ROMol::clearComputedProps
    // RDKit❗❌:   RDProps::clearComputedProps();
    // RDKit❗❌:   for (auto atom : atoms()) {
    // RDKit❗❌:     atom->clearComputedProps();
    // RDKit❗❌:   }
    // RDKit❗❌:   for (auto bond : bonds()) {
    // RDKit❗❌:     bond->clearComputedProps();
    // RDKit❗❌:   }
    // END RDKIT CPP FUNCTION ROMol::clearComputedProps
    // BEGIN RDKIT CPP FUNCTION RDProps::clearComputedProps
    // RDKit❗❌:   STR_VECT compLst;
    // RDKit❗❌:   if (getPropIfPresent(RDKit::detail::computedPropName, compLst) &&
    // RDKit❗❌:       !compLst.empty()) {
    // RDKit❗❌:     for (const auto &sv : compLst) {
    // RDKit❗❌:       d_props.clearVal(sv);
    // RDKit❗❌:     }
    // RDKit❗❌:     compLst.clear();
    // RDKit❗❌:     d_props.setVal(RDKit::detail::computedPropName, compLst);
    // RDKit❗❌:   }
    // END RDKIT CPP FUNCTION RDProps::clearComputedProps
    // Behavior review: the single-component clone preserves the marker;
    // subset extraction starts from an empty molecule, while clone/prune
    // clears computed but retains ordinary molecule properties at batch
    // commit. The source guard survives only when the extracted fragment still
    // has `_StereochemDone`.
    // Complexity review: both implementations scan the source graph per
    // fragment. The Rust preserving branch clones source and working topology,
    // then materializes mapped output rows, while RDKit clones once and prunes
    // in place; that branch has a clear extra full-topology copy.
    for (component_index, component) in components.iter().enumerate() {
        let mut inside = vec![false; topology.atoms.len()];
        for atom in component {
            inside[*atom] = true;
        }
        let clears_computed_props = components.len() > 1
            && (component.len() == 1
                || (components.len() > 3
                    && !fragment_has_challenging_features(topology, component, &inside)));
        let marker_survives_extraction = components.len() == 1
            || (!clears_computed_props && stereochem_done_marker_is_computed == Some(false));
        if stereochem_done_marker_is_computed.is_some() && marker_survives_extraction {
            continue;
        }

        let fragment = if components.len() == 1 {
            // BEGIN RDKIT CPP FUNCTION MolOps::getTheFrags one-component path
            // RDKit❗❌:   if (nFrags == 1) {
            // RDKit❗❌:     res.emplace_back(new RWMol(mol));
            // RDKit❗❌:     if (fragsMolAtomMapping) {
            // RDKit❗❌:       INT_VECT comp;
            // RDKit❗❌:       for (unsigned int idx = 0; idx < mol.getNumAtoms(); ++idx) {
            // RDKit❗❌:         comp.push_back(idx);
            // RDKit❗❌:       }
            // RDKit❗❌:       (*fragsMolAtomMapping).push_back(comp);
            // RDKit❗❌:     }
            // RDKit❗❌:   }
            // END RDKIT CPP FUNCTION MolOps::getTheFrags one-component path
            // The source deep-clones the sole component and maps each atom to
            // its unchanged source index. This topology clone also retains
            // ordinary/computed atom and bond properties; the writer's
            // detached adapter adds an atom/bond transport map and merge pass.
            WriterStereoFragment {
                topology: topology.clone(),
                source_atoms: (0..topology.atoms.len()).map(AtomId::new).collect(),
                source_bonds: (0..topology.bonds.len()).map(BondId::new).collect(),
            }
        } else if clears_computed_props {
            extract_writer_subset_fragment(topology, component, &inside)?
        } else {
            extract_writer_preserving_fragment(topology, &inside)?
        };

        let valence = cosmolkit_core::assign_valence_with_options_for_topology(
            &fragment.topology,
            ValenceModel::RdkitLike,
            false,
        )
        .map_err(|error| SmilesParseError::WriterValence(error.to_string()))?;
        // MolOps::getTheFrags preserves the sole component's RingInfo:
        // RDKit✔️✔️:   if (nFrags == 1) {
        // RDKit✔️✔️:     res.emplace_back(new RWMol(mol));
        // legacyStereoPerception must not replace a Fast-or-better cache:
        // RDKit✔️✔️:   if (!mol.getRingInfo()->isFindFastOrBetter()) {
        // RDKit✔️✔️:     MolOps::fastFindRings(mol);
        // RDKit✔️✔️:   }
        // Behavior: single-component indices are unchanged. Extracted/pruned
        // components reset RingInfo in the source and use the normal fallback.
        // The core owner handles an uninitialized/Other input cache itself.
        // Complexity: reuse a borrowed cache in O(1), without a clone or search;
        // absent/extracted caches retain the existing O(V+E) fast search.
        let fallback_rings;
        let rings = if let Some(rings) = source_rings.filter(|_| components.len() == 1) {
            rings
        } else {
            fallback_rings = fast_find_rings_from_parts(
                fragment.topology.atoms.len(),
                &fragment.topology.bonds,
                &fragment.topology.adjacency,
            )
            .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
            &fallback_rings
        };
        let assigned = if clean_it {
            // RDKit❗✔️: the writer passes cleanStereo while the omitted
            // flagPossibleStereoCenters argument defaults to false; the current
            // cleanup wrapper selects true for that flag. The serialized-output
            // regression below checks this wrapper independently.
            cosmolkit_core::assign_legacy_stereochemistry(fragment.topology, &valence, rings)
        } else {
            cosmolkit_core::assign_legacy_stereochemistry_for_depiction(
                fragment.topology,
                &valence,
                rings,
            )
        }
        .map_err(|error| match error {
            cosmolkit_core::LegacyStereoError::StereoGroup(cause) => {
                SmilesParseError::StereoGroup(cause)
            }
            other => SmilesParseError::WriterStereo(other.to_string()),
        })?;

        merge_writer_stereo_fragment(
            topology,
            &assigned,
            &fragment.source_atoms,
            &fragment.source_bonds,
        )?;

        // Keep the source's component walk visible in this local anchor. The
        // index is intentionally consumed by debug assertions rather than
        // changing the serialized traversal order.
        debug_assert!(component_index < components.len());
    }
    Ok(())
}

fn fragment_has_challenging_features(
    topology: &TopologyBlock,
    component: &[usize],
    inside: &[bool],
) -> bool {
    // BEGIN RDKIT CPP FUNCTION MolOps::getTheFrags fragmentHasChallengingFeatures
    // RDKit✔️✔️:       auto fragmentHasChallengingFeatures =
    // RDKit✔️✔️:           [&](const INT_VECT &comp,
    // RDKit✔️✔️:               const boost::dynamic_bitset<> &atomsInFrag) -> bool {
    // RDKit✔️✔️:         for (auto idx : comp) {
    // RDKit✔️✔️:           const auto atom = mol.getAtomWithIdx(idx);
    // RDKit✔️✔️:           if (atom->getChiralTag() != Atom::ChiralType::CHI_UNSPECIFIED &&
    // RDKit✔️✔️:               atom->getChiralTag() != Atom::ChiralType::CHI_OTHER) {
    // RDKit✔️✔️:             return true;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           for (auto bnd : mol.atomBonds(atom)) {
    // RDKit✔️✔️:             if (atomsInFrag[bnd->getOtherAtomIdx(idx)]) {
    // RDKit✔️✔️:               if (bnd->getStereo() != Bond::BondStereo::STEREONONE &&
    // RDKit✔️✔️:                   bnd->getStereo() != Bond::BondStereo::STEREOANY) {
    // RDKit✔️✔️:                 return true;
    // RDKit✔️✔️:               }
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         for (auto sgroup : getSubstanceGroups(mol)) {
    // RDKit✔️✔️:           for (auto aid : sgroup.getAtoms()) {
    // RDKit✔️✔️:             if (atomsInFrag[aid]) {
    // RDKit✔️✔️:               return true;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           for (auto aid : sgroup.getParentAtoms()) {
    // RDKit✔️✔️:             if (atomsInFrag[aid]) {
    // RDKit✔️✔️:               return true;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         for (auto stereoGroup : mol.getStereoGroups()) {
    // RDKit✔️✔️:           for (auto atom : stereoGroup.getAtoms()) {
    // RDKit✔️✔️:             if (atomsInFrag[atom->getIdx()]) {
    // RDKit✔️✔️:               return true;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           for (auto bond : stereoGroup.getBonds()) {
    // RDKit✔️✔️:             if (atomsInFrag[bond->getBeginAtomIdx()] &&
    // RDKit✔️✔️:                 atomsInFrag[bond->getEndAtomIdx()]) {
    // RDKit✔️✔️:               return true;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       };
    // END RDKIT CPP FUNCTION MolOps::getTheFrags fragmentHasChallengingFeatures
    // Behavior review: each modeled challenge category is checked on the
    // selected component. The source and Rust loops inspect atom tags,
    // internal bond stereo, SGroup atom/parent membership, and stereo-group
    // atom plus fully internal bond membership.
    // Complexity review: one component traversal and source-like scans of
    // incident bonds and group memberships; no nested whole-graph rescan is
    // added beyond the source lambda.
    for atom_index in component {
        let atom = &topology.atoms[*atom_index];
        if !matches!(atom.chiral_tag(), ChiralTag::Unspecified | ChiralTag::Other) {
            return true;
        }
        for neighbor in topology.adjacency.neighbors_of(*atom_index) {
            if inside[neighbor.atom_index] {
                let stereo = topology.bonds[neighbor.bond.index()].stereo();
                if !matches!(stereo, BondStereo::None | BondStereo::Any) {
                    return true;
                }
            }
        }
    }
    if topology.substance_groups.iter().any(|group| {
        group
            .atoms()
            .iter()
            .chain(group.parent_atoms())
            .any(|atom| inside[atom.index()])
    }) {
        return true;
    }
    topology.stereo_groups.iter().any(|group| {
        group.atoms().iter().any(|atom| inside[atom.index()])
            || group.bonds().iter().any(|bond| {
                let bond = &topology.bonds[bond.index()];
                inside[bond.begin().index()] && inside[bond.end().index()]
            })
    })
}

fn writer_component_uses_subset_stereo_copy(
    topology: &TopologyBlock,
    component: &[usize],
    component_count: usize,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION MolOps::getTheFrags component copy dispatch
    // RDKit❗✔️:   if (nFrags == 1) {
    // RDKit❗✔️:     res.emplace_back(new RWMol(mol));
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     if (comp.size() == 1 ||
    // RDKit❗✔️:         (nFrags > 3 && !fragmentHasChallengingFeatures(comp, atomsInFrag))) {
    // RDKit❗✔️:       auto submol = copyMolSubset(mol, atoms, info, opts);
    // RDKit❗✔️:       res.push_back(std::move(submol));
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       res.emplace_back(new RWMol(mol));
    // RDKit❗✔️:       auto &frag = res.back();
    // RDKit❗✔️:       frag->beginBatchEdit();
    // RDKit❗✔️:       for (unsigned int idx = 0; idx < mol.getNumAtoms(); ++idx) {
    // RDKit❗✔️:         if (!atomsInFrag[idx]) {
    // RDKit❗✔️:           frag->removeAtom(idx);
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:       frag->commitBatchEdit();
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION MolOps::getTheFrags component copy dispatch
    // Behavior review: the sole-fragment path is a full copy; singleton
    // fragments and sufficiently numerous uncomplicated fragments use subset
    // copying, while other multi-fragment cases clone and prune.
    // Complexity review: the branch decision makes one V-sized membership
    // vector and reuses the existing feature scan, matching the source's
    // component-local O(V+E) decision work without constructing a topology.
    if component_count == 1 {
        return false;
    }
    let mut inside = vec![false; topology.atoms.len()];
    for atom in component {
        inside[*atom] = true;
    }
    component.len() == 1
        || (component_count > 3 && !fragment_has_challenging_features(topology, component, &inside))
}

fn writer_component_stereo_groups(
    topology: &TopologyBlock,
    component: &[usize],
    subset_copy: bool,
) -> Result<Vec<StereoGroup>, SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION Subset::copySelectedStereoGroups selection
    // RDKit❗✔️:   auto is_selected_component = [](auto &objects, auto &selected_indices) {
    // RDKit❗✔️:     return objects.empty() ||
    // RDKit❗✔️:            std::any_of(objects.begin(), objects.end(), [&](auto &object) {
    // RDKit❗✔️:              return selected_indices[object->getIdx()];
    // RDKit❗✔️:            });
    // RDKit❗✔️:   };
    // RDKit❗✔️:   auto is_selected_stereo_group = [&](const auto &stereo_group) {
    // RDKit❗✔️:     return is_selected_component(stereo_group.getAtoms(),
    // RDKit❗✔️:                                  selection_info.selectedAtoms) &&
    // RDKit❗✔️:            is_selected_component(stereo_group.getBonds(),
    // RDKit❗✔️:                                  selection_info.selectedBonds);
    // RDKit❗✔️:   };
    // RDKit❗✔️:   std::vector<Atom *> extracted_atoms(extracted_mol.getNumAtoms());
    // RDKit❗✔️:   for (const auto &atom : extracted_mol.atoms()) {
    // RDKit❗✔️:     extracted_atoms[atom->getIdx()] = atom;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   std::vector<Bond *> extracted_bonds(extracted_mol.getNumBonds());
    // RDKit❗✔️:   for (const auto &bond : extracted_mol.bonds()) {
    // RDKit❗✔️:     extracted_bonds[bond->getIdx()] = bond;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   const auto &[selectedAtoms, selectedBonds, atomMapping, bondMapping] =
    // RDKit❗✔️:       selection_info;
    // RDKit❗✔️:   std::vector<StereoGroup> extracted_stereo_groups;
    // RDKit❗✔️:   for (const auto &stereo_group : reference_mol.getStereoGroups()) {
    // RDKit❗✔️:     if (!is_selected_stereo_group(stereo_group)) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     std::vector<Atom *> atoms;
    // RDKit❗✔️:     for (const auto &atom : stereo_group.getAtoms()) {
    // RDKit❗✔️:       auto mapping = atomMapping.find(atom->getIdx());
    // RDKit❗✔️:       if (mapping != atomMapping.end()) {
    // RDKit❗✔️:         atoms.push_back(extracted_atoms[mapping->second]);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     std::vector<Bond *> bonds;
    // RDKit❗✔️:     for (const auto &bond : stereo_group.getBonds()) {
    // RDKit❗✔️:       auto mapping = bondMapping.find(bond->getIdx());
    // RDKit❗✔️:       if (mapping != bondMapping.end()) {
    // RDKit❗✔️:         bonds.push_back(extracted_bonds[mapping->second]);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     extracted_stereo_groups.push_back({stereo_group.getGroupType(),
    // RDKit❗✔️:                                        std::move(atoms), std::move(bonds),
    // RDKit❗✔️:                                        stereo_group.getReadId()});
    // RDKit❗✔️:     extracted_stereo_groups.back().setWriteId(stereo_group.getWriteId());
    // RDKit❗✔️:   }
    // RDKit❗✔️:   extracted_mol.setStereoGroups(std::move(extracted_stereo_groups));
    // END RDKIT CPP FUNCTION Subset::copySelectedStereoGroups selection
    // BEGIN RDKIT CPP FUNCTION StereoGroup removeAtomFromGroups/removeBondFromGroups
    // RDKit❗✔️: void removeAtomFromGroups(const Atom *atom, std::vector<StereoGroup> &groups) {
    // RDKit❗✔️:   auto findAtom = [atom](StereoGroup &group) {
    // RDKit❗✔️:     return std::find(group.getAtoms().begin(), group.getAtoms().end(), atom);
    // RDKit❗✔️:   };
    // RDKit❗✔️:   for (auto &group : groups) {
    // RDKit❗✔️:     auto atomPos = findAtom(group);
    // RDKit❗✔️:     if (atomPos != group.d_atoms.end()) {
    // RDKit❗✔️:       group.d_atoms.erase(atomPos);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   groups.erase(std::remove_if(groups.begin(), groups.end(),
    // RDKit❗✔️:                               [](const auto &gp) {
    // RDKit❗✔️:                                 return gp.getAtoms().empty() &&
    // RDKit❗✔️:                                        gp.getBonds().empty();
    // RDKit❗✔️:                               }),
    // RDKit❗✔️:                groups.end());
    // RDKit❗✔️: }
    // RDKit❗✔️: void removeBondFromGroups(const Bond *bond, std::vector<StereoGroup> &groups) {
    // RDKit❗✔️:   auto findBond = [bond](StereoGroup &group) {
    // RDKit❗✔️:     return std::find(group.getBonds().begin(), group.getBonds().end(), bond);
    // RDKit❗✔️:   };
    // RDKit❗✔️:   for (auto &group : groups) {
    // RDKit❗✔️:     auto bondPos = findBond(group);
    // RDKit❗✔️:     if (bondPos != group.d_bonds.end()) {
    // RDKit❗✔️:       group.d_bonds.erase(bondPos);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   groups.erase(std::remove_if(groups.begin(), groups.end(),
    // RDKit❗✔️:                               [](const auto &gp) {
    // RDKit❗✔️:                                 return gp.getAtoms().empty() &&
    // RDKit❗✔️:                                        gp.getBonds().empty();
    // RDKit❗✔️:                               }),
    // RDKit❗✔️:                groups.end());
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION StereoGroup removeAtomFromGroups/removeBondFromGroups
    // Behavior review: source order, read IDs and write IDs are retained. Clone
    // pruning retains a group if any member survives; subset copying retains a
    // group only when every nonempty member category overlaps the selected set.
    // Pre-existing empty groups survive both paths. Component bonds are selected
    // exactly when both full-topology endpoints belong to this component.
    // Complexity review: one atom-membership mask and one pass over groups and
    // their members, with O(G+A+B) output allocation; source extraction also
    // scans membership rows and copies each retained group member once.
    let mut inside = vec![false; topology.atoms.len()];
    for atom in component {
        inside[*atom] = true;
    }

    topology
        .stereo_groups
        .iter()
        .map(|group| {
            let atoms = group
                .atoms()
                .iter()
                .copied()
                .filter(|atom| inside[atom.index()])
                .collect::<Vec<_>>();
            let bonds = group
                .bonds()
                .iter()
                .copied()
                .filter(|bond_id| {
                    let bond = &topology.bonds[bond_id.index()];
                    inside[bond.begin().index()] && inside[bond.end().index()]
                })
                .collect::<Vec<_>>();
            let atoms_selected = !atoms.is_empty();
            let bonds_selected = !bonds.is_empty();
            let keep_group = if subset_copy {
                (group.atoms().is_empty() || atoms_selected)
                    && (group.bonds().is_empty() || bonds_selected)
            } else {
                (group.atoms().is_empty() && group.bonds().is_empty())
                    || atoms_selected
                    || bonds_selected
            };
            if !keep_group {
                return Ok(None);
            }

            let mut selected_group = StereoGroup::new(group.kind(), atoms, bonds)?;
            if let Some(read_id) = group.id() {
                selected_group = selected_group.with_id(read_id);
            }
            set_stereo_group_write_id(&mut selected_group, stereo_group_write_id(group));
            Ok(Some(selected_group))
        })
        .collect::<Result<Vec<_>, SmilesParseError>>()
        .map(|rows| rows.into_iter().flatten().collect())
}

fn extract_writer_subset_fragment(
    topology: &TopologyBlock,
    component: &[usize],
    inside: &[bool],
) -> Result<WriterStereoFragment, SmilesParseError> {
    // Behavior review: for a multi-atom connected component, `component_bonds`
    // below contains every bond with both endpoints in `inside`; this is the
    // source `BONDS_BETWEEN_ATOMS` selection. `subtopology_from_path` creates
    // source-ordered rows, clears atom/bond computed properties, clears cis or
    // trans stereo when either stereo reference is unmappable, and returns
    // validated `new_to_old` maps. The singleton branch below supplies the
    // no-bond case directly and copies source-selected groups. The model has
    // additional typed SGroup references beyond RDKit's three selected index
    // lists, so those edge cases remain explicitly incomplete.
    // Complexity review: each component allocates a V-sized mask, scans all
    // source bonds to collect its internal edges, then the core subset path
    // scans atom and bond rows and allocates old/new maps. RDKit's
    // `copyMolSubset` consumes selected bitsets and scans its source rows once;
    // the local edge-inventory pass and mapping vectors are additional work.
    if component.len() == 1 {
        // BEGIN RDKIT CPP FUNCTION Subset::copySelectedAtomsAndBonds singleton
        // RDKit❗❗:   for (const auto &ref_atom : reference_mol.atoms()) {
        // RDKit❗❗:     if (!selectedAtoms[ref_atom->getIdx()]) {
        // RDKit❗❗:       continue;
        // RDKit❗❗:     }
        // RDKit❗❗:     std::unique_ptr<Atom> extracted_atom{
        // RDKit❗❗:         options.copyAsQuery ? new QueryAtom(*ref_atom) : ref_atom->copy()};
        // RDKit❗❗:     extracted_atom->clearComputedProps();
        // RDKit❗❗:   }
        // END RDKIT CPP FUNCTION Subset::copySelectedAtomsAndBonds singleton
        // BEGIN RDKIT CPP FUNCTION Subset::isSelectedSGroup/copySelectedSubstanceGroups
        // RDKit❗❗:   auto is_selected_component = [](auto &indices, auto &selection_test) {
        // RDKit❗❗:     return indices.empty() ||
        // RDKit❗❗:            std::all_of(indices.begin(), indices.end(), selection_test);
        // RDKit❗❗:   };
        // RDKit❗❗:   return is_selected_component(sgroup.getAtoms(), atom_test) &&
        // RDKit❗❗:          is_selected_component(sgroup.getBonds(), bond_test) &&
        // RDKit❗❗:          is_selected_component(sgroup.getParentAtoms(), atom_test);
        // RDKit❗❗:     if (!isSelectedSGroup(sgroup, selection_info)) {
        // RDKit❗❗:       continue;
        // RDKit❗❗:     }
        // RDKit❗❗:     SubstanceGroup extracted_sgroup(sgroup);
        // RDKit❗❗:     extracted_sgroup.setOwningMol(&extracted_mol);
        // RDKit❗❗:     update_indices(extracted_sgroup, std::mem_fn(&SubstanceGroup::getAtoms),
        // RDKit❗❗:                    std::mem_fn(&SubstanceGroup::setAtoms), atomMapping);
        // RDKit❗❗:     update_indices(extracted_sgroup,
        // RDKit❗❗:                    std::mem_fn(&SubstanceGroup::getParentAtoms),
        // RDKit❗❗:                    std::mem_fn(&SubstanceGroup::setParentAtoms), atomMapping);
        // RDKit❗❗:     update_indices(extracted_sgroup, std::mem_fn(&SubstanceGroup::getBonds),
        // RDKit❗❗:                    std::mem_fn(&SubstanceGroup::setBonds), bondMapping);
        // RDKit❗❗:     addSubstanceGroup(extracted_mol, std::move(extracted_sgroup));
        // END RDKIT CPP FUNCTION Subset::isSelectedSGroup/copySelectedSubstanceGroups
        // BEGIN RDKIT CPP FUNCTION Subset::copySelectedStereoGroups
        // RDKit❗❗:   auto is_selected_component = [](auto &objects, auto &selected_indices) {
        // RDKit❗❗:     return objects.empty() ||
        // RDKit❗❗:            std::any_of(objects.begin(), objects.end(), [&](auto &object) {
        // RDKit❗❗:              return selected_indices[object->getIdx()];
        // RDKit❗❗:            });
        // RDKit❗❗:   };
        // RDKit❗❗:   auto is_selected_stereo_group = [&](const auto &stereo_group) {
        // RDKit❗❗:     return is_selected_component(stereo_group.getAtoms(),
        // RDKit❗❗:                                  selection_info.selectedAtoms) &&
        // RDKit❗❗:            is_selected_component(stereo_group.getBonds(),
        // RDKit❗❗:                                  selection_info.selectedBonds);
        // RDKit❗❗:   };
        // RDKit❗❗:   if (!is_selected_stereo_group(stereo_group)) { continue; }
        // END RDKIT CPP FUNCTION Subset::copySelectedStereoGroups
        // BEGIN RDKIT CPP FUNCTION Subset::copyMolSubset computed properties
        // RDKit❗❗:   if (options.clearComputedProps) {
        // RDKit❗❗:     extracted_mol->clearComputedProps();
        // RDKit❗❗:   } else {
        // RDKit❗❗:     copyComputedProps(mol, *extracted_mol);
        // RDKit❗❗:   }
        // END RDKIT CPP FUNCTION Subset::copyMolSubset computed properties
        let source_atom = AtomId::new(component[0]);
        let mut atom_old_to_new = vec![None; topology.atoms.len()];
        atom_old_to_new[source_atom.index()] = Some(AtomId::new(0));
        let mut atom = topology.atoms[source_atom.index()]
            .clone()
            .with_id(AtomId::new(0));
        atom.clear_computed_props()?;
        atom.remap_template_attachment_order(&atom_old_to_new)
            .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
        let stereo_groups = topology
            .stereo_groups
            .iter()
            .map(|group| {
                let atom_side_selected = group.atoms().is_empty()
                    || group.atoms().iter().any(|member| inside[member.index()]);
                let bond_side_selected = group.bonds().is_empty();
                if !atom_side_selected || !bond_side_selected {
                    return Ok(None);
                }
                let atoms = group
                    .atoms()
                    .iter()
                    .filter(|member| inside[member.index()])
                    .map(|_| AtomId::new(0))
                    .collect();
                let local = StereoGroup::new(group.kind(), atoms, Vec::new())?;
                Ok(Some(match group.id() {
                    Some(id) => local.with_id(id),
                    None => local,
                }))
            })
            .collect::<Result<Vec<_>, SmilesParseError>>()?
            .into_iter()
            .flatten()
            .collect();
        let bond_old_to_new = vec![None; topology.bonds.len()];
        let mut selected_sgroups = topology
            .substance_groups
            .iter()
            .map(|group| {
                let source_members_selected = group
                    .atoms()
                    .iter()
                    .all(|member| atom_old_to_new[member.index()].is_some())
                    && group
                        .bonds()
                        .iter()
                        .all(|member| bond_old_to_new[member.index()].is_some())
                    && group
                        .parent_atoms()
                        .iter()
                        .all(|member| atom_old_to_new[member.index()].is_some());
                source_members_selected
                    && group.can_remap_without_parent(&atom_old_to_new, &bond_old_to_new)
            })
            .collect::<Vec<_>>();
        loop {
            let mut changed = false;
            for (index, group) in topology.substance_groups.iter().enumerate() {
                if selected_sgroups[index]
                    && group.parent().is_some_and(|parent| {
                        !selected_sgroups
                            .get(parent.index())
                            .copied()
                            .unwrap_or(false)
                    })
                {
                    selected_sgroups[index] = false;
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
        let mut sgroup_old_to_new = vec![None; topology.substance_groups.len()];
        let mut next_sgroup_index = 0;
        for (old_index, keep) in selected_sgroups.iter().copied().enumerate() {
            if keep {
                sgroup_old_to_new[old_index] =
                    Some(cosmolkit_model::SubstanceGroupId::new(next_sgroup_index));
                next_sgroup_index += 1;
            }
        }
        let substance_groups = topology
            .substance_groups
            .iter()
            .enumerate()
            .filter(|(index, _)| selected_sgroups[*index])
            .map(|(index, group)| {
                group
                    .remapped(
                        sgroup_old_to_new[index].expect("selected SGroup has a new id"),
                        &atom_old_to_new,
                        &bond_old_to_new,
                        &sgroup_old_to_new,
                    )
                    .ok_or_else(|| {
                        SmilesParseError::WriterStereo(
                            "selected singleton SGroup could not be remapped".into(),
                        )
                    })
            })
            .collect::<Result<Vec<_>, _>>()?;
        // The model's typed SGroup state can reference attach/crossing atoms
        // and CState bonds outside the three source membership lists. Such a
        // group is excluded above because its references cannot survive this
        // subset. Source-selected memberships and all representable references
        // are remapped into the singleton row before topology validation.
        let fragment =
            TopologyBlock::try_from_parts(vec![atom], Vec::new(), substance_groups, stereo_groups)
                .map_err(|error| match error {
                    cosmolkit_model::TopologyValidationError::StereoGroup(cause) => {
                        SmilesParseError::StereoGroup(cause)
                    }
                    other => SmilesParseError::WriterStereo(other.to_string()),
                })?;
        return Ok(WriterStereoFragment {
            topology: fragment,
            source_atoms: vec![source_atom],
            source_bonds: Vec::new(),
        });
    }

    // BEGIN RDKIT CPP FUNCTION Subset::getSubsetInfo BONDS_BETWEEN_ATOMS
    // RDKit✔️❌:   if (options.method == SubsetMethod::BONDS_BETWEEN_ATOMS) {
    // RDKit✔️❌:     for (const auto &atom_idx : path) {
    // RDKit✔️❌:       if (atom_idx < num_atoms) {
    // RDKit✔️❌:         selectedAtoms.set(atom_idx);
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:     for (const auto &bond : mol.bonds()) {
    // RDKit✔️❌:       if (selectedAtoms[bond->getBeginAtomIdx()] &&
    // RDKit✔️❌:           selectedAtoms[bond->getEndAtomIdx()]) {
    // RDKit✔️❌:         selectedBonds.set(bond->getIdx());
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // END RDKIT CPP FUNCTION Subset::getSubsetInfo BONDS_BETWEEN_ATOMS
    // BEGIN RDKIT CPP FUNCTION Subset::copySelectedAtomsAndBonds atom old/new rows
    // RDKit✔️❌:   for (const auto &ref_atom : reference_mol.atoms()) {
    // RDKit✔️❌:     if (!selectedAtoms[ref_atom->getIdx()]) {
    // RDKit✔️❌:       continue;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     std::unique_ptr<Atom> extracted_atom{
    // RDKit✔️❌:         options.copyAsQuery ? new QueryAtom(*ref_atom) : ref_atom->copy()};
    // RDKit✔️❌:     extracted_atom->clearComputedProps();
    // RDKit✔️❌:     atomMapping[ref_atom->getIdx()] = extracted_mol.addAtom(
    // RDKit✔️❌:         extracted_atom.release(), updateLabel, takeOwnership);
    // RDKit✔️❌:   }
    // END RDKIT CPP FUNCTION Subset::copySelectedAtomsAndBonds atom old/new rows
    // BEGIN RDKIT CPP FUNCTION Subset::copySelectedAtomsAndBonds bond old/new rows
    // RDKit✔️❌:     auto num_bonds =
    // RDKit✔️❌:         extracted_mol.addBond(extracted_bond.release(), takeOwnership);
    // RDKit✔️❌:     bondMapping[ref_bond->getIdx()] = num_bonds - 1;
    // END RDKIT CPP FUNCTION Subset::copySelectedAtomsAndBonds bond old/new rows
    // BEGIN RDKIT CPP FUNCTION Subset::copySelectedAtomsAndBonds mapped bonds
    // RDKit❗❗:   for (const auto &ref_bond : reference_mol.bonds()) {
    // RDKit❗❗:     if (!selectedBonds[ref_bond->getIdx()]) {
    // RDKit❗❗:       continue;
    // RDKit❗❗:     }
    // RDKit❗❗:     if (atomMapping.find(ref_bond->getBeginAtomIdx()) == atomMapping.end() ||
    // RDKit❗❗:         atomMapping.find(ref_bond->getEndAtomIdx()) == atomMapping.end()) {
    // RDKit❗❗:       throw ValueErrorException("copyMolSubset: subset bonds contain atoms not contained in subset atoms");
    // RDKit❗❗:     }
    // RDKit❗❗:     std::unique_ptr<Bond> extracted_bond{
    // RDKit❗❗:         options.copyAsQuery ? new QueryBond(*ref_bond) : ref_bond->copy()};
    // RDKit❗❗:     auto &atoms = extracted_bond->getStereoAtoms();
    // RDKit❗❗:     if (atoms.size() == 2) {
    // RDKit❗❗:       auto map1 = atomMapping.find(atoms[0]);
    // RDKit❗❗:       auto map2 = atomMapping.find(atoms[1]);
    // RDKit❗❗:       if (map1 != atomMapping.end() && map2 != atomMapping.end()) {
    // RDKit❗❗:         atoms[0] = map1->second;
    // RDKit❗❗:         atoms[1] = map2->second;
    // RDKit❗❗:       } else {
    // RDKit❗❗:         atoms.clear();  // We couldn't map the stereo atoms
    // RDKit❗❗:       }
    // RDKit❗❗:     }
    // RDKit❗❗:     for (auto &atomidx : atoms) {
    // RDKit❗❗:       auto map = atomMapping.find(atomidx);
    // RDKit❗❗:       if (map != atomMapping.end()) {
    // RDKit❗❗:         atomidx = map->second;
    // RDKit❗❗:       }
    // RDKit❗❗:     }
    // RDKit❗❗:   }
    // END RDKIT CPP FUNCTION Subset::copySelectedAtomsAndBonds mapped bonds
    let component_bonds = topology
        .bonds
        .iter()
        .filter(|bond| inside[bond.begin().index()] && inside[bond.end().index()])
        .map(Bond::id)
        .collect::<Vec<_>>();
    let result = cosmolkit_core::subtopology_from_path(
        topology,
        &component_bonds,
        &cosmolkit_core::SubtopologyParams::default(),
    )
    .map_err(|error| match error {
        cosmolkit_core::PathError::StereoGroup(cause) => SmilesParseError::StereoGroup(cause),
        other => SmilesParseError::WriterStereo(other.to_string()),
    })?;
    let cosmolkit_core::DetachedPathSubgraph::Concrete(fragment) = result.subgraph else {
        return Err(SmilesParseError::WriterStereo(
            "writer component extraction returned a query graph".into(),
        ));
    };
    Ok(WriterStereoFragment {
        topology: fragment,
        source_atoms: result
            .mapping
            .atoms
            .new_to_old
            .into_iter()
            .flatten()
            .collect(),
        source_bonds: result
            .mapping
            .bonds
            .new_to_old
            .into_iter()
            .flatten()
            .collect(),
    })
}

fn extract_writer_preserving_fragment(
    topology: &TopologyBlock,
    inside: &[bool],
) -> Result<WriterStereoFragment, SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::getTheFrags clone-and-prune branch
    // RDKit❗❌:         res.emplace_back(new RWMol(mol));
    // RDKit❗❌:         auto &frag = res.back();
    // RDKit❗❌:         frag->beginBatchEdit();
    // RDKit❗❌:         for (unsigned int idx = 0; idx < mol.getNumAtoms(); ++idx) {
    // RDKit❗❌:           if (!atomsInFrag[idx]) {
    // RDKit❗❌:             frag->removeAtom(idx);
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:         frag->commitBatchEdit();
    // END RDKIT CPP FUNCTION MolOps::getTheFrags clone-and-prune branch
    // Behavior review: row order and atom/bond mappings follow source order;
    // ordinary atom/bond properties and typed groups survive the validated
    // mapping. RDKit's `commitBatchEdit` clears computed properties on the
    // molecule, atoms, and bonds, so this topology-only adapter clears the
    // computed atom/bond properties after `finish`. The source removes
    // nonmembers in ascending source atom order and remaps surviving rows.
    // Complexity review: clearing computed properties scans retained atom
    // and bond rows, matching the source's linear property-clear pass;
    // `begin_batch_edit` still clones both source and working blocks, then
    // `finish` creates remapped rows. These two full topology copies are a
    // confirmed local cost increase over RDKit's one-copy in-place branch.
    let mut edit = topology
        .begin_batch_edit()
        .map_err(SmilesParseError::from)?;
    for (index, selected) in inside.iter().copied().enumerate() {
        if !selected {
            edit.remove_atom(AtomId::new(index))
                .map_err(SmilesParseError::from)?;
        }
    }
    let (mut fragment, mapping) = edit.finish().map_err(SmilesParseError::from)?;
    // BEGIN RDKIT CPP FUNCTION RWMol::commitBatchEdit computed-property clearing
    // RDKit❗❌:   // fix properties
    // RDKit❗❌:   clearComputedProps(true);
    // END RDKIT CPP FUNCTION RWMol::commitBatchEdit computed-property clearing
    // BEGIN RDKIT CPP FUNCTION ROMol::clearComputedProps
    // RDKit❗❌:   RDProps::clearComputedProps();
    // RDKit❗❌:   for (auto atom : atoms()) {
    // RDKit❗❌:     atom->clearComputedProps();
    // RDKit❗❌:   }
    // RDKit❗❌:   for (auto bond : bonds()) {
    // RDKit❗❌:     bond->clearComputedProps();
    // RDKit❗❌:   }
    // END RDKIT CPP FUNCTION ROMol::clearComputedProps
    // TopologyBlock has no molecule-level property block or persisted ring
    // cache; those source-owned fields remain handled by their caller. Its
    // retained atom and bond rows carry the source property-clearing effect.
    for atom in &mut fragment.atoms {
        atom.clear_computed_props()?;
    }
    for bond in &mut fragment.bonds {
        bond.clear_computed_props()?;
    }
    Ok(WriterStereoFragment {
        topology: fragment,
        source_atoms: mapping.atoms.new_to_old.into_iter().flatten().collect(),
        source_bonds: mapping.bonds.new_to_old.into_iter().flatten().collect(),
    })
}

fn merge_writer_stereo_fragment(
    topology: &mut TopologyBlock,
    assigned: &TopologyBlock,
    source_atoms: &[AtomId],
    source_bonds: &[BondId],
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles assignment producer
    // RDKit❗❌:       if (!tmol->hasProp(common_properties::_StereochemDone)) {
    // RDKit❗❌:         MolOps::assignStereochemistry(*tmol, params.cleanStereo);
    // RDKit❗❌:       }
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles assignment producer
    // Rust-only detached transport: RDKit mutates each temporary fragment in
    // place and has no merge function. The core assignment returns an owned
    // topology, so this adapter maps each assigned row through the extraction
    // `new_to_old` vectors: atoms carry chiral tag, explicit-H/no-implicit
    // state, and computed properties; bonds carry direction, stereo enum,
    // stereo-atom references, and computed properties. Ordinary target props
    // and the original typed groups stay on the full topology. `_ringStereoAtoms`
    // is the one computed property whose atom indices are fragment-relative;
    // signed one-based local IDs are converted to signed one-based source IDs.
    // D10's regression covers both signs and ordinary/computed property
    // preservation. The extraction mapping is validated before these rows are
    // returned by the subset/batch-edit owners.
    // Complexity review: source mutation has no copy-back pass. This adapter
    // scans each mapped atom/bond and clones computed property strings; that
    // is an extra O(V+E) transport and allocation cost.
    for (local_index, source_id) in source_atoms.iter().copied().enumerate() {
        let prepared = &assigned.atoms[local_index];
        let target = &mut topology.atoms[source_id.index()];
        target.set_chiral_tag(prepared.chiral_tag());
        target.set_explicit_hydrogens(prepared.explicit_hydrogens());
        target.set_no_implicit(prepared.no_implicit());
        target.clear_computed_props()?;
        for (key, value) in prepared.props() {
            if !prepared.is_prop_computed(key)? {
                continue;
            }
            let value = if key.as_bytes() == b"_ringStereoAtoms" {
                parse_ring_stereo_atoms(value, assigned.atoms.len())?
                    .into_iter()
                    .map(|(same_orientation, local_atom)| {
                        let source_atom = source_atoms[local_atom].index() + 1;
                        let signed = i32::try_from(source_atom).map_err(|_| {
                            SmilesParseError::WriterStereo(
                                "`_ringStereoAtoms` index is out of range".into(),
                            )
                        })?;
                        Ok(if same_orientation { signed } else { -signed })
                    })
                    .collect::<Result<Vec<i32>, SmilesParseError>>()?
                    .into()
            } else {
                value.clone()
            };
            target.set_computed_prop(key.clone(), value)?;
        }
    }
    for (local_index, source_id) in source_bonds.iter().copied().enumerate() {
        let prepared = &assigned.bonds[local_index];
        let target = &mut topology.bonds[source_id.index()];
        target.set_direction(prepared.direction());
        let stereo_atoms = prepared.stereo_atoms().map(|atoms| {
            [
                source_atoms[atoms[0].index()],
                source_atoms[atoms[1].index()],
            ]
        });
        target.set_stereo_atoms(stereo_atoms);
        target.set_stereo(prepared.stereo())?;
        target.clear_computed_props()?;
        for (key, value) in prepared.props() {
            if prepared.is_prop_computed(key)? {
                target
                    .set_computed_prop(key.clone(), value.clone())
                    .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
            }
        }
    }
    Ok(())
}

fn connected_components(topology: &TopologyBlock) -> Vec<Vec<usize>> {
    // BEGIN RDKIT CPP FUNCTION MolOps::getMolFrags component mapping
    // RDKit✔️❌: unsigned int getMolFrags(const ROMol &mol, INT_VECT &mapping) {
    // RDKit✔️❌:   unsigned int natms = mol.getNumAtoms();
    // RDKit✔️❌:   mapping.resize(natms);
    // RDKit✔️❌:   return natms ? boost::connected_components(mol.getTopology(), &mapping[0])
    // RDKit✔️❌:                    : 0;
    // RDKit✔️❌: };
    // RDKit✔️❌: unsigned int getMolFrags(const ROMol &mol, VECT_INT_VECT &frags) {
    // RDKit✔️❌:   frags.clear();
    // RDKit✔️❌:   INT_VECT mapping;
    // RDKit✔️❌:   getMolFrags(mol, mapping);
    // RDKit✔️❌:
    // RDKit✔️❌:   INT_INT_VECT_MAP comMap;
    // RDKit✔️❌:   for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
    // RDKit✔️❌:     int mi = mapping[i];
    // RDKit✔️❌:     if (comMap.find(mi) == comMap.end()) {
    // RDKit✔️❌:       INT_VECT comp;
    // RDKit✔️❌:       comMap[mi] = comp;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     comMap[mi].push_back(i);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   for (INT_INT_VECT_MAP_CI mci = comMap.begin(); mci != comMap.end(); mci++) {
    // RDKit✔️❌:     frags.push_back((*mci).second);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return rdcast<unsigned int>(frags.size());
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION MolOps::getMolFrags component mapping
    // Complexity review: local DFS is O(V + E), but sorting each component adds
    // O(sum(|component| log |component|)); RDKit emits each component's atom rows
    // in source order while grouping the connected-component labels.
    let mut seen = vec![false; topology.atoms.len()];
    let mut components = Vec::new();
    for start in 0..topology.atoms.len() {
        if seen[start] {
            continue;
        }
        seen[start] = true;
        let mut component = Vec::new();
        let mut pending = vec![start];
        while let Some(atom) = pending.pop() {
            component.push(atom);
            for neighbor in topology.adjacency.neighbors_of(atom) {
                if !seen[neighbor.atom_index] {
                    seen[neighbor.atom_index] = true;
                    pending.push(neighbor.atom_index);
                }
            }
        }
        component.sort_unstable();
        components.push(component);
    }
    components
}

fn reject_unmodeled_stereochemical_writing(
    topology: &TopologyBlock,
    doing_cx_smiles: bool,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION GetBondSmiles base token dispatch
    // RDKit❗❌:   switch (bond->getBondType()) {
    // RDKit❗❌:     case Bond::SINGLE:
    // RDKit❗❌:       if (dir != Bond::NONE && dir != Bond::UNKNOWN) {
    // END RDKIT CPP FUNCTION GetBondSmiles base token dispatch
    // Behavior review: the base token branch consumes bond kind and direction;
    // CX Atrop stereo stays on the prepared record for the later
    // CX_BOND_ATROPISOMER extension. Keep ordinary SMILES' existing gate.
    // Complexity review: this existing preflight scans all bonds before the
    // output traversal, an O(E) extra pass; these match alternatives add no
    // new scan, allocation or per-bond data structure.
    if topology.atoms.iter().any(|atom| {
        !matches!(
            atom.chiral_tag(),
            ChiralTag::Unspecified
                | ChiralTag::TetrahedralCw
                | ChiralTag::TetrahedralCcw
                | ChiralTag::SquarePlanar
                | ChiralTag::TrigonalBipyramidal
                | ChiralTag::Octahedral
        ) || atom.unknown_stereo()
    }) {
        return Err(SmilesParseError::UnsupportedWriter(
            "allene, generic, or unknown atom stereochemistry is not modeled by the detached writer",
        ));
    }
    if topology.bonds.iter().any(|bond| {
        (!doing_cx_smiles
            && (!matches!(
                bond.direction(),
                BondDirection::None | BondDirection::EndDownRight | BondDirection::EndUpRight
            ) || !matches!(
                bond.stereo(),
                // RDKit✔️✔️:         if (bond->getStereo() == Bond::BondStereo::STEREOANY) {
                // RDKit✔️✔️:           bond->setStereo(Bond::BondStereo::STEREONONE);
                // RDKit✔️✔️:         }
                // The existing writer-owned copy executes this source cleanup
                // before ranking. ANY is supported input; preserve its source
                // value while serializing the cleaned copy, O(1) per bond.
                BondStereo::None
                    | BondStereo::Any
                    | BondStereo::E
                    | BondStereo::Z
                    | BondStereo::Cis
                    | BondStereo::Trans
            ) || bond.unknown_stereo()))
            || (doing_cx_smiles
                && (!matches!(
                    bond.direction(),
                    BondDirection::None
                        | BondDirection::BeginWedge
                        | BondDirection::BeginDash
                        | BondDirection::EndDownRight
                        | BondDirection::EndUpRight
                        | BondDirection::EitherDouble
                        | BondDirection::Unknown
                ) || !matches!(
                    bond.stereo(),
                    BondStereo::None
                        | BondStereo::Any
                        | BondStereo::E
                        | BondStereo::Z
                        | BondStereo::Cis
                        | BondStereo::Trans
                        | BondStereo::AtropCw
                        | BondStereo::AtropCcw
                )))
    }) {
        return Err(SmilesParseError::UnsupportedWriter(
            "unknown, atropisomeric, or non-directional bond stereochemistry is not modeled by the detached writer",
        ));
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn dfs_find_cycles(
    topology: &TopologyBlock,
    atom: usize,
    incoming_bond: Option<BondId>,
    colors: &mut [AtomColor],
    ranks: &[i64],
    ring_bonds: &[bool],
    bonds_in_play: Option<&[bool]>,
    bond_symbols: Option<&[String]>,
    atom_ring_closures: &mut [Vec<BondId>],
    random_stream: Option<&mut cosmolkit_core::RdkitRandomGenerator<'_>>,
) -> Result<(), SmilesParseError> {
    let ranks = ranks.iter().map(|r| *r as u32).collect::<Vec<_>>();
    dfs_find_cycles_source(
        topology,
        atom,
        incoming_bond,
        colors,
        &ranks,
        ring_bonds,
        bonds_in_play,
        bond_symbols.map(DfsBondSymbols::Utf8),
        atom_ring_closures,
        random_stream,
    )
}

#[allow(clippy::too_many_arguments)]
fn dfs_find_cycles_source<G: SourceTraversalGraph>(
    graph: &G,
    atom: usize,
    incoming_bond: Option<BondId>,
    colors: &mut [AtomColor],
    ranks: &[u32],
    ring_bonds: &[bool],
    bonds_in_play: Option<&[bool]>,
    bond_symbols: Option<DfsBondSymbols<'_>>,
    atom_ring_closures: &mut [Vec<BondId>],
    mut random_stream: Option<&mut cosmolkit_core::RdkitRandomGenerator<'_>>,
) -> Result<(), SmilesParseError> {
    // RDKit❗❌: void dfsFindCycles(ROMol &mol, int atomIdx, int inBondIdx,
    // RDKit❗❌:                    std::vector<AtomColors> &colors, const UINT_VECT &ranks,
    // RDKit❗❌:                    VECT_INT_VECT &atomRingClosures,
    // RDKit❗❌:                    const boost::dynamic_bitset<> *bondsInPlay,
    // RDKit❗❌:                    const std::vector<std::string> *bondSymbols, bool doRandom) {
    // RDKit❗❌:   Atom *atom = mol.getAtomWithIdx(atomIdx);
    // RDKit❗❌:
    // RDKit❗❌:   colors[atomIdx] = GREY_NODE;
    // RDKit❗❌:
    // RDKit❗❌:   // ---------------------
    // RDKit❗❌:   //
    // RDKit❗❌:   //  Build the list of possible destinations from here
    // RDKit❗❌:   //
    // RDKit❗❌:   // ---------------------
    // RDKit❗❌:   std::vector<PossibleType> possibles;
    // RDKit❗❌:   auto bondsPair = mol.getAtomBonds(atom);
    // RDKit❗❌:   possibles.reserve(bondsPair.second - bondsPair.first);
    // RDKit❗❌:
    // RDKit❗❌:   while (bondsPair.first != bondsPair.second) {
    // RDKit❗❌:     Bond *theBond = mol[*(bondsPair.first)];
    // RDKit❗❌:     ++bondsPair.first;
    // RDKit❗❌:     if (bondsInPlay && !(*bondsInPlay)[theBond->getIdx()]) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:     if (inBondIdx < 0 ||
    // RDKit❗❌:         theBond->getIdx() != static_cast<unsigned int>(inBondIdx)) {
    // RDKit❗❌:       int otherIdx = theBond->getOtherAtomIdx(atomIdx);
    // RDKit❗❌:       auto rank = ranks[otherIdx];
    // RDKit❗❌:       // ---------------------
    // RDKit❗❌:       //
    // RDKit❗❌:       // things are a bit more complicated if we are sitting on a
    // RDKit❗❌:       // ring atom. we would like to traverse first to the
    // RDKit❗❌:       // ring-closure atoms, then to atoms outside the ring first,
    // RDKit❗❌:       // then to atoms in the ring that haven't already been visited
    // RDKit❗❌:       // (non-ring-closure atoms).
    // RDKit❗❌:       //
    // RDKit❗❌:       //  Here's how the black magic works:
    // RDKit❗❌:       //   - non-ring atom neighbors have their original ranks
    // RDKit❗❌:       //   - ring atom neighbors have this added to their ranks:
    // RDKit❗❌:       //       (MAX_BONDTYPE - bondOrder)*MAX_NATOMS*MAX_NATOMS
    // RDKit❗❌:       //   - ring-closure neighbors lose a factor of:
    // RDKit❗❌:       //       (MAX_BONDTYPE+1)*MAX_NATOMS*MAX_NATOMS
    // RDKit❗❌:       //
    // RDKit❗❌:       //  This tactic biases us to traverse to non-ring neighbors first,
    // RDKit❗❌:       //  original ordering if bond orders are all equal... crafty, neh?
    // RDKit❗❌:       //
    // RDKit❗❌:       // ---------------------
    // RDKit❗❌:       if (!doRandom) {
    // RDKit❗❌:         if (colors[otherIdx] == GREY_NODE) {
    // RDKit❗❌:           rank -= static_cast<int>(MAX_BONDTYPE + 1) * MAX_NATOMS * MAX_NATOMS;
    // RDKit❗❌:           if (!bondSymbols) {
    // RDKit❗❌:             rank += static_cast<int>(MAX_BONDTYPE - theBond->getBondType()) *
    // RDKit❗❌:                     MAX_NATOMS;
    // RDKit❗❌:           } else {
    // RDKit❗❌:             const std::string &symb = (*bondSymbols)[theBond->getIdx()];
    // RDKit❗❌:             std::uint32_t hsh = gboost::hash_range(symb.begin(), symb.end());
    // RDKit❗❌:             rank += (hsh % MAX_NATOMS) * MAX_NATOMS;
    // RDKit❗❌:           }
    // RDKit❗❌:         } else if (theBond->getOwningMol().getRingInfo()->numBondRings(
    // RDKit❗❌:                        theBond->getIdx())) {
    // RDKit❗❌:           if (!bondSymbols) {
    // RDKit❗❌:             rank += static_cast<int>(MAX_BONDTYPE - theBond->getBondType()) *
    // RDKit❗❌:                     MAX_NATOMS * MAX_NATOMS;
    // RDKit❗❌:           } else {
    // RDKit❗❌:             const std::string &symb = (*bondSymbols)[theBond->getIdx()];
    // RDKit❗❌:             std::uint32_t hsh = gboost::hash_range(symb.begin(), symb.end());
    // RDKit❗❌:             rank += (hsh % MAX_NATOMS) * MAX_NATOMS * MAX_NATOMS;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       } else {
    // RDKit❗❌:         // randomize the rank
    // RDKit❗❌:         rank = getRandomGenerator()();
    // RDKit❗❌:       }
    // RDKit❗❌:       // std::cerr << "            " << atomIdx << ": " << otherIdx << " " <<
    // RDKit❗❌:       // rank
    // RDKit❗❌:       //           << std::endl;
    // RDKit❗❌:       // std::cerr<<"aIdx: "<< atomIdx <<"   p: "<<otherIdx<<" Rank:
    // RDKit❗❌:       // "<<ranks[otherIdx] <<" "<<colors[otherIdx]<<"
    // RDKit❗❌:       // "<<theBond->getBondType()<<" "<<rank<<std::endl;
    // RDKit❗❌:       possibles.emplace_back(rank, otherIdx, theBond);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // ---------------------
    // RDKit❗❌:   //
    // RDKit❗❌:   //  Sort on ranks
    // RDKit❗❌:   //
    // RDKit❗❌:   // ---------------------
    // RDKit❗❌:   std::sort(possibles.begin(), possibles.end(), _possibleCompare);
    // RDKit❗❌:   // if (possibles.size())
    // RDKit❗❌:   //   std::cerr << " aIdx1: " << atomIdx
    // RDKit❗❌:   //             << " first: " << possibles.front()std:std::get<0>() << " "
    // RDKit❗❌:   //             << possibles.front()std:std::get<1>() << std::endl;
    // RDKit❗❌:   // // ---------------------
    // RDKit❗❌:   //
    // RDKit❗❌:   //  Now work the children
    // RDKit❗❌:   //
    // RDKit❗❌:   // ---------------------
    // RDKit❗❌:   for (auto &possible : possibles) {
    // RDKit❗❌:     int possibleIdx = std::get<1>(possible);
    // RDKit❗❌:     Bond *bond = std::get<2>(possible);
    // RDKit❗❌:     switch (colors[possibleIdx]) {
    // RDKit❗❌:       case WHITE_NODE:
    // RDKit❗❌:         // -----
    // RDKit❗❌:         // we haven't seen this node at all before, traverse
    // RDKit❗❌:         // -----
    // RDKit❗❌:         dfsFindCycles(mol, possibleIdx, bond->getIdx(), colors, ranks,
    // RDKit❗❌:                       atomRingClosures, bondsInPlay, bondSymbols, doRandom);
    // RDKit❗❌:         break;
    // RDKit❗❌:       case GREY_NODE:
    // RDKit❗❌:         // -----
    // RDKit❗❌:         // we've seen this, but haven't finished it (we're finishing a ring)
    // RDKit❗❌:         // -----
    // RDKit❗❌:         atomRingClosures[possibleIdx].push_back(bond->getIdx());
    // RDKit❗❌:         atomRingClosures[atomIdx].push_back(bond->getIdx());
    // RDKit❗❌:         break;
    // RDKit❗❌:       default:
    // RDKit❗❌:         // -----
    // RDKit❗❌:         // this node has been finished. don't do anything.
    // RDKit❗❌:         // -----
    // RDKit❗❌:         break;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   colors[atomIdx] = BLACK_NODE;
    // RDKit❗❌: }
    // Native uint32 arithmetic precedes the signed PossibleType conversion.
    // Preserve mask-before-incoming order and draws for every non-incoming
    // neighbor including already Black nodes. Equal-rank sort permutations
    // remain a final comparison gap; no index tie rule is added.
    // The full-degree candidate vector/sort/recursive passes match source;
    // existing Topology callers pay one rank-width projection allocation.
    if atom >= graph.source_atom_count() {
        return Err(SmilesParseError::TraversalStateIndex {
            state: "atom",
            index: atom,
            count: graph.source_atom_count(),
        });
    }
    let count = colors.len();
    *colors
        .get_mut(atom)
        .ok_or(SmilesParseError::TraversalStateIndex {
            state: "colors",
            index: atom,
            count,
        })? = AtomColor::Grey;
    let mut possibles = Vec::with_capacity(graph.source_neighbors(atom).len());
    for (other, bond) in graph.source_neighbors(atom) {
        if let Some(mask) = bonds_in_play {
            if !*mask
                .get(bond.index())
                .ok_or(SmilesParseError::TraversalStateIndex {
                    state: "bondsInPlay",
                    index: bond.index(),
                    count: mask.len(),
                })?
            {
                continue;
            }
        }
        if Some(bond) == incoming_bond {
            continue;
        }
        let mut rank = *ranks
            .get(other)
            .ok_or(SmilesParseError::TraversalStateIndex {
                state: "ranks",
                index: other,
                count: ranks.len(),
            })?;
        if let Some(random) = random_stream.as_deref_mut() {
            rank = random.next_u32();
        } else {
            let color = *colors
                .get(other)
                .ok_or(SmilesParseError::TraversalStateIndex {
                    state: "colors",
                    index: other,
                    count: colors.len(),
                })?;
            if color == AtomColor::Grey {
                rank = rank.wrapping_sub(((MAX_BONDTYPE + 1) * MAX_NATOMS * MAX_NATOMS) as u32);
                let salt = if let Some(symbols) = bond_symbols {
                    (boost_hash_range(symbols.get(bond)?) % MAX_NATOMS as u32)
                        .wrapping_mul(MAX_NATOMS as u32)
                } else {
                    ((MAX_BONDTYPE - rdkit_bond_type_code(graph.source_bond(bond)?.order()))
                        * MAX_NATOMS) as u32
                };
                rank = rank.wrapping_add(salt);
            } else if *ring_bonds.get(bond.index()).ok_or(
                SmilesParseError::TraversalStateIndex {
                    state: "ringInfo",
                    index: bond.index(),
                    count: ring_bonds.len(),
                },
            )? {
                let salt = if let Some(symbols) = bond_symbols {
                    (boost_hash_range(symbols.get(bond)?) % MAX_NATOMS as u32)
                        .wrapping_mul(MAX_NATOMS as u32)
                        .wrapping_mul(MAX_NATOMS as u32)
                } else {
                    ((MAX_BONDTYPE - rdkit_bond_type_code(graph.source_bond(bond)?.order()))
                        * MAX_NATOMS
                        * MAX_NATOMS) as u32
                };
                rank = rank.wrapping_add(salt);
            }
        }
        possibles.push(Possible {
            rank: i64::from(rank as i32),
            atom: other,
            bond,
        });
    }
    possibles.sort_unstable_by_key(|p| p.rank);
    for possible in possibles {
        match colors[possible.atom] {
            AtomColor::White => dfs_find_cycles_source(
                graph,
                possible.atom,
                Some(possible.bond),
                colors,
                ranks,
                ring_bonds,
                bonds_in_play,
                bond_symbols,
                atom_ring_closures,
                random_stream.as_deref_mut(),
            )?,
            AtomColor::Grey => {
                let count = atom_ring_closures.len();
                atom_ring_closures
                    .get_mut(possible.atom)
                    .ok_or(SmilesParseError::TraversalStateIndex {
                        state: "atomRingClosures",
                        index: possible.atom,
                        count,
                    })?
                    .push(possible.bond);
                atom_ring_closures
                    .get_mut(atom)
                    .ok_or(SmilesParseError::TraversalStateIndex {
                        state: "atomRingClosures",
                        index: atom,
                        count,
                    })?
                    .push(possible.bond);
            }
            AtomColor::Black => {}
        }
    }
    colors[atom] = AtomColor::Black;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn canonical_dfs_traversal_source<G: SourceTraversalGraph>(
    graph: &mut G,
    atom: usize,
    incoming_bond: Option<BondId>,
    colors: &mut [AtomColor],
    ranks: &[u32],
    ring_bonds: &[bool],
    stack: &mut Vec<MolStackElem>,
    atom_ring_closures: &mut [Vec<BondId>],
    atom_traversal_bond_order: &mut [Vec<BondId>],
    bonds_in_play: Option<&[bool]>,
    bond_symbols: Option<DfsBondSymbols<'_>>,
    mut random_stream: Option<&mut cosmolkit_core::RdkitRandomGenerator<'_>>,
) -> Result<Vec<bool>, SmilesParseError> {
    // RDKit❗❌: void canonicalDFSTraversal(ROMol &mol, int atomIdx, int inBondIdx,
    // RDKit❗❌:                            std::vector<AtomColors> &colors,
    // RDKit❗❌:                            const UINT_VECT &ranks, MolStack &molStack,
    // RDKit❗❌:                            VECT_INT_VECT &atomRingClosures,
    // RDKit❗❌:                            std::vector<INT_LIST> &atomTraversalBondOrder,
    // RDKit❗❌:                            const boost::dynamic_bitset<> *bondsInPlay,
    // RDKit❗❌:                            const std::vector<std::string> *bondSymbols,
    // RDKit❗❌:                            bool doRandom) {
    // RDKit❗❌:   PRECONDITION(colors.size() >= mol.getNumAtoms(), "vector too small");
    // RDKit❗❌:   PRECONDITION(ranks.size() >= mol.getNumAtoms(), "vector too small");
    // RDKit❗❌:   PRECONDITION(atomRingClosures.size() >= mol.getNumAtoms(),
    // RDKit❗❌:                "vector too small");
    // RDKit❗❌:   PRECONDITION(atomTraversalBondOrder.size() >= mol.getNumAtoms(),
    // RDKit❗❌:                "vector too small");
    // RDKit❗❌:   PRECONDITION(!bondsInPlay || bondsInPlay->size() >= mol.getNumBonds(),
    // RDKit❗❌:                "bondsInPlay too small");
    // RDKit❗❌:   PRECONDITION(!bondSymbols || bondSymbols->size() >= mol.getNumBonds(),
    // RDKit❗❌:                "bondSymbols too small");
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<AtomColors> tcolors(colors.begin(), colors.end());
    // RDKit❗❌:   dfsFindCycles(mol, atomIdx, inBondIdx, tcolors, ranks, atomRingClosures,
    // RDKit❗❌:                 bondsInPlay, bondSymbols, doRandom);
    // RDKit❗❌:
    // RDKit❗❌:   boost::dynamic_bitset<> cyclesAvailable(MAX_CYCLES);
    // RDKit❗❌:   cyclesAvailable.set();
    // RDKit❗❌:   dfsBuildStack(mol, atomIdx, inBondIdx, colors, ranks, cyclesAvailable,
    // RDKit❗❌:                 molStack, atomRingClosures, atomTraversalBondOrder, bondsInPlay,
    // RDKit❗❌:                 bondSymbols, doRandom);
    // RDKit❗❌: }
    // Same full color copy before cycle discovery and fresh1024-slot bitmap;
    // no whole graph clone. Byte bitmap/rank-view projection costs and native
    // pointer/allocator/equal-key-sort gaps are explicit shared DFS limitations.
    let n = graph.source_atom_count();
    let m = graph.source_bond_count();
    for (name, len, minimum) in [
        ("colors", colors.len(), n),
        ("ranks", ranks.len(), n),
        ("atomRingClosures", atom_ring_closures.len(), n),
        ("atomTraversalBondOrder", atom_traversal_bond_order.len(), n),
    ] {
        if len < minimum {
            return Err(SmilesParseError::TraversalStateIndex {
                state: name,
                index: minimum,
                count: len,
            });
        }
    }
    if let Some(mask) = bonds_in_play {
        if mask.len() < m {
            return Err(SmilesParseError::TraversalStateIndex {
                state: "bondsInPlay",
                index: m,
                count: mask.len(),
            });
        }
    }
    if let Some(symbols) = bond_symbols {
        if symbols.len() < m {
            return Err(SmilesParseError::TraversalStateIndex {
                state: "bondSymbols",
                index: m,
                count: symbols.len(),
            });
        }
    }
    let mut cycle_colors = colors.to_vec();
    dfs_find_cycles_source(
        graph,
        atom,
        incoming_bond,
        &mut cycle_colors,
        ranks,
        ring_bonds,
        bonds_in_play,
        bond_symbols,
        atom_ring_closures,
        random_stream.as_deref_mut(),
    )?;
    let mut available = vec![true; MAX_CYCLES];
    let mut ids = vec![None; m];
    let mut opened = vec![false; m];
    dfs_build_stack(
        graph,
        atom,
        incoming_bond,
        colors,
        ranks,
        ring_bonds,
        atom_ring_closures,
        &mut ids,
        &mut available,
        stack,
        atom_traversal_bond_order,
        &mut opened,
        bonds_in_play,
        bond_symbols,
        random_stream,
    )?;
    Ok(opened)
}

// One notation-owned DFS body supports both existing detached graph carriers.
// The closed private trait exposes only fields reached by the pinned function;
// it grants no live Molecule/runtime authority and clones no graph or queries.
trait SourceTraversalGraph {
    fn source_atom_count(&self) -> usize;
    fn source_bond_count(&self) -> usize;
    fn source_neighbors(&self, atom: usize) -> impl ExactSizeIterator<Item = (usize, BondId)>;
    fn source_bond(&self, bond: BondId) -> Result<&Bond, SmilesParseError>;
    fn source_bond_mut(&mut self, bond: BondId) -> Result<&mut Bond, SmilesParseError>;
    fn clear_source_traversal_order(&mut self, atom: usize) -> Result<(), SmilesParseError>;
}

impl SourceTraversalGraph for TopologyBlock {
    fn source_atom_count(&self) -> usize {
        self.atoms.len()
    }
    fn source_bond_count(&self) -> usize {
        self.bonds.len()
    }
    fn source_neighbors(&self, atom: usize) -> impl ExactSizeIterator<Item = (usize, BondId)> {
        self.adjacency
            .neighbors_of(atom)
            .iter()
            .map(|n| (n.atom_index, n.bond))
    }
    fn source_bond(&self, bond: BondId) -> Result<&Bond, SmilesParseError> {
        self.bonds
            .get(bond.index())
            .ok_or(SmilesParseError::TraversalStateIndex {
                state: "bond",
                index: bond.index(),
                count: self.bonds.len(),
            })
    }
    fn source_bond_mut(&mut self, bond: BondId) -> Result<&mut Bond, SmilesParseError> {
        let count = self.bonds.len();
        self.bonds
            .get_mut(bond.index())
            .ok_or(SmilesParseError::TraversalStateIndex {
                state: "bond",
                index: bond.index(),
                count,
            })
    }
    fn clear_source_traversal_order(&mut self, atom: usize) -> Result<(), SmilesParseError> {
        let count = self.atoms.len();
        self.atoms
            .get_mut(atom)
            .ok_or(SmilesParseError::TraversalStateIndex {
                state: "atom",
                index: atom,
                count,
            })?
            .clear_prop("_TraversalBondIndexOrder")?;
        Ok(())
    }
}

impl SourceTraversalGraph for cosmolkit_model::QueryGraph {
    fn source_atom_count(&self) -> usize {
        self.num_atoms()
    }
    fn source_bond_count(&self) -> usize {
        self.num_bonds()
    }
    fn source_neighbors(&self, atom: usize) -> impl ExactSizeIterator<Item = (usize, BondId)> {
        self.adjacency()[atom]
            .iter()
            .map(|(a, b)| (*a, BondId::new(*b)))
    }
    fn source_bond(&self, bond: BondId) -> Result<&Bond, SmilesParseError> {
        self.bond(bond.index())
            .map(|b| b.bond())
            .ok_or(SmilesParseError::TraversalStateIndex {
                state: "bond",
                index: bond.index(),
                count: self.num_bonds(),
            })
    }
    fn source_bond_mut(&mut self, bond: BondId) -> Result<&mut Bond, SmilesParseError> {
        let count = self.num_bonds();
        self.bonds_mut()
            .get_mut(bond.index())
            .map(|b| b.bond_mut())
            .ok_or(SmilesParseError::TraversalStateIndex {
                state: "bond",
                index: bond.index(),
                count,
            })
    }
    fn clear_source_traversal_order(&mut self, atom: usize) -> Result<(), SmilesParseError> {
        let count = self.num_atoms();
        self.atom_mut(atom)
            .ok_or(SmilesParseError::TraversalStateIndex {
                state: "atom",
                index: atom,
                count,
            })?
            .clear_prop("_TraversalBondIndexOrder")?;
        Ok(())
    }
}

/// Borrowed native byte symbols; the existing UTF-8 caller is one valid subset.
#[doc(hidden)]
#[derive(Clone, Copy)]
pub enum DfsBondSymbols<'a> {
    Utf8(&'a [String]),
    Bytes(&'a [PropertyText]),
}

impl<'a> DfsBondSymbols<'a> {
    fn len(self) -> usize {
        match self {
            Self::Utf8(v) => v.len(),
            Self::Bytes(v) => v.len(),
        }
    }
    fn get(self, bond: BondId) -> Result<&'a [u8], SmilesParseError> {
        let (value, count) = match self {
            Self::Utf8(s) => (s.get(bond.index()).map(|v| v.as_bytes()), s.len()),
            Self::Bytes(s) => (s.get(bond.index()).map(|v| v.as_bytes()), s.len()),
        };
        value.ok_or(SmilesParseError::TraversalStateIndex {
            state: "bondSymbols",
            index: bond.index(),
            count,
        })
    }
}

/// Query projection of the same source DFS; all supplied scratch/output state
/// and graph properties retain native prefix mutation on a later failure.
#[doc(hidden)]
#[allow(clippy::too_many_arguments)]
pub fn dfs_build_query_stack(
    graph: &mut cosmolkit_model::QueryGraph,
    atom: usize,
    incoming_bond: Option<BondId>,
    colors: &mut [AtomColor],
    ranks: &[u32],
    ring_bonds: &[bool],
    atom_ring_closures: &[Vec<BondId>],
    ring_ids: &mut [Option<usize>],
    available_ring_ids: &mut [bool],
    stack: &mut Vec<MolStackElem>,
    atom_traversal_bond_order: &mut [Vec<BondId>],
    traversal_ring_closure_bonds: &mut [bool],
    bonds_in_play: Option<&[bool]>,
    bond_symbols: Option<DfsBondSymbols<'_>>,
    random_stream: Option<&mut cosmolkit_core::RdkitRandomGenerator<'_>>,
) -> Result<(), SmilesParseError> {
    dfs_build_stack(
        graph,
        atom,
        incoming_bond,
        colors,
        ranks,
        ring_bonds,
        atom_ring_closures,
        ring_ids,
        available_ring_ids,
        stack,
        atom_traversal_bond_order,
        traversal_ring_closure_bonds,
        bonds_in_play,
        bond_symbols,
        random_stream,
    )
}

#[allow(clippy::too_many_arguments)]
fn dfs_build_stack<G: SourceTraversalGraph>(
    graph: &mut G,
    atom: usize,
    incoming_bond: Option<BondId>,
    colors: &mut [AtomColor],
    ranks: &[u32],
    ring_bonds: &[bool],
    atom_ring_closures: &[Vec<BondId>],
    ring_ids: &mut [Option<usize>],
    available_ring_ids: &mut [bool],
    stack: &mut Vec<MolStackElem>,
    atom_traversal_bond_order: &mut [Vec<BondId>],
    traversal_ring_closure_bonds: &mut [bool],
    bonds_in_play: Option<&[bool]>,
    bond_symbols: Option<DfsBondSymbols<'_>>,
    mut random_stream: Option<&mut cosmolkit_core::RdkitRandomGenerator<'_>>,
) -> Result<(), SmilesParseError> {
    // RDKit❗❌: void dfsBuildStack(ROMol &mol, int atomIdx, int inBondIdx,
    // RDKit❗❌:                    std::vector<AtomColors> &colors, const UINT_VECT &ranks,
    // RDKit❗❌:                    boost::dynamic_bitset<> &cyclesAvailable, MolStack &molStack,
    // RDKit❗❌:                    VECT_INT_VECT &atomRingClosures,
    // RDKit❗❌:                    std::vector<INT_LIST> &atomTraversalBondOrder,
    // RDKit❗❌:                    const boost::dynamic_bitset<> *bondsInPlay,
    // RDKit❗❌:                    const std::vector<std::string> *bondSymbols, bool doRandom) {
    // RDKit❗❌:   Atom *atom = mol.getAtomWithIdx(atomIdx);
    // RDKit❗❌:   boost::dynamic_bitset<> seenFromHere(mol.getNumAtoms());
    // RDKit❗❌:
    // RDKit❗❌:   seenFromHere.set(atomIdx);
    // RDKit❗❌:   molStack.push_back(MolStackElem(atom));
    // RDKit❗❌:   colors[atomIdx] = GREY_NODE;
    // RDKit❗❌:
    // RDKit❗❌:   INT_LIST travList;
    // RDKit❗❌:   if (inBondIdx >= 0) {
    // RDKit❗❌:     travList.push_back(inBondIdx);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // ---------------------
    // RDKit❗❌:   //
    // RDKit❗❌:   //  Add any ring closures
    // RDKit❗❌:   //
    // RDKit❗❌:   // ---------------------
    // RDKit❗❌:   if (!atomRingClosures[atomIdx].empty()) {
    // RDKit❗❌:     std::vector<unsigned int> ringsClosed;
    // RDKit❗❌:     for (auto bIdx : atomRingClosures[atomIdx]) {
    // RDKit❗❌:       travList.push_back(bIdx);
    // RDKit❗❌:       Bond *bond = mol.getBondWithIdx(bIdx);
    // RDKit❗❌:       seenFromHere.set(bond->getOtherAtomIdx(atomIdx));
    // RDKit❗❌:       unsigned int ringIdx = std::numeric_limits<unsigned int>::max();
    // RDKit❗❌:       if (bond->getPropIfPresent(common_properties::_TraversalRingClosureBond,
    // RDKit❗❌:                                  ringIdx)) {
    // RDKit❗❌:         // this is end of the ring closure
    // RDKit❗❌:         // we can just pull the ring index from the bond itself:
    // RDKit❗❌:         molStack.push_back(MolStackElem(bond, atomIdx));
    // RDKit❗❌:         molStack.push_back(MolStackElem(ringIdx));
    // RDKit❗❌:         // don't make the ring digit immediately available again: we don't want
    // RDKit❗❌:         // to have the same
    // RDKit❗❌:         // ring digit opening and closing rings on an atom.
    // RDKit❗❌:         ringsClosed.push_back(ringIdx - 1);
    // RDKit❗❌:       } else {
    // RDKit❗❌:         // this is the beginning of the ring closure, we need to come up with a
    // RDKit❗❌:         // ring index:
    // RDKit❗❌:         auto lowestRingIdx = cyclesAvailable.find_first();
    // RDKit❗❌:         if (lowestRingIdx == boost::dynamic_bitset<>::npos) {
    // RDKit❗❌:           throw ValueErrorException(
    // RDKit❗❌:               "Too many rings open at once. SMILES cannot be generated.");
    // RDKit❗❌:         }
    // RDKit❗❌:         cyclesAvailable.set(lowestRingIdx, false);
    // RDKit❗❌:         ++lowestRingIdx;
    // RDKit❗❌:         bond->setProp(common_properties::_TraversalRingClosureBond,
    // RDKit❗❌:                       static_cast<unsigned int>(lowestRingIdx));
    // RDKit❗❌:         molStack.push_back(MolStackElem(lowestRingIdx));
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     for (auto ringIdx : ringsClosed) {
    // RDKit❗❌:       cyclesAvailable.set(ringIdx);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // ---------------------
    // RDKit❗❌:   //
    // RDKit❗❌:   //  Build the list of possible destinations from here
    // RDKit❗❌:   //
    // RDKit❗❌:   // ---------------------
    // RDKit❗❌:   std::vector<PossibleType> possibles;
    // RDKit❗❌:   possibles.reserve(atom->getDegree());
    // RDKit❗❌:   for (auto theBond : mol.atomBonds(atom)) {
    // RDKit❗❌:     if (bondsInPlay && !(*bondsInPlay)[theBond->getIdx()]) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:     if (inBondIdx < 0 ||
    // RDKit❗❌:         theBond->getIdx() != static_cast<unsigned int>(inBondIdx)) {
    // RDKit❗❌:       int otherIdx = theBond->getOtherAtomIdx(atomIdx);
    // RDKit❗❌:       // ---------------------
    // RDKit❗❌:       //
    // RDKit❗❌:       // This time we skip the ring-closure atoms (we did them
    // RDKit❗❌:       // above); we want to traverse first to atoms outside the ring
    // RDKit❗❌:       // then to atoms in the ring that haven't already been visited
    // RDKit❗❌:       // (non-ring-closure atoms).
    // RDKit❗❌:       //
    // RDKit❗❌:       // otherwise it's the same ranking logic as above
    // RDKit❗❌:       //
    // RDKit❗❌:       // ---------------------
    // RDKit❗❌:       if (colors[otherIdx] != WHITE_NODE || seenFromHere[otherIdx]) {
    // RDKit❗❌:         // ring closure or finished atom... skip it.
    // RDKit❗❌:         continue;
    // RDKit❗❌:       }
    // RDKit❗❌:       auto rank = ranks[otherIdx];
    // RDKit❗❌:       if (!doRandom) {
    // RDKit❗❌:         if (theBond->getOwningMol().getRingInfo()->numBondRings(
    // RDKit❗❌:                 theBond->getIdx())) {
    // RDKit❗❌:           if (!bondSymbols) {
    // RDKit❗❌:             rank += static_cast<int>(MAX_BONDTYPE - theBond->getBondType()) *
    // RDKit❗❌:                     MAX_NATOMS * MAX_NATOMS;
    // RDKit❗❌:           } else {
    // RDKit❗❌:             const std::string &symb = (*bondSymbols)[theBond->getIdx()];
    // RDKit❗❌:             std::uint32_t hsh = gboost::hash_range(symb.begin(), symb.end());
    // RDKit❗❌:             rank += (hsh % MAX_NATOMS) * MAX_NATOMS * MAX_NATOMS;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       } else {
    // RDKit❗❌:         // randomize the rank
    // RDKit❗❌:         rank = getRandomGenerator()();
    // RDKit❗❌:       }
    // RDKit❗❌:
    // RDKit❗❌:       possibles.emplace_back(rank, otherIdx, theBond);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // ---------------------
    // RDKit❗❌:   //
    // RDKit❗❌:   //  Sort on ranks
    // RDKit❗❌:   //
    // RDKit❗❌:   // ---------------------
    // RDKit❗❌:   std::sort(possibles.begin(), possibles.end(), _possibleCompare);
    // RDKit❗❌:   // if (possibles.size())
    // RDKit❗❌:   //   std::cerr << " aIdx2: " << atomIdx
    // RDKit❗❌:   //             << " first: " << possibles.front()std:std::get<0>() << " "
    // RDKit❗❌:   //             << possibles.front()std:std::get<1>() << std::endl;
    // RDKit❗❌:
    // RDKit❗❌:   // ---------------------
    // RDKit❗❌:   //
    // RDKit❗❌:   //  Now work the children
    // RDKit❗❌:   //
    // RDKit❗❌:   // ---------------------
    // RDKit❗❌:   for (auto possiblesIt = possibles.begin(); possiblesIt != possibles.end();
    // RDKit❗❌:        ++possiblesIt) {
    // RDKit❗❌:     int possibleIdx = std::get<1>(*possiblesIt);
    // RDKit❗❌:     if (colors[possibleIdx] != WHITE_NODE) {
    // RDKit❗❌:       // we're either done or it's a ring-closure, which we already processed...
    // RDKit❗❌:       // this test isn't strictly required, because we only added WHITE notes to
    // RDKit❗❌:       // the possibles list, but it seems logical to document it
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:     Bond *bond = std::get<2>(*possiblesIt);
    // RDKit❗❌:     Atom *otherAtom = mol.getAtomWithIdx(possibleIdx);
    // RDKit❗❌:     // ww might have some residual data from earlier calls, clean that up:
    // RDKit❗❌:     otherAtom->clearProp(common_properties::_TraversalBondIndexOrder);
    // RDKit❗❌:     travList.push_back(bond->getIdx());
    // RDKit❗❌:     if (possiblesIt + 1 != possibles.end()) {
    // RDKit❗❌:       // we're branching
    // RDKit❗❌:       molStack.push_back(
    // RDKit❗❌:           MolStackElem("(", rdcast<int>(possiblesIt - possibles.begin())));
    // RDKit❗❌:     }
    // RDKit❗❌:     molStack.push_back(MolStackElem(bond, atomIdx));
    // RDKit❗❌:     dfsBuildStack(mol, possibleIdx, bond->getIdx(), colors, ranks,
    // RDKit❗❌:                   cyclesAvailable, molStack, atomRingClosures,
    // RDKit❗❌:                   atomTraversalBondOrder, bondsInPlay, bondSymbols, doRandom);
    // RDKit❗❌:     if (possiblesIt + 1 != possibles.end()) {
    // RDKit❗❌:       molStack.push_back(
    // RDKit❗❌:           MolStackElem(")", rdcast<int>(possiblesIt - possibles.begin())));
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   atomTraversalBondOrder[atom->getIdx()] = travList;
    // RDKit❗❌:   colors[atomIdx] = BLACK_NODE;
    // RDKit❗❌: }
    // RDKit❗❌: typedef std::tuple<int, int, Bond *> PossibleType;
    // RDKit❗❌: auto _possibleCompare = [](const PossibleType &arg1, const PossibleType &arg2) {
    // RDKit❗❌:   return (std::get<0>(arg1) < std::get<0>(arg2));
    // RDKit❗❌: };
    // Behavior: actual dictionaries own ring closure UInt values; no cached
    // ring-ID proxy guesses absence or swallows conversion errors. Source
    // atom/bond/ring/branch tokens, property clear/write order, unsigned rank
    // arithmetic then signed tuple ordering, masks, random draws and prefix
    // mutation are reproduced in one body for both detached carriers.
    // Safe state-index errors translate native invalid array/bitset access;
    // pointer identity/allocator errors and equal-key std::sort permutations
    // are not modeled/proven and retain the behavior gap for final review.
    // Complexity: same source full-size per-node seen bitmap and per-degree
    // possible vector plus recursion/sort. Vec<bool> uses bytes and find-first
    // scans bits individually rather than Boost packed-word scanning: this
    // known memory/search overhead gets a cost gap. Dict tree writes and old
    // SMILES rank-view conversion add costs, without cloning a whole graph.
    if atom >= graph.source_atom_count() {
        return Err(SmilesParseError::TraversalStateIndex {
            state: "atom",
            index: atom,
            count: graph.source_atom_count(),
        });
    }
    let mut seen = vec![false; graph.source_atom_count()];
    seen[atom] = true;
    stack.push(MolStackElem::Atom(atom));
    let color_count = colors.len();
    *colors
        .get_mut(atom)
        .ok_or(SmilesParseError::TraversalStateIndex {
            state: "colors",
            index: atom,
            count: color_count,
        })? = AtomColor::Grey;
    let mut traversal_order = Vec::new();
    if let Some(bond) = incoming_bond {
        traversal_order.push(bond);
    }
    let closures = atom_ring_closures
        .get(atom)
        .ok_or(SmilesParseError::TraversalStateIndex {
            state: "atomRingClosures",
            index: atom,
            count: atom_ring_closures.len(),
        })?;
    let mut closed = Vec::new();
    for &bond_id in closures {
        traversal_order.push(bond_id);
        let bond = graph.source_bond(bond_id)?;
        let other = other_atom(bond, atom)?;
        let seen_count = seen.len();
        *seen
            .get_mut(other)
            .ok_or(SmilesParseError::TraversalStateIndex {
                state: "seenFromHere",
                index: other,
                count: seen_count,
            })? = true;
        let ring = bond
            .prop("_TraversalRingClosureBond")
            .map(cosmolkit_core::property_value_to_uint)
            .transpose()
            .map_err(SmilesParseError::WriterNumeric)?;
        if let Some(ring) = ring {
            stack.push(MolStackElem::Bond {
                bond: bond_id,
                atom_to_left: atom,
            });
            stack.push(MolStackElem::Ring(ring as i32));
            closed.push(ring.wrapping_sub(1) as usize);
            let count = ring_ids.len();
            *ring_ids
                .get_mut(bond_id.index())
                .ok_or(SmilesParseError::TraversalStateIndex {
                    state: "ringIds",
                    index: bond_id.index(),
                    count,
                })? = Some(ring as usize);
        } else {
            let slot = available_ring_ids
                .iter()
                .position(|available| *available)
                .ok_or(SmilesParseError::TraversalTooManyOpenRings)?;
            available_ring_ids[slot] = false;
            let ring = (slot + 1) as u32;
            graph
                .source_bond_mut(bond_id)?
                .set_prop("_TraversalRingClosureBond", PropertyValue::UInt(ring))?;
            let count = ring_ids.len();
            *ring_ids
                .get_mut(bond_id.index())
                .ok_or(SmilesParseError::TraversalStateIndex {
                    state: "ringIds",
                    index: bond_id.index(),
                    count,
                })? = Some(ring as usize);
            let count = traversal_ring_closure_bonds.len();
            *traversal_ring_closure_bonds
                .get_mut(bond_id.index())
                .ok_or(SmilesParseError::TraversalStateIndex {
                    state: "traversalRingClosureBonds",
                    index: bond_id.index(),
                    count,
                })? = true;
            stack.push(MolStackElem::Ring(ring as i32));
        }
    }
    for slot in closed {
        let count = available_ring_ids.len();
        *available_ring_ids
            .get_mut(slot)
            .ok_or(SmilesParseError::TraversalStateIndex {
                state: "cyclesAvailable",
                index: slot,
                count,
            })? = true;
    }
    let mut possibles = Vec::with_capacity(graph.source_neighbors(atom).len());
    for (other, bond_id) in graph.source_neighbors(atom) {
        if let Some(mask) = bonds_in_play {
            if !*mask
                .get(bond_id.index())
                .ok_or(SmilesParseError::TraversalStateIndex {
                    state: "bondsInPlay",
                    index: bond_id.index(),
                    count: mask.len(),
                })?
            {
                continue;
            }
        }
        if Some(bond_id) == incoming_bond {
            continue;
        }
        let color = colors
            .get(other)
            .ok_or(SmilesParseError::TraversalStateIndex {
                state: "colors",
                index: other,
                count: colors.len(),
            })?;
        if *color != AtomColor::White || seen[other] {
            continue;
        }
        let mut rank = *ranks
            .get(other)
            .ok_or(SmilesParseError::TraversalStateIndex {
                state: "ranks",
                index: other,
                count: ranks.len(),
            })?;
        if let Some(random) = random_stream.as_deref_mut() {
            rank = random.next_u32();
        } else if *ring_bonds
            .get(bond_id.index())
            .ok_or(SmilesParseError::TraversalStateIndex {
                state: "ringInfo",
                index: bond_id.index(),
                count: ring_bonds.len(),
            })?
        {
            let salt = if let Some(symbols) = bond_symbols {
                (boost_hash_range(symbols.get(bond_id)?) % MAX_NATOMS as u32)
                    .wrapping_mul(MAX_NATOMS as u32)
                    .wrapping_mul(MAX_NATOMS as u32)
            } else {
                ((MAX_BONDTYPE - rdkit_bond_type_code(graph.source_bond(bond_id)?.order())) as u32)
                    .wrapping_mul(MAX_NATOMS as u32)
                    .wrapping_mul(MAX_NATOMS as u32)
            };
            rank = rank.wrapping_add(salt);
        }
        possibles.push(Possible {
            rank: i64::from(rank as i32),
            atom: other,
            bond: bond_id,
        });
    }
    possibles.sort_unstable_by_key(|possible| possible.rank);
    for (position, possible) in possibles.iter().copied().enumerate() {
        if colors[possible.atom] != AtomColor::White {
            continue;
        }
        graph.clear_source_traversal_order(possible.atom)?;
        traversal_order.push(possible.bond);
        let branch = position + 1 != possibles.len();
        if branch {
            stack.push(MolStackElem::BranchOpen(position as i32));
        }
        stack.push(MolStackElem::Bond {
            bond: possible.bond,
            atom_to_left: atom,
        });
        dfs_build_stack(
            graph,
            possible.atom,
            Some(possible.bond),
            colors,
            ranks,
            ring_bonds,
            atom_ring_closures,
            ring_ids,
            available_ring_ids,
            stack,
            atom_traversal_bond_order,
            traversal_ring_closure_bonds,
            bonds_in_play,
            bond_symbols,
            random_stream.as_deref_mut(),
        )?;
        if branch {
            stack.push(MolStackElem::BranchClose(position as i32));
        }
    }
    let count = atom_traversal_bond_order.len();
    *atom_traversal_bond_order
        .get_mut(atom)
        .ok_or(SmilesParseError::TraversalStateIndex {
            state: "atomTraversalBondOrder",
            index: atom,
            count,
        })? = traversal_order;
    colors[atom] = AtomColor::Black;
    Ok(())
}

/// Native overload derives atom participation only from actual selected bond
/// endpoints, then dispatches the one complete canonical fragment pipeline.
#[doc(hidden)]
#[allow(clippy::too_many_arguments)]
pub fn canonicalize_fragment_from_bond_mask_source(
    topology: &mut TopologyBlock,
    properties: &mut cosmolkit_model::MoleculeProperties,
    valence: &ValenceAssignment,
    rings: &mut cosmolkit_core::RingInfo,
    atom: usize,
    colors: &mut [AtomColor],
    ranks: &[u32],
    stack: &mut Vec<MolStackElem>,
    bonds_in_play: Option<&[bool]>,
    bond_symbols: Option<DfsBondSymbols<'_>>,
    isomeric_smiles: bool,
    random_stream: Option<&mut cosmolkit_core::RdkitRandomGenerator<'_>>,
    do_chiral_inversions: bool,
    groups: &BTreeMap<usize, usize>,
) -> Result<(), SmilesParseError> {
    // RDKit❗❌: void canonicalizeFragment(ROMol &mol, int atomIdx,
    // RDKit❗❌:                           std::vector<AtomColors> &colors,
    // RDKit❗❌:                           const UINT_VECT &ranks, MolStack &molStack,
    // RDKit❗❌:                           const boost::dynamic_bitset<> *bondsInPlay,
    // RDKit❗❌:                           const std::vector<std::string> *bondSymbols,
    // RDKit❗❌:                           bool doIsomericSmiles, bool doRandom,
    // RDKit❗❌:                           bool doChiralInversions) {
    // RDKit❗❌:   boost::dynamic_bitset<> atomsInPlay(mol.getNumAtoms());
    // RDKit❗❌:   if (!bondsInPlay) {
    // RDKit❗❌:     // if we weren't given a bondsInPlay, then all bonds are in play, so we need
    // RDKit❗❌:     // to set both those and the atomsInPlay here:
    // RDKit❗❌:     atomsInPlay.set();
    // RDKit❗❌:   } else {
    // RDKit❗❌:     for (const auto bnd : mol.bonds()) {
    // RDKit❗❌:       if ((*bondsInPlay)[bnd->getIdx()]) {
    // RDKit❗❌:         atomsInPlay.set(bnd->getBeginAtomIdx());
    // RDKit❗❌:         atomsInPlay.set(bnd->getEndAtomIdx());
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   canonicalizeFragment(mol, atomIdx, colors, ranks, molStack, &atomsInPlay,
    // RDKit❗❌:                        bondsInPlay, bondSymbols, doIsomericSmiles, doRandom,
    // RDKit❗❌:                        doChiralInversions);
    // RDKit❗❌: }
    // Same O(V+E) mask construction, no graph/property copy and no rank,
    // color, reachability, or bond-order heuristic. Byte mask versus packed
    // native bitset remains an explicit storage/word-operation cost gap.
    let n = topology.atoms.len();
    let mut atoms_in_play = vec![false; n];
    if let Some(mask) = bonds_in_play {
        for bond in &topology.bonds {
            let i = bond.id().index();
            if *mask.get(i).ok_or(SmilesParseError::TraversalStateIndex {
                state: "bondsInPlay",
                index: i,
                count: mask.len(),
            })? {
                for endpoint in [bond.begin().index(), bond.end().index()] {
                    *atoms_in_play.get_mut(endpoint).ok_or(
                        SmilesParseError::TraversalStateIndex {
                            state: "atomsInPlay",
                            index: endpoint,
                            count: n,
                        },
                    )? = true;
                }
            }
        }
    } else {
        atoms_in_play.fill(true);
    }
    canonicalize_fragment_source(
        topology,
        properties,
        valence,
        rings,
        atom,
        colors,
        ranks,
        stack,
        Some(&atoms_in_play),
        bonds_in_play,
        bond_symbols,
        isomeric_smiles,
        random_stream,
        do_chiral_inversions,
        groups,
    )
}

/// Canonical fragment traversal over actual detached molecule state. Query
/// carriers use the shared DFS; their source writer adapter is a separate
/// planned source function, not silently projected to ordinary element atoms.
#[doc(hidden)]
#[allow(clippy::too_many_arguments)]
pub fn canonicalize_fragment_source(
    topology: &mut TopologyBlock,
    properties: &mut cosmolkit_model::MoleculeProperties,
    valence: &ValenceAssignment,
    rings: &mut cosmolkit_core::RingInfo,
    atom: usize,
    colors: &mut [AtomColor],
    ranks: &[u32],
    stack: &mut Vec<MolStackElem>,
    atoms_in_play: Option<&[bool]>,
    bonds_in_play: Option<&[bool]>,
    bond_symbols: Option<DfsBondSymbols<'_>>,
    isomeric_smiles: bool,
    random_stream: Option<&mut cosmolkit_core::RdkitRandomGenerator<'_>>,
    do_chiral_inversions: bool,
    stereo_group_references: &BTreeMap<usize, usize>,
) -> Result<(), SmilesParseError> {
    canonicalize_fragment_impl(
        topology,
        properties,
        valence,
        rings,
        atom,
        colors,
        ranks,
        stack,
        atoms_in_play,
        bonds_in_play,
        bond_symbols,
        isomeric_smiles,
        random_stream,
        do_chiral_inversions,
        stereo_group_references,
        None,
    )
}

#[allow(clippy::too_many_arguments)]
fn canonicalize_fragment_impl(
    topology: &mut TopologyBlock,
    properties: &mut cosmolkit_model::MoleculeProperties,
    valence: &ValenceAssignment,
    rings: &mut cosmolkit_core::RingInfo,
    atom: usize,
    colors: &mut [AtomColor],
    ranks: &[u32],
    stack: &mut Vec<MolStackElem>,
    atoms_in_play: Option<&[bool]>,
    bonds_in_play: Option<&[bool]>,
    bond_symbols: Option<DfsBondSymbols<'_>>,
    isomeric_smiles: bool,
    random_stream: Option<&mut cosmolkit_core::RdkitRandomGenerator<'_>>,
    do_chiral_inversions: bool,
    stereo_group_references: &BTreeMap<usize, usize>,
    source_queries: Option<cosmolkit_model::QueryStateRef<'_>>,
) -> Result<(), SmilesParseError> {
    // RDKit❗❌: RDKIT_GRAPHMOL_EXPORT void canonicalizeFragment(
    // RDKit❗❌:     ROMol &mol, int atomIdx, std::vector<AtomColors> &colors,
    // RDKit❗❌:     const std::vector<unsigned int> &ranks, MolStack &molStack,
    // RDKit❗❌:     const boost::dynamic_bitset<> *atomsInPlay,
    // RDKit❗❌:     const boost::dynamic_bitset<> *bondsInPlay,
    // RDKit❗❌:     const std::vector<std::string> *bondSymbols, bool doIsomericSmiles,
    // RDKit❗❌:     bool doRandom, bool doChiralInversions) {
    // RDKit❗❌:   PRECONDITION(colors.size() >= mol.getNumAtoms(), "vector too small");
    // RDKit❗❌:   PRECONDITION(ranks.size() >= mol.getNumAtoms(), "vector too small");
    // RDKit❗❌:   PRECONDITION(!atomsInPlay || atomsInPlay->size() >= mol.getNumAtoms(),
    // RDKit❗❌:                "atomsInPlay too small");
    // RDKit❗❌:   PRECONDITION(!bondsInPlay || bondsInPlay->size() >= mol.getNumBonds(),
    // RDKit❗❌:                "bondsInPlay too small");
    // RDKit❗❌:   PRECONDITION(!bondSymbols || bondSymbols->size() >= mol.getNumBonds(),
    // RDKit❗❌:                "bondSymbols too small");
    // RDKit❗❌:   unsigned int nAtoms = mol.getNumAtoms();
    // RDKit❗❌:
    // RDKit❗❌:   // make sure that we've done the stereo perception:
    // RDKit❗❌:   if (!mol.hasProp(common_properties::_StereochemDone)) {
    // RDKit❗❌:     MolOps::assignStereochemistry(mol, false);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // we need ring information; make sure findSSSR has been called before
    // RDKit❗❌:   // if not call now
    // RDKit❗❌:   // NOTE: if called from the SMARTS code, the ring info will be set to SSSR,
    // RDKit❗❌:   // but no ring infor in actually set
    // RDKit❗❌:   if (!mol.getRingInfo()->isSymmSssr()) {
    // RDKit❗❌:     MolOps::findSSSR(mol);
    // RDKit❗❌:   }
    // RDKit❗❌:   mol.getAtomWithIdx(atomIdx)->setProp(common_properties::_TraversalStartPoint,
    // RDKit❗❌:                                        true);
    // RDKit❗❌:
    // RDKit❗❌:   VECT_INT_VECT atomRingClosures(nAtoms);
    // RDKit❗❌:   std::vector<INT_LIST> atomTraversalBondOrder(nAtoms);
    // RDKit❗❌:   Canon::canonicalDFSTraversal(mol, atomIdx, -1, colors, ranks, molStack,
    // RDKit❗❌:                                atomRingClosures, atomTraversalBondOrder,
    // RDKit❗❌:                                bondsInPlay, bondSymbols, doRandom);
    // RDKit❗❌:
    // RDKit❗❌:   CHECK_INVARIANT(!molStack.empty(), "Empty stack.");
    // RDKit❗❌:   CHECK_INVARIANT(molStack.begin()->type == MOL_STACK_ATOM,
    // RDKit❗❌:                   "Corrupted stack. First element should be an atom.");
    // RDKit❗❌:
    // RDKit❗❌:   // collect some information about traversal order on chiral atoms
    // RDKit❗❌:   boost::dynamic_bitset<> numSwapsChiralAtoms(nAtoms);
    // RDKit❗❌:   std::vector<int> atomPermutationIndices(nAtoms, 0);
    // RDKit❗❌:   if (doIsomericSmiles) {
    // RDKit❗❌:     for (const auto atom : mol.atoms()) {
    // RDKit❗❌:       if (atomsInPlay && !(*atomsInPlay)[atom->getIdx()]) {
    // RDKit❗❌:         continue;
    // RDKit❗❌:       }
    // RDKit❗❌:       if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
    // RDKit❗❌:         // check if all of this atom's bonds are in play
    // RDKit❗❌:         for (const auto bnd : mol.atomBonds(atom)) {
    // RDKit❗❌:           if (bondsInPlay && !(*bondsInPlay)[bnd->getIdx()]) {
    // RDKit❗❌:             atom->setProp(common_properties::_brokenChirality, true);
    // RDKit❗❌:             break;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:         if (atom->hasProp(common_properties::_brokenChirality)) {
    // RDKit❗❌:           continue;
    // RDKit❗❌:         }
    // RDKit❗❌:
    // RDKit❗❌:         // Extra check needed if/when @AL1/@AL2 supported
    // RDKit❗❌:         if (Chirality::detail::isAtomPotentialTetrahedralCenter(atom) ||
    // RDKit❗❌:             Chirality::hasNonTetrahedralStereo(atom)) {
    // RDKit❗❌:           int perm = 0;
    // RDKit❗❌:           if (Chirality::hasNonTetrahedralStereo(atom)) {
    // RDKit❗❌:             atom->getPropIfPresent(common_properties::_chiralPermutation, perm);
    // RDKit❗❌:           }
    // RDKit❗❌:
    // RDKit❗❌:           const unsigned int firstIdx = molStack.begin()->obj.atom->getIdx();
    // RDKit❗❌:           const bool firstInPart = atom->getIdx() == firstIdx;
    // RDKit❗❌:
    // RDKit❗❌:           // Check if the atom can be chiral, and if chirality needs inversion
    // RDKit❗❌:           const INT_LIST &trueOrder = atomTraversalBondOrder[atom->getIdx()];
    // RDKit❗❌:
    // RDKit❗❌:           // We have to make sure that trueOrder contains all the
    // RDKit❗❌:           // bonds, even if they won't be written to the SMILES
    // RDKit❗❌:           int nSwaps = 0;
    // RDKit❗❌:           if (trueOrder.size() < atom->getDegree()) {
    // RDKit❗❌:             INT_LIST tOrder = trueOrder;
    // RDKit❗❌:             for (const auto bnd : mol.atomBonds(atom)) {
    // RDKit❗❌:               int bndIdx = bnd->getIdx();
    // RDKit❗❌:               if (std::find(trueOrder.begin(), trueOrder.end(), bndIdx) ==
    // RDKit❗❌:                   trueOrder.end()) {
    // RDKit❗❌:                 tOrder.push_back(bndIdx);
    // RDKit❗❌:               }
    // RDKit❗❌:             }
    // RDKit❗❌:             if (!perm) {
    // RDKit❗❌:               nSwaps = atom->getPerturbationOrder(tOrder);
    // RDKit❗❌:             } else {
    // RDKit❗❌:               insertImplicitNbors(tOrder, atom->getChiralTag(), firstInPart);
    // RDKit❗❌:               perm = Chirality::getChiralPermutation(atom, tOrder);
    // RDKit❗❌:             }
    // RDKit❗❌:           } else {
    // RDKit❗❌:             if (!perm) {
    // RDKit❗❌:               nSwaps = atom->getPerturbationOrder(trueOrder);
    // RDKit❗❌:             } else {
    // RDKit❗❌:               INT_LIST tOrder = trueOrder;
    // RDKit❗❌:               insertImplicitNbors(tOrder, atom->getChiralTag(), firstInPart);
    // RDKit❗❌:               perm = Chirality::getChiralPermutation(atom, tOrder);
    // RDKit❗❌:             }
    // RDKit❗❌:           }
    // RDKit❗❌:
    // RDKit❗❌:           // in future this should be moved up and simplified, there should not
    // RDKit❗❌:           // be an option to not do chiral inversions
    // RDKit❗❌:           if (doChiralInversions &&
    // RDKit❗❌:               chiralAtomNeedsTagInversion(
    // RDKit❗❌:                   mol, atom, firstInPart,
    // RDKit❗❌:                   atomRingClosures[atom->getIdx()].size())) {
    // RDKit❗❌:             // This is a special case. Here's an example:
    // RDKit❗❌:             //   Our internal representation of a chiral center is equivalent
    // RDKit❗❌:             //   to:
    // RDKit❗❌:             //     [C@](F)(O)(C)[H]
    // RDKit❗❌:             //   we'll be dumping it without the H, which entails a
    // RDKit❗❌:             //   reordering:
    // RDKit❗❌:             //     [C@@H](F)(O)C
    // RDKit❗❌:             ++nSwaps;
    // RDKit❗❌:           }
    // RDKit❗❌:           if (nSwaps % 2) {
    // RDKit❗❌:             numSwapsChiralAtoms.set(atom->getIdx());
    // RDKit❗❌:           }
    // RDKit❗❌:           atomPermutationIndices[atom->getIdx()] = perm;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<unsigned int> atomVisitOrders(mol.getNumAtoms());
    // RDKit❗❌:   std::vector<unsigned int> bondVisitOrders(mol.getNumBonds());
    // RDKit❗❌:
    // RDKit❗❌:   unsigned int pos = 0;
    // RDKit❗❌:   for (const auto &msI : molStack) {
    // RDKit❗❌:     if (msI.type == MOL_STACK_ATOM) {
    // RDKit❗❌:       atomVisitOrders[msI.obj.atom->getIdx()] = pos;
    // RDKit❗❌:     } else if (msI.type == MOL_STACK_BOND) {
    // RDKit❗❌:       bondVisitOrders[msI.obj.bond->getIdx()] = pos;
    // RDKit❗❌:       auto dir = msI.obj.bond->getBondDir();
    // RDKit❗❌:       if (dir == Bond::ENDDOWNRIGHT || dir == Bond::ENDUPRIGHT) {
    // RDKit❗❌:         msI.obj.bond->setBondDir(Bond::NONE);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     ++pos;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<int8_t> bondDirCounts(mol.getNumBonds(), 0);
    // RDKit❗❌:   std::vector<int8_t> atomDirCounts(nAtoms, 0);
    // RDKit❗❌:   canonicalizeDoubleBonds(mol, bondVisitOrders, atomVisitOrders, bondDirCounts,
    // RDKit❗❌:                           atomDirCounts, molStack);
    // RDKit❗❌:
    // RDKit❗❌:   // traverse the stack and canonicalize atoms with (ring) stereochemistry
    // RDKit❗❌:   if (doIsomericSmiles) {
    // RDKit❗❌:     boost::dynamic_bitset<> ringStereoChemAdjusted(nAtoms);
    // RDKit❗❌:     for (auto &msI : molStack) {
    // RDKit❗❌:       if (msI.type == MOL_STACK_ATOM &&
    // RDKit❗❌:           msI.obj.atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit❗❌:           !msI.obj.atom->hasProp(common_properties::_brokenChirality)) {
    // RDKit❗❌:         if (msI.obj.atom->hasProp(common_properties::_ringStereoAtoms)) {
    // RDKit❗❌:           // FIX: handle stereogroups here too
    // RDKit❗❌:           if (!ringStereoChemAdjusted[msI.obj.atom->getIdx()]) {
    // RDKit❗❌:             msI.obj.atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
    // RDKit❗❌:             ringStereoChemAdjusted.set(msI.obj.atom->getIdx());
    // RDKit❗❌:           }
    // RDKit❗❌:           const INT_VECT &ringStereoAtoms = msI.obj.atom->getProp<INT_VECT>(
    // RDKit❗❌:               common_properties::_ringStereoAtoms);
    // RDKit❗❌:           for (auto nbrV : ringStereoAtoms) {
    // RDKit❗❌:             int nbrIdx = abs(nbrV) - 1;
    // RDKit❗❌:             // Adjust the chirality flag of the ring stereo atoms according to
    // RDKit❗❌:             // the first one
    // RDKit❗❌:             if (!ringStereoChemAdjusted[nbrIdx] &&
    // RDKit❗❌:                 atomVisitOrders[nbrIdx] >
    // RDKit❗❌:                     atomVisitOrders[msI.obj.atom->getIdx()]) {
    // RDKit❗❌:               mol.getAtomWithIdx(nbrIdx)->setChiralTag(
    // RDKit❗❌:                   msI.obj.atom->getChiralTag());
    // RDKit❗❌:               if (nbrV < 0) {
    // RDKit❗❌:                 mol.getAtomWithIdx(nbrIdx)->invertChirality();
    // RDKit❗❌:               }
    // RDKit❗❌:               // Odd number of swaps for first chiral ring atom --> needs to be
    // RDKit❗❌:               // swapped but we want to retain chirality
    // RDKit❗❌:               if (numSwapsChiralAtoms[msI.obj.atom->getIdx()]) {
    // RDKit❗❌:                 // Odd number of swaps for chiral ring neighbor --> needs to be
    // RDKit❗❌:                 // swapped but we want to retain chirality
    // RDKit❗❌:                 if (!numSwapsChiralAtoms[nbrIdx]) {
    // RDKit❗❌:                   mol.getAtomWithIdx(nbrIdx)->invertChirality();
    // RDKit❗❌:                 }
    // RDKit❗❌:               }
    // RDKit❗❌:               // Even number of swaps for first chiral ring atom --> don't need
    // RDKit❗❌:               // to be swapped
    // RDKit❗❌:               else {
    // RDKit❗❌:                 // Odd number of swaps for chiral ring neighbor --> needs to be
    // RDKit❗❌:                 // swapped
    // RDKit❗❌:                 if (numSwapsChiralAtoms[nbrIdx]) {
    // RDKit❗❌:                   mol.getAtomWithIdx(nbrIdx)->invertChirality();
    // RDKit❗❌:                 }
    // RDKit❗❌:               }
    // RDKit❗❌:               ringStereoChemAdjusted.set(nbrIdx);
    // RDKit❗❌:             }
    // RDKit❗❌:           }
    // RDKit❗❌:         } else if (size_t sgidx;
    // RDKit❗❌:                    msI.obj.atom->getPropIfPresent("_stereoGroup", sgidx) &&
    // RDKit❗❌:                    mol.getStereoGroups().size() > sgidx) {
    // RDKit❗❌:           // make sure that the reference atom in the stereogroup is CCW
    // RDKit❗❌:           auto &sg = mol.getStereoGroups()[sgidx];
    // RDKit❗❌:           bool swapIt =
    // RDKit❗❌:               msI.obj.atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW;
    // RDKit❗❌:           if (swapIt) {
    // RDKit❗❌:             msI.obj.atom->invertChirality();
    // RDKit❗❌:           }
    // RDKit❗❌:           if (swapIt || numSwapsChiralAtoms[msI.obj.atom->getIdx()]) {
    // RDKit❗❌:             for (auto at : sg.getAtoms()) {
    // RDKit❗❌:               if (at == msI.obj.atom) {
    // RDKit❗❌:                 continue;
    // RDKit❗❌:               }
    // RDKit❗❌:               at->invertChirality();
    // RDKit❗❌:             }
    // RDKit❗❌:           }
    // RDKit❗❌:
    // RDKit❗❌:         } else {
    // RDKit❗❌:           if (msI.obj.atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit❗❌:               msI.obj.atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit❗❌:             if ((numSwapsChiralAtoms[msI.obj.atom->getIdx()])) {
    // RDKit❗❌:               msI.obj.atom->invertChirality();
    // RDKit❗❌:             }
    // RDKit❗❌:           } else if (atomPermutationIndices[msI.obj.atom->getIdx()]) {
    // RDKit❗❌:             msI.obj.atom->setProp(
    // RDKit❗❌:                 common_properties::_chiralPermutation,
    // RDKit❗❌:                 atomPermutationIndices[msI.obj.atom->getIdx()]);
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   Canon::removeUnwantedBondDirSpecs(mol, molStack, bondDirCounts, atomDirCounts,
    // RDKit❗❌:                                     bondVisitOrders);
    // RDKit❗❌:
    // RDKit❗❌:   Canon::removeRedundantBondDirSpecs(mol, molStack, bondDirCounts,
    // RDKit❗❌:                                      atomDirCounts);
    // RDKit❗❌:
    // RDKit❗❌: #if ENABLE_EXTRA_CHECKS
    // RDKit❗❌:   checkDirCounts(mol, bondDirCounts, atomDirCounts);
    // RDKit❗❌: #endif
    // RDKit❗❌: }
    // The source pipeline uses the one borrowed CORE stereo/ring owner and
    // shared notation DFS and stack phases. Preparation errors retain actual
    // property/topology/cache prefix mutations. No whole graph clone or
    // synthetic cache is used. Native comparator permutations, existing helper
    // comparison gaps, packed-bitset costs, and query adapter work remain ❗/❌.
    let n = topology.atoms.len();
    let m = topology.bonds.len();
    for (name, len, min) in [("colors", colors.len(), n), ("ranks", ranks.len(), n)] {
        if len < min {
            return Err(SmilesParseError::TraversalStateIndex {
                state: name,
                index: min,
                count: len,
            });
        }
    }
    for (name, value, min) in [
        ("atomsInPlay", atoms_in_play.map(<[bool]>::len), n),
        ("bondsInPlay", bonds_in_play.map(<[bool]>::len), m),
        ("bondSymbols", bond_symbols.map(DfsBondSymbols::len), m),
    ] {
        if let Some(len) = value {
            if len < min {
                return Err(SmilesParseError::TraversalStateIndex {
                    state: name,
                    index: min,
                    count: len,
                });
            }
        }
    }
    if properties.prop("_StereochemDone").is_none() {
        let mut update = None;
        let mut atom_valence_updates = Vec::new();
        let result = cosmolkit_core::assign_legacy_stereochemistry_source(
            topology,
            valence,
            rings,
            source_queries,
            false,
            false,
            &mut update,
            &mut atom_valence_updates,
        );
        // cleanIt=false cannot enter the explicit-H removal branch.
        debug_assert!(atom_valence_updates.is_empty());
        if let Some(update) = update {
            *rings = update;
        }
        result.map_err(SmilesParseError::WriterLegacyStereo)?;
        properties.set_computed_prop("_StereochemDone", PropertyValue::Int(1))?;
    }
    if !rings.is_symm_sssr() {
        cosmolkit_core::find_sssr_with_source_outputs_from_parts(
            n,
            &topology.bonds,
            &topology.adjacency,
            rings,
            Some(properties),
            None,
            false,
            false,
        )
        .map_err(SmilesParseError::WriterRings)?;
    }
    let start = topology
        .atoms
        .get_mut(atom)
        .ok_or(SmilesParseError::TraversalStateIndex {
            state: "atom",
            index: atom,
            count: n,
        })?;
    start.set_prop("_TraversalStartPoint", true)?;
    let mut closures = vec![vec![]; n];
    let mut orders = vec![vec![]; n];
    let ring_bonds = (0..m)
        .map(|i| rings.num_bond_rings(BondId::new(i)) != 0)
        .collect::<Vec<_>>();
    let opened = canonical_dfs_traversal_source(
        topology,
        atom,
        None,
        colors,
        ranks,
        &ring_bonds,
        stack,
        &mut closures,
        &mut orders,
        bonds_in_play,
        bond_symbols,
        random_stream,
    )?;
    canonicalize_source_stack(
        topology,
        valence,
        Some(rings),
        stack,
        &closures,
        &orders,
        &opened,
        atoms_in_play,
        bonds_in_play,
        isomeric_smiles,
        do_chiral_inversions,
        stereo_group_references,
        source_queries,
    )
}

fn source_atom_permutation(atom: &Atom) -> Result<i32, SmilesParseError> {
    // RDKit❗✔️: atom->getPropIfPresent(common_properties::_chiralPermutation, perm);
    // Raw dictionary data takes precedence over the explicitly modeled source
    // fact. A modeled UInt follows the same source checked int conversion;
    // only actual property absence produces the caller-initialized zero.
    if let Some(value) = atom.prop("_chiralPermutation") {
        return cosmolkit_core::property_value_to_int(value).map_err(SmilesParseError::WriterInt);
    }
    match atom.chiral_permutation() {
        Some(value) => cosmolkit_core::property_value_to_int(&PropertyValue::UInt(value))
            .map_err(SmilesParseError::WriterInt),
        None => Ok(0),
    }
}

#[allow(clippy::too_many_arguments)]
fn canonicalize_source_stack(
    topology: &mut TopologyBlock,
    valence: &ValenceAssignment,
    rings: Option<&cosmolkit_core::RingInfo>,
    stack: &[MolStackElem],
    closures: &[Vec<BondId>],
    orders: &[Vec<BondId>],
    opened: &[bool],
    atoms_in_play: Option<&[bool]>,
    bonds_in_play: Option<&[bool]>,
    isomeric_smiles: bool,
    do_chiral_inversions: bool,
    stereo_group_references: &BTreeMap<usize, usize>,
    source_queries: Option<cosmolkit_model::QueryStateRef<'_>>,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKit 2026.03.6 COMPLETE Canon::canonicalizeFragment
    // RDKit❗❌: RDKIT_GRAPHMOL_EXPORT void canonicalizeFragment(
    // RDKit❗❌:     ROMol &mol, int atomIdx, std::vector<AtomColors> &colors,
    // RDKit❗❌:     const std::vector<unsigned int> &ranks, MolStack &molStack,
    // RDKit❗❌:     const boost::dynamic_bitset<> *atomsInPlay,
    // RDKit❗❌:     const boost::dynamic_bitset<> *bondsInPlay,
    // RDKit❗❌:     const std::vector<std::string> *bondSymbols, bool doIsomericSmiles,
    // RDKit❗❌:     bool doRandom, bool doChiralInversions) {
    // RDKit❗❌:   PRECONDITION(colors.size() >= mol.getNumAtoms(), "vector too small");
    // RDKit❗❌:   PRECONDITION(ranks.size() >= mol.getNumAtoms(), "vector too small");
    // RDKit❗❌:   PRECONDITION(!atomsInPlay || atomsInPlay->size() >= mol.getNumAtoms(),
    // RDKit❗❌:                "atomsInPlay too small");
    // RDKit❗❌:   PRECONDITION(!bondsInPlay || bondsInPlay->size() >= mol.getNumBonds(),
    // RDKit❗❌:                "bondsInPlay too small");
    // RDKit❗❌:   PRECONDITION(!bondSymbols || bondSymbols->size() >= mol.getNumBonds(),
    // RDKit❗❌:                "bondSymbols too small");
    // RDKit❗❌:   unsigned int nAtoms = mol.getNumAtoms();
    // RDKit❗❌:
    // RDKit❗❌:   // make sure that we've done the stereo perception:
    // RDKit❗❌:   if (!mol.hasProp(common_properties::_StereochemDone)) {
    // RDKit❗❌:     MolOps::assignStereochemistry(mol, false);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // we need ring information; make sure findSSSR has been called before
    // RDKit❗❌:   // if not call now
    // RDKit❗❌:   // NOTE: if called from the SMARTS code, the ring info will be set to SSSR,
    // RDKit❗❌:   // but no ring infor in actually set
    // RDKit❗❌:   if (!mol.getRingInfo()->isSymmSssr()) {
    // RDKit❗❌:     MolOps::findSSSR(mol);
    // RDKit❗❌:   }
    // RDKit❗❌:   mol.getAtomWithIdx(atomIdx)->setProp(common_properties::_TraversalStartPoint,
    // RDKit❗❌:                                        true);
    // RDKit❗❌:
    // RDKit❗❌:   VECT_INT_VECT atomRingClosures(nAtoms);
    // RDKit❗❌:   std::vector<INT_LIST> atomTraversalBondOrder(nAtoms);
    // RDKit❗❌:   Canon::canonicalDFSTraversal(mol, atomIdx, -1, colors, ranks, molStack,
    // RDKit❗❌:                                atomRingClosures, atomTraversalBondOrder,
    // RDKit❗❌:                                bondsInPlay, bondSymbols, doRandom);
    // RDKit❗❌:
    // RDKit❗❌:   CHECK_INVARIANT(!molStack.empty(), "Empty stack.");
    // RDKit❗❌:   CHECK_INVARIANT(molStack.begin()->type == MOL_STACK_ATOM,
    // RDKit❗❌:                   "Corrupted stack. First element should be an atom.");
    // RDKit❗❌:
    // RDKit❗❌:   // collect some information about traversal order on chiral atoms
    // RDKit❗❌:   boost::dynamic_bitset<> numSwapsChiralAtoms(nAtoms);
    // RDKit❗❌:   std::vector<int> atomPermutationIndices(nAtoms, 0);
    // RDKit❗❌:   if (doIsomericSmiles) {
    // RDKit❗❌:     for (const auto atom : mol.atoms()) {
    // RDKit❗❌:       if (atomsInPlay && !(*atomsInPlay)[atom->getIdx()]) {
    // RDKit❗❌:         continue;
    // RDKit❗❌:       }
    // RDKit❗❌:       if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
    // RDKit❗❌:         // check if all of this atom's bonds are in play
    // RDKit❗❌:         for (const auto bnd : mol.atomBonds(atom)) {
    // RDKit❗❌:           if (bondsInPlay && !(*bondsInPlay)[bnd->getIdx()]) {
    // RDKit❗❌:             atom->setProp(common_properties::_brokenChirality, true);
    // RDKit❗❌:             break;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:         if (atom->hasProp(common_properties::_brokenChirality)) {
    // RDKit❗❌:           continue;
    // RDKit❗❌:         }
    // RDKit❗❌:
    // RDKit❗❌:         // Extra check needed if/when @AL1/@AL2 supported
    // RDKit❗❌:         if (Chirality::detail::isAtomPotentialTetrahedralCenter(atom) ||
    // RDKit❗❌:             Chirality::hasNonTetrahedralStereo(atom)) {
    // RDKit❗❌:           int perm = 0;
    // RDKit❗❌:           if (Chirality::hasNonTetrahedralStereo(atom)) {
    // RDKit❗❌:             atom->getPropIfPresent(common_properties::_chiralPermutation, perm);
    // RDKit❗❌:           }
    // RDKit❗❌:
    // RDKit❗❌:           const unsigned int firstIdx = molStack.begin()->obj.atom->getIdx();
    // RDKit❗❌:           const bool firstInPart = atom->getIdx() == firstIdx;
    // RDKit❗❌:
    // RDKit❗❌:           // Check if the atom can be chiral, and if chirality needs inversion
    // RDKit❗❌:           const INT_LIST &trueOrder = atomTraversalBondOrder[atom->getIdx()];
    // RDKit❗❌:
    // RDKit❗❌:           // We have to make sure that trueOrder contains all the
    // RDKit❗❌:           // bonds, even if they won't be written to the SMILES
    // RDKit❗❌:           int nSwaps = 0;
    // RDKit❗❌:           if (trueOrder.size() < atom->getDegree()) {
    // RDKit❗❌:             INT_LIST tOrder = trueOrder;
    // RDKit❗❌:             for (const auto bnd : mol.atomBonds(atom)) {
    // RDKit❗❌:               int bndIdx = bnd->getIdx();
    // RDKit❗❌:               if (std::find(trueOrder.begin(), trueOrder.end(), bndIdx) ==
    // RDKit❗❌:                   trueOrder.end()) {
    // RDKit❗❌:                 tOrder.push_back(bndIdx);
    // RDKit❗❌:               }
    // RDKit❗❌:             }
    // RDKit❗❌:             if (!perm) {
    // RDKit❗❌:               nSwaps = atom->getPerturbationOrder(tOrder);
    // RDKit❗❌:             } else {
    // RDKit❗❌:               insertImplicitNbors(tOrder, atom->getChiralTag(), firstInPart);
    // RDKit❗❌:               perm = Chirality::getChiralPermutation(atom, tOrder);
    // RDKit❗❌:             }
    // RDKit❗❌:           } else {
    // RDKit❗❌:             if (!perm) {
    // RDKit❗❌:               nSwaps = atom->getPerturbationOrder(trueOrder);
    // RDKit❗❌:             } else {
    // RDKit❗❌:               INT_LIST tOrder = trueOrder;
    // RDKit❗❌:               insertImplicitNbors(tOrder, atom->getChiralTag(), firstInPart);
    // RDKit❗❌:               perm = Chirality::getChiralPermutation(atom, tOrder);
    // RDKit❗❌:             }
    // RDKit❗❌:           }
    // RDKit❗❌:
    // RDKit❗❌:           // in future this should be moved up and simplified, there should not
    // RDKit❗❌:           // be an option to not do chiral inversions
    // RDKit❗❌:           if (doChiralInversions &&
    // RDKit❗❌:               chiralAtomNeedsTagInversion(
    // RDKit❗❌:                   mol, atom, firstInPart,
    // RDKit❗❌:                   atomRingClosures[atom->getIdx()].size())) {
    // RDKit❗❌:             // This is a special case. Here's an example:
    // RDKit❗❌:             //   Our internal representation of a chiral center is equivalent
    // RDKit❗❌:             //   to:
    // RDKit❗❌:             //     [C@](F)(O)(C)[H]
    // RDKit❗❌:             //   we'll be dumping it without the H, which entails a
    // RDKit❗❌:             //   reordering:
    // RDKit❗❌:             //     [C@@H](F)(O)C
    // RDKit❗❌:             ++nSwaps;
    // RDKit❗❌:           }
    // RDKit❗❌:           if (nSwaps % 2) {
    // RDKit❗❌:             numSwapsChiralAtoms.set(atom->getIdx());
    // RDKit❗❌:           }
    // RDKit❗❌:           atomPermutationIndices[atom->getIdx()] = perm;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<unsigned int> atomVisitOrders(mol.getNumAtoms(), 0);
    // RDKit❗❌:   std::vector<unsigned int> bondVisitOrders(mol.getNumBonds(), 0);
    // RDKit❗❌:
    // RDKit❗❌:   unsigned int pos =
    // RDKit❗❌:       1;  // start at 1 since we use 0 to detect unvisited atoms/bonds
    // RDKit❗❌:   for (const auto &msI : molStack) {
    // RDKit❗❌:     if (msI.type == MOL_STACK_ATOM) {
    // RDKit❗❌:       atomVisitOrders[msI.obj.atom->getIdx()] = pos;
    // RDKit❗❌:     } else if (msI.type == MOL_STACK_BOND) {
    // RDKit❗❌:       bondVisitOrders[msI.obj.bond->getIdx()] = pos;
    // RDKit❗❌:       auto dir = msI.obj.bond->getBondDir();
    // RDKit❗❌:       if (dir == Bond::ENDDOWNRIGHT || dir == Bond::ENDUPRIGHT) {
    // RDKit❗❌:         msI.obj.bond->setBondDir(Bond::NONE);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     ++pos;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<int8_t> bondDirCounts(mol.getNumBonds(), 0);
    // RDKit❗❌:   std::vector<int8_t> atomDirCounts(nAtoms, 0);
    // RDKit❗❌:   canonicalizeDoubleBonds(mol, bondVisitOrders, atomVisitOrders, bondDirCounts,
    // RDKit❗❌:                           atomDirCounts, molStack);
    // RDKit❗❌:
    // RDKit❗❌:   // traverse the stack and canonicalize atoms with (ring) stereochemistry
    // RDKit❗❌:   if (doIsomericSmiles) {
    // RDKit❗❌:     boost::dynamic_bitset<> ringStereoChemAdjusted(nAtoms);
    // RDKit❗❌:     for (auto &msI : molStack) {
    // RDKit❗❌:       if (msI.type == MOL_STACK_ATOM &&
    // RDKit❗❌:           msI.obj.atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit❗❌:           !msI.obj.atom->hasProp(common_properties::_brokenChirality)) {
    // RDKit❗❌:         if (msI.obj.atom->hasProp(common_properties::_ringStereoAtoms)) {
    // RDKit❗❌:           // FIX: handle stereogroups here too
    // RDKit❗❌:           if (!ringStereoChemAdjusted[msI.obj.atom->getIdx()]) {
    // RDKit❗❌:             msI.obj.atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
    // RDKit❗❌:             ringStereoChemAdjusted.set(msI.obj.atom->getIdx());
    // RDKit❗❌:           }
    // RDKit❗❌:           const INT_VECT &ringStereoAtoms = msI.obj.atom->getProp<INT_VECT>(
    // RDKit❗❌:               common_properties::_ringStereoAtoms);
    // RDKit❗❌:           for (auto nbrV : ringStereoAtoms) {
    // RDKit❗❌:             int nbrIdx = abs(nbrV) - 1;
    // RDKit❗❌:             // Adjust the chirality flag of the ring stereo atoms according to
    // RDKit❗❌:             // the first one
    // RDKit❗❌:             if (!ringStereoChemAdjusted[nbrIdx] &&
    // RDKit❗❌:                 atomVisitOrders[nbrIdx] >
    // RDKit❗❌:                     atomVisitOrders[msI.obj.atom->getIdx()]) {
    // RDKit❗❌:               mol.getAtomWithIdx(nbrIdx)->setChiralTag(
    // RDKit❗❌:                   msI.obj.atom->getChiralTag());
    // RDKit❗❌:               if (nbrV < 0) {
    // RDKit❗❌:                 mol.getAtomWithIdx(nbrIdx)->invertChirality();
    // RDKit❗❌:               }
    // RDKit❗❌:               // Odd number of swaps for first chiral ring atom --> needs to be
    // RDKit❗❌:               // swapped but we want to retain chirality
    // RDKit❗❌:               if (numSwapsChiralAtoms[msI.obj.atom->getIdx()]) {
    // RDKit❗❌:                 // Odd number of swaps for chiral ring neighbor --> needs to be
    // RDKit❗❌:                 // swapped but we want to retain chirality
    // RDKit❗❌:                 if (!numSwapsChiralAtoms[nbrIdx]) {
    // RDKit❗❌:                   mol.getAtomWithIdx(nbrIdx)->invertChirality();
    // RDKit❗❌:                 }
    // RDKit❗❌:               }
    // RDKit❗❌:               // Even number of swaps for first chiral ring atom --> don't need
    // RDKit❗❌:               // to be swapped
    // RDKit❗❌:               else {
    // RDKit❗❌:                 // Odd number of swaps for chiral ring neighbor --> needs to be
    // RDKit❗❌:                 // swapped
    // RDKit❗❌:                 if (numSwapsChiralAtoms[nbrIdx]) {
    // RDKit❗❌:                   mol.getAtomWithIdx(nbrIdx)->invertChirality();
    // RDKit❗❌:                 }
    // RDKit❗❌:               }
    // RDKit❗❌:               ringStereoChemAdjusted.set(nbrIdx);
    // RDKit❗❌:             }
    // RDKit❗❌:           }
    // RDKit❗❌:         } else if (size_t sgidx;
    // RDKit❗❌:                    msI.obj.atom->getPropIfPresent("_stereoGroup", sgidx) &&
    // RDKit❗❌:                    mol.getStereoGroups().size() > sgidx) {
    // RDKit❗❌:           // make sure that the reference atom in the stereogroup is CCW
    // RDKit❗❌:           auto &sg = mol.getStereoGroups()[sgidx];
    // RDKit❗❌:           bool swapIt =
    // RDKit❗❌:               msI.obj.atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW;
    // RDKit❗❌:           if (swapIt) {
    // RDKit❗❌:             msI.obj.atom->invertChirality();
    // RDKit❗❌:           }
    // RDKit❗❌:           if (swapIt || numSwapsChiralAtoms[msI.obj.atom->getIdx()]) {
    // RDKit❗❌:             for (auto at : sg.getAtoms()) {
    // RDKit❗❌:               if (at == msI.obj.atom) {
    // RDKit❗❌:                 continue;
    // RDKit❗❌:               }
    // RDKit❗❌:               at->invertChirality();
    // RDKit❗❌:             }
    // RDKit❗❌:           }
    // RDKit❗❌:
    // RDKit❗❌:         } else {
    // RDKit❗❌:           if (msI.obj.atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit❗❌:               msI.obj.atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit❗❌:             if ((numSwapsChiralAtoms[msI.obj.atom->getIdx()])) {
    // RDKit❗❌:               msI.obj.atom->invertChirality();
    // RDKit❗❌:             }
    // RDKit❗❌:           } else if (atomPermutationIndices[msI.obj.atom->getIdx()]) {
    // RDKit❗❌:             msI.obj.atom->setProp(
    // RDKit❗❌:                 common_properties::_chiralPermutation,
    // RDKit❗❌:                 atomPermutationIndices[msI.obj.atom->getIdx()]);
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   Canon::removeUnwantedBondDirSpecs(mol, molStack, bondDirCounts, atomDirCounts,
    // RDKit❗❌:                                     bondVisitOrders);
    // RDKit❗❌:
    // RDKit❗❌:   Canon::removeRedundantBondDirSpecs(mol, molStack, bondDirCounts,
    // RDKit❗❌:                                      atomDirCounts);
    // RDKit❗❌:
    // RDKit❗❌: #if ENABLE_EXTRA_CHECKS
    // RDKit❗❌:   checkDirCounts(mol, bondDirCounts, atomDirCounts);
    // RDKit❗❌: #endif
    // RDKit❗❌: }
    // END RDKit 2026.03.6 COMPLETE Canon::canonicalizeFragment
    // The shared phases mutate the actual source fields, in source order.
    // O(V+E) scratch and degree-local perturbation work; source INT_LIST views
    // need Vec projections and the byte bitmap costs are retained as ❌.
    // ENABLE_EXTRA_CHECKS is pinned to 0 in Canon.cpp: no check is synthesized.
    let first = match stack.first() {
        Some(MolStackElem::Atom(i)) => *i,
        _ => {
            return Err(SmilesParseError::WriterCanonicalInvariant(
                "Empty or corrupted stack. First element should be an atom.",
            ));
        }
    };
    let n = topology.atoms.len();
    let m = topology.bonds.len();
    let mut swaps = vec![false; n];
    let mut permutations = vec![0i32; n];
    if isomeric_smiles {
        for i in 0..n {
            if atoms_in_play.is_some_and(|mask| !mask[i])
                || topology.atoms[i].chiral_tag() == ChiralTag::Unspecified
            {
                continue;
            }
            if topology
                .adjacency
                .neighbors_of(i)
                .iter()
                .any(|nb| bonds_in_play.is_some_and(|mask| !mask[nb.bond.index()]))
            {
                topology.atoms[i].set_prop("_brokenChirality", true)?;
            }
            if topology.atoms[i].prop("_brokenChirality").is_some() {
                continue;
            }
            // RDKit❗✔️: bool hasNonTetrahedralStereo(const Atom *cen) {
            // RDKit❗✔️:   PRECONDITION(cen, "bad center pointer");
            // RDKit❗✔️:   if (!cen->hasOwningMol()) {
            // RDKit❗✔️:     return false;
            // RDKit❗✔️:   }
            // RDKit❗✔️:   auto tag = cen->getChiralTag();
            // RDKit❗✔️:   return tag == Atom::ChiralType::CHI_SQUAREPLANAR ||
            // RDKit❗✔️:          tag == Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL ||
            // RDKit❗✔️:          tag == Atom::ChiralType::CHI_OCTAHEDRAL;
            // RDKit❗✔️: }
            // These are actual graph members of the canonicalization input;
            // source owning-molecule membership is supplied by that boundary.
            let non_tetra = matches!(
                topology.atoms[i].chiral_tag(),
                ChiralTag::SquarePlanar | ChiralTag::TrigonalBipyramidal | ChiralTag::Octahedral
            );
            if !cosmolkit_core::potential_tetrahedral_center_from_source(
                topology,
                valence,
                rings,
                AtomId::new(i),
            )
            .map_err(SmilesParseError::WriterPotentialStereo)?
                && !non_tetra
            {
                continue;
            }
            let mut perm = if non_tetra {
                source_atom_permutation(&topology.atoms[i])?
            } else {
                0
            };
            let true_order = &orders[i];
            let incident = topology.adjacency.neighbors_of(i);
            let mut padded = Vec::new();
            let traversal = if true_order.len() < incident.len() {
                padded.extend_from_slice(true_order);
                for nb in incident {
                    if !true_order.contains(&nb.bond) {
                        padded.push(nb.bond);
                    }
                }
                &padded[..]
            } else {
                &true_order[..]
            };
            let mut count = 0i32;
            if perm == 0 {
                let probe = traversal
                    .iter()
                    .map(|b| b.index() as u32 as i32)
                    .collect::<Vec<_>>();
                count = cosmolkit_core::atom_perturbation_order(
                    &probe,
                    incident.iter().map(|nb| nb.bond.index()),
                )
                .map_err(SmilesParseError::WriterStereoOrder)?;
            } else {
                let mut probe = traversal.iter().copied().map(Some).collect::<Vec<_>>();
                stereo::insert_implicit_nontetrahedral_neighbors(
                    &mut probe,
                    topology.atoms[i].chiral_tag(),
                    i == first,
                );
                let incident = incident.iter().map(|nb| nb.bond).collect::<Vec<_>>();
                // getChiralPermutation returns zero for a non-positive source int.
                perm = if perm <= 0 {
                    0
                } else {
                    stereo::nontetrahedral_chiral_permutation(
                        perm as u32,
                        topology.atoms[i].chiral_tag(),
                        m,
                        &incident,
                        &probe,
                        false,
                    )
                    .map_err(SmilesParseError::WriterParserStereoOrder)? as i32
                };
            }
            let atom = &topology.atoms[i];
            if do_chiral_inversions
                && stereo::chiral_atom_needs_tag_inversion(
                    incident.len(),
                    atom.explicit_hydrogens(),
                    i == first,
                    stereo::atom_has_fourth_valence(
                        atom.explicit_hydrogens(),
                        valence.implicit_hydrogens[i] == 1,
                    ) || source_queries.is_some_and(|q| {
                        q.atom_has_query(AtomId::new(i))
                            && source_query_has_single_h(q.atom_predicate(AtomId::new(i)))
                    }),
                    closures[i].len(),
                    incident.iter().any(|nb| {
                        stereo::bond_order_as_double(topology.bonds[nb.bond.index()].order()) > 1.0
                    }),
                )
            {
                count = count.wrapping_add(1);
            }
            swaps[i] = count % 2 != 0;
            permutations[i] = perm;
        }
    }
    let mut atom_visits = vec![0usize; n];
    let mut bond_visits = vec![0usize; m];
    let mut pos = 1u32;
    for item in stack {
        match *item {
            MolStackElem::Atom(i) => {
                let len = atom_visits.len();
                *atom_visits
                    .get_mut(i)
                    .ok_or(SmilesParseError::TraversalStateIndex {
                        state: "atomVisitOrders",
                        index: i,
                        count: len,
                    })? = pos as usize;
            }
            MolStackElem::Bond { bond, .. } => {
                let i = bond.index();
                let len = bond_visits.len();
                *bond_visits
                    .get_mut(i)
                    .ok_or(SmilesParseError::TraversalStateIndex {
                        state: "bondVisitOrders",
                        index: i,
                        count: len,
                    })? = pos as usize;
                if matches!(
                    topology.bonds[i].direction(),
                    BondDirection::EndDownRight | BondDirection::EndUpRight
                ) {
                    topology.bonds[i].set_direction(BondDirection::None);
                }
            }
            _ => {}
        }
        pos = pos.wrapping_add(1);
    }
    let mut bond_counts = vec![0i8; m];
    let mut atom_counts = vec![0i8; n];
    direction::canonicalize_double_bonds_for_writer(
        topology,
        &bond_visits,
        &atom_visits,
        opened,
        &mut bond_counts,
        &mut atom_counts,
        stack,
    );
    if isomeric_smiles {
        let mut adjusted = vec![false; n];
        for item in stack {
            let MolStackElem::Atom(i) = *item else {
                continue;
            };
            if topology.atoms[i].chiral_tag() == ChiralTag::Unspecified
                || topology.atoms[i].prop("_brokenChirality").is_some()
            {
                continue;
            }
            if topology.atoms[i].prop("_ringStereoAtoms").is_some() {
                if !adjusted[i] {
                    topology.atoms[i].set_chiral_tag(ChiralTag::TetrahedralCcw);
                    adjusted[i] = true;
                }
                // Borrow the actual INT_VECT once after changing the source tag.
                // Disjoint slices expose the center and the other real atoms: no
                // vector clone, repeated property lookup, or prevalidation scan.
                let (before, center_and_after) = topology.atoms.split_at_mut(i);
                let (center, after) = center_and_after.split_first_mut().unwrap();
                let relations = center.prop("_ringStereoAtoms").unwrap().as_int_vector()?;
                for &rel in relations {
                    if rel == 0 || rel == i32::MIN {
                        return Err(SmilesParseError::WriterCanonicalInvariant(
                            "invalid source ring stereo neighbor index",
                        ));
                    }
                    let j = (rel.abs() - 1) as usize;
                    if j >= n {
                        return Err(SmilesParseError::TraversalStateIndex {
                            state: "ringStereoChemAdjusted",
                            index: j,
                            count: n,
                        });
                    }
                    if !adjusted[j] && atom_visits[j] > atom_visits[i] {
                        // This strict visit-order test excludes the center itself.
                        let neighbor = if j < i {
                            &mut before[j]
                        } else {
                            &mut after[j - i - 1]
                        };
                        neighbor.set_chiral_tag(center.chiral_tag());
                        if rel < 0 {
                            cosmolkit_core::invert_atom_chirality(neighbor)
                                .map_err(SmilesParseError::WriterStereoOrder)?;
                        }
                        if swaps[i] != swaps[j] {
                            cosmolkit_core::invert_atom_chirality(neighbor)
                                .map_err(SmilesParseError::WriterStereoOrder)?;
                        }
                        adjusted[j] = true;
                    }
                }
            } else {
                // The sparse map represents the actual source transient Any<size_t>
                // property set by canonicalizeEnhancedStereo. Otherwise the raw
                // source property getter is checked without narrowing or fallback.
                let group = if let Some(&g) = stereo_group_references.get(&i) {
                    Some(g as u64)
                } else {
                    topology.atoms[i]
                        .prop("_stereoGroup")
                        .map(cosmolkit_core::property_value_to_ulong)
                        .transpose()
                        .map_err(SmilesParseError::WriterULong)?
                };
                if let Some(g) = group.filter(|&g| g < (topology.stereo_groups.len() as u64)) {
                    let flip = topology.atoms[i].chiral_tag() == ChiralTag::TetrahedralCw;
                    if flip {
                        cosmolkit_core::invert_atom_chirality(&mut topology.atoms[i])
                            .map_err(SmilesParseError::WriterStereoOrder)?;
                    }
                    if flip || swaps[i] {
                        for k in 0..topology.stereo_groups[g as usize].atoms().len() {
                            let j = topology.stereo_groups[g as usize].atoms()[k].index();
                            if j == i {
                                continue;
                            }
                            let atom = topology.atoms.get_mut(j).ok_or(
                                SmilesParseError::TraversalStateIndex {
                                    state: "stereoGroupAtom",
                                    index: j,
                                    count: n,
                                },
                            )?;
                            cosmolkit_core::invert_atom_chirality(atom)
                                .map_err(SmilesParseError::WriterStereoOrder)?;
                        }
                    }
                } else if matches!(
                    topology.atoms[i].chiral_tag(),
                    ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
                ) {
                    if swaps[i] {
                        cosmolkit_core::invert_atom_chirality(&mut topology.atoms[i])
                            .map_err(SmilesParseError::WriterStereoOrder)?;
                    }
                } else if permutations[i] != 0 {
                    topology.atoms[i]
                        .set_prop("_chiralPermutation", PropertyValue::Int(permutations[i]))?;
                }
            }
        }
    }
    direction::remove_unwanted_bond_dir_specs_for_writer(
        topology,
        stack,
        &mut bond_counts,
        &mut atom_counts,
        &bond_visits,
    );
    direction::remove_redundant_bond_dir_specs_for_writer(
        topology,
        stack,
        &mut bond_counts,
        &mut atom_counts,
    );
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn compute_chiral_adjustments(
    topology: &mut TopologyBlock,
    valence: &ValenceAssignment,
    rings: Option<&cosmolkit_core::RingInfo>,
    isomeric_smiles: bool,
    _start_atom: usize,
    closures: &[Vec<BondId>],
    orders: &[Vec<BondId>],
    stack: &[MolStackElem],
    groups: &BTreeMap<usize, usize>,
    atoms_in_play: Option<&[bool]>,
    bonds_in_play: Option<&[bool]>,
    opened: &[bool],
) -> Result<Vec<ChiralAdjustment>, SmilesParseError> {
    canonicalize_source_stack(
        topology,
        valence,
        rings,
        stack,
        closures,
        orders,
        opened,
        atoms_in_play,
        bonds_in_play,
        isomeric_smiles,
        true,
        groups,
        None,
    )?;
    // Emission reads actual mutated fields; no parallel chiral algorithm.
    Ok(vec![ChiralAdjustment::default(); topology.atoms.len()])
}

fn parse_ring_stereo_atoms(
    encoded: &cosmolkit_model::PropertyValue,
    atom_count: usize,
) -> Result<Vec<(bool, usize)>, SmilesParseError> {
    // BEGIN RDKIT CPP TYPE RDGeneral/types.h INT_VECT
    // RDKit❗❌: typedef std::vector<int> INT_VECT;
    // END RDKIT CPP TYPE RDGeneral/types.h INT_VECT
    // BEGIN RDKIT CPP FUNCTION Canon::canonicalizeFragment ring-relative references
    // RDKit❗❌:           const INT_VECT &ringStereoAtoms = atom->getProp<INT_VECT>(
    // RDKit❗❌:               common_properties::_ringStereoAtoms);
    // RDKit❗❌:           for (auto nbrV : ringStereoAtoms) {
    // RDKit❗❌:             int nbrIdx = abs(nbrV) - 1;
    // RDKit❗❌:               if (nbrV < 0) {
    // RDKit❗❌:                 mol.getAtomWithIdx(nbrIdx)->invertChirality();
    // END RDKIT CPP FUNCTION Canon::canonicalizeFragment ring-relative references
    // Source getProp<INT_VECT> is a strict tag cast, not string projection.
    // Preserve signed entries and duplicates. Borrow the stored vector and
    // scan once; no text parsing or compatibility encoding remains.
    let values = encoded.as_int_vector().map_err(|_| {
        SmilesParseError::WriterStereo(
            "`_ringStereoAtoms` is not a signed source INT_VECT value (bad_any_cast)".into(),
        )
    })?;
    let mut result = Vec::with_capacity(values.len());
    for &value in values {
        if value == 0 {
            return Err(SmilesParseError::WriterStereo(
                "`_ringStereoAtoms` cannot contain zero".into(),
            ));
        }
        let index = usize::try_from(value.unsigned_abs() - 1).map_err(|_| {
            SmilesParseError::WriterStereo("`_ringStereoAtoms` index is out of range".into())
        })?;
        if index >= atom_count {
            return Err(SmilesParseError::WriterStereo(
                "`_ringStereoAtoms` index is out of range".into(),
            ));
        }
        result.push((value > 0, index));
    }
    Ok(result)
}

fn write_mol_stack(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    stack: &[MolStackElem],
    chiral_adjustments: &[ChiralAdjustment],
    params: &SmilesWriteParams,
    atom_symbols: Option<&[String]>,
    bond_symbols: Option<&[String]>,
) -> Result<PropertyText, SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION FragmentSmilesConstruct MolStack emission section
    // RDKit✔️❌:   for (auto &mSE : molStack) {
    // RDKit✔️❌:     switch (mSE.type) {
    // RDKit✔️❌:       case Canon::MOL_STACK_ATOM:
    // RDKit✔️❌:         for (auto rclosure : ringClosuresToErase) {
    // RDKit✔️❌:           ringClosureMap.erase(rclosure);
    // RDKit✔️❌:         }
    // RDKit✔️❌:         ringClosuresToErase.clear();
    // RDKit✔️❌:         res << GetAtomSmiles(mSE.obj.atom, params);
    // RDKit✔️❌:         atomOrdering.push_back(mSE.obj.atom->getIdx());
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       case Canon::MOL_STACK_BOND:
    // RDKit✔️❌:         bond = mSE.obj.bond;
    // RDKit✔️❌:         res << GetBondSmiles(bond, params, mSE.number);
    // RDKit✔️❌:         bondOrdering.push_back(bond->getIdx());
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       case Canon::MOL_STACK_RING:
    // RDKit✔️❌:         ringIdx = mSE.number;
    // RDKit✔️❌:         if (ringClosureMap.count(ringIdx)) {
    // RDKit✔️❌:           closureVal = ringClosureMap[ringIdx];
    // RDKit✔️❌:           ringClosuresToErase.push_back(ringIdx);
    // RDKit✔️❌:         } else {
    // RDKit✔️❌:           closureVal = 1;
    // RDKit✔️❌:           bool done = false;
    // RDKit✔️❌:           while (!done) {
    // RDKit✔️❌:             std::map<int, int>::iterator mapIt;
    // RDKit✔️❌:             for (mapIt = ringClosureMap.begin();
    // RDKit✔️❌:                  mapIt != ringClosureMap.end(); ++mapIt) {
    // RDKit✔️❌:               if (mapIt->second == closureVal) {
    // RDKit✔️❌:                 break;
    // RDKit✔️❌:               }
    // RDKit✔️❌:             }
    // RDKit✔️❌:             if (mapIt == ringClosureMap.end()) {
    // RDKit✔️❌:               done = true;
    // RDKit✔️❌:             } else {
    // RDKit✔️❌:               closureVal += 1;
    // RDKit✔️❌:             }
    // RDKit✔️❌:           }
    // RDKit✔️❌:           ringClosureMap[ringIdx] = closureVal;
    // RDKit✔️❌:         }
    // RDKit✔️❌:         if (closureVal < 10) {
    // RDKit✔️❌:           res << (char)(closureVal + '0');
    // RDKit✔️❌:         } else if (closureVal < 100) {
    // RDKit✔️❌:           res << '%' << closureVal;
    // RDKit✔️❌:         } else {
    // RDKit✔️❌:           res << "%(" << closureVal << ')';
    // RDKit✔️❌:         }
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       case Canon::MOL_STACK_BRANCH_OPEN:
    // RDKit✔️❌:         res << "(";
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       case Canon::MOL_STACK_BRANCH_CLOSE:
    // RDKit✔️❌:         res << ")";
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       default:
    // RDKit✔️❌:         break;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // END RDKIT CPP FUNCTION FragmentSmilesConstruct MolStack emission section
    // `write_smiles_output` collects atom and bond output-order rows in two
    // separate scans after this emission pass; RDKit records both in this loop.
    let mut output = PropertyText::new();
    let mut display_digits = BTreeMap::<i32, usize>::new();
    let mut closures_to_erase = Vec::new();
    for element in stack {
        match *element {
            MolStackElem::Atom(atom) => {
                for ring_id in closures_to_erase.drain(..) {
                    display_digits.remove(&ring_id);
                }
                if let Some(symbols) = atom_symbols {
                    output.extend_bytes((&symbols[atom]).as_ref());
                } else {
                    output.extend_bytes(
                        (&atom_text(topology, valence, atom, chiral_adjustments[atom], params)?)
                            .as_ref(),
                    );
                }
            }
            MolStackElem::Bond { bond, atom_to_left } => {
                if let Some(symbols) = bond_symbols {
                    output.extend_bytes((&symbols[bond.index()]).as_ref());
                } else {
                    output.extend_bytes(
                        (&bond_text(
                            topology,
                            &topology.bonds[bond.index()],
                            atom_to_left,
                            params,
                        )?)
                            .as_ref(),
                    );
                }
            }
            MolStackElem::Ring(ring_id) => {
                let display_digit = if let Some(&digit) = display_digits.get(&ring_id) {
                    closures_to_erase.push(ring_id);
                    digit
                } else {
                    let digit = (1..)
                        .find(|candidate| !display_digits.values().any(|used| used == candidate))
                        .expect("positive ring-label space is not exhaustible");
                    display_digits.insert(ring_id, digit);
                    digit
                };
                write_ring_label(&mut output, display_digit);
            }
            MolStackElem::BranchOpen(_) => output.push_byte(b'('),
            MolStackElem::BranchClose(_) => output.push_byte(b')'),
        }
    }
    Ok(output)
}

fn atom_text(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    atom_index: usize,
    chiral_adjustment: ChiralAdjustment,
    params: &SmilesWriteParams,
) -> Result<PropertyText, SmilesParseError> {
    let atom = &topology.atoms[atom_index];
    // BEGIN RDKIT CPP FUNCTION atomNeedsBracket
    // RDKit❗✔️: bool atomNeedsBracket(const Atom *atom, const std::string &atString,
    // RDKit❗✔️:                       const SmilesWriteParams &params) {
    // RDKit❗✔️:   PRECONDITION(atom, "null atom");
    // RDKit❗✔️:   auto num = atom->getAtomicNum();
    // RDKit❗✔️:   if (!inOrganicSubset(num)) {
    // RDKit❗✔️:     return true;
    // RDKit❗✔️:   }
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
    // RDKit❗✔️:   if (atom->hasOwningMol()) {
    // RDKit❗✔️:     for (const auto bond : atom->getOwningMol().atomBonds(atom)) {
    // RDKit❗✔️:       auto oatom = bond->getOtherAtom(atom);
    // RDKit❗✔️:       if (QueryOps::isMetal(*oatom)) {
    // RDKit❗✔️:         return true;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return false;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION atomNeedsBracket
    // BEGIN RDKIT CPP FUNCTION GetAtomSmiles modeled non-stereo atom emission
    // RDKit❗✔️:   bool needsBracket = true;
    // RDKit❗✔️:   if (!hasCustomSymbol && !params.allHsExplicit) {
    // RDKit❗✔️:     needsBracket = atomNeedsBracket(atom, atString, params);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (needsBracket) {
    // RDKit❗✔️:     res += "[";
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (isotope && params.doIsomericSmiles) {
    // RDKit❗✔️:     res += std::to_string(isotope);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   res += symb;
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
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION GetAtomSmiles modeled non-stereo atom emission
    // BEGIN RDKIT CPP FUNCTION GetAtomSmiles custom symbol selection and aromatic spelling
    // RDKit❗✔️:   std::string symb;
    // RDKit❗✔️:   bool hasCustomSymbol =
    // RDKit❗✔️:       atom->getPropIfPresent(common_properties::smilesSymbol, symb);
    // RDKit❗✔️:   if (!hasCustomSymbol) {
    // RDKit❗✔️:     symb = PeriodicTable::getTable()->getElementSymbol(num);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // this was originally only done for the organic subset,
    // RDKit❗✔️:   // applying it to other atom-types is a fix for Issue 3152751:
    // RDKit❗✔️:   // Only accept for atom->getAtomicNum() in [5, 6, 7, 8, 14, 15, 16, 33, 34,
    // RDKit❗✔️:   // 52]
    // RDKit❗✔️:   if (!params.doKekule && atom->getIsAromatic() && symb[0] >= 'A' &&
    // RDKit❗✔️:       symb[0] <= 'Z') {
    // RDKit❗✔️:     switch (atom->getAtomicNum()) {
    // RDKit❗✔️:       case 5:
    // RDKit❗✔️:       case 6:
    // RDKit❗✔️:       case 7:
    // RDKit❗✔️:       case 8:
    // RDKit❗✔️:       case 14:
    // RDKit❗✔️:       case 15:
    // RDKit❗✔️:       case 16:
    // RDKit❗✔️:       case 33:
    // RDKit❗✔️:       case 34:
    // RDKit❗✔️:       case 52:
    // RDKit❗✔️:         symb[0] -= ('A' - 'a');
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION GetAtomSmiles custom symbol selection and aromatic spelling
    let custom_symbol = atom
        .prop("smilesSymbol")
        .map(cosmolkit_core::property_value_to_string)
        .transpose()
        .map_err(SmilesParseError::WriterProperty)?;
    // RDKit Dict's std::string overload uses rdvalue_tostring, including
    // scalar tags. This is not the INT_VECT cast used by _ringStereoAtoms.
    let raw_symbol = custom_symbol
        .as_ref()
        .map(PropertyText::as_bytes)
        .unwrap_or_else(|| atom.element().symbol().as_bytes());
    let mut symbol = raw_symbol.to_vec();
    if !params.kekule
        && atom.is_aromatic()
        && raw_symbol.first().is_some_and(u8::is_ascii_uppercase)
        && matches!(
            atom.atomic_number(),
            5 | 6 | 7 | 8 | 14 | 15 | 16 | 33 | 34 | 52
        )
    {
        symbol[0] = symbol[0].to_ascii_lowercase();
    }
    let symbol = PropertyText::from(symbol);
    // BEGIN RDKIT CPP FUNCTION GetAtomSmiles chirality selection
    // RDKit❗✔️:   if (params.doIsomericSmiles) {
    // RDKit❗✔️:     if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit❗✔️:         !atom->hasProp(common_properties::_brokenChirality)) {
    // RDKit❗✔️:       atString = getAtomChiralityInfo(atom);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION GetAtomSmiles chirality selection
    let mut chirality = String::new();
    if params.isomeric_smiles && atom.prop("_brokenChirality").is_none() {
        let base_chiral_tag = chiral_adjustment
            .chiral_tag_override
            .unwrap_or(atom.chiral_tag());
        let chiral_tag = if chiral_adjustment.invert_tetrahedral {
            stereo::invert_tetrahedral_tag(base_chiral_tag)
        } else {
            base_chiral_tag
        };
        chirality = match chiral_tag {
            ChiralTag::TetrahedralCw => "@@".to_owned(),
            ChiralTag::TetrahedralCcw => "@".to_owned(),
            ChiralTag::SquarePlanar => "@SP".to_owned(),
            ChiralTag::TrigonalBipyramidal => "@TB".to_owned(),
            ChiralTag::Octahedral => "@OH".to_owned(),
            _ => String::new(),
        };
        if matches!(
            chiral_tag,
            ChiralTag::SquarePlanar | ChiralTag::TrigonalBipyramidal | ChiralTag::Octahedral
        ) {
            // RDKit❗✔️: std::string getAtomChiralityInfo(const Atom *atom) {
            // RDKit❗✔️:   auto allowNontet = Chirality::getAllowNontetrahedralChirality();
            // RDKit❗✔️:   std::string atString;
            // RDKit❗✔️:   switch (atom->getChiralTag()) {
            // RDKit❗✔️:     case Atom::CHI_TETRAHEDRAL_CW:
            // RDKit❗✔️:       atString = "@@";
            // RDKit❗✔️:       break;
            // RDKit❗✔️:     case Atom::CHI_TETRAHEDRAL_CCW:
            // RDKit❗✔️:       atString = "@";
            // RDKit❗✔️:       break;
            // RDKit❗✔️:     default:
            // RDKit❗✔️:       break;
            // RDKit❗✔️:   }
            // RDKit❗✔️:   if (atString.empty() && allowNontet) {
            // RDKit❗✔️:     switch (atom->getChiralTag()) {
            // RDKit❗✔️:       case Atom::CHI_SQUAREPLANAR:
            // RDKit❗✔️:         atString = "@SP";
            // RDKit❗✔️:         break;
            // RDKit❗✔️:       case Atom::CHI_TRIGONALBIPYRAMIDAL:
            // RDKit❗✔️:         atString = "@TB";
            // RDKit❗✔️:         break;
            // RDKit❗✔️:       case Atom::CHI_OCTAHEDRAL:
            // RDKit❗✔️:         atString = "@OH";
            // RDKit❗✔️:         break;
            // RDKit❗✔️:       default:
            // RDKit❗✔️:         break;
            // RDKit❗✔️:     }
            // RDKit❗✔️:     if (!atString.empty()) {
            // RDKit❗✔️:       // we added info about non-tetrahedral stereo, so check whether or not
            // RDKit❗✔️:       // we need to also add permutation info
            // RDKit❗✔️:       int permutation = 0;
            // RDKit❗✔️:       if (atom->getChiralTag() > Atom::ChiralType::CHI_OTHER &&
            // RDKit❗✔️:           atom->getPropIfPresent(common_properties::_chiralPermutation,
            // RDKit❗✔️:                                  permutation) &&
            // RDKit❗✔️:           !SmilesParseOps::checkChiralPermutation(atom->getChiralTag(),
            // RDKit❗✔️:                                                   permutation)) {
            // RDKit❗✔️:         throw ValueErrorException("bad chirality spec");
            // RDKit❗✔️:       } else if (permutation) {
            // RDKit❗✔️:         atString += std::to_string(permutation);
            // RDKit❗✔️:       }
            // RDKit❗✔️:     }
            // RDKit❗✔️:   }
            // RDKit❗✔️:   return atString;
            // RDKit❗✔️: }
            // RDKit❗✔️: bool checkChiralPermutation(int chiralTag, int permutation) {
            // RDKit❗✔️:   if (chiralTag > RDKit::Atom::ChiralType::CHI_OTHER &&
            // RDKit❗✔️:       permutationLimits.find(chiralTag) != permutationLimits.end() &&
            // RDKit❗✔️:       (permutation < 0 || permutation > permutationLimits.at(chiralTag))) {
            // RDKit❗✔️:     return false;
            // RDKit❗✔️:   }
            // RDKit❗✔️:   return true;
            // RDKit❗✔️: }
            // Read the actual dictionary committed by canonicalizeFragment.
            // Explicit detached source facts use the same checked int getter.
            let permutation = if let Some(value) = chiral_adjustment.nontetrahedral_permutation {
                cosmolkit_core::property_value_to_int(&PropertyValue::UInt(value))
                    .map_err(SmilesParseError::WriterInt)?
            } else {
                source_atom_permutation(atom)?
            };
            let limit = match chiral_tag {
                ChiralTag::SquarePlanar => 3,
                ChiralTag::TrigonalBipyramidal => 20,
                ChiralTag::Octahedral => 30,
                _ => unreachable!(),
            };
            if permutation < 0 || permutation > limit {
                return Err(SmilesParseError::WriterStereo(format!(
                    "invalid {} permutation {permutation}; maximum is {limit}",
                    chiral_tag.rdkit_name()
                )));
            }
            if permutation != 0 {
                chirality.push_str(&permutation.to_string());
            }
        }
    }
    let needs_bracket = if custom_symbol.is_some() || params.all_hydrogens_explicit {
        true
    } else {
        if !in_organic_subset(i32::from(atom.atomic_number())) {
            true
        } else if atom.formal_charge() != 0 {
            true
        } else if params.isomeric_smiles && (atom.isotope().is_some() || !chirality.is_empty()) {
            true
        } else if atom.atom_map().is_some() {
            true
        } else {
            let default_valences = cosmolkit_core::rdkit_valence_list(atom.atomic_number())
                .map_err(|error| SmilesParseError::WriterValence(error.to_string()))?
                .ok_or_else(|| {
                    SmilesParseError::WriterValence(format!(
                        "RDKit valence list is unavailable for atomic number {}",
                        atom.atomic_number()
                    ))
                })?;
            let default_valence = default_valences.first().copied().ok_or_else(|| {
                SmilesParseError::WriterValence(format!(
                    "RDKit valence list is empty for atomic number {}",
                    atom.atomic_number()
                ))
            })?;
            let total_valence = valence.explicit_valence[atom_index]
                + valence.implicit_hydrogens[atom_index].max(0);
            let total_num_hydrogens = usize::from(atom.explicit_hydrogens())
                + usize::try_from(valence.implicit_hydrogens[atom_index].max(0))
                    .unwrap_or(usize::MAX);
            let nonstandard_valence = if atom.radical_electrons() != 0
                || (matches!(atom.atomic_number(), 7 | 15)
                    && atom.is_aromatic()
                    && atom.explicit_hydrogens() > 0)
            {
                true
            } else {
                total_valence != default_valence && total_num_hydrogens > 0
            };
            nonstandard_valence
                || topology
                    .adjacency
                    .neighbors_of(atom_index)
                    .iter()
                    .any(|neighbor| {
                        rdkit_query_ops_is_metal(
                            topology.atoms[neighbor.atom_index].atomic_number(),
                        )
                    })
        }
    };
    if !needs_bracket {
        return append_supplemental_label(atom, symbol);
    }

    let mut output = PropertyText::from("[");
    if params.isomeric_smiles {
        if let Some(isotope) = atom.isotope() {
            output.extend_bytes((&isotope.to_string()).as_ref());
        }
    }
    output.extend_bytes((&symbol).as_ref());
    output.extend_bytes((&chirality).as_ref());
    let total_num_hydrogens = usize::from(atom.explicit_hydrogens())
        + usize::try_from(valence.implicit_hydrogens[atom_index].max(0)).unwrap_or(usize::MAX);
    if total_num_hydrogens > 0 {
        output.push_byte(b'H');
        if total_num_hydrogens > 1 {
            output.extend_bytes((&total_num_hydrogens.to_string()).as_ref());
        }
    }
    match atom.formal_charge() {
        0 => {}
        1 => output.push_byte(b'+'),
        -1 => output.push_byte(b'-'),
        charge if charge > 1 => {
            output.push_byte(b'+');
            output.extend_bytes((&charge.to_string()).as_ref());
        }
        charge => output.extend_bytes((&charge.to_string()).as_ref()),
    }
    if let Some(atom_map) = atom.atom_map() {
        output.push_byte(b':');
        output.extend_bytes((&atom_map.to_string()).as_ref());
    }
    output.push_byte(b']');
    append_supplemental_label(atom, output)
}

fn append_supplemental_label(
    atom: &Atom,
    mut text: PropertyText,
) -> Result<PropertyText, SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION GetAtomSmiles supplemental label
    // RDKit❗✔️:   // If the atom has this property, the contained string will
    // RDKit❗✔️:   // be inserted directly in the SMILES:
    // RDKit❗✔️:   std::string label;
    // RDKit❗✔️:   if (atom->getPropIfPresent(common_properties::_supplementalSmilesLabel,
    // RDKit❗✔️:                              label)) {
    // RDKit❗✔️:     res += label;
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION GetAtomSmiles supplemental label
    if let Some(label) = atom.prop("_supplementalSmilesLabel") {
        text.extend_bytes(
            (&cosmolkit_core::property_value_to_string(label)
                .map_err(SmilesParseError::WriterProperty)?)
                .as_ref(),
        );
    }
    Ok(text)
}

fn bond_text(
    topology: &TopologyBlock,
    bond: &Bond,
    atom_to_left: usize,
    params: &SmilesWriteParams,
) -> Result<&'static str, SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION GetBondSmiles aromatic-context selection
    // RDKit❗✔️:   bool aromatic = false;
    // RDKit❗✔️:   if (!params.doKekule && (bond->getBondType() == Bond::SINGLE ||
    // RDKit❗✔️:                            bond->getBondType() == Bond::DOUBLE ||
    // RDKit❗✔️:                            bond->getBondType() == Bond::AROMATIC)) {
    // RDKit❗✔️:     if (bond->hasOwningMol()) {
    // RDKit❗✔️:       auto a1 = bond->getOwningMol().getAtomWithIdx(atomToLeftIdx);
    // RDKit❗✔️:       auto a2 = bond->getOwningMol().getAtomWithIdx(
    // RDKit❗✔️:           bond->getOtherAtomIdx(atomToLeftIdx));
    // RDKit❗✔️:       if ((a1->getIsAromatic() && a2->getIsAromatic()) &&
    // RDKit❗✔️:           (a1->getAtomicNum() || a2->getAtomicNum())) {
    // RDKit❗✔️:         aromatic = true;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       aromatic = false;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION GetBondSmiles aromatic-context selection
    // BEGIN RDKIT CPP FUNCTION GetBondSmiles direction selection
    // RDKit❗✔️:   Bond::BondDir dir = bond->getBondDir();
    // RDKit❗✔️:   switch (bond->getBondType()) {
    // RDKit❗✔️:     case Bond::SINGLE:
    // RDKit❗✔️:       if (dir != Bond::NONE && dir != Bond::UNKNOWN) {
    // RDKit❗✔️:         switch (dir) {
    // RDKit❗✔️:           case Bond::ENDDOWNRIGHT:
    // RDKit❗✔️:             if (params.allBondsExplicit || params.doIsomericSmiles) {
    // RDKit❗✔️:               res = "\\";
    // RDKit❗✔️:             }
    // RDKit❗✔️:             break;
    // RDKit❗✔️:           case Bond::ENDUPRIGHT:
    // RDKit❗✔️:             if (params.allBondsExplicit || params.doIsomericSmiles) {
    // RDKit❗✔️:               res = "/";
    // RDKit❗✔️:             }
    // RDKit❗✔️:             break;
    // RDKit❗✔️:           default:
    // RDKit❗✔️:             if (params.allBondsExplicit) {
    // RDKit❗✔️:               res = "-";
    // RDKit❗✔️:             }
    // RDKit❗✔️:             break;
    // RDKit❗✔️:         }
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
    // RDKit❗✔️:       if (dir != Bond::NONE && dir != Bond::UNKNOWN) {
    // RDKit❗✔️:         switch (dir) {
    // RDKit❗✔️:           case Bond::ENDDOWNRIGHT:
    // RDKit❗✔️:             if (params.allBondsExplicit || params.doIsomericSmiles) {
    // RDKit❗✔️:               res = "\\";
    // RDKit❗✔️:             }
    // RDKit❗✔️:             break;
    // RDKit❗✔️:           case Bond::ENDUPRIGHT:
    // RDKit❗✔️:             if (params.allBondsExplicit || params.doIsomericSmiles) {
    // RDKit❗✔️:               res = "/";
    // RDKit❗✔️:             }
    // RDKit❗✔️:             break;
    // RDKit❗✔️:           default:
    // RDKit❗✔️:             if (params.allBondsExplicit || !aromatic) {
    // RDKit❗✔️:               res = ":";
    // RDKit❗✔️:             }
    // RDKit❗✔️:             break;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       } else if (params.allBondsExplicit || !aromatic) {
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
    // END RDKIT CPP FUNCTION GetBondSmiles direction selection
    let other = other_atom(bond, atom_to_left)?;
    let aromatic_context = !params.kekule
        && matches!(
            bond.order(),
            BondOrder::Single | BondOrder::Double | BondOrder::Aromatic
        )
        && topology.atoms[atom_to_left].is_aromatic()
        && topology.atoms[other].is_aromatic()
        && (topology.atoms[atom_to_left].atomic_number() != 0
            || topology.atoms[other].atomic_number() != 0);
    let direction = bond.direction();
    let direction_is_specified = !matches!(direction, BondDirection::None | BondDirection::Unknown);
    match bond.order() {
        BondOrder::Single => {
            if direction_is_specified {
                Ok(match direction {
                    BondDirection::EndDownRight
                        if params.all_bonds_explicit || params.isomeric_smiles =>
                    {
                        "\\"
                    }
                    BondDirection::EndUpRight
                        if params.all_bonds_explicit || params.isomeric_smiles =>
                    {
                        "/"
                    }
                    _ if params.all_bonds_explicit => "-",
                    _ => "",
                })
            } else if params.all_bonds_explicit || (aromatic_context && !bond.is_aromatic()) {
                Ok("-")
            } else {
                Ok("")
            }
        }
        BondOrder::Double => Ok(
            if !aromatic_context || !bond.is_aromatic() || params.all_bonds_explicit {
                "="
            } else {
                ""
            },
        ),
        BondOrder::Triple => Ok("#"),
        BondOrder::Quadruple => Ok("$"),
        BondOrder::Aromatic => {
            if direction_is_specified {
                Ok(match direction {
                    BondDirection::EndDownRight
                        if params.all_bonds_explicit || params.isomeric_smiles =>
                    {
                        "\\"
                    }
                    BondDirection::EndUpRight
                        if params.all_bonds_explicit || params.isomeric_smiles =>
                    {
                        "/"
                    }
                    _ if params.all_bonds_explicit || !aromatic_context => ":",
                    _ => "",
                })
            } else {
                Ok(if params.all_bonds_explicit || !aromatic_context {
                    ":"
                } else {
                    ""
                })
            }
        }
        BondOrder::Dative if bond.begin().index() == atom_to_left => Ok("->"),
        BondOrder::Dative => Ok("<-"),
        _ => Ok("~"),
    }
}

fn other_atom(bond: &Bond, atom: usize) -> Result<usize, SmilesParseError> {
    if bond.begin().index() == atom {
        Ok(bond.end().index())
    } else if bond.end().index() == atom {
        Ok(bond.begin().index())
    } else {
        Err(SmilesParseError::Model(format!(
            "bond {} is not incident to atom {atom}",
            bond.id().index()
        )))
    }
}

/// Exact source organic-subset membership over the native signed input domain.
#[doc(hidden)]
pub fn in_organic_subset(atomic_number: i32) -> bool {
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
    // Behavior: the source array is ascending until its -1 sentinel. Negative
    // input stops at index0 and returns false; nonnegative input returns true
    // exactly for the eleven entries before the sentinel, including dummy0.
    // Direct membership in those exact constants reproduces every signed int
    // input without narrowing, saturation, value guesses, errors, or fallback.
    // Complexity: constant bounded membership and no allocations/clones/maps;
    // the source scans at most twelve constant entries. Both have O(1) cost,
    // and this direct existing-owner membership has no added hot-path buffer.
    matches!(
        atomic_number,
        0 | 5 | 6 | 7 | 8 | 9 | 15 | 16 | 17 | 35 | 53
    )
}

fn rdkit_query_ops_is_metal(atomic_number: u8) -> bool {
    // BEGIN RDKIT CPP FUNCTION QueryOps::makeMAtomQuery / QueryOps::isMetal
    // RDKit✔️✔️: // !#0!#1!#2!#5!#6!#7!#8!#9!#10!#14!#15!#16!#17!#18!#33!#34!#35!#36!#52!#53!#54!#85!#86
    // RDKit✔️✔️: bool isMetal(const Atom &atom) {
    // RDKit✔️✔️:   static const std::unique_ptr<ATOM_OR_QUERY> q(makeMAtomQuery());
    // RDKit✔️✔️:   return q->Match(&atom);
    // RDKit✔️✔️: }
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

fn write_ring_label(output: &mut PropertyText, label: usize) {
    // RDKit✔️❌:         if (closureVal < 10) {
    // RDKit✔️❌:           res << (char)(closureVal + '0');
    // RDKit✔️❌:         } else if (closureVal < 100) {
    // RDKit✔️❌:           res << '%' << closureVal;
    // RDKit✔️❌:         } else {
    // RDKit✔️❌:           res << "%(" << closureVal << ')';
    // RDKit✔️❌:         }
    // The multi-digit branches allocate a temporary with `to_string`; RDKit
    // streams the integer into its existing output buffer.
    if label < 10 {
        output.push_byte(b'0' + label as u8);
    } else if label < 100 {
        output.push_byte(b'%');
        output.extend_bytes((&label.to_string()).as_ref());
    } else {
        output.extend_bytes(("%(").as_ref());
        output.extend_bytes((&label.to_string()).as_ref());
        output.push_byte(b')');
    }
}

fn rdkit_bond_type_code(order: BondOrder) -> i64 {
    match order {
        BondOrder::Unspecified => 0,
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

fn boost_hash_range(value: impl AsRef<[u8]>) -> u32 {
    // BEGIN RDKIT CPP FUNCTION Canon::dfsFindCycles custom bond-symbol hash
    // RDKit❗✔️: std::uint32_t hsh = gboost::hash_range(symb.begin(), symb.end());
    // END RDKIT CPP FUNCTION Canon::dfsFindCycles custom bond-symbol hash
    // BEGIN RDKIT CPP FUNCTION RDGeneral/hash/hash.hpp::hash_range
    // RDKit❗✔️: std::hash_result_t seed = 0;
    // RDKit❗✔️: for (; first != last; ++first) {
    // RDKit❗✔️:   hash_combine(seed, *first);
    // RDKit❗✔️: }
    // RDKit❗✔️: seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    // RDKit❗✔️: inline std::hash_result_t hash_value(char v) {
    // RDKit❗✔️:   return static_cast<std::hash_result_t>(v);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION RDGeneral/hash/hash.hpp::hash_range
    // RDKit pins a 32-bit hash result. Iterating UTF-8 bytes and explicitly
    // sign-extending each platform char matches its Linux C++ char hashing;
    // this is O(L) time, O(1) extra space, and does not allocate.
    let mut seed = 0_u32;
    for byte in value.as_ref().iter().copied() {
        let hashed_char = (byte as i8 as i32) as u32;
        seed ^= hashed_char
            .wrapping_add(0x9e37_79b9)
            .wrapping_add(seed.wrapping_shl(6))
            .wrapping_add(seed >> 2);
    }
    seed
}

fn writer_ring_bonds(
    topology: &TopologyBlock,
    source_rings: Option<&cosmolkit_core::RingInfo>,
    component_count: usize,
) -> Vec<bool> {
    // RDKit✔️✔️:   if (nFrags == 1) {
    // RDKit✔️✔️:     res.emplace_back(new RWMol(mol));
    // BEGIN RDKIT CPP FUNCTION Canon::canonicalizeFragment ring-cache selection
    // RDKit✔️✔️:   if (!mol.getRingInfo()->isSymmSssr()) {
    // RDKit✔️✔️:     MolOps::findSSSR(mol);
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION Canon::canonicalizeFragment ring-cache selection
    // Behavior: the sole fragment is an RWMol copy retaining SymmSSSR. Changing
    // DATIVE to SINGLE does not reset that cache; freshly searching the changed
    // graph would incorrectly mark new metal cycles as existing source rings.
    // Complexity: O(E) membership projection, avoiding an O(V+E) bridge search.
    if let Some(rings) = source_rings.filter(|rings| component_count == 1 && rings.is_symm_sssr()) {
        return topology
            .bonds
            .iter()
            .map(|bond| rings.num_bond_rings(bond.id()) != 0)
            .collect();
    }
    find_ring_bonds(topology)
}

#[cfg(test)]
mod dative_ring_cache_tests {
    use super::*;

    #[test]
    fn converting_dative_bonds_preserves_symm_sssr_membership_on_the_writer_copy() {
        let mut record =
            crate::parse_smiles("C1CN->[Cu+2]1", &crate::SmilesParseParams::default()).unwrap();
        let source_rings = cosmolkit_core::symmetrized_sssr(
            &record.topology,
            &cosmolkit_core::RingSearchParams::default(),
        )
        .unwrap();
        assert!(source_rings.is_symm_sssr());
        assert!(source_rings.bond_rings().is_empty());
        for bond in &mut record.topology.bonds {
            if bond.order() == BondOrder::Dative {
                bond.set_order(BondOrder::Single);
            }
        }
        assert!(
            find_ring_bonds(&record.topology)
                .iter()
                .all(|in_ring| *in_ring)
        );
        assert!(
            writer_ring_bonds(&record.topology, Some(&source_rings), 1)
                .iter()
                .all(|in_ring| !*in_ring)
        );
        assert!(
            writer_ring_bonds(&record.topology, Some(&source_rings), 2)
                .iter()
                .all(|in_ring| *in_ring)
        );
    }
}

fn find_ring_bonds(topology: &TopologyBlock) -> Vec<bool> {
    // BEGIN RDKIT CPP FUNCTION findSSSR active bond selection
    // RDKit✔️🔝:   // Zero-order bonds are not candidates for rings, and dative bonds and
    // RDKit✔️🔝:   // hydrogen bonds may also be out
    // RDKit✔️🔝:   boost::dynamic_bitset<> activeBonds(nbnds);
    // RDKit✔️🔝:   activeBonds.set();
    // RDKit✔️🔝:   for (auto bond : mol.bonds()) {
    // RDKit✔️🔝:     if (auto bt = bond->getBondType();
    // RDKit✔️🔝:         bt == Bond::ZERO || (!includeDativeBonds && isDative(bt)) ||
    // RDKit✔️🔝:         (!includeHydrogenBonds && bt == Bond::HYDROGEN)) {
    // RDKit✔️🔝:       activeBonds[bond->getIdx()] = 0;
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:   }
    // BEGIN RDKIT CPP FUNCTION RingInfo::numBondRings
    // RDKit✔️🔝: unsigned int RingInfo::numBondRings(unsigned int idx) const {
    // RDKit✔️🔝:   PRECONDITION(df_init, "RingInfo not initialized");
    // RDKit✔️🔝:
    // RDKit✔️🔝:   if (idx < d_bondMembers.size()) {
    // RDKit✔️🔝:     return rdcast<unsigned int>(d_bondMembers[idx].size());
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:   return 0;
    // RDKit✔️🔝: }
    // END RDKIT CPP FUNCTION RingInfo::numBondRings
    // END RDKIT CPP FUNCTION findSSSR active bond selection
    // A cycle basis spans the source-selected active graph's cycle space, so
    // every bond in a cycle has a stored ring membership. Tarjan's bridge test
    // computes that same `numBondRings() != 0` predicate directly in O(V+E),
    // avoiding SSSR ring-vector construction without changing the result.
    let mut discovery = vec![usize::MAX; topology.atoms.len()];
    let mut low = vec![usize::MAX; topology.atoms.len()];
    let mut bridges = vec![false; topology.bonds.len()];
    let mut next_time = 0;

    fn visit(
        topology: &TopologyBlock,
        atom: usize,
        parent_bond: Option<BondId>,
        discovery: &mut [usize],
        low: &mut [usize],
        bridges: &mut [bool],
        next_time: &mut usize,
    ) {
        discovery[atom] = *next_time;
        low[atom] = *next_time;
        *next_time += 1;
        for neighbor in topology.adjacency.neighbors_of(atom) {
            if Some(neighbor.bond) == parent_bond {
                continue;
            }
            if !ring_perception_eligible(topology.bonds[neighbor.bond.index()].order()) {
                continue;
            }
            if discovery[neighbor.atom_index] == usize::MAX {
                visit(
                    topology,
                    neighbor.atom_index,
                    Some(neighbor.bond),
                    discovery,
                    low,
                    bridges,
                    next_time,
                );
                low[atom] = low[atom].min(low[neighbor.atom_index]);
                if low[neighbor.atom_index] > discovery[atom] {
                    bridges[neighbor.bond.index()] = true;
                }
            } else {
                low[atom] = low[atom].min(discovery[neighbor.atom_index]);
            }
        }
    }

    for atom in 0..topology.atoms.len() {
        if discovery[atom] == usize::MAX {
            visit(
                topology,
                atom,
                None,
                &mut discovery,
                &mut low,
                &mut bridges,
                &mut next_time,
            );
        }
    }
    bridges
        .into_iter()
        .enumerate()
        .map(|(index, is_bridge)| {
            ring_perception_eligible(topology.bonds[index].order()) && !is_bridge
        })
        .collect()
}

fn ring_perception_eligible(order: BondOrder) -> bool {
    !matches!(
        order,
        BondOrder::Zero
            | BondOrder::Dative
            | BondOrder::DativeOne
            | BondOrder::DativeLeft
            | BondOrder::DativeRight
            | BondOrder::Hydrogen
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{SmilesParseParams, parse_smiles};

    #[test]
    fn writer_kekulized_fragment_retains_existing_explicit_cache_and_recalculates_implicit() {
        let record = parse_smiles("CC", &SmilesParseParams::default()).expect("fixed C-C graph");
        // Source Atom::calcImplicitValence sees cached explicit=3 and neutral
        // carbon default valence=4, so implicit=1. It must not reconstruct
        // explicit=1 from the current single bond and then produce implicit=3.
        let retained = ValenceAssignment {
            explicit_valence: vec![3, 3],
            implicit_hydrogens: vec![0, 0],
        };
        let assignment = kekulize_writer_fragment(
            &record.topology,
            &retained,
            &[AtomId::new(0), AtomId::new(1)],
        )
        .expect("source cache-preserving component input");
        let final_valence = assignment.final_valence.expect("selected rows retained");
        assert_eq!(final_valence.explicit_valence, [3, 3]);
        assert_eq!(final_valence.implicit_hydrogens, [1, 1]);
        assert_eq!(assignment.topology, record.topology);
        assert_eq!(retained.explicit_valence, [3, 3]);
        assert_eq!(retained.implicit_hydrogens, [0, 0]);
    }

    #[test]
    fn writer_kekulized_valence_maps_both_complete_cache_rows_without_recalculation() {
        let mut writer = ValenceAssignment {
            explicit_valence: vec![90, 91, 92, 93, 94],
            implicit_hydrogens: vec![80, 81, 82, 83, 84],
        };
        // Nonidentity/noncontiguous source order; these are opaque cache rows,
        // deliberately not reconstructed from topology or restricted to N/P.
        let fragment = ValenceAssignment {
            explicit_valence: vec![2, -1, 4],
            implicit_hydrogens: vec![1, 0, -1],
        };
        project_kekulized_writer_valence(
            &mut writer,
            Some(&fragment),
            &[AtomId::new(3), AtomId::new(0), AtomId::new(4)],
        )
        .expect("complete mapped cache");
        assert_eq!(writer.explicit_valence, [-1, 91, 92, 2, 4]);
        assert_eq!(writer.implicit_hydrogens, [0, 81, 82, 1, -1]);
    }

    #[test]
    fn writer_kekulized_valence_rejects_missing_or_invalid_cache_without_fallback() {
        let original = ValenceAssignment {
            explicit_valence: vec![3, 2],
            implicit_hydrogens: vec![0, 1],
        };
        for (cache, sources) in [
            (None, vec![AtomId::new(0)]),
            (
                Some(ValenceAssignment {
                    explicit_valence: vec![2],
                    implicit_hydrogens: vec![],
                }),
                vec![AtomId::new(0)],
            ),
            (
                Some(ValenceAssignment {
                    explicit_valence: vec![2],
                    implicit_hydrogens: vec![1],
                }),
                vec![AtomId::new(2)],
            ),
        ] {
            let mut writer = original.clone();
            assert!(matches!(
                project_kekulized_writer_valence(&mut writer, cache.as_ref(), &sources),
                Err(SmilesParseError::WriterValence(_))
            ));
            assert_eq!(writer, original);
        }
    }

    #[test]
    fn random_root_uses_source_modulo_for_one_component_and_advances_once() {
        let record = parse_smiles("CCCC", &SmilesParseParams::default()).expect("parse");
        let before = record.clone();
        let params = SmilesWriteParams {
            canonical: false,
            ..SmilesWriteParams::default()
        };

        // Pinned Boost1.85 x1=2_027_382 selects root 2. The three edges draw
        // x2..x4 in cycle DFS and x5..x7 in stack DFS; x8 is the continuation.
        let (output, continuation) = cosmolkit_core::with_rdkit_random_generator(42, |rng| {
            let output =
                write_smiles_output_with_random_stream(&record, &params, false, Some(&mut *rng))
                    .expect("random root write");
            (output, rng.next_u32())
        });

        assert_eq!(2_027_382_u32 % 4, 2);
        assert_eq!(output.atom_order[0], AtomId::new(2));
        assert_eq!(continuation, 1_538_354_858);
        assert_eq!(record, before);
    }

    #[test]
    fn random_root_draws_once_per_unrooted_component_in_source_order() {
        let record = parse_smiles("CCCC.CCC", &SmilesParseParams::default()).expect("parse");
        let before = record.clone();
        let params = SmilesWriteParams {
            canonical: false,
            ..SmilesWriteParams::default()
        };

        // x1=2_027_382 selects the first root; its three edges consume x2..x7
        // across cycle/stack DFS. x8=1_538_354_858 selects atom 6 in the second
        // component, whose edges consume x9..x12; x13 is the continuation.
        let (output, continuation) = cosmolkit_core::with_rdkit_random_generator(42, |rng| {
            let output =
                write_smiles_output_with_random_stream(&record, &params, false, Some(&mut *rng))
                    .expect("random root write");
            (output, rng.next_u32())
        });

        assert_eq!(2_027_382_u32 % 4, 2);
        assert_eq!(1_538_354_858_u32 % 3, 2);
        assert_eq!(output.atom_order[0], AtomId::new(2));
        assert_eq!(output.atom_order[4], AtomId::new(6));
        assert_eq!(continuation, 974_199_846);
        assert_eq!(record, before);
    }

    #[test]
    fn random_root_explicit_component_suppresses_only_its_draw() {
        let record = parse_smiles("CCCC.CCC", &SmilesParseParams::default()).expect("parse");
        let before = record.clone();
        let params = SmilesWriteParams {
            canonical: false,
            rooted_at_atom: Some(AtomId::new(1)),
            ..SmilesWriteParams::default()
        };

        // The first component keeps its explicit root without drawing; its
        // cycle/stack passes consume x1..x6. The second root uses x7=1_350_734_175
        // modulo 3, then its cycle/stack passes consume x8..x11; x12 continues.
        let (output, continuation) = cosmolkit_core::with_rdkit_random_generator(42, |rng| {
            let output =
                write_smiles_output_with_random_stream(&record, &params, false, Some(&mut *rng))
                    .expect("random root write");
            (output, rng.next_u32())
        });

        assert_eq!(1_350_734_175_u32 % 3, 0);
        assert_eq!(output.atom_order[0], AtomId::new(1));
        assert_eq!(output.atom_order[4], AtomId::new(4));
        assert_eq!(continuation, 1_151_860_813);
        assert_eq!(record, before);
    }

    #[test]
    fn random_cycle_uses_each_adjacency_draw_and_records_seed42_closure() {
        let record = parse_smiles("C1CC1", &SmilesParseParams::default()).expect("parse");
        let before = record.clone();
        assert_eq!(record.topology.bonds.len(), 3);
        assert_eq!(
            record
                .topology
                .adjacency
                .neighbors_of(0)
                .iter()
                .map(|neighbor| neighbor.bond.index())
                .collect::<Vec<_>>(),
            vec![0, 2]
        );

        // Pinned Boost1.85 seed42 draws are x1=2_027_382, x2=1_226_992_407,
        // x3=551_494_037, x4=961_371_815, x5=1_404_753_842. At atom 0,
        // x1 sorts bond 0 before bond 2; recursive eligible edges draw x3/x4
        // before the exact x5 continuation.
        let (colors, closures, continuation) =
            cosmolkit_core::with_rdkit_random_generator(42, |rng| {
                let mut colors = vec![AtomColor::White; 3];
                let mut closures = vec![Vec::new(); 3];
                dfs_find_cycles(
                    &record.topology,
                    0,
                    None,
                    &mut colors,
                    &[0, 1, 2],
                    &[true, true, true],
                    None,
                    None,
                    &mut closures,
                    Some(&mut *rng),
                );
                (colors, closures, rng.next_u32())
            });

        assert!(
            colors
                .iter()
                .all(|color| matches!(*color, AtomColor::Black))
        );
        assert_eq!(
            closures,
            vec![
                vec![BondId::new(2)],
                Vec::<BondId>::new(),
                vec![BondId::new(2)],
            ]
        );
        assert_eq!(continuation, 1_404_753_842);
        assert_eq!(record, before);
    }

    #[test]
    fn random_cycle_rank_reorders_adjacency_and_records_seed247_closure() {
        let record = parse_smiles("C1CC1", &SmilesParseParams::default()).expect("parse");
        let before = record.clone();
        assert_eq!(
            record
                .topology
                .adjacency
                .neighbors_of(0)
                .iter()
                .map(|neighbor| neighbor.bond.index())
                .collect::<Vec<_>>(),
            vec![0, 2]
        );

        // Independent recurrence values for seed247 are
        // [11_922_937, 6_474_531, 1_146_957_086, 489_594_999, 182_661_494].
        // The first two random ranks reverse the source adjacency order; the
        // later eligible edges close bond 0, then x5 is the continuation.
        let (colors, closures, continuation) =
            cosmolkit_core::with_rdkit_random_generator(247, |rng| {
                let mut colors = vec![AtomColor::White; 3];
                let mut closures = vec![Vec::new(); 3];
                dfs_find_cycles(
                    &record.topology,
                    0,
                    None,
                    &mut colors,
                    &[0, 1, 2],
                    &[false, false, false],
                    None,
                    None,
                    &mut closures,
                    Some(&mut *rng),
                );
                (colors, closures, rng.next_u32())
            });

        assert!(
            colors
                .iter()
                .all(|color| matches!(*color, AtomColor::Black))
        );
        assert_eq!(
            closures,
            vec![
                vec![BondId::new(0)],
                vec![BondId::new(0)],
                Vec::<BondId>::new(),
            ]
        );
        assert_eq!(continuation, 182_661_494);
        assert_eq!(record, before);
    }

    #[test]
    fn random_cycle_skips_masked_edges_without_draws() {
        let record = parse_smiles("C1CC1", &SmilesParseParams::default()).expect("parse");
        let before = record.clone();

        // With bond 1 out of play, only root bonds 0 and 2 draw x1/x2. The
        // incoming edges and the masked bond draw nothing; x3 is unchanged.
        let (colors, closures, continuation) =
            cosmolkit_core::with_rdkit_random_generator(42, |rng| {
                let mut colors = vec![AtomColor::White; 3];
                let mut closures = vec![Vec::new(); 3];
                dfs_find_cycles(
                    &record.topology,
                    0,
                    None,
                    &mut colors,
                    &[0, 1, 2],
                    &[true, true, true],
                    Some(&[true, false, true]),
                    None,
                    &mut closures,
                    Some(&mut *rng),
                );
                (colors, closures, rng.next_u32())
            });

        assert!(
            colors
                .iter()
                .all(|color| matches!(*color, AtomColor::Black))
        );
        assert_eq!(closures, vec![Vec::<BondId>::new(); 3]);
        assert_eq!(continuation, 551_494_037);
        assert_eq!(record, before);
    }

    fn string_property(value: Option<&PropertyValue>) -> Option<&str> {
        match value {
            Some(PropertyValue::String(value)) => Some(
                std::str::from_utf8(value.as_bytes()).expect("original fixture String is UTF-8"),
            ),
            _ => None,
        }
    }

    fn roundtrip(input: &str) -> String {
        let record = parse_smiles(input, &SmilesParseParams::default()).expect("parse");
        String::from_utf8(
            write_smiles_with_params(
                &record,
                &SmilesWriteParams {
                    canonical: false,
                    ..Default::default()
                },
            )
            .expect("write")
            .into_bytes(),
        )
        .expect("original fixture SMILES text is UTF-8")
    }

    #[test]
    fn source_order_writer_places_all_but_last_child_in_branches() {
        assert_eq!(roundtrip("C(O)N"), "C(O)N");
        assert_eq!(roundtrip("OC.C"), "OC.C");
    }

    #[test]
    fn writer_serializes_aliphatic_and_aromatic_cycles() {
        assert_eq!(roundtrip("C1CCCCC1"), "C1CCCCC1");
        assert_eq!(roundtrip("c1ccccc1"), "c1ccccc1");
        assert_eq!(roundtrip("c1ccccc-1"), "c1ccccc-1");
    }

    #[test]
    fn writer_matches_source_order_for_fused_bridged_and_spiro_cycles() {
        for smiles in [
            "C1CCC2CCCCC2C1",
            "C12(CCCCC1)CCCCC2",
            "C1C2C3C1C2C3",
            "C1CC2CCC1C2",
            "C1(C2CC2)CC1",
            "C1CC2(CC1)CCC2",
        ] {
            assert_eq!(roundtrip(smiles), smiles);
        }
    }

    #[test]
    fn writer_preserves_dative_ring_orientation() {
        assert_eq!(roundtrip("N1CC->1"), "N1CC->1");
        assert_eq!(roundtrip("N1CC<-1"), "N1CC<-1");
    }

    #[test]
    fn writer_emits_bracket_fields_without_loss() {
        assert_eq!(roundtrip("[13CH3+:7]C"), "[13CH3+:7]C");
        assert_eq!(roundtrip("[Na+]Cl"), "[Na+][Cl]");
    }

    #[test]
    fn writer_matches_rdkit_tetrahedral_traversal_inversions() {
        for (input, expected) in [
            ("[C@H](F)(Cl)Br", "F[C@@H](Cl)Br"),
            ("[C@@H](F)(Cl)Br", "F[C@H](Cl)Br"),
            ("F[C@H](Cl)Br", "F[C@H](Cl)Br"),
            ("Br[C@@H](Cl)F", "F[C@H](Cl)Br"),
            ("N[C@](F)(Cl)Br", "N[C@](F)(Cl)Br"),
            ("C[C@H](F)Cl", "C[C@H](F)Cl"),
            ("F[C@]1(Br)CCO1", "F[C@]1(Br)CCO1"),
            ("F[C@]1(CCO1)Br", "F[C@@]1(Br)CCO1"),
            ("[C@H](F)(Cl)Br.[C@@H](I)(N)O", "F[C@@H](Cl)Br.N[C@H](O)I"),
        ] {
            let record = parse_smiles(input, &Default::default()).unwrap();
            assert_eq!(write_smiles(&record).unwrap(), expected.into(), "{input}");
        }

        let bond_stereo = parse_smiles("F/C=C/F", &Default::default()).unwrap();
        assert_eq!(write_smiles(&bond_stereo).unwrap(), "F/C=C/F".into());
    }

    #[test]
    fn writer_matches_rdkit_nontetrahedral_permutation_reordering() {
        for (input, expected_noncanonical, expected_canonical) in [
            (
                "[Pt@SP1](F)(Cl)(Br)I",
                "[Pt@SP1]([F])([Cl])([Br])[I]",
                "[F][Pt@SP1]([Cl])([Br])[I]",
            ),
            (
                "[Pt@SP2](F)(Cl)(Br)I",
                "[Pt@SP2]([F])([Cl])([Br])[I]",
                "[F][Pt@SP2]([Cl])([Br])[I]",
            ),
            (
                "I[Pt@SP1](Br)(Cl)F",
                "[I][Pt@SP1]([Br])([Cl])[F]",
                "[F][Pt@SP1]([Cl])([Br])[I]",
            ),
            (
                "[Pt@SP1](F)(Cl)Br",
                "[Pt@SP1]([F])([Cl])[Br]",
                "[F][Pt@SP3]([Cl])[Br]",
            ),
            (
                "[P@TB1](F)(Cl)(Br)(I)N",
                "[P@TB1](F)(Cl)(Br)(I)N",
                "N[P@TB8](F)(Cl)(Br)I",
            ),
            (
                "[P@TB20](F)(Cl)(Br)(I)N",
                "[P@TB20](F)(Cl)(Br)(I)N",
                "N[P@TB3](F)(Cl)(Br)I",
            ),
            (
                "N[P@TB1](I)(Br)(Cl)F",
                "N[P@TB1](I)(Br)(Cl)F",
                "N[P@TB8](F)(Cl)(Br)I",
            ),
            (
                "[P@TB1](F)(Cl)(Br)I",
                "[P@TB1](F)(Cl)(Br)I",
                "F[P@TB9](Cl)(Br)I",
            ),
            (
                "[Co@OH1](F)(Cl)(Br)(I)(N)O",
                "[Co@OH1]([F])([Cl])([Br])([I])([NH2])[OH]",
                "[NH2][Co@OH9]([OH])([F])([Cl])([Br])[I]",
            ),
            (
                "[Co@OH30](F)(Cl)(Br)(I)(N)O",
                "[Co@OH30]([F])([Cl])([Br])([I])([NH2])[OH]",
                "[NH2][Co@OH15]([OH])([F])([Cl])([Br])[I]",
            ),
            (
                "O[Co@OH1](N)(I)(Br)(Cl)F",
                "[OH][Co@OH1]([NH2])([I])([Br])([Cl])[F]",
                "[NH2][Co@OH9]([OH])([F])([Cl])([Br])[I]",
            ),
            (
                "[Co@OH1](F)(Cl)(Br)(I)N",
                "[Co@OH1]([F])([Cl])([Br])([I])[NH2]",
                "[NH2][Co@OH30]([F])([Cl])([Br])[I]",
            ),
        ] {
            let record = parse_smiles(input, &Default::default()).unwrap();
            assert_eq!(
                write_smiles_with_params(
                    &record,
                    &SmilesWriteParams {
                        canonical: false,
                        ..Default::default()
                    }
                )
                .unwrap(),
                expected_noncanonical.into(),
                "non-canonical {input}"
            );
            assert_eq!(
                write_smiles(&record).unwrap(),
                expected_canonical.into(),
                "canonical {input}"
            );
        }
    }

    #[test]
    fn writer_canonicalizes_ring_relative_stereo_like_rdkit() {
        for (input, relation, expected) in [
            (
                "C1[C@H](F)CC[C@H](Cl)C1",
                (6_i32, 2_i32),
                "F[C@H]1CC[C@@H](Cl)CC1",
            ),
            (
                "C1[C@H](F)CC[C@@H](Cl)C1",
                (-6_i32, -2_i32),
                "F[C@H]1CC[C@H](Cl)CC1",
            ),
        ] {
            let mut record = parse_smiles(input, &Default::default()).unwrap();
            record.topology.atoms[1].set_prop("_ringStereoAtoms", vec![relation.0]);
            record.topology.atoms[5].set_prop("_ringStereoAtoms", vec![relation.1]);
            assert_eq!(write_smiles(&record).unwrap(), expected.into(), "{input}");
        }
    }

    #[test]
    fn writer_rejects_malformed_ring_relative_stereo_state() {
        let mut record = parse_smiles("C1[C@H](F)CC[C@H](Cl)C1", &Default::default()).unwrap();
        record.topology.atoms[1].set_prop("_ringStereoAtoms", "0");
        assert!(matches!(
            write_smiles(&record),
            Err(SmilesParseError::WriterStereo(_))
        ));
    }

    #[test]
    fn canonical_writer_matches_rdkit_double_bond_directions() {
        let mut mismatches = Vec::new();
        for (input, expected) in [
            ("F/C=C/F", "F/C=C/F"),
            ("F/C=C\\F", "F/C=C\\F"),
            ("Br\\C=C/F", "F/C=C\\Br"),
            ("F\\C=C\\Br", "F/C=C/Br"),
            ("F\\C=C(/Cl)\\Br", "F/C=C(\\Cl)Br"),
            ("F/C=C(\\Cl)/Br", "F/C=C(\\Cl)Br"),
            ("C/C=C/C", "C/C=C/C"),
            ("C/C=C\\C", "C/C=C\\C"),
            ("C/C=C/C=C\\C", "C/C=C\\C=C\\C"),
            ("Cl/C=C(/C=C/C)\\C=C\\Br", "C/C=C/C(=C/Cl)/C=C/Br"),
            ("C(\\C/C=C/Cl)=C/O", "O/C=C/C/C=C/Cl"),
            ("O=C\\C=C/F", "O=C/C=C\\F"),
            ("C(=O)\\C=C/Br", "O=C/C=C\\Br"),
            ("CC(=O)\\C=C/Br", "CC(=O)/C=C\\Br"),
            ("C/C=C(/C)C", "CC=C(C)C"),
            ("C/C=C(/C(F))C(Cl)", "C/C=C(/CF)CCl"),
            ("C/C=C(/C(F))C(F)", "CC=C(CF)CF"),
            ("C/C=C(/CO)CN", "C/C=C(\\CN)CO"),
            ("F/C=C(/F)Cl", "F/C=C(/F)Cl"),
            ("F/C(Cl)=C(/Br)I", "F/C(Cl)=C(/Br)I"),
            ("C/C=C(/C)\\C", "CC=C(C)C"),
            ("C1COC/C=C\\CCC1", "C1=C\\COCCCCC/1"),
            ("C1COC/C=C/CCC1", "C1=C/COCCCCC/1"),
            ("C1CC/C=C/C=C/CCC1", "C1=C/CCCCCC/C=C/1"),
            ("C/1=C/C=C/CCCCCC1", "C1=C\\CCCCCC/C=C/1"),
            ("C1COC/C=C/C=C/C1", "C1=C/CCCOC/C=C/1"),
            ("C1=C/OCC/C=C\\CC\\1", "C1=C\\CCO/C=C\\CC/1"),
            ("CO/C1=C/C=C\\C=C/C=N\\1", "COC1=C/C=C\\C=C/C=N\\1"),
            ("C1C/C=C/CCCCCCCC1", "C1=C/CCCCCCCCCC/1"),
            ("C1/C=C/C=C/CCCCCCCCC1", "C1=C/CCCCCCCCCC/C=C/1"),
            ("C1/C=C\\CCCCC1", "C1=C\\CCCCCC/1"),
            // This crate's parse boundary is deliberately unsanitized. RDKit
            // with sanitize=false/removeHs=false preserves these directions;
            // its later full sanitize pipeline removes them as non-potential
            // fused-ring stereo.
            ("C1/C=C\\C2CCCCC2C1", "C1=C\\C2CCCCC2CC/1"),
            ("C1C/C=C/CC2CCCCC12", "C1=C/CC2CCCCC2CC/1"),
            ("F/C=C/F.C/C=C\\C", "C/C=C\\C.F/C=C/F"),
        ] {
            let record = parse_smiles(input, &Default::default()).unwrap();
            let actual = write_smiles(&record).unwrap();
            if actual != expected.into() {
                mismatches.push((input, expected, actual));
            }
        }
        assert!(mismatches.is_empty(), "{mismatches:#?}");
    }

    #[test]
    fn detached_legacy_cip_ranks_match_rdkit_compute_atom_cip_ranks() {
        for (input, expected) in [
            ("CCO", vec![0, 1, 2]),
            ("Cl/C=C(/C=C/C)\\C=C\\Br", vec![7, 5, 4, 2, 1, 0, 3, 6, 8]),
            ("C1CCCCC1", vec![0, 0, 0, 0, 0, 0]),
            ("c1ccccc1", vec![0, 0, 0, 0, 0, 0]),
            ("CC(C)O", vec![0, 1, 0, 2]),
            ("N#CC(=O)O", vec![2, 0, 1, 4, 3]),
            ("F[C@H](Cl)Br", vec![1, 0, 2, 3]),
            ("OP(=O)(O)O", vec![1, 2, 0, 1, 1]),
            ("[Na+].[O-]C=O", vec![3, 1, 0, 2]),
            ("[CH3:7]CO", vec![1, 0, 2]),
            ("[13CH3]C", vec![1, 0]),
            ("[12CH3][13CH3]", vec![0, 1]),
            ("[35Cl]C[37Cl]", vec![1, 0, 2]),
            ("[127I]C[129I]", vec![1, 0, 2]),
            ("[128Te]C[130Te]", vec![1, 0, 2]),
            ("C1/C=C\\C2CCCCC2C1", vec![5, 8, 9, 7, 4, 1, 0, 2, 6, 3]),
        ] {
            let record = parse_smiles(input, &Default::default()).unwrap();
            let valence = cosmolkit_core::assign_valence_with_options_for_topology(
                &record.topology,
                ValenceModel::RdkitLike,
                false,
            )
            .unwrap();
            assert_eq!(
                cosmolkit_core::assign_atom_cip_ranks(&record.topology, &valence).unwrap(),
                expected,
                "{input}"
            );
        }
    }

    #[test]
    fn ring_label_format_matches_smiles_writer_extensions() {
        let mut output = cosmolkit_model::PropertyText::new();
        write_ring_label(&mut output, 9);
        write_ring_label(&mut output, 10);
        write_ring_label(&mut output, 100);
        assert_eq!(output.as_bytes(), b"9%10%(100)");
    }

    #[test]
    fn canonical_writer_matches_rdkit_for_non_stereo_components_and_cycles() {
        for (input, expected) in [
            ("OCC", "CCO"),
            ("C(C)O", "CCO"),
            ("OC.C", "C.CO"),
            ("C.O.CC", "C.CC.O"),
            ("C1CCCCC1", "C1CCCCC1"),
            ("C1CCC2CCCCC2C1", "C1CCC2CCCCC2C1"),
            ("C12(CCCCC1)CCCCC2", "C1CCC2(CC1)CCCCC2"),
            ("C1C2C3C1C2C3", "C1C2C3CC2C13"),
            ("OC(C)C", "CC(C)O"),
            ("n1ccccc1", "c1ccncc1"),
            ("N#C", "C#N"),
            ("[CH3:7]CO", "OC[CH3:7]"),
            ("[Na+]Cl", "[Na+][Cl]"),
            ("CCN(CC)CC", "CCN(CC)CC"),
            ("CC(C)(C)O", "CC(C)(C)O"),
            ("O=C(O)C", "CC(=O)O"),
            ("C#CC=C", "C#CC=C"),
            ("N=C=O", "N=C=O"),
            ("[NH4+]", "[NH4+]"),
            ("[O-]C=O", "O=C[O-]"),
            ("[13CH3]C", "C[13CH3]"),
            ("[*:1]CC", "CC[*:1]"),
            ("B(O)O", "OBO"),
            ("OP(=O)(O)O", "O=P(O)(O)O"),
            ("c1cc[nH]c1", "c1cc[nH]c1"),
            ("[nH]1cccc1", "c1cc[nH]c1"),
            ("c1ccccc1O", "Oc1ccccc1"),
            ("C1=CC2=CC=CC=C2C=C1", "C1=CC=C2C=CC=CC2=C1"),
            ("C1CC1C2CC2", "C1CC1C1CC1"),
            ("N->B", "B<-N"),
            ("N1CC->1", "C1C->N1"),
            ("BrCCl", "ClCBr"),
            ("[SiH4]", "[SiH4]"),
            ("[AsH]", "[AsH]"),
            ("[se]1cccc1", "c1cc[se]c1"),
            ("[C]", "C"),
            ("[CH3]", "[CH3]"),
            ("[NH2-]", "[NH2-]"),
            ("[O:2]=[C:1]O", "O[C:1]=[O:2]"),
        ] {
            let record = parse_smiles(input, &Default::default()).unwrap();
            assert_eq!(write_smiles(&record).unwrap(), expected.into(), "{input}");
        }
    }

    #[test]
    fn writer_subset_and_clone_prune_preserve_source_property_and_group_rules() {
        let mut record =
            parse_smiles("C[C@H](F)Cl.C[C@H](O)Br", &SmilesParseParams::default()).unwrap();
        record.topology.atoms[1]
            .set_prop("ordinary_atom", "kept")
            .unwrap();
        record.topology.atoms[1]
            .set_computed_prop("_computed_atom", "cleared")
            .unwrap();
        record.topology.bonds[0]
            .set_prop("ordinary_bond", "kept")
            .unwrap();
        record.topology.bonds[0]
            .set_computed_prop("_computed_bond", "cleared")
            .unwrap();
        record.topology.stereo_groups.push(
            StereoGroup::new(
                cosmolkit_model::StereoGroupKind::Or,
                vec![AtomId::new(1), AtomId::new(5)],
                vec![BondId::new(0), BondId::new(3)],
            )
            .expect("valid distinct stereo members")
            .with_id(17),
        );
        let before = record.clone();
        let inside = [true, true, true, true, false, false, false, false];

        let subset =
            extract_writer_subset_fragment(&record.topology, &[0, 1, 2, 3], &inside).unwrap();
        assert_eq!(
            subset
                .source_atoms
                .iter()
                .map(|atom| atom.index())
                .collect::<Vec<_>>(),
            [0, 1, 2, 3]
        );
        assert_eq!(
            subset
                .source_bonds
                .iter()
                .map(|bond| bond.index())
                .collect::<Vec<_>>(),
            [0, 1, 2]
        );
        assert_eq!(
            string_property(subset.topology.atoms[1].prop("ordinary_atom")),
            Some("kept")
        );
        assert_eq!(subset.topology.atoms[1].prop("_computed_atom"), None);
        assert_eq!(
            string_property(subset.topology.bonds[0].prop("ordinary_bond")),
            Some("kept")
        );
        assert_eq!(subset.topology.bonds[0].prop("_computed_bond"), None);
        assert_eq!(subset.topology.stereo_groups.len(), 1);
        assert_eq!(subset.topology.stereo_groups[0].id(), Some(17));
        assert_eq!(
            subset.topology.stereo_groups[0].kind(),
            cosmolkit_model::StereoGroupKind::Or
        );
        assert_eq!(subset.topology.stereo_groups[0].atoms(), &[AtomId::new(1)]);
        assert_eq!(subset.topology.stereo_groups[0].bonds(), &[BondId::new(0)]);

        let clone_prune = extract_writer_preserving_fragment(&record.topology, &inside).unwrap();
        assert_eq!(
            clone_prune
                .source_atoms
                .iter()
                .map(|atom| atom.index())
                .collect::<Vec<_>>(),
            [0, 1, 2, 3]
        );
        assert_eq!(
            clone_prune
                .source_bonds
                .iter()
                .map(|bond| bond.index())
                .collect::<Vec<_>>(),
            [0, 1, 2]
        );
        assert_eq!(
            string_property(clone_prune.topology.atoms[1].prop("ordinary_atom")),
            Some("kept")
        );
        assert_eq!(clone_prune.topology.atoms[1].prop("_computed_atom"), None);
        assert!(
            !clone_prune.topology.atoms[1]
                .is_prop_computed("_computed_atom")
                .unwrap()
        );
        assert_eq!(
            string_property(clone_prune.topology.bonds[0].prop("ordinary_bond")),
            Some("kept")
        );
        assert_eq!(clone_prune.topology.bonds[0].prop("_computed_bond"), None);
        assert!(
            !clone_prune.topology.bonds[0]
                .is_prop_computed("_computed_bond")
                .unwrap()
        );
        assert_eq!(clone_prune.topology.stereo_groups.len(), 1);
        assert_eq!(clone_prune.topology.stereo_groups[0].id(), Some(17));
        assert_eq!(
            clone_prune.topology.stereo_groups[0].kind(),
            cosmolkit_model::StereoGroupKind::Or
        );
        assert_eq!(
            clone_prune.topology.stereo_groups[0].atoms(),
            &[AtomId::new(1)]
        );
        assert_eq!(
            clone_prune.topology.stereo_groups[0].bonds(),
            &[BondId::new(0)]
        );
        assert_eq!(
            record, before,
            "fragment extraction leaves its source unchanged"
        );
    }

    #[test]
    fn writer_singleton_subset_clears_computed_props_and_keeps_partial_group_membership() {
        let mut record = parse_smiles("C.CC", &SmilesParseParams::default()).unwrap();
        record.topology.atoms[0]
            .set_prop("ordinary_atom", "kept")
            .unwrap();
        record.topology.atoms[0]
            .set_computed_prop("_computed_atom", "cleared")
            .unwrap();
        record.topology.substance_groups.push(
            cosmolkit_model::SubstanceGroup::new(
                cosmolkit_model::SubstanceGroupId::new(0),
                cosmolkit_model::SubstanceGroupKind::Data,
            )
            .with_atoms(vec![AtomId::new(0)])
            .with_label("singleton"),
        );
        record.topology.stereo_groups.push(
            StereoGroup::new(
                cosmolkit_model::StereoGroupKind::Or,
                vec![AtomId::new(0), AtomId::new(1)],
                Vec::new(),
            )
            .expect("valid distinct stereo members")
            .with_id(23),
        );
        let singleton =
            extract_writer_subset_fragment(&record.topology, &[0], &[true, false, false]).unwrap();

        assert_eq!(singleton.source_atoms, [AtomId::new(0)]);
        assert!(singleton.source_bonds.is_empty());
        assert_eq!(
            string_property(singleton.topology.atoms[0].prop("ordinary_atom")),
            Some("kept")
        );
        assert_eq!(singleton.topology.atoms[0].prop("_computed_atom"), None);
        assert_eq!(singleton.topology.substance_groups.len(), 1);
        assert_eq!(
            singleton.topology.substance_groups[0].id(),
            cosmolkit_model::SubstanceGroupId::new(0)
        );
        assert_eq!(
            singleton.topology.substance_groups[0].atoms(),
            &[AtomId::new(0)]
        );
        assert_eq!(
            singleton.topology.substance_groups[0].label(),
            Some(&"singleton".into())
        );
        assert_eq!(singleton.topology.stereo_groups.len(), 1);
        assert_eq!(singleton.topology.stereo_groups[0].id(), Some(23));
        assert_eq!(
            singleton.topology.stereo_groups[0].kind(),
            cosmolkit_model::StereoGroupKind::Or
        );
        assert_eq!(
            singleton.topology.stereo_groups[0].atoms(),
            &[AtomId::new(0)]
        );
        assert!(singleton.topology.stereo_groups[0].bonds().is_empty());
    }

    #[test]
    fn writer_stereo_merge_remaps_signed_ring_relative_ids_and_only_computed_props() {
        let mut target = parse_smiles("CCCCCC", &SmilesParseParams::default()).unwrap();
        let mut assigned = parse_smiles("CCCC", &SmilesParseParams::default()).unwrap();
        target.topology.atoms[2]
            .set_prop("target_atom", "source")
            .unwrap();
        target.topology.atoms[2]
            .set_computed_prop("_stale_atom", "remove")
            .unwrap();
        target.topology.bonds[2]
            .set_prop("target_bond", "source")
            .unwrap();
        target.topology.bonds[2]
            .set_computed_prop("_stale_bond", "remove")
            .unwrap();
        assigned.topology.atoms[0]
            .set_prop("incoming_atom_ordinary", "do-not-copy")
            .unwrap();
        assigned.topology.atoms[0]
            .set_computed_prop("_ringStereoAtoms", vec![-2_i32, 4])
            .unwrap();
        assigned.topology.atoms[0]
            .set_computed_prop("_incoming_atom_computed", "copy")
            .unwrap();
        assigned.topology.bonds[0]
            .set_prop("incoming_bond_ordinary", "do-not-copy")
            .unwrap();
        assigned.topology.bonds[0]
            .set_computed_prop("_incoming_bond_computed", "copy")
            .unwrap();
        assigned.topology.bonds[0].set_direction(BondDirection::EndUpRight);

        merge_writer_stereo_fragment(
            &mut target.topology,
            &assigned.topology,
            &[
                AtomId::new(2),
                AtomId::new(3),
                AtomId::new(4),
                AtomId::new(5),
            ],
            &[BondId::new(2), BondId::new(3), BondId::new(4)],
        )
        .unwrap();

        assert_eq!(
            string_property(target.topology.atoms[2].prop("target_atom")),
            Some("source")
        );
        assert_eq!(
            target.topology.atoms[2].prop("incoming_atom_ordinary"),
            None
        );
        assert_eq!(target.topology.atoms[2].prop("_stale_atom"), None);
        assert_eq!(
            target.topology.atoms[2].prop("_ringStereoAtoms"),
            Some(&PropertyValue::IntVector(vec![-4_i32, 6]))
        );
        assert!(
            target.topology.atoms[2]
                .is_prop_computed("_ringStereoAtoms")
                .unwrap()
        );
        assert_eq!(
            string_property(target.topology.atoms[2].prop("_incoming_atom_computed")),
            Some("copy")
        );
        assert!(
            target.topology.atoms[2]
                .is_prop_computed("_incoming_atom_computed")
                .unwrap()
        );
        assert_eq!(
            string_property(target.topology.bonds[2].prop("target_bond")),
            Some("source")
        );
        assert_eq!(
            target.topology.bonds[2].prop("incoming_bond_ordinary"),
            None
        );
        assert_eq!(target.topology.bonds[2].prop("_stale_bond"), None);
        assert_eq!(
            string_property(target.topology.bonds[2].prop("_incoming_bond_computed")),
            Some("copy")
        );
        assert!(
            target.topology.bonds[2]
                .is_prop_computed("_incoming_bond_computed")
                .unwrap()
        );
        assert_eq!(
            target.topology.bonds[2].direction(),
            BondDirection::EndUpRight
        );
    }

    #[test]
    fn canonical_nonisomeric_stereo_fallback_is_deferred_until_after_ranking() {
        let mut record = parse_smiles("F[C@H](Cl)Br", &SmilesParseParams::default()).unwrap();
        for key in [
            "_CIPCode",
            "_CIPRank",
            "_ChiralityPossible",
            "_ringStereoAtoms",
        ] {
            record.topology.atoms[1].clear_prop(key);
        }
        let params = SmilesWriteParams {
            isomeric_smiles: false,
            canonical: true,
            ..Default::default()
        };
        let components = connected_components(&record.topology);

        prepare_writer_stereochemistry(&mut record.topology, None, None, &params, &components)
            .unwrap();

        assert!(
            record
                .topology
                .atoms
                .iter()
                .all(|atom| atom.prop("_CIPCode").is_none()),
            "canonical non-isomeric stereo assignment belongs after canonical ranking"
        );
    }

    fn assert_s33_standard_component_order_row(
        record: &SmilesRecord,
        params: &SmilesWriteParams,
        expected_text: &str,
        expected_atom_order: &[usize],
        expected_bond_order: &[usize],
        case: &str,
    ) {
        let before = record.clone();
        let output = write_smiles_output(record, params, false).expect(case);
        assert_eq!(output.text, expected_text.into(), "{case}: text");
        assert_eq!(
            output.atom_order,
            expected_atom_order
                .iter()
                .copied()
                .map(AtomId::new)
                .collect::<Vec<_>>(),
            "{case}: atom_order"
        );
        assert_eq!(
            output.bond_order,
            expected_bond_order
                .iter()
                .copied()
                .map(BondId::new)
                .collect::<Vec<_>>(),
            "{case}: bond_order"
        );
        assert_eq!(record, &before, "{case}: input mutation");
    }

    fn s33_standard_params(isomeric: bool) -> SmilesWriteParams {
        SmilesWriteParams {
            isomeric_smiles: isomeric,
            kekule: false,
            canonical: true,
            clean_stereo: true,
            rooted_at_atom: None,
            all_bonds_explicit: false,
            all_hydrogens_explicit: false,
            include_dative_bonds: true,
            ignore_atom_map_numbers: false,
        }
    }

    fn parse_s33_standard_input(input: &str) -> SmilesRecord {
        let parse_params = SmilesParseParams::default();
        let parsed = parse_smiles(input, &parse_params).unwrap();
        let remove_params = cosmolkit_core::RemoveHsParams {
            update_explicit_count: true,
            sanitize: parse_params.sanitize,
            ..cosmolkit_core::RemoveHsParams::default()
        };
        let prepared = cosmolkit_core::remove_hydrogens_with_params(
            parsed.topology,
            parsed.coordinates,
            parsed.properties,
            &remove_params,
        )
        .unwrap();
        // Move the real RH-prepared ring state into the finalizer carrier
        // (fourth argument); the root constructor installs the final carrier
        // in its validated live cache, while this writer helper only borrows
        // the prepared state for its detached record.
        let mut prepared_rings = prepared.final_rings;
        crate::finalize_smiles_stereo(
            SmilesRecord {
                topology: prepared.topology,
                coordinates: prepared.coordinates,
                properties: prepared.properties,
            },
            &parse_params,
            &mut None,
            &mut prepared_rings,
        )
        .unwrap()
    }

    fn s33_interleave_oc_co_components(record: &SmilesRecord) -> SmilesRecord {
        assert_eq!(record.topology.atoms.len(), 4);
        assert_eq!(record.topology.bonds.len(), 2);
        let new_to_old = [0, 2, 1, 3];
        let mut old_to_new = [0; 4];
        for (new_index, old_index) in new_to_old.into_iter().enumerate() {
            old_to_new[old_index] = new_index;
        }

        let mut interleaved = record.clone();
        interleaved.topology.atoms = new_to_old
            .into_iter()
            .enumerate()
            .map(|(new_index, old_index)| {
                record.topology.atoms[old_index]
                    .clone()
                    .with_id(AtomId::new(new_index))
            })
            .collect();
        for bond in &mut interleaved.topology.bonds {
            bond.set_endpoints(
                AtomId::new(old_to_new[bond.begin().index()]),
                AtomId::new(old_to_new[bond.end().index()]),
            );
        }
        interleaved.topology.adjacency = cosmolkit_model::AdjacencyList::from_topology(
            interleaved.topology.atoms.len(),
            &interleaved.topology.bonds,
        );
        interleaved.topology.validate().unwrap();
        interleaved
    }

    #[test]
    fn s33_standard_component_order_matches_all_pinned_text_and_map_rows() {
        // The ten rows are the pinned standard-writer oracle in
        // S33_component_order_oracle.md. The CX profile is deliberately not used.
        let unmapped = [
            (
                "C[C@H](F)Cl.[13CH3]O",
                true,
                "C[C@H](F)Cl.[13CH3]O",
                &[0, 1, 2, 3, 4, 5][..],
                &[0, 1, 2, 3][..],
                "unmapped isomeric source order",
            ),
            (
                "[13CH3]O.C[C@H](F)Cl",
                true,
                "C[C@H](F)Cl.[13CH3]O",
                &[2, 3, 4, 5, 0, 1][..],
                &[1, 2, 3, 0][..],
                "unmapped isomeric reversed components",
            ),
            (
                "C[C@H](F)Cl.[13CH3]O",
                false,
                "CC(F)Cl.CO",
                &[0, 1, 2, 3, 4, 5][..],
                &[0, 1, 2, 3][..],
                "unmapped nonisomeric source order",
            ),
            (
                "[13CH3]O.C[C@H](F)Cl",
                false,
                "CC(F)Cl.CO",
                &[2, 3, 4, 5, 0, 1][..],
                &[1, 2, 3, 0][..],
                "unmapped nonisomeric reversed components",
            ),
        ];
        for (input, isomeric, text, atoms, bonds, case) in unmapped {
            let record = parse_s33_standard_input(input);
            assert_s33_standard_component_order_row(
                &record,
                &s33_standard_params(isomeric),
                text,
                atoms,
                bonds,
                case,
            );
        }

        let mapped = [
            (
                "[CH3:10][C@H:11]([F:12])[Cl:13].[13CH3:20][O:21]",
                true,
                "[13CH3:20][O:21].[CH3:10][C@H:11]([F:12])[Cl:13]",
                &[4, 5, 0, 1, 2, 3][..],
                &[3, 0, 1, 2][..],
                "mapped isomeric source order",
            ),
            (
                "[13CH3:20][O:21].[CH3:10][C@H:11]([F:12])[Cl:13]",
                true,
                "[13CH3:20][O:21].[CH3:10][C@H:11]([F:12])[Cl:13]",
                &[0, 1, 2, 3, 4, 5][..],
                &[0, 1, 2, 3][..],
                "mapped isomeric reversed components",
            ),
            (
                "[CH3:10][C@H:11]([F:12])[Cl:13].[13CH3:20][O:21]",
                false,
                "[CH3:10][CH:11]([F:12])[Cl:13].[CH3:20][O:21]",
                &[0, 1, 2, 3, 4, 5][..],
                &[0, 1, 2, 3][..],
                "mapped nonisomeric source order",
            ),
            (
                "[13CH3:20][O:21].[CH3:10][C@H:11]([F:12])[Cl:13]",
                false,
                "[CH3:10][CH:11]([F:12])[Cl:13].[CH3:20][O:21]",
                &[2, 3, 4, 5, 0, 1][..],
                &[1, 2, 3, 0][..],
                "mapped nonisomeric reversed components",
            ),
        ];
        for (input, isomeric, text, atoms, bonds, case) in mapped {
            let record = parse_s33_standard_input(input);
            assert_s33_standard_component_order_row(
                &record,
                &s33_standard_params(isomeric),
                text,
                atoms,
                bonds,
                case,
            );
        }

        let source_order = parse_s33_standard_input("OC.CO");
        assert_s33_standard_component_order_row(
            &source_order,
            &s33_standard_params(true),
            "CO.CO",
            &[1, 0, 2, 3],
            &[0, 1],
            "equal-text components in source order",
        );
        let interleaved = s33_interleave_oc_co_components(&source_order);
        assert_eq!(
            interleaved
                .topology
                .bonds
                .iter()
                .map(|bond| (bond.begin().index(), bond.end().index()))
                .collect::<Vec<_>>(),
            [(0, 2), (1, 3)],
            "the [0,2,1,3] input permutation interleaves component rows"
        );
        assert_s33_standard_component_order_row(
            &interleaved,
            &s33_standard_params(true),
            "CO.CO",
            &[1, 3, 2, 0],
            &[1, 0],
            "equal-text interleaved components use atom/bond map tie-break",
        );
    }
}

#[cfg(test)]
mod enhanced_stereo_canonical_tests {
    use std::collections::{BTreeMap, BTreeSet};

    use cosmolkit_model::{
        AtomId, BondId, StereoGroup, StereoGroupKind, TopologyBlock, set_stereo_group_write_id,
        stereo_group_write_id,
    };
    use cosmolkit_types::{BondStereo, ChiralTag};

    use super::{SmilesWriteParams, canonicalize_enhanced_stereo, write_fragment_smiles_output};
    use crate::{SmilesParseParams, parse_smiles};

    fn topology(input: &str) -> TopologyBlock {
        parse_smiles(input, &SmilesParseParams::default())
            .expect("parse fixed enhanced-stereo topology")
            .topology
    }

    #[test]
    fn enhanced_stereo_canonical_no_groups_and_absolute_group_are_unchanged() {
        let mut topology = topology("CCC");
        topology.atoms[0].set_chiral_tag(ChiralTag::TetrahedralCw);
        topology.bonds[0]
            .set_stereo(BondStereo::AtropCw)
            .expect("atrop state is valid without reference atoms");

        let no_groups_before = topology.clone();
        assert_eq!(
            canonicalize_enhanced_stereo(&mut topology, &[0, 1, 2]).unwrap(),
            BTreeMap::new()
        );
        assert_eq!(topology, no_groups_before);

        let mut absolute = StereoGroup::new(
            StereoGroupKind::Absolute,
            vec![AtomId::new(0)],
            vec![BondId::new(0)],
        )
        .expect("valid distinct stereo members")
        .with_id(37);
        set_stereo_group_write_id(&mut absolute, 91);
        topology.stereo_groups.push(absolute.clone());
        let atoms_before = topology.atoms.clone();
        let bonds_before = topology.bonds.clone();

        assert_eq!(
            canonicalize_enhanced_stereo(&mut topology, &[0, 1, 2]).unwrap(),
            BTreeMap::new()
        );
        assert_eq!(topology.stereo_groups, [absolute]);
        assert_eq!(topology.atoms, atoms_before);
        assert_eq!(topology.bonds, bonds_before);
        assert_eq!(topology.stereo_groups[0].id(), Some(37));
        assert_eq!(stereo_group_write_id(&topology.stereo_groups[0]), 91);
    }

    #[test]
    fn enhanced_stereo_canonical_atom_groups_sort_invert_and_write_source_indices() {
        let mut topology = topology("CCCC");
        topology.atoms[0].set_chiral_tag(ChiralTag::TetrahedralCcw);
        topology.atoms[1].set_chiral_tag(ChiralTag::TetrahedralCw);
        topology.atoms[2].set_chiral_tag(ChiralTag::TetrahedralCcw);
        topology.atoms[3].set_chiral_tag(ChiralTag::TetrahedralCw);

        let mut absolute = StereoGroup::new(StereoGroupKind::Absolute, Vec::new(), Vec::new())
            .expect("valid distinct stereo members")
            .with_id(8);
        set_stereo_group_write_id(&mut absolute, 15);
        let mut or_group = StereoGroup::new(
            StereoGroupKind::Or,
            vec![AtomId::new(0), AtomId::new(1)],
            Vec::new(),
        )
        .expect("valid distinct stereo members")
        .with_id(21);
        set_stereo_group_write_id(&mut or_group, 22);
        let mut and_group = StereoGroup::new(
            StereoGroupKind::And,
            vec![AtomId::new(2), AtomId::new(3)],
            Vec::new(),
        )
        .expect("valid distinct stereo members")
        .with_id(31);
        set_stereo_group_write_id(&mut and_group, 32);
        topology.stereo_groups = vec![absolute.clone(), or_group, and_group];

        let references = canonicalize_enhanced_stereo(&mut topology, &[8, 1, 7, 2]).unwrap();

        assert_eq!(references, BTreeMap::from([(1, 1), (3, 2)]));
        assert_eq!(topology.stereo_groups[0], absolute);
        assert_eq!(topology.stereo_groups[0].id(), Some(8));
        assert_eq!(stereo_group_write_id(&topology.stereo_groups[0]), 15);
        assert_eq!(topology.stereo_groups[1].kind(), StereoGroupKind::Or);
        assert_eq!(
            topology.stereo_groups[1].atoms(),
            &[AtomId::new(1), AtomId::new(0)]
        );
        assert_eq!(topology.stereo_groups[1].id(), None);
        assert_eq!(stereo_group_write_id(&topology.stereo_groups[1]), 0);
        assert_eq!(topology.stereo_groups[2].kind(), StereoGroupKind::And);
        assert_eq!(
            topology.stereo_groups[2].atoms(),
            &[AtomId::new(3), AtomId::new(2)]
        );
        assert_eq!(topology.stereo_groups[2].id(), None);
        assert_eq!(stereo_group_write_id(&topology.stereo_groups[2]), 0);
        assert_eq!(topology.atoms[0].chiral_tag(), ChiralTag::TetrahedralCw);
        assert_eq!(topology.atoms[1].chiral_tag(), ChiralTag::TetrahedralCcw);
        assert_eq!(topology.atoms[2].chiral_tag(), ChiralTag::TetrahedralCw);
        assert_eq!(topology.atoms[3].chiral_tag(), ChiralTag::TetrahedralCcw);
    }

    #[test]
    fn enhanced_stereo_canonical_atom_rank_ties_remain_source_equivalent() {
        let mut topology = topology("CCCC");
        for atom in &mut topology.atoms {
            atom.set_chiral_tag(ChiralTag::TetrahedralCcw);
        }
        topology.stereo_groups.push(
            StereoGroup::new(
                StereoGroupKind::Or,
                vec![AtomId::new(3), AtomId::new(1), AtomId::new(2)],
                Vec::new(),
            )
            .expect("valid distinct stereo members"),
        );

        let references = canonicalize_enhanced_stereo(&mut topology, &[9, 4, 4, 1]).unwrap();
        let members = topology.stereo_groups[0].atoms();

        assert_eq!(members[0], AtomId::new(3));
        assert_eq!(
            members[1..].iter().copied().collect::<BTreeSet<_>>(),
            BTreeSet::from([AtomId::new(1), AtomId::new(2)])
        );
        assert_eq!(
            members
                .windows(2)
                .map(|pair| [9, 4, 4, 1][pair[0].index()] <= [9, 4, 4, 1][pair[1].index()])
                .all(|ordered| ordered),
            true,
            "rank-equivalent atoms have no source-defined secondary order"
        );
        assert_eq!(references, BTreeMap::from([(3, 0)]));
        assert!(
            topology
                .atoms
                .iter()
                .all(|atom| atom.chiral_tag() == ChiralTag::TetrahedralCcw)
        );
    }

    #[test]
    fn enhanced_stereo_canonical_already_ccw_and_empty_groups_write_only_atom_marker() {
        let mut topology = topology("CC");
        topology.atoms[0].set_chiral_tag(ChiralTag::TetrahedralCcw);
        topology.atoms[1].set_chiral_tag(ChiralTag::TetrahedralCw);
        let mut atom_group = StereoGroup::new(
            StereoGroupKind::Or,
            vec![AtomId::new(1), AtomId::new(0)],
            Vec::new(),
        )
        .expect("valid distinct stereo members")
        .with_id(6);
        set_stereo_group_write_id(&mut atom_group, 7);
        let empty_group = StereoGroup::new(StereoGroupKind::And, Vec::new(), Vec::new())
            .expect("valid distinct stereo members")
            .with_id(8);
        topology.stereo_groups = vec![atom_group, empty_group];

        let references = canonicalize_enhanced_stereo(&mut topology, &[0, 1]).unwrap();

        assert_eq!(references, BTreeMap::from([(0, 0)]));
        assert_eq!(topology.atoms[0].chiral_tag(), ChiralTag::TetrahedralCcw);
        assert_eq!(topology.atoms[1].chiral_tag(), ChiralTag::TetrahedralCw);
        assert_eq!(
            topology.stereo_groups[0].atoms(),
            &[AtomId::new(0), AtomId::new(1)]
        );
        assert_eq!(topology.stereo_groups[0].id(), None);
        assert_eq!(stereo_group_write_id(&topology.stereo_groups[0]), 0);
        assert!(topology.stereo_groups[1].is_empty());
        assert_eq!(topology.stereo_groups[1].kind(), StereoGroupKind::And);
        assert_eq!(topology.stereo_groups[1].id(), None);
        assert_eq!(stereo_group_write_id(&topology.stereo_groups[1]), 0);
    }

    #[test]
    fn enhanced_stereo_canonical_bond_pairs_sort_with_ties_and_invert_full_group() {
        let mut topology = topology("CCCCC");
        for (bond, stereo) in [
            (0, BondStereo::AtropCcw),
            (1, BondStereo::None),
            (2, BondStereo::AtropCcw),
            (3, BondStereo::AtropCw),
        ] {
            topology.bonds[bond]
                .set_stereo(stereo)
                .expect("source bond state is valid");
        }
        let mut group = StereoGroup::new(
            StereoGroupKind::Or,
            Vec::new(),
            (0..4).map(BondId::new).collect(),
        )
        .expect("valid distinct stereo members")
        .with_id(54);
        set_stereo_group_write_id(&mut group, 55);
        topology.stereo_groups.push(group);

        let references = canonicalize_enhanced_stereo(&mut topology, &[10, 0, 8, 0, 5]).unwrap();
        let ordered_bonds = topology.stereo_groups[0].bonds();
        let rank_pair = |bond_id: BondId| {
            let bond = &topology.bonds[bond_id.index()];
            let first = [10, 0, 8, 0, 5][bond.begin().index()];
            let second = [10, 0, 8, 0, 5][bond.end().index()];
            (first.max(second), first.min(second))
        };

        assert!(
            references.is_empty(),
            "bond-only groups write no atom marker"
        );
        assert_eq!(ordered_bonds[0], BondId::new(3));
        assert_eq!(ordered_bonds[3], BondId::new(0));
        assert_eq!(
            ordered_bonds[1..3].iter().copied().collect::<BTreeSet<_>>(),
            BTreeSet::from([BondId::new(1), BondId::new(2)])
        );
        assert!(
            ordered_bonds
                .windows(2)
                .all(|pair| rank_pair(pair[0]) <= rank_pair(pair[1]))
        );
        assert_eq!(topology.bonds[0].stereo(), BondStereo::AtropCw);
        assert_eq!(topology.bonds[1].stereo(), BondStereo::None);
        assert_eq!(topology.bonds[2].stereo(), BondStereo::AtropCw);
        assert_eq!(topology.bonds[3].stereo(), BondStereo::AtropCcw);
        assert_eq!(topology.stereo_groups[0].id(), None);
        assert_eq!(stereo_group_write_id(&topology.stereo_groups[0]), 0);
    }

    #[test]
    fn enhanced_stereo_canonical_atom_reference_precedes_bond_reference() {
        let mut topology = topology("CC");
        topology.atoms[0].set_chiral_tag(ChiralTag::TetrahedralCw);
        topology.atoms[1].set_chiral_tag(ChiralTag::TetrahedralCcw);
        topology.bonds[0]
            .set_stereo(BondStereo::AtropCcw)
            .expect("atrop state is valid without reference atoms");
        topology.stereo_groups.push(
            StereoGroup::new(
                StereoGroupKind::And,
                vec![AtomId::new(0)],
                vec![BondId::new(0)],
            )
            .expect("valid distinct stereo members"),
        );

        let references = canonicalize_enhanced_stereo(&mut topology, &[0, 1]).unwrap();

        assert_eq!(references, BTreeMap::from([(0, 0)]));
        assert_eq!(topology.atoms[0].chiral_tag(), ChiralTag::TetrahedralCcw);
        assert_eq!(topology.atoms[1].chiral_tag(), ChiralTag::TetrahedralCcw);
        assert_eq!(topology.bonds[0].stereo(), BondStereo::AtropCw);
        assert_eq!(topology.stereo_groups[0].atoms(), &[AtomId::new(0)]);
        assert_eq!(topology.stereo_groups[0].bonds(), &[BondId::new(0)]);
    }

    #[test]
    fn enhanced_stereo_canonical_non_ccw_non_atrop_bond_reference_uses_cw_fallback() {
        let mut topology = topology("CCC");
        topology.bonds[0]
            .set_stereo(BondStereo::None)
            .expect("ordinary bond stereo is valid");
        topology.bonds[1]
            .set_stereo(BondStereo::AtropCcw)
            .expect("atrop state is valid without reference atoms");
        topology.stereo_groups.push(
            StereoGroup::new(
                StereoGroupKind::Or,
                Vec::new(),
                vec![BondId::new(0), BondId::new(1)],
            )
            .expect("valid distinct stereo members"),
        );

        let references = canonicalize_enhanced_stereo(&mut topology, &[0, 1, 2]).unwrap();

        assert!(references.is_empty());
        assert_eq!(topology.bonds[0].stereo(), BondStereo::None);
        assert_eq!(topology.bonds[1].stereo(), BondStereo::AtropCw);
    }

    #[test]
    fn enhanced_stereo_canonical_nontetrahedral_reference_uses_source_group_state() {
        let record = crate::parse_smiles(
            "[P@TB1](F)(Cl)(Br)(I)N |o1:0|",
            &SmilesParseParams::default(),
        )
        .expect("parse pinned non-tetrahedral enhanced-stereo example");

        let output = super::write_smiles_for_cx(&record, &SmilesWriteParams::default())
            .expect("write canonical CX base SMILES");

        assert_eq!(output.text, "N[P@TB2](F)(Cl)(Br)I".into());

        let atoms = (0..record.topology.atoms.len())
            .map(AtomId::new)
            .collect::<Vec<_>>();
        let fragment = write_fragment_smiles_output(
            &record,
            &SmilesWriteParams::default(),
            &atoms,
            None,
            None,
            None,
            None,
            None,
        )
        .expect("write canonical detached fragment");

        assert_eq!(fragment.text, "N[P@TB2](F)(Cl)(Br)I".into());
    }

    #[test]
    fn enhanced_stereo_canonical_fragment_entry_preserves_input_record() {
        let mut record = parse_smiles("C[C@H](N)C[C@@H](N)C", &SmilesParseParams::default())
            .expect("parse fixed detached writer input");
        record.topology.stereo_groups.push(
            StereoGroup::new(
                StereoGroupKind::Or,
                vec![AtomId::new(1), AtomId::new(4)],
                Vec::new(),
            )
            .expect("valid distinct stereo members"),
        );
        let before = record.clone();
        let atoms = (0..record.topology.atoms.len())
            .map(AtomId::new)
            .collect::<Vec<_>>();

        let output = write_fragment_smiles_output(
            &record,
            &SmilesWriteParams::default(),
            &atoms,
            None,
            None,
            None,
            None,
            None,
        )
        .expect("write complete detached fragment");

        assert_eq!(record, before, "detached writer must preserve its input");
        assert_eq!(output.atom_order.len(), atoms.len());
        assert_eq!(
            output
                .text
                .as_bytes()
                .iter()
                .filter(|&&byte| byte == b'.')
                .count(),
            0
        );
    }
}

#[cfg(test)]
mod uint_complete_source_condition_cells {
    use super::*;
    use cosmolkit_model::PropertyValue;
    // FROZEN UINT CONDITION: STRICT_BOOL_VECTOR
    #[test]
    fn uint_cell_strict_bool_vector_writer() {
        for n in [0_u32, 1, 4294967295] {
            assert!(
                matches!(parse_ring_stereo_atoms(&PropertyValue::UInt(n),2),Err(SmilesParseError::WriterStereo(message)) if message=="`_ringStereoAtoms` is not a signed source INT_VECT value (bad_any_cast)")
            );
        }
    }
}

#[cfg(test)]
fn fixed_property_text(value: &cosmolkit_model::PropertyText) -> &str {
    std::str::from_utf8(value.as_bytes()).expect("original fixed fixture text is UTF8")
}

#[cfg(test)]
mod computed_clear_error_tests {
    use super::*;
    use crate::{SmilesParseParams, parse_smiles};

    fn assert_writer_failure(input: &str, bad_bond: bool) {
        for clean_stereo in [false, true] {
            let mut record = parse_smiles(input, &SmilesParseParams::default()).unwrap();
            record.properties.clear_prop("_StereochemDone").unwrap();
            if bad_bond {
                record.topology.bonds[0]
                    .set_prop("__computedProps", PropertyValue::Int(7))
                    .unwrap();
            } else {
                record.topology.atoms[0]
                    .set_prop("__computedProps", PropertyValue::Int(7))
                    .unwrap();
            }
            let before = record.clone();
            let params = SmilesWriteParams {
                clean_stereo,
                ..Default::default()
            };
            let error = write_smiles_with_params(&record, &params).unwrap_err();
            if bad_bond {
                assert!(matches!(
                    &error,
                    SmilesParseError::BondProperty(
                        cosmolkit_model::BondValueError::ComputedListKind(_)
                    )
                ));
                assert!(
                    std::error::Error::source(&error)
                        .unwrap()
                        .downcast_ref::<cosmolkit_model::BondValueError>()
                        .is_some()
                );
            } else {
                assert!(matches!(
                    &error,
                    SmilesParseError::AtomProperty(
                        cosmolkit_model::AtomPropertyError::ComputedListKind(_)
                    )
                ));
                assert!(
                    std::error::Error::source(&error)
                        .unwrap()
                        .downcast_ref::<cosmolkit_model::AtomPropertyError>()
                        .is_some()
                );
            }
            assert_eq!(
                record, before,
                "failed serialization preserves all input blocks"
            );

            // A genuine empty StringVector is valid source state. It must
            // keep the same writer options, succeed and preserve its input.
            if bad_bond {
                record.topology.bonds[0]
                    .set_prop("__computedProps", PropertyValue::StringVector(Vec::new()))
                    .unwrap();
            } else {
                record.topology.atoms[0]
                    .set_prop("__computedProps", PropertyValue::StringVector(Vec::new()))
                    .unwrap();
            }
            let valid_before = record.clone();
            assert_eq!(
                write_smiles_with_params(&record, &params)
                    .unwrap()
                    .as_bytes(),
                input.as_bytes()
            );
            assert_eq!(record, valid_before);
        }
    }
    #[test]
    fn singleton_subset_writer_retains_computed_list_error_and_input() {
        assert_writer_failure("C.C", false);
    }
    #[test]
    fn preserving_fragment_writer_retains_atom_computed_list_error_and_input() {
        assert_writer_failure("CC.CC", false);
    }
    #[test]
    fn preserving_fragment_writer_retains_bond_computed_list_error_and_input() {
        assert_writer_failure("CC.CC", true);
    }
}

#[cfg(test)]
mod source578_organic_membership_tests {
    use super::in_organic_subset;

    #[test]
    fn source578_native_signed_domain_zero_and_sentinel_boundaries() {
        assert!(in_organic_subset(0));
        assert!(!in_organic_subset(-1));
        assert!(!in_organic_subset(i32::MIN));
        assert!(!in_organic_subset(i32::MAX));
        assert!(!in_organic_subset(256));
        for number in 0..=255 {
            let expected = [0, 5, 6, 7, 8, 9, 15, 16, 17, 35, 53].contains(&number);
            assert_eq!(
                in_organic_subset(number),
                expected,
                "native atomicNumber={number}"
            );
        }
    }
}

#[cfg(test)]
mod source582_dfs_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, BondSpec, Element, QueryAtom, QueryBond, QueryGraph};

    fn graph(n: usize, edges: &[(usize, usize)]) -> QueryGraph {
        QueryGraph::from_parts(
            (0..n)
                .map(|i| QueryAtom::new(AtomId::new(i), AtomSpec::new(Element::C)))
                .collect(),
            edges
                .iter()
                .enumerate()
                .map(|(i, &(a, b))| {
                    QueryBond::new(
                        BondId::new(i),
                        BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                    )
                })
                .collect(),
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }

    struct Scratch {
        colors: Vec<AtomColor>,
        closures: Vec<Vec<BondId>>,
        ids: Vec<Option<usize>>,
        available: Vec<bool>,
        stack: Vec<MolStackElem>,
        orders: Vec<Vec<BondId>>,
        opened: Vec<bool>,
    }
    impl Scratch {
        fn new(n: usize, m: usize) -> Self {
            Self {
                colors: vec![AtomColor::White; n],
                closures: vec![vec![]; n],
                ids: vec![None; m],
                available: vec![true; 1024],
                stack: vec![],
                orders: vec![vec![]; n],
                opened: vec![false; m],
            }
        }
        fn run(
            &mut self,
            g: &mut QueryGraph,
            ranks: &[u32],
            rings: &[bool],
        ) -> Result<(), SmilesParseError> {
            dfs_build_query_stack(
                g,
                0,
                None,
                &mut self.colors,
                ranks,
                rings,
                &self.closures,
                &mut self.ids,
                &mut self.available,
                &mut self.stack,
                &mut self.orders,
                &mut self.opened,
                None,
                None,
                None,
            )
        }
    }

    #[test]
    fn source582_unsigned_rank_converts_to_signed_native_tuple_and_branches() {
        let mut g = graph(3, &[(0, 1), (0, 2)]);
        let mut x = Scratch::new(3, 2);
        x.run(&mut g, &[0, u32::MAX, 0], &[false; 2]).unwrap();
        assert_eq!(
            x.stack,
            vec![
                MolStackElem::Atom(0),
                MolStackElem::BranchOpen(0),
                MolStackElem::Bond {
                    bond: BondId::new(0),
                    atom_to_left: 0
                },
                MolStackElem::Atom(1),
                MolStackElem::BranchClose(0),
                MolStackElem::Bond {
                    bond: BondId::new(1),
                    atom_to_left: 0
                },
                MolStackElem::Atom(2)
            ]
        );
        assert_eq!(
            x.orders,
            vec![
                vec![BondId::new(0), BondId::new(1)],
                vec![BondId::new(0)],
                vec![BondId::new(1)]
            ]
        );
        assert_eq!(x.colors, vec![AtomColor::Black; 3]);
    }

    #[test]
    fn source582_invalid_ring_property_preserves_native_prefix_and_outputs() {
        let mut g = graph(3, &[(0, 1), (1, 2), (2, 0)]);
        let mut x = Scratch::new(3, 3);
        g.bonds_mut()[2]
            .bond_mut()
            .set_prop("_TraversalRingClosureBond", "not_uint")
            .unwrap();
        x.closures[0] = vec![BondId::new(2)];
        x.orders[0] = vec![BondId::new(1)];
        assert!(matches!(
            x.run(&mut g, &[0, 1, 2], &[true; 3]),
            Err(SmilesParseError::WriterNumeric(_))
        ));
        assert_eq!(x.stack, vec![MolStackElem::Atom(0)]);
        assert_eq!(
            x.colors,
            vec![AtomColor::Grey, AtomColor::White, AtomColor::White]
        );
        assert_eq!(x.orders[0], vec![BondId::new(1)]);
        assert!(x.available.iter().all(|v| *v));
        assert_eq!(
            g.bond(2).unwrap().bond().prop("_TraversalRingClosureBond"),
            Some(&PropertyValue::String("not_uint".into()))
        );
    }

    #[test]
    fn source582_child_property_clear_fails_before_branch_bond_and_traversal_append() {
        let mut g = graph(2, &[(0, 1)]);
        let mut x = Scratch::new(2, 1);
        g.atom_mut(1)
            .unwrap()
            .set_prop("_TraversalBondIndexOrder", PropertyValue::UInt(9))
            .unwrap();
        g.atom_mut(1)
            .unwrap()
            .set_prop("__computedProps", PropertyValue::Bool(false))
            .unwrap();
        assert!(matches!(
            x.run(&mut g, &[0, 1], &[false]),
            Err(SmilesParseError::AtomProperty(_))
        ));
        assert_eq!(x.stack, vec![MolStackElem::Atom(0)]);
        assert!(x.orders[0].is_empty());
        assert_eq!(x.colors, vec![AtomColor::Grey, AtomColor::White]);
        assert_eq!(
            g.atom(1).unwrap().prop("_TraversalBondIndexOrder"),
            Some(&PropertyValue::UInt(9))
        );
    }

    #[test]
    fn source582_ring_digits_closed_at_atom_reuse_only_after_all_closures() {
        let mut g = graph(3, &[(0, 1), (1, 2), (2, 0)]);
        let mut x = Scratch::new(3, 3);
        g.bonds_mut()[0]
            .bond_mut()
            .set_prop("_TraversalRingClosureBond", PropertyValue::UInt(1))
            .unwrap();
        x.closures[0] = vec![BondId::new(0), BondId::new(2)];
        x.available[0] = false;
        x.run(&mut g, &[0, 1, 2], &[true; 3]).unwrap();
        assert_eq!(
            x.stack,
            vec![
                MolStackElem::Atom(0),
                MolStackElem::Bond {
                    bond: BondId::new(0),
                    atom_to_left: 0
                },
                MolStackElem::Ring(1),
                MolStackElem::Ring(2)
            ]
        );
        assert!(x.available[0]);
        assert!(!x.available[1]);
        assert_eq!(
            g.bond(2).unwrap().bond().prop("_TraversalRingClosureBond"),
            Some(&PropertyValue::UInt(2))
        );
        assert_eq!(x.orders[0], vec![BondId::new(0), BondId::new(2)]);
    }

    #[test]
    fn source582_ring_capacity_is_source_value_error_without_property_write() {
        let mut g = graph(2, &[(0, 1)]);
        let mut x = Scratch::new(2, 1);
        x.closures[0] = vec![BondId::new(0)];
        x.available.fill(false);
        let error = x.run(&mut g, &[0, 1], &[true]).unwrap_err();
        assert_eq!(error, SmilesParseError::TraversalTooManyOpenRings);
        assert_eq!(
            error.to_string(),
            "Too many rings open at once. SMILES cannot be generated."
        );
        assert_eq!(x.stack, vec![MolStackElem::Atom(0)]);
        assert_eq!(x.colors[0], AtomColor::Grey);
        assert_eq!(
            g.bond(0).unwrap().bond().prop("_TraversalRingClosureBond"),
            None
        );
    }

    #[test]
    fn source582_random_branch_reads_actual_stream_and_skips_ring_and_symbol_inputs() {
        let mut g = graph(3, &[(0, 1), (0, 2)]);
        let mut x = Scratch::new(3, 2);
        let mut expected = cosmolkit_core::RdkitRandomEngine::from_seed(1);
        let a = expected.next_u32();
        let b = expected.next_u32();
        let next = expected.next_u32();
        let empty_symbols: Vec<PropertyText> = vec![];
        let actual = cosmolkit_core::with_rdkit_random_generator(1, |rng| {
            dfs_build_query_stack(
                &mut g,
                0,
                None,
                &mut x.colors,
                &[0; 3],
                &[],
                &x.closures,
                &mut x.ids,
                &mut x.available,
                &mut x.stack,
                &mut x.orders,
                &mut x.opened,
                None,
                Some(DfsBondSymbols::Bytes(&empty_symbols)),
                Some(rng),
            )
            .unwrap();
            rng.next_u32()
        });
        assert_eq!(actual, next);
        assert_eq!(
            x.stack
                .iter()
                .filter_map(|v| if let MolStackElem::Atom(a) = v {
                    Some(*a)
                } else {
                    None
                })
                .collect::<Vec<_>>(),
            if (a as i32) < (b as i32) {
                vec![0, 1, 2]
            } else {
                vec![0, 2, 1]
            }
        );
    }
}

#[cfg(test)]
mod source590_cycle_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, BondSpec, Element, QueryAtom, QueryBond, QueryGraph};
    fn graph(n: usize, edges: &[(usize, usize)]) -> QueryGraph {
        QueryGraph::from_parts(
            (0..n)
                .map(|i| QueryAtom::new(AtomId::new(i), AtomSpec::new(Element::C)))
                .collect(),
            edges
                .iter()
                .enumerate()
                .map(|(i, &(a, b))| {
                    QueryBond::new(
                        BondId::new(i),
                        BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                    )
                })
                .collect(),
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    #[test]
    fn source590_cycle_discovery_preserves_caller_colors_and_builds_reciprocal_closures() {
        let mut g = graph(3, &[(0, 1), (1, 2), (2, 0)]);
        let mut colors = vec![AtomColor::White; 3];
        let mut closures = vec![vec![]; 3];
        let mut orders = vec![vec![]; 3];
        let mut stack = vec![];
        let opened = canonical_dfs_traversal_source(
            &mut g,
            0,
            None,
            &mut colors,
            &[0, 1, 2],
            &[true; 3],
            &mut stack,
            &mut closures,
            &mut orders,
            None,
            None,
            None,
        )
        .unwrap();
        assert_eq!(
            closures,
            vec![vec![BondId::new(2)], vec![], vec![BondId::new(2)]]
        );
        assert_eq!(colors, vec![AtomColor::Black; 3]);
        assert!(opened[2]);
        assert_eq!(
            stack
                .iter()
                .filter_map(|x| if let MolStackElem::Atom(i) = x {
                    Some(*i)
                } else {
                    None
                })
                .collect::<Vec<_>>(),
            vec![0, 1, 2]
        );
        assert_eq!(
            stack
                .iter()
                .filter_map(|x| if let MolStackElem::Ring(i) = x {
                    Some(*i)
                } else {
                    None
                })
                .collect::<Vec<_>>(),
            vec![1, 1]
        );
        assert_eq!(
            g.bond(2).unwrap().bond().prop("_TraversalRingClosureBond"),
            Some(&PropertyValue::UInt(1))
        );
    }
    #[test]
    fn source590_random_cycle_pass_draws_for_black_neighbors_after_mask_and_incoming_filters() {
        let g = graph(4, &[(0, 1), (0, 2), (0, 3)]);
        let mut colors = vec![
            AtomColor::White,
            AtomColor::Black,
            AtomColor::White,
            AtomColor::Black,
        ];
        let mut closures = vec![vec![]; 4];
        let mut expected = cosmolkit_core::RdkitRandomEngine::from_seed(42);
        expected.next_u32();
        let next = expected.next_u32();
        let actual = cosmolkit_core::with_rdkit_random_generator(42, |rng| {
            dfs_find_cycles_source(
                &g,
                0,
                Some(BondId::new(1)),
                &mut colors,
                &[0; 4],
                &[],
                Some(&[true, true, false]),
                Some(DfsBondSymbols::Bytes(&[])),
                &mut closures,
                Some(rng),
            )
            .unwrap();
            rng.next_u32()
        });
        assert_eq!(actual, next);
        assert_eq!(
            colors,
            vec![
                AtomColor::Black,
                AtomColor::Black,
                AtomColor::White,
                AtomColor::Black
            ]
        );
        assert!(closures.iter().all(Vec::is_empty));
    }
    #[test]
    fn source590_dfs_preconditions_precede_color_and_property_mutation() {
        let mut g = graph(2, &[(0, 1)]);
        let mut colors = vec![AtomColor::White; 2];
        let mut closures = vec![vec![]; 2];
        let mut orders = vec![vec![]; 2];
        let mut stack = vec![];
        let result = canonical_dfs_traversal_source(
            &mut g,
            0,
            None,
            &mut colors,
            &[0; 2],
            &[false],
            &mut stack,
            &mut closures,
            &mut orders,
            None,
            Some(DfsBondSymbols::Bytes(&[])),
            None,
        );
        assert!(matches!(
            result,
            Err(SmilesParseError::TraversalStateIndex {
                state: "bondSymbols",
                ..
            })
        ));
        assert_eq!(colors, vec![AtomColor::White; 2]);
        assert!(stack.is_empty());
        assert_eq!(
            g.bond(0).unwrap().bond().prop("_TraversalRingClosureBond"),
            None
        );
    }
}

#[cfg(test)]
mod source590_fragment_tests {
    use super::*;
    use cosmolkit_model::StereoGroupKind;
    use cosmolkit_model::{AtomSpec, BondSpec, Element, MoleculeProperties};
    fn graph(n: usize, edges: &[(usize, usize)]) -> TopologyBlock {
        TopologyBlock::try_from_parts(
            (0..n)
                .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
                .collect(),
            edges
                .iter()
                .enumerate()
                .map(|(i, &(a, b))| {
                    Bond::from_spec(
                        BondId::new(i),
                        BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                    )
                })
                .collect(),
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn valence(g: &TopologyBlock) -> ValenceAssignment {
        cosmolkit_core::assign_valence_with_options_for_topology(g, ValenceModel::RdkitLike, false)
            .unwrap()
    }
    fn rings(n: usize, m: usize) -> cosmolkit_core::RingInfo {
        cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::SymmSssr, n, m)
    }
    fn phase(
        g: &mut TopologyBlock,
        v: &ValenceAssignment,
        r: &cosmolkit_core::RingInfo,
        stack: &[MolStackElem],
        mask: Option<&[bool]>,
        iso: bool,
    ) -> Result<(), SmilesParseError> {
        let n = g.atoms.len();
        let m = g.bonds.len();
        canonicalize_source_stack(
            g,
            v,
            Some(r),
            stack,
            &vec![vec![]; n],
            &vec![vec![]; n],
            &vec![false; m],
            None,
            mask,
            iso,
            true,
            &BTreeMap::new(),
            None,
        )
    }
    #[test]
    fn source590_fragment_preconditions_precede_all_preparation_effects() {
        let mut g = graph(2, &[(0, 1)]);
        let v = valence(&g);
        let mut props = MoleculeProperties::default();
        props.set_prop("__computedProps", false).unwrap();
        let mut r = rings(2, 1);
        let before = g.clone();
        let prop_before = props.clone();
        let ring_before = r.clone();
        let mut stack = vec![];
        assert!(matches!(
            canonicalize_fragment_source(
                &mut g,
                &mut props,
                &v,
                &mut r,
                0,
                &mut [AtomColor::White],
                &[0; 2],
                &mut stack,
                None,
                None,
                None,
                false,
                None,
                true,
                &BTreeMap::new()
            ),
            Err(SmilesParseError::TraversalStateIndex {
                state: "colors",
                ..
            })
        ));
        assert_eq!(g, before);
        assert_eq!(props, prop_before);
        assert_eq!(r, ring_before);
        assert!(stack.is_empty());
    }
    #[test]
    fn source590_stereo_done_false_is_presence_guard_and_traversal_changes_actual_child() {
        let mut g = graph(2, &[(0, 1)]);
        let v = valence(&g);
        g.atoms[1]
            .set_prop("_TraversalBondIndexOrder", 9u32)
            .unwrap();
        g.atoms[1].set_prop("_CIPCode", "kept").unwrap();
        let mut props = MoleculeProperties::default();
        props.set_prop("_StereochemDone", false).unwrap();
        let mut r = rings(2, 1);
        let mut stack = vec![];
        canonicalize_fragment_source(
            &mut g,
            &mut props,
            &v,
            &mut r,
            0,
            &mut [AtomColor::White; 2],
            &[0, 1],
            &mut stack,
            None,
            None,
            None,
            false,
            None,
            true,
            &BTreeMap::new(),
        )
        .unwrap();
        assert_eq!(
            props.prop("_StereochemDone"),
            Some(&PropertyValue::Bool(false))
        );
        assert_eq!(
            g.atoms[0].prop("_TraversalStartPoint"),
            Some(&PropertyValue::Bool(true))
        );
        assert!(g.atoms[1].prop("_TraversalBondIndexOrder").is_none());
        assert_eq!(
            g.atoms[1].prop("_CIPCode"),
            Some(&PropertyValue::String("kept".into()))
        );
        assert!(r.is_symm_sssr());
    }
    #[test]
    fn source590_missing_stereo_done_commits_computed_marker_and_actual_sssr_cache() {
        let mut g = graph(2, &[(0, 1)]);
        let v = valence(&g);
        let mut props = MoleculeProperties::default();
        let mut r =
            cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::OtherOrUnknown, 0, 0);
        let mut stack = vec![];
        canonicalize_fragment_source(
            &mut g,
            &mut props,
            &v,
            &mut r,
            0,
            &mut [AtomColor::White; 2],
            &[0, 1],
            &mut stack,
            None,
            None,
            None,
            false,
            None,
            true,
            &BTreeMap::new(),
        )
        .unwrap();
        assert_eq!(props.prop("_StereochemDone"), Some(&PropertyValue::Int(1)));
        assert!(
            props
                .computed_prop_names()
                .unwrap()
                .unwrap()
                .iter()
                .any(|v| v.as_bytes() == b"_StereochemDone")
        );
        assert!(r.find_type() == cosmolkit_core::RingFindType::Sssr);
        assert!(r.is_initialized());
        assert_eq!(r.atom_row_count(), 0);
    }
    #[test]
    fn source590_ring_error_retains_reset_cache_before_startpoint_write() {
        let mut g = graph(2, &[(0, 1)]);
        let v = valence(&g);
        let mut props = MoleculeProperties::default();
        props.set_prop("_StereochemDone", 0i32).unwrap();
        props.set_prop("__computedProps", false).unwrap();
        let mut r = cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::Fast, 2, 1);
        let mut stack = vec![];
        assert!(matches!(
            canonicalize_fragment_source(
                &mut g,
                &mut props,
                &v,
                &mut r,
                0,
                &mut [AtomColor::White; 2],
                &[0; 2],
                &mut stack,
                None,
                None,
                None,
                false,
                None,
                true,
                &BTreeMap::new()
            ),
            Err(SmilesParseError::WriterRings(_))
        ));
        assert!(r.find_type() == cosmolkit_core::RingFindType::Sssr);
        assert_eq!(r.atom_row_count(), 0);
        assert!(g.atoms[0].prop("_TraversalStartPoint").is_none());
        assert!(stack.is_empty());
    }
    #[test]
    fn source590_ring_bad_any_cast_follows_actual_tag_change() {
        let mut g = graph(1, &[]);
        g.atoms[0].set_chiral_tag(ChiralTag::TetrahedralCw);
        g.atoms[0].set_prop("_ringStereoAtoms", false).unwrap();
        let v = valence(&g);
        assert!(matches!(
            phase(
                &mut g,
                &v,
                &rings(1, 0),
                &[MolStackElem::Atom(0)],
                None,
                true
            ),
            Err(SmilesParseError::PropertyKind(_))
        ));
        assert_eq!(g.atoms[0].chiral_tag(), ChiralTag::TetrahedralCcw);
    }
    #[test]
    fn source590_ring_neighbor_prefix_survives_later_invalid_relation() {
        let mut g = graph(2, &[]);
        g.atoms[0].set_chiral_tag(ChiralTag::TetrahedralCw);
        g.atoms[0]
            .set_prop("_ringStereoAtoms", PropertyValue::IntVector(vec![2, 0]))
            .unwrap();
        let v = valence(&g);
        assert!(matches!(
            phase(
                &mut g,
                &v,
                &rings(2, 0),
                &[MolStackElem::Atom(0), MolStackElem::Atom(1)],
                None,
                true
            ),
            Err(SmilesParseError::WriterCanonicalInvariant(_))
        ));
        assert_eq!(g.atoms[0].chiral_tag(), ChiralTag::TetrahedralCcw);
        assert_eq!(g.atoms[1].chiral_tag(), ChiralTag::TetrahedralCcw);
    }
    #[test]
    fn source590_stereo_group_inverts_unvisited_members_and_checks_source_size_t_tag() {
        let mut g = graph(2, &[]);
        g.atoms[0].set_chiral_tag(ChiralTag::TetrahedralCw);
        g.atoms[1].set_chiral_tag(ChiralTag::TetrahedralCw);
        g.atoms[0].set_prop("_stereoGroup", 0u32).unwrap();
        g.stereo_groups.push(
            StereoGroup::new(
                StereoGroupKind::And,
                vec![AtomId::new(0), AtomId::new(1)],
                vec![],
            )
            .expect("valid distinct stereo members"),
        );
        let v = valence(&g);
        phase(
            &mut g,
            &v,
            &rings(2, 0),
            &[MolStackElem::Atom(0)],
            None,
            true,
        )
        .unwrap();
        assert_eq!(g.atoms[0].chiral_tag(), ChiralTag::TetrahedralCcw);
        assert_eq!(g.atoms[1].chiral_tag(), ChiralTag::TetrahedralCcw);
        g.atoms[0].set_prop("_stereoGroup", false).unwrap();
        assert!(matches!(
            phase(
                &mut g,
                &v,
                &rings(2, 0),
                &[MolStackElem::Atom(0)],
                None,
                true
            ),
            Err(SmilesParseError::WriterULong(_))
        ));
        g.atoms[0]
            .set_prop("_stereoGroup", "18446744073709551615")
            .unwrap();
        phase(
            &mut g,
            &v,
            &rings(2, 0),
            &[MolStackElem::Atom(0)],
            None,
            true,
        )
        .unwrap();
        assert_eq!(g.atoms[1].chiral_tag(), ChiralTag::TetrahedralCcw);
    }
    #[test]
    fn source590_nonisomeric_skips_broken_chirality_and_permutation_property_getters() {
        let mut g = graph(2, &[(0, 1)]);
        g.atoms[0].set_chiral_tag(ChiralTag::SquarePlanar);
        g.atoms[0].set_prop("_chiralPermutation", false).unwrap();
        let v = valence(&g);
        phase(
            &mut g,
            &v,
            &rings(2, 1),
            &[MolStackElem::Atom(0)],
            Some(&[false]),
            false,
        )
        .unwrap();
        assert!(g.atoms[0].prop("_brokenChirality").is_none());
        phase(
            &mut g,
            &v,
            &rings(2, 1),
            &[MolStackElem::Atom(0)],
            Some(&[false]),
            true,
        )
        .unwrap();
        assert_eq!(
            g.atoms[0].prop("_brokenChirality"),
            Some(&PropertyValue::Bool(true))
        );
    }

    #[test]
    fn source590_chiral_inversion_flag_controls_first_explicit_hydrogen_case() {
        let mut original = graph(4, &[(0, 1), (0, 2), (0, 3)]);
        original.atoms[0].set_chiral_tag(ChiralTag::TetrahedralCw);
        original.atoms[0].set_explicit_hydrogens(1);
        let v = valence(&original);
        for invert in [false, true] {
            let mut g = original.clone();
            let stack = (0..4).map(MolStackElem::Atom).collect::<Vec<_>>();
            let mut orders = vec![vec![]; 4];
            orders[0] = vec![BondId::new(0), BondId::new(1), BondId::new(2)];
            canonicalize_source_stack(
                &mut g,
                &v,
                Some(&rings(4, 3)),
                &stack,
                &vec![vec![]; 4],
                &orders,
                &[false; 3],
                None,
                None,
                true,
                invert,
                &BTreeMap::new(),
                None,
            )
            .unwrap();
            assert_eq!(
                g.atoms[0].chiral_tag(),
                if invert {
                    ChiralTag::TetrahedralCcw
                } else {
                    ChiralTag::TetrahedralCw
                }
            );
        }
    }
    #[test]
    fn source590_non_tetrahedral_native_swap_table_writes_raw_int_and_emission_reads_it() {
        // Pinned swap_squareplanar_table[SP1][swap(0,1)] is SP3.
        let mut g = graph(5, &[(0, 1), (0, 2), (0, 3), (0, 4)]);
        g.atoms[0].set_chiral_tag(ChiralTag::SquarePlanar);
        g.atoms[0].set_chiral_permutation(Some(1));
        let v = valence(&g);
        let mut orders = vec![vec![]; 5];
        orders[0] = vec![
            BondId::new(1),
            BondId::new(0),
            BondId::new(2),
            BondId::new(3),
        ];
        canonicalize_source_stack(
            &mut g,
            &v,
            Some(&rings(5, 4)),
            &[MolStackElem::Atom(0)],
            &vec![vec![]; 5],
            &orders,
            &[false; 4],
            None,
            None,
            true,
            false,
            &BTreeMap::new(),
            None,
        )
        .unwrap();
        assert_eq!(
            g.atoms[0].prop("_chiralPermutation"),
            Some(&PropertyValue::Int(3))
        );
        assert_eq!(g.atoms[0].chiral_permutation(), Some(1));
        let text = atom_text(
            &g,
            &v,
            0,
            ChiralAdjustment::default(),
            &SmilesWriteParams::default(),
        )
        .unwrap();
        assert!(text.as_bytes().windows(4).any(|w| w == b"@SP3"));
    }
    #[test]
    fn source590_non_tetrahedral_getter_propagates_raw_numeric_errors_and_negative_zero_result() {
        let mut g = graph(1, &[]);
        g.atoms[0].set_chiral_tag(ChiralTag::SquarePlanar);
        g.atoms[0].set_prop("_chiralPermutation", false).unwrap();
        let v = valence(&g);
        assert!(matches!(
            phase(
                &mut g,
                &v,
                &rings(1, 0),
                &[MolStackElem::Atom(0)],
                None,
                true
            ),
            Err(SmilesParseError::WriterInt(_))
        ));
        g.atoms[0].set_prop("_chiralPermutation", -1i32).unwrap();
        phase(
            &mut g,
            &v,
            &rings(1, 0),
            &[MolStackElem::Atom(0)],
            None,
            true,
        )
        .unwrap();
        assert_eq!(
            g.atoms[0].prop("_chiralPermutation"),
            Some(&PropertyValue::Int(-1))
        );
    }
}

#[cfg(test)]
mod source594_mask_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, BondSpec, Element, MoleculeProperties};
    fn graph() -> TopologyBlock {
        TopologyBlock::try_from_parts(
            (0..4)
                .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
                .collect(),
            [(0, 1), (1, 2)]
                .into_iter()
                .enumerate()
                .map(|(i, (a, b))| {
                    Bond::from_spec(
                        BondId::new(i),
                        BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                    )
                })
                .collect(),
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn run(g: &mut TopologyBlock, mask: Option<&[bool]>) -> Result<(), SmilesParseError> {
        let v = cosmolkit_core::assign_valence_with_options_for_topology(
            g,
            ValenceModel::RdkitLike,
            false,
        )
        .unwrap();
        let mut p = MoleculeProperties::default();
        p.set_prop("_StereochemDone", false).unwrap();
        let mut rings = cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::SymmSssr, 4, 2);
        canonicalize_fragment_from_bond_mask_source(
            g,
            &mut p,
            &v,
            &mut rings,
            3,
            &mut [AtomColor::Black; 4],
            &[0; 4],
            &mut vec![],
            mask,
            None,
            true,
            None,
            true,
            &BTreeMap::new(),
        )
    }
    #[test]
    fn source594_mask_selects_bond_endpoints_including_unvisited_black_atoms() {
        let mut g = graph();
        g.atoms[1].set_chiral_tag(ChiralTag::TetrahedralCw);
        g.atoms[3].set_chiral_tag(ChiralTag::SquarePlanar);
        g.atoms[3].set_prop("_chiralPermutation", false).unwrap();
        run(&mut g, Some(&[false, true])).unwrap();
        assert_eq!(
            g.atoms[1].prop("_brokenChirality"),
            Some(&PropertyValue::Bool(true))
        );
        assert!(g.atoms[3].prop("_brokenChirality").is_none());
        assert_eq!(
            g.atoms[3].prop("_chiralPermutation"),
            Some(&PropertyValue::Bool(false))
        );
    }
    #[test]
    fn source594_absent_mask_includes_isolated_atoms_but_explicit_empty_selection_does_not() {
        let mut g = graph();
        g.atoms[3].set_chiral_tag(ChiralTag::SquarePlanar);
        g.atoms[3].set_prop("_chiralPermutation", false).unwrap();
        let mut selected = g.clone();
        assert!(matches!(
            run(&mut selected, None),
            Err(SmilesParseError::WriterInt(_))
        ));
        assert_eq!(
            selected.atoms[3].prop("_TraversalStartPoint"),
            Some(&PropertyValue::Bool(true))
        );
        run(&mut g, Some(&[false, false])).unwrap();
        assert!(g.atoms.iter().all(|a| a.prop("_brokenChirality").is_none()));
    }
    #[test]
    fn source594_mask_bounds_are_reached_before_main_preconditions_or_state_changes() {
        let mut g = graph();
        let before = g.clone();
        let v = cosmolkit_core::assign_valence_with_options_for_topology(
            &g,
            ValenceModel::RdkitLike,
            false,
        )
        .unwrap();
        let mut p = MoleculeProperties::default();
        let oldp = p.clone();
        let mut r = cosmolkit_core::RingInfo::new(cosmolkit_core::RingFindType::SymmSssr, 4, 2);
        let oldr = r.clone();
        let error = canonicalize_fragment_from_bond_mask_source(
            &mut g,
            &mut p,
            &v,
            &mut r,
            3,
            &mut [],
            &[],
            &mut vec![],
            Some(&[false]),
            None,
            true,
            None,
            true,
            &BTreeMap::new(),
        )
        .unwrap_err();
        assert!(matches!(
            error,
            SmilesParseError::TraversalStateIndex {
                state: "bondsInPlay",
                index: 1,
                count: 1
            }
        ));
        assert_eq!(g, before);
        assert_eq!(p, oldp);
        assert_eq!(r, oldr);
    }
}

fn source_query_has_single_h(
    query: &cosmolkit_model::QueryNode<cosmolkit_model::AtomQueryPredicate>,
) -> bool {
    // RDKit❗❌: bool hasSingleHQuery(const Atom::QUERYATOM_QUERY *q) {
    // RDKit❗❌:   // list queries are series of nested ors of AtomAtomicNum queries
    // RDKit❗❌:   PRECONDITION(q, "bad query");
    // RDKit❗❌:   bool res = false;
    // RDKit❗❌:   const auto &descr = q->getDescription();
    // RDKit❗❌:   if (descr == "AtomAnd") {
    // RDKit❗❌:     for (auto cIt = q->beginChildren(); cIt != q->endChildren(); ++cIt) {
    // RDKit❗❌:       const auto &cDescr = (*cIt)->getDescription();
    // RDKit❗❌:       if (cDescr == "AtomHCount") {
    // RDKit❗❌:         return !(*cIt)->getNegation() &&
    // RDKit❗❌:                ((ATOM_EQUALS_QUERY *)(*cIt).get())->getVal() == 1;
    // RDKit❗❌:       } else if (cDescr == "AtomAnd") {
    // RDKit❗❌:         res = hasSingleHQuery((*cIt).get());
    // RDKit❗❌:         if (res) {
    // RDKit❗❌:           return true;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    use cosmolkit_model::{AtomQueryPredicate, QueryNode};
    fn own_flag(mut q: &QueryNode<AtomQueryPredicate>) -> (&QueryNode<AtomQueryPredicate>, bool) {
        let mut n = false;
        while let QueryNode::Not(c) = q {
            n = !n;
            q = c;
        }
        (q, n)
    }
    let (query, _) = own_flag(query);
    if let QueryNode::And(children) = query {
        for child in children {
            let (leaf, negated) = own_flag(child);
            match leaf {
                QueryNode::Predicate(AtomQueryPredicate::HydrogenCount(n)) => {
                    return !negated && *n == 1;
                }
                QueryNode::And(_) => {
                    if source_query_has_single_h(leaf) {
                        return true;
                    }
                }
                _ => {}
            }
        }
    }
    false
}

/// Query-aware carrier projection of the sole full Canon pipeline. Query trees
/// and source hasQuery flags remain authoritative borrowed QueryStateRef data;
/// no predicate or ordinary-atom replacement is inferred from chemical values.
#[doc(hidden)]
#[allow(clippy::too_many_arguments)]
pub fn canonicalize_query_fragment_source(
    query: &mut cosmolkit_model::QueryGraph,
    rings: &mut cosmolkit_core::RingInfo,
    atom: usize,
    colors: &mut [AtomColor],
    ranks: &[u32],
    stack: &mut Vec<MolStackElem>,
    atoms_in_play: Option<&[bool]>,
    bonds_in_play: Option<&[bool]>,
    isomeric_smiles: bool,
) -> Result<(), SmilesParseError> {
    // Source cache update already ran on every actual carrier before this
    // boundary, including non-element identities. A failed table lookup never
    // reaches the Element-only structural projection or substitutes a dummy.
    // Extra common-property and bond copies are a known performance ❌. Query
    // trees themselves are borrowed, not cloned; all completed effects return
    // on errors too, as required by Native mutable ROMol prefix semantics.
    let atoms = query
        .atoms()
        .iter()
        .map(|a| {
            a.try_to_atom().map_err(|_| {
                SmilesParseError::WriterCanonicalInvariant(
                    "query cache identity was not representable after source cache update",
                )
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let bonds = query
        .bonds()
        .iter()
        .map(|b| b.bond().clone())
        .collect::<Vec<_>>();
    let adjacency = cosmolkit_model::AdjacencyList::from_topology(atoms.len(), &bonds);
    let mut topology = TopologyBlock {
        atoms,
        bonds,
        adjacency,
        substance_groups: vec![],
        stereo_groups: query.stereo_groups().to_vec(),
    };
    let valence = ValenceAssignment {
        explicit_valence: topology
            .atoms
            .iter()
            .map(|a| i32::from(a.source_valence_facts().explicit_valence))
            .collect(),
        implicit_hydrogens: topology
            .atoms
            .iter()
            .map(|a| {
                if a.no_implicit() {
                    0
                } else {
                    i32::from(a.source_valence_facts().implicit_valence)
                }
            })
            .collect(),
    };
    let mut properties = query.source_molecule_properties();
    let source_queries =
        cosmolkit_model::QueryStateRef::try_for_topology(query.atoms(), query.bonds(), &topology)
            .map_err(|_| {
            SmilesParseError::WriterCanonicalInvariant("source query/carrier row alignment")
        })?;
    let result = canonicalize_fragment_impl(
        &mut topology,
        &mut properties,
        &valence,
        rings,
        atom,
        colors,
        ranks,
        stack,
        atoms_in_play,
        bonds_in_play,
        None,
        isomeric_smiles,
        None,
        true,
        &BTreeMap::new(),
        Some(source_queries),
    );
    for (out, carrier) in query.atoms_mut().iter_mut().zip(&topology.atoms) {
        out.replace_source_carrier_members_from(carrier);
    }
    for (out, carrier) in query.bonds_mut().iter_mut().zip(topology.bonds) {
        *out.bond_mut() = carrier;
    }
    query.replace_source_molecule_properties(&properties);
    result
}

#[cfg(test)]
mod source1110_query_h_tests {
    use super::*;
    use cosmolkit_model::{AtomQueryPredicate, QueryNode};
    fn h(n: i32) -> QueryNode<AtomQueryPredicate> {
        QueryNode::predicate(AtomQueryPredicate::HydrogenCount(n))
    }
    #[test]
    fn bare_h_leaf_and_or_are_not_the_source_and_description() {
        assert!(!source_query_has_single_h(&h(1)));
        assert!(!source_query_has_single_h(&QueryNode::or(vec![h(1), h(0)])));
    }
    #[test]
    fn direct_child_reads_own_negation_and_exact_integer() {
        for n in [-1, 0, 1, 2, i32::MAX] {
            assert_eq!(
                source_query_has_single_h(&QueryNode::and(vec![h(n)])),
                n == 1
            );
            assert!(!source_query_has_single_h(&QueryNode::and(vec![
                QueryNode::not(h(n))
            ])));
        }
    }
    #[test]
    fn first_h_child_returns_immediately_even_when_later_child_is_one() {
        assert!(!source_query_has_single_h(&QueryNode::and(vec![
            h(0),
            h(1)
        ])));
        assert!(source_query_has_single_h(&QueryNode::and(vec![h(1), h(0)])));
    }
    #[test]
    fn nested_and_false_result_continues_to_later_siblings() {
        assert!(source_query_has_single_h(&QueryNode::and(vec![
            QueryNode::and(vec![h(0)]),
            h(1)
        ])));
        assert!(source_query_has_single_h(&QueryNode::and(vec![
            QueryNode::and(vec![h(1)]),
            h(0)
        ])));
    }
    #[test]
    fn composite_own_negation_does_not_propagate_to_h_child() {
        assert!(source_query_has_single_h(&QueryNode::not(QueryNode::and(
            vec![h(1)]
        ))));
        assert!(source_query_has_single_h(&QueryNode::and(vec![
            QueryNode::not(QueryNode::and(vec![h(1)]))
        ])));
    }
    #[test]
    fn unrelated_and_or_branches_are_ignored_and_input_is_immutable() {
        let q = QueryNode::and(vec![
            QueryNode::or(vec![h(1), h(1)]),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            QueryNode::and(vec![h(0)]),
        ]);
        let before = q.clone();
        assert!(!source_query_has_single_h(&q));
        assert_eq!(q, before);
    }
}

/// Known causes from source SMARTS traversal preparation.
#[derive(Debug, Clone)]
pub enum SmartsTraversalError {
    Valence(cosmolkit_core::ValenceError),
    Stereo(std::sync::Arc<cosmolkit_core::LegacyStereoError>),
    Traversal(SmilesParseError),
}
impl std::fmt::Display for SmartsTraversalError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Valence(error) => write!(f, "SMARTS traversal valence preparation: {error}"),
            Self::Stereo(error) => write!(f, "SMARTS traversal stereo preparation: {error}"),
            Self::Traversal(error) => std::fmt::Display::fmt(error, f),
        }
    }
}
impl std::error::Error for SmartsTraversalError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(match self {
            Self::Valence(error) => error,
            Self::Stereo(error) => error.as_ref(),
            Self::Traversal(error) => error,
        })
    }
}
impl PartialEq for SmartsTraversalError {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Valence(a), Self::Valence(b)) => a == b,
            (Self::Stereo(a), Self::Stereo(b)) => std::sync::Arc::ptr_eq(a, b),
            (Self::Traversal(a), Self::Traversal(b)) => a == b,
            _ => false,
        }
    }
}
impl Eq for SmartsTraversalError {}
impl From<SmilesParseError> for SmartsTraversalError {
    fn from(error: SmilesParseError) -> Self {
        Self::Traversal(error)
    }
}

/// Reuse the Canon traversal and stereo owners for concrete SMARTS output.
/// The input is detached; no live molecule or runtime cache is accessible.
#[doc(hidden)]
pub fn prepare_smarts_serialization_topology(
    mut topology: TopologyBlock,
    properties: &cosmolkit_model::MoleculeProperties,
    rooted_at_atom: Option<usize>,
    do_isomeric: bool,
) -> Result<TopologyBlock, SmartsTraversalError> {
    // RDKit❗❌:   ROMol mol(inmol);
    // RDKit❗❌:   mol.getRingInfo()->reset();
    // RDKit❗❌:   mol.getRingInfo()->initialize(FIND_RING_TYPE_SYMM_SSSR);
    // RDKit❗❌:   for (auto &atom : mol.atoms()) {
    // RDKit❗❌:     atom->updatePropertyCache(false);
    // RDKit❗❌:   }
    // RDKit❗❌:   bool doRandom = false;
    // RDKit❗❌:   bool doChiralInversions = true;
    // RDKit❗❌:   Canon::canonicalizeFragment(
    // RDKit❗❌:       mol, atomIdx, colors, ranks, molStack, &atomsInPlay, bondsInPlay, nullptr,
    // RDKit❗❌:       params.doIsomericSmiles, doRandom, doChiralInversions);
    // This adapter selects the source SMARTS profile of the existing traversal,
    // permutation, relative-stereo and double-bond direction owners. SEARCH
    // retains its serializer; no SMARTS or stereo algorithm is duplicated here.
    // Cost review: traversal is linear plus existing neighbor sorts. SEARCH
    // currently traverses the prepared carrier again to emit tokens, an extra
    // O(V+E) pass compared with SOURCE's single shared stack, hence ❌ complexity.
    let count = topology.atoms.len();
    if count == 0 {
        return Ok(topology);
    }
    let valence = cosmolkit_core::assign_valence_with_options_for_topology(
        &topology,
        ValenceModel::RdkitLike,
        false,
    )
    .map_err(SmartsTraversalError::Valence)?;
    let mut rings = cosmolkit_core::RingInfo::new(
        cosmolkit_core::RingFindType::SymmSssr,
        count,
        topology.bonds.len(),
    );
    // RDKit❗✔️:   if (!mol.hasProp(common_properties::_StereochemDone)) {
    // RDKit❗✔️:     MolOps::assignStereochemistry(mol, false);
    // RDKit❗✔️:   }
    if properties.prop("_StereochemDone").is_none() {
        topology = cosmolkit_core::assign_legacy_stereochemistry_with_flags(
            topology, &valence, &rings, false, false,
        )
        .map_err(|error| SmartsTraversalError::Stereo(std::sync::Arc::new(error)))?;
    }
    // RDKit✔️✔️:   for (const auto &atom : mol.atoms()) {
    // RDKit✔️✔️:     ranks.push_back(atom->getIdx());
    // RDKit✔️✔️:   }
    let ranks = (0..count).map(|index| index as u32).collect::<Vec<_>>();
    let mut working_properties = properties.clone();
    let mut colors = vec![AtomColor::White; count];
    let mut root = rooted_at_atom;
    while colors.contains(&AtomColor::White) {
        // RDKit✔️✔️:       // Try to find a non-chiral atom we have not processed yet.
        // RDKit✔️✔️:       // If we can't find non-chiral atom, use the chiral atom with
        // RDKit✔️✔️:       // the lowest rank (we are guaranteed to find an unprocessed atom).
        let start = root.take().unwrap_or_else(|| {
            (0..count)
                .find(|index| {
                    colors[*index] == AtomColor::White
                        && !matches!(
                            topology.atoms[*index].chiral_tag(),
                            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw,
                        )
                })
                .unwrap_or_else(|| {
                    (0..count)
                        .find(|index| colors[*index] == AtomColor::White)
                        .expect("white atom exists")
                })
        });
        if start >= count {
            return Err(SmilesParseError::Model("bad atom index".into()).into());
        }
        let mut stack = Vec::with_capacity(count + topology.bonds.len());
        canonicalize_fragment_source(
            &mut topology,
            &mut working_properties,
            &valence,
            &mut rings,
            start,
            &mut colors,
            &ranks,
            &mut stack,
            None,
            None,
            None,
            do_isomeric,
            None,
            true,
            &BTreeMap::new(),
        )?;
    }
    Ok(topology)
}

#[cfg(test)]
mod recovery_chem09 {
    use super::*;

    fn chem27_controlled_graph(
        n: usize,
        edges: &[(usize, usize, BondOrder)],
        target: usize,
        stereo: BondStereo,
        pair: Option<[usize; 2]>,
    ) -> TopologyBlock {
        let atoms = (0..n)
            .map(|i| {
                Atom::from_spec(
                    AtomId::new(i),
                    cosmolkit_model::AtomSpec::new(cosmolkit_types::Element::C),
                )
            })
            .collect();
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(i, &(a, b, o))| {
                Bond::from_spec(
                    BondId::new(i),
                    cosmolkit_model::BondSpec::new(AtomId::new(a), AtomId::new(b), o),
                )
            })
            .collect();
        let mut t = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap();
        t.bonds[target].set_stereo_atoms(pair.map(|p| p.map(AtomId::new)));
        t.bonds[target].set_stereo(stereo).unwrap();
        t
    }
    fn chem27_assert_target(
        t: &TopologyBlock,
        b: usize,
        stereo: BondStereo,
        pair: Option<[usize; 2]>,
    ) {
        assert_eq!(t.bonds[b].stereo(), stereo);
        assert_eq!(t.bonds[b].stereo_atoms(), pair.map(|p| p.map(AtomId::new)));
    }
    fn chem09_chain() -> TopologyBlock {
        let mut m = chem27_controlled_graph(
            6,
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Double),
                (2, 3, BondOrder::Single),
                (3, 4, BondOrder::Double),
                (4, 5, BondOrder::Single),
            ],
            1,
            BondStereo::E,
            Some([0, 3]),
        );
        let b = &mut m.bonds[3];
        b.set_stereo(BondStereo::E).unwrap();
        b.set_stereo_atoms(Some([AtomId::new(2), AtomId::new(5)]));
        assert!(m.bonds.iter().all(|b| b.direction() == BondDirection::None));
        m
    }

    fn chem09_stack_through(last_atom: usize) -> Vec<MolStackElem> {
        let mut stack = vec![MolStackElem::Atom(0)];
        for i in 0..last_atom {
            stack.push(MolStackElem::Bond {
                bond: BondId::new(i),
                atom_to_left: i,
            });
            stack.push(MolStackElem::Atom(i + 1));
        }
        stack
    }

    fn chem09_run_queue(
        m: &mut TopologyBlock,
        stack: &[MolStackElem],
        atom_orders: &[usize],
        bond_orders: &[usize],
    ) -> (Vec<i8>, Vec<i8>) {
        // Direct queue entry has the actual outer caller's precleared directions.
        assert!(m.bonds.iter().all(|b| b.direction() == BondDirection::None));
        let mut bc = vec![0; m.bonds.len()];
        let mut ac = vec![0; m.atoms.len()];
        let rings = vec![false; m.bonds.len()];
        direction::canonicalize_double_bonds_for_writer(
            m,
            bond_orders,
            atom_orders,
            &rings,
            &mut bc,
            &mut ac,
            stack,
        )
        .unwrap();
        (bc, ac)
    }

    #[test]
    fn rdkit_2026_03_6_chem09_first_position_one_is_live_and_each_missing_controller_cleans_pair() {
        use direction::chem09_trace::{Event, capture};
        let stack = chem09_stack_through(3);
        for missing in [None, Some(0), Some(3)] {
            let mut m = chem09_chain();
            let mut visits = [1, 3, 5, 7, 0, 0];
            if let Some(a) = missing {
                visits[a] = 0;
            }
            let trace = capture();
            let (bc, ac) = chem09_run_queue(&mut m, &stack, &visits, &[2, 4, 6, 0, 0]);
            let events = trace.events();
            if missing.is_some() {
                chem27_assert_target(&m, 1, BondStereo::None, None);
                assert!(events.contains(&Event::InitialCleanup(BondId::new(1))));
                assert!(!events.contains(&Event::Canonicalize(BondId::new(1))));
                assert!(bc.iter().all(|&x| x == 0));
                assert!(ac.iter().all(|&x| x == 0));
            } else {
                chem27_assert_target(&m, 1, BondStereo::E, Some([0, 3]));
                assert!(events.contains(&Event::Canonicalize(BondId::new(1))));
                assert!(bc[0] > 0 && bc[2] > 0);
                assert!(ac[1] > 0 && ac[2] > 0);
                assert_ne!(m.bonds[0].direction(), BondDirection::None);
            }
            // Missing-controller arrays are controlled predicate masks, not natural
            // traversal claims; the first-position-one control is a consistent real stack.
        }
    }

    #[test]
    fn rdkit_2026_03_6_chem09_omitted_neighbor_is_enqueued_then_skipped_without_seen() {
        use direction::chem09_trace::{Event, capture};
        let mut m = chem09_chain();
        let trace = capture();
        let (bc, _) = chem09_run_queue(
            &mut m,
            &chem09_stack_through(3),
            &[1, 3, 5, 7, 0, 0],
            &[2, 4, 6, 0, 0],
        );
        let events = trace.events();
        let a = BondId::new(1);
        let b = BondId::new(3);
        let expected = [
            Event::Neighbors(a, vec![b]),
            Event::Canonicalize(a),
            Event::SeenMarked(a),
            Event::Enqueue(a, b),
            Event::Dequeue(b, false, true),
            Event::SkipOmitted(b),
        ];
        let mut previous = None;
        for event in expected {
            let pos = events
                .iter()
                .position(|e| *e == event)
                .unwrap_or_else(|| panic!("missing {event:?}: {events:?}"));
            if let Some(prev) = previous {
                assert!(prev < pos);
            }
            previous = Some(pos);
        }
        assert!(!events.contains(&Event::SeenMarked(b)));
        assert!(!events.contains(&Event::InvalidCleanup(b)));
        assert!(
            !events
                .iter()
                .any(|e| matches!(e,Event::Enqueue(from,_) if *from==b))
        );
        chem27_assert_target(&m, 3, BondStereo::E, Some([2, 5]));
        assert_eq!(bc[4], 0);
        assert_eq!(m.bonds[4].direction(), BondDirection::None);
        assert!(bc[2] > 0); // Shared edge belongs to valid A, so it must be allowed to change.
        drop(trace);
        let fresh = capture();
        assert!(fresh.events().is_empty());
    }

    #[test]
    fn rdkit_2026_03_6_chem09_initial_invalid_tags_clear_all_pairs_and_outer_slashes() {
        // Injected non-double/atrop/pairless cases isolate the source predicate;
        // they are not claims that parsing naturally produces these stereo states.
        for (order, stereo, pair) in [
            (BondOrder::Double, BondStereo::Any, Some([0, 3])),
            (BondOrder::Single, BondStereo::E, Some([0, 3])),
            (BondOrder::Dative, BondStereo::E, Some([0, 3])),
            (BondOrder::Double, BondStereo::AtropCw, Some([0, 3])),
            (BondOrder::Double, BondStereo::AtropCcw, Some([0, 3])),
            (BondOrder::Double, BondStereo::E, None),
            (BondOrder::Double, BondStereo::None, Some([0, 3])),
        ] {
            let mut m = chem27_controlled_graph(
                4,
                &[
                    (0, 1, BondOrder::Single),
                    (1, 2, order),
                    (2, 3, BondOrder::Single),
                ],
                1,
                stereo,
                pair,
            );
            for b in &mut m.bonds {
                b.set_direction(BondDirection::EndUpRight);
            }
            direction::canonicalize_double_bond_directions_for_writer(
                &mut m,
                &chem09_stack_through(3),
                &[false; 3],
            )
            .unwrap();
            chem27_assert_target(&m, 1, BondStereo::None, None);
            assert!(m.bonds.iter().all(|b| b.direction() == BondDirection::None));
        }
    }

    #[test]
    fn rdkit_2026_03_6_chem09_live_invalid_bfs_cleanup_marks_seen_without_propagation() {
        use direction::chem09_trace::{Event, capture};
        let mut m = chem09_chain();
        let trace = capture();
        let (bc, _) = chem09_run_queue(
            &mut m,
            &chem09_stack_through(4),
            &[1, 3, 5, 7, 9, 0],
            &[2, 4, 6, 8, 0],
        );
        let events = trace.events();
        let a = BondId::new(1);
        let b = BondId::new(3);
        let expected = [
            Event::Neighbors(a, vec![b]),
            Event::InitialCleanup(b),
            Event::Canonicalize(a),
            Event::Enqueue(a, b),
            Event::Dequeue(b, false, false),
            Event::InvalidCleanup(b),
            Event::SeenMarked(b),
        ];
        let mut previous = None;
        for event in expected {
            let pos = events
                .iter()
                .position(|e| *e == event)
                .unwrap_or_else(|| panic!("missing {event:?}: {events:?}"));
            if let Some(prev) = previous {
                assert!(prev < pos);
            }
            previous = Some(pos);
        }
        assert!(!events.contains(&Event::Canonicalize(b)));
        assert!(!events.contains(&Event::SkipOmitted(b)));
        assert!(
            !events
                .iter()
                .any(|e| matches!(e,Event::Enqueue(from,_) if *from==b))
        );
        chem27_assert_target(&m, 3, BondStereo::None, None);
        assert_eq!(bc[4], 0);
    }

    #[test]
    fn rdkit_2026_03_6_chem09_duplicate_dequeue_skips_seen_before_reprocessing() {
        use direction::chem09_trace::{Event, capture};
        let mut m = chem27_controlled_graph(
            8,
            &[
                (0, 1, BondOrder::Double),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Double),
                (3, 4, BondOrder::Single),
                (4, 5, BondOrder::Double),
                (5, 6, BondOrder::Single),
                (6, 7, BondOrder::Double),
                (7, 0, BondOrder::Single),
            ],
            0,
            BondStereo::E,
            Some([7, 2]),
        );
        for (bond, pair) in [(2, [1, 4]), (4, [3, 6]), (6, [5, 0])] {
            let b = &mut m.bonds[bond];
            b.set_stereo(BondStereo::E).unwrap();
            b.set_stereo_atoms(Some(pair.map(AtomId::new)));
        }
        // Controlled complete ring graph gives two queued paths to C before C is popped.
        // All initial directions are clear; no natural ring-stereo assignment claimed.
        let mut stack = chem09_stack_through(7);
        stack.push(MolStackElem::Bond {
            bond: BondId::new(7),
            atom_to_left: 7,
        });
        let trace = capture();
        let _ = chem09_run_queue(
            &mut m,
            &stack,
            &[1, 3, 5, 7, 9, 11, 13, 15],
            &[2, 4, 6, 8, 10, 12, 14, 16],
        );
        let events = trace.events();
        let c = BondId::new(4);
        assert_eq!(
            events
                .iter()
                .filter(|e| **e == Event::Canonicalize(c))
                .count(),
            1
        );
        let first = events
            .iter()
            .position(|e| *e == Event::Dequeue(c, false, false))
            .unwrap();
        let seen = events
            .iter()
            .position(|e| *e == Event::SeenMarked(c))
            .unwrap();
        let second = events
            .iter()
            .position(|e| *e == Event::Dequeue(c, true, false))
            .unwrap();
        let skip = events
            .iter()
            .position(|e| *e == Event::SkipSeen(c))
            .unwrap();
        assert!(first < seen && seen < second && second < skip);
        assert!(
            events
                .iter()
                .any(|e| matches!(e, Event::PrioritySkipSeen(_)))
        );
    }

    #[test]
    fn rdkit_2026_03_6_chem09_source_stack_first_controller_and_unvisited_ring_relative() {
        use direction::chem09_trace::{Event, capture};
        let mut g = chem09_chain();
        let v = cosmolkit_core::assign_valence_with_options_for_topology(
            &g,
            ValenceModel::RdkitLike,
            false,
        )
        .unwrap();
        let stack = chem09_stack_through(3);
        let trace = capture();
        canonicalize_source_stack(
            &mut g,
            &v,
            None,
            &stack,
            &vec![vec![]; 6],
            &vec![vec![]; 6],
            &[false; 5],
            None,
            None,
            false,
            false,
            &BTreeMap::new(),
            None,
        )
        .unwrap();
        assert!(
            trace
                .events()
                .contains(&Event::Canonicalize(BondId::new(1)))
        );
        assert!(trace.events().contains(&Event::SkipOmitted(BondId::new(3))));
        chem27_assert_target(&g, 1, BondStereo::E, Some([0, 3]));
        drop(trace);
        // The same source0/positive arrays feed the ring-relative comparison:
        // later visited member changes; unvisited member retains its prior tag.
        for visited in [false, true] {
            let atoms = (0..2)
                .map(|i| {
                    Atom::from_spec(
                        AtomId::new(i),
                        cosmolkit_model::AtomSpec::new(cosmolkit_types::Element::C),
                    )
                })
                .collect();
            let mut g = TopologyBlock::try_from_parts(atoms, vec![], vec![], vec![]).unwrap();
            g.atoms[0].set_chiral_tag(ChiralTag::TetrahedralCw);
            g.atoms[1].set_chiral_tag(ChiralTag::TetrahedralCw);
            g.atoms[0]
                .set_prop("_ringStereoAtoms", PropertyValue::IntVector(vec![2]))
                .unwrap();
            let v = cosmolkit_core::assign_valence_with_options_for_topology(
                &g,
                ValenceModel::RdkitLike,
                false,
            )
            .unwrap();
            let stack = if visited {
                vec![MolStackElem::Atom(0), MolStackElem::Atom(1)]
            } else {
                vec![MolStackElem::Atom(0)]
            };
            canonicalize_source_stack(
                &mut g,
                &v,
                None,
                &stack,
                &[vec![], vec![]],
                &[vec![], vec![]],
                &[],
                None,
                None,
                true,
                false,
                &BTreeMap::new(),
                None,
            )
            .unwrap();
            assert_eq!(g.atoms[0].chiral_tag(), ChiralTag::TetrahedralCcw);
            assert_eq!(
                g.atoms[1].chiral_tag(),
                if visited {
                    ChiralTag::TetrahedralCcw
                } else {
                    ChiralTag::TetrahedralCw
                }
            );
        }
    }
}

#[cfg(test)]
mod recovery_chem27 {
    use super::*;
    #[test]
    fn official_9368_ten_cases_plain_borrowed_context_and_cx_share_prepass() {
        for (text, selection, expected) in [
            ("[*:1]/C=C/C=C/c1ccc(OC)cc1", vec![0, 1, 2], "C=C[*:1]"),
            (
                "[*:1]/C=C/C=C/c1ccc(OC)cc1",
                vec![0, 1, 2, 3],
                "C/C=C/[*:1]",
            ),
            ("C/C=C/C", vec![0, 1, 2], "C=CC"),
            ("C/C=C/C", vec![1, 2, 3], "C=CC"),
            ("C/C=C(F)/C", vec![0, 1, 2, 3], "C/C=C\\F"),
            ("C/C=C(F)/C", vec![0, 1, 2, 4], "C/C=C/C"),
            ("C/C(F)=C/C", vec![0, 1, 3, 4], "C\\C=C\\C"),
            ("C/C(F)=C/C", vec![2, 1, 3, 4], "C/C=C\\F"),
            ("C/C(F)=C/C", vec![0, 1, 3], "C=CC"),
            ("C/C(F)=C/C", vec![2, 1, 3], "C=CF"),
        ] {
            let r = crate::parse_smiles_complete_source(text, &Default::default()).unwrap();
            let before = r.clone();
            let atoms = selection.into_iter().map(AtomId::new).collect::<Vec<_>>();
            let params = SmilesWriteParams::default();
            let plain =
                write_fragment_smiles_output(&r, &params, &atoms, None, None, None, None, None)
                    .unwrap();
            assert_eq!(plain.text.as_bytes(), expected.as_bytes(), "plain {text}");
            let view = crate::SmilesRecordView::from(&r);
            let context =
                write_fragment_smiles_output(view, &params, &atoms, None, None, None, None, None)
                    .unwrap();
            assert_eq!(context, plain, "borrowed {text}");
            let cxparams = CxSmilesWriteParams {
                fields: crate::CxSmilesFields::NONE,
                ..Default::default()
            };
            let cx = write_fragment_cx_smiles(&r, &cxparams, &atoms, None, None, None, None, None)
                .unwrap();
            assert_eq!(cx.as_bytes(), expected.as_bytes(), "CX {text}");
            assert_eq!(r, before);
        }
    }
}
