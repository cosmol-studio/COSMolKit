use std::borrow::Cow;

use cosmolkit_core::{
    CanonicalRankError, CanonicalRankParams, LegacyStereoError, RingFindType, RingFindingError,
    RingInfo, ValenceAssignment, ValenceError, ValenceModel,
};
use cosmolkit_model::{AtomId, AtomPropertyError, BondId, BondValueError, TopologyBlock};

use crate::SmilesRecord;
use crate::writer::SmilesWriteParams;

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum FragmentWriteInputError {
    #[error("no atoms provided")]
    NoAtomsProvided,
    #[error("root atom index {atom_index} must be less than the molecule atom count {atom_count}")]
    RootAtomOutOfRange {
        atom_index: usize,
        atom_count: usize,
    },
    #[error("root atom index {atom_index} was not found in the selected atoms")]
    RootAtomNotSelected { atom_index: usize },
    #[error("atom symbol vector length {len} is less than the molecule atom count {expected}")]
    AtomSymbolsTooShort { len: usize, expected: usize },
    #[error("bond symbol vector length {len} is less than the molecule bond count {expected}")]
    BondSymbolsTooShort { len: usize, expected: usize },
    #[error("selected atom index {atom_index} is out of range for {atom_count} atoms")]
    AtomOutOfRange {
        atom_index: usize,
        atom_count: usize,
    },
    #[error("selected bond index {bond_index} is out of range for {bond_count} bonds")]
    BondOutOfRange {
        bond_index: usize,
        bond_count: usize,
    },
    #[error(
        "root atom index {atom_index} requires a single-fragment molecule, found {fragment_count} fragments"
    )]
    RootAtomRequiresSingleFragment {
        atom_index: usize,
        fragment_count: usize,
    },
    #[error("invalid fragment topology: {0}")]
    InvalidTopology(String),
    #[error(
        "existing valence rows have lengths explicit={explicit_len}, implicit={implicit_len}; expected {atom_count} atoms"
    )]
    ValenceAssignmentLengthMismatch {
        atom_count: usize,
        explicit_len: usize,
        implicit_len: usize,
    },
    #[error("fragment valence preparation failed: {0}")]
    Valence(#[from] ValenceError),
    #[error("fragment ring preparation failed: {0}")]
    Ring(#[from] RingFindingError),
    #[error("fragment legacy stereochemistry failed: {0}")]
    Stereo(#[from] LegacyStereoError),
    #[error("fragment canonical ranking failed: {0}")]
    Rank(#[from] CanonicalRankError),
    #[error("fragment SMILES emission failed: {0}")]
    Writer(#[from] crate::SmilesParseError),
    #[error("fragment CIP property conversion failed: {0}")]
    CipProperty(#[from] cosmolkit_core::PropertyStringError),
    #[error("fragment CIP property assignment failed: {0}")]
    AtomProperty(#[from] AtomPropertyError),
    #[error("fragment enhanced-stereo bond inversion failed: {0}")]
    StereoBond(#[from] BondValueError),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct FragmentSelectionMasks {
    pub(super) atoms_in_play: Vec<bool>,
    pub(super) bonds_in_play: Vec<bool>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct FragmentRingRows {
    pub(super) initialized: bool,
    pub(super) atom_rings: Vec<Vec<AtomId>>,
    pub(super) bond_rings: Vec<Vec<BondId>>,
}

pub(super) struct PreparedFragmentStereo<'a> {
    pub(super) topology: TopologyBlock,
    pub(super) retained_rings: Option<RingInfo>,
    pub(super) ranking_rings: Option<RingInfo>,
    pub(super) valence: Cow<'a, ValenceAssignment>,
    pub(super) stereochem_done: bool,
}

pub(super) fn prepare_fragment_stereo<'record, 'a>(
    record: impl Into<crate::SmilesRecordView<'record>>,
    params: &SmilesWriteParams,
    atoms_to_use: &[AtomId],
    masks: &FragmentSelectionMasks,
    source_rings: Option<&RingInfo>,
    existing_valence: Option<&'a ValenceAssignment>,
) -> Result<PreparedFragmentStereo<'a>, FragmentWriteInputError> {
    let record = record.into();
    // BEGIN RDKIT CPP FUNCTION MolFragmentToSmiles stereo preparation
    // RDKit✔️✔️:   ROMol tmol(mol, true);
    // RDKit✔️✔️:   // copy over the rings that only involve atoms/bonds in this fragment:
    // RDKit✔️✔️:   if (mol.getRingInfo()->isInitialized()) {
    // RDKit✔️✔️:     tmol.getRingInfo()->reset();
    // RDKit✔️✔️:     tmol.getRingInfo()->initialize();
    // RDKit✔️✔️:   if (tmol.needsUpdatePropertyCache()) {
    // RDKit✔️✔️:     for (auto atom : tmol.atoms()) {
    // RDKit✔️✔️:       atom->updatePropertyCache(false);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (params.doIsomericSmiles) {
    // RDKit✔️✔️:     if (!mol.hasProp(common_properties::_StereochemDone)) {
    // RDKit✔️✔️:       MolOps::assignStereochemistry(tmol, true);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       tmol.setProp(common_properties::_StereochemDone, 1);
    // RDKit✔️✔️:       for (auto aidx : atomsToUse) {
    // RDKit✔️✔️:         const Atom *oAt = mol.getAtomWithIdx(aidx);
    // RDKit✔️✔️:         std::string cipCode;
    // RDKit✔️✔️:         if (oAt->getPropIfPresent(common_properties::_CIPCode, cipCode)) {
    // RDKit✔️✔️:           tmol.getAtomWithIdx(aidx)->setProp(common_properties::_CIPCode,
    // RDKit✔️✔️:                                              cipCode);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION MolFragmentToSmiles stereo preparation
    // Behavior review: the optional source ring state remains absent or is
    // materialized from retained source rows before valence and stereo. The
    // existing core legacy owner performs its source fast-ring prepass when
    // OtherOrUnknown is supplied. A private empty unknown argument serves
    // only the absent-input dispatch and does not change retained state.
    // Complexity review: the quick-copy topology clone and retained-ring
    // construction each copy their source-shaped rows once. Valence borrows
    // valid supplied rows and otherwise computes all atoms; the legacy owner
    // performs at most one source-required fast-ring search.
    let mut topology = record.topology.clone();
    let retained_rings = source_rings
        .filter(|rings| rings.is_initialized())
        .map(|rings| {
            let rows = retain_fragment_rings(rings, masks);
            cosmolkit_core::ring_info_from_selected_rows(
                topology.atoms.len(),
                topology.bonds.len(),
                &rows.atom_rings,
                &rows.bond_rings,
            )
        })
        .transpose()?;
    let valence = prepare_fragment_valence(&topology, existing_valence)?;
    let mut stereochem_done = false;
    let mut ranking_rings = None;
    if params.isomeric_smiles {
        if record.properties.prop("_StereochemDone").is_none() {
            let absent_rings;
            let rings = if let Some(rings) = retained_rings.as_ref() {
                rings
            } else {
                absent_rings = RingInfo::new(
                    RingFindType::OtherOrUnknown,
                    topology.atoms.len(),
                    topology.bonds.len(),
                );
                &absent_rings
            };
            // RDKit✔️✔️:   if (!mol.getRingInfo()->isFindFastOrBetter()) {
            // RDKit✔️✔️:     MolOps::fastFindRings(mol);
            // The legacy owner and the later canonical owner both consume the
            // source Fast state. Materialize it once at this source boundary,
            // while retaining the copied OtherOrUnknown rows separately.
            // Cost: one O(V+E) search only when the legacy branch needs it;
            // both later owners borrow the result instead of repeating it.
            if !rings.is_find_fast_or_better() {
                ranking_rings = Some(cosmolkit_core::fast_find_rings(&topology)?);
            }
            let legacy_rings = ranking_rings.as_ref().unwrap_or(rings);
            topology = cosmolkit_core::assign_legacy_stereochemistry_with_flags(
                topology,
                &valence,
                legacy_rings,
                true,
                false,
            )?;
        } else {
            for &atom_id in atoms_to_use {
                if let Some(value) = record.topology.atoms[atom_id.index()].prop("_CIPCode") {
                    let code = cosmolkit_core::property_value_to_string(value)?;
                    topology.atoms[atom_id.index()].set_prop("_CIPCode", code)?;
                }
            }
        }
        stereochem_done = true;
    }
    Ok(PreparedFragmentStereo {
        topology,
        retained_rings,
        ranking_rings,
        valence,
        stereochem_done,
    })
}

pub(super) fn rank_prepared_fragment(
    prepared: &PreparedFragmentStereo<'_>,
    masks: &FragmentSelectionMasks,
    atom_symbols: Option<&[String]>,
    bond_symbols: Option<&[String]>,
    params: &SmilesWriteParams,
) -> Result<Vec<usize>, FragmentWriteInputError> {
    // BEGIN RDKIT CPP FUNCTION MolFragmentToSmiles ranking dispatch
    // RDKit✔️✔️:   if (params.canonical) {
    // RDKit✔️✔️:     bool breakTies = true;
    // RDKit✔️✔️:     bool includeChiralPresence = false;
    // RDKit✔️✔️:     bool includeRingStereo = true;
    // RDKit✔️✔️:     Canon::rankFragmentAtoms(
    // RDKit✔️✔️:         tmol, ranks, atomsInPlay, bondsInPlay, atomSymbols, bondSymbols,
    // RDKit✔️✔️:         breakTies, params.doIsomericSmiles, params.doIsomericSmiles,
    // RDKit✔️✔️:         !params.ignoreAtomMapNumbers, includeChiralPresence, includeRingStereo);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     for (unsigned int i = 0; i < tmol.getNumAtoms(); ++i) {
    // RDKit✔️✔️:       ranks[i] = i;
    // Only canonical dispatch checks exact symbol lengths in the core owner.
    // The index branch is O(V), while canonical dispatch borrows the S59
    // valence/Fast-ring preparation and performs the owner's rank refinement.
    if !params.canonical {
        return Ok((0..prepared.topology.atoms.len()).collect());
    }
    let mut rank_params = CanonicalRankParams::default();
    rank_params.break_ties = true;
    rank_params.include_chirality = params.isomeric_smiles;
    rank_params.include_isotopes = params.isomeric_smiles;
    rank_params.include_atom_maps = !params.ignore_atom_map_numbers;
    rank_params.include_chiral_presence = false;
    rank_params.include_ring_stereo = true;
    Ok(cosmolkit_core::rank_fragment_atoms_with_prepared_state(
        &prepared.topology,
        prepared.valence.as_ref(),
        prepared
            .ranking_rings
            .as_ref()
            .or(prepared.retained_rings.as_ref()),
        &masks.atoms_in_play,
        &masks.bonds_in_play,
        atom_symbols,
        bond_symbols,
        &rank_params,
    )?)
}

pub(super) fn prepare_fragment_valence<'a>(
    topology: &TopologyBlock,
    existing: Option<&'a ValenceAssignment>,
) -> Result<Cow<'a, ValenceAssignment>, FragmentWriteInputError> {
    // BEGIN RDKIT CPP FUNCTION ROMol::needsUpdatePropertyCache
    // RDKit✔️✔️: bool ROMol::needsUpdatePropertyCache() const {
    // RDKit✔️✔️:   for (const auto atom : atoms()) {
    // RDKit✔️✔️:     if (atom->needsUpdatePropertyCache()) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // there is no test for bonds yet since they do not obtain a valence property
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION ROMol::needsUpdatePropertyCache
    // BEGIN RDKIT CPP FUNCTION Atom::needsUpdatePropertyCache
    // RDKit✔️✔️: bool Atom::needsUpdatePropertyCache() const {
    // RDKit✔️✔️:   return !(this->d_explicitValence >= 0 &&
    // RDKit✔️✔️:            (this->df_noImplicit || this->d_implicitValence >= 0));
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Atom::needsUpdatePropertyCache
    // BEGIN RDKIT CPP FUNCTION MolFragmentToSmiles cache preparation
    // RDKit✔️✔️:   if (tmol.needsUpdatePropertyCache()) {
    // RDKit✔️✔️:     for (auto atom : tmol.atoms()) {
    // RDKit✔️✔️:       atom->updatePropertyCache(false);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION MolFragmentToSmiles cache preparation
    // A supplied assignment is the only detached cache evidence. The O(A)
    // validity scan borrows all-valid rows without allocation; if any row
    // needs an update, the existing core owner calculates all rows with
    // strict=false, matching source all-atom scope. None has no supplied
    // cache and takes that same complete calculation path.
    if let Some(existing) = existing {
        if existing.explicit_valence.len() != topology.atoms.len()
            || existing.implicit_hydrogens.len() != topology.atoms.len()
        {
            return Err(FragmentWriteInputError::ValenceAssignmentLengthMismatch {
                atom_count: topology.atoms.len(),
                explicit_len: existing.explicit_valence.len(),
                implicit_len: existing.implicit_hydrogens.len(),
            });
        }
        if topology.atoms.iter().enumerate().all(|(index, atom)| {
            existing.explicit_valence[index] >= 0
                && (atom.no_implicit() || existing.implicit_hydrogens[index] >= 0)
        }) {
            return Ok(Cow::Borrowed(existing));
        }
    }
    Ok(Cow::Owned(
        cosmolkit_core::assign_valence_with_options_for_topology(
            topology,
            ValenceModel::RdkitLike,
            false,
        )?,
    ))
}

pub(super) fn validate_fragment_write_inputs(
    topology: &TopologyBlock,
    params: &SmilesWriteParams,
    atoms_to_use: &[AtomId],
    atom_symbols: Option<&[String]>,
    bond_symbols: Option<&[String]>,
) -> Result<(), FragmentWriteInputError> {
    // BEGIN RDKIT CPP FUNCTION MolFragmentToSmiles preconditions
    // RDKit✔️✔️:   PRECONDITION(atomsToUse.size(), "no atoms provided");
    // RDKit✔️✔️:   PRECONDITION(
    // RDKit✔️✔️:       params.rootedAtAtom < 0 ||
    // RDKit✔️✔️:           static_cast<unsigned int>(params.rootedAtAtom) < mol.getNumAtoms(),
    // RDKit✔️✔️:       "rootedAtomAtom must be less than the number of atoms");
    // RDKit✔️✔️:   PRECONDITION(params.rootedAtAtom < 0 ||
    // RDKit✔️✔️:                    std::find(atomsToUse.begin(), atomsToUse.end(),
    // RDKit✔️✔️:                              params.rootedAtAtom) != atomsToUse.end(),
    // RDKit✔️✔️:                "rootedAtAtom not found in atomsToUse");
    // RDKit✔️✔️:   PRECONDITION(!atomSymbols || atomSymbols->size() >= mol.getNumAtoms(),
    // RDKit✔️✔️:                "bad atomSymbols vector");
    // RDKit✔️✔️:   PRECONDITION(!bondSymbols || bondSymbols->size() >= mol.getNumBonds(),
    // RDKit✔️✔️:                "bad bondSymbols vector");
    // END RDKIT CPP FUNCTION MolFragmentToSmiles preconditions
    // The checks preserve the source order and cost shape: constant-time
    // size/range checks, one linear root-membership scan, and no allocation.
    if atoms_to_use.is_empty() {
        return Err(FragmentWriteInputError::NoAtomsProvided);
    }
    if let Some(root) = params.rooted_at_atom {
        if root.index() >= topology.atoms.len() {
            return Err(FragmentWriteInputError::RootAtomOutOfRange {
                atom_index: root.index(),
                atom_count: topology.atoms.len(),
            });
        }
        if !atoms_to_use.contains(&root) {
            return Err(FragmentWriteInputError::RootAtomNotSelected {
                atom_index: root.index(),
            });
        }
    }
    if let Some(symbols) = atom_symbols
        && symbols.len() < topology.atoms.len()
    {
        return Err(FragmentWriteInputError::AtomSymbolsTooShort {
            len: symbols.len(),
            expected: topology.atoms.len(),
        });
    }
    if let Some(symbols) = bond_symbols
        && symbols.len() < topology.bonds.len()
    {
        return Err(FragmentWriteInputError::BondSymbolsTooShort {
            len: symbols.len(),
            expected: topology.bonds.len(),
        });
    }
    Ok(())
}

pub(super) fn build_fragment_selection_masks(
    topology: &TopologyBlock,
    params: &SmilesWriteParams,
    atoms_to_use: &[AtomId],
    bonds_to_use: Option<&[BondId]>,
) -> Result<FragmentSelectionMasks, FragmentWriteInputError> {
    // BEGIN RDKIT CPP FUNCTION MolFragmentToSmiles selection masks
    // RDKit✔️✔️:   boost::dynamic_bitset<> atomsInPlay(mol.getNumAtoms(), 0);
    // RDKit✔️✔️:   for (auto aidx : atomsToUse) {
    // RDKit✔️✔️:     atomsInPlay.set(aidx);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   boost::dynamic_bitset<> bondsInPlay(mol.getNumBonds(), 0);
    // RDKit✔️✔️:   if (bondsToUse) {
    // RDKit✔️✔️:     for (auto bidx : *bondsToUse) {
    // RDKit✔️✔️:       bondsInPlay.set(bidx);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     PRECONDITION(
    // RDKit✔️✔️:         params.rootedAtAtom < 0 || MolOps::getMolFrags(mol).size() == 1,
    // RDKit✔️✔️:         "rootedAtAtom can only be used with molecules that have a single fragment");
    // RDKit✔️✔️:     for (auto aidx : atomsToUse) {
    // RDKit✔️✔️:       for (const auto &bndi : boost::make_iterator_range(
    // RDKit✔️✔️:                mol.getAtomBonds(mol.getAtomWithIdx(aidx)))) {
    // RDKit✔️✔️:         const Bond *bond = mol[bndi];
    // RDKit✔️✔️:         if (atomsInPlay[bond->getOtherAtomIdx(aidx)]) {
    // RDKit✔️✔️:           bondsInPlay.set(bond->getIdx());
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION MolFragmentToSmiles selection masks
    // Typed range errors replace the source bitset assertions. For valid
    // inputs the explicit path is O(A+B), while the induced path reuses the
    // O(V+E) component owner and scans only selected-atom adjacency rows.
    let mut atoms_in_play = vec![false; topology.atoms.len()];
    for atom in atoms_to_use {
        let Some(selected) = atoms_in_play.get_mut(atom.index()) else {
            return Err(FragmentWriteInputError::AtomOutOfRange {
                atom_index: atom.index(),
                atom_count: topology.atoms.len(),
            });
        };
        *selected = true;
    }

    let mut bonds_in_play = vec![false; topology.bonds.len()];
    if let Some(bonds_to_use) = bonds_to_use {
        for bond in bonds_to_use {
            let Some(selected) = bonds_in_play.get_mut(bond.index()) else {
                return Err(FragmentWriteInputError::BondOutOfRange {
                    bond_index: bond.index(),
                    bond_count: topology.bonds.len(),
                });
            };
            *selected = true;
        }
    } else {
        if let Some(root) = params.rooted_at_atom {
            let components = cosmolkit_core::connected_components(topology)
                .map_err(|error| FragmentWriteInputError::InvalidTopology(error.to_string()))?;
            if components.components.len() != 1 {
                return Err(FragmentWriteInputError::RootAtomRequiresSingleFragment {
                    atom_index: root.index(),
                    fragment_count: components.components.len(),
                });
            }
        }
        for atom in atoms_to_use {
            for neighbor in topology.adjacency.neighbors_of(atom.index()) {
                if atoms_in_play[neighbor.atom_index] {
                    bonds_in_play[neighbor.bond.index()] = true;
                }
            }
        }
    }

    Ok(FragmentSelectionMasks {
        atoms_in_play,
        bonds_in_play,
    })
}

pub(super) fn retain_fragment_rings(
    rings: &RingInfo,
    masks: &FragmentSelectionMasks,
) -> FragmentRingRows {
    // BEGIN RDKIT CPP FUNCTION MolFragmentToSmiles ring transport
    // RDKit✔️✔️:   // copy over the rings that only involve atoms/bonds in this fragment:
    // RDKit✔️✔️:   if (mol.getRingInfo()->isInitialized()) {
    // RDKit✔️✔️:     tmol.getRingInfo()->reset();
    // RDKit✔️✔️:     tmol.getRingInfo()->initialize();
    // RDKit✔️✔️:     for (unsigned int ridx = 0; ridx < mol.getRingInfo()->numRings(); ++ridx) {
    // RDKit✔️✔️:       const INT_VECT &aring = mol.getRingInfo()->atomRings()[ridx];
    // RDKit✔️✔️:       bool keepIt = true;
    // RDKit✔️✔️:       for (auto aidx : aring) {
    // RDKit✔️✔️:         if (!atomsInPlay[aidx]) {
    // RDKit✔️✔️:           keepIt = false;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (keepIt) {
    // RDKit✔️✔️:         const INT_VECT &bring = mol.getRingInfo()->bondRings()[ridx];
    // RDKit✔️✔️:         for (auto bidx : bring) {
    // RDKit✔️✔️:           if (!bondsInPlay[bidx]) {
    // RDKit✔️✔️:             keepIt = false;
    // RDKit✔️✔️:             break;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         if (keepIt) {
    // RDKit✔️✔️:           tmol.getRingInfo()->addRing(aring, bring);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION MolFragmentToSmiles ring transport
    // RingInfo rows and fragment masks share the validated source topology's
    // ID domain. Each retained source row is cloned once, as addRing does;
    // iteration remains linear in the total number of inspected ring members.
    if !rings.is_initialized() {
        return FragmentRingRows {
            initialized: false,
            atom_rings: Vec::new(),
            bond_rings: Vec::new(),
        };
    }

    let mut atom_rings = Vec::new();
    let mut bond_rings = Vec::new();
    for (atom_ring, bond_ring) in rings.atom_rings().iter().zip(rings.bond_rings()) {
        if atom_ring
            .iter()
            .all(|atom| masks.atoms_in_play[atom.index()])
            && bond_ring
                .iter()
                .all(|bond| masks.bonds_in_play[bond.index()])
        {
            atom_rings.push(atom_ring.clone());
            bond_rings.push(bond_ring.clone());
        }
    }

    FragmentRingRows {
        initialized: true,
        atom_rings,
        bond_rings,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse_smiles;
    use cosmolkit_types::BondStereo;

    fn three_atom_topology() -> TopologyBlock {
        parse_smiles("CCC", &Default::default())
            .expect("fixed fragment-preflight molecule parses")
            .topology
    }

    fn topology(input: &str) -> TopologyBlock {
        parse_smiles(input, &Default::default())
            .expect("fixed fragment-selection molecule parses")
            .topology
    }

    fn s59_record(input: &str) -> SmilesRecord {
        let mut record =
            parse_smiles(input, &Default::default()).expect("fixed S59 molecule parses");
        record.properties.clear_prop("_StereochemDone");
        record
    }

    #[test]
    fn fragment_ranking_custom_symbols_use_original_masks_and_prepared_state() {
        let record = s59_record("CCO");
        let before = record.clone();
        let params = SmilesWriteParams::default();
        let selected = [AtomId::new(0), AtomId::new(1), AtomId::new(2)];
        let masks =
            build_fragment_selection_masks(&record.topology, &params, &selected, None).unwrap();
        let prepared =
            prepare_fragment_stereo(&record, &params, &selected, &masks, None, None).unwrap();
        assert!(
            prepared
                .ranking_rings
                .as_ref()
                .unwrap()
                .is_find_fast_or_better()
        );
        let forward = ["N".to_owned(), "O".to_owned(), "S".to_owned()];
        let reverse = ["S".to_owned(), "O".to_owned(), "N".to_owned()];
        assert_eq!(
            rank_prepared_fragment(&prepared, &masks, Some(&forward), None, &params).unwrap(),
            vec![0, 2, 1]
        );
        assert_eq!(
            rank_prepared_fragment(&prepared, &masks, Some(&reverse), None, &params).unwrap(),
            vec![1, 2, 0]
        );
        let selected_masks = FragmentSelectionMasks {
            atoms_in_play: vec![true, true, false],
            bonds_in_play: vec![true, false],
        };
        assert_eq!(
            rank_prepared_fragment(&prepared, &selected_masks, None, None, &params).unwrap(),
            cosmolkit_core::rank_fragment_atoms_with_prepared_state(
                &prepared.topology,
                prepared.valence.as_ref(),
                prepared.ranking_rings.as_ref(),
                &selected_masks.atoms_in_play,
                &selected_masks.bonds_in_play,
                None,
                None,
                &CanonicalRankParams::default(),
            )
            .unwrap()
        );
        assert_eq!(record, before);
    }

    #[test]
    fn fragment_ranking_writer_flags_control_isotope_and_map_policy() {
        let mut record = s59_record("[13CH3]CC");
        record.topology.atoms[2].set_atom_map(Some(9));
        let selected = [AtomId::new(0), AtomId::new(1), AtomId::new(2)];
        let plain_params = SmilesWriteParams {
            isomeric_smiles: false,
            ..Default::default()
        };
        let masks =
            build_fragment_selection_masks(&record.topology, &plain_params, &selected, None)
                .unwrap();
        let plain =
            prepare_fragment_stereo(&record, &plain_params, &selected, &masks, None, None).unwrap();
        assert_eq!(
            rank_prepared_fragment(&plain, &masks, None, None, &plain_params).unwrap(),
            vec![0, 1, 2]
        );
        let isomeric_params = SmilesWriteParams::default();
        let isomeric =
            prepare_fragment_stereo(&record, &isomeric_params, &selected, &masks, None, None)
                .unwrap();
        assert_eq!(
            rank_prepared_fragment(&isomeric, &masks, None, None, &isomeric_params).unwrap(),
            vec![0, 1, 2]
        );
        let ignore_maps = SmilesWriteParams {
            ignore_atom_map_numbers: true,
            ..Default::default()
        };
        let mut expected_params = CanonicalRankParams::default();
        expected_params.include_atom_maps = false;
        let expected_without_maps = cosmolkit_core::rank_fragment_atoms_with_prepared_state(
            &isomeric.topology,
            isomeric.valence.as_ref(),
            isomeric.ranking_rings.as_ref(),
            &masks.atoms_in_play,
            &masks.bonds_in_play,
            None,
            None,
            &expected_params,
        )
        .unwrap();
        assert_eq!(
            rank_prepared_fragment(&isomeric, &masks, None, None, &ignore_maps).unwrap(),
            expected_without_maps
        );
        assert_eq!(expected_without_maps, vec![1, 2, 0]);
    }

    #[test]
    fn fragment_ranking_chirality_policy_matches_pinned_opposite_centers() {
        let mut record =
            parse_smiles("[C@H](F)(Cl)Br.[C@@H](F)(Cl)Br", &Default::default()).unwrap();
        record.properties.set_prop("_StereochemDone", "1").unwrap();
        let selected = (0..record.topology.atoms.len())
            .map(AtomId::new)
            .collect::<Vec<_>>();
        let plain_params = SmilesWriteParams {
            isomeric_smiles: false,
            ..Default::default()
        };
        let masks =
            build_fragment_selection_masks(&record.topology, &plain_params, &selected, None)
                .unwrap();
        let plain =
            prepare_fragment_stereo(&record, &plain_params, &selected, &masks, None, None).unwrap();
        assert_eq!(
            rank_prepared_fragment(&plain, &masks, None, None, &plain_params).unwrap(),
            vec![6, 0, 2, 4, 7, 1, 3, 5]
        );
        let isomeric_params = SmilesWriteParams::default();
        let isomeric =
            prepare_fragment_stereo(&record, &isomeric_params, &selected, &masks, None, None)
                .unwrap();
        assert_eq!(
            rank_prepared_fragment(&isomeric, &masks, None, None, &isomeric_params).unwrap(),
            vec![7, 1, 3, 5, 6, 0, 2, 4]
        );
    }

    #[test]
    fn fragment_ranking_noncanonical_skips_canonical_symbol_checks() {
        let record = s59_record("CCO");
        let selected = [AtomId::new(0), AtomId::new(1), AtomId::new(2)];
        let params = SmilesWriteParams {
            canonical: false,
            ..Default::default()
        };
        let masks =
            build_fragment_selection_masks(&record.topology, &params, &selected, None).unwrap();
        let prepared =
            prepare_fragment_stereo(&record, &params, &selected, &masks, None, None).unwrap();
        let oversized = vec!["C".to_owned(); 4];
        assert_eq!(
            validate_fragment_write_inputs(
                &record.topology,
                &params,
                &selected,
                Some(&oversized),
                None,
            ),
            Ok(())
        );
        assert_eq!(
            rank_prepared_fragment(&prepared, &masks, Some(&oversized), None, &params).unwrap(),
            vec![0, 1, 2]
        );
        let canonical = SmilesWriteParams {
            canonical: true,
            ..params
        };
        assert!(matches!(
            rank_prepared_fragment(&prepared, &masks, Some(&oversized), None, &canonical),
            Err(FragmentWriteInputError::Rank(
                CanonicalRankError::AtomSymbolLength {
                    expected: 3,
                    actual: 4,
                }
            ))
        ));
    }

    #[test]
    fn fragment_stereo_preparation_retains_rows_before_legacy_and_recomputes_at_dispatch() {
        let mut record = s59_record("C1=CCCCC1.C1CCCCC1");
        record.topology.bonds[0]
            .set_stereo(BondStereo::Any)
            .unwrap();
        let before = record.clone();
        let source_rings = cosmolkit_core::fast_find_rings(&record.topology).unwrap();
        assert_eq!(source_rings.num_rings(), 2);
        let source_rings_before = source_rings.clone();
        let selected = (6..12).map(AtomId::new).collect::<Vec<_>>();
        let masks = build_fragment_selection_masks(
            &record.topology,
            &SmilesWriteParams::default(),
            &selected,
            None,
        )
        .unwrap();
        let prepared = prepare_fragment_stereo(
            &record,
            &SmilesWriteParams::default(),
            &selected,
            &masks,
            Some(&source_rings),
            None,
        )
        .unwrap();
        let retained = prepared.retained_rings.as_ref().unwrap();
        assert!(retained.is_initialized());
        assert_eq!(retained.find_type(), RingFindType::OtherOrUnknown);
        assert_eq!(retained.num_rings(), 1);
        assert_eq!(retained.num_bond_rings(BondId::new(0)), 0);
        assert_eq!(
            prepared.topology.bonds[0].stereo(),
            BondStereo::None,
            "the pinned legacy consumer searches the full topology after retention"
        );
        assert!(prepared.stereochem_done);
        assert!(matches!(prepared.valence, Cow::Owned(_)));
        assert_eq!(record, before);
        assert_eq!(source_rings, source_rings_before);
    }

    #[test]
    fn fragment_stereo_preparation_keeps_initialized_empty_distinct_from_absent() {
        let mut record = s59_record("C1=CCCCC1");
        record.topology.bonds[0]
            .set_stereo(BondStereo::Any)
            .unwrap();
        let source_rings = cosmolkit_core::fast_find_rings(&record.topology).unwrap();
        let masks = FragmentSelectionMasks {
            atoms_in_play: vec![false; record.topology.atoms.len()],
            bonds_in_play: vec![false; record.topology.bonds.len()],
        };
        let selected = [AtomId::new(0)];
        let initialized = prepare_fragment_stereo(
            &record,
            &SmilesWriteParams::default(),
            &selected,
            &masks,
            Some(&source_rings),
            None,
        )
        .unwrap();
        assert_eq!(initialized.retained_rings.as_ref().unwrap().num_rings(), 0);
        assert!(
            initialized
                .retained_rings
                .as_ref()
                .unwrap()
                .is_initialized()
        );
        assert_eq!(initialized.topology.bonds[0].stereo(), BondStereo::None);

        let absent = prepare_fragment_stereo(
            &record,
            &SmilesWriteParams::default(),
            &selected,
            &masks,
            None,
            None,
        )
        .unwrap();
        assert!(absent.retained_rings.is_none());
        assert_eq!(absent.topology.bonds[0].stereo(), BondStereo::None);
    }

    #[test]
    fn fragment_stereo_preparation_nonisomeric_borrows_or_recalculates_valence() {
        let record = s59_record("CCC");
        let before = record.clone();
        let masks = FragmentSelectionMasks {
            atoms_in_play: vec![true; 3],
            bonds_in_play: vec![true; 2],
        };
        let params = SmilesWriteParams {
            isomeric_smiles: false,
            ..Default::default()
        };
        let existing = cosmolkit_core::assign_valence_with_options_for_topology(
            &record.topology,
            ValenceModel::RdkitLike,
            false,
        )
        .unwrap();
        let selected = [AtomId::new(0), AtomId::new(1), AtomId::new(2)];
        let borrowed =
            prepare_fragment_stereo(&record, &params, &selected, &masks, None, Some(&existing))
                .unwrap();
        assert!(matches!(borrowed.valence, Cow::Borrowed(_)));
        assert!(std::ptr::eq(borrowed.valence.as_ref(), &existing));
        assert!(!borrowed.stereochem_done);
        assert!(borrowed.retained_rings.is_none());
        assert_eq!(borrowed.topology, record.topology);

        let mut invalid = existing.clone();
        invalid.explicit_valence[0] = -1;
        invalid.explicit_valence[2] = 77;
        let recalculated =
            prepare_fragment_stereo(&record, &params, &selected, &masks, None, Some(&invalid))
                .unwrap();
        assert!(matches!(recalculated.valence, Cow::Owned(_)));
        assert_eq!(*recalculated.valence, existing);
        assert_eq!(record, before);
    }

    #[test]
    fn fragment_stereo_preparation_marker_presence_copies_selected_cip_without_assignment() {
        let mut record = s59_record("C=CC");
        record.topology.bonds[0]
            .set_stereo(BondStereo::Any)
            .unwrap();
        record.topology.atoms[0]
            .set_prop("_CIPCode", 7_i32)
            .unwrap();
        record.topology.atoms[1].set_prop("_CIPCode", "S").unwrap();
        record.properties.set_prop("_StereochemDone", "0").unwrap();
        let before = record.clone();
        let masks = FragmentSelectionMasks {
            atoms_in_play: vec![true, false, false],
            bonds_in_play: vec![false; 2],
        };
        let expected_valence = cosmolkit_core::assign_valence_with_options_for_topology(
            &record.topology,
            ValenceModel::RdkitLike,
            false,
        )
        .unwrap();
        let mut stale_valence = expected_valence.clone();
        stale_valence.explicit_valence[0] = -1;
        stale_valence.explicit_valence[2] = 77;
        let stale_before = stale_valence.clone();
        let prepared = prepare_fragment_stereo(
            &record,
            &SmilesWriteParams::default(),
            &[AtomId::new(0)],
            &masks,
            None,
            Some(&stale_valence),
        )
        .unwrap();
        assert!(prepared.stereochem_done);
        assert!(matches!(prepared.valence, Cow::Owned(_)));
        assert_eq!(*prepared.valence, expected_valence);
        assert_eq!(prepared.topology.bonds[0].stereo(), BondStereo::Any);
        assert_eq!(
            prepared.topology.atoms[0]
                .prop("_CIPCode")
                .unwrap()
                .as_string()
                .unwrap()
                .as_bytes(),
            b"7".as_slice()
        );
        assert_eq!(
            prepared.topology.atoms[1]
                .prop("_CIPCode")
                .unwrap()
                .as_string()
                .unwrap()
                .as_bytes(),
            b"S".as_slice()
        );
        assert_eq!(record, before);
        assert_eq!(stale_valence, stale_before);
    }

    #[test]
    fn fragment_valence_preparation_runs_for_nonisomeric_absent_cache() {
        let topology = topology("CCC");
        let before = topology.clone();
        let params = SmilesWriteParams {
            isomeric_smiles: false,
            ..Default::default()
        };
        assert!(!params.isomeric_smiles);
        let expected = cosmolkit_core::assign_valence_with_options_for_topology(
            &topology,
            ValenceModel::RdkitLike,
            false,
        )
        .unwrap();

        let prepared = prepare_fragment_valence(&topology, None).unwrap();
        assert!(matches!(prepared, Cow::Owned(_)));
        assert_eq!(*prepared, expected);
        assert_eq!(topology, before);
    }

    #[test]
    fn fragment_valence_preparation_borrows_valid_rows_unchanged() {
        let topology = topology("CCC");
        let existing = cosmolkit_core::assign_valence_with_options_for_topology(
            &topology,
            ValenceModel::RdkitLike,
            false,
        )
        .unwrap();
        let before = existing.clone();
        let prepared = prepare_fragment_valence(&topology, Some(&existing)).unwrap();
        assert!(matches!(prepared, Cow::Borrowed(_)));
        assert!(std::ptr::eq(prepared.as_ref(), &existing));
        assert_eq!(existing, before);
    }

    #[test]
    fn fragment_valence_preparation_recalculates_every_row_for_invalid_cache() {
        let topology = topology("CCC");
        let expected = cosmolkit_core::assign_valence_with_options_for_topology(
            &topology,
            ValenceModel::RdkitLike,
            false,
        )
        .unwrap();
        let mut bad_explicit = expected.clone();
        bad_explicit.explicit_valence[0] = -1;
        bad_explicit.explicit_valence[2] = 77;
        let before_explicit = bad_explicit.clone();
        let recalculated = prepare_fragment_valence(&topology, Some(&bad_explicit)).unwrap();
        assert!(matches!(recalculated, Cow::Owned(_)));
        assert_eq!(*recalculated, expected);
        assert_eq!(bad_explicit, before_explicit);

        let mut bad_implicit = expected.clone();
        bad_implicit.implicit_hydrogens[1] = -1;
        bad_implicit.explicit_valence[2] = 77;
        let before_implicit = bad_implicit.clone();
        let recalculated = prepare_fragment_valence(&topology, Some(&bad_implicit)).unwrap();
        assert!(matches!(recalculated, Cow::Owned(_)));
        assert_eq!(*recalculated, expected);
        assert_eq!(bad_implicit, before_implicit);
    }

    #[test]
    fn fragment_valence_preparation_honors_no_implicit_and_rejects_bad_lengths() {
        let mut topology = topology("CCC");
        topology.atoms[0].set_no_implicit(true);
        let mut existing = cosmolkit_core::assign_valence_with_options_for_topology(
            &topology,
            ValenceModel::RdkitLike,
            false,
        )
        .unwrap();
        existing.implicit_hydrogens[0] = -1;
        assert!(matches!(
            prepare_fragment_valence(&topology, Some(&existing)).unwrap(),
            Cow::Borrowed(_)
        ));

        let mut wrong_explicit = existing.clone();
        wrong_explicit.explicit_valence.pop();
        assert_eq!(
            prepare_fragment_valence(&topology, Some(&wrong_explicit)),
            Err(FragmentWriteInputError::ValenceAssignmentLengthMismatch {
                atom_count: 3,
                explicit_len: 2,
                implicit_len: 3,
            })
        );
        let mut wrong_implicit = existing;
        wrong_implicit.implicit_hydrogens.pop();
        assert_eq!(
            prepare_fragment_valence(&topology, Some(&wrong_implicit)),
            Err(FragmentWriteInputError::ValenceAssignmentLengthMismatch {
                atom_count: 3,
                explicit_len: 3,
                implicit_len: 2,
            })
        );
    }

    #[test]
    fn fragment_preflight_rejects_empty_selection_before_other_errors() {
        let topology = three_atom_topology();
        let params = SmilesWriteParams {
            rooted_at_atom: Some(AtomId::new(9)),
            ..Default::default()
        };
        let atom_symbols = vec!["C".to_owned()];
        let bond_symbols = vec!["-".to_owned()];

        assert_eq!(
            validate_fragment_write_inputs(
                &topology,
                &params,
                &[],
                Some(&atom_symbols),
                Some(&bond_symbols),
            ),
            Err(FragmentWriteInputError::NoAtomsProvided)
        );
    }

    #[test]
    fn fragment_preflight_checks_root_range_before_membership_and_symbols() {
        let topology = three_atom_topology();
        let params = SmilesWriteParams {
            rooted_at_atom: Some(AtomId::new(3)),
            ..Default::default()
        };
        let atom_symbols = vec!["C".to_owned()];
        let bond_symbols = vec!["-".to_owned()];

        assert_eq!(
            validate_fragment_write_inputs(
                &topology,
                &params,
                &[AtomId::new(0)],
                Some(&atom_symbols),
                Some(&bond_symbols),
            ),
            Err(FragmentWriteInputError::RootAtomOutOfRange {
                atom_index: 3,
                atom_count: 3,
            })
        );
    }

    #[test]
    fn fragment_preflight_checks_root_membership_before_symbols() {
        let topology = three_atom_topology();
        let params = SmilesWriteParams {
            rooted_at_atom: Some(AtomId::new(2)),
            ..Default::default()
        };
        let atom_symbols = vec!["C".to_owned()];
        let bond_symbols = vec!["-".to_owned()];

        assert_eq!(
            validate_fragment_write_inputs(
                &topology,
                &params,
                &[AtomId::new(0), AtomId::new(1)],
                Some(&atom_symbols),
                Some(&bond_symbols),
            ),
            Err(FragmentWriteInputError::RootAtomNotSelected { atom_index: 2 })
        );
    }

    #[test]
    fn fragment_preflight_checks_atom_symbols_before_bond_symbols() {
        let topology = three_atom_topology();
        let params = SmilesWriteParams {
            rooted_at_atom: Some(AtomId::new(0)),
            ..Default::default()
        };
        let short_atom_symbols = vec!["C".to_owned(), "C".to_owned()];
        let short_bond_symbols = vec!["-".to_owned()];

        assert_eq!(
            validate_fragment_write_inputs(
                &topology,
                &params,
                &[AtomId::new(0)],
                Some(&short_atom_symbols),
                Some(&short_bond_symbols),
            ),
            Err(FragmentWriteInputError::AtomSymbolsTooShort {
                len: 2,
                expected: 3,
            })
        );

        let exact_atom_symbols = vec!["C".to_owned(); 3];
        assert_eq!(
            validate_fragment_write_inputs(
                &topology,
                &params,
                &[AtomId::new(0)],
                Some(&exact_atom_symbols),
                Some(&short_bond_symbols),
            ),
            Err(FragmentWriteInputError::BondSymbolsTooShort {
                len: 1,
                expected: 2,
            })
        );
    }

    #[test]
    fn fragment_preflight_accepts_exact_and_oversized_symbol_tables() {
        let topology = three_atom_topology();
        let params = SmilesWriteParams {
            rooted_at_atom: Some(AtomId::new(1)),
            ..Default::default()
        };
        let before = topology.clone();
        let exact_atom_symbols = vec!["C".to_owned(); 3];
        let exact_bond_symbols = vec!["-".to_owned(); 2];
        let oversized_atom_symbols = vec!["C".to_owned(); 4];
        let oversized_bond_symbols = vec!["-".to_owned(); 3];

        for (atom_symbols, bond_symbols) in [
            (&exact_atom_symbols[..], &exact_bond_symbols[..]),
            (&oversized_atom_symbols[..], &oversized_bond_symbols[..]),
        ] {
            assert_eq!(
                validate_fragment_write_inputs(
                    &topology,
                    &params,
                    &[AtomId::new(0), AtomId::new(1)],
                    Some(atom_symbols),
                    Some(bond_symbols),
                ),
                Ok(())
            );
        }
        assert_eq!(topology, before, "preflight leaves the input unchanged");
    }

    #[test]
    fn fragment_selection_masks_distinguish_explicit_and_induced_ring_bonds() {
        let topology = topology("C1CCC1");
        let before = topology.clone();
        let atoms = [
            AtomId::new(0),
            AtomId::new(1),
            AtomId::new(2),
            AtomId::new(3),
        ];
        let explicit = [BondId::new(0), BondId::new(2)];

        let explicit_masks = build_fragment_selection_masks(
            &topology,
            &SmilesWriteParams::default(),
            &atoms,
            Some(&explicit),
        )
        .expect("explicit typed bond rows are valid");
        assert_eq!(explicit_masks.atoms_in_play, vec![true; 4]);
        assert_eq!(explicit_masks.bonds_in_play, vec![true, false, true, false]);

        let induced_masks =
            build_fragment_selection_masks(&topology, &SmilesWriteParams::default(), &atoms, None)
                .expect("all four ring bonds are induced");
        assert_eq!(induced_masks.atoms_in_play, vec![true; 4]);
        assert_eq!(induced_masks.bonds_in_play, vec![true; 4]);
        assert_eq!(topology, before, "mask construction leaves input unchanged");
    }

    #[test]
    fn fragment_selection_masks_apply_root_guard_only_without_explicit_bonds() {
        let topology = topology("CC.CC");
        let params = SmilesWriteParams {
            rooted_at_atom: Some(AtomId::new(0)),
            ..Default::default()
        };
        let atoms = [AtomId::new(0), AtomId::new(1)];

        assert_eq!(
            build_fragment_selection_masks(&topology, &params, &atoms, None),
            Err(FragmentWriteInputError::RootAtomRequiresSingleFragment {
                atom_index: 0,
                fragment_count: 2,
            })
        );

        let explicit_empty = build_fragment_selection_masks(&topology, &params, &atoms, Some(&[]))
            .expect("a present empty bond list bypasses the source guard");
        assert_eq!(explicit_empty.atoms_in_play, vec![true, true, false, false]);
        assert_eq!(explicit_empty.bonds_in_play, vec![false, false]);
    }

    #[test]
    fn fragment_selection_masks_report_typed_index_errors_in_source_order() {
        let topology = three_atom_topology();
        let invalid_bond = [BondId::new(2)];

        assert_eq!(
            build_fragment_selection_masks(
                &topology,
                &SmilesWriteParams::default(),
                &[AtomId::new(3)],
                Some(&invalid_bond),
            ),
            Err(FragmentWriteInputError::AtomOutOfRange {
                atom_index: 3,
                atom_count: 3,
            })
        );
        assert_eq!(
            build_fragment_selection_masks(
                &topology,
                &SmilesWriteParams::default(),
                &[AtomId::new(0)],
                Some(&invalid_bond),
            ),
            Err(FragmentWriteInputError::BondOutOfRange {
                bond_index: 2,
                bond_count: 2,
            })
        );
    }

    #[test]
    fn fragment_ring_transport_requires_every_ring_atom_and_bond() {
        let topology = topology("C1CCC1");
        let rings = cosmolkit_core::fast_find_rings(&topology)
            .expect("fixed four-membered ring has valid ring information");
        let complete = FragmentSelectionMasks {
            atoms_in_play: vec![true; topology.atoms.len()],
            bonds_in_play: vec![true; topology.bonds.len()],
        };
        let retained = retain_fragment_rings(&rings, &complete);
        assert!(retained.initialized);
        assert_eq!(retained.atom_rings, rings.atom_rings());
        assert_eq!(retained.bond_rings, rings.bond_rings());

        let mut missing_atom = complete.clone();
        missing_atom.atoms_in_play[rings.atom_rings()[0][0].index()] = false;
        assert_eq!(
            retain_fragment_rings(&rings, &missing_atom),
            FragmentRingRows {
                initialized: true,
                atom_rings: Vec::new(),
                bond_rings: Vec::new(),
            }
        );

        let mut missing_bond = complete;
        missing_bond.bonds_in_play[rings.bond_rings()[0][0].index()] = false;
        assert_eq!(
            retain_fragment_rings(&rings, &missing_bond),
            FragmentRingRows {
                initialized: true,
                atom_rings: Vec::new(),
                bond_rings: Vec::new(),
            }
        );
    }

    #[test]
    fn fragment_ring_transport_preserves_original_ids_and_source_order() {
        let topology = topology("CC.C1CC1");
        let before = topology.clone();
        let rings = cosmolkit_core::fast_find_rings(&topology)
            .expect("fixed offset triangle has valid ring information");
        assert_eq!(
            rings.atom_rings(),
            &[vec![AtomId::new(4), AtomId::new(3), AtomId::new(2)]]
        );
        assert_eq!(
            rings.bond_rings(),
            &[vec![BondId::new(2), BondId::new(1), BondId::new(3)]]
        );

        let masks = build_fragment_selection_masks(
            &topology,
            &SmilesWriteParams::default(),
            &[AtomId::new(2), AtomId::new(3), AtomId::new(4)],
            Some(&[BondId::new(1), BondId::new(2), BondId::new(3)]),
        )
        .expect("offset ring selection uses valid original IDs");
        assert_eq!(
            retain_fragment_rings(&rings, &masks),
            FragmentRingRows {
                initialized: true,
                atom_rings: vec![vec![AtomId::new(4), AtomId::new(3), AtomId::new(2)]],
                bond_rings: vec![vec![BondId::new(2), BondId::new(1), BondId::new(3)]],
            }
        );
        assert_eq!(topology, before, "ring transport leaves input unchanged");
    }
}
