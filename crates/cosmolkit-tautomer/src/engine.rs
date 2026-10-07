//! Detached TAU chemistry: existing source-backed transforms and stereo.
use crate::ordered::{TautomerExpandedProduct, TautomerExpansionAttempt};
use crate::{TautomerParams, TautomerTransform};
use cosmolkit_core::*;
use cosmolkit_model::{AtomId, BondId, CoordinateBlock, MoleculeProperties, TopologyBlock};
use cosmolkit_search::{
    CompiledQuery, SearchTarget, SubstructMatchParams, SubstructMatchResult,
    build_prepared_query_match_context,
    try_get_substruct_atom_matches_with_compiled_query_and_context,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo, ChiralTag, Hybridization};
use std::{collections::BTreeSet, sync::Arc};

/// Explicit detached source data; optional assignments retain source cache state.
#[derive(Clone, Copy)]
pub struct TautomerRecordView<'a> {
    pub topology: &'a TopologyBlock,
    pub coordinates: &'a CoordinateBlock,
    pub properties: &'a MoleculeProperties,
    pub valence: Option<&'a ValenceAssignment>,
    pub rings: Option<&'a RingInfo>,
}

/// Detached prepared candidate. It has no live cache or commit authority.
#[derive(Debug, Clone, PartialEq)]
pub struct TautomerRecord {
    pub topology: TopologyBlock,
    pub properties: MoleculeProperties,
    pub valence: ValenceAssignment,
    pub rings: RingInfo,
}
impl TautomerRecord {
    pub fn view<'a>(&'a self, coordinates: &'a CoordinateBlock) -> TautomerRecordView<'a> {
        TautomerRecordView {
            topology: &self.topology,
            properties: &self.properties,
            coordinates,
            valence: Some(&self.valence),
            rings: Some(&self.rings),
        }
    }
}

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum TautomerRunError {
    #[error(transparent)]
    Smiles(#[from] cosmolkit_smiles::SmilesParseError),
    #[error(transparent)]
    Valence(#[from] ValenceError),
    #[error(transparent)]
    Rings(#[from] RingFindingError),
    #[error(transparent)]
    Kekulize(#[from] KekulizeError),
    #[error(transparent)]
    Sanitize(#[from] SanitizeError),
    #[error(transparent)]
    LegacyStereo(#[from] LegacyStereoError),
    #[error(transparent)]
    QueryCompile(#[from] cosmolkit_search::QueryCompileError),
    #[error(transparent)]
    Match(#[from] cosmolkit_search::SubstructMatchError),
    #[error(transparent)]
    MatchContext(#[from] cosmolkit_search::QueryMatchContextError),
    #[error(transparent)]
    Topology(#[from] cosmolkit_model::TopologyValidationError),
    #[error(transparent)]
    AtomProperty(#[from] cosmolkit_model::AtomPropertyError),
    #[error(transparent)]
    MoleculeProperty(#[from] cosmolkit_model::MoleculePropertyError),
    #[error(transparent)]
    BondValue(#[from] cosmolkit_model::BondValueError),
    #[error(transparent)]
    PropertyString(#[from] PropertyStringError),
    #[error("tautomer match has {actual} atom rows, expected {expected}")]
    AtomMappingCount { expected: usize, actual: usize },
    #[error("tautomer match has {actual} bond rows, expected {expected}")]
    BondMappingCount { expected: usize, actual: usize },
    #[error("tautomer atom table has {actual} rows, expected {expected}")]
    AtomCountMismatch { expected: usize, actual: usize },
    #[error("tautomer bond table has {actual} rows, expected {expected}")]
    BondCountMismatch { expected: usize, actual: usize },
    #[error("tautomer match atom {atom} is outside {atom_count} atoms")]
    AtomOutOfRange { atom: usize, atom_count: usize },
    #[error("tautomer match bond {bond} is outside {bond_count} bonds")]
    BondOutOfRange { bond: usize, bond_count: usize },
    #[error("source transform has no donor/acceptor match endpoints")]
    EmptyTransformMatch,
    #[error("source transform has {actual} {field} edits but accesses {expected}")]
    EditCount {
        field: &'static str,
        expected: usize,
        actual: usize,
    },
    #[error("hydrogen count {count} on atom {atom} exceeds modeled storage")]
    HydrogenCountOutOfRange { atom: AtomId, count: u32 },
    #[error("charge {charge} on atom {atom} exceeds modeled storage")]
    FormalChargeOutOfRange { atom: AtomId, charge: i32 },
    #[error("query bond {bond} has no matched target bond")]
    MissingMappedBond { bond: usize },
    #[error("tautomer source has no hydrogen cache for atom {atom}")]
    MissingValence { atom: AtomId },
    #[error("no tautomer candidate can be selected canonically")]
    NoCanonicalTautomer,
    #[error("tautomer control invariant: {0}")]
    Control(String),
    #[error("tautomer callback/scorer failed: {0}")]
    Callback(String),
}

pub(crate) fn canonical_smiles(
    view: TautomerRecordView<'_>,
) -> Result<cosmolkit_model::PropertyText, TautomerRunError> {
    Ok(cosmolkit_smiles::write_smiles(
        cosmolkit_smiles::SmilesRecordView {
            topology: view.topology,
            coordinates: view.coordinates,
            properties: view.properties,
        },
    )?)
}
pub(crate) fn total_hydrogens(
    view: TautomerRecordView<'_>,
    atom: AtomId,
) -> Result<u32, TautomerRunError> {
    let empty = ValenceAssignment {
        explicit_valence: Vec::new(),
        implicit_hydrogens: Vec::new(),
    };
    let valence = match view.valence {
        Some(value) => value,
        None if view
            .topology
            .atoms
            .get(atom.index())
            .is_some_and(|a| a.no_implicit()) =>
        {
            &empty
        }
        None => return Err(ValenceError::ImplicitValenceCacheNotInitialized { atom }.into()),
    };
    Ok(total_hydrogen_count_from_validated(
        view.topology,
        valence,
        atom,
        false,
    )?)
}
pub(crate) fn prepared(view: TautomerRecordView<'_>) -> Result<TautomerRecord, TautomerRunError> {
    view.topology.validate()?;
    let valence = prepared_valence(view)?;
    let rings = match view.rings {
        Some(value) if value.is_symm_sssr() => value.clone(),
        _ => symmetrized_sssr(view.topology, &Default::default())?,
    };
    Ok(TautomerRecord {
        topology: view.topology.clone(),
        properties: view.properties.clone(),
        valence,
        rings,
    })
}
pub(crate) fn kekulized(candidate: &TautomerRecord) -> Result<TautomerRecord, TautomerRunError> {
    let assignment = kekulize_with_query_state_and_ring_info(
        &candidate.topology,
        &KekulizeParams {
            mark_atoms_bonds: false,
            canonical: true,
            max_backtracks: 100,
        },
        None,
        Some(&candidate.rings),
        Some(&candidate.valence),
    )?;
    let mut result = candidate.clone();
    result.topology = assignment.topology;
    if let Some(rings) = assignment.ring_update {
        result.rings = rings;
    }
    if let Some(valence) = assignment.final_valence {
        result.valence = valence;
    }

    Ok(result)
}
pub(crate) fn query_matches(
    view: TautomerRecordView<'_>,
    query: &CompiledQuery,
) -> Result<Vec<Vec<usize>>, TautomerRunError> {
    let owned_valence;
    let valence = match view.valence {
        Some(value) => value,
        None => {
            owned_valence = assign_valence(
                view.topology,
                &ValenceParams {
                    model: ValenceModel::RdkitLike,
                    strict: false,
                },
            )?;
            &owned_valence
        }
    };
    let owned_rings;
    let rings = match view.rings {
        Some(value) if value.is_initialized() => value,
        None | Some(_) => {
            owned_rings = symmetrized_sssr(view.topology, &Default::default())?;
            &owned_rings
        }
    };
    let target = SearchTarget::new(
        view.topology,
        view.coordinates,
        &view.topology.stereo_groups,
        Some(rings),
        Some(valence),
    );
    let context = build_prepared_query_match_context(view.topology, rings, valence)?;
    Ok(
        try_get_substruct_atom_matches_with_compiled_query_and_context(
            &target,
            query,
            &SubstructMatchParams::default(),
            &context,
        )?,
    )
}
pub(crate) fn transform_matches(
    candidate: &TautomerRecord,
    coordinates: &CoordinateBlock,
    transform: &TautomerTransform,
) -> Result<Vec<SubstructMatchResult>, TautomerRunError> {
    let rows = query_matches(candidate.view(coordinates), &transform.query)?;
    rows.into_iter()
        .map(|atom_mapping| {
            let mut bond_mapping = Vec::with_capacity(transform.query().num_bonds());
            for (index, bond) in transform.query().bonds().iter().enumerate() {
                let begin = atom_mapping[bond.begin().index()];
                let end = atom_mapping[bond.end().index()];
                let target = candidate
                    .topology
                    .adjacency
                    .neighbors_of(begin)
                    .iter()
                    .find(|n| n.atom_index == end)
                    .ok_or(TautomerRunError::MissingMappedBond { bond: index })?;
                bond_mapping.push(target.bond.index());
            }
            Ok(SubstructMatchResult {
                atom_mapping,
                bond_mapping,
            })
        })
        .collect()
}
pub(crate) fn assign_stereo(candidate: &mut TautomerRecord) -> Result<(), TautomerRunError> {
    let assignment = assign_legacy_stereochemistry_with_assignments(
        candidate.topology.clone(),
        &candidate.valence,
        &candidate.rings,
        true,
        false,
    )?;
    candidate.topology = assignment.topology;
    if let Some(rings) = assignment.ring_update {
        candidate.rings = rings;
    }
    // RDKit❗✔️: mol.setProp(common_properties::_StereochemDone, 1, true);
    candidate
        .properties
        .set_computed_prop("_StereochemDone", 1_i32)?;
    Ok(())
}
fn is_stereo_beyond_any(stereo: BondStereo) -> bool {
    matches!(
        stereo,
        BondStereo::Z
            | BondStereo::E
            | BondStereo::Cis
            | BondStereo::Trans
            | BondStereo::AtropCw
            | BondStereo::AtropCcw
    )
}

fn tautomer_bond_is_ring(
    molecule: &mut TautomerRecord,
    bond_id: BondId,
) -> Result<bool, RingFindingError> {
    // RDKit✔️✔️:       if (tautBond->getBondType() == Bond::DOUBLE) {
    // RDKit✔️✔️:         auto ringInfo = taut.getRingInfo();
    // RDKit✔️✔️:         if (!ringInfo || !ringInfo->isFindFastOrBetter()) {
    // RDKit✔️✔️:           MolOps::fastFindRings(taut);
    // RDKit✔️✔️:           ringInfo = taut.getRingInfo();
    // RDKit✔️✔️:         }
    // Source-required O(V+E) preparation is retained in the detached candidate;
    // every later query reads the same result, without another ring search.
    let bond = &molecule.topology.bonds[bond_id.index()];
    if bond.order() != BondOrder::Double {
        return Ok(false);
    }
    if !molecule.rings.is_find_fast_or_better() {
        molecule.rings = fast_find_rings(&molecule.topology)?;
    }
    let rings = &molecule.rings;
    Ok(rings.num_bond_rings(bond_id) != 0
        || (rings.num_atom_rings(bond.begin()) != 0 && rings.num_atom_rings(bond.end()) != 0))
}

pub(crate) fn set_tautomer_stereo_and_isotopic_hydrogens(
    source: &TautomerRecord,
    tautomer: &mut TautomerRecord,
    modified_atoms: &BTreeSet<AtomId>,
    modified_bonds: &BTreeSet<BondId>,
    options: TautomerParams,
    coordinates: &CoordinateBlock,
) -> Result<bool, TautomerRunError> {
    // RDKit✔️❌: bool TautomerEnumerator::setTautomerStereoAndIsoHs(
    // RDKit✔️❌:     const ROMol &mol, ROMol &taut, const TautomerEnumeratorResult &res) const {
    // RDKit✔️❌:   bool modified = false;
    // The source transition is reproduced below and validated through the
    // complete tautomer stereo matrix. BTreeSet membership is O(log n), versus
    // the source dynamic_bitset's O(1), so the complexity axis is intentionally
    // not claimed equivalent.
    if tautomer.topology.atoms.len() != source.topology.atoms.len() {
        return Err(TautomerRunError::AtomCountMismatch {
            expected: source.topology.atoms.len(),
            actual: tautomer.topology.atoms.len(),
        });
    }
    if tautomer.topology.bonds.len() != source.topology.bonds.len() {
        return Err(TautomerRunError::BondCountMismatch {
            expected: source.topology.bonds.len(),
            actual: tautomer.topology.bonds.len(),
        });
    }
    let mut modified = false;

    // RDKit✔️❌:   for (auto atom : mol.atoms()) {
    // RDKit✔️❌:     auto atomIdx = atom->getIdx();
    // RDKit✔️❌:     if (!res.d_modifiedAtoms.test(atomIdx)) {
    // RDKit✔️❌:       continue;
    // RDKit✔️❌:     }
    for atom_id in source.topology.atoms.iter().map(cosmolkit_model::Atom::id) {
        if !modified_atoms.contains(&atom_id) {
            continue;
        }
        let source_atom = &source.topology.atoms[atom_id.index()];
        let clear_isotopic_hydrogens = {
            let tautomer_atom = &tautomer.topology.atoms[atom_id.index()];
            !tautomer_atom.tracked_isotopic_hydrogens().is_empty()
                && (options.remove_isotopic_hydrogens()
                    || total_hydrogens(tautomer.view(coordinates), atom_id)? == 0)
        };
        let tautomer_atom = &mut tautomer.topology.atoms[atom_id.index()];

        // RDKit✔️❌:     auto tautAtom = taut.getAtomWithIdx(atomIdx);
        // RDKit✔️❌:     // clear chiral tag on sp2 atoms (also sp3 if d_removeSp3Stereo is true)
        // RDKit✔️❌:     if (tautAtom->getHybridization() == Atom::SP2 || d_removeSp3Stereo) {
        if tautomer_atom.hybridization() == Hybridization::Sp2 || options.remove_sp3_stereo() {
            // RDKit✔️❌:       modified |= (tautAtom->getChiralTag() != Atom::CHI_UNSPECIFIED);
            // RDKit✔️❌:       tautAtom->setChiralTag(Atom::CHI_UNSPECIFIED);
            modified |= tautomer_atom.chiral_tag() != ChiralTag::Unspecified;
            tautomer_atom.set_chiral_tag(ChiralTag::Unspecified);
            // RDKit✔️❌:       if (tautAtom->hasProp(common_properties::_CIPCode)) {
            // RDKit✔️❌:         tautAtom->clearProp(common_properties::_CIPCode);
            // RDKit✔️❌:       }
            if tautomer_atom.prop("_CIPCode").is_some() {
                tautomer_atom.clear_prop("_CIPCode")?;
            }
        } else {
            // RDKit✔️❌:     } else {
            // RDKit✔️❌:       modified |= (tautAtom->getChiralTag() != atom->getChiralTag());
            // RDKit✔️❌:       tautAtom->setChiralTag(atom->getChiralTag());
            modified |= tautomer_atom.chiral_tag() != source_atom.chiral_tag();
            tautomer_atom.set_chiral_tag(source_atom.chiral_tag());
            // RDKit✔️❌:       if (atom->hasProp(common_properties::_CIPCode)) {
            // RDKit✔️❌:         tautAtom->setProp(
            // RDKit✔️❌:             common_properties::_CIPCode,
            // RDKit✔️❌:             atom->getProp<std::string>(common_properties::_CIPCode));
            // RDKit✔️❌:       }
            if let Some(cip_code) = source_atom.prop("_CIPCode") {
                tautomer_atom.set_prop("_CIPCode", property_value_to_string(cip_code)?)?;
            }
        }
        // RDKit✔️❌:     // remove isotopic Hs if present (and if d_removeIsotopicHs is true)
        // RDKit✔️❌:     if (tautAtom->hasProp(common_properties::_isotopicHs) &&
        // RDKit✔️❌:         (d_removeIsotopicHs || !tautAtom->getTotalNumHs())) {
        // RDKit✔️❌:       tautAtom->clearProp(common_properties::_isotopicHs);
        // RDKit✔️❌:     }
        if clear_isotopic_hydrogens {
            tautomer_atom.set_tracked_isotopic_hydrogens(Vec::new());
        }
        // RDKit✔️❌:   }
    }

    // RDKit✔️❌:   // remove stereochemistry on bonds that are part of a tautomeric path
    // RDKit✔️❌:   for (auto bond : mol.bonds()) {
    // RDKit✔️❌:     auto bondIdx = bond->getIdx();
    // RDKit✔️❌:     if (!res.d_modifiedBonds.test(bondIdx)) {
    // RDKit✔️❌:       continue;
    // RDKit✔️❌:     }
    for bond_id in source.topology.bonds.iter().map(cosmolkit_model::Bond::id) {
        if !modified_bonds.contains(&bond_id) {
            continue;
        }
        let source_bond = &source.topology.bonds[bond_id.index()];
        // RDKit✔️❌:     std::vector<unsigned int> bondsToClearDirs;
        let mut bonds_to_clear_directions = Vec::new();
        // RDKit✔️❌:     if (bond->getBondType() == Bond::DOUBLE &&
        // RDKit✔️❌:         bond->getStereo() > Bond::STEREOANY) {
        if source_bond.order() == BondOrder::Double && is_stereo_beyond_any(source_bond.stereo()) {
            // RDKit✔️❌:       for (auto atom : {bond->getBeginAtom(), bond->getEndAtom()}) {
            // RDKit✔️❌:         for (const auto &nbri :
            // RDKit✔️❌:              boost::make_iterator_range(mol.getAtomBonds(atom))) {
            for atom_id in [source_bond.begin(), source_bond.end()] {
                for neighbor in source.topology.adjacency.neighbors_of(atom_id.index()) {
                    let adjacent_bond = &source.topology.bonds[neighbor.bond.index()];
                    // RDKit✔️❌:           const auto &obnd = mol[nbri];
                    // RDKit✔️❌:           if (obnd->getBondDir() == Bond::ENDDOWNRIGHT ||
                    // RDKit✔️❌:               obnd->getBondDir() == Bond::ENDUPRIGHT) {
                    // RDKit✔️❌:             bondsToClearDirs.push_back(obnd->getIdx());
                    // RDKit✔️❌:           }
                    if matches!(
                        adjacent_bond.direction(),
                        BondDirection::EndDownRight | BondDirection::EndUpRight
                    ) {
                        bonds_to_clear_directions.push(adjacent_bond.id());
                    }
                    // RDKit✔️❌:         }
                }
                // RDKit✔️❌:       }
            }
        }

        let tautomer_bond = &tautomer.topology.bonds[bond_id.index()];
        let candidate_order = tautomer_bond.order();
        let candidate_stereo = tautomer_bond.stereo();
        let candidate_stereo_atoms = tautomer_bond.stereo_atoms();
        let remove_stereo =
            tautomer_bond.order() != BondOrder::Double || options.remove_bond_stereo();
        let target_stereo = if remove_stereo {
            let is_ring_bond = tautomer_bond_is_ring(tautomer, bond_id)?;
            if candidate_order == BondOrder::Double && !is_ring_bond {
                BondStereo::Any
            } else {
                BondStereo::None
            }
        } else {
            source_bond.stereo()
        };
        let target_stereo_atoms = if remove_stereo {
            None
        } else {
            source_bond.stereo_atoms().or(candidate_stereo_atoms)
        };
        modified |= candidate_stereo != target_stereo
            || (!remove_stereo
                && candidate_stereo_atoms.is_some() != source_bond.stereo_atoms().is_some());

        // RDKit✔️❌:     auto tautBond = taut.getBondWithIdx(bondIdx);
        // RDKit✔️❌:     if (tautBond->getBondType() != Bond::DOUBLE || d_removeBondStereo) {
        // RDKit✔️❌:       tautBond->setStereo(targetStereo);
        // RDKit✔️❌:       tautBond->getStereoAtoms().clear();
        // RDKit✔️❌:     } else {
        // RDKit✔️❌:       const INT_VECT &sa = bond->getStereoAtoms();
        // RDKit✔️❌:       if (sa.size() == 2) {
        // RDKit✔️❌:         tautBond->setStereoAtoms(sa.front(), sa.back());
        // RDKit✔️❌:       }
        // RDKit✔️❌:       tautBond->setStereo(bond->getStereo());
        // RDKit✔️❌:     }
        let tautomer_bond = &mut tautomer.topology.bonds[bond_id.index()];
        tautomer_bond.set_stereo_atoms(target_stereo_atoms);
        tautomer_bond.set_stereo(target_stereo)?;
        for adjacent_bond_id in bonds_to_clear_directions {
            // RDKit✔️❌:       for (auto bi : bondsToClearDirs) {
            // RDKit✔️❌:         taut.getBondWithIdx(bi)->setBondDir(
            // RDKit✔️❌:             mol.getBondWithIdx(bi)->getBondDir());
            // RDKit✔️❌:       }
            let direction = if remove_stereo {
                BondDirection::None
            } else {
                source.topology.bonds[adjacent_bond_id.index()].direction()
            };
            tautomer.topology.bonds[adjacent_bond_id.index()].set_direction(direction);
        }
        // RDKit✔️❌:   }
    }

    // RDKit✔️❌:   if (d_reassignStereo) {
    if options.reassign_stereo() {
        // RDKit✔️❌:     static const bool cleanIt = true;
        // RDKit✔️❌:     static const bool force = true;
        // RDKit✔️❌:     MolOps::assignStereochemistry(taut, cleanIt, force);
        assign_stereo(tautomer)?;

        // RDKit✔️❌:     if (d_removeBondStereo) {
        if options.remove_bond_stereo() {
            // RDKit✔️✔️:       auto ringInfo = taut.getRingInfo();
            // RDKit✔️✔️:       if (!ringInfo || !ringInfo->isFindFastOrBetter()) {
            // RDKit✔️✔️:         MolOps::fastFindRings(taut);
            // RDKit✔️✔️:         ringInfo = taut.getRingInfo();
            // RDKit✔️✔️:       }
            if !tautomer.rings.is_find_fast_or_better() {
                tautomer.rings = fast_find_rings(&tautomer.topology)?;
            }
            for bond_id in modified_bonds.iter().copied() {
                let bond = &tautomer.topology.bonds[bond_id.index()];
                let target_stereo = if bond.order() != BondOrder::Double {
                    BondStereo::None
                } else if tautomer_bond_is_ring(tautomer, bond_id)? {
                    BondStereo::None
                } else {
                    BondStereo::Any
                };
                // RDKit✔️❌:       for (auto bond : taut.bonds()) {
                // RDKit✔️❌:         const auto bondIdx = bond->getIdx();
                // RDKit✔️❌:         if (!res.d_modifiedBonds.test(bondIdx)) {
                // RDKit✔️❌:           continue;
                // RDKit✔️❌:         }
                // RDKit✔️❌:         if (bond->getBondType() != Bond::DOUBLE) {
                // RDKit✔️❌:           bond->setStereo(Bond::STEREONONE);
                // RDKit✔️❌:           bond->getStereoAtoms().clear();
                // RDKit✔️❌:           continue;
                // RDKit✔️❌:         }
                // RDKit✔️❌:         bond->setStereo(isRingBond ? Bond::STEREONONE : Bond::STEREOANY);
                // RDKit✔️❌:         bond->getStereoAtoms().clear();
                let bond = &mut tautomer.topology.bonds[bond_id.index()];
                bond.set_stereo_atoms(None);
                bond.set_stereo(target_stereo)?;
                // RDKit✔️❌:       }
            }
            // RDKit✔️❌:     }
        }
    } else {
        // RDKit✔️❌:   } else {
        // RDKit✔️❌:     taut.setProp(common_properties::_StereochemDone, 1);
        tautomer.properties.set_prop("_StereochemDone", 1_i32)?;
        // RDKit✔️❌:   }
    }
    // RDKit✔️❌:   return modified;
    // RDKit✔️❌: }
    Ok(modified)
}
pub(crate) fn apply_tautomer_transform_match(
    source: &TautomerRecord,
    candidate: &TautomerRecord,
    coordinates: &CoordinateBlock,
    transform: &TautomerTransform,
    match_result: &SubstructMatchResult,
    current_modified_atoms: &BTreeSet<AtomId>,
    current_modified_bonds: &BTreeSet<BondId>,
    contains_smiles: &dyn Fn(&cosmolkit_model::PropertyText) -> bool,
    options: TautomerParams,
) -> Result<TautomerExpansionAttempt<Arc<TautomerRecord>>, TautomerRunError> {
    // RDKit✔️❌:           RWMOL_SPTR product(new RWMol(*kmol));
    // Only topology is copied for editing. Unchanged coordinates stay borrowed;
    // prepared ring/H assignments remain detached source values.
    let mut topology = candidate.topology.clone();
    if match_result.atom_mapping.len() != transform.query().num_atoms() {
        return Err(TautomerRunError::AtomMappingCount {
            expected: transform.query().num_atoms(),
            actual: match_result.atom_mapping.len(),
        });
    }
    if match_result.bond_mapping.len() != transform.query().num_bonds() {
        return Err(TautomerRunError::BondMappingCount {
            expected: transform.query().num_bonds(),
            actual: match_result.bond_mapping.len(),
        });
    }
    for &atom in &match_result.atom_mapping {
        if atom >= topology.atoms.len() {
            return Err(TautomerRunError::AtomOutOfRange {
                atom,
                atom_count: topology.atoms.len(),
            });
        }
    }
    for &bond in &match_result.bond_mapping {
        if bond >= topology.bonds.len() {
            return Err(TautomerRunError::BondOutOfRange {
                bond,
                bond_count: topology.bonds.len(),
            });
        }
    }

    if match_result.atom_mapping.is_empty() {
        return Err(TautomerRunError::EmptyTransformMatch);
    }
    for (field, actual, expected) in [
        (
            "bond",
            transform.bond_types().len(),
            transform.query().num_bonds(),
        ),
        (
            "charge",
            transform.charges().len(),
            transform.query().num_atoms(),
        ),
    ] {
        if actual != 0 && actual < expected {
            return Err(TautomerRunError::EditCount {
                field,
                actual,
                expected,
            });
        }
    }
    // RDKit✔️❌:           // Remove a hydrogen from the first matched atom and add one to the
    // RDKit✔️❌:           // last
    // RDKit✔️❌:           int firstIdx = match.front().second;
    // RDKit✔️❌:           int lastIdx = match.back().second;
    // RDKit✔️❌:           Atom *first = product->getAtomWithIdx(firstIdx);
    // RDKit✔️❌:           Atom *last = product->getAtomWithIdx(lastIdx);
    // RDKit✔️❌:           res.d_modifiedAtoms.set(firstIdx);
    // RDKit✔️❌:           res.d_modifiedAtoms.set(lastIdx);
    let first = AtomId::new(match_result.atom_mapping[0]);
    let last = AtomId::new(match_result.atom_mapping[match_result.atom_mapping.len() - 1]);
    let mut modified_atoms = current_modified_atoms.clone();
    let mut modified_bonds = current_modified_bonds.clone();
    modified_atoms.insert(first);
    modified_atoms.insert(last);

    // RDKit✔️❌:           first->setNumExplicitHs(
    // RDKit✔️❌:               std::max(0, static_cast<int>(first->getTotalNumHs()) - 1));
    // RDKit✔️❌:           last->setNumExplicitHs(last->getTotalNumHs() + 1);
    let first_hydrogens = total_hydrogens(candidate.view(coordinates), first)?;
    let first_explicit = u8::try_from(first_hydrogens.saturating_sub(1)).map_err(|_| {
        TautomerRunError::HydrogenCountOutOfRange {
            atom: first,
            count: first_hydrogens.saturating_sub(1),
        }
    })?;
    topology.atoms[first.index()].set_explicit_hydrogens(first_explicit);
    // Source reads the acceptor after updating the donor, even when a
    // custom one-atom query maps both endpoints to the same atom. The cached
    // implicit-H row remains unchanged until the subsequent sanitize stage.
    let last_hydrogens =
        total_hydrogen_count_from_validated(&topology, &candidate.valence, last, false)?;
    let last_total =
        last_hydrogens
            .checked_add(1)
            .ok_or(TautomerRunError::HydrogenCountOutOfRange {
                atom: last,
                count: u32::MAX,
            })?;
    let last_explicit =
        u8::try_from(last_total).map_err(|_| TautomerRunError::HydrogenCountOutOfRange {
            atom: last,
            count: last_total,
        })?;
    topology.atoms[last.index()].set_explicit_hydrogens(last_explicit);

    // RDKit✔️❌:           // Remove any implicit hydrogens from the first and last atoms
    // RDKit✔️❌:           // now we have set the count explicitly
    // RDKit✔️❌:           first->setNoImplicit(true);
    // RDKit✔️❌:           last->setNoImplicit(true);
    topology.atoms[first.index()].set_no_implicit(true);
    topology.atoms[last.index()].set_no_implicit(true);

    // RDKit✔️❌:           // Adjust bond orders
    // RDKit✔️❌:           unsigned int bi = 0;
    // RDKit✔️❌:           for (size_t i = 0; i < transform.Mol->getNumBonds(); ++i) {
    // RDKit✔️❌:             const auto tbond = transform.Mol->getBondWithIdx(i);
    // RDKit✔️❌:             Bond *bond = product->getBondBetweenAtoms(
    // RDKit✔️❌:                 match[tbond->getBeginAtomIdx()].second,
    // RDKit✔️❌:                 match[tbond->getEndAtomIdx()].second);
    // RDKit✔️❌:             ASSERT_INVARIANT(bond, "required bond not found");
    for (query_bond, &target_bond) in match_result.bond_mapping.iter().enumerate() {
        let bond_id = BondId::new(target_bond);
        let bond = &mut topology.bonds[target_bond];
        // RDKit✔️❌:             // check if bonds is specified in tautomer.in file
        // RDKit✔️❌:             if (!transform.BondTypes.empty()) {
        // RDKit✔️❌:               bond->setBondType(transform.BondTypes[bi]);
        // RDKit✔️❌:               ++bi;
        // RDKit✔️❌:             } else {
        if let Some(&bond_type) = transform.bond_types().get(query_bond) {
            bond.set_order(bond_type);
        } else {
            // RDKit✔️❌:               Bond::BondType bondtype = bond->getBondType();
            // RDKit✔️❌:               if (bondtype == Bond::SINGLE) {
            // RDKit✔️❌:                 bond->setBondType(Bond::DOUBLE);
            // RDKit✔️❌:               }
            // RDKit✔️❌:               if (bondtype == Bond::DOUBLE) {
            // RDKit✔️❌:                 bond->setBondType(Bond::SINGLE);
            // RDKit✔️❌:               }
            let bond_type = bond.order();
            if bond_type == BondOrder::Single {
                bond.set_order(BondOrder::Double);
            }
            if bond_type == BondOrder::Double {
                bond.set_order(BondOrder::Single);
            }
        }
        // RDKit✔️❌:             }
        // RDKit✔️❌:             res.d_modifiedBonds.set(bond->getIdx());
        modified_bonds.insert(bond_id);
        // RDKit✔️❌:           }
    }

    // RDKit✔️❌:           // TODO adjust charges
    // RDKit✔️❌:           if (!transform.Charges.empty()) {
    // RDKit✔️❌:             unsigned int ci = 0;
    // RDKit✔️❌:             for (const auto &pair : match) {
    // RDKit✔️❌:               Atom *atom = product->getAtomWithIdx(pair.second);
    // RDKit✔️❌:               atom->setFormalCharge(atom->getFormalCharge() +
    // RDKit✔️❌:                                     transform.Charges[ci++]);
    // RDKit✔️❌:             }
    // RDKit✔️❌:           }
    for (&target_atom, &delta) in match_result.atom_mapping.iter().zip(transform.charges()) {
        let atom_id = AtomId::new(target_atom);
        let charge = i32::from(topology.atoms[target_atom].formal_charge()) + delta;
        let charge =
            i8::try_from(charge).map_err(|_| TautomerRunError::FormalChargeOutOfRange {
                atom: atom_id,
                charge,
            })?;
        topology.atoms[target_atom].set_formal_charge(charge);
    }

    let sanitize_ops = SanitizeOperations::KEKULIZE
        | SanitizeOperations::SET_AROMATICITY
        | SanitizeOperations::SET_CONJUGATION
        | SanitizeOperations::SET_HYBRIDIZATION
        | SanitizeOperations::ADJUST_HS;
    // RDKit✔️❌:           unsigned int failedOp;
    // RDKit✔️❌:           try {
    // RDKit✔️❌:             MolOps::sanitizeMol(*product, failedOp,
    // RDKit✔️❌:                                 MolOps::SANITIZE_KEKULIZE |
    // RDKit✔️❌:                                     MolOps::SANITIZE_SETAROMATICITY |
    // RDKit✔️❌:                                     MolOps::SANITIZE_SETCONJUGATION |
    // RDKit✔️❌:                                     MolOps::SANITIZE_SETHYBRIDIZATION |
    // RDKit✔️❌:                                     MolOps::SANITIZE_ADJUSTHS);
    // RDKit✔️❌:           } catch (const KekulizeException &) {
    // RDKit✔️❌:             continue;
    // RDKit✔️❌:           }
    // BEGIN RDKIT CPP FUNCTION MolOps::sanitizeMol entry property clearing
    // RDKit✔️❌:   // clear out any cached properties
    // RDKit✔️❌:   mol.clearComputedProps();
    // END RDKIT CPP FUNCTION MolOps::sanitizeMol entry property clearing
    // The topology-only sanitizer owns atom/bond clearing. Its detached
    // molecule-property companion must run first, as ROMol clears RDProps
    // before the atom and bond loops. Keep one existing property-block clone;
    // moving it here adds no scan/allocation or chemistry fallback.
    let mut properties = candidate.properties.clone();
    properties.clear_computed_props()?;
    let assignment = match sanitize_topology(
        &topology,
        &SanitizeParams {
            operations: sanitize_ops,
        },
    ) {
        Ok(value) => value,
        Err(SanitizeError::Kekulize {
            source: KekulizeError::NotKekulizable { .. },
            ..
        }) => {
            return Ok(TautomerExpansionAttempt::RecoverableKekulizeFailure {
                modified_atoms,
                modified_bonds,
            });
        }
        Err(error) => return Err(error.into()),
    };
    let valence = assignment
        .final_valence
        .or(assignment.non_strict_valence)
        .expect("sanitize always produces source property-cache rows");
    let rings = assignment
        .final_rings
        .expect("SET_AROMATICITY initializes the source ring state");
    // RDKit❗✔️: int narom = 0;
    // RDKit❗✔️: mol.setProp(common_properties::numArom, narom, true);
    if let Some(count) = assignment.aromatic_ring_count {
        properties.set_computed_prop(
            "numArom",
            i32::try_from(count).map_err(|_| cosmolkit_core::SanitizeError::Aromaticity {
                stage: cosmolkit_core::SanitizeStage::SetAromaticity,
                source: cosmolkit_core::AromaticityError::IntegerOverflow {
                    field: "source numArom int",
                },
            })?,
        )?;
    }
    let mut product = TautomerRecord {
        topology: assignment.topology,
        properties,
        valence,
        rings,
    };
    // RDKit✔️❌:           setTautomerStereoAndIsoHs(mol, *product, res);
    set_tautomer_stereo_and_isotopic_hydrogens(
        source,
        &mut product,
        &modified_atoms,
        &modified_bonds,
        options,
        coordinates,
    )?;
    // RDKit✔️❌:           tsmiles = MolToSmiles(*product, true);
    let canonical_smiles = canonical_smiles(product.view(coordinates))?;

    // RDKit✔️❌:           if (res.d_tautomers.find(tsmiles) != res.d_tautomers.end()) {
    // RDKit✔️❌:             continue;
    // RDKit✔️❌:           }
    if contains_smiles(&canonical_smiles) {
        return Ok(TautomerExpansionAttempt::Duplicate {
            canonical_smiles,
            modified_atoms,
            modified_bonds,
        });
    }

    // RDKit✔️❌:           // in addition to the above transformations, sanitization may modify
    // RDKit✔️❌:           // bonds, e.g. Cc1nc2ccccc2[nH]1
    // RDKit✔️❌:           for (size_t i = 0; i < mol.getNumBonds(); i++) {
    // RDKit✔️❌:             auto molBondType = mol.getBondWithIdx(i)->getBondType();
    // RDKit✔️❌:             auto tautBondType = product->getBondWithIdx(i)->getBondType();
    // RDKit✔️❌:             if (molBondType != tautBondType && !res.d_modifiedBonds.test(i)) {
    // RDKit✔️❌:               res.d_modifiedBonds.set(i);
    // RDKit✔️❌:             }
    // RDKit✔️❌:           }
    for (index, (source_bond, product_bond)) in source
        .topology
        .bonds
        .iter()
        .zip(&product.topology.bonds)
        .enumerate()
    {
        if source_bond.order() != product_bond.order() {
            modified_bonds.insert(BondId::new(index));
        }
    }

    // RDKit✔️❌:           RWMOL_SPTR kekulized_product(new RWMol(*product));
    // RDKit✔️❌:           // canonical=true for order-independent tautomer deduplication
    // RDKit✔️❌:           MolOps::Kekulize(*kekulized_product, false, true);
    let kekulized_product = kekulized(&product)?;
    // RDKit✔️❌:           res.d_tautomers[tsmiles] = Tautomer(
    // RDKit✔️❌:               std::move(product), std::move(kekulized_product),
    // RDKit✔️❌:               res.d_modifiedAtoms.count(), res.d_modifiedBonds.count());
    Ok(TautomerExpansionAttempt::Product(TautomerExpandedProduct {
        tautomer: Arc::new(product),
        kekulized: Arc::new(kekulized_product),
        canonical_smiles,
        modified_atoms,
        modified_bonds,
    }))
}

#[cfg(test)]
#[path = "stereo_tests.rs"]
pub(crate) mod stereo_tests;

fn prepared_valence(view: TautomerRecordView<'_>) -> Result<ValenceAssignment, TautomerRunError> {
    Ok(match view.valence {
        Some(value)
            if value.explicit_valence.len() == view.topology.atoms.len()
                && value.implicit_hydrogens.len() == view.topology.atoms.len()
                && value.explicit_valence.iter().all(|v| *v >= 0)
                && view
                    .topology
                    .atoms
                    .iter()
                    .zip(&value.implicit_hydrogens)
                    .all(|(a, v)| a.no_implicit() || *v >= 0) =>
        {
            value.clone()
        }
        _ => assign_valence(
            view.topology,
            &ValenceParams {
                model: ValenceModel::RdkitLike,
                strict: false,
            },
        )?,
    })
}

pub(crate) fn copy_for_canonical_assignment(
    view: TautomerRecordView<'_>,
) -> Result<TautomerRecord, TautomerRunError> {
    // RDKit✔️❌:   ROMol *res = new ROMol(*bestMol);
    // Preserve existing ring quality before the unique legacy owner performs
    // source-required fast-ring and cleanIt SymmSSSR preparation.
    view.topology.validate()?;
    Ok(TautomerRecord {
        topology: view.topology.clone(),
        properties: view.properties.clone(),
        valence: prepared_valence(view)?,
        rings: match view.rings {
            Some(rings) => rings.clone(),
            None => RingInfo::new(
                RingFindType::OtherOrUnknown,
                view.topology.atoms.len(),
                view.topology.bonds.len(),
            ),
        },
    })
}

#[cfg(test)]
#[path = "application_tests.rs"]
mod application_tests;
