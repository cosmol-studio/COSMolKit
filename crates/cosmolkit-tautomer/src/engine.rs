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
    pub coordinates: Option<CoordinateBlock>,
}
impl TautomerRecord {
    pub fn view<'a>(&'a self, coordinates: &'a CoordinateBlock) -> TautomerRecordView<'a> {
        TautomerRecordView {
            topology: &self.topology,
            properties: &self.properties,
            coordinates: self.coordinates.as_ref().unwrap_or(coordinates),
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
    Coordinates(#[from] cosmolkit_model::CoordinateValidationError),
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
    #[error("source ring-cache commit failed; owning host retains the typed cause")]
    SourceCacheCommit,
    #[error("tautomer score snapshot does not match its cache-loan target")]
    ScoreCacheTargetMismatch,
    #[error("tautomer control invariant: {0}")]
    Control(String),
    #[error("tautomer callback/scorer failed: {0}")]
    Callback(String),
}

#[doc(hidden)]
pub fn canonical_smiles(
    view: TautomerRecordView<'_>,
) -> Result<cosmolkit_model::PropertyText, TautomerRunError> {
    Ok(cosmolkit_smiles::write_smiles(
        cosmolkit_smiles::SmilesRecordView {
            topology: view.topology,
            coordinates: view.coordinates,
            properties: view.properties,
            rings: view.rings,
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
        coordinates: None,
        valence,
        rings,
    })
}
pub(crate) fn kekulized(candidate: &TautomerRecord) -> Result<TautomerRecord, TautomerRunError> {
    let mut result = candidate.clone();
    cosmolkit_core::source_kekulize_attempt(
        &mut result.topology,
        &mut result.valence,
        &mut result.rings,
        &KekulizeParams {
            mark_atoms_bonds: false,
            canonical: true,
            max_backtracks: 100,
        },
    )?;
    Ok(result)
}
pub(crate) fn query_matches(
    view: TautomerRecordView<'_>,
    query: &CompiledQuery,
) -> Result<Vec<Vec<usize>>, TautomerRunError> {
    query_matches_with_params(view, query, &Default::default())
}

pub(crate) fn query_matches_with_params(
    view: TautomerRecordView<'_>,
    query: &CompiledQuery,
    params: &SubstructMatchParams,
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
            &target, query, params, &context,
        )?,
    )
}

pub(crate) fn query_match_count(
    view: TautomerRecordView<'_>,
    query: &CompiledQuery,
    params: &SubstructMatchParams,
) -> Result<u32, TautomerRunError> {
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
        cosmolkit_search::try_get_substruct_match_count_with_compiled_query_and_context(
            &target, query, params, &context,
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
    for (atom, facts) in assignment.atom_valence_updates {
        candidate.valence.explicit_valence[atom.index()] = i32::from(facts.explicit_valence);
        candidate.valence.implicit_hydrogens[atom.index()] = i32::from(facts.implicit_valence);
    }
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
    // BEGIN RDKIT CPP FUNCTION RDKit::MolStandardize::TautomerEnumerator::setTautomerStereoAndIsoHs
    // RDKit❗❌: bool TautomerEnumerator::setTautomerStereoAndIsoHs(
    // RDKit❗❌:     const ROMol &mol, ROMol &taut, const TautomerEnumeratorResult &res) const {
    // RDKit❗❌:   bool modified = false;
    // RDKit❗❌:   // Iterate only the atoms/bonds actually modified by transforms.
    // RDKit❗❌:   for (auto atomIdx = res.d_modifiedAtoms.find_first();
    // RDKit❗❌:        atomIdx != boost::dynamic_bitset<>::npos;
    // RDKit❗❌:        atomIdx = res.d_modifiedAtoms.find_next(atomIdx)) {
    // RDKit❗❌:     const auto atom = mol.getAtomWithIdx(static_cast<unsigned int>(atomIdx));
    // RDKit❗❌:     auto tautAtom = taut.getAtomWithIdx(atomIdx);
    // RDKit❗❌:     // clear chiral tag on sp2 atoms (also sp3 if d_removeSp3Stereo is true)
    // RDKit❗❌:     if (tautAtom->getHybridization() == Atom::SP2 || d_removeSp3Stereo) {
    // RDKit❗❌:       modified |= (tautAtom->getChiralTag() != Atom::CHI_UNSPECIFIED);
    // RDKit❗❌:       tautAtom->setChiralTag(Atom::CHI_UNSPECIFIED);
    // RDKit❗❌:       if (tautAtom->hasProp(common_properties::_CIPCode)) {
    // RDKit❗❌:         tautAtom->clearProp(common_properties::_CIPCode);
    // RDKit❗❌:       }
    // RDKit❗❌:     } else {
    // RDKit❗❌:       modified |= (tautAtom->getChiralTag() != atom->getChiralTag());
    // RDKit❗❌:       tautAtom->setChiralTag(atom->getChiralTag());
    // RDKit❗❌:       if (atom->hasProp(common_properties::_CIPCode)) {
    // RDKit❗❌:         tautAtom->setProp(
    // RDKit❗❌:             common_properties::_CIPCode,
    // RDKit❗❌:             atom->getProp<std::string>(common_properties::_CIPCode));
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     // remove isotopic Hs if present (and if d_removeIsotopicHs is true)
    // RDKit❗❌:     if (tautAtom->hasProp(common_properties::_isotopicHs) &&
    // RDKit❗❌:         (d_removeIsotopicHs || !tautAtom->getTotalNumHs())) {
    // RDKit❗❌:       tautAtom->clearProp(common_properties::_isotopicHs);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   // remove stereochemistry on bonds that are part of a tautomeric path
    // RDKit❗❌:   // Build bond lookup caches to avoid O(n) getBondWithIdx calls.
    // RDKit❗❌:   // getBondWithIdx iterates through all bonds to find the one with the given
    // RDKit❗❌:   // index, which is expensive when called multiple times in a loop.
    // RDKit❗❌:   std::vector<const Bond *> molBonds;
    // RDKit❗❌:   std::vector<Bond *> tautBonds;
    // RDKit❗❌:   if (res.d_modifiedBonds.any()) {
    // RDKit❗❌:     const auto numBonds = mol.getNumBonds();
    // RDKit❗❌:     molBonds.resize(numBonds);
    // RDKit❗❌:     tautBonds.resize(numBonds);
    // RDKit❗❌:     for (auto bond : mol.bonds()) {
    // RDKit❗❌:       molBonds[bond->getIdx()] = bond;
    // RDKit❗❌:     }
    // RDKit❗❌:     for (auto bond : taut.bonds()) {
    // RDKit❗❌:       tautBonds[bond->getIdx()] = bond;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   for (auto bondIdx = res.d_modifiedBonds.find_first();
    // RDKit❗❌:        bondIdx != boost::dynamic_bitset<>::npos;
    // RDKit❗❌:        bondIdx = res.d_modifiedBonds.find_next(bondIdx)) {
    // RDKit❗❌:     const auto bond = molBonds[bondIdx];
    // RDKit❗❌:     std::vector<unsigned int> bondsToClearDirs;
    // RDKit❗❌:     if (bond->getBondType() == Bond::DOUBLE &&
    // RDKit❗❌:         bond->getStereo() > Bond::STEREOANY) {
    // RDKit❗❌:       // look around the beginning and end atoms and check for bonds with
    // RDKit❗❌:       // direction set
    // RDKit❗❌:       for (auto atom : {bond->getBeginAtom(), bond->getEndAtom()}) {
    // RDKit❗❌:         for (const auto &nbri :
    // RDKit❗❌:              boost::make_iterator_range(mol.getAtomBonds(atom))) {
    // RDKit❗❌:           const auto &obnd = mol[nbri];
    // RDKit❗❌:           if (obnd->getBondDir() == Bond::ENDDOWNRIGHT ||
    // RDKit❗❌:               obnd->getBondDir() == Bond::ENDUPRIGHT) {
    // RDKit❗❌:             bondsToClearDirs.push_back(obnd->getIdx());
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     auto tautBond = tautBonds[bondIdx];
    // RDKit❗❌:     if (tautBond->getBondType() != Bond::DOUBLE || d_removeBondStereo ||
    // RDKit❗❌:         !hasValidSpecifiedDoubleBondStereo(*bond)) {
    // RDKit❗❌:       // When bond stereo is being removed for bonds involved in tautomerism,
    // RDKit❗❌:       // use STEREOANY (for double bonds not in rings or connecting two ring atoms)
    // RDKit❗❌:       // instead of STEREONONE.
    // RDKit❗❌:       // This prevents downstream tools (notably InChI) from inferring a specific E/Z
    // RDKit❗❌:       // assignment from 2D coordinates after bond orders have been changed.
    // RDKit❗❌:       const RingInfo *ringInfo =
    // RDKit❗❌:           tautBond->getBondType() == Bond::DOUBLE ? getFastRingInfo(taut)
    // RDKit❗❌:                                                   : nullptr;
    // RDKit❗❌:       const auto targetStereo = getClearedTautomerBondStereo(ringInfo, *tautBond);
    // RDKit❗❌:       modified |= (tautBond->getStereo() != targetStereo);
    // RDKit❗❌:       tautBond->setStereo(targetStereo);
    // RDKit❗❌:       tautBond->getStereoAtoms().clear();
    // RDKit❗❌:       for (auto bi : bondsToClearDirs) {
    // RDKit❗❌:         tautBonds[bi]->setBondDir(Bond::NONE);
    // RDKit❗❌:       }
    // RDKit❗❌:     } else {
    // RDKit❗❌:       const INT_VECT &sa = bond->getStereoAtoms();
    // RDKit❗❌:       modified |= (tautBond->getStereo() != bond->getStereo() ||
    // RDKit❗❌:                    sa.size() != tautBond->getStereoAtoms().size());
    // RDKit❗❌:       if (sa.size() == 2) {
    // RDKit❗❌:         tautBond->setStereoAtoms(sa.front(), sa.back());
    // RDKit❗❌:       }
    // RDKit❗❌:       tautBond->setStereo(bond->getStereo());
    // RDKit❗❌:       for (auto bi : bondsToClearDirs) {
    // RDKit❗❌:         tautBonds[bi]->setBondDir(molBonds[bi]->getBondDir());
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   if (d_reassignStereo) {
    // RDKit❗❌:     static const bool cleanIt = true;
    // RDKit❗❌:     static const bool force = true;
    // RDKit❗❌:     MolOps::assignStereochemistry(taut, cleanIt, force);
    // RDKit❗❌:
    // RDKit❗❌:     // assignStereochemistry() can overwrite the explicit "undefined" bond
    // RDKit❗❌:     // stereo (STEREOANY) that we set above in order to prevent downstream
    // RDKit❗❌:     // coordinate-based E/Z inference. If bond stereo removal is enabled,
    // RDKit❗❌:     // re-apply our contract to the bonds involved in tautomerism.
    // RDKit❗❌:     if (d_removeBondStereo) {
    // RDKit❗❌:       const auto ringInfo = getFastRingInfo(taut);
    // RDKit❗❌:       for (auto bond : taut.bonds()) {
    // RDKit❗❌:         const auto bondIdx = bond->getIdx();
    // RDKit❗❌:         if (!res.d_modifiedBonds.test(bondIdx)) {
    // RDKit❗❌:           continue;
    // RDKit❗❌:         }
    // RDKit❗❌:         if (bond->getBondType() != Bond::DOUBLE) {
    // RDKit❗❌:           bond->setStereo(Bond::STEREONONE);
    // RDKit❗❌:           bond->getStereoAtoms().clear();
    // RDKit❗❌:           continue;
    // RDKit❗❌:         }
    // RDKit❗❌:         bond->setStereo(getClearedTautomerBondStereo(ringInfo, *bond));
    // RDKit❗❌:         bond->getStereoAtoms().clear();
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   } else {
    // RDKit❗❌:     taut.setProp(common_properties::_StereochemDone, 1);
    // RDKit❗❌:   }
    // RDKit❗❌:   return modified;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION RDKit::MolStandardize::TautomerEnumerator::setTautomerStereoAndIsoHs
    // Sparse BTreeSet iteration preserves ascending source IDs. Native bond
    // arrays already give O(1) indexed lookup: no redundant O(E) pointer caches.
    // Existing assignment clone/error/state adaptations remain baseline costs;
    // the complete source owner is not a blanket behavior/performance upgrade.
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

    for atom_id in modified_atoms.iter().copied() {
        let source_atom = &source.topology.atoms[atom_id.index()];
        let clear_isotopic_hydrogens = {
            let tautomer_atom = &tautomer.topology.atoms[atom_id.index()];
            !tautomer_atom.tracked_isotopic_hydrogens().is_empty()
                && (options.remove_isotopic_hydrogens()
                    || total_hydrogens(tautomer.view(coordinates), atom_id)? == 0)
        };
        let tautomer_atom = &mut tautomer.topology.atoms[atom_id.index()];

        if tautomer_atom.hybridization() == Hybridization::Sp2 || options.remove_sp3_stereo() {
            modified |= tautomer_atom.chiral_tag() != ChiralTag::Unspecified;
            tautomer_atom.set_chiral_tag(ChiralTag::Unspecified);
            if tautomer_atom.prop("_CIPCode").is_some() {
                tautomer_atom.clear_prop("_CIPCode")?;
            }
        } else {
            modified |= tautomer_atom.chiral_tag() != source_atom.chiral_tag();
            tautomer_atom.set_chiral_tag(source_atom.chiral_tag());
            if let Some(cip_code) = source_atom.prop("_CIPCode") {
                tautomer_atom.set_prop("_CIPCode", property_value_to_string(cip_code)?)?;
            }
        }
        if clear_isotopic_hydrogens {
            tautomer_atom.set_tracked_isotopic_hydrogens(Vec::new());
        }
    }

    for bond_id in modified_bonds.iter().copied() {
        let source_bond = &source.topology.bonds[bond_id.index()];
        let mut bonds_to_clear_directions = Vec::new();
        if source_bond.order() == BondOrder::Double && is_stereo_beyond_any(source_bond.stereo()) {
            for atom_id in [source_bond.begin(), source_bond.end()] {
                for neighbor in source.topology.adjacency.neighbors_of(atom_id.index()) {
                    let adjacent_bond = &source.topology.bonds[neighbor.bond.index()];
                    if matches!(
                        adjacent_bond.direction(),
                        BondDirection::EndDownRight | BondDirection::EndUpRight
                    ) {
                        bonds_to_clear_directions.push(adjacent_bond.id());
                    }
                }
            }
        }

        let tautomer_bond = &tautomer.topology.bonds[bond_id.index()];
        let candidate_order = tautomer_bond.order();
        let candidate_stereo = tautomer_bond.stereo();
        let candidate_stereo_reference_count = tautomer_bond.stereo_atom_references().len();
        let remove_stereo = tautomer_bond.order() != BondOrder::Double
            || options.remove_bond_stereo()
            || !has_valid_specified_double_bond_stereo(source_bond);
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
            source_bond.stereo_atoms()
        };
        modified |= candidate_stereo != target_stereo
            || (!remove_stereo
                && candidate_stereo_reference_count != source_bond.stereo_atom_references().len());

        let tautomer_bond = &mut tautomer.topology.bonds[bond_id.index()];
        tautomer_bond.set_stereo_atoms(target_stereo_atoms);
        tautomer_bond.set_stereo(target_stereo)?;
        for adjacent_bond_id in bonds_to_clear_directions {
            let direction = if remove_stereo {
                BondDirection::None
            } else {
                source.topology.bonds[adjacent_bond_id.index()].direction()
            };
            tautomer.topology.bonds[adjacent_bond_id.index()].set_direction(direction);
        }
    }

    if options.reassign_stereo() {
        assign_stereo(tautomer)?;

        if options.remove_bond_stereo() {
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
                let bond = &mut tautomer.topology.bonds[bond_id.index()];
                bond.set_stereo_atoms(None);
                bond.set_stereo(target_stereo)?;
            }
        }
    } else {
        tautomer.properties.set_prop("_StereochemDone", 1_i32)?;
    }
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
    // BEGIN RDKIT CPP FUNCTION RDKit::ROMol::initFromOther quickCopy source dependency
    // RDKit❗❌: void ROMol::initFromOther(const ROMol &other, bool quickCopy, int confId) {
    // RDKit❗❌:   if (this == &other) {
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:   numBonds = 0;
    // RDKit❗❌:   // std::cerr<<"    init from other: "<<this<<" "<<&other<<std::endl;
    // RDKit❗❌:   // copy over the atoms
    // RDKit❗❌:   // Avoid repeated reallocations when copying: for MolGraph's vecS vertex
    // RDKit❗❌:   // container, reserving upfront can reduce allocation churn.
    // RDKit❗❌:   d_graph.m_vertices.reserve(other.getNumAtoms());
    // RDKit❗❌:   for (const auto oatom : other.atoms()) {
    // RDKit❗❌:     constexpr bool updateLabel = false;
    // RDKit❗❌:     constexpr bool takeOwnership = true;
    // RDKit❗❌:     addAtom(oatom->copy(), updateLabel, takeOwnership);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // and the bonds:
    // RDKit❗❌:   for (const auto obond : other.bonds()) {
    // RDKit❗❌:     addBond(obond->copy(), true);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // ring information
    // RDKit❗❌:   delete dp_ringInfo;
    // RDKit❗❌:   if (other.dp_ringInfo) {
    // RDKit❗❌:     dp_ringInfo = new RingInfo(*(other.dp_ringInfo));
    // RDKit❗❌:   } else {
    // RDKit❗❌:     dp_ringInfo = new RingInfo();
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // enhanced stereochemical information
    // RDKit❗❌:   d_stereo_groups.clear();
    // RDKit❗❌:   d_stereo_groups.reserve(other.d_stereo_groups.size());
    // RDKit❗❌:   for (auto &otherGroup : other.d_stereo_groups) {
    // RDKit❗❌:     std::vector<Atom *> atoms;
    // RDKit❗❌:     for (auto &otherAtom : otherGroup.getAtoms()) {
    // RDKit❗❌:       atoms.push_back(getAtomWithIdx(otherAtom->getIdx()));
    // RDKit❗❌:     }
    // RDKit❗❌:     std::vector<Bond *> bonds;
    // RDKit❗❌:     for (auto &otherBond : otherGroup.getBonds()) {
    // RDKit❗❌:       bonds.push_back(getBondWithIdx(otherBond->getIdx()));
    // RDKit❗❌:     }
    // RDKit❗❌:     d_stereo_groups.emplace_back(otherGroup.getGroupType(), std::move(atoms),
    // RDKit❗❌:                                  std::move(bonds), otherGroup.getReadId());
    // RDKit❗❌:     d_stereo_groups.back().setWriteId(otherGroup.getWriteId());
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (other.dp_delAtoms) {
    // RDKit❗❌:     dp_delAtoms.reset(new boost::dynamic_bitset<>(*other.dp_delAtoms));
    // RDKit❗❌:   } else {
    // RDKit❗❌:     dp_delAtoms.reset(nullptr);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (other.dp_delBonds) {
    // RDKit❗❌:     dp_delBonds.reset(new boost::dynamic_bitset<>(*other.dp_delBonds));
    // RDKit❗❌:   } else {
    // RDKit❗❌:     dp_delBonds.reset(nullptr);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (!quickCopy) {
    // RDKit❗❌:     // copy conformations
    // RDKit❗❌:     for (const auto &conf : other.d_confs) {
    // RDKit❗❌:       if (confId < 0 || rdcast<int>(conf->getId()) == confId) {
    // RDKit❗❌:         this->addConformer(new Conformer(*conf));
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     // Copy sgroups
    // RDKit❗❌:     for (const auto &sg : getSubstanceGroups(other)) {
    // RDKit❗❌:       addSubstanceGroup(*this, sg);
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     d_props = other.d_props;
    // RDKit❗❌:
    // RDKit❗❌:     // Bookmarks should be copied as well:
    // RDKit❗❌:     for (auto abmI : other.d_atomBookmarks) {
    // RDKit❗❌:       for (const auto *aptr : abmI.second) {
    // RDKit❗❌:         setAtomBookmark(getAtomWithIdx(aptr->getIdx()), abmI.first);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     for (auto bbmI : other.d_bondBookmarks) {
    // RDKit❗❌:       for (const auto *bptr : bbmI.second) {
    // RDKit❗❌:         setBondBookmark(getBondWithIdx(bptr->getIdx()), bbmI.first);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   } else {
    // RDKit❗❌:     d_props.reset();
    // RDKit❗❌:     STR_VECT computed;
    // RDKit❗❌:     d_props.setVal(RDKit::detail::computedPropName, computed);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // std::cerr<<"---------    done init from other: "<<this<<"
    // RDKit❗❌:   // "<<&other<<std::endl;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION RDKit::ROMol::initFromOther quickCopy source dependency
    let mut topology = candidate.topology.clone();
    topology.substance_groups.clear();
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

    // BEGIN RECOVERY SEARCH04 SOURCE TautomerProductTryCatch
    // RDKit❗❌:           try {
    // RDKit❗❌:             // We only change bond orders/H counts/charges; the molecular graph
    // RDKit❗❌:             // (and therefore ring topology) is unchanged.
    // RDKit❗❌:             // `sanitizeMol()` always calls `clearComputedProps()` which resets
    // RDKit❗❌:             // ring info and forces ring-finding for each generated tautomer.
    // RDKit❗❌:             // Avoid that by clearing computed props without touching rings,
    // RDKit❗❌:             // then running the specific sanitize steps we need.
    // RDKit❗❌:             product->clearComputedProps(false);
    // RDKit❗❌:             product->updatePropertyCache(false);
    // RDKit❗❌:             MolOps::Kekulize(*product);
    // RDKit❗❌:             MolOps::setAromaticity(*product);
    // RDKit❗❌:             MolOps::setConjugation(*product);
    // RDKit❗❌:             MolOps::setHybridization(*product);
    // RDKit❗❌:             MolOps::adjustHs(*product);
    // RDKit❗❌:           } catch (const KekulizeException &) {
    // RDKit❗❌:             continue;
    // RDKit❗❌:           }
    // END RECOVERY SEARCH04 SOURCE TautomerProductTryCatch
    let mut properties = MoleculeProperties::default();
    properties.set_prop(
        "__computedProps",
        cosmolkit_model::PropertyValue::StringVector(Vec::new()),
    )?;
    let mut valence = candidate.valence.clone();
    let mut rings = candidate.rings.clone();
    // quickCopy retained the actual source ring and atom cache rows. The
    // source product sequence keeps rings and performs exactly one canonical
    // attempt; sanitizeMol would reset them and add a different retry policy.
    match cosmolkit_core::source_sanitize_tautomer_product(
        &mut topology,
        &mut properties,
        &mut valence,
        &mut rings,
    ) {
        Ok(()) => {}
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
    }
    let mut product = TautomerRecord {
        topology,
        coordinates: Some(CoordinateBlock::default()),
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

    // BEGIN RDKIT CPP FUNCTION RDKit::TautomerEnumerator::enumerate lazy product publication
    // RDKit❗❌:           res.d_tautomers[tsmiles] = Tautomer(
    // RDKit❗❌:               std::move(product),
    // RDKit❗❌:               numModifiedAtoms, numModifiedBonds);
    // END RDKIT CPP FUNCTION RDKit::TautomerEnumerator::enumerate lazy product publication
    Ok(TautomerExpansionAttempt::Product(TautomerExpandedProduct {
        tautomer: Arc::new(product),
        kekulized: None,
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
        coordinates: Some(view.coordinates.clone()),
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

pub(crate) fn get_cached_kekulized(
    candidate: &mut crate::ordered::TautomerCandidate<Arc<TautomerRecord>>,
) -> Result<Arc<TautomerRecord>, TautomerRunError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::MolStandardize::Tautomer::getKekulized
    // RDKit❗❌:   const ROMOL_SPTR &getKekulized() const {
    // RDKit❗❌:     if (!kekulized && tautomer) {
    // RDKit❗❌:       kekulized.reset(new RWMol(*tautomer));
    // RDKit❗❌:       MolOps::Kekulize(static_cast<RWMol &>(*kekulized), false, true);
    // RDKit❗❌:     }
    // RDKit❗❌:     return kekulized;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // END RDKIT CPP FUNCTION RDKit::MolStandardize::Tautomer::getKekulized
    if candidate.kekulized.is_none() {
        if let Some(tautomer) = candidate.tautomer.as_ref() {
            candidate.kekulized = Some(Arc::new(tautomer.as_ref().clone()));
            let copied = Arc::make_mut(candidate.kekulized.as_mut().expect("copy published above"));
            // Shared CHEM17 one-attempt algorithm mutates borrowed state and
            // retains completed partial changes on Err. No catch or fallback.
            #[cfg(test)]
            search04_failed_trace::record("published-copy", copied);
            let outcome = cosmolkit_core::source_kekulize_attempt(
                &mut copied.topology,
                &mut copied.valence,
                &mut copied.rings,
                &KekulizeParams {
                    mark_atoms_bonds: false,
                    canonical: true,
                    max_backtracks: 100,
                },
            );
            #[cfg(test)]
            if outcome.is_err() {
                search04_failed_trace::record("failed-copy-still-installed", copied);
            }
            outcome?;
        }
    }
    candidate.kekulized.clone().ok_or_else(|| {
        TautomerRunError::Control("tautomer has neither source nor cached Kekule branch".to_owned())
    })
}

#[cfg(test)]
mod search04_lazy_tests {
    use super::*;
    use crate::ordered::TautomerCandidate;
    fn candidate(record: TautomerRecord) -> TautomerCandidate<Arc<TautomerRecord>> {
        TautomerCandidate {
            tautomer: Some(Arc::new(record)),
            kekulized: None,
            num_modified_atoms: 0,
            num_modified_bonds: 0,
            done: false,
        }
    }
    #[test]
    fn search04_lazy_success_is_materialized_once_and_reused() {
        let source = stereo_tests::fixture_from_smiles("c1ccccc1").unwrap();
        let before = source.clone();
        let mut entry = candidate(source);
        assert!(entry.kekulized.is_none());
        let first = get_cached_kekulized(&mut entry).unwrap();
        let second = get_cached_kekulized(&mut entry).unwrap();
        assert!(Arc::ptr_eq(&first, &second));
        assert_eq!(entry.tautomer.as_deref(), Some(&before));
        assert!(first.topology.bonds.iter().all(|b| b.is_aromatic()));
        assert_eq!(
            first
                .topology
                .bonds
                .iter()
                .filter(|b| b.order() == BondOrder::Double)
                .count(),
            3
        );
    }
    #[test]
    fn search04_lazy_failed_copy_is_published_then_reused_without_retry() {
        let parsed = cosmolkit_smiles::parse_smiles(
            "c1cccc1",
            &cosmolkit_smiles::SmilesParseParams {
                sanitize: false,
                ..Default::default()
            },
        )
        .unwrap();
        let source = prepared(TautomerRecordView {
            topology: &parsed.topology,
            coordinates: &parsed.coordinates,
            properties: &parsed.properties,
            valence: None,
            rings: None,
        })
        .unwrap();
        let before = source.clone();
        let mut entry = candidate(source);
        let error = get_cached_kekulized(&mut entry).unwrap_err();
        assert!(matches!(
            error,
            TautomerRunError::Kekulize(KekulizeError::NotKekulizable { .. })
        ));
        let failed = entry.kekulized.as_ref().unwrap().clone();
        let failed_state = failed.as_ref().clone();
        let second = get_cached_kekulized(&mut entry).unwrap();
        assert!(Arc::ptr_eq(&failed, &second));
        assert_eq!(*second, failed_state);
        assert_eq!(entry.tautomer.as_deref(), Some(&before));
    }
    #[test]
    fn search04_empty_transform_catalog_never_requests_kekule_branch() {
        let parsed = cosmolkit_smiles::parse_smiles(
            "c1cccc1",
            &cosmolkit_smiles::SmilesParseParams {
                sanitize: false,
                ..Default::default()
            },
        )
        .unwrap();
        let source = TautomerRecordView {
            topology: &parsed.topology,
            coordinates: &parsed.coordinates,
            properties: &parsed.properties,
            valence: None,
            rings: None,
        };
        let catalog = crate::TautomerCatalog::from_data(&[]).unwrap();
        // No transform means no getKekulized call, even with a source that would
        // fail a later first attempt; normal stereo-pruning still executes.
        let mut source_rings = RingInfo::new(
            RingFindType::OtherOrUnknown,
            source.topology.atoms.len(),
            source.topology.bonds.len(),
        );
        let result = crate::enumerate_with_catalog(
            crate::TautomerScoreView::new(source, &mut source_rings),
            &catalog,
            TautomerParams::default().with_reassign_stereo(false),
            None,
        )
        .unwrap();
        assert_eq!(result.entries.len(), 1);
    }
    #[test]
    fn search04_product_is_lazy_and_quickcopy_drops_only_molecule_metadata() {
        let mut source = stereo_tests::fixture_from_smiles("CC=O").unwrap();
        source.properties.set_prop("retained-input", 7_i32).unwrap();
        source.topology.atoms[0]
            .set_prop("atom-local", 9_i32)
            .unwrap();
        source.topology.substance_groups.push(
            cosmolkit_model::SubstanceGroup::new(
                cosmolkit_model::SubstanceGroupId::new(0),
                cosmolkit_model::SubstanceGroupKind::Data,
            )
            .with_atoms(vec![AtomId::new(0)]),
        );
        let input = source.clone();
        let kmol = kekulized(&source).unwrap();
        let transform = crate::TautomerCatalog::current()
            .unwrap()
            .transforms()
            .iter()
            .find(|t| t.name().as_bytes() == b"1,3 (thio)keto/enol f")
            .unwrap()
            .clone();
        let coordinates = CoordinateBlock::default();
        let matched = transform_matches(&kmol, &coordinates, &transform)
            .unwrap()
            .remove(0);
        let attempt = apply_tautomer_transform_match(
            &source,
            &kmol,
            &coordinates,
            &transform,
            &matched,
            &BTreeSet::new(),
            &BTreeSet::new(),
            &|_| false,
            TautomerParams::default().with_reassign_stereo(false),
        )
        .unwrap();
        let TautomerExpansionAttempt::Product(product) = attempt else {
            panic!("new product expected")
        };
        assert!(product.kekulized.is_none());
        assert!(product.tautomer.topology.substance_groups.is_empty());
        assert_eq!(source.topology.substance_groups.len(), 1);
        assert!(product.tautomer.properties.prop("retained-input").is_none());
        assert_eq!(
            product.tautomer.topology.atoms[0].prop("atom-local"),
            input.topology.atoms[0].prop("atom-local")
        );
        assert_eq!(
            product.tautomer.coordinates.as_ref(),
            Some(&CoordinateBlock::default())
        );
        assert_eq!(source, input);
    }
}

#[cfg(test)]
mod search04_failed_trace {
    use super::TautomerRecord;
    use std::cell::RefCell;
    std::thread_local! {
        static EVENTS: RefCell<Option<Vec<(&'static str, TautomerRecord)>>> = const { RefCell::new(None) };
    }
    pub(super) fn record(event: &'static str, record: &TautomerRecord) {
        EVENTS.with(|events| {
            if let Some(rows) = events.borrow_mut().as_mut() {
                rows.push((event, record.clone()));
            }
        });
    }
    pub(super) struct Capture;
    impl Capture {
        pub(super) fn new() -> Self {
            EVENTS.with(|v| {
                assert!(v.borrow().is_none());
                *v.borrow_mut() = Some(Vec::new());
            });
            Self
        }
        pub(super) fn rows(&self) -> Vec<(&'static str, TautomerRecord)> {
            EVENTS.with(|v| v.borrow().as_ref().unwrap().clone())
        }
    }
    impl Drop for Capture {
        fn drop(&mut self) {
            EVENTS.with(|v| {
                v.borrow_mut().take();
            });
        }
    }
}

#[cfg(test)]
mod search04_registered_detached_tests {
    use super::*;
    #[test]
    fn search04_public_detached_enumeration_keeps_failed_cache_before_propagation() {
        let parsed = cosmolkit_smiles::parse_smiles(
            "c1cccc1",
            &cosmolkit_smiles::SmilesParseParams {
                sanitize: false,
                ..Default::default()
            },
        )
        .unwrap();
        let before = parsed.clone();
        let capture = search04_failed_trace::Capture::new();
        let view = TautomerRecordView {
            topology: &parsed.topology,
            coordinates: &parsed.coordinates,
            properties: &parsed.properties,
            valence: None,
            rings: None,
        };
        let mut source_rings = RingInfo::new(
            RingFindType::OtherOrUnknown,
            view.topology.atoms.len(),
            view.topology.bonds.len(),
        );
        let error = crate::enumerate_with_catalog(
            crate::TautomerScoreView::new(view, &mut source_rings),
            &crate::TautomerCatalog::current().unwrap(),
            TautomerParams::default(),
            None,
        )
        .unwrap_err();
        assert!(matches!(
            error,
            TautomerRunError::Kekulize(KekulizeError::NotKekulizable { .. })
        ));
        let rows = capture.rows();
        assert_eq!(
            rows.iter().map(|r| r.0).collect::<Vec<_>>(),
            ["published-copy", "failed-copy-still-installed"]
        );
        assert_eq!(rows[0].1.topology, parsed.topology);
        assert!(rows[1].1.topology.atoms.iter().all(|a| a.is_aromatic()));
        assert!(
            rows[1]
                .1
                .topology
                .bonds
                .iter()
                .all(|b| b.order() == BondOrder::Single && b.is_aromatic())
        );
        assert_eq!(rows[0].1.rings, rows[1].1.rings);
        assert_eq!(parsed, before);
    }
}

fn has_valid_specified_double_bond_stereo(bond: &cosmolkit_model::Bond) -> bool {
    // BEGIN RDKIT CPP FUNCTION RDKit::MolStandardize::hasValidSpecifiedDoubleBondStereo
    // RDKit✔️✔️: bool hasValidSpecifiedDoubleBondStereo(const Bond &bond) {
    // RDKit✔️✔️:   const auto stereo = bond.getStereo();
    // RDKit✔️✔️:   return bond.getBondType() == Bond::DOUBLE && stereo >= Bond::STEREOZ &&
    // RDKit✔️✔️:          stereo <= Bond::STEREOTRANS && bond.getStereoAtoms().size() == 2;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDKit::MolStandardize::hasValidSpecifiedDoubleBondStereo
    bond.order() == BondOrder::Double
        && matches!(
            bond.stereo(),
            BondStereo::Z | BondStereo::E | BondStereo::Cis | BondStereo::Trans
        )
        && bond.stereo_atom_references().len() == 2
}

/// Exclusive authority for the scoreRings source cache write. Graph and all
/// other molecular blocks remain readonly. It is neither Clone nor Copy.
pub struct TautomerScoreView<'a> {
    pub topology: &'a TopologyBlock,
    pub coordinates: &'a CoordinateBlock,
    pub properties: &'a MoleculeProperties,
    pub valence: Option<&'a ValenceAssignment>,
    pub(crate) rings: &'a mut RingInfo,
}
impl<'a> TautomerScoreView<'a> {
    #[doc(hidden)]
    pub fn new(view: TautomerRecordView<'a>, rings: &'a mut RingInfo) -> Self {
        Self {
            topology: view.topology,
            coordinates: view.coordinates,
            properties: view.properties,
            valence: view.valence,
            rings,
        }
    }
    pub fn as_record_view(&self) -> TautomerRecordView<'_> {
        TautomerRecordView {
            topology: self.topology,
            coordinates: self.coordinates,
            properties: self.properties,
            valence: self.valence,
            rings: Some(self.rings),
        }
    }
    pub fn reborrow(&mut self) -> TautomerScoreView<'_> {
        TautomerScoreView {
            topology: self.topology,
            coordinates: self.coordinates,
            properties: self.properties,
            valence: self.valence,
            rings: self.rings,
        }
    }
    #[doc(hidden)]
    pub fn ring_cache(&self) -> &RingInfo {
        self.rings
    }
    /// Transport only a completed cache from an invocation-owned score snapshot.
    /// Project ownership seam: the graph stays readonly and is checked exactly.
    #[doc(hidden)]
    pub fn retain_scored_cache_from(
        &mut self,
        snapshot: TautomerRecordView<'_>,
    ) -> Result<(), TautomerRunError> {
        if self.topology.atoms != snapshot.topology.atoms
            || self.topology.bonds != snapshot.topology.bonds
            || self.topology.adjacency != snapshot.topology.adjacency
        {
            return Err(TautomerRunError::ScoreCacheTargetMismatch);
        }
        if !self.rings.is_symm_sssr() {
            if let Some(rings) = snapshot.rings.filter(|rings| rings.is_symm_sssr()) {
                *self.rings = rings.clone();
            }
        }
        Ok(())
    }
}
impl TautomerRecord {
    pub fn score_view<'a>(&'a mut self, coordinates: &'a CoordinateBlock) -> TautomerScoreView<'a> {
        TautomerScoreView {
            topology: &self.topology,
            coordinates: self.coordinates.as_ref().unwrap_or(coordinates),
            properties: &self.properties,
            valence: Some(&self.valence),
            rings: &mut self.rings,
        }
    }
}
