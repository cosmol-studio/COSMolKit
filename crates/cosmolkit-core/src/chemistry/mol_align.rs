//! Read-only molecular alignment and RMSD calculations.

use thiserror::Error;

use crate::chemistry::numerics::alignment::{Transform3D, align_points};
use crate::read_parts::MoleculeReadParts;
use crate::{Conformer3D, Molecule, QueryGraph, SubstructMatchParams};

const DEFAULT_MAX_ITERATIONS: u32 = 50;
const DEFAULT_MAX_MATCHES: i32 = 1_000_000;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AlignmentAtomMap {
    pub probe_atom: usize,
    pub reference_atom: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AlignmentParameters {
    pub probe_conformer_id: i32,
    pub reference_conformer_id: i32,
    pub atom_map: Option<Vec<AlignmentAtomMap>>,
    pub weights: Option<Vec<f64>>,
    pub reflect: bool,
    pub max_iterations: u32,
}

impl Default for AlignmentParameters {
    fn default() -> Self {
        Self {
            probe_conformer_id: -1,
            reference_conformer_id: -1,
            atom_map: None,
            weights: None,
            reflect: false,
            max_iterations: DEFAULT_MAX_ITERATIONS,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct BestAlignmentParameters {
    pub probe_conformer_id: i32,
    pub reference_conformer_id: i32,
    pub atom_maps: Vec<Vec<AlignmentAtomMap>>,
    pub weights: Option<Vec<f64>>,
    pub reflect: bool,
    pub max_iterations: u32,
    pub max_matches: i32,
    pub symmetrize_conjugated_terminal_groups: bool,
    pub ignore_hydrogens: bool,
    pub num_threads: i32,
}

impl Default for BestAlignmentParameters {
    fn default() -> Self {
        Self {
            probe_conformer_id: -1,
            reference_conformer_id: -1,
            atom_maps: Vec::new(),
            weights: None,
            reflect: false,
            max_iterations: DEFAULT_MAX_ITERATIONS,
            max_matches: DEFAULT_MAX_MATCHES,
            symmetrize_conjugated_terminal_groups: true,
            ignore_hydrogens: true,
            num_threads: 1,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct AllConformerRmsdParameters {
    pub atom_maps: Vec<Vec<AlignmentAtomMap>>,
    pub weights: Option<Vec<f64>>,
    pub max_matches: i32,
    pub symmetrize_conjugated_terminal_groups: bool,
    pub ignore_hydrogens: bool,
    pub num_threads: i32,
}

impl Default for AllConformerRmsdParameters {
    fn default() -> Self {
        Self {
            atom_maps: Vec::new(),
            weights: None,
            max_matches: DEFAULT_MAX_MATCHES,
            symmetrize_conjugated_terminal_groups: true,
            ignore_hydrogens: true,
            num_threads: 1,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct CoordinateRmsdParameters {
    pub probe_conformer_id: i32,
    pub reference_conformer_id: i32,
    pub atom_maps: Vec<Vec<AlignmentAtomMap>>,
    pub weights: Option<Vec<f64>>,
    pub max_matches: i32,
    pub symmetrize_conjugated_terminal_groups: bool,
}

impl Default for CoordinateRmsdParameters {
    fn default() -> Self {
        Self {
            probe_conformer_id: -1,
            reference_conformer_id: -1,
            atom_maps: Vec::new(),
            weights: None,
            max_matches: DEFAULT_MAX_MATCHES,
            symmetrize_conjugated_terminal_groups: true,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AlignmentTransform {
    pub matrix: Transform3D,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AlignmentResult {
    pub rmsd: f64,
    pub transform: AlignmentTransform,
    pub atom_map: Vec<AlignmentAtomMap>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConformerRmsd {
    pub probe_conformer_id: usize,
    pub reference_conformer_id: usize,
    pub rmsd: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ConformerAlignmentReport {
    pub rmsds: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum AlignmentError {
    #[error("molecule has no 3D conformers")]
    NoConformers,
    #[error("conformer id {id} was not found")]
    ConformerNotFound { id: i32 },
    #[error("alignment atom map is empty")]
    EmptyAtomMap,
    #[error("probe atom index {index} is out of range for {atom_count} atoms")]
    ProbeAtomOutOfRange { index: usize, atom_count: usize },
    #[error("reference atom index {index} is out of range for {atom_count} atoms")]
    ReferenceAtomOutOfRange { index: usize, atom_count: usize },
    #[error("alignment has {map_len} mapped atoms but {weight_len} weights")]
    WeightCountMismatch { map_len: usize, weight_len: usize },
    #[error("alignment weight at index {index} must be positive")]
    NonPositiveWeight { index: usize },
    #[error("no substructure match found between the reference and probe molecules")]
    NoSubstructureMatch,
    #[error("terminal-group symmetrization failed: {message}")]
    TerminalGroupSymmetrization { message: &'static str },
    #[error("alignment numerical precondition failed: {message}")]
    NumericalPrecondition { message: &'static str },
    #[error("alignment worker terminated unexpectedly")]
    WorkerTerminated,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ConformerAlignmentParameters {
    pub atom_indices: Option<Vec<usize>>,
    pub conformer_ids: Option<Vec<usize>>,
    pub weights: Option<Vec<f64>>,
    pub reflect: bool,
    pub max_iterations: u32,
}

impl Default for ConformerAlignmentParameters {
    fn default() -> Self {
        Self {
            atom_indices: None,
            conformer_ids: None,
            weights: None,
            reflect: false,
            max_iterations: DEFAULT_MAX_ITERATIONS,
        }
    }
}

fn conformer_by_id(molecule: &Molecule, id: i32) -> Result<&Conformer3D, AlignmentError> {
    // RDKit✔️✔️: const Conformer &prbCnf = prbMol.getConformer(prbCid);
    // RDKit✔️✔️: const Conformer &refCnf = refMol.getConformer(refCid);
    let conformers = molecule.conformers_3d();
    let index = crate::chemistry::conformer_selection::resolve_3d_conformer_index(conformers, i64::from(id))
        .ok_or_else(|| {
            if conformers.is_empty() {
                AlignmentError::NoConformers
            } else {
                AlignmentError::ConformerNotFound { id }
            }
        })?;
    Ok(&conformers[index])
}

fn first_alignment_map(
    probe: &Molecule,
    reference: &Molecule,
    params: &AlignmentParameters,
) -> Result<Vec<AlignmentAtomMap>, AlignmentError> {
    if let Some(map) = &params.atom_map {
        return Ok(map.clone());
    }
    // BEGIN RDKIT CPP FUNCTION RDKit::MolAlign::getAlignmentTransform automatic map
    // RDKit✔️✔️:   const bool recursionPossible = true;
    // RDKit✔️✔️:   const bool useChirality = false;
    // RDKit✔️✔️:   const bool useQueryQueryMatches = true;
    // RDKit✔️✔️:   if (SubstructMatch(refMol, prbMol, match, recursionPossible,
    // RDKit✔️✔️:                      useChirality, useQueryQueryMatches)) {
    // RDKit✔️✔️:     atomMap = &match;
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION RDKit::MolAlign::getAlignmentTransform automatic map
    let matcher_params = SubstructMatchParams {
        max_matches: 1,
        use_chirality: false,
        use_query_query_matches: true,
        ..SubstructMatchParams::default()
    };
    let probe_query = QueryGraph::from_concrete_molecule(probe).map_err(|_| AlignmentError::NoSubstructureMatch)?;
    let matched = crate::try_get_substruct_matches_with_params(reference, &probe_query, &matcher_params)
        .map_err(|_| AlignmentError::NoSubstructureMatch)?
        .into_iter()
        .next()
        .ok_or(AlignmentError::NoSubstructureMatch)?;
    Ok(matched
        .atom_mapping
        .into_iter()
        .enumerate()
        .map(|(probe_atom, reference_atom)| AlignmentAtomMap {
            probe_atom,
            reference_atom,
        })
        .collect())
}

fn validate_map(
    map: &[AlignmentAtomMap],
    probe: &Conformer3D,
    reference: &Conformer3D,
    weights: Option<&[f64]>,
) -> Result<(), AlignmentError> {
    if map.is_empty() {
        return Err(AlignmentError::EmptyAtomMap);
    }
    for entry in map {
        if entry.probe_atom >= probe.coordinates().len() {
            return Err(AlignmentError::ProbeAtomOutOfRange {
                index: entry.probe_atom,
                atom_count: probe.coordinates().len(),
            });
        }
        if entry.reference_atom >= reference.coordinates().len() {
            return Err(AlignmentError::ReferenceAtomOutOfRange {
                index: entry.reference_atom,
                atom_count: reference.coordinates().len(),
            });
        }
    }
    if let Some(weights) = weights {
        if weights.len() != map.len() {
            return Err(AlignmentError::WeightCountMismatch {
                map_len: map.len(),
                weight_len: weights.len(),
            });
        }
        if let Some(index) = weights.iter().position(|weight| !(*weight > 0.0)) {
            return Err(AlignmentError::NonPositiveWeight { index });
        }
    }
    Ok(())
}

fn automatic_maps(
    probe: &Molecule,
    reference: &Molecule,
    max_matches: i32,
    symmetrize_conjugated_terminal_groups: bool,
    ignore_hydrogens: bool,
) -> Result<Vec<Vec<AlignmentAtomMap>>, AlignmentError> {
    // RDKit✔️✔️: bool uniquify = false;
    // RDKit✔️✔️: bool recursionPossible = true;
    // RDKit✔️✔️: bool useChirality = false;
    // RDKit✔️✔️: bool useQueryQueryMatches = false;
    // RDKit✔️✔️: SubstructMatch(refMol, prbMolForMatch, matches, uniquify,
    // RDKit✔️✔️:                recursionPossible, useChirality,
    // RDKit✔️✔️:                useQueryQueryMatches, maxMatches);
    if max_matches == 0 {
        return Err(AlignmentError::NoSubstructureMatch);
    }
    let symmetrized_probe;
    let probe_for_match = if symmetrize_conjugated_terminal_groups {
        symmetrized_probe = crate::chemistry::mol_align_support::symmetrize_terminal_atoms(probe)
            .map_err(|message| AlignmentError::TerminalGroupSymmetrization { message })?;
        &symmetrized_probe
    } else {
        probe
    };
    let matcher_params = SubstructMatchParams {
        uniquify: false,
        // RDKit's legacy overload accepts unsigned int; BestAlignmentParams
        // stores int and relies on the same modulo-2^32 conversion here.
        max_matches: max_matches as u32 as usize,
        use_chirality: false,
        ..Default::default()
    };
    let probe_query =
        QueryGraph::from_concrete_molecule(probe_for_match).map_err(|_| AlignmentError::NoSubstructureMatch)?;
    let mut maps: Vec<_> = crate::try_get_substruct_matches_with_params(reference, &probe_query, &matcher_params)
        .map_err(|_| AlignmentError::NoSubstructureMatch)?
        .into_iter()
        .map(|matched| {
            matched
                .atom_mapping
                .into_iter()
                .enumerate()
                .filter(|(probe_atom, _)| !ignore_hydrogens || probe.atoms()[*probe_atom].atomic_number() != 1)
                .map(|(probe_atom, reference_atom)| AlignmentAtomMap {
                    probe_atom,
                    reference_atom,
                })
                .collect::<Vec<_>>()
        })
        .collect();
    maps.retain(|map| !map.is_empty());
    if maps.is_empty() {
        Err(AlignmentError::NoSubstructureMatch)
    } else {
        Ok(maps)
    }
}

fn maps_for(
    probe: &Molecule,
    reference: &Molecule,
    atom_maps: &[Vec<AlignmentAtomMap>],
    max_matches: i32,
    symmetrize_conjugated_terminal_groups: bool,
    ignore_hydrogens: bool,
) -> Result<Vec<Vec<AlignmentAtomMap>>, AlignmentError> {
    if atom_maps.is_empty() {
        automatic_maps(
            probe,
            reference,
            max_matches,
            symmetrize_conjugated_terminal_groups,
            ignore_hydrogens,
        )
    } else {
        Ok(atom_maps.to_vec())
    }
}

fn aligned_result(
    probe: &Conformer3D,
    reference: &Conformer3D,
    map: &[AlignmentAtomMap],
    weights: Option<&[f64]>,
    reflect: bool,
    max_iterations: u32,
) -> Result<AlignmentResult, AlignmentError> {
    // RDKit✔️✔️: for (const auto &mi : atomMap) {
    // RDKit✔️✔️:   prbPoints.push_back(&prbCnf.getAtomPos(mi.first));
    // RDKit✔️✔️:   refPoints.push_back(&refCnf.getAtomPos(mi.second));
    // RDKit✔️✔️: }
    // RDKit✔️✔️: double ssr = RDNumeric::Alignments::AlignPoints(
    // RDKit✔️✔️:     refPoints, prbPoints, trans, weights, reflect, maxIterations);
    // RDKit✔️✔️: return ssr / static_cast<double>(prbPoints.size());
    validate_map(map, probe, reference, weights)?;
    let probe_points: Vec<_> = map.iter().map(|entry| probe.coordinates()[entry.probe_atom]).collect();
    let reference_points: Vec<_> = map
        .iter()
        .map(|entry| reference.coordinates()[entry.reference_atom])
        .collect();
    let (ssr, matrix) = align_points(
        &reference_points,
        &probe_points,
        weights,
        reflect,
        max_iterations as usize,
    )
    .map_err(|message| AlignmentError::NumericalPrecondition { message })?;
    Ok(AlignmentResult {
        rmsd: (ssr / map.len() as f64).sqrt(),
        transform: AlignmentTransform { matrix },
        atom_map: map.to_vec(),
    })
}

fn num_threads_to_use(target: i32) -> usize {
    // BEGIN RDKIT CPP FUNCTION RDKit::getNumThreadsToUse
    // RDKit✔️✔️: if (target >= 1) {
    // RDKit✔️✔️:   return static_cast<unsigned int>(target);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: unsigned int res = std::thread::hardware_concurrency();
    // RDKit✔️✔️: if (res > rdcast<unsigned int>(-target)) {
    // RDKit✔️✔️:   return res + target;
    // RDKit✔️✔️: } else {
    // RDKit✔️✔️:   return 1;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDKit::getNumThreadsToUse
    if target >= 1 {
        return target as usize;
    }
    let available = std::thread::available_parallelism().map_or(1, std::num::NonZeroUsize::get);
    let offset = target.unsigned_abs() as usize;
    if available > offset { available - offset } else { 1 }
}

fn best_aligned_result(
    probe: &Conformer3D,
    reference: &Conformer3D,
    maps: &[Vec<AlignmentAtomMap>],
    weights: Option<&[f64]>,
    reflect: bool,
    max_iterations: u32,
    num_threads: i32,
) -> Result<AlignmentResult, AlignmentError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::MolAlign::getBestRMSInternal
    // RDKit✔️✔️: if (numThreads == 1) {
    // RDKit✔️✔️:   for (const auto &matche : matches) {
    // RDKit✔️✔️:     double msd = alignConfsOnAtomMap(prbCnf, refCnf, matche,
    // RDKit✔️✔️:                                             tmpTrans, weights, reflect, maxIters);
    // RDKit✔️✔️:     if (msd < msdBest) {
    // RDKit✔️✔️:       msdBest = msd;
    // RDKit✔️✔️:       bestMatchPtr = &matche;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: for (auto midx = tidx; midx < matches.size(); midx += numThreads) {
    // RDKit✔️✔️:   rmsds[tidx].emplace_back(msd, midx, tmpTrans);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: for (const auto &rv : rmsds) {
    // RDKit✔️✔️:   for (const auto &res : rv) {
    // RDKit✔️✔️:     if (msd < msdBest) {
    // RDKit✔️✔️:       msdBest = msd;
    // RDKit✔️✔️:       bestMatchPtr = &matches[midx];
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDKit::MolAlign::getBestRMSInternal
    if maps.is_empty() {
        return Err(AlignmentError::NoSubstructureMatch);
    }
    for map in maps {
        validate_map(map, probe, reference, weights)?;
    }
    let thread_count = num_threads_to_use(num_threads);
    if thread_count == 1 {
        let mut best: Option<AlignmentResult> = None;
        for map in maps {
            let candidate = aligned_result(probe, reference, map, weights, reflect, max_iterations)?;
            if best.as_ref().is_none_or(|current| candidate.rmsd < current.rmsd) {
                best = Some(candidate);
            }
        }
        return best.ok_or(AlignmentError::NoSubstructureMatch);
    }

    let buckets = std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(thread_count);
        for thread_index in 0..thread_count {
            handles.push(scope.spawn(move || {
                let mut bucket = Vec::new();
                let mut map_index = thread_index;
                while map_index < maps.len() {
                    bucket.push((
                        map_index,
                        aligned_result(probe, reference, &maps[map_index], weights, reflect, max_iterations),
                    ));
                    map_index += thread_count;
                }
                bucket
            }));
        }
        handles.into_iter().map(|handle| handle.join()).collect::<Vec<_>>()
    });
    let mut best: Option<AlignmentResult> = None;
    for bucket in buckets {
        let bucket = bucket.map_err(|_| AlignmentError::WorkerTerminated)?;
        for (_, candidate) in bucket {
            let candidate = candidate?;
            if best.as_ref().is_none_or(|current| candidate.rmsd < current.rmsd) {
                best = Some(candidate);
            }
        }
    }
    best.ok_or(AlignmentError::NoSubstructureMatch)
}

fn coordinate_rmsd(
    probe: &Conformer3D,
    reference: &Conformer3D,
    map: &[AlignmentAtomMap],
    weights: Option<&[f64]>,
) -> Result<f64, AlignmentError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::MolAlign::calcMSDInternal
    // RDKit✔️✔️:   unsigned int npt = atomMap.size();
    // RDKit✔️✔️:   if (!weights) {
    // RDKit✔️✔️:     unitWeights.reset(new RDNumeric::DoubleVector(npt, 1.0));
    // RDKit✔️✔️:     weights = unitWeights.get();
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     PRECONDITION(npt == weights->size(), "Mismatch in number of weights");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (const auto &mi : atomMap) {
    // RDKit✔️✔️:     prbPoints.push_back(&prbCnf.getAtomPos(mi.first));
    // RDKit✔️✔️:     refPoints.push_back(&refCnf.getAtomPos(mi.second));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (unsigned int i = 0; i < npt; ++i) {
    // RDKit✔️✔️:     ssr += (*weights)[i] * (*ppt - *rpt).lengthSq();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return ssr / static_cast<double>(npt);
    // END RDKIT CPP FUNCTION RDKit::MolAlign::calcMSDInternal
    validate_map(map, probe, reference, weights)?;
    let mut sum = 0.0;
    for (index, entry) in map.iter().enumerate() {
        let p = probe.coordinates()[entry.probe_atom];
        let r = reference.coordinates()[entry.reference_atom];
        let weight = weights.map_or(1.0, |weights| weights[index]);
        sum += weight * ((p[0] - r[0]).powi(2) + (p[1] - r[1]).powi(2) + (p[2] - r[2]).powi(2));
    }
    Ok((sum / map.len() as f64).sqrt())
}

#[cfg(test)]
std::thread_local! {
    static ALIGN_CONFORMERS_CALL_COUNT: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

#[cfg(test)]
pub(crate) fn reset_align_conformers_call_count() {
    ALIGN_CONFORMERS_CALL_COUNT.set(0);
}

#[cfg(test)]
pub(crate) fn align_conformers_call_count() -> usize {
    ALIGN_CONFORMERS_CALL_COUNT.get()
}

pub(crate) fn align_conformers_in_coordinate_block(
    coordinates: &mut crate::molecule::CoordinateBlock,
    atom_count: usize,
    params: &ConformerAlignmentParameters,
) -> Result<Vec<f64>, AlignmentError> {
    #[cfg(test)]
    ALIGN_CONFORMERS_CALL_COUNT.set(ALIGN_CONFORMERS_CALL_COUNT.get() + 1);

    // BEGIN RDKIT CPP FUNCTION RDKit::MolAlign::alignMolConformers
    // RDKit✔️✔️: if (mol.getNumConformers() == 0) {
    // RDKit✔️✔️:   return;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if ((confIds != nullptr) && (confIds->size() > 0)) {
    // RDKit✔️✔️:   cid = confIds->front();
    // RDKit✔️✔️: }
    // RDKit✔️✔️: const Conformer &refCnf = mol.getConformer(cid);
    // RDKit✔️✔️: _fillAtomPositions(refPoints, refCnf, atomIds);
    // RDKit✔️✔️: for (cai = confIds->begin(); cai != confIds->end(); cai++) {
    // RDKit✔️✔️:   if (i == 1) { continue; }
    // RDKit✔️✔️:   ssd = RDNumeric::Alignments::AlignPoints(refPoints, prbPoints, trans,
    // RDKit✔️✔️:                                            weights, reflect, maxIters);
    // RDKit✔️✔️:   ssd /= (prbPoints.size());
    // RDKit✔️✔️:   RMSlist->push_back(sqrt(ssd));
    // RDKit✔️✔️:   MolTransforms::transformConformer(conf, trans);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDKit::MolAlign::alignMolConformers
    if coordinates.conformers_3d.is_empty() {
        return Ok(Vec::new());
    }
    let atoms: Vec<usize> = params.atom_indices.clone().unwrap_or_else(|| (0..atom_count).collect());
    if atoms.is_empty() {
        return Err(AlignmentError::EmptyAtomMap);
    }
    if let Some(&index) = atoms.iter().find(|&&index| index >= atom_count) {
        return Err(AlignmentError::ProbeAtomOutOfRange { index, atom_count });
    }
    if let Some(weights) = &params.weights {
        if weights.len() != atoms.len() {
            return Err(AlignmentError::WeightCountMismatch {
                map_len: atoms.len(),
                weight_len: weights.len(),
            });
        }
        if let Some(index) = weights.iter().position(|weight| !(*weight > 0.0)) {
            return Err(AlignmentError::NonPositiveWeight { index });
        }
    }
    let selected: Vec<usize> = match &params.conformer_ids {
        None => (0..coordinates.conformers_3d.len()).collect(),
        Some(ids) if ids.is_empty() => (0..coordinates.conformers_3d.len()).collect(),
        Some(ids) => ids
            .iter()
            .map(|&id| {
                coordinates
                    .conformers_3d
                    .iter()
                    .position(|conformer| conformer.id() == id)
                    .ok_or(AlignmentError::ConformerNotFound { id: id as i32 })
            })
            .collect::<Result<_, _>>()?,
    };
    if selected.is_empty() {
        return Ok(Vec::new());
    }
    let reference_points: Vec<_> = atoms
        .iter()
        .map(|&atom| coordinates.conformers_3d[selected[0]].coordinates()[atom])
        .collect();
    let mut rmsds = Vec::with_capacity(selected.len().saturating_sub(1));
    for &conformer_index in selected.iter().skip(1) {
        let probe_points: Vec<_> = atoms
            .iter()
            .map(|&atom| coordinates.conformers_3d[conformer_index].coordinates()[atom])
            .collect();
        let (ssr, transform) = align_points(
            &reference_points,
            &probe_points,
            params.weights.as_deref(),
            params.reflect,
            params.max_iterations as usize,
        )
        .map_err(|message| AlignmentError::NumericalPrecondition { message })?;
        rmsds.push((ssr / atoms.len() as f64).sqrt());
        for point in coordinates.conformers_3d[conformer_index].coordinates_mut() {
            *point = crate::chemistry::numerics::alignment::transform_point(&transform, *point);
        }
    }
    Ok(rmsds)
}

pub(crate) fn alignment_result_from_read_parts(
    read: MoleculeReadParts<'_>,
    reference: &Molecule,
    params: &AlignmentParameters,
) -> Result<AlignmentResult, AlignmentError> {
    let probe = Molecule::from_operation_blocks(
        read.topology().clone(),
        read.coordinates().clone(),
        read.properties().clone(),
        read.derived_cache().clone(),
        read.capabilities(),
    )
    .expect("an operation read view must preserve molecule invariants");
    probe.alignment_transform_to(reference, params)
}

pub(crate) fn apply_alignment_result_to_coordinate_block(
    coordinates: &mut crate::molecule::CoordinateBlock,
    probe_conformer_id: i32,
    result: &AlignmentResult,
) -> Result<(), AlignmentError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::MolAlign::alignMol
    // RDKit✔️✔️: double res = getAlignmentTransform(prbMol, refMol, trans, prbCid, refCid,
    // RDKit✔️✔️:                                    atomMap, weights, reflect, maxIters);
    // RDKit✔️✔️: Conformer &prbCnf = prbMol.getConformer(prbCid);
    // RDKit✔️✔️: MolTransforms::transformConformer(prbCnf, trans);
    // RDKit✔️✔️: return res;
    // END RDKIT CPP FUNCTION RDKit::MolAlign::alignMol
    let conformer_index = crate::chemistry::conformer_selection::resolve_3d_conformer_index(
        &coordinates.conformers_3d,
        i64::from(probe_conformer_id),
    )
    .ok_or_else(|| {
        if coordinates.conformers_3d.is_empty() {
            AlignmentError::NoConformers
        } else {
            AlignmentError::ConformerNotFound { id: probe_conformer_id }
        }
    })?;
    for point in coordinates.conformers_3d[conformer_index].coordinates_mut() {
        *point = crate::chemistry::numerics::alignment::transform_point(&result.transform.matrix, *point);
    }
    Ok(())
}

impl Molecule {
    pub fn alignment_transform_to(
        &self,
        reference: &Self,
        params: &AlignmentParameters,
    ) -> Result<AlignmentResult, AlignmentError> {
        let probe = conformer_by_id(self, params.probe_conformer_id)?;
        let reference_conformer = conformer_by_id(reference, params.reference_conformer_id)?;
        let map = first_alignment_map(self, reference, params)?;
        aligned_result(
            probe,
            reference_conformer,
            &map,
            params.weights.as_deref(),
            params.reflect,
            params.max_iterations,
        )
    }

    pub fn best_alignment_to(
        &self,
        reference: &Self,
        params: &BestAlignmentParameters,
    ) -> Result<AlignmentResult, AlignmentError> {
        let probe = conformer_by_id(self, params.probe_conformer_id)?;
        let reference_conformer = conformer_by_id(reference, params.reference_conformer_id)?;
        let maps = maps_for(
            self,
            reference,
            &params.atom_maps,
            params.max_matches,
            params.symmetrize_conjugated_terminal_groups,
            params.ignore_hydrogens,
        )?;
        best_aligned_result(
            probe,
            reference_conformer,
            &maps,
            params.weights.as_deref(),
            params.reflect,
            params.max_iterations,
            params.num_threads,
        )
    }

    pub fn best_rmsd_to(&self, reference: &Self, params: &BestAlignmentParameters) -> Result<f64, AlignmentError> {
        Ok(self.best_alignment_to(reference, params)?.rmsd)
    }

    pub fn coordinate_rmsd_to(
        &self,
        reference: &Self,
        params: &CoordinateRmsdParameters,
    ) -> Result<f64, AlignmentError> {
        let probe = conformer_by_id(self, params.probe_conformer_id)?;
        let reference_conformer = conformer_by_id(reference, params.reference_conformer_id)?;
        let maps = maps_for(
            self,
            reference,
            &params.atom_maps,
            params.max_matches,
            params.symmetrize_conjugated_terminal_groups,
            false,
        )?;
        let mut best = f64::INFINITY;
        for map in maps {
            let value = coordinate_rmsd(probe, reference_conformer, &map, params.weights.as_deref())?;
            if value < best {
                best = value;
            }
        }
        Ok(best)
    }

    pub fn all_conformer_best_rmsds(
        &self,
        params: &AllConformerRmsdParameters,
    ) -> Result<Vec<ConformerRmsd>, AlignmentError> {
        // BEGIN RDKIT CPP FUNCTION RDKit::MolAlign::getAllConformerBestRMS
        // RDKit✔️✔️: std::vector<int> cids;
        // RDKit✔️✔️: for (auto cit = mol.beginConformers(); cit != mol.endConformers(); ++cit) {
        // RDKit✔️✔️:   cids.push_back((*cit)->getId());
        // RDKit✔️✔️: }
        // RDKit✔️✔️: for (auto ci = 0u; ci < mol.getNumConformers(); ++ci) {
        // RDKit✔️✔️:   for (auto cj = 0u; cj < ci; ++cj) {
        // RDKit✔️✔️:     res.push_back(getBestRMSInternal(mol, mol, cids[ci], cids[cj],
        // RDKit✔️✔️:                                      matches, &trans, nullptr, params.weights,
        // RDKit✔️✔️:                                      reflect, maxIters, 1));
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // RDKit✔️✔️: for (auto i = tidx; i < pairs.size(); i += numThreads) {
        // RDKit✔️✔️:   rmsds[tidx].emplace_back(i, rms);
        // RDKit✔️✔️: }
        // RDKit✔️✔️: res[v.first] = v.second;
        // END RDKIT CPP FUNCTION RDKit::MolAlign::getAllConformerBestRMS
        let conformers = self.conformers_3d();
        if conformers.len() < 2 {
            return Ok(Vec::new());
        }
        let maps = maps_for(
            self,
            self,
            &params.atom_maps,
            params.max_matches,
            params.symmetrize_conjugated_terminal_groups,
            params.ignore_hydrogens,
        )?;
        let mut pairs = Vec::with_capacity(conformers.len() * (conformers.len() - 1) / 2);
        for probe_index in 0..conformers.len() {
            for reference_index in 0..probe_index {
                pairs.push((probe_index, reference_index));
            }
        }
        let evaluate = |pair_index: usize| -> Result<ConformerRmsd, AlignmentError> {
            let (probe_index, reference_index) = pairs[pair_index];
            let aligned = best_aligned_result(
                &conformers[probe_index],
                &conformers[reference_index],
                &maps,
                params.weights.as_deref(),
                false,
                DEFAULT_MAX_ITERATIONS,
                1,
            )?;
            Ok(ConformerRmsd {
                probe_conformer_id: conformers[probe_index].id(),
                reference_conformer_id: conformers[reference_index].id(),
                rmsd: aligned.rmsd,
            })
        };
        let thread_count = num_threads_to_use(params.num_threads);
        if thread_count == 1 {
            return (0..pairs.len()).map(evaluate).collect();
        }
        let pair_count = pairs.len();
        let buckets = std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(thread_count);
            for thread_index in 0..thread_count {
                let evaluate = &evaluate;
                handles.push(scope.spawn(move || {
                    let mut bucket = Vec::new();
                    let mut pair_index = thread_index;
                    while pair_index < pair_count {
                        bucket.push((pair_index, evaluate(pair_index)));
                        pair_index += thread_count;
                    }
                    bucket
                }));
            }
            handles.into_iter().map(|handle| handle.join()).collect::<Vec<_>>()
        });
        let mut ordered = vec![None; pair_count];
        for bucket in buckets {
            for (pair_index, result) in bucket.map_err(|_| AlignmentError::WorkerTerminated)? {
                ordered[pair_index] = Some(result);
            }
        }
        ordered
            .into_iter()
            .map(|entry| entry.ok_or(AlignmentError::WorkerTerminated)?)
            .collect()
    }
}
