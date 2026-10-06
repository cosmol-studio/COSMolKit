//! Source-backed AtomPair fingerprints over explicitly prepared detached blocks.
use crate::generator::{
    FingerprintArguments, FingerprintEnvironment, FingerprintFuncArguments,
    accumulate_sparse_counts_into, project_count_fingerprint, project_fingerprint,
    project_sparse_fingerprint, with_fingerprint_environment_inputs_and_output,
};
use crate::{
    AtomCodeError, AtomCodeInput, AtomCodeOptions, Fingerprint, FingerprintAdditionalOutput,
    FingerprintError, MorganError, SparseBitFingerprint, SparseCountFingerprint,
    SparseCountFingerprint32, atom_code, atom_pair_code, hash::hash_combine,
};
use cosmolkit_core::{
    DistanceMatrix3dParams, MatrixError, RingInfo, TopologicalDistanceMatrixParams,
    ValenceAssignment, distance_matrix_3d, topological_distance_matrix,
};
use cosmolkit_model::{
    CoordinateBlock, MoleculeProperties, TopologyBlock, TopologyValidationError,
};
use std::fmt;

#[derive(Debug)]
pub enum AtomPairError {
    Fingerprint(FingerprintError),
    AtomCode(AtomCodeError),
    Matrix(MatrixError),
    Topology(TopologyValidationError),
    Preparation(MorganError),
    AtomInvariantLength { length: usize, atom_index: usize },
}
impl fmt::Display for AtomPairError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Fingerprint(e) => e.fmt(f),
            Self::AtomCode(e) => e.fmt(f),
            Self::Matrix(e) => e.fmt(f),
            Self::Topology(e) => e.fmt(f),
            Self::Preparation(e) => e.fmt(f),
            Self::AtomInvariantLength { .. } => f.write_str("bad atom invariants size"),
        }
    }
}
impl std::error::Error for AtomPairError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Fingerprint(e) => Some(e),
            Self::AtomCode(e) => Some(e),
            Self::Matrix(e) => Some(e),
            Self::Topology(e) => Some(e),
            Self::Preparation(e) => Some(e),
            Self::AtomInvariantLength { .. } => None,
        }
    }
}
impl From<FingerprintError> for AtomPairError {
    fn from(e: FingerprintError) -> Self {
        Self::Fingerprint(e)
    }
}
impl From<MorganError> for AtomPairError {
    fn from(e: MorganError) -> Self {
        Self::Preparation(e)
    }
}
impl From<MatrixError> for AtomPairError {
    fn from(e: MatrixError) -> Self {
        Self::Matrix(e)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomPairParams {
    pub min_distance: u32,
    pub max_distance: u32,
    pub include_chirality: bool,
    pub use_2d: bool,
    pub count_simulation: bool,
    pub fp_size: u32,
    pub count_bounds: Vec<u32>,
    pub bits_per_feature: u32,
}
impl Default for AtomPairParams {
    fn default() -> Self {
        Self {
            min_distance: 1,
            max_distance: 30,
            include_chirality: false,
            use_2d: true,
            count_simulation: true,
            fp_size: 2048,
            count_bounds: vec![1, 2, 4, 8],
            bits_per_feature: 1,
        }
    }
}
/// Preserve the pinned project-native batch wrapper's preflight ordering.
/// Modern scalar RDKit entrypoints retain their independent zero-size semantics.
pub fn validate_atom_pair_params(params: &AtomPairParams) -> Result<(), FingerprintError> {
    AtomPairBatchArguments::new(params).map(|_| ())
}

fn batch_common(params: &AtomPairParams) -> Result<FingerprintArguments, FingerprintError> {
    // COSMolKit❗✔️: d892ec3507c5b568c5ed5d86ae44e466f7d03855
    // properties/fingerprint/atom_pair.rs:812-842. The detached parameters
    // already carry a u32 width; the Python projection performs the source
    // checked usize conversion before invoking this narrow owner seam.
    //     pub fn new(params: &AtomPairFingerprintParams) -> Result<Self, FingerprintError> {
    //         if params.n_bits == 0 {
    //             return Err(FingerprintError::EmptyFingerprint);
    //         }
    //         let fp_size =
    //             u32::try_from(params.n_bits).map_err(|_| FingerprintError::InvalidArguments {
    //                 reason: "AtomPair n_bits exceeds the source uint32 range",
    //             })?;
    //         let mut generator = atom_pair_generator_with_parameters(
    //             params.min_distance,
    //             params.max_distance,
    //             params.use_chirality,
    //             params.use_2d,
    //             None,
    //             params.count_simulation,
    //             fp_size,
    //             params.count_bounds.clone(),
    //             false,
    //         )?;
    //         if params.num_bits_per_feature == 0 {
    //             return Err(FingerprintError::InvalidArguments {
    //                 reason: "num_bits_per_feature must be greater than zero",
    //             });
    //         }
    //         generator
    //             .arguments
    //             .fingerprint_arguments
    //             .d_num_bits_per_feature = params.num_bits_per_feature;
    //         Ok(generator)
    //     }
    if params.fp_size == 0 {
        return Err(FingerprintError::EmptyFingerprint);
    }
    let mut common = params.common_with_bits_per_feature(1)?;
    if params.bits_per_feature == 0 {
        return Err(FingerprintError::InvalidArguments {
            reason: "num_bits_per_feature must be greater than zero",
        });
    }
    common.bits_per_feature = params.bits_per_feature;
    Ok(common)
}

/// Immutable construction state for the existing batch owner boundary.
/// Never exposed by the public chemistry facade or language projections.
#[doc(hidden)]
pub struct AtomPairBatchArguments<'a> {
    params: &'a AtomPairParams,
    common: FingerprintArguments,
}

impl<'a> AtomPairBatchArguments<'a> {
    pub fn new(params: &'a AtomPairParams) -> Result<Self, FingerprintError> {
        // COSMolKit❗✔️: pinned batch::atom_pair_generator_for_batch constructs
        // one generator before collect_optional_values_with_options.
        // The source constructor anchor and ordered guards live in batch_common.
        // Cost: one O(count_bounds.len()) copy, then shared immutable borrows.
        Ok(Self {
            params,
            common: batch_common(params)?,
        })
    }
    pub fn sparse_count(
        &self,
        input: &AtomPairPreparedInput<'_>,
        call: &AtomPairCall<'_>,
        output: Option<&mut FingerprintAdditionalOutput>,
    ) -> Result<SparseCountFingerprint, AtomPairError> {
        atom_pair_sparse_count_with_common(input, self.params, &self.common, call, output)
    }
    pub fn sparse_bits(
        &self,
        input: &AtomPairPreparedInput<'_>,
        call: &AtomPairCall<'_>,
        output: Option<&mut FingerprintAdditionalOutput>,
    ) -> Result<SparseBitFingerprint, AtomPairError> {
        atom_pair_sparse_bits_with_common(input, self.params, &self.common, call, output)
    }
    pub fn count(
        &self,
        input: &AtomPairPreparedInput<'_>,
        call: &AtomPairCall<'_>,
        output: Option<&mut FingerprintAdditionalOutput>,
    ) -> Result<SparseCountFingerprint32, AtomPairError> {
        atom_pair_count_with_common(input, self.params, &self.common, call, output)
    }
    pub fn bits(
        &self,
        input: &AtomPairPreparedInput<'_>,
        call: &AtomPairCall<'_>,
        output: Option<&mut FingerprintAdditionalOutput>,
    ) -> Result<Fingerprint, AtomPairError> {
        atom_pair_bits_with_common(input, self.params, &self.common, call, output)
    }
}

impl AtomPairParams {
    fn common(&self) -> Result<FingerprintArguments, FingerprintError> {
        self.common_with_bits_per_feature(self.bits_per_feature)
    }
    fn common_with_bits_per_feature(
        &self,
        bits_per_feature: u32,
    ) -> Result<FingerprintArguments, FingerprintError> {
        // RDKit❗✔️: AtomPairArguments::AtomPairArguments(
        // RDKit❗✔️:     const bool countSimulation, const bool includeChirality, const bool use2D,
        // RDKit❗✔️:     const unsigned int minDistance, const unsigned int maxDistance,
        // RDKit❗✔️:     const std::vector<std::uint32_t> countBounds, const std::uint32_t fpSize)
        // RDKit❗✔️:     : FingerprintArguments(countSimulation, countBounds, fpSize, 1,
        // RDKit❗✔️:                            includeChirality),
        // RDKit❗✔️:       df_use2D(use2D),
        // RDKit❗✔️:       d_minDistance(minDistance),
        // RDKit❗✔️:       d_maxDistance(maxDistance) {
        // RDKit❗✔️:   PRECONDITION(minDistance <= maxDistance, "bad distances provided");
        // RDKit❗✔️: }
        let common = FingerprintArguments::new(
            self.count_simulation,
            self.count_bounds.clone(),
            self.fp_size,
            bits_per_feature,
            self.include_chirality,
        )?;
        if self.min_distance > self.max_distance {
            return Err(FingerprintError::PreconditionViolation {
                what: "bad distances provided",
            });
        }
        Ok(common)
    }
    fn result_size(&self) -> u64 {
        // RDKit❗✔️: OutputType AtomPairEnvGenerator<OutputType>::getResultSize() const {
        // RDKit❗✔️:   OutputType result = 1;
        // RDKit❗✔️:   return (result << (numAtomPairFingerprintBits +
        // RDKit❗✔️:                      2 * (this->dp_fingerprintArguments->df_includeChirality
        // RDKit❗✔️:                               ? numChiralBits
        // RDKit❗✔️:                               : 0)));
        // RDKit❗✔️: }
        1u64 << (23 + if self.include_chirality { 4 } else { 0 })
    }
}
#[derive(Debug, Clone, Copy)]
pub struct AtomPairPreparedInput<'a> {
    pub topology: &'a TopologyBlock,
    pub coordinates: &'a CoordinateBlock,
    pub properties: &'a MoleculeProperties,
    pub valence: &'a ValenceAssignment,
    pub rings: &'a RingInfo,
    pub use_legacy_stereo_perception: bool,
}
#[derive(Debug, Clone, Copy)]
pub struct AtomPairCall<'a> {
    pub from_atoms: Option<&'a [u32]>,
    pub ignore_atoms: Option<&'a [u32]>,
    pub conformer_id: i32,
    pub custom_atom_invariants: Option<&'a [u32]>,
    pub custom_bond_invariants: Option<&'a [u32]>,
    pub atom_invariants_generator: Option<AtomPairAtomInvariantsGenerator>,
}
impl Default for AtomPairCall<'_> {
    fn default() -> Self {
        Self {
            from_atoms: None,
            ignore_atoms: None,
            conformer_id: -1,
            custom_atom_invariants: None,
            custom_bond_invariants: None,
            atom_invariants_generator: None,
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct AtomPairAtomInvariantsGenerator {
    pub include_chirality: bool,
    pub topological_torsion_correction: bool,
}
impl AtomPairAtomInvariantsGenerator {
    pub(crate) fn from_json(
        &mut self,
        json: &str,
    ) -> Result<(), crate::metadata::FingerprintJsonError> {
        if json.trim().is_empty() {
            return Ok(());
        }
        self.from_json_value(&crate::metadata::parse_object(json)?)
    }
    pub(crate) fn from_json_value(
        &mut self,
        value: &crate::metadata::SourceNode,
    ) -> Result<(), crate::metadata::FingerprintJsonError> {
        // RDKit source: AtomPairGenerator.cpp lines 56-61
        // RDKit✔️✔️: void AtomPairAtomInvGenerator::fromJSON(const boost::property_tree::ptree &pt) {
        // RDKit✔️✔️:   df_includeChirality = pt.get<bool>("includeChirality", df_includeChirality);
        self.include_chirality =
            crate::metadata::bool_or(&value, "includeChirality", self.include_chirality);
        // RDKit✔️✔️:   df_topologicalTorsionCorrection = pt.get<bool>(
        // RDKit✔️✔️:       "topologicalTorsionCorrection", df_topologicalTorsionCorrection);
        self.topological_torsion_correction = crate::metadata::bool_or(
            &value,
            "topologicalTorsionCorrection",
            self.topological_torsion_correction,
        );
        // RDKit✔️✔️:   AtomInvariantsGenerator::fromJSON(pt);
        // RDKit✔️✔️: }
        Ok(())
    }

    #[must_use]
    pub fn info_string(&self) -> String {
        // RDKit source: AtomPairGenerator.cpp lines 45-48
        // RDKit✔️✔️: std::string AtomPairAtomInvGenerator::infoString() const {
        // RDKit✔️✔️:   return "AtomPairInvariantGenerator topologicalTorsionCorrection=" +
        // RDKit✔️✔️:          std::to_string(df_topologicalTorsionCorrection);
        // RDKit✔️✔️: }
        format!(
            "AtomPairInvariantGenerator topologicalTorsionCorrection={}",
            self.topological_torsion_correction as u8
        )
    }

    #[must_use]
    pub fn to_json(&self) -> String {
        // RDKit source: AtomPairGenerator.cpp lines 50-55
        // RDKit✔️✔️: void AtomPairAtomInvGenerator::toJSON(boost::property_tree::ptree &pt) const {
        // RDKit✔️✔️:   pt.put("type", "AtomPairAtomInvGenerator");
        // RDKit✔️✔️:   pt.put("includeChirality", df_includeChirality);
        // RDKit✔️✔️:   pt.put("topologicalTorsionCorrection", df_topologicalTorsionCorrection);
        // RDKit✔️✔️:   AtomInvariantsGenerator::toJSON(pt);
        // RDKit✔️✔️: }
        format!(
            "{{\"type\":\"AtomPairAtomInvGenerator\",\"includeChirality\":\"{}\",\"topologicalTorsionCorrection\":\"{}\"}}",
            self.include_chirality, self.topological_torsion_correction
        )
    }

    pub fn atom_invariants(
        &self,
        input: &AtomPairPreparedInput<'_>,
    ) -> Result<Vec<u32>, AtomPairError> {
        // RDKit source: AtomPairGenerator.cpp lines 31-43
        // RDKit❗✔️: std::vector<std::uint32_t> *AtomPairAtomInvGenerator::getAtomInvariants(
        // RDKit❗✔️:     const ROMol &mol) const {
        // RDKit❗✔️:   auto *atomInvariants = new std::vector<std::uint32_t>(mol.getNumAtoms());
        let mut atom_invariants = vec![0; input.topology.atoms.len()];
        let mut state = AtomCodeInput::from_ref(input.topology, input.properties)
            .map_err(AtomPairError::Topology)?;

        // RDKit❗✔️:   for (ROMol::ConstAtomIterator atomItI = mol.beginAtoms();
        // RDKit❗✔️:        atomItI != mol.endAtoms(); ++atomItI) {
        for atom in &input.topology.atoms {
            // RDKit❗✔️:     (*atomInvariants)[(*atomItI)->getIdx()] =
            // RDKit❗✔️:         getAtomCode(*atomItI, 0, df_includeChirality) -
            // RDKit❗✔️:         (df_topologicalTorsionCorrection ? 2 : 0);
            let correction = if self.topological_torsion_correction {
                2
            } else {
                0
            };
            let assignment = atom_code(
                state,
                Some(atom.id()),
                input
                    .valence
                    .explicit_valence
                    .get(atom.id().index())
                    .copied()
                    .map(|v| v as i8),
                &AtomCodeOptions {
                    branch_subtract: 0,
                    include_chirality: self.include_chirality,
                    use_legacy_stereo_perception: input.use_legacy_stereo_perception,
                },
            )
            .map_err(AtomPairError::AtomCode)?;
            atom_invariants[atom.id().index()] = assignment.code().wrapping_sub(correction);
            state = assignment.into_input();
            // RDKit❗✔️:   }
        }

        // RDKit❗✔️:   return atomInvariants;
        Ok(atom_invariants)
        // RDKit❗✔️: }
    }
}
#[allow(clippy::too_many_arguments)]
fn environments(
    input: &AtomPairPreparedInput<'_>,
    arguments: &AtomPairParams,
    from_atoms: Option<&[u32]>,
    ignore_atoms: Option<&[u32]>,
    conf_id: i32,
) -> Result<Vec<AtomPairEnvironment>, AtomPairError> {
    // RDKit source file: AtomPairGenerator.cpp
    // RDKit source: AtomPairGenerator.cpp lines 177-237
    // RDKit❗✔️: std::vector<AtomEnvironment<OutputType> *>
    // RDKit❗✔️: AtomPairEnvGenerator<OutputType>::getEnvironments(
    // RDKit❗✔️:     const ROMol &mol, FingerprintArguments *arguments,
    // RDKit❗✔️:     const std::vector<std::uint32_t> *fromAtoms,
    // RDKit❗✔️:     const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
    // RDKit❗✔️:     const AdditionalOutput *additionalOutput,
    // RDKit❗✔️:     const std::vector<std::uint32_t> *,  // atomInvariants
    // RDKit❗✔️:     const std::vector<std::uint32_t> *,  // bondInvariants,
    // RDKit❗✔️:     const bool                           // hashResults
    // RDKit❗✔️: ) const {
    // RDKit❗✔️:   const unsigned int atomCount = mol.getNumAtoms();
    let atom_count = input.topology.atoms.len();
    // RDKit❗✔️:   PRECONDITION(!additionalOutput || !additionalOutput->atomToBits ||
    // RDKit❗✔️:                    additionalOutput->atomToBits->size() == atomCount,
    // RDKit❗✔️:                "bad atomToBits size in AdditionalOutput");
    // The shared projector resets every allocated atomToBits vector to the
    // molecule's atom count before calling this family generator.

    // RDKit❗✔️:   auto *atomPairArguments = dynamic_cast<AtomPairArguments *>(arguments);
    // Rust's concrete argument type makes the dynamic-cast precondition structural.
    // RDKit❗✔️:   std::vector<AtomEnvironment<OutputType> *> result =
    // RDKit❗✔️:       std::vector<AtomEnvironment<OutputType> *>();
    let mut result = Vec::new();
    // RDKit❗✔️:   const double *distanceMatrix;
    // RDKit❗✔️:   if (atomPairArguments->df_use2D) {
    let distance_matrix = if arguments.use_2d {
        // RDKit❗✔️:     distanceMatrix = MolOps::getDistanceMat(mol);
        topological_distance_matrix(input.topology, &TopologicalDistanceMatrixParams::default())?
    // RDKit❗✔️:   } else {
    } else {
        // RDKit❗✔️:     distanceMatrix = MolOps::get3DDistanceMat(mol, confId);
        distance_matrix_3d(
            input.topology,
            input.coordinates,
            &DistanceMatrix3dParams {
                conformer_id: (conf_id >= 0).then_some(conf_id as usize),
                use_atom_weights: false,
            },
        )
        .map_err(AtomPairError::Matrix)?
        // RDKit❗✔️:   }
    };

    // RDKit❗✔️:   for (ROMol::ConstAtomIterator atomItI = mol.beginAtoms();
    // RDKit❗✔️:        atomItI != mol.endAtoms(); ++atomItI) {
    for atom_id_first in 0..atom_count {
        // RDKit❗✔️:     unsigned int i = (*atomItI)->getIdx();
        // RDKit❗✔️:     if (ignoreAtoms && std::find(ignoreAtoms->begin(), ignoreAtoms->end(), i) !=
        // RDKit❗✔️:                            ignoreAtoms->end()) {
        if ignore_atoms.is_some_and(|atoms| atoms.contains(&(atom_id_first as u32))) {
            // RDKit❗✔️:       continue;
            continue;
            // RDKit❗✔️:     }
        }

        // RDKit❗✔️:     for (ROMol::ConstAtomIterator atomItJ = atomItI + 1;
        // RDKit❗✔️:          atomItJ != mol.endAtoms(); ++atomItJ) {
        for atom_id_second in (atom_id_first + 1)..atom_count {
            // RDKit❗✔️:       unsigned int j = (*atomItJ)->getIdx();
            // RDKit❗✔️:       if (ignoreAtoms && std::find(ignoreAtoms->begin(), ignoreAtoms->end(),
            // RDKit❗✔️:                                    j) != ignoreAtoms->end()) {
            if ignore_atoms.is_some_and(|atoms| atoms.contains(&(atom_id_second as u32))) {
                // RDKit❗✔️:         continue;
                continue;
                // RDKit❗✔️:       }
            }

            // RDKit❗✔️:       if (fromAtoms &&
            // RDKit❗✔️:           (std::find(fromAtoms->begin(), fromAtoms->end(), i) ==
            // RDKit❗✔️:            fromAtoms->end()) &&
            // RDKit❗✔️:           (std::find(fromAtoms->begin(), fromAtoms->end(), j) ==
            // RDKit❗✔️:            fromAtoms->end())) {
            if from_atoms.is_some_and(|atoms| {
                !atoms.contains(&(atom_id_first as u32))
                    && !atoms.contains(&(atom_id_second as u32))
            }) {
                // RDKit❗✔️:         continue;
                continue;
                // RDKit❗✔️:       }
            }
            // RDKit❗✔️:       auto distance =
            // RDKit❗✔️:           static_cast<unsigned int>(floor(distanceMatrix[i * atomCount + j]));
            let distance = distance_matrix.values()[atom_id_first * atom_count + atom_id_second]
                .floor() as u32;

            // RDKit❗✔️:       if (distance >= atomPairArguments->d_minDistance &&
            // RDKit❗✔️:           distance <= atomPairArguments->d_maxDistance) {
            if distance >= arguments.min_distance && distance <= arguments.max_distance {
                // RDKit❗✔️:         result.push_back(new AtomPairAtomEnv<OutputType>(i, j, distance));
                result.push(AtomPairEnvironment::new(
                    atom_id_first,
                    atom_id_second,
                    distance,
                ));
                // RDKit❗✔️:       }
            }
            // RDKit❗✔️:     }
        }
        // RDKit❗✔️:   }
    }

    // RDKit❗✔️:   return result;
    Ok(result)
    // RDKit❗✔️: }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct AtomPairEnvironment {
    atom_id_first: usize,
    atom_id_second: usize,
    distance: u32,
}

impl AtomPairEnvironment {
    #[must_use]
    pub(crate) const fn new(atom_id_first: usize, atom_id_second: usize, distance: u32) -> Self {
        // RDKit source file: AtomPairGenerator.cpp
        // RDKit source: AtomPairGenerator.cpp lines 170-175
        // RDKit❗✔️: AtomPairAtomEnv<OutputType>::AtomPairAtomEnv(const unsigned int atomIdFirst,
        // RDKit❗✔️:                                              const unsigned int atomIdSecond,
        // RDKit❗✔️:                                              const unsigned int distance)
        // RDKit❗✔️:     : d_atomIdFirst(atomIdFirst),
        // RDKit❗✔️:       d_atomIdSecond(atomIdSecond),
        // RDKit❗✔️:       d_distance(distance) {}
        Self {
            atom_id_first,
            atom_id_second,
            distance,
        }
    }

    pub(crate) fn bit_id(
        &self,
        arguments: &FingerprintArguments,
        atom_invariants: &[u32],
        hash_results: bool,
    ) -> Result<u32, AtomPairError> {
        // RDKit source file: AtomPairGenerator.cpp
        // RDKit source: AtomPairGenerator.cpp lines 131-167
        // RDKit❗✔️: OutputType AtomPairAtomEnv<OutputType>::getBitId(
        // RDKit❗✔️:     FingerprintArguments *arguments,
        // RDKit❗✔️:     const std::vector<std::uint32_t> *atomInvariants,
        // RDKit❗✔️:     const std::vector<std::uint32_t> *,  // bondInvariants
        // RDKit❗✔️:     AdditionalOutput *,                  // additionalOutput,
        // RDKit❗✔️:     const bool hashResults,
        // RDKit❗✔️:     const std::uint64_t  // fpSize
        // RDKit❗✔️: ) const {
        // RDKit❗✔️:   PRECONDITION((atomInvariants->size() >= d_atomIdFirst) &&
        // RDKit❗✔️:                    (atomInvariants->size() >= d_atomIdSecond),
        // RDKit❗✔️:                "bad atom invariants size");
        // The source spells this precondition with `>=`, while the indexed
        // access below requires strict containment. Preserve the effective
        // source boundary without reproducing an out-of-bounds access.
        for atom_index in [self.atom_id_first, self.atom_id_second] {
            if atom_index >= atom_invariants.len() {
                return Err(AtomPairError::AtomInvariantLength {
                    length: atom_invariants.len(),
                    atom_index,
                });
            }
        }

        // RDKit❗✔️:   auto *atomPairArguments = dynamic_cast<AtomPairArguments *>(arguments);
        // The concrete Rust type makes the source dynamic-cast precondition structural.

        // RDKit❗✔️:   std::uint32_t codeSizeLimit =
        // RDKit❗✔️:       (1 << (codeSize +
        // RDKit❗✔️:              (atomPairArguments->df_includeChirality ? numChiralBits : 0))) -
        // RDKit❗✔️:       1;
        let code_size_limit = (1u32 << (9 + if arguments.include_chirality { 2 } else { 0 })) - 1;

        // RDKit❗✔️:   std::uint32_t atomCodeFirst =
        // RDKit❗✔️:       (*atomInvariants)[d_atomIdFirst] % codeSizeLimit;
        let atom_code_first = atom_invariants[self.atom_id_first] % code_size_limit;

        // RDKit❗✔️:   std::uint32_t atomCodeSecond =
        // RDKit❗✔️:       (*atomInvariants)[d_atomIdSecond] % codeSizeLimit;
        let atom_code_second = atom_invariants[self.atom_id_second] % code_size_limit;

        // RDKit❗✔️:   std::uint32_t bitId = 0;
        let mut bit_id = 0u32;
        // RDKit❗✔️:   if (hashResults) {
        if hash_results {
            // RDKit❗✔️:     gboost::hash_combine(bitId, std::min(atomCodeFirst, atomCodeSecond));
            hash_combine(&mut bit_id, atom_code_first.min(atom_code_second));
            // RDKit❗✔️:     gboost::hash_combine(bitId, d_distance);
            hash_combine(&mut bit_id, self.distance);
            // RDKit❗✔️:     gboost::hash_combine(bitId, std::max(atomCodeFirst, atomCodeSecond));
            hash_combine(&mut bit_id, atom_code_first.max(atom_code_second));
        // RDKit❗✔️:   } else {
        } else {
            // RDKit❗✔️:     bitId = getAtomPairCode(atomCodeFirst, atomCodeSecond, d_distance,
            // RDKit❗✔️:                             atomPairArguments->df_includeChirality);
            bit_id = atom_pair_code(
                atom_code_first,
                atom_code_second,
                self.distance,
                arguments.include_chirality,
            )?;
            // RDKit❗✔️:   }
        }

        // RDKit❗✔️:   return bitId;
        Ok(bit_id)
        // RDKit❗✔️: }
    }

    pub(crate) fn update_additional_output(
        &self,
        additional_output: &mut FingerprintAdditionalOutput,
        bit_id: u64,
    ) {
        // RDKit source file: AtomPairGenerator.cpp
        // RDKit source: AtomPairGenerator.cpp lines 109-128
        // RDKit❗✔️: void AtomPairAtomEnv<OutputType>::updateAdditionalOutput(
        // RDKit❗✔️:     AdditionalOutput *additionalOutput, size_t bitId) const {
        // RDKit❗✔️:   PRECONDITION(additionalOutput, "bad output pointer");
        // The Rust reference makes the non-null precondition structural.
        // RDKit❗✔️:   if (additionalOutput->bitInfoMap) {
        if let Some(bit_info_map) = additional_output.bit_info_map.as_mut() {
            // RDKit❗✔️:     (*additionalOutput->bitInfoMap)[bitId].emplace_back(d_atomIdFirst,
            // RDKit❗✔️:                                                         d_atomIdSecond);
            bit_info_map
                .entry(bit_id)
                .or_default()
                .push((self.atom_id_first as u32, self.atom_id_second as u32));
            // RDKit❗✔️:   }
        }
        // RDKit❗✔️:   if (additionalOutput->atomToBits) {
        if let Some(atom_to_bits) = additional_output.atom_to_bits.as_mut() {
            // RDKit❗✔️:     additionalOutput->atomToBits->at(d_atomIdFirst).push_back(bitId);
            atom_to_bits[self.atom_id_first].push(bit_id);
            // RDKit❗✔️:     additionalOutput->atomToBits->at(d_atomIdSecond).push_back(bitId);
            atom_to_bits[self.atom_id_second].push(bit_id);
            // RDKit❗✔️:   }
        }
        // RDKit❗✔️:   if (additionalOutput->atomCounts) {
        if let Some(atom_counts) = additional_output.atom_counts.as_mut() {
            // RDKit❗✔️:     additionalOutput->atomCounts->at(d_atomIdFirst)++;
            atom_counts[self.atom_id_first] = atom_counts[self.atom_id_first].wrapping_add(1);
            // RDKit❗✔️:     additionalOutput->atomCounts->at(d_atomIdSecond)++;
            atom_counts[self.atom_id_second] = atom_counts[self.atom_id_second].wrapping_add(1);
            // RDKit❗✔️:   }
        }
        // RDKit❗✔️:   if (additionalOutput->atomsPerBit) {
        if let Some(atoms_per_bit) = additional_output.atoms_per_bit.as_mut() {
            // RDKit❗✔️:     (*additionalOutput->atomsPerBit)[bitId].push_back(std::vector<int>{
            // RDKit❗✔️:         static_cast<int>(d_atomIdFirst), static_cast<int>(d_atomIdSecond)});
            atoms_per_bit
                .entry(bit_id)
                .or_default()
                .push(vec![self.atom_id_first as i32, self.atom_id_second as i32]);
            // RDKit❗✔️:   }
        }
        // RDKit❗✔️: }
    }
}

impl FingerprintEnvironment<AtomPairError> for AtomPairEnvironment {
    type Output = u32;
    type State = ();
    fn bit_id(
        &self,
        args: &FingerprintArguments,
        atoms: &[u32],
        _bonds: &[u32],
        _output: Option<&mut FingerprintAdditionalOutput>,
        hashed: bool,
        _fp_size: u64,
    ) -> Result<u32, AtomPairError> {
        AtomPairEnvironment::bit_id(self, args, atoms, hashed)
    }
    fn update_output(
        &self,
        output: &mut FingerprintAdditionalOutput,
        bit: u64,
        _state: &mut (),
    ) -> Result<(), AtomPairError> {
        self.update_additional_output(output, bit);
        Ok(())
    }
}
fn count_helper(
    input: &AtomPairPreparedInput<'_>,
    params: &AtomPairParams,
    common: &FingerprintArguments,
    call: &AtomPairCall<'_>,
    fp_size: u64,
    mut output: Option<&mut FingerprintAdditionalOutput>,
) -> Result<SparseCountFingerprint, AtomPairError> {
    let args = FingerprintFuncArguments {
        from_atoms: call.from_atoms,
        ignore_atoms: call.ignore_atoms,
        custom_atom_invariants: call.custom_atom_invariants,
        custom_bond_invariants: call.custom_bond_invariants,
        conformer_id: call.conformer_id,
    };
    with_fingerprint_environment_inputs_and_output(
        input.topology,
        input.properties,
        input.valence,
        input.rings,
        common,
        &args,
        output,
        || {
            call.atom_invariants_generator
                .unwrap_or(AtomPairAtomInvariantsGenerator {
                    include_chirality: params.include_chirality,
                    topological_torsion_correction: false,
                })
                .atom_invariants(input)
        },
        || Ok(Vec::new()),
        |topology, properties, atoms, bonds, output| {
            let prepared_input = AtomPairPreparedInput {
                topology,
                properties,
                ..*input
            };
            let envs = environments(
                &prepared_input,
                params,
                call.from_atoms,
                call.ignore_atoms,
                call.conformer_id,
            )?;
            let mut result = SparseCountFingerprint::new(if fp_size == 0 {
                params.result_size()
            } else {
                fp_size
            });
            accumulate_sparse_counts_into(
                envs,
                common,
                atoms,
                bonds,
                fp_size,
                output,
                &mut result,
            )?;
            Ok(result)
        },
    )
}

pub fn atom_pair_sparse_count(
    input: &AtomPairPreparedInput<'_>,
    params: &AtomPairParams,
    call: &AtomPairCall<'_>,
    output: Option<&mut FingerprintAdditionalOutput>,
) -> Result<SparseCountFingerprint, AtomPairError> {
    let common = params.common()?;
    atom_pair_sparse_count_with_common(input, params, &common, call, output)
}
fn atom_pair_sparse_count_with_common(
    input: &AtomPairPreparedInput<'_>,
    params: &AtomPairParams,
    common: &FingerprintArguments,
    call: &AtomPairCall<'_>,
    output: Option<&mut FingerprintAdditionalOutput>,
) -> Result<SparseCountFingerprint, AtomPairError> {
    count_helper(input, params, common, call, 0, output)
}
pub fn atom_pair_sparse_bits(
    input: &AtomPairPreparedInput<'_>,
    params: &AtomPairParams,
    call: &AtomPairCall<'_>,
    output: Option<&mut FingerprintAdditionalOutput>,
) -> Result<SparseBitFingerprint, AtomPairError> {
    let common = params.common()?;
    atom_pair_sparse_bits_with_common(input, params, &common, call, output)
}
fn atom_pair_sparse_bits_with_common(
    input: &AtomPairPreparedInput<'_>,
    params: &AtomPairParams,
    common: &FingerprintArguments,
    call: &AtomPairCall<'_>,
    output: Option<&mut FingerprintAdditionalOutput>,
) -> Result<SparseBitFingerprint, AtomPairError> {
    project_sparse_fingerprint(
        common,
        params.result_size(),
        input.topology.atoms.len(),
        output,
        |fp_size, output| count_helper(input, params, common, call, fp_size, output),
    )
}
pub fn atom_pair_count(
    input: &AtomPairPreparedInput<'_>,
    params: &AtomPairParams,
    call: &AtomPairCall<'_>,
    output: Option<&mut FingerprintAdditionalOutput>,
) -> Result<SparseCountFingerprint32, AtomPairError> {
    let common = params.common()?;
    atom_pair_count_with_common(input, params, &common, call, output)
}
fn atom_pair_count_with_common(
    input: &AtomPairPreparedInput<'_>,
    params: &AtomPairParams,
    common: &FingerprintArguments,
    call: &AtomPairCall<'_>,
    output: Option<&mut FingerprintAdditionalOutput>,
) -> Result<SparseCountFingerprint32, AtomPairError> {
    project_count_fingerprint(
        common,
        input.topology.atoms.len(),
        output,
        |fp_size, output| count_helper(input, params, common, call, fp_size, output),
    )
}
pub fn atom_pair_bits(
    input: &AtomPairPreparedInput<'_>,
    params: &AtomPairParams,
    call: &AtomPairCall<'_>,
    output: Option<&mut FingerprintAdditionalOutput>,
) -> Result<Fingerprint, AtomPairError> {
    let common = params.common()?;
    atom_pair_bits_with_common(input, params, &common, call, output)
}
fn atom_pair_bits_with_common(
    input: &AtomPairPreparedInput<'_>,
    params: &AtomPairParams,
    common: &FingerprintArguments,
    call: &AtomPairCall<'_>,
    output: Option<&mut FingerprintAdditionalOutput>,
) -> Result<Fingerprint, AtomPairError> {
    project_fingerprint(
        common,
        input.topology.atoms.len(),
        output,
        |fp_size, output| count_helper(input, params, common, call, fp_size, output),
    )
}

#[cfg(test)]
mod tests {
    #[test]
    fn project_batch_preflight_preserves_legacy_rejection_order_and_scalar_zero_size() {
        let zero = AtomPairParams {
            fp_size: 0,
            ..Default::default()
        };
        assert_eq!(
            validate_atom_pair_params(&zero).unwrap_err().to_string(),
            "fingerprint requires n_bits > 0"
        );
        assert_eq!(zero.common().unwrap().fp_size, 0);
        let distance_and_bits = AtomPairParams {
            min_distance: 31,
            max_distance: 30,
            bits_per_feature: 0,
            ..Default::default()
        };
        assert_eq!(
            validate_atom_pair_params(&distance_and_bits),
            Err(FingerprintError::PreconditionViolation {
                what: "bad distances provided"
            })
        );
        let bounds_and_distance = AtomPairParams {
            count_bounds: vec![],
            min_distance: 31,
            max_distance: 30,
            bits_per_feature: 0,
            ..Default::default()
        };
        assert_eq!(
            validate_atom_pair_params(&bounds_and_distance),
            Err(FingerprintError::PreconditionViolation {
                what: "bad count bounds provided"
            })
        );
        let bits = AtomPairParams {
            bits_per_feature: 0,
            ..Default::default()
        };
        assert_eq!(
            validate_atom_pair_params(&bits),
            Err(FingerprintError::InvalidArguments {
                reason: "num_bits_per_feature must be greater than zero"
            })
        );
        assert!(
            validate_atom_pair_params(&AtomPairParams {
                count_simulation: false,
                count_bounds: vec![],
                ..Default::default()
            })
            .is_ok()
        );
        let input = TestBuilder::default().build().unwrap();
        // getFingerprint's source count-simulation branch rejects bounds.len() >= fpSize;
        // the defined zero-size no-simulation branch returns empty for no environments.
        assert!(matches!(
            atom_pair_bits(&input.input(), &zero, &AtomPairCall::default(), None),
            Err(AtomPairError::Fingerprint(
                FingerprintError::InvalidArguments {
                    reason: "Count bounds size is >= fingerprint size"
                }
            ))
        ));
        let scalar_zero = AtomPairParams {
            count_simulation: false,
            ..zero.clone()
        };
        assert_eq!(
            atom_pair_bits(&input.input(), &scalar_zero, &AtomPairCall::default(), None)
                .unwrap()
                .n_bits(),
            0
        );
        assert_eq!(
            atom_pair_count(&input.input(), &zero, &AtomPairCall::default(), None)
                .unwrap()
                .length(),
            0
        );
    }

    use super::*;
    #[test]
    fn environment_bit_id_exact_mode_matches_pair_packing_and_endpoint_reversal() {
        let arguments = AtomPairParams::default();
        let forward = AtomPairEnvironment::new(0, 1, 1);
        let reverse = AtomPairEnvironment::new(1, 0, 1);
        let invariants = [10, 20];
        assert_eq!(
            forward
                .bit_id(&arguments.common().unwrap(), &invariants, false)
                .unwrap(),
            328_001
        );
        assert_eq!(
            reverse
                .bit_id(&arguments.common().unwrap(), &invariants, false)
                .unwrap(),
            328_001
        );
    }

    #[test]
    fn environment_bit_id_hashed_mode_uses_canonical_three_stage_gboost_order() {
        let arguments = AtomPairParams::default();
        let forward = AtomPairEnvironment::new(0, 1, 1);
        let reverse = AtomPairEnvironment::new(1, 0, 1);
        let invariants = [10, 20];
        assert_eq!(
            forward
                .bit_id(&arguments.common().unwrap(), &invariants, true)
                .unwrap(),
            4_217_127_294
        );
        assert_eq!(
            reverse
                .bit_id(&arguments.common().unwrap(), &invariants, true)
                .unwrap(),
            4_217_127_294
        );
    }

    #[test]
    fn environment_bit_id_uses_source_modulo_limit_not_a_bit_mask() {
        let unchiral = AtomPairParams::default();
        let environment = AtomPairEnvironment::new(0, 1, 1);
        assert_eq!(
            environment
                .bit_id(&unchiral.common().unwrap(), &[510, 511], false)
                .unwrap(),
            atom_pair_code(510, 0, 1, false).unwrap()
        );
        assert_eq!(
            environment
                .bit_id(&unchiral.common().unwrap(), &[511, 0], false)
                .unwrap(),
            environment
                .bit_id(&unchiral.common().unwrap(), &[0, 0], false)
                .unwrap()
        );
        assert_eq!(
            environment
                .bit_id(&unchiral.common().unwrap(), &[2047, u32::MAX], false,)
                .unwrap(),
            atom_pair_code(3, 31, 1, false).unwrap()
        );

        let mut chiral = AtomPairParams::default();
        chiral.include_chirality = true;
        assert_eq!(
            environment
                .bit_id(&chiral.common().unwrap(), &[2047, u32::MAX], false)
                .unwrap(),
            atom_pair_code(0, 1023, 1, true).unwrap()
        );
    }

    #[test]
    fn environment_bit_id_preserves_modulo_collisions_before_hashing() {
        let arguments = AtomPairParams::default();
        let environment = AtomPairEnvironment::new(0, 1, 1);
        assert_eq!(
            environment
                .bit_id(&arguments.common().unwrap(), &[0, 0], true)
                .unwrap(),
            environment
                .bit_id(&arguments.common().unwrap(), &[511, 511], true)
                .unwrap()
        );
        assert_eq!(
            environment
                .bit_id(&arguments.common().unwrap(), &[0, 0], true)
                .unwrap(),
            4_216_857_148
        );
    }

    #[test]
    fn additional_output_each_supported_allocation_receives_the_source_shape() {
        let environment = AtomPairEnvironment::new(0, 2, 3);

        let mut bit_info = FingerprintAdditionalOutput::default();
        bit_info.allocate_bit_info_map();
        bit_info.reinitialize(3);
        environment.update_additional_output(&mut bit_info, 17);
        assert_eq!(bit_info.bit_info_map.unwrap().get(&17).unwrap(), &[(0, 2)]);

        let mut atom_to_bits = FingerprintAdditionalOutput::default();
        atom_to_bits.allocate_atom_to_bits();
        atom_to_bits.reinitialize(3);
        environment.update_additional_output(&mut atom_to_bits, 17);
        assert_eq!(
            atom_to_bits.atom_to_bits.unwrap(),
            [vec![17], vec![], vec![17]]
        );

        let mut atom_counts = FingerprintAdditionalOutput::default();
        atom_counts.allocate_atom_counts();
        atom_counts.reinitialize(3);
        environment.update_additional_output(&mut atom_counts, 17);
        assert_eq!(atom_counts.atom_counts.unwrap(), [1, 0, 1]);

        let mut atoms_per_bit = FingerprintAdditionalOutput::default();
        atoms_per_bit.allocate_atoms_per_bit();
        atoms_per_bit.reinitialize(3);
        environment.update_additional_output(&mut atoms_per_bit, 17);
        assert_eq!(
            atoms_per_bit.atoms_per_bit.unwrap().get(&17).unwrap(),
            &[vec![0, 2]]
        );
    }

    #[test]
    fn additional_output_combined_allocations_preserve_pair_and_call_order() {
        let mut output = FingerprintAdditionalOutput::default();
        output.allocate_bit_info_map();
        output.allocate_atom_to_bits();
        output.allocate_atom_counts();
        output.allocate_atoms_per_bit();
        output.allocate_bit_paths();
        output.reinitialize(4);

        AtomPairEnvironment::new(0, 2, 3).update_additional_output(&mut output, 17);
        AtomPairEnvironment::new(1, 3, 2).update_additional_output(&mut output, 19);
        AtomPairEnvironment::new(2, 0, 3).update_additional_output(&mut output, 17);

        assert_eq!(
            output.bit_info_map.as_ref().unwrap().get(&17).unwrap(),
            &[(0, 2), (2, 0)]
        );
        assert_eq!(
            output.bit_info_map.as_ref().unwrap().get(&19).unwrap(),
            &[(1, 3)]
        );
        assert_eq!(
            output.atom_to_bits.as_ref().unwrap(),
            &[vec![17, 17], vec![19], vec![17, 17], vec![19]]
        );
        assert_eq!(output.atom_counts.as_ref().unwrap(), &[2, 1, 2, 1]);
        assert_eq!(
            output.atoms_per_bit.as_ref().unwrap().get(&17).unwrap(),
            &[vec![0, 2], vec![2, 0]]
        );
        assert!(output.bit_paths.as_ref().unwrap().is_empty());
    }

    #[test]
    fn invalid_invariant_lengths_retain_selected_index_and_length() {
        let args = AtomPairParams::default().common().unwrap();
        for (first, second, atoms, expected) in [
            (0, 0, Vec::new(), (0, 0)),
            (0, 1, vec![10], (1, 1)),
            (1, 0, vec![10], (1, 1)),
        ] {
            let error = AtomPairEnvironment::new(first, second, 1)
                .bit_id(&args, &atoms, false)
                .unwrap_err();
            assert!(
                matches!(error,AtomPairError::AtomInvariantLength { length, atom_index } if (length,atom_index)==expected)
            );
        }
    }

    use crate::test_support::{TestBuilder, TestMolecule};
    use cosmolkit_model::{AtomSpec, BondOrder, BondSpec, Conformer3D, Element};
    fn environments_2d(
        molecule: &TestMolecule,
        arguments: &AtomPairParams,
        from_atoms: Option<&[u32]>,
        ignore_atoms: Option<&[u32]>,
    ) -> Vec<AtomPairEnvironment> {
        environments(&molecule.input(), arguments, from_atoms, ignore_atoms, -1).unwrap()
    }
    fn three_atom_3d_molecule() -> TestMolecule {
        let mut builder = TestBuilder::default();
        for _ in 0..3 {
            builder.add_atom(AtomSpec::new(Element::C));
        }
        builder
            .add_conformer(Conformer3D::new(
                7,
                vec![[0.0, 0.0, 0.0], [1.999, 0.0, 0.0], [2.001, 0.0, 0.0]],
                true,
            ))
            .unwrap();
        builder
            .add_conformer(Conformer3D::new(
                2,
                vec![[0.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]],
                true,
            ))
            .unwrap();
        builder.build().unwrap()
    }

    #[test]
    fn additional_output_preserves_duplicate_provenance_after_hash_collision() {
        let arguments = AtomPairParams::default();
        let first = AtomPairEnvironment::new(0, 1, 1);
        let second = AtomPairEnvironment::new(2, 3, 1);
        let invariants = [0, 0, 511, 511];
        let first_bit = first
            .bit_id(&arguments.common().unwrap(), &invariants, true)
            .unwrap();
        let second_bit = second
            .bit_id(&arguments.common().unwrap(), &invariants, true)
            .unwrap();
        assert_eq!(first_bit, second_bit);

        let mut output = FingerprintAdditionalOutput::default();
        output.allocate_bit_info_map();
        output.allocate_atom_to_bits();
        output.allocate_atom_counts();
        output.allocate_atoms_per_bit();
        output.reinitialize(4);
        first.update_additional_output(&mut output, u64::from(first_bit));
        second.update_additional_output(&mut output, u64::from(second_bit));

        let bit_id = u64::from(first_bit);
        assert_eq!(
            output.bit_info_map.as_ref().unwrap().get(&bit_id).unwrap(),
            &[(0, 1), (2, 3)]
        );
        assert_eq!(
            output.atoms_per_bit.as_ref().unwrap().get(&bit_id).unwrap(),
            &[vec![0, 1], vec![2, 3]]
        );
        assert_eq!(output.atom_counts.as_ref().unwrap(), &[1, 1, 1, 1]);
    }

    #[test]
    fn environment_generation_empty_single_chain_branch_ring_and_fused_order() {
        let arguments = AtomPairParams::default();
        assert!(environments_2d(&TestMolecule::new(), &arguments, None, None).is_empty());
        assert!(
            environments_2d(
                &TestMolecule::from_smiles("[He]").unwrap(),
                &arguments,
                None,
                None
            )
            .is_empty()
        );

        let chain = TestMolecule::from_smiles("CCCC").unwrap();
        assert_eq!(
            environments_2d(&chain, &arguments, None, None),
            [
                AtomPairEnvironment::new(0, 1, 1),
                AtomPairEnvironment::new(0, 2, 2),
                AtomPairEnvironment::new(0, 3, 3),
                AtomPairEnvironment::new(1, 2, 1),
                AtomPairEnvironment::new(1, 3, 2),
                AtomPairEnvironment::new(2, 3, 1),
            ]
        );

        let branch = TestMolecule::from_smiles("CC(C)C").unwrap();
        assert_eq!(
            environments_2d(&branch, &arguments, None, None),
            [
                AtomPairEnvironment::new(0, 1, 1),
                AtomPairEnvironment::new(0, 2, 2),
                AtomPairEnvironment::new(0, 3, 2),
                AtomPairEnvironment::new(1, 2, 1),
                AtomPairEnvironment::new(1, 3, 1),
                AtomPairEnvironment::new(2, 3, 2),
            ]
        );

        let ring = TestMolecule::from_smiles("C1CC1").unwrap();
        assert_eq!(
            environments_2d(&ring, &arguments, None, None),
            [
                AtomPairEnvironment::new(0, 1, 1),
                AtomPairEnvironment::new(0, 2, 1),
                AtomPairEnvironment::new(1, 2, 1),
            ]
        );

        let fused = TestMolecule::from_smiles("c1ccc2ccccc2c1").unwrap();
        let fused_environments = environments_2d(&fused, &arguments, None, None);
        assert_eq!(fused_environments.len(), 45);
        assert!(fused_environments.windows(2).all(|pair| {
            (pair[0].atom_id_first, pair[0].atom_id_second)
                < (pair[1].atom_id_first, pair[1].atom_id_second)
        }));
    }

    #[test]
    fn environment_generation_includes_explicit_hydrogens_and_handles_disconnected_distance() {
        let arguments = AtomPairParams::default();
        let mut explicit_h_builder = TestBuilder::default();
        let carbon = explicit_h_builder.add_atom(AtomSpec::new(Element::C));
        for _ in 0..4 {
            let hydrogen = explicit_h_builder.add_atom(AtomSpec::new(Element::H));
            explicit_h_builder
                .add_bond(BondSpec::new(carbon, hydrogen, BondOrder::Single))
                .unwrap();
        }
        let explicit_h = explicit_h_builder.build().unwrap();
        assert_eq!(explicit_h.num_atoms(), 5);
        assert_eq!(
            environments_2d(&explicit_h, &arguments, None, None).len(),
            10
        );

        let disconnected = TestMolecule::from_smiles("C.C").unwrap();
        assert!(environments_2d(&disconnected, &arguments, None, None).is_empty());
        let mut wide = arguments;
        wide.max_distance = 100_000_000;
        assert_eq!(
            environments_2d(&disconnected, &wide, None, None),
            [AtomPairEnvironment::new(0, 1, 100_000_000)]
        );
    }

    #[test]
    fn environment_generation_root_and_ignore_filters_match_source_precedence() {
        let molecule = TestMolecule::from_smiles("CCCC").unwrap();
        let arguments = AtomPairParams::default();
        assert_eq!(
            environments_2d(&molecule, &arguments, Some(&[1]), None),
            [
                AtomPairEnvironment::new(0, 1, 1),
                AtomPairEnvironment::new(1, 2, 1),
                AtomPairEnvironment::new(1, 3, 2),
            ]
        );
        assert_eq!(
            environments_2d(&molecule, &arguments, Some(&[1, 1, 99]), None),
            environments_2d(&molecule, &arguments, Some(&[1]), None)
        );
        assert!(environments_2d(&molecule, &arguments, Some(&[]), None).is_empty());
        assert!(environments_2d(&molecule, &arguments, Some(&[99]), None).is_empty());

        assert_eq!(
            environments_2d(&molecule, &arguments, None, Some(&[1])),
            [
                AtomPairEnvironment::new(0, 2, 2),
                AtomPairEnvironment::new(0, 3, 3),
                AtomPairEnvironment::new(2, 3, 1),
            ]
        );
        assert_eq!(
            environments_2d(&molecule, &arguments, None, Some(&[1, 1, 99])),
            environments_2d(&molecule, &arguments, None, Some(&[1]))
        );
        assert!(environments_2d(&molecule, &arguments, Some(&[1]), Some(&[1])).is_empty());
        assert_eq!(
            environments_2d(&molecule, &arguments, Some(&[0, 2]), Some(&[0])),
            [
                AtomPairEnvironment::new(1, 2, 1),
                AtomPairEnvironment::new(2, 3, 1),
            ]
        );
    }

    #[test]
    fn environment_generation_distance_bounds_are_inclusive_and_allow_equality() {
        let molecule = TestMolecule::from_smiles("CCCCC").unwrap();
        let exactly_two = AtomPairParams {
            count_simulation: false,
            min_distance: 2,
            max_distance: 2,
            count_bounds: Vec::new(),
            ..Default::default()
        };
        assert_eq!(
            environments_2d(&molecule, &exactly_two, None, None),
            [
                AtomPairEnvironment::new(0, 2, 2),
                AtomPairEnvironment::new(1, 3, 2),
                AtomPairEnvironment::new(2, 4, 2),
            ]
        );

        let one_through_two = AtomPairParams {
            count_simulation: false,
            max_distance: 2,
            count_bounds: Vec::new(),
            ..Default::default()
        };
        let distances = environments_2d(&molecule, &one_through_two, None, None)
            .into_iter()
            .map(|environment| environment.distance)
            .collect::<Vec<_>>();
        assert_eq!(distances, [1, 2, 1, 2, 1, 2, 1]);
    }

    #[test]
    fn environment_generation_3d_floors_fractional_distances_and_selects_conformer_ids() {
        let molecule = three_atom_3d_molecule();
        let arguments = AtomPairParams {
            count_simulation: false,
            use_2d: false,
            min_distance: 0,
            max_distance: 10,
            count_bounds: Vec::new(),
            ..Default::default()
        };
        assert_eq!(
            environments(&molecule.input(), &arguments, None, None, -1).unwrap(),
            [
                AtomPairEnvironment::new(0, 1, 1),
                AtomPairEnvironment::new(0, 2, 2),
                AtomPairEnvironment::new(1, 2, 0),
            ]
        );
        assert_eq!(
            environments(&molecule.input(), &arguments, None, None, 2).unwrap(),
            [
                AtomPairEnvironment::new(0, 1, 3),
                AtomPairEnvironment::new(0, 2, 4),
                AtomPairEnvironment::new(1, 2, 5),
            ]
        );
    }
}
