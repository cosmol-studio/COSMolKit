//! RDKit AtomPair fingerprint source reproduction.
//!
//! This is the sole owner of AtomPair-specific encoding, invariant, and
//! environment behavior. Generic vector projection remains in `generator`.

use crate::chemistry::atom_properties::{AtomPropertyError, num_pi_electrons};
use crate::chemistry::ciplabeler::assign_cip_labels;
use crate::chemistry::matrices::{DistanceMatrixError, distance_matrix_3d, topological_distance_matrix};
use crate::{AtomId, ChiralTag, Molecule};

use super::{
    AdditionalOutput, Fingerprint, FingerprintArguments, FingerprintError, FingerprintFuncArguments,
    SparseBitFingerprint, SparseCountFingerprint, generator, hash_combine, json_value_as_bool, json_value_as_u32,
    rdkit_use_legacy_stereo_perception,
};
use serde_json::Value;

// RDKit source file: FingerprintUtil.h
// RDKit source: FingerprintUtil.h lines 27-41
// RDKit✔️✔️: namespace AtomPairs {
// RDKit✔️✔️: const unsigned int numTypeBits = 4;
pub(crate) const NUM_TYPE_BITS: u32 = 4;
// RDKit✔️✔️: const unsigned int atomNumberTypes[1 << numTypeBits] = {
// RDKit✔️✔️:     5, 6, 7, 8, 9, 14, 15, 16, 17, 33, 34, 35, 51, 52, 53};
//
// C++ aggregate initialization zero-initializes the unlisted sixteenth slot.
pub(crate) const ATOM_NUMBER_TYPES: [u32; 1usize << NUM_TYPE_BITS] =
    [5, 6, 7, 8, 9, 14, 15, 16, 17, 33, 34, 35, 51, 52, 53, 0];
// RDKit✔️✔️: const unsigned int numPiBits = 2;
pub(crate) const NUM_PI_BITS: u32 = 2;
// RDKit✔️✔️: const unsigned int maxNumPi = (1 << numPiBits) - 1;
pub(crate) const MAX_NUM_PI: u32 = (1u32 << NUM_PI_BITS) - 1;
// RDKit✔️✔️: const unsigned int numBranchBits = 3;
pub(crate) const NUM_BRANCH_BITS: u32 = 3;
// RDKit✔️✔️: const unsigned int maxNumBranches = (1 << numBranchBits) - 1;
pub(crate) const MAX_NUM_BRANCHES: u32 = (1u32 << NUM_BRANCH_BITS) - 1;
// RDKit✔️✔️: const unsigned int numChiralBits = 2;
pub(crate) const NUM_CHIRAL_BITS: u32 = 2;
// RDKit✔️✔️: const unsigned int codeSize = numTypeBits + numPiBits + numBranchBits;
pub(crate) const CODE_SIZE: u32 = NUM_TYPE_BITS + NUM_PI_BITS + NUM_BRANCH_BITS;
// RDKit✔️✔️: const unsigned int numPathBits = 5;
pub(crate) const NUM_PATH_BITS: u32 = 5;
// RDKit✔️✔️: const unsigned int maxPathLen = (1 << numPathBits) - 1;
pub(crate) const MAX_PATH_LEN: u32 = (1u32 << NUM_PATH_BITS) - 1;
// RDKit✔️✔️: const unsigned int numAtomPairFingerprintBits =
// RDKit✔️✔️:     numPathBits + 2 * codeSize;  // note that this is only accurate if chirality
// RDKit✔️✔️:                                  // is not included
pub(crate) const NUM_ATOM_PAIR_FINGERPRINT_BITS: u32 = NUM_PATH_BITS + 2 * CODE_SIZE;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub(crate) struct AtomPairEnvironmentGenerator;

impl AtomPairEnvironmentGenerator {
    #[must_use]
    pub(crate) fn result_size(&self, arguments: &AtomPairArguments) -> u64 {
        // RDKit source file: AtomPairGenerator.cpp
        // RDKit source: AtomPairGenerator.cpp lines 68-75
        // RDKit✔️✔️: OutputType AtomPairEnvGenerator<OutputType>::getResultSize() const {
        // RDKit✔️✔️:   OutputType result = 1;
        let result = 1u64;
        // RDKit✔️✔️:   return (result << (numAtomPairFingerprintBits +
        // RDKit✔️✔️:                      2 * (this->dp_fingerprintArguments->df_includeChirality
        // RDKit✔️✔️:                               ? numChiralBits
        // RDKit✔️✔️:                               : 0)));
        result
            << (NUM_ATOM_PAIR_FINGERPRINT_BITS
                + 2 * if arguments.fingerprint_arguments.df_include_chirality {
                    NUM_CHIRAL_BITS
                } else {
                    0
                })
        // RDKit✔️✔️: }
    }

    #[must_use]
    pub(crate) const fn info_string(&self) -> &'static str {
        // RDKit source: AtomPairGenerator.cpp lines 239-242
        // RDKit✔️✔️: std::string AtomPairEnvGenerator<OutputType>::infoString() const {
        // RDKit✔️✔️:   return "AtomPairEnvironmentGenerator";
        // RDKit✔️✔️: }
        "AtomPairEnvironmentGenerator"
    }

    #[must_use]
    pub(crate) const fn to_json(&self) -> &'static str {
        // RDKit source: AtomPairGenerator.cpp lines 244-249
        // RDKit✔️✔️: void AtomPairEnvGenerator<OutputType>::toJSON(
        // RDKit✔️✔️:     boost::property_tree::ptree &pt) const {
        // RDKit✔️✔️:   pt.put("type", "AtomPairEnvGenerator");
        // RDKit✔️✔️:   AtomEnvironmentGenerator<OutputType>::toJSON(pt);
        // RDKit✔️✔️: }
        r#"{"type":"AtomPairEnvGenerator"}"#
    }

    pub(crate) fn from_json(&mut self, json: &str) -> Result<(), FingerprintError> {
        // AtomPairEnvGenerator inherits the source base implementation, whose
        // state is empty. Restoration therefore validates the property-tree
        // object while retaining the same stateless generator value.
        if json.trim().is_empty() {
            return Ok(());
        }
        let value: Value =
            serde_json::from_str(json).map_err(|error| FingerprintError::InvalidArgumentsJson(error.to_string()))?;
        value
            .as_object()
            .ok_or_else(|| FingerprintError::InvalidArgumentsJson("expected JSON object".to_string()))?;
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn environments(
        &self,
        molecule: &Molecule,
        arguments: &AtomPairArguments,
        from_atoms: Option<&[usize]>,
        ignore_atoms: Option<&[usize]>,
        conf_id: i32,
    ) -> Result<Vec<AtomPairEnvironment>, AtomPairCodeError> {
        // RDKit source file: AtomPairGenerator.cpp
        // RDKit source: AtomPairGenerator.cpp lines 177-237
        // RDKit✔️✔️: std::vector<AtomEnvironment<OutputType> *>
        // RDKit✔️✔️: AtomPairEnvGenerator<OutputType>::getEnvironments(
        // RDKit✔️✔️:     const ROMol &mol, FingerprintArguments *arguments,
        // RDKit✔️✔️:     const std::vector<std::uint32_t> *fromAtoms,
        // RDKit✔️✔️:     const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
        // RDKit✔️✔️:     const AdditionalOutput *additionalOutput,
        // RDKit✔️✔️:     const std::vector<std::uint32_t> *,  // atomInvariants
        // RDKit✔️✔️:     const std::vector<std::uint32_t> *,  // bondInvariants,
        // RDKit✔️✔️:     const bool                           // hashResults
        // RDKit✔️✔️: ) const {
        // RDKit✔️✔️:   const unsigned int atomCount = mol.getNumAtoms();
        let atom_count = molecule.num_atoms();
        // RDKit✔️✔️:   PRECONDITION(!additionalOutput || !additionalOutput->atomToBits ||
        // RDKit✔️✔️:                    additionalOutput->atomToBits->size() == atomCount,
        // RDKit✔️✔️:                "bad atomToBits size in AdditionalOutput");
        // The shared projector resets every allocated atomToBits vector to the
        // molecule's atom count before calling this family generator.

        // RDKit✔️✔️:   auto *atomPairArguments = dynamic_cast<AtomPairArguments *>(arguments);
        // Rust's concrete argument type makes the dynamic-cast precondition structural.
        // RDKit✔️✔️:   std::vector<AtomEnvironment<OutputType> *> result =
        // RDKit✔️✔️:       std::vector<AtomEnvironment<OutputType> *>();
        let mut result = Vec::new();
        // RDKit✔️✔️:   const double *distanceMatrix;
        // RDKit✔️✔️:   if (atomPairArguments->df_use2D) {
        let distance_matrix = if arguments.use_2d {
            // RDKit✔️✔️:     distanceMatrix = MolOps::getDistanceMat(mol);
            topological_distance_matrix(molecule)
        // RDKit✔️✔️:   } else {
        } else {
            // RDKit✔️✔️:     distanceMatrix = MolOps::get3DDistanceMat(mol, confId);
            distance_matrix_3d(molecule, conf_id, false).map_err(AtomPairCodeError::DistanceMatrix)?
            // RDKit✔️✔️:   }
        };

        // RDKit✔️✔️:   for (ROMol::ConstAtomIterator atomItI = mol.beginAtoms();
        // RDKit✔️✔️:        atomItI != mol.endAtoms(); ++atomItI) {
        for atom_id_first in 0..atom_count {
            // RDKit✔️✔️:     unsigned int i = (*atomItI)->getIdx();
            // RDKit✔️✔️:     if (ignoreAtoms && std::find(ignoreAtoms->begin(), ignoreAtoms->end(), i) !=
            // RDKit✔️✔️:                            ignoreAtoms->end()) {
            if ignore_atoms.is_some_and(|atoms| atoms.contains(&atom_id_first)) {
                // RDKit✔️✔️:       continue;
                continue;
                // RDKit✔️✔️:     }
            }

            // RDKit✔️✔️:     for (ROMol::ConstAtomIterator atomItJ = atomItI + 1;
            // RDKit✔️✔️:          atomItJ != mol.endAtoms(); ++atomItJ) {
            for atom_id_second in (atom_id_first + 1)..atom_count {
                // RDKit✔️✔️:       unsigned int j = (*atomItJ)->getIdx();
                // RDKit✔️✔️:       if (ignoreAtoms && std::find(ignoreAtoms->begin(), ignoreAtoms->end(),
                // RDKit✔️✔️:                                    j) != ignoreAtoms->end()) {
                if ignore_atoms.is_some_and(|atoms| atoms.contains(&atom_id_second)) {
                    // RDKit✔️✔️:         continue;
                    continue;
                    // RDKit✔️✔️:       }
                }

                // RDKit✔️✔️:       if (fromAtoms &&
                // RDKit✔️✔️:           (std::find(fromAtoms->begin(), fromAtoms->end(), i) ==
                // RDKit✔️✔️:            fromAtoms->end()) &&
                // RDKit✔️✔️:           (std::find(fromAtoms->begin(), fromAtoms->end(), j) ==
                // RDKit✔️✔️:            fromAtoms->end())) {
                if from_atoms.is_some_and(|atoms| !atoms.contains(&atom_id_first) && !atoms.contains(&atom_id_second)) {
                    // RDKit✔️✔️:         continue;
                    continue;
                    // RDKit✔️✔️:       }
                }
                // RDKit✔️✔️:       auto distance =
                // RDKit✔️✔️:           static_cast<unsigned int>(floor(distanceMatrix[i * atomCount + j]));
                let distance = distance_matrix[atom_id_first * atom_count + atom_id_second].floor() as u32;

                // RDKit✔️✔️:       if (distance >= atomPairArguments->d_minDistance &&
                // RDKit✔️✔️:           distance <= atomPairArguments->d_maxDistance) {
                if distance >= arguments.min_distance && distance <= arguments.max_distance {
                    // RDKit✔️✔️:         result.push_back(new AtomPairAtomEnv<OutputType>(i, j, distance));
                    result.push(AtomPairEnvironment::new(atom_id_first, atom_id_second, distance));
                    // RDKit✔️✔️:       }
                }
                // RDKit✔️✔️:     }
            }
            // RDKit✔️✔️:   }
        }

        // RDKit✔️✔️:   return result;
        Ok(result)
        // RDKit✔️✔️: }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomPairArguments {
    pub fingerprint_arguments: FingerprintArguments,
    pub use_2d: bool,
    pub min_distance: u32,
    pub max_distance: u32,
}

impl Default for AtomPairArguments {
    fn default() -> Self {
        Self {
            fingerprint_arguments: FingerprintArguments {
                df_count_simulation: true,
                df_include_chirality: false,
                d_count_bounds: vec![1, 2, 4, 8],
                d_fp_size: 2048,
                d_num_bits_per_feature: 1,
            },
            use_2d: true,
            min_distance: 1,
            max_distance: MAX_PATH_LEN - 1,
        }
    }
}

impl AtomPairArguments {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        count_simulation: bool,
        include_chirality: bool,
        use_2d: bool,
        min_distance: u32,
        max_distance: u32,
        count_bounds: Vec<u32>,
        fp_size: u32,
    ) -> Result<Self, FingerprintError> {
        // RDKit source file: AtomPairGenerator.cpp
        // RDKit source: AtomPairGenerator.cpp lines 77-87
        // RDKit✔️✔️: AtomPairArguments::AtomPairArguments(
        // RDKit✔️✔️:     const bool countSimulation, const bool includeChirality, const bool use2D,
        // RDKit✔️✔️:     const unsigned int minDistance, const unsigned int maxDistance,
        // RDKit✔️✔️:     const std::vector<std::uint32_t> countBounds, const std::uint32_t fpSize)
        // RDKit✔️✔️:     : FingerprintArguments(countSimulation, countBounds, fpSize, 1,
        // RDKit✔️✔️:                            includeChirality),
        let fingerprint_arguments =
            FingerprintArguments::new(count_simulation, count_bounds, fp_size, 1, include_chirality)?;
        // RDKit✔️✔️:       df_use2D(use2D),
        // RDKit✔️✔️:       d_minDistance(minDistance),
        // RDKit✔️✔️:       d_maxDistance(maxDistance) {
        // RDKit✔️✔️:   PRECONDITION(minDistance <= maxDistance, "bad distances provided");
        if min_distance > max_distance {
            return Err(FingerprintError::InvalidArguments {
                reason: "minDistance must be less than or equal to maxDistance",
            });
        }
        // RDKit✔️✔️: }
        Ok(Self {
            fingerprint_arguments,
            use_2d,
            min_distance,
            max_distance,
        })
    }

    #[must_use]
    pub fn info_string(&self) -> String {
        // RDKit source: AtomPairGenerator.cpp lines 89-93
        // RDKit✔️✔️: std::string AtomPairArguments::infoString() const {
        // RDKit✔️✔️:   return "AtomPairArguments use2D=" + std::to_string(df_use2D) +
        // RDKit✔️✔️:          " minDistance=" + std::to_string(d_minDistance) +
        // RDKit✔️✔️:          " maxDistance=" + std::to_string(d_maxDistance);
        // RDKit✔️✔️: }
        format!(
            "AtomPairArguments use2D={} minDistance={} maxDistance={}",
            self.use_2d as u8, self.min_distance, self.max_distance
        )
    }

    #[must_use]
    pub fn to_json(&self) -> String {
        // RDKit source: AtomPairGenerator.cpp lines 94-100
        // RDKit✔️✔️: void AtomPairArguments::toJSON(boost::property_tree::ptree &pt) const {
        // RDKit✔️✔️:   pt.put("type", "AtomPairArguments");
        // RDKit✔️✔️:   pt.put("use2D", df_use2D);
        // RDKit✔️✔️:   pt.put("minDistance", d_minDistance);
        // RDKit✔️✔️:   pt.put("maxDistance", d_maxDistance);
        // RDKit✔️✔️:   FingerprintArguments::toJSON(pt);
        // RDKit✔️✔️: }
        let common = self.fingerprint_arguments.to_json();
        let common_body = &common[1..common.len().saturating_sub(1)];
        format!(
            "{{\"type\":\"AtomPairArguments\",\"use2D\":\"{}\",\"minDistance\":\"{}\",\"maxDistance\":\"{}\",{common_body}}}",
            self.use_2d, self.min_distance, self.max_distance
        )
    }

    pub fn from_json(&mut self, json: &str) -> Result<(), FingerprintError> {
        if json.trim().is_empty() {
            return Ok(());
        }
        let value: Value =
            serde_json::from_str(json).map_err(|error| FingerprintError::InvalidArgumentsJson(error.to_string()))?;
        self.from_json_value(&value)
    }

    fn from_json_value(&mut self, value: &Value) -> Result<(), FingerprintError> {
        // RDKit source: AtomPairGenerator.cpp lines 101-106
        // RDKit✔️✔️: void AtomPairArguments::fromJSON(const boost::property_tree::ptree &pt) {
        let object = value
            .as_object()
            .ok_or_else(|| FingerprintError::InvalidArgumentsJson("expected JSON object".to_string()))?;
        // RDKit✔️✔️:   df_use2D = pt.get<bool>("use2D", df_use2D);
        if let Some(field) = object.get("use2D") {
            self.use_2d = json_value_as_bool("use2D", field)?;
        }
        // RDKit✔️✔️:   d_minDistance = pt.get<unsigned int>("minDistance", d_minDistance);
        if let Some(field) = object.get("minDistance") {
            self.min_distance = json_value_as_u32("minDistance", field)?;
        }
        // RDKit✔️✔️:   d_maxDistance = pt.get<unsigned int>("maxDistance", d_maxDistance);
        if let Some(field) = object.get("maxDistance") {
            self.max_distance = json_value_as_u32("maxDistance", field)?;
        }
        // RDKit✔️✔️:   FingerprintArguments::fromJSON(pt);
        self.fingerprint_arguments.from_json_value(value)?;
        // RDKit✔️✔️: }
        Ok(())
    }
}

/// Parameters for the project-native AtomPair fingerprint API.
///
/// The defaults are the pinned RDKit AtomPair generator defaults. Per-call
/// selectors and custom atom invariants live here for the convenient molecule
/// method; reusable generators accept [`FingerprintFuncArguments`] directly.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomPairFingerprintParams {
    pub n_bits: usize,
    pub min_distance: u32,
    pub max_distance: u32,
    pub use_2d: bool,
    pub use_chirality: bool,
    pub count_simulation: bool,
    pub count_bounds: Vec<u32>,
    pub num_bits_per_feature: u32,
    pub from_atoms: Option<Vec<usize>>,
    pub ignore_atoms: Option<Vec<usize>>,
    pub conformer_id: i32,
    pub custom_atom_invariants: Option<Vec<u32>>,
    pub collect_additional_output: bool,
}

impl Default for AtomPairFingerprintParams {
    fn default() -> Self {
        let arguments = AtomPairArguments::default();
        Self {
            n_bits: arguments.fingerprint_arguments.d_fp_size as usize,
            min_distance: arguments.min_distance,
            max_distance: arguments.max_distance,
            use_2d: arguments.use_2d,
            use_chirality: arguments.fingerprint_arguments.df_include_chirality,
            count_simulation: arguments.fingerprint_arguments.df_count_simulation,
            count_bounds: arguments.fingerprint_arguments.d_count_bounds,
            num_bits_per_feature: arguments.fingerprint_arguments.d_num_bits_per_feature,
            from_atoms: None,
            ignore_atoms: None,
            conformer_id: -1,
            custom_atom_invariants: None,
            collect_additional_output: false,
        }
    }
}

/// Explicit-bit AtomPair fingerprint together with optional provenance.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomPairFingerprintOutput {
    pub fingerprint: Fingerprint,
    pub additional_output: Option<AdditionalOutput>,
}

#[derive(Debug, PartialEq, Eq)]
pub struct AtomPairAtomInvariantsGenerator {
    pub include_chirality: bool,
    pub topological_torsion_correction: bool,
}

impl AtomPairAtomInvariantsGenerator {
    #[must_use]
    pub const fn new(include_chirality: bool, topological_torsion_correction: bool) -> Self {
        // RDKit source file: AtomPairGenerator.cpp
        // RDKit source: AtomPairGenerator.cpp lines 26-29
        // RDKit✔️✔️: AtomPairAtomInvGenerator::AtomPairAtomInvGenerator(
        // RDKit✔️✔️:     bool includeChirality, bool topologicalTorsionCorrection)
        // RDKit✔️✔️:     : df_includeChirality(includeChirality),
        // RDKit✔️✔️:       df_topologicalTorsionCorrection(topologicalTorsionCorrection) {}
        Self {
            include_chirality,
            topological_torsion_correction,
        }
    }

    pub(crate) fn atom_invariants(&self, molecule: &Molecule) -> Result<Vec<u32>, AtomPairCodeError> {
        // RDKit source: AtomPairGenerator.cpp lines 31-43
        // RDKit✔️✔️: std::vector<std::uint32_t> *AtomPairAtomInvGenerator::getAtomInvariants(
        // RDKit✔️✔️:     const ROMol &mol) const {
        // RDKit✔️✔️:   auto *atomInvariants = new std::vector<std::uint32_t>(mol.getNumAtoms());
        let mut atom_invariants = vec![0; molecule.num_atoms()];

        // RDKit✔️✔️:   for (ROMol::ConstAtomIterator atomItI = mol.beginAtoms();
        // RDKit✔️✔️:        atomItI != mol.endAtoms(); ++atomItI) {
        for atom in molecule.atoms() {
            // RDKit✔️✔️:     (*atomInvariants)[(*atomItI)->getIdx()] =
            // RDKit✔️✔️:         getAtomCode(*atomItI, 0, df_includeChirality) -
            // RDKit✔️✔️:         (df_topologicalTorsionCorrection ? 2 : 0);
            let correction = if self.topological_torsion_correction { 2 } else { 0 };
            atom_invariants[atom.id().index()] =
                get_atom_code(molecule, atom.id(), 0, self.include_chirality)?.wrapping_sub(correction);
            // RDKit✔️✔️:   }
        }

        // RDKit✔️✔️:   return atomInvariants;
        Ok(atom_invariants)
        // RDKit✔️✔️: }
    }

    #[allow(non_snake_case)]
    pub fn getAtomInvariants(&self, molecule: &Molecule) -> Result<Vec<u32>, FingerprintError> {
        self.atom_invariants(molecule).map_err(FingerprintError::from)
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

    #[allow(non_snake_case)]
    #[must_use]
    pub fn infoString(&self) -> String {
        self.info_string()
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

    #[allow(non_snake_case)]
    #[must_use]
    pub fn toJSON(&self) -> String {
        self.to_json()
    }

    pub fn from_json(&mut self, json: &str) -> Result<(), FingerprintError> {
        if json.trim().is_empty() {
            return Ok(());
        }
        let value: Value =
            serde_json::from_str(json).map_err(|error| FingerprintError::InvalidArgumentsJson(error.to_string()))?;
        let object = value
            .as_object()
            .ok_or_else(|| FingerprintError::InvalidArgumentsJson("expected JSON object".to_string()))?;
        // RDKit source: AtomPairGenerator.cpp lines 56-61
        // RDKit✔️✔️: void AtomPairAtomInvGenerator::fromJSON(const boost::property_tree::ptree &pt) {
        // RDKit✔️✔️:   df_includeChirality = pt.get<bool>("includeChirality", df_includeChirality);
        if let Some(field) = object.get("includeChirality") {
            self.include_chirality = json_value_as_bool("includeChirality", field)?;
        }
        // RDKit✔️✔️:   df_topologicalTorsionCorrection = pt.get<bool>(
        // RDKit✔️✔️:       "topologicalTorsionCorrection", df_topologicalTorsionCorrection);
        if let Some(field) = object.get("topologicalTorsionCorrection") {
            self.topological_torsion_correction = json_value_as_bool("topologicalTorsionCorrection", field)?;
        }
        // RDKit✔️✔️:   AtomInvariantsGenerator::fromJSON(pt);
        // RDKit✔️✔️: }
        Ok(())
    }

    #[allow(non_snake_case)]
    pub fn fromJSON(&mut self, json: &str) -> Result<(), FingerprintError> {
        self.from_json(json)
    }
}

impl Default for AtomPairAtomInvariantsGenerator {
    fn default() -> Self {
        Self::new(false, false)
    }
}

impl Clone for AtomPairAtomInvariantsGenerator {
    fn clone(&self) -> Self {
        // RDKit source: AtomPairGenerator.cpp lines 63-66
        // RDKit✔️✔️: AtomPairAtomInvGenerator *AtomPairAtomInvGenerator::clone() const {
        // RDKit✔️✔️:   return new AtomPairAtomInvGenerator(df_includeChirality,
        // RDKit✔️✔️:                                       df_topologicalTorsionCorrection);
        // RDKit✔️✔️: }
        Self::new(self.include_chirality, self.topological_torsion_correction)
    }
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
        // RDKit✔️✔️: AtomPairAtomEnv<OutputType>::AtomPairAtomEnv(const unsigned int atomIdFirst,
        // RDKit✔️✔️:                                              const unsigned int atomIdSecond,
        // RDKit✔️✔️:                                              const unsigned int distance)
        // RDKit✔️✔️:     : d_atomIdFirst(atomIdFirst),
        // RDKit✔️✔️:       d_atomIdSecond(atomIdSecond),
        // RDKit✔️✔️:       d_distance(distance) {}
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
    ) -> Result<u32, AtomPairCodeError> {
        // RDKit source file: AtomPairGenerator.cpp
        // RDKit source: AtomPairGenerator.cpp lines 131-167
        // RDKit✔️✔️: OutputType AtomPairAtomEnv<OutputType>::getBitId(
        // RDKit✔️✔️:     FingerprintArguments *arguments,
        // RDKit✔️✔️:     const std::vector<std::uint32_t> *atomInvariants,
        // RDKit✔️✔️:     const std::vector<std::uint32_t> *,  // bondInvariants
        // RDKit✔️✔️:     AdditionalOutput *,                  // additionalOutput,
        // RDKit✔️✔️:     const bool hashResults,
        // RDKit✔️✔️:     const std::uint64_t  // fpSize
        // RDKit✔️✔️: ) const {
        // RDKit✔️✔️:   PRECONDITION((atomInvariants->size() >= d_atomIdFirst) &&
        // RDKit✔️✔️:                    (atomInvariants->size() >= d_atomIdSecond),
        // RDKit✔️✔️:                "bad atom invariants size");
        // The source spells this precondition with `>=`, while the indexed
        // access below requires strict containment. Preserve the effective
        // source boundary without reproducing an out-of-bounds access.
        for atom_index in [self.atom_id_first, self.atom_id_second] {
            if atom_index >= atom_invariants.len() {
                return Err(AtomPairCodeError::AtomInvariantLength {
                    length: atom_invariants.len(),
                    atom_index,
                });
            }
        }

        // RDKit✔️✔️:   auto *atomPairArguments = dynamic_cast<AtomPairArguments *>(arguments);
        // The concrete Rust type makes the source dynamic-cast precondition structural.

        // RDKit✔️✔️:   std::uint32_t codeSizeLimit =
        // RDKit✔️✔️:       (1 << (codeSize +
        // RDKit✔️✔️:              (atomPairArguments->df_includeChirality ? numChiralBits : 0))) -
        // RDKit✔️✔️:       1;
        let code_size_limit = (1u32
            << (CODE_SIZE
                + if arguments.df_include_chirality {
                    NUM_CHIRAL_BITS
                } else {
                    0
                }))
            - 1;

        // RDKit✔️✔️:   std::uint32_t atomCodeFirst =
        // RDKit✔️✔️:       (*atomInvariants)[d_atomIdFirst] % codeSizeLimit;
        let atom_code_first = atom_invariants[self.atom_id_first] % code_size_limit;

        // RDKit✔️✔️:   std::uint32_t atomCodeSecond =
        // RDKit✔️✔️:       (*atomInvariants)[d_atomIdSecond] % codeSizeLimit;
        let atom_code_second = atom_invariants[self.atom_id_second] % code_size_limit;

        // RDKit✔️✔️:   std::uint32_t bitId = 0;
        let mut bit_id = 0u32;
        // RDKit✔️✔️:   if (hashResults) {
        if hash_results {
            // RDKit✔️✔️:     gboost::hash_combine(bitId, std::min(atomCodeFirst, atomCodeSecond));
            hash_combine(&mut bit_id, atom_code_first.min(atom_code_second));
            // RDKit✔️✔️:     gboost::hash_combine(bitId, d_distance);
            hash_combine(&mut bit_id, self.distance);
            // RDKit✔️✔️:     gboost::hash_combine(bitId, std::max(atomCodeFirst, atomCodeSecond));
            hash_combine(&mut bit_id, atom_code_first.max(atom_code_second));
        // RDKit✔️✔️:   } else {
        } else {
            // RDKit✔️✔️:     bitId = getAtomPairCode(atomCodeFirst, atomCodeSecond, d_distance,
            // RDKit✔️✔️:                             atomPairArguments->df_includeChirality);
            bit_id = get_atom_pair_code(
                atom_code_first,
                atom_code_second,
                self.distance,
                arguments.df_include_chirality,
            )?;
            // RDKit✔️✔️:   }
        }

        // RDKit✔️✔️:   return bitId;
        Ok(bit_id)
        // RDKit✔️✔️: }
    }

    pub(crate) fn update_additional_output(&self, additional_output: &mut AdditionalOutput, bit_id: u64) {
        // RDKit source file: AtomPairGenerator.cpp
        // RDKit source: AtomPairGenerator.cpp lines 109-128
        // RDKit✔️✔️: void AtomPairAtomEnv<OutputType>::updateAdditionalOutput(
        // RDKit✔️✔️:     AdditionalOutput *additionalOutput, size_t bitId) const {
        // RDKit✔️✔️:   PRECONDITION(additionalOutput, "bad output pointer");
        // The Rust reference makes the non-null precondition structural.
        // RDKit✔️✔️:   if (additionalOutput->bitInfoMap) {
        if let Some(bit_info_map) = additional_output.bit_info_map.as_mut() {
            // RDKit✔️✔️:     (*additionalOutput->bitInfoMap)[bitId].emplace_back(d_atomIdFirst,
            // RDKit✔️✔️:                                                         d_atomIdSecond);
            bit_info_map
                .entry(bit_id)
                .or_default()
                .push((self.atom_id_first as u32, self.atom_id_second as u32));
            // RDKit✔️✔️:   }
        }
        // RDKit✔️✔️:   if (additionalOutput->atomToBits) {
        if let Some(atom_to_bits) = additional_output.atom_to_bits.as_mut() {
            // RDKit✔️✔️:     additionalOutput->atomToBits->at(d_atomIdFirst).push_back(bitId);
            atom_to_bits[self.atom_id_first].push(bit_id);
            // RDKit✔️✔️:     additionalOutput->atomToBits->at(d_atomIdSecond).push_back(bitId);
            atom_to_bits[self.atom_id_second].push(bit_id);
            // RDKit✔️✔️:   }
        }
        // RDKit✔️✔️:   if (additionalOutput->atomCounts) {
        if let Some(atom_counts) = additional_output.atom_counts.as_mut() {
            // RDKit✔️✔️:     additionalOutput->atomCounts->at(d_atomIdFirst)++;
            atom_counts[self.atom_id_first] += 1;
            // RDKit✔️✔️:     additionalOutput->atomCounts->at(d_atomIdSecond)++;
            atom_counts[self.atom_id_second] += 1;
            // RDKit✔️✔️:   }
        }
        // RDKit✔️✔️:   if (additionalOutput->atomsPerBit) {
        if let Some(atoms_per_bit) = additional_output.atoms_per_bit.as_mut() {
            // RDKit✔️✔️:     (*additionalOutput->atomsPerBit)[bitId].push_back(std::vector<int>{
            // RDKit✔️✔️:         static_cast<int>(d_atomIdFirst), static_cast<int>(d_atomIdSecond)});
            atoms_per_bit
                .entry(bit_id)
                .or_default()
                .push(vec![self.atom_id_first, self.atom_id_second]);
            // RDKit✔️✔️:   }
        }
        // RDKit✔️✔️: }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomPairFingerprintGenerator {
    atom_environment_generator: AtomPairEnvironmentGenerator,
    arguments: AtomPairArguments,
    atom_invariants_generator: AtomPairAtomInvariantsGenerator,
    owns_atom_invariants_generator: bool,
}

impl generator::FingerprintEnvironment for AtomPairEnvironment {
    fn bit_id(
        &self,
        arguments: &FingerprintArguments,
        atom_invariants: &[u32],
        _bond_invariants: &[u32],
        hash_results: bool,
        _fp_size: u64,
    ) -> Result<u64, FingerprintError> {
        AtomPairEnvironment::bit_id(self, arguments, atom_invariants, hash_results)
            .map(u64::from)
            .map_err(FingerprintError::from)
    }

    fn update_additional_output(&self, output: &mut AdditionalOutput, bit_id: u64) {
        AtomPairEnvironment::update_additional_output(self, output, bit_id);
    }
}

impl generator::FingerprintFamily for AtomPairFingerprintGenerator {
    type Environment = AtomPairEnvironment;

    fn common_arguments(&self) -> &FingerprintArguments {
        &self.arguments.fingerprint_arguments
    }

    fn result_size(&self) -> Result<u64, FingerprintError> {
        Ok(self.atom_environment_generator.result_size(&self.arguments))
    }

    fn arguments_info_string(&self) -> String {
        self.arguments.info_string()
    }

    fn environment_info_string(&self) -> String {
        self.atom_environment_generator.info_string().to_string()
    }

    fn atom_invariants_info_string(&self) -> Option<String> {
        Some(self.atom_invariants_generator.info_string())
    }

    fn bond_invariants_info_string(&self) -> Option<String> {
        None
    }

    fn arguments_json(&self) -> String {
        self.arguments.to_json()
    }

    fn environment_json(&self) -> String {
        self.atom_environment_generator.to_json().to_string()
    }

    fn atom_invariants_json(&self) -> Option<String> {
        Some(self.atom_invariants_generator.to_json())
    }

    fn bond_invariants_json(&self) -> Option<String> {
        None
    }

    fn atom_invariants(&self, molecule: &Molecule) -> Result<Vec<u32>, FingerprintError> {
        self.atom_invariants_generator
            .atom_invariants(molecule)
            .map_err(FingerprintError::from)
    }

    fn bond_invariants(&self, _molecule: &Molecule) -> Result<Vec<u32>, FingerprintError> {
        Ok(Vec::new())
    }

    fn environments(
        &self,
        molecule: &Molecule,
        from_atoms: Option<&[usize]>,
        ignore_atoms: Option<&[usize]>,
        conf_id: i32,
        _atom_invariants: &[u32],
        _bond_invariants: &[u32],
        _hash_results: bool,
    ) -> Result<Vec<Self::Environment>, FingerprintError> {
        self.atom_environment_generator
            .environments(molecule, &self.arguments, from_atoms, ignore_atoms, conf_id)
            .map_err(FingerprintError::from)
    }
}

impl AtomPairFingerprintGenerator {
    /// Construct a reusable AtomPair generator from project-native parameters.
    pub fn new(params: &AtomPairFingerprintParams) -> Result<Self, FingerprintError> {
        if params.n_bits == 0 {
            return Err(FingerprintError::EmptyFingerprint);
        }
        let fp_size = u32::try_from(params.n_bits).map_err(|_| FingerprintError::InvalidArguments {
            reason: "AtomPair n_bits exceeds the source uint32 range",
        })?;
        let mut generator = atom_pair_generator_with_parameters(
            params.min_distance,
            params.max_distance,
            params.use_chirality,
            params.use_2d,
            None,
            params.count_simulation,
            fp_size,
            params.count_bounds.clone(),
            false,
        )?;
        if params.num_bits_per_feature == 0 {
            return Err(FingerprintError::InvalidArguments {
                reason: "num_bits_per_feature must be greater than zero",
            });
        }
        generator.arguments.fingerprint_arguments.d_num_bits_per_feature = params.num_bits_per_feature;
        Ok(generator)
    }

    /// Restore an AtomPair generator through the common fingerprint JSON
    /// dispatcher. JSON for another fingerprint family is rejected.
    pub fn from_json(json: &str) -> Result<Self, FingerprintError> {
        match generator::generator_from_json(json)? {
            generator::RestoredFingerprintGenerator::AtomPair(generator) => Ok(generator),
            generator::RestoredFingerprintGenerator::Morgan(_) => Err(FingerprintError::InvalidArgumentsJson(
                "serialized generator is not an AtomPair generator".to_string(),
            )),
            generator::RestoredFingerprintGenerator::TopologicalTorsion(_) => Err(
                FingerprintError::InvalidArgumentsJson("serialized generator is not an AtomPair generator".to_string()),
            ),
        }
    }

    #[must_use]
    pub(crate) fn arguments(&self) -> &AtomPairArguments {
        &self.arguments
    }

    #[must_use]
    pub(crate) const fn owns_atom_invariants_generator(&self) -> bool {
        self.owns_atom_invariants_generator
    }

    #[must_use]
    pub fn info_string(&self) -> String {
        generator::FingerprintGenerator::new(self).info_string()
    }

    #[must_use]
    pub fn to_json(&self) -> String {
        generator::FingerprintGenerator::new(self).to_json()
    }

    pub fn sparse_count_fingerprint(
        &self,
        molecule: &Molecule,
        arguments: &mut FingerprintFuncArguments,
    ) -> Result<SparseCountFingerprint, FingerprintError> {
        generator::FingerprintGenerator::new(self).get_sparse_count_fingerprint(molecule, arguments)
    }

    pub fn sparse_bit_fingerprint(
        &self,
        molecule: &Molecule,
        arguments: &mut FingerprintFuncArguments,
    ) -> Result<SparseBitFingerprint, FingerprintError> {
        generator::FingerprintGenerator::new(self).get_sparse_fingerprint(molecule, arguments)
    }

    pub fn count_fingerprint(
        &self,
        molecule: &Molecule,
        arguments: &mut FingerprintFuncArguments,
    ) -> Result<SparseCountFingerprint, FingerprintError> {
        generator::FingerprintGenerator::new(self).get_count_fingerprint(molecule, arguments)
    }

    pub fn fingerprint(
        &self,
        molecule: &Molecule,
        arguments: &mut FingerprintFuncArguments,
    ) -> Result<Fingerprint, FingerprintError> {
        generator::FingerprintGenerator::new(self).get_fingerprint(molecule, arguments)
    }
}

/// Generate the explicit-bit AtomPair result without mutating the molecule.
pub fn atom_pair_fingerprint(
    molecule: &Molecule,
    params: &AtomPairFingerprintParams,
) -> Result<Fingerprint, FingerprintError> {
    Ok(atom_pair_fingerprint_with_output(molecule, params)?.fingerprint)
}

/// Generate an explicit-bit AtomPair fingerprint and, when requested, all
/// source-supported provenance containers.
pub fn atom_pair_fingerprint_with_output(
    molecule: &Molecule,
    params: &AtomPairFingerprintParams,
) -> Result<AtomPairFingerprintOutput, FingerprintError> {
    let generator = AtomPairFingerprintGenerator::new(params)?;
    let mut arguments = atom_pair_function_arguments(params);
    let fingerprint = generator.fingerprint(molecule, &mut arguments)?;
    Ok(AtomPairFingerprintOutput {
        fingerprint,
        additional_output: arguments.additional_output,
    })
}

pub(crate) fn atom_pair_function_arguments(params: &AtomPairFingerprintParams) -> FingerprintFuncArguments {
    let mut arguments = FingerprintFuncArguments {
        from_atoms: params.from_atoms.clone(),
        ignore_atoms: params.ignore_atoms.clone(),
        conf_id: params.conformer_id,
        custom_atom_invariants: params.custom_atom_invariants.clone(),
        ..Default::default()
    };
    if params.collect_additional_output {
        let mut additional_output = AdditionalOutput::new();
        additional_output.allocate_atom_counts();
        additional_output.allocate_atom_to_bits();
        additional_output.allocate_bit_info_map();
        additional_output.allocate_atoms_per_bit();
        arguments.additional_output = Some(additional_output);
    }
    arguments
}

/// Arguments retained only for source-compatibility verification of RDKit's
/// deprecated `AtomPairs.cpp` entry points. These adapters are intentionally
/// private and do not create a second public fingerprint implementation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct LegacyAtomPairArguments {
    pub(crate) min_length: u32,
    pub(crate) max_length: u32,
    pub(crate) from_atoms: Option<Vec<usize>>,
    pub(crate) ignore_atoms: Option<Vec<usize>>,
    pub(crate) atom_invariants: Option<Vec<u32>>,
    pub(crate) include_chirality: bool,
    pub(crate) use_2d: bool,
    pub(crate) conf_id: i32,
}

impl Default for LegacyAtomPairArguments {
    fn default() -> Self {
        Self {
            min_length: 1,
            max_length: MAX_PATH_LEN - 1,
            from_atoms: None,
            ignore_atoms: None,
            atom_invariants: None,
            include_chirality: false,
            use_2d: true,
            conf_id: -1,
        }
    }
}

fn validate_legacy_atom_invariants(
    molecule: &Molecule,
    arguments: &LegacyAtomPairArguments,
) -> Result<(), FingerprintError> {
    // RDKit source: AtomPairs.cpp lines 61-62 and 121-122.
    // RDKit✔️✔️: PRECONDITION(!atomInvariants ||
    // RDKit✔️✔️:              atomInvariants->size() >= mol.getNumAtoms(),
    // RDKit✔️✔️:              "bad atomInvariants size");
    if arguments
        .atom_invariants
        .as_ref()
        .is_some_and(|invariants| invariants.len() < molecule.num_atoms())
    {
        return Err(FingerprintError::InvalidArguments {
            reason: "legacy AtomPair atom invariants are shorter than the atom count",
        });
    }
    Ok(())
}

fn legacy_call_arguments(arguments: &LegacyAtomPairArguments) -> FingerprintFuncArguments {
    // RDKit source: AtomPairs.cpp lines 72-77.
    // RDKit✔️✔️: FingerprintFuncArguments args;
    // RDKit✔️✔️: args.fromAtoms = fromAtoms;
    // RDKit✔️✔️: args.ignoreAtoms = ignoreAtoms;
    // RDKit✔️✔️: args.customAtomInvariants = atomInvariants;
    // RDKit✔️✔️: args.confId = confId;
    FingerprintFuncArguments {
        from_atoms: arguments.from_atoms.clone(),
        ignore_atoms: arguments.ignore_atoms.clone(),
        conf_id: arguments.conf_id,
        custom_atom_invariants: arguments.atom_invariants.clone(),
        ..Default::default()
    }
}

fn legacy_generator(
    n_bits: u32,
    arguments: &LegacyAtomPairArguments,
    count_bounds: Vec<u32>,
) -> Result<AtomPairFingerprintGenerator, FingerprintError> {
    // RDKit source: AtomPairs.cpp lines 78-81.
    // RDKit✔️✔️: RDKit::AtomPair::getAtomPairGenerator<std::uint32_t>(
    // RDKit✔️✔️:     minLength, maxLength, includeChirality, use2D, nullptr,
    // RDKit✔️✔️:     true, nBits)
    atom_pair_generator_with_parameters(
        arguments.min_length,
        arguments.max_length,
        arguments.include_chirality,
        arguments.use_2d,
        None,
        true,
        n_bits,
        count_bounds,
        false,
    )
}

pub(crate) fn legacy_sparse_count_fingerprint(
    molecule: &Molecule,
    arguments: &LegacyAtomPairArguments,
) -> Result<SparseCountFingerprint, FingerprintError> {
    validate_legacy_atom_invariants(molecule, arguments)?;
    // RDKit source: AtomPairs.cpp lines 89-109.
    // RDKit✔️✔️: getAtomPairFingerprintInternal(mol, 0, minLength,
    // RDKit✔️✔️:     maxLength, fromAtoms, ignoreAtoms, atomInvariants,
    // RDKit✔️✔️:     includeChirality, use2D, confId, true)
    let generator = legacy_generator(0, arguments, vec![1, 2, 4, 8])?;
    generator.sparse_count_fingerprint(molecule, &mut legacy_call_arguments(arguments))
}

pub(crate) fn legacy_hashed_count_fingerprint(
    molecule: &Molecule,
    n_bits: u32,
    arguments: &LegacyAtomPairArguments,
) -> Result<SparseCountFingerprint, FingerprintError> {
    validate_legacy_atom_invariants(molecule, arguments)?;
    // RDKit source: AtomPairs.cpp lines 111-118.
    // RDKit✔️✔️: getAtomPairFingerprintInternal(mol, nBits, minLength,
    // RDKit✔️✔️:     maxLength, fromAtoms, ignoreAtoms, atomInvariants,
    // RDKit✔️✔️:     includeChirality, use2D, confId, false)
    let generator = legacy_generator(n_bits, arguments, vec![1, 2, 4, 8])?;
    generator.count_fingerprint(molecule, &mut legacy_call_arguments(arguments))
}

pub(crate) fn legacy_hashed_bit_fingerprint(
    molecule: &Molecule,
    n_bits: u32,
    n_bits_per_entry: u32,
    arguments: &LegacyAtomPairArguments,
) -> Result<Fingerprint, FingerprintError> {
    validate_legacy_atom_invariants(molecule, arguments)?;
    if n_bits_per_entry == 0 {
        return Err(FingerprintError::InvalidArguments {
            reason: "legacy AtomPair n_bits_per_entry must be greater than zero",
        });
    }

    // RDKit source: AtomPairs.cpp lines 124-151.
    // RDKit✔️✔️: static int bounds[4] = {1, 2, 4, 8};
    // RDKit✔️✔️: unsigned int blockLength = nBits / nBitsPerEntry;
    // RDKit✔️✔️: if (nBitsPerEntry != 4) count thresholds are i + 1;
    // RDKit✔️✔️: else the thresholds are bounds[i].
    // The sole modern generator core performs the same block-length division
    // and bit projection when given these source-derived count bounds.
    let count_bounds = if n_bits_per_entry == 4 {
        vec![1, 2, 4, 8]
    } else {
        (1..=n_bits_per_entry).collect()
    };
    let generator = legacy_generator(n_bits, arguments, count_bounds)?;
    generator.fingerprint(molecule, &mut legacy_call_arguments(arguments))
}

pub(crate) fn atom_pair_generator(
    arguments: &AtomPairArguments,
    atom_invariants_generator: Option<AtomPairAtomInvariantsGenerator>,
    owns_atom_invariants_generator: bool,
) -> AtomPairFingerprintGenerator {
    // RDKit source file: AtomPairGenerator.cpp
    // RDKit source: AtomPairGenerator.cpp lines 252-268
    // RDKit✔️✔️: FingerprintGenerator<OutputType> *getAtomPairGenerator(
    // RDKit✔️✔️:     const AtomPairArguments &args,
    // RDKit✔️✔️:     AtomInvariantsGenerator *atomInvariantsGenerator,
    // RDKit✔️✔️:     const bool ownsAtomInvGen) {
    // RDKit✔️✔️:   AtomEnvironmentGenerator<OutputType> *atomPairEnvGenerator =
    // RDKit✔️✔️:       new AtomPair::AtomPairEnvGenerator<OutputType>();
    let atom_environment_generator = AtomPairEnvironmentGenerator;
    // RDKit✔️✔️:   bool ownsAtomInvGenerator = ownsAtomInvGen;
    let mut owns_atom_invariants_generator = owns_atom_invariants_generator;
    // RDKit✔️✔️:   if (!atomInvariantsGenerator) {
    // RDKit✔️✔️:     atomInvariantsGenerator =
    // RDKit✔️✔️:         new AtomPairAtomInvGenerator(args.df_includeChirality);
    // RDKit✔️✔️:     ownsAtomInvGenerator = true;
    // RDKit✔️✔️:   }
    let atom_invariants_generator = match atom_invariants_generator {
        Some(generator) => generator,
        None => {
            owns_atom_invariants_generator = true;
            AtomPairAtomInvariantsGenerator::new(arguments.fingerprint_arguments.df_include_chirality, false)
        }
    };

    // RDKit✔️✔️:   return new FingerprintGenerator<OutputType>(
    // RDKit✔️✔️:       atomPairEnvGenerator, new AtomPairArguments(args),
    // RDKit✔️✔️:       atomInvariantsGenerator, nullptr, ownsAtomInvGenerator, false);
    // RDKit✔️✔️: }
    AtomPairFingerprintGenerator {
        atom_environment_generator,
        arguments: arguments.clone(),
        atom_invariants_generator,
        owns_atom_invariants_generator,
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn atom_pair_generator_with_parameters(
    min_distance: u32,
    max_distance: u32,
    include_chirality: bool,
    use_2d: bool,
    atom_invariants_generator: Option<AtomPairAtomInvariantsGenerator>,
    use_count_simulation: bool,
    fp_size: u32,
    count_bounds: Vec<u32>,
    owns_atom_invariants_generator: bool,
) -> Result<AtomPairFingerprintGenerator, FingerprintError> {
    // RDKit source: AtomPairGenerator.cpp lines 270-283
    // RDKit✔️✔️: FingerprintGenerator<OutputType> *getAtomPairGenerator(
    // RDKit✔️✔️:     const unsigned int minDistance, const unsigned int maxDistance,
    // RDKit✔️✔️:     const bool includeChirality, const bool use2D,
    // RDKit✔️✔️:     AtomInvariantsGenerator *atomInvariantsGenerator,
    // RDKit✔️✔️:     const bool useCountSimulation, const std::uint32_t fpSize,
    // RDKit✔️✔️:     const std::vector<std::uint32_t> countBounds, const bool ownsAtomInvGen) {
    // RDKit✔️✔️:   AtomPair::AtomPairArguments arguments(useCountSimulation, includeChirality,
    // RDKit✔️✔️:                                         use2D, minDistance, maxDistance,
    // RDKit✔️✔️:                                         countBounds, fpSize);
    let arguments = AtomPairArguments::new(
        use_count_simulation,
        include_chirality,
        use_2d,
        min_distance,
        max_distance,
        count_bounds,
        fp_size,
    )?;

    // RDKit✔️✔️:   return getAtomPairGenerator<OutputType>(arguments, atomInvariantsGenerator,
    // RDKit✔️✔️:                                           ownsAtomInvGen);
    // RDKit✔️✔️: }
    Ok(atom_pair_generator(
        &arguments,
        atom_invariants_generator,
        owns_atom_invariants_generator,
    ))
}

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub(crate) enum AtomPairCodeError {
    #[error("atom index {atom} is outside the molecule")]
    AtomIndex { atom: AtomId },
    #[error(transparent)]
    AtomProperty(#[from] AtomPropertyError),
    #[error("CIPLabeler failed while preparing AtomPair atom chirality: {reason}")]
    CipLabeler { reason: String },
    #[error("AtomPair atom code {code} exceeds its {width}-bit representation")]
    CodeWidth { code: u32, width: u32 },
    #[error("AtomPair distance {distance} must be less than {maximum}")]
    DistanceTooLong { distance: u32, maximum: u32 },
    #[error("AtomPair atom invariants length {length} does not contain endpoint index {atom_index}")]
    AtomInvariantLength { length: usize, atom_index: usize },
    #[error("AtomPair distance-matrix construction failed: {0:?}")]
    DistanceMatrix(DistanceMatrixError),
}

impl From<crate::chemistry::ciplabeler::CipLabelerError> for AtomPairCodeError {
    fn from(error: crate::chemistry::ciplabeler::CipLabelerError) -> Self {
        Self::CipLabeler {
            reason: error.to_string(),
        }
    }
}

impl From<AtomPairCodeError> for FingerprintError {
    fn from(error: AtomPairCodeError) -> Self {
        match error {
            AtomPairCodeError::AtomProperty(AtomPropertyError::MissingExplicitValence { atom }) => {
                Self::Valence(crate::ValenceError::ExplicitValenceCacheNotInitialized { atom })
            }
            AtomPairCodeError::AtomProperty(AtomPropertyError::Valence(error)) => Self::Valence(error),
            other => Self::AtomPair {
                reason: other.to_string(),
            },
        }
    }
}

pub(crate) fn get_atom_code(
    molecule: &Molecule,
    atom_id: AtomId,
    branch_subtract: u32,
    include_chirality: bool,
) -> Result<u32, AtomPairCodeError> {
    get_atom_code_with_stereo_mode(
        molecule,
        atom_id,
        branch_subtract,
        include_chirality,
        rdkit_use_legacy_stereo_perception(),
    )
}

fn get_atom_code_with_stereo_mode(
    molecule: &Molecule,
    atom_id: AtomId,
    branch_subtract: u32,
    include_chirality: bool,
    use_legacy_stereo_perception: bool,
) -> Result<u32, AtomPairCodeError> {
    // RDKit source file: FingerprintUtil.cpp
    // RDKit source: FingerprintUtil.cpp lines 45-97
    // RDKit✔️✔️: std::uint32_t getAtomCode(const Atom *atom, unsigned int branchSubtract,
    // RDKit✔️✔️:                           bool includeChirality) {
    // RDKit✔️✔️:   PRECONDITION(atom, "no atom");
    let atom = molecule
        .atoms()
        .get(atom_id.index())
        .ok_or(AtomPairCodeError::AtomIndex { atom: atom_id })?;
    // RDKit✔️✔️:   std::uint32_t code;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   unsigned int numBranches = 0;
    let degree = molecule.topology_block().adjacency.neighbors_of(atom_id.index()).len() as u32;
    // RDKit✔️✔️:   if (atom->getDegree() > branchSubtract) {
    // RDKit✔️✔️:     numBranches = atom->getDegree() - branchSubtract;
    // RDKit✔️✔️:   }
    let num_branches = if degree > branch_subtract {
        degree - branch_subtract
    } else {
        0
    };

    // RDKit✔️✔️:   code = numBranches % maxNumBranches;
    let mut code = num_branches % MAX_NUM_BRANCHES;
    // RDKit✔️✔️:   unsigned int nPi = numPiElectrons(*atom) % maxNumPi;
    let num_pi = num_pi_electrons(molecule, atom_id)? % MAX_NUM_PI;
    // RDKit✔️✔️:   code |= nPi << numBranchBits;
    code |= num_pi << NUM_BRANCH_BITS;

    // RDKit✔️✔️:   unsigned int typeIdx = 0;
    let mut type_index = 0usize;
    // RDKit✔️✔️:   unsigned int nTypes = 1 << numTypeBits;
    let num_types = 1usize << NUM_TYPE_BITS;
    // RDKit✔️✔️:   while (typeIdx < nTypes) {
    while type_index < num_types {
        // RDKit✔️✔️:     if (atomNumberTypes[typeIdx] ==
        // RDKit✔️✔️:         static_cast<unsigned int>(atom->getAtomicNum())) {
        if ATOM_NUMBER_TYPES[type_index] == u32::from(atom.atomic_number()) {
            // RDKit✔️✔️:       break;
            break;
        // RDKit✔️✔️:     } else if (atomNumberTypes[typeIdx] >
        // RDKit✔️✔️:                static_cast<unsigned int>(atom->getAtomicNum())) {
        } else if ATOM_NUMBER_TYPES[type_index] > u32::from(atom.atomic_number()) {
            // RDKit✔️✔️:       typeIdx = nTypes;
            type_index = num_types;
            // RDKit✔️✔️:       break;
            break;
        }
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     ++typeIdx;
        type_index += 1;
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   if (typeIdx == nTypes) {
    if type_index == num_types {
        // RDKit✔️✔️:     --typeIdx;
        type_index -= 1;
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   code |= typeIdx << (numBranchBits + numPiBits);
    code |= (type_index as u32) << (NUM_BRANCH_BITS + NUM_PI_BITS);
    // RDKit✔️✔️:   if (includeChirality) {
    if include_chirality {
        // RDKit✔️✔️:     // if we aren't using legacy stereo, we need to compute the CIP codes
        // RDKit✔️✔️:     if (!Chirality::getUseLegacyStereoPerception() &&
        // RDKit✔️✔️:         atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
        // RDKit✔️✔️:         !atom->getOwningMol().hasProp(common_properties::_CIPComputed)) {
        let labeled;
        let chirality_atom = if !use_legacy_stereo_perception
            && atom.chiral_tag() != ChiralTag::Unspecified
            && molecule.prop("_CIPComputed").is_none()
        {
            // RDKit✔️✔️:       CIPLabeler::assignCIPLabels(atom->getOwningMol());
            labeled = assign_cip_labels(molecule, 0)?;
            &labeled.atoms()[atom_id.index()]
        // RDKit✔️✔️:     }
        } else {
            atom
        };
        // RDKit✔️✔️:     std::string cipCode;
        // RDKit✔️✔️:     if (atom->getPropIfPresent(common_properties::_CIPCode, cipCode)) {
        if let Some(cip_code) = chirality_atom.prop("_CIPCode") {
            // RDKit✔️✔️:       std::uint32_t offset = numBranchBits + numPiBits + numTypeBits;
            let offset = NUM_BRANCH_BITS + NUM_PI_BITS + NUM_TYPE_BITS;
            // RDKit✔️✔️:       if (cipCode == "R") {
            if cip_code == "R" {
                // RDKit✔️✔️:         code |= 1 << offset;
                code |= 1u32 << offset;
            // RDKit✔️✔️:       } else if (cipCode == "S") {
            } else if cip_code == "S" {
                // RDKit✔️✔️:         code |= 2 << offset;
                code |= 2u32 << offset;
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   POSTCONDITION(code < static_cast<std::uint32_t>(
    // RDKit✔️✔️:                            1 << (codeSize + (includeChirality ? 2 : 0))),
    // RDKit✔️✔️:                 "code exceeds number of bits");
    let width = CODE_SIZE + if include_chirality { NUM_CHIRAL_BITS } else { 0 };
    if code >= (1u32 << width) {
        return Err(AtomPairCodeError::CodeWidth { code, width });
    }
    // RDKit✔️✔️:   return code;
    Ok(code)
    // RDKit✔️✔️: };
}

pub(crate) fn get_atom_pair_code(
    code_i: u32,
    code_j: u32,
    distance: u32,
    include_chirality: bool,
) -> Result<u32, AtomPairCodeError> {
    // RDKit source file: FingerprintUtil.cpp
    // RDKit source: FingerprintUtil.cpp lines 99-107
    // RDKit✔️✔️: std::uint32_t getAtomPairCode(std::uint32_t codeI, std::uint32_t codeJ,
    // RDKit✔️✔️:                               unsigned int dist, bool includeChirality) {
    // RDKit✔️✔️:   PRECONDITION(dist < maxPathLen, "dist too long");
    if distance >= MAX_PATH_LEN {
        return Err(AtomPairCodeError::DistanceTooLong {
            distance,
            maximum: MAX_PATH_LEN,
        });
    }
    // RDKit✔️✔️:   std::uint32_t res = dist;
    let mut result = distance;
    // RDKit✔️✔️:   res |= std::min(codeI, codeJ) << numPathBits;
    result |= code_i.min(code_j) << NUM_PATH_BITS;
    // RDKit✔️✔️:   res |= std::max(codeI, codeJ)
    // RDKit✔️✔️:          << (numPathBits + codeSize + (includeChirality ? numChiralBits : 0));
    result |= code_i.max(code_j) << (NUM_PATH_BITS + CODE_SIZE + if include_chirality { NUM_CHIRAL_BITS } else { 0 });
    // RDKit✔️✔️:   return res;
    Ok(result)
    // RDKit✔️✔️: }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{AtomSpec, BondOrder, BondSpec, Conformer3D, Element, Hybridization};

    fn isolated_atom(atomic_number: u8) -> Molecule {
        let mut builder = Molecule::builder();
        builder.add_atom(
            AtomSpec::new(Element::from_atomic_number(atomic_number).unwrap()).with_hybridization(Hybridization::Sp3),
        );
        builder.build().unwrap()
    }

    #[test]
    fn atom_code_covers_every_source_element_bucket_and_fallback_boundaries() {
        for (type_index, &atomic_number) in ATOM_NUMBER_TYPES[..15].iter().enumerate() {
            let molecule = isolated_atom(atomic_number as u8);
            assert_eq!(
                get_atom_code(&molecule, AtomId::new(0), 0, false),
                Ok((type_index as u32) << (NUM_BRANCH_BITS + NUM_PI_BITS)),
                "atomic number {atomic_number}"
            );
        }

        for atomic_number in [1, 4, 10, 13, 18, 32, 36, 50, 54, 118] {
            let molecule = isolated_atom(atomic_number);
            assert_eq!(
                get_atom_code(&molecule, AtomId::new(0), 0, false),
                Ok(15 << (NUM_BRANCH_BITS + NUM_PI_BITS)),
                "fallback atomic number {atomic_number}"
            );
        }
    }

    #[test]
    fn atom_code_preserves_degree_subtraction_and_modulo_seven() {
        let mut builder = Molecule::builder();
        let center = builder.add_atom(AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3));
        for _ in 0..8 {
            let neighbor = builder.add_atom(AtomSpec::new(Element::H).with_hybridization(Hybridization::Sp3));
            builder
                .add_bond(BondSpec::new(center, neighbor, BondOrder::Single))
                .unwrap();
        }
        let molecule = builder.build().unwrap();

        assert_eq!(get_atom_code(&molecule, center, 0, false), Ok(33));
        assert_eq!(get_atom_code(&molecule, center, 1, false), Ok(32));
        assert_eq!(get_atom_code(&molecule, center, 8, false), Ok(32));
        assert_eq!(get_atom_code(&molecule, center, 9, false), Ok(32));
    }

    #[test]
    fn atom_code_preserves_aromatic_multiple_bond_and_explicit_hydrogen_pi_bits() {
        let aromatic = Molecule::from_smiles("c1ccccc1").unwrap();
        assert_eq!(get_atom_code(&aromatic, AtomId::new(0), 0, false), Ok(42));

        let carbon_dioxide = Molecule::from_smiles("O=C=O").unwrap();
        assert_eq!(get_atom_code(&carbon_dioxide, AtomId::new(1), 0, false), Ok(50));
        assert_eq!(get_atom_code(&carbon_dioxide, AtomId::new(0), 0, false), Ok(105));

        let explicit_hydrogens = Molecule::from_smiles("[CH2]=C").unwrap();
        assert_eq!(get_atom_code(&explicit_hydrogens, AtomId::new(0), 0, false), Ok(41));
    }

    #[test]
    fn atom_code_preserves_r_s_absent_and_disabled_cip_bits() {
        fn chiral_atom(cip: Option<&str>) -> Molecule {
            let mut spec = AtomSpec::new(Element::C)
                .with_hybridization(Hybridization::Sp3)
                .with_chiral_tag(ChiralTag::TetrahedralCw);
            if let Some(cip) = cip {
                spec = spec.with_prop("_CIPCode", cip);
            }
            let mut builder = Molecule::builder();
            builder.add_atom(spec);
            builder.build().unwrap()
        }

        let r = chiral_atom(Some("R"));
        let s = chiral_atom(Some("S"));
        let absent = chiral_atom(None);
        assert_eq!(get_atom_code(&r, AtomId::new(0), 0, true), Ok(544));
        assert_eq!(get_atom_code(&s, AtomId::new(0), 0, true), Ok(1056));
        assert_eq!(get_atom_code(&absent, AtomId::new(0), 0, true), Ok(32));
        assert_eq!(get_atom_code(&r, AtomId::new(0), 0, false), Ok(32));
    }

    #[test]
    fn atom_code_new_stereo_mode_computes_missing_cip_labels_without_mutating_input() {
        let molecule = Molecule::from_smiles("C[C@H](O)F").unwrap();
        let center = AtomId::new(1);
        let original_cip = molecule.atoms()[center.index()].prop("_CIPCode");
        let code = get_atom_code_with_stereo_mode(&molecule, center, 0, true, false).unwrap();
        assert!(matches!(code >> CODE_SIZE, 1 | 2));
        assert_eq!(molecule.atoms()[center.index()].prop("_CIPCode"), original_cip);
    }

    #[test]
    fn atom_code_reports_invalid_atom_and_missing_valence_state() {
        let molecule = isolated_atom(6);
        assert_eq!(
            get_atom_code(&molecule, AtomId::new(1), 0, false),
            Err(AtomPairCodeError::AtomIndex { atom: AtomId::new(1) })
        );

        let mut builder = Molecule::builder();
        let atom = builder.add_atom(AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2));
        let uncached = builder.build().unwrap();
        assert_eq!(
            get_atom_code(&uncached, atom, 0, false),
            Err(AtomPairCodeError::AtomProperty(
                AtomPropertyError::MissingExplicitValence { atom }
            ))
        );
    }

    #[test]
    fn pair_code_canonicalizes_reversed_endpoints_in_both_width_modes() {
        assert_eq!(get_atom_pair_code(10, 20, 1, false), Ok(328_001));
        assert_eq!(get_atom_pair_code(20, 10, 1, false), Ok(328_001));
        assert_eq!(get_atom_pair_code(10, 20, 1, true), Ok(1_311_041));
        assert_eq!(get_atom_pair_code(20, 10, 1, true), Ok(1_311_041));
    }

    #[test]
    fn pair_code_preserves_distance_zero_and_last_valid_distance() {
        assert_eq!(get_atom_pair_code(0, 0, 0, false), Ok(0));
        assert_eq!(get_atom_pair_code(0, 0, 30, false), Ok(30));
        assert_eq!(get_atom_pair_code(1, 2, 0, false), Ok(32_800));
        assert_eq!(get_atom_pair_code(1, 2, 30, false), Ok(32_830));
    }

    #[test]
    fn pair_code_preserves_maximum_valid_codes_and_u32_packing() {
        let max_unchiral_code = (1u32 << CODE_SIZE) - 1;
        assert_eq!(
            get_atom_pair_code(max_unchiral_code, max_unchiral_code, 30, false),
            Ok(8_388_606)
        );

        let max_chiral_code = (1u32 << (CODE_SIZE + NUM_CHIRAL_BITS)) - 1;
        assert_eq!(
            get_atom_pair_code(max_chiral_code, max_chiral_code, 30, true),
            Ok(134_217_726)
        );

        assert_eq!(get_atom_pair_code(u32::MAX, u32::MAX, 0, false), Ok(u32::MAX - 31));
    }

    #[test]
    fn pair_code_rejects_source_distance_precondition_boundaries() {
        assert_eq!(
            get_atom_pair_code(0, 0, 31, false),
            Err(AtomPairCodeError::DistanceTooLong {
                distance: 31,
                maximum: 31,
            })
        );
        assert_eq!(
            get_atom_pair_code(0, 0, u32::MAX, true),
            Err(AtomPairCodeError::DistanceTooLong {
                distance: u32::MAX,
                maximum: 31,
            })
        );
    }

    #[test]
    fn arguments_defaults_and_metadata_match_the_source_constructor() {
        let arguments = AtomPairArguments::default();
        assert!(arguments.fingerprint_arguments.df_count_simulation);
        assert!(!arguments.fingerprint_arguments.df_include_chirality);
        assert_eq!(arguments.fingerprint_arguments.d_count_bounds, [1, 2, 4, 8]);
        assert_eq!(arguments.fingerprint_arguments.d_fp_size, 2048);
        assert_eq!(arguments.fingerprint_arguments.d_num_bits_per_feature, 1);
        assert!(arguments.use_2d);
        assert_eq!(arguments.min_distance, 1);
        assert_eq!(arguments.max_distance, 30);
        assert_eq!(
            arguments.info_string(),
            "AtomPairArguments use2D=1 minDistance=1 maxDistance=30"
        );
    }

    #[test]
    fn arguments_constructor_enforces_source_common_and_distance_preconditions() {
        assert!(matches!(
            AtomPairArguments::new(true, false, true, 1, 30, Vec::new(), 2048),
            Err(FingerprintError::InvalidArguments {
                reason: "countSimulation requires non-empty countBounds"
            })
        ));
        assert!(matches!(
            AtomPairArguments::new(false, false, true, 4, 3, Vec::new(), 2048),
            Err(FingerprintError::InvalidArguments {
                reason: "minDistance must be less than or equal to maxDistance"
            })
        ));
        let arguments = AtomPairArguments::new(false, true, false, 0, 0, Vec::new(), 17).unwrap();
        assert!(!arguments.fingerprint_arguments.df_count_simulation);
        assert!(arguments.fingerprint_arguments.df_include_chirality);
        assert!(!arguments.use_2d);
        assert_eq!(arguments.min_distance, 0);
        assert_eq!(arguments.max_distance, 0);
        assert_eq!(arguments.fingerprint_arguments.d_fp_size, 17);
    }

    #[test]
    fn arguments_json_has_exact_family_then_common_field_shape_and_round_trips() {
        let original = AtomPairArguments::default();
        assert_eq!(
            original.to_json(),
            r#"{"type":"AtomPairArguments","use2D":"true","minDistance":"1","maxDistance":"30","countSimulation":"true","fpSize":"2048","numBitsPerFeature":"1","includeChirality":"false","countBounds":["1","2","4","8"]}"#
        );

        let mut restored = AtomPairArguments::new(false, true, false, 0, 0, Vec::new(), 32).unwrap();
        restored.from_json(&original.to_json()).unwrap();
        assert_eq!(restored, original);
    }

    #[test]
    fn arguments_partial_json_preserves_unmentioned_family_fields_and_source_clears_bounds() {
        let mut arguments = AtomPairArguments::default();
        arguments.from_json(r#"{"minDistance":4,"fpSize":4096}"#).unwrap();
        assert!(arguments.use_2d);
        assert_eq!(arguments.min_distance, 4);
        assert_eq!(arguments.max_distance, 30);
        assert_eq!(arguments.fingerprint_arguments.d_fp_size, 4096);
        assert!(arguments.fingerprint_arguments.df_count_simulation);
        assert!(arguments.fingerprint_arguments.d_count_bounds.is_empty());

        let before = arguments.clone();
        arguments.from_json("").unwrap();
        assert_eq!(arguments, before);
    }

    #[test]
    fn arguments_json_rejects_invalid_document_and_field_types() {
        for invalid in [
            "[1,2,3]",
            r#"{"use2D":"yes"}"#,
            r#"{"minDistance":-1}"#,
            r#"{"maxDistance":1.5}"#,
            r#"{"fpSize":4294967296}"#,
            r#"{"countBounds":{}}"#,
            "{",
        ] {
            let mut arguments = AtomPairArguments::default();
            assert!(arguments.from_json(invalid).is_err(), "{invalid}");
        }
    }

    #[test]
    fn atom_invariant_generator_preserves_atom_index_placement_and_exact_codes() {
        let mut builder = Molecule::builder();
        for element in [Element::O, Element::C, Element::B] {
            builder.add_atom(AtomSpec::new(element).with_hybridization(Hybridization::Sp3));
        }
        let molecule = builder.build().unwrap();
        assert_eq!(
            AtomPairAtomInvariantsGenerator::default()
                .atom_invariants(&molecule)
                .unwrap(),
            [96, 32, 0]
        );
    }

    #[test]
    fn atom_invariant_generator_preserves_chiral_codes() {
        let mut builder = Molecule::builder();
        builder.add_atom(
            AtomSpec::new(Element::C)
                .with_hybridization(Hybridization::Sp3)
                .with_chiral_tag(ChiralTag::TetrahedralCw)
                .with_prop("_CIPCode", "R"),
        );
        builder.add_atom(
            AtomSpec::new(Element::C)
                .with_hybridization(Hybridization::Sp3)
                .with_chiral_tag(ChiralTag::TetrahedralCcw)
                .with_prop("_CIPCode", "S"),
        );
        let molecule = builder.build().unwrap();

        assert_eq!(
            AtomPairAtomInvariantsGenerator::new(true, false)
                .atom_invariants(&molecule)
                .unwrap(),
            [544, 1056]
        );
    }

    #[test]
    fn atom_invariant_generator_uses_source_unsigned_torsion_correction() {
        let boron = isolated_atom(5);
        assert_eq!(
            AtomPairAtomInvariantsGenerator::new(false, true)
                .atom_invariants(&boron)
                .unwrap(),
            [u32::MAX - 1]
        );

        let mut builder = Molecule::builder();
        let center = builder.add_atom(AtomSpec::new(Element::B).with_hybridization(Hybridization::Sp3));
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H).with_hybridization(Hybridization::Sp3));
        builder
            .add_bond(BondSpec::new(center, hydrogen, BondOrder::Single))
            .unwrap();
        let bonded = builder.build().unwrap();
        assert_eq!(
            AtomPairAtomInvariantsGenerator::new(false, true)
                .atom_invariants(&bonded)
                .unwrap()[center.index()],
            u32::MAX
        );
    }

    #[test]
    fn atom_invariant_generator_metadata_json_and_partial_restore_match_source() {
        let original = AtomPairAtomInvariantsGenerator::new(true, true);
        assert_eq!(
            original.info_string(),
            "AtomPairInvariantGenerator topologicalTorsionCorrection=1"
        );
        assert_eq!(
            original.to_json(),
            r#"{"type":"AtomPairAtomInvGenerator","includeChirality":"true","topologicalTorsionCorrection":"true"}"#
        );

        let mut restored = AtomPairAtomInvariantsGenerator::default();
        restored.from_json(&original.to_json()).unwrap();
        assert_eq!(restored, original);
        restored.from_json(r#"{"includeChirality":false}"#).unwrap();
        assert!(!restored.include_chirality);
        assert!(restored.topological_torsion_correction);
        restored.from_json(r#"{"includeChirality":1}"#).unwrap();
        assert!(restored.include_chirality);
        assert!(restored.from_json(r#"{"includeChirality":2}"#).is_err());
        assert!(restored.from_json("[]").is_err());
    }

    #[test]
    fn atom_invariant_generator_clone_has_independent_owned_state() {
        let original = AtomPairAtomInvariantsGenerator::new(true, false);
        let mut cloned = original.clone();
        cloned.include_chirality = false;
        cloned.topological_torsion_correction = true;
        assert_eq!(original, AtomPairAtomInvariantsGenerator::new(true, false));
        assert_eq!(cloned, AtomPairAtomInvariantsGenerator::new(false, true));
    }

    #[test]
    fn environment_bit_id_exact_mode_matches_pair_packing_and_endpoint_reversal() {
        let arguments = AtomPairArguments::default();
        let forward = AtomPairEnvironment::new(0, 1, 1);
        let reverse = AtomPairEnvironment::new(1, 0, 1);
        let invariants = [10, 20];
        assert_eq!(
            forward.bit_id(&arguments.fingerprint_arguments, &invariants, false),
            Ok(328_001)
        );
        assert_eq!(
            reverse.bit_id(&arguments.fingerprint_arguments, &invariants, false),
            Ok(328_001)
        );
    }

    #[test]
    fn environment_bit_id_hashed_mode_uses_canonical_three_stage_gboost_order() {
        let arguments = AtomPairArguments::default();
        let forward = AtomPairEnvironment::new(0, 1, 1);
        let reverse = AtomPairEnvironment::new(1, 0, 1);
        let invariants = [10, 20];
        assert_eq!(
            forward.bit_id(&arguments.fingerprint_arguments, &invariants, true),
            Ok(4_217_127_294)
        );
        assert_eq!(
            reverse.bit_id(&arguments.fingerprint_arguments, &invariants, true),
            Ok(4_217_127_294)
        );
    }

    #[test]
    fn environment_bit_id_uses_source_modulo_limit_not_a_bit_mask() {
        let unchiral = AtomPairArguments::default();
        let environment = AtomPairEnvironment::new(0, 1, 1);
        assert_eq!(
            environment.bit_id(&unchiral.fingerprint_arguments, &[510, 511], false),
            get_atom_pair_code(510, 0, 1, false)
        );
        assert_eq!(
            environment.bit_id(&unchiral.fingerprint_arguments, &[511, 0], false),
            environment.bit_id(&unchiral.fingerprint_arguments, &[0, 0], false)
        );
        assert_eq!(
            environment.bit_id(&unchiral.fingerprint_arguments, &[2047, u32::MAX], false,),
            get_atom_pair_code(3, 31, 1, false)
        );

        let mut chiral = AtomPairArguments::default();
        chiral.fingerprint_arguments.df_include_chirality = true;
        assert_eq!(
            environment.bit_id(&chiral.fingerprint_arguments, &[2047, u32::MAX], false),
            get_atom_pair_code(0, 1023, 1, true)
        );
    }

    #[test]
    fn environment_bit_id_preserves_modulo_collisions_before_hashing() {
        let arguments = AtomPairArguments::default();
        let environment = AtomPairEnvironment::new(0, 1, 1);
        assert_eq!(
            environment.bit_id(&arguments.fingerprint_arguments, &[0, 0], true),
            environment.bit_id(&arguments.fingerprint_arguments, &[511, 511], true)
        );
        assert_eq!(
            environment.bit_id(&arguments.fingerprint_arguments, &[0, 0], true),
            Ok(4_216_857_148)
        );
    }

    #[test]
    fn environment_bit_id_reports_every_invalid_invariant_length_boundary() {
        let arguments = AtomPairArguments::default();
        assert_eq!(
            AtomPairEnvironment::new(0, 0, 1).bit_id(&arguments.fingerprint_arguments, &[], false,),
            Err(AtomPairCodeError::AtomInvariantLength {
                length: 0,
                atom_index: 0,
            })
        );
        assert_eq!(
            AtomPairEnvironment::new(0, 1, 1).bit_id(&arguments.fingerprint_arguments, &[10], true,),
            Err(AtomPairCodeError::AtomInvariantLength {
                length: 1,
                atom_index: 1,
            })
        );
        assert_eq!(
            AtomPairEnvironment::new(1, 0, 1).bit_id(&arguments.fingerprint_arguments, &[10], true,),
            Err(AtomPairCodeError::AtomInvariantLength {
                length: 1,
                atom_index: 1,
            })
        );
    }

    #[test]
    fn additional_output_each_supported_allocation_receives_the_source_shape() {
        let environment = AtomPairEnvironment::new(0, 2, 3);

        let mut bit_info = AdditionalOutput::new();
        bit_info.allocate_bit_info_map();
        bit_info.reset_for_atom_count(3);
        environment.update_additional_output(&mut bit_info, 17);
        assert_eq!(bit_info.bit_info_map.unwrap().get(&17).unwrap(), &[(0, 2)]);

        let mut atom_to_bits = AdditionalOutput::new();
        atom_to_bits.allocate_atom_to_bits();
        atom_to_bits.reset_for_atom_count(3);
        environment.update_additional_output(&mut atom_to_bits, 17);
        assert_eq!(atom_to_bits.atom_to_bits.unwrap(), [vec![17], vec![], vec![17]]);

        let mut atom_counts = AdditionalOutput::new();
        atom_counts.allocate_atom_counts();
        atom_counts.reset_for_atom_count(3);
        environment.update_additional_output(&mut atom_counts, 17);
        assert_eq!(atom_counts.atom_counts.unwrap(), [1, 0, 1]);

        let mut atoms_per_bit = AdditionalOutput::new();
        atoms_per_bit.allocate_atoms_per_bit();
        atoms_per_bit.reset_for_atom_count(3);
        environment.update_additional_output(&mut atoms_per_bit, 17);
        assert_eq!(atoms_per_bit.atoms_per_bit.unwrap().get(&17).unwrap(), &[vec![0, 2]]);
    }

    #[test]
    fn additional_output_combined_allocations_preserve_pair_and_call_order() {
        let mut output = AdditionalOutput::new();
        output.allocate_bit_info_map();
        output.allocate_atom_to_bits();
        output.allocate_atom_counts();
        output.allocate_atoms_per_bit();
        output.allocate_bit_paths();
        output.reset_for_atom_count(4);

        AtomPairEnvironment::new(0, 2, 3).update_additional_output(&mut output, 17);
        AtomPairEnvironment::new(1, 3, 2).update_additional_output(&mut output, 19);
        AtomPairEnvironment::new(2, 0, 3).update_additional_output(&mut output, 17);

        assert_eq!(
            output.bit_info_map.as_ref().unwrap().get(&17).unwrap(),
            &[(0, 2), (2, 0)]
        );
        assert_eq!(output.bit_info_map.as_ref().unwrap().get(&19).unwrap(), &[(1, 3)]);
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
    fn additional_output_preserves_duplicate_provenance_after_hash_collision() {
        let arguments = AtomPairArguments::default();
        let first = AtomPairEnvironment::new(0, 1, 1);
        let second = AtomPairEnvironment::new(2, 3, 1);
        let invariants = [0, 0, 511, 511];
        let first_bit = first
            .bit_id(&arguments.fingerprint_arguments, &invariants, true)
            .unwrap();
        let second_bit = second
            .bit_id(&arguments.fingerprint_arguments, &invariants, true)
            .unwrap();
        assert_eq!(first_bit, second_bit);

        let mut output = AdditionalOutput::new();
        output.allocate_bit_info_map();
        output.allocate_atom_to_bits();
        output.allocate_atom_counts();
        output.allocate_atoms_per_bit();
        output.reset_for_atom_count(4);
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

    fn environments_2d(
        molecule: &Molecule,
        arguments: &AtomPairArguments,
        from_atoms: Option<&[usize]>,
        ignore_atoms: Option<&[usize]>,
    ) -> Vec<AtomPairEnvironment> {
        AtomPairEnvironmentGenerator
            .environments(molecule, arguments, from_atoms, ignore_atoms, -1)
            .unwrap()
    }

    #[test]
    fn environment_generation_empty_single_chain_branch_ring_and_fused_order() {
        let arguments = AtomPairArguments::default();
        assert!(environments_2d(&Molecule::new(), &arguments, None, None).is_empty());
        assert!(environments_2d(&Molecule::from_smiles("[He]").unwrap(), &arguments, None, None).is_empty());

        let chain = Molecule::from_smiles("CCCC").unwrap();
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

        let branch = Molecule::from_smiles("CC(C)C").unwrap();
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

        let ring = Molecule::from_smiles("C1CC1").unwrap();
        assert_eq!(
            environments_2d(&ring, &arguments, None, None),
            [
                AtomPairEnvironment::new(0, 1, 1),
                AtomPairEnvironment::new(0, 2, 1),
                AtomPairEnvironment::new(1, 2, 1),
            ]
        );

        let fused = Molecule::from_smiles("c1ccc2ccccc2c1").unwrap();
        let fused_environments = environments_2d(&fused, &arguments, None, None);
        assert_eq!(fused_environments.len(), 45);
        assert!(fused_environments.windows(2).all(|pair| {
            (pair[0].atom_id_first, pair[0].atom_id_second) < (pair[1].atom_id_first, pair[1].atom_id_second)
        }));
    }

    #[test]
    fn environment_generation_includes_explicit_hydrogens_and_handles_disconnected_distance() {
        let arguments = AtomPairArguments::default();
        let mut explicit_h_builder = Molecule::builder();
        let carbon = explicit_h_builder.add_atom(AtomSpec::new(Element::C));
        for _ in 0..4 {
            let hydrogen = explicit_h_builder.add_atom(AtomSpec::new(Element::H));
            explicit_h_builder
                .add_bond(BondSpec::new(carbon, hydrogen, BondOrder::Single))
                .unwrap();
        }
        let explicit_h = explicit_h_builder.build().unwrap();
        assert_eq!(explicit_h.num_atoms(), 5);
        assert_eq!(environments_2d(&explicit_h, &arguments, None, None).len(), 10);

        let disconnected = Molecule::from_smiles("C.C").unwrap();
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
        let molecule = Molecule::from_smiles("CCCC").unwrap();
        let arguments = AtomPairArguments::default();
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
            [AtomPairEnvironment::new(1, 2, 1), AtomPairEnvironment::new(2, 3, 1),]
        );
    }

    #[test]
    fn environment_generation_distance_bounds_are_inclusive_and_allow_equality() {
        let molecule = Molecule::from_smiles("CCCCC").unwrap();
        let exactly_two = AtomPairArguments::new(false, false, true, 2, 2, Vec::new(), 2048).unwrap();
        assert_eq!(
            environments_2d(&molecule, &exactly_two, None, None),
            [
                AtomPairEnvironment::new(0, 2, 2),
                AtomPairEnvironment::new(1, 3, 2),
                AtomPairEnvironment::new(2, 4, 2),
            ]
        );

        let one_through_two = AtomPairArguments::new(false, false, true, 1, 2, Vec::new(), 2048).unwrap();
        let distances = environments_2d(&molecule, &one_through_two, None, None)
            .into_iter()
            .map(|environment| environment.distance)
            .collect::<Vec<_>>();
        assert_eq!(distances, [1, 2, 1, 2, 1, 2, 1]);
    }

    fn three_atom_3d_molecule() -> Molecule {
        let mut builder = Molecule::builder();
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
    fn environment_generation_3d_floors_fractional_distances_and_selects_conformer_ids() {
        let molecule = three_atom_3d_molecule();
        let arguments = AtomPairArguments::new(false, false, false, 0, 10, Vec::new(), 2048).unwrap();
        assert_eq!(
            AtomPairEnvironmentGenerator
                .environments(&molecule, &arguments, None, None, -1)
                .unwrap(),
            [
                AtomPairEnvironment::new(0, 1, 1),
                AtomPairEnvironment::new(0, 2, 2),
                AtomPairEnvironment::new(1, 2, 0),
            ]
        );
        assert_eq!(
            AtomPairEnvironmentGenerator
                .environments(&molecule, &arguments, None, None, 2)
                .unwrap(),
            [
                AtomPairEnvironment::new(0, 1, 3),
                AtomPairEnvironment::new(0, 2, 4),
                AtomPairEnvironment::new(1, 2, 5),
            ]
        );
    }

    #[test]
    fn environment_generation_3d_reports_missing_conformers_and_preserves_nonfinite_cast() {
        let mut no_conformer_builder = Molecule::builder();
        no_conformer_builder.add_atom(AtomSpec::new(Element::C));
        let no_conformer = no_conformer_builder.build().unwrap();
        let arguments = AtomPairArguments::new(false, false, false, 0, 30, Vec::new(), 2048).unwrap();
        assert_eq!(
            AtomPairEnvironmentGenerator.environments(&no_conformer, &arguments, None, None, -1),
            Err(AtomPairCodeError::DistanceMatrix(
                DistanceMatrixError::MissingConformer { conf_id: -1 }
            ))
        );

        let molecule = three_atom_3d_molecule();
        assert_eq!(
            AtomPairEnvironmentGenerator.environments(&molecule, &arguments, None, None, 99),
            Err(AtomPairCodeError::DistanceMatrix(
                DistanceMatrixError::MissingConformer { conf_id: 99 }
            ))
        );

        let mut nonfinite_builder = Molecule::builder();
        nonfinite_builder.add_atom(AtomSpec::new(Element::C));
        nonfinite_builder.add_atom(AtomSpec::new(Element::C));
        nonfinite_builder
            .add_conformer(Conformer3D::new(0, vec![[0.0, 0.0, 0.0], [f64::NAN, 0.0, 0.0]], true))
            .unwrap();
        let nonfinite = nonfinite_builder.build().unwrap();
        assert_eq!(
            AtomPairEnvironmentGenerator
                .environments(&nonfinite, &arguments, None, None, -1)
                .unwrap(),
            [AtomPairEnvironment::new(0, 1, 0)]
        );
    }

    #[test]
    fn environment_generator_metadata_result_size_matches_both_source_widths() {
        let generator = AtomPairEnvironmentGenerator;
        let unchiral = AtomPairArguments::default();
        assert_eq!(generator.result_size(&unchiral), 1u64 << 23);

        let mut chiral = unchiral;
        chiral.fingerprint_arguments.df_include_chirality = true;
        assert_eq!(generator.result_size(&chiral), 1u64 << 27);
        assert_eq!(
            generator.result_size(&chiral),
            u64::from(u32::MAX.min((1u32 << 27) - 1)) + 1
        );
    }

    #[test]
    fn environment_generator_metadata_strings_and_json_restore_are_exact_and_stateless() {
        let mut generator = AtomPairEnvironmentGenerator;
        assert_eq!(generator.info_string(), "AtomPairEnvironmentGenerator");
        assert_eq!(generator.to_json(), r#"{"type":"AtomPairEnvGenerator"}"#);
        generator.from_json("").unwrap();
        generator.from_json(generator.to_json()).unwrap();
        generator.from_json(r#"{"type":"OtherGenerator"}"#).unwrap();
        assert_eq!(generator, AtomPairEnvironmentGenerator);
        assert!(generator.from_json("[]").is_err());
        assert!(generator.from_json("{").is_err());
    }

    #[test]
    fn generator_factory_defaults_and_parameter_forwarding_match_source() {
        let default_arguments = AtomPairArguments::default();
        let default_generator = atom_pair_generator(&default_arguments, None, false);
        assert_eq!(default_generator.arguments(), &default_arguments);
        assert!(default_generator.owns_atom_invariants_generator());
        assert_eq!(
            default_generator.atom_invariants_generator,
            AtomPairAtomInvariantsGenerator::new(false, false)
        );

        let parameterized =
            atom_pair_generator_with_parameters(2, 9, true, false, None, false, 4096, vec![2, 3, 7], false).unwrap();
        assert_eq!(parameterized.arguments.min_distance, 2);
        assert_eq!(parameterized.arguments.max_distance, 9);
        assert!(!parameterized.arguments.use_2d);
        assert!(parameterized.arguments.fingerprint_arguments.df_include_chirality);
        assert!(!parameterized.arguments.fingerprint_arguments.df_count_simulation);
        assert_eq!(parameterized.arguments.fingerprint_arguments.d_count_bounds, [2, 3, 7]);
        assert_eq!(parameterized.arguments.fingerprint_arguments.d_fp_size, 4096);
        assert_eq!(
            parameterized.atom_invariants_generator,
            AtomPairAtomInvariantsGenerator::new(true, false)
        );
        assert!(parameterized.owns_atom_invariants_generator());
    }

    #[test]
    fn generator_factory_owns_a_logically_independent_custom_invariant_generator() {
        let arguments = AtomPairArguments::default();
        let mut source_generator = AtomPairAtomInvariantsGenerator::new(true, true);
        let fingerprint_generator = atom_pair_generator(&arguments, Some(source_generator.clone()), false);
        source_generator.include_chirality = false;
        source_generator.topological_torsion_correction = false;

        assert_eq!(
            fingerprint_generator.atom_invariants_generator,
            AtomPairAtomInvariantsGenerator::new(true, true)
        );
        assert!(!fingerprint_generator.owns_atom_invariants_generator());
    }

    #[test]
    fn generator_factory_metadata_and_json_envelope_are_exact() {
        let generator = atom_pair_generator(&AtomPairArguments::default(), None, false);
        assert_eq!(
            generator.info_string(),
            "Common arguments : countSimulation=1 fpSize=2048 bitsPerFeature=1 includeChirality=0 --- AtomPairArguments use2D=1 minDistance=1 maxDistance=30 --- AtomPairEnvironmentGenerator --- AtomPairInvariantGenerator topologicalTorsionCorrection=0 --- No bond invariants generator"
        );
        assert_eq!(
            generator.to_json(),
            r#"{"name":"FingerprintGenerator","fingerprintArguments":{"type":"AtomPairArguments","use2D":"true","minDistance":"1","maxDistance":"30","countSimulation":"true","fpSize":"2048","numBitsPerFeature":"1","includeChirality":"false","countBounds":["1","2","4","8"]},"atomEnvironmentGenerator":{"type":"AtomPairEnvGenerator"},"atomInvariantsGenerator":{"type":"AtomPairAtomInvGenerator","includeChirality":"false","topologicalTorsionCorrection":"false"}}"#
        );
        let core = generator::FingerprintGenerator::new(&generator);
        core.from_json("").unwrap();
        core.from_json(r#"{"ignored":true}"#).unwrap();
        assert!(core.from_json("[]").is_err());
        assert!(core.from_json("{").is_err());
    }

    #[test]
    fn generator_factory_full_and_partial_json_restore_atom_pair_dispatch() {
        let original = atom_pair_generator(
            &AtomPairArguments::new(false, true, true, 2, 8, vec![1, 3], 1024).unwrap(),
            Some(AtomPairAtomInvariantsGenerator::new(false, true)),
            false,
        );
        let restored = generator::generator_from_json(&original.to_json()).unwrap();
        let generator::RestoredFingerprintGenerator::AtomPair(restored) = restored else {
            panic!("AtomPair JSON restored a different fingerprint family");
        };
        assert_eq!(restored.arguments, original.arguments);
        assert_eq!(restored.atom_invariants_generator, original.atom_invariants_generator);
        assert_eq!(restored.to_json(), original.to_json());

        let partial = r#"{
          "fingerprintArguments": {
            "type": "AtomPairArguments",
            "minDistance": 4,
            "countBounds": [1, 2, 4, 8]
          },
          "atomEnvironmentGenerator": {"type": "AtomPairEnvGenerator"}
        }"#;
        let partial = generator::generator_from_json(partial).unwrap();
        let generator::RestoredFingerprintGenerator::AtomPair(partial) = partial else {
            panic!("partial AtomPair JSON restored a different fingerprint family");
        };
        assert_eq!(partial.arguments.min_distance, 4);
        assert_eq!(partial.arguments.max_distance, 30);
        assert!(partial.arguments.use_2d);
        assert_eq!(
            partial.atom_invariants_generator,
            AtomPairAtomInvariantsGenerator::default()
        );
        assert!(partial.owns_atom_invariants_generator());

        for invalid in [
            r#"{"atomEnvironmentGenerator":{"type":"AtomPairEnvGenerator"}}"#,
            r#"{"fingerprintArguments":{"type":"AtomPairArguments"},"atomEnvironmentGenerator":{"type":"MorganEnvGenerator"}}"#,
            r#"{"fingerprintArguments":{"type":"Unknown"},"atomEnvironmentGenerator":{"type":"AtomPairEnvGenerator"}}"#,
            r#"{"fingerprintArguments":{"type":"AtomPairArguments"},"atomEnvironmentGenerator":{"type":"AtomPairEnvGenerator"},"atomInvariantsGenerator":{"type":"Unknown"}}"#,
            r#"{"fingerprintArguments":{"type":"AtomPairArguments"},"atomEnvironmentGenerator":{"type":"AtomPairEnvGenerator"},"bondInvariantsGenerator":{"type":"MorganBondInvGenerator"}}"#,
        ] {
            assert!(generator::generator_from_json(invalid).is_err(), "{invalid}");
        }
    }

    #[test]
    fn generator_factory_restored_generator_runs_all_four_shared_output_forms() {
        let molecule = Molecule::from_smiles("CCCO").unwrap();
        let original = atom_pair_generator(&AtomPairArguments::default(), None, false);
        let restored = generator::generator_from_json(&original.to_json()).unwrap();

        assert_eq!(
            original
                .sparse_count_fingerprint(&molecule, &mut FingerprintFuncArguments::default())
                .unwrap(),
            restored
                .sparse_count_fingerprint(&molecule, &mut FingerprintFuncArguments::default())
                .unwrap()
        );
        assert_eq!(
            original
                .count_fingerprint(&molecule, &mut FingerprintFuncArguments::default())
                .unwrap(),
            restored
                .folded_count_fingerprint(&molecule, &mut FingerprintFuncArguments::default())
                .unwrap()
        );
        assert_eq!(
            original
                .sparse_bit_fingerprint(&molecule, &mut FingerprintFuncArguments::default())
                .unwrap(),
            restored
                .sparse_bit_fingerprint(&molecule, &mut FingerprintFuncArguments::default())
                .unwrap()
        );
        assert_eq!(
            original
                .fingerprint(&molecule, &mut FingerprintFuncArguments::default())
                .unwrap(),
            restored
                .explicit_bit_fingerprint(&molecule, &mut FingerprintFuncArguments::default())
                .unwrap()
        );
        assert_eq!(restored.info_string(), original.info_string());
        assert_eq!(restored.to_json(), original.to_json());
    }

    #[test]
    fn legacy_adapters_delegate_default_and_both_count_simulation_modes_to_modern_core() {
        let molecule = Molecule::from_smiles("CCCO").unwrap();
        let legacy_arguments = LegacyAtomPairArguments::default();

        let legacy_sparse = legacy_sparse_count_fingerprint(&molecule, &legacy_arguments).unwrap();
        let modern_default = AtomPairFingerprintGenerator::new(&AtomPairFingerprintParams::default()).unwrap();
        let modern_sparse = modern_default
            .sparse_count_fingerprint(&molecule, &mut FingerprintFuncArguments::default())
            .unwrap();
        assert_eq!(legacy_sparse, modern_sparse);

        let legacy_count = legacy_hashed_count_fingerprint(&molecule, 128, &legacy_arguments).unwrap();
        let modern_count = AtomPairFingerprintGenerator::new(&AtomPairFingerprintParams {
            n_bits: 128,
            ..Default::default()
        })
        .unwrap()
        .count_fingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .unwrap();
        assert_eq!(legacy_count, modern_count);

        for (n_bits, n_bits_per_entry, count_bounds) in [(128, 4, vec![1, 2, 4, 8]), (96, 3, vec![1, 2, 3])] {
            let legacy = legacy_hashed_bit_fingerprint(&molecule, n_bits, n_bits_per_entry, &legacy_arguments).unwrap();
            let modern = AtomPairFingerprintGenerator::new(&AtomPairFingerprintParams {
                n_bits: n_bits as usize,
                count_bounds,
                ..Default::default()
            })
            .unwrap()
            .fingerprint(&molecule, &mut FingerprintFuncArguments::default())
            .unwrap();
            assert_eq!(legacy, modern);
        }
    }

    #[test]
    fn legacy_adapters_delegate_custom_chiral_rooted_and_ignored_pairs_to_modern_core() {
        let molecule = Molecule::from_smiles("C[C@H](O)F").unwrap();
        let legacy_arguments = LegacyAtomPairArguments {
            min_length: 1,
            max_length: 2,
            from_atoms: Some(vec![1]),
            ignore_atoms: Some(vec![3]),
            atom_invariants: Some(vec![11, 22, 33, 44]),
            include_chirality: true,
            ..Default::default()
        };
        let modern_generator = AtomPairFingerprintGenerator::new(&AtomPairFingerprintParams {
            n_bits: 64,
            min_distance: 1,
            max_distance: 2,
            use_chirality: true,
            ..Default::default()
        })
        .unwrap();
        let modern_arguments = || FingerprintFuncArguments {
            from_atoms: Some(vec![1]),
            ignore_atoms: Some(vec![3]),
            custom_atom_invariants: Some(vec![11, 22, 33, 44]),
            ..Default::default()
        };

        assert_eq!(
            legacy_hashed_count_fingerprint(&molecule, 64, &legacy_arguments).unwrap(),
            modern_generator
                .count_fingerprint(&molecule, &mut modern_arguments())
                .unwrap()
        );
        assert_eq!(
            legacy_hashed_bit_fingerprint(&molecule, 64, 4, &legacy_arguments).unwrap(),
            modern_generator
                .fingerprint(&molecule, &mut modern_arguments())
                .unwrap()
        );
    }

    #[test]
    fn legacy_adapters_delegate_three_dimensional_conformer_selection_to_modern_core() {
        let molecule = Molecule::from_mol_block(
            "\n     RDKit          3D\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.6000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\nM  END\n",
        )
        .unwrap();
        let legacy_arguments = LegacyAtomPairArguments {
            min_length: 2,
            max_length: 3,
            use_2d: false,
            conf_id: 0,
            ..Default::default()
        };
        let modern = AtomPairFingerprintGenerator::new(&AtomPairFingerprintParams {
            n_bits: 32,
            min_distance: 2,
            max_distance: 3,
            use_2d: false,
            ..Default::default()
        })
        .unwrap();
        let args = || FingerprintFuncArguments {
            conf_id: 0,
            ..Default::default()
        };

        assert_eq!(
            legacy_hashed_count_fingerprint(&molecule, 32, &legacy_arguments).unwrap(),
            modern.count_fingerprint(&molecule, &mut args()).unwrap()
        );
        assert_eq!(
            legacy_hashed_bit_fingerprint(&molecule, 32, 4, &legacy_arguments).unwrap(),
            modern.fingerprint(&molecule, &mut args()).unwrap()
        );
    }

    #[test]
    fn legacy_adapters_reject_source_precondition_violations_as_typed_errors() {
        let molecule = Molecule::from_smiles("CCO").unwrap();
        let too_few_invariants = LegacyAtomPairArguments {
            atom_invariants: Some(vec![1, 2]),
            ..Default::default()
        };
        assert!(matches!(
            legacy_sparse_count_fingerprint(&molecule, &too_few_invariants),
            Err(FingerprintError::InvalidArguments { .. })
        ));
        assert!(matches!(
            legacy_hashed_bit_fingerprint(&molecule, 128, 0, &LegacyAtomPairArguments::default()),
            Err(FingerprintError::InvalidArguments { .. })
        ));
    }
}
