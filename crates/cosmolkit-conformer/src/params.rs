//! Complete source-backed embedding options, factories, JSON and failure vocabulary.
use crate::ConformerError;
use crate::bounds::BoundsMatrix;
use std::collections::BTreeMap;
use std::sync::Arc;

// BEGIN RDKIT CPP ENUM DGeomHelpers::EmbedFailureCauses (Embedder.h:25-39)
// RDKit✔️✔️: enum EmbedFailureCauses {
// RDKit✔️✔️:   INITIAL_COORDS = 0,
// RDKit✔️✔️:   FIRST_MINIMIZATION = 1,
// RDKit✔️✔️:   CHECK_TETRAHEDRAL_CENTERS = 2,
// RDKit✔️✔️:   CHECK_CHIRAL_CENTERS = 3,
// RDKit✔️✔️:   MINIMIZE_FOURTH_DIMENSION = 4,
// RDKit✔️✔️:   ETK_MINIMIZATION = 5,
// RDKit✔️✔️:   FINAL_CHIRAL_BOUNDS = 6,
// RDKit✔️✔️:   FINAL_CENTER_IN_VOLUME = 7,
// RDKit✔️✔️:   LINEAR_DOUBLE_BOND = 8,
// RDKit✔️✔️:   BAD_DOUBLE_BOND_STEREO = 9,
// RDKit✔️✔️:   CHECK_CHIRAL_CENTERS2 = 10,
// RDKit✔️✔️:   EXCEEDED_TIMEOUT = 11,
// RDKit✔️✔️:   END_OF_ENUM = 12,
// RDKit✔️✔️: };
// END RDKIT CPP ENUM DGeomHelpers::EmbedFailureCauses
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u32)]
pub enum EmbedFailureCause {
    InitialCoords = 0,
    FirstMinimization = 1,
    CheckTetrahedralCenters = 2,
    CheckChiralCenters = 3,
    MinimizeFourthDimension = 4,
    EtkMinimization = 5,
    FinalChiralBounds = 6,
    FinalCenterInVolume = 7,
    LinearDoubleBond = 8,
    BadDoubleBondStereo = 9,
    CheckChiralCenters2 = 10,
    ExceededTimeout = 11,
    EndOfEnum = 12,
}

impl EmbedFailureCause {
    pub const ALL: [Self; 13] = [
        Self::InitialCoords,
        Self::FirstMinimization,
        Self::CheckTetrahedralCenters,
        Self::CheckChiralCenters,
        Self::MinimizeFourthDimension,
        Self::EtkMinimization,
        Self::FinalChiralBounds,
        Self::FinalCenterInVolume,
        Self::LinearDoubleBond,
        Self::BadDoubleBondStereo,
        Self::CheckChiralCenters2,
        Self::ExceededTimeout,
        Self::EndOfEnum,
    ];

    #[must_use]
    pub const fn rdkit_ordinal(self) -> u32 {
        self as u32
    }

    #[must_use]
    pub const fn from_rdkit_ordinal(value: u32) -> Option<Self> {
        match value {
            0 => Some(Self::InitialCoords),
            1 => Some(Self::FirstMinimization),
            2 => Some(Self::CheckTetrahedralCenters),
            3 => Some(Self::CheckChiralCenters),
            4 => Some(Self::MinimizeFourthDimension),
            5 => Some(Self::EtkMinimization),
            6 => Some(Self::FinalChiralBounds),
            7 => Some(Self::FinalCenterInVolume),
            8 => Some(Self::LinearDoubleBond),
            9 => Some(Self::BadDoubleBondStereo),
            10 => Some(Self::CheckChiralCenters2),
            11 => Some(Self::ExceededTimeout),
            12 => Some(Self::EndOfEnum),
            _ => None,
        }
    }

    #[must_use]
    pub const fn rdkit_name(self) -> &'static str {
        match self {
            Self::InitialCoords => "INITIAL_COORDS",
            Self::FirstMinimization => "FIRST_MINIMIZATION",
            Self::CheckTetrahedralCenters => "CHECK_TETRAHEDRAL_CENTERS",
            Self::CheckChiralCenters => "CHECK_CHIRAL_CENTERS",
            Self::MinimizeFourthDimension => "MINIMIZE_FOURTH_DIMENSION",
            Self::EtkMinimization => "ETK_MINIMIZATION",
            Self::FinalChiralBounds => "FINAL_CHIRAL_BOUNDS",
            Self::FinalCenterInVolume => "FINAL_CENTER_IN_VOLUME",
            Self::LinearDoubleBond => "LINEAR_DOUBLE_BOND",
            Self::BadDoubleBondStereo => "BAD_DOUBLE_BOND_STEREO",
            Self::CheckChiralCenters2 => "CHECK_CHIRAL_CENTERS2",
            Self::ExceededTimeout => "EXCEEDED_TIMEOUT",
            Self::EndOfEnum => "END_OF_ENUM",
        }
    }
}

// BEGIN RDKIT CPP STRUCT DGeomHelpers::EmbedParameters (Embedder.h:122-191)
// RDKit✔️✔️: struct RDKIT_DISTGEOMHELPERS_EXPORT EmbedParameters {
// RDKit✔️✔️:   unsigned int maxIterations{0};
// RDKit✔️✔️:   int numThreads{1};
// RDKit✔️✔️:   int randomSeed{-1};
// RDKit✔️✔️:   bool clearConfs{true};
// RDKit✔️✔️:   bool useRandomCoords{false};
// RDKit✔️✔️:   double boxSizeMult{2.0};
// RDKit✔️✔️:   bool randNegEig{true};
// RDKit✔️✔️:   unsigned int numZeroFail{1};
// RDKit✔️✔️:   const std::map<int, RDGeom::Point3D> *coordMap{nullptr};
// RDKit✔️✔️:   double optimizerForceTol{1e-3};
// RDKit✔️✔️:   bool ignoreSmoothingFailures{false};
// RDKit✔️✔️:   bool enforceChirality{true};
// RDKit✔️✔️:   bool useExpTorsionAnglePrefs{false};
// RDKit✔️✔️:   bool useBasicKnowledge{false};
// RDKit✔️✔️:   bool verbose{false};
// RDKit✔️✔️:   double basinThresh{5.0};
// RDKit✔️✔️:   double pruneRmsThresh{-1.0};
// RDKit✔️✔️:   bool onlyHeavyAtomsForRMS{true};
// RDKit✔️✔️:   unsigned int ETversion{2};
// RDKit✔️✔️:   boost::shared_ptr<const DistGeom::BoundsMatrix> boundsMat;
// RDKit✔️✔️:   bool embedFragmentsSeparately{true};
// RDKit✔️✔️:   bool useSmallRingTorsions{false};
// RDKit✔️✔️:   bool useMacrocycleTorsions{false};
// RDKit✔️✔️:   bool useMacrocycle14config{false};
// RDKit✔️✔️:   unsigned int timeout{0};
// RDKit✔️✔️:   std::shared_ptr<std::map<std::pair<unsigned int, unsigned int>, double>> CPCI;
// RDKit✔️✔️:   void (*callback)(unsigned int);
// RDKit✔️✔️:   bool forceTransAmides{true};
// RDKit✔️✔️:   bool useSymmetryForPruning{true};
// RDKit✔️✔️:   double boundsMatForceScaling{1.0};
// RDKit✔️✔️:   bool trackFailures{false};
// RDKit✔️✔️:   std::vector<unsigned int> failures;
// RDKit✔️✔️:   bool enableSequentialRandomSeeds{false};
// RDKit✔️✔️:   bool symmetrizeConjugatedTerminalGroupsForPruning{true};
#[derive(Clone)]
pub struct EmbedParams {
    pub max_iterations: u32,
    pub num_threads: i32,
    pub random_seed: i32,
    pub clear_conformers: bool,
    pub use_random_coords: bool,
    pub box_size_mult: f64,
    pub rand_neg_eig: bool,
    pub num_zero_fail: u32,
    pub coord_map: Option<BTreeMap<i32, [f64; 3]>>,
    pub optimizer_force_tol: f64,
    pub ignore_smoothing_failures: bool,
    pub enforce_chirality: bool,
    pub use_exp_torsion_angle_prefs: bool,
    pub use_basic_knowledge: bool,
    pub verbose: bool,
    pub basin_thresh: f64,
    pub prune_rms_thresh: f64,
    pub only_heavy_atoms_for_rms: bool,
    pub et_version: u32,
    pub(crate) bounds_mat: Option<Arc<BoundsMatrix>>,
    pub embed_fragments_separately: bool,
    pub use_small_ring_torsions: bool,
    pub use_macrocycle_torsions: bool,
    pub use_macrocycle14config: bool,
    pub timeout: u32,
    pub cpci: Option<BTreeMap<(u32, u32), f64>>,
    pub callback: Option<fn(u32)>,
    pub force_trans_amides: bool,
    pub use_symmetry_for_pruning: bool,
    pub bounds_mat_force_scaling: f64,
    pub track_failures: bool,
    pub failures: Vec<u32>,
    pub enable_sequential_random_seeds: bool,
    pub symmetrize_conjugated_terminal_groups_for_pruning: bool,
}

impl Default for EmbedParams {
    fn default() -> Self {
        // RDKit✔️✔️:   EmbedParameters() : boundsMat(nullptr), CPCI(nullptr), callback(nullptr) {}
        Self {
            max_iterations: 0,
            num_threads: 1,
            random_seed: -1,
            clear_conformers: true,
            use_random_coords: false,
            box_size_mult: 2.0,
            rand_neg_eig: true,
            num_zero_fail: 1,
            coord_map: None,
            optimizer_force_tol: 1e-3,
            ignore_smoothing_failures: false,
            enforce_chirality: true,
            use_exp_torsion_angle_prefs: false,
            use_basic_knowledge: false,
            verbose: false,
            basin_thresh: 5.0,
            prune_rms_thresh: -1.0,
            only_heavy_atoms_for_rms: true,
            et_version: 2,
            bounds_mat: None,
            embed_fragments_separately: true,
            use_small_ring_torsions: false,
            use_macrocycle_torsions: false,
            use_macrocycle14config: false,
            timeout: 0,
            cpci: None,
            callback: None,
            force_trans_amides: true,
            use_symmetry_for_pruning: true,
            bounds_mat_force_scaling: 1.0,
            track_failures: false,
            failures: Vec::new(),
            enable_sequential_random_seeds: false,
            symmetrize_conjugated_terminal_groups_for_pruning: true,
        }
    }
}

impl EmbedParams {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    pub fn set_bounds_matrix(&mut self, bounds: Vec<Vec<f64>>) -> Result<(), ConformerError> {
        let n = bounds.len();
        if bounds.iter().any(|row| row.len() != n) {
            return Err(ConformerError::GenerationFailed(
                "bounds matrix must be square".to_string(),
            ));
        }
        for row in &bounds {
            for &value in row {
                if !value.is_finite() || value < 0.0 {
                    return Err(ConformerError::GenerationFailed(
                        "bounds matrix values must be finite and non-negative".to_string(),
                    ));
                }
            }
        }
        self.bounds_mat = Some(Arc::new(
            BoundsMatrix::from_data(n, bounds.into_iter().flatten().collect())
                .map_err(|error| ConformerError::GenerationFailed(format!("{error:?}")))?,
        ));
        Ok(())
    }

    #[must_use]
    pub fn has_bounds_matrix(&self) -> bool {
        self.bounds_mat.is_some()
    }

    #[allow(clippy::too_many_arguments)]
    fn from_rdkit_constructor(
        max_iterations: u32,
        num_threads: i32,
        random_seed: i32,
        clear_conformers: bool,
        use_random_coords: bool,
        box_size_mult: f64,
        rand_neg_eig: bool,
        num_zero_fail: u32,
        coord_map: Option<BTreeMap<i32, [f64; 3]>>,
        optimizer_force_tol: f64,
        ignore_smoothing_failures: bool,
        enforce_chirality: bool,
        use_exp_torsion_angle_prefs: bool,
        use_basic_knowledge: bool,
        verbose: bool,
        basin_thresh: f64,
        prune_rms_thresh: f64,
        only_heavy_atoms_for_rms: bool,
        et_version: u32,
        bounds_mat: Option<Arc<BoundsMatrix>>,
        embed_fragments_separately: bool,
        use_small_ring_torsions: bool,
        use_macrocycle_torsions: bool,
        use_macrocycle14config: bool,
        timeout: u32,
        cpci: Option<BTreeMap<(u32, u32), f64>>,
        callback: Option<fn(u32)>,
    ) -> Self {
        // RDKit✔️✔️:       : maxIterations(maxIterations),
        // RDKit✔️✔️:         numThreads(numThreads),
        // RDKit✔️✔️:         randomSeed(randomSeed),
        // RDKit✔️✔️:         clearConfs(clearConfs),
        // RDKit✔️✔️:         useRandomCoords(useRandomCoords),
        // RDKit✔️✔️:         boxSizeMult(boxSizeMult),
        // RDKit✔️✔️:         randNegEig(randNegEig),
        // RDKit✔️✔️:         numZeroFail(numZeroFail),
        // RDKit✔️✔️:         coordMap(coordMap),
        // RDKit✔️✔️:         optimizerForceTol(optimizerForceTol),
        // RDKit✔️✔️:         ignoreSmoothingFailures(ignoreSmoothingFailures),
        // RDKit✔️✔️:         enforceChirality(enforceChirality),
        // RDKit✔️✔️:         useExpTorsionAnglePrefs(useExpTorsionAnglePrefs),
        // RDKit✔️✔️:         useBasicKnowledge(useBasicKnowledge),
        // RDKit✔️✔️:         verbose(verbose),
        // RDKit✔️✔️:         basinThresh(basinThresh),
        // RDKit✔️✔️:         pruneRmsThresh(pruneRmsThresh),
        // RDKit✔️✔️:         onlyHeavyAtomsForRMS(onlyHeavyAtomsForRMS),
        // RDKit✔️✔️:         ETversion(ETversion),
        // RDKit✔️✔️:         boundsMat(boundsMat),
        // RDKit✔️✔️:         embedFragmentsSeparately(embedFragmentsSeparately),
        // RDKit✔️✔️:         useSmallRingTorsions(useSmallRingTorsions),
        // RDKit✔️✔️:         useMacrocycleTorsions(useMacrocycleTorsions),
        // RDKit✔️✔️:         useMacrocycle14config(useMacrocycle14config),
        // RDKit✔️✔️:         timeout(timeout),
        // RDKit✔️✔️:         CPCI(std::move(CPCI)),
        // RDKit✔️✔️:         callback(callback) {}
        let mut params = Self {
            max_iterations,
            num_threads,
            random_seed,
            clear_conformers,
            use_random_coords,
            box_size_mult,
            rand_neg_eig,
            num_zero_fail,
            coord_map,
            optimizer_force_tol,
            ignore_smoothing_failures,
            enforce_chirality,
            use_exp_torsion_angle_prefs,
            use_basic_knowledge,
            verbose,
            basin_thresh,
            prune_rms_thresh,
            only_heavy_atoms_for_rms,
            et_version,
            bounds_mat,
            embed_fragments_separately,
            use_small_ring_torsions,
            use_macrocycle_torsions,
            use_macrocycle14config,
            timeout,
            cpci,
            callback,
            ..Self::default()
        };
        params.failures = Vec::new();
        params
    }

    #[must_use]
    pub fn kdg() -> Self {
        // RDKit✔️✔️: const EmbedParameters KDG(0,        // maxIterations
        // RDKit✔️✔️:                           1,        // numThreads
        // RDKit✔️✔️:                           -1,       // randomSeed
        // RDKit✔️✔️:                           true,     // clearConfs
        // RDKit✔️✔️:                           false,    // useRandomCoords
        // RDKit✔️✔️:                           2.0,      // boxSizeMult
        // RDKit✔️✔️:                           true,     // randNegEig
        // RDKit✔️✔️:                           1,        // numZeroFail
        // RDKit✔️✔️:                           nullptr,  // coordMap
        // RDKit✔️✔️:                           1e-3,     // optimizerForceTol
        // RDKit✔️✔️:                           false,    // ignoreSmoothingFailures
        // RDKit✔️✔️:                           true,     // enforceChirality
        // RDKit✔️✔️:                           false,    // useExpTorsionAnglePrefs
        // RDKit✔️✔️:                           true,     // useBasicKnowledge
        // RDKit✔️✔️:                           false,    // verbose
        // RDKit✔️✔️:                           5.0,      // basinThresh
        // RDKit✔️✔️:                           -1.0,     // pruneRmsThresh
        // RDKit✔️✔️:                           true,     // onlyHeavyAtomsForRMS
        // RDKit✔️✔️:                           1,        // ETversion
        // RDKit✔️✔️:                           nullptr,  // boundsMat
        // RDKit✔️✔️:                           true,     // embedFragmentsSeparately
        // RDKit✔️✔️:                           false,    // useSmallRingTorsions
        // RDKit✔️✔️:                           false,    // useMacrocycleTorsions
        // RDKit✔️✔️:                           false,    // useMacrocycle14config
        // RDKit✔️✔️:                           0,        // timeout
        // RDKit✔️✔️:                           nullptr,  // CPCI
        // RDKit✔️✔️:                           nullptr   // callback
        Self::from_rdkit_constructor(
            0, 1, -1, true, false, 2.0, true, 1, None, 1e-3, false, true, false, true, false, 5.0,
            -1.0, true, 1, None, true, false, false, false, 0, None, None,
        )
    }

    #[must_use]
    pub fn etdg() -> Self {
        // RDKit✔️✔️: const EmbedParameters ETDG(0,        // maxIterations
        // RDKit✔️✔️:                            1,        // numThreads
        // RDKit✔️✔️:                            -1,       // randomSeed
        // RDKit✔️✔️:                            true,     // clearConfs
        // RDKit✔️✔️:                            false,    // useRandomCoords
        // RDKit✔️✔️:                            2.0,      // boxSizeMult
        // RDKit✔️✔️:                            true,     // randNegEig
        // RDKit✔️✔️:                            1,        // numZeroFail
        // RDKit✔️✔️:                            nullptr,  // coordMap
        // RDKit✔️✔️:                            1e-3,     // optimizerForceTol
        // RDKit✔️✔️:                            false,    // ignoreSmoothingFailures
        // RDKit✔️✔️:                            false,    // enforceChirality
        // RDKit✔️✔️:                            true,     // useExpTorsionAnglePrefs
        // RDKit✔️✔️:                            false,    // useBasicKnowledge
        // RDKit✔️✔️:                            false,    // verbose
        // RDKit✔️✔️:                            5.0,      // basinThresh
        // RDKit✔️✔️:                            -1.0,     // pruneRmsThresh
        // RDKit✔️✔️:                            true,     // onlyHeavyAtomsForRMS
        // RDKit✔️✔️:                            1,        // ETversion
        // RDKit✔️✔️:                            nullptr,  // boundsMat
        // RDKit✔️✔️:                            true,     // embedFragmentsSeparately
        // RDKit✔️✔️:                            false,    // useSmallRingTorsions
        // RDKit✔️✔️:                            false,    // useMacrocycleTorsions
        // RDKit✔️✔️:                            false,    // useMacrocycle14config
        // RDKit✔️✔️:                            0,        // timeout
        // RDKit✔️✔️:                            nullptr,  // CPCI
        // RDKit✔️✔️:                            nullptr   // callback
        Self::from_rdkit_constructor(
            0, 1, -1, true, false, 2.0, true, 1, None, 1e-3, false, false, true, false, false, 5.0,
            -1.0, true, 1, None, true, false, false, false, 0, None, None,
        )
    }

    #[must_use]
    pub fn etdg_v2() -> Self {
        // RDKit✔️✔️: const EmbedParameters ETDGv2(0,        // maxIterations
        // RDKit✔️✔️:                              1,        // numThreads
        // RDKit✔️✔️:                              -1,       // randomSeed
        // RDKit✔️✔️:                              true,     // clearConfs
        // RDKit✔️✔️:                              false,    // useRandomCoords
        // RDKit✔️✔️:                              2.0,      // boxSizeMult
        // RDKit✔️✔️:                              true,     // randNegEig
        // RDKit✔️✔️:                              1,        // numZeroFail
        // RDKit✔️✔️:                              nullptr,  // coordMap
        // RDKit✔️✔️:                              1e-3,     // optimizerForceTol
        // RDKit✔️✔️:                              false,    // ignoreSmoothingFailures
        // RDKit✔️✔️:                              false,    // enforceChirality
        // RDKit✔️✔️:                              true,     // useExpTorsionAnglePrefs
        // RDKit✔️✔️:                              false,    // useBasicKnowledge
        // RDKit✔️✔️:                              false,    // verbose
        // RDKit✔️✔️:                              5.0,      // basinThresh
        // RDKit✔️✔️:                              -1.0,     // pruneRmsThresh
        // RDKit✔️✔️:                              true,     // onlyHeavyAtomsForRMS
        // RDKit✔️✔️:                              2,        // ETversion
        // RDKit✔️✔️:                              nullptr,  // boundsMat
        // RDKit✔️✔️:                              true,     // embedFragmentsSeparately
        // RDKit✔️✔️:                              false,    // useSmallRingTorsions
        // RDKit✔️✔️:                              false,    // useMacrocycleTorsions
        // RDKit✔️✔️:                              false,    // useMacrocycle14config
        // RDKit✔️✔️:                              0,        // timeout
        // RDKit✔️✔️:                              nullptr,  // CPCI
        // RDKit✔️✔️:                              nullptr   // callback
        Self::from_rdkit_constructor(
            0, 1, -1, true, false, 2.0, true, 1, None, 1e-3, false, false, true, false, false, 5.0,
            -1.0, true, 2, None, true, false, false, false, 0, None, None,
        )
    }

    #[must_use]
    pub fn etkdg() -> Self {
        // RDKit✔️✔️: const EmbedParameters ETKDG(0,        // maxIterations
        // RDKit✔️✔️:                             1,        // numThreads
        // RDKit✔️✔️:                             -1,       // randomSeed
        // RDKit✔️✔️:                             true,     // clearConfs
        // RDKit✔️✔️:                             false,    // useRandomCoords
        // RDKit✔️✔️:                             2.0,      // boxSizeMult
        // RDKit✔️✔️:                             true,     // randNegEig
        // RDKit✔️✔️:                             1,        // numZeroFail
        // RDKit✔️✔️:                             nullptr,  // coordMap
        // RDKit✔️✔️:                             1e-3,     // optimizerForceTol
        // RDKit✔️✔️:                             false,    // ignoreSmoothingFailures
        // RDKit✔️✔️:                             true,     // enforceChirality
        // RDKit✔️✔️:                             true,     // useExpTorsionAnglePrefs
        // RDKit✔️✔️:                             true,     // useBasicKnowledge
        // RDKit✔️✔️:                             false,    // verbose
        // RDKit✔️✔️:                             5.0,      // basinThresh
        // RDKit✔️✔️:                             -1.0,     // pruneRmsThresh
        // RDKit✔️✔️:                             true,     // onlyHeavyAtomsForRMS
        // RDKit✔️✔️:                             1,        // ETversion
        // RDKit✔️✔️:                             nullptr,  // boundsMat
        // RDKit✔️✔️:                             true,     // embedFragmentsSeparately
        // RDKit✔️✔️:                             false,    // useSmallRingTorsions
        // RDKit✔️✔️:                             false,    // useMacrocycleTorsions
        // RDKit✔️✔️:                             false,    // useMacrocycle14config
        // RDKit✔️✔️:                             0,        // timeout
        // RDKit✔️✔️:                             nullptr,  // CPCI
        // RDKit✔️✔️:                             nullptr   // callback
        Self::from_rdkit_constructor(
            0, 1, -1, true, false, 2.0, true, 1, None, 1e-3, false, true, true, true, false, 5.0,
            -1.0, true, 1, None, true, false, false, false, 0, None, None,
        )
    }

    #[must_use]
    pub fn etkdg_v2() -> Self {
        // RDKit✔️✔️: const EmbedParameters ETKDGv2(0,        // maxIterations
        // RDKit✔️✔️:                               1,        // numThreads
        // RDKit✔️✔️:                               -1,       // randomSeed
        // RDKit✔️✔️:                               true,     // clearConfs
        // RDKit✔️✔️:                               false,    // useRandomCoords
        // RDKit✔️✔️:                               2.0,      // boxSizeMult
        // RDKit✔️✔️:                               true,     // randNegEig
        // RDKit✔️✔️:                               1,        // numZeroFail
        // RDKit✔️✔️:                               nullptr,  // coordMap
        // RDKit✔️✔️:                               1e-3,     // optimizerForceTol
        // RDKit✔️✔️:                               false,    // ignoreSmoothingFailures
        // RDKit✔️✔️:                               true,     // enforceChirality
        // RDKit✔️✔️:                               true,     // useExpTorsionAnglePrefs
        // RDKit✔️✔️:                               true,     // useBasicKnowledge
        // RDKit✔️✔️:                               false,    // verbose
        // RDKit✔️✔️:                               5.0,      // basinThresh
        // RDKit✔️✔️:                               -1.0,     // pruneRmsThresh
        // RDKit✔️✔️:                               true,     // onlyHeavyAtomsForRMS
        // RDKit✔️✔️:                               2,        // ETversion
        // RDKit✔️✔️:                               nullptr,  // boundsMat
        // RDKit✔️✔️:                               true,     // embedFragmentsSeparately
        // RDKit✔️✔️:                               false,    // useSmallRingTorsions
        // RDKit✔️✔️:                               false,    // useMacrocycleTorsions
        // RDKit✔️✔️:                               false,    // useMacrocycle14config
        // RDKit✔️✔️:                               0,        // timeout
        // RDKit✔️✔️:                               nullptr,  // CPCI
        // RDKit✔️✔️:                               nullptr   // callback
        Self::from_rdkit_constructor(
            0, 1, -1, true, false, 2.0, true, 1, None, 1e-3, false, true, true, true, false, 5.0,
            -1.0, true, 2, None, true, false, false, false, 0, None, None,
        )
    }

    #[must_use]
    pub fn etkdg_v3() -> Self {
        // RDKit✔️✔️: const EmbedParameters ETKDGv3(0,        // maxIterations
        // RDKit✔️✔️:                               1,        // numThreads
        // RDKit✔️✔️:                               -1,       // randomSeed
        // RDKit✔️✔️:                               true,     // clearConfs
        // RDKit✔️✔️:                               false,    // useRandomCoords
        // RDKit✔️✔️:                               2.0,      // boxSizeMult
        // RDKit✔️✔️:                               true,     // randNegEig
        // RDKit✔️✔️:                               1,        // numZeroFail
        // RDKit✔️✔️:                               nullptr,  // coordMap
        // RDKit✔️✔️:                               1e-3,     // optimizerForceTol
        // RDKit✔️✔️:                               false,    // ignoreSmoothingFailures
        // RDKit✔️✔️:                               true,     // enforceChirality
        // RDKit✔️✔️:                               true,     // useExpTorsionAnglePrefs
        // RDKit✔️✔️:                               true,     // useBasicKnowledge
        // RDKit✔️✔️:                               false,    // verbose
        // RDKit✔️✔️:                               5.0,      // basinThresh
        // RDKit✔️✔️:                               -1.0,     // pruneRmsThresh
        // RDKit✔️✔️:                               true,     // onlyHeavyAtomsForRMS
        // RDKit✔️✔️:                               2,        // ETversion
        // RDKit✔️✔️:                               nullptr,  // boundsMat
        // RDKit✔️✔️:                               true,     // embedFragmentsSeparately
        // RDKit✔️✔️:                               false,    // useSmallRingTorsions
        // RDKit✔️✔️:                               true,     // useMacrocycleTorsions
        // RDKit✔️✔️:                               true,     // useMacrocycle14config
        // RDKit✔️✔️:                               0,        // timeout
        // RDKit✔️✔️:                               nullptr,  // CPCI
        // RDKit✔️✔️:                               nullptr   // callback
        Self::from_rdkit_constructor(
            0, 1, -1, true, false, 2.0, true, 1, None, 1e-3, false, true, true, true, false, 5.0,
            -1.0, true, 2, None, true, false, true, true, 0, None, None,
        )
    }

    #[must_use]
    pub fn sr_etkdg_v3() -> Self {
        // RDKit✔️✔️: const EmbedParameters srETKDGv3(0,        // maxIterations
        // RDKit✔️✔️:                                 1,        // numThreads
        // RDKit✔️✔️:                                 -1,       // randomSeed
        // RDKit✔️✔️:                                 true,     // clearConfs
        // RDKit✔️✔️:                                 false,    // useRandomCoords
        // RDKit✔️✔️:                                 2.0,      // boxSizeMult
        // RDKit✔️✔️:                                 true,     // randNegEig
        // RDKit✔️✔️:                                 1,        // numZeroFail
        // RDKit✔️✔️:                                 nullptr,  // coordMap
        // RDKit✔️✔️:                                 1e-3,     // optimizerForceTol
        // RDKit✔️✔️:                                 false,    // ignoreSmoothingFailures
        // RDKit✔️✔️:                                 true,     // enforceChirality
        // RDKit✔️✔️:                                 true,     // useExpTorsionAnglePrefs
        // RDKit✔️✔️:                                 true,     // useBasicKnowledge
        // RDKit✔️✔️:                                 false,    // verbose
        // RDKit✔️✔️:                                 5.0,      // basinThresh
        // RDKit✔️✔️:                                 -1.0,     // pruneRmsThresh
        // RDKit✔️✔️:                                 true,     // onlyHeavyAtomsForRMS
        // RDKit✔️✔️:                                 2,        // ETversion
        // RDKit✔️✔️:                                 nullptr,  // boundsMat
        // RDKit✔️✔️:                                 true,     // embedFragmentsSeparately
        // RDKit✔️✔️:                                 true,     // useSmallRingTorsions
        // RDKit✔️✔️:                                 false,    // useMacrocycleTorsions
        // RDKit✔️✔️:                                 false,    // useMacrocycle14config
        // RDKit✔️✔️:                                 0,        // timeout
        // RDKit✔️✔️:                                 nullptr,  // CPCI
        // RDKit✔️✔️:                                 nullptr   // callback
        Self::from_rdkit_constructor(
            0, 1, -1, true, false, 2.0, true, 1, None, 1e-3, false, true, true, true, false, 5.0,
            -1.0, true, 2, None, true, true, false, false, 0, None, None,
        )
    }

    pub(crate) fn update_from_json(&mut self, json: &str) -> Result<(), ConformerError> {
        // BEGIN RDKIT CPP FUNCTION DGeomHelpers::updateEmbedParametersFromJSON (EmbedderUtils.cpp:56-87)
        // RDKit✔️✔️: void updateEmbedParametersFromJSON(EmbedParameters &params,
        // RDKit✔️✔️:                                    const std::string &json) {
        // RDKit✔️✔️:   if (json.empty()) {
        // RDKit✔️✔️:     return;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   std::istringstream ss(json);
        // RDKit✔️✔️:   boost::property_tree::ptree pt;
        // RDKit✔️✔️:   boost::property_tree::read_json(ss, pt);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   EMBED_PARAMS_FIELDS(PT_OPT_GET)
        // RDKit✔️✔️:
        // RDKit✔️✔️:   std::map<int, RDGeom::Point3D> *cmap = nullptr;
        // RDKit✔️✔️:   const auto coordMap = pt.get_child_optional("coordMap");
        // RDKit✔️✔️:   if (coordMap) {
        // RDKit✔️✔️:     // NOTE: this leaks since EmbedParameters uses a naked pointer and we don't
        // RDKit✔️✔️:     // have any way to tie the lifetime of the memory we allocate here to the
        // RDKit✔️✔️:     // EmbedParameters object itself.
        // RDKit✔️✔️:     cmap = new std::map<int, RDGeom::Point3D>();
        // RDKit✔️✔️:     for (const auto &entry : *coordMap) {
        // RDKit✔️✔️:       RDGeom::Point3D pt;
        // RDKit✔️✔️:
        // RDKit✔️✔️:       auto itm = entry.second.begin();
        // RDKit✔️✔️:       pt.x = itm->second.get_value<float>();
        // RDKit✔️✔️:       ++itm;
        // RDKit✔️✔️:       pt.y = itm->second.get_value<float>();
        // RDKit✔️✔️:       ++itm;
        // RDKit✔️✔️:       pt.z = itm->second.get_value<float>();
        // RDKit✔️✔️:
        // RDKit✔️✔️:       (*cmap)[boost::lexical_cast<int>(entry.first)] = pt;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     params.coordMap = cmap;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION DGeomHelpers::updateEmbedParametersFromJSON
        if json.is_empty() {
            return Ok(());
        }
        let value: serde_json::Value = serde_json::from_str(json)
            .map_err(|err| ConformerError::InvalidEmbedParametersJson(err.to_string()))?;

        update_f64_field(&value, "basinThresh", &mut self.basin_thresh)?;
        update_f64_field(
            &value,
            "boundsMatForceScaling",
            &mut self.bounds_mat_force_scaling,
        )?;
        update_f64_field(&value, "boxSizeMult", &mut self.box_size_mult)?;
        update_bool_field(&value, "clearConfs", &mut self.clear_conformers)?;
        update_bool_field(
            &value,
            "embedFragmentsSeparately",
            &mut self.embed_fragments_separately,
        )?;
        update_bool_field(
            &value,
            "enableSequentialRandomSeeds",
            &mut self.enable_sequential_random_seeds,
        )?;
        update_bool_field(&value, "enforceChirality", &mut self.enforce_chirality)?;
        update_u32_field(&value, "ETversion", &mut self.et_version)?;
        update_bool_field(&value, "forceTransAmides", &mut self.force_trans_amides)?;
        update_bool_field(
            &value,
            "ignoreSmoothingFailures",
            &mut self.ignore_smoothing_failures,
        )?;
        update_u32_field(&value, "maxIterations", &mut self.max_iterations)?;
        update_i32_field(&value, "numThreads", &mut self.num_threads)?;
        update_u32_field(&value, "numZeroFail", &mut self.num_zero_fail)?;
        update_bool_field(
            &value,
            "onlyHeavyAtomsForRMS",
            &mut self.only_heavy_atoms_for_rms,
        )?;
        update_f64_field(&value, "optimizerForceTol", &mut self.optimizer_force_tol)?;
        update_f64_field(&value, "pruneRmsThresh", &mut self.prune_rms_thresh)?;
        update_bool_field(&value, "randNegEig", &mut self.rand_neg_eig)?;
        update_i32_field(&value, "randomSeed", &mut self.random_seed)?;
        update_bool_field(
            &value,
            "symmetrizeConjugatedTerminalGroupsForPruning",
            &mut self.symmetrize_conjugated_terminal_groups_for_pruning,
        )?;
        update_u32_field(&value, "timeout", &mut self.timeout)?;
        update_bool_field(&value, "trackFailures", &mut self.track_failures)?;
        update_bool_field(&value, "useBasicKnowledge", &mut self.use_basic_knowledge)?;
        update_bool_field(
            &value,
            "useExpTorsionAnglePrefs",
            &mut self.use_exp_torsion_angle_prefs,
        )?;
        update_bool_field(
            &value,
            "useMacrocycle14config",
            &mut self.use_macrocycle14config,
        )?;
        update_bool_field(
            &value,
            "useMacrocycleTorsions",
            &mut self.use_macrocycle_torsions,
        )?;
        update_bool_field(&value, "useRandomCoords", &mut self.use_random_coords)?;
        update_bool_field(
            &value,
            "useSmallRingTorsions",
            &mut self.use_small_ring_torsions,
        )?;
        update_bool_field(
            &value,
            "useSymmetryForPruning",
            &mut self.use_symmetry_for_pruning,
        )?;
        update_bool_field(&value, "verbose", &mut self.verbose)?;

        if let Some(coord_map) = value.get("coordMap") {
            self.coord_map = Some(parse_embed_parameters_coord_map(coord_map)?);
        }

        Ok(())
    }

    /// Returns an independent parameter value updated with source JSON fields.
    pub fn with_json(&self, json: &str) -> Result<Self, ConformerError> {
        let mut updated = self.clone();
        updated.update_from_json(json)?;
        Ok(updated)
    }

    #[must_use]
    pub fn to_json(&self) -> String {
        // BEGIN RDKIT CPP FUNCTION DGeomHelpers::embedParametersToJSON (EmbedderUtils.cpp:90-126)
        // RDKit✔️✔️: std::string embedParametersToJSON(const EmbedParameters &params) {
        // RDKit✔️✔️:   boost::property_tree::ptree pt;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   EMBED_PARAMS_FIELDS(PT_OPT_PUT)
        // RDKit✔️✔️:
        // RDKit✔️✔️:   if (params.coordMap) {
        // RDKit✔️✔️:     boost::property_tree::ptree coordMapPT;
        // RDKit✔️✔️:
        // RDKit✔️✔️:     for (const auto &kv : *params.coordMap) {
        // RDKit✔️✔️:       boost::property_tree::ptree pointPT;
        // RDKit✔️✔️:       pointPT.push_back(
        // RDKit✔️✔️:           {"", boost::property_tree::ptree(std::to_string(kv.second.x))});
        // RDKit✔️✔️:       pointPT.push_back(
        // RDKit✔️✔️:           {"", boost::property_tree::ptree(std::to_string(kv.second.y))});
        // RDKit✔️✔️:       pointPT.push_back(
        // RDKit✔️✔️:           {"", boost::property_tree::ptree(std::to_string(kv.second.z))});
        // RDKit✔️✔️:
        // RDKit✔️✔️:       coordMapPT.add_child(std::to_string(kv.first), pointPT);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:
        // RDKit✔️✔️:     pt.add_child("coordMap", coordMapPT);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   if (params.boundsMat) {
        // RDKit✔️✔️:     boost::property_tree::ptree matrixPT;
        // RDKit✔️✔️:     const unsigned int N = params.boundsMat->numCols();
        // RDKit✔️✔️:     for (unsigned i = 0; i < N; ++i) {
        // RDKit✔️✔️:       boost::property_tree::ptree rowPT;
        // RDKit✔️✔️:
        // RDKit✔️✔️:       for (unsigned j = 0; j < N; ++j) {
        // RDKit✔️✔️:         boost::property_tree::ptree v;
        // RDKit✔️✔️:         v.put("", params.boundsMat->getVal(i, j));
        // RDKit✔️✔️:         rowPT.push_back({"", v});
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:
        // RDKit✔️✔️:       matrixPT.push_back({"", rowPT});
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:
        // RDKit✔️✔️:     pt.add_child("boundsMatrix", matrixPT);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   std::ostringstream ss;
        // RDKit✔️✔️:   boost::property_tree::write_json(ss, pt, false);
        // RDKit✔️✔️:   auto str = ss.str();
        // RDKit✔️✔️:   boost::algorithm::trim(str);
        // RDKit✔️✔️:   return str;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION DGeomHelpers::embedParametersToJSON
        let mut fields = Vec::with_capacity(31);
        push_json_field(&mut fields, "basinThresh", self.basin_thresh);
        push_json_field(
            &mut fields,
            "boundsMatForceScaling",
            self.bounds_mat_force_scaling,
        );
        push_json_field(&mut fields, "boxSizeMult", self.box_size_mult);
        push_json_field(&mut fields, "clearConfs", self.clear_conformers);
        push_json_field(
            &mut fields,
            "embedFragmentsSeparately",
            self.embed_fragments_separately,
        );
        push_json_field(
            &mut fields,
            "enableSequentialRandomSeeds",
            self.enable_sequential_random_seeds,
        );
        push_json_field(&mut fields, "enforceChirality", self.enforce_chirality);
        push_json_field(&mut fields, "ETversion", self.et_version);
        push_json_field(&mut fields, "forceTransAmides", self.force_trans_amides);
        push_json_field(
            &mut fields,
            "ignoreSmoothingFailures",
            self.ignore_smoothing_failures,
        );
        push_json_field(&mut fields, "maxIterations", self.max_iterations);
        push_json_field(&mut fields, "numThreads", self.num_threads);
        push_json_field(&mut fields, "numZeroFail", self.num_zero_fail);
        push_json_field(
            &mut fields,
            "onlyHeavyAtomsForRMS",
            self.only_heavy_atoms_for_rms,
        );
        push_json_field(&mut fields, "optimizerForceTol", self.optimizer_force_tol);
        push_json_field(&mut fields, "pruneRmsThresh", self.prune_rms_thresh);
        push_json_field(&mut fields, "randNegEig", self.rand_neg_eig);
        push_json_field(&mut fields, "randomSeed", self.random_seed);
        push_json_field(
            &mut fields,
            "symmetrizeConjugatedTerminalGroupsForPruning",
            self.symmetrize_conjugated_terminal_groups_for_pruning,
        );
        push_json_field(&mut fields, "timeout", self.timeout);
        push_json_field(&mut fields, "trackFailures", self.track_failures);
        push_json_field(&mut fields, "useBasicKnowledge", self.use_basic_knowledge);
        push_json_field(
            &mut fields,
            "useExpTorsionAnglePrefs",
            self.use_exp_torsion_angle_prefs,
        );
        push_json_field(
            &mut fields,
            "useMacrocycle14config",
            self.use_macrocycle14config,
        );
        push_json_field(
            &mut fields,
            "useMacrocycleTorsions",
            self.use_macrocycle_torsions,
        );
        push_json_field(&mut fields, "useRandomCoords", self.use_random_coords);
        push_json_field(
            &mut fields,
            "useSmallRingTorsions",
            self.use_small_ring_torsions,
        );
        push_json_field(
            &mut fields,
            "useSymmetryForPruning",
            self.use_symmetry_for_pruning,
        );
        push_json_field(&mut fields, "verbose", self.verbose);

        if let Some(coord_map) = &self.coord_map {
            let entries = coord_map
                .iter()
                .map(|(atom_idx, point)| {
                    format!(
                        "\"{}\":[\"{:.6}\",\"{:.6}\",\"{:.6}\"]",
                        atom_idx, point[0], point[1], point[2]
                    )
                })
                .collect::<Vec<_>>()
                .join(",");
            fields.push(format!("\"coordMap\":{{{entries}}}"));
        }

        if let Some(bounds_mat) = &self.bounds_mat {
            let rows = (0..bounds_mat.dimension())
                .map(|i| {
                    let row = (0..bounds_mat.dimension())
                        .map(|j| {
                            format!(
                                "\"{}\"",
                                bounds_mat
                                    .get_val(i, j)
                                    .expect("bounds loop indices are in range")
                            )
                        })
                        .collect::<Vec<_>>()
                        .join(",");
                    format!("[{row}]")
                })
                .collect::<Vec<_>>()
                .join(",");
            fields.push(format!("\"boundsMatrix\":[{rows}]"));
        }

        format!("{{{}}}", fields.join(","))
    }
}

fn embed_parameters_json_field<'a>(
    value: &'a serde_json::Value,
    name: &str,
) -> Option<&'a serde_json::Value> {
    value.as_object().and_then(|object| object.get(name))
}

fn push_json_field<T: ToString>(fields: &mut Vec<String>, name: &str, value: T) {
    fields.push(format!("\"{}\":\"{}\"", name, value.to_string()));
}

fn embed_parameters_invalid_json(name: &str, expected: &str) -> ConformerError {
    ConformerError::InvalidEmbedParametersJson(format!("{name} must be {expected}"))
}

fn json_value_as_f64(name: &str, value: &serde_json::Value) -> Result<f64, ConformerError> {
    if let Some(number) = value.as_f64() {
        Ok(number)
    } else if let Some(text) = value.as_str() {
        text.parse::<f64>()
            .map_err(|_| embed_parameters_invalid_json(name, "a floating-point number"))
    } else {
        Err(embed_parameters_invalid_json(
            name,
            "a floating-point number",
        ))
    }
}

fn json_value_as_i32(name: &str, value: &serde_json::Value) -> Result<i32, ConformerError> {
    if let Some(number) = value.as_i64() {
        i32::try_from(number).map_err(|_| embed_parameters_invalid_json(name, "a 32-bit integer"))
    } else if let Some(text) = value.as_str() {
        text.parse::<i32>()
            .map_err(|_| embed_parameters_invalid_json(name, "a 32-bit integer"))
    } else {
        Err(embed_parameters_invalid_json(name, "a 32-bit integer"))
    }
}

fn json_value_as_u32(name: &str, value: &serde_json::Value) -> Result<u32, ConformerError> {
    if let Some(number) = value.as_u64() {
        u32::try_from(number)
            .map_err(|_| embed_parameters_invalid_json(name, "an unsigned 32-bit integer"))
    } else if let Some(text) = value.as_str() {
        text.parse::<u32>()
            .map_err(|_| embed_parameters_invalid_json(name, "an unsigned 32-bit integer"))
    } else {
        Err(embed_parameters_invalid_json(
            name,
            "an unsigned 32-bit integer",
        ))
    }
}

fn json_value_as_bool(name: &str, value: &serde_json::Value) -> Result<bool, ConformerError> {
    if let Some(value) = value.as_bool() {
        Ok(value)
    } else if let Some(number) = value.as_u64() {
        match number {
            0 => Ok(false),
            1 => Ok(true),
            _ => Err(embed_parameters_invalid_json(name, "a boolean")),
        }
    } else if let Some(text) = value.as_str() {
        match text {
            "0" => Ok(false),
            "1" => Ok(true),
            "false" => Ok(false),
            "true" => Ok(true),
            _ => Err(embed_parameters_invalid_json(name, "a boolean")),
        }
    } else {
        Err(embed_parameters_invalid_json(name, "a boolean"))
    }
}

fn update_f64_field(
    value: &serde_json::Value,
    name: &str,
    target: &mut f64,
) -> Result<(), ConformerError> {
    if let Some(field) = embed_parameters_json_field(value, name) {
        *target = json_value_as_f64(name, field)?;
    }
    Ok(())
}

fn update_i32_field(
    value: &serde_json::Value,
    name: &str,
    target: &mut i32,
) -> Result<(), ConformerError> {
    if let Some(field) = embed_parameters_json_field(value, name) {
        *target = json_value_as_i32(name, field)?;
    }
    Ok(())
}

fn update_u32_field(
    value: &serde_json::Value,
    name: &str,
    target: &mut u32,
) -> Result<(), ConformerError> {
    if let Some(field) = embed_parameters_json_field(value, name) {
        *target = json_value_as_u32(name, field)?;
    }
    Ok(())
}

fn update_bool_field(
    value: &serde_json::Value,
    name: &str,
    target: &mut bool,
) -> Result<(), ConformerError> {
    if let Some(field) = embed_parameters_json_field(value, name) {
        *target = json_value_as_bool(name, field)?;
    }
    Ok(())
}

fn parse_embed_parameters_coord_map(
    value: &serde_json::Value,
) -> Result<BTreeMap<i32, [f64; 3]>, ConformerError> {
    let object = value
        .as_object()
        .ok_or_else(|| embed_parameters_invalid_json("coordMap", "an object"))?;
    let mut coord_map = BTreeMap::new();
    for (key, point_value) in object {
        let atom_idx = key
            .parse::<i32>()
            .map_err(|_| embed_parameters_invalid_json("coordMap key", "a 32-bit integer"))?;
        let point = point_value
            .as_array()
            .ok_or_else(|| embed_parameters_invalid_json("coordMap value", "an array"))?;
        if point.len() < 3 {
            return Err(embed_parameters_invalid_json(
                "coordMap value",
                "an array with at least three coordinates",
            ));
        }
        coord_map.insert(
            atom_idx,
            [
                json_value_as_f64("coordMap x", &point[0])?,
                json_value_as_f64("coordMap y", &point[1])?,
                json_value_as_f64("coordMap z", &point[2])?,
            ],
        );
    }
    Ok(coord_map)
}

// RDKit✔️✔️:   EmbedParameters(
// RDKit✔️✔️:       unsigned int maxIterations, int numThreads, int randomSeed,
// RDKit✔️✔️:       bool clearConfs, bool useRandomCoords, double boxSizeMult,
// RDKit✔️✔️:       bool randNegEig, unsigned int numZeroFail,
// RDKit✔️✔️:       const std::map<int, RDGeom::Point3D> *coordMap, double optimizerForceTol,
// RDKit✔️✔️:       bool ignoreSmoothingFailures, bool enforceChirality,
// RDKit✔️✔️:       bool useExpTorsionAnglePrefs, bool useBasicKnowledge, bool verbose,
// RDKit✔️✔️:       double basinThresh, double pruneRmsThresh, bool onlyHeavyAtomsForRMS,
// RDKit✔️✔️:       unsigned int ETversion = 2,
// RDKit✔️✔️:       const DistGeom::BoundsMatrix *boundsMat = nullptr,
// RDKit✔️✔️:       bool embedFragmentsSeparately = true, bool useSmallRingTorsions = false,
// RDKit✔️✔️:       bool useMacrocycleTorsions = false, bool useMacrocycle14config = false,
// RDKit✔️✔️:       unsigned int timeout = 0,
// RDKit✔️✔️:       std::shared_ptr<std::map<std::pair<unsigned int, unsigned int>, double>>
// RDKit✔️✔️:           CPCI = nullptr,
// RDKit✔️✔️:       void (*callback)(unsigned int) = nullptr)
// RDKit✔️✔️:       : maxIterations(maxIterations),
// RDKit✔️✔️:         numThreads(numThreads),
// RDKit✔️✔️:         randomSeed(randomSeed),
// RDKit✔️✔️:         clearConfs(clearConfs),
// RDKit✔️✔️:         useRandomCoords(useRandomCoords),
// RDKit✔️✔️:         boxSizeMult(boxSizeMult),
// RDKit✔️✔️:         randNegEig(randNegEig),
// RDKit✔️✔️:         numZeroFail(numZeroFail),
// RDKit✔️✔️:         coordMap(coordMap),
// RDKit✔️✔️:         optimizerForceTol(optimizerForceTol),
// RDKit✔️✔️:         ignoreSmoothingFailures(ignoreSmoothingFailures),
// RDKit✔️✔️:         enforceChirality(enforceChirality),
// RDKit✔️✔️:         useExpTorsionAnglePrefs(useExpTorsionAnglePrefs),
// RDKit✔️✔️:         useBasicKnowledge(useBasicKnowledge),
// RDKit✔️✔️:         verbose(verbose),
// RDKit✔️✔️:         basinThresh(basinThresh),
// RDKit✔️✔️:         pruneRmsThresh(pruneRmsThresh),
// RDKit✔️✔️:         onlyHeavyAtomsForRMS(onlyHeavyAtomsForRMS),
// RDKit✔️✔️:         ETversion(ETversion),
// RDKit✔️✔️:         boundsMat(boundsMat),
// RDKit✔️✔️:         embedFragmentsSeparately(embedFragmentsSeparately),
// RDKit✔️✔️:         useSmallRingTorsions(useSmallRingTorsions),
// RDKit✔️✔️:         useMacrocycleTorsions(useMacrocycleTorsions),
// RDKit✔️✔️:         useMacrocycle14config(useMacrocycle14config),
// RDKit✔️✔️:         timeout(timeout),
// RDKit✔️✔️:         CPCI(std::move(CPCI)),
// RDKit✔️✔️:         callback(callback) {}
// RDKit✔️✔️: };
// END RDKIT CPP STRUCT DGeomHelpers::EmbedParameters

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn embed_failure_causes_match_rdkit_ordinals_and_names() {
        let expected = [
            (EmbedFailureCause::InitialCoords, 0, "INITIAL_COORDS"),
            (
                EmbedFailureCause::FirstMinimization,
                1,
                "FIRST_MINIMIZATION",
            ),
            (
                EmbedFailureCause::CheckTetrahedralCenters,
                2,
                "CHECK_TETRAHEDRAL_CENTERS",
            ),
            (
                EmbedFailureCause::CheckChiralCenters,
                3,
                "CHECK_CHIRAL_CENTERS",
            ),
            (
                EmbedFailureCause::MinimizeFourthDimension,
                4,
                "MINIMIZE_FOURTH_DIMENSION",
            ),
            (EmbedFailureCause::EtkMinimization, 5, "ETK_MINIMIZATION"),
            (
                EmbedFailureCause::FinalChiralBounds,
                6,
                "FINAL_CHIRAL_BOUNDS",
            ),
            (
                EmbedFailureCause::FinalCenterInVolume,
                7,
                "FINAL_CENTER_IN_VOLUME",
            ),
            (EmbedFailureCause::LinearDoubleBond, 8, "LINEAR_DOUBLE_BOND"),
            (
                EmbedFailureCause::BadDoubleBondStereo,
                9,
                "BAD_DOUBLE_BOND_STEREO",
            ),
            (
                EmbedFailureCause::CheckChiralCenters2,
                10,
                "CHECK_CHIRAL_CENTERS2",
            ),
            (EmbedFailureCause::ExceededTimeout, 11, "EXCEEDED_TIMEOUT"),
            (EmbedFailureCause::EndOfEnum, 12, "END_OF_ENUM"),
        ];

        assert_eq!(EmbedFailureCause::ALL.len(), expected.len());
        for (idx, (cause, ordinal, name)) in expected.iter().copied().enumerate() {
            assert_eq!(EmbedFailureCause::ALL[idx], cause);
            assert_eq!(cause.rdkit_ordinal(), ordinal);
            assert_eq!(EmbedFailureCause::from_rdkit_ordinal(ordinal), Some(cause));
            assert_eq!(cause.rdkit_name(), name);
        }
        assert_eq!(EmbedFailureCause::from_rdkit_ordinal(13), None);
    }

    #[test]
    fn embed_parameters_defaults_match_rdkit_constructor_defaults() {
        let params = EmbedParams::default();

        assert_eq!(params.max_iterations, 0);
        assert_eq!(params.num_threads, 1);
        assert_eq!(params.random_seed, -1);
        assert!(params.clear_conformers);
        assert!(!params.use_random_coords);
        assert_eq!(params.box_size_mult, 2.0);
        assert!(params.rand_neg_eig);
        assert_eq!(params.num_zero_fail, 1);
        assert!(params.coord_map.is_none());
        assert_eq!(params.optimizer_force_tol, 1e-3);
        assert!(!params.ignore_smoothing_failures);
        assert!(params.enforce_chirality);
        assert!(!params.use_exp_torsion_angle_prefs);
        assert!(!params.use_basic_knowledge);
        assert!(!params.verbose);
        assert_eq!(params.basin_thresh, 5.0);
        assert_eq!(params.prune_rms_thresh, -1.0);
        assert!(params.only_heavy_atoms_for_rms);
        assert_eq!(params.et_version, 2);
        assert!(params.bounds_mat.is_none());
        assert!(params.embed_fragments_separately);
        assert!(!params.use_small_ring_torsions);
        assert!(!params.use_macrocycle_torsions);
        assert!(!params.use_macrocycle14config);
        assert_eq!(params.timeout, 0);
        assert!(params.cpci.is_none());
        assert!(params.callback.is_none());
        assert!(params.force_trans_amides);
        assert!(params.use_symmetry_for_pruning);
        assert_eq!(params.bounds_mat_force_scaling, 1.0);
        assert!(!params.track_failures);
        assert!(params.failures.is_empty());
        assert!(!params.enable_sequential_random_seeds);
        assert!(params.symmetrize_conjugated_terminal_groups_for_pruning);
    }

    #[test]
    fn embed_parameters_new_matches_default_constructor() {
        let from_new = EmbedParams::new();
        let from_default = EmbedParams::default();

        assert_eq!(from_new.max_iterations, from_default.max_iterations);
        assert_eq!(from_new.num_threads, from_default.num_threads);
        assert_eq!(from_new.random_seed, from_default.random_seed);
        assert_eq!(from_new.clear_conformers, from_default.clear_conformers);
        assert_eq!(from_new.use_random_coords, from_default.use_random_coords);
        assert_eq!(from_new.box_size_mult, from_default.box_size_mult);
        assert_eq!(from_new.rand_neg_eig, from_default.rand_neg_eig);
        assert_eq!(from_new.num_zero_fail, from_default.num_zero_fail);
        assert_eq!(
            from_new.optimizer_force_tol,
            from_default.optimizer_force_tol
        );
        assert_eq!(
            from_new.ignore_smoothing_failures,
            from_default.ignore_smoothing_failures
        );
        assert_eq!(from_new.enforce_chirality, from_default.enforce_chirality);
        assert_eq!(
            from_new.use_exp_torsion_angle_prefs,
            from_default.use_exp_torsion_angle_prefs
        );
        assert_eq!(
            from_new.use_basic_knowledge,
            from_default.use_basic_knowledge
        );
        assert_eq!(from_new.verbose, from_default.verbose);
        assert_eq!(from_new.basin_thresh, from_default.basin_thresh);
        assert_eq!(from_new.prune_rms_thresh, from_default.prune_rms_thresh);
        assert_eq!(
            from_new.only_heavy_atoms_for_rms,
            from_default.only_heavy_atoms_for_rms
        );
        assert_eq!(from_new.et_version, from_default.et_version);
        assert_eq!(
            from_new.embed_fragments_separately,
            from_default.embed_fragments_separately
        );
        assert_eq!(
            from_new.use_small_ring_torsions,
            from_default.use_small_ring_torsions
        );
        assert_eq!(
            from_new.use_macrocycle_torsions,
            from_default.use_macrocycle_torsions
        );
        assert_eq!(
            from_new.use_macrocycle14config,
            from_default.use_macrocycle14config
        );
        assert_eq!(from_new.timeout, from_default.timeout);
        assert_eq!(from_new.force_trans_amides, from_default.force_trans_amides);
        assert_eq!(
            from_new.use_symmetry_for_pruning,
            from_default.use_symmetry_for_pruning
        );
        assert_eq!(
            from_new.bounds_mat_force_scaling,
            from_default.bounds_mat_force_scaling
        );
        assert_eq!(from_new.track_failures, from_default.track_failures);
        assert_eq!(
            from_new.enable_sequential_random_seeds,
            from_default.enable_sequential_random_seeds
        );
        assert_eq!(
            from_new.symmetrize_conjugated_terminal_groups_for_pruning,
            from_default.symmetrize_conjugated_terminal_groups_for_pruning
        );
    }

    #[allow(clippy::too_many_arguments)]
    fn assert_embed_parameters_preset(
        params: &EmbedParams,
        max_iterations: u32,
        num_threads: i32,
        random_seed: i32,
        clear_conformers: bool,
        use_random_coords: bool,
        box_size_mult: f64,
        rand_neg_eig: bool,
        num_zero_fail: u32,
        optimizer_force_tol: f64,
        ignore_smoothing_failures: bool,
        enforce_chirality: bool,
        use_exp_torsion_angle_prefs: bool,
        use_basic_knowledge: bool,
        verbose: bool,
        basin_thresh: f64,
        prune_rms_thresh: f64,
        only_heavy_atoms_for_rms: bool,
        et_version: u32,
        embed_fragments_separately: bool,
        use_small_ring_torsions: bool,
        use_macrocycle_torsions: bool,
        use_macrocycle14config: bool,
        timeout: u32,
    ) {
        assert_eq!(params.max_iterations, max_iterations);
        assert_eq!(params.num_threads, num_threads);
        assert_eq!(params.random_seed, random_seed);
        assert_eq!(params.clear_conformers, clear_conformers);
        assert_eq!(params.use_random_coords, use_random_coords);
        assert_eq!(params.box_size_mult, box_size_mult);
        assert_eq!(params.rand_neg_eig, rand_neg_eig);
        assert_eq!(params.num_zero_fail, num_zero_fail);
        assert!(params.coord_map.is_none());
        assert_eq!(params.optimizer_force_tol, optimizer_force_tol);
        assert_eq!(params.ignore_smoothing_failures, ignore_smoothing_failures);
        assert_eq!(params.enforce_chirality, enforce_chirality);
        assert_eq!(
            params.use_exp_torsion_angle_prefs,
            use_exp_torsion_angle_prefs
        );
        assert_eq!(params.use_basic_knowledge, use_basic_knowledge);
        assert_eq!(params.verbose, verbose);
        assert_eq!(params.basin_thresh, basin_thresh);
        assert_eq!(params.prune_rms_thresh, prune_rms_thresh);
        assert_eq!(params.only_heavy_atoms_for_rms, only_heavy_atoms_for_rms);
        assert_eq!(params.et_version, et_version);
        assert!(params.bounds_mat.is_none());
        assert_eq!(
            params.embed_fragments_separately,
            embed_fragments_separately
        );
        assert_eq!(params.use_small_ring_torsions, use_small_ring_torsions);
        assert_eq!(params.use_macrocycle_torsions, use_macrocycle_torsions);
        assert_eq!(params.use_macrocycle14config, use_macrocycle14config);
        assert_eq!(params.timeout, timeout);
        assert!(params.cpci.is_none());
        assert!(params.callback.is_none());
        assert!(params.force_trans_amides);
        assert!(params.use_symmetry_for_pruning);
        assert_eq!(params.bounds_mat_force_scaling, 1.0);
        assert!(!params.track_failures);
        assert!(params.failures.is_empty());
        assert!(!params.enable_sequential_random_seeds);
        assert!(params.symmetrize_conjugated_terminal_groups_for_pruning);
    }

    #[test]
    fn embed_parameter_presets_match_rdkit_global_parameters() {
        assert_embed_parameters_preset(
            &EmbedParams::kdg(),
            0,
            1,
            -1,
            true,
            false,
            2.0,
            true,
            1,
            1e-3,
            false,
            true,
            false,
            true,
            false,
            5.0,
            -1.0,
            true,
            1,
            true,
            false,
            false,
            false,
            0,
        );
        assert_embed_parameters_preset(
            &EmbedParams::etdg(),
            0,
            1,
            -1,
            true,
            false,
            2.0,
            true,
            1,
            1e-3,
            false,
            false,
            true,
            false,
            false,
            5.0,
            -1.0,
            true,
            1,
            true,
            false,
            false,
            false,
            0,
        );
        assert_embed_parameters_preset(
            &EmbedParams::etdg_v2(),
            0,
            1,
            -1,
            true,
            false,
            2.0,
            true,
            1,
            1e-3,
            false,
            false,
            true,
            false,
            false,
            5.0,
            -1.0,
            true,
            2,
            true,
            false,
            false,
            false,
            0,
        );
        assert_embed_parameters_preset(
            &EmbedParams::etkdg(),
            0,
            1,
            -1,
            true,
            false,
            2.0,
            true,
            1,
            1e-3,
            false,
            true,
            true,
            true,
            false,
            5.0,
            -1.0,
            true,
            1,
            true,
            false,
            false,
            false,
            0,
        );
        assert_embed_parameters_preset(
            &EmbedParams::etkdg_v2(),
            0,
            1,
            -1,
            true,
            false,
            2.0,
            true,
            1,
            1e-3,
            false,
            true,
            true,
            true,
            false,
            5.0,
            -1.0,
            true,
            2,
            true,
            false,
            false,
            false,
            0,
        );
        assert_embed_parameters_preset(
            &EmbedParams::etkdg_v3(),
            0,
            1,
            -1,
            true,
            false,
            2.0,
            true,
            1,
            1e-3,
            false,
            true,
            true,
            true,
            false,
            5.0,
            -1.0,
            true,
            2,
            true,
            false,
            true,
            true,
            0,
        );
        assert_embed_parameters_preset(
            &EmbedParams::sr_etkdg_v3(),
            0,
            1,
            -1,
            true,
            false,
            2.0,
            true,
            1,
            1e-3,
            false,
            true,
            true,
            true,
            false,
            5.0,
            -1.0,
            true,
            2,
            true,
            true,
            false,
            false,
            0,
        );
    }

    #[test]
    fn update_embed_parameters_from_json_empty_string_is_noop() {
        let mut params = EmbedParams::etkdg_v3();
        let before = params.clone();

        params.update_from_json("").unwrap();

        assert_eq!(params.max_iterations, before.max_iterations);
        assert_eq!(params.et_version, before.et_version);
        assert_eq!(
            params.use_macrocycle_torsions,
            before.use_macrocycle_torsions
        );
        assert_eq!(params.coord_map, before.coord_map);
    }

    #[test]
    fn update_embed_parameters_from_json_updates_only_present_scalar_fields() {
        let mut params = EmbedParams::default();

        params
            .update_from_json(
                r#"{
              "maxIterations": 23,
              "randomSeed": 17,
              "useRandomCoords": true,
              "optimizerForceTol": 0.25,
              "useMacrocycleTorsions": true,
              "unknownIgnored": false
            }"#,
            )
            .unwrap();

        assert_eq!(params.max_iterations, 23);
        assert_eq!(params.random_seed, 17);
        assert!(params.use_random_coords);
        assert_eq!(params.optimizer_force_tol, 0.25);
        assert!(params.use_macrocycle_torsions);
        assert_eq!(params.num_threads, 1);
        assert!(params.clear_conformers);
        assert!(!params.use_macrocycle14config);
    }

    #[test]
    fn update_embed_parameters_from_json_updates_all_rdkit_macro_fields() {
        let mut params = EmbedParams::default();

        params
            .update_from_json(
                r#"{
              "basinThresh": 1.25,
              "boundsMatForceScaling": 2.5,
              "boxSizeMult": 3.5,
              "clearConfs": false,
              "embedFragmentsSeparately": false,
              "enableSequentialRandomSeeds": true,
              "enforceChirality": false,
              "ETversion": 3,
              "forceTransAmides": false,
              "ignoreSmoothingFailures": true,
              "maxIterations": 101,
              "numThreads": 4,
              "numZeroFail": 5,
              "onlyHeavyAtomsForRMS": false,
              "optimizerForceTol": 0.125,
              "pruneRmsThresh": 0.75,
              "randNegEig": false,
              "randomSeed": 99,
              "symmetrizeConjugatedTerminalGroupsForPruning": false,
              "timeout": 44,
              "trackFailures": true,
              "useBasicKnowledge": true,
              "useExpTorsionAnglePrefs": true,
              "useMacrocycle14config": true,
              "useMacrocycleTorsions": true,
              "useRandomCoords": true,
              "useSmallRingTorsions": true,
              "useSymmetryForPruning": false,
              "verbose": true
            }"#,
            )
            .unwrap();

        assert_eq!(params.basin_thresh, 1.25);
        assert_eq!(params.bounds_mat_force_scaling, 2.5);
        assert_eq!(params.box_size_mult, 3.5);
        assert!(!params.clear_conformers);
        assert!(!params.embed_fragments_separately);
        assert!(params.enable_sequential_random_seeds);
        assert!(!params.enforce_chirality);
        assert_eq!(params.et_version, 3);
        assert!(!params.force_trans_amides);
        assert!(params.ignore_smoothing_failures);
        assert_eq!(params.max_iterations, 101);
        assert_eq!(params.num_threads, 4);
        assert_eq!(params.num_zero_fail, 5);
        assert!(!params.only_heavy_atoms_for_rms);
        assert_eq!(params.optimizer_force_tol, 0.125);
        assert_eq!(params.prune_rms_thresh, 0.75);
        assert!(!params.rand_neg_eig);
        assert_eq!(params.random_seed, 99);
        assert!(!params.symmetrize_conjugated_terminal_groups_for_pruning);
        assert_eq!(params.timeout, 44);
        assert!(params.track_failures);
        assert!(params.use_basic_knowledge);
        assert!(params.use_exp_torsion_angle_prefs);
        assert!(params.use_macrocycle14config);
        assert!(params.use_macrocycle_torsions);
        assert!(params.use_random_coords);
        assert!(params.use_small_ring_torsions);
        assert!(!params.use_symmetry_for_pruning);
        assert!(params.verbose);
    }

    #[test]
    fn update_embed_parameters_from_json_updates_coord_map() {
        let mut params = EmbedParams::default();

        params
            .update_from_json(
                r#"{
              "coordMap": {
                "2": [1.0, 2.0, 3.0],
                "5": ["4.5", "5.5", "6.5"]
              }
            }"#,
            )
            .unwrap();

        let coord_map = params.coord_map.as_ref().unwrap();
        assert_eq!(coord_map.len(), 2);
        assert_eq!(coord_map[&2], [1.0, 2.0, 3.0]);
        assert_eq!(coord_map[&5], [4.5, 5.5, 6.5]);
    }

    #[test]
    fn update_embed_parameters_from_json_rejects_invalid_json_and_field_types() {
        let mut params = EmbedParams::default();

        assert!(params.update_from_json("{").is_err());
        assert!(params.update_from_json(r#"{"maxIterations": -1}"#).is_err());
        assert!(
            params
                .update_from_json(r#"{"coordMap": {"x": [1.0, 2.0, 3.0]}}"#)
                .is_err()
        );
        assert!(
            params
                .update_from_json(r#"{"coordMap": {"1": [1.0, 2.0]}}"#)
                .is_err()
        );
    }

    #[test]
    fn embed_parameters_to_json_matches_rdkit_without_maps() {
        let json = EmbedParams::kdg().to_json();
        let expected = r#"{"basinThresh":"5","boundsMatForceScaling":"1","boxSizeMult":"2","clearConfs":"true","embedFragmentsSeparately":"true","enableSequentialRandomSeeds":"false","enforceChirality":"true","ETversion":"1","forceTransAmides":"true","ignoreSmoothingFailures":"false","maxIterations":"0","numThreads":"1","numZeroFail":"1","onlyHeavyAtomsForRMS":"true","optimizerForceTol":"0.001","pruneRmsThresh":"-1","randNegEig":"true","randomSeed":"-1","symmetrizeConjugatedTerminalGroupsForPruning":"true","timeout":"0","trackFailures":"false","useBasicKnowledge":"true","useExpTorsionAnglePrefs":"false","useMacrocycle14config":"false","useMacrocycleTorsions":"false","useRandomCoords":"false","useSmallRingTorsions":"false","useSymmetryForPruning":"true","verbose":"false"}"#;

        assert_eq!(json, expected);
    }

    #[test]
    fn embed_parameters_to_json_matches_rdkit_coord_map_shape() {
        let mut params = EmbedParams::kdg();
        params
            .coord_map
            .get_or_insert_with(BTreeMap::new)
            .insert(3, [1.1, 2.2, 3.3]);

        let json = params.to_json();
        let expected = r#"{"basinThresh":"5","boundsMatForceScaling":"1","boxSizeMult":"2","clearConfs":"true","embedFragmentsSeparately":"true","enableSequentialRandomSeeds":"false","enforceChirality":"true","ETversion":"1","forceTransAmides":"true","ignoreSmoothingFailures":"false","maxIterations":"0","numThreads":"1","numZeroFail":"1","onlyHeavyAtomsForRMS":"true","optimizerForceTol":"0.001","pruneRmsThresh":"-1","randNegEig":"true","randomSeed":"-1","symmetrizeConjugatedTerminalGroupsForPruning":"true","timeout":"0","trackFailures":"false","useBasicKnowledge":"true","useExpTorsionAnglePrefs":"false","useMacrocycle14config":"false","useMacrocycleTorsions":"false","useRandomCoords":"false","useSmallRingTorsions":"false","useSymmetryForPruning":"true","verbose":"false","coordMap":{"3":["1.100000","2.200000","3.300000"]}}"#;

        assert_eq!(json, expected);
    }

    #[test]
    fn embed_parameters_to_json_includes_bounds_matrix_values() {
        let mut params = EmbedParams::kdg();
        let mut bounds = BoundsMatrix::new(2).unwrap();
        bounds.set_val(0, 0, 0.0).unwrap();
        bounds.set_val(0, 1, 2.5).unwrap();
        bounds.set_val(1, 0, 1.25).unwrap();
        bounds.set_val(1, 1, 0.0).unwrap();
        params.bounds_mat = Some(Arc::new(bounds));

        let json = params.to_json();

        assert!(json.ends_with(r#","boundsMatrix":[["0","2.5"],["1.25","0"]]}"#));
    }

    #[test]
    fn embed_parameters_to_json_round_trips_through_update_from_json() {
        let params = EmbedParams::etkdg_v3();
        let json = params.to_json();
        let mut round_trip = EmbedParams::default();

        round_trip.update_from_json(&json).unwrap();

        assert_eq!(round_trip.to_json(), json);
    }
}

impl EmbedParams {
    pub fn dg() -> Self {
        Self::default()
    }
    pub fn max_iterations(&self) -> u32 {
        self.max_iterations
    }
    pub fn num_threads(&self) -> i32 {
        self.num_threads
    }
    pub fn random_seed(&self) -> i32 {
        self.random_seed
    }
    pub fn clear_conformers(&self) -> bool {
        self.clear_conformers
    }
    pub fn use_random_coords(&self) -> bool {
        self.use_random_coords
    }
    pub fn box_size_mult(&self) -> f64 {
        self.box_size_mult
    }
    pub fn rand_neg_eig(&self) -> bool {
        self.rand_neg_eig
    }
    pub fn num_zero_fail(&self) -> u32 {
        self.num_zero_fail
    }
    pub fn coord_map(&self) -> &Option<BTreeMap<i32, [f64; 3]>> {
        &self.coord_map
    }
    pub fn optimizer_force_tol(&self) -> f64 {
        self.optimizer_force_tol
    }
    pub fn ignore_smoothing_failures(&self) -> bool {
        self.ignore_smoothing_failures
    }
    pub fn enforce_chirality(&self) -> bool {
        self.enforce_chirality
    }
    pub fn use_exp_torsion_angle_prefs(&self) -> bool {
        self.use_exp_torsion_angle_prefs
    }
    pub fn use_basic_knowledge(&self) -> bool {
        self.use_basic_knowledge
    }
    pub fn verbose(&self) -> bool {
        self.verbose
    }
    pub fn basin_thresh(&self) -> f64 {
        self.basin_thresh
    }
    pub fn prune_rms_thresh(&self) -> f64 {
        self.prune_rms_thresh
    }
    pub fn only_heavy_atoms_for_rms(&self) -> bool {
        self.only_heavy_atoms_for_rms
    }
    pub fn et_version(&self) -> u32 {
        self.et_version
    }
    pub fn embed_fragments_separately(&self) -> bool {
        self.embed_fragments_separately
    }
    pub fn use_small_ring_torsions(&self) -> bool {
        self.use_small_ring_torsions
    }
    pub fn use_macrocycle_torsions(&self) -> bool {
        self.use_macrocycle_torsions
    }
    pub fn use_macrocycle14config(&self) -> bool {
        self.use_macrocycle14config
    }
    pub fn timeout(&self) -> u32 {
        self.timeout
    }
    pub fn cpci(&self) -> &Option<BTreeMap<(u32, u32), f64>> {
        &self.cpci
    }
    pub fn force_trans_amides(&self) -> bool {
        self.force_trans_amides
    }
    pub fn use_symmetry_for_pruning(&self) -> bool {
        self.use_symmetry_for_pruning
    }
    pub fn bounds_mat_force_scaling(&self) -> f64 {
        self.bounds_mat_force_scaling
    }
    pub fn track_failures(&self) -> bool {
        self.track_failures
    }
    pub fn failures(&self) -> &Vec<u32> {
        &self.failures
    }
    pub fn enable_sequential_random_seeds(&self) -> bool {
        self.enable_sequential_random_seeds
    }
    pub fn symmetrize_conjugated_terminal_groups_for_pruning(&self) -> bool {
        self.symmetrize_conjugated_terminal_groups_for_pruning
    }
}
