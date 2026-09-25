//! Detached structural metadata values.

/// Classification of the software entry recorded in structural metadata.
///
/// The explicit representation and declaration order preserve Gemmi's
/// `SoftwareItem::Classification` integer values.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum BioSoftwareClassification {
    DataCollection,
    DataExtraction,
    DataProcessing,
    DataReduction,
    DataScaling,
    ModelBuilding,
    Phasing,
    Refinement,
    Unspecified,
}

/// One source-ordered software record from structural metadata.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BioSoftwareItem {
    pub name: String,
    pub version: String,
    pub date: String,
    pub description: String,
    pub contact_author: String,
    pub contact_author_email: String,
    pub classification: BioSoftwareClassification,
}

/// Reflection statistics retained from diffraction metadata.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BioReflectionsInfo {
    pub resolution_high: f64,
    pub resolution_low: f64,
    pub completeness: f64,
    pub redundancy: f64,
    pub r_merge: f64,
    pub r_sym: f64,
    pub mean_i_over_sigma: f64,
}

impl Default for BioReflectionsInfo {
    fn default() -> Self {
        // Gemmi❗✔️: struct ReflectionsInfo {
        // Gemmi❗✔️:   double resolution_high = NAN; // _reflns.d_resolution_high
        // Gemmi❗✔️:                                 // (or _reflns_shell.d_res_high)
        // Gemmi❗✔️:   double resolution_low = NAN;  // _reflns.d_resolution_low
        // Gemmi❗✔️:   double completeness = NAN;    // _reflns.percent_possible_obs
        // Gemmi❗✔️:   double redundancy = NAN;      // _reflns.pdbx_redundancy
        // Gemmi❗✔️:   double r_merge = NAN;         // _reflns.pdbx_Rmerge_I_obs
        // Gemmi❗✔️:   double r_sym = NAN;           // _reflns.pdbx_Rsym_value
        // Gemmi❗✔️:   double mean_I_over_sigma = NAN; // _reflns.pdbx_netI_over_sigmaI
        // Gemmi❗✔️: };
        // Behavior review: each source-missing reflection statistic is retained
        // as NaN rather than a numeric default or an absent optional field.
        // Complexity review: construction assigns seven scalar fields in O(1)
        // time and allocates no memory, matching the source value initialization.
        Self {
            resolution_high: f64::NAN,
            resolution_low: f64::NAN,
            completeness: f64::NAN,
            redundancy: f64::NAN,
            r_merge: f64::NAN,
            r_sym: f64::NAN,
            mean_i_over_sigma: f64::NAN,
        }
    }
}

/// Refinement statistics shared by overall refinement and per-bin values.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BioBasicRefinementInfo {
    pub resolution_high: f64,
    pub resolution_low: f64,
    pub completeness: f64,
    pub reflection_count: i32,
    pub work_set_count: i32,
    pub rfree_set_count: i32,
    pub r_all: f64,
    pub r_work: f64,
    pub r_free: f64,
    pub cc_fo_fc_work: f64,
    pub cc_fo_fc_free: f64,
    pub fsc_work: f64,
    pub fsc_free: f64,
    pub cc_intensity_work: f64,
    pub cc_intensity_free: f64,
}

impl Default for BioBasicRefinementInfo {
    fn default() -> Self {
        // Gemmi❗✔️: // BasicRefinementInfo is used for both total and per-bin statistics.
        // Gemmi❗✔️: // For per-bin data, each values corresponds to one _refine_ls_shell.* tag.
        // Gemmi❗✔️: struct BasicRefinementInfo {
        // Gemmi❗✔️:   double resolution_high = NAN; // _refine.ls_d_res_high,         _refine_ls_shell.d_res_high
        // Gemmi❗✔️:   double resolution_low = NAN;  // _refine.ls_d_res_low,          _refine_ls_shell.d_res_low
        // Gemmi❗✔️:   double completeness = NAN;    // _refine.ls_percent_reflns_obs, _refine_ls_shell.percent...
        // Gemmi❗✔️:   int reflection_count = -1;    // _refine.ls_number_reflns_obs, _refine_ls_shell.number...
        // Gemmi❗✔️:   int work_set_count = -1;      // _refine.ls_number_reflns_R_work, _refine_ls_shell.number...
        // Gemmi❗✔️:   int rfree_set_count = -1;     // _refine.ls_number_reflns_R_free, _refine_ls_shell.number...
        // Gemmi❗✔️:   double r_all = NAN;           // _refine.ls_R_factor_obs,       _refine_ls_shell.R_factor_obs
        // Gemmi❗✔️:   double r_work = NAN;          // _refine.ls_R_factor_R_work,    _refine_ls_shell.R_factor_R_work
        // Gemmi❗✔️:   double r_free = NAN;          // _refine.ls_R_factor_R_free,    _refine_ls_shell.R_factor_R_free
        // Gemmi❗✔️:   double cc_fo_fc_work = NAN;   // _refine.correlation_coeff_Fo_to_Fc, _refine_ls_shell.corr...
        // Gemmi❗✔️:   double cc_fo_fc_free = NAN;   // _refine.correlation_coeff_Fo_to_Fc_free, _refine_ls_shell.c...
        // Gemmi❗✔️:   double fsc_work = NAN;        // _refine.pdbx_average_fsc_work, _refine_ls_shell.pdbx_fsc_work
        // Gemmi❗✔️:   double fsc_free = NAN;        // _refine.pdbx_average_fsc_free, _refine_ls_shell.pdbx_fsc_free
        // Gemmi❗✔️:   double cc_intensity_work = NAN;  // _refine.correlation_coeff_I_to_Fcsqd_work, ...
        // Gemmi❗✔️:   double cc_intensity_free = NAN;  // _refine.correlation_coeff_I_to_Fcsqd_free, ...
        // Gemmi❗✔️: };
        // Behavior review: retain every source field independently; the three
        // counts default to -1 and all twelve floating statistics to NaN. This
        // value is suitable for both total and per-bin use without conversion.
        // Complexity review: fifteen scalar assignments are fixed O(1) work,
        // with no allocation.
        Self {
            resolution_high: f64::NAN,
            resolution_low: f64::NAN,
            completeness: f64::NAN,
            reflection_count: -1,
            work_set_count: -1,
            rfree_set_count: -1,
            r_all: f64::NAN,
            r_work: f64::NAN,
            r_free: f64::NAN,
            cc_fo_fc_work: f64::NAN,
            cc_fo_fc_free: f64::NAN,
            fsc_work: f64::NAN,
            fsc_free: f64::NAN,
            cc_intensity_work: f64::NAN,
            cc_intensity_free: f64::NAN,
        }
    }
}

/// One named restraint statistic in refinement metadata.
#[derive(Debug, Clone, PartialEq)]
pub struct BioRefinementRestraint {
    pub name: String,
    pub count: i32,
    pub weight: f64,
    pub function: String,
    pub dev_ideal: f64,
}

impl Default for BioRefinementRestraint {
    fn default() -> Self {
        // Gemmi❗✔️:   struct Restr {
        // Gemmi❗✔️:     std::string name;
        // Gemmi❗✔️:     int count = -1;
        // Gemmi❗✔️:     double weight = NAN;
        // Gemmi❗✔️:     std::string function;
        // Gemmi❗✔️:     double dev_ideal = NAN;
        // Gemmi❗✔️:     Restr() = default;
        // Gemmi❗✔️:   };
        // Behavior review: empty source strings remain empty; the source count
        // sentinel remains -1 and both absent numeric statistics remain NaN.
        // Complexity review: five field initializations are constant-time and
        // empty String values allocate no heap storage.
        Self {
            name: String::new(),
            count: -1,
            weight: f64::NAN,
            function: String::new(),
            dev_ideal: f64::NAN,
        }
    }
}

impl BioRefinementRestraint {
    /// Creates a restraint statistic with its source-defined name constructor.
    pub fn new(name: impl Into<String>) -> Self {
        // Gemmi❗✔️:     explicit Restr(const std::string& name_) : name(name_) {}
        // Behavior review: the supplied name is retained while count, weight,
        // function, and ideal deviation keep their member-initializer defaults.
        // Complexity review: borrowed strings are copied once into the owned
        // field; an owned String is moved without cloning.
        Self {
            name: name.into(),
            ..Self::default()
        }
    }
}

/// Experimental method and its associated reflection statistics.
#[derive(Debug, Clone, PartialEq)]
pub struct BioExperimentInfo {
    pub method: String,
    pub number_of_crystals: i32,
    pub unique_reflections: i32,
    pub reflections: BioReflectionsInfo,
    pub b_wilson: f64,
    pub shells: Vec<BioReflectionsInfo>,
    pub diffraction_ids: Vec<String>,
}

impl Default for BioExperimentInfo {
    fn default() -> Self {
        // Gemmi❗✔️: struct ExperimentInfo {
        // Gemmi❗✔️:   std::string method;             // _exptl.method
        // Gemmi❗✔️:   int number_of_crystals = -1;    // _exptl.crystals_number
        // Gemmi❗✔️:   int unique_reflections = -1;    // _reflns.number_obs
        // Gemmi❗✔️:   ReflectionsInfo reflections;
        // Gemmi❗✔️:   double b_wilson = NAN;          // _reflns.B_iso_Wilson_estimate
        // Gemmi❗✔️:   std::vector<ReflectionsInfo> shells;
        // Gemmi❗✔️:   std::vector<std::string> diffraction_ids;
        // Gemmi❗✔️: };
        // Behavior review: the method and vectors default empty, both source
        // count sentinels remain -1, Wilson B remains NaN, and the nested
        // reflection record uses its own source-defined NaN defaults.
        // Complexity review: initialization is constant work with no heap
        // allocation for empty strings/vectors and nested scalar defaults.
        Self {
            method: String::new(),
            number_of_crystals: -1,
            unique_reflections: -1,
            reflections: BioReflectionsInfo::default(),
            b_wilson: f64::NAN,
            shells: Vec::new(),
            diffraction_ids: Vec::new(),
        }
    }
}

/// Source-ordered diffraction metadata for an experiment.
#[derive(Debug, Clone, PartialEq)]
pub struct BioDiffractionInfo {
    pub id: String,
    pub temperature: f64,
    pub source: String,
    pub source_type: String,
    pub synchrotron: String,
    pub beamline: String,
    pub wavelengths: String,
    pub scattering_type: String,
    /// The source's single-byte monochromatic/Laue code; NUL means unset.
    pub mono_or_laue: u8,
    pub monochromator: String,
    pub collection_date: String,
    pub optics: String,
    pub detector: String,
    pub detector_make: String,
}

impl Default for BioDiffractionInfo {
    fn default() -> Self {
        // Gemmi❗✔️: struct DiffractionInfo {
        // Gemmi❗✔️:   std::string id;                // _diffrn.id
        // Gemmi❗✔️:   double temperature = NAN;      // _diffrn.ambient_temp
        // Gemmi❗✔️:   std::string source;            // _diffrn_source.source
        // Gemmi❗✔️:   std::string source_type;       // _diffrn_source.type
        // Gemmi❗✔️:   std::string synchrotron;       // _diffrn_source.pdbx_synchrotron_site
        // Gemmi❗✔️:   std::string beamline;          // _diffrn_source.pdbx_synchrotron_beamline
        // Gemmi❗✔️:   std::string wavelengths;       // _diffrn_source.pdbx_wavelength
        // Gemmi❗✔️:   std::string scattering_type;   // _diffrn_radiation.pdbx_scattering_type
        // Gemmi❗✔️:   char mono_or_laue = '\0'; // _diffrn_radiation.pdbx_monochromatic_or_laue_m_l
        // Gemmi❗✔️:   std::string monochromator;     // _diffrn_radiation.monochromator
        // Gemmi❗✔️:   std::string collection_date;   // _diffrn_detector.pdbx_collection_date
        // Gemmi❗✔️:   std::string optics;            // _diffrn_detector.details
        // Gemmi❗✔️:   std::string detector;          // _diffrn_detector.detector
        // Gemmi❗✔️:   std::string detector_make;     // _diffrn_detector.type
        // Gemmi❗✔️: };
        // Behavior review: default-constructed source strings are empty,
        // temperature is NaN, and the single-byte code starts at NUL. Rust
        // stores the C++ char as u8 so arbitrary source bytes remain representable.
        // Complexity review: this fixed-size value initializes in O(1), with
        // empty String values and no heap allocation, matching the source defaults.
        Self {
            id: String::new(),
            temperature: f64::NAN,
            source: String::new(),
            source_type: String::new(),
            synchrotron: String::new(),
            beamline: String::new(),
            wavelengths: String::new(),
            scattering_type: String::new(),
            mono_or_laue: 0,
            monochromator: String::new(),
            collection_date: String::new(),
            optics: String::new(),
            detector: String::new(),
            detector_make: String::new(),
        }
    }
}

/// Experimental crystal metadata, distinct from lattice/coordinate crystal state.
#[derive(Debug, Clone, PartialEq)]
pub struct BioExperimentalCrystalInfo {
    pub id: String,
    pub description: String,
    pub ph: f64,
    pub ph_range: String,
    pub diffractions: Vec<BioDiffractionInfo>,
}

impl Default for BioExperimentalCrystalInfo {
    fn default() -> Self {
        // Gemmi❗✔️: struct CrystalInfo {
        // Gemmi❗✔️:   std::string id;                 // _exptl_crystal.id
        // Gemmi❗✔️:   std::string description;        // _exptl_crystal.description
        // Gemmi❗✔️:   double ph = NAN;                // _exptl_crystal_grow.pH
        // Gemmi❗✔️:   std::string ph_range;           // _exptl_crystal_grow.pdbx_pH_range
        // Gemmi❗✔️:   std::vector<DiffractionInfo> diffractions;
        // Gemmi❗✔️: };
        // Behavior review: strings and the ordered diffraction collection
        // default empty; pH uses the source NaN sentinel. The explicit Rust
        // name distinguishes experimental metadata from lattice crystal state.
        // Complexity review: fixed scalar/default construction is O(1), with
        // no heap allocation for empty String and Vec values.
        Self {
            id: String::new(),
            description: String::new(),
            ph: f64::NAN,
            ph_range: String::new(),
            diffractions: Vec::new(),
        }
    }
}

/// One TLS residue-range selection from refinement metadata.
///
/// The chain uses COSMolKit's approved four-byte ASCII source-chain boundary.
/// Sequence numbers and insertion codes retain the Gemmi `SeqId` value; an
/// absent sequence number is represented by `i32::MIN`, and no insertion code
/// by `None` (the source blank character).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BioTlsSelection {
    pub chain: crate::PdbChainId,
    pub res_begin: crate::PdbSeqId,
    pub res_end: crate::PdbSeqId,
    pub details: String,
}

impl Default for BioTlsSelection {
    fn default() -> Self {
        // Gemmi❗✔️: struct Selection {
        // Gemmi❗✔️:   std::string chain;
        // Gemmi❗✔️:   SeqId res_begin;
        // Gemmi❗✔️:   SeqId res_end;
        // Gemmi❗✔️:   std::string details;  // _pdbx_refine_tls_group.selection_details
        // Gemmi❗✔️: };
        // Gemmi❗✔️: struct SeqId {
        // Gemmi❗✔️:   using OptionalNum = OptionalInt<INT_MIN>;
        // Gemmi❗✔️:   OptionalNum num;   // sequence number
        // Gemmi❗✔️:   char icode = ' ';  // insertion code
        // Gemmi❗✔️:   SeqId() = default;
        // Gemmi❗✔️: };
        // Gemmi❗✔️: template<int N> struct OptionalInt {
        // Gemmi❗✔️:   enum { None=N };
        // Gemmi❗✔️:   int value = None;
        // Gemmi❗✔️:   OptionalInt() = default;
        // Gemmi❗✔️: };
        // Behavior review: empty strings map to empty chain/details; the
        // bounded chain type preserves the approved ASCII-width contract.
        // `INT_MIN` is the source absent-number sentinel and `None` maps the
        // source blank insertion character. Complexity review: only fixed
        // scalar values and empty strings are initialized, in O(1) time with
        // no heap allocation for default strings.
        Self {
            chain: crate::PdbChainId::from_ascii(b"")
                .expect("empty chain id is within the approved representation"),
            res_begin: crate::PdbSeqId::new(i32::MIN, None),
            res_end: crate::PdbSeqId::new(i32::MIN, None),
            details: String::new(),
        }
    }
}

/// Tensor, origin, and ordered selections for one TLS group.
#[derive(Debug, Clone, PartialEq)]
pub struct BioTlsGroup {
    /// Source's small numeric ID; `-1` means no numeric ID was assigned.
    pub num_id: i16,
    pub id: String,
    pub selections: Vec<BioTlsSelection>,
    pub origin: [f64; 3],
    /// Symmetric components in source order: u11, u22, u33, u12, u13, u23.
    pub t: [f64; 6],
    /// Symmetric components in source order: u11, u22, u33, u12, u13, u23.
    pub l: [f64; 6],
    /// Full row-major 3x3 S matrix.
    pub s: [[f64; 3]; 3],
}

impl Default for BioTlsGroup {
    fn default() -> Self {
        // Gemmi❗✔️: struct TlsGroup {
        // Gemmi❗✔️:   struct Selection {
        // Gemmi❗✔️:     std::string chain;
        // Gemmi❗✔️:     SeqId res_begin;
        // Gemmi❗✔️:     SeqId res_end;
        // Gemmi❗✔️:     std::string details;  // _pdbx_refine_tls_group.selection_details
        // Gemmi❗✔️:   };
        // Gemmi❗✔️:   short num_id = -1;      // id stored as number (optimization)
        // Gemmi❗✔️:   std::string id;         // _pdbx_refine_tls.id
        // Gemmi❗✔️:   std::vector<Selection> selections;
        // Gemmi❗✔️:   Position origin;        // _pdbx_refine_tls.origin_x/y/z
        // Gemmi❗✔️:   SMat33<double> T = {NAN, NAN, NAN, NAN, NAN, NAN};  // _pdbx_refine_tls.T[][]
        // Gemmi❗✔️:   SMat33<double> L = {NAN, NAN, NAN, NAN, NAN, NAN};  // _pdbx_refine_tls.L[][]
        // Gemmi❗✔️:   Mat33 S = Mat33{NAN};   // _pdbx_refine_tls.S[][]
        // Gemmi❗✔️: };
        // Gemmi❗✔️: template<typename T> struct SMat33 {
        // Gemmi❗✔️:   T u11, u22, u33, u12, u13, u23;
        // Gemmi❗✔️: struct Position : Vec3 {
        // Gemmi❗✔️:   using Vec3::Vec3;
        // Gemmi❗✔️:   Position() = default;
        // Gemmi❗✔️: template <typename Real>
        // Gemmi❗✔️: struct Vec3_ {
        // Gemmi❗✔️:   Real x, y, z;
        // Gemmi❗✔️:   Vec3_() : x(0), y(0), z(0) {}
        // Gemmi❗✔️:   double a[3][3] = { {1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.} };
        // Gemmi❗✔️:   Mat33() = default;
        // Gemmi❗✔️:   explicit Mat33(double d) : a{{d, d, d}, {d, d, d}, {d, d, d}} {}
        // Behavior review: map the source short to i16; preserve source
        // missing/default IDs, ordered empty selections, zero Cartesian origin,
        // six NaN symmetric components in declared aggregate order, and all
        // nine NaN S entries. Complexity review: fixed arrays initialize in
        // O(1); the empty String and Vec allocate no backing storage.
        Self {
            num_id: -1,
            id: String::new(),
            selections: Vec::new(),
            origin: [0.0; 3],
            t: [f64::NAN; 6],
            l: [f64::NAN; 6],
            s: [[f64::NAN; 3]; 3],
        }
    }
}

/// Overall refinement metadata, including the source-defined per-bin and
/// restraint records.
#[derive(Debug, Clone, PartialEq)]
pub struct BioRefinementInfo {
    /// The inherited overall refinement statistics.
    pub basic: BioBasicRefinementInfo,
    pub id: String,
    pub cross_validation_method: String,
    pub rfree_selection_method: String,
    pub bin_count: i32,
    pub bins: Vec<BioBasicRefinementInfo>,
    pub mean_b: f64,
    /// Symmetric components in Gemmi's source order: u11, u22, u33, u12, u13, u23.
    pub aniso_b: [f64; 6],
    pub luzzati_error: f64,
    pub dpi_blow_r: f64,
    pub dpi_blow_rfree: f64,
    pub dpi_cruickshank_r: f64,
    pub dpi_cruickshank_rfree: f64,
    pub restr_stats: Vec<BioRefinementRestraint>,
    pub tls_groups: Vec<BioTlsGroup>,
    pub remarks: String,
}

impl Default for BioRefinementInfo {
    fn default() -> Self {
        // Gemmi✔️✔️: struct RefinementInfo : BasicRefinementInfo {
        // Gemmi✔️✔️:   struct Restr {
        // Gemmi✔️✔️:     std::string name;
        // Gemmi✔️✔️:     int count = -1;
        // Gemmi✔️✔️:     double weight = NAN;
        // Gemmi✔️✔️:     std::string function;
        // Gemmi✔️✔️:     double dev_ideal = NAN;
        // Gemmi✔️✔️:
        // Gemmi✔️✔️:     Restr() = default;
        // Gemmi✔️✔️:     explicit Restr(const std::string& name_) : name(name_) {}
        // Gemmi✔️✔️:   };
        // Gemmi✔️✔️:   std::string id;
        // Gemmi✔️✔️:   std::string cross_validation_method; // _refine.pdbx_ls_cross_valid_method
        // Gemmi✔️✔️:   std::string rfree_selection_method;  // _refine.pdbx_R_Free_selection_details
        // Gemmi✔️✔️:   int bin_count = -1;        // _refine_ls_shell.pdbx_total_number_of_bins_used
        // Gemmi✔️✔️:   std::vector<BasicRefinementInfo> bins;
        // Gemmi✔️✔️:   double mean_b = NAN;                // _refine.B_iso_mean
        // Gemmi✔️✔️:   SMat33<double> aniso_b{NAN, NAN, NAN, NAN, NAN, NAN};  // _refine.aniso_B[][]
        // Gemmi✔️✔️:   double luzzati_error = NAN; // _refine_analyze.Luzzati_coordinate_error_obs
        // Gemmi✔️✔️:   double dpi_blow_r = NAN;            // _refine.pdbx_overall_SU_R_Blow_DPI
        // Gemmi✔️✔️:   double dpi_blow_rfree = NAN;        // _refine.pdbx_overall_SU_R_free_Blow_DPI
        // Gemmi✔️✔️:   double dpi_cruickshank_r = NAN;     // _refine.overall_SU_R_Cruickshank_DPI
        // Gemmi✔️✔️:   double dpi_cruickshank_rfree = NAN; // _refine.pdbx_overall_SU_R_free_Cruickshank_DPI
        // Gemmi✔️✔️:   std::vector<Restr> restr_stats;     // _refine_ls_restr
        // Gemmi✔️✔️:   std::vector<TlsGroup> tls_groups;   // _pdbx_refine_tls
        // Gemmi✔️✔️:   std::string remarks;
        // Gemmi✔️✔️: };
        // Behavior review: the Rust `basic` field represents the source base
        // value without duplicating its fifteen fields; bins, restraint rows,
        // and TLS groups retain their declared element types and vector order.
        // The six anisotropic-B values keep SMat33's declared component order.
        // Empty text and vectors, -1 bin count, and NaN statistics match the
        // source member initializers. Complexity review: construction delegates
        // the fixed-size base defaults and initializes fixed-size scalars/array
        // in O(1); empty String and Vec values allocate no backing storage.
        Self {
            basic: BioBasicRefinementInfo::default(),
            id: String::new(),
            cross_validation_method: String::new(),
            rfree_selection_method: String::new(),
            bin_count: -1,
            bins: Vec::new(),
            mean_b: f64::NAN,
            aniso_b: [f64::NAN; 6],
            luzzati_error: f64::NAN,
            dpi_blow_r: f64::NAN,
            dpi_blow_rfree: f64::NAN,
            dpi_cruickshank_r: f64::NAN,
            dpi_cruickshank_rfree: f64::NAN,
            restr_stats: Vec::new(),
            tls_groups: Vec::new(),
            remarks: String::new(),
        }
    }
}

/// Source-ordered metadata attached to a biological structure.
#[derive(Debug, Clone, PartialEq)]
pub struct BioMetadata {
    pub authors: Vec<String>,
    pub experiments: Vec<BioExperimentInfo>,
    pub crystals: Vec<BioExperimentalCrystalInfo>,
    pub refinement: Vec<BioRefinementInfo>,
    pub software: Vec<BioSoftwareItem>,
    pub solved_by: String,
    pub starting_model: String,
    pub remark_300_detail: String,
}

impl Default for BioMetadata {
    fn default() -> Self {
        // Gemmi✔️✔️: struct Metadata {
        // Gemmi✔️✔️:   std::vector<std::string> authors;  // _audit_author.name
        // Gemmi✔️✔️:   std::vector<ExperimentInfo> experiments;
        // Gemmi✔️✔️:   std::vector<CrystalInfo> crystals;
        // Gemmi✔️✔️:   std::vector<RefinementInfo> refinement;
        // Gemmi✔️✔️:   std::vector<SoftwareItem> software;
        // Gemmi✔️✔️:   std::string solved_by;       // _refine.pdbx_method_to_determine_struct
        // Gemmi✔️✔️:   std::string starting_model;  // _refine.pdbx_starting_model
        // Gemmi✔️✔️:   std::string remark_300_detail; // _struct_biol.details
        // Behavior review: all eight source members are retained using their
        // existing canonical BIO value types; each ordered collection and
        // free-text field receives its source default of empty. Complexity
        // review: eight default-empty String/Vec values require fixed O(1)
        // setup and allocate no element/string backing storage.
        // The pinned struct continues with presence and TLS-selection helpers;
        // those writer/heuristic behaviors are outside this passive-value packet.
        Self {
            authors: Vec::new(),
            experiments: Vec::new(),
            crystals: Vec::new(),
            refinement: Vec::new(),
            software: Vec::new(),
            solved_by: String::new(),
            starting_model: String::new(),
            remark_300_detail: String::new(),
        }
    }
}

impl Default for BioSoftwareItem {
    fn default() -> Self {
        // Gemmi❗✔️: struct SoftwareItem {
        // Gemmi❗✔️:   enum Classification {
        // Gemmi❗✔️:     DataCollection, DataExtraction, DataProcessing, DataReduction,
        // Gemmi❗✔️:     DataScaling, ModelBuilding, Phasing, Refinement, Unspecified
        // Gemmi❗✔️:   };
        // Gemmi❗✔️:   std::string name;
        // Gemmi❗✔️:   std::string version;
        // Gemmi❗✔️:   std::string date;
        // Gemmi❗✔️:   std::string description;
        // Gemmi❗✔️:   std::string contact_author;
        // Gemmi❗✔️:   std::string contact_author_email;
        // Gemmi❗✔️:   Classification classification = Unspecified;
        // Gemmi❗✔️: };
        // Behavior review: the six strings default to empty, and the source's
        // classification member initializer maps to `Unspecified`; the enum
        // declaration above retains the pinned source order and integer values.
        // Complexity review: constructing the fixed set of empty String values
        // is constant work with no heap allocation, matching default-constructed
        // std::string members and the enum initializer.
        Self {
            name: String::new(),
            version: String::new(),
            date: String::new(),
            description: String::new(),
            contact_author: String::new(),
            contact_author_email: String::new(),
            classification: BioSoftwareClassification::Unspecified,
        }
    }
}
