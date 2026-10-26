use crate::bio_values::BioTransform;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::collections::BTreeMap;

// Native readonly scalar records. Omission/None selects the exact Rust field
// default, including source-defined NaN and negative sentinel values. Borrowed
// getters avoid cloning Rust strings before Python creates its own value.
macro_rules! scalar_record {
    ($(#[$doc:meta])* $name:ident, {$($(#[$field_doc:meta])* $field:ident: $ty:ty),+ $(,)?}) => {
        crate::configuration::type_projection! {
            #[pyclass(module="cosmolkit", frozen, from_py_object)]
            #[derive(Clone)]
            $(#[$doc])*
            pub struct $name {
                pub(crate) inner: ck::$name
            }
        }
        #[pymethods]
        impl $name {
            #[new]
            #[pyo3(signature=(*, $($field=None),+))]
            fn new($($field: Option<$ty>),+) -> Self {
                Self { inner: ck::$name { $($field: $field.unwrap_or_else(|| ck::$name::default().$field)),+ } }
            }
            $($(#[$field_doc])* #[getter]
            fn $field(&self) -> &$ty { &self.inner.$field })+
        }
    };
}
scalar_record!(
    /// Readonly diffraction-reflection statistics. Omitted values retain native missing-value sentinels.
    BioReflectionsInfo, {
    /// High-resolution limit in angstroms.
    resolution_high: f64,
    /// Low-resolution limit in angstroms.
    resolution_low: f64,
    /// Reported reflection completeness.
    completeness: f64,
    /// Reported measurement redundancy.
    redundancy: f64,
    /// Merging R factor.
    r_merge: f64,
    /// Symmetry-related R factor.
    r_sym: f64,
    /// Mean intensity divided by its standard uncertainty.
    mean_i_over_sigma: f64,
});
scalar_record!(
    /// Readonly refinement statistics for reflection work and free sets; unavailable values retain native sentinels.
    BioBasicRefinementInfo, {
    /// High-resolution refinement limit in angstroms.
    resolution_high: f64,
    /// Low-resolution refinement limit in angstroms.
    resolution_low: f64,
    /// Reported completeness of the refinement data.
    completeness: f64,
    /// Total reflection count.
    reflection_count: i32,
    /// Number of reflections in the working set.
    work_set_count: i32,
    /// Number of reflections in the R-free set.
    rfree_set_count: i32,
    /// R factor over all reflections.
    r_all: f64,
    /// R factor over the working reflections.
    r_work: f64,
    /// R factor over the free reflections.
    r_free: f64,
    /// Observed/calculated structure-factor correlation for the working set.
    cc_fo_fc_work: f64,
    /// Observed/calculated structure-factor correlation for the free set.
    cc_fo_fc_free: f64,
    /// Fourier-shell correlation reported for the working set.
    fsc_work: f64,
    /// Fourier-shell correlation reported for the free set.
    fsc_free: f64,
    /// Intensity correlation for the working set.
    cc_intensity_work: f64,
    /// Intensity correlation for the free set.
    cc_intensity_free: f64,
});
scalar_record!(
    /// Readonly statistics for one refinement restraint category.
    BioRefinementRestraint, {
    /// Restraint category name.
    name: String,
    /// Number of restraints in this category.
    count: i32,
    /// Reported restraint weight.
    weight: f64,
    /// Restraint function description from the source.
    function: String,
    /// Reported deviation from the ideal restraint target.
    dev_ideal: f64,
});
scalar_record!(
    /// Readonly diffraction-experiment metadata preserved from the structural file.
    BioDiffractionInfo, {
    /// Source diffraction-experiment identifier.
    id: String,
    /// Reported experiment temperature in kelvin.
    temperature: f64,
    /// Radiation source description.
    source: String,
    /// Radiation source type.
    source_type: String,
    /// Synchrotron facility name.
    synchrotron: String,
    /// Beamline identifier.
    beamline: String,
    /// Source wavelength description, preserved as text.
    wavelengths: String,
    /// Scattering experiment type.
    scattering_type: String,
    /// Native single-byte monochromatic/Laue classification marker.
    mono_or_laue: u8,
    /// Monochromator description.
    monochromator: String,
    /// Data collection date as recorded in the source.
    collection_date: String,
    /// Beam optics description.
    optics: String,
    /// Detector technology description.
    detector: String,
    /// Detector manufacturer or model description.
    detector_make: String,
});

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Source-format metadata and parsing state retained with a biological structure.
    pub struct BioStructureSourceState {
        pub(crate) inner: ck::BioStructureSourceState,
    }
}
#[pymethods]
impl BioStructureSourceState {
    #[new]
    #[pyo3(signature=(*, name=ck::BioStructureSourceState::default().name, resolution=ck::BioStructureSourceState::default().resolution, conect_map=None, has_d_fraction=ck::BioStructureSourceState::default().has_d_fraction, non_ascii_line=ck::BioStructureSourceState::default().non_ascii_line, ter_status=ck::BioStructureSourceState::default().ter_status, has_origx=ck::BioStructureSourceState::default().has_origx, origx=None, info=None, raw_remarks=ck::BioStructureSourceState::default().raw_remarks))]
    fn new(
        name: String,
        resolution: f64,
        conect_map: Option<BTreeMap<i32, Vec<i32>>>,
        has_d_fraction: bool,
        non_ascii_line: i32,
        ter_status: u8,
        has_origx: bool,
        origx: Option<BioTransform>,
        info: Option<BTreeMap<String, String>>,
        raw_remarks: Vec<String>,
    ) -> Self {
        Self {
            inner: ck::BioStructureSourceState {
                name,
                resolution,
                conect_map: conect_map
                    .unwrap_or_else(|| ck::BioStructureSourceState::default().conect_map),
                has_d_fraction,
                non_ascii_line,
                ter_status,
                has_origx,
                origx: origx
                    .map(|value| value.inner)
                    .unwrap_or_else(|| ck::BioStructureSourceState::default().origx),
                info: info.unwrap_or_else(|| ck::BioStructureSourceState::default().info),
                raw_remarks,
            },
        }
    }
    #[getter]
    /// Structure name retained in the source metadata.
    fn name(&self) -> &str {
        &self.inner.name
    }
    #[getter]
    /// Resolution reported by the source structure.
    fn resolution(&self) -> f64 {
        self.inner.resolution
    }
    #[getter]
    /// PDB CONECT relationships keyed by atom serial number.
    fn conect_map(&self) -> &BTreeMap<i32, Vec<i32>> {
        &self.inner.conect_map
    }
    #[getter]
    /// Whether the input supplied deuterium-fraction information.
    fn has_d_fraction(&self) -> bool {
        self.inner.has_d_fraction
    }
    #[getter]
    /// Native source line indicator for non-ASCII input.
    fn non_ascii_line(&self) -> i32 {
        self.inner.non_ascii_line
    }
    #[getter]
    /// Stored native status of PDB TER record handling.
    fn ter_status(&self) -> u8 {
        self.inner.ter_status
    }
    #[getter]
    /// Whether the source supplied an ORIGX transformation.
    fn has_origx(&self) -> bool {
        self.inner.has_origx
    }
    #[getter]
    /// Original-coordinate transformation retained from the source.
    fn origx(&self) -> BioTransform {
        BioTransform {
            inner: self.inner.origx,
        }
    }
    #[getter]
    /// Source metadata key/value pairs.
    fn info(&self) -> &BTreeMap<String, String> {
        &self.inner.info
    }
    #[getter]
    /// Retained raw PDB REMARK records in source order.
    fn raw_remarks(&self) -> &[String] {
        &self.inner.raw_remarks
    }
}
