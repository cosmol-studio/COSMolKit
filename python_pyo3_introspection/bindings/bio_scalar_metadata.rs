use crate::bio_values::BioTransform;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::collections::BTreeMap;

// Native readonly scalar records. Omission/None selects the exact Rust field
// default, including source-defined NaN and negative sentinel values. Borrowed
// getters avoid cloning Rust strings before Python creates its own value.
macro_rules! scalar_record {
    ($name:ident, {$($field:ident: $ty:ty),+ $(,)?}) => {
        #[pyclass(module="cosmolkit", frozen, from_py_object)]
        #[derive(Clone)]
        pub struct $name { pub(crate) inner: ck::$name }
        #[pymethods]
        impl $name {
            #[new]
            #[pyo3(signature=(*, $($field=None),+))]
            fn new($($field: Option<$ty>),+) -> Self {
                Self { inner: ck::$name { $($field: $field.unwrap_or_else(|| ck::$name::default().$field)),+ } }
            }
            $(#[getter]
            fn $field(&self) -> &$ty { &self.inner.$field })+
        }
    };
}
scalar_record!(BioReflectionsInfo, {
    resolution_high: f64, resolution_low: f64, completeness: f64, redundancy: f64,
    r_merge: f64, r_sym: f64, mean_i_over_sigma: f64,
});
scalar_record!(BioBasicRefinementInfo, {
    resolution_high: f64, resolution_low: f64, completeness: f64,
    reflection_count: i32, work_set_count: i32, rfree_set_count: i32,
    r_all: f64, r_work: f64, r_free: f64, cc_fo_fc_work: f64, cc_fo_fc_free: f64,
    fsc_work: f64, fsc_free: f64, cc_intensity_work: f64, cc_intensity_free: f64,
});
scalar_record!(BioRefinementRestraint, {
    name: String, count: i32, weight: f64, function: String, dev_ideal: f64,
});
scalar_record!(BioDiffractionInfo, {
    id: String, temperature: f64, source: String, source_type: String,
    synchrotron: String, beamline: String, wavelengths: String, scattering_type: String,
    mono_or_laue: u8, monochromator: String, collection_date: String,
    optics: String, detector: String, detector_make: String,
});

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub struct BioStructureSourceState {
    pub(crate) inner: ck::BioStructureSourceState,
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
    fn name(&self) -> &str {
        &self.inner.name
    }
    #[getter]
    fn resolution(&self) -> f64 {
        self.inner.resolution
    }
    #[getter]
    fn conect_map(&self) -> &BTreeMap<i32, Vec<i32>> {
        &self.inner.conect_map
    }
    #[getter]
    fn has_d_fraction(&self) -> bool {
        self.inner.has_d_fraction
    }
    #[getter]
    fn non_ascii_line(&self) -> i32 {
        self.inner.non_ascii_line
    }
    #[getter]
    fn ter_status(&self) -> u8 {
        self.inner.ter_status
    }
    #[getter]
    fn has_origx(&self) -> bool {
        self.inner.has_origx
    }
    #[getter]
    fn origx(&self) -> BioTransform {
        BioTransform {
            inner: self.inner.origx,
        }
    }
    #[getter]
    fn info(&self) -> &BTreeMap<String, String> {
        &self.inner.info
    }
    #[getter]
    fn raw_remarks(&self) -> &[String] {
        &self.inner.raw_remarks
    }
}
