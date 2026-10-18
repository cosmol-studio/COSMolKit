use crate::bio_scalar_metadata::{
    BioBasicRefinementInfo, BioDiffractionInfo, BioRefinementRestraint, BioReflectionsInfo,
};
use crate::bio_source_values::PdbSeqId;
use crate::chemical_values::{EnumValue, enum_value};
use cosmolkit as ck;
use pyo3::prelude::*;

// Typed transport only: every field belongs to the exhaustive Rust literal.
// Nested values are copied into/out of detached value objects, never rebuilt
// from guessed fields. The source Default owns every missing-value policy.
macro_rules! metadata_record {
    ($name:ident, scalars {$($field:ident: $ty:ty),+ $(,)?},
        values {$($value:ident: $value_ty:ident),* $(,)?},
        lists {$($list:ident: $list_ty:ident),* $(,)?}) => {
        #[pyclass(module="fixed_projection._native", frozen, from_py_object)]
        #[derive(Clone)]
        pub struct $name { pub(crate) inner: ck::$name }
        #[pymethods]
        impl $name {
            #[new]
            #[pyo3(signature=(*, $($field=None,)+ $($value=None,)* $($list=None),*))]
            fn new($($field: Option<$ty>,)+ $($value: Option<$value_ty>,)* $($list: Option<Vec<$list_ty>>),*) -> Self {
                let defaults = ck::$name::default();
                Self { inner: ck::$name {
                    $($field: $field.unwrap_or(defaults.$field),)+
                    $($value: $value.map(|value| value.inner).unwrap_or(defaults.$value),)*
                    $($list: $list.map(|values| values.into_iter().map(|value| value.inner).collect()).unwrap_or(defaults.$list)),*
                } }
            }
            $(#[getter]
            fn $field(&self) -> &$ty { &self.inner.$field })+
            $(#[getter]
            fn $value(&self) -> $value_ty { $value_ty { inner: self.inner.$value.clone() } })*
            $(#[getter]
            fn $list(&self) -> Vec<$list_ty> {
                self.inner.$list.iter().cloned().map(|inner| $list_ty { inner }).collect()
            })*
        }
    };
}

enum_value!(BioSoftwareClassification, { DataCollection, DataExtraction, DataProcessing,
    DataReduction, DataScaling, ModelBuilding, Phasing, Refinement, Unspecified });
#[derive(FromPyObject, IntoPyObject)]
pub enum BioSoftwareClassificationInput {
    Enum(BioSoftwareClassification),
    Name(String),
}
impl BioSoftwareClassificationInput {
    fn into_rust(self) -> PyResult<ck::BioSoftwareClassification> {
        match self {
            Self::Enum(value) => Ok(value.into()),
            Self::Name(name) => BioSoftwareClassification::from_name(&name).map(Into::into),
        }
    }
}
#[pyclass(module = "fixed_projection._native", frozen, from_py_object)]
#[derive(Clone)]
pub struct BioSoftwareItem {
    pub(crate) inner: ck::BioSoftwareItem,
}
#[pymethods]
impl BioSoftwareItem {
    #[new]
    #[pyo3(signature=(*, name=ck::BioSoftwareItem::default().name, version=ck::BioSoftwareItem::default().version, date=ck::BioSoftwareItem::default().date, description=ck::BioSoftwareItem::default().description, contact_author=ck::BioSoftwareItem::default().contact_author, contact_author_email=ck::BioSoftwareItem::default().contact_author_email, classification=None))]
    fn new(
        name: String,
        version: String,
        date: String,
        description: String,
        contact_author: String,
        contact_author_email: String,
        classification: Option<BioSoftwareClassificationInput>,
    ) -> PyResult<Self> {
        let classification = classification
            .map(BioSoftwareClassificationInput::into_rust)
            .transpose()?
            .unwrap_or_else(|| ck::BioSoftwareItem::default().classification);
        Ok(Self {
            inner: ck::BioSoftwareItem {
                name,
                version,
                date,
                description,
                contact_author,
                contact_author_email,
                classification,
            },
        })
    }
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }
    #[getter]
    fn version(&self) -> &str {
        &self.inner.version
    }
    #[getter]
    fn date(&self) -> &str {
        &self.inner.date
    }
    #[getter]
    fn description(&self) -> &str {
        &self.inner.description
    }
    #[getter]
    fn contact_author(&self) -> &str {
        &self.inner.contact_author
    }
    #[getter]
    fn contact_author_email(&self) -> &str {
        &self.inner.contact_author_email
    }
    #[getter]
    fn classification(&self) -> BioSoftwareClassification {
        self.inner.classification.into()
    }
}

#[pyclass(module = "fixed_projection._native", frozen, from_py_object)]
#[derive(Clone)]
pub struct BioTlsSelection {
    pub(crate) inner: ck::BioTlsSelection,
}
#[pymethods]
impl BioTlsSelection {
    #[new]
    #[pyo3(signature=(*, chain=None, res_begin=None, res_end=None, details=ck::BioTlsSelection::default().details))]
    fn new(
        chain: Option<&str>,
        res_begin: Option<PdbSeqId>,
        res_end: Option<PdbSeqId>,
        details: String,
    ) -> PyResult<Self> {
        let defaults = ck::BioTlsSelection::default();
        let chain = chain
            .map(|value| {
                ck::PdbChainId::from_ascii(value.as_bytes()).ok_or_else(|| {
                    pyo3::exceptions::PyValueError::new_err("chain is not a valid PdbChainId")
                })
            })
            .transpose()?
            .unwrap_or(defaults.chain);
        Ok(Self {
            inner: ck::BioTlsSelection {
                chain,
                res_begin: res_begin
                    .map(|value| value.inner)
                    .unwrap_or(defaults.res_begin),
                res_end: res_end.map(|value| value.inner).unwrap_or(defaults.res_end),
                details,
            },
        })
    }
    #[getter]
    fn chain(&self) -> &str {
        self.inner.chain.as_str()
    }
    #[getter]
    fn res_begin(&self) -> PdbSeqId {
        PdbSeqId {
            inner: self.inner.res_begin,
        }
    }
    #[getter]
    fn res_end(&self) -> PdbSeqId {
        PdbSeqId {
            inner: self.inner.res_end,
        }
    }
    #[getter]
    fn details(&self) -> &str {
        &self.inner.details
    }
}

metadata_record!(BioExperimentInfo,
    scalars { method: String, number_of_crystals: i32, unique_reflections: i32,
        b_wilson: f64, diffraction_ids: Vec<String> },
    values { reflections: BioReflectionsInfo }, lists { shells: BioReflectionsInfo });
metadata_record!(
    BioExperimentalCrystalInfo,
    scalars {
        id: String,
        description: String,
        ph: f64,
        ph_range: String
    },
    values {},
    lists {
        diffractions: BioDiffractionInfo
    }
);
metadata_record!(
    BioTlsGroup,
    scalars {
        num_id: i16,
        id: String,
        origin: [f64; 3],
        t: [f64; 6],
        l: [f64; 6],
        s: [[f64; 3]; 3]
    },
    values {},
    lists {
        selections: BioTlsSelection
    }
);
metadata_record!(
    BioRefinementInfo,
    scalars {
        id: String,
        cross_validation_method: String,
        rfree_selection_method: String,
        bin_count: i32,
        mean_b: f64,
        aniso_b: [f64; 6],
        luzzati_error: f64,
        dpi_blow_r: f64,
        dpi_blow_rfree: f64,
        dpi_cruickshank_r: f64,
        dpi_cruickshank_rfree: f64,
        remarks: String
    },
    values {
        basic: BioBasicRefinementInfo
    },
    lists {
        bins: BioBasicRefinementInfo,
        restr_stats: BioRefinementRestraint,
        tls_groups: BioTlsGroup
    }
);
metadata_record!(BioMetadata,
    scalars { authors: Vec<String>, solved_by: String, starting_model: String, remark_300_detail: String },
    values {}, lists { experiments: BioExperimentInfo, crystals: BioExperimentalCrystalInfo,
        refinement: BioRefinementInfo, software: BioSoftwareItem });
