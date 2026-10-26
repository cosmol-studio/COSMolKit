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
    ($(#[$class_doc:meta])* $name:ident, scalars {$($(#[$field_doc:meta])* $field:ident: $ty:ty),+ $(,)?},
        values {$($(#[$value_doc:meta])* $value:ident: $value_ty:ident),* $(,)?},
        lists {$($(#[$list_doc:meta])* $list:ident: $list_ty:ident),* $(,)?}) => {
        crate::configuration::type_projection! {
            #[pyclass(module="cosmolkit", frozen, from_py_object)]
            #[derive(Clone)]
            $(#[$class_doc])*
            pub struct $name {
                pub(crate) inner: ck::$name
            }
        }
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
            $(#[$field_doc])*
            fn $field(&self) -> &$ty { &self.inner.$field })+
            $(#[getter]
            $(#[$value_doc])*
            fn $value(&self) -> $value_ty { $value_ty { inner: self.inner.$value.clone() } })*
            $(#[getter]
            $(#[$list_doc])*
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
crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Software provenance reported during structure determination or refinement.
    pub struct BioSoftwareItem {
        pub(crate) inner: ck::BioSoftwareItem,
    }
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
    /// Name of the software reported in the source.
    fn name(&self) -> &str {
        &self.inner.name
    }
    #[getter]
    /// Reported software version.
    fn version(&self) -> &str {
        &self.inner.version
    }
    #[getter]
    /// Reported software date.
    fn date(&self) -> &str {
        &self.inner.date
    }
    #[getter]
    /// Source description of the software.
    fn description(&self) -> &str {
        &self.inner.description
    }
    #[getter]
    /// Contact author recorded for the software.
    fn contact_author(&self) -> &str {
        &self.inner.contact_author
    }
    #[getter]
    /// Contact author's email address.
    fn contact_author_email(&self) -> &str {
        &self.inner.contact_author_email
    }
    #[getter]
    /// Role of the software in structure determination.
    fn classification(&self) -> BioSoftwareClassification {
        self.inner.classification.into()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Chain and residue range to which a TLS refinement group applies.
    pub struct BioTlsSelection {
        pub(crate) inner: ck::BioTlsSelection,
    }
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
    /// Chain identifier selected by this TLS range.
    fn chain(&self) -> &str {
        self.inner.chain.as_str()
    }
    #[getter]
    /// First residue sequence identifier of the TLS selection.
    fn res_begin(&self) -> PdbSeqId {
        PdbSeqId {
            inner: self.inner.res_begin,
        }
    }
    #[getter]
    /// Last residue sequence identifier of the TLS selection.
    fn res_end(&self) -> PdbSeqId {
        PdbSeqId {
            inner: self.inner.res_end,
        }
    }
    #[getter]
    /// Additional source description of this TLS selection.
    fn details(&self) -> &str {
        &self.inner.details
    }
}

metadata_record!(
/// Experimental method and diffraction/reflection statistics reported for a structure.
BioExperimentInfo,
scalars {
    /// Experimental method recorded by the source.
    method: String,
    /// Reported number of crystals used.
    number_of_crystals: i32,
    /// Reported number of unique reflections.
    unique_reflections: i32,
    /// Wilson B factor reported for the experiment.
    b_wilson: f64,
    /// Identifiers linking this experiment to diffraction records.
    diffraction_ids: Vec<String>
},
values {
    /// Overall reflection statistics.
    reflections: BioReflectionsInfo
}, lists {
    /// Reflection statistics by resolution shell.
    shells: BioReflectionsInfo
});
metadata_record!(
    /// Crystal identifier, preparation description, pH and diffraction records.
    BioExperimentalCrystalInfo,
    scalars {
        /// Crystal identifier in the source.
        id: String,
        /// Source description of crystal preparation or properties.
        description: String,
        /// Reported crystal preparation pH.
        ph: f64,
        /// Source text describing the pH range.
        ph_range: String
    },
    values {},
    lists {
        /// Diffraction measurements associated with this crystal.
        diffractions: BioDiffractionInfo
    }
);
metadata_record!(
    /// Translation/libration/screw refinement group, its origin and atom selections.
    BioTlsGroup,
    scalars {
        /// Numeric source identifier for this TLS group.
        num_id: i16,
        /// Textual source identifier for this TLS group.
        id: String,
        /// Cartesian origin of the TLS model.
        origin: [f64; 3],
        /// Six independent components of the symmetric translation tensor.
        t: [f64; 6],
        /// Six independent components of the symmetric libration tensor.
        l: [f64; 6],
        /// Three-by-three screw tensor.
        s: [[f64; 3]; 3]
    },
    values {},
    lists {
        /// Chain/residue selections belonging to this TLS group.
        selections: BioTlsSelection
    }
);
metadata_record!(
    /// Refinement methods, quality statistics and retained restraint/TLS records.
    BioRefinementInfo,
    scalars {
        /// Refinement identifier from the source.
        id: String,
        /// Reported cross-validation method.
        cross_validation_method: String,
        /// Method used to select the free-R reflection set.
        rfree_selection_method: String,
        /// Reported number of resolution bins.
        bin_count: i32,
        /// Mean atomic B factor.
        mean_b: f64,
        /// Six components of the reported anisotropic B tensor.
        aniso_b: [f64; 6],
        /// Coordinate error estimated by the Luzzati method.
        luzzati_error: f64,
        /// Blow diffraction precision index calculated from R.
        dpi_blow_r: f64,
        /// Blow diffraction precision index calculated from free R.
        dpi_blow_rfree: f64,
        /// Cruickshank diffraction precision index calculated from R.
        dpi_cruickshank_r: f64,
        /// Cruickshank diffraction precision index calculated from free R.
        dpi_cruickshank_rfree: f64,
        /// Additional source refinement remarks.
        remarks: String
    },
    values {
        /// Overall refinement statistics.
        basic: BioBasicRefinementInfo
    },
    lists {
        /// Refinement statistics for each resolution bin.
        bins: BioBasicRefinementInfo,
        /// Reported restraint statistics.
        restr_stats: BioRefinementRestraint,
        /// Translation/libration/screw refinement groups.
        tls_groups: BioTlsGroup
    }
);
metadata_record!(
/// Detached experimental, authorship and refinement metadata for a biological structure.
BioMetadata,
scalars {
    /// Source author names in retained order.
    authors: Vec<String>,
    /// Reported structure solution method.
    solved_by: String,
    /// Description of the starting structural model.
    starting_model: String,
    /// Retained PDB REMARK 300 details.
    remark_300_detail: String
},
values {}, lists {
    /// Experimental records.
    experiments: BioExperimentInfo,
    /// Crystal preparation and diffraction records.
    crystals: BioExperimentalCrystalInfo,
    /// Refinement records.
    refinement: BioRefinementInfo,
    /// Software records reported by the source.
    software: BioSoftwareItem
});
