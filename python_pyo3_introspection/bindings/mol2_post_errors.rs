//! The concrete MOL2 finalization record and its existing native source graph.
use crate::error_framework::Source;
use crate::hydrogen_errors::HydrogenError;
use crate::io_source_errors::ErrorRef;
use crate::model_edit_errors::MoleculePropertyError;
use crate::model_errors::CoordinateValidationError;
use crate::mol_post_errors::MolPostError;
use crate::retained_error::RetainedError;
use crate::ring_errors::RingFindingError;
use crate::sanitize_errors::SanitizeError;
use crate::stereo_errors::LegacyStereoError;
use crate::stereo_errors::StereoError;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::{exceptions::PyValueError, prelude::*};
type Mol2Source = Source<
    Py<ValenceError>,
    Source<
        Py<StereoError>,
        Source<
            Py<SanitizeError>,
            Source<
                Py<MolPostError>,
                Source<
                    Py<HydrogenError>,
                    Source<
                        Py<RingFindingError>,
                        Source<
                            Py<CoordinateValidationError>,
                            Source<Py<LegacyStereoError>, Py<MoleculePropertyError>>,
                        >,
                    >,
                >,
            >,
        >,
    >,
>;

// The private native constructor is used only at the nine concrete stage types
// in mol2_post.rs. Match the actual Error type, never stage or Display text.
fn source(node: ErrorRef) -> PyErr {
    if let Some(value) = node.error().downcast_ref::<ck::ValenceError>() {
        return ValenceError::from_rust(value.clone());
    }
    if let Some(value) = node.error().downcast_ref::<ck::StereoError>() {
        return StereoError::from_rust(value.clone());
    }
    if let Some(value) = node.error().downcast_ref::<ck::SanitizeError>() {
        return SanitizeError::from_rust(value.clone());
    }
    if let Some(value) = node.error().downcast_ref::<ck::MolPostError>() {
        return MolPostError::from_rust(value.clone());
    }
    if let Some(value) = node.error().downcast_ref::<ck::HydrogenError>() {
        return HydrogenError::from_rust(value.clone());
    }
    if let Some(value) = node.error().downcast_ref::<ck::RingFindingError>() {
        return RingFindingError::from_rust(value.clone());
    }
    if let Some(value) = node.error().downcast_ref::<ck::CoordinateValidationError>() {
        return CoordinateValidationError::from_rust(value.clone());
    }
    if let Some(value) = node.error().downcast_ref::<ck::LegacyStereoError>() {
        return LegacyStereoError::from_rust(value.clone());
    }
    let value = node
        .error()
        .downcast_ref::<ck::MoleculePropertyError>()
        .expect("actual MOL2 StereoMetadata source");
    MoleculePropertyError::from_rust(value.clone())
}

#[doc = r#"{"role":"error_record"}"#]
#[pyclass(module="cosmolkit",extends=PyValueError)]
pub struct Mol2PostError {
    inner: RetainedError<ck::Mol2PostError>,
}
impl Mol2PostError {
    pub(crate) fn from_rust(inner: ck::Mol2PostError) -> PyErr {
        Self::from_retained(RetainedError::new(inner))
    }
    pub(crate) fn from_ref(node: ErrorRef) -> PyErr {
        Self::from_retained(RetainedError::from_ref(node))
    }
    fn from_retained(inner: RetainedError<ck::Mol2PostError>) -> PyErr {
        let cause = source(
            inner
                .node()
                .source()
                .expect("actual Mol2PostError boxed source"),
        );
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => {
                let error = PyErr::from_value(value.into_bound(py).into_any());
                error.set_cause(py, Some(cause));
                error
            }
            Err(error) => error,
        })
    }
}
#[pymethods]
impl Mol2PostError {
    #[getter]
    fn domain(&self) -> &'static str {
        "io"
    }
    #[getter]
    fn stage(&self) -> &'static str {
        self.inner.stage
    }
    #[getter]
    #[doc = r#"{"role":"error_source"}"#]
    fn source(slf: Bound<'_, Self>) -> PyResult<Mol2Source> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::error::Error;
    #[test]
    fn real_mol2_finalization_retains_native_stage_exact_cause_and_owned_descendants() {
        Python::initialize();
        Python::attach(|py|->PyResult<()> {
            let text="@<TRIPOS>MOLECULE\ninvalid_valence\n6 5 0 0 0\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n1 C1 0.0 0.0 0.0 C.3 1 MOL 0.0\n2 C2 1.0 0.0 0.0 C.3 1 MOL 0.0\n3 C3 -1.0 0.0 0.0 C.3 1 MOL 0.0\n4 C4 0.0 1.0 0.0 C.3 1 MOL 0.0\n5 C5 0.0 -1.0 0.0 C.3 1 MOL 0.0\n6 C6 0.0 0.0 1.0 C.3 1 MOL 0.0\n@<TRIPOS>BOND\n1 1 2 1\n2 1 3 1\n3 1 4 1\n4 1 5 1\n5 1 6 1\n";
            let error=ck::Molecule::from_mol2_with_params(text,&ck::Mol2ReadParams{remove_hs:false,..Default::default()}).unwrap_err();
            let ck::MolecularIoError::Mol2Post(native)=error else {panic!("actual MOL2 finalization failure expected, got {error:?}")};
            let display=native.to_string();
            let stage=native.stage;
            let native_source=native.source().expect("real boxed finalization source");
            let source_display=native_source.to_string();
            let expected=if native_source.is::<ck::ValenceError>(){"ValenceError"}
                else if native_source.is::<ck::StereoError>(){"StereoError"}
                else if native_source.is::<ck::SanitizeError>(){"SanitizeError"}
                else if native_source.is::<ck::MolPostError>(){"MolPostError"}
                else if native_source.is::<ck::HydrogenError>(){"HydrogenError"}
                else if native_source.is::<ck::RingFindingError>(){"RingFindingError"}
                else if native_source.is::<ck::CoordinateValidationError>(){"CoordinateValidationError"}
                else if native_source.is::<ck::LegacyStereoError>(){"LegacyStereoError"}
                else if native_source.is::<ck::MoleculePropertyError>(){"MoleculePropertyError"}
                else{panic!("the actual closed native source type must be preserved")};
            let error=Mol2PostError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?,display);
            assert_eq!(error.value(py).getattr("stage")?.extract::<String>()?,stage);
            assert!(error.value(py).setattr("stage","different").is_err());
            assert!(!error.value(py).hasattr("variant")?);
            let cause=error.cause(py).expect("native cause must remain");
            assert_eq!(cause.value(py).get_type().name()?.to_str()?,expected);
            assert_eq!(cause.value(py).str()?.to_str()?,source_display);
            assert!(error.value(py).getattr("source")?.is(cause.value(py)));
            let mut cursor=cause.clone_ref(py);
            let mut observed=Vec::new();
            loop {
                observed.push(cursor.value(py).str()?.to_str()?.to_owned());
                match cursor.cause(py){Some(next)=>cursor=next,None=>break}
            }
            drop(error);
            assert_eq!(cause.value(py).str()?.to_str()?,source_display);
            assert!(!observed.is_empty());
            assert!(!cause.value(py).getattr("domain")?.extract::<String>()?.is_empty());
            assert!(cause.value(py).setattr("domain","changed").is_err());
            Ok(())
        }).unwrap();
    }
}
