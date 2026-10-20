use crate::callbacks::{CallbackErrors, typed_callback};
use cosmolkit as ck;
use pyo3::prelude::*;
use std::sync::Arc;

typed_callback!(
    ProgressCallback,
    [],
    pyo3::inspect::PyStaticExpr::Constant {
        value: pyo3::inspect::PyStaticConstant::None
    }
);

#[pyclass(module = "cosmolkit")]
pub struct BatchQueryParams {
    inner: ck::BatchQueryParams,
    progress_callback: Option<ProgressCallback>,
}
macro_rules! batch_query_configuration {
    ($($field:ident: $ty:ty),+ $(,)?) => {
    paste::paste! {
    #[pymethods]
    impl BatchQueryParams {
    fn __traverse__(&self, visit: pyo3::gc::PyVisit<'_>) -> Result<(), pyo3::PyTraverseError> {
        if let Some(callback) = &self.progress_callback {
            callback.traverse(&visit)?;
        }
        Ok(())
    }
    fn __clear__(&mut self) {
        self.progress_callback = None;
    }
    #[new]
    #[pyo3(signature=(*, $($field=ck::BatchQueryParams::default().$field,)+ progress_callback=None))]
    fn new(
        $($field: $ty,)+
        progress_callback: Option<ProgressCallback>,
    ) -> Self {
        Self {
            inner: ck::BatchQueryParams {
                $($field,)+
                progress_callback: None,
            },
            progress_callback,
        }
    }
    $(#[getter] fn $field(&self) -> $ty { self.inner.$field }
    #[setter] fn [<set_ $field>](&mut self, value: $ty) { self.inner.$field = value; })+
    #[getter]
    fn progress_callback(&self, py: Python<'_>) -> Option<ProgressCallback> {
        self.progress_callback
            .as_ref()
            .map(|value| value.clone_ref(py))
    }
    #[setter]
    fn set_progress_callback(&mut self, value: Option<ProgressCallback>) {
        self.progress_callback = value;
    }
    }
    }
    };
}
ck::project_configuration!(ck::BatchQueryParams, using batch_query_configuration);
impl BatchQueryParams {
    pub(crate) fn with_callback<T: Send, E: Send>(
        &self,
        py: Python<'_>,
        run: impl FnOnce(&ck::BatchQueryParams) -> Result<T, E> + Send,
    ) -> PyResult<Result<T, E>> {
        let mut params = self.inner.clone();
        let errors = Arc::new(CallbackErrors::new());
        if let Some(callback) = &self.progress_callback {
            let callback = callback.clone_ref(py);
            let errors = Arc::clone(&errors);
            params.progress_callback = Some(Arc::new(move || {
                errors.notify(|py| {
                    callback.bind(py).call0()?;
                    Ok(())
                });
            }));
        }
        // The native batch owns parallel scheduling and row order. Worker
        // callbacks need Python, so release the caller's interpreter while it runs.
        let result = py.detach(|| run(&params));
        match errors.take() {
            Some(error) => Err(error),
            None => Ok(result),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pyo3::ffi::c_str;
    use pyo3::types::PyModule;

    #[test]
    fn native_parallel_batch_progress_preserves_errors_and_recovers() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let callbacks = PyModule::from_code(
                py,
                c_str!(
                    r#"
seen = []
failed = []
marker = ValueError("original progress failure")
def progress():
    seen.append(True)
def fail():
    failed.append(True)
    raise marker
"#
                ),
                c_str!("batch_progress_probe.py"),
                c_str!("batch_progress_probe"),
            )?;
            let batch = ck::MoleculeBatch::from_smiles_list(&[
                "CCO".to_owned(),
                "C".to_owned(),
                "O".to_owned(),
            ])
            .unwrap();
            let mut params = BatchQueryParams::new(
                Some(2),
                Some(false),
                Some(callbacks.getattr("progress")?.extract()?),
            );
            let run = |native: &ck::BatchQueryParams| {
                batch.to_smiles_with_params(&ck::SmilesWriteParams::default(), native)
            };
            let rows = params.with_callback(py, run)?.unwrap();
            assert_eq!(rows.len(), 3);
            assert_eq!(callbacks.getattr("seen")?.len()?, 3);
            assert_eq!(rows[0].as_ref().unwrap().as_bytes(), b"CCO");
            params.progress_callback = Some(callbacks.getattr("fail")?.extract()?);
            let error = params.with_callback(py, run).unwrap_err();
            assert!(error.value(py).is(&callbacks.getattr("marker")?));
            assert_eq!(callbacks.getattr("failed")?.len()?, 3);
            params.progress_callback = Some(callbacks.getattr("progress")?.extract()?);
            assert_eq!(params.with_callback(py, run)?.unwrap().len(), 3);
            assert_eq!(callbacks.getattr("seen")?.len()?, 6);
            Ok(())
        })
        .unwrap();
    }
}
