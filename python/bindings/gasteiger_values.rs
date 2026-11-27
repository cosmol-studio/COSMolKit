//! Checked Python projections for Gasteiger options and detached atom rows.

use cosmolkit as ck;
use pyo3::{exceptions::PyValueError, prelude::*};

#[pyclass(module = "cosmolkit", from_py_object)]
#[derive(Clone)]
/// Checked settings for the Gasteiger partial-charge iteration.
pub struct GasteigerChargeParams {
    pub(crate) inner: ck::GasteigerChargeParams,
}

// Checked native accessors use the existing explicit configuration adapter hook.
// The native declaration supplies the complete field list and storage types.
macro_rules! checked_fields {
    ($(#[$count_doc:meta])* $count:ident: $count_ty:ty,
     $(#[$strict_doc:meta])* $strict:ident: $strict_ty:ty) => {
        paste::paste! {
            #[pymethods]
            impl GasteigerChargeParams {
                #[new]
                #[pyo3(signature = ($count=i64::from(ck::GasteigerChargeParams::default().$count()), $strict=ck::GasteigerChargeParams::default().$strict()))]
                pub(crate) fn new($count: i64, $strict: $strict_ty) -> PyResult<Self> {
                    Ok(Self { inner: ck::GasteigerChargeParams::new($count, $strict)
                        .map_err(|error| PyValueError::new_err(error.to_string()))? })
                }
                $(#[$count_doc])*
                #[getter]
                fn $count(&self) -> $count_ty { self.inner.$count() }
                #[setter]
                fn [<set_ $count>](&mut self, value: i64) -> PyResult<()> {
                    self.inner.[<set_ $count>](value)
                        .map_err(|error| PyValueError::new_err(error.to_string()))
                }
                $(#[$strict_doc])*
                #[getter]
                fn $strict(&self) -> $strict_ty { self.inner.$strict() }
                #[setter]
                fn [<set_ $strict>](&mut self, value: $strict_ty) {
                    self.inner.[<set_ $strict>](value);
                }
            }
        }
    };
}
ck::project_configuration!(ck::GasteigerChargeParams, using checked_fields);

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Detached Gasteiger charge values for one graph atom.
    pub struct GasteigerCharge {
        pub(crate) inner: ck::GasteigerCharge,
    }
}

#[pymethods]
impl GasteigerCharge {
    #[getter]
    /// Charge on this graph atom, in elementary-charge units.
    fn atom_charge(&self) -> f64 {
        self.inner.atom_charge()
    }

    #[getter]
    /// Sum of charges on this atom's attached hydrogens that are not graph vertices.
    fn attached_hydrogen_charge(&self) -> f64 {
        self.inner.attached_hydrogen_charge()
    }
}
