//! Complete detached MMFF property queries and their concrete table failures.
use crate::chemical_values::enum_value;
use crate::error_framework::{Source, error_projection};
use crate::kekulize_errors::{KekulizeError, KekulizeSource};
use crate::leaf_values::{MmffAtomProperties, MmffVariant};
use crate::ring_errors::RingFindingError;
use crate::sanitize_support_errors::AromaticityError;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;

enum_value!(MmffParamTable, { Def, Prop, Pbci, Chg, Bond, Bndk, Stbn, Dfsb, Angle, Oop, Tor, Vdw, HerschbachLaurie, CovRadPauEle });

crate::enums::binding_enum!(MmffParamParseCauseKind, {
    InvalidUnsigned,
    InvalidFloat,
    MissingToken,
    EmptyProcessedLine,
});
#[pyclass(module = "cosmolkit", frozen)]
#[doc = r#"{"role":"error_value"}"#]
pub struct MmffParamParseCause {
    inner: ck::MmffParamParseCause,
}
#[pymethods]
impl MmffParamParseCause {
    #[getter]
    fn variant(&self) -> MmffParamParseCauseKind {
        match &self.inner {
            ck::MmffParamParseCause::InvalidUnsigned { cell: _ } => {
                MmffParamParseCauseKind::InvalidUnsigned
            }
            ck::MmffParamParseCause::InvalidFloat { cell: _ } => {
                MmffParamParseCauseKind::InvalidFloat
            }
            ck::MmffParamParseCause::MissingToken => MmffParamParseCauseKind::MissingToken,
            ck::MmffParamParseCause::EmptyProcessedLine => {
                MmffParamParseCauseKind::EmptyProcessedLine
            }
        }
    }
    #[getter]
    fn cell(&self) -> Option<&str> {
        match &self.inner {
            ck::MmffParamParseCause::InvalidUnsigned { cell }
            | ck::MmffParamParseCause::InvalidFloat { cell } => Some(cell),
            ck::MmffParamParseCause::MissingToken | ck::MmffParamParseCause::EmptyProcessedLine => {
                None
            }
        }
    }
}

#[pyclass(module="cosmolkit",extends=PyValueError)]
#[doc = r#"{"role":"error_record"}"#]
pub struct MmffParamParseError {
    inner: ck::MmffParamParseError,
}
impl MmffParamParseError {
    pub(crate) fn from_rust(inner: ck::MmffParamParseError) -> PyErr {
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => PyErr::from_value(value.into_bound(py).into_any()),
            Err(error) => error,
        })
    }
}
#[pymethods]
impl MmffParamParseError {
    #[getter]
    fn domain(&self) -> &'static str {
        "mmff"
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn table(&self) -> MmffParamTable {
        self.inner.table().into()
    }
    #[getter]
    fn line(&self) -> usize {
        self.inner.line()
    }
    #[getter]
    fn column(&self) -> usize {
        self.inner.column()
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn cause(&self) -> MmffParamParseCause {
        MmffParamParseCause {
            inner: self.inner.cause().clone(),
        }
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}

error_projection!(MmffAngleLookupError, MmffAngleLookupErrorKind, ck::MmffAngleLookupError,"mmff", {
    MissingDefinition { atom_type:u32 },
},cause=|_:&ck::MmffAngleLookupError| None::<PyErr>);
error_projection!(MmffOopLookupError, MmffOopLookupErrorKind, ck::MmffOopLookupError,"mmff", {
    MissingDefinition { atom_type:u32 },
},cause=|_:&ck::MmffOopLookupError| None::<PyErr>);
error_projection!(MmffTorLookupError, MmffTorLookupErrorKind, ck::MmffTorLookupError,"mmff", {
    MissingDefinition { atom_type:u32 }, InvalidEquivalentLevel { level:usize },
},cause=|_:&ck::MmffTorLookupError| None::<PyErr>);
error_projection!(MmffVdwParseError, MmffVdwParseErrorKind, ck::MmffVdwParseError,"mmff", {
    Table(_),MissingHeader,
},cause=|inner:&ck::MmffVdwParseError| match inner {
    ck::MmffVdwParseError::Table(value) => Some(MmffParamParseError::from_rust(value.clone())),
    ck::MmffVdwParseError::MissingHeader => None,
},source=Option<Py<MmffParamParseError>>);

error_projection!(MmffParamError, MmffParamErrorKind, ck::MmffParamError,"mmff", {
    Table(_), Vdw(_), Angle(_), Torsion(_), Oop(_), MissingDefaultData {symbol:&'static str},
},cause=|inner:&ck::MmffParamError| Python::attach(|py| match inner {
    ck::MmffParamError::Table(value) => MmffParamParseError::from_rust(value.clone()).cause(py),
    ck::MmffParamError::Vdw(value) => MmffVdwParseError::from_rust(value.clone()).cause(py),
    ck::MmffParamError::Angle(value) => MmffAngleLookupError::from_rust(*value).cause(py),
    ck::MmffParamError::Torsion(value) => MmffTorLookupError::from_rust(*value).cause(py),
    ck::MmffParamError::Oop(value) => MmffOopLookupError::from_rust(*value).cause(py),
    ck::MmffParamError::MissingDefaultData {..} => None,
}),source=Option<Py<MmffParamParseError>>,methods={
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn table(&self,py:Python<'_>)->PyResult<Option<Py<MmffParamParseError>>> {match &self.inner {
        ck::MmffParamError::Table(value)=>MmffParamParseError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _=>Ok(None),
    }}
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn vdw(&self,py:Python<'_>)->PyResult<Option<Py<MmffVdwParseError>>> {match &self.inner {
        ck::MmffParamError::Vdw(value)=>MmffVdwParseError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _=>Ok(None),
    }}
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn angle(&self,py:Python<'_>)->PyResult<Option<Py<MmffAngleLookupError>>> {match &self.inner {
        ck::MmffParamError::Angle(value)=>MmffAngleLookupError::from_rust(*value).value(py).extract().map(Some).map_err(Into::into), _=>Ok(None),
    }}
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn torsion(&self,py:Python<'_>)->PyResult<Option<Py<MmffTorLookupError>>> {match &self.inner {
        ck::MmffParamError::Torsion(value)=>MmffTorLookupError::from_rust(*value).value(py).extract().map(Some).map_err(Into::into), _=>Ok(None),
    }}
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn oop(&self,py:Python<'_>)->PyResult<Option<Py<MmffOopLookupError>>> {match &self.inner {
        ck::MmffParamError::Oop(value)=>MmffOopLookupError::from_rust(*value).value(py).extract().map(Some).map_err(Into::into), _=>Ok(None),
    }}
});

type MmffSource = Source<KekulizeSource, Py<MmffParamParseError>>;
error_projection!(MmffMolPropertiesError,MmffMolPropertiesErrorKind,ck::MmffMolPropertiesError,"mmff",{
    Params(_),Kekulize(_),Aromaticity(_),RingFinding(_),Valence(_),
    AtomIndexOutOfRange {atom_index:usize,atoms:usize}, AtomTypePropertiesMissing {atom_type},AtomTypePbciMissing {atom_type},
},cause=|inner:&ck::MmffMolPropertiesError| Python::attach(|py| match inner {
    ck::MmffMolPropertiesError::Params(value)=>MmffParamError::from_rust(value.clone()).cause(py),
    ck::MmffMolPropertiesError::Kekulize(value)=>KekulizeError::from_rust(value.clone()).cause(py),
    ck::MmffMolPropertiesError::Aromaticity(value)=>AromaticityError::from_rust(value.clone()).cause(py),
    ck::MmffMolPropertiesError::RingFinding(value)=>RingFindingError::from_rust(value.clone()).cause(py),
    ck::MmffMolPropertiesError::Valence(value)=>ValenceError::from_rust(value.clone()).cause(py),
    ck::MmffMolPropertiesError::AtomIndexOutOfRange {..}|ck::MmffMolPropertiesError::AtomTypePropertiesMissing {..}|ck::MmffMolPropertiesError::AtomTypePbciMissing {..}=>None,
}),source=Option<MmffSource>,methods={
    #[getter]
    fn atom_type(&self)->Option<u8>{match self.inner {ck::MmffMolPropertiesError::AtomTypePropertiesMissing {atom_type}|ck::MmffMolPropertiesError::AtomTypePbciMissing {atom_type}=>Some(atom_type),_=>None}}
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn params(&self,py:Python<'_>)->PyResult<Option<Py<MmffParamError>>> {match &self.inner {
        ck::MmffMolPropertiesError::Params(value)=>MmffParamError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _=>Ok(None),
    }}
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn kekulize(&self,py:Python<'_>)->PyResult<Option<Py<KekulizeError>>> {match &self.inner {
        ck::MmffMolPropertiesError::Kekulize(value)=>KekulizeError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _=>Ok(None),
    }}
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn aromaticity(&self,py:Python<'_>)->PyResult<Option<Py<AromaticityError>>> {match &self.inner {
        ck::MmffMolPropertiesError::Aromaticity(value)=>AromaticityError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _=>Ok(None),
    }}
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn ring_finding(&self,py:Python<'_>)->PyResult<Option<Py<RingFindingError>>> {match &self.inner {
        ck::MmffMolPropertiesError::RingFinding(value)=>RingFindingError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _=>Ok(None),
    }}
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn valence(&self,py:Python<'_>)->PyResult<Option<Py<ValenceError>>> {match &self.inner {
        ck::MmffMolPropertiesError::Valence(value)=>ValenceError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _=>Ok(None),
    }}
});

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen)]
    pub struct MmffProperties {
        pub(crate) inner: ck::MmffProperties,
    }
}
#[pymethods]
impl MmffProperties {
    #[getter]
    fn is_valid(&self) -> bool {
        self.inner.is_valid()
    }
    #[getter]
    fn variant(&self) -> MmffVariant {
        self.inner.variant().into()
    }
    #[getter]
    fn atoms(&self) -> Vec<MmffAtomProperties> {
        self.inner
            .atoms()
            .iter()
            .copied()
            .map(|inner| MmffAtomProperties { inner })
            .collect()
    }
    fn atom_type(&self, atom_index: usize) -> PyResult<u8> {
        self.inner
            .atom_type(atom_index)
            .map_err(MmffMolPropertiesError::from_rust)
    }
    fn formal_charge(&self, atom_index: usize) -> PyResult<f64> {
        self.inner
            .formal_charge(atom_index)
            .map_err(MmffMolPropertiesError::from_rust)
    }
    fn partial_charge(&self, atom_index: usize) -> PyResult<f64> {
        self.inner
            .partial_charge(atom_index)
            .map_err(MmffMolPropertiesError::from_rust)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::native::Molecule;
    use crate::params::MmffPropertiesParams;
    use pyo3::types::PyDict;

    #[test]
    fn complete_mmff_properties_queries_forward_variants_forms_values_and_concrete_failures() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native = ck::Molecule::from_smiles("CCO").unwrap();
            let before = native.to_binary().unwrap();
            let expected_default = native.mmff_has_all_molecule_params().unwrap();
            let molecule = Py::new(py, Molecule { inner: native })?;
            assert_eq!(
                molecule
                    .bind(py)
                    .call_method0("mmff_has_all_molecule_params")?
                    .extract::<bool>()?,
                expected_default
            );
            let default_result = molecule.bind(py).call_method0("mmff_properties")?;
            assert_eq!(
                default_result
                    .call_method0("variant")?
                    .extract::<MmffVariant>()? as isize,
                MmffVariant::Mmff94 as isize
            );
            let mut retained = Vec::new();
            let invalid = PyDict::new(py);
            invalid.set_item("mmff_variant", "source-defined-unknown")?;
            assert!(
                molecule
                    .bind(py)
                    .call_method("mmff_properties", (), Some(&invalid))
                    .unwrap_err()
                    .is_instance_of::<pyo3::exceptions::PyValueError>(py)
            );
            for variant in ["MMFF94", "MMFF94s"] {
                let params = Py::new(
                    py,
                    MmffPropertiesParams {
                        inner: crate::parameters::NativeField::Owned(ck::MmffPropertiesParams {
                            mmff_variant: ck::MmffVariant::from_rdkit_constructor_arg(variant),
                        }),
                    },
                )?;
                let expected = molecule
                    .borrow(py)
                    .inner
                    .mmff_properties_with_params(&params.borrow(py).to_rust(py)?)
                    .unwrap();
                let keywords = PyDict::new(py);
                keywords.set_item("mmff_variant", variant)?;
                let values = [
                    molecule
                        .bind(py)
                        .call_method1("mmff_properties", (params.clone_ref(py),))?,
                    molecule
                        .bind(py)
                        .call_method("mmff_properties", (), Some(&keywords))?,
                    molecule
                        .bind(py)
                        .call_method1("mmff_properties_with_params", (params.clone_ref(py),))?,
                ];
                for value in values {
                    assert_eq!(
                        value.call_method0("is_valid")?.extract::<bool>()?,
                        expected.is_valid()
                    );
                    assert!(
                        value.call_method0("variant")?.extract::<MmffVariant>()?
                            == MmffVariant::from(expected.variant())
                    );
                    let rows = value.call_method0("atoms")?;
                    assert_eq!(rows.len()?, expected.atoms().len());
                    for (index, atom) in expected.atoms().iter().enumerate() {
                        assert_eq!(
                            value.call_method1("atom_type", (index,))?.extract::<u8>()?,
                            atom.atom_type
                        );
                        assert_eq!(
                            value
                                .call_method1("formal_charge", (index,))?
                                .extract::<f64>()?
                                .to_bits(),
                            atom.formal_charge.to_bits()
                        );
                        assert_eq!(
                            value
                                .call_method1("partial_charge", (index,))?
                                .extract::<f64>()?
                                .to_bits(),
                            atom.partial_charge.to_bits()
                        );
                        let row = rows
                            .get_item(index)?
                            .extract::<PyRef<MmffAtomProperties>>()?;
                        assert_eq!(row.inner, *atom);
                    }
                    let error = value.call_method1("atom_type", (usize::MAX,)).unwrap_err();
                    assert!(error.is_instance_of::<MmffMolPropertiesError>(py));
                    assert_eq!(
                        error.value(py).getattr("atom_index")?.extract::<usize>()?,
                        usize::MAX
                    );
                    assert_eq!(
                        error.value(py).getattr("atoms")?.extract::<usize>()?,
                        expected.atoms().len()
                    );
                    assert!(error.value(py).getattr("source")?.is_none());
                    assert!(error.value(py).setattr("atom_index", 0).is_err());
                    retained.push(value);
                }
                assert!(
                    molecule
                        .bind(py)
                        .call_method("mmff_properties", (params,), Some(&keywords))
                        .unwrap_err()
                        .is_instance_of::<pyo3::exceptions::PyTypeError>(py)
                );
            }
            assert_eq!(molecule.borrow(py).inner.to_binary().unwrap(), before);
            drop(molecule);
            for value in retained {
                assert_eq!(value.call_method0("atoms")?.len()?, 3);
            }

            let native = ck::MmffMolPropertiesError::Params(ck::MmffParamError::Angle(
                ck::MmffAngleLookupError::MissingDefinition {
                    atom_type: u32::MAX,
                },
            ));
            assert!(native.source().is_none());
            let error = MmffMolPropertiesError::from_rust(native);
            assert!(error.value(py).getattr("source")?.is_none());
            let leaf = error.value(py).getattr("params")?.getattr("angle")?;
            assert!(leaf.is_instance_of::<MmffAngleLookupError>());
            assert_eq!(leaf.getattr("atom_type")?.extract::<u32>()?, u32::MAX);
            let cause = Py::new(
                py,
                MmffParamParseCause {
                    inner: ck::MmffParamParseCause::InvalidFloat {
                        cell: "1\0x".into(),
                    },
                },
            )?;
            assert!(
                cause
                    .bind(py)
                    .getattr("variant")?
                    .extract::<MmffParamParseCauseKind>()?
                    == MmffParamParseCauseKind::InvalidFloat
            );
            assert_eq!(cause.bind(py).getattr("cell")?.extract::<String>()?, "1\0x");
            assert!(cause.bind(py).setattr("cell", "1").is_err());
            Ok(())
        })
        .unwrap();
    }
}
