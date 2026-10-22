//! Typed atom producers forward to the native owner without changing COW.
use crate::native::Molecule;
use crate::operation_errors::{OperationError, OperationErrorKind};
use crate::ordered_results::AtomPairAtomCodeResult;
use crate::property_values::AtomPropertyReadError;
use cosmolkit as ck;
use pyo3::IntoPyObjectExt;
use pyo3::prelude::*;

#[test]
fn atom_properties_and_code_preserve_native_types_value_inplace_cow_and_errors() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let original = ck::Molecule::from_smiles("C[C@H](F)Cl").unwrap();
        let before = original.to_binary().unwrap();
        let atom = ck::AtomId::new(1);
        let source = Py::new(
            py,
            Molecule {
                inner: original.clone(),
            },
        )?;
        for (suffix, native_value, value) in [
            (
                "string",
                ck::PropertyValue::String("tag\0text".into()),
                "tag\0text".into_bound_py_any(py)?,
            ),
            (
                "int",
                ck::PropertyValue::Int(i32::MIN),
                i32::MIN.into_bound_py_any(py)?,
            ),
            (
                "uint",
                ck::PropertyValue::UInt(u32::MAX),
                u32::MAX.into_bound_py_any(py)?,
            ),
            (
                "bool",
                ck::PropertyValue::Bool(true),
                true.into_bound_py_any(py)?,
            ),
            (
                "double",
                ck::PropertyValue::Double(-0.0),
                (-0.0_f64).into_bound_py_any(py)?,
            ),
            (
                "int_vector",
                ck::PropertyValue::IntVector(vec![i32::MIN, i32::MAX]),
                vec![i32::MIN, i32::MAX].into_bound_py_any(py)?,
            ),
            (
                "string_vector",
                ck::PropertyValue::StringVector(vec!["a\0b".into(), "".into()]),
                vec!["a\0b", ""].into_bound_py_any(py)?,
            ),
        ] {
            let expected = original
                .with_atom_property(atom, "tracking_id", &native_value)
                .unwrap();
            let output = source
                .bind(py)
                .call_method1(
                    format!("with_atom_property_{suffix}"),
                    (1, "tracking_id", &value),
                )?
                .extract::<Py<Molecule>>()?;
            assert_eq!(
                output.borrow(py).inner.to_binary().unwrap(),
                expected.to_binary().unwrap()
            );
            let read = output
                .bind(py)
                .call_method1(format!("atom_property_{suffix}"), (1, "tracking_id"))?;
            assert!(read.eq(&value)?);
            if suffix == "double" {
                assert_eq!(read.extract::<f64>()?.to_bits(), (-0.0_f64).to_bits());
            }
            assert!(
                source
                    .bind(py)
                    .call_method1(format!("atom_property_{suffix}"), (1, "tracking_id"))?
                    .is_none()
            );
            let mutable = Py::new(
                py,
                Molecule {
                    inner: original.clone(),
                },
            )?;
            assert!(
                mutable
                    .bind(py)
                    .call_method1(
                        format!("set_atom_property_{suffix}_"),
                        (1, "tracking_id", &value)
                    )?
                    .is_none()
            );
            assert_eq!(
                mutable.borrow(py).inner.to_binary().unwrap(),
                expected.to_binary().unwrap()
            );
            assert_eq!(source.borrow(py).inner.to_binary().unwrap(), before);
        }
        for (method, index, key) in [
            ("with_atom_property_int", 1, "_CIPCode"),
            ("set_atom_property_int_", 1, "_CIPCode"),
            ("with_atom_property_int", usize::MAX, "tracking_id"),
            ("set_atom_property_int_", usize::MAX, "tracking_id"),
        ] {
            let native = original
                .with_atom_property_int(ck::AtomId::new(index), key, 42)
                .unwrap_err();
            let error = source
                .bind(py)
                .call_method1(method, (index, key, 42))
                .unwrap_err();
            assert!(error.is_instance_of::<OperationError>(py));
            assert_eq!(error.value(py).str()?.to_str()?, native.to_string());
            assert_eq!(source.borrow(py).inner.to_binary().unwrap(), before);
        }
        let error = source
            .bind(py)
            .call_method1("atom_property_int", (usize::MAX, "tracking_id"))
            .unwrap_err();
        assert!(error.is_instance_of::<AtomPropertyReadError>(py));
        let cause = error.cause(py).unwrap();
        assert!(
            cause
                .value(py)
                .getattr("variant")?
                .extract::<OperationErrorKind>()?
                == OperationErrorKind::AtomPropertyIndex
        );
        assert_eq!(
            cause.value(py).getattr("atom_count")?.extract::<usize>()?,
            original.num_atoms()
        );
        for (branch, chirality, legacy) in [(0, false, true), (0, true, false), (1, true, true)] {
            let want = original
                .with_atom_pair_atom_code(atom, branch, chirality, legacy)
                .unwrap();
            let got = source
                .bind(py)
                .call_method1("with_atom_pair_atom_code", (1, branch, chirality, legacy))?
                .extract::<Py<AtomPairAtomCodeResult>>()?;
            assert_eq!(got.borrow(py).inner.code, want.code);
            assert_eq!(
                got.borrow(py).inner.molecule.to_binary().unwrap(),
                want.molecule.to_binary().unwrap()
            );
            assert_eq!(source.borrow(py).inner.to_binary().unwrap(), before);
        }
        let want = original
            .with_atom_pair_atom_code(atom, 0, false, true)
            .unwrap();
        let got = source
            .bind(py)
            .call_method1("with_atom_pair_atom_code", (1,))?
            .extract::<Py<AtomPairAtomCodeResult>>()?;
        assert_eq!(got.borrow(py).inner, want);
        assert_eq!(original.to_binary().unwrap(), before);
        Ok(())
    })
    .unwrap();
}
