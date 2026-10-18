//! Concentrated native forwarding checks for the remaining typed atom producers.
use crate::native::Molecule;
use crate::operation_errors::{OperationError, OperationErrorKind};
use crate::ordered_results::AtomPairAtomCodeResult;
use crate::property_values::PropertyValue;
use cosmolkit as ck;
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
        assert!(
            source
                .bind(py)
                .call_method1("atom_property", (1, "tracking_id"))?
                .is_none()
        );
        for native_value in [
            ck::PropertyValue::String("tag\0text".into()),
            ck::PropertyValue::Int(i32::MIN),
            ck::PropertyValue::UInt(u32::MAX),
            ck::PropertyValue::Bool(true),
            ck::PropertyValue::Double(-0.0),
            ck::PropertyValue::IntVector(vec![i32::MIN, i32::MAX]),
            ck::PropertyValue::StringVector(vec!["a\0b".into(), "".into()]),
        ] {
            let expected = original
                .with_atom_property(atom, "tracking_id", &native_value)
                .unwrap();
            let value = Py::new(
                py,
                PropertyValue {
                    inner: native_value.clone(),
                },
            )?;
            let output = source
                .bind(py)
                .call_method1("with_atom_property", (1, "tracking_id", &value))?
                .extract::<Py<Molecule>>()?;
            assert_eq!(
                output.borrow(py).inner.to_binary().unwrap(),
                expected.to_binary().unwrap()
            );
            let read = output
                .bind(py)
                .call_method1("atom_property", (1, "tracking_id"))?
                .extract::<Py<PropertyValue>>()?;
            assert_eq!(read.borrow(py).inner, native_value);
            if let ck::PropertyValue::Double(want) = native_value {
                assert_eq!(
                    read.borrow(py).inner.as_double().unwrap().to_bits(),
                    want.to_bits()
                );
            }
            assert_eq!(source.borrow(py).inner.to_binary().unwrap(), before);
            let mutable = Py::new(
                py,
                Molecule {
                    inner: original.clone(),
                },
            )?;
            assert!(
                mutable
                    .bind(py)
                    .call_method1("set_atom_property_", (1, "tracking_id", &value))?
                    .is_none()
            );
            assert_eq!(
                mutable.borrow(py).inner.to_binary().unwrap(),
                expected.to_binary().unwrap()
            );
            drop(mutable);
            drop(output);
            // Reading returns an independent detached property value.
            assert_eq!(read.borrow(py).inner.kind(), value.borrow(py).inner.kind());
        }
        let value = Py::new(
            py,
            PropertyValue {
                inner: ck::PropertyValue::Int(42),
            },
        )?;
        for (method, index, key) in [
            ("with_atom_property", 1, "_CIPCode"),
            ("set_atom_property_", 1, "_CIPCode"),
            ("with_atom_property", usize::MAX, "tracking_id"),
            ("set_atom_property_", usize::MAX, "tracking_id"),
        ] {
            let native = original
                .with_atom_property(ck::AtomId::new(index), key, &value.borrow(py).inner)
                .unwrap_err();
            let error = source
                .bind(py)
                .call_method1(method, (index, key, &value))
                .unwrap_err();
            assert!(error.is_instance_of::<OperationError>(py));
            assert_eq!(error.value(py).str()?.to_str()?, native.to_string());
            assert_eq!(source.borrow(py).inner.to_binary().unwrap(), before);
        }
        let error = source
            .bind(py)
            .call_method1("atom_property", (usize::MAX, "tracking_id"))
            .unwrap_err();
        assert!(
            error
                .value(py)
                .getattr("variant")?
                .extract::<OperationErrorKind>()?
                == OperationErrorKind::AtomPropertyIndex
        );
        assert_eq!(
            error.value(py).getattr("atom_count")?.extract::<usize>()?,
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
