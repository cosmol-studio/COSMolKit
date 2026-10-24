//! Actual atom/bond property lexical queries; no formatting in the binding.
use crate::native::Molecule;
use cosmolkit as ck;
use pyo3::prelude::*;

#[test]
fn native_property_strings_preserve_all_kinds_absence_bytes_and_source() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let values = [
            ("string", ck::PropertyValue::String("a\0b\u{e9}".into())),
            ("int", ck::PropertyValue::Int(i32::MIN)),
            ("uint", ck::PropertyValue::UInt(u32::MAX)),
            ("ints", ck::PropertyValue::IntVector(vec![-7, 0, i32::MAX])),
            (
                "names\0sequence",
                ck::PropertyValue::StringVector(vec!["A\0B".into(), "\u{e9}".into()]),
            ),
            ("double", ck::PropertyValue::Double(-0.0)),
            ("bool", ck::PropertyValue::Bool(true)),
            (
                "bytes",
                ck::PropertyValue::String(ck::PropertyText::from_bytes(&[b'x', 0xff])),
            ),
        ];
        let mut atom = ck::AtomSpec::new(ck::Element::C);
        let mut bond = ck::BondSpec::new(
            ck::AtomId::new(0),
            ck::AtomId::new(1),
            ck::BondOrder::Single,
        );
        for (key, value) in &values {
            atom = atom.with_prop(*key, value.clone()).unwrap();
            bond = bond.with_prop(*key, value.clone()).unwrap();
        }
        let topology = ck::TopologyBlock::try_from_parts(
            vec![
                ck::Atom::from_spec(ck::AtomId::new(0), atom),
                ck::Atom::from_spec(ck::AtomId::new(1), ck::AtomSpec::new(ck::Element::C)),
            ],
            vec![ck::Bond::from_spec(ck::BondId::new(0), bond)],
            vec![],
            vec![],
        )
        .unwrap();
        let native =
            ck::Molecule::from_parts(topology, Default::default(), Default::default()).unwrap();
        let before = native.to_binary().unwrap();
        let molecule = Py::new(py, Molecule { inner: native })?;
        for method in ["bond_property_string"] {
            for (key, value) in &values {
                let expected = ck::property_value_to_text(value).unwrap();
                let actual = molecule.bind(py).call_method1(method, (0, *key));
                if *key == "bytes" {
                    let error = actual.unwrap_err();
                    assert!(error.is_instance_of::<pyo3::exceptions::PyUnicodeDecodeError>(py));
                    assert_eq!(
                        error.value(py).getattr("object")?.extract::<Vec<u8>>()?,
                        expected.as_bytes()
                    );
                } else {
                    assert_eq!(actual?.extract::<String>()?.as_bytes(), expected.as_bytes());
                }
            }
            for (index, key) in [(0, "missing"), (usize::MAX, "string")] {
                assert!(
                    molecule
                        .bind(py)
                        .call_method1(method, (index, key))?
                        .is_none()
                );
            }
            assert!(
                molecule
                    .bind(py)
                    .call_method1(method, (-1, "string"))
                    .unwrap_err()
                    .is_instance_of::<pyo3::exceptions::PyOverflowError>(py)
            );
        }
        assert_eq!(
            molecule
                .bind(py)
                .call_method1("atom_property_string", (0, "string"))?
                .extract::<String>()?,
            "a\0b\u{e9}"
        );
        let wrong = molecule
            .bind(py)
            .call_method1("atom_property_string", (0, "int"))
            .unwrap_err();
        assert!(wrong.is_instance(py, <crate::property_values::AtomPropertyReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
        assert!(
            wrong
                .cause(py)
                .unwrap()
                .is_instance(py, <crate::property_values::PropertyValueError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())
        );
        let invalid = molecule
            .bind(py)
            .call_method1("atom_property_string", (0, "bytes"))
            .unwrap_err();
        assert!(
            invalid
                .cause(py)
                .unwrap()
                .is_instance_of::<pyo3::exceptions::PyUnicodeError>(py)
        );
        assert_eq!(
            invalid
                .value(py)
                .getattr("valid_up_to")?
                .extract::<usize>()?,
            1
        );
        assert!(
            molecule
                .bind(py)
                .call_method1("atom_property_string", (0, "missing"))?
                .is_none()
        );
        assert!(
            molecule
                .bind(py)
                .call_method1("atom_property_string", (usize::MAX, "string"))
                .is_err()
        );
        assert_eq!(molecule.borrow(py).inner.to_binary().unwrap(), before);
        Ok(())
    })
    .unwrap();
}
