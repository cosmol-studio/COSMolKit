//! Exact-type atom properties delegate to the Rust facade.
use crate::Molecule;
use cosmolkit as ck;
impl Molecule {
    pub fn atom_property_int(
        &self,
        atom: usize,
        key: &str,
    ) -> Result<Option<i32>, ck::AtomPropertyReadError> {
        self.inner
            .borrow()
            .atom_property_int(ck::AtomId::new(atom), key)
    }
    #[cfg(feature = "cap-transforms")]
    pub fn with_atom_property_int(
        &self,
        atom: usize,
        key: &str,
        value: i32,
    ) -> Result<Self, ck::OperationError> {
        self.inner
            .borrow()
            .with_atom_property_int(ck::AtomId::new(atom), key, value)
            .map(|inner| Self {
                inner: inner.into(),
            })
    }
    #[cfg(feature = "cap-transforms")]
    pub fn set_atom_property_int_(
        &self,
        atom: usize,
        key: &str,
        value: i32,
    ) -> Result<(), ck::OperationError> {
        self.inner
            .borrow_mut()
            .set_atom_property_int_(ck::AtomId::new(atom), key, value)
    }
    pub fn atom_property_uint(
        &self,
        atom: usize,
        key: &str,
    ) -> Result<Option<u32>, ck::AtomPropertyReadError> {
        self.inner
            .borrow()
            .atom_property_uint(ck::AtomId::new(atom), key)
    }
    #[cfg(feature = "cap-transforms")]
    pub fn with_atom_property_uint(
        &self,
        atom: usize,
        key: &str,
        value: u32,
    ) -> Result<Self, ck::OperationError> {
        self.inner
            .borrow()
            .with_atom_property_uint(ck::AtomId::new(atom), key, value)
            .map(|inner| Self {
                inner: inner.into(),
            })
    }
    #[cfg(feature = "cap-transforms")]
    pub fn set_atom_property_uint_(
        &self,
        atom: usize,
        key: &str,
        value: u32,
    ) -> Result<(), ck::OperationError> {
        self.inner
            .borrow_mut()
            .set_atom_property_uint_(ck::AtomId::new(atom), key, value)
    }
    pub fn atom_property_double(
        &self,
        atom: usize,
        key: &str,
    ) -> Result<Option<f64>, ck::AtomPropertyReadError> {
        self.inner
            .borrow()
            .atom_property_double(ck::AtomId::new(atom), key)
    }
    #[cfg(feature = "cap-transforms")]
    pub fn with_atom_property_double(
        &self,
        atom: usize,
        key: &str,
        value: f64,
    ) -> Result<Self, ck::OperationError> {
        self.inner
            .borrow()
            .with_atom_property_double(ck::AtomId::new(atom), key, value)
            .map(|inner| Self {
                inner: inner.into(),
            })
    }
    #[cfg(feature = "cap-transforms")]
    pub fn set_atom_property_double_(
        &self,
        atom: usize,
        key: &str,
        value: f64,
    ) -> Result<(), ck::OperationError> {
        self.inner
            .borrow_mut()
            .set_atom_property_double_(ck::AtomId::new(atom), key, value)
    }
    pub fn atom_property_bool(
        &self,
        atom: usize,
        key: &str,
    ) -> Result<Option<bool>, ck::AtomPropertyReadError> {
        self.inner
            .borrow()
            .atom_property_bool(ck::AtomId::new(atom), key)
    }
    #[cfg(feature = "cap-transforms")]
    pub fn with_atom_property_bool(
        &self,
        atom: usize,
        key: &str,
        value: bool,
    ) -> Result<Self, ck::OperationError> {
        self.inner
            .borrow()
            .with_atom_property_bool(ck::AtomId::new(atom), key, value)
            .map(|inner| Self {
                inner: inner.into(),
            })
    }
    #[cfg(feature = "cap-transforms")]
    pub fn set_atom_property_bool_(
        &self,
        atom: usize,
        key: &str,
        value: bool,
    ) -> Result<(), ck::OperationError> {
        self.inner
            .borrow_mut()
            .set_atom_property_bool_(ck::AtomId::new(atom), key, value)
    }
    pub fn atom_property_string(
        &self,
        atom: usize,
        key: &str,
    ) -> Result<Option<String>, ck::AtomPropertyReadError> {
        self.inner
            .borrow()
            .atom_property_string(ck::AtomId::new(atom), key)
    }
    #[cfg(feature = "cap-transforms")]
    pub fn with_atom_property_string(
        &self,
        atom: usize,
        key: &str,
        value: &str,
    ) -> Result<Self, ck::OperationError> {
        self.inner
            .borrow()
            .with_atom_property_string(ck::AtomId::new(atom), key, value)
            .map(|inner| Self {
                inner: inner.into(),
            })
    }
    #[cfg(feature = "cap-transforms")]
    pub fn set_atom_property_string_(
        &self,
        atom: usize,
        key: &str,
        value: &str,
    ) -> Result<(), ck::OperationError> {
        self.inner
            .borrow_mut()
            .set_atom_property_string_(ck::AtomId::new(atom), key, value)
    }
    pub fn atom_property_int_vector(
        &self,
        atom: usize,
        key: &str,
    ) -> Result<Option<Vec<i32>>, ck::AtomPropertyReadError> {
        self.inner
            .borrow()
            .atom_property_int_vector(ck::AtomId::new(atom), key)
    }
    #[cfg(feature = "cap-transforms")]
    pub fn with_atom_property_int_vector(
        &self,
        atom: usize,
        key: &str,
        value: &[i32],
    ) -> Result<Self, ck::OperationError> {
        self.inner
            .borrow()
            .with_atom_property_int_vector(ck::AtomId::new(atom), key, value)
            .map(|inner| Self {
                inner: inner.into(),
            })
    }
    #[cfg(feature = "cap-transforms")]
    pub fn set_atom_property_int_vector_(
        &self,
        atom: usize,
        key: &str,
        value: &[i32],
    ) -> Result<(), ck::OperationError> {
        self.inner
            .borrow_mut()
            .set_atom_property_int_vector_(ck::AtomId::new(atom), key, value)
    }
    pub fn atom_property_string_vector(
        &self,
        atom: usize,
        key: &str,
    ) -> Result<Option<Vec<String>>, ck::AtomPropertyReadError> {
        self.inner
            .borrow()
            .atom_property_string_vector(ck::AtomId::new(atom), key)
    }
    #[cfg(feature = "cap-transforms")]
    pub fn with_atom_property_string_vector(
        &self,
        atom: usize,
        key: &str,
        value: &[String],
    ) -> Result<Self, ck::OperationError> {
        self.inner
            .borrow()
            .with_atom_property_string_vector(ck::AtomId::new(atom), key, value)
            .map(|inner| Self {
                inner: inner.into(),
            })
    }
    #[cfg(feature = "cap-transforms")]
    pub fn set_atom_property_string_vector_(
        &self,
        atom: usize,
        key: &str,
        value: &[String],
    ) -> Result<(), ck::OperationError> {
        self.inner
            .borrow_mut()
            .set_atom_property_string_vector_(ck::AtomId::new(atom), key, value)
    }
    #[cfg(feature = "cap-io")]
    pub fn bond_property_string(
        &self,
        id: usize,
        key: &str,
    ) -> Result<Option<ck::PropertyText>, ck::PropertyStringError> {
        self.inner
            .borrow()
            .bond_property_string(ck::BondId::new(id), key)
    }
}
#[cfg(all(test, feature = "cap-io"))]
mod tests {
    use super::*;
    #[test]
    fn property_strings_all_six_value_kinds_are_exact_read_only_and_absence_is_none() {
        let cases = [
            (ck::PropertyValue::String("a\0b".into()), "a\0b"),
            (ck::PropertyValue::Int(i32::MIN), "-2147483648"),
            (ck::PropertyValue::UInt(u32::MAX), "4294967295"),
            (
                ck::PropertyValue::IntVector(vec![i32::MIN, 0, i32::MAX]),
                "[-2147483648,0,2147483647]",
            ),
            (ck::PropertyValue::Double(0.1), "0.10000000000000001"),
            (ck::PropertyValue::Double(-0.0), "-0"),
            (ck::PropertyValue::Bool(true), "1"),
            (ck::PropertyValue::Bool(false), "0"),
        ];
        let mut atom = ck::AtomSpec::new(ck::Element::C);
        let mut bond = ck::BondSpec::new(
            ck::AtomId::new(0),
            ck::AtomId::new(1),
            ck::BondOrder::Single,
        );
        for (i, (value, _)) in cases.iter().enumerate() {
            atom = atom.with_prop(format!("v{i}"), value.clone()).unwrap();
            bond = bond.with_prop(format!("v{i}"), value.clone()).unwrap();
        }
        let mut builder = ck::MoleculeBuilder::new();
        builder.add_atom(atom);
        builder.add_atom(ck::AtomSpec::new(ck::Element::N));
        builder.add_bond(bond).unwrap();
        let m = Molecule {
            inner: builder.build().unwrap().into(),
        };
        let original = m.inner.borrow().clone();
        for (i, (_, expected)) in cases.iter().enumerate() {
            let key = format!("v{i}");
            assert_eq!(
                m.inner
                    .borrow()
                    .atom_property(ck::AtomId::new(0), &key)
                    .unwrap(),
                Some(&cases[i].0)
            );
            assert_eq!(
                m.bond_property_string(0, &key).unwrap(),
                Some((*expected).into())
            );
            if i == 0 {
                assert_eq!(
                    m.atom_property_string(0, &key).unwrap(),
                    Some((*expected).into())
                );
            } else {
                assert!(m.atom_property_string(0, &key).is_err());
            }
        }
        for id in [0, 1, usize::MAX] {
            if id < 2 {
                assert_eq!(m.atom_property_string(id, "missing"), Ok(None));
            } else {
                assert!(m.atom_property_string(id, "missing").is_err());
            }
            assert_eq!(m.bond_property_string(id, "missing"), Ok(None));
        }
        assert_eq!(*m.inner.borrow(), original);
        let e = ck::PropertyStringError::UnsupportedKind {
            kind: ck::PropertyValueKind::UInt,
        };
        assert_eq!(e.kind(), ck::PropertyValueKind::UInt);
    }
}
