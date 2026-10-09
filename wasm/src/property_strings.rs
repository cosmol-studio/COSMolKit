//! Property string queries delegate source formatting to the canonical owner.
use crate::Molecule;
use cosmolkit as ck;
impl Molecule {
    pub fn atom_property(
        &self,
        atom: usize,
        key: &str,
    ) -> Result<Option<ck::PropertyValue>, ck::OperationError> {
        self.inner
            .borrow()
            .atom_property(ck::AtomId::new(atom), key)
            .map(|v| v.cloned())
    }
    #[cfg(feature = "cap-transforms")]
    pub fn with_atom_property(
        &self,
        atom: usize,
        key: &str,
        value: &ck::PropertyValue,
    ) -> Result<Self, ck::OperationError> {
        self.inner
            .borrow()
            .with_atom_property(ck::AtomId::new(atom), key, value)
            .map(|inner| Self {
                inner: inner.into(),
            })
    }
    #[cfg(feature = "cap-transforms")]
    pub fn set_atom_property_(
        &self,
        atom: usize,
        key: &str,
        value: &ck::PropertyValue,
    ) -> Result<(), ck::OperationError> {
        self.inner
            .borrow_mut()
            .set_atom_property_(ck::AtomId::new(atom), key, value)
    }
    pub fn atom_property_string(
        &self,
        id: usize,
        key: &str,
    ) -> Result<Option<ck::PropertyText>, ck::PropertyStringError> {
        // COSMolKit❗✔️: self.inner.atom_property_string(ck::AtomId::new(id),key)
        self.inner
            .borrow()
            .atom_property_string(ck::AtomId::new(id), key)
    }
    pub fn bond_property_string(
        &self,
        id: usize,
        key: &str,
    ) -> Result<Option<ck::PropertyText>, ck::PropertyStringError> {
        // COSMolKit❗✔️: self.inner.bond_property_string(ck::BondId::new(id),key)
        self.inner
            .borrow()
            .bond_property_string(ck::BondId::new(id), key)
    }
}

#[cfg(test)]
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
                m.atom_property_string(0, &key).unwrap(),
                Some((*expected).into())
            );
            assert_eq!(
                m.bond_property_string(0, &key).unwrap(),
                Some((*expected).into())
            );
            assert_eq!(
                m.atom_property_string(0, &key),
                original.atom_property_string(ck::AtomId::new(0), &key)
            );
        }
        for id in [0, 1, usize::MAX] {
            assert_eq!(m.atom_property_string(id, "missing"), Ok(None));
            assert_eq!(m.bond_property_string(id, "missing"), Ok(None));
        }
        assert_eq!(*m.inner.borrow(), original);
        let e = ck::PropertyStringError::UnsupportedKind {
            kind: ck::PropertyValueKind::UInt,
        };
        assert_eq!(e.kind(), ck::PropertyValueKind::UInt);
    }
}
