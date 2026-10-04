//! Canonical detached atom and bond property values.

use std::collections::{BTreeMap, BTreeSet};

/// The modeled source value kinds supported by atom and bond properties.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PropertyValueKind {
    String,
    Int,
    Double,
    Bool,
}

/// A canonical detached atom or bond property value.
#[derive(Debug, Clone)]
pub enum PropertyValue {
    String(String),
    Int(i32),
    Double(f64),
    Bool(bool),
}

impl PartialEq for PropertyValue {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::String(left), Self::String(right)) => left == right,
            (Self::Int(left), Self::Int(right)) => left == right,
            (Self::Double(left), Self::Double(right)) => left.to_bits() == right.to_bits(),
            (Self::Bool(left), Self::Bool(right)) => left == right,
            _ => false,
        }
    }
}

impl Eq for PropertyValue {}

impl From<String> for PropertyValue {
    fn from(value: String) -> Self {
        Self::String(value)
    }
}

impl From<&str> for PropertyValue {
    fn from(value: &str) -> Self {
        Self::String(value.to_owned())
    }
}

impl From<&String> for PropertyValue {
    fn from(value: &String) -> Self {
        Self::String(value.clone())
    }
}

impl From<&PropertyValue> for PropertyValue {
    fn from(value: &PropertyValue) -> Self {
        value.clone()
    }
}

impl From<i32> for PropertyValue {
    fn from(value: i32) -> Self {
        Self::Int(value)
    }
}

impl From<f64> for PropertyValue {
    fn from(value: f64) -> Self {
        Self::Double(value)
    }
}

impl From<bool> for PropertyValue {
    fn from(value: bool) -> Self {
        Self::Bool(value)
    }
}

/// A property value was read using an accessor for a different value kind.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[error("property value has kind {actual:?}, expected {expected:?}")]
pub struct PropertyValueError {
    expected: PropertyValueKind,
    actual: PropertyValueKind,
}

impl PropertyValueError {
    #[must_use]
    pub const fn expected(&self) -> PropertyValueKind {
        self.expected
    }

    #[must_use]
    pub const fn actual(&self) -> PropertyValueKind {
        self.actual
    }
}

impl PropertyValue {
    #[must_use]
    pub const fn kind(&self) -> PropertyValueKind {
        // BEGIN RDKIT CPP FUNCTION RDValue::getTag / RDTypeTag modeled subset
        // RDKit✔️✔️: case RDTypeTag::StringTag:
        // RDKit✔️✔️: case RDTypeTag::IntTag:
        // RDKit✔️✔️: case RDTypeTag::DoubleTag: {
        // RDKit✔️✔️: case RDTypeTag::BoolTag:
        // END RDKIT CPP FUNCTION RDValue::getTag / RDTypeTag modeled subset
        // Behavior review: each modeled detached variant has one source tag;
        // no text inspection or numeric coercion changes the stored kind.
        // Complexity review: one enum discriminant match is constant time and
        // allocation free, equivalent to the source tag switch.
        match self {
            Self::String(_) => PropertyValueKind::String,
            Self::Int(_) => PropertyValueKind::Int,
            Self::Double(_) => PropertyValueKind::Double,
            Self::Bool(_) => PropertyValueKind::Bool,
        }
    }

    fn wrong_kind(&self, expected: PropertyValueKind) -> PropertyValueError {
        PropertyValueError {
            expected,
            actual: self.kind(),
        }
    }

    pub fn as_string(&self) -> Result<&str, PropertyValueError> {
        match self {
            Self::String(value) => Ok(value),
            _ => Err(self.wrong_kind(PropertyValueKind::String)),
        }
    }

    pub fn as_int(&self) -> Result<i32, PropertyValueError> {
        match self {
            Self::Int(value) => Ok(*value),
            _ => Err(self.wrong_kind(PropertyValueKind::Int)),
        }
    }

    pub fn as_double(&self) -> Result<f64, PropertyValueError> {
        match self {
            Self::Double(value) => Ok(*value),
            _ => Err(self.wrong_kind(PropertyValueKind::Double)),
        }
    }

    pub fn as_bool(&self) -> Result<bool, PropertyValueError> {
        match self {
            Self::Bool(value) => Ok(*value),
            _ => Err(self.wrong_kind(PropertyValueKind::Bool)),
        }
    }
}

/// One canonical typed value map with source insertion order and computed state.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct PropertyStore {
    values: BTreeMap<String, PropertyValue>,
    order: Vec<String>,
    computed: BTreeSet<String>,
}

impl PropertyStore {
    pub(crate) const fn new() -> Self {
        Self {
            values: BTreeMap::new(),
            order: Vec::new(),
            computed: BTreeSet::new(),
        }
    }

    pub(crate) fn values(&self) -> &BTreeMap<String, PropertyValue> {
        &self.values
    }

    pub(crate) fn get(&self, key: &str) -> Option<&PropertyValue> {
        self.values.get(key)
    }

    pub(crate) fn computed_names(&self) -> &BTreeSet<String> {
        &self.computed
    }

    pub(crate) fn is_computed(&self, key: &str) -> bool {
        self.computed.contains(key)
    }

    pub(crate) fn ordered(&self) -> impl ExactSizeIterator<Item = (&str, &PropertyValue)> + '_ {
        self.order.iter().map(|key| {
            let value = self
                .values
                .get(key)
                .expect("private property order must match its canonical value map");
            (key.as_str(), value)
        })
    }

    pub(crate) fn set(&mut self, key: String, value: PropertyValue) {
        // BEGIN RDKIT CPP FUNCTION Dict::setVal
        // RDKit✔️🔝: for (auto &&data : _data) {
        // RDKit✔️🔝:   if (data.key == what) {
        // RDKit✔️🔝:     RDValue::cleanup_rdvalue(data.val);
        // RDKit✔️🔝:     data.val = val;
        // RDKit✔️🔝:     return;
        // RDKit✔️🔝:   }
        // RDKit✔️🔝: }
        // RDKit✔️🔝: _data.push_back(Pair(what, val));
        // END RDKIT CPP FUNCTION Dict::setVal
        // The tree replaces the source linear value lookup. The order vector
        // stores only keys, so overwrites preserve position without duplicating
        // any String/Int/Double/Bool value.
        if !self.values.contains_key(&key) {
            self.order.push(key.clone());
        }
        self.values.insert(key, value);
    }

    pub(crate) fn set_computed(&mut self, key: String, value: PropertyValue) {
        self.set(key.clone(), value);
        self.computed.insert(key);
    }

    pub(crate) fn clear(&mut self, key: &str) {
        // BEGIN RDKIT CPP FUNCTION Dict::clearVal
        // RDKit✔️🔝: for (auto it = _data.begin(); it < _data.end(); ++it) {
        // RDKit✔️🔝:   if (it->key == what) {
        // RDKit✔️🔝:     if (_hasNonPodData) {
        // RDKit✔️🔝:       RDValue::cleanup_rdvalue(it->val);
        // RDKit✔️🔝:     }
        // RDKit✔️🔝:     _data.erase(it);
        // RDKit✔️🔝:     return;
        // RDKit✔️🔝:   }
        // RDKit✔️🔝: }
        // END RDKIT CPP FUNCTION Dict::clearVal
        // Tree/set removal plus one linear order-key erase preserves the source
        // transition with no value cloning.
        if self.values.remove(key).is_some()
            && let Some(position) = self.order.iter().position(|name| name == key)
        {
            self.order.remove(position);
        }
        self.computed.remove(key);
    }

    pub(crate) fn clear_computed(&mut self) {
        for key in std::mem::take(&mut self.computed) {
            self.values.remove(&key);
        }
        self.order.retain(|key| self.values.contains_key(key));
    }

    #[cfg(test)]
    pub(crate) fn ordered_keys(&self) -> &[String] {
        &self.order
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Element};

    #[test]
    fn typed_property_value_preserves_all_variants_and_exact_access_errors() {
        let values = [
            PropertyValue::String("7".to_owned()),
            PropertyValue::Int(7),
            PropertyValue::Double(7.0),
            PropertyValue::Bool(true),
        ];
        assert_eq!(values[0].kind(), PropertyValueKind::String);
        assert_eq!(values[1].kind(), PropertyValueKind::Int);
        assert_eq!(values[2].kind(), PropertyValueKind::Double);
        assert_eq!(values[3].kind(), PropertyValueKind::Bool);
        assert_eq!(values[0].as_string(), Ok("7"));
        assert_eq!(values[1].as_int(), Ok(7));
        assert_eq!(values[2].as_double(), Ok(7.0));
        assert_eq!(values[3].as_bool(), Ok(true));
        assert_ne!(values[0], values[1]);
        assert_ne!(values[1], values[2]);
        assert_eq!(
            values[1].as_string(),
            Err(PropertyValueError {
                expected: PropertyValueKind::String,
                actual: PropertyValueKind::Int,
            })
        );
        assert_eq!(
            values[0].as_bool(),
            Err(PropertyValueError {
                expected: PropertyValueKind::Bool,
                actual: PropertyValueKind::String,
            })
        );
    }

    #[test]
    fn typed_property_value_double_equality_preserves_bits_and_signed_zero() {
        assert_ne!(PropertyValue::Double(0.0), PropertyValue::Double(-0.0));
        assert_eq!(
            PropertyValue::Double(f64::from_bits(0x7ff8_0000_0000_0042)),
            PropertyValue::Double(f64::from_bits(0x7ff8_0000_0000_0042))
        );
        assert_ne!(
            PropertyValue::Double(f64::from_bits(0x7ff8_0000_0000_0042)),
            PropertyValue::Double(f64::from_bits(0x7ff8_0000_0000_0043))
        );
    }

    #[test]
    fn typed_property_value_order_and_lifecycle_preserve_type_transitions() {
        let mut store = PropertyStore::new();
        store.set("z".to_owned(), PropertyValue::String("7".to_owned()));
        store.set("a".to_owned(), PropertyValue::Int(7));
        assert_eq!(store.ordered_keys(), &["z", "a"]);
        store.set("z".to_owned(), PropertyValue::Double(-0.0));
        assert_eq!(store.ordered_keys(), &["z", "a"]);
        assert_eq!(store.get("z"), Some(&PropertyValue::Double(-0.0)));

        store.clear("z");
        assert_eq!(store.ordered_keys(), &["a"]);
        store.set("z".to_owned(), PropertyValue::Bool(false));
        assert_eq!(store.ordered_keys(), &["a", "z"]);
        store.set_computed("c".to_owned(), PropertyValue::Double(1.25));
        store.set_computed("a".to_owned(), PropertyValue::String("seven".to_owned()));
        assert_eq!(store.ordered_keys(), &["a", "z", "c"]);
        assert!(store.is_computed("a"));
        assert!(store.is_computed("c"));
        store.clear_computed();
        assert_eq!(store.ordered_keys(), &["z"]);
        assert_eq!(store.get("z"), Some(&PropertyValue::Bool(false)));
        assert!(store.computed_names().is_empty());
    }

    #[test]
    fn typed_property_value_invalid_keys_are_failure_atomic_for_atom_and_bond() {
        let mut atom = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_prop("kept", PropertyValue::Int(7))
                .unwrap(),
        );
        let atom_before = atom.clone();
        assert!(atom.set_prop("", PropertyValue::Bool(true)).is_err());
        assert!(
            atom.set_computed_prop("", PropertyValue::Double(-0.0))
                .is_err()
        );
        assert_eq!(atom, atom_before);

        let mut bond = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                .with_prop("kept", PropertyValue::Bool(false))
                .unwrap(),
        );
        let bond_before = bond.clone();
        assert!(bond.set_prop("", PropertyValue::Int(4)).is_err());
        assert!(
            bond.set_computed_prop("", PropertyValue::String("bad".to_owned()))
                .is_err()
        );
        assert_eq!(bond, bond_before);
    }
}
