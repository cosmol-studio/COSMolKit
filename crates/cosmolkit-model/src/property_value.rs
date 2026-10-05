//! Canonical detached atom and bond property values.

use std::collections::{BTreeMap, BTreeSet};

/// The modeled source value kinds supported by atom and bond properties.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PropertyValueKind {
    String,
    Int,
    UInt,
    IntVector,
    Double,
    Bool,
}

/// A canonical detached atom or bond property value.
#[derive(Debug, Clone)]
pub enum PropertyValue {
    String(String),
    Int(i32),
    UInt(u32),
    IntVector(Vec<i32>),
    Double(f64),
    Bool(bool),
}

impl PartialEq for PropertyValue {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::String(left), Self::String(right)) => left == right,
            (Self::Int(left), Self::Int(right)) => left == right,
            (Self::UInt(left), Self::UInt(right)) => left == right,
            (Self::IntVector(left), Self::IntVector(right)) => left == right,
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

impl From<u32> for PropertyValue {
    fn from(value: u32) -> Self {
        // BEGIN RDKIT COMPLETE PROPOSED CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue-taggedunion.h:168-168
        // RDKit❗✔️:   inline Value(unsigned int v) : u(v) {}
        // END RDKIT COMPLETE PROPOSED CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue-taggedunion.h:168-168
        // BEGIN RDKIT COMPLETE PROPOSED CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue-taggedunion.h:227-227
        // RDKit❗✔️:   inline RDValue(unsigned v) : value(v), type(RDTypeTag::UnsignedIntTag) {}
        // END RDKIT COMPLETE PROPOSED CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue-taggedunion.h:227-227
        // Behavior: exact unsigned source constructor tag and full u32 POD value.
        // Complexity: constant-time discriminant construction, no allocation.
        Self::UInt(value)
    }
}

impl From<Vec<i32>> for PropertyValue {
    fn from(value: Vec<i32>) -> Self {
        Self::IntVector(value)
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
        // BEGIN RDKIT COMPLETE PROPOSED CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue-taggedunion.h:253-253
        // RDKit❗✔️:   short getTag() const { return type; }
        // END RDKIT COMPLETE PROPOSED CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue-taggedunion.h:253-253
        // BEGIN RDKIT COMPLETE PROPOSED CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue-taggedunion.h:112-115
        // RDKit❗✔️: template <>
        // RDKit❗✔️: inline short GetTag<unsigned int>() {
        // RDKit❗✔️:   return UnsignedIntTag;
        // RDKit❗✔️: }
        // END RDKIT COMPLETE PROPOSED CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue-taggedunion.h:112-115
        // RDKit❗✔️: const short UnsignedIntTag = 6;
        // UInt is proposed/unrun. Existing other source tag mappings unchanged.
        // Behavior review: each modeled detached variant has one source tag;
        // no text inspection or numeric coercion changes the stored kind.
        // Complexity review: one enum discriminant match is constant time and
        // allocation free, equivalent to the source tag switch.
        match self {
            Self::String(_) => PropertyValueKind::String,
            Self::Int(_) => PropertyValueKind::Int,
            Self::UInt(_) => PropertyValueKind::UInt,
            Self::IntVector(_) => PropertyValueKind::IntVector,
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

    /// Read only the exact unsigned tag; this never coerces signed values.
    pub fn as_uint(&self) -> Result<u32, PropertyValueError> {
        // BEGIN RDKIT COMPLETE PROPOSED CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue-taggedunion.h:112-115
        // RDKit❗✔️: template <>
        // RDKit❗✔️: inline short GetTag<unsigned int>() {
        // RDKit❗✔️:   return UnsignedIntTag;
        // RDKit❗✔️: }
        // END RDKIT COMPLETE PROPOSED CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue-taggedunion.h:112-115
        // BEGIN RDKIT COMPLETE PROPOSED CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue-taggedunion.h:365-382
        // RDKit❗✔️: template <class T>
        // RDKit❗✔️: inline bool rdvalue_is(RDValue_cast_t v) {
        // RDKit❗✔️:   const short tag =
        // RDKit❗✔️:       RDTypeTag::GetTag<typename boost::remove_reference<T>::type>();
        // RDKit❗✔️:
        // RDKit❗✔️:   // If we are an Any tag, check the any type info
        // RDKit❗✔️:   //  see the template specialization below if we are
        // RDKit❗✔️:   //  looking for a boost any directly
        // RDKit❗✔️:   if (v.getTag() == RDTypeTag::AnyTag) {
        // RDKit❗✔️:     return v.value.a->type() == typeid(T);
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   if (v.getTag() == tag) {
        // RDKit❗✔️:     return true;
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   return false;
        // RDKit❗✔️: }
        // END RDKIT COMPLETE PROPOSED CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue-taggedunion.h:365-382
        // BEGIN RDKIT COMPLETE PROPOSED CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue-taggedunion.h:488-497
        // RDKit❗✔️: template <>
        // RDKit❗✔️: inline unsigned int rdvalue_cast<unsigned int>(RDValue_cast_t v) {
        // RDKit❗✔️:   if (rdvalue_is<unsigned int>(v)) {
        // RDKit❗✔️:     return v.value.u;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   if (rdvalue_is<int>(v)) {
        // RDKit❗✔️:     return boost::numeric_cast<unsigned int>(v.value.i);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   throw std::bad_any_cast();
        // RDKit❗✔️: }
        // END RDKIT COMPLETE PROPOSED CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue-taggedunion.h:488-497
        // Canonical policy: dev/public_api_design.md exact-kind accessors reject
        // every non-UInt variant. This uses source unsigned tag identity and its
        // direct POD read; it does not expose the source signed numeric_cast
        // coercion. Coercive reads remain in the owning source algorithm helpers.
        // Any type introspection is independently unmodeled; no Any variant here.
        // Behavior status: proposal only; no Rust execution or acceptance.
        // Complexity: one tag branch and copy; constant time, no allocation.
        match self {
            Self::UInt(value) => Ok(*value),
            _ => Err(self.wrong_kind(PropertyValueKind::UInt)),
        }
    }

    /// Borrow the canonical signed integer vector without conversion.
    pub fn as_int_vector(&self) -> Result<&[i32], PropertyValueError> {
        // RDKit✔️✔️: typedef std::vector<int> INT_VECT;
        // RDKit✔️✔️:   return rdvalue_cast<T>(arg);
        // Behavior: exact vector tag; borrowed elements retain order and duplicates.
        // Complexity: constant-time tag check and borrow, with no allocation.
        match self {
            Self::IntVector(value) => Ok(value),
            _ => Err(self.wrong_kind(PropertyValueKind::IntVector)),
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
#[cfg(test)]
mod q01_b1_tests {
    use super::*;
    use crate::{
        Atom, AtomId, AtomQueryPredicate, AtomSpec, Bond, BondId, BondOrder, BondQueryPredicate,
        BondSpec, Element, QueryAtom, QueryBond, QueryNode,
    };

    #[test]
    fn q01_b1_int_vector_strict_borrow_equality_and_width() {
        for values in [vec![], vec![1], vec![1, -2, 1], vec![i32::MIN, i32::MAX]] {
            let value = PropertyValue::from(values.clone());
            assert_eq!(value.kind(), PropertyValueKind::IntVector);
            let PropertyValue::IntVector(stored) = &value else {
                panic!("canonical variant");
            };
            assert_eq!(value.as_int_vector().unwrap(), values);
            assert_eq!(value.as_int_vector().unwrap().as_ptr(), stored.as_ptr());
            assert_eq!(value.clone(), value);
            let errors = [
                value.as_string().unwrap_err(),
                value.as_int().unwrap_err(),
                value.as_double().unwrap_err(),
                value.as_bool().unwrap_err(),
            ];
            for (error, expected) in errors.into_iter().zip([
                PropertyValueKind::String,
                PropertyValueKind::Int,
                PropertyValueKind::Double,
                PropertyValueKind::Bool,
            ]) {
                assert_eq!(error.expected(), expected);
                assert_eq!(error.actual(), PropertyValueKind::IntVector);
            }
        }
        for value in [
            PropertyValue::String("[1]".into()),
            PropertyValue::Int(1),
            PropertyValue::Double(1.0),
            PropertyValue::Bool(true),
        ] {
            let error = value.as_int_vector().unwrap_err();
            assert_eq!(error.expected(), PropertyValueKind::IntVector);
            assert_eq!(error.actual(), value.kind());
        }
        assert_ne!(
            PropertyValue::from(vec![1, -2, 1]),
            PropertyValue::from(vec![1, 1, -2])
        );
        let original = PropertyValue::from(vec![1, 1]);
        let mut cloned = original.clone();
        let PropertyValue::IntVector(v) = &mut cloned else {
            panic!()
        };
        v[0] = -1;
        assert_eq!(original.as_int_vector().unwrap(), [1, 1]);
        assert_ne!(cloned, original);
    }

    #[test]
    fn q01_b1_int_vector_store_lifecycle_is_typed_and_atomic() {
        let mut store = PropertyStore::new();
        store.set("a".into(), vec![1, -2, 1].into());
        store.set_computed("b".into(), vec![i32::MIN].into());
        store.set("b".into(), vec![i32::MAX].into());
        assert!(store.is_computed("b"));
        assert_eq!(store.ordered_keys(), ["a", "b"]);
        store.clear("a");
        store.set("a".into(), vec![].into());
        assert_eq!(store.ordered_keys(), ["b", "a"]);
        store.clear_computed();
        assert_eq!(store.ordered_keys(), ["a"]);
        assert_eq!(store.get("a").unwrap().as_int_vector().unwrap(), []);
        let mut atom = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_prop("v", vec![1, -2, 1])
                .unwrap(),
        );
        let before = atom.clone();
        assert!(atom.set_prop("", vec![1]).is_err());
        assert_eq!(atom, before);
        atom.set_computed_prop("v", vec![7]).unwrap();
        assert_eq!(
            before.prop("v").unwrap().as_int_vector().unwrap(),
            [1, -2, 1]
        );
        let qa = QueryAtom::from_carrier_parts(
            before.clone(),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
        );
        assert_eq!(qa.prop("v"), before.prop("v"));
        let mut bond = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                .with_prop("v", vec![-1, 2])
                .unwrap(),
        );
        let before = bond.clone();
        assert!(bond.set_computed_prop("", vec![3]).is_err());
        assert_eq!(bond, before);
        let qb = QueryBond::from_carrier_parts(
            before.clone(),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
        );
        assert_eq!(qb.bond().prop("v"), before.prop("v"));
    }
}

#[cfg(test)]
mod uint_dependency_proposed_tests {
    use super::*;
    #[test]
    fn proposed_uint_identity_width_and_strict_accessors() {
        for value in [
            0_u32,
            1,
            i32::MAX as u32 - 1,
            i32::MAX as u32,
            i32::MAX as u32 + 1,
            u32::MAX,
        ] {
            let property = PropertyValue::from(value);
            assert_eq!(property.kind(), PropertyValueKind::UInt);
            assert_eq!(property.as_uint(), Ok(value));
            assert_eq!(property.clone(), property);
            assert_eq!(
                property.as_int().unwrap_err().actual(),
                PropertyValueKind::UInt
            );
            assert!(property.as_string().is_err());
            assert!(property.as_bool().is_err());
            assert!(property.as_double().is_err());
            assert!(property.as_int_vector().is_err());
        }
        assert_ne!(PropertyValue::UInt(1), PropertyValue::Int(1));
        assert!(PropertyValue::Int(1).as_uint().is_err());
    }
}

#[cfg(test)]
mod uint_source_transport_proposed_tests {
    use super::*;
    use crate::{
        Atom, AtomId, AtomQueryPredicate, AtomSpec, Bond, BondId, BondOrder, BondQueryPredicate,
        BondSpec, Element, QueryAtom, QueryBond, QueryNode,
    };
    #[test]
    fn proposed_uint_computed_order_clone_and_query_transport() {
        for number in [0_u32, 1, 2147483646, 2147483647, 2147483648, 4294967295] {
            let mut store = PropertyStore::new();
            store.set("a".into(), PropertyValue::Int(1));
            store.set_computed("rank".into(), PropertyValue::UInt(number));
            store.set("rank".into(), PropertyValue::UInt(number));
            assert!(store.is_computed("rank"));
            assert_eq!(store.ordered_keys(), ["a", "rank"]);
            let saved = store.clone();
            store.clear_computed();
            assert_eq!(saved.get("rank"), Some(&PropertyValue::UInt(number)));
            assert_eq!(store.ordered_keys(), ["a"]);
            let mut atom = Atom::from_spec(
                AtomId::new(0),
                AtomSpec::new(Element::C)
                    .with_computed_prop("rank", PropertyValue::UInt(number))
                    .unwrap(),
            );
            let original = atom.clone();
            assert!(atom.set_prop("", PropertyValue::UInt(number)).is_err());
            assert_eq!(atom, original);
            let mut query = QueryAtom::from_carrier_parts(
                original.clone(),
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            );
            assert_eq!(query.prop("rank"), Some(&PropertyValue::UInt(number)));
            assert!(query.is_prop_computed("rank"));
            query.clear_computed_props();
            assert_eq!(query.prop("rank"), None);
            assert_eq!(original.prop("rank"), Some(&PropertyValue::UInt(number)));
            let bond = Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                    .with_computed_prop("rank", PropertyValue::UInt(number))
                    .unwrap(),
            );
            let query = QueryBond::from_carrier_parts(
                bond.clone(),
                QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
            );
            assert_eq!(query.bond().prop("rank"), bond.prop("rank"));
            assert!(query.bond().is_prop_computed("rank"));
        }
    }
}

#[cfg(test)]
mod uint_complete_source_condition_cells {
    use super::*;
    use crate::{
        Atom, AtomId, AtomQueryPredicate, AtomSpec, Bond, BondId, BondQueryPredicate, BondSpec,
        QueryAtom, QueryBond, QueryNode, SdfPropertyList, SdfPropertyListTarget, TopologyBlock,
    };
    use cosmolkit_types::{BondOrder, Element};

    // FROZEN UINT CONDITION: STRICT_AND_WIDTH_0
    #[test]
    fn uint_cell_strict_and_width_0() {
        let v = PropertyValue::UInt(0_u32);
        assert_eq!(v.kind(), PropertyValueKind::UInt);
        assert_eq!(v.as_uint(), Ok(0_u32));
        assert_eq!(v.clone(), PropertyValue::UInt(0_u32));
        assert_ne!(v, PropertyValue::Int(0));
        let e = v.as_int().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::Int);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        let e = v.as_string().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::String);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        let e = v.as_double().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::Double);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        let e = v.as_bool().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::Bool);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        let e = v.as_int_vector().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::IntVector);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        assert_eq!(
            PropertyValue::Int(1).as_uint().unwrap_err().expected(),
            PropertyValueKind::UInt
        );
    }
    // FROZEN UINT CONDITION: STRICT_AND_WIDTH_1
    #[test]
    fn uint_cell_strict_and_width_1() {
        let v = PropertyValue::UInt(1_u32);
        assert_eq!(v.kind(), PropertyValueKind::UInt);
        assert_eq!(v.as_uint(), Ok(1_u32));
        assert_eq!(v.clone(), PropertyValue::UInt(1_u32));
        assert_ne!(v, PropertyValue::Int(1));
        let e = v.as_int().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::Int);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        let e = v.as_string().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::String);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        let e = v.as_double().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::Double);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        let e = v.as_bool().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::Bool);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        let e = v.as_int_vector().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::IntVector);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        assert_eq!(
            PropertyValue::Int(1).as_uint().unwrap_err().expected(),
            PropertyValueKind::UInt
        );
    }
    // FROZEN UINT CONDITION: STRICT_AND_WIDTH_2147483646
    #[test]
    fn uint_cell_strict_and_width_2147483646() {
        let v = PropertyValue::UInt(2147483646_u32);
        assert_eq!(v.kind(), PropertyValueKind::UInt);
        assert_eq!(v.as_uint(), Ok(2147483646_u32));
        assert_eq!(v.clone(), PropertyValue::UInt(2147483646_u32));
        assert_ne!(v, PropertyValue::Int(2147483646));
        let e = v.as_int().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::Int);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        let e = v.as_string().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::String);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        let e = v.as_double().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::Double);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        let e = v.as_bool().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::Bool);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        let e = v.as_int_vector().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::IntVector);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        assert_eq!(
            PropertyValue::Int(1).as_uint().unwrap_err().expected(),
            PropertyValueKind::UInt
        );
    }
    // FROZEN UINT CONDITION: STRICT_AND_WIDTH_2147483647
    #[test]
    fn uint_cell_strict_and_width_2147483647() {
        let v = PropertyValue::UInt(2147483647_u32);
        assert_eq!(v.kind(), PropertyValueKind::UInt);
        assert_eq!(v.as_uint(), Ok(2147483647_u32));
        assert_eq!(v.clone(), PropertyValue::UInt(2147483647_u32));
        assert_ne!(v, PropertyValue::Int(2147483647));
        let e = v.as_int().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::Int);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        let e = v.as_string().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::String);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        let e = v.as_double().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::Double);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        let e = v.as_bool().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::Bool);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        let e = v.as_int_vector().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::IntVector);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        assert_eq!(
            PropertyValue::Int(1).as_uint().unwrap_err().expected(),
            PropertyValueKind::UInt
        );
    }
    // FROZEN UINT CONDITION: STRICT_AND_WIDTH_2147483648
    #[test]
    fn uint_cell_strict_and_width_2147483648() {
        let v = PropertyValue::UInt(2147483648_u32);
        assert_eq!(v.kind(), PropertyValueKind::UInt);
        assert_eq!(v.as_uint(), Ok(2147483648_u32));
        assert_eq!(v.clone(), PropertyValue::UInt(2147483648_u32));
        assert_ne!(v, PropertyValue::Int(2147483647));
        let e = v.as_int().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::Int);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        let e = v.as_string().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::String);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        let e = v.as_double().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::Double);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        let e = v.as_bool().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::Bool);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        let e = v.as_int_vector().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::IntVector);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        assert_eq!(
            PropertyValue::Int(1).as_uint().unwrap_err().expected(),
            PropertyValueKind::UInt
        );
    }
    // FROZEN UINT CONDITION: STRICT_AND_WIDTH_4294967295
    #[test]
    fn uint_cell_strict_and_width_4294967295() {
        let v = PropertyValue::UInt(4294967295_u32);
        assert_eq!(v.kind(), PropertyValueKind::UInt);
        assert_eq!(v.as_uint(), Ok(4294967295_u32));
        assert_eq!(v.clone(), PropertyValue::UInt(4294967295_u32));
        assert_ne!(v, PropertyValue::Int(2147483647));
        let e = v.as_int().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::Int);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        let e = v.as_string().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::String);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        let e = v.as_double().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::Double);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        let e = v.as_bool().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::Bool);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        let e = v.as_int_vector().unwrap_err();
        assert_eq!(e.expected(), PropertyValueKind::IntVector);
        assert_eq!(e.actual(), PropertyValueKind::UInt);
        assert_eq!(
            PropertyValue::Int(1).as_uint().unwrap_err().expected(),
            PropertyValueKind::UInt
        );
    }
    // FROZEN UINT CONDITION: CLONE_COMPUTED_0
    #[test]
    fn uint_cell_clone_computed_0_property_value() {
        let mut store = PropertyStore::new();
        store.set("a".into(), PropertyValue::Int(1));
        store.set_computed("rank".into(), PropertyValue::UInt(0_u32));
        store.set("rank".into(), PropertyValue::UInt(0_u32));
        assert_eq!(store.ordered_keys(), ["a", "rank"]);
        assert!(store.is_computed("rank"));
        let saved = store.clone();
        store.clear_computed();
        assert_eq!(saved.get("rank"), Some(&PropertyValue::UInt(0_u32)));
        assert_eq!(store.ordered_keys(), ["a"]);
        let atom = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_computed_prop("rank", PropertyValue::UInt(0_u32))
                .unwrap(),
        );
        let mut query = QueryAtom::from_carrier_parts(
            atom.clone(),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
        );
        assert_eq!(query.prop("rank"), Some(&PropertyValue::UInt(0_u32)));
        assert!(query.is_prop_computed("rank"));
        query.clear_computed_props();
        assert_eq!(query.prop("rank"), None);
        assert_eq!(atom.prop("rank"), Some(&PropertyValue::UInt(0_u32)));
        let bond = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single)
                .with_computed_prop("rank", PropertyValue::UInt(0_u32))
                .unwrap(),
        );
        let querybond = QueryBond::from_carrier_parts(
            bond.clone(),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
        );
        assert_eq!(
            querybond.bond().prop("rank"),
            Some(&PropertyValue::UInt(0_u32))
        );
        assert!(querybond.bond().is_prop_computed("rank"));
        let g = TopologyBlock::try_from_parts(
            vec![
                atom.clone(),
                Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
                Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::C)),
            ],
            vec![bond],
            vec![],
            vec![],
        )
        .unwrap();
        let before = g.clone();
        let (reordered, _) = g
            .reordered_atoms(&[AtomId::new(2), AtomId::new(1), AtomId::new(0)])
            .unwrap();
        assert_eq!(
            reordered.atoms[2].prop("rank"),
            Some(&PropertyValue::UInt(0_u32))
        );
        assert_eq!(
            reordered.bonds[0].prop("rank"),
            Some(&PropertyValue::UInt(0_u32))
        );
        assert!(reordered.atoms[2].is_prop_computed("rank"));
        assert!(reordered.bonds[0].is_prop_computed("rank"));
        let mut edit = g.begin_batch_edit().unwrap();
        edit.remove_atom(AtomId::new(1)).unwrap();
        let (fragment, _) = edit.finish().unwrap();
        assert_eq!(
            fragment.atoms[0].prop("rank"),
            Some(&PropertyValue::UInt(0_u32))
        );
        assert_eq!(
            fragment.bonds[0].prop("rank"),
            Some(&PropertyValue::UInt(0_u32))
        );
        assert_eq!(g, before);
        let mut props = crate::MoleculeProperties::default()
            .with_sdf_property_list(SdfPropertyList::new(
                SdfPropertyListTarget::Atom,
                "rank",
                vec![
                    Some(PropertyValue::UInt(0_u32)),
                    None,
                    Some(PropertyValue::UInt(0_u32)),
                ],
            ))
            .with_sdf_property_list(SdfPropertyList::new(
                SdfPropertyListTarget::Bond,
                "rank",
                vec![Some(PropertyValue::UInt(0_u32))],
            ));
        let original = props.clone();
        props.remap_topology(
            &[Some(AtomId::new(2)), Some(AtomId::new(0))],
            &[Some(BondId::new(0))],
        );
        assert_eq!(
            props.sdf_property_lists()[0].values(),
            &[
                Some(PropertyValue::UInt(0_u32)),
                Some(PropertyValue::UInt(0_u32))
            ]
        );
        assert_eq!(
            props.sdf_property_lists()[1].values(),
            &[Some(PropertyValue::UInt(0_u32))]
        );
        assert_eq!(
            original.sdf_property_lists()[0].values(),
            &[
                Some(PropertyValue::UInt(0_u32)),
                None,
                Some(PropertyValue::UInt(0_u32))
            ]
        );
    }
    // FROZEN UINT CONDITION: CLONE_COMPUTED_1
    #[test]
    fn uint_cell_clone_computed_1_property_value() {
        let mut store = PropertyStore::new();
        store.set("a".into(), PropertyValue::Int(1));
        store.set_computed("rank".into(), PropertyValue::UInt(1_u32));
        store.set("rank".into(), PropertyValue::UInt(1_u32));
        assert_eq!(store.ordered_keys(), ["a", "rank"]);
        assert!(store.is_computed("rank"));
        let saved = store.clone();
        store.clear_computed();
        assert_eq!(saved.get("rank"), Some(&PropertyValue::UInt(1_u32)));
        assert_eq!(store.ordered_keys(), ["a"]);
        let atom = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_computed_prop("rank", PropertyValue::UInt(1_u32))
                .unwrap(),
        );
        let mut query = QueryAtom::from_carrier_parts(
            atom.clone(),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
        );
        assert_eq!(query.prop("rank"), Some(&PropertyValue::UInt(1_u32)));
        assert!(query.is_prop_computed("rank"));
        query.clear_computed_props();
        assert_eq!(query.prop("rank"), None);
        assert_eq!(atom.prop("rank"), Some(&PropertyValue::UInt(1_u32)));
        let bond = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single)
                .with_computed_prop("rank", PropertyValue::UInt(1_u32))
                .unwrap(),
        );
        let querybond = QueryBond::from_carrier_parts(
            bond.clone(),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
        );
        assert_eq!(
            querybond.bond().prop("rank"),
            Some(&PropertyValue::UInt(1_u32))
        );
        assert!(querybond.bond().is_prop_computed("rank"));
        let g = TopologyBlock::try_from_parts(
            vec![
                atom.clone(),
                Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
                Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::C)),
            ],
            vec![bond],
            vec![],
            vec![],
        )
        .unwrap();
        let before = g.clone();
        let (reordered, _) = g
            .reordered_atoms(&[AtomId::new(2), AtomId::new(1), AtomId::new(0)])
            .unwrap();
        assert_eq!(
            reordered.atoms[2].prop("rank"),
            Some(&PropertyValue::UInt(1_u32))
        );
        assert_eq!(
            reordered.bonds[0].prop("rank"),
            Some(&PropertyValue::UInt(1_u32))
        );
        assert!(reordered.atoms[2].is_prop_computed("rank"));
        assert!(reordered.bonds[0].is_prop_computed("rank"));
        let mut edit = g.begin_batch_edit().unwrap();
        edit.remove_atom(AtomId::new(1)).unwrap();
        let (fragment, _) = edit.finish().unwrap();
        assert_eq!(
            fragment.atoms[0].prop("rank"),
            Some(&PropertyValue::UInt(1_u32))
        );
        assert_eq!(
            fragment.bonds[0].prop("rank"),
            Some(&PropertyValue::UInt(1_u32))
        );
        assert_eq!(g, before);
        let mut props = crate::MoleculeProperties::default()
            .with_sdf_property_list(SdfPropertyList::new(
                SdfPropertyListTarget::Atom,
                "rank",
                vec![
                    Some(PropertyValue::UInt(1_u32)),
                    None,
                    Some(PropertyValue::UInt(1_u32)),
                ],
            ))
            .with_sdf_property_list(SdfPropertyList::new(
                SdfPropertyListTarget::Bond,
                "rank",
                vec![Some(PropertyValue::UInt(1_u32))],
            ));
        let original = props.clone();
        props.remap_topology(
            &[Some(AtomId::new(2)), Some(AtomId::new(0))],
            &[Some(BondId::new(0))],
        );
        assert_eq!(
            props.sdf_property_lists()[0].values(),
            &[
                Some(PropertyValue::UInt(1_u32)),
                Some(PropertyValue::UInt(1_u32))
            ]
        );
        assert_eq!(
            props.sdf_property_lists()[1].values(),
            &[Some(PropertyValue::UInt(1_u32))]
        );
        assert_eq!(
            original.sdf_property_lists()[0].values(),
            &[
                Some(PropertyValue::UInt(1_u32)),
                None,
                Some(PropertyValue::UInt(1_u32))
            ]
        );
    }
    // FROZEN UINT CONDITION: CLONE_COMPUTED_2147483646
    #[test]
    fn uint_cell_clone_computed_2147483646_property_value() {
        let mut store = PropertyStore::new();
        store.set("a".into(), PropertyValue::Int(1));
        store.set_computed("rank".into(), PropertyValue::UInt(2147483646_u32));
        store.set("rank".into(), PropertyValue::UInt(2147483646_u32));
        assert_eq!(store.ordered_keys(), ["a", "rank"]);
        assert!(store.is_computed("rank"));
        let saved = store.clone();
        store.clear_computed();
        assert_eq!(
            saved.get("rank"),
            Some(&PropertyValue::UInt(2147483646_u32))
        );
        assert_eq!(store.ordered_keys(), ["a"]);
        let atom = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_computed_prop("rank", PropertyValue::UInt(2147483646_u32))
                .unwrap(),
        );
        let mut query = QueryAtom::from_carrier_parts(
            atom.clone(),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
        );
        assert_eq!(
            query.prop("rank"),
            Some(&PropertyValue::UInt(2147483646_u32))
        );
        assert!(query.is_prop_computed("rank"));
        query.clear_computed_props();
        assert_eq!(query.prop("rank"), None);
        assert_eq!(
            atom.prop("rank"),
            Some(&PropertyValue::UInt(2147483646_u32))
        );
        let bond = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single)
                .with_computed_prop("rank", PropertyValue::UInt(2147483646_u32))
                .unwrap(),
        );
        let querybond = QueryBond::from_carrier_parts(
            bond.clone(),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
        );
        assert_eq!(
            querybond.bond().prop("rank"),
            Some(&PropertyValue::UInt(2147483646_u32))
        );
        assert!(querybond.bond().is_prop_computed("rank"));
        let g = TopologyBlock::try_from_parts(
            vec![
                atom.clone(),
                Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
                Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::C)),
            ],
            vec![bond],
            vec![],
            vec![],
        )
        .unwrap();
        let before = g.clone();
        let (reordered, _) = g
            .reordered_atoms(&[AtomId::new(2), AtomId::new(1), AtomId::new(0)])
            .unwrap();
        assert_eq!(
            reordered.atoms[2].prop("rank"),
            Some(&PropertyValue::UInt(2147483646_u32))
        );
        assert_eq!(
            reordered.bonds[0].prop("rank"),
            Some(&PropertyValue::UInt(2147483646_u32))
        );
        assert!(reordered.atoms[2].is_prop_computed("rank"));
        assert!(reordered.bonds[0].is_prop_computed("rank"));
        let mut edit = g.begin_batch_edit().unwrap();
        edit.remove_atom(AtomId::new(1)).unwrap();
        let (fragment, _) = edit.finish().unwrap();
        assert_eq!(
            fragment.atoms[0].prop("rank"),
            Some(&PropertyValue::UInt(2147483646_u32))
        );
        assert_eq!(
            fragment.bonds[0].prop("rank"),
            Some(&PropertyValue::UInt(2147483646_u32))
        );
        assert_eq!(g, before);
        let mut props = crate::MoleculeProperties::default()
            .with_sdf_property_list(SdfPropertyList::new(
                SdfPropertyListTarget::Atom,
                "rank",
                vec![
                    Some(PropertyValue::UInt(2147483646_u32)),
                    None,
                    Some(PropertyValue::UInt(2147483646_u32)),
                ],
            ))
            .with_sdf_property_list(SdfPropertyList::new(
                SdfPropertyListTarget::Bond,
                "rank",
                vec![Some(PropertyValue::UInt(2147483646_u32))],
            ));
        let original = props.clone();
        props.remap_topology(
            &[Some(AtomId::new(2)), Some(AtomId::new(0))],
            &[Some(BondId::new(0))],
        );
        assert_eq!(
            props.sdf_property_lists()[0].values(),
            &[
                Some(PropertyValue::UInt(2147483646_u32)),
                Some(PropertyValue::UInt(2147483646_u32))
            ]
        );
        assert_eq!(
            props.sdf_property_lists()[1].values(),
            &[Some(PropertyValue::UInt(2147483646_u32))]
        );
        assert_eq!(
            original.sdf_property_lists()[0].values(),
            &[
                Some(PropertyValue::UInt(2147483646_u32)),
                None,
                Some(PropertyValue::UInt(2147483646_u32))
            ]
        );
    }
    // FROZEN UINT CONDITION: CLONE_COMPUTED_2147483647
    #[test]
    fn uint_cell_clone_computed_2147483647_property_value() {
        let mut store = PropertyStore::new();
        store.set("a".into(), PropertyValue::Int(1));
        store.set_computed("rank".into(), PropertyValue::UInt(2147483647_u32));
        store.set("rank".into(), PropertyValue::UInt(2147483647_u32));
        assert_eq!(store.ordered_keys(), ["a", "rank"]);
        assert!(store.is_computed("rank"));
        let saved = store.clone();
        store.clear_computed();
        assert_eq!(
            saved.get("rank"),
            Some(&PropertyValue::UInt(2147483647_u32))
        );
        assert_eq!(store.ordered_keys(), ["a"]);
        let atom = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_computed_prop("rank", PropertyValue::UInt(2147483647_u32))
                .unwrap(),
        );
        let mut query = QueryAtom::from_carrier_parts(
            atom.clone(),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
        );
        assert_eq!(
            query.prop("rank"),
            Some(&PropertyValue::UInt(2147483647_u32))
        );
        assert!(query.is_prop_computed("rank"));
        query.clear_computed_props();
        assert_eq!(query.prop("rank"), None);
        assert_eq!(
            atom.prop("rank"),
            Some(&PropertyValue::UInt(2147483647_u32))
        );
        let bond = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single)
                .with_computed_prop("rank", PropertyValue::UInt(2147483647_u32))
                .unwrap(),
        );
        let querybond = QueryBond::from_carrier_parts(
            bond.clone(),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
        );
        assert_eq!(
            querybond.bond().prop("rank"),
            Some(&PropertyValue::UInt(2147483647_u32))
        );
        assert!(querybond.bond().is_prop_computed("rank"));
        let g = TopologyBlock::try_from_parts(
            vec![
                atom.clone(),
                Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
                Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::C)),
            ],
            vec![bond],
            vec![],
            vec![],
        )
        .unwrap();
        let before = g.clone();
        let (reordered, _) = g
            .reordered_atoms(&[AtomId::new(2), AtomId::new(1), AtomId::new(0)])
            .unwrap();
        assert_eq!(
            reordered.atoms[2].prop("rank"),
            Some(&PropertyValue::UInt(2147483647_u32))
        );
        assert_eq!(
            reordered.bonds[0].prop("rank"),
            Some(&PropertyValue::UInt(2147483647_u32))
        );
        assert!(reordered.atoms[2].is_prop_computed("rank"));
        assert!(reordered.bonds[0].is_prop_computed("rank"));
        let mut edit = g.begin_batch_edit().unwrap();
        edit.remove_atom(AtomId::new(1)).unwrap();
        let (fragment, _) = edit.finish().unwrap();
        assert_eq!(
            fragment.atoms[0].prop("rank"),
            Some(&PropertyValue::UInt(2147483647_u32))
        );
        assert_eq!(
            fragment.bonds[0].prop("rank"),
            Some(&PropertyValue::UInt(2147483647_u32))
        );
        assert_eq!(g, before);
        let mut props = crate::MoleculeProperties::default()
            .with_sdf_property_list(SdfPropertyList::new(
                SdfPropertyListTarget::Atom,
                "rank",
                vec![
                    Some(PropertyValue::UInt(2147483647_u32)),
                    None,
                    Some(PropertyValue::UInt(2147483647_u32)),
                ],
            ))
            .with_sdf_property_list(SdfPropertyList::new(
                SdfPropertyListTarget::Bond,
                "rank",
                vec![Some(PropertyValue::UInt(2147483647_u32))],
            ));
        let original = props.clone();
        props.remap_topology(
            &[Some(AtomId::new(2)), Some(AtomId::new(0))],
            &[Some(BondId::new(0))],
        );
        assert_eq!(
            props.sdf_property_lists()[0].values(),
            &[
                Some(PropertyValue::UInt(2147483647_u32)),
                Some(PropertyValue::UInt(2147483647_u32))
            ]
        );
        assert_eq!(
            props.sdf_property_lists()[1].values(),
            &[Some(PropertyValue::UInt(2147483647_u32))]
        );
        assert_eq!(
            original.sdf_property_lists()[0].values(),
            &[
                Some(PropertyValue::UInt(2147483647_u32)),
                None,
                Some(PropertyValue::UInt(2147483647_u32))
            ]
        );
    }
    // FROZEN UINT CONDITION: CLONE_COMPUTED_2147483648
    #[test]
    fn uint_cell_clone_computed_2147483648_property_value() {
        let mut store = PropertyStore::new();
        store.set("a".into(), PropertyValue::Int(1));
        store.set_computed("rank".into(), PropertyValue::UInt(2147483648_u32));
        store.set("rank".into(), PropertyValue::UInt(2147483648_u32));
        assert_eq!(store.ordered_keys(), ["a", "rank"]);
        assert!(store.is_computed("rank"));
        let saved = store.clone();
        store.clear_computed();
        assert_eq!(
            saved.get("rank"),
            Some(&PropertyValue::UInt(2147483648_u32))
        );
        assert_eq!(store.ordered_keys(), ["a"]);
        let atom = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_computed_prop("rank", PropertyValue::UInt(2147483648_u32))
                .unwrap(),
        );
        let mut query = QueryAtom::from_carrier_parts(
            atom.clone(),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
        );
        assert_eq!(
            query.prop("rank"),
            Some(&PropertyValue::UInt(2147483648_u32))
        );
        assert!(query.is_prop_computed("rank"));
        query.clear_computed_props();
        assert_eq!(query.prop("rank"), None);
        assert_eq!(
            atom.prop("rank"),
            Some(&PropertyValue::UInt(2147483648_u32))
        );
        let bond = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single)
                .with_computed_prop("rank", PropertyValue::UInt(2147483648_u32))
                .unwrap(),
        );
        let querybond = QueryBond::from_carrier_parts(
            bond.clone(),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
        );
        assert_eq!(
            querybond.bond().prop("rank"),
            Some(&PropertyValue::UInt(2147483648_u32))
        );
        assert!(querybond.bond().is_prop_computed("rank"));
        let g = TopologyBlock::try_from_parts(
            vec![
                atom.clone(),
                Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
                Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::C)),
            ],
            vec![bond],
            vec![],
            vec![],
        )
        .unwrap();
        let before = g.clone();
        let (reordered, _) = g
            .reordered_atoms(&[AtomId::new(2), AtomId::new(1), AtomId::new(0)])
            .unwrap();
        assert_eq!(
            reordered.atoms[2].prop("rank"),
            Some(&PropertyValue::UInt(2147483648_u32))
        );
        assert_eq!(
            reordered.bonds[0].prop("rank"),
            Some(&PropertyValue::UInt(2147483648_u32))
        );
        assert!(reordered.atoms[2].is_prop_computed("rank"));
        assert!(reordered.bonds[0].is_prop_computed("rank"));
        let mut edit = g.begin_batch_edit().unwrap();
        edit.remove_atom(AtomId::new(1)).unwrap();
        let (fragment, _) = edit.finish().unwrap();
        assert_eq!(
            fragment.atoms[0].prop("rank"),
            Some(&PropertyValue::UInt(2147483648_u32))
        );
        assert_eq!(
            fragment.bonds[0].prop("rank"),
            Some(&PropertyValue::UInt(2147483648_u32))
        );
        assert_eq!(g, before);
        let mut props = crate::MoleculeProperties::default()
            .with_sdf_property_list(SdfPropertyList::new(
                SdfPropertyListTarget::Atom,
                "rank",
                vec![
                    Some(PropertyValue::UInt(2147483648_u32)),
                    None,
                    Some(PropertyValue::UInt(2147483648_u32)),
                ],
            ))
            .with_sdf_property_list(SdfPropertyList::new(
                SdfPropertyListTarget::Bond,
                "rank",
                vec![Some(PropertyValue::UInt(2147483648_u32))],
            ));
        let original = props.clone();
        props.remap_topology(
            &[Some(AtomId::new(2)), Some(AtomId::new(0))],
            &[Some(BondId::new(0))],
        );
        assert_eq!(
            props.sdf_property_lists()[0].values(),
            &[
                Some(PropertyValue::UInt(2147483648_u32)),
                Some(PropertyValue::UInt(2147483648_u32))
            ]
        );
        assert_eq!(
            props.sdf_property_lists()[1].values(),
            &[Some(PropertyValue::UInt(2147483648_u32))]
        );
        assert_eq!(
            original.sdf_property_lists()[0].values(),
            &[
                Some(PropertyValue::UInt(2147483648_u32)),
                None,
                Some(PropertyValue::UInt(2147483648_u32))
            ]
        );
    }
    // FROZEN UINT CONDITION: CLONE_COMPUTED_4294967295
    #[test]
    fn uint_cell_clone_computed_4294967295_property_value() {
        let mut store = PropertyStore::new();
        store.set("a".into(), PropertyValue::Int(1));
        store.set_computed("rank".into(), PropertyValue::UInt(4294967295_u32));
        store.set("rank".into(), PropertyValue::UInt(4294967295_u32));
        assert_eq!(store.ordered_keys(), ["a", "rank"]);
        assert!(store.is_computed("rank"));
        let saved = store.clone();
        store.clear_computed();
        assert_eq!(
            saved.get("rank"),
            Some(&PropertyValue::UInt(4294967295_u32))
        );
        assert_eq!(store.ordered_keys(), ["a"]);
        let atom = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_computed_prop("rank", PropertyValue::UInt(4294967295_u32))
                .unwrap(),
        );
        let mut query = QueryAtom::from_carrier_parts(
            atom.clone(),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
        );
        assert_eq!(
            query.prop("rank"),
            Some(&PropertyValue::UInt(4294967295_u32))
        );
        assert!(query.is_prop_computed("rank"));
        query.clear_computed_props();
        assert_eq!(query.prop("rank"), None);
        assert_eq!(
            atom.prop("rank"),
            Some(&PropertyValue::UInt(4294967295_u32))
        );
        let bond = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single)
                .with_computed_prop("rank", PropertyValue::UInt(4294967295_u32))
                .unwrap(),
        );
        let querybond = QueryBond::from_carrier_parts(
            bond.clone(),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
        );
        assert_eq!(
            querybond.bond().prop("rank"),
            Some(&PropertyValue::UInt(4294967295_u32))
        );
        assert!(querybond.bond().is_prop_computed("rank"));
        let g = TopologyBlock::try_from_parts(
            vec![
                atom.clone(),
                Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
                Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::C)),
            ],
            vec![bond],
            vec![],
            vec![],
        )
        .unwrap();
        let before = g.clone();
        let (reordered, _) = g
            .reordered_atoms(&[AtomId::new(2), AtomId::new(1), AtomId::new(0)])
            .unwrap();
        assert_eq!(
            reordered.atoms[2].prop("rank"),
            Some(&PropertyValue::UInt(4294967295_u32))
        );
        assert_eq!(
            reordered.bonds[0].prop("rank"),
            Some(&PropertyValue::UInt(4294967295_u32))
        );
        assert!(reordered.atoms[2].is_prop_computed("rank"));
        assert!(reordered.bonds[0].is_prop_computed("rank"));
        let mut edit = g.begin_batch_edit().unwrap();
        edit.remove_atom(AtomId::new(1)).unwrap();
        let (fragment, _) = edit.finish().unwrap();
        assert_eq!(
            fragment.atoms[0].prop("rank"),
            Some(&PropertyValue::UInt(4294967295_u32))
        );
        assert_eq!(
            fragment.bonds[0].prop("rank"),
            Some(&PropertyValue::UInt(4294967295_u32))
        );
        assert_eq!(g, before);
        let mut props = crate::MoleculeProperties::default()
            .with_sdf_property_list(SdfPropertyList::new(
                SdfPropertyListTarget::Atom,
                "rank",
                vec![
                    Some(PropertyValue::UInt(4294967295_u32)),
                    None,
                    Some(PropertyValue::UInt(4294967295_u32)),
                ],
            ))
            .with_sdf_property_list(SdfPropertyList::new(
                SdfPropertyListTarget::Bond,
                "rank",
                vec![Some(PropertyValue::UInt(4294967295_u32))],
            ));
        let original = props.clone();
        props.remap_topology(
            &[Some(AtomId::new(2)), Some(AtomId::new(0))],
            &[Some(BondId::new(0))],
        );
        assert_eq!(
            props.sdf_property_lists()[0].values(),
            &[
                Some(PropertyValue::UInt(4294967295_u32)),
                Some(PropertyValue::UInt(4294967295_u32))
            ]
        );
        assert_eq!(
            props.sdf_property_lists()[1].values(),
            &[Some(PropertyValue::UInt(4294967295_u32))]
        );
        assert_eq!(
            original.sdf_property_lists()[0].values(),
            &[
                Some(PropertyValue::UInt(4294967295_u32)),
                None,
                Some(PropertyValue::UInt(4294967295_u32))
            ]
        );
    }
    // FROZEN UINT CONDITION: STRICT_BOOL_VECTOR
    #[test]
    fn uint_cell_strict_bool_vector_property_value() {
        for n in [0_u32, 1, 4294967295] {
            let v = PropertyValue::UInt(n);
            let e = v.as_bool().unwrap_err();
            assert_eq!(e.expected(), PropertyValueKind::Bool);
            assert_eq!(e.actual(), PropertyValueKind::UInt);
            let e = v.as_int_vector().unwrap_err();
            assert_eq!(e.expected(), PropertyValueKind::IntVector);
            assert_eq!(e.actual(), PropertyValueKind::UInt);
        }
    }
}
