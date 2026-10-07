//! Canonical detached atom and bond property values.

use crate::PropertyText;
use std::collections::BTreeMap;

/// The modeled source value kinds supported by atom and bond properties.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PropertyValueKind {
    String,
    Int,
    UInt,
    IntVector,
    StringVector,
    Double,
    Bool,
}

/// A canonical detached atom or bond property value.
#[derive(Debug)]
pub enum PropertyValue {
    String(PropertyText),
    Int(i32),
    UInt(u32),
    IntVector(Vec<i32>),
    StringVector(Vec<PropertyText>),
    Double(f64),
    Bool(bool),
}

impl Clone for PropertyValue {
    fn clone(&self) -> Self {
        // RDKit✔️✔️: inline void copy_rdvalue(RDValue &dest, const RDValue &src) {
        // RDKit✔️✔️:   if (&dest == &src) {  // don't copy over yourself
        // RDKit✔️✔️:     return;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   dest.destroy();
        // RDKit✔️✔️:   dest.type = src.type;
        // RDKit✔️✔️:   switch (src.type) {
        // RDKit✔️✔️:     case RDTypeTag::StringTag:
        // RDKit✔️✔️:       dest.value.s = new std::string(*src.value.s);
        // RDKit✔️✔️:       break;
        // RDKit❌❌:     case RDTypeTag::AnyTag:
        // RDKit❌❌:       dest.value.a = new std::any(*src.value.a);
        // RDKit❌❌:       break;
        // RDKit❌❌:     case RDTypeTag::VecDoubleTag:
        // RDKit❌❌:       dest.value.vd = new std::vector<double>(*src.value.vd);
        // RDKit❌❌:       break;
        // RDKit❌❌:     case RDTypeTag::VecFloatTag:
        // RDKit❌❌:       dest.value.vf = new std::vector<float>(*src.value.vf);
        // RDKit❌❌:       break;
        // RDKit✔️✔️:     case RDTypeTag::VecIntTag:
        // RDKit✔️✔️:       dest.value.vi = new std::vector<int>(*src.value.vi);
        // RDKit✔️✔️:       break;
        // RDKit❌❌:     case RDTypeTag::VecUnsignedIntTag:
        // RDKit❌❌:       dest.value.vu = new std::vector<unsigned int>(*src.value.vu);
        // RDKit❌❌:       break;
        // RDKit✔️✔️:     case RDTypeTag::VecStringTag:
        // RDKit✔️✔️:       dest.value.vs = new std::vector<std::string>(*src.value.vs);
        // RDKit✔️✔️:       break;
        // RDKit✔️✔️:     default:
        // RDKit✔️✔️:       dest = src;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // RDKit✔️✔️:   void destroy() {
        // RDKit✔️✔️:     switch (type) {
        // RDKit✔️✔️:       case RDTypeTag::StringTag:
        // RDKit✔️✔️:         delete value.s;
        // RDKit✔️✔️:         break;
        // RDKit❌❌:       case RDTypeTag::AnyTag:
        // RDKit❌❌:         delete value.a;
        // RDKit❌❌:         break;
        // RDKit❌❌:       case RDTypeTag::VecDoubleTag:
        // RDKit❌❌:         delete value.vd;
        // RDKit❌❌:         break;
        // RDKit❌❌:       case RDTypeTag::VecFloatTag:
        // RDKit❌❌:         delete value.vf;
        // RDKit❌❌:         break;
        // RDKit✔️✔️:       case RDTypeTag::VecIntTag:
        // RDKit✔️✔️:         delete value.vi;
        // RDKit✔️✔️:         break;
        // RDKit❌❌:       case RDTypeTag::VecUnsignedIntTag:
        // RDKit❌❌:         delete value.vu;
        // RDKit❌❌:         break;
        // RDKit✔️✔️:       case RDTypeTag::VecStringTag:
        // RDKit✔️✔️:         delete value.vs;
        // RDKit✔️✔️:         break;
        // RDKit✔️✔️:       default:
        // RDKit✔️✔️:         break;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     type = RDTypeTag::EmptyTag;
        // RDKit✔️✔️:   }
        // Behavior: the new destination keeps the modeled tag and exact
        // payload. Source counted strings and vectors deep-copy their
        // bytes/elements, including each StringVector byte buffer. Scalars copy
        // bits. Rust ownership destroys a
        // replaced destination. Safe mutable and shared references cannot alias
        // that destination. Independent unmodeled source tags remain absent.
        // Complexity: O(elements + total bytes) vector/string-buffer deep copy, O(1)
        // for scalars; no decoding, validation pass or second property store.
        match self {
            Self::String(value) => Self::String(value.clone()),
            Self::Int(value) => Self::Int(*value),
            Self::UInt(value) => Self::UInt(*value),
            Self::IntVector(value) => Self::IntVector(value.clone()),
            Self::StringVector(value) => Self::StringVector(value.clone()),
            Self::Double(value) => Self::Double(*value),
            Self::Bool(value) => Self::Bool(*value),
        }
    }
}

impl PartialEq for PropertyValue {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::String(left), Self::String(right)) => left == right,
            (Self::Int(left), Self::Int(right)) => left == right,
            (Self::UInt(left), Self::UInt(right)) => left == right,
            (Self::IntVector(left), Self::IntVector(right)) => left == right,
            (Self::StringVector(left), Self::StringVector(right)) => left == right,
            (Self::Double(left), Self::Double(right)) => left.to_bits() == right.to_bits(),
            (Self::Bool(left), Self::Bool(right)) => left == right,
            _ => false,
        }
    }
}

impl Eq for PropertyValue {}

impl From<String> for PropertyValue {
    fn from(value: String) -> Self {
        Self::String(value.into())
    }
}

impl From<&str> for PropertyValue {
    fn from(value: &str) -> Self {
        Self::String(value.into())
    }
}

impl From<&String> for PropertyValue {
    fn from(value: &String) -> Self {
        Self::String(value.into())
    }
}

impl From<PropertyText> for PropertyValue {
    fn from(value: PropertyText) -> Self {
        Self::String(value)
    }
}

impl From<&PropertyText> for PropertyValue {
    fn from(value: &PropertyText) -> Self {
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

impl From<&[PropertyText]> for PropertyValue {
    fn from(value: &[PropertyText]) -> Self {
        // RDKit✔️✔️: inline RDValue(const std::vector<std::string> &v)
        // RDKit✔️✔️:       : value(new std::vector<std::string>(v)), type(RDTypeTag::VecStringTag) {}
        // Behavior: construct the exact VecStringTag payload, retaining vector
        // order, duplicates, empty strings and every counted string byte.
        // Each element owns a deep copy, independent of the borrowed input.
        // Complexity: O(elements + total bytes), one vector allocation plus
        // each nonempty string buffer, as in vector<string>'s copy constructor.
        // The enum owns all buffers; ordinary Rust destruction frees every
        // element and the vector once, without an auxiliary property store.
        Self::StringVector(value.to_vec())
    }
}

impl From<Vec<PropertyText>> for PropertyValue {
    fn from(value: Vec<PropertyText>) -> Self {
        // Owning input adapter transfers the same canonical vector and bytes.
        Self::StringVector(value)
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
            Self::StringVector(_) => PropertyValueKind::StringVector,
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

    pub fn as_string(&self) -> Result<&PropertyText, PropertyValueError> {
        // RDKit✔️✔️: inline std::string &rdvalue_cast<std::string &>(RDValue_cast_t v) {
        // RDKit✔️✔️:   if (rdvalue_is<std::string>(v)) {
        // RDKit✔️✔️:     return *v.ptrCast<std::string>();
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   throw std::bad_any_cast();
        // RDKit✔️✔️: }
        // RDKit✔️✔️: template <>
        // RDKit✔️✔️: inline short GetTag<std::string>() {
        // RDKit✔️✔️:   return StringTag;
        // RDKit✔️✔️: }
        // RDKit✔️✔️:   short getTag() const { return type; }
        // RDKit✔️✔️:   template <class T>
        // RDKit✔️✔️:   inline T *ptrCast() const {
        // RDKit✔️✔️:     return RDTypeTag::detail::valuePtrCast<T>(value);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: template <>
        // RDKit✔️✔️: inline std::string *valuePtrCast<std::string>(Value value) {
        // RDKit✔️✔️:   return value.s;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: template <class T>
        // RDKit✔️✔️: inline bool rdvalue_is(RDValue_cast_t v) {
        // RDKit✔️✔️:   const short tag =
        // RDKit✔️✔️:       RDTypeTag::GetTag<typename boost::remove_reference<T>::type>();
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // If we are an Any tag, check the any type info
        // RDKit✔️✔️:   //  see the template specialization below if we are
        // RDKit✔️✔️:   //  looking for a boost any directly
        // RDKit❌❌:   if (v.getTag() == RDTypeTag::AnyTag) {
        // RDKit❌❌:     return v.value.a->type() == typeid(T);
        // RDKit❌❌:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   if (v.getTag() == tag) {
        // RDKit✔️✔️:     return true;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   return false;
        // RDKit✔️✔️: }
        // Behavior: the enum's String discriminant is the source String tag;
        // borrow that counted payload unchanged, including arbitrary bytes.
        // Every other modeled tag yields the project exact-kind error. The
        // independent Any introspection capability is not modeled. This local
        // accessor grants a shared value borrow, no live storage mutation.
        // Complexity: constant discriminant dispatch and borrow, O(1), without
        // allocation, payload scan, scalar conversion or UTF8 validation.
        match self {
            Self::String(value) => Ok(value),
            _ => Err(self.wrong_kind(PropertyValueKind::String)),
        }
    }

    /// Copy the exact String-tag payload into a new owning text value.
    #[doc(hidden)]
    pub fn string_clone(&self) -> Result<PropertyText, PropertyValueError> {
        // RDKit✔️✔️: inline std::string rdvalue_cast<std::string>(RDValue_cast_t v) {
        // RDKit✔️✔️:   if (rdvalue_is<std::string>(v)) {
        // RDKit✔️✔️:     return *v.ptrCast<std::string>();
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   throw std::bad_any_cast();
        // RDKit✔️✔️: }
        // RDKit✔️✔️: template <>
        // RDKit✔️✔️: inline short GetTag<std::string>() {
        // RDKit✔️✔️:   return StringTag;
        // RDKit✔️✔️: }
        // RDKit✔️✔️:   short getTag() const { return type; }
        // RDKit✔️✔️:   template <class T>
        // RDKit✔️✔️:   inline T *ptrCast() const {
        // RDKit✔️✔️:     return RDTypeTag::detail::valuePtrCast<T>(value);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: template <>
        // RDKit✔️✔️: inline std::string *valuePtrCast<std::string>(Value value) {
        // RDKit✔️✔️:   return value.s;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: template <class T>
        // RDKit✔️✔️: inline bool rdvalue_is(RDValue_cast_t v) {
        // RDKit✔️✔️:   const short tag =
        // RDKit✔️✔️:       RDTypeTag::GetTag<typename boost::remove_reference<T>::type>();
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // If we are an Any tag, check the any type info
        // RDKit✔️✔️:   //  see the template specialization below if we are
        // RDKit✔️✔️:   //  looking for a boost any directly
        // RDKit❌❌:   if (v.getTag() == RDTypeTag::AnyTag) {
        // RDKit❌❌:     return v.value.a->type() == typeid(T);
        // RDKit❌❌:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   if (v.getTag() == tag) {
        // RDKit✔️✔️:     return true;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   return false;
        // RDKit✔️✔️: }
        // Behavior: reuse the canonical checked String borrow, then its one
        // counted deep copy. The source return-by-value copies std::string;
        // arbitrary bytes are preserved. Wrong modeled kinds retain their
        // structured expected/actual error, not conversion/default/Unsupported.
        // Complexity: O(payload) one byte-buffer copy on success, O(1) on
        // wrong-kind failure. No second getter algorithm, map or UTF8 decode.
        self.as_string().cloned()
    }

    /// Borrow the exact source string-vector payload without conversion.
    #[doc(hidden)]
    pub fn as_string_vector(&self) -> Result<&[PropertyText], PropertyValueError> {
        // RDKit✔️✔️: inline std::vector<std::string> &rdvalue_cast<std::vector<std::string> &>(
        // RDKit✔️✔️:     RDValue_cast_t v) {
        // RDKit✔️✔️:   if (rdvalue_is<std::vector<std::string>>(v)) {
        // RDKit✔️✔️:     return *v.ptrCast<std::vector<std::string>>();
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   throw std::bad_any_cast();
        // RDKit✔️✔️: }
        // RDKit✔️✔️: template <>
        // RDKit✔️✔️: inline short GetTag<std::vector<std::string>>() {
        // RDKit✔️✔️:   return VecStringTag;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: short getTag() const { return type; }
        // RDKit✔️✔️: template <class T>
        // RDKit✔️✔️:   inline T *ptrCast() const {
        // RDKit✔️✔️:     return RDTypeTag::detail::valuePtrCast<T>(value);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: template <>
        // RDKit✔️✔️: inline std::vector<std::string> *valuePtrCast<std::vector<std::string>>(
        // RDKit✔️✔️:     Value value) {
        // RDKit✔️✔️:   return value.vs;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: template <class T>
        // RDKit✔️✔️: inline bool rdvalue_is(RDValue_cast_t v) {
        // RDKit✔️✔️:   const short tag =
        // RDKit✔️✔️:       RDTypeTag::GetTag<typename boost::remove_reference<T>::type>();
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // If we are an Any tag, check the any type info
        // RDKit✔️✔️:   //  see the template specialization below if we are
        // RDKit✔️✔️:   //  looking for a boost any directly
        // RDKit❌❌:   if (v.getTag() == RDTypeTag::AnyTag) {
        // RDKit❌❌:     return v.value.a->type() == typeid(T);
        // RDKit❌❌:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   if (v.getTag() == tag) {
        // RDKit✔️✔️:     return true;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   return false;
        // RDKit✔️✔️: }
        // Behavior: the exact StringVector tag borrows the original ordered
        // byte strings. Every other modeled tag returns expected/actual kind;
        // independent Any type introspection is absent from the model.
        // Complexity: O(1) discriminant check and shared borrow. No allocation,
        // per-element traversal, decoding, coercion or auxiliary authority.
        match self {
            Self::StringVector(value) => Ok(value),
            _ => Err(self.wrong_kind(PropertyValueKind::StringVector)),
        }
    }

    /// Copy the exact source string-vector payload into independent owners.
    #[doc(hidden)]
    pub fn string_vector_clone(&self) -> Result<Vec<PropertyText>, PropertyValueError> {
        // RDKit✔️✔️: inline std::vector<std::string> rdvalue_cast<std::vector<std::string>>(
        // RDKit✔️✔️:     RDValue_cast_t v) {
        // RDKit✔️✔️:   if (rdvalue_is<std::vector<std::string>>(v)) {
        // RDKit✔️✔️:     return *v.ptrCast<std::vector<std::string>>();
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   throw std::bad_any_cast();
        // RDKit✔️✔️: }
        // RDKit✔️✔️: template <>
        // RDKit✔️✔️: inline short GetTag<std::vector<std::string>>() {
        // RDKit✔️✔️:   return VecStringTag;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: short getTag() const { return type; }
        // RDKit✔️✔️: template <class T>
        // RDKit✔️✔️:   inline T *ptrCast() const {
        // RDKit✔️✔️:     return RDTypeTag::detail::valuePtrCast<T>(value);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: template <>
        // RDKit✔️✔️: inline std::vector<std::string> *valuePtrCast<std::vector<std::string>>(
        // RDKit✔️✔️:     Value value) {
        // RDKit✔️✔️:   return value.vs;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: template <class T>
        // RDKit✔️✔️: inline bool rdvalue_is(RDValue_cast_t v) {
        // RDKit✔️✔️:   const short tag =
        // RDKit✔️✔️:       RDTypeTag::GetTag<typename boost::remove_reference<T>::type>();
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // If we are an Any tag, check the any type info
        // RDKit✔️✔️:   //  see the template specialization below if we are
        // RDKit✔️✔️:   //  looking for a boost any directly
        // RDKit❌❌:   if (v.getTag() == RDTypeTag::AnyTag) {
        // RDKit❌❌:     return v.value.a->type() == typeid(T);
        // RDKit❌❌:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   if (v.getTag() == tag) {
        // RDKit✔️✔️:     return true;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   return false;
        // RDKit✔️✔️: }
        // Behavior: reuse the canonical checked StringVector borrow, then
        // deep-copy its vector and each counted string. Order, duplicates,
        // empty strings and arbitrary bytes survive without decoding. Wrong
        // modeled tags retain the original structured expected/actual error.
        // Complexity: O(elements + total bytes), one vector and per-string
        // buffers as source return-by-value; O(1) failure with no allocation.
        self.as_string_vector().map(<[PropertyText]>::to_vec)
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

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum PropertyStoreError {
    EmptyKey,
    ComputedListKind(PropertyValueError),
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[doc(hidden)]
pub struct MissingPropertyError {
    key: PropertyText,
}

impl MissingPropertyError {
    #[must_use]
    pub fn key(&self) -> &PropertyText {
        &self.key
    }
}

impl std::fmt::Display for MissingPropertyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "missing property key bytes {:?}", self.key.as_bytes())
    }
}

impl std::error::Error for MissingPropertyError {}

/// One canonical typed value map with source insertion order and computed state.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct PropertyStore {
    values: BTreeMap<PropertyText, PropertyValue>,
    order: Vec<PropertyText>,
}

impl PropertyStore {
    pub(crate) const fn new() -> Self {
        Self {
            values: BTreeMap::new(),
            order: Vec::new(),
        }
    }

    pub(crate) fn from_records(
        records: impl IntoIterator<Item = (PropertyText, PropertyValue)>,
    ) -> Self {
        // Detached transport consumes explicit record order, preserving a source
        // clone's ordered dictionary rather than inferring order from keys.
        let mut store = Self::new();
        for (key, value) in records {
            match store.values.entry(key) {
                std::collections::btree_map::Entry::Occupied(mut entry) => {
                    entry.insert(value);
                }
                std::collections::btree_map::Entry::Vacant(entry) => {
                    store.order.push(entry.key().clone());
                    entry.insert(value);
                }
            }
        }
        store
    }

    pub(crate) fn values(&self) -> &BTreeMap<PropertyText, PropertyValue> {
        &self.values
    }

    pub(crate) fn get(&self, key: &[u8]) -> Option<&PropertyValue> {
        // BEGIN COMPLETE PINNED SF383
        // RDKit✔️🔝: bool getPropIfPresent(const std::string_view key, T &res) const {
        // RDKit✔️🔝:     return d_props.getValIfPresent(key, res);
        // RDKit✔️🔝:   }
        // END COMPLETE PINNED SF383
        // BEGIN COMPLETE REACHED Dict::getValIfPresent
        // RDKit✔️🔝:   template <typename T>
        // RDKit✔️🔝:   bool getValIfPresent(const std::string_view what, T &res) const {
        // RDKit✔️🔝:     for (const auto &data : _data) {
        // RDKit✔️🔝:       if (data.key == what) {
        // RDKit✔️🔝:         res = from_rdvalue<T>(data.val);
        // RDKit✔️🔝:         return true;
        // RDKit✔️🔝:       }
        // RDKit✔️🔝:     }
        // RDKit✔️🔝:     return false;
        // RDKit✔️🔝:   }
        // RDKit✔️🔝:
        // END COMPLETE REACHED Dict::getValIfPresent
        // The canonical counted-byte lookup distinguishes absence from every
        // present tagged value. Reached typed conversions remain in the one
        // CORE owner and execute only on present values, in caller order.
        // Lookup O(log P) replaces the source Dict O(P) scan without copies.
        // BEGIN RDKIT CPP FUNCTION Dict::hasVal
        // RDKit✔️🔝: bool hasVal(const std::string_view what) const {
        // RDKit✔️🔝:     for (const auto &data : _data) {
        // RDKit✔️🔝:       if (data.key == what) {
        // RDKit✔️🔝:         return true;
        // RDKit✔️🔝:       }
        // RDKit✔️🔝:     }
        // RDKit✔️🔝:     return false;
        // RDKit✔️🔝:   }
        // END RDKIT CPP FUNCTION Dict::hasVal
        // Behavior: get(...).is_some() is precisely source key presence.
        // Presence never converts/validates the payload or uses its truth
        // value; zero/false/empty String values remain present.
        // Complexity: one O(log P) byte-key lookup, no clone or allocation,
        // including on absence. The borrow also serves canonical value reads.
        self.values.get(key)
    }

    pub(crate) fn get_required(
        &self,
        key: impl AsRef<[u8]>,
    ) -> Result<&PropertyValue, MissingPropertyError> {
        // BEGIN RDKIT CPP FUNCTION Dict::getRDValue
        // RDKit✔️🔝: const RDValue &getRDValue(const std::string_view what) const {
        // RDKit✔️🔝:     for (const auto &data : _data) {
        // RDKit✔️🔝:       if (data.key == what) {
        // RDKit✔️🔝:         return data.val;
        // RDKit✔️🔝:       }
        // RDKit✔️🔝:     }
        // RDKit✔️🔝:     throw KeyErrorException(what);
        // RDKit✔️🔝:   }
        // END RDKIT CPP FUNCTION Dict::getRDValue
        // Behavior: return the same tagged borrowed payload, or a structural
        // missing-key error retaining the exact requested bytes. Optional
        // source getValIfPresent/hasVal lookup remains the separate get seam.
        // Complexity: O(log P) tree lookup replaces the source O(P) scan.
        // The key is copied only on failure, as in the owning source error.
        let key = key.as_ref();
        self.get(key).ok_or_else(|| MissingPropertyError {
            key: PropertyText::from_bytes(key),
        })
    }

    pub(crate) fn computed_names(&self) -> Result<Option<&[PropertyText]>, PropertyValueError> {
        self.get(b"__computedProps")
            .map(PropertyValue::as_string_vector)
            .transpose()
    }

    pub(crate) fn is_computed(&self, key: impl AsRef<[u8]>) -> Result<bool, PropertyValueError> {
        let key = key.as_ref();
        Ok(self
            .computed_names()?
            .is_some_and(|names| names.iter().any(|name| name.as_bytes() == key)))
    }

    pub(crate) fn ordered(
        &self,
    ) -> impl ExactSizeIterator<Item = (&PropertyText, &PropertyValue)> + '_ {
        // RDKit✔️❌: STR_VECT keys() const {
        // RDKit✔️❌:     STR_VECT res;
        // RDKit✔️❌:     res.reserve(_data.size());
        // RDKit✔️❌:     for (const auto &item : _data) {
        // RDKit✔️❌:       res.push_back(item.key);
        // RDKit✔️❌:     }
        // RDKit✔️❌:     return res;
        // RDKit✔️❌:   }
        // Behavior: keys follow source dictionary insertion order, including
        // the actual __computedProps entry. Overwrite retains position; erase
        // removes it; reinsert appends. No lexical sort, computed filtering,
        // text decoding or hidden reserved-entry omission occurs here.
        // This canonical read seam additionally borrows each corresponding
        // value. Rust's shared borrow keeps the dictionary stable throughout
        // iteration, rather than building the source's copied key list.
        // Complexity: no temporary key-vector/buffer copies, but each value
        // is resolved through the tree: O(P log P), versus source keys' linear
        // traversal/copy. That asymptotic value-lookup cost is recorded as a
        // known gap, not a claim of source complexity equivalence.
        self.order.iter().map(|key| {
            let value = self
                .values
                .get(key)
                .expect("private property order must match its canonical value map");
            (key, value)
        })
    }

    pub(crate) fn filtered_ordered(
        &self,
        include_private: bool,
        include_computed: bool,
    ) -> Result<impl Iterator<Item = (&PropertyText, &PropertyValue)> + '_, PropertyStoreError>
    {
        // RDKit✔️❌: STR_VECT getPropList(bool includePrivate = true,
        // RDKit✔️❌:                        bool includeComputed = true) const {
        // RDKit✔️❌:     const STR_VECT &tmp = d_props.keys();
        // RDKit✔️❌:     STR_VECT res, computed;
        // RDKit✔️❌:     if (!includeComputed &&
        // RDKit✔️❌:         getPropIfPresent(RDKit::detail::computedPropName, computed)) {
        // RDKit✔️❌:       computed.emplace_back(RDKit::detail::computedPropName);
        // RDKit✔️❌:     }
        // RDKit✔️❌:
        // RDKit✔️❌:     auto pos = tmp.begin();
        // RDKit✔️❌:     while (pos != tmp.end()) {
        // RDKit✔️❌:       if ((includePrivate || (*pos)[0] != '_') &&
        // RDKit✔️❌:           std::find(computed.begin(), computed.end(), *pos) == computed.end()) {
        // RDKit✔️❌:         res.push_back(*pos);
        // RDKit✔️❌:       }
        // RDKit✔️❌:       ++pos;
        // RDKit✔️❌:     }
        // RDKit✔️❌:     return res;
        // RDKit✔️❌:   }
        // RDKit✔️✔️: template <typename T>
        // RDKit✔️✔️:   bool getPropIfPresent(const std::string_view key, T &res) const {
        // RDKit✔️✔️:     return d_props.getValIfPresent(key, res);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: template <typename T>
        // RDKit✔️✔️:   bool getValIfPresent(const std::string_view what, T &res) const {
        // RDKit✔️✔️:     for (const auto &data : _data) {
        // RDKit✔️✔️:       if (data.key == what) {
        // RDKit✔️✔️:         res = from_rdvalue<T>(data.val);
        // RDKit✔️✔️:         return true;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     return false;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: template <class T>
        // RDKit✔️✔️: typename boost::disable_if<boost::is_arithmetic<T>, T>::type from_rdvalue(
        // RDKit✔️✔️:     RDValue_cast_t arg) {
        // RDKit✔️✔️:   return rdvalue_cast<T>(arg);
        // RDKit✔️✔️: }
        // Behavior: skip the reserved-vector read entirely when computed
        // values are included. Otherwise copy/cast its real value before
        // traversing/filtering any key: wrong kind must fail even if every
        // available key is private. A present vector also excludes its own
        // reserved name, including when empty. Exact byte equality and first
        // underscore byte implement the source membership/private conditions.
        // Complexity: the source computed-vector copy and linear membership
        // scan are retained. The shared ordered seam omits temporary key/output
        // copies, but O(P log P) value resolution remains a known cost gap.
        let mut computed = Vec::new();
        if !include_computed && let Some(value) = self.get(b"__computedProps") {
            computed = value
                .string_vector_clone()
                .map_err(PropertyStoreError::ComputedListKind)?;
            computed.push("__computedProps".into());
        }
        Ok(self.ordered().filter(move |(key, _)| {
            (include_private || key.as_bytes().first() != Some(&b'_')) && !computed.contains(key)
        }))
    }

    pub(crate) fn set(
        &mut self,
        key: PropertyText,
        value: PropertyValue,
    ) -> Result<(), PropertyStoreError> {
        // BEGIN RDKIT CPP FUNCTION Dict::Pair::Pair(std::string_view,const RDValue&)
        // RDKit✔️❌: Pair(std::string_view s, const RDValue &v) : key(std::string(s)), val(v) {}
        // END RDKIT CPP FUNCTION Dict::Pair::Pair(std::string_view,const RDValue&)
        // Behavior: canonical keys own the exact counted bytes, including NUL
        // and non-UTF-8. The tagged payload moves without conversion.
        // Complexity: key ownership is O(key bytes), but the existing tree plus
        // insertion index owns an additional copy of each new key compared
        // with the source Pair vector. This extra allocation is recorded as
        // a performance gap; no second property value is created.
        // BEGIN RDKIT CPP FUNCTION Dict::setVal<T&>
        // RDKit✔️🔝: void setVal(const std::string_view what, T &val) {
        // RDKit✔️✔️:     static_assert(!std::is_same_v<T, std::string_view>,
        // RDKit✔️✔️:                   "T cannot be string_view");
        // RDKit✔️🔝:     if (what.empty()) {
        // RDKit✔️🔝:       throw ValueErrorException("Cannot set value with empty key");
        // RDKit✔️🔝:     }
        // RDKit✔️✔️:     _hasNonPodData = true;
        // RDKit✔️🔝:     for (auto &&data : _data) {
        // RDKit✔️🔝:       if (data.key == what) {
        // RDKit✔️🔝:         RDValue::cleanup_rdvalue(data.val);
        // RDKit✔️🔝:         data.val = val;
        // RDKit✔️🔝:         return;
        // RDKit✔️🔝:       }
        // RDKit✔️🔝:     }
        // RDKit✔️🔝:     _data.push_back(Pair(what, val));
        // RDKit✔️🔝:   }
        // END RDKIT CPP FUNCTION Dict::setVal<T&>
        // The tree replaces the source linear value lookup. The order vector
        // stores only keys, so overwrites preserve position without duplicating
        // any String/Int/Double/Bool value.
        // BEGIN RDKIT CPP HELPERS RDValue::destroy / cleanup_rdvalue
        // RDKit✔️✔️:   void destroy() {
        // RDKit✔️✔️:     switch (type) {
        // RDKit✔️✔️:       case RDTypeTag::StringTag:
        // RDKit✔️✔️:         delete value.s;
        // RDKit✔️✔️:         break;
        // RDKit❌❌:       case RDTypeTag::AnyTag:
        // RDKit❌❌:         delete value.a;
        // RDKit❌❌:         break;
        // RDKit❌❌:       case RDTypeTag::VecDoubleTag:
        // RDKit❌❌:         delete value.vd;
        // RDKit❌❌:         break;
        // RDKit❌❌:       case RDTypeTag::VecFloatTag:
        // RDKit❌❌:         delete value.vf;
        // RDKit❌❌:         break;
        // RDKit✔️✔️:       case RDTypeTag::VecIntTag:
        // RDKit✔️✔️:         delete value.vi;
        // RDKit✔️✔️:         break;
        // RDKit❌❌:       case RDTypeTag::VecUnsignedIntTag:
        // RDKit❌❌:         delete value.vu;
        // RDKit❌❌:         break;
        // RDKit❌❌:       case RDTypeTag::VecStringTag:
        // RDKit❌❌:         delete value.vs;
        // RDKit❌❌:         break;
        // RDKit✔️✔️:       default:
        // RDKit✔️✔️:         break;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     type = RDTypeTag::EmptyTag;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   static  // Given a type and an RDAnyValue - delete the appropriate structure
        // RDKit✔️✔️:       inline void
        // RDKit✔️✔️:       cleanup_rdvalue(RDValue &rdvalue) {
        // RDKit✔️✔️:     rdvalue.destroy();
        // RDKit✔️✔️:   }
        // END RDKIT CPP HELPERS RDValue::destroy / cleanup_rdvalue
        // Rust's owning payload destroys replaced non-POD values by tag;
        // no mutable source cleanup flag is needed. PropertyValue cannot
        // borrow a string_view, so the source static assertion holds by type.
        // Reject before changing either ordered keys or the canonical map.
        if key.is_empty() {
            return Err(PropertyStoreError::EmptyKey);
        }
        match self.values.entry(key) {
            std::collections::btree_map::Entry::Occupied(mut entry) => {
                entry.insert(value);
            }
            std::collections::btree_map::Entry::Vacant(entry) => {
                self.order.push(entry.key().clone());
                entry.insert(value);
            }
        }
        Ok(())
    }

    pub(crate) fn set_computed(
        &mut self,
        key: PropertyText,
        value: PropertyValue,
    ) -> Result<(), PropertyStoreError> {
        // RDKit✔️❌: void setProp(const std::string_view key, T val, bool computed = false) const {
        // RDKit✔️❌:     if(key.empty()) {
        // RDKit✔️❌:       throw ValueErrorException("Cannot set property with empty key");
        // RDKit✔️❌:     }
        // RDKit✔️❌:     if (computed) {
        // RDKit✔️❌:       STR_VECT compLst;
        // RDKit✔️❌:       getPropIfPresent(RDKit::detail::computedPropName, compLst);
        // RDKit✔️❌:       if (std::find(compLst.begin(), compLst.end(), key) == compLst.end()) {
        // RDKit✔️❌:         compLst.emplace_back(key);
        // RDKit✔️❌:         d_props.setVal(RDKit::detail::computedPropName, compLst);
        // RDKit✔️❌:       }
        // RDKit✔️❌:     }
        // RDKit✔️❌:     d_props.setVal(key, val);
        // RDKit✔️❌:   }
        // RDKit✔️✔️: template <typename T>
        // RDKit✔️✔️:   bool getPropIfPresent(const std::string_view key, T &res) const {
        // RDKit✔️✔️:     return d_props.getValIfPresent(key, res);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: template <typename T>
        // RDKit✔️✔️:   bool getValIfPresent(const std::string_view what, T &res) const {
        // RDKit✔️✔️:     for (const auto &data : _data) {
        // RDKit✔️✔️:       if (data.key == what) {
        // RDKit✔️✔️:         res = from_rdvalue<T>(data.val);
        // RDKit✔️✔️:         return true;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     return false;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: template <class T>
        // RDKit✔️✔️: typename boost::disable_if<boost::is_arithmetic<T>, T>::type from_rdvalue(
        // RDKit✔️✔️:     RDValue_cast_t arg) {
        // RDKit✔️✔️:   return rdvalue_cast<T>(arg);
        // RDKit✔️✔️: }
        // Behavior: reject an empty key first. A present reserved entry is
        // actually read as StringVector and copied before any write, so wrong
        // kind propagates with no ordinary-value mutation. Append a missing
        // membership and store the reserved vector before the requested value.
        // Ordinary writes retain membership; a write to __computedProps itself
        // may overwrite its vector with another modeled kind, as in source.
        // Complexity: source linear membership scan and deep vector/string
        // copy retained; tree lookups O(log P). The insertion-order index owns
        // extra key bytes versus source Pair vector, recorded as a cost gap.
        if key.is_empty() {
            return Err(PropertyStoreError::EmptyKey);
        }
        let mut names = match self.get(b"__computedProps") {
            Some(value) => value
                .string_vector_clone()
                .map_err(PropertyStoreError::ComputedListKind)?,
            None => Vec::new(),
        };
        if !names.contains(&key) {
            names.push(key.clone());
            // setVal copies this borrowed vector. Keep that second source copy.
            self.set(
                "__computedProps".into(),
                PropertyValue::from(names.as_slice()),
            )?;
        }
        self.set(key, value)
    }

    pub(crate) fn clear_value(&mut self, key: impl AsRef<[u8]>) {
        // RDKit✔️❌: void clearVal(const std::string_view what) {
        // RDKit✔️❌:     for (auto it = _data.begin(); it < _data.end(); ++it) {
        // RDKit✔️❌:       if (it->key == what) {
        // RDKit✔️❌:         if (_hasNonPodData) {
        // RDKit✔️❌:           RDValue::cleanup_rdvalue(it->val);
        // RDKit✔️❌:         }
        // RDKit✔️❌:         _data.erase(it);
        // RDKit✔️❌:         return;
        // RDKit✔️❌:       }
        // RDKit✔️❌:     }
        // RDKit✔️❌:   }
        // RDKit✔️✔️:   void destroy() {
        // RDKit✔️✔️:     switch (type) {
        // RDKit✔️✔️:       case RDTypeTag::StringTag:
        // RDKit✔️✔️:         delete value.s;
        // RDKit✔️✔️:         break;
        // RDKit❌❌:       case RDTypeTag::AnyTag:
        // RDKit❌❌:         delete value.a;
        // RDKit❌❌:         break;
        // RDKit❌❌:       case RDTypeTag::VecDoubleTag:
        // RDKit❌❌:         delete value.vd;
        // RDKit❌❌:         break;
        // RDKit❌❌:       case RDTypeTag::VecFloatTag:
        // RDKit❌❌:         delete value.vf;
        // RDKit❌❌:         break;
        // RDKit✔️✔️:       case RDTypeTag::VecIntTag:
        // RDKit✔️✔️:         delete value.vi;
        // RDKit✔️✔️:         break;
        // RDKit❌❌:       case RDTypeTag::VecUnsignedIntTag:
        // RDKit❌❌:         delete value.vu;
        // RDKit❌❌:         break;
        // RDKit✔️✔️:       case RDTypeTag::VecStringTag:
        // RDKit✔️✔️:         delete value.vs;
        // RDKit✔️✔️:         break;
        // RDKit✔️✔️:       default:
        // RDKit✔️✔️:         break;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     type = RDTypeTag::EmptyTag;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   static  // Given a type and an RDAnyValue - delete the appropriate structure
        // RDKit✔️✔️:       inline void
        // RDKit✔️✔️:       cleanup_rdvalue(RDValue &rdvalue) {
        // RDKit✔️✔️:     rdvalue.destroy();
        // RDKit✔️✔️:   }
        // Behavior: the dictionary primitive removes only the matching key,
        // owning payload and its insertion-order position. Missing keys do
        // nothing; empty keys are permitted; no computed-vector read/update
        // belongs here. Rust drop destroys the removed modeled payload once.
        // Complexity: logarithmic tree removal followed by one ordered-key
        // search/erase is O(P), as source linear find/erase. The separate key
        // index requires additional key storage/destruction, a known cost gap.
        let key = key.as_ref();
        if self.values.remove(key).is_some()
            && let Some(position) = self.order.iter().position(|name| name.as_bytes() == key)
        {
            self.order.remove(position);
        }
    }

    pub(crate) fn clear(&mut self, key: impl AsRef<[u8]>) -> Result<(), PropertyStoreError> {
        // RDKit✔️❌: void clearProp(const std::string_view key) const {
        // RDKit✔️❌:     STR_VECT compLst;
        // RDKit✔️❌:     if (getPropIfPresent(RDKit::detail::computedPropName, compLst)) {
        // RDKit✔️❌:       auto svi = std::find(compLst.begin(), compLst.end(), key);
        // RDKit✔️❌:       if (svi != compLst.end()) {
        // RDKit✔️❌:         compLst.erase(svi);
        // RDKit✔️❌:         d_props.setVal(RDKit::detail::computedPropName, compLst);
        // RDKit✔️❌:       }
        // RDKit✔️❌:     }
        // RDKit✔️❌:     d_props.clearVal(key);
        // RDKit✔️❌:   }
        // RDKit✔️✔️: template <typename T>
        // RDKit✔️✔️:   bool getPropIfPresent(const std::string_view key, T &res) const {
        // RDKit✔️✔️:     return d_props.getValIfPresent(key, res);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: template <typename T>
        // RDKit✔️✔️:   bool getValIfPresent(const std::string_view what, T &res) const {
        // RDKit✔️✔️:     for (const auto &data : _data) {
        // RDKit✔️✔️:       if (data.key == what) {
        // RDKit✔️✔️:         res = from_rdvalue<T>(data.val);
        // RDKit✔️✔️:         return true;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     return false;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: template <class T>
        // RDKit✔️✔️: typename boost::disable_if<boost::is_arithmetic<T>, T>::type from_rdvalue(
        // RDKit✔️✔️:     RDValue_cast_t arg) {
        // RDKit✔️✔️:   return rdvalue_cast<T>(arg);
        // RDKit✔️✔️: }
        // Behavior: read/copy a present real StringVector before any mutation;
        // wrong-kind failure also occurs for an absent requested key. Remove
        // only the first matching membership, rewrite its vector only then,
        // and finally clear the requested dictionary entry. An empty key is
        // not rejected; clearing __computedProps removes its entire entry.
        // Complexity: source vector copy, linear find/erase and setter copy
        // retained. Tree lookup plus an ordered-key erase has the source's
        // linear erase shape, but its extra key index retains a known cost gap.
        let key = key.as_ref();
        if let Some(value) = self.get(b"__computedProps") {
            let mut names = value
                .string_vector_clone()
                .map_err(PropertyStoreError::ComputedListKind)?;
            if let Some(position) = names.iter().position(|name| name.as_bytes() == key) {
                names.remove(position);
                self.set(
                    "__computedProps".into(),
                    PropertyValue::from(names.as_slice()),
                )?;
            }
        }
        self.clear_value(key);
        Ok(())
    }

    pub(crate) fn clear_computed(&mut self) -> Result<(), PropertyStoreError> {
        // RDKit✔️❌: void clearComputedProps() const {
        // RDKit✔️❌:     STR_VECT compLst;
        // RDKit✔️❌:     if (getPropIfPresent(RDKit::detail::computedPropName, compLst) &&
        // RDKit✔️❌:         !compLst.empty()) {
        // RDKit✔️❌:       for (const auto &sv : compLst) {
        // RDKit✔️❌:         d_props.clearVal(sv);
        // RDKit✔️❌:       }
        // RDKit✔️❌:       compLst.clear();
        // RDKit✔️❌:       d_props.setVal(RDKit::detail::computedPropName, compLst);
        // RDKit✔️❌:     }
        // RDKit✔️❌:   }
        // RDKit✔️✔️: template <typename T>
        // RDKit✔️✔️:   bool getPropIfPresent(const std::string_view key, T &res) const {
        // RDKit✔️✔️:     return d_props.getValIfPresent(key, res);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: template <typename T>
        // RDKit✔️✔️:   bool getValIfPresent(const std::string_view what, T &res) const {
        // RDKit✔️✔️:     for (const auto &data : _data) {
        // RDKit✔️✔️:       if (data.key == what) {
        // RDKit✔️✔️:         res = from_rdvalue<T>(data.val);
        // RDKit✔️✔️:         return true;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     return false;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: template <class T>
        // RDKit✔️✔️: typename boost::disable_if<boost::is_arithmetic<T>, T>::type from_rdvalue(
        // RDKit✔️✔️:     RDValue_cast_t arg) {
        // RDKit✔️✔️:   return rdvalue_cast<T>(arg);
        // RDKit✔️✔️: }
        // Behavior: a present reserved value must genuinely copy/cast as
        // StringVector. Absence and present-empty remain distinct and untouched.
        // Clear names in copied vector order, including duplicates and the
        // reserved key itself, before writing a present empty reserved vector.
        // If that key was erased in the loop, its reinsertion is at the end.
        // Wrong-kind errors precede every dictionary mutation.
        // Complexity: the same source ordered copy and per-key clear loop;
        // no sorted set or bulk retain changes source effects. Empty copy into
        // the stored value needs no element allocation. Additional map/index
        // key ownership remains a known cost gap versus source Pair storage.
        if let Some(value) = self.get(b"__computedProps") {
            let names = value
                .string_vector_clone()
                .map_err(PropertyStoreError::ComputedListKind)?;
            if !names.is_empty() {
                for name in &names {
                    self.clear_value(name);
                }
                self.set(
                    "__computedProps".into(),
                    PropertyValue::StringVector(Vec::new()),
                )?;
            }
        }
        Ok(())
    }

    pub(crate) fn update_from(&mut self, source: &Self, preserve_existing: bool) {
        // RDKit✔️❌: void update(const Dict &other, bool preserveExisting = false) {
        // RDKit✔️❌:     if (!preserveExisting) {
        // RDKit✔️❌:       *this = other;
        // RDKit✔️❌:     } else {
        // RDKit✔️❌:       if (other._hasNonPodData) {
        // RDKit✔️❌:         _hasNonPodData = true;
        // RDKit✔️❌:       }
        // RDKit✔️❌:       for (const auto &opair : other._data) {
        // RDKit✔️❌:         Pair *target = nullptr;
        // RDKit✔️❌:         for (auto &dpair : _data) {
        // RDKit✔️❌:           if (dpair.key == opair.key) {
        // RDKit✔️❌:             target = &dpair;
        // RDKit✔️❌:             break;
        // RDKit✔️❌:           }
        // RDKit✔️❌:         }
        // RDKit✔️❌:
        // RDKit✔️❌:         if (!target) {
        // RDKit✔️❌:           // need to create blank entry and copy
        // RDKit✔️❌:           _data.push_back(Pair(opair.key));
        // RDKit✔️❌:           copy_rdvalue(_data.back().val, opair.val);
        // RDKit✔️❌:         } else {
        // RDKit✔️❌:           // just copy
        // RDKit✔️❌:           copy_rdvalue(target->val, opair.val);
        // RDKit✔️❌:         }
        // RDKit✔️❌:       }
        // RDKit✔️❌:     }
        // RDKit✔️❌:   }
        // RDKit✔️❌: Dict &operator=(const Dict &other) {
        // RDKit✔️❌:     if (this == &other) {
        // RDKit✔️❌:       return *this;
        // RDKit✔️❌:     }
        // RDKit✔️❌:     if (_hasNonPodData) {
        // RDKit✔️❌:       reset();
        // RDKit✔️❌:     }
        // RDKit✔️❌:
        // RDKit✔️❌:     if (other._hasNonPodData) {
        // RDKit✔️❌:       std::vector<Pair> data(other._data.size());
        // RDKit✔️❌:       _data.swap(data);
        // RDKit✔️❌:       for (size_t i = 0; i < _data.size(); ++i) {
        // RDKit✔️❌:         _data[i].key = other._data[i].key;
        // RDKit✔️❌:         copy_rdvalue(_data[i].val, other._data[i].val);
        // RDKit✔️❌:       }
        // RDKit✔️❌:     } else {
        // RDKit✔️❌:       _data = other._data;
        // RDKit✔️❌:     }
        // RDKit✔️❌:     _hasNonPodData = other._hasNonPodData;
        // RDKit✔️❌:     return *this;
        // RDKit✔️❌:   }
        // RDKit✔️❌: void reset() {
        // RDKit✔️❌:     if (_hasNonPodData) {
        // RDKit✔️❌:       for (auto &&data : _data) {
        // RDKit✔️❌:         RDValue::cleanup_rdvalue(data.val);
        // RDKit✔️❌:       }
        // RDKit✔️❌:     }
        // RDKit✔️❌:     DataType data;
        // RDKit✔️❌:     _data.swap(data);
        // RDKit✔️❌:   }
        // RDKit✔️❌: inline void copy_rdvalue(RDValue &dest, const RDValue &src) {
        // RDKit✔️❌:   if (&dest == &src) {  // don't copy over yourself
        // RDKit✔️❌:     return;
        // RDKit✔️❌:   }
        // RDKit✔️❌:   dest.destroy();
        // RDKit✔️❌:   dest.type = src.type;
        // RDKit✔️❌:   switch (src.type) {
        // RDKit✔️✔️:     case RDTypeTag::StringTag:
        // RDKit✔️✔️:       dest.value.s = new std::string(*src.value.s);
        // RDKit✔️✔️:       break;
        // RDKit❌❌:     case RDTypeTag::AnyTag:
        // RDKit❌❌:       dest.value.a = new std::any(*src.value.a);
        // RDKit❌❌:       break;
        // RDKit❌❌:     case RDTypeTag::VecDoubleTag:
        // RDKit❌❌:       dest.value.vd = new std::vector<double>(*src.value.vd);
        // RDKit❌❌:       break;
        // RDKit❌❌:     case RDTypeTag::VecFloatTag:
        // RDKit❌❌:       dest.value.vf = new std::vector<float>(*src.value.vf);
        // RDKit❌❌:       break;
        // RDKit✔️✔️:     case RDTypeTag::VecIntTag:
        // RDKit✔️✔️:       dest.value.vi = new std::vector<int>(*src.value.vi);
        // RDKit✔️✔️:       break;
        // RDKit❌❌:     case RDTypeTag::VecUnsignedIntTag:
        // RDKit❌❌:       dest.value.vu = new std::vector<unsigned int>(*src.value.vu);
        // RDKit❌❌:       break;
        // RDKit✔️✔️:     case RDTypeTag::VecStringTag:
        // RDKit✔️✔️:       dest.value.vs = new std::vector<std::string>(*src.value.vs);
        // RDKit✔️✔️:       break;
        // RDKit✔️✔️:     default:
        // RDKit✔️✔️:       dest = src;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // Behavior: false replaces the complete ordered dictionary by deep
        // copy. True walks source records in their insertion order, overwrites
        // existing values in place, and appends new keys. Despite its name,
        // source preserveExisting does not keep a same-key destination value.
        // The real reserved computed entry copies normally: absent source
        // leaves it untouched in the merge branch; present-empty overwrites it.
        // Direct pair copying has no setVal empty-key precondition or typed
        // computed-list read, and Rust ownership cleans replaced values.
        // Safe mutable/shared borrows cannot alias the source dictionary.
        // Complexity: indexed merge lookup is O(S log D), improving source's
        // nested O(S*D) key search. Deep payload copies are source-required;
        // separate ordered key ownership and tree nodes add a known storage
        // and allocation cost, hence the overall cost-gap marker.
        if !preserve_existing {
            *self = source.clone();
            return;
        }
        for (key, value) in source.ordered() {
            // Dict::update copies pairs directly, without setVal's empty-key
            // precondition. Keep that distinction while transporting bytes.
            match self.values.entry(key.clone()) {
                std::collections::btree_map::Entry::Occupied(mut entry) => {
                    entry.insert(value.clone());
                }
                std::collections::btree_map::Entry::Vacant(entry) => {
                    self.order.push(entry.key().clone());
                    entry.insert(value.clone());
                }
            }
        }
    }

    #[cfg(test)]
    pub(crate) fn ordered_keys(&self) -> &[PropertyText] {
        &self.order
    }
}

impl Default for PropertyStore {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    // These original fixtures contain UTF-8 literals. Decode only their
    // borrowed test projection; raw-byte controls assert as_bytes directly.
    // Invalid bytes fail this assertion, never change a chemistry outcome.
    fn fixture_text(value: &crate::PropertyText) -> &str {
        std::str::from_utf8(value.as_bytes()).expect("unchanged UTF-8 fixture bytes")
    }

    use super::*;
    use crate::{Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Element};

    #[test]
    fn typed_property_value_preserves_all_variants_and_exact_access_errors() {
        let values = [
            PropertyValue::String("7".into()),
            PropertyValue::Int(7),
            PropertyValue::Double(7.0),
            PropertyValue::Bool(true),
        ];
        assert_eq!(values[0].kind(), PropertyValueKind::String);
        assert_eq!(values[1].kind(), PropertyValueKind::Int);
        assert_eq!(values[2].kind(), PropertyValueKind::Double);
        assert_eq!(values[3].kind(), PropertyValueKind::Bool);
        assert_eq!(values[0].as_string().map(fixture_text), Ok("7"));
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
        store.set("z".into(), PropertyValue::String("7".into()));
        store.set("a".into(), PropertyValue::Int(7));
        assert_eq!(
            store
                .ordered_keys()
                .iter()
                .map(fixture_text)
                .collect::<Vec<_>>(),
            &["z", "a"]
        );
        store.set("z".into(), PropertyValue::Double(-0.0));
        assert_eq!(
            store
                .ordered_keys()
                .iter()
                .map(fixture_text)
                .collect::<Vec<_>>(),
            &["z", "a"]
        );
        assert_eq!(
            store.get("z".as_bytes()),
            Some(&PropertyValue::Double(-0.0))
        );

        store.clear("z");
        assert_eq!(
            store
                .ordered_keys()
                .iter()
                .map(fixture_text)
                .collect::<Vec<_>>(),
            &["a"]
        );
        store.set("z".into(), PropertyValue::Bool(false));
        assert_eq!(
            store
                .ordered_keys()
                .iter()
                .map(fixture_text)
                .collect::<Vec<_>>(),
            &["a", "z"]
        );
        store.set_computed("c".into(), PropertyValue::Double(1.25));
        store.set_computed("a".into(), PropertyValue::String("seven".into()));
        assert_eq!(
            store
                .ordered_keys()
                .iter()
                .map(fixture_text)
                .collect::<Vec<_>>(),
            &["a", "z", "__computedProps", "c"]
        );
        assert!(store.is_computed("a").unwrap());
        assert!(store.is_computed("c").unwrap());
        store.clear_computed();
        assert_eq!(
            store
                .ordered_keys()
                .iter()
                .map(fixture_text)
                .collect::<Vec<_>>(),
            &["z", "__computedProps"]
        );
        assert_eq!(store.get("z".as_bytes()), Some(&PropertyValue::Bool(false)));
        assert!(
            store
                .computed_names()
                .unwrap()
                .expect("computed list retained")
                .is_empty()
        );
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
            bond.set_computed_prop("", PropertyValue::String("bad".into()))
                .is_err()
        );
        assert_eq!(bond, bond_before);
    }
}
#[cfg(test)]
mod q01_b1_tests {
    fn fixture_text(value: &crate::PropertyText) -> &str {
        std::str::from_utf8(value.as_bytes()).expect("unchanged UTF-8 fixture bytes")
    }

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
        assert!(store.is_computed("b").unwrap());
        assert_eq!(
            store
                .ordered_keys()
                .iter()
                .map(fixture_text)
                .collect::<Vec<_>>(),
            ["a", "__computedProps", "b"]
        );
        store.clear("a");
        store.set("a".into(), Vec::<i32>::new().into());
        assert_eq!(
            store
                .ordered_keys()
                .iter()
                .map(fixture_text)
                .collect::<Vec<_>>(),
            ["__computedProps", "b", "a"]
        );
        store.clear_computed();
        assert_eq!(
            store
                .ordered_keys()
                .iter()
                .map(fixture_text)
                .collect::<Vec<_>>(),
            ["__computedProps", "a"]
        );
        assert_eq!(
            store.get("a".as_bytes()).unwrap().as_int_vector().unwrap(),
            [] as [i32; 0]
        );
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
    fn fixture_text(value: &crate::PropertyText) -> &str {
        std::str::from_utf8(value.as_bytes()).expect("unchanged UTF-8 fixture bytes")
    }

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
            assert!(store.is_computed("rank").unwrap());
            assert_eq!(
                store
                    .ordered_keys()
                    .iter()
                    .map(fixture_text)
                    .collect::<Vec<_>>(),
                ["a", "__computedProps", "rank"]
            );
            let saved = store.clone();
            store.clear_computed();
            assert_eq!(
                saved.get("rank".as_bytes()),
                Some(&PropertyValue::UInt(number))
            );
            assert_eq!(
                store
                    .ordered_keys()
                    .iter()
                    .map(fixture_text)
                    .collect::<Vec<_>>(),
                ["a", "__computedProps"]
            );
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
            assert!(query.is_prop_computed("rank").unwrap());
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
            assert!(query.bond().is_prop_computed("rank").unwrap());
        }
    }
}

#[cfg(test)]
mod uint_complete_source_condition_cells {
    fn fixture_text(value: &crate::PropertyText) -> &str {
        std::str::from_utf8(value.as_bytes()).expect("unchanged UTF-8 fixture bytes")
    }

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
        assert_eq!(
            store
                .ordered_keys()
                .iter()
                .map(fixture_text)
                .collect::<Vec<_>>(),
            ["a", "__computedProps", "rank"]
        );
        assert!(store.is_computed("rank").unwrap());
        let saved = store.clone();
        store.clear_computed();
        assert_eq!(
            saved.get("rank".as_bytes()),
            Some(&PropertyValue::UInt(0_u32))
        );
        assert_eq!(
            store
                .ordered_keys()
                .iter()
                .map(fixture_text)
                .collect::<Vec<_>>(),
            ["a", "__computedProps"]
        );
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
        assert!(query.is_prop_computed("rank").unwrap());
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
        assert!(querybond.bond().is_prop_computed("rank").unwrap());
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
        assert!(reordered.atoms[2].is_prop_computed("rank").unwrap());
        assert!(reordered.bonds[0].is_prop_computed("rank").unwrap());
        let mut edit = g.begin_batch_edit().unwrap();
        edit.remove_atom(AtomId::new(1)).unwrap();
        let (fragment, _) = edit.finish().unwrap();
        assert_eq!(fragment.atoms[0].prop("rank"), None);
        assert_eq!(fragment.bonds[0].prop("rank"), None);
        assert!(!fragment.atoms[0].is_prop_computed("rank").unwrap());
        assert!(!fragment.bonds[0].is_prop_computed("rank").unwrap());
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
        assert_eq!(
            store
                .ordered_keys()
                .iter()
                .map(fixture_text)
                .collect::<Vec<_>>(),
            ["a", "__computedProps", "rank"]
        );
        assert!(store.is_computed("rank").unwrap());
        let saved = store.clone();
        store.clear_computed();
        assert_eq!(
            saved.get("rank".as_bytes()),
            Some(&PropertyValue::UInt(1_u32))
        );
        assert_eq!(
            store
                .ordered_keys()
                .iter()
                .map(fixture_text)
                .collect::<Vec<_>>(),
            ["a", "__computedProps"]
        );
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
        assert!(query.is_prop_computed("rank").unwrap());
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
        assert!(querybond.bond().is_prop_computed("rank").unwrap());
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
        assert!(reordered.atoms[2].is_prop_computed("rank").unwrap());
        assert!(reordered.bonds[0].is_prop_computed("rank").unwrap());
        let mut edit = g.begin_batch_edit().unwrap();
        edit.remove_atom(AtomId::new(1)).unwrap();
        let (fragment, _) = edit.finish().unwrap();
        assert_eq!(fragment.atoms[0].prop("rank"), None);
        assert_eq!(fragment.bonds[0].prop("rank"), None);
        assert!(!fragment.atoms[0].is_prop_computed("rank").unwrap());
        assert!(!fragment.bonds[0].is_prop_computed("rank").unwrap());
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
        assert_eq!(
            store
                .ordered_keys()
                .iter()
                .map(fixture_text)
                .collect::<Vec<_>>(),
            ["a", "__computedProps", "rank"]
        );
        assert!(store.is_computed("rank").unwrap());
        let saved = store.clone();
        store.clear_computed();
        assert_eq!(
            saved.get("rank".as_bytes()),
            Some(&PropertyValue::UInt(2147483646_u32))
        );
        assert_eq!(
            store
                .ordered_keys()
                .iter()
                .map(fixture_text)
                .collect::<Vec<_>>(),
            ["a", "__computedProps"]
        );
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
        assert!(query.is_prop_computed("rank").unwrap());
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
        assert!(querybond.bond().is_prop_computed("rank").unwrap());
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
        assert!(reordered.atoms[2].is_prop_computed("rank").unwrap());
        assert!(reordered.bonds[0].is_prop_computed("rank").unwrap());
        let mut edit = g.begin_batch_edit().unwrap();
        edit.remove_atom(AtomId::new(1)).unwrap();
        let (fragment, _) = edit.finish().unwrap();
        assert_eq!(fragment.atoms[0].prop("rank"), None);
        assert_eq!(fragment.bonds[0].prop("rank"), None);
        assert!(!fragment.atoms[0].is_prop_computed("rank").unwrap());
        assert!(!fragment.bonds[0].is_prop_computed("rank").unwrap());
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
        assert_eq!(
            store
                .ordered_keys()
                .iter()
                .map(fixture_text)
                .collect::<Vec<_>>(),
            ["a", "__computedProps", "rank"]
        );
        assert!(store.is_computed("rank").unwrap());
        let saved = store.clone();
        store.clear_computed();
        assert_eq!(
            saved.get("rank".as_bytes()),
            Some(&PropertyValue::UInt(2147483647_u32))
        );
        assert_eq!(
            store
                .ordered_keys()
                .iter()
                .map(fixture_text)
                .collect::<Vec<_>>(),
            ["a", "__computedProps"]
        );
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
        assert!(query.is_prop_computed("rank").unwrap());
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
        assert!(querybond.bond().is_prop_computed("rank").unwrap());
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
        assert!(reordered.atoms[2].is_prop_computed("rank").unwrap());
        assert!(reordered.bonds[0].is_prop_computed("rank").unwrap());
        let mut edit = g.begin_batch_edit().unwrap();
        edit.remove_atom(AtomId::new(1)).unwrap();
        let (fragment, _) = edit.finish().unwrap();
        assert_eq!(fragment.atoms[0].prop("rank"), None);
        assert_eq!(fragment.bonds[0].prop("rank"), None);
        assert!(!fragment.atoms[0].is_prop_computed("rank").unwrap());
        assert!(!fragment.bonds[0].is_prop_computed("rank").unwrap());
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
        assert_eq!(
            store
                .ordered_keys()
                .iter()
                .map(fixture_text)
                .collect::<Vec<_>>(),
            ["a", "__computedProps", "rank"]
        );
        assert!(store.is_computed("rank").unwrap());
        let saved = store.clone();
        store.clear_computed();
        assert_eq!(
            saved.get("rank".as_bytes()),
            Some(&PropertyValue::UInt(2147483648_u32))
        );
        assert_eq!(
            store
                .ordered_keys()
                .iter()
                .map(fixture_text)
                .collect::<Vec<_>>(),
            ["a", "__computedProps"]
        );
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
        assert!(query.is_prop_computed("rank").unwrap());
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
        assert!(querybond.bond().is_prop_computed("rank").unwrap());
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
        assert!(reordered.atoms[2].is_prop_computed("rank").unwrap());
        assert!(reordered.bonds[0].is_prop_computed("rank").unwrap());
        let mut edit = g.begin_batch_edit().unwrap();
        edit.remove_atom(AtomId::new(1)).unwrap();
        let (fragment, _) = edit.finish().unwrap();
        assert_eq!(fragment.atoms[0].prop("rank"), None);
        assert_eq!(fragment.bonds[0].prop("rank"), None);
        assert!(!fragment.atoms[0].is_prop_computed("rank").unwrap());
        assert!(!fragment.bonds[0].is_prop_computed("rank").unwrap());
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
        assert_eq!(
            store
                .ordered_keys()
                .iter()
                .map(fixture_text)
                .collect::<Vec<_>>(),
            ["a", "__computedProps", "rank"]
        );
        assert!(store.is_computed("rank").unwrap());
        let saved = store.clone();
        store.clear_computed();
        assert_eq!(
            saved.get("rank".as_bytes()),
            Some(&PropertyValue::UInt(4294967295_u32))
        );
        assert_eq!(
            store
                .ordered_keys()
                .iter()
                .map(fixture_text)
                .collect::<Vec<_>>(),
            ["a", "__computedProps"]
        );
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
        assert!(query.is_prop_computed("rank").unwrap());
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
        assert!(querybond.bond().is_prop_computed("rank").unwrap());
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
        assert!(reordered.atoms[2].is_prop_computed("rank").unwrap());
        assert!(reordered.bonds[0].is_prop_computed("rank").unwrap());
        let mut edit = g.begin_batch_edit().unwrap();
        edit.remove_atom(AtomId::new(1)).unwrap();
        let (fragment, _) = edit.finish().unwrap();
        assert_eq!(fragment.atoms[0].prop("rank"), None);
        assert_eq!(fragment.bonds[0].prop("rank"), None);
        assert!(!fragment.atoms[0].is_prop_computed("rank").unwrap());
        assert!(!fragment.bonds[0].is_prop_computed("rank").unwrap());
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
