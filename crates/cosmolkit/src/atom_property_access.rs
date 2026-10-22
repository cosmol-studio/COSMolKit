//! Exact-type atom metadata access; all writes use the registered property operation.
use crate::{AtomId, Molecule, OperationError, PropertyValue, PropertyValueError};

/// A typed atom-property read failed; absence alone is returned as None.
#[derive(Debug, Clone, PartialEq)]
pub enum AtomPropertyReadError {
    Operation(OperationError),
    Type(PropertyValueError),
    Utf8(std::str::Utf8Error),
}
impl std::fmt::Display for AtomPropertyReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Operation(e) => e.fmt(f),
            Self::Type(e) => e.fmt(f),
            Self::Utf8(e) => e.fmt(f),
        }
    }
}
impl std::error::Error for AtomPropertyReadError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(match self {
            Self::Operation(e) => e,
            Self::Type(e) => e,
            Self::Utf8(e) => e,
        })
    }
}
impl From<OperationError> for AtomPropertyReadError {
    fn from(e: OperationError) -> Self {
        Self::Operation(e)
    }
}
impl From<PropertyValueError> for AtomPropertyReadError {
    fn from(e: PropertyValueError) -> Self {
        Self::Type(e)
    }
}
impl From<std::str::Utf8Error> for AtomPropertyReadError {
    fn from(e: std::str::Utf8Error) -> Self {
        Self::Utf8(e)
    }
}

impl Molecule {
    /// Read exactly the stored Int kind; missing keys return None.
    /// Wrong kinds, invalid atom IDs and invalid UTF-8 are errors.
    pub fn atom_property_int(
        &self,
        atom: AtomId,
        key: &str,
    ) -> Result<Option<i32>, AtomPropertyReadError> {
        self.atom_property(atom, key)?
            .map(|v| Ok(v.as_int()?))
            .transpose()
    }

    /// Return a new molecule with a Int property; leave the receiver unchanged.
    #[cfg(feature = "cap-transforms")]
    pub fn with_atom_property_int(
        &self,
        atom: AtomId,
        key: &str,
        value: i32,
    ) -> Result<Self, OperationError> {
        self.with_atom_property(atom, key, &PropertyValue::Int(value))
    }

    /// Set a Int property through the existing COW operation.
    #[cfg(feature = "cap-transforms")]
    pub fn set_atom_property_int_(
        &mut self,
        atom: AtomId,
        key: &str,
        value: i32,
    ) -> Result<(), OperationError> {
        self.set_atom_property_(atom, key, &PropertyValue::Int(value))
    }

    /// Read exactly the stored UInt kind; missing keys return None.
    /// Wrong kinds, invalid atom IDs and invalid UTF-8 are errors.
    pub fn atom_property_uint(
        &self,
        atom: AtomId,
        key: &str,
    ) -> Result<Option<u32>, AtomPropertyReadError> {
        self.atom_property(atom, key)?
            .map(|v| Ok(v.as_uint()?))
            .transpose()
    }

    /// Return a new molecule with a UInt property; leave the receiver unchanged.
    #[cfg(feature = "cap-transforms")]
    pub fn with_atom_property_uint(
        &self,
        atom: AtomId,
        key: &str,
        value: u32,
    ) -> Result<Self, OperationError> {
        self.with_atom_property(atom, key, &PropertyValue::UInt(value))
    }

    /// Set a UInt property through the existing COW operation.
    #[cfg(feature = "cap-transforms")]
    pub fn set_atom_property_uint_(
        &mut self,
        atom: AtomId,
        key: &str,
        value: u32,
    ) -> Result<(), OperationError> {
        self.set_atom_property_(atom, key, &PropertyValue::UInt(value))
    }

    /// Read exactly the stored Double kind; missing keys return None.
    /// Wrong kinds, invalid atom IDs and invalid UTF-8 are errors.
    pub fn atom_property_double(
        &self,
        atom: AtomId,
        key: &str,
    ) -> Result<Option<f64>, AtomPropertyReadError> {
        self.atom_property(atom, key)?
            .map(|v| Ok(v.as_double()?))
            .transpose()
    }

    /// Return a new molecule with a Double property; leave the receiver unchanged.
    #[cfg(feature = "cap-transforms")]
    pub fn with_atom_property_double(
        &self,
        atom: AtomId,
        key: &str,
        value: f64,
    ) -> Result<Self, OperationError> {
        self.with_atom_property(atom, key, &PropertyValue::Double(value))
    }

    /// Set a Double property through the existing COW operation.
    #[cfg(feature = "cap-transforms")]
    pub fn set_atom_property_double_(
        &mut self,
        atom: AtomId,
        key: &str,
        value: f64,
    ) -> Result<(), OperationError> {
        self.set_atom_property_(atom, key, &PropertyValue::Double(value))
    }

    /// Read exactly the stored Bool kind; missing keys return None.
    /// Wrong kinds, invalid atom IDs and invalid UTF-8 are errors.
    pub fn atom_property_bool(
        &self,
        atom: AtomId,
        key: &str,
    ) -> Result<Option<bool>, AtomPropertyReadError> {
        self.atom_property(atom, key)?
            .map(|v| Ok(v.as_bool()?))
            .transpose()
    }

    /// Return a new molecule with a Bool property; leave the receiver unchanged.
    #[cfg(feature = "cap-transforms")]
    pub fn with_atom_property_bool(
        &self,
        atom: AtomId,
        key: &str,
        value: bool,
    ) -> Result<Self, OperationError> {
        self.with_atom_property(atom, key, &PropertyValue::Bool(value))
    }

    /// Set a Bool property through the existing COW operation.
    #[cfg(feature = "cap-transforms")]
    pub fn set_atom_property_bool_(
        &mut self,
        atom: AtomId,
        key: &str,
        value: bool,
    ) -> Result<(), OperationError> {
        self.set_atom_property_(atom, key, &PropertyValue::Bool(value))
    }

    /// Read exactly the stored String kind; missing keys return None.
    /// Wrong kinds, invalid atom IDs and invalid UTF-8 are errors.
    pub fn atom_property_string(
        &self,
        atom: AtomId,
        key: &str,
    ) -> Result<Option<String>, AtomPropertyReadError> {
        self.atom_property(atom, key)?
            .map(|v| Ok(std::str::from_utf8(v.as_string()?.as_bytes())?.to_owned()))
            .transpose()
    }

    /// Return a new molecule with a String property; leave the receiver unchanged.
    #[cfg(feature = "cap-transforms")]
    pub fn with_atom_property_string(
        &self,
        atom: AtomId,
        key: &str,
        value: &str,
    ) -> Result<Self, OperationError> {
        self.with_atom_property(atom, key, &PropertyValue::String(value.into()))
    }

    /// Set a String property through the existing COW operation.
    #[cfg(feature = "cap-transforms")]
    pub fn set_atom_property_string_(
        &mut self,
        atom: AtomId,
        key: &str,
        value: &str,
    ) -> Result<(), OperationError> {
        self.set_atom_property_(atom, key, &PropertyValue::String(value.into()))
    }

    /// Read exactly the stored IntVector kind; missing keys return None.
    /// Wrong kinds, invalid atom IDs and invalid UTF-8 are errors.
    pub fn atom_property_int_vector(
        &self,
        atom: AtomId,
        key: &str,
    ) -> Result<Option<Vec<i32>>, AtomPropertyReadError> {
        self.atom_property(atom, key)?
            .map(|v| Ok(v.as_int_vector()?.to_vec()))
            .transpose()
    }

    /// Return a new molecule with a IntVector property; leave the receiver unchanged.
    #[cfg(feature = "cap-transforms")]
    pub fn with_atom_property_int_vector(
        &self,
        atom: AtomId,
        key: &str,
        value: &[i32],
    ) -> Result<Self, OperationError> {
        self.with_atom_property(atom, key, &PropertyValue::IntVector(value.to_vec()))
    }

    /// Set a IntVector property through the existing COW operation.
    #[cfg(feature = "cap-transforms")]
    pub fn set_atom_property_int_vector_(
        &mut self,
        atom: AtomId,
        key: &str,
        value: &[i32],
    ) -> Result<(), OperationError> {
        self.set_atom_property_(atom, key, &PropertyValue::IntVector(value.to_vec()))
    }

    /// Read exactly the stored StringVector kind; missing keys return None.
    /// Wrong kinds, invalid atom IDs and invalid UTF-8 are errors.
    pub fn atom_property_string_vector(
        &self,
        atom: AtomId,
        key: &str,
    ) -> Result<Option<Vec<String>>, AtomPropertyReadError> {
        self.atom_property(atom, key)?
            .map(|v| {
                Ok(v.as_string_vector()?
                    .iter()
                    .map(|text| std::str::from_utf8(text.as_bytes()).map(str::to_owned))
                    .collect::<Result<Vec<_>, _>>()?)
            })
            .transpose()
    }

    /// Return a new molecule with a StringVector property; leave the receiver unchanged.
    #[cfg(feature = "cap-transforms")]
    pub fn with_atom_property_string_vector(
        &self,
        atom: AtomId,
        key: &str,
        value: &[String],
    ) -> Result<Self, OperationError> {
        self.with_atom_property(
            atom,
            key,
            &PropertyValue::StringVector(value.iter().map(Into::into).collect()),
        )
    }

    /// Set a StringVector property through the existing COW operation.
    #[cfg(feature = "cap-transforms")]
    pub fn set_atom_property_string_vector_(
        &mut self,
        atom: AtomId,
        key: &str,
        value: &[String],
    ) -> Result<(), OperationError> {
        self.set_atom_property_(
            atom,
            key,
            &PropertyValue::StringVector(value.iter().map(Into::into).collect()),
        )
    }
}
