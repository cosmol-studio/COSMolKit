//! Canonical projection of detached property values to source lexical bytes.
use crate::{PropertyText, PropertyValue};
pub use cosmolkit_core::PropertyStringError;

/// Return the source lexical bytes without changing the stored value or tag.
pub fn property_value_to_text(value: &PropertyValue) -> Result<PropertyText, PropertyStringError> {
    cosmolkit_core::property_value_to_string(value)
}
