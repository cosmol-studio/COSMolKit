//! Facade crate that re-exports COSMolKit core modules.

pub use cosmolkit_core as core;
pub use cosmolkit_core::bio;
pub use cosmolkit_core::io;
pub use cosmolkit_core::*;

/// Returns the crate version at compile time.
#[must_use]
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}
