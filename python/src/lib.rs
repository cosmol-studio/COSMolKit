//! Selected Python projections of the public COSMolKit facade.
//!
//! The explicit drawing selection is usable independently of the still
//! unmigrated historical default surface. It does not duplicate chemistry.

#[cfg(feature = "drawing-bindings")]
mod drawing_binding;

#[cfg(all(feature = "drawing-bindings", feature = "stubgen"))]
pyo3_stub_gen::define_stub_info_gatherer!(stub_info);

#[cfg(not(feature = "drawing-bindings"))]
include!("unmigrated_bindings.rs");
