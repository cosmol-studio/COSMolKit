//! Canonical Python projections of the public COSMolKit facade.
//!
//! Default and explicit drawing selections share one live Molecule class.
//! Historical adapters remain as source evidence, outside module compilation.

mod canonical_element_metadata;
mod canonical_fingerprint_values;
mod canonical_values;
mod drawing_binding;

#[cfg(feature = "stubgen")]
pyo3_stub_gen::define_stub_info_gatherer!(stub_info);
