// The source tree is a line-preserving Rust port of the official C engine.
// Keep its mechanical translation lints local to this module so they do not
// hide diagnostics in the public adapter or other workspace crates.
#![allow(
    dead_code,
    non_camel_case_types,
    non_snake_case,
    non_upper_case_globals,
    unused_assignments,
    unused_imports,
    unused_mut,
    unused_variables
)]

pub(crate) mod api;
pub(crate) mod base;
pub(crate) mod rdkit;
