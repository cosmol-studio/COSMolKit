//! Shared CXSMILES extension syntax for COSMolKit parsers.
//!
//! This crate owns representation-independent CX records. SMILES and SMARTS
//! parsers own the lowering of those records into their respective graph
//! states. No molecule, query graph, parser state, or operation-runtime type
//! belongs here.

#[cfg(test)]
mod float_input;
mod parse;
mod records;
mod scan;

pub use parse::{parse_cx_extensions, parse_cx_extensions_progress};
#[doc(hidden)]
pub use parse::{
    parse_cx_extensions_progress_with_atom_window, parse_cx_extensions_with_atom_window,
};
pub use records::*;
