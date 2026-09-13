//! Shared CXSMILES extension syntax for COSMolKit parsers.
//!
//! This crate owns representation-independent CX records. SMILES and SMARTS
//! parsers own the lowering of those records into their respective graph
//! states. No molecule, query graph, parser state, or operation-runtime type
//! belongs here.

mod parse;
mod records;
mod scan;

pub use parse::parse_cx_extensions;
pub use records::*;
