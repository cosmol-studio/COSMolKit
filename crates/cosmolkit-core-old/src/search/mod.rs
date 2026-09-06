//! Query, SMARTS, and substructure search.

mod generic_groups;
pub mod query;
pub mod query_graph;
pub mod smarts_parse;
mod smarts_write;
pub mod substruct;
mod target;

pub use query_graph::{
    CompiledQuery, MatchResult, McsResult, QueryAtom, QueryBond, QueryGraph, QueryGraphError,
    QueryGraphOperator,
};
pub use smarts_parse::{SmartsParseParams, parse_smarts};
pub use smarts_write::{
    SmartsWriteError, get_atom_smarts, get_bond_smarts, mol_fragment_to_cx_smarts,
    mol_fragment_to_smarts, mol_to_cx_smarts, mol_to_smarts,
};
