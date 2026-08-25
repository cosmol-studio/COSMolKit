//! Query, SMARTS, and substructure search.

mod generic_groups;
pub mod query;
pub mod smarts_parse;
mod smarts_write;
pub mod substruct;

pub use smarts_parse::{SmartsParseParams, mol_from_smarts};
pub use smarts_write::{
    SmartsWriteError, get_atom_smarts, get_bond_smarts, mol_fragment_to_cx_smarts,
    mol_fragment_to_smarts, mol_to_cx_smarts, mol_to_smarts,
};
