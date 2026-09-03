//! Core aliases for the shared concrete bond model.
//!
//! The concrete value and its local mutation methods are owned by
//! `cosmolkit-model`; core only supplies the search query payload specialization.

pub use cosmolkit_model::BondId;
pub use cosmolkit_model::{BondDirection, BondOrder, BondStereo};

pub type BondSpec = cosmolkit_model::BondSpec<crate::QueryNode<crate::BondQueryPredicate>>;
pub type Bond = cosmolkit_model::Bond<crate::QueryNode<crate::BondQueryPredicate>>;
