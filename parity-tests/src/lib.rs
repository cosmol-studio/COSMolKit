//! Bounded, runtime-selected reference preparation and exact Cargo comparisons.
pub mod bio;
pub mod chemistry;
pub mod cli;
pub mod corpus;
pub mod descriptors;
pub mod fingerprints;
pub mod fuzzy_fingerprints;
pub mod geometry;
pub mod pair_fingerprints;
pub mod path_fingerprints;
pub mod random_profiles;
pub mod runner;
pub mod search_tautomer;
pub mod special;
pub mod special_geometry;
pub mod storage;
pub mod structure_tags;
pub mod tasks;
pub mod value;

pub type Result<T> = std::result::Result<T, Box<dyn std::error::Error + Send + Sync>>;
