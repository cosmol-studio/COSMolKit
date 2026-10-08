//! Corpus and designated source-regression preparation/comparison only.
pub mod bio_mmcif;
mod descriptors;
mod execute;
pub mod fingerprints;
pub mod mmff;
pub mod molalign;
pub mod molecular;
pub mod persistent_forcefields;
mod reference;
pub mod registry;
pub mod search;
pub mod smiles_write;
pub mod special_regression;
pub mod tautomer;
pub mod testing;
pub mod uff;
mod workflow;

pub(crate) use registry::{Corpus, Input, Record, Task};
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::path::{Path, PathBuf};
pub use workflow::{Preparation, Selection, prepare};
pub type Result<T> = std::result::Result<T, String>;

pub fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .to_owned()
}
pub fn directory() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}
pub fn expected() -> PathBuf {
    directory().join("expected")
}
pub(crate) fn read(path: &Path) -> Result<Vec<u8>> {
    std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))
}
pub(crate) fn encode<T: Serialize + ?Sized>(value: &T) -> Result<Vec<u8>> {
    serde_json::to_vec(value).map_err(|e| e.to_string())
}
pub(crate) fn digest(bytes: &[u8]) -> String {
    Sha256::digest(bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}
