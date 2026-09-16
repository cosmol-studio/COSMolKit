//! Optional, separately compiled documentation search target in the docs package.

#[cfg(feature = "search-engine")]
mod engine;

#[cfg(all(feature = "search-engine", target_arch = "wasm32"))]
mod browser;

#[cfg(all(feature = "search-engine", target_arch = "wasm32"))]
pub use browser::mount_search;
