//! Compile-time selectors for checks owned by the live molecule runtime.
//!
//! The leaf selectors are intentionally independent. Runtime invariant and
//! operation-contract call sites use their respective leaf value, while the
//! combined value proves that the development/CI feature selected both.

pub(crate) const RUNTIME_INVARIANTS_ENABLED: bool = cfg!(feature = "runtime-invariants");
pub(crate) const OPERATION_CONTRACTS_ENABLED: bool = cfg!(feature = "op-contracts");
pub(crate) const STRICT_RUNTIME_ENABLED: bool =
    RUNTIME_INVARIANTS_ENABLED && OPERATION_CONTRACTS_ENABLED;

#[cfg(all(
    feature = "op-contracts-strict",
    not(all(feature = "runtime-invariants", feature = "op-contracts"))
))]
compile_error!("op-contracts-strict must enable both parent runtime check features");
