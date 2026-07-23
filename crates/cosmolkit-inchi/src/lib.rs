//! Reusable Rust source-port boundary for the IUPAC InChI engine.
//!
//! The crate deliberately exposes no generation or parsing API until the
//! pinned official InChI source call graph has been ported and tested. It is
//! independent of `cosmolkit-core`, RDKit, Python, and external executables so
//! completed functionality can be reused by other Rust projects.
//!
//! See `dev/rdkit_inchi_full_port_plan.md` in the workspace for the required
//! source inventory, marker, and parity contract.

mod source_types;

mod source;

#[cfg(test)]
mod test_support;
