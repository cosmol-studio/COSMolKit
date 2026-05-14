//! QUARANTINED: pdb_parser — not compiled from lib.rs until properly re-ported.
//! See dev/bio_structure_io_policy.md, dev/tracked_code_policy_audit.md, and
//! dev/ordinary_agent_correction_guide.md.
//!
//! This module was previously tracked but contained heuristic placeholder
//! bodies that violated project policy. It is quarantined until exact
//! source reproduction under dev/source_reproduction_protocol.md is complete.
//!
//! Do not expose this as a structural BioStructure parser. Structural PDB IO
//! belongs in the Gemmi-aligned io::bio path. Any future RDKit-compatible PDB
//! molecule input must be implemented as Gemmi structural parse plus explicit
//! BioStructure-to-Molecule conversion plus RDKit-compatible post-processing,
//! not as an independent parser.
#![allow(dead_code)]
pub(crate) fn quarantined_pdb_parser() -> ! {
    panic!("pdb_parser is quarantined: heuristic code without source reproduction is forbidden");
}
