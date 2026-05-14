//! QUARANTINED: mmcif_parser — not compiled from lib.rs until properly re-ported.
//!
//! This module previously contained summarized source comments and overlaps the
//! Gemmi-aligned BioStructure IO policy. See dev/bio_structure_io_policy.md and
//! dev/tracked_code_policy_audit.md.
//!
//! Do not expose this as a structural BioStructure parser. Structural mmCIF IO
//! belongs in the Gemmi-aligned io::bio path. RDKit-derived mmCIF parser work is
//! deferred unless a concrete Molecule compatibility need is approved, and must
//! not bypass the Gemmi-primary BioStructure parse without human approval.
#![allow(dead_code)]
pub(crate) fn quarantined_mmcif_parser() -> ! {
    panic!("mmcif_parser is quarantined until a source-backed Molecule compatibility need exists");
}
