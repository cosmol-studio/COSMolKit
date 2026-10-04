//! Ring-count descriptor owner (RDKit `Lipinski.cpp`).

use cosmolkit_core::RingInfo;
use cosmolkit_model::{AtomId, Element};

use crate::{DescriptorError, DescriptorInput, DescriptorResult};

/// Version string the source exports next to `calcNumRings`
/// (`Lipinski.cpp:205`).
pub const NUM_RINGS_VERSION: &str = "1.0.1";

/// The ONE RingInfo-only count kernel shared by the prepared and cold
/// entrypoints.
///
/// Behavior review: reproduces `calcNumRings` exactly — the count of the
/// SUPPLIED `RingInfo` rows via the core owner
/// (`cosmolkit_core::RingInfo::num_rings`, which carries the verbatim
/// `RingInfo::numRings` anchors: initialized/length-match preconditions and
/// the atom-ring row count). The kernel NEVER recalculates a ring set; it
/// reads the rows its caller hands it. The source's unchecked
/// `rdcast<unsigned int>` is projected as a typed
/// [`DescriptorError::CountOverflow`] (same value on every realizable
/// input; fail-closed instead of wrapping).
///
/// Complexity review: one row-count read plus one integer conversion — the
/// same single-access cost class as the source's pointer fetch + size read;
/// no allocation beyond the `Result`.
pub(crate) fn num_rings_kernel(ring_info: &RingInfo) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:206-209):
    //   unsigned int calcNumRings(const ROMol &mol) {
    //     return mol.getRingInfo()->numRings();
    //   }
    // RDKit✔️✔️: unsigned int calcNumRings(const ROMol &mol) {
    // RDKit✔️✔️:   return mol.getRingInfo()->numRings();
    // RDKit✔️✔️: }
    //
    // `mol.getRingInfo()` is whatever final RingInfo the caller supplied
    // (prepared input rows, or the cold wrapper's one SSSR assignment);
    // `numRings()` is the core owner's `RingInfo::num_rings()` with its own
    // in-function anchors.
    u32::try_from(ring_info.num_rings()).map_err(|_| DescriptorError::CountOverflow {
        function: "num_rings",
        field: "ring_rows",
    })
}

/// Number of rings over prepared FINAL input, reading the SUPPLIED ring
/// rows.
///
/// Behavior review: thin delegate to [`num_rings_kernel`] over
/// `input.ring_info()` — the rows the caller supplied, whatever ring
/// finder produced them (not universally SSSR). Complexity review: the
/// kernel's single row-count read; no ring recomputation of any kind.
pub fn num_rings_prepared(input: &DescriptorInput<'_>) -> DescriptorResult<u32> {
    num_rings_kernel(input.ring_info())
}

/// Version string the source exports next to `calcNumHeterocycles`
/// (`Lipinski.cpp:232`).
pub const NUM_HETEROCYCLES_VERSION: &str = "1.0.0";

/// The ONE hetero-ring count kernel shared by the prepared and cold
/// entrypoints.
///
/// Behavior review: reproduces `calcNumHeterocycles` exactly — one pass over
/// the SUPPLIED `RingInfo::atom_rings()` rows; a row counts ONCE if ANY
/// member atom is not carbon (`element() != Element::C` is the typed exact
/// projection of the source's `getAtomicNum() != 6`, so every non-carbon
/// element — dummy atomicNum 0 and hydrogen included — counts as hetero);
/// the inner scan stops at the first non-carbon member, matching the
/// source `break`. The kernel NEVER recalculates a ring set. The final
/// count is widened fail-closed via `u32::try_from` (same value on every
/// realizable input).
///
/// Complexity review: early-exit per row — the inner scan stops at the
/// first non-carbon member, so a hetero row costs O(1) after its first
/// member check (best case O(number of rows) when the first member of
/// each row is non-carbon) versus O(total supplied ring membership)
/// worst case when all-carbon rows scan fully — exactly the source loop
/// shape; no allocation beyond the `Result`.
pub(crate) fn num_heterocycles_kernel(
    ring_info: &RingInfo,
    atoms: &[cosmolkit_model::Atom],
) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:233-242):
    //   unsigned int calcNumHeterocycles(const ROMol &mol) {
    //     unsigned int res = 0;
    //     for (const auto &iv : mol.getRingInfo()->atomRings()) {
    //       for (auto i : iv) {
    //         if (mol.getAtomWithIdx(i)->getAtomicNum() != 6) {
    //           ++res;
    //           break;
    //         }
    //       }
    //     }
    //     return res;
    //   }
    // RDKit✔️✔️: unsigned int calcNumHeterocycles(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->atomRings()) {
    // RDKit✔️✔️:     for (auto i : iv) {
    // RDKit✔️✔️:       if (mol.getAtomWithIdx(i)->getAtomicNum() != 6) {
    // RDKit✔️✔️:         ++res;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // `atomRings()` are the SUPPLIED rows; `getAtomicNum() != 6` is the
    // typed `element() != Element::C` predicate over the caller's atoms.
    let mut hetero_rows = 0usize;
    for row in ring_info.atom_rings() {
        for atom_id in row {
            if atoms[atom_id.index()].element() != Element::C {
                hetero_rows += 1;
                break;
            }
        }
    }
    u32::try_from(hetero_rows).map_err(|_| DescriptorError::CountOverflow {
        function: "num_heterocycles",
        field: "hetero_ring_rows",
    })
}

/// Number of heterocycles over prepared FINAL input, reading the SUPPLIED
/// ring rows.
///
/// Thin delegate to [`num_heterocycles_kernel`] over `input.ring_info()`
/// and `input.topology().atoms` — the rows and atoms the caller supplied
/// (not universally SSSR). Complexity: the kernel's single pass.
pub fn num_heterocycles_prepared(input: &DescriptorInput<'_>) -> DescriptorResult<u32> {
    num_heterocycles_kernel(input.ring_info(), &input.topology().atoms)
}

/// Version string the source exports next to `calcNumAromaticRings`
/// (`Lipinski.cpp:246`).
pub const NUM_AROMATIC_RINGS_VERSION: &str = "1.0.0";

/// The ONE aromatic-ring count kernel shared by the prepared and cold
/// entrypoints.
///
/// Behavior review: reproduces `calcNumAromaticRings` exactly — one pass
/// over the SUPPLIED `RingInfo::bond_rings()` rows; each row OPTIMISTICALLY
/// increments the count, then retracts once at the FIRST member bond whose
/// typed `is_aromatic()` flag is false and breaks — the net source predicate
/// counts exactly the rows whose EVERY member bond carries the aromatic
/// flag. The predicate reads the boolean flag ONLY (never `BondOrder`, ring
/// size or any chemical classification), so any partially-flagged row
/// retracts regardless of which member position is the first false one. The
/// kernel NEVER recalculates a ring set; it reads the rows its caller hands
/// it. The final count is widened fail-closed via `u32::try_from` (same
/// value on every realizable input).
///
/// Complexity review: early-exit per row — the inner scan stops at the
/// first non-aromatic member bond, so a retracting row costs O(1) after its
/// first false flag (best case O(number of rows)) versus O(total supplied
/// bond-ring membership) worst case when rows scan fully — exactly the
/// source loop shape; no allocation beyond the `Result`.
pub(crate) fn num_aromatic_rings_kernel(
    ring_info: &RingInfo,
    bonds: &[cosmolkit_model::Bond],
) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:247-257):
    //   unsigned int calcNumAromaticRings(const ROMol &mol) {
    //     unsigned int res = 0;
    //     for (const auto &iv : mol.getRingInfo()->bondRings()) {
    //       ++res;
    //       for (auto i : iv) {
    //         if (!mol.getBondWithIdx(i)->getIsAromatic()) {
    //           --res;
    //           break;
    //         }
    //       }
    //     }
    //     return res;
    //   }
    // RDKit✔️✔️: unsigned int calcNumAromaticRings(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->bondRings()) {
    // RDKit✔️✔️:     ++res;
    // RDKit✔️✔️:     for (auto i : iv) {
    // RDKit✔️✔️:       if (!mol.getBondWithIdx(i)->getIsAromatic()) {
    // RDKit✔️✔️:         --res;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // `bondRings()` are the SUPPLIED rows; `getBondWithIdx(i)` is the
    // caller's `bonds` slice indexed by `BondId`; `getIsAromatic()` is the
    // typed `Bond::is_aromatic()` flag read.
    let mut aromatic_rows = 0usize;
    for row in ring_info.bond_rings() {
        aromatic_rows += 1;
        for bond_id in row {
            if !bonds[bond_id.index()].is_aromatic() {
                aromatic_rows -= 1;
                break;
            }
        }
    }
    u32::try_from(aromatic_rows).map_err(|_| DescriptorError::CountOverflow {
        function: "num_aromatic_rings",
        field: "aromatic_ring_rows",
    })
}

/// Number of aromatic rings over prepared FINAL input, reading the SUPPLIED
/// ring rows.
///
/// Thin delegate to [`num_aromatic_rings_kernel`] over `input.ring_info()`
/// and `input.topology().bonds` — the rows and bonds the caller supplied
/// (not universally SSSR). Complexity: the kernel's single pass.
pub fn num_aromatic_rings_prepared(input: &DescriptorInput<'_>) -> DescriptorResult<u32> {
    num_aromatic_rings_kernel(input.ring_info(), &input.topology().bonds)
}

/// Version string the source exports next to `calcNumSaturatedRings`
/// (`Lipinski.cpp:260`).
pub const NUM_SATURATED_RINGS_VERSION: &str = "1.0.0";

/// The ONE saturated-ring count kernel shared by the prepared and cold
/// entrypoints.
///
/// Behavior review: reproduces `calcNumSaturatedRings` exactly — one pass
/// over the SUPPLIED `RingInfo::bond_rings()` rows; each row OPTIMISTICALLY
/// increments the count, then retracts once at the FIRST member bond
/// failing the source conjunct `getBondType() != Bond::SINGLE ||
/// getIsAromatic()` and breaks — the net predicate counts exactly the rows
/// whose EVERY member bond has order Single AND aromatic flag false. The
/// two arms are independent: a Single-ORDERED bond that still carries the
/// aromatic flag retracts via the flag arm exactly as in the source (and a
/// Double/Triple/Aromatic-ordered bond retracts via the order arm). The
/// kernel NEVER recalculates a ring set. The final count is widened
/// fail-closed via `u32::try_from` (same value on every realizable input).
///
/// Complexity review: early-exit per row — the inner scan stops at the
/// first failing member bond, so a retracting row costs O(1) after its
/// first failure (best case O(number of rows)) versus O(total supplied
/// bond-ring membership) worst case when rows scan fully — exactly the
/// source loop shape; no allocation beyond the `Result`.
pub(crate) fn num_saturated_rings_kernel(
    ring_info: &RingInfo,
    bonds: &[cosmolkit_model::Bond],
) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:261-272):
    //   unsigned int calcNumSaturatedRings(const ROMol &mol) {
    //     unsigned int res = 0;
    //     for (const auto &iv : mol.getRingInfo()->bondRings()) {
    //       ++res;
    //       for (int i : iv) {
    //         if (mol.getBondWithIdx(i)->getBondType() != Bond::SINGLE ||
    //             mol.getBondWithIdx(i)->getIsAromatic()) {
    //           --res;
    //           break;
    //         }
    //       }
    //     }
    //     return res;
    //   }
    // RDKit✔️✔️: unsigned int calcNumSaturatedRings(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->bondRings()) {
    // RDKit✔️✔️:     ++res;
    // RDKit✔️✔️:     for (int i : iv) {
    // RDKit✔️✔️:       if (mol.getBondWithIdx(i)->getBondType() != Bond::SINGLE ||
    // RDKit✔️✔️:           mol.getBondWithIdx(i)->getIsAromatic()) {
    // RDKit✔️✔️:         --res;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // `bondRings()` are the SUPPLIED rows; `getBondWithIdx(i)` is the
    // caller's `bonds` slice indexed by `BondId`; `getBondType()` is the
    // typed `Bond::order()` read and `getIsAromatic()` the typed
    // `Bond::is_aromatic()` flag read (both carry their own anchors in the
    // model crate).
    let mut saturated_rows = 0usize;
    for row in ring_info.bond_rings() {
        saturated_rows += 1;
        for bond_id in row {
            let bond = &bonds[bond_id.index()];
            if bond.order() != cosmolkit_model::BondOrder::Single || bond.is_aromatic() {
                saturated_rows -= 1;
                break;
            }
        }
    }
    u32::try_from(saturated_rows).map_err(|_| DescriptorError::CountOverflow {
        function: "num_saturated_rings",
        field: "saturated_ring_rows",
    })
}

/// Number of saturated rings over prepared FINAL input, reading the
/// SUPPLIED ring rows.
///
/// Thin delegate to [`num_saturated_rings_kernel`] over `input.ring_info()`
/// and `input.topology().bonds` — the rows and bonds the caller supplied
/// (not universally SSSR). Complexity: the kernel's single pass.
pub fn num_saturated_rings_prepared(input: &DescriptorInput<'_>) -> DescriptorResult<u32> {
    num_saturated_rings_kernel(input.ring_info(), &input.topology().bonds)
}

/// Version string the source exports next to `calcNumAliphaticRings`
/// (`Lipinski.cpp:275`).
pub const NUM_ALIPHATIC_RINGS_VERSION: &str = "1.0.0";

/// The ONE aliphatic-ring count kernel shared by the prepared and cold
/// entrypoints.
///
/// Behavior review: reproduces `calcNumAliphaticRings` exactly — one pass
/// over the SUPPLIED `RingInfo::bond_rings()` rows with NO optimistic
/// increment: the ++ sits INSIDE the inner loop at the FIRST member bond
/// whose typed `is_aromatic()` flag is false, followed by `break`. Net
/// predicate: a row counts ONCE iff it has AT LEAST ONE non-aromatic
/// member bond — the set-complement of the N03 all-aromatic predicate over
/// the same supplied rows (a fully aromatic row scans to its end and adds
/// nothing). The predicate reads the boolean flag ONLY, never `BondOrder`.
/// The kernel NEVER recalculates a ring set. The final count is widened
/// fail-closed via `u32::try_from` (same value on every realizable input).
///
/// Complexity review: early-exit per row — the inner scan stops at the
/// first non-aromatic member bond, so a counting row costs O(1) after its
/// first member check (best case O(number of rows)) versus O(total
/// supplied bond-ring membership) worst case when fully-aromatic rows scan
/// to their end — exactly the source loop shape; no allocation beyond the
/// `Result`.
pub(crate) fn num_aliphatic_rings_kernel(
    ring_info: &RingInfo,
    bonds: &[cosmolkit_model::Bond],
) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:276-286):
    //   unsigned int calcNumAliphaticRings(const ROMol &mol) {
    //     unsigned int res = 0;
    //     for (const auto &iv : mol.getRingInfo()->bondRings()) {
    //       for (auto i : iv) {
    //         if (!mol.getBondWithIdx(i)->getIsAromatic()) {
    //           ++res;
    //           break;
    //         }
    //       }
    //     }
    //     return res;
    //   }
    // RDKit✔️✔️: unsigned int calcNumAliphaticRings(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->bondRings()) {
    // RDKit✔️✔️:     for (auto i : iv) {
    // RDKit✔️✔️:       if (!mol.getBondWithIdx(i)->getIsAromatic()) {
    // RDKit✔️✔️:         ++res;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // `bondRings()` are the SUPPLIED rows; `getBondWithIdx(i)` is the
    // caller's `bonds` slice indexed by `BondId`; `getIsAromatic()` is the
    // typed `Bond::is_aromatic()` flag read.
    let mut aliphatic_rows = 0usize;
    for row in ring_info.bond_rings() {
        for bond_id in row {
            if !bonds[bond_id.index()].is_aromatic() {
                aliphatic_rows += 1;
                break;
            }
        }
    }
    u32::try_from(aliphatic_rows).map_err(|_| DescriptorError::CountOverflow {
        function: "num_aliphatic_rings",
        field: "aliphatic_ring_rows",
    })
}

/// Number of aliphatic rings (ring rows with AT LEAST ONE non-aromatic
/// member bond) over prepared FINAL input, reading the SUPPLIED ring rows.
///
/// Thin delegate to [`num_aliphatic_rings_kernel`] over `input.ring_info()`
/// and `input.topology().bonds` — the rows and bonds the caller supplied
/// (not universally SSSR). Complexity: the kernel's single pass.
pub fn num_aliphatic_rings_prepared(input: &DescriptorInput<'_>) -> DescriptorResult<u32> {
    num_aliphatic_rings_kernel(input.ring_info(), &input.topology().bonds)
}

/// Version string the source exports next to `calcNumAromaticHeterocycles`
/// (`Lipinski.cpp:288`).
pub const NUM_AROMATIC_HETEROCYCLES_VERSION: &str = "1.0.0";

/// The ONE aromatic-heterocycle count kernel shared by the prepared and
/// cold entrypoints.
///
/// Behavior review: reproduces `calcNumAromaticHeterocycles` exactly — per
/// SUPPLIED `RingInfo::bond_rings()` row a `count_it` flag starts false;
/// the inner scan forces `count_it = false` and BREAKS at the FIRST
/// non-aromatic member bond; otherwise, while still unset, a member bond
/// with a non-carbon endpoint (begin OR end `element() != Element::C`, the
/// typed `getAtomicNum() != 6` projection) sets it true (sticky — the
/// source's `!countIt` guard is a perf nicety, not a semantic difference;
/// the source itself notes each hetero atom is checked twice via its two
/// bonds, "kind of doofy"). AFTER the inner loop, `count_it` increments
/// the count. Net predicate: a row counts iff EVERY member bond is
/// aromatic AND at least one member bond has a non-carbon endpoint. The
/// kernel NEVER recalculates a ring set. Final count widened fail-closed
/// via `u32::try_from` (same value on every realizable input).
///
/// Complexity review: early break per row at the first non-aromatic
/// member; the sticky flag skips endpoint checks after the first hetero
/// bond — exactly the source loop shape; no allocation beyond the
/// `Result`.
pub(crate) fn num_aromatic_heterocycles_kernel(
    ring_info: &RingInfo,
    bonds: &[cosmolkit_model::Bond],
    atoms: &[cosmolkit_model::Atom],
) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:289-311):
    //   unsigned int calcNumAromaticHeterocycles(const ROMol &mol) {
    //     unsigned int res = 0;
    //     for (const auto &iv : mol.getRingInfo()->bondRings()) {
    //       bool countIt = false;
    //       for (auto i : iv) {
    //         if (!mol.getBondWithIdx(i)->getIsAromatic()) {
    //           countIt = false;
    //           break;
    //         }
    //         // we're checking each atom twice, which is kind of doofy, but this
    //         // function is hopefully not going to be a big time sink.
    //         if (!countIt &&
    //             (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
    //              mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6)) {
    //           countIt = true;
    //         }
    //       }
    //       if (countIt) {
    //         ++res;
    //       }
    //     }
    //     return res;
    //   }
    // RDKit✔️✔️: unsigned int calcNumAromaticHeterocycles(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->bondRings()) {
    // RDKit✔️✔️:     bool countIt = false;
    // RDKit✔️✔️:     for (auto i : iv) {
    // RDKit✔️✔️:       if (!mol.getBondWithIdx(i)->getIsAromatic()) {
    // RDKit✔️✔️:         countIt = false;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // we're checking each atom twice, which is kind of doofy, but this
    // RDKit✔️✔️:       // function is hopefully not going to be a big time sink.
    // RDKit✔️✔️:       if (!countIt &&
    // RDKit✔️✔️:           (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
    // RDKit✔️✔️:            mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6)) {
    // RDKit✔️✔️:         countIt = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (countIt) {
    // RDKit✔️✔️:       ++res;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // `bondRings()` are the SUPPLIED rows; `getBondWithIdx(i)` is the
    // caller's `bonds` slice indexed by `BondId`; `getBeginAtom()` /
    // `getEndAtom()` are the caller's `atoms` slice indexed by the bond's
    // typed `begin()` / `end()` AtomIds.
    let mut aromatic_hetero_rows = 0usize;
    for row in ring_info.bond_rings() {
        let mut count_it = false;
        'members: for bond_id in row {
            let bond = &bonds[bond_id.index()];
            if !bond.is_aromatic() {
                count_it = false;
                break 'members;
            }
            if !count_it
                && (atoms[bond.begin().index()].element() != Element::C
                    || atoms[bond.end().index()].element() != Element::C)
            {
                count_it = true;
            }
        }
        if count_it {
            aromatic_hetero_rows += 1;
        }
    }
    u32::try_from(aromatic_hetero_rows).map_err(|_| DescriptorError::CountOverflow {
        function: "num_aromatic_heterocycles",
        field: "aromatic_heterocycle_rows",
    })
}

/// Number of aromatic heterocycles over prepared FINAL input, reading the
/// SUPPLIED ring rows.
///
/// Thin delegate to [`num_aromatic_heterocycles_kernel`] over
/// `input.ring_info()`, `input.topology().bonds` and
/// `input.topology().atoms` — the rows the caller supplied (not
/// universally SSSR). Complexity: the kernel's single pass.
pub fn num_aromatic_heterocycles_prepared(input: &DescriptorInput<'_>) -> DescriptorResult<u32> {
    num_aromatic_heterocycles_kernel(
        input.ring_info(),
        &input.topology().bonds,
        &input.topology().atoms,
    )
}

/// Version string the source exports next to `calcNumAromaticCarbocycles`
/// (`Lipinski.cpp:312`).
pub const NUM_AROMATIC_CARBOCYCLES_VERSION: &str = "1.0.0";

/// The ONE aromatic-carbocycle count kernel shared by the prepared and
/// cold entrypoints.
///
/// Behavior review: reproduces `calcNumAromaticCarbocycles` exactly — per
/// SUPPLIED `RingInfo::bond_rings()` row `count_it` starts TRUE; the inner
/// scan breaks with `count_it = false` at EITHER the first non-aromatic
/// member bond OR the first member bond with a non-carbon endpoint (begin
/// OR end `element() != Element::C`, the typed `getAtomicNum() != 6`
/// projection — the source again notes each atom is checked twice, here
/// "kind of doofy ... big time sync" wording verbatim); after the loop
/// `count_it` increments the count. Net predicate: counts fully-aromatic
/// ALL-CARBON rows (the carbocycle complement of the N06 heterocycle
/// predicate over the same rows). The kernel NEVER recalculates a ring
/// set. Final count widened fail-closed via `u32::try_from` (same value
/// on every realizable input).
///
/// Complexity review: early break per row at the first failing member
/// (either arm) — exactly the source loop shape; no allocation beyond
/// the `Result`.
pub(crate) fn num_aromatic_carbocycles_kernel(
    ring_info: &RingInfo,
    bonds: &[cosmolkit_model::Bond],
    atoms: &[cosmolkit_model::Atom],
) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:313-333):
    //   unsigned int calcNumAromaticCarbocycles(const ROMol &mol) {
    //     unsigned int res = 0;
    //     for (const auto &iv : mol.getRingInfo()->bondRings()) {
    //       bool countIt = true;
    //       for (auto i : iv) {
    //         if (!mol.getBondWithIdx(i)->getIsAromatic()) {
    //           countIt = false;
    //           break;
    //         }
    //         // we're checking each atom twice, which is kind of doofy, but this
    //         // function is hopefully not going to be a big time sync.
    //         if (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
    //             mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6) {
    //           countIt = false;
    //           break;
    //         }
    //       }
    //       if (countIt) {
    //         ++res;
    //       }
    //     }
    //     return res;
    //   }
    // RDKit✔️✔️: unsigned int calcNumAromaticCarbocycles(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->bondRings()) {
    // RDKit✔️✔️:     bool countIt = true;
    // RDKit✔️✔️:     for (auto i : iv) {
    // RDKit✔️✔️:       if (!mol.getBondWithIdx(i)->getIsAromatic()) {
    // RDKit✔️✔️:         countIt = false;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // we're checking each atom twice, which is kind of doofy, but this
    // RDKit✔️✔️:       // function is hopefully not going to be a big time sync.
    // RDKit✔️✔️:       if (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
    // RDKit✔️✔️:           mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6) {
    // RDKit✔️✔️:         countIt = false;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (countIt) {
    // RDKit✔️✔️:       ++res;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // `bondRings()` are the SUPPLIED rows; `getBondWithIdx(i)` is the
    // caller's `bonds` slice indexed by `BondId`; `getBeginAtom()` /
    // `getEndAtom()` are the caller's `atoms` slice indexed by the bond's
    // typed `begin()` / `end()` AtomIds.
    let mut aromatic_carbo_rows = 0usize;
    for row in ring_info.bond_rings() {
        let mut count_it = true;
        'members: for bond_id in row {
            let bond = &bonds[bond_id.index()];
            if !bond.is_aromatic() {
                count_it = false;
                break 'members;
            }
            if atoms[bond.begin().index()].element() != Element::C
                || atoms[bond.end().index()].element() != Element::C
            {
                count_it = false;
                break 'members;
            }
        }
        if count_it {
            aromatic_carbo_rows += 1;
        }
    }
    u32::try_from(aromatic_carbo_rows).map_err(|_| DescriptorError::CountOverflow {
        function: "num_aromatic_carbocycles",
        field: "aromatic_carbocycle_rows",
    })
}

/// Number of aromatic carbocycles over prepared FINAL input, reading the
/// SUPPLIED ring rows.
///
/// Thin delegate to [`num_aromatic_carbocycles_kernel`] over
/// `input.ring_info()`, `input.topology().bonds` and
/// `input.topology().atoms` — the rows the caller supplied (not
/// universally SSSR). Complexity: the kernel's single pass.
pub fn num_aromatic_carbocycles_prepared(input: &DescriptorInput<'_>) -> DescriptorResult<u32> {
    num_aromatic_carbocycles_kernel(
        input.ring_info(),
        &input.topology().bonds,
        &input.topology().atoms,
    )
}

/// Version string the source exports next to
/// `calcNumAliphaticHeterocycles` (`Lipinski.cpp:334`).
pub const NUM_ALIPHATIC_HETEROCYCLES_VERSION: &str = "1.0.0";

/// The ONE aliphatic-heterocycle count kernel shared by the prepared and
/// cold entrypoints.
///
/// Behavior review: reproduces `calcNumAliphaticHeterocycles` exactly —
/// per SUPPLIED `RingInfo::bond_rings()` row TWO independent flags start
/// false: `has_aliph` is set by ANY non-aromatic member bond, and
/// `has_hetero` (sticky, guarded `!hasHetero`) is set by any member bond
/// with a non-carbon endpoint (begin OR end `element() != Element::C`,
/// the typed `getAtomicNum() != 6` projection; the source's "checking
/// each atom twice, kind of doofy ... big time sink" comment preserved
/// verbatim). The inner scan runs over ALL members — NO break anywhere
/// (unlike the N03-N07 kernels). After the loop, BOTH flags => ++. Net
/// predicate: a row counts iff it contains at least one non-aromatic
/// member bond AND at least one hetero-endpoint member bond, at any
/// positions. The kernel NEVER recalculates a ring set. Final count
/// widened fail-closed via `u32::try_from` (same value on every
/// realizable input).
///
/// Complexity review: full-row scan with two sticky flags — O(total
/// supplied bond-ring membership), exactly the source loop shape (no
/// early exit in the source either); no allocation beyond the `Result`.
pub(crate) fn num_aliphatic_heterocycles_kernel(
    ring_info: &RingInfo,
    bonds: &[cosmolkit_model::Bond],
    atoms: &[cosmolkit_model::Atom],
) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:335-358):
    //   unsigned int calcNumAliphaticHeterocycles(const ROMol &mol) {
    //     unsigned int res = 0;
    //     for (const auto &iv : mol.getRingInfo()->bondRings()) {
    //       bool hasAliph = false;
    //       bool hasHetero = false;
    //       for (auto i : iv) {
    //         if (!mol.getBondWithIdx(i)->getIsAromatic()) {
    //           hasAliph = true;
    //         }
    //         // we're checking each atom twice, which is kind of doofy, but this
    //         // function is hopefully not going to be a big time sink.
    //         if (!hasHetero &&
    //             (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
    //              mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6)) {
    //           hasHetero = true;
    //         }
    //       }
    //       if (hasHetero && hasAliph) {
    //         ++res;
    //       }
    //     }
    //     return res;
    //   }
    // RDKit✔️✔️: unsigned int calcNumAliphaticHeterocycles(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->bondRings()) {
    // RDKit✔️✔️:     bool hasAliph = false;
    // RDKit✔️✔️:     bool hasHetero = false;
    // RDKit✔️✔️:     for (auto i : iv) {
    // RDKit✔️✔️:       if (!mol.getBondWithIdx(i)->getIsAromatic()) {
    // RDKit✔️✔️:         hasAliph = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // we're checking each atom twice, which is kind of doofy, but this
    // RDKit✔️✔️:       // function is hopefully not going to be a big time sink.
    // RDKit✔️✔️:       if (!hasHetero &&
    // RDKit✔️✔️:           (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
    // RDKit✔️✔️:            mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6)) {
    // RDKit✔️✔️:         hasHetero = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (hasHetero && hasAliph) {
    // RDKit✔️✔️:       ++res;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // `bondRings()` are the SUPPLIED rows; `getBondWithIdx(i)` is the
    // caller's `bonds` slice indexed by `BondId`; `getBeginAtom()` /
    // `getEndAtom()` are the caller's `atoms` slice indexed by the bond's
    // typed `begin()` / `end()` AtomIds.
    let mut aliphatic_hetero_rows = 0usize;
    for row in ring_info.bond_rings() {
        let mut has_aliph = false;
        let mut has_hetero = false;
        for bond_id in row {
            let bond = &bonds[bond_id.index()];
            if !bond.is_aromatic() {
                has_aliph = true;
            }
            if !has_hetero
                && (atoms[bond.begin().index()].element() != Element::C
                    || atoms[bond.end().index()].element() != Element::C)
            {
                has_hetero = true;
            }
        }
        if has_hetero && has_aliph {
            aliphatic_hetero_rows += 1;
        }
    }
    u32::try_from(aliphatic_hetero_rows).map_err(|_| DescriptorError::CountOverflow {
        function: "num_aliphatic_heterocycles",
        field: "aliphatic_heterocycle_rows",
    })
}

/// Number of aliphatic heterocycles over prepared FINAL input, reading
/// the SUPPLIED ring rows.
///
/// Thin delegate to [`num_aliphatic_heterocycles_kernel`] over
/// `input.ring_info()`, `input.topology().bonds` and
/// `input.topology().atoms` — the rows the caller supplied (not
/// universally SSSR). Complexity: the kernel's full-row pass.
pub fn num_aliphatic_heterocycles_prepared(input: &DescriptorInput<'_>) -> DescriptorResult<u32> {
    num_aliphatic_heterocycles_kernel(
        input.ring_info(),
        &input.topology().bonds,
        &input.topology().atoms,
    )
}

/// Version string the source exports next to
/// `calcNumAliphaticCarbocycles` (`Lipinski.cpp:360`).
pub const NUM_ALIPHATIC_CARBOCYCLES_VERSION: &str = "1.0.0";

/// The ONE aliphatic-carbocycle count kernel shared by the prepared and
/// cold entrypoints.
///
/// Behavior review: reproduces `calcNumAliphaticCarbocycles` exactly —
/// per SUPPLIED `RingInfo::bond_rings()` row two flags start false: a
/// non-aromatic member bond sets `has_aliph` WITHOUT breaking (the scan
/// continues); a member bond with a non-carbon endpoint (begin OR end
/// `element() != Element::C`, the typed `getAtomicNum() != 6` projection)
/// sets `has_hetero` AND BREAKS immediately — unlike the N08 kernel there
/// is NO sticky guard, the source sets-and-breaks directly (its inline
/// comment reads "big time sync", preserved verbatim). The break is an
/// outcome-preserving early exit: once `has_hetero` is true the row is
/// already disqualified, so skipping later member bonds cannot change the
/// result. After the loop `has_aliph && !has_hetero` increments the count.
/// Net predicate: a row counts iff it has at least one non-aromatic member
/// bond AND no hetero-endpoint member bond at all. The kernel NEVER
/// recalculates a ring set. Final count widened fail-closed via
/// `u32::try_from` (same value on every realizable input).
///
/// Complexity review: early break at the first hetero-endpoint member;
/// otherwise full-row scan for the all-carbon rows — exactly the source
/// loop shape; no allocation beyond the `Result`.
pub(crate) fn num_aliphatic_carbocycles_kernel(
    ring_info: &RingInfo,
    bonds: &[cosmolkit_model::Bond],
    atoms: &[cosmolkit_model::Atom],
) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:361-382):
    //   unsigned int calcNumAliphaticCarbocycles(const ROMol &mol) {
    //     unsigned int res = 0;
    //     for (const auto &iv : mol.getRingInfo()->bondRings()) {
    //       bool hasAliph = false;
    //       bool hasHetero = false;
    //       for (auto i : iv) {
    //         if (!mol.getBondWithIdx(i)->getIsAromatic()) {
    //           hasAliph = true;
    //         }
    //         // we're checking each atom twice, which is kind of doofy, but this
    //         // function is hopefully not going to be a big time sync.
    //         if (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
    //             mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6) {
    //           hasHetero = true;
    //           break;
    //         }
    //       }
    //       if (hasAliph && !hasHetero) {
    //         ++res;
    //       }
    //     }
    //     return res;
    //   }
    // RDKit✔️✔️: unsigned int calcNumAliphaticCarbocycles(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->bondRings()) {
    // RDKit✔️✔️:     bool hasAliph = false;
    // RDKit✔️✔️:     bool hasHetero = false;
    // RDKit✔️✔️:     for (auto i : iv) {
    // RDKit✔️✔️:       if (!mol.getBondWithIdx(i)->getIsAromatic()) {
    // RDKit✔️✔️:         hasAliph = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // we're checking each atom twice, which is kind of doofy, but this
    // RDKit✔️✔️:       // function is hopefully not going to be a big time sync.
    // RDKit✔️✔️:       if (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
    // RDKit✔️✔️:           mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6) {
    // RDKit✔️✔️:         hasHetero = true;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (hasAliph && !hasHetero) {
    // RDKit✔️✔️:       ++res;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // `bondRings()` are the SUPPLIED rows; `getBondWithIdx(i)` is the
    // caller's `bonds` slice indexed by `BondId`; `getBeginAtom()` /
    // `getEndAtom()` are the caller's `atoms` slice indexed by the bond's
    // typed `begin()` / `end()` AtomIds.
    let mut aliphatic_carbo_rows = 0usize;
    for row in ring_info.bond_rings() {
        let mut has_aliph = false;
        let mut has_hetero = false;
        'members: for bond_id in row {
            let bond = &bonds[bond_id.index()];
            if !bond.is_aromatic() {
                has_aliph = true;
            }
            if atoms[bond.begin().index()].element() != Element::C
                || atoms[bond.end().index()].element() != Element::C
            {
                has_hetero = true;
                break 'members;
            }
        }
        if has_aliph && !has_hetero {
            aliphatic_carbo_rows += 1;
        }
    }
    u32::try_from(aliphatic_carbo_rows).map_err(|_| DescriptorError::CountOverflow {
        function: "num_aliphatic_carbocycles",
        field: "aliphatic_carbocycle_rows",
    })
}

/// Number of aliphatic carbocycles over prepared FINAL input, reading the
/// SUPPLIED ring rows.
///
/// Thin delegate to [`num_aliphatic_carbocycles_kernel`] over
/// `input.ring_info()`, `input.topology().bonds` and
/// `input.topology().atoms` — the rows the caller supplied (not
/// universally SSSR). Complexity: the kernel's pass with hetero early
/// break.
pub fn num_aliphatic_carbocycles_prepared(input: &DescriptorInput<'_>) -> DescriptorResult<u32> {
    num_aliphatic_carbocycles_kernel(
        input.ring_info(),
        &input.topology().bonds,
        &input.topology().atoms,
    )
}

/// Version string the source exports next to
/// `calcNumSaturatedHeterocycles` (`Lipinski.cpp:384`).
pub const NUM_SATURATED_HETEROCYCLES_VERSION: &str = "1.0.0";

/// The ONE saturated-heterocycle count kernel shared by the prepared and
/// cold entrypoints.
///
/// Behavior review: reproduces `calcNumSaturatedHeterocycles` exactly —
/// per SUPPLIED `RingInfo::bond_rings()` row `count_it` starts false; the
/// inner scan breaks with `count_it = false` at the FIRST member bond
/// failing the saturated conjunct `getBondType() != Bond::SINGLE ||
/// getIsAromatic()` (identical arms to the N04 kernel — both
/// independent); otherwise the N06-style sticky hetero arm
/// (`!countIt &&` begin-or-end endpoint `element() != Element::C`, the
/// typed `getAtomicNum() != 6` projection, sets `count_it = true`). After
/// the loop `count_it` increments the count. Net predicate: a row counts
/// iff EVERY member bond is Single-ordered AND non-aromatic AND at least
/// one member bond has a hetero endpoint — the hetero arm only ever fires
/// among all-Single non-aromatic bonds since any other bond breaks first
/// (source inline comment "big time sync", preserved verbatim). The
/// kernel NEVER recalculates a ring set. Final count widened fail-closed
/// via `u32::try_from` (same value on every realizable input).
///
/// Complexity review: early break at the first non-saturated member; the
/// sticky flag skips endpoint checks after the first hetero bond —
/// exactly the source loop shape; no allocation beyond the `Result`.
pub(crate) fn num_saturated_heterocycles_kernel(
    ring_info: &RingInfo,
    bonds: &[cosmolkit_model::Bond],
    atoms: &[cosmolkit_model::Atom],
) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:385-407):
    //   unsigned int calcNumSaturatedHeterocycles(const ROMol &mol) {
    //     unsigned int res = 0;
    //     for (const auto &iv : mol.getRingInfo()->bondRings()) {
    //       bool countIt = false;
    //       for (auto i : iv) {
    //         if (mol.getBondWithIdx(i)->getBondType() != Bond::SINGLE ||
    //             mol.getBondWithIdx(i)->getIsAromatic()) {
    //           countIt = false;
    //           break;
    //         }
    //         // we're checking each atom twice, which is kind of doofy, but this
    //         // function is hopefully not going to be a big time sync.
    //         if (!countIt &&
    //             (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
    //              mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6)) {
    //           countIt = true;
    //         }
    //       }
    //       if (countIt) {
    //         ++res;
    //       }
    //     }
    //     return res;
    //   }
    // RDKit✔️✔️: unsigned int calcNumSaturatedHeterocycles(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->bondRings()) {
    // RDKit✔️✔️:     bool countIt = false;
    // RDKit✔️✔️:     for (auto i : iv) {
    // RDKit✔️✔️:       if (mol.getBondWithIdx(i)->getBondType() != Bond::SINGLE ||
    // RDKit✔️✔️:           mol.getBondWithIdx(i)->getIsAromatic()) {
    // RDKit✔️✔️:         countIt = false;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // we're checking each atom twice, which is kind of doofy, but this
    // RDKit✔️✔️:       // function is hopefully not going to be a big time sync.
    // RDKit✔️✔️:       if (!countIt &&
    // RDKit✔️✔️:           (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
    // RDKit✔️✔️:            mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6)) {
    // RDKit✔️✔️:         countIt = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (countIt) {
    // RDKit✔️✔️:       ++res;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // `bondRings()` are the SUPPLIED rows; `getBondWithIdx(i)` is the
    // caller's `bonds` slice indexed by `BondId`; `getBeginAtom()` /
    // `getEndAtom()` are the caller's `atoms` slice indexed by the bond's
    // typed `begin()` / `end()` AtomIds.
    let mut saturated_hetero_rows = 0usize;
    for row in ring_info.bond_rings() {
        let mut count_it = false;
        'members: for bond_id in row {
            let bond = &bonds[bond_id.index()];
            if bond.order() != cosmolkit_model::BondOrder::Single || bond.is_aromatic() {
                count_it = false;
                break 'members;
            }
            if !count_it
                && (atoms[bond.begin().index()].element() != Element::C
                    || atoms[bond.end().index()].element() != Element::C)
            {
                count_it = true;
            }
        }
        if count_it {
            saturated_hetero_rows += 1;
        }
    }
    u32::try_from(saturated_hetero_rows).map_err(|_| DescriptorError::CountOverflow {
        function: "num_saturated_heterocycles",
        field: "saturated_heterocycle_rows",
    })
}

/// Number of saturated heterocycles over prepared FINAL input, reading the
/// SUPPLIED ring rows.
///
/// Thin delegate to [`num_saturated_heterocycles_kernel`] over
/// `input.ring_info()`, `input.topology().bonds` and
/// `input.topology().atoms` — the rows the caller supplied (not
/// universally SSSR). Complexity: the kernel's pass with early break.
pub fn num_saturated_heterocycles_prepared(input: &DescriptorInput<'_>) -> DescriptorResult<u32> {
    num_saturated_heterocycles_kernel(
        input.ring_info(),
        &input.topology().bonds,
        &input.topology().atoms,
    )
}

/// Version string the source exports next to
/// `calcNumSaturatedCarbocycles` (`Lipinski.cpp:409`).
pub const NUM_SATURATED_CARBOCYCLES_VERSION: &str = "1.0.0";

/// The ONE saturated-carbocycle count kernel shared by the prepared and
/// cold entrypoints.
///
/// Behavior review: reproduces `calcNumSaturatedCarbocycles` exactly —
/// per SUPPLIED `RingInfo::bond_rings()` row `count_it` starts TRUE; the
/// inner scan breaks with `count_it = false` at EITHER the first member
/// bond failing the saturated conjunct `getBondType() != Bond::SINGLE ||
/// getIsAromatic()` (identical independent arms to the N04 kernel) OR the
/// first member bond with a non-carbon endpoint (begin OR end
/// `element() != Element::C`, the typed `getAtomicNum() != 6` projection;
/// direct set+break, NO sticky guard — the N07-style shape; "time sync"
/// inline comment preserved verbatim). After the loop `count_it`
/// increments the count. Net predicate: a row counts iff EVERY member
/// bond is Single-ordered, non-aromatic AND has two carbon endpoints —
/// saturated CARBOcycles, DISTINCT from the N04 count which ignores
/// endpoints (piperidine: N04 counts 1, this kernel counts 0). The
/// kernel NEVER recalculates a ring set. Final count widened fail-closed
/// via `u32::try_from` (same value on every realizable input).
///
/// Complexity review: early break at the first failing member on either
/// arm — exactly the source loop shape; no allocation beyond the `Result`.
pub(crate) fn num_saturated_carbocycles_kernel(
    ring_info: &RingInfo,
    bonds: &[cosmolkit_model::Bond],
    atoms: &[cosmolkit_model::Atom],
) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:410-432):
    //   unsigned int calcNumSaturatedCarbocycles(const ROMol &mol) {
    //     unsigned int res = 0;
    //     for (const auto &iv : mol.getRingInfo()->bondRings()) {
    //       bool countIt = true;
    //       for (auto i : iv) {
    //         if (mol.getBondWithIdx(i)->getBondType() != Bond::SINGLE ||
    //             mol.getBondWithIdx(i)->getIsAromatic()) {
    //           countIt = false;
    //           break;
    //         }
    //         // we're checking each atom twice, which is kind of doofy, but this
    //         // function is hopefully not going to be a big time sync.
    //         if (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
    //             mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6) {
    //           countIt = false;
    //           break;
    //         }
    //       }
    //       if (countIt) {
    //         ++res;
    //       }
    //     }
    //     return res;
    //   }
    // RDKit✔️✔️: unsigned int calcNumSaturatedCarbocycles(const ROMol &mol) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   for (const auto &iv : mol.getRingInfo()->bondRings()) {
    // RDKit✔️✔️:     bool countIt = true;
    // RDKit✔️✔️:     for (auto i : iv) {
    // RDKit✔️✔️:       if (mol.getBondWithIdx(i)->getBondType() != Bond::SINGLE ||
    // RDKit✔️✔️:           mol.getBondWithIdx(i)->getIsAromatic()) {
    // RDKit✔️✔️:         countIt = false;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // we're checking each atom twice, which is kind of doofy, but this
    // RDKit✔️✔️:       // function is hopefully not going to be a big time sync.
    // RDKit✔️✔️:       if (mol.getBondWithIdx(i)->getBeginAtom()->getAtomicNum() != 6 ||
    // RDKit✔️✔️:           mol.getBondWithIdx(i)->getEndAtom()->getAtomicNum() != 6) {
    // RDKit✔️✔️:         countIt = false;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (countIt) {
    // RDKit✔️✔️:       ++res;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // `bondRings()` are the SUPPLIED rows; `getBondWithIdx(i)` is the
    // caller's `bonds` slice indexed by `BondId`; `getBeginAtom()` /
    // `getEndAtom()` are the caller's `atoms` slice indexed by the bond's
    // typed `begin()` / `end()` AtomIds.
    let mut saturated_carbo_rows = 0usize;
    for row in ring_info.bond_rings() {
        let mut count_it = true;
        'members: for bond_id in row {
            let bond = &bonds[bond_id.index()];
            if bond.order() != cosmolkit_model::BondOrder::Single || bond.is_aromatic() {
                count_it = false;
                break 'members;
            }
            if atoms[bond.begin().index()].element() != Element::C
                || atoms[bond.end().index()].element() != Element::C
            {
                count_it = false;
                break 'members;
            }
        }
        if count_it {
            saturated_carbo_rows += 1;
        }
    }
    u32::try_from(saturated_carbo_rows).map_err(|_| DescriptorError::CountOverflow {
        function: "num_saturated_carbocycles",
        field: "saturated_carbocycle_rows",
    })
}

/// Number of saturated carbocycles over prepared FINAL input, reading the
/// SUPPLIED ring rows.
///
/// Thin delegate to [`num_saturated_carbocycles_kernel`] over
/// `input.ring_info()`, `input.topology().bonds` and
/// `input.topology().atoms` — the rows the caller supplied (not
/// universally SSSR). Complexity: the kernel's pass with early break.
pub fn num_saturated_carbocycles_prepared(input: &DescriptorInput<'_>) -> DescriptorResult<u32> {
    num_saturated_carbocycles_kernel(
        input.ring_info(),
        &input.topology().bonds,
        &input.topology().atoms,
    )
}

/// Version string the source exports next to `calcNumSpiroAtoms`
/// (`Lipinski.cpp:435`).
pub const NUM_SPIRO_ATOMS_VERSION: &str = "1.0.0";

/// The ONE spiro-atom kernel shared by the prepared and cold entrypoints.
///
/// Behavior review: reproduces the pair-scan body of `calcNumSpiroAtoms`
/// exactly — for every unordered pair (i < j) of the RingInfo's
/// `atom_rings()` rows, compute the intersection; a row pair whose
/// intersection has EXACTLY ONE member contributes that atom as a spiro
/// atom, appended to `atoms` iff not already present (linear `contains`
/// mirroring the source's `std::find` — a prepopulated vector keeps its
/// entries, dedups against them, and inflates the returned count exactly
/// as the source's out-parameter does). Fused pairs (|i∩j| = 2) and
/// cubane-style larger shares do NOT contribute. The kernel is PURE: the
/// source's ensure-SSSR acquisition precondition lives in the
/// entrypoints (see [`num_spiro_atoms_prepared`]); this mirrors the
/// source's structure where `rInfo` is read AFTER the ensure step. The
/// returned count is `atoms.len()` (INCLUDING preexisting entries)
/// widened fail-closed via `u32::try_from`.
///
/// Complexity review: R^2 unordered row pairs with a per-pair linear
/// membership intersect O(len_i * len_j) — the same cost class as the
/// source's `Intersect` (whose own inline EFF comment concedes it does
/// more work than required; comment preserved verbatim below); linear
/// `contains` dedup mirrors `std::find`.
pub(crate) fn spiro_atom_ids_kernel(
    ring_info: &RingInfo,
    atoms: &mut Vec<AtomId>,
) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:448-462):
    //   for (unsigned int i = 0; i < rInfo->atomRings().size(); ++i) {
    //     const INT_VECT &ri = rInfo->atomRings()[i];
    //     for (unsigned int j = i + 1; j < rInfo->atomRings().size(); ++j) {
    //       const INT_VECT &rj = rInfo->atomRings()[j];
    //       // EFF: using intersect here does more work and memory allocation than is
    //       // required
    //       INT_VECT inter;
    //       Intersect(ri, rj, inter);
    //       if (inter.size() == 1) {
    //         if (std::find(atoms->begin(), atoms->end(), inter[0]) == atoms->end()) {
    //           atoms->push_back(inter[0]);
    //         }
    //       }
    //     }
    //   }
    // RDKit✔️✔️: for (unsigned int i = 0; i < rInfo->atomRings().size(); ++i) {
    // RDKit✔️✔️:   const INT_VECT &ri = rInfo->atomRings()[i];
    // RDKit✔️✔️:   for (unsigned int j = i + 1; j < rInfo->atomRings().size(); ++j) {
    // RDKit✔️✔️:     const INT_VECT &rj = rInfo->atomRings()[j];
    // RDKit✔️✔️:     // EFF: using intersect here does more work and memory allocation than is
    // RDKit✔️✔️:     // required
    // RDKit✔️✔️:     INT_VECT inter;
    // RDKit✔️✔️:     Intersect(ri, rj, inter);
    // RDKit✔️✔️:     if (inter.size() == 1) {
    // RDKit✔️✔️:       if (std::find(atoms->begin(), atoms->end(), inter[0]) == atoms->end()) {
    // RDKit✔️✔️:         atoms->push_back(inter[0]);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    //
    // `atomRings()` are the rows of the RingInfo handed in by the caller
    // AFTER the entrypoint's ensure step; `atoms` mirrors the source's
    // optional out-parameter (dedup + count include its prior contents).
    //
    // --- Cross-file helper closure (RDGeneral/types.cpp:37-45), copied
    // verbatim per the protocol's cross-file helper rule beside the Rust
    // intersection that implements it (the `shared` collection below):
    //
    // RDKit source (RDGeneral/types.cpp:37-45):
    //   void Intersect(const INT_VECT &r1, const INT_VECT &r2, INT_VECT &res) {
    //     res.resize(0);
    //     INT_VECT_CI ri;
    //     for (ri = r1.begin(); ri != r1.end(); ri++) {
    //       if (std::find(r2.begin(), r2.end(), (*ri)) != r2.end()) {
    //         res.push_back(*ri);
    //       }
    //     }
    //   }
    // RDKit✔️✔️: void Intersect(const INT_VECT &r1, const INT_VECT &r2, INT_VECT &res) {
    // RDKit✔️✔️:   res.resize(0);
    // RDKit✔️✔️:   INT_VECT_CI ri;
    // RDKit✔️✔️:   for (ri = r1.begin(); ri != r1.end(); ri++) {
    // RDKit✔️✔️:     if (std::find(r2.begin(), r2.end(), (*ri)) != r2.end()) {
    // RDKit✔️✔️:       res.push_back(*ri);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    //
    // Helper behavior review: `Intersect` clears the output length then
    // keeps EVERY member of r1 (in r1's own order) that linear membership
    // search finds in r2 — the intersection preserves the FIRST row's
    // member order and may contain duplicates if r1 itself does. The Rust
    // `shared` collection below reproduces exactly this: first-row
    // iteration order and `contains` as the typed `std::find` linear
    // membership test, with the caller's local `inter` / `shared`
    // per-pair temporary container (`res.resize(0)` only clears the
    // length; it neither constructs nor allocates). The ID-order products
    // pin the PAIR-TRAVERSAL order and the OUTPUT-ID order; they do not
    // pin a multi-member intersection's internal member order beyond the
    // first-row-order rule shared with the source.
    //
    // Helper complexity review: O(len(r1) * len(r2)) per pair via the
    // nested linear membership scan; the per-pair temporary container is
    // the caller's local `inter` (source) / `shared` (Rust) — an empty
    // intersection needs no allocation and nonempty growth may allocate
    // or reallocate; no exactly-one-allocation claim is made. Identical
    // scan shape — no sort, no hash set, no asymptotic divergence; output
    // dedup is the linear `contains`/`std::find` scan.
    let rows = ring_info.atom_rings();
    for i in 0..rows.len() {
        for j in (i + 1)..rows.len() {
            let shared: Vec<_> = rows[i]
                .iter()
                .filter(|atom_id| rows[j].contains(atom_id))
                .collect();
            if shared.len() == 1 {
                let spiro = *shared[0];
                if !atoms.contains(&spiro) {
                    atoms.push(spiro);
                }
            }
        }
    }
    u32::try_from(atoms.len()).map_err(|_| DescriptorError::CountOverflow {
        function: "num_spiro_atoms",
        field: "spiro_atoms",
    })
}

/// Number of spiro atoms over prepared FINAL input, with the source's
/// ensure-SSSR acquisition branch.
///
/// Behavior review: mirrors `calcNumSpiroAtoms`'s precondition exactly —
/// the SUPPLIED rows are used as-is when they are SSSR-or-better
/// (`RingInfo::is_sssr_or_better`, the typed `isSssrOrBetter()`); ONLY a
/// weaker-or-absent ring set triggers the recomputation (one canonical
/// `find_sssr` on the input's topology, the `MolOps::findSSSR(mol)`
/// branch). No heuristic and no unconditional recompute.
///
/// Complexity review: at most ONE `find_sssr` acquisition plus the
/// kernel's R^2 pair scan.
pub fn num_spiro_atoms_prepared(input: &DescriptorInput<'_>) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:438-441):
    //   if (!mol.getRingInfo() || !mol.getRingInfo()->isSssrOrBetter()) {
    //     MolOps::findSSSR(mol);
    //   }
    //   const RingInfo *rInfo = mol.getRingInfo();
    // RDKit✔️✔️: if (!mol.getRingInfo() || !mol.getRingInfo()->isSssrOrBetter()) {
    // RDKit✔️✔️:   MolOps::findSSSR(mol);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: const RingInfo *rInfo = mol.getRingInfo();
    //
    // The prepared input always carries a RingInfo (never absent), so only
    // the weaker-than-SSSR arm can fire here; the recompute is the ONE
    // canonical cold owner on the input's own topology.
    let mut atoms = Vec::new();
    if input.ring_info().is_sssr_or_better() {
        spiro_atom_ids_kernel(input.ring_info(), &mut atoms)
    } else {
        let ensured = crate::ring_info(input.topology(), "num_spiro_atoms")?;
        spiro_atom_ids_kernel(&ensured, &mut atoms)
    }
}

/// Version string the source exports next to `calcNumBridgeheadAtoms`
/// (`Lipinski.cpp:465`).
pub const NUM_BRIDGEHEAD_ATOMS_VERSION: &str = "2.0.0";

/// The ONE bridgehead-atom kernel shared by the prepared and cold
/// entrypoints.
///
/// Behavior review: reproduces the pair-scan body of
/// `calcNumBridgeheadAtoms` exactly — for every unordered pair (i < j) of
/// the RingInfo's `bond_rings()` rows, compute the shared-BOND
/// intersection; a pair qualifies ONLY when it shares MORE THAN ONE bond
/// (fused pairs sharing a single bond never qualify). For each qualifying
/// pair a FRESH per-pair incidence array over all atoms counts endpoint
/// touches of the SHARED bonds (both endpoints of every shared bond,
/// mirroring the source's per-pair `atomCounts` constructor); an atom
/// with incidence EXACTLY ONE is a bridgehead (the shared bonds form a
/// chain — interior chain atoms touch two shared bonds, the two chain
/// ends touch one) and is appended to `atoms` iff absent (linear
/// `contains` mirroring `std::find` — preexisting contents are kept,
/// deduped against, and INCLUDED in the returned count). The kernel is
/// PURE: the ensure-SSSR acquisition precondition lives in the
/// entrypoints. The returned count is `atoms.len()` widened fail-closed
/// via `u32::try_from`.
///
/// Complexity review: R^2 unordered row pairs with per-pair linear
/// membership intersect; each qualifying pair costs O(inter_size)
/// endpoint increments plus an O(num_atoms) incidence scan and a fresh
/// zeroed allocation — the source allocates a fresh numAtoms-sized array
/// per qualifying pair and we keep that exact shape (fresh per pair, no
/// cross-pair accumulation, which is also semantically load-bearing).
pub(crate) fn bridgehead_atom_ids_kernel(
    ring_info: &RingInfo,
    bonds: &[cosmolkit_model::Bond],
    atom_count: usize,
    atoms: &mut Vec<AtomId>,
) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:479-500):
    //   for (unsigned int i = 0; i < rInfo->bondRings().size(); ++i) {
    //     const INT_VECT &ri = rInfo->bondRings()[i];
    //     for (unsigned int j = i + 1; j < rInfo->bondRings().size(); ++j) {
    //       const INT_VECT &rj = rInfo->bondRings()[j];
    //       // EFF: using intersect here does more work and memory allocation than is
    //       // required
    //       INT_VECT inter;
    //       Intersect(ri, rj, inter);
    //       if (inter.size() > 1) {
    //         INT_VECT atomCounts(mol.getNumAtoms(), 0);
    //         for (auto ii : inter) {
    //           atomCounts[mol.getBondWithIdx(ii)->getBeginAtomIdx()] += 1;
    //           atomCounts[mol.getBondWithIdx(ii)->getEndAtomIdx()] += 1;
    //         }
    //         for (unsigned int ti = 0; ti < atomCounts.size(); ++ti) {
    //           if (atomCounts[ti] == 1) {
    //             if (std::find(atoms->begin(), atoms->end(), ti) == atoms->end()) {
    //               atoms->push_back(ti);
    //             }
    //           }
    //         }
    //       }
    //     }
    //   }
    // RDKit✔️✔️: for (unsigned int i = 0; i < rInfo->bondRings().size(); ++i) {
    // RDKit✔️✔️:   const INT_VECT &ri = rInfo->bondRings()[i];
    // RDKit✔️✔️:   for (unsigned int j = i + 1; j < rInfo->bondRings().size(); ++j) {
    // RDKit✔️✔️:     const INT_VECT &rj = rInfo->bondRings()[j];
    // RDKit✔️✔️:     // EFF: using intersect here does more work and memory allocation than is
    // RDKit✔️✔️:     // required
    // RDKit✔️✔️:     INT_VECT inter;
    // RDKit✔️✔️:     Intersect(ri, rj, inter);
    // RDKit✔️✔️:     if (inter.size() > 1) {
    // RDKit✔️✔️:       INT_VECT atomCounts(mol.getNumAtoms(), 0);
    // RDKit✔️✔️:       for (auto ii : inter) {
    // RDKit✔️✔️:         atomCounts[mol.getBondWithIdx(ii)->getBeginAtomIdx()] += 1;
    // RDKit✔️✔️:         atomCounts[mol.getBondWithIdx(ii)->getEndAtomIdx()] += 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       for (unsigned int ti = 0; ti < atomCounts.size(); ++ti) {
    // RDKit✔️✔️:         if (atomCounts[ti] == 1) {
    // RDKit✔️✔️:           if (std::find(atoms->begin(), atoms->end(), ti) == atoms->end()) {
    // RDKit✔️✔️:             atoms->push_back(ti);
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    //
    // `bondRings()` are the rows of the RingInfo handed in by the caller
    // AFTER the entrypoint's ensure step; `bonds` is the caller's bond
    // slice (typed begin()/end() endpoint reads); `atoms` mirrors the
    // source's optional out-parameter.
    let rows = ring_info.bond_rings();
    for i in 0..rows.len() {
        for j in (i + 1)..rows.len() {
            let shared: Vec<_> = rows[i]
                .iter()
                .filter(|bond_id| rows[j].contains(bond_id))
                .collect();
            if shared.len() > 1 {
                let mut atom_counts = vec![0usize; atom_count];
                for bond_id in &shared {
                    let bond = &bonds[bond_id.index()];
                    atom_counts[bond.begin().index()] += 1;
                    atom_counts[bond.end().index()] += 1;
                }
                for (atom_index, count) in atom_counts.iter().enumerate() {
                    if *count == 1 {
                        let bridgehead = AtomId::new(atom_index);
                        if !atoms.contains(&bridgehead) {
                            atoms.push(bridgehead);
                        }
                    }
                }
            }
        }
    }
    u32::try_from(atoms.len()).map_err(|_| DescriptorError::CountOverflow {
        function: "num_bridgehead_atoms",
        field: "bridgehead_atoms",
    })
}

/// Number of bridgehead atoms over prepared FINAL input, with the
/// source's ensure-SSSR acquisition branch (identical shape to
/// [`num_spiro_atoms_prepared`]).
pub fn num_bridgehead_atoms_prepared(input: &DescriptorInput<'_>) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:468-471):
    //   if (!mol.getRingInfo() || !mol.getRingInfo()->isSssrOrBetter()) {
    //     MolOps::findSSSR(mol);
    //   }
    //   const RingInfo *rInfo = mol.getRingInfo();
    // RDKit✔️✔️: if (!mol.getRingInfo() || !mol.getRingInfo()->isSssrOrBetter()) {
    // RDKit✔️✔️:   MolOps::findSSSR(mol);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: const RingInfo *rInfo = mol.getRingInfo();
    let mut atoms = Vec::new();
    if input.ring_info().is_sssr_or_better() {
        bridgehead_atom_ids_kernel(
            input.ring_info(),
            &input.topology().bonds,
            input.topology().atoms.len(),
            &mut atoms,
        )
    } else {
        let ensured = crate::ring_info(input.topology(), "num_bridgehead_atoms")?;
        bridgehead_atom_ids_kernel(
            &ensured,
            &input.topology().bonds,
            input.topology().atoms.len(),
            &mut atoms,
        )
    }
}
