//! Stereo-count descriptor owner (RDKit `Lipinski.cpp` stereo functions).

use std::borrow::Cow;

use cosmolkit_model::ChiralTag;
use cosmolkit_model::TopologyBlock;

use crate::{DescriptorError, DescriptorInput, DescriptorResult};

/// The `_ChiralityPossible` atom computed-property key the legacy stereo
/// owner sets (`legacy_stereo.rs` mirrors RDKit
/// `common_properties::_ChiralityPossible`).
const CHIRALITY_POSSIBLE: &str = "_ChiralityPossible";

/// The molecule-level stereo-assigned marker key (RDKit
/// `common_properties::_StereochemDone`).
const STEREOCHEM_DONE: &str = "_StereochemDone";

/// Ensures stereochemistry is assigned, mirroring the source dispatch.
///
/// Behavior review: reproduces the anonymous-namespace
/// `hasStereoAssigned` helper and the copy-and-assign branch of
/// `numAtomStereoCenters` / `numUnspecifiedAtomStereoCenters` exactly —
/// the guard is a PRESENCE check on the molecule-level `_StereochemDone`
/// property (present with ANY value, including "0", counts as assigned;
/// RDKit writes it as 1 when done). When the property is ABSENT the
/// source CLONES the molecule and runs
/// `MolOps::assignStereochemistry(copy, cleanIt=true, force=true,
/// flagPossible=true)`; the COSMolKit mapping passes the input topology
/// BY VALUE to the legacy owner (the owned parameter IS the copy) with
/// both independent flags true — `clean_it = true` and
/// `flag_possible_stereo_centers = true` — and `force = true` is exactly
/// the unconditional run this caller-side guard models (the detached
/// topology carries no `_StereochemDone`, so the owner's documented
/// caller-owned force guard is discharged here). The input's SUPPLIED
/// valence rows and ring rows are reused, never recomputed.
///
/// Complexity review: O(1) presence check on the assigned arm; the absent
/// arm costs ONE topology clone plus ONE legacy stereo assignment (the
/// source's own copy + assign shape — no additional passes).
pub(crate) fn ensured_stereo_topology<'a>(
    input: &'a DescriptorInput<'_>,
    function: &'static str,
) -> Result<Cow<'a, TopologyBlock>, DescriptorError> {
    // RDKit source (Lipinski.cpp:504-506):
    //   namespace {
    //   bool hasStereoAssigned(const ROMol &mol) {
    //     return mol.hasProp(common_properties::_StereochemDone);
    //   }
    //   }  // namespace
    // RDKit✔️✔️: bool hasStereoAssigned(const ROMol &mol) {
    // RDKit✔️✔️:   return mol.hasProp(common_properties::_StereochemDone);
    // RDKit✔️✔️: }
    //
    // Presence semantics: `hasProp` is true for any stored value; a
    // `_StereochemDone` of "0" still means "assigned" exactly as in the
    // source.
    if input.properties().prop(STEREOCHEM_DONE).is_some() {
        return Ok(Cow::Borrowed(input.topology()));
    }
    // RDKit source (Lipinski.cpp:510-517):
    //   std::unique_ptr<ROMol> tmol;
    //   const ROMol *mptr = &mol;
    //   if (!hasStereoAssigned(mol)) {
    //     tmol.reset(new ROMol(mol));
    //     constexpr bool cleanIt = true;
    //     constexpr bool force = true;
    //     constexpr bool flagPossible = true;
    //     MolOps::assignStereochemistry(*tmol, cleanIt, force, flagPossible);
    //     mptr = tmol.get();
    //   }
    // RDKit✔️✔️: if (!hasStereoAssigned(mol)) {
    // RDKit✔️✔️:   tmol.reset(new ROMol(mol));
    // RDKit✔️✔️:   constexpr bool cleanIt = true;
    // RDKit✔️✔️:   constexpr bool force = true;
    // RDKit✔️✔️:   constexpr bool flagPossible = true;
    // RDKit✔️✔️:   MolOps::assignStereochemistry(*tmol, cleanIt, force, flagPossible);
    // RDKit✔️✔️:   mptr = tmol.get();
    // RDKit✔️✔️: }
    //
    // `new ROMol(mol)` maps to passing the topology BY VALUE to the
    // legacy owner (its owned parameter is the copy); `force = true` is
    // the unconditional run this dispatch performs.
    let assigned = cosmolkit_core::assign_legacy_stereochemistry_with_flags(
        input.topology().clone(),
        input.valence(),
        input.ring_info(),
        true,
        true,
    )
    .map_err(|source| DescriptorError::Stereo { function, source })?;
    Ok(Cow::Owned(assigned))
}

/// Version string the source exports next to `numAtomStereoCenters`
/// (`Lipinski.cpp:509`; the S01 record cited 507 — comment-only
/// correction, the value "1.0.1" was always exact).
pub const NUM_ATOM_STEREO_CENTERS_VERSION: &str = "1.0.1";

/// The ONE atom-stereo-center count kernel.
///
/// Behavior review: reproduces the counting loop of
/// `numAtomStereoCenters` exactly — one pass over the (possibly
/// stereo-ensured) topology's atoms; an atom counts iff it carries the
/// `_ChiralityPossible` computed property (set by the legacy stereo
/// owner, mirroring `common_properties::_ChiralityPossible`). The count
/// is widened fail-closed via `u32::try_from` (same value on every
/// realizable input).
///
/// Complexity review: O(atoms) with an O(log n) per-atom property
/// lookup (BTreeMap-backed props) — the source's linear atom scan over a
/// property-carrying molecule; no allocation beyond the `Result`.
pub(crate) fn num_atom_stereo_centers_kernel(topology: &TopologyBlock) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:519-524):
    //   unsigned int res = 0;
    //   for (const auto &atom : mptr->atoms()) {
    //     if (atom->hasProp(common_properties::_ChiralityPossible)) {
    //       ++res;
    //     }
    //   }
    //   return res;
    // RDKit✔️✔️: unsigned int res = 0;
    // RDKit✔️✔️: for (const auto &atom : mptr->atoms()) {
    // RDKit✔️✔️:   if (atom->hasProp(common_properties::_ChiralityPossible)) {
    // RDKit✔️✔️:     ++res;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: return res;
    let mut centers = 0usize;
    for atom in &topology.atoms {
        if atom.prop(CHIRALITY_POSSIBLE).is_some() {
            centers += 1;
        }
    }
    u32::try_from(centers).map_err(|_| DescriptorError::CountOverflow {
        function: "num_atom_stereo_centers",
        field: "atom_stereo_centers",
    })
}

/// Number of potential atom stereo centers over prepared FINAL input,
/// with the source's assigned-state dispatch.
pub fn num_atom_stereo_centers_prepared(input: &DescriptorInput<'_>) -> DescriptorResult<u32> {
    let topology = ensured_stereo_topology(input, "num_atom_stereo_centers")?;
    num_atom_stereo_centers_kernel(&topology)
}

/// Version string the source exports next to
/// `numUnspecifiedAtomStereoCenters` (`Lipinski.cpp:530`; the S01 record
/// cited 527 — comment-only correction, the value "1.0.1" was always
/// exact).
pub const NUM_UNSPECIFIED_ATOM_STEREO_CENTERS_VERSION: &str = "1.0.1";

/// The ONE unspecified-atom-stereo-center count kernel.
///
/// Behavior review: reproduces the counting loop of
/// `numUnspecifiedAtomStereoCenters` exactly — identical dispatch, and an
/// atom counts iff it carries `_ChiralityPossible` AND its
/// `chiral_tag()` is `ChiralTag::Unspecified` (the typed
/// `Atom::CHI_UNSPECIFIED` projection). Widened fail-closed via
/// `u32::try_from`.
///
/// Complexity review: O(atoms) with the same per-atom property lookup
/// plus an O(1) tag read — the source's linear scan; no allocation
/// beyond the `Result`.
pub(crate) fn num_unspecified_atom_stereo_centers_kernel(
    topology: &TopologyBlock,
) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:546-551):
    //   unsigned int res = 0;
    //   for (const auto &atom : mptr->atoms()) {
    //     if (atom->hasProp(common_properties::_ChiralityPossible) &&
    //         atom->getChiralTag() == Atom::CHI_UNSPECIFIED) {
    //       ++res;
    //     }
    //   }
    //   return res;
    // RDKit✔️✔️: unsigned int res = 0;
    // RDKit✔️✔️: for (const auto &atom : mptr->atoms()) {
    // RDKit✔️✔️:   if (atom->hasProp(common_properties::_ChiralityPossible) &&
    // RDKit✔️✔️:     atom->getChiralTag() == Atom::CHI_UNSPECIFIED) {
    // RDKit✔️✔️:     ++res;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: return res;
    let mut centers = 0usize;
    for atom in &topology.atoms {
        if atom.prop(CHIRALITY_POSSIBLE).is_some() && atom.chiral_tag() == ChiralTag::Unspecified {
            centers += 1;
        }
    }
    u32::try_from(centers).map_err(|_| DescriptorError::CountOverflow {
        function: "num_unspecified_atom_stereo_centers",
        field: "unspecified_atom_stereo_centers",
    })
}

/// Number of potential-but-unspecified atom stereo centers over prepared
/// FINAL input, with the source's assigned-state dispatch.
pub fn num_unspecified_atom_stereo_centers_prepared(
    input: &DescriptorInput<'_>,
) -> DescriptorResult<u32> {
    let topology = ensured_stereo_topology(input, "num_unspecified_atom_stereo_centers")?;
    num_unspecified_atom_stereo_centers_kernel(&topology)
}
