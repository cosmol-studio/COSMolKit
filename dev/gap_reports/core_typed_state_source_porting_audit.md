# Core Typed State Source-Porting Audit

This is the Stage 1 typed-state audit for
[`porting_plan.md`](../archive/roadmaps/porting_plan.md).

No runtime behavior is changed by this document. It records current redesigned
Rust core state and the typed-state gaps that should be resolved before marking
checked README source-porting items complete.

## Evidence Snapshot

- Normative documents reread before this audit:
  `dev/README.md`, `dev/policy_invariants.md`,
  `dev/topology_operations.md`,
  `dev/operation_system_standard.md`, and
  `dev/testing_contract.md`.
- Inspected modules:
  `atom.rs`, `bond.rs`, `molecule.rs`, `derived.rs`, `stereo.rs`,
  `sgroup.rs`, `query.rs`, `ops.rs`, `smiles.rs`, `io/sdf.rs`,
  `hydrogens.rs`, `builder.rs`.
- Scope is Rust-only and source-porting only. No new parity tests are proposed
  here.

## Current Typed State

| Area | Current typed state | Remap / invalidation status |
|---|---|---|
| Atom identity | `AtomId`, `Element`, atom table rows | Atom ids are reassigned through `TopologyBlock::remove_atoms_with_mapping`; appended atoms use `TopologyMapping::with_appended`. |
| Atom chemical state | formal charge, explicit H count, chiral tag, chiral permutation, unknown-stereo flag, molfile parity, molfile inversion flag, implicit-hydrogen origin flag, tracked isotopic hydrogens, aromatic flag, isotope, atom map, no-implicit, radical electrons, hybridization | Stored on `Atom`; cloned/remapped with atom rows. Weak operations mutate through `OpParts::topology_mut()`. |
| Atom query state | `QueryNode<AtomQueryPredicate>` | Stored on atom rows; remaps with atom rows. Query semantics remain subset/preserved-only. |
| Atom properties | `BTreeMap<String, String>` | Stored on atom rows; remaps with atom rows. Some RDKit compatibility state still lives here and may need typed migration. |
| PDB residue info | `AtomPdbResidueInfo` | Typed atom state; remaps with atom rows. AddHs now plans and applies typed updates through operation assignment. |
| Bond identity | `BondId`, begin/end atoms | Bond ids/endpoints remapped in `TopologyBlock::remove_atoms_with_mapping`. |
| Bond chemical state | order, aromatic, conjugated, direction, stereo, stereo atom refs, unknown-stereo flag | Stored on bond rows; stereo atom refs remap when bonds survive compaction. |
| Bond query state | `QueryNode<BondQueryPredicate>` | Stored on bond rows; remaps with bond rows. Query semantics remain subset/preserved-only. |
| Bond properties | `BTreeMap<String, String>` | Stored on bond rows; remaps with bond rows. Some stereo/source flags still live here as props. |
| Coordinates | optional 2D rows, 3D conformer rows, source coordinate dimension | `CoordinateBlock::remap_topology` remaps coordinates for registered strong ops with `auto_remap: [coordinates]`. |
| Molecule properties | name, SDF data fields, SDF property lists, molecule props | Stored in `MoleculeProperties`; molecule-level props preserve by Arc clone or property COW. Atom/bond property lists need explicit policy for topology edits. |
| SGroups | typed `SubstanceGroup` rows | `SubstanceGroup::remapped` remaps atoms/bonds/parent atoms/attach points/cstates. Parent SGroup remap needs audit, see gaps. |
| Stereo groups | typed `StereoGroup` rows | Remapped through atom and bond maps; groups that lose referenced rows are dropped. |
| Derived cache | adjacency, rings, ring families, valence, aromaticity-valid, stereo-valid | Registered operations invalidate/recompute through `DerivedState`. Drawing/fingerprint bits exist but no cache blocks exist yet. |

## Typed-State Gaps

These are not permission to add fields immediately. Each addition still requires
rereading the normative documents and adding fields only when the source port
needs them.

| Gap | Current representation | Why it matters | Required policy-compliant path |
|---|---|---|---|
| Chiral permutation | Typed atom field `chiral_permutation: Option<u32>` | SMILES stereo parsing, non-tetrahedral permutation adjustment, and sanitize cleanup now read/write typed atom state. Existing RDKit source markers still name `_chiralPermutation` because that is the upstream property. | Continue using typed state for operation behavior. Reintroduce a compatibility prop only if a source-aligned external preservation path requires it. |
| Molfile parity / inversion flags | Typed atom fields `mol_parity` and `mol_inversion_flag` | V2000 atom-line and V3000 atom-property source paths now preserve these values as typed atom state. Existing RDKit markers still name `molParity` and `molInversionFlag` because that is upstream storage. | Continue using typed state for SDF/stereo behavior. Add compatibility serialization only when source-ported writer branches need external molfile property names. |
| Unknown stereo flags | Typed atom and bond `unknown_stereo` flags | RemoveHs writes atom unknown-stereo state; SDF single-bond direction cleanup writes bond unknown-stereo state. Existing RDKit markers still name `_UnknownStereo` because that is upstream storage. | Continue using typed state for operation and IO behavior. Add compatibility serialization only if a source-aligned writer path requires preserving the upstream property name externally. |
| AddHs implicit marker | Typed atom field `implicit_hydrogen` | AddHs marks generated implicit hydrogens as typed atom state; RemoveHs predicates read the typed flag instead of a string prop. Existing RDKit markers still name `isImplicit` because that is upstream storage. | Continue using typed state for AddHs/RemoveHs behavior. Add compatibility serialization only if a source-aligned external format requires it. |
| Isotopic hydrogen tracking | Typed atom field `tracked_isotopic_hydrogens: Vec<u16>` | RemoveHs records isotopic hydrogens on heavy atoms as typed state; AddHs replays and clears that typed state. Existing RDKit markers still name `_isotopicHs` because that is upstream storage. | Continue using typed state for AddHs/RemoveHs behavior. Add compatibility serialization only if a source-aligned external format requires it. |
| SDF atom/bond property lists | `MoleculeProperties::sdf_property_lists` with target and row-aligned values | Strong topology edits now auto-remap properties for registered hydrogen operations. Deleted atom/bond rows are filtered; appended atom/bond rows receive `None` values. | Keep extending tests as new strong topology operations are registered. |
| SGroup parent remap | `TopologyBlock::remove_atoms_with_mapping` now builds an old-to-new SGroup map for surviving groups before remapping | Parent links are preserved when both parent and child survive atom/bond compaction. Broader SGroup remap scenarios still need source audit. | Keep this behavior covered by operation-contract tests and extend source audit for nested/partially removed SGroups. |
| Drawing cache | `DerivedState::DRAWING` bit only | Drawing is a checked README feature but has no typed cache/prepared-state storage yet. | Add drawing cache only if source port needs cached state. Otherwise keep drawing as recomputed output and use invalidation bit for future compatibility. |
| Fingerprint cache | `DerivedState::FINGERPRINT` bit only | Fingerprint is a checked README feature but generation is unsupported. | Add cache only after Morgan source port establishes what can become stale; otherwise compute directly and keep invalidation bit reserved. |
| DG bounds feature separation | `dg_bounds_matrix` returns a coordinate-2d unsupported error | DG bounds is a checked README feature distinct from 2D coordinate generation. | Add or split a dedicated support feature before marking DG bounds complete. |
| PDB residue info completeness | `AtomPdbResidueInfo` models the subset needed by AddHs | RDKit `AtomPDBResidueInfo` also has secondary structure, segment number, and monomer class. | Add fields only when Rust source ports for PDB/sequence/writer behavior need them. Do not add broad compatibility fields without source need. |

## Remap Risks To Resolve Before Completion

1. Strong topology operations must account for molecule-level row-aligned data.
   Registered hydrogen strong operations now auto-remap SDF atom/bond property
   lists through the properties block. Future strong operations must declare
   and implement the same policy before claiming SDF roundtrip support.
2. SGroup parent relationships now use a real SGroup mapping during atom
   compaction. Broader nested and partially removed SGroup cases still need
   source-level audit before SDF roundtrip claims are complete.
3. Stereo state is split between typed fields and compatibility props. Chiral
   permutation, unknown-stereo flags, molfile parity, and molfile inversion are
   typed. Remaining stereo/source props should stay compatibility props until a
   source branch needs operation-visible state.
4. Drawing, fingerprint, and DG bounds are checked README areas but have no
   implementation state beyond public data types or reserved derived bits.

## Recommended Next Implementation Targets

1. Continue hydrogens source porting beyond typed state: `AddHsParams` and
   `RemoveHsParams` are now public Rust parameters for registered operations,
   including `add_coords`, `add_residue_info`, and `remove_and_track_isotopes`.
   Remaining work should focus on deeper source-aligned coordinate geometry,
   stereo side effects, and sanitize-after-removal coverage.
