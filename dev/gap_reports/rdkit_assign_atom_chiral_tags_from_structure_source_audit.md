# RDKit `AssignAtomChiralTagsFromStructure` Source Audit

## Audit Scope

This audit covers the pinned RDKit `2026.3.1` behavior exposed to Python as
`Chem.AssignAtomChiralTagsFromStructure`. The Python wrapper at
`third_party/rdkit/Code/GraphMol/Wrap/MolOps.cpp:2390-2409` binds that name
directly to `RDKit::MolOps::assignChiralTypesFrom3D`; it is not a wrapper for a
separate C++ function named `assignAtomChiralTagsFromStructure`.

The declaration is `third_party/rdkit/Code/GraphMol/MolOps.h:1089-1105` and the
active implementation is
`third_party/rdkit/Code/GraphMol/Chirality.cpp:3377-3500`. The public parameters
are `confId=-1` and `replaceExistingTags=true`. The function deliberately does
not check whether four substituents are distinct, so it may assign a tag to an
atom that is not an actual stereocenter.

No compile-time branch changes this function in the pinned Linux build. The
non-tetrahedral branch is controlled at run time by
`RDK_ENABLE_NONTETRAHEDRAL_STEREO` and defaults to enabled.

## Verified Active Call Graph

The following closure was read from the definitions and call expressions, not
inferred from names.

```text
Python Chem.AssignAtomChiralTagsFromStructure
  MolOps::assignChiralTypesFrom3D
    ROMol::getNumConformers
    ROMol::getConformer
    Conformer::is3D
    RDProps::{hasProp,clearProp,setProp,getPropIfPresent}
    Chirality::getAllowNontetrahedralChirality
      getValFromEnvironment
        std::getenv
        std::strcmp
    Chirality::detail::getAtomNonzeroDegree
      ROMol::atomBonds
      Chirality::detail::bondAffectsAtomChirality
    Atom::getTotalNumHs(false)
      Atom::getNumExplicitHs
      Atom::getNumImplicitHs
        Atom::getValence(IMPLICIT)
    assignNontetrahedralChiralTypeFrom3D
      ROMol::atomBonds
      isWigglyBond
        Bond accessors
        RDProps::getPropIfPresent<int>
      ROMol::getAtomNeighbors
      Conformer::getAtomPos
      RDGeom::Point3D::directionVector
        RDGeom::Point3D::normalize
          RDGeom::Point3D::length
            sqrt
          zero_tolerance
      RDGeom::Point3D::dotProduct
      RDGeom::Point3D::angleTo
        RDGeom::Point3D::lengthSq
        RDGeom::Point3D::dotProduct
        sqrt
        acos
      VOLTEST
        RDGeom::Point3D::crossProduct
        RDGeom::Point3D::dotProduct
      OctahedralPermFrom3D
        VOLTEST
    tetrahedral signed-volume branch
      ROMol::atomBonds
      isWigglyBond
      Chirality::detail::bondAffectsAtomChirality
      Conformer::getAtomPos
      RDGeom::Point3D::operator-
      RDGeom::Point3D::crossProduct
      RDGeom::Point3D::dotProduct
```

`std::getenv`, `std::strcmp`, `sqrt`, and `acos` are libc/compiler math or
environment primitives, not chemistry functions to port. Atom, bond,
conformer, property, and graph iterator accessors map to existing COSMolKit
typed storage, but their ordering and exception behavior are part of the port.

### Source locations

| Behavior | Active source |
|---|---|
| environment constants/default | `GraphMol/Chirality.h:33-36` |
| `getValFromEnvironment` | `GraphMol/Chirality.cpp:45-54` |
| `getAllowNontetrahedralChirality` | `GraphMol/Chirality.cpp:860-863` |
| `bondAffectsAtomChirality` | `GraphMol/Chirality.cpp:878-889` |
| `getAtomNonzeroDegree` | `GraphMol/Chirality.cpp:890-900` |
| `VOLTEST` | `GraphMol/Chirality.cpp:3114` |
| `OctahedralPermFrom3D` | `GraphMol/Chirality.cpp:3116-3171` |
| `isWigglyBond` | `GraphMol/Chirality.cpp:3173-3184` |
| non-tetrahedral assignment | `GraphMol/Chirality.cpp:3188-3375` |
| top-level assignment | `GraphMol/Chirality.cpp:3377-3500` |
| total hydrogen count | `GraphMol/Atom.cpp:286-294` |
| conformer selection | `GraphMol/ROMol.cpp:630-649` |
| graph adjacency access | `GraphMol/ROMol.h:283-305`, `ROMol.cpp:352-364` |
| `Point3D` normalization/direction | `Geometry/point.h:147-161`, `210-220` |
| `Point3D` dot/angle/cross | `Geometry/point.h:169-193`, `223-234` |
| `zero_tolerance` | `Numerics/Vector.h:24` (`1.e-16`) |

## Active Branch Inventory

### Entry and conformer selection

- Zero conformers is an unchanged, successful no-op.
- Negative `confId` selects the first conformer. A nonnegative id selects by
  conformer id, not vector position. A missing id throws
  `ConformerException("Can't find conformation with ID: N")`.
- A selected non-3D conformer is an unchanged, successful no-op. In particular,
  `_StereochemDone` is not cleared on this path.
- For a 3D conformer, `_StereochemDone` is cleared before atom processing.
- Coordinate lookup is by atom index. Invalid model alignment is a source
  precondition and must become a structured Rust invariant/stereo error.

### Environment behavior

- Missing `RDK_ENABLE_NONTETRAHEDRAL_STEREO` returns the default `true`.
- The exact C string `"0"` returns `false`.
- Every other present value, including empty, `"00"`, and `"false"`, returns
  `true`.
- The value is sampled once per top-level call. Production code must read the
  process environment directly and must not launch a subprocess.

### Bond and degree behavior

- `UNSPECIFIED` and `ZERO` never contribute to nonzero degree.
- A `DATIVE` bond does not contribute at its begin atom and does contribute at
  its end atom. Other modeled bond types contribute.
- A wiggly bond suppresses assignment only when the center is the begin atom,
  the bond is `SINGLE`, and either direction is `UNKNOWN` or integer property
  `_UnknownStereo` exists with a nonzero value.
- Incident bonds and neighbors are visited in the owning graph's adjacency
  order. COSMolKit's CSR adjacency is built in bond-table insertion order and
  can reproduce the modeled order without repeated full-table scans.

### Non-tetrahedral behavior

- Atomic numbers below 15 are rejected before geometry work.
- Every topology neighbor contributes a normalized direction vector; this
  helper does not filter zero or dative bonds after its wiggly check.
- More than six neighbors and fewer than three neighbors return false.
- Duplicate center/neighbor coordinates throw
  `runtime_error("Cannot normalize a zero length vector")`. Lengths strictly
  below `1.e-16` throw; equality does not.
- Antiparallel pairing uses strict `dot < -0.9`; equality is not paired. Reuse
  of either endpoint by a second pair rejects the assignment.
- The active shape branches are: three-coordinate T-shaped square planar;
  four-coordinate see-saw split at strict angle `< 100*pi/180`; five-coordinate
  trigonal bipyramidal; four-coordinate/two-pair square planar;
  five-coordinate/two-pair square pyramidal represented as octahedral; and
  six-coordinate/three-pair octahedral.
- `OctahedralPermFrom3D` selects permutations 1 through 30 through its nested
  pair-table switches. Its trailing `return 0` is source-unreachable for valid
  tables and must not be exposed as a fallback.
- `VOLTEST` is `dot(v[X], cross(v[Y], v[Z])) >= 0.0`. Consequently both signed
  zero values select true and NaN selects false.

### Tetrahedral and mutation behavior

- Explicitness is snapshotted before any tags are cleared. It includes begin
  atoms of `BEGINWEDGE`/`BEGINDASH` bonds and every atom initially carrying a
  non-unspecified tag.
- With `replaceExistingTags=false`, an initially tagged atom is skipped without
  mutation. Otherwise each atom's tag is set to unspecified immediately,
  before all degree and geometry rejection branches.
- Nonzero degree must be at least three and nonzero degree plus explicit and
  implicit H must not exceed six.
- If non-tetrahedral assignment succeeds it wins before tetrahedral handling.
- Tetrahedral handling requires total degree at most four. Except for sulfur
  (16) and selenium (34), total degree must equal four.
- The first three affected neighbors determine signed volume. Strict values
  below `-0.1` assign clockwise and above `0.1` assign counterclockwise.
  Equality and the interval between remain ambiguous.
- When the first volume is ambiguous and a fourth explicit neighbor exists,
  the source computes the negated volume using that fourth neighbor. It does
  not use a tolerance-based approximation.
- A successful geometry assignment writes `_NonExplicit3DChirality=1` only
  when the atom was not in the explicitness snapshot.
- Existing `_chiralPermutation`, `_NonExplicit3DChirality`, and unrelated
  properties are not globally cleared by this function. Their exact
  preservation or overwrite follows source mutation order.

## Pre-Port Baseline Gaps (Closed)

This section records the implementation state found by the initial audit. It is
historical evidence for the work selected by the plan; every gap listed here
was closed by the completed kernel, operation, parser, and oracle work described
below.

### Fabricated public helper

`crates/cosmolkit-core/src/chemistry/stereo.rs:3760-3834` defines
`assign_atom_chiral_tags_from_structure(&Molecule) -> Vec<(usize, ChiralTag)>`.
Its abbreviated C++ frame names a C++ function that does not exist at the cited
location. It infers tags from 2D wedge/dash state, never reads a conformer, and
uses an `is_opposite_bonds` helper that always returns true. Its current
`RDKit✔️✔️` marker is invalid. It has no production callers in the repository,
but is publicly reachable through the public `chemistry::stereo` module and
must be removed or replaced, not retained as a compatibility fallback.

### Substantial internal port

`crates/cosmolkit-core/src/notation/smiles/stereo.rs:2230-3067` contains most
numeric and shape branches. It is not yet exact:

| Source behavior | Pre-port Rust mismatch |
|---|---|
| signed `confId`, default `-1` | accepts `usize`; cannot express default selection |
| missing nonnegative id throws | returns `Ok(())` |
| public `replaceExistingTags` | parameter absent; always replaces |
| non-3D no-op before cache clear | clears `_StereochemDone` first |
| environment-gated non-tetra branch | always enabled |
| begin-relative wiggly predicate | treats either endpoint as wiggly |
| `_UnknownStereo` integer property | typed boolean collapses missing and zero/nonzero representation |
| center-relative dative predicate | helper excludes dative aliases without atom direction and fails ordinary `Dative` begin/end semantics |
| normalized length `< 1.e-16` throws | only exact zero becomes `None` |
| zero-length exception | silently converts to no assignment |
| Point3D `angleTo` arithmetic | substitutes a zero-vector result and uses Rust `clamp`, changing NaN behavior |
| source atom-by-atom mutation | computes all tags, then writes in a second phase |
| source entry uses cached H state | invokes whole-molecule valence assignment and may return a new `ValenceError` |
| source property preservation | delayed arrays do not reproduce every partial mutation/error path |
| adjacency iteration | repeatedly scans the whole bond table per atom |
| fallible operation | wrapper discards every error |

The helper frames for `VOLTEST`, `OctahedralPermFrom3D`, and most shape
branches are structurally complete, but their first-axis complete markers are
premature because their numeric dependencies and a direct all-branch oracle
have not been closed.

### Parser adapter duplicate

`crates/cosmolkit-core/src/io/sdf/postprocess.rs:1218-1365` originally repeated
the top-level C++ frame and delegated to the notation helper. The final state
keeps only parser sequencing/error adapters: SDF, MOL2, PDB-molecule, and 3D
SMILES all call `chemistry::stereo::assign_chiral_types_from_3d_molecule`, and
all propagate a parser-specific structured error instead of discarding it.

## Operation-Contract Audit

The completed public operation is registered through `molecule_ops!` as
`with_chiral_tags_from_structure`, with trailing-underscore in-place form
`assign_chiral_tags_from_structure_`. It is a weak topology-state operation
with `RequiredNow` parity and a dedicated operation-contract profile. The body
uses `OpParts`, performs all source work on detached topology/property blocks,
and commits only after successful completion.

The required operation is a weak topology-state operation: atom and bond
indices, topology, coordinates, and bond state remain unchanged; atom chiral
tags and stereo properties plus molecule `_StereochemDone` may change. It must
read topology, coordinates, properties, and the existing valence/H cache, and
write topology atom stereo state and properties. It requires trusted topology,
coordinate-row alignment, and source-defined conformer selection. It preserves
coordinates, atom/bond ordering, unrelated properties, and all non-stereo
derived state. It must provide value-style and trailing-underscore in-place
forms and use `STEREO_FEATURE` with required pinned-RDKit parity.

`StereoError` now has explicit `ConformerNotFound`,
`ConformerAtomCountMismatch`, `ZeroLengthVector`,
`MissingImplicitHydrogenState`, and `ImplicitHydrogenCountMismatch` paths. The
public operation and parser adapters preserve these as visible errors rather
than translating them into no-ops or silently discarding them.

## Pre-Port Test Coverage

- `crates/cosmolkit-core/src/notation/smiles/tests.rs` covers several ordinary
  tetrahedral, sulfur/selenium, non-3D, non-tetrahedral, and wiggly examples at
  internal/parser boundaries. It does not exercise the public parameter and
  full state matrix.
- `crates/cosmolkit-core/src/io/sdf.rs` contains parser-focused 3D assignment
  tests, including tetrahedral and square-planar examples. These verify SDF
  sequencing, not the standalone operation.
- `crates/cosmolkit-core/tests/tetrahedral_stereo_geometry.rs` compares ordered
  tetrahedral ligand reporting against RDKit ETKDG coordinates. It does not
  call or compare `AssignAtomChiralTagsFromStructure`.
- `tools/testdata/rdkit/_generate_tetrahedral_stereo_geometry.py` currently
  emits only SMILES, embedding status, initially tagged centers, and positions.
  It has no operation options, before/after state, errors, non-tetrahedral
  permutations, or property fields.
- RDKit's own tests cover ordinary R/S geometry (`testChirality.cpp:1100+`),
  sulfur/selenium (`molopstest.cpp:6814+`), wiggly bonds
  (`catch_chirality.cpp:2617+`), explicit hydrogens
  (`catch_chirality.cpp:2810+`), and additional parser/non-tetrahedral cases.
  They do not by themselves compare COSMolKit or cover every boundary below.

## Required Fixed Oracle Matrix

The new fixture and generated expected rows must include, without filtering:

- no conformer; negative default selection; first conformer; selected later
  conformer by id; missing id exception; selected non-3D conformer;
- `replaceExistingTags` true and false with existing tetrahedral and
  non-tetrahedral tags;
- absent, `"0"`, empty, `"00"`, `"false"`, and ordinary enabled environment
  values in isolated processes;
- ordinary tetrahedral CW/CCW; implicit-H and explicit-H centers; S and Se;
  degree below three, total degree above four and six; exact and adjacent
  signed-volume boundaries `-0.1` and `0.1`; coplanar first three with decisive
  fourth; fully degenerate geometry;
- zero, unspecified, ordinary, and both center perspectives of dative bonds;
  begin/end wiggly bonds; unknown direction; missing/zero/nonzero
  `_UnknownStereo`; non-single unknown-direction bond;
- non-tetra atomic-number rejection; neighbor counts 2 through 7; dot-product
  values below/equal/above `-0.9`; reused antiparallel endpoints; all
  square-planar, see-saw angle, trigonal-bipyramidal, square-pyramidal, and
  octahedral pair-table branches; both `VOLTEST` signs;
- normalization lengths below/equal/above `1.e-16`, duplicate coordinates,
  signed zero, subnormal, maximum finite, NaN, and infinities where RDKit has
  defined observable behavior;
- preexisting `_StereochemDone`, `_chiralPermutation`,
  `_NonExplicit3DChirality`, and unrelated atom/molecule properties, including
  partial mutation before an exception.

Each expected row must record operation status, error type/text category,
selected conformer id, all atom tags, all chiral permutations, all
`_NonExplicit3DChirality` values, molecule `_StereochemDone`, atom and bond
counts/order/types/directions/unknown-stereo state, coordinates, conformer
dimensionality, and unrelated property preservation.

## Explicitly Out of Scope

The following functions are not transitive callees of this public entry and
must not be folded into the port:

- `MolOps::assignStereochemistryFrom3D`;
- double-bond direction detection or E/Z assignment;
- CIP rank/code assignment or filtering;
- validation that substituents are distinct;
- atropisomer perception;
- parser-specific cleanup and subsequent stereo sequencing, except for later
  migration of existing call sites to the completed operation.

## Completed Kernel And Oracle Result

The completed source closure is implemented in
`crates/cosmolkit-core/src/chemistry/stereo.rs`. The first-axis markers for
`getValFromEnvironment`, `getAllowNontetrahedralChirality`,
`bondAffectsAtomChirality`, `getAtomNonzeroDegree`, `isWigglyBond`, the used
`Point3D` arithmetic, `VOLTEST`, `OctahedralPermFrom3D`,
`assignNontetrahedralChiralTypeFrom3D`, the tetrahedral branch, and
`assignChiralTypesFrom3D` are now `✔️` only after their focused tests and the
complete public-operation oracle passed.

The exact parity command was:

```bash
cargo test -p cosmolkit-core --release --features op-contracts-strict \
  --test rdkit_assign_atom_chiral_tags_from_structure_parity -- --nocapture
```

It executed one integration test with no filtering or ignored tests:

```text
1 passed; 0 failed; 0 ignored; 0 measured; 0 filtered out
```

That test loaded and asserted all 77 fixed oracle records. The 73 successful
records matched the pinned RDKit result exactly for status, selected conformer,
complete atom and bond state, chiral tags, chiral permutations,
`_NonExplicit3DChirality`, `_StereochemDone`, conformer coordinates and
dimensionality, and unrelated typed properties. The four error records also
matched the official error surface exactly:

| Case | RDKit oracle | Structured Rust error | Source object after error |
|---|---|---|---|
| `missing_specific_conformer` | `ValueError: Bad Conformer Id` | `StereoError::ConformerNotFound` mapped to the same type/text | unchanged |
| `direction_length_below_zero_tolerance` | `RuntimeError: Cannot normalize a zero length vector` | `StereoError::ZeroLengthVector` mapped to the same type/text | unchanged |
| `direction_duplicate_coordinate` | `RuntimeError: Cannot normalize a zero length vector` | `StereoError::ZeroLengthVector` mapped to the same type/text | unchanged |
| `partial_mutation_before_later_exception` | `RuntimeError: Cannot normalize a zero length vector` | `StereoError::ZeroLengthVector` mapped to the same type/text | unchanged by the public Rust operation |

The final row exposes an intentional policy boundary. RDKit mutates its working
`ROMol` before the later geometry exception, whereas COSMolKit's registered
value-style/in-place operation computes into detached blocks and commits only
on success. COSMolKit therefore preserves the caller's source molecule on
error, as required by the operation transaction contract. This is not a silent
fallback and does not alter any successful result.

The CSR adjacency implementation now visits incident bonds in insertion order
and has the same `O(degree)` helper traversal as RDKit's adjacency iterator, so
the bond-iteration helpers use `✔️✔️`. The top-level operation retains
`✔️❌`: transactionality clones topology and property blocks before execution,
which is additional `O(atoms + bonds + properties)` work compared with direct
RDKit mutation. Environment lookup also retains `✔️❌` because `var_os`
materializes an owned value while `getenv` returns a borrowed process pointer.

No active branch in the fixed 77-record kernel matrix remains unclosed. Parser
call-site consolidation, Python exposure, and repository-wide validation are
also complete and are recorded below.

## Final Source-Frame Status

The verified call graph at the start of this report is the complete active
closure for the selected RDKit entry point. The corresponding Rust source
frames and final markers are:

| Source unit | Official source | Rust status | Reason for second axis |
|---|---|---|---|
| environment constants, `getValFromEnvironment`, `getAllowNontetrahedralChirality` | `Chirality.h:33-36`, `Chirality.cpp:45-54`, `860-863` | `RDKit✔️❌` | `var_os` materializes an owned value; C `getenv` returns a borrowed pointer |
| `bondAffectsAtomChirality` | `Chirality.cpp:878-889` | `RDKit✔️✔️` | CSR adjacency preserves insertion order and `O(degree)` traversal |
| `getAtomNonzeroDegree` | `Chirality.cpp:890-900` | `RDKit✔️✔️` | same adjacency complexity and no per-call allocation |
| `Point3D` normalize/direction, dot, angle, cross | `Geometry/point.h:147-194`, `210-234`; `Vector.h:24` | `RDKit✔️✔️` | same scalar operation order and constant-space arithmetic |
| macro-only `VOLTEST` | `Chirality.cpp:3114` | `RDKit✔️✔️` | direct dot/cross expression with source comparison semantics |
| `OctahedralPermFrom3D` | `Chirality.cpp:3116-3171` | `RDKit✔️✔️` | fixed switch table and constant-time predicates |
| `isWigglyBond` | `Chirality.cpp:3173-3184` | `RDKit✔️✔️` | direct bond/property lookup |
| `assignNontetrahedralChiralTypeFrom3D` | `Chirality.cpp:3188-3375` | `RDKit✔️✔️` | same bounded neighbor matrix and nested branch structure |
| tetrahedral signed-volume branch | `Chirality.cpp:3434-3498` | `RDKit✔️✔️` | same source-order incident-neighbor traversal and constant-size arithmetic |
| `assignChiralTypesFrom3D` | `Chirality.cpp:3377-3500` | `RDKit✔️❌` | transactional operation clones topology/property blocks before committing, adding `O(atoms + bonds + properties)` work |

All corresponding C++ lines are present verbatim beside the Rust behavior in
`crates/cosmolkit-core/src/chemistry/stereo.rs`. The former constant-true
`is_opposite_bonds` helper and fabricated 2D
`assign_atom_chiral_tags_from_structure` production function no longer exist.
The notation-level function with a similar name is now only a parameter adapter
to the canonical chemistry kernel.

## Focused Validation

Each source-focused unit command executed exactly one named test and passed
with `1 passed`, `0 failed`, and `0 ignored`:

```text
cargo test -p cosmolkit-core --features op-contracts-strict --lib chemistry::stereo::tests::get_allow_nontetrahedral_chirality_matches_rdkit_environment_semantics -- --exact
cargo test -p cosmolkit-core --features op-contracts-strict --lib chemistry::stereo::tests::bond_affects_atom_chirality_matches_rdkit_directional_matrix -- --exact
cargo test -p cosmolkit-core --features op-contracts-strict --lib chemistry::stereo::tests::atom_nonzero_degree_matches_rdkit_bond_matrix -- --exact
cargo test -p cosmolkit-core --features op-contracts-strict --lib chemistry::stereo::tests::is_wiggly_bond_matches_rdkit_atom_relative_matrix -- --exact
cargo test -p cosmolkit-core --features op-contracts-strict --lib chemistry::stereo::tests::rdkit_direction_vector_matches_source_boundaries -- --exact
cargo test -p cosmolkit-core --features op-contracts-strict --lib chemistry::stereo::tests::rdkit_angle_to_matches_source_boundaries -- --exact
cargo test -p cosmolkit-core --features op-contracts-strict --lib chemistry::stereo::tests::voltest_matches_rdkit_signed_zero_and_nonfinite_matrix -- --exact
cargo test -p cosmolkit-core --features op-contracts-strict --lib chemistry::stereo::tests::octahedral_perm_from_3d_matches_all_rdkit_switch_branches -- --exact
cargo test -p cosmolkit-core --features op-contracts-strict --lib chemistry::stereo::tests::assign_nontetrahedral_chiral_type_from_3d_covers_all_active_rdkit_branches -- --exact
cargo test -p cosmolkit-core --features op-contracts-strict --lib chemistry::stereo::tests::tetrahedral_chiral_type_from_3d_matches_rdkit_branch_boundaries -- --exact
cargo test -p cosmolkit-core --features op-contracts-strict --lib chemistry::stereo::tests::assign_chiral_types_from_3d_covers_complete_rdkit_control_flow -- --exact
```

The final dedicated boundary runs were:

| Command | Actual result |
|---|---|
| `cargo test -p cosmolkit-core --features op-contracts-strict --test assign_atom_chiral_tags_from_structure_operation -- --nocapture` | 4 passed; 0 failed; 0 ignored; 0 filtered |
| `cargo test -p cosmolkit-core --features op-contracts-strict --test assign_atom_chiral_tags_from_structure_parser_paths -- --nocapture` | 4 passed; 0 failed; 0 ignored; 0 filtered |
| `cargo test -p cosmolkit-core --features op-contracts-strict --lib notation::smiles::tests::assign_chiral_types_from_3d` | 13 passed; 0 failed; 0 ignored; 2593 filtered |
| `.venv/bin/pytest python/tests/test_assign_atom_chiral_tags_from_structure.py python/tests/test_stub_surface.py` | 8 passed; 0 failed; 0 skipped; 0 deselected |

No focused command was accepted with zero executed tests.

## Exact RDKit Oracle

The generated file
`testdata/stereo/expected/rdkit/smiles_small/assign_atom_chiral_tags_from_structure.jsonl`
contains exactly 77 records produced with pinned RDKit `2026.3.1`. The final
command was:

```bash
cargo test -p cosmolkit-core --release --features op-contracts-strict \
  --test rdkit_assign_atom_chiral_tags_from_structure_parity -- --nocapture
```

It executed `1 passed`, `0 failed`, `0 ignored`, and `0 filtered`. The test
asserts both fixture and oracle lengths are 77 before iterating, then compares
every row without tolerance, filtering, fallback, or field omission. All 73
successful rows and four error rows matched the pinned RDKit status and defined
observable fields described above. Re-running the public generator reported:

```text
[expected] reuse stereo/smiles_small: exact identity and output checksums match
Validated 2 RDKit expected outputs for profile smiles_small suite stereo
```

## Public API Mapping

| Layer | API | Mapping |
|---|---|---|
| RDKit Python | `Chem.AssignAtomChiralTagsFromStructure(mol, confId=-1, replaceExistingTags=True)` | wrapper calls `MolOps::assignChiralTypesFrom3D` |
| Rust value API | `Molecule::with_chiral_tags_from_structure(i32, bool)` | registered weak topology-state operation; returns a new molecule |
| Rust in-place API | `Molecule::assign_chiral_tags_from_structure_(i32, bool)` | same registered operation; transactional trailing-underscore mutation |
| Python value API | `Molecule.with_chiral_tags_from_structure(conf_id=-1, replace_existing_tags=True)` | calls the Rust value API and returns a new molecule |
| Python in-place API | `Molecule.assign_chiral_tags_from_structure_(conf_id=-1, replace_existing_tags=True)` | calls the Rust in-place API and returns `None` |

Both Python declarations are generated into `python/cosmolkit.pyi`; there is no
hand-written duplicate stub declaration. Missing conformers map to `ValueError`
and zero-length direction vectors map to `RuntimeError` with the fixed RDKit
oracle text.

## Parser Call Sites

There is one chemistry implementation. Format-specific functions preserve
their original parser sequencing but only adapt parameters and error types:

| Path | Canonical call | Sequencing status |
|---|---|---|
| SDF/V2000/V3000 postprocess | `assign_chiral_types_from_3d_molecule(&mut molecule, -1, true)` | retained at the source-defined post-parse position |
| MOL2 | same kernel with `-1, true` | retained before sanitization and hydrogen removal |
| PDB-to-Molecule compatibility path | same kernel with `-1, true` | retained before standard residue chirality handling |
| 3D CXSMILES | same kernel with selected 3D conformer id and `true` | retained only when no 2D conformer precedes the first 3D conformer |

The four-path parser integration target compares tag, permutation,
`_NonExplicit3DChirality`, coordinate, and conformer-dimensionality state and
passes all four tests. No parser contains an operation-local copy of the
assignment algorithm, and no parser silently discards a kernel error.

## Structured Error Differences

- RDKit directly mutates its `ROMol`; the
  `partial_mutation_before_later_exception` oracle row therefore observes
  partial RDKit working-state mutation before the zero-length-vector exception.
  COSMolKit deliberately computes on detached operation blocks and commits only
  on success. Both Rust public forms leave the caller molecule unchanged on
  error. This is the operation transaction policy, is explicitly tested, and
  does not alter successful chemistry results.
- Invalid coordinate-row alignment is a C++ precondition violation. COSMolKit
  exposes it as `StereoError::ConformerAtomCountMismatch` rather than allowing
  unchecked indexing or a panic.
- The public Rust operation requires modeled implicit-hydrogen cache state and
  reports `StereoError::MissingImplicitHydrogenState` if it is absent. Parser
  adapters construct that state before invoking the kernel. This is a
  structured model-boundary error, not a fallback or guessed valence result.

## Repository Validation

- `cargo fmt --all -- --check`: passed after the final source-comment and test
  cleanup.
- `cargo check -p cosmolkit-core --features op-contracts-strict`: passed.
- `cargo test -p cosmolkit-core --release --features op-contracts-strict`:
  core unit target reported 2561 passed, 0 failed, and 45 explicitly ignored;
  every integration and doc target completed successfully.
- `cargo test --workspace --release --features
  cosmolkit-core/op-contracts-strict`: 46 test/doc targets completed
  successfully. The InChI library target reported 1197 passed and 74 ignored;
  the ignored set is 43 official-C and 31 RDKit live-oracle tests that require
  pinned vendored source trees and must be run explicitly with `--ignored`.
  The provenance integration target executed its metadata test and retained
  one explicit live-C ignored test. None of those ignored tests is claimed as
  executed.
- `.venv/bin/pytest`: 538 passed, 0 failed, 37 skipped, and 0 deselected from
  575 collected tests.
- `.venv/bin/basedpyright python/tests python/examples`: 0 errors, 416 existing
  warnings, and 0 notes. The command exits nonzero under the repository warning
  policy; this result is not described as warning-free or zero-exit.

## Unclosed Active Branches

None within the selected `AssignAtomChiralTagsFromStructure` /
`assignChiralTypesFrom3D` modeled input state space. The following remain
explicitly outside this call graph and are not completion blockers for this
plan: `assignStereochemistryFrom3D`, double-bond E/Z assignment, CIP
assignment/filtering, atropisomer perception, and distinct-substituent
validation.

## Audit Conclusion

The selected active source closure, registered Rust operation, Rust and Python
public mappings, and four parser call sites are complete. Every selected source
branch has a verbatim frame and an honest final two-axis marker. The focused
tests execute nonzero cases, the fixed 77-row pinned-RDKit oracle passes exact
field comparison, expected output regeneration reuses identical checksums, and
there is no unclosed active branch in this plan's scope.
