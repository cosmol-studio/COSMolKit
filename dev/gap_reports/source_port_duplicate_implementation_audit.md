# Source Port Duplicate-Implementation Audit

## Scope

This audit checks source-backed production Rust under:

- `crates/cosmolkit-core/src/`
- `crates/cosmolkit-inchi/src/`
- `crates/cosmolkit-ringdecomposer/src/`

Tests, oracle programs, `third_party/`, `tmp/`, generated artifacts, and archived
plans are not treated as production implementations. Test code is used only as
evidence about whether a production path is exercised.

The audit distinguishes four cases:

1. the same upstream function has two independent Rust behavior bodies;
2. identical Rust logic or data was copied without a shared kernel;
3. one upstream function was intentionally split into coordinator, parser, or
   wrapper phases which delegate instead of reimplementing behavior;
4. functions have the same unqualified name but are different upstream
   definitions because their source files or static scopes differ.

Textual source-frame equality alone is not evidence of duplicate behavior.
The classification below was made by inspecting callers and Rust bodies, not by
counting function names alone.

## Inventory Result

| Area | Standard source frames | Duplicate upstream identities | Result |
|---|---:|---:|---|
| InChI | 1,157 | 0 | No duplicate `source file + function` identity found. |
| RingDecomposerLib | 36 | 0 | No duplicate source-function identity found. |
| Core RDKit/Gemmi ports | not normalized to one frame format | 34 repeated frame-header groups | Mixture of real duplication and legitimate decomposition; classified below. |

In InChI, `PrepareSaveOptBits` appears in both `ichiread.c` and `runichi.c`.
Those are distinct upstream definitions and are not duplicate ports. No evidence
was found that the stepwise InChI port implemented one exact upstream function
twice.

## Consolidation Boundary

Only the following groups are approved for consolidation by the current
same-source rule:

| Group | Upstream identity | Allowed shared portion |
|---|---|---|
| Swap counting | `countSwapsToInterconvert` | Complete generic swap-count kernel; caller error projection remains separate. |
| Ring angle | `DGeomHelpers::_setRingAngle` | Hybridization/ring-size calculation; molecule access remains separate. |
| Periodic table | `PeriodicTable::getAtomicNumber` and `getElementSymbol` | Canonical lookup tables; format parsing and error types remain separate. |
| Force-field conformers | `ForceFieldsHelper::OptimizeMoleculeConfs`, active non-threaded `getNumThreadsToUse`, and `OptimizeMoleculeConfsST` | Shared conformer iteration/minimization driver; MMFF/UFF construction and result types remain separate. |
| Fused-ring utilities | `RingUtils::makeRingNeighborMap` and `pickFusedRings` | Overlap map and traversal kernel; caller filtering remains separate. |
| Canon stereo helpers | `chiralAtomNeedsTagInversion`, `atomHasFourthValence`, `hasSingleHQuery`, and `getChiralPermutation` | Source-identical decisions parameterized by explicit read adapters. |
| Isotope map | `getIsoMap` | Source traversal; builder/molecule access and errors remain separate. |

The following candidates are frozen and must not be mechanically consolidated:

- the forward, indexed, and multithreaded SDF record readers come from
  different supplier functions with different delimiter and state behavior;
- atropisomer format paths contain different source functions and modeled input
  state despite sharing local vector operations;
- the two current `assignChiralTypesFromBondDirs` bodies require a branch-level
  source audit before any common kernel can be proven;
- `atom_nonzero_degree` currently has differing dative semantics and cannot be
  unified before resolving that difference against source;
- SGroup membership, bond-order conversion/ranking, and minstd RNG transition
  copies do not currently carry enough common source identity in every location
  to satisfy the same-source rule.

## Completed Same-Source Consolidations

The approved groups now have one behavior owner for the proven
source-identical portion. Caller adapters remain separate where ownership,
input representation, or error projection differs.

| Group | Canonical Rust owner | Retained adapters | Focused evidence |
|---|---|---|---|
| `countSwapsToInterconvert` | `crates/cosmolkit-core/src/source_port_helpers.rs::count_swaps_to_interconvert` | Canon ranking, substructure matching, remove-H stereo, and SMILES writing retain their source-defined panic/`None`/domain-error projection. | `shared_count_swaps`: 7 passed, 0 failed. |
| `DGeomHelpers::_setRingAngle` | `crates/cosmolkit-core/src/source_port_helpers.rs::rdkit_set_ring_angle` | ConfSeq and distance-geometry callers retain their molecule/ring access. | `shared_set_ring_angle`: 2 passed, 0 failed. |
| `PeriodicTable::getAtomicNumber` / `getElementSymbol` | `chemistry/valence.rs::rdkit_atomic_number_from_symbol` and `rdkit_element_symbol` | MOL2 retains parse errors; PDB and SMILES retain their source-defined invalid-symbol fallback. | `shared_periodic_table_delegation`: 3 passed, 0 failed. |
| `ForceFieldsHelper::OptimizeMoleculeConfs` active non-threaded path | `chemistry/forcefield/mod.rs::optimize_molecule_confs_non_threaded` | MMFF and UFF retain parameter construction, missing-parameter status, and result types. | `shared_forcefield_conformer_driver`: 9 passed, 0 failed. |
| `RingUtils::makeRingNeighborMap` / `pickFusedRings` | `source_port_helpers.rs::rdkit_make_ring_neighbor_map` and `rdkit_pick_fused_rings` | Aromaticity and kekulization retain caller-specific ring selection. | `shared_ring_utils`: 4 passed, 0 failed. |
| Canon stereo helpers | `notation/smiles.rs::{canon_chiral_atom_needs_tag_inversion, canon_atom_has_fourth_valence, atom_query_has_single_h_count, nontetrahedral_chiral_permutation_kernel}` | Builder and immutable-molecule access, cache lookup, and caller errors remain separate. | `shared_canon_helpers`: 9 passed, 0 failed. |
| `AddHs.cpp::getIsoMap` bond traversal | `source_port_helpers.rs::rdkit_get_iso_map` | SMILES builder mutation and remove-H operation assignment retain separate atom-property clearing and structured index errors. | `shared_get_iso_map`: 5 passed, 0 failed. |

The `getIsoMap` consolidation also removes a pre-existing ordering difference:
the remove-H copy accumulated keys by first bond encounter, while official
`std::map` returns keys in atom-index order. The shared kernel uses `BTreeMap`,
preserves isotope value order within each key, and is covered by a focused
ordering test.

## High-Risk Findings

### 1. `assignAtomChiralTagsFromStructure` has a divergent public implementation

The public implementation at
`crates/cosmolkit-core/src/chemistry/stereo.rs:3769` describes itself as
simplified, carries only an abbreviated source frame, and delegates an essential
geometry decision to `is_opposite_bonds()` at line 3830. That helper currently
returns `true` for every input. This is not merely duplicated utility code: it
is a separate, behaviorally incomplete implementation of the same upstream
capability.

The source-backed 3D implementation already developed for parsing lives in
`crates/cosmolkit-core/src/notation/smiles/stereo.rs` around line 2843. The SDF
path in `crates/cosmolkit-core/src/io/sdf/postprocess.rs` around line 1218 is a
wrapper around that implementation rather than a third geometry engine.

Required resolution:

- complete the existing source-level port plan for
  `assignAtomChiralTagsFromStructure`;
- establish one canonical geometry kernel with explicit input-state adapters;
- make the public API and SDF/SMILES paths call that kernel;
- remove the unconditional `is_opposite_bonds()` behavior only after focused
  source-branch and RDKit parity tests pass.

This is the highest-priority issue because a public API can currently return a
chemically meaningful-looking result from behavior explicitly documented as
simplified.

### 2. Atropisomer handling has multiple behavior centers

There are three independently evolved areas:

- public APIs in `crates/cosmolkit-core/src/chemistry/atropisomer.rs`;
- SMILES/conformer assignment in
  `crates/cosmolkit-core/src/notation/smiles/stereo.rs`;
- the most complete SDF source-backed implementation in
  `crates/cosmolkit-core/src/io/sdf/postprocess.rs:1367` through approximately
  line 2180.

The public `chemistry::atropisomer` module is re-exported by
`crates/cosmolkit-core/src/lib.rs:83`. Its
`assign_atropisomer_stereo()` at line 484 explicitly says that full RDKit
behavior is unfinished and returns `(BondId, ChiralTag)` even though the source
capability assigns bond stereo. Its detection path also has its own candidate
scan and structural filters. Searches found no production caller of these
public assignment/detection functions outside that module's tests.

Across the SMILES and SDF paths, local copies also exist for candidate
selection, neighbor-bond collection, opposite-atom lookup, wedge direction,
projected vectors, and vector arithmetic. The input adapters legitimately
differ, but the source-defined geometric and parity decisions should not.

Required resolution:

- choose the complete `BondStereo`-based source behavior as the canonical
  internal kernel;
- keep format-specific coordinate and error adapters at the boundary;
- make public APIs delegate to the canonical kernel or return structured
  unsupported errors for state that the canonical kernel cannot model;
- do not preserve the current `ChiralTag` projection as a compatibility
  implementation unless it is deliberately specified as a different API.

### 3. `assignChiralTypesFromBondDirs` is independently ported twice

Active implementations exist at:

- `crates/cosmolkit-core/src/notation/smiles/stereo.rs:2041`;
- `crates/cosmolkit-core/src/io/sdf/postprocess.rs:1075`.

They reproduce the same upstream operation but have independently accumulated
2D fallback, error, mutation, and cache behavior. The format boundary explains
different adapters, not two separate stereo assignment algorithms.

Required resolution: extract a source-backed assignment kernel whose inputs
explicitly represent the available conformer and replacement policy. Preserve
format-specific error conversion and output application in thin callers.

The related `atom_nonzero_degree()` helpers at
`chemistry/stereo.rs:3324` and `notation/smiles/stereo.rs:2818` are also not
semantically identical: one has atom-relative dative handling while the other
excludes dative orders broadly. Resolve this against the upstream definition
before deduplicating; selecting either current helper by convenience could
change behavior.

### 4. Periodic-table conversion data was copied into format code

Canonical conversion APIs exist at:

- `crates/cosmolkit-core/src/chemistry/valence.rs:2325`
  (`rdkit_atomic_number_from_symbol`);
- `crates/cosmolkit-core/src/chemistry/valence.rs:2623`
  (`rdkit_element_symbol`).

Independent tables were found in:

- `crates/cosmolkit-core/src/io/mol2.rs:288`, symbol to atomic number;
- `crates/cosmolkit-core/src/io/pdb_writer.rs:74`, atomic number to symbol;
- `crates/cosmolkit-core/src/notation/smiles_write.rs:1867`, atomic number to
  symbol.

This finding is resolved. MOL2, PDB, and SMILES now delegate to the canonical
helpers while retaining their format-specific invalid-input behavior. No
format API changed.

## Confirmed Duplicate Kernels

### 5. SDF forward and indexed raw-record loops are intentionally separate

`extract_next_forward_sdf_record()` at `io/sdf.rs:1096` and
`extract_next_indexed_sdf_record()` at `io/sdf.rs:1203` have line-for-line
equivalent Rust loops for reading through `$$$$`, normalizing line endings, and
constructing `RawSdfRecord`.

Their upstream source functions and supplier-level contracts differ, so the
current same-source policy does not permit extracting their textual overlap as
a shared source-port implementation. `extract_next_raw_sdf_record()` also has
different source-defined `M  END`, blank-line, and end-state conditions. All
three paths remain independent.

### 6. `countSwapsToInterconvert` has four active copies

The same swap-count algorithm appears in:

- `crates/cosmolkit-core/src/notation/canon_rank.rs:567`;
- `crates/cosmolkit-core/src/search/substruct.rs:1602`;
- `crates/cosmolkit-core/src/chemistry/hydrogens.rs:1949`;
- `crates/cosmolkit-core/src/notation/smiles_write/stereo.rs:1938`.

The copies differ mainly in element type and in whether invalid input panics,
returns `None`, or returns a domain error. A generic private kernel should
return a neutral result that distinguishes length mismatch and missing
elements; each caller should preserve its current source-defined error mapping.

Resolved by `source_port_helpers.rs::count_swaps_to_interconvert`; the four
caller projections remain separate.

### 7. MMFF and UFF copy the force-field conformer optimization driver

The following upstream helpers are independently embedded in both public paths:

- `ForceFieldsHelper::OptimizeMoleculeConfs`;
- the non-threaded `getNumThreadsToUse` branch;
- `ForceFieldsHelper::detail::OptimizeMoleculeConfsST`.

Locations:

- MMFF: `chemistry/forcefield/mmff/mol_properties.rs:396-454`;
- UFF: `chemistry/forcefield/uff/atom_typer.rs:798-866`.

Force-field construction and result types legitimately differ. Iteration over
conformers, force-field initialization, minimization, and status collection do
not. Introduce a shared internal driver parameterized by force-field setup and
result projection. Keep MMFF/UFF source wrappers and focused tests separate.

Resolved by `chemistry/forcefield/mod.rs::optimize_molecule_confs_non_threaded`.

### 8. Fused-ring utility algorithms are ported in aromaticity and kekulization

`RingUtils::makeRingNeighborMap` and `RingUtils::pickFusedRings` are implemented
in both:

- `chemistry/aromaticity.rs:533` and `:626`;
- `chemistry/kekulize.rs:1807` and `:1837`.

The caller-specific ring filtering differs, but overlap-map construction and
the source traversal algorithm are common. They should be generic over the ring
representation and selection predicate rather than copied. Tests must preserve
RDKit's traversal/order behavior, not just set equality.

Resolved by the two RingUtils kernels in `source_port_helpers.rs`; tests cover
duplicate-sensitive `Intersect` behavior and recursive traversal order.

### 9. Canonical SMILES helpers are copied across parser and writer paths

Repeated source functions include:

- `Canon::chiralAtomNeedsTagInversion` in `notation/smiles.rs:3109` and
  `notation/smiles_write/stereo.rs:1976`;
- `Canon::details::atomHasFourthValence` in `notation/smiles.rs:3135` and
  `notation/smiles_write/stereo.rs:1996`;
- `Canon::details::hasSingleHQuery` in `notation/smiles.rs:4760` and
  `notation/smiles_write/stereo.rs:2064`;
- `Chirality::getChiralPermutation` in `notation/smiles.rs:3193` and the
  molecule probe implementation at `notation/smiles.rs:4934`.

The parser operates on builder state while the writer operates on a molecule
and may consume derived valence state. In particular,
`atomHasFourthValence` cannot be replaced by either current body unchanged.
The source-defined pure decisions should be parameterized over a small read
view; builder/molecule cache access remains in adapters.

Resolved for the listed pure decisions. The builder and molecule adapters were
not merged because they read different modeled state and project different
errors.

### 10. Smaller copied production helpers

The following lower-risk candidates remain frozen because common upstream
function identity has not been established for every location:

- identical SGroup membership predicates in `chemistry/hydrogens.rs:1166` and
  `notation/smiles.rs:4685`;
- identical bond-order numeric conversion in `notation/smiles.rs:4732` and
  `notation/smiles_write/stereo.rs:2034`;
- the RDKit bond-order rank table in `notation/canon_rank.rs:2940` and
  `notation/smiles_write/stereo.rs:1166` (the return integer type differs);
- the minstd LCG state transition in `chemistry/coordinates.rs:144` and
  `chemistry/distgeom.rs:6573`.

The two formerly listed same-source items, `DGeomHelpers::_setRingAngle` and
`AddHs.cpp::getIsoMap`, are resolved in `source_port_helpers.rs`. For the RNG
case, seed handling and real-number distributions differ; even the common state
transition remains separate until both locations are proven to reproduce the
same upstream function rather than merely the same formula.

## Repeated Source Frames That Are Not Duplicate Behavior

The following repeated headers are intentional decomposition or delegation and
should not be removed as if they were duplicate ports:

| Repeated source function | Classification |
|---|---|
| `gemmi::make_structure` | The document adapter and structure builder are separate stages; the adapter delegates. |
| `Atom::calcImplicitValence` | Public assignment wrapper delegates to the canonical per-atom kernel. |
| `Bond::setStereo` | Internal setter invariant and fallible builder-boundary validation enforce the same contract at different boundaries. |
| `MolOps::addHs` | Assignment planning and topology application are deliberately split phases. |
| `details::KekulizeFragment` | SMILES writer calls the canonical chemistry assignment; it is not a second Kekulize engine. |
| `FindRings::removeExtraRings` | The sorting code is a called kernel used to reproduce source-library ordering, not a parallel `removeExtraRings`. |
| `PeriodicTable::getElementSymbol` in drawing | Drawing delegates to the canonical valence helper. Other copied tables remain a finding above. |
| UFF parameter helpers in `confseq` | Conformer-sequence wrappers delegate to canonical UFF parameter functions. |
| `GetMolFileBondStereoInfo` and `BondGetDirCode` | V2000/V3000 wrappers delegate to one molfile implementation. |
| `GetMolFileQueryInfo` | Value-query and atom-list output are different source branches, not duplicate bodies. |
| SDF `finishMolProcessing`, `MolFromMolDataStream`, `ParseV3000CTAB`, `ParseV3000BondBlock`, `ParseMolBlockProperties`, and `readMolProps` | Large upstream procedures were split into coordinator and focused parsing/application phases. No duplicate complete behavior body was found from the repeated header alone. |
| `ROMol::atomBonds` | Small adapters return different domain representations and error semantics; consolidation is optional and lower value. |

Source frames should remain complete in the canonical kernel. Thin wrappers may
retain the exact upstream call-site lines needed to explain the delegation, but
should not copy the full callee body after the behavior has one owner.

## Process Artifact

`dev/source_reproduction_protocol.md` contains the complete "Example D" section
twice, at lines 259 and 288. This does not affect production behavior, but it is
a stepwise-edit artifact in the normative protocol and should be reduced to one
copy in a separate documentation-only change.

## Recommended Repair Order

1. Complete and centralize `assignAtomChiralTagsFromStructure`; do not preserve
   the unconditional geometry result.
2. Define one source-backed atropisomer assignment/detection kernel and migrate
   format adapters and public APIs to it.
3. Unify `assignChiralTypesFromBondDirs` and resolve dative-degree semantics
   against the source.
4. Replace copied periodic-table data with the canonical conversion APIs.
5. Keep the different-source SDF supplier loops separate; do not consolidate
   them from textual similarity.
6. Keep the unresolved smaller helpers separate until exact upstream identity
   and active branches are established.

Each repair should be behavior-neutral unless the audit identifies an existing
incomplete/divergent implementation. For those divergent paths, source-level
focused tests and the pinned external parity oracle must establish the corrected
behavior before the old body is removed. No wrapper should gain a fallback to
the old implementation.

## Audit Conclusion

All seven groups approved by the same-source matrix have been consolidated.
Their focused filters executed 39 tests in total, all passing, with no failures.
This count excludes filtered-only test binaries.

No InChI or RingDecomposerLib production function was merged because the audit
found no duplicate exact upstream identity there. Remaining core repetition is
either tied to different upstream functions, intentionally split architecture,
or blocked on a source/branch audit. In particular, the public
`assignAtomChiralTagsFromStructure`, atropisomer paths, and the two
`assignChiralTypesFromBondDirs` bodies remain behavioral repair work; they were
not treated as safe deduplication targets in this pass.

The stepwise port did not duplicate any standardized InChI or RingDecomposerLib
source-function identity. It did leave meaningful duplication in the core,
including copied tables and helpers, two stereo bond-direction engines, several
atropisomer behavior centers, and one explicitly incomplete public 3D chirality
implementation alongside a more complete internal port.

The main risk is behavioral drift, not repository size. Consolidation should
therefore follow source-backed parity boundaries, not broad mechanical deletion.
