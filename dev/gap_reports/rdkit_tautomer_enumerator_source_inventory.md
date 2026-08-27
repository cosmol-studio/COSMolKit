# RDKit Tautomer Enumerator Source Inventory

Status: pre-port source inventory. This report fixes the selected upstream
call graph, test branches, COSMolKit prerequisite ownership, and target API
mapping before implementation begins. It is an audit baseline, not a support
claim.

## Pinned source and evidence boundary

- Oracle: vendored RDKit `2026.03.1`; the execution plan pins revision
  `351f8f378f8ad6bbd517980c38896e66bf907af8`.
- Core source: `Code/GraphMol/MolStandardize/Tautomer.cpp` and `Tautomer.h`.
- Catalog source: `TautomerCatalogParams.*`, `TautomerCatalogUtils.*`, the
  current 37-transform table, and the V1 36-transform table.
- Wrapper source: `Code/GraphMol/MolStandardize/Wrap/Tautomer.cpp`.
- Focused source tests: all 14 test functions in `testTautomer.cpp` and the
  tautomer-reachable sections of `catch_tests.cpp`.
- Canonical corpora: the vendored 1,000-row PCS file and the 99,991-row file
  named `100kPCS_tautomer.csv.gz`.
- Exact comparison includes ordered canonical isomeric SMILES, complete
  molecular state, status, modified atom/bond sets, scoring components,
  option branches, cancellation, custom scoring, and source immutability.
- Catalog-entry deserialization and the wider cleanup/parent pipeline are not
  reachable from the selected enumerator API and remain outside this port.

## Normative upstream call graph

1. `TautomerEnumerator` copies the cleanup options and constructs one catalog
   from a file or in-memory transform definitions.
2. `enumerate(const ROMol&)` produces canonical isomeric SMILES, updates the
   property cache non-strictly when required, computes SymmSSSR, and creates a
   canonical-kekulized working form.
3. The state machine visits a `std::map` in SMILES order, transforms in catalog
   order, and substructure matches in source order.
4. Each accepted match moves one hydrogen, sets endpoint `noImplicit`, edits
   bonds and charges, then executes exactly `KEKULIZE | SETAROMATICITY |
   SETCONJUGATION | SETHYBRIDIZATION | ADJUSTHS`.
5. Only `KekulizeException` rejects a product and continues. Other failures are
   not silently converted into a missing tautomer.
6. Stereo/isotopic-H handling is applied before SMILES deduplication and again
   when the global modified sets expand. Rekeying and duplicate merging remain
   ordered-map operations.
7. Canonical selection maximizes signed score, breaks ties by lexical canonical
   SMILES, copies the selected molecule, and performs clean forced stereo
   assignment.

## Function ledger

| Source function or type | Pinned location | Required COSMolKit representation | Pre-port state |
| --- | --- | --- | --- |
| `TautomerTransform` ctor/copy/assignment/destructor | `TautomerCatalogParams.h:34-62` | One owned value with one compiled query and ordered bond/charge edits | Missing |
| anonymous `getTautomer(name,...)` | `TautomerCatalogUtils.cpp:24-38` | Sole fields-to-transform constructor | Missing |
| anonymous `getTautomer(line)` | `TautomerCatalogUtils.cpp:40-86` | Sole line parser preserving tab/token/space behavior | Missing |
| `stringToBondType` | `TautomerCatalogUtils.cpp:91-110` | Bond-token parser | Missing |
| `stringToCharge` | `TautomerCatalogUtils.cpp:112-130` | Charge-token parser with structured invalid-symbol error | Missing |
| three `readTautomers` overloads | `TautomerCatalogUtils.cpp:132-181` | File, stream/text, and definition-list adapters into one parser | Missing |
| catalog constructors/copy/accessors | `TautomerCatalogParams.cpp:21-57` | `TautomerCatalog` value semantics and indexed access | Missing |
| catalog `toStream`/`Serialize` | `TautomerCatalogParams.cpp:59-67` | Deterministic source-equivalent count serialization | Missing |
| catalog `initFromStream`/`initFromString` | `TautomerCatalogParams.cpp:69-75` | Explicit upstream-under-construction error | Missing |
| catalog-entry stream methods | `TautomerCatalogEntry.cpp:21-67` | Reachability note only; transform serialization is absent upstream | Deliberate non-goal |
| `SubstructTerm` constructor | `Tautomer.cpp:75-81` | `TautomerScoreTerm` with compiled SMARTS | Missing |
| `getDefaultTautomerScoreSubstructs` | `Tautomer.cpp:83-99` | One static 12-term table | Missing |
| `scoreRings` | `Tautomer.cpp:33-72` | Ring component in `TautomerScore` | Missing |
| `scoreSubstructs` | `Tautomer.cpp:101-119` | Ordered term matching and signed accumulation | Missing |
| `scoreHeteroHs` | `Tautomer.cpp:121-130` | P/S/Se/Te total-H penalty | Missing |
| inline `scoreTautomer` | `Tautomer.h:89-91` | Aggregate score over the three components | Missing |
| private `Tautomer` record | `Tautomer.h:101-119` | Candidate molecule, kekulized form, counters, done flag | Missing |
| result constructors/iteration/accessors | `Tautomer.h:121-223` | One ordered `TautomerEnumeration` representation | Missing |
| callback interface | `Tautomer.h:225-232` | Borrowed Rust/Python callable per invocation | Missing |
| enumerator constructors/copy/options | `Tautomer.h:235-344`, `Tautomer.cpp:133-147` | `TautomerOptions` plus clonable `TautomerEnumerator` | Missing |
| `setTautomerStereoAndIsoHs` | `Tautomer.cpp:149-297` | One private source transition using shared typed state | Missing |
| deprecated vector `enumerate` | `Tautomer.cpp:299-310` | Private adapter over rich result | Missing |
| rich `enumerate` | `Tautomer.cpp:312-597` | Sole state machine through multi-molecule ops | Missing |
| result-aware `pickCanonical` | `Tautomer.cpp:600-633` | Fast path using retained SMILES keys | Missing |
| iterable `pickCanonical` | `Tautomer.h:374-417` | Public iterable path sharing scoring/tie-break logic | Missing |
| `canonicalize` | `Tautomer.cpp:635-646` | Value-returning delegate with reassign-stereo disabled during enumeration | Missing |
| `canonicalizeInPlace` | `Tautomer.cpp:648-688` | Private correspondence test only | Missing |
| params and V1 factories | `Tautomer.h:452-466` | `from_options()` and `v1()` delegates | Missing |

The wrapper helpers at `Wrap/Tautomer.cpp:23-274` cover result sequence
materialization, negative indexing, callback validation, constructors,
iterable dispatch, custom score callables, enumeration, and the default score
table. The bindings at lines 289-560 expose status, result properties,
enumerator options, scoring, V1 construction, and callbacks. COSMolKit must
provide idiomatic snake-case surfaces over the same Rust core rather than
copying this CamelCase API.

## Focused upstream branch inventory

### `testTautomer.cpp`

| Function and source range | Observable branches retained for parity |
| --- | --- |
| `testEnumerator`, 23-384 | ordinary enumeration, charged/aromatic/long-path examples, lexical output, completed and max-transform statuses |
| `testEnumeratorParams`, 385-815 | tautomer/transform caps, SP3 stereo, bond stereo in legacy/current modes, modified bonds, isotopic-H remove/retain paths |
| `testEnumeratorCallback`, 816-970 | cancellation timing/order, completed-or-canceled status, callback preservation through copy |
| `testCanonicalize`, 971-995 | 37-transform default, value and private in-place correspondence |
| `testPickCanonical`, 996-1019 | result-aware canonical choice over the canonical dataset |
| `testCustomScoreFunc`, 1020-1088 | canonicalize/result/vector/set custom scoring equivalence |
| `testEnumerationProblems`, 1089-1122 | known enumeration regressions and invalid product rejection |
| `testPickCanonical2`, 1123-1155 | score ordering and lexical tie behavior |
| `testEnumerateDetails`, 1156-1210 | exact modified atom/bond sets, map, SMILES and result access |
| `testGithub2990`, 1211-1285 | E/Z removal, `STEREOANY`, ring/non-ring handling and option preservation |
| `testPickCanonicalCIPChangeOnChiralCenter`, 1286-1522 | canonical selection with CIP/tag changes under every stereo option combination |
| `testTautomerEnumeratorResult_const_iterator`, 1523-1611 | bidirectional iteration, random indexing, ordered map consistency |
| `testGithub3430`, 1612-1636 | scoring regression for hydrate/aminal/hemiaminal-like forms |
| `testGithub3755`, 1637-1676 | canonical and private in-place equivalence on retained endpoints |

### `catch_tests.cpp` and wrapper branches

- `provide tautomer parameters as JSON`, lines 831-878: custom data with
  omitted fields, bond tokens, charge tokens, and canonical helper behavior.
- GitHub 5008, lines 959-981: phosphorus non-transforming, transforming, and
  canonical branches.
- `asymmetric imine tautomer generation`, lines 1002-1021: five asymmetric
  transform-count cases.
- GitHub 5318 tautomer section, lines 1067-1073: unsanitized input plus private
  in-place correspondence.
- GitHub 5402, lines 1096-1141: atom-order independence with default and custom
  transform ordering.
- GitHub 5784, lines 1143-1154: three recoverable kekulization failures.
- `Custom Scoring Functions`, lines 1736-1757: all score components, 12-term
  table size, and a replacement term list.
- Tautomer-parent and bulk parent tests exercise a wider MolStandardize
  pipeline and are recorded as downstream composition evidence, not as an
  authorization to port cleanup/fragment/reionization in this plan.
- Wrapper-only branches retained: negative indices, out-of-range `IndexError`,
  non-iterable input, callback `None`, invalid callback type/method, custom
  scorer conversion, result-aware vs iterable dispatch, and copy construction.

## Existing COSMolKit prerequisite audit

| Prerequisite | Existing source-backed implementation | Gap that must be closed on the tautomer mainline |
| --- | --- | --- |
| Multiple-molecule operation contract | `MultiMoleculeOpParts`, typed `result_type`/`assemble_fn`, per-branch validation and ordered emission | Add the tautomer registry entry and use it as the only mutation route |
| SMARTS parsing | Lists, negation, recursion, `z/Z`, `R/r/x`, ranges, valence, degree and bond queries are represented | Prove all 73 built-in transform SMARTS compile; fix shared parser only on demonstrated divergence |
| Substructure matching | Structured matcher and default parameter path exist | Add a narrow `MoleculeReadParts` entry path; operation bodies may not recover raw `Molecule`; current public convenience methods also erase matcher errors |
| Canonical isomeric SMILES | Source-backed writer/canonical ranking exists | Add a narrow read-parts writer path so candidates can be keyed inside operation capability boundaries |
| Kekulization | Canonical assignment and read-parts helpers exist | Reuse exact `clearAromaticFlags=false, canonical=true`; validate tautomer-specific failure filtering |
| SymmSSSR/ring state | `symmetrize_sssr_from_read_parts` and typed `RingInfo` exist | Reuse and validate cache lifecycle per emitted branch |
| Partial sanitization | All five required flags exist in the shared sanitize pipeline | Expose a narrow sibling operation helper under the tautomer spec; do not call a nested public operation or duplicate sanitization |
| Aromaticity/conjugation/hybridization/valence | Shared source-backed assignments exist | No local substitutes; any mismatch is fixed in the shared owner |
| Stereo and CIP | Typed atom/bond stereo, tracked isotopic H, CIP properties and assignment implementations exist | Tautomer-specific transition is missing; some shared stereo markers remain behaviorally provisional and require focused proof before being treated as exact prerequisites |
| Total hydrogen semantics | `valence::total_num_hydrogens()` is the single shared RDKit `Atom::getTotalNumHs()` reproduction, including `noImplicit`, cache precondition, and optional H-neighbor counting | Tautomer scoring calls this shared owner; no tautomer-local approximation remains |
| Atom/bond local edits | `OpParts::with_topology_mut` and cache/effect recording exist | Implement source-ordered endpoint, bond and charge edits within one branch lifecycle |
| Errors | `OperationError` has structured chemistry/sanitize/unsupported variants | Add tautomer/catalog-specific structured context without string-only error erasure |
| Rust/Python public surface | Molecule value semantics, PyO3, stub generator, docs patterns exist | Entire tautomer surface is missing |
| Evidence tooling | Pinned RDKit, PCS files, 5,000 corpus, ChEMBL 37 pipeline patterns exist | Tautomer generators, manifests, schemas and aggregate runner are missing |

The critical architecture gaps are adapters, not permission to fork shared
chemistry. Canonical SMILES and matching currently require whole `Molecule`
references, while a registered multi-output body may read candidates only via
`MoleculeReadParts`. Those components must gain narrow capability-compatible
entry points backed by the same implementation.

## Tautomer CIP transition audit

The pinned `setTautomerStereoAndIsoHs()` transition cannot truthfully use any
existing `CipStatePolicy`:

| Existing policy | Why it does not describe the source transition |
| --- | --- |
| `Preserve` | Modified atoms may clear or restore chiral tags and `_CIPCode`; modified bonds may clear/restore stereo atoms and adjacent bond directions; forced stereo assignment may replace CIP state. |
| `ClearComputed` | RDKit does not perform a uniform `clearComputedProps()` transition. It selectively preserves, restores, clears, and then optionally recomputes state according to the modified atom/bond sets and six enumerator options. |
| `Assign` | This policy is intentionally reserved for the modern `with_cip_labels_with_options` producer. Tautomer enumeration calls legacy `assignStereochemistry(cleanIt=true, force=true)` only when `reassignStereo` is enabled and also has a source branch that marks `_StereochemDone` without assignment. |

The operation therefore requires one narrowly named
`TautomerSourceTransition` policy. It delegates only validation of this
source-defined CIP/stereo lifecycle. It does not grant block access, raw
storage access, a new derived-cache effect, or authority to call ordinary
mutation APIs. Compile-time macro validation and strict runtime validation
must both reject that policy unless the registry method is exactly
`enumerate_tautomers_with_options`, and a registry test must prove it is the
only operation carrying the policy.

The registered operation requires topology and molecule-properties read/write
capability because the source transition touches atom/bond chiral, CIP, stereo,
direction, tracked isotopic-H, and `_StereochemDone` state. Derived ring/stereo
state remains governed separately by declared access and derived effects. Each
candidate must still finish its own `OpParts` lifecycle and molecule invariant
validation before `MultiMoleculeOpParts::emit()` can retain it.

Correctness remains an operation-specific obligation: focused strict tests
must cover selective `_CIPCode` clear/restore, computed CIP cleanup or
reassignment, the `reassignStereo=false` marker branch, SP2/SP3 option branches,
modified/non-modified atoms and bonds, and rejection of the policy on every
non-tautomer spec. The policy is not a general-purpose escape hatch.

## Catalog serialization reachability boundary

`TautomerCatalogParams::toStream()` and `Serialize()` are reachable catalog
APIs, but the pinned source writes only the decimal transform count followed by
one newline. It does not serialize transform definitions. COSMolKit therefore
reproduces that count-only output exactly and does not present it as a
round-trippable catalog format.

`TautomerCatalogParams::initFromStream()` and `initFromString()` both consist
solely of `UNDER_CONSTRUCTION("not implemented")` in the pinned source. The
Rust catalog returns `TautomerCatalogError::DeserializationUnderConstruction`
for the corresponding boundary; it does not infer transforms from the stored
count or invent a private serialization format.

`TautomerCatalogEntry` is not constructed, queried, or serialized anywhere in
the enumerator call graph. Its `toStream()` deliberately leaves the transform
pickling call commented out and stores only a bit id plus description. Its
`initFromStream()` constructs a transform with a null molecule and empty edit
vectors while the reaction-unpickling call is also commented out. Consequently
this hierarchy entry cannot preserve a tautomer transform and is excluded from
the production Rust API. Porting its bit-id/description shell would create an
apparently usable but chemically empty duplicate representation, so it remains
a written upstream-under-construction boundary rather than an implementation
gap in the enumerator.

## Target ownership and API mapping

- `chemistry::tautomer` owns catalog values, scoring values, candidate metadata,
  the exact source state machine, and public enumerator/result types.
- `operations::ops` owns only the registered multi-output wrapper, branch
  mutation application, cache/effect recording, and rich-result assembly.
- Existing SMARTS, substructure, canonical SMILES, rings, kekulization,
  sanitization, stereo and CIP modules remain the sole shared implementations.
- Rust public types are `TautomerTransform`, `TautomerCatalog`,
  `TautomerOptions`, `TautomerEnumerator`, `TautomerEnumeration`,
  `TautomerEnumerationStatus`, `TautomerScoreTerm`, and `TautomerScore`.
- `Molecule::enumerate_tautomers()` and
  `Molecule::enumerate_tautomers_with_options(...)` are value-style registered
  entry points. `canonicalize()` returns a new molecule; no duplicate public
  in-place canonicalizer is introduced.
- Python exposes the same ownership and behavior with snake-case methods,
  integer indices for modified sets, negative sequence indexing, structured
  exceptions, and no access to internal branch handles.

## Pre-port conclusion

No tautomer enumerator implementation currently exists, so there is no
historical algorithm to preserve or reconcile. The operation system now has
the required multi-molecule shape, and most chemistry prerequisites exist as
single shared implementations. The confirmed prerequisite work is limited to
narrow operation-compatible entry points for matching and canonical SMILES, a
shared exact total-H helper, the tautomer-specific stereo/isotopic-H lifecycle,
and registry/effect integration. Every chemistry discrepancy found during the
port must be repaired at its shared source owner and locked with a local
regression before tautomer evidence is accepted.
