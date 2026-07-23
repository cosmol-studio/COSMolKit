# cosmolkit-inchi

`cosmolkit-inchi` is the independent, toolkit-neutral Rust boundary for the
source-level port of the official IUPAC InChI engine and the narrow RDKit
molecule-adapter behavior required by COSMolKit parity.

This document fixes the architecture for the port. Described modules and
types are a required destination, not a claim that they are already
implemented. No production generation, parsing, conversion, AuxInfo, or
InChIKey API may be made public until its selected source call graph is fully
ported and tested.

## Source Alignment

The official engine source is pinned to:

- upstream project: IUPAC InChI;
- upstream repository: <https://github.com/IUPAC-InChI/InChI>;
- release/tag: `v1.07.5`, released `2026-02-17`;
- annotated tag object: `fa28fcb6fbba554c952ad321b5d7797a5fb001a5`;
- peeled commit: `11a87982bb518f57ac013f0b258c283655e1ea1d`;
- commit tree: `214d60b92a94b03508a6b262590262bb4f89e275`;
- source archive SHA-256:
  `88532b3f599d125940e91af5d3135b31b0392b4c5a6e3f25e6418d9b56c5d5e3`;
- aggregate checked-out source-tree SHA-256 (excluding submodule `.git`
  metadata):
  `4f903503c370a6dd21de0a2e6beb77e076a1cd6063b34b6bde4b57d91f9c4edd`;
- local source path: `third_party/InChI`;
- repository integration: pinned Git submodule;
- upstream license: MIT, preserved at `third_party/InChI/LICENSE`;
- approved production target:
  `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/CMakeLists.txt`.

The human author approved `v1.07.5` and MIT source redistribution on
`2026-07-19`. The complete 643-file upstream release tree is checked out by
the `third_party/InChI` submodule at the peeled commit above. Initialize it
with `git submodule update --init third_party/InChI`. The production port
covers the complete official `libinchi` target: 60 C files, 1,590 C-source
function definitions, and 13 production header-defined functions. CLI and
demo sources remain reference and fixture material, not production modules.

The independently traceable RDKit adapter source is pinned to the RDKit
submodule's `2026.03.1` baseline:

- upstream repository: <https://github.com/rdkit/rdkit>;
- release: `Release_2026.03.1`;
- local source path: `third_party/rdkit/External/INCHI-API`;
- adapter sources: `inchi.cpp` and `inchi.h`;
- submodule license: `third_party/rdkit/license.txt`.

The official engine source and RDKit adapter are separate authorities. The
adapter cannot substitute for the engine, and their copied source frames and
markers must never be combined.

## Crate Boundary

The crate owns toolkit-neutral atom, bond, coordinate, isotope, radical,
hydrogen, stereo, polymer, V3000, status, diagnostics, option, AuxInfo,
InChIKey, engine, and adapter data and behavior. It must remain independent of:

- `cosmolkit-core` and COSMolKit `Molecule`;
- the molecule operation registry and molecule caches;
- RDKit runtime libraries or binaries;
- Python and PyO3;
- Open Babel;
- native or system InChI libraries in production;
- subprocesses, external executables, and network access;
- SMILES or MolBlock regeneration used as an InChI substitute;
- `crates/cosmolkit-core-old`.

The conversion between COSMolKit `Molecule` and this crate's neutral graph
belongs in `cosmolkit-core`. Official C may be compiled only by independent,
non-production oracle tooling under `tests/`; it is never linked, loaded, or
invoked by the production crate.

## Module Map

### Facade modules

The eventual reusable facade is split by owned behavior:

| Module | Ownership |
|---|---|
| `model` | Validated neutral atoms, bonds, coordinates, isotopes, radicals, hydrogens, 0D stereo, polymers, and V3000 collections |
| `options` | Parsed engine and adapter options, standard/non-standard mode, timeouts, and exact option diagnostics |
| `output` | InChI text, AuxInfo, generation status, warning/error/log bytes, and reconstructed neutral structures |
| `key` | InChIKey input, output, XHash fields, and exact source return status |
| `status` | Structured engine, parser, adapter, allocation, and unsupported errors without erased source status codes |
| `engine` | Toolkit-neutral generation, validation, InChI-to-InChI, structure reconstruction, and AuxInfo entry points |
| `adapter` | Source-derived RDKit option normalization, cleanup ordering, graph mapping, stereo mapping, and diagnostics over neutral data |

Until the complete source call graph and its tests are complete, these modules
and their entry points remain `pub(crate)` or absent. A partial source port is
not exposed as an experimental production API.

### Shared internal modules

| Module | Responsibility |
|---|---|
| `source_types` | Layout-independent owned equivalents of every production C typedef, struct, union, enum, constant, callback shape, and configured integer semantic |
| `memory` | C allocation/free/reallocation behavior translated into explicit Rust ownership and fallible allocation |
| `bytes` | NUL-terminated strings, fixed C character arrays, bounded byte fields, cursors, and source-compatible byte classification |
| `source` | One private Rust module per production C translation unit, preserving source names and function ownership |

Types shared by source modules live in `source_types`; they are not duplicated
to break dependency cycles. General helpers live in `memory` or `bytes` only
when they model a shared source primitive. Chemistry logic remains in the
module corresponding to the C translation unit that defines it.

### Official source modules

Every production translation unit has exactly one private Rust module. Module
names retain the C basename so inventory entries, source frames, tests, and
call-graph edges remain mechanically traceable.

`source::api` maps the 12 `INCHI_API/libinchi/src` production files:

| Rust module | Official C file |
|---|---|
| `source::api::ichilnct` | `ichilnct.c` |
| `source::api::inchi_dll` | `inchi_dll.c` |
| `source::api::inchi_dll_a` | `inchi_dll_a.c` |
| `source::api::inchi_dll_a2` | `inchi_dll_a2.c` |
| `source::api::inchi_dll_b` | `inchi_dll_b.c` |
| `source::api::inchi_dll_main` | `inchi_dll_main.c` |
| `source::api::ixa::ixa_builder` | `ixa/ixa_builder.c` |
| `source::api::ixa::ixa_inchikey_builder` | `ixa/ixa_inchikey_builder.c` |
| `source::api::ixa::ixa_mol` | `ixa/ixa_mol.c` |
| `source::api::ixa::ixa_read_inchi` | `ixa/ixa_read_inchi.c` |
| `source::api::ixa::ixa_read_mol` | `ixa/ixa_read_mol.c` |
| `source::api::ixa::ixa_status` | `ixa/ixa_status.c` |

`source::base` maps the 48 `INCHI_BASE/src` production files:

| Rust modules | Official C files |
|---|---|
| `bcf_s`, `ichi_bns`, `ichi_io` | `bcf_s.c`, `ichi_bns.c`, `ichi_io.c` |
| `ichican2`, `ichicano`, `ichicans` | `ichican2.c`, `ichicano.c`, `ichicans.c` |
| `ichierr`, `ichiisot` | `ichierr.c`, `ichiisot.c` |
| `ichimak2`, `ichimake` | `ichimak2.c`, `ichimake.c` |
| `ichimap1`, `ichimap2`, `ichimap4` | `ichimap1.c`, `ichimap2.c`, `ichimap4.c` |
| `ichinorm`, `ichiparm` | `ichinorm.c`, `ichiparm.c` |
| `ichiprt1`, `ichiprt2`, `ichiprt3` | `ichiprt1.c`, `ichiprt2.c`, `ichiprt3.c` |
| `ichiqueu`, `ichiread`, `ichiring` | `ichiqueu.c`, `ichiread.c`, `ichiring.c` |
| `ichirvr1`, `ichirvr2`, `ichirvr3`, `ichirvr4` | `ichirvr1.c`, `ichirvr2.c`, `ichirvr3.c`, `ichirvr4.c` |
| `ichirvr5`, `ichirvr6`, `ichirvr7` | `ichirvr5.c`, `ichirvr6.c`, `ichirvr7.c` |
| `ichisort`, `ichister`, `ichitaut` | `ichisort.c`, `ichister.c`, `ichitaut.c` |
| `ikey_base26`, `ikey_dll` | `ikey_base26.c`, `ikey_dll.c` |
| `inchi_gui`, `mol2atom` | `inchi_gui.c`, `mol2atom.c` |
| `mol_fmt1`, `mol_fmt2`, `mol_fmt3`, `mol_fmt4` | `mol_fmt1.c`, `mol_fmt2.c`, `mol_fmt3.c`, `mol_fmt4.c` |
| `permutation_util`, `readinch` | `permutation_util.c`, `readinch.c` |
| `runichi`, `runichi2`, `runichi3`, `runichi4` | `runichi.c`, `runichi2.c`, `runichi3.c`, `runichi4.c` |
| `sha2`, `strutil`, `util` | `sha2.c`, `strutil.c`, `util.c` |

Header-defined production functions are placed in the private module matching
the defining header basename under `source::base` unless that header is owned
by the API tree. Their source frame records the header path and line exactly;
they are never silently rewritten as unrelated generic helpers.

The RDKit translation is separate under `adapter::rdkit_2026_03` with source
ownership split into `inchi_cpp` and `inchi_h`. It may call completed engine
facade behavior but may not move official engine functions into adapter
modules.

## Owned Type Policy

### Scalar semantics

The source layer uses explicit-width Rust types and never substitutes
`usize`/`isize` for source arithmetic. The audited GCC LP64 production profile
is represented as follows:

| C semantic | Rust source semantic |
|---|---|
| `signed char` / `S_CHAR` | `i8` |
| `unsigned char` / `U_CHAR` | `u8` |
| `signed short` / `S_SHORT` | `i16` |
| `unsigned short` / `U_SHORT` | `u16` |
| `int` / `unsigned int` | `i32` / `u32` |
| `long` / `unsigned long` | `i64` / `u64` in the pinned LP64 profile |
| `AT_NUM` | `i16` |
| `AT_NUMB`, `AT_NUMBR`, `AT_RANK` | `u16` |
| `NUM_H`, `NUM_HS`, `ST_CAP_FLOW` | `i16` |
| `INCHI_MODE`, `INCHI_MODES` | explicit source-mode bit type with `u64` LP64 storage |
| `float`, `double` | `f32`, `f64` with tested IEEE-754 assumptions |

Platform-conditioned code must be translated from the selected GCC/Linux
preprocessor branch recorded by the inventory. A function may not inherit
host-dependent Rust widths. Integer casts, sign extension, truncation, masks,
and source-permitted unsigned wrap use explicit operations. Signed overflow,
out-of-range conversion, NaN behavior, and shifts are reproduced only from
source evidence and focused tests; they are never guessed.

C enums and flag sets have raw source representations that preserve every
integer value accepted by the C call graph, including unknown or invalid
values needed for exact error paths. Validated facade enums are separate
types. Converting a raw value to a facade enum is fallible and cannot collapse
unknown values into a default.

### Text and byte containers

C text is bytes first:

- `char *`, ASCIIZ fields, Molfile text, option strings, InChI text, AuxInfo,
  logs, and error buffers use owned `Vec<u8>`, `Box<[u8]>`, fixed `[u8; N]`,
  or an explicit NUL-terminated byte wrapper in the source layer;
- fixed arrays preserve their exact source length, zero fill, truncation, and
  terminator behavior;
- `MOL_COORD` is an exact fixed byte field, not a floating coordinate tuple;
- cursors and slices carry explicit bounds and current offsets;
- UTF-8 `String` is used only after a source-backed validation boundary says
  the bytes are valid text. Lossy conversion is not production behavior.

The neutral facade may expose validated UTF-8 where the source contract makes
that valid, while retaining raw diagnostic bytes when the source can emit
non-UTF-8 data.

### Owned graph and aggregate containers

Production structs, unions, callbacks, and arrays are translated into owned,
layout-independent Rust values:

- fixed C arrays become `[T; N]` when `N` is a source constant;
- pointer-plus-count arrays become `Vec<T>` or `Box<[T]>`, retaining a
  separate logical length/capacity field when the C algorithm observes it;
- nullable unique allocations become `Option<Box<T>>`, `Option<Vec<T>>`, or
  an equivalent owned wrapper;
- pointer identity, interior pointers, aliasing, and pointer subtraction
  become typed indices, IDs, ranges, or arena handles with the same ordering
  and lifetime semantics;
- pointer-to-pointer mutation becomes explicit ownership transfer or a typed
  mutable slot, not shared raw mutable storage;
- tagged unions become Rust enums only when the source has a complete active
  tag. Untagged or partially tagged unions use a dedicated representation
  that preserves every source-observable bit/value state;
- callback typedefs become exact internal Rust callback traits or function
  shapes, including nullability and return status. They are not `extern "C"`
  declarations;
- source global work state becomes engine-owned context passed explicitly.
  There is no unsynchronized mutable global state.

Atom and bond ordering, adjacency order, component order, queue/stack order,
canonical ranks, stereo neighbor order, and iteration tie breaks are stored
explicitly. Hash-based containers may not replace ordered source structures
when iteration is observable. Convenience cloning, sorting, deduplication, or
normalization is forbidden unless the copied source performs the same action.

Public values own their data. No public borrow points into temporary engine
state, and no native allocation crosses the crate boundary.

## Allocation Policy

The production crate uses the Rust allocator only. It contains no libc
allocator calls, FFI, dynamic loading, system InChI handle, or external
process. Source allocation functions remain individually ported functions
because their overflow, zeroing, null/no-op, ownership, and failure behavior
are part of the call graph.

The translation rules are:

1. `malloc`-style source allocation creates an owned buffer with the exact
   requested logical size. Bytes are not assumed initialized unless the source
   initializes them before reading.
2. `calloc` checks count-times-size overflow and creates exactly zeroed source
   values.
3. `realloc` preserves the old allocation and contents on failure, preserves
   the source-defined prefix on success, and implements the source's zero-size
   case exactly.
4. `free(NULL)` and source free helpers are explicit no-ops; freeing a present
   allocation consumes or takes its Rust owner exactly once.
5. Fallible growth uses `try_reserve` or an equivalent checked constructor.
   Allocation failure is carried as a structured internal allocation status
   and mapped to the exact official return code and diagnostics at the same
   boundary as C. It is not converted to an empty result or ignored.
6. Size multiplication, addition, alignment, and index conversion are checked
   before allocation. Source-defined unsigned wrapping is used only where the
   source actually relies on it.
7. Capacity and allocation count follow the C algorithm where they affect
   pointer/index validity or complexity. Rust containers may use a different
   capacity only after local complexity and allocation review proves no
   observable behavior or material performance regression.
8. `MaybeUninit` or `unsafe` is allowed only when exact source behavior cannot
   be represented reasonably with initialized safe values. Any use must be
   localized, document its soundness invariants, and have focused tests.

Rust drop replaces cleanup only where it is behaviorally equivalent. Named C
free functions still exist as private source-port boundaries when inventoried,
with complete source frames and focused tests. Ownership wrappers must permit
the same partial-construction cleanup paths without double-free, leaks, or
stale aliases.

## Source-Frame Policy

Every Rust function corresponding to official production source contains the
complete original C definition inside that Rust function body. The frame
includes the defining path, source line, exact C signature/body text, and
two-axis marker on every copied line:

```rust
fn example(/* owned Rust arguments */) {
    // BEGIN INCHI C FUNCTION: INCHI_BASE/src/example.c:123 example
    // INCHI❌❌: int example(int value)
    // INCHI❌❌: {
    // INCHI❌❌:     return value + 1;
    // INCHI❌❌: }
    // END INCHI C FUNCTION: example

    // Rust implementation follows the frame.
}
```

The text after `INCHI<behavior><complexity>:` is copied verbatim. Frames are
permanent review anchors and are not moved to appendices, generated files, or
unrelated functions. The complete definition runs from its declaration through
its closing brace. Relevant macro-controlled and implicit C behavior is copied
or anchored beside the operation that reproduces it.

Official engine frames use only the `INCHI` prefix. RDKit adapter frames use
`BEGIN RDKIT CXX FUNCTION` / `END RDKIT CXX FUNCTION` and only the `RDKit`
prefix. The official MIT notice and RDKit license attribution remain
traceable; source must not be relabeled between projects.

Marker axes are independent:

- the first symbol records behavioral reproduction;
- the second records locally reviewed performance and algorithmic complexity;
- `❌❌` remains until behavior exists;
- a behavioral `✔️` requires implementation and source-derived tests for the
  currently modeled source state space;
- a complexity `✔️` requires local review of asymptotics, loop nesting,
  scanning, allocation, temporary values, cloning, container choice, lookup,
  buffering, cache locality, and hot-path branching;
- `❗` is not a default for unbenchmarked code;
- `🔝` requires a nearby explanation of the behavior-preserving improvement.

Copying a frame, adding a dispatch stub, declaring FFI, returning a
placeholder, or calling an unported fallback never completes a port. If a
function depends on unmodeled state, it remains `INCHI❌❌` and returns a
structured unsupported result where reachable; it may not guess.

Cross-file helper behavior is framed from the source file that actually
defines the helper. For audited mutual recursion, the first port may contain
complete source-framed operation-local companion implementations exactly as
specified by the generated plan; the later companion step promotes them and
removes the local copies. No other duplication is used to bypass callee-first
ordering.

Each function receives its uniquely named source-derived focused test from
`dev/rdkit_inchi_full_port_plan.md`. A marker is upgraded only after that
function's implementation, edge behavior, allocation shape, and focused test
have been inspected. Aggregate parity suites cannot replace these tests.

## Unsupported Behavior

Before its source call graph is complete, behavior stays private or returns a
structured unsupported error. The crate must never produce chemically
meaningful-looking output through a formula, SMILES, MolBlock regeneration,
similarity rule, external command, native library, Open Babel call, or silent
fallback. Source ambiguity is resolved from pinned source and oracle fixtures,
not by heuristic approximation.
