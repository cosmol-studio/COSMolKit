# Native Molecule Archive 2.0

Status: implemented by `src/molecule_binary/archive_v2.rs`. The filename is
retained from the approved draft. Archive 2.0 is a CK-native format, not RDKit
pickle compatibility; its version is independent of the CK package version.

This is the first V2 definition, using Müsli 0.1.9. There are no earlier
supported V2 schemas or codecs: no V2-to-V2 compatibility branch, migration,
codec ID, or schema version is introduced for this update.
Legacy raw and archive 1.x inputs still use their existing reader branches.

## 1. Goals and public API

Keep `Molecule::to_binary()` and `Molecule::from_binary()` unchanged. The writer
emits archive 2.0; the reader automatically selects the new or legacy format
from the file header. Users do not rename APIs, pass a format version, or
manually convert their existing files.

Use `musli::storage` for block payloads. CK owns the envelope, schema, field
identities, migration rules, resource limits, and validation. Müsli is the
encoding implementation, not the compatibility policy.

Store complete molecule state once. Do not carry forward the proposed
`raw 4 + canonical 3` duplicate complete payloads.

## 2. File hierarchy and version ownership

Implemented hierarchy:

```text
CK file header: magic + archive major/minor + block count
|
+-- metadata block header: id + schema version + codec + flags + length
|   +-- Müsli payload: producer/schema metadata
|
+-- molecule block header: id + schema version + codec + flags + length
|   +-- Müsli payload: complete detached molecule wire record
|
+-- derived block header: id + schema version + codec + flags + length
    +-- Müsli payload: cache payloads and their authoritative validity state
```

The fixed envelope reuses the existing section-header layout. Integers in the
envelope are little-endian, independent of the Müsli payload encoding:

```text
file: 8 bytes b"COSMOL\0\0" + u16 major(2) + u16 minor(0) + u16 block count
block: u16 id + u16 schema + u8 flags + u8 codec + u32 payload length + payload
```

IDs 1 (metadata), 2 (molecule), and 3 (derived) occur exactly once and are all
required (`flags = 1`), schema 1, codec 2. No raw or canonical companion is
written. Unknown optional blocks are skipped by length without decoding their
payloads. Unknown required blocks, duplicate IDs, missing required blocks,
unknown flag bits, and trailing bytes are errors.

The new magic selects the archive 2 reader. The old `CSMOLPKL` magic selects
only the legacy archive 1.x reader; legacy raw inputs retain their own branch.
Neither archive reader falls back to the other after an error.

Codec 2 is CK's `ck-storage-v1`: Müsli storage-compatible encoding, numeric
field/variant IDs and default Binary-mode options. Its current implementation
is pinned to `musli = "=0.1.9"` (`default-features = false`, features `std`,
`alloc`, `storage`); the lockfile resolves `musli`, `musli-core`, and
`musli-macros` to 0.1.9 together. The library version is an implementation
dependency, not a permanent on-disk identity. Metadata records
`codec_contract = "ck-storage-v1"` and molecule/derived schema versions, which
must agree with their respective block headers. The producer package version
is informational.

Metadata, molecule and derived blocks have independent schema constants and
dispatch on `(block ID, block schema)`. All three start at 1; equality of these
initial values does not couple their future evolution. For example, a future
`metadata = 1, molecule = 2, derived = 1` requires a molecule-schema-2 reader,
not new readers for the other two blocks. No schema-2 reader is introduced now.

| Version or identifier | Controls | Does not control |
|---|---|---|
| Archive version `2.0` | Envelope layout, required blocks, allowed combinations and dispatch | CK package release number |
| Block schema version | That block's field types, meanings and interpretation | Other blocks' schemas |
| Block codec identifier | Encoding implementation and frozen configuration | Chemical meaning or migration rules |
| Stable field/variant identifier | Identity inside a Müsli record or enum | Source declaration position |

There is no new raw version and no canonical block in archive 2.0. Legacy
raw/canonical version numbers belong only to legacy decoding. A new block
schema version does not automatically require an archive-version bump.
Changing envelope layout, required-block rules, or an incompatible block
combination does require an archive-version decision.

After this initial V2 definition is finalized, verify the frozen codec bytes
and decoding behavior before upgrading the library or its configuration.
An implementation upgrade that preserves the CK encoding contract retains
codec 2 and `ck-storage-v1`. An incompatible encoding requires a distinct
codec ID/contract, for example codec 3 / `ck-storage-v2`, and explicit decoder
dispatch; it is not represented by changing a third-party version string in
the existing protocol. A package-version change or different producer string
alone does not imply a new archive schema.

## 3. Complete state, without duplicate payloads

The molecule wire record owns topology, ordered typed properties and their
source-defined computed representation, conformers and IDs/order, SGroups,
stereo groups, names, SDF fields/lists, and other represented molecule state.
Property keys/values that admit raw bytes must use byte-preserving wire types,
not lossy UTF-8 conversions. Preserve scalar/vector distinctions, presence
versus absence, sparse records, and floating-point bits where required.

The current model admits UTF-8 `String` values, including embedded NUL, not
arbitrary invalid UTF-8. Atom/bond wire keys and string values are byte vectors;
conversion rejects invalid UTF-8 rather than replacing bytes. This implementation
does not introduce a new runtime text/property model. Molecule properties remain
the current model's string values; atom/bond/SDF-list values retain the modeled
`String`, `Int(i32)`, `UInt(u32)`, `IntVector(Vec<i32>)`, `Double(f64)`, and
`Bool` variants. Floating-point payloads are stored as integer bit patterns.
An extension of the model's admitted value types needs explicit wire conversion
and tests; changing this format alone does not implement those model types.

The derived record owns cache payloads, cache extents, initialization/quality,
and validity metadata. Put each authoritative datum in exactly one block;
cross-block references and dimensions are checked, not duplicated as a second
complete molecule copy. Preserve absent, stale and valid cache states as
distinct states; reading must not silently sanitize or recompute them.

For a molecule consisting of carbon atom 0 and oxygen atom 1:

```text
molecule block:
    atoms: atomic numbers 6 and 8
    bonds: one single bond from 0 to 1
    atom 0 properties: x = Int(123)

derived block:
    existing ring/valence cache payloads and validity, if present
```

Atomic numbers and `x` are not written again into a canonical companion.

Do not derive the archive codec directly on runtime `Molecule`, Arc-backed
blocks, locks, or mutable runtime implementation types. IO owns explicit wire
records and conversion to validated detached values. `cosmolkit` retains the
existing checked runtime construction and cache-restoration boundary.

## 4. Changes Müsli can handle without a new schema reader

Use explicit stable numeric field and enum-variant identifiers. Never infer
persistent identity from declaration order or reuse a retired identifier.

| Change | Conditions for using the existing schema reader |
|---|---|
| Reorder Rust field declarations | Explicit field IDs, types and meanings remain unchanged |
| Add an optional field | Explicit decoding default; absence has the defined meaning `None` |
| Add a list | Explicit decoding default; historical absence genuinely means an empty list |
| Add another defaultable field | Its exact default is a faithful interpretation of historical absence |
| Change runtime storage layout | Explicit conversion preserves the same wire record and semantics |

Example: adding a new optional annotation does not require a special old
reader solely because old files lack that field:

```rust
// Illustrative wire field, not a production type declaration.
#[musli(Binary, name = 20)]
#[musli(default)]
annotation: Option<String>,
```

The new reader interprets a missing field as `None`. Existing fields still
decode by their fixed IDs. Add compatibility tests even when no migration
branch is needed. Never choose a default merely to make decoding succeed:
missing information is not proof of zero, false, valid cache state, or an
inferred conformer order.

Production records use `#[derive(Encode, Decode)]`, explicit numeric
`#[musli(Binary, name = ...)]` identities, and `#[musli(default)]` for faithful
missing-field defaults. Müsli owns field dispatch, order, missing-field handling,
sequence/string decoding, and enum tags. CK does not maintain a parallel map
decoder or add a repeated-field rejection policy. Private wire enums carry
their own fixed variant IDs; explicit model/wire conversions do not encode or
decode numeric tags. The regression suite checks defaults and field reordering
independently of round-trip encoding.

`storage` cannot skip unknown fields inside a recognized record. Consequently,
this policy does not promise that old software can read newly added fields.
The user requirement is one-way: new software reads old files; old software
need not read new files. Optional envelope blocks can still be skipped because
their lengths are known independently of the payload codec.

## 5. Changes requiring explicit versioned decoding/migration

These rules govern future changes after the initial V2 definition is finalized;
they do not introduce older V2 schema readers into the current implementation.

| Change | Required handling |
|---|---|
| Field type changes, e.g. `String` to `Vec<String>` | New block schema version; decode the old schema using its old type, then explicitly convert |
| Integer width or signedness changes | Versioned decoding and checked conversion; preserve or explicitly handle the old range |
| Existing field meaning or unit changes | Versioned interpretation even if the Rust type stays identical |
| New mandatory field with no faithful default | Explicit migration deriving it from preserved information, or a documented unresolved compatibility requirement |
| Remove a stored field | Explicit migration/retention decision; `storage` cannot simply ignore the old field |
| Change an existing enum variant ID, payload type, or meaning | Versioned decoding; never reinterpret old bytes under the new meaning |
| Add an enum variant with a new stable ID | Extend the new reader; no old-schema branch solely for the additive variant if all old IDs/payloads remain unchanged |
| Change computed-property representation | Explicit state conversion and collision policy, not a generic default |
| Incompatible codec/configuration change | A distinct codec contract and decoder dispatch; preserve the old decoder |

A schema reader need not be a separate complete file-reader implementation.
Use one envelope reader and version-specific block decoding/conversion only
where necessary. Decode old layouts into old wire records rather than making
the current record interpret old fields as their new types.

Example, only if the approved migration says an old string is one list item:

```text
block schema v1: decode String("ABC")
    -> explicit v1-to-current conversion
    -> current Vec<String>(["ABC"])

block schema v2: decode Vec<String>(["ABC", "DEF"])
    -> current Vec<String>(["ABC", "DEF"])
```

Do not split `"ABC,DEF"` on commas unless the old field's contract defines
that delimiter. Müsli does not perform this semantic conversion automatically.

An additive enum variant can keep the same schema when the old meaning is
unchanged, because new software still recognizes every old tag. Old software
may reject the new tag; that is allowed. Archive/block version numbers are
not counters incremented for every source-code edit.

## 6. Legacy read policy

Version dispatch is outside the archive-2 Müsli schema. The existing 1.x/raw
branch retains its own reader; archive-2 records do not decode legacy payloads,
emulate legacy layouts, or add fields for legacy encoding compatibility.

Retain reading of previously supported, valid standalone raw 1/2/3 and
sectioned archive 1.0/1.1/1.2, including the supported canonical section
versions. They are legacy inputs, not new-write targets. Preserve their actual
wire definitions, version-specific interpretation, integrity checks and
validation; use the actual legacy version when checking a legacy companion.

```text
from_binary(bytes)
    -> recognize standalone legacy raw or archive header
    -> legacy: old block decoding + explicit conversion
    -> archive 2.0: Müsli 0.1.9 block decoding + explicit model conversion
    -> validated detached values
    -> checked Molecule construction
```

Do not reject an otherwise valid file simply because it is legacy. Do not
require users to rewrite it first. Corrupt files, unsupported future versions,
invalid indices/counts, or failed integrity checks still return structured
errors; "legacy reads without errors" is not permission to accept corruption.
Do not try multiple decoders after a recognized format fails validation.

### Legacy computed-property collision: compatibility requirement

Old CK files can represent an ordinary `__computedProps = "opaque"` entry
and a separate nonempty computed-name set at the same time. A sole new store
using `__computedProps` as the computed-name vector cannot simply overwrite
one with the other.

The user-approved legacy import policy is to read this conflicting slot as an
empty computed-name vector, without reporting an error. The old ordinary value
and independent computed flags are discarded for that property store; all other
property values and molecule data remain intact. This is an explicit lossy
legacy normalization, not a claim of lossless preservation of the conflict.
Nonconflicting legacy properties keep their existing conversion behavior.

The decoder retains the original conflicting value temporarily when validating
legacy raw/canonical companion bytes, so normalization does not create a false
integrity failure. That provenance is not molecule state and is not written to
archive 2.0. The normalized empty vector survives archive-2 round trips. The
archive-2 Müsli codec needs no legacy-specific decoding branch or compatibility
field; format dispatch and normalization stay in the legacy decoder.

### Read limits and failure behavior

Archive-2 input is limited to 256 MiB and at most 1,024 envelope blocks.
Atom/bond counts and ring-cache tables have a 1,000,000-row boundary, checked
on both writing and reading. On reads these row checks run after Müsli has
decoded the DTO, not before its vectors allocate. Other sequences and strings
are framed and decoded by Müsli. The archive byte limit bounds encoded input,
not total decoded heap allocation. Müsli owns allocation and sequence framing;
there are no CK sequence decoders or a second generic validation framework.
These are archive-2 limits, not new limits imposed on legacy readers.

All indices, topology/group references, coordinate dimensions/order, ring-cache
extents, and valence row counts are validated during detached conversion. The
existing runtime constructor validates cache validity/payload relationships.
No sanitation or cache recalculation runs merely because an archive is read.
Errors are returned through the existing `PickleError` API.

#### Tracked hardening: decode-time allocation budgets

Status: follow-up, not implemented and not a 0.5.0 release gate. Investigate
Müsli's allocator/context/decoder controls for per-sequence admission before
allocation and an aggregate allocation budget for nested data. Keep derive-based
decoding rather than adding a CK sequence decoder. Cover the exact limit,
limit-plus-one, and truncated payloads claiming huge sequence lengths. Until
then, do not claim that the post-decode row limit bounds allocations when
reading untrusted uploads or internet data. This does not change legacy readers.

## 7. Validation required before delivery

- Decode fixed valid legacy examples for every supported raw/archive/section
  version into the current model, including typed-property conversion and the
  computed-property collision obligation above.
- Round-trip archive 2.0 complete state, including ordered raw-byte properties,
  every supported value/enum variant, computed records, sparse SGroup roles,
  conformer IDs/order, floating-point bits, cache payloads/extents/validity.
- Exercise defaultable fields absent from schema-1 payloads. Add type/meaning
  migration comparisons when a future schema is actually introduced.
- Reject malformed/truncated payloads, duplicate blocks,
  unknown required blocks, invalid cross-block references and resource-limit
  violations; check trailing-byte handling explicitly.
- Freeze codec byte examples separately from semantic compatibility tests.
  Producer metadata may change without a schema change. A byte difference
  alone is neither proof of incompatibility nor permission to rewrite goldens.
- Keep focused owner tests in their owner; shared fixtures in repository
  `testdata/`; corpus preparation/comparison in `parity-tests/`. Tests must not
  generate their own expectations or invoke reference implementations.

## 8. References and limits

- [Müsli storage documentation](https://docs.rs/musli/latest/musli/storage/):
  missing defaulted fields are tolerated; unknown record fields are not.
- [Müsli format and field-identity guidance](https://github.com/udoprog/musli).
- [Current CK implementation](src/molecule_binary.rs): legacy wire dispatch,
  section validation, canonical reconstruction and checked detached records.
- [RDKit pickle field inventory](rdkit_molecule_pickle_fields.md): a reference
  inventory, not the schema or wire format of this CK-native archive.

No performance improvement or cross-release byte stability is claimed. The
codec regression freezes a `ck-storage-v1` molecule-schema-1 byte example,
currently verified with Müsli 0.1.9. Future implementation upgrades must verify
that contract; a library version change alone does not change the protocol or
authorize regenerating goldens.

Owner validation commands:

```bash
cargo test -p cosmolkit-io --locked --release
cargo test -p cosmolkit-model --locked --release --lib
cargo test -p cosmolkit --locked --release --features op-contracts-strict --test native_archive_v2
cargo test -p cosmolkit --locked --release --features op-contracts-strict --lib original_binary_live_tests
```

Fixed legacy byte fixtures cover raw 1/2/3, archives 1.0/1.1, and archive 1.2
canonical schemas 1/2. Rich legacy-upgrade tests additionally retain typed
properties, computed collisions, PDB metadata, coordinates and companion
checks against each actual raw version. New-format tests cover full modeled
state, sparse SGroup roles, float bits, schema defaults, field ordering,
malformed payloads and resource limits. Public tests verify unchanged APIs and
input sharing. These local regressions are not a claim of a full corpus or
workspace acceptance.
