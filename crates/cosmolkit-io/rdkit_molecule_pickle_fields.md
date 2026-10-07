# RDKit molecule pickle: fields, types and ownership

This is a source inventory of RDKit molecule binary serialization, not a
COSMolKit archive proposal. It describes the vendored `MolPickler` writer with
pickle version **16.3.0**, inspected on 2026-10-06. It does not claim that a
pickle is a complete snapshot of RDKit's in-memory objects, or that COSMolKit
already supports these fields. Reaction serialization is a separate format.

Every field emitted by the normal molecule writer and its reached query,
property and built-in custom-property writers is listed below. The alternate
`OLD_PICKLE` writer is listed separately. Historical reader branches are not
mistaken for fields emitted by the current writer. Arbitrary property keys and
application-registered custom handlers are extension points, not a finite
catalogue of fixed RDKit fields.

## 1. Ownership hierarchy

The hierarchy is semantic ownership, not the order of all bytes in the stream.
For example, atom property tables are written together after the atom records,
but each table still belongs to one atom, not to the molecule property table.

```text
Molecule pickle
├── Format header and molecule counts
├── Molecule.properties: property table
├── Atoms[i]
│   ├── Atomic number, aromaticity, implicit-H policy, charge, chirality,
│   │   hybridization, explicit H count, valences, radicals, isotope
│   ├── Atom map number and dummy label
│   ├── Query: recursive query-node tree
│   ├── MonomerInfo / PDBResidueInfo: atom-owned record
│   └── properties: this atom's property table and compact special properties
├── Bonds[j]
│   ├── Begin/end atom indices, type, aromaticity, conjugation, direction
│   ├── Stereo kind and stereo-reference atom indices
│   ├── Query: recursive query-node tree
│   ├── Molfile endpoint/attachment strings
│   └── properties: this bond's property table and compact special properties
├── RingInfo
│   └── Ring-finding kind and ordered atom-index rows
├── SubstanceGroups[k]
│   ├── properties: this SGroup's property table
│   ├── Atom, parent-atom and bond indices
│   ├── Brackets: three 3D points per bracket
│   ├── CStates: bond index and conditional 3D vector
│   └── AttachPoints: atom index, leaving-atom index and string ID
├── StereoGroups[k]
│   └── Group kind, conditional output ID, atom and bond indices
└── Conformers[c]
    ├── is3D, conformer ID and atom count
    ├── Positions[i]: x, y, z for this conformer and atom
    └── properties: this conformer's property table
```

| Property owner | Address/association in the pickle | Selection for its generic table |
|---|---|---|
| Molecule | One molecule table | `MolProps` |
| Atom | One table for each atom, in atom order | `AtomProps` |
| Bond | One table for each bond, in bond order | `BondProps` |
| Conformer | One table for each conformer, in conformer order | `MolProps`, and conformers not excluded |
| SGroup | A table inside each SGroup record | Always visited when that SGroup is written; its own filters apply |
| Query value predicate | One key/type/value record inside a query node | Structural query payload, not an optional owner property table |

Property keys are **not atom-index strings**. An atom's index is determined by
its position in the atom sequence. Its property table is associated by the
same sequence position. Bond endpoint and membership references are explicit
integer indices. Coordinates are owned by a conformer and indexed by atom,
not stored as an ordinary atom property.

Sources: [main writer](../../third_party/rdkit/Code/GraphMol/MolPickler.cpp#L1185),
[property writers](../../third_party/rdkit/Code/GraphMol/MolPickler.cpp#L110).

## 2. Encoding vocabulary

| Notation | Exact source/write type | Encoding |
|---|---|---|
| `tag` | `MolPickler::Tags` converted to `unsigned char` | One byte, except the initial `VERSION` header field below |
| `char` | C++ `char` | One byte; signedness is platform-dependent, so it is not automatically an unsigned count |
| `signed char` / `uint8_t` / `int8_t` | Explicit source type | One byte |
| `int32_t`, `uint32_t`, `int16_t`, `uint16_t` | Explicit fixed-width integer | Respectively 4, 4, 2, 2 bytes, little-endian |
| `int`, `unsigned int` | Native C++ source type | `sizeof(T)` bytes, little-endian; normally 4 bytes on supported builds, not a source-level `int32_t` typedef |
| `bool` | Native C++ `bool` | `sizeof(bool)` bytes; normally one byte; distinct from flags packed into `char` |
| `float`, `double` | C++ floating type | `sizeof(T)` bytes with endian conversion; normally binary32/binary64 (4/8 bytes) |
| `String` | `std::string` | `unsigned int` byte length, then exactly that many raw bytes; no terminator |
| `Vector<V>` | `std::vector<V>` | `boost::uint64_t` element count, followed by each element's encoding |
| `T` | Writer's index/count template parameter | `unsigned char` if molecule atom count is at most 255; otherwise `int32_t` |
| `C` | Conformer coordinate template parameter | `float` by default; `double` with `CoordsAsDouble` |
| Framed block | Temporary `stringstream` payload | `int32_t` byte length, followed by those raw bytes |

`T` is selected from **molecule atom count**, not separately from each list's
size or bond count. The writer casts several counts and indices to `T` or
`char`; this document does not replace those casts with wider hypothetical
types. Normal molecule counts are fixed-width writes, not varints. The nested
bit-vector payload in Section 11 is the exception using packed integers.

Strings may contain embedded NUL and non-UTF-8 bytes. Neither `streamWrite`
nor its string reader imposes UTF-8 validation. Ordinary scalar writes perform
little-endian conversion and use `sizeof(T)`; there is no general portable
"all integers are u32" rule.

Sources: [tag encoding](../../third_party/rdkit/Code/GraphMol/MolPickler.cpp#L42),
[scalar/string/vector encoding](../../third_party/rdkit/Code/RDGeneral/StreamOps.h#L263),
[template selection](../../third_party/rdkit/Code/GraphMol/MolPickler.cpp#L1043).

## 3. Options and complete outer stream

### 3.1 Writer options

| Option | Type/value | Effect |
|---|---|---|
| `NoProps` | `unsigned int`, `0` | Initial global default: no optional molecule/atom/bond/conformer property tables; structural fields remain |
| `MolProps` | `0x1` | Molecule and conformer property tables |
| `AtomProps` | `0x2` | Atom property tables and compact atom special properties |
| `BondProps` | `0x4` | Bond property tables and compact bond special properties |
| `QueryAtomData` | `0x2` | Deprecated alias of `AtomProps`, not a second independent category |
| `PrivateProps` | `0x10` | Include private keys in selected generic tables |
| `ComputedProps` | `0x20` | Include computed keys in selected generic tables |
| `AllProps` | `0x0000FFFF` | Includes property category/filter bits; does not include the two higher coordinate-option bits |
| `CoordsAsDouble` | `0x00010000` | Use `double` conformer coordinates |
| `NoConformers` | `0x00020000` | Omit conformers and their properties |

The global default can be changed through `setDefaultPickleProperties`.
`PrivateProps`/`ComputedProps` alone do not select a property owner category.
The complete options word is **not written as a field**; section presence and
tags encode its effects. SGroup properties use separate fixed filtering
arguments, explained in Section 8.

Source: [options](../../third_party/rdkit/Code/GraphMol/MolPickler.h#L50).

### 3.2 All outer fields, in write order

| Level | Field or framing marker | Write type | Condition/content |
|---|---|---|---|
| Format | Endian sentinel | `int32_t` | `0xDEADBEEF` bit pattern |
| Format | `VERSION` marker | `int` | Explicit cast of `VERSION == 0`; unlike following tags, not a one-byte write |
| Format | Version major | `int32_t` | `16` |
| Format | Version minor | `int32_t` | `3` |
| Format | Version patch | `int32_t` | `0` |
| Molecule | Atom count | `int32_t` | Always |
| Molecule | Bond count | `int32_t` | Always |
| Molecule | Flags | `unsigned char` | Current writer always writes `0x80`; reader names bit 7 `includeCoords` |
| Molecule | `BEGINATOM` | `tag` | Followed by exactly atom-count records, Section 4 |
| Molecule | `BEGINBOND` | `tag` | Followed by exactly bond-count records, Section 5 |
| RingInfo | Ring-kind tag and ring records | `tag`, then Section 7 | Only if RingInfo exists and is initialized |
| Molecule | `BEGINSGROUP` | `tag` | Only if SGroup list is nonempty |
| Molecule | SGroup count | `int32_t` | Then that many Section 8 records |
| Molecule | `BEGINSTEREOGROUP` | `tag` | Only if stereo-group list is nonempty; then Section 9 |
| Molecule | `BEGINCONFS` or `BEGINCONFS_DOUBLE` | `tag` | Unless `NoConformers`; selects `C` |
| Molecule | Conformer-block byte length | `int32_t` | Framed payload, even when conformer count is zero |
| Molecule | Conformer count inside that block | `int32_t` | Then Section 6 records |
| Conformer block | `BEGINCONFPROPS` | `tag` | If `MolProps`; follows all conformer coordinate records |
| Conformer block | Conformer-properties byte length | `int32_t` | One nested frame containing one table per conformer |
| Molecule | `BEGINPROPS`, property-block length, payload, `ENDPROPS` | `tag`, `int32_t`, Section 10, `tag` | If `MolProps` and temporary stream nonempty; a zero-count table still occupies two bytes |
| Atoms | `BEGINATOMPROPS`, block length, payload, `ENDPROPS` | `tag`, `int32_t`, Section 10, `tag` | If `AtomProps` and at least one atom writes a generic or compact property |
| Bonds | `BEGINBONDPROPS`, block length, payload, `ENDPROPS` | `tag`, `int32_t`, Section 10, `tag` | If `BondProps` and at least one bond writes a generic or compact property |
| Molecule | `ENDMOL` | `tag` | Always |

Atom/bond property blocks, when present, contain tables for **every** atom/bond,
including empty ones, in the corresponding owner order. There is no separate
index field in each table. There is no per-record `ENDATOM` or `ENDBOND` in the
normal writer; those belong to the alternate legacy writer.

Sources: [header](../../third_party/rdkit/Code/GraphMol/MolPickler.cpp#L1025),
[outer layout](../../third_party/rdkit/Code/GraphMol/MolPickler.cpp#L1185).

## 4. Atom-owned fields

### 4.1 Fixed prefix and boolean mask

| Atom field | C++ semantic/source type | Write type | Condition |
|---|---|---|---|
| Atomic number | `int`, from `getAtomicNum()` | `unsigned char` | Always; source takes `% 256` before assignment |
| Atom flags | Packed mask | `char` | Always; bits below |
| Atom-data presence mask (`propFlags`) | Packed mask | `int32_t` | Always; governs Section 4.2, not the optional property table |

| Atom flag bit | Atom-owned meaning | Semantic type |
|---|---|---|
| 6 | `isAromatic` | `bool` |
| 5 | `noImplicit` (implicit-H policy) | `bool` |
| 4 | Has query | `bool`, followed by query payload |
| 3 | Atom map number successfully retrieved | Presence of `int`-interpretable map property |
| 2 | Has `dummyLabel` property | Presence flag |
| 1 | Has monomer info | Pointer presence, not pointer bytes |
| 0, 7 | Not set by this writer | No current payload |

### 4.2 Conditional atom-data scalars, in write order

| Atom field | C++ semantic/source type | Write type | Presence-mask bit and source condition | Reader value when omitted |
|---|---|---|---|---|
| Formal charge | `int` getter | `signed char` | 1; converted byte is nonzero | `0` |
| Chiral tag | `Atom::ChiralType` enum | `char` | 2; converted byte is nonzero | Zero enum value |
| Hybridization | `Atom::HybridizationType` enum | `char` | 3; converted byte differs from `Atom::SP3` | `SP3` |
| Explicit H count | `unsigned int` getter | `char` | 4; converted byte is nonzero | `0` |
| Cached explicit valence | `d_explicitValence: std::int8_t` | `char` | 5; original member is greater than zero | `0` |
| Cached implicit valence | `d_implicitValence: std::int8_t` | `char` | 6; original member is greater than zero | `0` |
| Radical electron count | `unsigned int` getter | `char` | 7; converted byte is nonzero | `0` |
| Isotope | `unsigned int` | `unsigned int` | 8; greater than zero | `0` |

Bit 0's mass-delta `float` writer is commented out. It is not an additional
current field. A reader branch for that older field is not evidence that the
current writer preserves an independent mass value. Positive cached valences
are saved, but omitted negative/uninitialized values do not retain their
original sentinel. This is not a full cache-validity snapshot.

### 4.3 Remaining atom payloads, in write order

| Atom field | C++ type | Write representation | Condition |
|---|---|---|---|
| Query | `Query<int, Atom const *, true>` tree | `BEGINQUERY`, Section 12 node, `ENDQUERY` | Has query |
| Atom map number | Retrieved as `int`; string values can be interpreted by the helper | `ATOM_MAPNUMBER`, `char` for values 0..127; otherwise byte `0xFF` followed by `int` | Successful map retrieval; independent of optional property flags |
| Dummy label | `std::string` via property getter | `ATOM_DUMMYLABEL`, `String` | Property exists |
| PDB monomer info | `AtomPDBResidueInfo` | `BEGIN_PDB_RESIDUE`, Section 13.1, `END_PDB_RESIDUE` | Monomer type is `PDBRESIDUE` |
| Other monomer info | `AtomMonomerInfo` | `BEGIN_ATOM_MONOMER_INFO`, Section 13.2, `END_ATOM_MONOMER_INFO` | Present and not `PDBRESIDUE` |
| Ordinary/special properties | `RDProps` dictionary | Section 10.3 in the later atom-property block | `AtomProps`, subject to block/filter rules |

Source: [atom data](../../third_party/rdkit/Code/GraphMol/MolPickler.cpp#L1607),
[atom writer](../../third_party/rdkit/Code/GraphMol/MolPickler.cpp#L1732),
[member/getter types](../../third_party/rdkit/Code/GraphMol/Atom.h).

## 5. Bond-owned fields

| Bond field | C++ semantic/source type | Write type | Condition/order |
|---|---|---|---|
| Begin atom index | `unsigned int` | `T` after atom-index mapping | Always, first |
| End atom index | `unsigned int` | `T` after atom-index mapping | Always, second |
| Bond flags | Packed mask | `char` | Always, third |
| Bond type | `Bond::BondType` enum | `char` | If not `SINGLE` |
| Bond direction | `Bond::BondDir` enum | `char` | If not `NONE` |
| Bond stereo | `Bond::BondStereo` enum | `char` | If not `STEREONONE` |
| Stereo-reference atom count | `INT_VECT` length | `char` after size conversion | Only with non-`STEREONONE` stereo |
| Stereo-reference atom indices | `INT_VECT == std::vector<int>` | One `T` per entry | In stored vector order; source casts each entry directly, without `atomIdxMap` at this site |
| Query | `Query<int, Bond const *, true>` tree | `BEGINQUERY`, Section 12 node, `ENDQUERY` | Has query |
| `_MolFileBondEndPts` | `std::string` property | `String` | Retrieved endpoint string is nonempty |
| `_MolFileBondAttach` | `std::string` property | `String` | Immediately after nonempty endpoint string; fallback is source-defined `"ALL"` |
| Ordinary/special properties | `RDProps` dictionary | Section 10.4 in later bond-property block | `BondProps`, subject to block/filter rules |

| Bond flag bit | Bond-owned meaning | Semantic type |
|---|---|---|
| 6 | Aromatic | `bool` |
| 5 | Conjugated | `bool` |
| 4 | Has query | `bool` |
| 3 | Non-default bond type follows | Presence flag |
| 2 | Non-default direction follows | Presence flag |
| 1 | Non-default stereo and its reference list follow | Presence flag |
| 0 | `_MolFileBondEndPts` was successfully retrieved | Presence flag |
| 7 | Not set by this writer | No current payload |

Do not conflate the bit-0 condition with the subsequent string-write condition:
the former tests successful retrieval, the latter tests **nonempty** content.
This inventory records both source conditions, without asserting that an
empty-but-present endpoint string round-trips. Structural endpoint strings can
also be eligible for the generic bond property table; only Section 10.4's five
compact keys are excluded from that generic table.

Source: [bond writer](../../third_party/rdkit/Code/GraphMol/MolPickler.cpp#L2040),
[bond types](../../third_party/rdkit/Code/GraphMol/Bond.h).

## 6. Conformer-owned fields

| Level | Field | C++ type | Write type | Condition |
|---|---|---|---|---|
| Conformer | `is3D` | `bool` | `char` | Always per included conformer |
| Conformer | ID | `unsigned int` | `int32_t` | Always |
| Conformer | Atom count | `unsigned int` | `T` | Always |
| Conformer / atom position | `x` | `Point3D` component, `double` | `C` | For every stored position |
| Conformer / atom position | `y` | `double` | `C` | After `x` |
| Conformer / atom position | `z` | `double` | `C` | After `y`, even for a 2D conformer |
| Conformer | Properties | `RDProps` dictionary | Section 10 table | `MolProps`; in nested conformer-properties frame |

`NoConformers` removes this entire collection. Default coordinate serialization
casts from double to float; double bitwise preservation is not implied by
default pickling. `CoordsAsDouble` affects conformer coordinates, **not** SGroup
bracket/CState coordinates. There is no per-position atom-index field; order
associates each coordinate triple with its atom.

Source: [conformer writer](../../third_party/rdkit/Code/GraphMol/MolPickler.cpp#L1802).

## 7. Molecule-owned RingInfo

| Level | Field | C++ type | Write type | Condition |
|---|---|---|---|---|
| RingInfo | Finding kind | `FIND_RING_TYPE` enum | One tag: `BEGINFASTFIND`, `BEGINSSSR`, `BEGINSYMMSSSR`, or `BEGINFINDOTHERORUNKNOWN` | Initialized RingInfo only |
| RingInfo | Ring count | `numRings()` converted to `std::uint32_t` | `uint32_t` | After finding-kind tag |
| Ring row | Atom count | `INT_VECT` size | `T` | Per ring |
| Ring row | Ordered atom indices | `INT_VECT` | One mapped `T` per atom | Per ring, preserving row order |

The normal writer does not separately emit bond-ring rows, per-atom/per-bond
membership arrays, fused-ring caches or ring-family payloads. The reader derives
bond indices from atom rows and the molecule's edges. Initialized empty ring
state still has a kind tag and zero ring count; uninitialized state has no ring
section. `OtherOrUnknown` is one fallback tag, not an exact arbitrary enum value.

Sources: [kind dispatch](../../third_party/rdkit/Code/GraphMol/MolPickler.cpp#L1229),
[ring writer/reader](../../third_party/rdkit/Code/GraphMol/MolPickler.cpp#L2233).

## 8. Molecule-owned SGroups and their child fields

The outer SGroup count is `int32_t`. Each SGroup has the following fields in
this order; all listed membership indices refer to molecule atoms/bonds.

| Level | Field | C++ type | Write type |
|---|---|---|---|
| SGroup | Properties | Inherited `RDProps` dictionary | Section 10 table, `uint16_t` property count |
| SGroup | Atom count | `std::vector<unsigned int>` size | `T` |
| SGroup | Atoms | `std::vector<unsigned int>` | Mapped atom indices, one `T` each |
| SGroup | Parent-atom count | `std::vector<unsigned int>` size | `T` |
| SGroup | Parent atoms | `std::vector<unsigned int>` | Mapped atom indices, one `T` each |
| SGroup | Bond count | `std::vector<unsigned int>` size | `T` |
| SGroup | Bonds | `std::vector<unsigned int>` | Mapped bond indices, one `T` each |
| SGroup | Bracket count | `std::vector<Bracket>` size | `T` |
| Bracket | Three points, each with `x`, `y`, `z` | `std::array<Point3D, 3>`; double components | Nine `float` values per bracket, in point order |
| SGroup | CState count | `std::vector<CState>` size | `T` |
| CState | Bond index | `unsigned int` | Mapped bond index as `T` |
| CState | Vector `x`, `y`, `z` | `Point3D`, double components | Three `float` values **only if SGroup property `TYPE` is `"SUP"`** |
| SGroup | Attach-point count | `std::vector<AttachPoint>` size | `T` |
| Attach point | Attached atom `aIdx` | `unsigned int` | Mapped atom index as `T` |
| Attach point | Leaving atom `lvIdx` | `int` | `signed int`; `-1` when absent, otherwise mapped atom index |
| Attach point | ID | `std::string` | `String`, including an empty ID |

The SGroup property call is exactly `streamWriteProps<uint16_t>(ss, sgroup)`:
`savePrivate=false`, `saveComputed=false`, **no custom handlers**, regardless of
the molecule's property flags. Thus `AllProps` does not change these filters.

SGroup names such as `TYPE` and other metadata keys live in its property table,
not in a second fixed name/value schema in `MolPickler`. The table supports the
ordinary serializable types in Section 10.2. Keys are open-ended, so listing
only familiar SDF keys would not enumerate all writable properties.
`Contained`/`Crossing` (`CBOND`/`XBOND`) is not a separate enum field here;
the SGroup bond membership is written. Neither its owner pointer nor
`d_isValid` is emitted as an independent field.

Source: [SGroup writer](../../third_party/rdkit/Code/GraphMol/MolPickler.cpp#L2327),
[SGroup value types](../../third_party/rdkit/Code/GraphMol/SubstanceGroup.h).

## 9. Molecule-owned enhanced stereo groups

| Level | Field | C++ type | Write type | Condition |
|---|---|---|---|---|
| Stereo-group collection | Count | `std::vector<StereoGroup>` size | `T` | Nonempty collection section |
| Stereo group | Kind | `StereoGroupType` enum | `T` | `STEREO_ABSOLUTE=0`, `STEREO_OR=1`, `STEREO_AND=2` |
| Stereo group | Output/write ID | `unsigned` | `T` | Only non-absolute groups |
| Stereo group | Atom count | `std::vector<Atom *>` size | `T` | Always |
| Stereo group | Atom references | Atom pointers in memory | Mapped atom indices as `T` | Per member, no pointer bytes |
| Stereo group | Bond count | `std::vector<Bond *>` size | `T` | Always |
| Stereo group | Bond references | Bond pointers in memory | Mapped bond indices as `T` | Per member, no pointer bytes |

The writer calls `assignStereoGroupIds` on its by-value group vector before
writing IDs. It does not independently write both `d_readId` and `d_writeId`.

Source: [stereo writer](../../third_party/rdkit/Code/GraphMol/MolPickler.cpp#L2501),
[stereo-group types](../../third_party/rdkit/Code/GraphMol/StereoGroup.h).

## 10. Owner-local property dictionaries

### 10.1 Complete generic property-table fields

This same encoding is used for molecule, atom, bond, conformer and SGroup
tables, with their different selection rules. Query value predicates reuse
the individual property record, **without a table count**.

| Level | Field | C++ type | Write type |
|---|---|---|---|
| Owner property table | Serializable selected property count | Counter | `uint16_t` in all current MolPickler table calls |
| Property record | Key | `std::string` | `String` |
| Property record | Value type discriminator | `DTags` constant | `unsigned char` |
| Property record | Value | `RDValue`, actual stored type | One of Section 10.2's payloads |
| Custom property record | Handler name | Handler `const char *`, converted to `std::string` | `String`, after discriminator `0xFE` |
| Custom property record | Handler payload | Handler-defined | Section 11 for the built-in handler; arbitrary for registered extensions |

Selection uses `RDProps::getPropList(savePrivate, saveComputed)` and the
caller's ignore set, then excludes nonserializable values. Records follow the
underlying dictionary iteration order. There is no end-tag per ordinary
property and no independent owner index inside a property record.

### 10.2 Every ordinary writable property value type

| Discriminator | Semantic C++ value type | Exact payload |
|---|---|---|
| `0` / `StringTag` | `std::string` | `String` |
| `1` / `IntTag` | `int` | `int` |
| `2` / `UnsignedIntTag` | `unsigned int` | `unsigned int` |
| `3` / `BoolTag` | `bool` | `bool`, not a string or packed property flag |
| `4` / `FloatTag` | `float` | `float` |
| `5` / `DoubleTag` | `double` | `double` |
| `6` / `VecStringTag` | `std::vector<std::string>` | `uint64` element count; then one `String` per element |
| `7` / `VecIntTag` | `std::vector<int>` | `uint64` element count; then one `int` per element |
| `8` / `VecUIntTag` | `std::vector<unsigned int>` | `uint64` element count; then one `unsigned int` per element |
| `10` / `VecFloatTag` | `std::vector<float>` | `uint64` element count; then one `float` per element |
| `11` / `VecDoubleTag` | `std::vector<double>` | `uint64` element count; then one `double` per element |
| `0xFE` / `CustomTag` | Handler-supported `RDValue` object | Handler name and handler-specific payload |

`VecBoolTag=9` is declared, but has **no ordinary serialization admission,
writer case or reader case** here. It must not be claimed as supported merely
because the constant exists. `EndTag=0xFF` is also declared, but is not written
as a property-table terminator. Arbitrary `AnyTag` values without a matching
custom handler are skipped by table serialization, not automatically encoded.

### 10.3 Atom compact property fields

After **each atom's** generic table, a `uint8_t` presence mask is written,
then the present values in the following order as **`int16_t`**:

| Atom property key | Mask | Getter/write type | Reader property type |
|---|---|---|---|
| `molStereoCare` | `0x01` | `int16_t` via `getPropIfPresent` | `int` |
| `molParity` | `0x02` | `int16_t` | `int` |
| `molInversionFlag` | `0x04` | `int16_t` | `int` |
| `_ChiralityPossible` | `0x08` | `int16_t` | `int` |

These four keys are excluded from the atom's generic table. `molAtomMapNumber`
and `dummyLabel` are also excluded there because Section 4.3 writes them in
the structural atom record. The compact lookup does not apply the generic
private/computed filters; successful value retrieval governs its mask.

### 10.4 Bond compact property fields

After **each bond's** generic table, a `uint8_t` presence mask is written,
then present values in the following order as **`int8_t`**:

| Bond property key | Mask | Getter/write type | Reader property type |
|---|---|---|---|
| `_MolFileBondType` | `0x01` | `int8_t` via `getPropIfPresent` | `int` |
| `_MolFileBondStereo` | `0x02` | `int8_t` | `int` |
| `_MolFileBondCfg` | `0x04` | `int8_t` | `int` |
| `_MolFileBondQuery` | `0x08` | `int8_t` | `int` |
| `molStereoCare` | `0x10` | `int8_t` | `int` |

These five keys are excluded from the generic bond table. As for atom compact
properties, successful retrieval rather than generic private/computed filtering
governs their inclusion. Retrieval/conversion is not proof that the original
`RDValue` type is preserved: these compact fields restore as `int`.

### 10.5 Private and computed properties: where the metadata lives

Each `RDProps` owner has its **own** dictionary and computed-key list.
`__computedProps` is the reserved property key whose conventional value is
`STR_VECT == std::vector<std::string>` containing computed property names.
It is encoded as a normal key/type/value record (type `VecStringTag`) when
selected. There is **no additional per-property `computed: bool` wire field**.

Private keys start with `_`. Without `PrivateProps`, they are filtered out of
the applicable generic tables. Without `ComputedProps`, `getPropList` reads
the owner's computed-key list, excludes those keys, and also excludes the
reserved list key when that list was retrieved. Preserving that list therefore
normally requires both private and computed inclusion, plus the owner's
category bit. Compact and structural special fields have their separate rules.

Thus a computed atom property is **atom-local**, a computed molecule property
is **molecule-local**, and identical key names on different owners are not one
shared slot. There is no finite list of all allowable property names: keys are
user-extensible. Their types come from actual values, not from the spelling of
the key. An arbitrary ordinary value under `__computedProps` is not a second
independently serialized computed-marking channel.

Sources: [generic property admission/writer](../../third_party/rdkit/Code/RDGeneral/StreamOps.h#L388),
[compact properties](../../third_party/rdkit/Code/GraphMol/MolPickler.cpp#L123),
[owner-local filtering](../../third_party/rdkit/Code/RDGeneral/RDProps.h#L48).

## 11. Built-in custom property: ExplicitBitVect

`MolPickler` initializes its handler list with
`DataStructsExplicitBitVecPropHandler`. This handles an actual `ExplicitBitVect`
stored as an `RDValue`, wherever that table/query call passes the handler list.
This is **not** an automatically stored molecule fingerprint field.

| Level | Field | Write type/encoding |
|---|---|---|
| Property | Key | `String` |
| Property | Type tag | `unsigned char`, `0xFE` |
| Property | Handler name | `String`, exactly `"ExplicitBVProp"` |
| Property | Nested bit-vector pickle | `String` length and bytes from `ExplicitBitVect::toString()` |
| Bit-vector payload | Version marker | `int32_t`, `-32` (`-ci_BITVECT_VERSION`) |
| Bit-vector payload | Total bit count | `int32_t`, converted from bit-vector size |
| Bit-vector payload | Set-bit count | `int32_t` |
| Bit-vector payload | Zero run before each set bit | Packed unsigned 32-bit integer input; variable 1..4-byte encoding |
| Bit-vector payload | Trailing zero run | Same packed encoding; written even after the last set bit, or for an all-zero vector |

Packed integers use RDKit's `appendPackedIntToStream`, not a generic LEB128
assumption: low bits encode byte count, with range-dependent offsets. On-bit
positions are traversed in ascending order; each run is `i - prev - 1`, with
initial `prev=-1`. The final run is `size - prev - 1`.

Applications may register additional `CustomPropHandler`s. Their names and
payloads cannot be exhaustively specified from `MolPickler` alone; each actual
registered handler is another source dependency. SGroup's call does not pass
this handler list.

Sources: [handler registration](../../third_party/rdkit/Code/GraphMol/MolPickler.cpp#L291),
[built-in handler](../../third_party/rdkit/Code/DataStructs/DatastructsStreamOps.h#L39),
[nested writer](../../third_party/rdkit/Code/DataStructs/ExplicitBitVect.cpp#L207),
[packed integer encoding](../../third_party/rdkit/Code/RDGeneral/StreamOps.h#L114).

## 12. Atom/bond query-node fields

Query ownership follows the enclosing atom or bond. Recursive children belong
to their parent node; a recursive-structure node additionally owns a nested
query molecule pickle. A query property predicate's expected value is not the
target atom/bond's property dictionary.

### 12.1 Common node fields, in write order

| Query-node field | C++ type | Write encoding | Condition |
|---|---|---|---|
| Description | `std::string` | `String` | Always |
| Type label | `std::string` | `QUERY_TYPELABEL`, `String` | Only nonempty label |
| Negation | `bool` | `QUERY_ISNEGATED` tag only | Only when true; absence means false |
| Query kind | Classified C++ query subtype | One `tag` | Always; kinds below |
| Kind-specific operands | Subtype-specific | Section 12.2 | Per kind |
| Child count | Iterator distance | `QUERY_NUMCHILDREN`, `unsigned char` cast | Always, including zero |
| Children | Query child vector | One complete node encoding per child | Recursively, in vector order |

Only atom/bond top-level query records have the external `BEGINQUERY` and
`ENDQUERY` wrapper. Children do not get an extra wrapper. The count is byte-sized
but the source visits all children; this inventory does not invent an overflow
admission check. Unknown query subtypes throw `MolPicklerException`.

### 12.2 Every currently classified query kind and operand

| Kind tag | C++ query subtype/value type | Operand encoding after kind tag |
|---|---|---|
| `QUERY_AND` | `AndQuery<int, owner const *, true>` | No scalar operand; common children follow |
| `QUERY_OR` | `OrQuery<int, owner const *, true>` | No scalar operand |
| `QUERY_XOR` | `XOrQuery<int, owner const *, true>` | No scalar operand |
| `QUERY_NULL` | Base `Query<int, owner const *, true>` | No scalar operand |
| `QUERY_EQUALS` | `EqualityQuery`, value and tolerance | `QUERY_VALUE`, `int32_t` value, `int32_t` tolerance |
| `QUERY_GREATER` | `GreaterQuery` | `QUERY_VALUE`, `int32_t` value, `int32_t` tolerance |
| `QUERY_GREATEREQUAL` | `GreaterEqualQuery` | `QUERY_VALUE`, `int32_t` value, `int32_t` tolerance |
| `QUERY_LESS` | `LessQuery` | `QUERY_VALUE`, `int32_t` value, `int32_t` tolerance |
| `QUERY_LESSEQUAL` | `LessEqualQuery` | `QUERY_VALUE`, `int32_t` value, `int32_t` tolerance |
| `QUERY_ATOMRING` | `AtomRingQuery` | `QUERY_VALUE`, `int32_t` value, `int32_t` tolerance |
| `QUERY_RANGE` | `RangeQuery`, bounds/tolerance/open-end booleans | `QUERY_VALUE`, `int32_t` lower, `int32_t` upper, `int32_t` tolerance, `char` open-end mask |
| `QUERY_SET` | `SetQuery`, converted to `std::set<int32_t>` | `QUERY_VALUE`, `int32_t` member count, one `int32_t` per member in set order |
| `QUERY_PROPERTY` | `HasPropQuery`, property name | `QUERY_VALUE`, `String` name |
| `QUERY_PROPERTY_WITH_VALUE` | `HasPropWithValueQueryBase`, `PairHolder` and `double` tolerance | `QUERY_VALUE`, `double` tolerance **first**, then one Section 10 property record (key/type/value), using custom handlers |
| `QUERY_RECURSIVE` | `RecursiveStructureQuery` | `QUERY_VALUE`, complete nested molecule pickle, including its header |

Range open-end mask: bit 1 is `lowerOpen`; bit 0 is `upperOpen`. The recursive
molecule uses the overload with **global default property flags**, not automatic
inheritance of the enclosing molecule's explicitly supplied flags.

The writer also has a `QueryDetails` variant arm for `(tag, int32_t)` producing
`tag`, `QUERY_VALUE`, one `int32_t`. The current `getQueryDetails` classifier
does not produce that arm. Likewise the declared `QUERY_BOOL` tag is not a
currently classified extra boolean-query payload. A query property value uses
`streamWriteProp` directly, not table filtering; unsupported value types do not
gain an invented encoding by being placed in a query.

Source: [classifier and complete node writer](../../third_party/rdkit/Code/GraphMol/MolPickler.cpp#L313).

## 13. Atom-owned monomer and PDB residue records

These are attached to **each atom**, even when many atoms describe the same
residue. They are not a molecule-level normalized residue table. They are
structural payloads and are not gated by `AtomProps`, `PrivateProps` or
`ComputedProps`.

### 13.1 AtomPDBResidueInfo

Name and monomer type are unconditional. Each later field has its own named
tag followed by the value, only when the source condition holds.

| Atom / PDB-info field | C++ semantic/source type | Write type | Condition/tag |
|---|---|---|---|
| Atom name | `std::string` | `String` | Always, first |
| Monomer type | `AtomMonomerInfo::AtomMonomerType` | `unsigned int` | Always; `PDBRESIDUE` for this branch |
| Serial number | Member `unsigned int`, getter returns `int` | `int` | Nonzero; `ATOM_PDB_RESIDUE_SERIALNUMBER` |
| Alternate location | `std::string` | `String` | Nonempty; `ATOM_PDB_RESIDUE_ALTLOC` |
| Residue name | `std::string` | `String` | Nonempty; `ATOM_PDB_RESIDUE_RESIDUENAME` |
| Residue number | `int` | `int` | Nonzero; `ATOM_PDB_RESIDUE_RESIDUENUMBER` |
| Chain ID | `std::string` | `String` | Nonempty; `ATOM_PDB_RESIDUE_CHAINID` |
| Insertion code | `std::string` | `String` | Nonempty; `ATOM_PDB_RESIDUE_INSERTIONCODE` |
| Occupancy | `double` | `double` | Nonzero; `ATOM_PDB_RESIDUE_OCCUPANCY` |
| Temperature factor | `double` | `double` | Nonzero; `ATOM_PDB_RESIDUE_TEMPFACTOR` |
| Hetero-atom flag | `bool` | `char` | True; `ATOM_PDB_RESIDUE_ISHETEROATOM` |
| Secondary structure | `unsigned int` | `unsigned int` | Nonzero; `ATOM_PDB_RESIDUE_SECONDARYSTRUCTURE` |
| Segment number | `unsigned int` | `unsigned int` | Nonzero; `ATOM_PDB_RESIDUE_SEGMENTNUMBER` |
| Monomer class | `std::string` | `String` | Nonempty; `ATOM_PDB_RESIDUE_MONOMERCLASS` |

The outer `END_PDB_RESIDUE` terminates the tagged record. The occupancy test is
**nonzero**, not "different from 1.0": the reader constructs an info object
whose default occupancy is 1.0, so a source occupancy of zero is omitted and
is not independently preserved by this writer. Other defaults follow the
constructor; absent strings are empty, numeric fields zero and hetero false.
Element and formal charge belong to the atom fields, not duplicated PDB fields.

### 13.2 Non-PDB AtomMonomerInfo

| Atom / monomer-info field | C++ type | Write type | Condition/tag |
|---|---|---|---|
| Atom name | `std::string` | `String` | Always |
| Monomer type | `AtomMonomerType` enum | `uint8_t` | Always; enum declares `UNKNOWN=0`, `PDBRESIDUE=1`, `OTHER=2`; PDB branch is separately dispatched |
| Residue name | `std::string` | `String` | Nonempty; `ATOM_MONOMER_INFO_RESIDUENAME` |
| Residue number | `int` | `int` | Nonzero; `ATOM_MONOMER_INFO_RESIDUENUMBER` |
| Chain ID | `std::string` | `String` | Nonempty; `ATOM_MONOMER_INFO_CHAINID` |
| Monomer class | `std::string` | `String` | Nonempty; `ATOM_MONOMER_INFO_MONOMERCLASS` |

The outer `END_ATOM_MONOMER_INFO` terminates this record. The base-class fields
were extended in pickle version 16.3. The different PDB/base monomer-type wire
widths above are deliberate descriptions of the two actual writer functions,
not a proposed unified encoding.

Sources: [PDB writer/reader](../../third_party/rdkit/Code/GraphMol/MolPickler.cpp#L793),
[base monomer writer](../../third_party/rdkit/Code/GraphMol/MolPickler.cpp#L930),
[member types/defaults](../../third_party/rdkit/Code/GraphMol/MonomerInfo.h).

## 14. Alternate legacy writer and version boundaries

### 14.1 `OLD_PICKLE` / `_pickleV1` body

This is a compile-time alternate writer, **not additional fields appended to
the normal 16.3 layout**. Its caller still writes the common header; the legacy
body itself writes the following tagged atom/bond records instead of the normal
counted blocks. The present compile-time alternate path is not asserted to
produce a reader-compatible historical header/body combination.

| Level | Tag/field | Value write type | Condition |
|---|---|---|---|
| Atom | `BEGINATOM` | `tag` only | Each atom |
| Atom | `ATOM_NUMBER` | `int` atomic number | Always |
| Atom | `ATOM_INDEX` | `unsigned int` | Always; explicit in this alternate body |
| Atom position | `ATOM_POS`, `x`, `y`, `z` | Three `double` values | Selected/default conformer if present, otherwise default `Point3D` |
| Atom | `ATOM_CHARGE` | `int` | Formal charge nonzero |
| Atom | `ATOM_NEXPLICIT` | `unsigned int` | Explicit H count nonzero |
| Atom | `ATOM_CHIRALTAG` | Native `Atom::ChiralType` enum via generic `sizeof(T)` write | Nonzero enum |
| Atom | `ATOM_ISAROMATIC` | `char` | True |
| Atom | `ENDATOM` | `tag` only | Each atom |
| Bond | `BEGINBOND` | `tag` only | Each bond |
| Bond | `BOND_INDEX` | `unsigned int` | Always |
| Bond | `BOND_BEGATOMIDX` | `unsigned int` | Always |
| Bond | `BOND_ENDATOMIDX` | `unsigned int` | Always |
| Bond | `BOND_TYPE` | Native `Bond::BondType` enum via generic `sizeof(T)` write | Always |
| Bond | `BOND_DIR` | Native `Bond::BondDir` enum via generic `sizeof(T)` write | Nonzero enum |
| Bond | `ENDBOND` | `tag` only | Each bond |
| Molecule | `ENDMOL` | `tag` only | End of body |

### 14.2 Do not infer old schemas from current declarations

- The reader uses 32-bit tags for versions below 7.0; the normal current tag
  writer uses one byte. The initial header `VERSION` is a separate exception.
- Reader compatibility paths before version 14.0 use `unsigned int` property
  counts in the corresponding tables; the current writer uses `uint16_t` and
  the compact atom/bond property packs.
- Current framed property/conformer blocks are read with byte lengths from
  version 13.0 onward; the existence of older unframed paths does not change
  the current layout above.
- `ATOM_MASS`, older inline atom coordinates, `BEGINQUERYATOMDATA` and
  `QUERY_BOOL` compatibility/declaration paths must not be advertised as extra
  current writer fields. `ENDSSSR` is not emitted after current counted ring
  rows. `INVALID_TAG=255` is a rejection sentinel, not serialized data.
- Any complete migration of **historical reading** needs its own per-version
  reader inventory. This document enumerates current writes and the explicit
  alternate legacy writer, not every byte schema from every RDKit release.

Sources: [alternate writer](../../third_party/rdkit/Code/GraphMol/MolPickler.cpp#L2596),
[tag declarations](../../third_party/rdkit/Code/GraphMol/MolPickler.h#L83),
[normal reader](../../third_party/rdkit/Code/GraphMol/MolPickler.cpp#L1348).

## 15. What this pickle does not independently store

| State | Actual treatment |
|---|---|
| Atom/bond indices in normal records | Implied by record order; references use indices, not stable UUIDs |
| Adjacency structure | Reconstructed from bond endpoint fields; no separate adjacency cache dump |
| C++ owner pointers, locks, allocations, caches' object identities | Not serialized as fields |
| Entire cache initialization/invalidation state | Not one preserved validity structure; atom valences and RingInfo have the specific limited encodings above |
| Computed-property classification | Owner-local `__computedProps` value when selected, not a second per-property wire flag |
| Mass | No independent current atom mass field; isotope/atomic number are written |
| Ring bond rows, membership/fusion/ring-family caches | Not separate current writer payloads |
| SGroup owner pointer and validity flag | Not separate current writer payloads |
| Stereo-group read ID | Not separately written alongside assigned write ID |
| Original file's complete text/records | Only explicitly modeled fields/properties; no automatic verbatim SDF/PDB/CIF capture |
| Fingerprint, descriptors, arbitrary user objects | Only if stored as selected serializable properties; no automatic computation or implicit fixed field |
| Arbitrary query implementation/function pointers | Supported subtype/description/operand tree only; unknown subtypes fail |
| Reaction templates/execution state | Not part of this molecule pickle inventory |

## 16. Source coverage and reproducible identities

The field inventory follows every `streamWrite`, `streamWriteProps`,
`streamWriteProp`, raw payload write and framing helper reached from the normal
writer: header, `_pickle`, `_pickleAtomData`, `_pickleAtom`, `_pickleBond`,
`_pickleConformer`, `_pickleSSSR`, `_pickleSubstanceGroup`, `_pickleStereo`, both
monomer writers, generic/compact property writers, query classification/node
writer, and the built-in ExplicitBitVect handler and nested writer. `_pickleV1`
is accounted for separately. This is source inspection, not a runtime parity
test or a claim of format support in COSMolKit.

The following SHA-256 values identify the inspected local files; they are not
claims about a Git commit or the latest upstream release. Paths are relative
to `third_party/rdkit/Code/`.

| Source | SHA-256 |
|---|---|
| `GraphMol/MolPickler.cpp` | `48d233675fd3f8c4b62120da5d7da4c201ed3bad50fa83c9fdfe4e1e622b50c0` |
| `GraphMol/MolPickler.h` | `3bdfeb72dfac983ec3640d415d289d886bb8710c71d05043ee8d7ae5c4c8dbd2` |
| `RDGeneral/StreamOps.h` | `c8f1c61b2b8f82a81779eaa35794ac2cd293b3870a5036855ce8c601178d2377` |
| `RDGeneral/RDProps.h` | `08a12ae1d2dba8aae46d290a24d2f25e8b741d908ccb68045c72007f2c4eb270` |
| `RDGeneral/types.h` | `24e361adf251ee265c012a1f64d203342f28a24c863d529a76aecf98381b8e30` |
| `GraphMol/Atom.h` | `3b8c91e5c86004fe12e85c2d56056a5ae809cebfcaf20807b2bbe81684a4b4ea` |
| `GraphMol/Bond.h` | `15f873f32dac816103b77a28e6f45e651deffca89637fc8b869f51bf4f1c35e1` |
| `GraphMol/SubstanceGroup.h` | `1c26d4fcb14779c12a5eb8e4b2280ffc74225f5991a26ff72cf95c19fcc8502a` |
| `GraphMol/StereoGroup.h` | `2deedd35cd3fa6931deb33ad97f6d7e5ae2da64779e20076c7fc136356c0b0e1` |
| `GraphMol/MonomerInfo.h` | `ce2291eb6f634da8cd59ba86c945ea7db9388adcd7531b7a2131a4937ebe6253` |
| `DataStructs/DatastructsStreamOps.h` | `b04ec36cadbe395ea918aa3aab4bd42934dd66cae7e1e84f32828db9db8206f6` |
| `DataStructs/ExplicitBitVect.cpp` | `4965982db56c59ae884efd0257cb0df57fd8e9f7710b204e84db5d358fec3af3` |
| `DataStructs/BitVect.h` | `8702f119a5e7c92406a93c6514a14b24fe5b3294d22ea4e0f23079ced696177a` |
