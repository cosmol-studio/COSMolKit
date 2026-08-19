# Avalon Fingerprint Source Inventory

Status: Step 55 audit complete. This is the source boundary for the selected
RDKit `GetAvalonFP(ROMol, ...)` explicit-bit-vector path. It is not a port and
does not change the current fail-closed behavior.

## Top-level source path

The active adapter call graph is:

```text
pyAvalonTools.cpp:158-165
  GetAvalonFP(ROMol, nBits, isQuery, resetVect, bitFlags)
AvalonTools.cpp:204-215
  warn on nBits % 8; nBytes = nBits / 8
  molToReaccs(ROMol)
AvalonTools.cpp:147-153
  RDKit MolToMolBlock(mol, true); LocaleSwitcher; MolStr2Mol()
reaccsio.c:1428-1448
  FortranStringOpen(); ReadREACCSMolecule(); FortranClose()
AvalonTools.cpp:63-77
  round internal nBytes up to a multiple of four; SetFingerprintBits twice
ssmatch.c:1743-1765
  allocate count array, call SetFingerprintCountsWithFocus, pack positive counts
ssmatch.c:1815-4232
  CountFingerprintPatterns: preprocessing, all selected feature families, cleanup
hashcode.c:62-76
  next_hash() and hash_position()
```

The source-level conversion boundary is therefore not a SMILES parser. It is
the exact RDKit MolBlock writer followed by the Avalon REACCS MolBlock parser.
The corresponding RDKit writer declarations are in
`third_party/rdkit/Code/GraphMol/FileParsers/FileWriters.h:40-91`; the active
writer implementation is `MolFileWriter.cpp:1247-1450`.

## Adapter and packing ranges

| Source | Range | Responsibility |
|---|---:|---|
| `AvalonTools.cpp` | 48-61 | count-buffer allocation and two-pass count call |
| `AvalonTools.cpp` | 63-77 | byte-buffer allocation, four-byte rounding, two-pass bit call |
| `AvalonTools.cpp` | 79-136 | reset, explicit bit-vector packing, little-endian word packing |
| `AvalonTools.cpp` | 147-153 | RDKit MolBlock conversion, locale guard, REACCS parse |
| `AvalonTools.cpp` | 195-215 | molecule count/bit public overloads and nBits/8 conversion |
| `pyAvalonTools.cpp` | 40-58 | Python ROMol overload and words overload |
| `pyAvalonTools.cpp` | 157-165 | Python `GetAvalonFP` defaults and return allocation |
| `pyAvalonTools.cpp` | 193-201 | Python `GetAvalonFPAsWords` defaults |
| `reaccsio.c` | 1428-1448 | MolBlock string to `reaccs_molecule_t` |
| `reaccsio.c` | 709-1300 | complete REACCS molecule/atom/bond/query field parser |
| `forio.c` | 153-198 | locale-independent in-memory Fortran stream and close |

The Python default profile is `nBits=512`, `isQuery=false`,
`resetVect=false`, `bitFlags=0xF07FFF`. The C++ header default is a separate
overload contract (`resetVect=true`, `bitFlags=0x007FFF`) and must not be
silently substituted for the Python API.

## Shared fingerprint engine

### Hash and recursive path helpers

All source bodies below are active for the selected bit-vector call graph and
must receive verbatim C comments with two-axis markers when ported:

| Source | Range | Responsibility |
|---|---:|---|
| `ssmatch.c` | 1223-1296 | query hydrogen inference |
| `ssmatch.c` | 1298-1314 | family seeds |
| `ssmatch.c` | 1321-1326 | `NEXT_SEED`, modulo position and byte-bit macros |
| `ssmatch.c` | 1347-1379 | path-length matrix recursion |
| `ssmatch.c` | 1382-1424 | special-neighbour recursion |
| `ssmatch.c` | 1458-1580 | path traversal and ring closure hashing |
| `ssmatch.c` | 1595-1680 | typed feature-pair hashing and count expansion |
| `ssmatch.c` | 1682-1734 | ring-size-pair hashing |
| `hashcode.c` | 62-76 | 64-bit `next_hash` and slot hashing |
| `hashcode.c` | 82-88 | shortcut-label string hashing |

The old `__OLDFP` branch in `hashcode.c:50-60` is not selected by the pinned
RDKit build. It remains source context and must not be accidentally used as a
compatibility fallback.

### Preprocessing and state mutation

`CountFingerprintPatterns` is `ssmatch.c:1815-4232`. Its ordered setup is:

1. zero or preserve the output count array according to `ACCUMULATE_BITS`;
2. allocate and fill `neighbourhood_t` using `utilities.c:2423-2467`;
3. snapshot bond types, allocate touched/path state, and calculate H counts;
4. for query molecules call `GuessHCountsFromSubstitution` and count explicit
   H/D/T bonds; for ordinary molecules call `ComputeImplicitH` then add explicit
   hydrogen bonds (`utilities.c:2774-2851`);
5. use `PerceiveDYAromaticity` (`perceive.c:398-662`) when
   `USE_DY_AROMATICITY` is set, otherwise `PerceiveAromaticBonds`
   (`perceive.c:263-382`);
6. compute `RingState` (`utilities.c:2905-2948`) and `SetRingSizeFlags`
   (`perceive.c:141-178`);
7. assign atom colors, shortcut-label colors, bond colors, degree, carbon
   degree, unsaturation and ring-fusion counts.

The engine mutates the REACCS working molecule while doing this. At cleanup
(`ssmatch.c:4213-4224`) it restores bond types and frees working arrays. The
port must reproduce this source mutation boundary and restoration behavior;
it must not assume that a generic immutable traversal is equivalent.

## Flag-family line ledger

The `which_bits` mask is defined in `include/ssmatch.h:127-157`. The active
branches in `CountFingerprintPatterns` are:

| Flag family | Mask | Source range | Prerequisites / notes |
|---|---:|---:|---|
| `USE_ATOM_COUNT` | `0x000010` | `2038-2210` | atom colors, H counts, ring/fusion counters |
| `USE_ATOM_SYMBOL_PATH` | `0x000004` | `2215-2403` | path matrix and atom color seeds |
| `USE_AUGMENTED_ATOM` | `0x000020` | `2406-2534` | atom colors, degree, H/ring state |
| `USE_AUGMENTED_BOND` | `0x000400` | `2537-2608` | bond colors, endpoint atom colors |
| `USE_HCOUNT_PAIR` | `0x000100` | `2611-2668` | H counts, bond and atom colors |
| `USE_HCOUNT_PATH` | `0x000040` | `2671-2805` | path matrix, H counts |
| `USE_RING_PATH` | `0x000002` | `2808-2870` | ring status and path matrix |
| `USE_BOND_PATH` | `0x000200` | `2872-3030` | path matrix, bond colors |
| `USE_HCOUNT_CLASS_PATH` | `0x000080` | `3032-3073` | path matrix, H-count classes |
| `USE_ATOM_CLASS_PATH` | `0x000008` | `3076-3307` | path matrix, atom classes |
| `USE_RING_PATTERN` | `0x000001` | `3310-3414` | ring status, ring bond types, closure recursion |
| `USE_RING_SIZE_COUNTS` | `0x000800` | `3418-3523` | ring-size flags and ring counters |
| `USE_DEGREE_PATH` | `0x001000` | `3527-3625` | degree arrays and path matrix |
| `USE_CLASS_SPIDERS` | `0x002000` | `3628-4021` | distance matrix, ring classes, `SetFeatureBits` |
| `USE_FEATURE_PAIRS` | `0x004000` | `3628-4021` | distance matrix, ring classes, `SetFeatureBits` |
| `USE_SCAFFOLD_IDS` | `0x100000` | `4025-4204` | non-query only; ring ext-connectivity |
| `USE_SCAFFOLD_COLORS` | `0x200000` | `4025-4204` | non-query only; ring ext-connectivity |
| `USE_SCAFFOLD_LINKS` | `0x400000` | `4025-4204` | non-query only; ring distance matrix |
| `USE_SHORTCUT_LABELS` | `0x800000` | setup `1968-1971`, downstream colors | `R` atom text and `hash_string` |

`USE_ALL_FEATURES` is `0x007FFF`; `USE_NON_SSS_BITS` is `0xF00000`.
The non-SSS branch is gated by `!as_query` (`4031-4033`), so query mode is a
behavioral branch rather than merely an option annotation.

## Public option mapping

| COSMolKit public field | Pinned source field | Mapping status |
|---|---|---|
| `n_bits` | Python `nBits` / C++ `nBits` | source-backed; preserve integer conversion and byte boundary behavior |
| `is_query` | `isQuery` | source-backed; controls H inference, second pass, and non-SSS gating |
| `bit_flags` | `bitFlags` | source-backed typed mask; all defined flag families remain addressable |
| fresh result allocation | wrapper `new ExplicitBitVect(nBits)` | `resetVect` is internal, not a public COSMolKit value option |
| old `min_path`, `max_path`, `n_bits_per_hash`, `use_bond_order`, `use_hs`, `tautomeric_fingerprint`, `from_atoms` | none | invented placeholder fields; remove before enabling support |

## Error and cleanup boundary

`molToReaccs` uses a postcondition after `MolStr2Mol`; parser failure is a
structured source failure in the adapter boundary, not permission to emit an
empty fingerprint. `SetFingerprintBits` returns zero for a null REACCS
molecule and otherwise owns a temporary `int` count array. The port must
return a structured COSMolKit error for conversion or unsupported molecule
states and must not silently fall back to another fingerprint family.

## Next source-backed units

The ordered plan now proceeds through public type correction and focused
conversion tests before copying the `MolToReaccs`/REACCS parser and then the
preprocessing/traversal body. No flag branch may be marked behaviorally
complete until its complete helper closure and exact-bit tests are present.
