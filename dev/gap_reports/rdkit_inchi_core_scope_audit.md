# RDKit InChI Core Scope Audit

## Scope Decision

This audit replaces the all-export `libinchi` port objective with the exact engine and adapter closure required by four scalar RDKit APIs:

- `Chem.MolToInchi(mol)`;
- `Chem.MolToInchiKey(mol)`;
- `InchiToInchiKey(inchi)`;
- `Chem.MolFromInchi(inchi)`.

MolBlock, SDF, V3000 Molfile, IXA, AuxInfo conversion, incremental INCHIGEN, version-query, CLI/demo, and unrelated export parity are outside this scope. Already ported outside-scope code remains private and frozen; this plan neither deletes it nor spends further parity work on it.

## Compiler Evidence

- Official target: GCC/Linux `libinchi` with `COMPILE_ANSI_ONLY`, `TARGET_API_LIB`, and `fPIC` from the official CMake target.
- Definition and direct-edge evidence: the configured 60-entry `compile_commands.json`, GCC `-aux-info`, and GCC `-fdump-ipa-cgraph` artifacts.
- Callback evidence: active GCC IPA `References: ... (addr)` records from selected callers; no function-name-only inference is used.
- Active configured definitions: `1364`.
- Selected official roots: `5`.
- Compiler direct-edge selected closure: `1014` definitions.
- Callback-complete selected closure: `1050` definitions.
- Compiler address-taken callback edges in the selected closure: `57`.
- Address-taken callback target definitions added to the direct closure: `36`.
- Selected source frames present: `1050`.
- Behavior-closed selected source locations: `1049`.
- Selected source locations with an open first axis: `1` (`NormalizeAndCompare`).
- Remaining source-implementation Port units: `0`; remaining official exact-parity closure units: `1`.
- Active completed functions outside scope and frozen: `108`.
- Active functions outside scope and not scheduled: `206`.

## Root Mapping

| RDKit API | Official roots | RDKit evidence |
|---|---|---|
| `Chem.MolToInchi` | `GetINCHI`, `FreeINCHI` | `inchi.cpp:2082`, `inchi.cpp:2100` |
| `Chem.MolToInchiKey` | `GetINCHI`, `FreeINCHI`, `GetINCHIKeyFromINCHI` | `inchi.h:107` composition |
| `InchiToInchiKey` | `GetINCHIKeyFromINCHI` | `inchi.cpp:2149` |
| `Chem.MolFromInchi` | `GetStructFromINCHI`, `FreeStructFromINCHI` | `inchi.cpp:1273`, `inchi.cpp:1371`, `inchi.cpp:1648` |

## State Constraints

- `GetINCHI` wraps plain `inchi_Input`, sets `extended_input.polymer = NULL` and `extended_input.v3000 = NULL`, then calls `GetINCHI1`.
- `GetStructFromINCHI` uses the InChI string path and forces `INPUT_INCHI`; it does not consume a Molfile.
- User option strings remain part of supported RDKit behavior, so their reachable normalization, canonicalization, reconstruction, stereo, FixedH, reconnected, warning, and error paths are retained.
- Global `UserAction` and `ConsoleQuit` callbacks are initialized to null and have no configured RDKit API setter; no invented callback target is added for them.
- Macro-only allocation behavior remains attached to callers and does not create fake functions for inactive `util.c` allocator bodies.
- The abandoned import-only staging for the unimplemented `ReadMolfileToInpAtoms` step is not a completed port and is scheduled for removal before the first remaining core function; no completed Molfile/V3000 implementation is removed.

## Step 5269 Closure Reaudit

The source-frame marker scan found two stale/open first-axis summaries before this reaudit. `MolToInchiKey` was a stale marker: its focused Rust test and pinned RDKit C++ oracle had already passed exactly and its audit row already recorded `RDKit✔️❌`. The source-frame marker is now corrected to `RDKit✔️❌`, preserving the reviewed performance axis.

`NormalizeAndCompare` is different and remains open. The configured compiler direct-edge path from a selected root is:

```text
GetStructFromINCHI
  -> GetStructFromINCHIEx
  -> ReadWriteInChI
  -> ConvertInChI2Struct
  -> AllInchiToStructure
  -> InChI2Atom
  -> OneInChI2Atom
  -> RestoreAtomMakeBNS
  -> RunBnsRestore1
  -> NormalizeAndCompare
```

Every edge above is present in the GCC `-fdump-ipa-cgraph`-derived direct-edge table in `official_inchi_active_call_graph_audit.md`; no function-name inference is used. The function therefore belongs to the narrowed five-root selected closure and cannot be frozen as outside scope.

| Function | Source frame | Focused Rust test | Official C oracle | Exact result | Marker | Unclosed branches |
|---|---|---|---|---|---|---|
| `MolToInchiKey` | `third_party/rdkit/External/INCHI-API/inchi.h:107-110`; complete verbatim header-defined inline frame | `cargo test -p cosmolkit-inchi source_port__inchi__moltoinchikey__line_107 -- --nocapture` (`1` main-harness test executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__moltoinchikey__exact -- --nocapture` (`1` main-harness test executed; `9` complete pinned C++ records) | exact pass for every field recorded in the existing row below | `RDKit✔️❌` | none within the modeled toolkit-neutral composition boundary |
| `NormalizeAndCompare` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:1274-1672`; `diff -u` against the in-function Rust comment frame returned zero differences | `cargo test -p cosmolkit-inchi source_port__ichirvr4__normalizeandcompare__line_1274 -- --nocapture` (`1` main-harness test executed and passed; the separate provenance harness executed `0` and is not counted) | `cargo test -p cosmolkit-inchi official_c_oracle__normalizeandcompare__exact -- --nocapture` (`1` main-harness test executed; `241` provenance-checked official-C records; zero field mismatches) | exact pass for all 240 source-defined records; the separate undefined-path evidence record is reproduced as process-termination evidence and is not treated as a defined return/cleanup result | `INCHI❗❌` | the active initial `inchi_strbuf_init` allocation-failure path reaches official-C undefined behavior before a source-defined return or cleanup sequence; see the blocker below |

The aggregate scalar RDKit test does not close this internal source function: it verifies selected end-to-end adapter records, not every observable `NormalizeAndCompare` branch. The dedicated official-C oracle now executes and passes every compared field for all source-defined records, but public API exposure remains blocked by the active undefined allocation path described below.

## Step 5271 NormalizeAndCompare Branch And Field Matrix

### Active Configuration

- Source frame: `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:1274-1672`; the in-function Rust comment frame has zero `diff -u` differences from these 399 official lines.
- Target: GCC/Linux `libinchi`, `COMPILE_ANSI_ONLY=1`, `TARGET_API_LIB=1`, `bRELEASE_VERSION=1`, `FIX_ADD_PROTON_FOR_ADP=0`.
- Inactive and excluded from production cases: both `bRELEASE_VERSION == 0 && !TARGET_API_LIB` diagnostic `fprintf` blocks and the complete `FIX_ADD_PROTON_FOR_ADP == 1` block.
- The source has no direct diagnostic stream output in the selected configuration. Diagnostics produced inside callees remain observable through caller-owned `STRUCT_DATA` and callee-mutated state.

### Required Branch Cases

| Case family | Active source decisions that must be distinguished exactly |
|---|---|
| initial rebuild | `MakeOneInChIOutOfStrFromINChI2` returns negative or nonnegative; negative return must still enter the common cleanup block |
| original/reversed layer selection | `bMobileH == TAUT_NON` and the non-fixed-H alternative; independently for original and reversed layer 1: null pointer, zero atoms, deleted, and live nondeleted atoms |
| Zz/Zy gate | each false operand of `n_zy && n_pzz`, true gate with null formula, true gate with ordinary formula, and true gate with mergeable `Zz` formula |
| Zz/Zy string buffer | `inchi_strbuf_init > 0` and `<= 0`; `MergeZzInHillFormula` success and allocation failure; source-reachable non-growth and source-unreachable growth predicate; `inchi_realloc` null/non-null if a controlled direct-callee oracle forces growth; close after every defined path |
| first fixed-H enrichment | `iMobileH != TAUT_NON`, `FillOutExtraFixedHDataRestr == 0`, and every nonzero signed return with immediate cleanup |
| first comparison | `IDIF_PROBLEM` precedence, nonzero `err`, and ordinary flags; serialize the full `ICR` and exact direct-callee call arguments |
| less-H repair | each short-circuit operand false; positive delta; fixer negative, zero, and positive; rebuild negative; repeated fixed-H enrichment nonzero; loop exits for cleared flag, null prepared atoms, zero delta, nondecreasing delta, and continues for strictly decreasing positive delta |
| more-H repair | each short-circuit operand false; positive delta; fixer negative, zero, and positive; rebuild negative; repeated fixed-H enrichment nonzero; loop exits for cleared flag, null prepared atoms, zero delta, nondecreasing delta, and continues for strictly decreasing positive delta |
| extra-endpoint repair | each short-circuit operand false; positive endpoint delta; fixer negative, zero, and positive; rebuild negative; repeated fixed-H enrichment nonzero; loop exits for cleared flag, null normalized atoms, zero delta, nondecreasing delta, and continues for strictly decreasing positive delta |
| fixed-H restoration | `bMobileH == TAUT_NON`; fixer negative, zero, and positive; exact post-increment loop behavior for one, two, and three positive returns before exit |
| mobile-H restoration | `bMobileH == TAUT_YES`; fixer negative, zero, and positive; a source-defined non-enum value that takes neither fixed nor mobile branch |
| final primary comparison | `IDIF_PROBLEM` precedence, nonzero `err`, and ordinary flags |
| optional layer comparison | both indexes zero versus either nonzero; ordinary result and nonzero `err`; preserve the official line-1576 behavior that tests primary `cmpInChI`, not `cmpInChI2`, for `IDIF_PROBLEM` |
| stereo repair | negative, zero, and positive `FixRestoredStructureStereo` returns |
| endpoint counting | zero/nonzero `pTCGroups->num_tgroups`; each `num_edges` contribution and source-defined signed sum boundary; null/non-null `t_group_info`; zero/nonzero `num_t_groups`; each `num[0]` false/true branch and `nNumEndpoints` contribution |
| common cleanup | for both `TAUT_NUM` slots, every null/non-null combination of `pOneINChI`, `pOneINChI_Aux`, and `pOne_norm_data`; exact `Free_INChI`, `Free_INChI_Aux`, `FreeInpAtomData`, holder `inchi_free`, pointer-null assignment, and final `free_t_group_info` order after success and every defined early return |
| legal integer boundaries | negative/zero/positive helper returns, `LONG_MIN`/`LONG_MAX` `num_inp` under inactive logging, wrapping-sensitive Rust deltas only where the corresponding C subtraction is defined, valid zero/maximum representable counts that do not cause C signed overflow, and exact `num_tries++ < 2` ordering |

### Required Observable Record

Every oracle record must contain the case id, active macro tuple, direct-callee script, exact call trace with argument-layer indexes, return value, `pnNumRunBNS`, `pnTotalDelta`, complete caller-owned `StrFromINChI` scalar state, both cleaned pointer slots, pre-free snapshots of each pointed `INChI`, `INChI_Aux`, and `INP_ATOM_DATA`, `One_ti` before and after cleanup, complete `BN_STRUCT` vertices/edges/totals, `BN_DATA`, `at`/`at2`/`at3`, `VAL_AT`, `ALL_TC_GROUPS` and `TC_GROUP` arrays, preservation of both input `pInChI` objects, formula bytes before cleanup, allocation/free event order, outstanding allocation count, and process termination status. Rust must construct the same record independently and compare the complete JSON value; shared expected constants or omitted fields are not permitted.

### Undefined Official Allocation Path

The Zz/Zy branch contains an active official-C undefined-behavior path:

1. `inchi_strbuf_init` zeroes the complete `INCHI_IOS_STRING` and returns `-1` when its `inchi_calloc` fails (`ichi_io.c:1370-1395`).
2. `NormalizeAndCompare` does not branch on that negative result except to skip `inchi_strbuf_printf`.
3. `MergeZzInHillFormula` returns `0` immediately because `strbuf->pStr == NULL` (`ichiprt1.c:5435-5443`).
4. Official line 1359 then executes `strcpy(existing_formula, strbuf->pStr)` with a null source pointer.

The C language defines no return value, cleanup sequence, mutation state, or portable signal for this execution. A pinned GCC build may terminate with a particular signal, but that observation is not source-defined exact behavior and cannot justify a behavioral `✔️` marker. The Rust implementation currently returns `SourceHeapError::AllocationFailed` before the common cleanup on this path, which is also not exact C behavior. Because Step 5273 requires allocation-failure and cleanup-order exact parity and forbids skipping active exceptional input, this branch is a concrete blocker. `NormalizeAndCompare` remains `INCHI❗❌`; no public API exposure is permitted under the current completion rule.

### Step 5297 Zero-Mismatch Artifact

- Command: `cargo test -p cosmolkit-inchi official_c_oracle__normalizeandcompare__exact -- --nocapture`.
- Main harness: `running 1 test`; `official_c_oracle__normalizeandcompare__exact ... ok`; `1 passed; 0 failed`.
- Oracle source: pinned official GCC/Linux `libinchi` with `COMPILE_ANSI_ONLY` and `TARGET_API_LIB`; the source-tree checksum and instrumented `ichirvr4.c` token restoration are verified before execution.
- Record count: `241`, comprising `128` null/cleanup records, `32` layer-selection records, `1` common-success record, `9` defined Zz/Zy records, `1` undefined-path evidence record, `18` first-less records, `26` more/extra records, and `26` final restoration/comparison/stereo/endpoint records.
- Exact result: zero JSON field mismatches across all records. The final family includes `LONG_MIN`, `LONG_MAX`, Fixed-H negative/zero and one/two/three positive-return iterations, mobile-H negative/zero/positive, non-enum `bMobileH`, final primary precedence/error/ordinary flags, skipped/original-only/reversed-only/both optional selection, the official line-1576 primary-flag check, stereo negative/zero/positive, zero/nonzero TC groups, null/non-null `t_group_info`, false/true `num[0]`, signed endpoint-count boundaries, and early/success cleanup order.
- Marker decision: preserve `INCHI❗❌`. The exact source-defined records require no Rust repair, but the undefined allocation path above remains unclosed and prohibits a behavioral-axis upgrade.

Step 5301 repeated the complete exact oracle after this artifact update. The main harness again executed one matching test and all `241` records matched with zero field mismatches, so no production Rust repair was indicated. The first axis remains `❗` and the independently reviewed second axis remains `❌` because the undefined path is unchanged.

### Step 5307 Final NormalizeAndCompare Closure Decision

- Source frame: `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:1274-1672`; the complete active function remains copied verbatim in `crates/cosmolkit-inchi/src/source/base/ichirvr4.rs`.
- Focused Rust command: `cargo test -p cosmolkit-inchi source_port__ichirvr4__normalizeandcompare__line_1274 -- --nocapture`; the main harness executed `1` matching test and passed with zero failures.
- Official C oracle command: `cargo test -p cosmolkit-inchi official_c_oracle__normalizeandcompare__exact -- --nocapture`; the main harness executed `1` matching test and passed with zero failures.
- Exact source-defined result: all `240` records for which official C defines observable return, mutation, allocation, cleanup, integer, comparison, and call-order behavior match Rust field-for-field with zero mismatches.
- Undefined-path evidence: the remaining record forces initial `inchi_strbuf_init` allocation failure and confirms that official C reaches `strcpy(existing_formula, NULL)`; it is process-termination evidence only because C defines no exact return, signal, mutation, or cleanup sequence for this path.
- Final marker: preserve `INCHI❗❌`; the first axis cannot be upgraded and the independently reviewed second-axis performance gap remains unchanged.
- Exposure decision: the selected call graph retains one open active branch, so the four scalar public APIs remain blocked under the current completion rule.

## Closed During Core Execution

| Function | Source | Focused Rust test | Official C oracle | Exact result | Marker | Unclosed branches |
|---|---|---|---|---|---|---|
| `cmp_components` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4247` | `cargo test -p cosmolkit-inchi source_port__strutil__cmp_components__line_4247` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__cmp_components__exact` (`1` executed; `10` complete official C records) | exact pass | `INCHI✔️✔️` | none |
| `cmp_rad_endpoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:7784` | `cargo test -p cosmolkit-inchi source_port__ichi_bns__cmp_rad_endpoints__line_7784` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__cmp_rad_endpoints__exact` (`1` executed; `5` complete official C records) | exact pass | `INCHI✔️✔️` | none |
| `CompTGroupNumber` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9195` | `cargo test -p cosmolkit-inchi source_port__ichi_bns__comptgroupnumber__line_9195` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__comptgroupnumber__exact` (`1` executed; `5` complete official C records) | exact pass | `INCHI✔️✔️` | none |
| `CompCGroupNumber` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9202` | `cargo test -p cosmolkit-inchi source_port__ichi_bns__compcgroupnumber__line_9202` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__compcgroupnumber__exact` (`1` executed; `5` complete official C records) | exact pass | `INCHI✔️✔️` | none |
| `cmp_iso_atw_diff_component_no` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:166` | `cargo test -p cosmolkit-inchi source_port__strutil__cmp_iso_atw_diff_component_no__line_166` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__cmp_iso_atw_diff_component_no__exact` (`1` executed; `7` complete official C records) | exact pass | `INCHI✔️✔️` | none |
| `CompNeighLists` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:632` | `cargo test -p cosmolkit-inchi source_port__ichisort__compneighlists__line_632` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__compneighlists__exact` (`1` executed; `7` complete official C records) | exact pass | `INCHI✔️❌` | none |
| `CompNeighListsUpToMaxRank` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:644` | `cargo test -p cosmolkit-inchi source_port__ichisort__compneighlistsuptomaxrank__line_644` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__compneighlistsuptomaxrank__exact` (`1` executed; `8` complete official C records) | exact pass | `INCHI✔️❌` | none |
| `cmp_charge_val` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:732` | `cargo test -p cosmolkit-inchi source_port__ichirvr1__cmp_charge_val__line_732` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__cmp_charge_val__exact` (`1` executed; `14` complete official C records) | exact pass | `INCHI✔️✔️` | none |
| `comp_cc_cand` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4584` | `cargo test -p cosmolkit-inchi source_port__ichirvr1__comp_cc_cand__line_4584` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__comp_cc_cand__exact` (`1` executed; `12` complete official C records) | exact pass | `INCHI✔️✔️` | none |
| `base26_triplet_1` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1173` | `cargo test -p cosmolkit-inchi source_port__ikey_base26__base26_triplet_1__line_1173` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__base26_triplet_1__exact` (`1` executed; `16388` complete official C records) | exact pass | `INCHI✔️✔️` | none |
| `base26_triplet_2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1191` | `cargo test -p cosmolkit-inchi source_port__ikey_base26__base26_triplet_2__line_1191` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__base26_triplet_2__exact` (`1` executed; `16388` complete official C records) | exact pass | `INCHI✔️✔️` | none |
| `base26_triplet_3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1212` | `cargo test -p cosmolkit-inchi source_port__ikey_base26__base26_triplet_3__line_1212` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__base26_triplet_3__exact` (`1` executed; `16388` complete official C records) | exact pass | `INCHI✔️✔️` | none |
| `base26_triplet_4` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1233` | `cargo test -p cosmolkit-inchi source_port__ikey_base26__base26_triplet_4__line_1233` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__base26_triplet_4__exact` (`1` executed; `16388` complete official C records) | exact pass | `INCHI✔️✔️` | none |
| `base26_dublet_for_bits_28_to_36` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1262` | `cargo test -p cosmolkit-inchi source_port__ikey_base26__base26_dublet_for_bits_28_to_36__line_1262` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__base26_dublet_for_bits_28_to_36__exact` (`1` executed; `516` complete official C records) | exact pass | `INCHI✔️✔️` | none |
| `base26_dublet_for_bits_56_to_64` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1287` | `cargo test -p cosmolkit-inchi source_port__ikey_base26__base26_dublet_for_bits_56_to_64__line_1287` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__base26_dublet_for_bits_56_to_64__exact` (`1` executed; `516` complete official C records) | exact pass | `INCHI✔️✔️` | none |
| `get_xtra_hash_major_hex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1307` | `cargo test -p cosmolkit-inchi source_port__ikey_base26__get_xtra_hash_major_hex__line_1307` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__get_xtra_hash_major_hex__exact` (`1` executed; `256` complete official C records) | exact pass | `INCHI✔️✔️` | none |
| `get_xtra_hash_minor_hex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1328` | `cargo test -p cosmolkit-inchi source_port__ikey_base26__get_xtra_hash_minor_hex__line_1328` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__get_xtra_hash_minor_hex__exact` (`1` executed; `256` complete official C records) | exact pass | `INCHI✔️✔️` | none |
| `sha2_process` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:72` | `cargo test -p cosmolkit-inchi source_port__sha2__sha2_process__line_72` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__sha2_process__exact` (`1` executed; `256` complete official C records) | exact pass | `INCHI✔️✔️` | none |
| `sha2_update` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:206` | `cargo test -p cosmolkit-inchi source_port__sha2__sha2_update__line_206` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__sha2_update__exact` (`1` executed; `36` complete official C records) | exact pass | `INCHI✔️❌` | none; the second axis remains `❌` because safe partial-buffer compression copies one fixed 64-byte block before processing |
| `sha2_finish` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:257` | `cargo test -p cosmolkit-inchi source_port__sha2__sha2_finish__line_257` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__sha2_finish__exact` (`1` executed; `12` complete official C records) | exact pass | `INCHI✔️❌` | none; the second axis retains the completed `sha2_update` fixed-block copy cost |
| `sha2_starts` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:57` | `cargo test -p cosmolkit-inchi source_port__sha2__sha2_starts__line_57` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__sha2_starts__exact` (`1` executed; `32` complete official C records) | exact pass | `INCHI✔️✔️` | none |
| `sha2_csum` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:312` | `cargo test -p cosmolkit-inchi source_port__sha2__sha2_csum__line_312` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__sha2_csum__exact` (`1` executed; `26` complete official C records) | exact pass | `INCHI✔️❌` | none; the second axis retains the completed `sha2_update` fixed-block copy cost |
| `GetINCHIKeyFromINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_dll.c:113` | `cargo test -p cosmolkit-inchi source_port__ikey_dll__getinchikeyfrominchi__line_113` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__getinchikeyfrominchi__exact` (`1` executed; `38` complete official C records) | exact pass | `INCHI✔️❌` | none; the second axis remains `❌` because checked `SourceHeap` access and byte-buffer modeling add known overhead |

## Remaining Official Port Units

None. The selected official five-root callback-complete closure contains `1050`
active definitions, and all `1050` source locations have complete source frames
and Rust implementations. This does not mean that all `1050` behavior axes are
closed: `1049` are behavior-closed and `NormalizeAndCompare` remains the single
official exact-parity closure blocker described above.

No remaining source-implementation unit belongs to an unfinished SCC. The
completed `CheckINCHI`/`GetINCHIfromINCHI` SCC remains in the checked prefix and
is not rescheduled.

## Closed During RDKit Adapter Execution

| Function | Source | Focused Rust test | Pinned RDKit C++ oracle | Exact result | Marker | Unclosed branches |
|---|---|---|---|---|---|---|
| `assignBondDirs` | `third_party/rdkit/External/INCHI-API/inchi.cpp:91` | `cargo test -p cosmolkit-inchi source_port__inchi__assignbonddirs__line_91` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__assignbonddirs__exact` (`1` executed; `13` complete pinned C++ records extracted from the provenance-checked official function) | exact pass for return, range exception, diagnostic, rule-vector preservation, all bond-direction fields, ordered-set/FIFO effects, conflicts, and partial mutation | `RDKit✔️✔️` | none within the modeled `RWMol` bond-direction boundary; atom fields and properties are not accessed by this source function |
| `findAlternatingBonds` | `third_party/rdkit/External/INCHI-API/inchi.cpp:207` | `cargo test -p cosmolkit-inchi source_port__inchi__findalternatingbonds__line_195` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__findalternatingbonds__exact` (`1` executed; `20` complete pinned C++ records extracted from the provenance-checked official function) | exact pass for target pointer/index, complete reverse path stack, complete global visited set, root reset, shortest/equal path behavior, insertion-order effects, depth and unsigned boundaries, all active bond branches, every atom/bond/stereo field, diagnostics, properties, and graph immutability | `RDKit✔️✔️` | none within the modeled valid `ROMol` graph boundary; the source exposes no diagnostic or exception path for valid atom/bond pointers |
| `getNumDoubleBondedNegativelyChargedNeighboringSi` | `third_party/rdkit/External/INCHI-API/inchi.cpp:306` | `cargo test -p cosmolkit-inchi source_port__inchi__getnumdoublebondednegativelychargedneighboringsi__line_306` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__getnumdoublebondednegativelychargedneighboringsi__exact` (`1` executed; `5` complete pinned C++ records extracted from the provenance-checked official function) | exact pass for isolated atoms, reversed bond endpoints, every predicate failure, multiple matching silicon neighbors, complete atom/bond/stereo fields, graph immutability, diagnostics, properties, and exception fields | `RDKit✔️✔️` | none within the modeled valid `ROMol` graph boundary; the source exposes no diagnostic or exception path for valid atom and adjacency entries |
| `_Valence4NCleanUp1` | `third_party/rdkit/External/INCHI-API/inchi.cpp:324` | `cargo test -p cosmolkit-inchi source_port__inchi__valence4ncleanup1__line_324` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__valence4ncleanup1__exact` (`1` executed; `9` complete pinned C++ records extracted from the provenance-checked official function) | exact pass for return and exception paths, diagnostics, eligibility predicates, zero/unique/multiple/remote substructure matches, temporary atom mutation and restoration, bond flips, complete modeled atom/bond/stereo fields, explicit hydrogens, counts, properties, and unchanged graph state | `RDKit✔️✔️` | none within the modeled `InchiToMol` graph boundary; unsupported RDKit bond types produce the source-defined `ValueErrorException` equivalent before any mutation |
| `_Valence4NCleanUp2` | `third_party/rdkit/External/INCHI-API/inchi.cpp:376` | `cargo test -p cosmolkit-inchi source_port__inchi__valence4ncleanup2__line_376` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__valence4ncleanup2__exact` (`1` executed; `6` complete pinned C++ records extracted from the provenance-checked official function) | exact pass for null-target and selected-target returns, every direct-target predicate, reversed endpoints, first-candidate adjacency ordering, miss-before-hit traversal, graph mutation order, complete modeled atom/bond/stereo fields, explicit hydrogens, counts, diagnostics, exceptions, properties, and unchanged graph state | `RDKit✔️✔️` | none within the modeled valid `InchiToMol` `RWMol` boundary; the source exposes no exception path for valid atom and graph pointers |
| `_Valence5NCleanUp1` | `third_party/rdkit/External/INCHI-API/inchi.cpp:393` | `cargo test -p cosmolkit-inchi source_port__inchi__valence5ncleanup1__line_393` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__valence5ncleanup1__exact` (`1` executed; `5` complete pinned C++ records extracted from the provenance-checked official function) | exact pass for null-target return, reversed endpoints and field preservation, alternating-path stack order and bond flips, source depth limit, target-charge mutation before explicit-valence exception, complete modeled atom/bond/stereo fields, explicit hydrogens, diagnostics, properties, and graph mutation order | `RDKit✔️✔️` | none within the modeled valid `InchiToMol` `RWMol` boundary; unsupported RDKit bond types produce the source-defined `ValueErrorException` equivalent after the source-defined target-charge mutation |
| `_Valence5NCleanUp2` | `third_party/rdkit/External/INCHI-API/inchi.cpp:417` | `cargo test -p cosmolkit-inchi source_port__inchi__valence5ncleanup2__line_417` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__valence5ncleanup2__exact` (`1` executed; `6` complete pinned C++ records extracted from the provenance-checked official function) | exact pass for null-target and depth-limit returns, both root-bond endpoint branches, equal-length target adjacency ordering, triple/single bond mutation order, intermediate and target charge updates, target-before-root explicit-valence exceptions, complete modeled atom/bond/stereo fields, explicit hydrogens, diagnostics, properties, and partial mutation | `RDKit✔️✔️` | none within the modeled valid `InchiToMol` `RWMol` boundary; unsupported RDKit bond types produce the source-defined `ValueErrorException` equivalent after all source-defined graph mutations, at the target or root check in source order |
| `_Valence5NCleanUp3` | `third_party/rdkit/External/INCHI-API/inchi.cpp:443` | `cargo test -p cosmolkit-inchi source_port__inchi__valence5ncleanup3__line_443` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__valence5ncleanup3__exact` (`1` executed; `5` complete pinned C++ records extracted from the provenance-checked official function) | exact pass for null-target return, neutral-oxygen guard without mutation, charged-oxygen mutation path, target charge before target explicit-valence exception, bond and both charge mutations before root explicit-valence exception, complete modeled atom/bond/stereo fields, explicit hydrogens, diagnostics, properties, and partial mutation | `RDKit✔️✔️` | none within the modeled valid `InchiToMol` `RWMol` boundary; unsupported RDKit bond types produce the source-defined `ValueErrorException` equivalent after the source-defined target charge mutation or after the source-defined bond and root charge mutations |
| `_Valence5NCleanUp4` | `third_party/rdkit/External/INCHI-API/inchi.cpp:477` | `cargo test -p cosmolkit-inchi source_port__inchi__valence5ncleanup4__line_477` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__valence5ncleanup4__exact` (`1` executed; `5` complete pinned C++ records extracted from the provenance-checked official function) | exact pass for zero and one match returns, all three silicon predicates, exactly two mutations in adjacency order, reversed endpoints, third-match early return without mutation, late nonmatch, complete modeled atom/bond/stereo fields, explicit hydrogens, diagnostics, properties, and graph preservation | `RDKit✔️✔️` | none within the modeled valid `InchiToMol` `RWMol` boundary; the source `#if 0` sulfur/carbon block is preprocessor-inactive in the pinned GCC/Linux build and no active exception path exists for valid graph pointers |
| `_Valence5NCleanUp5` | `third_party/rdkit/External/INCHI-API/inchi.cpp:544` | `cargo test -p cosmolkit-inchi source_port__inchi__valence5ncleanup5__line_544` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__valence5ncleanup5__exact` (`1` executed; `14` complete pinned C++ records extracted from the provenance-checked official function) | exact pass for invalid atomic-number precondition diagnostics, O/S/F/Cl searches, absent/neutral/charged/both targets, shared visited reset, reversed endpoints, depth 7 and over-depth behavior, selected stack mutation, hydrogen transfer, charge mutation, charged-before-neutral explicit-valence exception order, partial mutation, complete modeled atom/bond/stereo fields, counts, diagnostics, properties, and graph preservation | `RDKit✔️✔️` | none within the modeled valid `InchiToMol` `RWMol` boundary; `CHECK_INVARIANT` is source-unreachable after a successful +1-charge target search, while active precondition and value-exception paths are exact |
| `_Valence5NCleanUp6` | `third_party/rdkit/External/INCHI-API/inchi.cpp:625` | `cargo test -p cosmolkit-inchi source_port__inchi__valence5ncleanup6__line_625` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__valence5ncleanup6__exact` (`1` executed; `40` complete pinned C++ records extracted from the provenance-checked official function) | exact pass for atomic-number and charge short-circuiting, entry explicit-valence exception before mutation, temporary nitrogen-to-tin mutation, zero/one/multiple default atom-set-unique whole-molecule substructure matches, unrelated unique tin match, every modeled wildcard bond type, atom and bond predicate misses, restoration on non-unique match, exact four-bond mutation order, final atomic-number restoration and charge update, complete modeled atom/bond/stereo fields, explicit hydrogens, diagnostics, properties, graph preservation, and exception state | `RDKit✔️✔️` | none within the modeled valid `InchiToMol` `RWMol` boundary; default `SubstructMatch` searches the complete molecule and `Bond::UNSPECIFIED` accepts every modeled RDKit bond type exactly as exercised by the pinned oracle |
| `_Valence5NCleanUp7` | `third_party/rdkit/External/INCHI-API/inchi.cpp:679` | `cargo test -p cosmolkit-inchi source_port__inchi__valence5ncleanup7__line_679` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__valence5ncleanup7__exact` (`1` executed; `39` complete pinned C++ records extracted from the provenance-checked official function) | exact pass for alternating-path search before all atom preconditions, null-target and source depth-limit returns, atomic-number and charge short-circuiting, explicit-valence mismatch and exception order, temporary nitrogen-to-tin mutation and restoration, zero/one/multiple atom-set-unique whole-molecule matches, unrelated unique match excluding both argument and target, every modeled wildcard bond type, atom and nonpath bond predicate misses, query-bond mutation before stack-top path flips, double and non-double path branches, target oxygen charge mutation, complete modeled atom/bond/stereo fields, explicit hydrogens, diagnostics, properties, graph preservation, and exception state | `RDKit✔️✔️` | none within the modeled valid `InchiToMol` `RWMol` boundary; source frame mechanically matches all lines 679-746, while the pinned oracle verifies complete observable state for every active branch |
| `_Valence5NCleanUp8` | `third_party/rdkit/External/INCHI-API/inchi.cpp:750` | `cargo test -p cosmolkit-inchi source_port__inchi__valence5ncleanup8__line_750` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__valence5ncleanup8__exact` (`1` executed; `18` complete pinned C++ records extracted from the provenance-checked official function) | exact pass for atomic-number, charge, and explicit-valence short-circuiting, entry explicit-valence exception, temporary nitrogen-to-tin mutation and restoration, every query atom and bond predicate, zero/one/multiple default atom-set-unique whole-molecule matches, unrelated unique tin match, exact five-bond mutation order, both charge updates, complete modeled atom/bond/stereo fields, explicit hydrogens, diagnostics, properties, graph preservation, and exception state | `RDKit✔️✔️` | none within the modeled valid `InchiToMol` `RWMol` boundary; source frame mechanically matches all lines 750-800, while the pinned oracle verifies complete observable state for every active branch |
| `_Valence5NCleanUp9` | `third_party/rdkit/External/INCHI-API/inchi.cpp:804` | `cargo test -p cosmolkit-inchi source_port__inchi__valence5ncleanup9__line_804` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__valence5ncleanup9__exact` (`1` executed; `18` complete pinned C++ records extracted from the provenance-checked official function) | exact pass for atomic-number, charge, and explicit-valence short-circuiting, entry explicit-valence exception, temporary nitrogen-to-tin mutation and restoration, every query atom and bond predicate, zero/one/multiple default atom-set-unique whole-molecule matches, unrelated unique tin match, exact three-bond mutation order, both charge updates, reversed endpoints, complete modeled atom/bond/stereo fields, explicit hydrogens, diagnostics, properties, graph preservation, and exception state | `RDKit✔️✔️` | none within the modeled valid `InchiToMol` `RWMol` boundary; source frame mechanically matches all lines 804-852, while the pinned oracle verifies complete observable state for every active branch |
| `_Valence5NCleanUpA` | `third_party/rdkit/External/INCHI-API/inchi.cpp:855` | `cargo test -p cosmolkit-inchi source_port__inchi__valence5ncleanupa__line_855` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__valence5ncleanupa__exact` (`1` executed; `12` complete pinned C++ records extracted from the provenance-checked official function) | exact pass for atomic-number and charge short-circuiting, explicit-valence mismatch and exception order, zero N=N matches, current-atom match exclusion, disconnected match, temporary N-to-Sn mutation and restoration, direct and alternating paths, path lengths 9 and 11, shortest-candidate replacement, equal-length first-candidate retention, stack-order bond flips, final charge and explicit-valence update, complete modeled atom/bond/stereo fields, explicit hydrogens, diagnostics, properties, graph preservation, and exception state | `RDKit✔️✔️` | none within the modeled valid `InchiToMol` `RWMol` boundary; source frame mechanically matches all lines 855-913, while the pinned oracle verifies complete observable state for every active branch |
| `_Valence5NCleanUpB` | `third_party/rdkit/External/INCHI-API/inchi.cpp:916` | `cargo test -p cosmolkit-inchi source_port__inchi__valence5ncleanupb__line_916` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__valence5ncleanupb__exact` (`1` executed; `7` complete pinned C++ records extracted from the provenance-checked official function) | exact pass for null target, exact target charge, one-bond depth limit, arbitrary source atom fields, adjacency-first equal-path selection, direction and explicit-hydrogen preservation, target charge mutation before target valence exception, bond and source charge mutation before source valence exception, complete modeled atom/bond/stereo fields, diagnostics, properties, graph preservation, and exception state | `RDKit✔️✔️` | none within the modeled valid `InchiToMol` `RWMol` boundary; source frame mechanically matches all lines 916-938, while the pinned oracle verifies complete observable state for every active branch |
| `_Valence7SCleanUp1` | `third_party/rdkit/External/INCHI-API/inchi.cpp:941` | `cargo test -p cosmolkit-inchi source_port__inchi__valence7scleanup1__line_941` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__valence7scleanup1__exact` (`1` executed; `12` complete pinned C++ records extracted from the provenance-checked official function) | exact pass for atomic-number and charge short-circuiting, explicit-valence mismatch and entry exception, oxygen/carbon/other-neighbor bond predicates, no-oxygen and selection-criteria failures, one-carbon and three-oxygen success criteria, last-valid-oxygen adjacency selection, exact bond/charge mutation order, selected-oxygen explicit-valence exception after partial mutation, complete modeled atom/bond/stereo fields, explicit hydrogens, diagnostics, properties, graph preservation, and exception state | `RDKit✔️✔️` | none within the modeled valid `InchiToMol` `RWMol` boundary; source frame mechanically matches all lines 941-986, while the pinned oracle verifies complete observable state for every active branch |
| `_Valence7SCleanUp2` | `third_party/rdkit/External/INCHI-API/inchi.cpp:989` | `cargo test -p cosmolkit-inchi source_port__inchi__valence7scleanup2__line_989` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__valence7scleanup2__exact` (`1` executed; `10` complete pinned C++ records extracted from the provenance-checked official function) | exact pass for atomic-number and charge short-circuiting, explicit-valence mismatch and entry exception, absent neutral-nitrogen target, source depth limit, direct and one-/two-intermediate paths, single/double/triple stack mutations, first equal-length adjacency path retention, root charge update, explicit-hydrogen and bond-direction preservation, complete modeled atom/bond/stereo fields, diagnostics, properties, graph preservation, and exception state | `RDKit✔️✔️` | none within the modeled valid `InchiToMol` `RWMol` boundary; source frame mechanically matches all lines 989-1018, while the pinned oracle verifies complete observable state for every active branch |
| `_Valence7SCleanUp3` | `third_party/rdkit/External/INCHI-API/inchi.cpp:1021` | `cargo test -p cosmolkit-inchi source_port__inchi__valence7scleanup3__line_1021` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__valence7scleanup3__exact` (`1` executed; `10` complete pinned C++ records extracted from the provenance-checked official function) | exact pass for atomic-number and charge short-circuiting, explicit-valence mismatch and entry exception, absent and charged nitrogen targets, one-bond source depth limit, direct success, first equal-length adjacency path retention, exact bond/target/source mutation order, source explicit-valence recomputation, deliberate target-valence non-recomputation, complete modeled atom/bond/stereo fields, explicit hydrogens, diagnostics, properties, graph preservation, and exception state | `RDKit✔️✔️` | none within the modeled valid `InchiToMol` `RWMol` boundary; source frame mechanically matches all lines 1021-1041, while the pinned oracle verifies complete observable state for every active branch |
| `_Valence8SCleanUp1` | `third_party/rdkit/External/INCHI-API/inchi.cpp:1044` | `cargo test -p cosmolkit-inchi source_port__inchi__valence8scleanup1__line_1044` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__valence8scleanup1__exact` (`1` executed; `12` complete pinned C++ records extracted from the provenance-checked official function) | exact pass for atomic-number and charge short-circuiting, explicit-valence mismatch and entry exception, absent and charged nitrogen targets, direct and alternating paths, accepted nine-bond and rejected eleven-bond depth boundaries, double/other stack branches, first equal-length adjacency path retention, exact path/target/source mutation order, target explicit-hydrogen reset, target-valence exception after partial mutation, complete modeled atom/bond/stereo fields, diagnostics, properties, graph preservation, and exception state | `RDKit✔️✔️` | none within the modeled valid `InchiToMol` `RWMol` boundary; source frame mechanically matches all lines 1044-1074, while the pinned oracle verifies complete observable state for every active branch |
| `_Valence8ClCleanUp1` | `third_party/rdkit/External/INCHI-API/inchi.cpp:1079` | `cargo test -p cosmolkit-inchi source_port__inchi__valence8clcleanup1__line_1079` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__valence8clcleanup1__exact` (`1` executed; `7` complete pinned C++ records extracted from the provenance-checked official function) | exact pass for explicit-valence and charge short-circuiting, valence-before-charge exception order, all-oxygen rejection, empty-neighbor vacuous success without an atomic-number guard, double- versus non-double-bond handling, exact root/bond/oxygen mutation order, ordered partial mutation before a later oxygen valence exception, complete modeled atom/bond/stereo fields, explicit hydrogens, diagnostics, properties, graph preservation, and exception state | `RDKit✔️✔️` | none within the modeled valid `InchiToMol` `RWMol` boundary; source frame mechanically matches all lines 1079-1111, while the pinned oracle verifies complete observable state for every active branch |
| `_Valence5ClCleanUp1` | `third_party/rdkit/External/INCHI-API/inchi.cpp:1114` | `cargo test -p cosmolkit-inchi source_port__inchi__valence5clcleanup1__line_1114` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__valence5clcleanup1__exact` (`1` executed; `8` complete pinned C++ records extracted from the provenance-checked official function) | exact pass for explicit-valence and charge short-circuiting, valence-before-charge exception order, oxygen element/charge and single-bond target predicates, absence of an atomic-number guard, first equal-length adjacency target retention, exact bond/root/target charge mutation order, deliberate target-valence non-recomputation, root-valence recomputation, complete modeled atom/bond/stereo fields, explicit hydrogens, diagnostics, properties, graph preservation, and exception state | `RDKit✔️✔️` | none within the modeled valid `InchiToMol` `RWMol` boundary; source frame mechanically matches all lines 1114-1131, while the pinned oracle verifies complete observable state for every active branch |
| `_Valence3ClCleanUp1` | `third_party/rdkit/External/INCHI-API/inchi.cpp:1134` | `cargo test -p cosmolkit-inchi source_port__inchi__valence3clcleanup1__line_1134` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__valence3clcleanup1__exact` (`1` executed; `7` complete pinned C++ records extracted from the provenance-checked official function) | exact pass for valence-before-charge exception order, explicit-valence and charge short-circuiting, sulfur element/charge and triple-bond target predicates, absence of an atomic-number guard, exact selected-bond mutation before root-valence recomputation, deliberate target-valence non-recomputation, complete modeled atom/bond/stereo fields, explicit hydrogens, diagnostics, properties, graph preservation, and exception state | `RDKit✔️✔️` | none within the modeled valid `InchiToMol` `RWMol` boundary; source frame mechanically matches all lines 1134-1149, while two equal triple-bond targets are unreachable after the required entry valence 3 check |
| `cleanUp` | `third_party/rdkit/External/INCHI-API/inchi.cpp:1151` | `cargo test -p cosmolkit-inchi source_port__inchi__cleanup__line_1151` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__cleanup__exact` (`1` executed; `37` complete pinned C++ records extracted from the provenance-checked official function) | exact pass for empty and ignored atoms, stored atom iteration order, nitrogen valence-4 routing, charged-nitrogen short-circuiting, aromatic clearing/restoration and exception preservation, every ordered valence-5 short-circuit operand, chlorine valence 8/5/3 routing, sulfur valence 7/8 routing, long alternating sulfur path, bromine degree and selenium predicates, complete atom/bond/aromatic/direction fields, diagnostics, properties, graph preservation, and partial mutation before exceptions | `RDKit✔️✔️` | none within the modeled valid `InchiToMol` `RWMol` boundary; source frame mechanically matches all lines 1151-1251, including the source-defined valence-8 sulfur call that is rejected by `_Valence8SCleanUp1`'s internal valence-7 guard |
| `InchiToMol` | `third_party/rdkit/External/INCHI-API/inchi.cpp:1254` | `cargo test -p cosmolkit-inchi source_port__inchi__inchitomol__line_1254` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__inchitomol__exact` (`1` executed; `33` complete pinned C++ records extracted from the provenance-checked official function) | exact pass for return and exception state, all `ExtraInchiReturnValues` fields, complete atom/bond/stereo fields, CIP ranks, warning/error bytes, toolkit call order, engine get/free counts, source allocation state, input/option bytes, success/warning/failure returns, empty structures, isotope-H, charge, radical, every bond/stereo mapping, duplicate and illegal bonds, double/tetrahedral/allene/unknown stereo, parity filtering and E/Z switching, unsigned-rank conversion, direction conflicts, sanitize/remove-H dispatch, and injected engine/toolkit exceptions | `RDKit✔️❌` | none within the modeled `Chem.MolFromInchi` adapter boundary; source-output cloning and checked graph access retain the documented performance gap |
| `fixOptionSymbol` | `third_party/rdkit/External/INCHI-API/inchi.cpp:1674` | `cargo test -p cosmolkit-inchi source_port__inchi__fixoptionsymbol__line_1674` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__fixoptionsymbol__exact` (`1` executed; `5` complete pinned C++ records extracted from the provenance-checked official function, including all non-NUL byte values) | exact pass for the GCC/Linux slash replacement and ordinary-byte branches, empty input, first-NUL termination, input preservation, output NUL position, exact-capacity output, and every untouched output-tail byte | `RDKit✔️✔️` | none within the production distinct-buffer, NUL-terminated caller boundary; the `_WIN32` option-prefix branch is preprocessor-inactive for the selected GCC/Linux target |
| `rCleanUp` | `third_party/rdkit/External/INCHI-API/inchi.cpp:1694` | `cargo test -p cosmolkit-inchi source_port__inchi__rcleanup__line_1694` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__rcleanup__exact` (`1` executed; `15` complete pinned C++ records extracted from the provenance-checked official function and charge-aware fixed SMILES query dependencies) | exact pass for no-match and topology-miss paths, non-default query charge matching, neutral-query charge wildcard behavior, all-negative and arbitrary wildcard oxygen charges, every neutral target position, precollected overlapping and disconnected matches, mutation and early-return ordering, reversed bond endpoints, and every modeled atom/bond/stereo field before and after mutation | `RDKit✔️❌` | none within the modeled `MolToInchi` adapter graph boundary; the dedicated safe fixed-query matcher preserves observable source behavior but retains extra temporary allocation and repeated sorting compared with RDKit's substructure matcher |
| `MolToInchi` | `third_party/rdkit/External/INCHI-API/inchi.cpp:1747` | `cargo test -p cosmolkit-inchi source_port__inchi__moltoinchi__line_1747` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__moltoinchi__exact` (`1` executed; `65` complete pinned C++ records extracted from the provenance-checked official function) | exact pass for return and exception state, all `ExtraInchiReturnValues` fields, property-cache and toolkit call order, clone-only `rCleanUp`, no/first conformer behavior, complete modeled `inchi_Atom` and `inchi_Stereo0D` fields, element/isotope/charge/radical/hydrogen conversion, tetrahedral degree/parity branches, all 22 bond types, all 7 bond directions with both endpoint orders, double-bond stereo/parity/swap/short-list branches, STEREOANY coordinate collapse, option conversion and embedded NUL termination, MAXVAL early return, nullable output preservation, GetINCHI/FreeINCHI counts and cleanup, original graph preservation, and injected toolkit/GetINCHI/FreeINCHI exceptions | `RDKit✔️❌` | none within the modeled toolkit-neutral `Chem.MolToInchi` adapter boundary; full molecule cloning and additional owned atom/stereo/option vectors retain the documented allocation and performance gap |
| `InchiToInchiKey` | `third_party/rdkit/External/INCHI-API/inchi.cpp:2145` | `cargo test -p cosmolkit-inchi source_port__inchi__inchitoinchikey__line_2145` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__inchitoinchikey__exact` (`1` executed; `9` complete pinned C++ records extracted from the provenance-checked official function) | exact pass for NUL-terminated and first-NUL-truncated successful keys, all six named return-code mappings, unknown return status, empty input, embedded input NUL, exact error-stream text and newline, direct-callee input prefix, fixed zero xtra arguments, call count, in-tree source-engine methane key, cleanup state, and structured engine/invalid-source-output errors | `RDKit✔️❌` | none within the modeled toolkit-neutral `InchiToInchiKey` adapter boundary; additional owned input/key/xtra buffers and diagnostic construction retain the documented allocation and performance gap |
| `MolToInchiKey` | `third_party/rdkit/External/INCHI-API/inchi.h:107` | `cargo test -p cosmolkit-inchi source_port__inchi__moltoinchikey__line_107` (`1` executed, exact pass) | `cargo test -p cosmolkit-inchi rdkit_cpp_oracle__moltoinchikey__exact` (`1` executed; `9` complete pinned C++ records extracted from the provenance-checked header-defined inline function and both pinned dependency frames) | exact pass for complete atom/bond/stereo/coordinate/options forwarding, original graph preservation, generation warning before key error diagnostics, successful first-NUL key construction, empty and fatal-code generation outputs, MAXVAL empty-return continuation, exact cross-callee call order, GetINCHI/FreeINCHI/key counts and outstanding output state, and toolkit/GetINCHI/FreeINCHI/key exception short-circuit and cleanup order | `RDKit✔️❌` | none within the modeled toolkit-neutral `Chem.MolToInchiKey` composition boundary; separate owned generation/key results and diagnostic-vector merging retain the documented allocation and performance gap |

## Remaining RDKit Adapter Functions

None. All 30 adapter functions in the source-level helper and public-entry
closure of the four selected APIs have complete source frames, focused tests
that execute, pinned exact C++ oracles, and closed first-axis markers.

`MolBlockToInchi` and `getInchiVersion` remain intentionally absent.

## Step 5309 Selected Closure Reaudit

- Compiler-derived scope: the selected five-root callback-complete official closure remains `1050` active definitions; the RDKit adapter closure remains `30` functions. No raw function-name inventory was used to expand or shrink this scope.
- Official completed rows: the `Closed During Core Execution` table contains `23` rows. Every row contains its specific focused Rust command, official C oracle command, exact result, final first-axis marker, preserved second-axis marker, and unclosed-branch field; automated row validation found `0` missing focused commands and `0` missing oracle commands.
- RDKit completed rows: the `Closed During RDKit Adapter Execution` table contains `30` rows. Every row contains its specific focused Rust command, pinned RDKit C++ oracle command, exact result, final first-axis marker, preserved second-axis marker, and unclosed-branch field; automated row validation found `0` missing focused commands and `0` missing oracle commands.
- Existing checked prefix: the previously completed official source locations remain covered by `dev/gap_reports/rdkit_inchi_completed_prefix_parity_audit.md` and the compiler-derived selected-closure accounting above; no completed marker or historical step was reopened by this audit.
- Function-level open-marker scan: `rg -n '^\s*// (INCHI|RDKit)❗' crates/cosmolkit-inchi/src --glob '*.rs'` returns exactly one row, `crates/cosmolkit-inchi/src/source/base/ichirvr4.rs:8364`, for `NormalizeAndCompare` with marker `INCHI❗❌`.
- Exact result: all `1049` source-defined, behavior-closed selected official locations and all `30` selected RDKit adapter functions retain their recorded exact focused/oracle results. The second-axis markers are preserved independently and are not promoted by this audit.
- Remaining branch: `NormalizeAndCompare` initial string-buffer allocation failure reaches official `strcpy(existing_formula, NULL)`, so official C defines no exact return, mutation, signal, or cleanup sequence. Its `240` source-defined records match exactly; the separate undefined-path record is termination evidence only.
- Decision: the selected call graph still has one open first axis. Stop before Step 5310 and do not expose `mol_to_inchi`, `mol_to_inchi_key`, `inchi_to_inchi_key`, or `mol_from_inchi` under the current completion rule.

## Step 5311 Human-Authorized Undefined-Behavior Policy

On 2026-07-31, the human author explicitly authorized Rust to return a
structured error for official-C undefined behavior. This authorization is
narrowly applied to the audited `NormalizeAndCompare` initial string-buffer
allocation-failure path and does not authorize inferred behavior for any other
undefined, unported, or unverified source path.

- Source-defined behavior remains unchanged: all `240` source-defined `NormalizeAndCompare` records must continue to match official C field-for-field.
- Official undefined-path evidence remains unchanged: the official oracle may establish that execution reaches `strcpy(existing_formula, NULL)`, but no observed signal, exit code, mutation, or cleanup sequence is an exact expected result.
- Rust behavior for this path is deterministic `SourceHeapError::AllocationFailed`. It must be verified independently for caller-visible mutation order, allocation state, retained ownership, cleanup state, and absence of fallback.
- The source marker remains `INCHI❗❌`. The first axis is not upgraded because the complete active C function includes a path with no source-defined behavior; the second-axis performance marker is independent and unchanged.
- Exact-parity claims remain limited to source-defined behavior. The structured Rust error is an explicit project behavior, not an official-C parity result.
- Production remains Rust-only. This decision does not permit production FFI, external commands, process-signal emulation, SMILES/MolBlock transit, heuristics, placeholders, or silent fallback.
- The Step 5309 exposure stop is superseded only after the dedicated structured-error focused test and the unchanged official-C oracle both execute and pass. Until then, public API exposure remains blocked.

## Step 5319 Authorized Undefined-Path Closure

The two conditions in Step 5311 now execute and pass:

- Structured-error focused command: `cargo test -p cosmolkit-inchi source_policy__normalizeandcompare__undefined_initial_buffer_allocation_returns_structured_error -- --nocapture`. The main harness executed exactly `1` matching test and reported `1 passed; 0 failed; 1265 filtered out`. The later provenance harness executed `0` tests and is not counted. The test verifies the exact `SourceHeapError::AllocationFailed`, `pnNumRunBNS` and `pnTotalDelta` preservation, allocation-call and outstanding-allocation state, generated formula and reversed-InChI cleanup, original caller-input preservation, complete caller state, exact cleanup-event order, and absence of comparison, repair, or fallback calls.
- Repeated official-C command: `cargo test -p cosmolkit-inchi official_c_oracle__normalizeandcompare__exact -- --nocapture`. The main harness executed exactly `1` matching test and reported `1 passed; 0 failed; 1265 filtered out`. The oracle processed exactly `241` records. All `240` source-defined records matched Rust field-for-field with zero mismatches. The single `zz-zy-undefined` record remains separate evidence that official C terminates after reaching the null-source `strcpy`; it has no asserted exact result, signal, mutation, or cleanup sequence.
- Final marker: `NormalizeAndCompare` remains `INCHI❗❌`. The first axis is intentionally open because official C has no behavior to reproduce on the authorized path, and the independently reviewed second-axis performance gap remains unchanged.
- Revised exposure gate: the Step 5309 stop is superseded for the four audited toolkit-neutral scalar APIs only. Integration may proceed because every source-defined selected behavior is closed exactly and the one audited official-C undefined path has the human-authorized deterministic structured Rust error. This is not an exact-parity result for that path and does not authorize any other undefined, unverified, or unsupported path.
- Remaining selected active call graph: no source-implementation Port unit remains. One selected source location retains an open first-axis marker solely for the authorized undefined path; it does not permit a heuristic, placeholder, partial result, silent fallback, production FFI, external production command, or SMILES/MolBlock transit.

## Schematic Issue 11 Pointer Regression Reaudit

Eight valid rows in `tests/smiles_5000.smi` previously failed Rust InChI
generation with `SourceHeapError::PointerOutOfBounds`: `440`, `904`, `2040`,
`2556`, `3248`, `3620`, `3811`, and `4944`. The failing active call path was:

```text
RunBalancedNetworkSearch
  -> BalancedNetworkSearch
  -> bIgnoreVertexNonTACN_group
  -> GetPrevVertex
```

The source defect was in the existing `BalancedNetworkSearch` port, not a
missing function. Official `ichi_bns.c:10951-10957` evaluates
`bIgnoreVertexNonTACN_group(pBNS, prim(v), u, SwitchEdge)` only after
`TREE_IS_S_REACHABLE(prim(v))`, the T-prime edge exclusion, and `b_u != b_v`
have all succeeded. Rust had evaluated that helper before the enclosing
condition, so it followed an official-C-unreachable `SwitchEdge[-6]` path.
The Rust expression at `crates/cosmolkit-inchi/src/source/base/ichi_bns.rs:13540`
now keeps the fallible helper call in the same lazy `&&` position as the
official source. No pointer bound was widened, no error was swallowed, and no
input-specific branch, heuristic, fallback, FFI, or external production call
was added.

The focused command was:

```bash
cargo test -q -p cosmolkit-core \
  --features op-contracts-strict \
  --lib \
  inchi::tests::inchi_generation_matches_chematic_issue_11_pointer_regressions \
  -- --exact --nocapture
```

The main harness executed exactly one matching test and reported `1 passed; 0
failed; 2573 filtered out`. The test covers all eight original failures and
asserts the complete pinned-RDKit InChI string, official warning return code
`1`, and complete InChIKey for every row.

After rebuilding the Python extension from the repaired Rust source, a
fail-closed public-entry comparison parsed every corpus row independently with
pinned RDKit `2026.03.1` and COSMolKit, then compared
`rdkit.Chem.MolToInchi` with `cosmolkit.Chem.MolToInchi` by string equality.
The result was exactly `5000` equal InChI strings, `0` mismatches, `0` RDKit
parse errors, `0` COSMolKit parse errors, `0` RDKit generation errors, and `0`
COSMolKit generation errors. Any nonzero failure category made the harness
exit nonzero. This regression closes no new function and changes no existing
source or performance marker.

## Frozen And Excluded Active Definitions

Every configured active definition outside the five-root callback-complete closure is frozen or unscheduled below. A `completed-frozen` row keeps its existing Rust code but receives no further Port or parity step. An `excluded-unported` row is not part of the new plan.

| Function | Source | Disposition |
|---|---|---|
| `Free_std_inchi_Input` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ichilnct.c:182` | completed-frozen |
| `Free_inchi_Input` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ichilnct.c:187` | completed-frozen |
| `InchiToInchi_Input` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ichilnct.c:204` | completed-frozen |
| `Get_std_inchi_Input_FromAuxInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ichilnct.c:79` | completed-frozen |
| `Get_inchi_Input_FromAuxInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ichilnct.c:89` | completed-frozen |
| `FreeStdINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:183` | completed-frozen |
| `FreeStructFromStdINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:195` | completed-frozen |
| `GetStringLength` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:2086` | completed-frozen |
| `GetStdINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:242` | completed-frozen |
| `GetStructFromStdINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:2461` | completed-frozen |
| `FreeStructFromINCHIEx` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:2892` | completed-frozen |
| `GetINCHIEx` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:325` | completed-frozen |
| `STDINCHIGEN_Reset` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1125` | completed-frozen |
| `INCHIGEN_Reset` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1134` | completed-frozen |
| `STDINCHIGEN_Create` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:122` | completed-frozen |
| `STDINCHIGEN_Destroy` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1308` | completed-frozen |
| `INCHIGEN_Create` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:130` | completed-frozen |
| `INCHIGEN_Destroy` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1315` | completed-frozen |
| `STDINCHIGEN_Setup` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:185` | completed-frozen |
| `INCHIGEN_Setup` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:260` | completed-frozen |
| `STDINCHIGEN_DoNormalization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:435` | completed-frozen |
| `INCHIGEN_DoNormalization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:442` | completed-frozen |
| `STDINCHIGEN_DoCanonicalization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:689` | completed-frozen |
| `INCHIGEN_DoCanonicalization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:697` | completed-frozen |
| `STDINCHIGEN_DoSerialization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:859` | completed-frozen |
| `INCHIGEN_DoSerialization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:869` | completed-frozen |
| `Normalization_step` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:1138` | completed-frozen |
| `Canonicalization_step` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:1612` | completed-frozen |
| `NormOneStructureINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:188` | completed-frozen |
| `FillOutINChIReducedWarn` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:2254` | completed-frozen |
| `make_norm_atoms_from_inp_atoms` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:2937` | completed-frozen |
| `CanonOneStructureINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:503` | completed-frozen |
| `NormOneComponentINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:654` | completed-frozen |
| `CanonOneComponentINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:923` | completed-frozen |
| `FreeInchi_Atom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:106` | completed-frozen |
| `CreateInchiAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:117` | completed-frozen |
| `MakeINCHIFromMolfileText` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:127` | excluded-unported |
| `PrepareToMakeINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:285` | completed-frozen |
| `PostMakeINCHICleanup` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:391` | completed-frozen |
| `FreeInchi_Input` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:439` | completed-frozen |
| `is_in_the_slist` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:531` | completed-frozen |
| `is_element_a_metal` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:543` | completed-frozen |
| `InchiToInchiAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:582` | completed-frozen |
| `GetSingleStereoCode` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:145` | excluded-unported |
| `IXA_INCHIBUILDER_SetOption_Stereo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1583` | excluded-unported |
| `IXA_INCHIBUILDER_SetOption_Timeout` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1597` | excluded-unported |
| `IXA_INCHIBUILDER_SetOption_Timeout_MilliSeconds` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1612` | excluded-unported |
| `IXA_INCHIBUILDER_SetOption` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1626` | excluded-unported |
| `IXA_INCHIBUILDER_CheckOption` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1776` | excluded-unported |
| `GetDoubleStereoCode` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:181` | excluded-unported |
| `IXA_INCHIBUILDER_CheckOption_Stereo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1971` | excluded-unported |
| `BUILDER_Unpack` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:198` | excluded-unported |
| `IXA_INCHIBUILDER_GetOption_Timeout_MilliSeconds` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1993` | excluded-unported |
| `IXA_INCHIBUILDER_GetInChIVersion` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:2008` | excluded-unported |
| `IXA_INCHIBUILDER_GetInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:2022` | excluded-unported |
| `IXA_INCHIBUILDER_GetInChIEx` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:2037` | excluded-unported |
| `IXA_INCHIBUILDER_GetAuxInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:2052` | excluded-unported |
| `IXA_INCHIBUILDER_GetLog` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:2067` | excluded-unported |
| `BUILDER_Pack` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:212` | excluded-unported |
| `TranslateTetrahedralVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:219` | excluded-unported |
| `ExtendAllene` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:246` | excluded-unported |
| `ExtendCumulene` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:313` | excluded-unported |
| `IsRectangularVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:361` | excluded-unported |
| `IsRectOrAntiRectCentre` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:392` | excluded-unported |
| `ClearMolecule` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:422` | excluded-unported |
| `AppendOption` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:441` | excluded-unported |
| `BUILDER_ClearOptions` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:453` | excluded-unported |
| `BUILDER_Update` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:513` | excluded-unported |
| `IXA_INCHIBUILDER_Create` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:802` | excluded-unported |
| `IXA_INCHIBUILDER_Destroy` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:823` | excluded-unported |
| `IXA_INCHIBUILDER_SetMolecule` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:841` | excluded-unported |
| `IXA_INCHIKEYBUILDER_SetInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_inchikey_builder.c:107` | excluded-unported |
| `IXA_INCHIKEYBUILDER_GetInChIKey` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_inchikey_builder.c:131` | excluded-unported |
| `KEYBUILDER_Unpack` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_inchikey_builder.c:59` | excluded-unported |
| `KEYBUILDER_Pack` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_inchikey_builder.c:73` | excluded-unported |
| `IXA_INCHIKEYBUILDER_Create` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_inchikey_builder.c:80` | excluded-unported |
| `IXA_INCHIKEYBUILDER_Destroy` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_inchikey_builder.c:95` | excluded-unported |
| `IXA_MOL_GetAtomY` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1000` | excluded-unported |
| `IXA_MOL_SetAtomZ` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1023` | excluded-unported |
| `IXA_MOL_GetAtomZ` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1047` | excluded-unported |
| `IXA_MOL_SetAtomElement` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1070` | excluded-unported |
| `IXA_MOL_GetAtomElement` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1110` | excluded-unported |
| `IXA_MOL_SetAtomAtomicNumber` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1133` | excluded-unported |
| `IXA_MOL_GetAtomAtomicNumber` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1163` | excluded-unported |
| `IXA_MOL_SetAtomHydrogens` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1186` | excluded-unported |
| `GetVertexCount` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:118` | excluded-unported |
| `IXA_MOL_GetAtomHydrogens` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1236` | excluded-unported |
| `IXA_MOL_SetAtomMass` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1271` | excluded-unported |
| `IXA_MOL_GetAtomMass` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1315` | excluded-unported |
| `IXA_MOL_SetAtomRadical` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1338` | excluded-unported |
| `MOL_PackAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:135` | excluded-unported |
| `IXA_MOL_GetAtomRadical` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1362` | excluded-unported |
| `IXA_MOL_SetAtomCharge` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1385` | excluded-unported |
| `IXA_MOL_GetAtomCharge` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1409` | excluded-unported |
| `MOL_PackBond` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:142` | excluded-unported |
| `IXA_MOL_ReserveSpace` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1432` | excluded-unported |
| `IXA_MOL_CreateBond` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1497` | excluded-unported |
| `MOL_PackStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:149` | excluded-unported |
| `MOL_PackPolymerUnit` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:156` | excluded-unported |
| `MOL_UnpackAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:163` | excluded-unported |
| `IXA_MOL_GetNumBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1643` | excluded-unported |
| `IXA_MOL_GetBondId` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1657` | excluded-unported |
| `IXA_MOL_GetBondIndex` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1684` | excluded-unported |
| `IXA_MOL_GetBondAtom1` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1707` | excluded-unported |
| `IXA_MOL_GetBondAtom2` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1732` | excluded-unported |
| `IXA_MOL_GetBondOtherAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1756` | excluded-unported |
| `IXA_MOL_SetBondType` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1796` | excluded-unported |
| `MOL_UnpackBond` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:180` | excluded-unported |
| `IXA_MOL_GetBondType` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1822` | excluded-unported |
| `IXA_MOL_SetBondWedge` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1844` | excluded-unported |
| `IXA_MOL_GetBondWedge` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1881` | excluded-unported |
| `IXA_MOL_SetDblBondConfig` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1917` | excluded-unported |
| `IXA_MOL_GetDblBondConfig` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1941` | excluded-unported |
| `IXA_MOL_GetCommonBond` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1964` | excluded-unported |
| `MOL_UnpackStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:197` | excluded-unported |
| `IXA_MOL_CreateStereoTetrahedron` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2010` | excluded-unported |
| `IXA_MOL_CreateStereoRectangle` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2045` | excluded-unported |
| `IXA_MOL_CreateStereoAntiRectangle` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2079` | excluded-unported |
| `IXA_MOL_GetNumStereos` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2113` | excluded-unported |
| `IXA_MOL_GetStereoId` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2127` | excluded-unported |
| `IXA_MOL_GetStereoIndex` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2148` | excluded-unported |
| `MOL_UnpackPolymerUnit` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:214` | excluded-unported |
| `IXA_MOL_GetStereoTopology` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2164` | excluded-unported |
| `IXA_MOL_GetStereoCentralAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2186` | excluded-unported |
| `IXA_MOL_GetStereoCentralBond` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2219` | excluded-unported |
| `IXA_MOL_GetStereoNumVertices` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2252` | excluded-unported |
| `IXA_MOL_GetStereoVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2274` | excluded-unported |
| `IXA_MOL_SetStereoParity` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2310` | excluded-unported |
| `IXA_MOL_GetStereoParity` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2331` | excluded-unported |
| `IXA_MOL_SetPolymerUnit` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2353` | excluded-unported |
| `MOL_GetAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:238` | excluded-unported |
| `IXA_MOL_SetExtMolDataByInChIExtInput` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2429` | excluded-unported |
| `MOL_GetBond` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:249` | excluded-unported |
| `MOL_GetStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:260` | excluded-unported |
| `IXA_MOL_CreatePolymerUnit` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2681` | excluded-unported |
| `IXA_MOL_GetPolymerUnitId` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2696` | excluded-unported |
| `IXA_MOL_GetPolymerUnitIndex` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2717` | excluded-unported |
| `MOL_GetSGroup` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:271` | excluded-unported |
| `MOL_ClearExtMolData` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:285` | excluded-unported |
| `MOL_Clear` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:395` | excluded-unported |
| `MOL_GuessNewSize` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:408` | excluded-unported |
| `MOL_CreateAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:428` | excluded-unported |
| `MOL_CreateStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:520` | excluded-unported |
| `MOL_CreatePolymerUnit` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:598` | excluded-unported |
| `MOL_Pack` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:671` | excluded-unported |
| `MOL_Unpack` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:678` | excluded-unported |
| `MOL_GetBondOtherAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:692` | excluded-unported |
| `IXA_MOL_Create` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:729` | excluded-unported |
| `IXA_MOL_Destroy` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:744` | excluded-unported |
| `IXA_MOL_Clear` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:759` | excluded-unported |
| `IXA_MOL_SetChiral` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:773` | excluded-unported |
| `IXA_MOL_GetChiral` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:788` | excluded-unported |
| `IXA_MOL_CreateAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:801` | excluded-unported |
| `IXA_MOL_GetNumAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:816` | excluded-unported |
| `IXA_MOL_GetAtomId` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:830` | excluded-unported |
| `IXA_MOL_GetAtomIndex` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:857` | excluded-unported |
| `IXA_MOL_GetAtomNumBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:879` | excluded-unported |
| `IXA_MOL_GetAtomBond` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:902` | excluded-unported |
| `IXA_MOL_SetAtomX` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:929` | excluded-unported |
| `IXA_MOL_GetAtomX` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:953` | excluded-unported |
| `IXA_MOL_SetAtomY` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:976` | excluded-unported |
| `FindCumuleneCentre` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_read_inchi.c:112` | excluded-unported |
| `IXA_MOL_ReadInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_read_inchi.c:189` | excluded-unported |
| `AnalyseInternalVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_read_inchi.c:55` | excluded-unported |
| `IXA_MOL_ReadAuxInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_read_inchi.c:686` | excluded-unported |
| `IXA_MOL_ReadMolfile` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_read_mol.c:78` | excluded-unported |
| `STATUS_Clear` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:102` | excluded-unported |
| `STATUS_Pack` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:123` | excluded-unported |
| `STATUS_Unpack` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:131` | excluded-unported |
| `STATUS_PushMessage` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:138` | excluded-unported |
| `INCHISTATUS_TestSeverity` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:211` | excluded-unported |
| `IXA_STATUS_Create` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:233` | excluded-unported |
| `IXA_STATUS_Destroy` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:246` | excluded-unported |
| `IXA_STATUS_Clear` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:259` | excluded-unported |
| `IXA_STATUS_HasError` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:271` | excluded-unported |
| `IXA_STATUS_HasWarning` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:278` | excluded-unported |
| `IXA_STATUS_GetCount` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:285` | excluded-unported |
| `IXA_STATUS_GetMessage` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:298` | excluded-unported |
| `IXA_STATUS_GetSeverity` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:316` | excluded-unported |
| `BLOCK_clear` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:81` | excluded-unported |
| `STATUS_init` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:93` | excluded-unported |
| `dbl2int_f` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/bcf_s.c:117` | excluded-unported |
| `dbl2int_e` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/bcf_s.c:188` | excluded-unported |
| `dbl2int_g` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/bcf_s.c:276` | excluded-unported |
| `max_3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/bcf_s.c:58` | excluded-unported |
| `memcpy_custom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/bcf_s.c:68` | excluded-unported |
| `dbl2int` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/bcf_s.c:92` | excluded-unported |
| `fix_explicitly_indicated_bonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:2704` | excluded-unported |
| `update_some_attype_totals` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:2846` | excluded-unported |
| `inchi_ios_create_copy` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:106` | excluded-unported |
| `inchi_ios_flush` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:144` | excluded-unported |
| `inchi_strbuf_printf_from` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1544` | excluded-unported |
| `inchi_strbuf_getline` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1584` | excluded-unported |
| `inchi_strbuf_addline` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1635` | completed-frozen |
| `_inchi_trace` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1699` | excluded-unported |
| `Output_RecordInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1709` | excluded-unported |
| `inchi_ios_flush2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:183` | excluded-unported |
| `inchi_ios_free_str` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:272` | excluded-unported |
| `inchi_ios_str_getc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:294` | completed-frozen |
| `inchi_ios_str_gets` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:324` | completed-frozen |
| `inchi_ios_str_getsTab` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:354` | completed-frozen |
| `inchi_ios_gets` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:386` | completed-frozen |
| `inchi_ios_getsTab` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:420` | completed-frozen |
| `inchi_ios_getsTab1` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:451` | completed-frozen |
| `push_to_winchi_text_window` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:576` | excluded-unported |
| `inchi_ios_flush_not_displayed` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:679` | excluded-unported |
| `inchi_fprintf` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:800` | excluded-unported |
| `inchi_fgetsLfTab` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:893` | excluded-unported |
| `inchi_fgetsTab` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:966` | excluded-unported |
| `ProcessStructError` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2511` | excluded-unported |
| `is_Aryl` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1486` | excluded-unported |
| `is_Saturated_C` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1509` | excluded-unported |
| `is_C_Alk` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1517` | excluded-unported |
| `is_P_TB_N` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1740` | excluded-unported |
| `get_CO_opposite` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1763` | excluded-unported |
| `set_R2C_el_numbers` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:712` | excluded-unported |
| `has_atom_pair` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:835` | excluded-unported |
| `OutputINChIPlainError` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:514` | excluded-unported |
| `OutputInChIOutOfStrFromINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:5888` | excluded-unported |
| `CheckINCHIKey` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_dll.c:513` | excluded-unported |
| `fprint_digest` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_dll.c:601` | excluded-unported |
| `GetStdINCHIKeyFromStdINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_dll.c:92` | excluded-unported |
| `CreateOrigInpDataFromMolfile` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:108` | excluded-unported |
| `SetExtOrigAtDataByMolfileExtInput` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1340` | completed-frozen |
| `ReadMolfileToInpAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:309` | excluded-unported |
| `MakeInpAtomsFromMolfileData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:508` | completed-frozen |
| `calculate_valences` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:807` | completed-frozen |
| `SetInpAtomsXYZ` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:975` | completed-frozen |
| `MolfileReadSgroupOfPolymer` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:1425` | completed-frozen |
| `MolfileReadDataLines` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:165` | completed-frozen |
| `MolfileTreatPseudoElementAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:1860` | completed-frozen |
| `MolfileReadHeaderLines` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:451` | completed-frozen |
| `MolfileReadCountsLine` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:569` | completed-frozen |
| `MolfileReadAtomsBlock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:661` | completed-frozen |
| `MolfileReadBondsBlock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:820` | completed-frozen |
| `MolfileReadSTextBlock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:917` | completed-frozen |
| `ReadMolfile` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:92` | completed-frozen |
| `MolfileReadPropBlock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:960` | completed-frozen |
| `MolfileReadField` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:105` | completed-frozen |
| `MolfileExtractStrucNum` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:333` | completed-frozen |
| `MolfileHasNoChemStruc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:362` | completed-frozen |
| `MolfileGetXYZDimAndNormFactors` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:474` | completed-frozen |
| `MolfileStrnread` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:65` | completed-frozen |
| `FreeMolfileData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:664` | completed-frozen |
| `MolfileV3000ReadBondsBlock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:1262` | completed-frozen |
| `get_actual_atom_number` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:1585` | completed-frozen |
| `MolfileV3000ReadTailOfCTAB` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:1602` | completed-frozen |
| `DeleteMolfileV3000Info` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:165` | completed-frozen |
| `MolfileV3000ReadHapticBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:1732` | completed-frozen |
| `MolfileV3000ReadStereoCollection` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:1817` | completed-frozen |
| `get_V3000_input_line_to_strbuf` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:1889` | completed-frozen |
| `inchi_fgetsLf_V3000` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:216` | excluded-unported |
| `MolfileV3000ReadField` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:257` | completed-frozen |
| `MolfileV3000ReadKeyword` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:410` | completed-frozen |
| `MolfileV3000ReadCTABBeginAndCountsLine` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:444` | completed-frozen |
| `MolfileV3000ReadSGroup` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:561` | completed-frozen |
| `MolfileV3000Read3DBlock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:608` | completed-frozen |
| `MolfileV3000ReadCollections` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:647` | completed-frozen |
| `MolfileV3000Init` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:68` | completed-frozen |
| `MolfileV3000ReadAtomsBlock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:862` | completed-frozen |
| `SDFileSkipExtraData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:161` | completed-frozen |
| `SDFileIdentifyLabel` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:346` | completed-frozen |
| `SDFileExtractCASNo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:407` | completed-frozen |
| `NumLists_Alloc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:433` | completed-frozen |
| `NumLists_ReAlloc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:449` | completed-frozen |
| `NumLists_Append` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:471` | completed-frozen |
| `NumLists_Free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:491` | completed-frozen |
| `MolFmtSgroup_Create` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:634` | completed-frozen |
| `MolFmtSgroup_Free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:658` | completed-frozen |
| `MolFmtSgroups_Alloc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:675` | completed-frozen |
| `MolFmtSgroups_ReAlloc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:693` | completed-frozen |
| `MolFmtSgroups_Append` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:715` | completed-frozen |
| `MolFmtSgroups_Free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:751` | completed-frozen |
| `MolFmtSgroups_GetIndexBySgroupId` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:769` | completed-frozen |
| `CreateInchi_Stereo0D` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:118` | completed-frozen |
| `FreeInchi_Stereo0D` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:125` | completed-frozen |
| `find_and_interpret_structure_header` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:1742` | excluded-unported |
| `FindToken` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:445` | completed-frozen |
| `LoadLine` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:489` | completed-frozen |
| `InchiToInpAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:549` | excluded-unported |
| `InchiToOrigAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1067` | excluded-unported |
| `GetOneStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:134` | excluded-unported |
| `POSEContext_DebugPrint` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1672` | excluded-unported |
| `ReadTheStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:396` | excluded-unported |
| `OAD_PolymerUnit_CompareAtomLists` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1400` | excluded-unported |
| `OrigAtData_RemoveAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2564` | excluded-unported |
| `OAD_CollectFragmentBondsAndAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2805` | excluded-unported |
| `OAD_CollectBackboneAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2946` | excluded-unported |
| `DisplayStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1222` | excluded-unported |
| `winchi_calc_inchikey` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:421` | excluded-unported |
| `sha2_file` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:288` | excluded-unported |
| `sha2_hmac` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:325` | excluded-unported |
| `sha2_self_test` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:429` | excluded-unported |
| `stbsp_sprintf` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:1376` | excluded-unported |
| `stbsp__clamp_callback` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:1393` | excluded-unported |
| `stbsp__count_clamp_callback` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:1421` | excluded-unported |
| `stbsp_vsnprintf` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:1430` | excluded-unported |
| `stbsp_snprintf` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:1460` | excluded-unported |
| `stbsp_vsprintf` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:1472` | excluded-unported |
| `stbsp__real_to_parts` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:1491` | excluded-unported |
| `stbsp__raise_to_power10` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:1638` | excluded-unported |
| `stbsp__real_to_str` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:1705` | excluded-unported |
| `stbsp_set_separators` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:274` | excluded-unported |
| `stbsp__lead_sign` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:294` | excluded-unported |
| `stbsp__strlen_limited` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:309` | excluded-unported |
| `stbsp_vsprintfcb` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:349` | excluded-unported |
| `subgraf_debug_trace` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5192` | excluded-unported |
| `get_atomic_mass` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1040` | completed-frozen |
| `normalize_string` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1589` | completed-frozen |
| `dotify_non_printable_chars` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1630` | completed-frozen |
| `read_upto_delim` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1658` | completed-frozen |
| `is_matching_any_delim` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1710` | completed-frozen |
| `remove_trailing_spaces` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1728` | completed-frozen |
| `inchi__strnset` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:2022` | excluded-unported |

## Historical Non-Function Anchors

The following checked historical source locations do not correspond to active GCC function bodies under the selected target. They remain historical records and are not rescheduled:

- `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1552`
- `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1561`
- `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1570`

## Completion Rule

No selected public API may be exposed until all 23 official functions and all
30 RDKit adapter functions have complete source frames, focused tests that
actually execute, and exact non-production oracles for every source-defined
behavior. Closed first-axis markers are required except for the single audited
`NormalizeAndCompare` undefined initial string-buffer allocation path authorized
in Step 5311. That exception requires a deterministic structured Rust error, a
dedicated focused test that executes, preserved `INCHI❗❌`, and explicit
documentation that the path is not an exact-parity result. It does not apply to
future undefined or unverified paths without separate human-author approval.
Performance markers remain independent. Exact parity does not permit ignored
fields, tolerances, compatibility assertions, skipped malformed inputs, or
molecule-specific exceptions.
