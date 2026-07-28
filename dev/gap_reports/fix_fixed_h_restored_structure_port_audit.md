# `FixFixedHRestoredStructure` source-port audit

## Scope and active configuration

- Plan step: 3725.
- Source frame: `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr3.c:333-6534`.
- Active target: official `libinchi` CMake target with `COMPILE_ANSI_ONLY`,
  `TARGET_API_LIB`, GCC/Linux, and `READ_INCHI_STRING == 1`.
- `INCHI_HEAPCHK` expands to no production operation.
- `INC_ADD_EDGE` is the active integer macro value `64`.
- `MAX_DIFF_FIXH` is `256`.
- The function is active and export-root reachable. It is not CLI/demo-only code.

Step 3725 is integrated as one production function. The complete active source
frame and all ordered repair cases are present in the Rust function; no case is
implemented by a companion copy, fallback, or test-only production path.

## Validation result

The focused command is:

```bash
cargo test -p cosmolkit-inchi source_port__ichirvr3__fixfixedhrestoredstructure__line_333
```

The main `src/lib.rs` harness ran exactly one matching test and passed: `1
passed; 0 failed; 908 filtered out`. The separate `official_c_oracle` harness
ran zero matching tests and is not counted as evidence for this step.

The focused test directly invokes the production function and validates the
no-fixed-H exit, both initial comparison exits, charge-list allocation failure,
successful temporary-list allocation followed by the source `RI_ERR_PROGR`
exit, both canonical-map allocation failure ordinals, outer temporary-list
cleanup, retained caller-owned allocations, exact return values, and intentional
non-modification of the four unused source parameters.

## Atomic source boundary

The function is one source-defined operation containing ordered repair cases.
Each successful case may mutate BNS flow, reconstruct InChI, refill canonical
maps, recompute `CMP2FHINCHI`, and thereby change the eligibility of every later
case. The cases therefore cannot be implemented as independent public behavior,
reordered, or short-circuited from stale comparison data.

| Ordered block | Source start | Local exit label |
|---|---:|---|
| Initial edge collection and comparison | 333 | `exit_function` |
| case 01 | 516 | none |
| case 02 | 649 | none |
| case 03 | 799 | none |
| case 03a | 948 | none |
| case 04 | 1084 | none |
| case 05 | 1228 | none |
| case 06c | 1361 | `exit_case_06c` at 1488 |
| case 06d | 1526 | `exit_case_06d` at 1649 |
| case 06 | 1683 | none |
| case 06a | 1816 | `exit_case_06a` at 1913 |
| case 06b | 1950 | `exit_case_06b` at 2060 |
| case 06e | 2097 | `exit_case_06e` at 2294 |
| case 07 | 2329 | none |
| case 07a | 2470 | none |
| case 08 | 2613 | none |
| case 09 | 2763 | none |
| case 10 | 2920 | none |
| case 11 | 3093 | none |
| case 12 | 3304 | none |
| case 13 | 3519 | none |
| case 14 | 3701 | `exit_case_14` at 3989 |
| case 15 | 4033 | none |
| case 16 | 4210 | none |
| case 17 | 4370 | none |
| case 18 | 4545 | `exit_case_18` at 4649 |
| case 19 | 4688 | none |
| case 20 | 4842 | none |
| case 21 | 5016 | none |
| case 21a | 5224 | `exit_case_21a` at 5428 |
| case 22 | 5465 | none |
| case 23 | 5609 | `exit_case_23` at 5933 |
| case 24 | 5978 | `exit_case_24` at 6222 |
| case 25 | 6272 | `exit_case_25` at 6480 |
| Mandatory shared cleanup and return mapping | 6523 | `exit_function` |

The source labels above are intra-function control-flow targets, not callable C
functions and not candidates for companion Rust ports.

## Direct callee closure status

All compiler-verified direct callees already have source-backed Rust functions.
Their call multiplicities in the active C body are retained here to detect
accidental omission during translation.

| Callee | Calls | Rust source status before Step 3725 |
|---|---:|---|
| `AddToEdgeList` | 105 | present |
| `AllocEdgeList` | 44 | present |
| `FillOutCMP2FHINCHI` | 34 | present; Step 3721 |
| `FillOutExtraFixedHDataRestr` | 34 | present; Step 3717 |
| `FillTgDiffHChgFH` | 1 | present |
| `FindInEdgeList` | 18 | present |
| `GetChargeFlowerUpperEdge` | 7 | present |
| `MakeOneInChIOutOfStrFromINChI2` | 33 | present |
| `RemoveForbiddenEdgeMask` | 75 | present |
| `RemoveFromEdgeListByValue` | 1 | present |
| `RunBnsRestoreOnce` | 33 | present |
| `RunBnsTestOnce` | 59 | present |
| `SetForbiddenEdgeMask` | 41 | present |
| `bHas_N_V` | 1 | present |
| `is_bond_in_Nmax_memb_ring` | 1 | present |

`inchi_min` is header macro behavior and must remain an inline comparison, not
a fabricated source function. `memset` is libc behavior. `INCHI_HEAPCHK` is the
active empty macro and must not become observable Rust behavior.

## Mutable state and cleanup contract

- Seven outer edge lists are initialized with `EDGE_LIST_CLEAR` before any
  chemistry branch and are always released at `exit_function` in source order:
  `AllChargeEdges`, `CurrEdges`, `NFlowerEdges`, `SFlowerEdges`,
  `OtherNFlowerEdges`, `FixedLargeRingStereoEdges`, `AllBondEdges`.
- Cases 06c, 06a, 06b, 23, 24, and 25 allocate additional local edge-list
  arrays. Their case-local labels free every initialized list before later
  cases run or before the outer exit is taken.
- Source globals `iat_*[MAX_DIFF_FIXH]` and `tdhc[MAX_DIFF_FIXH]` are mutable
  scratch buffers. A Rust representation must preserve capacity, signed-short
  writes, and reuse order without introducing cross-call data races.
- Direct BNS mutations include edge `flow`, `forbidden`, vertex `st_edge.flow`,
  and `tot_st_flow`. Every failed search restores exactly the fields changed
  before that search and in the same order.
- A successful `RunBnsRestoreOnce` is followed by counter changes and, for each
  successful case block, source-ordered reconstruction through
  `MakeOneInChIOutOfStrFromINChI2`, `FillOutExtraFixedHDataRestr`, and
  `FillOutCMP2FHINCHI`.
- `CurrEdges.num_edges = 0` clears logical membership without freeing storage.
  It must not be replaced by allocation replacement.
- The final return is `ret < 0 ? ret : (tot_succes && pc2i->bHasDifference)`.
  Positive intermediate return values are not returned directly.

## Parameter audit

`num_inp`, `bHasSomeFixedH`, `pnNumRunBNS`, and `pnTotalDelta` occur only in
the active C signature. The function instead updates its local
`nNumRunBNS`; the caller-provided counter pointers are not modified. Rust must
preserve that source behavior rather than infer intended counter updates.

All other pointer parameters participate in reads, mutations, reconstruction,
or ring/BNS work. Nullability may only be modeled where the C source performs a
null test. A checked Rust pointer failure outside the valid C input contract
must remain a structured `SourceHeapError`; it must not trigger fallback
chemistry.

The prior Rust port of `MakeOneInChIOutOfStrFromINChI2` represented
`T_GROUP_INFO **` as a borrowed `Option<&mut T_GROUP_INFO>`. That representation
could not be reassigned safely across this function's 33 reconstruction calls.
It has been replaced by `SourceTGroupInfoPointer`, which preserves the three
observable source states: `NULL`, `&pStruct->One_ti`, and an unchanged external
pointer identity on pre-assignment error. Its focused test command is:

```bash
cargo test -p cosmolkit-inchi source_port__ichirvr1__makeoneinchioutofstrfrominchi2__line_5087
```

The main test harness ran one test and passed with 907 tests filtered out.

## Source-review and focused-test matrix

The complete implementation was reviewed against the following source-defined
matrix. The focused test exercises the stable externally constructible prefix
and cleanup boundary; the ordered case bodies are retained and reviewed as one
atomic operation because the not-yet-ported caller is the source component that
constructs their full internal BNS state.

- no-fixed-H early exit, no-difference exit, and the secondary empty-difference
  exit;
- initial collection of minus, plus, N-flower, S-flower, other-N-flower, bond,
  and fixed-large-ring-stereo edges;
- every ordered case 01 through 25, including 03a, 06a-e, 07a, and 21a;
- every case-local `goto` cleanup path and every outer `goto exit_function`;
- `AddToEdgeList` allocation failure at each distinct ownership phase;
- BNS test negative, zero, wrong endpoint, wrong delta-H/charge, successful
  restore, and restore-error outcomes;
- exact rollback of edge flow, forbidden bits, vertex flow, and total flow;
- all 33 reconstruction sites, including negative reconstruction return and
  nonzero refill/recomparison return;
- fixed-H disappearance after reconstruction and no-difference early exit;
- signed `char`, unsigned `char`, `short`, edge-index, vertex-index, mask
  complement, and wrapping counter behavior at source-defined boundaries;
- all local and outer edge-list allocations are released on every return path;
- final `ret` mapping and the intentional non-modification of the four unused
  parameters identified above.

The first behavior marker is `✔️`; the second remains `❌` because SourceHeap
lookup and ownership checks add known overhead compared with direct C pointer
access.

## Translation ledger

The non-production translation draft is kept outside the crate so an
incomplete chemical operation cannot be linked accidentally.

| Source range | Status | Reviewed behavior |
|---|---|---|
| 333-513 | integrated; line-reviewed | seven-list initialization, fixed-H early return, charge/flower/bond edge collection, stereo-ring collection, first canonical-map refill and comparison, both no-work exits |
| 514-647 | integrated; line-reviewed | complete case 01 classification, mask changes, BNS test/restore conditions, exact failed-test rollback, success reconstruction, refill and recomparison |
| 648-795 | integrated; line-reviewed | complete case 02 canonical/remapped classification, mobile-H pointer fallback, cumulative zero-return count, one-time N-flower release and retry, rollback, reconstruction and recomparison |
| 796-944 | integrated; line-reviewed | complete case 03 difference/canonical classification, reversed mobile-H endpoint lookup, charge-edge mask changes, BNS endpoint and delta-charge checks, exact failed-test rollback, success reconstruction, refill and recomparison |
| 945-1083 | integrated; line-reviewed | complete case 03a canonical classification, input fixed/mobile-H endpoint-or-acceptor test, reversed endpoint lookup, selective current-edge collection, BNS mutation and exact rollback, success reconstruction, refill and recomparison |
| 1084-1227 | integrated; line-reviewed | complete case 04 canonical classification with input endpoint/fixed/mobile-H state, charged versus neutral exclusion, positive-edge collection, negative delta-charge acceptance, exact rollback, success reconstruction, refill and recomparison |
| 1228-1360 | integrated; line-reviewed | complete case 05 difference classification with exclusive positive-NH and negative-O paths, current charge-edge collection, positive delta-charge acceptance, exact rollback, success reconstruction, refill and recomparison |
| 1361-1525 | integrated; line-reviewed | complete case 06c fixed-H/charge guard, current positive-edge collection, adjacent-charge discovery with duplicate suppression, case-local edge-list allocation and label cleanup, BNS mutation without a forbidden-bit change, exact failed-test rollback, success reconstruction, refill and recomparison |
| 1526-1682 | integrated; line-reviewed | complete case 06d difference classification with C integer promotion, case-local charge-edge ownership, one- or two-pass sulfur-flower masking, BNS mutation without a forbidden-bit change, exact rollback, case-label cleanup, success reconstruction, refill and recomparison |
| 1683-1815 | integrated; line-reviewed | complete case 06 mutually exclusive neutral and positive fixed-H classification, shared charge-edge collection, temporary candidate-edge forbidden bit, negative delta-charge acceptance, exact rollback, success reconstruction, refill and recomparison |
| 1816-1949 | integrated; line-reviewed | complete case 06a t-group/proton and normalized-fixed-bond guard, canonical neutral and charged classification, source-order Add-result control flow, case-local charge-edge cleanup, BNS mutation without forbidden-bit change, exact rollback, success reconstruction, refill and recomparison |
| 1950-2096 | integrated; line-reviewed | complete case 06b relaxed proton guard, source-defined chalcogen-or-period-row-7 charged classification, fixed-bond charge requirement, source-order Add control flow, one- or two-pass sulfur-flower masking, case-local cleanup, exact rollback, success reconstruction, refill and recomparison |
| 2097-2328 | integrated; line-reviewed | complete case 06e normalized fixed-bond guard and fallback expression, t-group difference fill, donor/receptor counts and endpoint lists, source-defined max-like removal count, reverse adjacency scans, three-edge/two-vertex flow mutation, delta-charge 0-or-1 acceptance, exact rollback, post-attempt donor-edge insertion, allocation-failure mask retention, local-list cleanup, success reconstruction, refill and recomparison |
| 2329-2469 | integrated; line-reviewed | complete case 07 difference and canonical negative-oxygen classification, reversed endpoint and input H exclusion, shared minus-edge collection, temporary candidate-edge forbidden bit, positive delta-charge acceptance, exact rollback, success reconstruction, refill and recomparison |
| 2470-2612 | integrated; line-reviewed | complete case 07a non-taut neutral double-bond classification and canonical negative-chalcogen classification with source-ordered single N(V) neighbor lookup, input endpoint/H exclusion, temporary minus-edge forbidden bit, positive delta-charge acceptance, exact rollback, success reconstruction, refill and recomparison |
| 2613-2762 | integrated; line-reviewed | complete case 08 single-input-t-group guard, exclusive endpoint negative-chalcogen and non-stereogenic neutral-nitrogen classification, shared minus-edge collection, temporary candidate forbidden bit, independent fixed-large-ring stereo mask release and restoration, positive delta-charge acceptance, exact rollback, success reconstruction, refill and recomparison |
| 2763-2919 | integrated; line-reviewed | complete case 09 charged nitrogen/double-bond/central-carbon/second-nitrogen nested discovery, C integer promotion, bond-edge lookup by atom neighbor order, charge-edge masking, zero-flow continue that intentionally retains masks and current edges, bond-flow mutation, zero delta-charge acceptance, unconditional post-test forbidden-bit removal, exact failed-test rollback, success reconstruction, refill and recomparison |
| 2920-3092 | integrated; line-reviewed | complete case 10 used-entry skip, neutral-N/positive-center/second-neutral-N nested search, full comparison lookup, first-then-second charge-edge selection using corresponding bond forbidden bits and charge flows, center-edge unmasking, selected charge-edge mutation without forbidden-bit change, negative delta-charge acceptance, exact failed-test rollback, positive-restore-only in-place `nValue` marking and loop termination, success reconstruction, refill and recomparison |
| 3093-3303 | integrated; line-reviewed | complete case 11 charged N/chalcogen and neutral-N canonical classification, neutral-edge-only current list, zero-flow charged-edge requirement, asymmetric reverse adjacency scans, three-edge/two-remote-vertex mutation, matching remote endpoints with delta-charge 0-or-1 acceptance, exact five-field failed-test rollback, success reconstruction, refill and recomparison |
| 3304-3518 | integrated; line-reviewed | complete case 12 source guard and minus-H boundary, N(V)-only then general retry, short-circuit `bN_V` assignment with C integer promotion, first-N(V) candidate reset and flower-edge gating, donor minus-edge permanent unmask, `iedge[0]` bond mutation, positive-two delta-charge acceptance, unconditional temporary bond-mask clearing, exact failed-test rollback, success reconstruction, refill and recomparison |
| 3519-3700 | integrated; line-reviewed | complete case 13 paired t-group comparison with C-promoted unsigned-rank subtraction, canonical/atom H-array indexing, source-defined mixed `pVA[i]` donor classification versus `pVA[iat]` receiver and BNS use, donor minus-edge permanent unmask, `iedge[0]` bond mutation, positive-two delta-charge acceptance, exact failed-test rollback, unconditional temporary bond-mask clearing, successful reconstruction and t-group-loop break |
| 3701-4032 | integrated; line-reviewed | complete case 14 short-circuit `bHas_N_V` guard, canonical reversed-endpoint indexing, local candidate-list ownership, source edge-zero exclusions, invalid-candidate marking, reverse adjacency scans excluding N(V)-plus, test #1 mask/mutation and negative-error no-rollback path, test #1 normalization and restoration, test #2 alternate charge unmasking, source-defined positive-mismatch no-rollback/mark behavior, candidate retirement, local-label cleanup, error propagation, success reconstruction, refill and recomparison |
| 4029-4205 | integrated; line-reviewed | complete case 15 fixed/mobile-H-free double-bond positive-chalcogen versus neutral single-bond nitrogen classification, neutral-N plus-edge collection, zero-flow charged plus-edge requirement, symmetric reverse adjacency scans, three-edge/two-remote-vertex mutation, matching remote endpoints with delta-charge 0-or-1 acceptance, exact five-field failed-test rollback, success reconstruction, refill and recomparison |
| 4207-4365 | integrated; line-reviewed | complete case 16 t-group endpoint traversal, source-order neutral double-bond chalcogen versus fixed-H negative nitrogen classification, current minus-edge collection, flowing nitrogen minus-edge mutation, positive-one delta-charge acceptance, source-defined forbidden-bit clearing on failed tests, exact flow rollback, per-group success reconstruction, refill and recomparison |
| 4367-4538 | integrated; line-reviewed | complete case 17 positive-H double-bond pnictogen/chalcogen versus non-endpoint neutral single-bond heteroatom classification, neutral plus-edge collection, removed-H-difference attempt cap, symmetric reverse adjacency scans, three-edge/two-remote-vertex mutation, delta-charge 0-or-1 acceptance, exact failed-test rollback, success reconstruction, refill and recomparison |
| 4540-4683 | integrated; line-reviewed | complete case 18 source-defined final qualifying t-group candidate retention, mixed canonical/atom valence indexing, non-taut positive-edge search, per-candidate charge and other-N-flower masking, plus-edge mutation, positive-one delta-charge acceptance, case-label mask cleanup and negative-return propagation, success reconstruction, refill and recomparison |
| 4685-4837 | integrated; line-reviewed | complete case 19 comparison-entry classification with source short-circuit access order, atom-number storage in the shared edge list, metal plus/minus and O-metal bond unmasking, C-promoted metal-charge calculation, expected delta-charge selection, positive-O edge mutation and rollback, source continue path mask behavior, global mask cleanup, success reconstruction, refill and recomparison |
| 4839-5009 | integrated; line-reviewed | complete case 20 primary comparison-entry donor/receiver classification, original t-group multiplicity restriction, no-primary-receiver endpoint fallback without an invented element condition, receiver minus-edge collection, flowing donor minus-edge mutation, positive-one delta-charge acceptance, exact failed-test rollback, global mask cleanup, success reconstruction, refill and recomparison |
| 5011-5217 | integrated; line-reviewed | complete case 21 comparison receiver collection, canonical endpoint donor scan, source-defined `at2[i]` center lookup and `at[iS]` loop bound, duplicate-center suppression, full tautomer/other/negative-endpoint neighbor counters, donor minus-edge masking and mutation, positive-one delta-charge acceptance, global mask cleanup, success reconstruction, refill and recomparison |
| 5219-5460 | integrated; line-reviewed | complete case 21a four-list ownership and local cleanup label, central chalcogen candidate discovery, source mixed difference-index endpoint comparison, full neighbor counters, fixed/changeable sulfur-oxygen edge partition, removal of fixed acceptors, donor minus-edge mutation and positive-one delta-charge acceptance, exact rollback, list cleanup, success reconstruction, refill and recomparison |
| 5462-5604 | integrated; line-reviewed | complete case 22 comparison receiver collection, canonical N(-)=N(+)=C chain recognition with source boolean neighbor subscript, donor minus-edge collection, target first-bond mutation, target minus-edge unmasking, zero delta-charge acceptance, exact failed-test rollback, per-candidate mask restoration, global mask cleanup, success reconstruction, refill and recomparison |
| 5606-5973 | integrated; line-reviewed | complete case 23 NO2 target discovery, four changeable-edge sets, source-order fixed/NOOH/wrong-taut/taut edge collection with duplicate fixed-edge insertion, local case cleanup, target bond mutation, ordered three-set retries without intermediate rollback, set-specific delta-charge 2-or-0 acceptance, final failed-target rollback, success reconstruction, refill and recomparison |
| 5975-6262 | integrated; line-reviewed | complete case 24 both NN candidate forms, source boolean neighbor selection and raw unshifted minus-edge avoidance, NN edge/flower/expected-delta triples, five-list ownership, missing/O/N tautomer acceptor partition, ordered three-set retries without intermediate rollback, optional flower unmasking, exact failed-target rollback, success reconstruction, refill and recomparison |
| 6264-6520 | integrated; line-reviewed | complete case 25 six-list ownership, source `bFirst` condition, canonical/atom endpoint mapping, underconstrained `nNumO != 1 && nNumOthers != 1` acceptance, N-plus/O-minus pair collection, required-list gate, ordered four-set retries with constant expected delta-charge 3, exact failed-pair rollback, success reconstruction, refill and recomparison |
| 6523-6534 | integrated; line-reviewed | shared function cleanup frees seven outer edge lists in source order and returns negative `ret` unchanged or the C logical-and result of total success and remaining difference |


The integrated function parses under `rustfmt` and is available to its source-defined caller.
