# RDKit InChI Completed-Prefix Exact-Parity Audit

## Scope And Gate

This audit covers the 18 source functions whose historical Port steps were
checked while their source frames used `INCHI❗❌` at repair entry. The selected
production configuration is the official GCC/Linux `libinchi` target with
`COMPILE_ANSI_ONLY` and `TARGET_API_LIB`. The first behavioral marker may be
upgraded only after complete source and active-macro review, a complete focused
Rust test, and an independent official C behavior oracle exact match.

`SetAtomAndBondProperties` remains frozen. Its renumbered Port step must not run
until every row below has marker `INCHI✔❌`, exact focused and oracle results,
and no unclosed branch.

The second axis remains `❌`: this repair phase does not claim removal of the
known or unresolved performance differences represented by the historical
markers.

## Baseline Findings

- All 18 source frames exist, have been reviewed against the selected active C
  macro branches, and use `INCHI✔❌`.
- All 18 functions have completed their function-specific focused Rust test and
  independent official C behavior oracle. Every recorded observable field
  compares exactly for the currently modeled C-defined input domain.
- `InchiToInchiAtom` retains the original eight official C behavior cases and
  now has 193 additional fixed cases. All 201 cases compare every emitted
  status, error, text, atom, bond, coordinate, stereo, isotope-H, charge,
  label/value, caller-capacity, pointer-state, allocation, cleanup, and input
  position field exactly. The source frame remains line-verbatim after
  ignoring source-only trailing whitespace on blank lines.
- `WriteCoord`, `parse_options_string`, and
  `el_number_in_internal_ref_table`, `get_periodic_table_number`, and
  `get_el_valence`, `detect_unusual_el_valence`, and
  `extract_charges_and_radicals`, `extract_H_atoms`, `is_el_a_metal`, and
  `get_atomic_mass_from_elnum`, `is_in_the_list`, `nBondsValToMetal`, and
  `SetAtomProperties`, `SetBondProperties`, `InchiToInchiAtom`,
  `InchiToInchi_Input`, `Get_inchi_Input_FromAuxInfo`, and
  `Get_std_inchi_Input_FromAuxInfo` have complete focused and official C oracle
  evidence.
- `cargo fmt --all -- --check` passed. `cargo test -p cosmolkit-inchi` passed
  with 62 Rust unit tests, one provenance integration test, and zero failures.

## Function Ledger

| Order | Function | Source path/line | Active macro configuration | Focused test command | Official C oracle command | Exact parity | Marker | Unclosed branches |
|---:|---|---|---|---|---|---|---|---|
| 1 | `WriteCoord` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:890` | GCC/Linux, `COMPILE_ANSI_ONLY`, `TARGET_API_LIB`; `GHI100_FIX` undefined, `SPRINTF_FLAG=2` inert, eight active leaves use libc `sprintf` | `cargo test -p cosmolkit-inchi source_port__ichimak2__writecoord__line_890` (exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__ichimak2__writecoord__line_890` (4,172 complete 64-byte records, exact pass) | exact pass | `INCHI✔❌` | none |
| 2 | `parse_options_string` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:1037` | GCC/Linux, `COMPILE_ANSI_ONLY`, `TARGET_API_LIB`; no function-local conditional branch, preprocessed body matches source | `cargo test -p cosmolkit-inchi source_port__inchi_dll__parse_options_string__line_1037` (exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__inchi_dll__parse_options_string__line_1037` (18 full-state records, exact pass) | exact pass | `INCHI✔❌` | none |
| 3 | `el_number_in_internal_ref_table` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:347` | GCC/Linux, `COMPILE_ANSI_ONLY`, `TARGET_API_LIB`; `INCHI_ZFRAG` undefined gives 122 active rows plus terminal empty sentinel | `cargo test -p cosmolkit-inchi source_port__util__el_number_in_internal_ref_table__line_347` (exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__util__el_number_in_internal_ref_table__line_347` (718 records exhaustively discovering all 122 indices, exact pass) | exact pass | `INCHI✔❌` | none |
| 4 | `get_periodic_table_number` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:364` | GCC/Linux, `COMPILE_ANSI_ONLY`, `TARGET_API_LIB`; active single-letter macros plus 122-row lookup, `INCHI_ZFRAG` undefined | `cargo test -p cosmolkit-inchi source_port__util__get_periodic_table_number__line_364` (exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__util__get_periodic_table_number__line_364` (719 exhaustive records, exact pass) | exact pass | `INCHI✔❌` | none |
| 5 | `get_el_valence` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:439` | GCC/Linux, `COMPILE_ANSI_ONLY`, `TARGET_API_LIB`; no function-local conditional branch; `INCHI_ZFRAG` undefined; oracle covers the H row selected by every `nPeriodicNum <= 1`, active periodic rows 2..120, terminal sentinel 121, charges -3..3, and slots 0..5 | `cargo test -p cosmolkit-inchi source_port__util__get_el_valence__line_439` (exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__util__get_el_valence__line_439` (5,208 complete C records, exact pass) | exact pass | `INCHI✔❌` | none |
| 6 | `detect_unusual_el_valence` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:620` | GCC/Linux, `COMPILE_ANSI_ONLY`, `TARGET_API_LIB`; no function-local conditional branch; direct callee `get_el_valence` is closed; oracle covers every legal table selector, charges -3..3, all radical branches/default boundaries, chemical valences -2..16, all early returns, and nonoverflowing `i32` sum boundaries | `cargo test -p cosmolkit-inchi source_port__util__detect_unusual_el_valence__line_620` (exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__util__detect_unusual_el_valence__line_620` (131,956 complete C records, exact pass) | exact pass | `INCHI✔❌` | none; C signed-addition overflow is outside the source-defined input domain and Rust reports it structurally |
| 7 | `extract_charges_and_radicals` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:700` | GCC/Linux, `COMPILE_ANSI_ONLY`, `TARGET_API_LIB`; no function-local conditional branch; active glibc LP64 `strtol`, `strpbrk`, `strrchr`, `strlen`, and `memmove`; oracle covers 52 fixed edge cases plus all 1..6-length sign runs with eight decimal/malformed suffixes | `cargo test -p cosmolkit-inchi source_port__util__extract_charges_and_radicals__line_700` (exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__util__extract_charges_and_radicals__line_700` (1,060 full-state C records, exact pass) | exact pass | `INCHI✔❌` | none; null pointers, missing NUL, and subsequent C signed-arithmetic overflow are outside the C-defined domain and remain structured Rust errors |
| 8 | `extract_H_atoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:774` | GCC/Linux, `COMPILE_ANSI_ONLY`, `TARGET_API_LIB`; no function-local conditional branch; active glibc `islower`, `isdigit`, LP64 `strtol`, `strlen`, and `memmove`; oracle covers H/D/T token classes, lowercase exclusions, repeated/mixed tokens, decimal stops and narrowing, `S_CHAR` wrapping, alias correction, and full-buffer mutation | `cargo test -p cosmolkit-inchi source_port__util__extract_h_atoms__line_774` (exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__util__extract_h_atoms__line_774` (56 full-state C records, exact pass) | exact pass | `INCHI✔❌` | none; null pointers, missing NUL, undersized isotope storage, and subsequent C signed-arithmetic overflow are outside the C-defined domain and remain structured Rust errors |
| 9 | `is_el_a_metal` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:688` | GCC/Linux, `COMPILE_ANSI_ONLY`, `TARGET_API_LIB`; no function-local conditional branch; `INCHI_ZFRAG` undefined; oracle covers every C-defined selector -1..121, including H/D/T mapping, pseudo-elements, and terminal sentinel | `cargo test -p cosmolkit-inchi source_port__util__is_el_a_metal__line_688` (exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__util__is_el_a_metal__line_688` (123 complete C records, exact pass) | exact pass | `INCHI✔❌` | none; selectors outside -1..121 are C out-of-bounds/overflow undefined behavior and remain structured Rust errors |
| 10 | `get_atomic_mass_from_elnum` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1007` | GCC/Linux, `COMPILE_ANSI_ONLY`, `TARGET_API_LIB`; no function-local conditional branch; `INCHI_ZFRAG` undefined; oracle covers every atomic number -10..130 plus the two extreme nonoverflowing arithmetic boundaries; the source deliberately bypasses D/T rows | `cargo test -p cosmolkit-inchi source_port__util__get_atomic_mass_from_elnum__line_1007` (exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__util__get_atomic_mass_from_elnum__line_1007` (143 complete C records, exact pass) | exact pass | `INCHI✔❌` | none; `INT_MIN` and `INT_MAX` invoke C signed overflow and remain structured Rust errors |
| 11 | `is_in_the_list` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1059` | GCC/Linux, `COMPILE_ANSI_ONLY`, `TARGET_API_LIB`; no function-local conditional branch; oracle exhausts lengths 0..5 over a three-value `AT_NUMB` alphabet, four targets, and every legal prefix length, plus NULL with zero length | `cargo test -p cosmolkit-inchi source_port__util__is_in_the_list__line_1059` (exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__util__is_in_the_list__line_1059` (8,021 pointer-result C records, exact pass) | exact pass | `INCHI✔❌` | none; negative lengths, NULL with positive length, and lengths beyond allocation are C undefined behavior and remain structured Rust errors |
| 12 | `nBondsValToMetal` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1100` | GCC/Linux, `COMPILE_ANSI_ONLY`, `TARGET_API_LIB`; no function-local conditional branch; direct callee `is_el_a_metal` is closed; oracle exhausts element selectors 0..121 against bond codes 0..255 and covers negative/zero/MAXVAL, mixed neighbors, accumulation, and early invalid-order returns | `cargo test -p cosmolkit-inchi source_port__util__nbondsvaltometal__line_1100` (exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__util__nbondsvaltometal__line_1100` (31,239 complete C records, exact pass) | exact pass | `INCHI✔❌` | none; null/base/index, invalid neighbor, valence beyond MAXVAL, and element selector beyond sentinel are C undefined behavior and remain structured Rust errors |
| 13 | `SetAtomProperties` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:1139` | GCC/Linux, `COMPILE_ANSI_ONLY`, `TARGET_API_LIB`; `SINGLET_IS_TRIPLET == 1`; direct `WriteCoord`, `fabs`, `strcpy`, `memcpy`, `sprintf`, `TREAT_ERR`/`AddErrorMessage`; oracle exhausts all 256 radical values and covers coordinate threshold branches, signed zero, NaN/Inf, optional coordinate storage, index 1, element/charge limits, seeded errors, full atom state, and mutation order | `cargo test -p cosmolkit-inchi source_port__inchi_dll__setatomproperties__line_1139` (exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__inchi_dll__setatomproperties__line_1139` (264 full-state C records, exact pass) | exact pass | `INCHI✔❌` | none; null/invalid pointers, unterminated element text, and atom-index arithmetic overflow are C undefined behavior and remain structured Rust errors |
| 14 | `SetBondProperties` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:1235` | GCC/Linux, `COMPILE_ANSI_ONLY`, `TARGET_API_LIB`; no function-local conditional branch; direct `is_in_the_list`, `sprintf`, `TREAT_ERR`/`AddErrorMessage`; oracle exhausts all 256 bond-type and 256 stereo byte values plus 15 nonexistent/self/duplicate/multiple/one-sided/MAXVAL/partial-mutation states | `cargo test -p cosmolkit-inchi source_port__inchi_dll__setbondproperties__line_1235` (exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__inchi_dll__setbondproperties__line_1235` (527 full-state C records, exact pass) | exact pass | `INCHI✔❌` | none; null/invalid pointers and indices outside allocated arrays are C undefined behavior and remain structured Rust errors |
| 15 | `InchiToInchiAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:582` | GCC/Linux, `COMPILE_ANSI_ONLY`, `TARGET_API_LIB`; no function-local conditional compilation; active `inchi_malloc -> malloc`, `inchi_calloc -> calloc`, and macro `inchi_free`; all 2,518 source lines present verbatim modulo unrepresentable trailing whitespace on blank source lines | `cargo test -p cosmolkit-inchi source_port__inchi_dll_b__inchitoinchiatom__line_582` (exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle_inchi_to_inchi_atom_exact` (201 complete full-state records, exact pass) | exact pass | `INCHI✔❌` | none; all remaining untaken compiler edges and unexecuted defensive lines are proven unreachable below |
| 16 | `InchiToInchi_Input` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ichilnct.c:204` | GCC/Linux, `COMPILE_ANSI_ONLY`, `TARGET_API_LIB`; function lines 204-403 have no conditional-compilation split; active allocation macros are the closed `CreateInchiAtom`/`CreateInchi_Stereo0D` callees and `inchi_free`; all 200 source lines are verbatim | `cargo test -p cosmolkit-inchi source_port__ichilnct__inchitoinchi_input__line_204` (exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__ichilnct__inchitoinchi_input__line_204` (50 full-state records, exact pass) | exact pass | `INCHI✔❌` | none; all remaining zero compiler edges and two unexecuted defensive cleanup lines are proven unreachable below |
| 17 | `Get_inchi_Input_FromAuxInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ichilnct.c:89` | GCC/Linux, `COMPILE_ANSI_ONLY`, `TARGET_API_LIB`; lines 89-180 have no function-local conditional compilation; direct callees `inchi_ios_init`, `InchiToInchi_Input`, and `Free_inchi_Input` are closed; all 92 lines are verbatim; stack `INCHI_IOSTREAM`/label/value storage is excluded from C allocation ordinals | `cargo test -p cosmolkit-inchi source_port__ichilnct__get_inchi_input_fromauxinfo__line_89` (exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__ichilnct__get_inchi_input_fromauxinfo__line_89` (28 full-state C records, exact pass) | exact pass | `INCHI✔❌` | none in the C-defined target domain; `err == 9` is unreachable through fixed `INPUT_INCHI_PLAIN`, and `err == 0 && szErrMsg[0]` contradicts the closed plain-parser error-write invariant; both remain represented in Rust return mapping |
| 18 | `Get_std_inchi_Input_FromAuxInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ichilnct.c:79` | GCC/Linux, `COMPILE_ANSI_ONLY`, `TARGET_API_LIB`; lines 79-87 have no function-local conditional compilation; all 9 lines are verbatim; sole direct callee `Get_inchi_Input_FromAuxInfo` is closed and the wrapper fixes `bDiffUnkUndfStereo` to zero | `cargo test -p cosmolkit-inchi source_port__ichilnct__get_std_inchi_input_fromauxinfo__line_79` (exact pass) | `cargo test -p cosmolkit-inchi official_c_oracle__ichilnct__get_std_inchi_input_fromauxinfo__line_79` (27 full-state C records, exact pass) | exact pass | `INCHI✔❌` | none in the C-defined target domain; every observable delegated branch is represented by the standard-wrapper oracle suite |

## `InchiToInchiAtom` Defensive-Branch Closure

The official C coverage build used `--coverage -O0 -fno-inline`,
`COMPILE_ANSI_ONLY`, `TARGET_API_LIB`, and the GCC/Linux libinchi objects. Its
final `target/inchi-official-c-coverage/report/inchi_dll_b.c.gcov` result is
87.59% of 1,458 lines executed, 94.78% of 882 branches executed, and 87.98% of
882 branches taken at least once. Coverage below 100% is not being treated as
proof by percentage; every remaining source line and zero compiler edge was
classified against the active source:

- Lines 1351-1359 and 2605-2614: `strtod(p, &q)` receives a non-NULL `endptr`.
  The C library contract always stores a pointer in `q`, including conversion
  failure, so the `q == NULL` body is unreachable. The same contract closes
  the `q == NULL` short-circuit edges after `strtol` at lines 679, 812, 1891,
  and 2008.
- Lines 1363 and 2618: `k` is controlled by `for (k = 0; k < 3; k++)`; the
  coordinate switch default edge is unreachable.
- Lines 1438-1441 and 2697-2700: every accepted bond character is mapped to a
  value in `MIN_INPUT_BOND_TYPE..=MAX_INPUT_BOND_TYPE`; unknown characters
  leave through the earlier input-error path. The post-parse range repair is
  therefore unreachable for source-defined parser state.
- Lines 1486 and 1493-1496, plus 2745 and 2752-2755: each parsed bond writes
  both adjacency directions before aromatic postprocessing. A reciprocal
  aromatic neighbor cannot disappear without prior out-of-bounds undefined
  behavior, which is outside the C-defined input domain.
- Lines 1593-1601 and 2859-2867: double-bond stereo records are created only
  from a bond immediately inserted in both directions, so both reciprocal
  lookups exist. The `INCHI_StereoType_None` mutation is defensive and
  unreachable.
- Lines 1667-1668 and 1691-1696, plus 2929-2930 and 2955-2960: entering the
  one-endpoint cumulene branch initializes a two-node chain and advances once
  before its terminal lookup, so `len > 2`; every advance chooses the other
  neighbor while retaining `prev`, which remains a neighbor of `cur`.
- Lines 1571, 1576, 1711-1712, 2829, 2834, and 2978-2979: stereo rows are
  constructed only as tetrahedral or double-bond rows and only for nonzero
  parity. `None`, zero parity, and negative stereobond orders require one of
  the already-proven-unreachable defensive failures above.
- Lines 885, 1008, 2073, 2118, and 2228 expose zero branches from the active
  `inchi_free` macro. Each occurs only after the corresponding allocation has
  succeeded, so the macro's NULL arm is unreachable in that block.
- Lines 1002, 1040, 1260, 1312, 1394, 1820, 1982, 2093, 2222, 2499, 2649,
  2814, 3058, and 3061 contain compiler-expanded short-circuit edges. Their
  source-level outcomes are covered, including no terminator, trailing data,
  premature section end, cross-buffer input, NULL outputs, and caller-owned
  cleanup. Remaining zero edges require contradictory `bItemIsOver`/`s`/`p`,
  section-length, allocation, or ownership state established immediately
  before each condition.
- Lines 1612, 1622, 2878, and 2887 are nested short-circuit edges reached only
  after one endpoint is already a degree-two double-bond middle atom. A parsed
  stereobond's opposite endpoint supplies its non-stereobond neighbor; the
  omitted degree combination contradicts that entry state. Both middle-atom
  orientations, isotope-H rejection, and later-chain isotope-H rejection are
  covered in plain and XML cases.
- Lines 1757 and 3025 are success cleanup after allocation, so `atom == NULL`
  is unreachable there. Lines 1774-1776 are after natural termination of the
  plain outer scanner; every path that allocated stereo state jumps to
  `bypass_end_of_INChI_plain`, so natural termination has no stereo allocation.
  Lines 1782 and 3045 follow loops whose only fallthrough condition is
  `res <= 0`. Lines 1792 and 3056 are labels, not executable behavior.

No source-defined active branch remains unclosed for this function. NULL or
invalid required pointers, unterminated caller buffers, and states requiring
prior C out-of-bounds access remain outside the C-defined domain and are
reported structurally by Rust.

## `InchiToInchi_Input` Defensive-Branch Closure

The independent oracle has 18 fixed semantic records and two 16-record
allocation-ordinal sweeps, for 50 full-state records total. It compares return
status, error code/text, stream position, label/value/id/flags, complete atom,
bond, coordinate, and stereo fields, atom/stereo/options pointer state, old
allocation release state, output reset, and exact C allocation-call count.
The fixed records cover string/file input, single-structure stop, merge through
EOF, existing atom/stereo replacement, atom/bond and stereo index remapping,
count-only output, empty input, malformed input, MAX_ATOMS, seeded error,
warning text, fatal cleanup, XML input, NULL optional outputs, zero-count
non-NULL old storage, and double-bond stereo with negative central index.

The GCC coverage object was compiled directly from `ichilnct.c` with
`--coverage -O0 -fno-inline`, `COMPILE_ANSI_ONLY`, and `TARGET_API_LIB`, then
linked against the official CMake target objects. The report is
`target/inchi-official-c-input-coverage/report/ichilnct.c.gcov`. Every
executable line in the active function body is covered except lines 378 and
382, and every remaining zero edge is source-invariant unreachable:

- Line 257 exact-zero EOF guard: whenever `InchiToInchiAtom` returns exactly
  zero with a positive error after the prior terms pass, its error is in
  11..19. Non-EOF parser errors return a negative status, so the `err >= 20`
  compiler edge is unreachable. Reaching the final merge term with prior
  atoms also proves the loop was entered with merge enabled, so its false edge
  is contradictory.
- Line 305: a positive parser return with atom output requested guarantees
  non-NULL `at_new`; allocation failure returns a fatal negative status.
- Line 308: positive `num_inp_0D_new` is produced only with non-NULL
  `stereo0D_new`; allocation failure resets the count and returns fatal.
- Lines 376-378: `szCoordNew` is initialized to NULL and never assigned in the
  active function, so its cleanup body cannot execute.
- Lines 380-382: every loop path adopts, frees, or NULLs `at_new` before the
  loop exits; the post-loop non-NULL cleanup is defensive and unreachable.

No source-defined active branch remains unclosed for this function. The Rust
model's stack/static temporary storage is excluded from C allocator-failure
ordinals; only source `malloc/calloc` calls participate, matching the official
C oracle exactly.

## Current Gate Result

The completed-prefix repair gate is closed for all 18 functions. Every ledger
row has a complete source frame, selected-macro review, passing focused Rust
test, passing independent official C behavior oracle with exact field
comparison, final `INCHI✔❌` marker, and no unclosed source-defined branch in
the currently modeled input domain. The second marker remains `❌` because the
known performance differences have not been removed.

This result does not claim exact parity for the broader active InChI call graph.
The frozen suffix may resume only after the mandatory protocol reading in Step
323; `SetAtomAndBondProperties` remains unexecuted at Step 324 until then.
