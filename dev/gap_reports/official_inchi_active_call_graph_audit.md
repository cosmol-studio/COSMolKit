# Official InChI Active Call Graph Audit

## Audit Basis

- Generated: `2026-07-19`; callback-closure correction: `2026-07-21`.
- Selected target: official `libinchi` shared-library CMake target from `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/CMakeLists.txt`.
- Platform branch: GCC 13.3.0 on Linux, C11.
- Active target definitions: `COMPILE_ANSI_ONLY`, `TARGET_API_LIB`, `fPIC`.
- Compiler evidence: the target's 60-entry `compile_commands.json`, GCC `-aux-info`, GCC `-fdump-ipa-cgraph`, and the linked ELF version-script exports.
- Raw tree-sitter inventory is not used to decide activity or call edges.
- Preprocessed active vendored function definitions: `1364`.
- Export-root reachable definitions through compiler direct edges plus compiler-verified callback-target edges: `1274`.
- Active definitions not reachable from an exported API after callback closure: `90`.
- Exported production roots: `127`.
- Compiler direct source edges: `3052`.
- Compiler indirect-call sites: `31`.
- Compiler-verified standard-library callback-target edges added to reachability closure: `2`.

Activity means that GCC saw a function body after preprocessing under the selected target configuration. Reachability is the recursive closure of the ELF version-script exports over compiler-recorded direct edges plus callback targets whose function addresses and callback arguments are both present in GCC's configured tree and IPA output. Unresolved indirect callback sites remain separate and are not converted into guessed name-based edges.

## Required Classification

- **Active source functions:** listed in the active-definition appendix below; each has a GCC `NF/OF/IF` definition record and an IPA function body.
- **Inactive conditional functions:** raw bodies excluded by preprocessing have no GCC definition/body and cannot receive a production Port step.
- **Macro-only behavior:** allocation names under this target are preprocessor aliases, not C function bodies; they are documented in the allocation section below.
- **Header-defined behavior:** active header bodies are included only when GCC emits an active vendored body; plain macros are represented as behavior attached to their caller.
- **External declarations:** `NC` declaration records without a matching active vendored definition are not functions to port.
- **libc/compiler intrinsics:** unresolved compiler edges such as `malloc`, `free`, `memcpy`, ctype accessors, and checked builtins are runtime semantics, not vendored Port steps.
- **CLI/demo and non-production behavior:** a compiled body outside the exported-root closure is active but is classified as target-unreachable and does not justify a public API.

## `SortAndPrintINChI` Callback-Closure Correction

The original extraction treated only ordinary call expressions as call-graph edges. That omitted callbacks passed as function-pointer arguments to known standard-library algorithms. Under the selected GCC/Linux `COMPILE_ANSI_ONLY` and `TARGET_API_LIB` configuration, this omission incorrectly classified four comparator functions as export-unreachable and left the `SortAndPrintINChI` closure incomplete.

The correction is compiler-verified rather than name-inferred:

- GCC `-fdump-tree-original` for configured `runichi4.c` contains active calls at source lines `236-245` equivalent to `qsort(..., CompINChINonTaut2)` and `qsort(..., CompINChITaut2)`.
- GCC `-fdump-ipa-cgraph` records both comparator symbols as address-taken references from `SortAndPrintINChI` and records both `qsort` calls from that function.
- GCC configured trees for `ichimake.c` contain complete active bodies for both wrappers. Because `mode.h:542` defines `CANON_FIXH_TRANS` as `1`, each wrapper contains both calls to `CompINChI2` before the stable `ord_number` tie-break.
- GCC IPA output records `CompINChI2` calling `CompareHillFormulas`, `CompareHillFormulasNoH`, `CompareTautNonIsoPartOfINChI`, and `CompareInchiStereo`; it records `CompareHillFormulas` calling `GetElementAndCount` plus libc `strcmp`.

The corrected vendored closure introduced before `SortAndPrintINChI` is therefore:

```text
SortAndPrintINChI
|- qsort callback -> CompINChINonTaut2
|  `- CompINChI2
|     |- CompareHillFormulas -> GetElementAndCount
|     |- CompareHillFormulasNoH -> GetElementAndCount
|     |- CompareTautNonIsoPartOfINChI
|     `- CompareInchiStereo
`- qsort callback -> CompINChITaut2
   `- CompINChI2 (same closure)
```

`GetElementAndCount`, `CompareHillFormulasNoH`, `CompareTautNonIsoPartOfINChI`, and `CompareInchiStereo` were already export-root reachable and already precede this point in the port plan. `CompareHillFormulas`, `CompINChI2`, `CompINChINonTaut2`, and `CompINChITaut2` must be ported before `SortAndPrintINChI`; they are production-reachable through the configured serialization export path `INCHIGEN_DoSerialization -> SortAndPrintINChI`.

## Allocation Branch

Under `TARGET_API_LIB` on GCC/Linux, `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mode.h:1172-1181` selects:

```c
#define inchi_malloc   malloc
#define inchi_calloc   calloc
#define inchi_realloc  realloc
#define inchi_free(X)  do{ if(X) free(X); }while(0)
```

Therefore `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1552`, `:1561`, and `:1570` are preprocessed out. In particular, the raw self-referential body at `util.c:1570` is not the active GCC/Linux `TARGET_API_LIB` implementation. Rust must reproduce libc allocation sizing/zeroing/reallocation and the null-checked `free` macro behavior, anchored to `mode.h`, without inventing active C functions named `inchi_malloc`, `inchi_calloc`, `inchi_realloc`, or `inchi_free`.

The same rule applies to `qmalloc`, `qfree`, `qzfree`, `fast_alloc`, and `fast_free`: they are macro-expanded behavior. Under `USE_ALLOCA == 0` and GCC, they ultimately use the active allocation macros above; no standalone production function exists for them.

## Required Function Checks

| Name | Configured status | In `InchiToInchiAtom` closure | Definition |
|---|---|---|---|
| `is_in_the_slist` | export-reachable | yes | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:531` |
| `is_element_a_metal` | export-reachable | yes | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:543` |
| `InchiToInchiAtom` | export-reachable | yes | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:582` |
| `FindToken` | export-reachable | yes | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:445` |
| `LoadLine` | export-reachable | yes | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:489` |
| `inchi_ios_gets` | export-reachable | yes | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:386` |
| `inchi_ios_getsTab` | export-reachable | yes | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:420` |
| `inchi_ios_getsTab1` | export-reachable | yes | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:451` |
| `lrtrim` | export-reachable | yes | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1804` |
| `mystrncpy` | export-reachable | yes | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1760` |
| `CreateInchiAtom` | export-reachable | yes | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:117` |
| `CreateInchi_Stereo0D` | export-reachable | yes | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:118` |
| `FreeInchi_Stereo0D` | export-reachable | yes | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:125` |

Both `is_in_the_slist` and `is_element_a_metal` are active GCC function bodies and direct callees of `InchiToInchiAtom`; they must be ported before that caller.

## `InchiToInchiAtom` Direct Callees

| Callee | Definition | Classification |
|---|---|---|
| `CreateInchiAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:117` | vendored active definition |
| `is_in_the_slist` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:531` | vendored active definition |
| `is_element_a_metal` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:543` | vendored active definition |
| `inchi_ios_gets` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:386` | vendored active definition |
| `inchi_ios_getsTab` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:420` | vendored active definition |
| `inchi_ios_getsTab1` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:451` | vendored active definition |
| `AddErrorMessage` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichierr.c:106` | vendored active definition |
| `CreateInchi_Stereo0D` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:118` | vendored active definition |
| `FreeInchi_Stereo0D` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:125` | vendored active definition |
| `FindToken` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:445` | vendored active definition |
| `LoadLine` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:489` | vendored active definition |
| `mystrncpy` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1760` | vendored active definition |
| `lrtrim` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1804` | vendored active definition |
| `__ctype_b_loc` | runtime | libc/compiler intrinsic or macro expansion |
| `free` | runtime | libc/compiler intrinsic or macro expansion |
| `malloc` | runtime | libc/compiler intrinsic or macro expansion |
| `memcmp` | runtime | libc/compiler intrinsic or macro expansion |
| `memcpy` | runtime | libc/compiler intrinsic or macro expansion |
| `memmove` | runtime | libc/compiler intrinsic or macro expansion |
| `memset` | runtime | libc/compiler intrinsic or macro expansion |
| `sprintf` | runtime | libc/compiler intrinsic or macro expansion |
| `strchr` | runtime | libc/compiler intrinsic or macro expansion |
| `strlen` | runtime | libc/compiler intrinsic or macro expansion |
| `strspn` | runtime | libc/compiler intrinsic or macro expansion |
| `strstr` | runtime | libc/compiler intrinsic or macro expansion |
| `strtod` | runtime | libc/compiler intrinsic or macro expansion |
| `strtol` | runtime | libc/compiler intrinsic or macro expansion |

## Recursive Transitive Callees

| Minimum depth | Function | Definition | Direct vendored callees |
|---|---|---|---|
| 2 | `inchi_ios_str_gets` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:324` | `inchi_ios_str_getc` |
| 2 | `inchi_ios_str_getsTab` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:354` | `inchi_ios_str_getc` |
| 2 | `already_have_this_message` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichierr.c:160` | none |
| 3 | `inchi_ios_str_getc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:294` | none |

The complete non-vendored closure leaves are:

`__ctype_b_loc`, `calloc`, `ferror`, `fgetc`, `free`, `malloc`, `memchr`, `memcmp`, `memcpy`, `memmove`, `memset`, `sprintf`, `strcat`, `strchr`, `strlen`, `strrchr`, `strspn`, `strstr`, `strtod`, `strtol`.

No function-name inference was used for these edges. The direct and transitive tables are generated from the configured GCC IPA graph; unresolved indirect calls remain explicitly indirect.

## Active Definition Appendix

| Function | Source | Linkage | Production classification | Direct active callees | Indirect sites |
|---|---|---|---|---|---|
| `Get_std_inchi_Input_FromAuxInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ichilnct.c:79` | external | export root | `Get_inchi_Input_FromAuxInfo` | 0 |
| `Get_inchi_Input_FromAuxInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ichilnct.c:89` | external | export root | `Free_inchi_Input`, `InchiToInchi_Input`, `inchi_ios_init` | 0 |
| `Free_std_inchi_Input` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ichilnct.c:182` | external | export root | `Free_inchi_Input` | 0 |
| `Free_inchi_Input` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ichilnct.c:187` | external | export root | `FreeInchi_Atom`, `FreeInchi_Stereo0D` | 0 |
| `InchiToInchi_Input` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ichilnct.c:204` | external | reachable | `FreeInchi_Atom`, `CreateInchiAtom`, `FreeInchi_Input`, `InchiToInchiAtom`, `AddErrorMessage`, `CreateInchi_Stereo0D`, `FreeInchi_Stereo0D` | 0 |
| `FreeINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:153` | external | export root | none | 0 |
| `FreeStdINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:183` | external | export root | `FreeINCHI` | 0 |
| `FreeStructFromStdINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:195` | external | export root | `FreeStructFromINCHI` | 0 |
| `FreeStructFromINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:208` | external | export root | none | 0 |
| `GetStdINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:242` | external | export root | `input_erroneously_contains_pseudoatoms`, `GetINCHI1` | 0 |
| `GetINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:270` | external | export root | `input_erroneously_contains_pseudoatoms`, `GetINCHI1` | 0 |
| `input_erroneously_contains_pseudoatoms` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:293` | external | reachable | none | 0 |
| `GetINCHIEx` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:325` | external | export root | `GetINCHI1` | 0 |
| `GetINCHI1` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:345` | internal | reachable | `produce_generation_output`, `parse_options_string`, `ExtractOneStructure`, `inchi_ios_init`, `inchi_ios_close`, `inchi_ios_eprint`, `inchi_strbuf_init`, `inchi_strbuf_close`, `SetBitFree`, `ReadCommandLineParms`, `PrintInputParms`, `HelpCommandLineParms`, `FreeOrigAtData`, `ProcessOneStructureEx`, `FreeAllINChIArrays`, `inchi_stricmp` | 0 |
| `produce_generation_output` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:741` | external | reachable | `copy_corrected_log_tail` | 0 |
| `copy_corrected_log_tail` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:787` | external | reachable | none | 0 |
| `CheckINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:834` | external | export root | `GetINCHIfromINCHI`, `extract_inchi_substring` | 0 |
| `SetNumImplicitH` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:994` | external | reachable | `is_el_a_metal`, `get_num_H` | 0 |
| `parse_options_string` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:1037` | external | reachable | none | 0 |
| `SetAtomProperties` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:1139` | external | reachable | `AddErrorMessage`, `WriteCoord` | 0 |
| `SetBondProperties` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:1235` | external | reachable | `AddErrorMessage`, `is_in_the_list` | 0 |
| `SetAtomAndBondProperties` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:1439` | external | reachable | `AddErrorMessage`, `get_periodic_table_number`, `detect_unusual_el_valence`, `extract_charges_and_radicals`, `extract_H_atoms`, `get_atomic_mass_from_elnum`, `nBondsValToMetal`, `mystrncpy` | 0 |
| `InpAtom0DToInchiAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:1670` | external | reachable | none | 0 |
| `ExtractOneStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:1863` | external | reachable | `SetNumImplicitH`, `SetAtomProperties`, `SetBondProperties`, `SetAtomAndBondProperties`, `SetExtOrigAtDataByInChIExtInput`, `AddErrorMessage`, `FreeOrigAtData`, `Extract0DParities`, `TreatErrorsInReadTheStructure` | 0 |
| `GetStringLength` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:2086` | external | export root | none | 0 |
| `GetINCHIfromINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:2114` | external | export root | `CheckINCHI`, `parse_options_string`, `inchi_ios_init`, `inchi_ios_close`, `inchi_ios_reset`, `inchi_ios_eprint`, `SetBitFree`, `ReadCommandLineParms`, `PrintInputParms`, `HelpCommandLineParms`, `ReadWriteInChI`, `inchi_stricmp` | 0 |
| `GetStructFromStdINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:2461` | external | export root | `GetStructFromINCHI` | 0 |
| `GetStructFromINCHIEx` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:2485` | external | export root | `CheckINCHI`, `parse_options_string`, `InpAtom0DToInchiAtom`, `SetInChIExtInputByExtOrigAtData`, `inchi_ios_init`, `inchi_ios_close`, `inchi_ios_print_nodisplay`, `inchi_ios_eprint`, `SetBitFree`, `ReadCommandLineParms`, `PrintInputParms`, `HelpCommandLineParms`, `ReadWriteInChI`, `FreeExtOrigAtData`, `OAD_Polymer_SmartReopenCyclizedUnits`, `inchi_stricmp` | 0 |
| `GetStructFromINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:2856` | external | export root | `GetStructFromINCHIEx` | 0 |
| `FreeStructFromINCHIEx` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:2892` | external | export root | `FreeInChIExtInput` | 0 |
| `FreeInChIExtInput` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:2923` | external | reachable | none | 0 |
| `SetExtOrigAtDataByInChIExtInput` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3021` | external | reachable | `FreeExtOrigAtData` | 0 |
| `SetInChIExtInputByExtOrigAtData` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3258` | external | reachable | `FreeInChIExtInput` | 0 |
| `STDINCHIGEN_Create` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:122` | external | export root | `INCHIGEN_Create` | 0 |
| `INCHIGEN_Create` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:130` | external | export root | `inchi_ios_init`, `inchi_strbuf_init` | 0 |
| `STDINCHIGEN_Setup` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:185` | external | export root | `INCHIGEN_Setup` | 0 |
| `INCHIGEN_Setup` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:260` | external | export root | `input_erroneously_contains_pseudoatoms`, `parse_options_string`, `ExtractOneStructure`, `AddErrorMessage`, `ReadCommandLineParms`, `PrintInputParms`, `HelpCommandLineParms`, `inchi_stricmp` | 0 |
| `STDINCHIGEN_DoNormalization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:435` | external | export root | `INCHIGEN_DoNormalization` | 0 |
| `INCHIGEN_DoNormalization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:442` | external | export root | `NormOneStructureINChI`, `make_norm_atoms_from_inp_atoms`, `inchi_ios_init`, `AddErrorMessage`, `OrigStruct_FillOut`, `OrigAtData_WriteToSDfile` | 0 |
| `STDINCHIGEN_DoCanonicalization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:689` | external | export root | `INCHIGEN_DoCanonicalization` | 0 |
| `INCHIGEN_DoCanonicalization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:697` | external | export root | `CanonOneStructureINChI`, `inchi_ios_init`, `AddErrorMessage`, `bIsStructChiral`, `TreatCreateINChIWarning` | 0 |
| `STDINCHIGEN_DoSerialization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:859` | external | export root | `INCHIGEN_DoSerialization` | 0 |
| `INCHIGEN_DoSerialization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:869` | external | export root | `inchi_ios_init`, `SetBitFree`, `AddErrorMessage`, `FreeCompAtomData`, `MolfileSaveCopy`, `SortAndPrintINChI` | 0 |
| `STDINCHIGEN_Reset` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1125` | external | export root | `INCHIGEN_Reset` | 0 |
| `INCHIGEN_Reset` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1134` | external | export root | `inchi_ios_init`, `inchi_ios_close`, `inchi_strbuf_reset`, `OrigStruct_Free`, `free_t_group_info`, `FreeInpAtomData`, `FreeCompAtomData`, `FreeOrigAtData`, `FreeAllINChIArrays` | 0 |
| `STDINCHIGEN_Destroy` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1308` | external | export root | `INCHIGEN_Destroy` | 0 |
| `INCHIGEN_Destroy` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1315` | external | export root | `inchi_ios_close`, `inchi_strbuf_close` | 0 |
| `NormOneStructureINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:188` | external | reachable | `NormOneComponentINChI`, `CreateCompositeNormAtom`, `inchi_ios_init`, `InchiTimeGet`, `InchiTimeElapsed`, `AddErrorMessage`, `FreeInpAtomData`, `GetOneComponent`, `TreatErrorsInReadTheStructure`, `PreprocessOneStructure`, `TreatErrorsInCreateOneComponentINChI` | 0 |
| `CanonOneStructureINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:503` | external | reachable | `CanonOneComponentINChI`, `inchi_ios_init`, `InchiTimeGet`, `InchiTimeElapsed`, `FreeInpAtomData`, `GetOneComponent`, `TreatErrorsInCreateOneComponentINChI` | 0 |
| `NormOneComponentINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:654` | external | reachable | `Normalization_step`, `InchiTimeGet`, `InchiTimeElapsed`, `InchiTimeAddMsec`, `FreeInpAtomData`, `CreateInpAtomData`, `GetProcessingWarningsOneComponentInChI`, `SetConnectedComponentNumber`, `Alloc_INChI`, `Alloc_INChI_Aux` | 0 |
| `CanonOneComponentINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:923` | external | reachable | `Canonicalization_step`, `InchiTimeGet`, `InchiTimeElapsed`, `InchiTimeAddMsec`, `GetProcessingWarningsOneComponentInChI`, `SetConnectedComponentNumber` | 0 |
| `Normalization_step` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:1138` | external | reachable | `mark_alt_bonds_and_taut_groups`, `GetCanonLengths`, `set_atom_iso_sort_keys`, `inp2spATOM`, `MarkRingSystemsInp`, `set_stereo_parity`, `make_a_copy_of_t_group_info`, `set_tautomer_iso_sort_keys`, `CountTautomerGroups`, `FreeInpAtom`, `add_DT_to_num_H`, `remove_terminal_HDT` | 0 |
| `Canonicalization_step` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:1612` | external | reachable | `FillOutINChIReducedWarn`, `DeAllocBCN`, `GetBaseCanonRanking`, `DeAllocateCS`, `AllocateCS`, `Canon_INChI`, `CheckCanonNumberingCorrectness`, `CreateNeighList`, `FreeNeighList`, `free_t_group_info` | 0 |
| `CreateCompositeNormAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:1973` | external | reachable | `CreateCompAtomData` | 0 |
| `CreateCompAtomData` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:2232` | external | reachable | `CreateInpAtom`, `FreeCompAtomData` | 0 |
| `FillOutINChIReducedWarn` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:2254` | external | reachable | `AddErrorMessage`, `AllocateAndFillHillFormula`, `CopyLinearCTStereoToINChIStereo`, `MarkAmbiguousStereo`, `UnmarkAllUndefinedUnknownStereo`, `switch_ptrs`, `inchi_qsort`, `SortTautomerGroupsAndEndpoints`, `get_unusual_el_valence` | 0 |
| `make_norm_atoms_from_inp_atoms` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:2937` | external | reachable | none | 0 |
| `FreeInchi_Atom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:106` | external | reachable | none | 0 |
| `CreateInchiAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:117` | external | reachable | none | 0 |
| `MakeINCHIFromMolfileText` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:127` | external | export root | `produce_generation_output`, `copy_corrected_log_tail`, `PrepareToMakeINCHI`, `PostMakeINCHICleanup`, `AddErrorMessage`, `ProcessOneStructureEx`, `GetOneStructure` | 0 |
| `PrepareToMakeINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:285` | internal | reachable | `parse_options_string`, `inchi_ios_init`, `inchi_ios_eprint`, `inchi_strbuf_init`, `ReadCommandLineParms` | 0 |
| `PostMakeINCHICleanup` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:391` | internal | reachable | `inchi_ios_close`, `inchi_strbuf_close`, `SetBitFree`, `FreeOrigAtData`, `FreeAllINChIArrays` | 0 |
| `FreeInchi_Input` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:439` | external | reachable | `FreeInchi_Atom`, `FreeInchi_Stereo0D` | 0 |
| `is_in_the_slist` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:531` | external | reachable | none | 0 |
| `is_element_a_metal` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:543` | external | reachable | none | 0 |
| `InchiToInchiAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:582` | external | reachable | `CreateInchiAtom`, `is_in_the_slist`, `is_element_a_metal`, `inchi_ios_gets`, `inchi_ios_getsTab`, `inchi_ios_getsTab1`, `AddErrorMessage`, `CreateInchi_Stereo0D`, `FreeInchi_Stereo0D`, `FindToken`, `LoadLine`, `mystrncpy`, `lrtrim` | 0 |
| `GetSingleStereoCode` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:145` | internal | reachable | `STATUS_PushMessage` | 0 |
| `GetDoubleStereoCode` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:181` | internal | reachable | `STATUS_PushMessage` | 0 |
| `BUILDER_Unpack` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:198` | internal | reachable | `STATUS_PushMessage` | 0 |
| `BUILDER_Pack` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:212` | internal | reachable | none | 0 |
| `TranslateTetrahedralVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:219` | internal | reachable | `IXA_MOL_GetAtomIndex`, `IXA_MOL_GetStereoVertex`, `IXA_STATUS_HasError` | 0 |
| `ExtendAllene` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:246` | internal | reachable | `MOL_GetBondOtherAtom`, `IXA_MOL_GetAtomNumBonds`, `IXA_MOL_GetAtomBond`, `IXA_MOL_GetBondType`, `STATUS_PushMessage`, `IXA_STATUS_HasError` | 0 |
| `ExtendCumulene` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:313` | internal | reachable | `MOL_GetBondOtherAtom`, `IXA_MOL_GetAtomNumBonds`, `IXA_MOL_GetAtomBond`, `IXA_MOL_GetBondType`, `STATUS_PushMessage`, `IXA_STATUS_HasError` | 0 |
| `IsRectangularVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:361` | internal | reachable | `IXA_MOL_GetCommonBond`, `IXA_STATUS_HasError` | 0 |
| `IsRectOrAntiRectCentre` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:392` | internal | reachable | `IsRectangularVertex` | 0 |
| `ClearMolecule` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:422` | internal | reachable | `FreeInChIExtInput` | 0 |
| `AppendOption` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:441` | internal | reachable | none | 0 |
| `BUILDER_ClearOptions` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:453` | internal | reachable | none | 0 |
| `BUILDER_Update` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:513` | internal | reachable | `FreeINCHI`, `GetINCHIEx`, `AppendOption`, `STATUS_PushMessage` | 0 |
| `IXA_INCHIBUILDER_Create` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:802` | external | export root | `BUILDER_Pack`, `BUILDER_ClearOptions`, `STATUS_PushMessage` | 0 |
| `IXA_INCHIBUILDER_Destroy` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:823` | external | export root | `BUILDER_Unpack`, `ClearMolecule` | 0 |
| `IXA_INCHIBUILDER_SetMolecule` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:841` | external | export root | `GetSingleStereoCode`, `GetDoubleStereoCode`, `BUILDER_Unpack`, `TranslateTetrahedralVertex`, `ExtendAllene`, `ExtendCumulene`, `IsRectOrAntiRectCentre`, `ClearMolecule`, `MOL_Unpack`, `IXA_MOL_GetChiral`, `IXA_MOL_GetNumAtoms`, `IXA_MOL_GetAtomId`, `IXA_MOL_GetAtomIndex`, `IXA_MOL_GetAtomNumBonds`, `IXA_MOL_GetAtomBond`, `IXA_MOL_GetAtomX`, `IXA_MOL_GetAtomY`, `IXA_MOL_GetAtomZ`, `IXA_MOL_GetAtomElement`, `IXA_MOL_GetAtomHydrogens`, `IXA_MOL_GetAtomMass`, `IXA_MOL_GetAtomRadical`, `IXA_MOL_GetAtomCharge`, `IXA_MOL_GetBondAtom1`, `IXA_MOL_GetBondAtom2`, `IXA_MOL_GetBondType`, `IXA_MOL_GetBondWedge`, `IXA_MOL_GetDblBondConfig`, `IXA_MOL_GetNumStereos`, `IXA_MOL_GetStereoId`, `IXA_MOL_GetStereoTopology`, `IXA_MOL_GetStereoCentralAtom`, `IXA_MOL_GetStereoCentralBond`, `IXA_MOL_GetStereoVertex`, `IXA_MOL_GetStereoParity`, `STATUS_PushMessage`, `IXA_STATUS_HasError` | 0 |
| `IXA_INCHIBUILDER_SetOption_Stereo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1583` | external | export root | `BUILDER_Unpack` | 0 |
| `IXA_INCHIBUILDER_SetOption_Timeout` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1597` | external | export root | `BUILDER_Unpack`, `IXA_INCHIBUILDER_SetOption_Timeout_MilliSeconds` | 0 |
| `IXA_INCHIBUILDER_SetOption_Timeout_MilliSeconds` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1612` | external | export root | `BUILDER_Unpack` | 0 |
| `IXA_INCHIBUILDER_SetOption` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1626` | external | export root | `BUILDER_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_INCHIBUILDER_CheckOption` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1776` | external | export root | `BUILDER_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_INCHIBUILDER_CheckOption_Stereo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1971` | external | export root | `BUILDER_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_INCHIBUILDER_GetOption_Timeout_MilliSeconds` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1993` | external | export root | `BUILDER_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_INCHIBUILDER_GetInChIVersion` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:2008` | external | export root | none | 0 |
| `IXA_INCHIBUILDER_GetInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:2022` | external | export root | `BUILDER_Unpack`, `BUILDER_Update` | 0 |
| `IXA_INCHIBUILDER_GetInChIEx` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:2037` | external | export root | `BUILDER_Unpack`, `BUILDER_Update` | 0 |
| `IXA_INCHIBUILDER_GetAuxInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:2052` | external | export root | `BUILDER_Unpack`, `BUILDER_Update` | 0 |
| `IXA_INCHIBUILDER_GetLog` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:2067` | external | export root | `BUILDER_Unpack`, `BUILDER_Update` | 0 |
| `KEYBUILDER_Unpack` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_inchikey_builder.c:59` | internal | reachable | `STATUS_PushMessage` | 0 |
| `KEYBUILDER_Pack` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_inchikey_builder.c:73` | internal | reachable | none | 0 |
| `IXA_INCHIKEYBUILDER_Create` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_inchikey_builder.c:80` | external | export root | `KEYBUILDER_Pack`, `STATUS_PushMessage` | 0 |
| `IXA_INCHIKEYBUILDER_Destroy` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_inchikey_builder.c:95` | external | export root | `KEYBUILDER_Unpack` | 0 |
| `IXA_INCHIKEYBUILDER_SetInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_inchikey_builder.c:107` | external | export root | `KEYBUILDER_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_INCHIKEYBUILDER_GetInChIKey` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_inchikey_builder.c:131` | external | export root | `KEYBUILDER_Unpack`, `GetINCHIKeyFromINCHI` | 0 |
| `GetVertexCount` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:118` | internal | reachable | `STATUS_PushMessage` | 0 |
| `MOL_PackAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:135` | internal | reachable | none | 0 |
| `MOL_PackBond` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:142` | internal | reachable | none | 0 |
| `MOL_PackStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:149` | internal | reachable | none | 0 |
| `MOL_PackPolymerUnit` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:156` | internal | reachable | none | 0 |
| `MOL_UnpackAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:163` | internal | reachable | `STATUS_PushMessage` | 0 |
| `MOL_UnpackBond` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:180` | internal | reachable | `STATUS_PushMessage` | 0 |
| `MOL_UnpackStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:197` | internal | reachable | `STATUS_PushMessage` | 0 |
| `MOL_UnpackPolymerUnit` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:214` | internal | reachable | `STATUS_PushMessage` | 0 |
| `MOL_GetAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:238` | internal | reachable | `MOL_UnpackAtom` | 0 |
| `MOL_GetBond` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:249` | internal | reachable | `MOL_UnpackBond` | 0 |
| `MOL_GetStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:260` | internal | reachable | `MOL_UnpackStereo` | 0 |
| `MOL_GetSGroup` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:271` | internal | reachable | `MOL_UnpackPolymerUnit` | 0 |
| `MOL_ClearExtMolData` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:285` | external | reachable | none | 0 |
| `MOL_Clear` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:395` | external | reachable | `MOL_ClearExtMolData` | 0 |
| `MOL_GuessNewSize` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:408` | internal | reachable | none | 0 |
| `MOL_CreateAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:428` | internal | reachable | `MOL_PackAtom`, `MOL_GuessNewSize`, `STATUS_PushMessage` | 0 |
| `MOL_CreateStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:520` | internal | reachable | `MOL_GuessNewSize`, `STATUS_PushMessage` | 0 |
| `MOL_CreatePolymerUnit` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:598` | internal | reachable | `MOL_PackPolymerUnit`, `MOL_GuessNewSize`, `STATUS_PushMessage` | 0 |
| `MOL_Pack` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:671` | internal | reachable | none | 0 |
| `MOL_Unpack` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:678` | external | reachable | `STATUS_PushMessage` | 0 |
| `MOL_GetBondOtherAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:692` | external | reachable | `IXA_MOL_GetBondAtom1`, `IXA_MOL_GetBondAtom2`, `STATUS_PushMessage`, `IXA_STATUS_HasError` | 0 |
| `IXA_MOL_Create` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:729` | external | export root | `MOL_Pack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_Destroy` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:744` | external | export root | `MOL_Clear`, `MOL_Unpack` | 0 |
| `IXA_MOL_Clear` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:759` | external | export root | `MOL_Clear`, `MOL_Unpack` | 0 |
| `IXA_MOL_SetChiral` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:773` | external | export root | `MOL_Unpack` | 0 |
| `IXA_MOL_GetChiral` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:788` | external | export root | `MOL_Unpack` | 0 |
| `IXA_MOL_CreateAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:801` | external | export root | `MOL_CreateAtom`, `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_GetNumAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:816` | external | export root | `MOL_Unpack` | 0 |
| `IXA_MOL_GetAtomId` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:830` | external | export root | `MOL_PackAtom`, `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_GetAtomIndex` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:857` | external | export root | `MOL_UnpackAtom`, `MOL_Unpack` | 0 |
| `IXA_MOL_GetAtomNumBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:879` | external | export root | `MOL_GetAtom`, `MOL_Unpack` | 0 |
| `IXA_MOL_GetAtomBond` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:902` | external | export root | `MOL_GetAtom`, `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_SetAtomX` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:929` | external | export root | `MOL_GetAtom`, `MOL_Unpack` | 0 |
| `IXA_MOL_GetAtomX` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:953` | external | export root | `MOL_GetAtom`, `MOL_Unpack` | 0 |
| `IXA_MOL_SetAtomY` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:976` | external | export root | `MOL_GetAtom`, `MOL_Unpack` | 0 |
| `IXA_MOL_GetAtomY` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1000` | external | export root | `MOL_GetAtom`, `MOL_Unpack` | 0 |
| `IXA_MOL_SetAtomZ` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1023` | external | export root | `MOL_GetAtom`, `MOL_Unpack` | 0 |
| `IXA_MOL_GetAtomZ` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1047` | external | export root | `MOL_GetAtom`, `MOL_Unpack` | 0 |
| `IXA_MOL_SetAtomElement` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1070` | external | export root | `MOL_GetAtom`, `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_GetAtomElement` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1110` | external | export root | `MOL_GetAtom`, `MOL_Unpack` | 0 |
| `IXA_MOL_SetAtomAtomicNumber` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1133` | external | export root | `MOL_GetAtom`, `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_GetAtomAtomicNumber` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1163` | external | export root | `MOL_GetAtom`, `MOL_Unpack` | 0 |
| `IXA_MOL_SetAtomHydrogens` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1186` | external | export root | `MOL_GetAtom`, `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_GetAtomHydrogens` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1236` | external | export root | `MOL_GetAtom`, `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_SetAtomMass` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1271` | external | export root | `MOL_GetAtom`, `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_GetAtomMass` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1315` | external | export root | `MOL_GetAtom`, `MOL_Unpack` | 0 |
| `IXA_MOL_SetAtomRadical` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1338` | external | export root | `MOL_GetAtom`, `MOL_Unpack` | 0 |
| `IXA_MOL_GetAtomRadical` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1362` | external | export root | `MOL_GetAtom`, `MOL_Unpack` | 0 |
| `IXA_MOL_SetAtomCharge` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1385` | external | export root | `MOL_GetAtom`, `MOL_Unpack` | 0 |
| `IXA_MOL_GetAtomCharge` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1409` | external | export root | `MOL_GetAtom`, `MOL_Unpack` | 0 |
| `IXA_MOL_ReserveSpace` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1432` | external | export root | `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_CreateBond` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1497` | external | export root | `MOL_PackBond`, `MOL_GetAtom`, `MOL_GetBond`, `MOL_GuessNewSize`, `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_GetNumBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1643` | external | export root | `MOL_Unpack` | 0 |
| `IXA_MOL_GetBondId` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1657` | external | export root | `MOL_PackBond`, `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_GetBondIndex` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1684` | external | export root | `MOL_UnpackBond`, `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_GetBondAtom1` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1707` | external | export root | `MOL_GetBond`, `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_GetBondAtom2` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1732` | external | export root | `MOL_GetBond`, `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_GetBondOtherAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1756` | external | export root | `MOL_GetBond`, `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_SetBondType` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1796` | external | export root | `MOL_GetBond`, `MOL_Unpack` | 0 |
| `IXA_MOL_GetBondType` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1822` | external | export root | `MOL_GetBond`, `MOL_Unpack` | 0 |
| `IXA_MOL_SetBondWedge` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1844` | external | export root | `MOL_GetBond`, `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_GetBondWedge` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1881` | external | export root | `MOL_GetBond`, `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_SetDblBondConfig` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1917` | external | export root | `MOL_GetBond`, `MOL_Unpack` | 0 |
| `IXA_MOL_GetDblBondConfig` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1941` | external | export root | `MOL_GetBond`, `MOL_Unpack` | 0 |
| `IXA_MOL_GetCommonBond` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1964` | external | export root | `MOL_GetAtom`, `MOL_Unpack` | 0 |
| `IXA_MOL_CreateStereoTetrahedron` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2010` | external | export root | `MOL_PackStereo`, `MOL_CreateStereo`, `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_CreateStereoRectangle` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2045` | external | export root | `MOL_PackStereo`, `MOL_CreateStereo`, `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_CreateStereoAntiRectangle` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2079` | external | export root | `MOL_PackStereo`, `MOL_CreateStereo`, `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_GetNumStereos` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2113` | external | export root | `MOL_Unpack` | 0 |
| `IXA_MOL_GetStereoId` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2127` | external | export root | `MOL_PackStereo`, `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_GetStereoIndex` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2148` | external | export root | `MOL_UnpackStereo`, `MOL_Unpack` | 0 |
| `IXA_MOL_GetStereoTopology` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2164` | external | export root | `MOL_GetStereo`, `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_GetStereoCentralAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2186` | external | export root | `MOL_GetStereo`, `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_GetStereoCentralBond` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2219` | external | export root | `MOL_GetStereo`, `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_GetStereoNumVertices` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2252` | external | export root | `GetVertexCount`, `MOL_GetStereo`, `MOL_Unpack` | 0 |
| `IXA_MOL_GetStereoVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2274` | external | export root | `GetVertexCount`, `MOL_GetStereo`, `MOL_Unpack`, `STATUS_PushMessage`, `IXA_STATUS_HasError` | 0 |
| `IXA_MOL_SetStereoParity` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2310` | external | export root | `MOL_GetStereo`, `MOL_Unpack` | 0 |
| `IXA_MOL_GetStereoParity` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2331` | external | export root | `MOL_GetStereo`, `MOL_Unpack` | 0 |
| `IXA_MOL_SetPolymerUnit` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2353` | external | export root | `MOL_GetSGroup`, `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_SetExtMolDataByInChIExtInput` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2429` | external | reachable | `MOL_ClearExtMolData`, `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_CreatePolymerUnit` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2681` | external | export root | `MOL_CreatePolymerUnit`, `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_GetPolymerUnitId` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2696` | external | export root | `MOL_PackPolymerUnit`, `MOL_Unpack`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_GetPolymerUnitIndex` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2717` | external | export root | `MOL_UnpackPolymerUnit`, `MOL_Unpack` | 0 |
| `AnalyseInternalVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_read_inchi.c:55` | internal | reachable | `MOL_GetBondOtherAtom`, `IXA_MOL_GetAtomNumBonds`, `IXA_MOL_GetAtomBond`, `STATUS_PushMessage`, `IXA_STATUS_HasError` | 0 |
| `FindCumuleneCentre` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_read_inchi.c:112` | internal | reachable | `MOL_GetBondOtherAtom`, `IXA_MOL_GetAtomNumBonds`, `IXA_MOL_GetAtomBond`, `IXA_MOL_GetBondType`, `IXA_MOL_GetCommonBond`, `STATUS_PushMessage`, `IXA_STATUS_HasError` | 0 |
| `IXA_MOL_ReadInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_read_inchi.c:189` | external | export root | `GetStructFromINCHIEx`, `FreeStructFromINCHIEx`, `IXA_MOL_Clear`, `IXA_MOL_CreateAtom`, `IXA_MOL_GetAtomId`, `IXA_MOL_GetAtomNumBonds`, `IXA_MOL_SetAtomX`, `IXA_MOL_SetAtomY`, `IXA_MOL_SetAtomZ`, `IXA_MOL_SetAtomElement`, `IXA_MOL_SetAtomHydrogens`, `IXA_MOL_SetAtomMass`, `IXA_MOL_SetAtomRadical`, `IXA_MOL_SetAtomCharge`, `IXA_MOL_CreateBond`, `IXA_MOL_GetBondAtom1`, `IXA_MOL_GetBondAtom2`, `IXA_MOL_SetBondType`, `IXA_MOL_SetBondWedge`, `IXA_MOL_SetDblBondConfig`, `IXA_MOL_GetCommonBond`, `IXA_MOL_CreateStereoTetrahedron`, `IXA_MOL_CreateStereoRectangle`, `IXA_MOL_CreateStereoAntiRectangle`, `IXA_MOL_SetStereoParity`, `IXA_MOL_SetExtMolDataByInChIExtInput`, `AnalyseInternalVertex`, `FindCumuleneCentre`, `STATUS_PushMessage`, `IXA_STATUS_HasError`, `get_atomic_mass` | 0 |
| `IXA_MOL_ReadAuxInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_read_inchi.c:686` | external | export root | `Get_inchi_Input_FromAuxInfo`, `Free_inchi_Input`, `FreeINCHI`, `GetINCHI`, `IXA_MOL_ReadInChI`, `STATUS_PushMessage` | 0 |
| `IXA_MOL_ReadMolfile` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_read_mol.c:78` | external | export root | `IXA_MOL_Clear`, `IXA_MOL_SetChiral`, `IXA_MOL_CreateAtom`, `IXA_MOL_GetAtomId`, `IXA_MOL_GetAtomNumBonds`, `IXA_MOL_GetAtomBond`, `IXA_MOL_SetAtomX`, `IXA_MOL_SetAtomY`, `IXA_MOL_SetAtomZ`, `IXA_MOL_SetAtomElement`, `IXA_MOL_SetAtomHydrogens`, `IXA_MOL_SetAtomMass`, `IXA_MOL_SetAtomRadical`, `IXA_MOL_SetAtomCharge`, `IXA_MOL_ReserveSpace`, `IXA_MOL_CreateBond`, `IXA_MOL_SetBondType`, `IXA_MOL_GetBondType`, `IXA_MOL_SetBondWedge`, `IXA_MOL_SetDblBondConfig`, `IXA_MOL_SetPolymerUnit`, `IXA_MOL_CreatePolymerUnit`, `IXA_MOL_GetPolymerUnitId`, `STATUS_PushMessage`, `IXA_STATUS_HasError`, `inchi_ios_init`, `ReadMolfile`, `FreeMolfileData` | 0 |
| `BLOCK_clear` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:81` | internal | reachable | `BLOCK_clear` | 0 |
| `STATUS_init` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:93` | internal | reachable | none | 0 |
| `STATUS_Clear` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:102` | internal | reachable | `BLOCK_clear`, `STATUS_init` | 0 |
| `STATUS_Pack` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:123` | external | reachable | none | 0 |
| `STATUS_Unpack` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:131` | internal | reachable | none | 0 |
| `STATUS_PushMessage` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:138` | external | reachable | `STATUS_Unpack` | 0 |
| `INCHISTATUS_TestSeverity` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:211` | internal | reachable | `STATUS_Unpack` | 0 |
| `IXA_STATUS_Create` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:233` | external | export root | `STATUS_init`, `STATUS_Pack` | 0 |
| `IXA_STATUS_Destroy` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:246` | external | export root | `STATUS_Clear`, `STATUS_Unpack` | 0 |
| `IXA_STATUS_Clear` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:259` | external | export root | `STATUS_Clear`, `STATUS_Unpack` | 0 |
| `IXA_STATUS_HasError` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:271` | external | export root | `INCHISTATUS_TestSeverity` | 0 |
| `IXA_STATUS_HasWarning` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:278` | external | export root | `INCHISTATUS_TestSeverity` | 0 |
| `IXA_STATUS_GetCount` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:285` | external | export root | `STATUS_Unpack` | 0 |
| `IXA_STATUS_GetMessage` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:298` | external | export root | `STATUS_Unpack` | 0 |
| `IXA_STATUS_GetSeverity` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:316` | external | export root | `STATUS_Unpack` | 0 |
| `max_3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/bcf_s.c:58` | external | active, export-unreachable | none | 0 |
| `memcpy_custom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/bcf_s.c:68` | external | active, export-unreachable | none | 0 |
| `dbl2int` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/bcf_s.c:92` | external | active, export-unreachable | `dbl2int_f`, `dbl2int_e`, `dbl2int_g` | 0 |
| `dbl2int_f` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/bcf_s.c:117` | internal | active, export-unreachable | none | 0 |
| `dbl2int_e` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/bcf_s.c:188` | internal | active, export-unreachable | none | 0 |
| `dbl2int_g` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/bcf_s.c:276` | internal | active, export-unreachable | none | 0 |
| `RestoreEdgeFlow` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:442` | external | reachable | none | 0 |
| `SetAtomBondType` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:480` | external | reachable | none | 0 |
| `RunBalancedNetworkSearch` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:689` | external | reachable | `ReInitBnData`, `BalancedNetworkSearch`, `bInchiTimeIsOver` | 0 |
| `SetAtomRadAndChemValFromVertexCapFlow` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:726` | external | reachable | none | 0 |
| `AddChangedAtHChargeBNS` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:775` | external | reachable | `GetAtomChargeType` | 0 |
| `EliminatePlusMinusChargeAmbiguity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:812` | external | reachable | none | 0 |
| `AddOrRemoveExplOrImplH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:902` | external | reachable | `get_opposite_sb_atom` | 0 |
| `SubtractOrChangeAtHChargeBNS` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:1057` | external | reachable | `AddOrRemoveExplOrImplH`, `GetAtomChargeType` | 0 |
| `SetBondsFromBnStructFlow` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:1165` | external | reachable | `SetAtomBondType`, `SetAtomRadAndChemValFromVertexCapFlow` | 0 |
| `MarkAtomsAtTautGroups` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:1331` | external | reachable | none | 0 |
| `RestoreBnStructFlow` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:1525` | external | reachable | `RestoreEdgeFlow` | 0 |
| `bNeedToTestTheFlow` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:1590` | external | reachable | none | 0 |
| `nBondsValenceInpAt` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:1673` | external | reachable | none | 0 |
| `BnsAdjustFlowBondsRad` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:1725` | external | reachable | `RunBalancedNetworkSearch`, `SetBondsFromBnStructFlow`, `RestoreBnStructFlow`, `nBondsValenceInpAt`, `ReInitBnStructAltPaths` | 0 |
| `BnsTestAndMarkAltBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:1825` | external | reachable | `RunBalancedNetworkSearch`, `SetBondsFromBnStructFlow`, `RestoreBnStructFlow`, `bNeedToTestTheFlow`, `nMaxFlow2Check`, `nCurFlow2Check`, `nMinFlow2Check`, `bSetBondsAfterCheckOneBond`, `bRestoreFlowAfterCheckOneBond`, `bSetFlowToCheckOneBond`, `ReInitBnStructAltPaths` | 0 |
| `remove_alt_bond_marks` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:1987` | internal | reachable | none | 0 |
| `SetForbiddenEdges` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:2002` | external | reachable | `fix_special_bonds` | 0 |
| `fix_special_bonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:2181` | external | reachable | `get_el_valence`, `is_el_a_metal`, `num_of_H`, `ion_el_group`, `nNoMetalNumBonds`, `nNoMetalBondsValence`, `nNoMetalNeighIndex`, `nNoMetalOtherNeighIndex`, `nNoMetalOtherNeighIndex2` | 0 |
| `fix_explicitly_indicated_bonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:2704` | external | active, export-unreachable | none | 0 |
| `is_Z_atom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:2781` | external | reachable | none | 0 |
| `IsZOX` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:2814` | external | reachable | none | 0 |
| `update_some_attype_totals` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:2846` | external | active, export-unreachable | none | 0 |
| `GetAtomChargeType` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:2873` | external | reachable | `is_Z_atom`, `IsZOX`, `detect_unusual_el_valence`, `is_el_a_metal` | 0 |
| `SimpleRemoveHplusNPO` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:3432` | external | reachable | `AddOrRemoveExplOrImplH`, `GetAtomChargeType` | 0 |
| `bIsAtomTypeHard` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:3480` | external | reachable | `GetAtomChargeType` | 0 |
| `bIsHDonorAccAtomType` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:3501` | external | reachable | `bIsAtomTypeHard` | 0 |
| `bIsNegAtomType` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:3539` | external | reachable | `bIsAtomTypeHard` | 0 |
| `bIsHardRemHCandidate` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:3582` | external | reachable | `bIsHDonorAccAtomType`, `bIsNegAtomType` | 0 |
| `CreateCGroupInBnStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:3599` | external | reachable | `GetAtomChargeType` | 0 |
| `CreateTGroupInBnStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:3813` | external | reachable | `GetAtomChargeType` | 0 |
| `RemoveLastGroupFromBnStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:3988` | external | reachable | none | 0 |
| `SetInitCapFlowToCurrent` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:4094` | external | reachable | none | 0 |
| `SimpleRemoveAcidicProtons` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:4129` | external | reachable | `AddOrRemoveExplOrImplH`, `GetAtomChargeType` | 0 |
| `bHasAcidicHydrogen` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:4200` | external | reachable | `GetAtomChargeType` | 0 |
| `bHasOtherExchangableH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:4220` | external | reachable | `GetAtomChargeType` | 0 |
| `SimpleAddAcidicProtons` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:4253` | external | reachable | `AddOrRemoveExplOrImplH`, `GetAtomChargeType` | 0 |
| `bHasAcidicMinus` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:4324` | external | reachable | `GetAtomChargeType` | 0 |
| `HardRemoveAcidicProtons` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:4351` | external | reachable | `CreateCGroupInBnStruct`, `CreateTGroupInBnStruct`, `RemoveLastGroupFromBnStruct`, `bExistsAltPath` | 0 |
| `HardAddAcidicProtons` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:4556` | external | reachable | `CreateCGroupInBnStruct`, `CreateTGroupInBnStruct`, `RemoveLastGroupFromBnStruct`, `bExistsAltPath` | 0 |
| `HardRemoveHplusNP` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:4770` | external | reachable | `CreateCGroupInBnStruct`, `CreateTGroupInBnStruct`, `RemoveLastGroupFromBnStruct`, `bExistsAltPath` | 0 |
| `mark_at_type` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:5040` | external | reachable | `GetAtomChargeType` | 0 |
| `RemoveNPProtonsAndAcidCharges` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:5084` | external | reachable | `SimpleRemoveHplusNPO`, `SimpleRemoveAcidicProtons`, `SimpleAddAcidicProtons`, `HardRemoveAcidicProtons`, `HardAddAcidicProtons`, `HardRemoveHplusNP` | 0 |
| `mark_alt_bonds_and_taut_groups` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:5242` | external | reachable | `BnsAdjustFlowBondsRad`, `BnsTestAndMarkAltBonds`, `remove_alt_bond_marks`, `SetForbiddenEdges`, `SetInitCapFlowToCurrent`, `mark_at_type`, `RemoveNPProtonsAndAcidCharges`, `AllocateAndInitBnStruct`, `DeAllocateBnStruct`, `ReInitBnStructAddGroups`, `AddCGroups2BnStruct`, `DeAllocateBnData`, `AllocateAndInitBnData`, `MarkRingSystemsAltBns`, `ReInitBnStructForAltBns`, `MarkNonStereoAltBns`, `MarkChargeGroups`, `MarkSaltChargeGroups2`, `MarkSaltChargeGroups`, `MergeSaltTautGroups`, `MakeIsotopicHGroup`, `MarkTautomerGroups` | 0 |
| `nMaxFlow2Check` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:6010` | external | reachable | none | 0 |
| `nCurFlow2Check` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:6024` | external | reachable | none | 0 |
| `nMinFlow2Check` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:6033` | external | reachable | none | 0 |
| `bSetBondsAfterCheckOneBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:6071` | external | reachable | `SetAtomBondType` | 0 |
| `bRestoreFlowAfterCheckOneBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:6180` | external | reachable | none | 0 |
| `bSetFlowToCheckOneBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:6222` | external | reachable | none | 0 |
| `bAddNewVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:6482` | external | reachable | none | 0 |
| `AddNewEdge` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:6556` | external | reachable | none | 0 |
| `GetEdgeToGroupVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:6610` | external | reachable | none | 0 |
| `GetGroupVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:6644` | external | reachable | none | 0 |
| `bAddStCapToAVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:6703` | external | reachable | none | 0 |
| `bSetBnsToCheckAltPath` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:6778` | external | reachable | `bSetFlowToCheckOneBond`, `bAddNewVertex`, `GetEdgeToGroupVertex`, `GetGroupVertex`, `bAddStCapToAVertex` | 0 |
| `bRestoreBnsAfterCheckAltPath` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:7398` | external | reachable | none | 0 |
| `bExistsAnyAltPath` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:7518` | external | reachable | `bExistsAltPath` | 0 |
| `bIsBnsEndpoint` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:7561` | external | reachable | none | 0 |
| `bRadChangesAtomType` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:7588` | external | reachable | `GetPrevVertex` | 0 |
| `RegisterRadEndpoint` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:7639` | external | reachable | `bRadChangesAtomType`, `GetPrevVertex` | 0 |
| `cmp_rad_endpoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:7784` | external | active, export-unreachable | none | 0 |
| `RemoveRadEndpoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:7812` | external | reachable | none | 0 |
| `RestoreRadicalsOnly` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:7907` | external | reachable | none | 0 |
| `SetRadEndpoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:7969` | external | reachable | `bAddNewVertex`, `AddNewEdge`, `RemoveRadEndpoints`, `ReInitBnStructAltPaths`, `ReInitBnData`, `BalancedNetworkSearch` | 0 |
| `SetRadEndpoints2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:8064` | external | reachable | `bAddNewVertex`, `AddNewEdge`, `RemoveRadEndpoints`, `ReInitBnStructAltPaths`, `ReInitBnData`, `BalancedNetworkSearch`, `NodeSetCreate`, `NodeSetFree`, `NodeSetFromRadEndpoints`, `RemoveFromNodeSet`, `DoNodeSetsIntersect`, `IsNodeSetEmpty`, `AddNodeSet2ToNodeSet1`, `AddNodesToRadEndpoints`, `SetBitCreate` | 0 |
| `bExistsAltPath` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:8393` | external | reachable | `RunBalancedNetworkSearch`, `AddChangedAtHChargeBNS`, `EliminatePlusMinusChargeAmbiguity`, `SubtractOrChangeAtHChargeBNS`, `SetBondsFromBnStructFlow`, `MarkAtomsAtTautGroups`, `RestoreBnStructFlow`, `bSetBondsAfterCheckOneBond`, `bRestoreFlowAfterCheckOneBond`, `bSetBnsToCheckAltPath`, `bRestoreBnsAfterCheckAltPath`, `bIsBnsEndpoint`, `RemoveRadEndpoints`, `RestoreRadicalsOnly`, `SetRadEndpoints2`, `ReInitBnStructAltPaths`, `nGetEndpointInfo`, `nGetEndpointInfo_PT_22_00`, `nGetEndpointInfo_PT_16_00`, `nGetEndpointInfo_PT_06_00`, `nGetEndpointInfo_PT_39_00`, `nGetEndpointInfo_PT_13_00`, `nGetEndpointInfo_PT_18_00`, `nGetEndpointInfo_KET` | 0 |
| `AllocateAndInitBnStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:8783` | external | reachable | `DeAllocateBnStruct`, `is_centerpoint_elem`, `get_endpoint_valence` | 0 |
| `DeAllocateBnStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9010` | external | reachable | none | 0 |
| `ReInitBnStructAltPaths` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9042` | external | reachable | none | 0 |
| `ReInitBnStructAddGroups` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9062` | external | reachable | `ReInitBnStruct`, `AddTGroups2BnStruct`, `AddCGroups2BnStruct` | 0 |
| `ReInitBnStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9102` | external | reachable | `ReInitBnStructAltPaths` | 0 |
| `CompTGroupNumber` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9195` | external | active, export-unreachable | none | 0 |
| `CompCGroupNumber` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9202` | external | active, export-unreachable | none | 0 |
| `AddTGroups2BnStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9209` | external | reachable | `insertions_sort`, `nGetEndpointInfo`, `nGetEndpointInfo_PT_22_00`, `nGetEndpointInfo_PT_16_00`, `nGetEndpointInfo_PT_06_00`, `nGetEndpointInfo_PT_39_00`, `nGetEndpointInfo_PT_13_00`, `nGetEndpointInfo_PT_18_00`, `nGetEndpointInfo_KET` | 0 |
| `AddCGroups2BnStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9406` | external | reachable | `insertions_sort` | 0 |
| `ClearAllBnDataVertices` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9609` | external | reachable | none | 0 |
| `ClearAllBnDataEdges` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9620` | external | reachable | none | 0 |
| `DeAllocateBnData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9631` | external | reachable | none | 0 |
| `AllocateAndInitBnData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9677` | external | reachable | `ClearAllBnDataVertices`, `ClearAllBnDataEdges`, `DeAllocateBnData` | 0 |
| `ReInitBnData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9718` | external | reachable | none | 0 |
| `GetVertexDegree` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9774` | external | reachable | none | 0 |
| `Get2ndEdgeVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9796` | external | reachable | none | 0 |
| `GetVertexNeighbor` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9834` | external | reachable | none | 0 |
| `GetEdgePointer` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9879` | external | reachable | none | 0 |
| `AugmentEdge` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9934` | external | reachable | `GetEdgePointer` | 0 |
| `rescap_mark` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:10102` | external | reachable | `GetEdgePointer` | 0 |
| `GetPrevVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:10163` | external | reachable | `Get2ndEdgeVertex` | 0 |
| `bIgnoreVertexNonTACN_atom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:10291` | external | reachable | `GetVertexDegree`, `GetVertexNeighbor`, `rescap` | 0 |
| `bIgnoreVertexNonTACN_group` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:10419` | external | reachable | `GetPrevVertex` | 0 |
| `rescap` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:10644` | external | reachable | `GetEdgePointer` | 0 |
| `BalancedNetworkSearch` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:10779` | external | reachable | `RegisterRadEndpoint`, `GetVertexDegree`, `Get2ndEdgeVertex`, `GetVertexNeighbor`, `bIgnoreVertexNonTACN_atom`, `bIgnoreVertexNonTACN_group`, `rescap`, `FindBase`, `MakeBlossom`, `PullFlow`, `FindPathCap` | 0 |
| `FindBase` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:11068` | external | reachable | `FindBase` | 0 |
| `FindPathToVertex_s` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:11092` | external | reachable | `FindBase` | 0 |
| `MakeBlossom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:11122` | external | reachable | `rescap`, `FindPathToVertex_s` | 0 |
| `PullFlow` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:11290` | external | reachable | `Get2ndEdgeVertex`, `AugmentEdge`, `PullFlow` | 0 |
| `FindPathCap` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:11367` | external | reachable | `Get2ndEdgeVertex`, `rescap_mark`, `FindPathCap` | 0 |
| `MarkRingSystemsAltBns` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:11459` | external | reachable | none | 0 |
| `ReInitBnStructForAltBns` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:11800` | external | reachable | `ReInitBnStruct` | 0 |
| `MarkNonStereoAltBns` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:11909` | external | reachable | none | 0 |
| `bHasChargedNeighbor` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:12005` | external | reachable | none | 0 |
| `AddRemoveProtonsRestr` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:12056` | external | reachable | `GetAtomChargeType`, `mark_at_type`, `bHasChargedNeighbor`, `is_centerpoint_elem`, `get_endpoint_valence` | 0 |
| `AddRemoveIsoProtonsRestr` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:12360` | external | reachable | `bHeteroAtomMayHaveXchgIsoH` | 0 |
| `inchi_ios_init` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:85` | external | reachable | none | 0 |
| `inchi_ios_create_copy` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:106` | external | active, export-unreachable | none | 0 |
| `inchi_ios_flush` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:144` | external | active, export-unreachable | none | 0 |
| `inchi_ios_flush2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:183` | external | active, export-unreachable | none | 0 |
| `inchi_ios_close` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:229` | external | reachable | none | 0 |
| `inchi_ios_reset` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:255` | external | reachable | none | 0 |
| `inchi_ios_free_str` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:272` | external | active, export-unreachable | none | 0 |
| `inchi_ios_str_getc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:294` | external | reachable | none | 0 |
| `inchi_ios_str_gets` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:324` | external | reachable | `inchi_ios_str_getc` | 0 |
| `inchi_ios_str_getsTab` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:354` | external | reachable | `inchi_ios_str_getc` | 0 |
| `inchi_ios_gets` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:386` | external | reachable | `inchi_ios_str_gets`, `lrtrim` | 0 |
| `inchi_ios_getsTab` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:420` | external | reachable | `inchi_ios_str_getsTab`, `lrtrim` | 0 |
| `inchi_ios_getsTab1` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:451` | external | reachable | `inchi_ios_str_getsTab`, `lrtrim` | 0 |
| `inchi_ios_print` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:477` | external | reachable | `GetMaxPrintfLength` | 0 |
| `push_to_winchi_text_window` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:576` | external | active, export-unreachable | none | 0 |
| `inchi_ios_print_nodisplay` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:604` | external | reachable | `inchi_print_nodisplay`, `GetMaxPrintfLength` | 0 |
| `inchi_ios_flush_not_displayed` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:679` | external | active, export-unreachable | `inchi_ios_print` | 0 |
| `inchi_ios_eprint` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:708` | external | reachable | `inchi_vfprintf`, `GetMaxPrintfLength` | 0 |
| `inchi_fprintf` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:800` | external | active, export-unreachable | `inchi_vfprintf` | 0 |
| `inchi_vfprintf` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:828` | external | reachable | none | 0 |
| `inchi_print_nodisplay` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:867` | external | reachable | none | 0 |
| `inchi_fgetsLfTab` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:893` | external | active, export-unreachable | `inchi_fgetsTab`, `lrtrim` | 0 |
| `inchi_fgetsTab` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:966` | external | active, export-unreachable | none | 0 |
| `inchi_fgetsLf` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:996` | external | reachable | `inchi_sgets` | 0 |
| `GetMaxPrintfLength` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1065` | external | reachable | none | 0 |
| `inchi_sgets` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1322` | external | reachable | none | 0 |
| `inchi_strbuf_init` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1370` | external | reachable | none | 0 |
| `inchi_strbuf_reset` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1403` | external | reachable | none | 0 |
| `inchi_strbuf_close` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1422` | external | reachable | none | 0 |
| `inchi_strbuf_create_copy` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1438` | external | reachable | none | 0 |
| `inchi_strbuf_update` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1459` | external | reachable | none | 0 |
| `inchi_strbuf_printf` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1507` | external | reachable | `GetMaxPrintfLength`, `inchi_strbuf_update` | 0 |
| `inchi_strbuf_printf_from` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1544` | external | active, export-unreachable | `GetMaxPrintfLength`, `inchi_strbuf_update` | 0 |
| `inchi_strbuf_getline` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1584` | external | active, export-unreachable | `inchi_strbuf_reset`, `inchi_strbuf_printf` | 0 |
| `inchi_strbuf_addline` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1635` | external | reachable | `inchi_ios_str_getc`, `inchi_strbuf_printf` | 0 |
| `_inchi_trace` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1699` | external | active, export-unreachable | none | 0 |
| `Output_RecordInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1709` | external | active, export-unreachable | `inchi_ios_print_nodisplay` | 0 |
| `TranspositionCreate` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:694` | external | reachable | none | 0 |
| `TranspositionFree` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:707` | external | reachable | none | 0 |
| `NodeSetCreate` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:718` | external | reachable | none | 0 |
| `NodeSetFree` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:754` | external | reachable | none | 0 |
| `CTableCreate` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:769` | external | reachable | none | 0 |
| `CTableFree` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:851` | external | reachable | none | 0 |
| `UnorderedPartitionCreate` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:904` | external | reachable | none | 0 |
| `UnorderedPartitionFree` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:919` | external | reachable | none | 0 |
| `UnorderedPartitionMakeDiscrete` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:933` | external | reachable | none | 0 |
| `PartitionCreate` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:946` | external | reachable | none | 0 |
| `PartitionFree` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:959` | external | reachable | none | 0 |
| `PartitionIsDiscrete` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:978` | external | reachable | none | 0 |
| `PartitionGetFirstCell` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:998` | external | reachable | none | 0 |
| `CellMakeEmpty` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1040` | external | reachable | none | 0 |
| `NodeSetFromVertices` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1052` | external | reachable | none | 0 |
| `AllNodesAreInSet` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1071` | external | reachable | none | 0 |
| `PartitionGetMcrAndFixSet` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1098` | external | reachable | none | 0 |
| `NodeSetFromRadEndpoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1138` | external | reachable | none | 0 |
| `RemoveFromNodeSet` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1159` | external | reachable | none | 0 |
| `DoNodeSetsIntersect` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1178` | external | reachable | none | 0 |
| `IsNodeSetEmpty` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1201` | external | reachable | none | 0 |
| `AddNodeSet2ToNodeSet1` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1223` | external | reachable | none | 0 |
| `AddNodesToRadEndpoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1241` | external | reachable | none | 0 |
| `PartitionGetTransposition` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1285` | external | reachable | none | 0 |
| `nGetMcr2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1305` | external | reachable | none | 0 |
| `nJoin2Mcrs2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1342` | external | reachable | `nGetMcr2` | 0 |
| `GetUnorderedPartitionMcrNode` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1372` | external | reachable | `nGetMcr2` | 0 |
| `UnorderedPartitionJoin` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1385` | external | reachable | `nJoin2Mcrs2` | 0 |
| `PartitionSatisfiesLemma_2_25` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1407` | external | reachable | none | 0 |
| `PartitionCopy` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1445` | external | reachable | none | 0 |
| `PartitionColorVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1464` | external | reachable | `PartitionCopy`, `DifferentiateRanks3`, `DifferentiateRanks4` | 0 |
| `CellGetMinNode` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1586` | external | reachable | none | 0 |
| `CellGetNumberOfNodes` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1748` | external | reachable | none | 0 |
| `CellIntersectWithSet` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1768` | external | reachable | none | 0 |
| `CtPartClear` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1798` | external | reachable | none | 0 |
| `insertions_sort_NeighList_AT_NUMBERS2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1819` | external | reachable | none | 0 |
| `CtPartFill` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1849` | external | reachable | `insertions_sort_NeighList_AT_NUMBERS2` | 0 |
| `CtPartINCHI_CANON_INFINITY` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:2065` | external | reachable | none | 0 |
| `CtPartCompare` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:2123` | external | reachable | none | 0 |
| `CtFullCompare` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:2623` | external | reachable | none | 0 |
| `CtFullCompareLayers` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:2863` | external | reachable | none | 0 |
| `CtCompareLayersGetFirstDiff` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:2880` | external | reachable | none | 0 |
| `CtPartCompareLayers` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:2929` | external | reachable | `CtCompareLayersGetFirstDiff` | 0 |
| `UpdateCompareLayers` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:2947` | external | reachable | none | 0 |
| `CtPartCopy` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:2965` | external | reachable | none | 0 |
| `CtFullCopy` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:3113` | external | reachable | `CtPartCopy` | 0 |
| `TranspositionGetMcrAndFixSetAndUnorderedPartition` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:3125` | external | reachable | none | 0 |
| `SetBitCreate` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:3205` | external | reachable | none | 0 |
| `SetBitFree` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:3259` | external | reachable | none | 0 |
| `GetOneAdditionalLayer` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:3566` | external | reachable | none | 0 |
| `CanonGraph` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:3662` | external | reachable | `TranspositionCreate`, `TranspositionFree`, `NodeSetCreate`, `NodeSetFree`, `CTableCreate`, `CTableFree`, `UnorderedPartitionCreate`, `UnorderedPartitionFree`, `UnorderedPartitionMakeDiscrete`, `PartitionCreate`, `PartitionFree`, `PartitionIsDiscrete`, `PartitionGetFirstCell`, `CellMakeEmpty`, `NodeSetFromVertices`, `AllNodesAreInSet`, `PartitionGetMcrAndFixSet`, `PartitionGetTransposition`, `GetUnorderedPartitionMcrNode`, `UnorderedPartitionJoin`, `PartitionSatisfiesLemma_2_25`, `PartitionCopy`, `PartitionColorVertex`, `CellGetMinNode`, `CellGetNumberOfNodes`, `CellIntersectWithSet`, `CtPartClear`, `CtPartFill`, `CtPartINCHI_CANON_INFINITY`, `CtPartCompare`, `CtFullCompare`, `CtFullCompareLayers`, `CtCompareLayersGetFirstDiff`, `CtPartCompareLayers`, `UpdateCompareLayers`, `CtPartCopy`, `CtFullCopy`, `TranspositionGetMcrAndFixSetAndUnorderedPartition`, `SetBitCreate`, `GetOneAdditionalLayer`, `bInchiTimeIsOver` | 2 |
| `SetInitialRanks2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:4674` | external | reachable | `inchi_qsort`, `CompAtomInvariants2Only` | 0 |
| `FillOutAtomInvariant2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:4719` | external | reachable | none | 0 |
| `CleanNumH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:4881` | external | reachable | none | 0 |
| `CleanCt` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:4902` | external | reachable | none | 0 |
| `CleanIsoSortKeys` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:4915` | external | reachable | none | 0 |
| `DeAllocBCN` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:4960` | external | reachable | `PartitionFree`, `FreeNeighList` | 0 |
| `GetBaseCanonRanking` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:5105` | external | reachable | `CTableCreate`, `CTableFree`, `CtFullCopy`, `CanonGraph`, `SetInitialRanks2`, `FillOutAtomInvariant2`, `CleanNumH`, `CleanCt`, `CleanIsoSortKeys`, `FixCanonEquivalenceInfo`, `make_iso_sort_key`, `DifferentiateRanks2`, `DifferentiateRanks4`, `CreateNeighList`, `FreeNeighList` | 0 |
| `InchiClock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:111` | internal | reachable | none | 0 |
| `FillMaxMinClock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:128` | internal | reachable | none | 0 |
| `InchiTimeGet` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:146` | external | reachable | `InchiClock` | 0 |
| `InchiTimeMsecDiff` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:151` | external | reachable | `FillMaxMinClock` | 0 |
| `InchiTimeElapsed` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:223` | external | reachable | `InchiTimeGet`, `InchiTimeMsecDiff` | 0 |
| `InchiTimeAddMsec` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:234` | external | reachable | `FillMaxMinClock` | 0 |
| `bInchiTimeIsOver` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:257` | external | reachable | `InchiClock`, `FillMaxMinClock` | 0 |
| `GetCanonLengths` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:418` | external | reachable | none | 0 |
| `DeAllocateCS` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:491` | external | reachable | `FreeNeighList` | 0 |
| `AllocateCS` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:575` | external | reachable | `DeAllocateCS` | 0 |
| `FillIsotopicAtLinearCT` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:780` | external | reachable | none | 0 |
| `FillTautLinearCT2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:858` | external | reachable | `inchi_qsort` | 0 |
| `UpdateFullLinearCT` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:1050` | external | reachable | `insertions_sort` | 0 |
| `FixCanonEquivalenceInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:1278` | external | reachable | `SortedEquInfoToRanks`, `SortedRanksToEquInfo`, `inchi_qsort` | 0 |
| `Canon_INChI3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:1370` | external | reachable | `InchiTimeGet`, `InchiTimeMsecDiff`, `FillIsotopicAtLinearCT`, `FillTautLinearCT2`, `FixCanonEquivalenceInfo`, `SwitchAtomStereoAndIsotopicStereo`, `SetCtToIsotopicStereo`, `SetCtToNonIsotopicStereo`, `InvertStereo`, `FillOutStereoParities`, `CompareLinCtStereo`, `CurTreeAlloc`, `CurTreeFree`, `CurTreeSetPos`, `map_stereo_bonds4`, `inchi_swap`, `SortTautomerGroupsAndEndpoints` | 0 |
| `Canon_INChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:2727` | external | reachable | `Canon_INChI3` | 0 |
| `find_atoms_with_parity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:155` | external | reachable | `find_atoms_with_parity` | 0 |
| `SetHalfStereoBondIllDefPariy` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:189` | external | reachable | none | 0 |
| `RemoveHalfStereoBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:208` | external | reachable | none | 0 |
| `SetOneStereoBondIllDefParity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:242` | external | reachable | `SetHalfStereoBondIllDefPariy` | 0 |
| `RemoveOneStereoBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:271` | external | reachable | `RemoveHalfStereoBond` | 0 |
| `RemoveOneStereoCenter` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:300` | external | reachable | none | 0 |
| `UnmarkNonStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:322` | external | reachable | `find_atoms_with_parity`, `RemoveHalfStereoBond`, `insertions_sort` | 0 |
| `FillSingleStereoDescriptors` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:525` | external | reachable | `HalfStereoBondParity`, `insertions_sort` | 0 |
| `SwitchAtomStereoAndIsotopicStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:705` | external | reachable | `inchi_swap` | 0 |
| `SetCtToIsotopicStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:729` | external | reachable | none | 0 |
| `SetCtToNonIsotopicStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:746` | external | reachable | none | 0 |
| `FillAllStereoDescriptors` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:766` | external | reachable | `FillSingleStereoDescriptors` | 0 |
| `SetKnownStereoBondParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:801` | external | reachable | `insertions_sort` | 0 |
| `MarkKnownEqualStereoBondParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:1093` | external | reachable | none | 0 |
| `GetNextNeighborAndRank` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:1419` | external | reachable | none | 0 |
| `GetAndCheckNextNeighbors` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:1456` | external | reachable | `GetNextNeighborAndRank` | 0 |
| `PathsHaveIdenticalKnownParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:1530` | external | reachable | `GetAndCheckNextNeighbors`, `PathsHaveIdenticalKnownParities` | 0 |
| `RemoveKnownNonStereoBondParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:1608` | external | reachable | `RemoveOneStereoBond`, `PathsHaveIdenticalKnownParities` | 0 |
| `SetKnownStereoCenterParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:1714` | external | reachable | `insertions_sort` | 0 |
| `RemoveKnownNonStereoCenterParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:1833` | external | reachable | `PathsHaveIdenticalKnownParities`, `insertions_sort` | 0 |
| `MarkKnownEqualStereoCenterParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:1942` | external | reachable | none | 0 |
| `InvertStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:2016` | external | reachable | none | 0 |
| `FillOutStereoParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:2127` | external | reachable | `UnmarkNonStereo`, `FillAllStereoDescriptors`, `SetKnownStereoBondParities`, `MarkKnownEqualStereoBondParities`, `RemoveKnownNonStereoBondParities`, `SetKnownStereoCenterParities`, `RemoveKnownNonStereoCenterParities`, `MarkKnownEqualStereoCenterParities` | 0 |
| `GetStereoNeighborPos` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:2197` | external | reachable | none | 0 |
| `GetStereoBondParity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:2219` | external | reachable | `HalfStereoBondParity` | 0 |
| `GetPermutationParity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:2296` | external | reachable | `insertions_sort` | 0 |
| `GetStereoCenterParity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:2340` | external | reachable | `insertions_sort` | 0 |
| `ErrMsg` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichierr.c:54` | external | reachable | none | 0 |
| `AddErrorMessage` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichierr.c:106` | external | reachable | `already_have_this_message` | 0 |
| `already_have_this_message` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichierr.c:160` | internal | reachable | none | 0 |
| `make_iso_sort_key` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiisot.c:47` | external | reachable | none | 0 |
| `set_atom_iso_sort_keys` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiisot.c:67` | external | reachable | `make_iso_sort_key` | 0 |
| `MakeHillFormulaString` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:121` | external | reachable | `inchi_strbuf_printf` | 0 |
| `GetHillFormulaIndexLength` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:144` | external | reachable | none | 0 |
| `GetHillFormulaCounts` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:159` | external | reachable | `GetHillFormulaIndexLength`, `get_element_or_pseudoelement_symbol` | 0 |
| `AddElementAndCount` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:279` | external | reachable | none | 0 |
| `MakeHillFormula` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:316` | external | reachable | `AddElementAndCount`, `get_element_or_pseudoelement_symbol` | 0 |
| `AllocateAndFillHillFormula` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:402` | external | reachable | `GetHillFormulaCounts`, `MakeHillFormula` | 0 |
| `Copy2StereoBondOrAllene` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:451` | external | reachable | none | 0 |
| `CopyLinearCTStereoToINChIStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:566` | external | reachable | `Copy2StereoBondOrAllene`, `CompareLinCtStereoDble` | 0 |
| `MarkAmbiguousStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:713` | external | reachable | none | 0 |
| `UnmarkAllUndefinedUnknownStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:823` | external | reachable | none | 0 |
| `WriteCoord` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:890` | external | reachable | none | 0 |
| `FillOutINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:1072` | external | reachable | `AddErrorMessage`, `AllocateAndFillHillFormula`, `CopyLinearCTStereoToINChIStereo`, `MarkAmbiguousStereo`, `UnmarkAllUndefinedUnknownStereo`, `switch_ptrs`, `inchi_qsort`, `SortTautomerGroupsAndEndpoints`, `get_unusual_el_valence` | 0 |
| `inp2spATOM` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:119` | external | reachable | `get_periodic_table_number` | 0 |
| `GetElementAndCount` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:173` | external | reachable | none | 0 |
| `CompareHillFormulas` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:241` | external | reachable | `GetElementAndCount` | 0 |
| `CompareHillFormulasNoH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:272` | external | reachable | `GetElementAndCount` | 0 |
| `CompareTautNonIsoPartOfINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:316` | external | reachable | none | 0 |
| `CompINChITautVsNonTaut` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:341` | external | reachable | `CompareHillFormulasNoH`, `CompareTautNonIsoPartOfINChI`, `CompareInchiStereo` | 0 |
| `GetSp3RelRacAbs` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:593` | external | reachable | none | 0 |
| `CompINChILayers` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:645` | external | reachable | `CompareHillFormulasNoH`, `GetSp3RelRacAbs`, `Eql_INChI_Stereo` | 0 |
| `INChI_SegmentAction` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:1486` | external | reachable | none | 0 |
| `MarkUnusedAndEmptyLayers` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:1525` | external | reachable | none | 0 |
| `CompareInchiStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:1607` | external | reachable | none | 0 |
| `CompINChI2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:1712` | external | reachable | `CompareHillFormulas`, `CompareHillFormulasNoH`, `CompareTautNonIsoPartOfINChI`, `CompareInchiStereo` | 0 |
| `CompINChINonTaut2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2042` | external | reachable through `SortAndPrintINChI` `qsort` callback | `CompINChI2` | 0 |
| `CompINChITaut2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2064` | external | reachable through `SortAndPrintINChI` `qsort` callback | `CompINChI2` | 0 |
| `mystrrev` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2090` | external | reachable | none | 0 |
| `CompareDfsDescendants4CT` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2119` | internal | active, export-unreachable | none | 0 |
| `GetDfsOrder4CT` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2154` | external | reachable | `insertions_sort`, `CreateNeighListFromLinearCT`, `FreeNeighList` | 0 |
| `GetInpStructErrorType` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2480` | external | reachable | none | 0 |
| `ProcessStructError` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2511` | external | active, export-unreachable | `inchi_ios_print`, `inchi_ios_eprint`, `OutputINChIPlainError` | 0 |
| `CompareReversedStereoINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2555` | external | reachable | none | 0 |
| `CompareReversedStereoINChI2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2635` | external | reachable | none | 0 |
| `CompareReversedINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2936` | external | reachable | `CompareReversedStereoINChI` | 0 |
| `CompareIcr` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:3189` | external | reachable | none | 0 |
| `CompareReversedINChI2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:3262` | external | reachable | `CompareHillFormulasNoH`, `CompareReversedStereoINChI2`, `insertions_sort_AT_RANK` | 0 |
| `Create_INChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:3707` | external | reachable | `mark_alt_bonds_and_taut_groups`, `DeAllocBCN`, `GetBaseCanonRanking`, `GetCanonLengths`, `DeAllocateCS`, `AllocateCS`, `Canon_INChI`, `set_atom_iso_sort_keys`, `FillOutINChI`, `inp2spATOM`, `CheckCanonNumberingCorrectness`, `MarkRingSystemsInp`, `CreateNeighList`, `FreeNeighList`, `set_stereo_parity`, `free_t_group_info`, `make_a_copy_of_t_group_info`, `set_tautomer_iso_sort_keys`, `CountTautomerGroups`, `FreeInpAtom`, `add_DT_to_num_H`, `remove_terminal_HDT` | 0 |
| `CheckCanonNumberingCorrectness` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:6230` | external | reachable | `UpdateFullLinearCT` | 0 |
| `All_SC_Same` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:53` | external | reachable | none | 0 |
| `Next_SC_At_CanonRank2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:99` | external | reachable | none | 0 |
| `CompareLinCtStereoDble` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:178` | external | reachable | none | 0 |
| `CompareLinCtStereoCarb` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:223` | external | reachable | none | 0 |
| `CompareLinCtStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:262` | external | reachable | `CompareLinCtStereoDble`, `CompareLinCtStereoCarb` | 0 |
| `CompareLinCtStereoAtomToValues` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:287` | external | reachable | none | 0 |
| `bUniqueAtNbrFromMappingRank` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:316` | external | reachable | none | 0 |
| `nGetMcr` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:336` | external | reachable | none | 0 |
| `nJoin2Mcrs` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:366` | external | reachable | `nGetMcr` | 0 |
| `All_SB_Same` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:393` | external | reachable | none | 0 |
| `Next_SB_At_CanonRanks2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:518` | external | reachable | none | 0 |
| `NextStereoParity2Test` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:659` | external | reachable | none | 0 |
| `CompareLinCtStereoDoubleToValues` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:764` | external | reachable | none | 0 |
| `SetUseAtomForStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:804` | external | reachable | none | 0 |
| `CurTreeAlloc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:823` | external | reachable | none | 0 |
| `CurTreeReAlloc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:850` | external | reachable | none | 0 |
| `CurTreeFree` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:871` | external | reachable | none | 0 |
| `CurTreeAddRank` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:881` | external | reachable | `CurTreeReAlloc` | 0 |
| `CurTreeIsLastRank` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:902` | external | reachable | none | 0 |
| `CurTreeRemoveLastRankIfNoAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:920` | external | reachable | `CurTreeRemoveLastRank` | 0 |
| `CurTreeAddAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:935` | external | reachable | `CurTreeReAlloc` | 0 |
| `CurTreeKeepLastAtomsOnly` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:960` | external | reachable | `CurTreeKeepLastAtomsOnly` | 0 |
| `CurTreeRemoveIfLastAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:993` | external | reachable | none | 0 |
| `CurTreeGetPos` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:1011` | external | reachable | none | 0 |
| `CurTreeSetPos` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:1023` | external | reachable | none | 0 |
| `CurTreeRemoveLastRank` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:1036` | external | reachable | none | 0 |
| `CurTreeIsLastAtomEqu` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:1054` | external | reachable | none | 0 |
| `SortedEquInfoToRanks` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:148` | external | reachable | none | 0 |
| `SortedRanksToEquInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:199` | external | reachable | none | 0 |
| `switch_ptrs` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:230` | external | reachable | none | 0 |
| `SetNewRanksFromNeighLists3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:241` | external | reachable | `insertions_sort_AT_NUMBERS`, `CompareNeighListLex` | 0 |
| `SetNewRanksFromNeighLists4` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:308` | external | reachable | `insertions_sort_AT_NUMBERS`, `CompareNeighListLexUpToMaxRank` | 0 |
| `SetNewRanksFromNeighLists` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:380` | external | reachable | `inchi_qsort`, `insertions_sort`, `CompNeighListRanks` | 0 |
| `SortNeighListsBySymmAndCanonRank` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:434` | external | reachable | `insertions_sort_NeighListBySymmAndCanonRank` | 0 |
| `SortNeighLists2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:448` | external | reachable | `insertions_sort_NeighList_AT_NUMBERS` | 0 |
| `SortNeighLists3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:482` | external | reachable | `insertions_sort_NeighList_AT_NUMBERS3` | 0 |
| `DifferentiateRanks2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:518` | external | reachable | `switch_ptrs`, `SetNewRanksFromNeighLists`, `SortNeighLists2`, `inchi_qsort`, `insertions_sort` | 0 |
| `DifferentiateRanks3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:561` | external | reachable | `switch_ptrs`, `SetNewRanksFromNeighLists3`, `SortNeighLists3` | 0 |
| `DifferentiateRanks4` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:602` | external | reachable | `switch_ptrs`, `SetNewRanksFromNeighLists4`, `SortNeighLists3` | 0 |
| `DifferentiateRanksBasic` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:637` | external | reachable | `switch_ptrs`, `SetNewRanksFromNeighLists`, `SortNeighLists2`, `inchi_qsort`, `insertions_sort` | 0 |
| `NumberOfTies` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:680` | external | reachable | none | 0 |
| `HalfStereoBondParity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:802` | external | reachable | none | 0 |
| `parity_of_mapped_half_bond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:958` | external | reachable | none | 0 |
| `parity_of_mapped_atom2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:1168` | external | reachable | `insertions_sort` | 0 |
| `ClearPreviousMappings` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:1334` | external | reachable | none | 0 |
| `map_an_atom2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:1350` | external | reachable | `DifferentiateRanks2`, `NumberOfTies` | 0 |
| `might_change_other_atom_parity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:1461` | external | reachable | none | 0 |
| `DeAllocateForNonStereoRemoval` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:1504` | external | reachable | `FreeNeighList` | 0 |
| `AllocateForNonStereoRemoval` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:1551` | external | reachable | `SortNeighListsBySymmAndCanonRank`, `DeAllocateForNonStereoRemoval`, `CreateNeighList` | 0 |
| `GetMinNewRank` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:1587` | external | reachable | none | 0 |
| `BreakNeighborsTie` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:1609` | external | reachable | `SortedEquInfoToRanks`, `SortNeighListsBySymmAndCanonRank`, `DifferentiateRanksBasic`, `GetMinNewRank`, `insertions_sort` | 0 |
| `CheckNextSymmNeighborsAndBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:2060` | external | reachable | none | 0 |
| `CreateCheckSymmPaths` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:2191` | external | reachable | `CheckNextSymmNeighborsAndBonds`, `CreateCheckSymmPaths` | 0 |
| `CalculatedPathsParitiesAreIdentical` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:2399` | external | reachable | `GetStereoNeighborPos`, `GetStereoBondParity`, `GetPermutationParity`, `GetStereoCenterParity` | 0 |
| `RemoveCalculatedNonStereoBondParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:3099` | external | reachable | `SetOneStereoBondIllDefParity`, `RemoveOneStereoBond`, `BreakNeighborsTie`, `CreateCheckSymmPaths`, `CalculatedPathsParitiesAreIdentical` | 0 |
| `RemoveCalculatedNonStereoCenterParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:3356` | external | reachable | `RemoveOneStereoCenter`, `BreakNeighborsTie`, `CreateCheckSymmPaths`, `CalculatedPathsParitiesAreIdentical` | 0 |
| `RemoveCalculatedNonStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:3672` | external | reachable | `DeAllocateForNonStereoRemoval`, `AllocateForNonStereoRemoval`, `RemoveCalculatedNonStereoBondParities`, `RemoveCalculatedNonStereoCenterParities` | 0 |
| `map_stereo_bonds4` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap4.c:83` | external | reachable | `All_SB_Same`, `Next_SB_At_CanonRanks2`, `NextStereoParity2Test`, `CompareLinCtStereoDoubleToValues`, `SetUseAtomForStereo`, `CurTreeAddRank`, `CurTreeIsLastRank`, `CurTreeRemoveLastRankIfNoAtoms`, `CurTreeAddAtom`, `CurTreeKeepLastAtomsOnly`, `CurTreeRemoveIfLastAtom`, `CurTreeGetPos`, `CurTreeSetPos`, `CurTreeRemoveLastRank`, `CurTreeIsLastAtomEqu`, `parity_of_mapped_half_bond`, `ClearPreviousMappings`, `map_an_atom2`, `map_stereo_bonds4`, `map_stereo_atoms4` | 0 |
| `map_stereo_atoms4` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap4.c:1126` | external | reachable | `bInchiTimeIsOver`, `All_SC_Same`, `Next_SC_At_CanonRank2`, `CompareLinCtStereoAtomToValues`, `bUniqueAtNbrFromMappingRank`, `nGetMcr`, `nJoin2Mcrs`, `NextStereoParity2Test`, `CurTreeAddRank`, `CurTreeIsLastRank`, `CurTreeRemoveLastRankIfNoAtoms`, `CurTreeAddAtom`, `CurTreeKeepLastAtomsOnly`, `CurTreeRemoveIfLastAtom`, `CurTreeGetPos`, `CurTreeSetPos`, `CurTreeRemoveLastRank`, `CurTreeIsLastAtomEqu`, `parity_of_mapped_atom2`, `ClearPreviousMappings`, `map_an_atom2`, `might_change_other_atom_parity`, `RemoveCalculatedNonStereo`, `map_stereo_atoms4`, `BreakAllTies` | 2 |
| `MarkRingSystemsInp` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:59` | external | reachable | none | 0 |
| `UnMarkDisconnectedComponents` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:652` | external | reachable | none | 0 |
| `UnMarkOtherIndicators` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:681` | external | reachable | none | 0 |
| `UnMarkOneComponent` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:698` | external | reachable | none | 0 |
| `set_R2C_el_numbers` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:712` | external | active, export-unreachable | none | 0 |
| `subtract_DT_from_num_H` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:734` | external | reachable | none | 0 |
| `add_inp_ATOM` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:749` | external | reachable | none | 0 |
| `mark_arom_bonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:784` | external | reachable | `mark_alt_bonds_and_taut_groups` | 0 |
| `cmp_r2c_atpair` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:800` | external | reachable | none | 0 |
| `has_atom_pair_seq` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:815` | external | reachable | none | 0 |
| `has_atom_pair` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:835` | external | active, export-unreachable | `cmp_r2c_atpair` | 0 |
| `mark_atoms_ap` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:871` | external | reachable | `has_atom_pair_seq`, `mark_atoms_ap` | 0 |
| `mark_atoms_cFlags` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1368` | external | reachable | `mark_atoms_cFlags` | 0 |
| `unmark_atoms_cFlags` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1388` | external | reachable | `unmark_atoms_cFlags` | 0 |
| `is_C_or_S_DB_O` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1410` | external | reachable | none | 0 |
| `is_C_DB_O` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1434` | external | reachable | none | 0 |
| `is_C_unsat_not_arom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1457` | external | reachable | none | 0 |
| `is_Aryl` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1486` | external | active, export-unreachable | none | 0 |
| `is_Saturated_C` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1509` | external | active, export-unreachable | none | 0 |
| `is_C_Alk` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1517` | external | active, export-unreachable | none | 0 |
| `is_Phenyl` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1559` | external | reachable | none | 0 |
| `is_PentaFluoroPhenyl` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1615` | external | reachable | none | 0 |
| `is_Methyl` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1679` | external | reachable | none | 0 |
| `is_Ethyl` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1694` | external | reachable | `is_Methyl` | 0 |
| `is_Methyl_or_Etyl` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1709` | external | reachable | `is_Methyl`, `is_Ethyl` | 0 |
| `is_Si_IV` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1727` | external | reachable | none | 0 |
| `is_P_TB_N` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1740` | external | active, export-unreachable | none | 0 |
| `get_CO_opposite` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1763` | external | active, export-unreachable | none | 0 |
| `is_DERIV_RING_DMOX_DEOX_O` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1802` | external | reachable | `is_Methyl_or_Etyl` | 0 |
| `is_DERIV_RING_DMOX_DEOX_N` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1962` | external | reachable | `is_Methyl_or_Etyl` | 0 |
| `is_DERIV_RING2_PRRLDD_PPRDN` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:2147` | external | reachable | none | 0 |
| `check_arom_chain` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:2281` | external | reachable | none | 0 |
| `is_Dansyl` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:2334` | external | reachable | `is_Methyl`, `check_arom_chain` | 0 |
| `is_possibly_deriv_neigh` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:2610` | external | reachable | `is_C_or_S_DB_O`, `is_C_unsat_not_arom`, `is_Si_IV`, `is_deriv_chain2` | 0 |
| `get_traversed_deriv_type` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:2677` | external | reachable | `mark_atoms_cFlags`, `unmark_atoms_cFlags`, `is_C_or_S_DB_O`, `is_C_DB_O`, `is_Phenyl`, `is_PentaFluoroPhenyl`, `is_Methyl`, `is_Ethyl`, `is_Si_IV`, `is_DERIV_RING_DMOX_DEOX_O`, `is_DERIV_RING_DMOX_DEOX_N`, `is_DERIV_RING2_PRRLDD_PPRDN`, `is_Dansyl`, `is_possibly_deriv_neigh`, `is_silyl2`, `is_CF3_or_linC3F7a`, `is_el_a_metal` | 0 |
| `add_to_da` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:3245` | external | reachable | none | 0 |
| `mark_atoms_deriv` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:3332` | external | reachable | `get_traversed_deriv_type`, `add_to_da`, `mark_atoms_deriv` | 0 |
| `count_one_bond_atoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:3536` | external | reachable | `mark_atoms_deriv` | 0 |
| `is_silyl` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:3633` | external | reachable | none | 0 |
| `is_silyl2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:3738` | external | reachable | none | 0 |
| `is_nButyl` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:3831` | external | reachable | none | 0 |
| `is_Me_or_Et` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:3862` | external | reachable | none | 0 |
| `is_CF3_or_linC3F7a` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:4000` | external | reachable | `is_CF3_or_linC3F7` | 0 |
| `is_CF3_or_linC3F7` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:4017` | external | reachable | `is_in_the_list` | 0 |
| `is_phenyl` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:4103` | external | reachable | none | 0 |
| `is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:4136` | external | reachable | `is_nButyl`, `is_Me_or_Et`, `is_phenyl`, `underiv_list_add`, `underiv_list_get_last`, `is_in_the_list` | 0 |
| `is_deriv_chain2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:4483` | external | reachable | `is_silyl`, `is_CF3_or_linC3F7`, `underiv_list_add`, `underiv_list_get_last`, `is_in_the_list` | 0 |
| `is_deriv_chain` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:4980` | external | reachable | `is_deriv_chain2` | 0 |
| `is_deriv_chain_or_ring` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:4996` | external | reachable | `is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR` | 0 |
| `remove_deriv` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5046` | external | reachable | none | 0 |
| `remove_deriv_mark` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5108` | external | reachable | none | 0 |
| `underiv_buf_clear` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5147` | external | reachable | none | 0 |
| `underiv_list_add` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5156` | external | reachable | none | 0 |
| `underiv_list_get_last` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5179` | external | reachable | none | 0 |
| `underiv_compare` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5197` | external | reachable | none | 0 |
| `underiv_list_add_two_cuts` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5204` | external | reachable | `underiv_list_add` | 0 |
| `sort_merge_underiv` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5281` | external | reachable | `underiv_list_add`, `underiv_compare` | 0 |
| `eliminate_deriv_not_in_list` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5351` | external | reachable | `is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR`, `is_deriv_chain`, `underiv_buf_clear`, `underiv_list_add`, `underiv_list_add_two_cuts` | 0 |
| `make_single_cut` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5495` | external | reachable | `DisconnectInpAtBond`, `is_in_the_list` | 0 |
| `fill_out_bond_cuts` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5628` | external | reachable | `cmp_r2c_atpair` | 0 |
| `mark_deriv_agents` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5743` | external | reachable | `UnMarkOtherIndicators`, `mark_atoms_ap` | 0 |
| `replace_arom_bonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5809` | external | reachable | `is_in_the_list` | 0 |
| `add_explicit_H` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5868` | external | reachable | none | 0 |
| `OAD_Edit_Underivatize` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5934` | external | reachable | `MarkRingSystemsInp`, `UnMarkDisconnectedComponents`, `UnMarkOtherIndicators`, `UnMarkOneComponent`, `mark_arom_bonds`, `mark_atoms_ap`, `count_one_bond_atoms`, `is_silyl2`, `is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR`, `is_deriv_chain`, `is_deriv_chain_or_ring`, `remove_deriv`, `remove_deriv_mark`, `sort_merge_underiv`, `eliminate_deriv_not_in_list`, `make_single_cut`, `fill_out_bond_cuts`, `mark_deriv_agents`, `replace_arom_bonds`, `add_explicit_H`, `OAD_Edit_MergeComponentsAndRecreateOAD`, `free_underiv_temp_data`, `remove_cut_derivs`, `CreateInpAtomData`, `UnMarkRingSystemsInp`, `add_DT_to_num_H`, `remove_terminal_HDT`, `MarkDisconnectedComponents`, `ExtractConnectedComponent` | 0 |
| `detect_r2c_Zatom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:7159` | external | reachable | none | 0 |
| `cut_ring_to_chain` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:7296` | external | reachable | `DisconnectInpAtBond`, `is_in_the_list` | 0 |
| `Ring2Chain` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:7362` | external | reachable | `MarkRingSystemsInp`, `UnMarkDisconnectedComponents`, `UnMarkOtherIndicators`, `UnMarkOneComponent`, `mark_arom_bonds`, `mark_atoms_ap`, `detect_r2c_Zatom`, `cut_ring_to_chain`, `OAD_Edit_MergeComponentsAndRecreateOAD`, `FreeInpAtomData`, `CreateInpAtomData`, `UnMarkRingSystemsInp`, `add_DT_to_num_H`, `remove_terminal_HDT`, `MarkDisconnectedComponents`, `ExtractConnectedComponent` | 0 |
| `OAD_Edit_MergeComponentsAndRecreateOAD` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:7563` | external | reachable | `UnMarkDisconnectedComponents`, `UnMarkOtherIndicators`, `UnMarkOneComponent`, `subtract_DT_from_num_H`, `add_inp_ATOM`, `UnMarkRingSystemsInp` | 0 |
| `free_underiv_temp_data` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:7632` | external | reachable | `FreeInpAtomData` | 0 |
| `remove_cut_derivs` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:7666` | external | reachable | `UnMarkDisconnectedComponents`, `UnMarkOtherIndicators`, `UnMarkOneComponent`, `add_inp_ATOM`, `FreeInpAtomData`, `CreateInpAtomData`, `UnMarkRingSystemsInp`, `MarkDisconnectedComponents`, `ExtractConnectedComponent` | 0 |
| `set_common_options_by_parg` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiparm.c:131` | external | reachable | `mystrncpy`, `lrtrim`, `inchi_memicmp`, `inchi_stricmp` | 0 |
| `ReadCommandLineParms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiparm.c:602` | external | reachable | `inchi_ios_eprint`, `set_common_options_by_parg`, `mystrncpy`, `lrtrim`, `inchi_memicmp`, `inchi_stricmp`, `inchi__strdup` | 0 |
| `PrintInputParms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiparm.c:2130` | external | reachable | `inchi_ios_eprint` | 0 |
| `HelpCommandLineParms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiparm.c:2700` | external | reachable | `inchi_ios_print_nodisplay` | 0 |
| `OutputINChIPlainError` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:514` | external | active, export-unreachable | `inchi_ios_print` | 0 |
| `EquString` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:564` | external | reachable | none | 0 |
| `OutputINChI2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:966` | external | reachable | `OutputINChI1` | 0 |
| `OutputINChI1` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:1043` | external | reachable | `inchi_ios_print`, `inchi_ios_print_nodisplay`, `inchi_ios_eprint`, `inchi_strbuf_reset`, `CompINChITautVsNonTaut`, `CompINChILayers`, `MarkUnusedAndEmptyLayers`, `OutputINChI1`, `GetSaveOptLetters`, `set_line_separators`, `OutputINCHI_VersionAndKind`, `OutputINCHI_MainLayerFormula`, `OutputINCHI_MainLayerConnections`, `OutputINCHI_MainLayerHydrogens`, `OutputINCHI_ChargeAndRemovedAddedProtonsLayers`, `OutputINCHI_StereoLayer`, `OutputINCHI_IsotopicLayer`, `OutputINCHI_FixedHLayerWithSublayers`, `OutputINCHI_PolymerLayer`, `OutputAUXINFO_HeaderAndNormalization_type`, `OutputAUXINFO_OriginalNumbersAndEquivalenceClasses`, `OutputAUXINFO_TautomericGroupsEquivalence`, `OutputAUXINFO_Stereo`, `OutputAUXINFO_IsotopicInfo`, `OutputAUXINFO_ChargesRadicalsAndUnusualValences`, `OutputAUXINFO_ReversibilityInfo`, `OutputAUXINFO_PolymerInfo`, `Eql_INChI_Stereo`, `bHasOrigInfo`, `bHasEquString` | 0 |
| `szGetTag` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:2087` | external | reachable | none | 0 |
| `str_LineEnd` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:2183` | internal | reachable | `inchi_strbuf_update` | 0 |
| `CleanOrigCoord` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:2228` | internal | reachable | `lrtrim` | 0 |
| `WriteOrigCoord` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:2340` | internal | reachable | `CleanOrigCoord` | 0 |
| `WriteOrigAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:2410` | internal | reachable | `nBondsValenceInpAt`, `insertions_sort`, `needed_unusual_el_valence`, `get_atomic_mass_from_elnum`, `is_in_the_list` | 0 |
| `WriteOrigBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:2610` | internal | reachable | `AddErrorMessage`, `insertions_sort`, `get_opposite_sb_atom`, `is_el_a_metal`, `is_in_the_list` | 0 |
| `OrigStruct_FillOut` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:2933` | external | reachable | `WriteOrigCoord`, `WriteOrigAtoms`, `WriteOrigBonds` | 0 |
| `OrigStruct_Free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3051` | external | reachable | none | 0 |
| `GetSaveOptLetters` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3100` | internal | reachable | none | 0 |
| `set_line_separators` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3112` | external | reachable | none | 0 |
| `OutputINCHI_VersionAndKind` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3137` | internal | reachable | `inchi_ios_print_nodisplay`, `inchi_strbuf_reset`, `inchi_strbuf_printf` | 0 |
| `OutputINCHI_MainLayerFormula` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3170` | internal | reachable | `inchi_ios_print_nodisplay`, `inchi_strbuf_reset`, `szGetTag`, `str_LineEnd`, `MergeZzInHillFormula`, `str_HillFormula` | 0 |
| `OutputINCHI_MainLayerConnections` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3213` | internal | reachable | `inchi_ios_print_nodisplay`, `inchi_strbuf_reset`, `szGetTag`, `str_LineEnd`, `str_Connections` | 0 |
| `OutputINCHI_MainLayerHydrogens` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3251` | internal | reachable | `inchi_ios_print_nodisplay`, `inchi_strbuf_reset`, `INChI_SegmentAction`, `szGetTag`, `str_LineEnd`, `str_H_atoms` | 0 |
| `OutputINCHI_ChargeAndRemovedAddedProtonsLayers` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3293` | internal | reachable | `inchi_ios_print_nodisplay`, `inchi_strbuf_reset`, `inchi_strbuf_printf`, `INChI_SegmentAction`, `szGetTag`, `str_LineEnd`, `str_Charge2` | 0 |
| `OutputINCHI_StereoLayer` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3354` | internal | reachable | `inchi_ios_print_nodisplay`, `inchi_strbuf_reset`, `INChI_SegmentAction`, `szGetTag`, `str_LineEnd`, `MakeDelim`, `str_Sp2`, `str_Sp3`, `str_StereoAbsInv` | 0 |
| `OutputINCHI_IsotopicLayer` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3502` | internal | reachable | `inchi_ios_print_nodisplay`, `inchi_strbuf_reset`, `INChI_SegmentAction`, `szGetTag`, `str_LineEnd`, `MakeDelim`, `MakeIsoHString`, `str_IsoAtoms`, `str_IsoSp2`, `str_IsoSp3`, `str_IsoStereoAbsInv`, `bin_AuxTautTrans`, `str_AuxTautTrans` | 0 |
| `OutputINCHI_FixedHLayerWithSublayers` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3750` | internal | reachable | `inchi_ios_print_nodisplay`, `inchi_strbuf_reset`, `INChI_SegmentAction`, `szGetTag`, `str_LineEnd`, `MergeZzInHillFormula`, `str_HillFormula2`, `str_FixedH_atoms` | 0 |
| `OutputINCHI_PolymerLayer` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3880` | internal | reachable | `inchi_ios_print_nodisplay`, `inchi_strbuf_printf`, `OutputINCHI_PolymerLayer_SingleUnit`, `InternallyGetCanoNumsAndComponentNums`, `OAD_PolymerUnit_CreateCopy`, `OAD_PolymerUnit_Free`, `OAD_Polymer_PrepareWorkingSet`, `OAD_Polymer_SetAtProps` | 0 |
| `OutputINCHI_PolymerLayer_SingleUnit` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:4111` | internal | reachable | `inchi_strbuf_printf`, `IsBondAtomNumsLesser`, `inchi_sort_int_pair_ascending`, `print_sequence_of_nums_compressing_ranges`, `OAD_PolymerUnit_SetReopeningDetails`, `OAD_PolymerUnit_SortBackboneBondsAndSetSeniors`, `OAD_Polymer_IsFirstAtomRankLower`, `is_in_the_ilist` | 0 |
| `OutputAUXINFO_HeaderAndNormalization_type` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:4339` | internal | reachable | `inchi_ios_print`, `inchi_strbuf_reset`, `inchi_strbuf_printf`, `szGetTag`, `str_LineEnd` | 0 |
| `OutputAUXINFO_OriginalNumbersAndEquivalenceClasses` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:4396` | internal | reachable | `inchi_ios_print`, `inchi_strbuf_reset`, `szGetTag`, `str_LineEnd`, `str_AuxEqu`, `str_AuxNumb` | 0 |
| `OutputAUXINFO_TautomericGroupsEquivalence` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:4461` | internal | reachable | `inchi_ios_print`, `inchi_strbuf_reset`, `szGetTag`, `str_LineEnd`, `str_AuxTgroupEqu` | 0 |
| `OutputAUXINFO_Stereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:4500` | internal | reachable | `inchi_ios_print`, `inchi_strbuf_reset`, `szGetTag`, `str_LineEnd`, `str_AuxInvSp3`, `str_AuxInvSp3Numb` | 0 |
| `OutputAUXINFO_IsotopicInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:4560` | internal | reachable | `inchi_ios_print`, `inchi_strbuf_reset`, `szGetTag`, `str_LineEnd`, `str_AuxIsoNumb`, `str_AuxIsoEqu`, `str_AuxInvIsoSp3`, `str_AuxInvIsoSp3Numb`, `str_AuxIsoTgroupEqu` | 0 |
| `OutputAUXINFO_ChargesRadicalsAndUnusualValences` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:4717` | internal | reachable | `inchi_ios_print`, `inchi_strbuf_reset`, `szGetTag`, `str_LineEnd`, `str_AuxChargeRadVal` | 0 |
| `OutputAUXINFO_ReversibilityInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:4753` | internal | reachable | `inchi_ios_print`, `inchi_strbuf_reset`, `szGetTag` | 0 |
| `OutputAUXINFO_PolymerInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:4923` | internal | reachable | `inchi_ios_print`, `inchi_strbuf_reset`, `inchi_strbuf_printf`, `print_sequence_of_nums_compressing_ranges` | 0 |
| `IsBondAtomNumsLesser` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:5007` | internal | reachable | none | 0 |
| `EditINCHI_HidePolymerZz` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:5039` | external | reachable | `inchi_ios_print_nodisplay`, `inchi_strtol` | 0 |
| `InternallyGetCanoNumsAndComponentNums` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:5362` | internal | reachable | `inchi_strbuf_reset`, `str_AuxNumb` | 0 |
| `MergeZzInHillFormula` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:5435` | external | reachable | `inchi_strbuf_reset`, `inchi_strbuf_printf`, `MergeZzInStrHillFormulaComponent` | 0 |
| `MergeZzInStrHillFormulaComponent` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:5495` | internal | reachable | none | 0 |
| `inchi_sort_int_pair_ascending` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:5525` | internal | reachable | none | 0 |
| `Eql_INChI_Stereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:54` | external | reachable | none | 0 |
| `Eql_INChI_Isotopic` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:203` | external | reachable | none | 0 |
| `Eql_INChI_Aux_Equ` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:227` | external | reachable | `bHasEquString` | 0 |
| `Eql_INChI_Aux_Num` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:305` | external | reachable | none | 0 |
| `bHasOrigInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:369` | external | reachable | none | 0 |
| `EqlOrigInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:387` | external | reachable | `bHasOrigInfo` | 0 |
| `bHasEquString` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:398` | external | reachable | none | 0 |
| `MakeMult` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:431` | external | reachable | `inchi_strbuf_printf`, `MakeAbcNumber`, `MakeDecNumber` | 0 |
| `MakeDelim` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:476` | external | reachable | `inchi_strbuf_printf` | 0 |
| `MakeEqStr` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:506` | external | reachable | `inchi_strbuf_printf`, `MakeDecNumber` | 0 |
| `MakeCtStringNew` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:551` | external | reachable | `inchi_strbuf_printf`, `GetDfsOrder4CT`, `MakeAbcNumber`, `MakeDecNumber` | 0 |
| `MakeCtStringOld` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:703` | external | reachable | `inchi_strbuf_printf`, `MakeAbcNumber`, `MakeDecNumber` | 0 |
| `MakeHString` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:789` | external | reachable | `inchi_strbuf_printf`, `MakeAbcNumber`, `MakeDecNumber` | 0 |
| `MakeCtString` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:1069` | external | reachable | `MakeCtStringNew`, `MakeCtStringOld` | 0 |
| `MakeTautString` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:1111` | external | reachable | `inchi_strbuf_printf`, `MakeAbcNumber`, `MakeDecNumber` | 0 |
| `MakeCRVString` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:1308` | external | reachable | `inchi_strbuf_printf`, `MakeAbcNumber`, `MakeDecNumber` | 0 |
| `MakeEquString` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:1559` | external | reachable | `inchi_strbuf_printf`, `MakeAbcNumber`, `MakeDecNumber` | 0 |
| `MakeIsoAtomString` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:1645` | external | reachable | `inchi_strbuf_printf`, `MakeDecNumber` | 1 |
| `MakeIsoTautString` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:1803` | external | reachable | `inchi_strbuf_printf`, `MakeDecNumber` | 1 |
| `MakeIsoHString` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:1924` | external | reachable | `inchi_strbuf_printf`, `MakeDecNumber` | 0 |
| `MakeStereoString` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:2019` | external | reachable | `inchi_strbuf_printf`, `MakeDecNumber` | 1 |
| `MakeAbcNumber` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:2143` | external | reachable | `mystrrev` | 0 |
| `abctol` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:2203` | internal | reachable | none | 0 |
| `inchi_strtol` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:2269` | external | reachable | `abctol` | 0 |
| `inchi_strtod` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:2288` | external | reachable | none | 0 |
| `MakeDecNumber` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:2302` | external | reachable | `mystrrev` | 0 |
| `print_sequence_of_nums_compressing_ranges` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:2376` | external | reachable | `inchi_strbuf_printf` | 0 |
| `str_HillFormula` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:60` | external | reachable | `MakeHillFormulaString`, `MakeMult`, `MakeDelim` | 0 |
| `str_HillFormula2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:137` | external | reachable | `MakeHillFormulaString`, `MakeMult`, `MakeDelim` | 0 |
| `str_Connections` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:224` | external | reachable | `MakeMult`, `MakeDelim`, `MakeCtStringNew` | 0 |
| `str_H_atoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:309` | external | reachable | `CompareTautNonIsoPartOfINChI`, `MakeMult`, `MakeDelim`, `MakeHString`, `MakeTautString` | 0 |
| `str_Charge2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:427` | external | reachable | `inchi_strbuf_printf`, `EquString`, `MakeMult`, `MakeDelim`, `MakeEqStr` | 0 |
| `str_FixedH_atoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:618` | external | reachable | `MakeMult`, `MakeDelim`, `MakeHString` | 0 |
| `str_Sp2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:725` | external | reachable | `EquString`, `Eql_INChI_Stereo`, `MakeMult`, `MakeDelim`, `MakeEqStr`, `MakeStereoString` | 0 |
| `str_Sp3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:957` | external | reachable | `EquString`, `Eql_INChI_Stereo`, `MakeMult`, `MakeDelim`, `MakeEqStr`, `MakeStereoString` | 0 |
| `str_StereoAbsInv` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:1191` | external | reachable | `MakeDelim` | 0 |
| `str_IsoAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:1229` | external | reachable | `EquString`, `Eql_INChI_Isotopic`, `MakeMult`, `MakeDelim`, `MakeEqStr`, `MakeIsoAtomString`, `MakeIsoTautString` | 0 |
| `str_IsoSp2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:1479` | external | reachable | `EquString`, `Eql_INChI_Stereo`, `MakeMult`, `MakeDelim`, `MakeEqStr`, `MakeStereoString` | 0 |
| `str_IsoSp3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:1768` | external | reachable | `EquString`, `Eql_INChI_Stereo`, `MakeMult`, `MakeDelim`, `MakeEqStr`, `MakeStereoString` | 0 |
| `str_IsoStereoAbsInv` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:2058` | external | reachable | `MakeDelim` | 0 |
| `str_AuxEqu` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:2106` | external | reachable | `EquString`, `Eql_INChI_Aux_Equ`, `bHasEquString`, `MakeMult`, `MakeDelim`, `MakeEqStr`, `MakeEquString` | 0 |
| `str_AuxInvSp3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:2315` | external | reachable | `EquString`, `Eql_INChI_Stereo`, `MakeMult`, `MakeDelim`, `MakeEqStr`, `MakeStereoString` | 0 |
| `str_AuxInvSp3Numb` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:2584` | external | reachable | `EquString`, `Eql_INChI_Aux_Num`, `MakeDelim`, `MakeEqStr`, `MakeCtString` | 0 |
| `str_AuxIsoNumb` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:2769` | external | reachable | `EquString`, `Eql_INChI_Aux_Num`, `MakeDelim`, `MakeEqStr`, `MakeCtString` | 0 |
| `str_AuxIsoEqu` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:2962` | external | reachable | `EquString`, `Eql_INChI_Aux_Equ`, `bHasEquString`, `MakeMult`, `MakeDelim`, `MakeEqStr`, `MakeEquString` | 0 |
| `str_AuxInvIsoSp3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:3213` | external | reachable | `EquString`, `Eql_INChI_Stereo`, `MakeMult`, `MakeDelim`, `MakeEqStr`, `MakeStereoString` | 0 |
| `str_AuxInvIsoSp3Numb` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:3575` | external | reachable | `EquString`, `Eql_INChI_Aux_Num`, `MakeDelim`, `MakeEqStr`, `MakeCtString` | 0 |
| `str_AuxNumb` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:3818` | external | reachable | `EquString`, `Eql_INChI_Aux_Num`, `MakeDelim`, `MakeEqStr`, `MakeCtString` | 0 |
| `str_AuxTgroupEqu` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:3953` | external | reachable | `Eql_INChI_Aux_Equ`, `bHasEquString`, `MakeMult`, `MakeDelim`, `MakeEquString` | 0 |
| `str_AuxChargeRadVal` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:4013` | external | reachable | `bHasOrigInfo`, `EqlOrigInfo`, `MakeMult`, `MakeDelim`, `MakeCRVString` | 0 |
| `bin_AuxTautTrans` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:4091` | external | reachable | none | 0 |
| `str_AuxTautTrans` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:4187` | external | reachable | `MakeDelim`, `MakeCtString` | 0 |
| `str_AuxIsoTgroupEqu` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:4236` | external | reachable | `EquString`, `Eql_INChI_Aux_Equ`, `bHasEquString`, `MakeMult`, `MakeDelim`, `MakeEqStr`, `MakeEquString` | 0 |
| `bIsCenterPointStrict` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:215` | external | reachable | `is_centerpoint_elem_strict`, `get_endpoint_valence` | 0 |
| `nGet14TautIn7MembAltRing` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:240` | external | reachable | `DFS_FindTautInARing` | 0 |
| `nGet14TautIn5MembAltRing` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:283` | external | reachable | `DFS_FindTautInARing` | 0 |
| `nGet12TautIn5MembAltRing` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:322` | external | reachable | `DFS_FindTautInARing` | 0 |
| `nGet15TautIn6MembAltRing` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:361` | external | reachable | `DFS_FindTautInARing` | 0 |
| `nGet15TautInAltPath` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:407` | external | reachable | `DFS_FindTautAltPath` | 0 |
| `DFS_FindTautInARing` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:458` | external | reachable | none | 2 |
| `DFS_FindTautAltPath` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:616` | external | reachable | none | 2 |
| `are_alt_bonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:778` | external | active, export-unreachable | none | 0 |
| `AddBondsPos` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:846` | external | active, export-unreachable | none | 0 |
| `AddEndPoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:897` | external | active, export-unreachable | none | 0 |
| `Check7MembTautRing` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:951` | external | active, export-unreachable | `bExistsAnyAltPath`, `are_alt_bonds`, `AddBondsPos`, `AddEndPoints`, `AddAtom2num`, `AddAtom2DA`, `nGetEndpointInfo` | 0 |
| `Check6MembTautRing` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:1142` | external | active, export-unreachable | `bExistsAnyAltPath`, `are_alt_bonds`, `AddBondsPos`, `AddEndPoints`, `AddAtom2num`, `AddAtom2DA`, `nGetEndpointInfo` | 0 |
| `Check15TautPathCenterpoint` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:1409` | external | active, export-unreachable | `bIsCenterPointStrict` | 0 |
| `Check15TautPath` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:1434` | external | active, export-unreachable | `bExistsAnyAltPath`, `AddBondsPos`, `AddEndPoints`, `AddAtom2num`, `AddAtom2DA`, `nGetEndpointInfo` | 0 |
| `Check5MembTautRing` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:1666` | external | active, export-unreachable | `bExistsAnyAltPath`, `are_alt_bonds`, `AddBondsPos`, `AddEndPoints`, `AddAtom2num`, `AddAtom2DA`, `nGetEndpointInfo` | 0 |
| `getInchiStateReadErr` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:536` | internal | reachable | none | 0 |
| `getInchiErrName` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:562` | internal | reachable | none | 0 |
| `SetHillFormFromInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:583` | internal | reachable | `AllocateAndFillHillFormula` | 0 |
| `ReadWriteInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:615` | external | reachable | `inchi_ios_init`, `inchi_ios_close`, `inchi_ios_reset`, `inchi_ios_eprint`, `InchiTimeGet`, `InchiTimeElapsed`, `inchi_strtol`, `SetProtonsAndXchgIsoH`, `InChILine2Data`, `PrepareSaveOptBits`, `TreatErrorsInReadInChIString`, `ConvertInChI2InChI`, `ConvertInChI2Struct`, `DetectAndExposePolymerInternals`, `FreeInpInChI` | 0 |
| `OutputInChIAsRequested` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1213` | internal | reachable | `inchi_strbuf_init`, `inchi_strbuf_close`, `CompareReversedINChI`, `bInChIHasReconnectedMetal`, `SortAndPrintINChI`, `FreeAllINChIArrays` | 0 |
| `GetNumNeighborsFromInchi` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1646` | external | reachable | none | 0 |
| `CountStereoTypes` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1723` | external | reachable | `GetNumNeighborsFromInchi` | 0 |
| `bInpInchiComponentExists` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1800` | external | reachable | none | 0 |
| `bInpInchiComponentDeleted` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1819` | external | reachable | none | 0 |
| `bRevInchiComponentExists` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1838` | external | reachable | none | 0 |
| `bRevInchiComponentDeleted` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1859` | external | reachable | none | 0 |
| `DetectInpInchiCreationOptions` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1880` | external | reachable | `CountStereoTypes`, `bInChIHasReconnectedMetal` | 0 |
| `bInChIHasReconnectedMetal` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1993` | internal | reachable | `is_el_a_metal` | 0 |
| `SetProtonsAndXchgIsoH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:2015` | internal | reachable | `nFillOutProtonMobileH`, `Free_INChI_Stereo`, `Free_INChI_Members` | 0 |
| `GetInChIFormulaNumH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:2215` | internal | reachable | `inchi_strtol` | 0 |
| `GetInChINumH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:2254` | internal | reachable | none | 0 |
| `GetInChIIsoH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:2296` | internal | reachable | none | 0 |
| `InChILine2Data` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:2342` | external | reachable | `GetInChIFormulaNumH`, `GetInChINumH`, `GetInChIIsoH`, `ReadInChICoord`, `ReadInChILine`, `nFillOutProtonMobileH`, `nProtonCopyIsotopicInfo`, `CopySegment` | 0 |
| `ParseAuxSegmentVersion` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:3783` | internal | reachable | `inchi_strtol` | 0 |
| `CopyAtomNumbers` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:3806` | internal | reachable | none | 0 |
| `ParseAuxSegmentNumbers` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:3847` | internal | reachable | `inchi_strtol`, `CopyAtomNumbers` | 0 |
| `ParseAuxSegmentAtomEqu` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:4154` | internal | reachable | none | 0 |
| `ParseAuxSegmentGroupEqu` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:4213` | internal | reachable | none | 0 |
| `ParseAuxSegmentSp3Inv` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:4252` | internal | reachable | none | 0 |
| `ParseAuxSegmentSp3InvNumbers` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:4312` | internal | reachable | none | 0 |
| `ParseAuxSegmentReverseCRV` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:4371` | internal | reachable | none | 0 |
| `ParseAuxSegmentReverseAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:4392` | internal | reachable | none | 0 |
| `ParseAuxSegmentReverseBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:4413` | internal | reachable | none | 0 |
| `ParseAuxSegmentReverseXYZ` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:4432` | internal | reachable | `inchi_strtod` | 0 |
| `AddAuxSegmentCoord` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:4491` | internal | reachable | `CopyAtomNumbers` | 0 |
| `ReadInChICoord` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:4683` | external | reachable | `ParseAuxSegmentVersion`, `ParseAuxSegmentNumbers`, `ParseAuxSegmentAtomEqu`, `ParseAuxSegmentGroupEqu`, `ParseAuxSegmentSp3Inv`, `ParseAuxSegmentSp3InvNumbers`, `ParseAuxSegmentReverseCRV`, `ParseAuxSegmentReverseAtoms`, `ParseAuxSegmentReverseBonds`, `ParseAuxSegmentReverseXYZ`, `AddAuxSegmentCoord`, `getInChIChar`, `nGetInChISegment`, `inchi_memicmp` | 0 |
| `ReadInChILine` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:5016` | internal | reachable | `ParseSegmentIsoExchgH`, `ParseSegmentPerm`, `ParseSegmentIsoAtoms`, `ParseSegmentSp3s`, `ParseSegmentSp3m`, `ParseSegmentSp3`, `ParseSegmentSp2`, `ParseSegmentProtons`, `ParseSegmentPolymer`, `ParseSegmentCharge`, `ParseSegmentMobileH`, `ParseSegmentConnections`, `ParseSegmentFormula`, `getInChIChar`, `nGetInChISegment` | 0 |
| `ParseSegmentIsoExchgH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:5541` | internal | reachable | `inchi_strtol` | 0 |
| `ParseSegmentPerm` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:5663` | internal | reachable | `inchi_strtol` | 0 |
| `ParseSegmentIsoAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:5806` | internal | reachable | `inchi_strtol`, `CopySegment` | 0 |
| `ParseSegmentSp3s` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:6235` | internal | reachable | `inchi_strtol` | 0 |
| `bIsSp3LayerNotEmpty` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:6346` | internal | reachable | none | 0 |
| `ParseSegmentSp3m` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:6378` | internal | reachable | `bIsSp3LayerNotEmpty` | 0 |
| `ParseSegmentSp3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:6600` | internal | reachable | `inchi_strtol`, `SegmentSp3CreateEmpty`, `SegmentSp3StoreStereoCenters`, `SegmentSp3CopyMultiplierCovered`, `SegmentSp3ProcessAbbreviation` | 0 |
| `ParseSegmentSp2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:6866` | internal | reachable | `inchi_strtol`, `CopySegment` | 0 |
| `ParseSegmentProtons` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:7426` | internal | reachable | `inchi_strtol` | 0 |
| `ParseSegmentPolymer` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:7494` | internal | reachable | `inchi_strtol`, `ParseSegmentReadDelimitedNumbers`, `IntArray_Alloc`, `IntArray_Append`, `IntArray_Reset`, `IntArray_Free`, `OAD_PolymerUnit_New`, `OAD_Polymer_Free`, `imat_new`, `imat_free` | 0 |
| `ParseSegmentReadDelimitedNumbers` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:7966` | external | reachable | `inchi_strtol`, `IntArray_Append` | 0 |
| `ParseSegmentCharge` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:8032` | internal | reachable | `inchi_strtol` | 0 |
| `ParseSegmentMobileH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:8230` | internal | reachable | `inchi_strtol`, `GetInChIFormulaNumH`, `GetInChINumH` | 0 |
| `ParseSegmentConnections` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:9084` | internal | reachable | `inchi_strtol`, `insertions_sort_AT_NUMB`, `AddLinkedBond` | 0 |
| `nFillOutProtonMobileH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:9718` | internal | reachable | none | 0 |
| `nProtonCopyIsotopicInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:9770` | internal | reachable | none | 0 |
| `ParseSegmentFormula` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:9810` | internal | reachable | `inchi_strtol`, `nFillOutProtonMobileH`, `get_periodic_table_number` | 0 |
| `CopySegment` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:10263` | internal | reachable | none | 0 |
| `insertions_sort_AT_NUMB` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:10546` | external | reachable | none | 0 |
| `getInChIChar` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:10565` | internal | reachable | none | 0 |
| `AddInChIChar` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:10603` | internal | reachable | `getInChIChar` | 0 |
| `nGetInChISegment` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:10676` | internal | reachable | `AddInChIChar` | 0 |
| `AddLinkedBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:10699` | internal | reachable | none | 0 |
| `PrepareSaveOptBits` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:10771` | external | reachable | `inchi_ios_eprint` | 0 |
| `TreatErrorsInReadInChIString` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:10901` | external | reachable | `inchi_ios_eprint`, `getInchiStateReadErr`, `getInchiErrName`, `FreeInpInChI` | 0 |
| `ConvertInChI2InChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:11003` | external | reachable | `InchiTimeGet`, `InchiTimeElapsed`, `SetHillFormFromInChI`, `OutputInChIAsRequested` | 0 |
| `ConvertInChI2Struct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:11203` | external | reachable | `inchi_ios_eprint`, `InchiTimeGet`, `InchiTimeElapsed`, `DetectInpInchiCreationOptions`, `RemoveFixHInChIIdentical2MobH`, `MarkDisconectedIdenticalToReconnected`, `SetUpSrm`, `MergeStructureComponents`, `AllInchiToStructure`, `AddProtonAndIsoHBalanceToMobHStruct`, `FreeStrFromINChI`, `FreeInpInChI`, `CompareAllOrigInchiToRevInChI`, `CompareAllDisconnectedOrigInchiToRevInChI`, `AddOneMsg`, `FillOutCompareMessage` | 0 |
| `DetectAndExposePolymerInternals` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:11629` | external | reachable | `inchi_ios_print`, `inchi_strbuf_close`, `inchi_strtol`, `DetectHiddenPolymerStuff`, `get_periodic_table_number` | 0 |
| `DetectHiddenPolymerStuff` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:12006` | external | reachable | none | 0 |
| `SegmentSp3CreateEmpty` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:12060` | internal | reachable | none | 0 |
| `SegmentSp3StoreStereoCenters` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:12123` | internal | reachable | `inchi_strtol` | 0 |
| `SegmentSp3CopyMultiplierCovered` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:12203` | internal | reachable | `CopySegment` | 0 |
| `SegmentSp3ProcessAbbreviation` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:12239` | internal | reachable | `CopySegment` | 0 |
| `extract_from_inchi_string` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:12416` | internal | reachable | `inchi_ios_init`, `inchi_ios_close`, `inchi_ios_print`, `InChILine2Data`, `DetectAndExposePolymerInternals` | 0 |
| `extract_stereo_info_from_inchi_string` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:12499` | external | reachable | `extract_from_inchi_string`, `FreeInpInChI` | 0 |
| `extract_all_backbone_bonds_from_inchi_string` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:12547` | external | reachable | `extract_from_inchi_string`, `FreeInpInChI` | 0 |
| `QueueCreate` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:67` | external | reachable | none | 0 |
| `QueueAdd` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:87` | external | reachable | none | 0 |
| `QueueGet` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:101` | external | reachable | none | 0 |
| `QueueGetAny` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:122` | external | reachable | none | 0 |
| `QueueDelete` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:203` | external | reachable | none | 0 |
| `QueueReinit` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:216` | external | reachable | none | 0 |
| `QueueLength` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:231` | external | reachable | none | 0 |
| `QueueWrittenLength` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:245` | external | reachable | none | 0 |
| `GetMinRingSize` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:262` | external | reachable | `QueueAdd`, `QueueGet`, `QueueLength` | 0 |
| `is_bond_in_Nmax_memb_ring` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:362` | external | reachable | `QueueAdd`, `QueueGetAny`, `QueueReinit`, `QueueWrittenLength`, `GetMinRingSize` | 0 |
| `is_atom_in_3memb_ring` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:420` | external | reachable | none | 0 |
| `clear_t_group_info` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:509` | external | reachable | none | 0 |
| `GetTgroupInfoFromInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:575` | external | reachable | `clear_t_group_info` | 0 |
| `FillOutpStructEndpointFromInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:692` | external | reachable | none | 0 |
| `cmp_charge_val` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:732` | external | active, export-unreachable | none | 0 |
| `bMayBeACationInMobileHLayer` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:749` | external | reachable | none | 0 |
| `clean_charge_val` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:792` | external | reachable | `bMayBeACationInMobileHLayer`, `get_sp_element_type`, `insertions_sort`, `if_skip_add_H` | 0 |
| `GetAtomRestoreInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:912` | external | reachable | `clean_charge_val`, `if_skip_add_H`, `get_el_valence` | 0 |
| `get_sp_element_type` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:1414` | external | reachable | none | 0 |
| `ReallocTCGroups` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:1502` | external | reachable | none | 0 |
| `RegisterTCGroup` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:1527` | external | reachable | `ReallocTCGroups` | 0 |
| `nTautEndpointEdgeCap` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:1574` | external | reachable | none | 0 |
| `BondFlowMaxcapMinorder` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:1690` | external | reachable | none | 0 |
| `AtomStcapStflow` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:1790` | external | reachable | `BondFlowMaxcapMinorder` | 0 |
| `nCountBnsSizes` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:1870` | external | reachable | `RegisterTCGroup`, `nTautEndpointEdgeCap`, `AtomStcapStflow` | 0 |
| `nAddSuperCGroups` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:2279` | external | reachable | `RegisterTCGroup` | 0 |
| `AddTGroups2TCGBnStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:2425` | external | reachable | `ConnectTwoVertices` | 0 |
| `ConnectTwoVertices` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:2627` | external | reachable | none | 0 |
| `AddRadicalToMetal` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:2706` | external | reachable | none | 0 |
| `ConnectMetalFlower` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:2743` | external | reachable | `ConnectTwoVertices`, `SetEdgeCapFlow`, `SetStCapFlow` | 0 |
| `SetEdgeCapFlow` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:2911` | external | reachable | none | 0 |
| `AddEdgeFlow` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:2923` | external | reachable | none | 0 |
| `ConnectSuperCGroup` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:3058` | external | reachable | `ConnectTwoVertices`, `AddEdgeFlow` | 0 |
| `AddStCapFlow` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:3218` | external | reachable | none | 0 |
| `SetStCapFlow` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:3231` | external | reachable | none | 0 |
| `AddCGroups2TCGBnStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:3249` | external | reachable | `AtomStcapStflow`, `ConnectTwoVertices`, `AddRadicalToMetal`, `ConnectMetalFlower`, `ConnectSuperCGroup`, `AddStCapFlow` | 0 |
| `nNumEdgesToCnVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:3824` | external | reachable | none | 0 |
| `AllocateAndInitTCGBnStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:3879` | external | reachable | `DeAllocateBnStruct`, `BondFlowMaxcapMinorder`, `AtomStcapStflow`, `nNumEdgesToCnVertex` | 0 |
| `IncrZeroBondsAndClearEndpts` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4152` | external | reachable | none | 0 |
| `IncrZeroBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4172` | external | reachable | none | 0 |
| `ClearEndpts` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4191` | external | reachable | none | 0 |
| `GetDeltaChargeFromVF` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4229` | external | reachable | none | 0 |
| `EvaluateChargeChanges` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4320` | external | reachable | `GetDeltaChargeFromVF` | 0 |
| `RunBnsTestOnce` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4507` | external | reachable | `ReInitBnStructAltPaths`, `ReInitBnData`, `BalancedNetworkSearch`, `EvaluateChargeChanges` | 0 |
| `RunBnsRestoreOnce` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4546` | external | reachable | `RunBalancedNetworkSearch`, `ReInitBnStructAltPaths`, `ReInitBnData` | 0 |
| `comp_cc_cand` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4584` | external | active, export-unreachable | none | 0 |
| `get_pVA_atom_type` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4628` | external | reachable | `is_el_a_metal`, `nNoMetalNumBonds`, `nNoMetalBondsValence`, `get_endpoint_valence` | 0 |
| `AllocEdgeList` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4690` | external | reachable | none | 0 |
| `AddToEdgeList` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4738` | external | reachable | `AllocEdgeList` | 0 |
| `RemoveFromEdgeListByIndex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4759` | external | reachable | none | 0 |
| `FindInEdgeList` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4777` | external | reachable | none | 0 |
| `RemoveFromEdgeListByValue` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4794` | external | reachable | `RemoveFromEdgeListByIndex` | 0 |
| `AllocBfsQueue` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4815` | external | reachable | `QueueCreate`, `QueueDelete`, `AllocBfsQueue` | 0 |
| `RemoveForbiddenEdgeMask` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4869` | external | reachable | none | 0 |
| `SetForbiddenEdgeMask` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4880` | external | reachable | none | 0 |
| `RemoveForbiddenBondFlowBits` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4891` | external | reachable | none | 0 |
| `GetChargeFlowerUpperEdge` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4915` | external | reachable | none | 0 |
| `MakeOneInChIOutOfStrFromINChI2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:5087` | external | reachable | `MakeOneInChIOutOfStrFromINChI`, `CopyBnsToAtom` | 0 |
| `MakeOneInChIOutOfStrFromINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:5168` | external | reachable | `Create_INChI`, `IncrZeroBondsAndClearEndpts`, `CopySt2At`, `FixUnkn0DStereoBonds`, `ReconcileAllCmlBondParities`, `free_t_group_info`, `FreeInpAtomData`, `CreateInpAtomData`, `fix_odd_things`, `bNumHeterAtomHasIsotopicH`, `SetConnectedComponentNumber`, `Free_INChI`, `Alloc_INChI`, `Free_INChI_Aux`, `Alloc_INChI_Aux` | 0 |
| `ConnectDisconnectedH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:5480` | external | reachable | none | 0 |
| `DisconnectedConnectedH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:5594` | external | reachable | none | 0 |
| `MakeInChIOutOfStrFromINChI2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:5669` | external | reachable | `inchi_strbuf_init`, `inchi_strbuf_close`, `bRevInchiComponentExists`, `IncrZeroBonds`, `ClearEndpts`, `ConnectDisconnectedH`, `CopySt2At`, `FixUnkn0DStereoBonds`, `ReconcileAllCmlBondParities`, `FreeInpAtom`, `FreeOrigAtData`, `ProcessOneStructure` | 0 |
| `OutputInChIOutOfStrFromINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:5888` | external | active, export-unreachable | `inchi_strbuf_init`, `inchi_strbuf_close`, `ClearEndpts`, `FixUnkn0DStereoBonds`, `ReconcileAllCmlBondParities`, `FreeOrigAtData`, `ProcessOneStructure`, `FreeAllINChIArrays` | 0 |
| `CopyAt2St` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:57` | external | reachable | none | 0 |
| `CopySt2At` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:79` | external | reachable | none | 0 |
| `RestoreAtomConnectionsSetStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:106` | external | reachable | `GetTgroupInfoFromInChI`, `CopyAt2St`, `AddExplicitDeletedH`, `bFindCumuleneChain`, `set_cumulene_0D_parity`, `set_atom_0D_parity`, `bCanAtomBeMiddleAllene`, `bCanAtomBeTerminalAllene`, `get_element_chemical_symbol`, `is_in_the_list` | 0 |
| `SetStereoBondTypeFor0DParity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:790` | external | reachable | `set_bond_type`, `bCanAtomBeMiddleAllene` | 0 |
| `SetStereoBondTypesFrom0DStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:844` | external | reachable | `SetStereoBondTypeFor0DParity`, `set_bond_type` | 0 |
| `CopyBnsToAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:979` | external | reachable | `BondFlowMaxcapMinorder` | 0 |
| `CheckBnsConsistency` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:1155` | external | reachable | none | 0 |
| `AddExplicitDeletedH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:1470` | external | reachable | none | 0 |
| `bFindCumuleneChain` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:1565` | external | reachable | `bCanAtomBeMiddleAllene` | 0 |
| `set_bond_type` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:1614` | external | reachable | `is_in_the_list` | 0 |
| `set_cumulene_0D_parity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:1639` | external | reachable | `bFindCumuleneChain`, `is_in_the_list` | 0 |
| `set_atom_0D_parity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:1773` | external | reachable | none | 0 |
| `MoveRadToAtomsAddCharges` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:1841` | external | reachable | `RemoveRadEndpoints`, `SetRadEndpoints`, `RunBnsRestoreOnce`, `CopyBnsToAtom` | 0 |
| `AdjustTgroupsToForbiddenEdges2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:2134` | external | reachable | `get_pVA_atom_type`, `is_centerpoint_elem`, `get_el_valence`, `is_el_a_metal`, `nNoMetalNumBonds`, `nNoMetalBondsValence`, `get_endpoint_valence` | 0 |
| `RearrangePlusMinusEdgesFlow` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:3403` | external | reachable | `RunBnsRestoreOnce`, `AllocEdgeList`, `AddToEdgeList`, `RemoveForbiddenEdgeMask` | 0 |
| `IncrementZeroOrderBondsToHeteroat` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:3520` | external | reachable | `RunBnsTestOnce`, `RunBnsRestoreOnce`, `AllocEdgeList`, `AddToEdgeList`, `RemoveForbiddenEdgeMask`, `CopyBnsToAtom`, `ForbidCarbonChargeEdges` | 0 |
| `MovePlusFromS2DiaminoCarbon` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:3736` | external | reachable | `RunBnsTestOnce`, `RunBnsRestoreOnce`, `AllocEdgeList`, `AddToEdgeList`, `RemoveForbiddenEdgeMask`, `SetForbiddenEdgeMask`, `CopyBnsToAtom` | 0 |
| `EliminateChargeSeparationOnHeteroatoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:3884` | external | reachable | `is_bond_in_Nmax_memb_ring`, `RunBnsTestOnce`, `RunBnsRestoreOnce`, `AllocEdgeList`, `AddToEdgeList`, `RemoveForbiddenEdgeMask`, `SetForbiddenEdgeMask`, `CopyBnsToAtom`, `ForbidCarbonChargeEdges` | 0 |
| `RestoreCyanoGroup` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:4554` | external | reachable | `RunBnsTestOnce`, `RunBnsRestoreOnce`, `AllocEdgeList`, `RemoveForbiddenEdgeMask`, `CopyBnsToAtom`, `ForbidCarbonChargeEdges` | 0 |
| `RestoreIsoCyanoGroup` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:4663` | external | reachable | `RunBnsTestOnce`, `RunBnsRestoreOnce`, `AllocEdgeList`, `AddToEdgeList`, `RemoveForbiddenEdgeMask`, `SetForbiddenEdgeMask`, `GetChargeFlowerUpperEdge`, `CopyBnsToAtom` | 0 |
| `FixMetal_Nminus_Ominus` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:4965` | external | reachable | `RunBnsTestOnce`, `RunBnsRestoreOnce`, `AllocEdgeList`, `AddToEdgeList`, `RemoveForbiddenEdgeMask`, `SetForbiddenEdgeMask`, `GetChargeFlowerUpperEdge`, `CopyBnsToAtom` | 0 |
| `RestoreNNNgroup` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:5115` | external | reachable | `RunBnsTestOnce`, `RunBnsRestoreOnce`, `AllocEdgeList`, `AddToEdgeList`, `RemoveFromEdgeListByValue`, `RemoveForbiddenEdgeMask`, `SetForbiddenEdgeMask`, `GetChargeFlowerUpperEdge`, `CopyBnsToAtom`, `ForbidCarbonChargeEdges` | 0 |
| `EliminateNitrogen5Val3Bonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:5810` | external | reachable | `RunBnsTestOnce`, `RunBnsRestoreOnce`, `AllocEdgeList`, `RemoveForbiddenEdgeMask`, `GetChargeFlowerUpperEdge`, `CopyBnsToAtom`, `ForbidCarbonChargeEdges` | 0 |
| `Convert_SIV_to_SVI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:6054` | external | reachable | `RunBnsTestOnce`, `RunBnsRestoreOnce`, `AllocEdgeList`, `AddToEdgeList`, `RemoveForbiddenEdgeMask`, `GetChargeFlowerUpperEdge`, `CopyBnsToAtom`, `ForbidCarbonChargeEdges` | 0 |
| `PlusFromDB_N_DB_O_to_Metal` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:6270` | external | reachable | `RunBnsTestOnce`, `RunBnsRestoreOnce`, `AllocEdgeList`, `AddToEdgeList`, `RemoveForbiddenEdgeMask`, `SetForbiddenEdgeMask`, `CopyBnsToAtom` | 0 |
| `MoveMobileHToAvoidFixedBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:6446` | external | reachable | `SetForbiddenEdges`, `MarkRingSystemsInp`, `RunBnsRestoreOnce`, `RemoveForbiddenBondFlowBits`, `CopyBnsToAtom`, `AdjustTgroupsToForbiddenEdges2` | 0 |
| `RemoveRadFromMobileHEndpoint` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:6525` | external | reachable | `RunBnsTestOnce`, `RunBnsRestoreOnce`, `CopyBnsToAtom` | 0 |
| `RemoveRadFromMobileHEndpointFixH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:6963` | external | reachable | `RunBnsTestOnce`, `RunBnsRestoreOnce`, `AllocEdgeList`, `AddToEdgeList`, `RemoveForbiddenEdgeMask`, `SetForbiddenEdgeMask`, `CopyBnsToAtom` | 0 |
| `MoveChargeToMakeCenerpoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:7658` | external | reachable | `RunBnsTestOnce`, `RunBnsRestoreOnce`, `CopyBnsToAtom`, `is_centerpoint_elem` | 0 |
| `bHas_N_V` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr3.c:92` | external | reachable | none | 0 |
| `FillTgDiffHChgFH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr3.c:113` | internal | reachable | `AllocEdgeList`, `AddToEdgeList` | 0 |
| `FixFixedHRestoredStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr3.c:333` | external | reachable | `is_bond_in_Nmax_memb_ring`, `RunBnsTestOnce`, `RunBnsRestoreOnce`, `AllocEdgeList`, `AddToEdgeList`, `FindInEdgeList`, `RemoveFromEdgeListByValue`, `RemoveForbiddenEdgeMask`, `SetForbiddenEdgeMask`, `GetChargeFlowerUpperEdge`, `MakeOneInChIOutOfStrFromINChI2`, `bHas_N_V`, `FillTgDiffHChgFH`, `FillOutExtraFixedHDataRestr`, `FillOutCMP2FHINCHI` | 0 |
| `ForbidCarbonChargeEdges` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:57` | external | reachable | `AllocEdgeList`, `AddToEdgeList` | 0 |
| `ForbidNintrogenPlus2BondsInSmallRings` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:116` | external | reachable | `AddToEdgeList` | 0 |
| `FixLessHydrogenInFormula` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:177` | external | reachable | `RunBnsTestOnce`, `RunBnsRestoreOnce`, `AllocEdgeList`, `AddToEdgeList`, `RemoveForbiddenEdgeMask`, `GetChargeFlowerUpperEdge` | 0 |
| `FixMoreHydrogenInFormula` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:398` | external | reachable | `RunBnsTestOnce`, `RunBnsRestoreOnce`, `AllocEdgeList`, `AddToEdgeList`, `RemoveForbiddenEdgeMask` | 0 |
| `FixRemoveExtraTautEndpoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:589` | external | reachable | `RunBnsTestOnce`, `RunBnsRestoreOnce`, `AllocEdgeList`, `AddToEdgeList`, `RemoveForbiddenEdgeMask`, `is_centerpoint_elem` | 0 |
| `FillOutExtraFixedHDataRestr` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:754` | external | reachable | none | 0 |
| `FillOutExtraFixedHDataInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:831` | external | reachable | `GetTgroupInfoFromInChI` | 0 |
| `FillOutCMP2FHINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:870` | external | reachable | none | 0 |
| `FillOutCMP2MHINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:1021` | external | reachable | none | 0 |
| `NormalizeAndCompare` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:1274` | external | reachable | `inchi_strbuf_init`, `inchi_strbuf_close`, `inchi_strbuf_printf`, `CompareReversedINChI2`, `MergeZzInHillFormula`, `MakeOneInChIOutOfStrFromINChI2`, `FixFixedHRestoredStructure`, `FixLessHydrogenInFormula`, `FixMoreHydrogenInFormula`, `FixRemoveExtraTautEndpoints`, `FillOutExtraFixedHDataRestr`, `FixMobileHRestoredStructure`, `FixRestoredStructureStereo`, `free_t_group_info`, `FreeInpAtomData`, `Free_INChI`, `Free_INChI_Aux` | 0 |
| `CheckAndRefixStereobonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:1678` | external | reachable | `RunBnsRestoreOnce`, `AllocEdgeList`, `AddToEdgeList`, `RemoveForbiddenEdgeMask`, `SetForbiddenEdgeMask`, `CopyBnsToAtom`, `ForbidCarbonChargeEdges` | 0 |
| `MoveChargeToRemoveCenerpoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:1860` | external | reachable | `SetForbiddenEdges`, `MarkRingSystemsInp`, `RunBnsTestOnce`, `RunBnsRestoreOnce`, `AllocEdgeList`, `AddToEdgeList`, `RemoveForbiddenEdgeMask`, `SetForbiddenEdgeMask`, `RemoveForbiddenBondFlowBits`, `CopyBnsToAtom`, `is_centerpoint_elem`, `get_endpoint_valence` | 0 |
| `MakeSingleBondsMetal2ChargedHeteroat` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:2194` | external | reachable | `RunBnsRestoreOnce`, `CopyBnsToAtom` | 0 |
| `SaltBondsToCoordBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:2373` | external | reachable | `RunBnsTestOnce`, `RunBnsRestoreOnce`, `AllocEdgeList`, `AddToEdgeList`, `RemoveForbiddenEdgeMask`, `SetForbiddenEdgeMask`, `CopyBnsToAtom`, `GetPlusMinusVertex`, `bIsMetalSalt` | 0 |
| `RunBnsRestore1` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:2740` | external | reachable | `RunBnsRestoreOnce`, `AllocEdgeList`, `RemoveForbiddenEdgeMask`, `MoveRadToAtomsAddCharges`, `RearrangePlusMinusEdgesFlow`, `IncrementZeroOrderBondsToHeteroat`, `MovePlusFromS2DiaminoCarbon`, `EliminateChargeSeparationOnHeteroatoms`, `RestoreCyanoGroup`, `RestoreIsoCyanoGroup`, `FixMetal_Nminus_Ominus`, `RestoreNNNgroup`, `EliminateNitrogen5Val3Bonds`, `Convert_SIV_to_SVI`, `PlusFromDB_N_DB_O_to_Metal`, `MoveMobileHToAvoidFixedBonds`, `RemoveRadFromMobileHEndpoint`, `RemoveRadFromMobileHEndpointFixH`, `MoveChargeToMakeCenerpoints`, `ForbidCarbonChargeEdges`, `ForbidNintrogenPlus2BondsInSmallRings`, `FillOutExtraFixedHDataInChI`, `NormalizeAndCompare`, `CheckAndRefixStereobonds`, `MoveChargeToRemoveCenerpoints`, `MakeSingleBondsMetal2ChargedHeteroat`, `SaltBondsToCoordBonds` | 0 |
| `RestoreAtomMakeBNS` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:3198` | external | reachable | `DeAllocateBnStruct`, `DeAllocateBnData`, `AllocateAndInitBnData`, `is_bond_in_Nmax_memb_ring`, `FillOutpStructEndpointFromInChI`, `GetAtomRestoreInfo`, `get_sp_element_type`, `nCountBnsSizes`, `nAddSuperCGroups`, `AddTGroups2TCGBnStruct`, `AddCGroups2TCGBnStruct`, `AllocateAndInitTCGBnStruct`, `AllocBfsQueue`, `MakeOneInChIOutOfStrFromINChI`, `CopyBnsToAtom`, `CheckBnsConsistency`, `RunBnsRestore1`, `FreeInpAtomData`, `Free_INChI`, `Free_INChI_Aux`, `is_el_a_metal` | 0 |
| `OneInChI2Atom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:3477` | external | reachable | `MakeInChIOutOfStrFromINChI2`, `RestoreAtomConnectionsSetStereo`, `SetStereoBondTypesFrom0DStereo`, `RestoreAtomMakeBNS`, `ReconcileAllCmlBondParities` | 0 |
| `MakeProtonComponent` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:3574` | external | reachable | none | 0 |
| `AddRemProtonsInRestrStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:3615` | external | reachable | `AddRemoveProtonsRestr`, `bRevInchiComponentExists`, `ConnectDisconnectedH`, `DisconnectedConnectedH`, `MakeInChIOutOfStrFromINChI2`, `MakeProtonComponent`, `FreeAllINChIArrays` | 0 |
| `AddRemIsoProtonsInRestrStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:3816` | external | reachable | `AddRemoveIsoProtonsRestr`, `DisconnectedConnectedH`, `MakeInChIOutOfStrFromINChI2`, `FreeAllINChIArrays` | 0 |
| `GetPlusMinusVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr5.c:58` | external | reachable | none | 0 |
| `bIsUnsatCarbonInASmallRing` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr5.c:92` | external | reachable | `is_bond_in_Nmax_memb_ring` | 0 |
| `FixMobileHRestoredStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr5.c:141` | external | reachable | `is_bond_in_Nmax_memb_ring`, `RunBnsTestOnce`, `RunBnsRestoreOnce`, `AllocEdgeList`, `AddToEdgeList`, `FindInEdgeList`, `RemoveForbiddenEdgeMask`, `SetForbiddenEdgeMask`, `GetChargeFlowerUpperEdge`, `MakeOneInChIOutOfStrFromINChI2`, `FillOutExtraFixedHDataRestr`, `FillOutCMP2MHINCHI`, `GetPlusMinusVertex`, `bIsUnsatCarbonInASmallRing` | 0 |
| `FixRestoredStructureStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr6.c:59` | external | reachable | `CompareIcr`, `CompareReversedINChI2`, `RunBnsTestOnce`, `RunBnsRestoreOnce`, `AllocEdgeList`, `AddToEdgeList`, `RemoveForbiddenEdgeMask`, `SetForbiddenEdgeMask`, `GetChargeFlowerUpperEdge`, `MakeOneInChIOutOfStrFromINChI2`, `FillOutExtraFixedHDataRestr`, `GetPlusMinusVertex` | 0 |
| `InChI2Atom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:101` | external | reachable | `OneInChI2Atom` | 0 |
| `RemoveFixHInChIIdentical2MobH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:256` | external | reachable | `CompareReversedINChI`, `Free_INChI_Members` | 0 |
| `MarkDisconectedIdenticalToReconnected` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:283` | external | reachable | `CompareReversedINChI` | 0 |
| `SetUpSrm` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:364` | external | reachable | none | 0 |
| `MergeStructureComponents` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:406` | external | reachable | none | 0 |
| `AllInchiToStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:1042` | external | reachable | `InchiTimeGet`, `InchiTimeElapsed`, `InChI2Atom` | 0 |
| `AddProtonAndIsoHBalanceToMobHStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:1180` | external | reachable | `bInpInchiComponentExists`, `AddRemProtonsInRestrStruct`, `AddRemIsoProtonsInRestrStruct` | 0 |
| `FreeStrFromINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:1353` | external | reachable | `free_t_group_info`, `FreeAllINChIArrays` | 0 |
| `FreeInpInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:1438` | external | reachable | `FreeExtOrigAtData`, `Free_INChI_Members` | 0 |
| `CompareAllOrigInchiToRevInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:1548` | external | reachable | `CompareOneOrigInchiToRevInChI` | 0 |
| `CompareAllDisconnectedOrigInchiToRevInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:1683` | external | reachable | `CompareReversedINChI`, `bInpInchiComponentExists`, `bInpInchiComponentDeleted`, `bRevInchiComponentExists`, `bRevInchiComponentDeleted`, `CompareTwoPairsOfInChI` | 0 |
| `CompareTwoPairsOfInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:2203` | external | reachable | `CompareReversedINChI3` | 0 |
| `CompareOneOrigInchiToRevInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:2240` | external | reachable | `CompareReversedINChI3` | 0 |
| `CompareReversedStereoINChI3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:2362` | external | reachable | none | 0 |
| `CompareReversedINChI3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:2644` | external | reachable | `CompareHillFormulasNoH`, `insertions_sort_AT_NUMB`, `CompareReversedStereoINChI3` | 0 |
| `AddOneMsg` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:3138` | external | reachable | none | 0 |
| `FillOutCompareMessage` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:3179` | external | reachable | `AddOneMsg` | 0 |
| `inchi_qsort` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:66` | external | reachable | `inchi_swap` | 8 |
| `inchi_swap` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:286` | external | reachable | none | 0 |
| `insertions_sort` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:304` | external | reachable | `inchi_swap` | 1 |
| `insertions_sort_AT_NUMBERS` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:331` | external | reachable | none | 1 |
| `insertions_sort_NeighList_AT_NUMBERS` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:355` | external | reachable | none | 0 |
| `insertions_sort_AT_RANK` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:375` | external | reachable | none | 0 |
| `insertions_sort_NeighList_AT_NUMBERS3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:396` | external | reachable | none | 0 |
| `insertions_sort_NeighListBySymmAndCanonRank` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:421` | external | reachable | none | 0 |
| `CompNeighborsAT_NUMBER` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:453` | external | active, export-unreachable | none | 0 |
| `comp_AT_RANK` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:467` | external | active, export-unreachable | none | 0 |
| `CompRank` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:475` | external | active, export-unreachable | none | 0 |
| `CompRanksOrd` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:486` | external | active, export-unreachable | none | 0 |
| `CompAtomInvariants2Only` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:502` | external | reachable | none | 0 |
| `CompAtomInvariants2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:536` | external | active, export-unreachable | `CompAtomInvariants2Only` | 0 |
| `CompChemElemLex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:551` | external | active, export-unreachable | none | 0 |
| `CompareNeighListLex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:560` | external | reachable | none | 0 |
| `CompareNeighListLexUpToMaxRank` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:582` | external | reachable | none | 0 |
| `compare_NeighLists` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:607` | external | reachable | `CompareNeighListLex` | 0 |
| `CompNeighListRanks` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:615` | external | reachable | `compare_NeighLists` | 0 |
| `CompNeighLists` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:632` | external | active, export-unreachable | `compare_NeighLists` | 0 |
| `CompNeighListsUpToMaxRank` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:644` | external | active, export-unreachable | `CompareNeighListLexUpToMaxRank` | 0 |
| `CompNeighListRanksOrd` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:657` | external | active, export-unreachable | `CompNeighListRanks` | 0 |
| `CompRanksInvOrd` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:670` | external | active, export-unreachable | none | 0 |
| `CompNeighborsRanksCountEql` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:677` | external | active, export-unreachable | none | 0 |
| `CreateNeighListFromLinearCT` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:701` | external | reachable | none | 0 |
| `CreateNeighList` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:810` | external | reachable | none | 0 |
| `FreeNeighList` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:941` | external | reachable | none | 0 |
| `BreakAllTies` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:956` | external | reachable | `DifferentiateRanks2` | 0 |
| `iisort` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:1014` | external | reachable | none | 0 |
| `comp_AT_NUMB` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:140` | external | active, export-unreachable | none | 0 |
| `get_z_coord` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:147` | internal | reachable | none | 0 |
| `len3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:221` | internal | reachable | none | 0 |
| `len2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:234` | internal | reachable | none | 0 |
| `diff3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:247` | internal | reachable | none | 0 |
| `add3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:259` | internal | reachable | none | 0 |
| `mult3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:270` | internal | reachable | none | 0 |
| `change_sign3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:292` | internal | reachable | none | 0 |
| `dot_prod3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:303` | external | reachable | none | 0 |
| `dot_prodchar3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:310` | internal | reachable | none | 0 |
| `cross_prod3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:330` | external | reachable | none | 0 |
| `triple_prod` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:347` | internal | reachable | `len3`, `dot_prod3`, `cross_prod3` | 0 |
| `CompDble` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:379` | internal | active, export-unreachable | none | 0 |
| `Get2DTetrahedralAmbiguity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:402` | internal | reachable | `inchi_swap`, `insertions_sort` | 0 |
| `triple_prod_and_min_abs_sine2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:926` | internal | reachable | `len3`, `diff3`, `cross_prod3`, `triple_prod` | 0 |
| `triple_prod_and_min_abs_sine` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1145` | internal | reachable | `triple_prod` | 0 |
| `are_3_vect_in_one_plane` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1176` | internal | reachable | `triple_prod_and_min_abs_sine` | 0 |
| `are_4at_in_one_plane` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1190` | internal | reachable | `diff3`, `triple_prod_and_min_abs_sine` | 0 |
| `triple_prod_char` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1218` | internal | reachable | `len3`, `add3`, `mult3`, `triple_prod` | 0 |
| `bInpAtomHasRequirdNeigh` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1317` | internal | reachable | `get_endpoint_valence` | 0 |
| `bCanInpAtomBeAStereoCenter` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1432` | external | reachable | `is_atom_in_3memb_ring`, `bInpAtomHasRequirdNeigh` | 0 |
| `bAtomHasValence3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1574` | external | reachable | none | 0 |
| `bCanAtomHaveAStereoBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1594` | external | reachable | none | 0 |
| `bCanAtomBeMiddleAllene` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1615` | external | reachable | none | 0 |
| `bIsSuitableHeteroInpAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1634` | internal | reachable | `get_endpoint_valence` | 0 |
| `bIsOxide` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1667` | internal | reachable | `get_endpoint_valence` | 0 |
| `bCanAtomBeTerminalAllene` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1705` | external | reachable | none | 0 |
| `GetHalfStereobond0DParity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1725` | external | reachable | none | 0 |
| `FixSb0DParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1824` | external | reachable | `len3`, `mult3`, `cross_prod3` | 0 |
| `FixUnkn0DStereoBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:2006` | external | reachable | none | 0 |
| `half_stereo_bond_parity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:2121` | internal | reachable | `get_z_coord`, `len3`, `len2`, `diff3`, `mult3`, `dot_prod3`, `cross_prod3`, `bAtomHasValence3`, `bCanAtomHaveAStereoBond`, `GetHalfStereobond0DParity` | 0 |
| `save_a_stereo_bond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:2542` | internal | reachable | none | 0 |
| `get_allowed_stereo_bond_type` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:2588` | internal | reachable | none | 0 |
| `can_be_a_stereo_bond_with_isotopic_H` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:2638` | internal | reachable | `bCanAtomHaveAStereoBond`, `bIsSuitableHeteroInpAtom`, `bIsOxide`, `get_allowed_stereo_bond_type`, `get_endpoint_valence` | 0 |
| `half_stereo_bond_action` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:2901` | internal | reachable | none | 0 |
| `set_stereo_bonds_parity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:3009` | internal | reachable | `is_bond_in_Nmax_memb_ring`, `dot_prodchar3`, `triple_prod_char`, `bCanAtomHaveAStereoBond`, `bCanAtomBeMiddleAllene`, `bIsSuitableHeteroInpAtom`, `bIsOxide`, `bCanAtomBeTerminalAllene`, `FixSb0DParities`, `half_stereo_bond_parity`, `save_a_stereo_bond`, `get_allowed_stereo_bond_type`, `half_stereo_bond_action`, `get_endpoint_valence` | 0 |
| `can_be_a_stereo_atom_with_isotopic_H` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:3690` | internal | reachable | `bCanInpAtomBeAStereoCenter` | 0 |
| `GetStereocenter0DParity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:3734` | external | reachable | `insertions_sort` | 0 |
| `set_stereo_atom_parity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:3790` | internal | reachable | `get_z_coord`, `len3`, `len2`, `diff3`, `add3`, `mult3`, `change_sign3`, `Get2DTetrahedralAmbiguity`, `triple_prod_and_min_abs_sine2`, `are_3_vect_in_one_plane`, `are_4at_in_one_plane`, `bCanInpAtomBeAStereoCenter`, `GetStereocenter0DParity` | 0 |
| `set_stereo_parity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:4398` | external | reachable | `QueueCreate`, `QueueDelete`, `can_be_a_stereo_bond_with_isotopic_H`, `set_stereo_bonds_parity`, `can_be_a_stereo_atom_with_isotopic_H`, `set_stereo_atom_parity` | 0 |
| `ReconcileAllCmlBondParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:4663` | external | reachable | `ReconcileCmlIncidentBondParities`, `is_el_a_metal` | 0 |
| `ReconcileCmlIncidentBondParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:4690` | external | reachable | `ReconcileCmlIncidentBondParities`, `get_opposite_sb_atom` | 0 |
| `get_opposite_sb_atom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:4861` | external | reachable | none | 0 |
| `is_centerpoint_elem` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:157` | external | reachable | none | 0 |
| `is_centerpoint_elem_KET` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:182` | external | reachable | none | 0 |
| `is_centerpoint_elem_strict` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:190` | external | reachable | none | 0 |
| `AddAtom2num` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:211` | external | reachable | none | 0 |
| `AddAtom2DA` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:250` | external | reachable | none | 0 |
| `AddEndPoint` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:330` | external | reachable | `AddAtom2num`, `AddAtom2DA` | 0 |
| `nGetEndpointInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:359` | external | reachable | `GetChargeType`, `get_endpoint_valence` | 0 |
| `nGetEndpointInfo_PT_22_00` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:452` | external | reachable | `GetChargeType` | 0 |
| `nGetEndpointInfo_PT_16_00` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:524` | external | reachable | `GetChargeType`, `get_endpoint_valence_KET` | 0 |
| `nGetEndpointInfo_PT_06_00` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:600` | external | reachable | `GetChargeType`, `get_endpoint_valence` | 0 |
| `nGetEndpointInfo_PT_39_00` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:677` | external | reachable | `GetChargeType` | 0 |
| `nGetEndpointInfo_PT_13_00` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:756` | external | reachable | `GetChargeType` | 0 |
| `nGetEndpointInfo_PT_18_00` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:832` | external | reachable | `GetChargeType` | 0 |
| `nGetEndpointInfo_KET` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:916` | external | reachable | `GetChargeType`, `get_endpoint_valence_KET` | 0 |
| `RegisterEndPoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:1021` | external | reachable | `ReInitBnStruct`, `AddTGroups2BnStruct`, `AddCGroups2BnStruct` | 0 |
| `SetTautomericBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:1523` | external | reachable | none | 0 |
| `GetNeutralRepsIfNeeded` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:1565` | external | reachable | none | 0 |
| `FindAccessibleEndPoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:1712` | external | reachable | `bExistsAnyAltPath`, `GetNeutralRepsIfNeeded` | 0 |
| `bCanBeACPoint` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2047` | external | reachable | none | 0 |
| `GetChargeType` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2181` | external | reachable | `bCanBeACPoint`, `get_endpoint_valence` | 0 |
| `CmpCCandidates` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2229` | external | active, export-unreachable | none | 0 |
| `RegisterCPoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2249` | external | reachable | none | 0 |
| `MarkChargeGroups` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2397` | external | reachable | `bExistsAltPath`, `GetChargeType`, `RegisterCPoints` | 0 |
| `GetSaltChargeType` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2565` | external | reachable | `get_el_valence` | 0 |
| `bDoNotMergeNonTautAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2691` | external | reachable | none | 0 |
| `GetOtherSaltChargeType` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2702` | external | reachable | `is_centerpoint_elem`, `nGetEndpointInfo` | 0 |
| `GetOtherSaltType` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2828` | external | reachable | `nGetEndpointInfo` | 0 |
| `comp_candidates` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2930` | external | active, export-unreachable | none | 0 |
| `MarkSaltChargeGroups2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2961` | external | reachable | `bIsHardRemHCandidate`, `bExistsAltPath`, `AddEndPoint`, `RegisterEndPoints`, `GetSaltChargeType`, `GetOtherSaltChargeType`, `GetOtherSaltType` | 0 |
| `MarkSaltChargeGroups` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:3483` | external | reachable | `bExistsAltPath`, `AddEndPoint`, `RegisterEndPoints`, `GetSaltChargeType`, `GetOtherSaltChargeType`, `GetOtherSaltType` | 0 |
| `MergeSaltTautGroups` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:3953` | external | reachable | `bHasAcidicHydrogen`, `bHasAcidicMinus`, `AddEndPoint`, `RegisterEndPoints`, `GetSaltChargeType`, `bDoNotMergeNonTautAtom`, `GetOtherSaltChargeType`, `GetOtherSaltType` | 0 |
| `MakeIsotopicHGroup` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:4156` | external | reachable | `bHasAcidicHydrogen`, `bHasOtherExchangableH`, `bHasAcidicMinus`, `GetSaltChargeType`, `GetOtherSaltChargeType`, `GetOtherSaltType` | 0 |
| `MarkTautomerGroups` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:4336` | external | reachable | `bIsCenterPointStrict`, `nGet14TautIn7MembAltRing`, `nGet14TautIn5MembAltRing`, `nGet12TautIn5MembAltRing`, `nGet15TautIn6MembAltRing`, `nGet15TautInAltPath`, `is_centerpoint_elem`, `is_centerpoint_elem_KET`, `AddAtom2num`, `AddAtom2DA`, `nGetEndpointInfo`, `nGetEndpointInfo_PT_22_00`, `nGetEndpointInfo_PT_16_00`, `nGetEndpointInfo_PT_06_00`, `nGetEndpointInfo_PT_39_00`, `nGetEndpointInfo_PT_13_00`, `nGetEndpointInfo_PT_18_00`, `nGetEndpointInfo_KET`, `RegisterEndPoints`, `SetTautomericBonds`, `FindAccessibleEndPoints` | 0 |
| `free_t_group_info` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:6336` | external | reachable | none | 0 |
| `make_a_copy_of_t_group_info` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:6364` | external | reachable | `free_t_group_info` | 0 |
| `set_tautomer_iso_sort_keys` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:6477` | external | reachable | none | 0 |
| `CountTautomerGroups` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:6519` | external | reachable | none | 0 |
| `CompRankTautomer` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:7008` | external | active, export-unreachable | none | 0 |
| `SortTautomerGroupsAndEndpoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:7019` | external | reachable | `insertions_sort` | 0 |
| `base26_triplet_1` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1173` | external | reachable | none | 0 |
| `base26_triplet_2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1191` | external | reachable | none | 0 |
| `base26_triplet_3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1212` | external | reachable | none | 0 |
| `base26_triplet_4` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1233` | external | reachable | none | 0 |
| `base26_dublet_for_bits_28_to_36` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1262` | external | reachable | none | 0 |
| `base26_dublet_for_bits_56_to_64` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1287` | external | reachable | none | 0 |
| `get_xtra_hash_major_hex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1307` | external | reachable | none | 0 |
| `get_xtra_hash_minor_hex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1328` | external | reachable | none | 0 |
| `GetStdINCHIKeyFromStdINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_dll.c:92` | external | export root | `GetINCHIKeyFromINCHI` | 0 |
| `GetINCHIKeyFromINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_dll.c:113` | external | export root | `base26_triplet_1`, `base26_triplet_2`, `base26_triplet_3`, `base26_triplet_4`, `base26_dublet_for_bits_28_to_36`, `base26_dublet_for_bits_56_to_64`, `get_xtra_hash_major_hex`, `get_xtra_hash_minor_hex`, `sha2_csum`, `extract_inchi_substring` | 0 |
| `CheckINCHIKey` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_dll.c:513` | external | export root | none | 0 |
| `fprint_digest` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_dll.c:601` | external | active, export-unreachable | none | 0 |
| `CreateOrigInpDataFromMolfile` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:108` | external | reachable | `AddErrorMessage`, `ReadMolfileToInpAtoms`, `FreeOrigAtData` | 0 |
| `ReadMolfileToInpAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:309` | external | reachable | `AddErrorMessage`, `MakeInpAtomsFromMolfileData`, `SetInpAtomsXYZ`, `SetExtOrigAtDataByMolfileExtInput`, `ReadMolfile`, `MolfileExtractStrucNum`, `FreeMolfileData`, `mystrncpy`, `lrtrim`, `inchi_stricmp` | 0 |
| `MakeInpAtomsFromMolfileData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:508` | external | reachable | `AddErrorMessage`, `calculate_valences`, `CreateInpAtom`, `MolfileHasNoChemStruc`, `get_periodic_table_number`, `extract_H_atoms`, `is_in_the_list`, `mystrncpy` | 0 |
| `calculate_valences` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:807` | external | reachable | `AddErrorMessage`, `detect_unusual_el_valence`, `is_el_a_metal`, `get_num_H`, `is_in_the_list`, `nBondsValToMetal` | 0 |
| `SetInpAtomsXYZ` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:975` | external | reachable | `MolfileGetXYZDimAndNormFactors` | 0 |
| `CreateInpAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1033` | external | reachable | none | 0 |
| `FreeInpAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1046` | external | reachable | none | 0 |
| `FreeInpAtomData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1058` | external | reachable | `FreeInpAtom` | 0 |
| `CreateInpAtomData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1071` | external | reachable | `CreateInpAtom`, `FreeInpAtomData` | 0 |
| `FreeCompAtomData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1090` | external | reachable | `FreeInpAtom` | 0 |
| `FreeOrigAtData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1224` | external | reachable | `FreeInpAtom`, `FreeExtOrigAtData` | 0 |
| `FreeExtOrigAtData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1269` | external | reachable | `OAD_Polymer_Free` | 0 |
| `SetExtOrigAtDataByMolfileExtInput` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1340` | external | reachable | `AddErrorMessage`, `FreeExtOrigAtData` | 0 |
| `ReadMolfile` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:92` | external | reachable | `MolfileReadDataLines`, `MolfileTreatPseudoElementAtoms`, `SDFileSkipExtraData` | 0 |
| `MolfileReadDataLines` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:165` | internal | reachable | `AddErrorMessage`, `MolfileReadHeaderLines`, `MolfileReadCountsLine`, `MolfileReadAtomsBlock`, `MolfileReadBondsBlock`, `MolfileReadSTextBlock`, `MolfileReadPropBlock`, `FreeMolfileData`, `MolfileV3000Init`, `MolfileV3000ReadCTABBeginAndCountsLine`, `MolfileV3000ReadAtomsBlock`, `MolfileV3000ReadBondsBlock`, `MolfileV3000ReadTailOfCTAB` | 0 |
| `MolfileReadHeaderLines` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:451` | internal | reachable | `inchi_fgetsLf`, `MolfileReadField`, `remove_one_lf` | 0 |
| `MolfileReadCountsLine` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:569` | internal | reachable | `inchi_fgetsLf`, `AddErrorMessage`, `MolfileReadField`, `MolFmtSgroups_Alloc`, `dotify_non_printable_chars`, `remove_one_lf` | 0 |
| `MolfileReadAtomsBlock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:661` | internal | reachable | `inchi_fgetsLf`, `AddErrorMessage`, `MolfileReadField`, `dotify_non_printable_chars`, `remove_one_lf`, `mystrncpy` | 0 |
| `MolfileReadBondsBlock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:820` | internal | reachable | `inchi_fgetsLf`, `AddErrorMessage`, `MolfileReadField`, `dotify_non_printable_chars`, `remove_one_lf` | 0 |
| `MolfileReadSTextBlock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:917` | internal | reachable | `inchi_fgetsLf`, `AddErrorMessage` | 0 |
| `MolfileReadPropBlock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:960` | internal | reachable | `inchi_fgetsLf`, `AddErrorMessage`, `MolfileReadSgroupOfPolymer`, `MolfileReadField`, `IntArray_DebugPrint`, `extract_charges_and_radicals`, `get_atomic_mass`, `normalize_string`, `dotify_non_printable_chars`, `remove_one_lf` | 0 |
| `MolfileReadSgroupOfPolymer` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:1425` | internal | reachable | `MolfileReadField`, `IntArray_Append`, `MolFmtSgroups_Append`, `MolFmtSgroups_Free`, `MolFmtSgroups_GetIndexBySgroupId`, `lrtrim` | 0 |
| `MolfileTreatPseudoElementAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:1860` | internal | reachable | `AddErrorMessage`, `mystrncpy` | 0 |
| `MolfileStrnread` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:65` | external | reachable | none | 0 |
| `MolfileReadField` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:105` | external | reachable | `MolfileStrnread` | 0 |
| `MolfileExtractStrucNum` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:333` | external | reachable | `inchi_memicmp` | 0 |
| `MolfileHasNoChemStruc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:362` | external | reachable | none | 0 |
| `MolfileSaveCopy` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:385` | external | reachable | `inchi_fgetsLf`, `mystrncpy`, `lrtrim` | 0 |
| `MolfileGetXYZDimAndNormFactors` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:474` | external | reachable | `AddErrorMessage`, `MolfileHasNoChemStruc` | 0 |
| `FreeMolfileData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:664` | external | reachable | `DeleteMolfileV3000Info`, `MolFmtSgroups_Free` | 0 |
| `MolfileV3000Init` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:68` | external | reachable | `AddErrorMessage`, `NumLists_Alloc` | 0 |
| `DeleteMolfileV3000Info` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:165` | external | reachable | `NumLists_Free` | 0 |
| `inchi_fgetsLf_V3000` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:216` | external | active, export-unreachable | `inchi_fgetsLf`, `normalize_string` | 0 |
| `MolfileV3000ReadField` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:257` | external | reachable | `read_upto_delim`, `mystrncpy` | 0 |
| `MolfileV3000ReadKeyword` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:410` | external | reachable | `read_upto_delim`, `mystrncpy` | 0 |
| `MolfileV3000ReadCTABBeginAndCountsLine` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:444` | external | reachable | `inchi_ios_init`, `inchi_strbuf_reset`, `inchi_strbuf_close`, `AddErrorMessage`, `MolfileV3000ReadField`, `get_V3000_input_line_to_strbuf`, `dotify_non_printable_chars`, `remove_one_lf` | 0 |
| `MolfileV3000ReadSGroup` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:561` | external | reachable | `inchi_ios_init`, `inchi_ios_close`, `get_V3000_input_line_to_strbuf`, `remove_one_lf` | 0 |
| `MolfileV3000Read3DBlock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:608` | external | reachable | `inchi_ios_init`, `inchi_ios_close`, `AddErrorMessage`, `get_V3000_input_line_to_strbuf`, `remove_one_lf` | 0 |
| `MolfileV3000ReadCollections` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:647` | external | reachable | `inchi_ios_init`, `inchi_ios_close`, `inchi_strbuf_reset`, `AddErrorMessage`, `MolfileV3000ReadField`, `MolfileV3000ReadKeyword`, `get_actual_atom_number`, `MolfileV3000ReadStereoCollection`, `get_V3000_input_line_to_strbuf`, `NumLists_Append`, `dotify_non_printable_chars`, `read_upto_delim`, `remove_one_lf` | 0 |
| `MolfileV3000ReadAtomsBlock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:862` | external | reachable | `inchi_ios_init`, `inchi_strbuf_reset`, `inchi_strbuf_close`, `AddErrorMessage`, `MolfileV3000ReadField`, `MolfileV3000ReadKeyword`, `get_V3000_input_line_to_strbuf`, `get_atomic_mass`, `dotify_non_printable_chars`, `remove_one_lf`, `mystrncpy` | 0 |
| `MolfileV3000ReadBondsBlock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:1262` | external | reachable | `inchi_ios_init`, `inchi_ios_close`, `inchi_strbuf_reset`, `AddErrorMessage`, `MolfileV3000ReadField`, `MolfileV3000ReadKeyword`, `get_actual_atom_number`, `MolfileV3000ReadHapticBond`, `get_V3000_input_line_to_strbuf`, `NumLists_Append`, `dotify_non_printable_chars`, `remove_one_lf` | 0 |
| `get_actual_atom_number` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:1585` | internal | reachable | none | 0 |
| `MolfileV3000ReadTailOfCTAB` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:1602` | external | reachable | `inchi_ios_init`, `inchi_strbuf_reset`, `inchi_strbuf_close`, `AddErrorMessage`, `MolfileV3000ReadSGroup`, `MolfileV3000Read3DBlock`, `MolfileV3000ReadCollections`, `get_V3000_input_line_to_strbuf`, `remove_one_lf` | 0 |
| `MolfileV3000ReadHapticBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:1732` | external | reachable | `MolfileV3000ReadField`, `read_upto_delim` | 0 |
| `MolfileV3000ReadStereoCollection` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:1817` | external | reachable | `MolfileV3000ReadField`, `read_upto_delim` | 0 |
| `get_V3000_input_line_to_strbuf` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:1889` | external | reachable | `inchi_strbuf_reset`, `inchi_strbuf_addline`, `remove_trailing_spaces` | 0 |
| `SDFileSkipExtraData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:161` | external | reachable | `inchi_fgetsLf`, `AddErrorMessage`, `SDFileIdentifyLabel`, `SDFileExtractCASNo`, `normalize_string`, `dotify_non_printable_chars`, `remove_trailing_spaces`, `mystrncpy`, `inchi_memicmp` | 0 |
| `SDFileIdentifyLabel` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:346` | external | reachable | `inchi_memicmp` | 0 |
| `SDFileExtractCASNo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:407` | external | reachable | none | 0 |
| `NumLists_Alloc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:433` | external | reachable | none | 0 |
| `NumLists_ReAlloc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:449` | external | reachable | none | 0 |
| `NumLists_Append` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:471` | external | reachable | `NumLists_ReAlloc` | 0 |
| `NumLists_Free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:491` | external | reachable | none | 0 |
| `IntArray_Alloc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:510` | external | reachable | none | 0 |
| `IntArray_ReAlloc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:525` | external | reachable | none | 0 |
| `IntArray_Append` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:550` | external | reachable | `IntArray_ReAlloc` | 0 |
| `IntArray_AppendIfAbsent` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:572` | external | reachable | `IntArray_Append`, `is_in_the_ilist` | 0 |
| `IntArray_DebugPrint` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:582` | external | reachable | none | 0 |
| `IntArray_Reset` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:603` | external | reachable | none | 0 |
| `IntArray_Free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:610` | external | reachable | none | 0 |
| `MolFmtSgroup_Create` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:634` | external | reachable | `IntArray_Alloc`, `MolFmtSgroup_Free` | 0 |
| `MolFmtSgroup_Free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:658` | external | reachable | `IntArray_Free` | 0 |
| `MolFmtSgroups_Alloc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:675` | external | reachable | none | 0 |
| `MolFmtSgroups_ReAlloc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:693` | external | reachable | none | 0 |
| `MolFmtSgroups_Append` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:715` | external | reachable | `MolFmtSgroup_Create`, `MolFmtSgroup_Free`, `MolFmtSgroups_ReAlloc` | 0 |
| `MolFmtSgroups_Free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:751` | external | reachable | `MolFmtSgroup_Free` | 0 |
| `MolFmtSgroups_GetIndexBySgroupId` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:769` | external | reachable | none | 0 |
| `OrigAtData_WriteToSDfile` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:783` | external | reachable | `inchi_ios_print_nodisplay`, `IntArray_Alloc`, `IntArray_Free`, `OrigAtData_WriteToSDfileHeaderAndCountThings`, `OrigAtData_WriteToSDfileAtomsBlock`, `OrigAtData_WriteToSDfileBondsBlock`, `OrigAtData_WriteToSDfileAdditionalLines` | 0 |
| `OrigAtData_WriteToSDfileHeaderAndCountThings` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:852` | external | reachable | `inchi_ios_print_nodisplay` | 0 |
| `OrigAtData_WriteToSDfileAtomsBlock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:989` | external | reachable | `nBondsValenceInpAt`, `inchi_ios_print_nodisplay`, `needed_unusual_el_valence` | 0 |
| `OrigAtData_WriteToSDfileBondsBlock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:1122` | external | reachable | `inchi_ios_print_nodisplay`, `IntArray_Append` | 0 |
| `OrigAtData_WriteToSDfileAdditionalLines` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:1182` | external | reachable | `inchi_ios_print_nodisplay`, `OrigAtData_WriteToSDfilePolymerData`, `get_atomic_mass_from_elnum` | 0 |
| `OrigAtData_WriteToSDfilePolymerData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:1390` | external | reachable | `inchi_ios_print_nodisplay` | 0 |
| `CreateInchi_Stereo0D` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:118` | external | reachable | none | 0 |
| `FreeInchi_Stereo0D` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:125` | external | reachable | none | 0 |
| `Extract0DParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:136` | external | reachable | `AddErrorMessage`, `FixUnkn0DStereoBonds`, `ReconcileAllCmlBondParities`, `is_in_the_list` | 0 |
| `FindToken` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:445` | external | reachable | `inchi_ios_getsTab1` | 0 |
| `LoadLine` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:489` | external | reachable | `inchi_ios_getsTab1` | 0 |
| `InchiToInpAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:549` | external | reachable | `inchi_ios_getsTab`, `inchi_ios_getsTab1`, `AddErrorMessage`, `CreateInpAtom`, `CreateInchi_Stereo0D`, `FreeInchi_Stereo0D`, `Extract0DParities`, `FindToken`, `LoadLine`, `find_and_interpret_structure_header`, `get_periodic_table_number`, `detect_unusual_el_valence`, `is_el_a_metal`, `get_num_H`, `get_atomic_mass_from_elnum`, `is_in_the_list`, `nBondsValToMetal` | 0 |
| `find_and_interpret_structure_header` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:1742` | external | reachable | `mystrncpy`, `lrtrim` | 0 |
| `ProcessOneStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:218` | external | reachable | `AddErrorMessage`, `OrigStruct_Free`, `FreeCompAtomData`, `DoOneStructureEarlyPreprocessing`, `OrigAtData_SaveMolfile`, `OrigAtData_StoreNativeInput`, `PrepareSaveOptBits`, `DisplayOrigAndResultStructuresAndComponents`, `SaveOkProcessedMolfile`, `CreateOneStructureINChI`, `OAD_Polymer_FindBackbones`, `OAD_Polymer_GetRepresentation`, `OAD_Polymer_SmartReopenCyclizedUnits`, `SortAndPrintINChI`, `bIsStructChiral`, `TreatCreateINChIWarning` | 0 |
| `DoOneStructureEarlyPreprocessing` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:493` | internal | reachable | `AddErrorMessage`, `OAD_Edit_Underivatize`, `Ring2Chain`, `TreatErrorsInReadTheStructure` | 0 |
| `OrigAtData_SaveMolfile` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:564` | external | reachable | `OrigAtData_WriteToSDfile` | 0 |
| `OrigAtData_StoreNativeInput` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:593` | internal | reachable | `AddErrorMessage`, `OrigStruct_FillOut` | 0 |
| `PrepareSaveOptBits` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:621` | internal | reachable | none | 0 |
| `DisplayOrigAndResultStructuresAndComponents` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:679` | internal | reachable | none | 0 |
| `SaveOkProcessedMolfile` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:785` | internal | reachable | `MolfileSaveCopy` | 0 |
| `CreateOneStructureINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:809` | external | reachable | `CreateCompositeNormAtom`, `InchiTimeGet`, `InchiTimeElapsed`, `AddErrorMessage`, `FreeInpAtomData`, `CreateOneComponentINChI`, `GetOneComponent`, `TreatErrorsInReadTheStructure`, `PreprocessOneStructure`, `DisplayTheWholeStructure`, `TreatErrorsInCreateOneComponentINChI` | 0 |
| `CreateOneComponentINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:1747` | external | reachable | `InchiTimeGet`, `InchiTimeElapsed`, `InchiTimeAddMsec`, `Create_INChI`, `FreeInpAtomData`, `CreateInpAtomData`, `GetProcessingWarningsOneComponentInChI`, `SetConnectedComponentNumber`, `Alloc_INChI`, `Alloc_INChI_Aux` | 0 |
| `ProcessOneStructureEx` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:2029` | external | reachable | `inchi_ios_eprint`, `PreprocessPolymerCRUData`, `ProcessOneStructureExCore`, `remove_one_lf` | 0 |
| `PreprocessPolymerCRUData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:2147` | external | reachable | `AddErrorMessage`, `OAD_StructureEdits_Apply`, `OAD_ProcessOneStructureNoEdits`, `OAD_ProcessOneStructure105Plus`, `OAD_StructureEdits_Init`, `OAD_StructureEdits_Clear`, `OAD_StructureEdits_DebugPrint`, `OAD_Polymer_PrepareFoldCRUEdits`, `OAD_Polymer_PrepareFrameShiftEdits` | 0 |
| `swap_atoms_xyz` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:2401` | internal | reachable | none | 0 |
| `OAD_StructureEdits_Apply` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:2427` | external | reachable | `swap_atoms_xyz`, `set_renumbered_or_delete`, `mark_atoms_to_delete_or_renumber`, `bIsSameBond`, `OrigAtData_DebugTrace`, `OrigAtData_RemoveBond`, `OrigAtData_AddBond`, `OAD_Polymer_DebugTrace` | 0 |
| `set_renumbered_or_delete` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:2775` | internal | reachable | none | 0 |
| `ProcessOneStructureExCore` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:2803` | external | reachable | `EditINCHI_HidePolymerZz`, `ProcessOneStructure`, `ValidateAndPreparePolymerAndPseudoatoms` | 0 |
| `ValidateAndPreparePolymerAndPseudoatoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:2879` | external | reachable | `inchi_ios_eprint`, `AddErrorMessage`, `OAD_ValidatePolymerAndPseudoElementData`, `OAD_Polymer_CyclizeCloseableUnits` | 0 |
| `OAD_ProcessOneStructureNoEdits` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:2967` | external | reachable | `ProcessOneStructureExCore`, `POSEContext_Init`, `POSEContext_Free`, `extract_inchi_substring`, `extract_auxinfo_substring` | 0 |
| `OAD_ProcessOneStructure105Plus` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:3062` | external | reachable | `ProcessOneStructureExCore`, `POSEContext_Init`, `POSEContext_Free`, `extract_inchi_substring`, `extract_auxinfo_substring` | 0 |
| `mark_atoms_to_delete_or_renumber` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:3155` | internal | reachable | `IntArray_AppendIfAbsent`, `subgraf_new`, `subgraf_free`, `subgraf_pathfinder_new`, `subgraf_pathfinder_free`, `subgraf_pathfinder_collect_all`, `is_in_the_ilist` | 0 |
| `GetOneStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:134` | external | reachable | `FreeOrigAtData`, `ReadTheStructure`, `TreatErrorsInReadTheStructure` | 0 |
| `GetOneComponent` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:345` | external | reachable | `inchi_ios_eprint`, `InchiTimeGet`, `InchiTimeElapsed`, `AddErrorMessage`, `CreateInpAtomData`, `ExtractConnectedComponent` | 0 |
| `ReadTheStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:396` | external | reachable | `InchiTimeGet`, `InchiTimeElapsed`, `CreateOrigInpDataFromMolfile`, `InchiToOrigAtom` | 0 |
| `TreatErrorsInReadTheStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:716` | external | reachable | `inchi_ios_eprint`, `GetInpStructErrorType`, `MolfileSaveCopy` | 0 |
| `InchiToOrigAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1067` | external | reachable | `AddErrorMessage`, `ReconcileAllCmlBondParities`, `FreeOrigAtData`, `InchiToInpAtom` | 0 |
| `bIsSameBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1293` | external | reachable | none | 0 |
| `GetFrameShiftInfoFrom105PlusInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1305` | internal | reachable | `inchi_strtol` | 0 |
| `extract_orig_nums_from_auxinfo_string` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1389` | external | reachable | `inchi_strtol` | 0 |
| `extract_nonstereo_eq_classes_from_auxinfo_string` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1427` | external | reachable | `inchi_strtol` | 0 |
| `POSEContext_Init` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1503` | external | reachable | `inchi_ios_init`, `inchi_strbuf_init`, `inchi_strbuf_create_copy`, `OrigAtData_Duplicate` | 0 |
| `POSEContext_Free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1636` | external | reachable | `inchi_ios_close`, `inchi_strbuf_close`, `FreeOrigAtData`, `FreeAllINChIArrays` | 0 |
| `POSEContext_DebugPrint` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1672` | external | active, export-unreachable | none | 0 |
| `OAD_StructureEdits_Init` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1735` | external | reachable | `IntArray_Alloc`, `OAD_StructureEdits_Clear` | 0 |
| `OAD_StructureEdits_Clear` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1769` | external | reachable | `IntArray_Free` | 0 |
| `OAD_StructureEdits_DebugPrint` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1807` | external | reachable | `IntArray_DebugPrint` | 0 |
| `OAD_Polymer_PrepareFoldCRUEdits` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1838` | external | reachable | `extract_stereo_info_from_inchi_string`, `extract_all_backbone_bonds_from_inchi_string`, `extract_orig_nums_from_auxinfo_string`, `extract_nonstereo_eq_classes_from_auxinfo_string`, `analyze_CRU_folding`, `OAD_ValidatePolymerAndPseudoElementData`, `Inp_Atom_GetBondType` | 0 |
| `DiylFrag_New` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2007` | internal | reachable | `inchi_strbuf_printf`, `DiylFrag_Free` | 0 |
| `DiylFrag_Free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2049` | internal | reachable | `inchi_strbuf_close` | 0 |
| `DiylFrag_MakeSignature` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2069` | internal | reachable | `inchi_strbuf_printf`, `count_colors_in_sequence` | 0 |
| `DiylFrag_Diff` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2098` | internal | reachable | none | 0 |
| `DiylFrag_DebugTrace` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2119` | internal | reachable | none | 0 |
| `analyze_CRU_folding` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2146` | internal | reachable | `IntArray_Append`, `IntArray_AppendIfAbsent`, `bIsSameBond`, `DiylFrag_New`, `DiylFrag_Free`, `DiylFrag_MakeSignature`, `DiylFrag_Diff`, `DiylFrag_DebugTrace`, `len_repeating_subsequence`, `OAD_CollectReachableAtoms`, `OAD_CollectBackboneBonds`, `OAD_PolymerUnit_DebugTrace` | 0 |
| `count_colors_in_sequence` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2463` | internal | reachable | none | 0 |
| `len_repeating_subsequence` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2489` | internal | reachable | none | 0 |
| `OAD_Polymer_PrepareFrameShiftEdits` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2523` | external | reachable | `IntArray_Append`, `bIsSameBond`, `GetFrameShiftInfoFrom105PlusInChI`, `extract_orig_nums_from_auxinfo_string`, `ModSCenter_Init`, `ModSCenter_AddTo`, `ModSCenter_DelFrom`, `ModSCenter_IsChanged`, `OAD_PolymerUnit_FindEndsAndCaps` | 0 |
| `ModSCenter_Init` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2760` | internal | reachable | `NDefStereoBonds` | 0 |
| `NDefStereoBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2775` | internal | reachable | none | 0 |
| `ModSCenter_AddTo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2803` | internal | reachable | `is_in_the_ilist` | 0 |
| `ModSCenter_DelFrom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2815` | internal | reachable | none | 0 |
| `ModSCenter_IsChanged` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2836` | internal | reachable | `iisort`, `dot_prod3`, `cross_prod3`, `NDefStereoBonds`, `is_in_the_ilist` | 0 |
| `OrigAtData_bCheckUnusualValences` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:89` | internal | reachable | `AddErrorMessage`, `detect_unusual_el_valence` | 0 |
| `OrigAtData_Duplicate` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:148` | external | reachable | `OAD_PolymerUnit_CreateCopy` | 0 |
| `PreprocessOneStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:488` | external | reachable | `AddErrorMessage`, `ReconcileAllCmlBondParities`, `OrigAtData_bCheckUnusualValences`, `OrigAtData_Duplicate`, `fix_odd_things`, `post_fix_odd_things`, `remove_ion_pairs`, `DisconnectSalts`, `bMayDisconnectMetals`, `DisconnectMetals`, `bNumHeterAtomHasIsotopicH`, `MarkDisconnectedComponents` | 0 |
| `OrigAtData_DebugTrace` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1108` | external | reachable | none | 0 |
| `OAD_PolymerUnit_New` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1162` | external | reachable | `OAD_PolymerUnit_Free` | 0 |
| `OAD_PolymerUnit_CreateCopy` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1260` | external | reachable | `OAD_PolymerUnit_Free`, `imat_new` | 0 |
| `OAD_PolymerUnit_Free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1342` | external | reachable | `OAD_PolymerUnit_DebugTrace`, `imat_free` | 0 |
| `OAD_PolymerUnit_CompareAtomListsMod` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1377` | external | reachable | none | 0 |
| `OAD_PolymerUnit_CompareAtomLists` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1400` | external | active, export-unreachable | none | 0 |
| `OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1437` | external | reachable | `is_in_the_ilist` | 0 |
| `OAD_ValidatePolymerAndPseudoElementData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1516` | external | reachable | `AddErrorMessage`, `OAD_PolymerUnit_SetEndsAndCaps`, `OAD_Polymer_GetRepresentation`, `OAD_ValidateAndSortOutPseudoElementAtoms`, `imat_new`, `imat_free`, `is_in_the_ilist`, `is_ilist_inside` | 0 |
| `UnMarkRingSystemsInp` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1961` | external | reachable | none | 0 |
| `OAD_Polymer_CyclizeCloseableUnits` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1979` | external | reachable | `AddErrorMessage`, `OAD_PolymerUnit_HasMetal`, `OAD_PolymerUnit_UnlinkCapsAndConnectEndAtoms`, `OAD_PolymerUnit_SetEndsAndCaps` | 0 |
| `OAD_PolymerUnit_HasMetal` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2055` | external | reachable | `is_el_a_metal` | 0 |
| `OAD_Polymer_Free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2071` | external | reachable | `OAD_PolymerUnit_Free` | 0 |
| `OAD_PolymerUnit_UnlinkCapsAndConnectEndAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2101` | external | reachable | `OrigAtData_RemoveBond`, `OrigAtData_AddSingleStereolessBond`, `OrigAtData_IncreaseBondOrder` | 0 |
| `OAD_PolymerUnit_FindEndsAndCaps` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2166` | external | reachable | `AddErrorMessage`, `is_in_the_ilist` | 0 |
| `OAD_PolymerUnit_SetEndsAndCaps` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2266` | external | reachable | `OAD_PolymerUnit_FindEndsAndCaps` | 0 |
| `OAD_Polymer_PrepareWorkingSet` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2342` | external | reachable | `iisort`, `OAD_PolymerUnit_CompareAtomListsMod`, `OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves` | 0 |
| `OrigAtData_RemoveHalfBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2517` | external | reachable | none | 0 |
| `OrigAtData_RemoveAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2564` | external | active, export-unreachable | none | 0 |
| `OrigAtData_RemoveBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2577` | external | reachable | `OrigAtData_RemoveHalfBond` | 0 |
| `OrigAtData_AddBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2607` | external | reachable | none | 0 |
| `OrigAtData_AddSingleStereolessBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2682` | external | reachable | `OrigAtData_AddBond` | 0 |
| `OrigAtData_IncreaseBondOrder` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2692` | external | reachable | none | 0 |
| `OrigAtData_DecreaseBondOrder` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2750` | external | reachable | none | 0 |
| `OAD_CollectFragmentBondsAndAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2805` | external | active, export-unreachable | `AddErrorMessage`, `subgraf_new`, `subgraf_free`, `subgraf_pathfinder_new`, `subgraf_pathfinder_free`, `subgraf_pathfinder_run` | 0 |
| `OAD_Polymer_FindBackbones` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2872` | external | reachable | `OAD_CollectBackboneBonds`, `OAD_PolymerUnit_DelistIntraRingBackboneBonds`, `OAD_PolymerUnit_DelistHighOrderBackboneBonds` | 0 |
| `OAD_CollectBackboneAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2946` | external | active, export-unreachable | `AddErrorMessage`, `imat_new`, `imat_free`, `subgraf_new`, `subgraf_free`, `subgraf_pathfinder_new`, `subgraf_pathfinder_free`, `subgraf_pathfinder_run` | 0 |
| `OAD_CollectReachableAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3019` | external | reachable | `subgraf_new`, `subgraf_free`, `subgraf_pathfinder_new`, `subgraf_pathfinder_free`, `subgraf_pathfinder_collect_all` | 0 |
| `OAD_CollectBackboneBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3109` | external | reachable | `AddErrorMessage`, `subgraf_new`, `subgraf_free`, `subgraf_pathfinder_new`, `subgraf_pathfinder_free`, `subgraf_pathfinder_run` | 0 |
| `OAD_PolymerUnit_DelistIntraRingBackboneBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3168` | external | reachable | `OAD_Polymer_FindRingSystems`, `OAD_PolymerUnit_RemoveLinkFromCRUChain` | 0 |
| `OAD_Polymer_FindRingSystems` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3236` | external | reachable | `MarkRingSystemsInp`, `UnMarkRingSystemsInp`, `OrigAtData_RemoveBond`, `OrigAtData_AddSingleStereolessBond` | 0 |
| `OAD_Polymer_SetAtProps` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3312` | external | reachable | `OrigAtData_RemoveBond`, `OrigAtData_AddBond`, `OAD_Polymer_FindRingSystems` | 0 |
| `OAD_PolymerUnit_DelistHighOrderBackboneBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3509` | external | reachable | `OAD_PolymerUnit_RemoveLinkFromCRUChain`, `CompAtomData_GetNumMapping` | 0 |
| `OAD_PolymerUnit_RemoveLinkFromCRUChain` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3607` | internal | reachable | none | 0 |
| `OAD_PolymerUnit_DebugTrace` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3642` | external | reachable | none | 0 |
| `OAD_Polymer_DebugTrace` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3757` | external | reachable | `OAD_PolymerUnit_DebugTrace` | 0 |
| `OAD_Polymer_GetRepresentation` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3788` | external | reachable | none | 0 |
| `OAD_Polymer_SmartReopenCyclizedUnits` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3878` | external | reachable | `OAD_Polymer_SetAtProps`, `OAD_PolymerUnit_ReopenCyclized`, `OAD_PolymerUnit_SetReopeningDetails`, `OAD_PolymerUnit_SortBackboneBondsAndSetSeniors` | 0 |
| `OAD_PolymerUnit_ReopenCyclized` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3953` | external | reachable | `OrigAtData_RemoveBond`, `OrigAtData_AddSingleStereolessBond`, `OrigAtData_DecreaseBondOrder` | 0 |
| `OAD_PolymerUnit_SetReopeningDetails` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4007` | external | reachable | none | 0 |
| `OAD_PolymerUnit_SortBackboneBondsAndSetSeniors` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4056` | external | reachable | `OAD_PolymerUnit_SortBackboneBonds`, `OAD_Polymer_IsFirstAtomRankLower` | 0 |
| `OAD_PolymerUnit_SortBackboneBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4097` | external | reachable | `OAD_Polymer_CompareBackboneBondsSeniority` | 0 |
| `OAD_Polymer_CompareBackboneBondsSeniority` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4128` | external | reachable | `OAD_Polymer_CompareRanksOfTwoAtoms`, `OAD_Polymer_IsFirstAtomRankLower` | 0 |
| `OAD_Polymer_CompareRanksOfTwoAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4209` | external | reachable | none | 0 |
| `OAD_Polymer_IsFirstAtomRankLower` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4369` | external | reachable | `OAD_Polymer_CompareRanksOfTwoAtoms` | 0 |
| `OAD_ValidateAndSortOutPseudoElementAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4396` | external | reachable | `AddErrorMessage`, `mystrncpy` | 0 |
| `Inp_Atom_GetBondType` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4507` | external | reachable | none | 0 |
| `SortAndPrintINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:102` | external | reachable | `inchi_ios_print`, `OutputINChI2`, `CompINChINonTaut2` (`qsort` callback), `CompINChITaut2` (`qsort` callback) | 0 |
| `winchi_calc_inchikey` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:421` | external | active, export-unreachable | `inchi_ios_flush`, `inchi_ios_print`, `GetINCHIKeyFromINCHI`, `extract_inchi_substring` | 0 |
| `DisplayTheWholeStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1204` | external | reachable | none | 0 |
| `DisplayStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1222` | external | active, export-unreachable | none | 0 |
| `bIsStructChiral` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1271` | external | reachable | none | 0 |
| `FreeAllINChIArrays` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1322` | external | reachable | `FreeINChIArrays` | 0 |
| `FreeINChIArrays` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1357` | external | reachable | `Free_INChI`, `Free_INChI_Aux` | 0 |
| `TreatErrorsInCreateOneComponentINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1387` | external | reachable | `inchi_ios_eprint`, `ErrMsg`, `AddErrorMessage` | 0 |
| `TreatCreateINChIWarning` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1453` | external | reachable | `inchi_ios_eprint`, `MolfileSaveCopy` | 0 |
| `GetProcessingWarningsOneComponentInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1550` | external | reachable | `GetProcessingWarningsOneInChI` | 0 |
| `GetProcessingWarningsOneInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1574` | external | reachable | `AddErrorMessage` | 0 |
| `sha2_starts` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:57` | external | reachable | none | 0 |
| `sha2_process` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:72` | internal | reachable | none | 0 |
| `sha2_update` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:206` | external | reachable | `sha2_process` | 0 |
| `sha2_finish` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:257` | external | reachable | `sha2_update` | 0 |
| `sha2_file` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:288` | external | active, export-unreachable | `sha2_starts`, `sha2_update`, `sha2_finish` | 0 |
| `sha2_csum` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:312` | external | reachable | `sha2_starts`, `sha2_update`, `sha2_finish` | 0 |
| `sha2_hmac` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:325` | external | active, export-unreachable | `sha2_starts`, `sha2_update`, `sha2_finish` | 0 |
| `sha2_self_test` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:429` | external | active, export-unreachable | none | 0 |
| `stbsp_set_separators` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:274` | external | active, export-unreachable | none | 0 |
| `stbsp__lead_sign` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:294` | internal | active, export-unreachable | none | 0 |
| `stbsp__strlen_limited` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:309` | internal | active, export-unreachable | none | 0 |
| `stbsp_vsprintfcb` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:349` | external | active, export-unreachable | `stbsp__lead_sign`, `stbsp__strlen_limited`, `stbsp__real_to_parts`, `stbsp__real_to_str` | 10 |
| `stbsp_sprintf` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:1376` | external | active, export-unreachable | `stbsp_vsprintfcb` | 0 |
| `stbsp__clamp_callback` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:1393` | internal | active, export-unreachable | none | 0 |
| `stbsp__count_clamp_callback` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:1421` | internal | active, export-unreachable | none | 0 |
| `stbsp_vsnprintf` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:1430` | external | active, export-unreachable | `stbsp_vsprintfcb`, `stbsp__clamp_callback` | 0 |
| `stbsp_snprintf` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:1460` | external | active, export-unreachable | `stbsp_vsnprintf` | 0 |
| `stbsp_vsprintf` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:1472` | external | active, export-unreachable | `stbsp_vsprintfcb` | 0 |
| `stbsp__real_to_parts` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:1491` | internal | active, export-unreachable | none | 0 |
| `stbsp__raise_to_power10` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:1638` | internal | active, export-unreachable | none | 0 |
| `stbsp__real_to_str` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:1705` | internal | active, export-unreachable | `stbsp__raise_to_power10` | 0 |
| `cmp_iso_atw_diff_component_no` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:166` | external | active, export-unreachable | none | 0 |
| `the_only_doublet_neigh` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:179` | external | reachable | none | 0 |
| `fix_non_uniform_drawn_oxoanions` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:227` | external | reachable | none | 0 |
| `fix_non_uniform_drawn_amidiniums` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:452` | external | reachable | none | 0 |
| `fix_odd_things` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:603` | external | reachable | `the_only_doublet_neigh`, `fix_non_uniform_drawn_oxoanions`, `fix_non_uniform_drawn_amidiniums`, `remove_ion_pairs`, `get_el_valence` | 0 |
| `post_fix_odd_things` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:922` | external | reachable | none | 0 |
| `nFindOneOM` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:932` | external | reachable | none | 0 |
| `remove_ion_pairs` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:1049` | external | reachable | `nFindOneOM`, `get_el_valence`, `is_in_the_list`, `num_of_H`, `ion_el_group`, `has_other_ion_neigh`, `has_other_ion_in_sphere_2`, `nNoMetalNumBonds`, `nNoMetalBondsValence`, `nNoMetalNeighIndex`, `nNoMetalOtherNeighIndex`, `nNoMetalOtherNeighIndex2` | 0 |
| `RemoveInpAtBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2046` | external | reachable | `get_opposite_sb_atom` | 0 |
| `DisconnectInpAtBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2261` | external | reachable | `RemoveInpAtBond` | 0 |
| `bIsAmmoniumSalt` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2300` | external | reachable | none | 0 |
| `DisconnectAmmoniumSalt` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2398` | external | reachable | `RemoveInpAtBond` | 0 |
| `bIsMetalSalt` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2511` | external | reachable | `get_el_valence`, `get_el_type` | 0 |
| `DisconnectMetalSalt` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2612` | external | reachable | none | 0 |
| `DisconnectSalts` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2668` | external | reachable | `bIsAmmoniumSalt`, `DisconnectAmmoniumSalt`, `bIsMetalSalt`, `DisconnectMetalSalt` | 0 |
| `bIsMetalToDisconnect` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2719` | external | reachable | `get_el_valence`, `get_el_type` | 0 |
| `bMayDisconnectMetals` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2762` | external | reachable | `bIsAmmoniumSalt`, `bIsMetalSalt`, `bIsMetalToDisconnect` | 0 |
| `DisconnectMetals` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2865` | external | reachable | `bIsAmmoniumSalt`, `bIsMetalSalt`, `bIsMetalToDisconnect`, `DisconnectOneLigand`, `move_explicit_Hcation`, `get_periodic_table_number` | 0 |
| `DisconnectOneLigand` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:3112` | external | reachable | `DisconnectInpAtBond`, `get_el_valence` | 0 |
| `dist3D` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:3245` | external | reachable | none | 0 |
| `GetMinDistDistribution` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:3264` | external | reachable | `inchi_swap` | 0 |
| `move_explicit_Hcation` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:3480` | external | reachable | `RemoveInpAtBond`, `dist3D`, `GetMinDistDistribution` | 0 |
| `add_DT_to_num_H` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:3705` | external | reachable | none | 0 |
| `remove_terminal_HDT` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:3723` | external | reachable | `insertions_sort_AT_RANK` | 0 |
| `get_iat_number` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4024` | external | reachable | none | 0 |
| `bHeteroAtomMayHaveXchgIsoH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4047` | external | reachable | `get_iat_number` | 0 |
| `bNumHeterAtomHasIsotopicH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4133` | external | reachable | `get_iat_number` | 0 |
| `cmp_components` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4247` | external | active, export-unreachable | none | 0 |
| `MarkDisconnectedComponents` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4277` | external | reachable | none | 0 |
| `ExtractConnectedComponent` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4558` | external | reachable | none | 0 |
| `SetConnectedComponentNumber` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4598` | external | reachable | none | 0 |
| `Free_INChI_Stereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4611` | external | reachable | none | 0 |
| `Alloc_INChI_Stereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4629` | external | reachable | `Free_INChI_Stereo` | 0 |
| `Free_INChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4675` | external | reachable | `Free_INChI_Members` | 0 |
| `Free_INChI_Members` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4698` | external | reachable | `Free_INChI_Stereo` | 0 |
| `Alloc_INChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4722` | external | reachable | `Alloc_INChI_Stereo`, `Free_INChI` | 0 |
| `Free_INChI_Aux` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4843` | external | reachable | none | 0 |
| `Alloc_INChI_Aux` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4884` | external | reachable | `Free_INChI_Aux` | 0 |
| `CompAtomData_GetNumMapping` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4986` | external | reachable | none | 0 |
| `imat_new` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5005` | external | reachable | `imat_free` | 0 |
| `imat_free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5036` | external | reachable | none | 0 |
| `subgraf_new` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5064` | external | reachable | `subgraf_free` | 0 |
| `subgraf_free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5155` | external | reachable | none | 0 |
| `subgraf_debug_trace` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5192` | external | active, export-unreachable | none | 0 |
| `subgraf_pathfinder_new` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5221` | external | reachable | none | 0 |
| `subgraf_pathfinder_free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5252` | external | reachable | none | 0 |
| `subgraf_pathfinder_run` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5273` | external | reachable | `bIsSameBond`, `subgraf_pathfinder_run`, `add_bond_if_unseen`, `is_in_the_ilist` | 0 |
| `add_bond_if_unseen` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5380` | external | reachable | `bIsSameBond` | 0 |
| `subgraf_pathfinder_collect_all` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5422` | external | reachable | `bIsSameBond`, `subgraf_pathfinder_collect_all`, `is_in_the_ilist` | 0 |
| `get_element_chemical_symbol` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:289` | external | reachable | none | 0 |
| `get_element_or_pseudoelement_symbol` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:316` | external | reachable | none | 0 |
| `el_number_in_internal_ref_table` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:347` | external | reachable | none | 0 |
| `get_periodic_table_number` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:364` | external | reachable | `el_number_in_internal_ref_table` | 0 |
| `if_skip_add_H` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:428` | external | reachable | none | 0 |
| `get_el_valence` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:439` | external | reachable | none | 0 |
| `get_unusual_el_valence` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:454` | external | reachable | `get_el_valence` | 0 |
| `needed_unusual_el_valence` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:518` | external | reachable | `get_element_chemical_symbol`, `if_skip_add_H`, `get_el_valence`, `get_num_H` | 0 |
| `detect_unusual_el_valence` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:620` | external | reachable | `get_el_valence` | 0 |
| `get_el_type` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:679` | external | reachable | none | 0 |
| `is_el_a_metal` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:688` | external | reachable | none | 0 |
| `extract_charges_and_radicals` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:700` | external | reachable | none | 0 |
| `extract_H_atoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:774` | external | reachable | none | 0 |
| `get_num_H` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:862` | external | reachable | `el_number_in_internal_ref_table` | 0 |
| `get_atomic_mass_from_elnum` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1007` | external | reachable | none | 0 |
| `get_atomic_mass` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1040` | external | reachable | `el_number_in_internal_ref_table` | 0 |
| `is_in_the_list` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1059` | external | reachable | none | 0 |
| `is_in_the_ilist` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1072` | external | reachable | none | 0 |
| `is_ilist_inside` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1085` | external | reachable | `is_in_the_ilist` | 0 |
| `nBondsValToMetal` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1100` | external | reachable | `is_el_a_metal` | 0 |
| `num_of_H` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1129` | external | reachable | none | 0 |
| `ion_el_group` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1151` | external | reachable | none | 0 |
| `has_other_ion_neigh` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1175` | external | reachable | `ion_el_group` | 0 |
| `has_other_ion_in_sphere_2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1201` | external | reachable | `ion_el_group` | 0 |
| `nNoMetalNumBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1253` | external | reachable | `get_el_valence`, `is_el_a_metal`, `get_endpoint_valence` | 0 |
| `nNoMetalBondsValence` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1320` | external | reachable | `get_el_valence`, `is_el_a_metal`, `get_endpoint_valence` | 0 |
| `nNoMetalNeighIndex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1386` | external | reachable | `is_el_a_metal` | 0 |
| `nNoMetalOtherNeighIndex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1405` | external | reachable | `is_el_a_metal` | 0 |
| `nNoMetalOtherNeighIndex2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1426` | external | reachable | `is_el_a_metal` | 0 |
| `get_endpoint_valence` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1509` | external | reachable | none | 0 |
| `get_endpoint_valence_KET` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1530` | external | reachable | none | 0 |
| `normalize_string` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1589` | external | reachable | none | 0 |
| `dotify_non_printable_chars` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1630` | external | reachable | none | 0 |
| `read_upto_delim` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1658` | external | reachable | `is_matching_any_delim`, `mystrncpy` | 0 |
| `is_matching_any_delim` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1710` | external | reachable | none | 0 |
| `remove_trailing_spaces` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1728` | external | reachable | none | 0 |
| `remove_one_lf` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1740` | external | reachable | none | 0 |
| `mystrncpy` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1760` | external | reachable | none | 0 |
| `lrtrim` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1804` | external | reachable | none | 0 |
| `extract_inchi_substring` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1860` | external | reachable | none | 0 |
| `extract_auxinfo_substring` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1921` | external | reachable | none | 0 |
| `inchi_memicmp` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1971` | external | reachable | none | 0 |
| `inchi_stricmp` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1995` | external | reachable | none | 0 |
| `inchi__strnset` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:2022` | external | active, export-unreachable | none | 0 |
| `inchi__strdup` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:2035` | external | reachable | none | 0 |

## Plan Input

The machine-readable companion is `target/inchi-active-call-graph-audit/official_inchi_active_call_graph.json`. Only active vendored definitions may become function Port steps. Macro behavior is attached to source-backed caller/behavior steps, inactive raw bodies are excluded, and no production API may be exposed until its complete reachable closure is ported and exact official-C behavior fixtures pass.
