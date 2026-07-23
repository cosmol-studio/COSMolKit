# Official InChI Function Inventory

## Generation Contract

- Generated: `2026-07-19`
- Upstream release: `v1.07.5`
- Upstream commit: `11a87982bb518f57ac013f0b258c283655e1ea1d`
- Parser: `tree-sitter-c` over raw vendored source definitions
- Parsed C files: `77`
- C-source function definitions: `1868`
- Production C files: `60`
- Production C-source functions: `1590`
- Production header-defined functions: `13`
- GCC-configured production function locations: `1365`
- GCC-recovered production definitions: `5`
- GCC-corrected production symbol names: `6`
- GCC subset check: `complete`
- CLI C files: `2`
- Demo C files: `15`
- C-source aggregate SHA-256: `bc3fbe28c2abe261a81011a80a116fcd76abe51c2d62702d9fcef9e104332df2`

The production classification is the complete official `libinchi` target:
`INCHI_BASE/src` plus `INCHI_API/libinchi/src`. CLI and demo functions
are inventoried because the bootstrap plan requires every vendored C source
file, but they are not selected for the reusable production Rust crate.

A function is one raw C `function_definition` node. Conditional definitions
remain visible because extraction is performed before C preprocessing.
Production header definitions are listed separately so macro-controlled and
header-implemented behavior cannot disappear from the port call graph.

## Parse Integrity

| Source | Scope | Functions | Parse errors |
|---|---|---:|---:|
| `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c` | demo | 37 | 11 |
| `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_io.c` | demo | 22 | 4 |
| `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_parms.c` | demo | 17 | 210 |
| `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichimain.c` | demo | 6 | 4 |
| `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichimain_a.c` | demo | 6 | 6 |
| `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_inchi_atom.c` | demo | 6 | 0 |
| `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_mol2atom.c` | demo | 5 | 0 |
| `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_readinch.c` | demo | 0 | 0 |
| `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_readmol.c` | demo | 17 | 27 |
| `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_readstru.c` | demo | 3 | 0 |
| `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c` | demo | 29 | 3 |
| `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/mol2inchi.c` | demo | 14 | 0 |
| `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/moreutil.c` | demo | 16 | 4 |
| `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/moreutil.c` | demo | 7 | 0 |
| `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c` | demo | 35 | 59 |
| `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ichilnct.c` | production | 5 | 4 |
| `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c` | production | 57 | 52 |
| `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c` | production | 43 | 51 |
| `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c` | production | 10 | 0 |
| `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c` | production | 11 | 7 |
| `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_main.c` | production | 1 | 1 |
| `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c` | production | 58 | 44 |
| `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_inchikey_builder.c` | production | 6 | 4 |
| `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c` | production | 85 | 61 |
| `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_read_inchi.c` | production | 4 | 2 |
| `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_read_mol.c` | production | 1 | 2 |
| `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c` | production | 15 | 8 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/bcf_s.c` | production | 6 | 0 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c` | production | 111 | 33 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c` | production | 39 | 4 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c` | production | 84 | 56 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c` | production | 22 | 5 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c` | production | 27 | 3 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichierr.c` | production | 3 | 0 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiisot.c` | production | 3 | 0 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c` | production | 12 | 0 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c` | production | 32 | 9 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c` | production | 30 | 5 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c` | production | 30 | 1 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap4.c` | production | 2 | 7 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c` | production | 77 | 80 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiparm.c` | production | 10 | 156 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c` | production | 39 | 18 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c` | production | 27 | 0 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c` | production | 26 | 11 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c` | production | 16 | 2 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c` | production | 85 | 87 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c` | production | 15 | 0 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c` | production | 54 | 38 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c` | production | 32 | 8 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr3.c` | production | 6 | 15 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c` | production | 23 | 10 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr5.c` | production | 5 | 15 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr6.c` | production | 1 | 2 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c` | production | 22 | 9 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c` | production | 29 | 0 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c` | production | 47 | 22 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c` | production | 41 | 29 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c` | production | 8 | 0 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_dll.c` | production | 8 | 8 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/inchi_gui.c` | production | 3 | 1 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c` | production | 20 | 1 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c` | production | 10 | 21 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c` | production | 7 | 2 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c` | production | 16 | 2 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c` | production | 27 | 3 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/permutation_util.c` | production | 3 | 0 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c` | production | 7 | 1 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c` | production | 21 | 1 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c` | production | 31 | 1 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c` | production | 51 | 8 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c` | production | 15 | 11 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c` | production | 9 | 0 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c` | production | 53 | 4 |
| `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c` | production | 49 | 13 |
| `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c` | cli | 43 | 6 |
| `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/ichimain.c` | cli | 20 | 17 |

Production headers with parser errors:
- `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.h`: 2
- `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.h`: 3
- `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.h`: 2
- `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.h`: 2
- `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicomn.h`: 2
- `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicomp.h`: 2
- `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichierr.h`: 3
- `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimain.h`: 2
- `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.h`: 2
- `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.h`: 2
- `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.h`: 2
- `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvrs.h`: 34
- `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.h`: 2
- `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.h`: 3
- `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitime.h`: 2
- `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/inchi_api.h`: 39
- `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/incomdef.h`: 2
- `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/inpdef.h`: 2
- `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ixa.h`: 89
- `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mode.h`: 2
- `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h`: 19
- `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.h`: 2
- `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.h`: 2

## Production C Functions

### `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ichilnct.c`

Parse errors: `4`. Function definitions: `5`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `Get_std_inchi_Input_FromAuxInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ichilnct.c:79` | external | tree-sitter+gcc-aux | `extern int Get_std_inchi_Input_FromAuxInfo (char *szInchiAuxInfo, int bDoNotAddH, InchiInpData *pInchiInp)` |
| `Get_inchi_Input_FromAuxInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ichilnct.c:89` | external | tree-sitter+gcc-aux | `extern int Get_inchi_Input_FromAuxInfo (char *szInchiAuxInfo, int bDoNotAddH, int bDiffUnkUndfStereo, InchiInpData *pInchiInp)` |
| `Free_std_inchi_Input` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ichilnct.c:182` | external | tree-sitter+gcc-aux | `extern void Free_std_inchi_Input (inchi_Input *pInp)` |
| `Free_inchi_Input` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ichilnct.c:187` | external | tree-sitter+gcc-aux | `extern void Free_inchi_Input (inchi_Input *pInp)` |
| `InchiToInchi_Input` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ichilnct.c:204` | external | tree-sitter+gcc-aux | `extern int InchiToInchi_Input (INCHI_IOSTREAM *inp_molfile, inchi_Input *orig_at_data, int bMergeAllInputStructures, int bDoNotAddH, int vABParityUnknown, INPUT_TYPE nInputType, char *pSdfLabel, char *pSdfValue, long int *lSdfId, INCHI_MODE *pInpAtomFlags, int *err, char *pStrErr)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c`

Parse errors: `52`. Function definitions: `57`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `FreeINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:153` | external | tree-sitter+gcc-aux | `extern void FreeINCHI (inchi_Output *out)` |
| `FreeStdINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:183` | external | tree-sitter+gcc-aux | `extern void FreeStdINCHI (inchi_Output *out)` |
| `FreeStructFromStdINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:195` | external | tree-sitter+gcc-aux | `extern void FreeStructFromStdINCHI (inchi_OutputStruct *out)` |
| `FreeStructFromINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:208` | external | tree-sitter+gcc-aux | `extern void FreeStructFromINCHI (inchi_OutputStruct *out)` |
| `GetStdINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:242` | external | tree-sitter+gcc-aux | `extern int GetStdINCHI (inchi_Input *inp, inchi_Output *out)` |
| `GetINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:270` | external | tree-sitter+gcc-aux | `extern int GetINCHI (inchi_Input *inp, inchi_Output *out)` |
| `input_erroneously_contains_pseudoatoms` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:293` | external | tree-sitter+gcc-aux | `extern int input_erroneously_contains_pseudoatoms (inchi_Input *inp, inchi_Output *out)` |
| `GetINCHIEx` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:325` | external | tree-sitter+gcc-aux | `extern int GetINCHIEx (inchi_InputEx *inp, inchi_Output *out)` |
| `GetINCHI1` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:345` | internal | tree-sitter+gcc-aux | `static int GetINCHI1 (inchi_InputEx *extended_input, inchi_Output *out, int enforce_std_format)` |
| `produce_generation_output` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:741` | external | tree-sitter+gcc-aux | `extern void produce_generation_output (inchi_Output *out, STRUCT_DATA *sd, INPUT_PARMS *ip, INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *out_file)` |
| `copy_corrected_log_tail` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:787` | external | tree-sitter+gcc-aux | `extern void copy_corrected_log_tail (inchi_Output *out, INCHI_IOSTREAM *log_file)` |
| `CheckINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:834` | external | tree-sitter+gcc-aux | `extern int CheckINCHI (const char *szINCHI, const const int strict)` |
| `SetNumImplicitH` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:994` | external | tree-sitter+gcc-aux | `extern void SetNumImplicitH (inp_ATOM *at, int num_atoms)` |
| `parse_options_string` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:1037` | external | tree-sitter+gcc-aux | `extern int parse_options_string (char *cmd, const char **argv, int maxargs)` |
| `SetAtomProperties` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:1139` | external | tree-sitter+gcc-aux | `extern int SetAtomProperties (inp_ATOM *at, MOL_COORD (*szCoord), inchi_Atom *ati, int a1, int *nDim, char *pStrErr, int *err)` |
| `SetBondProperties` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:1235` | external | tree-sitter+gcc-aux | `extern int SetBondProperties (inp_ATOM *at, inchi_Atom *ati, int a1, int j, int nNumAtoms, int *nNumBonds, char *pStrErr, int *err)` |
| `SetAtomAndBondProperties` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:1439` | external | tree-sitter+gcc-aux | `extern int SetAtomAndBondProperties (inp_ATOM *at, inchi_Atom *ati, int a1, int bDoNotAddH, char *pStrErr, int *err)` |
| `InpAtom0DToInchiAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:1670` | external | tree-sitter+gcc-aux | `extern int InpAtom0DToInchiAtom (inp_ATOM *at, int num_inp_atoms, AT_NUM *num_atoms, inchi_Atom **atom, AT_NUM *num_stereo0D, inchi_Stereo0D **stereo0D)` |
| `ExtractOneStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:1863` | external | tree-sitter+gcc-aux | `extern int ExtractOneStructure (STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle, inchi_InputEx *inp, INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *out_file, INCHI_IOSTREAM *prb_file, ORIG_ATOM_DATA *orig_inp_data, long int *num_inp)` |
| `GetStringLength` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:2086` | external | tree-sitter+gcc-aux | `extern int GetStringLength (char *p)` |
| `GetINCHIfromINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:2114` | external | tree-sitter+gcc-aux | `extern int GetINCHIfromINCHI (inchi_InputINCHI *inpInChI, inchi_Output *out)` |
| `GetStructFromStdINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:2461` | external | tree-sitter+gcc-aux | `extern int GetStructFromStdINCHI (inchi_InputINCHI *inpInChI, inchi_OutputStruct *outStruct)` |
| `GetStructFromINCHIEx` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:2485` | external | tree-sitter+gcc-aux | `extern int GetStructFromINCHIEx (inchi_InputINCHI *inpInChI, inchi_OutputStructEx *outStruct)` |
| `GetStructFromINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:2856` | external | tree-sitter+gcc-aux | `extern int GetStructFromINCHI (inchi_InputINCHI *inpInChI, inchi_OutputStruct *out)` |
| `FreeStructFromINCHIEx` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:2892` | external | tree-sitter+gcc-aux | `extern void FreeStructFromINCHIEx (inchi_OutputStructEx *out)` |
| `FreeInChIExtInput` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:2923` | external | tree-sitter+gcc-aux | `extern void FreeInChIExtInput (inchi_Input_Polymer *polymer, inchi_Input_V3000 *v3000)` |
| `SetExtOrigAtDataByInChIExtInput` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3021` | external | tree-sitter+gcc-aux | `extern int SetExtOrigAtDataByInChIExtInput (OAD_Polymer **ppPolymer, OAD_V3000 **ppV3000, inchi_Input_Polymer *iep, inchi_Input_V3000 *iev, int nat)` |
| `SetInChIExtInputByExtOrigAtData` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3258` | external | tree-sitter+gcc-aux | `extern int SetInChIExtInputByExtOrigAtData (OAD_Polymer *orp, OAD_V3000 *orv, inchi_Input_Polymer **iip, inchi_Input_V3000 **iiv, int nat)` |
| `cdecl_GetINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3550` | external | tree-sitter | `int cdecl_GetINCHI( inchi_Input *inp, inchi_Output *out )` |
| `cdecl_GetStdINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3557` | external | tree-sitter | `int cdecl_GetStdINCHI( inchi_Input *inp, inchi_Output *out )` |
| `cdecl_FreeINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3564` | external | tree-sitter | `void cdecl_FreeINCHI( inchi_Output *out )` |
| `cdecl_FreeStdINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3571` | external | tree-sitter | `void cdecl_FreeStdINCHI( inchi_Output *out )` |
| `cdecl_GetStringLength` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3578` | external | tree-sitter | `int cdecl_GetStringLength( char *p )` |
| `cdecl_Get_inchi_Input_FromAuxInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3585` | external | tree-sitter | `int cdecl_Get_inchi_Input_FromAuxInfo( char *szInchiAuxInfo, int bDoNotAddH, int bDiffUnkUndfStereo, InchiInpData *pInchiInp )` |
| `cdecl_Get_std_inchi_Input_FromAuxInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3598` | external | tree-sitter | `int cdecl_Get_std_inchi_Input_FromAuxInfo( char *szInchiAuxInfo, int bDoNotAddH, InchiInpData *pInchiInp )` |
| `cdecl_Free_std_inchi_Input` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3609` | external | tree-sitter | `void cdecl_Free_std_inchi_Input( inchi_Input *pInp )` |
| `cdecl_Free_inchi_Input` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3616` | external | tree-sitter | `void cdecl_Free_inchi_Input( inchi_Input *pInp )` |
| `cdecl_GetStructFromINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3623` | external | tree-sitter | `int cdecl_GetStructFromINCHI( inchi_InputINCHI *inpInChI, inchi_OutputStruct *outStruct )` |
| `cdecl_GetStructFromStdINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3631` | external | tree-sitter | `int cdecl_GetStructFromStdINCHI( inchi_InputINCHI *inpInChI, inchi_OutputStruct *outStruct )` |
| `cdecl_FreeStructFromINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3638` | external | tree-sitter | `void cdecl_FreeStructFromINCHI( inchi_OutputStruct *outStruct )` |
| `cdecl_GetINCHIfromINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3645` | external | tree-sitter | `int cdecl_GetINCHIfromINCHI( inchi_InputINCHI *inpInChI, inchi_Output *out )` |
| `cdecl_FreeStructFromStdINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3653` | external | tree-sitter | `void cdecl_FreeStructFromStdINCHI( inchi_OutputStruct *outStruct )` |
| `cdecl_CheckINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3660` | external | tree-sitter | `int cdecl_CheckINCHI( const char *szINCHI, const int strict )` |
| `pasc_GetINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3709` | external | tree-sitter | `int PASCAL pasc_GetINCHI( inchi_Input *inp, inchi_Output *out )` |
| `pasc_GetStdINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3716` | external | tree-sitter | `int PASCAL pasc_GetStdINCHI( inchi_Input *inp, inchi_Output *out )` |
| `pasc_FreeINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3723` | external | tree-sitter | `void PASCAL pasc_FreeINCHI( inchi_Output *out )` |
| `pasc_FreeStdINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3730` | external | tree-sitter | `void PASCAL pasc_FreeStdINCHI( inchi_Output *out )` |
| `pasc_GetStringLength` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3737` | external | tree-sitter | `int PASCAL pasc_GetStringLength( char *p )` |
| `pasc_Get_inchi_Input_FromAuxInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3744` | external | tree-sitter | `int PASCAL pasc_Get_inchi_Input_FromAuxInfo( char *szInchiAuxInfo, int bDoNotAddH, int bDiffUnkUndfStereo, InchiInpData *pInchiInp )` |
| `pasc_Get_std_inchi_Input_FromAuxInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3755` | external | tree-sitter | `int PASCAL pasc_Get_std_inchi_Input_FromAuxInfo( char *szInchiAuxInfo, int bDoNotAddH, InchiInpData *pInchiInp )` |
| `pasc_Free_inchi_Input` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3764` | external | tree-sitter | `void PASCAL pasc_Free_inchi_Input( inchi_Input *pInp )` |
| `pasc_Free_std_inchi_Input` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3771` | external | tree-sitter | `void PASCAL pasc_Free_std_inchi_Input( inchi_Input *pInp )` |
| `pasc_FreeStructFromINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3778` | external | tree-sitter | `void PASCAL pasc_FreeStructFromINCHI( inchi_OutputStruct *out )` |
| `pasc_FreeStructFromStdINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3785` | external | tree-sitter | `void PASCAL pasc_FreeStructFromStdINCHI( inchi_OutputStruct *out )` |
| `pasc_GetStructFromINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3792` | external | tree-sitter | `int PASCAL pasc_GetStructFromINCHI( inchi_InputINCHI *inp, inchi_OutputStruct *out )` |
| `pasc_GetStructFromStdINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3799` | external | tree-sitter | `int PASCAL pasc_GetStructFromStdINCHI( inchi_InputINCHI *inp, inchi_OutputStruct *out )` |
| `pasc_CheckINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll.c:3806` | external | tree-sitter | `int PASCAL pasc_CheckINCHI( const char *szINCHI, const int strict )` |

### `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c`

Parse errors: `51`. Function definitions: `43`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `STDINCHIGEN_Create` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:122` | external | tree-sitter+gcc-aux | `extern INCHIGEN_HANDLE STDINCHIGEN_Create (void)` |
| `INCHIGEN_Create` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:130` | external | tree-sitter+gcc-aux | `extern INCHIGEN_HANDLE INCHIGEN_Create (void)` |
| `STDINCHIGEN_Setup` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:185` | external | tree-sitter+gcc-aux | `extern int STDINCHIGEN_Setup (INCHIGEN_HANDLE _HGen, INCHIGEN_DATA *pGenData, inchi_Input *pInp)` |
| `INCHIGEN_Setup` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:260` | external | tree-sitter+gcc-aux | `extern int INCHIGEN_Setup (INCHIGEN_HANDLE _HGen, INCHIGEN_DATA *pGenData, inchi_Input *pInp)` |
| `STDINCHIGEN_DoNormalization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:435` | external | tree-sitter+gcc-aux | `extern int STDINCHIGEN_DoNormalization (INCHIGEN_HANDLE HGen, INCHIGEN_DATA *pGenData)` |
| `INCHIGEN_DoNormalization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:442` | external | tree-sitter+gcc-aux | `extern int INCHIGEN_DoNormalization (INCHIGEN_HANDLE _HGen, INCHIGEN_DATA *pGenData)` |
| `STDINCHIGEN_DoCanonicalization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:689` | external | tree-sitter+gcc-aux | `extern int STDINCHIGEN_DoCanonicalization (INCHIGEN_HANDLE HGen, INCHIGEN_DATA *pGenData)` |
| `INCHIGEN_DoCanonicalization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:697` | external | tree-sitter+gcc-aux | `extern int INCHIGEN_DoCanonicalization (INCHIGEN_HANDLE _HGen, INCHIGEN_DATA *pGenData)` |
| `STDINCHIGEN_DoSerialization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:859` | external | tree-sitter+gcc-aux | `extern int STDINCHIGEN_DoSerialization (INCHIGEN_HANDLE HGen, INCHIGEN_DATA *pGenData, inchi_Output *pResults)` |
| `INCHIGEN_DoSerialization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:869` | external | tree-sitter+gcc-aux | `extern int INCHIGEN_DoSerialization (INCHIGEN_HANDLE _HGen, INCHIGEN_DATA *pGenData, inchi_Output *pResults)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:999` | external | tree-sitter | `else if (prb_file->f && (bSortPrintINChIFlags & (FLAG_SORT_PRINT_TRANSPOS_BAS \| FLAG_SORT_PRINT_TRANSPOS_REC)) )` |
| `STDINCHIGEN_Reset` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1125` | external | tree-sitter+gcc-aux | `extern void STDINCHIGEN_Reset (INCHIGEN_HANDLE HGen, INCHIGEN_DATA *pGenData, inchi_Output *pResults)` |
| `INCHIGEN_Reset` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1134` | external | tree-sitter+gcc-aux | `extern void INCHIGEN_Reset (INCHIGEN_HANDLE _HGen, INCHIGEN_DATA *pGenData, inchi_Output *pResults)` |
| `STDINCHIGEN_Destroy` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1308` | external | tree-sitter+gcc-aux | `extern void STDINCHIGEN_Destroy (INCHIGEN_HANDLE HGen)` |
| `INCHIGEN_Destroy` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1315` | external | tree-sitter+gcc-aux | `extern void INCHIGEN_Destroy (INCHIGEN_HANDLE _HGen)` |
| `cdecl_INCHIGEN_Create` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1369` | external | tree-sitter | `INCHIGEN_HANDLE cdecl_INCHIGEN_Create( void )` |
| `cdecl_STDINCHIGEN_Create` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1376` | external | tree-sitter | `INCHIGEN_HANDLE cdecl_STDINCHIGEN_Create( void )` |
| `cdecl_INCHIGEN_Setup` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1383` | external | tree-sitter | `int cdecl_INCHIGEN_Setup( INCHIGEN_HANDLE HGen, INCHIGEN_DATA * pGenData, inchi_Input * pInp )` |
| `cdecl_STDINCHIGEN_Setup` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1390` | external | tree-sitter | `int cdecl_STDINCHIGEN_Setup( INCHIGEN_HANDLE HGen, INCHIGEN_DATA * pGenData, inchi_Input * pInp )` |
| `cdecl_INCHIGEN_DoNormalization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1397` | external | tree-sitter | `int cdecl_INCHIGEN_DoNormalization( INCHIGEN_HANDLE HGen, INCHIGEN_DATA *pGenData )` |
| `cdecl_STDINCHIGEN_DoNormalization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1404` | external | tree-sitter | `int cdecl_STDINCHIGEN_DoNormalization( INCHIGEN_HANDLE HGen, INCHIGEN_DATA *pGenData )` |
| `cdecl_INCHIGEN_DoCanonicalization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1410` | external | tree-sitter | `int cdecl_INCHIGEN_DoCanonicalization( INCHIGEN_HANDLE HGen, INCHIGEN_DATA *pGenData )` |
| `cdecl_STDINCHIGEN_DoCanonicalization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1417` | external | tree-sitter | `int cdecl_STDINCHIGEN_DoCanonicalization( INCHIGEN_HANDLE HGen, INCHIGEN_DATA *pGenData )` |
| `cdecl_INCHIGEN_DoSerialization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1424` | external | tree-sitter | `int cdecl_INCHIGEN_DoSerialization( INCHIGEN_HANDLE HGen, INCHIGEN_DATA * pGenData, inchi_Output * pResults )` |
| `cdecl_STDINCHIGEN_DoSerialization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1431` | external | tree-sitter | `int cdecl_STDINCHIGEN_DoSerialization( INCHIGEN_HANDLE HGen, INCHIGEN_DATA * pGenData, inchi_Output * pResults )` |
| `cdecl_INCHIGEN_Reset` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1438` | external | tree-sitter | `void cdecl_INCHIGEN_Reset( INCHIGEN_HANDLE HGen, INCHIGEN_DATA *pGenData, inchi_Output *pResults )` |
| `cdecl_STDINCHIGEN_Reset` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1445` | external | tree-sitter | `void cdecl_STDINCHIGEN_Reset( INCHIGEN_HANDLE HGen, INCHIGEN_DATA *pGenData, inchi_Output *pResults )` |
| `cdecl_INCHIGEN_Destroy` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1452` | external | tree-sitter | `void cdecl_INCHIGEN_Destroy( INCHIGEN_HANDLE HGen )` |
| `cdecl_STDINCHIGEN_Destroy` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1459` | external | tree-sitter | `void cdecl_STDINCHIGEN_Destroy( INCHIGEN_HANDLE HGen )` |
| `pasc_INCHIGEN_Create` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1502` | external | tree-sitter | `INCHIGEN_HANDLE PASCAL pasc_INCHIGEN_Create( void )` |
| `pasc_STDINCHIGEN_Create` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1506` | external | tree-sitter | `INCHIGEN_HANDLE PASCAL pasc_STDINCHIGEN_Create( void )` |
| `pasc_INCHIGEN_Setup` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1510` | external | tree-sitter | `int PASCAL pasc_INCHIGEN_Setup( INCHIGEN_HANDLE HGen, INCHIGEN_DATA * pGenData, inchi_Input * pInp )` |
| `pasc_STDINCHIGEN_Setup` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1514` | external | tree-sitter | `int PASCAL pasc_STDINCHIGEN_Setup( INCHIGEN_HANDLE HGen, INCHIGEN_DATA * pGenData, inchi_Input * pInp )` |
| `pasc_INCHIGEN_DoNormalization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1518` | external | tree-sitter | `int PASCAL pasc_INCHIGEN_DoNormalization( INCHIGEN_HANDLE HGen, INCHIGEN_DATA *pGenData )` |
| `pasc_STDINCHIGEN_DoNormalization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1522` | external | tree-sitter | `int PASCAL pasc_STDINCHIGEN_DoNormalization( INCHIGEN_HANDLE HGen, INCHIGEN_DATA *pGenData )` |
| `pasc_INCHIGEN_DoCanonicalization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1526` | external | tree-sitter | `int PASCAL pasc_INCHIGEN_DoCanonicalization( INCHIGEN_HANDLE HGen, INCHIGEN_DATA *pGenData )` |
| `pasc_STDINCHIGEN_DoCanonicalization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1530` | external | tree-sitter | `int PASCAL pasc_STDINCHIGEN_DoCanonicalization( INCHIGEN_HANDLE HGen, INCHIGEN_DATA *pGenData )` |
| `pasc_INCHIGEN_DoSerialization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1534` | external | tree-sitter | `int PASCAL pasc_INCHIGEN_DoSerialization( INCHIGEN_HANDLE HGen, INCHIGEN_DATA * pGenData, inchi_Output * pResults )` |
| `pasc_STDINCHIGEN_DoSerialization` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1538` | external | tree-sitter | `int PASCAL pasc_STDINCHIGEN_DoSerialization( INCHIGEN_HANDLE HGen, INCHIGEN_DATA * pGenData, inchi_Output * pResults )` |
| `pasc_INCHIGEN_Reset` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1542` | external | tree-sitter | `void PASCAL pasc_INCHIGEN_Reset( INCHIGEN_HANDLE HGen, INCHIGEN_DATA *pGenData, inchi_Output *pResults )` |
| `pasc_STDINCHIGEN_Reset` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1546` | external | tree-sitter | `void PASCAL pasc_STDINCHIGEN_Reset( INCHIGEN_HANDLE HGen, INCHIGEN_DATA *pGenData, inchi_Output *pResults )` |
| `pasc_INCHIGEN_Destroy` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1550` | external | tree-sitter | `void PASCAL pasc_INCHIGEN_Destroy( INCHIGEN_HANDLE HGen )` |
| `pasc_STDINCHIGEN_Destroy` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a.c:1554` | external | tree-sitter | `void PASCAL pasc_STDINCHIGEN_Destroy( INCHIGEN_HANDLE HGen )` |

### `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c`

Parse errors: `0`. Function definitions: `10`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `NormOneStructureINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:188` | external | tree-sitter+gcc-aux | `extern int NormOneStructureINChI (CANON_GLOBALS *pCG, INCHI_CLOCK *ic, INCHIGEN_DATA *gendata, INCHIGEN_CONTROL *genctl, int iINChI, INCHI_IOSTREAM *inp_file)` |
| `CanonOneStructureINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:503` | external | tree-sitter+gcc-aux | `extern int CanonOneStructureINChI (CANON_GLOBALS *pCG, INCHI_CLOCK *ic, INCHIGEN_CONTROL *genctl, int iINChI, INCHI_IOSTREAM *inp_file)` |
| `NormOneComponentINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:654` | external | tree-sitter+gcc-aux | `extern int NormOneComponentINChI (CANON_GLOBALS *pCG, INCHI_CLOCK *ic, INCHIGEN_CONTROL *genctl, int iINChI, int i)` |
| `CanonOneComponentINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:923` | external | tree-sitter+gcc-aux | `extern int CanonOneComponentINChI (CANON_GLOBALS *pCG, INCHI_CLOCK *ic, INCHIGEN_CONTROL *genctl, int iINChI, int i)` |
| `Normalization_step` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:1138` | external | tree-sitter+gcc-aux | `extern int Normalization_step (CANON_GLOBALS *pCG, INCHI_CLOCK *ic, INChI **ppINChI, INChI_Aux **ppINChI_Aux, inp_ATOM *inp_at, INP_ATOM_DATA **out_norm_data, int num_inp_at, struct tagInchiTime *ulMaxTime, INCHI_MODE *pbTautFlags, INCHI_MODE *pbTautFlagsDone, COMPONENT_TREAT_INFO *z)` |
| `Canonicalization_step` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:1612` | external | tree-sitter+gcc-aux | `extern int Canonicalization_step (CANON_GLOBALS *pCG, INCHI_CLOCK *ic, INChI **ppINChI, INChI_Aux **ppINChI_Aux, INP_ATOM_DATA **out_norm_data, struct tagInchiTime *ulMaxTime, T_GROUP_INFO *ti_out, char *pStrErrStruct, COMPONENT_TREAT_INFO *z, int LargeMolecules)` |
| `CreateCompositeNormAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:1973` | external | tree-sitter+gcc-aux | `extern int CreateCompositeNormAtom (COMP_ATOM_DATA *composite_norm_data, INP_ATOM_DATA2 (*all_inp_norm_data), int num_components)` |
| `CreateCompAtomData` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:2232` | external | tree-sitter+gcc-aux | `extern int CreateCompAtomData (COMP_ATOM_DATA *inp_at_data, int num_atoms, int num_components, int bIntermediateTaut)` |
| `FillOutINChIReducedWarn` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:2254` | external | tree-sitter+gcc-aux | `extern int FillOutINChIReducedWarn (INChI *pINChI, INChI_Aux *pINChI_Aux, int num_atoms, int num_at_tg, int num_removed_H, sp_ATOM *at, inp_ATOM *norm_at, CANON_STAT *pCS, CANON_GLOBALS *pCG, int bTautomeric, INCHI_MODE nUserMode, char *pStrErrStruct)` |
| `make_norm_atoms_from_inp_atoms` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_a2.c:2937` | external | tree-sitter+gcc-aux | `extern void make_norm_atoms_from_inp_atoms (INCHIGEN_DATA *gendata, INCHIGEN_CONTROL *genctl)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c`

Parse errors: `7`. Function definitions: `11`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `FreeInchi_Atom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:106` | external | tree-sitter+gcc-aux | `extern void FreeInchi_Atom (inchi_Atom **at)` |
| `CreateInchiAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:117` | external | tree-sitter+gcc-aux | `extern inchi_Atom *CreateInchiAtom (int num_atoms)` |
| `MakeINCHIFromMolfileText` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:127` | external | tree-sitter+gcc-aux | `extern int MakeINCHIFromMolfileText (const char *moltext, char *szOptions, inchi_Output *result)` |
| `PrepareToMakeINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:285` | internal | tree-sitter+gcc-aux | `static int PrepareToMakeINCHI (STRUCT_DATA *sd, INPUT_PARMS *ip, ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data, PINChI2 (**pINChI), PINChI_Aux2 (**pINChI_Aux), INCHI_IOSTREAM *pout, INCHI_IOSTREAM *plog, INCHI_IOSTREAM *pprb, INCHI_IOSTREAM *inp_file, const char *moltext, char *options, INCHI_IOS_STRING *strbuf)` |
| `PostMakeINCHICleanup` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:391` | internal | tree-sitter+gcc-aux | `static int PostMakeINCHICleanup (struct tagCANON_GLOBALS *pCG, STRUCT_DATA *sd, INPUT_PARMS *ip, ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data, PINChI2 (**pINChI), PINChI_Aux2 (**pINChI_Aux), INCHI_IOSTREAM *pout, INCHI_IOSTREAM *plog, INCHI_IOSTREAM *pprb, INCHI_IOSTREAM *inp_file, const char *moltext, INCHI_IOS_STRING *strbuf)` |
| `FreeInchi_Input` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:439` | external | tree-sitter+gcc-aux | `extern void FreeInchi_Input (inchi_Input *inp_at_data)` |
| `cdecl_MakeINCHIFromMolfileText` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:477` | external | tree-sitter | `int cdecl_MakeINCHIFromMolfileText( const char *moltext, char *options, inchi_Output *out )` |
| `pasc_MakeINCHIFromMolfileText` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:516` | external | tree-sitter | `int PASCAL pasc_MakeINCHIFromMolfileText( const char *moltext, char *options, inchi_Output *out )` |
| `is_in_the_slist` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:531` | external | tree-sitter+gcc-aux | `extern S_SHORT *is_in_the_slist (S_SHORT *pathAtom, S_SHORT nNextAtom, int nPathLen)` |
| `is_element_a_metal` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:543` | external | tree-sitter+gcc-aux | `extern int is_element_a_metal (char *szEl)` |
| `InchiToInchiAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_b.c:582` | external | tree-sitter+gcc-aux | `extern int InchiToInchiAtom (INCHI_IOSTREAM *inp_file, inchi_Stereo0D **stereo0D, int *num_stereo0D, int bDoNotAddH, int vABParityUnknown, INPUT_TYPE nInputType, inchi_Atom **at, int max_num_at, int *num_dimensions, int *num_bonds, char *pSdfLabel, char *pSdfValue, long int *Id, INCHI_MODE *pInpAtomFlags, int *err, char *pStrErr)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_main.c`

Parse errors: `1`. Function definitions: `1`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `DllMain` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/inchi_dll_main.c:52` | external | tree-sitter | `int INCHI_DLLMAIN_TYPE DllMain( HANDLE hModule, DWORD ul_reason_for_call, LPVOID lpReserved )` |

### `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c`

Parse errors: `44`. Function definitions: `58`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `GetSingleStereoCode` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:145` | internal | tree-sitter+gcc-aux | `static int GetSingleStereoCode (IXA_STATUS_HANDLE hStatus, IXA_BOND_WEDGE direction, IXA_BOND_WEDGE reverse_direction)` |
| `GetDoubleStereoCode` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:181` | internal | tree-sitter+gcc-aux | `static int GetDoubleStereoCode (IXA_STATUS_HANDLE hStatus, IXA_DBLBOND_CONFIG vConfig)` |
| `BUILDER_Unpack` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:198` | internal | tree-sitter+gcc-aux | `static INCHIBUILDER *BUILDER_Unpack (IXA_STATUS_HANDLE hStatus, IXA_INCHIBUILDER_HANDLE hBuilder)` |
| `BUILDER_Pack` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:212` | internal | tree-sitter+gcc-aux | `static IXA_INCHIBUILDER_HANDLE BUILDER_Pack (INCHIBUILDER *pBuilder)` |
| `TranslateTetrahedralVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:219` | internal | tree-sitter+gcc-aux | `static void TranslateTetrahedralVertex (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_STEREOID vSourceStereo, inchi_Stereo0D *pTargetStereo, int vVertexIndex)` |
| `ExtendAllene` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:246` | internal | tree-sitter+gcc-aux | `static void ExtendAllene (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vCentralAtom, IXA_ATOMID *pAtom1, IXA_ATOMID *pAtom2)` |
| `ExtendCumulene` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:313` | internal | tree-sitter+gcc-aux | `static IXA_ATOMID ExtendCumulene (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_BONDID vBond, IXA_ATOMID vAtom)` |
| `IsRectangularVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:361` | internal | tree-sitter+gcc-aux | `static IXA_BOOL IsRectangularVertex (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vVertex, IXA_ATOMID vInternal)` |
| `IsRectOrAntiRectCentre` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:392` | internal | tree-sitter+gcc-aux | `static IXA_BOOL IsRectOrAntiRectCentre (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vVertex1a, IXA_ATOMID vVertex1b, IXA_ATOMID vInternal1, IXA_ATOMID vVertex2a, IXA_ATOMID vVertex2b, IXA_ATOMID vInternal2)` |
| `ClearMolecule` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:422` | internal | tree-sitter+gcc-aux | `static void ClearMolecule (BuilderMolecule *pMolecule)` |
| `AppendOption` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:441` | internal | tree-sitter+gcc-aux | `static void AppendOption (char *pString, const char *pOption)` |
| `BUILDER_ClearOptions` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:453` | internal | tree-sitter+gcc-aux | `static void BUILDER_ClearOptions (INCHIBUILDER *pBuilder)` |
| `BUILDER_Update` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:513` | internal | tree-sitter+gcc-aux | `static void BUILDER_Update (IXA_STATUS_HANDLE hStatus, INCHIBUILDER *pBuilder)` |
| `IXA_INCHIBUILDER_Create` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:802` | external | tree-sitter+gcc-aux | `extern IXA_INCHIBUILDER_HANDLE IXA_INCHIBUILDER_Create (IXA_STATUS_HANDLE hStatus)` |
| `IXA_INCHIBUILDER_Destroy` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:823` | external | tree-sitter+gcc-aux | `extern void IXA_INCHIBUILDER_Destroy (IXA_STATUS_HANDLE hStatus, IXA_INCHIBUILDER_HANDLE hBuilder)` |
| `IXA_INCHIBUILDER_SetMolecule` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:841` | external | tree-sitter+gcc-aux | `extern void IXA_INCHIBUILDER_SetMolecule (IXA_STATUS_HANDLE hStatus, IXA_INCHIBUILDER_HANDLE hBuilder, IXA_MOL_HANDLE hMolecule)` |
| `IXA_INCHIBUILDER_SetOption_Stereo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1583` | external | tree-sitter+gcc-aux | `extern void IXA_INCHIBUILDER_SetOption_Stereo (IXA_STATUS_HANDLE hStatus, IXA_INCHIBUILDER_HANDLE hBuilder, IXA_INCHIBUILDER_STEREOOPTION vValue)` |
| `IXA_INCHIBUILDER_SetOption_Timeout` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1597` | external | tree-sitter+gcc-aux | `extern void IXA_INCHIBUILDER_SetOption_Timeout (IXA_STATUS_HANDLE hStatus, IXA_INCHIBUILDER_HANDLE hBuilder, int vValue)` |
| `IXA_INCHIBUILDER_SetOption_Timeout_MilliSeconds` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1612` | external | tree-sitter+gcc-aux | `extern void IXA_INCHIBUILDER_SetOption_Timeout_MilliSeconds (IXA_STATUS_HANDLE hStatus, IXA_INCHIBUILDER_HANDLE hBuilder, long int vValue)` |
| `IXA_INCHIBUILDER_SetOption` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1626` | external | tree-sitter+gcc-aux | `extern void IXA_INCHIBUILDER_SetOption (IXA_STATUS_HANDLE hStatus, IXA_INCHIBUILDER_HANDLE hBuilder, IXA_INCHIBUILDER_OPTION vOption, IXA_BOOL vValue)` |
| `IXA_INCHIBUILDER_CheckOption` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1776` | external | tree-sitter+gcc-aux | `extern IXA_BOOL IXA_INCHIBUILDER_CheckOption (IXA_STATUS_HANDLE hStatus, IXA_INCHIBUILDER_HANDLE hBuilder, IXA_INCHIBUILDER_OPTION vOption)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1823` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_PT_22_00)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1829` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_PT_16_00)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1835` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_PT_06_00)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1841` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_PT_39_00)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1847` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_PT_13_00)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1853` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_PT_18_00)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1858` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_SaveOpt)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1862` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_AuxNone)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1866` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_WarnOnEmptyStructure)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1870` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_LargeMolecules)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1874` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_Polymers)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1878` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_Polymers105)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1882` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_NPZZ)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1886` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_SATZZ)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1890` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_NoFrameShift)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1894` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_FoldCRU)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1898` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_NoEdits)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1902` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_LooseTSACheck)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1906` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_OutErrInChI)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1910` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_NoWarnings)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1916` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_DoDrv)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1920` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_DoDrvReport)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1924` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_DoR2C)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1928` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_DoneOnly)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1932` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_OnlyRecSalt)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1936` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_OnlyExact)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1940` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_OnlyRecMet)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1944` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_Polymers105)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1948` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_FilterSS)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1952` | external | tree-sitter | `else if (vOption == IXA_INCHIBUILDER_OPTION_InvFilterSS)` |
| `IXA_INCHIBUILDER_CheckOption_Stereo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1971` | external | tree-sitter+gcc-aux | `extern IXA_BOOL IXA_INCHIBUILDER_CheckOption_Stereo (IXA_STATUS_HANDLE hStatus, IXA_INCHIBUILDER_HANDLE hBuilder, IXA_INCHIBUILDER_STEREOOPTION vValue)` |
| `IXA_INCHIBUILDER_GetOption_Timeout_MilliSeconds` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:1993` | external | tree-sitter+gcc-aux | `extern long int IXA_INCHIBUILDER_GetOption_Timeout_MilliSeconds (IXA_STATUS_HANDLE hStatus, IXA_INCHIBUILDER_HANDLE hBuilder)` |
| `IXA_INCHIBUILDER_GetInChIVersion` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:2008` | external | tree-sitter+gcc-aux | `extern const char *IXA_INCHIBUILDER_GetInChIVersion (IXA_BOOL vFullDescription)` |
| `IXA_INCHIBUILDER_GetInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:2022` | external | tree-sitter+gcc-aux | `extern const char *IXA_INCHIBUILDER_GetInChI (IXA_STATUS_HANDLE hStatus, IXA_INCHIBUILDER_HANDLE hBuilder)` |
| `IXA_INCHIBUILDER_GetInChIEx` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:2037` | external | tree-sitter+gcc-aux | `extern const char *IXA_INCHIBUILDER_GetInChIEx (IXA_STATUS_HANDLE hStatus, IXA_INCHIBUILDER_HANDLE hBuilder)` |
| `IXA_INCHIBUILDER_GetAuxInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:2052` | external | tree-sitter+gcc-aux | `extern const char *IXA_INCHIBUILDER_GetAuxInfo (IXA_STATUS_HANDLE hStatus, IXA_INCHIBUILDER_HANDLE hBuilder)` |
| `IXA_INCHIBUILDER_GetLog` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_builder.c:2067` | external | tree-sitter+gcc-aux | `extern const char *IXA_INCHIBUILDER_GetLog (IXA_STATUS_HANDLE hStatus, IXA_INCHIBUILDER_HANDLE hBuilder)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_inchikey_builder.c`

Parse errors: `4`. Function definitions: `6`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `KEYBUILDER_Unpack` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_inchikey_builder.c:59` | internal | tree-sitter+gcc-aux | `static INCHIKEYBUILDER *KEYBUILDER_Unpack (IXA_STATUS_HANDLE hStatus, IXA_INCHIKEYBUILDER_HANDLE hKeyBuilder)` |
| `KEYBUILDER_Pack` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_inchikey_builder.c:73` | internal | tree-sitter+gcc-aux | `static IXA_INCHIKEYBUILDER_HANDLE KEYBUILDER_Pack (INCHIKEYBUILDER *pKeyBuilder)` |
| `IXA_INCHIKEYBUILDER_Create` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_inchikey_builder.c:80` | external | tree-sitter+gcc-aux | `extern IXA_INCHIKEYBUILDER_HANDLE IXA_INCHIKEYBUILDER_Create (IXA_STATUS_HANDLE hStatus)` |
| `IXA_INCHIKEYBUILDER_Destroy` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_inchikey_builder.c:95` | external | tree-sitter+gcc-aux | `extern void IXA_INCHIKEYBUILDER_Destroy (IXA_STATUS_HANDLE hStatus, IXA_INCHIKEYBUILDER_HANDLE hKeyBuilder)` |
| `IXA_INCHIKEYBUILDER_SetInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_inchikey_builder.c:107` | external | tree-sitter+gcc-aux | `extern void IXA_INCHIKEYBUILDER_SetInChI (IXA_STATUS_HANDLE hStatus, IXA_INCHIKEYBUILDER_HANDLE hKeyBuilder, const char *pInChI)` |
| `IXA_INCHIKEYBUILDER_GetInChIKey` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_inchikey_builder.c:131` | external | tree-sitter+gcc-aux | `extern const char *IXA_INCHIKEYBUILDER_GetInChIKey (IXA_STATUS_HANDLE hStatus, IXA_INCHIKEYBUILDER_HANDLE hKeyBuilder)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c`

Parse errors: `61`. Function definitions: `85`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `GetVertexCount` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:118` | internal | tree-sitter+gcc-aux | `static int GetVertexCount (IXA_STATUS_HANDLE hStatus, IXA_STEREO_TOPOLOGY vTopology)` |
| `MOL_PackAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:135` | internal | tree-sitter+gcc-aux | `static IXA_ATOMID MOL_PackAtom (int vAtomIndex)` |
| `MOL_PackBond` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:142` | internal | tree-sitter+gcc-aux | `static IXA_BONDID MOL_PackBond (int vBondIndex)` |
| `MOL_PackStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:149` | internal | tree-sitter+gcc-aux | `static IXA_STEREOID MOL_PackStereo (int vStereoIndex)` |
| `MOL_PackPolymerUnit` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:156` | internal | tree-sitter+gcc-aux | `static IXA_POLYMERUNITID MOL_PackPolymerUnit (int vPunitIndex)` |
| `MOL_UnpackAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:163` | internal | tree-sitter+gcc-aux | `static IXA_BOOL MOL_UnpackAtom (IXA_STATUS_HANDLE hStatus, INCHIMOL *pMolecule, IXA_ATOMID vAtom, int *pAtomIndex)` |
| `MOL_UnpackBond` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:180` | internal | tree-sitter+gcc-aux | `static IXA_BOOL MOL_UnpackBond (IXA_STATUS_HANDLE hStatus, INCHIMOL *pMolecule, IXA_BONDID vBond, int *pBondIndex)` |
| `MOL_UnpackStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:197` | internal | tree-sitter+gcc-aux | `static IXA_BOOL MOL_UnpackStereo (IXA_STATUS_HANDLE hStatus, INCHIMOL *pMolecule, IXA_STEREOID vStereo, int *pStereoIndex)` |
| `MOL_UnpackPolymerUnit` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:214` | internal | tree-sitter+gcc-aux | `static IXA_BOOL MOL_UnpackPolymerUnit (IXA_STATUS_HANDLE hStatus, INCHIMOL *pMolecule, IXA_POLYMERUNITID vPunit, int *pPunitIndex)` |
| `MOL_GetAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:238` | internal | tree-sitter+gcc-aux | `static INCHIMOL_ATOM *MOL_GetAtom (IXA_STATUS_HANDLE hStatus, INCHIMOL *pMolecule, IXA_ATOMID vAtom)` |
| `MOL_GetBond` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:249` | internal | tree-sitter+gcc-aux | `static INCHIMOL_BOND *MOL_GetBond (IXA_STATUS_HANDLE hStatus, INCHIMOL *pMolecule, IXA_BONDID vBond)` |
| `MOL_GetStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:260` | internal | tree-sitter+gcc-aux | `static INCHIMOL_STEREO *MOL_GetStereo (IXA_STATUS_HANDLE hStatus, INCHIMOL *pMolecule, IXA_STEREOID vStereo)` |
| `MOL_GetSGroup` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:271` | internal | tree-sitter+gcc-aux | `static INCHIMOL_SGROUP *MOL_GetSGroup (IXA_STATUS_HANDLE hStatus, INCHIMOL *pMolecule, IXA_POLYMERUNITID vPunit)` |
| `MOL_ClearExtMolData` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:285` | external | tree-sitter+gcc-aux | `extern void MOL_ClearExtMolData (INCHIMOL_POLYMER *pd, INCHIMOL_V3000 *v3k)` |
| `MOL_Clear` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:395` | external | tree-sitter+gcc-aux | `extern void MOL_Clear (INCHIMOL *pMolecule)` |
| `MOL_GuessNewSize` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:408` | internal | tree-sitter+gcc-aux | `static int MOL_GuessNewSize (int reserved, int start_size, int max_size)` |
| `MOL_CreateAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:428` | internal | tree-sitter+gcc-aux | `static IXA_ATOMID MOL_CreateAtom (IXA_STATUS_HANDLE hStatus, INCHIMOL *pMolecule)` |
| `MOL_CreateStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:520` | internal | tree-sitter+gcc-aux | `static int MOL_CreateStereo (IXA_STATUS_HANDLE hStatus, INCHIMOL *pMolecule)` |
| `MOL_CreatePolymerUnit` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:598` | internal | tree-sitter+gcc-aux | `static IXA_POLYMERUNITID MOL_CreatePolymerUnit (IXA_STATUS_HANDLE hStatus, INCHIMOL *pMolecule)` |
| `MOL_Pack` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:671` | internal | tree-sitter+gcc-aux | `static IXA_MOL_HANDLE MOL_Pack (INCHIMOL *pMolecule)` |
| `MOL_Unpack` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:678` | external | tree-sitter+gcc-aux | `extern INCHIMOL *MOL_Unpack (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule)` |
| `MOL_GetBondOtherAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:692` | external | tree-sitter+gcc-aux | `extern IXA_ATOMID MOL_GetBondOtherAtom (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_BONDID vBond, IXA_ATOMID vAtom)` |
| `IXA_MOL_Create` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:729` | external | tree-sitter+gcc-aux | `extern IXA_MOL_HANDLE IXA_MOL_Create (IXA_STATUS_HANDLE hStatus)` |
| `IXA_MOL_Destroy` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:744` | external | tree-sitter+gcc-aux | `extern void IXA_MOL_Destroy (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule)` |
| `IXA_MOL_Clear` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:759` | external | tree-sitter+gcc-aux | `extern void IXA_MOL_Clear (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule)` |
| `IXA_MOL_SetChiral` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:773` | external | tree-sitter+gcc-aux | `extern void IXA_MOL_SetChiral (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_BOOL vChiral)` |
| `IXA_MOL_GetChiral` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:788` | external | tree-sitter+gcc-aux | `extern IXA_BOOL IXA_MOL_GetChiral (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule)` |
| `IXA_MOL_CreateAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:801` | external | tree-sitter+gcc-aux | `extern IXA_ATOMID IXA_MOL_CreateAtom (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule)` |
| `IXA_MOL_GetNumAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:816` | external | tree-sitter+gcc-aux | `extern int IXA_MOL_GetNumAtoms (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule)` |
| `IXA_MOL_GetAtomId` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:830` | external | tree-sitter+gcc-aux | `extern IXA_ATOMID IXA_MOL_GetAtomId (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, int vAtomIndex)` |
| `IXA_MOL_GetAtomIndex` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:857` | external | tree-sitter+gcc-aux | `extern int IXA_MOL_GetAtomIndex (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vAtom)` |
| `IXA_MOL_GetAtomNumBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:879` | external | tree-sitter+gcc-aux | `extern int IXA_MOL_GetAtomNumBonds (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vAtom)` |
| `IXA_MOL_GetAtomBond` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:902` | external | tree-sitter+gcc-aux | `extern IXA_BONDID IXA_MOL_GetAtomBond (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vAtom, int vBondIndex)` |
| `IXA_MOL_SetAtomX` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:929` | external | tree-sitter+gcc-aux | `extern void IXA_MOL_SetAtomX (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vAtom, double vX)` |
| `IXA_MOL_GetAtomX` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:953` | external | tree-sitter+gcc-aux | `extern double IXA_MOL_GetAtomX (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vAtom)` |
| `IXA_MOL_SetAtomY` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:976` | external | tree-sitter+gcc-aux | `extern void IXA_MOL_SetAtomY (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vAtom, double vY)` |
| `IXA_MOL_GetAtomY` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1000` | external | tree-sitter+gcc-aux | `extern double IXA_MOL_GetAtomY (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vAtom)` |
| `IXA_MOL_SetAtomZ` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1023` | external | tree-sitter+gcc-aux | `extern void IXA_MOL_SetAtomZ (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vAtom, double vZ)` |
| `IXA_MOL_GetAtomZ` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1047` | external | tree-sitter+gcc-aux | `extern double IXA_MOL_GetAtomZ (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vAtom)` |
| `IXA_MOL_SetAtomElement` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1070` | external | tree-sitter+gcc-aux | `extern void IXA_MOL_SetAtomElement (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vAtom, const char *pElement)` |
| `IXA_MOL_GetAtomElement` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1110` | external | tree-sitter+gcc-aux | `extern const char *IXA_MOL_GetAtomElement (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vAtom)` |
| `IXA_MOL_SetAtomAtomicNumber` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1133` | external | tree-sitter+gcc-aux | `extern void IXA_MOL_SetAtomAtomicNumber (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vAtom, int vAtomicNumber)` |
| `IXA_MOL_GetAtomAtomicNumber` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1163` | external | tree-sitter+gcc-aux | `extern int IXA_MOL_GetAtomAtomicNumber (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vAtom)` |
| `IXA_MOL_SetAtomHydrogens` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1186` | external | tree-sitter+gcc-aux | `extern void IXA_MOL_SetAtomHydrogens (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vAtom, int vHydrogenIsotope, int vHydrogenCount)` |
| `IXA_MOL_GetAtomHydrogens` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1236` | external | tree-sitter+gcc-aux | `extern int IXA_MOL_GetAtomHydrogens (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vAtom, int vHydrogenIsotope)` |
| `IXA_MOL_SetAtomMass` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1271` | external | tree-sitter+gcc-aux | `extern void IXA_MOL_SetAtomMass (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vAtom, int vMassNumber)` |
| `IXA_MOL_GetAtomMass` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1315` | external | tree-sitter+gcc-aux | `extern int IXA_MOL_GetAtomMass (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vAtom)` |
| `IXA_MOL_SetAtomRadical` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1338` | external | tree-sitter+gcc-aux | `extern void IXA_MOL_SetAtomRadical (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vAtom, IXA_ATOM_RADICAL vRadical)` |
| `IXA_MOL_GetAtomRadical` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1362` | external | tree-sitter+gcc-aux | `extern IXA_ATOM_RADICAL IXA_MOL_GetAtomRadical (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vAtom)` |
| `IXA_MOL_SetAtomCharge` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1385` | external | tree-sitter+gcc-aux | `extern void IXA_MOL_SetAtomCharge (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vAtom, int vCharge)` |
| `IXA_MOL_GetAtomCharge` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1409` | external | tree-sitter+gcc-aux | `extern int IXA_MOL_GetAtomCharge (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vAtom)` |
| `IXA_MOL_ReserveSpace` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1432` | external | tree-sitter+gcc-aux | `extern int IXA_MOL_ReserveSpace (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, int num_atoms, int num_bonds, int num_stereos)` |
| `IXA_MOL_CreateBond` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1497` | external | tree-sitter+gcc-aux | `extern IXA_BONDID IXA_MOL_CreateBond (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vAtom1, IXA_ATOMID vAtom2)` |
| `IXA_MOL_GetNumBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1643` | external | tree-sitter+gcc-aux | `extern int IXA_MOL_GetNumBonds (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule)` |
| `IXA_MOL_GetBondId` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1657` | external | tree-sitter+gcc-aux | `extern IXA_BONDID IXA_MOL_GetBondId (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, int vBondIndex)` |
| `IXA_MOL_GetBondIndex` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1684` | external | tree-sitter+gcc-aux | `extern int IXA_MOL_GetBondIndex (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_BONDID vBond)` |
| `IXA_MOL_GetBondAtom1` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1707` | external | tree-sitter+gcc-aux | `extern IXA_ATOMID IXA_MOL_GetBondAtom1 (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_BONDID vBond)` |
| `IXA_MOL_GetBondAtom2` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1732` | external | tree-sitter+gcc-aux | `extern IXA_ATOMID IXA_MOL_GetBondAtom2 (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_BONDID vBond)` |
| `IXA_MOL_GetBondOtherAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1756` | external | tree-sitter+gcc-aux | `extern IXA_ATOMID IXA_MOL_GetBondOtherAtom (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_BONDID vBond, IXA_ATOMID vAtom)` |
| `IXA_MOL_SetBondType` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1796` | external | tree-sitter+gcc-aux | `extern void IXA_MOL_SetBondType (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_BONDID vBond, IXA_BOND_TYPE vType)` |
| `IXA_MOL_GetBondType` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1822` | external | tree-sitter+gcc-aux | `extern IXA_BOND_TYPE IXA_MOL_GetBondType (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_BONDID vBond)` |
| `IXA_MOL_SetBondWedge` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1844` | external | tree-sitter+gcc-aux | `extern void IXA_MOL_SetBondWedge (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_BONDID vBond, IXA_ATOMID vRefAtom, IXA_BOND_WEDGE vDirection)` |
| `IXA_MOL_GetBondWedge` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1880` | external | tree-sitter | `IXA_BOND_WEDGE INCHI_DECL IXA_MOL_GetBondWedge( IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_BONDID vBond, IXA_ATOMID vRefAtom )` |
| `IXA_MOL_GetBondWedge` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1881` | external | gcc-aux | `extern IXA_BOND_WEDGE IXA_MOL_GetBondWedge (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_BONDID vBond, IXA_ATOMID vRefAtom)` |
| `IXA_MOL_SetDblBondConfig` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1917` | external | tree-sitter+gcc-aux | `extern void IXA_MOL_SetDblBondConfig (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_BONDID vBond, IXA_DBLBOND_CONFIG vConfig)` |
| `IXA_MOL_GetDblBondConfig` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1941` | external | tree-sitter+gcc-aux | `extern IXA_DBLBOND_CONFIG IXA_MOL_GetDblBondConfig (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_BONDID vBond)` |
| `IXA_MOL_GetCommonBond` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:1964` | external | tree-sitter+gcc-aux | `extern IXA_BONDID IXA_MOL_GetCommonBond (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vAtom1, IXA_ATOMID vAtom2)` |
| `IXA_MOL_CreateStereoTetrahedron` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2010` | external | tree-sitter+gcc-aux | `extern IXA_STEREOID IXA_MOL_CreateStereoTetrahedron (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vCentralAtom, IXA_ATOMID vVertex1, IXA_ATOMID vVertex2, IXA_ATOMID vVertex3, IXA_ATOMID vVertex4)` |
| `IXA_MOL_CreateStereoRectangle` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2045` | external | tree-sitter+gcc-aux | `extern IXA_STEREOID IXA_MOL_CreateStereoRectangle (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_BONDID vCentralBond, IXA_ATOMID vVertex1, IXA_ATOMID vVertex2, IXA_ATOMID vVertex3, IXA_ATOMID vVertex4)` |
| `IXA_MOL_CreateStereoAntiRectangle` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2079` | external | tree-sitter+gcc-aux | `extern IXA_STEREOID IXA_MOL_CreateStereoAntiRectangle (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vCentralAtom, IXA_ATOMID vVertex1, IXA_ATOMID vVertex2, IXA_ATOMID vVertex3, IXA_ATOMID vVertex4)` |
| `IXA_MOL_GetNumStereos` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2113` | external | tree-sitter+gcc-aux | `extern int IXA_MOL_GetNumStereos (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule)` |
| `IXA_MOL_GetStereoId` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2127` | external | tree-sitter+gcc-aux | `extern IXA_STEREOID IXA_MOL_GetStereoId (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, int vStereoIndex)` |
| `IXA_MOL_GetStereoIndex` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2148` | external | tree-sitter+gcc-aux | `extern int IXA_MOL_GetStereoIndex (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_STEREOID vStereo)` |
| `IXA_MOL_GetStereoTopology` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2164` | external | tree-sitter+gcc-aux | `extern IXA_STEREO_TOPOLOGY IXA_MOL_GetStereoTopology (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_STEREOID vStereo)` |
| `IXA_MOL_GetStereoCentralAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2186` | external | tree-sitter+gcc-aux | `extern IXA_ATOMID IXA_MOL_GetStereoCentralAtom (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_STEREOID vStereo)` |
| `IXA_MOL_GetStereoCentralBond` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2219` | external | tree-sitter+gcc-aux | `extern IXA_BONDID IXA_MOL_GetStereoCentralBond (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_STEREOID vStereo)` |
| `IXA_MOL_GetStereoNumVertices` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2252` | external | tree-sitter+gcc-aux | `extern int IXA_MOL_GetStereoNumVertices (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_STEREOID vStereo)` |
| `IXA_MOL_GetStereoVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2274` | external | tree-sitter+gcc-aux | `extern IXA_ATOMID IXA_MOL_GetStereoVertex (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_STEREOID vStereo, int vVertexIndex)` |
| `IXA_MOL_SetStereoParity` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2310` | external | tree-sitter+gcc-aux | `extern void IXA_MOL_SetStereoParity (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_STEREOID vStereo, IXA_STEREO_PARITY vParity)` |
| `IXA_MOL_GetStereoParity` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2331` | external | tree-sitter+gcc-aux | `extern IXA_STEREO_PARITY IXA_MOL_GetStereoParity (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_STEREOID vStereo)` |
| `IXA_MOL_SetPolymerUnit` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2353` | external | tree-sitter+gcc-aux | `extern void IXA_MOL_SetPolymerUnit (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_POLYMERUNITID vPunit, int vid, int vtype, int vsubtype, int vconn, int vlabel, int vna, int vnb, double *vxbr1, double *vxbr2, char *vsmt, int *valist, int *vblist)` |
| `IXA_MOL_SetExtMolDataByInChIExtInput` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2429` | external | tree-sitter+gcc-aux | `extern int IXA_MOL_SetExtMolDataByInChIExtInput (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, inchi_Output_Polymer *polymer, inchi_Output_V3000 *v3000, int nat)` |
| `IXA_MOL_CreatePolymerUnit` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2681` | external | tree-sitter+gcc-aux | `extern IXA_POLYMERUNITID IXA_MOL_CreatePolymerUnit (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule)` |
| `IXA_MOL_GetPolymerUnitId` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2696` | external | tree-sitter+gcc-aux | `extern IXA_POLYMERUNITID IXA_MOL_GetPolymerUnitId (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, int vPolymerUnitIndex)` |
| `IXA_MOL_GetPolymerUnitIndex` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_mol.c:2717` | external | tree-sitter+gcc-aux | `extern int IXA_MOL_GetPolymerUnitIndex (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_POLYMERUNITID vPolymerUnit)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_read_inchi.c`

Parse errors: `2`. Function definitions: `4`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `AnalyseInternalVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_read_inchi.c:55` | internal | tree-sitter+gcc-aux | `static IXA_ATOMID AnalyseInternalVertex (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vVertex, IXA_ATOMID vInternal, IXA_ATOMID vIgnore1, IXA_ATOMID vIgnore2)` |
| `FindCumuleneCentre` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_read_inchi.c:112` | internal | tree-sitter+gcc-aux | `static IXA_BONDID FindCumuleneCentre (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, IXA_ATOMID vInternal1, IXA_ATOMID vInternal2)` |
| `IXA_MOL_ReadInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_read_inchi.c:189` | external | tree-sitter+gcc-aux | `extern void IXA_MOL_ReadInChI (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, const char *pInChI)` |
| `IXA_MOL_ReadAuxInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_read_inchi.c:686` | external | tree-sitter+gcc-aux | `extern void IXA_MOL_ReadAuxInfo (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, char *pAuxInfo, int bDoNotAddH, int bDiffUnkUndfStereo)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_read_mol.c`

Parse errors: `2`. Function definitions: `1`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `IXA_MOL_ReadMolfile` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_read_mol.c:78` | external | tree-sitter+gcc-aux | `extern void IXA_MOL_ReadMolfile (IXA_STATUS_HANDLE hStatus, IXA_MOL_HANDLE hMolecule, const char *pBytes)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c`

Parse errors: `8`. Function definitions: `15`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `BLOCK_clear` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:81` | internal | tree-sitter+gcc-aux | `static void BLOCK_clear (StatusBlock *pBlock)` |
| `STATUS_init` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:93` | internal | tree-sitter+gcc-aux | `static void STATUS_init (INCHISTATUS *pStatus)` |
| `STATUS_Clear` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:102` | internal | tree-sitter+gcc-aux | `static void STATUS_Clear (INCHISTATUS *pStatus)` |
| `STATUS_Pack` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:123` | external | tree-sitter+gcc-aux | `extern IXA_STATUS_HANDLE STATUS_Pack (INCHISTATUS *pStatus)` |
| `STATUS_Unpack` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:131` | internal | tree-sitter+gcc-aux | `static INCHISTATUS *STATUS_Unpack (IXA_STATUS_HANDLE hStatus)` |
| `STATUS_PushMessage` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:138` | external | tree-sitter+gcc-aux | `extern void STATUS_PushMessage (IXA_STATUS_HANDLE hStatus, IXA_STATUS vSeverity, char *pFormat, ...)` |
| `INCHISTATUS_TestSeverity` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:211` | internal | tree-sitter+gcc-aux | `static IXA_BOOL INCHISTATUS_TestSeverity (IXA_STATUS_HANDLE hStatus, IXA_STATUS vSeverity)` |
| `IXA_STATUS_Create` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:233` | external | tree-sitter+gcc-aux | `extern IXA_STATUS_HANDLE IXA_STATUS_Create (void)` |
| `IXA_STATUS_Destroy` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:246` | external | tree-sitter+gcc-aux | `extern void IXA_STATUS_Destroy (IXA_STATUS_HANDLE hStatus)` |
| `IXA_STATUS_Clear` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:259` | external | tree-sitter+gcc-aux | `extern void IXA_STATUS_Clear (IXA_STATUS_HANDLE hStatus)` |
| `IXA_STATUS_HasError` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:271` | external | tree-sitter+gcc-aux | `extern IXA_BOOL IXA_STATUS_HasError (IXA_STATUS_HANDLE hStatus)` |
| `IXA_STATUS_HasWarning` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:278` | external | tree-sitter+gcc-aux | `extern IXA_BOOL IXA_STATUS_HasWarning (IXA_STATUS_HANDLE hStatus)` |
| `IXA_STATUS_GetCount` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:285` | external | tree-sitter+gcc-aux | `extern int IXA_STATUS_GetCount (IXA_STATUS_HANDLE hStatus)` |
| `IXA_STATUS_GetMessage` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:298` | external | tree-sitter+gcc-aux | `extern const char *IXA_STATUS_GetMessage (IXA_STATUS_HANDLE hStatus, int vIndex)` |
| `IXA_STATUS_GetSeverity` | `third_party/InChI/INCHI-1-SRC/INCHI_API/libinchi/src/ixa/ixa_status.c:316` | external | tree-sitter+gcc-aux | `extern IXA_STATUS IXA_STATUS_GetSeverity (IXA_STATUS_HANDLE hStatus, int vIndex)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/bcf_s.c`

Parse errors: `0`. Function definitions: `6`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `max_3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/bcf_s.c:58` | external | tree-sitter+gcc-aux | `extern int max_3 (int a, int b, int c)` |
| `memcpy_custom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/bcf_s.c:68` | external | tree-sitter+gcc-aux | `extern int memcpy_custom (char **dst, char *src, long long unsigned int len)` |
| `dbl2int` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/bcf_s.c:92` | external | tree-sitter+gcc-aux | `extern int dbl2int (char *str, int fwidth, int ndecpl, char dbl_flag, double dblinp)` |
| `dbl2int_f` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/bcf_s.c:117` | internal | tree-sitter+gcc-aux | `static int dbl2int_f (double dblinp, int fwidth, int ndecpl, char *str)` |
| `dbl2int_e` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/bcf_s.c:188` | internal | tree-sitter+gcc-aux | `static int dbl2int_e (double dblinp, int fwidth, int ndecpl, char *str)` |
| `dbl2int_g` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/bcf_s.c:276` | internal | tree-sitter+gcc-aux | `static int dbl2int_g (double dblinp, int fwidth, int ndecpl, char *str)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c`

Parse errors: `33`. Function definitions: `111`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `RestoreEdgeFlow` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:442` | external | tree-sitter+gcc-aux | `extern int RestoreEdgeFlow (BNS_EDGE *edge, int delta, int bChangeFlow)` |
| `SetAtomBondType` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:480` | external | tree-sitter+gcc-aux | `extern int SetAtomBondType (BNS_EDGE *edge, U_CHAR *bond_type12, U_CHAR *bond_type21, int delta, int bChangeFlow)` |
| `RunBalancedNetworkSearch` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:689` | external | tree-sitter+gcc-aux | `extern int RunBalancedNetworkSearch (BN_STRUCT *pBNS, BN_DATA *pBD, int bChangeFlow)` |
| `SetAtomRadAndChemValFromVertexCapFlow` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:726` | external | tree-sitter+gcc-aux | `extern int SetAtomRadAndChemValFromVertexCapFlow (BN_STRUCT *pBNS, inp_ATOM *atom, int v1)` |
| `AddChangedAtHChargeBNS` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:775` | external | tree-sitter+gcc-aux | `extern int AddChangedAtHChargeBNS (inp_ATOM *at, int num_atoms, int *nAtTypeTotals, S_CHAR *mark)` |
| `EliminatePlusMinusChargeAmbiguity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:812` | external | tree-sitter+gcc-aux | `extern int EliminatePlusMinusChargeAmbiguity (BN_STRUCT *pBNS, int num_atoms)` |
| `AddOrRemoveExplOrImplH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:902` | external | tree-sitter+gcc-aux | `extern int AddOrRemoveExplOrImplH (int nDelta, inp_ATOM *at, int num_atoms, AT_NUMB at_no, T_GROUP_INFO *t_group_info)` |
| `SubtractOrChangeAtHChargeBNS` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:1057` | external | tree-sitter+gcc-aux | `extern int SubtractOrChangeAtHChargeBNS (BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, int *nAtTypeTotals, S_CHAR *mark, T_GROUP_INFO *t_group_info, int bSubtract)` |
| `SetBondsFromBnStructFlow` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:1165` | external | tree-sitter+gcc-aux | `extern int SetBondsFromBnStructFlow (BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, int bChangeFlow0)` |
| `MarkAtomsAtTautGroups` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:1331` | external | tree-sitter+gcc-aux | `extern int MarkAtomsAtTautGroups (BN_STRUCT *pBNS, int num_atoms, BN_AATG *pAATG, int nEnd1, int nEnd2)` |
| `RestoreBnStructFlow` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:1525` | external | tree-sitter+gcc-aux | `extern int RestoreBnStructFlow (BN_STRUCT *pBNS, int bChangeFlow)` |
| `bNeedToTestTheFlow` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:1590` | external | tree-sitter+gcc-aux | `extern int bNeedToTestTheFlow (int bond_type, int nTestFlow, int bTestForNonStereoBond)` |
| `nBondsValenceInpAt` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:1673` | external | tree-sitter+gcc-aux | `extern int nBondsValenceInpAt (const inp_ATOM *at, int *nNumAltBonds, int *nNumWrongBonds)` |
| `BnsAdjustFlowBondsRad` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:1725` | external | tree-sitter+gcc-aux | `extern int BnsAdjustFlowBondsRad (BN_STRUCT *pBNS, BN_DATA *pBD, inp_ATOM *at, int num_atoms)` |
| `BnsTestAndMarkAltBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:1825` | external | tree-sitter+gcc-aux | `extern int BnsTestAndMarkAltBonds (BN_STRUCT *pBNS, BN_DATA *pBD, inp_ATOM *at, int num_atoms, BNS_FLOW_CHANGES *fcd, int bChangeFlow, int nBondTypeToTest)` |
| `remove_alt_bond_marks` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:1987` | internal | tree-sitter+gcc-aux | `static void remove_alt_bond_marks (inp_ATOM *at, int num_atoms)` |
| `SetForbiddenEdges` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:2002` | external | tree-sitter+gcc-aux | `extern int SetForbiddenEdges (BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, int forbidden_mask, int nebend, int *ebend)` |
| `TempFix_NH_NH_Bonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:2108` | external | tree-sitter | `int TempFix_NH_NH_Bonds( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms )` |
| `CorrectFixing_NH_NH_Bonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:2143` | external | tree-sitter | `int CorrectFixing_NH_NH_Bonds( BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms )` |
| `fix_special_bonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:2181` | external | tree-sitter+gcc-aux | `extern int fix_special_bonds (BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, int forbidden_mask)` |
| `fix_explicitly_indicated_bonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:2704` | external | tree-sitter+gcc-aux | `extern int fix_explicitly_indicated_bonds (int nebend, int *ebend, BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms)` |
| `is_Z_atom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:2781` | external | tree-sitter+gcc-aux | `extern int is_Z_atom (U_CHAR el_number)` |
| `IsZOX` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:2814` | external | tree-sitter+gcc-aux | `extern int IsZOX (inp_ATOM *atom, int at_x, int ord)` |
| `update_some_attype_totals` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:2846` | external | tree-sitter+gcc-aux | `extern void update_some_attype_totals (int *nAtTypeTotals, int mask, int delta, S_CHAR at_charge)` |
| `GetAtomChargeType` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:2873` | external | tree-sitter+gcc-aux | `extern int GetAtomChargeType (inp_ATOM *atom, int at_no, int *nAtTypeTotals, int *pMask, int bSubtract)` |
| `SimpleRemoveHplusNPO` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:3432` | external | tree-sitter+gcc-aux | `extern int SimpleRemoveHplusNPO (inp_ATOM *at, int num_atoms, int *nAtTypeTotals, T_GROUP_INFO *t_group_info)` |
| `bIsAtomTypeHard` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:3480` | external | tree-sitter+gcc-aux | `extern int bIsAtomTypeHard (inp_ATOM *at, int endpoint, int nType, int nMask, int nCharge)` |
| `bIsHDonorAccAtomType` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:3501` | external | tree-sitter+gcc-aux | `extern int bIsHDonorAccAtomType (inp_ATOM *at, int endpoint, int *cSubType)` |
| `bIsNegAtomType` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:3539` | external | tree-sitter+gcc-aux | `extern int bIsNegAtomType (inp_ATOM *at, int endpoint, int *cSubType)` |
| `bIsHardRemHCandidate` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:3582` | external | tree-sitter+gcc-aux | `extern int bIsHardRemHCandidate (inp_ATOM *at, int i, int *cSubType)` |
| `CreateCGroupInBnStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:3599` | external | tree-sitter+gcc-aux | `extern int CreateCGroupInBnStruct (inp_ATOM *at, int num_atoms, BN_STRUCT *pBNS, int nType, int nMask, int nCharge)` |
| `CreateTGroupInBnStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:3813` | external | tree-sitter+gcc-aux | `extern int CreateTGroupInBnStruct (inp_ATOM *at, int num_atoms, BN_STRUCT *pBNS, int nType, int nMask)` |
| `RemoveLastGroupFromBnStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:3988` | external | tree-sitter+gcc-aux | `extern int RemoveLastGroupFromBnStruct (inp_ATOM *at, int num_atoms, int tg, BN_STRUCT *pBNS)` |
| `SetInitCapFlowToCurrent` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:4094` | external | tree-sitter+gcc-aux | `extern int SetInitCapFlowToCurrent (BN_STRUCT *pBNS)` |
| `SimpleRemoveAcidicProtons` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:4129` | external | tree-sitter+gcc-aux | `extern int SimpleRemoveAcidicProtons (inp_ATOM *at, int num_atoms, BN_AATG *pAATG, int num2remove)` |
| `bHasAcidicHydrogen` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:4200` | external | tree-sitter+gcc-aux | `extern int bHasAcidicHydrogen (inp_ATOM *at, int i)` |
| `bHasOtherExchangableH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:4220` | external | tree-sitter+gcc-aux | `extern int bHasOtherExchangableH (inp_ATOM *at, int i)` |
| `SimpleAddAcidicProtons` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:4253` | external | tree-sitter+gcc-aux | `extern int SimpleAddAcidicProtons (inp_ATOM *at, int num_atoms, BN_AATG *pAATG, int num2add)` |
| `bHasAcidicMinus` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:4324` | external | tree-sitter+gcc-aux | `extern int bHasAcidicMinus (inp_ATOM *at, int i)` |
| `HardRemoveAcidicProtons` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:4351` | external | tree-sitter+gcc-aux | `extern int HardRemoveAcidicProtons (CANON_GLOBALS *pCG, inp_ATOM *at, int num_atoms, BN_AATG *pAATG, int num2remove, int *nNumCanceledCharges, BN_STRUCT *pBNS, BN_DATA *pBD)` |
| `HardAddAcidicProtons` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:4556` | external | tree-sitter+gcc-aux | `extern int HardAddAcidicProtons (CANON_GLOBALS *pCG, inp_ATOM *at, int num_atoms, BN_AATG *pAATG, int num2add, int *nNumCanceledCharges, BN_STRUCT *pBNS, BN_DATA *pBD)` |
| `HardRemoveHplusNP` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:4770` | external | tree-sitter+gcc-aux | `extern int HardRemoveHplusNP (CANON_GLOBALS *pCG, inp_ATOM *at, int num_atoms, int bCancelChargesAlways, int *nNumCanceledCharges, BN_AATG *pAATG, BN_STRUCT *pBNS, BN_DATA *pBD)` |
| `mark_at_type` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:5040` | external | tree-sitter+gcc-aux | `extern int mark_at_type (inp_ATOM *atom, int num_atoms, int *nAtTypeTotals)` |
| `RemoveNPProtonsAndAcidCharges` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:5084` | external | tree-sitter+gcc-aux | `extern int RemoveNPProtonsAndAcidCharges (CANON_GLOBALS *pCG, inp_ATOM *at, int num_atoms, BN_AATG *pAATG, BN_STRUCT *pBNS, BN_DATA *pBD)` |
| `mark_alt_bonds_and_taut_groups` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:5242` | external | tree-sitter+gcc-aux | `extern int mark_alt_bonds_and_taut_groups (struct tagINCHI_CLOCK *ic, struct tagCANON_GLOBALS *pCG, inp_ATOM *at, inp_ATOM *at_fixed_bonds_out, int num_atoms, struct tagInchiTime *ulTimeOutTime, T_GROUP_INFO *t_group_info, INCHI_MODE *inpbTautFlags, INCHI_MODE *inpbTautFlagsDone, int nebend, int *ebend)` |
| `nMaxFlow2Check` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:6010` | external | tree-sitter+gcc-aux | `extern int nMaxFlow2Check (BN_STRUCT *pBNS, int iedge)` |
| `nCurFlow2Check` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:6024` | external | tree-sitter+gcc-aux | `extern int nCurFlow2Check (BN_STRUCT *pBNS, int iedge)` |
| `nMinFlow2Check` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:6033` | external | tree-sitter+gcc-aux | `extern int nMinFlow2Check (BN_STRUCT *pBNS, int iedge)` |
| `bSetBondsAfterCheckOneBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:6071` | external | tree-sitter+gcc-aux | `extern int bSetBondsAfterCheckOneBond (BN_STRUCT *pBNS, BNS_FLOW_CHANGES *fcd, int nTestFlow, inp_ATOM *at, int num_atoms, int bChangeFlow0)` |
| `bRestoreFlowAfterCheckOneBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:6180` | external | tree-sitter+gcc-aux | `extern int bRestoreFlowAfterCheckOneBond (BN_STRUCT *pBNS, BNS_FLOW_CHANGES *fcd)` |
| `bSetFlowToCheckOneBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:6222` | external | tree-sitter+gcc-aux | `extern int bSetFlowToCheckOneBond (BN_STRUCT *pBNS, int iedge, int flow, BNS_FLOW_CHANGES *fcd)` |
| `bAddNewVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:6482` | external | tree-sitter+gcc-aux | `extern int bAddNewVertex (BN_STRUCT *pBNS, int nVertDoubleBond, int nCap, int nFlow, int nMaxAdjEdges, int *nDots)` |
| `AddNewEdge` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:6556` | external | tree-sitter+gcc-aux | `extern int AddNewEdge (BNS_VERTEX *p1, BNS_VERTEX *p2, BN_STRUCT *pBNS, int nEdgeCap, int nEdgeFlow)` |
| `GetEdgeToGroupVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:6610` | external | tree-sitter+gcc-aux | `extern BNS_IEDGE GetEdgeToGroupVertex (BN_STRUCT *pBNS, Vertex v1, AT_NUMB type)` |
| `GetGroupVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:6644` | external | tree-sitter+gcc-aux | `extern Vertex GetGroupVertex (BN_STRUCT *pBNS, Vertex v1, AT_NUMB type)` |
| `bAddStCapToAVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:6703` | external | tree-sitter+gcc-aux | `extern int bAddStCapToAVertex (BN_STRUCT *pBNS, Vertex v1, Vertex v2, VertexFlow *nOldCapVertSingleBond, int *nDots, int bAdjacentDonors)` |
| `bSetBnsToCheckAltPath` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:6778` | external | tree-sitter+gcc-aux | `extern int bSetBnsToCheckAltPath (BN_STRUCT *pBNS, int nVertDoubleBond, int nVertSingleBond, AT_NUMB type, int path_type, ALT_PATH_CHANGES *apc, BNS_FLOW_CHANGES *fcd, int *nDots)` |
| `bRestoreBnsAfterCheckAltPath` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:7398` | external | tree-sitter+gcc-aux | `extern int bRestoreBnsAfterCheckAltPath (BN_STRUCT *pBNS, ALT_PATH_CHANGES *apc, int bChangeFlow)` |
| `bExistsAnyAltPath` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:7518` | external | tree-sitter+gcc-aux | `extern int bExistsAnyAltPath (CANON_GLOBALS *pCG, BN_STRUCT *pBNS, BN_DATA *pBD, inp_ATOM *at, int num_atoms, int nVert2, int nVert1, int path_type)` |
| `bIsBnsEndpoint` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:7561` | external | tree-sitter+gcc-aux | `extern int bIsBnsEndpoint (BN_STRUCT *pBNS, int v)` |
| `bRadChangesAtomType` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:7588` | external | tree-sitter+gcc-aux | `extern int bRadChangesAtomType (BN_STRUCT *pBNS, BN_DATA *pBD, Vertex v, Vertex v_1, Vertex v_2)` |
| `RegisterRadEndpoint` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:7639` | external | tree-sitter+gcc-aux | `extern int RegisterRadEndpoint (BN_STRUCT *pBNS, BN_DATA *pBD, Vertex u)` |
| `cmp_rad_endpoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:7784` | external | tree-sitter+gcc-aux | `extern int cmp_rad_endpoints (const void *a1, const void *a2)` |
| `RemoveRadEndpoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:7812` | external | tree-sitter+gcc-aux | `extern int RemoveRadEndpoints (BN_STRUCT *pBNS, BN_DATA *pBD, inp_ATOM *at)` |
| `RestoreRadicalsOnly` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:7907` | external | tree-sitter+gcc-aux | `extern int RestoreRadicalsOnly (BN_STRUCT *pBNS, BN_DATA *pBD, inp_ATOM *at)` |
| `SetRadEndpoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:7969` | external | tree-sitter+gcc-aux | `extern int SetRadEndpoints (BN_STRUCT *pBNS, BN_DATA *pBD, BRS_MODE bRadSrchMode)` |
| `SetRadEndpoints2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:8064` | external | tree-sitter+gcc-aux | `extern int SetRadEndpoints2 (CANON_GLOBALS *pCG, BN_STRUCT *pBNS, BN_DATA *pBD, BRS_MODE bRadSrchMode)` |
| `SetRadEndpoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:8369` | external | tree-sitter | `int SetRadEndpoints( BN_STRUCT *pBNS, BN_DATA *pBD, BRS_MODE bRadSrchMode )` |
| `RemoveRadEndpoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:8373` | external | tree-sitter | `int RemoveRadEndpoints( BN_STRUCT *pBNS, BN_DATA *pBD, inp_ATOM *at )` |
| `SetRadEndpoints2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:8377` | external | tree-sitter | `int SetRadEndpoints2( CANON_GLOBALS *pCG, BN_STRUCT *pBNS, BN_DATA *pBD, BRS_MODE bRadSrchMode )` |
| `bExistsAltPath` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:8393` | external | tree-sitter+gcc-aux | `extern int bExistsAltPath (CANON_GLOBALS *pCG, BN_STRUCT *pBNS, BN_DATA *pBD, BN_AATG *pAATG, inp_ATOM *at, int num_atoms, int nVertDoubleBond, int nVertSingleBond, int path_type)` |
| `AllocateAndInitBnStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:8783` | external | tree-sitter+gcc-aux | `extern BN_STRUCT *AllocateAndInitBnStruct (inp_ATOM *at, int num_atoms, int nMaxAddAtoms, int nMaxAddEdges, int max_altp, int *pNum_changed_bonds)` |
| `DeAllocateBnStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9010` | external | tree-sitter+gcc-aux | `extern BN_STRUCT *DeAllocateBnStruct (BN_STRUCT *pBNS)` |
| `ReInitBnStructAltPaths` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9042` | external | tree-sitter+gcc-aux | `extern int ReInitBnStructAltPaths (BN_STRUCT *pBNS)` |
| `ReInitBnStructAddGroups` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9062` | external | tree-sitter+gcc-aux | `extern int ReInitBnStructAddGroups (CANON_GLOBALS *pCG, BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, T_GROUP_INFO *tgi, C_GROUP_INFO *cgi)` |
| `ReInitBnStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9102` | external | tree-sitter+gcc-aux | `extern int ReInitBnStruct (BN_STRUCT *pBNS, inp_ATOM *at, int num_at, int bRemoveGroupsFromAtoms)` |
| `CompTGroupNumber` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9195` | external | tree-sitter+gcc-aux | `extern int CompTGroupNumber (const void *tg1, const void *tg2, void *p)` |
| `CompCGroupNumber` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9202` | external | tree-sitter+gcc-aux | `extern int CompCGroupNumber (const void *cg1, const void *cg2, void *p)` |
| `AddTGroups2BnStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9209` | external | tree-sitter+gcc-aux | `extern int AddTGroups2BnStruct (CANON_GLOBALS *pCG, BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, T_GROUP_INFO *tgi)` |
| `AddCGroups2BnStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9406` | external | tree-sitter+gcc-aux | `extern int AddCGroups2BnStruct (CANON_GLOBALS *pCG, BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, C_GROUP_INFO *cgi)` |
| `ClearAllBnDataVertices` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9609` | external | tree-sitter+gcc-aux | `extern void ClearAllBnDataVertices (Vertex *v, Vertex value, int size)` |
| `ClearAllBnDataEdges` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9620` | external | tree-sitter+gcc-aux | `extern void ClearAllBnDataEdges (Edge (*e), Vertex value, int size)` |
| `DeAllocateBnData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9631` | external | tree-sitter+gcc-aux | `extern BN_DATA *DeAllocateBnData (BN_DATA *pBD)` |
| `AllocateAndInitBnData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9677` | external | tree-sitter+gcc-aux | `extern BN_DATA *AllocateAndInitBnData (int max_num_vertices)` |
| `ReInitBnData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9718` | external | tree-sitter+gcc-aux | `extern int ReInitBnData (BN_DATA *pBD)` |
| `GetVertexDegree` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9774` | external | tree-sitter+gcc-aux | `extern int GetVertexDegree (BN_STRUCT *pBNS, Vertex v)` |
| `Get2ndEdgeVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9796` | external | tree-sitter+gcc-aux | `extern Vertex Get2ndEdgeVertex (BN_STRUCT *pBNS, int *uv)` |
| `GetVertexNeighbor` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9834` | external | tree-sitter+gcc-aux | `extern Vertex GetVertexNeighbor (BN_STRUCT *pBNS, Vertex v, int neigh, EdgeIndex *iedge)` |
| `GetEdgePointer` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9879` | external | tree-sitter+gcc-aux | `extern int GetEdgePointer (BN_STRUCT *pBNS, Vertex u, Vertex v, EdgeIndex iuv, BNS_EDGE **uv, S_CHAR *s_or_t)` |
| `AugmentEdge` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:9934` | external | tree-sitter+gcc-aux | `extern int AugmentEdge (BN_STRUCT *pBNS, Vertex u, Vertex v, EdgeIndex iuv, int delta, S_CHAR bReverse, int bChangeFlow)` |
| `rescap_mark` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:10102` | external | tree-sitter+gcc-aux | `extern int rescap_mark (BN_STRUCT *pBNS, Vertex u, Vertex v, EdgeIndex iuv)` |
| `GetPrevVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:10163` | external | tree-sitter+gcc-aux | `extern Vertex GetPrevVertex (BN_STRUCT *pBNS, Vertex y, Edge (*SwitchEdge), EdgeIndex *iuv)` |
| `bIgnoreVertexNonTACN_atom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:10234` | external | tree-sitter | `int bIgnoreVertexNonTACN_atom( BN_STRUCT* pBNS, Vertex u, Vertex v )` |
| `bIgnoreVertexNonTACN_atom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:10291` | external | tree-sitter+gcc-aux | `extern int bIgnoreVertexNonTACN_atom (BN_STRUCT *pBNS, Vertex u, Vertex v)` |
| `bIgnoreVertexNonTACN_group` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:10360` | external | tree-sitter | `int bIgnoreVertexNonTACN_group( BN_STRUCT* pBNS, Vertex v, Vertex w, Edge *SwitchEdge )` |
| `bIgnoreVertexNonTACN_group` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:10419` | external | tree-sitter+gcc-aux | `extern int bIgnoreVertexNonTACN_group (BN_STRUCT *pBNS, Vertex v, Vertex w, Edge (*SwitchEdge))` |
| `bIsRemovedHfromNHaion` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:10483` | external | tree-sitter | `int bIsRemovedHfromNHaion( BN_STRUCT* pBNS, Vertex u, Vertex v )` |
| `bIsAggressiveDeprotonation` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:10550` | external | tree-sitter | `int bIsAggressiveDeprotonation( BN_STRUCT* pBNS, Vertex v, Vertex w, Edge *SwitchEdge )` |
| `rescap` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:10644` | external | tree-sitter+gcc-aux | `extern int rescap (BN_STRUCT *pBNS, Vertex u, Vertex v, EdgeIndex iuv)` |
| `BalancedNetworkSearch` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:10779` | external | tree-sitter+gcc-aux | `extern int BalancedNetworkSearch (BN_STRUCT *pBNS, BN_DATA *pBD, int bChangeFlow)` |
| `FindBase` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:11068` | external | tree-sitter+gcc-aux | `extern Vertex FindBase (Vertex u, Vertex *BasePtr)` |
| `FindPathToVertex_s` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:11092` | external | tree-sitter+gcc-aux | `extern int FindPathToVertex_s (Vertex x, Edge (*SwitchEdge), Vertex *BasePtr, Vertex *Path, int MaxPathLen)` |
| `MakeBlossom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:11122` | external | tree-sitter+gcc-aux | `extern Vertex MakeBlossom (BN_STRUCT *pBNS, Vertex *ScanQ, int *pQSize, Vertex *Pu, Vertex *Pv, int max_len_Pu_Pv, Edge (*SwitchEdge), Vertex *BasePtr, Vertex u, Vertex v, EdgeIndex iuv, Vertex b_u, Vertex b_v, S_CHAR *Tree)` |
| `PullFlow` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:11290` | external | tree-sitter+gcc-aux | `extern int PullFlow (BN_STRUCT *pBNS, Edge (*SwitchEdge), Vertex x, Vertex y, int delta, S_CHAR bReverse, int bChangeFlow)` |
| `FindPathCap` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:11367` | external | tree-sitter+gcc-aux | `extern int FindPathCap (BN_STRUCT *pBNS, Edge (*SwitchEdge), Vertex x, Vertex y, int delta)` |
| `MarkRingSystemsAltBns` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:11459` | external | tree-sitter+gcc-aux | `extern int MarkRingSystemsAltBns (BN_STRUCT *pBNS, int bUnknAltAsNoStereo)` |
| `ReInitBnStructForAltBns` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:11800` | external | tree-sitter+gcc-aux | `extern int ReInitBnStructForAltBns (BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, int bUnknAltAsNoStereo)` |
| `MarkNonStereoAltBns` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:11909` | external | tree-sitter+gcc-aux | `extern int MarkNonStereoAltBns (BN_STRUCT *pBNS, inp_ATOM *at, int num_atoms, int bUnknAltAsNoStereo)` |
| `bHasChargedNeighbor` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:12005` | external | tree-sitter+gcc-aux | `extern int bHasChargedNeighbor (inp_ATOM *at, int iat)` |
| `AddRemoveProtonsRestr` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:12056` | external | tree-sitter+gcc-aux | `extern int AddRemoveProtonsRestr (inp_ATOM *at, int num_atoms, int *num_protons_to_add, int nNumProtAddedByRestr, INCHI_MODE bNormalizationFlags, int num_tg, int nChargeRevrs, int nChargeInChI)` |
| `AddRemoveIsoProtonsRestr` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_bns.c:12360` | external | tree-sitter+gcc-aux | `extern int AddRemoveIsoProtonsRestr (inp_ATOM *at, int num_atoms, NUM_H *num_protons_to_add, int num_tg)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c`

Parse errors: `4`. Function definitions: `39`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `inchi_ios_init` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:85` | external | tree-sitter+gcc-aux | `extern void inchi_ios_init (INCHI_IOSTREAM *ios, int io_type, FILE *f)` |
| `inchi_ios_create_copy` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:106` | external | tree-sitter+gcc-aux | `extern int inchi_ios_create_copy (INCHI_IOSTREAM *ios, INCHI_IOSTREAM *ios0)` |
| `inchi_ios_flush` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:144` | external | tree-sitter+gcc-aux | `extern void inchi_ios_flush (INCHI_IOSTREAM *ios)` |
| `inchi_ios_flush2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:183` | external | tree-sitter+gcc-aux | `extern void inchi_ios_flush2 (INCHI_IOSTREAM *ios, FILE *f2)` |
| `inchi_ios_close` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:229` | external | tree-sitter+gcc-aux | `extern void inchi_ios_close (INCHI_IOSTREAM *ios)` |
| `inchi_ios_reset` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:255` | external | tree-sitter+gcc-aux | `extern void inchi_ios_reset (INCHI_IOSTREAM *ios)` |
| `inchi_ios_free_str` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:272` | external | tree-sitter+gcc-aux | `extern void inchi_ios_free_str (INCHI_IOSTREAM *ios)` |
| `inchi_ios_str_getc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:294` | external | tree-sitter+gcc-aux | `extern int inchi_ios_str_getc (INCHI_IOSTREAM *ios)` |
| `inchi_ios_str_gets` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:324` | external | tree-sitter+gcc-aux | `extern char *inchi_ios_str_gets (char *szLine, int len, INCHI_IOSTREAM *f)` |
| `inchi_ios_str_getsTab` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:354` | external | tree-sitter+gcc-aux | `extern char *inchi_ios_str_getsTab (char *szLine, int len, INCHI_IOSTREAM *f)` |
| `inchi_ios_gets` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:386` | external | tree-sitter+gcc-aux | `extern int inchi_ios_gets (char *szLine, int len, INCHI_IOSTREAM *f, int *bTooLongLine)` |
| `inchi_ios_getsTab` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:420` | external | tree-sitter+gcc-aux | `extern int inchi_ios_getsTab (char *szLine, int len, INCHI_IOSTREAM *f, int *bTooLongLine)` |
| `inchi_ios_getsTab1` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:451` | external | tree-sitter+gcc-aux | `extern int inchi_ios_getsTab1 (char *szLine, int len, INCHI_IOSTREAM *f, int *bTooLongLine)` |
| `inchi_ios_print` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:477` | external | tree-sitter+gcc-aux | `extern int inchi_ios_print (INCHI_IOSTREAM *ios, const char *lpszFormat, ...)` |
| `push_to_winchi_text_window` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:576` | external | tree-sitter+gcc-aux | `extern int push_to_winchi_text_window (INCHI_IOSTREAM *ios)` |
| `inchi_ios_print_nodisplay` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:604` | external | tree-sitter+gcc-aux | `extern int inchi_ios_print_nodisplay (INCHI_IOSTREAM *ios, const char *lpszFormat, ...)` |
| `inchi_ios_flush_not_displayed` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:679` | external | tree-sitter+gcc-aux | `extern int inchi_ios_flush_not_displayed (INCHI_IOSTREAM *ios)` |
| `inchi_ios_eprint` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:708` | external | tree-sitter+gcc-aux | `extern int inchi_ios_eprint (INCHI_IOSTREAM *ios, const char *lpszFormat, ...)` |
| `inchi_fprintf` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:800` | external | tree-sitter+gcc-aux | `extern int inchi_fprintf (FILE *f, const char *lpszFormat, ...)` |
| `inchi_vfprintf` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:828` | external | tree-sitter+gcc-aux | `extern int inchi_vfprintf (FILE *f, const char *lpszFormat, __va_list_tag *argList)` |
| `inchi_print_nodisplay` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:867` | external | tree-sitter+gcc-aux | `extern int inchi_print_nodisplay (FILE *f, const char *lpszFormat, ...)` |
| `inchi_fgetsLfTab` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:893` | external | tree-sitter+gcc-aux | `extern int inchi_fgetsLfTab (char *szLine, int len, FILE *f)` |
| `inchi_fgetsLfTab` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:926` | external | tree-sitter | `int inchi_fgetsLfTab(char* szLine, int len, FILE* f)` |
| `inchi_fgetsTab` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:966` | external | tree-sitter+gcc-aux | `extern char *inchi_fgetsTab (char *szLine, int len, FILE *f)` |
| `inchi_fgetsLf` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:996` | external | tree-sitter+gcc-aux | `extern char *inchi_fgetsLf (char *line, int line_len, INCHI_IOSTREAM *inp_stream)` |
| `GetMaxPrintfLength` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1065` | external | tree-sitter+gcc-aux | `extern int GetMaxPrintfLength (const char *lpszFormat, __va_list_tag *argList)` |
| `inchi_sgets` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1322` | external | tree-sitter+gcc-aux | `extern char *inchi_sgets (char *s, int n, INCHI_IOSTREAM *ios)` |
| `inchi_strbuf_init` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1370` | external | tree-sitter+gcc-aux | `extern int inchi_strbuf_init (INCHI_IOS_STRING *buf, int start_size, int incr_size)` |
| `inchi_strbuf_reset` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1403` | external | tree-sitter+gcc-aux | `extern void inchi_strbuf_reset (INCHI_IOS_STRING *buf)` |
| `inchi_strbuf_close` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1422` | external | tree-sitter+gcc-aux | `extern void inchi_strbuf_close (INCHI_IOS_STRING *buf)` |
| `inchi_strbuf_create_copy` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1438` | external | tree-sitter+gcc-aux | `extern int inchi_strbuf_create_copy (INCHI_IOS_STRING *buf2, INCHI_IOS_STRING *buf)` |
| `inchi_strbuf_update` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1459` | external | tree-sitter+gcc-aux | `extern int inchi_strbuf_update (INCHI_IOS_STRING *buf, int new_addition_size)` |
| `inchi_strbuf_printf` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1507` | external | tree-sitter+gcc-aux | `extern int inchi_strbuf_printf (INCHI_IOS_STRING *buf, const char *lpszFormat, ...)` |
| `inchi_strbuf_printf_from` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1544` | external | tree-sitter+gcc-aux | `extern int inchi_strbuf_printf_from (INCHI_IOS_STRING *buf, int npos, const char *lpszFormat, ...)` |
| `inchi_strbuf_getline` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1584` | external | tree-sitter+gcc-aux | `extern int inchi_strbuf_getline (INCHI_IOS_STRING *buf, FILE *f, int crlf2lf, int preserve_lf)` |
| `inchi_strbuf_addline` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1635` | external | tree-sitter+gcc-aux | `extern int inchi_strbuf_addline (INCHI_IOS_STRING *buf, INCHI_IOSTREAM *inp_stream, int crlf2lf, int preserve_lf)` |
| `_inchi_trace` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1682` | external | tree-sitter | `int _inchi_trace(char* format, ...)` |
| `_inchi_trace` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1699` | external | tree-sitter+gcc-aux | `extern int _inchi_trace (char *format, ...)` |
| `Output_RecordInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichi_io.c:1709` | external | tree-sitter+gcc-aux | `extern int Output_RecordInfo (INCHI_IOSTREAM *out_file, int num_input_struct, int bNoStructLabels, const char *szSdfLabel, const char *szSdfValue, long unsigned int lSdfId, char *pLF, char *pTAB)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c`

Parse errors: `56`. Function definitions: `84`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `CanonGraph01` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:348` | external | tree-sitter | `int CanonGraph01( int n, int n_tg, int n_max, int bDigraph, Graph *G, Partition pi[], AT_RANK *nSymmRank, AT_RANK *nCanonRank, AT_NUMB *nAtomNumberCanon, CANON_DATA *pCD, CANON_COUNTS *pCC, ConTable **pp_zb_rho_inp, ConTable **pp_zb_rho_out, int LargeMolecules )` |
| `CanonGraph02` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:371` | external | tree-sitter | `int CanonGraph02( int n, int n_tg, int n_max, int bDigraph, Graph *G, Partition pi[], AT_RANK *nSymmRank, AT_RANK *nCanonRank, AT_NUMB *nAtomNumberCanon, CANON_DATA *pCD, CANON_COUNTS *pCC, ConTable **pp_zb_rho_inp, ConTable **pp_zb_rho_out, int LargeMolecules )` |
| `CanonGraph03` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:395` | external | tree-sitter | `int CanonGraph03( int n, int n_tg, int n_max, int bDigraph, Graph *G, Partition pi[], AT_RANK *nSymmRank, AT_RANK *nCanonRank, AT_NUMB *nAtomNumberCanon, CANON_DATA *pCD, CANON_COUNTS *pCC, ConTable **pp_zb_rho_inp, ConTable **pp_zb_rho_out, int LargeMolecules )` |
| `CanonGraph04` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:418` | external | tree-sitter | `int CanonGraph04( int n, int n_tg, int n_max, int bDigraph, Graph *G, Partition pi[], AT_RANK *nSymmRank, AT_RANK *nCanonRank, AT_NUMB *nAtomNumberCanon, CANON_DATA *pCD, CANON_COUNTS *pCC, ConTable **pp_zb_rho_inp, ConTable **pp_zb_rho_out, int LargeMolecules )` |
| `CanonGraph05` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:442` | external | tree-sitter | `int CanonGraph05( int n, int n_tg, int n_max, int bDigraph, Graph *G, Partition pi[], AT_RANK *nSymmRank, AT_RANK *nCanonRank, AT_NUMB *nAtomNumberCanon, CANON_DATA *pCD, CANON_COUNTS *pCC, ConTable **pp_zb_rho_inp, ConTable **pp_zb_rho_out, int LargeMolecules )` |
| `CanonGraph06` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:466` | external | tree-sitter | `int CanonGraph06( int n, int n_tg, int n_max, int bDigraph, Graph *G, Partition pi[], AT_RANK *nSymmRank, AT_RANK *nCanonRank, AT_NUMB *nAtomNumberCanon, CANON_DATA *pCD, CANON_COUNTS *pCC, ConTable **pp_zb_rho_inp, ConTable **pp_zb_rho_out, int LargeMolecules )` |
| `CanonGraph07` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:490` | external | tree-sitter | `int CanonGraph07( int n, int n_tg, int n_max, int bDigraph, Graph *G, Partition pi[], AT_RANK *nSymmRank, AT_RANK *nCanonRank, AT_NUMB *nAtomNumberCanon, CANON_DATA *pCD, CANON_COUNTS *pCC, ConTable **pp_zb_rho_inp, ConTable **pp_zb_rho_out, int LargeMolecules )` |
| `CanonGraph08` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:514` | external | tree-sitter | `int CanonGraph08( int n, int n_tg, int n_max, int bDigraph, Graph *G, Partition pi[], AT_RANK *nSymmRank, AT_RANK *nCanonRank, AT_NUMB *nAtomNumberCanon, CANON_DATA *pCD, CANON_COUNTS *pCC, ConTable **pp_zb_rho_inp, ConTable **pp_zb_rho_out, int LargeMolecules )` |
| `CanonGraph09` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:538` | external | tree-sitter | `int CanonGraph09( int n, int n_tg, int n_max, int bDigraph, Graph *G, Partition pi[], AT_RANK *nSymmRank, AT_RANK *nCanonRank, AT_NUMB *nAtomNumberCanon, CANON_DATA *pCD, CANON_COUNTS *pCC, ConTable **pp_zb_rho_inp, ConTable **pp_zb_rho_out, int LargeMolecules )` |
| `CanonGraph10` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:562` | external | tree-sitter | `int CanonGraph10( int n, int n_tg, int n_max, int bDigraph, Graph *G, Partition pi[], AT_RANK *nSymmRank, AT_RANK *nCanonRank, AT_NUMB *nAtomNumberCanon, CANON_DATA *pCD, CANON_COUNTS *pCC, ConTable **pp_zb_rho_inp, ConTable **pp_zb_rho_out, int LargeMolecules )` |
| `CanonGraph11` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:586` | external | tree-sitter | `int CanonGraph11( int n, int n_tg, int n_max, int bDigraph, Graph *G, Partition pi[], AT_RANK *nSymmRank, AT_RANK *nCanonRank, AT_NUMB *nAtomNumberCanon, CANON_DATA *pCD, CANON_COUNTS *pCC, ConTable **pp_zb_rho_inp, ConTable **pp_zb_rho_out, int LargeMolecules )` |
| `CanonGraph12` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:610` | external | tree-sitter | `int CanonGraph12( int n, int n_tg, int n_max, int bDigraph, Graph *G, Partition pi[], AT_RANK *nSymmRank, AT_RANK *nCanonRank, AT_NUMB *nAtomNumberCanon, CANON_DATA *pCD, CANON_COUNTS *pCC, ConTable **pp_zb_rho_inp, ConTable **pp_zb_rho_out, int LargeMolecules )` |
| `fill_crc32_data` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:657` | external | tree-sitter | `void fill_crc32_data( )` |
| `add2crc32` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:676` | external | tree-sitter | `unsigned long add2crc32( unsigned long crc32, AT_NUMB n )` |
| `TranspositionCreate` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:694` | external | tree-sitter+gcc-aux | `extern int TranspositionCreate (Transposition *p, int n)` |
| `TranspositionFree` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:707` | external | tree-sitter+gcc-aux | `extern void TranspositionFree (Transposition *p)` |
| `NodeSetCreate` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:718` | external | tree-sitter+gcc-aux | `extern int NodeSetCreate (struct tagCANON_GLOBALS *pCG, NodeSet *pSet, int n, int L)` |
| `NodeSetFree` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:754` | external | tree-sitter+gcc-aux | `extern void NodeSetFree (struct tagCANON_GLOBALS *pCG, NodeSet *pSet)` |
| `CTableCreate` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:769` | external | tree-sitter+gcc-aux | `extern int CTableCreate (ConTable *Ct, int n, CANON_DATA *pCD)` |
| `CTableFree` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:851` | external | tree-sitter+gcc-aux | `extern void CTableFree (ConTable *Ct)` |
| `UnorderedPartitionCreate` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:904` | external | tree-sitter+gcc-aux | `extern int UnorderedPartitionCreate (UnorderedPartition *p, int n)` |
| `UnorderedPartitionFree` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:919` | external | tree-sitter+gcc-aux | `extern void UnorderedPartitionFree (UnorderedPartition *p)` |
| `UnorderedPartitionMakeDiscrete` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:933` | external | tree-sitter+gcc-aux | `extern void UnorderedPartitionMakeDiscrete (UnorderedPartition *p, int n)` |
| `PartitionCreate` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:946` | external | tree-sitter+gcc-aux | `extern int PartitionCreate (Partition *p, int n)` |
| `PartitionFree` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:959` | external | tree-sitter+gcc-aux | `extern void PartitionFree (Partition *p)` |
| `PartitionIsDiscrete` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:978` | external | tree-sitter+gcc-aux | `extern int PartitionIsDiscrete (Partition *p, int n)` |
| `PartitionGetFirstCell` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:998` | external | tree-sitter+gcc-aux | `extern int PartitionGetFirstCell (Partition *p, Cell *baseW, int k, int n)` |
| `CellMakeEmpty` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1040` | external | tree-sitter+gcc-aux | `extern void CellMakeEmpty (Cell *baseW, int k)` |
| `NodeSetFromVertices` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1052` | external | tree-sitter+gcc-aux | `extern void NodeSetFromVertices (CANON_GLOBALS *pCG, NodeSet *cur_nodes, int l, Node *v, int num_v)` |
| `AllNodesAreInSet` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1071` | external | tree-sitter+gcc-aux | `extern int AllNodesAreInSet (NodeSet *cur_nodes, int lcur_nodes, NodeSet *set, int lset)` |
| `PartitionGetMcrAndFixSet` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1098` | external | tree-sitter+gcc-aux | `extern void PartitionGetMcrAndFixSet (CANON_GLOBALS *pCG, Partition *p, NodeSet *Mcr, NodeSet *Fix, int n, int l)` |
| `NodeSetFromRadEndpoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1138` | external | tree-sitter+gcc-aux | `extern void NodeSetFromRadEndpoints (CANON_GLOBALS *pCG, NodeSet *cur_nodes, int k, Vertex *RadEndpoints, int num_v)` |
| `RemoveFromNodeSet` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1159` | external | tree-sitter+gcc-aux | `extern void RemoveFromNodeSet (CANON_GLOBALS *pCG, NodeSet *cur_nodes, int k, Vertex *v, int num_v)` |
| `DoNodeSetsIntersect` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1178` | external | tree-sitter+gcc-aux | `extern int DoNodeSetsIntersect (NodeSet *cur_nodes, int k1, int k2)` |
| `IsNodeSetEmpty` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1201` | external | tree-sitter+gcc-aux | `extern int IsNodeSetEmpty (NodeSet *cur_nodes, int k)` |
| `AddNodeSet2ToNodeSet1` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1223` | external | tree-sitter+gcc-aux | `extern void AddNodeSet2ToNodeSet1 (NodeSet *cur_nodes, int k1, int k2)` |
| `AddNodesToRadEndpoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1241` | external | tree-sitter+gcc-aux | `extern int AddNodesToRadEndpoints (CANON_GLOBALS *pCG, NodeSet *cur_nodes, int k, Vertex *RadEndpoints, Vertex vRad, int nStart, int nLen)` |
| `PartitionGetTransposition` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1285` | external | tree-sitter+gcc-aux | `extern void PartitionGetTransposition (Partition *pFrom, Partition *pTo, int n, Transposition *gamma)` |
| `nGetMcr2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1305` | external | tree-sitter+gcc-aux | `extern AT_RANK nGetMcr2 (AT_RANK *nEqArray, AT_RANK n)` |
| `nJoin2Mcrs2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1342` | external | tree-sitter+gcc-aux | `extern int nJoin2Mcrs2 (AT_RANK *nEqArray, AT_RANK n1, AT_RANK n2)` |
| `GetUnorderedPartitionMcrNode` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1372` | external | tree-sitter+gcc-aux | `extern Node GetUnorderedPartitionMcrNode (UnorderedPartition *p1, Node v)` |
| `UnorderedPartitionJoin` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1385` | external | tree-sitter+gcc-aux | `extern int UnorderedPartitionJoin (UnorderedPartition *p1, UnorderedPartition *p2, int n)` |
| `PartitionSatisfiesLemma_2_25` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1407` | external | tree-sitter+gcc-aux | `extern int PartitionSatisfiesLemma_2_25 (Partition *p, int n)` |
| `PartitionCopy` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1445` | external | tree-sitter+gcc-aux | `extern void PartitionCopy (Partition *To, Partition *From, int n)` |
| `PartitionColorVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1464` | external | tree-sitter+gcc-aux | `extern int PartitionColorVertex (CANON_GLOBALS *pCG, Graph *G, Partition *p, Node v, int n, int n_tg, int n_max, int bDigraph, int nNumPrevRanks)` |
| `CellGetMinNode` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1586` | external | tree-sitter+gcc-aux | `extern Node CellGetMinNode (Partition *p, Cell *W, Node v, CANON_DATA *pCD)` |
| `CellGetNumberOfNodes` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1748` | external | tree-sitter+gcc-aux | `extern int CellGetNumberOfNodes (Partition *p, Cell *W)` |
| `CellIntersectWithSet` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1768` | external | tree-sitter+gcc-aux | `extern int CellIntersectWithSet (CANON_GLOBALS *pCG, Partition *p, Cell *W, NodeSet *Mcr, int l)` |
| `CtPartClear` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1798` | external | tree-sitter+gcc-aux | `extern void CtPartClear (ConTable *Ct, int k)` |
| `insertions_sort_NeighList_AT_NUMBERS2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1819` | external | tree-sitter+gcc-aux | `extern void insertions_sort_NeighList_AT_NUMBERS2 (NEIGH_LIST base, AT_RANK *nRank, AT_RANK max_rj)` |
| `CtPartFill` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1849` | external | tree-sitter+gcc-aux | `extern void CtPartFill (Graph *G, CANON_DATA *pCD, Partition *p, ConTable *Ct, int k, int n, int n_tg, int n_max)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1934` | external | tree-sitter | `INCHI_HEAPCHK /****************** Well-defined part of base hydrogen atoms *******************/ if (Ct)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1963` | external | tree-sitter | `INCHI_HEAPCHK /****************** Well-defined part of fixed hydrogen atoms *******************/ if (pCD->NumHfixed && Ct->NumHfixed)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1982` | external | tree-sitter | `INCHI_HEAPCHK /****************** Well-defined part of isotopic keys ***************************/ if (pCD->iso_sort_key && Ct->iso_sort_key)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:1998` | external | tree-sitter | `INCHI_HEAPCHK /****************** Well-defined part of isotopic iso_exchg_atnos ***************************/ if (pCD->iso_exchg_atnos && Ct->iso_exchg_atnos)` |
| `CtPartINCHI_CANON_INFINITY` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:2065` | external | tree-sitter+gcc-aux | `extern void CtPartINCHI_CANON_INFINITY (ConTable *Ct, S_CHAR *cmp, int k)` |
| `CtPartCompare` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:2123` | external | tree-sitter+gcc-aux | `extern int CtPartCompare (ConTable *Ct1, ConTable *Ct2, S_CHAR *cmp, kLeast *kLeastForLayer, int k, int bOnlyCommon, int bSplitTautCompare)` |
| `CtFullCompare` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:2623` | external | tree-sitter+gcc-aux | `extern int CtFullCompare (ConTable *Ct1, ConTable *Ct2, int bOnlyCommon, int bSplitTautCompare)` |
| `CtFullCompareLayers` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:2863` | external | tree-sitter+gcc-aux | `extern int CtFullCompareLayers (kLeast *kLeastForLayer)` |
| `CtCompareLayersGetFirstDiff` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:2880` | external | tree-sitter+gcc-aux | `extern int CtCompareLayersGetFirstDiff (kLeast *kLeast_rho, int nOneAdditionalLayer, int *L_rho, int *I_rho, int *k_rho)` |
| `CtPartCompareLayers` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:2929` | external | tree-sitter+gcc-aux | `extern int CtPartCompareLayers (kLeast *kLeast_rho, int L_rho_fix_prev, int nOneAdditionalLayer)` |
| `UpdateCompareLayers` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:2947` | external | tree-sitter+gcc-aux | `extern void UpdateCompareLayers (kLeast *kLeastForLayer, int hzz)` |
| `CtPartCopy` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:2965` | external | tree-sitter+gcc-aux | `extern void CtPartCopy (ConTable *Ct1, ConTable *Ct2, int k)` |
| `CtFullCopy` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:3113` | external | tree-sitter+gcc-aux | `extern void CtFullCopy (ConTable *Ct1, ConTable *Ct2)` |
| `TranspositionGetMcrAndFixSetAndUnorderedPartition` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:3125` | external | tree-sitter+gcc-aux | `extern void TranspositionGetMcrAndFixSetAndUnorderedPartition (CANON_GLOBALS *pCG, Transposition *gamma, NodeSet *McrSet, NodeSet *FixSet, int n, int l, UnorderedPartition *p)` |
| `SetBitCreate` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:3205` | external | tree-sitter+gcc-aux | `extern int SetBitCreate (CANON_GLOBALS *pCG)` |
| `SetBitFree` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:3259` | external | tree-sitter+gcc-aux | `extern int SetBitFree (CANON_GLOBALS *pCG)` |
| `RenumberTheGraph` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:3278` | external | tree-sitter | `void RenumberTheGraph( int n, NEIGH_LIST *NeighList, AT_NUMB *old2new, AT_NUMB *new2old, S_CHAR *mark, int bDone )` |
| `RearrangeAtRankArray` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:3356` | external | tree-sitter | `void RearrangeAtRankArray( int n, AT_RANK *nRank, AT_NUMB *new2old, S_CHAR *mark, int bDone )` |
| `RenumberAtNumbArray` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:3394` | external | tree-sitter | `void RenumberAtNumbArray( int n, AT_NUMB *nAtNumb, AT_NUMB *old2new )` |
| `GetCanonRanking2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:3405` | external | tree-sitter | `int GetCanonRanking2( int num_atoms, int num_at_tg, int num_max, int bDigraph, sp_ATOM* at, AT_RANK **pRankStack, int nNumPrevRanks, AT_RANK *nSymmRank, AT_RANK *nCanonRank, NEIGH_LIST *NeighList, AT_RANK *nTempRank, CANON_STAT* pCS )` |
| `GetOneAdditionalLayer` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:3566` | external | tree-sitter+gcc-aux | `extern int GetOneAdditionalLayer (CANON_DATA *pCD, ConTable *pzb_rho_fix)` |
| `CanonGraph` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:3662` | external | tree-sitter+gcc-aux | `extern int CanonGraph (INCHI_CLOCK *ic, CANON_GLOBALS *pCG, int n, int n_tg, int n_max, int bDigraph, Graph *G, Partition *pi, AT_RANK *nSymmRank, AT_RANK *nCanonRank, AT_NUMB *nAtomNumberCanon, CANON_DATA *pCD, CANON_COUNTS *pCC, ConTable **pp_zb_rho_inp, ConTable **pp_zb_rho_out, int LargeMolecules)` |
| `while` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:3920` | external | tree-sitter | `INCHI_HEAPCHK /* L2: reach the first leaf and save it in zeta and rho */ while (k)` |
| `SetInitialRanks2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:4674` | external | tree-sitter+gcc-aux | `extern int SetInitialRanks2 (int num_atoms, ATOM_INVARIANT2 *pAtomInvariant2, AT_RANK *nNewRank, AT_RANK *nAtomNumber, CANON_GLOBALS *pCG)` |
| `FillOutAtomInvariant2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:4719` | external | tree-sitter+gcc-aux | `extern void FillOutAtomInvariant2 (sp_ATOM *at, int num_atoms, int num_at_tg, ATOM_INVARIANT2 *pAtomInvariant, int bIgnoreIsotopic, int bHydrogensInRanks, int bHydrogensFixedInRanks, int bDigraph, int bTautGroupsOnly, T_GROUP_INFO *t_group_info)` |
| `CleanNumH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:4881` | external | tree-sitter+gcc-aux | `extern void CleanNumH (NUM_H *NumH, int len)` |
| `CleanCt` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:4902` | external | tree-sitter+gcc-aux | `extern int CleanCt (AT_RANK *Ct, int len)` |
| `CleanIsoSortKeys` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:4915` | external | tree-sitter+gcc-aux | `extern void CleanIsoSortKeys (AT_ISO_SORT_KEY *isk, int len)` |
| `MergeCleanIsoSortKeys` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:4934` | external | tree-sitter | `void MergeCleanIsoSortKeys( AT_ISO_SORT_KEY * isk1, AT_ISO_SORT_KEY * isk2, int len )` |
| `DeAllocBCN` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:4960` | external | tree-sitter+gcc-aux | `extern void DeAllocBCN (BCN *pBCN)` |
| `bCanonIsFinerThanEquitablePartition` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:5011` | external | tree-sitter | `int bCanonIsFinerThanEquitablePartition( int num_atoms, sp_ATOM* at, AT_RANK *nSymmRank )` |
| `GetBaseCanonRanking` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:5105` | external | gcc-aux | `extern int GetBaseCanonRanking (INCHI_CLOCK *ic, int num_atoms, int num_at_tg, sp_ATOM **at, T_GROUP_INFO *t_group_info, ATOM_SIZES *s, BCN *pBCN, struct tagInchiTime *ulTimeOutTime, CANON_GLOBALS *pCG, int bFixIsoFixedH, int LargeMolecules)` |
| `FREE_CONTABLE` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c:6836` | external | tree-sitter | `FREE_CONTABLE( Ct_NoH ) FREE_CONTABLE( Ct_NoTautH ) FREE_CONTABLE( Ct_Base ) FREE_CONTABLE( Ct_FixH ) FREE_CONTABLE( Ct_Temp ) /* isotopic canonicalization */ FREE_CONTABLE( Ct_NoTautHIso ) FREE_CONTABLE( Ct_BaseIso ) FREE_CONTABLE( Ct_FixHIso ) /* free the first two pointers from pBCN->pRankStack */ FREE_ARRAY( nRank ) FREE_ARRAY( nAtomNumber ) if (pBCN->pRankStack)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c`

Parse errors: `5`. Function definitions: `22`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `InchiClock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:100` | internal | tree-sitter | `static clock_t InchiClock( void )` |
| `InchiClock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:111` | internal | tree-sitter+gcc-aux | `static clock_t InchiClock (void)` |
| `FillMaxMinClock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:128` | internal | tree-sitter+gcc-aux | `static void FillMaxMinClock (INCHI_CLOCK *ic)` |
| `InchiTimeGet` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:146` | external | tree-sitter+gcc-aux | `extern void InchiTimeGet (inchiTime *TickEnd)` |
| `InchiTimeMsecDiff` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:151` | external | tree-sitter+gcc-aux | `extern long int InchiTimeMsecDiff (INCHI_CLOCK *ic, inchiTime *TickEnd, inchiTime *TickStart)` |
| `InchiTimeElapsed` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:223` | external | tree-sitter+gcc-aux | `extern long int InchiTimeElapsed (INCHI_CLOCK *ic, inchiTime *TickStart)` |
| `InchiTimeAddMsec` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:234` | external | tree-sitter+gcc-aux | `extern void InchiTimeAddMsec (INCHI_CLOCK *ic, inchiTime *TickEnd, long unsigned int nNumMsec)` |
| `bInchiTimeIsOver` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:257` | external | tree-sitter+gcc-aux | `extern int bInchiTimeIsOver (INCHI_CLOCK *ic, inchiTime *TickStart)` |
| `InchiTimeGet` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:328` | external | tree-sitter | `void InchiTimeGet( inchiTime *TickEnd )` |
| `InchiTimeMsecDiff` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:341` | external | tree-sitter | `long InchiTimeMsecDiff( INCHI_CLOCK *ic, inchiTime *TickEnd, inchiTime *TickStart )` |
| `InchiTimeElapsed` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:365` | external | tree-sitter | `long InchiTimeElapsed( INCHI_CLOCK *ic, inchiTime *TickStart )` |
| `InchiTimeAddMsec` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:378` | external | tree-sitter | `void InchiTimeAddMsec( INCHI_CLOCK *ic, inchiTime *TickEnd, unsigned long nNumMsec )` |
| `bInchiTimeIsOver` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:393` | external | tree-sitter | `int bInchiTimeIsOver( INCHI_CLOCK *ic, inchiTime *TickEnd )` |
| `GetCanonLengths` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:418` | external | tree-sitter+gcc-aux | `extern int GetCanonLengths (int num_at, sp_ATOM *at, ATOM_SIZES *s, T_GROUP_INFO *t_group_info)` |
| `DeAllocateCS` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:491` | external | tree-sitter+gcc-aux | `extern int DeAllocateCS (CANON_STAT *pCS)` |
| `AllocateCS` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:575` | external | tree-sitter+gcc-aux | `extern int AllocateCS (CANON_STAT *pCS, int num_at, int num_at_tg, int nLenCT, int nLenCTAtOnly, int nLenLinearCTStereoDble, int nLenLinearCTIsotopicStereoDble, int nLenLinearCTStereoCarb, int nLenLinearCTIsotopicStereoCarb, int nLenLinearCTTautomer, int nLenLinearCTIsotopicTautomer, int nLenIsotopic, INCHI_MODE nMode, BCN *pBCN)` |
| `FillIsotopicAtLinearCT` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:780` | external | tree-sitter+gcc-aux | `extern int FillIsotopicAtLinearCT (int num_atoms, sp_ATOM *at, const AT_RANK *nAtomNumber, AT_ISOTOPIC *LinearCTIsotopic, int nMaxLenLinearCTIsotopic, int *pnLenLinearCTIsotopic)` |
| `FillTautLinearCT2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:858` | external | tree-sitter+gcc-aux | `extern int FillTautLinearCT2 (CANON_GLOBALS *pCG, int num_atoms, int num_at_tg, int bIsoTaut, const AT_RANK *nRank, const AT_RANK *nAtomNumber, const AT_RANK *nSymmRank, const AT_RANK *nRankIso, const AT_RANK *nAtomNumberIso, const AT_RANK *nSymmRankIso, AT_TAUTOMER *LinearCTTautomer, int nMaxLenLinearCTTautomer, int *pnLenLinearCTTautomer, AT_ISO_TGROUP *LinearCTIsotopicTautomer, int nMaxLenLinearCTIsotopicTautomer, int *pnLenLinearCTIsotopicTautomer, T_GROUP_INFO *t_group_info)` |
| `UpdateFullLinearCT` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:1050` | external | tree-sitter+gcc-aux | `extern int UpdateFullLinearCT (int num_atoms, int num_at_tg, sp_ATOM *at, AT_RANK *nRank, AT_RANK *nAtomNumber, CANON_STAT *pCS, CANON_GLOBALS *pCG, int bFirstTime)` |
| `FixCanonEquivalenceInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:1278` | external | tree-sitter+gcc-aux | `extern int FixCanonEquivalenceInfo (CANON_GLOBALS *pCG, int num_at_tg, AT_RANK *nSymmRank, AT_RANK *nCurrRank, AT_RANK *nTempRank, AT_NUMB *nAtomNumber, int *bChanged)` |
| `Canon_INChI3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:1370` | external | tree-sitter+gcc-aux | `extern int Canon_INChI3 (INCHI_CLOCK *ic, int num_atoms, int num_at_tg, sp_ATOM *at, CANON_STAT *pCS, CANON_GLOBALS *pCG, INCHI_MODE nMode, int bTautFtcn)` |
| `Canon_INChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicano.c:2727` | external | tree-sitter+gcc-aux | `extern int Canon_INChI (INCHI_CLOCK *ic, int num_atoms, int num_at_tg, sp_ATOM *at, CANON_STAT *pCS, CANON_GLOBALS *pCG, INCHI_MODE nMode, int bTautFtcn)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c`

Parse errors: `3`. Function definitions: `27`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `find_atoms_with_parity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:155` | external | tree-sitter+gcc-aux | `extern int find_atoms_with_parity (sp_ATOM *at, S_CHAR *visited, int from_atom, int cur_atom)` |
| `SetHalfStereoBondIllDefPariy` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:189` | external | tree-sitter+gcc-aux | `extern int SetHalfStereoBondIllDefPariy (sp_ATOM *at, int jn, int k1, int new_parity)` |
| `RemoveHalfStereoBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:208` | external | tree-sitter+gcc-aux | `extern int RemoveHalfStereoBond (sp_ATOM *at, int jn, int k1)` |
| `SetOneStereoBondIllDefParity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:242` | external | tree-sitter+gcc-aux | `extern int SetOneStereoBondIllDefParity (sp_ATOM *at, int jc, int k, int new_parity)` |
| `RemoveOneStereoBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:271` | external | tree-sitter+gcc-aux | `extern int RemoveOneStereoBond (sp_ATOM *at, int jc, int k)` |
| `RemoveOneStereoCenter` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:300` | external | tree-sitter+gcc-aux | `extern int RemoveOneStereoCenter (sp_ATOM *at, int jc)` |
| `UnmarkNonStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:322` | external | tree-sitter+gcc-aux | `extern int UnmarkNonStereo (CANON_GLOBALS *pCG, sp_ATOM *at, int num_atoms, const AT_RANK *nRank, const AT_RANK *nAtomNumber, int bIsotopic)` |
| `FillSingleStereoDescriptors` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:525` | external | tree-sitter+gcc-aux | `extern int FillSingleStereoDescriptors (CANON_GLOBALS *pCG, sp_ATOM *at, int i, int num_trans, const AT_RANK *nRank, AT_STEREO_CARB *LinearCTStereoCarb, int *nStereoCarbLen, int nMaxStereoCarbLen, AT_STEREO_DBLE *LinearCTStereoDble, int *nStereoDbleLen, int nMaxStereoDbleLen, int bAllene)` |
| `SwitchAtomStereoAndIsotopicStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:705` | external | tree-sitter+gcc-aux | `extern void SwitchAtomStereoAndIsotopicStereo (sp_ATOM *at, int num_atoms, int *bSwitched)` |
| `SetCtToIsotopicStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:729` | external | tree-sitter+gcc-aux | `extern void SetCtToIsotopicStereo (CANON_STAT *pCS, CANON_STAT *pCS2)` |
| `SetCtToNonIsotopicStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:746` | external | tree-sitter+gcc-aux | `extern void SetCtToNonIsotopicStereo (CANON_STAT *pCS, CANON_STAT *pCS2)` |
| `FillAllStereoDescriptors` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:766` | external | tree-sitter+gcc-aux | `extern int FillAllStereoDescriptors (CANON_GLOBALS *pCG, sp_ATOM *at, int num_atoms, const AT_RANK *nCanonRank, const AT_RANK *nAtomNumberCanon, CANON_STAT *pCS)` |
| `SetKnownStereoBondParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:801` | external | tree-sitter+gcc-aux | `extern int SetKnownStereoBondParities (CANON_GLOBALS *pCG, sp_ATOM *at, int num_atoms, const AT_RANK *nCanonRank, const AT_RANK *nRank, const AT_RANK *nAtomNumber)` |
| `MarkKnownEqualStereoBondParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:1093` | external | tree-sitter+gcc-aux | `extern int MarkKnownEqualStereoBondParities (sp_ATOM *at, int num_atoms, const AT_RANK *nRank, const AT_RANK *nAtomNumber)` |
| `GetNextNeighborAndRank` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:1419` | external | tree-sitter+gcc-aux | `extern int GetNextNeighborAndRank (sp_ATOM *at, AT_RANK cur, AT_RANK prev, AT_RANK *n, AT_RANK *cr, const AT_RANK *nCanonRank)` |
| `GetAndCheckNextNeighbors` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:1456` | external | tree-sitter+gcc-aux | `extern int GetAndCheckNextNeighbors (sp_ATOM *at, AT_RANK cur1, AT_RANK prev1, AT_RANK cur2, AT_RANK prev2, AT_RANK *n1, AT_RANK *n2, AT_RANK *nVisited1, AT_RANK *nVisited2, const AT_RANK *nRank, const AT_RANK *nCanonRank)` |
| `PathsHaveIdenticalKnownParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:1530` | external | tree-sitter+gcc-aux | `extern AT_RANK PathsHaveIdenticalKnownParities (sp_ATOM *at, AT_RANK prev1, AT_RANK cur1, AT_RANK prev2, AT_RANK cur2, AT_RANK *nVisited1, AT_RANK *nVisited2, const AT_RANK *nRank, const AT_RANK *nCanonRank, AT_RANK nLength)` |
| `RemoveKnownNonStereoBondParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:1608` | external | tree-sitter+gcc-aux | `extern int RemoveKnownNonStereoBondParities (sp_ATOM *at, int num_atoms, const AT_RANK *nCanonRank, const AT_RANK *nRank, CANON_STAT *pCS)` |
| `SetKnownStereoCenterParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:1714` | external | tree-sitter+gcc-aux | `extern int SetKnownStereoCenterParities (CANON_GLOBALS *pCG, sp_ATOM *at, int num_atoms, const AT_RANK *nCanonRank, const AT_RANK *nRank, const AT_RANK *nAtomNumber)` |
| `RemoveKnownNonStereoCenterParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:1833` | external | tree-sitter+gcc-aux | `extern int RemoveKnownNonStereoCenterParities (CANON_GLOBALS *pCG, sp_ATOM *at, int num_atoms, const AT_RANK *nCanonRank, const AT_RANK *nRank, CANON_STAT *pCS)` |
| `MarkKnownEqualStereoCenterParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:1942` | external | tree-sitter+gcc-aux | `extern int MarkKnownEqualStereoCenterParities (sp_ATOM *at, int num_atoms, const AT_RANK *nRank, const AT_RANK *nAtomNumber)` |
| `InvertStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:2016` | external | tree-sitter+gcc-aux | `extern int InvertStereo (sp_ATOM *at, int num_at_tg, AT_RANK *nCanonRank, AT_RANK *nAtomNumberCanon, CANON_STAT *pCS, int bInvertLinearCTStereo)` |
| `FillOutStereoParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:2127` | external | tree-sitter+gcc-aux | `extern int FillOutStereoParities (sp_ATOM *at, int num_atoms, const AT_RANK *nCanonRank, const AT_RANK *nAtomNumberCanon, const AT_RANK *nRank, const AT_RANK *nAtomNumber, CANON_STAT *pCS, CANON_GLOBALS *pCG, int bIsotopic)` |
| `GetStereoNeighborPos` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:2197` | external | tree-sitter+gcc-aux | `extern int GetStereoNeighborPos (sp_ATOM *at, int iAt1, int iAt2)` |
| `GetStereoBondParity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:2219` | external | tree-sitter+gcc-aux | `extern int GetStereoBondParity (sp_ATOM *at, int i, int n, AT_RANK *nRank)` |
| `GetPermutationParity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:2296` | external | tree-sitter+gcc-aux | `extern int GetPermutationParity (CANON_GLOBALS *pCG, sp_ATOM *at, AT_RANK nAvoidNeighbor, AT_RANK *nCanonRank)` |
| `GetStereoCenterParity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichicans.c:2340` | external | tree-sitter+gcc-aux | `extern int GetStereoCenterParity (CANON_GLOBALS *pCG, sp_ATOM *at, int i, AT_RANK *nRank)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichierr.c`

Parse errors: `0`. Function definitions: `3`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `ErrMsg` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichierr.c:54` | external | tree-sitter+gcc-aux | `extern const char *ErrMsg (int nErrorCode)` |
| `AddErrorMessage` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichierr.c:106` | external | tree-sitter+gcc-aux | `extern int AddErrorMessage (char *all_messages, const char *new_message)` |
| `already_have_this_message` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichierr.c:160` | internal | tree-sitter+gcc-aux | `static int already_have_this_message (char *prev_messages, const char *new_message)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiisot.c`

Parse errors: `0`. Function definitions: `3`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `make_iso_sort_key` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiisot.c:47` | external | tree-sitter+gcc-aux | `extern AT_ISO_SORT_KEY make_iso_sort_key (int iso_atw_diff, int num_1H, int num_2H, int num_3H)` |
| `set_atom_iso_sort_keys` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiisot.c:67` | external | tree-sitter+gcc-aux | `extern int set_atom_iso_sort_keys (int num_at, sp_ATOM *at, T_GROUP_INFO *t_group_info, int *bHasIsotopicInTautomerGroups)` |
| `unpack_iso_sort_key` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiisot.c:112` | external | tree-sitter | `int unpack_iso_sort_key( AT_ISO_SORT_KEY iso_sort_key, S_CHAR *num_1H, S_CHAR *num_2H, S_CHAR *num_3H, S_CHAR *iso_atw_diff )` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c`

Parse errors: `0`. Function definitions: `12`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `MakeHillFormulaString` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:121` | external | tree-sitter+gcc-aux | `extern int MakeHillFormulaString (char *szHillFormula, INCHI_IOS_STRING *strbuf, int *bOverflow)` |
| `GetHillFormulaIndexLength` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:144` | external | tree-sitter+gcc-aux | `extern int GetHillFormulaIndexLength (int count)` |
| `GetHillFormulaCounts` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:159` | external | tree-sitter+gcc-aux | `extern int GetHillFormulaCounts (U_CHAR *nAtom, S_CHAR *nNum_H, int num_atoms, AT_NUMB *nTautomer, int lenTautomer, int *pnum_C, int *pnum_H, int *pnLen, int *pnNumNonHAtoms)` |
| `AddElementAndCount` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:279` | external | tree-sitter+gcc-aux | `extern int AddElementAndCount (const char *szElement, int mult, char *szLinearCT, int nLenLinearCT, int *bOverflow)` |
| `MakeHillFormula` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:316` | external | tree-sitter+gcc-aux | `extern int MakeHillFormula (U_CHAR *nAtom, int num_atoms, char *szLinearCT, int nLen_szLinearCT, int num_C, int num_H, int *bOverflow)` |
| `AllocateAndFillHillFormula` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:402` | external | tree-sitter+gcc-aux | `extern char *AllocateAndFillHillFormula (INChI *pINChI)` |
| `Copy2StereoBondOrAllene` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:451` | external | tree-sitter+gcc-aux | `extern int Copy2StereoBondOrAllene (INChI_Stereo *Stereo, int *nNumberOfStereoCenters, int *nNumberOfStereoBonds, AT_STEREO_DBLE *LinearCTStereoDble, AT_NUMB *pCanonOrd, AT_RANK *pCanonRank, sp_ATOM *at, int bIsotopic)` |
| `CopyLinearCTStereoToINChIStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:566` | external | tree-sitter+gcc-aux | `extern int CopyLinearCTStereoToINChIStereo (INChI_Stereo *Stereo, AT_STEREO_CARB *LinearCTStereoCarb, int nLenLinearCTStereoCarb, AT_STEREO_DBLE *LinearCTStereoDble, int nLenLinearCTStereoDble, AT_NUMB *pCanonOrd, AT_RANK *pCanonRank, sp_ATOM *at, int bIsotopic, AT_STEREO_CARB *LinearCTStereoCarbInv, AT_STEREO_DBLE *LinearCTStereoDbleInv, AT_NUMB *pCanonOrdInv, AT_RANK *pCanonRankInv)` |
| `MarkAmbiguousStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:713` | external | tree-sitter+gcc-aux | `extern int MarkAmbiguousStereo (sp_ATOM *at, inp_ATOM *norm_at, int bIsotopic, AT_NUMB *pCanonOrd, AT_STEREO_CARB *LinearCTStereoCarb, int nLenLinearCTStereoCarb, AT_STEREO_DBLE *LinearCTStereoDble, int nLenLinearCTStereoDble)` |
| `UnmarkAllUndefinedUnknownStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:823` | external | tree-sitter+gcc-aux | `extern INCHI_MODE UnmarkAllUndefinedUnknownStereo (INChI_Stereo *Stereo, INCHI_MODE nUserMode)` |
| `WriteCoord` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:890` | external | tree-sitter+gcc-aux | `extern void WriteCoord (char *str, double x)` |
| `FillOutINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimak2.c:1072` | external | tree-sitter+gcc-aux | `extern int FillOutINChI (INChI *pINChI, INChI_Aux *pINChI_Aux, int num_atoms, int num_at_tg, int num_removed_H, sp_ATOM *at, inp_ATOM *norm_at, CANON_STAT *pCS, CANON_GLOBALS *pCG, int bTautomeric, INCHI_MODE nUserMode, char *pStrErrStruct, int bNoWarnings)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c`

Parse errors: `9`. Function definitions: `32`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `inp2spATOM` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:119` | external | tree-sitter+gcc-aux | `extern int inp2spATOM (inp_ATOM *inp_at, int num_inp_at, sp_ATOM *at)` |
| `GetElementAndCount` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:173` | external | tree-sitter+gcc-aux | `extern int GetElementAndCount (const char **f, char *szEl, int *count)` |
| `CompareHillFormulas` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:241` | external | tree-sitter+gcc-aux | `extern int CompareHillFormulas (const char *f1, const char *f2)` |
| `CompareHillFormulasNoH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:272` | external | tree-sitter+gcc-aux | `extern int CompareHillFormulasNoH (const char *f1, const char *f2, int *num_H1, int *num_H2)` |
| `CompareTautNonIsoPartOfINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:316` | external | tree-sitter+gcc-aux | `extern int CompareTautNonIsoPartOfINChI (const INChI *i1, const INChI *i2)` |
| `CompINChITautVsNonTaut` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:341` | external | tree-sitter+gcc-aux | `extern int CompINChITautVsNonTaut (const INCHI_SORT *p1, const INCHI_SORT *p2, int bCompareIsotopic)` |
| `GetSp3RelRacAbs` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:593` | external | tree-sitter+gcc-aux | `extern int GetSp3RelRacAbs (const INChI *pINChI, INChI_Stereo *Stereo)` |
| `CompINChILayers` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:645` | external | tree-sitter+gcc-aux | `extern int CompINChILayers (const INCHI_SORT *p1, const INCHI_SORT *p2, char (*sDifSegs)[11], int bFixTranspChargeBug)` |
| `INChI_SegmentAction` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:1486` | external | tree-sitter+gcc-aux | `extern int INChI_SegmentAction (char cDifSegs)` |
| `MarkUnusedAndEmptyLayers` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:1525` | external | tree-sitter+gcc-aux | `extern int MarkUnusedAndEmptyLayers (char (*sDifSegs)[11])` |
| `CompareInchiStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:1607` | external | tree-sitter+gcc-aux | `extern int CompareInchiStereo (INChI_Stereo *Stereo1, INCHI_MODE nFlags1, INChI_Stereo *Stereo2, INCHI_MODE nFlags2)` |
| `CompINChI2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:1712` | external | tree-sitter+gcc-aux | `extern int CompINChI2 (const INCHI_SORT *p1, const INCHI_SORT *p2, int bTaut, int bCompareIsotopic)` |
| `CompINChINonTaut2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2042` | external | tree-sitter+gcc-aux | `extern int CompINChINonTaut2 (const void *p1, const void *p2)` |
| `CompINChITaut2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2064` | external | tree-sitter+gcc-aux | `extern int CompINChITaut2 (const void *p1, const void *p2)` |
| `mystrrev` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2090` | external | tree-sitter+gcc-aux | `extern void mystrrev (char *p)` |
| `CompareDfsDescendants4CT` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2119` | internal | tree-sitter+gcc-aux | `static int CompareDfsDescendants4CT (const void *a1, const void *a2, void *p)` |
| `GetDfsOrder4CT` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2154` | external | tree-sitter+gcc-aux | `extern AT_NUMB *GetDfsOrder4CT (CANON_GLOBALS *pCG, AT_NUMB *LinearCT, int nLenCT, S_CHAR *nNum_H, int num_atoms, int nCtMode)` |
| `GetInpStructErrorType` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2480` | external | tree-sitter+gcc-aux | `extern int GetInpStructErrorType (INPUT_PARMS *ip, int err, char *pStrErrStruct, int num_inp_atoms)` |
| `ProcessStructError` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2511` | external | tree-sitter+gcc-aux | `extern int ProcessStructError (INCHI_IOSTREAM *out_file, INCHI_IOSTREAM *log_file, char *pStrErrStruct, int nErrorType, long int num_inp, INPUT_PARMS *ip)` |
| `CompareReversedStereoINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2555` | external | tree-sitter+gcc-aux | `extern int CompareReversedStereoINChI (INChI_Stereo *s1, INChI_Stereo *s2)` |
| `CompareReversedStereoINChI2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2635` | external | tree-sitter+gcc-aux | `extern int CompareReversedStereoINChI2 (INChI_Stereo *s1, INChI_Stereo *s2, ICR *picr)` |
| `CompareReversedINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:2936` | external | tree-sitter+gcc-aux | `extern int CompareReversedINChI (INChI *i1, INChI *i2, INChI_Aux *a1, INChI_Aux *a2)` |
| `CompareIcr` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:3189` | external | tree-sitter+gcc-aux | `extern int CompareIcr (ICR *picr1, ICR *picr2, INCHI_MODE *pin1, INCHI_MODE *pin2, INCHI_MODE mask)` |
| `CompareReversedINChI2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:3262` | external | tree-sitter+gcc-aux | `extern INCHI_MODE CompareReversedINChI2 (INChI *i1, INChI *i2, INChI_Aux *a1, INChI_Aux *a2, ICR *picr, int *err)` |
| `Create_INChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:3707` | external | tree-sitter+gcc-aux | `extern int Create_INChI (CANON_GLOBALS *pCG, INCHI_CLOCK *ic, INPUT_PARMS *ip, INChI **ppINChI, INChI_Aux **ppINChI_Aux, ORIG_ATOM_DATA *orig_inp_data, inp_ATOM *inp_at, INP_ATOM_DATA **out_norm_data, int num_inp_at, INCHI_MODE nUserMode, INCHI_MODE *pbTautFlags, INCHI_MODE *pbTautFlagsDone, struct tagInchiTime *ulMaxTime, T_GROUP_INFO *ti_out, char *pStrErrStruct)` |
| `GetAtomOrdNbrInCanonOrd` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:4569` | external | tree-sitter | `int GetAtomOrdNbrInCanonOrd(struct tagCANON_GLOBALS* pCG, inp_ATOM* norm_at, AT_NUMB* nAtomOrdNbr, AT_NUMB* nOrigAtNosInCanonOrd, int num_at)` |
| `FillOutCanonInfAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:4640` | external | tree-sitter | `int FillOutCanonInfAtom(struct tagCANON_GLOBALS* pCG, inp_ATOM* norm_at, INF_ATOM_DATA* inf_norm_at_data, int init_num_at, int bIsotopic, INChI* pINChI, INChI_Aux* pINChI_Aux, int bAbcNumbers, INCHI_MODE nMode)` |
| `FillOutOneCanonInfAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:5255` | external | tree-sitter | `int FillOutOneCanonInfAtom(struct tagCANON_GLOBALS* pCG, inp_ATOM* inp_norm_at, INF_ATOM_DATA* inf_norm_at_data, AT_NUMB* pStereoFlags, int init_num_at, int offset, int offset_H, int bIsotopic, INChI* pINChI, INChI_Aux* pINChI_Aux, int bAbcNumbers, INCHI_MODE nMode)` |
| `FillOutInputInfAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:5882` | external | tree-sitter | `int FillOutInputInfAtom(inp_ATOM* inp_at, INF_ATOM_DATA* inf_at_data, int init_num_at, int num_removed_H, int bAdd_DT_to_num_H, int nNumRemovedProtons, NUM_H* nNumRemovedProtonsIsotopic, int bIsotopic, int bAbcNumbers) /* djb-rwth: this whole function seems to be useless as it returns the value of a function argument -- init_num_at */` |
| `FillOutInfAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:6098` | external | tree-sitter | `int FillOutInfAtom(struct tagCANON_GLOBALS* pCG, inp_ATOM* norm_at, INF_ATOM_DATA* inf_norm_at_data, int init_num_at, int num_removed_H, int bAdd_DT_to_num_H, int nNumRemovedProtons, NUM_H* nNumRemovedProtonsIsotopic, int bIsotopic, INChI* pINChI, INChI_Aux* pINChI_Aux, int bAbcNumbers, INCHI_MODE nMode)` |
| `FillOutCompositeCanonInfAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:6131` | external | tree-sitter | `int FillOutCompositeCanonInfAtom(struct tagCANON_GLOBALS* pCG, COMP_ATOM_DATA* composite_norm_data, INF_ATOM_DATA* inf_norm_at_data, int bIsotopic, int bTautomeric, PINChI2* pINChI2, PINChI_Aux2* pINChI_Aux2, int bAbcNumbers, INCHI_MODE nMode)` |
| `CheckCanonNumberingCorrectness` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimake.c:6230` | external | tree-sitter+gcc-aux | `extern int CheckCanonNumberingCorrectness (int num_atoms, int num_at_tg, sp_ATOM *at, CANON_STAT *pCS, CANON_GLOBALS *pCG, int bTautomeric, char *pStrErrStruct)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c`

Parse errors: `5`. Function definitions: `30`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `All_SC_Same` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:53` | external | tree-sitter+gcc-aux | `extern int All_SC_Same (AT_RANK canon_rank1, const const ppAT_RANK pRankStack1, const const ppAT_RANK pRankStack2, const AT_RANK *nAtomNumberCanonFrom, const sp_ATOM *at)` |
| `Next_SC_At_CanonRank2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:99` | external | tree-sitter+gcc-aux | `extern int Next_SC_At_CanonRank2 (AT_RANK *canon_rank1, AT_RANK *canon_rank1_min, int *bFirstTime, S_CHAR *bAtomUsedForStereo, const const ppAT_RANK pRankStack1, const const ppAT_RANK pRankStack2, const AT_RANK *nAtomNumberCanonFrom, int num_atoms)` |
| `CompareLinCtStereoDble` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:178` | external | tree-sitter+gcc-aux | `extern int CompareLinCtStereoDble (AT_STEREO_DBLE *LinearCTStereoDble1, int nLenLinearCTStereoDble1, AT_STEREO_DBLE *LinearCTStereoDble2, int nLenLinearCTStereoDble2)` |
| `CompareLinCtStereoCarb` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:223` | external | tree-sitter+gcc-aux | `extern int CompareLinCtStereoCarb (AT_STEREO_CARB *LinearCTStereoCarb1, int nLenLinearCTStereoCarb1, AT_STEREO_CARB *LinearCTStereoCarb2, int nLenLinearCTStereoCarb2)` |
| `CompareLinCtStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:262` | external | tree-sitter+gcc-aux | `extern int CompareLinCtStereo (AT_STEREO_DBLE *LinearCTStereoDble1, int nLenLinearCTStereoDble1, AT_STEREO_CARB *LinearCTStereoCarb1, int nLenLinearCTStereoCarb1, AT_STEREO_DBLE *LinearCTStereoDble2, int nLenLinearCTStereoDble2, AT_STEREO_CARB *LinearCTStereoCarb2, int nLenLinearCTStereoCarb2)` |
| `CompareLinCtStereoAtomToValues` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:287` | external | tree-sitter+gcc-aux | `extern int CompareLinCtStereoAtomToValues (AT_STEREO_CARB *LinearCTStereoCarb, AT_RANK at_rank_canon1, U_CHAR parity)` |
| `bUniqueAtNbrFromMappingRank` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:316` | external | tree-sitter+gcc-aux | `extern int bUniqueAtNbrFromMappingRank (AT_RANK **pRankStack, AT_RANK nAtRank, AT_NUMB *nAtNumber)` |
| `nGetMcr` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:336` | external | tree-sitter+gcc-aux | `extern AT_RANK nGetMcr (AT_RANK *nEqArray, AT_RANK n)` |
| `nJoin2Mcrs` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:366` | external | tree-sitter+gcc-aux | `extern int nJoin2Mcrs (AT_RANK *nEqArray, AT_RANK n1, AT_RANK n2)` |
| `All_SB_Same` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:393` | external | tree-sitter+gcc-aux | `extern int All_SB_Same (AT_RANK canon_rank1, AT_RANK canon_rank2, const const ppAT_RANK pRankStack1, const const ppAT_RANK pRankStack2, const AT_RANK *nAtomNumberCanonFrom, sp_ATOM *at)` |
| `Next_SB_At_CanonRanks2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:518` | external | tree-sitter+gcc-aux | `extern int Next_SB_At_CanonRanks2 (AT_RANK *canon_rank1, AT_RANK *canon_rank2, AT_RANK *canon_rank1_min, AT_RANK *canon_rank2_min, int *bFirstTime, S_CHAR *bAtomUsedForStereo, const const ppAT_RANK pRankStack1, const const ppAT_RANK pRankStack2, const AT_RANK *nCanonRankFrom, const AT_RANK *nAtomNumberCanonFrom, const sp_ATOM *at, int num_atoms, int bAllene)` |
| `NextStereoParity2Test` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:659` | external | tree-sitter+gcc-aux | `extern int NextStereoParity2Test (int *stereo_bond_parity, int *sb_parity_calc, int nNumBest, int nNumWorse, int nNumUnkn, int nNumUndf, int nNumCalc, int vABParityUnknown)` |
| `CompareLinCtStereoDoubleToValues` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:764` | external | tree-sitter+gcc-aux | `extern int CompareLinCtStereoDoubleToValues (AT_STEREO_DBLE *LinearCTStereoDble, AT_RANK at_rank_canon1, AT_RANK at_rank_canon2, U_CHAR bond_parity)` |
| `SetUseAtomForStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:804` | external | tree-sitter+gcc-aux | `extern void SetUseAtomForStereo (S_CHAR *bAtomUsedForStereo, sp_ATOM *at, int num_atoms)` |
| `CurTreeAlloc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:823` | external | tree-sitter+gcc-aux | `extern int CurTreeAlloc (CUR_TREE *cur_tree, int num_atoms)` |
| `CurTreeReAlloc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:850` | external | tree-sitter+gcc-aux | `extern int CurTreeReAlloc (CUR_TREE *cur_tree)` |
| `CurTreeFree` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:871` | external | tree-sitter+gcc-aux | `extern void CurTreeFree (CUR_TREE *cur_tree)` |
| `CurTreeAddRank` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:881` | external | tree-sitter+gcc-aux | `extern int CurTreeAddRank (CUR_TREE *cur_tree, AT_NUMB rank)` |
| `CurTreeIsLastRank` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:902` | external | tree-sitter+gcc-aux | `extern int CurTreeIsLastRank (CUR_TREE *cur_tree, AT_NUMB rank)` |
| `CurTreeRemoveLastRankIfNoAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:920` | external | tree-sitter+gcc-aux | `extern int CurTreeRemoveLastRankIfNoAtoms (CUR_TREE *cur_tree)` |
| `CurTreeAddAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:935` | external | tree-sitter+gcc-aux | `extern int CurTreeAddAtom (CUR_TREE *cur_tree, int at_no)` |
| `CurTreeKeepLastAtomsOnly` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:960` | external | tree-sitter+gcc-aux | `extern void CurTreeKeepLastAtomsOnly (CUR_TREE *cur_tree, int tpos, int shift)` |
| `CurTreeRemoveIfLastAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:993` | external | tree-sitter+gcc-aux | `extern int CurTreeRemoveIfLastAtom (CUR_TREE *cur_tree, int at_no)` |
| `CurTreeGetPos` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:1011` | external | tree-sitter+gcc-aux | `extern int CurTreeGetPos (CUR_TREE *cur_tree)` |
| `CurTreeSetPos` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:1023` | external | tree-sitter+gcc-aux | `extern int CurTreeSetPos (CUR_TREE *cur_tree, int len)` |
| `CurTreeRemoveLastRank` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:1036` | external | tree-sitter+gcc-aux | `extern int CurTreeRemoveLastRank (CUR_TREE *cur_tree)` |
| `CurTreeIsLastAtomEqu` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:1054` | external | tree-sitter+gcc-aux | `extern int CurTreeIsLastAtomEqu (CUR_TREE *cur_tree, int at_no, AT_NUMB *nSymmStereo)` |
| `CurTreeRemoveLastAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:1075` | external | tree-sitter | `int CurTreeRemoveLastAtom( CUR_TREE *cur_tree )` |
| `CurTreeReplaceLastRank` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:1091` | external | tree-sitter | `int CurTreeReplaceLastRank( CUR_TREE *cur_tree, AT_NUMB rank )` |
| `CurTreeFindTheRankPos` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap1.c:1104` | external | tree-sitter | `int CurTreeFindTheRankPos( CUR_TREE *cur_tree, AT_NUMB rank )` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c`

Parse errors: `1`. Function definitions: `30`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `SortedEquInfoToRanks` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:148` | external | tree-sitter+gcc-aux | `extern int SortedEquInfoToRanks (const AT_RANK *nSymmRank, AT_RANK *nRank, const AT_RANK *nAtomNumber, int num_atoms, int *bChanged)` |
| `SortedRanksToEquInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:199` | external | tree-sitter+gcc-aux | `extern int SortedRanksToEquInfo (AT_RANK *nSymmRank, const AT_RANK *nRank, const AT_RANK *nAtomNumber, int num_atoms)` |
| `switch_ptrs` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:230` | external | tree-sitter+gcc-aux | `extern void switch_ptrs (AT_RANK **p1, AT_RANK **p2)` |
| `SetNewRanksFromNeighLists3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:241` | external | tree-sitter+gcc-aux | `extern int SetNewRanksFromNeighLists3 (CANON_GLOBALS *pCG, int num_atoms, NEIGH_LIST *NeighList, AT_RANK *nRank, AT_RANK *nNewRank, AT_RANK *nAtomNumber)` |
| `SetNewRanksFromNeighLists4` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:308` | external | tree-sitter+gcc-aux | `extern int SetNewRanksFromNeighLists4 (CANON_GLOBALS *pCG, int num_atoms, NEIGH_LIST *NeighList, AT_RANK *nRank, AT_RANK *nNewRank, AT_RANK *nAtomNumber, AT_RANK nMaxAtRank)` |
| `SetNewRanksFromNeighLists` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:380` | external | tree-sitter+gcc-aux | `extern int SetNewRanksFromNeighLists (CANON_GLOBALS *pCG, int num_atoms, NEIGH_LIST *NeighList, AT_RANK *nRank, AT_RANK *nNewRank, AT_RANK *nAtomNumber, int bUseAltSort, int (*comp) (const void *, const void *, void *))` |
| `SortNeighListsBySymmAndCanonRank` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:434` | external | tree-sitter+gcc-aux | `extern void SortNeighListsBySymmAndCanonRank (int num_atoms, NEIGH_LIST *NeighList, const AT_RANK *nSymmRank, const AT_RANK *nCanonRank)` |
| `SortNeighLists2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:448` | external | tree-sitter+gcc-aux | `extern int SortNeighLists2 (int num_atoms, AT_RANK *nRank, NEIGH_LIST *NeighList, AT_RANK *nAtomNumber)` |
| `SortNeighLists3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:482` | external | tree-sitter+gcc-aux | `extern int SortNeighLists3 (int num_atoms, AT_RANK *nRank, NEIGH_LIST *NeighList, AT_RANK *nAtomNumber)` |
| `DifferentiateRanks2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:518` | external | tree-sitter+gcc-aux | `extern int DifferentiateRanks2 (CANON_GLOBALS *pCG, int num_atoms, NEIGH_LIST *NeighList, int nNumCurrRanks, AT_RANK *pnCurrRank, AT_RANK *pnPrevRank, AT_RANK *nAtomNumber, long int *lNumIter, int bUseAltSort)` |
| `DifferentiateRanks3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:561` | external | tree-sitter+gcc-aux | `extern int DifferentiateRanks3 (CANON_GLOBALS *pCG, int num_atoms, NEIGH_LIST *NeighList, int nNumCurrRanks, AT_RANK *pnCurrRank, AT_RANK *pnPrevRank, AT_RANK *nAtomNumber, long int *lNumIter)` |
| `DifferentiateRanks4` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:602` | external | tree-sitter+gcc-aux | `extern int DifferentiateRanks4 (CANON_GLOBALS *pCG, int num_atoms, NEIGH_LIST *NeighList, int nNumCurrRanks, AT_RANK *pnCurrRank, AT_RANK *pnPrevRank, AT_RANK *nAtomNumber, AT_RANK nMaxAtRank, long int *lNumIter)` |
| `DifferentiateRanksBasic` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:637` | external | tree-sitter+gcc-aux | `extern int DifferentiateRanksBasic (CANON_GLOBALS *pCG, int num_atoms, NEIGH_LIST *NeighList, int nNumCurrRanks, AT_RANK *pnCurrRank, AT_RANK *pnPrevRank, AT_RANK *nAtomNumber, long int *lNumIter, int bUseAltSort)` |
| `NumberOfTies` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:680` | external | tree-sitter+gcc-aux | `extern int NumberOfTies (AT_RANK **pRankStack1, AT_RANK **pRankStack2, int length, int at_no1, int at_no2, AT_RANK *nNewRank, int *bAddStack, int *bMapped1)` |
| `HalfStereoBondParity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:802` | external | tree-sitter+gcc-aux | `extern int HalfStereoBondParity (sp_ATOM *at, int at_no1, int i_sb_neigh, const AT_RANK *nRank)` |
| `parity_of_mapped_half_bond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:958` | external | tree-sitter+gcc-aux | `extern int parity_of_mapped_half_bond (int from_at, int to_at, int from_neigh, int to_neigh, sp_ATOM *at, EQ_NEIGH *pEN, const AT_RANK *nCanonRankFrom, const AT_RANK *nRankFrom, const AT_RANK *nRankTo)` |
| `parity_of_mapped_atom2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:1168` | external | tree-sitter+gcc-aux | `extern int parity_of_mapped_atom2 (CANON_GLOBALS *pCG, int from_at, int to_at, const sp_ATOM *at, EQ_NEIGH *pEN, const AT_RANK *nCanonRankFrom, const AT_RANK *nRankFrom, const AT_RANK *nRankTo)` |
| `ClearPreviousMappings` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:1334` | external | tree-sitter+gcc-aux | `extern int ClearPreviousMappings (AT_RANK **pRankStack1)` |
| `map_an_atom2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:1350` | external | tree-sitter+gcc-aux | `extern int map_an_atom2 (CANON_GLOBALS *pCG, int num_atoms, int num_max, int at_no1, int at_no2, AT_RANK *nTempRank, int nNumMappedRanks, int *pnNewNumMappedRanks, CANON_STAT *pCS, NEIGH_LIST *NeighList, AT_RANK **pRankStack1, AT_RANK **pRankStack2, int *bAddStack)` |
| `might_change_other_atom_parity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:1461` | external | tree-sitter+gcc-aux | `extern int might_change_other_atom_parity (sp_ATOM *at, int num_atoms, int at_no, AT_RANK *nRank2, AT_RANK *nRank1)` |
| `DeAllocateForNonStereoRemoval` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:1504` | external | tree-sitter+gcc-aux | `extern void DeAllocateForNonStereoRemoval (AT_RANK **nAtomNumberCanon1, AT_RANK **nAtomNumberCanon2, NEIGH_LIST **nl, NEIGH_LIST **nl1, NEIGH_LIST **nl2, AT_RANK **nVisited1, AT_RANK **nVisited2)` |
| `AllocateForNonStereoRemoval` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:1551` | external | tree-sitter+gcc-aux | `extern int AllocateForNonStereoRemoval (sp_ATOM *at, int num_atoms, const AT_RANK *nSymmRank, AT_RANK *nCanonRank, AT_RANK **nAtomNumberCanon1, AT_RANK **nAtomNumberCanon2, NEIGH_LIST **nl, NEIGH_LIST **nl1, NEIGH_LIST **nl2, AT_RANK **nVisited1, AT_RANK **nVisited2)` |
| `GetMinNewRank` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:1587` | external | tree-sitter+gcc-aux | `extern AT_RANK GetMinNewRank (AT_RANK *nAtomRank, AT_RANK *nAtomNumb, AT_RANK nRank1)` |
| `BreakNeighborsTie` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:1609` | external | tree-sitter+gcc-aux | `extern int BreakNeighborsTie (CANON_GLOBALS *pCG, sp_ATOM *at, int num_atoms, int num_at_tg, int ib, int ia, AT_RANK *neigh_num, int in1, int in2, int mode, AT_RANK **pRankStack1, AT_RANK **pRankStack2, AT_RANK *nTempRank, NEIGH_LIST *NeighList, const AT_RANK *nSymmRank, AT_RANK *nCanonRank, NEIGH_LIST *nl1, NEIGH_LIST *nl2, long int *lNumIter)` |
| `CheckNextSymmNeighborsAndBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:2060` | external | tree-sitter+gcc-aux | `extern int CheckNextSymmNeighborsAndBonds (sp_ATOM *at, AT_RANK cur1, AT_RANK cur2, AT_RANK n1, AT_RANK n2, AT_RANK *nAvoidCheckAtom, AT_RANK *nVisited1, AT_RANK *nVisited2, AT_RANK *nVisitOrd1, AT_RANK *nVisitOrd2, const AT_RANK *nRank1, const AT_RANK *nRank2)` |
| `CreateCheckSymmPaths` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:2191` | external | tree-sitter+gcc-aux | `extern int CreateCheckSymmPaths (sp_ATOM *at, AT_RANK prev1, AT_RANK cur1, AT_RANK prev2, AT_RANK cur2, AT_RANK *nAvoidCheckAtom, AT_RANK *nVisited1, AT_RANK *nVisited2, AT_RANK *nVisitOrd1, AT_RANK *nVisitOrd2, NEIGH_LIST *nl1, NEIGH_LIST *nl2, const AT_RANK *nRank1, const AT_RANK *nRank2, AT_RANK *nCanonRank, AT_RANK *nLength, int *bParitiesInverted, int mode)` |
| `CalculatedPathsParitiesAreIdentical` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:2399` | external | tree-sitter+gcc-aux | `extern int CalculatedPathsParitiesAreIdentical (CANON_GLOBALS *pCG, sp_ATOM *at, int num_atoms, const AT_RANK *nSymmRank, AT_RANK *nCanonRank, AT_RANK *nAtomNumberCanon, AT_RANK *nAtomNumberCanon1, AT_RANK *nAtomNumberCanon2, AT_RANK *nVisited1, AT_RANK *nVisited2, AT_RANK prev_sb_neigh, AT_RANK cur, AT_RANK next1, AT_RANK next2, int nNeighMode, int bParitiesInverted, int mode, CANON_STAT *pCS, int vABParityUnknown)` |
| `RemoveCalculatedNonStereoBondParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:3099` | external | tree-sitter+gcc-aux | `extern int RemoveCalculatedNonStereoBondParities (CANON_GLOBALS *pCG, sp_ATOM *at, int num_atoms, int num_at_tg, AT_RANK **pRankStack1, AT_RANK **pRankStack2, AT_RANK *nTempRank, NEIGH_LIST *NeighList, AT_RANK *nCanonRank, const AT_RANK *nSymmRank, AT_RANK *nAtomNumberCanon, AT_RANK *nAtomNumberCanon1, AT_RANK *nAtomNumberCanon2, NEIGH_LIST *nl, NEIGH_LIST *nl1, NEIGH_LIST *nl2, AT_RANK *nVisited1, AT_RANK *nVisited2, CANON_STAT *pCS, int vABParityUnknown)` |
| `RemoveCalculatedNonStereoCenterParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:3356` | external | tree-sitter+gcc-aux | `extern int RemoveCalculatedNonStereoCenterParities (CANON_GLOBALS *pCG, sp_ATOM *at, int num_atoms, int num_at_tg, AT_RANK **pRankStack1, AT_RANK **pRankStack2, AT_RANK *nTempRank, NEIGH_LIST *NeighList, AT_RANK *nCanonRank, const AT_RANK *nSymmRank, AT_RANK *nAtomNumberCanon, AT_RANK *nAtomNumberCanon1, AT_RANK *nAtomNumberCanon2, NEIGH_LIST *nl, NEIGH_LIST *nl1, NEIGH_LIST *nl2, AT_RANK *nVisited1, AT_RANK *nVisited2, CANON_STAT *pCS, int vABParityUnknown)` |
| `RemoveCalculatedNonStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap2.c:3672` | external | tree-sitter+gcc-aux | `extern int RemoveCalculatedNonStereo (CANON_GLOBALS *pCG, sp_ATOM *at, int num_atoms, int num_at_tg, AT_RANK **pRankStack1, AT_RANK **pRankStack2, AT_RANK *nTempRank, NEIGH_LIST *NeighList, const AT_RANK *nSymmRank, AT_RANK *nCanonRank, AT_RANK *nAtomNumberCanon, CANON_STAT *pCS, int vABParityUnknown)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap4.c`

Parse errors: `7`. Function definitions: `2`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `map_stereo_bonds4` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap4.c:83` | external | tree-sitter+gcc-aux | `extern int map_stereo_bonds4 (struct tagINCHI_CLOCK *ic, CANON_GLOBALS *pCG, sp_ATOM *at, int num_atoms, int num_at_tg, int num_max, int bAllene, const AT_RANK *nCanonRankFrom, const AT_RANK *nAtomNumberCanonFrom, AT_RANK *nCanonRankTo, const AT_RANK *nSymmRank, AT_RANK **pRankStack1, AT_RANK **pRankStack2, AT_RANK *nTempRank, int nNumMappedRanksInput, AT_RANK *nSymmStereo, NEIGH_LIST *NeighList, CANON_STAT *pCS, CUR_TREE *cur_tree, int nNumMappedBonds, int vABParityUnknown)` |
| `map_stereo_atoms4` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichimap4.c:1126` | external | gcc-aux | `extern int map_stereo_atoms4 (struct tagINCHI_CLOCK *ic, CANON_GLOBALS *pCG, sp_ATOM *at, int num_atoms, int num_at_tg, int num_max, const AT_RANK *nCanonRankFrom, const AT_RANK *nAtomNumberCanonFrom, AT_RANK *nCanonRankTo, const AT_RANK *nSymmRank, AT_RANK **pRankStack1, AT_RANK **pRankStack2, AT_RANK *nTempRank, int nNumMappedRanksInput, AT_RANK *nSymmStereo, NEIGH_LIST *NeighList, CANON_STAT *pCS, CUR_TREE *cur_tree, int nNumMappedAtoms, int vABParityUnknown)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c`

Parse errors: `80`. Function definitions: `77`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `MarkRingSystemsInp` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:59` | external | tree-sitter+gcc-aux | `extern int MarkRingSystemsInp (inp_ATOM *at, int num_atoms, int start)` |
| `UnMarkDisconnectedComponents` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:652` | external | tree-sitter+gcc-aux | `extern int UnMarkDisconnectedComponents (ORIG_ATOM_DATA *orig_inp_data)` |
| `UnMarkOtherIndicators` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:681` | external | tree-sitter+gcc-aux | `extern int UnMarkOtherIndicators (inp_ATOM *at, int num_atoms)` |
| `UnMarkOneComponent` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:698` | external | tree-sitter+gcc-aux | `extern int UnMarkOneComponent (inp_ATOM *at, int num_atoms)` |
| `set_R2C_el_numbers` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:712` | external | tree-sitter+gcc-aux | `extern void set_R2C_el_numbers (void)` |
| `subtract_DT_from_num_H` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:734` | external | tree-sitter+gcc-aux | `extern int subtract_DT_from_num_H (int num_atoms, inp_ATOM *at)` |
| `add_inp_ATOM` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:749` | external | tree-sitter+gcc-aux | `extern int add_inp_ATOM (inp_ATOM *at, int len_at, int len_cur, inp_ATOM *add, int len_add)` |
| `mark_arom_bonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:784` | external | tree-sitter+gcc-aux | `extern int mark_arom_bonds (struct tagINCHI_CLOCK *ic, struct tagCANON_GLOBALS *pCG, inp_ATOM *at, int num_atoms)` |
| `cmp_r2c_atpair` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:800` | external | tree-sitter+gcc-aux | `extern int cmp_r2c_atpair (const void *p1, const void *p2)` |
| `has_atom_pair_seq` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:815` | external | tree-sitter+gcc-aux | `extern int has_atom_pair_seq (R2C_ATPAIR *ap, int num_ap, AT_NUMB at1, AT_NUMB at2)` |
| `has_atom_pair` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:835` | external | tree-sitter+gcc-aux | `extern int has_atom_pair (R2C_ATPAIR *ap, int num_ap, AT_NUMB at1, AT_NUMB at2)` |
| `mark_atoms_ap` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:871` | external | tree-sitter+gcc-aux | `extern int mark_atoms_ap (inp_ATOM *at, AT_NUMB start, R2C_ATPAIR *ap, int num_ap, int num, AT_NUMB cFlags)` |
| `mark_atoms_cFlags` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1368` | external | tree-sitter+gcc-aux | `extern int mark_atoms_cFlags (inp_ATOM *at, int start, int num, char cFlags)` |
| `unmark_atoms_cFlags` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1388` | external | tree-sitter+gcc-aux | `extern int unmark_atoms_cFlags (inp_ATOM *at, int start, int num, char cFlags, char cInvFlags)` |
| `is_C_or_S_DB_O` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1410` | external | tree-sitter+gcc-aux | `extern int is_C_or_S_DB_O (inp_ATOM *at, int i)` |
| `is_C_DB_O` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1434` | external | tree-sitter+gcc-aux | `extern int is_C_DB_O (inp_ATOM *at, int i)` |
| `is_C_unsat_not_arom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1457` | external | tree-sitter+gcc-aux | `extern int is_C_unsat_not_arom (inp_ATOM *at, int i)` |
| `is_Aryl` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1486` | external | tree-sitter+gcc-aux | `extern int is_Aryl (inp_ATOM *at, int outside_point, int attachment_pont)` |
| `is_Saturated_C` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1509` | external | tree-sitter+gcc-aux | `extern int is_Saturated_C (inp_ATOM *at, int attachment_pont)` |
| `is_C_Alk` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1517` | external | tree-sitter+gcc-aux | `extern int is_C_Alk (inp_ATOM *at, int i, char cFlags)` |
| `is_Phenyl` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1559` | external | tree-sitter+gcc-aux | `extern int is_Phenyl (inp_ATOM *at, int outside_point, int attachment_point)` |
| `is_PentaFluoroPhenyl` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1615` | external | tree-sitter+gcc-aux | `extern int is_PentaFluoroPhenyl (inp_ATOM *at, int outside_point, int attachment_point)` |
| `is_Methyl` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1679` | external | tree-sitter+gcc-aux | `extern int is_Methyl (inp_ATOM *at, int attachment_point)` |
| `is_Ethyl` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1694` | external | tree-sitter+gcc-aux | `extern int is_Ethyl (inp_ATOM *at, int outside_point, int attachment_point)` |
| `is_Methyl_or_Etyl` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1709` | external | tree-sitter+gcc-aux | `extern int is_Methyl_or_Etyl (inp_ATOM *at, int outside_point, int attachment_point)` |
| `is_Si_IV` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1727` | external | tree-sitter+gcc-aux | `extern int is_Si_IV (inp_ATOM *at, int i)` |
| `is_P_TB_N` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1740` | external | tree-sitter+gcc-aux | `extern int is_P_TB_N (inp_ATOM *at, int i)` |
| `get_CO_opposite` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1763` | external | tree-sitter+gcc-aux | `extern int get_CO_opposite (inp_ATOM *at, int iat, int iord, int *iat_opposite, int *iord_opposite)` |
| `is_DERIV_RING_DMOX_DEOX_O` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1802` | external | tree-sitter+gcc-aux | `extern int is_DERIV_RING_DMOX_DEOX_O (inp_ATOM *at, int cur_atom, int from_ord, DERIV_AT *da, DERIV_AT *da1)` |
| `is_DERIV_RING_DMOX_DEOX_N` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:1962` | external | tree-sitter+gcc-aux | `extern int is_DERIV_RING_DMOX_DEOX_N (inp_ATOM *at, int cur_atom, int from_ord, DERIV_AT *da, DERIV_AT *da1)` |
| `is_DERIV_RING2_PRRLDD_PPRDN` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:2147` | external | tree-sitter+gcc-aux | `extern int is_DERIV_RING2_PRRLDD_PPRDN (inp_ATOM *at, int cur_atom, int from_ord, DERIV_AT *da, DERIV_AT *da1)` |
| `check_arom_chain` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:2281` | external | tree-sitter+gcc-aux | `extern int check_arom_chain (inp_ATOM *at, int cur, int from, int last, int len)` |
| `is_Dansyl` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:2334` | external | tree-sitter+gcc-aux | `extern int is_Dansyl (inp_ATOM *at, int cur_atom, int to_ord, DERIV_AT *da, DERIV_AT *da1)` |
| `is_possibly_deriv_neigh` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:2610` | external | tree-sitter+gcc-aux | `extern int is_possibly_deriv_neigh (inp_ATOM *at, int iat, int iord, int type, char cFlags)` |
| `get_traversed_deriv_type` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:2677` | external | tree-sitter+gcc-aux | `extern int get_traversed_deriv_type (inp_ATOM *at, DERIV_AT *da, int k, DERIV_AT *da1, char cFlags)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:2937` | external | tree-sitter | `else if (is_Methyl( at, n2 ))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:2948` | external | tree-sitter | `else if (is_Ethyl( at, n1, n2 ))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:2958` | external | tree-sitter | `else if (is_Phenyl( at, n1, n2 ))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:2968` | external | tree-sitter | `else if (is_PentaFluoroPhenyl( at, n1, n2 ))` |
| `add_to_da` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:3245` | external | tree-sitter+gcc-aux | `extern int add_to_da (DERIV_AT *da, DERIV_AT *add)` |
| `mark_atoms_deriv` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:3332` | external | tree-sitter+gcc-aux | `extern int mark_atoms_deriv (inp_ATOM *at, DERIV_AT *da, int start, int num, char cFlags, int *pbFound)` |
| `count_one_bond_atoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:3536` | external | tree-sitter+gcc-aux | `extern int count_one_bond_atoms (inp_ATOM *at, DERIV_AT *da, int start, int ord, char cFlags, int *bFound)` |
| `is_silyl` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:3633` | external | tree-sitter+gcc-aux | `extern int is_silyl (inp_ATOM *at, int start, int ord_prev)` |
| `is_silyl2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:3738` | external | tree-sitter+gcc-aux | `extern int is_silyl2 (inp_ATOM *at, int start, int from_at)` |
| `is_nButyl` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:3831` | external | tree-sitter+gcc-aux | `extern int is_nButyl (inp_ATOM *at, int start, int ord_prev)` |
| `is_Me_or_Et` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:3862` | external | tree-sitter+gcc-aux | `extern int is_Me_or_Et (inp_ATOM *at, int start, int ord_prev)` |
| `is_CF3_or_isoC3F7` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:3904` | external | tree-sitter | `int is_CF3_or_isoC3F7( inp_ATOM *at, int start, int ord_prev )` |
| `is_CF3_or_linC3F7a` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:4000` | external | tree-sitter+gcc-aux | `extern int is_CF3_or_linC3F7a (inp_ATOM *at, int start, int iat_prev)` |
| `is_CF3_or_linC3F7` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:4017` | external | tree-sitter+gcc-aux | `extern int is_CF3_or_linC3F7 (inp_ATOM *at, int start, int ord_prev)` |
| `is_phenyl` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:4103` | external | tree-sitter+gcc-aux | `extern int is_phenyl (inp_ATOM *at, int start, int ord_prev)` |
| `is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:4136` | external | tree-sitter+gcc-aux | `extern int is_DERIV_RING_O_or_NH_OUTSIDE_PRECURSOR (inp_ATOM *at, int start, int num_atoms, DERIV_AT *da1, int idrv, char *szUnderiv, int lenUnderiv, char *szUnderiv2, int lenUnderiv2, BIT_UNDERIV *bitUnderiv)` |
| `is_deriv_chain2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:4483` | external | tree-sitter+gcc-aux | `extern int is_deriv_chain2 (inp_ATOM *at, int start, int type, int num, int ord, int idrv, char *szUnderiv, int lenUnderiv, char *szUnderiv2, int lenUnderiv2, BIT_UNDERIV *bitUnderiv)` |
| `is_deriv_chain` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:4980` | external | tree-sitter+gcc-aux | `extern int is_deriv_chain (inp_ATOM *at, int start, int num_atoms, DERIV_AT *da1, int idrv, char *szUnderiv, int lenUnderiv, char *szUnderiv2, int lenUnderiv2, BIT_UNDERIV *bitUnderiv)` |
| `is_deriv_chain_or_ring` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:4996` | external | tree-sitter+gcc-aux | `extern int is_deriv_chain_or_ring (inp_ATOM *at, int start, int num_atoms, DERIV_AT *da1, int *idrv)` |
| `remove_deriv` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5046` | external | tree-sitter+gcc-aux | `extern int remove_deriv (DERIV_AT *da1, int idrv)` |
| `remove_deriv_mark` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5108` | external | tree-sitter+gcc-aux | `extern int remove_deriv_mark (DERIV_AT *da1, int idrv)` |
| `underiv_buf_clear` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5147` | external | tree-sitter+gcc-aux | `extern void underiv_buf_clear (char *szUnderiv)` |
| `underiv_list_add` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5156` | external | tree-sitter+gcc-aux | `extern int underiv_list_add (char *szUnderivList, int lenUnderivList, const char *szUnderiv, char cDelimiter)` |
| `underiv_list_get_last` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5179` | external | tree-sitter+gcc-aux | `extern const char *underiv_list_get_last (const char *szUnderivList, char cDelimiter)` |
| `underiv_compare` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5197` | external | tree-sitter+gcc-aux | `extern int underiv_compare (const void *p1, const void *p2)` |
| `underiv_list_add_two_cuts` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5204` | external | tree-sitter+gcc-aux | `extern int underiv_list_add_two_cuts (char *szUnderivList, int lenUnderivList, char *szUnderiv, const const char cDelim)` |
| `sort_merge_underiv` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5281` | external | tree-sitter+gcc-aux | `extern int sort_merge_underiv (char *pSdfValue, int bOutputSdf, char *szUnderivList, char cDerivSeparator, const char *pszUnderivPrefix, const char *pszUnderivPostfix)` |
| `eliminate_deriv_not_in_list` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5351` | external | tree-sitter+gcc-aux | `extern int eliminate_deriv_not_in_list (inp_ATOM *at, DERIV_AT *da, int num_atoms, char *szUnderivList, int lenUnderivList, char *szUnderivList2, int lenUnderivList2, BIT_UNDERIV *bitUnderivList)` |
| `make_single_cut` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5495` | external | tree-sitter+gcc-aux | `extern int make_single_cut (inp_ATOM *at, DERIV_AT *da, int iat, int icut)` |
| `fill_out_bond_cuts` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5628` | external | tree-sitter+gcc-aux | `extern int fill_out_bond_cuts (inp_ATOM *at, DERIV_AT *da, int num_atoms, R2C_ATPAIR *ap, int num_cuts_to_check)` |
| `mark_deriv_agents` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5743` | external | tree-sitter+gcc-aux | `extern int mark_deriv_agents (inp_ATOM *at, DERIV_AT *da, int num_atoms, R2C_ATPAIR *ap, int num_cuts_to_check, AT_NUMB *pnum_comp, int *pcur_num_at)` |
| `replace_arom_bonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5809` | external | tree-sitter+gcc-aux | `extern int replace_arom_bonds (inp_ATOM *at, int num_atoms, inp_ATOM *at2, int num_atoms2)` |
| `add_explicit_H` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5868` | external | tree-sitter+gcc-aux | `extern int add_explicit_H (INP_ATOM_DATA *inp_cur_data)` |
| `OAD_Edit_Underivatize` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:5934` | external | tree-sitter+gcc-aux | `extern int OAD_Edit_Underivatize (struct tagINCHI_CLOCK *ic, struct tagCANON_GLOBALS *pCG, ORIG_ATOM_DATA *orig_inp_data, int bOutputSdf, int bOutputReport, char *pSdfValue)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:6213` | external | tree-sitter | `else if (da[i].typ[0] == da[i].typ[1])` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:6605` | external | tree-sitter | `ALLOC_AP if (num_cuts_to_check >= 2)` |
| `detect_r2c_Zatom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:7159` | external | tree-sitter+gcc-aux | `extern int detect_r2c_Zatom (inp_ATOM *at, R2C_AT *da, int iZ)` |
| `cut_ring_to_chain` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:7296` | external | tree-sitter+gcc-aux | `extern int cut_ring_to_chain (inp_ATOM *at, R2C_AT *da, int iZ)` |
| `Ring2Chain` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:7362` | external | tree-sitter+gcc-aux | `extern int Ring2Chain (struct tagINCHI_CLOCK *ic, struct tagCANON_GLOBALS *pCG, ORIG_ATOM_DATA *orig_inp_data)` |
| `OAD_Edit_MergeComponentsAndRecreateOAD` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:7563` | external | tree-sitter+gcc-aux | `extern void OAD_Edit_MergeComponentsAndRecreateOAD (ORIG_ATOM_DATA *orig_OrigAtomData, INP_ATOM_DATA *curr_InpAtomData, int num_components, int *errcode)` |
| `free_underiv_temp_data` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:7632` | external | tree-sitter+gcc-aux | `extern void free_underiv_temp_data (R2C_ATPAIR *ap, DERIV_AT *da, inp_ATOM *at2, INP_ATOM_DATA *inp_cur_data, int num_components)` |
| `remove_cut_derivs` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichinorm.c:7666` | external | tree-sitter+gcc-aux | `extern void remove_cut_derivs (int num_atoms, inp_ATOM *at, INP_ATOM_DATA *inp_cur_data, int i_component, int *errcode)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiparm.c`

Parse errors: `156`. Function definitions: `10`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `set_common_options_by_parg` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiparm.c:131` | external | tree-sitter+gcc-aux | `extern int set_common_options_by_parg (const char *pArg, int developer_options, INPUT_PARMS *ip, INCHI_MODE *pbVer1DefaultMode, int *pnMode, int *pbINChIOutputOptions, int *pbINChIOutputOptions2, int *pbStdFormat, int *pbHashKey, int *pbHashXtra1, int *pbHashXtra2, int *pbFixSp3bug, int *pbFixFB2, int *pbAddPhosphineStereo, int *pbAddArsineStereo, int *pbNoStructLabels, int *pbPointedEdgeStereo, int *pbDoNotAddH, int *pbForcedChiralFlag, int *pbReconnectCoord, int *pbKetoEnolTaut, int *pb15TautNonRing, int *pbPT_06_00_Taut, int *pbPT_13_00_Taut, int *pbPT_16_00_Taut, int *pbPT_18_00_Taut, int *pbPT_22_00_Taut, int *pbPT_39_00_Taut, int *pbLooseTSACheck, int *pbLargeMolecules, int *pbPolymers, int *pbFoldPolymerSRU, int *pbFrameShiftScheme, int *pbStereoAtZz, int *pbNPZz, int *pbNoWarnings, int *pbMergeHash, int *pbHideInChI)` |
| `ReadCommandLineParms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiparm.c:602` | external | tree-sitter+gcc-aux | `extern int ReadCommandLineParms (int argc, const char **argv, INPUT_PARMS *ip, char *szSdfDataValue, long unsigned int *ulDisplTime, int bReleaseVersion, INCHI_IOSTREAM *log_file)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiparm.c:1221` | external | tree-sitter | `else if ((bVer1Options & 1) && INCHI_OPTION_PREFX == argv[i][0] && argv[i][1]) { /* Parsing stand-alone executable/libinchi options */ pArg = argv[i] + 1; bRecognizedOption = 2; bVer1Options += 2; /* always on: REQ_MODE_TAUT \| REQ_MODE_ISO \| REQ_MODE_STEREO */ got = set_common_options_by_parg(pArg, developer_options, ip, &bVer1DefaultMode, &nMode, &bINChIOutputOptions, &bINChIOutputOptions2, &bStdFormat, &bHashKey, &bHashXtra1, &bHashXtra2, &bFixSp3bug, &bFixFB2, &bAddPhosphineStereo, &bAddArsineStereo, &bNoStructLabels, &bPointedEdgeStereo, &bDoNotAddH, &bForcedChiralFlag, &bReconnectCoord, &bKetoEnolTaut, &b15TautNonRing, &bPT_06_00_Taut, &bPT_13_00_Taut, &bPT_16_00_Taut, &bPT_18_00_Taut, &bPT_22_00_Taut, &bPT_39_00_Taut, &bLooseTSACheck, &bLargeMolecules, &bPolymers, &bFoldPolymerSRU, &bFrameShiftScheme, &bStereoAtZz, &bNPZz, &bNoWarnings, &bMergeHash, &bHideInChI); if (got)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiparm.c:1665` | external | tree-sitter | `else if (ip->num_paths < MAX_NUM_PATHS)` |
| `PrintInputParms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiparm.c:2130` | external | tree-sitter+gcc-aux | `extern int PrintInputParms (INCHI_IOSTREAM *log_file, INPUT_PARMS *ip)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiparm.c:2214` | external | tree-sitter | `else if (bInChI2Struct)` |
| `HelpCommandLineParms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiparm.c:2700` | external | tree-sitter+gcc-aux | `extern void HelpCommandLineParms (INCHI_IOSTREAM *f)` |
| `OpenFiles` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiparm.c:2880` | external | tree-sitter | `int OpenFiles(FILE** inp_file, FILE** out_file, FILE** log_file, FILE** prb_file, INPUT_PARMS* ip)` |
| `bMatchOnePrefix` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiparm.c:3043` | internal | tree-sitter | `static int bMatchOnePrefix(int len, char* str, int lenPrefix[], char strPrefix[][LEN_VERSIONS], int numPrefix)` |
| `DetectInputINChIFileType` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiparm.c:3062` | external | tree-sitter | `int DetectInputINChIFileType(FILE** inp_file, INPUT_PARMS* ip, const char* fmode)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c`

Parse errors: `18`. Function definitions: `39`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `OutputINChIPlainError` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:514` | external | tree-sitter+gcc-aux | `extern int OutputINChIPlainError (INCHI_IOSTREAM *out_file, char *pErrorText, int bError)` |
| `EquString` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:564` | external | tree-sitter+gcc-aux | `extern const char *EquString (int EquVal)` |
| `OutputINChI2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:966` | external | tree-sitter+gcc-aux | `extern int OutputINChI2 (CANON_GLOBALS *pCG, INCHI_IOS_STRING *strbuf, INCHI_SORT *(*pINChISortTautAndNonTaut2)[2], int INCHI_basic_or_INCHI_reconnected, ORIG_ATOM_DATA *orig_inp_data, ORIG_STRUCT *pOrigStruct, INPUT_PARMS *ip, int bDisconnectedCoord, int bOutputType, int bINChIOutputOptions, int *num_components2, int *num_non_taut2, int *num_taut2, INCHI_IOSTREAM *output_file, INCHI_IOSTREAM *log_file, int num_input_struct, int *pSortPrintINChIFlags, unsigned char save_opt_bits)` |
| `OutputINChI1` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:1043` | external | gcc-aux | `extern int OutputINChI1 (CANON_GLOBALS *pCG, INCHI_IOS_STRING *strbuf, INCHI_SORT *(*pINChISortTautAndNonTaut2)[2], int INCHI_basic_or_INCHI_reconnected, ORIG_ATOM_DATA *orig_inp_data, ORIG_STRUCT *pOrigStruct, INPUT_PARMS *ip, int bDisconnectedCoord, int bOutputType, int bINChIOutputOptions, int *num_components2, int *num_non_taut2, int *num_taut2, INCHI_IOSTREAM *out_file, INCHI_IOSTREAM *log_file, int num_input_struct, int *pSortPrintINChIFlags, unsigned char save_opt_bits)` |
| `szGetTag` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:2087` | external | tree-sitter+gcc-aux | `extern char *szGetTag (const INCHI_TAG *Tag, int nTag, int bTag, char *szTag, int *bAlways, short int tag_flag)` |
| `str_LineEnd` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:2183` | internal | tree-sitter+gcc-aux | `static int str_LineEnd (const char *tag, int *bOverflow, INCHI_IOS_STRING *buf, int ind, int bPlainTextTags)` |
| `CleanOrigCoord` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:2228` | internal | tree-sitter+gcc-aux | `static int CleanOrigCoord (char *szCoord, int delim)` |
| `WriteOrigCoord` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:2340` | internal | tree-sitter+gcc-aux | `static int WriteOrigCoord (int num_inp_atoms, MOL_COORD (*szMolCoord), int *i, char *szBuf, int buf_len)` |
| `WriteOrigAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:2410` | internal | tree-sitter+gcc-aux | `static int WriteOrigAtoms (CANON_GLOBALS *pCG, int num_inp_atoms, inp_ATOM *at, int *i, char *szBuf, int buf_len, STRUCT_DATA *sd)` |
| `WriteOrigBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:2610` | internal | tree-sitter+gcc-aux | `static int WriteOrigBonds (CANON_GLOBALS *pCG, int num_inp_atoms, inp_ATOM *at, int *i, char *szBuf, int buf_len, STRUCT_DATA *sd)` |
| `OrigStruct_FillOut` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:2933` | external | tree-sitter+gcc-aux | `extern int OrigStruct_FillOut (CANON_GLOBALS *pCG, ORIG_ATOM_DATA *orig_inp_data, ORIG_STRUCT *pOrigStruct, STRUCT_DATA *sd)` |
| `OrigStruct_Free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3051` | external | tree-sitter+gcc-aux | `extern void OrigStruct_Free (ORIG_STRUCT *pOrigStruct)` |
| `GetSaveOptLetters` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3100` | internal | tree-sitter+gcc-aux | `static void GetSaveOptLetters (unsigned char save_opt_bits, char *let1, char *let2)` |
| `set_line_separators` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3112` | external | tree-sitter+gcc-aux | `extern void set_line_separators (int bINChIOutputOptions, char **pLF, char **pTAB)` |
| `OutputINCHI_VersionAndKind` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3137` | internal | tree-sitter+gcc-aux | `static int OutputINCHI_VersionAndKind (INCHI_IOSTREAM *out_file, INCHI_IOS_STRING *strbuf, int bINChIOutputOptions, int is_beta, char *pLF, char *pTAB)` |
| `OutputINCHI_MainLayerFormula` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3170` | internal | tree-sitter+gcc-aux | `static int OutputINCHI_MainLayerFormula (CANON_GLOBALS *pCG, INCHI_IOSTREAM *out_file, INCHI_IOS_STRING *strbuf, int *num_components2, int *INCHI_basic_or_INCHI_reconnected, INCHI_OUT_CTL *io, char *pLF, char *pTAB)` |
| `OutputINCHI_MainLayerConnections` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3213` | internal | tree-sitter+gcc-aux | `static int OutputINCHI_MainLayerConnections (CANON_GLOBALS *pCG, INCHI_IOSTREAM *out_file, INCHI_IOS_STRING *strbuf, int *num_components2, int *INCHI_basic_or_INCHI_reconnected, INCHI_OUT_CTL *io, char *pLF, char *pTAB)` |
| `OutputINCHI_MainLayerHydrogens` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3251` | internal | tree-sitter+gcc-aux | `static int OutputINCHI_MainLayerHydrogens (CANON_GLOBALS *pCG, INCHI_IOSTREAM *out_file, INCHI_IOS_STRING *strbuf, int *num_components2, int *INCHI_basic_or_INCHI_reconnected, INCHI_OUT_CTL *io, char *pLF, char *pTAB)` |
| `OutputINCHI_ChargeAndRemovedAddedProtonsLayers` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3293` | internal | tree-sitter+gcc-aux | `static int OutputINCHI_ChargeAndRemovedAddedProtonsLayers (CANON_GLOBALS *pCG, INCHI_IOSTREAM *out_file, INCHI_IOS_STRING *strbuf, INCHI_OUT_CTL *io, char *pLF, char *pTAB)` |
| `OutputINCHI_StereoLayer` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3354` | internal | tree-sitter+gcc-aux | `static int OutputINCHI_StereoLayer (CANON_GLOBALS *pCG, INCHI_IOSTREAM *out_file, INCHI_IOS_STRING *strbuf, INCHI_OUT_CTL *io, char *pLF, char *pTAB)` |
| `OutputINCHI_IsotopicLayer` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3502` | internal | tree-sitter+gcc-aux | `static int OutputINCHI_IsotopicLayer (CANON_GLOBALS *pCG, INCHI_IOSTREAM *out_file, INCHI_IOS_STRING *strbuf, int *INCHI_basic_or_INCHI_reconnected, INCHI_OUT_CTL *io, char *pLF, char *pTAB)` |
| `OutputINCHI_FixedHLayerWithSublayers` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3750` | internal | tree-sitter+gcc-aux | `static int OutputINCHI_FixedHLayerWithSublayers (CANON_GLOBALS *pCG, INCHI_IOSTREAM *out_file, INCHI_IOS_STRING *strbuf, int *INCHI_basic_or_INCHI_reconnected, INCHI_OUT_CTL *io, char *pLF, char *pTAB, int *then_goto_repeat)` |
| `OutputINCHI_PolymerLayer` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:3880` | internal | tree-sitter+gcc-aux | `static int OutputINCHI_PolymerLayer (CANON_GLOBALS *pCG, INCHI_IOSTREAM *out_file, INCHI_IOS_STRING *strbuf, int *INCHI_basic_or_INCHI_reconnected, ORIG_ATOM_DATA *orig_inp_data, ORIG_STRUCT *pOrigStruct, INCHI_OUT_CTL *io, char *pLF, char *pTAB)` |
| `OutputINCHI_PolymerLayer_SingleUnit` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:4111` | internal | tree-sitter+gcc-aux | `static int OutputINCHI_PolymerLayer_SingleUnit (OAD_PolymerUnit *u, int bPolymers, int total_star_atoms, int *n_used_stars, OAD_AtProps *aprops, int *cano_nums, ORIG_ATOM_DATA *orig_inp_data, ORIG_STRUCT *pOrigStruct, INCHI_IOS_STRING *strbuf)` |
| `OutputAUXINFO_HeaderAndNormalization_type` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:4339` | internal | tree-sitter+gcc-aux | `static int OutputAUXINFO_HeaderAndNormalization_type (CANON_GLOBALS *pCG, INCHI_IOSTREAM *out_file, INCHI_IOS_STRING *strbuf, int bINChIOutputOptions, int *INCHI_basic_or_INCHI_reconnected, int *num_components2, INCHI_OUT_CTL *io, char *pLF, char *pTAB)` |
| `OutputAUXINFO_OriginalNumbersAndEquivalenceClasses` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:4396` | internal | tree-sitter+gcc-aux | `static int OutputAUXINFO_OriginalNumbersAndEquivalenceClasses (CANON_GLOBALS *pCG, INCHI_IOSTREAM *out_file, INCHI_IOS_STRING *strbuf, int *num_components2, INCHI_OUT_CTL *io, char *pLF, char *pTAB)` |
| `OutputAUXINFO_TautomericGroupsEquivalence` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:4461` | internal | tree-sitter+gcc-aux | `static int OutputAUXINFO_TautomericGroupsEquivalence (CANON_GLOBALS *pCG, INCHI_IOSTREAM *out_file, INCHI_IOS_STRING *strbuf, INCHI_OUT_CTL *io)` |
| `OutputAUXINFO_Stereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:4500` | internal | tree-sitter+gcc-aux | `static int OutputAUXINFO_Stereo (CANON_GLOBALS *pCG, INCHI_IOSTREAM *out_file, INCHI_IOS_STRING *strbuf, INCHI_OUT_CTL *io, char *pLF, char *pTAB)` |
| `OutputAUXINFO_IsotopicInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:4560` | internal | tree-sitter+gcc-aux | `static int OutputAUXINFO_IsotopicInfo (CANON_GLOBALS *pCG, INCHI_IOSTREAM *out_file, INCHI_IOS_STRING *strbuf, int *INCHI_basic_or_INCHI_reconnected, INCHI_OUT_CTL *io, char *pLF, char *pTAB)` |
| `OutputAUXINFO_ChargesRadicalsAndUnusualValences` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:4717` | internal | tree-sitter+gcc-aux | `static int OutputAUXINFO_ChargesRadicalsAndUnusualValences (CANON_GLOBALS *pCG, INCHI_IOSTREAM *out_file, INCHI_IOS_STRING *strbuf, INCHI_OUT_CTL *io, char *pLF, char *pTAB)` |
| `OutputAUXINFO_ReversibilityInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:4753` | internal | tree-sitter+gcc-aux | `static int OutputAUXINFO_ReversibilityInfo (CANON_GLOBALS *pCG, INCHI_IOSTREAM *out_file, INCHI_IOS_STRING *strbuf, ORIG_STRUCT *pOrigStruct, INCHI_OUT_CTL *io, char *pLF, char *pTAB)` |
| `OutputAUXINFO_PolymerInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:4923` | internal | tree-sitter+gcc-aux | `static int OutputAUXINFO_PolymerInfo (CANON_GLOBALS *pCG, INCHI_IOSTREAM *out_file, INCHI_IOS_STRING *strbuf, ORIG_STRUCT *pOrigStruct, INCHI_OUT_CTL *io, char *pLF, char *pTAB)` |
| `IsBondAtomNumsLesser` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:5007` | internal | tree-sitter+gcc-aux | `static int IsBondAtomNumsLesser (int *bond1, int *bond2)` |
| `EditINCHI_HidePolymerZz` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:5039` | external | tree-sitter+gcc-aux | `extern void EditINCHI_HidePolymerZz (INCHI_IOSTREAM *out, int n_pzz, int n_zy)` |
| `CountPseudoElementInFormula` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:5279` | internal | tree-sitter+gcc-aux | `static int CountPseudoElementInFormula (const char *pseudo, char *s)` |
| `InternallyGetCanoNumsAndComponentNums` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:5362` | internal | tree-sitter+gcc-aux | `static int InternallyGetCanoNumsAndComponentNums (CANON_GLOBALS *pCG, INCHI_IOS_STRING *strbuf, INCHI_OUT_CTL *io, int nat, int *cano_nums, int *compnt_nums)` |
| `MergeZzInHillFormula` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:5435` | external | tree-sitter+gcc-aux | `extern int MergeZzInHillFormula (INCHI_IOS_STRING *strbuf)` |
| `MergeZzInStrHillFormulaComponent` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:5495` | internal | tree-sitter+gcc-aux | `static void MergeZzInStrHillFormulaComponent (char *s)` |
| `inchi_sort_int_pair_ascending` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt1.c:5525` | internal | tree-sitter+gcc-aux | `static void inchi_sort_int_pair_ascending (int *a, int *b)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c`

Parse errors: `0`. Function definitions: `27`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `Eql_INChI_Stereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:54` | external | tree-sitter+gcc-aux | `extern int Eql_INChI_Stereo (INChI_Stereo *s1, int eql1, INChI_Stereo *s2, int eql2, int bRelRac)` |
| `Eql_INChI_Isotopic` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:203` | external | tree-sitter+gcc-aux | `extern int Eql_INChI_Isotopic (INChI *i1, INChI *i2)` |
| `Eql_INChI_Aux_Equ` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:227` | external | tree-sitter+gcc-aux | `extern int Eql_INChI_Aux_Equ (INChI_Aux *a1, int eql1, INChI_Aux *a2, int eql2)` |
| `Eql_INChI_Aux_Num` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:305` | external | tree-sitter+gcc-aux | `extern int Eql_INChI_Aux_Num (INChI_Aux *a1, int eql1, INChI_Aux *a2, int eql2)` |
| `bHasOrigInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:369` | external | tree-sitter+gcc-aux | `extern int bHasOrigInfo (ORIG_INFO *OrigInfo, int num_atoms)` |
| `EqlOrigInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:387` | external | tree-sitter+gcc-aux | `extern int EqlOrigInfo (INChI_Aux *a1, INChI_Aux *a2)` |
| `bHasEquString` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:398` | external | tree-sitter+gcc-aux | `extern int bHasEquString (AT_NUMB *LinearCT, int nLenCT)` |
| `MakeMult` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:431` | external | tree-sitter+gcc-aux | `extern int MakeMult (int mult, const char *szTailingDelim, INCHI_IOS_STRING *buf, int nCtMode, int *bOverflow)` |
| `MakeDelim` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:476` | external | tree-sitter+gcc-aux | `extern int MakeDelim (const char *szTailingDelim, INCHI_IOS_STRING *buf, int *bOverflow)` |
| `MakeEqStr` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:506` | external | tree-sitter+gcc-aux | `extern int MakeEqStr (const char *szTailingDelim, int mult, INCHI_IOS_STRING *buf, int *bOverflow)` |
| `MakeCtStringNew` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:551` | external | tree-sitter+gcc-aux | `extern int MakeCtStringNew (CANON_GLOBALS *pCG, AT_NUMB *LinearCT, int nLenCT, int bAddDelim, S_CHAR *nNum_H, int num_atoms, INCHI_IOS_STRING *strbuf, int nCtMode, int *bOverflow)` |
| `MakeCtStringOld` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:703` | external | tree-sitter+gcc-aux | `extern int MakeCtStringOld (AT_NUMB *LinearCT, int nLenCT, int bAddDelim, INCHI_IOS_STRING *strbuf, int nCtMode, int *bOverflow)` |
| `MakeHString` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:789` | external | tree-sitter+gcc-aux | `extern int MakeHString (int bAddDelim, S_CHAR *LinearCT, int nLenCT, INCHI_IOS_STRING *strbuf, int nCtMode, int *bOverflow)` |
| `MakeCtString` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:1069` | external | tree-sitter+gcc-aux | `extern int MakeCtString (CANON_GLOBALS *pCG, AT_NUMB *LinearCT, int nLenCT, int bAddDelim, S_CHAR *nNum_H, int num_atoms, INCHI_IOS_STRING *strbuf, int nCtMode, int *bOverflow)` |
| `MakeTautString` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:1111` | external | tree-sitter+gcc-aux | `extern int MakeTautString (AT_NUMB *LinearCT, int nLenCT, int bAddDelim, INCHI_IOS_STRING *strbuf, int nCtMode, int *bOverflow)` |
| `MakeCRVString` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:1308` | external | tree-sitter+gcc-aux | `extern int MakeCRVString (ORIG_INFO *OrigInfo, int nLenCT, int bAddDelim, INCHI_IOS_STRING *strbuf, int nCtMode, int *bOverflow)` |
| `MakeEquString` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:1559` | external | tree-sitter+gcc-aux | `extern int MakeEquString (AT_NUMB *LinearCT, int nLenCT, int bAddDelim, INCHI_IOS_STRING *strbuf, int nCtMode, int *bOverflow)` |
| `MakeIsoAtomString` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:1645` | external | tree-sitter+gcc-aux | `extern int MakeIsoAtomString (INChI_IsotopicAtom *IsotopicAtom, int nNumberOfIsotopicAtoms, INCHI_IOS_STRING *strbuf, int nCtMode, int *bOverflow)` |
| `MakeIsoTautString` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:1803` | external | tree-sitter+gcc-aux | `extern int MakeIsoTautString (INChI_IsotopicTGroup *IsotopicTGroup, int nNumberOfIsotopicTGroups, INCHI_IOS_STRING *strbuf, int nCtMode, int *bOverflow)` |
| `MakeIsoHString` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:1924` | external | tree-sitter+gcc-aux | `extern int MakeIsoHString (int *num_iso_H, INCHI_IOS_STRING *strbuf, int nCtMode, int *bOverflow)` |
| `MakeStereoString` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:2019` | external | tree-sitter+gcc-aux | `extern int MakeStereoString (AT_NUMB *at1, AT_NUMB *at2, S_CHAR *parity, int bAddDelim, int nLenCT, INCHI_IOS_STRING *strbuf, int nCtMode, int *bOverflow)` |
| `MakeAbcNumber` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:2143` | external | tree-sitter+gcc-aux | `extern int MakeAbcNumber (char *szString, int nStringLen, const char *szLeadingDelim, int nValue)` |
| `abctol` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:2203` | internal | tree-sitter+gcc-aux | `static long int abctol (const char *szString, char **q)` |
| `inchi_strtol` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:2269` | external | tree-sitter+gcc-aux | `extern long int inchi_strtol (const char *str, const char **p, int base)` |
| `inchi_strtod` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:2288` | external | tree-sitter+gcc-aux | `extern double inchi_strtod (const char *str, const char **p)` |
| `MakeDecNumber` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:2302` | external | tree-sitter+gcc-aux | `extern int MakeDecNumber (char *szString, int nStringLen, const char *szLeadingDelim, int nValue)` |
| `print_sequence_of_nums_compressing_ranges` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt2.c:2376` | external | tree-sitter+gcc-aux | `extern void print_sequence_of_nums_compressing_ranges (int n, int *num, INCHI_IOS_STRING *strbuf)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c`

Parse errors: `11`. Function definitions: `26`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `str_HillFormula` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:60` | external | tree-sitter+gcc-aux | `extern int str_HillFormula (INCHI_SORT *pINChISort, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int num_components, int bUseMulipliers)` |
| `str_HillFormula2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:137` | external | tree-sitter+gcc-aux | `extern int str_HillFormula2 (INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int num_components, int bUseMulipliers)` |
| `str_Connections` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:224` | external | tree-sitter+gcc-aux | `extern int str_Connections (CANON_GLOBALS *pCG, INCHI_SORT *pINChISort, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int ATOM_MODE, int num_components, int bUseMulipliers)` |
| `str_H_atoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:309` | external | tree-sitter+gcc-aux | `extern int str_H_atoms (INCHI_SORT *pINChISort, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int ATOM_MODE, int TAUT_MODE, int num_components, int bUseMulipliers)` |
| `str_Charge2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:427` | external | tree-sitter+gcc-aux | `extern int str_Charge2 (INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int num_components, int bSecondNonTautPass, int bOmitRepetitions, int bUseMulipliers)` |
| `str_FixedH_atoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:618` | external | tree-sitter+gcc-aux | `extern int str_FixedH_atoms (INCHI_SORT *pINChISort, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int ATOM_MODE, int num_components, int bUseMulipliers)` |
| `str_Sp2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:725` | external | tree-sitter+gcc-aux | `extern int str_Sp2 (INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int TAUT_MODE, int num_components, int bSecondNonTautPass, int bOmitRepetitions, int bUseMulipliers)` |
| `str_Sp3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:957` | external | tree-sitter+gcc-aux | `extern int str_Sp3 (INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int TAUT_MODE, int num_components, int bRelRac, int bSecondNonTautPass, int bOmitRepetitions, int bUseMulipliers)` |
| `str_StereoAbsInv` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:1191` | external | tree-sitter+gcc-aux | `extern int str_StereoAbsInv (INCHI_SORT *pINChISort, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int num_components)` |
| `str_IsoAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:1229` | external | tree-sitter+gcc-aux | `extern int str_IsoAtoms (INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int TAUT_MODE, int num_components, int bAbcNumbers, int bSecondNonTautPass, int bOmitRepetitions, int bUseMulipliers)` |
| `str_IsoSp2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:1479` | external | tree-sitter+gcc-aux | `extern int str_IsoSp2 (INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int TAUT_MODE, int num_components, int bSecondNonTautPass, int bOmitRepetitions, int bUseMulipliers)` |
| `str_IsoSp3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:1768` | external | tree-sitter+gcc-aux | `extern int str_IsoSp3 (INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int TAUT_MODE, int num_components, int bRelRac, int bSecondNonTautPass, int bOmitRepetitions, int bUseMulipliers)` |
| `str_IsoStereoAbsInv` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:2058` | external | tree-sitter+gcc-aux | `extern int str_IsoStereoAbsInv (INCHI_SORT *pINChISort, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int num_components)` |
| `str_AuxEqu` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:2106` | external | tree-sitter+gcc-aux | `extern int str_AuxEqu (INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int TAUT_MODE, int num_components, int bSecondNonTautPass, int bOmitRepetitions, int bUseMulipliers)` |
| `str_AuxInvSp3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:2315` | external | tree-sitter+gcc-aux | `extern int str_AuxInvSp3 (INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int TAUT_MODE, int num_components, int bSecondNonTautPass, int bOmitRepetitions, int bUseMulipliers)` |
| `str_AuxInvSp3Numb` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:2584` | external | tree-sitter+gcc-aux | `extern int str_AuxInvSp3Numb (CANON_GLOBALS *pCG, INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int TAUT_MODE, int num_components, int bSecondNonTautPass, int bOmitRepetitions)` |
| `str_AuxIsoNumb` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:2769` | external | tree-sitter+gcc-aux | `extern int str_AuxIsoNumb (CANON_GLOBALS *pCG, INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int TAUT_MODE, int num_components, int bSecondNonTautPass, int bOmitRepetitions)` |
| `str_AuxIsoEqu` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:2962` | external | tree-sitter+gcc-aux | `extern int str_AuxIsoEqu (INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int TAUT_MODE, int num_components, int bSecondNonTautPass, int bOmitRepetitions, int bUseMulipliers)` |
| `str_AuxInvIsoSp3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:3213` | external | tree-sitter+gcc-aux | `extern int str_AuxInvIsoSp3 (INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int TAUT_MODE, int num_components, int bSecondNonTautPass, int bOmitRepetitions, int bUseMulipliers)` |
| `str_AuxInvIsoSp3Numb` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:3575` | external | tree-sitter+gcc-aux | `extern int str_AuxInvIsoSp3Numb (CANON_GLOBALS *pCG, INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int TAUT_MODE, int num_components, int bSecondNonTautPass, int bOmitRepetitions)` |
| `str_AuxNumb` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:3818` | external | tree-sitter+gcc-aux | `extern int str_AuxNumb (CANON_GLOBALS *pCG, INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int TAUT_MODE, int num_components, int bSecondNonTautPass, int bOmitRepetitions)` |
| `str_AuxTgroupEqu` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:3953` | external | tree-sitter+gcc-aux | `extern int str_AuxTgroupEqu (INCHI_SORT *pINChISort, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int TAUT_MODE, int num_components, int bUseMulipliers)` |
| `str_AuxChargeRadVal` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:4013` | external | tree-sitter+gcc-aux | `extern int str_AuxChargeRadVal (INCHI_SORT *pINChISort, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int TAUT_MODE, int num_components, int bUseMulipliers)` |
| `bin_AuxTautTrans` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:4091` | external | tree-sitter+gcc-aux | `extern int bin_AuxTautTrans (INCHI_SORT *pINChISort, INCHI_SORT *pINChISort2, AT_NUMB **pTrans_n, AT_NUMB **pTrans_s, int bOutType, int num_components)` |
| `str_AuxTautTrans` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:4187` | external | tree-sitter+gcc-aux | `extern int str_AuxTautTrans (CANON_GLOBALS *pCG, AT_NUMB *nTrans_n, AT_NUMB *nTrans_s, INCHI_IOS_STRING *strbuf, int *bOverflow, int TAUT_MODE, int num_components)` |
| `str_AuxIsoTgroupEqu` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiprt3.c:4236` | external | tree-sitter+gcc-aux | `extern int str_AuxIsoTgroupEqu (INCHI_SORT *pINChISort, INCHI_IOS_STRING *strbuf, int *bOverflow, int bOutType, int TAUT_MODE, int num_components, int bOmitRepetitions, int bUseMulipliers)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c`

Parse errors: `2`. Function definitions: `16`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `bIsCenterPointStrict` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:215` | external | tree-sitter+gcc-aux | `extern int bIsCenterPointStrict (inp_ATOM *atom, int iat)` |
| `nGet14TautIn7MembAltRing` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:240` | external | tree-sitter+gcc-aux | `extern int nGet14TautIn7MembAltRing (struct tagCANON_GLOBALS *pCG, inp_ATOM *atom, int nStartAtom, int nStartAtomNeighbor, int nStartAtomNeighborEndpoint, int nStartAtomNeighborNeighborEndpoint, AT_RANK *nDfsPathPos, DFS_PATH *DfsPath, int nMaxLenDfsPath, T_ENDPOINT *EndPoint, int nMaxNumEndPoint, T_BONDPOS *BondPos, int nMaxNumBondPos, int *pnNumEndPoint, int *pnNumBondPos, struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD, int num_atoms)` |
| `nGet14TautIn5MembAltRing` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:283` | external | tree-sitter+gcc-aux | `extern int nGet14TautIn5MembAltRing (struct tagCANON_GLOBALS *pCG, inp_ATOM *atom, int nStartAtom, int nStartAtomNeighbor, int nStartAtomNeighborEndpoint, int nStartAtomNeighborNeighborEndpoint, AT_RANK *nDfsPathPos, DFS_PATH *DfsPath, int nMaxLenDfsPath, T_ENDPOINT *EndPoint, int nMaxNumEndPoint, T_BONDPOS *BondPos, int nMaxNumBondPos, int *pnNumEndPoint, int *pnNumBondPos, struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD, int num_atoms)` |
| `nGet12TautIn5MembAltRing` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:322` | external | tree-sitter+gcc-aux | `extern int nGet12TautIn5MembAltRing (struct tagCANON_GLOBALS *pCG, inp_ATOM *atom, int nStartAtom, int nStartAtomNeighbor, AT_RANK *nDfsPathPos, DFS_PATH *DfsPath, int nMaxLenDfsPath, T_ENDPOINT *EndPoint, int nMaxNumEndPoint, T_BONDPOS *BondPos, int nMaxNumBondPos, int *pnNumEndPoint, int *pnNumBondPos, struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD, int num_atoms)` |
| `nGet15TautIn6MembAltRing` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:361` | external | tree-sitter+gcc-aux | `extern int nGet15TautIn6MembAltRing (struct tagCANON_GLOBALS *pCG, inp_ATOM *atom, int nStartAtom, AT_RANK *nDfsPathPos, DFS_PATH *DfsPath, int nMaxLenDfsPath, T_ENDPOINT *EndPoint, int nMaxNumEndPoint, T_BONDPOS *BondPos, int nMaxNumBondPos, int *pnNumEndPoint, int *pnNumBondPos, struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD, int num_atoms)` |
| `nGet15TautInAltPath` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:407` | external | tree-sitter+gcc-aux | `extern int nGet15TautInAltPath (struct tagCANON_GLOBALS *pCG, inp_ATOM *atom, int nStartAtom, AT_RANK *nDfsPathPos, DFS_PATH *DfsPath, int nMaxLenDfsPath, T_ENDPOINT *EndPoint, int nMaxNumEndPoint, T_BONDPOS *BondPos, int nMaxNumBondPos, int *pnNumEndPoint, int *pnNumBondPos, struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD, int num_atoms)` |
| `DFS_FindTautInARing` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:458` | external | tree-sitter+gcc-aux | `extern int DFS_FindTautInARing (struct tagCANON_GLOBALS *pCG, inp_ATOM *atom, int nStartAtom, int nStartAtomNeighbor, int nStartAtomNeighbor2, int nStartAtomNeighborNeighbor, int nCycleLen, AT_RANK *nDfsPathPos, DFS_PATH *DfsPath, CHECK_DFS_RING (*CheckDfsRing), CHECK_CENTERPOINT (*CheckCenterPoint), T_ENDPOINT *EndPoint, int nMaxNumEndPoint, T_BONDPOS *BondPos, int nMaxNumBondPos, int *pnNumEndPoint, int *pnNumBondPos, struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD, int num_atoms)` |
| `DFS_FindTautAltPath` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:616` | external | tree-sitter+gcc-aux | `extern int DFS_FindTautAltPath (struct tagCANON_GLOBALS *pCG, inp_ATOM *atom, int nStartAtom, int nStartAtomNeighbor, int nStartAtomNeighbor2, int nStartAtomNeighborNeighbor, int nCycleLen, AT_RANK *nDfsPathPos, DFS_PATH *DfsPath, CHECK_DFS_PATH (*CheckDfsPath), CHECK_DFS_CENTERPOINT (*CheckCenterPointPath), T_ENDPOINT *EndPoint, int nMaxNumEndPoint, T_BONDPOS *BondPos, int nMaxNumBondPos, int *pnNumEndPoint, int *pnNumBondPos, struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD, int num_atoms)` |
| `are_alt_bonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:778` | external | tree-sitter+gcc-aux | `extern int are_alt_bonds (U_CHAR *bonds, int len)` |
| `AddBondsPos` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:846` | external | tree-sitter+gcc-aux | `extern int AddBondsPos (inp_ATOM *atom, T_BONDPOS *BondPosTmp, int nNumBondPosTmp, T_BONDPOS *BondPos, int nMaxNumBondPos, int nNumBondPos)` |
| `AddEndPoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:897` | external | tree-sitter+gcc-aux | `extern int AddEndPoints (T_ENDPOINT *EndPointTmp, int nNumNewEndPoint, T_ENDPOINT *EndPoint, int nMaxNumEndPoint, int nNumEndPoint)` |
| `Check7MembTautRing` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:951` | external | tree-sitter+gcc-aux | `extern int Check7MembTautRing (struct tagCANON_GLOBALS *pCG, inp_ATOM *atom, DFS_PATH *DfsPath, int nLenDfsPath, int nStartAtomNeighbor, int nStartAtomNeighbor2, int nStartAtomNeighborNeighbor, T_ENDPOINT *EndPoint, int nMaxNumEndPoint, T_BONDPOS *BondPos, int nMaxNumBondPos, int *pnNumEndPoint, int *pnNumBondPos, struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD, int num_atoms)` |
| `Check6MembTautRing` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:1142` | external | tree-sitter+gcc-aux | `extern int Check6MembTautRing (struct tagCANON_GLOBALS *pCG, inp_ATOM *atom, DFS_PATH *DfsPath, int nLenDfsPath, int nStartAtomNeighbor, int nStartAtomNeighbor2, int nStartAtomNeighborNeighbor, T_ENDPOINT *EndPoint, int nMaxNumEndPoint, T_BONDPOS *BondPos, int nMaxNumBondPos, int *pnNumEndPoint, int *pnNumBondPos, struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD, int num_atoms)` |
| `Check15TautPathCenterpoint` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:1409` | external | tree-sitter+gcc-aux | `extern int Check15TautPathCenterpoint (inp_ATOM *atom, DFS_PATH *DfsPath, int nLenDfsPath, int jNxtNeigh, struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD, int num_atoms)` |
| `Check15TautPath` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:1434` | external | tree-sitter+gcc-aux | `extern int Check15TautPath (struct tagCANON_GLOBALS *pCG, inp_ATOM *atom, DFS_PATH *DfsPath, int nLenDfsPath, int jNxtNeigh, int nStartAtomNeighbor, int nStartAtomNeighbor2, int nStartAtomNeighborNeighbor, T_ENDPOINT *EndPoint, int nMaxNumEndPoint, T_BONDPOS *BondPos, int nMaxNumBondPos, int *pnNumEndPoint, int *pnNumBondPos, struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD, int num_atoms)` |
| `Check5MembTautRing` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiqueu.c:1666` | external | tree-sitter+gcc-aux | `extern int Check5MembTautRing (struct tagCANON_GLOBALS *pCG, inp_ATOM *atom, DFS_PATH *DfsPath, int nLenDfsPath, int nStartAtomNeighbor, int nStartAtomNeighbor2, int nStartAtomNeighborNeighbor, T_ENDPOINT *EndPoint, int nMaxNumEndPoint, T_BONDPOS *BondPos, int nMaxNumBondPos, int *pnNumEndPoint, int *pnNumBondPos, struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD, int num_atoms)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c`

Parse errors: `87`. Function definitions: `85`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `getInchiStateReadErr` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:536` | internal | tree-sitter+gcc-aux | `static void getInchiStateReadErr (int stat, char *szMsg)` |
| `getInchiErrName` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:562` | internal | tree-sitter+gcc-aux | `static const char *getInchiErrName (int nErr)` |
| `SetHillFormFromInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:583` | internal | tree-sitter+gcc-aux | `static int SetHillFormFromInChI (InpInChI *OneInput)` |
| `ReadWriteInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:615` | external | tree-sitter+gcc-aux | `extern int ReadWriteInChI (INCHI_CLOCK *ic, struct tagCANON_GLOBALS *pCG, INCHI_IOSTREAM *pInp, INCHI_IOSTREAM *pOut, INCHI_IOSTREAM *pLog, INPUT_PARMS *ip_inp, STRUCT_DATA *sd_inp, inp_ATOM **at, int *num_at, int *num_bonds, OAD_Polymer **polymer, OAD_V3000 **v3000, char *szMsg, int nMsgLen, long unsigned int (*WarningFlags)[2])` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1003` | external | tree-sitter | `INCHI_HEAPCHK if (ret < 0)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1195` | external | tree-sitter | `INCHI_HEAPCHK if (sd_inp)` |
| `OutputInChIAsRequested` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1213` | internal | tree-sitter+gcc-aux | `static int OutputInChIAsRequested (struct tagCANON_GLOBALS *pCG, INCHI_IOSTREAM *pOut, INCHI_IOSTREAM *pLog, const INPUT_PARMS *ip_inp, STRUCT_DATA *sd_inp, InpInChI *OneInput, int *num_components, MODE_PIXH *nModeProtonIsoExchgH, long int num_inp, unsigned char save_opt_bits)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1264` | external | tree-sitter | `INCHI_HEAPCHK if (num_components[INCHI_BAS])` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1335` | external | tree-sitter | `INCHI_HEAPCHK if (k1 \|\| k2 /*\|\| !pStr*/)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1392` | external | tree-sitter | `INCHI_HEAPCHK /* take care of protons in AuxInfo */ if (nModeProtonIsoExchgH[iINChI] == MODE_PIXH_ADD_TO_EACH && j == TAUT_YES)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1630` | external | tree-sitter | `INCHI_HEAPCHK if (nRet1 == _IS_FATAL \|\| nRet1 == _IS_ERROR)` |
| `GetNumNeighborsFromInchi` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1646` | external | tree-sitter+gcc-aux | `extern int GetNumNeighborsFromInchi (INChI *pInChI, AT_NUMB nAtNumber)` |
| `CountStereoTypes` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1723` | external | tree-sitter+gcc-aux | `extern int CountStereoTypes (INChI *pInChI, int *num_known_SB, int *num_known_SC, int *num_unk_und_SB, int *num_unk_und_SC, int *num_SC_PIII, int *num_SC_AsIII)` |
| `bInpInchiComponentExists` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1800` | external | tree-sitter+gcc-aux | `extern int bInpInchiComponentExists (InpInChI *pOneInput, int iInChI, int bMobileH, int k)` |
| `bInpInchiComponentDeleted` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1819` | external | tree-sitter+gcc-aux | `extern int bInpInchiComponentDeleted (InpInChI *pOneInput, int iInChI, int bMobileH, int k)` |
| `bRevInchiComponentExists` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1838` | external | tree-sitter+gcc-aux | `extern int bRevInchiComponentExists (StrFromINChI *pStruct, int iInChI, int bMobileH, int k)` |
| `bRevInchiComponentDeleted` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1859` | external | tree-sitter+gcc-aux | `extern int bRevInchiComponentDeleted (StrFromINChI *pStruct, int iInChI, int bMobileH, int k)` |
| `DetectInpInchiCreationOptions` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1880` | external | tree-sitter+gcc-aux | `extern int DetectInpInchiCreationOptions (InpInChI *pOneInput, int *bHasReconnected, int *bHasMetal, int *bHasFixedH, int *nModeFlagsStereo, int *bTautFlagsStereo)` |
| `bInChIHasReconnectedMetal` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:1993` | internal | tree-sitter+gcc-aux | `static int bInChIHasReconnectedMetal (INChI *pInChI)` |
| `SetProtonsAndXchgIsoH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:2015` | internal | tree-sitter+gcc-aux | `static int SetProtonsAndXchgIsoH (int bInChI2Structure, int bReqSplitOutputInChI, int bReqProtonsForEachComponent, int bReqNonTaut, int bReqStereo, int *num_components, MODE_PIXH *nModeProtonIsoExchgH, InpInChI *OneInput)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:2176` | external | tree-sitter | `INCHI_HEAPCHK /* remove unneeded Stereo and/or Fixed H */ if (!bReqStereo)` |
| `GetInChIFormulaNumH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:2215` | internal | tree-sitter+gcc-aux | `static int GetInChIFormulaNumH (INChI *pInChI, int *nNumH)` |
| `GetInChINumH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:2254` | internal | tree-sitter+gcc-aux | `static int GetInChINumH (INChI *pInChI, int *nNumH)` |
| `GetInChIIsoH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:2296` | internal | tree-sitter+gcc-aux | `static int GetInChIIsoH (INChI *pInChI, int *nNumIsotopicH)` |
| `InChILine2Data` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:2342` | external | tree-sitter+gcc-aux | `extern int InChILine2Data (INCHI_IOSTREAM *pInp, SEGM_LINE *pLine, char **pStr, int *pState, int *nErr, INChI *(*pInpInChI)[2], int (*nNumComponents)[2], REM_PROTONS (*nNumProtons)[2], int (*s)[2][2], int bReadCoord, int bInchi2Struct, INCHI_MODE nMode, int *bStdFormat, int *input_has_save_opt, unsigned char *input_save_opt_bits, OAD_Polymer **ppolymer, OAD_V3000 **pv3000)` |
| `bIsoMayBeArranged` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:3633` | external | tree-sitter | `int bIsoMayBeArranged(int bInchi2Struct, int iso_diff[NUM_H_ISOTOPES], REM_PROTONS nNumProtons[INCHI_NUM][TAUT_NUM], INChI* pInpInChI[INCHI_NUM][TAUT_NUM], int nNumComponents[INCHI_NUM][TAUT_NUM], int iINChI)` |
| `ParseAuxSegmentVersion` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:3783` | internal | tree-sitter+gcc-aux | `static int ParseAuxSegmentVersion (const char *str, int bMobileH, INChI **pInpInChI, int *ppnNumComponents, int state)` |
| `CopyAtomNumbers` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:3806` | internal | tree-sitter+gcc-aux | `static int CopyAtomNumbers (INChI *pInChI_To, int bIsoTo, INChI *pInChI_From, int bIsoFrom)` |
| `ParseAuxSegmentNumbers` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:3847` | internal | tree-sitter+gcc-aux | `static int ParseAuxSegmentNumbers (const char *str, int bMobileH, INChI **pInpInChI, int *ppnNumComponents, int state, int *pbAbc)` |
| `ParseAuxSegmentAtomEqu` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:4154` | internal | tree-sitter+gcc-aux | `static int ParseAuxSegmentAtomEqu (const char *str, int bMobileH, INChI **pInpInChI, int *ppnNumComponents, int state)` |
| `ParseAuxSegmentGroupEqu` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:4213` | internal | tree-sitter+gcc-aux | `static int ParseAuxSegmentGroupEqu (const char *str, int bMobileH, INChI **pInpInChI, int *ppnNumComponents, int state)` |
| `ParseAuxSegmentSp3Inv` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:4252` | internal | tree-sitter+gcc-aux | `static int ParseAuxSegmentSp3Inv (const char *str, int bMobileH, INChI **pInpInChI, int *ppnNumComponents, int state)` |
| `ParseAuxSegmentSp3InvNumbers` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:4312` | internal | tree-sitter+gcc-aux | `static int ParseAuxSegmentSp3InvNumbers (const char *str, int bMobileH, INChI **pInpInChI, int *ppnNumComponents, int state)` |
| `ParseAuxSegmentReverseCRV` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:4371` | internal | tree-sitter+gcc-aux | `static int ParseAuxSegmentReverseCRV (const char *str, int state)` |
| `ParseAuxSegmentReverseAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:4392` | internal | tree-sitter+gcc-aux | `static int ParseAuxSegmentReverseAtoms (const char *str, int state)` |
| `ParseAuxSegmentReverseBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:4413` | internal | tree-sitter+gcc-aux | `static int ParseAuxSegmentReverseBonds (const char *str, int state)` |
| `ParseAuxSegmentReverseXYZ` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:4432` | internal | tree-sitter+gcc-aux | `static int ParseAuxSegmentReverseXYZ (const char *str, XYZ_COORD **ppXYZ, int state)` |
| `AddAuxSegmentCoord` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:4491` | internal | tree-sitter+gcc-aux | `static int AddAuxSegmentCoord (int nRet, XYZ_COORD *pXYZ, int nLenXYZ, INChI *(*pInpInChI)[2], int (*nNumComponents)[2])` |
| `ReadInChICoord` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:4683` | external | tree-sitter+gcc-aux | `extern int ReadInChICoord (INCHI_IOSTREAM *pInp, SEGM_LINE *pLine, int *pState, INChI *(*pInpInChI)[2], int (*nNumComponents)[2])` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:4704` | external | tree-sitter | `INCHI_HEAPCHK /* Get "InChI=1/" */ if (pLine->len)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:4739` | external | tree-sitter | `INCHI_HEAPCHK if (ret < 0)` |
| `ReadInChILine` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:5016` | internal | tree-sitter+gcc-aux | `static int ReadInChILine (INCHI_IOSTREAM *pInp, SEGM_LINE *pLine, char **pStr, int *pState, INChI *(*pInpInChI)[2], int (*nNumComponents)[2], REM_PROTONS (*nNumProtons)[2], int (*s)[2][2], int *bStdFormat, int *input_has_save_opt, unsigned char *input_save_opt_bits, int bInchi2Struct, OAD_Polymer **ppPolymer, OAD_V3000 **ppV3000)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:5071` | external | tree-sitter | `INCHI_HEAPCHK if (pLine->str && (pLine->len == 0 \|\| (c != SEG_END && c != RI_ERR_EOF) \|\| !(p = strstr(pLine->str, "InChI=1")))) /* djb-rwth: fixing a NULL pointer dereference; addressing LLVM warning; ignoring LLVM warning: value used */` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:5117` | external | tree-sitter | `INCHI_HEAPCHK if (ret < 0)` |
| `ParseSegmentIsoExchgH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:5541` | internal | tree-sitter+gcc-aux | `static int ParseSegmentIsoExchgH (const char *str, int bMobileH, REM_PROTONS *nNumProtons, int *pnNumComponents, int state, int *pbAbc)` |
| `ParseSegmentPerm` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:5663` | internal | tree-sitter+gcc-aux | `static int ParseSegmentPerm (const char *str, int bMobileH, INChI **pInpInChI, int *ppnNumComponents, int state, int *pbAbc)` |
| `ParseSegmentIsoAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:5806` | internal | tree-sitter+gcc-aux | `static int ParseSegmentIsoAtoms (const char *str, int bMobileH, INChI **pInpInChI, int *ppnNumComponents, int state, int *pbAbc)` |
| `ParseSegmentSp3s` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:6235` | internal | tree-sitter+gcc-aux | `static int ParseSegmentSp3s (const char *str, int bMobileH, INChI **pInpInChI, int (*s)[2], int *ppnNumComponents, int state)` |
| `bIsSp3LayerNotEmpty` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:6346` | internal | tree-sitter+gcc-aux | `static int bIsSp3LayerNotEmpty (INChI **pInpInChI, int bMobileH, int bIso, int nNumComponents)` |
| `ParseSegmentSp3m` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:6378` | internal | tree-sitter+gcc-aux | `static int ParseSegmentSp3m (const char *str, int bMobileH, INChI **pInpInChI, int *ppnNumComponents, int state)` |
| `ParseSegmentSp3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:6600` | internal | tree-sitter+gcc-aux | `static int ParseSegmentSp3 (const char *str, int bMobileH, INChI **pInpInChI, int *ppnNumComponents, int state, int *pbAbc)` |
| `ParseSegmentSp2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:6866` | internal | tree-sitter+gcc-aux | `static int ParseSegmentSp2 (const char *str, int bMobileH, INChI **pInpInChI, int *ppnNumComponents, int state, int *pbAbc)` |
| `ParseSegmentProtons` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:7426` | internal | tree-sitter+gcc-aux | `static int ParseSegmentProtons (const char *str, int bMobileH, REM_PROTONS *nNumProtons, int *ppnNumComponents)` |
| `ParseSegmentPolymer` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:7494` | internal | tree-sitter+gcc-aux | `static int ParseSegmentPolymer (const char *str, int bMobileH, REM_PROTONS *nNumProtons, int *ppnNumComponents, int na_total, int nb_total, int bInchi2Struct, OAD_Polymer **ppPolymer, OAD_V3000 **ppV3000)` |
| `ParseSegmentReadDelimitedNumbers` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:7966` | external | tree-sitter+gcc-aux | `extern const char *ParseSegmentReadDelimitedNumbers (const char *str, const char *pEnd, INT_ARRAY *numlist, char c_delim, char c_stop, int *ret)` |
| `ParseSegmentCharge` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:8032` | internal | tree-sitter+gcc-aux | `static int ParseSegmentCharge (const char *str, int bMobileH, INChI **pInpInChI, int *ppnNumComponents)` |
| `ParseSegmentMobileH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:8230` | internal | tree-sitter+gcc-aux | `static int ParseSegmentMobileH (const char *str, int bMobileH, INChI **pInpInChI, int *pnNumComponents, int *pbAbc)` |
| `ParseSegmentConnections` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:9084` | internal | tree-sitter+gcc-aux | `static int ParseSegmentConnections (const char *str, int bMobileH, INChI **pInpInChI, int *pnNumComponents, int *pbAbc, int *nb_total)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:9387` | external | tree-sitter | `INCHI_HEAPCHK if (pInChI[iComponent + i].nAtom)` |
| `nFillOutProtonMobileH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:9718` | internal | tree-sitter+gcc-aux | `static int nFillOutProtonMobileH (INChI *pInChI)` |
| `nProtonCopyIsotopicInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:9770` | internal | tree-sitter+gcc-aux | `static int nProtonCopyIsotopicInfo (INChI *pInChI_to, INChI *pInChI_from)` |
| `ParseSegmentFormula` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:9810` | internal | tree-sitter+gcc-aux | `static int ParseSegmentFormula (const char *str, int bMobileH, INChI **pInpInChI, int *pnNumComponents, int *na_total)` |
| `CopySegment` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:10263` | internal | tree-sitter+gcc-aux | `static int CopySegment (INChI *pInChITo, INChI *pInChIFrom, int SegmentType, int bIsotopicTo, int bIsotopicFrom)` |
| `insertions_sort_AT_NUMB` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:10546` | external | tree-sitter+gcc-aux | `extern int insertions_sort_AT_NUMB (AT_NUMB *base, int num)` |
| `getInChIChar` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:10565` | internal | tree-sitter+gcc-aux | `static int getInChIChar (INCHI_IOSTREAM *pInp)` |
| `AddInChIChar` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:10603` | internal | tree-sitter+gcc-aux | `static int AddInChIChar (INCHI_IOSTREAM *pInp, SEGM_LINE *Line, const char *pszToken)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:10614` | external | tree-sitter | `INCHI_HEAPCHK if (Line->len + 2 >= Line->len_alloc)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:10619` | external | tree-sitter | `INCHI_HEAPCHK if (str)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:10641` | external | tree-sitter | `INCHI_HEAPCHK if (c < 0)` |
| `nGetInChISegment` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:10676` | internal | tree-sitter+gcc-aux | `static int nGetInChISegment (INCHI_IOSTREAM *pInp, SEGM_LINE *Line, const char *pszToken)` |
| `AddLinkedBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:10699` | internal | tree-sitter+gcc-aux | `static int AddLinkedBond (AT_NUMB at1, AT_NUMB at2, AT_NUMB num_at, LINKED_BONDS *pLB)` |
| `PrepareSaveOptBits` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:10771` | external | tree-sitter+gcc-aux | `extern void PrepareSaveOptBits (INPUT_PARMS *ip, INCHI_IOSTREAM *pLog, const const long int num_inp, const char *szCurHdr, int input_has_save_opt, unsigned char input_save_opt_bits, unsigned char *save_opt_bits)` |
| `TreatErrorsInReadInChIString` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:10901` | external | tree-sitter+gcc-aux | `extern void TreatErrorsInReadInChIString (int nReadStatus, int nErr, int pState, INPUT_PARMS *ip, INCHI_IOSTREAM *pOut, INCHI_IOSTREAM *pLog, long int *num_inp, long int *num_errors, long int *num_processed, char **pstrHdr, char **pszCurHdr, InpInChI *pOneInput)` |
| `ConvertInChI2InChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:11003` | external | tree-sitter+gcc-aux | `extern int ConvertInChI2InChI (INPUT_PARMS *ip, InpInChI *pOneInput, INCHI_IOSTREAM *pOut, INCHI_IOSTREAM *pLog, STRUCT_DATA *sd, int *num_components, MODE_PIXH *nModeProtonIsoExchgH, char **pszCurHdr, long int num_inp, long int *num_errors, unsigned char save_opt_bits, inchiTime *pulTStart, long int *ulProcessingTime, struct tagINCHI_CLOCK *ic, struct tagCANON_GLOBALS *pCG)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:11176` | external | tree-sitter | `else if (*pszCurHdr && (*pszCurHdr)[0])` |
| `ConvertInChI2Struct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:11203` | external | tree-sitter+gcc-aux | `extern int ConvertInChI2Struct (const INPUT_PARMS *ip_inp, INPUT_PARMS *ip, InpInChI *pOneInput, inp_ATOM **at, int *num_at, OAD_Polymer **polymer, OAD_V3000 **v3000, INCHI_IOSTREAM *pOut, INCHI_IOSTREAM *pLog, STRUCT_DATA *sd, int *num_components, MODE_PIXH *nModeProtonIsoExchgH, char **pszCurHdr, char *szMsg, int nMsgLen, char *szMessage, int nInitLenMessage, int nMessageLen, int input_is_stdinchi, int bHasSomeReconnected, int bHasSomeFixedH, int bHasMetal, int nModeFlagsStereo, int bTautFlags, int bReqNonTaut, long unsigned int (*WarningFlags)[2], long int num_inp, long int *num_errors, unsigned char save_opt_bits, inchiTime *pulTStart, long int *ulProcessingTime, struct tagINCHI_CLOCK *ic, struct tagCANON_GLOBALS *pCG)` |
| `DetectAndExposePolymerInternals` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:11629` | external | tree-sitter+gcc-aux | `extern int DetectAndExposePolymerInternals (INCHI_IOSTREAM *is)` |
| `DetectHiddenPolymerStuff` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:12006` | external | tree-sitter+gcc-aux | `extern int DetectHiddenPolymerStuff (char *tmpstr, int tmpstrlen, int *ninsert, int *insert_pos, int insert_lead_offset, int *nstars)` |
| `SegmentSp3CreateEmpty` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:12060` | internal | tree-sitter+gcc-aux | `static int SegmentSp3CreateEmpty (const char *str, int bMobileH, INChI **pInpInChI, int nNumComponents, int state, int *pbAbc)` |
| `SegmentSp3StoreStereoCenters` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:12123` | internal | tree-sitter+gcc-aux | `static int SegmentSp3StoreStereoCenters (int *pbAbc, const char *pStart, const char *pEnd, int pInChI_iComponent_nNumberOfAtoms, INChI_Stereo *PStereo_0)` |
| `SegmentSp3CopyMultiplierCovered` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:12203` | internal | tree-sitter+gcc-aux | `static int SegmentSp3CopyMultiplierCovered (int mpy_component, int iComponent, INChI *pInChI, int bIso, int nCpyType)` |
| `SegmentSp3ProcessAbbreviation` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:12239` | internal | tree-sitter+gcc-aux | `static int SegmentSp3ProcessAbbreviation (int *mpy_component, int iComponent, int nNumComponents, int val, const char *q, int state, int *pbAbc, int bMobileH, int nCpyType, INChI *pInChI, INChI *pInpInChI_ALT_TAUT_bMobileH)` |
| `extract_from_inchi_string` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:12416` | internal | tree-sitter+gcc-aux | `static int extract_from_inchi_string (char *sinchi, InpInChI *OneInput)` |
| `extract_stereo_info_from_inchi_string` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:12499` | external | tree-sitter+gcc-aux | `extern int extract_stereo_info_from_inchi_string (char *sinchi, int nat, int *orig, int *at_stereo_mark_orig)` |
| `extract_all_backbone_bonds_from_inchi_string` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiread.c:12547` | external | tree-sitter+gcc-aux | `extern int extract_all_backbone_bonds_from_inchi_string (char *sinchi, int *n_all_bkb_orig, int *orig, int *all_bkb_orig)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c`

Parse errors: `0`. Function definitions: `15`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `QueueCreate` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:67` | external | tree-sitter+gcc-aux | `extern QUEUE *QueueCreate (int nTotLength, int nSize)` |
| `QueueAdd` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:87` | external | tree-sitter+gcc-aux | `extern int QueueAdd (QUEUE *q, qInt *Val)` |
| `QueueGet` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:101` | external | tree-sitter+gcc-aux | `extern int QueueGet (QUEUE *q, qInt *Val)` |
| `QueueGetAny` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:122` | external | tree-sitter+gcc-aux | `extern int QueueGetAny (QUEUE *q, qInt *Val, int ord)` |
| `QueueCreate` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:139` | external | tree-sitter | `QUEUE *QueueCreate( int nTotLength, int nSize )` |
| `QueueAdd` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:155` | external | tree-sitter | `int QueueAdd( QUEUE *q, QINT_TYPE *Val )` |
| `QueueGet` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:169` | external | tree-sitter | `int QueueGet( QUEUE *q, QINT_TYPE *Val )` |
| `QueueGetAny` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:186` | external | tree-sitter | `int QueueGetAny( QUEUE *q, QINT_TYPE *Val, int ord )` |
| `QueueDelete` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:203` | external | tree-sitter+gcc-aux | `extern QUEUE *QueueDelete (QUEUE *q)` |
| `QueueReinit` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:216` | external | tree-sitter+gcc-aux | `extern int QueueReinit (QUEUE *q)` |
| `QueueLength` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:231` | external | tree-sitter+gcc-aux | `extern int QueueLength (QUEUE *q)` |
| `QueueWrittenLength` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:245` | external | tree-sitter+gcc-aux | `extern int QueueWrittenLength (QUEUE *q)` |
| `GetMinRingSize` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:262` | external | tree-sitter+gcc-aux | `extern int GetMinRingSize (inp_ATOM *atom, QUEUE *q, AT_RANK *nAtomLevel, S_CHAR *cSource, AT_RANK nMaxRingSize)` |
| `is_bond_in_Nmax_memb_ring` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:362` | external | tree-sitter+gcc-aux | `extern int is_bond_in_Nmax_memb_ring (inp_ATOM *atom, int at_no, int neigh_ord, QUEUE *q, AT_RANK *nAtomLevel, S_CHAR *cSource, AT_RANK nMaxRingSize)` |
| `is_atom_in_3memb_ring` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichiring.c:420` | external | tree-sitter+gcc-aux | `extern int is_atom_in_3memb_ring (inp_ATOM *atom, int at_no)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c`

Parse errors: `38`. Function definitions: `54`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `clear_t_group_info` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:509` | external | tree-sitter+gcc-aux | `extern void clear_t_group_info (T_GROUP_INFO *ti)` |
| `GetTgroupInfoFromInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:575` | external | tree-sitter+gcc-aux | `extern int GetTgroupInfoFromInChI (T_GROUP_INFO *ti, inp_ATOM *at, AT_NUMB *endpoint, INChI *pInChI)` |
| `FillOutpStructEndpointFromInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:692` | external | tree-sitter+gcc-aux | `extern int FillOutpStructEndpointFromInChI (INChI *pInChI, AT_NUMB **pEndpoint)` |
| `cmp_charge_val` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:732` | external | tree-sitter+gcc-aux | `extern int cmp_charge_val (const void *a1, const void *a2, void *p)` |
| `bMayBeACationInMobileHLayer` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:749` | external | tree-sitter+gcc-aux | `extern int bMayBeACationInMobileHLayer (inp_ATOM *at, VAL_AT *pVA, int iat, int bMobileH)` |
| `clean_charge_val` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:792` | external | tree-sitter+gcc-aux | `extern int clean_charge_val (struct tagCANON_GLOBALS *pCG, CHARGE_VAL *pChargeVal, int len, inp_ATOM *atom, VAL_AT *pVA, int iat, int bIsMetal, int bMobileH, AT_NUMB *endpoint)` |
| `GetAtomRestoreInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:912` | external | tree-sitter+gcc-aux | `extern int GetAtomRestoreInfo (struct tagCANON_GLOBALS *pCG, inp_ATOM *atom, int iat, VAL_AT *pVArray, const SRM *pSrm, int bMobileH, AT_NUMB *endpoint)` |
| `get_bonds_valences` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:1382` | external | tree-sitter | `int get_bonds_valences( int nPeriodicNum, int bonds_valence, int num_H, VAL_AT *pVA )` |
| `get_sp_element_type` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:1414` | external | tree-sitter+gcc-aux | `extern int get_sp_element_type (int nPeriodicNumber, int *nRow)` |
| `ReallocTCGroups` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:1502` | external | tree-sitter+gcc-aux | `extern int ReallocTCGroups (ALL_TC_GROUPS *pTCGroups, int nAdd)` |
| `RegisterTCGroup` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:1527` | external | tree-sitter+gcc-aux | `extern int RegisterTCGroup (ALL_TC_GROUPS *pTCGroups, int nGroupType, int nGroupOrdNum, int nVertexCap, int nVertexFlow, int nEdgeCap, int nEdgeFlow, int nNumEdges)` |
| `nTautEndpointEdgeCap` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:1574` | external | tree-sitter+gcc-aux | `extern int nTautEndpointEdgeCap (inp_ATOM *at, VAL_AT *pVA, int i)` |
| `BondFlowMaxcapMinorder` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:1690` | external | tree-sitter+gcc-aux | `extern int BondFlowMaxcapMinorder (inp_ATOM *atom, VAL_AT *pVA, const SRM *pSrm, int iat, int ineigh, int *pnMaxcap, int *pnMinorder, int *pbNeedsFlower)` |
| `AtomStcapStflow` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:1790` | external | tree-sitter+gcc-aux | `extern int AtomStcapStflow (inp_ATOM *atom, VAL_AT *pVA, const SRM *pSrm, int iat, int *pnStcap, int *pnStflow, EdgeFlow *pnMGroupEdgeCap, EdgeFlow *pnMGroupEdgeFlow)` |
| `nCountBnsSizes` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:1870` | external | tree-sitter+gcc-aux | `extern int nCountBnsSizes (inp_ATOM *at, int num_at, int nAddEdges2eachAtom, int nAddVertices, T_GROUP_INFO *ti, VAL_AT *pVA, const SRM *pSrm, ALL_TC_GROUPS *pTCGroups)` |
| `nAddSuperCGroups` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:2279` | external | tree-sitter+gcc-aux | `extern int nAddSuperCGroups (ALL_TC_GROUPS *pTCGroups)` |
| `AddTGroups2TCGBnStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:2425` | external | tree-sitter+gcc-aux | `extern int AddTGroups2TCGBnStruct (BN_STRUCT *pBNS, StrFromINChI *pStruct, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int nMaxAddEdges)` |
| `ConnectTwoVertices` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:2627` | external | tree-sitter+gcc-aux | `extern int ConnectTwoVertices (BNS_VERTEX *p1, BNS_VERTEX *p2, BNS_EDGE *e, BN_STRUCT *pBNS, int bClearEdge)` |
| `AddRadicalToMetal` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:2706` | external | tree-sitter+gcc-aux | `extern int AddRadicalToMetal (int *tot_st_cap, int *tot_st_flow, const SRM *pSrm, BN_STRUCT *pBNS, ALL_TC_GROUPS *pTCGroups)` |
| `ConnectMetalFlower` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:2743` | external | tree-sitter+gcc-aux | `extern int ConnectMetalFlower (int *pcur_num_vertices, int *pcur_num_edges, int *tot_st_cap, int *tot_st_flow, const SRM *pSrm, BN_STRUCT *pBNS, ALL_TC_GROUPS *pTCGroups)` |
| `SetEdgeCapFlow` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:2911` | external | tree-sitter+gcc-aux | `extern void SetEdgeCapFlow (BNS_EDGE *e, int edge_cap, int edge_flow)` |
| `AddEdgeFlow` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:2923` | external | tree-sitter+gcc-aux | `extern int AddEdgeFlow (int edge_cap, int edge_flow, BNS_EDGE *e01, BNS_VERTEX *pSrc, BNS_VERTEX *pDst, int *tot_st_cap, int *tot_st_flow)` |
| `ConnectSuperCGroup` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:3058` | external | tree-sitter+gcc-aux | `extern int ConnectSuperCGroup (int nSuperCGroup, int *nAddGroups, int num_add, int *pcur_num_vertices, int *pcur_num_edges, int *tot_st_cap, int *tot_st_flow, BN_STRUCT *pBNS, ALL_TC_GROUPS *pTCGroups)` |
| `AddStCapFlow` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:3218` | external | tree-sitter+gcc-aux | `extern void AddStCapFlow (BNS_VERTEX *vert_ficpoint, int *tot_st_flow, int *tot_st_cap, int cap, int flow)` |
| `SetStCapFlow` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:3231` | external | tree-sitter+gcc-aux | `extern void SetStCapFlow (BNS_VERTEX *vert_ficpoint, int *tot_st_flow, int *tot_st_cap, int cap, int flow)` |
| `AddCGroups2TCGBnStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:3249` | external | tree-sitter+gcc-aux | `extern int AddCGroups2TCGBnStruct (BN_STRUCT *pBNS, StrFromINChI *pStruct, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int nMaxAddEdges)` |
| `nNumEdgesToCnVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:3824` | external | tree-sitter+gcc-aux | `extern int nNumEdgesToCnVertex (const C_NODE *pCN, int len, int v)` |
| `AllocateAndInitTCGBnStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:3879` | external | tree-sitter+gcc-aux | `extern BN_STRUCT *AllocateAndInitTCGBnStruct (StrFromINChI *pStruct, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int nMaxAddAtoms, int nMaxAddEdges, int max_altp, int *pNum_changed_bonds)` |
| `IncrZeroBondsAndClearEndpts` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4152` | external | tree-sitter+gcc-aux | `extern void IncrZeroBondsAndClearEndpts (inp_ATOM *at, int num_at, int iComponent)` |
| `IncrZeroBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4172` | external | tree-sitter+gcc-aux | `extern void IncrZeroBonds (inp_ATOM *at, int num_at, int iComponent)` |
| `ClearEndpts` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4191` | external | tree-sitter+gcc-aux | `extern void ClearEndpts (inp_ATOM *at, int num_at)` |
| `GetDeltaChargeFromVF` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4229` | external | tree-sitter+gcc-aux | `extern int GetDeltaChargeFromVF (BN_STRUCT *pBNS, VAL_AT *pVA, VF *vf)` |
| `EvaluateChargeChanges` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4320` | external | tree-sitter+gcc-aux | `extern int EvaluateChargeChanges (BN_STRUCT *pBNS, VAL_AT *pVA, int *pnDeltaH, int *pnDeltaCharge, int *pnNumVisitedAtoms)` |
| `RunBnsTestOnce` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4507` | external | tree-sitter+gcc-aux | `extern int RunBnsTestOnce (BN_STRUCT *pBNS, BN_DATA *pBD, VAL_AT *pVA, Vertex *pvFirst, Vertex *pvLast, int *pPathLen, int *pnDeltaH, int *pnDeltaCharge, int *pnNumVisitedAtoms)` |
| `RunBnsRestoreOnce` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4546` | external | tree-sitter+gcc-aux | `extern int RunBnsRestoreOnce (BN_STRUCT *pBNS, BN_DATA *pBD, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups)` |
| `comp_cc_cand` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4584` | external | tree-sitter+gcc-aux | `extern int comp_cc_cand (const void *a1, const void *a2)` |
| `get_pVA_atom_type` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4628` | external | tree-sitter+gcc-aux | `extern int get_pVA_atom_type (VAL_AT *pVA, inp_ATOM *at, int iat, int bond_type)` |
| `AllocEdgeList` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4690` | external | tree-sitter+gcc-aux | `extern int AllocEdgeList (EDGE_LIST *pEdges, int nLen)` |
| `AddToEdgeList` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4738` | external | tree-sitter+gcc-aux | `extern int AddToEdgeList (EDGE_LIST *pEdges, int iedge, int nAddLen)` |
| `RemoveFromEdgeListByIndex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4759` | external | tree-sitter+gcc-aux | `extern int RemoveFromEdgeListByIndex (EDGE_LIST *pEdges, int index)` |
| `FindInEdgeList` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4777` | external | tree-sitter+gcc-aux | `extern int FindInEdgeList (EDGE_LIST *pEdges, int iedge)` |
| `RemoveFromEdgeListByValue` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4794` | external | tree-sitter+gcc-aux | `extern int RemoveFromEdgeListByValue (EDGE_LIST *pEdges, int iedge)` |
| `AllocBfsQueue` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4815` | external | tree-sitter+gcc-aux | `extern int AllocBfsQueue (BFS_Q *pQ, int num_at, int min_ring_size)` |
| `RemoveForbiddenEdgeMask` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4869` | external | tree-sitter+gcc-aux | `extern void RemoveForbiddenEdgeMask (BN_STRUCT *pBNS, EDGE_LIST *pEdges, int forbidden_edge_mask)` |
| `SetForbiddenEdgeMask` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4880` | external | tree-sitter+gcc-aux | `extern void SetForbiddenEdgeMask (BN_STRUCT *pBNS, EDGE_LIST *pEdges, int forbidden_edge_mask)` |
| `RemoveForbiddenBondFlowBits` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4891` | external | tree-sitter+gcc-aux | `extern void RemoveForbiddenBondFlowBits (BN_STRUCT *pBNS, int forbidden_edge_mask_int)` |
| `GetChargeFlowerUpperEdge` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:4915` | external | tree-sitter+gcc-aux | `extern int GetChargeFlowerUpperEdge (BN_STRUCT *pBNS, VAL_AT *pVA, int nChargeEdge)` |
| `NormalizeStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:5013` | external | tree-sitter | `int NormalizeStructure( ICHICONST INPUT_PARMS *ip, STRUCT_DATA *sd, BN_STRUCT *pBNS, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, inp_ATOM *at_norm, inp_ATOM *at_fixed_bonds_out, T_GROUP_INFO *t_group_info )` |
| `MakeOneInChIOutOfStrFromINChI2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:5087` | external | tree-sitter+gcc-aux | `extern int MakeOneInChIOutOfStrFromINChI2 (struct tagCANON_GLOBALS *pCG, INCHI_CLOCK *ic, const INPUT_PARMS *ip_inp, STRUCT_DATA *sd_inp, BN_STRUCT *pBNS, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, inp_ATOM *at3, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, T_GROUP_INFO **t_group_info, inp_ATOM **at_norm, inp_ATOM **at_prep)` |
| `MakeOneInChIOutOfStrFromINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:5168` | external | tree-sitter+gcc-aux | `extern int MakeOneInChIOutOfStrFromINChI (struct tagCANON_GLOBALS *pCG, INCHI_CLOCK *ic, const INPUT_PARMS *ip, STRUCT_DATA *sd, StrFromINChI *pStruct, inp_ATOM *at2, inp_ATOM *at3, ALL_TC_GROUPS *pTCGroups)` |
| `ConnectDisconnectedH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:5480` | external | tree-sitter+gcc-aux | `extern int ConnectDisconnectedH (inp_ATOM *at, int num_atoms, int num_deleted_H)` |
| `DisconnectedConnectedH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:5594` | external | tree-sitter+gcc-aux | `extern int DisconnectedConnectedH (inp_ATOM *at, int num_atoms, int num_deleted_H)` |
| `MakeInChIOutOfStrFromINChI2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:5669` | external | tree-sitter+gcc-aux | `extern int MakeInChIOutOfStrFromINChI2 (INCHI_CLOCK *ic, CANON_GLOBALS *pCG, const INPUT_PARMS *ip_inp, STRUCT_DATA *sd_inp, StrFromINChI *pStruct, int iComponent, int iAtNoOffset, long int num_inp)` |
| `OutputInChIOutOfStrFromINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr1.c:5888` | external | tree-sitter+gcc-aux | `extern int OutputInChIOutOfStrFromINChI (struct tagINCHI_CLOCK *ic, struct tagCANON_GLOBALS *pCG, const INPUT_PARMS *ip_inp, STRUCT_DATA *sd_inp, long int num_inp, int bINChIOutputOptions, INCHI_IOSTREAM *pout, INCHI_IOSTREAM *plog, InpInChI *pOneInput, int bHasSomeFixedH, unsigned char save_opt_bits)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c`

Parse errors: `8`. Function definitions: `32`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `CopyAt2St` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:57` | external | tree-sitter+gcc-aux | `extern void CopyAt2St (inp_ATOM *at, inp_ATOM_STEREO *st, int num_atoms)` |
| `CopySt2At` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:79` | external | tree-sitter+gcc-aux | `extern void CopySt2At (inp_ATOM *at, inp_ATOM_STEREO *st, int num_atoms)` |
| `RestoreAtomConnectionsSetStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:106` | external | tree-sitter+gcc-aux | `extern int RestoreAtomConnectionsSetStereo (StrFromINChI *pStruct, int iComponent, int iAtNoOffset, INChI *pInChI, INChI *pInChIMobH)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:196` | external | tree-sitter | `INCHI_HEAPCHK /* isotopic atoms */ if (pInChI->IsotopicAtom && pInChI->nNumberOfIsotopicAtoms)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:414` | external | tree-sitter | `INCHI_HEAPCHK if (nNumDeletedH)` |
| `SetStereoBondTypeFor0DParity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:790` | external | tree-sitter+gcc-aux | `extern int SetStereoBondTypeFor0DParity (inp_ATOM *at, int i1, int m1)` |
| `SetStereoBondTypesFrom0DStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:844` | external | tree-sitter+gcc-aux | `extern int SetStereoBondTypesFrom0DStereo (StrFromINChI *pStruct, INChI *pInChI)` |
| `CopyBnsToAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:979` | external | tree-sitter+gcc-aux | `extern int CopyBnsToAtom (StrFromINChI *pStruct, BN_STRUCT *pBNS, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int bAllowZeroBondOrder)` |
| `CheckBnsConsistency` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:1155` | external | tree-sitter+gcc-aux | `extern int CheckBnsConsistency (StrFromINChI *pStruct, BN_STRUCT *pBNS, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int bNoRad)` |
| `AddExplicitDeletedH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:1470` | external | tree-sitter+gcc-aux | `extern int AddExplicitDeletedH (inp_ATOM *at, int jv, int num_at, int *iDeletedH, int *iH, int nNumDeletedH, int bTwoStereo)` |
| `bFindCumuleneChain` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:1565` | external | tree-sitter+gcc-aux | `extern int bFindCumuleneChain (inp_ATOM *at, AT_NUMB i1, AT_NUMB i2, AT_NUMB *nCumulene, int nMaxLen)` |
| `set_bond_type` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:1614` | external | tree-sitter+gcc-aux | `extern int set_bond_type (inp_ATOM *at, AT_NUMB i1, AT_NUMB i2, int bType)` |
| `set_cumulene_0D_parity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:1639` | external | tree-sitter+gcc-aux | `extern int set_cumulene_0D_parity (inp_ATOM *at, inp_ATOM_STEREO *st, int num_at, int idelH1, int i1, int i2, int idelH2, int parity, int len)` |
| `set_atom_0D_parity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:1773` | external | tree-sitter+gcc-aux | `extern int set_atom_0D_parity (inp_ATOM *at, inp_ATOM_STEREO *st, int num_at, int num_deleted_H, int i1, int parity)` |
| `MoveRadToAtomsAddCharges` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:1841` | external | tree-sitter+gcc-aux | `extern int MoveRadToAtomsAddCharges (BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int forbidden_mask)` |
| `AdjustTgroupsToForbiddenEdges2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:2134` | external | tree-sitter+gcc-aux | `extern int AdjustTgroupsToForbiddenEdges2 (BN_STRUCT *pBNS, inp_ATOM *at, VAL_AT *pVA, int num_atoms, int forbidden_mask)` |
| `RearrangePlusMinusEdgesFlow` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:3403` | external | tree-sitter+gcc-aux | `extern int RearrangePlusMinusEdgesFlow (BN_STRUCT *pBNS, BN_DATA *pBD, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int forbidden_edge_mask)` |
| `IncrementZeroOrderBondsToHeteroat` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:3520` | external | tree-sitter+gcc-aux | `extern int IncrementZeroOrderBondsToHeteroat (BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)` |
| `MovePlusFromS2DiaminoCarbon` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:3736` | external | tree-sitter+gcc-aux | `extern int MovePlusFromS2DiaminoCarbon (BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)` |
| `EliminateChargeSeparationOnHeteroatoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:3884` | external | tree-sitter+gcc-aux | `extern int EliminateChargeSeparationOnHeteroatoms (BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask, int forbidden_stereo_edge_mask)` |
| `MoveChargeFromHeteroatomsToMetals` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:4138` | external | tree-sitter | `int MoveChargeFromHeteroatomsToMetals( BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask ) /********* Avoid charge separation on heteroatoms ******************/` |
| `RestoreCyanoGroup` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:4554` | external | tree-sitter+gcc-aux | `extern int RestoreCyanoGroup (BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)` |
| `RestoreIsoCyanoGroup` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:4663` | external | tree-sitter+gcc-aux | `extern int RestoreIsoCyanoGroup (BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)` |
| `FixMetal_Nminus_Ominus` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:4965` | external | tree-sitter+gcc-aux | `extern int FixMetal_Nminus_Ominus (BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)` |
| `RestoreNNNgroup` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:5115` | external | tree-sitter+gcc-aux | `extern int RestoreNNNgroup (BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)` |
| `EliminateNitrogen5Val3Bonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:5810` | external | tree-sitter+gcc-aux | `extern int EliminateNitrogen5Val3Bonds (BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)` |
| `Convert_SIV_to_SVI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:6054` | external | tree-sitter+gcc-aux | `extern int Convert_SIV_to_SVI (BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)` |
| `PlusFromDB_N_DB_O_to_Metal` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:6270` | external | tree-sitter+gcc-aux | `extern int PlusFromDB_N_DB_O_to_Metal (BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)` |
| `MoveMobileHToAvoidFixedBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:6446` | external | tree-sitter+gcc-aux | `extern int MoveMobileHToAvoidFixedBonds (BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)` |
| `RemoveRadFromMobileHEndpoint` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:6525` | external | tree-sitter+gcc-aux | `extern int RemoveRadFromMobileHEndpoint (BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)` |
| `RemoveRadFromMobileHEndpointFixH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:6963` | external | tree-sitter+gcc-aux | `extern int RemoveRadFromMobileHEndpointFixH (BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)` |
| `MoveChargeToMakeCenerpoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr2.c:7658` | external | tree-sitter+gcc-aux | `extern int MoveChargeToMakeCenerpoints (BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr3.c`

Parse errors: `15`. Function definitions: `6`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `bHas_N_V` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr3.c:92` | external | tree-sitter+gcc-aux | `extern int bHas_N_V (inp_ATOM *at2, int num_atoms)` |
| `FillTgDiffHChgFH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr3.c:113` | internal | tree-sitter+gcc-aux | `static int FillTgDiffHChgFH (TgDiffHChgFH *tdhc, int max_tdhc, inp_ATOM *at2, inp_ATOM *atf, AT_NUMB *nCanon2AtnoRevrs, VAL_AT *pVA, T_GROUP_INFO *ti, EDGE_LIST *pAtomIndList)` |
| `FixFixedHRestoredStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr3.c:333` | external | tree-sitter+gcc-aux | `extern int FixFixedHRestoredStructure (CANON_GLOBALS *pCG, INCHI_CLOCK *ic, const INPUT_PARMS *ip, STRUCT_DATA *sd, BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, inp_ATOM *at3, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, T_GROUP_INFO **ppt_group_info, inp_ATOM **ppat_norm, inp_ATOM **ppat_prep, INChI **pInChI, long int num_inp, int bHasSomeFixedH, int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask, int forbidden_stereo_edge_mask)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr3.c:489` | external | tree-sitter | `INCHI_HEAPCHK if ((ret = FillOutExtraFixedHDataRestr( pStruct ))) /* djb-rwth: addressing LLVM warning */` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr3.c:494` | external | tree-sitter | `INCHI_HEAPCHK if ((ret = FillOutCMP2FHINCHI( pStruct, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr3.c:3897` | external | tree-sitter | `INCHI_HEAPCHK if (ret < 0)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c`

Parse errors: `10`. Function definitions: `23`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `ForbidCarbonChargeEdges` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:57` | external | tree-sitter+gcc-aux | `extern int ForbidCarbonChargeEdges (BN_STRUCT *pBNS, ALL_TC_GROUPS *pTCGroups, EDGE_LIST *pCarbonChargeEdges, int forbidden_edge_mask)` |
| `ForbidNintrogenPlus2BondsInSmallRings` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:116` | external | tree-sitter+gcc-aux | `extern int ForbidNintrogenPlus2BondsInSmallRings (BN_STRUCT *pBNS, inp_ATOM *at, int num_at, VAL_AT *pVA, int min_ring_size, ALL_TC_GROUPS *pTCGroups, EDGE_LIST *pNplus2BondsEdges, int forbidden_edge_mask)` |
| `FixLessHydrogenInFormula` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:177` | external | tree-sitter+gcc-aux | `extern int FixLessHydrogenInFormula (BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, inp_ATOM *atf, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)` |
| `FixMoreHydrogenInFormula` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:398` | external | tree-sitter+gcc-aux | `extern int FixMoreHydrogenInFormula (BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, inp_ATOM *atf, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)` |
| `FixAddProtonForADP` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:541` | external | tree-sitter | `int FixAddProtonForADP(BN_STRUCT* pBNS, BN_DATA* pBD, StrFromINChI* pStruct, inp_ATOM* at, inp_ATOM* at2, inp_ATOM* atf, VAL_AT* pVA, ALL_TC_GROUPS* pTCGroups, ICR* picr, int* pnNumRunBNS, int* pnTotalDelta, int forbidden_edge_mask)` |
| `FixRemoveExtraTautEndpoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:589` | external | tree-sitter+gcc-aux | `extern int FixRemoveExtraTautEndpoints (BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, inp_ATOM *atf, inp_ATOM *atn, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, ICR *picr, int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)` |
| `FillOutExtraFixedHDataRestr` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:754` | external | tree-sitter+gcc-aux | `extern int FillOutExtraFixedHDataRestr (StrFromINChI *pStruct)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:813` | external | tree-sitter | `INCHI_HEAPCHK if (pStruct->nAtno2Canon[i])` |
| `FillOutExtraFixedHDataInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:831` | external | tree-sitter+gcc-aux | `extern int FillOutExtraFixedHDataInChI (StrFromINChI *pStruct, INChI **pInChI)` |
| `FillOutCMP2FHINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:870` | external | tree-sitter+gcc-aux | `extern int FillOutCMP2FHINCHI (StrFromINChI *pStruct, inp_ATOM *at2, VAL_AT *pVA, INChI **pInChI, CMP2FHINCHI *pc2i)` |
| `FillOutCMP2MHINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:1021` | external | tree-sitter+gcc-aux | `extern int FillOutCMP2MHINCHI (StrFromINChI *pStruct, ALL_TC_GROUPS *pTCGroups, inp_ATOM *at2, VAL_AT *pVA, INChI **pInChI, CMP2MHINCHI *pc2i)` |
| `NormalizeAndCompare` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:1274` | external | tree-sitter+gcc-aux | `extern int NormalizeAndCompare (CANON_GLOBALS *pCG, INCHI_CLOCK *ic, const INPUT_PARMS *ip, STRUCT_DATA *sd, BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, inp_ATOM *at3, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, INChI **pInChI, long int num_inp, int bHasSomeFixedH, int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask, int forbidden_stereo_edge_mask)` |
| `CheckAndRefixStereobonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:1678` | external | tree-sitter+gcc-aux | `extern int CheckAndRefixStereobonds (BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)` |
| `MoveChargeToRemoveCenerpoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:1860` | external | tree-sitter+gcc-aux | `extern int MoveChargeToRemoveCenerpoints (BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)` |
| `MakeSingleBondsMetal2ChargedHeteroat` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:2194` | external | tree-sitter+gcc-aux | `extern int MakeSingleBondsMetal2ChargedHeteroat (BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)` |
| `SaltBondsToCoordBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:2373` | external | tree-sitter+gcc-aux | `extern int SaltBondsToCoordBonds (BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask)` |
| `ForbidMetalCarbonEdges` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:2666` | external | tree-sitter | `int ForbidMetalCarbonEdges(BN_STRUCT* pBNS, inp_ATOM* at, int num_at, VAL_AT* pVA, ALL_TC_GROUPS* pTCGroups, EDGE_LIST* pMetalCarbonEdges, int forbidden_edge_mask)` |
| `RunBnsRestore1` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:2740` | external | tree-sitter+gcc-aux | `extern int RunBnsRestore1 (CANON_GLOBALS *pCG, INCHI_CLOCK *ic, const INPUT_PARMS *ip, STRUCT_DATA *sd, BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, INChI **pInChI, long int num_inp, int bHasSomeFixedH)` |
| `RestoreAtomMakeBNS` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:3198` | external | tree-sitter+gcc-aux | `extern int RestoreAtomMakeBNS (INCHI_CLOCK *ic, CANON_GLOBALS *pCG, const INPUT_PARMS *ip, STRUCT_DATA *sd, StrFromINChI *pStruct, int iComponent, int iAtNoOffset, INChI **pInChI, const char *szCurHdr, long int num_inp, int bHasSomeFixedH)` |
| `OneInChI2Atom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:3477` | external | tree-sitter+gcc-aux | `extern int OneInChI2Atom (INCHI_CLOCK *ic, CANON_GLOBALS *pCG, const INPUT_PARMS *ip_inp, STRUCT_DATA *sd, const char *szCurHdr, long int num_inp, StrFromINChI *pStruct, int iComponent, int iAtNoOffset, int bHasSomeFixedH, INChI **pInChI)` |
| `MakeProtonComponent` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:3574` | external | tree-sitter+gcc-aux | `extern int MakeProtonComponent (StrFromINChI *pStruct, int iComponent, int num_prot)` |
| `AddRemProtonsInRestrStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:3615` | external | tree-sitter+gcc-aux | `extern int AddRemProtonsInRestrStruct (INCHI_CLOCK *ic, CANON_GLOBALS *pCG, const INPUT_PARMS *ip_inp, STRUCT_DATA *sd, long int num_inp, int bHasSomeFixedH, StrFromINChI *pStruct, int num_components, StrFromINChI *pStructR, int num_componentsR, NUM_H *nProtonsToBeRemovedByNormFromRevrs, int *recmet_change_balance)` |
| `AddRemIsoProtonsInRestrStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr4.c:3816` | external | tree-sitter+gcc-aux | `extern int AddRemIsoProtonsInRestrStruct (INCHI_CLOCK *ic, CANON_GLOBALS *pCG, const INPUT_PARMS *ip_inp, STRUCT_DATA *sd, long int num_inp, int bHasSomeFixedH, StrFromINChI *pStruct, int num_components, StrFromINChI *pStructR, int num_componentsR, NUM_H *pProtonBalance, NUM_H *recmet_change_balance)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr5.c`

Parse errors: `15`. Function definitions: `5`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `GetPlusMinusVertex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr5.c:58` | external | tree-sitter+gcc-aux | `extern int GetPlusMinusVertex (BN_STRUCT *pBNS, ALL_TC_GROUPS *pTCGroups, int bCheckForbiddenPlus, int bCheckForbiddenMinus)` |
| `bIsUnsatCarbonInASmallRing` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr5.c:92` | external | tree-sitter+gcc-aux | `extern int bIsUnsatCarbonInASmallRing (inp_ATOM *at2, VAL_AT *pVA, int iat, BFS_Q *pbfsq, int min_ring_size)` |
| `FixMobileHRestoredStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr5.c:141` | external | tree-sitter+gcc-aux | `extern int FixMobileHRestoredStructure (CANON_GLOBALS *pCG, INCHI_CLOCK *ic, const INPUT_PARMS *ip, STRUCT_DATA *sd, BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, inp_ATOM *at3, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, T_GROUP_INFO **ppt_group_info, inp_ATOM **ppat_norm, inp_ATOM **ppat_prep, INChI **pInChI, long int num_inp, int bHasSomeFixedH, int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask, int forbidden_stereo_edge_mask)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr5.c:276` | external | tree-sitter | `INCHI_HEAPCHK if ((ret = FillOutExtraFixedHDataRestr( pStruct ))) /* djb-rwth: addressing LLVM warning */` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr5.c:283` | external | tree-sitter | `INCHI_HEAPCHK if ((ret = FillOutCMP2MHINCHI( pStruct, pTCGroups, at2, pVA, pInChI, pc2i ))) /* djb-rwth: addressing LLVM warning */` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr6.c`

Parse errors: `2`. Function definitions: `1`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `FixRestoredStructureStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr6.c:59` | external | tree-sitter+gcc-aux | `extern int FixRestoredStructureStereo (struct tagCANON_GLOBALS *pCG, INCHI_CLOCK *ic, INCHI_MODE cmpInChI, ICR *icr, INCHI_MODE cmpInChI2, ICR *icr2, const INPUT_PARMS *ip, STRUCT_DATA *sd, BN_STRUCT *pBNS, BN_DATA *pBD, StrFromINChI *pStruct, inp_ATOM *at, inp_ATOM *at2, inp_ATOM *at3, VAL_AT *pVA, ALL_TC_GROUPS *pTCGroups, T_GROUP_INFO **ppt_group_info, inp_ATOM **ppat_norm, inp_ATOM **ppat_prep, INChI **pInChI, long int num_inp, int *pnNumRunBNS, int *pnTotalDelta, int forbidden_edge_mask, int forbidden_stereo_edge_mask)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c`

Parse errors: `9`. Function definitions: `22`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `DisplayStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:78` | external | tree-sitter | `int DisplayStructure(struct tagCANON_GLOBALS* pCG, inp_ATOM* at, int num_at, OAD_Polymer* polymer, int num_removed_H, int bAdd_DT_to_num_H, int nNumRemovedProtons, NUM_H* nNumRemovedProtonsIsotopic, int bIsotopic, int j /*bTautomeric*/, INChI** cur_INChI, INChI_Aux** cur_INChI_Aux, int bAbcNumbers, DRAW_PARMS* dp, INCHI_MODE nMode, char* szTitle)` |
| `InChI2Atom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:101` | external | tree-sitter+gcc-aux | `extern int InChI2Atom (INCHI_CLOCK *ic, CANON_GLOBALS *pCG, const INPUT_PARMS *ip, STRUCT_DATA *sd, const char *szCurHdr, long int num_inp, StrFromINChI *pStruct, int iComponent, int iAtNoOffset, int bI2A_Flag, int bHasSomeFixedH, InpInChI *OneInput)` |
| `RemoveFixHInChIIdentical2MobH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:256` | external | tree-sitter+gcc-aux | `extern void RemoveFixHInChIIdentical2MobH (InpInChI *pOneInput)` |
| `MarkDisconectedIdenticalToReconnected` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:283` | external | tree-sitter+gcc-aux | `extern int MarkDisconectedIdenticalToReconnected (InpInChI *pOneInput)` |
| `SetUpSrm` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:364` | external | tree-sitter+gcc-aux | `extern void SetUpSrm (SRM *pSrm)` |
| `MergeStructureComponents` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:406` | external | tree-sitter+gcc-aux | `extern int MergeStructureComponents (const INPUT_PARMS *ip, STRUCT_DATA *sd, long int num_inp, char *szCurHdr, const SRM *pSrm, int bReqNonTaut, StrFromINChI *(*pStruct)[2], InpInChI *pOneInput)` |
| `DisplayAllRestoredComponents` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:622` | external | tree-sitter | `int DisplayAllRestoredComponents(struct tagCANON_GLOBALS* pCG, inp_ATOM* at, int num_at, const char* szCurHdr)` |
| `DisplayOneRestoredComponent` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:655` | external | tree-sitter | `int DisplayOneRestoredComponent(struct tagCANON_GLOBALS* pCG, StrFromINChI* pStruct, inp_ATOM* at, int iComponent, int nNumComponents, int bMobileH, const char* szCurHdr)` |
| `DisplayRestoredComponent` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:716` | external | tree-sitter | `int DisplayRestoredComponent(struct tagCANON_GLOBALS* pCG, StrFromINChI* pStruct, int iComponent, int iAtNoOffset, INChI* pInChI, const char* szCurHdr)` |
| `DisplayStructureComponents` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:771` | external | tree-sitter | `int DisplayStructureComponents(struct tagCANON_GLOBALS* pCG, ICHICONST INPUT_PARMS* ip, STRUCT_DATA* sd, long num_inp, char* szCurHdr, ICHICONST SRM* pSrm, int bReqNonTaut, StrFromINChI* pStruct[INCHI_NUM][TAUT_NUM], InpInChI* pOneInput)` |
| `AllInchiToStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:1042` | external | tree-sitter+gcc-aux | `extern int AllInchiToStructure (INCHI_CLOCK *ic, CANON_GLOBALS *pCG, const INPUT_PARMS *ip_inp, STRUCT_DATA *sd_inp, long int num_inp, char *szCurHdr, const SRM *pSrm, int bHasSomeFixedH, StrFromINChI *(*pStruct)[2], InpInChI *pOneInput)` |
| `AddProtonAndIsoHBalanceToMobHStruct` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:1180` | external | tree-sitter+gcc-aux | `extern int AddProtonAndIsoHBalanceToMobHStruct (INCHI_CLOCK *ic, CANON_GLOBALS *pCG, const INPUT_PARMS *ip, STRUCT_DATA *sd, long int num_inp, int bHasSomeFixedH, char *szCurHdr, StrFromINChI *(*pStruct)[2], InpInChI *pOneInput)` |
| `FreeStrFromINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:1353` | external | tree-sitter+gcc-aux | `extern void FreeStrFromINChI (StrFromINChI *(*pStruct)[2], int (*nNumComponents)[2])` |
| `FreeInpInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:1438` | external | tree-sitter+gcc-aux | `extern void FreeInpInChI (InpInChI *pOneInput)` |
| `CompareAllOrigInchiToRevInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:1548` | external | tree-sitter+gcc-aux | `extern int CompareAllOrigInchiToRevInChI (StrFromINChI *(*pStruct)[2], InpInChI *pOneInput, int bReqNonTaut, long int num_inp, char *szCurHdr)` |
| `CompareAllDisconnectedOrigInchiToRevInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:1683` | external | tree-sitter+gcc-aux | `extern int CompareAllDisconnectedOrigInchiToRevInChI (StrFromINChI *(*pStruct)[2], InpInChI *pOneInput, int bHasSomeFixedH, long int num_inp, char *szCurHdr)` |
| `CompareTwoPairsOfInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:2203` | external | tree-sitter+gcc-aux | `extern int CompareTwoPairsOfInChI (INChI **pInChI1, INChI **pInChI2, int bMobileH, INCHI_MODE *CompareInchiFlags)` |
| `CompareOneOrigInchiToRevInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:2240` | external | tree-sitter+gcc-aux | `extern int CompareOneOrigInchiToRevInChI (StrFromINChI *pStruct, INChI **pInChI, int bMobileH, int iComponent, long int num_inp, char *szCurHdr, COMPONENT_REM_PROTONS *nCurRemovedProtons, INCHI_MODE *CompareInchiFlags)` |
| `CompareReversedStereoINChI3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:2362` | external | tree-sitter+gcc-aux | `extern INCHI_MODE CompareReversedStereoINChI3 (INChI_Stereo *s1, INChI_Stereo *s2, ICR *picr)` |
| `CompareReversedINChI3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:2644` | external | tree-sitter+gcc-aux | `extern INCHI_MODE CompareReversedINChI3 (INChI *i1, INChI *i2, INChI_Aux *a1, INChI_Aux *a2, int *err)` |
| `AddOneMsg` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:3138` | external | tree-sitter+gcc-aux | `extern int AddOneMsg (char *szMsg, int used_len, int tot_len, const char *szAddMsg, const char *szDelim)` |
| `FillOutCompareMessage` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichirvr7.c:3179` | external | tree-sitter+gcc-aux | `extern int FillOutCompareMessage (char *szMsg, int nLenMsg, INCHI_MODE *bits)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c`

Parse errors: `0`. Function definitions: `29`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `inchi_qsort` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:66` | external | tree-sitter+gcc-aux | `extern void inchi_qsort (void *pParam, void *base, size_t num, size_t width, int (*comp) (const void *, const void *, void *))` |
| `inchi_swap` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:286` | external | tree-sitter+gcc-aux | `extern void inchi_swap (char *a, char *b, size_t width)` |
| `insertions_sort` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:304` | external | tree-sitter+gcc-aux | `extern int insertions_sort (void *pCG, void *base, size_t num, size_t width, int (*compare) (const void *, const void *, void *))` |
| `insertions_sort_AT_NUMBERS` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:331` | external | tree-sitter+gcc-aux | `extern int insertions_sort_AT_NUMBERS (void *pCG, AT_NUMB *base, int num, int (*compare) (const void *, const void *, void *))` |
| `insertions_sort_NeighList_AT_NUMBERS` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:355` | external | tree-sitter+gcc-aux | `extern void insertions_sort_NeighList_AT_NUMBERS (NEIGH_LIST base, AT_RANK *nRank)` |
| `insertions_sort_AT_RANK` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:375` | external | tree-sitter+gcc-aux | `extern int insertions_sort_AT_RANK (AT_RANK *base, int num)` |
| `insertions_sort_NeighList_AT_NUMBERS3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:396` | external | tree-sitter+gcc-aux | `extern int insertions_sort_NeighList_AT_NUMBERS3 (NEIGH_LIST base, AT_RANK *nRank)` |
| `insertions_sort_NeighListBySymmAndCanonRank` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:421` | external | tree-sitter+gcc-aux | `extern void insertions_sort_NeighListBySymmAndCanonRank (NEIGH_LIST base, const AT_RANK *nSymmRank, const AT_RANK *nCanonRank)` |
| `CompNeighborsAT_NUMBER` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:453` | external | tree-sitter+gcc-aux | `extern int CompNeighborsAT_NUMBER (const void *a1, const void *a2, void *p)` |
| `comp_AT_RANK` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:467` | external | tree-sitter+gcc-aux | `extern int comp_AT_RANK (const void *a1, const void *a2, void *p)` |
| `CompRank` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:475` | external | tree-sitter+gcc-aux | `extern int CompRank (const void *a1, const void *a2, void *p)` |
| `CompRanksOrd` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:486` | external | tree-sitter+gcc-aux | `extern int CompRanksOrd (const void *a1, const void *a2, void *p)` |
| `CompAtomInvariants2Only` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:502` | external | tree-sitter+gcc-aux | `extern int CompAtomInvariants2Only (const void *a1, const void *a2, void *p)` |
| `CompAtomInvariants2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:536` | external | tree-sitter+gcc-aux | `extern int CompAtomInvariants2 (const void *a1, const void *a2, void *p)` |
| `CompChemElemLex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:551` | external | tree-sitter+gcc-aux | `extern int CompChemElemLex (const void *a1, const void *a2)` |
| `CompareNeighListLex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:560` | external | tree-sitter+gcc-aux | `extern int CompareNeighListLex (NEIGH_LIST pp1, NEIGH_LIST pp2, const AT_RANK *nRank)` |
| `CompareNeighListLexUpToMaxRank` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:582` | external | tree-sitter+gcc-aux | `extern int CompareNeighListLexUpToMaxRank (NEIGH_LIST pp1, NEIGH_LIST pp2, const AT_RANK *nRank, AT_RANK nMaxAtNeighRank)` |
| `compare_NeighLists` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:607` | external | tree-sitter+gcc-aux | `extern int compare_NeighLists (const NEIGH_LIST *op1, const NEIGH_LIST *op2, void *p)` |
| `CompNeighListRanks` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:615` | external | tree-sitter+gcc-aux | `extern int CompNeighListRanks (const void *a1, const void *a2, void *p)` |
| `CompNeighLists` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:632` | external | tree-sitter+gcc-aux | `extern int CompNeighLists (const void *a1, const void *a2, void *p)` |
| `CompNeighListsUpToMaxRank` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:644` | external | tree-sitter+gcc-aux | `extern int CompNeighListsUpToMaxRank (const void *a1, const void *a2, void *p)` |
| `CompNeighListRanksOrd` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:657` | external | tree-sitter+gcc-aux | `extern int CompNeighListRanksOrd (const void *a1, const void *a2, void *p)` |
| `CompRanksInvOrd` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:670` | external | tree-sitter+gcc-aux | `extern int CompRanksInvOrd (const void *a1, const void *a2, void *p)` |
| `CompNeighborsRanksCountEql` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:677` | external | tree-sitter+gcc-aux | `extern int CompNeighborsRanksCountEql (const void *a1, const void *a2, void *p)` |
| `CreateNeighListFromLinearCT` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:701` | external | tree-sitter+gcc-aux | `extern NEIGH_LIST *CreateNeighListFromLinearCT (AT_NUMB *LinearCT, int nLenCT, int num_atoms)` |
| `CreateNeighList` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:810` | external | tree-sitter+gcc-aux | `extern NEIGH_LIST *CreateNeighList (int num_atoms, int num_at_tg, sp_ATOM *at, int bDoubleBondSquare, T_GROUP_INFO *t_group_info)` |
| `FreeNeighList` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:941` | external | tree-sitter+gcc-aux | `extern void FreeNeighList (NEIGH_LIST *pp)` |
| `BreakAllTies` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:956` | external | tree-sitter+gcc-aux | `extern int BreakAllTies (CANON_GLOBALS *pCG, int num_atoms, int num_max, AT_RANK **pRankStack, NEIGH_LIST *NeighList, AT_RANK *nTempRank, CANON_STAT *pCS)` |
| `iisort` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichisort.c:1014` | external | tree-sitter+gcc-aux | `extern int *iisort (int *list, int num)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c`

Parse errors: `22`. Function definitions: `47`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `comp_AT_NUMB` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:140` | external | tree-sitter+gcc-aux | `extern int comp_AT_NUMB (const void *a1, const void *a2, void *p)` |
| `get_z_coord` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:147` | internal | tree-sitter+gcc-aux | `static double get_z_coord (inp_ATOM *at, int cur_atom, int neigh_no, int *nType, int bPointedEdgeStereo)` |
| `len3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:221` | internal | tree-sitter+gcc-aux | `static double len3 (const double *c)` |
| `len2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:234` | internal | tree-sitter+gcc-aux | `static double len2 (const double *c)` |
| `diff3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:247` | internal | tree-sitter+gcc-aux | `static void *diff3 (const double *a, const double *b, double *result)` |
| `add3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:259` | internal | tree-sitter+gcc-aux | `static void add3 (const double *a, const double *b, double *result)` |
| `mult3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:270` | internal | tree-sitter+gcc-aux | `static void mult3 (const double *a, double b, double *result)` |
| `change_sign3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:292` | internal | tree-sitter+gcc-aux | `static void change_sign3 (const double *a, double *result)` |
| `dot_prod3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:303` | external | tree-sitter+gcc-aux | `extern double dot_prod3 (const double *a, const double *b)` |
| `dot_prodchar3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:310` | internal | tree-sitter+gcc-aux | `static int dot_prodchar3 (const S_CHAR *a, const S_CHAR *b)` |
| `cross_prod3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:330` | external | tree-sitter+gcc-aux | `extern void *cross_prod3 (const double *a, const double *b, double *result)` |
| `triple_prod` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:347` | internal | tree-sitter+gcc-aux | `static double triple_prod (double *a, double *b, double *c, double *sine_value)` |
| `CompDble` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:379` | internal | tree-sitter+gcc-aux | `static int CompDble (const void *a1, const void *a2, void *p)` |
| `Get2DTetrahedralAmbiguity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:402` | internal | tree-sitter+gcc-aux | `static int Get2DTetrahedralAmbiguity (CANON_GLOBALS *pCG, double (*at_coord)[3], int bAddExplicitNeighbor, int bFix2DstereoBorderCase, double vMinAngle)` |
| `triple_prod_and_min_abs_sine2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:926` | internal | tree-sitter+gcc-aux | `static double triple_prod_and_min_abs_sine2 (double (*at_coord)[3], double *central_at_coord, int bAddedExplicitNeighbor, double *min_sine, int *bAmbiguous, double vMinSine)` |
| `triple_prod_and_min_abs_sine` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1145` | internal | tree-sitter+gcc-aux | `static double triple_prod_and_min_abs_sine (double (*at_coord)[3], double *min_sine)` |
| `are_3_vect_in_one_plane` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1176` | internal | tree-sitter+gcc-aux | `static int are_3_vect_in_one_plane (double (*at_coord)[3], double min_sine)` |
| `are_4at_in_one_plane` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1190` | internal | tree-sitter+gcc-aux | `static int are_4at_in_one_plane (double (*at_coord)[3], double min_sine)` |
| `triple_prod_char` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1218` | internal | tree-sitter+gcc-aux | `static int triple_prod_char (inp_ATOM *at, int at_1, int i_next_at_1, S_CHAR *z_dir1, int at_2, int i_next_at_2, S_CHAR *z_dir2)` |
| `bInpAtomHasRequirdNeigh` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1317` | internal | tree-sitter+gcc-aux | `static int bInpAtomHasRequirdNeigh (inp_ATOM *at, int cur_at, int RequirdNeighType, int NumDbleBonds, int bStereoAtZz)` |
| `bCanInpAtomBeAStereoCenter` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1432` | external | tree-sitter+gcc-aux | `extern int bCanInpAtomBeAStereoCenter (inp_ATOM *at, int cur_at, int bPointedEdgeStereo, int bStereoAtZz)` |
| `bCanAtomBeAStereoCenter` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1553` | external | tree-sitter | `int bCanAtomBeAStereoCenter( char *elname, S_CHAR charge, S_CHAR radical )` |
| `bAtomHasValence3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1574` | external | tree-sitter+gcc-aux | `extern int bAtomHasValence3 (char *elname, S_CHAR charge, S_CHAR radical)` |
| `bCanAtomHaveAStereoBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1594` | external | tree-sitter+gcc-aux | `extern int bCanAtomHaveAStereoBond (char *elname, S_CHAR charge, S_CHAR radical)` |
| `bCanAtomBeMiddleAllene` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1615` | external | tree-sitter+gcc-aux | `extern int bCanAtomBeMiddleAllene (char *elname, S_CHAR charge, S_CHAR radical)` |
| `bIsSuitableHeteroInpAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1634` | internal | tree-sitter+gcc-aux | `static int bIsSuitableHeteroInpAtom (inp_ATOM *at)` |
| `bIsOxide` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1667` | internal | tree-sitter+gcc-aux | `static int bIsOxide (inp_ATOM *at, int cur_at)` |
| `bCanAtomBeTerminalAllene` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1705` | external | tree-sitter+gcc-aux | `extern int bCanAtomBeTerminalAllene (char *elname, S_CHAR charge, S_CHAR radical)` |
| `GetHalfStereobond0DParity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1725` | external | tree-sitter+gcc-aux | `extern int GetHalfStereobond0DParity (inp_ATOM *at, int cur_at, AT_NUMB *nSbNeighOrigAtNumb, int nNumExplictAttachments, int bond_parity, int nFlag)` |
| `FixSb0DParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:1824` | external | tree-sitter+gcc-aux | `extern int FixSb0DParities (inp_ATOM *at, int chain_length, int at_1, int i_next_at_1, S_CHAR *z_dir1, int at_2, int i_next_at_2, S_CHAR *z_dir2, int *pparity1, int *pparity2)` |
| `FixUnkn0DStereoBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:2006` | external | tree-sitter+gcc-aux | `extern int FixUnkn0DStereoBonds (inp_ATOM *at, int num_at)` |
| `half_stereo_bond_parity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:2121` | internal | tree-sitter+gcc-aux | `static int half_stereo_bond_parity (inp_ATOM *at, int cur_at, inp_ATOM *at_removed_H, int num_removed_H, S_CHAR *z_dir, int bPointedEdgeStereo, int vABParityUnknown)` |
| `save_a_stereo_bond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:2542` | internal | tree-sitter+gcc-aux | `static int save_a_stereo_bond (int z_prod, int result_action, int at1, int ord1, AT_NUMB *stereo_bond_neighbor1, S_CHAR *stereo_bond_ord1, S_CHAR *stereo_bond_z_prod1, S_CHAR *stereo_bond_parity1, int at2, int ord2, AT_NUMB *stereo_bond_neighbor2, S_CHAR *stereo_bond_ord2, S_CHAR *stereo_bond_z_prod2, S_CHAR *stereo_bond_parity2)` |
| `get_allowed_stereo_bond_type` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:2588` | internal | tree-sitter+gcc-aux | `static int get_allowed_stereo_bond_type (int bond_type)` |
| `can_be_a_stereo_bond_with_isotopic_H` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:2638` | internal | tree-sitter+gcc-aux | `static int can_be_a_stereo_bond_with_isotopic_H (inp_ATOM *at, int cur_at, INCHI_MODE nMode)` |
| `half_stereo_bond_action` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:2901` | internal | tree-sitter+gcc-aux | `static int half_stereo_bond_action (int nParity, int bUnknown, int bIsotopic, int vABParityUnknown)` |
| `set_stereo_bonds_parity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:3009` | internal | tree-sitter+gcc-aux | `static int set_stereo_bonds_parity (sp_ATOM *out_at, inp_ATOM *at, int at_1, inp_ATOM *at_removed_H, int num_removed_H, INCHI_MODE nMode, QUEUE *q, AT_RANK *nAtomLevel, S_CHAR *cSource, AT_RANK min_sb_ring_size, int bPointedEdgeStereo, int vABParityUnknown)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:3261` | external | tree-sitter | `else if (at[at_2].c_point)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:3266` | external | tree-sitter | `else if (num_2s_2 > 2)` |
| `can_be_a_stereo_atom_with_isotopic_H` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:3690` | internal | tree-sitter+gcc-aux | `static int can_be_a_stereo_atom_with_isotopic_H (inp_ATOM *at, int cur_at, int bPointedEdgeStereo, int bStereoAtZz)` |
| `can_be_a_stereo_atom_with_isotopic_H` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:3711` | external | tree-sitter | `int can_be_a_stereo_atom_with_isotopic_H( inp_ATOM *at, int cur_at )` |
| `GetStereocenter0DParity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:3734` | external | tree-sitter+gcc-aux | `extern int GetStereocenter0DParity (CANON_GLOBALS *pCG, inp_ATOM *at, int cur_at, int j1, AT_NUMB *nSbNeighOrigAtNumb, int nFlag)` |
| `set_stereo_atom_parity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:3790` | internal | tree-sitter+gcc-aux | `static int set_stereo_atom_parity (CANON_GLOBALS *pCG, sp_ATOM *out_at, inp_ATOM *at, int cur_at, inp_ATOM *at_removed_H, int num_removed_H, int bPointedEdgeStereo, int vABParityUnknown, int LooseTSACheck, int bStereoAtZz)` |
| `set_stereo_parity` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:4398` | external | tree-sitter+gcc-aux | `extern int set_stereo_parity (CANON_GLOBALS *pCG, inp_ATOM *at, sp_ATOM *at_output, int num_at, int num_removed_H, int *nMaxNumStereoAtoms, int *nMaxNumStereoBonds, INCHI_MODE nMode, int bPointedEdgeStereo, int vABParityUnknown, int bLooseTSACheck, int bStereoAtZz)` |
| `ReconcileAllCmlBondParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:4663` | external | tree-sitter+gcc-aux | `extern int ReconcileAllCmlBondParities (inp_ATOM *at, int num_atoms, int bDisconnected)` |
| `ReconcileCmlIncidentBondParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:4690` | external | tree-sitter+gcc-aux | `extern int ReconcileCmlIncidentBondParities (inp_ATOM *at, int cur_atom, int prev_atom, S_CHAR *visited, int bDisconnected)` |
| `get_opposite_sb_atom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichister.c:4861` | external | tree-sitter+gcc-aux | `extern int get_opposite_sb_atom (inp_ATOM *at, int cur_atom, int icur2nxt, int *pnxt_atom, int *pinxt2cur, int *pinxt_sb_parity_ord)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c`

Parse errors: `29`. Function definitions: `41`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `is_centerpoint_elem` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:157` | external | tree-sitter+gcc-aux | `extern int is_centerpoint_elem (U_CHAR el_number)` |
| `is_centerpoint_elem_KET` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:182` | external | tree-sitter+gcc-aux | `extern int is_centerpoint_elem_KET (U_CHAR el_number)` |
| `is_centerpoint_elem_strict` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:190` | external | tree-sitter+gcc-aux | `extern int is_centerpoint_elem_strict (U_CHAR el_number)` |
| `AddAtom2num` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:211` | external | tree-sitter+gcc-aux | `extern int AddAtom2num (AT_RANK *num, inp_ATOM *atom, int at_no, int bSubtract)` |
| `AddAtom2DA` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:250` | external | tree-sitter+gcc-aux | `extern void AddAtom2DA (AT_RANK *num_DA, inp_ATOM *atom, int at_no, int bSubtract)` |
| `AddEndPoint` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:330` | external | tree-sitter+gcc-aux | `extern int AddEndPoint (T_ENDPOINT *pEndPoint, inp_ATOM *at, int iat)` |
| `nGetEndpointInfo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:359` | external | tree-sitter+gcc-aux | `extern int nGetEndpointInfo (inp_ATOM *atom, int iat, ENDPOINT_INFO *eif)` |
| `nGetEndpointInfo_PT_22_00` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:452` | external | tree-sitter+gcc-aux | `extern int nGetEndpointInfo_PT_22_00 (inp_ATOM *atom, int iat, ENDPOINT_INFO *eif)` |
| `nGetEndpointInfo_PT_16_00` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:524` | external | tree-sitter+gcc-aux | `extern int nGetEndpointInfo_PT_16_00 (inp_ATOM *atom, int iat, ENDPOINT_INFO *eif)` |
| `nGetEndpointInfo_PT_06_00` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:600` | external | tree-sitter+gcc-aux | `extern int nGetEndpointInfo_PT_06_00 (inp_ATOM *atom, int iat, ENDPOINT_INFO *eif)` |
| `nGetEndpointInfo_PT_39_00` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:677` | external | tree-sitter+gcc-aux | `extern int nGetEndpointInfo_PT_39_00 (inp_ATOM *atom, int iat, ENDPOINT_INFO *eif)` |
| `nGetEndpointInfo_PT_13_00` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:756` | external | tree-sitter+gcc-aux | `extern int nGetEndpointInfo_PT_13_00 (inp_ATOM *atom, int iat, ENDPOINT_INFO *eif)` |
| `nGetEndpointInfo_PT_18_00` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:832` | external | tree-sitter+gcc-aux | `extern int nGetEndpointInfo_PT_18_00 (inp_ATOM *atom, int iat, ENDPOINT_INFO *eif)` |
| `nGetEndpointInfo_KET` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:916` | external | tree-sitter+gcc-aux | `extern int nGetEndpointInfo_KET (inp_ATOM *atom, int iat, ENDPOINT_INFO *eif)` |
| `RegisterEndPoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:1021` | external | tree-sitter+gcc-aux | `extern int RegisterEndPoints (CANON_GLOBALS *pCG, T_GROUP_INFO *t_group_info, T_ENDPOINT *EndPoint, int nNumEndPoints, inp_ATOM *at, int num_atoms, C_GROUP_INFO *cgi, struct BalancedNetworkStructure *pBNS)` |
| `SetTautomericBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:1523` | external | tree-sitter+gcc-aux | `extern int SetTautomericBonds (inp_ATOM *at, int nNumBondPos, T_BONDPOS *BondPos)` |
| `GetNeutralRepsIfNeeded` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:1565` | external | tree-sitter+gcc-aux | `extern int GetNeutralRepsIfNeeded (AT_NUMB *pri, AT_NUMB *prj, inp_ATOM *at, int num_atoms, T_ENDPOINT *EndPoint, int nNumEndPoints, C_GROUP_INFO *cgi)` |
| `FindAccessibleEndPoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:1712` | external | tree-sitter+gcc-aux | `extern int FindAccessibleEndPoints (CANON_GLOBALS *pCG, T_ENDPOINT *EndPoint, int *nNumEndPoints, T_BONDPOS *BondPos, int *nNumBondPos, struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD, inp_ATOM *at, int num_atoms, C_GROUP_INFO *cgi, int taut_mode)` |
| `bCanBeACPoint` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2047` | external | tree-sitter+gcc-aux | `extern int bCanBeACPoint (inp_ATOM *at, S_CHAR cCharge, S_CHAR cChangeValence, S_CHAR neutral_bonds_valence, S_CHAR neutral_valence, S_CHAR nEndpointValence, S_CHAR *cChargeSubtype)` |
| `GetChargeType` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2181` | external | tree-sitter+gcc-aux | `extern int GetChargeType (inp_ATOM *atom, int iat, S_CHAR *cChargeSubtype)` |
| `CmpCCandidates` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2229` | external | tree-sitter+gcc-aux | `extern int CmpCCandidates (const void *a1, const void *a2)` |
| `RegisterCPoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2249` | external | tree-sitter+gcc-aux | `extern int RegisterCPoints (C_GROUP *c_group, int *pnum_c, int max_num_c, T_GROUP_INFO *t_group_info, int point1, int point2, int ctype, inp_ATOM *at, int num_atoms)` |
| `MarkChargeGroups` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2397` | external | tree-sitter+gcc-aux | `extern int MarkChargeGroups (struct tagCANON_GLOBALS *pCG, inp_ATOM *at, int num_atoms, C_GROUP_INFO *c_group_info, T_GROUP_INFO *t_group_info, struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD)` |
| `GetSaltChargeType` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2565` | external | tree-sitter+gcc-aux | `extern int GetSaltChargeType (inp_ATOM *at, int at_no, T_GROUP_INFO *t_group_info, int *s_subtype)` |
| `bDoNotMergeNonTautAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2691` | external | tree-sitter+gcc-aux | `extern int bDoNotMergeNonTautAtom (inp_ATOM *at, int at_no)` |
| `GetOtherSaltChargeType` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2702` | external | tree-sitter+gcc-aux | `extern int GetOtherSaltChargeType (inp_ATOM *at, int at_no, T_GROUP_INFO *t_group_info, int *s_subtype, int bAccept_O)` |
| `GetOtherSaltType` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2828` | external | tree-sitter+gcc-aux | `extern int GetOtherSaltType (inp_ATOM *at, int at_no, int *s_subtype)` |
| `comp_candidates` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2930` | external | tree-sitter+gcc-aux | `extern int comp_candidates (const void *a1, const void *a2)` |
| `MarkSaltChargeGroups2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:2961` | external | tree-sitter+gcc-aux | `extern int MarkSaltChargeGroups2 (CANON_GLOBALS *pCG, inp_ATOM *at, int num_atoms, S_GROUP_INFO *s_group_info, T_GROUP_INFO *t_group_info, C_GROUP_INFO *c_group_info, struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD)` |
| `MarkSaltChargeGroups` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:3483` | external | tree-sitter+gcc-aux | `extern int MarkSaltChargeGroups (CANON_GLOBALS *pCG, inp_ATOM *at, int num_atoms, S_GROUP_INFO *s_group_info, T_GROUP_INFO *t_group_info, C_GROUP_INFO *c_group_info, struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD)` |
| `MarkSaltChargeGroups` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:3753` | external | tree-sitter | `int MarkSaltChargeGroups( CANON_GLOBALS *pCG, inp_ATOM *at, int num_atoms, S_GROUP_INFO *s_group_info, T_GROUP_INFO *t_group_info, C_GROUP_INFO *c_group_info, struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD )` |
| `MergeSaltTautGroups` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:3953` | external | tree-sitter+gcc-aux | `extern int MergeSaltTautGroups (CANON_GLOBALS *pCG, inp_ATOM *at, int num_atoms, S_GROUP_INFO *s_group_info, T_GROUP_INFO *t_group_info, C_GROUP_INFO *c_group_info, struct BalancedNetworkStructure *pBNS)` |
| `MakeIsotopicHGroup` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:4156` | external | tree-sitter+gcc-aux | `extern int MakeIsotopicHGroup (inp_ATOM *at, int num_atoms, S_GROUP_INFO *s_group_info, T_GROUP_INFO *t_group_info)` |
| `MarkTautomerGroups` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:4336` | external | gcc-aux | `extern int MarkTautomerGroups (CANON_GLOBALS *pCG, inp_ATOM *at, int num_atoms, T_GROUP_INFO *t_group_info, C_GROUP_INFO *c_group_info, struct BalancedNetworkStructure *pBNS, struct BalancedNetworkData *pBD)` |
| `free_t_group_info` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:6336` | external | tree-sitter+gcc-aux | `extern int free_t_group_info (T_GROUP_INFO *t_group_info)` |
| `make_a_copy_of_t_group_info` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:6364` | external | tree-sitter+gcc-aux | `extern int make_a_copy_of_t_group_info (T_GROUP_INFO *t_group_info, T_GROUP_INFO *t_group_info_orig)` |
| `set_tautomer_iso_sort_keys` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:6477` | external | tree-sitter+gcc-aux | `extern int set_tautomer_iso_sort_keys (T_GROUP_INFO *t_group_info)` |
| `CountTautomerGroups` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:6519` | external | tree-sitter+gcc-aux | `extern int CountTautomerGroups (sp_ATOM *at, int num_atoms, T_GROUP_INFO *t_group_info)` |
| `CountTautomerGroupsInpAt` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:6772` | external | tree-sitter | `int CountTautomerGroupsInpAt( inp_ATOM *at, int num_atoms, T_GROUP_INFO *t_group_info )` |
| `CompRankTautomer` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:7008` | external | tree-sitter+gcc-aux | `extern int CompRankTautomer (const void *a1, const void *a2, void *p)` |
| `SortTautomerGroupsAndEndpoints` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ichitaut.c:7019` | external | tree-sitter+gcc-aux | `extern int SortTautomerGroupsAndEndpoints (CANON_GLOBALS *pCG, T_GROUP_INFO *t_group_info, int num_atoms, int num_at_tg, AT_RANK *nRank)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c`

Parse errors: `0`. Function definitions: `8`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `base26_triplet_1` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1173` | external | tree-sitter+gcc-aux | `extern const char *base26_triplet_1 (const unsigned char *a)` |
| `base26_triplet_2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1191` | external | tree-sitter+gcc-aux | `extern const char *base26_triplet_2 (const unsigned char *a)` |
| `base26_triplet_3` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1212` | external | tree-sitter+gcc-aux | `extern const char *base26_triplet_3 (const unsigned char *a)` |
| `base26_triplet_4` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1233` | external | tree-sitter+gcc-aux | `extern const char *base26_triplet_4 (const unsigned char *a)` |
| `base26_dublet_for_bits_28_to_36` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1262` | external | tree-sitter+gcc-aux | `extern const char *base26_dublet_for_bits_28_to_36 (unsigned char *a)` |
| `base26_dublet_for_bits_56_to_64` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1287` | external | tree-sitter+gcc-aux | `extern const char *base26_dublet_for_bits_56_to_64 (unsigned char *a)` |
| `get_xtra_hash_major_hex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1307` | external | tree-sitter+gcc-aux | `extern void get_xtra_hash_major_hex (const unsigned char *a, char *szXtra)` |
| `get_xtra_hash_minor_hex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_base26.c:1328` | external | tree-sitter+gcc-aux | `extern void get_xtra_hash_minor_hex (const unsigned char *a, char *szXtra)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_dll.c`

Parse errors: `8`. Function definitions: `8`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `GetStdINCHIKeyFromStdINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_dll.c:92` | external | tree-sitter+gcc-aux | `extern int GetStdINCHIKeyFromStdINCHI (const char *szINCHISource, char *szINCHIKey)` |
| `GetINCHIKeyFromINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_dll.c:113` | external | tree-sitter+gcc-aux | `extern int GetINCHIKeyFromINCHI (const char *szINCHISource, const const int xtra1, const const int xtra2, char *szINCHIKey, char *szXtra1, char *szXtra2)` |
| `CheckINCHIKey` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_dll.c:513` | external | tree-sitter+gcc-aux | `extern int CheckINCHIKey (const char *szINCHIKey)` |
| `fprint_digest` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_dll.c:601` | external | tree-sitter+gcc-aux | `extern void fprint_digest (FILE *fw, const char *header, unsigned char *a)` |
| `cdecl_GetStdINCHIKeyFromStdINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_dll.c:643` | external | tree-sitter | `int cdecl_GetStdINCHIKeyFromStdINCHI( const char* szINCHISource, char* szINCHIKey )` |
| `cdecl_CheckINCHIKey` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_dll.c:650` | external | tree-sitter | `int cdecl_CheckINCHIKey( const char *szINCHIKey )` |
| `cdecl_GetINCHIKeyFromINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_dll.c:657` | external | tree-sitter | `int cdecl_GetINCHIKeyFromINCHI( const char* szINCHISource, const int xtra1, const int xtra2, char* szINCHIKey, char* szXtra1, char* szXtra2 )` |
| `pasc_CheckINCHIKey` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/ikey_dll.c:691` | external | tree-sitter | `int PASCAL pasc_CheckINCHIKey( const char *szINCHIKey )` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/inchi_gui.c`

Parse errors: `1`. Function definitions: `3`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `DisplayStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/inchi_gui.c:61` | external | tree-sitter | `int DisplayStructure( struct tagCANON_GLOBALS *pCG, inp_ATOM *at, int num_at, OAD_Polymer *polymer, int num_removed_H, int bAdd_DT_to_num_H, int nNumRemovedProtons, NUM_H *nNumRemovedProtonsIsotopic, int bIsotopic, int j /*bTautomeric*/, INChI **cur_INChI, INChI_Aux **cur_INChI_Aux, int bAbcNumbers, DRAW_PARMS *dp, INCHI_MODE nMode, char *szTitle )` |
| `DisplayCompositeStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/inchi_gui.c:104` | external | tree-sitter | `int DisplayCompositeStructure( struct tagCANON_GLOBALS *pCG, COMP_ATOM_DATA *composite_norm_data, OAD_Polymer *polymer, int bIsotopic, int bTautomeric, PINChI2 *pINChI2, PINChI_Aux2 *pINChI_Aux2, int bAbcNumbers, DRAW_PARMS *dp, INCHI_MODE nMode, char *szTitle )` |
| `FillCompositeTableParms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/inchi_gui.c:319` | external | tree-sitter | `void FillCompositeTableParms( SET_DRAW_PARMS *sdp, AT_NUMB StereoFlags, INCHI_MODE nMode, int bShowIsotopic, int bShowTaut )` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c`

Parse errors: `1`. Function definitions: `20`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `CreateOrigInpDataFromMolfile` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:108` | external | tree-sitter+gcc-aux | `extern int CreateOrigInpDataFromMolfile (INCHI_IOSTREAM *inp_file, ORIG_ATOM_DATA *orig_at_data, int bMergeAllInputStructures, int bGetOrigCoord, int bDoNotAddH, int treat_polymers, int treat_NPZz, const char *pSdfLabel, char *pSdfValue, long unsigned int *lSdfId, long int *lMolfileNumber, INCHI_MODE *pInpAtomFlags, int *err, char *pStrErr, int bNoWarnings)` |
| `ReadMolfileToInpAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:309` | external | tree-sitter+gcc-aux | `extern int ReadMolfileToInpAtoms (INCHI_IOSTREAM *inp_file, int bDoNotAddH, inp_ATOM **at, MOL_COORD (**szCoord), OAD_Polymer **polymer, OAD_V3000 **v3000, int treat_polymers, int treat_NPZz, int max_num_at, int *num_dimensions, int *num_bonds, const char *pSdfLabel, char *pSdfValue, long unsigned int *Id, long int *lMolfileNumber, INCHI_MODE *pInpAtomFlags, int *err, char *pStrErr, int bNoWarnings)` |
| `MakeInpAtomsFromMolfileData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:508` | external | tree-sitter+gcc-aux | `extern inp_ATOM *MakeInpAtomsFromMolfileData (MOL_FMT_DATA *mfdata, int *num_atoms, int *num_bonds, inp_ATOM *at_inp, int bDoNotAddH, int *err, char *pStrErr)` |
| `calculate_valences` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:807` | external | tree-sitter+gcc-aux | `extern void calculate_valences (MOL_FMT_DATA *mfdata, inp_ATOM *at, int *num_atoms, int bDoNotAddH, int *err, char *pStrErr)` |
| `SetInpAtomsXYZ` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:975` | external | tree-sitter+gcc-aux | `extern int SetInpAtomsXYZ (MOL_FMT_DATA *mfdata, int num_atoms, inp_ATOM *at, int *err, char *pStrErr)` |
| `CreateInpAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1033` | external | tree-sitter+gcc-aux | `extern inp_ATOM *CreateInpAtom (int num_atoms)` |
| `FreeInpAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1046` | external | tree-sitter+gcc-aux | `extern void FreeInpAtom (inp_ATOM **at)` |
| `FreeInpAtomData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1058` | external | tree-sitter+gcc-aux | `extern void FreeInpAtomData (INP_ATOM_DATA *inp_at_data)` |
| `CreateInpAtomData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1071` | external | tree-sitter+gcc-aux | `extern int CreateInpAtomData (INP_ATOM_DATA *inp_at_data, int num_atoms, int create_at_fixed_bonds)` |
| `FreeCompAtomData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1090` | external | tree-sitter+gcc-aux | `extern void FreeCompAtomData (COMP_ATOM_DATA *inp_at_data)` |
| `CreateCompAtomData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1104` | external | tree-sitter | `int CreateCompAtomData(COMP_ATOM_DATA *inp_at_data, int num_atoms, int num_components, int bIntermediateTaut)` |
| `FreeInfAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1129` | external | tree-sitter | `void FreeInfAtom(inf_ATOM **at)` |
| `CreateInfAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1140` | external | tree-sitter | `inf_ATOM *CreateInfAtom(int num_atoms)` |
| `FreeInfoAtomData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1146` | external | tree-sitter | `void FreeInfoAtomData(INF_ATOM_DATA *inf_at_data)` |
| `CreateInfoAtomData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1159` | external | tree-sitter | `int CreateInfoAtomData(INF_ATOM_DATA *inf_at_data, int num_atoms, int num_components)` |
| `AllocateInfoAtomData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1183` | external | tree-sitter | `int AllocateInfoAtomData(INF_ATOM_DATA *inf_at_data, int num_atoms, int num_components)` |
| `DuplicateInfoAtomData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1202` | external | tree-sitter | `int DuplicateInfoAtomData(INF_ATOM_DATA *inf_at_data_to, const INF_ATOM_DATA *inf_at_data_from)` |
| `FreeOrigAtData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1224` | external | tree-sitter+gcc-aux | `extern void FreeOrigAtData (ORIG_ATOM_DATA *orig_at_data)` |
| `FreeExtOrigAtData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1269` | external | tree-sitter+gcc-aux | `extern void FreeExtOrigAtData (OAD_Polymer *pd, OAD_V3000 *v3k)` |
| `SetExtOrigAtDataByMolfileExtInput` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol2atom.c:1340` | external | tree-sitter+gcc-aux | `extern int SetExtOrigAtDataByMolfileExtInput (MOL_FMT_DATA *mfdata, OAD_Polymer **ppPolymer, OAD_V3000 **ppV3000, char *pStrErr)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c`

Parse errors: `21`. Function definitions: `10`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `ReadMolfile` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:92` | external | tree-sitter+gcc-aux | `extern MOL_FMT_DATA *ReadMolfile (INCHI_IOSTREAM *inp_file, MOL_FMT_HEADER_BLOCK *OnlyHeaderBlock, MOL_FMT_CTAB *OnlyCTab, int bGetOrigCoord, int treat_polymers, int treat_NPZz, char *pname, int lname, long unsigned int *Id, const char *pSdfLabel, char *pSdfValue, int *err, char *pStrErr, int bNoWarnings)` |
| `MolfileReadDataLines` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:165` | internal | tree-sitter+gcc-aux | `static MOL_FMT_DATA *MolfileReadDataLines (INCHI_IOSTREAM *inp_file, MOL_FMT_HEADER_BLOCK *OnlyHeaderBlock, MOL_FMT_CTAB *OnlyCTab, int bGetOrigCoord, int treat_polymers, int *err, char *pStrErr, int bNoWarnings)` |
| `MolfileReadHeaderLines` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:451` | internal | tree-sitter+gcc-aux | `static int MolfileReadHeaderLines (MOL_FMT_HEADER_BLOCK *hdr, INCHI_IOSTREAM *inp_file, char *pStrErr)` |
| `MolfileReadCountsLine` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:569` | internal | tree-sitter+gcc-aux | `static int MolfileReadCountsLine (MOL_FMT_CTAB *ctab, INCHI_IOSTREAM *inp_file, char *pStrErr)` |
| `MolfileReadAtomsBlock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:661` | internal | tree-sitter+gcc-aux | `static int MolfileReadAtomsBlock (MOL_FMT_CTAB *ctab, INCHI_IOSTREAM *inp_file, int err, char *pStrErr)` |
| `MolfileReadBondsBlock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:820` | internal | tree-sitter+gcc-aux | `static int MolfileReadBondsBlock (MOL_FMT_CTAB *ctab, INCHI_IOSTREAM *inp_file, int err, char *pStrErr)` |
| `MolfileReadSTextBlock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:917` | internal | tree-sitter+gcc-aux | `static int MolfileReadSTextBlock (MOL_FMT_CTAB *ctab, INCHI_IOSTREAM *inp_file, int err, char *pStrErr)` |
| `MolfileReadPropBlock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:960` | internal | tree-sitter+gcc-aux | `static int MolfileReadPropBlock (MOL_FMT_CTAB *ctab, MOL_FMT_HEADER_BLOCK *pHdr, INCHI_IOSTREAM *inp_file, int treat_polymers, int err, char *pStrErr, int bNoWarnings)` |
| `MolfileReadSgroupOfPolymer` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:1425` | internal | tree-sitter+gcc-aux | `static int MolfileReadSgroupOfPolymer (MOL_FMT_CTAB *ctab, MOL_FMT_HEADER_BLOCK *pHdr, INCHI_IOSTREAM *inp_file, char *line, char *szType, char *p, int err, char *pStrErr)` |
| `MolfileTreatPseudoElementAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt1.c:1860` | internal | tree-sitter+gcc-aux | `static int MolfileTreatPseudoElementAtoms (MOL_FMT_CTAB *ctab, int pseudos_allowed, int *err, char *pStrErr)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c`

Parse errors: `2`. Function definitions: `7`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `MolfileStrnread` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:65` | external | tree-sitter+gcc-aux | `extern int MolfileStrnread (char *dest, char *source, int len, char **first_space)` |
| `MolfileReadField` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:105` | external | tree-sitter+gcc-aux | `extern int MolfileReadField (void *data, int field_len, int data_type, char **line_ptr)` |
| `MolfileExtractStrucNum` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:333` | external | tree-sitter+gcc-aux | `extern long int MolfileExtractStrucNum (MOL_FMT_HEADER_BLOCK *pHdr)` |
| `MolfileHasNoChemStruc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:362` | external | tree-sitter+gcc-aux | `extern int MolfileHasNoChemStruc (MOL_FMT_DATA *mfdata)` |
| `MolfileSaveCopy` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:385` | external | tree-sitter+gcc-aux | `extern int MolfileSaveCopy (INCHI_IOSTREAM *inp_file, long int fPtrStart, long int fPtrEnd, FILE *outfile, long int num)` |
| `MolfileGetXYZDimAndNormFactors` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:474` | external | tree-sitter+gcc-aux | `extern int MolfileGetXYZDimAndNormFactors (MOL_FMT_DATA *mfdata, int find_norm_factors, double *x0, double *y0, double *z0, double *xmin, double *ymin, double *zmin, double *scaler, int *err, char *pStrErr)` |
| `FreeMolfileData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt2.c:664` | external | tree-sitter+gcc-aux | `extern MOL_FMT_DATA *FreeMolfileData (MOL_FMT_DATA *mfdata)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c`

Parse errors: `2`. Function definitions: `16`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `MolfileV3000Init` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:68` | external | tree-sitter+gcc-aux | `extern int MolfileV3000Init (MOL_FMT_CTAB *ctab, char *pStrErr)` |
| `DeleteMolfileV3000Info` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:165` | external | tree-sitter+gcc-aux | `extern int DeleteMolfileV3000Info (MOL_FMT_v3000 *v3000)` |
| `inchi_fgetsLf_V3000` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:216` | external | tree-sitter+gcc-aux | `extern char *inchi_fgetsLf_V3000 (char *line, INCHI_IOSTREAM *inp_stream)` |
| `MolfileV3000ReadField` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:257` | external | tree-sitter+gcc-aux | `extern int MolfileV3000ReadField (void *data, int data_type, char **line_ptr)` |
| `MolfileV3000ReadKeyword` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:410` | external | tree-sitter+gcc-aux | `extern int MolfileV3000ReadKeyword (char *key, char **line_ptr)` |
| `MolfileV3000ReadCTABBeginAndCountsLine` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:444` | external | tree-sitter+gcc-aux | `extern int MolfileV3000ReadCTABBeginAndCountsLine (MOL_FMT_CTAB *ctab, INCHI_IOSTREAM *inp_file, char *pStrErr)` |
| `MolfileV3000ReadSGroup` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:561` | external | tree-sitter+gcc-aux | `extern int MolfileV3000ReadSGroup (MOL_FMT_CTAB *ctab, INCHI_IOSTREAM *inp_file, int err, char *pStrErr)` |
| `MolfileV3000Read3DBlock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:608` | external | tree-sitter+gcc-aux | `extern int MolfileV3000Read3DBlock (MOL_FMT_CTAB *ctab, INCHI_IOSTREAM *inp_file, int err, char *pStrErr)` |
| `MolfileV3000ReadCollections` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:647` | external | tree-sitter+gcc-aux | `extern int MolfileV3000ReadCollections (MOL_FMT_CTAB *ctab, INCHI_IOSTREAM *inp_file, int err, char *pStrErr)` |
| `MolfileV3000ReadAtomsBlock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:862` | external | tree-sitter+gcc-aux | `extern int MolfileV3000ReadAtomsBlock (MOL_FMT_CTAB *ctab, INCHI_IOSTREAM *inp_file, int err, char *pStrErr)` |
| `MolfileV3000ReadBondsBlock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:1262` | external | tree-sitter+gcc-aux | `extern int MolfileV3000ReadBondsBlock (MOL_FMT_CTAB *ctab, INCHI_IOSTREAM *inp_file, int err, char *pStrErr)` |
| `get_actual_atom_number` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:1585` | internal | tree-sitter+gcc-aux | `static int get_actual_atom_number (int index, int n, int *orig, int *fin)` |
| `MolfileV3000ReadTailOfCTAB` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:1602` | external | tree-sitter+gcc-aux | `extern int MolfileV3000ReadTailOfCTAB (MOL_FMT_CTAB *ctab, INCHI_IOSTREAM *inp_file, int err, char *pStrErr)` |
| `MolfileV3000ReadHapticBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:1732` | external | tree-sitter+gcc-aux | `extern int MolfileV3000ReadHapticBond (MOL_FMT_CTAB *ctab, char **line_ptr, int **num_list, char *pStrErr)` |
| `MolfileV3000ReadStereoCollection` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:1817` | external | tree-sitter+gcc-aux | `extern int MolfileV3000ReadStereoCollection (MOL_FMT_CTAB *ctab, char **line_ptr, int **num_list, char *pStrErr)` |
| `get_V3000_input_line_to_strbuf` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt3.c:1889` | external | tree-sitter+gcc-aux | `extern int get_V3000_input_line_to_strbuf (INCHI_IOS_STRING *buf, INCHI_IOSTREAM *inp_stream)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c`

Parse errors: `3`. Function definitions: `27`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `SDFileSkipExtraData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:161` | external | tree-sitter+gcc-aux | `extern int SDFileSkipExtraData (INCHI_IOSTREAM *inp_file, long unsigned int *CAS_num, char *comment, int lcomment, char *name, int lname, int prev_err, const char *pSdfLabel, char *pSdfValue, char *pStrErr, int bNoWarnings)` |
| `SDFileIdentifyLabel` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:346` | external | tree-sitter+gcc-aux | `extern int SDFileIdentifyLabel (char *inp_line, const char *pSdfLabel)` |
| `SDFileExtractCASNo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:407` | external | tree-sitter+gcc-aux | `extern long unsigned int SDFileExtractCASNo (char *line)` |
| `NumLists_Alloc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:433` | external | tree-sitter+gcc-aux | `extern int NumLists_Alloc (NUM_LISTS *num_lists, int nlists)` |
| `NumLists_ReAlloc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:449` | external | tree-sitter+gcc-aux | `extern int NumLists_ReAlloc (NUM_LISTS *num_lists)` |
| `NumLists_Append` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:471` | external | tree-sitter+gcc-aux | `extern int NumLists_Append (NUM_LISTS *num_lists, int *list)` |
| `NumLists_Free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:491` | external | tree-sitter+gcc-aux | `extern void NumLists_Free (NUM_LISTS *num_lists)` |
| `IntArray_Alloc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:510` | external | tree-sitter+gcc-aux | `extern int IntArray_Alloc (INT_ARRAY *items, int nitems)` |
| `IntArray_ReAlloc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:525` | external | tree-sitter+gcc-aux | `extern int IntArray_ReAlloc (INT_ARRAY *items)` |
| `IntArray_Append` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:550` | external | tree-sitter+gcc-aux | `extern int IntArray_Append (INT_ARRAY *items, int new_item)` |
| `IntArray_AppendIfAbsent` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:572` | external | tree-sitter+gcc-aux | `extern int IntArray_AppendIfAbsent (INT_ARRAY *items, int new_item)` |
| `IntArray_DebugPrint` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:582` | external | tree-sitter+gcc-aux | `extern void IntArray_DebugPrint (INT_ARRAY *items)` |
| `IntArray_Reset` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:603` | external | tree-sitter+gcc-aux | `extern void IntArray_Reset (INT_ARRAY *items)` |
| `IntArray_Free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:610` | external | tree-sitter+gcc-aux | `extern void IntArray_Free (INT_ARRAY *items)` |
| `MolFmtSgroup_Create` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:634` | external | tree-sitter+gcc-aux | `extern int MolFmtSgroup_Create (MOL_FMT_SGROUP **sgroup, int id, int type)` |
| `MolFmtSgroup_Free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:658` | external | tree-sitter+gcc-aux | `extern void MolFmtSgroup_Free (MOL_FMT_SGROUP *sgroup)` |
| `MolFmtSgroups_Alloc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:675` | external | tree-sitter+gcc-aux | `extern int MolFmtSgroups_Alloc (MOL_FMT_SGROUPS *sgroups, int nsgroups)` |
| `MolFmtSgroups_ReAlloc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:693` | external | tree-sitter+gcc-aux | `extern int MolFmtSgroups_ReAlloc (MOL_FMT_SGROUPS *sgroups)` |
| `MolFmtSgroups_Append` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:715` | external | tree-sitter+gcc-aux | `extern int MolFmtSgroups_Append (MOL_FMT_SGROUPS *sgroups, int id, int type)` |
| `MolFmtSgroups_Free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:751` | external | tree-sitter+gcc-aux | `extern void MolFmtSgroups_Free (MOL_FMT_SGROUPS *sgroups)` |
| `MolFmtSgroups_GetIndexBySgroupId` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:769` | external | tree-sitter+gcc-aux | `extern int MolFmtSgroups_GetIndexBySgroupId (int id, MOL_FMT_SGROUPS *sgroups)` |
| `OrigAtData_WriteToSDfile` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:783` | external | tree-sitter+gcc-aux | `extern int OrigAtData_WriteToSDfile (const ORIG_ATOM_DATA *inp_at_data, INCHI_IOSTREAM *fcb, const char *name, const char *comment, int bChiralFlag, int bAtomsDT, const char *szLabel, const char *szValue)` |
| `OrigAtData_WriteToSDfileHeaderAndCountThings` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:852` | external | tree-sitter+gcc-aux | `extern int OrigAtData_WriteToSDfileHeaderAndCountThings (const ORIG_ATOM_DATA *inp_at_data, INCHI_IOSTREAM *fcb, const char *name, const char *comment, int bChiralFlag, int bAtomsDT, const char *szLabel, const char *szValue, int *nNumAliasLines, int *nNumChargeLines, int *nNumRadicalLines, int *nNumIsoLines, int *nNumAddLines, int *num_bonds)` |
| `OrigAtData_WriteToSDfileAtomsBlock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:989` | external | tree-sitter+gcc-aux | `extern int OrigAtData_WriteToSDfileAtomsBlock (const ORIG_ATOM_DATA *inp_at_data, INCHI_IOSTREAM *fcb, const char *name, const char *comment, int bAtomsDT, const char *szLabel, const char *szValue)` |
| `OrigAtData_WriteToSDfileBondsBlock` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:1122` | external | tree-sitter+gcc-aux | `extern int OrigAtData_WriteToSDfileBondsBlock (const ORIG_ATOM_DATA *inp_at_data, INCHI_IOSTREAM *fcb, const char *name, const char *comment, const char *szLabel, const char *szValue, INT_ARRAY *written_bond_ends)` |
| `OrigAtData_WriteToSDfileAdditionalLines` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:1182` | external | tree-sitter+gcc-aux | `extern int OrigAtData_WriteToSDfileAdditionalLines (const ORIG_ATOM_DATA *inp_at_data, INCHI_IOSTREAM *fcb, const char *name, const char *comment, int bAtomsDT, const char *szLabel, const char *szValue, int nNumAliasLines, int nNumChargeLines, int nNumRadicalLines, int nNumIsoLines, INT_ARRAY *written_bond_ends)` |
| `OrigAtData_WriteToSDfilePolymerData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/mol_fmt4.c:1390` | external | tree-sitter+gcc-aux | `extern int OrigAtData_WriteToSDfilePolymerData (const ORIG_ATOM_DATA *inp_at_data, INCHI_IOSTREAM *fcb, const char *name, const char *comment, const char *szLabel, const char *szValue, INT_ARRAY *written_bond_ends)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/permutation_util.c`

Parse errors: `0`. Function definitions: `3`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `rrand` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/permutation_util.c:74` | external | tree-sitter | `int rrand(int m)` |
| `shuffle` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/permutation_util.c:80` | external | tree-sitter | `void shuffle(void* obj, size_t nmemb, size_t size)` |
| `OrigAtData_Permute` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/permutation_util.c:103` | external | tree-sitter | `void OrigAtData_Permute(ORIG_ATOM_DATA* permuted, ORIG_ATOM_DATA* saved, int* numbers)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c`

Parse errors: `1`. Function definitions: `7`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `CreateInchi_Stereo0D` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:118` | external | tree-sitter+gcc-aux | `extern inchi_Stereo0D *CreateInchi_Stereo0D (int num_stereo0D)` |
| `FreeInchi_Stereo0D` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:125` | external | tree-sitter+gcc-aux | `extern void FreeInchi_Stereo0D (inchi_Stereo0D **stereo0D)` |
| `Extract0DParities` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:136` | external | tree-sitter+gcc-aux | `extern int Extract0DParities (inp_ATOM *at, int nNumAtoms, inchi_Stereo0D *stereo0D, int num_stereo0D, char *pStrErr, int *err, int vABParityUnknown)` |
| `FindToken` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:445` | external | tree-sitter+gcc-aux | `extern char *FindToken (INCHI_IOSTREAM *inp_file, int *bTooLongLine, const char *sToken, int lToken, char *szLine, int nLenLine, char *p, int *res)` |
| `LoadLine` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:489` | external | tree-sitter+gcc-aux | `extern char *LoadLine (INCHI_IOSTREAM *inp_file, int *bTooLongLine, int *bItemIsOver, char **s, char *szLine, int nLenLine, int nMinLen2Load, char *p, int *res)` |
| `InchiToInpAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:549` | external | tree-sitter+gcc-aux | `extern int InchiToInpAtom (INCHI_IOSTREAM *inp_file, MOL_COORD (**szCoord), int bDoNotAddH, int vABParityUnknown, INPUT_TYPE nInputType, inp_ATOM **at, int max_num_at, int *num_dimensions, int *num_bonds, char *pSdfLabel, char *pSdfValue, long unsigned int *Id, INCHI_MODE *pInpAtomFlags, int *err, char *pStrErr)` |
| `find_and_interpret_structure_header` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/readinch.c:1742` | external | tree-sitter+gcc-aux | `extern void find_and_interpret_structure_header (char *szLine, char *pSdfLabel, char *pSdfValue, long unsigned int *Id, int hlen, ReadINCHI_CtlData *ir)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c`

Parse errors: `1`. Function definitions: `21`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `ProcessOneStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:218` | external | tree-sitter+gcc-aux | `extern int ProcessOneStructure (INCHI_CLOCK *ic, CANON_GLOBALS *pCG, STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle, PINChI2 (**pINChI), PINChI_Aux2 (**pINChI_Aux), INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *out_file, INCHI_IOSTREAM *prb_file, ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data, long int num_inp, INCHI_IOS_STRING *strbuf, unsigned char save_opt_bits)` |
| `DoOneStructureEarlyPreprocessing` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:493` | internal | tree-sitter+gcc-aux | `static int DoOneStructureEarlyPreprocessing (INCHI_CLOCK *ic, CANON_GLOBALS *pCG, long int num_inp, STRUCT_DATA *sd, INPUT_PARMS *ip, INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *out_file, INCHI_IOSTREAM *prb_file, ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data)` |
| `OrigAtData_SaveMolfile` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:564` | external | tree-sitter+gcc-aux | `extern int OrigAtData_SaveMolfile (ORIG_ATOM_DATA *orig_inp_data, STRUCT_DATA *sd, INPUT_PARMS *ip, long int num_inp, INCHI_IOSTREAM *out_file)` |
| `OrigAtData_StoreNativeInput` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:593` | internal | tree-sitter+gcc-aux | `static ORIG_STRUCT *OrigAtData_StoreNativeInput (CANON_GLOBALS *pCG, int *nRet, STRUCT_DATA *sd, INPUT_PARMS *ip, ORIG_ATOM_DATA *orig_inp_data, ORIG_STRUCT *pOrigStruct)` |
| `PrepareSaveOptBits` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:621` | internal | tree-sitter+gcc-aux | `static void PrepareSaveOptBits (unsigned char *save_opt_bits, INPUT_PARMS *ip)` |
| `DisplayOrigAndResultStructuresAndComponents` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:679` | internal | tree-sitter+gcc-aux | `static void DisplayOrigAndResultStructuresAndComponents (int nRet, INCHI_CLOCK *ic, CANON_GLOBALS *pCG, STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle, PINChI2 (**pINChI), PINChI_Aux2 (**pINChI_Aux), INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *out_file, ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data, long int num_inp, int maxINChI, COMP_ATOM_DATA (*composite_norm_data)[3])` |
| `SaveOkProcessedMolfile` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:785` | internal | tree-sitter+gcc-aux | `static void SaveOkProcessedMolfile (int nRet, STRUCT_DATA *sd, INPUT_PARMS *ip, INCHI_IOSTREAM *prb_file, INCHI_IOSTREAM *inp_file)` |
| `CreateOneStructureINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:809` | external | tree-sitter+gcc-aux | `extern int CreateOneStructureINChI (CANON_GLOBALS *pCG, INCHI_CLOCK *ic, STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle, PINChI2 (**pINChI2), PINChI_Aux2 (**pINChI_Aux2), int iINChI, INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *out_file, INCHI_IOSTREAM *prb_file, ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data, COMP_ATOM_DATA (*composite_norm_data2)[3], long int num_inp, INCHI_IOS_STRING *strbuf, NORM_CANON_FLAGS *pncFlags)` |
| `CreateOneComponentINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:1747` | external | tree-sitter+gcc-aux | `extern int CreateOneComponentINChI (CANON_GLOBALS *pCG, INCHI_CLOCK *ic, STRUCT_DATA *sd, INPUT_PARMS *ip, INP_ATOM_DATA *inp_cur_data, ORIG_ATOM_DATA *orig_inp_data, PINChI2 (*pINChI), PINChI_Aux2 (*pINChI_Aux), int iINChI, int i, long int num_inp, INP_ATOM_DATA **inp_norm_data, NORM_CANON_FLAGS *pncFlags, INCHI_IOSTREAM *log_file)` |
| `ProcessOneStructureEx` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:2029` | external | tree-sitter+gcc-aux | `extern int ProcessOneStructureEx (struct tagINCHI_CLOCK *ic, struct tagCANON_GLOBALS *CG, STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle, PINChI2 (**pINChI2), PINChI_Aux2 (**pINChI_Aux2), INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *out_file, INCHI_IOSTREAM *prb_file, ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data, long int num_inp, INCHI_IOS_STRING *strbuf, unsigned char save_opt_bits)` |
| `PreprocessPolymerCRUData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:2147` | external | tree-sitter+gcc-aux | `extern int PreprocessPolymerCRUData (struct tagINCHI_CLOCK *ic, struct tagCANON_GLOBALS *CG, STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle, PINChI2 (**pINChI2), PINChI_Aux2 (**pINChI_Aux2), INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *out_file, INCHI_IOSTREAM *prb_file, ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data, long int num_inp, INCHI_IOS_STRING *strbuf, unsigned char save_opt_bits, char **sinchi_noedits, char **saux_noedits)` |
| `swap_atoms_xyz` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:2401` | internal | tree-sitter+gcc-aux | `static void swap_atoms_xyz (ORIG_ATOM_DATA *orig_at_data, int ia1, int ia2)` |
| `OAD_StructureEdits_Apply` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:2427` | external | tree-sitter+gcc-aux | `extern int OAD_StructureEdits_Apply (STRUCT_DATA *sd, INPUT_PARMS *ip, ORIG_ATOM_DATA *orig_at_data, OAD_StructureEdits *ed, int *ret)` |
| `set_renumbered_or_delete` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:2775` | internal | tree-sitter+gcc-aux | `static int set_renumbered_or_delete (int *number, int *buf, int nelems, int *renum, int base)` |
| `ProcessOneStructureExCore` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:2803` | external | tree-sitter+gcc-aux | `extern int ProcessOneStructureExCore (struct tagINCHI_CLOCK *ic, struct tagCANON_GLOBALS *CG, STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle, PINChI2 (**pINChI2), PINChI_Aux2 (**pINChI_Aux2), INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *out_file, INCHI_IOSTREAM *prb_file, ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data, long int num_inp, INCHI_IOS_STRING *strbuf, unsigned char save_opt_bits)` |
| `ValidateAndPreparePolymerAndPseudoatoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:2879` | external | tree-sitter+gcc-aux | `extern int ValidateAndPreparePolymerAndPseudoatoms (struct tagINCHI_CLOCK *ic, struct tagCANON_GLOBALS *CG, STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle, PINChI2 (**pINChI2), PINChI_Aux2 (**pINChI_Aux2), INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *out_file, INCHI_IOSTREAM *prb_file, ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data, long int num_inp, INCHI_IOS_STRING *strbuf, unsigned char save_opt_bits, int *mind_polymers)` |
| `OAD_ProcessOneStructureNoEdits` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:2967` | external | tree-sitter+gcc-aux | `extern int OAD_ProcessOneStructureNoEdits (struct tagINCHI_CLOCK *ic, struct tagCANON_GLOBALS *CG, STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle, PINChI2 (**pINChI2), PINChI_Aux2 (**pINChI_Aux2), INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *out_file, INCHI_IOSTREAM *prb_file, ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data, long int num_inp, INCHI_IOS_STRING *strbuf, unsigned char save_opt_bits, int *n_pzz, char **sinchi, char **saux)` |
| `OAD_ProcessOneStructure105Plus` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:3062` | external | tree-sitter+gcc-aux | `extern int OAD_ProcessOneStructure105Plus (struct tagINCHI_CLOCK *ic, struct tagCANON_GLOBALS *CG, STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle, PINChI2 (**pINChI2), PINChI_Aux2 (**pINChI_Aux2), INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *out_file, INCHI_IOSTREAM *prb_file, ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data, long int num_inp, INCHI_IOS_STRING *strbuf, unsigned char save_opt_bits, char **sinchi, char **saux)` |
| `mark_atoms_to_delete_or_renumber` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:3155` | internal | tree-sitter+gcc-aux | `static int mark_atoms_to_delete_or_renumber (ORIG_ATOM_DATA *orig_at_data, OAD_StructureEdits *ed, int *at_renum)` |
| `OrigAtData_CheckForSubstructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:3261` | external | tree-sitter | `int OrigAtData_CheckForSubstructure(ORIG_ATOM_DATA *orig_inp_data)` |
| `check_presence_of_the_encoded_substructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi.c:3269` | external | tree-sitter | `int check_presence_of_the_encoded_substructure(ORIG_ATOM_DATA *oad)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c`

Parse errors: `1`. Function definitions: `31`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `GetOneStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:134` | external | tree-sitter+gcc-aux | `extern int GetOneStructure (INCHI_CLOCK *ic, STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle, INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *out_file, INCHI_IOSTREAM *prb_file, ORIG_ATOM_DATA *orig_inp_data, long int *num_inp, STRUCT_FPTRS *struct_fptrs)` |
| `GetOneComponent` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:345` | external | tree-sitter+gcc-aux | `extern int GetOneComponent (INCHI_CLOCK *ic, STRUCT_DATA *sd, INPUT_PARMS *ip, INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *out_file, INP_ATOM_DATA *inp_cur_data, ORIG_ATOM_DATA *orig_inp_data, int i, long int num_inp)` |
| `ReadTheStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:396` | external | tree-sitter+gcc-aux | `extern int ReadTheStructure (struct tagINCHI_CLOCK *ic, STRUCT_DATA *sd, INPUT_PARMS *ip, INCHI_IOSTREAM *inp_file, ORIG_ATOM_DATA *orig_inp_data, int inp_index, int *out_index)` |
| `TreatErrorsInReadTheStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:716` | external | tree-sitter+gcc-aux | `extern int TreatErrorsInReadTheStructure (STRUCT_DATA *sd, INPUT_PARMS *ip, int nLogMask, INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *out_file, INCHI_IOSTREAM *prb_file, ORIG_ATOM_DATA *orig_inp_data, long int *num_inp)` |
| `InchiToInchi_Input` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:839` | external | tree-sitter | `int InchiToInchi_Input( INCHI_IOSTREAM *inp_file, inchi_Input *orig_at_data, int bMergeAllInputStructures, int bDoNotAddH, int vABParityUnknown, INPUT_TYPE nInputType, char *pSdfLabel, char *pSdfValue, unsigned long *lSdfId, INCHI_MODE *pInpAtomFlags, int *err, char *pStrErr )` |
| `InchiToOrigAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1067` | external | tree-sitter+gcc-aux | `extern int InchiToOrigAtom (INCHI_IOSTREAM *inp_molfile, ORIG_ATOM_DATA *orig_at_data, int bMergeAllInputStructures, int bGetOrigCoord, int bDoNotAddH, int vABParityUnknown, INPUT_TYPE nInputType, char *pSdfLabel, char *pSdfValue, long unsigned int *lSdfId, INCHI_MODE *pInpAtomFlags, int *err, char *pStrErr)` |
| `bIsSameBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1293` | external | tree-sitter+gcc-aux | `extern int bIsSameBond (int a1, int a2, int b1, int b2)` |
| `GetFrameShiftInfoFrom105PlusInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1305` | internal | tree-sitter+gcc-aux | `static int GetFrameShiftInfoFrom105PlusInChI (char *sinchi, int *frame_shift_info, int max_crossing)` |
| `extract_orig_nums_from_auxinfo_string` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1389` | external | tree-sitter+gcc-aux | `extern int extract_orig_nums_from_auxinfo_string (char *saux, int *orig)` |
| `extract_nonstereo_eq_classes_from_auxinfo_string` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1427` | external | tree-sitter+gcc-aux | `extern int extract_nonstereo_eq_classes_from_auxinfo_string (char *saux, int nat, int *orig, int *nclasses, int *eclass, int *eclass_by_origs)` |
| `POSEContext_Init` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1503` | external | tree-sitter+gcc-aux | `extern int POSEContext_Init (POSEContext *context, STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle, PINChI2 (**pINChI2), PINChI_Aux2 (**pINChI_Aux2), INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *out_file, INCHI_IOSTREAM *prb_file, ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data, long int num_inp, INCHI_IOS_STRING *strbuf, unsigned char save_opt_bits)` |
| `POSEContext_Free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1636` | external | tree-sitter+gcc-aux | `extern void POSEContext_Free (POSEContext *context)` |
| `POSEContext_DebugPrint` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1672` | external | tree-sitter+gcc-aux | `extern void POSEContext_DebugPrint (POSEContext *context)` |
| `OAD_StructureEdits_Init` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1735` | external | tree-sitter+gcc-aux | `extern int OAD_StructureEdits_Init (OAD_StructureEdits *ed)` |
| `OAD_StructureEdits_Clear` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1769` | external | tree-sitter+gcc-aux | `extern void OAD_StructureEdits_Clear (OAD_StructureEdits *ed)` |
| `OAD_StructureEdits_DebugPrint` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1807` | external | tree-sitter+gcc-aux | `extern void OAD_StructureEdits_DebugPrint (OAD_StructureEdits *ed)` |
| `OAD_Polymer_PrepareFoldCRUEdits` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:1838` | external | tree-sitter+gcc-aux | `extern int OAD_Polymer_PrepareFoldCRUEdits (ORIG_ATOM_DATA *orig_at_data, char *sinchi_noedits, char *saux_noedits, char *sinchi, char *saux, OAD_StructureEdits *ed)` |
| `DiylFrag_New` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2007` | internal | tree-sitter+gcc-aux | `static DiylFrag *DiylFrag_New (int na, int end1, int end2, char *s)` |
| `DiylFrag_Free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2049` | internal | tree-sitter+gcc-aux | `static void DiylFrag_Free (DiylFrag *pfrag)` |
| `DiylFrag_MakeSignature` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2069` | internal | tree-sitter+gcc-aux | `static void DiylFrag_MakeSignature (DiylFrag *pfrag, int nxc, int *xc, int *cnt)` |
| `DiylFrag_Diff` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2098` | internal | tree-sitter+gcc-aux | `static int DiylFrag_Diff (DiylFrag *pfrag1, DiylFrag *pfrag2)` |
| `DiylFrag_DebugTrace` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2119` | internal | tree-sitter+gcc-aux | `static void DiylFrag_DebugTrace (DiylFrag *pfrag)` |
| `analyze_CRU_folding` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2146` | internal | tree-sitter+gcc-aux | `static int analyze_CRU_folding (ORIG_ATOM_DATA *orig_at_data, int iunit, int n_all_bkb, int *all_bkb, int nxclasses, int *xc, OAD_StructureEdits *ed)` |
| `count_colors_in_sequence` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2463` | internal | tree-sitter+gcc-aux | `static int count_colors_in_sequence (int *color, int n, int maxcol, int *counts)` |
| `len_repeating_subsequence` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2489` | internal | tree-sitter+gcc-aux | `static int len_repeating_subsequence (int *color, int *color2, int n)` |
| `OAD_Polymer_PrepareFrameShiftEdits` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2523` | external | tree-sitter+gcc-aux | `extern int OAD_Polymer_PrepareFrameShiftEdits (ORIG_ATOM_DATA *orig_at_data, char *sinchi, char *saux, OAD_StructureEdits *ed)` |
| `ModSCenter_Init` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2760` | internal | tree-sitter+gcc-aux | `static void ModSCenter_Init (ModSCenterInfo *scinfo, inp_ATOM *at, int iatom)` |
| `NDefStereoBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2775` | internal | tree-sitter+gcc-aux | `static int NDefStereoBonds (inp_ATOM *at, int iatom, int bOnlyPointedEndMatters)` |
| `ModSCenter_AddTo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2803` | internal | tree-sitter+gcc-aux | `static void ModSCenter_AddTo (ModSCenterInfo *scinfo, int iadd)` |
| `ModSCenter_DelFrom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2815` | internal | tree-sitter+gcc-aux | `static void ModSCenter_DelFrom (ModSCenterInfo *scinfo, int idel)` |
| `ModSCenter_IsChanged` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi2.c:2836` | internal | tree-sitter+gcc-aux | `static int ModSCenter_IsChanged (ModSCenterInfo *scinfo, inp_ATOM *at)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c`

Parse errors: `8`. Function definitions: `51`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `OrigAtData_bCheckUnusualValences` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:89` | internal | tree-sitter+gcc-aux | `static int OrigAtData_bCheckUnusualValences (ORIG_ATOM_DATA *orig_at_data, int bAddIsoH, char *pStrErrStruct, int bNoWarnings)` |
| `OrigAtData_Duplicate` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:148` | external | tree-sitter+gcc-aux | `extern int OrigAtData_Duplicate (ORIG_ATOM_DATA *new_orig_atom, ORIG_ATOM_DATA *orig_atom)` |
| `PreprocessOneStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:488` | external | tree-sitter+gcc-aux | `extern int PreprocessOneStructure (struct tagINCHI_CLOCK *ic, STRUCT_DATA *sd, INPUT_PARMS *ip, ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data)` |
| `CreateCompositeNormAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:846` | external | tree-sitter | `int CreateCompositeNormAtom( COMP_ATOM_DATA *composite_norm_data, INP_ATOM_DATA2 *all_inp_norm_data, int num_components )` |
| `OrigAtData_DebugTrace` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1108` | external | tree-sitter+gcc-aux | `extern void OrigAtData_DebugTrace (ORIG_ATOM_DATA *d)` |
| `OAD_PolymerUnit_New` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1162` | external | tree-sitter+gcc-aux | `extern OAD_PolymerUnit *OAD_PolymerUnit_New (int maxatoms, int maxbonds, int id, int label, int type, int subtype, int conn, char *smt, int na, INT_ARRAY *alist, int nb, INT_ARRAY *blist, int nbkbonds, int **bkbonds)` |
| `OAD_PolymerUnit_CreateCopy` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1260` | external | tree-sitter+gcc-aux | `extern OAD_PolymerUnit *OAD_PolymerUnit_CreateCopy (OAD_PolymerUnit *u)` |
| `OAD_PolymerUnit_Free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1342` | external | tree-sitter+gcc-aux | `extern void OAD_PolymerUnit_Free (OAD_PolymerUnit *unit)` |
| `OAD_PolymerUnit_CompareAtomListsMod` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1377` | external | tree-sitter+gcc-aux | `extern int OAD_PolymerUnit_CompareAtomListsMod (OAD_PolymerUnit *u1, OAD_PolymerUnit *u2)` |
| `OAD_PolymerUnit_CompareAtomLists` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1400` | external | tree-sitter+gcc-aux | `extern int OAD_PolymerUnit_CompareAtomLists (OAD_PolymerUnit *u1, OAD_PolymerUnit *u2)` |
| `OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1437` | external | tree-sitter+gcc-aux | `extern int OAD_PolymerUnit_OrderBondAtomsAndBondsThemselves (OAD_PolymerUnit *u, int n_star_atoms, int *star_atoms)` |
| `OAD_ValidatePolymerAndPseudoElementData` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1516` | external | tree-sitter+gcc-aux | `extern int OAD_ValidatePolymerAndPseudoElementData (ORIG_ATOM_DATA *orig_at_data, int treat_polymers, int bNPZz, char *pStrErr, int bNoWarnings)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1777` | external | tree-sitter | `else if (representation == POLYMER_REPRESENTATION_STRUCTURE_BASED) #endif` |
| `UnMarkRingSystemsInp` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1961` | external | tree-sitter+gcc-aux | `extern int UnMarkRingSystemsInp (inp_ATOM *at, int num_atoms)` |
| `OAD_Polymer_CyclizeCloseableUnits` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:1979` | external | tree-sitter+gcc-aux | `extern int OAD_Polymer_CyclizeCloseableUnits (ORIG_ATOM_DATA *orig_at_data, int treat_polymers, char *pStrErr, int bNoWarnings)` |
| `OAD_PolymerUnit_HasMetal` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2055` | external | tree-sitter+gcc-aux | `extern int OAD_PolymerUnit_HasMetal (OAD_PolymerUnit *u, inp_ATOM *at)` |
| `OAD_Polymer_Free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2071` | external | tree-sitter+gcc-aux | `extern void OAD_Polymer_Free (OAD_Polymer *pd)` |
| `OAD_PolymerUnit_UnlinkCapsAndConnectEndAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2101` | external | tree-sitter+gcc-aux | `extern void OAD_PolymerUnit_UnlinkCapsAndConnectEndAtoms (OAD_PolymerUnit *unit, ORIG_ATOM_DATA *orig_inp_data, int *err, char *pStrErr)` |
| `OAD_PolymerUnit_FindEndsAndCaps` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2166` | external | tree-sitter+gcc-aux | `extern void OAD_PolymerUnit_FindEndsAndCaps (OAD_PolymerUnit *unit, ORIG_ATOM_DATA *orig_at_data, int *end1, int *cap1, int *cap1_is_star, int *end2, int *cap2, int *cap2_is_star, int *err, char *pStrErr)` |
| `OAD_PolymerUnit_SetEndsAndCaps` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2266` | external | tree-sitter+gcc-aux | `extern void OAD_PolymerUnit_SetEndsAndCaps (OAD_PolymerUnit *unit, ORIG_ATOM_DATA *orig_at_data, int *err, char *pStrErr)` |
| `OAD_Polymer_PrepareWorkingSet` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2342` | external | tree-sitter+gcc-aux | `extern int OAD_Polymer_PrepareWorkingSet (OAD_Polymer *p, int *cano_nums, int *compnt_nums, OAD_PolymerUnit **units2, int *unum)` |
| `OrigAtData_RemoveHalfBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2517` | external | tree-sitter+gcc-aux | `extern int OrigAtData_RemoveHalfBond (int this_atom, int other_atom, inp_ATOM *at, int *bond_type, int *bond_stereo)` |
| `OrigAtData_RemoveAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2564` | external | tree-sitter+gcc-aux | `extern int OrigAtData_RemoveAtom (ORIG_ATOM_DATA *orig_at_data, int iatom)` |
| `OrigAtData_RemoveBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2577` | external | tree-sitter+gcc-aux | `extern int OrigAtData_RemoveBond (int this_atom, int other_atom, inp_ATOM *at, int *bond_type, int *bond_stereo, int *num_inp_bonds)` |
| `OrigAtData_AddBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2607` | external | tree-sitter+gcc-aux | `extern int OrigAtData_AddBond (int this_atom, int other_atom, inp_ATOM *at, int bond_type, int bond_stereo, int *num_bonds)` |
| `OrigAtData_AddSingleStereolessBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2682` | external | tree-sitter+gcc-aux | `extern int OrigAtData_AddSingleStereolessBond (int this_atom, int other_atom, inp_ATOM *at, int *num_bonds)` |
| `OrigAtData_IncreaseBondOrder` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2692` | external | tree-sitter+gcc-aux | `extern int OrigAtData_IncreaseBondOrder (int this_atom, int other_atom, inp_ATOM *at)` |
| `OrigAtData_DecreaseBondOrder` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2750` | external | tree-sitter+gcc-aux | `extern int OrigAtData_DecreaseBondOrder (int this_atom, int other_atom, inp_ATOM *at)` |
| `OAD_CollectFragmentBondsAndAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2805` | external | tree-sitter+gcc-aux | `extern void OAD_CollectFragmentBondsAndAtoms (ORIG_ATOM_DATA *orig_at_data, int nforbidden, int *forbidden_orig, int *n_fragbonds, int **fragbonds, int *n_fragatoms, int *fragatoms, int *err, char *pStrErr)` |
| `OAD_Polymer_FindBackbones` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2872` | external | tree-sitter+gcc-aux | `extern void OAD_Polymer_FindBackbones (ORIG_ATOM_DATA *at_data, COMP_ATOM_DATA *composite_norm_data, int *err, char *pStrErr)` |
| `OAD_CollectBackboneAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:2946` | external | tree-sitter+gcc-aux | `extern void OAD_CollectBackboneAtoms (ORIG_ATOM_DATA *at_data, int na, int *alist, int end_atom1, int end_atom2, int *nbkatoms, int *bkatoms, int *err, char *pStrErr)` |
| `OAD_CollectReachableAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3019` | external | tree-sitter+gcc-aux | `extern int OAD_CollectReachableAtoms (ORIG_ATOM_DATA *orig_at_data, int start_atom, int nforbidden_bonds, int *forbidden_bond_atoms, int *n_reachable, int *reachable, int *err, char *pStrErr)` |
| `OAD_CollectBackboneBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3109` | external | tree-sitter+gcc-aux | `extern void OAD_CollectBackboneBonds (ORIG_ATOM_DATA *at_data, int na, int *alist, int end_atom1, int end_atom2, int *nbkbonds, int **bkbonds, int *err, char *pStrErr)` |
| `OAD_PolymerUnit_DelistIntraRingBackboneBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3168` | external | tree-sitter+gcc-aux | `extern void OAD_PolymerUnit_DelistIntraRingBackboneBonds (OAD_PolymerUnit *unit, ORIG_ATOM_DATA *at_data, int *err, char *pStrErr)` |
| `OAD_Polymer_FindRingSystems` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3236` | external | tree-sitter+gcc-aux | `extern int OAD_Polymer_FindRingSystems (OAD_Polymer *pd, inp_ATOM *at, int nat, int *num_inp_bonds, int *num_ring_sys, int *size_ring_sys, int start)` |
| `OAD_Polymer_SetAtProps` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3312` | external | tree-sitter+gcc-aux | `extern void OAD_Polymer_SetAtProps (OAD_Polymer *pd, inp_ATOM *at, int nat, int *num_inp_bonds, OAD_AtProps *aprops, int *cano_nums)` |
| `OAD_PolymerUnit_DelistHighOrderBackboneBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3509` | external | tree-sitter+gcc-aux | `extern void OAD_PolymerUnit_DelistHighOrderBackboneBonds (OAD_PolymerUnit *unit, ORIG_ATOM_DATA *orig_at_data, COMP_ATOM_DATA *composite_norm_data, int *err, char *pStrErr)` |
| `OAD_PolymerUnit_RemoveLinkFromCRUChain` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3607` | internal | tree-sitter+gcc-aux | `static void OAD_PolymerUnit_RemoveLinkFromCRUChain (int at1, int at2, int *nbonds, int **bonds)` |
| `OAD_PolymerUnit_DebugTrace` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3642` | external | tree-sitter+gcc-aux | `extern void OAD_PolymerUnit_DebugTrace (OAD_PolymerUnit *u)` |
| `OAD_Polymer_DebugTrace` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3757` | external | tree-sitter+gcc-aux | `extern void OAD_Polymer_DebugTrace (OAD_Polymer *p)` |
| `OAD_Polymer_GetRepresentation` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3788` | external | tree-sitter+gcc-aux | `extern int OAD_Polymer_GetRepresentation (OAD_Polymer *p)` |
| `OAD_Polymer_SmartReopenCyclizedUnits` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3878` | external | tree-sitter+gcc-aux | `extern void OAD_Polymer_SmartReopenCyclizedUnits (OAD_Polymer *p, inp_ATOM *at, int nat, int *num_inp_bonds)` |
| `OAD_PolymerUnit_ReopenCyclized` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:3953` | external | tree-sitter+gcc-aux | `extern void OAD_PolymerUnit_ReopenCyclized (OAD_PolymerUnit *u, inp_ATOM *at, OAD_AtProps *aprops, int nat, int *num_inp_bonds)` |
| `OAD_PolymerUnit_SetReopeningDetails` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4007` | external | tree-sitter+gcc-aux | `extern int OAD_PolymerUnit_SetReopeningDetails (OAD_PolymerUnit *u, inp_ATOM *at)` |
| `OAD_PolymerUnit_SortBackboneBondsAndSetSeniors` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4056` | external | tree-sitter+gcc-aux | `extern void OAD_PolymerUnit_SortBackboneBondsAndSetSeniors (OAD_PolymerUnit *u, inp_ATOM *at, OAD_AtProps *aprops, int *senior_bond)` |
| `OAD_PolymerUnit_SortBackboneBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4097` | external | tree-sitter+gcc-aux | `extern void OAD_PolymerUnit_SortBackboneBonds (OAD_PolymerUnit *u, OAD_AtProps *aprops, int *bnum)` |
| `OAD_Polymer_CompareBackboneBondsSeniority` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4128` | external | tree-sitter+gcc-aux | `extern int OAD_Polymer_CompareBackboneBondsSeniority (int *b1, int *b2, OAD_AtProps *aprops)` |
| `OAD_Polymer_CompareRanksOfTwoAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4209` | external | tree-sitter+gcc-aux | `extern int OAD_Polymer_CompareRanksOfTwoAtoms (int atom1, int atom2, OAD_AtProps *aprops)` |
| `OAD_Polymer_IsFirstAtomRankLower` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4369` | external | tree-sitter+gcc-aux | `extern int OAD_Polymer_IsFirstAtomRankLower (int atom1, int atom2, OAD_AtProps *aprops)` |
| `OAD_ValidateAndSortOutPseudoElementAtoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4396` | external | tree-sitter+gcc-aux | `extern void OAD_ValidateAndSortOutPseudoElementAtoms (ORIG_ATOM_DATA *orig_at_data, int treat_polymers, int bNPZz, int *err, char *pStrErr)` |
| `Inp_Atom_GetBondType` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi3.c:4507` | external | tree-sitter+gcc-aux | `extern int Inp_Atom_GetBondType (inp_ATOM *at, int iatom1, int iatom2)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c`

Parse errors: `11`. Function definitions: `15`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `SortAndPrintINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:102` | external | tree-sitter+gcc-aux | `extern int SortAndPrintINChI (CANON_GLOBALS *pCG, INCHI_IOSTREAM *out_file, INCHI_IOS_STRING *strbuf, INCHI_IOSTREAM *log_file, INPUT_PARMS *ip, ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data, COMP_ATOM_DATA (*composite_norm_data)[3], ORIG_STRUCT *pOrigStruct, int *num_components, int *num_non_taut, int *num_taut, INCHI_MODE *bTautFlags, INCHI_MODE *bTautFlagsDone, NORM_CANON_FLAGS *pncFlags, long int num_inp, PINChI2 (**pINChI), PINChI_Aux2 (**pINChI_Aux), int *pSortPrintINChIFlags, unsigned char save_opt_bits)` |
| `winchi_calc_inchikey` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:421` | external | tree-sitter+gcc-aux | `extern void winchi_calc_inchikey (int ret, int *ikflag, INPUT_PARMS *ip, INCHI_IOSTREAM *out_file, INCHI_IOSTREAM *log_file)` |
| `SaveEquComponentsInfoAndSortOrder` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:503` | external | tree-sitter | `int SaveEquComponentsInfoAndSortOrder( int iINChI, INCHI_SORT *pINChISort[TAUT_NUM], int *num_components, ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data, COMP_ATOM_DATA composite_norm_data[TAUT_NUM + 1], int bCompareComponents )` |
| `DisplayTheWholeCompositeStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:703` | external | tree-sitter | `int DisplayTheWholeCompositeStructure( struct tagCANON_GLOBALS *pCG, struct tagINCHI_CLOCK *ic, INPUT_PARMS *ip, struct tagStructData *sd, long num_inp, int iINChI, PINChI2 *pINChI2, PINChI_Aux2 *pINChI_Aux2, ORIG_ATOM_DATA *orig_inp_data, ORIG_ATOM_DATA *prep_inp_data, COMP_ATOM_DATA composite_norm_data[TAUT_NUM + 1] )` |
| `DisplayTheWholeStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:994` | external | tree-sitter | `int DisplayTheWholeStructure( struct tagCANON_GLOBALS *pCG, struct tagINCHI_CLOCK *ic, STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle, INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file, ORIG_ATOM_DATA *orig_inp_data, long num_inp, int iINChI, int bShowStruct, int bINCHI_LIB_Flag )` |
| `DisplayTheWholeStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1204` | external | tree-sitter+gcc-aux | `extern int DisplayTheWholeStructure (struct tagCANON_GLOBALS *pCG, struct tagINCHI_CLOCK *ic, STRUCT_DATA *sd, INPUT_PARMS *ip, char *szTitle, INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file, ORIG_ATOM_DATA *orig_inp_data, long int num_inp, int iINChI, int bShowStruct, int bINCHI_LIB_Flag)` |
| `DisplayStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1222` | external | tree-sitter+gcc-aux | `extern int DisplayStructure (struct tagCANON_GLOBALS *pCG, inp_ATOM *at, int num_at, OAD_Polymer *polymer, int num_removed_H, int bAdd_DT_to_num_H, int nNumRemovedProtons, NUM_H *nNumRemovedProtonsIsotopic, int bIsotopic, int j, INChI **cur_INChI, INChI_Aux **cur_INChI_Aux, int bAbcNumbers, DRAW_PARMS *dp, INCHI_MODE nMode, char *szTitle)` |
| `SplitTime` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1252` | external | tree-sitter | `void SplitTime( unsigned long ulTotalTime, int *hours, int *minutes, int *seconds, int *mseconds )` |
| `bIsStructChiral` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1271` | external | tree-sitter+gcc-aux | `extern int bIsStructChiral (PINChI2 (**pINChI2), int *num_components)` |
| `FreeAllINChIArrays` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1322` | external | tree-sitter+gcc-aux | `extern void FreeAllINChIArrays (PINChI2 (**pINChI), PINChI_Aux2 (**pINChI_Aux), int *num_components)` |
| `FreeINChIArrays` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1357` | external | tree-sitter+gcc-aux | `extern void FreeINChIArrays (PINChI2 (*pINChI), PINChI_Aux2 (*pINChI_Aux), int num_components)` |
| `TreatErrorsInCreateOneComponentINChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1387` | external | tree-sitter+gcc-aux | `extern int TreatErrorsInCreateOneComponentINChI (STRUCT_DATA *sd, INPUT_PARMS *ip, ORIG_ATOM_DATA *orig_inp_data, int i, long int num_inp, INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *out_file, INCHI_IOSTREAM *prb_file)` |
| `TreatCreateINChIWarning` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1453` | external | tree-sitter+gcc-aux | `extern int TreatCreateINChIWarning (STRUCT_DATA *sd, INPUT_PARMS *ip, ORIG_ATOM_DATA *orig_inp_data, long int num_inp, INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *out_file, INCHI_IOSTREAM *prb_file)` |
| `GetProcessingWarningsOneComponentInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1550` | external | tree-sitter+gcc-aux | `extern int GetProcessingWarningsOneComponentInChI (INChI **cur_INChI, INP_ATOM_DATA **inp_norm_data, STRUCT_DATA *sd, int bNoWarnings)` |
| `GetProcessingWarningsOneInChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/runichi4.c:1574` | external | tree-sitter+gcc-aux | `extern int GetProcessingWarningsOneInChI (INChI *pINChI, INP_ATOM_DATA *inp_norm_data, char *pStrErrStruct, int bNoWarnings)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c`

Parse errors: `0`. Function definitions: `9`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `sha2_starts` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:57` | external | tree-sitter+gcc-aux | `extern void sha2_starts (sha2_context *ctx)` |
| `sha2_process` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:72` | internal | tree-sitter+gcc-aux | `static void sha2_process (sha2_context *ctx, unsigned char *data)` |
| `sha2_update` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:206` | external | tree-sitter+gcc-aux | `extern void sha2_update (sha2_context *ctx, unsigned char *input, int ilen)` |
| `sha2_finish` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:257` | external | tree-sitter+gcc-aux | `extern void sha2_finish (sha2_context *ctx, unsigned char *output)` |
| `sha2_file` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:288` | external | tree-sitter+gcc-aux | `extern int sha2_file (char *path, unsigned char *output)` |
| `sha2_csum` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:312` | external | tree-sitter+gcc-aux | `extern void sha2_csum (unsigned char *input, int ilen, unsigned char *output)` |
| `sha2_hmac` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:325` | external | tree-sitter+gcc-aux | `extern void sha2_hmac (unsigned char *key, int keylen, unsigned char *input, int ilen, unsigned char *output)` |
| `sha2_self_test` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:391` | external | tree-sitter | `int sha2_self_test(void)` |
| `sha2_self_test` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/sha2.c:429` | external | tree-sitter+gcc-aux | `extern int sha2_self_test (void)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c`

Parse errors: `4`. Function definitions: `53`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `cmp_iso_atw_diff_component_no` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:166` | external | tree-sitter+gcc-aux | `extern int cmp_iso_atw_diff_component_no (const void *a1, const void *a2)` |
| `the_only_doublet_neigh` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:179` | external | tree-sitter+gcc-aux | `extern int the_only_doublet_neigh (inp_ATOM *at, int i1, int *ineigh1, int *ineigh2)` |
| `fix_non_uniform_drawn_oxoanions` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:227` | external | tree-sitter+gcc-aux | `extern int fix_non_uniform_drawn_oxoanions (int num_atoms, inp_ATOM *at, int *num_changes)` |
| `fix_non_uniform_drawn_amidiniums` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:452` | external | tree-sitter+gcc-aux | `extern int fix_non_uniform_drawn_amidiniums (int num_atoms, inp_ATOM *at, int *num_changes)` |
| `fix_odd_things` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:603` | external | tree-sitter+gcc-aux | `extern int fix_odd_things (int num_atoms, inp_ATOM *at, int bFixBug, int bFixNonUniformDraw)` |
| `post_fix_odd_things` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:922` | external | tree-sitter+gcc-aux | `extern int post_fix_odd_things (int num_atoms, inp_ATOM *at)` |
| `nFindOneOM` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:932` | external | tree-sitter+gcc-aux | `extern int nFindOneOM (inp_ATOM *at, int at_no, int *ord_OM, int num_OM)` |
| `remove_ion_pairs` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:1049` | external | tree-sitter+gcc-aux | `extern int remove_ion_pairs (int num_atoms, inp_ATOM *at)` |
| `RemoveInpAtBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2046` | external | tree-sitter+gcc-aux | `extern int RemoveInpAtBond (inp_ATOM *atom, int iat, int k)` |
| `DisconnectInpAtBond` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2261` | external | tree-sitter+gcc-aux | `extern int DisconnectInpAtBond (inp_ATOM *at, AT_NUMB *nOldCompNumber, int iat, int neigh_ord)` |
| `bIsAmmoniumSalt` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2300` | external | tree-sitter+gcc-aux | `extern int bIsAmmoniumSalt (inp_ATOM *at, int i, int *piO, int *pk, S_CHAR *num_explicit_H)` |
| `DisconnectAmmoniumSalt` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2398` | external | tree-sitter+gcc-aux | `extern int DisconnectAmmoniumSalt (inp_ATOM *at, int iN, int iO, int k, S_CHAR *num_explicit_H)` |
| `bIsMetalSalt` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2511` | external | tree-sitter+gcc-aux | `extern int bIsMetalSalt (inp_ATOM *at, int i)` |
| `DisconnectMetalSalt` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2612` | external | tree-sitter+gcc-aux | `extern int DisconnectMetalSalt (inp_ATOM *at, int i)` |
| `DisconnectSalts` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2668` | external | tree-sitter+gcc-aux | `extern int DisconnectSalts (ORIG_ATOM_DATA *orig_inp_data, int bDisconnect)` |
| `bIsMetalToDisconnect` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2719` | external | tree-sitter+gcc-aux | `extern int bIsMetalToDisconnect (inp_ATOM *at, int i, int bCheckMetalValence)` |
| `bMayDisconnectMetals` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2762` | external | tree-sitter+gcc-aux | `extern int bMayDisconnectMetals (ORIG_ATOM_DATA *orig_inp_data, int bCheckMetalValence, INCHI_MODE *bTautFlagsDone)` |
| `bHasMetalAtom` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2820` | external | tree-sitter | `int bHasMetalAtom( ORIG_ATOM_DATA *orig_inp_data )` |
| `DisconnectMetals` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:2865` | external | tree-sitter+gcc-aux | `extern int DisconnectMetals (ORIG_ATOM_DATA *orig_inp_data, int bCheckMetalValence, INCHI_MODE *bTautFlagsDone)` |
| `DisconnectOneLigand` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:3112` | external | tree-sitter+gcc-aux | `extern int DisconnectOneLigand (inp_ATOM *at, AT_NUMB *nOldCompNumber, S_CHAR *bMetal, char *elnumber_Heteroat, int num_halogens, int num_atoms, int iMetal, int jLigand, INCHI_MODE *bTautFlagsDone)` |
| `dist3D` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:3245` | external | tree-sitter+gcc-aux | `extern double dist3D (inp_ATOM *at1, inp_ATOM *at2)` |
| `GetMinDistDistribution` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:3264` | external | tree-sitter+gcc-aux | `extern double GetMinDistDistribution (inp_ATOM *at, int num_at, int iat, int iat_H, int bInAllComponents, double *min_dist, int num_segm)` |
| `move_explicit_Hcation` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:3480` | external | tree-sitter+gcc-aux | `extern int move_explicit_Hcation (inp_ATOM *at, int num_at, int iat, int iat_H, int bInAllComponents)` |
| `add_DT_to_num_H` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:3705` | external | tree-sitter+gcc-aux | `extern int add_DT_to_num_H (int num_atoms, inp_ATOM *at)` |
| `remove_terminal_HDT` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:3723` | external | tree-sitter+gcc-aux | `extern int remove_terminal_HDT (int num_atoms, inp_ATOM *at, int bFixTermHChrg)` |
| `get_iat_number` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4024` | external | tree-sitter+gcc-aux | `extern int get_iat_number (int el_number)` |
| `bHeteroAtomMayHaveXchgIsoH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4047` | external | tree-sitter+gcc-aux | `extern int bHeteroAtomMayHaveXchgIsoH (inp_ATOM *atom, int iat)` |
| `bNumHeterAtomHasIsotopicH` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4133` | external | tree-sitter+gcc-aux | `extern int bNumHeterAtomHasIsotopicH (inp_ATOM *atom, int num_atoms)` |
| `cmp_components` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4247` | external | tree-sitter+gcc-aux | `extern int cmp_components (const void *a1, const void *a2)` |
| `MarkDisconnectedComponents` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4277` | external | tree-sitter+gcc-aux | `extern int MarkDisconnectedComponents (ORIG_ATOM_DATA *orig_at_data, int bProcessOldCompNumbers)` |
| `ExtractConnectedComponent` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4558` | external | tree-sitter+gcc-aux | `extern int ExtractConnectedComponent (inp_ATOM *at, int num_at, int component_number, inp_ATOM *component_at)` |
| `SetConnectedComponentNumber` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4598` | external | tree-sitter+gcc-aux | `extern int SetConnectedComponentNumber (inp_ATOM *at, int num_at, int component_number)` |
| `Free_INChI_Stereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4611` | external | tree-sitter+gcc-aux | `extern int Free_INChI_Stereo (INChI_Stereo *pINChI_Stereo)` |
| `Alloc_INChI_Stereo` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4629` | external | tree-sitter+gcc-aux | `extern INChI_Stereo *Alloc_INChI_Stereo (int num_at, int num_bonds)` |
| `Free_INChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4675` | external | tree-sitter+gcc-aux | `extern int Free_INChI (INChI **ppINChI)` |
| `Free_INChI_Members` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4698` | external | tree-sitter+gcc-aux | `extern int Free_INChI_Members (INChI *pINChI)` |
| `Alloc_INChI` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4722` | external | tree-sitter+gcc-aux | `extern INChI *Alloc_INChI (inp_ATOM *at, int num_at, int *found_num_bonds, int *found_num_isotopic, int nAllocMode)` |
| `Free_INChI_Aux` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4843` | external | tree-sitter+gcc-aux | `extern int Free_INChI_Aux (INChI_Aux **ppINChI_Aux)` |
| `Alloc_INChI_Aux` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4884` | external | tree-sitter+gcc-aux | `extern INChI_Aux *Alloc_INChI_Aux (int num_at, int num_isotopic_atoms, int nAllocMode, int bOrigCoord)` |
| `CompAtomData_GetNumMapping` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:4986` | external | tree-sitter+gcc-aux | `extern void CompAtomData_GetNumMapping (COMP_ATOM_DATA *adata, int *orig_num, int *curr_num)` |
| `imat_new` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5005` | external | tree-sitter+gcc-aux | `extern int imat_new (int m, int n, int ***a)` |
| `imat_free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5036` | external | tree-sitter+gcc-aux | `extern void imat_free (int m, int **a)` |
| `subgraf_new` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5064` | external | tree-sitter+gcc-aux | `extern subgraf *subgraf_new (ORIG_ATOM_DATA *orig_inp_data, int nnodes, int *nodes)` |
| `subgraf_free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5155` | external | tree-sitter+gcc-aux | `extern void subgraf_free (subgraf *sg)` |
| `subgraf_debug_trace` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5192` | external | tree-sitter+gcc-aux | `extern void subgraf_debug_trace (subgraf *sg)` |
| `subgraf_pathfinder_new` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5221` | external | tree-sitter+gcc-aux | `extern subgraf_pathfinder *subgraf_pathfinder_new (subgraf *sg, ORIG_ATOM_DATA *orig_inp_data, int start, int end)` |
| `subgraf_pathfinder_free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5252` | external | tree-sitter+gcc-aux | `extern void subgraf_pathfinder_free (subgraf_pathfinder *spf)` |
| `subgraf_pathfinder_run` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5273` | external | tree-sitter+gcc-aux | `extern void subgraf_pathfinder_run (subgraf_pathfinder *spf, int nforbidden, int *forbidden, int *nbonds, int **bonds, int *natoms, int *atoms)` |
| `add_bond_if_unseen` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5380` | external | tree-sitter+gcc-aux | `extern void add_bond_if_unseen (subgraf_pathfinder *spf, int node0, int node, int *nbonds, int **bonds)` |
| `subgraf_pathfinder_collect_all` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5422` | external | tree-sitter+gcc-aux | `extern int subgraf_pathfinder_collect_all (subgraf_pathfinder *spf, int nforbidden, int *forbidden, int *atnums)` |
| `FixNextRadicals` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5471` | external | tree-sitter | `int FixNextRadicals( int cur_at, inp_ATOM *at )` |
| `FixAdjacentRadicals` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5491` | external | tree-sitter | `int FixAdjacentRadicals( int num_inp_atoms, inp_ATOM *at )` |
| `PrintFileName` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/strutil.c:5533` | external | tree-sitter | `void PrintFileName( const char *fmt, FILE *out_file, /* INCHI_IOSTREAM *out_file, */ const char *szFname )` |

### `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c`

Parse errors: `13`. Function definitions: `49`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `get_element_chemical_symbol` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:289` | external | tree-sitter+gcc-aux | `extern int get_element_chemical_symbol (int nAtNum, char *szElement)` |
| `get_element_or_pseudoelement_symbol` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:316` | external | tree-sitter+gcc-aux | `extern int get_element_or_pseudoelement_symbol (int nAtNum, char *szElement)` |
| `el_number_in_internal_ref_table` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:347` | external | tree-sitter+gcc-aux | `extern int el_number_in_internal_ref_table (const char *elname)` |
| `get_periodic_table_number` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:364` | external | tree-sitter+gcc-aux | `extern int get_periodic_table_number (const char *elname)` |
| `if_skip_add_H` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:428` | external | tree-sitter+gcc-aux | `extern int if_skip_add_H (int nPeriodicNum)` |
| `get_el_valence` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:439` | external | tree-sitter+gcc-aux | `extern int get_el_valence (int nPeriodicNum, int charge, int val_num)` |
| `get_unusual_el_valence` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:454` | external | tree-sitter+gcc-aux | `extern int get_unusual_el_valence (int nPeriodicNum, int charge, int radical, int bonds_valence, int num_H, int num_bonds)` |
| `needed_unusual_el_valence` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:518` | external | tree-sitter+gcc-aux | `extern int needed_unusual_el_valence (int nPeriodicNum, int charge, int radical, int bonds_valence, int actual_bonds_valence, int num_H, int num_bonds)` |
| `detect_unusual_el_valence` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:620` | external | tree-sitter+gcc-aux | `extern int detect_unusual_el_valence (int nPeriodicNum, int charge, int radical, int bonds_valence, int num_H, int num_bonds)` |
| `get_el_type` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:679` | external | tree-sitter+gcc-aux | `extern int get_el_type (int nPeriodicNum)` |
| `is_el_a_metal` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:688` | external | tree-sitter+gcc-aux | `extern int is_el_a_metal (int nPeriodicNum)` |
| `extract_charges_and_radicals` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:700` | external | tree-sitter+gcc-aux | `extern int extract_charges_and_radicals (char *elname, int *pnRadical, int *pnCharge)` |
| `extract_H_atoms` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:774` | external | tree-sitter+gcc-aux | `extern int extract_H_atoms (char *elname, S_CHAR *num_iso_H)` |
| `get_num_H` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:862` | external | tree-sitter+gcc-aux | `extern int get_num_H (const char *elname, int inp_num_H, S_CHAR *inp_num_iso_H, int charge, int radical, int chem_bonds_valence, int atom_input_valence, int bAliased, int bDoNotAddH, int bHasMetalNeighbor)` |
| `get_atomic_mass_from_elnum` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1007` | external | tree-sitter+gcc-aux | `extern int get_atomic_mass_from_elnum (int nAtNum)` |
| `get_atomic_mass` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1040` | external | tree-sitter+gcc-aux | `extern int get_atomic_mass (const char *elname)` |
| `is_in_the_list` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1059` | external | tree-sitter+gcc-aux | `extern AT_NUMB *is_in_the_list (AT_NUMB *pathAtom, AT_NUMB nNextAtom, int nPathLen)` |
| `is_in_the_ilist` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1072` | external | tree-sitter+gcc-aux | `extern int *is_in_the_ilist (int *pathAtom, int nNextAtom, int nPathLen)` |
| `is_ilist_inside` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1085` | external | tree-sitter+gcc-aux | `extern int is_ilist_inside (int *ilist, int nlist, int *ilist2, int nlist2)` |
| `nBondsValToMetal` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1100` | external | tree-sitter+gcc-aux | `extern int nBondsValToMetal (inp_ATOM *at, int iat)` |
| `num_of_H` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1129` | external | tree-sitter+gcc-aux | `extern int num_of_H (inp_ATOM *at, int iat)` |
| `ion_el_group` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1151` | external | tree-sitter+gcc-aux | `extern U_CHAR ion_el_group (int el)` |
| `has_other_ion_neigh` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1175` | external | tree-sitter+gcc-aux | `extern int has_other_ion_neigh (inp_ATOM *at, int iat, int iat_ion_neigh)` |
| `has_other_ion_in_sphere_2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1201` | external | tree-sitter+gcc-aux | `extern int has_other_ion_in_sphere_2 (inp_ATOM *at, int iat, int iat_ion_neigh)` |
| `nNoMetalNumBonds` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1253` | external | tree-sitter+gcc-aux | `extern int nNoMetalNumBonds (inp_ATOM *at, int at_no)` |
| `nNoMetalBondsValence` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1320` | external | tree-sitter+gcc-aux | `extern int nNoMetalBondsValence (inp_ATOM *at, int at_no)` |
| `nNoMetalNeighIndex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1386` | external | tree-sitter+gcc-aux | `extern int nNoMetalNeighIndex (inp_ATOM *at, int at_no)` |
| `nNoMetalOtherNeighIndex` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1405` | external | tree-sitter+gcc-aux | `extern int nNoMetalOtherNeighIndex (inp_ATOM *at, int at_no, int cur_neigh)` |
| `nNoMetalOtherNeighIndex2` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1426` | external | tree-sitter+gcc-aux | `extern int nNoMetalOtherNeighIndex2 (inp_ATOM *at, int at_no, int cur_neigh, int cur_neigh2)` |
| `MakeRemovedProtonsString` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1454` | external | tree-sitter | `int MakeRemovedProtonsString( int nNumRemovedProtons, NUM_H *nNumExchgIsotopicH, NUM_H *nNumRemovedProtonsIsotopic, int bIsotopic, char *szRemovedProtons, int *num_removed_iso_H )` |
| `get_endpoint_valence` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1509` | external | tree-sitter+gcc-aux | `extern int get_endpoint_valence (U_CHAR el_number)` |
| `get_endpoint_valence_KET` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1530` | external | tree-sitter+gcc-aux | `extern int get_endpoint_valence_KET (U_CHAR el_number)` |
| `inchi_malloc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1552` | external | tree-sitter | `void *inchi_malloc( size_t c )` |
| `inchi_calloc` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1561` | external | tree-sitter | `void *inchi_calloc( size_t c, size_t n )` |
| `inchi_free` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1570` | external | tree-sitter | `void inchi_free( void *p )` |
| `normalize_string` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1589` | external | tree-sitter+gcc-aux | `extern int normalize_string (char *name)` |
| `dotify_non_printable_chars` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1630` | external | tree-sitter+gcc-aux | `extern int dotify_non_printable_chars (char *line)` |
| `read_upto_delim` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1658` | external | tree-sitter+gcc-aux | `extern int read_upto_delim (char **pstring, char *field, int maxlen, char *delims)` |
| `is_matching_any_delim` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1710` | external | tree-sitter+gcc-aux | `extern int is_matching_any_delim (char c, char *delims)` |
| `remove_trailing_spaces` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1728` | external | tree-sitter+gcc-aux | `extern void remove_trailing_spaces (char *p)` |
| `remove_one_lf` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1740` | external | tree-sitter+gcc-aux | `extern void remove_one_lf (char *p)` |
| `mystrncpy` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1760` | external | tree-sitter+gcc-aux | `extern int mystrncpy (char *target, const char *source, unsigned int maxlen)` |
| `lrtrim` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1804` | external | tree-sitter+gcc-aux | `extern char *lrtrim (char *p, int *nLen)` |
| `extract_inchi_substring` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1860` | external | tree-sitter+gcc-aux | `extern void extract_inchi_substring (char **buf, const char *str, size_t slen)` |
| `extract_auxinfo_substring` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1921` | external | tree-sitter+gcc-aux | `extern void extract_auxinfo_substring (char **buf, const char *str, size_t slen)` |
| `inchi_memicmp` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1971` | external | tree-sitter+gcc-aux | `extern int inchi_memicmp (const void *p1, const void *p2, size_t length)` |
| `inchi_stricmp` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:1995` | external | tree-sitter+gcc-aux | `extern int inchi_stricmp (const char *s1, const char *s2)` |
| `inchi__strnset` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:2022` | external | tree-sitter+gcc-aux | `extern char *inchi__strnset (char *s, int val, size_t length)` |
| `inchi__strdup` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/util.c:2035` | external | tree-sitter+gcc-aux | `extern char *inchi__strdup (const char *string)` |


## CLI C Functions (Inventory Only)

### `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c`

Parse errors: `6`. Function definitions: `43`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `PrintFileName` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:154` | external | tree-sitter | `void PrintFileName( const char *fmt, FILE *out_file, const char *szFname )` |
| `ReallySetForegroundWindow` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:177` | external | tree-sitter | `BOOL ReallySetForegroundWindow( HWND hWnd )` |
| `GetConsoleHwnd` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:195` | external | tree-sitter | `HWND GetConsoleHwnd( void )` |
| `GetFontHeight` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:233` | external | tree-sitter | `int GetFontHeight( HDC pDC )` |
| `GetFontAscent` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:244` | external | tree-sitter | `int GetFontAscent( HDC pDC )` |
| `GetFontDescent` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:255` | external | tree-sitter | `int GetFontDescent( HDC pDC )` |
| `GetFontAveWidth` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:266` | external | tree-sitter | `int GetFontAveWidth( HDC pDC )` |
| `GetSubstringWidth` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:277` | external | tree-sitter | `int GetSubstringWidth( HDC pDC,int len, char *pString )` |
| `GetTextSize` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:307` | external | tree-sitter | `void GetTextSize( HDC pDC, int len, char *pString,int *width, int *height )` |
| `GetVertTextSize` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:318` | external | tree-sitter | `void GetVertTextSize( HDC pDC, int len, char *pString, int *width, int *height )` |
| `TextOutVert` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:337` | external | tree-sitter | `BOOL TextOutVert( HDC pDC, /* handle to DC */ int nXStart, /* x-coordinate of starting position */ int nYStart, /* y-coordinate of starting position */ LPCTSTR lpString, /* character string */ int cbString, /* number of characters */ int cell_width /* width for center alignment */ )` |
| `TextOutHoriz` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:362` | external | tree-sitter | `BOOL TextOutHoriz( HDC pDC, /* handle to DC */ int nXStart, /* x-coordinate of starting position*/ int nYStart, /* y-coordinate of starting position*/ LPCTSTR lpString, /* character string */ int cbString, /* number of characters */ int cell_width )` |
| `GetStringWidth` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:378` | external | tree-sitter | `int GetStringWidth( HDC pDC, char *pString )` |
| `GetOneCharInStringWidth` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:389` | external | tree-sitter | `int GetOneCharInStringWidth( HDC pDC, const char *pString )` |
| `DrawStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:400` | external | tree-sitter | `int DrawStructure( HDC pDC, inp_ATOM *at, INF_ATOM_DATA *inf_at_data, int num_at, int xoff, int yoff, COLORREF clrPen, int nPenWidth )` |
| `DrawBond` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:618` | external | tree-sitter | `int DrawBond( HDC pDC, int x1, int y1, int x2, int y2, int b_type, int b_stereo, int b_parity, int bInvertBonds, COLORREF clrPen, int nPenWidth )` |
| `RoundDouble` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:677` | external | tree-sitter | `double RoundDouble( double X )` |
| `nRound` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:684` | external | tree-sitter | `int nRound( double X )` |
| `roundoff_coord` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:691` | external | tree-sitter | `void roundoff_coord( double dx1, double dx2, int *new_ix1, int *new_ix2 )` |
| `DrawPenColorFilledPolygon` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:719` | external | tree-sitter | `void DrawPenColorFilledPolygon( HDC pDC, const POINT *pnt, int num )` |
| `DrawTextColorDot` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:746` | external | tree-sitter | `int DrawTextColorDot( HDC pDC )` |
| `DrawBondStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:812` | external | tree-sitter | `int DrawBondStereo( HDC pDC, int x1, int y1, int x2, int y2, int b_stereo, int b_highlight, int bInvertBonds, COLORREF clrPen, int nPenWidth )` |
| `DrawBondParity` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:993` | external | tree-sitter | `int DrawBondParity( HDC pDC, int x1, int y1, int x2, int y2, int parity_mark0 )` |
| `DrawBondCrossingCRU` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:1081` | external | tree-sitter | `int DrawBondCrossingCRU(HDC pDC, int x1, int y1, int x2, int y2)` |
| `DrawBondNoStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:1122` | external | tree-sitter | `int DrawBondNoStereo( HDC pDC, int x1, int y1, int x2, int y2, int b_type, int b_highlight, COLORREF clrPen, int nPenWidth )` |
| `MoveHydrogenAtomToTheLeft` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:1359` | external | tree-sitter | `int MoveHydrogenAtomToTheLeft( char *s,int start, int H )` |
| `MyTextOutABC` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:1391` | external | tree-sitter | `int MyTextOutABC( const char *p, int iFst, int iLst, HDC pDC )` |
| `DrawColorString` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:1431` | external | tree-sitter | `int DrawColorString( HDC pDC,const char *st, int xs, int ys, int bHighlightTheAtom )` |
| `DrawPreparedString` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:1675` | external | tree-sitter | `int DrawPreparedString( HDC pDC, char *st1, int shift, int x, int y, int bHighlightTheAtom )` |
| `DrawString` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:1688` | external | tree-sitter | `int DrawString( HDC pDC, char *st1, int shift, int x, int y )` |
| `DrawLine` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:1745` | external | tree-sitter | `void DrawLine( HDC pDC, int x1, int y1, int x2, int y2 )` |
| `nGetNumLegendOptions` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:1753` | external | tree-sitter | `int nGetNumLegendOptions( inf_ATOM *inf_at, int num_at )` |
| `DrawTheInputStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:1771` | external | tree-sitter | `int DrawTheInputStructure( inp_ATOM *at, INF_ATOM_DATA *inf_at_data, int num_at, HDC pDC, int tx_off, int ty_off, int xoff, int yoff, int width_pix, int height_pix, int bDraw, int bOrigAtom, COLORREF clrPen, int nPenWidth )` |
| `CalcTblParms` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:1887` | external | tree-sitter | `void CalcTblParms( HDC hMemoryDC, TBL_PARMS *tp, TBL_DRAW_PARMS *tdp, int *xStructOffs, int *yStructOffs, int *xStructSize, int *yStructSize, int yoffs1 )` |
| `DrawTheTable` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:1949` | external | tree-sitter | `int DrawTheTable( HDC hDC, TBL_PARMS *tp, TBL_DRAW_PARMS *tdp, int x_offs, int y_offs )` |
| `GetStructSizes` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:2045` | external | tree-sitter | `void GetStructSizes(HDC hDC, inf_ATOM* inf_at, inp_ATOM* at0, inp_ATOM* at1, int num_at, int* xoffs1, int* xoffs2, INT_DRAW_PARMS* idp)` |
| `ResizeAtomForDrawing` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:2223` | external | tree-sitter | `void ResizeAtomForDrawing( inf_ATOM *inf_at, inp_ATOM *at0, inp_ATOM *at1, int num_at, INT_DRAW_PARMS *idp, int width, int height, int nFontWidth, int *xoffs1, int *xoffs2, int *draw_width, int *draw_height, int *xdraw_offs, int *ydraw_offs )` |
| `InpStructureMarkEquComponents` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:2357` | external | tree-sitter | `void InpStructureMarkEquComponents( MY_WINDOW_DATA *pWinData, AT_NUMB nNewEquLabel, inp_ATOM *at0, inp_ATOM *at1, inf_ATOM *inf_at, int num_at )` |
| `CreateInputStructPicture` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:2448` | external | tree-sitter | `int CreateInputStructPicture( HDC hDC, MY_WINDOW_DATA *pWinData, RECT *rc, int bPrint, AT_NUMB nNewEquLabel )` |
| `WndProcDisplayInputStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:2919` | external | tree-sitter | `LRESULT CALLBACK WndProcDisplayInputStructure( HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam )` |
| `FreeWinData` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:3074` | external | tree-sitter | `void FreeWinData( MY_WINDOW_DATA* pWinData )` |
| `DisplayInputStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:3108` | external | tree-sitter | `int DisplayInputStructure( char *szOutputString, inp_ATOM *at, INF_ATOM_DATA *inf_at_data, int num_at, DRAW_PARMS *dp )` |
| `MySleep` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/dispstru.c:3321` | external | tree-sitter | `void MySleep( unsigned long ms )` |

### `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/ichimain.c`

Parse errors: `17`. Function definitions: `20`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `user_quit` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/ichimain.c:96` | external | tree-sitter | `int user_quit(struct tagINCHI_CLOCK* ic, const char* msg, unsigned long ulMaxTime)` |
| `eat_keyboard_input` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/ichimain.c:153` | external | tree-sitter | `void eat_keyboard_input(void)` |
| `MyHandlerRoutine` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/ichimain.c:190` | external | tree-sitter | `BOOL WINAPI MyHandlerRoutine(DWORD dwCtrlType /* control signal type */)` |
| `WasInterrupted` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/ichimain.c:205` | external | tree-sitter | `int WasInterrupted(void)` |
| `main` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/ichimain.c:224` | external | tree-sitter | `int main(int argc, char* argv[])` |
| `ProcessMultipleInputFiles` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/ichimain.c:266` | external | tree-sitter | `int ProcessMultipleInputFiles(int argc, char* argv[])` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/ichimain.c:397` | external | tree-sitter | `else if (pOutPath)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/ichimain.c:419` | external | tree-sitter | `else if (pOutPath)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/ichimain.c:440` | external | tree-sitter | `else if (pOutPath)` |
| `ProcessSingleInputFile` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/ichimain.c:523` | external | tree-sitter | `int ProcessSingleInputFile(int argc, char* argv[])` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/ichimain.c:930` | external | tree-sitter | `else if (orig_inp_data->num_inp_atoms == 4)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/ichimain.c:939` | external | tree-sitter | `else if (orig_inp_data->num_inp_atoms == 4)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/ichimain.c:943` | external | tree-sitter | `else if (orig_inp_data->num_inp_atoms == 5)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/ichimain.c:947` | external | tree-sitter | `else if (orig_inp_data->num_inp_atoms == 6)` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/ichimain.c:951` | external | tree-sitter | `else if (orig_inp_data->num_inp_atoms == 7)` |
| `save_command_line` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/ichimain.c:1044` | external | tree-sitter | `void save_command_line(int argc, char* argv[], INCHI_IOSTREAM* plog)` |
| `emit_empty_inchi` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/ichimain.c:1065` | external | tree-sitter | `void emit_empty_inchi(INPUT_PARMS* ip, long num_inp, char* pLF, char* pTAB, INCHI_IOSTREAM* pout)` |
| `GetTheNextRecordOfInputFile` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/ichimain.c:1085` | external | tree-sitter | `int GetTheNextRecordOfInputFile(struct tagINCHI_CLOCK* ic, STRUCT_DATA* sd, INPUT_PARMS* ip, char* szTitle, INCHI_IOSTREAM* inp_file, INCHI_IOSTREAM* plog, INCHI_IOSTREAM* pout, INCHI_IOSTREAM* pprb, ORIG_ATOM_DATA* orig_inp_data, long* num_inp, STRUCT_FPTRS* pStructPtrs, int* nRet, int* have_err_in_GetOneStructure, long* num_err, int output_error_inchi)` |
| `CalcAndPrintINCHIAndINCHIKEY` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/ichimain.c:1158` | external | tree-sitter | `int CalcAndPrintINCHIAndINCHIKEY(struct tagINCHI_CLOCK* ic, CANON_GLOBALS* CG, STRUCT_DATA* sd, INPUT_PARMS* ip, char* szTitle, PINChI2* pINChI[INCHI_NUM], PINChI_Aux2* pINChI_Aux[INCHI_NUM], INCHI_IOSTREAM* inp_file, INCHI_IOSTREAM* plog, INCHI_IOSTREAM* pout, INCHI_IOSTREAM* pprb, ORIG_ATOM_DATA* orig_inp_data, ORIG_ATOM_DATA* prep_inp_data, long* num_inp, STRUCT_FPTRS* pStructPtrs, int* nRet, int have_err_in_GetOneStructure, long* num_err, int output_error_inchi, INCHI_IOS_STRING* strbuf, unsigned long* pulTotalProcessingTime, char* pLF, char* pTAB, char* ikey, int silent)` |
| `RepeatedlyRenumberAtomsAndRecalcINCHI` | `third_party/InChI/INCHI-1-SRC/INCHI_EXE/inchi-1/src/ichimain.c:1470` | external | tree-sitter | `int RepeatedlyRenumberAtomsAndRecalcINCHI(struct tagINCHI_CLOCK* ic, CANON_GLOBALS* CG, STRUCT_DATA* sd, INPUT_PARMS* ip, char* szTitle, PINChI2* pINChI[INCHI_NUM], PINChI_Aux2* pINChI_Aux[INCHI_NUM], INCHI_IOSTREAM* inp_file, INCHI_IOSTREAM* plog, INCHI_IOSTREAM* pout, INCHI_IOSTREAM* pprb, ORIG_ATOM_DATA* orig_inp_data, ORIG_ATOM_DATA* prep_inp_data, long* num_inp, STRUCT_FPTRS* pStructPtrs, int* nRet, int have_err_in_GetOneStructure, long* num_err, int output_error_inchi, INCHI_IOS_STRING* strbuf, unsigned long* pulTotalProcessingTime, char* pLF, char* pTAB, long int nrepeat)` |


## Demo C Functions (Inventory Only)

### `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c`

Parse errors: `11`. Function definitions: `37`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `ee_extract_ChargeRadical` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:194` | external | tree-sitter | `int ee_extract_ChargeRadical( char *elname, int *pnRadical, int *pnCharge )` |
| `ee_extract_H_atoms` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:262` | external | tree-sitter | `int ee_extract_H_atoms( char *elname, S_CHAR num_iso_H[] )` |
| `e_GetElType` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:323` | external | tree-sitter | `int e_GetElType( inchi_Atom *at, int cur_atom )` |
| `e_inchi_swap` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:860` | external | tree-sitter | `void e_inchi_swap( char *a, char *b, size_t width )` |
| `e_insertions_sort` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:876` | external | tree-sitter | `int e_insertions_sort( void *base, size_t num, size_t width, int( *compare )( const void*, const void* ) ) /* djb-rwth: types of variables are sufficient */` |
| `e_get_z_coord` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:972` | external | tree-sitter | `double e_get_z_coord( inchi_Atom* at, int cur_atom, int neigh_no, int *nType, int bPointedEdgeStereo )` |
| `e_len3` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:1043` | external | tree-sitter | `double e_len3( const double c[] ) /* djb-rwth: avoiding uninitialised values */` |
| `e_len2` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:1056` | external | tree-sitter | `double e_len2( const double c[] ) /* djb-rwth: avoiding uninitialised values */` |
| `e_diff3` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:1069` | external | tree-sitter | `void* e_diff3( const double a[], const double b[], double result[] ) /* djb-rwth: changed function type */` |
| `e_add3` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:1081` | external | tree-sitter | `void e_add3( const double a[], const double b[], double result[] ) /* djb-rwth: changed function type */` |
| `e_mult3` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:1092` | external | tree-sitter | `void e_mult3( const double a[], double b, double result[] ) /* djb-rwth: changed function type */` |
| `e_change_sign3` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:1114` | external | tree-sitter | `void e_change_sign3( const double a[], double result[] ) /* djb-rwth: changed function type */` |
| `e_dot_prod3` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:1125` | external | tree-sitter | `double e_dot_prod3( const double a[], const double b[] )` |
| `e_dot_prodchar3` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:1132` | external | tree-sitter | `int e_dot_prodchar3( const S_CHAR a[], const S_CHAR b[] )` |
| `e_cross_prod3` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:1150` | external | tree-sitter | `double* e_cross_prod3( const double a[], const double b[], double result[] )` |
| `e_triple_prod` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:1167` | external | tree-sitter | `double e_triple_prod( double a[], double b[], double c[], double *sine_value )` |
| `e_CompDble` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:1202` | external | tree-sitter | `int e_CompDble( const void *a1, const void *a2 )` |
| `e_Get2DTetrahedralAmbiguity` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:1225` | external | tree-sitter | `int e_Get2DTetrahedralAmbiguity( double at_coord[][3], int bAddExplicitNeighbor )` |
| `e_triple_prod_and_min_abs_sine2` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:1665` | external | tree-sitter | `double e_triple_prod_and_min_abs_sine2( double at_coord[][3], double central_at_coord[], int bAddedExplicitNeighbor, double *min_sine, int *bAmbiguous )` |
| `e_triple_prod_and_min_abs_sine` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:1884` | external | tree-sitter | `double e_triple_prod_and_min_abs_sine( double at_coord[][3], double *min_sine )` |
| `are_3_vect_in_one_plane` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:1914` | external | tree-sitter | `int are_3_vect_in_one_plane( double at_coord[][3], double min_sine )` |
| `e_are_4at_in_one_plane` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:1926` | external | tree-sitter | `int e_are_4at_in_one_plane( double at_coord[][3], double min_sine )` |
| `e_triple_prod_char` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:1954` | external | tree-sitter | `int e_triple_prod_char( inchi_Atom *at, int at_1, int i_next_at_1, S_CHAR *z_dir1, int at_2, int i_next_at_2, S_CHAR *z_dir2 )` |
| `bCanInpAtomBeAStereoCenter` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:2037` | external | tree-sitter | `int bCanInpAtomBeAStereoCenter( inchi_Atom *at, int cur_at )` |
| `e_bCanInpAtomBeAStereoCenter` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:2128` | external | tree-sitter | `int e_bCanInpAtomBeAStereoCenter( int cur_at, S_CHAR *cAtType )` |
| `e_bCanAtomHaveAStereoBond` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:2164` | external | tree-sitter | `int e_bCanAtomHaveAStereoBond( inchi_Atom *at, int cur_at, S_CHAR *cAtType )` |
| `e_bCanAtomBeMiddleAllene` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:2207` | external | tree-sitter | `int e_bCanAtomBeMiddleAllene( int cur_at, S_CHAR *cAtType )` |
| `e_bCanAtomBeTerminalAllene` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:2223` | external | tree-sitter | `int e_bCanAtomBeTerminalAllene( int cur_at, S_CHAR *cAtType )` |
| `e_FixSb0DParities` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:2237` | external | tree-sitter | `int e_FixSb0DParities( inchi_Atom *at, Stereo0D *pStereo, int chain_length, AT_NUM at_middle, int at_1, int i_next_at_1, S_CHAR z_dir1[], S_CHAR z_dir1NM[], int bOnlyNM1, int bAnomaly1NM, int parity1, int parity1NM, int at_2, int i_next_at_2, S_CHAR z_dir2[], S_CHAR z_dir2NM[], int bOnlyNM2, int bAnomaly2NM, int parity2, int parity2NM )` |
| `e_nNumNonMetalNeigh` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:2487` | external | tree-sitter | `int e_nNumNonMetalNeigh( inchi_Atom *atom, int cur_at, Stereo0D *pStereo, int *i_ord_LastMetal )` |
| `e_half_stereo_bond_parity` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:2586` | external | tree-sitter | `int e_half_stereo_bond_parity( inchi_Atom *at, int cur_at, S_CHAR *z_dir, int *bOnlyNonMetal, int bPointedEdgeStereo, Stereo0D *pStereo )` |
| `e_set_stereo_bonds_parity` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:2899` | external | tree-sitter | `int e_set_stereo_bonds_parity( Stereo0D *pStereo, inchi_Atom *at, int at_1, int bPointedEdgeStereo )` |
| `e_set_stereo_atom_parity` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:3289` | external | tree-sitter | `int e_set_stereo_atom_parity( Stereo0D *pStereo, inchi_Atom *at, int cur_at, int bPointedEdgeStereo )` |
| `e_GetNewStereo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:3682` | external | tree-sitter | `inchi_Stereo0D *e_GetNewStereo( Stereo0D *pStereo )` |
| `set_0D_stereo_parities` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:3709` | external | tree-sitter | `int set_0D_stereo_parities( inchi_Input *pInp, int bPointedEdgeStereo )` |
| `Clear3D2Dstereo` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:3870` | external | tree-sitter | `int Clear3D2Dstereo( inchi_Input *pInp )` |
| `e_GetChiralFlagString` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_0dstereo.c:3891` | external | tree-sitter | `char *e_GetChiralFlagString( int bChiralFlagOn )` |

### `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_io.c`

Parse errors: `4`. Function definitions: `22`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `e_PrintFileName` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_io.c:79` | external | tree-sitter | `void e_PrintFileName( const char *fmt, FILE *output_file, const char *szFname )` |
| `inchi_ios_init` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_io.c:117` | external | tree-sitter | `void inchi_ios_init( INCHI_IOSTREAM* ios, int io_type, FILE *f )` |
| `inchi_ios_flush` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_io.c:140` | external | tree-sitter | `void inchi_ios_flush( INCHI_IOSTREAM* ios )` |
| `inchi_ios_flush2` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_io.c:178` | external | tree-sitter | `void inchi_ios_flush2( INCHI_IOSTREAM* ios, FILE *f2 )` |
| `inchi_ios_close` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_io.c:222` | external | tree-sitter | `void inchi_ios_close( INCHI_IOSTREAM* ios )` |
| `inchi_ios_reset` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_io.c:244` | external | tree-sitter | `void inchi_ios_reset( INCHI_IOSTREAM* ios )` |
| `inchi_ios_str_getc` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_io.c:259` | external | tree-sitter | `int inchi_ios_str_getc( INCHI_IOSTREAM *ios )` |
| `inchi_ios_str_gets` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_io.c:281` | external | tree-sitter | `char *inchi_ios_str_gets( char *szLine, int len, INCHI_IOSTREAM *f )` |
| `inchi_ios_str_getsTab` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_io.c:309` | external | tree-sitter | `char *inchi_ios_str_getsTab( char *szLine, int len, INCHI_IOSTREAM *f )` |
| `inchi_ios_gets` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_io.c:339` | external | tree-sitter | `int inchi_ios_gets( char *szLine, int len, INCHI_IOSTREAM *f, int *bTooLongLine )` |
| `inchi_ios_getsTab` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_io.c:372` | external | tree-sitter | `int inchi_ios_getsTab( char *szLine, int len, INCHI_IOSTREAM *f, int *bTooLongLine )` |
| `inchi_ios_getsTab1` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_io.c:400` | external | tree-sitter | `int inchi_ios_getsTab1( char *szLine, int len, INCHI_IOSTREAM *f, int *bTooLongLine )` |
| `inchi_ios_print` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_io.c:430` | external | tree-sitter | `int inchi_ios_print( INCHI_IOSTREAM * ios, const char* lpszFormat, ... )` |
| `inchi_ios_print_nodisplay` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_io.c:520` | external | tree-sitter | `int inchi_ios_print_nodisplay( INCHI_IOSTREAM * ios, const char* lpszFormat, ... )` |
| `inchi_ios_eprint` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_io.c:600` | external | tree-sitter | `int inchi_ios_eprint( INCHI_IOSTREAM * ios, const char* lpszFormat, ... )` |
| `inchi_fprintf` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_io.c:699` | external | tree-sitter | `int inchi_fprintf( FILE* f, const char* lpszFormat, ... )` |
| `inchi_vfprintf` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_io.c:730` | external | tree-sitter | `int inchi_vfprintf( FILE* f, const char* lpszFormat, va_list argList )` |
| `inchi_fgetsLfTab` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_io.c:766` | external | tree-sitter | `int inchi_fgetsLfTab( char *szLine, int len, FILE *f )` |
| `inchi_fgetsLfTab` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_io.c:798` | external | tree-sitter | `int inchi_fgetsLfTab( char *szLine, int len, FILE *f )` |
| `inchi_fgetsTab` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_io.c:838` | external | tree-sitter | `char *inchi_fgetsTab( char *szLine, int len, FILE *f )` |
| `inchi_fgetsLf` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_io.c:868` | external | tree-sitter | `char* inchi_fgetsLf( char* line, int line_len, FILE* inp )` |
| `GetMaxPrintfLength` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_io.c:905` | external | tree-sitter | `int GetMaxPrintfLength( const char *lpszFormat, va_list argList )` |

### `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_parms.c`

Parse errors: `210`. Function definitions: `17`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `ReadCommandLineParms` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_parms.c:72` | external | tree-sitter | `int ReadCommandLineParms( int argc, const char *argv[], INPUT_PARMS *ip, char *szSdfDataValue, unsigned long *ulDisplTime, int bReleaseVersion, INCHI_IOSTREAM *log_file )` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_parms.c:320` | external | tree-sitter | `else if (!inchi_memicmp( pArg, "DISCONSALT:", 11 ))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_parms.c:325` | external | tree-sitter | `else if (!inchi_memicmp( pArg, "DISCONMETAL:", 12 ))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_parms.c:330` | external | tree-sitter | `else if (!inchi_memicmp( pArg, "RECONMETAL:", 11 ))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_parms.c:335` | external | tree-sitter | `else if (!inchi_memicmp( pArg, "DISCONMETALCHKVAL:", 18 ))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_parms.c:340` | external | tree-sitter | `else if (!inchi_memicmp( pArg, "MOVEPOS:", 8 ))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_parms.c:345` | external | tree-sitter | `else if (!inchi_memicmp( pArg, "MERGESALTTG:", 12 ))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_parms.c:350` | external | tree-sitter | `else if (!inchi_memicmp( pArg, "UNCHARGEDACIDS:", 15 ))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_parms.c:355` | external | tree-sitter | `else if (!inchi_memicmp( pArg, "ACIDTAUT:", 9 ))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_parms.c:969` | external | tree-sitter | `else if (( bVer1Options & 1 ) && INCHI_OPTION_PREFX == argv[i][0] && argv[i][1]) { /*=== parsing stand-alone/library InChI options ===*/ pArg = argv[i] + 1; bRecognizedOption = 2; bVer1Options += 2; /* always on: REQ_MODE_TAUT \| REQ_MODE_ISO \| REQ_MODE_STEREO */ /*--- Input options ---*/ if (!inchi_stricmp( pArg, "STDIO" )) { bNameSuffix = 0; } else if (!inchi_stricmp( pArg, "INPAUX" ))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_parms.c:1681` | external | tree-sitter | `else if (ip->num_paths < MAX_NUM_PATHS)` |
| `PrintInputParms` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_parms.c:2132` | external | tree-sitter | `int PrintInputParms( INCHI_IOSTREAM *log_file, INPUT_PARMS *ip )` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_parms.c:2192` | external | tree-sitter | `else if (bInChI2Struct)` |
| `HelpCommandLineParms` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_parms.c:2642` | external | tree-sitter | `void HelpCommandLineParms(INCHI_IOSTREAM* f)` |
| `OpenFiles` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_parms.c:2830` | external | tree-sitter | `int OpenFiles( FILE **inp_file, FILE **out_file, FILE **log_file, FILE **prb_file, INPUT_PARMS *ip )` |
| `bMatchOnePrefix` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_parms.c:3003` | internal | tree-sitter | `static int bMatchOnePrefix( int len, char *str, int lenPrefix[], char strPrefix[][LEN_VERSIONS], int numPrefix )` |
| `DetectInputINChIFileType` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichi_parms.c:3022` | external | tree-sitter | `int DetectInputINChIFileType( FILE **inp_file, INPUT_PARMS *ip, const char *fmode )` |

### `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichimain.c`

Parse errors: `4`. Function definitions: `6`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `MyHandlerRoutine` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichimain.c:136` | external | tree-sitter | `BOOL WINAPI MyHandlerRoutine( DWORD dwCtrlType /* control signal type */ )` |
| `e_WasInterrupted` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichimain.c:153` | external | tree-sitter | `int e_WasInterrupted( void )` |
| `e_stristr` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichimain.c:173` | internal | tree-sitter | `static char *e_stristr( const char * string1, const char * string2 )` |
| `e_bEnableCmdlineOption` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichimain.c:202` | internal | tree-sitter | `static int e_bEnableCmdlineOption( char *szCmdLine, const char *szOption, int bEnable )` |
| `main` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichimain.c:240` | external | tree-sitter | `int main( int argc, char *argv[] )` |
| `e_MakeOutputHeader` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichimain.c:1089` | external | tree-sitter | `int e_MakeOutputHeader(char* pSdfLabel, char* pSdfValue, long lSdfId, long num_inp, char* pStr1, char* pStr2, int pStr1_size, int pStr2_size)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichimain_a.c`

Parse errors: `6`. Function definitions: `6`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `MyHandlerRoutine` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichimain_a.c:133` | external | tree-sitter | `BOOL WINAPI MyHandlerRoutine( DWORD dwCtrlType ) /* control signal type */` |
| `e_WasInterrupted` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichimain_a.c:146` | external | tree-sitter | `int e_WasInterrupted( void )` |
| `main` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichimain_a.c:172` | external | tree-sitter | `int main( int argc, char *argv[] )` |
| `e_MakeOutputHeader` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichimain_a.c:1016` | external | tree-sitter | `int e_MakeOutputHeader( char *pSdfLabel, char *pSdfValue, long lSdfId, long num_inp, char *pStr1, char *pStr2 )` |
| `e_HelpCommandLineParmsReduced` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichimain_a.c:1044` | external | tree-sitter | `void e_HelpCommandLineParmsReduced( INCHI_IOSTREAM *f )` |
| `e_set_inchi_input_by_extended_inchi_input` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_ichimain_a.c:1156` | external | tree-sitter | `int e_set_inchi_input_by_extended_inchi_input(inchi_Input *pInp, inchi_InputEx *pInpEx)` |

### `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_inchi_atom.c`

Parse errors: `0`. Function definitions: `6`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `e_FreeInchi_Atom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_inchi_atom.c:58` | external | tree-sitter | `void e_FreeInchi_Atom( inchi_Atom **at )` |
| `e_FreeInchi_Stereo0D` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_inchi_atom.c:67` | external | tree-sitter | `void e_FreeInchi_Stereo0D( inchi_Stereo0D **stereo0D )` |
| `e_CreateInchi_Atom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_inchi_atom.c:76` | external | tree-sitter | `inchi_Atom *e_CreateInchi_Atom( int num_atoms )` |
| `e_CreateInchi_Stereo0D` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_inchi_atom.c:82` | external | tree-sitter | `inchi_Stereo0D *e_CreateInchi_Stereo0D( int num_stereo0D )` |
| `e_FreeInchi_Input` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_inchi_atom.c:90` | external | tree-sitter | `void e_FreeInchi_Input( inchi_InputEx *inp_at_data )` |
| `e_RemoveRedundantNeighbors` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_inchi_atom.c:200` | external | tree-sitter | `int e_RemoveRedundantNeighbors( inchi_Input *inp_at_data )` |

### `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_mol2atom.c`

Parse errors: `0`. Function definitions: `5`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `mol_to_inchi_Atom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_mol2atom.c:96` | external | tree-sitter | `inchi_Atom* mol_to_inchi_Atom( MOL_DATA* mol_data, int *num_atoms, int *num_bonds, inchi_Atom* at_inp, int bDoNotAddH, int *err, char *pStrErr )` |
| `mol_to_inchi_Atom_xyz` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_mol2atom.c:331` | external | tree-sitter | `int mol_to_inchi_Atom_xyz( MOL_DATA* mol_data, int num_atoms, inchi_Atom* at, int *err, char *pStrErr )` |
| `GetMolfileNumber` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_mol2atom.c:505` | external | tree-sitter | `long GetMolfileNumber( MOL_HEADER_BLOCK *pHdr )` |
| `MolfileToInchi_Atom` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_mol2atom.c:530` | external | tree-sitter | `int MolfileToInchi_Atom( FILE *inp_molfile, int bDoNotAddH, inchi_Atom **at, int max_num_at, int *num_dimensions, int *num_bonds, const char *pSdfLabel, char *pSdfValue, long *Id, long *lMolfileNumber, INCHI_MODE *pInpAtomFlags, int *err, char *pStrErr )` |
| `e_MolfileToInchi_Input` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_mol2atom.c:659` | external | tree-sitter | `int e_MolfileToInchi_Input( FILE *inp_molfile, inchi_InputEx *orig_at_data, int bMergeAllInputStructures, int bDoNotAddH, int bAllowEmptyStructure, const char *pSdfLabel, char *pSdfValue, long *lSdfId, long *lMolfileNumber, INCHI_MODE *pInpAtomFlags, int *err, char *pStrErr )` |

### `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_readinch.c`

Parse errors: `0`. Function definitions: `0`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|

### `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_readmol.c`

Parse errors: `27`. Function definitions: `17`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `AddMOLfileError` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_readmol.c:131` | external | tree-sitter | `int AddMOLfileError( char *pStrErr, const char *szMsg )` |
| `mol_copy_check_empty` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_readmol.c:175` | external | tree-sitter | `int mol_copy_check_empty( char* dest, char* source, int len, char **first_space )` |
| `mol_read_datum` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_readmol.c:195` | external | tree-sitter | `int mol_read_datum( void* data, int field_len, int data_type, char** line_ptr )` |
| `mol_read_hdr` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_readmol.c:404` | external | tree-sitter | `int mol_read_hdr( MOL_HEADER_BLOCK *hdr, FILE* inp, char *pStrErr )` |
| `RemoveNonPrintable` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_readmol.c:487` | external | tree-sitter | `int RemoveNonPrintable( char *line )` |
| `mol_read_counts_line` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_readmol.c:508` | external | tree-sitter | `int mol_read_counts_line( MOL_CTAB* ctab, FILE *inp, char *pStrErr )` |
| `read_atom_block` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_readmol.c:564` | external | tree-sitter | `int read_atom_block( MOL_CTAB* ctab, FILE *inp, int err, char *pStrErr )` |
| `read_bonds_block` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_readmol.c:709` | external | tree-sitter | `int read_bonds_block( MOL_CTAB* ctab, FILE *inp, int err, char *pStrErr )` |
| `read_stext_block` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_readmol.c:789` | external | tree-sitter | `int read_stext_block( MOL_CTAB* ctab, FILE *inp, int err, char *pStrErr )` |
| `read_properties_block` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_readmol.c:825` | external | tree-sitter | `int read_properties_block( MOL_CTAB* ctab, MOL_HEADER_BLOCK *pHdr, FILE *inp, int err, char *pStrErr )` |
| `delete_mol_data` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_readmol.c:1188` | external | tree-sitter | `MOL_DATA* delete_mol_data( MOL_DATA* mol_data )` |
| `read_mol_file` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_readmol.c:1215` | external | tree-sitter | `MOL_DATA* read_mol_file( FILE* inp, MOL_HEADER_BLOCK *OnlyHeaderBlock, MOL_CTAB *OnlyCtab, int bGetOrigCoord, int *err, char *pStrErr )` |
| `extract_cas_rn` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_readmol.c:1357` | external | tree-sitter | `long extract_cas_rn( char *line )` |
| `identify_sdf_label` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_readmol.c:1377` | external | tree-sitter | `int identify_sdf_label( char* inp_line, const char *pSdfLabel )` |
| `bypass_sdf_data_items` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_readmol.c:1421` | external | tree-sitter | `int bypass_sdf_data_items( FILE* inp, long *cas_reg_no, char* comment, int lcomment, char *name, int lname, int prev_err, const char *pSdfLabel, char *pSdfValue, char *pStrErr )` |
| `read_sdfile_segment` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_readmol.c:1598` | external | tree-sitter | `MOL_DATA* read_sdfile_segment( FILE* inp, MOL_HEADER_BLOCK *OnlyHeaderBlock, MOL_CTAB *OnlyCtab, int bGetOrigCoord, char *pname, int lname, long *Id, const char *pSdfLabel, char *pSdfValue, int *err, char *pStrErr )` |
| `CopyMOLfile` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_readmol.c:1640` | external | tree-sitter | `int CopyMOLfile( FILE *inp_file, long fPtrStart, long fPtrEnd, FILE *prb_file, long lNumb )` |

### `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_readstru.c`

Parse errors: `0`. Function definitions: `3`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `e_ReadStructure` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_readstru.c:78` | external | tree-sitter | `int e_ReadStructure( STRUCT_DATA *sd, INPUT_PARMS *ip, INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM *log_file, INCHI_IOSTREAM *output_file, INCHI_IOSTREAM *prb_file, inchi_InputEx *pInp, long num_inp, int inp_index, int *out_index )` |
| `e_TreatReadTheStructureErrors` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_readstru.c:223` | external | tree-sitter | `int e_TreatReadTheStructureErrors( STRUCT_DATA *sd, INPUT_PARMS *ip, INCHI_IOSTREAM *inp_file, INCHI_IOSTREAM* log_file, INCHI_IOSTREAM* output_file, INCHI_IOSTREAM *prb_file, inchi_InputEx *pInp, long *num_inp )` |
| `e_GetInpStructErrorType` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_readstru.c:287` | external | tree-sitter | `int e_GetInpStructErrorType( INPUT_PARMS *ip, int err, char *pStrErrStruct, int num_inp_atoms )` |

### `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c`

Parse errors: `3`. Function definitions: `29`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `extract_ChargeRadical` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:76` | external | tree-sitter | `int extract_ChargeRadical( char *elname, int *pnRadical, int *pnCharge )` |
| `normalize_name` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:140` | external | tree-sitter | `int normalize_name( char* name )` |
| `e_inchi_malloc` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:176` | external | tree-sitter | `void *e_inchi_malloc( size_t c )` |
| `e_inchi_calloc` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:182` | external | tree-sitter | `void *e_inchi_calloc( size_t c, size_t n )` |
| `e_inchi_free` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:188` | external | tree-sitter | `void e_inchi_free( void *p )` |
| `e_mystrncpy` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:201` | external | tree-sitter | `int e_mystrncpy( char *target, const char *source, unsigned maxlen )` |
| `e_LtrimRtrim` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:224` | external | tree-sitter | `char* e_LtrimRtrim( char *p, int* nLen )` |
| `e_remove_trailing_spaces` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:246` | external | tree-sitter | `void e_remove_trailing_spaces( char* p )` |
| `e_remove_one_lf` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:254` | external | tree-sitter | `void e_remove_one_lf( char* p )` |
| `e_is_in_the_slist` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:265` | external | tree-sitter | `S_SHORT *e_is_in_the_slist( S_SHORT *pathAtom, S_SHORT nNextAtom, int nPathLen )` |
| `e_is_element_a_metal` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:272` | external | tree-sitter | `int e_is_element_a_metal( char szEl[] )` |
| `InchiClock` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:301` | internal | tree-sitter | `static clock_t InchiClock( void )` |
| `InchiClock` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:312` | internal | tree-sitter | `static clock_t InchiClock( void )` |
| `FillMaxMinClock` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:340` | internal | tree-sitter | `static void FillMaxMinClock( void )` |
| `InchiTimeMsecDiff` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:358` | external | tree-sitter | `long InchiTimeMsecDiff( e_inchiTime *TickEnd, e_inchiTime *TickStart )` |
| `InchiTimeElapsed` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:429` | external | tree-sitter | `long InchiTimeElapsed( e_inchiTime *TickStart )` |
| `InchiTimeAddMsec` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:438` | external | tree-sitter | `void InchiTimeAddMsec( e_inchiTime *TickEnd, unsigned long nNumMsec )` |
| `bInchiTimeIsOver` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:459` | external | tree-sitter | `int bInchiTimeIsOver( e_inchiTime *TickStart )` |
| `e_inchiTimeGet` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:529` | external | tree-sitter | `void e_inchiTimeGet( e_inchiTime *TickEnd )` |
| `InchiTimeMsecDiff` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:540` | external | tree-sitter | `long InchiTimeMsecDiff( e_inchiTime *TickEnd, e_inchiTime *TickStart )` |
| `InchiTimeElapsed` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:562` | external | tree-sitter | `long InchiTimeElapsed( e_inchiTime *TickStart )` |
| `InchiTimeAddMsec` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:571` | external | tree-sitter | `void InchiTimeAddMsec( e_inchiTime *TickEnd, unsigned long nNumMsec )` |
| `bInchiTimeIsOver` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:582` | external | tree-sitter | `int bInchiTimeIsOver( e_inchiTime *TickEnd )` |
| `e_inchiTimeMsecDiff` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:603` | external | tree-sitter | `long e_inchiTimeMsecDiff( e_inchiTime *TickEnd, e_inchiTime *TickStart )` |
| `e_inchiTimeElapsed` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:674` | external | tree-sitter | `long e_inchiTimeElapsed( e_inchiTime *TickStart )` |
| `e_inchiTimeAddMsec` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:683` | external | tree-sitter | `void e_inchiTimeAddMsec( e_inchiTime *TickEnd, unsigned long nNumMsec )` |
| `e_inchiTimeGet` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:706` | external | tree-sitter | `void e_inchiTimeGet( e_inchiTime *TickEnd )` |
| `e_inchi_memicmp` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:717` | external | tree-sitter | `int e_inchi_memicmp( const void * p1, const void * p2, size_t length )` |
| `e_inchi_stricmp` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/inchi_main/src/e_util.c:738` | external | tree-sitter | `int e_inchi_stricmp( const char *s1, const char *s2 )` |

### `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/mol2inchi.c`

Parse errors: `0`. Function definitions: `14`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `main` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/mol2inchi.c:77` | external | tree-sitter | `int main( int argc, char *argv[] )` |
| `m2i_WorkPool_wait_and_print_all` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/mol2inchi.c:369` | external | tree-sitter | `int m2i_WorkPool_wait_and_print_all( m2i_WorkPool *pool, THREAD_PTR *pthreads )` |
| `m2i_Worker_run` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/mol2inchi.c:415` | external | tree-sitter | `M2I_THREADFUNC m2i_Worker_run( void *arg ) /* djb-rwth: ignoring LLVM warning */` |
| `print_help` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/mol2inchi.c:497` | external | tree-sitter | `void print_help( void )` |
| `m2i_WorkDetails_init` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/mol2inchi.c:520` | external | tree-sitter | `int m2i_WorkDetails_init( m2i_WorkDetails *wd, int argc, char *argv[] )` |
| `m2i_WorkPool_init` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/mol2inchi.c:637` | external | tree-sitter | `int m2i_WorkPool_init( m2i_WorkPool *pool, m2i_WorkDetails *wd )` |
| `m2i_WorkPool_close` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/mol2inchi.c:667` | external | tree-sitter | `void m2i_WorkPool_close( m2i_WorkPool *pool )` |
| `m2i_WorkPool_reset` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/mol2inchi.c:679` | external | tree-sitter | `void m2i_WorkPool_reset( m2i_WorkPool *pool ) /* do not touch already allocd pointers to workers */` |
| `m2i_Worker_init` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/mol2inchi.c:695` | external | tree-sitter | `int m2i_Worker_init( m2i_Worker *w, m2i_WorkDetails *wd )` |
| `m2i_Worker_reset` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/mol2inchi.c:723` | external | tree-sitter | `void m2i_Worker_reset( m2i_Worker *w ) /* do not touch already allocd pointers to tasks */` |
| `m2i_Worker_close` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/mol2inchi.c:735` | external | tree-sitter | `void m2i_Worker_close( m2i_Worker *w )` |
| `m2i_Task_reset` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/mol2inchi.c:752` | external | tree-sitter | `void m2i_Task_reset( m2i_Task *task )` |
| `m2i_Task_close` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/mol2inchi.c:761` | external | tree-sitter | `void m2i_Task_close( m2i_Task *task )` |
| `m2i_Task_init` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/mol2inchi.c:770` | external | tree-sitter | `int m2i_Task_init( m2i_Task *task, m2i_WorkDetails *wd )` |

### `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/moreutil.c`

Parse errors: `4`. Function definitions: `16`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `print_time` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/moreutil.c:70` | external | tree-sitter | `void print_time( void )` |
| `own_stricmp` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/moreutil.c:84` | external | tree-sitter | `int own_stricmp( const char *s1, const char *s2 )` |
| `own_memicmp` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/moreutil.c:110` | external | tree-sitter | `int own_memicmp( const void *p1, const void *p2, size_t length )` |
| `get_substr_in_between` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/moreutil.c:142` | external | tree-sitter | `char* get_substr_in_between( char *s, char *pat1, char *pat2, char *buf, size_t max_symbols, size_t *copied )` |
| `fgets_lf` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/moreutil.c:192` | external | tree-sitter | `char* fgets_lf( char* line, int line_len, FILE *f )` |
| `GrowStr_init` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/moreutil.c:230` | external | tree-sitter | `int GrowStr_init( GrowStr *buf, int start_size, int incr_size )` |
| `GrowStr_copy` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/moreutil.c:258` | external | tree-sitter | `int GrowStr_copy( GrowStr *dest, GrowStr *src )` |
| `GrowStr_reset` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/moreutil.c:282` | external | tree-sitter | `int GrowStr_reset( GrowStr *buf )` |
| `GrowStr_close` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/moreutil.c:303` | external | tree-sitter | `void GrowStr_close( GrowStr *buf )` |
| `GrowStr_update` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/moreutil.c:320` | external | tree-sitter | `int GrowStr_update( GrowStr *buf, int requested_addition )` |
| `GrowStr_putc` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/moreutil.c:376` | external | tree-sitter | `int GrowStr_putc( GrowStr *buf, char c )` |
| `GrowStr_printf` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/moreutil.c:400` | external | tree-sitter | `int GrowStr_printf( GrowStr *buf, const char* lpszFormat, ... )` |
| `GrowStr_readline` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/moreutil.c:438` | external | tree-sitter | `int GrowStr_readline( GrowStr *buf, FILE *f )` |
| `GetMaxPrintfLength` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/moreutil.c:471` | external | tree-sitter | `int GetMaxPrintfLength( const char *lpszFormat, va_list argList )` |
| `get_next_molfile_as_growing_str_buffer` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/moreutil.c:718` | external | tree-sitter | `int get_next_molfile_as_growing_str_buffer( FILE *f, GrowStr *buf, int *at_eof )` |
| `get_msec_timer` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/mol2inchi/src/moreutil.c:767` | external | tree-sitter | `unsigned int get_msec_timer( void )` |

### `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/moreutil.c`

Parse errors: `0`. Function definitions: `7`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `print_time` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/moreutil.c:61` | external | tree-sitter | `void print_time( void )` |
| `own_stricmp` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/moreutil.c:74` | external | tree-sitter | `int own_stricmp( const char *s1, const char *s2 )` |
| `own_memicmp` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/moreutil.c:101` | external | tree-sitter | `int own_memicmp( const void * p1, const void * p2, size_t length )` |
| `get_substr_in_between` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/moreutil.c:133` | external | tree-sitter | `char* get_substr_in_between( char *s, char *pat1, char *pat2, char *buf, size_t max_symbols, size_t *copied )` |
| `fgets_lf` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/moreutil.c:186` | external | tree-sitter | `char* fgets_lf( char* line, int line_len, FILE *f )` |
| `get_next_molfile_as_text` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/moreutil.c:221` | external | tree-sitter | `int get_next_molfile_as_text( FILE *f, char *buf, size_t buflen )` |
| `is_empty_text` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/moreutil.c:269` | external | tree-sitter | `int is_empty_text( char *buf )` |

### `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c`

Parse errors: `59`. Function definitions: `35`.

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `CheckStatus` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:65` | internal | tree-sitter | `static int CheckStatus( IXA_STATUS_HANDLE hStatus, long nrecord )` |
| `OptionCompare` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:99` | internal | tree-sitter | `static int OptionCompare( const char *pArg, const char *pOption )` |
| `ReadOptions` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:132` | internal | tree-sitter | `static int ReadOptions( int argc, const char *argv[], IXA_BOOL *pKey, IXA_BOOL *pRoundTrip, IXA_BOOL *pGenerateAuxinfo, IXA_BOOL *pVerbose, char *pOptions, IXA_STATUS_HANDLE hStatus, IXA_INCHIBUILDER_HANDLE hBuilder )` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:218` | external | tree-sitter | `else if (OptionCompare(argv[index], "PT_22_00"))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:224` | external | tree-sitter | `else if (OptionCompare(argv[index], "PT_16_00"))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:230` | external | tree-sitter | `else if (OptionCompare(argv[index], "PT_06_00"))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:236` | external | tree-sitter | `else if (OptionCompare(argv[index], "PT_39_00"))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:242` | external | tree-sitter | `else if (OptionCompare(argv[index], "PT_13_00"))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:248` | external | tree-sitter | `else if (OptionCompare(argv[index], "PT_18_00"))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:253` | external | tree-sitter | `else if (OptionCompare( argv[index], "AuxNone" ))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:258` | external | tree-sitter | `else if (OptionCompare( argv[index], "WarnOnEmptyStructure" ))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:262` | external | tree-sitter | `else if (OptionCompare( argv[index], "LargeMolecules" ))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:267` | external | tree-sitter | `else if (OptionCompare( argv[index], "DoDrv" ))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:271` | external | tree-sitter | `else if (OptionCompare( argv[index], "DoDrvReport" ))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:275` | external | tree-sitter | `else if (OptionCompare( argv[index], "DoR2C" ))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:279` | external | tree-sitter | `else if (OptionCompare( argv[index], "DoneOnly" ))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:283` | external | tree-sitter | `else if (OptionCompare( argv[index], "OnlyRecSalt" ))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:287` | external | tree-sitter | `else if (OptionCompare( argv[index], "OnlyExact" ))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:291` | external | tree-sitter | `else if (OptionCompare( argv[index], "OnlyRecMet" ))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:295` | external | tree-sitter | `else if (OptionCompare(argv[index], "Polymers105+"))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:299` | external | tree-sitter | `else if (OptionCompare(argv[index], "FilterSS"))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:303` | external | tree-sitter | `else if (OptionCompare(argv[index], "InvFilterSS"))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:309` | external | tree-sitter | `else if (OptionCompare( argv[index], "SaveOpt" ))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:314` | external | tree-sitter | `else if (OptionCompare(argv[index], "OutErrInChI"))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:318` | external | tree-sitter | `else if (OptionCompare(argv[index], "LooseTSACheck"))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:322` | external | tree-sitter | `else if (OptionCompare(argv[index], "NoWarnings"))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:326` | external | tree-sitter | `else if (OptionCompare(argv[index], "Polymers"))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:330` | external | tree-sitter | `else if (OptionCompare(argv[index], "Polymers105"))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:334` | external | tree-sitter | `else if (OptionCompare(argv[index], "NPZz"))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:338` | external | tree-sitter | `else if (OptionCompare(argv[index], "SAtZz"))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:342` | external | tree-sitter | `else if (OptionCompare(argv[index], "FoldSRU"))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:346` | external | tree-sitter | `else if (OptionCompare(argv[index], "FoldCRU"))` |
| `if` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:350` | external | tree-sitter | `else if (OptionCompare(argv[index], "NoEdits"))` |
| `main` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:397` | external | tree-sitter | `int main( int argc, const char *argv[] ) /* djb-rwth: main function needs to be int as LLVM/Clang otherwise reports error */` |
| `print_help` | `third_party/InChI/INCHI-1-SRC/INCHI_API/demos/test_ixa/src/test_ixa.c:810` | external | tree-sitter | `void print_help( void )` |


## Production Header-Defined Functions

| Function | Source | Linkage | Extraction | Signature |
|---|---|---|---|---|
| `stbsp_set_separators` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:274` | external | tree-sitter+gcc-aux | `extern void stbsp_set_separators (char pcomma, char pperiod)` |
| `stbsp__lead_sign` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:294` | internal | tree-sitter+gcc-aux | `static void stbsp__lead_sign (unsigned int fl, char *sign)` |
| `stbsp__strlen_limited` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:309` | internal | tree-sitter+gcc-aux | `static unsigned int stbsp__strlen_limited (const char *s, unsigned int limit)` |
| `stbsp_vsprintfcb` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:349` | external | tree-sitter+gcc-aux | `extern int stbsp_vsprintfcb (STBSP_SPRINTFCB (*callback), void *user, char *buf, const char *fmt, __va_list_tag *va)` |
| `stbsp_sprintf` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:1376` | external | tree-sitter+gcc-aux | `extern int stbsp_sprintf (char *buf, const char *fmt, ...)` |
| `stbsp__clamp_callback` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:1393` | internal | tree-sitter+gcc-aux | `static char *stbsp__clamp_callback (const char *buf, void *user, int len)` |
| `stbsp__count_clamp_callback` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:1421` | internal | tree-sitter+gcc-aux | `static char *stbsp__count_clamp_callback (const char *buf, void *user, int len)` |
| `stbsp_vsnprintf` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:1430` | external | tree-sitter+gcc-aux | `extern int stbsp_vsnprintf (char *buf, int count, const char *fmt, __va_list_tag *va)` |
| `stbsp_snprintf` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:1460` | external | tree-sitter+gcc-aux | `extern int stbsp_snprintf (char *buf, int count, const char *fmt, ...)` |
| `stbsp_vsprintf` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:1472` | external | tree-sitter+gcc-aux | `extern int stbsp_vsprintf (char *buf, const char *fmt, __va_list_tag *va)` |
| `stbsp__real_to_parts` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:1491` | internal | tree-sitter+gcc-aux | `static int stbsp__real_to_parts (long long int *bits, int *expo, double value)` |
| `stbsp__raise_to_power10` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:1638` | internal | tree-sitter+gcc-aux | `static void stbsp__raise_to_power10 (double *ohi, double *olo, double d, int power)` |
| `stbsp__real_to_str` | `third_party/InChI/INCHI-1-SRC/INCHI_BASE/src/stb_sprintf.h:1705` | internal | tree-sitter+gcc-aux | `static int stbsp__real_to_str (const char **start, unsigned int *len, char *out, int *decimal_pos, double value, unsigned int frac_digits)` |

## Regeneration

Run from the repository root:

```bash
uv run --no-project --with tree-sitter --with tree-sitter-c python dev/generate_inchi_function_inventory.py
```
