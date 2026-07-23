# RDKit InChI Adapter Source Inventory

## Audit Boundary

This report inventories the pinned RDKit `2026.03.1` adapter under
`third_party/rdkit/External/INCHI-API`. It is not an inventory of the official
IUPAC InChI engine. The adapter translates between RDKit molecule state and
the official C API; it cannot generate or parse InChI without that engine.

The RDKit version is recorded by `third_party/rdkit/CMakeLists.txt` as year
`2026`, month `03`, revision `1`, and by `third_party/rdkit/ReleaseNotes.md` as
`Release_2026.03.1`. The vendored RDKit license is
`third_party/rdkit/license.txt`.

## Primary Artifacts

| Path | Lines | SHA-256 | Role |
|---|---:|---|---|
| `inchi.cpp` | 2178 | `104c1ee0c7978f92077c30d5f7a4566af791aa254bbf8d53a4a4bb4c590bad3f` | C++ graph adapter and public API implementation |
| `inchi.h` | 115 | `27b4eaef1714869c42dbf2998807018a03389b4f9ce40438e843248ebfc3614e` | Public declarations and inline `MolToInchiKey` |
| `test.cpp` | 988 | `d3faba4a81f82c03647c275f0c0978cd7887e883b0fd9a3ca3e473252d1a9b01` | Native adapter fixtures and concurrency tests |
| `CMakeLists.txt` | 115 | `bcec4e6642d4cb2ec5455792a2ab9f97c78649e50d00582ad0da4338f17f6a74` | Adapter and official-engine build boundary |
| `Wrap/pyInchi.cpp` |  | `550894c79f6596d35772b02884dec9c5570589339d9c897b73fe489ef31a2bad` | Python tuple-return wrapper |
| `python/inchi.py` |  | `6276b5a3dba0d4cc4fe9830dba7cc066d9714ac2decff2037d1165c630a2408d` | Python convenience/error policy and AuxInfo handling |
| `python/noinchi.py` |  | `b51d80a81037d93d0f6caef0298d47db2269afb58862aa3d7e8a1b7ee35c319d` | Explicit unavailable sentinel |
| `download-inchi.sh` |  | `092486654b55a3e45309d99e6a112b7373c7a30703bfe63397a7ea7b62332e60` | Obsolete 1.05 downloader, not a source-version authority |
| `README` |  | `17b9c1d6215ec03915f1d70ab34f830aec2930b420945debe5767f54effdd3bf` | Older manual engine acquisition instructions |
| `third_party/rdkit/license.txt` |  | `daeb8d194502cbcf34c05c39541a0d02be65bc9bada5b891c1974cd24e9fca30` | RDKit redistribution license |

The three native molecule fixtures are `test_data/github296.mol`,
`test_data/github3.mol`, and `test_data/github8_extra.mol`, with SHA-256 values
`6a56e7d17b94a7a008a76490eb4accb6f65d5a16b3f735fe5911ff29e25c34e9`,
`65850917e63dc83265796148ebecb280167a7223358c7dbf2dcf882497db3803`,
and `153ed38be62fdd579371d39b68030aa55066f3720868f6ca6cf9ea995a2697e0`.

## Adapter Function Inventory

Every row below requires its own named `Port` step in the regenerated plan.
Anonymous lambdas inside `InchiToMol` remain part of that function's copied
source body and test surface; they are not independent source functions.

| Function | Source | Linkage | Responsibility |
|---|---|---|---|
| `assignBondDirs` | `inchi.cpp:91` | anonymous namespace | Propagate same/opposite double-bond direction constraints with deterministic set and queue order |
| `findAlternatingBonds` | `inchi.cpp:195` | anonymous namespace | Recursive shortest alternating-bond path search with shared visited/path state |
| `getNumDoubleBondedNegativelyChargedNeighboringSi` | `inchi.cpp:306` | anonymous namespace | Count double-bonded, negatively charged silicon neighbors |
| `_Valence4NCleanUp1` | `inchi.cpp:324` | anonymous namespace | N(-1), valence-4 ring cleanup via temporary tin substitution and substructure matching |
| `_Valence4NCleanUp2` | `inchi.cpp:376` | anonymous namespace | N(-1), valence-4 alternating-bond cleanup |
| `_Valence5NCleanUp1` | `inchi.cpp:393` | anonymous namespace | Neutral valence-5 N to N(+1) alternating-path cleanup |
| `_Valence5NCleanUp2` | `inchi.cpp:417` | anonymous namespace | Triple/single path cleanup to N(-1) |
| `_Valence5NCleanUp3` | `inchi.cpp:443` | anonymous namespace | Direct/alternating neutral-N cleanup branches |
| `_Valence5NCleanUp4` | `inchi.cpp:477` | anonymous namespace | Two negatively charged silicon-neighbor cleanup |
| `_Valence5NCleanUp5` | `inchi.cpp:544` | anonymous namespace | Parameterized O/S/F/Cl alternating-path cleanup |
| `_Valence5NCleanUp6` | `inchi.cpp:625` | anonymous namespace | Source-specific valence-5 N ring cleanup |
| `_Valence5NCleanUp7` | `inchi.cpp:679` | anonymous namespace | Source-specific N/O alternating-ring cleanup |
| `_Valence5NCleanUp8` | `inchi.cpp:750` | anonymous namespace | Source-specific valence-5 N ring cleanup |
| `_Valence5NCleanUp9` | `inchi.cpp:804` | anonymous namespace | Source-specific valence-5 N ring cleanup |
| `_Valence5NCleanUpA` | `inchi.cpp:855` | anonymous namespace | Tin-substitution/substructure cleanup branch |
| `_Valence5NCleanUpB` | `inchi.cpp:916` | anonymous namespace | Last-resort carbon alternating-path cleanup |
| `_Valence7SCleanUp1` | `inchi.cpp:941` | anonymous namespace | Valence-7 sulfur cleanup branch 1 |
| `_Valence7SCleanUp2` | `inchi.cpp:989` | anonymous namespace | Valence-7 sulfur N-path cleanup |
| `_Valence7SCleanUp3` | `inchi.cpp:1021` | anonymous namespace | Valence-7 sulfur alternating N-path cleanup |
| `_Valence8SCleanUp1` | `inchi.cpp:1044` | anonymous namespace | Valence-8 sulfur alternating-path cleanup |
| `_Valence8ClCleanUp1` | `inchi.cpp:1079` | anonymous namespace | Valence-8 chlorine cleanup |
| `_Valence5ClCleanUp1` | `inchi.cpp:1114` | anonymous namespace | Valence-5 chlorine oxygen-path cleanup |
| `_Valence3ClCleanUp1` | `inchi.cpp:1134` | anonymous namespace | Valence-3 chlorine sulfur-path cleanup |
| `cleanUp` | `inchi.cpp:1151` | anonymous namespace | Ordered post-InChI cleanup dispatcher, including bromine/selenium special case |
| `InchiToMol` | `inchi.cpp:1254` | public | Convert `inchi_OutputStruct` to RDKit graph, isotopes, H, bond/stereo state, sanitize/remove-H, and diagnostics |
| `fixOptionSymbol` | `inchi.cpp:1674` | public translation unit | Convert option prefixes for Windows versus non-Windows |
| `rCleanUp` | `inchi.cpp:1694` | public translation unit | Reverse RDKit perchlorate cleanup before engine input |
| `MolToInchi` | `inchi.cpp:1747` | public | Clone/kekulize graph, build `inchi_Input`, map coordinates/isotopes/H/radicals/stereo, and collect output diagnostics |
| `MolBlockToInchi` | `inchi.cpp:2106` | public | Pass original MolBlock text and normalized options to the official parser/engine |
| `InchiToInchiKey` | `inchi.cpp:2145` | public | Map official InChIKey return codes to key or logged empty result |
| `getInchiVersion` | `inchi.cpp:2177` | public | Return official engine `CURRENT_VER` |
| `MolToInchiKey` | `inchi.h:107` | public inline | Compose `MolToInchi` and `InchiToInchiKey` |

## Official Engine Dependencies

`inchi.cpp` includes `inchi_api.h` and `bcf_s.h`. Its direct official API
surface is:

- `GetStructFromINCHI` and `FreeStructFromINCHI`;
- `GetINCHI` and `FreeINCHI`;
- `MakeINCHIFromMolfileText`;
- `GetINCHIKeyFromINCHI`;
- `CURRENT_VER`.

It also depends on exact definitions and numeric mappings for `inchi_Atom`,
`inchi_Stereo0D`, `inchi_InputINCHI`, `inchi_OutputStruct`, `inchi_Input`,
`inchi_Output`, return codes, bond types, bond stereo, stereo types, parity,
`ISOTOPIC_SHIFT_FLAG`, `MAXVAL`, and `NO_ATOM`.

`External/INCHI-API/src` is absent. Consequently none of the direct API calls,
types, constants, allocation/release behavior, canonicalization, parsing,
printing, or InChIKey hashing is available in this RDKit snapshot. The adapter
alone is not executable and must not be used as an engine substitute.

The CMake file contains a production source list for a separately acquired
official engine, but does not pin a release, commit, checksum, or license. The
legacy downloader says 1.05 while its opening comment says 1.04; neither is
acceptable provenance for the selected full source port.

## RDKit State Dependencies

The adapter relies on behavior outside this directory that must be mapped to
toolkit-neutral graph operations before the adapter can be reused:

- periodic-table symbol, mass, and isotope conversions;
- atom explicit/implicit valence, total H, formal charge, radical, aromatic,
  and chiral state;
- bond direction, bond stereo, stereo-atom ordering, and conformer access;
- graph cloning, atom/bond insertion, adjacency order, and property caches;
- kekulization, sanitization, hydrogen removal, and stereo assignment;
- CIP ranking and perturbation order;
- substructure matching and SMILES query construction used by cleanup code;
- deterministic set, queue, stack, neighbor, atom, and bond iteration order;
- RDKit logging and sanitize exceptions.

These dependencies may be represented by owned neutral data and adapter-local
helpers, but their observable order and failure behavior cannot be replaced by
heuristics.

## Existing Test Surface

`test.cpp` provides roundtrip, InChIKey, MolBlock, multithreading, charge,
radical, isotope/H, stereo, cleanup, sanitize-error, and regression coverage.
Its named tests are `testMultiThread`, `testMultiThread2`,
`testGithubIssue3`, `testGithubIssue8`, `testGithubIssue40`,
`testGithubIssue67`, `testGithubIssue68`, `testGithubIssue296`,
`testGithubIssue437`, `testGithubIssue562`, `testGithubIssue614`,
`testGithubIssue1572`, `testMolBlockToInchi`, `testGithub3365`,
`testGithub3645`, `test_clean_up_on_kekulization_error`,
`testGithub6172`, `testGithub5311`, `testGithub8123`, and
`testGithub8239`.

The Python wrapper additionally fixes tuple shapes and defaults:

- `InchiToMol(inchi, sanitize=true, removeHs=true)` returns molecule or null,
  return code, message, and log;
- `MolToInchi` and `MolBlockToInchi` return identifier, return code, message,
  log, and AuxInfo;
- `MolToInchiKey`, `InchiToInchiKey`, and `GetInchiVersion` return strings;
- `python/inchi.py` adds warning/error policy and AuxInfo atom-order/coordinate
  reconstruction that is separate from `inchi.cpp` and must be inventoried if
  selected for COSMolKit Python parity.

## Source-Declared Gaps

The adapter source itself declares unresolved behavior for allene/extended
double-bond stereochemistry, advanced options such as FixedH, broken input
molecules, metals on InChI read, large-ring stereo without coordinates, and
some radical restoration. Existing code logs and ignores several illegal or
unrecognized radical, bond, and stereo values. These are source behaviors to
reproduce for RDKit adapter parity; they are not permission to invent fallback
chemistry or claim broader official-engine support.

## Required Closure

1. Vendor and pin the approved official InChI engine before any production API
   is exposed.
2. Generate an independent function inventory for the complete selected
   production engine call graph.
3. Give all 32 adapter functions above individual source-backed `Port` steps.
4. Preserve original C++ bodies with `RDKit` markers and official C bodies
   separately with `INCHI` markers.
5. Build exact focused tests before aggregated RDKit parity tests.
6. Keep InChI explicitly unsupported until official engine, adapter, Rust core,
   and Python exposed branches all close without heuristic fallback.
