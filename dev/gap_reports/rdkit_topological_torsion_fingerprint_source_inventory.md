# RDKit Topological Torsion Fingerprint Source Inventory

Status: completed pre-port inventory retained as an audit baseline. The tables
below describe the tree before the port and are intentionally not rewritten as
a current implementation ledger. Current implementation status, exact focused
and maintained-corpus gates, and the complete ChEMBL 37 result are recorded in
[`rdkit_topological_torsion_fingerprint_full_port_validation.md`](rdkit_topological_torsion_fingerprint_full_port_validation.md).

## Audit Contract

This report is the concrete artifact for Step 3 of
`dev/plans/rdkit_topological_torsion_fingerprint_full_port_plan.md`.
It compares the pinned RDKit Topological Torsion production and wrapper call
paths with the current COSMolKit tree. It does not implement any behavior.

The audit follows `dev/source_reproduction_protocol.md`: every behavior-bearing
source function listed here must eventually be copied verbatim into its
corresponding Rust function as two-axis source-marker comments. Dispatch,
ownership, and wrapper lines do not substitute for the helper bodies they call.

## Source Pin

- RDKit distribution version: `2026.3.1`.
- RDKit source version: `2026.03.1`.
- Canonical revision used by this plan and by the verified source URLs:
  `351f8f378f8ad6bbd517980c38896e66bf907af8`.
- Repository metadata: `testdata/reference/rdkit.json` currently records
  `351f8f378f8ad6bbd517980c38896e66bf907af8c`, a 41-character value. GitHub's
  raw endpoint resolves it to the same content as the canonical 40-character
  revision, but the metadata is not a canonical Git SHA-1 and must be
  normalized before final source-version closure.
- License: RDKit BSD-3-Clause.
- Source root: `https://github.com/rdkit/rdkit/tree/351f8f378f8ad6bbd517980c38896e66bf907af8`.

Changing the source revision invalidates the line ranges, source comments, and
oracle outputs in this report.

## Public Calling Paths

### Modern generator path

```text
Python GetTopologicalTorsionGenerator / Rust factory
└─ TopologicalTorsion::getTopologicalTorsionGenerator<uint64_t>
   ├─ TopologicalTorsionArguments
   ├─ TopologicalTorsionEnvGenerator<uint64_t>
   └─ default AtomPairAtomInvGenerator(includeChirality, true)
      └─ FingerprintGenerator<uint64_t>
         ├─ getSparseCountFingerprint
         ├─ getSparseFingerprint
         ├─ getCountFingerprint
         ├─ getFingerprint
         └─ four ordered multi-molecule variants
            └─ getFingerprintHelper
               ├─ optional private ROMol copy + MolOps::assignStereochemistry
               ├─ custom or generated atom invariants
               │  └─ AtomPairAtomInvGenerator::getAtomInvariants
               │     └─ AtomPairs::getAtomCode
               │        └─ Atom::getNumPiElectrons
               └─ TopologicalTorsionEnvGenerator::getEnvironments
                  ├─ findAllPathsOfLengthN
                  │  └─ findAllPathsOfLengthsMtoN
                  │     ├─ pathFinderHelper
                  │     │  └─ extendPaths
                  │     ├─ optional MolOps::getDistanceMat
                  │     │  └─ FloydWarshall
                  │     └─ bond-set deduplication
                  ├─ fromAtoms/ignoreAtoms/cycle filtering
                  ├─ endpoint/internal atom-code correction
                  ├─ getTopologicalTorsionCode for sparse ids
                  └─ getTopologicalTorsionHash for folded ids
                     └─ boost::hash_combine-compatible shared hash
```

The pinned `getFingerprintHelper()` uses two molecule pointers deliberately:
it builds `lmol` as a private stereochemistry-prepared copy when necessary,
but lines 346-359 obtain atom and bond invariants from the original `mol`;
only environment generation at lines 364-367 receives `*lmol`. The Rust port
must reproduce and test this source-object selection rather than assuming that
the private copy supplies atom invariants.

### Legacy paths

```text
getTopologicalTorsionFingerprint
└─ direct invariant/path/code loop
   └─ compatibility sparse length = (1 << totalWidth) - 1

getHashedTopologicalTorsionFingerprint
└─ TorsionFpCalc
   └─ modern generator getCountFingerprint

getHashedTopologicalTorsionFingerprintAsBitVect
└─ TorsionFpCalc(nBits / nBitsPerEntry)
   └─ legacy threshold expansion into nBits
```

The legacy un-hashed function is not replaced by the modern result wrapper
because its signed sparse-vector length compatibility behavior is public. The
legacy hashed-count path already delegates upstream to the modern generator and
must do the same in COSMolKit. The legacy bit-vector threshold expansion remains
an adapter over that shared count result.

### Python paths

```text
rdFingerprintGenerator.GetTopologicalTorsionGenerator
└─ TopologicalTorsionWrapper.cpp
   └─ shared FingerprintGeneratorWrapper scalar/bulk/options/JSON surface

rdMolDescriptors.GetTopologicalTorsionFingerprint*
└─ rdMolDescriptors.cpp legacy wrappers
   └─ AtomPairs.cpp legacy entry points

Chem.AtomPairs.Torsions helpers
├─ pyScorePath -> Utils.GetAtomCode -> GetAtomPairAtomCode
├─ ExplainPathScore -> Utils.ExplainAtomCode
└─ GetTopologicalTorsionFingerprintAsIds -> legacy sparse fingerprint
```

COSMolKit-native public names remain snake case, but every exposed behavior and
default is tied to the same Rust atom-code, torsion-code, environment, and
vector core.

## Production Source Inventory

| Source range at the pinned revision | Behavior reached by Topological Torsion | Port disposition |
|---|---|---|
| `Code/GraphMol/Fingerprints/FingerprintUtil.h:28-41` | Atom type table and bit-width constants | One shared typed constant set; no Python-only duplicate. |
| `Code/GraphMol/Atom.cpp:934-953` | `Atom::getNumPiElectrons()` | Add one shared atom helper with full source body. |
| `Code/GraphMol/Fingerprints/FingerprintUtil.cpp:45-97` | `AtomPairs::getAtomCode()` | Add one shared atom-code kernel. |
| `Code/GraphMol/Fingerprints/AtomPairGenerator.h:22-47` | invariant-generator contract/defaults/clone | One reusable typed invariant-generator unit. |
| `Code/GraphMol/Fingerprints/AtomPairGenerator.cpp:26-66` | constructor, `getAtomInvariants`, info, JSON, clone | Port complete unit; correction flag is required by torsions. |
| `Code/GraphMol/Fingerprints/FingerprintUtil.cpp:109-138` | canonical packed 64-bit torsion code | Add one shared un-hashed code function. |
| `Code/GraphMol/Fingerprints/FingerprintUtil.cpp:140-167` | canonical 32-bit torsion hash | Add function delegating to the existing shared hash combine. |
| `Code/GraphMol/Subgraphs/Subgraphs.h:39-48,107-130` | ordered path types and call defaults | Preserve ordering and argument semantics internally. |
| `Code/GraphMol/Subgraphs/Subgraphs.cpp:189-242` | `extendPaths()` | Port exact atom-path extension and final ring closure. |
| `Code/GraphMol/Subgraphs/Subgraphs.cpp:244-282` | `pathFinderHelper()` | Port exact seed and length-map ordering. |
| `Code/GraphMol/Subgraphs/Subgraphs.cpp:442-549` | `findAllPathsOfLengthsMtoN()` | Port exact adjacency, optional distances, conversion, deduplication. |
| `Code/GraphMol/Subgraphs/Subgraphs.cpp:550-555` | `findAllPathsOfLengthN()` | Single adapter over the range function. |
| `Code/GraphMol/Matrices.cpp:39-100` | `FloydWarshall()` | Move the existing Rust port to one shared private module. |
| `Code/GraphMol/Matrices.cpp:166-230` | `MolOps::getDistanceMat()` topology branch | Reuse shared matrix kernel for `onlyShortestPaths`. |
| `Code/GraphMol/Chirality.cpp:2889-2905` | default `MolOps::assignStereochemistry()` call | Reuse shared stereo machinery only on the private copy; audit exact observable boundary. |
| `Code/GraphMol/MolOps.h:1185-1187` | default stereo-call flags | Preserve `cleanIt=false`, `force=false`, `flagPossibleStereoCenters=false`. |
| `Code/GraphMol/Fingerprints/FingerprintGenerator.h:27-82` | `AdditionalOutput` allocations | Retain the existing Rust container after reset-semantics correction. |
| `Code/GraphMol/Fingerprints/FingerprintGenerator.h:90-123` | common arguments and lifetime contract | Retain typed Rust arguments. |
| `Code/GraphMol/Fingerprints/FingerprintGenerator.h:292-425` | call arguments, generator ownership, options, scalar overloads, JSON declarations | Extract one family-neutral Rust core; Rust ownership replaces raw pointers. |
| `Code/GraphMol/Fingerprints/FingerprintGenerator.cpp:36-87` | common arguments, information string, JSON | Existing implementation is the mainline; regression-check and correct only proven gaps. |
| `Code/GraphMol/Fingerprints/FingerprintGenerator.cpp:90-173` | constructor/destructor/info and `reinitAdditionalOutput()` | Consolidate ownership and correct exact reset behavior. |
| `Code/GraphMol/Fingerprints/FingerprintGenerator.cpp:182-321` | generator JSON serialization/deserialization | Extend shared typed dispatch with torsion types. |
| `Code/GraphMol/Fingerprints/FingerprintGenerator.cpp:323-435` | environment-to-count-vector driver | Extract from Morgan as the sole shared driver. |
| `Code/GraphMol/Fingerprints/FingerprintGenerator.cpp:438-500` | provenance duplication and temporary output | Retain one shared implementation. |
| `Code/GraphMol/Fingerprints/FingerprintGenerator.cpp:503-652` | four scalar vector conversions | Retain one shared implementation. |
| `Code/GraphMol/Fingerprints/FingerprintGenerator.cpp:654-756` | ordered multi-molecule generation | Add one generic batch path over scalar methods. |
| `Code/GraphMol/Fingerprints/FingerprintGenerator.cpp:826-989` | default scalar/bulk `FPType` dispatch | Extend shared dispatch with torsions. |
| `Code/GraphMol/Fingerprints/TopologicalTorsionGenerator.h:21-136` | arguments, environment storage, contracts, factory defaults | Add distinct typed Rust torsion surface. |
| `Code/GraphMol/Fingerprints/TopologicalTorsionGenerator.cpp:25-65` | argument constructor/info/JSON | Port complete argument source unit. |
| `Code/GraphMol/Fingerprints/TopologicalTorsionGenerator.cpp:33-45` | result size | Port exact valid shifts; reject source-undefined widths structurally. |
| `Code/GraphMol/Fingerprints/TopologicalTorsionGenerator.cpp:67-99` | atom environment bit id and provenance | Port all four observed output branches. |
| `Code/GraphMol/Fingerprints/TopologicalTorsionGenerator.cpp:101-194` | environment generation | Port exact selection, repeat, correction, code/hash, and order behavior. |
| `Code/GraphMol/Fingerprints/TopologicalTorsionGenerator.cpp:196-210` | environment info/JSON | Port complete metadata unit. |
| `Code/GraphMol/Fingerprints/TopologicalTorsionGenerator.cpp:212-252` | two factories and 64-bit instantiation | One safe Rust factory; reject unsupported output widths. |
| `Code/GraphMol/Fingerprints/AtomPairs.h:197-271` | legacy declarations/defaults | Add distinct compatibility surface. |
| `Code/GraphMol/Fingerprints/AtomPairs.cpp:159-265` | legacy un-hashed fingerprint | Direct adapter over shared chemistry kernels with legacy size. |
| `Code/GraphMol/Fingerprints/AtomPairs.cpp:267-297` | `TorsionFpCalc` | Sole legacy-to-modern hashed-count adapter. |
| `Code/GraphMol/Fingerprints/AtomPairs.cpp:298-310` | legacy hashed count | Thin result adapter only. |
| `Code/GraphMol/Fingerprints/AtomPairs.cpp:312-347` | legacy hashed bit vector | Preserve block division and threshold branches. |
| `Code/GraphMol/Fingerprints/Wrap/TopologicalTorsionWrapper.cpp:21-92` | generator construction and typed options | Add Python construction/options over the Rust generator. |
| `Code/GraphMol/Fingerprints/Wrap/FingerprintGeneratorWrapper.cpp:43-320,327-772` | scalar/bulk methods, options, output, JSON | Reuse one Python wrapper pattern. |
| `Code/GraphMol/Descriptors/Wrap/rdMolDescriptors.cpp:250-309,1035-1070` | legacy torsion wrappers | Add native compatibility functions with source defaults. |
| `Code/GraphMol/Descriptors/Wrap/rdMolDescriptors.cpp:964-985` | AtomPairs constants and atom-code wrapper | Bind the shared Rust constants/kernel. |
| `rdkit/Chem/AtomPairs/Utils.py:18-107` | `ExplainAtomCode` and `GetAtomCode` alias | Port pure helper behavior without a second decoder table. |
| `rdkit/Chem/AtomPairs/Torsions.py:31-162` | aliases, score/explain/id helpers | Port pure helpers over shared code/vector functions. |

## Reference-Test Inventory

The source tests are a minimum branch seed, not a substitute for the planned
focused and 5,000-row oracle profiles.

| Source range | Covered upstream behavior |
|---|---|
| `Code/GraphMol/Fingerprints/test1.cpp:1510-1598` | legacy codes, counts, target sizes, hashed determinism |
| `test1.cpp:1729-1821` | rooted and ignored atoms |
| `test1.cpp:2128-2210` | custom invariants for all legacy forms |
| `test1.cpp:2537-2631` | chirality across un-hashed and hashed results |
| `test1.cpp:3051-3095` | explicit-H regression |
| `testFingerprintGenerators.cpp:1574-2218` | modern/legacy relationship, four vector forms, counts, roots, ignores, chirality |
| `testFingerprintGenerators.cpp:2257-2366` | explicit-H, long-path, and cycle regressions |
| `testFingerprintGenerators.cpp:2380-2420` | generic `FPType::TopologicalTorsionFP` bulk dispatch |
| `catch_tests.cpp:404-447` | `AdditionalOutput` atom counts/to-bits/bit paths |
| `catch_tests.cpp:755-772` | `onlyShortestPaths` behavior |
| `catch_tests.cpp:850-873` | generator JSON roundtrip |
| `fpgen_catch_tests.cpp:49-85` | mutable common generator options |
| `fpgen_catch_tests.cpp:523-539` | `atomsPerBit` output |
| `Wrap/testGenerators.py:107-114,142-203` | Python scalar and generic bulk dispatch |
| `Wrap/testGenerators.py:252-257` | custom count bounds |
| `Wrap/testGenerators.py:333-343` | mutable `onlyShortestPaths` |
| `Wrap/testGenerators.py:394-401` | Python `atomsPerBit` |

## Public Option Mapping

### Modern construction and mutable options

| RDKit option | Pinned default | COSMolKit planned field | Output effect |
|---|---:|---|---|
| `includeChirality` | `false` | `include_chirality` | Adds two chiral bits to atom-code width and reads R/S CIP labels. |
| `torsionAtomCount` | `4` | `torsion_atom_count` | Atom count in each enumerated path and packed result width. |
| `atomInvariantsGenerator` | null | typed optional generator | Null creates `AtomPairAtomInvGenerator(includeChirality, true)`. |
| `countSimulation` | `true` | common `count_simulation` | Expands counts into bound-index bits for sparse/explicit bit forms. |
| `countBounds` | `[1,2,4,8]` | common `count_bounds` | Count-simulation thresholds; comparison is `count >= bound`. |
| `fpSize` | `2048` | common `fp_size` | Folded count/explicit-bit size; sparse-count ignores it. |
| `ownsAtomInvGen` | `false` | Rust ownership enum/value | Lifetime-only difference; safe ownership must preserve clone semantics. |
| `onlyShortestPaths` | `false` | torsion option | Prunes a path when endpoint shortest distance is shorter than its length. |
| `numBitsPerFeature` | `1` | common mutable option | Adds deterministic RNG-derived bits per environment. |

### Per-call arguments

| RDKit argument | Exact torsion interpretation |
|---|---|
| `fromAtoms=nullptr` | No root restriction. |
| `fromAtoms=[]` | Non-null empty set; no path has an allowed endpoint, so no environments. |
| non-empty `fromAtoms` | Keep a path when its first or last atom is selected; internal selection alone is insufficient. |
| `ignoreAtoms=nullptr` | No exclusion. |
| `ignoreAtoms=[]` | Non-null empty set with no excluded atom. |
| non-empty `ignoreAtoms` | Reject a path when any atom in it is selected. |
| `confId` | Accepted by common API and ignored by torsion environment generation. |
| `customAtomInvariants` | Copied verbatim, then each path element is normalized by the torsion environment logic. |
| `customBondInvariants` | Common API accepts it; torsion generation does not consume bond invariants. |
| `AdditionalOutput` | Supports `atomToBits`, `atomCounts`, `bitPaths`, and observed `atomsPerBit`; not `bitInfoMap`. |

### Four scalar outputs

| Output | Shared source behavior |
|---|---|
| sparse count | Full direction-canonical packed torsion ids, no `fpSize` hashing, integer counts. |
| sparse bit | Result size capped at `u32::MAX`; uses folded environment ids and optional count simulation. |
| hashed count | Hashes environments into `fpSize`; stores exact collision-merged counts. |
| explicit bit | Uses `fpSize`, optional count simulation, and deterministic extra bits. |

All four bulk forms must preserve input order and null positions and produce the
same value as calling the corresponding scalar function on each valid record.

### Legacy defaults

| Entry point | Defaults and compatibility rules |
|---|---|
| un-hashed | `targetSize=4`, null selections/invariants, `includeChirality=false`; sparse length is `2^width - 1`. |
| hashed count | `nBits=2048`, `targetSize=4`, null selections/invariants, `includeChirality=false`. |
| hashed bit | same plus `nBitsPerEntry=4`; block length is integer `nBits/nBitsPerEntry`. |

For legacy hashed bits, `nBitsPerEntry==4` uses thresholds `[1,2,4,8]`.
Otherwise bit `i` is set when `count > i`, which is intentionally different
from the modern `count >= bound` simulation rule.

## COSMolKit Integration and Same-Source Decisions

| Current COSMolKit artifact | Audit status | Required mainline decision |
|---|---|---|
| `AdditionalOutput` at `fingerprint.rs:134` | reusable with one proven reset mismatch | Keep one type; correct shared reset semantics. |
| `FingerprintArguments` at `fingerprint.rs:199` | reusable common arguments | Keep one type and extend typed dispatch only. |
| `FingerprintFuncArguments` at `fingerprint.rs:417` | reusable per-call value holder | Keep one type; validate selection/invariant lengths before indexed access. |
| `SparseCountFingerprint` at `fingerprint.rs:427` | reusable `u64` key/count vector | Keep one implementation. |
| `SparseBitFingerprint` at `fingerprint.rs:466` | reusable sparse bit vector | Keep one implementation. |
| `Fingerprint` at `fingerprint.rs:2541` | reusable explicit bit vector | Keep one implementation. |
| `MorganFingerprintGenerator` methods at `fingerprint.rs:1914-2410` | generic source driver is family-bound | Extract the generic driver; Morgan delegates unchanged. |
| `duplicate_additional_output_bit` and `setup_temp_additional_output` | active shared-source logic but stored in Morgan-oriented module flow | Retain once behind family-neutral generator core. |
| `RdkitFingerprintMtRng` at `fingerprint.rs:495` | reusable architecture; exact typedef parameter mapping unresolved | Keep one RNG only after pinned sequence tests close the discrepancy. |
| `hash_combine` at `fingerprint.rs:3309` | already shared and used by other fingerprint/hash code | Call directly from torsion hash; never duplicate formula. |
| `build_fingerprint`/`fold_invariant`/`atom_is_excluded` at `fingerprint.rs:3320-3679` | unreachable historical Morgan branch; only internal calls are within that branch | Delete before torsion integration. |
| `compute_topological_distances` at `distgeom.rs:6170` | source-backed but private to distance geometry | Move into one private graph helper and share with both callers. |
| `enumerate_rdkit_fp_paths` at `fingerprint.rs:4092` | ports RDKFingerprint-specific path behavior | Keep separate; it is not a torsion path-enumeration shortcut. |
| `TopologicalFingerprintParams` and `topological_fingerprint` | existing `RDKFingerprintMol` family | Preserve API and semantics; Topological Torsion gets distinct names. |
| `FINGERPRINT_FEATURE` | describes Morgan; not an adequate torsion status boundary | Add a distinct Topological Torsion `FeatureSpec` only after parity closes. |
| Python Morgan/RDK fingerprint bindings | useful wrapper patterns, not a torsion chemistry implementation | Reuse conversion/provenance patterns; call the Rust torsion core. |
| Morgan batch API | reusable order/error design | Add torsion batch adapters over scalar torsion generation. |
| existing fingerprint oracle tools | reusable deterministic schema style | Add a new pinned torsion generator/profile, not fixture lookup in production. |

### Proven shared-core discrepancies

1. Pinned `reinitAdditionalOutput()` clears `atomCounts`, `atomToBits`,
   `bitInfoMap`, and `bitPaths`, but does not clear `atomsPerBit`.
   `AdditionalOutput::reset_for_atom_count()` currently clears
   `atoms_per_bit`. Repeated-use parity must decide and reproduce the pinned
   behavior in the single shared reset implementation.
2. The pinned generator typedef at `FingerprintGenerator.cpp:378-381` ends in
   parameter `3346425566U`; the current Rust RNG declares `F=1812433253` and
   comments that the Boost adapter forwards that value. Existing RDKFingerprint
   parity suggests adapter details may matter. A pinned output-sequence test,
   not inspection alone, must determine the exact shared engine behavior before
   torsion `numBitsPerFeature>1` is claimed.
3. The active Morgan-bound `getFingerprintHelper()` returns
   `UnsupportedOption` when chirality is requested without `_StereochemDone`,
   while pinned RDKit uses a private molecule copy. The shared refactor must
   close this branch once, preserve source immutability, and reproduce the
   original `mol` versus `lmol` call flow.

## Prerequisite Closure

| Prerequisite | Current state | Closure required before dependent function |
|---|---|---|
| generic generator core | embedded in `MorganFingerprintGenerator` | Extract and regression-test before torsion generator factory. |
| private-copy stereo preparation | explicit unsupported branch | Close in shared generator core before chirality parity. |
| Atom `numPiElectrons` | no source-backed helper found by symbol audit | Port before `getAtomCode`. |
| AtomPairs constants and atom code | absent | Port once before invariant generator and Python explanation helper. |
| ordered torsion path enumeration | absent | Port exact Subgraphs functions; do not reuse RDK path enumeration. |
| shortest-path matrix | exact implementation exists privately | Move before `onlyShortestPaths` path port. |
| torsion code/hash | absent | Port after atom constants; hash calls existing `hash_combine`. |
| generator JSON type dispatch | Morgan/RDK families exist | Extend after torsion arguments/env/invariant types exist. |
| Python common vectors/provenance | exists for current families | Reuse after Rust public torsion types exist. |
| deterministic reference generator | no torsion-specific artifact | Add before focused/corpus parity tests. |
| distinct feature status | absent | Add only after complete parity evidence. |

The broader Atom Pair fingerprint family is not a prerequisite. Only its
source constants, atom code, and invariant generator reached by Topological
Torsion are included.

## Error Ledger

| Input or source failure | Pinned source behavior | Rust policy mapping |
|---|---|---|
| `countSimulation=true`, empty bounds at construction | precondition/value error | `FingerprintError::InvalidArguments`. |
| explicit bit with empty bounds after option mutation | `ValueErrorException` | structured invalid-arguments error. |
| count-bound count `>= fpSize` | `ValueErrorException` | structured invalid-arguments error. |
| `numBitsPerFeature==0` | constructor precondition | structured invalid-arguments error. |
| `fpSize==0` for folded/explicit output | source precondition or impossible vector | structured invalid-arguments error; sparse-count remains valid. |
| custom invariant length shorter than atom count | unchecked source indexing | reject with expected/actual length context; never panic. |
| selection index out of range | source bitset range/precondition failure | reject with index and atom-count context. |
| torsion width shift at or beyond output width | C++ undefined behavior | structured invalid-width error; no fabricated output. |
| `torsionAtomCount==0` or other path-length boundary | depends on Subgraphs source branch | oracle-classify before support; structured error only for undefined branch. |
| legacy `nBitsPerEntry==0` | divide by zero in source adapter | structured invalid-arguments error. |
| legacy `nBits<nBitsPerEntry` | zero block size reaches downstream code | oracle/source-classify and reject any undefined branch. |
| non-divisible legacy `nBits` | integer division leaves a partial tail | reproduce source-defined valid behavior exactly. |
| malformed/unknown generator JSON | property-tree/type dispatch failure | `InvalidArgumentsJson` or a new typed generator-kind error retaining cause. |
| unsupported 32-bit torsion output generator | no upstream instantiation | structured unsupported-output-width error. |
| query or unmodeled atom state needed by atom code | source branch must be inspected | exact modeled branch or structured unsupported state; no heuristic bucket. |
| stereo/CIP preparation failure | source operation failure/precondition | retain structured stereo/CIP cause and source immutability. |
| null molecule in bulk input | upstream preserves a null output position | Rust batch result preserves input index with typed null/invalid record. |
| worker failure | upstream thread path preserves ordering | indexed structured error; no silent dropping. |
| allocation size cannot fit Rust address space | possible source allocation failure | structured size/capacity error, never wrap or panic. |

`FingerprintError` currently has broad `InvalidArguments`, JSON, CIP,
unsupported-option, and bit-length variants. The port should add narrow typed
variants where an index, expected length, width, batch position, or capacity is
known instead of erasing that information into static strings.

## Line-Coverage Ledger

Status here records the pre-port state observed during this audit.

| Source unit | Current status | Planned closure step(s) |
|---|---|---|
| common vector and argument types | partial/existing | 11, 14 |
| generic generator helper and scalar outputs | existing but Morgan-bound with open branches | 11, 14 |
| graph distance matrix | implemented but incorrectly scoped for reuse | 19, 22 |
| `numPiElectrons` | absent | 27, 36 |
| `getAtomCode` | absent | 30, 36 |
| `AtomPairAtomInvGenerator` | absent | 33, 36 |
| torsion code/hash | absent | 41, 44, 47 |
| `extendPaths` | absent | 52, 64 |
| `pathFinderHelper` | absent | 55, 64 |
| `findAllPathsOfLengthsMtoN` | no torsion-compatible implementation | 58, 64 |
| `findAllPathsOfLengthN` | absent | 61, 64 |
| torsion arguments/result size | absent | 69, 72, 75 |
| torsion atom environment | absent | 80, 89 |
| torsion environment generator | absent | 83, 86, 89 |
| modern generator factories | absent | 94, 103 |
| shared JSON/FPType dispatch | no torsion variant | 97, 100, 103 |
| legacy un-hashed | absent | 108, 120 |
| legacy hashed adapter/count/bit | absent | 111, 114, 117, 120 |
| public Rust surface | absent | 125, 128 |
| Python generator/wrappers/helpers | absent | 132, 134, 136 |
| focused deterministic oracle/parity | absent | 140-147 |
| 5,000-row parity | absent | 149-156 |
| ordered batch API | absent | 159, 162 |
| docs/stubs/feature status | absent | 166-171 |
| final source/performance closure | absent | 173-191 |

No Topological Torsion source block may be marked behaviorally complete merely
because an existing RDKFingerprint or Morgan helper looks similar. Every row
above closes only when its own source lines and edge tests are present.

## Mainline Decision

The executable mainline is:

1. remove the unreachable historical Morgan constructor branch;
2. extract and correct the one shared generator output core;
3. share the existing exact topological-distance kernel;
4. port the shared atom invariant and torsion code/hash prerequisites;
5. port the exact ordered Subgraphs path functions;
6. port torsion arguments, environments, factories, and shared dispatch;
7. layer legacy, Rust, Python, and batch adapters over those kernels;
8. close focused and maintained-corpus parity before publishing a distinct
   supported feature status.

This ordering prevents the historical Morgan branch, existing RDKFingerprint
path implementation, or wrapper convenience code from becoming a competing
chemistry implementation. The completed tree must contain exactly one active
implementation of common vector storage, generator assembly, count simulation,
provenance duplication, fingerprint RNG, Boost-compatible hashing, and
topological distance matrices.
