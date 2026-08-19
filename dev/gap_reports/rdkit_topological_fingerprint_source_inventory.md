# RDKit Topological Fingerprint Source Inventory

Status: source inventory retained as provenance for the completed source-backed
port. The selected `RDKFingerprintMol` call graph is implemented and validated;
the public support declaration is upgraded only by the release/support-matrix
step after this report.

## Pinned source

- RDKit source tree: `third_party/rdkit`.
- Pinned release/revision: RDKit `2026.03.1`,
  `351f8f378f8ad6bbd517980c38896e66bf907af`.
- Public boundary: `RDKFingerprintMol()` in
  `Code/GraphMol/Fingerprints/Fingerprints.cpp:178-248` and the Python wrapper
  in `Code/GraphMol/Wrap/MolOps.cpp:2513-2527`.

## Active call graph

1. `RDKFingerprintMol()` validates `minPath`, `maxPath`, `fpSize`,
   `nBitsPerHash`, `atomInvariants`, and `atomBits`, constructs an
   `RDKitFPEnvGenerator`, sets `d_fpSize` and `d_numBitsPerFeature`, maps
   `atomInvariants`/`fromAtoms` into `FingerprintFuncArguments`, allocates
   requested `AdditionalOutput` holders, invokes `getFingerprint`, copies
   `atomToBits` and `bitPaths` to the legacy output containers, then repeatedly
   folds by two while the requested density is not reached
   (`Fingerprints.cpp:178-248`).
2. `getRDKitFPGenerator()` constructs `RDKitFPArguments` and an
   `RDKitFPEnvGenerator`, installs the default `RDKitFPAtomInvGenerator`, and
   creates the generic `FingerprintGenerator`
   (`RDKitFPGenerator.cpp:244-295`).
3. `FingerprintGenerator::getFingerprintHelper()` optionally clones and
   stereochemically labels the molecule, reinitializes additional output,
   obtains atom/bond invariants, asks the environment generator for all
   environments, allocates the sparse result, and emits one base bit plus
   `numBitsPerFeature - 1` deterministic random bits per environment
   (`FingerprintGenerator.cpp:325-435`).
4. `RDKitFPEnvGenerator::getEnvironments()` calls
   `RDKitFPUtils::enumerateAllPaths()`, identifies query bonds, generates bond
   hashes for every ordered path, sorts multi-bond hashes, appends the number
   of distinct atoms, hashes the range, and creates `RDKitFPAtomEnv` objects
   (`RDKitFPGenerator.cpp:180-242`).
5. `enumerateAllPaths()` dispatches to
   `findAllSubgraphsOfLengthsMtoN()` for branched paths or
   `findAllPathsOfLengthsMtoN()` for linear paths, preserving the map and
   insertion ordering, including repeated `fromAtoms`
   (`FingerprintUtil.cpp:283-321`). The implementation lives in
   `Code/GraphMol/Subgraphs/Subgraphs.cpp:347-552`; those helper bodies are in
   scope under the source-reproduction protocol.
6. `generateBondHashes()` resets the atom membership bitset, counts path-local
   atom degrees, rejects a path containing a query bond, canonicalizes endpoint
   invariant/degree ordering, selects aromatic or ordinary bond order, and
   applies `gboost::hash_combine` in the exact field order
   (`FingerprintUtil.cpp:345-434`).
7. `RDKitFPAtomEnv::updateAdditionalOutput()` records bit paths and updates
   `atomsPerBit`, deduplicated `atomToBits`, and atom counts in source order
   (`RDKitFPGenerator.cpp:120-149`).
8. Generic output uses the customized Boost Mersenne Twister and
   `boost::uniform_int<>(0, INT_MAX)`, seeded with `42` and reseeded with the
   environment seed (`FingerprintGenerator.cpp:373-430`).

## Public option mapping

| Python/RDKit option | Rust field required | Default | Current state |
| --- | --- | ---: | --- |
| `minPath` | `min_path: u32` | 1 | present, not executed |
| `maxPath` | `max_path: u32` | 7 | present, not executed |
| `fpSize` | `fp_size: u32` (or checked `usize` adapter) | 2048 | current `n_bits: usize`; type correction required |
| `nBitsPerHash` | `num_bits_per_feature: u32` | 2 | current name `n_bits_per_hash`; rename/mapping required |
| `useHs` | `use_hs: bool` | true | missing from current placeholder |
| `tgtDensity` | `target_density: f64` | 0.0 | missing |
| `minSize` | `min_size: u32` | 128 | missing |
| `branchedPaths` | `branched_paths: bool` | true | missing |
| `useBondOrder` | `use_bond_order: bool` | true | current `use_bond_types`; source naming/semantics correction required |
| `atomInvariants` | checked `Vec<u32>` | `None` | missing typed field |
| `fromAtoms` | checked `Vec<u32>` | `None` | current `Option<Vec<usize>>`; conversion and source bounds required |
| `atomBits` | typed `Vec<Vec<u32>>` provenance output | omitted | result type required |
| `bitInfo` | typed ordered map `u32 -> Vec<Vec<i32>>` | omitted | result type required |

`ignore_atoms` is not a parameter of `RDKFingerprintMol`; it must not remain
in the public topological fingerprint API. The generic generator has an
internal `ignoreAtoms` argument, but the selected legacy boundary never passes
one and the port must not expose it as a compatibility invention.

The Python wrapper defaults and keyword names are pinned by
`MolOps.cpp:2513-2523`. C++ overload defaults are the same for the selected
boundary except that the returned legacy output containers are allocated only
when the caller supplies `atomBits` or `bitInfo`.

## Output and error branches

### Preconditions

The source rejects `minPath == 0`, `maxPath < minPath`, `fpSize == 0`, and
`nBitsPerHash == 0` with `PRECONDITION` failures
(`Fingerprints.cpp:186-193`; constructor checks also occur at
`RDKitFPGenerator.cpp:103-118`). A Rust public wrapper must convert these to
the existing structured invalid-argument error, without panic or fallback.
`atomInvariants` and `atomBits` must contain at least `mol.getNumAtoms()`
entries. `fromAtoms` values are consumed by the source subgraph helpers and
therefore need checked atom-index behavior matching those helpers.

### Query paths

`generateBondHashes()` returns no hashes when any bond in a candidate path, or
either endpoint atom, is a complex query (`identifyQueryBonds`,
`FingerprintUtil.cpp:323-343`; hash rejection at `:355-370`). This is not a
generic "skip unsupported molecule" heuristic; it is an exact source branch
and must be retained with structured behavior at the public boundary.

### Density folding

The result is folded by exactly a factor of two while
`on_bits / num_bits < tgtDensity` and `num_bits >= 2 * minSize`
(`Fingerprints.cpp:236-247`). Folding happens after provenance is copied and
therefore provenance bit identifiers remain those emitted before folding; this
ordering must be tested explicitly.

### Provenance

`atomBits` is cleared and rebuilt from `AdditionalOutput::atomToBits`, while
`bitInfo` is cleared and populated from `AdditionalOutput::bitPaths`
(`Fingerprints.cpp:217-234`). The generic generator invokes provenance updates
for the base bit and every random expansion bit (`FingerprintGenerator.cpp:400-430`).
`RDKitFPAtomEnv::updateAdditionalOutput()` deduplicates an atom's bit list but
does not deduplicate `bitInfo` paths (`RDKitFPGenerator.cpp:120-149`).

## Source-line coverage ledger

| Source unit | Pinned ranges | Required Rust unit | Status |
| --- | --- | --- | --- |
| Legacy boundary and preconditions | `Fingerprints.cpp:178-248` | public `topological_fingerprint` adapter | ✔️✔️ |
| RDKitFP argument defaults/JSON | `RDKitFPGenerator.h:20-58`, `RDKitFPGenerator.cpp:77-118` | typed params and validation | ✔️✔️ |
| Default atom invariants | `FingerprintUtil.cpp:271-280`, `RDKitFPGenerator.h:60-71` | atom invariant generator | ✔️✔️ |
| Path/subgraph dispatch | `FingerprintUtil.cpp:283-321`; `Subgraphs.cpp:347-552` | ordered path enumerator | ✔️✔️ |
| Query-bond cache | `FingerprintUtil.cpp:323-343` | path-local query rejection | ✔️✔️ |
| Bond hash preparation | `FingerprintUtil.cpp:345-434` | canonical bond hash inputs | ✔️✔️ |
| Environment seed and provenance | `RDKitFPGenerator.cpp:120-242` | environment collection/output | ✔️✔️ |
| Generic result/RNG | `FingerprintGenerator.cpp:325-435` | bit expansion and counts | ✔️✔️ |
| Density fold | `Fingerprints.cpp:236-247` | explicit-bit fold loop | ✔️✔️ |
| Python surface | `MolOps.cpp:2513-2527` | Python wrapper and typed outputs | ✔️✔️ |

## Reusable COSMolKit machinery and gaps

`fingerprint.rs` already contains a source-backed generic Morgan
`FingerprintGenerator`, `AdditionalOutput`, checked JSON parsing, sparse and
explicit projections, and a deterministic implementation of the customized
RDKit RNG used by the Morgan path. These components are reusable only where
their source contracts match. In particular:

- the RDK environment generator needs a separate atom/path environment type;
- the RDK atom invariant formula is not Morgan's connectivity invariant;
- RDK uses Boost `hash_combine` over sorted per-bond fields and a distinct
  path enumerator;
- RDK's `atomToBits`/`bitPaths` shapes and legacy copying order differ from
  Morgan's `bitInfoMap` output;
- density folding is a legacy RDK-only post-processing branch.

No local path hash, fallback fingerprint family, approximation, filtering, or
corpus-specific special case is source-equivalent and may be used to close the
gap.

## Required focused evidence

Before enabling support, the plan must add exact tests for invalid ranges,
default and every public option, linear/branched/cyclic/fused/disconnected
paths, explicit hydrogens, repeated and restricted `fromAtoms`, custom
invariants, query rejection, random-bit collisions, density folding and
minimum-size boundaries, typed `atomBits`/`bitInfo` provenance, deterministic
repeated calls, and all active 5,000-row profiles. The validation report must
record source revision, corpus checksum, commands, timing, and machine
artifacts under ignored `tmp/parity-audit/`.
