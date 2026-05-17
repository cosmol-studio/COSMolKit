# SMILES Writer Port-Rigor Repair

Date: 2026-05-17

## Scope

This note records the remaining non-rigorous RDKit writer/parser alignment gaps
found during the SMILES writer parity repair pass and the source-backed repair
plan for each one.

Relevant Rust files:

- `crates/cosmolkit-core/src/notation/smiles.rs`
- `crates/cosmolkit-core/src/notation/smiles_write.rs`

Relevant RDKit sources:

- `third_party/rdkit/Code/GraphMol/Chirality.cpp`
- `third_party/rdkit/Code/GraphMol/Canon.cpp`

## Findings Before Repair

### 1. Parser did not persist ring-stereo computed props

Rust location:

- `smiles.rs:3870-3904`

RDKit source basis:

- `Chirality.cpp:1521-1641` `findChiralAtomSpecialCases()`

Problem:

- COSMolKit used `find_chiral_atom_special_cases()` only as an internal cleanup
  decision input.
- RDKit also persists computed `_ringStereochemCand` and `_ringStereoAtoms`
  props on the molecule.
- Because parser output did not carry those props, the writer later lacked the
  same state RDKit canonicalization expects to reuse.

Required repair:

- Persist `_ringStereochemCand` and `_ringStereoAtoms` during parser stereo
  finalization on the same source-backed special-case path.

### 2. Writer used local reconstruction when `_ringStereoAtoms` was absent

Rust location:

- `smiles_write.rs:4744-4894`

RDKit source basis:

- `Canon.cpp:1544-1564`

Problem:

- RDKit canonical fragment writing only reuses `_ringStereoAtoms` if the prop
  is already present.
- COSMolKit added a single-fragment fallback that recomputed ring stereo
  neighbors locally when the prop was absent.
- That fallback was not source-native behavior; it existed only because parser
  output had failed to persist the computed prop.

Required repair:

- Remove the local reconstruction fallback.
- Require the prop on supported writer paths.
- Treat corrupted ring-stereo state as an explicit error instead of silently
  synthesizing replacement data.

### 3. Double-bond canonicalization markers lagged the actual implementation

Rust location:

- `smiles_write.rs:5289-5658`

RDKit source basis:

- `Canon.cpp:1187-1380`

Problem:

- The code path already closely matched RDKit structure and ordering, but the
  copied-source markers still advertised it as `RDKit❗❌`.
- That made the implementation look heuristic even where the Rust body now
  follows the RDKit control flow directly.

Required repair:

- Re-review against RDKit source and upgrade markers only where the local Rust
  implementation is actually source-aligned.

### 4. Parser stereochemistry finalization still had two non-source-native CIP paths

Rust locations:

- `smiles.rs:3710-3908`
- `stereo.rs:859-931`
- `stereo.rs:3469-3567`

RDKit source basis:

- `Chirality.cpp:1651-1821` `isAtomPotentialChiralCenter()` /
  `assignAtomChiralCodes()`
- `Chirality.cpp:2067-2117` `rerankAtoms()`
- `Chirality.cpp:2267-2445` `legacyStereoPerception()`

Problem:

- COSMolKit's parser finalization had started to model RDKit's fixed-point
  legacy stereo loop, but two details were still off the source path:
  1. `setBondStereoFromDirections()` had been injected into the coordinate-free
     SMILES parser finalization even though RDKit's
     `legacyStereoPerception()` does not pre-assign double-bond stereo there.
  2. `assign_atom_chiral_codes()` still used a simplified local route:
     neighbor records carried atom indices instead of bond indices, and
     `nSwaps` came from `chiral_permutation.unwrap_or(0)` instead of
     RDKit's `Atom::getPerturbationOrder()` over the current rank-ordered bond
     list.
  3. `rerank_atoms()` only recognized `BondStereo::E/Z`, while the parser path
     records the same stereo information as `BondStereo::Trans/Cis`.

Required repair:

- Keep coordinate-free parser finalization on the
  `legacyStereoPerception()` path and remove the premature
  `setBondStereoFromDirections()` call.
- Make `assign_atom_chiral_codes()` use bond-index neighbor records plus an
  explicit `getPerturbationOrder()`-equivalent swap count.
- Make `rerank_atoms()` consume the parser's `Trans/Cis` bond stereo state as
  the same RDKit E/Z ordering signal.

## Repair Rules For This Pass

- No heuristic output-path widening.
- No parity-preserving hacks without direct RDKit source support.
- Any counterexample created only by corrupted internal ring-stereo props must
  fail explicitly instead of producing guessed SMILES.
- Any counterexample created only by corrupted internal perturbation-order
  inputs must fail explicitly instead of silently falling back to a guessed
  swap parity.
