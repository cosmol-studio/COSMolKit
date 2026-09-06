# cosmolkit-core source layout

This crate is organized by ownership boundary, not by porting order.

- `model/`: value-style molecule storage, row ids, builders, read-only views,
  source invariants, and core error/state types.
- `operations/`: registered mutation/finalization machinery. Existing
  molecules must be changed here, not by exposing mutable model storage.
- `chemistry/`: chemistry algorithms over the model: valence, aromaticity,
  rings, hydrogens, stereo, coordinates, transforms, and distance geometry.
- `notation/`: textual notations and sequence/fragment construction: SMILES,
  canonical ranking/writing, sequence, and fragment helpers.
- `search/`: query trees, SMARTS parsing, and substructure matching.
- `properties/`: fingerprints, drawing, hashing, pickling, and batch helpers.
- `bio/`: `BioStructure`, biomolecular invariants, and BioStructure operations.
- `io/`: file-format readers/writers. Structural PDB/mmCIF goes through
  `io::bio`; molecule molblock/SDF code stays in `io::molblock` and `io::sdf`.
- Historical `unported/` quarantine modules have been deleted. Do not
  reintroduce parallel inactive stubs or use deleted quarantine paths as
  implementation shortcuts.

`lib.rs` keeps root-level compatibility re-exports such as `crate::atom`,
`crate::ops`, and `crate::smiles_write`. New code should prefer the domain
directory when navigating, but it should not bypass operation or porting policy
just because a root alias exists.
