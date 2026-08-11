# Feature Parity Scope

This page summarizes how rigorously the main COSMolKit features are checked
against RDKit. It is organized by user-visible capability: what the feature
does, how strict the comparison is, and where the current boundary is.

RDKit `2026.03.1` is the current reference for chemistry parity. A feature is
described as parity-covered only at the boundary listed here; open branches stay
explicitly unfinished rather than being counted as compatible behavior.

## How To Read This Table

- **Strict parity** means 100% agreement on every known covered case for that
  feature boundary.
- We do not report pass percentages for strict parity surfaces. Any mismatch is
  a failure to investigate, not an acceptable error rate.
- **RNG-sensitive parity** means the test fixes seed and execution conditions and
  compares deterministic coordinates or optimizer state.
- **Final-output parity** means the serialized user-visible artifact is compared
  directly.
- **Unfinished** means the API may exist, but it must not be presented as
  RDKit-compatible until the listed exact parity boundary passes.
- Default parity tests use the small `smiles_small` corpus for daily feedback.
  Most stricter profiles use the same assertions with separate expected data
  under `testdata/<domain>/expected/rdkit/smiles_5000/`. Force fields use an
  explicit two-tier strict boundary: all 5,000 rows cover parameter
  availability, explicit-H behavior, MMFF atom types, and MMFF charges, while
  all 150 curated rows cover initial energy, full gradient, and single- and
  multi-conformer optimization. The tiers are separate tests and neither is
  reported as coverage of the other's behavioral surface.
- Golden baselines are condition-bound artifacts. A corpus, RDKit version,
  branch-matrix, normalization, RNG/coordinate setting, or compared-field change
  requires a new profile directory or file name. Reusing a name for incompatible
  conditions is not allowed.

| Feature | Current rigor | What is compared | Important boundary |
|---|---:|---|---|
| Molecular graph state | Strict parity | Atom identity, charge, degree, hydrogens, radicals, hybridization, aromaticity, ring membership, bond type, bond stereo, conjugation, CIP rank/code, direct molecules, and explicit-H molecules. | This is the base state other chemistry features depend on. |
| SMILES writing | Strict branch-matrix parity | RDKit `MolToSmiles()` output across the full writer parameter matrix: isomeric on/off, kekule on/off, canonical on/off, clean-stereo on/off, explicit-bond on/off, explicit-H on/off, dative-bond inclusion on/off, atom-map ignoring on/off, and rooted output for none/first/last atom. Canonicalized SMILES are compared directly as part of the canonical branches. | The exhaustive all-combination matrix is the parity target; everyday runs cover common branches and the full matrix is run explicitly when auditing writer parity. Unsupported writer branches must remain explicit. |
| Distance-geometry bounds | Strict numeric parity | Full bounds-matrix shape and every lower/upper-bound entry, including topology, ring, macrocycle, amide, VDW, and triangle-smoothing effects. | Entry tolerance is `1e-8`; this is matrix equality, not a smoke test. |
| 3D conformer generation | RNG-sensitive coordinate parity | DG/KDG/ETDG/ETKDG presets, fixed seeds, single-thread execution, random-coordinate embedding, coordMap, custom pairwise terms, bounds-matrix handling, timeout/max-iteration controls, multi-conformer generation, and sequential-seed policy. | Successful seeded embeddings must reproduce RDKit coordinates within `1e-6`. |
| UFF/MMFF force fields | Tiered state and optimizer parity | All 5,000 strict-corpus rows compare UFF/MMFF parameter availability, explicit-H behavior, MMFF atom types, and MMFF formal/partial charges. All 150 curated rows compare seeded initial coordinates, initial energy, every gradient component, single-conformer optimization, multi-conformer optimization, result codes, final energies, and final coordinates. | Energy, gradient, and coordinate tolerances remain `1e-6`. The 5,000-row coverage tier is not presented as optimizer parity; optimizer parity is exhaustive on the curated corpus and all locked regressions. |
| Molecular descriptors | Strict source-backed descriptor parity | Molecular weight, exact molecular weight, formula variants, hydrogen-bond donor/acceptor counts, fraction Csp3, Crippen LogP/MR, TPSA with and without S/P contributions, aromatic-ring count, rotatable-bond counts for default/non-strict/strict/strict-linkages modes, and QED. | Float descriptors are compared by exact RDKit bit pattern from the golden baseline; integer and string descriptors are exact equality. QED is pinned to the complete reference runtime `RDKit 2026.03.1 + CPython 3.13.12` because RDKit delegates its reductions to Python `sum()`, whose float algorithm differs before CPython 3.12. Unsupported or failed-closed descriptor calls are parity failures for the covered fields. |
| 2D depiction preparation | Strict prepared-state parity | RDKit-style preparation before drawing: kekulization for depiction, chiral-H insertion, wedge-bond assignment, prepared atom order, prepared 2D coordinates, bond order, aromatic flags, and bond directions. | This is the prepared drawing boundary, not a blanket claim that every standalone `Compute2DCoords` branch is complete. |
| SVG drawing | Final-output parity | Final SVG text after normalizing only tool namespace/prefix metadata. The comparison covers preparation, scaling, atom labels, bond geometry, text metrics, metadata, CSS/data attributes, and rendered paths present in the covered surface. | This is final-SVG parity for the covered depiction surface. Link-node extraction, StereoGroup masking, atomRegions, and other marker-open drawing branches remain unfinished. |
| PNG drawing | Rasterization consistency | PNG export is checked as local SVG rasterization through the Rust rendering stack. | PNG is not claimed to be RDKit Cairo/Qt bit parity. |
| MolBlock/SDF writing | Strict output parity for supported branches | V2000/V3000 output, 2D and 3D sources, stereochemistry, kekulization, SGroups, RGroups, aliases, value lines, and aromatic-bond bookkeeping in the covered writer branches. | Unsupported writer branches must fail closed or stay explicitly marked. |
| MolBlock/SDF reading | Strict parse/state parity for supported branches | Topology, atom fields, bond fields, 2D coordinates, 3D coordinates, marker-derived chirality, coordinate-derived chirality, delayed sanitize/remove-H paths, and SMILES round-trip checks. | Reader coverage is broad but not a blanket claim for every CTAB/SDF extension. |
| MOL2 reading | Source-backed parser parity | Parser parameters, topology, atom and bond fields, 3D coordinates, chirality, CORINA variant behavior, cleanup-substructure behavior, and SMILES output. | Source-ported for the exposed Tripos MOL2 profile; broader marker audit work remains tracked separately. |
| InChI scalar APIs | Strict source-defined field and branch parity | `Chem.MolToInchi`, `Chem.MolToInchiKey`, root `InchiToInchiKey`, and `Chem.MolFromInchi`; official return fields, diagnostics, graph state, atom/bond/stereo/isotope/H/charge/radical fields, options, cleanup, malformed input, source preservation, and concurrent calls. | Exact parity applies only where pinned official InChI v1.07.5 and RDKit 2026.03.1 define behavior. The official `NormalizeAndCompare` initial-buffer allocation-failure path is undefined C behavior; Rust returns a deterministic structured `allocation_failed` error instead. MolBlock, SDF/V3000, IXA, AuxInfo, INCHIGEN, version query, and extended-polymer APIs are frozen or unsupported. |
| XYZ reading | Strict simple-format parity | Atom identities, one 3D conformer, coordinates, and absence of inferred bonds. | XYZ is coordinate input, not a chemistry perception format. |
| Atom chiral tags from 3D structure | Strict full-state parity | Public `with_chiral_tags_from_structure` behavior over 77 fixed RDKit oracle records: operation status/errors, selected conformer, every atom chiral tag and permutation, `_NonExplicit3DChirality`, `_StereochemDone`, complete atom/bond state, conformer coordinates/dimensionality, replacement behavior, environment-controlled non-tetrahedral assignment, and unrelated properties. | This is exact parity with `assignChiralTypesFrom3D` for the modeled state space. It does not claim `assignStereochemistryFrom3D`, 3D double-bond direction/E-Z assignment, CIP orchestration, or distinct-substituent validation. COSMolKit intentionally preserves caller state when RDKit throws after partial mutation. |
| Tetrahedral stereo geometry | Geometry parity | Ordered-ligand orientation and signed-volume behavior against RDKit ETKDG-generated chiral geometries. | Checked against seeded RDKit ETKDG geometry; it validates stereo geometry, not all conformer-generation behavior by itself. |
| Batch surfaces | Determinism and order parity | Parallel batch paths preserve record order and reproduce the same outputs as scalar paths for bounds, prepared drawing, SVG, graph features, and supported transformations. | Batch checks are value/order checks, not separate chemistry definitions. |
| Fingerprints | Exact only for enumerated branches | Morgan sparse-count/sparse-bit/hashed-count/explicit-bit/AdditionalOutput branches and MACCS raw-167/public-166 vectors are checked for exact RDKit equality. `RDKFingerprintMol`/topological and Avalon are not implemented and return structured unsupported errors. | Passing Morgan/MACCS tests do not imply coverage of unmodeled RDKit APIs or preparation branches. No approximate vector is returned; similarity correlation and 99.9% bit agreement are failures. |
| Substructure matching | Partially checked, unfinished | Molecule-query matching is checked against RDKit query-to-target atom mappings. | Direct SMARTS-query parity remains a separate unfinished surface. |
