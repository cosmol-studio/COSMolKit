# RDKit MOL2 Source Map And Gap Report

## Step 3 Current COSMolKit MOL2 Support Audit

Date: 2026-05-28

Scope:

- Current COSMolKit Rust core MOL2 parser support.
- Current public MOL2 API exposure.
- Current support-matrix feature declaration.
- Current test, fixture, and golden coverage.

Evidence commands executed during this step:

- `rg -n "MOL2|Mol2|mol2|Tripos|SYBYL" crates tests dev -g '!target'`
- `find crates/cosmolkit-core/src/io -maxdepth 2 -type f | sort`
- `sed -n '1,220p' crates/cosmolkit-core/src/io/mod.rs`
- `sed -n '1,220p' crates/cosmolkit-core/src/support.rs`
- `find tests -maxdepth 3 -type f | sort | rg -i 'mol2|tripos|sybyl'`
- `find tests -maxdepth 3 -type d | sort | rg -i 'mol2|tripos|sybyl'`
- `rg -n "MOL2_READ_FEATURE|mol2.read|read_mol2|mol_from_mol2|Mol2" crates/cosmolkit-core/src crates/cosmolkit/src python tests -g '!target'`

Findings:

- `crates/cosmolkit-core/src/io/` currently contains `bio.rs`, `gemmi_spacegroup_table.rs`, `mod.rs`, `molblock.rs`, `molfile.rs`, `pdb_molecule.rs`, `pdb_writer.rs`, `sdf.rs`, and `xyz.rs`; it does not contain `mol2.rs`.
- `crates/cosmolkit-core/src/io/mod.rs` does not declare or re-export a MOL2 module.
- `crates/cosmolkit-core/src/support.rs` does not define `MOL2_READ_FEATURE`, `mol2.read`, or any MOL2 support spec.
- `crates/cosmolkit-core/src/lib.rs`, `crates/cosmolkit/src`, and `python/` do not expose `read_mol2`, `Mol2`, or `mol_from_mol2` APIs.
- `tests/` currently has no MOL2, Tripos, or SYBYL fixture directory or golden file.
- Existing `mol2` matches outside `dev/rdkit_mol2_full_port_plan.md` are variable names such as `let mol2 = ...`; they are not MOL2 file-format implementation evidence.

Current support status:

```text
MOL2 parser module: absent
MOL2 public Rust API: absent
MOL2 public Python API: absent
MOL2 FeatureSpec: absent
MOL2 fixtures: absent
MOL2 golden parity data: absent
MOL2 source markers: absent
```

Step 3 conclusion:

COSMolKit currently has no implemented MOL2 parsing surface. The next MOL2 work must start from source-backed RDKit mapping and explicit unsupported feature declaration before exposing any public parser behavior.

## Step 5 RDKit `FileParsers.h` MOL2 Declaration Audit

Source audited:

- `third_party/rdkit/Code/GraphMol/FileParsers/FileParsers.h:257-363`

Evidence commands executed during this step:

- `sed -n '250,370p' third_party/rdkit/Code/GraphMol/FileParsers/FileParsers.h`
- `rg -n "Mol2Type|Mol2ParserParams|MolFromMol2DataStream|MolFromMol2Block|MolFromMol2File|Mol2FileToMol|Mol2DataStreamToMol|Mol2BlockToMol" third_party/rdkit/Code/GraphMol/FileParsers/FileParsers.h`

RDKit declaration surface:

```text
Mol2Type::CORINA = 0
Mol2ParserParams::sanitize = true
Mol2ParserParams::removeHs = true
Mol2ParserParams::variant = Mol2Type::CORINA
Mol2ParserParams::cleanupSubstructures = true
MolFromMol2DataStream(std::istream&, Mol2ParserParams)
MolFromMol2Block(std::string, Mol2ParserParams)
MolFromMol2File(std::string, Mol2ParserParams)
Mol2FileToMol(fName, sanitize, removeHs, variant, cleanupSubstructures)
Mol2DataStreamToMol(inStream, sanitize, removeHs, variant, cleanupSubstructures)
Mol2BlockToMol(molBlock, sanitize, removeHs, variant, cleanupSubstructures)
```

Rust API mapping:

| RDKit declaration | COSMolKit target |
|---|---|
| `Mol2Type::CORINA` | `Mol2Type::Corina` |
| `Mol2ParserParams` | `Mol2ReadParams` |
| `Mol2ParserParams::sanitize` | `Mol2ReadParams::sanitize` |
| `Mol2ParserParams::removeHs` | `Mol2ReadParams::remove_hs` |
| `Mol2ParserParams::variant` | `Mol2ReadParams::variant` |
| `Mol2ParserParams::cleanupSubstructures` | `Mol2ReadParams::cleanup_substructures` |
| `MolFromMol2DataStream(std::istream&, params)` | `mol_from_mol2_data_stream_like_rdkit(input, params)` internal parser over buffered text state |
| `MolFromMol2Block(std::string, params)` | `mol_from_mol2_block_like_rdkit(text, params)` internal parser and `read_mol2_from_str_with_params(text, params)` public convenience |
| `MolFromMol2File(std::string, params)` | `mol_from_mol2_file_like_rdkit(path, params)` internal parser and `read_mol2_file_with_params(path, params)` public convenience |
| `Mol2BlockToMol(molBlock, sanitize, removeHs, variant, cleanupSubstructures)` | `read_mol2_from_str(text)` and `read_mol2_from_str_with_params(text, params)` |
| `Mol2FileToMol(fName, sanitize, removeHs, variant, cleanupSubstructures)` | `read_mol2_file(path)` and `read_mol2_file_with_params(path, params)` |
| `Mol2DataStreamToMol(...)` | No public streaming API in initial Rust surface; behavior is reproduced in the internal data-stream parser used by block and file entrypoints. |

Default policy required by source:

```text
sanitize: true
remove_hs: true
variant: Corina
cleanup_substructures: true
```

Step 5 conclusion:

The public Rust surface can mirror RDKit block and file convenience APIs while keeping stream parsing as a private implementation detail. The parameter defaults must match RDKit exactly before any parser behavior is exposed.

## Step 7 RDKit `Mol2FileParser.cpp` Function-Body Audit

Source audited:

- `third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp`

Evidence commands executed during this step:

- `rg -n "^(void|unsigned int|bool|Atom \\*|Bond \\*|std::unique_ptr<RWMol>) [A-Za-z0-9_]+|^std::unique_ptr<RWMol> MolFromMol2" third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp`
- `sed -n '52,180p' third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp`
- `sed -n '180,360p' third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp`
- `sed -n '360,528p' third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp`
- `sed -n '528,812p' third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp`
- `sed -n '820,1030p' third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp`
- `rg -n "MolOps::|SanitizeFlags|removeHs|assignChiralTypesFrom3D|detectBondStereochemistry|cleanUpMol2Substructures|guessFormalCharges|readFormalChargesFromAttr|ParseMol2" third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp`

Function-by-function source map:

| RDKit function | Source lines | COSMolKit target | Required behavior |
|---|---:|---|---|
| `fixNitroSubstructureAndCharge(RWMol*, unsigned int)` | 55-75 | `fix_nitro_substructure_and_charge_like_rdkit` | Count double-bond oxygen neighbors on a nitrogen, convert the last matching double bond to single, assign nitrogen `+1`, and assign that oxygen `-1` only when exactly two double-bond oxygen neighbors exist. |
| `readFormalChargesFromAttr(std::istream*, RWMol*)` | 77-137 | `read_formal_charges_from_attr_like_rdkit` | Parse `UNITY_ATOM_ATTR` records, accept only `AtomExpr` values without `=`, assign one-based atom formal charges, and stop on EOF, blank line, section line, or comment line. |
| `guessFormalCharges(RWMol*)` | 139-264 | `guess_formal_charges_like_rdkit` | Guess charges for non-carbon, non-query, zero-charge atoms from explicit valence, aromatic bond count, ring membership, Tripos atom type, periodic valence list, outer electron count, aromatic charge clamping, and nitro repair handoff. |
| `chkNoHNeighbNOx(RWMol*, ROMol::ADJ_ITER, int&)` | 267-284 | `check_no_h_neighbors_n_oxide_like_rdkit` | Count hydrogen neighbors of a candidate nitrogen and set `toModIdx` when it has a terminal oxygen neighbor used by the N-oxide precedence path. |
| `cleanUpMol2Substructures(RWMol*)` | 287-526 | `clean_up_mol2_substructures_like_rdkit` | Apply RDKit Tripos cleanup for `N.4`, `O.co2` phosphate/carboxylate/sulfonate cases, and `C.cat` two- and three-nitrogen guanidinium/amidine cases with `isFixed` ordering. |
| `ParseMol2FileAtomLine(std::string, Point3D&)` | 528-643 | `parse_mol2_atom_line_like_rdkit` | Tokenize atom lines, parse coordinates, store `_TriposAtomName`, store `_TriposAtomType`, map SYBYL symbols, remove `LP`, create ANY/Du/HEV/HET/HAL query atoms, set no implicit hydrogens, skip subst fields, and preserve `_TriposPartialCharge`. |
| `ParseMol2FileBondLine(std::string, INT_VECT const&)` | 646-720 | `parse_mol2_bond_line_like_rdkit` | Tokenize bond lines, convert one-based atom indices through `idxCorresp`, skip LP endpoints, map bond types `1`, `am`, `2`, `3`, `ar`, `du`, and `un`, reject index mismatch, and ignore unsupported strings including `nc`. |
| `ParseMol2AtomBlock(std::istream*, RWMol*, unsigned int, INT_VECT&)` | 723-776 | `parse_mol2_atom_block_like_rdkit` | Read exactly `nAtoms` atom lines, preserve LP-skipped index correspondence, collect 3D coordinates, detect explicit hydrogens, warn when absent, and add a 3D conformer sized to non-LP atoms. |
| `ParseMol2BondBlock(std::istream*, RWMol*, unsigned int, INT_VECT const&)` | 779-810 | `parse_mol2_bond_block_like_rdkit` | Read exactly `nBonds` bond lines, skip null bonds, set aromatic flags on aromatic bonds and endpoint atoms, and enforce the resulting bond-count postcondition. |
| `MolFromMol2DataStream(std::istream&, Mol2ParserParams const&)` | 820-997 | `mol_from_mol2_data_stream_like_rdkit`, `scan_mol2_sections_like_rdkit`, `parse_mol2_molecule_header_like_rdkit`, `assign_chiral_types_from_3d_for_mol2_like_rdkit`, `sanitize_mol2_molecule_like_rdkit`, `finish_mol2_read_like_rdkit` | Scan sections, reject missing required blocks, parse name/counts/charge type, parse atom and bond blocks, choose formal charge source, run 3D chirality before sanitization, and apply RDKit MOL2-specific sanitize/remove-H/stereo ordering. |
| `MolFromMol2Block(std::string const&, Mol2ParserParams const&)` | 1000-1005 | `mol_from_mol2_block_like_rdkit` | Wrap the input string in a stream and delegate to the data-stream parser. |
| `MolFromMol2File(std::string const&, Mol2ParserParams const&)` | 1011-1028 | `mol_from_mol2_file_like_rdkit` | Open file in binary mode, report bad input files, return parsed data for non-empty streams, and return null for immediate EOF. |

Cross-function dependencies:

- `MolFromMol2DataStream` is the required top-level semantic source for ordering; helper behavior must not be implemented independently of that ordering.
- `guessFormalCharges` depends on `fixNitroSubstructureAndCharge` and ring perception semantics.
- `cleanUpMol2Substructures` depends on `chkNoHNeighbNOx` and ring perception semantics.
- `ParseMol2AtomBlock` depends on `ParseMol2FileAtomLine`.
- `ParseMol2BondBlock` depends on `ParseMol2FileBondLine`.

Step 7 conclusion:

The MOL2 port must be implemented as a function-level RDKit reproduction. The source functions are available locally, and each future Rust function must carry copied source markers in the corresponding Rust body before behavior is marked reproduced.

## Step 9 RDKit `testMol2ToMol.cpp` Parity Fixture Audit

Source audited:

- `third_party/rdkit/Code/GraphMol/FileParsers/testMol2ToMol.cpp`
- `third_party/rdkit/Code/GraphMol/FileParsers/test_data/*.mol2`

Evidence commands executed during this step:

- `rg -n "testMol2|Mol2FileToMol|MolFromMol2File|test_data/.*\\.mol2|BOOST_LOG|TEST_ASSERT|Mol2Type|cleanupSubstructures|EZ_mol2" third_party/rdkit/Code/GraphMol/FileParsers/testMol2ToMol.cpp`
- `sed -n '1,430p' third_party/rdkit/Code/GraphMol/FileParsers/testMol2ToMol.cpp`
- `find third_party/rdkit/Code/GraphMol/FileParsers/test_data -maxdepth 1 -type f -name '*.mol2' -printf '%f\\n' | sort`

RDKit tested fixture matrix:

| Fixture | RDKit test scope |
|---|---|
| `nonExistFile.mol2` | Bad-file error path. |
| `pyrazole_pyridine.mol2` | Basic parse, atom count, one 3D conformer, and first atom coordinate regression. |
| `benzene.mol2` | Basic aromatic MOL2 parse and atom count. |
| `mol_noatoms.mol2` | File parse exception for zero atoms. |
| `mol_nomol.mol2` | File parse exception for missing `MOLECULE` block. |
| `lonePairMol.mol2` | `LP` atom removal and resulting atom/bond counts. |
| `symmetricGuanidine.mol2` | Unsanitized guanidine charge assignment. |
| `highlySymmetricGuanidine.mol2` | Sanitized guanidine charge assignment across repeated symmetric groups. |
| `Noxide.mol2` | N-oxide cleanup charges. |
| `fusedRing.mol2` | Aromatic fused-ring charge-guess exclusions. |
| `pyridiniumPhenyl.mol2` | Pyridinium charge assignment. |
| `sulfonAmide.mol2` | Sulfonamide charge retention. |
| `chargedAmidineRWH.mol2` | Charged amidine charge assignment with hydrogen precedence. |
| `chargedAmidineEC.mol2` | Charged amidine charge assignment with electron/cleanup precedence. |
| `chargedAmidine.mol2` | Charged amidine charge assignment baseline. |
| `dbtranslateCharged.mol2` | `UNITY_ATOM_ATTR` or dbtranslate charge behavior yielding charged atom. |
| `dbtranslateUncharged.mol2` | dbtranslate charge behavior yielding uncharged atom. |
| `dbtranslateUnchargedRing.mol2` | dbtranslate aromatic-ring no-charge behavior. |
| `badSubstPyridine.mol2` | Substituted aromatic group case accepted by test even when sanitize behavior is noted as problematic. |
| `Issue3399798.mol2` | 3D chirality remains unspecified for selected atoms. |
| `Issue3399798.2.mol2` | 3D chirality is assigned for one selected atom and remains unspecified for another. |
| `EZ_mol2_issue114.mol2` | Double-bond stereo assignment, including legacy and non-legacy stereo expectations. |
| `github438_1.mol2` | Metal MOL2 formal charge `+1`. |
| `github438_2.mol2` | Metal MOL2 formal charge `+2`. |
| `3505.mol2` | `cleanupSubstructures=true` versus `false` changes bond type and formal charge. |

Additional local RDKit `.mol2` fixtures present but not directly referenced by `testMol2ToMol.cpp`:

```text
Canion.mol2
Sulfonate.mol2
```

COSMolKit parity coverage required from this audit:

- Error parity for missing files, missing molecule section, and zero atoms.
- Topology parity for atom count, bond count, LP removal, aromatic flags, and bond typing.
- Coordinate parity for parsed 3D conformers.
- Formal-charge parity for cleanup, charge guessing, `UNITY_ATOM_ATTR`, metal cases, and cleanup toggling.
- Chirality and double-bond stereo parity for the RDKit issue fixtures.

Step 9 conclusion:

The initial COSMolKit MOL2 golden should be generated from the RDKit test fixture set rather than from ad hoc examples. The current plan's unit and parity steps cover the fixture categories found in RDKit's own MOL2 parser tests.
