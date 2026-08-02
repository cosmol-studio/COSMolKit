# RDKit MolBlock Reader Source Map

Step: `dev/archive/plans/rdkit_molblock_sdf_full_port_plan.md` Step 181 final update.

Scope audited:

- RDKit `third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp::MolFromMolDataStream`
- RDKit `third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp::MolFromMolBlock`
- RDKit `third_party/rdkit/Code/GraphMol/FileParsers/MolFileParser.cpp::MolFromMolFile`
- RDKit `third_party/rdkit/Code/GraphMol/FileParsers/FileParsers.h::MolBlockToMol`
- RDKit `third_party/rdkit/Code/GraphMol/FileParsers/FileParsers.h::MolFileToMol`
- COSMolKit `crates/cosmolkit-core/src/io/sdf.rs::SdfReadParams`
- COSMolKit `crates/cosmolkit-core/src/io/sdf.rs::mol_from_mol_block`
- COSMolKit `crates/cosmolkit-core/src/io/sdf.rs::mol_from_mol_data_stream`
- COSMolKit `crates/cosmolkit-core/src/io/molfile.rs::read_mol_record_from_str_with_params`
- COSMolKit `crates/cosmolkit-core/src/io/molfile.rs::read_mol_file_with_params`

## Final Scope Statement

The public MolBlock/molfile reader scope covered by this plan has no remaining
known behavioral gap against the mapped RDKit reader entry points for the
modeled parameter space:

```text
sanitize in {true,false}
removeHs in {true,false}
strictParsing in {true,false}
public MolBlock and molfile inputs, including V2000 and V3000 records
```

RDKit `MolFileParserParams::parsingSCSRMol` remains outside COSMolKit's public
MolBlock/SDF reader model. It is not silently approximated: no public
`SdfReadParams` field claims SCSR parsing support, and this report does not
count SCSR-only `expectMEND=false` / macro-atom behavior as supported.

Some copied source markers in the large shared CTAB reader code still use
second-axis `❗` or `❌` markers for performance/complexity status in adjacent
SDF supplier or parser internals. Those are not remaining MolBlock reader
behavior gaps. They record unresolved or known non-equivalent implementation
shape such as buffered SDF record handling or indexed SDF rescans.

## Source Parameter Map

| RDKit field or wrapper argument | COSMolKit target | Final status |
|---|---|---|
| `MolFileParserParams::sanitize` | `SdfReadParams::sanitize` | Source-backed and publicly propagated through Rust and Python MolBlock/molfile readers. |
| `MolFileParserParams::removeHs` | `SdfReadParams::remove_hs` | Source-backed and publicly propagated through Rust and Python MolBlock/molfile readers. |
| `MolFileParserParams::strictParsing` | `SdfReadParams::strict_parsing` | Source-backed for counts-line, CTAB-version, V2000, and V3000 reader behavior in the modeled scope. |
| `MolFileParserParams::expandAttachmentPoints` | `SdfReadParams::expand_attachment_points` | Source-backed through finalization; detailed finalization coverage is tracked in `rdkit_molblock_finalization_source_map.md`. |
| `MolFileParserParams::parsingSCSRMol` | No public `SdfReadParams` field | Out of public scope and not claimed supported. |
| v1 `MolBlockToMol(... sanitize, removeHs, strictParsing)` | `read_mol_record_from_str_with_params` and Python wrappers | Covered through the same source-backed parameter carrier and parity tests. |
| v1 `MolFileToMol(... sanitize, removeHs, strictParsing)` | `read_mol_file_with_params` and Python wrappers | Covered through the same source-backed parameter carrier and parity tests. |

## Function Map

| RDKit function | COSMolKit function | Final finding |
|---|---|---|
| `MolFromMolBlock` | `mol_from_mol_block`; `read_mol_record_from_str_with_params` | No known behavioral gap in modeled scope. RDKit source is copied in both Rust entry points. Both construct a stream-like reader and delegate to `mol_from_mol_data_stream`; unread trailing text after `M  END` is ignored like RDKit. |
| `MolFromMolFile` | `read_mol_file_with_params` | No known behavioral gap in modeled scope. RDKit source is copied in the Rust function. The implementation opens a file, checks immediate empty input, streams into `read_mol_data_stream_molecule_with_params`, and maps RDKit's null-molecule branch to a structured parse error because COSMolKit's public API returns `Result<MolFileRecord, SdfReadError>`, not nullable molecules. |
| `MolFromMolDataStream` header reads | `parse_mol_header` | No known behavioral gap in modeled scope. Name, info line, comments line, line counting, `_MolFileInfoLine`, `_MolFileComments`, and 3D label handling are source-backed and tested. |
| `MolFromMolDataStream` counts-line parse | `parse_counts_line` | No known behavioral gap in modeled scope. Fixed-width atom/bond counts, ignored optional field conversion failures, chiral flag, CTAB version selection, strict/non-strict CTAB version behavior, and V3000 counts are source-backed and tested. |
| `MolFromMolDataStream` V2000 dispatch | `parse_v2000_ctab` | No known reader-dispatch behavioral gap in modeled scope. Detailed V2000 atom, bond, property, query, SGroup, stereo, and finalization behavior is covered by the focused tests and parity corpus created later in the plan. |
| `MolFromMolDataStream` V3000 dispatch | `parse_v3000_ctab` | No known reader-dispatch behavioral gap in modeled public scope. V3000 dispatch, initial-count strictness, `END CTAB`, `M  END`, atom/bond/SGroup/stereo collection coverage, and non-strict recovery paths are source-backed and tested where modeled. SCSR-specific macro-atom parsing remains out of scope. |
| `MolFromMolDataStream` missing `M  END` handling | `mol_from_mol_data_stream`, `parse_v2000_ctab`, `parse_v3000_ctab` | No known behavioral gap in modeled scope. Missing completion returns structured parse errors with RDKit-aligned message categories and line behavior covered by focused tests. |
| `MolFromMolDataStream` finalization call | `finish_mol_processing` | No known reader-entry behavioral gap. Detailed post-CTAB finalization is owned by `rdkit_molblock_finalization_source_map.md`. |
| `FileParsers.h::MolBlockToMol` | Public MolBlock wrappers | No known behavioral gap in modeled scope. Wrapper parameter propagation is covered by Rust and Python parity/sanity tests. |
| `FileParsers.h::MolFileToMol` | Public molfile wrappers | No known behavioral gap in modeled scope. Wrapper parameter propagation is covered by Rust and Python parity/sanity tests. |

## Evidence

Source anchors:

- `crates/cosmolkit-core/src/io/molfile.rs::read_mol_file_with_params`
- `crates/cosmolkit-core/src/io/molfile.rs::read_mol_record_from_str_with_params`
- `crates/cosmolkit-core/src/io/sdf.rs::mol_from_mol_block`
- `crates/cosmolkit-core/src/io/sdf.rs::mol_from_mol_data_stream`
- `crates/cosmolkit-core/src/io/sdf.rs::parse_mol_header`
- `crates/cosmolkit-core/src/io/sdf.rs::parse_counts_line`

Focused tests:

- `crates/cosmolkit-core/src/io/molfile.rs` covers MolBlock unread trailing text, empty input, missing counts line, `strictParsing`, `sanitize`, `removeHs`, file-open failure, file-read failure, empty file, and single-record molfile behavior.
- `crates/cosmolkit-core/src/io/sdf.rs` covers `mol_from_mol_data_stream` header reads, counts-line parsing, strict/non-strict CTAB-version handling, V2000 and V3000 dispatch, initial V3000 counts, missing `M  END`, and completion behavior.
- `crates/cosmolkit-core/tests/rdkit_molfile_read_parity.rs` compares generated RDKit MolBlock/molfile golden rows against COSMolKit atom fields, bond fields, conformers, molecule properties, SGroups, SMILES, delayed sanitize, delayed hydrogen removal, and error rows.

Validation performed after the full port:

```bash
cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_molfile_read_parity -- --nocapture
cargo test -p cosmolkit-core --features op-contracts-strict
cargo test --workspace --features cosmolkit-core/op-contracts-strict
.venv/bin/pytest
```

The final full validation completed successfully:

```text
cargo test -p cosmolkit-core --features op-contracts-strict: passed
cargo test --workspace --features cosmolkit-core/op-contracts-strict: passed
.venv/bin/pytest: 491 passed, 37 skipped
```

## Residual Non-Claims

- COSMolKit does not claim public SCSR molfile parsing parity.
- COSMolKit does not expose RDKit warning logs as a public warning stream; strict/non-strict observable success or structured error behavior is what is covered.
- Performance/complexity markers in nearby SDF supplier indexing code are not reader behavior gaps, but they should not be batch-upgraded without a separate local performance review.
