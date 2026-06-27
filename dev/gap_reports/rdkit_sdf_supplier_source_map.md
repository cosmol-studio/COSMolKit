# RDKit SDF Supplier Source Map

Step: `dev/rdkit_molblock_sdf_full_port_plan.md` Step 185 final update.

Scope audited:

- RDKit `third_party/rdkit/Code/GraphMol/FileParsers/ForwardSDMolSupplier.cpp::ForwardSDMolSupplier::next`
- RDKit `third_party/rdkit/Code/GraphMol/FileParsers/ForwardSDMolSupplier.cpp::ForwardSDMolSupplier::_next`
- RDKit `third_party/rdkit/Code/GraphMol/FileParsers/ForwardSDMolSupplier.cpp::ForwardSDMolSupplier::readMolProps`
- RDKit `third_party/rdkit/Code/GraphMol/FileParsers/ForwardSDMolSupplier.cpp::ForwardSDMolSupplier::checkForEnd`
- RDKit `third_party/rdkit/Code/GraphMol/FileParsers/SDMolSupplier.cpp::SDMolSupplier`
- RDKit `third_party/rdkit/Code/GraphMol/FileParsers/SDMolSupplier.cpp::SDMolSupplier::buildIndexTo`
- RDKit `third_party/rdkit/Code/GraphMol/FileParsers/SDMolSupplier.cpp::SDMolSupplier::operator[]`
- RDKit `third_party/rdkit/Code/GraphMol/FileParsers/MolSupplier.v1API.h`
- RDKit `third_party/rdkit/Code/GraphMol/FileParsers/MultithreadedSDMolSupplier.cpp::processMoleculeRecord`
- COSMolKit `crates/cosmolkit-core/src/io/sdf.rs`
- COSMolKit `crates/cosmolkit-core/src/properties/batch.rs`

## Final Scope Statement

The public SDF supplier scope covered by this plan has no remaining known
behavioral gap against the mapped RDKit supplier functions for the modeled
parameter space:

```text
sanitize in {true,false}
removeHs in {true,false}
strictParsing in {true,false}
processPropertyLists in {true,false}
ForwardSDMolSupplier-like streaming reads
SDMolSupplier-like indexed reads
batch reads over reader and dataset paths
valid, invalid, malformed, missing-delimiter, and EOF records
```

This is a behavioral finding for the modeled public supplier surface. It is
not a performance-equivalence finding. The indexed dataset and forward-reader
paths intentionally retain `RDKit✔️❌` second-axis markers because COSMolKit
buffers complete records and `SdfDataset::record_with_params` reopens/rescans
instead of seeking RDKit's cached stream offsets. Those are documented
performance/complexity gaps, not known public behavior gaps.

The RDKit `_next()` generic C++ `catch (...)` block remains a residual
non-claim. COSMolKit exposes structured `Result` errors and does not claim to
recover from arbitrary Rust panics as if they were RDKit catch-all exceptions.
Modeled parse, sanitize, data-field, missing-delimiter, and EOF recovery paths
are source-backed and tested.

## Parameter Propagation Map

| RDKit source | COSMolKit target | Final finding |
|---|---|---|
| `ForwardSDMolSupplier(... MolFileParserParams)` stores `d_params` | `SdfReader::with_params` stores `SdfReadParams` | No known behavioral gap. Streaming reads carry `sanitize`, `remove_hs`, `strict_parsing`, `coordinate_mode`, and `process_property_lists` through record parsing. |
| v1 `ForwardSDMolSupplier(... sanitize, removeHs, strictParsing)` | Rust `SdfReader`; Python `PySdfReader::open` | No known behavioral gap in modeled public scope. Python propagation was completed and validated by Python tests and full pytest. |
| v1 `SDMolSupplier(file, sanitize, removeHs, strictParsing)` | `SdfDataset::open_with_params`; `SdfDataset::record_with_params`; Python `SdfDataset.open` | No known behavioral gap. Open parameters are retained, record-level override parameters are supported, and out-of-range behavior matches RDKit's error category/message. |
| `MultithreadedSDMolSupplier::d_parseParams` | `MoleculeBatch::read_sdf_records_from_reader_with_params_and_options`; dataset batch APIs | No known behavioral gap in modeled scope. Batch reads now route through `SdfReadParams` and preserve recovered no-molecule records as error records unless the reader is truly at end. |
| `df_processPropertyLists` | `SdfReadParams::process_property_lists` | No known behavioral gap in modeled scope. Atom and bond property lists are processed when enabled and preserved only as raw data fields when disabled. |

## Function Map

| RDKit function | COSMolKit function | Final finding |
|---|---|---|
| `ForwardSDMolSupplier::next` | `SdfReader::next_record`; `read_next_sdf_record` | No known behavioral gap in modeled scope. EOF, empty input, `df_end`-like state, and no-record return behavior are represented by `Ok(None)` plus `SdfReader::is_end()`. |
| `ForwardSDMolSupplier::_next` normal record path | `read_next_sdf_record`; `extract_next_forward_sdf_record`; `parse_forward_sdf_record_text` | No known behavioral gap in modeled scope. COSMolKit consumes the delimiter, updates line/byte/index state, returns the first record before the second, and handles final records without trailing newline. |
| `_next` parse-exception recovery | `parse_forward_sdf_record_text`; `extract_next_forward_sdf_record` | No known behavioral gap in modeled scope. A malformed mol block returns `Ok(None)` for the current record, consumes through the next `$$$$`, keeps the stream aligned, and allows the next valid record to be read. |
| `_next` sanitize-exception recovery | `parse_forward_sdf_record_text` | No known behavioral gap in modeled scope. Sanitization failures return `Ok(None)` for the current record and leave the reader positioned at the following record, matching RDKit supplier recovery behavior. |
| `_next` malformed data-field recovery | `parse_forward_sdf_record_text`; `parse_sdf_data_fields` | No known behavioral gap in modeled scope. Forward supplier behavior keeps the molecule and drops malformed data fields when RDKit would recover at supplier level. Direct `parse_sdf_record_text` keeps strict structured errors for record-local parsing. |
| `_next` missing delimiter at EOF | `extract_next_forward_sdf_record`; `read_next_sdf_record` | No known behavioral gap in modeled scope. A final record without `$$$$` is returned and marks the reader ended. Empty input returns no record and marks EOF. |
| `_next` null molecule before EOF | `read_next_sdf_record`; batch reader no-molecule handling | No known behavioral gap in modeled scope. `Ok(None)` before end is preserved as a per-record error in `MoleculeBatch`, while `Ok(None)` with `is_end()` terminates iteration. |
| `_next` generic `catch (...)` | No public panic-recovery claim | Residual non-claim. COSMolKit does not model arbitrary C++ catch-all exception recovery as public Rust behavior. |
| `ForwardSDMolSupplier::readMolProps` | `parse_sdf_data_fields`; `process_sdf_property_lists` | No known behavioral gap in modeled scope. Named fields, multiline values, indented blank continuation lines, terminal `\r` stripping, repeated fields, invalid/empty headers, strict spurious-data errors, non-strict truncation/ignore behavior, and missing-delimiter behavior are source-backed and tested. |
| `FileParserUtils::processMolPropertyList` | `process_sdf_property_list`; property-list helpers | No known behavioral gap in public modeled scope. Prefix matching, atom/bond targeting, missing-value handling, bool/int/double/string token boundaries, and `process_property_lists=false` behavior match RDKit parity tests. COSMolKit stores public atom/bond property payloads as strings rather than RDKit typed `RDProps`; this is not a known public behavior gap for exposed comparison data. |
| `ForwardSDMolSupplier::checkForEnd` | raw-record extraction and reader EOF state | No known behavioral gap in modeled scope. Empty EOF, final-record EOF, and records without trailing newline are covered. |
| `SDMolSupplier` constructor/open | `SdfDataset::open_with_params`; `open_sdf_dataset`; `build_sdf_index` | No known behavioral gap in modeled scope. Metadata construction records raw record index, byte offset, byte length, line offset, line length, and title, including invalid records. |
| `SDMolSupplier::buildIndexTo` | `build_sdf_index`; `extract_next_indexed_sdf_record` | No known behavioral gap in modeled scope. Behavior is source-backed for record discovery, offsets, invalid-record inclusion, trailing EOF, and delimiter handling. Complexity is intentionally marked worse than RDKit's chunked indexed seek path. |
| `SDMolSupplier::operator[]` | `SdfDataset::record`; `SdfDataset::record_with_params`; `read_indexed_sdf_record` | No known behavioral gap in modeled scope. Random access, out-of-range errors, invalid-record errors, open parameter reuse, record parameter override, and sliced batch reads are covered. Indexed access remains materially slower because it rescans to the requested record. |
| `MolSupplier.v1API.h` wrapper defaults | Rust and Python SDF reader/dataset/batch APIs | No known behavioral gap in modeled public scope. Parameter propagation was completed through Rust and Python wrappers and validated after stub/docs/example updates. |
| `MultithreadedSDMolSupplier::processMoleculeRecord` | `parse_sdf_record_text`; batch reader/dataset paths | No known behavioral gap in modeled scope. Record-local parse behavior and default batch output are covered by SDF parity tests. The implementation is sequential despite accepting `n_jobs`; this is a performance/concurrency non-claim, not a behavior gap in returned records. |

## Evidence

Source anchors:

- `crates/cosmolkit-core/src/io/sdf.rs::SdfReader::with_params`
- `crates/cosmolkit-core/src/io/sdf.rs::SdfReader::next_record`
- `crates/cosmolkit-core/src/io/sdf.rs::read_next_sdf_record`
- `crates/cosmolkit-core/src/io/sdf.rs::extract_next_forward_sdf_record`
- `crates/cosmolkit-core/src/io/sdf.rs::parse_forward_sdf_record_text`
- `crates/cosmolkit-core/src/io/sdf.rs::parse_sdf_record_text`
- `crates/cosmolkit-core/src/io/sdf.rs::parse_sdf_data_fields`
- `crates/cosmolkit-core/src/io/sdf.rs::process_sdf_property_lists`
- `crates/cosmolkit-core/src/io/sdf.rs::process_sdf_property_list`
- `crates/cosmolkit-core/src/io/sdf.rs::SdfDataset::open_with_params`
- `crates/cosmolkit-core/src/io/sdf.rs::build_sdf_index`
- `crates/cosmolkit-core/src/io/sdf.rs::read_indexed_sdf_record`
- `crates/cosmolkit-core/src/properties/batch.rs::read_sdf_records_from_reader_with_params_and_progress`
- `crates/cosmolkit-core/src/properties/batch.rs::read_sdf_records_from_reader_with_params_and_options`

Focused tests:

- `crates/cosmolkit-core/src/io/sdf.rs` covers forward supplier delimiter consumption, EOF state, empty EOF, parse recovery, sanitize recovery, malformed data-field recovery, missing delimiter at EOF, index metadata, random access, invalid indexed records, out-of-range indexed records, open-parameter reuse, parameter overrides, and sliced dataset batches.
- `crates/cosmolkit-core/src/io/sdf.rs` also covers SDF data fields and property lists: named fields, blank values, repeated fields, invalid headers, strict/non-strict spurious data, atom and bond property lists, disabled property-list processing, bool/token boundaries, and terminal carriage return handling.
- `crates/cosmolkit-core/tests/rdkit_sdf_read_parity.rs` compares generated RDKit `ForwardSDMolSupplier` and `SDMolSupplier` golden rows against direct string reads, `SdfReader::with_params`, `SdfDataset::record_with_params`, reader batch paths, dataset batch paths, data fields, property lists, atom/bond/conformer/finalization state, SMILES, and errors.
- `python/tests/test_python_api_sanity.py` and `python/tests/test_value_semantics.py` cover Python SDF reader parameter propagation and batch/dataset/reader behavior.

Validation performed after the full port:

```bash
cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_sdf_read_parity -- --nocapture
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

- COSMolKit does not claim performance parity for SDF indexed access or stream parsing. Current source markers correctly retain `RDKit✔️❌` where behavior is reproduced through buffering or rescanning.
- COSMolKit does not claim multithreaded supplier concurrency parity. `n_jobs` is accepted in batch APIs, but supplier record behavior is validated independently of parallel execution.
- COSMolKit does not claim RDKit warning-log stream parity; strict/non-strict observable success, recovered no-molecule records, structured errors, and preserved data are the covered public contract.
- COSMolKit does not claim recovery from arbitrary Rust panics as RDKit's C++ `catch (...)` branch. Modeled parse and sanitize failures recover through structured `Result` paths.
