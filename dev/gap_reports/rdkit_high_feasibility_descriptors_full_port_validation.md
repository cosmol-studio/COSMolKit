# RDKit High-Feasibility Descriptor Full-Port Validation

## Validation Boundary

- Reference: RDKit `2026.03.1`, revision
  `351f8f378f8ad6bbd517980c38896e66bf907af8`.
- Families: connectivity/Chi, Hall-Kier/Kappa/Phi, Lipinski count
  extensions, MQN, LabuteASA, SlogP-VSA, and SMR-VSA.
- Numeric policy: every finite floating-point leaf is compared by its exact
  IEEE-754 binary64 bits. Integer counts, all 42 MQN entries, every VSA bin,
  every atom contribution, and every cache-sequence result are separate
  observations.
- The matrix includes orders 0 through 6, fixed Chi projections, cold/warm/
  forced cache branches, both Labute hydrogen branches, default and custom
  VSA bins, and ordered contribution arrays.

## Complete Evidence

| Corpus | Input records | Mutually parseable records | Paired rejections | Exact descriptor observations | Mismatches |
|---|---:|---:|---:|---:|---:|
| Focused behavioral fixtures | 69 | 69 | 0 | 17,265 | 0 |
| `smiles_small` | 152 | 140 | 12 | 40,224 | 0 |
| Maintained 5,000-row corpus | 5,000 | 5,000 | 0 | 1,588,749 | 0 |
| ChEMBL 37 | 2,897,819 | 2,897,804 | 15 | 903,237,331 | 0 |

The paired rejection rows are retained parse-disposition checks. They do not
produce descriptor values and are therefore not included in the descriptor
observation column.

## ChEMBL 37 Execution

The repository-owned streaming phase processed all 128 deterministic corpus
shards and satisfied its configured record boundary.

| Field | Result |
|---|---:|
| Source records processed | 2,897,819 |
| Successful shard tasks | 128 / 128 |
| Mutually parseable records | 2,897,804 |
| Paired parser rejections | 15 |
| High-feasibility descriptor observations | 903,237,331 |
| Blocking mismatches | 0 |
| Retained finding examples | 0 |
| Complete | Yes |
| Passed | Yes |

The fixed ChEMBL source SHA-256 is
`ea6181ce8dc7af41974e35b92e1febb0c9dcbe2c62f7ccc4a5d983ac19f696e7`.
The complete local run manifest and aggregate are retained outside tracked
data under:

```text
tmp/parity-audit/chembl37-high-feasibility-descriptors-20260827/manifest.json
tmp/parity-audit/chembl37-high-feasibility-descriptors-20260827/descriptors/aggregate.json
```

The ChEMBL observation total is the sum of only
`match.descriptor.high_feasibility.*` leaves. It excludes the older descriptor
surface already present in the same phase, parser counters, and aggregate
status counters.

### Current-Code Stereo-Count Closure

The complete descriptor run above was produced before the later shared
potential-stereo implementation landed. That implementation is read by only
two descriptors in this selected family, so both affected fields were rerun
over every ChEMBL 37 shard against the current worktree:

| Field | Result |
|---|---:|
| Source records processed | 2,897,819 |
| Successful shard tasks | 128 / 128 |
| Mutually parseable records | 2,897,804 |
| Paired parser rejections | 15 |
| `num_atom_stereo_centers` exact comparisons | 2,897,804 |
| `num_unspecified_atom_stereo_centers` exact comparisons | 2,897,804 |
| Blocking mismatches | 0 |
| Retained finding examples | 0 |

The rerun manifest identifies Git `494fe2812d7def1ca3d99d0c6d506f8c3310e2d1`,
the complete staged descriptor diff, the installed extension binary, the
temporary targeted runner, and the unchanged ChEMBL corpus by SHA-256. Its
5,795,608 comparisons reconfirm fields already included in the
903,237,331-observation descriptor total and are therefore not added to that
total a second time. The ignored run evidence is retained under:

```text
tmp/parity-audit/chembl37-descriptor-stereo-counts-current/manifest.json
tmp/parity-audit/chembl37-descriptor-stereo-counts-current/aggregate.json
```

## Reproduction Commands

```bash
.venv/bin/python tools/testdata/rdkit/generate_all.py \
  --python .venv/bin/python \
  --profile smiles_5000 \
  --suite descriptors \
  --jobs 16 \
  --clean

COSMOLKIT_PARITY_PROFILE=smiles_5000 \
cargo test -p cosmolkit-core --release --features op-contracts-strict \
  --test rdkit_high_feasibility_descriptor_parity \
  active_corpus_high_feasibility_descriptors_match_pinned_rdkit_exactly \
  -- --exact --nocapture

.venv/bin/python dev/tools/chembl_parity/run.py \
  --corpus-dir /path/to/chembl37-128-v1 \
  --run-dir /path/to/run \
  --workers 96 \
  --phase descriptors
```

Generated expected data and complete ChEMBL run artifacts remain ignored or
external. The generator, schema, focused input fixtures, comparison code, and
streaming harness are repository-owned.

## Performance Evidence

The isolated-worker benchmark verifies exact output before timing and compares
the same public descriptor call on molecules prepared outside the timer. Cold
measurements use molecules that have not entered the selected descriptor;
warm measurements repeat the descriptor on one warmed molecule. Results are
machine observations, not substitutes for the source-level complexity audit.

| Family | COSMolKit/RDKit cold ratio range | COSMolKit/RDKit warm ratio range | Disposition |
|---|---:|---:|---|
| Chi order 2/6 | 0.58–1.23x | 0.56–1.80x | Source-equivalent path and typed delta-cache route; comparable complexity and call cost. |
| LabuteASA | 0.22–0.41x | 0.52–0.94x | Source-equivalent arithmetic and retained typed cache. |
| Ring-complexity counts | 0.44–0.72x | 0.38–0.64x | Shared initialized ring state is retained and reused. |
| MQN | 0.46–0.99x | 0.46–1.57x | Shared ring state, hydrogen cache, and parsed rotatable queries are retained. |
| SlogP-VSA | 1.19–2.04x | 0.55–0.86x | Cold atom typing is source-equivalent; warm typed Crippen/Labute caches are retained. |
| SMR-VSA | 0.77–2.02x | 0.55–0.85x | Cold atom typing is source-equivalent; warm typed Crippen/Labute caches are retained. |

The benchmark output is retained at:

```text
tmp/high-feasibility-descriptor-benchmark.json
```

All 35 case/family preflights matched exactly. The measurement includes the
shared ring-state provider, parsed Crippen and rotatable-query flyweights,
typed Crippen/Labute/connectivity caches, and Python descriptor wrappers that
borrow the retained `Molecule` through `PyRef`. Cold paths remain linear in the
same graph, path, ring, or SMARTS work as the source; warm paths retain the
same observational state instead of writing it to an extracted clone. No
selected family has a remaining material performance- or complexity-axis gap.

## Source Marker Audit

Every source function selected by the plan has a line-level RDKit source
anchor. The behavioral axis is closed for the modeled descriptor state: the
focused matrix, maintained 5,000-row corpus, and complete ChEMBL 37 execution
all compare the ordered scalar, vector, contribution, custom-bin, and cache-
sequence outputs without a mismatch.

The performance axis is closed for the selected function and wrapper paths.
The final local review covered loop nesting, allocation, cache ownership,
query reuse, path enumeration, ring-state acquisition, and Python extraction.

| Source-function group | Behavioral axis | Performance-axis disposition |
|---|---|---|
| Connectivity, Hall-Kier, Kappa, Phi, Lipinski scalar counts, stereo counts, Labute arithmetic, and VSA bin assignment | Exact for the modeled state | `RDKit✔️✔️`: source-equivalent work, allocation shape, and cache access; final ratios are comparable. |
| `calcChiNv`, `calcChiNn`, and path-backed Kappa/Phi calls | Exact | `RDKit✔️✔️`: one shared source-order path kernel and retained typed delta vectors; final Chi ratios are 0.56–1.80x. |
| Lipinski ring-family projections, spiro, and bridgehead | Exact | `RDKit✔️✔️`: one provider borrows initialized SSSR-or-better state and computes the source-required state only on a cold molecule. |
| `calcMQNs` | Exact | `RDKit✔️✔️`: direct source MQN polarity rules plus shared O(1) total-H, ring-state, and default-rotor paths. |
| `getCrippenAtomContribs`, `calcSlogP_VSA`, and `calcSMR_VSA` | Exact | `RDKit✔️✔️`: parsed parameter queries and typed atom contributions are retained; warm ratios are 0.55–0.86x. |
| `getLabuteAtomContribs`, `calcLabuteASA`, and `_LabuteHelper` | Exact | `RDKit✔️✔️`: the typed atom/H/total cache and borrowed Python wrapper reproduce the source lifecycle; warm ratios are 0.52–0.94x. |

The remaining `RDKit✔️❌` block in `descriptors.rs` belongs to the pre-existing
`calc_crippen_descriptors(include_hs=true)` working-molecule branch. The selected
VSA functions call the atom-contribution core directly and never enter that
branch, so the marker remains accurate without limiting this port's two-axis
claim. No selected source block retains an unresolved marker.

## Architecture Audit

| Boundary | Result | Disposition |
|---|---|---|
| Core ownership | One connectivity module, one MQN core, one Labute core, one VSA binning loop, one Crippen contribution engine, one neutral path enumerator, and one authoritative Rb0 table were found. The existing aromatic-ring and rotatable-bond implementations are reused rather than copied. | Complete and enforced by architecture tests. |
| Public API paths | The crate-root exports delegate to the descriptor facade; the Python wrappers delegate to the crate-root functions; fixed VSA scalar functions index the vector core. No alternate public calculation path was found. | Complete and guarded as projection-only. |
| Mutable storage | Descriptor production code accepts `&Molecule` and writes only through typed observational computed-cache methods. No public mutable topology, coordinate, property, or derived-cache storage is exposed by these descriptor APIs. Direct derived-cache mutation found in `lipinski.rs` is confined to a module test fixture. | Complete. |
| Descriptor-local chemistry | Chi/Kappa use the shared source-reproduced subgraph kernel; SMARTS counts use the shared parser/matcher; MQN reuses shared HBA/HBD, total-H, ring, and rotatable-bond behavior; Hall-Kier/Labute use the authoritative periodic table. No second valence, SMARTS, path, or Rb0 implementation was found. | Complete and guarded against family-local replacements. |
| Ring-state acquisition | `descriptor_ring_info()` is the single selected-family provider. It borrows initialized SSSR-or-better state and owns only cold source-required perception output. | Complete. |
| Crippen state | `rdkit_crippen_params()` retains parsed queries in one `OnceLock`; `rdkit_crippen_atom_contribs()` uses a typed topology-invalidated contribution cache. | Complete for the selected VSA path. |
| Labute state | One typed cache stores atom, hydrogen, and total ASA contributions; Python wrappers borrow the retained molecule. | Complete. |
| Static architecture guard | The guard checks the actual VSA core spelling, sole family owners, facade delegation, scalar projection, and borrowed Python descriptor extraction. | Complete. |

The architecture is coherent at the public API, chemistry-ownership, cache,
and Python-wrapper layers.

## Final Signoff

**Decision: complete and accepted.** Every completion criterion in the
execution plan is satisfied for the declared descriptor boundary.

| Completion criterion | Final evidence | Result |
|---|---|---|
| Source-ledger closure | The completed source inventory assigns every selected function to one concrete implementation or reuse proof. Focused tests cover the source branches; no selected path contains an approximate, placeholder, heuristic, or silent fallback. | Pass |
| Source reproduction markers | Every selected C++ or Python source function has an in-function source anchor and a closed behavioral/performance disposition. The only remaining `RDKit✔️❌` descriptor block is the pre-existing, excluded `calc_crippen_descriptors(include_hs=true)` working-molecule branch. | Pass |
| Single-core architecture | Static guards confirm one path enumerator, Rb0 table, selected-family ring provider, standard HBA/HBD and rotor cores, Crippen contribution engine, MQN engine, Labute engine, and VSA binning engine. The SMARTS inventory now preserves production code after interleaved test modules and explicitly tracks the one descriptor flyweight helper. | Pass |
| Cache and ownership semantics | Focused cold, warm, forced, clone, invalidation, mixed-order, and concurrent-read tests pass for connectivity, Crippen, Labute, ring, and Python `PyRef` paths. | Pass |
| Exact parity evidence | Focused, `smiles_small`, maintained 5,000-row, and complete ChEMBL 37 gates report zero mismatch across finite float bits, counts, contribution arrays, 42 MQN entries, VSA bins, custom boundaries, and cache sequences. | Pass |
| Baseline disposition | The five pre-port failures were recorded as stale expected-data identity failures. Expected data were regenerated only through the repository entrypoint; the final strict core and workspace runs pass without weakening manifest or checksum validation. | Pass |
| Public surface and documentation | Rust exports, Python bindings, generated stub, the descriptor example, README, Sphinx docs, support metadata, `VALIDATION.md`, and the unreleased changelog describe the validated surface. The example executes successfully against the installed extension. | Pass |
| Final validation gates | `cargo fmt --all -- --check`; strict release core; strict release workspace; `maturin develop`; Python tests (`653 passed`, `37 skipped`); and Sphinx HTML complete successfully. Basedpyright reports `0 errors`; its repository-wide warning set remains a separate non-clean gate. | Pass for this descriptor boundary; project-wide warning cleanup remains open. |

Generated expected data and benchmark/ChEMBL artifacts remain in their
policy-defined ignored or external locations. No Git operation was used while
executing this plan.
