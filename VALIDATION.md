# COSMolKit Validation

## Current status: COSMolKit 0.5.1-rc.2 / RDKit 2026.09.1

Current validation covers **180,000 ChEMBL 37 source records** across **894
completed task runs**, with **220,258,878 matching observations**, including
matching source-defined errors and verified corrections. The identified
conformer differences have passed bit-exact revalidation, and CIP node-limit
rejections now pair correctly. **No unresolved chemical mismatches remain in
these completed observations.**

### Reference and comparison rules

- Chemistry reference: RDKit **2026.09.1**, source revision
  `fece8caa860bdf6c9c82bdb9253f22a0cf87150c`.
- Corpus: ChEMBL 37, 2,897,819 available source records, processed in
  20,000-record batches through `parity-tests` prepare and Cargo comparison.
- Execution: release optimization with `cosmolkit/op-contracts-strict`.
- All numerical observations compare their encoded IEEE-754 bits, including
  descriptors, coordinates, force-field energies/gradients and alignment
  matrices. **No numerical tolerance is used.** Other outputs follow their
  declared exact representations; source-defined errors pair with concrete
  Rust error kinds.
- Counts combine completed measurements and verified corrections. Duplicate
  attempts and interrupted task groups are excluded. A comparison is one
  input/profile observation, not one additional distinct molecule.

### Completed ChEMBL coverage

| Result | Count |
|---|---:|
| Source records with completed task coverage | 180,000 |
| Records covered by all 100 scheduled corpus tasks | 60,000 |
| Records covered by 99 tasks, excluding incomplete tautomer preparation | 120,000 |
| Completed task runs | 894 |
| Exact matching values, including verified corrections | 220,015,543 |
| Matching source-defined errors, including corrected pairing | 243,335 |
| Unresolved mismatches | 0 |
| Reference timeout observations | 1,072 |
| COSMolKit timeout observations | 633 |
| Timeout observations on either side, not compared | 1,122 |
| Upstream crashes / undefined outcomes | 0 / 0 |

Coverage includes SMILES parsing and the 768-profile writer matrix; CIP,
valence and potential stereo; Morgan, MACCS, Avalon, topological, Layered,
Pattern, AtomPair and Torsion fingerprints and fingerprint operations;
descriptors and ring counts; sanitization, kekulization and hydrogen transforms;
fragments and scaffolds; substructure matching; tautomer enumeration; 2D
coordinates and drawing; bounds and conformer generation; alignment/RMSD;
MMFF/UFF availability, optimization and persistent force fields.

Tautomer coverage comprises 60,000 records. Its preparation in six other
batches exceeded the 128 MiB reference-record bound and remains incomplete.
Timeouts and unfinished groups are not counted as matches. Million-scale
ChEMBL coverage remains pending.

### Verified corrections and maintained-suite results

| Check | Current result | Evidence |
|---|---|---|
| ETKDGv3 single-/multi-conformer coordinates | **Pass — bit-exact** | Retained-reference reproduction matches all ten profiles, covering count 1 and count 2; all 4,233 instrumented checkpoints agree |
| CIP digraph node-limit rejection | **Pass — matching source error** | Retained-reference replay matches all four profiles; typed error pairing retains wrong-kind, wrong-limit and wrong-task negative checks |
| Standard full suite with `smiles_5000` | **Pass — zero mismatches** | 113 task reports; 6,139,696 matching values and 8,083 matching source errors; 44 timeout observations separately recorded |
| Post-correction CIP `smiles_5000` | **Pass — zero mismatches** | Four profiles: 16,918 matching values and 3,082 matching source errors; no timeouts |
| Post-correction conformer `smiles_5000` | **Pass — bit-exact, zero mismatches** | Ten profiles: 49,953 matching values and one matching source error; 46 timeout observations separately recorded |

The full maintained suite ran before the larger ChEMBL queue. After correction,
the affected CIP and conformer tasks were compared once more using their
existing references. These checks are separate evidence sets, not additional
ChEMBL input breadth. Fixed-input BIO and special-regression groups in the
standard suite retain their own recipe boundaries. Corrections are protected
by local behavior regressions, not new corpus fixtures.

Completed batch receipts and original reports remain under
`target/parity-20k-20261013/`. Source-localization and correction evidence is in
`target/parity-known-fixes-20261014/EVIDENCE.md`; the post-correction 5000 reports
are under `parity-tests/reports/run-7e9eo4/`. Original measurements are retained
alongside correction results. Generated inputs, references and reports are not
committed.

## Previous validation: COSMolKit 0.5.0 / RDKit 2026.03.6

<details>
<summary>Expand the previous ChEMBL 37 results and comparison boundaries</summary>

### ChEMBL 37 parity results

Completed validation comprises nine batches covering **18,000 ChEMBL 37
source records** and **20,505,775 parameterized comparisons**, with **zero
mismatches**. Million-scale validation remains pending.

These results and comparison rules apply to **COSMolKit 0.5.0 against RDKit
2026.03.6**, not the current development version targeting RDKit 2026.09.1.

### Reference and execution scope

- Chemistry reference: RDKit **2026.03.6**, source revision
  `0e0d85f4ca34aeae15dfc0f7cf5503bdb0a8e985`, CPython **3.13.12**.
- Corpus: ChEMBL 37, 2,897,819 available source records; source archive SHA-256
  `ea6181ce8dc7af41974e35b92e1febb0c9dcbe2c62f7ccc4a5d983ac19f696e7`.
- Selection: nine consecutive 2,000-record batches, retaining source order and
  ChEMBL IDs. The completed boundary is source offsets `0..18000`.
- Execution: the then-current `parity-tests_fixed`
  preparation and Cargo comparison, release optimization and
  `cosmolkit/op-contracts-strict`. Every completed batch's Cargo command exited
  successfully; together they produced 1,046 task reports.
- Only completed comparisons are included. The additional 2,000-record batch
  with unfinished reference preparation is excluded from the results below.

| Result | Count |
|---|---:|
| Completed source records | 18,000 |
| Actual comparisons | 20,505,775 |
| Mismatches | 0 |
| Reference timeout observations, not compared | 216 |
| Upstream crash observations, not compared | 27 |
| Upstream undefined observations, not compared | 0 |

A comparison is one registered input/profile observation, potentially
containing a vector, molecular state or provenance collection. Parameter
profiles multiply comparisons, not the number of distinct input molecules.
Matching parser rejections and defined error outcomes are included where the
task contract compares them; they are not successful molecule calculations.

### Completed functional coverage

All counts below come from the nine completed batch receipts. Each row groups
related registered tasks; the comparison column sums to the overall total.

| Surface | Actual comparisons | Compared behavior |
|---|---:|---|
| SMILES writing | 13,824,000 | Complete 768-profile writer matrix |
| SMILES parsing | 72,000 | Registered parser profiles and rejection outcomes |
| Molecular descriptors | 3,528,000 | 78 tasks, including weights/formula, Chi, Kappa/Phi, QED, TPSA, Crippen, MQN, ASA/contributions, VSA bins and molecular counts |
| Substructure matching | 324,000 | Registered matching parameter profiles |
| Morgan fingerprints | 1,152,000 | Four output forms, 16 ordinary profiles per form, including provenance |
| MACCS fingerprints | 18,000 | Registered raw/public fingerprint observations |
| Avalon fingerprints | 238,000 | Ordinary profiles, then the full 23-profile matrix |
| RDK topological fingerprints | 148,000 | Ordinary profiles, then the full 14-profile matrix |
| Layered fingerprints | 119,973 | Ordinary profiles, then the full 18-profile matrix |
| Pattern fingerprints | 72,000 | Ordinary profiles, then the full ten-profile matrix |
| AtomPair fingerprints | 80,000 | Ten profiles across four output forms, with provenance |
| Topological Torsion fingerprints | 72,000 | Nine profiles across four output forms, with provenance |
| Fuzzy fingerprint operations | 36,000 | Registered AND/OR pair observations |
| Tautomer enumeration and canonicalization | 264,000 | Ordered results, state and registered catalog/option profiles |
| Persistent MMFF/UFF force fields | 36,000 | Initial/final energy, gradients and positions with two-step minimization |
| MMFF single/conformer optimization | 71,856 | Prepared optimizer observations and defined outcomes |
| UFF single/conformer optimization | 35,928 | Prepared optimizer observations and defined outcomes |
| MMFF/UFF parameter availability | 72,000 | Registered parameter-availability branches |
| Molecular alignment and RMSD | 18,000 | Six rotated call families, including transforms, maps and source preservation |
| Direct 2D coordinates | 18,000 | Topology and coordinate outputs |
| SVG drawing | 18,000 | Final artifacts under the registered SVG comparison rules |
| Fragments and largest fragment | 36,000 | Registered fragment extraction and selection |
| Murcko decomposition/scaffold and net scaffold | 54,000 | Registered scaffold outputs |
| Sanitize, kekulize, hydrogen transforms and distance matrices | 198,000 | Registered graph/state and matrix observations |
| Batch SMILES | 18 | Two whole-batch observations per completed batch |

Historical parameter matrices were added at batch boundaries without rerunning
earlier completed batches. Their completed input breadth is therefore:

| Parameter matrix | Completed ChEMBL source records |
|---|---:|
| SMILES writer, 768 profiles | 18,000 |
| Tautomers, eight enumeration and eight canonicalization profiles | 16,000 |
| Avalon, 23 profiles; RDK topological, 14 profiles | 10,000 each |
| Layered, 18 profiles; Pattern, ten profiles | 6,000 each |
| AtomPair, 40 observations; Topological Torsion, 36 observations | 2,000 each |

The additional 15 historical dense-Morgan recipes passed a separate focused
standard prepare/Cargo comparison: **186/186 observations matched**, with no
timeouts or upstream crashes. This includes 96 ordinary requests, 90 historical
requests, 31 paired parser rejections and 75 actual historical fingerprint
vectors. Source-defined atom/bond eligibility is retained. This focused check
does not add ChEMBL input breadth; the corresponding large-batch comparison
is not included in the completed results.

### Comparison rules and reference exceptions

Fingerprint vectors, integer/string/state outputs and provenance use exact
comparison. Descriptor binary64 values, persistent MMFF/UFF numerical outputs
and ordinary UFF optimization compare their IEEE-754 bits. Ordinary MMFF
optimization compares energies and coordinates with absolute `1e-6` tolerance;
statuses and conformer IDs are exact. Ordinary MolAlign results and direct 2D
coordinates retain their declared absolute `1e-8` tolerance. Error and parser
outcomes follow each registered task's structural comparison contract; exception
message spelling alone is not a chemical mismatch. SVG uses the registered
artifact normalization rules. The run did not weaken these comparisons.

The 27 upstream crashes occurred in Layered fingerprint reference calls. The
216 reference timeouts were recorded as 72 MMFF single, 72 MMFF conformer,
36 UFF single and 36 UFF conformer observations. These are observation counts,
not necessarily distinct molecules; none is counted as a match. One earlier
memory-limit interruption was recovered using valid saved references and is
recorded separately from chemical mismatches and upstream crashes.

### Evidence and remaining coverage

Inputs, reference manifests, exact parameters, comparison reports and immutable
per-batch source identities are retained locally. The run's consolidated record
is `target/overnight-parity-20261010/morning-report.json`; batch receipts are
`completed.json` in the same directory. Native outcomes remain under
`parity-tests_fixed/expected/`, and comparison/exception reports under
`parity-tests_fixed/reports/`. Generated data are not committed.

The test-runner library validation recorded for this run passed 77 tests. Focused
adapter checks and replayed existing-reference comparisons are additional
runner evidence, not additional distinct ChEMBL records.

Remaining work includes million-scale breadth, uniform historical profiles
across that breadth, the complete descriptor/MolAlign historical boundaries,
modern CIP and potential-stereo/stereoisomer matrices, InChI diagnostics and
options, distance-geometry/conformer matrices, full MOL/SDF/binary-state
coverage, unsanitized/in-place hydrogen/cache observations, and full
composition/batch/concurrent-read matrices. The functional coverage above
records what was executed; it is not a complete 0.3.0-equivalent audit.

</details>

## Historical evidence: 0.3.0

<details>
<summary>Expand the previous full-corpus results and comparison boundaries</summary>

The historical results below used RDKit **2026.03.1**. They retain their
original corpus and API boundaries and are separate from the 0.5.0 run above.

## Corpus Tiers

The three corpus tiers are separate evidence sets. Their names must not be used
interchangeably in tests or documentation.

| Tier | Corpus | Records | Primary use |
|---|---|---:|---|
| Project small | `testdata/smiles/corpus/smiles_small.smi` | 152 | Fast daily parity, focused source-port regressions, and high-detail state checks |
| Maintained strict | `testdata/smiles/corpus/smiles_5000.smi` | 5,000 | Committed exhaustive profile matrices for surfaces not yet audited over ChEMBL 37 |
| Large stress | ChEMBL 37 structure table | 2,897,819 source; 2,897,804 mutually parseable | Large-scale parity, operation-order, batch, and concurrent-read stress auditing |

The complete ChEMBL 37 evidence inventory contains the 22-phase chemistry,
composition, batch, and concurrency audit; the complete modern CIPLabeler,
ordinary MolAlign, RDKFingerprint/Avalon, AtomPair, Pattern, Topological
Torsion, Layered fingerprint, tautomer, and potential-stereo/stereoisomer
enumeration audits; and the 15-phase execution recorded below. The
chemistry/composition audit contains 2,816 tasks over 128 corpus shards, while
the nine complete independent audit phases add 1,152 tasks over the same
shards.
Phases with a configured atom-count boundary evaluate 2,854,362 eligible
records. Large subset phases select their stated number of records from the
same ChEMBL 37 source; they are not separate corpora.

The repository-owned preparation, phase definitions, resumable runner, result
identity, and acceptance procedure are documented in
[`dev/tools/chembl_parity/README.md`](dev/tools/chembl_parity/README.md). The
current `complete.json` profile preserves the completed phases, adds four
configured phases for topology operations, fragments, direct 2D coordinates,
and binary roundtrips, and includes modern CIPLabeler, ordinary MolAlign, the
five complete fingerprint audit phases, tautomer enumeration, and
stereoisomer enumeration. The extended profile contains 35 phases and 4,480
shard tasks. The consolidated validation covers all additions over their
complete configured ChEMBL boundaries. Modern CIPLabeler, ordinary MolAlign,
fragments, direct 2D coordinates, topology operations, binary roundtrips,
tautomer enumeration, and stereoisomer enumeration each have a complete
accepted result with zero mismatch. The summary below reports this final
evidence boundary.
Corpus shards and run outputs remain uncommitted; code, profiles, source and
reference pins, and pass/fail rules are version controlled.

## Tautomer Validation

The complete ChEMBL 37 tautomer matrix processes all 2,897,819 inputs across
128 shards and covers 211,539,707 exact output and molecular-state
observations. The accepted result has zero mismatch across parse behavior and
all eight configured enumeration profiles.

The complete-corpus result is closed by exhaustive post-correction verification
of every source row implicated during alignment and a deterministic
524,288-record regression selected across all 128 corpus shards. The
cross-shard regression executes the same complete branch matrix and verifies
that the source-level corrections preserve the previously matching corpus
surface. Focused fixtures independently lock every discovered source boundary.

## Complete ChEMBL 37 Validation Summary

The consolidated validation record combines the 22-phase audit, the complete
modern CIPLabeler, ordinary MolAlign, RDKFingerprint/Avalon, AtomPair, Pattern,
Topological Torsion, Layered fingerprint, tautomer, and stereoisomer-enumeration
audits, and the 15-phase execution. The 15-phase execution used COSMolKit
`0.2.12`, pinned RDKit `2026.03.1`, the same 2,897,819-record ChEMBL 37 source,
128 corpus shards, and 112 workers. All 15 phases and all 1,920 shard tasks
completed; every result artifact and checksum was independently revalidated.
Together with the nine independent complete audits, the consolidated evidence
table records 4,176,585,416 matching checks and zero blocking mismatch. Bounds
additionally traversed 2,757,910,995 matrix entries.

The first 15 rows summarize the consolidated chemistry and composition
execution. The final nine rows are independent complete modern CIPLabeler,
ordinary MolAlign, fingerprint, tautomer, and stereoisomer-enumeration audits
over the same ChEMBL 37 source.

| Phase or independent audit | Records evaluated | Matching checks | Blocking mismatches | Status |
|---|---:|---:|---:|---|
| Descriptors | 2,897,819 source; 2,897,804 mutually parseable | 987,273,647 | 0 | Pass |
| Morgan/MACCS | 2,897,819 | 89,831,720 | 0 | Pass |
| Explicit hydrogen | 2,854,362 | 11,417,448 | 0 | Pass |
| Topology operations | 2,854,376 | 45,669,848 | 0 | Pass |
| Connected fragments | 2,854,362 | 5,852,676 | 0 | Pass |
| Direct 2D coordinates | 2,854,362 | 14,271,810 | 0 | Pass |
| Binary roundtrip | 524,288 | 11,534,336 | 0 | Pass |
| Bounds | 2,854,362 | 2,854,362 | 0 | Pass |
| InChI parse branches | 2,854,362 | 11,417,448 | 0 | Pass |
| SMILES matrix | 2,854,362 | 2,192,150,016 | 0 | Pass |
| Final SVG | 2,854,362 | 2,854,362 | 0 | Pass |
| Molecular I/O | 524,288 | 12,559,920 | 0 | Pass |
| Scalar/batch composition | 1,027,660 | 20,561,265 | 0 | Pass |
| Fixed-seed conformer outcome | 524,288 | 524,288 | 0 | Pass |
| UFF/MMFF force fields | 524,288 | 3,139,761 | 0 | Pass |
| Modern CIPLabeler | 2,854,362 compared; 15 parser rejects; 43,442 over 80 atoms | 11,417,448 | 0 | Pass |
| Ordinary MolAlign and RMSD | 2,854,362 compared; 15 parser rejects; 43,442 over 80 atoms | 11,417,207 | 0 | Pass |
| RDKFingerprint/Avalon fingerprints | 2,897,819 source; 2,897,804 compared | 113,014,356 | 0 | Pass |
| AtomPair fingerprints | 2,897,819 source; 2,897,804 compared | 118,809,964 | 0 | Pass |
| Pattern fingerprints | 2,897,819 source; 2,897,804 compared | 28,978,040 | 0 | Pass |
| Topological Torsion fingerprints | 2,897,819 source; 2,897,804 compared | 127,503,376 | 0 | Pass |
| Layered fingerprints | 2,897,819 source; 2,897,804 compared | 52,160,472 | 0 | Pass |
| Tautomer enumeration | 2,897,819 source | 211,539,707 | 0 | Pass |
| Potential stereo and stereoisomer enumeration | 2,897,819 source; 2,897,804 accepted by both parsers | 89,831,939 | 0 | Pass |

The topology phase covers the `RemoveHs(sanitize=false)` value and in-place
branches, including RDKit's non-strict property-cache update before removal
and the surviving explicit-valence and implicit-H fields afterward. COSMolKit
migrates those cache rows through the declared topology mapping. The complete
128-shard result processed 2,854,376 records and produced 45,669,848 exact
matches with no mismatch or finding.

The binary phase requires validated valence, ring, aromaticity, and stereo
state to survive both public deserialization entrypoints while preserving
legacy and archive-v1.0 read compatibility. The complete configured boundary
produced 11,534,336 matches and no mismatch or finding across 524,288 records.

The immutable identity and detailed acceptance record are in
[`dev/gap_reports/chembl37_complete_validation_summary_20260820.md`](dev/gap_reports/chembl37_complete_validation_summary_20260820.md).
The run validates the recorded extension and start-of-run repository snapshot.
Later uncommitted working-tree changes require a rebuilt extension and a new
identity before they can inherit this evidence.

## Numeric Coverage Summary

The comparison count is the number of output values, states, matrices, or
matrix entries actually compared, not merely the number of input molecules.

| Surface | Corpus | Records evaluated | Profiles or output branches | Comparisons |
|---|---|---:|---:|---:|
| Core graph and default scalar APIs | ChEMBL 37 | 2,897,804 | 9 | 26,080,236 |
| Molecular descriptors | ChEMBL 37 | 2,897,804 | 29 established outputs + complete high-feasibility scalar/vector/contribution/cache matrix | 987,273,647 |
| Morgan and MACCS fingerprints | ChEMBL 37 | 2,897,804 | 31 configured | 89,831,720 completed |
| Topological Torsion fingerprints | ChEMBL 37 | 2,897,804 | 36 vectors + 8 provenance outputs | 127,503,376 |
| Layered fingerprints | ChEMBL 37 | 2,897,804 | 18 complete bit/count result profiles | 52,160,472 |
| Explicit-hydrogen graph/state | ChEMBL 37 | 2,854,362 | 4 | 11,417,448 |
| Distance-geometry bounds | ChEMBL 37 | 2,854,362 matrices | every matrix entry | 2,757,910,995 entries |
| InChI parse branches | ChEMBL 37 | 2,854,362 | 4 | 11,417,448 |
| InChI writer options and keys | ChEMBL 37 | 2,854,362 | 10 options x 2 outputs | 57,087,240 |
| SMILES writer matrix | ChEMBL 37 | 2,854,362 | 768 | 2,192,150,016 per run; 3 runs |
| Final SVG | ChEMBL 37 | 2,854,362 | 1 | 2,854,362 |
| MolBlock/SDF I/O | ChEMBL 37 subset | 524,288 | 16 writer + 8 state/coordinate | 12,559,920 completed |
| Scalar operation-order composition | ChEMBL 37 subset | 516,650 | 4 orders x 12 surfaces | 24,799,200 |
| Scalar-versus-batch execution | ChEMBL 37 subset | 1,027,660 | 2 job counts x 6 outputs | 12,331,920 |
| Shared-object concurrent reads | ChEMBL 37 subset | 64,524 | 16 values + 16 graph states | 2,064,768 per run; 2 runs |
| Fixed-seed conformer outcome | ChEMBL 37 subset | 524,288 | coordinate or matching failure status | 524,288 |
| UFF/MMFF parameter and optimizer paths | ChEMBL 37 subset | 516,096 | 4 parameter + 2 optimizer | 3,090,706 completed |
| RDKFingerprint/topological | ChEMBL 37 | 2,897,804 | 14 vectors + 2 provenance outputs | 46,364,864 |
| Avalon explicit-bit | ChEMBL 37 | 2,897,804 | 23 | 66,649,492 |
| AtomPair | ChEMBL 37 | 2,897,804 | 40 vectors + 1 provenance output | 118,809,964 |
| Pattern | ChEMBL 37 | 2,897,804 | 10 complete vectors | 28,978,040 |
| Modern CIPLabeler | ChEMBL 37 | 2,854,362 | full, selected atom, selected bond, empty selection | 11,417,448 |
| Ordinary MolAlign and RMSD | ChEMBL 37 | 2,854,362 | 6 rotated call families; transforms, maps, RMSD, coordinates, order, immutability | 11,417,207 |
| Tautomer enumeration | ChEMBL 37 | 2,897,819 | parse + 8 complete enumeration profiles | 211,539,707 |
| Potential stereo and stereoisomer enumeration | ChEMBL 37 | 2,897,819 source; 2,897,804 compared | 4 potential-stereo + 4 bounded enumeration profiles | 89,831,939 |

Except for the separately identified modern CIPLabeler, ordinary MolAlign,
RDKFingerprint, Avalon, AtomPair, Pattern, Topological Torsion, Layered,
tautomer, and stereoisomer-enumeration rows,
the ChEMBL 37 figures above are the complete chemistry/composition audit
counters. The acceptance model combines that full-corpus breadth with exact
source-regression matrices for every covered branch: 10 descriptor rows,
3,092 fingerprint rows, one explicit-hydrogen row, 180 bounds matrices with
308,734 entries, 50 InChI rows with 200 parse-branch comparisons, 108 SMILES
rows with 82,944 writer comparisons, 441 SVG rows, 160 I/O rows with 2,560
writer comparisons, seven conformer rows, and 262 force-field rows. Every
listed regression matrix has exact pinned-RDKit agreement; no percentage
threshold or unexplained mismatch is accepted. RDKFingerprint, Avalon,
AtomPair, Pattern, Topological Torsion, Layered, and tautomer enumeration
additionally have complete current full-corpus evidence.

### Potential Stereo And Stereoisomer Enumeration

The project-native value API covers typed, source-ordered potential-stereo
records and the pinned RDKit Python ``EnumerateStereoisomers`` behavior. Four
potential-stereo profiles span both cleanup and possible-center switches. Four
enumeration profiles cover default bounded output, assigned-center inclusion,
non-unique output, and deterministic seeded random generation. Comparisons
include parse acceptance, ordered records, affected molecule state,
theoretical counts, bounded flags, ordered full output state, source
preservation, and parser-rejection outcomes.

The complete 128-shard ChEMBL 37 execution processed all 2,897,819 records and
fully compared 2,897,804 rows accepted by both parsers. The remaining 15 rows
were rejected by both parsers at parse entry. It produced 89,831,939 exact
observations with no blocking, informational, or other mismatch. Focused
source fixtures, the 152-record project corpus, and the maintained 5,000-record
corpus additionally cover option defaults, enhanced groups, random draw
sequences, arbitrary-width counts, uniqueness, iterator errors, embedding
outcomes, and full-state regressions. This boundary is the Python enumerator;
it does not claim the separate newer C++ enumerator or atropisomer enumeration.

The 3,480 archived ChEMBL reference InChIs and 3,480 archived reference
InChIKeys that differ from the current result are corpus-version differences:
COSMolKit agrees with pinned RDKit on those rows. They are not COSMolKit parity
failures.

## Comparison Details

### Molecular Graph And Hydrogen State

The core graph comparison covers atom identity, charge, degree, hydrogens,
radicals, hybridization, aromaticity, ring membership, bond type, bond stereo,
conjugation, CIP rank/code, atom count, and bond count. The ChEMBL core phase
also compares default SMILES and the four default InChI scalar surfaces. The
explicit-hydrogen phase separately compares atom state, bond state, atom count,
and bond count after the transform.

### SMILES Writing

The 768-profile matrix is the Cartesian expansion of isomeric, kekule,
canonical, clean-stereo, explicit-bond, explicit-H, dative-bond, and
atom-map-ignore switches, with rooted output for no root, first atom, and last
atom. Canonicalized SMILES are compared directly in canonical branches. The
matrix was run three times to distinguish deterministic source differences from
execution instability. Upstream writer options not enumerated by this matrix
are separate capabilities outside the current parity claim; every enumerated
branch remains subject to the zero-mismatch rule.

### InChI

The four public scalar APIs are `Molecule.to_inchi()`,
`Molecule.to_inchi_key()`, top-level `inchi_to_key()`, and
`Molecule.from_inchi()`. Comparisons cover official return fields, diagnostics,
graph atom/bond/stereo/isotope/H/charge/radical state, ten writer option
profiles, cleanup, malformed input, source preservation, and concurrent calls.

Exact parity applies where official InChI v1.07.5 and RDKit `2026.03.1` define
behavior. The official `NormalizeAndCompare` initial-buffer allocation-failure
path is undefined C behavior; Rust returns a deterministic structured
`allocation_failed` error. MolBlock, SDF/V3000, IXA, AuxInfo, INCHIGEN, version
query, and extended-polymer APIs are separate upstream surfaces outside the
four-scalar-API claim.

### Molecular Descriptors

The established 29 ChEMBL branches cover molecular weight, exact molecular
weight, formula variants, H-bond donor/acceptor counts, fraction Csp3,
Crippen LogP/MR option branches, TPSA with and without S/P contributions and
force recomputation, aromatic-ring count, four rotatable-bond modes, and QED.

The complete high-feasibility matrix adds connectivity Chi orders 0 through
6, fixed valence and principal-quantum Chi projections, Hall-Kier alpha and
atom contributions, Kappa/Phi, direct Lipinski and ring/stereo counts, every
entry of the 42-component MQN vector, Labute ASA and atom/H contributions,
and every default or custom SlogP/SMR VSA bin. It separately compares cold,
warm, forced, and mixed-order cache results. These additions contribute
903,237,331 exact ChEMBL observations beyond the 84,036,316 established
descriptor comparisons.

| Descriptor evidence | Input records | Mutually parseable records | Exact observations | Mismatches |
|---|---:|---:|---:|---:|
| Focused behavioral fixtures | 69 | 69 | 17,265 | 0 |
| Project small corpus | 152 | 140 | 40,224 | 0 |
| Maintained 5,000-row corpus | 5,000 | 5,000 | 1,588,749 | 0 |
| ChEMBL 37 | 2,897,819 | 2,897,804 | 903,237,331 | 0 |

Every finite floating-point leaf uses exact IEEE-754 binary64-bit comparison;
integer and string outputs use exact equality. Contribution arrays preserve
atom order, vectors compare every entry, and paired parser rejections do not
enter the descriptor-observation count. QED is pinned to
`RDKit 2026.03.1 + CPython 3.13.12` because RDKit delegates reductions to
Python `sum()` and older CPython versions use a different floating-point
reduction algorithm. The implementation boundary and reproduction commands
are recorded in the [descriptor validation report](dev/gap_reports/rdkit_high_feasibility_descriptors_full_port_validation.md).

### Tautomer Enumeration

The tautomer audit compares parser acceptance and eight source-defined
enumeration profiles: default and V1 catalogs, one-tautomer and one-transform
limits, retained tetrahedral stereo, retained double-bond stereo, retained
isotopic hydrogens, and disabled stereo reassignment. Every successful branch
compares the ordered tautomer SMILES, complete molecule states, enumeration
status, modified atom and bond sets, per-tautomer scores, canonical tautomer
SMILES, and canonical molecule state; error outcomes are compared exactly.

The final source alignment preserves RDKit's separate legacy and modern
tetrahedral-center predicates, ring-stereo cleanup order, canonical-writer
traversal gate, and connected-component computed-property lifecycle. These are
shared chemistry and writer paths rather than corpus-local tautomer branches.
The focused fixture set locks each boundary, and the deterministic cross-shard
regression covers all eight profiles over 524,288 ChEMBL 37 records.

### Modern CIPLabeler

The modern CIPLabeler audit uses the reproducible ChEMBL runner and the same
128 deterministic corpus shards as the other complete audits. Each of the
2,854,362 eligible records is compared after full assignment, selected-atom
assignment, selected-bond assignment, and empty-selection dispatch. All
11,417,448 comparisons match pinned RDKit exactly; all 128 shards completed
with zero mismatch or retained finding. Fifteen source records are rejected by
both parsers, and 43,442 accepted records exceed the phase's configured
80-atom boundary.

Each comparison includes success or exact error state, molecule
`_CIPComputed`, atom and bond `_CIPCode` and `_CIPNeighborOrder`, atom
`_CIPRank`, computed-property membership, chiral tags, bond stereo, and stereo
atoms. The pinned source boundary covers tetrahedral, double-bond, and
atropisomeric configurations constructed by `findConfigs`; descriptor families
outside that dispatcher are separate capabilities. The immutable result and
reproduction command are recorded in the [modern CIPLabeler validation
report](dev/gap_reports/rdkit_ciplabeler_full_port_validation.md).

### Molecular Alignment And RMSD

The ordinary MolAlign audit rotates alignment-transform, best-alignment,
coordinate-frame RMSD, explicit value-style alignment, all-conformer best-RMS,
and conformer-alignment calls across stable ChEMBL row numbers. It compares
RMSD values and complete 4x4 transforms to `1e-8`, atom maps and conformer-pair
order exactly, changed coordinates to `1e-8`, and complete source/reference
coordinate snapshots for immutability.

All 2,854,362 eligible molecules produced 11,417,207 matching checks with zero
mismatch across 128/128 shards. Fifteen source rows were rejected by both
parsers and 43,442 accepted rows exceeded the configured 80-atom boundary.
Focused, 152-row, and explicit exhaustive 5,000-row gates additionally cover
weights, reflection, iteration limits, sparse conformer IDs, map validation,
terminal-group symmetry, structured errors, and exact mutation separation.
The [MolAlign validation report](dev/gap_reports/rdkit_molalign_full_port_validation.md)
records the source boundary and reproduction commands. O3A and its MMFF/Crippen
scoring path are separate capabilities and are not included in ordinary
MolAlign.

### Morgan And MACCS Fingerprints

Morgan compares 15 parameter profiles as both the main vector and
`AdditionalOutput`; MACCS adds one output branch. The large-run counter omits
102 rows from each of two custom-bond-invariant output counters because of an
audit-harness accounting gap, not an accepted mismatch. Exact source
regressions cover both affected branches. The selected boundary is documented
in the [Morgan/MACCS validation report](dev/gap_reports/rdkit_morgan_maccs_full_port_validation.md).

### RDKFingerprint And Avalon Fingerprints

RDKFingerprint comparisons include vector width, the complete on-bit set,
`atomBits`, and `bitInfo`. Avalon compares complete explicit-bit vectors across
23 profiles. Their combined 128-shard audit completed with zero mismatch or
failed task. Family-specific boundaries and reproducibility identities are in
the [RDKFingerprint](dev/gap_reports/rdkit_topological_fingerprint_full_port_validation.md)
and [Avalon](dev/gap_reports/avalon_fingerprint_full_port_validation.md) validation
reports.

### AtomPair Fingerprints

AtomPair compares sparse-count, folded-count, sparse-bit, and explicit-bit
outputs. Its ten profiles cover chirality, count simulation and bounds,
distance bounds, atom selection and exclusion, custom atom invariants, and
`numBitsPerFeature`; the provenance profile compares atom-to-bit, atom-count,
and bit-info mappings. The source-aligned implementation also covers 2D and 3D
distance-cache reuse, conformer keys, invalidation, and concurrent cold/warm
reads. The full audit has zero mismatch, invalid profile, or retained finding;
details are in the [AtomPair validation report](dev/gap_reports/rdkit_atom_pair_fingerprint_full_port_validation.md).

### Pattern Fingerprints

The Pattern boundary is pinned to RDKit `2026.03.1` revision
`351f8f378f8ad6bbd517980c38896e66bf907af8` and its 13 built-in Pattern SMARTS
queries in source order. The ChEMBL 37 execution completed all 128/128
deterministic shards against corpus SHA-256
`ea6181ce8dc7af41974e35b92e1febb0c9dcbe2c62f7ccc4a5d983ac19f696e7`.

Pattern compares exact width and the complete ordered on-bit set for default,
tautomeric, widths 1 and 127, a collision-heavy seven-bit tautomeric profile,
zeroed and pre-populated atom-count inputs, and all-zero, all-one, and sparse
set-only masks. The latter two source arguments are size-validated but
observably inert in the pinned ordinary overload; the port preserves that
behavior instead of implementing the stale header description.

All 2,897,804 mutually parseable ChEMBL 37 records match across all ten
profiles, totaling 28,978,040 exact comparisons with zero mismatch. Focused,
small, and maintained 5,000-row gates use 11 branches because they additionally
execute the deterministic wrong-width `setOnlyBits` precondition error. They
also cover query-bearing graphs, validation errors, source screening
regressions, exact match-derived intermediates, scalar/batch composition, and
concurrent compile-once cache reuse. RDKit's `MolBundle` intersection overload
is outside COSMolKit's ordinary molecule and ordered batch model. Behavioral
parity is complete for the stated boundary, while the measured
`1.508-3.104x` warm runtime ratio is reported without claiming performance
parity. Details are in the [Pattern validation
report](dev/gap_reports/rdkit_pattern_fingerprint_full_port_validation.md).

### Topological Torsion Fingerprints

Topological Torsion is distinct from the RDKFingerprint path/subgraph family.
Its modern and legacy boundary covers live options, selections, chirality,
path behavior, JSON, collisions, structured errors, and all four vector forms.
Requested provenance compares atom counts, atom-to-bits, bit paths,
atoms-per-bit, and the source-defined empty bit-info map.

Focused state tests additionally cover dynamic unfolded widths, zero live
`numBitsPerFeature`, empty live count bounds, zero live folded size, query and
dative graphs, hypervalent sulfur, explicit hydrogens, and the legacy boron
endpoint correction. Fingerprinting reads the source explicit-valence cache;
`with_hydrogens()` commits the source-defined transition needed by both modern
and legacy calls. The complete audit has zero mismatch, invalid profile, or
retained finding; details are in the [Topological Torsion validation
report](dev/gap_reports/rdkit_topological_torsion_fingerprint_full_port_validation.md).

### Layered Fingerprints

Layered fingerprint comparisons cover all six active source layers, active,
default, substructure, zero, and inactive-high layer masks, two width/path
configurations, branched and rooted-linear paths, rooted branched selections,
bit masks, and seeded atom counts. Every successful profile compares the exact
fingerprint width, complete ordered on-bit set, and complete optional atom-count
vector. The complete 128-shard ChEMBL 37 audit evaluated all 18 profiles for
every mutually parseable row, totaling 52,160,472 exact comparisons with zero
mismatch, invalid profile, or retained finding; details and immutable run
identity are in the [Layered validation
report](dev/gap_reports/rdkit_layered_fingerprint_full_port_validation.md).

Pinned RDKit's unrooted-linear call passes atom paths into code that consumes
bond paths and can terminate the process. COSMolKit does not reproduce that
crash: it returns the bond-path result documented by the source API and used by
the valid rooted-linear branch. This process-safety difference is excluded
from, rather than counted as part of, the zero-mismatch total.

### Fingerprint Acceptance Rules

Fingerprint families remain separate parity boundaries even when they share
vector, hashing, distance, or generator infrastructure. Complete fingerprint
audits process all 2,897,819 ChEMBL 37 source records: 2,897,804 are mutually
parseable and 15 are rejected by both libraries.

The maintained 5,000-row gates consume all 14 RDKFingerprint, 23 Avalon, ten
AtomPair, 11 Pattern branches, nine Topological Torsion, and 18 Layered
profiles. Each active profile defines the vector forms and provenance it
claims. Comparisons use exact equality, deterministic ordered collection, and
validated manifests; there is no sampling, tolerance, skipped mismatch,
approximate vector, similarity threshold, fallback fingerprint, or runtime
RDKit dependency. Run identities and aggregate checksums remain in the linked
validation reports instead of being repeated in this scope summary.

### Bounds, Conformers, And Force Fields

Bounds comparisons cover matrix shape and every lower/upper entry, including
topology, ring, macrocycle, amide, VDW, and triangle-smoothing effects. Entry
tolerance is `1e-8`.

Conformer comparisons fix seeds and execution conditions. They cover
DG/KDG/ETDG/ETKDG presets, single-thread execution, random-coordinate
embedding, coordMap, custom pairwise terms, bounds-matrix handling,
timeout/max-iteration controls, multi-conformer generation, and sequential seed
policy. The ChEMBL subset produced 421,217 coordinate-bearing outcomes and
103,071 matching failure statuses. Successful seeded coordinates use a `1e-6`
tolerance.

The ChEMBL force-field subset compares UFF/MMFF parameter availability with and
without explicit hydrogens plus optimizer outcomes. The project small tier also
contains 150 curated cases that compare seeded initial coordinates, initial
energy, every gradient component, single- and multi-conformer optimization,
result codes, final energies, and final coordinates. Energy, gradient, and
coordinate tolerance is `1e-6`. Parameter availability over a large corpus is
not presented as equivalent to the curated full-gradient boundary.

### Molecular I/O And Depiction

MolBlock/SDF writer comparisons cover 16 V2000/V3000 branches over 2D and 3D
sources, including stereochemistry, kekulization, SGroups, RGroups, aliases,
value lines, and aromatic-bond bookkeeping. Reader comparisons cover topology,
atom fields, bond fields, coordinates, marker- and coordinate-derived
chirality, delayed sanitize/remove-H paths, and parsed-SMILES checks. This is
broad CTAB/SDF coverage, not a claim for every extension.

MOL2 reading is source-ported for the exposed Tripos parser parameters,
topology, atom/bond state, 3D coordinates, chirality, CORINA behavior,
cleanup-substructure behavior, and SMILES output. XYZ reading compares atom
identity, one 3D conformer, coordinates, and the absence of inferred bonds.

Prepared 2D drawing state covers depiction kekulization, chiral-H insertion,
wedge bonds, atom order, 2D coordinates, bond order, aromatic flags, and bond
directions. Final SVG text is compared after normalizing only tool
namespace/prefix metadata. Link-node extraction, StereoGroup masking,
atomRegions, and other separately named drawing extensions are outside this
SVG boundary. PNG is checked as local SVG rasterization and is not claimed to
have RDKit Cairo/Qt bit parity.

### Composition, Batch, And Concurrent Reads

Four scalar operation orders each compare 12 named output/state surfaces,
including graph state, explicit hydrogens, bounds, descriptors, Morgan, MACCS,
SMILES, InChI/InChIKey, and SVG. Batch comparisons run add/remove-H, bounds,
Morgan, Morgan `AdditionalOutput`, SMILES, SVG, raw-graph sanitize, both
kekulize flag branches, and direct 2D coordinates with one and eight jobs while
checking scalar equality and input order. Shared-object tests read the same
molecule concurrently and compare both returned values and graph state; they
do not define separate chemistry semantics.

### Configured Additional ChEMBL Coverage

The extended profile adds strict structural-state comparisons for raw parse
and full-state comparisons after sanitize, value and in-place kekulization
with both aromatic-flag branches, AddHs followed by RemoveHs with both sanitize
branches, connected fragments, and direct 2D coordinate generation. Graph
comparisons retain atom valence and hydrogen-cache fields so unsanitized
intermediate-state differences cannot be hidden by an equal final SMILES.
Value-style calls additionally require the source molecule to remain
byte-identical.

The raw topology entrance processes 2,854,376 records with at most 80 atoms.
It deliberately includes 14 records outside the ordinary 2,854,362-record
sanitized boundary so raw-success/sanitize-rejection parity is exercised.

The binary-roundtrip phase creates 2D coordinates and compares graph state,
SMILES, descriptors, hash, fingerprint, coordinates, conformer count, and
deterministic reserialization across both public deserialization entrypoints.
This is a COSMolKit serialization invariant exercised on ChEMBL structures,
not a claim that COSMolKit and RDKit share a binary representation.

The consolidated 15-phase validation covers these four phases across every
configured shard with no blocking mismatch. The topology phase passes every
value-style and in-place branch, including unsanitized valence/H state. The
binary phase passes graph, SMILES, Kekule SMILES, hash, Morgan, descriptor,
coordinate, conformer-count, and deterministic-byte observations through both
deserialization entrypoints. Focused source regressions independently lock the
corresponding state-transition behavior.

### Focused And Scope-Bound Surfaces

The project regression tier contains 77 fixed full-state oracle records for
`with_chiral_tags_from_structure`. They compare status/errors, selected
conformer, every atom chiral tag and permutation, `_NonExplicit3DChirality`,
`_StereochemDone`, complete atom/bond state, coordinates/dimensionality,
replacement behavior, environment-controlled non-tetrahedral assignment, and
unrelated properties. This covers `assignChiralTypesFrom3D`, not
`assignStereochemistryFrom3D`, 3D double-bond E/Z assignment, CIP orchestration,
or distinct-substituent validation.

Tetrahedral stereo geometry separately checks ordered-ligand orientation and
signed volume against seeded RDKit ETKDG geometries. SMARTS parsing, writing,
and direct query matching use one canonical ordinary-molecule query graph;
the pinned RDKit 2026.03.1 source corpus contains 91 distinct accepted parser
inputs, nine rejected inputs, and 62 matcher rows. It compares parser
acceptance, atom/bond counts, exact molecule SMARTS output, and exact ordered
atom mappings. Separate strict tests exercise atom/bond/fragment/CXSMARTS
composition, parser parameters, matcher parameters and callbacks, Python
conversion, and existing consumers; those tests are not conflated with the
162-row golden corpus. Reaction and database/container SMARTS remain outside
the claim. The Bison-specific `debug_parse=true` diagnostic stream is an
explicitly unsupported diagnostic mode, not a chemistry result inside this
parity boundary.

## Interpretation Rules

- Strict parity means 100% agreement on every covered case. A mismatch is work
  to investigate, not an accepted error rate.
- RNG-sensitive parity fixes seed and execution conditions before comparing
  coordinates or optimizer state.
- Final-output parity compares the serialized user-visible artifact directly.
- Outside scope means a separately identified upstream capability is not part
  of the claim. It never means that a failing input inside a covered branch is
  accepted or relabeled unsupported.
- Golden baselines are condition-bound artifacts. Changing corpus, RDKit
  version, branch matrix, normalization, RNG/coordinate settings, or compared
  fields requires a new profile directory or file name.

</details>
