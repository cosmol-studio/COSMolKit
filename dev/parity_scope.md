# Feature Parity Scope

This page records the exact corpus and comparison boundary behind COSMolKit's
RDKit parity claims. RDKit `2026.03.1` is the current chemistry reference. A
feature is parity-covered only at the boundary listed here. Separate upstream
APIs outside that boundary are not claimed as implemented.
This scope distinction never permits a mismatching row inside a covered
boundary to be reclassified as unsupported.

## Corpus Tiers

The three corpus tiers are separate evidence sets. Their names must not be used
interchangeably in tests or documentation.

| Tier | Corpus | Records | Primary use |
|---|---|---:|---|
| Project small | `testdata/smiles/corpus/smiles_small.smi` | 152 | Fast daily parity, focused source-port regressions, and high-detail state checks |
| Maintained strict | `testdata/smiles/corpus/smiles_5000.smi` | 5,000 | Committed exhaustive profile matrices for surfaces not yet audited over ChEMBL 37 |
| Large stress | ChEMBL 37 structure table | 2,897,819 source; 2,897,804 mutually parseable | Large-scale parity, operation-order, batch, and concurrent-read stress auditing |

The complete ChEMBL 37 evidence inventory contains the 22-phase chemistry,
composition, batch, and concurrency audit; the later complete
RDKFingerprint/Avalon audit; and the final 15-phase execution recorded below.
The original audit contains 2,816 tasks over 128 corpus shards, while the
fingerprint audit adds 128 tasks over the same shards.
Phases with a configured atom-count boundary evaluate 2,854,362 eligible
records. Large subset phases select their stated number of records from the
same ChEMBL 37 source; they are not separate corpora.

The repository-owned preparation, phase definitions, resumable runner, result
identity, and acceptance procedure are documented in
[`dev/tools/chembl_parity/README.md`](tools/chembl_parity/README.md). The
current `complete.json` profile preserves the completed phases and adds four
configured phases for topology operations, fragments, direct 2D coordinates,
and binary roundtrips. The extended profile contains 27 phases and 3,456 shard
tasks. The final 15-phase execution ran all four additions over their complete
configured ChEMBL boundaries. Fragments and direct 2D
coordinates were clean; topology operations and binary roundtrips exposed
blocking systemic differences. The consolidated summary records both the clean
and failed boundaries while keeping failed phases outside accepted parity
evidence.
Corpus shards and run outputs remain uncommitted; code, profiles, source and
reference pins, and pass/fail rules are version controlled.

## Complete ChEMBL 37 Validation Summary

The consolidated validation record combines the original 22-phase audit, the
complete RDKFingerprint/Avalon audit, and the final 15-phase execution. The
final execution used COSMolKit `0.2.12`, pinned RDKit `2026.03.1`, the same
2,897,819-record ChEMBL 37 source, 128 corpus shards, and 112 workers. All 15
phases and all 1,920 shard tasks completed; every result artifact and checksum
was independently revalidated. It performed 31,755,360 phase-record executions
and recorded 2,500,870,740 passing checks plus 7,804,836 blocking mismatch
checks. Bounds additionally traversed 2,757,910,995 matrix entries.

| Evidence boundary | Result |
|---|---|
| Original 22-phase chemistry/composition audit | Complete; source-fixed retained-branch replays are recorded without claiming a complete post-fix corpus rerun |
| Complete RDKFingerprint/Avalon audit | 113,014,356 exact comparisons, zero mismatch |
| Final 15-phase execution | 13 phase aggregates pass; topology operations and binary roundtrip fail |

| Phase group | Records evaluated | Matching checks | Blocking mismatches | Status |
|---|---:|---:|---:|---|
| Descriptors | 2,897,819 | 84,036,316 | 0 | Pass |
| Morgan/MACCS | 2,897,819 | 89,831,720 | 0 | Pass |
| Explicit hydrogen | 2,854,362 | 11,417,448 | 0 | Pass |
| Topology operations | 2,854,376 | 39,962,026 | 5,707,822 | **Fail** |
| Connected fragments | 2,854,362 | 5,852,676 | 0 | Pass |
| Direct 2D coordinates | 2,854,362 | 14,271,810 | 0 | Pass |
| Binary roundtrip | 524,288 | 9,437,322 | 2,097,014 | **Fail** |
| Bounds | 2,854,362 | 2,854,362 | 0 | Pass |
| InChI parse branches | 2,854,362 | 11,417,448 | 0 | Pass |
| SMILES matrix | 2,854,362 | 2,192,150,016 | 0 | Pass |
| Final SVG | 2,854,362 | 2,854,362 | 0 | Pass |
| Molecular I/O | 524,288 | 12,559,920 | 0 | Pass |
| Scalar/batch composition | 1,027,660 | 20,561,265 | 0 | Pass |
| Fixed-seed conformer outcome | 524,288 | 524,288 | 0 | Pass |
| UFF/MMFF force fields | 524,288 | 3,139,761 | 0 | Pass |

The topology differences are confined to the `RemoveHs(sanitize=false)`
value and in-place branches: each has 2,853,911 mismatches and 451 matches. All
retained examples differ only in explicit valence, implicit-H, and total-H
intermediate-state fields. The binary differences affect hash behavior on all
524,288 rows and Morgan on 524,219 rows through each deserialization entrypoint;
graph state, conformer count, and deterministic reserialization still match.
Both patterns were reproduced after completion and remain implementation
defects, not harness, shard, or aggregation failures.

The immutable identity and detailed acceptance record are in
[`dev/gap_reports/chembl37_complete_validation_summary_20260820.md`](gap_reports/chembl37_complete_validation_summary_20260820.md).
The run validates the recorded extension and start-of-run repository snapshot.
Later uncommitted working-tree changes require a rebuilt extension and a new
identity before they can inherit this evidence.

## Numeric Coverage Summary

The comparison count is the number of output values, states, matrices, or
matrix entries actually compared, not merely the number of input molecules.

| Surface | Corpus | Records evaluated | Profiles or output branches | Comparisons |
|---|---|---:|---:|---:|
| Core graph and default scalar APIs | ChEMBL 37 | 2,897,804 | 9 | 26,080,236 |
| Molecular descriptors | ChEMBL 37 | 2,897,804 | 29 | 84,036,316 |
| Morgan and MACCS fingerprints | ChEMBL 37 | 2,897,804 | 31 configured | 89,831,720 completed |
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

Except for the separately identified RDKFingerprint and Avalon rows, the
ChEMBL 37 figures above are the original complete audit's counters. That audit
deliberately retained every observed difference instead of accepting a
percentage threshold. Source-aligned fixes were then replayed over the retained
records and their full affected branch matrices: 10 descriptor records, 3,092
fingerprint records, one explicit-hydrogen record, 180 bounds matrices with
308,734 entries, 50 InChI records with 200 parse-branch comparisons, 108 SMILES
records with 82,944 writer comparisons, 441 SVG records, 160 I/O records with
2,560 writer comparisons, seven conformer records, and 262 force-field records.
Those retained-set replays have no unexplained mismatch. They are not described
as a post-fix rerun of the entire ChEMBL 37 corpus. RDKFingerprint and Avalon,
by contrast, were subsequently rerun across the entire corpus after their
source ports were completed.

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

The 29 ChEMBL branches cover molecular weight, exact molecular weight, formula
variants, H-bond donor/acceptor counts, fraction Csp3, Crippen LogP/MR option
branches, TPSA with and without S/P contributions and force recomputation,
aromatic-ring count, four rotatable-bond modes, and QED. Float descriptors use
exact RDKit bit-pattern comparison; integer and string outputs use exact
equality. QED is pinned to `RDKit 2026.03.1 + CPython 3.13.12` because RDKit
delegates reductions to Python `sum()` and older CPython versions use a
different floating-point reduction algorithm.

### Fingerprints

Morgan coverage includes 15 parameter profiles. Each profile compares the main
vector and `AdditionalOutput`, giving 30 configured output branches; MACCS adds
one. The large-run counter contains 89,831,720 completed comparisons because
102 records in each of the two custom-bond-invariant outputs were not recorded
as completed branch comparisons by the audit harness. The source fix was
validated across the retained fingerprint set and all affected branches.

The dedicated RDKFingerprint/Avalon full-corpus audit processed all 2,897,819
source records: 2,897,804 were mutually parseable and 15 were rejected by both
libraries. Every mutually parseable molecule matched exactly across 14
topological vector profiles (40,569,256 comparisons), 23 Avalon explicit-bit
profiles (66,649,492 comparisons), and two topological provenance profiles
(5,795,608 comparisons). The provenance comparison includes vector width,
complete on-bit sets, `atomBits`, and `bitInfo`. All 128 shards completed with
zero mismatches and no failed task, for 113,014,356 exact comparisons total.
The 2026-08-19 run used COSMolKit `0.2.12` at commit `ea7c581c`, pinned RDKit
`2026.03.1`, and ChEMBL 37 source SHA-256
`ea6181ce8dc7af41974e35b92e1febb0c9dcbe2c62f7ccc4a5d983ac19f696e7`.

The maintained 5,000-row corpus continues to provide committed exhaustive
regression gates for the same 14 topological and 23 Avalon profiles. Sparse-
count, sparse-bit, hashed-count, explicit-bit, raw/public projection, and
provenance output are claimed only where enumerated by each family's active
profile. No approximate vector, similarity correlation, or 99.9% bit
agreement is accepted as parity.

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

The final 2026-08-20 execution completed these four phases across every
configured shard. Fragments and direct 2D coordinates had no blocking
mismatch. Topology operations exposed the unsanitized `RemoveHs` intermediate
valence/H-state difference, and binary roundtrips exposed post-deserialization
hash/fingerprint differences. Because the run-level gate failed, the results
remain an audit record rather than accepted zero-mismatch evidence. They can be
promoted only after source-level fixes, focused regressions, complete affected-
phase reruns, and a passing complete-profile identity.

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
signed volume against seeded RDKit ETKDG geometries. Substructure matching
currently compares molecule-query atom mappings; direct SMARTS-query parity is
outside that specific claim.

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
