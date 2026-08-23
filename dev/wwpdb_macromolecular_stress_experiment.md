# wwPDB Macromolecular Stress Experiment

This document defines the reproducible large-macromolecule validation boundary
for COSMolKit. It is an experiment protocol, not evidence that the experiment
has already passed. Results may become release evidence only after the complete
profile finishes and its immutable run manifest is accepted.

## Objective

The experiment validates the public structural-biology and PDB/mmCIF molecule
workflows against a large, heterogeneous, single-source corpus. It is designed
to detect:

- Gemmi-port differences in PDB and PDBx/mmCIF parsing;
- hierarchy, source-identifier, coordinate, metadata, connection, assembly,
  and protein-projection differences;
- PDB/mmCIF path differences that disappear in small hand-written fixtures;
- RDKit-parity differences in the public PDB-to-`Molecule` and PDB writer
  profiles;
- composition failures across parse, projection, molecule conversion, PDB
  writing, reparsing, and water removal;
- panics, aborts, unbounded per-entry memory growth, pathological long tails,
  nondeterminism, and concurrency-only differences.

This is a valid-archive stress experiment. Malformed-input testing and parser
fuzzing are separate boundaries.

## Current Public Boundary

The experiment must test what COSMolKit currently exposes and keep the
structural writer boundary explicit.

Supported compositions in scope are:

```text
PDB text   -> BioStructure -> Protein
mmCIF text -> BioStructure -> Protein
PDB text   -> BioStructure -> Molecule
mmCIF text -> BioStructure -> Molecule
Molecule   -> PDB text
BioStructure -> mmCIF text -> BioStructure
BioStructure.without_waters() -> BioStructure
```

The project does not expose `BioStructure.to_pdb()`, `Protein.to_pdb()`, or
`Protein.to_mmcif()`. `BioStructure.to_mmcif()` is implemented and belongs to
the complete structural model. Consequently, “PDB/mmCIF conversion” in this
experiment has three precise meanings:

1. PDB and mmCIF representations of the same wwPDB entry are parsed and
   compared through their common structural projection, relative to pinned
   Gemmi.
2. The implemented `PDB -> BioStructure -> mmCIF -> BioStructure` composition
   is compared through complete structural state and the pinned Gemmi writer.
3. The implemented `mmCIF -> Molecule -> PDB -> Molecule` and
   `PDB -> Molecule -> PDB -> Molecule` compositions are exercised.

Direct `BioStructure` PDB-to-mmCIF writing is exercised through the public
composition above. PDB output remains a `Molecule` compatibility writer; there
is no `BioStructure.to_pdb()`.

## Corpus Identity

The corpus is the official **wwPDB 20260101 annual archive snapshot**. The
[wwPDB archive documentation](https://www.wwpdb.org/ftp/pdb-ftp-sites)
describes annual snapshots as stable datasets for research, and the
[RCSB snapshot inventory](https://s3snapshots.rcsb.org/) records this snapshot
as 246,905 experimentally determined coordinate entries current on
2026-01-01. The complete snapshot is approximately 1,583 GB; the experiment
prepares only the coordinate and holdings data it needs.

Canonical source prefixes are:

```text
s3://pdbsnapshots/20260101/pub/pdb/data/structures/divided/mmCIF/
s3://pdbsnapshots/20260101/pub/pdb/data/structures/divided/pdb/
s3://pdbsnapshots/20260101/pub/pdb/holdings/
```

The equivalent HTTPS snapshot endpoint is
`https://pdbsnapshots.s3.us-west-2.amazonaws.com/`. Current archive URLs such
as `https://files.rcsb.org/download/4hhb.cif` are convenient but must not be
used as corpus identity because current entries can be revised. The
[wwPDB versioned archive](https://www.wwpdb.org/ftp/pdb-versioned-ftp-site)
may be used to independently recover a specific entry revision.

PDBx/mmCIF is present for the complete coordinate corpus. Legacy PDB is a
paired subset. This distinction is required because the legacy format cannot
represent all current structures; the
[official RCSB format-limit documentation](https://www.rcsb.org/docs/general-help/structures-without-legacy-pdb-format-files)
lists chain-count, atom-count, identifier, and numeric-field boundaries.
Absence of a legacy PDB file is therefore corpus metadata, not a failed or
skipped mmCIF case.

The run manifest must separately count mmCIF entries rejected because a chain
or subchain identifier exceeds the current four-byte `PdbChainId` model. That
declared model boundary is not a matching result and cannot be silently
truncated or excluded from corpus accounting. Entity source identifiers are
variable-length strings and are not part of this rejection class.

## Reproducible Preparation

Preparation must run outside tracked repository data and must be atomic.
Repository code may populate a temporary sibling directory and rename it only
after every check passes.

The preparation manifest must contain:

- snapshot name and source prefixes;
- the expected 246,905 current coordinate entries;
- every entry identifier and available input formats;
- compressed and uncompressed byte counts and SHA-256 for every selected file;
- a sorted manifest SHA-256 and a Merkle root over all file identities;
- exact Gemmi and RDKit reference identities;
- the preparation-script and profile SHA-256 values;
- every derived diversity label and the reference values used to assign it;
- deterministic subset membership and shard assignment;
- rejected, missing, corrupt, and duplicate input counts.

Preparation must fail if an expected mmCIF entry is missing, a checksum
changes, an identifier occurs twice, decompression fails, or the discovered
entry count differs from the snapshot identity. It must not download from the
current archive to repair a snapshot hole.

Inputs, prepared shards, Gemmi expected data, and run results remain outside
Git. The repository must commit the preparation code, schemas, profile,
reference pins, and final aggregate report before the experiment can be
called reproducible.

## Reference Implementations

Two references have distinct ownership boundaries:

- Gemmi commit `5cc1c23c6007e0e6cbd69289c6f7c0bff50e943e`, matching the
  repository submodule pin, is the structural PDB/mmCIF parser, hierarchy, and
  protein-projection reference.
- RDKit `2026.03.1` is the PDB-to-`Molecule` and `Molecule`-to-PDB reference.
  RDKit has no direct `MolFromMMCIFBlock` oracle, so mmCIF molecule conversion
  must not be labeled RDKit parity.

The reference adapter must emit a versioned, documented JSON schema. It must
not normalize a field merely because COSMolKit differs. Adapter changes alter
the run identity and invalidate cached expected data.

The run manifest also records the COSMolKit Git commit and status, installed
Rust and Python artifact hashes, Cargo feature set, interpreter, reference
package hashes, operating system, CPU, worker count, and profile hash. A dirty
tree requires retaining the exact tracked patch.

## Diversity Model

All complete structural parsing phases use the complete applicable corpus.
The diversity model is used to define smaller preflight suites and the
expensive molecule/writer/composition boundary; it is not a filter that can
remove a failing full-parser row.

Labels are derived by pinned Gemmi from the snapshot mmCIF, never from
COSMolKit output. Each entry may have multiple labels.

Content labels include:

- protein-only, nucleic-acid-only, mixed polymer, branched entity, and
  non-polymer content;
- ligands, waters, metals, modified amino acids, and unknown components;
- single-model and multi-model structures;
- alternate locations, insertion codes, deuterium, formal charges, and
  missing coordinates;
- explicit connections, cis-peptides, modified-residue records, crystal
  metadata, NCS operators, and assemblies;
- X-ray, electron microscopy, NMR, neutron, and other deposited methods.

Stress labels include bins for compressed size, uncompressed size, atom count,
residue count, chain count, model count, entity count, connection count, and
assembly count. Bin boundaries are stored in the profile rather than inferred
from a particular run.

The deterministic expensive-path set is the union of:

1. the 4,096 entries with the smallest
   `BLAKE2b-128("wwpdb-20260101\0" + entry_id)` value;
2. the 128 lowest-digest entries for every content label;
3. the 128 lowest-digest entries for every populated stress bin;
4. the largest 512 entries by each of atom count, uncompressed bytes, model
   count, chain count, connection count, and assembly count; and
5. the corresponding paired-format entries needed by a selected composition.

Ties are resolved by the digest and then the ASCII entry identifier. The
manifest records the final union and count. A smaller preflight profile uses
the same algorithm with 256 base entries, 16 per label/bin, and 32 per extreme
metric. Only the complete profile may update public validation claims.

## Sharding And Scheduling

Entries are assigned by
`BLAKE2b-64("wwpdb-20260101-shard\0" + entry_id) modulo 256`. A shard is the
atomic resumable task. Within a shard, entries are processed in ASCII entry-ID
order.

Every worker processes one entry at a time. This prevents a worker from
retaining an archive-sized collection and makes peak memory attributable to a
specific entry. Reference and COSMolKit measurements run in separate worker
processes so allocator reuse in one implementation does not contaminate the
other measurement.

An entry timeout is an experiment failure, not an exclusion. An interrupted
run is incomplete even if every finished shard is clean. Resume is allowed
only when corpus, profile, code, binary, reference, and environment identities
are identical.

## Comparison Schema

### Structural state

The Gemmi differential schema compares ordered models, chains, residues,
atoms, coordinates, entities, metadata, crystal state, resolution,
connections, cis-peptides, modified residues, assemblies, and every source
identifier exposed by the corresponding COSMolKit public rows.

Integer, enum, optional, string, ordering, and cardinality fields compare
exactly. For source values stored by COSMolKit as `f32`, the reference value is
converted once to IEEE-754 binary32 and the resulting bits compare exactly.
Values stored as `f64` compare as binary64 after the same source-defined
conversion. No decimal tolerance is introduced to hide representation
differences.

PDB and mmCIF are each compared with Gemmi independently. Cross-format state
uses a documented common projection keyed by deposited model, author chain,
author residue sequence/insertion code, residue name, atom name, and alternate
location. A direct PDB-versus-mmCIF mismatch is blocking only when the same
relation is not present in Gemmi. This reference-differential rule preserves
legitimate precision and representability differences in the legacy PDB
format.

### Protein projection

The `Protein` schema compares:

- file and in-memory constructors for PDB and mmCIF;
- model, chain, residue, and atom counts;
- chain indexing, negative Python indexing, iteration order, and lengths;
- residue raw name, `ResidueCode`, kind, source-derived information, FASTA
  code, canonical one-letter code, parent standard code, and modified status;
- atom index, name, element, and optional Cartesian position;
- equality of flattened traversal with nested chain/residue traversal; and
- exclusion of ligands, nucleic acids, and waters according to the documented
  protein-only projection.

The Python list-returning APIs are exercised on the extreme-size set as well
as ordinary entries so clone or ownership costs are visible rather than hidden
by scalar-count-only tests.

### Molecule state

PDB-to-`Molecule` profiles compare COSMolKit with RDKit across the Cartesian
product of:

```text
sanitize:          false, true
remove_hs:         false, true
proximity_bonding: false, true
reader flavor:     0, 1, 8, 9
```

The comparison retains success/failure class, complete atom and bond state,
PDB residue information, conformer count, coordinates, and structured
diagnostics. A sanitizer rejection from both implementations is a matching
outcome; an exception, panic, missing result, or different graph is not.

mmCIF-to-`Molecule` uses the same profile matrix and state schema, but its
oracle is the Gemmi-derived structural state plus COSMolKit's documented
conversion invariants. Paired PDB/mmCIF entries additionally compare the
relation between both conversion paths. These results must be called
cross-format consistency, not RDKit mmCIF parity.

### PDB writing and compositions

For every selected structural PDB profile, `BioStructure.to_mmcif()` output is
compared with the pinned Gemmi writer at exact bytes, then reparsed by both
implementations and compared through complete represented structural state.
Category suppression and every public `MmcifWriteOptions` formatting branch
are included in the deterministic writer phase.

Every successful selected PDB molecule profile is written with `conf_id=-1`
and `conf_id=0` across all 64 combinations of the six documented RDKit PDB
writer flavor bits. COSMolKit and RDKit text compare byte-for-byte. The output
is then reparsed by the same implementation and compared through complete
molecule state.

Composition phases are:

```text
PDB   -> Molecule -> PDB -> Molecule
mmCIF -> Molecule -> PDB -> Molecule
PDB   -> BioStructure -> Protein traversal
mmCIF -> BioStructure -> Protein traversal
BioStructure -> without_waters -> hierarchy validation
BioStructure -> mmCIF -> BioStructure -> hierarchy and metadata validation
```

The value-style `without_waters()` result must preserve the source object,
remove exactly the Gemmi-classified water residues, compact every hierarchy
span and coordinate row consistently, and satisfy BioStructure invariants.

## Execution Phases

| Phase | Corpus boundary | Required oracle or invariant |
|---|---|---|
| `inventory` | all 246,905 entries | checksums, formats, labels, and shards complete |
| `mmcif_biostructure` | all mmCIF entries | complete Gemmi structural state |
| `pdb_biostructure` | every paired legacy PDB entry | complete Gemmi structural state |
| `paired_common_projection` | every paired entry | COSMolKit relation equals Gemmi relation |
| `protein_projection` | all mmCIF plus all paired PDB | complete public Protein traversal schema |
| `without_waters` | all mmCIF entries | source preservation, exact selection, valid compacted hierarchy |
| `pdb_molecule` | deterministic expensive-path paired set | 32 RDKit reader profiles |
| `mmcif_molecule` | deterministic expensive-path set | 32 documented conversion profiles |
| `pdb_writer` | every successful selected PDB molecule outcome | 128 writer outputs per outcome and RDKit exact text |
| `biostructure_mmcif_writer` | every selected structural PDB/mmCIF outcome | exact Gemmi writer bytes, semantic reparse, option matrix |
| `format_composition` | deterministic expensive-path set | complete post-roundtrip state relations |
| `concurrency_determinism` | every diversity label and extreme metric | 1-worker and configured-worker digests identical |
| `resource_profile` | every phase | per-entry time, CPU, allocation/RSS, and output bytes retained |

Phase counters count actual compared fields, rows, coordinates, graphs, and
writer outputs. A record count alone is not a comparison-strength claim.

## Performance And Resource Evidence

The first complete run establishes the measured performance baseline; it must
not invent a pass threshold after seeing outliers. It records, by phase and
size bin:

- total and per-entry wall and CPU time;
- compressed and uncompressed input throughput;
- median, p90, p95, p99, p99.9, and maximum latency;
- peak worker RSS and maximum single-entry RSS delta;
- atoms and coordinate rows processed per second;
- COSMolKit/Gemmi and COSMolKit/RDKit ratios on the same inputs;
- the complete slowest and highest-memory entry lists; and
- one-worker versus configured-worker throughput and output digests.

No valid entry may disappear from correctness accounting because it exceeds a
resource budget. Timeout, process termination, panic, abort, or allocation
failure remains a named failing outcome. After an accepted baseline exists, a
later profile may add predeclared regression thresholds bound to that baseline
identity.

## Findings And Acceptance

Every entry in every configured phase must end in exactly one terminal class:

```text
matched
reference_rejected
cosmolkit_rejected
mismatched
timed_out
worker_failed
```

The valid archive is expected to parse; a reference rejection is retained and
investigated, not silently skipped. A capability advertised as supported may
not turn a failing row into “unsupported”. Unsupported may be reported only
for a separately named public capability outside the experiment boundary.

A complete profile passes only when:

1. all configured shards and phases finish;
2. every input and output checksum validates;
3. every covered correctness counter has zero unexplained mismatch;
4. no worker panics, aborts, times out, or disappears;
5. single-worker and concurrent output identities are equal;
6. every retained finding names the entry, phase, profile, first differing
   field, COSMolKit value, reference value, and artifact hashes; and
7. the aggregate report distinguishes complete-corpus phases from
   deterministic expensive-path phases.

Finding examples may be bounded per mismatch signature, but the complete list
of failing entry IDs and counters must be retained. Fixes require source-level
analysis and focused regressions. A retained-case replay is not a replacement
for rerunning the complete affected phase.

## Required Artifacts

The eventual repository-owned implementation belongs under
`dev/tools/wwpdb_stress/` and must provide one documented preparation command,
one profile-driven resumable runner, versioned schemas, reference adapters,
and unit tests for corpus identity, selection, sharding, resume, aggregation,
and failure propagation.

An accepted run publishes a concise report under `dev/gap_reports/` containing:

- immutable corpus, code, binary, reference, and environment identities;
- phase completeness and exact comparison counters;
- cross-format coverage and paired-entry count;
- all mismatch signatures and terminal-class counts;
- performance distributions and the complete long-tail explanation; and
- an explicit statement of which public claims the run does and does not
  support.
