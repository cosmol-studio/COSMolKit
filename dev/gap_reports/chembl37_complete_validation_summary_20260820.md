# ChEMBL 37 Complete Validation Summary — 2026-08-20

## Overall Verdict

This document is the consolidated ChEMBL 37 validation record for COSMolKit.
It combines the original complete chemistry audit, the subsequent complete
RDKFingerprint/Avalon audit, and the final 15-phase execution that closes the
new topology, fragment, direct-2D, binary-roundtrip, composition, conformer,
force-field, and I/O boundaries.

All recorded corpora, tasks, and artifacts passed integrity validation. The
overall parity gate is nevertheless **not yet accepted**: unsanitized hydrogen
removal and binary deserialization expose two reproducible systemic gaps. A
complete summary must retain those failures rather than report only the clean
surfaces.

## Evidence Sets

The executions below share the pinned ChEMBL 37 source and RDKit reference but
have separate immutable identities. Counts from overlapping executions are not
added together and presented as distinct molecules.

| Evidence set | Corpus/task boundary | Result |
|---|---|---|
| Complete chemistry, composition, batch, and concurrency audit | 2,897,819 source records; 22 phases; 2,816 shard tasks | Findings were source-traced and replayed over every affected retained branch; this is not represented as a complete post-fix corpus rerun |
| Complete RDKFingerprint/Avalon audit | 2,897,819 source records; 128 shard tasks; 113,014,356 exact comparisons | Pass, zero mismatches and zero failed tasks |
| Complete 2026-08-20 validation execution | 15 phases; 1,920 shard tasks; 31,755,360 phase-record executions | 13 phases pass; topology operations and binary roundtrip fail |

The final table below gives the complete current surface status:

| Surface group | Current evidence status |
|---|---|
| Molecular descriptors, Morgan/MACCS, explicit hydrogen, fragments, direct 2D coordinates, bounds, InChI parse, SMILES matrix, SVG, I/O, scalar/batch composition, conformers, and UFF/MMFF | Complete configured 2026-08-20 phase passed with zero blocking mismatch |
| RDKFingerprint and Avalon | Separate complete ChEMBL 37 audit passed with zero mismatch |
| Earlier core, InChI-option, operation-order, and shared-read matrices | Complete original audit plus source-fixed retained-branch replay; no complete post-fix corpus rerun is claimed |
| `RemoveHs(sanitize=false)` intermediate state | Blocking source-parity gap |
| COSMolKit binary roundtrip behavior | Blocking serialization-state gap |

## Final Execution Identity

| Field | Recorded value |
|---|---|
| Run directory | `tmp/parity-audit/chembl37-postfix-full-20260820` |
| Started | `2026-08-20T06:30:02.813488+00:00` |
| Finished | `2026-08-20T09:21:51.580389+00:00` |
| COSMolKit version | `0.2.12` |
| Git commit | `8d285f9d9dbfd453a2fb91bc5b33c60a8aa1fffd` |
| Installed extension SHA-256 | `e5fb253cae3291e32c67e0efc973251281560a96e2285322bf39af80b621736d` |
| Recorded tracked-diff SHA-256 | `5bcf70f5957ece60154044f3cf42f5106be77606a3c620adb8c9aeed96ec63c6` |
| RDKit | `2026.03.1` |
| CPython | `3.13.12` |
| ChEMBL source records | `2,897,819` |
| Mutually parseable sanitized records | `2,897,804` |
| ChEMBL source SHA-256 | `ea6181ce8dc7af41974e35b92e1febb0c9dcbe2c62f7ccc4a5d983ac19f696e7` |
| Corpus manifest SHA-256 | `f8b1a516a7794c8d4428a8309c80c049b5c4f34c6c9c54381f0aebccd9ccb976` |
| Workers | `112` |

The identity retains the repository status, tracked diff, installed extension,
profile, scripts, reference profiles, interpreter, NumPy, RDKit, and corpus.
The installed extension still hashes to the recorded value. The working tree
changed after execution began, so this identity validates the recorded
extension and repository snapshot, not later uncommitted Element, Residue,
dependency, or XYZ changes.

## Completeness And Integrity

| Check | Result |
|---|---:|
| Phases completed | 15 / 15 |
| Successful shard tasks | 1,920 / 1,920 |
| Missing or failed tasks | 0 |
| Processed phase-record executions | 31,755,360 / 31,755,360 expected |
| Output summaries with valid schema | 1,920 / 1,920 |
| Output/finding checksum failures | 0 |
| Missing output/finding artifacts | 0 |
| Passing `match.*` counters | 2,500,870,740 |
| Blocking `mismatch.*` counters | 7,804,836 |
| Distance-bounds matrix entries traversed | 2,757,910,995 |

The corpus was independently revalidated with per-shard checksums. Every task
output and bounded finding file was independently rehashed against the run
manifest after completion.

## Complete Phase Results

`Checks matched` and `checks mismatched` are sums of each phase's named
`match.*` and `mismatch.*` counters. They count branch/state observations, not
distinct molecules or distinct implementation defects.

| Phase | Records processed | Checks matched | Checks mismatched | Result |
|---|---:|---:|---:|---|
| Molecular descriptors | 2,897,819 | 84,036,316 | 0 | Pass |
| Morgan/MACCS | 2,897,819 | 89,831,720 | 0 | Pass |
| Explicit hydrogen | 2,854,362 | 11,417,448 | 0 | Pass |
| Topology operations | 2,854,376 | 39,962,026 | 5,707,822 | **Fail** |
| Connected fragments | 2,854,362 | 5,852,676 | 0 | Pass |
| Direct 2D coordinates | 2,854,362 | 14,271,810 | 0 | Pass |
| Binary roundtrip | 524,288 | 9,437,322 | 2,097,014 | **Fail** |
| Distance bounds | 2,854,362 | 2,854,362 | 0 | Pass |
| InChI parse branches | 2,854,362 | 11,417,448 | 0 | Pass |
| SMILES writer matrix | 2,854,362 | 2,192,150,016 | 0 | Pass |
| Final SVG | 2,854,362 | 2,854,362 | 0 | Pass |
| Molecular I/O | 524,288 | 12,559,920 | 0 | Pass |
| Scalar/batch composition | 1,027,660 | 20,561,265 | 0 | Pass |
| Fixed-seed conformer outcome | 524,288 | 524,288 | 0 | Pass |
| UFF/MMFF force fields | 524,288 | 3,139,761 | 0 | Pass |

The conformer phase contains 421,217 matching coordinate outcomes and 103,071
matching failure outcomes. The force-field phase covers UFF/MMFF parameter
availability with and without explicit hydrogens and both optimizers. These
phases completed without a blocking mismatch.

## Blocking Gap: Unsanitized Hydrogen Removal

The topology mismatch is confined to `AddHs` followed by
`RemoveHs(sanitize=false)`:

| Branch | Matches | Mismatches |
|---|---:|---:|
| Value-style graph state | 451 | 2,853,911 |
| In-place graph state | 451 | 2,853,911 |

Every one of the 3,072 retained examples differs only in the atom-level
`explicit_valence`, `implicit_hydrogens`, and `total_hydrogens` fields. The
retained examples do not show a separate bond, atom-count, element, charge,
aromaticity, or chirality difference. The same field set was reproduced after
the completed execution.

This is one systemic intermediate-state parity gap expressed over many rows,
not 5.7 million unrelated bugs. It remains blocking because `sanitize=false`
is an explicitly tested public branch and the comparison contract preserves
unsanitized valence/H state rather than normalizing it away.

## Blocking Gap: Binary Roundtrip

The binary phase serializes a COSMolKit molecule with 2D coordinates, restores
it through each public deserialization entrypoint, and compares behavior before
and after. It does not compare COSMolKit bytes with RDKit bytes.

| Observation | Class method | Module function |
|---|---:|---:|
| Hash mismatches | 524,288 / 524,288 | 524,288 / 524,288 |
| Morgan mismatches | 524,219 / 524,288 | 524,219 / 524,288 |
| Graph-state matches | 524,288 / 524,288 | 524,288 / 524,288 |
| Deterministic-byte matches | 524,288 / 524,288 | 524,288 / 524,288 |
| Conformer-count matches | 524,288 / 524,288 | 524,288 / 524,288 |

All retained hash examples are successful before serialization and return the
same structured valence/CIP error after restoration. Retained Morgan examples
run successfully on both sides but produce different on-bit sets. A post-run
reproduction confirmed that graph state and reserialized bytes remain equal
while hash fails and Morgan changes.

This is a serialization-state completeness defect. Byte determinism and
user-visible graph equality do not make a roundtrip valid when supported
downstream operations change behavior.

## Final Acceptance State

- The complete evidence inventory is recorded here; no clean phase is omitted
  and no failed phase is hidden.
- The 13 clean final-execution phases and the separate complete fingerprint
  audit retain their exact counters and identities.
- No mismatch is reclassified as informational, unsupported, or an accepted
  error rate.
- The two systemic gaps require source-level fixes and focused regressions that
  preserve the failing branch/state fields.
- After rebuilding the Python extension, both affected phases must be rerun
  over all configured shards. A passing complete-profile identity is required
  before the overall zero-mismatch claim is accepted.
