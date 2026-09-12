# RDKit MolBlock Strict And Legacy-State Migration

## Outcome

`cosmolkit-io` now owns the detached syntax policy used by MolBlock and SDF
readers:

- `MolBlockReadParams { strict_parsing }` is public and defaults to RDKit's
  strict mode;
- `read_mol_block_detached_with_params()`,
  `read_v2000_detached_with_params()`, and
  `read_v3000_detached_with_params()` expose that policy without adding
  live-molecule chemistry options; and
- `SdfDataReadParams::strict_parsing` is passed through to the embedded
  MolBlock parse instead of controlling only trailing SDF data fields.

The `cosmolkit` block, file, forward-SDF, and indexed-SDF production readers
all call that detached owner and then install a validated concrete graph before
running the live-molecule chemistry-finalization policy. The former runtime
`MolFromMolDataStream` reproduction is test-only and is no longer reachable
from production IO.

## Source-Backed Coverage

The detached parser now follows the pinned RDKit `MolFileParser.cpp` branches
for:

- CTAB version validation and strict/non-strict fallback;
- the V2000 34-character strict versus 32-character non-strict atom-line
  boundary;
- non-strict unknown atom symbols as labeled dummy atoms;
- strict V3000 zero-valued outer atom/bond counts;
- non-strict blank lines before `M  END`;
- V2000 atom aliases, atom values, deprecated group-line skipping, and
  `S  SKP` records;
- old-style V2000 atom-list records, including negation and replacement by a
  later `M  ALS` record while preserving non-atomic query components;
- `M  PXA`, `M  ZCH`, `M  HYD`, `M  ZBO`, `M  APO`, and `M  LIN` lowering; and
- strict duplicate-attachment rejection versus non-strict first-value
  preservation;
- V3000 `RGROUPS` syntax, labels, isotope projection, and null-query state in
  `MolBlockRecord::Query`;
- V2000 field-level SGroup recovery across STY/SST/SLB/SCN/SDS,
  SAL/SBL/SPA, SMT/SDI/SBV/SDT/SDD/SCD/SED, SPL/SNC/SAP/SCL/SBT, with an
  invalid group dropped without discarding unaffected groups;
- V2000 short-SAP recovery by inferring the leaving atom from the sole
  crossing bond, or dropping the group when that inference is impossible;
  and
- V3000 invalid-group cleanup, crossing-bond validation, duplicate/count
  policy, and strict/non-strict handling of missing or unexpected SGROUP
  blocks and terminators.

V2000 and V3000 link nodes use COSMolKit's established
`_MolFileLinkNodes` property spelling so MolBlock, CXSMILES, and the runtime
adapter share one internal key. The pinned RDKit symbolic constant resolves to
`_molLinkNodes`; this internal spelling difference does not change the encoded
MolBlock/CX record.

The ZBO path preserves RDKit's exact bond-type transition: zero maps to the
RDKit zero bond, nonzero values map through the source bond-type enumeration,
and existing aromatic/query flags are not silently cleared. Property updates
move `AtomSpec` values rather than cloning their property maps.

## Explicit Residual Boundaries

- `M  MRV SMA` remains structured unsupported because recursive SMARTS belongs
  to `cosmolkit-search`, while `cosmolkit-io` intentionally has no dependency
  on that crate.
- Warning-log text is not modeled. Non-strict accepted-state results are
  preserved, and affected source markers remain behavior-axis partial where
  the warning side effect differs.
- The runtime intentionally owns live-molecule compatibility finalization
  (`sanitize`, `remove_hs`, and stereo handling), but not MolBlock syntax or
  graph assembly.
- Query-bearing records cannot be installed into the concrete runtime
  `Molecule`; this adapter returns a structured unsupported error. The same
  applies when attachment-point expansion would create RDKit null-query dummy
  atoms. Invalid attachment values remain preserved without expansion.

## Validation

```text
cargo test -p cosmolkit-io
75 passed, 0 failed

cargo test -p cosmolkit --features op-contracts-strict io::sdf::tests::
108 passed, 0 failed

cargo test -p cosmolkit --features op-contracts-strict io::molfile::tests::
11 passed, 0 failed
```

The regression set includes the pinned RDKit `H3BNH3.mol` ZCH/HYD/ZBO fixture
and focused strict/non-strict cases for every boundary listed above, including
old-style-only, negated, and conflicting old/new atom lists, malformed V2000
SGroup fields and attachment points, and malformed or structurally inconsistent
V3000 SGROUP blocks.
