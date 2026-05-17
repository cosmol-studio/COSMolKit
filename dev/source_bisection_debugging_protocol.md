# Source Bisection Debugging Protocol

This document defines the mandatory debugging workflow for source-level
port mismatches, especially RDKit parity or source-reproduction gaps.

This protocol exists to prevent heuristic patching. The goal is to locate
the **first divergent state boundary** between upstream C++ and COSMolKit
Rust, then patch only the corresponding source-backed logic.

---

## 1. Golden Rule

Do **not** debug from the final wrong output backward by intuition.

Do **not** patch the first Rust function that "looks suspicious".

Do **not** make a local fix before identifying the first upstream/Rust
state divergence.

The required workflow is:

1. Define the upstream call chain.
2. Choose observable state boundaries on that call chain.
3. Instrument upstream and Rust at the same boundaries.
4. Compare boundary states using bisection.
5. Find the **first** boundary that differs.
6. Narrow that boundary to a specific source block or source line region.
7. Patch only with explicit source evidence.
8. Add a permanent test at the smallest stable boundary.

---

## 2. What Counts As a Boundary

A valid debugging boundary is a point where upstream and Rust can expose the
same state in a comparable form.

Examples:

- raw parse result
- property-cache-updated temporary molecule
- ring info after ring finding
- bond direction state after direction normalization
- bond stereo state after stereo assignment
- CIP ranks after rank assignment
- canonical ranks before traversal
- traversal start atom
- canonical traversal stack
- final serialized output

Bad boundaries:

- "somewhere in canonicalization"
- "around the writer"
- "probably traversal-related"

Each boundary must be named concretely.

---

## 3. Required Debugging Shape

For each bug, define a state chain like:

1. raw input state
2. stage A output
3. stage B output
4. stage C output
5. final output

Then perform bisection on these states:

- if stage A matches and stage C differs, inspect stage B
- if stage B already differs, stop searching later stages
- never patch stage C when stage B is the first divergence

This is state bisection, not guess-and-check.

---

## 4. Upstream Instrumentation First

If the upstream C++ source is available and buildable, instrumentation should
start there.

Reason:

- the upstream call chain is the source of truth
- one upstream run can expose multiple boundary states
- it prevents Rust-side debugging from inventing nonexistent invariants

Preferred order:

1. instrument upstream at multiple checkpoints in one pass
2. capture all relevant boundary states
3. instrument Rust with the same checkpoint names
4. compare states

For RDKit-style issues, one instrumentation pass should usually include:

- raw parsed bond/atom state
- post-property-cache state
- post-stereo-preparation state
- post-ring-finding state when relevant
- post-canonical-ranking ranks
- traversal start atom
- traversal stack
- final output

---

## 5. Three-Probe Rule

For a single bug, instrumentation must converge quickly.

One probe means:

- one upstream/Rust comparison round with a defined checkpoint set
- or one refined follow-up probe after the previous probe reduced the range

The expected discipline is:

1. Probe 1
   Produce a coarse call-chain state map and reduce the bug to a bounded
   stage interval.

2. Probe 2
   Split that interval and reduce the bug to a specific function, helper, or
   source block.

3. Probe 3
   Reduce the bug to a concrete source-backed logic difference or line region.

If Probe 3 still does not identify the first divergent source block, the
debugging plan is wrong and must be rewritten before any further code changes.

Do **not** accumulate five, six, or ten speculative probes.

---

## 6. Mandatory Reporting After Each Probe

After every instrumentation round, report:

1. the checkpoints inspected
2. the last checkpoint that still matched
3. the first checkpoint that differed
4. the reduced problem range
5. the next narrower probe to run

Required format:

```text
Probe N
- Checked: <checkpoint list>
- Last matching boundary: <boundary>
- First differing boundary: <boundary>
- Reduced range: <from X to Y>
- Next probe: <specific narrower instrumentation plan>
```

Bad report:

```text
Maybe writer-prep is wrong.
```

Good report:

```text
Probe 1
- Checked: raw parse, post-setBondStereoFromDirections, post-updateDoubleBondStereo, canonical ranks
- Last matching boundary: post-setBondStereoFromDirections
- First differing boundary: post-updateDoubleBondStereo
- Reduced range: RDKit updateDoubleBondStereo / FindStereo cleanup vs Rust writer-prep stereo assignment
- Next probe: instrument controlling atoms and candidate-bond filtering for the target double bond
```

---

## 7. No Heuristic Patches

A code change is forbidden if it is justified by:

- "this seems likely"
- "this probably matches the output better"
- "this branch looks related"
- "this makes the failing test pass"

A code change is allowed only when justified by:

- an identified first divergent boundary
- an upstream source block corresponding to that boundary
- a concrete mismatch between upstream state transition and Rust state transition

If the change does not have a source-backed boundary explanation, it is
heuristic and must not be made.

---

## 8. Instrument Multiple Checkpoints Per Run

Do not waste probes by inspecting only one point when the same run can expose
an entire sub-chain.

Prefer one probe that prints:

- boundary A
- boundary B
- boundary C
- boundary D

instead of four isolated probes.

The purpose of instrumentation is not "more logs". The purpose is maximum
range reduction per run.

---

## 9. Patch Only After the Range Is Narrow

Patching is allowed only after the bug has been narrowed to:

- one function
- one helper cluster
- or one copied source block

If the current range still spans multiple conceptual stages, keep instrumenting.

Examples:

- acceptable: "difference is inside `updateDoubleBondStereo()` cleanup of
  non-potential stereo bonds"
- not acceptable: "difference is somewhere in stereochemistry or ranking"

---

## 10. Permanent Test Requirement

Every fixed bug must leave behind:

1. one smallest-boundary unit test for the newly understood state transition
2. one higher-level regression test if the bug was user-visible

The smallest-boundary test is mandatory because it prevents future regressions
from hiding behind a large end-to-end failure surface.

---

## 11. RDKit-Specific Recommended Checkpoint Vocabulary

Use stable checkpoint names so future debugging stays comparable:

- `raw_parse`
- `post_update_property_cache`
- `post_set_bond_stereo_from_directions`
- `post_find_potential_stereo`
- `post_update_double_bond_stereo`
- `post_assign_stereochemistry`
- `post_canonical_rank`
- `post_choose_start_atom`
- `post_canonicalize_fragment`
- `final_smiles`

When possible, use the same names in:

- upstream C++ instrumentation
- Rust debug output
- probe reports

---

## 12. Required Attitude

For source-level ports, debugging is a localization problem before it is an
implementation problem.

If localization is weak, implementation quality will drift into heuristic
repair.

Therefore:

- first localize
- then explain
- then patch
- then lock with tests

Never reverse that order.
