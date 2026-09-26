# InChI static source audit assignment

## Authority and workspace

Work only in the checkout on the dedicated `audit/inchi-source` branch.
This is a static audit plus deliberate panic insertion, not implementation,
repair or validation execution.
The audited tree becomes a quarantined, intentionally non-runnable audit tree.
Do not merge, synchronize, publish, commit, stash, reset, or operate on siblings.
No Git commands, including read-only Git commands, are authorized.

Read AGENTS.md, dev/README.md, dev/crate_architecture.md,
dev/policy_invariants.md, dev/source_reproduction_protocol.md,
dev/source_bisection_debugging_protocol.md, dev/agent_plan_standard.md and
the InChI-specific instructions/port documentation before beginning.
Respect approved behavioral differences backed by explicit evidence.
Do not treat old completion claims or inline source comments as proof.

## First deliverable: complete function-level plan

Start in the terminal client's actual /plan mode. During that mode, inspect
read-only and propose the complete function-level plan; do not insert panics
or assume a chat instruction overrides Plan mode's mutation restrictions.
Persist the proposed inventory/plan when the client permits it; otherwise
return the complete plan for the mode's normal approval transition. The
execution requirements below apply only after implementation mode is enabled.

Before scanning or changing implementation, write:
- dev/audits/inchi/plan.md: complete numbered Read/action plan.
- dev/audits/inchi/inventory.md: full function map.
- dev/audits/inchi/findings.md: findings ledger.

Add only a short audit delegation pointer to this workspace's
dev/plans/crate_architecture_completion_plan.md; do not alter implementation
checkmarks or create an alternative migration architecture/queue.

Inventory the complete already-ported production InChI implementation:
all crates/cosmolkit-inchi/src files recursively, public entrypoints,
private helpers, impl methods, trait defaults, closures containing algorithms,
macro-produced implementations, generated source types, source constants/tables,
and any existing toolkit adapter or InChI-specific runtime body elsewhere.
Inspect actual code/reference mappings; do not assume shallow directory lists
or a grep of public fn is a complete inventory. Tests/examples are evidence,
not production implementation; missing unported functionality is separate.

For each function, record a stable ID, exact Rust path/symbol/range, upstream
file/symbol/range and pinned version identity, relevant defines/platform
assumptions, dependencies, and audit status. One-to-many and many-to-one
mappings must be explicit. Distinguish Rust-only plumbing from missing anchors.
Record table/type coverage separately so algorithmic state cannot escape review.
Freeze total function count and plan denominator before the first audit action;
newly discovered functions require an explicit inventory/plan amendment.

Explicit task-specific override (2026-09-26): read dev/policy_invariants.md and
dev/source_reproduction_protocol.md fully at startup and after context
compaction, not before every function, region, finding, or panic insertion.
This overrides the repeated Read/action pattern for this audit only; all
substantive source and evidence requirements remain. Remove pending repetitive
policy Read steps from the execution plan, preserving completed history and
stable function/finding IDs, and recompute step totals and active pointers.
Read the actual relevant upstream/Rust function bodies for every audit action.
Do not replace policy rereads with thousands of synthetic acknowledgment steps.
Display progress
as e.g. INCHI-AUDIT Step 3/N [x], replacing N with the actual total. Every
function has its own audit action and source mapping; never bundle an entire
module or dozens of helpers into one step. For large functions, additionally
enumerate bounded source control-flow regions, then a whole-function closure
check. Order callees/shared definitions before consumers where possible and
identify recursive groups. A read refresh is not audited function completion.
Do not pre-check future steps or batch-complete an inventory as an audit.

After proposing the complete plan, briefly restate scope and counts and use
the client's normal Plan-mode approval transition. Do not leave Plan mode
or start code changes merely because the earlier assignment said to begin.

## Line-by-line audit

Read full actual pinned upstream bodies and their relevant helpers/macros,
not only pasted anchors, summaries or an output example. Confirm official
InChI identity from available pin/version evidence without Git. Inspect the
pinned RDKit adapter only where the CK adapter actually derives from it.
Do not fetch a different source version or modify third_party.

For each function examine every source branch and Rust counterpart: predicates,
iteration and ordering, indices and sentinel values, integer widths/signedness/
overflow, pointer/aliasing translations, allocations and failure propagation,
state updates, lifetime/ownership, error/status/message/log/AuxInfo fields,
stereo/isotope/H/charge behavior, options, conditional compilation and helper
dispatch. Review complexity/allocations independently; matching-looking syntax
does not establish source parity. Distinguish upstream undefined behavior and
already approved departures from newly suspected deviations.

Both axes are mandatory acceptance gates: source behavior/control-flow/state
correspondence AND performance/complexity no worse than the original in the
modeled input space. Review worst-case and amortized time, auxiliary/peak
space, allocations/reallocations, temporary buffers, clones, repeated scans,
lookup/data-structure cost, recursion and hot-path constants. A matching
Big-O label alone does not justify extra full copies or avoidable allocations.
Concrete suspected regressions on either axis are findings requiring the same
panic quarantine during execution. No benchmarks may be run: distinguish
statically demonstrated regression, source-supported suspicion and unresolved
performance evidence; do not claim measured speed or equivalence without basis.

Write evidence for every audited function, including no-difference outcomes and
unresolved dependencies. A missing source prevents a verified result; record it
and continue independently auditable functions. Never infer parity from tests
not executed, or from tests that merely exist.

## Findings and required panic insertion

On a concrete source-supported POSSIBLE implementation difference, first record
a stable INCHI-AUDIT-xxxx ID in findings.md, then immediately insert a direct
panic! at the corresponding implementation location. Do not wait for a batch.
Record suspected vs confirmed-by-static-analysis separately, source excerpts,
Rust excerpts, trigger/path, expected source behavior, actual Rust behavior,
affected callers, behavior vs performance classification and exact panic site.
No runtime reproduction is claimed.

Use a literal marker such as:
panic!("INCHI-AUDIT-0001: suspected source divergence; see dev/audits/inchi/findings.md");

Keep original code and source anchors intact underneath/around the inserted
statement. No algorithm repair, rewritten branch, replacement return/default,
new error variant, weakened checker, cfg gate, debug_assert, or swallowed error.
Insert at the first divergent point within the EXISTING affected branch where
possible. Do not invent an approximate trigger predicate or panic an entire
public entrypoint when the suspect branch can be isolated. If divergence has
no safe existing branch boundary, guard the narrowest affected function entry
and explicitly record the broader quarantine scope. Inspect const/extern/unsafe
contexts statically and report special restrictions rather than changing ABI.
Do not insert into upstream sources, tests, signatures, or generated tables.
A table/type discrepancy is recorded at the declaration and guarded at the
narrowest existing consuming code site; record coverage limitations explicitly.

Panic insertion is the ONLY permitted production code modification. Restrict
it to InChI-owned implementation and actual InChI-specific adapters/bodies.
Other shared algorithms may be read as dependencies but not edited; if a
difference originates there, record it and place the quarantine at the
InChI-owned callsite. Missing tests/anchors alone are documentation/evidence
gaps, not automatically an algorithm defect deserving a fabricated trigger.
Already-approved differences must not be reclassified as unapproved defects.
Retain existing changes and all original behavior for later review.

## Absolutely no execution

Static file inspection/search/edit tools are allowed. Do not run project code,
upstream code, examples, binaries, benchmarks, parity/oracle generators, test
runners, cargo check/build/test/run, Python project scripts, build scripts,
code-generation tools or compilation. Do not run cargo fmt --all or rewrite
unrelated code; this task explicitly overrides normal post-edit execution
requirements. Match local formatting manually for inserted statements.
Do not install dependencies or change manifests/lockfiles/features.
Only static correspondence checks of the inventory, findings and panic IDs
are allowed. Report: compilation/tests NOT RUN by instruction.

## Continuation and reporting

Continue through all functions until the complete audit is done, a genuine
task-wide blocker occurs, or interrupted. A finding is not a reason to stop:
record it, insert the panic, and proceed. Persist the current function/region,
completed/remaining function counts, next exact source location, and outstanding
dependency evidence in the plan so compaction cannot erase the assignment.
Do not confuse read-step progress with function coverage. Final handoff must
report audited/total functions, finding and inserted-panic counts, unresolved
source gaps, all touched paths and the deliberate non-runnable status.
For a genuine blocker notify the supervising session with the function/step,
precise missing authority/source and evidence path; otherwise do not wait idle.
