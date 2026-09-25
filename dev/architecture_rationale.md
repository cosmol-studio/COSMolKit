# COSMolKit: Constrained Agent-Assisted Source Reproduction

COSMolKit aims to be a Rust-native cheminformatics toolkit with source-backed
behavior, efficient execution, and a coherent public API. Its architecture also
addresses a production problem: how to build and extend such a toolkit when
much of the implementation work is performed by agents operating at different
times, with different context, across many interacting domains.

The answer is not simply to require agents to be more careful. It is to arrange
the system so that many incorrect implementations are difficult to express,
unauthorized state changes are difficult to perform, and the remaining
structurally legal but semantically incorrect results can be investigated
through source-backed validation.

This explains why the model layer, algorithm crates, operation contracts,
private runtime, copy-on-write, strict checks, and parity tests belong in one
design. Individually, some can look redundant or unnecessarily elaborate.
Together, they address different failure modes while keeping that complexity
out of the user's everyday API.

This is a self-contained explanation of the design and its tradeoffs, not an
implementation-status report or an additional execution standard. The
authoritative documents are linked at the end.

## 1. Why a straightforward Rust rewrite is not enough

The simplest plausible port architecture puts a molecule in a foundational
crate and lets every algorithm receive it. Conceptually:

```text
                   core::Molecule
                          ^
       smiles / search / descriptors / conformers / force fields
```

An algorithm reads `&Molecule`; a transformation receives `&mut Molecule`.
Developers translate upstream code, run tests, and repair differential
mismatches. This has real advantages: low initial framework cost, familiar
interfaces, and a short path from source code to a working feature.

Its hidden assumption is that every implementer understands the entire
molecule's state model and reliably changes only the parts they should.
A topology edit may also require coordinate remapping, stereo-reference
updates, property handling, and cache invalidation. Each algorithm author
must remember the relevant obligations, including every error path.

A small team with stable ownership can sometimes sustain that arrangement.
Agent-assisted implementation makes the assumption much more fragile. A local
shortcut can look reasonable in one context while violating a decision made
elsewhere. Parallel work and context compression make that risk cumulative:
each function may appear defensible while the global state model gradually
loses consistency.

The problem is therefore larger than whether an agent can translate a
particular function correctly. It is whether hundreds of independently
implemented operations can coexist without each author understanding and
manually preserving the entire system.

COSMolKit optimizes for a smaller space of globally invalid implementations,
not merely the shortest path to each local feature.

## 2. Define authority before choosing crates

Asking where a function belongs is easier after asking who may do what.
COSMolKit separates four kinds of authority.

**Semantic authority** determines what a source-backed operation means.
Its basis is pinned upstream behavior, together with explicit project policies
and approved differences. An output that looks chemically plausible is not
enough, and a test expectation is evidence to investigate rather than a license
to invent a matching implementation.

**State authority** determines which representation is the live, user-visible
molecule. There is one `Molecule` lifecycle, owned by the public runtime.
Algorithm crates do not each keep their own competing live molecule model.

**Mutation authority** determines which state a particular operation may read
or change. Permission to alter topology does not imply permission to replace
coordinates, inspect arbitrary runtime caches, or rewrite unrelated properties.

**Commit authority** determines when a computed result becomes authoritative.
An algorithm may calculate a candidate, but successful calculation is not
itself permission to install that candidate into a live object.

These distinctions separate two questions that are often collapsed into one:
“Can this function calculate the result?” and “May this function change the
application's authoritative state?” An algorithm needs the first ability
without automatically receiving the second.

## 3. The crate graph makes that separation concrete

The main dependency direction follows those responsibilities:

```text
language bindings
       |
       v
cosmolkit: public Molecule + private runtime
       |
       +----> core and domain algorithm owners
                         |
                         v
              detached model values
                         |
                         v
                foundational vocabulary
```

The public runtime calls algorithms; algorithms operate on lower-level values.
No dependency back into `cosmolkit` is needed to run the chemistry. This makes
it possible to withhold live-state and commit authority from domain crates,
rather than merely asking them not to use it.

The split is not an argument for one crate per feature. Closely coupled
foundational chemistry can share `cosmolkit-core`; separately reusable domains
have their assigned owners. What matters is the ownership and dependency
boundary, not the number of packages.

This is also why `cosmolkit-model` exists. Algorithms need a common language
for atoms, bonds, adjacency, query graphs, coordinates, properties, and detached
transformation results. That language must not require the runtime that
eventually installs those results.

The model is therefore closer to a shared intermediate representation than a
collection of passive business objects. Local structural validation belongs
there; molecule lifecycle, cache validity authority, and commit policy do not.
Foundational vocabulary such as elements and bond-order enums belongs in
`cosmolkit-types`; domain-specific values retain their assigned owners.

The distinction is not “nouns go in model, verbs go elsewhere.” It is shared
representation versus domain interpretation versus live-state authority.

## 4. Why hiding a large internal object is insufficient

A natural alternative is a public wrapper around a complete internal object:

```rust
// Alternative design, not the COSMolKit algorithm boundary.
pub struct Molecule {
    inner: MoleculeData,
}
```

The public method delegates to an algorithm using `&self.inner` or
`&mut self.inner`. This gives users one attractive API, hides implementation
details, and avoids requiring extension traits for every domain.

Those are useful properties, but they solve API aggregation and information
hiding, not algorithm authority control. If `MoleculeData` exposes the whole
state to its algorithms, the wrapper protects it from users while leaving
implementation code broadly privileged.

COSMolKit needs to constrain that implementation code too. A helper that needs
only topology should not receive coordinates, caches, and runtime metadata
because handing over the complete object is convenient.

Splitting such helpers into separate crates does not fix the problem by
itself. An independent crate receiving a broadly mutable domain object can
still have too much authority. Crate boundaries and state-level capability
boundaries solve different problems.

## 5. Operation contracts turn intent into an interface

Every state-changing operation already has an implicit contract: what it reads,
what it changes, what it preserves, what it invalidates, and what relationship
must hold between input and output. Without an explicit representation, that
contract lives in comments, review habits, and the implementer's memory.

An illustrative topology operation might require:

```text
input:        writable topology, readable ordinary properties
preserve:     atom identity and coordinates
update:       selected derived state
output:       detached topology result and required effect evidence
```

The operation contract makes those obligations available to the framework.
Its value is not another checklist for an agent to remember. It is that the
declaration can drive the interface the operation body receives, the wrapper
that executes it, and the metadata against which it is validated.

Generated access markers expose operation-specific capabilities. If coordinates
are not authorized, an operation body does not get a general coordinate method
and a request to avoid calling it. The corresponding capability is unavailable.
Runtime-private unrestricted helpers remain outside the body's semantic module
boundary.

This is stronger than rejecting an unauthorized call at runtime: the call
should not compile in the real module layout. Compile-fail tests check that
distinction, while compile-pass tests ensure legitimate work remains possible.

Generating wrappers and contract matrices from the same declaration serves a
related purpose. A handwritten permission list, support table, and wrapper can
each be reasonable yet disagree. One declaration reduces those independently
maintained versions of the operation's meaning.

## 6. Narrow signatures bound the algorithm's world

Capabilities protect the operation body inside the runtime crate. Detached
function signatures extend that separation across the algorithm boundary.

For example, this schematic signature carries more information than an
ordinary calling convention:

```rust
fn calculate(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
) -> Result<Descriptor, Error>;
```

The function can inspect the supplied topology and properties. The interface
does not hand it coordinates, mutable topology, live cache storage, or commit
authority. A transformation can instead receive owned detached blocks and
return a typed transformation result.

In this sense, a signature defines the algorithm's available world. The
implementation works within explicit inputs rather than searching a large
object for whatever happens to be useful.

`OpParts` is the runtime-side working-state machinery, not a public molecule
decomposition API or an object passed to domain algorithms. Operation bodies
see its generated marker-specific surface; algorithm crates receive only the
authorized detached values. Keeping those two boundaries distinct prevents
a “narrow helper” from quietly carrying runtime authority into a lower crate.

Unavailable capability is more reliable than a convention not to use an
available one. This does not prove that an algorithm is correct, but it removes
whole categories of unrelated state access from ordinary implementation work.

## 7. A candidate result is not yet an authoritative result

Restricted access alone is not enough. An algorithm can have legitimate
permission to change topology and still return broken atom indices, an invalid
mapping, or coordinates that no longer correspond to its atom table.

If it works directly on the live object, an error halfway through can leave
that object in an intermediate state. Each implementation then needs its own
careful recovery logic, and every early return becomes part of the global
correctness argument.

Detached computation introduces a distinct candidate boundary:

```text
live input
    -> authorized working values
    -> domain computation
    -> candidate result
    -> applicable contract and invariant validation
    -> runtime-controlled installation
```

The domain algorithm computes; the runtime owns the transition to live state.
Mappings, effects, and structural relationships are explicit parts of that
transition, rather than unrelated cleanup that each algorithm may remember
or forget.

This is transaction-like, not a claim that every operation implements database
rollback. Value-returning operations preserve their source on failure.
Explicit in-place operations follow their documented error and panic policy,
which can permit completed partial changes while requiring complete, usable
storage. The common architectural point is that lifecycle and recovery policy
belong to the runtime, not to arbitrary domain helpers.

The candidate boundary also explains why merely returning `Ok` cannot be the
whole completion criterion. A successful computation can still violate the
operation's declared relationship to its input.

## 8. Isolation does not require copying the entire molecule

The obvious implementation of detached computation would deep-clone everything
before each operation. That would provide simple isolation, but at potentially
unacceptable cost for large molecules, many conformers, or repeated transforms.

COSMolKit separates logical isolation from physical copying through block-level
copy-on-write. Two value-style molecules can share unchanged storage. If an
operation modifies topology, it detaches the affected block as needed while
unmodified coordinate and property blocks can remain shared.

For the caller:

```rust
let next = molecule.with_hydrogens()?;
```

means that `molecule` is unchanged. It does not mean every byte of its state
must be copied immediately.

The location of COW matters. Returning an unrestricted `&mut Molecule` to an
algorithm in the name of efficiency would undo the authority boundary. Storage
sharing and detachment therefore remain inside runtime-controlled block
handling; algorithms still operate on their declared values.

COW is not just a performance feature added after the architecture. It makes
strong logical separation practical without forcing users to pay for full
copies whenever they choose predictable value semantics.

## 9. Value-style and in-place APIs share one design

Value semantics are easy to reason about, but some workloads benefit from
reusing uniquely owned storage. Maintaining separate algorithms for those
two cases would introduce another source of semantic drift.

Instead, eligible public forms use the same registered implementation:

```rust
let next = molecule.with_hydrogens()?;
molecule.add_hydrogens_()?;
```

The difference is the storage/lifecycle mode, not a second hydrogen algorithm.
The value form preserves its input and shares unchanged blocks. The in-place
form may reuse unique storage; if storage is shared, it still needs COW to
protect other values.

The trailing underscore makes mutation visible to the user. It does not expose
raw mutable storage or remove contract obligations. Nor does sharing the
algorithm mean the two forms promise identical failure guarantees: those are
part of their explicit lifecycle contracts.

Together, the candidate boundary and COW make value semantics, in-place
efficiency, and controlled authority compatible rather than competing goals.

## 10. Strict validation separates verification cost from use cost

Candidate validation raises another question: if each operation repeatedly
scans mappings, caches, coordinates, and preserved state, does the framework
become too expensive to use?

There are two different kinds of checking. Required production checks protect
the public contract and safe state handling. Additional strict checks provide
diagnostics, redundant verification, and evidence during development and CI.
They need not all have the same deployment policy.

Strict release-mode tests exercise the optimized implementation with those
additional checks enabled. Ordinary distribution builds can omit designated
diagnostic costs while retaining type boundaries, private visibility, COW,
required checks, and the same chemistry algorithm.

That last condition is essential. If disabling strict checks chooses a
different chemistry implementation, validation has tested a different program.
The purpose of strict mode is to observe and verify behavior, not to repair
results behind the scenes.

“Zero-cost abstraction” here is an aspiration about avoidable overhead, not a
claim that allocation, COW, validation, or runtime machinery costs nothing.
The useful distinction is between the cost of establishing confidence during
software production and the cost every user must pay during software use.

## 11. Three defenses against three different errors

The preceding layers cannot replace one another because they answer different
questions.

| Failure | Main defense | What it cannot establish |
|---|---|---|
| An implementation reaches state outside its authority | Crate direction, visibility, narrow inputs, generated capabilities | Correct chemistry within the permitted state |
| An authorized operation produces an invalid state transition | Mapping, effect, invariant and lifecycle checks | Equivalence to the reference algorithm |
| A structurally valid result has the wrong chemical meaning | Source review and pinned-reference parity | The absence of every untested defect |

For example, a coordinate algorithm should not gain authority to rewrite
topology merely because both belong to a molecule. A topology algorithm that
is authorized to remove an atom must still maintain the required mapping and
dependent state. Even after those checks pass, it might have selected the
wrong atom according to the reference semantics.

Rust can prevent many access and ownership mistakes. It cannot establish that
an aromaticity model or stereo assignment reproduces the selected upstream
behavior. Structural validity is necessary but not sufficient.

This is why source-backed reproduction and parity remain essential even in
a heavily constrained implementation.

### Verification follows ownership

Architecture assigns guarantees to shared mechanisms rather than asking each
algorithm to implement and prove them again. Tests follow the same ownership:

| Owner | Responsibility | Verification |
|---|---|---|
| Model | Local structural validity of detached values | Fixed validation and rejection cases |
| Macros and module boundaries | Generate declared capabilities and prevent unauthorized access | Real-module compile-pass and compile-fail tests |
| Runtime | Apply mappings and effects, enforce commit rules, manage COW and cache state | Focused mechanism regressions, including sharing and failure paths |
| Operation declaration and integration | Specify the correct access, effects and semantic obligations; connect the intended algorithm | Source review of declarations and focused integration checks |
| Domain algorithm | Produce the correct values and mappings for the operation's semantics | Fixed behavior regressions and public-API corpus parity |

The runtime can apply a deletion mapping correctly while an algorithm supplies
a structurally valid mapping that deletes the wrong atom. The former is a
runtime mechanism defect; the latter is an algorithm or declaration defect.
Neither layer's tests replace the other's.

Shared guarantees are tested at their owning boundary. An operation adds a
specific integration regression when its wiring or behavior needs one; it
does not copy a universal cache/COW/mapping/roundtrip suite. Registry metadata
alone is not proof that a declaration matches the reference semantics.

Compile-time restrictions and runtime checks are distinct defenses. Strict
checks verify declared invariants without becoming a second chemistry path.
Tests must respect the actual value-style versus in-place failure guarantees
and each operation's cache contract, not impose universal rollback or cache
recomputation promises. Unsafe/FFI owners likewise document their safety
assumptions and test the relevant boundaries; chemistry tests do not prove
memory safety.

Concrete test selection and comparison methods belong in
[test boundaries](./test_boundaries.md), not in another architecture checklist.

## 12. Parity verifies semantics; it does not invent them

A tempting porting loop is to observe a differential mismatch, add a special
case, and repeat until the corpus passes. That can produce an implementation
which fits known outputs but has no coherent explanation for an unseen input.

The source-backed loop starts elsewhere: determine the pinned source's control
flow and state transitions, reproduce them, then test the reproduction.
When parity fails, locate the first semantic divergence rather than asking
which local patch would make the final value match.

The defect can be in the algorithm, the operation contract, the runtime, or
the reference setup. Treating every mismatch as an algorithm bug can conceal
an incorrect contract; treating an invariant failure as inconvenient can
conceal an incorrect port.

Source anchors and separate behavioral/complexity review make a different
kind of audit possible. A reviewer can ask not only whether outputs agree,
but why a branch exists, which upstream behavior it represents, and whether
the Rust translation introduced an avoidable algorithmic cost.

Pinned versions make that reasoning reproducible. Approved project differences
remain explicit, and undefined upstream behavior does not acquire a meaningful
value merely because one reference run produced one. Neither case justifies
inventing a silent fallback.

Validation is evidence about an implementation, not a replacement for its
semantic basis.

## 13. Why the public Molecule belongs at the top

So far, the architecture has protected implementation boundaries. It must also
avoid making users pay for those boundaries in every call.

Rust inherent methods belong to the crate that defines their type. If the
public `Molecule` lived in the lowest-level model crate, higher-level search,
fingerprint, and conformer crates could not independently add inherent methods
to it. The usual alternatives are free functions, extension traits, or another
umbrella wrapper.

Those alternatives can be appropriate, but they make users more aware of the
internal decomposition. They may need to know which trait enables a method
or which crate owns a particular representation.

Owning the live `Molecule` in `cosmolkit` lets that crate provide a coherent
method-based API while calling each algorithm's unique owner. In this document,
“facade” means that public aggregation role; it does not mean an empty crate.
The same crate also owns the private runtime and lifecycle authority.

The important distinction from the large `MoleculeData` wrapper is what
crosses the boundary. Public methods do not hand the entire internal object
to every algorithm. They expose convenient user operations while the runtime
projects narrow detached inputs.

A simple public interface and a strongly constrained implementation are thus
two sides of the same arrangement, not opposing architectural choices.

## 14. One canonical value can support multiple behavior views

Not every useful type should follow `Molecule`'s ownership pattern.
`QueryGraph` is shared query data: predicates, graph structure, properties,
coordinates, and locally validated values. Moving parsing, matching, writing,
and fingerprint algorithms into its model crate merely to obtain method
syntax would reverse the intended dependency direction.

Creating a separate query model for each domain is no better. That would
introduce conversions and competing notions of query identity.

A borrowing behavior facade, such as `QueryGraphOperator`, offers another
arrangement: keep one canonical query value and interpret it through a
domain-owned view. The view does not become a second graph, and the model does
not need to depend on every algorithm that can consume it.

Construction stays distinct from interpretation. A parser produces a new owned
query; it does not need a pre-existing query for an operator to borrow:

```text
text -> parser -> owned QueryGraph -> borrowed behavior view
```

This distinction is more than naming preference. It keeps construction,
representation, interpretation, and runtime authority from being combined
solely to make a call look like a method.

Likewise, preserving query semantics does not require turning every concrete
molecule into a container for arbitrary query state. A file-result container
can distinguish concrete and query payloads while each model retains its own
meaning. Convenience belongs at the boundary, not in a duplicated or ambiguous
chemical representation.

## 15. Internal dependencies are not public feature promises

The same separation applies to features. A fingerprint implementation may use
a lower-level chemistry primitive without implying that users enabling that
fingerprint also requested every public operation associated with the
primitive's crate.

Cargo implementation dependencies answer what code is needed to compute a
result. Public feature bundles answer what capabilities the library exposes.
Conflating them makes an internal refactor unexpectedly change the user's API.

Explicit dependency wiring enables reuse; explicit public bundles define
exposure. The goal is not to eliminate either graph but to prevent the
implementation graph from silently becoming the public contract.

This completes a recurring theme: internal decomposition should support the
user's workflow rather than become another concept the user must reconstruct.

## 16. Why the pieces form a chain rather than a collection

Each design choice closes a gap left by another:

- Source fidelity does not prevent unauthorized state changes, so source ports
  also need operation contracts and authority boundaries.
- A written contract does not enforce itself, so declarations drive generated
  capabilities, narrow inputs, and private runtime access.
- Legal inputs do not guarantee a valid output, so candidates need mapping,
  effect, and invariant checks.
- Valid structure does not guarantee correct chemistry, so the result still
  needs source-backed parity.
- Detached computation can be expensive if implemented as whole-object copying,
  so isolation needs block-level sharing and COW.
- COW does not control algorithm authority, so it stays inside the restricted
  runtime lifecycle rather than becoming a mutable-object escape hatch.
- A safe internal design can still produce a fragmented public API, so the
  top-level public owner aggregates operations.
- A convenient public wrapper alone does not constrain its implementation, so
  its delegation crosses the same narrow authority boundary.

The resulting feedback loop is:

```text
pinned source behavior
        |
operation contract
        |
generated capabilities and detached inputs
        |
source-backed algorithm
        |
candidate state
        |
applicable contract checks and controlled installation
        |
reference parity validation
        |
first-divergence analysis -> source / contract / runtime / implementation
```

Parity and strict testing are development-time evidence in this loop, not
external oracle calls made during each user's operation.

The combination reduces both the opportunity to create certain errors and
the chance that remaining errors go undetected. These are distinct objectives.
A high test pass rate alone says little about how much unrelated authority
each implementation possesses, just as strong type boundaries alone say
little about whether the chemistry is right.

## 17. Architecture as external memory for agents

The central agent-era benefit is that more project knowledge can survive
outside an agent's current context.

An instruction to “leave coordinates alone” can be forgotten; an unavailable
coordinate capability remains unavailable. A note to preserve the required
mapping can disappear during compaction; a declared mapping obligation and
its validation still exist. A reminder not to construct live molecules in
a domain crate is stronger when that crate has no dependency or interface
providing the authority.

This does not remove the need to read standards, understand source code,
or design correct contracts. It relocates repeatable global obligations into
structures that do not depend on every implementer remembering them each time.

The desired scaling property is that project complexity can grow faster than
the amount of global knowledge required for one well-scoped operation. An
implementer should increasingly be able to work from a bounded set of facts:
the source behavior, the contract, the inputs available, the complete outputs
required, and the validation that establishes completion.

The goal is not less rigorous implementation. It is less need for each local
implementation to reconstruct the entire system before doing rigorous work.

## 18. The tradeoff: more framework work, less distributed complexity

This design is not necessarily the fastest way to produce the first working
feature. Contracts, generators, runtime lifecycles, typed mappings, and
validation machinery require substantial design and maintenance effort.

The intended return is a higher sustainable complexity ceiling. Instead of
spreading state-management obligations across every algorithm and binding,
the project concentrates them in a smaller set of shared mechanisms. Those
mechanisms can receive deeper review and repeated validation.

This is a form of complexity inversion:

```text
framework:       shared lifecycle and verification complexity
implementation: bounded source behavior and explicit inputs/outputs
user:           coherent chemistry operations
```

The framework is successful only if that separation is real. If every new
operation still requires a bespoke capability wrapper, another parallel
registry, or detailed knowledge of all runtime internals, the abstractions
have not delivered the promised reduction in local complexity.

The benefit should be judged by complete, auditable functionality and stable
public semantics, not by crate count, generated-code volume, or agent
throughput in isolation.

## 19. Constraints concentrate risk; they do not eliminate it

Shared machinery becomes a small but important trusted base.

**Contracts can be wrong.** A declaration that permits only one block when
the reference behavior requires two can enforce the wrong specification very
consistently. Source review must examine the declaration as well as the body.

**Generators have broad impact.** One wrapper or capability-generation defect
can affect many operations. Deterministic generation, focused structural tests,
and real-module privacy tests matter precisely because so much code relies
on these mechanisms.

**COW and failure semantics need exact reasoning.** Sharing, uniqueness,
derived-state ownership, returned errors, and unwinding can interact subtly.
Storage reuse cannot become a hidden aliasing or authority bypass, and the
documented distinction between value and in-place guarantees must remain clear.

**Strict and production builds must execute the same chemistry.** Additional
verification can expose errors; it must not make correctness depend on a
diagnostic-only computation or fallback.

**The model layer can become a dumping ground.** A shared representation belongs
there because of its role and dependencies, not merely because it is a struct.
Putting every convenient type into model recreates coupling at a lower level.

These risks explain why framework code deserves disproportionate scrutiny.
Constraining agents is not a substitute for engineering judgment. It focuses
that judgment on the contracts and mechanisms whose correctness many local
implementations rely on.

## 20. The project vision

COW, contracts, facades, intermediate representations, dependency DAGs, and
parity tests are not individually new. The value lies in connecting them so
that each addresses a limitation the others leave unresolved.

COSMolKit treats agent-assisted source reproduction as an engineering process
that needs enforceable boundaries, not merely as a faster way to generate code.
It seeks source-level behavioral fidelity without sacrificing a Rust-native
state model, production performance, or a usable API.

The long-term objective is a toolkit that can gain substantial new capability
without granting each contributor broader authority over everything already
built. The framework carries shared lifecycle obligations; domain owners carry
complete source-backed behavior; users receive a coherent library rather than
the internal machinery needed to produce it.

Higher agent throughput makes this separation more important, not less.
Producing more code faster increases the need for boundaries that prevent
local decisions from silently becoming global policy.

Constrained agent-assisted source reproduction therefore has two inseparable
goals: reduce the space in which invalid implementations can arise, and make
the remaining semantic errors traceable. The desired outcome is not simply
more automation. It is reliable growth in functionality while preserving
behavioral clarity, performance, and user-facing simplicity.

## Authoritative references

This essay explains why the design fits together. Exact requirements remain in
[crate architecture](./crate_architecture.md),
[operation contracts](./operation_system_standard.md),
[public API design](./public_api_design.md),
[source reproduction](./source_reproduction_protocol.md), and
[test boundaries](./test_boundaries.md).
The [sole split-crate plan](./plans/crate_architecture_completion_plan.md)
records execution and acceptance; this rationale grants no implementation
permission and makes no completion claim.
