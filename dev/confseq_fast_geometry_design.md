# COSMolKit ConfSeq FastGeometry Plan

This document defines the target FastGeometry design. It should describe only
the intended architecture and unfinished work. Historical implementations,
temporary benchmark snapshots, and abandoned intermediate routes do not belong
here.

## Decoder Boundary

COSMolKit keeps two ConfSeq template backends with different responsibilities.

`ConfSeqTemplateBackend::DistanceGeometry` is the broad reference backend. It
uses the RDKit-derived distance-geometry implementation and may use UFF
relaxation according to decode options.

`ConfSeqTemplateBackend::FastGeometry` is the fast backend. It must not
run whole-molecule UFF/MMFF optimization and must not silently fall back to the
distance-geometry backend. Its allowed conditioning target is the minimum
distance-geometry work already present inside RDKit `EmbedMolecule()`, applied
only to rigid fragments so the work is smaller and fragment templates can be
cached.

## Target Pipeline

FastGeometry is a rigid-fragment template builder plus a rigid assembly step:

```text
parse molecule
-> identify ConfSeq angle-token centers
-> cut ConfSeq rotatable single bonds
-> extract rigid heavy-atom fragments
-> build/cache fragment templates
-> condition each fragment with EmbedMolecule-like local distance constraints
-> assemble fragments by making reciprocal connector-stub bond segments coincide
-> validate local geometry and stereo
-> run the same full ConfSeq angle, dihedral, and ring-deferred decode used by
   the DistanceGeometry backend
```

FastGeometry replaces only the initial coordinate source. It must not redefine
ConfSeq recovery semantics. After the base coordinates exist, the decoder must
apply the same token stream, same angle transforms, same acyclic dihedral
transforms, and same ring-deferred dihedral logic as
`ConfSeqTemplateBackend::DistanceGeometry`.

Fragment cut bonds are an initialization and assembly concept only. They may be
used to construct and attach cheap rigid fragments, but they must not filter or
reinterpret final ConfSeq token application.

Implementation work is split into exactly two framework categories:

1. Fragment construction: identify rigid fragments, build reusable templates,
   and realize fragment-local coordinates.
2. Fragment assembly: attach realized fragments through cut-bond connector
   stubs by shared-bond rigid transforms only.

Issues should be classified into one of these categories before code changes.
Functional-group-specific fixes are out of scope for FastGeometry unless they
are introduced as source-backed generic parameter data used by the fragment
template model.

## Fragment Template Model

Each rigid fragment should be represented before coordinate realization:

```text
RigidFragmentTemplate
  atoms
  realization_atoms: owned atoms plus external connector stubs
  bonds:          atom-atom target lengths
  angles:         atom-center-atom target angles not delegated to ConfSeq tokens
  connectors:     cut rotatable bonds represented as external stub atoms
  shape_class:    Atom | Bond | Angle | CenterPlanar | CenterTetrahedral |
                  PlanarPi | RingPlanar | RingPuckered |
                  FragmentInternalRingSystem | Unsupported
  cache_key:      canonical local topology and atom-state key
```

Template data sources are local and source-backed:

- Bond lengths: UFF `r0` from the ported UFF parameter API, with explicit
  source-anchored static defaults only where the source model defines them.
- Non-tokenized bond angles: UFF `theta0` interpreted as degrees and converted
  to radians, RDKit ring-angle rules, and other already ported RDKit
  bounds-builder priors.
- Planar/pi shape: RDKit 2D scaffold shape may provide an initial layout or
  planarity prior, but average scaffold scaling is not a final geometry
  solution.
- 1-4 local preferences: RDKit cis/trans/other priors, used only to select a
  finite local template when needed.
- Stereo validation: RDKit-derived double-bond and chiral-volume checks.

## Fragment Conditioning Rules

Fragment conditioning is the target local-coordinate path. Connector stubs
participate in fragment conditioning exactly like local atoms, but final
coordinate commit writes only fragment-owned atoms.

The conditioning boundary is:

- build a fragment-local bounds matrix from the same RDKit-derived 1-2, 1-3,
  1-4, ring, planar, stereo, and local prior data used by `EmbedMolecule()`;
- triangle-smooth the fragment-local bounds matrix;
- generate a fragment-local coordinate realization using the already ported
  distance-geometry embedding machinery;
- run only the distance-geometry constraint minimization that `EmbedMolecule()`
  uses to satisfy bounds and local stereochemistry;
- cache the resulting fragment template by canonical local topology and atom
  state.

This is not `UFFOptimizeMolecule()`. It must not add full-molecule UFF/MMFF,
nonbonded packing, or a later relaxation step after fragment assembly.

Closed-form realizers remain allowed only for trivial fragments where they are
exact by construction:

- `Atom`: one atom at the local origin.
- `Bond`: two atoms on the local x-axis at target bond length.
- `Angle`: three atoms constructed directly from two bond lengths and one angle.
- single-center star fragments where all local bond lengths and non-tokenized
  angles are satisfied exactly.

Nontrivial planar/pi, ring, polycyclic, and multi-center fragments should use
the fragment-local conditioner. Average bond-length scaling of scaffold
coordinates is not a target mechanism.

Every realized template must pass cheap local validation before assembly:

- every template bond, including connector-stub bonds, must be near its target
  length;
- every recorded local template angle must be finite and near its target angle;
- every connector core and external stub atom must be present in realized local
  coordinates.

Unfinished realizers:

- generic fragment-internal templates for non-simple nonplanar ring systems.
  These must produce one realized rigid fragment with connector stubs; they
  must not introduce a separate ring scheduler in the assembly layer.
- more accurate nonplanar ring template families beyond the first simple-ring
  finite template, especially larger saturated rings where the current
  alternating/envelope pucker is locally valid but not close enough to the
  distance-geometry reference in rigid RMSD.
- multi-center rigid fragments that are not decomposable into a single center,
  planar/pi scaffold, or finite ring template.
- batch-global fragment-template cache and cache-hit diagnostics.

These unfinished cases must return structured unsupported errors until they
have explicit templates or the fragment-local conditioner supports them. They
must not fall back to whole-molecule distance geometry or force-field
relaxation.

## Force-Field Boundary

UFFOptimizeMolecule is not part of FastGeometry.

Allowed:

- read UFF bond equilibrium lengths;
- read UFF angle equilibrium values;
- use those values when building fragment templates.
- use the already ported RDKit distance-geometry force-field terms that
  `EmbedMolecule()` uses internally to satisfy bounds, but only inside one
  rigid fragment template.

Not allowed:

- whole-molecule UFF or MMFF optimization;
- force-field relaxation after fragment attachment;
- nonbonded packing optimization;
- using force-field energy to repair global RMSD;
- hidden fallback from unsupported fragments to force-field fitting.
- any fragment conditioner more expensive or more global than the corresponding
  `EmbedMolecule()` work on that fragment.

The intent of cutting into fragments is speed and reuse. The fragment
conditioner should do less work than full-molecule `EmbedMolecule()` by using
smaller bounds matrices and cached results, not more work.

## Structural Diagnostics And Operability

Topology-only structural classes are diagnostics, not routing decisions. They
may describe rings, shared ring centers, fused edges, or nonplanar multi-ring
components, but they must not route a molecule away from FastGeometry by
themselves.

FastGeometry unsupported status is based on initializer/decode operability:

- collapsed or near-zero bond lengths in the initialized coordinates;
- unusable angle triplets or dihedral frames needed by ConfSeq transforms;
- ring-deferred decode cannot apply required ring single-bond dihedral tokens
  after temporary ring-bond removal;
- the cheap initializer cannot produce connected, finite, nondegenerate
  coordinates for the molecule.

Not allowed:

- topology-only rejection for nonplanar rings, fused rings, shared ring atoms,
  or pucker classes before attempting cheap initialization;
- rejecting merely because a molecule contains a nonplanar five- or
  six-member ring; ConfSeq has ring-deferred dihedral recovery for ring
  single-bond tokens;
- rejecting from global RMSD, rigid-fragment RMSD, or any reference/DG
  coordinate comparison;
- functional-group-specific rejection rules;
- expanding this into a force-field repair path.

The production-facing topology result is `ConfSeqBaseStructuralRiskPrecheck`.
For now, its `classes` field is diagnostic and its topology-derived
`fallback_candidate` remains false for API compatibility with earlier
diagnostics. New code should treat it as a legacy diagnostic field; production
FastGeometry routing must come from explicit initializer/decode degeneracy
checks, not from topology class alone.

## Assembly Rules

Fragment interiors are rigid after realization. Assembly may translate and
rotate whole fragments, but must not deform internal bond lengths or
non-tokenized local angles.

The preferred assembly unit is the ConfSeq rotatable cut bond, not a single
atom anchor. Each side of a cut bond keeps an external connector stub in its
local template:

```text
left template:  ... i - j - k(stub)
right template: j(stub) - k - l ...
```

When one side has been placed, it records the global position predicted for the
other side of the cut bond. The unplaced side then aligns its local
`j(stub)-k` connector segment to the placed global `j-k` connector segment.
Rotation around `j-k` remains arbitrary at this stage because ConfSeq dihedral
tokens own that degree of freedom.

Assembly is defined by reciprocal connector stubs:

- Connector-stub attachment is the default for ConfSeq rotatable cuts when both
  sides expose the shared connector bond. It aligns the shared connector bond
  and must not write stub atoms into the final molecule coordinates.
- Rotation around the shared connector bond is not optimized by FastGeometry.
  ConfSeq dihedral tokens own that degree of freedom.
- Single-anchor and multi-anchor placement are temporary diagnostic incomplete
  paths for missing connector information, not the target path.
- Fragment interiors stay rigid after realization.

Assembly should progress in this order:

1. reciprocal connector-stub shared-bond alignment for ConfSeq rotatable cuts;
2. structured unsupported error when the required reciprocal stubs are absent;
3. diagnostic incomplete paths only while implementation is incomplete,
   reported separately from the target shared-bond path.

## Token Application Boundary

FastGeometry must keep ConfSeq token parsing compatible with the reference
DistanceGeometry backend. The token stream is not shortened and token counts
are not reinterpreted by fragment construction.

Execution is shared, not narrowed. FastGeometry and DistanceGeometry must use
the same final decode routine after their respective initial templates are
built:

```text
DistanceGeometry base -> full ConfSeq decode -> final
FastGeometry template -> full ConfSeq decode -> final
```

In particular, ring-internal single-bond dihedral tokens must still reach the
ring-deferred decode path. FastGeometry must not disable them merely because a
ring was represented as one rigid fragment during cheap initialization.

The intended division of labor is:

```text
FastGeometry: cheap nondegenerate initial coordinates
ConfSeq tokens: all encoded recovery semantics
```

## Metrics

Development metrics must compare FastGeometry against
`ConfSeqTemplateBackend::DistanceGeometry`.

Required report fields:

- distance-geometry reference success count;
- FastGeometry success where reference succeeds;
- final decoded global heavy-atom RMSD;
- automorphism-aware final decoded global heavy-atom RMSD;
- initializer/decode unsupported count and rate, reported separately from pass
  rates;
- rigid-fragment RMSD after cutting ConfSeq rotatable bonds as a diagnostic;
- rigid-fragment thresholds at 0.1 A, 0.2 A, and 0.3 A;
- rigid-fragment RMSD grouped by framework key, including nonplanar ring
  topology;
- worst rigid-fragment diagnostics with atom and bond details;
- failure classes by structured error.
- cut-bond count, reciprocal connector-stub count, shared-bond assembly count,
  and temporary incomplete-path assembly count.
- mirror-branch-like rigid-fragment counts: proper rigid RMSD fails while
  mirror-invariant shape RMSD passes. This keeps physical local-shape errors
  separate from distance-geometry branch choices that are not fixed by tuning
  local bond lengths or angles.

Final decoded RMSD against the DistanceGeometry decoded reference is the
primary metric. Pre-token geometry, rigid-fragment RMSD, shape RMSD, connector
RMSD, and non-token local constraint error are diagnostics only. They identify
whether remaining failures come from cheap initialization, fragment assembly,
or ConfSeq token recovery, but they must not become the optimization target.

Proper rigid-fragment RMSD remains required and must not be hidden. When
proper RMSD is high but mirror-invariant shape RMSD is low, the issue should be
classified as fragment branch consistency, not fragment geometry construction.
The acceptable solutions are framework-level:

- preserve enough source-backed stereo/context atoms during fragment
  conditioning to make the local branch deterministic when the full molecule
  contains such information;
- use graph-automorphism or self-match logic only to diagnose symmetry-equivalent
  labels, not to arbitrarily permute final atom identities;
- leave genuinely unconstrained mirror-equivalent nonchiral fragments as
  physically valid local shapes unless a cut-boundary or stereo constraint can
  source-back the branch choice.

Unacceptable solutions are local angle tuning, functional-group-specific
patches, or whole-molecule force-field repair.

## Implementation Progress

Completed target pieces:

- ConfSeq rotatable-bond based rigid-fragment extraction.
- Rigid-fragment classification.
- Explicit fragment template data structure separate from realized coordinates.
- Per-build fragment template cache keyed by local atom and bond state.
- Rotatable cut-bond connector stubs stored separately from fragment-owned
  atoms.
- UFF-backed bond length targets.
- UFF/RDKit-backed local angle targets, including explicit UFF degree-to-radian
  conversion.
- Closed-form realizers for atom, bond, three-atom angle, single-center
  two/three/four-ligand fragments, and acyclic simple-path/tree fragments where
  target bond lengths are exact.
- Non-simple nonplanar ring topology classification:
  `edge_fused_chain`, `edge_fused_polycyclic`, `spiro`, `bridged_or_cage`,
  `single_non_simple`, and `unknown`.
- Non-simple nonplanar ring dispatch now has an explicit fragment-internal
  support boundary. Unsupported ring-system fragments return structured errors;
  assembly does not special-case their topology.
- Connector stubs are included in `realization_atoms`; post-hoc guessed stub
  directions have been removed from the target path.
- Realized template bond, angle, and connector-stub validation.
- Shared connector-bond alignment for attaching fragments across ConfSeq
  rotatable cuts is the default assembly path when a connector target is
  available and the unplaced fragment has the corresponding local stub segment.
- FastGeometry and DistanceGeometry share full ConfSeq token execution after
  their initial templates are built; fragment cut bonds do not filter decode.
- Structured errors for unsupported FastGeometry cases.
- Corpus diagnostics against the DistanceGeometry backend, including rigid
  fragment RMSD, shape RMSD, terminal-symmetry RMSD, connector-inclusive RMSD,
  nonplanar ring topology counts, and framework-key RMSD buckets.
- Separate stereo context atoms for fragment conditioning. These atoms may
  constrain local branch choice during fragment-local embedding, but they do
  not become fragment-owned realization atoms and do not pollute template bond,
  angle, cache, or assembly state.

Unfinished target pieces:

- Batch-global fragment-template cache and cache-hit reporting.
- Keep the default initializer cheap and nondegenerate. Fragment-local
  EmbedMolecule-like conditioning is the target for nontrivial fragments where
  closed-form placement is not exact; it is not a route for making fragment
  shapes match DistanceGeometry by global repair.
- Enforce reciprocal connector-stub coverage as the normal assembly contract,
  then retire single-anchor and multi-anchor temporary diagnostic paths from
  the target path.
- Branch-aware fragment conditioning and assembly diagnostics. This is not a
  blind mirror-candidate search; branch changes are allowed only when justified
  by stereo context, connector-stub constraints, or graph self-match evidence.
- Generic fragment-internal realization for non-simple ring systems when they
  appear inside a single ConfSeq rigid fragment.

The next implementation step is to separate final decoded failures into cheap
initializer degeneracy, fragment assembly failure, and ConfSeq token recovery
failure. Pre-token geometry may be measured for diagnosis, but it is not a
success criterion.

Current diagnostic boundary:

- Fragment construction is monitored by shape-class counts, nonplanar ring
  topology counts, structured template failures, and framework-key rigid RMSD
  buckets.
- Fragment assembly is monitored by cut-bond count, reciprocal connector-stub
  count, shared-bond assembly count, and explicit incomplete-path count.
- Token recovery is monitored by comparing pre-token FastGeometry geometry
  with final decoded geometry and DistanceGeometry reference RMSD.
- Fused/spiro/bridged labels may explain fragment-internal realization
  failures, but they must not appear as assembly-layer scheduling concepts.
- Mirror-branch-like failures are reported separately from shape failures.
  They are a framework-level branch-consistency problem, not evidence that the
  fragment should be locally angle-optimized or expanded with arbitrary
  functional-group rules.

## Usage Examples

Rust callers should select the template backend explicitly when comparing the
two paths. `DistanceGeometry` is the reference path; `FastGeometry` is the
fast path that only replaces the initial coordinate source and then runs
the same ConfSeq token decode.

```rust
use cosmolkit::{
    ConfSeqDecodeOptions, ConfSeqTemplateBackend, decode_confseq_record_with_options,
};

// Tokenized corpus record: dude_aa2ar__L021087 candidate 0.
// Spaces are ConfSeq token boundaries and must be preserved.
let confseq = "N C <112> | ( = O ) <84> C <111> | <-45> n 1 c c ( <21> N <123> | <173> C <116> | ( = O ) <120> c 2 c c c c c 2 <112> C <112> | <-66> N 2 <-172> C ( = O ) <-2> C <0> N <3> C <174> 2 = O ) c n 1";

let reference = decode_confseq_record_with_options(
    confseq,
    &ConfSeqDecodeOptions {
        optimize_with_uff: false,
        template_backend: ConfSeqTemplateBackend::DistanceGeometry,
        ..ConfSeqDecodeOptions::default()
    },
)?;

let fast = decode_confseq_record_with_options(
    confseq,
    &ConfSeqDecodeOptions {
        optimize_with_uff: false,
        template_backend: ConfSeqTemplateBackend::FastGeometry,
        ..ConfSeqDecodeOptions::default()
    },
)?;
```

Python callers use the same explicit backend names:

```python
import cosmolkit

# Tokenized corpus record: dude_aa2ar__L021087 candidate 0.
# Spaces are ConfSeq token boundaries and must be preserved.
confseq = (
    "N C <112> | ( = O ) <84> C <111> | <-45> n 1 c c ( <21> N <123> | "
    "<173> C <116> | ( = O ) <120> c 2 c c c c c 2 <112> C <112> | <-66> "
    "N 2 <-172> C ( = O ) <-2> C <0> N <3> C <174> 2 = O ) c n 1"
)

reference = cosmolkit.confseq.decode(
    confseq,
    optimize_with_uff=False,
    template_backend="distance_geometry",
)

fast = cosmolkit.confseq.decode(
    confseq,
    optimize_with_uff=False,
    template_backend="fast_geometry",
)
```

## Validation

After FastGeometry code changes, run:

```bash
cargo fmt --all
cargo check -p cosmolkit-core --features op-contracts-strict
cargo test -p cosmolkit-core confseq_base_corpus_pass_rate_snapshot --features op-contracts-strict -- --ignored --nocapture
```

For broader core changes, also run the full strict core test suite required by
the repository guidelines.
