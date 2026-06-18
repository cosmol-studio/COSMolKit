# COSMolKit ConfSeq Decoder Status

COSMolKit currently has two native ConfSeq decoding paths. They serve different
purposes and should not be conflated.

## Distance-Geometry Path

`ConfSeqTemplateBackend::DistanceGeometry` is the closer-to-original path. It
parses ConfSeq tokens in Rust, prepares the molecule with COSMolKit, builds the
initial conformer through the RDKit-derived distance-geometry implementation,
and can optionally run the ported UFF relaxation.

This path is the default because it currently covers the broadest input space
and is the reference behavior for the native decoder.

## Base-Conformer Path

`ConfSeqTemplateBackend::BaseConformer` is a custom lightweight base conformer
for ConfSeq. It does not run distance geometry and does not silently fall back to
UFF or the distance-geometry path.

The intended contract is strict:

- Build a deterministic heavy-atom base geometry from topology.
- Preserve source-backed bond lengths and local geometry where implemented.
- Apply ConfSeq angles and dihedrals after the base geometry exists.
- Return structured unsupported errors for topology or stereochemistry that the
  base builder cannot model rigorously yet.

The existing Base-Conformer implementation should now be treated as a legacy
baseline for diagnostics and regression comparison, not as the architecture to
extend. Its current direct-coordinate construction has useful pieces, but it
mixes topology classification, local geometry, ring placement, stereo repair,
and fallback heuristics in one path. Further work should replace the base
generation path with a source-backed local-constraint model instead of adding
more case-specific placement fixes.

Current modeled coverage includes acyclic heavy-atom graphs, several local
functional-group geometries, planar aromatic rings, supported saturated rings,
edge-fused aromatic systems, explicit double-bond stereo planes, and explicit
tetrahedral stereo volume sign checks. Ring substituent directions are
constrained from the ring-neighbor local angle instead of a centroid-only radial
direction; aromatic and sp2 ring centers remain in-plane, while supported sp3
ring centers may use out-of-plane local geometry. Tetrahedral stereo correction
now prefers independent ligand-side moves and uses whole cyclic-side moves only
as a strict fallback after explicit chiral-volume revalidation.

## Current Base-Conformer Test Progress

The active corpus benchmark is:

```text
~/sh4090/confseq_test_strings_100x10.jsonl
```

The current comparison uses the completed distance-geometry/UFF output as the
dirty reference and measures heavy-atom RMSD for the base-conformer result with a
0.3 A pass threshold.

Latest recorded result:

```text
total candidates: 1000
dirty reference succeeded: 989
base succeeded where dirty succeeded: 829
RMSD pass count at 0.3 A: 415
total pass rate: 41.5%
pass rate among comparable candidates: 50.1%
p50 RMSD: 0.298 A
p75 RMSD: 0.683 A
p90 RMSD: 1.344 A
p95 RMSD: 2.443 A
p99 RMSD: 7.266 A
```

This is a net gain of 139 candidates over the initial 0.3 A baseline. The base
builder now has 160 strict base errors on this corpus, down from 194 before the
cyclic-side tetrahedral stereo fallback. The same work still exposes
non-monotonic tail behavior, so further base-conformer work should treat RMSD
regressions as diagnostic data, not noise.

The largest current strict gaps are tetrahedral stereo correction for ring-bound
centers that have no independent movable ligand side, and unsupported complex
ring systems with shared constraints that cannot be placed by the current
one-ring-at-a-time scaffold. These are base-conformer gaps, not ConfSeq parser
gaps.

## Base-Conformer Rebuild Direction

The next Base-Conformer should be a deterministic constructive solver driven by
RDKit-derived local prior knowledge, not a force-field-like coordinate
optimization loop.

Conceptually:

```text
RDKit-derived local constraints
  -> molecule component classification
  -> deterministic economical placement
  -> constraint validation
  -> ConfSeq angle/dihedral application
```

This differs from the distance-geometry backend:

```text
DistanceGeometry: constraints -> randomized embedding/minimization
BaseConformer:   constraints -> deterministic constructive placement
```

Both paths may share the same source-backed local chemistry knowledge, but the
Base-Conformer must not call DG embedding, UFF/MMFF relaxation, or a general
iterative optimizer. The goal is a broadly useful heavy-atom base conformation
for most small molecules, not exact strain modeling for every ring system.

### Rebuild Contract

- Use RDKit bounds-builder prior knowledge as the main source of local geometry.
- Extract exact target values or discrete preferences from the bounds logic
  where possible: bond lengths, 1-3 angle-derived distances, 1-4 cis/trans/other
  preferences, double-bond planarity, and chiral volume signs.
- Prefer closed-form construction, local frame placement, rigid transforms, and
  small finite template selection over iterative coordinate fitting.
- Keep all work deterministic for a fixed molecule and ConfSeq input.
- Keep computation close to linear or low-order local graph traversal for common
  small molecules.
- Return structured unsupported errors when a component requires a global
  constrained solve, rather than silently falling back to DG/UFF.
- Preserve diagnostics from the legacy path so coverage and failure classes can
  be compared during the rebuild.

### Non-Goals And Prohibited Shortcuts

- Do not tune coordinates by repeatedly minimizing a bounds or angle objective.
- Do not implement a hidden force field under the Base-Conformer name.
- Do not use DG embedding as an internal fallback.
- Do not chase corpus failures one by one by adding unprincipled geometry
  patches.
- Do not require exact reproduction of aromatic or fused-ring strain when the
  intended small-molecule base scaffold is chemically reasonable and the
  unsupported boundary is explicit.

### Source-Backed Constraint Model

Introduce a `ConfSeqBaseConstraintModel` before rewriting placement. It should
be the only source of local geometry decisions for the new Base-Conformer.

Suggested model:

```text
ConfSeqBaseConstraintModel
  bond_targets       # 1-2 bond-length targets
  angle_targets      # 1-3 target angles or equivalent target distances
  torsion_priors     # 1-4 cis/trans/other preferences
  planar_groups      # aromatic, sp2, double-bond, amide-like groups
  chiral_sets        # RDKit-style chiral volume constraints
  ring_components    # ring sharing graph and placement class
```

Primary RDKit-derived sources already ported in
`crates/cosmolkit-core/src/chemistry/distgeom.rs`:

- `set12Bounds`: source-backed 1-2 bond targets.
- `_set13BoundsHelper` and `compute13Dist`: angle-to-distance local geometry.
- `_setRingAngle`: ring angle rules by hybridization and ring size.
- `set14Bounds`, `_setInRing14Bounds`, and `_setChain14Bounds`: 1-4 local
  torsion preferences, ring cis/trans decisions, and amide-like handling.
- `set15Bounds`: short-range extended local constraints useful for validation,
  not for expensive fitting.
- `findChiralSets`, `_volumeTest`, and final chiral checks: chiral volume
  expression and validation.
- `doubleBondGeometryChecks` and `doubleBondStereoChecks`: double-bond
  planarity and stereo validation.

The model should extract central constructive intent from these routines rather
than copying the full DG process. For example, a 1-3 bound derived from a bond
angle can be converted back to an angle target for local frame placement, while
a 1-4 cis/trans preference can choose one of a finite set of deterministic
dihedral templates.

### Constructive Placement Plan

1. Build `ConfSeqBaseConstraintModel` from the molecule.
2. Classify connected components and ring components:
   - acyclic tree / branch components,
   - isolated aromatic or sp2 planar components,
   - simple saturated 5/6-member rings,
   - fused or spiro ring forests,
   - global constrained systems that remain unsupported.
3. Place acyclic heavy-atom regions with a Z-matrix-like local-frame builder:
   bond target + angle target + torsion prior determine each new atom.
4. Place planar aromatic/sp2 components as rigid local frames with source-backed
   bond and angle targets.
5. Place simple saturated rings from small finite templates selected by local
   targets and stereo priors, not by iterative minimization.
6. Attach ring and acyclic components by rigid transforms over shared atoms,
   shared bonds, or explicit anchor frames.
7. Enforce chiral sets during placement when possible; use post-placement
   correction only as a validation-backed fallback, not as the primary stereo
   strategy.
8. Validate local constraints cheaply:
   - bond target tolerance,
   - angle target tolerance,
   - selected torsion-prior class,
   - double-bond stereo,
   - chiral volume sign,
   - all atoms placed exactly once.
9. Apply ConfSeq angle and dihedral tokens using the existing decode stage.

### Development Steps

1. Keep the current implementation available only as a legacy baseline while
   the new path is built behind internal helpers.
2. Add the constraint model and tests for extracted constraints without changing
   coordinate generation.
3. Move existing bond-length, ring-angle, torsion-prior, and chiral-set logic to
   consume the constraint model.
4. Replace acyclic placement first, because it should benefit most from a
   Z-matrix-like deterministic builder.
5. Replace planar/aromatic ring placement next, using the same angle and bond
   targets rather than ad hoc polygon rules.
6. Replace saturated simple-ring placement with finite source-backed templates.
7. Reintroduce fused/spiro placement only after ring components expose shared
   atom/bond constraints explicitly.
8. Delete the legacy coordinate-generation path once the new path covers the
   same supported contract and diagnostics show no loss in failure
   transparency.

### Acceptance Criteria

- Default ConfSeq decoding remains `DistanceGeometry`.
- `BaseConformer` remains opt-in and strict.
- Unsupported systems return typed errors, not hidden fallback results.
- The new path is deterministic and avoids iterative minimization.
- Diagnostics identify whether failures come from constraint extraction,
  component classification, placement, validation, or ConfSeq token application.
- Benchmark comparisons use RMSD and failure-class movement as diagnostics, but
  implementation decisions must be justified by source-backed local chemistry
  constraints rather than by single-case overfitting.
