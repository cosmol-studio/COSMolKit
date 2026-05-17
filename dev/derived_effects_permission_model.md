# Derived Effects Permission Model

## Core Idea

Each category in `derived_effects` is no longer just a "side-effect label". It defines a set of **read/write permissions** for that derived state.

## Permission Table

```
Category     | Read | Write | Meaning
requires     |  ✅  |  ✅  | Prerequisite: the operation may read this cache, and may rebuild it if needed
recompute    |  ❌  |  ✅  | Old value is stale: write-only, no read access. Must compute from raw data
preserve     |  ✅  |  ❌  | Old value is still valid: read-only. Can record a preservation proof, but cannot mutate
invalidate   |  ❌  |  ❌  | Neither read nor write. Simply declares "this cache is now dirty"
unsupported  |  ❌  |  ❌  | The system knows this cache exists, but the operation cannot handle it
```

### Per-category Rationale

**requires**: The operation body needs to read this cache to perform its chemistry logic. For example,
`assigned_valence` needs `adjacency` to iterate each atom's incident bonds. The body may also
rebuild the cache if it finds it is `None`, so write is also allowed.

**recompute**: The operation commits to computing this cache from raw data (e.g., the bond table).
The old cache value is untrustworthy, so reading it is forbidden. After the operation completes,
this cache is guaranteed fresh.

**preserve**: The operation can prove that topological changes did not affect this cache.
Reading the old value is allowed. Writing the cache directly is forbidden —
only `prove_preserved()` may be called.

**invalidate**: The operation marks this cache as dirty. It cannot read or write it.
Downstream operations will detect the staleness via `needs_update()` and handle it.

**unsupported**: The system acknowledges this cache exists, but the current operation
cannot handle it correctly. No read or write access — explicitly recorded as "known gap."

## Dev-mode Verification

Under `op-contracts`, every cache getter/setter in `OpParts` checks the trace:

- Reading `adjacency` → verify that `requires` or `preserve` contains `adjacency`; reject if undeclared
- Writing `valence` → verify that `requires` or `recompute` contains `valence`; reject if undeclared
- At `finish()`: confirm every cache in `recompute` was written; confirm every cache in `requires` is `Some`

## Non-op Entry Points

Some chemistry functions (e.g., `atom_has_valence_violation`) are called directly from `&Molecule`,
bypassing the operation system entirely. In this case they manage their own cache fallback
("use if available, build if not"). This does not violate the model — no ops path means no permission check.

Operation bodies that go through ops **must** obey the permission model.
Direct API calls that bypass ops are the caller's own responsibility.

## Per-operation Analysis

### 1. `with_hydrogens` / `without_hydrogens` / `without_hydrogens_with_params`

These are strong topology operations (append/compact atoms and bonds).

**Actual cache access:**
- No direct cache reads or `set_*_cache` writes
- `clear_cache`: adjacency, valence, aromaticity, stereo, drawing, fingerprint (and rings/ring_families for without_hydrogens)
- `prove_preserved` (with_hydrogens only): rings, ring_families
- `mark_states_handled`: valence, aromaticity, stereo

**New model mapping:**
```
requires:        []            (no cache reads)
recompute:       []            (no cache writes)
preserve:        [rings, ring_families]   (with_hydrogens only, read-only + proof)
invalidate:      [adjacency, valence, aromaticity, stereo, drawing, fingerprint]
```

`preserve` fits perfectly — the body proves rings are untouched and never writes them.
`invalidate` fits perfectly — the body clears these and never touches them again.
`mark_states_handled` is not a cache permission — it's a **contract obligation** that must
remain a separate mechanism (already enforced by `validate_contract` via `must_handle`).

**Verdict**: ✅ Fits cleanly.

---

### 2. `with_kekulized_bonds`

**Actual cache access:**
- Reads `derived_cache().rings`                            → requires read
- `set_rings_cache(rings)` (conditional, when None)       → requires or recompute write
- `set_valence_cache(valence)` (always)                   → recompute write
- `clear_cache`: aromaticity, drawing, fingerprint        → invalidate

**New model mapping:**
```
requires:        [rings]
recompute:       [valence]
invalidate:      [aromaticity, drawing, fingerprint]
```

The conditional `set_rings_cache` is interesting — it only writes rings when the cache
was missing. Under `requires: [rings]`, this is allowed (requires permits both read and write).
The write-once-when-missing pattern is valid: the body reads rings, finds None, rebuilds them.

**Verdict**: ✅ Fits. `requires` handles the "read, and write only if absent" pattern.

---

### 3. `sanitized`

This is the most complex operation. It conditionally rebuilds many caches depending on which
`SanitizeOps` flags are set.

**Actual cache access:**
- Reads `derived_cache().rings` (multiple conditional paths)   → requires read
- `set_adjacency_cache` (via `sanitize_recompute_property_cache`, called up to 7 times)
- `set_valence_cache` (via same helper, up to 7 times)
- `set_rings_cache` (conditional, in SYMMRINGS, KEKULIZE, SET_AROMATICITY, CLEANUP_ATROPISOMERS)
- `clear_cache` for ring_families, aromaticity, stereo, drawing, fingerprint (initial clear)
- Additional conditional `clear_cache` per sanitize step

**New model mapping:**
```
requires:        [rings]
recompute:       [adjacency, valence]     (rebuilds these from raw data as needed)
preserve:        []
invalidate:      [ring_families, aromaticity, stereo, drawing, fingerprint]
```

The conditional `set_rings_cache` is declared in requires (because the body reads rings
to decide whether to rebuild them). Under `requires: [rings]`, both the reads and the
conditional write are legal.

But there's a subtlety: adjacency is rebuilt up to 7 times (each time `sanitize_recompute_property_cache` is called). Under `recompute: [adjacency]`, each `set_adjacency_cache` call is a valid write. The framework doesn't care how many times you write — it only checks that you don't read.

**Verdict**: ✅ Fits, though verbose. The "rebuild adjacency + valence as dependency" pattern
occurs in every sanitize sub-step. This is inherent to the sanitize pipeline, not a model flaw.

---

### 4. `assigned_valence`

**Current (problematic) declaration:**
```
recompute: [adjacency, valence]
```

**Actual cache access:**
- `set_adjacency_cache` (always)
- `set_valence_cache` (always)

The body does NOT read adjacency from the cache — it rebuilds it from topology:
```
let adjacency = AdjacencyList::try_from_topology(read.num_atoms(), read.bonds())?;
```

If the body were refactored to reuse cached adjacency:
```
requires:        [adjacency]      (needs it, may read if present)
recompute:       [valence]        (produces this)
```

If the body always rebuilds adjacency from scratch:
```
recompute:       [adjacency, valence]
```

Both are valid under the new model. The choice depends on whether the implementation
wants to check the cache or always rebuild. **The model enforces consistency**:
declare requires → you may read; declare recompute → you must not read the old value.

**Verdict**: ✅ Fits both implementation strategies. Eliminates the old semantic lie
where `recompute: [adjacency]` was declared but the value was silently reused.

---

### 5. `assigned_rings`

**Current:**
```
recompute: [adjacency, rings]
```

Same pattern as `assigned_valence` — adjacency is always rebuilt from topology,
then rings are computed from adjacency + topology.

**New model:**
```
requires:        [adjacency]      (needs it for ring perception; or keep in recompute if always rebuild)
recompute:       [rings]
```

**Verdict**: ✅ Same as assigned_valence.

---

### 6. `assigned_ring_families`

**Current:**
```
recompute: [adjacency, ring_families]
```

Same pattern — adjacency rebuilt, then ring_families computed.

**New model:**
```
requires:        [adjacency]
recompute:       [ring_families]
```

**Verdict**: ✅ Same.

---

### 7. `assigned_aromaticity`

**Current:**
```
recompute: [adjacency, rings, valence, aromaticity]
invalidate: [drawing, fingerprint]
```

This operation always rebuilds adjacency, rings, valence, and aromaticity from scratch.
It does not read any existing cache — everything is computed fresh from topology data.

**New model:**
```
requires:        []
recompute:       [adjacency, rings, valence, aromaticity]
invalidate:      [drawing, fingerprint]
```

This is a pure "recompute everything" operation. No requires needed.

**Verdict**: ✅ Fits cleanly. Demonstrates that an operation can have zero `requires`.

---

### 8. `assigned_radicals`

**Current:**
```
recompute: [adjacency, valence]
```

Same pattern as `assigned_valence`.

**New model:**
```
requires:        [adjacency]
recompute:       [valence]
```

**Verdict**: ✅ Same as assigned_valence.

---

### 9. `with_2d_coordinates`

**Current:**
```
invalidate: [drawing]
```

No derived cache reads or writes. Only clears drawing.

**New model:**
```
requires:        []
recompute:       []
invalidate:      [drawing]
```

**Verdict**: ✅ Simplest case.

---

## Cross-cutting Observations

### `require_handle` / `must_handle` is orthogonal
The current `require_handle: [valence, aromaticity, stereo]` on hydrogen operations
is a CONTRACT OBLIGATION, not a cache permission. The permission model does not replace it.
It should remain as a separate field.

### `preserve` semantics
Only `with_hydrogens` uses `preserve`. The new model's read-only constraint matches
the current behavior — the body calls `prove_preserved`, not `set_rings_cache`.

### adjacency dependency pattern
A common pattern across operations:
| Operation | Adjacency role |
|-----------|---------------|
| `assigned_valence` | needs it (requires or recompute) |
| `assigned_rings` | needs it (requires or recompute) |
| `assigned_aromaticity` | always rebuilds (recompute) |
| `sanitized` | always rebuilds (recompute) |
| `with_kekulized_bonds` | does NOT touch adjacency |

Under the old model, all of these except `with_kekulized_bonds` stuffed adjacency into
`recompute` — including `assigned_valence` which COULD reuse the existing cache but
was declared as "recompute." The new model forces a choice and enforces consistency.

### `clear_cache` validity
Clearing a cache is allowed under `invalidate` (or `recompute`).
It is NOT allowed under `preserve` (read-only) or `requires` alone (unless also in invalidate).
The model makes this constraint explicit.

