# Source-Reproduction Protocol for C/C++ Library Ports

This document defines the marking and review discipline for source-level
porting from any external C/C++ library (e.g., RDKit, Gemmi, OpenBabel) into
COSMolKit's Rust core. It is the single source of truth for the marker
convention used across the entire repository.

Every file that contains copied C++ source lines must include a reference to
this protocol in its header comment. The canonical file is:

```
dev/source_reproduction_protocol.md
```

If a specific module uses a library-specific variant (e.g.,
`dev/rdkit_porting_notes.md`), it must defer to this document for the marker
definitions and only add library-specific guidance.

---

## 1. Copy Rule

The original C++ source lines **must** be copied verbatim into the
corresponding Rust function body as comments. Do not summarise, paraphrase,
or move them to an unrelated reference appendix when a corresponding Rust
function exists.

**Rationale:** Copied lines serve as an in-line diff anchor. Reviewers and
agents must be able to compare the original and the Rust implementation
side-by-side without switching files.

## 2. Two-Axis Status Marker

Each copied C++ line carries a **two-symbol** status marker placed immediately
before or after the copied block:

```
// <Library>❌❌: <description>
// <Library>❗✔️: <description>
```

The format is: `<LibraryPrefix><firstSymbol><secondSymbol>:`

- The **first** symbol records **behavioral reproduction** status.
- The **second** symbol records **performance and algorithmic-complexity**
  status within the currently modeled input state space.

These two axes **must not be mixed**.

### First Symbol — Behavioral Reproduction Status

| Symbol | Meaning |
|--------|---------|
| `❌`   | The line has **not** been behaviorally reproduced. |
| `❗`   | The Rust code implements an approximate, partial, or structurally
         | similar behavior, but exact behavioral equivalence has not yet been
         | proven, or the behavior may still differ on edge cases. |
| `✔️`   | The Rust code below the copied C++ block reproduces the source
         | library behavior for the currently modeled input state space. |

### Second Symbol — Performance and Complexity Status

| Symbol | Meaning |
|--------|---------|
| `❌`   | The Rust implementation is known or expected to be materially
         | worse than the original in performance or algorithmic complexity.
         |
         | Use when local code inspection reveals an obvious cost, such as
         | repeated scanning where the original uses indexed access, extra
         | allocation or cloning in a hot path, a less suitable data structure,
         | worse asymptotic complexity, avoidable buffering overhead, or an
         | intentionally simplified early-port implementation. |
| `❗`   | Performance status remains **unresolved** after local performance
         | review.
         |
         | ⚠️ **Not a default** for different-looking Rust code. Must **not**
         | be used merely because the implementation has not been benchmarked.
         | Before using `❗`, the reviewer or agent **must** inspect the
         | implementation shape, including data structures, allocation
         | behavior, loop nesting, lookup complexity, parsing/buffering
         | strategy, ownership and cloning behavior, temporary object creation,
         | and hot-path branching. Only use when local inspection does not
         | reveal a clear material win or loss. |
| `✔️`   | Local code inspection supports performance- and
         | complexity-equivalence with the original for the currently modeled
         | input state space.
         |
         | Does not require identical source form. Requires comparable
         | asymptotic complexity, allocation behavior, data-structure cost,
         | and no obvious hot-path regression. |
| `🔝`   | The Rust implementation is expected or measured to be materially
         | **better** than the original in performance or algorithmic
         | complexity, while preserving the original behavior for the
         | currently modeled input state space.
         |
         | Only use when the improvement is explainable from the code or
         | measured by benchmarks (better asymptotic complexity, fewer
         | allocations, direct indexing, simpler buffering, more cache-friendly
         | representation). |

## 3. Valid Marker Combinations

| Marker              | Meaning |
|---------------------|---------|
| `Lib❌❌:`           | Not implemented. |
| `Lib❗❗:`           | Approximately implemented; behavior is not fully
                      | confirmed and performance status remains unresolved. |
| `Lib❗✔️:`           | Approximately implemented; behavior may still differ,
                      | but local inspection supports perf/complexity
                      | equivalence. |
| `Lib✔️❗:`           | Behavior is confirmed equivalent, but performance
                      | status remains unresolved after local review.
                      |
                      | Requires prior inspection of: data structures,
                      | allocation behavior, loop structure, lookup
                      | complexity, parsing/buffering strategy, ownership/
                      | cloning costs, temporary object creation, and
                      | hot-path branching.
                      |
                      | If a clear material win or loss is visible, use
                      | `Lib✔️🔝:` or `Lib✔️❌:` instead. |
| `Lib✔️✔️:`           | Behavior, performance, and complexity are considered
                      | equivalent after local code inspection. |
| `Lib✔️❌:`           | Behavior is confirmed equivalent, but the Rust
                      | implementation is known or expected to be slower or
                      | algorithmically worse. |
| `Lib✔️🔝:`           | Behavior is confirmed equivalent, and the Rust
                      | implementation is expected or measured to be faster
                      | or algorithmically better.
                      |
                      | **Must** have a nearby Rust comment explaining why
                      | the improvement does not change semantics and why it
                      | is expected or measured to be better. |

## 4. Behavioral Upgrade Rules

- A `✔️` first marker is only allowed after the Rust code **actually
  implements** the original behavior for the currently modeled input state
  space. Do not batch-upgrade behavior markers without checking the
  corresponding Rust behavior.

- A `✔️` second marker is only allowed after **local performance review**
  supports performance- and complexity-equivalence. Review must cover at
  least: asymptotic complexity, loop nesting, repeated scanning, allocation
  and temporary object creation, cloning and ownership movement,
  map/set/vector choice and lookup complexity, parsing/buffering strategy,
  cache locality, and hot-path branching.

## 5. Special Cases

- **`🔝` second marker**: Requires a nearby Rust comment explaining why the
  improvement preserves semantics and why it is better.

- **`❌` second marker**: Use when the Rust implementation preserves behavior
  but is known or expected to be materially worse in performance/algorithmic
  complexity (e.g., a simple implementation intentionally replacing an
  optimized original path during an early porting phase).

- **Implicit original behavior**: If C++ behavior is implicit, surprising, or
  depends on language rules (unsigned overflow, substring bounds, stream
  state, object lifetime, exception behavior), reproduce the observed
  behavior explicitly or keep the first marker as `❗` until the edge behavior
  is understood and tested.

- **Unsupported features**: If a copied C++ function depends on source
  library features that COSMolKit has not modeled yet, keep the relevant
  lines as `❌❌` until the required state, operation, or policy exists.
  Do not hide unsupported behavior behind a heuristic.

- **Cross-file helper functions**: Source library logic may be split across
  helper files. If a Rust function implements behavior from another source
  file, the corresponding helper function body from that file must be
  copied into the Rust function before behavior is expanded. For example,
  if top-level dispatch lives in `parser.cpp` but the real logic is in
  `parser_helpers.cpp`, future work must align against the helper bodies
  directly, not only against the dispatcher call sites.

## 6. Testing

Tests should be added at the smallest stable boundary where a reproduced
behavior can be observed. Boundary cases from C++ language behavior are part
of the port, not optional cleanup.

## 7. Plain-Text Fallback

If rendering in a plain-text environment without emoji support, use these
substitutes:

| Unicode | Fallback |
|---------|----------|
| ❌      | `X`      |
| ❗      | `!`      |
| ✔️      | `V`      |
| 🔝      | `UP`     |

---

## 8. Examples

### Example A — RDKit (current primary use)

```cpp
// RDKit source:
//   if (!params.doKekule && atom->getIsAromatic() && symb[0] >= 'A' &&
//       symb[0] <= 'Z') {
//     symb[0] = tolower(symb[0]);
//   }
```

```rust
// RDKit✔️✔️: if (!params.doKekule && atom->getIsAromatic() && symb[0] >= 'A' &&
// RDKit✔️✔️:     symb[0] <= 'Z') {
// RDKit✔️✔️:     symb[0] = tolower(symb[0]);
// RDKit✔️✔️:   }
//
// Rust implementation: lowercase first character of element symbol when
// the atom is aromatic and kekulization is off.
let symbol: &str = if !params.do_kekule
    && atom.is_aromatic()
    && raw_symbol.as_bytes().first().is_some_and(u8::is_ascii_uppercase)
{
    // ... lowercase logic ...
    &lowered_symbol
} else {
    raw_symbol
};
```

### Example B — Hypothetical Gemmi port

```cpp
// Gemmi source:
//   if (model.findChain(name)) {
//     chain = model.findChain(name);
//     chain->setEntityId(entity_id);
//   }
```

```rust
// Gemmi❌❌: if (model.findChain(name)) {
// Gemmi❌❌:   chain = model.findChain(name);
// Gemmi❌❌:   chain->setEntityId(entity_id);
// Gemmi❌❌: }
//
// Not yet ported — COSMolKit has no Chain entity or setEntityId equivalent.
// Blocked on: CRATE-42 (Chain/entity model).
```

### Example C — Early-port with intentional performance gap

```cpp
// RDKit source:
//   for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
//     if (colors[i] == Canon::WHITE_NODE && ranks[i] < nextRank) {
//       nextRank = ranks[i];
//       nextAtomIdx = i;
//     }
//   }
```

```rust
// RDKit✔️❌: for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
// RDKit✔️❌:   if (colors[i] == Canon::WHITE_NODE && ranks[i] < nextRank) {
// RDKit✔️❌:     nextRank = ranks[i];
// RDKit✔️❌:     nextAtomIdx = i;
// RDKit✔️❌:   }
// RDKit✔️❌: }
//
// Rust implementation uses .iter().enumerate().min_by_key() instead of
// manual loop. Behavior is equivalent but introduces an extra closure
// allocation per call. Replace with indexed loop if profiling shows
// this on a hot path.
let (idx, _) = ranks
    .iter()
    .enumerate()
    .min_by_key(|(_, rank)| **rank)
    .ok_or(...)?;
```

### Example D — Hydrogenes with typed state (simplified)

```cpp
// RDKit source (simplified):
//   int max_serial = 0;
//   for (unsigned int aidx = 0; aidx < stopIdx; ++aidx) {
//     auto *info = (AtomPDBResidueInfo *)(mol.getAtomWithIdx(aidx)->getMonomerInfo());
//     if (info) max_serial = max(max_serial, info->getSerialNumber());
//   }
```

```rust
// RDKit✔️✔️: int max_serial = 0;
// RDKit✔️✔️: for (unsigned int aidx = 0; aidx < stopIdx; ++aidx) {
// RDKit✔️✔️:   auto *info = (AtomPDBResidueInfo *)(mol.getAtomWithIdx(aidx)->getMonomerInfo());
// RDKit✔️✔️:   if (info) max_serial = max(max_serial, info->getSerialNumber());
// RDKit✔️✔️: }
//
// COSMolKit uses typed `pdb_residue_info()` accessors instead of raw
// `getMonomerInfo()` casts. The max_serial computation maps directly.
let mut max_serial = 0;
for atom in read_parts.atoms() {
    if let Some(info) = atom.pdb_residue_info() {
        max_serial = max_serial.max(info.serial_number());
    }
}
```


### Example D — Hydrogenes with typed state (simplified)

```cpp
// RDKit source (simplified):
//   int max_serial = 0;
//   for (unsigned int aidx = 0; aidx < stopIdx; ++aidx) {
//     auto *info = (AtomPDBResidueInfo *)(mol.getAtomWithIdx(aidx)->getMonomerInfo());
//     if (info) max_serial = max(max_serial, info->getSerialNumber());
//   }
```

```rust
// RDKit✔️✔️: int max_serial = 0;
// RDKit✔️✔️: for (unsigned int aidx = 0; aidx < stopIdx; ++aidx) {
// RDKit✔️✔️:   auto *info = (AtomPDBResidueInfo *)(mol.getAtomWithIdx(aidx)->getMonomerInfo());
// RDKit✔️✔️:   if (info) max_serial = max(max_serial, info->getSerialNumber());
// RDKit✔️✔️: }
//
// COSMolKit uses typed `pdb_residue_info()` accessors instead of raw
// `getMonomerInfo()` casts. The max_serial computation maps directly.
let mut max_serial = 0;
for atom in read_parts.atoms() {
    if let Some(info) = atom.pdb_residue_info() {
        max_serial = max_serial.max(info.serial_number());
    }
}
```
