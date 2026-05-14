# Agent Workflow Rules

This file documents recurring mistakes made by Hermes agent during
RDKit porting work, and the correct workflow to prevent them.

**This file's contents must be reproduced verbatim (not by path reference)
in every subagent task instruction that involves RDKit source reproduction.**

---

## 1. Golden Rule: Copy First, Code Second

The correct order of operations is NEVER:

```
写 Rust → 贴 RDKit 标记当装饰
```

The ALWAYS correct order is:

```
① 打开 RDKit C++ 源文件
② 将目标函数按行复制为 Rust 注释
③ 给每行加上 // RDKit<beh><perf>: 标记
④ 在注释块下方写等价的 Rust 实现
```

After writing the Rust, verify that the C++ block above it is COMPLETE
(no `...`, no `// see file X` placeholders).

**Verification:** At the end, grep for `\.\.\.` in RDKit marker comments.
If any exist, the protocol was violated.

## 2. Subagent Tasks Must Be Self-Contained

Never write:

```
Follow dev/source_reproduction_protocol.md
```

Always write the key RULES verbatim:

```
# CRITICAL - These rules are MANDATORY:
#
# 1. Copy Rule: C++ source lines MUST be copied VERBATIM as comments.
#    No summarising, no paraphrasing, no "see file X" references.
#
# 2. Each C++ line gets a TWO-AXIS marker (per dev/source_reproduction_protocol.md):
#
#    First axis (behavior):
#      ❌ = not reproduced
#      ❗ = approximate / unverified
#      ✔️ = behavior verified against RDKit
#
#    Second axis (perf/complexity):
#      ❌ = known to be worse
#      ❗ = unresolved after local inspection
#      ✔️ = comparable after local inspection
#      🔝 = known to be better
#
#    Valid marker combinations (choose the BEST fit):
#      // RDKit❌❌: <line>   — Not implemented at all
#      // RDKit❗❗: <line>   — Approximate, perf not reviewed
#      // RDKit❗✔️: <line>  — Approximate, perf looks comparable
#      // RDKit✔️❗: <line>  — Verified behavior, perf not reviewed
#      // RDKit✔️✔️: <line>  — Verified behavior, perf comparable
#      // RDKit✔️❌: <line>  — Verified behavior, known perf regression
#      // RDKit✔️🔝: <line>  — Verified behavior, known perf improvement
#
#    NEVER default to ✔️✔️ for a first port. Choose ❗❗ or ❗✔️
#    unless you specifically verified against RDKit output.
#
# 3. Write the Rust equivalent code BELOW each C++ comment block.
#
# 4. After writing, grep for 'RDKit.*\.\.\.' — zero matches required.
#    Also grep for 'RDKit✔️✔️' — verify each one is justified.
#
# 5. File header MUST reference dev/source_reproduction_protocol.md.
```

## 3. Do Not Trust Subagent Summaries

When a subagent reports "all functions carry RDKit status markers",
THIS IS NOT ENOUGH. You must verify by:

```
grep -n "\.\.\." <file>        # Check for summarized C++ blocks
grep -c "// RDKit✔️✔️:" <file> # Count verbatim markers
grep -c "// RDKit" <file>      # Count ALL RDKit references
```

If `...` appears in RDKit marker comments, the subagent violated the
protocol. Reject and rewrite.

## 4. What a Correct Implementation Looks Like

CORRECT (verbatim C++ → Rust):
```
    // RDKit✔️✔️: for (unsigned int layer = 0; layer < numLayers; ++layer) {
    // RDKit✔️✔️:   for (unsigned int atomIdx = 0; atomIdx < nAtoms; ++atomIdx) {
    // RDKit✔️✔️:     std::uint32_t invar = layer;
    // RDKit✔️✔️:     gboost::hash_combine(invar, currentInvariants[atomIdx]);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    for layer in 0..num_layers {
        for atom_idx in 0..n_atoms {
            let mut invar = layer as u32;
            hash_combine(&mut invar, current_invariants[atom_idx]);
        }
    }
```

WRONG (summarized):
```
    // RDKit✔️❌: for (unsigned int layer = 0; layer < numLayers; ++layer) {
    // RDKit✔️❌:   ...
    // RDKit✔️❌: }
```

## 5. Source File Attribution

Each function block must name the EXACT RDKit source file and line range:

```
// RDKit source: MorganGenerator.cpp lines 384-459
```

Not just the module name. This allows future reviewers to find the
original.

## 6. Check All Related Source Files

When porting a function, check ALL source files in the RDKit module
directory for helper functions. Example:

Drawing SVG support → check `MolDraw2DSVG.cpp` AND `MolDraw2DDetails.h`
Morgan fingerprints → check `MorganGenerator.cpp`, `FingerprintUtil.cpp`

## 7. Default Markers: ❗❗ Not ✔️✔️

When a function is first ported, the default two-axis status MUST be:

```
// RDKit❗❗: <C++ line>
```

NOT `RDKit✔️✔️`. The ❗❗ means:

- **Behavior** (first axis) = ❗ — "approximately implemented, exact behavioral
  equivalence has not been proven"
- **Performance** (second axis) = ❗ — "performance status unresolved"

Only upgrade a marker to ✔️ on either axis AFTER:

- **Behavior → ✔️**: You have verified against RDKit reference output for
  the function's full input state space.
- **Performance → ✔️**: You have done local code inspection confirming
  comparable algorithmic complexity, OR the code is structurally identical
  to the C++ original (same loop nesting, same data structure choice).

Do NOT assume equivalence just because the Rust code looks similar. Every
first-port marker is ❗❗ unless you have specific evidence for a higher
status.

### Upgrade progression

```
❗❗ → (behavior verified) → ✔️❗ or ✔️❌ or ✔️🔝
❗❗ → (perf verified)      → ❗✔️
❗❗ → (both verified)      → ✔️✔️
❗❗ → (known regression)   → ❗❌
```

### Exceptions

You may use a higher marker immediately if the Rust function is:
- A direct 1:1 translation with identical algorithm, data structures, and
  no architecture differences (e.g., pure math helpers, lookup tables)
- A thin delegate/wrapper that calls another already-✔️ function

When in doubt: default to ❗❗.

## 8. Self-Correction Checklist

Before submitting any RDKit porting work, run:

```
grep -n "RDKit.*\.\.\.\|RDKit✔️✔️" <output file> | head -20
grep -n "// see " <output file>
```

Check that RDKit✔️✔️ markers are ONLY used when behavior has been verified.
Default should be RDKit❗❗ for first ports.
```

If any match, STOP and fix.
