# Crate Dependency Comparison

This note explains the difference between the workspace as it exists today
and the target architecture in [`crate_architecture.md`](./crate_architecture.md).
It is intentionally a learning document: the trees show ownership and
dependency direction, not every Cargo feature or every internal module.

## 1. The Three Layers

The target design separates three different responsibilities:

```text
cosmolkit-types
    public vocabulary: Element, BondOrder, Hybridization, ...

cosmolkit-model
    detached values: Atom, Bond, TopologyBlock, QueryGraph, Conformer, ...

cosmolkit
    live Molecule, OpParts, operation runtime, public API assembly
```

Algorithms sit below the runtime and above the values:

```text
algorithm crate
    receives model values
    returns a value, assignment, report, or transformed detached block

runtime crate
    extracts values from Molecule
    calls the algorithm
    validates and commits the result
```

The important distinction is between a detached value and live molecule
state. An algorithm may edit a detached `TopologyBlock`; it must not receive
the authority to replace the state of a live `Molecule` or edit its caches.

## 2. Current Workspace (Actual Code)

The current Cargo-level dependency graph is approximately:

```text
                         python / language adapters
                                  |
                                  v
                             cosmolkit
                         public re-export facade
                          /          |          \
                         v           v           v
               cosmolkit-core  cosmolkit-model  cosmolkit-types
               (large crate)       |             ^
                 |                  |             |
                 |                  v             |
                 |            cosmolkit-types ----+
                 |
                 +--> cosmolkit-model
                 +--> cosmolkit-types
                 +--> cosmolkit-macros
                 +--> cosmolkit-ringdecomposer
                 +--> cosmolkit-inchi (optional)
```

Inside `cosmolkit-core`, many roles are still colocated:

```text
cosmolkit-core
├── model::molecule       Molecule and live state
├── operations::ops       OpParts, registry, contracts, runtime
├── model::*               aliases and shared value access
├── chemistry::*           valence, hydrogens, sanitize, stereo, 3D, forcefields
├── search::*              SMARTS parser, matcher, writer, QueryGraph operator
├── properties::*          descriptors, fingerprints, batch, hashes
├── notation::*            SMILES and writers
└── io::*                  SDF, PDB, mmCIF, XYZ, ...
```

So the current situation is not:

```text
cosmolkit -> core -> model
```

only. It is also:

```text
cosmolkit-core contains both the runtime owner and the algorithm owner
```

The query extraction slice has already changed one part of this picture:

```text
cosmolkit-model::query
    owns QueryNode, predicates, QueryAtom, QueryBond, QueryGraph

cosmolkit-core::search
    owns parsing, matching, writing, compilation, and operator behaviour
```

The query implementation now has one canonical model owner. `QueryAtom` and
`QueryBond` explicitly compose their concrete local values with predicates, and
search behavior consumes the graph without projecting it back into a concrete
`Molecule`.

## 3. Target Workspace

The target dependency direction is:

```text
                         python / language adapters
                                  |
                                  v
                             cosmolkit
                  Molecule + OpParts + public API/runtime
                         /          |           \
                        v           v            v
             cosmolkit-core   cosmolkit-search   other domain crates
             foundational       (future split)   descriptors, fp, io, ...
                    |              |             |
                    +--------------+-------------+
                                   v
                          cosmolkit-model
                                   |
                                   v
                          cosmolkit-types
```

The arrows mean “depends on”, so every arrow points downward toward simpler
values or lower-level algorithms. No child crate points back to `cosmolkit`.

The target roles are:

```text
cosmolkit-types
    shared public vocabulary and enums

cosmolkit-model
    Atom, Bond, IDs, adjacency, conformers, blocks, QueryGraph data
    local structural validation only

cosmolkit-core
    reusable foundational chemistry algorithms over model values
    valence, sanitization, hydrogen transforms, and shared primitives
    no Molecule, OpParts, registry, or runtime cache authority

cosmolkit-search       (future domain split)
    SMARTS parser, QueryGraphOperator, matcher, writer, MCS
    depends on model and any genuinely reused lower-level algorithm crate

cosmolkit
    live Molecule and MoleculeBuilder
    OpParts and operation registry/contracts
    extract -> call algorithm -> validate -> commit
```

## 4. One Concrete Call Compared

### Current SMARTS path

The current path is approximately:

```text
public cosmolkit API
        |
        v
cosmolkit-core::search::parse_smarts
        |
        +--> parser state / QueryGraphBuilder
        v
cosmolkit-model::QueryGraph
```

Writers, parsers, and matchers use `QueryGraph` directly. Concrete-molecule
fingerprints retain their separate `&Molecule` input contract; query-target
fingerprinting, if added later, must define an explicit query input rather than
reintroducing a projection adapter.

### Target SMARTS path

The target path is:

```text
SMARTS text
    |
    v
search parser state + QueryGraphBuilder
    |
    | one controlled lowering
    v
cosmolkit-model::QueryGraph
    |
    v
QueryGraphOperator / matcher / writer
```

The parser and writer operate on the query value directly. Algorithms that
need a concrete target explicitly receive `Molecule`; no generic
No `QueryGraph -> Molecule` adapter is retained.

### Target topology operation path

For a mutating operation such as hydrogen removal:

```text
cosmolkit::Molecule
    |
    | runtime checks registry access
    v
cosmolkit::OpParts
    |
    | extract authorized detached blocks
    v
cosmolkit-core::remove_hydrogens(
    TopologyBlock,
    CoordinateBlock,
    MoleculeProperties,
)
    |
    v
transformed detached blocks
    |
    | runtime validates invariants/effects and commits
    v
cosmolkit::Molecule
```

The child algorithm never sees `Molecule`, `OpParts`, derived caches, or the
operation registry. The parent runtime remains the only place that can make
the result authoritative.

## 5. Why `core` Is Not the Runtime in the Target

`cosmolkit-core` can still be the home of tightly coupled foundational
chemistry such as valence, sanitization, and hydrogen transforms. “Core” here
means reusable chemistry implementation, not ownership of the live object.

The ownership split is therefore:

```text
model  = what a detached molecular value is
core   = how foundational chemistry algorithms transform or inspect it
runtime = who is allowed to install a result into a live Molecule
```

Keeping these meanings separate prevents every algorithm from depending on a
large runtime object merely because it needs atom or bond data.

## 6. Feature Flags Versus Crate Boundaries

Feature flags select capabilities inside a crate; they do not redefine data
ownership.

```text
fine-grained feature:
    cosmolkit/fingerprints -> exposes fingerprint API

crate boundary:
    cosmolkit-fingerprints -> depends on model/core values
```

An implementation dependency must not automatically expose an unrelated
public capability. For example, a hashing implementation may reuse a
fingerprint primitive, but enabling `hashing` must not implicitly expose the
Morgan fingerprint API unless that bundle explicitly requests it.

## 7. QueryGraph 的逻辑依赖树

这一节描述的是“谁依赖谁的数据和语义”，不是 crate 之间的 Cargo
依赖。可以把 `QueryGraph` 看成根节点，把它内部的结构和外围行为展开：

### 7.1 当前逻辑结构

当前 `QueryGraph` 的数据关系大致是：

```text
QueryGraph
├── atoms: Vec<QueryAtom>
│   ├── atom: Atom<QueryNode<AtomQueryPredicate>>
│   │   ├── concrete atom fields
│   │   └── query payload（旧模型遗留）
│   └── predicate: QueryNode<AtomQueryPredicate>
│       └── AtomQueryPredicate
│           ├── AtomicNumber / IsAromatic / FormalCharge / ...
│           ├── Range / RecursiveSmarts / ...
│           └── And / Or / Xor / Not（由 QueryNode 组织）
├── bonds: Vec<QueryBond>
│   ├── bond: Bond<QueryNode<BondQueryPredicate>>
│   │   ├── endpoint AtomId
│   │   ├── representative BondOrder
│   │   └── query payload（旧模型遗留）
│   └── predicate: QueryNode<BondQueryPredicate>
│       └── BondQueryPredicate
│           ├── Order / IsAromatic / Direction / Stereo / ...
│           └── And / Or / Xor / Not
├── adjacency: Vec<Vec<(atom_index, bond_index)>>
├── props: BTreeMap<String, String>
├── conformers_2d / conformers_3d
└── stereo_groups
```

当前外围行为依赖关系是：

```text
SMARTS 文本
    ↓
parser state + QueryGraphBuilder
    ↓ finish()
QueryGraph
    ├── query predicate matching
    │       ├── QueryNode evaluation
    │       ├── AtomQueryPredicate evaluation
    │       └── BondQueryPredicate evaluation
    │
    ├── substructure matcher
    │       └── 读取 QueryGraph.atoms / bonds / adjacency
    │
    ├── SMARTS writer
    │       └── 读取 predicate tree、方向、映射和 properties
    │
    ├── CompiledQuery
    │       └── 当前只是保存 QueryGraph 的 wrapper
    │
    └── pattern fingerprint
            └── concrete Molecule -> concrete fingerprint path
```

Pattern fingerprinting intentionally remains a concrete-molecule operation;
query data is not materialized as a fake molecule to reach it.

### 7.2 目标逻辑结构

目标结构把 value data 和解释这些 data 的行为分开：

```text
QueryGraph                         // 纯 query value
├── QueryAtom
│   ├── query atom attributes
│   └── AtomQueryPredicate
├── QueryBond
│   ├── endpoint IDs / bond attributes
│   └── BondQueryPredicate
├── adjacency
├── properties
├── conformers / stereo groups
└── local graph validation
```

目标行为关系是：

```text
SMARTS 文本
    ↓
parser state + QueryGraphBuilder
    ↓ controlled lowering
QueryGraph
    ├── QueryGraphOperator
    │   ├── to_smarts()
    │   ├── to_cx_smarts()
    │   ├── compile()
    │   └── matches(Molecule)
    │
    ├── matcher
    │   └── QueryGraph + concrete Molecule -> MatchResult
    │
    ├── writer
    │   └── QueryGraph -> SMARTS/CXSMARTS
    │
    ├── compiler
    │   └── QueryGraph -> CompiledQuery/MatchPlan
    │
    └── MCS
        └── input molecules -> QueryGraph/McsResult
```

这里的关键变化是：

```text
QueryGraph -> matcher/writer/compiler
```

而不是：

```text
QueryGraph -> concrete Molecule -> matcher/writer/compiler
```

### 7.3 具体类型之间的依赖方向

可以按下面的顺序理解各类型：

```text
QueryNode<T>
    ↓ 组织递归布尔逻辑
AtomQueryPredicate / BondQueryPredicate
    ↓ 描述单个 atom / bond 的匹配条件
QueryAtom / QueryBond
    ↓ 携带 query 条件、ID、局部属性
QueryGraph
    ↓ 携带节点、边、adjacency 和图级属性
QueryGraphOperator / matcher / writer / compiler
    ↓ 解释 QueryGraph
MatchResult / CompiledQuery / McsResult
    ↓ 输出操作结果
```

反方向不应存在：

```text
QueryGraph
    X 不应要求 concrete Molecule 才能表达自身
QueryNode / predicate
    X 不应依赖 parser、matcher 或 writer
QueryGraph data
    X 不应拥有 runtime cache 或 OpParts
```

### 7.4 `Molecule` 在逻辑树中的位置

`Molecule` 与 `QueryGraph` 是两个并列的图 value，不是父子关系：

```text
Molecule                         QueryGraph
├── concrete Atom/Bond           ├── QueryAtom/QueryBond
├── concrete chemistry state     ├── predicate semantics
├── runtime-owned caches         ├── query adjacency
└── OpParts/runtime authority    └── local query validation
          \                         /
           \                       /
            matcher / alignment / search behaviour
```

matcher 可以同时读取二者，但不通过把 `QueryGraph` 转成 `Molecule`
来完成匹配。需要 concrete target 的算法明确接收 `Molecule`，需要 query
semantics 的算法明确接收 `QueryGraph`。

## 8. Reading the Migration Status

The current query slice has completed the first ownership step:

```text
[done] query data has one canonical owner: cosmolkit-model
[done] core/search no longer defines a second QueryGraph
[done] no-default core compilation passes
[done] remove the broad QueryGraph -> Molecule projection
[done] make SMARTS parsing direct-to-QueryGraph
[open] remove query payload compatibility baggage from concrete Atom/Bond
[future] move live Molecule and OpParts from core into cosmolkit
```

这就是为什么下一步应按 vertical slice 推进：先清理 query value 的逻辑
路径，再让一个低耦合算法通过 extract/call/validate/commit 边界，最后再
迁移更大的 runtime owner。
