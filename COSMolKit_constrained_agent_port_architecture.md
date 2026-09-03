# COSMolKit 的约束式 Agent Port 架构：从语义权威到运行时权威的闭环设计

> 本文不是对某一个 crate 划分方案的说明，而是试图从整体上解释 COSMolKit 当前架构背后的设计逻辑：为什么 `Molecule`、`cosmolkit-model`、算法 crate、Ops Contract、COW、strict checks、source-backed port 与 parity validation 必须被放在同一个体系中理解；以及为什么这些看似“偏重”的约束，实际上是在主动降低大规模 agent-assisted source port 的错误空间，同时把复杂度留在框架内部，而不是暴露给最终用户。
>
> 核心观点可以浓缩为一句话：
>
> **COSMolKit 不是要求 agent“足够聪明地写对代码”，而是在构造一个系统，使大量错误实现难以表达、越权修改难以发生、非法中间态难以污染真实对象，而剩余的合法但语义错误又能够被 source parity validation 系统性发现。**

---

## 1. 这不是一个普通的“把 RDKit 改写成 Rust”项目

如果目标只是把已有功能翻译到 Rust，最直接的架构其实很简单：

```text
core::Molecule
      ↑
      ├── smiles
      ├── search
      ├── fingerprint
      ├── descriptors
      ├── conformer
      └── forcefield
```

所有算法都可以直接接受：

```rust
fn algorithm(mol: &Molecule) -> Result<...>;

fn transform(mol: &mut Molecule) -> Result<...>;
```

然后通过单元测试、差分测试和 corpus parity 修复不一致。

这种架构的优势非常明显：

- 初期实现速度快；
- agent 很容易理解；
- crate 划分直观；
- 算法 port 时上下文成本低；
- 几乎不需要额外 runtime framework。

但它有一个隐藏前提：

> **每一个算法实现者都必须正确理解整个 `Molecule` 的状态语义，并自觉只读取、修改和失效自己应该处理的部分。**

对于少量由核心作者长期维护的代码，这个前提也许可以成立。

对于大规模 agent port，它不够可靠。

COSMolKit 面对的问题不是“agent 能不能把一个函数翻译出来”，而是：

```text
当几十、几百、上千个 operation
被不同时间、不同上下文、不同 agent 并行迁移时，

如何保证每一个局部实现
不会因为一个方便的捷径
逐渐破坏全局 state model？
```

因此 COSMolKit 的架构目标从一开始就不是：

> maximize local coding convenience

而是：

> **minimize the space of globally invalid implementations**

这决定了后面几乎所有设计。

---

## 2. 第一个关键决定：先定义“权威”，再定义 crate

很多多-crate设计从“功能应该放在哪个 crate”开始。

COSMolKit 更深的一层是先问：

```text
谁拥有真实 Molecule？
谁能看到哪些状态？
谁能修改哪些状态？
谁能声明 cache 失效？
谁能安装一个 transformation 的结果？
谁能决定一个中间结果是否成为 authoritative state？
```

这实际上是在设计 **authority model**，而不仅仅是 dependency graph。

当前体系中可以区分至少四种不同的权威。

### 2.1 Semantic authority

对于 source-backed port，语义来自 pinned upstream source。

不是：

```text
“这个输出看起来合理”
```

也不是：

```text
“这个 corpus case 需要得到 X，所以我们 patch 到 X”
```

而是：

```text
upstream source semantics
        ↓
Rust reproduction
```

Validation 是对 reproduction 的审计，而不是新的语义来源。

---

### 2.2 State authority

真正对用户可观察的 `Molecule` 状态由顶层 runtime 持有。

算法 crate 不拥有一个“自己的 Molecule”。

它们处理的是：

- immutable views；
- stable model values；
- detached blocks；
- assignments；
- typed transformation results。

因此系统中只有一个 authoritative molecule lifecycle。

---

### 2.3 Mutation authority

“能拿到一个值”和“有权修改整个分子”是两件不同的事情。

一个 operation 如果 contract 只允许：

```text
read:
    topology
    properties

write:
    topology
```

那么算法实现就不应该获得：

```rust
&mut Molecule
```

甚至不应该获得一个包含 coordinates、cache authority、runtime metadata 的万能内部对象。

它应该只得到与 contract 对应的 working set。

---

### 2.4 Commit authority

即使算法成功计算出一个候选结果，它也不意味着有权让这个结果成为真实 `Molecule` 状态。

算法做的是：

```text
input values
    ↓
candidate result
```

runtime 做的是：

```text
candidate result
    ↓
contract checks
    ↓
mapping/effect/invariant checks
    ↓
commit
```

**算法能够计算结果，但不能自行宣布结果有效。**

这是整个体系里非常重要的一条边界。

---

## 3. 因此 crate graph 只是 authority graph 的物理表达

一旦上述权威关系确定，crate 结构就不再是任意的。

概念上：

```text
                  ┌───────────────────────┐
                  │   Public API Layer    │
                  │ cosmolkit::Molecule   │
                  └───────────┬───────────┘
                              │
                              ▼
                  ┌───────────────────────┐
                  │ Operation Runtime     │
                  │ contracts / effects   │
                  │ validation / commit   │
                  └───────────┬───────────┘
                              │
                      authorized values
                              │
                              ▼
        ┌─────────────────────┼─────────────────────┐
        ▼                     ▼                     ▼
┌──────────────┐      ┌──────────────┐      ┌──────────────┐
│ search       │      │ fingerprint  │      │ core algos   │
│ SMARTS/MCS   │      │ descriptors  │      │ chemistry    │
└──────┬───────┘      └──────┬───────┘      └──────┬───────┘
       │                     │                     │
       └─────────────────────┼─────────────────────┘
                             ▼
                  ┌───────────────────────┐
                  │ cosmolkit-model       │
                  │ stable value model    │
                  └───────────────────────┘
```

依赖方向与权威方向是相反的：

```text
dependency:
    facade/runtime → algorithms → model

authority:
    model values have no runtime authority
    algorithms have no commit authority
    runtime owns authoritative lifecycle
```

因此，crate 拆分不是为了让代码“看起来模块化”。

它把一条架构规则固化成 Rust 无法轻易绕开的 dependency DAG：

> **lower-level implementation must not climb upward and obtain runtime authority.**

---

# 4. 为什么需要 `cosmolkit-model`

`cosmolkit-model` 不是“把所有 struct 扔进去”的传统贫血 model crate。

它存在的真正理由是：

> **为多个算法 crate 提供一个稳定、runtime-independent、可验证的共享 value vocabulary。**

它适合承载：

- atom / bond / element / identifier 等 value types；
- query representation；
- adjacency / graph value structures；
- topology blocks；
- coordinate / conformer values；
-普通 molecule properties；
- runtime-independent mappings / assignments。

但它不应该拥有：

- public authoritative `Molecule`；
- operation runtime；
- `OpParts`；
- registry / operation spec；
- capability tracking；
- cache lifecycle authority；
- commit / rollback policy；
- facade compatibility policy。

这是一个非常关键的区别：

```text
错误理解：

    nouns → model
    verbs → other crates
```

COSMolKit 的理解应该是：

```text
stable shared representations → model
runtime identity / lifecycle  → facade/runtime
domain interpretation         → capability crates
```

所以 `model` 不是“贫血业务对象仓库”，而更接近 shared IR / interchange representation layer。

---

# 5. 为什么不能简单用一个内部 `MoleculeData` wrapper

另一个看起来很自然的方案是：

```rust
// internal crate
pub struct MoleculeData {
    // everything
}

// public facade
pub struct Molecule {
    inner: MoleculeData,
}
```

然后：

```rust
impl Molecule {
    pub fn fingerprint(&self) -> Fingerprint {
        fingerprint::calculate(&self.inner)
    }

    pub fn add_hydrogens_(&mut self) -> Result<()> {
        hydrogen::transform(&mut self.inner)
    }
}
```

这个方案已经能解决很多问题：

- 用户只看到一个 `Molecule`；
- 所有功能都可以是 inherent methods；
- 用户不用 import 一堆 extension traits；
- 下层 crate 可以继续独立；
- 内部 representation 可以隐藏。

对于很多 Rust library，这已经是非常好的设计。

但它解决的是：

> **API aggregation 和 implementation hiding**

它没有解决：

> **algorithm authority control**

因为当算法拿到：

```rust
&MoleculeData
```

它理论上能观察所有内部字段。

当算法拿到：

```rust
&mut MoleculeData
```

它理论上能修改所有内部字段。

也就是说：

```text
wrapper protects inner from the user,
but not from implementation code.
```

而 COSMolKit 需要限制的恰恰包括 implementation code 本身。

---

# 6. Ops Contract：整个设计真正的中心

从这个角度看，COSMolKit 最重要的抽象不是 crate，也不是 wrapper。

而是 **Operation Contract**。

一个 operation 不只是一个函数。

它隐含着一份权限和状态转换声明：

```text
reads:
    topology
    properties

writes:
    topology

preserves:
    atom mapping
    coordinates

invalidates:
    selected derived state

requires:
    preconditions...

returns:
    transformed block / assignment / report
```

传统实现可能把这份规则写在：

- 文档；
- 注释；
- code review checklist；
- 作者脑子里。

COSMolKit 的目标则是：

> **尽量让 contract 变成 executable architecture。**

不是告诉 agent：

```text
“不要修改 coordinates”
```

而是让它的算法函数根本没有 coordinates。

不是告诉 agent：

```text
“不要直接清 runtime cache”
```

而是算法 crate 根本拿不到 cache authority。

不是告诉 agent：

```text
“失败时不要污染 Molecule”
```

而是算法运行在 detached / isolated working state 上，没有 commit 能力。

这是一种非常强的设计思想：

> **Do not rely on implementation discipline where the type graph can remove the capability entirely.**

---

# 7. 函数签名本身就是权限声明

因此，一个好的 algorithm boundary 类似：

```rust
pub fn calculate(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
) -> Result<Descriptor, Error>;
```

它不仅是在描述数据依赖。

它实际上同时声明：

```text
algorithm CAN:
    read topology
    read properties

algorithm CANNOT:
    read coordinates
    mutate topology
    mutate properties
    inspect runtime caches
    commit state
```

同样，一个 transformation：

```rust
pub fn transform(
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
    params: &Params,
) -> Result<TopologyTransform, Error>;
```

也不仅是普通 API。

它说明：

```text
这个 operation 的计算宇宙
只包含 contract 授权给它的 working values。
```

这让 Rust 的函数签名从“调用约定”升级成：

> **compile-time data-visibility boundary**

---

# 8. `OpParts` 的真正含义：contract-scoped working set

在这个体系下，`OpParts` 不应该被理解成“方便把 Molecule 拆开传参数”的公开对象。

它更接近：

> **runtime 根据 operation contract 构造的私有 working set。**

概念上：

```text
Authoritative Molecule
        │
        │ operation contract
        │
        ▼
┌─────────────────────────┐
│ OpParts / working set   │
│                         │
│ topology:    writable   │
│ properties:  readable   │
│ coordinates: absent     │
│ caches:      absent     │
│ runtime:     absent     │
└─────────────┬───────────┘
              │
              ▼
          algorithm
```

这里最有价值的状态不是：

```text
coordinates = read-only
```

而是：

```text
coordinates = not present
```

因为：

> **不存在的 capability 比“请勿使用”的 capability 更可靠。**

这对于 agent-assisted implementation 尤其重要。

---

# 9. 为什么 detached state 很重要

如果算法直接在 authoritative `Molecule` 上工作：

```text
valid Molecule
    ↓
temporary invalid state
    ↓
half-transformed state
    ↓
more mutation
    ↓
valid Molecule
```

那么中间任何：

- error；
- panic；
- early return；
- unexpected branch；
- agent port bug；

都有机会把真实对象停留在部分修改状态。

COSMolKit 把它变成：

```text
authoritative state
      │
      │ extract / project
      ▼
working state
      │
      │ algorithm may transform
      ▼
candidate result
      │
      │ validate
      ▼
verified result
      │
      │ controlled commit
      ▼
new authoritative state
```

这是一种 transaction-like state transition。

核心不在于是否真的实现数据库意义上的事务，而在于：

> **候选状态与已提交状态有明确边界。**

算法产生 candidate。

runtime 决定 candidate 是否有资格成为 authority。

---

# 10. COW：逻辑隔离不等于物理深拷贝

如果上述设计要求每个 operation 都：

```text
clone whole molecule
→ operate
→ validate
→ replace
```

那么安全性虽然很好，但性能代价可能不可接受。

COSMolKit 的关键设计之一，是把：

```text
logical isolation
```

与：

```text
physical copying
```

彻底解耦。

概念上，多个 value-style molecule 可以共享未变化的 storage：

```text
mol A                         mol B
  │                             │
  ├──── topology ───────────────┤
  ├──── coordinates ────────────┤
  └──── properties ─────────────┘
```

当 `mol B` 的某个 operation 真正写 topology：

```text
mol A                         mol B
  │                             │
old topology                new topology
  │                             │
  └──── coordinates ────────────┤
  └──── properties ─────────────┘
```

于是：

> **value semantics 并不意味着 eager deep-copy semantics。**

这使外部 API 可以非常容易推理：

```rust
let mol2 = mol.with_hydrogens()?;
```

用户知道 `mol` 不会变化。

内部则可以：

```text
unchanged blocks → shared
written blocks   → copy on write
```

---

# 11. In-place API 与 value-style API 不是两套架构

COSMolKit 对用户同时提供：

```rust
let new_mol = mol.with_hydrogens()?;
```

和显式的：

```rust
mol.with_hydrogens_()?;
```

这种设计的关键不是维护两套完全不同的算法。

真正理想的模型是：

```text
same operation contract
same algorithm semantics
same validation rules

storage strategy differs
```

对于 value-style operation：

```text
share unchanged state
copy on actual write
return a new value
```

对于 in-place operation：

```text
if storage is uniquely owned:
    reuse allocation / mutate authorized working block

if storage is shared:
    COW first

then validate and commit
```

因此：

```text
value semantics
in-place efficiency
authority control
```

并不矛盾。

这也是为什么 COW 在整个设计里不是一个孤立的性能优化，而是 **让强隔离架构具备生产可行性** 的关键组成部分。

---

# 12. Strict checks：把昂贵证明留在软件生产过程，而不是用户生产环境

重约束 architecture 的另一个常见问题是：

```text
如果每次 operation 都扫描大量 invariant，
production runtime 会不会很重？
```

COSMolKit 的答案是区分两类检查。

### 12.1 必须存在的 runtime correctness

真正关系到 API 正确性和内存/状态安全的必要检查，production 仍然保留。

### 12.2 Strict / redundant / diagnostic validation

为了 port、CI、parity、架构验证而存在的强检查，可以通过 strict feature 在开发和 release-mode validation 中打开，而在普通分发 build 中编译掉。

于是形成：

```text
                  Development / CI       Production

type boundary             ✓                  ✓
crate authority           ✓                  ✓
COW semantics             ✓                  ✓
required checks           ✓                  ✓
strict diagnostics        ✓                  optional/off
parity oracle             ✓                  ✕
large corpus              ✓                  ✕
```

这是一个非常成熟的权衡：

> **让最终二进制保留被证明过的结构，而不是保留证明过程的全部成本。**

因此“heavy architecture”不等于“heavy runtime”。

---

# 13. 三层错误防线

把 Ops Contract、strict checks 和 parity validation 放在一起后，COSMolKit 实际上建立了三种性质完全不同的防线。

## 第一层：结构性越权错误

例如：

```text
fingerprint implementation
尝试修改 coordinates
```

理想情况下它根本拿不到 writable coordinates。

这类错误由：

```text
crate boundary
type boundary
function signature
ownership/capability boundary
```

解决。

目标是：

> **让错误很难写出来。**

---

## 第二层：合法权限内的非法状态错误

例如算法确实允许修改 topology，但产生：

- broken mapping；
- 不满足 local structural invariant 的 block；
- contract 声明需要 preserve 的信息被破坏；
- derived effects 与声明不一致。

这类问题由：

```text
strict operation contract checks
intermediate/final invariant checks
runtime commit validation
```

发现。

目标是：

> **即使代码能编译，也不能把不合格的候选状态 commit。**

---

## 第三层：结构合法但 chemistry semantics 错误

最难的一类错误是：

```text
所有类型都正确
所有 invariant 都满足
状态转换也合法

但结果和 RDKit source semantics 不一致。
```

这无法靠 Rust 类型系统证明。

因此需要：

```text
source-backed reproduction
+
pinned parity oracle
+
focused regression
+
large corpus validation
```

目标是：

> **验证“合法的实现”是否仍然是“正确的 chemistry 实现”。**

三层组合起来：

```text
1. Type / capability constraints
       ↓
   防止越权

2. Contract / strict validation
       ↓
   防止非法状态

3. Source parity validation
       ↓
   防止合法但语义错误
```

任何一层单独存在都不够。

---

# 14. 为什么 source parity 不能反过来定义实现

Agent-assisted port 中一个非常危险的开发模式是：

```text
run differential test
    ↓
observe mismatch
    ↓
add special case
    ↓
test passes
    ↓
repeat
```

当 corpus 足够大，这很容易形成一个：

> observed-output-fitting implementation

它可能在已见样本上获得非常高的 parity，却越来越难回答：

```text
为什么这个 branch 存在？
它对应 upstream 哪段 semantics？
下一个没出现过的 state 应该怎样处理？
```

COSMolKit 更合理的顺序是：

```text
upstream source
      ↓
determine semantics / call path
      ↓
Rust reproduction
      ↓
focused source-backed tests
      ↓
parity validation
```

因此：

> **validation verifies the port; validation does not invent the port.**

当 mismatch 出现，正确的问题不是：

```text
“怎样 patch 这个 molecule？”
```

而是：

```text
“first semantic divergence 在哪里？”
```

然后回到：

- source behavior；
- call graph；
- parameter semantics；
- state transition；
- operation contract；
- Rust reproduction。

这显著降低了 agent 用 heuristic “修测试”的空间。

---

# 15. Agent 不是架构的可信主体

传统软件工程经常隐含依赖：

```text
资深开发者知道这个函数不能碰什么。
```

大规模 agent port 中，这不是足够强的保证。

Agent 很擅长：

- 局部完成任务；
- 根据可见接口寻找最短路径；
- 为让测试通过做合理-looking 修改；
- 使用当前 scope 中最方便的数据；
- 在缺少硬约束时引入隐式 coupling。

因此 COSMolKit 的设计不是：

```text
write a longer AGENTS.md
and tell the agent to be careful
```

而是：

> **Move architectural knowledge out of agent memory and into executable structure.**

例如：

```text
不要要求 agent 记住“不准改 coordinates”
→ 不给它 coordinates

不要要求 agent 记住“失败不能污染 molecule”
→ 不让它直接拥有 authoritative state

不要要求 agent 记住“cache 要由 runtime invalidation”
→ 不给它 cache authority

不要要求 agent 记住“结果必须满足 mapping contract”
→ commit path 自动验证

不要要求 agent 从 oracle 猜 chemistry
→ source semantics 是 implementation authority
```

这种设计尤其适合 high-throughput autonomous implementation。

---

# 16. 真正的目标：让局部实现者越来越“笨”

这是整个体系最反直觉、也最深的一点。

大型系统通常随着功能增加，会要求实现者理解越来越多的全局知识：

```text
System complexity ↑
Developer-required global knowledge ↑
```

这最终限制并行开发。

COSMolKit 试图获得另一种 scaling property：

```text
System complexity ↑↑↑

but

Agent-visible problem
    ≈ source behavior
    + explicit contract
    + narrow inputs
    + expected outputs
```

换句话说：

> **系统可以越来越复杂，但单个 port task 暴露给实现者的局部世界应该保持有限。**

理想状态下，一个 agent 不需要理解：

- 所有 cache；
- 所有 operation；
- Python API；
- 整个 feature graph；
- 所有其他 domain crate；
- Molecule 完整 lifecycle；

它只需要理解：

```text
这是 source；
这是 contract；
这是允许读取的数据；
这是允许产生的结果；
这是 validation。
```

如果系统能够长期维持这个性质，它的可持续复杂度上限会比依赖“实现者理解全局状态”的 architecture 高很多。

---

# 17. 为什么 `Molecule` 必须位于顶层 facade/runtime

Rust 的一个现实是：

> inherent methods 只能由定义类型的 crate 提供。

假设 canonical public `Molecule` 位于最低层：

```text
core::Molecule
      ↑
      ├── fingerprint
      ├── search
      ├── smarts
      ├── descriptors
      └── 3d
```

那么这些上层 crate 无法直接给 `core::Molecule` 添加 inherent methods。

只能选择：

### Free functions

```rust
fingerprint::morgan(&mol);
search::matches(&mol, &query);
```

### Extension traits

```rust
use fp::MoleculeFingerprintExt;
use search::MoleculeSearchExt;

mol.fingerprint();
mol.matches();
```

### 再增加一个 umbrella wrapper

```rust
pub struct PublicMolecule {
    inner: core::Molecule,
}
```

COSMolKit 从根上绕开了这个问题：

```text
                cosmolkit::Molecule
                       │
        ┌──────────────┼──────────────┐
        ▼              ▼              ▼
      search      fingerprint        core
        │              │              │
        └──────────────┼──────────────┘
                       ▼
                  model values
```

因为 public `Molecule` 位于 dependency graph 顶层，所以它天然可以聚合所有 lower capability：

```rust
mol.to_smiles();
mol.fingerprint_morgan(...);
mol.with_hydrogens();
mol.with_3d_conformer(...);
mol.with_uff_optimized(...);
```

这带来一个重要性质：

> **内部 crate topology 不需要泄露到用户 API。**

用户不应该为了使用一个分子库先理解：

```text
哪个方法属于哪个 crate？
哪个 extension trait 要 import？
哪个 lower-level representation 是 canonical？
```

复杂度应该由框架承担，而不是用户承担。

---

# 18. Complexity inversion：内部越重，外部反而越简单

COSMolKit 的设计可以被理解成一种 **complexity inversion**。

内部：

```text
operation specs
registries
contract matrices
OpParts
COW
working blocks
effects
mapping checks
strict validation
commit rules
feature wiring
source markers
parity suites
```

用户：

```python
mol2 = mol.with_hydrogens()
fp = mol.fingerprint_morgan(radius=2, n_bits=2048)
mol3d = mol.with_3d_conformer(params)
```

或者 Rust：

```rust
let mol2 = mol.with_hydrogens(&params)?;
let fp = mol.fingerprint_morgan(&params)?;
```

这不是“为了架构漂亮增加复杂度”。

而是有意把复杂度放在：

> **最少数、最可测试、最容易系统审查的位置**

从而让：

- 算法 port 更局部；
- 用户 API 更简单；
- crate graph 更稳定；
- runtime authority 更清晰。

一个优秀的 heavy architecture 不应该让每个调用者都变重。

它应该：

> **内部约束很重，外部使用很轻。**

---

# 19. QueryGraph：canonical identity 与 behavior facade

`QueryGraph` 是这套思想在跨-crate行为组织上的一个典型案例。

假设：

```text
model::QueryGraph
```

是多个 domain 共享的 canonical query representation。

它应该承担：

- query data；
- local construction；
- accessors；
- local structural validation。

它不应该为了 method syntax 被迫依赖：

- parser；
- matcher；
- SMARTS writer；
- fingerprint；
- planner；
- runtime。

因此不能简单把所有方法写回 `model::QueryGraph`。

但如果全变成：

```rust
search::matches(&query, &mol);
smarts::write(&query, &params);
fingerprint::pattern(&query, &params);
```

public API 又会变得割裂。

于是引入：

```rust
pub struct QueryGraphOperator<'a> {
    inner: &'a model::QueryGraph,
}
```

把 interpretation behavior 聚合：

```rust
operator.matches(...);
operator.to_smarts(...);
operator.pattern_fingerprint(...);
```

这不是为了制造第二个 QueryGraph。

恰恰相反，它遵守：

> **one canonical value identity, multiple capability views**

底层永远只有一个 `QueryGraph`。

behavior facade 只是借用它。

这避免了：

```text
model::QueryGraph
search::QueryGraph
```

两个 public domain identity 带来的 conversion、ownership 和 API 分裂。

---

# 20. Parser 为什么应该保持 construction entrypoint

`QueryGraphOperator<'a>` 的语义是：

```text
对一个已经存在的 QueryGraph 进行解释/操作。
```

而：

```rust
parse_smarts(...)
```

产生的是一个新的 owned `QueryGraph`。

它没有 existing `inner`。

因此更自然的边界是：

```rust
let query = parse_smarts(text, &params)?;
let op = QueryGraphOperator::new(&query);
let matches = op.matches(&molecule);
```

也就是：

```text
parser
  ↓ constructs

canonical value
  ↓ borrowed by

behavior facade
```

这体现了 COSMolKit 设计里一个反复出现的原则：

> **construction、representation、interpretation、runtime authority 不因为“调用语法漂亮”而混在同一层。**

---

# 21. Public feature graph 不应该等于 implementation dependency graph

大型 workspace 另一个很容易出现的泄露是：

```text
Feature A implementation needs crate B
→ public Feature A automatically exposes Feature B
```

这样长期下来：

```text
implementation dependency
=
public capability dependency
```

feature graph 会越来越难控制。

COSMolKit 的思路应该明确区分：

```text
implementation graph
```

和：

```text
public capability graph
```

如果 fingerprint 内部复用了某个 lower-level chemistry algorithm：

```text
fingerprint → lower-level implementation
```

不意味着用户开启 fingerprint feature 时必须得到另一个完全不同 domain 的 public API。

共享实现应该通过：

- lower-level internal crate；
- explicit algorithm dependency；

解决，而不是通过：

```text
“顺便把另一个 public feature 打开”
```

这与前面的设计哲学完全一致：

> **不要把内部拓扑泄露成外部语义。**

---

# 22. 为什么这一整套设计是一环扣一环的

现在可以看到，每一项设计实际上都在补另一项设计留下的缺口。

---

## 22.1 只有 source-backed port，不够

因为 agent 仍然可以在实现时：

- 越权修改 state；
- 污染 cache；
- 破坏 mapping；
- 在 error path 留下 partial mutation。

所以需要 Ops Contract 和 authority boundary。

---

## 22.2 只有 Ops Contract 文档，不够

因为 agent 可以无意违反约定。

所以 contract 必须通过：

- narrow signatures；
- unavailable state；
- crate dependency；
- private runtime authority；

尽量变成结构约束。

---

## 22.3 只有 narrow algorithm inputs，不够

算法仍然可能在允许修改的 block 内生成非法候选状态。

所以需要 strict invariant / contract validation。

---

## 22.4 只有 strict checks，不够

一个结果可以完全满足 structural invariants，却 chemistry semantics 错误。

所以需要 pinned source parity validation。

---

## 22.5 只有 detached state，不够

如果每次都 deep clone 整个 Molecule，安全性会换来不可接受的运行时成本。

所以需要 COW。

---

## 22.6 只有 COW，不够

如果为了 COW 性能又给算法完整 `&mut Molecule`，authority boundary 会重新被打穿。

所以 COW 必须发生在 contract-scoped block / runtime-controlled working state 层。

---

## 22.7 只有多 crate，不够

crate dependency isolation 不等于 mutation authority isolation。

一个独立 crate 拿到 `&mut Molecule`，依然拥有几乎整个系统的状态权限。

所以需要 state-level capability boundary。

---

## 22.8 只有安全的内部结构，不够

如果用户必须：

```rust
use cosmolkit_fp::MoleculeFingerprintExt;
use cosmolkit_search::MoleculeSearchExt;
use cosmolkit_smarts::MoleculeSmartsExt;
```

才能获得完整 API，那么内部拆分复杂度泄露给了用户。

所以需要顶层 `cosmolkit::Molecule` facade/runtime 聚合 inherent API。

---

## 22.9 只有漂亮 facade，不够

如果 facade 只是：

```rust
Molecule { inner: MoleculeData }
```

并把完整 `MoleculeData` 交给算法，API 虽然漂亮，但 authority control 仍不存在。

所以 facade 必须同时是 runtime authority boundary，而不是单纯 UI wrapper。

---

# 23. 整体闭环

因此最终得到的是一个闭环，而不是一组技巧：

```text
                     ┌────────────────────┐
                     │ Upstream Source    │
                     │ Semantic Authority │
                     └─────────┬──────────┘
                               │
                               ▼
                     ┌────────────────────┐
                     │ Operation Contract │
                     │ declared behavior  │
                     └─────────┬──────────┘
                               │
                               ▼
                ┌─────────────────────────────┐
                │ Capability-scoped Working   │
                │ Set / Narrow Model Values   │
                └─────────────┬───────────────┘
                              │
                              ▼
                     ┌────────────────────┐
                     │ Agent Port         │
                     │ Algorithm          │
                     └─────────┬──────────┘
                               │
                               ▼
                     ┌────────────────────┐
                     │ Candidate State    │
                     │ COW / detached     │
                     └─────────┬──────────┘
                               │
                               ▼
                     ┌────────────────────┐
                     │ Strict Checks      │
                     │ Contract/Invariants│
                     └─────────┬──────────┘
                               │
                               ▼
                     ┌────────────────────┐
                     │ Controlled Commit  │
                     │ Runtime Authority  │
                     └─────────┬──────────┘
                               │
                               ▼
                     ┌────────────────────┐
                     │ Parity Validation  │
                     │ Semantic Evidence  │
                     └─────────┬──────────┘
                               │ mismatch
                               └──────────────► source / contract /
                                                implementation analysis
```

这就是所谓“一环扣一环”。

任意抽掉其中某一层，系统都会退化：

```text
没有 contract
→ agent 权限边界变成约定

没有 narrow state
→ contract 只是文档

没有 detached/COW
→ error path 有污染 authority 的机会，或复制成本过高

没有 strict checks
→ structural contract 违规更晚暴露

没有 parity
→ 合法但错误的 chemistry 无法系统发现

没有 source authority
→ parity 容易变成 heuristic fitting

没有 top-level facade
→ 内部 crate complexity 泄露给用户
```

---

# 24. 与常见 Agent Port 模式的本质区别

一个典型的 agent-driven port pipeline 可能是：

```text
source/API description
      ↓
agent implementation
      ↓
compile
      ↓
tests
      ↓
differential mismatch
      ↓
patch
```

COSMolKit 更接近：

```text
source semantics
      ↓
explicit operation contract
      ↓
architecturally restricted implementation space
      ↓
agent implementation
      ↓
strict structural validation
      ↓
controlled state transition
      ↓
source-backed parity validation
```

可以对比：

| 维度 | 常见 Agent Port | COSMolKit 目标 |
|---|---|---|
| semantic authority | tests / observed output 容易成为事实标准 | pinned upstream source |
| algorithm input | 整个 domain object 很常见 | contract-scoped values |
| mutation authority | 依赖实现者自觉 | 结构化限制 |
| commit authority | algorithm 常直接修改对象 | runtime 独占 |
| error isolation | 依赖 careful mutation | detached/COW candidate state |
| structural checks | 普通 unit tests | strict operation/invariant checks |
| semantic checks | differential tests | source-backed parity |
| mismatch 修复 | 容易 case-driven patch | first semantic divergence |
| user API | 容易暴露 crate/trait topology | single facade |
| scale assumption | agent 需要理解越来越多上下文 | agent-visible context 尽量有界 |

真正的差异不是“用了更多测试”。

而是：

> **COSMolKit 同时压缩错误产生空间，并扩大错误检测空间。**

---

# 25. 一个有用的风险模型

可以把 correctness confidence 粗略理解为四个乘数，而不是一个 pass rate：

```text
Correctness Confidence
    ≈
    Semantic Constraint
    × Implementation Constraint
    × State-transition Constraint
    × Validation Strength
```

其中：

### Semantic Constraint

```text
实现是否真正来自 source semantics？
```

### Implementation Constraint

```text
算法是否只能访问 contract 允许的世界？
```

### State-transition Constraint

```text
错误候选状态是否能污染 authoritative state？
```

### Validation Strength

```text
合法结果是否真的与 reference semantics 一致？
```

很多 port 主要强化最后一个。

COSMolKit 的特殊之处在于：

> **四个方向同时加约束。**

---

# 26. “Zero-cost abstraction” 在这里真正意味着什么

这里的 zero-cost 不是说整个 architecture 没有成本。

它真正应该表达的是：

> **开发期为了证明正确性增加的冗余检查，不必全部转化为 production runtime overhead。**

production 仍然保留：

- Rust ownership；
- type boundary；
- crate DAG；
- public/private visibility；
-必要 runtime invariants；
- COW state mechanics。

而：

- expensive strict scans；
- redundant contract assertions；
- large parity corpus；
- reference comparison；

可以留在 validation builds / CI。

因此最理想的状态是：

```text
high verification cost during software production
+
low unnecessary verification overhead during software use
```

这比单纯追求“所有检查 production 都开着”更合理。

---

# 27. 这套架构真正“巧妙”的地方

它的巧妙不在于任何一个技术特别新。

COW、facade、contract、DAG、parity test、IR、transaction-like commit 都是已有思想。

真正巧妙的是 **组合方式**。

### 27.1 用 crate graph 约束 dependency authority

不是仅仅整理文件。

### 27.2 用 function signature 约束 data authority

不是仅仅描述输入。

### 27.3 用 detached/COW working state 约束 mutation blast radius

不是简单 clone。

### 27.4 用 runtime commit 约束 state authority

算法不拥有最终真相。

### 27.5 用 strict mode 提高开发期证明强度

不强迫用户支付全部成本。

### 27.6 用 source parity 约束 chemistry semantics

不让“类型正确”被误认为“化学正确”。

### 27.7 用 top-level facade 吸收内部复杂度

不让用户为内部 architecture 付 UX 成本。

这些约束之间不是重复关系，而是正交关系。

每一层都覆盖上一层无法证明的错误。

---

# 28. 设计最深的一点：把“人/Agent 应该记住的规则”变成“系统不允许忘记的规则”

这是我认为 COSMolKit 当前设计最值得长期保留的思想。

传统规则：

```text
sanitize 不应该改 coordinates
fingerprint 不应该清 cache
某个 transform 失败必须恢复 state
某个 operation 必须 preserve mapping
某个 feature 不应该暴露另一个 public capability
```

如果只存在于：

- README；
- AGENTS.md；
- code review；
- 作者记忆；

随着代码规模和 agent 数量增长，它们一定会衰减。

更强的方法是尽可能把这些规则映射为：

```text
crate dependency
type ownership
private visibility
generated wrapper
contract matrix
working-set projection
return type
strict validation
commit gate
feature wiring
```

因此架构本身成为一种“外部记忆”。

可以把原则写成：

> **Architectural knowledge should live in executable structure whenever possible, not only in the memory of the implementer.**

对 agent-heavy development 来说，这可能比任何 prompt engineering 都更重要。

---

# 29. 这种架构提高的是“可持续复杂度上限”

这种设计不一定在第一个月实现最快。

甚至相反，早期需要投入大量精力设计：

- operation contracts；
- model boundaries；
- runtime flow；
- generated registries；
- COW semantics；
- validation modes。

但它换来的不是某一个 feature 的速度。

它换来的是：

> **当系统规模扩大时，单个新 operation 不必获得同比增长的全局认知。**

因此它优化的是：

```text
long-term sustainable complexity
```

而不是：

```text
short-term local velocity
```

当未来出现：

```text
1000+ operations
dozens of capability domains
complex cache/mapping semantics
multiple language bindings
many parallel agent port tasks
```

理想的新增任务仍然应该能够被描述为：

```text
1. 这是 upstream source path
2. 这是 operation contract
3. 这是你能看到的 input
4. 这是允许产生的 output
5. 这是 strict verification
6. 这是 parity evidence
```

如果能够长期维持这个性质，架构就真正证明了自己。

---

# 30. 这套体系仍然有哪些真正需要警惕的地方

强约束 architecture 不等于自动正确。

它只是把 risk 移动到了更少数、更核心的位置。

这些位置应该被视为类似 trusted computing base。

---

## 30.1 Operation Contract 可能写错

如果 upstream behavior 实际上：

```text
writes A + B
```

而 contract 写成：

```text
writes A
```

那么系统可能非常严格地执行一个错误 specification。

因此 mismatch taxonomy 必须允许：

```text
algorithm bug
contract bug
runtime bug
reference misunderstanding
```

而不是默认所有 mismatch 都是算法 bug。

---

## 30.2 Registry / generator 变得非常关键

如果 operation specs、contract matrices 和 wrappers 大量自动生成，那么 generator 本身需要：

- 小；
- 稳定；
- 可 review；
- snapshot tested；
- structural tested；
- 尽量 deterministic。

因为一旦所有 operation 都依赖 generator，它实际上进入了 architecture TCB。

---

## 30.3 COW aliasing semantics 必须极其清楚

必须明确：

- 哪些 block 可以共享；
- 什么时候 ensure uniqueness；
- 哪些 derived state 跟随 block；
- in-place 与 value-style 是否严格保持同一语义；
- error/panic 时 authoritative state 是否保持合同要求。

COW 解决性能问题，但不能成为绕开 contract 的“隐形 mutable backdoor”。

---

## 30.4 Strict checks 被关闭后不能改变 semantics

Strict mode 应该主要增加：

```text
verification
diagnostics
redundant assertions
```

而不能让：

```text
strict build
```

和：

```text
production build
```

走本质不同的 chemistry algorithm path。

否则验证的是另一个程序。

---

## 30.5 Model crate 必须防止变成 dumping ground

一个类型进入 `cosmolkit-model` 的理由应该是：

> **它是 stable、runtime-independent、且被多个下层算法/domain crate meaningful shared 的 value representation。**

而不是：

```text
“它是一个 struct，所以放 model。”
```

随着项目扩大，应根据真实 dependency clusters 决定是否进一步拆分，而不是为了 crate 数量而拆分。

---

# 31. 可以写进项目文档的核心架构公理

下面这些原则几乎可以作为长期 architectural invariants。

## Axiom 1 — One authoritative runtime identity

```text
There is one authoritative public Molecule lifecycle.
```

---

## Axiom 2 — Algorithms do not own runtime authority

```text
Algorithm crates compute values;
they do not own Molecule lifecycle, cache authority, or commit authority.
```

---

## Axiom 3 — Capability follows contract

```text
An implementation should receive no more state capability
than its operation contract requires.
```

---

## Axiom 4 — Candidate state is not authoritative state

```text
A computed transformation becomes authoritative
only after runtime validation and commit.
```

---

## Axiom 5 — Logical isolation does not require eager copying

```text
COW and structural sharing are storage strategies;
they must not weaken authority boundaries.
```

---

## Axiom 6 — Validation verifies semantics; it does not invent them

```text
Pinned upstream source defines supported source-backed behavior.
Parity evidence audits the reproduction.
```

---

## Axiom 7 — One canonical value identity

```text
A shared domain representation should have one canonical value type.
Behavior facades may borrow it; they should not duplicate its identity.
```

---

## Axiom 8 — Internal topology is not public API

```text
Users should consume cosmolkit,
not reconstruct the workspace dependency graph in their imports.
```

---

## Axiom 9 — Expensive verification belongs primarily to development

```text
Retain structural guarantees in production;
compile out purely diagnostic strict checks when safe.
```

---

## Axiom 10 — Architecture should constrain agents, not trust them

```text
If an invariant can be encoded in ownership, visibility, types,
contracts, or generated structure, do not rely only on prompt discipline.
```

---

# 32. 最终总结

COSMolKit 当前架构最值得强调的，不是：

```text
我们拆了很多 crate。
```

也不是：

```text
我们做了大量 parity tests。
```

更不是：

```text
我们使用 agent port RDKit。
```

真正独特的设计是：

> **把 agent-assisted source reproduction 视为一个需要被约束的生产系统，而不是单纯的代码生成过程。**

整个体系围绕两个问题展开：

```text
How do we reduce the probability of producing an invalid port?

How do we maximize the probability of detecting a remaining invalid port?
```

前者通过：

```text
authority separation
crate DAG
Ops Contract
narrow inputs
private runtime
detached/COW state
controlled commit
```

实现。

后者通过：

```text
strict contract checks
invariant validation
source traceability
focused regression
pinned parity
large corpus validation
first-divergence debugging
```

实现。

而顶层 facade 又把所有内部复杂度吸收掉：

```text
内部：
    heavy constraints

算法实现者：
    narrow local world

用户：
    simple Molecule API
```

因此这套 architecture 的真正价值可以概括为：

> **复杂度被集中到框架；权限被压缩到 operation；语义被锚定到 source；错误被限制在 candidate state；验证被推到最大；用户则只看到一个稳定而自然的 Molecule。**

它不是为了证明“agent 可以替代工程设计”。

恰恰相反。

它说明：

> **agent 的吞吐量越高，工程设计越不能依赖 agent 自己具备全局判断力。**

最理想的 agent-native architecture，不是让 agent 拥有更多自由。

而是让它在一个足够狭窄、足够明确、足够可验证的局部世界里高速工作。

从这个角度看，COSMolKit 当前正在形成的不只是一个 cheminformatics crate architecture。

它更接近一套可以复用的方法论：

# **Constrained Agent-Assisted Source Reproduction**

其目标不是“生成更多代码”，而是：

> **在大规模自动化实现中，系统性压缩错误自由度，同时保留 source-level semantic fidelity、production performance 与最终用户 API 的简洁性。**

这才是这些设计一环扣一环之后真正形成的整体。
