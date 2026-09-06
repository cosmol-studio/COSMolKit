# COSMolKit 最终目标架构

状态：目标设计，尚未全部实现。

本文把 [`crate_architecture.md`](./crate_architecture.md) 中的所有权和
算法边界具体化为最终 crate 划分。`cosmolkit-core` 现在承载 detached 算法迁移实现；
原先的完整 runtime 混合实现已改名为 `cosmolkit-core-old`，仅作为迁移期间的兼容实现保留。

## 1. 总体分层

最终 workspace 分为四层：

```text
┌─────────────────────────────────────────────────────────────────────┐
│  Language adapters                                                  │
│  cosmolkit-py   cosmolkit-wasm   TypeScript/JavaScript declarations  │
└───────────────────────────────┬─────────────────────────────────────┘
                                │ only public facade dependency
                                v
┌─────────────────────────────────────────────────────────────────────┐
│  cosmolkit                                                         │
│  唯一 public Rust entrypoint                                        │
│  Molecule / MoleculeBuilder                                         │
│  OpParts / MultiOutputOpParts                                       │
│  operation registry / contracts / cache authority                   │
│  feature bundles / binding contract registry                        │
└───────────────┬───────────────────────────────┬─────────────────────┘
                │ extracts and commits          │ calls detached APIs
                v                               v
┌───────────────────────────────┐   ┌─────────────────────────────────┐
│  Domain algorithm crates       │   │  Foundational implementation     │
│                               │   │  cosmolkit-core                  │
│  smiles / search / io          │   │  valence / sanitize / hydrogens │
│  descriptors / fingerprints   │   │  aromaticity / ring / stereo    │
│  tautomer / stereo / conformer│   │  graph and chemistry primitives │
│  forcefields / alignment      │   └────────────────┬────────────────┘
│  depict / batch / bio          │                    │
└───────────────┬───────────────┘                    │
                └──────────────────┬────────────────┘
                                   v
┌───────────────────────────────┐   ┌─────────────────────────────────┐
│  cosmolkit-cx                  │   │  cosmolkit-model                │
│  representation-independent   │   │  detached canonical values      │
│  CX extension parser/records   │   │  Atom / Bond / IDs              │
└───────────────────────────────┘   │  Topology / Coordinates / Props  │
                                    │  Adjacency / Conformer           │
                                    │  QueryGraph / QueryAtom / Bond  │
                                    └────────────────┬────────────────┘
                                                     v
                                    ┌─────────────────────────────────┐
                                    │  cosmolkit-types                 │
                                    │  public chemical vocabulary      │
                                    │  Element / BondOrder / stereo    │
                                    │  hybridization and base enums   │
                                    └─────────────────────────────────┘
```

箭头表示 Cargo 依赖方向。所有算法和模型依赖都向下，禁止任何子 crate 反向
依赖 `cosmolkit`。

## 2. 最终 crate 数量

最终化学 runtime 和算法层共 **20 个 crate**。另外有 4 个非化学运行时或开发
crate。`cosmolkit-core-old` 不计入最终数量；它在迁移完成后删除。`cosmolkit-core`
是最终 foundational algorithm crate。

### 2.1 公共类型和共享语法：4 个

| crate | 最终职责 |
|---|---|
| `cosmolkit-types` | 可复用的公共化学词汇：`Element`、`BondOrder`、`BondStereo`、`ChiralTag`、`Hybridization` 等 |
| `cosmolkit-model` | 唯一 detached value model：`Atom`、`Bond`、ID、`AdjacencyList`、`TopologyBlock`、`CoordinateBlock`、`MoleculeProperties`、`Conformer`、`QueryGraph` 数据 |
| `cosmolkit-cx` | 与目标表示无关的 CX 扩展语法和 `ParsedCxExtensions`；不构造 `Molecule` 或 `QueryGraph` |
| `cosmolkit-ringdecomposer` | 独立的环分解图算法；只接收图/模型值，不拥有 runtime 状态 |

### 2.2 基础算法：2 个

| crate | 最终职责 |
|---|---|
| `cosmolkit-core` | 共享基础化学实现：化合价、sanitization、加氢/去氢、芳香性、基础环和立体化学、通用图变换 |
| `cosmolkit-inchi` | 官方 InChI source port 和标识符转换；不拥有 `Molecule` runtime |

`cosmolkit-core` 不是 runtime，也不是所有功能的垃圾桶。只有多个领域真正共享
且难以独立的基础化学实现才放入其中。

`cosmolkit-core` 不依赖 `cosmolkit-cx`。CX 是表示层语法，不是基础化学能力；
只有 `cosmolkit-smiles` 和 `cosmolkit-search` 需要解析 CX 扩展，因此它们各自
依赖 `cosmolkit-cx` 并负责目标表示的 lowering。当前旧 `cosmolkit-core` 中的
CX 依赖来自 SMILES/SMARTS 尚未完全拆出的迁移遗留，不能作为最终依赖方向。

### 2.3 领域算法：13 个

| crate | 最终职责 | 主要依赖 |
|---|---|---|
| `cosmolkit-smiles` | SMILES parser/writer、SMILES 到 model 的 lowering、CXSMILES lowering | `model`, `types`, `cx`, `core` |
| `cosmolkit-search` | SMARTS parser、`QueryGraphOperator`、compiler/planner、matcher、writer、MCS | `model`, `types`, `cx`, `core` |
| `cosmolkit-io` | SDF/MolBlock、PDB、XYZ、MOL2 等文件读写和 detached block 转换 | `model`, `types`, `core`, `bio`（需要时） |
| `cosmolkit-descriptors` | 分子量、TPSA、组成、拓扑描述符等 | `model`, `core` |
| `cosmolkit-fingerprints` | Morgan、MACCS、Layered、Pattern 等指纹 | `model`, `core` |
| `cosmolkit-stereo` | stereoisomer 枚举、CIP/立体结果的高层算法 | `model`, `core` |
| `cosmolkit-tautomer` | tautomer 枚举、规范化和结果报告 | `model`, `core` |
| `cosmolkit-conformer` | conformer 生成、distance geometry、坐标结果 | `model`, `core` |
| `cosmolkit-forcefields` | UFF、MMFF 及力场参数/优化 | `model`, `conformer`, `core` |
| `cosmolkit-alignment` | 分子/构象对齐、RMSD、原子映射结果 | `model`, `conformer` |
| `cosmolkit-depict` | SVG/PNG 等分子绘图和布局 | `model` |
| `cosmolkit-batch` | detached batch record、并行批处理和批量结果 | `model`，以及显式声明的领域算法依赖 |
| `cosmolkit-bio` | 蛋白质/大分子结构值、残基词汇、结构选择和基础转换 | `types`, `model` |

领域 crate 不接受 `Molecule`。单分子能力由 `cosmolkit` 的 `Molecule` 方法
提取 block 后调用；跨分子能力才以两个或多个 detached 输入作为函数参数。

### 2.4 唯一 runtime：1 个

| crate | 最终职责 |
|---|---|
| `cosmolkit` | 唯一公共 Rust facade；拥有 `Molecule`、`MoleculeBuilder`、`OpParts`、operation registry、contract checks、derived-cache authority，以及所有跨语言公共方法的 canonical surface |

`cosmolkit` 不再是薄 facade。它负责把 algorithm crate 的 detached 结果
验证后安装到 live `Molecule`，但不把算法实现复制到自身。

### 2.5 非化学支持和适配器：4 个

| crate/目录 | 最终定位 |
|---|---|
| `cosmolkit-macros` | 当前 operation/binding registry proc-macro；原实现保留为 `cosmolkit-macros-old`，待迁移完成后删除 |
| `cosmolkit-test-support` | workspace 内测试数据、fixture 和 parity 工具，`publish = false` |
| `python` / `cosmolkit-py` | PyO3 binding，只依赖 `cosmolkit`，不依赖算法 crate |
| `wasm` / `cosmolkit-wasm` | WASM/JavaScript binding，只依赖 `cosmolkit`，不重新实现 chemistry |

## 3. 逻辑所有权

```text
cosmolkit-types
  Element / BondOrder / BondStereo / ChiralTag / Hybridization

cosmolkit-model
  Atom / Bond / AtomId / BondId
  TopologyBlock / CoordinateBlock / MoleculeProperties
  AdjacencyList / Conformer
  QueryNode / QueryAtom / QueryBond / QueryGraph

cosmolkit-cx
  CxRecord / ParsedCxExtensions / CX syntax errors

cosmolkit-core and domain crates
  parser state, compiler, matcher, descriptors, transformations,
  assignments, reports, and algorithm-specific errors

cosmolkit
  Molecule / MoleculeBuilder
  OpParts / MultiOutputOpParts
  operation specs, registry, capability projection, cache authority
```

`QueryGraph` 只保存 query 数据和局部结构校验。SMARTS 的 parser、compiler、
matcher、writer、MCS 行为属于 `cosmolkit-search::QueryGraphOperator` 或其
内部函数。解析创建值的入口是 `search::parse_smarts`，不挂在 operator 上。

`ParsedCxExtensions` 只表达 CX 文本读到了什么。SMILES 和 SMARTS 各自拥有
自己的 lowering：

```text
CX text
  -> cosmolkit-cx::ParsedCxExtensions
       ├── cosmolkit-smiles  -> Molecule construction state
       └── cosmolkit-search  -> QueryGraph construction state
```

## 4. 典型运行流程

### 4.1 单分子值操作

```text
Molecule
  -> cosmolkit::OpParts checks operation access
  -> extract owned TopologyBlock / CoordinateBlock / Properties
  -> cosmolkit-core or domain crate
  -> detached blocks / assignment / report
  -> cosmolkit validates model and operation invariants
  -> commit into the authoritative Molecule
```

### 4.2 多输出操作

```text
Molecule
  -> OpParts / MultiOutputOpParts extracts read blocks
  -> tautomer/stereo/conformer algorithm
  -> Vec<detached candidate values>
  -> runtime validates every candidate
  -> runtime constructs Vec<Molecule>
```

算法 crate 不返回 `Molecule`、`MoleculeDraft`、分支句柄或 `OpParts`。

## 5. 依赖约束

允许的方向：

```text
bindings -> cosmolkit
cosmolkit -> domain crates
domain crates -> cosmolkit-core / cosmolkit-cx / cosmolkit-model
cosmolkit-core -> cosmolkit-ringdecomposer / cosmolkit-model
cosmolkit-model -> cosmolkit-types
```

其中 `cosmolkit-cx` 的允许消费者只有需要 CX 语法的表示层 crate：

```text
cosmolkit-smiles -> cosmolkit-cx
cosmolkit-search -> cosmolkit-cx
```

禁止的方向：

```text
任何算法 crate  -> cosmolkit
任何算法 crate  -> Molecule / OpParts / derived-cache authority
cosmolkit-model  -> search / parser / matcher / writer
cosmolkit-cx     -> Molecule / QueryGraph
cosmolkit-core   -> cosmolkit-cx
binding crate    -> cosmolkit-core 或其它领域实现 crate
```

所有实现 crate 共同形成无环图。某个 feature 使用另一个 crate 的内部 primitive
时，应增加明确的实现依赖，但不能因此开启另一个领域的公共 API。

## 6. Feature 与公共 API

每个细粒度 feature 只打开自己的 facade capability 和实现依赖：

```text
cosmolkit/smiles       -> cosmolkit-smiles
cosmolkit/search       -> cosmolkit-search
cosmolkit/fingerprints -> cosmolkit-fingerprints
cosmolkit/3d           -> cosmolkit-conformer + required primitives
```

`common_api`、`chemistry_api`、`3d_api`、`full` 是显式 feature bundle，不能
反向改变 crate 所有权。`hashing` 复用 fingerprint primitive 时，也不能自动
暴露 `Molecule::morgan_fingerprint`。

公共 API 只有一个语义来源：`cosmolkit` facade。Python 保持 snake_case，
JavaScript/TypeScript 将同一 semantic id 转为 camelCase：

```text
Molecule.molecular_weight
  Rust:       mol.molecular_weight()
  Python:     mol.molecular_weight()
  JavaScript: mol.molecularWeight()
```

每个迁移到 facade 的公共能力都必须先登记到 canonical naming registry
（当前实现为 `cosmolkit` 的 binding contract registry）。注册项至少包含
`semantic_id`、receiver/owner、API kind、Rust canonical name、feature 和
返回状态模型。旧实现如果使用了不符合本规范的名称，迁移时直接改为新名称；
不得为了保留旧名称而新增重复的 public alias，除非另有明确的兼容性决策。

## 7. 迁移完成判据

最终架构只有在以下条件同时满足时才算完成：

1. `cosmolkit` 是唯一 live `Molecule` owner 和 operation runtime。
2. `cosmolkit-core` 及所有领域 crate 的公开函数签名只出现 model values、
   参数、assignment、report 或 detached blocks。
3. `QueryGraph`、`Atom`、`Bond` 不再依赖 parser/matcher/writer 实现。
4. CX parser 与 SMILES/SMARTS lowering 分离。
5. `cosmolkit-core-old` 已删除，不作为并行实现长期存在。
6. 每个已迁移的 public capability 都已注册到 canonical naming registry，
   并且 Rust facade 存在与注册项一致的 owner、命名、receiver、feature 和
   state model。旧的不合理命名已在迁移时替换，不以兼容 alias 掩盖命名错误。
7. binding crate 只依赖 `cosmolkit`。Python/JavaScript/WASM 的导出检查、
   命名投影和跨语言语义向量属于后续 binding 阶段，暂不作为本阶段迁移完成的
   阻塞条件。
8. workspace 中不存在第二套 `Molecule`、第二套 `OpParts` 或算法到 runtime
   的反向依赖。
