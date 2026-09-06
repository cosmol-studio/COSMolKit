# COSMolKit 旧 `cosmolkit-core` 功能迁移总账

状态：执行中的迁移跟踪文档。目标架构以
[`final_target_architecture.md`](./final_target_architecture.md) 和
[`crate_architecture.md`](./crate_architecture.md) 为准。

这份文档回答一个问题：**旧 `cosmolkit-core` 中的每一类能力，是否已经
进入最终的 owner、是否已经由 `cosmolkit` facade 暴露、以及旧实现是否已经
被移除。** 它不是代码行数统计，也不把“复制到新目录”算作完成。

## 1. 完成判据

一个迁移单元只有同时满足以下条件才记为 `✅ 完成`：

1. 算法实现位于目标 crate，输入是 `cosmolkit-model` 的 detached value 或
   明确的 assignment/report；实现不接收 `Molecule`、`OpParts` 或 runtime cache。
2. `cosmolkit` facade 提供 canonical API；单分子操作由 facade 提取、校验并
   提交 detached 结果。
3. canonical naming/binding registry 已注册，名称和 receiver/state model
   与 facade 一致。
4. 旧 `cosmolkit-core` 路径已经删除，或者只保留明确的内部 forwarding layer，
   不存在第二套实际算法实现。
5. 至少有 focused parity 或行为回归验证；涉及操作时使用 strict contract
   检查。

以下状态不算最终完成：

- `🟡 垂直切片`：新边界和部分功能可用，但旧实现仍在，或只完成了一个功能子集。
- `🟠 实验结构`：目标 crate、类型或签名已建立，但尚未承载完整算法。
- `⚪ 未开始`：仍完全由旧 core 提供。
- `🔴 阻塞`：依赖 runtime/Molecule、model 不变量或另一个尚未迁移的基础层。

## 2. 当前总进度

进度按 **32 个迁移单元**计算，而不是按源文件数量计算。每个单元的权重相同；
只有 `✅` 计入最终完成，`🟡/🟠` 只计入“已建立迁移边界”。

```text
最终完成（旧实现已移除）       1 / 32   [█░░░░░░░░░░░░░░░░░░░]   3%
已建立 detached/目标边界       13 / 32   [████████░░░░░░░░░░░░░]  41%
```

当前唯一可视为完整基础单元的是 `cosmolkit-types` 的公共词汇；
`cosmolkit-model`、`cosmolkit-cx`、runtime 和 descriptor 都仍处于迁移中，
因为旧 core 中仍存在重复 owner 或未完成的行为实现。

## 3. 基础层与 runtime

| 单元 | 旧 core 位置 | 最终 owner | 状态 | 说明 |
|---|---|---|---|---|
| 公共化学词汇 | `model/atom.rs`、`model/bond.rs` 中的基础枚举 | `cosmolkit-types` | ✅ | `Element`、`BondOrder`、立体/杂化枚举已抽出；后续只需清理旧 re-export。 |
| detached model | `model/atom.rs`、`bond.rs`、`adjacency.rs`、`query.rs`、`coordinates.rs`、`topology.rs` 等 | `cosmolkit-model` | 🟡 | 类型已迁移，但旧 core 仍有完整 model/`Molecule` owner；必须完成唯一 owner 清理。 |
| CX 语法记录 | `notation/cx.rs` 及 SMILES/SMARTS 内部 CX 逻辑 | `cosmolkit-cx` | 🟠 | crate 已建立，`ParsedCxExtensions` 和目标 lowering 尚未完整迁出。 |
| 环分解基础 | 旧 core 环分解调用及 `cosmolkit-ringdecomposer` | `cosmolkit-ringdecomposer` | 🟡 | 独立 crate 已存在，但旧 core 仍保留部分环 perception/调用路径，尚未完成唯一实现收口。 |
| `Molecule` / builder | `model/molecule.rs`、`model/builder.rs` | `cosmolkit` | 🟠 | facade 中已有骨架，但尚未承载旧 core 的构造、解析和完整状态。 |
| `OpParts` / operation runtime | `operations/*`、`operations/ops.rs` | `cosmolkit` | 🟠 | 新 capability context 骨架存在；registry、invariants、effects 仍未完成迁移。 |
| operation macros | `cosmolkit-macros-old` | `cosmolkit-macros` | 🟠 | 新宏已接管 canonical 名称，目前只覆盖新的 operation body/binding registry 试验；旧宏保留为迁移兼容实现。 |
| feature/binding registry | 旧 support/feature 表和新 binding contract | `cosmolkit` | 🟡 | descriptor/hydrogen 已注册；其余 public capability 尚未注册。 |

```text
基础层 + runtime 最终完成       1 / 8   [██░░░░░░░░░░░░░░░░░░░░]  13%
已建立边界                     7 / 8   [█████████████████░░░░░]  88%
```

## 4. `cosmolkit-core` 基础算法迁移

这些能力仍属于 foundational implementation layer，不强行拆成“一功能一个
crate”。它们最终进入 `cosmolkit-core`，但必须接收 model values。

| 单元 | 旧 core 模块 | 最终 owner | 状态 | 下一验收点 |
|---|---|---|---|---|
| 化合价/隐式氢 | `chemistry/valence.rs` | `cosmolkit-core` | 🟡 | 已补齐 `assign_valence_for_topology`、`assign_radicals_for_topology`、`total_num_hydrogens_for_topology` 以及 detached 质量 primitive；旧 live/cache 入口仍保留，待 runtime 迁移后收口。 |
| sanitize/property cache | `chemistry/atom_properties.rs`、`operations/sanitize.rs`、`sanitize_pipeline.rs` | `cosmolkit-core` + `cosmolkit` runtime | 🟡 | 已将 property-cache 的纯 valence assignment 收口到 `chemistry::sanitize::assign_property_cache_for_topology`；runtime 仍负责 cache 安装、失效和 sanitize stage 编排，后续继续拆其它 assignment。 |
| 加氢/去氢 | `chemistry/hydrogens.rs` | `cosmolkit-core` | 🟡 | `add_hs_assignment_for_topology`、`remove_hs_candidates_for_topology`、`remove_hs_assignment_for_topology` 已公开为 detached assignment planner；runtime apply、坐标/属性安装和完整 mapping/invariant parity 仍待收口。 |
| 芳香性/kekulize | `chemistry/aromaticity.rs`、`kekulize.rs` | `cosmolkit-core` | ⚪ | 依赖 detached valence、ring 和 runtime commit。 |
| 基础环 perception | `chemistry/rings.rs`、`ringdecomposer` 调用 | `cosmolkit-core` | ⚪ | 需要把 ring assignment/report 与 derived-cache authority 分离。 |
| 基础立体/CIP | `chemistry/stereo.rs`、`cip.rs`、`potential_stereo.rs`、`atropisomer.rs` | `cosmolkit-core` / `cosmolkit-stereo` | ⚪ | 先保留共享 primitive，枚举类高层能力另迁 `stereo`。 |
| 通用图/子图 primitive | `chemistry/subgraph.rs`、`mol_transforms.rs`、`matrices.rs` | `cosmolkit-core` | ⚪ | 只接收 topology/assignment，不再接收 live `Molecule`。 |
| InChI source port | `cosmolkit-core/src/inchi.rs` | `cosmolkit-inchi` | 🟡 | 独立 crate 已存在；旧 core wrapper、facade ownership 和 warning policy 仍需收口。 |

sanitize/property-cache 当前阶段验收：

- detached API：`cosmolkit_core::chemistry::sanitize::assign_property_cache_for_topology`
  （顶层同时导出同名函数）。
- facade/runtime：`operations::sanitize::sanitize_properties_for_topology` 保留为
  旧路径的薄 forwarding；`sanitize_pipeline` 直接调用 chemistry detached API，
  只在 runtime 写入 valence cache。
- source/parity：复用 `valence::assign_valence_with_options_for_topology` 的
  RDKit source-backed assignment；现有 detached property-cache parity 测试覆盖
  `CCO`，并通过 strict core check。
- old-core removal：sanitize pipeline 的 property-cache 计算重复路径已删除；
  其余 sanitize stages 仍依赖 `MoleculeReadParts` 和 runtime cache，尚未移除。

hydrogens 当前阶段验收：

- detached API：`add_hs_assignment_for_topology`、
  `remove_hs_candidates_for_topology`、`remove_hs_assignment_for_topology`；
  assignment 及其 atom/bond/stereo/SGroup update records 已成为可跨 crate
  传递的公开结果值。
- facade/runtime：旧 `MoleculeReadParts` 入口仍保留为兼容包装；现有 operation
  body 继续负责将 assignment 应用到 topology、coordinates、properties，并执行
  cache/effect bookkeeping。
- source/parity：planner 复用原 RDKit addHs/removeHs source-port 代码和现有
  assignment 测试；本阶段仅提升 detached 可见性，不改变算法分支。
- old-core removal：尚未删除 live/runtime apply 路径；下一步应在 `cosmolkit`
  runtime 中接入这些 detached planners，再删除重复的 read-parts 规划入口。

```text
基础算法最终完成             0 / 8   [░░░░░░░░░░░░░░░░░░░░░]   0%
已建立 detached 边界          5 / 8   [███████████░░░░░░░░░░░░]  63%
```

## 5. 领域算法 crate

| 单元 | 旧 core 模块 | 目标 crate | 状态 | 迁移范围 |
|---|---|---|---|---|
| SMILES parser/writer | `notation/smiles.rs`、`smiles_write.rs`、`canon_smiles.rs`、`canon_rank.rs`、`fragment.rs` | `cosmolkit-smiles` | 🟠 | parser state、writer、canonicalization、CXSMILES lowering。 |
| SMARTS/search/MCS | `search/query.rs`、`query_graph.rs`、`smarts_parse.rs`、`smarts_write.rs`、`substruct.rs`、`generic_groups.rs` | `cosmolkit-search` | 🟠 | `QueryGraph` 数据留 model；parser/compiler/matcher/writer/MCS 全迁 operator。 |
| 文件 IO | `io/molblock.rs`、`molfile.rs`、`sdf.rs`、`mol2.rs`、`xyz.rs`、PDB 读写 | `cosmolkit-io` | 🟠 | parser 返回 detached blocks；不构造 runtime `Molecule`。 |
| 分子量/组成/描述符 | `properties/descriptors.rs` | `cosmolkit-descriptors` | 🟡 | 三个质量/公式 API 已完成 detached vertical slice；其余描述符仍在旧 core。 |
| 指纹 | `properties/fingerprint.rs`、`avalon_fingerprint.rs` | `cosmolkit-fingerprints` | 🟠 | Morgan、Pattern、Layered、MACCS、Avalon、topological torsion 等。 |
| 立体异构体枚举 | `chemistry/stereo_enumerate.rs` | `cosmolkit-stereo` | ⚪ | 多输出 `Vec<detached candidate>`，由 runtime 组装 `Vec<Molecule>`。 |
| 互变异构体 | `chemistry/tautomer.rs`、`tautomer_transforms.rs` | `cosmolkit-tautomer` | ⚪ | 多输出和 canonical result 都必须脱离 `MoleculeDraft`。 |
| 构象/距离几何 | `chemistry/coordinates.rs`、`distgeom.rs`、`conformer_selection.rs` | `cosmolkit-conformer` | ⚪ | 输入 topology/coordinates，输出 conformer 或 assignment。 |
| 力场/MMFF/UFF | `chemistry/forcefield/*` | `cosmolkit-forcefields` | ⚪ | 参数表、能量、梯度和优化；不拥有 Molecule。 |
| alignment/RMSD | `chemistry/mol_align.rs`、`mol_align_support.rs` | `cosmolkit-alignment` | ⚪ | 输出 alignment/mapping/report。 |
| depiction | `properties/draw.rs`、`data/rdkit_depictor_template_smarts.h` | `cosmolkit-depict` | ⚪ | 布局和 SVG/PNG 输出；不回写 runtime。 |
| batch | `properties/batch.rs` | `cosmolkit-batch` | ⚪ | `BatchRecord`、线程数、并行调度和 batch result。 |
| bio/macromolecule | `bio/*`、`io/bio.rs` | `cosmolkit-bio` + `cosmolkit-io` | ⚪ | 结构值、残基信息、选择和 PDB/mmCIF 转换。 |

```text
领域算法最终完成             0 / 13  [░░░░░░░░░░░░░░░░░░░░░]   0%
已建立目标 crate/scaffold      13 / 13  [████████████████████] 100%

这里的 `100%` 只表示目标目录和 Cargo package 已建立；除 descriptor 垂直切片
外，其余 crate 仍然没有承载旧 core 的完整算法，因此不计入最终完成度。
```

## 6. 其它旧 core 能力归属

| 旧模块/能力 | 目标归属 | 当前状态 |
|---|---|---|
| `properties/mol_hash.rs` | `cosmolkit-fingerprints` 或独立 hashing primitive | ⚪ 未迁移；不能因为内部复用 fingerprint 就打开 fingerprint public API。 |
| `properties/mol_pickler.rs` | `cosmolkit-io`/serialization 层 | ⚪ 未迁移；需先确定 binary format 的 public ownership。 |
| `confseq/*` | `cosmolkit-conformer`（必要时拆 sequence helper） | ⚪ 未迁移。 |
| `source_port_helpers.rs` | `cosmolkit-core` 内部支持 | 🟡 暂留；不属于 public domain crate。 |
| `support.rs` / feature matrices | `cosmolkit` registry/runtime | 🟡 迁移中；最终由 facade feature/capability registry 维护。 |
| `data/*` 参数和 source tables | 使用它们的目标算法 crate | 🟡 暂留旧 core；随算法逐单元移动，禁止复制两份。 |

## 7. 跨语言适配器

| 单元 | 目标 owner | 状态 | 说明 |
|---|---|---|---|
| Python binding | `python` / `cosmolkit-py` | 🟠 | 迁移期间暂不要求命名检查通过；最终只能依赖 `cosmolkit`。 |
| WASM/JS binding | `wasm` / `cosmolkit-wasm` | 🟠 | 当前已有 API smoke surface；最终只能依赖 `cosmolkit`。 |
| canonical naming parity | `cosmolkit` registry + language probes | 🟡 | Rust registry 已开始使用；Python/JS probe 属于后续阶段。 |

```text
适配器最终完成                0 / 3   [░░░░░░░░░░░░░░░░░░░░░]   0%
已建立 contract/adapter 骨架    2 / 3   [█████████████░░░░░░░]  67%
```

## 8. 推荐迁移顺序

接下来不按“旧目录从上到下复制”，而按依赖闭包推进：

1. **完成 valence detached boundary**：把质量、同位素、元素符号和隐式氢
   assignment 整理成 core 的稳定基础 API；这是 descriptors、hydrogens、sanitize
   的共同依赖。

当前该步骤已完成第一阶段：detached assignment、radical、total-H 及质量查找均
有公开的 model-only 入口，并通过 `CCO` parity 与 strict core 检查。剩余工作是
让后续 hydrogens/sanitize 直接使用这些入口，并在旧 live API 被替代后删除重复
实现；因此状态仍保持 `🟡`。
2. **迁移 sanitize 与 hydrogens**：property-cache 的第一段纯 assignment 已迁入
   `chemistry::sanitize`；继续先迁其它纯 topology assignments，再迁纯 topology
   transformation 和
   assignment，再在 `cosmolkit` 中接入 `OpParts`、invariant、derived effects。
3. **迁移 SMILES + SMARTS/search**：先完成 CX parser/lowering，再把 parser
   state、QueryGraph operator 和 writer 拆出旧 core。
4. **迁移 IO**：SDF/MolBlock/XYZ/PDB 等直接以 detached model 为输入输出，
   runtime 只负责组装 `Molecule`。
5. **迁移高层算法**：fingerprints、剩余 descriptors、tautomer、stereo、
   conformer、forcefields、alignment、depict、batch、bio。
6. **收尾 runtime 和 bindings**：清理旧 `Molecule`/`OpParts` owner、移除旧
   forwarding/重复实现、更新 feature bundles，最后做 Python/JS/WASM contract
   probes。

## 9. 每个单元的验收记录模板

新增或推进一个单元时，在对应行后补充：

```text
- detached API：<函数/类型>
- facade API：<Molecule 方法或 batch 入口>
- naming registry：<semantic_id>
- source/parity：<测试或 fixture>
- old-core removal：<删除/尚未删除及原因>
- strict validation：<命令和结果>
```

迁移单元只有在 `old-core removal` 不再是“尚未删除”时，才能把状态从
`🟡/🟠` 改为 `✅`。
