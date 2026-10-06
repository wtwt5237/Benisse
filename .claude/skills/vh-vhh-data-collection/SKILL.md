---
name: vh-vhh-data-collection
description: Specification for collecting and preprocessing antibody heavy-chain variable domain (VH/VHH, nanobody) data from the PDB — full submitted sequences (SEQRES/mmCIF, never reconstructed from ATOM records), AbNumber Chothia numbering, explicit sequence-to-structure residue mapping, 1-based renumbered structures, backbone [N, CA, C, O] extraction with missing-residue/atom masks, and QC outputs. Use this whenever the task involves building a VH/VHH/nanobody structure dataset, extracting heavy-chain variable domains from PDB/mmCIF, Chothia numbering, residue mapping, renumbering structures from 1, or extracting backbone coordinates with missing annotations — even if the user doesn't mention this spec. 抗体重链可变区 / VHH / 纳米抗体 数据收集、预处理、Chothia 编号、残基映射、主链提取、缺失标注。
---

# 抗体重链可变区数据收集与预处理规范

你负责按照本文件实现、运行并核验数据收集流程，产出可追溯的序列、结构、残基映射和缺失标注。先阅读当前项目已有代码与配置，复用可用实现；不要只给计划或伪代码。不要擅自扩展为结构预测、缺失环建模或模型训练任务。

## 1. 必须遵循的流程

**查找候选数据 → 提取重链及其完整提交序列 → 用 AbNumber 定位并提取 VH/VHH → Chothia 编号 → 建立序列—结构对应关系并按 Chothia 对齐 → 导出以 1 为基的结构 → 最后提取主链，保留所有缺失标注。**

最重要的约束：

1. **序列直接下载，或从 PDB 的 SEQRES / mmCIF 的聚合物序列字段读取。不能把 ATOM 中有坐标的残基拼起来当作完整序列。**
2. **完整序列是位置基准。结构缺失不意味着序列缺失；不能删除这些序列位置。**
3. **Chothia 编号、完整链序号、可变区内从 1 开始的序号、原始结构残基编号必须分别保存，并建立显式映射。**
4. **先完成映射，再重编号、裁剪和导出。不能靠两个列表长度相近、残基编号相同或 `zip(sequence, residues)` 推断对应关系。**
5. 缺失残基和缺失主链原子都要标注。不得插值、复制相邻坐标，或用预测坐标冒充实验观测。

## 2. 数据范围与检索

- 以用户指定的数据来源、PDB ID 列表和筛选条件为准。未指定来源时，可从 SAbDab 的抗体注释与 RCSB PDB 检索候选，再从 PDB 获取原始序列和结构。
- 按本流程收集重链可变区；截图中的目标称为 VHH。**AbNumber 判定为重链，只能确认 H 类型，不能单独证明该结构域是 VHH。** 记录 `domain_type=VHH/VH/unknown` 及证据；若任务严格限定 VHH，只把有明确来源注释或文献证据的 VHH 纳入正式集合，其余单列待核验。
- 不要假设重链的 PDB chain ID 一定是 `H`，也不要因为结构中没有轻链就直接认定为 VHH。链选择需结合数据库注释和序列识别。
- 优先使用原始 PDBx/mmCIF；PDB 格式可作为兼容输入。候选索引或第三方重编号结构不能替代原始数据来源。
- 记录 PDB ID、entity ID、`label_asym_id`、`auth_asym_id`、model ID、domain ID、实验方法、分辨率（若适用）、来源链接、获取日期及原始文件校验和。
- 保存查询条件、分页/返回总数、数据库快照信息和每个候选的处理状态；不得只保留成功样本。
- 未指定的物种、分辨率、序列长度、是否含抗原、完整度阈值等写入配置，注明未限定，不要擅自添加隐藏筛选条件。实验结构与计算模型分开，默认正式集合使用实验结构。
- 同一序列可能对应多个构象和实验结构，不能只因序列相同就删除结构。先记录重复组；具体保留策略服从项目配置。

## 3. 获取完整序列：禁止从坐标反推后替代

这里的“完整序列”指**该结构条目中实际提交的聚合物/实验构建体序列**，不是自动替换成天然蛋白、UniProt canonical 序列或推测的完整生物学序列。

优先顺序：

1. 直接下载目标 PDB polymer entity / chain 对应的完整 FASTA，并核对其 entity 与链的对应关系。
2. 若不能直接下载，从 mmCIF 的 `_entity_poly.pdbx_seq_one_letter_code_can` 读取，并用 `_entity_poly_seq` 的逐位残基记录核验；保留原始非标准残基身份。
3. 对 legacy PDB，从所选链的 `SEQRES` 读取。

RCSB 官方下载接口示例模板（执行时核对当前服务，不要使用下列占位符发请求）：

```text
https://www.rcsb.org/fasta/entity/{PDB_ID}_{ENTITY_ID}/download
https://www.rcsb.org/fasta/chain/{PDB_ID}.{LABEL_ASYM_ID}/download
https://files.rcsb.org/download/{PDB_ID}.cif
```

注意：FASTA entity 可能对应多条结构链；chain 下载接口使用 `label_asym_id`，不要混用作者链名 `auth_asym_id`。保存原始 FASTA header 和解析后的映射。

- `full_chain_sequence` 保存完整提交序列；AbNumber 提取的部分另存为 `variable_sequence`。
- 从坐标读取的 `observed_sequence` 只用于核验或回退对齐，不能成为完整序列的替代品。
- 仅有 ATOM 序列而没有可信的完整提交序列时，标记 `missing_full_sequence`，继续查找来源；未解决前不进入正式集合。
- 多个来源冲突时，先核对版本、链、entity、标签和构建体差异；不得悄悄挑一个或用参考蛋白序列补齐。
- 非标准残基按明确的母体氨基酸/化学组分映射处理，并记录映射依据；未知残基保留为 `X` 及原始身份，不能删掉或猜成具体氨基酸。对无法可靠编号或映射的样本单列。
- 构建体本身截短与“序列有此残基但结构无坐标”是两类情况，必须区分；不能凭参考序列补造提交序列。

## 4. 提取 VH/VHH 并进行 Chothia 编号

对完整提交序列执行 AbNumber，识别重链可变结构域，保留其在原始序列上的起止位置。不要先按结构中可见残基裁剪，也不要用固定长度或固定作者编号截取可变区。

单结构域、可直接解析的序列可使用以下入口；具体 API 需与实际安装版本核对：

```python
from abnumber import Chain

numbered = Chain(full_chain_sequence, scheme="chothia", cdr_definition="chothia")
assert numbered.chain_type == "H"
variable_sequence = numbered.seq
numbered_positions = list(numbered)  # 按序列顺序保存 (Position, amino_acid)
assert "".join(aa for _, aa in numbered_positions) == variable_sequence
```

- 带标签、信号肽、恒定区或融合片段的输入，要明确哪些部分被保留、哪些被剔除，并保存原始边界。检查 N 端前缀，不能只处理 C 端 tail。
- 对 scFv、多价 VHH、串联可变域等多结构域输入，使用经过验证的多域解析方式，给每个域独立的 `domain_id` 和边界；禁止静默只拿第一个域或把多个域拼成一个。
- 只有域序列在完整序列中具有唯一匹配时，才可用字符串定位边界；出现重复匹配必须借助域解析结果确定，不能无条件取 `.find()` 的首个结果。
- 默认编号方案和 CDR 定义均为 Chothia，记录工具及依赖版本。不要在不同样本间混用 IMGT/Kabat/Chothia 的位置或 CDR 边界。
- 保留完整 Chothia 位置，包括插入码，例如 `H52A`、`H100A`。它们都是独立序列位置，不能强转整数、合并、丢弃或按普通字符串排序。
- 保存 FR1、CDR1、FR2、CDR2、FR3、CDR3、FR4 的区域标注。AbNumber 的 raw 数值索引为 0 基，输出序号必须显式转换成 1 基。
- 编号失败、序列被意外截断、结构域边界不明确或非重链结果均进入错误/待核验表，不能强行补编号。

## 5. 建立不可丢失的残基映射

对 `variable_sequence` 的每个位置建立一行记录，包括没有坐标的位置。长度为 L 的可变区必须有 L 行。

| 字段 | 语义 |
| --- | --- |
| `sample_id` | 唯一标识结构、链、模型和结构域 |
| `full_seq_idx1` | 在完整提交序列中的位置，从 1 开始 |
| `seq_idx1` | 在所提取可变区内的位置，严格为 1…L |
| `aa` | 该位置的序列氨基酸 |
| `chothia_pos` | 完整 Chothia 位置，保留插入码 |
| `region` | Chothia 定义的 FR/CDR 区域 |
| `entity_id / label_asym_id / auth_asym_id` | 原始聚合物及链身份 |
| `label_seq_id` | 原始 mmCIF 的聚合物序列位置；不一定等于裁剪后的 `seq_idx1` |
| `auth_seq_id / insertion_code` | 作者残基编号及原始插入码；未知时留空，不猜测 |
| `model_id / original_resname` | 原始模型和残基身份 |
| `observed_idx1` | 可选：已出现坐标记录的聚合物残基次序，仅供追踪，不能用于序列索引 |
| `renumbered_res_id` | 本流程导出的结构编号，等于 `seq_idx1` |
| `coordinate_record_present` | 是否存在该残基的原始坐标记录 |
| `atom_mask / backbone_mask` | 四个主链原子的有效性及主链是否齐全 |
| `mapping_method / mapping_status` | 映射方法、是否唯一可靠、冲突原因 |

### 5.1 mmCIF 优先：使用文件内的显式对应关系

- 通过 entity、链和模型确定对象，用 `_atom_site.label_seq_id` 对接聚合物序列位置；通过 `_pdbx_poly_seq_scheme` 核对 label 编号、作者编号、链名和插入码。
- 原始坐标记录至少以 model、链、作者残基编号、插入码等联合定位，不能只用一个整数 residue ID 作为跨链或跨模型键。
- 原始作者插入码与 Chothia 插入码分别保存；两者不一定相同。
- 先建立完整链的序列—结构映射，再按可变区边界取子集，加入 Chothia 位置和 `seq_idx1`。
- 逐位检查序列残基与结构残基是否相容；有明确依据的修饰残基归一化可接受，但必须留下原始身份及转换记录。无法解释的 mismatch 不能当作正确匹配。

### 5.2 PDB 或已重编号结构：必要时回退到序列对齐

- 优先结合 SEQRES、原始结构、已有残基映射表及 REMARK 465/470 等信息建立关系。
- 若输入已经把可见残基连续 renumber 成 `1…N_observed`，这些数字不等于完整序列位置。首先寻找重编号前后的映射，不能直接据此配对。
- 没有显式映射时，允许用按链内顺序提取的 `observed_sequence` 对齐完整提交序列；支持内部缺失和端部缺失，保存算法、参数、alignment 以及逐位映射。
- 不能通过统一减去首残基编号、加一个 offset 或机械 zip 来处理内部 missing；offset 只适用于已证实没有内部缺口的特定编号转换。
- 同分对齐、多处重复片段、无法解释的结构侧插入或残基不匹配，应结合原始信息消歧；仍无法唯一定位就标记 `ambiguous_mapping`，不进入正式集合。
- 不要对缺失后的短 ATOM 序列重新跑 AbNumber，再把结果当成完整序列的 Chothia 编号。

## 6. 明确“以 1 为基的结构”是什么意思

**本流程主输出采用序列基准：完整可变区的第一个序列位置为 1，`renumbered_res_id = seq_idx1`。有坐标的残基使用其对应序列位置，缺失位置保留空位，不压缩后续编号。**

示例仅说明编号规则，不是实际样本：

| 可变区序列位置 | 坐标情况 | 导出结构 residue ID |
| --- | --- | --- |
| 29 | 有坐标 | 29 |
| 30 | 整个残基缺失 | 无坐标记录；映射表保留第 30 行 |
| 31 | 整个残基缺失 | 无坐标记录；映射表保留第 31 行 |
| 32 | 整个残基缺失 | 无坐标记录；映射表保留第 32 行 |
| 33 | 有坐标 | 33，不能变成 30 |

- 若第 1–3 位缺失，第一条可见残基记录可以编号为 4。这仍然是完整可变区从 1 开始的编号体系。
- 若下游工具必须要求可见残基连续编号，额外导出独立的兼容副本，并保存 `observed_idx1 ↔ seq_idx1 ↔ 原始结构编号` 映射；不能覆盖主输出或误称其为完整序列编号。
- 同时保存 Chothia 编号的对应关系；如项目需要 Chothia 编号结构，可额外导出。Chothia 视图与 1 基序列视图都从同一映射生成，不能各自推断。
- 不修改原始文件。重编号只改变标识，不改变原子坐标；导出后验证这一点。
- 若导出新的 mmCIF，保持聚合物序列、atom_site、链和编号相关类别一致；不能只改一个字段而留下相互矛盾的表。legacy PDB 表达不了的编号不能截断，应使用 mmCIF。

## 7. 区分三种“对齐”

### A. 同一样本的完整序列—结构映射（必须）

按第 5 节把每个有坐标的残基映射回正确序列位置。坐标缺失只改变 mask，不改变该位置的氨基酸、Chothia 位置或后续残基次序。

### B. 不同 VH/VHH 的序列位置对齐（数据集对齐）

按完整可变区的 Chothia 位置建立共同列，包括插入位置。可使用 AbNumber 的编号对齐能力，并保存公共列顺序与每条序列的映射。不要拿各结构当前 residue ID 直接当同源位置。

必须区分：

- `sequence_gap`：某条序列没有公共 Chothia 列对应的残基。
- `missing_residue`：该序列有残基，但没有该残基坐标记录。
- `missing_atom`：残基有坐标记录，但部分主链原子没有记录。
- `zero_occupancy`：原子有坐标记录，但占有率为零，按配置视为无有效观测，默认不纳入有效坐标。

若输出公共长度 M 的对齐数组，额外保存 `sequence_mask[M]` 和 `alignment_col ↔ seq_idx1` 映射。不能把 `sequence_gap` 与实验结构 missing 统一写成一种缺失原因。

### C. 三维结构叠合（仅在项目需要时）

序列对齐不自动等于坐标叠合。只有项目要求比较坐标/RMSD 时才执行叠合；使用明确的对应位置和两侧均有效的原子，记录参考结构、拟合区域和变换。默认可选择共同框架区 Cα 作为拟合集，CDR 单独评估；不要让缺失占位坐标参与拟合。

## 8. 最后提取主链，并显式保存缺失

映射、Chothia 对齐和 1 基编号确定后，再从目标结构域提取主链。固定原子顺序为 **`[N, CA, C, O]`**，`OXT` 不能替代 `O`。

每个样本至少输出：

```text
sequence:                 长度 L 的完整可变区序列
coords:                   float array [L, 4, 3]，单位 Å
atom_mask:                bool array  [L, 4]
backbone_mask:            bool array  [L]，等于 atom_mask.all(axis=1)
coordinate_record_present: bool array [L]
chothia_positions:        长度 L，含插入码
regions:                  长度 L
residue_mapping:          L 行
```

- 缺失或被判为无效的坐标在持久化数组中存为 `NaN`，对应 `atom_mask=False`；合法坐标可恰好等于 0，不能以数值是否为零推断存在性。
- 全残基无坐标记录：保留该氨基酸和映射行，四个主链 mask 均为 False，`coordinate_record_present=False`。
- 部分原子缺失：只关闭相应原子的 mask，其余坐标保留；不能因缺 O 就把整个残基从序列删除。
- 记录逐原子缺失/无效原因及原始 occupancy；结合 mmCIF 的未观测/零占有率类别或 PDB REMARK 信息交叉核验。注释缺失不等于坐标齐全，最终必须检查坐标记录本身。
- 若训练/运算接口需要有限数填充，只在接口层按 mask 临时转换；原始导出保留真实缺失信息，损失和几何计算必须使用 mask。
- 有多个模型时显式选择并记录 model ID，或把各模型独立处理；不得把不同模型的原子混在一个样本中。
- altloc 采用可复现的残基层面选择规则：兼容空白共用原子，选择一个一致的构象标签，并记录选择理由。不能逐原子拼出不存在的混合构象。默认可优先主链齐全且有效占有率较高的构象，平分时使用固定顺序。
- 聚合物中的修饰残基即使使用 HETATM 记录，也应按聚合物身份判断是否保留；不能把所有 HETATM 一律删除。水、配体、抗原链不进入目标主链数组。
- 主链结构文件仅写实际保留的有效原子；不写虚假的缺失占位原子。完整序列和缺失必须由配套 FASTA、映射表和 mask 表达。

## 9. 输出与可复现性

尊重项目已有目录结构；若尚无约定，按以下类型组织文件：

| 产物 | 最低内容 |
| --- | --- |
| 原始数据 | 未改动的 mmCIF/PDB、完整 FASTA、来源元数据、校验和 |
| 检索与配置记录 | 查询条件、筛选参数、版本、运行命令、下载时间和失败重试记录 |
| `metadata.tsv` | 每个候选一行，结构/链/域身份、序列来源、类型证据、质量统计及状态 |
| `full_chain.fasta` / `variable.fasta` | 完整提交序列与提取的可变区序列，header 能对应 sample ID |
| `residue_mapping.tsv` | 每个可变区序列位置一行，包含原始编号、Chothia 和 1 基映射 |
| `missing_atoms.tsv` | 缺失/无效位置、原子名、原因、原始注释及 occupancy |
| `backbone_1based.cif` 或 `.pdb` | 按可变区序列位置编号的主链坐标 |
| `backbone.npz` 或项目约定格式 | 固定轴顺序的坐标、mask、序列和位置标签 |
| Chothia 对齐产物 | 对齐序列、公共位置列、公共列到各序列位置的映射 |
| `failed_or_review.tsv` | 阶段、错误类型、具体原因、样本身份和已尝试的修复 |
| `qc_report.md` | 各阶段数量、纳入/排除理由、缺失统计及核验结果 |

不同样本使用独立目录或明确的 sample ID，避免同名覆盖。原始链 ID、导出链 ID 均记录在映射中。下载采用缓存、合理重试和断点续跑；一次失败不能导致整批静默中止或样本无记录地消失。

## 10. 验收条件与执行顺序

先完成少量真实样本的端到端试运行，核验后再按配置执行批量处理。不要把试运行样本数当成用户最终要求的数据规模。

必须满足的检查：

1. `L == len(variable_sequence) == 映射表行数 == coords.shape[0] == atom_mask.shape[0]`。
2. `seq_idx1` 严格为 1…L，无重复；可变区序列能够按已记录边界还原到完整提交序列。
3. 每个正式纳入的坐标残基映射唯一；插入码、原始链名、model 和域身份没有丢失。
4. 每个有效坐标残基与其序列身份一致，或具有明确记录的化学修饰解释；未解释 mismatch/歧义样本隔离。
5. 重编号前后，同一保留原子的坐标不变；内部或端部缺失不会令后续 `seq_idx1` 偏移。
6. `atom_mask=True` 的坐标全部有限；缺失位置不会参与几何计算。主链 mask 与四原子 mask 一致。
7. Chothia 对齐去掉 sequence gaps 后能恢复各样本的可变区序列，结构 missing 不会被删除成 sequence gap。
8. 每个检索候选都有最终状态；报告发现数、下载成功数、重链域识别数、映射成功数、正式纳入数及各类失败数。

为核心映射逻辑设置有针对性的回归样例：N 端缺失、内部连续 missing、只缺一个主链原子、原始或 Chothia 插入码、已有结构从 1 连续压缩编号、重复片段造成对齐歧义。可以使用明确标为合成的最小样例验证逻辑，但不能冒充真实数据结果。

完成后说明实际使用的数据来源、检索范围、样本数量、可复现运行命令、输出路径和未解决问题。没有执行的下载、验证或处理不得宣称已完成；网络或依赖失败应报告真实原因，并保留已完成结果。

## 官方参考（执行时核对当前接口与版本）

- RCSB 文件与完整序列下载：https://www.rcsb.org/docs/programmatic-access/file-download-services
- RCSB 链与残基标识：https://www.rcsb.org/docs/general-help/identifiers-in-pdb
- PDB-101 序列与缺失坐标说明：https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/primary-sequences
- wwPDB mmCIF 用户指南：https://mmcif.wwpdb.org/docs/user-guide/guide.html
- AbNumber 官方文档：https://abnumber.readthedocs.io/en/stable/
- AbNumber 官方代码：https://github.com/prihoda/AbNumber
- SAbDab 候选与注释入口：https://sabdab.opig.stats.ox.ac.uk/
