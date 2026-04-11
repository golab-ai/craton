# Craton 命令行参考

本文档描述 Craton 命令行入口 `craton` 下各子命令的用法。安装后可通过 `craton --help` 查看顶层命令列表，通过 `craton <命令> --help` 查看该命令的子命令与选项。

---

## 1. 总览

### 1.1 调用方式

```bash
craton [OPTIONS] COMMAND [ARGS]...
```

全局选项：`-h` / `--help` 显示帮助。

### 1.2 顶层命令速查

| 命令 | 说明 |
|------|------|
| `simulation` | 提交 MD / QM 模拟（含 FEP、溶液、蛋白等）；支持 YAML 任务与 autoqm |
| `prepare` | 分子与蛋白准备：分子信息、配体质子化、蛋白预处理、残基突变与修饰、UniProt/PDB 获取 |
| `md_check` | 检查 MD 作业运行状态 |
| `analyze` | 分析 MD 结果与自由能（BAR/ddG/RBFE 等）、性质与 ADMET、化学空间、训练/测试集划分、环闭合约简、FEP 全对检查 |
| `data` | 向/从数据库写入或读取化合物、QM 数据、项目虚拟化合物（prj_vcompound） |
| `ff` | 力场：原子类型、赋力场、AM1-BCC 电荷、拟合/验证/读参（含 yaml、autofit） |
| `mm` | 基于力场的计算：能量/力/Hessian/频率、优化、扭转扫描、多极矩、体积/表面积、质心、转动惯量 |
| `stru` | 分子拓扑与结构分析：环/手性/扭转/片段/官能团/粗粒化/作用位点等、结构参数测量与变更、RMSD |
| `tool` | 工具集：分子描述保存、Gaussian 扭转分析、化学信息、PubChem、粗粒化、性质查询、肽/RNA/DNA 生成、力场格式转换（atf/ff/rtp/Amber） |

---

## 2. simulation — 模拟提交

用于准备并提交 GROMACS MD 或 QM 计算，大部分默认设置在 `configure/configure.yaml` 中；需更细控制时可使用 YAML 任务文件。

### 2.1 基本使用方式

```bash
craton simulation [SIMULATION_TYPE] [OPTIONS]
```

模拟类型 (SIMULATION_TYPE) 包含如下选项：

| 类型 | 说明 | 示例 |
|------|------|------|
| `vacuum` | 气相 | `craton simulation vacuum --molecules mols.sdf -o out` |
| `solution` | 溶液相 | `craton simulation solution --molecules . -o out` |
| `liquid` | 纯液体 | `craton simulation liquid --molecules liquid.sdf -o out` |
| `complex` | 蛋白-配体复合物 | `craton simulation complex --protein prot.pdb --ligands lig.sdf -o out` |
| `protein` | 蛋白 | `craton simulation protein --protein prot.pdb -o out` |
| `rbfe` | 相对结合自由能 (FEP) | `craton simulation rbfe --protein prot.pdb --ligands series.sdf -o out` |
| `abfe` | 绝对结合自由能 | `craton simulation abfe --protein prot.pdb --ligands lig.sdf -o out` |
| `rhfe` / `ahfe` | 水合自由能（相对/绝对） | `craton simulation rhfe --molecules mols.sdf -o out` |
| `alogp` / `rlogp` | 分配系数（绝对/相对） | `craton simulation alogp --molecules mol.sdf -o out` |
| `alogs` / `rlogs` | 溶解度（绝对/相对） | `craton simulation alogs --molecules mol.sdf -o out` |
| `yaml` | 完全由 YAML 定义任务 | `craton simulation yaml -f task.yaml` |
| `Q0`–`Q10`, `Q100` | QM 阶段（构象搜索等） | `craton simulation Q0 --molecules . -o out` |
| `autoqm` | 自动 QM 计算流程 | `craton simulation autoqm -n 10` |

<!-- **仅通过 YAML 指定 `simulation_type` 时可用**（命令行无对应选项）：`hfe`、`mem-rbfe`、`cov-rbfe`、`mutation`、`pep-rbfe`、`rna-rbfe`、`bilayer`、`biomembrane`、`mem-protein`。 -->

基本的命令行选项如下：

| 选项 | 简写 | 默认 | 说明 |
|------|------|------|------|
| `--ligands` | — | — | 配体文件，如 `.sdf` |
| `--protein` | — | — | 蛋白文件，如 `.pdb` |
| `--molecules` | — | — | 分子（SMILES 或文件） |
| `--coligands` | — | — | 共配体文件 |
| `--repeat` | — | 1 | 重复计算次数 |
| `--output` | `-o` | output | 输出目录 |
| `--simulation_time` | `-t` | — | 模拟时间（与步长配合） |
| `--charge_method` | `-c` | — | 电荷方法 |
| `--molecule_number` | `-n` | — | autoqm 的分子数量 |
| `--yaml_file` | `-f` | — | 细控用 YAML 任务文件 |

### 2.4 使用 YAML 提交作业示例

通过 YAML 可指定 `simulation_type`、蛋白、配体、共配体等，适合复现或批量任务。下面示例：先查看 `task.yaml` 内容，再提交带共配体 (coligands) 的作业。

```bash
# 查看任务 YAML 内容
cat task.yaml
```

`task.yaml` 示例（蛋白-配体-共配体复合物）：

```yaml
simulation_type: complex
protein: protein.pdb
ligands: ligands.sdf
coligands: cofactor.sdf
output_directory: ./output_complex
```

```bash
# 使用该 YAML 提交作业
craton simulation yaml -f task.yaml
```

多任务时可在同一 YAML 中用 `task0`、`task1` 等键定义多个任务，例如：

```yaml
task0:
  simulation_type: complex
  protein: protein.pdb
  ligands: ligands.sdf
  coligands: cofactor.sdf
  output_directory: ./output_0
task1:
  simulation_type: rbfe
  protein: protein.pdb
  ligands: series.sdf
  output_directory: ./output_rbfe
```

然后同样执行 `craton simulation yaml -f task.yaml` 即可。

### 2.5 MD 模拟

以下仅针对 **MD 类** 的 `simulation_type`。每种类型对应一种模拟流程，并依赖不同的输入选项；未列出的选项（如 `--output`、`--simulation_time`）为通用可选参数。

| SIMULATION_TYPE | 模拟内容 | 必须提供的选项 |
|-----------------|----------|----------------|
| `vacuum` | 气相中的分子模拟：构建单分子或少量分子体系，无溶剂。 | `--molecules`（SMILES 或分子文件） |
| `solution` | 溶液相模拟：溶质置于溶剂盒子中，进行溶剂化 MD。 | `--molecules`（溶质分子） |
| `liquid` | 纯液体模拟：多分子堆积成液相盒子。 | `--molecules`（液体分子） |
| `complex` | 蛋白-配体复合物常规 MD：蛋白与配体一起建盒、跑 MD。 | `--protein`，`--ligands`；可选 `--coligands` |
| `protein` | 仅蛋白的 MD：无配体，仅蛋白建盒与模拟。 | `--molecules` (pdb文件传递给`--molecules`参数即可，勿传递给`--protein`) |
| `rbfe` | 相对结合自由能 (FEP)：双拓扑、配体对，在蛋白环境下计算 ΔΔG。 | `--protein`，`--ligands`（含多个配体的文件或目录） |
| `abfe` | 绝对结合自由能：单配体在蛋白结合位点的结合自由能。 | `--protein`，`--ligands` |
| `rhfe` | 相对水合自由能：多个分子之间的相对水合自由能差。 | `--molecules`（多个分子） |
| `ahfe` | 绝对水合自由能：单分子的水合自由能。 | `--molecules` |
| `alogp` | 绝对分配系数 (logP)：单分子在溶剂相间的分配。 | `--molecules`（或依配置提供配体） |
| `rlogp` | 相对分配系数：分子间相对 logP。 | `--molecules` 或 `--protein` + `--ligands`（依流程） |
| `alogs` | 绝对溶解度：单分子溶解度相关计算。 | `--molecules`（或依配置） |
| `rlogs` | 相对溶解度：分子间相对溶解度。 | `--molecules` 或 `--protein` + `--ligands`（依流程） |

**仅通过 YAML 指定 `simulation_type` 时可用**（命令行无对应子命令）：

| SIMULATION_TYPE | 模拟内容 | 必须提供的配置/选项（在 YAML 中） |
|-----------------|----------|-----------------------------------|
| `hfe` | 水合自由能（通用，同 ahfe/rhfe 流程） | `MoleculeFileSetting.molecules` 或相应分子输入 |
| `mem-rbfe` | 膜环境相对结合自由能 | 蛋白/膜/配体相关路径（见配置说明） |
| `cov-rbfe` | 共价结合相对自由能 | 蛋白、配体及共价信息 |
| `mutation` | 蛋白突变扫描（如点突变 FEP） | 蛋白、突变列表（如 `AlignmentSetting.mutation`） |
| `pep-rbfe` | 肽段相对结合自由能 | 蛋白（肽段）及序列/突变设置 |
| `rna-rbfe` | RNA 相关相对自由能 | 蛋白/RNA、配体设置 |
| `bilayer` | 双层膜体系 | 蛋白/膜、分子文件路径 |
| `biomembrane` | 生物膜体系 | 同上 |
| `mem-protein` | 膜蛋白体系 | 蛋白、膜、分子文件路径 |

以上 YAML 类型的必填项以 `configure/configure.yaml` 及任务 YAML 中的 `MoleculeFileSetting`、`AlignmentSetting` 等为准。

### 2.6 QM 作业（Q0–Q100）

通过 `craton simulation Q0`、`Q1`、…、`Q10`、`Q100` 提交的是 **QM 构象搜索与单点计算** 流程，用于力场拟合、构型采样或高精度能量计算。各阶段含义与配置均可在 `configure/configure.yaml` 中查看与修改。

#### 阶段含义（与代码/配置对应）

各阶段对应的 **qmjobs** 在 `craton/craton/software/gaussian.py` 中会写成 Gaussian 路由关键字（`# 方法 基组 关键字`），下表给出关键字部分（方法/基组由 configure 的 `opt_method_basisset` 或 `sp_method_basisset` 指定）。

| 阶段 | 说明 | Gaussian 关键字 |
|------|------|-----------------|
| `Q0` | 初始结构优化、频率、电荷计算 | `opt freq pop=chelpg`（或配置的 charge_model） |
| `Q1` | 扭转扫描（一维） | `opt=addred` |
| `Q2` | 键长/键角扫描，单点能 | `force` |
| `Q3` | 二维扭转扫描 | 见代码/扩展配置 |
| `Q4` | 扭转扫描得到的局部极小优化 | `opt freq pop=chelpg` |
| `Q5` | 随机构象搜索 | 见代码/扩展配置 |
| `Q6` | 寡聚体/分子间（intermolecular） | `opt freq pop=chelpg` |
| `Q7` | 过渡态 | 见代码/扩展配置 |
| `Q8` | 扫描（含 inter_val、可忽略环等） | `opt=addred` |
| `Q10` | 单点与频率/电荷（高精度方法） | `freq pop=chelpg` |
| `Q100` | CCSD(T)/CBS 等高精度单点 | 见代码/扩展配置 |


未在表中单独列出的阶段（如 Q3、Q5、Q7、Q100）若在命令行或 YAML 中使用，其行为以代码实现为准。

#### 输入与提交方式

- **分子来源**：通常由 `MoleculeFileSetting.molecules` 或 YAML 中的 `molecules` 指定（文件或目录）；Q6 阶段可从数据库读取（`get_from_db`）。  
- **命令行示例**：  
  `craton simulation Q0 --molecules . -o out`  
  `craton simulation Q1 --molecules ./mol_dir -o out`  
- **YAML 示例**：在任务 YAML 中设置 `simulation_type: Q0`（或 Q1、…、Q10、Q100）及 `molecules`、`output_directory` 等，再执行：  
  `craton simulation yaml -f task.yaml`

---

## 3. md_check — MD 作业状态检查

检查指定目录下 MD 作业的完成状态（基于 GROMACS 等引擎）。

### 3.1 用法

```bash
craton md_check [OPTIONS]
```

### 3.2 选项

| 选项 | 简写 | 默认 | 说明 |
|------|------|------|------|
| `--dir` | `-d` | . | 运行目录 |
| `--batchfile` | `-f` | — | 匹配作业 ID 的批处理文件 |
| `--mdengine` | `-e` | gmx | MD 引擎 |

### 3.3 示例

```bash
craton md_check -d ./rbfe_output
craton md_check -d ./output -f submit.sh
```

---

## 4. analyze — 分析

包含 GMX 轨迹/自由能分析、性质统计、化学空间、数据集划分等。

### 4.1 子命令速查

| 子命令 | 说明 |
|--------|------|
| `gmx` | 基于 GROMACS 结果的分析（BAR、ddG、RMSD、能量、扭转、相互作用、FEP 等） |
| `get-property` | 按温度等条件获取/整理实验性质 |
| `property-bin` | 性质分布分箱统计 |
| `admet-bin` | ADMET 数据分箱 |
| `property-info` | 在性质文件中加入分子信息（如重原子数、环等） |
| `property-block` | 将指定分子的实验数据标记为无效 |
| `property-result` | 性质预测结果汇总 |
| `admet-result` | ADMET 预测结果汇总 |
| `property-figure` | 重新生成性质相关图 |
| `property-script` | 生成微调用脚本 |
| `chem-space` | 化学空间分布分析 |
| `train-test` | 训练/测试集划分（random/arom/ring/halogen/element 等） |
| `property-add` | 向性质文件追加新实验数据 |
| `cc` | 环闭合约简（cycle closure）重算 |
| `all_check` | 检查 FEP 所有 pair 状态并可选生成重启脚本 |

### 4.2 analyze gmx — GMX 分析类型

```bash
craton analyze gmx <analyze_type> [OPTIONS]
```

**analyze_type** 可选：`bar`, `ddg`, `rmsd`, `energy`, `torsion`, `interaction`, `pair_interaction`/`pair-interaction`, `fep_exchange`/`fep-exchange`, `block_ddg`/`block-ddg`, `accum_ddg`/`accum-ddg`, `rbfe`, `abfe`, `complex`, `protein`, `normal`, `all-rbfe`/`all_rbfe`, `all-abfe`/`all_abfe`, `dimer`/`all_dimer`/`all-dimer`, `total-dimer`/`total_dimer`。

| 选项 | 简写 | 默认 | 说明 |
|------|------|------|------|
| `--input_directory` | `-i` | . | 任务/轨迹所在目录 |
| `--output_directory` | `-o` | md_result | 分析结果输出目录 |
| `--expt_file` | `-exp` | — | 实验 dG 文件（csv 或 gpickle） |
| `--two_stages` / `--no-two_stages` | — | False | 是否两阶段 FEP |
| `--pka-file` | `-pka` | — | pKa 文件 |
| `--molecule_dir` | `-mol` | — | 分子文件目录 |
| `--attributes` | `-attrs` | — | 分析属性（如 energy 的 Potential,Temperature,Pressure） |
| `--dimer_type` | — | — | 二聚体类型等 |

### 4.3 示例

```bash
# BAR 自由能
craton analyze gmx bar -i ./fep_pairs -o ./bar_result

# 单对 RBFE 分析
craton analyze gmx rbfe -i ./pair_A_to_B -o ./rbfe_analysis

# 全 RBFE 任务汇总
craton analyze gmx all-rbfe -i ./rbfe_task -o ./all_rbfe_result -exp ref_dG.csv
```

---

## 5. data — 数据库读写

向 MongoDB 写入或从 MongoDB 读取化合物、QM 数据等。

### 5.1 子命令

| 子命令 | 说明 |
|--------|------|
| `insert` | 插入数据到库（compound / qmdata / prj_vcompound） |
| `get` | 按条件从库中查询并导出（如 result.json） |

### 5.2 data insert

```bash
craton data insert [compound|qmdata|prj_vcompound] [OPTIONS]
```

| 选项 | 简写 | 默认 | 说明 |
|------|------|------|------|
| `--inputs` | `-i` | . | 输入文件或目录 |
| `--molecule_type` | `-t` | molecule | 分子类型（small molecule, peptide, protein 等） |
| `--element_flag` | `-elem` | F | 是否为每条记录写元素符号 (T/F) |

### 5.3 data get

```bash
craton data get [compound|qmdata] [OPTIONS]
```

| 选项 | 简写 | 说明 |
|------|------|------|
| `--inputs` | `-i` | 查询输入/条件 |
| `--yaml_files` | `-f` | 查询选择器 YAML（qmdata 时） |

### 5.4 示例

```bash
craton data insert compound -i ./molecules.sdf -t molecule
craton data get qmdata -i selector.yaml -f query.yaml
```

---

## 6. ff — 力场

原子类型划分、力场赋值、AM1-BCC 电荷、力场拟合与验证。

### 6.1 子命令速查

| 子命令 | 说明 |
|--------|------|
| `atom_type` | 对分子做原子类型划分 |
| `assign_ff` | 为分子分配力场参数 |
| `am1bcc` | 计算 AM1-BCC 电荷 |
| `fit` | 力场拟合/验证/读参（fitting, validation, read_parameter, yaml, autofit） |

### 6.2 选项（atom_type / assign_ff / am1bcc）

| 选项 | 简写 | 默认 | 说明 |
|------|------|------|------|
| `--input_files` | `-i` | . | 分子文件（SMILES, *.mol, *.sdf, *.pdb 等） |
| `--atom_type_file` | `-f`（atom_type）/ `-af`（assign_ff） | — | 原子类型定义文件 |
| `--force_field_file` | `-f`（assign_ff） | — | 力场参数文件 |
| `--output_directory` | `-o` | . | 输出目录 |

### 6.3 ff fit

```bash
craton ff fit [fitting|validation|read_parameter|yaml|autofit] [OPTIONS]
```

| 选项 | 简写 | 默认 | 说明 |
|------|------|------|------|
| `--input_files` | `-i` | . | 拟合分子或目录 |
| `--output_directory` | `-o` | . | 输出目录 |
| `--init_force_field_file` | `-if` | ./0.ff | 初始力场文件 |
| `--yaml_file` | `-f` | — | 细控 YAML（fit_type=yaml 时必填） |

YAML 可配置：`fitting_molecules`, `atom_type_file`, `force_field_file`, `qm_data_path`, `charge_method`, `fitting_terms`, `target_prop`, `out_put_dir` 等。

### 6.4 示例

```bash
craton ff atom_type -i ligands.sdf -o ./atom_types
craton ff assign_ff -i ligands.sdf -f params.ff -af types.atf -o ./assigned
craton ff am1bcc -i ligands.sdf -o ./charges
craton ff fit fitting -i ./train_mols -if 0.ff -o ./fit_result
craton ff fit yaml -f fit_config.yaml
```

---

## 7. mm — 力场计算

基于已赋力场的分子做单点、优化、扫描与性质计算。

### 7.1 子命令速查

| 子命令 | 说明 |
|--------|------|
| `calculate` | 能量、力、Hessian、频率 |
| `opt` | 结构优化 |
| `scan` | 扭转扫描 |
| `multipole` | 偶极、四极、八极等多极矩 |
| `volume` | 体积相关 |
| `surface` | 表面积相关 |
| `center` | 几何中心（cog/com/cob/size 等） |
| `inertia` | 转动惯量 |

### 7.2 选项（通用）

| 选项 | 简写 | 默认 | 说明 |
|------|------|------|------|
| `--inputs` | `-i` | . | 分子文件或目录 |
| `--output_dir` | `-o` | . | 输出目录 |
| `--prop` | `-p` | （见各子命令） | 计算类型或属性 |

- **calculate**：`prop` 取 `energy`, `force`, `hessian`, `freq`/`frequency`。
- **multipole**：`prop` 取 `dipole`, `quadrupole`, `octupole`, `multipole`，默认 `energy`。
- **center**：`prop` 取 `center`, `cog`, `com`, `cob`, `size`，默认 `cog`。

### 7.3 示例

```bash
craton mm calculate energy -i ./mols -o ./energies
craton mm opt -i initial.sdf -o ./optimized
craton mm scan -i mol.sdf -o ./torsion_scan
craton mm multipole -i mol.sdf -p dipole -o ./multipole
craton mm center -i mol.sdf -p cog -o ./centers
```

---

## 8. stru — 结构与拓扑分析

分析分子拓扑（环、手性、扭转、片段、作用位点等）、测量结构参数、改变结构参数、计算 RMSD。

### 8.1 子命令速查

| 子命令 | 说明 |
|--------|------|
| `topol` | 拓扑分析（类型见下） |
| `measure` | 键长/键角/二面角等结构参数测量 |
| `vary` | 改变结构参数（键长、角、二面角等） |
| `rmsd` | 两分子间 RMSD |

### 8.2 stru topol — 拓扑类型

```bash
craton stru topol [ring|chiral|torsion|hybrid|interaction_site|frag|cg|fg|image] [OPTIONS]
```

| 类型 | 说明 |
|------|------|
| `ring` | 环分析 |
| `chiral` | 手性中心 |
| `torsion` | 可旋转键/二面角 |
| `hybrid` | 杂化等 |
| `interaction_site` | 相互作用位点 |
| `frag` | 片段化 (fragmentation) |
| `cg` | 粗粒化 bead (atom_cluster) |
| `fg` | 官能团 (function_group) |
| `image` | 分子平面图 |

### 8.3 选项（通用）

| 选项 | 简写 | 默认 | 说明 |
|------|------|------|------|
| `--inputs` | `-i` | . | 输入文件或目录 |
| `--output_dir` | `-o` | . | 输出目录 |
| `--atoms` | `-a` | — | 原子 ID，如 `2-3`、`4-5-6`（measure/vary 必填） |
| `--value` | `-v` | — | vary 时的目标值 |
| `--del_value` | `-dv` | False | 是否删除该参数约束 |
| `--inputs2` | `-t` | — | rmsd 时的目标分子（参考为 `-i`） |

### 8.4 示例

```bash
craton stru topol ring -i molecules.sdf -o ./ring_info
craton stru topol torsion -i ligand.sdf -o ./torsions
craton stru measure -i mol.sdf -a 2-3 -o ./bond_length   # 键长
craton stru measure -i mol.sdf -a 4-5-6 -o ./angle        # 键角
craton stru vary -i mol.sdf -a 3-7-9-12 -v 90 -o ./rotated
craton stru rmsd -i ref.pdb -t target.pdb -o ./rmsd_out
```

---

## 9. tool — 工具集

分子描述保存、高斯扭转分析、PubChem 查询、粗粒化、肽/核酸生成、力场格式转换等。

### 9.1 子命令速查

| 子命令 | 说明 |
|--------|------|
| `save` | 生成并保存分子描述文件 |
| `torsion` | 基于 Gaussian log 的扭转扫描结果分析 |
| `chem-info` | 化学信息输出 |
| `pubchem` | PubChem 查询（smiles/name/cas_no 等） |
| `cg` | 粗粒化并输出 mtx 等 |
| `property` | 按温度等计算/查询性质（如 density） |
| `peptide` | 生成肽段（残基数、N/C 端帽、模板） |
| `rna` / `dna` | 生成 RNA/DNA 片段 |
| `atf2json` | 原子类型定义文件 → JSON |
| `fff2json` | 力场参数文件 → JSON |
| `aminoacid` | rtp → 氨基酸 JSON |
| `amberff` | Amber 力场 → 内部 ff 格式 |

### 9.2 常用选项与示例

**save**

```bash
craton tool save -i ./mols -o ./descriptors
```

**torsion**（Gaussian 扭转扫描分析）

```bash
craton tool torsion -i ./gaussian_logs -o ./torsion_curves
```

**pubchem**

```bash
# 输入类型：smiles, name, cas_no 等；输出类型：name 等
craton tool pubchem -i "CC(=O)OC1=CC=CC=C1C(=O)O" -it smiles -ot name
```

**property**

```bash
craton tool property -i "514-10-3" -p density -t 300.0
# 可用 -f 指定 YAML 配置
```

**peptide**

```bash
# -n 残基数，-r C 端帽，-l N 端帽，-t 模板
craton tool peptide -n 5 -r NME -l ACE -o ./pep_out
```

**rna / dna**

```bash
craton tool rna -n 3 -t rna -o ./rna_out
craton tool dna -n 3 -t dna -o ./dna_out
```

**格式转换**

```bash
craton tool atf2json -f atom_types.txt
craton tool fff2json -f forcefield.ff
craton tool aminoacid -f aminoacids.rtp
craton tool amberff -a types.atf -v nonbond.ff -b bond.ff
```

---

## 10. prepare — 分子与蛋白准备

小分子质子化、蛋白预处理、点突变、修饰、UniProt/PDB 数据获取。

### 10.1 子命令速查

| 子命令 | 说明 |
|--------|------|
| `mol_info` | 分子信息（IUPAC、CAS、SMILES 等，如通过 PubChem） |
| `ligand` | 小分子准备（如 pH 下质子化） |
| `uniprot` | UniProt 信息（序列、PDB、fasta） |
| `pdb` | 从 RCSB 下载 PDB |
| `protein` | 蛋白预处理（质子化、端基、加氢等） |
| `mutation` | 残基突变 |
| `modify` | 残基修饰（pho/suf/met/n-met） |

### 10.2 选项摘要

**mol_info**

| 选项 | 简写 | 默认 | 说明 |
|------|------|------|------|
| `--inputs` | `-i` | . | 分子输入 |
| `--in_type` | `-it` | smiles | 输入类型（name, smiles, cas-no 等） |
| `--output_dir` | `-o` | . | 输出目录 |

**ligand**

| 选项 | 简写 | 默认 | 说明 |
|------|------|------|------|
| `--inputs` | `-i` | — | 输入分子文件 |
| `--ph_min` | `-pi` | 7.4 | pH 下限 |
| `--ph_max` | `-pa` | 7.4 | pH 上限 |
| `--output_file` | `-of` | — | 输出文件名 |
| `--output_dir` | `-o` | . | 输出目录 |

**uniprot**

| 选项 | 简写 | 说明 |
|------|------|------|
| `--target` | `-t` | 目标蛋白名称 |
| `--uniprot_id` | `-i` | UniProt ID |
| `--output_dir` | `-o` | 输出目录 |

**pdb**

| 选项 | 简写 | 说明 |
|------|------|------|
| `--info_file` | `-i` | 含 PDB ID 列表的文件 |
| `--pdb_id` | `-id` | 单个 PDB ID |
| `--output_dir` | `-o` | 输出目录 |

**protein**

| 选项 | 简写 | 说明 |
|------|------|------|
| `--protein_file` | `-i` | 输入 PDB |
| `--output_file` | `-of` | 输出 PDB 文件名 |
| `--output_dir` | `-o` | 输出目录 |

**mutation / modify**

| 选项 | 简写 | 说明 |
|------|------|------|
| `--protein_file` | `-i` | 输入 PDB |
| `--residue` | `-r` | 目标残基，如 127-ARG-H |
| `--mutation` | `-m` | 突变目标（如 LEU）或修饰类型（pho/suf/met/n-met） |
| `--output_file` | `-of` | 输出文件名 |
| `--output_dir` | `-o` | 输出目录 |

### 10.3 示例

```bash
# 小分子在 pH 7.4 下准备
craton prepare ligand -i ligands.sdf -pi 7.4 -pa 7.4 -o ./prepared_ligands

# 从 UniProt 获取信息
craton prepare uniprot -t "MyProtein" -i P12345 -o ./uniprot_out

# 下载 PDB
craton prepare pdb -id 1ABC -o ./pdb_download

# 蛋白预处理
craton prepare protein -i raw.pdb -of prepared.pdb -o ./

# 残基突变
craton prepare mutation -i protein.pdb -r 127-ARG-H -m LEU -o ./mutant

# 残基修饰（如磷酸化 pho）
craton prepare modify -i protein.pdb -r 127-ARG-H -m pho -o ./modified
```

---

## 11. 使用建议

1. **优先查帮助**：`craton <command> --help`、`craton <command> <subcommand> --help` 可看到最新选项与默认值。
2. **细控用 YAML**：`simulation`、`ff fit` 等支持 `-f <file>.yaml`，适合复杂或多任务场景；示例可参考 `testing/` 下 YAML。
3. **路径**：未写绝对路径时，输入多为当前目录或 `configure` 中的 path；输出多用 `-o`/`--output_dir` 指定。
4. **数据库**：`data` 依赖 `configure/configure.yaml` 中 `mongodb` 配置，使用前确保连接与库名/集合正确。

若某子命令在本文未完全展开，请以 `craton <command> <subcommand> --help` 输出为准。
