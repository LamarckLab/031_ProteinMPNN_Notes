## Lamarck &nbsp; &nbsp; &nbsp; 2025-10-30
#### 该文档用于记录 server 上跑 ProteinMPNN 的命令
---

*环境 & 路径*
```bash
236 server上的环境: lmk_ProteinMPNN
236 server上的路径: /data/lmk/ProteinMPNN/protein_mpnn_run.py
```

*输入输出路径*
```bash
输入目录: /data/lmk/ProteinMPNN/mpnn_inputs    # 存放每次运行用的 PDB
输出目录: /data/lmk/ProteinMPNN/mpnn_outputs   # 解析得到的 jsonl 和生成的 FASTA
```

*GPU选择*
```bash
export CUDA_DEVICE_ORDER=PCI_BUS_ID
export CUDA_VISIBLE_DEVICES=3
```

> [!NOTE]
> 除「单 PDB 直接输入」(18) 外，所有功能都要先用 `parse_multiple_chains.py` 把一个文件夹里的 PDB 解析成 `parsed_pdbs.jsonl`，再喂给 `protein_mpnn_run.py --jsonl_path`。固定链 / 固定位点 / 绑定位点 / 偏置等额外约束都是先用对应 helper 脚本生成一个 `.jsonl`，再用相应 flag 传进去。所有命令统一用 `mpnn_inputs` / `mpnn_outputs` 一对固定目录，换功能时只需替换 `mpnn_inputs` 里的 PDB 再复跑（同名产物会被覆盖）。下面每个例子都在统一基线 `--num_seq_per_target 5 --sampling_temp "0.1" --seed 37 --batch_size 1` 上只改它要演示的那一项。

---

### 一、核心设计流程

> **01 简单单体设计**

**输入**：`mpnn_inputs/` 放单体 PDB（如 5L33.pdb、6MRR.pdb）  
**功能**：不加任何约束，固定骨架对整条链从头设计序列，每个骨架输出 5 条候选
```bash
python /data/lmk/ProteinMPNN/helper_scripts/parse_multiple_chains.py \
  --input_path=/data/lmk/ProteinMPNN/mpnn_inputs \
  --output_path=/data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl

python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
  --jsonl_path /data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --out_folder /data/lmk/ProteinMPNN/mpnn_outputs \
  --num_seq_per_target 5 \
  --sampling_temp "0.1" \
  --seed 37 \
  --batch_size 1
```

> **02 多链复合物 -- 指定设计链**

**输入**：`mpnn_inputs/` 放复合物 PDB（如 3HTN.pdb）  
**功能**：只重新设计多链复合物里的 A、B 链，其余链保持原序列作为固定环境

加 `assign_fixed_chains.py` 生成 `assigned_pdbs.jsonl` 标记哪些链参与设计，run 时用 `--chain_id_jsonl` 传入。下例只设计 A、B 链，其余链保持原序列
```bash
python /data/lmk/ProteinMPNN/helper_scripts/parse_multiple_chains.py \
  --input_path=/data/lmk/ProteinMPNN/mpnn_inputs \
  --output_path=/data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl

python /data/lmk/ProteinMPNN/helper_scripts/assign_fixed_chains.py \
  --input_path=/data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --output_path=/data/lmk/ProteinMPNN/mpnn_outputs/assigned_pdbs.jsonl \
  --chain_list "A B"

python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
  --jsonl_path /data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --chain_id_jsonl /data/lmk/ProteinMPNN/mpnn_outputs/assigned_pdbs.jsonl \
  --out_folder /data/lmk/ProteinMPNN/mpnn_outputs \
  --num_seq_per_target 5 \
  --sampling_temp "0.1" \
  --seed 37 \
  --batch_size 1
```

> **03 固定特定残基 -- 部分残基不参与设计**

**输入**：`mpnn_inputs/` 放复合物 PDB（如 3HTN.pdb）  
**功能**：在设计链上锁住指定残基（如关键活性位点），其余位置重新设计

`make_fixed_positions_dict.py` 生成 `fixed_pdbs.jsonl`，`--position_list` 中链间用逗号分隔、链内位点用空格。下例中 A 链锁 1-8/23/25，C 链锁 10-20/40

> [!IMPORTANT]
> `--position_list` 用的是 `parse_multiple_chains.py` 解析后**链内从 1 开始的顺序编号**（第几个残基），**不是 PDB 里的作者残基号**。链首编号 ≠ 1 或中间有缺失残基时，两者会错位。
```bash
python /data/lmk/ProteinMPNN/helper_scripts/parse_multiple_chains.py \
  --input_path=/data/lmk/ProteinMPNN/mpnn_inputs \
  --output_path=/data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl

python /data/lmk/ProteinMPNN/helper_scripts/assign_fixed_chains.py \
  --input_path=/data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --output_path=/data/lmk/ProteinMPNN/mpnn_outputs/assigned_pdbs.jsonl \
  --chain_list "A C"

python /data/lmk/ProteinMPNN/helper_scripts/make_fixed_positions_dict.py \
  --input_path=/data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --output_path=/data/lmk/ProteinMPNN/mpnn_outputs/fixed_pdbs.jsonl \
  --chain_list "A C" \
  --position_list "1 2 3 4 5 6 7 8 23 25, 10 11 12 13 14 15 16 17 18 19 20 40"

python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
  --jsonl_path /data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --chain_id_jsonl /data/lmk/ProteinMPNN/mpnn_outputs/assigned_pdbs.jsonl \
  --fixed_positions_jsonl /data/lmk/ProteinMPNN/mpnn_outputs/fixed_pdbs.jsonl \
  --out_folder /data/lmk/ProteinMPNN/mpnn_outputs \
  --num_seq_per_target 5 \
  --sampling_temp "0.1" \
  --seed 37 \
  --batch_size 1
```

> **04 仅设计特定残基 -- 反向约束**

**输入**：`mpnn_inputs/` 放复合物 PDB（如 3HTN.pdb）  
**功能**：与 03 相反，只重设计列出的少数位点，其余全部保持原样

在 `make_fixed_positions_dict.py` 末尾加 `--specify_non_fixed`，`--position_list` 列出的位点变成“需要被设计的位点”，其余全部固定
```bash
python /data/lmk/ProteinMPNN/helper_scripts/parse_multiple_chains.py \
  --input_path=/data/lmk/ProteinMPNN/mpnn_inputs \
  --output_path=/data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl

python /data/lmk/ProteinMPNN/helper_scripts/assign_fixed_chains.py \
  --input_path=/data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --output_path=/data/lmk/ProteinMPNN/mpnn_outputs/assigned_pdbs.jsonl \
  --chain_list "A C"

python /data/lmk/ProteinMPNN/helper_scripts/make_fixed_positions_dict.py \
  --input_path=/data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --output_path=/data/lmk/ProteinMPNN/mpnn_outputs/fixed_pdbs.jsonl \
  --chain_list "A C" \
  --position_list "1 2 3 4 5 6 7 8 23 25, 10 11 12 13 14 15 16 17 18 19 20 40" \
  --specify_non_fixed

python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
  --jsonl_path /data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --chain_id_jsonl /data/lmk/ProteinMPNN/mpnn_outputs/assigned_pdbs.jsonl \
  --fixed_positions_jsonl /data/lmk/ProteinMPNN/mpnn_outputs/fixed_pdbs.jsonl \
  --out_folder /data/lmk/ProteinMPNN/mpnn_outputs \
  --num_seq_per_target 5 \
  --sampling_temp "0.1" \
  --seed 37 \
  --batch_size 1
```

> **05 自定义残基绑定 -- tied_positions 协同设计**

**输入**：`mpnn_inputs/` 放复合物 PDB（如 3HTN.pdb）  
**功能**：把多条链的对应位点绑定，强制它们设计成相同氨基酸

`make_tied_positions_dict.py` 生成 `tied_pdbs.jsonl`，把不同链上对应的位点「绑定」，run 时用 `--tied_positions_jsonl`。下例把 A、B 两链的 1-8 位两两绑定
```bash
python /data/lmk/ProteinMPNN/helper_scripts/parse_multiple_chains.py \
  --input_path=/data/lmk/ProteinMPNN/mpnn_inputs \
  --output_path=/data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl

python /data/lmk/ProteinMPNN/helper_scripts/make_tied_positions_dict.py \
  --input_path=/data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --output_path=/data/lmk/ProteinMPNN/mpnn_outputs/tied_pdbs.jsonl \
  --chain_list "A B" \
  --position_list "1 2 3 4 5 6 7 8, 1 2 3 4 5 6 7 8"

python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
  --jsonl_path /data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --tied_positions_jsonl /data/lmk/ProteinMPNN/mpnn_outputs/tied_pdbs.jsonl \
  --out_folder /data/lmk/ProteinMPNN/mpnn_outputs \
  --num_seq_per_target 5 \
  --sampling_temp "0.1" \
  --seed 37 \
  --batch_size 1
```

> **06 同源寡聚体 -- C2/C3/C4 对称设计**

**输入**：`mpnn_inputs/` 放同源多聚体 PDB（如 4GYT.pdb、6EHB.pdb）  
**功能**：自动绑定所有链的对应位点，得到各链序列完全相同的同源寡聚体

同源寡聚体的各链等长、序列应一致，`make_tied_positions_dict.py` 加 `--homooligomer 1` 自动把所有链对应位点全部绑定，输出里每条链序列完全相同
```bash
python /data/lmk/ProteinMPNN/helper_scripts/parse_multiple_chains.py \
  --input_path=/data/lmk/ProteinMPNN/mpnn_inputs \
  --output_path=/data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl

python /data/lmk/ProteinMPNN/helper_scripts/make_tied_positions_dict.py \
  --input_path=/data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --output_path=/data/lmk/ProteinMPNN/mpnn_outputs/tied_pdbs.jsonl \
  --homooligomer 1

python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
  --jsonl_path /data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --tied_positions_jsonl /data/lmk/ProteinMPNN/mpnn_outputs/tied_pdbs.jsonl \
  --out_folder /data/lmk/ProteinMPNN/mpnn_outputs \
  --num_seq_per_target 5 \
  --sampling_temp "0.1" \
  --seed 37 \
  --batch_size 1
```

> **07 设置骨架噪声**

**输入**：此处直接复用 01 的 5L33.pdb、6MRR.pdb
**功能**：设计前给骨架坐标加高斯噪声，增加序列多样性

run 时加 `--backbone_noise 0.10`（单位 Å），注意这里的噪声是指的推理的骨架噪声，训练时默认是有噪声的（比如 v_48_020 有 0.2 的训练噪声）
```bash
python /data/lmk/ProteinMPNN/helper_scripts/parse_multiple_chains.py \
  --input_path=/data/lmk/ProteinMPNN/mpnn_inputs \
  --output_path=/data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl

python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
  --jsonl_path /data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --out_folder /data/lmk/ProteinMPNN/mpnn_outputs \
  --num_seq_per_target 5 \
  --sampling_temp "0.1" \
  --seed 37 \
  --batch_size 1 \
  --backbone_noise 0.10
```

> **08 氨基酸偏置 -- 鼓励或抑制特定氨基酸残基**

**输入**：此处直接复用 01 的 5L33.pdb、6MRR.pdb 
**功能**：在全局层面调整氨基酸偏好性

`make_bias_AA.py` 生成全局偏置表（log 空间，正值偏好 / 负值抑制），run 时用 `--bias_AA_jsonl`，下例整体偏向极性氨基酸
```bash
python /data/lmk/ProteinMPNN/helper_scripts/parse_multiple_chains.py \
  --input_path=/data/lmk/ProteinMPNN/mpnn_inputs \
  --output_path=/data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl

python /data/lmk/ProteinMPNN/helper_scripts/make_bias_AA.py \
  --output_path=/data/lmk/ProteinMPNN/mpnn_outputs/bias_AA.jsonl \
  --AA_list="D E H K N Q R S T W Y" \
  --bias_list="1.39 1.39 1.39 1.39 1.39 1.39 1.39 1.39 1.39 1.39 1.39"

python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
  --jsonl_path /data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --bias_AA_jsonl /data/lmk/ProteinMPNN/mpnn_outputs/bias_AA.jsonl \
  --out_folder /data/lmk/ProteinMPNN/mpnn_outputs \
  --num_seq_per_target 5 \
  --sampling_temp "0.1" \
  --seed 37 \
  --batch_size 1
```

---

### 二、进阶参数（均在 01 基线上只加一项）

> **09 全局排除氨基酸 -- |通用|禁用某些 AA|--omit_AAs|**

**输入**：同 01（单体骨架）；**功能**：完全禁止指定氨基酸出现在任何位点（这里禁用半胱氨酸 C）

`--omit_AAs` 禁止某些氨基酸出现在所有位点（默认 `"X"`）。下例禁用半胱氨酸 C，避免设计出游离巯基
```bash
python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
  --jsonl_path /data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --out_folder /data/lmk/ProteinMPNN/mpnn_outputs \
  --num_seq_per_target 5 \
  --sampling_temp "0.1" \
  --seed 37 \
  --batch_size 1 \
  --omit_AAs "C"
```

> [!TIP]
> **做抗体 CDR 设计时建议把 `--omit_AAs "C"` 当默认参数加上，** CDR 里多出的游离巯基会与之竞争错配，导致异构体、聚集、批次不一致。

> **10 仅打分（不设计）-- |通用|给原始序列打分|--score_only 1|**

**输入**：同 01（单体骨架）；**功能**：不做设计，只给输入 PDB 自带的天然序列打分

`--score_only 1` 不做设计，只对输入 PDB 自带序列算 score / global_score，结果存 `.npz`（在 `<out>/score_only/`）。用于评估天然序列或对比设计前后
```bash
python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
  --jsonl_path /data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --out_folder /data/lmk/ProteinMPNN/mpnn_outputs \
  --seed 37 \
  --batch_size 1 \
  --score_only 1
```

> **11 对外部序列打分 -- |通用|给定序列打分|--path_to_fasta|**

**输入**：骨架同 01，外加 `mpnn_inputs/` 里的待评序列 `seqs_to_score.fa`；**功能**：在该骨架上给一条外部序列打分（质检验证回来的序列）

配合 `--score_only 1`，用 `--path_to_fasta` 在同一 backbone 上给外部序列（如 AF2/AF3 验证回来的序列）打分。FASTA 里多链用 `/` 分隔，链顺序按字母序
```bash
python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
  --jsonl_path /data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --out_folder /data/lmk/ProteinMPNN/mpnn_outputs \
  --seed 37 \
  --batch_size 1 \
  --score_only 1 \
  --path_to_fasta /data/lmk/ProteinMPNN/mpnn_inputs/seqs_to_score.fa
```

> **12 soluble 模型 -- |通用|抑制表面疏水|--use_soluble_model|**

**输入**：同 01（单体骨架）；**功能**：换用 soluble 权重，设计可溶蛋白时倾向减少表面疏水残基

加 `--use_soluble_model` 切到只在可溶蛋白上训练的权重，设计可溶蛋白时减少暴露疏水残基（含义见 [FASTA_Format](./ProteinMPNN_FASTA_Format.md) 的 `model_name` 解码表）
```bash
python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
  --jsonl_path /data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --out_folder /data/lmk/ProteinMPNN/mpnn_outputs \
  --num_seq_per_target 5 \
  --sampling_temp "0.1" \
  --seed 37 \
  --batch_size 1 \
  --use_soluble_model
```

> **13 Cα-only 模型 -- |通用|仅用 Cα 坐标|--ca_only|**

**输入**：同 01（单体骨架，仅读 Cα 坐标）；**功能**：用 Cα-only 权重设计，适合骨架原子缺失或只有 Cα trace 的输入

加 `--ca_only` 切到 Cα-only 权重，只读 Cα 坐标，适合骨架原子缺失或只有 Cα trace 的输入
```bash
python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
  --jsonl_path /data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --out_folder /data/lmk/ProteinMPNN/mpnn_outputs \
  --num_seq_per_target 5 \
  --sampling_temp "0.1" \
  --seed 37 \
  --batch_size 1 \
  --ca_only
```

> **14 切换 vanilla 噪声档 -- |通用|挑选权重档位|--model_name|**

**输入**：同 01（单体骨架）；**功能**：切换 vanilla 权重的训练噪声档位，权衡稳健性与对骨架精度的依赖

`--model_name` 选 vanilla 权重档位（默认 `v_48_020`）。下例换 `v_48_002`（训练噪声 0.02 Å，依赖高精度骨架，适合晶体结构）
```bash
python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
  --jsonl_path /data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --out_folder /data/lmk/ProteinMPNN/mpnn_outputs \
  --num_seq_per_target 5 \
  --sampling_temp "0.1" \
  --seed 37 \
  --batch_size 1 \
  --model_name v_48_002
```

> **15 保存每位点概率 / 打分 -- |通用|导出 npz|--save_probs / --save_score|**

**输入**：同 01（单体骨架）；**功能**：设计的同时导出每位点的氨基酸概率（供保守性 / 熵分析）

`--save_probs 1` 把每位点 20 种氨基酸的概率存成 `.npz`（在 `<out>/probs/`），可做位点保守性 / 熵分析；`--save_score 1` 另存逐序列 score
```bash
python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
  --jsonl_path /data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --out_folder /data/lmk/ProteinMPNN/mpnn_outputs \
  --num_seq_per_target 5 \
  --sampling_temp "0.1" \
  --seed 37 \
  --batch_size 1 \
  --save_probs 1
```

> **16 条件概率分析 -- |通用|p(s_i | 骨架,其余序列)|--conditional_probs_only 1|**

**输入**：同 01（单体骨架）；**功能**：不出序列，只导出「给定骨架 + 其余序列」下的逐位点条件概率

`--conditional_probs_only 1` 不采样，只输出每位点在「给定骨架 + 其余序列」下的条件概率（存 `.npz`），用于分析某个位点对上下文的依赖
```bash
python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
  --jsonl_path /data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --out_folder /data/lmk/ProteinMPNN/mpnn_outputs \
  --seed 37 \
  --batch_size 1 \
  --conditional_probs_only 1
```

> **17 无条件概率分析 -- |通用|p(s_i | 骨架)|--unconditional_probs_only 1|**

**输入**：同 01（单体骨架）；**功能**：不出序列，只导出「仅给定骨架」下的逐位点无条件概率

`--unconditional_probs_only 1` 一次前向就给出每位点「仅给定骨架」的无条件概率，比 16 快，适合快速看骨架本身偏好哪些氨基酸
```bash
python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
  --jsonl_path /data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --out_folder /data/lmk/ProteinMPNN/mpnn_outputs \
  --seed 37 \
  --batch_size 1 \
  --unconditional_probs_only 1
```

> **18 单 PDB 直接输入 -- |单文件|跳过 parse|--pdb_path / --pdb_path_chains|**

**输入**：`mpnn_inputs/` 里单个 PDB 文件（如 5L33.pdb，无需 parse）；**功能**：跳过解析步骤，直接对一个 PDB 设计

只设计一个 PDB 时可跳过 `parse_multiple_chains.py`，用 `--pdb_path` 直接传文件、`--pdb_path_chains` 指定要设计的链（其余链固定）
```bash
python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
  --pdb_path /data/lmk/ProteinMPNN/mpnn_inputs/5L33.pdb \
  --pdb_path_chains "A" \
  --out_folder /data/lmk/ProteinMPNN/mpnn_outputs \
  --num_seq_per_target 5 \
  --sampling_temp "0.1" \
  --seed 37 \
  --batch_size 1
```

> **19 采样温度多值扫描 -- |通用|一次跑多个温度|--sampling_temp|**

**输入**：同 01（单体骨架）；**功能**：一次运行在多个采样温度下各采一批序列，方便比较多样性

`--sampling_temp` 接受多个空格分隔的温度，一次性在每个温度各采一批序列（温度越高越多样、置信度越低），FASTA header 里的 `T=` 会标明各自温度
```bash
python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
  --jsonl_path /data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --out_folder /data/lmk/ProteinMPNN/mpnn_outputs \
  --num_seq_per_target 5 \
  --sampling_temp "0.1 0.2 0.3" \
  --seed 37 \
  --batch_size 1
```

##### [ProteinMPNN官方文档](https://github.com/dauparas/ProteinMPNN)
