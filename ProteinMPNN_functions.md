## Lamarck &nbsp; &nbsp; &nbsp; 2025-10-30
#### 该文档用于记录 server 上跑 ProteinMPNN 的各种命令
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
9 个示例的样例 PDB 按 sample 分组放在仓库的 [inputs/](./inputs/) 目录下，跑某 sample 前先把对应文件传到 server 的 mpnn_inputs/ 即可

*GPU选择*
```bash
export CUDA_DEVICE_ORDER=PCI_BUS_ID
export CUDA_VISIBLE_DEVICES=3
```

---

### 01 Basic Design -- 基础设计

> **01.1 简单单体设计 -- 标准 jsonl 工作流**

最基础的入门脚本，先用 `parse_multiple_chains.py` 把 PDB 解析成 `parsed_pdbs.jsonl` (记录链 ID、残基编号、坐标/掩码等元数据)，再喂给 `protein_mpnn_run.py`
```bash
folder_with_pdbs="/data/lmk/ProteinMPNN/mpnn_inputs/"
output_dir="/data/lmk/ProteinMPNN/mpnn_outputs"

path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"

python /data/lmk/ProteinMPNN/helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains

python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
        --jsonl_path $path_for_parsed_chains \
        --out_folder $output_dir \
        --num_seq_per_target 5 \
        --sampling_temp "0.1" \
        --seed 37 \
        --batch_size 1
```

> **01.2 直接 PDB 路径输入 -- 跳过 jsonl 预解析**

不需要先生成 `.jsonl`，可以直接指定 PDB 文件路径运行，方便快速测试单个 PDB
```bash
path_to_PDB="/data/lmk/ProteinMPNN/mpnn_inputs/3HTN.pdb"
output_dir="/data/lmk/ProteinMPNN/mpnn_outputs"

chains_to_design="A B"

python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
        --pdb_path $path_to_PDB \
        --pdb_path_chains "$chains_to_design" \
        --out_folder $output_dir \
        --num_seq_per_target 2 \
        --sampling_temp "0.1" \
        --seed 37 \
        --batch_size 1
```

---

### 02 Complex & Residue Constraints -- 复合物与残基约束

> **02.1 多链复合物 -- 指定设计链**

抗原-抗体、受体-配体等多链复合物，用 `assign_fixed_chains.py` 生成 `assigned_pdbs.jsonl`，在 parsed 基础上写入"哪些链设计、哪些链固定"的分配信息
```bash
folder_with_pdbs="/data/lmk/ProteinMPNN/mpnn_inputs/"
output_dir="/data/lmk/ProteinMPNN/mpnn_outputs"

path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"
path_for_assigned_chains=$output_dir"/assigned_pdbs.jsonl"

chains_to_design="A B"  # 表示把链A和B作为待设计链，其余链自动视作固定链

python /data/lmk/ProteinMPNN/helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains
python /data/lmk/ProteinMPNN/helper_scripts/assign_fixed_chains.py --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design"

python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
        --jsonl_path $path_for_parsed_chains \
        --chain_id_jsonl $path_for_assigned_chains \
        --out_folder $output_dir \
        --num_seq_per_target 2 \
        --sampling_temp "0.1" \
        --seed 37 \
        --batch_size 1
```

> **02.2 固定特定残基 -- 部分残基不参与设计**

固定活性位点、金属配位残基、抗原表位等。用 `make_fixed_positions_dict.py` 生成 `fixed_pdbs.jsonl`，指定哪些残基位置需要固定 (即不被设计)
```bash
folder_with_pdbs="/data/lmk/ProteinMPNN/mpnn_inputs/"
output_dir="/data/lmk/ProteinMPNN/mpnn_outputs"

path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"
path_for_assigned_chains=$output_dir"/assigned_pdbs.jsonl"
path_for_fixed_positions=$output_dir"/fixed_pdbs.jsonl"

chains_to_design="A C"

fixed_positions="1 2 3 4 5 6 7 8 23 25, 10 11 12 13 14 15 16 17 18 19 20 40"  # 这里的序号是严格的残基排序，而不是pdb中的残基index

python /data/lmk/ProteinMPNN/helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains
python /data/lmk/ProteinMPNN/helper_scripts/assign_fixed_chains.py --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design"
python /data/lmk/ProteinMPNN/helper_scripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list "$chains_to_design" --position_list "$fixed_positions"

python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
        --jsonl_path $path_for_parsed_chains \
        --chain_id_jsonl $path_for_assigned_chains \
        --fixed_positions_jsonl $path_for_fixed_positions \
        --out_folder $output_dir \
        --num_seq_per_target 2 \
        --sampling_temp "0.1" \
        --seed 37 \
        --batch_size 1
```

> **02.3 仅设计特定残基 -- 反向约束**

与上一例相反，加 `--specify_non_fixed` 让指定的残基可被设计、其余全部保持固定，常用于探索性突变扫描
```bash
folder_with_pdbs="/data/lmk/ProteinMPNN/mpnn_inputs/"
output_dir="/data/lmk/ProteinMPNN/mpnn_outputs"

path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"
path_for_assigned_chains=$output_dir"/assigned_pdbs.jsonl"
path_for_fixed_positions=$output_dir"/fixed_pdbs.jsonl"

chains_to_design="A C"

design_only_positions="1 2 3 4 5 6 7 8 9 10, 3 4 5 6 7 8"  # 这里的序号是严格的残基排序，而不是pdb中的残基index

python /data/lmk/ProteinMPNN/helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains
python /data/lmk/ProteinMPNN/helper_scripts/assign_fixed_chains.py --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design"
python /data/lmk/ProteinMPNN/helper_scripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list "$chains_to_design" --position_list "$design_only_positions" --specify_non_fixed

python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
        --jsonl_path $path_for_parsed_chains \
        --chain_id_jsonl $path_for_assigned_chains \
        --fixed_positions_jsonl $path_for_fixed_positions \
        --out_folder $output_dir \
        --num_seq_per_target 2 \
        --sampling_temp "0.1" \
        --seed 37 \
        --batch_size 1
```

---

### 03 Symmetry -- 对称性设计

> **03.1 自定义残基绑定 -- tied_positions 协同设计**

把多条链中对应残基绑定，使它们使用相同的氨基酸类型。`tied_pdbs.jsonl` 保存绑定位点 (即必须一起设计、同步变动的残基位置)，常用于对称多聚体、重复结构、功能相关位点的协同设计
```bash
folder_with_pdbs="/data/lmk/ProteinMPNN/mpnn_inputs/"
output_dir="/data/lmk/ProteinMPNN/mpnn_outputs"

path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"
path_for_assigned_chains=$output_dir"/assigned_pdbs.jsonl"
path_for_fixed_positions=$output_dir"/fixed_pdbs.jsonl"
path_for_tied_positions=$output_dir"/tied_pdbs.jsonl"

chains_to_design="A C"

fixed_positions="9 10 11 12 13 14 15 16 17 18 19 20 21 22 23, 10 11 18 19 20 22"
tied_positions="1 2 3 4 5 6 7 8, 1 2 3 4 5 6 7 8"  # 指定哪些位置是绑定的，也就是必须一起设计和变动的残基，两个链中对应的残基会同步变动

python /data/lmk/ProteinMPNN/helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains
python /data/lmk/ProteinMPNN/helper_scripts/assign_fixed_chains.py --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design"
python /data/lmk/ProteinMPNN/helper_scripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list "$chains_to_design" --position_list "$fixed_positions"
python /data/lmk/ProteinMPNN/helper_scripts/make_tied_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_tied_positions --chain_list "$chains_to_design" --position_list "$tied_positions"

python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
        --jsonl_path $path_for_parsed_chains \
        --chain_id_jsonl $path_for_assigned_chains \
        --fixed_positions_jsonl $path_for_fixed_positions \
        --tied_positions_jsonl $path_for_tied_positions \
        --out_folder $output_dir \
        --num_seq_per_target 2 \
        --sampling_temp "0.1" \
        --seed 37 \
        --batch_size 1
```

> **03.2 同源寡聚体 -- C2/C3/C4 对称设计**

加 `--homooligomer 1` 直接对同源对称重复单元做对称序列生成，多条链共享序列信息
```bash
folder_with_pdbs="/data/lmk/ProteinMPNN/mpnn_inputs/"
output_dir="/data/lmk/ProteinMPNN/mpnn_outputs"

path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"
path_for_tied_positions=$output_dir"/tied_pdbs.jsonl"

python /data/lmk/ProteinMPNN/helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains
python /data/lmk/ProteinMPNN/helper_scripts/make_tied_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_tied_positions --homooligomer 1

python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
        --jsonl_path $path_for_parsed_chains \
        --tied_positions_jsonl $path_for_tied_positions \
        --out_folder $output_dir \
        --num_seq_per_target 2 \
        --sampling_temp "0.2" \
        --seed 37 \
        --batch_size 1
```

---

### 04 Probability & Bias -- 概率输出与氨基酸偏置

> **04.1 输出无条件概率 -- 模型内部的 PSSM**

加 `--unconditional_probs_only 1`，不生成确定序列，而输出每个位置上 20 个氨基酸的概率分布，可用于统计分析或后续能量计算
```bash
folder_with_pdbs="/data/lmk/ProteinMPNN/mpnn_inputs/"
output_dir="/data/lmk/ProteinMPNN/mpnn_outputs"

path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"

python /data/lmk/ProteinMPNN/helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains

python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
        --jsonl_path $path_for_parsed_chains \
        --out_folder $output_dir \
        --num_seq_per_target 1 \
        --sampling_temp "0.1" \
        --unconditional_probs_only 1 \
        --seed 37 \
        --batch_size 1
```

> **04.2 氨基酸偏置 -- 鼓励或抑制特定残基**

用 `make_bias_AA.py` 生成 `bias_pdbs.jsonl`，加到模型 logits 上做加法偏置 (正数=鼓励、负数=抑制、0=不变)。`1.39 ≈ ln(4)`，相对几率提高约 4 倍。例如想全局偏好芳香族就列出 `F W Y`；想强烈抑制 C 防止二硫键就给 `-2.0` 之类的负偏置；不改动就别把该氨基酸放入列表，或给 0
```bash
folder_with_pdbs="/data/lmk/ProteinMPNN/mpnn_inputs/"
output_dir="/data/lmk/ProteinMPNN/mpnn_outputs"

path_for_bias=$output_dir"/bias_pdbs.jsonl"
path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"

AA_list="D E H K N Q R S T W Y"  # 用单字母列出要偏置的氨基酸，这里是11种: D, E, H, K, N, Q, R, S, T, W, Y
bias_list="1.39 1.39 1.39 1.39 1.39 1.39 1.39 1.39 1.39 1.39 1.39"  # 与AA_list一一对应的偏置数值

python /data/lmk/ProteinMPNN/helper_scripts/make_bias_AA.py --output_path=$path_for_bias --AA_list="$AA_list" --bias_list="$bias_list"
python /data/lmk/ProteinMPNN/helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains

python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
        --jsonl_path $path_for_parsed_chains \
        --out_folder $output_dir \
        --bias_AA_jsonl $path_for_bias \
        --num_seq_per_target 2 \
        --sampling_temp "0.1" \
        --seed 37 \
        --batch_size 1
```

##### [ProteinMPNN官方文档](https://github.com/dauparas/ProteinMPNN)
