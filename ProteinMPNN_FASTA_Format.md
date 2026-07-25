## Lamarck &nbsp; &nbsp; &nbsp; 2025-11-01
#### 该文档用于展示 ProteinMPNN 输出 .fa 文件的字段含义
---

> **00 FASTA 基本结构**

MPNN 为每个输入 PDB 输出一个 `.fa` 文件（位于 `mpnn_outputs/seqs/`），文件内由若干条 (header, sequence) 对组成：
- **第 1 条**：输入 PDB 的原始序列（参考行）
- **第 2 条起**：模型设计的新序列（数量 = `--num_seq_per_target`）

---

> **01 第 1 条记录：输入参考行**

```
>monomer, score=1.3597, global_score=1.3597, fixed_chains=[], designed_chains=['A'], model_name=v_48_020, git_hash=8907e667..., seed=37
GGGGGGGGGGGGGGG...
```

| 字段              | 含义                                                                         |
| :---------------- | :--------------------------------------------------------------------------- |
| `monomer`         | 输入 PDB 文件名，同时也是这条 FASTA 的 ID                                    |
| `score`           | 模型对**输入序列**的负对数似然，**只在可设计位点上平均**；越低代表模型越认可 |
| `global_score`    | 同样是负对数似然，但在**复合物全部有坐标残基**上平均（含固定链）             |
| `fixed_chains`    | 设计时被固定（不参与重新设计）的链 ID 列表                                   |
| `designed_chains` | 设计时参与重新设计的链 ID 列表                                               |
| `model_name`      | 使用的模型权重，含义见下方解码表                                             |
| `git_hash`        | ProteinMPNN 仓库的 commit 哈希，用于复现                                     |
| `seed`            | 随机种子，对应 CLI 的 `--seed`                                               |

`model_name` 解码：

| 名称       | 含义                                                                                            |
| :--------- | :---------------------------------------------------------------------------------------------- |
| `v_48_020` | vanilla 模型，k-NN=48，训练时骨架噪声 0.20 Å（**默认**，最稳健，对 RFdiffusion 等粗糙骨架友好） |
| `v_48_010` | vanilla 模型，k-NN=48，训练噪声 0.10 Å                                                          |
| `v_48_002` | vanilla 模型，k-NN=48，训练噪声 0.02 Å（对骨架精度依赖较高，适合晶体结构）                      |
| `s_48_020` | soluble 模型，惩罚表面疏水残基，适合做可溶蛋白                                                  |
| `c_48_020` | Cα-only 模型，只用 Cα 坐标，适合骨架原子缺失或只有 Cα 的输入                                    |

> 上例中 `score=1.3597` 是输入序列的 score —— 因为这是一个 poly-G 骨架（来自 RFdiffusion 输出），模型自然不会很认可全是 G 的序列，所以这个 score 偏高。

---

> **02 第 2+ 条记录：模型设计的序列**

```
>T=0.1, sample=1, score=0.9813, global_score=0.9813, seq_recovery=0.0133
MEKEKIKEKLKEIREKIE...
```

| 字段           | 含义                                                                                                            |
| :------------- | :-------------------------------------------------------------------------------------------------------------- |
| `T=0.1`        | 采样温度，对应 CLI 的 `--sampling_temp`；值越大序列越多样但置信度下降                                           |
| `sample=N`     | 该设计是第 N 个样本（N 从 1 到 `--num_seq_per_target`）                                                         |
| `score`        | 模型对**这条设计序列**的负对数似然，范围同 01（仅可设计位点）；一般低于参考行的 score（模型更"认可"自己的设计） |
| `global_score` | 同 01，在复合物全部有坐标残基上平均（含固定链）                                                                 |
| `seq_recovery` | 与输入序列的位点重合率（0-1），分母同样只数**可设计位点**，不是整个复合物                                       |

---

> **03 多链复合物输出**

设计多链时，**每条序列内部用 `/` 分隔不同链**。下面是 3HTN 的真实输出（截断显示），A、B 两链参与设计，C 链固定：

```
>3HTN, score=1.1514, global_score=1.2018, fixed_chains=['C'], designed_chains=['A', 'B'], model_name=v_48_020, git_hash=8907e667..., seed=37
NMYSYKKIGNKYIVSINNHTEI...RTYNPDLGLNIYDFER/NMYSYKKIGNKYIVSINNHTEI...LRFFNPKXXXXDDKTFREQ...RTYNPDLGLNIYDFER

>T=0.1, sample=1, score=0.7382, global_score=0.9213, seq_recovery=0.5567
NMYKYKEIGNKYIVSINNNTDL...YKYDEELGLYLLDFNK/HMYSYKKIGNKYIVSINNGQDL...LSFFDPNXXXXTTKTFNDY...VKYNEETGLYLLDFDL
```

---

> **04 同源多聚体（tied positions）输出**

06 同源寡聚体设计 (`--homooligomer 1`) 时，多条链共享同一套残基决策，输出里所有链序列**完全相同**：

```
>homomer, score=1.40, global_score=1.40, fixed_chains=[], designed_chains=['A','B','C'], ...
GGGGG.../GGGGG.../GGGGG...

>T=0.2, sample=1, score=1.00, global_score=1.00, seq_recovery=0.0
MKLLV.../MKLLV.../MKLLV...
```


---

##### [ProteinMPNN官方文档](https://github.com/dauparas/ProteinMPNN)
