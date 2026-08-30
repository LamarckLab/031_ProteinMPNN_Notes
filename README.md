<p align="right">
  <a href="./docs_EN/README_EN.md">English</a> | <strong>中文</strong>
</p>

<h1 align="center">🧬 ProteinMPNN 部署与功能测试笔记</h1>

<p align="center"><em>—— 2025.11.01</em></p>

<p align="center">
  <img src="https://img.shields.io/badge/Tool-ProteinMPNN-blue?style=flat-square" />
  <img src="https://img.shields.io/badge/Field-Protein%20Design-orange?style=flat-square" />
  <img src="https://img.shields.io/badge/Platform-Linux-555?style=flat-square" />
  <img src="https://img.shields.io/badge/Status-Complete-brightgreen?style=flat-square" />
</p>

---

## 内容索引

| 文档                                                         | 说明                                                            |
| :----------------------------------------------------------- | :-------------------------------------------------------------- |
| [ProteinMPNN_Setup.md](./ProteinMPNN_Setup.md)               | 克隆源码、conda 环境创建、PyTorch 与 NumPy 安装                 |
| [ProteinMPNN_Functions.md](./ProteinMPNN_Functions.md)       | 各功能命令：链/位点固定、tied、同源寡聚、偏置、打分、概率导出等 |
| [ProteinMPNN_FASTA_Format.md](./ProteinMPNN_FASTA_Format.md) | 输出 .fa 文件的字段含义（score / seq_recovery / model_name 等） |
| [input_data/](./input_data/)                                 | 各功能测试用的示例 PDB，按 sample 1~8 分组                      |

---

##### [ProteinMPNN 官方文档](https://github.com/dauparas/ProteinMPNN)
