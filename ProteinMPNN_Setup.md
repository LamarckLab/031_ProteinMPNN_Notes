## Lamarck &nbsp; &nbsp; &nbsp; 2025-11-01
#### 该文档用于部署 ProteinMPNN
---

## 01  克隆官方的源码仓库
> **236 机子路径 /data/lmk/ProteinMPNN**
```bash
cd /data/lmk
git clone https://github.com/dauparas/ProteinMPNN.git
```

## 02  创建 conda 环境
> ProteinMPNN 依赖很轻，官方未提供 yml 文件，自己拉一个干净环境即可
```bash
conda create -n lmk_ProteinMPNN python=3.11 -y --solver=classic
conda activate lmk_ProteinMPNN
```

## 03  安装 PyTorch 与 NumPy
> 官方要求 PyTorch >= 1.8，amax 上 CUDA 12.8，用 PyTorch 官方 cu126 wheel 兼容
```bash
python -m pip install --index-url https://download.pytorch.org/whl/cu126 torch
python -m pip install numpy
```

## 04  验证安装
```bash
cd /data/lmk/ProteinMPNN
python protein_mpnn_run.py --help
```

## 05  ProteinMPNN 输入输出目录
**236 机子路径**
> 输入目录：/data/lmk/ProteinMPNN/mpnn_inputs  
> 输出目录：/data/lmk/ProteinMPNN/mpnn_outputs

##### [ProteinMPNN官方文档](https://github.com/dauparas/ProteinMPNN)
