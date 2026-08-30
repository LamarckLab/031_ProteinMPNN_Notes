<p align="left">
  <a href="./README_EN.md">Homepage</a>
</p>
<p align="right">
  <strong>English</strong> | <a href="../ProteinMPNN_Setup.md">中文</a>
</p>

## Lamarck &nbsp; &nbsp; &nbsp; 2025-11-01
#### Deploying ProteinMPNN
---

## 01  Cloning the official source repository
> **Path on the 236 machine: /data/lmk/ProteinMPNN**
```bash
cd /data/lmk
git clone https://github.com/dauparas/ProteinMPNN.git
```

## 02  Creating the conda environment
> ProteinMPNN has light dependencies and the authors provide no yml file, so create a clean environment directly
```bash
conda create -n lmk_ProteinMPNN python=3.11 -y --solver=classic
conda activate lmk_ProteinMPNN
```

## 03  Installing PyTorch and NumPy
> PyTorch >= 1.8 is required. amax runs CUDA 12.8, so the official cu126 wheel is the compatible build
```bash
python -m pip install --index-url https://download.pytorch.org/whl/cu126 torch
python -m pip install numpy
```

## 04  Verifying the installation
```bash
cd /data/lmk/ProteinMPNN
python protein_mpnn_run.py --help
```

## 05  ProteinMPNN input and output directories
**Paths on the 236 machine**
> Input dir:  /data/lmk/ProteinMPNN/mpnn_inputs  
> Output dir: /data/lmk/ProteinMPNN/mpnn_outputs

##### [ProteinMPNN official documentation](https://github.com/dauparas/ProteinMPNN)
