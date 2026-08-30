<p align="left">
  <a href="./README_EN.md">Homepage</a>
</p>
<p align="right">
  <strong>English</strong> | <a href="../ProteinMPNN_FASTA_Format.md">中文</a>
</p>

## Lamarck &nbsp; &nbsp; &nbsp; 2025-11-01
#### This document explains the fields of the .fa file ProteinMPNN writes
---

> **00 Basic FASTA structure**

MPNN writes one `.fa` file per input PDB (under `mpnn_outputs/seqs/`), holding a number of (header, sequence) pairs:
- **The 1st record**: the original sequence of the input PDB (the reference line)
- **From the 2nd record on**: the new sequences the model designed (as many as `--num_seq_per_target`)

---

> **01 The 1st record: the input reference line**

```
>monomer, score=1.3597, global_score=1.3597, fixed_chains=[], designed_chains=['A'], model_name=v_48_020, git_hash=8907e667..., seed=37
GGGGGGGGGGGGGGG...
```

| Field | Meaning |
| :--- | :--- |
| `monomer` | the input PDB filename, which is also the ID of this FASTA record |
| `score` | the model's negative log-likelihood of the **input sequence**, **averaged over the designable positions only**; lower means the model likes it more |
| `global_score` | the same negative log-likelihood, but averaged over **every residue with coordinates in the complex** (fixed chains included) |
| `fixed_chains` | the IDs of the chains held fixed (not redesigned) |
| `designed_chains` | the IDs of the chains that were redesigned |
| `model_name` | the model weights used; see the decoding table below |
| `git_hash` | the commit hash of the ProteinMPNN repo, for reproducibility |
| `seed` | the random seed, matching `--seed` on the CLI |

Decoding `model_name`:

| Name | Meaning |
| :--- | :--- |
| `v_48_020` | vanilla model, k-NN=48, 0.20 Å backbone noise during training (**the default**, the most robust, and forgiving of rough backbones such as RFdiffusion output) |
| `v_48_010` | vanilla model, k-NN=48, 0.10 Å training noise |
| `v_48_002` | vanilla model, k-NN=48, 0.02 Å training noise (depends more on backbone accuracy, suited to crystal structures) |
| `s_48_020` | soluble model, penalises surface hydrophobic residues, suited to soluble protein design |
| `c_48_020` | Cα-only model, uses Cα coordinates alone, suited to inputs missing backbone atoms or holding only Cα |

> In the example above, `score=1.3597` is the score of the input sequence — this is a poly-G backbone (from RFdiffusion output), and the model naturally does not think much of an all-G sequence, so the score is on the high side.

---

> **02 The 2nd and later records: the sequences the model designed**

```
>T=0.1, sample=1, score=0.9813, global_score=0.9813, seq_recovery=0.0133
MEKEKIKEKLKEIREKIE...
```

| Field | Meaning |
| :--- | :--- |
| `T=0.1` | the sampling temperature, matching `--sampling_temp` on the CLI; higher values give more diverse sequences but lower confidence |
| `sample=N` | this design is the Nth sample (N runs from 1 to `--num_seq_per_target`) |
| `score` | the model's negative log-likelihood of **this designed sequence**, over the same range as in 01 (designable positions only); usually lower than the reference line's score (the model "likes" its own design better) |
| `global_score` | as in 01, averaged over every residue with coordinates in the complex (fixed chains included) |
| `seq_recovery` | how often the design matches the input sequence per position (0-1); the denominator again counts **designable positions only**, not the whole complex |

---

> **03 Multi-chain complex output**

When several chains are designed, **`/` separates the chains inside each sequence**. Below is real 3HTN output (truncated), with chains A and B designed and chain C fixed:

```
>3HTN, score=1.1514, global_score=1.2018, fixed_chains=['C'], designed_chains=['A', 'B'], model_name=v_48_020, git_hash=8907e667..., seed=37
NMYSYKKIGNKYIVSINNHTEI...RTYNPDLGLNIYDFER/NMYSYKKIGNKYIVSINNHTEI...LRFFNPKXXXXDDKTFREQ...RTYNPDLGLNIYDFER

>T=0.1, sample=1, score=0.7382, global_score=0.9213, seq_recovery=0.5567
NMYKYKEIGNKYIVSINNNTDL...YKYDEELGLYLLDFNK/HMYSYKKIGNKYIVSINNGQDL...LSFFDPNXXXXTTKTFNDY...VKYNEETGLYLLDFDL
```

---

> **04 Homo-oligomer (tied positions) output**

In the homo-oligomer design of 06 (`--homooligomer 1`) the chains share one set of residue decisions, so every chain in the output carries **exactly the same** sequence:

```
>homomer, score=1.40, global_score=1.40, fixed_chains=[], designed_chains=['A','B','C'], ...
GGGGG.../GGGGG.../GGGGG...

>T=0.2, sample=1, score=1.00, global_score=1.00, seq_recovery=0.0
MKLLV.../MKLLV.../MKLLV...
```


---

##### [ProteinMPNN official documentation](https://github.com/dauparas/ProteinMPNN)
