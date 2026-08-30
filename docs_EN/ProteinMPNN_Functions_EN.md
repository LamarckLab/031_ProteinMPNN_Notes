<p align="left">
  <a href="./README_EN.md">Homepage</a>
</p>
<p align="right">
  <strong>English</strong> | <a href="../ProteinMPNN_Functions.md">中文</a>
</p>

## Lamarck &nbsp; &nbsp; &nbsp; 2025-10-30
#### This document records the commands for running ProteinMPNN on the server
---

*Environment & paths*
```bash
Env on the 236 server:  lmk_ProteinMPNN
Path on the 236 server: /data/lmk/ProteinMPNN/protein_mpnn_run.py
```

*Input & output paths*
```bash
Input dir:  /data/lmk/ProteinMPNN/mpnn_inputs    # PDB files used for each run
Output dir: /data/lmk/ProteinMPNN/mpnn_outputs   # the parsed jsonl and the generated FASTA
```

*GPU selection*
```bash
export CUDA_DEVICE_ORDER=PCI_BUS_ID
export CUDA_VISIBLE_DEVICES=3
```

> [!NOTE]
> Except for "single PDB straight in" (18), every feature first runs `parse_multiple_chains.py` to turn a folder of PDBs into `parsed_pdbs.jsonl`, which is then fed to `protein_mpnn_run.py --jsonl_path`. The extra constraints — fixed chains, fixed positions, tied positions, bias — all work the same way: the matching helper script writes a `.jsonl`, and the corresponding flag passes it in. Every command uses the same fixed pair of directories, `mpnn_inputs` / `mpnn_outputs`; to switch features, just swap the PDBs in `mpnn_inputs` and rerun (outputs with the same name are overwritten). Each example below changes only the one thing it demonstrates, on the shared baseline `--num_seq_per_target 5 --sampling_temp "0.1" --seed 37 --batch_size 1`.

---

### Part 1 — Core design workflow

> **01 Simple monomer design**

**Input**: put monomer PDBs in `mpnn_inputs/` (e.g. 5L33.pdb, 6MRR.pdb)  
**What it does**: with no constraints at all, design a sequence for the whole chain on the fixed backbone, giving 5 candidates per backbone
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

> **02 Multi-chain complex -- choosing which chains to design**

**Input**: put a complex PDB in `mpnn_inputs/` (e.g. 3HTN.pdb)  
**What it does**: redesign only chains A and B of the complex, leaving the rest at their original sequence as fixed context

Use `assign_fixed_chains.py` to write `assigned_pdbs.jsonl` marking which chains take part in the design, then pass it at run time with `--chain_id_jsonl`. The example below designs only chains A and B, leaving the others at their original sequence
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

> **03 Fixing specific residues -- keeping some residues out of the design**

**Input**: put a complex PDB in `mpnn_inputs/` (e.g. 3HTN.pdb)  
**What it does**: lock the named residues on the designed chains (a key active site, say) and redesign every other position

`make_fixed_positions_dict.py` writes `fixed_pdbs.jsonl`; inside `--position_list`, chains are separated by commas and positions within a chain by spaces. The example below locks 1-8/23/25 on chain A and 10-20/40 on chain C

> [!IMPORTANT]
> `--position_list` uses the **1-based ordinal position within the chain** after `parse_multiple_chains.py` has parsed it (the nth residue), **not the author residue number in the PDB**. Whenever the first residue is not numbered 1, or residues are missing in the middle, the two drift apart.
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

> **04 Designing only specific residues -- the inverse constraint**

**Input**: put a complex PDB in `mpnn_inputs/` (e.g. 3HTN.pdb)  
**What it does**: the inverse of 03 — redesign only the few positions listed and leave everything else untouched

Append `--specify_non_fixed` to `make_fixed_positions_dict.py` and the positions listed in `--position_list` become "the positions to design", with everything else fixed
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

> **05 Custom residue tying -- co-design with tied_positions**

**Input**: put a complex PDB in `mpnn_inputs/` (e.g. 3HTN.pdb)  
**What it does**: tie corresponding positions across several chains so they are forced to the same amino acid

`make_tied_positions_dict.py` writes `tied_pdbs.jsonl`, "tying" corresponding positions across chains; pass it at run time with `--tied_positions_jsonl`. The example below ties positions 1-8 of chains A and B pairwise
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

> **06 Homo-oligomer -- C2/C3/C4 symmetric design**

**Input**: put homo-oligomer PDBs in `mpnn_inputs/` (e.g. 4GYT.pdb, 6EHB.pdb)  
**What it does**: tie the corresponding positions of every chain automatically, yielding a homo-oligomer whose chains carry identical sequences

The chains of a homo-oligomer are the same length and should carry the same sequence; adding `--homooligomer 1` to `make_tied_positions_dict.py` ties every corresponding position across all chains automatically, so each chain in the output has exactly the same sequence
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

> **07 Setting the backbone noise**

**Input**: reuse 5L33.pdb and 6MRR.pdb from 01 directly
**What it does**: add Gaussian noise to the backbone coordinates before designing, for more sequence diversity

Add `--backbone_noise 0.10` at run time (in Å). Note that this is the backbone noise at inference — training carries noise by default anyway (v_48_020, for instance, was trained with 0.2)
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

> **08 Amino-acid bias -- encouraging or suppressing particular residues**

**Input**: reuse 5L33.pdb and 6MRR.pdb from 01 directly 
**What it does**: adjust the amino-acid preference globally

`make_bias_AA.py` writes a global bias table (in log space; positive favours, negative suppresses), passed at run time with `--bias_AA_jsonl`. The example below tilts the whole design towards polar amino acids
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

### Part 2 — Advanced options (each adding exactly one flag to the 01 baseline)

> **09 Globally excluding amino acids -- |general|ban certain AAs|--omit_AAs|**

**Input**: as in 01 (monomer backbone); **What it does**: forbid the named amino acids at every position outright (cysteine C here)

`--omit_AAs` bans amino acids from every position (default `"X"`). The example below bans cysteine C, so no free thiol gets designed in
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
> **When designing antibody CDRs, treat `--omit_AAs "C"` as a default parameter,** an extra free thiol in a CDR competes for mispairing and leads to isoforms, aggregation and batch-to-batch inconsistency.

> **10 Scoring only, no design -- |general|score the original sequence|--score_only 1|**

**Input**: as in 01 (monomer backbone); **What it does**: skip design entirely and just score the native sequence that came with the input PDB

`--score_only 1` skips design and computes score / global_score for the sequence already in the input PDB, saving the result as `.npz` (under `<out>/score_only/`). Use it to judge a native sequence, or to compare before and after design
```bash
python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
  --jsonl_path /data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --out_folder /data/lmk/ProteinMPNN/mpnn_outputs \
  --seed 37 \
  --batch_size 1 \
  --score_only 1
```

> **11 Scoring an external sequence -- |general|score a supplied sequence|--path_to_fasta|**

**Input**: the backbone as in 01, plus the sequence to be scored, `seqs_to_score.fa`, in `mpnn_inputs/`; **What it does**: score an external sequence on that backbone (one that came back from validation)

Together with `--score_only 1`, `--path_to_fasta` scores an external sequence (one returned from AF2/AF3 validation, say) on the same backbone. Inside the FASTA, chains are separated by `/` in alphabetical order
```bash
python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
  --jsonl_path /data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --out_folder /data/lmk/ProteinMPNN/mpnn_outputs \
  --seed 37 \
  --batch_size 1 \
  --score_only 1 \
  --path_to_fasta /data/lmk/ProteinMPNN/mpnn_inputs/seqs_to_score.fa
```

> **12 The soluble model -- |general|suppress surface hydrophobics|--use_soluble_model|**

**Input**: as in 01 (monomer backbone); **What it does**: switch to the soluble weights, which lean towards fewer surface hydrophobic residues when designing soluble proteins

Add `--use_soluble_model` to switch to weights trained on soluble proteins only, reducing exposed hydrophobic residues in soluble protein design (see the `model_name` decoding table in [FASTA_Format](./ProteinMPNN_FASTA_Format_EN.md))
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

> **13 The Cα-only model -- |general|Cα coordinates alone|--ca_only|**

**Input**: as in 01 (monomer backbone, only the Cα coordinates are read); **What it does**: design with the Cα-only weights, for inputs missing backbone atoms or holding only a Cα trace

Add `--ca_only` to switch to the Cα-only weights, which read the Cα coordinates alone — suited to inputs missing backbone atoms or holding only a Cα trace
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

> **14 Switching the vanilla noise tier -- |general|pick a weight tier|--model_name|**

**Input**: as in 01 (monomer backbone); **What it does**: switch the training-noise tier of the vanilla weights, trading robustness against reliance on backbone accuracy

`--model_name` picks the vanilla weight tier (default `v_48_020`). The example below switches to `v_48_002` (0.02 Å training noise, which relies on a high-accuracy backbone and suits crystal structures)
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

> **15 Saving per-position probabilities / scores -- |general|export npz|--save_probs / --save_score|**

**Input**: as in 01 (monomer backbone); **What it does**: export the per-position amino-acid probabilities alongside the design (for conservation / entropy analysis)

`--save_probs 1` saves the probabilities of all 20 amino acids at each position as `.npz` (under `<out>/probs/`), which supports per-position conservation / entropy analysis; `--save_score 1` additionally saves the per-sequence score
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

> **16 Conditional probability analysis -- |general|p(s_i | backbone, rest of sequence)|--conditional_probs_only 1|**

**Input**: as in 01 (monomer backbone); **What it does**: produce no sequences, only the per-position conditional probabilities given "the backbone plus the rest of the sequence"

`--conditional_probs_only 1` does not sample; it outputs only each position's probability given "the backbone plus the rest of the sequence" (saved as `.npz`), for analysing how much a position depends on its context
```bash
python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
  --jsonl_path /data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --out_folder /data/lmk/ProteinMPNN/mpnn_outputs \
  --seed 37 \
  --batch_size 1 \
  --conditional_probs_only 1
```

> **17 Unconditional probability analysis -- |general|p(s_i | backbone)|--unconditional_probs_only 1|**

**Input**: as in 01 (monomer backbone); **What it does**: produce no sequences, only the per-position unconditional probabilities given "the backbone alone"

`--unconditional_probs_only 1` gives each position's probability given "the backbone alone" in a single forward pass — faster than 16, and handy for a quick look at which amino acids the backbone itself prefers
```bash
python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
  --jsonl_path /data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --out_folder /data/lmk/ProteinMPNN/mpnn_outputs \
  --seed 37 \
  --batch_size 1 \
  --unconditional_probs_only 1
```

> **18 Single PDB straight in -- |one file|skip parsing|--pdb_path / --pdb_path_chains|**

**Input**: a single PDB file in `mpnn_inputs/` (5L33.pdb, say — no parsing needed); **What it does**: skip the parsing step and design straight from one PDB

When only one PDB is being designed, `parse_multiple_chains.py` can be skipped: `--pdb_path` passes the file directly and `--pdb_path_chains` names the chains to design (the rest stay fixed)
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

> **19 Sweeping the sampling temperature -- |general|several temperatures in one run|--sampling_temp|**

**Input**: as in 01 (monomer backbone); **What it does**: draw a batch of sequences at each of several sampling temperatures in one run, making diversity easy to compare

`--sampling_temp` accepts several space-separated temperatures and draws a batch at each in one go (higher temperature means more diversity and lower confidence); the `T=` in the FASTA header records which temperature each came from
```bash
python /data/lmk/ProteinMPNN/protein_mpnn_run.py \
  --jsonl_path /data/lmk/ProteinMPNN/mpnn_outputs/parsed_pdbs.jsonl \
  --out_folder /data/lmk/ProteinMPNN/mpnn_outputs \
  --num_seq_per_target 5 \
  --sampling_temp "0.1 0.2 0.3" \
  --seed 37 \
  --batch_size 1
```

##### [ProteinMPNN official documentation](https://github.com/dauparas/ProteinMPNN)
