# SB-OOB-100seed

100-seed simulation pipeline accompanying the paper:

> **Sequential Bootstrap for Out-of-Bag Error Estimation: A 100-Seed Replication Study and Variance-Structure Analysis**
> Cheng Peng. arXiv:2511.18065.

This repository provides the complete code, configuration, and aggregated outputs for reproducing every numerical result in the paper's Section 5 and Appendix A.

---

## Repository layout

```
SB-OOB-100seed/
├── README.md                        this file
├── Unified_funs.R                   bagging/OOB/SB-OOB helper functions
├── Unified_parallel.R               main per-seed simulation driver
├── run.slurm                        SLURM array job (general partition)
├── run_scavenger.slurm              SLURM array job (scavenger partition, with --requeue + idempotency)
├── aggregate_100seeds.R             cross-seed aggregation -> tables/
└── tables/                          paper-ready aggregated outputs
    ├── all_seeds_flat.csv           flat long table (6,200 rows)
    ├── agg_100seeds.csv             cross-seed paired statistics (57 cells × 18 stats)
    ├── summary_by_exp.csv           per-EXP overview
    ├── table_exp1_100seeds.csv
    ├── table_exp2_100seeds.csv
    ├── table_exp3_100seeds.csv
    ├── table_exp4_100seeds.csv
    └── table_exp5_100seeds.csv
```

The **raw per-seed RData outputs** (500 summary files + 500 raw-paired files, ~XX MB compressed) are too large for the main repository. They are attached to the **GitHub Release**:

- Release: **[v1.0-data](https://github.com/Cheng-Peng0718/SB-OOB-100seed/releases/tag/v1.0-data)**
  - `results_100seeds.tar.gz` — all 100 `results/seed_NNN/` directories
  - `logs_100seeds.tar.gz` — all SLURM stdout/stderr logs

Download and extract:

```bash
wget https://github.com/Cheng-Peng0718/SB-OOB-100seed/releases/download/v1.0-data/results_100seeds.tar.gz
tar -xzf results_100seeds.tar.gz
```

---

## Reproducing the paper's tables and figures

### Option A: re-aggregate from pre-existing RData (fast, ~1 minute)

If you have `results/seed_001/ ... seed_100/` (either by downloading the release tarball above, or by re-running the simulation, see Option B):

```bash
Rscript --vanilla aggregate_100seeds.R
```

This recreates `tables/agg_100seeds.csv` and the five per-EXP tables, and prints the cross-experiment summary plus the real-vs-synthetic permutation test referenced in Section 5.1 of the paper.

### Option B: re-run the simulation from scratch

Requires R ≥ 4.3 with `rpart`, `mlbench`, `MASS`, `dplyr`, `foreach`, `doParallel`, `doRNG`. On a single 50-core machine, one seed takes ~2 hours; on a SLURM cluster the full 100-seed array completes in 24–48 hours depending on backfill.

```bash
# single-seed run (for testing)
Rscript --vanilla Unified_parallel.R 1

# full SLURM array job
sbatch run.slurm
# or, on clusters with a scavenger partition:
sbatch run_scavenger.slurm
```

Both SLURM scripts write to `results/seed_NNN/` and are idempotent — already-completed seeds (those with `exp5_raw.RData` present) are skipped on requeue.

---

## File-by-file description of `tables/`

| File | Rows | What it contains |
|---|---:|---|
| `all_seeds_flat.csv` | 6,200 | One row per (seed, experiment, dataset, metric). Within-seed paired statistics computed across 50 internal replications. |
| `agg_100seeds.csv` | 62 | One row per (experiment, dataset, metric). Cross-seed paired statistics (Wilcoxon, paired-t, sign test, Pitman–Morgan, Cohen's d_z, sd-ratio, 95 % CI). **This is the master table for paper Section 5.** |
| `table_exp{1..5}_100seeds.csv` | 5–18 each | The same table partitioned by experiment family, exactly as used in Tables 2–7 of the paper. |
| `summary_by_exp.csv` | 5 | Per-EXP counts of significant cells (mean-level and variance-level), broken down by sign. Corresponds to Table 8 of the paper. |

---

## Method summary (for reviewers)

For each of 100 independent random seeds and each of 12 datasets, we:

1. Construct two bagged CART ensembles using identical hyperparameters, base learners, and aggregation logic; the **only** difference is the resampling scheme:
   - **OOB**: classical multinomial bootstrap (`bag_cart_OOB` in `Unified_funs.R`).
   - **SB-OOB**: Sequential Bootstrap with target distinct-count `k_n = ⌈0.632·n⌉` (`bag_cart_SB`).
2. Compute the 5 experimental-family diagnostics of Breiman (1996b): node-level classification accuracy (E1, E2), node-level regression accuracy (EB1, EB2), within-node stability (R1–R4), OOB-vs-test alignment (absdiff, eOB, eTS, ratio), and meta-prediction MSE.
3. Repeat 50 times within each seed (independent data generation or train-test resampling), so each (seed, dataset, metric) cell has 50 paired (OOB, SB-OOB) observations.
4. Aggregate up two levels: **within-seed** paired tests (in `Unified_parallel.R`'s `.summarise_runs`), then **cross-seed** paired tests (in `aggregate_100seeds.R`).

Cross-seed paired statistical comparison across n=100 seeds is the basis for every claim in Section 5 of the paper.

---

## Citation

If you use this code or data, please cite:

```bibtex
@article{peng2026sboob,
  title  = {Sequential Bootstrap for Out-of-Bag Error Estimation:
            A 100-Seed Replication Study and Variance-Structure Analysis},
  author = {Peng, Cheng},
  journal= {arXiv preprint arXiv:2511.18065},
  year   = {2026}
}
```

## License

MIT.
