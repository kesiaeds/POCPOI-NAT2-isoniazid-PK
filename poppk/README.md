# POCPOI popPK on Sherlock — step-by-step

Goal: fit the combined-cohort (CG + Manaus) two-compartment popPK model in
nlmixr2 on Sherlock and copy the artifacts back to this laptop.

The whole flow is **transfer → install (once, ~1 h) → fit (~30-60 min) → pull
results back**. Commands below use copy-pasteable blocks.

## Sherlock account assumed

You are already logged in via `ssh kesiaeds@login.sherlock.stanford.edu` in
your other window.

---

## 0. One-time choices

```bash
# On your laptop, edit POCPOI_WORKDIR if your scratch path differs.
export POCPOI_WORKDIR=/scratch/users/kesiaeds/POCPOI
```

## 1. Stage and transfer the package

From this laptop, in `/Users/kesiaeds/repos/POCPOI/`:

```bash
# Pack everything Sherlock needs
cd /Users/kesiaeds/repos/POCPOI
tar -czf sherlock_package.tar.gz -C sherlock data scripts

# Ship it
scp sherlock_package.tar.gz \
    kesiaeds@login.sherlock.stanford.edu:/scratch/users/kesiaeds/
```

Then, on the Sherlock login shell:

```bash
export POCPOI_WORKDIR=/scratch/users/$USER/POCPOI
mkdir -p "$POCPOI_WORKDIR"
cd "$POCPOI_WORKDIR"
tar -xzf /scratch/users/$USER/sherlock_package.tar.gz
mkdir -p output Rlib/4.3
ls
# expected: data/  scripts/  output/  Rlib/
```

## 2. Install nlmixr2 + dependencies (do this ONCE per Sherlock account)

`rxode2` compiles C++ for every model; the install is the slow part.
Run on an `sdev` node, NOT the login node:

```bash
# Inside Sherlock
sdev -t 03:00:00 -c 8 -m 32G
ml R/4.3 gcc/12 cmake openblas
cd "$POCPOI_WORKDIR"
export R_LIBS_USER="$POCPOI_WORKDIR/Rlib/4.3"
Rscript scripts/install_packages.R 2>&1 | tee output/install.log
exit   # leave sdev
```

If the sdev session times out before the install finishes, just re-enter
`sdev` and re-run — `install.packages` skips anything already installed.

Alternatively, submit it as a batch job (logs land in `output/install_*.log`):

```bash
cd "$POCPOI_WORKDIR"
export POCPOI_WORKDIR=/scratch/users/$USER/POCPOI
sbatch --export=POCPOI_WORKDIR scripts/install_packages.sbatch
```

Confirm everything loads:

```bash
ml R/4.3 gcc/12 openblas
export R_LIBS_USER="$POCPOI_WORKDIR/Rlib/4.3"
Rscript -e 'suppressPackageStartupMessages(library(nlmixr2)); cat("nlmixr2 OK\n")'
```

## 3. Submit the popPK fit

```bash
cd "$POCPOI_WORKDIR"
export POCPOI_WORKDIR=/scratch/users/$USER/POCPOI
sbatch --export=POCPOI_WORKDIR scripts/run_poppk.sbatch
squeue -u $USER
```

The fit writes:
- `output/fit_combined_focei.rds` — saved model object
- `output/popPK_parameters.csv` — typical-value table + RSE
- `output/individual_CL.csv` — empirical-Bayes CL per subject
- `output/exposures.csv` — NCA-style AUC/Cmax/C24 per subject-occasion
- `output/exposure_summary.csv` — by-phenotype summary
- `output/figures/vpc.pdf`
- `output/figures/profiles_by_phenotype.pdf`

Run-time estimate: 20-60 min depending on convergence.

## 4. Pull results back to the laptop

From this laptop, when the job is done:

```bash
cd /Users/kesiaeds/repos/POCPOI
mkdir -p sherlock/output
scp -r kesiaeds@login.sherlock.stanford.edu:/scratch/users/kesiaeds/POCPOI/output/ \
       sherlock/
```

Then tell me to update `POCPOI_Results.docx` with the new popPK parameters
and append the popPK paragraph + VPC figure.

---

## Quick troubleshooting

| Symptom | Fix |
|---|---|
| `Error: could not find function "nlmixr2"` after install | `export R_LIBS_USER="$POCPOI_WORKDIR/Rlib/4.3"` before launching R |
| `cc1plus: error: '-march=native'` etc. | Module mismatch — re-load `gcc/12` |
| Job dies with OOM | Bump `--mem=64G` in `run_poppk.sbatch` |
| FOCEI fails to converge | Loosen `eta.cl` initial estimate or switch to `saem` (`est = "saem"` in `pk_analysis_combined.R`) |
| Long compile blocks sdev session | Submit `install_packages.sbatch` as a real batch job |

---

## Files in this package

```
sherlock/
├── data/
│   ├── pk_nonmem_combined.csv   # 1033 rows, 105 IDs, CG+MA combined
│   └── pk_subset.csv            # 105 PK subject-level demographics
├── scripts/
│   ├── install_packages.R       # one-time package install (R)
│   ├── install_packages.sbatch  # SLURM job for #install (~1 h)
│   ├── pk_analysis_combined.R   # the popPK fit
│   └── run_poppk.sbatch         # SLURM job for #fit (~30-60 min)
└── output/                      # created on Sherlock; results land here
```
