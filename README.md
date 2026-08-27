# POCPOI — Genotype-guided isoniazid dosing harmonizes drug exposure in 3HP tuberculosis preventive therapy

**ClinicalTrials.gov:** [NCT05413551](https://clinicaltrials.gov/study/NCT05413551)

This repository contains the analysis code, anonymized participant data, and figures associated with:

> da Silva et al. *Genotype-guided isoniazid dosing harmonizes drug exposure in 3HP tuberculosis preventive therapy*. Under review, 2026.

---

## Study overview

A prospective pharmacokinetic (PK) study enrolling adults with HIV on dolutegravir-based ART initiating isoniazid preventive therapy (IPT) at two sites in Brazil:

- **Campo Grande (CG):** 163 participants
- **Manaus (MA):** 65 participants (all PLHIV)

Participants were genotyped for NAT2 (7-SNP panel: G191A, C282T, T341C, C481T, G590A, A803G, G857A) and classified as Slow, Intermediate, or Rapid acetylators. PK sampling was performed on Day 7 (standard dose, 900 mg guided) and Day 14 (flat 900 mg). A population PK model with Michaelis-Menten elimination was fitted using FOCEI (nlmixr2est v5.0.2, n=104 participants, 816 observations).

---

## Repository structure

```
data/
  all_participants.csv        Demographic and genotyping data, all enrolled (n=228)
  pk_subset.csv               Participants with PK data (n=104)
  pk_nonmem_combined.csv      Full PK dataset (concentration-time records, n=104, 816 obs)

scripts/
  harmonize_data.py                    Harmonize raw CG + Manaus source files into analysis-ready CSVs
  table1_demographics.py               Table 1 (demographics by phenotype and site)
  exposures_descriptive.py             NCA-based exposure summaries (AUC, Cmax, C24) → output/exposures_subject_occasion.csv
  genotyping_qc.py                     NAT2 genotyping QC (qPCR vs. nanopore concordance)
  build_figures.py                     Figures 1–4 and Figure S3
  build_consort_fig.py                 CONSORT/participant flow diagram
  build_structural_model_fig.py        Two-compartment MM structural model schematic
  dose_finding_simulation_empirical.py Empirical CL-based dose-finding (exploratory)
  dose_finding_simulation_proportional.py  Dose-proportional exploratory analysis
  dose_finding_bootstrap.py            Bootstrap uncertainty for NCA-based dose recommendations
  pk_analysis_combined.R               NCA summary statistics for PK exposures
  compare_nca_vs_model.py              Comparison of NCA vs. MM model-derived exposures

poppk/
  pk_analysis_MM.R            MM popPK model fitting (2-cmt + Michaelis-Menten, FOCEI) → produces fit_Km10.rds
  pk_analysis_combined.R      Combined linear popPK fit (for comparison; not primary model)
  vpc_gof.R                   Visual predictive check and goodness-of-fit diagnostics → Figures S1, S2
  dose_finding_simulation_MM.R  Monte Carlo dose-finding using MM model → Figure 5
  exposures_mm_model.R        Model-derived AUC/Cmax/C24 from individual EBE parameters → Figure S4
  run_poppk.sbatch            SLURM job script for popPK fit (Stanford Sherlock HPC)
  vpc_gof.sbatch              SLURM job script for VPC/GOF
  install_packages.R          R package installation for HPC environment
  README.md                   HPC setup and execution instructions

output/
  exposures_subject_occasion.csv  NCA-derived per-subject-occasion AUC, Cmax, C24 (input to figures)
  figures/
    Figure1.pdf/png           INH conc-time profiles by dosing strategy and NAT2 phenotype
    Figure2.pdf/png           AUC0-24 by NAT2 phenotype and dosing strategy
    Figure3.pdf/png           Within-subject paired AUC0-24 (slopegraph)
    Figure4.pdf/png           C24 by NAT2 phenotype and dosing strategy
    Figure5.pdf/png           Phenotype-specific dose recommendations (MM Monte Carlo simulation)
    Figure_S1.pdf/png         Visual predictive check by NAT2 phenotype (poppk/vpc_gof.R)
    Figure_S2.pdf/png         Goodness-of-fit diagnostics (poppk/vpc_gof.R)
    Figure_S3.pdf/png         Individual INH concentration-time profiles (scripts/build_figures.py)
    Figure_S4.pdf/png         Model-derived exposures by phenotype (poppk/exposures_mm_model.R)
  simulation/                 Dose-finding simulation results (CSV)
```

---

## Data dictionary

### all_participants.csv / pk_subset.csv

| Column | Description |
|---|---|
| ID | Pseudonymized participant ID (prefix: CG_ or MA_) |
| SITE | Study site: CG = Campo Grande, MA = Manaus |
| AGE_GROUP | Age group in years: 18-29, 30-39, 40-49, 50-59, 60-69, 70+ |
| SEX | M = Male, F = Female |
| RACE | Self-reported: Mixed/Pardo, White/Branco, Black/Preta, Yellow/Amarela, Indigenous/Indigena |
| WT | Body weight (kg) |
| ACE_PROFILE | NAT2 phenotype: Slow, Intermediate, Rapid |
| STD_DOSE | Standard (Day 7) isoniazid dose (mg) |
| MOD_DOSE | Modified (Day 14) dose (mg) |
| HAS_PK | Whether participant has PK concentration data |

### pk_nonmem_combined.csv

| Column | Description |
|---|---|
| ID | Pseudonymized participant ID |
| SITE | Study site |
| OCC | Occasion: 1 = Day 7, 2 = Day 14 |
| DAY | Nominal sampling day (7 or 14) |
| TIME | Accumulated time from first dose (h); Day 14 records start at 168 h |
| TAD | Time after dose for the current occasion (h) |
| EVID | NONMEM event ID: 1 = dose, 0 = observation |
| MDV | Missing DV flag (1 for dose records) |
| AMT | Dose amount (mg); populated for dose records |
| DV | Isoniazid plasma concentration (µg/mL) |
| CENS | Below-LOQ flag (1 = censored) |
| DOSE | Weight-normalized dose (mg/kg) |
| ACE_PROFILE | NAT2 phenotype |
| AGE_GROUP | Age group |
| SEX | M/F |
| RACE | Self-reported race/ethnicity |
| WT | Body weight (kg) |

---

## Dose recommendations

Based on the population PK model targeting AUC₀₋₂₄ ≥ 10.5 µg·h/mL (52nd percentile, consistent with clinical outcomes literature):

| NAT2 phenotype | Recommended dose |
|---|---|
| Slow acetylator | 600 mg |
| Intermediate acetylator | 900 mg |
| Rapid acetylator | 1,200 mg |

---

## Software requirements

**Python** (≥ 3.9): `pandas`, `numpy`, `matplotlib`, `seaborn`, `scipy`, `python-docx`

**R** (≥ 4.3): `nlmixr2est` (≥ 5.0.2), `rxode2`, `tidyvpc`, `ggplot2`, `dplyr`, `tidyr`, `vpc`

---

## Data availability and ethics

Deidentified participant data (demographics, genotyping, and pharmacokinetic concentration-time records) and all analysis code are provided in this repository. Age is reported in 10-year groups; exact dates, names, and other direct identifiers are not included. The population pharmacokinetic model fit object (`fit_Km10.rds`) is not included due to file size but is available upon reasonable request to the corresponding author (kesiaeds@stanford.edu).

The study was approved by the Research Ethics Committees at UFMS (Campo Grande) and FMT-HVD (Manaus), Brazil, and at Stanford University. All participants provided written informed consent.

---

## Citation

da Silva K.E., et al. *Genotype-guided isoniazid dosing harmonizes drug exposure in 3HP tuberculosis preventive therapy.* Under review, 2026.

---

## Contact

Kesia Esber da Silva — kesiaeds@stanford.edu  
Stanford University, Division of Infectious Diseases
