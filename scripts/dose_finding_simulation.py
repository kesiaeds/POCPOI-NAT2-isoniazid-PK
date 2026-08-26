"""Monte Carlo dose-finding simulation off the n=104 popPK model (M3 fit).

For each NAT2 phenotype, simulate AUC0-24 and C24 across a dose grid
(200-1300 mg, 50 mg steps) and identify the dose that best reproduces the
intermediate-metabolizer flat-900 reference exposure while keeping the
proportion of subjects with supratherapeutic trough (C24 > 0.15 µg/mL) low.

Approach
  - Parameters: typical values from sherlock/output/popPK_parameters.csv.
  - Per-subject variability: lognormal η on CL with ω = sqrt(ln(1 + CV^2))
    derived from the reported BSV(CV%); weight resampled from the observed
    PK cohort within each phenotype; residual error not added (we focus on
    structural exposure, not the assay).
  - Concentration / AUC computed analytically from the 2-compartment
    first-order absorption solution -- vectorized over virtual subjects.

Outputs (analysis/output/simulation/):
  dose_finding_results.csv       per-phenotype x dose summary stats
  dose_finding_recommended.csv   recommended dose per phenotype
analysis/output/figures/:
  Figure_dose_finding.{png,pdf}  A: dose vs AUC band, B: dose vs %supra
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path('/Users/kesiaeds/repos/POCPOI')
DATA = ROOT / 'analysis' / 'data'
ART  = ROOT / 'analysis' / 'output'
FIG  = ART / 'figures'
SIM  = ART / 'simulation'
SIM.mkdir(parents=True, exist_ok=True)

# ---------- popPK parameters (Sherlock M3 fit, n=104) ----------------------
pars = pd.read_csv(ROOT / 'sherlock' / 'output' / 'popPK_parameters.csv',
                   index_col=0)
tKA = pars.loc['tKA', 'Estimate']
tCL = pars.loc['tCL', 'Estimate']
tV2 = pars.loc['tV2', 'Estimate']
tQ  = pars.loc['tQ',  'Estimate']
tV3 = pars.loc['tV3', 'Estimate']
dCL = {'Slow':         pars.loc['dCldACEPROFILE1', 'Estimate'],
       'Intermediate': 0.0,
       'Rapid':        pars.loc['dCldACEPROFILE3', 'Estimate']}
bsv_cv = pars.loc['eta.cl', 'BSV(CV%)'] / 100.0
omega  = np.sqrt(np.log1p(bsv_cv ** 2))

EXP_CL_WT = 0.75   # dCldWT fixed
EXP_V_WT  = 1.0    # dVdWT fixed (applied to V2 and V3)
REF_WT    = 70.0
SUPRA     = 0.15   # µg/mL
SUPRA_CAP = 15.0   # % of simulated subjects allowed above SUPRA (C24-ceiling logic)

PHEN_ORDER  = ['Slow', 'Intermediate', 'Rapid']
PHEN_ID     = {'Slow': 1, 'Intermediate': 2, 'Rapid': 3}
PHEN_COLOUR = {'Slow': '#1f77b4', 'Intermediate': '#2ca02c', 'Rapid': '#d62728'}

# ---------- reference targets: intermediate × flat 900 mg ------------------
# Two clinically-motivated targets are computed:
#   AUC target  -- reproduce intermediate-metabolizer exposure (AUC0-24)
#   C24 target  -- match intermediate-metabolizer trough; more directly tied
#                  to hepatotoxicity risk under weekly 3HP for LTBI, where
#                  efficacy is not exposure-limited.
exps = pd.read_csv(ART / 'exposures_subject_occasion.csv')
ref_im_flat = exps[(exps['Phenotype'] == 'Intermediate') &
                   (exps['Dosing']    == 'Flat 900 mg')]
ref_auc = ref_im_flat['AUC24'].dropna()
ref_c24 = ref_im_flat['C24'].dropna()

AUC_TARGET_MED  = ref_auc.median()
AUC_TARGET_LO   = ref_auc.quantile(0.25)
AUC_TARGET_HI   = ref_auc.quantile(0.75)

C24_TARGET_MED  = ref_c24.median()
C24_TARGET_LO   = ref_c24.quantile(0.25)
C24_TARGET_HI   = ref_c24.quantile(0.75)

# ---------- weight pool (PK cohort, per phenotype) -------------------------
nm = pd.read_csv(DATA / 'pk_nonmem_combined.csv')
wt_pool = nm.drop_duplicates('ID')[['ACE_PROFILE', 'WT']]

# ---------- analytical 2-cmt FO absorption ---------------------------------
def _hybrid(cl, v2, q, v3):
    k20 = cl / v2
    k23 = q  / v2
    k32 = q  / v3
    s = k20 + k23 + k32
    disc = np.sqrt(s * s - 4 * k20 * k32)
    return k32, 0.5 * (s + disc), 0.5 * (s - disc)

def _auc24(dose, ka, cl, v2, q, v3):
    k32, a, b = _hybrid(cl, v2, q, v3)
    pre = ka * dose / v2
    ta = (k32 - a)  / ((ka - a) * (b - a))  * (1 - np.exp(-24 * a))  / a
    tb = (k32 - b)  / ((ka - b) * (a - b))  * (1 - np.exp(-24 * b))  / b
    tk = (k32 - ka) / ((a - ka) * (b - ka)) * (1 - np.exp(-24 * ka)) / ka
    return pre * (ta + tb + tk)

def _conc(t, dose, ka, cl, v2, q, v3):
    k32, a, b = _hybrid(cl, v2, q, v3)
    pre = ka * dose / v2
    ta = (k32 - a)  / ((ka - a) * (b - a))  * np.exp(-a  * t)
    tb = (k32 - b)  / ((ka - b) * (a - b))  * np.exp(-b  * t)
    tk = (k32 - ka) / ((a - ka) * (b - ka)) * np.exp(-ka * t)
    return pre * (ta + tb + tk)

# ---------- simulation -----------------------------------------------------
N_SUBJ = 2000
DOSES  = np.arange(200, 1801, 50)
rng    = np.random.default_rng(20260601)

records = []
for phen in PHEN_ORDER:
    pool = wt_pool[wt_pool['ACE_PROFILE'] == PHEN_ID[phen]]['WT'].values
    wts  = rng.choice(pool, size=N_SUBJ, replace=True)
    eta  = rng.normal(0.0, omega, size=N_SUBJ)
    cl_i = tCL * (wts / REF_WT) ** EXP_CL_WT * np.exp(dCL[phen] + eta)
    v2_i = tV2 * (wts / REF_WT) ** EXP_V_WT
    v3_i = tV3 * (wts / REF_WT) ** EXP_V_WT
    q_i  = np.full(N_SUBJ, tQ)
    for dose in DOSES:
        auc = _auc24(dose, tKA, cl_i, v2_i, q_i, v3_i)
        c24 = _conc(24.0, dose, tKA, cl_i, v2_i, q_i, v3_i)
        records.append({
            'Phenotype':    phen,
            'Dose_mg':      int(dose),
            'AUC_med':      float(np.median(auc)),
            'AUC_q1':       float(np.quantile(auc, 0.25)),
            'AUC_q3':       float(np.quantile(auc, 0.75)),
            'C24_med':      float(np.median(c24)),
            'C24_q1':       float(np.quantile(c24, 0.25)),
            'C24_q3':       float(np.quantile(c24, 0.75)),
            'pct_in_auc_target':  float(np.mean((auc >= AUC_TARGET_LO) &
                                                (auc <= AUC_TARGET_HI)) * 100),
            'pct_in_c24_target':  float(np.mean((c24 >= C24_TARGET_LO) &
                                                (c24 <= C24_TARGET_HI)) * 100),
            'pct_in_dual_target': float(np.mean((auc >= AUC_TARGET_LO) &
                                                (auc <= AUC_TARGET_HI) &
                                                (c24 >= C24_TARGET_LO) &
                                                (c24 <= C24_TARGET_HI)) * 100),
            'pct_supra':    float(np.mean(c24 > SUPRA) * 100),
        })

sim = pd.DataFrame(records)
sim.to_csv(SIM / 'dose_finding_results.csv', index=False)

# ---------- recommended dose per phenotype ---------------------------------
# Primary criterion: minimise |median(metric) - target median|.
# Tiebreaker: prefer doses with lower %supra.
def _pick(phen, metric_col, target_med):
    sub = sim[sim['Phenotype'] == phen].copy()
    sub['gap'] = (sub[metric_col] - target_med).abs()
    return sub.sort_values(['gap', 'pct_supra']).iloc[0]

opt_auc = (pd.DataFrame([_pick(p, 'AUC_med', AUC_TARGET_MED) for p in PHEN_ORDER])
             [['Phenotype', 'Dose_mg', 'AUC_med', 'C24_med',
               'pct_in_auc_target', 'pct_supra']]
             .rename(columns={'Dose_mg': 'AUC_match_dose_mg'}))
opt_c24 = (pd.DataFrame([_pick(p, 'C24_med', C24_TARGET_MED) for p in PHEN_ORDER])
             [['Phenotype', 'Dose_mg', 'AUC_med', 'C24_med',
               'pct_in_c24_target', 'pct_supra']]
             .rename(columns={'Dose_mg': 'C24_match_dose_mg'}))

# C24-ceiling: highest dose whose %supra stays at or below SUPRA_CAP.
# If even the lowest dose exceeds the cap, fall back to that lowest dose
# (and flag it via pct_supra_at_ceiling so the manuscript shows the violation).
def _pick_ceiling(phen):
    sub = sim[sim['Phenotype'] == phen].sort_values('Dose_mg')
    ok  = sub[sub['pct_supra'] <= SUPRA_CAP]
    row = ok.iloc[-1] if len(ok) else sub.iloc[0]
    return row

opt_ceil = (pd.DataFrame([_pick_ceiling(p) for p in PHEN_ORDER])
              [['Phenotype', 'Dose_mg', 'AUC_med', 'C24_med', 'pct_supra']]
              .rename(columns={'Dose_mg':   'Ceiling_dose_mg',
                               'pct_supra': 'pct_supra_at_ceiling'}))

opt = (opt_auc
       .merge(opt_c24[['Phenotype', 'C24_match_dose_mg', 'pct_in_c24_target']],
              on='Phenotype')
       .merge(opt_ceil[['Phenotype', 'Ceiling_dose_mg', 'pct_supra_at_ceiling']],
              on='Phenotype'))
opt.to_csv(SIM / 'dose_finding_recommended.csv', index=False)

print('Reference targets (Intermediate × Flat 900 mg)')
print(f'  AUC0-24 : median {AUC_TARGET_MED:.1f} '
      f'(IQR {AUC_TARGET_LO:.1f}-{AUC_TARGET_HI:.1f}) mg·h/L')
print(f'  C24     : median {C24_TARGET_MED:.4f} '
      f'(IQR {C24_TARGET_LO:.4f}-{C24_TARGET_HI:.4f}) µg/mL')
print()
print('Recommended dose per phenotype (both targets):')
print(opt.to_string(index=False))

# ---------- figure ---------------------------------------------------------
fig, axes = plt.subplots(2, 3, figsize=(16, 9.4))
(axA, axB, axC), (axD, axE, axF) = axes

for phen in PHEN_ORDER:
    sub = sim[sim['Phenotype'] == phen]
    c   = PHEN_COLOUR[phen]
    for ax in (axA, axD):
        ax.plot(sub['Dose_mg'], sub['AUC_med'], '-o', color=c, label=phen,
                linewidth=1.8, markersize=4)
        ax.fill_between(sub['Dose_mg'], sub['AUC_q1'], sub['AUC_q3'],
                        color=c, alpha=0.18)
    for ax in (axB, axE):
        ax.plot(sub['Dose_mg'], sub['C24_med'], '-o', color=c, label=phen,
                linewidth=1.8, markersize=4)
        ax.fill_between(sub['Dose_mg'], sub['C24_q1'], sub['C24_q3'],
                        color=c, alpha=0.18)
    for ax in (axC, axF):
        ax.plot(sub['Dose_mg'], sub['pct_supra'], '-o', color=c, label=phen,
                linewidth=1.8, markersize=4)

# Panel A — AUC vs dose (full range)
axA.axhspan(AUC_TARGET_LO, AUC_TARGET_HI, color='grey', alpha=0.25,
            label=f'AUC target IQR ({AUC_TARGET_LO:.0f}-{AUC_TARGET_HI:.0f})')
axA.axhline(AUC_TARGET_MED, color='black', linestyle='--', linewidth=1,
            alpha=0.6, label=f'AUC target median ({AUC_TARGET_MED:.0f})')
axA.set_xlabel('Dose (mg)')
axA.set_ylabel('Simulated AUC$_{0-24}$ (mg·h/L)')
axA.legend(loc='upper left', fontsize=8)
axA.grid(axis='y', alpha=.3)
axA.text(-0.10, 1.02, 'A', transform=axA.transAxes,
         fontsize=14, fontweight='bold', va='bottom', ha='left')

# Panel B — C24 vs dose (log scale, full range)
axB.axhspan(C24_TARGET_LO, C24_TARGET_HI, color='grey', alpha=0.25,
            label=f'C24 target IQR ({C24_TARGET_LO:.3f}-{C24_TARGET_HI:.3f})')
axB.axhline(C24_TARGET_MED, color='black', linestyle='--', linewidth=1,
            alpha=0.6, label=f'C24 target median ({C24_TARGET_MED:.3f})')
axB.axhline(SUPRA, color='red', linestyle='--', linewidth=1, alpha=0.6,
            label=f'Supratherapeutic ({SUPRA} µg/mL)')
axB.set_xlabel('Dose (mg)')
axB.set_ylabel('Simulated C$_{24}$ (µg/mL)')
axB.set_yscale('log')
axB.legend(loc='upper left', fontsize=8)
axB.grid(axis='y', alpha=.3)
axB.text(-0.10, 1.02, 'B', transform=axB.transAxes,
         fontsize=14, fontweight='bold', va='bottom', ha='left')

# Panel C — %supra vs dose (full range)
axC.axhline(SUPRA_CAP, color='red', linestyle='--', linewidth=1, alpha=0.6,
            label=f'{SUPRA_CAP:.0f}% supratherapeutic cap')
axC.set_xlabel('Dose (mg)')
axC.set_ylabel(f'% subjects with C$_{{24}}$ > {SUPRA} µg/mL')
axC.legend(loc='upper left', fontsize=8)
axC.grid(axis='y', alpha=.3)
axC.text(-0.10, 1.02, 'C', transform=axC.transAxes,
         fontsize=14, fontweight='bold', va='bottom', ha='left')

# Panel D — AUC zoom around target IQR
auc_pad = 0.2 * (AUC_TARGET_HI - AUC_TARGET_LO)
axD.axhspan(AUC_TARGET_LO, AUC_TARGET_HI, color='grey', alpha=0.25,
            label=f'AUC target IQR ({AUC_TARGET_LO:.0f}-{AUC_TARGET_HI:.0f})')
axD.axhline(AUC_TARGET_MED, color='black', linestyle='--', linewidth=1,
            alpha=0.6, label=f'AUC target median ({AUC_TARGET_MED:.0f})')
axD.set_xlabel('Dose (mg)')
axD.set_ylabel('Simulated AUC$_{0-24}$ (mg·h/L)')
axD.set_ylim(AUC_TARGET_LO - auc_pad, AUC_TARGET_HI + auc_pad)
axD.legend(loc='upper left', fontsize=8)
axD.grid(axis='y', alpha=.3)
axD.text(-0.10, 1.02, 'D', transform=axD.transAxes,
         fontsize=14, fontweight='bold', va='bottom', ha='left')
axD.set_title('Zoom: AUC target IQR', fontsize=10)

# Panel E — C24 zoom around target IQR (linear, narrow window)
c24_pad_lo = 0.6 * C24_TARGET_LO     # show a bit below LO
c24_top    = max(C24_TARGET_HI * 2, SUPRA * 1.05)
axE.axhspan(C24_TARGET_LO, C24_TARGET_HI, color='grey', alpha=0.25,
            label=f'C24 target IQR ({C24_TARGET_LO:.3f}-{C24_TARGET_HI:.3f})')
axE.axhline(C24_TARGET_MED, color='black', linestyle='--', linewidth=1,
            alpha=0.6, label=f'C24 target median ({C24_TARGET_MED:.3f})')
axE.axhline(SUPRA, color='red', linestyle='--', linewidth=1, alpha=0.6,
            label=f'Supratherapeutic ({SUPRA} µg/mL)')
axE.set_xlabel('Dose (mg)')
axE.set_ylabel('Simulated C$_{24}$ (µg/mL)')
axE.set_ylim(c24_pad_lo, c24_top)
axE.legend(loc='upper left', fontsize=8)
axE.grid(axis='y', alpha=.3)
axE.text(-0.10, 1.02, 'E', transform=axE.transAxes,
         fontsize=14, fontweight='bold', va='bottom', ha='left')
axE.set_title('Zoom: C$_{24}$ target IQR', fontsize=10)

# Panel F — %supra zoom (0-20%)
axF.axhline(SUPRA_CAP, color='red', linestyle='--', linewidth=1, alpha=0.6,
            label=f'{SUPRA_CAP:.0f}% supratherapeutic cap')
axF.set_xlabel('Dose (mg)')
axF.set_ylabel(f'% subjects with C$_{{24}}$ > {SUPRA} µg/mL')
axF.set_ylim(0, 20)
axF.legend(loc='upper left', fontsize=8)
axF.grid(axis='y', alpha=.3)
axF.text(-0.10, 1.02, 'F', transform=axF.transAxes,
         fontsize=14, fontweight='bold', va='bottom', ha='left')
axF.set_title('Zoom: 0-20% supratherapeutic', fontsize=10)

fig.tight_layout()
fig.savefig(FIG / 'Figure_dose_finding.png', dpi=300, bbox_inches='tight')
fig.savefig(FIG / 'Figure_dose_finding.pdf', bbox_inches='tight')
plt.close(fig)
print('Wrote', SIM, 'and Figure_dose_finding.{png,pdf}')

# ---------- target attainment figure --------------------------------------
# Three panels: %subjects in AUC IQR, in C24 IQR, in BOTH IQRs (dual target)
PROTOCOL = {'Slow': 300, 'Intermediate': 900, 'Rapid': 1500}

fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), sharey=True)
panel_specs = [
    (axes[0], 'pct_in_auc_target',
     f'AUC IQR ({AUC_TARGET_LO:.0f}-{AUC_TARGET_HI:.0f})',
     'A', 'In AUC$_{0-24}$ target IQR'),
    (axes[1], 'pct_in_c24_target',
     f'C24 IQR ({C24_TARGET_LO:.3f}-{C24_TARGET_HI:.3f})',
     'B', 'In C$_{24}$ target IQR'),
    (axes[2], 'pct_in_dual_target',
     'Both IQRs simultaneously', 'C', 'In both target IQRs (dual)'),
]

for ax, col, sub_title, label, ylab in panel_specs:
    for phen in PHEN_ORDER:
        sub = sim[sim['Phenotype'] == phen]
        ax.plot(sub['Dose_mg'], sub[col], '-o', color=PHEN_COLOUR[phen],
                label=phen, linewidth=1.8, markersize=4)
        ax.axvline(PROTOCOL[phen], color=PHEN_COLOUR[phen], linestyle=':',
                   linewidth=1, alpha=0.7)
    ax.axhline(50, color='grey', linestyle='--', linewidth=1, alpha=0.5,
               label='50% attainment')
    ax.set_xlabel('Dose (mg)')
    ax.set_ylabel(f'% subjects {ylab}')
    ax.set_ylim(0, 100)
    ax.grid(axis='y', alpha=.3)
    ax.set_title(sub_title, fontsize=10)
    ax.legend(loc='upper right', fontsize=8)
    ax.text(-0.10, 1.02, label, transform=ax.transAxes,
            fontsize=14, fontweight='bold', va='bottom', ha='left')

fig.suptitle('Target attainment vs dose '
             '(vertical dotted lines = protocol phenotype doses '
             '300 / 900 / 1500 mg)', fontsize=10.5, y=1.02)
fig.tight_layout()
fig.savefig(FIG / 'Figure_dose_finding_target_attainment.png',
            dpi=300, bbox_inches='tight')
fig.savefig(FIG / 'Figure_dose_finding_target_attainment.pdf',
            bbox_inches='tight')
plt.close(fig)
print('Wrote Figure_dose_finding_target_attainment.{png,pdf}')
