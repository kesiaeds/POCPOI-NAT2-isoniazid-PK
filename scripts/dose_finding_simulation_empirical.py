"""Empirical-anchored dose-finding simulation (sister to dose_finding_simulation.py).

Re-anchors the virtual cohort to the OBSERVED data instead of popPK typical
values. For each phenotype, draws (CL, WT) pairs with replacement from the
flat 900 mg arm, where CL is computed empirically as 900 mg / observed AUC24.
This bypasses the popPK shrinkage/bias toward typical CL (which over-predicts
clearance for the intermediate-metabolizer reference arm by ~30%).

Other PK parameters (KA, V2, Q, V3) still use popPK typical values because
no reliable subject-level estimates exist for them.

Targets remain the observed intermediate × flat 900 mg AUC and C24 medians,
so by construction the simulation reproduces the reference arm at 900 mg.

Outputs (analysis/output/simulation/):
  dose_finding_empirical_results.csv
  dose_finding_empirical_recommended.csv
analysis/output/figures/:
  Figure_dose_finding_empirical.{png,pdf}        (6 panels, same layout as
                                                   Figure_dose_finding.png)
  Figure_dose_finding_empirical_attainment.{png,pdf}  (3 panels, IQR attainment)
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

# ---------- popPK fixed structural params (still needed for KA/V/Q) -------
pars = pd.read_csv(ROOT / 'sherlock' / 'output' / 'popPK_parameters.csv',
                   index_col=0)
tKA = pars.loc['tKA', 'Estimate']
tV2 = pars.loc['tV2', 'Estimate']
tQ  = pars.loc['tQ',  'Estimate']
tV3 = pars.loc['tV3', 'Estimate']

EXP_V_WT = 1.0
REF_WT   = 70.0
SUPRA    = 0.15
SUPRA_CAP = 15.0

PHEN_ORDER  = ['Slow', 'Intermediate', 'Rapid']
PHEN_COLOUR = {'Slow': '#0072B2', 'Intermediate': '#009E73', 'Rapid': '#D55E00'}

# ---------- reference targets and per-phenotype empirical CL pool ---------
# Anchor: flat 900 mg arm (everyone got the same dose so CL_emp = 900 / AUC24)
exps = pd.read_csv(ART / 'exposures_subject_occasion.csv')
flat = exps[(exps['Dosing'] == 'Flat 900 mg') & exps['AUC24'].notna()].copy()
flat['CL_emp'] = 900.0 / flat['AUC24']

# Sanity-trim implausible CL values (one intermediate subject had near-zero
# AUC24 yielding CL ~ 1e4 L/h, which would distort the bootstrap pool)
flat = flat[flat['CL_emp'].between(0.5, 200)].copy()

# Reference targets (observed intermediate × 900 mg, before trim is identical)
ref = flat[flat['Phenotype'] == 'Intermediate']
AUC_TARGET_MED = ref['AUC24'].median()
AUC_TARGET_LO  = ref['AUC24'].quantile(0.25)
AUC_TARGET_HI  = ref['AUC24'].quantile(0.75)
C24_TARGET_MED = ref['C24'].median()
C24_TARGET_LO  = ref['C24'].quantile(0.25)
C24_TARGET_HI  = ref['C24'].quantile(0.75)

# Empirical (CL, WT) pool per phenotype
cl_pool = {p: flat[flat['Phenotype'] == p][['CL_emp', 'WT']].values
           for p in PHEN_ORDER}

print('Empirical CL pool sizes (subjects per phenotype):')
for p in PHEN_ORDER:
    a = cl_pool[p]
    print(f'  {p:12s}: n={len(a)}  CL median={np.median(a[:,0]):.2f}  '
          f'mean WT={np.mean(a[:,1]):.1f} kg')
print(f'\nReference targets (intermediate × 900 mg, observed):')
print(f'  AUC0-24 : median {AUC_TARGET_MED:.2f}  IQR {AUC_TARGET_LO:.2f}-{AUC_TARGET_HI:.2f}')
print(f'  C24     : median {C24_TARGET_MED:.4f}  IQR {C24_TARGET_LO:.4f}-{C24_TARGET_HI:.4f}')

# ---------- analytical 2-cmt FO absorption --------------------------------
def _hybrid(cl, v2, q, v3):
    k20 = cl / v2; k23 = q / v2; k32 = q / v3
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

def _c24(dose, ka, cl, v2, q, v3):
    k32, a, b = _hybrid(cl, v2, q, v3)
    pre = ka * dose / v2
    ta = (k32 - a)  / ((ka - a) * (b - a))  * np.exp(-24 * a)
    tb = (k32 - b)  / ((ka - b) * (a - b))  * np.exp(-24 * b)
    tk = (k32 - ka) / ((a - ka) * (b - ka)) * np.exp(-24 * ka)
    return pre * (ta + tb + tk)

# ---------- simulation -----------------------------------------------------
N_SUBJ = 2000
DOSES  = np.arange(200, 1801, 50)
rng    = np.random.default_rng(20260608)

records = []
for phen in PHEN_ORDER:
    pool = cl_pool[phen]
    idx  = rng.integers(0, len(pool), size=N_SUBJ)
    cl_i = pool[idx, 0]
    wts  = pool[idx, 1]
    v2_i = tV2 * (wts / REF_WT) ** EXP_V_WT
    v3_i = tV3 * (wts / REF_WT) ** EXP_V_WT
    q_i  = np.full(N_SUBJ, tQ)
    for dose in DOSES:
        auc = _auc24(dose, tKA, cl_i, v2_i, q_i, v3_i)
        c24 = _c24 (dose, tKA, cl_i, v2_i, q_i, v3_i)
        in_auc  = (auc >= AUC_TARGET_LO) & (auc <= AUC_TARGET_HI)
        in_c24  = (c24 >= C24_TARGET_LO) & (c24 <= C24_TARGET_HI)
        records.append({
            'Phenotype': phen, 'Dose_mg': int(dose),
            'AUC_med': float(np.median(auc)),
            'AUC_q1':  float(np.quantile(auc, 0.25)),
            'AUC_q3':  float(np.quantile(auc, 0.75)),
            'C24_med': float(np.median(c24)),
            'C24_q1':  float(np.quantile(c24, 0.25)),
            'C24_q3':  float(np.quantile(c24, 0.75)),
            'pct_in_auc_target':  float(np.mean(in_auc) * 100),
            'pct_in_c24_target':  float(np.mean(in_c24) * 100),
            'pct_in_dual_target': float(np.mean(in_auc & in_c24) * 100),
            'pct_supra':          float(np.mean(c24 > SUPRA) * 100),
        })

sim = pd.DataFrame(records)
sim.to_csv(SIM / 'dose_finding_empirical_results.csv', index=False)

# ---------- recommended doses (3 logics) ----------------------------------
def _pick(phen, col, target):
    s = sim[sim['Phenotype'] == phen].copy()
    s['gap'] = (s[col] - target).abs()
    return s.sort_values(['gap', 'pct_supra']).iloc[0]

def _pick_ceiling(phen):
    s = sim[sim['Phenotype'] == phen].sort_values('Dose_mg')
    ok = s[s['pct_supra'] <= SUPRA_CAP]
    return ok.iloc[-1] if len(ok) else s.iloc[0]

opt_auc = pd.DataFrame([_pick(p, 'AUC_med', AUC_TARGET_MED) for p in PHEN_ORDER]
                       )[['Phenotype', 'Dose_mg', 'AUC_med', 'C24_med',
                          'pct_in_auc_target', 'pct_supra']].rename(
                          columns={'Dose_mg': 'AUC_match_dose_mg'})
opt_c24 = pd.DataFrame([_pick(p, 'C24_med', C24_TARGET_MED) for p in PHEN_ORDER]
                       )[['Phenotype', 'Dose_mg', 'pct_in_c24_target']
                       ].rename(columns={'Dose_mg': 'C24_match_dose_mg'})
opt_ceil = pd.DataFrame([_pick_ceiling(p) for p in PHEN_ORDER]
                        )[['Phenotype', 'Dose_mg', 'pct_supra']].rename(
                          columns={'Dose_mg': 'Ceiling_dose_mg',
                                   'pct_supra': 'pct_supra_at_ceiling'})

opt = opt_auc.merge(opt_c24, on='Phenotype').merge(opt_ceil, on='Phenotype')
opt.to_csv(SIM / 'dose_finding_empirical_recommended.csv', index=False)
print('\nRecommended doses (empirical anchor):')
print(opt.to_string(index=False))

# ---------- 6-panel dose-response figure ----------------------------------
fig, axes = plt.subplots(2, 3, figsize=(16, 9.4))
(axA, axB, axC), (axD, axE, axF) = axes
for phen in PHEN_ORDER:
    sub = sim[sim['Phenotype'] == phen]; c = PHEN_COLOUR[phen]
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

# A — AUC full
axA.axhspan(AUC_TARGET_LO, AUC_TARGET_HI, color='grey', alpha=0.25,
            label=f'AUC target IQR ({AUC_TARGET_LO:.0f}-{AUC_TARGET_HI:.0f})')
axA.axhline(AUC_TARGET_MED, color='black', linestyle='--', linewidth=1,
            alpha=0.6, label=f'AUC target median ({AUC_TARGET_MED:.0f})')
axA.set_xlabel('Dose (mg)'); axA.set_ylabel('Simulated AUC$_{0-24}$ (mg·h/L)')
axA.legend(loc='upper left', fontsize=8); axA.grid(axis='y', alpha=.3)
axA.text(-0.10, 1.02, 'A', transform=axA.transAxes, fontsize=14,
         fontweight='bold', va='bottom', ha='left')

# B — C24 full log
axB.axhspan(C24_TARGET_LO, C24_TARGET_HI, color='grey', alpha=0.25,
            label=f'C24 target IQR ({C24_TARGET_LO:.3f}-{C24_TARGET_HI:.3f})')
axB.axhline(C24_TARGET_MED, color='black', linestyle='--', linewidth=1,
            alpha=0.6, label=f'C24 target median ({C24_TARGET_MED:.3f})')
axB.axhline(SUPRA, color='red', linestyle='--', linewidth=1, alpha=0.6,
            label=f'Supratherapeutic ({SUPRA} µg/mL)')
axB.set_xlabel('Dose (mg)'); axB.set_ylabel('Simulated C$_{24}$ (µg/mL)')
axB.set_yscale('log'); axB.legend(loc='upper left', fontsize=8)
axB.grid(axis='y', alpha=.3)
axB.text(-0.10, 1.02, 'B', transform=axB.transAxes, fontsize=14,
         fontweight='bold', va='bottom', ha='left')

# C — supra full
axC.axhline(SUPRA_CAP, color='red', linestyle='--', linewidth=1, alpha=0.6,
            label=f'{SUPRA_CAP:.0f}% supratherapeutic cap')
axC.set_xlabel('Dose (mg)')
axC.set_ylabel(f'% subjects with C$_{{24}}$ > {SUPRA} µg/mL')
axC.legend(loc='upper left', fontsize=8); axC.grid(axis='y', alpha=.3)
axC.text(-0.10, 1.02, 'C', transform=axC.transAxes, fontsize=14,
         fontweight='bold', va='bottom', ha='left')

# D — AUC zoom
auc_pad = 0.2 * (AUC_TARGET_HI - AUC_TARGET_LO)
axD.axhspan(AUC_TARGET_LO, AUC_TARGET_HI, color='grey', alpha=0.25,
            label=f'AUC target IQR ({AUC_TARGET_LO:.0f}-{AUC_TARGET_HI:.0f})')
axD.axhline(AUC_TARGET_MED, color='black', linestyle='--', linewidth=1,
            alpha=0.6, label=f'AUC target median ({AUC_TARGET_MED:.0f})')
axD.set_xlabel('Dose (mg)'); axD.set_ylabel('Simulated AUC$_{0-24}$ (mg·h/L)')
axD.set_ylim(AUC_TARGET_LO - auc_pad, AUC_TARGET_HI + auc_pad)
axD.legend(loc='upper left', fontsize=8); axD.grid(axis='y', alpha=.3)
axD.text(-0.10, 1.02, 'D', transform=axD.transAxes, fontsize=14,
         fontweight='bold', va='bottom', ha='left')
axD.set_title('Zoom: AUC target IQR', fontsize=10)

# E — C24 zoom linear
c24_pad_lo = 0.6 * C24_TARGET_LO
c24_top    = max(C24_TARGET_HI * 2, SUPRA * 1.05)
axE.axhspan(C24_TARGET_LO, C24_TARGET_HI, color='grey', alpha=0.25,
            label=f'C24 target IQR ({C24_TARGET_LO:.3f}-{C24_TARGET_HI:.3f})')
axE.axhline(C24_TARGET_MED, color='black', linestyle='--', linewidth=1,
            alpha=0.6, label=f'C24 target median ({C24_TARGET_MED:.3f})')
axE.axhline(SUPRA, color='red', linestyle='--', linewidth=1, alpha=0.6,
            label=f'Supratherapeutic ({SUPRA} µg/mL)')
axE.set_xlabel('Dose (mg)'); axE.set_ylabel('Simulated C$_{24}$ (µg/mL)')
axE.set_ylim(c24_pad_lo, c24_top)
axE.legend(loc='upper left', fontsize=8); axE.grid(axis='y', alpha=.3)
axE.text(-0.10, 1.02, 'E', transform=axE.transAxes, fontsize=14,
         fontweight='bold', va='bottom', ha='left')
axE.set_title('Zoom: C$_{24}$ target IQR', fontsize=10)

# F — supra zoom
axF.axhline(SUPRA_CAP, color='red', linestyle='--', linewidth=1, alpha=0.6,
            label=f'{SUPRA_CAP:.0f}% supratherapeutic cap')
axF.set_xlabel('Dose (mg)')
axF.set_ylabel(f'% subjects with C$_{{24}}$ > {SUPRA} µg/mL')
axF.set_ylim(0, 20)
axF.legend(loc='upper left', fontsize=8); axF.grid(axis='y', alpha=.3)
axF.text(-0.10, 1.02, 'F', transform=axF.transAxes, fontsize=14,
         fontweight='bold', va='bottom', ha='left')
axF.set_title('Zoom: 0-20% supratherapeutic', fontsize=10)

fig.suptitle('Empirical-anchored dose-finding (CL resampled from observed '
             'flat 900 mg arm; reference targets reproduced by construction)',
             y=1.01, fontsize=10.5)
fig.tight_layout()
fig.savefig(FIG / 'Figure_dose_finding_empirical.png', dpi=300,
            bbox_inches='tight')
fig.savefig(FIG / 'Figure_dose_finding_empirical.pdf', bbox_inches='tight')
plt.close(fig)

# ---------- target attainment figure (3 panels) --------------------------
PROTOCOL = {'Slow': 300, 'Intermediate': 900, 'Rapid': 1500}
fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), sharey=True)
specs = [
    (axes[0], 'pct_in_auc_target',
     f'AUC IQR ({AUC_TARGET_LO:.0f}-{AUC_TARGET_HI:.0f})',
     'A', 'In AUC$_{0-24}$ target IQR'),
    (axes[1], 'pct_in_c24_target',
     f'C24 IQR ({C24_TARGET_LO:.3f}-{C24_TARGET_HI:.3f})',
     'B', 'In C$_{24}$ target IQR'),
    (axes[2], 'pct_in_dual_target',
     'Both IQRs simultaneously', 'C', 'In both target IQRs (dual)'),
]
for ax, col, sub_title, lbl, ylab in specs:
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
    ax.set_ylim(0, 100); ax.grid(axis='y', alpha=.3)
    ax.set_title(sub_title, fontsize=10)
    ax.legend(loc='upper right', fontsize=8)
    ax.text(-0.10, 1.02, lbl, transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='bottom', ha='left')
fig.suptitle('Empirical-anchored target attainment vs dose '
             '(dotted lines = protocol doses 300 / 900 / 1500 mg)',
             y=1.02, fontsize=10.5)
fig.tight_layout()
fig.savefig(FIG / 'Figure_dose_finding_empirical_attainment.png',
            dpi=300, bbox_inches='tight')
fig.savefig(FIG / 'Figure_dose_finding_empirical_attainment.pdf',
            bbox_inches='tight')
plt.close(fig)

print('\nWrote dose_finding_empirical_results.csv,')
print('      dose_finding_empirical_recommended.csv,')
print('      Figure_dose_finding_empirical.{png,pdf},')
print('      Figure_dose_finding_empirical_attainment.{png,pdf}')
