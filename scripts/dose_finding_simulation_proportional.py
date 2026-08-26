"""Dose-proportional empirical dose-finding (Option 1, purely empirical).

For each subject in the flat 900 mg arm, both AUC0-24 and C24 are scaled
linearly with dose:
    AUC(D) = (D / 900) * AUC_obs
    C24(D) = (D / 900) * C24_obs
No popPK structural model is used. The only assumption is linear PK over
the dose range, which is well established for isoniazid.

By construction this reproduces both AUC AND C24 exactly at 900 mg for every
phenotype. Recommendations therefore inherit the empirical between-subject
variability of the trial cohort without any model bias.

Outputs (analysis/output/simulation/):
  dose_finding_proportional_results.csv
  dose_finding_proportional_recommended.csv
analysis/output/figures/:
  Figure_dose_finding_proportional.{png,pdf}              (6-panel dose-response)
  Figure_dose_finding_proportional_attainment.{png,pdf}   (3-panel attainment)
  Figure_dose_finding_methods_comparison.{png,pdf}        (3x3 method comparison)
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path('/Users/kesiaeds/repos/POCPOI')
ART  = ROOT / 'analysis' / 'output'
FIG  = ART / 'figures'
SIM  = ART / 'simulation'

SUPRA     = 0.15
SUPRA_CAP = 15.0
PHEN_ORDER  = ['Slow', 'Intermediate', 'Rapid']
PHEN_COLOUR = {'Slow': '#0072B2', 'Intermediate': '#009E73', 'Rapid': '#D55E00'}
METHOD_COLOUR = {
    'popPK-typical':       '#888888',
    'Empirical CL anchor': '#D55E00',
    'Dose-proportional':   '#0072B2',
}
PROTOCOL = {'Slow': 300, 'Intermediate': 900, 'Rapid': 1500}

# ---------- observed reference arm ----------------------------------------
exps = pd.read_csv(ART / 'exposures_subject_occasion.csv')
flat = exps[(exps['Dosing'] == 'Flat 900 mg') &
            exps['AUC24'].notna() & exps['C24'].notna()].copy()

# Observed intermediate × 900 mg defines the targets
ref = flat[flat['Phenotype'] == 'Intermediate']
AUC_TARGET_MED = ref['AUC24'].median()
AUC_TARGET_LO  = ref['AUC24'].quantile(0.25)
AUC_TARGET_HI  = ref['AUC24'].quantile(0.75)
C24_TARGET_MED = ref['C24'].median()
C24_TARGET_LO  = ref['C24'].quantile(0.25)
C24_TARGET_HI  = ref['C24'].quantile(0.75)

# Per-phenotype (AUC_per_mg, C24_per_mg) pool
pool = {p: np.column_stack([
            flat.loc[flat['Phenotype'] == p, 'AUC24'].values / 900.0,
            flat.loc[flat['Phenotype'] == p, 'C24'].values   / 900.0,
        ]) for p in PHEN_ORDER}

print('Dose-proportional empirical pool sizes (subjects per phenotype):')
for p in PHEN_ORDER:
    a = pool[p]
    print(f'  {p:12s}: n={len(a)}  AUC/mg median={np.median(a[:,0]):.5f}  '
          f'C24/mg median={np.median(a[:,1]):.6f}')
print(f'\nTargets (intermediate × 900 mg, observed):')
print(f'  AUC0-24 : median {AUC_TARGET_MED:.2f}  '
      f'IQR {AUC_TARGET_LO:.2f}-{AUC_TARGET_HI:.2f}')
print(f'  C24     : median {C24_TARGET_MED:.4f}  '
      f'IQR {C24_TARGET_LO:.4f}-{C24_TARGET_HI:.4f}')

# ---------- simulation -----------------------------------------------------
N_SUBJ = 2000
DOSES  = np.arange(200, 1801, 50)
rng    = np.random.default_rng(20260608)

records = []
for phen in PHEN_ORDER:
    p     = pool[phen]
    idx   = rng.integers(0, len(p), size=N_SUBJ)
    a_pm  = p[idx, 0]    # AUC per mg
    c_pm  = p[idx, 1]    # C24 per mg
    for dose in DOSES:
        auc = dose * a_pm
        c24 = dose * c_pm
        in_auc = (auc >= AUC_TARGET_LO) & (auc <= AUC_TARGET_HI)
        in_c24 = (c24 >= C24_TARGET_LO) & (c24 <= C24_TARGET_HI)
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
sim.to_csv(SIM / 'dose_finding_proportional_results.csv', index=False)

# ---------- recommendations (3 logics) ------------------------------------
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
opt.to_csv(SIM / 'dose_finding_proportional_recommended.csv', index=False)
print('\nDose-proportional recommended doses:')
print(opt.to_string(index=False))

# ---------- 6-panel dose-response figure ----------------------------------
def _make_6panel(sim_df, out_name, suptitle):
    fig, axes = plt.subplots(2, 3, figsize=(16, 9.4))
    (axA, axB, axC), (axD, axE, axF) = axes
    for phen in PHEN_ORDER:
        sub = sim_df[sim_df['Phenotype'] == phen]; c = PHEN_COLOUR[phen]
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
    axA.axhspan(AUC_TARGET_LO, AUC_TARGET_HI, color='grey', alpha=0.25,
                label=f'AUC target IQR ({AUC_TARGET_LO:.0f}-{AUC_TARGET_HI:.0f})')
    axA.axhline(AUC_TARGET_MED, color='black', linestyle='--', linewidth=1,
                alpha=0.6, label=f'AUC target median ({AUC_TARGET_MED:.0f})')
    axA.set_xlabel('Dose (mg)'); axA.set_ylabel('Simulated AUC$_{0-24}$ (mg·h/L)')
    axA.legend(loc='upper left', fontsize=8); axA.grid(axis='y', alpha=.3)
    axA.text(-0.10, 1.02, 'A', transform=axA.transAxes, fontsize=14,
             fontweight='bold', va='bottom', ha='left')
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
    axC.axhline(SUPRA_CAP, color='red', linestyle='--', linewidth=1, alpha=0.6,
                label=f'{SUPRA_CAP:.0f}% supratherapeutic cap')
    axC.set_xlabel('Dose (mg)')
    axC.set_ylabel(f'% subjects with C$_{{24}}$ > {SUPRA} µg/mL')
    axC.legend(loc='upper left', fontsize=8); axC.grid(axis='y', alpha=.3)
    axC.text(-0.10, 1.02, 'C', transform=axC.transAxes, fontsize=14,
             fontweight='bold', va='bottom', ha='left')
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
    axF.axhline(SUPRA_CAP, color='red', linestyle='--', linewidth=1, alpha=0.6,
                label=f'{SUPRA_CAP:.0f}% supratherapeutic cap')
    axF.set_xlabel('Dose (mg)')
    axF.set_ylabel(f'% subjects with C$_{{24}}$ > {SUPRA} µg/mL')
    axF.set_ylim(0, 20)
    axF.legend(loc='upper left', fontsize=8); axF.grid(axis='y', alpha=.3)
    axF.text(-0.10, 1.02, 'F', transform=axF.transAxes, fontsize=14,
             fontweight='bold', va='bottom', ha='left')
    axF.set_title('Zoom: 0-20% supratherapeutic', fontsize=10)
    fig.suptitle(suptitle, y=1.01, fontsize=10.5)
    fig.tight_layout()
    fig.savefig(FIG / f'{out_name}.png', dpi=300, bbox_inches='tight')
    fig.savefig(FIG / f'{out_name}.pdf', bbox_inches='tight')
    plt.close(fig)

_make_6panel(sim, 'Figure_dose_finding_proportional',
             'Dose-proportional empirical dose-finding '
             '(linear scaling of observed AUC and C24 from flat 900 mg arm)')

# ---------- target attainment figure --------------------------------------
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
    ax.set_xlabel('Dose (mg)'); ax.set_ylabel(f'% subjects {ylab}')
    ax.set_ylim(0, 100); ax.grid(axis='y', alpha=.3)
    ax.set_title(sub_title, fontsize=10)
    ax.legend(loc='upper right', fontsize=8)
    ax.text(-0.10, 1.02, lbl, transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='bottom', ha='left')
fig.suptitle('Dose-proportional target attainment vs dose',
             y=1.02, fontsize=10.5)
fig.tight_layout()
fig.savefig(FIG / 'Figure_dose_finding_proportional_attainment.png',
            dpi=300, bbox_inches='tight')
fig.savefig(FIG / 'Figure_dose_finding_proportional_attainment.pdf',
            bbox_inches='tight')
plt.close(fig)

# ---------- 3x3 method-comparison figure -----------------------------------
sim_pk  = pd.read_csv(SIM / 'dose_finding_results.csv')               # popPK
sim_emp = pd.read_csv(SIM / 'dose_finding_empirical_results.csv')     # empirical CL
sim_pp  = sim                                                         # dose-proportional
methods = [
    ('popPK-typical',       sim_pk),
    ('Empirical CL anchor', sim_emp),
    ('Dose-proportional',   sim_pp),
]
metrics = [
    ('AUC_med',   'AUC$_{0-24}$ (mg·h/L)',              False, AUC_TARGET_MED),
    ('C24_med',   'C$_{24}$ (µg/mL)',                    True,  C24_TARGET_MED),
    ('pct_supra', f'% C$_{{24}}$ > {SUPRA} µg/mL',       False, None),
]
fig, axes = plt.subplots(3, 3, figsize=(15, 11), sharex=True)
for r, phen in enumerate(PHEN_ORDER):
    for c, (col, ylab, logy, target) in enumerate(metrics):
        ax = axes[r, c]
        for m_name, df in methods:
            sub = df[df['Phenotype'] == phen]
            ax.plot(sub['Dose_mg'], sub[col], '-o',
                    color=METHOD_COLOUR[m_name],
                    label=m_name if r == 0 and c == 0 else None,
                    linewidth=1.6, markersize=3, alpha=0.9)
        if target is not None:
            ax.axhline(target, color='black', linestyle='--', linewidth=1,
                       alpha=0.6,
                       label='Observed target median'
                             if r == 0 and c == 0 else None)
        if col == 'C24_med':
            ax.axhline(SUPRA, color='red', linestyle=':', linewidth=1,
                       alpha=0.5,
                       label='0.15 µg/mL'
                             if r == 0 and c == 0 else None)
        if col == 'pct_supra':
            ax.axhline(SUPRA_CAP, color='red', linestyle=':', linewidth=1,
                       alpha=0.5,
                       label='15% cap'
                             if r == 0 and c == 0 else None)
        ax.axvline(PROTOCOL[phen], color=PHEN_COLOUR[phen], linestyle=':',
                   linewidth=1.2, alpha=0.7,
                   label=f'Protocol {PROTOCOL[phen]} mg'
                         if r == 0 and c == 0 else None)
        if logy:
            ax.set_yscale('log')
        ax.grid(axis='y', alpha=.3)
        if c == 0:
            ax.set_ylabel(f'{phen}\n{ylab}', fontsize=10)
        if r == 2:
            ax.set_xlabel('Dose (mg)')
        if r == 0:
            ax.set_title(ylab, fontsize=11)
# legend in upper-left panel
axes[0, 0].legend(loc='upper left', fontsize=7.5, ncol=1, framealpha=0.9)
fig.suptitle('Method comparison: popPK-typical vs empirical CL anchor vs '
             'dose-proportional\n(rows = phenotype, columns = exposure metric; '
             'dotted vertical = protocol dose)',
             y=1.01, fontsize=11)
fig.tight_layout()
fig.savefig(FIG / 'Figure_dose_finding_methods_comparison.png',
            dpi=300, bbox_inches='tight')
fig.savefig(FIG / 'Figure_dose_finding_methods_comparison.pdf',
            bbox_inches='tight')
plt.close(fig)

# ---------- side-by-side recommended-dose comparison table -----------------
rec_pk  = pd.read_csv(SIM / 'dose_finding_recommended.csv')
rec_emp = pd.read_csv(SIM / 'dose_finding_empirical_recommended.csv')
rec_pp  = opt

comp = pd.DataFrame({
    'Phenotype': PHEN_ORDER,
    'Protocol': [PROTOCOL[p] for p in PHEN_ORDER],
    'popPK_AUC_match': [int(rec_pk[rec_pk.Phenotype==p].AUC_match_dose_mg.iloc[0])
                       for p in PHEN_ORDER],
    'Emp_AUC_match':   [int(rec_emp[rec_emp.Phenotype==p].AUC_match_dose_mg.iloc[0])
                       for p in PHEN_ORDER],
    'DP_AUC_match':    [int(rec_pp[rec_pp.Phenotype==p].AUC_match_dose_mg.iloc[0])
                       for p in PHEN_ORDER],
    'popPK_Ceiling':   [int(rec_pk[rec_pk.Phenotype==p].Ceiling_dose_mg.iloc[0])
                       for p in PHEN_ORDER],
    'Emp_Ceiling':     [int(rec_emp[rec_emp.Phenotype==p].Ceiling_dose_mg.iloc[0])
                       for p in PHEN_ORDER],
    'DP_Ceiling':      [int(rec_pp[rec_pp.Phenotype==p].Ceiling_dose_mg.iloc[0])
                       for p in PHEN_ORDER],
})
comp.to_csv(SIM / 'dose_finding_method_comparison.csv', index=False)
print('\n========== Method comparison ==========')
print(comp.to_string(index=False))
print('\nWrote dose_finding_proportional_*, Figure_dose_finding_proportional*,')
print('      Figure_dose_finding_methods_comparison.{png,pdf}, '
      'dose_finding_method_comparison.csv')
