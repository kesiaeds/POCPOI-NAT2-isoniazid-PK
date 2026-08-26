"""Parameter-uncertainty bootstrap for the dose-finding simulation.

Extends dose_finding_simulation.py by sampling popPK fixed effects from
their estimated sampling distribution (parametric bootstrap from N(theta, SE^2)
with positivity truncation for physical parameters), rerunning the full Monte
Carlo per draw, and reporting recommended doses as point estimate + 90% CI.

Outputs (analysis/output/simulation/):
  dose_finding_bootstrap_draws.csv      one row per (phenotype, logic, draw)
  dose_finding_bootstrap_summary.csv    point estimate + 90% CI per logic
analysis/output/figures/:
  Figure_dose_finding_bootstrap.{png,pdf}
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

# ---------- popPK point estimates and SEs ---------------------------------
pars = pd.read_csv(ROOT / 'sherlock' / 'output' / 'popPK_parameters.csv',
                   index_col=0)

THETA = {
    'tKA':  (pars.loc['tKA', 'Estimate'],  pars.loc['tKA', 'SE']),
    'tCL':  (pars.loc['tCL', 'Estimate'],  pars.loc['tCL', 'SE']),
    'tV2':  (pars.loc['tV2', 'Estimate'],  pars.loc['tV2', 'SE']),
    'tQ':   (pars.loc['tQ',  'Estimate'],  pars.loc['tQ',  'SE']),
    'tV3':  (pars.loc['tV3', 'Estimate'],  pars.loc['tV3', 'SE']),
    'dCL_Slow':  (pars.loc['dCldACEPROFILE1', 'Estimate'],
                  pars.loc['dCldACEPROFILE1', 'SE']),
    'dCL_Rapid': (pars.loc['dCldACEPROFILE3', 'Estimate'],
                  pars.loc['dCldACEPROFILE3', 'SE']),
}
BSV_CV = pars.loc['eta.cl', 'BSV(CV%)'] / 100.0
OMEGA  = np.sqrt(np.log1p(BSV_CV ** 2))   # treated as fixed (no SE reported)

EXP_CL_WT = 0.75
EXP_V_WT  = 1.0
REF_WT    = 70.0
SUPRA     = 0.15
SUPRA_CAP = 15.0

PHEN_ORDER = ['Slow', 'Intermediate', 'Rapid']
PHEN_ID    = {'Slow': 1, 'Intermediate': 2, 'Rapid': 3}
PHEN_COLOUR = {'Slow': '#0072B2', 'Intermediate': '#009E73', 'Rapid': '#D55E00'}

# ---------- reference targets (observed intermediate * flat 900 mg) -------
exps = pd.read_csv(ART / 'exposures_subject_occasion.csv')
ref  = exps[(exps['Phenotype'] == 'Intermediate') &
            (exps['Dosing']    == 'Flat 900 mg')]
AUC_TARGET_MED  = ref['AUC24'].dropna().median()
C24_TARGET_MED  = ref['C24'].dropna().median()

# ---------- weight pool ---------------------------------------------------
nm = pd.read_csv(DATA / 'pk_nonmem_combined.csv')
wt_pool = nm.drop_duplicates('ID')[['ACE_PROFILE', 'WT']]

# ---------- analytical 2-cmt FO absorption --------------------------------
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

def _c24(dose, ka, cl, v2, q, v3):
    k32, a, b = _hybrid(cl, v2, q, v3)
    pre = ka * dose / v2
    ta = (k32 - a)  / ((ka - a) * (b - a))  * np.exp(-24 * a)
    tb = (k32 - b)  / ((ka - b) * (a - b))  * np.exp(-24 * b)
    tk = (k32 - ka) / ((a - ka) * (b - ka)) * np.exp(-24 * ka)
    return pre * (ta + tb + tk)

# ---------- bootstrap sampler ---------------------------------------------
def sample_theta(rng):
    """One draw of fixed effects. Truncate physical parameters at small +eps."""
    eps = 1e-3
    def pos(name):
        m, se = THETA[name]
        x = rng.normal(m, se)
        return max(x, eps)
    return {
        'tKA':       pos('tKA'),
        'tCL':       pos('tCL'),
        'tV2':       pos('tV2'),
        'tQ':        pos('tQ'),
        'tV3':       pos('tV3'),
        'dCL_Slow':  rng.normal(*THETA['dCL_Slow']),    # log-scale shift, unbounded
        'dCL_Rapid': rng.normal(*THETA['dCL_Rapid']),
    }

# ---------- single Monte Carlo simulation given one parameter draw --------
N_SUBJ = 2000
DOSES  = np.arange(200, 1801, 50)

def simulate_one_draw(th, rng):
    """Return dict[phen] -> DataFrame(Dose_mg, AUC_med, C24_med, pct_supra)."""
    dCL = {'Slow': th['dCL_Slow'], 'Intermediate': 0.0, 'Rapid': th['dCL_Rapid']}
    out = {}
    for phen in PHEN_ORDER:
        pool = wt_pool[wt_pool['ACE_PROFILE'] == PHEN_ID[phen]]['WT'].values
        wts  = rng.choice(pool, size=N_SUBJ, replace=True)
        eta  = rng.normal(0.0, OMEGA, size=N_SUBJ)
        cl_i = th['tCL'] * (wts / REF_WT) ** EXP_CL_WT * np.exp(dCL[phen] + eta)
        v2_i = th['tV2'] * (wts / REF_WT) ** EXP_V_WT
        v3_i = th['tV3'] * (wts / REF_WT) ** EXP_V_WT
        q_i  = np.full(N_SUBJ, th['tQ'])
        rows = []
        for dose in DOSES:
            auc = _auc24(dose, th['tKA'], cl_i, v2_i, q_i, v3_i)
            c24 = _c24 (dose, th['tKA'], cl_i, v2_i, q_i, v3_i)
            rows.append({
                'Dose_mg':   int(dose),
                'AUC_med':   float(np.median(auc)),
                'C24_med':   float(np.median(c24)),
                'pct_supra': float(np.mean(c24 > SUPRA) * 100),
            })
        out[phen] = pd.DataFrame(rows)
    return out

# ---------- recommendation extractors -------------------------------------
def pick_auc(df, target):
    return int(df.iloc[(df['AUC_med'] - target).abs().argmin()]['Dose_mg'])

def pick_c24(df, target):
    return int(df.iloc[(df['C24_med'] - target).abs().argmin()]['Dose_mg'])

def pick_ceiling(df, cap):
    ok = df[df['pct_supra'] <= cap].sort_values('Dose_mg')
    return int(ok.iloc[-1]['Dose_mg']) if len(ok) else int(df['Dose_mg'].min())

# ---------- bootstrap loop ------------------------------------------------
N_BOOT = 500
rng    = np.random.default_rng(20260605)

records = []
for b in range(N_BOOT):
    th  = sample_theta(rng)
    sim = simulate_one_draw(th, rng)
    for phen in PHEN_ORDER:
        df = sim[phen]
        records.append({'draw': b, 'Phenotype': phen,
                        'logic': 'AUC_match',
                        'dose':  pick_auc(df, AUC_TARGET_MED)})
        records.append({'draw': b, 'Phenotype': phen,
                        'logic': 'C24_match',
                        'dose':  pick_c24(df, C24_TARGET_MED)})
        records.append({'draw': b, 'Phenotype': phen,
                        'logic': 'C24_ceiling',
                        'dose':  pick_ceiling(df, SUPRA_CAP)})
    if (b + 1) % 50 == 0:
        print(f'  bootstrap draw {b+1}/{N_BOOT}')

draws = pd.DataFrame(records)
draws.to_csv(SIM / 'dose_finding_bootstrap_draws.csv', index=False)

# ---------- summarize ------------------------------------------------------
def _summ(g):
    return pd.Series({
        'median': int(np.median(g)),
        'q05':    int(np.quantile(g, 0.05)),
        'q95':    int(np.quantile(g, 0.95)),
        'mean':   float(np.mean(g)),
    })

summary = (draws.groupby(['Phenotype', 'logic'])['dose']
                .apply(_summ).unstack(-1).reset_index())
summary = summary[['Phenotype', 'logic', 'median', 'q05', 'q95', 'mean']]
summary['Phenotype'] = pd.Categorical(summary['Phenotype'], PHEN_ORDER, ordered=True)
summary = summary.sort_values(['Phenotype', 'logic']).reset_index(drop=True)
summary.to_csv(SIM / 'dose_finding_bootstrap_summary.csv', index=False)

print('\nReference targets (Intermediate × Flat 900 mg)')
print(f'  AUC0-24 : median {AUC_TARGET_MED:.1f} mg·h/L')
print(f'  C24     : median {C24_TARGET_MED:.4f} µg/mL')
print('\nBootstrap recommended doses (median [5th-95th percentile], mg):')
for phen in PHEN_ORDER:
    print(f'\n  {phen}:')
    for logic in ['AUC_match', 'C24_match', 'C24_ceiling']:
        row = summary[(summary['Phenotype'] == phen) &
                      (summary['logic']     == logic)].iloc[0]
        print(f'    {logic:12s}: {int(row["median"]):4d} mg '
              f'[{int(row["q05"]):4d} - {int(row["q95"]):4d}]')

# ---------- figure: recommended-dose distributions ------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), sharey=True)
logics = ['AUC_match', 'C24_match', 'C24_ceiling']
logic_label = {'AUC_match': 'AUC-matched',
               'C24_match': 'C$_{24}$-matched',
               'C24_ceiling': 'C$_{24}$-ceiling\n(≤15% supra)'}
PROTOCOL = {'Slow': 300, 'Intermediate': 900, 'Rapid': 1500}

for ax, logic in zip(axes, logics):
    parts_data = []
    for phen in PHEN_ORDER:
        parts_data.append(
            draws[(draws['Phenotype'] == phen) &
                  (draws['logic']     == logic)]['dose'].values)
    parts = ax.violinplot(parts_data, positions=[1, 2, 3], showmedians=True,
                          widths=0.7)
    for body, phen in zip(parts['bodies'], PHEN_ORDER):
        body.set_facecolor(PHEN_COLOUR[phen]); body.set_alpha(0.55)
        body.set_edgecolor(PHEN_COLOUR[phen])
    for key in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
        parts[key].set_color('black'); parts[key].set_linewidth(1)
    for i, phen in enumerate(PHEN_ORDER, start=1):
        ax.scatter([i], [PROTOCOL[phen]], marker='D', color='black',
                   s=42, zorder=4, label='Protocol' if i == 1 else None)
    ax.set_xticks([1, 2, 3]); ax.set_xticklabels(PHEN_ORDER)
    ax.set_title(logic_label[logic])
    ax.set_ylim(0, 1900)
    ax.grid(axis='y', alpha=.3)
axes[0].set_ylabel('Recommended dose (mg)')
axes[0].legend(loc='upper left', fontsize=9)
fig.suptitle(f'Bootstrap dose-recommendation distributions ({N_BOOT} draws of '
             'popPK fixed effects)', y=1.02, fontsize=11)
fig.tight_layout()
fig.savefig(FIG / 'Figure_dose_finding_bootstrap.png', dpi=300,
            bbox_inches='tight')
fig.savefig(FIG / 'Figure_dose_finding_bootstrap.pdf', bbox_inches='tight')
plt.close(fig)
print('\nWrote', SIM, 'and Figure_dose_finding_bootstrap.{png,pdf}')
