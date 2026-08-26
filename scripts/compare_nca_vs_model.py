"""Side-by-side comparison: NCA on observed DV vs MM-model-derived (IPRED)
exposure metrics by phenotype x dosing strategy.

Writes a single CSV that puts both methods next to each other so the user
can decide which to use in the manuscript.
"""
import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats

PHEN = {1: 'Slow', 2: 'Intermediate', 3: 'Rapid'}
OCCN = {1: 'Genotype-guided', 2: 'Flat 900 mg'}
SUPRA = 0.15

OUT = Path('/tmp/mm_exposures')

# ---- (a) NCA on observed DV ------------------------------------------------
raw = pd.read_csv('/Users/kesiaeds/Documents/POCPOI_dose_finding/data/pk_nonmem_combined.csv')
obs = raw[raw['EVID'] == 0].copy()

def per_occ(g):
    g = g.sort_values('TAD')
    auc = float(np.trapezoid(g['DV'].values, g['TAD'].values))
    cmax = float(g['DV'].max())
    c24 = float(g.loc[g['TAD'] == 24, 'DV'].iloc[0]) if (g['TAD'] == 24).any() else np.nan
    return pd.Series({'AUC24': auc, 'Cmax': cmax, 'C24': c24})

nca = (obs.groupby(['ID', 'OCC', 'ACE_PROFILE'])
          .apply(per_occ, include_groups=False).reset_index())
nca['Phenotype'] = nca['ACE_PROFILE'].map(PHEN)
nca['Dosing']    = nca['OCC'].map(OCCN)
nca['Method']    = 'NCA (observed DV)'

# ---- (b) MM model-derived (IPRED dense) -----------------------------------
mm = pd.read_csv(OUT / 'exposures_mm_model.csv')
mm['Method'] = 'MM model (IPRED dense)'
mm = mm.rename(columns={'dosenum': 'OCC'})

# ---- combine and summarize ------------------------------------------------
common_cols = ['Phenotype', 'Dosing', 'Method', 'AUC24', 'Cmax', 'C24']
nca_s = nca[common_cols].copy()
mm_s  = mm[common_cols].copy()
both  = pd.concat([nca_s, mm_s], ignore_index=True)

def fmt(x):
    return f"{np.median(x):.2f} ({np.percentile(x, 25):.2f}-{np.percentile(x, 75):.2f})"

rows = []
for (phen, dosing, method), g in both.groupby(['Phenotype', 'Dosing', 'Method'], sort=False):
    rows.append({
        'Phenotype': phen, 'Dosing': dosing, 'Method': method, 'N': len(g),
        'AUC24 median (IQR)': fmt(g['AUC24'].dropna()),
        'Cmax median (IQR)':  fmt(g['Cmax'].dropna()),
        'C24 median (IQR)':   fmt(g['C24'].dropna()),
        'C24 > 0.15 (%)':     f"{100*np.mean(g['C24'] > SUPRA):.1f}%",
    })
summary = pd.DataFrame(rows)

# Order rows: phenotype-dosing first, method second so NCA + MM stack visually
phen_order  = ['Slow', 'Intermediate', 'Rapid']
dose_order  = ['Genotype-guided', 'Flat 900 mg']
meth_order  = ['NCA (observed DV)', 'MM model (IPRED dense)']
summary['Phenotype'] = pd.Categorical(summary['Phenotype'], phen_order)
summary['Dosing']    = pd.Categorical(summary['Dosing'],    dose_order)
summary['Method']    = pd.Categorical(summary['Method'],    meth_order)
summary = summary.sort_values(['Phenotype', 'Dosing', 'Method']).reset_index(drop=True)
summary.to_csv(OUT / 'comparison_NCA_vs_MM_summary.csv', index=False)
print('\n=== Phenotype x dosing summary: NCA vs MM-IPRED ===')
print(summary.to_string(index=False))

# ---- Paired stats: by method, within phenotype (Day 7 vs Day 14) ----------
def paired_long(df, method_name):
    w = df.pivot(index=['ID', 'Phenotype'], columns='Dosing', values=['AUC24', 'Cmax', 'C24']).reset_index()
    w.columns = ['_'.join([c for c in tup if c]) for tup in w.columns]
    out = []
    for phen in phen_order:
        sub = w[w['Phenotype'] == phen].dropna()
        if sub.empty:
            continue
        for metric in ['AUC24', 'Cmax', 'C24']:
            g = sub[f'{metric}_Genotype-guided'].values
            f = sub[f'{metric}_Flat 900 mg'].values
            n = len(g)
            try:
                # Drop pairs with identical values (Wilcoxon ranks zero diffs)
                if np.all(g == f):
                    p = np.nan
                else:
                    p = stats.wilcoxon(g, f).pvalue
            except Exception:
                p = np.nan
            out.append({
                'Method': method_name, 'Metric': metric, 'Phenotype': phen,
                'N_pairs': n,
                'Guided median (IQR)': fmt(g),
                'Flat median (IQR)':   fmt(f),
                'p (Wilcoxon)':        f"{p:.3g}" if not np.isnan(p) else "all pairs identical",
            })
    return pd.DataFrame(out)

p_nca = paired_long(nca[['ID', 'Phenotype', 'Dosing', 'AUC24', 'Cmax', 'C24']],
                    'NCA (observed DV)')
p_mm  = paired_long(mm[['ID', 'Phenotype', 'Dosing', 'AUC24', 'Cmax', 'C24']]
                    .rename(columns={'ID': 'ID'}),
                    'MM model (IPRED dense)')
paired = pd.concat([p_nca, p_mm], ignore_index=True)
paired.to_csv(OUT / 'comparison_NCA_vs_MM_paired.csv', index=False)
print('\n=== Within-phenotype paired stats: NCA vs MM-IPRED ===')
print(paired.to_string(index=False))
print('\nWrote:')
print(' ', OUT / 'comparison_NCA_vs_MM_summary.csv')
print(' ', OUT / 'comparison_NCA_vs_MM_paired.csv')
