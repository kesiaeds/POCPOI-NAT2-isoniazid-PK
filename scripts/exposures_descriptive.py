"""Compute NCA-style exposures (AUC0-24, Cmax, C24) per subject-occasion,
then summarize by phenotype x dosing strategy, plus paired statistics.

Does NOT require nlmixr2 — pure pandas/scipy.
"""
import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats

OUT = Path('/Users/kesiaeds/repos/POCPOI/analysis/output')
DATA = Path('/Users/kesiaeds/repos/POCPOI/analysis/data')

PHEN = {1: 'Slow', 2: 'Intermediate', 3: 'Rapid'}
OCCN = {1: 'Genotype-guided', 2: 'Flat 900 mg'}

# Supratherapeutic threshold for C24 (µg/mL) often cited in the LTBI literature
SUPRA_C24 = 0.15

d = pd.read_csv(DATA / 'pk_nonmem_combined.csv')
obs = d[d['EVID'] == 0].copy()

def per_occ(g):
    g = g.sort_values('TAD')
    auc = float(np.trapezoid(g['DV'].values, g['TAD'].values))
    cmax = float(g['DV'].max())
    c24 = float(g.loc[g['TAD'] == 24, 'DV'].iloc[0]) if (g['TAD'] == 24).any() else np.nan
    return pd.Series({'AUC24': auc, 'Cmax': cmax, 'C24': c24})

exposures = (obs.groupby(['ID', 'OCC', 'ACE_PROFILE', 'SITE'])
                .apply(per_occ).reset_index())
exposures['Phenotype'] = exposures['ACE_PROFILE'].map(PHEN)
exposures['Dosing']    = exposures['OCC'].map(OCCN)

# attach per-subject weight and the actual mg dose for that occasion (mg/kg = AMT/WT)
dose_info = (d[d['EVID'] == 1]
             .drop_duplicates(['ID', 'OCC'])
             [['ID', 'OCC', 'WT', 'AMT']]
             .rename(columns={'AMT': 'DOSE_MG'}))
exposures = exposures.merge(dose_info, on=['ID', 'OCC'], how='left')
exposures['DOSE_MGKG'] = exposures['DOSE_MG'] / exposures['WT']
exposures.to_csv(OUT / 'exposures_subject_occasion.csv', index=False)

# --- Summary: Dosing x Phenotype ------------------------------------------
def med_iqr(s):
    s = s.dropna()
    if len(s) == 0:
        return 'NA'
    return f"{s.median():.2f} ({s.quantile(.25):.2f}-{s.quantile(.75):.2f})"

summ = []
for dosing in ['Genotype-guided', 'Flat 900 mg']:
    for phen in ['Slow', 'Intermediate', 'Rapid', 'All']:
        sub = exposures[exposures['Dosing'] == dosing]
        if phen != 'All':
            sub = sub[sub['Phenotype'] == phen]
        summ.append({
            'Dosing': dosing, 'Phenotype': phen, 'N': len(sub),
            'Weight (kg) median (IQR)': med_iqr(sub['WT']),
            'Dose (mg/kg) median (IQR)': med_iqr(sub['DOSE_MGKG']),
            'AUC24 median (IQR)': med_iqr(sub['AUC24']),
            'Cmax median (IQR)':  med_iqr(sub['Cmax']),
            'C24 median (IQR)':   med_iqr(sub['C24']),
        })
summary = pd.DataFrame(summ)
summary.to_csv(OUT / 'exposure_summary_descriptive.csv', index=False)
print('=== Exposure summary (descriptive, all PK IDs) ===')
print(summary.to_string(index=False))

# --- Paired stats: Day 7 (guided) vs Day 14 (flat) per subject ------------
wide = exposures.pivot_table(index=['ID', 'Phenotype'],
                              columns='Dosing',
                              values=['AUC24', 'Cmax', 'C24']).reset_index()
# flatten
wide.columns = ['_'.join([c for c in col if c]).strip('_')
                for col in wide.columns.values]

# Wilcoxon signed-rank (paired): flat vs guided
def wilcox(a, b):
    pairs = pd.concat([a, b], axis=1).dropna()
    if pairs.shape[0] < 5:
        return ('NA', pairs.shape[0])
    r = stats.wilcoxon(pairs.iloc[:, 0], pairs.iloc[:, 1])
    return (f"{r.pvalue:.3g}", pairs.shape[0])

paired_stats = []
for metric in ['AUC24', 'Cmax', 'C24']:
    a = wide[f'{metric}_Genotype-guided']
    b = wide[f'{metric}_Flat 900 mg']
    p, n = wilcox(a, b)
    paired_stats.append({'Metric': metric, 'N_pairs': n,
                         'Median guided': f"{a.dropna().median():.3f}",
                         'Median flat':   f"{b.dropna().median():.3f}",
                         'Wilcoxon p (paired)': p})
    # by phenotype
    for phen in ['Slow', 'Intermediate', 'Rapid']:
        sub = wide[wide['Phenotype'] == phen]
        p, n = wilcox(sub[f'{metric}_Genotype-guided'],
                      sub[f'{metric}_Flat 900 mg'])
        paired_stats.append({'Metric': f'{metric} ({phen})', 'N_pairs': n,
                             'Median guided': f"{sub[f'{metric}_Genotype-guided'].dropna().median():.3f}",
                             'Median flat':   f"{sub[f'{metric}_Flat 900 mg'].dropna().median():.3f}",
                             'Wilcoxon p (paired)': p})

paired_df = pd.DataFrame(paired_stats)
paired_df.to_csv(OUT / 'paired_stats.csv', index=False)
print('\n=== Paired statistics (Wilcoxon signed-rank) ===')
print(paired_df.to_string(index=False))

# --- Variability reduction: IQR-based --------------------------------------
def iqr(s): return s.dropna().quantile(.75) - s.dropna().quantile(.25)

auc_guided_iqr = iqr(exposures.loc[exposures['Dosing'] == 'Genotype-guided', 'AUC24'])
auc_flat_iqr   = iqr(exposures.loc[exposures['Dosing'] == 'Flat 900 mg',     'AUC24'])
auc_guided_cv  = exposures.loc[exposures['Dosing'] == 'Genotype-guided', 'AUC24'].std() / \
                 exposures.loc[exposures['Dosing'] == 'Genotype-guided', 'AUC24'].mean()
auc_flat_cv    = exposures.loc[exposures['Dosing'] == 'Flat 900 mg',     'AUC24'].std() / \
                 exposures.loc[exposures['Dosing'] == 'Flat 900 mg',     'AUC24'].mean()

print(f"\nAUC0-24 IQR  : guided {auc_guided_iqr:.2f} vs flat {auc_flat_iqr:.2f} "
      f"(reduction {auc_flat_iqr/auc_guided_iqr:.2f}-fold)")
print(f"AUC0-24 CV%  : guided {100*auc_guided_cv:.1f}% vs flat {100*auc_flat_cv:.1f}%")

# --- McNemar on C24 > supratherapeutic threshold (paired binary) ----------
mc = wide.copy()
mc['supra_g'] = mc['C24_Genotype-guided'] > SUPRA_C24
mc['supra_f'] = mc['C24_Flat 900 mg']     > SUPRA_C24
mc = mc.dropna(subset=['C24_Genotype-guided', 'C24_Flat 900 mg'])

# 2x2 table
b = ((mc['supra_g'] == False) & (mc['supra_f'] == True)).sum()
c = ((mc['supra_g'] == True)  & (mc['supra_f'] == False)).sum()
a = ((mc['supra_g'] == True)  & (mc['supra_f'] == True)).sum()
dd = ((mc['supra_g'] == False) & (mc['supra_f'] == False)).sum()
print(f"\nC24 > {SUPRA_C24} µg/mL  (paired n={len(mc)}):")
print(f"  Both supra   : {a}")
print(f"  Flat only    : {b}")
print(f"  Guided only  : {c}")
print(f"  Neither      : {dd}")
if (b + c) > 0:
    chi2 = (abs(b - c) - 1)**2 / (b + c)
    p = 1 - stats.chi2.cdf(chi2, df=1)
    print(f"  McNemar chi2 = {chi2:.2f}, p = {p:.3g}")

print('\nFiles written to', OUT)
