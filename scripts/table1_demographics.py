"""Table 1: baseline characteristics, stratified by NAT2 acetylator phenotype.

Produces two flavours per the user's request:
  - Table 1A : ALL enrolled participants  (CG + MA, n_total)
  - Table 1B : PK subset only             (CG + MA, n_pk)
Both tables also split by site (CG vs MA) and report combined totals.
"""
import pandas as pd
from pathlib import Path

OUT = Path('/Users/kesiaeds/repos/POCPOI/analysis/output')
OUT.mkdir(parents=True, exist_ok=True)
DATA = Path('/Users/kesiaeds/repos/POCPOI/analysis/data')

PHEN  = {1: 'Slow', 2: 'Intermediate', 3: 'Rapid'}
RACE  = {1: 'Mixed', 2: 'White', 3: 'Black', 4: 'Yellow', 5: 'Indigenous'}


def _fmt_n_pct(n, total):
    return f"{int(n)} ({100*n/total:.1f}%)" if total else "0 (0.0%)"


def _fmt_med_iqr(s):
    s = s.dropna()
    if len(s) == 0:
        return "NA"
    return f"{s.median():.1f} ({s.quantile(0.25):.1f}-{s.quantile(0.75):.1f})"


def build_table(df, label):
    rows = []
    n_total = len(df)
    groups = [('Overall', df)] + [
        (PHEN[k], df[df['ACE_PROFILE'] == k])
        for k in [1, 2, 3]
    ]

    def row(name, fn):
        rows.append([name] + [fn(g) for _, g in groups])

    rows.append(['Characteristic'] + [f"{lab} (n={len(g)})" for lab, g in groups])
    row('Age, median (IQR), years',
        lambda g: _fmt_med_iqr(g['AGE']))
    row('Weight, median (IQR), kg',
        lambda g: _fmt_med_iqr(g['WT']))
    row('Sex — Male, n (%)',
        lambda g: _fmt_n_pct((g['SEX']==1).sum(), len(g)))
    row('Sex — Female, n (%)',
        lambda g: _fmt_n_pct((g['SEX']==0).sum(), len(g)))
    for race_code, race_lab in RACE.items():
        row(f'Race — {race_lab}, n (%)',
            lambda g, c=race_code: _fmt_n_pct((g['RACE']==c).sum(), len(g)))
    row('Site — Campo Grande, n (%)',
        lambda g: _fmt_n_pct((g['SITE']=='CG').sum(), len(g)))
    row('Site — Manaus, n (%)',
        lambda g: _fmt_n_pct((g['SITE']=='MA').sum(), len(g)))
    out = pd.DataFrame(rows)
    out.to_csv(OUT / f'table1_{label}.csv', index=False, header=False)
    print(f'\n=== Table 1{label.upper()}: {label} ===')
    print(out.to_string(index=False, header=False))


if __name__ == '__main__':
    allp = pd.read_csv(DATA / 'all_participants.csv')
    pks  = pd.read_csv(DATA / 'pk_subset.csv')
    build_table(allp, 'all')
    build_table(pks,  'pk')
    print(f"\nWritten to {OUT}")
