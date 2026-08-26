"""POCPOI manuscript figures.

Outputs PNG + PDF to output/figures/:
  Figure 1  : INH concentration-time profiles by dosing strategy and NAT2 phenotype
  Figure 2  : AUC0-24 by NAT2 phenotype and dosing strategy
  Figure 3  : Within-subject paired AUC0-24 (slopegraph)
  Figure 4  : C24 by NAT2 phenotype and dosing strategy
  Figure S3 : Individual INH concentration-time profiles (supplementary)
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'data'
OUT  = ROOT / 'output'
FIG  = OUT / 'figures'
FIG.mkdir(parents=True, exist_ok=True)

# ---------- style -------------------------------------------------------------
plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size':   10,
    'axes.spines.top':   False,
    'axes.spines.right': False,
    'axes.labelsize': 11,
    'axes.titlesize': 11,
    'legend.frameon': False,
    'figure.dpi':     120,
})
PHEN_ORDER  = ['Slow', 'Intermediate', 'Rapid']
PHEN_COLOUR = {'Slow': '#0072B2', 'Intermediate': '#009E73', 'Rapid': '#D55E00'}
DOSE_COLOUR = {'Genotype-guided': '#4c72b0', 'Flat 900 mg': '#dd8452'}

def _save(fig, name):
    fig.savefig(FIG / f'{name}.png', dpi=300, bbox_inches='tight')
    fig.savefig(FIG / f'{name}.pdf', bbox_inches='tight')
    plt.close(fig)

# ---------- data --------------------------------------------------------------
nm  = pd.read_csv(DATA / 'pk_nonmem_combined.csv')
obs = nm[nm['EVID'] == 0].copy()
obs['Phenotype'] = obs['ACE_PROFILE']   # text values in public dataset
obs['Dosing']    = obs['OCC'].map({1: 'Genotype-guided', 2: 'Flat 900 mg'})

exps = pd.read_csv(OUT / 'exposures_subject_occasion.csv')

SUPRA       = 0.15
PHEN_GROUPS = ['All'] + PHEN_ORDER

legendA = [mpatches.Patch(facecolor=DOSE_COLOUR['Flat 900 mg'],       alpha=.7, label='Flat 900 mg'),
           mpatches.Patch(facecolor=DOSE_COLOUR['Genotype-guided'],   alpha=.7, label='Genotype-guided')]
legendB = legendA + [plt.Line2D([0], [0], color='red', linestyle='--', linewidth=1.2,
                                label=f'Supratherapeutic (>{SUPRA} µg/mL)')]

def _boxes(ax, metric):
    positions, data, colours = [], [], []
    pos = 1
    for group in PHEN_GROUPS:
        for dosing in ['Flat 900 mg', 'Genotype-guided']:
            if group == 'All':
                v = exps.loc[exps['Dosing'] == dosing, metric]
            else:
                v = exps.loc[(exps['Phenotype'] == group) & (exps['Dosing'] == dosing), metric]
            data.append(v.dropna().values)
            positions.append(pos); pos += 1
            colours.append(DOSE_COLOUR[dosing])
        pos += 1
    bp = ax.boxplot(data, positions=positions, widths=0.7, patch_artist=True,
                    medianprops=dict(color='black', linewidth=1.5),
                    flierprops=dict(marker='o', markersize=3, alpha=.6))
    for patch, c in zip(bp['boxes'], colours):
        patch.set_facecolor(c); patch.set_alpha(.7)
    for i, vals in enumerate(data):
        jitter = np.random.uniform(-0.12, 0.12, size=len(vals))
        ax.scatter(np.full(len(vals), positions[i]) + jitter, vals,
                   color='black', s=8, alpha=.35, zorder=3)
    group_centres = [np.mean(positions[i*2:i*2+2]) for i in range(len(PHEN_GROUPS))]
    ax.set_xticks(group_centres)
    ax.set_xticklabels(PHEN_GROUPS)
    ax.axvline(positions[1] + 0.65, color='grey', linewidth=0.8, linestyle='--', alpha=0.5)
    ax.grid(axis='y', alpha=.3)

# ============ Figure 1 : conc-time by dosing arm ==============================
fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharey=True)
for ax, dosing in zip(axes, ['Flat 900 mg', 'Genotype-guided']):
    sub = obs[obs['Dosing'] == dosing]
    for phen in PHEN_ORDER:
        ss = sub[sub['Phenotype'] == phen]
        summary = (ss.groupby('TAD')['DV']
                     .agg(median='median', q1=lambda x: x.quantile(.25),
                          q3=lambda x: x.quantile(.75)).reset_index())
        ax.fill_between(summary['TAD'], summary['q1'], summary['q3'],
                        alpha=.20, color=PHEN_COLOUR[phen])
        ax.plot(summary['TAD'], summary['median'], '-o',
                color=PHEN_COLOUR[phen], label=phen,
                linewidth=1.8, markersize=4)
    ax.set_title(dosing)
    ax.set_xlabel('Time after dose (h)')
    ax.set_xticks([1, 2, 8, 24])
    ax.set_yscale('log')
    ax.set_ylim(0.005, 30)
    ax.grid(axis='y', which='major', alpha=.3)
    ax.legend(loc='upper right', fontsize=9, title='NAT2 phenotype')
axes[0].set_ylabel('INH plasma conc. (µg/mL)')
fig.tight_layout()
_save(fig, 'Figure1')

# ============ Figure 2 : AUC0-24 boxplots =====================================
fig, ax = plt.subplots(figsize=(8, 4.8))
_boxes(ax, 'AUC24')
ax.set_ylabel('AUC$_{0-24}$ (mg·h/L)')
ax.legend(handles=legendA, loc='upper right')
fig.tight_layout()
_save(fig, 'Figure2')

# ============ Figure 3 : within-subject paired AUC0-24 slopegraph ============
wide = exps.pivot_table(index=['ID', 'Phenotype'], columns='Dosing', values='AUC24').reset_index()
wide = wide.dropna(subset=['Flat 900 mg', 'Genotype-guided'])
fig, axes = plt.subplots(1, 3, figsize=(11, 4.2), sharey=True)
for ax, phen in zip(axes, PHEN_ORDER):
    sub = wide[wide['Phenotype'] == phen]
    for _, row in sub.iterrows():
        ax.plot([0, 1], [row['Flat 900 mg'], row['Genotype-guided']],
                '-o', color=PHEN_COLOUR[phen], alpha=.45, linewidth=1, markersize=4)
    med_f = sub['Flat 900 mg'].median()
    med_g = sub['Genotype-guided'].median()
    ax.plot([0, 1], [med_f, med_g], '-', color='black', linewidth=3,
            label=f'median: {med_f:.1f} → {med_g:.1f}')
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Flat\n900 mg', 'Genotype-\nguided'])
    ax.set_title(phen)
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(axis='y', alpha=.3)
axes[0].set_ylabel('AUC$_{0-24}$ (mg·h/L)')
fig.tight_layout()
_save(fig, 'Figure3')

# ============ Figure 4 : C24 boxplots ========================================
fig, ax = plt.subplots(figsize=(8, 4.8))
_boxes(ax, 'C24')
ax.axhline(SUPRA, color='red', linestyle='--', linewidth=1.2)
ax.set_yscale('log')
ax.set_ylim(0.003, 2)
ax.set_ylabel('C$_{24}$ (µg/mL)')
ax.legend(handles=legendB, loc='upper right', fontsize=9)
fig.tight_layout()
_save(fig, 'Figure4')

# ============ Figure S3 : individual conc-time profiles (supplementary) ======
fig, axes = plt.subplots(2, 3, figsize=(13, 8), sharey=True)
for col, phen in enumerate(PHEN_ORDER):
    for row, dosing in enumerate(['Flat 900 mg', 'Genotype-guided']):
        ax = axes[row][col]
        ss = obs[(obs['Phenotype'] == phen) & (obs['Dosing'] == dosing)]
        for subj_id, grp in ss.groupby('ID'):
            grp_s = grp.sort_values('TAD')
            ax.plot(grp_s['TAD'], grp_s['DV'], '-o',
                    color=PHEN_COLOUR[phen], alpha=0.35, linewidth=0.8, markersize=3)
        summary = (ss.groupby('TAD')['DV']
                     .agg(median='median').reset_index())
        ax.plot(summary['TAD'], summary['median'], '-',
                color='black', linewidth=2.2, label='Median', zorder=5)
        ax.set_title(f'{phen}\n{dosing}', fontsize=9)
        ax.set_xlabel('Time after dose (h)')
        ax.set_xticks([1, 2, 8, 24])
        ax.set_yscale('log')
        ax.set_ylim(0.005, 30)
        ax.grid(axis='y', which='major', alpha=.3)
        if col == 0:
            ax.set_ylabel('INH plasma conc. (µg/mL)')
axes[0][0].legend(fontsize=8)
fig.tight_layout()
_save(fig, 'Figure_S3')

print('Figures written to', FIG)
for f in sorted(FIG.glob('*.png')):
    print(' ', f.name)
