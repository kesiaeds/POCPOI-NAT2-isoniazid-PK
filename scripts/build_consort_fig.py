"""Draw the CONSORT-style enrollment + PK flow figure for the manuscript.
   Output: analysis/output/Figure_CONSORT.pdf

   Flow:
       Enrolled (228, CG 163 + MA 65)
           |
       Genotyped (224)  -- pending sequencing (4 Manaus IDs)
           |
       Phenotype assigned (99 slow / 108 intermediate / 21 rapid)
           |
       Intensive PK (105)  -- no intensive PK (119)
           |
       Excluded (1): MA_3 (post-dose phenotype reclassification)
           |
       Analyzed: popPK fit (104 subj / 816 obs / M3 BLQ) + paired exposure
"""
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path('/Users/kesiaeds/repos/POCPOI')
OUT  = ROOT / 'analysis' / 'output' / 'Figure_CONSORT.pdf'

fig, ax = plt.subplots(figsize=(8.5, 9.5))
ax.set_xlim(0, 10); ax.set_ylim(0, 12)
ax.set_axis_off()

MAIN_X, MAIN_W = 1.6, 5.6           # main flow column
SIDE_X, SIDE_W = 7.6, 2.2           # side-branch (exclusions) column

def main_box(y, lines, fill='#e8f0fe'):
    h = 1.05
    ax.add_patch(FancyBboxPatch(
        (MAIN_X, y), MAIN_W, h, boxstyle='round,pad=0.05,rounding_size=0.18',
        linewidth=1.4, edgecolor='black', facecolor=fill))
    n = len(lines)
    for i, (txt, bold, sz) in enumerate(lines):
        yy = y + h - 0.20 - i * (h - 0.30) / max(n - 1, 1) if n > 1 else y + h / 2
        ax.text(MAIN_X + MAIN_W / 2, yy, txt,
                ha='center', va='center',
                fontsize=sz, fontweight=('bold' if bold else 'normal'))
    return y, y + h

def side_box(y, lines, fill='#fde0d3'):
    h = 0.95
    ax.add_patch(FancyBboxPatch(
        (SIDE_X, y), SIDE_W, h, boxstyle='round,pad=0.04,rounding_size=0.15',
        linewidth=1.0, edgecolor='black', facecolor=fill))
    n = len(lines)
    for i, (txt, bold, sz) in enumerate(lines):
        yy = y + h - 0.18 - i * (h - 0.28) / max(n - 1, 1) if n > 1 else y + h / 2
        ax.text(SIDE_X + SIDE_W / 2, yy, txt,
                ha='center', va='center',
                fontsize=sz, fontweight=('bold' if bold else 'normal'))
    return y, y + h

def down_arrow(y_top, y_bot):
    ax.add_patch(FancyArrowPatch(
        (MAIN_X + MAIN_W / 2, y_top), (MAIN_X + MAIN_W / 2, y_bot),
        arrowstyle='-|>', mutation_scale=18, linewidth=1.4, color='black'))

def side_arrow(y_main, y_side):
    # arrow from main box right edge to side box left edge
    ax.add_patch(FancyArrowPatch(
        (MAIN_X + MAIN_W, y_main), (SIDE_X, y_side),
        arrowstyle='-|>', mutation_scale=14, linewidth=1.0, color='black'))

# ---- main flow: top to bottom ---------------------------------------------
y = 10.6
b1 = main_box(y, [
    ('Enrolled and consented', True, 12),
    ('n = 228', False, 11),
    ('Campo Grande 163  •  Manaus 65', False, 10),
], fill='#cfe2ff')

y_next = b1[0] - 1.55
down_arrow(b1[0], y_next + 1.05)
b2 = main_box(y_next, [
    ('NAT2 genotyping', True, 12),
    ('n = 228', False, 11),
    ('7-SNP nanopore panel (n = 224); qPCR for n = 4 Manaus', False, 9),
])

y_next = b2[0] - 1.65
down_arrow(b2[0], y_next + 1.05)
b3 = main_box(y_next, [
    ('NAT2 phenotype assigned', True, 12),
    ('n = 228', False, 11),
    ('Slow 99  •  Intermediate 108  •  Rapid 21', False, 10),
])

y_next = b3[0] - 1.85
down_arrow(b3[0], y_next + 1.05)
b4 = main_box(y_next, [
    ('Completed intensive PK sampling', True, 12),
    ('n = 105', False, 11),
    ('Campo Grande 73  •  Manaus 32', False, 10),
], fill='#fff4cc')

# side-branch: no intensive PK
side_y = b4[0] + 0.05
side_box(side_y, [
    ('No intensive PK', True, 10),
    ('n = 119', False, 9.5),
    ('Sparse-sampling cohort', False, 9),
])
side_arrow(b4[0] + 0.55, side_y + 0.55)

y_next = b4[0] - 2.55
down_arrow(b4[0], y_next + 1.05)
b5 = main_box(y_next, [
    ('Analyzed', True, 12),
    ('n = 104', False, 11),
    ('Campo Grande 73  •  Manaus 31', False, 10),
    ('Slow 43  •  Intermediate 45  •  Rapid 16', False, 10),
    ('Population PK fit: 816 INH plasma conc. (Beal M3 BLQ)', False, 9.5),
], fill='#d4edda')

# side-branch on arrow from intensive-PK to Analyzed: protocol-deviation exclusion
excl_y = b4[0] - 1.55
side_box(excl_y, [
    ('Excluded from PK analysis', True, 10),
    ('n = 1 (MA_3)', False, 9.5),
    ('Post-dose phenotype', False, 9),
    ('reclassification', False, 9),
])
# arrow taps off the middle of the down-arrow column
ax.add_patch(FancyArrowPatch(
    (MAIN_X + MAIN_W / 2, excl_y + 0.50), (SIDE_X, excl_y + 0.50),
    arrowstyle='-|>', mutation_scale=14, linewidth=1.0, color='black'))

# ---- footer / caption -----------------------------------------------------
ax.text(5.0, 0.25,
        "OCC 1 = Day 7 (genotype-guided dosing); OCC 2 = Day 14 (flat 900 mg). "
        "Plasma INH measured at 1, 2, 8, 24 h post-dose on each occasion.",
        ha='center', va='center', fontsize=9, style='italic')

plt.tight_layout()
plt.savefig(OUT, bbox_inches='tight')
plt.close()
print(f'Wrote {OUT}')
