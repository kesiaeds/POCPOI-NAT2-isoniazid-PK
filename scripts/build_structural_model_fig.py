"""Draw the structural popPK model schematic (Sarkodie-style Figure 2 equivalent):
   depot --Ka--> central <--Q--> peripheral, central --CL--> elimination.
   Saves analysis/output/Figure_StructuralModel.pdf.
"""
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path('/Users/kesiaeds/repos/POCPOI')
OUT  = ROOT / 'analysis' / 'output' / 'Figure_StructuralModel.pdf'

fig, ax = plt.subplots(figsize=(9, 3.6))
ax.set_xlim(0, 10); ax.set_ylim(0, 4)
ax.set_axis_off()

def box(x, y, w, h, label, sublabel=None, fill='#f2f2f2'):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle='round,pad=0.05,rounding_size=0.15',
        linewidth=1.5, edgecolor='black', facecolor=fill))
    ax.text(x + w/2, y + h/2 + (0.18 if sublabel else 0), label,
            ha='center', va='center', fontsize=12, fontweight='bold')
    if sublabel:
        ax.text(x + w/2, y + h/2 - 0.25, sublabel,
                ha='center', va='center', fontsize=10, style='italic')

def arrow(x0, y0, x1, y1, label, label_above=True, mut=18, ls='-'):
    ax.add_patch(FancyArrowPatch(
        (x0, y0), (x1, y1), arrowstyle='-|>', mutation_scale=mut,
        linewidth=1.5, color='black', linestyle=ls))
    mx, my = (x0 + x1) / 2, (y0 + y1) / 2
    ax.text(mx, my + (0.28 if label_above else -0.28), label,
            ha='center', va='center', fontsize=11)

# ---- compartments ----------------------------------------------------------
box(0.4, 1.5, 1.8, 1.0, 'Depot',       sublabel='(oral dose)',  fill='#e8f0fe')
box(3.6, 1.5, 2.0, 1.0, 'Central',     sublabel='V2/F',         fill='#fff4cc')
box(7.2, 1.5, 2.0, 1.0, 'Peripheral',  sublabel='V3/F',         fill='#fde0d3')

# ---- dose-in arrow into Depot ---------------------------------------------
ax.add_patch(FancyArrowPatch(
    (0.4, 3.3), (1.3, 2.55), arrowstyle='-|>', mutation_scale=18,
    linewidth=1.5, color='black'))
ax.text(0.05, 3.45, 'Oral dose\n(900 mg or weight/genotype-guided)',
        ha='left', va='center', fontsize=10)

# ---- Ka arrow Depot -> Central --------------------------------------------
arrow(2.25, 2.0, 3.55, 2.0, r'$k_a$')

# ---- Q arrows Central <-> Peripheral --------------------------------------
arrow(5.65, 2.12, 7.15, 2.12, 'Q', label_above=True)
arrow(7.15, 1.88, 5.65, 1.88, 'Q', label_above=False)

# ---- CL arrow out of Central ----------------------------------------------
ax.add_patch(FancyArrowPatch(
    (4.6, 1.45), (4.6, 0.55), arrowstyle='-|>', mutation_scale=18,
    linewidth=1.5, color='black'))
ax.text(4.85, 1.0,
        'CL = $\\theta_{CL}\\cdot(WT/70)^{0.75}\\cdot e^{\\beta_{NAT2}}\\cdot e^{\\eta_{CL}}$',
        ha='left', va='center', fontsize=10)
ax.text(4.6, 0.32, 'elimination', ha='center', va='top', fontsize=9, style='italic')

# ---- V scaling note bottom right ------------------------------------------
ax.text(8.2, 0.7,
        'V2, V3 scaled by $(WT/70)^{1.0}$\nProportional residual error',
        ha='center', va='center', fontsize=9, style='italic')

plt.tight_layout()
plt.savefig(OUT, bbox_inches='tight')
plt.close()
print(f'Wrote {OUT}')
