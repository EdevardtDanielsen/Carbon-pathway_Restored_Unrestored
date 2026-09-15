import os
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory
from matplotlib.ticker import ScalarFormatter

np.random.seed(42)

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data' / 'field_survey.csv'
data = pd.read_csv(DATA, encoding='utf-8')
data.columns = data.columns.str.strip()

CONV = 86400.0 * 12.011 / 1e6
COL_CO2 = 'FCO₂ nmol/(m⁻²·s⁻¹)'
COL_CH4 = 'FCH₄ nmol/(m⁻²·s⁻¹)'
for _c in (COL_CO2, COL_CH4):
    data[_c] = pd.to_numeric(data[_c], errors='coerce') * CONV
data['Position'] = pd.to_numeric(data['Position'], errors='coerce')

ZONE_POS = {
    'R1': [1, 2, 3, 4],
    'R2': [9, 10, 11, 12],
    'U1': [5, 6, 7, 8],
    'U2': [13, 14, 15, 16],
}
ZONE_ORDER = ['R1', 'R2', 'U1', 'U2']
ZONE_TREAT = {'R1': 'Restored', 'R2': 'Restored', 'U1': 'Unrestored', 'U2': 'Unrestored'}
TREAT_COLOR = {'Restored': '#0072B2', 'Unrestored': '#009E73'}

XPOS = {'R1': 0.0, 'R2': 1.0, 'U1': 2.4, 'U2': 3.4}
RESTORED_CENTER = (XPOS['R1'] + XPOS['R2']) / 2.0
UNREST_CENTER = (XPOS['U1'] + XPOS['U2']) / 2.0

ZONE_LETTERS = {
    COL_CH4: {'R1': 'b', 'R2': 'a', 'U1': 'c', 'U2': 'c'},
    COL_CO2: {'R1': 'b', 'R2': 'a', 'U1': 'b', 'U2': 'b'},
}

PLOT_COLUMNS = [COL_CO2, COL_CH4]
AXIS_LABELS = {
    COL_CO2: r'$\mathit{F}\mathrm{CO_2}$ (mg C m$^{-2}$ d$^{-1}$)',
    COL_CH4: r'$\mathit{F}\mathrm{CH_4}$ (mg C m$^{-2}$ d$^{-1}$)',
}

OUTDIR = ROOT / 'figures'
OUTNAME = 'Figure_1'

def zone_values(col, zone):
    sub = data[data['Position'].isin(ZONE_POS[zone])]
    return sub[col].dropna().astype(float).values

def plot():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'mathtext.fontset': 'custom',
        'mathtext.default': 'rm',
        'mathtext.rm': 'Arial',
        'mathtext.it': 'Arial:italic',
        'mathtext.bf': 'Arial:bold',
        'font.size': 10,
        'axes.labelsize': 12,
        'axes.titlesize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'svg.fonttype': 'none',
    })

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 4.4))
    panel_labels = ['(a)', '(b)']

    for idx, col in enumerate(PLOT_COLUMNS):
        ax = axes[idx]
        is_log = 'FCH' in col
        all_vals = []

        for zone in ZONE_ORDER:
            vals = zone_values(col, zone)
            all_vals.append(vals)
            x = XPOS[zone]
            color = TREAT_COLOR[ZONE_TREAT[zone]]

            ax.boxplot(
                vals, positions=[x], widths=0.62, showfliers=False,
                patch_artist=True, manage_ticks=False,
                boxprops={'facecolor': 'none', 'edgecolor': color, 'linewidth': 1.2},
                whiskerprops={'linewidth': 1.2, 'color': color},
                capprops={'linewidth': 1.2, 'color': color},
                medianprops={'color': 'black', 'linewidth': 1.2},
            )
            xj = x + np.random.uniform(-0.14, 0.14, size=vals.size)
            ax.scatter(xj, vals, s=6, color=color, alpha=0.55, linewidths=0, zorder=3)

        ticks, labels = [], []
        for zone in ZONE_ORDER:
            n = zone_values(col, zone).size
            ticks.append(XPOS[zone])
            labels.append(f'{zone}\n$\\mathit{{n}}$={n}')
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels, rotation=0)
        ax.set_xlim(-0.7, 4.1)

        flat = np.concatenate(all_vals) if all_vals else np.array([0.0])
        vmax, vmin = np.nanmax(flat), np.nanmin(flat)
        if is_log:
            ax.set_yscale('symlog', linthresh=0.5, linscale=0.3)
            ax.set_ylim(bottom=-0.05, top=vmax * 4)
            ax.set_yticks([0, 1, 10, 100])
            ax.yaxis.set_major_formatter(ScalarFormatter())
        else:
            rng = vmax - vmin
            ax.set_ylim(bottom=vmin - rng * 0.12, top=vmax + rng * 0.32)

        letters = ZONE_LETTERS[col]
        yr = np.nanmax(np.concatenate(all_vals)) - np.nanmin(np.concatenate(all_vals))
        for zone in ZONE_ORDER:
            vals = zone_values(col, zone)
            if vals.size == 0:
                continue
            top = vals.max()
            ylab = top * 1.5 if is_log else top + yr * 0.04
            ax.text(XPOS[zone], ylab, f'$\\mathit{{{letters[zone]}}}$',
                    ha='center', va='bottom', fontsize=11, color='black')
        trans = blended_transform_factory(ax.transData, ax.transAxes)

        ax.text(RESTORED_CENTER, -0.15, 'Restored', transform=trans, ha='center',
                va='top', fontsize=11, color=TREAT_COLOR['Restored'])
        ax.text(UNREST_CENTER, -0.15, 'Unrestored', transform=trans, ha='center',
                va='top', fontsize=11, color=TREAT_COLOR['Unrestored'])

        ax.set_xlabel('')
        ax.set_ylabel(AXIS_LABELS.get(col, col), fontsize=12)
        ax.grid(ls='--', alpha=0.3, linewidth=0.6)
        ax.set_axisbelow(True)
        for spine in ax.spines.values():
            spine.set_linewidth(0.8)
        ax.set_box_aspect(1)
        ax.text(-0.02, 1.04, panel_labels[idx], transform=ax.transAxes,
                fontsize=12, va='bottom', ha='left')

    fig.tight_layout(pad=0.6, w_pad=2.0)
    os.makedirs(OUTDIR, exist_ok=True)
    png = os.path.join(OUTDIR, OUTNAME + '.png')
    svg = os.path.join(OUTDIR, OUTNAME + '.svg')
    pdf = os.path.join(OUTDIR, OUTNAME + '.pdf')
    fig.savefig(png, dpi=600, bbox_inches='tight')
    fig.savefig(svg, bbox_inches='tight')
    fig.savefig(pdf, bbox_inches='tight')
    print('Saved:\n  ' + png + '\n  ' + svg + '\n  ' + pdf)

if __name__ == '__main__':
    plot()
