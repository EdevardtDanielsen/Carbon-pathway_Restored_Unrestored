import re
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

ROOT = Path(__file__).resolve().parent
data_path = ROOT / 'data' / 'mesocosm_experiment.csv'
df = pd.read_csv(data_path, encoding='utf-8')

def extract_group(mesocosm):
    if pd.isna(mesocosm):
        return np.nan
    match = re.match(r'^(URT|URC|TH|CH)', str(mesocosm).strip().upper())
    return match.group(1) if match else np.nan

df['Group'] = df['Mesocosm'].apply(extract_group)

def rename_columns(df):
    new_names = {}
    for col in df.columns:
        col_lower = col.lower() if isinstance(col, str) else ''
        col_str = str(col)
        if 'δ13c-dic' in col_lower or col == 'δ13C-DIC':
            new_names[col] = 'd13C_DIC'
        elif 'δ13c-doc' in col_lower or col == 'δ13C-DOC':
            new_names[col] = 'd13C_DOC'
        elif 'ch4 conc' in col_lower:
            new_names[col] = 'd13C_CH4_Conc'
        elif 'co2 conc' in col_lower:
            new_names[col] = 'd13C_CO2_Conc'
        elif 'zooplankton' in col_lower:
            new_names[col] = 'd13C_Zooplankton'
        elif 'macrophyte' in col_lower:
            new_names[col] = 'd13C_Macrophytes'
        elif '16:1' in col_str and '7c' in col_lower:
            new_names[col] = 'PLFA_16_1w7c'
        elif '16:1' in col_str and '5c' in col_lower:
            new_names[col] = 'PLFA_16_1w5c'
        elif '18:1' in col_str and '7c' in col_lower:
            new_names[col] = 'PLFA_18_1w7c'
        elif '18:3' in col_str:
            new_names[col] = 'PLFA_18_3n3'
        elif '18:2' in col_str:
            new_names[col] = 'PLFA_18_2w6'
        elif 'i15:0' in col_str or 'i15.0' in col_str:
            new_names[col] = 'PLFA_i15_0'
        elif 'a15:0' in col_str or 'a15.0' in col_str:
            new_names[col] = 'PLFA_a15_0'
    return df.rename(columns=new_names)

df = rename_columns(df)

if 'PLFA_18_3n3' in df.columns and 'PLFA_18_2w6' in df.columns:
    df['Algal_PLFA'] = df[['PLFA_18_3n3', 'PLFA_18_2w6']].mean(axis=1, skipna=True)
elif 'PLFA_18_2w6' in df.columns:
    df['Algal_PLFA'] = df['PLFA_18_2w6']

if 'PLFA_i15_0' in df.columns and 'PLFA_a15_0' in df.columns:
    df['Bacterial_PLFA'] = df[['PLFA_i15_0', 'PLFA_a15_0']].mean(axis=1, skipna=True)

if 'PLFA_16_1w7c' in df.columns and 'PLFA_16_1w5c' in df.columns:
    df['Methanotroph_I'] = df[['PLFA_16_1w7c', 'PLFA_16_1w5c']].mean(axis=1, skipna=True)
elif 'PLFA_16_1w5c' in df.columns:
    df['Methanotroph_I'] = df['PLFA_16_1w5c']

if 'PLFA_18_1w7c' in df.columns:
    df['Methanotroph_II'] = df['PLFA_18_1w7c']

treatments = {
    'TH':  {'marker': '^', 'color': '#0072B2', 'facecolor': '#0072B2', 'label': 'Macrophyte + labelled'},
    'CH':  {'marker': '^', 'color': '#0072B2', 'facecolor': 'white',   'label': 'Macrophyte + control'},
    'URT': {'marker': 'o', 'color': '#009E73', 'facecolor': '#009E73', 'label': 'Phytoplankton + labelled'},
    'URC': {'marker': 'o', 'color': '#009E73', 'facecolor': 'white',   'label': 'Phytoplankton + control'},
}

offsets = {'TH': -0.15, 'CH': -0.05, 'URT': 0.05, 'URC': 0.15}

panels = [
    (0, 0, 'd13C_DIC',          'DIC'),
    (0, 1, 'd13C_CH4_Conc',     r'CH$_4$(aq)'),
    (0, 2, 'd13C_CO2_Conc',     r'CO$_2$(aq)'),
    (1, 0, 'd13C_DOC',          'DOC'),
    (1, 1, 'd13C_Macrophytes',  'Macrophyte'),
    (1, 2, 'd13C_Zooplankton',  'Zooplankton'),
    (2, 0, 'Algal_PLFA',        'Algal PLFA'),
    (2, 1, 'Bacterial_PLFA',    'Bacterial PLFA'),
    (2, 2, 'Methanotroph_I',    'Methanotroph I'),
    (2, 3, 'Methanotroph_II',   'Methanotroph II'),
]
LEGEND_CELLS = [(0, 3), (1, 3)]
PANEL_LETTERS = 'abcdefghij'

def calc_stats(data, var):
    stats = data.groupby(['Group', 'Day'])[var].agg(['mean', 'std', 'count']).reset_index()
    stats['se'] = stats['std'] / np.sqrt(stats['count'])
    return stats

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'mathtext.default': 'regular',
    'font.size': 10,
    'axes.labelsize': 12,
    'axes.titlesize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 9,
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'svg.fonttype': 'none',
})

YLABEL = r'$\delta^{13}$C (‰)'

def get_name_pos(row_idx):
    if row_idx == 0:
        return (0.95, 0.95, 'right', 'top')
    elif row_idx == 1:
        return (0.05, 0.95, 'left', 'top')
    else:
        return (0.05, 0.05, 'left', 'bottom')

bottom_row_of_col = {}
for (r, c, _, _) in panels:
    bottom_row_of_col[c] = max(bottom_row_of_col.get(c, -1), r)

PLFA_ROW = 2

def _var_extent(var):
    if var not in df.columns:
        return None
    s = calc_stats(df, var)
    lo = (s['mean'] - s['se'].fillna(0)).min()
    hi = (s['mean'] + s['se'].fillna(0)).max()
    return float(lo), float(hi)

_ext = [e for e in (_var_extent(v) for (r, c, v, n) in panels if r == PLFA_ROW) if e]
_ylo = min(e[0] for e in _ext)
_yhi = max(e[1] for e in _ext)
_pad = (_yhi - _ylo) * 0.08
PLFA_YLIM = (_ylo - _pad, _yhi + _pad)

fig, axes = plt.subplots(3, 4, figsize=(7.0, 5.1), sharex=True)

legend_handles = None

for i, (r, c, var, name) in enumerate(panels):
    ax = axes[r, c]
    letter = f'({PANEL_LETTERS[i]})'
    nx, ny, ha, va = get_name_pos(r)

    ax.text(-0.02, 1.04, letter, transform=ax.transAxes,
            fontsize=12, fontweight='normal', va='bottom', ha='left')

    if var not in df.columns:
        ax.text(nx, ny, name, transform=ax.transAxes, fontsize=9, va=va, ha=ha)
        ax.text(0.5, 0.5, 'No data', ha='center', va='center',
                transform=ax.transAxes, fontsize=9, color='gray')
    else:
        stats = calc_stats(df, var)
        for group, style in treatments.items():
            group_data = stats[stats['Group'] == group]
            if group_data.empty:
                continue
            x = group_data['Day'] + offsets[group]
            y = group_data['mean']
            yerr = group_data['se']

            ax.errorbar(x, y, yerr=yerr, color=style['color'], linewidth=1.0,
                        linestyle='-', capsize=2, capthick=0.8, alpha=0.55,
                        marker='none')
            ax.plot(x, y, marker=style['marker'], color=style['color'],
                    markerfacecolor=style['facecolor'],
                    markeredgecolor=style['color'], markeredgewidth=1.0,
                    markersize=5, linestyle='none',
                    label=style['label'] if i == 0 else None)

        ax.text(nx, ny, name, transform=ax.transAxes, fontsize=9, va=va, ha=ha)

    ax.set_xlim(0.5, 6.5)
    ax.set_xticks([1, 2, 3, 4, 5, 6])

    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
    ax.tick_params(axis='both', width=0.8, length=3)

    ax.grid(ls='--', alpha=0.3, linewidth=0.6)
    ax.set_axisbelow(True)

    ax.yaxis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=3))

    if c == 0:
        ax.set_ylabel(YLABEL)

    if r == PLFA_ROW:
        ax.set_ylim(PLFA_YLIM)
        if c != 0:
            ax.tick_params(labelleft=False)

    if r == bottom_row_of_col[c]:
        ax.set_xlabel('Day')
    else:
        ax.tick_params(labelbottom=False)

    if i == 0:
        legend_handles = ax.get_legend_handles_labels()

for (r, c) in LEGEND_CELLS:
    axes[r, c].axis('off')

plt.subplots_adjust(left=0.08, bottom=0.09, right=0.985, top=0.95,
                    wspace=0.33, hspace=0.22)

ax_top = axes[LEGEND_CELLS[0]]
ax_mid = axes[LEGEND_CELLS[1]]
p_top = ax_top.get_position()
p_mid = ax_mid.get_position()
cx = (p_top.x0 + p_top.x1) / 2
cy = (p_mid.y0 + p_top.y1) / 2
handles, labels = legend_handles
# right-align the legend frame to the column edge so it never extends past panel (j),
# which sits in the same grid column
fig.legend(handles, labels, loc='center right', bbox_to_anchor=(p_top.x1, cy),
           bbox_transform=fig.transFigure, frameon=True, edgecolor='0.8',
           framealpha=1.0, borderpad=0.4, labelspacing=0.5, handletextpad=0.4,
           handlelength=0.9, borderaxespad=0.0)

OUTDIR = ROOT / 'figures'
OUTDIR.mkdir(exist_ok=True)
OUTNAME = 'Figure_2'

png_path = OUTDIR / (OUTNAME + '.png')
svg_path = OUTDIR / (OUTNAME + '.svg')
pdf_path = OUTDIR / (OUTNAME + '.pdf')
fig.savefig(png_path, dpi=600, bbox_inches='tight', facecolor='white')
fig.savefig(svg_path, bbox_inches='tight', facecolor='white')
fig.savefig(pdf_path, bbox_inches='tight', facecolor='white')
print(f"Saved:\n  {png_path}\n  {svg_path}\n  {pdf_path}")

