import csv
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data'
OUTDIR = ROOT / 'tables'
OUTNAME = 'Table_S1_environment_medians.csv'

RESTORED_POS = [1, 2, 3, 4, 9, 10, 11, 12]
# field zones are restored or unrestored; the mesocosms are named by dominance
FIELD_SYSTEMS = ['Restored', 'Unrestored']
MESO_SYSTEMS = ['macrophyte', 'phytoplankton']

ROWS = [
    ('DIC (ppm)',       'DIC (ppm)',    'DIC (ppm)'),
    ('DOC (ppm)',       'DOC (ppm)',    'DOC (ppm)'),
    ('TN (mg/L)',       'TN (mg/L)',    'TN (mg/L)'),
    ('Chl a (µg/L)',    'Chl a (μg/L)', 'Chl a (μg/L)'),
    ('Temp (°C)',       'Temp (°C)',    'Temp (°C)'),
    ('DO (mg/L)',       'DO (mg/L)',    'DO (mg/L)'),
    ('TP (mg/L)',       'TP (mg/L)',    'TP (mg/L)'),
    ('Water depth (m)', None,           None),
    ('PO4-P (mg/L)',    'PO₄-P (mg/L)', 'PO₄-P (mg/L)'),
    ('NO2-N (mg/L)',    'NO₂-N (mg/L)', 'NO₂-N (mg/L)'),
    ('NO3-N (mg/L)',    'NO₃-N (mg/L)', 'NO₃-N (mg/L)'),
]

DEPTH_COL = 'Water depth (m)'


def load_field():
    df = pd.read_csv(DATA / 'field_survey.csv', encoding='utf-8')
    df['Position'] = pd.to_numeric(df['Position'], errors='coerce')
    df['System'] = np.where(df['Position'].isin(RESTORED_POS), 'Restored', 'Unrestored')
    return df


def load_mesocosm():
    df = pd.read_csv(DATA / 'mesocosm_experiment.csv', encoding='utf-8')
    df['System'] = df['Dominance']
    return df


def load_depth():
    df = pd.read_csv(DATA / 'field_water_depth.csv', encoding='utf-8')
    df['System'] = np.where(df['Position'].isin(RESTORED_POS), 'Restored', 'Unrestored')
    return df


def summarise(frame, column, system):
    if column is None or column not in frame.columns:
        return 'NA'
    s = pd.to_numeric(frame.loc[frame['System'] == system, column], errors='coerce').dropna()
    if s.empty:
        return 'NA'
    return f'{s.median():.2f} ({s.quantile(0.25):.2f}-{s.quantile(0.75):.2f}); n={s.size}'


def build():
    field, meso, depth = load_field(), load_mesocosm(), load_depth()
    table = [['Variable (unit)', 'Field Restored', 'Field Unrestored',
              'Mesocosm Macrophyte-dominated', 'Mesocosm Phytoplankton-dominated']]
    for label, field_col, meso_col in ROWS:
        if label == DEPTH_COL:
            row = [label] + [summarise(depth, DEPTH_COL, s) for s in FIELD_SYSTEMS] + ['NA', 'NA']
        else:
            row = ([label]
                   + [summarise(field, field_col, s) for s in FIELD_SYSTEMS]
                   + [summarise(meso, meso_col, s) for s in MESO_SYSTEMS])
        table.append(row)
    return table


if __name__ == '__main__':
    table = build()
    OUTDIR.mkdir(exist_ok=True)
    path = OUTDIR / OUTNAME
    with open(path, 'w', newline='', encoding='utf-8-sig') as fh:
        csv.writer(fh).writerows(table)
    for row in table:
        print(' | '.join(row))
    print('\nSaved: ' + str(path))
