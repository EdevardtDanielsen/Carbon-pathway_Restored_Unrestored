# Analysis and figure code

Code and input data behind the figures and statistical models of the Lake Xuanwu
restoration study (field flux survey + ¹³C-labelled mesocosm experiment).

## Contents

| File | Produces |
|---|---|
| `01_figure1_flux_by_zone.py` | **Fig. 1** — FCO₂ and FCH₄ boxplots by zone (`figures/Figure_1.*`) |
| `02_figure1_zone_letters.R` | Zone-level mixed model and Tukey compact letters used in Fig. 1 |
| `03_tableS2_treatment_model.R` | **Table S2** — nested restored-vs-unrestored mixed model |
| `04_figure2_stable_isotopes.py` | **Fig. 2** — δ¹³C time series, 10 panels (`figures/Figure_2.*`) |
| `05_figure3_partial_correlations.R` | **Fig. 3** — partial-correlation matrices, panels stacked (`figures/Figure_3.*`) |
| `05b_figure3_partial_correlations_side_by_side.R` | Fig. 3 alternate layout, panels side by side (`figures/Figure_3_side_by_side.*`) |
| `00_load_sem_data.R` | Shared data loader, sourced by `06` |
| `06_figure4_sem.R` | **Fig. 4** — piecewise SEM coefficients and Fisher's C (`sem_output/`) |
| `figure4_sem_diagram.drawio` | **Fig. 4** — the path diagram itself, drawn in draw.io from the `06` output |
| `07_tableS3_leave_one_mesocosm_out.R` | **Table S3** — leave-one-mesocosm-out path sensitivity |
| `08_tableS4_dsep_claims.R` | **Table S4** — full directed-separation claim list |
| `09_tableS1_environment_medians.py` | **Table S1** — median (IQR) environmental variables (`tables/`) |

## Data

| File | Used by |
|---|---|
| `data/field_survey.csv` | `01`, `02`, `03`, `09` — field flux survey, 16 sites × 24 sampling occasions, Sep 2023 – Feb 2024 |
| `data/mesocosm_experiment.csv` | `04`, `05`, `06`, `09` — mesocosm δ¹³C, PLFA and water chemistry |
| `data/field_water_depth.csv` | `09` — field water depth, 16 sites × 15 monthly visits |

Every data file has been reduced to the columns these scripts use.

| File | Columns |
|---|---|
| `field_survey.csv` | Position, Month, Temp (°C), DO (mg/L), FCH₄ and FCO₂ nmol m⁻² s⁻¹, Chl a (μg/L), NO₃-N, NO₂-N, PO₄-P, TN, TP (mg/L), DOC and DIC (ppm) |
| `mesocosm_experiment.csv` | Mesocosm, Day, Dominance, Label, δ¹³C of DIC, DOC, CH₄ and CO₂, zooplankton and macrophyte tissue, the seven PLFA biomarkers, and ten environmental variables |
| `field_water_depth.csv` | Position, Site, Month, Water depth (m) |

All three files are UTF-8 CSV with a single header row.

The water-depth series is monthly and runs from 2022-09 to 2024-01, which is a
longer and partly earlier window than the 24 weekly flux visits between 2023-09
and 2024-02. Table S1 pools all 220 depth measurements, so that row of the table
covers a different period from the rest of it.

## Running

Scripts resolve their own paths, so they can be run from anywhere:

```
python 01_figure1_flux_by_zone.py
Rscript 02_figure1_zone_letters.R
Rscript 03_tableS2_treatment_model.R
python 04_figure2_stable_isotopes.py
Rscript 05_figure3_partial_correlations.R
Rscript 05b_figure3_partial_correlations_side_by_side.R
Rscript 06_figure4_sem.R
Rscript 07_tableS3_leave_one_mesocosm_out.R
Rscript 08_tableS4_dsep_claims.R
python 09_tableS1_environment_medians.py
```

Figures are written to `figures/` as 600 dpi PNG, editable SVG and vector PDF.
`06` writes timestamped coefficient, Fisher's C and log files to `sem_output/`
and takes about five minutes (1000 within-cluster bootstrap refits plus MCMC per
equation, four model variants).

Fig. 1 hard-codes the compact letters that `02` computes; re-run `02` and check
them if the field data change.

Figure 3 exists in two layouts. `05` stacks the two panels vertically at
122 × 236 mm, which is the layout Su asked for in figure comment 18 so that each
panel spans the page and the 13 × 13 labels sit at a legible ~8 pt. `05b` keeps
the earlier side-by-side layout at 180 × 100 mm, where the same labels fall to
~5–6 pt. The statistics are identical in both; only the arrangement and font
sizes differ. `05` is the version to submit.

Fig. 4 is not produced by a script alone. `06` gives the coefficients, Fisher's
C and R² values; the diagram is drawn by hand in
[draw.io](https://app.diagrams.net) from `figure4_sem_diagram.drawio` and
exported to PDF/SVG.

## Requirements

Python 3.12 with `pandas`, `numpy` and `matplotlib` (see `requirements.txt`). Arial must be installed for the figures to typeset as
published.

R 4.5 with `lme4`, `lmerTest`, `emmeans`, `GGally`, `ggplot2`,
`cowplot`, `svglite`, `ragg`, `nlme`, `MuMIn`, `MCMCglmm`, `coda` and
`piecewiseSEM`.

## Licence

Code under the MIT licence, see LICENSE. Data under the Creative Commons CC0 1.0
Universal Public Domain Dedication, see data/LICENSE.

## Reproduction check

Re-run 2026-09-01, and again 2026-09-02 after the data files were reduced to the
columns actually used. Every result below is unchanged by that reduction.

Original against the figures embedded in the submitted manuscript:

- **Fig. 1** — regenerated PNG byte-identical to the published image.
- **Fig. 3** — both layouts regenerate byte-identically: the stacked `05` matches
  the 14 July figure set, and the side-by-side `05b` matches the image embedded
  in the 23 July manuscript.
- **Fig. 2** — same 415 vector paths and same 81 text elements as the published
  SVG, canvas 3.5 pt narrower. The unmodified original script reproduces this
  same shift, so it comes from font-metric changes in the current matplotlib,
  not from the code.
- **Table S2** — CH₄ p = 0.0643, CO₂ p = 0.7383, as reported.
- **Fig. 1 letters** — CH₄ R1 *b*, R2 *a*, U1 *c*, U2 *c* (F₃,₁₂.₄ = 42.95,
  p = 8.1 × 10⁻⁷); CO₂ R1 *b*, R2 *a*, U1 *b*, U2 *b* (F₃,₁₂.₂ = 12.73,
  p = 4.6 × 10⁻⁴). Both match the letters drawn in Fig. 1.
- **Fig. 4** — global fit reproduced for both panels: macrophyte-dominated
  (MACRO-zoo) Fisher's C = 74.13, df = 80, p = 0.664; phytoplankton-dominated
  (PHYTO-zoo) Fisher's C = 96.91, df = 92, p = 0.343, against 74.1/80/0.664 and
  96.9/92/0.343 as reported. Every standardized coefficient quoted in the Results
  also reproduces: δ¹³C-DIC → δ¹³C-CO₂ +1.116 macrophyte-dominated and +0.950
  phytoplankton-dominated, TN →
  δ¹³C-DIC −0.756, PO₄³⁻-P → δ¹³C-DIC −0.548, δ¹³C-zooplankton → δ¹³C-CH₄
  −0.280.
- **Table S3** — all four mesocosm-drop columns reproduce, including the
  δ¹³C-DOC → δ¹³C-CH₄ path at p = 0.028 with all mesocosms and p = 0.207 without
  TH1.
- **Table S4** — 40 claims / 0 failures (macrophyte-dominated) and 46 claims /
  4 failures (phytoplankton-dominated), matching the SI table headings.

- **Table S1** — all 55 cells reproduce the published table exactly, including
  the water-depth row at n = 100 restored and n = 120 unrestored. The original
  code for this table was never saved; `09` was reconstructed in 2026-09 and
  checked cell by cell against both the published SI table and the archived
  `Supplementary_Median_Table_Field_vs_Experiment.csv`.

Restored and unrestored refer to the FIELD zones only, positions 1–4 and 9–12 versus
5–8 and 13–16. The mesocosms were never restored and are named by dominance instead,
macrophyte-dominated (TH, CH) and phytoplankton-dominated (URT, URC).
