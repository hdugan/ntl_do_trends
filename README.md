# ntl_do_trends

Selected figures and their generating R scripts from the North Temperate Lakes
long-term record analysis (NTL-LTER, EDI packages 1, 29, 35). Split out from
the main `BuoyExploration_Mendota` project.

Run in order: `00_pull_data.R` pulls the newest revision of each needed EDI
data package (via `EDIutils`) into a local `data/` folder (gitignored, not
included here); `01_plot_profiles.R` plots every individual profile so
sensor glitches are visible; `02_clean_profiles.R` removes the ones that
are physically implausible, writing `data/profiles_clean.csv` (the raw
`data/profiles.csv` is left untouched). The fig*.R scripts then read
`data/profiles_clean.csv` and write their PNG(s) to `figures/`.
`Figure2_clarity_trends.R` must run before `Figure1_rate_of_change.R`,
which reads `figures/fig02_trends_table.csv` for its DOC/Secchi arrows.

| Script | Figure(s) |
|---|---|
| `00_pull_data.R` | (no figure) — downloads `data/profiles.csv`, `chem_north.csv`, `secchi.csv`, `color.csv`, `chl.csv`, `lake_levels.csv` from EDI (NTL-LTER packages 29, 1, 31, 87, 35, 30). `lake_levels.csv` covers the 7 northern lakes only (AL, BM, CB, CR, SP, TB, TR) — not currently used by any figure script |
| `01_plot_profiles.R` | `figures/diagnostics/diag_<lake>_<bin>_profiles.png` — every individual temp & DO(%sat) profile, per lake x half-month bin (same bins as fig12), baseline (pre-2016) vs recent (2016-2025) |
| `02_clean_profiles.R` | (no figure) — for lakes with zmax > 10 m, removes unrealistic bottom-point jumps (temp up >5°C or DO up >20 %sat at the single deepest reading vs the next-shallowest); plus manual removals: one bad-sensor DO profile (Sparkling 2004-08-16) and all temp/DO for Crystal Lake 2012-2013 (whole-lake mixing experiment); writes `data/profiles_clean.csv` and logs removed points to `data/profiles_removed_points.csv` |
| `FigureSX_do_heatmaps.R` | `figSX_northern_do/sat.png`, `figSX_southern_do/sat.png` — depth-time DO heatmaps, northern vs southern lakes |
| `Figure2_clarity_trends.R` | `fig02_clarity_trends.png` — standardized (robust z-score, median/MAD) water-clarity trends, all 11 lakes: Kd/DOC/Secchi for the 7 northern lakes, DOC/Secchi only for the 4 southern lakes (no Kd data there). CDOM colour (EDI 87, a440) is deliberately excluded. Also writes `figures/fig02_trends_table.csv` (lake x metric slopes with p-values) |
| `Figure1_rate_of_change.R` | `fig01_rate_mgL.png`, `fig01_rate_sat.png` — Theil-Sen trend heatmaps by depth & season, all 11 lakes. The DOC/Secchi arrows beside each lake name are read from `figures/fig02_trends_table.csv` rather than fit here — **run `Figure2_clarity_trends.R` first** |
| `Figure3_o2_profiles.R` | `fig03_o2_profiles.png` — August DO (% saturation) depth profile, all 11 lakes, four overlaid decade medians (1985-1995 / 1995-2005 / 2005-2015 / 2015-2025, blue->green->orange->red). No annotations (no anoxic-boundary markers, no metalimnetic-O2-max diamonds) — just the raw profile shape across decades. Panels built individually + assembled with patchwork (not facet_wrap) so the gap left by the 7th northern lake is a true blank space and the legend is extracted once and placed manually. Independent of the other two `Figure*.R` scripts |
| `Figure4_chl_profiles.R` | `fig04_chl_profiles.png` — August chlorophyll profile shift (1995-2005 vs 2016-2025), 4 northern deep oligotrophic lakes, one row |
| `FigureSD_BM_o2sat.R` | `figSD_BM_o2sat_1996-05.png` (+ legend) — Big Muskellunge O2 saturation panel |

`Figure1_rate_of_change.R` (6.5x8in), `Figure2_clarity_trends.R` (6.5x6in),
`Figure3_o2_profiles.R` (6.5x6in), and `Figure4_chl_profiles.R` (6.5x2.7in), all 500dpi, don't
draw an overall title/explainer on their PNGs — that text is written instead to the shared
`figures/captions.csv` (columns: file, title, caption; each script merge-writes its own rows
without clobbering the others').
