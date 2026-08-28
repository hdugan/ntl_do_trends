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

| Script | Figure(s) |
|---|---|
| `00_pull_data.R` | (no figure) — downloads `data/profiles.csv`, `chem_north.csv`, `secchi.csv`, `color.csv`, `chl.csv` from EDI (NTL-LTER packages 29, 1, 31, 87, 35) |
| `01_plot_profiles.R` | `figures/diagnostics/diag_<lake>_<bin>_profiles.png` — every individual temp & DO(%sat) profile, per lake x half-month bin (same bins as fig12), baseline (pre-2016) vs recent (2016-2025) |
| `02_clean_profiles.R` | (no figure) — for lakes with zmax > 10 m, removes unrealistic bottom-point jumps (temp up >5°C or DO up >20 %sat at the single deepest reading vs the next-shallowest), plus one manually-flagged bad-sensor DO profile (Sparkling 2004-08-16); writes `data/profiles_clean.csv` and logs removed points to `data/profiles_removed_points.csv` |
| `fig03_04_depthtime.R` | `fig03_northern_do/sat.png`, `fig04_southern_do/sat.png` — depth-time DO heatmaps, northern vs southern lakes |
| `fig07_alllakes_metO2_aug.R` | `fig07_alllakes_metO2_aug.png` — metalimnetic O2 maxima, August, all lakes |
| `fig11_light_attenuation.R` | `fig11_light_attenuation.png` — light attenuation vs DOC/Secchi |
| `fig12_rate_of_change.R` | `fig12_rate_mgL.png`, `fig12_rate_sat.png` — Theil-Sen trend heatmaps by depth & season, all 11 lakes |
| `fig15_chl_profiles_aug.R` | `fig15_chl_profiles_aug.png` — August chlorophyll profiles |
| `figS_BM_o2sat_panel.R` | `figS_BM_o2sat_1996-05.png` (+ legend) — Big Muskellunge O2 saturation panel |
