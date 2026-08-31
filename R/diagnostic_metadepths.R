## diagnostic_metadepths.R: sanity-check plot for the rLakeAnalyzer::meta.depths() thermocline
## detection used in TableS2_ZoneTrends.R. Faceted temperature-depth profiles for one lake/year,
## each panel with dashed lines at that cast's own metalimnion top/bottom (or none drawn if the
## cast was classified as isothermal/mixed).
suppressMessages({library(data.table); library(ggplot2); library(rLakeAnalyzer)})

LAKE <- "CR"
YEAR <- 2024

prof <- fread("data/profiles_clean.csv")[lakeid==LAKE & year4==YEAR & !is.na(wtemp)]
prof[, date := as.Date(sampledate)]

## same per-cast thermocline calc as TableS2_ZoneTrends.R
strat <- prof[, {
  o <- order(depth); dd <- depth[o]; ww <- wtemp[o]
  keep <- !duplicated(dd)
  dd <- dd[keep]; ww <- ww[keep]
  if (length(dd) < 3) {
    list(meta_top=NA_real_, meta_bot=NA_real_)
  } else {
    md <- tryCatch(meta.depths(ww, dd), error=function(e) c(NA_real_, NA_real_))
    list(meta_top=md[1], meta_bot=md[2])
  }
}, by=.(date)]

g <- ggplot(prof, aes(x=wtemp, y=depth)) +
  geom_hline(data=strat, aes(yintercept=meta_top), linetype="dashed", color="firebrick", linewidth=0.4) +
  geom_hline(data=strat, aes(yintercept=meta_bot), linetype="dashed", color="steelblue", linewidth=0.4) +
  geom_path(linewidth=0.5, color="grey30") +
  geom_point(size=0.9, color="grey30") +
  scale_y_reverse() +
  facet_wrap(~date) +
  labs(
    title=sprintf("%s %d: temperature profiles with rLakeAnalyzer::meta.depths() thermocline", LAKE, YEAR),
    subtitle="Dashed red = metalimnion top, dashed blue = metalimnion bottom (no lines = cast classified as mixed/isothermal)",
    x="Water temperature (°C)", y="Depth (m)") +
  theme_minimal(base_size=10) +
  theme(strip.text=element_text(face="bold"), plot.title=element_text(face="bold", size=12))

outfile <- sprintf("figures/diagnostics/diag_%s_%d_metadepths.png", LAKE, YEAR)
dir.create("figures/diagnostics", showWarnings=FALSE, recursive=TRUE)
ggsave(outfile, g, width=11, height=9, dpi=300, bg="white")
cat("wrote", outfile, "\n")
