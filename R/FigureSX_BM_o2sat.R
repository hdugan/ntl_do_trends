## Small standalone panel: Big Muskellunge, DO % saturation, 1996-05 era.
## Replicates fig03's DO-sat panel exactly -- same doy window, 0.5 m interpolation, half-month
## bins, MIN_YEARS=4 blanking rule and palette -- but for one lake and one era, sized 3 x 2 in.
## No significance outlines: in fig03 those are drawn on the FINAL era (2016-25) only.
suppressMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

MIN_YEARS <- 4
LAKE <- "BM"
ERA <- "1996–05"
brk <- c(1985, 1995, 2005, 2015, 2025) # northern-lake era breaks (fig03)
labs4 <- c("1986–95", "1996–05", "2006–15", "2016–25")

prof <- fread("data/profiles_clean.csv")[lakeid == LAKE]
prof[, `:=`(date = as.Date(sampledate), doy = yday(as.Date(sampledate)), year = year4)]
prof <- prof[doy >= 105 & doy <= 310]
prof <- prof[is.finite(o2sat)]

d <- prof[, .(val = mean(o2sat)), by = .(date, doy, year, depth)]
g <- d[,
  {
    o <- order(depth)
    dd <- depth[o]
    yy <- val[o]
    if (length(dd) >= 3 && diff(range(dd)) >= 1) {
      zg <- seq(0, max(dd), 0.5)
      zg <- zg[zg >= min(dd)]
      .(dbin = zg, val = approx(dd, yy, zg, rule = 2)$y)
    } else {
      .(dbin = numeric(0), val = numeric(0))
    }
  },
  by = .(date, doy, year)
]

hedges <- c(91, 106, 121, 136, 152, 167, 182, 197, 213, 228, 244, 259, 274, 289, 305)
g[, bi := findInterval(doy, hedges)]
g <- g[bi >= 1 & bi <= length(hedges) - 1]
g[, `:=`(tbin = (hedges[bi] + hedges[bi + 1]) / 2, twidth = hedges[bi + 1] - hedges[bi])]
g[, eidx := cut(year, breaks = brk, labels = FALSE, right = TRUE)]
g <- g[!is.na(eidx)][, era := factor(labs4[eidx], levels = labs4)]

## identical blanking rule to fig03: grey where this era has <MIN_YEARS, dropped where no era does
agg <- g[, .(m = mean(val), n_yr = uniqueN(year)), by = .(eidx, era, tbin, twidth, dbin)]
agg[, ever_ok := any(n_yr >= MIN_YEARS), by = .(tbin, twidth, dbin)]
agg <- agg[ever_ok == TRUE][, ever_ok := NULL]
agg[n_yr < MIN_YEARS, m := NA]
agg <- agg[era == ERA]

pal_sat <- scale_fill_gradientn(
  colours = c("#000000", "#7f0000", "#d7301f", "#fdae61", "#ffffbf", "#66bd63", "#1a9850"),
  values = rescale(c(0, 10, 20, 40, 70, 95, 140)), name = "DO (% sat)",
  limits = c(0, 140), oob = squish, na.value = "grey80"
)
mbound <- c(91, 121, 152, 182, 213, 244, 274, 305)

## NB: fig03 draws white month-boundary rules here; omitted in this small panel by request,
## so the seasonal progression reads as a continuous surface.
p <- ggplot(agg, aes(tbin, dbin)) +
  geom_tile(aes(fill = m, width = twidth), height = 0.5) +
  scale_y_reverse(expand = c(0, 0)) +
  scale_x_continuous(
    breaks = c(106, 136.5, 167, 197.5, 228.5, 259, 289.5),
    labels = c("Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"),
    minor_breaks = NULL, expand = c(0, 0)
  ) +
  pal_sat +
  labs(x = NULL, y = "Depth (m)") +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(), legend.position = "none",
    plot.margin = margin(2, 3, 1, 2)
  )

ggsave("figures/figSD_BM_o2sat_1996-05.png", p, width = 3, height = 2, dpi = 600, bg = "white")
cat("wrote figures/figSD_BM_o2sat_1996-05.png\n")

## same panel with a compact legend, in case it is needed standalone
p2 <- p + theme(
  legend.position = "right", legend.key.width = unit(0.18, "cm"),
  legend.key.height = unit(0.42, "cm"), legend.title = element_text(size = 5.4),
  legend.text = element_text(size = 5), legend.margin = margin(0, 0, 0, 1)
)
ggsave("figures/figSD_BM_o2sat_1996-05_legend.png", p2, width = 3, height = 2, dpi = 600, bg = "white")
cat("wrote figures/figSD_BM_o2sat_1996-05_legend.png\n")

cat(
  "\ncells drawn:", nrow(agg), " | greyed (<", MIN_YEARS, "yrs):", sum(is.na(agg$m)),
  " | depth range:", min(agg$dbin), "-", max(agg$dbin), "m",
  " | years contributing:", paste(range(g[era == ERA, year]), collapse = "-"), "\n"
)
