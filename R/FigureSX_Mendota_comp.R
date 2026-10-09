## FigureSX_Mendota_comp.R: Lake Mendota temperature and dissolved oxygen, 1906 vs 2016-2025.
##
## Four panels -- rows = era (1906, 2016-2025), columns = variable (temperature, DO mg/L).
## The 1906 data are digitized from Birge & Juday (1911), Wisconsin Geol. Nat. Hist. Surv.
## Bulletin XXII, Plates I (temperature) and II (dissolved oxygen):
##   data/Juday/Mendota_1906_temp.csv      -- 9 fixed depths
##   data/Juday/Mendota_1906_DO_plate.csv  -- 11 fixed depths; Plate II is in cc O2/L, so the
##                                            do_mgL column (x 1.429) is used here
## Repeated casts at a set of FIXED depths across the season, not continuous profiles. To keep
## the eras directly comparable, each variable's modern panel is built at exactly that
## variable's 1906 depths (each 2016-2025 cast interpolated to them) rather than the fine
## continuous depth grid FigureSX_do_heatmaps.R uses elsewhere -- showing modern data at a finer
## resolution than 1906 ever had would imply a precision the historical panel can't match.
## Same half-month time bins and DO colour scale as FigureSX_do_heatmaps.R.
suppressMessages({library(data.table); library(ggplot2); library(scales); library(patchwork)})

MIN_YEARS <- 4   # modern cells need data from at least this many distinct years to be averaged

## ---- 1906 (digitized plates) --------------------------------------------------------------------
t1906 <- fread("data/Juday/Mendota_1906_temp.csv")[, .(doy=yday, depth=depth_m, val=temp_C)]
o1906 <- fread("data/Juday/Mendota_1906_DO_plate.csv")[, .(doy=yday, depth=depth_m, val=do_mgL)]
DEPTHS <- list(temp=sort(unique(t1906$depth)), do=sort(unique(o1906$depth)))

## ---- 2016-2025: interpolate each cast to the 1906 depths for that variable -------------------
prof <- fread("data/profiles_clean.csv")[lakeid=="ME"]
prof[, `:=`(doy=yday(as.Date(sampledate)), year=year4)]
prof <- prof[year>=2016 & year<=2025]

## only interpolate a target depth if it's within 2 m of the cast's own sampled range -- don't
## extrapolate a 22 m estimate from a cast that only reached 12 m
interp_casts <- function(var, depths) {
  out <- prof[!is.na(get(var)), {
    o <- order(depth); dd <- depth[o]; yy <- get(var)[o]
    if (length(dd)>=2) {
      est <- approx(dd, yy, depths, rule=2, ties=mean)$y
      est[!(depths >= min(dd)-2 & depths <= max(dd)+2)] <- NA_real_
      .(depth=depths, val=est)
    } else .(depth=numeric(0), val=numeric(0))
  }, by=.(date=as.Date(sampledate), doy, year)]
  out[!is.na(val)]
}
## ---- % saturation, computed the SAME way for both eras -------------------------------------
## Benson & Krause (1984) freshwater solubility x barometric correction at Mendota's own
## elevation (258 m, EDI 434). Not NTL-LTER's o2sat column: back-solving it for Mendota shows a
## constant +2.9% offset, i.e. it carries the ~495 m correction of the northern lakes.
o2_solubility_mgL <- function(tempC) {
  Ts <- tempC + 273.15
  exp(-139.34411 + 1.575701e5/Ts - 6.642308e7/Ts^2 + 1.2438e10/Ts^3 - 8.621949e11/Ts^4)
}
ME_ELEV <- fread("data/lake_characteristics.csv")[waterbody_name=="Lake Mendota", elevation_m]
cs_local <- function(tempC) o2_solubility_mgL(tempC) * (1 - 2.25577e-5*ME_ELEV)^5.25588
prof[, o2sat_me := 100 * o2 / cs_local(wtemp)]

## 1906: temperature from Plate I at each Plate II cast -- interpolated in time at each
## temperature depth, then in depth (6 and 7 m sit between 5 and 8 m). Only within 1.5 d of the
## temperature record's span, so the November DO casts after Plate I ends get no %sat.
t_split <- split(t1906, by="depth", keep.by=FALSE)
t_depths <- as.numeric(names(t_split))
temp_at <- function(d0, z0) {
  Ts <- vapply(t_split, function(s) {
    if (d0 < min(s$doy) - 1.5 || d0 > max(s$doy) + 1.5) NA_real_
    else approx(s$doy, s$val, d0, rule=2)$y
  }, numeric(1))
  if (anyNA(Ts)) return(NA_real_)
  approx(t_depths, Ts, z0, rule=2)$y
}
s1906 <- copy(o1906)[, wtemp := mapply(temp_at, doy, depth)]
s1906 <- s1906[!is.na(wtemp), .(doy, depth, val = 100 * val / cs_local(wtemp))]
cat(sprintf("1906 surface %%sat, median over season: %.0f%%\n", median(s1906[depth==0, val])))

tmod <- interp_casts("wtemp", DEPTHS$temp)
omod <- interp_casts("o2", DEPTHS$do)
smod <- interp_casts("o2sat_me", DEPTHS$do)

## ---- half-month time bins, identical edges to FigureSX_do_heatmaps.R / Figure1 ---------------
hedges <- c(91,106,121,136,152,167,182,197,213,228,244,259,274,289,305,320)
mbound <- c(91,121,152,182,213,244,274,305,320)
xsc <- scale_x_continuous(breaks=c(106,136.5,167,197.5,228.5,259,289.5,312.5),
                          labels=c("Apr","May","Jun","Jul","Aug","Sep","Oct","Nov"),
                          minor_breaks=NULL, expand=c(0,0), limits=c(91,320))
bin_season <- function(d) {
  d <- copy(d)
  d[, bi := findInterval(doy, hedges)]
  d <- d[bi>=1 & bi<=length(hedges)-1]
  d[, `:=`(tbin=(hedges[bi]+hedges[bi+1])/2, twidth=hedges[bi+1]-hedges[bi])]
  d
}
## 1906: mean of however many casts fall in a cell (usually 0 or 1, occasionally 2)
agg_1906 <- function(d) bin_season(d)[, .(m=mean(val)), by=.(tbin, twidth, depth)]
## modern: mean across years, blanked out if fewer than MIN_YEARS distinct years contributed
agg_modern <- function(d) {
  a <- bin_season(d)[, .(m=mean(val), n_yr=uniqueN(year)), by=.(tbin, twidth, depth)]
  a[n_yr < MIN_YEARS, m := NA][, n_yr := NULL][]
}

## ---- absolute-value palettes (DO ones shared with FigureSX_do_heatmaps.R) ----------------------
pal_temp <- scale_fill_distiller(palette="RdYlBu", direction=-1, limits=c(4,28), oob=squish,
                                 name="Temp (°C)", na.value="grey90")
do_cols <- c("#000000","#7f0000","#d7301f","#fdae61","#ffffbf","#66bd63","#1a9850")
pal_do <- scale_fill_gradientn(colours=do_cols, values=rescale(c(0,1,2,4,6,9,14)), name="DO (mg/L)",
                               limits=c(0,14), oob=squish, na.value="grey90")
pal_sat <- scale_fill_gradientn(colours=do_cols, values=rescale(c(0,10,20,40,70,95,140)), name="DO (% sat)",
                                limits=c(0,140), oob=squish, na.value="grey90")

## ---- difference palettes: only cells beyond a threshold are coloured, everything else very light grey ----
## Thresholds are Figure 1's +/- values: 1 deg C, 1.5 mg/L, 15 % sat. Colours are the darkest end of
## the extended Figure 1 scales (temp red/blue, DO green/brown).
THRESH <- c(temp=1, mgL=1.5, sat=15)
thresh_pal <- function(hi, lo, lab_hi, lab_lo, name) {
  scale_fill_manual(values=c(hi=hi, lo=lo, none="#efefef"), breaks=c("hi","lo"),
                    labels=c(lab_hi, lab_lo), name=name, drop=FALSE)
}
pal_tdiff   <- thresh_pal("#b2182b", "#2166ac", "> +1 \u00b0C",   "< \u22121 \u00b0C",   "\u0394 Temp")
pal_odiff   <- thresh_pal("#01665e", "#8c510a", "> +1.5 mg/L", "< \u22121.5 mg/L", "\u0394 DO")
pal_satdiff <- thresh_pal("#01665e", "#8c510a", "> +15 % sat", "< \u221215 % sat", "\u0394 DO")
classify <- function(a, thr) {
  a <- copy(a); a[, m := factor(fifelse(m > thr, "hi", fifelse(m < -thr, "lo", "none")),
                                levels=c("hi","lo","none"))]; a
}

panel <- function(a, depths, pal, title, show_y, tile_col="white") {
  a <- copy(a)
  a[, depth_lab := factor(sprintf("%g m", depth), levels=sprintf("%g m", rev(depths)))]
  ggplot(a, aes(tbin, depth_lab, fill=m)) +
    geom_tile(aes(width=twidth), height=0.92, color=tile_col, linewidth=0.12) +
    ## month starts drawn dark + dashed over the tiles, so they read differently from the
    ## white mid-month tile seams; labels sit mid-month between them
    geom_vline(xintercept=mbound, color="grey15", linewidth=0.2, linetype="22") +
    scale_y_discrete(drop=FALSE) + xsc + pal +
    labs(title=title, x=NULL, y=if (show_y) "Depth" else NULL) +
    theme_minimal(base_size=6) +
    theme(panel.grid=element_blank(), plot.title=element_text(face="bold", size=6.5),
          axis.text=element_text(size=5),
          ## compact colorbars so the four stacked legends fit beside the panels
          legend.key.height=unit(0.22,"cm"), legend.key.width=unit(0.2,"cm"),
          legend.title=element_text(size=5.5), legend.text=element_text(size=5),
          legend.spacing.y=unit(0.1,"cm"), plot.margin=margin(2,2,2,2))
}

## difference = 2016-2025 minus 1906, only in cells where both eras have a value
diff_agg <- function(a_old, a_new) {
  d <- merge(a_old, a_new, by=c("tbin","twidth","depth"), suffixes=c("_old","_new"))
  d[, .(tbin, twidth, depth, m = m_new - m_old)][!is.na(m)]
}

a_t06 <- agg_1906(t1906); a_tmd <- agg_modern(tmod); a_tdf <- diff_agg(a_t06, a_tmd)  # continuous difference, classified at plot time
DO <- list(
  mgL = list(a06=agg_1906(o1906), amd=agg_modern(omod), pal=pal_do,  dpal=pal_odiff, thr=THRESH[["mgL"]],
             lab="DO (mg/L)", file="figures/figSX_mendota_1906_comp.png"),
  sat = list(a06=agg_1906(s1906), amd=agg_modern(smod), pal=pal_sat, dpal=pal_satdiff, thr=THRESH[["sat"]],
             lab="DO (% sat)", file="figures/figSX_mendota_1906_comp_sat.png"))
cat(sprintf("diff ranges -- temp %.1f to %.1f; mg/L %.1f to %.1f; %%sat %.0f to %.0f\n",
            min(a_tdf$m), max(a_tdf$m),
            min(diff_agg(DO$mgL$a06, DO$mgL$amd)$m), max(diff_agg(DO$mgL$a06, DO$mgL$amd)$m),
            min(diff_agg(DO$sat$a06, DO$sat$amd)$m), max(diff_agg(DO$sat$a06, DO$sat$amd)$m)))

## rows = variable, columns = 1906 | 2016-2025 | difference
p_t06 <- panel(a_t06, DEPTHS$temp, pal_temp,  "Temperature — 1906", TRUE)
p_tmd <- panel(a_tmd, DEPTHS$temp, pal_temp,  "Temperature — 2016–2025", FALSE)
p_tdf <- panel(classify(a_tdf, THRESH[["temp"]]), DEPTHS$temp, pal_tdiff, "Temperature — difference", FALSE, tile_col="white")
fig_height <- 0.6 + 0.135 * sum(lengths(DEPTHS))   # 6.5 in wide, 500 dpi (journal single-page width)
for (v in DO) {
  p_o06 <- panel(v$a06, DEPTHS$do, v$pal,  paste(v$lab, "— 1906"), TRUE)
  p_omd <- panel(v$amd, DEPTHS$do, v$pal,  paste(v$lab, "— 2016–2025"), FALSE)
  p_odf <- panel(classify(diff_agg(v$a06, v$amd), v$thr), DEPTHS$do, v$dpal, paste(v$lab, "— difference"), FALSE, tile_col="white")
  ## each row collects its OWN legends (the absolute scale + its difference key) at its right
  ## edge, so the keys sit beside the panels they describe instead of in one detached stack
  row_t <- (p_t06 | p_tmd | p_tdf) + plot_layout(guides="collect")
  row_o <- (p_o06 | p_omd | p_odf) + plot_layout(guides="collect")
  g <- (row_t / row_o) + plot_layout(heights=c(length(DEPTHS$temp), length(DEPTHS$do))) +
    plot_annotation(title="Lake Mendota: 1906 vs. 2016\u20132025",
                    theme=theme(plot.title=element_text(face="bold", size=8)))
  ggsave(v$file, g, width=6.5, height=fig_height, dpi=500, bg="white")
  cat("wrote", v$file, "\n")
}

## Title + caption off the PNG, into the shared captions.csv (same convention as the other fig*.R)
write_captions <- function(new_caps) {
  path <- "figures/captions.csv"
  old <- if (file.exists(path)) fread(path) else data.table(file=character(), title=character(), caption=character())
  fwrite(rbind(old[!file %in% new_caps$file], new_caps), path)
  cat("wrote", path, "\n")
}
dlist <- function(x) paste0(paste(head(x, -1), collapse=", "), " and ", tail(x, 1), " m")
base_cap <- paste0(
  "by half-month period: 1906 (left), 2016-2025 (middle), and the difference 2016-2025 minus 1906 ",
  "(right; only where both eras have a value, and coloured only where the change exceeds a threshold: ",
  "temperature +/-1 deg C, DO +/-1.5 mg/L or +/-15 % sat -- Figure 1's values -- red/blue for warmer/cooler, ",
  "green/brown for more/less oxygen; all other cells very light grey). Note Figure 1 is a rate per decade, this is the total ",
  "change across ~110 years. Dashed vertical lines mark the start of each month; labels sit mid-month (the Nov label covers Nov 1-15 only). ",
  "The 1906 values are digitized from Birge and Juday (1911), Wisconsin Geological and Natural History ",
  "Survey Bulletin XXII, Plates I and II -- repeated casts at fixed depths (temperature: ", dlist(DEPTHS$temp),
  "; oxygen: ", dlist(DEPTHS$do), "). Plate II reports oxygen in cc/L, converted to mg/L (x 1.429). Each ",
  "modern panel is built at the same fixed depths as its 1906 counterpart, interpolating every 2016-2025 ",
  "cast to those depths rather than the finer continuous grid used elsewhere in this project, to avoid ",
  "implying a precision the historical record never had. A modern cell is blank if fewer than ", MIN_YEARS,
  " distinct years contributed data; a 1906 cell is blank if no cast fell in that half-month window. ",
  "Absolute DO colour scales shared with FigureSX_do_heatmaps.R.")
write_captions(data.table(
  file=c(DO$mgL$file, DO$sat$file),
  title=c("Lake Mendota temperature and dissolved oxygen (mg/L), 1906 vs. 2016-2025",
          "Lake Mendota temperature and dissolved oxygen (% saturation), 1906 vs. 2016-2025"),
  caption=c(
    paste0("Water temperature (top row) and dissolved oxygen (mg/L, bottom row) ", base_cap),
    paste0("Water temperature (top row) and dissolved oxygen (% saturation, bottom row) ", base_cap,
           " % saturation is computed identically for both eras from DO (mg/L) and water temperature: ",
           "Benson and Krause (1984) freshwater solubility with a barometric correction for Mendota's ",
           "elevation (258 m). For 1906, temperature at each oxygen cast comes from Plate I, interpolated in ",
           "time and depth; oxygen casts after Plate I ends (31 Oct) have no %sat. The modern NTL-LTER o2sat ",
           "column is not used: for Mendota it carries a constant +2.9% offset (a ~495 m elevation correction)."))
))
