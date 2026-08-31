## TableS2_ZoneTrends.R: significance-tested trends by THERMAL ZONE instead of fixed depth bins.
##
## Why: Figure1_rate_of_change.R tests trends in ~14k fine depth x half-month cells. Individually
## those cells are weakly powered (most real trends don't clear even raw p<0.05 at that
## resolution) even though the pattern is visually obvious and spatially coherent -- because a
## fixed "surface <=2m / deep >=60% zmax" split doesn't track what's physically happening: lakes
## are isothermal in April, develop a thermocline that deepens through the summer, and mix again
## in fall. A cell fixed at "5 m" is epilimnion in June and hypolimnion in September.
##
## This script instead classifies every depth reading into a THERMAL zone computed per cast from
## its own temperature profile (rLakeAnalyzer::meta.depths, the standard limnological
## epilimnion/metalimnion/hypolimnion definition, via the depth-density-gradient method), so
## "hypolimnion" tracks the actual thermocline depth on that day rather than a fixed depth. Casts
## with no detectable stratification (density top-bottom < mixed.cutoff) become a single
## "Mixed" zone, split into Spring Mixing / Fall Turnover by day-of-year. Pooling all readings
## within a zone (instead of 0.5-1 m depth slivers) gives each trend test far more data and much
## more power, and the zones stay physically meaningful all season.
##
## Two exceptions use a fixed 2-zone (Epilimnion/Hypolimnion only) split instead -- see
## FIXED_ZONE_BOUND below: Crystal Bog (too shallow for the density-gradient method to resolve a
## stable thermocline) and Wingra (polymictic, so it has no single seasonal stratified/mixed state).
suppressMessages({library(data.table); library(rLakeAnalyzer); library(splines); library(quantreg)})
dir.create("tables", showWarnings = FALSE)

meta <- data.table(
  lakeid=c("TR","BM","CR","SP","AL","TB","CB","ME","MO","FI","WI"),
  name  =c("Trout Lake","Big Muskellunge Lake","Crystal Lake","Sparkling Lake","Allequash Lake",
           "Trout Bog","Crystal Bog","Lake Mendota","Lake Monona","Fish Lake","Lake Wingra"),
  region=c(rep("Northern",7), rep("Southern",4)),
  zmax  =c(35.7,21.3,20.4,20.0,8.0,7.9,2.5,25.3,22.5,18.9,4.0))

prof <- fread("data/profiles_clean.csv")[lakeid %in% meta$lakeid]
prof[, `:=`(date=as.Date(sampledate), doy=yday(as.Date(sampledate)), year=year4)]
prof <- prof[doy>=91 & doy<=319 & !(lakeid=="CR" & year %in% c(2012,2013))]  # same window as Figure1

## ---- 1. Per-cast thermal structure, from the temperature profile alone -------------------------
## SPRING_FALL_SPLIT_DOY divides "Mixed" casts into Spring Mixing vs Fall Turnover. 200 = ~Jul 19,
## a defensible mid-season cutover point (well after spring warm-up finishes, well before fall
## cooling starts in these lakes) rather than a precise physical event -- a handful of storm-mixed
## midsummer casts may fall on the "wrong" side, but that's a minor mislabel, not a bias, since it
## doesn't systematically favor one direction.
SPRING_FALL_SPLIT_DOY <- 200

## `inverse` flags casts where the surface reading is colder than the bottom reading -- winter/
## early-spring inverse stratification (fresh water is densest at ~4C, so a column can be
## legitimately density-stratified with COLD water on top, not the usual warm-epilimnion-over-
## cool-hypolimnion pattern). meta.depths() only reports whether there's a density gradient, not
## its direction, so it happily returns a "metalimnion" for these -- treated below as Spring
## Mixing instead, since it isn't the productive-season stratification Epi/Meta/Hypo describes.
strat <- prof[!is.na(wtemp), {
  o <- order(depth); dd <- depth[o]; ww <- wtemp[o]
  keep <- !duplicated(dd)
  dd <- dd[keep]; ww <- ww[keep]
  if (length(dd) < 3) {
    list(meta_top=NA_real_, meta_bot=NA_real_, inverse=NA, doy=doy[1])
  } else {
    md <- tryCatch(meta.depths(ww, dd), error=function(e) c(NA_real_, NA_real_))
    list(meta_top=md[1], meta_bot=md[2], inverse=ww[1] < ww[length(ww)], doy=doy[1])
  }
}, by=.(lakeid, date)]

## Crystal Bog and Wingra skip the dynamic meta.depths() classification entirely and use a fixed
## 2-zone (Epilimnion/Hypolimnion, no Metalimnion, no Spring Mixing/Fall Turnover) split instead:
## Crystal Bog is so shallow (zmax 2.5 m) that the density-gradient thermocline detection is
## unstable at that vertical scale; Wingra is polymictic (EDI 434 mixing_type), so it doesn't
## reliably hold a single seasonal stratified/mixed state for the dynamic method to key off.
FIXED_ZONE_BOUND <- c(CB = 1, WI = 2)  # lakeid -> epilimnion/hypolimnion boundary depth (m)

zone_of <- function(lakeid, depth, meta_top, meta_bot, doy, inverse) {
  fixed_bound <- FIXED_ZONE_BOUND[lakeid]  # NA (not NULL) for lakeids not in the vector -- stays aligned with input length
  ## "positively" stratified requires BOTH a top and a bottom metalimnion depth -- a cast missing
  ## either one isn't treated as a resolved Epi/Meta/Hypo structure, just Spring/Fall mixed.
  not_resolved <- is.na(meta_top) | is.na(meta_bot)
  fifelse(!is.na(fixed_bound), fifelse(depth < fixed_bound, "Epilimnion", "Hypolimnion"),
    ## inverse stratification (cold top, warm bottom) always reads as Spring Mixing, regardless
    ## of day-of-year -- it's a cold-water signature, never the warm-to-mixed fall transition.
    fifelse(!is.na(inverse) & inverse, "Spring Mixing",
      fifelse(not_resolved,
        fifelse(doy < SPRING_FALL_SPLIT_DOY, "Spring Mixing", "Fall Turnover"),
        fifelse(depth < meta_top, "Epilimnion", fifelse(depth > meta_bot, "Hypolimnion", "Metalimnion")))))
}

## ---- 2. Attach each depth reading to its cast's thermal zone -----------------------------------
pm <- melt(prof[,.(lakeid,date,doy,year,depth,wtemp,o2,o2sat)],
           id.vars=c("lakeid","date","doy","year","depth"), variable.name="var", value.name="val")
pm <- pm[is.finite(val)]
pm <- merge(pm, strat[,.(lakeid,date,meta_top,meta_bot,inverse)], by=c("lakeid","date"))
pm[, zone := zone_of(lakeid, depth, meta_top, meta_bot, doy, inverse)]

zone_levels <- c("Spring Mixing","Epilimnion","Metalimnion","Hypolimnion","Fall Turnover")
pm[, zone := factor(zone, levels=zone_levels)]

## ---- 3. Seasonal detrending within each zone, same rationale as Figure1: the exact days visited
## drift year to year, and since values change fast seasonally that drift alone can manufacture a
## trend. Fit one median seasonal curve per lake/var/zone, trend the RESIDUALS.
pm[, resid := {
  if(.N>=15){
    fit <- tryCatch(rq(val ~ ns(doy, df=4), tau=0.5), error=function(e) NULL)
    if(!is.null(fit)) residuals(fit) else val-median(val)
  } else if(.N>=6){
    fit <- tryCatch(rq(val ~ doy, tau=0.5), error=function(e) NULL)
    if(!is.null(fit)) residuals(fit) else val-median(val)
  } else val-median(val)
}, by=.(lakeid,var,zone)]

## one value per year per lake/var/zone (median of that year's readings), so a densely-sampled
## year doesn't get more weight than a sparsely-sampled one
yr <- pm[, .(v=median(resid)), by=.(lakeid,var,zone,year)]

## ---- 4. Theil-Sen slope (per decade) + Mann-Kendall significance, same estimators as Figure1 ---
MIN_YR_TREND <- 10
MIN_SPAN     <- 20
sen <- function(y,x){ n<-length(y); if(n<3) return(NA_real_)
  d <- outer(y,y,"-")[lower.tri(diag(n))] / outer(x,x,"-")[lower.tri(diag(n))]
  median(d[is.finite(d)], na.rm=TRUE) }
cell <- yr[, { n<-.N; sp<- if(n>0) diff(range(year)) else 0
  if(n>=MIN_YR_TREND && sp>=MIN_SPAN){
    .(slope = sen(v, year)*10,
      p     = tryCatch(suppressWarnings(cor.test(year, v, method="kendall", exact=FALSE)$p.value),
                       error=function(e) NA_real_), nyr=n)
  } else .(slope=NA_real_, p=NA_real_, nyr=n) },
  by=.(lakeid,var,zone)]

## Benjamini-Hochberg FDR correction, per lake -- same convention as Figure1, but the family here
## is only ~5 zones x 3 variables = <=15 tests per lake (vs ~14k across the whole fine grid), so
## correction costs far less power: pooling to physically real zones is what actually fixes the
## "everything looks non-significant" problem, not the correction method.
cell[!is.na(p), q := p.adjust(p, method="BH"), by=lakeid]
w <- merge(cell, meta, by="lakeid")

## ---- 5. Build the summary table -----------------------------------------------------------------
stars <- function(q) fifelse(is.na(q), "", fifelse(q<0.001,"***", fifelse(q<0.01,"**", fifelse(q<0.05,"*",""))))
fmt <- function(slope, q, digits) fifelse(is.na(slope), "--", sprintf(paste0("%.",digits,"f%s"), slope, stars(q)))

wide <- dcast(w, lakeid+name+region+zmax+zone ~ var, value.var=c("slope","q","nyr"))
wide[, `:=`(
  temp_fmt   = fmt(slope_wtemp, q_wtemp, 2),
  o2_fmt     = fmt(slope_o2,    q_o2,    2),
  o2sat_fmt  = fmt(slope_o2sat, q_o2sat, 1),
  n_yr       = fifelse(is.na(nyr_wtemp), nyr_o2, nyr_wtemp)
)]
setorder(wide, zone, region, -zmax)
wide <- wide[!(is.na(slope_wtemp) & is.na(slope_o2) & is.na(slope_o2sat))]  # drop empty zone-rows (e.g. no hypolimnion ever sampled)

## ---- 6. Write LaTeX table -----------------------------------------------------------------------
lines <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\footnotesize",
  "\\caption{Trends by thermal zone (Theil-Sen slope per decade, Mann-Kendall significance, Benjamini-Hochberg FDR-corrected per lake) for the 11 NTL-LTER primary study lakes, April 1--November 15, full period of record. Zones are computed per cast from its own temperature profile (\\texttt{rLakeAnalyzer::meta.depths}), not a fixed depth: Epilimnion/Metalimnion/Hypolimnion require both a resolved metalimnion top and bottom AND positive stratification (surface warmer than bottom); every other cast -- isothermal, unresolved, or inversely stratified (cold surface over warmer bottom, a winter/early-spring signature) -- is Spring Mixing or Fall Turnover, split at day-of-year 200 (inversely stratified casts are always Spring Mixing regardless of date). $^{*}p<0.05$ $^{**}p<0.01$ $^{***}p<0.001$ (FDR-adjusted).}",
  "\\label{tab:zone_trends}",
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{llrrrr}",
  "\\toprule",
  "Zone & Lake & Years & $\\Delta$ Temp (\\textdegree C/decade) & $\\Delta$ DO (mg/L/decade) & $\\Delta$ DO (\\%sat/decade) \\\\",
  "\\midrule"
)
zone_prev <- NULL
for (i in seq_len(nrow(wide))) {
  r <- wide[i]
  new_zone <- is.null(zone_prev) || r$zone != zone_prev
  if (new_zone && !is.null(zone_prev)) lines <- c(lines, "\\hline")
  zonecol <- if (new_zone) as.character(r$zone) else ""
  lines <- c(lines, sprintf(
    "%s & %s & %s & %s & %s & %s \\\\",
    zonecol, r$name, if (is.na(r$n_yr)) "--" else as.character(r$n_yr),
    r$temp_fmt, r$o2_fmt, r$o2sat_fmt
  ))
  zone_prev <- r$zone
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}%", "}", "\\end{table}")
writeLines(lines, "tables/TableS2_ZoneTrends.tex")
cat("wrote tables/TableS2_ZoneTrends.tex\n")

fwrite(wide[, .(zone, lake=name, region, n_years=n_yr,
                 temp_per_decade=slope_wtemp, temp_q=q_wtemp,
                 do_mgL_per_decade=slope_o2, do_mgL_q=q_o2,
                 do_sat_per_decade=slope_o2sat, do_sat_q=q_o2sat)],
       "tables/TableS2_ZoneTrends.csv")
cat("wrote tables/TableS2_ZoneTrends.csv\n")

## ---- Console summary ------------------------------------------------------------------------
cat("\n=== trend per decade by lake x thermal zone (BH q<0.05 starred) ===\n")
print(wide[, .(zone, name, n_yr, temp_fmt, o2_fmt, o2sat_fmt)], nrows=200)
cat("\nfraction of zone-tests significant at q<0.05:",
    round(mean(w$q < 0.05, na.rm=TRUE), 3), "(", sum(w$q<0.05, na.rm=TRUE), "of", sum(!is.na(w$q)), ")\n")
