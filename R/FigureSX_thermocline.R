## FigureSD_thermocline.R: thermocline depth and center of buoyancy through the open-water season,
## one panel per lake. Writes two figures: figSD_thermocline.png (rLakeAnalyzer::thermo.depth, depth of
## the steepest density gradient) and figSD_centerbuoyancy.png (rLakeAnalyzer::center.buoyancy, the
## buoyancy-frequency-weighted mean depth -- smoother, ~20-40% less cast-to-cast noise in most lakes).
## Both metrics are computed per cast using the same eligibility rule as
## TableS2_ZoneTrends.R / FigureSX_metalimnionwindow.R: the cast must have a resolved metalimnion
## (meta.depths top AND bottom) and be positively stratified (surface warmer than bottom); other
## casts have no thermocline. April 1-October 31, Crystal Lake 2012-2013 excluded.
## Each year is one thin line (day of year on x); colour = the four decade groups and palette of
## Figure3_o2_profiles.R (1985-1995 blue, 1995-2005 green, 2005-2015 orange, 2015-2025 red). As in
## Figure 3, years before 1985 (northern lakes 1981-1984) fall outside the groups and are not drawn.
## Thick lines = each group's median thermocline depth by half-month period (Figure 1's bin edges).
suppressMessages({library(data.table); library(ggplot2); library(rLakeAnalyzer)})

meta <- data.table(
  lakeid=c("TR","BM","CR","SP","AL","TB","CB","ME","MO"),
  name  =c("Trout Lake","Big Muskellunge Lake","Crystal Lake","Sparkling Lake","Allequash Lake",
           "Trout Bog","Crystal Bog","Lake Mendota","Lake Monona"),
  region=c(rep("Northern",7), rep("Southern",2)),
  zmax  =c(35.7,21.3,20.4,20.0,8.0,7.9,2.5,25.3,22.5))
## Crystal Bog (2.5 m, sampled every 0.5 m) is dropped: its thermocline can only fall at 0.5, 1.0 or
## 1.5 m, so the lines just jump between those values (Table S2 uses a fixed zone split for it instead).
meta <- meta[lakeid != "CB"]

prof <- fread("data/profiles_clean.csv")[lakeid %in% meta$lakeid & !is.na(wtemp)]
prof[, `:=`(date=as.Date(sampledate), doy=yday(as.Date(sampledate)), year=year4)]
prof <- prof[doy >= 91 & doy <= 304 & !(lakeid=="CR" & year %in% c(2012,2013))]

th <- prof[, {
  o <- order(depth); dd <- depth[o]; ww <- wtemp[o]
  keep <- !duplicated(dd); dd <- dd[keep]; ww <- ww[keep]
  if (length(dd) < 3) list(thermo=NA_real_, cb=NA_real_) else {
    md <- tryCatch(meta.depths(ww, dd), error=function(e) c(NA_real_, NA_real_))
    ok <- !any(is.na(md)) && ww[1] > ww[length(ww)]
    if (!ok) list(thermo=NA_real_, cb=NA_real_) else
      list(thermo = tryCatch(thermo.depth(ww, dd), error=function(e) NA_real_),
           cb     = tryCatch(center.buoyancy(ww, dd), error=function(e) NA_real_))
  }
}, by=.(lakeid, date, doy, year)]

breaks <- c(1985,1995,2005,2015,2026)
lv <- c("1985-1995","1995-2005","2005-2015","2015-2025")
th[, era := { idx <- findInterval(year, breaks)
              fifelse(idx >= 1 & idx <= 4, lv[pmin(pmax(idx,1L),4L)], NA_character_) }]
th <- th[!is.na(era)]
th[, era := factor(era, levels=lv)]
pal <- setNames(c("#2166ac","#1a9850","#e08214","#b2182b"), lv)

th <- merge(th, meta, by="lakeid")
th[, strip := factor(name, levels=meta[order(region, -zmax), name])]
setorder(th, lakeid, year, doy)

metrics <- list(
  list(col="thermo", file="figures/figSD_thermocline.png",    ylab="Thermocline depth (m)"),
  list(col="cb",     file="figures/figSD_centerbuoyancy.png", ylab="Center of buoyancy (m)"))
hedges <- c(91,106,121,136,152,167,182,197,213,228,244,259,274,289,305)
mbound <- c(91,121,152,182,213,244,274,305)
th[, bi := findInterval(doy, hedges)]
for (m in metrics) {
  d <- th[is.finite(get(m$col))][, val := get(m$col)]
  ## era medians by half-month bin (Figure 1's edges), drawn at bin centres
  med <- d[bi >= 1 & bi < length(hedges), .(val=median(val), n=.N), by=.(strip, era, bi)][n >= 3]
  med[, doy := (hedges[bi] + hedges[bi + 1]) / 2]
  setorder(med, strip, era, doy)
  g <- ggplot(d, aes(doy, val, color=era)) +
    geom_vline(xintercept=mbound, color="grey92", linewidth=0.3) +
    geom_line(aes(group=year), linewidth=0.2, alpha=0.35) +
    geom_line(data=med, linewidth=1) +
    scale_y_reverse() +
    scale_x_continuous(breaks=c(106,136.5,167,197.5,228.5,259,289.5),
                       labels=c("Apr","May","Jun","Jul","Aug","Sep","Oct"),
                       limits=c(91, 305), expand=c(0,0)) +
    scale_color_manual(values=pal, name=NULL) +
    facet_wrap(~strip, ncol=3, scales="free_y") +
    labs(x=NULL, y=m$ylab) +
    theme_minimal(base_size=10.5) +
    theme(panel.grid=element_blank(), strip.text=element_text(face="bold", size=9.5),
          axis.text.x=element_text(size=7.5), legend.position="bottom",
          panel.spacing=unit(0.7,"lines")) +
    guides(color=guide_legend(override.aes=list(linewidth=1.2, alpha=1)))
  ggsave(m$file, g, width=9.5, height=7.5, dpi=500, bg="white")
  cat("wrote", m$file, "\n")
}

write_captions <- function(new_caps) {
  path <- "figures/captions.csv"
  old <- if (file.exists(path)) fread(path) else data.table(file=character(), title=character(), caption=character())
  fwrite(rbind(old[!file %in% new_caps$file], new_caps), path)
  cat("wrote", path, "\n")
}
common <- paste0(
  "A cast contributes only if its metalimnion top and bottom are both resolved (meta.depths) and it is ",
  "positively stratified (surface warmer than bottom) -- the same rule as Table S2. Thin lines are individual ",
  "years; thick lines are each decade group's median by half-month period (shown where at least 3 casts ",
  "contribute). Colours follow Figure 3: 1985-1995 blue, 1995-2005 green, 2005-2015 orange, 2015-2025 red; ",
  "years before 1985 are not shown. Panels ordered north then south by max depth (zmax); each panel has its ",
  "own depth axis. Crystal Bog is not shown: at 2.5 m deep and 0.5 m sampling, its thermocline can only take ",
  "three values. Crystal Lake 2012-2013 excluded (whole-lake mixing experiment).")
write_captions(data.table(
  file=c("figures/figSD_thermocline.png", "figures/figSD_centerbuoyancy.png"),
  title=c("Seasonal thermocline depth by decade", "Seasonal center of buoyancy by decade"),
  caption=c(
    paste("Thermocline depth (rLakeAnalyzer::thermo.depth; depth of the steepest density gradient) for each cast, April 1-October 31.", common),
    paste("Center of buoyancy (rLakeAnalyzer::center.buoyancy; the depth of the buoyancy-frequency-weighted centroid of the density gradient) for each cast, April 1-October 31. It integrates the whole profile, so it is less sensitive than thermocline depth to the sampling interval and to competing gradients.", common))
))
