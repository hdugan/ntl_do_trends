## FigureSX_metalimnionwindow.R: the metalimnion (thermocline band) as a ribbon through time.
##
## Uses the exact same per-cast classification as TableS2_ZoneTrends.R: a cast only has a
## resolved metalimnion if rLakeAnalyzer::meta.depths() returns BOTH a top and a bottom AND the
## profile is positively stratified (surface warmer than bottom) -- isothermal, unresolved, or
## inversely-stratified (cold-water, winter/early-spring) casts contribute no ribbon.
##
## Instead of a seasonal (day-of-year) climatology, the x-axis here is real calendar time across
## the full 1981-2025 record. Each year's run of resolved casts is drawn as its own ribbon
## polygon (ymin=metalimnion top, ymax=metalimnion bottom), so the figure reads as a sequence of
## pulses -- the thermocline band forming in spring, widening/deepening through summer, and
## collapsing at fall turnover -- repeated every year, with NO ribbon drawn across the intervening
## unstratified winter gap. Ribbons are colored by year (viridis, not the warm/cool diverging
## palette Figure1 uses for temperature, to avoid implying a warming/cooling read here) so
## multi-decade drift in the band's depth is visible directly as a color trend, not just position.
suppressMessages({library(data.table); library(ggplot2); library(rLakeAnalyzer)})

meta <- data.table(
  lakeid=c("TR","BM","CR","SP","AL","TB","CB","ME","MO","FI","WI"),
  name  =c("Trout Lake","Big Muskellunge Lake","Crystal Lake","Sparkling Lake","Allequash Lake",
           "Trout Bog","Crystal Bog","Lake Mendota","Lake Monona","Fish Lake","Lake Wingra"),
  region=c(rep("Northern",7), rep("Southern",4)),
  zmax  =c(35.7,21.3,20.4,20.0,8.0,7.9,2.5,25.3,22.5,18.9,4.0))

prof <- fread("data/profiles_clean.csv")[lakeid %in% meta$lakeid]
prof[, `:=`(date=as.Date(sampledate), doy=yday(as.Date(sampledate)), year=year4)]
prof <- prof[doy>=91 & doy<=319 & !(lakeid=="CR" & year %in% c(2012,2013))]  # same window as TableS2

## ---- per-cast thermal structure, identical method to TableS2_ZoneTrends.R -----------------------
strat <- prof[!is.na(wtemp), {
  o <- order(depth); dd <- depth[o]; ww <- wtemp[o]
  keep <- !duplicated(dd)
  dd <- dd[keep]; ww <- ww[keep]
  if (length(dd) < 3) {
    list(meta_top=NA_real_, meta_bot=NA_real_, inverse=NA)
  } else {
    md <- tryCatch(meta.depths(ww, dd), error=function(e) c(NA_real_, NA_real_))
    list(meta_top=md[1], meta_bot=md[2], inverse=ww[1] < ww[length(ww)])
  }
}, by=.(lakeid, date, year)]

## resolved AND positively stratified only -- the same eligibility rule TableS2 uses for
## Epilimnion/Metalimnion/Hypolimnion (not the fixed-zone override it applies for Crystal Bog/
## Wingra specifically -- this figure shows the raw dynamic signal for all 11 lakes, including
## those two, since the point here is to see the actual thermocline behavior, unstable or not)
strat[, valid := !is.na(meta_top) & !is.na(meta_bot) & !is.na(inverse) & !inverse]
ribbon_d <- merge(strat[valid==TRUE], meta, by="lakeid")

ord <- meta[order(region, -zmax), name]
ribbon_d[, strip := factor(name, levels=ord)]

MIN_PTS_PER_YEAR <- 3  # need at least this many resolved casts in a year to draw that year's pulse
ribbon_d[, n_this_year := .N, by=.(lakeid, year)]
ribbon_d <- ribbon_d[n_this_year >= MIN_PTS_PER_YEAR]
setorder(ribbon_d, lakeid, date)

g <- ggplot(ribbon_d, aes(x=date, group=interaction(lakeid, year))) +
  geom_ribbon(aes(ymin=meta_top, ymax=meta_bot, fill=year), color="grey15", linewidth=0.06, alpha=0.85) +
  scale_y_reverse() +
  scale_x_date(date_breaks="5 years", date_labels="%Y", expand=c(0.01,0.01)) +
  scale_fill_viridis_c(name="Year", option="viridis") +
  facet_wrap(~strip, ncol=4, scales="free_y") +
  labs(
    title="The metalimnion, 1981–2025: a ribbon through time",
    subtitle="Each pulse is one year's resolved, positively-stratified metalimnion band (rLakeAnalyzer::meta.depths); gaps are winter/spring-mixed/inversely-stratified periods with no resolved thermocline",
    x=NULL, y="Depth (m)") +
  theme_minimal(base_size=11) +
  theme(
    panel.grid.minor=element_blank(),
    strip.text=element_text(face="bold", size=9.5),
    plot.title=element_text(face="bold", size=15),
    plot.subtitle=element_text(size=9.5, color="grey35"),
    legend.position="right",
    panel.spacing=unit(0.7,"lines"))

ggsave("figures/figSX_metalimnion_window.png", g, width=13, height=9, dpi=500, bg="white")
cat("wrote figures/figSX_metalimnion_window.png\n")

## Title + long explanation off the PNG, into the shared captions.csv (same convention as the
## other fig*.R scripts)
write_captions <- function(new_caps) {
  path <- "figures/captions.csv"
  old <- if (file.exists(path)) fread(path) else data.table(file=character(), title=character(), caption=character())
  fwrite(rbind(old[!file %in% new_caps$file], new_caps), path)
  cat("wrote", path, "\n")
}
write_captions(data.table(
  file="figures/figSX_metalimnion_window.png",
  title="The metalimnion as a ribbon through time",
  caption=paste0(
    "Metalimnion top and bottom depth for each cast, April 1--November 15, full period of record ",
    "(Crystal Lake 2012-2013 excluded; whole-lake mixing experiment), computed with rLakeAnalyzer::",
    "meta.depths() -- the same method and eligibility rule as TableS2_ZoneTrends.R: a cast only ",
    "contributes if the metalimnion top AND bottom are both resolved and the profile is positively ",
    "stratified (surface warmer than bottom); isothermal, unresolved, or inversely-stratified ",
    "(cold-water, winter/early-spring) casts leave no ribbon. Each year's run of eligible casts is ",
    "one ribbon polygon (a 'pulse'), colored by year (viridis), so long-term drift in the band's ",
    "depth appears as a positional AND color trend across the record, not just within one season. ",
    "Shown for all 11 lakes including Crystal Bog and Wingra, which TableS2_ZoneTrends.R instead ",
    "assigns a fixed shallow-depth Epilimnion/Hypolimnion split for trend testing (too shallow for ",
    "stable thermocline detection, and polymictic, respectively) -- their ribbons here are the raw, ",
    "often sparse or unstable signal that motivated that override."
  )
))
