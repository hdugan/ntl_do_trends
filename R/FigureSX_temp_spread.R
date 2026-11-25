## FigureSX_temp_spread.R: how much surface and bottom-water temperature VARY by month and lake,
## ignoring their absolute level. One figure, two panels sharing one colour scale:
##   (A) surface: mean of 0-1 m readings per profile
##   (B) bottom:  mean of readings in the deepest 1 m of each profile; a profile only counts if it reaches
##       at least 90% of that lake's median profile depth, so a cast that stopped short of the bottom is
##       not mistaken for a bottom reading
## Spread = 5th-95th percentile range (degC) across all profiles in that lake x month, full record.
## April-November; Crystal Lake 2012-2013 excluded (whole-lake mixing experiment).
## Lakes in the same order as Figures 2 and 6 (northern no-clarity-change, northern browning, southern;
## descending max depth within a group).
## (Replaces FigureSX_surface_temp_spread.R and FigureSX_bottom_temp.R.)
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})

meta <- data.table(
  lakeid=c("TR","BM","SP","AL","CR","TB","CB","ME","MO"),
  name  =c("Trout Lake","Big Muskellunge Lake","Sparkling Lake","Allequash Lake","Crystal Lake",
           "Trout Bog","Crystal Bog","Lake Mendota","Lake Monona"))   # Figure 2 order

prof <- fread("data/profiles_clean.csv")[lakeid %in% meta$lakeid & !is.na(wtemp)]
prof[, `:=`(date=as.Date(sampledate), year=year4)]
prof <- prof[!(lakeid=="CR" & year %in% c(2012,2013))]

## surface: mean of 0-1 m per cast
sfc <- prof[depth <= 1, .(temp=mean(wtemp)), by=.(lakeid, date)]
## bottom: mean of the deepest 1 m per cast, casts reaching >= 90% of the lake's median cast depth
btm <- prof[, .(zbot=max(depth), temp=mean(wtemp[depth >= max(depth) - 1])), by=.(lakeid, date)]
btm[, zref := median(zbot), by=lakeid]
btm <- btm[zbot >= 0.9 * zref, .(lakeid, date, temp)]

spread <- function(d){
  d <- d[month(date) %in% 4:11]
  sp <- d[, .(spread=unname(diff(quantile(temp, c(0.05, 0.95)))), n=.N), by=.(lakeid, month=month(date))]
  sp <- merge(sp, meta, by="lakeid")
  sp[, mlab := factor(month.abb[month], levels=month.abb[4:11])]
  sp[, lake := factor(name, levels=rev(meta$name))]
  sp
}
S <- spread(sfc); B <- spread(btm)
LIM <- range(c(S$spread, B$spread))
thr <- mean(LIM)                                   # white labels on the darker half of the scale

panel <- function(sp, ttl, top_months){
  ggplot(sp, aes(mlab, lake, fill=spread)) +
    geom_tile(color="white", linewidth=0.6) +
    geom_text(aes(label=sprintf("%.1f", spread), color=spread > thr), size=2.4) +
    scale_color_manual(values=c(`TRUE`="white", `FALSE`="grey15"), guide="none") +
    scale_fill_viridis_c(option="mako", direction=-1, limits=LIM, name="Spread (°C, 5th–95th percentile)") +
    scale_x_discrete(position="top", expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
    labs(x=NULL, y=NULL, title=ttl) +
    theme_minimal(base_size=8.5) +
    theme(panel.grid=element_blank(), axis.text.x=element_text(size=7.5),
          plot.title=element_text(face="bold", size=9, margin=margin(b=2)), plot.title.position="plot",
          legend.position="bottom", legend.key.width=unit(1,"cm"), legend.key.height=unit(0.25,"cm"),
          legend.title.position="top", legend.title=element_text(hjust=0.5, size=8)) +
    { if(!top_months) theme(axis.text.x=element_blank()) }
}
g <- (panel(S, "A  Surface water (0–1 m)", TRUE) / panel(B, "B  Bottom water (deepest 1 m)", FALSE)) +
  plot_layout(guides="collect") & theme(legend.position="bottom")
ggsave("figures/figSX_temp_spread.png", g, width=6.5, height=6.6, dpi=500, bg="white")
cat("wrote figures/figSX_temp_spread.png\n")
print(merge(S[, .(surface=round(mean(spread),1)), by=mlab], B[, .(bottom=round(mean(spread),1)), by=mlab])[order(mlab)])

write_captions <- function(new_caps) {
  path <- "figures/captions.csv"
  old <- if (file.exists(path)) fread(path) else data.table(file=character(), title=character(), caption=character())
  old <- old[!file %in% c("figures/figSX_surface_temp_spread.png", "figures/figSX_bottom_temp_spread.png")]  # replaced
  fwrite(rbind(old[!file %in% new_caps$file], new_caps), path)
  cat("wrote", path, "\n")
}
write_captions(data.table(
  file="figures/figSX_temp_spread.png",
  title="Spread of surface and bottom-water temperature by month and lake",
  caption=paste0(
    "Range of (A) surface water temperature (mean of 0-1 m readings per profile) and (B) bottom-water ",
    "temperature (mean of readings in the deepest 1 m of each profile; profiles that did not reach at least ",
    "90% of the lake's median profile depth were excluded) across all profiles in each lake and month, ",
    "April-November, full period of record. Spread is the difference between the 95th and 5th percentiles ",
    "(°C) and is independent of the absolute temperature level; darker cells vary more. Both panels share ",
    "one colour scale. Lakes ordered as in Figure 2. Crystal Lake 2012-2013 excluded (whole-lake mixing ",
    "experiment).")))
