## Fig 6: which parts of the open-water season are changing most? A seasonal summary of Figure 2
## (Figure1_rate_of_change.R), all 9 NTL lakes (Fish Lake and Wingra excluded).
## Reads figures/fig01_cell_trends.csv and figures/fig01_block_trends.csv -- run Figure1_rate_of_change.R first.
##
## For each lake x half-month period, collapsing over depth:
##   tile fill = mean Theil-Sen slope across that lake's depth bins (signed; degC or mg/L per decade)
##   dot size  = % of that lake's depth bins lying in a SIGNIFICANT DEPTH BLOCK (the block test of the
##               fig01_rate_block_*.png versions: BH q<0.05 per lake), so Figures 2 and 6 agree
## Top bars: that same % in each half-month, averaged across lakes so each lake counts
## equally (otherwise Trout, with the most depth bins, would dominate).
## Panel (a) weights every depth bin equally; panel (b) weights each depth bin by the lake AREA at that
## depth x bin thickness, i.e. the share of lake VOLUME it represents (data/Bathymetry/lakehypsometry.csv), so
## tiles become volume-weighted mean trends and dots the % of lake VOLUME in significant blocks. Shallow,
## broad layers then count for more than the narrow deep hole.
## Caveat: averaging over depth lets opposite trends cancel (e.g. Crystal's surface warming and deep
## cooling in Sep-Oct give a pale tile) -- the dot still shows that the period is active.
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(scales)})

path <- "figures/fig01_cell_trends.csv"
if(!file.exists(path)) stop(path, " not found -- run Figure1_rate_of_change.R first")
w <- fread(path)[var %in% c("wtemp","o2")]
bpath <- "figures/fig01_block_trends.csv"
if(!file.exists(bpath)) stop(bpath, " not found -- run Figure1_rate_of_change.R first")
bk <- fread(bpath)[var %in% c("wtemp","o2")]
## each 1-m (0.5-m) cell inherits the significance of the depth block it sits in
w[bk, on=.(lakeid, var, tbin, dbin >= depth_top, dbin <= depth_bottom), bsig := i.sig]
w[is.na(bsig), bsig := FALSE]

ord <- unique(w[, .(name, region, zmax)])[order(region, -zmax)]$name   # north first, then by zmax
w[, name := factor(name, levels=rev(ord))]
NSOUTH <- uniqueN(w[region=="Southern", name])   # southern lakes sit at the bottom of the reversed axis

## volume weights from data/Bathymetry/lakehypsometry.csv: hp_factor = fraction of lake volume in a layer
## (each lake sums to 1). Two layer conventions in that file: northern lakes -- `depth` is the layer's BOTTOM
## (0.5 or 1 m layers); Mendota and Monona -- `depth` is the layer's CENTRE (0.25 = 0-0.5 m, then d-0.5 to d+0.5).
## Cumulative volume V(z) at layer bottoms is interpolated, and a cell centred on dbin with thickness dstep
## gets weight V(dbin + dstep/2) - V(dbin - dstep/2).
hyp <- fread("data/Bathymetry/lakehypsometry.csv")
setorder(hyp, lakeid, depth)
hyp[, bot := if(.BY$lakeid %in% c("ME","MO")) fifelse(depth == min(depth), 2*depth, depth + 0.5) else depth, by=lakeid]
cumv <- hyp[, .(z=c(0, bot), V=c(0, cumsum(hp_factor))/sum(hp_factor)), by=lakeid][, .SD[!duplicated(z)], by=lakeid]
w[, dstep := fifelse(zmax < 10, 0.5, 1)]
w[, wt := { cv <- cumv[lakeid == .BY$lakeid]
  pmax(0, approx(cv$z, cv$V, dbin + dstep/2, rule=2)$y - approx(cv$z, cv$V, pmax(0, dbin - dstep/2), rule=2)$y) }, by=lakeid]

summ <- function(weighted){
  x <- copy(w)[, ww := if(weighted) wt else 1]
  s <- x[, .(net=sum(slope*ww)/sum(ww), pct=pmin(100, 100*sum(bsig*ww)/sum(ww))), by=.(var, name, tbin, twidth)]   # pmin: float rounding can give 100.0000001
  list(s=s, top=s[, .(pct=mean(pct)), by=.(var, tbin, twidth)])
}
A <- summ(FALSE); B <- summ(TRUE)
for(z in list(list("a: each depth bin equal", A), list("b: weighted by lake area at depth", B))){
  cat("\n", z[[1]], "-- % significant by half-month (mean across lakes):\n")
  print(dcast(z[[2]]$top[, .(var, tbin, pct=round(pct))], var ~ tbin, value.var="pct")) }

## month boundaries + tick centres, same convention as Figure1_rate_of_change.R
mbound <- c(91,121,152,182,213,244,274,305,320)
xsc <- scale_x_continuous(breaks=c(106,136.5,167,197.5,228.5,259,289.5,312.5),
         labels=c("Apr","May","Jun","Jul","Aug","Sep","Oct","Nov"), limits=c(91,320), expand=c(0,0))
th <- theme_minimal(base_size=7) + theme(
  panel.grid=element_blank(), plot.title=element_text(face="bold", size=7.5, hjust=0.5),
  legend.position="bottom", legend.key.height=unit(0.22,"cm"), legend.key.width=unit(0.5,"cm"),
  legend.title=element_text(size=6), legend.text=element_text(size=5.5), axis.text=element_text(size=6))

## same palettes and limits as Figure 2 so the colours mean the same thing in both figures
pal_T <- scale_fill_gradientn(colours=c("#2166ac","#4393c3","#92c5de","#d1e5f0","#f7f7f7","#fddbc7","#f4a582","#d6604d","#b2182b"),
           limits=c(-1,1), oob=squish, breaks=c(-1,0,1), name="Mean Δ temp\n(°C/decade)")
pal_O <- scale_fill_gradientn(colours=c("#8c510a","#bf812d","#dfc27d","#f6e8c3","#f5f5f5","#c7eae5","#80cdc1","#35978f","#01665e"),
           limits=c(-1.5,1.5), oob=squish, breaks=c(-1.5,0,1.5), name="Mean Δ DO\n(mg/L/decade)")
YMAX <- ceiling(max(c(A$top$pct, B$top$pct))/10)*10

column <- function(vv, ttl, pal, left, S, tag=NULL){
  s <- S$s; top <- S$top
  bars <- ggplot(top[var==vv], aes(tbin, pct)) +
    geom_col(aes(width=twidth*0.9), fill="grey40") +
    geom_vline(xintercept=mbound, color="white", linewidth=0.4) +
    xsc + scale_y_continuous(limits=c(0,YMAX), expand=c(0,0)) +
    labs(title=ttl, tag=tag, x=NULL, y=if(left) "% significant" else NULL) + th +
    theme(plot.tag=element_text(face="bold", size=9)) +
    theme(axis.text.x=element_blank(), legend.position="none") +
    { if(!left) theme(axis.text.y=element_blank()) }
  tiles <- ggplot(s[var==vv], aes(tbin, name)) +
    geom_tile(aes(width=twidth, fill=net), height=0.95) +
    geom_point(data=s[var==vv & pct>0], aes(size=pct), shape=21, fill="white", color="grey20", stroke=0.25) +
    geom_vline(xintercept=mbound, color="white", linewidth=0.4) +
    geom_hline(yintercept=NSOUTH + 0.5, linetype="dashed", color="grey20", linewidth=0.35) +   # north | south
    xsc + pal +
    scale_size_area(max_size=3, limits=c(0,100), breaks=c(10,30,60), name="% of depth (a) or volume (b)\nin significant blocks") +
    guides(fill=guide_colourbar(order=1, title.position="top"), size=guide_legend(order=2, title.position="top")) +
    labs(x=NULL, y=NULL) + th +
    { if(!left) theme(axis.text.y=element_blank()) }
  bars / tiles + plot_layout(heights=c(1,4))
}
pa <- column("wtemp", "Temperature", pal_T, TRUE, A, "a") | column("o2", "Dissolved oxygen", pal_O, FALSE, A)
pb <- column("wtemp", "Temperature, weighted by lake area at depth", pal_T, TRUE, B, "b") |
      column("o2", "Dissolved oxygen, weighted by lake area at depth", pal_O, FALSE, B)
p <- (pa / pb) + plot_layout(guides="collect") & theme(legend.position="bottom")
ggsave("figures/fig06_volatility.png", p, width=6.5, height=9, dpi=500, bg="white")
cat("wrote figures/fig06_volatility.png\n")

## Title + caption off the PNG, into the shared captions.csv (same convention as the other fig*.R)
write_captions <- function(new_caps){
  path <- "figures/captions.csv"
  old <- if(file.exists(path)) fread(path) else data.table(file=character(),title=character(),caption=character())
  fwrite(rbind(old[!file %in% new_caps$file], new_caps), path)
  cat("wrote", path, "\n")
}
write_captions(data.table(
  file="figures/fig06_volatility.png",
  title="When in the open-water season temperature and oxygen are changing",
  caption=paste0(
    "Seasonal summary of the depth-by-period trends in Figure 2 (Theil-Sen slopes of seasonally ",
    "detrended residuals, 1 April-15 November; northern lakes 1981-2025, southern 1995-2025). ",
    "Tiles: mean trend across all depth bins for each lake and half-month period (temperature in °C ",
    "per decade; dissolved oxygen in mg/L per decade). Dots: percentage of that lake's depth range lying in ",
    "a significant depth block (5 m in lakes deeper than 10 m, 2 m in lakes of 5-10 m, 1 m in Crystal Bog; Mann-Kendall ",
    "on block-mean residuals, Benjamini-Hochberg q<0.05 within each lake; as in the block version of Figure 2). ",
    "Top bars: the same percentage in each half-month, averaged across lakes so each lake counts ",
    "equally. (a) Every depth bin weighted equally. (b) Each depth bin weighted by the lake's area at that depth ",
    "(layer volume fractions from lake hypsometry), so tiles are ",
    "volume-weighted mean trends and dots the percentage of lake volume in significant blocks. Averaging over depth lets opposite trends cancel (e.g. surface warming and deep cooling ",
    "in Crystal Lake in September-October), so a pale tile with a large dot marks a period of change ",
    "in both directions. Crystal Lake 2012-2013 excluded (whole-lake mixing experiment)."
  )))
