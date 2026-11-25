## Fig 3: which parts of the open-water season are changing most? A seasonal summary of Figure 2
## (Figure2_rate_of_change.R), all 9 NTL lakes (Fish Lake and Wingra excluded).
## Reads figures/fig02_cell_trends.csv and figures/fig02_block_trends.csv -- run Figure2_rate_of_change.R first.
##
## For each lake x half-month period, collapsing over depth:
##   tile fill = mean Theil-Sen slope across that lake's depth bins (signed; degC or mg/L per decade)
##   dot size  = % of that lake's depth bins lying in a SIGNIFICANT DEPTH BLOCK (the block test of the
##               fig02_rate_block_*.png versions: BH q<0.05 per lake), so Figures 2 and 3 agree
## Top bars: that same %, split into significant INCREASING blocks (above 0) and DECREASING blocks
## (below 0), in each half-month, averaged across lakes so each lake counts equally (otherwise Trout,
## with the most depth bins, would dominate).
## Two outputs: figures/fig03_volatility.png weights every depth bin equally (main figure);
## figures/fig03_volatility_volume.png weights each depth bin by the lake AREA at that depth x bin thickness,
## i.e. the share of lake VOLUME it represents (data/Bathymetry/lakehypsometry.csv), so tiles become
## volume-weighted mean trends and dots the % of lake VOLUME in significant blocks (robustness check).
## Shallow, broad layers then count for more than the narrow deep hole.
## Caveat: averaging over depth lets opposite trends cancel (e.g. Crystal's surface warming and deep
## cooling in Sep-Oct give a pale tile) -- the dot still shows that the period is active.
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(scales)})

path <- "figures/fig02_cell_trends.csv"
if(!file.exists(path)) stop(path, " not found -- run Figure2_rate_of_change.R first")
w <- fread(path)[var %in% c("wtemp","o2")]
bpath <- "figures/fig02_block_trends.csv"
if(!file.exists(bpath)) stop(bpath, " not found -- run Figure2_rate_of_change.R first")
bk <- fread(bpath)[var %in% c("wtemp","o2")]
## each 1-m (0.5-m) cell inherits the significance of the depth block it sits in
w[bk, on=.(lakeid, var, tbin, dbin >= depth_top, dbin <= depth_bottom), `:=`(bsig = i.sig, bslope = i.slope)]
w[is.na(bsig), bsig := FALSE]

## same three groups and order as Figure 2 (Figure2_rate_of_change.R): northern lakes with no clarity
## change, northern lakes with a browning trend, southern lakes; descending max depth within a group
w[, grp := factor(fcase(lakeid %in% c("CR","TB","CB"), "brown", region == "Northern", "clear", default="south"),
                  levels=c("clear","brown","south"))]
ord <- unique(w[, .(name, grp, zmax)])[order(grp, -zmax)]$name
w[, name := factor(name, levels=rev(ord))]           # first lake at the top of the reversed axis
## y positions of the browning group's edges (counted from the bottom of the reversed axis)
GRP_BREAKS <- cumsum(rev(w[, uniqueN(name), keyby=grp]$V1))[1:2] + 0.5

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
  ## pct = % in significant blocks (dot size); up/down split it by the sign of the block's trend
  s <- x[, .(net=sum(slope*ww)/sum(ww), pct=pmin(100, 100*sum(bsig*ww)/sum(ww)),     # pmin: float rounding can give 100.0000001
             up=100*sum((bsig & bslope > 0)*ww)/sum(ww), down=100*sum((bsig & bslope < 0)*ww)/sum(ww)),
         by=.(var, name, tbin, twidth)]
  list(s=s, top=s[, .(pct=mean(pct), up=mean(up), down=mean(down)), by=.(var, tbin, twidth)])
}
A <- summ(FALSE); B <- summ(TRUE)
for(z in list(list("Fig 3 (fig03_volatility): each depth bin equal", A), list("Fig 3 (fig03_volatility_volume): weighted by lake volume", B))){
  cat("\n", z[[1]], "-- % in significant increasing / decreasing blocks by half-month (mean across lakes):\n")
  print(dcast(z[[2]]$top[, .(var, tbin, v=sprintf("+%d/-%d", round(up), round(down)))], var ~ tbin, value.var="v")) }

## month boundaries + tick centres, same convention as Figure2_rate_of_change.R
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
## bar colours: the increase / decrease ends of each variable's palette
BARCOL <- list(wtemp=c(up="#d6604d", down="#4393c3"), o2=c(up="#35978f", down="#bf812d"))
YMAX <- ceiling(max(c(A$top$up,   B$top$up))/10)*10
YMIN <- ceiling(max(c(A$top$down, B$top$down))/10)*10

column <- function(vv, ttl, pal, left, S, sizename, tag=NULL){
  s <- S$s; top <- S$top
  ## top bars: % of depth (or volume) in significant blocks that are increasing (above 0) or
  ## decreasing (below 0), mean across lakes
  bars <- ggplot(top[var==vv]) +
    geom_col(aes(tbin, up, width=twidth*0.9), fill=BARCOL[[vv]][["up"]]) +
    geom_col(aes(tbin, -down, width=twidth*0.9), fill=BARCOL[[vv]][["down"]]) +
    geom_hline(yintercept=0, color="grey30", linewidth=0.3) +
    ## y title drawn as text just left of the tick labels: a normal axis title would sit left of the
    ## (aligned) column holding the tiles' wide lake-name labels, far from this axis
    { if(left) annotate("text", x=-Inf, y=(YMAX - YMIN)/2, label="% significant", angle=90, vjust=-3.3, size=2.5) } +
    coord_cartesian(clip="off") +
    geom_vline(xintercept=mbound, color="white", linewidth=0.4) +
    annotate("text", x=113, y=0.7*YMAX,  label="Increase", hjust=0, size=2.1, fontface="italic", color="grey30") +
    annotate("text", x=113, y=-0.75*YMIN, label="Decrease", hjust=0, size=2.1, fontface="italic", color="grey30") +
    xsc + scale_y_continuous(limits=c(-YMIN, YMAX), expand=c(0,0), labels=function(b) abs(b)) +
    labs(title=ttl, tag=tag, x=NULL, y=NULL) + th +
    theme(plot.tag=element_text(face="bold", size=9)) +
    theme(axis.text.x=element_blank(), legend.position="none") +
    { if(!left) theme(axis.text.y=element_blank()) }
  tiles <- ggplot(s[var==vv], aes(tbin, name)) +
    geom_tile(aes(width=twidth, fill=net), height=0.95) +
    geom_point(data=s[var==vv & pct>0], aes(size=pct), shape=21, fill="white", color="grey20", stroke=0.25) +
    geom_vline(xintercept=mbound, color="white", linewidth=0.4) +
    ## dashed box around the browning lakes; labelled "Browning" outside the right-hand column
    annotate("rect", xmin=91, xmax=320, ymin=GRP_BREAKS[1], ymax=GRP_BREAKS[2],
             fill=NA, color="grey20", linetype="dashed", linewidth=0.35) +
    { if(!left) annotate("text", x=Inf, y=mean(GRP_BREAKS), label="Browning", angle=-90, vjust=-0.6,
                         size=2.3, fontface="italic", color="grey20") } +
    coord_cartesian(clip="off") +
    xsc + pal +
    scale_size_area(max_size=3, limits=c(0,100), breaks=c(10,30,60), name=sizename) +
    guides(fill=guide_colourbar(order=1, title.position="top"), size=guide_legend(order=2, title.position="top")) +
    labs(x=NULL, y=NULL) + th +
    { if(!left) theme(axis.text.y=element_blank()) }
  bars / tiles + plot_layout(heights=c(1,4))
}
## two separate figures (same bar axis limits, so they can be compared directly):
##   fig03_volatility.png         every depth bin weighted equally
##   fig03_volatility_volume.png  the same, weighted by lake volume
save_fig <- function(S, ttl_T, ttl_O, sizename, outfile){
  p <- (column("wtemp", ttl_T, pal_T, TRUE, S, sizename) | column("o2", ttl_O, pal_O, FALSE, S, sizename)) +
    plot_layout(guides="collect") & theme(legend.position="bottom")
  p <- p + plot_annotation(theme=theme(plot.margin=margin(5, 14, 5, 5)))   # room for the 'Browning' label
  ggsave(outfile, p, width=6.5, height=4.9, dpi=500, bg="white"); cat("wrote", outfile, "\n")
}
save_fig(A, "Temperature", "Dissolved oxygen", "% of depth range\nin significant blocks", "figures/fig03_volatility.png")
save_fig(B, "Temperature, weighted by lake volume", "Dissolved oxygen, weighted by lake volume",
         "% of lake volume\nin significant blocks", "figures/fig03_volatility_volume.png")

## Title + caption off the PNG, into the shared captions.csv (same convention as the other fig*.R)
write_captions <- function(new_caps){
  path <- "figures/captions.csv"
  old <- if(file.exists(path)) fread(path) else data.table(file=character(),title=character(),caption=character())
  fwrite(rbind(old[!file %in% new_caps$file], new_caps), path)
  cat("wrote", path, "\n")
}
common <- paste0(
  "Seasonal summary of the depth-by-period trends in Figure 2 (Theil-Sen slopes of seasonally ",
  "detrended residuals, 1 April-15 November; northern lakes 1981-2025, southern 1995-2025). ",
  "Lakes are grouped as in Figure 2; the dashed box marks the northern lakes with a browning trend. ")
dots <- paste0("lying in a significant depth block (5 m in lakes deeper than 10 m, 2 m in lakes of 5-10 m, 1 m in ",
  "Crystal Bog; Mann-Kendall on block-mean residuals, Benjamini-Hochberg q<0.05 within each lake; as in the ",
  "block version of Figure 2). ")
bars <- paste0("Top bars: the same percentage split by direction (significant increases above the line, decreases ",
  "below; colours as in the tile scales) in each half-month, averaged across lakes so each lake counts equally. ")
tail_txt <- paste0("Averaging over depth lets opposite trends cancel (e.g. surface warming and deep cooling in ",
  "Crystal Lake in September-October), so a pale tile with a large dot marks a period of change in both ",
  "directions. Crystal Lake 2012-2013 excluded (whole-lake mixing experiment).")
write_captions(rbind(
  data.table(file="figures/fig03_volatility.png",
    title="When in the open-water season temperature and oxygen are changing",
    caption=paste0(common,
      "Tiles: mean trend across all depth bins for each lake and half-month period (temperature in °C per ",
      "decade; dissolved oxygen in mg/L per decade), each depth bin weighted equally. Dots: percentage of that ",
      "lake's depth range ", dots, bars, tail_txt)),
  data.table(file="figures/fig03_volatility_volume.png",
    title="Seasonal summary of trends weighted by lake volume",
    caption=paste0("As Figure 3, but with each depth bin weighted by its share of lake volume (lake area at ",
      "that depth times bin thickness, from layer volume fractions in lake hypsometry; Fig. SX hypsography), so ",
      "shallow, broad layers count for more than the narrow deep basin. ", common,
      "Tiles: volume-weighted mean trend for each lake and half-month period. Dots: percentage of lake volume ",
      dots, bars, tail_txt))))
