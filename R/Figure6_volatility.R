## Fig 6: which parts of the open-water season are changing most? A seasonal summary of Figure 2
## (Figure1_rate_of_change.R), all 9 NTL lakes (Fish Lake and Wingra excluded).
## Reads the per-cell trend table figures/fig01_cell_trends.csv -- run Figure1_rate_of_change.R first.
##
## For each lake x half-month period, collapsing over depth:
##   tile fill = mean Theil-Sen slope across that lake's depth bins (signed; degC or mg/L per decade)
##   dot size  = % of that lake's depth bins whose trend is significant (BH q<0.05, as in Figure 2)
## Top bars: % of cells significant in each half-month, averaged across lakes so each lake counts
## equally (otherwise Trout, with the most depth bins, would dominate).
## Caveat: averaging over depth lets opposite trends cancel (e.g. Crystal's surface warming and deep
## cooling in Sep-Oct give a pale tile) -- the dot still shows that the period is active.
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(scales)})

path <- "figures/fig01_cell_trends.csv"
if(!file.exists(path)) stop(path, " not found -- run Figure1_rate_of_change.R first")
w <- fread(path)[var %in% c("wtemp","o2")]

ord <- unique(w[, .(name, region, zmax)])[order(region, -zmax)]$name   # north first, then by zmax
w[, name := factor(name, levels=rev(ord))]

s   <- w[, .(net=mean(slope), pct=100*mean(sig)), by=.(var, name, tbin, twidth)]
top <- s[, .(pct=mean(pct)), by=.(var, tbin, twidth)]
cat("% of cells significant by half-month (mean across lakes):\n")
print(dcast(top[, .(var, tbin, pct=round(pct))], var ~ tbin, value.var="pct"))

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
YMAX <- ceiling(max(top$pct)/10)*10

column <- function(vv, ttl, pal, left){
  bars <- ggplot(top[var==vv], aes(tbin, pct)) +
    geom_col(aes(width=twidth*0.9), fill="grey40") +
    geom_vline(xintercept=mbound, color="white", linewidth=0.4) +
    xsc + scale_y_continuous(limits=c(0,YMAX), expand=c(0,0)) +
    labs(title=ttl, x=NULL, y=if(left) "% of cells\nsignificant" else NULL) + th +
    theme(axis.text.x=element_blank(), legend.position="none") +
    { if(!left) theme(axis.text.y=element_blank()) }
  tiles <- ggplot(s[var==vv], aes(tbin, name)) +
    geom_tile(aes(width=twidth, fill=net), height=0.95) +
    geom_point(data=s[var==vv & pct>0], aes(size=pct), shape=21, fill="white", color="grey20", stroke=0.25) +
    geom_vline(xintercept=mbound, color="white", linewidth=0.4) +
    xsc + pal +
    scale_size_area(max_size=3, limits=c(0,100), breaks=c(10,30,60), name="% of depth bins\nsignificant") +
    guides(fill=guide_colourbar(order=1, title.position="top"), size=guide_legend(order=2, title.position="top")) +
    labs(x=NULL, y=NULL) + th +
    { if(!left) theme(axis.text.y=element_blank()) }
  bars / tiles + plot_layout(heights=c(1,4))
}
p <- (column("wtemp", "Temperature", pal_T, TRUE) | column("o2", "Dissolved oxygen", pal_O, FALSE)) +
  plot_layout(guides="collect") & theme(legend.position="bottom")
ggsave("figures/fig06_volatility.png", p, width=6.5, height=4.8, dpi=500, bg="white")
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
    "per decade; dissolved oxygen in mg/L per decade). Dots: percentage of that lake's depth bins with ",
    "a significant trend (Mann-Kendall, Benjamini-Hochberg q<0.05 within each lake). Top bars: ",
    "percentage of cells significant in each half-month, averaged across lakes so each lake counts ",
    "equally. Averaging over depth lets opposite trends cancel (e.g. surface warming and deep cooling ",
    "in Crystal Lake in September-October), so a pale tile with a large dot marks a period of change ",
    "in both directions. Crystal Lake 2012-2013 excluded (whole-lake mixing experiment)."
  )))
