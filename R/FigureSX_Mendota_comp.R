## FigureSX_Mendota_comp.R: Lake Mendota dissolved oxygen (mg/L), 1906 vs 2016-2025.
##
## data/Juday/Mendota_1906.csv is a hand-digitized historical record (Juday era UW limnology
## survey) -- repeated casts at four FIXED depths (0, 8, 10, 22 m) across the 1906 open-water
## season, not full continuous profiles. To keep the two eras directly comparable, BOTH panels
## here are built at those same four depths: the modern 2016-2025 panel interpolates every cast
## to exactly 0/8/10/22 m rather than the fine continuous depth grid FigureSX_do_heatmaps.R uses
## elsewhere -- showing modern data at a finer resolution than 1906 ever had would imply a
## precision the historical panel can't match and make the comparison look apples-to-oranges.
## Same half-month time bins and DO colour scale as FigureSX_do_heatmaps.R for visual consistency
## with the rest of the series.
suppressMessages({library(data.table); library(ggplot2); library(scales)})

TARGET_DEPTHS <- c(0, 8, 10, 22)
MIN_YEARS <- 4   # modern cells need data from at least this many distinct years to be averaged

## ---- 1906: four fixed depths, repeated casts across the season -------------------------------
d1906 <- fread("data/Juday/Mendota_1906.csv")[, .(doy=yday, depth=depth_m, val=do_mgl)]
d1906[, era := "1906"]

## ---- 2016-2025: interpolate each cast to the SAME four depths --------------------------------
prof <- fread("data/profiles_clean.csv")[lakeid=="ME" & !is.na(o2)]
prof[, `:=`(doy=yday(as.Date(sampledate)), year=year4)]
prof <- prof[year>=2016 & year<=2025]

## only interpolate a target depth if it's within 2 m of the cast's own sampled range -- don't
## extrapolate a 22 m estimate from a cast that only reached 12 m
modern <- prof[, {
  o <- order(depth); dd <- depth[o]; yy <- o2[o]
  if (length(dd)>=2) {
    est <- approx(dd, yy, TARGET_DEPTHS, rule=2)$y
    ok <- TARGET_DEPTHS >= (min(dd)-2) & TARGET_DEPTHS <= (max(dd)+2)
    est[!ok] <- NA_real_
    .(depth=TARGET_DEPTHS, val=est)
  } else .(depth=numeric(0), val=numeric(0))
}, by=.(date=as.Date(sampledate), doy, year)]
modern <- modern[!is.na(val)]
modern[, era := "2016-2025"]

## ---- half-month time bins, identical edges to FigureSX_do_heatmaps.R / Figure1 ---------------
hedges <- c(91,106,121,136,152,167,182,197,213,228,244,259,274,289,305,320)
mbound <- c(91,121,152,182,213,244,274,305,320)
xsc <- scale_x_continuous(breaks=c(106,136.5,167,197.5,228.5,259,289.5,312.5),
                          labels=c("Apr","May","Jun","Jul","Aug","Sep","Oct","Nov"),
                          minor_breaks=NULL, expand=c(0,0))

bin_season <- function(d) {
  d <- copy(d)
  d[, bi := findInterval(doy, hedges)]
  d <- d[bi>=1 & bi<=length(hedges)-1]
  d[, `:=`(tbin=(hedges[bi]+hedges[bi+1])/2, twidth=hedges[bi+1]-hedges[bi])]
  d
}
d1906 <- bin_season(d1906)
modern <- bin_season(modern)

## 1906: mean of however many casts fall in a cell (usually 0 or 1, occasionally 2)
agg1906 <- d1906[, .(m=mean(val)), by=.(era, tbin, twidth, depth)]

## modern: mean across years, blanked out if fewer than MIN_YEARS distinct years contributed
agg_modern <- modern[, .(m=mean(val), n_yr=uniqueN(year)), by=.(era, tbin, twidth, depth)]
agg_modern[n_yr < MIN_YEARS, m := NA]
agg_modern[, n_yr := NULL]

agg <- rbind(agg1906, agg_modern)
agg[, era := factor(era, levels=c("1906","2016-2025"))]
## discrete depth rows (0 m at top), not a continuous axis -- 1906 never resolved anything finer
agg[, depth_lab := factor(sprintf("%d m", depth), levels=sprintf("%d m", rev(TARGET_DEPTHS)))]

pal_do <- scale_fill_gradientn(colours=c("#000000","#7f0000","#d7301f","#fdae61","#ffffbf","#66bd63","#1a9850"),
            values=rescale(c(0,1,2,4,6,9,14)), name="DO (mg/L)", limits=c(0,14), oob=squish, na.value="grey90")

g <- ggplot(agg, aes(tbin, depth_lab, fill=m)) +
  geom_tile(aes(width=twidth), height=0.92, color="white", linewidth=0.3) +
  geom_vline(xintercept=mbound, color="white", linewidth=0.35, alpha=0.85) +
  facet_wrap(~era, ncol=1, strip.position="left") +
  xsc + pal_do +
  labs(title="Lake Mendota dissolved oxygen: 1906 vs. 2016–2025", x=NULL, y="Depth") +
  theme_minimal(base_size=12) +
  theme(panel.grid=element_blank(), plot.title=element_text(face="bold", size=13.5),
        strip.text=element_text(face="bold", size=11), strip.placement="outside",
        axis.text=element_text(size=10), panel.spacing=unit(0.8,"lines"))

ggsave("figures/figSX_mendota_1906_comp.png", g, width=7, height=4.2, dpi=500, bg="white")
cat("wrote figures/figSX_mendota_1906_comp.png\n")

## Title + caption off the PNG, into the shared captions.csv (same convention as the other fig*.R)
write_captions <- function(new_caps) {
  path <- "figures/captions.csv"
  old <- if (file.exists(path)) fread(path) else data.table(file=character(), title=character(), caption=character())
  fwrite(rbind(old[!file %in% new_caps$file], new_caps), path)
  cat("wrote", path, "\n")
}
write_captions(data.table(
  file="figures/figSX_mendota_1906_comp.png",
  title="Lake Mendota dissolved oxygen, 1906 vs. 2016-2025",
  caption=paste0(
    "Dissolved oxygen (mg/L) by half-month period at four fixed depths (0, 8, 10, 22 m), comparing ",
    "a historical 1906 record (data/Juday/Mendota_1906.csv, hand-digitized from repeated casts across ",
    "the open-water season) against the modern 2016-2025 record. Both panels are built at the same four ",
    "depths for a direct comparison: 1906 only resolved these four fixed depths (not full continuous ",
    "profiles), so the modern panel interpolates each cast to the same depths rather than the finer ",
    "continuous grid used elsewhere in this project, to avoid implying a precision the historical ",
    "record never had. A modern cell is blank if fewer than 4 distinct years contributed data; a 1906 ",
    "cell is blank if no cast fell in that half-month window. Color scale shared with ",
    "FigureSX_do_heatmaps.R."
  )
))
