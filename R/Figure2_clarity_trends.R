## Fig 2: standardized water-clarity trends, all 11 NTL lakes. Data: EDI 1 (DOC), 29 (light/PAR),
## 31 (Secchi).
## THREE independently-measured views of clarity for the 7 NORTHERN lakes (they cross-validate
## each other), but only TWO for the 4 SOUTHERN lakes:
##   Kd        -- light-extinction coefficient (m^-1) from PAR profiles;  higher = darker  [north only]
##   DOC       -- dissolved organic carbon (mg/L), the lab-measured driver;  higher = darker  [all 11]
##   Secchi    -- Secchi disk depth (m), read by eye;                        lower  = darker  [all 11]
## Kd is northern-only: southern PAR profiles (frlight) only start in 2019 -- too short to trend.
## (CDOM colour / a440, EDI 87, is deliberately excluded from this analysis.)
##
## Every metric is standardized per lake as a ROBUST z-score -- (annual value - its own MEDIAN) /
## its own MAD (median absolute deviation), not mean/SD -- so one wild bloom-day or post-storm
## cast can't inflate the scale the rest of the record is measured against. Sign-flipped so that,
## for every metric, POSITIVE = DARKER: Kd/DOC already have that sign; Secchi is negated (low
## Secchi = dark = should be positive). This "darkness anomaly" lets metrics in different units,
## and lakes with only 2 of the 3, share ONE axis and one fixed scale across every panel --
## panels are directly comparable, which the old min-max rescale-onto-Kd's-raw-units approach did
## not allow. Plotted with UP = darker.
##
## Secchi is recorded twice per visit -- secview (with viewing scope) and secnview (no scope,
## reads ~0.5 m shallower, r=0.97). secview is essentially northern-only (southern lakes have 0-1
## readings each), so PRIMARY is secview for the 7 northern lakes and secnview for the 4 southern
## lakes.
suppressMessages({library(data.table); library(ggplot2); library(ggh4x); library(ggtext)})

SEASON  <- c(152, 273)  # Jun 1 - Sep 30, stratified season
MIN_OBS <- 3            # a lake-year is plotted only if >=3 sampling DATES fall in the window

## Annual value = MEDIAN across sampling dates (not mean): robust to the odd bloom-day or
## post-storm cast, consistent with the median-based treatment used in fig08/fig09/fig12.
## Counting DATES not rows matters -- DOC has several depths per sample, so a raw row count would
## let one visit satisfy the 3-observation threshold.
annual <- function(d, min_obs=MIN_OBS)
  d[, .(v=median(v), n=uniqueN(date)), by=.(lakeid,year)][n>=min_obs, .(lakeid,year,v)]

meta <- data.table(
  lakeid=c("TR","BM","CR","SP","AL","TB","CB","ME","MO","FI","WI"),
  name  =c("Trout","Big Muskie","Crystal","Sparkling","Allequash","Trout Bog","Crystal Bog",
           "Mendota","Monona","Fish","Wingra"),
  region=c(rep("Northern",7), rep("Southern",4)),
  zmax  =c(35.7,21.3,20.4,20.0,8.0,7.9,2.5,25.3,22.5,18.9,4.0))
NORTH <- meta[region=="Northern", lakeid]

## Theil-Sen slope + median intercept + Mann-Kendall p. Returns the intercept too so the drawn
## line IS the reported fit.
senfit <- function(y, x){
  ok <- is.finite(y) & is.finite(x); y <- y[ok]; x <- x[ok]; n <- length(y)
  if(n < 8) return(list(slope=NA_real_, int=NA_real_, p=NA_real_, n=n))
  sl <- median(outer(y,y,"-")[lower.tri(diag(n))] / outer(x,x,"-")[lower.tri(diag(n))], na.rm=TRUE)
  list(slope=sl, int=median(y - sl*x),
       p=suppressWarnings(cor.test(x, y, method="kendall")$p.value), n=n)
}
## no p<0.1 tier: a "." marker on a p=0.059 null is what produced a spurious "weak darkening" read
star <- function(p) ifelse(is.na(p),"", ifelse(p<0.001,"***", ifelse(p<0.01,"**", ifelse(p<0.05,"*",""))))

## robust z-score: median/MAD instead of mean/SD, so outliers don't inflate the scale
robz <- function(d) d[, .(year,v,z={m<-mad(v); if(is.finite(m) && m>0) (v-median(v))/m else v-median(v)}), by=lakeid]

## ---------------- metric 1: Kd (m^-1), slope of ln(fraction surface PAR) vs depth, thru origin ----
## Northern only -- see file header.
p <- fread("data/profiles_clean.csv")[lakeid %in% NORTH & !is.na(frlight) & frlight>0 & depth>0.1 &
                                !(lakeid=="CR" & year4 %in% c(2012,2013))]
p[, `:=`(year=year4, doy=yday(as.Date(sampledate)), date=as.Date(sampledate))]
p <- p[doy>=SEASON[1] & doy<=SEASON[2]]
kd <- p[, { v <- NA_real_
  if(.N>=3){ f <- tryCatch(lm(log(frlight) ~ depth + 0), error=function(e) NULL)
    if(!is.null(f)) v <- as.numeric(-coef(f)[1]) }
  list(v=v) }, by=.(lakeid, date, year)]           # one Kd per cast
kdyr <- robz(annual(kd[is.finite(v) & v>0]))

## ---------------- metric 2: DOC (mg/L), surface (0 m exactly), all 11 lakes ------------------------
## depth==0, NOT depth<=2: discrete sampled depths differ by lake, so <=2m is not like-for-like
## (see fig08/fig12). Same change applied there.
doc <- fread("data/chem_north.csv")[lakeid %in% meta$lakeid & doc>0 & doc<80 & depth==0 &
                                    (is.na(flagdoc) | flagdoc=="") &
                                    !(lakeid=="CR" & year4 %in% c(2012,2013))]
doc[, `:=`(year=year4, date=as.Date(sampledate), doy=yday(as.Date(sampledate)))]
docdate <- doc[doy>=SEASON[1] & doy<=SEASON[2], .(v=median(doc)), by=.(lakeid,date,year)]
docyr   <- robz(annual(docdate))

## ---------------- metric 3: Secchi (m), all 11 lakes -------------------------------------------------
sec <- fread("data/secchi.csv")[lakeid %in% meta$lakeid &
                                !(lakeid=="CR" & year4 %in% c(2012,2013))]
sec[, `:=`(year=year4, date=as.Date(sampledate), doy=yday(as.Date(sampledate)))]
sec <- sec[doy>=SEASON[1] & doy<=SEASON[2]]
secyr_view  <- robz(annual(sec[!is.na(secview)  & secview>0,  .(lakeid,date,year,v=secview)]))
secyr_nview <- robz(annual(sec[!is.na(secnview) & secnview>0, .(lakeid,date,year,v=secnview)]))
## secview for northern lakes, secnview for southern (secview is essentially unmeasured there)
secyr <- rbind(secyr_view[lakeid %in% NORTH], secyr_nview[!lakeid %in% NORTH])

## ---------------- assemble: darkness anomaly (robust z, sign-flipped so + = darker) -----------------
mk <- function(d, nm, flip=FALSE) d[, .(lakeid, year, metric=nm, darkness = if(flip) -z else z)]
dat <- rbind(mk(kdyr,"Kd"), mk(docyr,"DOC"), mk(secyr,"Secchi", flip=TRUE))

## ---------------- trends, fit directly on the darkness anomaly (a per-lake affine transform of the
## raw value, so Kendall's p and the sign of "getting darker" are unchanged; slope is now in
## MAD/decade, directly comparable across metrics and lakes) ------------------------------------------
tr <- function(nm, flip=FALSE){ d <- switch(nm, Kd=kdyr, DOC=docyr, Secchi=secyr)
  d[, { zz <- if(flip) -z else z; f<-senfit(zz, year); .(metric=nm, slope=f$slope*10, p=f$p, nyr=f$n) }, by=lakeid] }
fits <- rbind(tr("Kd"), tr("DOC"), tr("Secchi", flip=TRUE))
fw <- dcast(fits, lakeid ~ metric, value.var=c("slope","p","nyr"))
fw <- merge(fw, meta, by="lakeid")

## ---------------- output: trend table with p-values, all lakes x metrics --------------------------
tbl <- fits[meta, on="lakeid", nomatch=0][
  , .(lakeid, name, region, metric, slope_mad_per_decade=round(slope,3), p=round(p,4),
      sig=!is.na(p) & p<0.05, nyr)]
setorder(tbl, region, name, metric)
fwrite(tbl, "figures/fig02_trends_table.csv")
cat("wrote figures/fig02_trends_table.csv\n")
print(tbl, nrows=100)

## ---------------- strip labels: 3 metrics for northern, 2 for southern lakes -----------------------
fw[, striplab_n := sprintf(
      "<span style='font-size:7pt;font-weight:bold'>%s</span> (northern)<br>Kd %+.2f%s &middot; DOC %+.2f%s &middot; Secchi %+.2f%s",
      name, slope_Kd, star(p_Kd), slope_DOC, star(p_DOC), slope_Secchi, star(p_Secchi))]
fw[, striplab_s := sprintf(
      "<span style='font-size:7pt;font-weight:bold'>%s</span> (southern)<br>DOC %+.2f%s &middot; Secchi %+.2f%s",
      name, slope_DOC, star(p_DOC), slope_Secchi, star(p_Secchi))]
fw[, striplab := fifelse(region=="Northern", striplab_n, striplab_s)]
ord <- fw[order(region, name), striplab]
fw[, strip := factor(striplab, levels=ord)]
dat <- merge(dat, fw[,.(lakeid,strip,region)], by="lakeid")

## trend segments, fitted on the plotted darkness anomaly so the line matches what's drawn
segs <- dat[, { f<-senfit(darkness, year)
  .(x=min(year), xend=max(year), y=f$int+f$slope*min(year), yend=f$int+f$slope*max(year), p=f$p) },
  by=.(lakeid,strip,metric)]
segs[, solid := !is.na(p) & p<0.05]

lv <- c("Kd  (light extinction)","DOC  (mg/L)","Secchi  (disk depth)")
dat[,  metric := factor(metric, levels=c("Kd","DOC","Secchi"), labels=lv)]
segs[, metric := factor(metric, levels=c("Kd","DOC","Secchi"), labels=lv)]
pal <- setNames(c("#2b6a3d","#6a4c93","#1f6f9e"), lv)

th <- theme_minimal(base_size=6.5) + theme(
  panel.grid.minor=element_blank(), panel.spacing=unit(0.4,"lines"),
  legend.position="bottom", legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.22,"cm"),
  legend.title=element_text(size=6.5), plot.title=element_text(face="bold",size=9.5),
  plot.subtitle=element_text(face="bold", size=6.5, color="grey30"),
  axis.text=element_text(size=5), axis.text.x=element_text(size=6),
  axis.title.y=element_text(size=7.5),
  strip.text=element_markdown(size=5, lineheight=1.2))

g <- ggplot(dat, aes(year, darkness, color=metric)) +
  geom_hline(yintercept=0, color="grey80", linewidth=0.3) +
  geom_line(alpha=0.3, linewidth=0.35) +
  geom_point(size=0.7, alpha=0.75) +
  geom_segment(data=segs, aes(x=x, xend=xend, y=y, yend=yend, color=metric, linetype=solid),
               linewidth=0.8, inherit.aes=FALSE) +
  facet_wrap(~strip, ncol=3) +
  scale_color_manual(values=pal, name=NULL) +
  scale_linetype_manual(values=c(`TRUE`="solid",`FALSE`="22"), guide="none") +
  labs(x=NULL, y="Clarity anomaly (MAD from lake median) — up = darker") +
  th
ggsave("figures/fig02_clarity_trends.png", g, width=6.5, height=8, dpi=500, bg="white")
cat("wrote figures/fig02_clarity_trends.png\n")

## Title + long explanatory paragraph are NOT drawn on the PNG -- written to the shared
## figures/captions.csv instead (see Figure1_rate_of_change.R for the same convention).
write_captions <- function(new_caps){
  path <- "figures/captions.csv"
  old <- if(file.exists(path)) fread(path) else data.table(file=character(),title=character(),caption=character())
  fwrite(rbind(old[!file %in% new_caps$file], new_caps), path)
  cat("wrote", path, "\n")
}
caption <- paste0(
  "June-September, standardized trends: each lake's own annual series is centred on its own MEDIAN and scaled by its own MAD (median absolute deviation, ",
  "robust to outliers), then oriented so POSITIVE = DARKER for every metric (Secchi is sign-flipped since low Secchi = dark water) -- this puts metrics in ",
  "different raw units, and lakes with only 2 of the 3 available, on one shared, fixed axis. Northern lakes (7): Kd (light extinction), DOC, and Secchi ",
  "(secview). Southern lakes (4): DOC and Secchi (secnview) only -- southern PAR profiles only start in 2019 (too short to trend).\n",
  "Each point is the annual MEDIAN across sampling dates; lake-years with <3 sampling dates are dropped.\n",
  "Solid trend = significant (Theil-Sen slope / Mann-Kendall p<0.05); dashed = not significant. Labels: slope in MAD/decade (*<0.05 **<0.01 ***<0.001). ",
  "Crystal 2012-13 excluded (whole-lake mixing experiment). Full trend table with p-values: figures/fig02_trends_table.csv.")
write_captions(data.table(file="figures/fig02_clarity_trends.png",
  title="Standardized water-clarity trends across Wisconsin's lakes", caption=caption))
