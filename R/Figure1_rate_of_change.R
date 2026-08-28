## Fig 12: RATE of change across depth x season, all 11 NTL lakes. Data: EDI 29.
## Companion to fig09, but instead of differencing two eras this fits a trend through the WHOLE
## record: per cell, a Theil-Sen slope (per decade) with a Mann-Kendall significance test.
## Why this exists: an era split (fig09) discards the time ordering -- a step change and a steady
## ramp look identical -- and throws away the years in between. A Sen slope uses every year, is
## robust to outliers (same motivation as the medians in fig09), and yields an interpretable RATE.
suppressMessages({library(data.table); library(ggplot2); library(scales); library(patchwork)
                  library(ggh4x); library(splines); library(quantreg); library(ggtext)})

meta <- data.table(
  lakeid=c("TR","BM","CR","SP","AL","TB","CB","ME","MO","FI","WI"),
  name  =c("Trout","Big Muskie","Crystal","Sparkling","Allequash","Trout Bog","Crystal Bog",
           "Mendota","Monona","Fish","Wingra"),
  region=c(rep("Northern",7), rep("Southern",4)),
  zmax  =c(35.7,21.3,20.4,20.0,8.0,7.9,2.5,25.3,22.5,18.9,4.0))

prof <- fread("data/profiles_clean.csv")[lakeid %in% meta$lakeid]
prof[, `:=`(date=as.Date(sampledate), doy=yday(as.Date(sampledate)), year=year4)]
prof <- prof[doy>=91 & doy<=319 & !(lakeid=="CR" & year %in% c(2012,2013))]  # Apr 1 - Nov 15

## long over the variables, interpolate each cast to 0.5 m, half-month bins (same as fig09)
pm <- melt(prof[,.(lakeid,date,doy,year,depth,wtemp,o2,o2sat)],
           id.vars=c("lakeid","date","doy","year","depth"), variable.name="var", value.name="val")
pm <- pm[is.finite(val)]
g <- pm[, { o<-order(depth); dd<-depth[o]; yy<-val[o]
  if(length(dd)>=3 && diff(range(dd))>=1){ zg<-seq(0,max(dd),0.5); zg<-zg[zg>=min(dd)]
    .(dbin=zg, val=approx(dd,yy,zg,rule=2)$y) } else .(dbin=numeric(0),val=numeric(0)) },
  by=.(lakeid,date,doy,year,var)]
hedges <- c(91,106,121,136,152,167,182,197,213,228,244,259,274,289,305,320)  # Apr 1 - Nov 16
g[, bi:=findInterval(doy,hedges)]; g<-g[bi>=1 & bi<=length(hedges)-1]
g[, `:=`(tbin=(hedges[bi]+hedges[bi+1])/2, twidth=hedges[bi+1]-hedges[bi])]

## Seasonal detrending, exactly as in fig09 and for the same reason: within a half-month window the
## days actually visited drift later/earlier over the record, and since temp & O2 change fast within
## a fortnight that drift alone can manufacture a trend. Fit one median (tau=0.5) seasonal curve per
## lake/var/depth across the whole Apr-Nov season, then trend the RESIDUALS from it.
## df=6 suits this ~229-day window; a narrower window would need fewer df to avoid the
## 'seasonal' curve absorbing real interannual signal.
g[, resid := {
  if(.N>=15){
    fit <- tryCatch(rq(val ~ ns(doy, df=6), tau=0.5), error=function(e) NULL)
    if(!is.null(fit)) residuals(fit) else val-median(val)
  } else if(.N>=6){
    fit <- tryCatch(rq(val ~ doy, tau=0.5), error=function(e) NULL)
    if(!is.null(fit)) residuals(fit) else val-median(val)
  } else val-median(val)
}, by=.(lakeid,var,dbin)]

## one value per YEAR per cell first (median of that year's casts), so a year sampled twice in a
## window doesn't get double weight in the slope
yr <- g[, .(v=median(resid)), by=.(lakeid,var,tbin,twidth,dbin,year)]

## Theil-Sen slope: median of all pairwise slopes. Mann-Kendall (Kendall tau) for significance --
## the natural non-parametric partner to Sen, and the same pairing already used in fig08/fig11.
MIN_YR_TREND <- 10   # distinct years required to fit a trend at all
MIN_SPAN     <- 20   # ...spread over at least this many years, so it's a trend not a short burst
sen <- function(y,x){ n<-length(y); if(n<3) return(NA_real_)
  d <- outer(y,y,"-")[lower.tri(diag(n))] / outer(x,x,"-")[lower.tri(diag(n))]
  median(d[is.finite(d)], na.rm=TRUE) }
cell <- yr[, { n<-.N; sp<- if(n>0) diff(range(year)) else 0
  if(n>=MIN_YR_TREND && sp>=MIN_SPAN){
    .(slope = sen(v, year)*10,      # per decade
      p     = tryCatch(suppressWarnings(cor.test(year, v, method="kendall", exact=FALSE)$p.value),
                       error=function(e) NA_real_), nyr=n)
  } else .(slope=NA_real_, p=NA_real_, nyr=n) },
  by=.(lakeid,var,tbin,twidth,dbin)]
w <- cell[!is.na(slope)]
w[, sig := !is.na(p) & p<0.05]
w <- merge(w, meta, by="lakeid")

## ---- context arrows beside each lake name: surface DOC trend and Secchi trend -----------------
## Pulled from figures/fig02_trends_table.csv (written by Figure2_clarity_trends.R) instead of
## re-fit here -- one source of truth for the DOC/Secchi trend, not two slightly-different fits.
## NOTE this changes the window from this figure's own Apr-mid Nov to fig02's June-Sept
## (stratified season), and Secchi source from secnview-for-all-11 to fig02's region split
## (secview north, secnview south) -- run Figure2_clarity_trends.R before this script.
trends_path <- "figures/fig02_trends_table.csv"
if(!file.exists(trends_path)) stop("figures/fig02_trends_table.csv not found -- run Figure2_clarity_trends.R first")
t2 <- fread(trends_path)[metric %in% c("DOC","Secchi")]
## fig02's slope is a DARKNESS-oriented z-score (sign-flipped for Secchi so + = darker); un-flip
## Secchi back to its own raw direction (+ = deeper/clearer) to match this figure's "up = that
## variable increased" glyph convention. DOC needs no flip (raw direction already = darker).
t2[, raw_slope := fifelse(metric=="Secchi", -slope_mad_per_decade, slope_mad_per_decade)]
arrows <- dcast(t2, lakeid ~ metric, value.var=c("raw_slope","p"))
setnames(arrows, c("raw_slope_DOC","p_DOC","raw_slope_Secchi","p_Secchi"),
         c("doc_s","doc_p","sec_s","sec_p"))
## brown up / blue down / grey X, coloured by DIRECTION of that variable (not by "darkening"):
## for DOC, up = more carbon; for Secchi, up = clearer water. Labelled so the two can't be confused.
glyph <- function(s,p) fifelse(!is.na(p) & p<0.05 & s>0,
    "<span style='color:#a1622f;font-weight:bold;font-size:12pt'>&#8593;</span>",
  fifelse(!is.na(p) & p<0.05 & s<0,
    "<span style='color:#1f6f9e;font-weight:bold;font-size:12pt'>&#8595;</span>",
    "<span style='color:#9a9a9a;font-weight:bold;font-size:9pt'>&#10005;</span>"))
arrows[, lab_md := sprintf("<span style='font-size:7pt'>DOC</span>%s&nbsp;&nbsp;<span style='font-size:7pt'>Sec</span>%s",
        glyph(doc_s,doc_p), glyph(sec_s,sec_p))]

w <- merge(w, arrows[, .(lakeid, lab_md)], by="lakeid", all.x=TRUE)
w[, striplab := sprintf("<span style='font-size:7pt;font-weight:bold'>%s</span><br>%s", name, fifelse(is.na(lab_md), "", lab_md))]
ord <- unique(w[, .(name, striplab, region, zmax)])[order(region,-zmax)]$striplab
w[, strip := factor(striplab, levels=ord)]

## Non-significant cells get a light diagonal stroke (one per tile). This is deliberately a plain
## geom_segment rather than ggpattern: ggpattern rasterizes every tile separately and takes many
## minutes at these tile counts (~14k cells), whereas segments stay vector and render instantly.
hatch <- function(d) d[sig==FALSE, .(x=tbin-twidth/2, xend=tbin+twidth/2, y=dbin-0.25, yend=dbin+0.25, strip)]

## diverging palettes centred on zero rate; limits from the observed 5-95% spread, extremes squished
pal_T <- scale_fill_gradientn(colours=c("#2166ac","#4393c3","#92c5de","#d1e5f0","#f7f7f7","#fddbc7","#f4a582","#d6604d","#b2182b"),
           limits=c(-1,1), oob=squish, name="Δ Temp (°C/decade)", breaks=c(-1,-0.5,0,0.5,1))
brbg <- c("#8c510a","#bf812d","#dfc27d","#f6e8c3","#f5f5f5","#c7eae5","#80cdc1","#35978f","#01665e")
pal_O    <- scale_fill_gradientn(colours=brbg, limits=c(-1.5,1.5), oob=squish, breaks=c(-1.5,-0.75,0,0.75,1.5),
              name="Δ DO (mg/L/decade)")
pal_Osat <- scale_fill_gradientn(colours=brbg, limits=c(-15,15), oob=squish, breaks=c(-15,-7.5,0,7.5,15),
              name="Δ DO (% sat/decade)")
mbound <- c(91,121,152,182,213,244,274,305,320)
xsc <- scale_x_continuous(breaks=c(106,136.5,167,197.5,228.5,259,289.5,312.5),
                          labels=c("Apr","May","Jun","Jul","Aug","Sep","Oct","Nov"),
                          minor_breaks=NULL, expand=c(0,0))
## text sizes scaled down for the 6.5in-wide output (was 12in) so panel titles like
## "Northern Forested Lakes" fit without clipping
th <- theme_minimal(base_size=6.5) + theme(
  panel.grid=element_blank(), panel.spacing=unit(0.25,"lines"),
  legend.position="bottom", legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.22,"cm"),
  legend.title=element_text(size=6.5), plot.title=element_text(face="bold",size=9.5),
  plot.subtitle=element_text(face="bold", size=6.5, color="grey30"),
  axis.text=element_text(size=5), axis.text.x=element_text(size=6),
  axis.title.y=element_text(size=7.5),
  strip.text.y.right=element_markdown(angle=0, hjust=0, size=5, lineheight=1.2))

block <- function(reg, vv, pal, right, ttl, sub){
  d <- w[region==reg & var==vv]
  ggplot(d, aes(tbin,dbin,fill=slope)) +
    geom_tile(aes(width=twidth), height=0.5) +
    geom_segment(data=hatch(d), aes(x=x,xend=xend,y=y,yend=yend), inherit.aes=FALSE,
                 color="grey25", linewidth=0.16, alpha=0.5) +
    geom_vline(xintercept=mbound, color="white", linewidth=0.3, alpha=0.8) +
    facet_grid2(strip ~ ., scales="free_y", switch=NULL, render_empty=FALSE) +
    scale_y_reverse(expand=c(0,0), breaks=scales::breaks_pretty(3)) + xsc + pal + th +
    labs(title=ttl, subtitle=sub, x=NULL, y=if(!right) "Depth (m)" else NULL) +
    { if(!right) theme(strip.text.y.right=element_blank()) } +
    { if(right) theme(axis.text.y=element_blank(), axis.title.y=element_blank()) }
}
## Figure title + long explanatory paragraph are NOT drawn on the PNG (kept out of the image so
## the figure can be captioned in a manuscript/report instead) -- collected here and written to
## the shared figures/captions.csv (one row per output PNG, across all fig*.R scripts) alongside
## the panel-level labels ("Northern Forested Lakes...", "Trend in temperature", etc.), which DO
## stay on the PNG since they're needed to tell the panels apart.
FIG_TITLE <- "Rate of change across Wisconsin's lakes: full-record trends by depth & season"
captions <- list()

## merge-write so this script's rows don't clobber captions written by other fig*.R scripts
write_captions <- function(new_caps){
  path <- "figures/captions.csv"
  old <- if(file.exists(path)) fread(path) else data.table(file=character(),title=character(),caption=character())
  fwrite(rbind(old[!file %in% new_caps$file], new_caps), path)
  cat("wrote", path, "\n")
}

make_fig <- function(dovar, dopal, dolab, subtitle, outfile){
  NT <- block("Northern","wtemp", pal_T, FALSE, "Northern Forested Lakes", "Trend in temperature")
  NO <- block("Northern", dovar,  dopal, TRUE,  " ", dolab)
  ST <- block("Southern","wtemp", pal_T, FALSE, "Southern Agricultural/Urban Lakes", "Trend in temperature")
  SO <- block("Southern", dovar,  dopal, TRUE,  " ", dolab)
  left  <- (NT / ST) + plot_layout(heights=c(7,4), guides="collect") & theme(legend.position="bottom")
  right <- (NO / SO) + plot_layout(heights=c(7,4), guides="collect") & theme(legend.position="bottom")
  design <- (left | right) + plot_layout(guides="keep")
  final <- wrap_elements(design)
  ggsave(outfile, final, width=6.5, height=8, dpi=500, bg="white"); cat("wrote", outfile, "\n")
  captions[[outfile]] <<- data.table(file=outfile, title=FIG_TITLE, caption=subtitle)
}

yrs <- range(prof$year)
base <- sprintf("April–mid-November. Theil-Sen slope per DECADE over the full %d–%d record (not a two-era difference): every year contributes, so a steady ramp and a step change are no longer confounded.\nCells are hatched where the trend is NOT significant (Mann-Kendall p<0.05); plain cells are significant. Blank = fewer than %d years or <%d-year span, too little to fit a trend.\nValues are trends in the seasonally-detrended residual, so drift in which days of a fortnight got sampled cannot masquerade as a trend.\nBeside each lake name: full-record trend in surface DOC and in Secchi depth over the same window — brown ↑ rising, blue ↓ falling, grey ✕ no significant trend (Mann-Kendall p<0.05).",
   yrs[1], yrs[2], MIN_YR_TREND, MIN_SPAN)
make_fig("o2",    pal_O,    "Trend in dissolved oxygen (mg/L)",
  paste0(base, "\nOxygen trends are strongly skewed negative — deep water is losing O₂ across nearly every lake, while surface/metalimnetic waters hold steadier."), "figures/fig01_rate_mgL.png")
make_fig("o2sat", pal_Osat, "Trend in dissolved oxygen (% sat)",
  paste0(base, "\nAs % saturation: warming lowers O₂ solubility, so the same mg/L loss reads as a larger saturation deficit in deep water."), "figures/fig01_rate_sat.png")

write_captions(rbindlist(captions))

## Console summary. NOTE: deliberately summarises ALL cells, not just the significant ones.
## Taking the median over only p<0.05 cells is selection-biased (winner's curse) -- conditioning on
## significance preferentially keeps large-magnitude slopes. For Wingra wtemp that inflated the
## median from +0.24 to +1.34 C/decade, ~5.6x. Report the unconditional rate, and report the
## fraction significant separately as the evidence-strength measure.
w[, zone := fifelse(dbin <= 2, "surface(<=2m)", fifelse(dbin >= 0.6*zmax, "deep(>=60% zmax)", "mid"))]
sm <- w[zone!="mid", .(slope_per_decade=round(median(slope),3), pct_sig=round(100*mean(sig)), n=.N),
        by=.(region,name,var,zone)]
cat("\n=== median trend per decade over ALL cells (unconditional), by lake & zone ===\n")
print(dcast(sm, region+name+zone ~ var, value.var="slope_per_decade")[order(region,name,zone)], nrows=60)
cat("\n=== % of cells with a significant trend (evidence strength, reported separately) ===\n")
print(dcast(sm, region+name+zone ~ var, value.var="pct_sig")[order(region,name,zone)], nrows=60)
