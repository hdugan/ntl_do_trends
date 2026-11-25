## Paper Figure 2: RATE of change in temperature and dissolved oxygen across depth x season, for the
## 9 NTL lakes (Fish Lake and Wingra excluded), in three groups: northern lakes with no clarity
## change (Trout, Big Muskie, Sparkling, Allequash), northern lakes with a browning trend (Crystal,
## Trout Bog, Crystal Bog), and southern lakes (Mendota, Monona); descending max depth within a group.
## Input: data/profiles_clean.csv (from 02_clean_profiles.R): NTL-LTER profiles (EDI knb-lter-ntl.29)
##   plus, for Mendota only, extra 2006-2019 MEMO casts (EDI knb-lter-ntl.415; NTL wins on shared
##   days). Window Apr 1 - Nov 15; Crystal 2012-2013 excluded (whole-lake mixing experiment).
##   Run Figure2_clarity_trends.R first (supplies the DOC/Secchi arrows beside each lake name).
## Method, per lake x variable (wtemp, o2, o2sat):
##   - each cast interpolated to depth bins (1 m; 0.5 m for lakes < 10 m deep) in half-month bins
##   - DO < 0.5 mg/L floored to 0 (anoxic floor, below)
##   - residuals from a median seasonal curve (quantile regression, ns df=6) per lake/var/depth
##   - one median residual per year per cell; a cell needs >= 10 years and a >= 20-year span
##   - Theil-Sen slope per decade, Mann-Kendall test, Benjamini-Hochberg FDR per lake
##   - block version: the same test on depth blocks (5 m for lakes >= 10 m, 2 m for 5-10 m, 1 m
##     below) within each half-month; significant blocks are outlined
## A Sen slope uses every year (unlike a two-era difference, which cannot tell a step from a ramp),
## is robust to outliers, and gives an interpretable RATE.
## Outputs (figures/):
##   fig01_rate_mgL.png, fig01_rate_sat.png              cell significance (hatched = not significant)
##   fig01_rate_block_mgL.png, fig01_rate_block_sat.png  block significance (outlined = significant)
##   fig01_cell_trends.csv, fig01_block_trends.csv       trend tables (read by Figure6_volatility.R)
##   captions.csv                                         rows for each PNG (merge-written)
## Set OUTDIR_OVERRIDE / BLK_DEEP_OVERRIDE before sourcing to write elsewhere / change block size.
suppressMessages({library(data.table); library(ggplot2); library(scales); library(patchwork)
                  library(ggh4x); library(splines); library(quantreg); library(ggtext)})

meta <- data.table(
  lakeid=c("TR","BM","CR","SP","AL","TB","CB","ME","MO","FI","WI"),
  name  =c("Trout","Big Muskie","Crystal","Sparkling","Allequash","Trout Bog","Crystal Bog",
           "Mendota","Monona","Fish","Wingra"),
  region=c(rep("Northern",7), rep("Southern",4)),
  zmax  =c(35.7,21.3,20.4,20.0,8.0,7.9,2.5,25.3,22.5,18.9,4.0))
## Fish Lake (FI) and Lake Wingra (WI) are dropped from every figure.
meta <- meta[!lakeid %in% c("FI","WI")]

prof <- fread("data/profiles_clean.csv")[lakeid %in% meta$lakeid]
prof[, `:=`(date=as.Date(sampledate), doy=yday(as.Date(sampledate)), year=year4)]
prof <- prof[doy>=91 & doy<=319 & !(lakeid=="CR" & year %in% c(2012,2013))]  # Apr 1 - Nov 15

## ANOXIC FLOOR: DO below 0.5 mg/L is set to 0 (and o2sat to 0). Near-zero readings in anoxic water
## depend on sensor zero/offset and on how values were recorded (e.g. early NTL values in 0.1 mg/L
## steps), so a drift from ~0.3 to ~0.05 mg/L in water that is anoxic either way otherwise shows up as
## a highly significant (low-variance) "trend" -- e.g. Mendota 20-24 m in late July. Flooring makes
## anoxic water trendless; only movement across or above the threshold counts.
DO_FLOOR <- 0.5
n_floor <- prof[is.finite(o2) & o2 < DO_FLOOR, .N]
prof[is.finite(o2) & o2 < DO_FLOOR, `:=`(o2 = 0, o2sat = 0)]
prof[!is.finite(o2) & is.finite(o2sat) & o2sat < 5, o2sat := 0]   # ~0.5 mg/L equivalent when only %sat exists
cat(sprintf("anoxic floor: %d DO readings < %.1f mg/L set to 0\n", n_floor, DO_FLOOR))

## Per-region year range for the panel titles -- Southern lakes join the record later
## (1995/96) than Northern (1981), so these must NOT share a single hardcoded range.
ryrs <- merge(prof[, .(lakeid, year)], meta[, .(lakeid, region)], by = "lakeid")[
  , .(lo = min(year), hi = max(year)), by = region]
nyr <- ryrs[region == "Northern"]; syr <- ryrs[region == "Southern"]

## long over the variables, interpolate each cast to depth bins, then assign half-month time bins.
## Bin width is 0.5 m for lakes under 10 m (zmax) and 1 m for deeper lakes -- at 0.5 m throughout,
## the deep lakes' water columns rendered far denser than the pattern actually resolves.
zmax_lookup <- setNames(meta$zmax, meta$lakeid)
pm <- melt(prof[,.(lakeid,date,doy,year,depth,wtemp,o2,o2sat)],
           id.vars=c("lakeid","date","doy","year","depth"), variable.name="var", value.name="val")
pm <- pm[is.finite(val)]
g <- pm[, { o<-order(depth); dd<-depth[o]; yy<-val[o]
  dstep <- if(zmax_lookup[[lakeid]] < 10) 0.5 else 1
  if(length(dd)>=3 && diff(range(dd))>=1){ zg<-seq(0,max(dd),dstep); zg<-zg[zg>=min(dd)]
    .(dbin=zg, val=approx(dd,yy,zg,rule=2)$y) } else .(dbin=numeric(0),val=numeric(0)) },
  by=.(lakeid,date,doy,year,var)]
hedges <- c(91,106,121,136,152,167,182,197,213,228,244,259,274,289,305,320)  # Apr 1 - Nov 16
g[, bi:=findInterval(doy,hedges)]; g<-g[bi>=1 & bi<=length(hedges)-1]
g[, `:=`(tbin=(hedges[bi]+hedges[bi+1])/2, twidth=hedges[bi+1]-hedges[bi])]

## Seasonal detrending: within a half-month window the
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
## the natural non-parametric partner to Sen.
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

## Benjamini-Hochberg FDR correction, per lake (pooling all its depth x time-bin x variable
## cells into one family, computed once here so wtemp's q-values don't differ between the
## mg/L and %sat output figures). With ~14k Mann-Kendall tests across the whole figure,
## uncorrected p<0.05 alone would flag hundreds of cells by chance; Bonferroni would be overly
## conservative given how correlated neighboring cells are (same casts interpolated to adjacent
## depths, overlapping years between time bins) -- BH controls the expected false-discovery
## proportion instead, the standard fix for many correlated tests shown as a significance grid.
cell[!is.na(p), q := p.adjust(p, method="BH"), by=lakeid]

## ---- BLOCK significance (for the fig01_rate_block_*.png versions) ------------------------------
## Colours stay the per-cell trends; significance is instead tested on depth BLOCKS within each
## half-month: each year's block value is the mean of its depth cells' yearly residuals, then
## Sen + Mann-Kendall, BH FDR per lake across all its block tests (all variables), as for cells.
## Pooling adjacent depths gives fewer, less noisy tests, so coherent features that fail cell by
## cell (e.g. autumn surface warming) can pass. Block height scales with lake depth so shallow
## lakes aren't one or two giant boxes: BLK_DEEP m for zmax >= 10 m, 2 m for 5-10 m, 1 m below.
BLK_DEEP <- if(exists("BLK_DEEP_OVERRIDE")) BLK_DEEP_OVERRIDE else 5
blk_size <- function(z) fifelse(z >= 10, BLK_DEEP, fifelse(z >= 5, 2, 1))
## block index; a partial block at the bottom (fewer depth bins than a full block, e.g. Trout's
## 30-32 m below 25-29 m) is merged into the block above, so there are no sliver boxes at the bottom
lake_maxd <- cell[!is.na(slope), .(maxd=max(dbin)), by=lakeid][, setNames(maxd, lakeid)]   # deepest bin with a trend, per lake
blk_of <- function(dbin, lakeid){
  z <- unname(zmax_lookup[lakeid]); sz <- blk_size(z); st <- fifelse(z < 10, 0.5, 1)
  b <- floor(dbin / sz)
  md <- unname(lake_maxd[lakeid])
  bmax <- floor(md / sz)
  full <- sz / st                                     # depth bins in a full block
  nlast <- (md - bmax*sz) / st + 1                    # depth bins in the bottom block
  fifelse(b == bmax & bmax > 0 & nlast < full, bmax - 1, b)
}
yr[, blk := blk_of(dbin, lakeid)]
bkyr <- yr[, .(v=mean(v), nd=.N), by=.(lakeid,var,tbin,twidth,blk,year)]
bkyr[, ndmax := max(nd), by=.(lakeid,var,tbin,blk)]
bkyr <- bkyr[nd >= 0.5*ndmax]          # a year counts only if it has at least half the block's depths
bcell <- bkyr[, { n<-.N; sp<- if(n>0) diff(range(year)) else 0
  if(n>=MIN_YR_TREND && sp>=MIN_SPAN)
    .(bslope=sen(v,year)*10,
      bp=tryCatch(suppressWarnings(cor.test(year, v, method="kendall", exact=FALSE)$p.value), error=function(e) NA_real_))
  else .(bslope=NA_real_, bp=NA_real_) }, by=.(lakeid,var,tbin,twidth,blk)]
bcell[!is.na(bp), bq := p.adjust(bp, method="BH"), by=lakeid]
bcell[, bsig := !is.na(bq) & bq<0.05]

w <- cell[!is.na(slope)]
w[, sig := !is.na(q) & q<0.05]
w <- merge(w, meta, by="lakeid")
w[, dstep := fifelse(zmax<10, 0.5, 1)]  # matches the depth-bin width used to build dbin above

## Per-cell trend table: one source of truth for anything summarising this figure
## (Figure6_volatility.R reads it rather than re-fitting).
fwrite(w[, .(lakeid, name, region, zmax, var, tbin, twidth, dbin, dstep, slope, p, q, sig, nyr)],
       "figures/fig01_cell_trends.csv")
cat("wrote figures/fig01_cell_trends.csv\n")

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
## n = number of profiles (distinct sampling dates with temperature or DO) per lake in the Apr 1-Nov 15 window
ncast <- prof[is.finite(wtemp) | is.finite(o2) | is.finite(o2sat), .(ncast=uniqueN(date)), by=lakeid]
w <- merge(w, ncast, by="lakeid", all.x=TRUE)
w[, striplab := sprintf("<span style='font-size:7pt;font-weight:bold'>%s</span><br>%s<br><span style='font-size:5pt;color:#555555'>n = %s profiles</span>",
                        name, fifelse(is.na(lab_md), "", lab_md), format(ncast, big.mark=","))]
## three panel groups: northern lakes split by clarity trend (the bogs and Crystal are browning),
## then the southern lakes; within a group, lakes run by descending max depth
GROUPS <- c(clear="A  Northern Forested Lakes (%d–%d): No clarity change",
            brown="B  Northern Forested Lakes (%d–%d): Browning trend",
            south="C  Southern Ag/Urban Lakes (%d–%d)")
w[, grp := factor(fcase(lakeid %in% c("CR","TB","CB"), "brown", region == "Northern", "clear", default="south"),
                  levels=names(GROUPS))]
ord <- unique(w[, .(name, striplab, grp, zmax)])[order(grp,-zmax)]$striplab
w[, strip := factor(striplab, levels=ord)]

## block outlines in plot coordinates: span of the cells actually drawn in each block
bx <- w[, .(ymin=min(dbin)-dstep[1]/2, ymax=max(dbin)+dstep[1]/2),
        by=.(lakeid, var, tbin, twidth, blk=blk_of(dbin, lakeid), strip, region, grp, name, zmax)]
bx <- merge(bx, bcell[, .(lakeid, var, tbin, blk, bslope, bq, bsig)], by=c("lakeid","var","tbin","blk"))
fwrite(bx[, .(lakeid, name, var, tbin, depth_top=ymin, depth_bottom=ymax, slope=bslope, q=bq, sig=bsig)],
       "figures/fig01_block_trends.csv")

## Non-significant cells get a light diagonal stroke (one per tile). This is deliberately a plain
## geom_segment rather than ggpattern: ggpattern rasterizes every tile separately and takes many
## minutes at these tile counts (~14k cells), whereas segments stay vector and render instantly.
hatch <- function(d) d[sig==FALSE, .(x=tbin-twidth/2, xend=tbin+twidth/2, y=dbin-dstep/2, yend=dbin+dstep/2, strip)]

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
  panel.grid=element_blank(), panel.spacing=unit(0.25,"lines"), plot.margin=margin(1, 3, 1, 3),
  legend.position="bottom", legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.22,"cm"),
  legend.title=element_text(size=6.5), plot.title=element_text(face="bold",size=9.5),
  plot.subtitle=element_text(face="bold", size=6.5, color="grey30"),
  axis.text=element_text(size=5), axis.text.x=element_text(size=6),
  axis.title.y=element_text(size=7.5),
  strip.text.y.right=element_markdown(angle=0, hjust=0, size=5, lineheight=1.2))

block <- function(gg, vv, pal, right, ttl, sub, mode="cell"){
  d <- w[grp==gg & var==vv]
  ggplot(d, aes(tbin,dbin,fill=slope)) +
    geom_tile(aes(width=twidth, height=dstep)) +
    { if(mode=="cell") geom_segment(data=hatch(d), aes(x=x,xend=xend,y=y,yend=yend), inherit.aes=FALSE,
                                    color="grey25", linewidth=0.16, alpha=0.5)
      else geom_rect(data=bx[grp==gg & var==vv & bsig==TRUE],
                     aes(xmin=tbin-twidth/2, xmax=tbin+twidth/2, ymin=ymin, ymax=ymax),
                     inherit.aes=FALSE, fill=NA, color="black", linewidth=0.35) } +
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

make_fig <- function(dovar, dopal, dolab, subtitle, outfile, mode="cell"){
  yrs_of <- function(gg) if(gg == "south") syr else nyr
  ## each group is a row: its title on a strip spanning both columns, then temperature (left) and
  ## DO (right) panels; the column subtitles only on the top group
  ttl_row <- function(txt) wrap_elements(full=grid::textGrob(txt, x=unit(0.02, "npc"), y=unit(0.2, "npc"), hjust=0, vjust=0,
                 gp=grid::gpar(fontface="bold", fontsize=9.5))) & theme(plot.margin=margin(0, 0, 0, 0))
  rows <- list(); hts <- c()
  for(gg in names(GROUPS)){ y <- yrs_of(gg); top <- gg == names(GROUPS)[1]
    rows <- c(rows, list(ttl_row(sprintf(GROUPS[[gg]], y$lo, y$hi)),
      block(gg, "wtemp", pal_T, FALSE, NULL, if(top) "Trend in temperature" else NULL, mode),
      block(gg, dovar,   dopal, TRUE,  NULL, if(top) dolab else NULL, mode)))
    hts <- c(hts, 0.32, w[grp == gg, uniqueN(lakeid)] + if(top) 0.3 else 0) }  # panel height ~ lakes
  lay <- paste(sprintf("%s%s\n%s%s", LETTERS[seq(1, 9, 3)], LETTERS[seq(1, 9, 3)],
                       LETTERS[seq(2, 9, 3)], LETTERS[seq(3, 9, 3)]), collapse="\n")
  design <- wrap_plots(rows, design=lay) + plot_layout(heights=hts, guides="collect") &
    theme(legend.position="bottom")
  final <- wrap_elements(design) + theme(plot.margin=margin(2, 2, 2, 2))
  ggsave(outfile, final, width=6.5, height=8, dpi=500, bg="white"); cat("wrote", outfile, "\n")
  captions[[outfile]] <<- data.table(file=outfile, title=FIG_TITLE, caption=subtitle)
}

yrs <- range(prof$year)
base <- sprintf("April–mid-November. Theil-Sen slope per DECADE over the full %d–%d record (not a two-era difference): every year contributes, so a steady ramp and a step change are no longer confounded.\nCells are hatched where the trend is NOT significant after Benjamini-Hochberg FDR correction (per lake, across all its depth/time/variable cells; Mann-Kendall q<0.05); plain cells are significant. Blank = fewer than %d years or <%d-year span, too little to fit a trend.\nValues are trends in the seasonally-detrended residual, so drift in which days of a fortnight got sampled cannot masquerade as a trend. DO below 0.5 mg/L is set to 0 (anoxic) so that drift in near-zero sensor readings cannot either.\nBeside each lake name: n = number of profiles used, and full-record trend in surface DOC and in Secchi depth over the same window — brown ↑ rising, blue ↓ falling, grey ✕ no significant trend (Mann-Kendall p<0.05).",
   yrs[1], yrs[2], MIN_YR_TREND, MIN_SPAN)
make_fig("o2",    pal_O,    "Trend in dissolved oxygen (mg/L)",
  paste0(base, "\nOxygen trends are strongly skewed negative — deep water is losing O₂ across nearly every lake, while surface/metalimnetic waters hold steadier."), "figures/fig01_rate_mgL.png")
make_fig("o2sat", pal_Osat, "Trend in dissolved oxygen (% sat)",
  paste0(base, "\nAs % saturation: warming lowers O₂ solubility, so the same mg/L loss reads as a larger saturation deficit in deep water."), "figures/fig01_rate_sat.png")

## block-significance versions: same colours, boxes around significant depth blocks instead of hatching
OUTDIR <- if(exists("OUTDIR_OVERRIDE")) OUTDIR_OVERRIDE else "figures"
bbase <- sprintf(paste0("April–mid-November. Colours: Theil-Sen slope per DECADE over the full %d–%d record for each depth x half-month cell ",
  "(trend in the seasonally-detrended residual; DO < 0.5 mg/L set to 0 as anoxic). Significance is tested on depth BLOCKS within each half-month -- %d m in lakes deeper than ",
  "10 m, 2 m in lakes of 5-10 m, 1 m in Crystal Bog -- using each year's block-mean residual (Mann-Kendall, Benjamini-Hochberg q<0.05 per lake ",
  "across all its block tests); significant blocks are outlined in black. Blank = fewer than %d years or <%d-year span.\n",
  "Beside each lake name: n = number of profiles used, and full-record trend in surface DOC and in Secchi depth — brown ↑ rising, blue ↓ falling, grey ✕ no significant trend (Mann-Kendall p<0.05)."),
  yrs[1], yrs[2], BLK_DEEP, MIN_YR_TREND, MIN_SPAN)
make_fig("o2",    pal_O,    "Trend in dissolved oxygen (mg/L)", bbase, file.path(OUTDIR, "fig01_rate_block_mgL.png"), mode="block")
make_fig("o2sat", pal_Osat, "Trend in dissolved oxygen (% sat)", paste0(bbase, "\nAs % saturation: warming lowers O₂ solubility, so the same mg/L loss reads as a larger saturation deficit in deep water."),
         file.path(OUTDIR, "fig01_rate_block_sat.png"), mode="block")
cat("\n=== % of drawn area in significant blocks, by lake & variable ===\n")
print(dcast(bx[, .(pct=round(100*mean(bsig))), by=.(name, var)], name ~ var, value.var="pct"))

if(OUTDIR == "figures") write_captions(rbindlist(captions))

## Console summary. NOTE: deliberately summarises ALL cells, not just the significant ones.
## Taking the median over only q<0.05 cells is selection-biased (winner's curse) -- conditioning on
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
