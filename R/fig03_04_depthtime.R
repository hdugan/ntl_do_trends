## Figs 8/9 with per-tile change outlines + %-saturation versions.
## Outline a tile YELLOW if DO significantly DECREASED vs the previous era, RED if it INCREASED (t-test p<0.05).
suppressMessages({library(data.table); library(ggplot2); library(scales)})

meta <- data.table(
  lakeid=c("TR","BM","CR","SP","AL","TB","CB","ME","MO","FI","WI"),
  name  =c("Trout","Big Muskellunge","Crystal","Sparkling","Allequash","Trout Bog","Crystal Bog",
           "Mendota","Monona","Fish","Wingra"),
  region=c(rep("Northern",7), rep("Southern",4)),
  trophic=c("oligotrophic","oligotrophic","oligotrophic","oligotrophic","mesotrophic","dystrophic","dystrophic",
            "eutrophic","eutrophic","mesotrophic","eutrophic"),
  zmax  =c(35.7,21.3,20.4,20.0,8.0,7.9,2.5,25.3,22.5,18.9,4.0))
tro <- c(oligotrophic="oligotrophic",mesotrophic="mesotrophic",dystrophic="dystrophic (bog)",eutrophic="eutrophic")

prof <- fread("data/profiles_clean.csv")[lakeid %in% meta$lakeid]
prof[, `:=`(date=as.Date(sampledate), doy=yday(as.Date(sampledate)), year=year4)]
prof <- prof[doy>=105 & doy<=310 & !(lakeid=="CR" & year %in% c(2012,2013))]  # drop Crystal mixing yrs
prof <- merge(prof, meta, by="lakeid")
prof[, strip := factor(sprintf("%s\n%s · %.0f m", name, tro[trophic], zmax))]

## palettes
MIN_YEARS <- 4   # per-era cells need data from at least this many distinct years to be averaged
pal_do  <- scale_fill_gradientn(colours=c("#000000","#7f0000","#d7301f","#fdae61","#ffffbf","#66bd63","#1a9850"),
             values=rescale(c(0,1,2,4,6,9,14)), name="DO (mg/L)", limits=c(0,14), oob=squish, na.value="grey80")
pal_sat <- scale_fill_gradientn(colours=c("#000000","#7f0000","#d7301f","#fdae61","#ffffbf","#66bd63","#1a9850"),
             values=rescale(c(0,10,20,40,70,95,140)), name="DO (% sat)", limits=c(0,140), oob=squish, na.value="grey80")
mbound <- c(91,121,152,182,213,244,274,305)                 # 1st of Apr…Nov (day-of-year)
hedges <- c(91,106,121,136,152,167,182,197,213,228,244,259,274,289,305)  # half-month edges (1st & 16th)
xsc <- scale_x_continuous(breaks=c(106,136.5,167,197.5,228.5,259,289.5),   # month centres
                          labels=c("Apr","May","Jun","Jul","Aug","Sep","Oct"),
                          minor_breaks=NULL, expand=c(0,0))
outline <- scale_color_manual(values=c(`DO decreased`="#ffe500", `DO increased`="#ff2a2a"),
             name="2016–25 vs\nearlier baseline\n(t-test p<0.05)", na.translate=FALSE)

build <- function(vv, reg, breaks, labels, title, sub, pal, outfile, W, H, min_years=MIN_YEARS){
  d <- prof[region==reg & is.finite(get(vv))]
  d <- d[, .(val=mean(get(vv))), by=.(lakeid, date, doy, year, strip, zmax, depth)]
  # vertical interpolation to 0.5 m grid
  g <- d[, { o<-order(depth); dd<-depth[o]; yy<-val[o]
    if(length(dd)>=3 && diff(range(dd))>=1){ zg<-seq(0,max(dd),0.5); zg<-zg[zg>=min(dd)]
      .(dbin=zg, val=approx(dd, yy, zg, rule=2)$y) } else .(dbin=numeric(0), val=numeric(0)) },
    by=.(lakeid, date, doy, year, strip, zmax)]
  g[, bi := findInterval(doy, hedges)]                       # half-month bins aligned to month edges
  g <- g[bi>=1 & bi<=length(hedges)-1]
  g[, `:=`(tbin=(hedges[bi]+hedges[bi+1])/2, twidth=hedges[bi+1]-hedges[bi])]
  g[, eidx := cut(year, breaks=breaks, labels=FALSE, right=TRUE)]
  g <- g[!is.na(eidx)]; g[, era := factor(labels[eidx], levels=labels)]

  ## per-era cell means (fill); cells with too few contributing years are blanked out (grey),
  ## unless that depth/time window never reaches min_years in ANY era, in which case it's dropped (white)
  agg <- g[, .(m=mean(val), n_yr=uniqueN(year)), by=.(strip, zmax, eidx, era, tbin, twidth, dbin)]
  agg[, ever_ok := any(n_yr>=min_years), by=.(strip, zmax, tbin, twidth, dbin)]
  agg <- agg[ever_ok==TRUE][, ever_ok := NULL]
  agg[n_yr < min_years, m := NA]
  ## a cell is eligible for significance testing only if every era with data there has >=min_years
  sig_elig <- agg[, .(sig_ok=all(n_yr>=min_years)), by=.(strip, zmax, tbin, twidth, dbin)]

  ## significance: FINAL era vs ALL earlier eras pooled; outline the final era only
  maxidx <- max(g$eidx)
  sig <- g[, {
    fin <- val[eidx==maxidx]; ear <- val[eidx<maxidx]
    if (length(fin)>=3 && length(ear)>=3) {
      pv <- tryCatch(t.test(fin, ear)$p.value, error=function(e) NA_real_); dd <- mean(fin)-mean(ear)
      .(chg=fifelse(!is.na(pv)&pv<0.05 & dd<0, "DO decreased",
             fifelse(!is.na(pv)&pv<0.05 & dd>0, "DO increased", NA_character_)))
    } else .(chg=NA_character_)
  }, by=.(strip, zmax, tbin, twidth, dbin)]
  sig <- merge(sig, sig_elig, by=c("strip","zmax","tbin","twidth","dbin"), all.x=TRUE)
  sig[is.na(sig_ok) | !sig_ok, chg := NA_character_]
  sig <- sig[!is.na(chg)][, sig_ok := NULL]; sig[, era := factor(labels[maxidx], levels=labels)]

  ord <- meta[region==reg][order(-zmax)][, sprintf("%s\n%s · %.0f m", name, tro[trophic], zmax)]
  agg[, strip := factor(strip, levels=ord)]; sig[, strip := factor(strip, levels=ord)]
  ggp <- ggplot(agg, aes(tbin, dbin)) +
    geom_tile(aes(fill=m, width=twidth), height=0.5) +
    geom_tile(data=sig, aes(tbin, dbin, color=chg, width=twidth), fill=NA, linewidth=0.32, height=0.5, inherit.aes=FALSE) +
    geom_vline(xintercept=mbound, color="white", linewidth=0.35, alpha=0.85) +   # month boundaries
    facet_grid(strip ~ era, scales="free_y") +
    scale_y_reverse(expand=c(0,0)) + xsc + pal + outline +
    guides(color=guide_legend(override.aes=list(fill="grey85", linewidth=1))) +
    labs(title=title, subtitle=sub, x=NULL, y="Depth (m)") +
    theme_minimal(base_size=12) + theme(panel.grid=element_blank(),
      plot.title=element_text(face="bold",size=14), plot.subtitle=element_text(color="grey35",size=10),
      strip.text.y=element_text(angle=0, hjust=0, size=10, lineheight=0.9), panel.spacing=unit(0.5,"lines"))
  ggsave(outfile, ggp, width=W, height=H, dpi=500, bg="white")
  cat("wrote", outfile, "\n")
}

subN <- function(my) sprintf("Outlines on 2016–25 only: yellow = DO significantly lower than the pooled 1986–2015 baseline, red = higher (t-test p<0.05), tested only where every era has ≥%d years of data. Crystal 2012–13 excluded.\nGrey = <%d years of data this era, though other eras were adequately sampled; blank = that window was never adequately sampled.", my, my)
subS <- function(my) sprintf("Outlines on 2016–25 only: yellow = DO significantly lower than the pooled 1996–2015 baseline, red = higher (t-test p<0.05), tested only where every era has ≥%d years of data. Crystal 2012–13 excluded.\nGrey = <%d years of data this era, though other eras were adequately sampled; blank = that window was never adequately sampled.", my, my)
nb <- c(1985,1995,2005,2015,2025); nl <- c("1986–95","1996–05","2006–15","2016–25")
sb <- c(1995,2005,2015,2025);      sl <- c("1996–05","2006–15","2016–25")

build("o2",    "Northern", nb, nl, "Northern NTL lakes: dissolved oxygen (mg/L) by depth & season", subN(MIN_YEARS), pal_do,  "figures/fig03_northern_do.png",  15, 13)
build("o2sat", "Northern", nb, nl, "Northern NTL lakes: oxygen saturation (%) by depth & season",    subN(MIN_YEARS), pal_sat, "figures/fig03_northern_sat.png", 15, 13)
build("o2",    "Southern", sb, sl, "Southern NTL lakes: dissolved oxygen (mg/L) by depth & season",   subS(MIN_YEARS), pal_do,  "figures/fig04_southern_do.png",  12, 8.5)
build("o2sat", "Southern", sb, sl, "Southern NTL lakes: oxygen saturation (%) by depth & season",     subS(MIN_YEARS), pal_sat, "figures/fig04_southern_sat.png", 12, 8.5)
