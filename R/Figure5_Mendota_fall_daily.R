## Figure5_Mendota_fall_daily.R: Lake Mendota autumn (1 Sep-9 Nov) temperature trends by depth x DAY.
## A zoom on Figure 2 (Figure2_rate_of_change.R) for the period when warming is strongest, using the
## Mendota buoy (EDI knb-lter-ntl.130, data/Buoy/ntl130_1_v14.csv, 2006-2025; plus 2026 from
## data/Buoy/ntl130_daily_2026.csv, built by R/01c_buoy_daily_2026.R): daily means, 0-20 m.
##
## Per year x 1-m depth: daily buoy means (flagged values removed); small depth gaps in a day's profile
## filled by linear interpolation (needs >= 10 of 21 depths); then a 7-day CENTRED running mean (>= 4 of
## 7 days), which removes single-day weather without the staircase of fixed weekly bins.
## Colour = Theil-Sen trend (degC/decade) of that smoothed value across years, per day x depth cell.
##
## NEAR-COMMON YEARS: only years with a smoothed value on >= 75% of the season's days are used, so
## (nearly) the same autumns underlie every day. Otherwise the set of years changes through the season
## (buoy pulled mid-Oct in 2006/2007/2013; partial years like 2010) and that composition change -- not
## the lake -- creates spurious day-to-day differences.
##
## SIGNIFICANCE is tested on blocks of 5 m (0-4, 5-9, 10-14, 15-20 m) x 2 weeks: each year's block value
## is the mean of its daily 1-m values, then Theil-Sen + Mann-Kendall, Benjamini-Hochberg FDR across the
## 20 blocks. Only significant blocks are outlined. (With ~14 autumns, year-to-year autumn weather
## dominates and no block is currently significant -- the figure is descriptive; significance for
## Mendota's autumn warming rests on the 31-year record in Figure 2.)
## Panel B: fall turnover date in each buoy year (end of stratification: bottom 18-20 m minus surface 0-2 m
## density < 0.3 kg/m3, sustained 30 days), read from figures/figSX_Mendota_turnover.csv -- run
## FigureSX_Mendota_turnover.R first. Theil-Sen trend + Mann-Kendall over the buoy years.
suppressMessages({library(data.table); library(ggplot2); library(scales); library(patchwork)})

ZMAX     <- 20
NDAY     <- 70      # pday 1 = 1 Sep ... 70 = 9 Nov (counted from 1 Sep, so leap years don't shift it)
MIN_COV  <- 0.75    # near-common-year rule: smoothed value on >= 75% of days 4..67
MIN_YR <- 10; MIN_SPAN <- 15
LIM      <- 2       # colour limits, degC/decade

## ---- buoy: daily -> 1-m bins -> gap-filled profiles -> 7-day running mean -------------------------
## 2006-2025 from EDI (ntl130); 2026 from the 1-minute buoy record, made daily by R/01c_buoy_daily_2026.R
b <- rbind(fread("data/Buoy/ntl130_1_v14.csv", colClasses=c(flag_wtemp="character")),
           fread("data/Buoy/ntl130_daily_2026.csv", colClasses=c(flag_wtemp="character")))[(is.na(flag_wtemp) | flag_wtemp=="") & is.finite(wtemp) & depth<=ZMAX]
b[, date := as.IDate(sampledate)]
b[, `:=`(year=year(date), pday=as.integer(date - as.IDate(sprintf("%d-09-01", year(date)))) + 1L, dbin=round(depth))]
b <- b[pday>=1 & pday<=NDAY]
d <- b[, .(t=mean(wtemp)), by=.(year, pday, dbin)]
d <- d[, { if(uniqueN(dbin) >= 10){ o <- order(dbin); .(dbin=0:ZMAX, t=approx(dbin[o], t[o], 0:ZMAX, rule=2)$y) }
           else .(dbin=integer(0), t=numeric(0)) }, by=.(year, pday)]
d <- d[CJ(year=unique(d$year), dbin=0:ZMAX, pday=1:NDAY), on=.(year, dbin, pday)]
setorder(d, year, dbin, pday)
d[, t7 := frollapply(t, 7, function(v) if(sum(!is.na(v)) >= 4) mean(v, na.rm=TRUE) else NA_real_, align="center"), by=.(year, dbin)]

## near-common years
inner <- 4:(NDAY-3)
cov <- d[pday %in% inner & dbin==0, .(cov=mean(is.finite(t7))), by=year]
keep <- cov[cov >= MIN_COV, year]
d <- d[year %in% keep & is.finite(t7)]
cat("years kept:", paste(sort(keep), collapse=","), "\n")
cat("years dropped (coverage):", cov[!year %in% keep, paste(sprintf("%d(%.0f%%)", year, 100*cov), collapse=" ")], "\n")

## ---- trends per day x depth ------------------------------------------------------------------------
sen <- function(y, x){ n <- length(y); dd <- outer(y,y,"-")[lower.tri(diag(n))] / outer(x,x,"-")[lower.tri(diag(n))]; median(dd[is.finite(dd)]) }
tr <- d[, if(.N >= MIN_YR && diff(range(year)) >= MIN_SPAN) .(slope=10*sen(t7, year), nyr=.N), by=.(pday, dbin)]
fwrite(tr, "figures/fig05_Mendota_fall_daily_trends.csv")

## ---- 5 m x 2-week block significance ---------------------------------------------------------------
BLK <- c(0, 5, 10, 15, ZMAX + 1)
d[, `:=`(blk=findInterval(dbin, BLK), per=(pday - 1L) %/% 14L + 1L)]          # per 1..5
bk <- d[, .(t=mean(t7)), by=.(year, per, blk)][, .(slope=10*sen(t, year),
        p=suppressWarnings(cor.test(year, t, method="kendall", exact=FALSE)$p.value), nyr=.N), by=.(per, blk)]
bk[, q := p.adjust(p, "BH")][, sig := q < 0.05]
bk[, `:=`(x0=(per-1)*14 + 0.5, x1=pmin(per*14, NDAY) + 0.5, ytop=BLK[blk]-0.5, ybot=BLK[blk+1]-0.5,
          depth=sprintf("%d-%d m", BLK[blk], BLK[blk+1]-1))]
fwrite(bk[, .(per, days=sprintf("%d-%d", (per-1)*14+1, pmin(per*14, NDAY)), depth, slope=round(slope,2),
              p=signif(p,3), q=signif(q,3), sig, nyr)], "figures/fig05_Mendota_fall_daily_blocks.csv")
cat("\nblock trends (degC/decade; * = BH q<0.05):\n")
print(dcast(bk[, .(depth, per, s=sprintf("%+.2f%s (p=%.3f)", slope, fifelse(sig,"*",""), p))], depth ~ per, value.var="s"))

## ---- plot ----------------------------------------------------------------------------------------
ga <- ggplot(tr, aes(pday, dbin, fill=slope)) +
  geom_raster() +
  geom_vline(xintercept=c(30.5, 61.5), color="white", linewidth=0.4) +           # 1 Oct, 1 Nov
  geom_rect(data=bk[sig==TRUE], aes(xmin=x0, xmax=x1, ymin=ytop, ymax=ybot), inherit.aes=FALSE,
            fill=NA, color="black", linewidth=0.7) +
  scale_fill_gradientn(colours=c("#2166ac","#4393c3","#92c5de","#d1e5f0","#f7f7f7","#fddbc7","#f4a582","#d6604d","#b2182b"),
                       limits=c(-LIM, LIM), oob=squish, breaks=c(-2,-1,0,1,2), name="Δ Temp (°C/decade)") +
  guides(fill=guide_colourbar(title.position="top", title.hjust=0.5)) +
  scale_x_continuous(breaks=c(1,8,15,22,31,38,45,52,62,69), labels=c("1 Sep","8","15","22","1 Oct","8","15","22","1 Nov","8"),
                     limits=c(0.5, NDAY+0.5), expand=c(0,0)) +
  scale_y_reverse(expand=c(0,0), breaks=seq(0, 20, 5)) +
  labs(x=NULL, y="Depth (m)", tag="A",
       title=sprintf("Temperature trend by depth and day, %d–%d", min(keep), max(keep))) +
  theme_minimal(base_size=7) +
  theme(panel.grid=element_blank(), plot.title=element_text(face="bold", size=7), legend.position="bottom",
        legend.key.width=unit(0.6,"cm"), legend.key.height=unit(0.22,"cm"), plot.tag=element_text(face="bold", size=10))

## ---- panel B: buoy turnover dates ----------------------------------------------------------------
tpath <- "figures/figSX_Mendota_turnover.csv"
if(!file.exists(tpath)) stop(tpath, " not found -- run FigureSX_Mendota_turnover.R first")
to <- fread(tpath)[is.finite(turn_buoy), .(year, turn=turn_buoy)]
to_sen <- 10 * sen(to$turn, to$year)
to_p <- suppressWarnings(cor.test(to$year, to$turn, method="kendall", exact=FALSE)$p.value)
to_int <- median(to$turn - to_sen/10 * to$year)
cat(sprintf("buoy turnover %d-%d (n=%d): %+.1f days/decade, Mann-Kendall p=%.3f\n",
            min(to$year), max(to$year), nrow(to), to_sen, to_p))
gb <- ggplot(to, aes(year, turn)) +
  annotate("segment", x=min(to$year), xend=max(to$year), y=to_int + to_sen/10*min(to$year),
           yend=to_int + to_sen/10*max(to$year), linewidth=0.5, color="grey30") +
  geom_point(size=1.6, color="#b2182b", shape=15) +
  annotate("text", x=min(to$year), y=Inf, hjust=0, vjust=1.3, size=2.2, color="grey20",
           label=sprintf("%+.1f days/decade (p = %.2f)", to_sen, to_p)) +
  scale_y_continuous(breaks=c(274, 281, 288, 295, 305, 312), labels=c("1 Oct","8","15","22","1 Nov","8"),
                     limits=c(274, 314)) +
  scale_x_continuous(breaks=seq(2010, 2025, 5)) +
  labs(x=NULL, y=NULL, tag="B", title="Fall turnover date") +
  theme_minimal(base_size=7) +
  theme(panel.grid.minor=element_blank(), plot.title=element_text(face="bold", size=7),
        plot.tag=element_text(face="bold", size=10))
g <- ga + gb + plot_layout(widths=c(2.6, 1))
ggsave("figures/fig05_Mendota_fall_daily.png", g, width=6.5, height=3.4, dpi=500, bg="white")
cat("wrote figures/fig05_Mendota_fall_daily.png\n")

write_captions <- function(new_caps){
  path <- "figures/captions.csv"
  old <- if(file.exists(path)) fread(path) else data.table(file=character(),title=character(),caption=character())
  fwrite(rbind(old[!file %in% new_caps$file], new_caps), path); cat("wrote", path, "\n")
}
write_captions(data.table(
  file="figures/fig05_Mendota_fall_daily.png",
  title="Lake Mendota autumn warming and fall turnover from the buoy record",
  caption=paste0(
    "(A) Theil-Sen trend in water temperature (°C per decade) for each day and 1-m depth, 1 September-9 November, ",
    "from the Lake Mendota buoy (EDI knb-lter-ntl.130; flagged values removed). Daily profiles were gap-filled by linear ",
    "interpolation over depth and smoothed with a centred 7-day running mean before fitting trends. Only autumns with ",
    "data on at least 75% of days are used (", length(keep), " years, ", min(keep), "-", max(keep), "), so the same years ",
    "underlie every day. Significance was tested on blocks of 5 m by two weeks using the block-mean temperature each year ",
    "(Mann-Kendall, Benjamini-Hochberg q < 0.05 across 20 blocks); significant blocks would be outlined in black, and none ",
    "was significant. The sloping boundary between warming (above) and cooling (below) follows the base of the deepening ",
    "autumn mixed layer. (B) Date of fall turnover in each buoy year, defined as the day the density difference ",
    "between bottom (18-20 m) and surface (0-2 m) water falls below 0.3 kg m-3 and stays below it for 30 days ",
    sprintf("(%d years, %d-%d; years when the buoy was removed before turnover are omitted). Line: Theil-Sen trend ",
            nrow(to), min(to$year), max(to$year)),
    sprintf("(%+.1f days per decade, Mann-Kendall p = %.2f).", to_sen, to_p))))
