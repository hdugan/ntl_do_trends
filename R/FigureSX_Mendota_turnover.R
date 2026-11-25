## FigureSX_Mendota_turnover.R: when does "fall" start in Lake Mendota? Fall turnover (the end of summer
## stratification) dated each year, and how it has shifted.
## Turnover = the day the density difference between bottom water (mean of 18-20 m) and surface water
## (mean of 0-2 m) falls below 0.3 kg/m3 after mid-July and stays below it for the next 30 days (the
## stratification-end rule of Jane et al. 2021 / data/Stetler/stratification_duration.R, with fixed layers
## so all sources measure the same thing). Density from temperature (rLakeAnalyzer::water.density);
## profiles are linearly interpolated to 0.5 m first.
## Sources:
##   buoy      data/Buoy/ntl130_1_v14.csv + ntl130_daily_2026.csv (R/01c_buoy_daily_2026.R), daily means, 2007-2026 (2006 only reaches 15 m; 2014-2015 and
##             some years are pulled before turnover)
##   profiles  data/profiles_clean.csv (NTL-LTER + 2006-2019 MEMO casts), 1995-2025; the crossing is linearly
##             interpolated between the last stratified cast and the first mixed cast, whose dates bound it
##   1906      data/Juday/Mendota_1906_temp.csv (Birge & Juday plate, digitized; casts every 2-10 days)
## Best estimate per year = buoy date where available, else profile date.
## Trend: Theil-Sen slope + Mann-Kendall, 1995-2025.
##   (A) turnover date by year, 1906 and 1995-2025
##   (B) bottom-minus-surface density difference through late summer and fall: each buoy year (grey),
##       1906 (black), threshold dashed
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(rLakeAnalyzer)})
THR <- 0.3        # kg/m3
START <- 196      # day of year (mid-July) after which a crossing counts as fall turnover

## bottom (18-20 m) minus surface (0-2 m) density from one profile
ddiff <- function(z, t){
  ok <- is.finite(t); z <- z[ok]; t <- t[ok]
  if(length(unique(z)) < 3 || min(z) > 1 || max(z) < 18) return(NA_real_)
  zg <- seq(0, 20, 0.5); tt <- approx(z, t, zg, rule=2, ties=mean)$y
  mean(water.density(tt[zg >= 18])) - mean(water.density(tt[zg <= 2])) }
## first crossing below THR after START that follows a stratified observation and stays below for 30 days;
## returns the interpolated crossing and the dates of the two observations that bracket it
turnover <- function(doy, d){
  o <- order(doy); doy <- doy[o]; d <- d[o]; k <- doy >= START; doy <- doy[k]; d <- d[k]
  na <- list(turn=NA_real_, lo=NA_real_, hi=NA_real_)
  if(length(d) < 2 || d[1] < THR) return(na)
  for(i in 2:length(d)) if(d[i] < THR && d[i-1] >= THR && all(d[doy > doy[i] & doy <= doy[i] + 30] < THR))
    return(list(turn=as.numeric(doy[i-1] + (doy[i] - doy[i-1]) * (d[i-1] - THR)/(d[i-1] - d[i])), lo=as.numeric(doy[i-1]), hi=as.numeric(doy[i])))
  na }

## ---- profiles ----
p <- fread("data/profiles_clean.csv")[lakeid == "ME" & is.finite(wtemp)]
p[, doy := yday(as.IDate(sampledate))]
pd <- p[, .(d=ddiff(depth, wtemp)), by=.(year=year4, doy)][is.finite(d)]
tp <- pd[, turnover(doy, d), by=year]
## ---- buoy (daily) ----
b <- rbind(fread("data/Buoy/ntl130_1_v14.csv", colClasses=c(flag_wtemp="character")),
           fread("data/Buoy/ntl130_daily_2026.csv", colClasses=c(flag_wtemp="character")))[is.na(flag_wtemp) | flag_wtemp == ""]
b[, doy := yday(as.IDate(sampledate))]
bd <- b[, .(d=ddiff(depth, wtemp)), by=.(year=year4, doy)][is.finite(d)]
tb <- bd[, turnover(doy, d), by=year][is.finite(turn), .(year, turn_buoy=turn)]
## ---- 1906 ----
j <- fread("data/Juday/Mendota_1906_temp.csv")[is.finite(temp_C)]
jd <- j[, .(d=ddiff(depth_m, temp_C)), by=.(doy=yday)][is.finite(d)]
t06 <- jd[, turnover(doy, d)]

tt <- merge(tp, tb, by="year", all=TRUE)
tt[, best := fifelse(is.finite(turn_buoy), turn_buoy, turn)]
tt[, src := fifelse(is.finite(turn_buoy), "Buoy (daily)", "Profiles (bracketed by casts)")]
fwrite(tt, "figures/figSX_Mendota_turnover.csv")
LAST <- max(tt[is.finite(best), year])                   # last year with a turnover date
BUOY_YRS <- sprintf("%d–%d", min(tb$year), max(tb$year))

sen <- function(x, y){ ok <- is.finite(y); x <- x[ok]; y <- y[ok]
  s <- outer(y, y, "-") / outer(x, x, "-"); s <- s[upper.tri(s)]; m <- median(s[is.finite(s)])
  list(slope=m, int=median(y - m*x), p=Kendall::MannKendall(y)$sl[1], n=length(y)) }
S_best <- sen(tt$year, tt$best); S_prof <- sen(tt$year, tt$turn)
cat(sprintf("turnover trend 1995-%d: best %.1f d/decade (p=%.3f, n=%d); profiles only %.1f d/decade (p=%.3f)\n",
            LAST, 10*S_best$slope, S_best$p, S_best$n, 10*S_prof$slope, S_prof$p))
cat(sprintf("1906 turnover: day %.1f (%s), bracketed by casts on days %.0f-%.0f\n", t06$turn,
            format(as.Date(t06$turn - 1, origin="1906-01-01"), "%d %b"), t06$lo, t06$hi))
med <- tt[, .(med=median(best, na.rm=TRUE)), by=.(decade=10*(year %/% 10))]; print(med)

## ---- plots ----
dlab <- function(x) format(as.Date(x - 1, origin="2001-01-01"), "%b %d")   # day of year -> date (non-leap)
YL <- c(270, 316)
th <- theme_minimal(base_size=8) + theme(panel.grid.minor=element_blank(),
  plot.tag=element_text(face="bold", size=10), legend.position="bottom", legend.title=element_blank(),
  legend.text=element_text(size=7))
pts <- rbind(tt[is.finite(best), .(year, best, lo=fifelse(is.finite(turn_buoy), NA_real_, lo),
                                    hi=fifelse(is.finite(turn_buoy), NA_real_, hi), src)],
             data.table(year=1906, best=t06$turn, lo=t06$lo, hi=t06$hi, src="1906 (Birge & Juday)"))
pts[, src := factor(src, levels=c("1906 (Birge & Juday)", "Profiles (bracketed by casts)", "Buoy (daily)"))]
## break the year axis between 1906 and 1995 by plotting 1906 at a stand-in x position
X06 <- 1990; pts[, x := fifelse(year == 1906, X06, year)]
pa <- ggplot(pts, aes(x, best, color=src, shape=src)) +
  annotate("rect", xmin=1992.5, xmax=1993.5, ymin=-Inf, ymax=Inf, fill="white") +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=0, linewidth=0.3, alpha=0.35, na.rm=TRUE) +
  annotate("segment", x=1995, xend=LAST, y=S_best$int + S_best$slope*1995, yend=S_best$int + S_best$slope*LAST,
           linewidth=0.5, color="grey30") +
  geom_point(size=1.6) +
  annotate("text", x=LAST, y=YL[2], hjust=1, vjust=1, size=2.4, color="grey20",
           label=sprintf("1995–%d: +%.1f days/decade (p = %.2f)", LAST, 10*S_best$slope, S_best$p)) +
  scale_x_continuous(breaks=c(X06, seq(1995, 2025, 10)), labels=c("1906", seq(1995, 2025, 10))) +
  scale_y_continuous(labels=dlab, breaks=seq(270, 310, 10)) + coord_cartesian(ylim=YL) +
  scale_color_manual(values=c("black", "#4393c3", "#b2182b")) + scale_shape_manual(values=c(17, 16, 15)) +
  labs(x=NULL, y="Fall turnover date", tag="A") + th +
  theme(panel.grid.major.x=element_blank())

dens <- rbind(bd[doy >= 182 & doy <= 335, .(year, doy, d, grp=as.character(year), src="buoy")],
              jd[doy >= 182, .(year=1906, doy, d, grp="1906", src="1906 (Birge & Juday)")])
pb <- ggplot(dens, aes(doy, d, group=grp)) +
  geom_line(data=dens[src != "1906 (Birge & Juday)"], color="grey60", linewidth=0.25, alpha=0.7) +
  geom_line(data=dens[src == "1906 (Birge & Juday)"], color="black", linewidth=0.6) +
  geom_point(data=dens[src == "1906 (Birge & Juday)"], color="black", size=0.8, shape=17) +
  geom_hline(yintercept=THR, linetype="dashed", linewidth=0.35) +
  annotate("text", x=330, y=THR, label="0.3 kg m⁻³", vjust=-0.5, hjust=1, size=2.3) +
  annotate("text", x=186, y=0.75, label="1906", hjust=0, size=2.5, fontface="bold") +
  annotate("text", x=186, y=0.5, label=sprintf("Buoy, %d–%d", min(bd$year), max(bd$year)), hjust=0, size=2.3, color="grey45") +
  scale_x_continuous(breaks=c(182, 213, 244, 274, 305, 335), labels=c("Jul","Aug","Sep","Oct","Nov","Dec"),
                     limits=c(182, 335)) +
  labs(x=NULL, y="Bottom − surface density (kg m⁻³)", tag="B") + th
g <- pa + pb + plot_layout(widths=c(1.15, 1))
ggsave("figures/figSX_Mendota_turnover.png", g, width=6.5, height=3.2, dpi=500, bg="white")
cat("wrote figures/figSX_Mendota_turnover.png\n")

write_captions <- function(new_caps){
  path <- "figures/captions.csv"
  old <- if(file.exists(path)) fread(path) else data.table(file=character(), title=character(), caption=character())
  fwrite(rbind(old[!file %in% new_caps$file], new_caps), path); cat("wrote", path, "\n")
}
write_captions(data.table(
  file="figures/figSX_Mendota_turnover.png",
  title="Timing of fall turnover in Lake Mendota",
  caption=sprintf(paste0(
    "(A) Date of fall turnover, defined as the day the density difference between bottom (18-20 m) and ",
    "surface (0-2 m) water falls below 0.3 kg m-3 after mid-July and remains below it for 30 days. Red squares: ",
    "daily buoy data (", BUOY_YRS, "); blue circles: biweekly profiles (NTL-LTER, with additional 2006-2019 profiles), ",
    "with bars spanning the last stratified and first mixed profile; black triangle: 1906 profiles of Birge and ",
    "Juday. Line: Theil-Sen trend 1995-", LAST, " in the best estimate per year (buoy where available, otherwise profiles) (+%.1f days per decade, Mann-Kendall p = %.2f). (B) Bottom-minus-surface ",
    "density difference (18-20 m minus 0-2 m) from July to November for each buoy year (grey) and for 1906 (black); dashed line, ",
    "turnover threshold."), 10*S_best$slope, S_best$p)))
