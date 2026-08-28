## August DO profile shift, all NTL lakes, 1995-2005 vs 2016-2025.
## Quantifies how the bottom-anoxic boundary (depth where DO permanently drops below 2 mg/L) moved and
## whether the anoxic layer expanded or contracted between eras.
suppressMessages({library(data.table); library(ggplot2); library(lubridate)})

meta <- data.table(
  lakeid=c("TR","BM","CR","SP","AL","TB","CB","ME","MO","FI","WI"),
  name  =c("Trout","Big Muskellunge","Crystal","Sparkling","Allequash","Trout Bog","Crystal Bog",
           "Mendota","Monona","Fish","Wingra"),
  region=c(rep("Northern",7), rep("Southern",4)),
  trophic=c("oligotrophic","oligotrophic","oligotrophic","oligotrophic","mesotrophic","dystrophic","dystrophic",
            "eutrophic","eutrophic","mesotrophic","eutrophic"),
  zmax  =c(35.7,21.3,20.4,20.0,8.0,7.9,2.5,25.3,22.5,18.9,4.0))
tro <- c(oligotrophic="oligotrophic",mesotrophic="mesotrophic",dystrophic="dystrophic (bog)",eutrophic="eutrophic")
early<-"#2c7fb8"; recent<-"#d7301f"

## Crystal 2012-13 dropped everywhere: whole-lake mixing experiment, not natural conditions
p <- fread("data/profiles_clean.csv")[lakeid %in% meta$lakeid & !is.na(o2) &
                                !(lakeid=="CR" & year4 %in% c(2012,2013))]
p[, `:=`(month=month(as.Date(sampledate)), year=year4)]
p <- p[month==8]                                            # August only
p[, era := fifelse(year<=2005,"1995–2005", fifelse(year>=2016,"2016–2025", NA_character_))]
do <- p[!is.na(era), .(o2=mean(o2)), by=.(lakeid, era, depth=round(depth))]
do <- merge(do, meta, by="lakeid")
do <- do[depth<=zmax]

## --- depth of the persistent bottom HYPOXIC boundary (DO permanently <3 mg/L to the bottom) ---
## 3 mg/L is the conventional hypoxia / fish-stress threshold; 2 mg/L (anoxia) was used previously.
anoxic_boundary <- function(depth, o2, thresh=3){
  o <- order(depth); depth <- depth[o]; o2 <- o2[o]
  below <- o2 < thresh; n <- length(below)
  if(!any(below) || !below[n]) return(NA_real_)   # never hypoxic, or bottom itself is oxygenated
  k <- n; while(k>1 && below[k-1]) k <- k-1        # walk up from the bottom while still hypoxic
  if(k==1) return(depth[1])
  approx(o2[c(k-1,k)], depth[c(k-1,k)], xout=thresh)$y
}
## --- metalimnetic O2 maximum: the most prominent subsurface local peak in the profile (if any) ---
meta_peak <- function(depth, o2, min_depth=1, min_prom=0.5){
  o <- order(depth); depth <- depth[o]; o2 <- o2[o]
  keep <- depth>=min_depth; d <- depth[keep]; y <- o2[keep]; n <- length(y)
  if(n<3) return(c(pkz=NA_real_, pko2=NA_real_))
  ismax <- c(FALSE, y[2:(n-1)]>y[1:(n-2)] & y[2:(n-1)]>y[3:n], FALSE)
  if(!any(ismax)) return(c(pkz=NA_real_, pko2=NA_real_))
  cand <- which(ismax)
  prom <- sapply(cand, function(i) y[i] - max(min(y[1:i]), min(y[i:n])))
  best <- cand[which.max(prom)]
  if(max(prom) < min_prom) return(c(pkz=NA_real_, pko2=NA_real_))
  c(pkz=d[best], pko2=y[best])
}
pk <- do[, { r<-meta_peak(depth,o2); .(pkz=r["pkz"], pko2=r["pko2"]) }, by=.(lakeid,era)]
pkw <- merge(pk[era=="1995–2005",.(lakeid,pkz_e=pkz,pko2_e=pko2)],
             pk[era=="2016–2025", .(lakeid,pkz_r=pkz,pko2_r=pko2)], by="lakeid", all=TRUE)
pkw[, peaklab := mapply(function(ze,zr,oe,or){
  if(is.na(ze) && is.na(zr)) return("")
  if(is.na(ze)) return(sprintf("\nO2max new: %.1f at %.0fm", or, zr))
  if(is.na(zr)) return("\nO2max lost")
  dz <- zr-ze; do2 <- or-oe
  mv  <- if(dz < -0.5) sprintf("up %.1fm", -dz) else if(dz>0.5) sprintf("down %.1fm", dz) else "depth held"
  mag <- if(do2 > 0.5) sprintf("+%.1f", do2) else if(do2 < -0.5) sprintf("%.1f", do2) else "~same"
  sprintf("\nO2max %s, %s", mv, mag)
}, pkz_e, pkz_r, pko2_e, pko2_r)]

zb <- do[, .(zb=anoxic_boundary(depth,o2)), by=.(lakeid,era)]
wide <- dcast(zb, lakeid~era, value.var="zb")
setnames(wide, c("1995–2005","2016–2025"), c("zb_early","zb_recent"))
wide <- merge(wide, meta[,.(lakeid,zmax)], by="lakeid")
wide[, `:=`(th_early = zmax-fifelse(is.na(zb_early),zmax,zb_early), th_recent = zmax-fifelse(is.na(zb_recent),zmax,zb_recent))]
wide[, shiftlab := mapply(function(ze, zr, the, thr){
  if(is.na(ze) && is.na(zr)) return("no persistent hypoxia")
  if(is.na(ze)) return(sprintf("hypoxia emerged <%.0fm", zr))
  if(is.na(zr)) return("hypoxia gone")
  dz <- zr-ze; dth <- thr-the
  moved <- if(dz < -0.5) sprintf("bdy up %.1fm", -dz) else
           if(dz >  0.5) sprintf("bdy down %.1fm", dz) else "bdy held"
  sizechg <- if(dth >  0.5) sprintf(", +%.1fm", dth) else
             if(dth < -0.5) sprintf(", %.1fm", dth) else ", ~same"
  paste0(moved, sizechg)
}, zb_early, zb_recent, th_early, th_recent)]

meta2 <- merge(meta, wide[,.(lakeid,zb_early,zb_recent,shiftlab)], by="lakeid")
meta2 <- merge(meta2, pkw[,.(lakeid,peaklab)], by="lakeid")
meta2[, lab := sprintf("%s · %.0f m\n%s%s", name, zmax, shiftlab, peaklab)]
ordtab <- meta2[order(region,-zmax)]
do   <- merge(do,   meta2[,.(lakeid,lab)], by="lakeid")[, lab := factor(lab, levels=ordtab$lab)]
segs <- merge(wide,  meta2[,.(lakeid,lab)], by="lakeid")[, lab := factor(lab, levels=ordtab$lab)]
peaksegs <- merge(pkw, meta2[,.(lakeid,lab)], by="lakeid")[, lab := factor(lab, levels=ordtab$lab)]
## geom_path joins points in ROW order; guard against a merge re-sorting them out of depth order
setorder(do, lakeid, era, depth)

g <- ggplot(do, aes(o2, depth, color=era)) +
  geom_vline(xintercept=3, linetype=3, color="grey65") +
  geom_path(linewidth=0.5) + geom_point(size=0.4) +
  geom_hline(data=segs[!is.na(zb_early)],  aes(yintercept=zb_early),  color=early,  linetype=2, linewidth=0.3) +
  geom_hline(data=segs[!is.na(zb_recent)], aes(yintercept=zb_recent), color=recent, linetype=2, linewidth=0.3) +
  geom_segment(data=segs[!is.na(zb_early) & !is.na(zb_recent)],
    aes(x=3, xend=3, y=zb_early, yend=zb_recent), inherit.aes=FALSE,
    arrow=arrow(length=unit(0.05,"in")), color="grey30", linewidth=0.35) +
  geom_point(data=peaksegs[!is.na(pkz_e)], aes(pko2_e, pkz_e), inherit.aes=FALSE,
    shape=23, size=1.3, fill=early, color="black", stroke=0.25) +
  geom_point(data=peaksegs[!is.na(pkz_r)], aes(pko2_r, pkz_r), inherit.aes=FALSE,
    shape=23, size=1.3, fill=recent, color="black", stroke=0.25) +
  geom_segment(data=peaksegs[!is.na(pkz_e) & !is.na(pkz_r)],
    aes(x=pko2_e, xend=pko2_r, y=pkz_e, yend=pkz_r), inherit.aes=FALSE,
    arrow=arrow(length=unit(0.05,"in")), color="grey30", linewidth=0.35, linetype=1) +
  facet_wrap(~lab, scales="free_y", ncol=3) +
  scale_y_reverse() +
  scale_color_manual(values=c("1995–2005"=early,"2016–2025"=recent), name=NULL) +
  labs(title="August DO profile shift, NTL lakes: 1995–2005 vs 2016–2025",
       subtitle="Dashed lines/arrow: bottom hypoxic boundary, where DO stays below 3 mg/L to the bottom (blue = 1995–2005, red = 2016–2025).\nDiamonds/arrow: metalimnetic O2 maximum. Northern lakes first, then Southern, by max depth. Crystal 2012–13 excluded.",
       x="Dissolved oxygen (mg/L)", y="Depth (m)") +
  theme_minimal(base_size=7) + theme(panel.grid.minor=element_blank(), legend.position="top",
    plot.title=element_text(face="bold",size=8.5), plot.subtitle=element_text(color="grey40",size=5.4, lineheight=1.15),
    legend.text=element_text(size=6), legend.key.size=unit(0.3,"cm"), legend.margin=margin(0,0,0,0),
    axis.text=element_text(size=5.4), axis.title=element_text(size=6.4),
    strip.text=element_text(size=5.4, lineheight=1.05), panel.spacing=unit(0.28,"lines"))
ggsave("figures/fig07_alllakes_metO2_aug.png", g, width=6, height=7.6, dpi=600, bg="white")
cat("wrote figures/fig07_alllakes_metO2_aug.png\n")
