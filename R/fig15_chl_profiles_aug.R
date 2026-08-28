## Fig 15: August chlorophyll profile shift, northern NTL lakes, 1995-2005 vs 2016-2025.
## Chlorophyll counterpart to fig07. Where fig07 tracks the hypoxic boundary and the metalimnetic
## O2 maximum, this tracks the DEEP CHLOROPHYLL MAXIMUM: its depth and its magnitude.
##
## DEPTH SELECTION: split into the two eras first, then within each era keep only depths sampled
## more than MIN_N times -- the routinely-occupied levels, not opportunistic one-offs. Each era is
## then drawn on its own surviving depth set; the two need not match.
suppressMessages({library(data.table); library(ggplot2)})

## Four deep oligotrophic lakes only. Allequash and the two bogs are excluded: Allequash's deepest
## level carries benthic material rather than plankton, and neither bog retains a usable recent-era
## profile (Trout Bog's recent August casts are surface-only), so their era comparison would report
## a sampling change rather than a chlorophyll change.
meta <- data.table(
  lakeid=c("TR","BM","CR","SP"),
  name  =c("Trout","Big Muskellunge","Crystal","Sparkling"),
  trophic=c("oligotrophic","oligotrophic","oligotrophic","oligotrophic"),
  zmax  =c(35.7,21.3,20.4,20.0))
early <- "#2c7fb8"; recent <- "#d7301f"
MIN_N <- 7         # keep a depth only if sampled MORE than this many times within the era

## QC FLAGS ARE NOT USED TO EXCLUDE DATA. An earlier version dropped every flagged value, which
## turned out to bias the result badly: the "V" flag occurs at 0% of samples above 12 m but 18-23%
## at 17-19 m, and only from the 2010s onward. Excluding it therefore strips recent DEEP samples
## specifically -- exactly the cells this figure is about. At Sparkling 18 m in August 2016-2024
## that left 2 of 17 measurements, and those 2 were the lowest of the set (median 1.2 vs 5.6),
## manufacturing an apparent collapse of deep chlorophyll that is not in the data.
## All finite values are kept; the negative-value rule below plus median aggregation do the work.
ch <- fread("data/chl.csv")[lakeid %in% meta$lakeid & is.finite(chlor)]
## Negative chlorophyll: small negatives are blank/baseline noise around a true value of ~0, so
## floor them at 0 rather than discard (discarding would bias the low end upward). Below -5 ug/L is
## not recoverable noise and is deleted.
cat("negative chl: ", sum(ch$chlor < 0), " total; ", sum(ch$chlor >= -5 & ch$chlor < 0),
    " floored to 0, ", sum(ch$chlor < -5), " deleted\n", sep="")
ch <- ch[chlor >= -5]
ch[chlor < 0, chlor := 0]
ch[, `:=`(mon=month(as.Date(sampledate)), year=year4, date=as.Date(sampledate), z=round(depth))]
ch <- ch[mon==8 & !(lakeid=="CR" & year %in% c(2012,2013))]     # August; Crystal mixing yrs out
ch[, era := fifelse(year>=1995 & year<=2005, "1995–2005",
             fifelse(year>=2016, "2016–2025", NA_character_))]
ch <- ch[!is.na(era)]

## per-era depth selection
keep <- ch[, .(n=.N), by=.(lakeid, era, z)][n > MIN_N]
ch   <- merge(ch, keep[, .(lakeid, era, z)], by=c("lakeid","era","z"))
cat("=== depths retained per era (>", MIN_N, "August measurements) ===\n")
print(merge(keep[, .(levels=.N, depths=paste(sort(z), collapse=",")), by=.(lakeid,era)],
            meta[, .(lakeid,name)], by="lakeid")[order(name,era), .(name, era, levels, depths)], nrows=20)

## era profile as the MEDIAN over casts, not the mean: one extreme reading otherwise dominates
## (Trout has a single 137.6 ug/L value at 15 m in 2019 against a median of 4.3 there).
prof <- ch[, .(chl=median(chlor), n=.N), by=.(lakeid, era, depth=z)]
prof <- merge(prof, meta, by="lakeid")

## ---- deep chlorophyll maximum: depth and magnitude of the profile peak ----------------------
dcm <- prof[, .SD[which.max(chl)], by=.(lakeid,era)][, .(lakeid, era, pz=depth, pchl=chl)]
w <- merge(dcm[era=="1995–2005", .(lakeid, pz_e=pz, pc_e=pchl)],
           dcm[era=="2016–2025", .(lakeid, pz_r=pz, pc_r=pchl)], by="lakeid", all=TRUE)
w[, dcmlab := mapply(function(ze,zr,ce,cr){
  if(is.na(ze)||is.na(zr)) return("")
  dz <- zr-ze; dc <- cr-ce
  mv  <- if(dz < -0.5) sprintf("up %.0f m", -dz) else if(dz > 0.5) sprintf("down %.0f m", dz) else "depth held"
  mag <- if(dc >  0.3) sprintf("+%.1f µg/L", dc) else
         if(dc < -0.3) sprintf("%.1f µg/L", dc) else "~same"
  sprintf("\nDCM %s, %s", mv, mag)
}, pz_e, pz_r, pc_e, pc_r)]

m2 <- merge(meta, w[, .(lakeid, dcmlab)], by="lakeid")
m2[, lab := sprintf("%s (%s) · %.0f m%s", name, trophic, zmax, dcmlab)]
ord <- m2[order(-zmax)]$lab
prof <- merge(prof, m2[, .(lakeid,lab)], by="lakeid")[, lab := factor(lab, levels=ord)]
pk   <- merge(w,    m2[, .(lakeid,lab)], by="lakeid")[, lab := factor(lab, levels=ord)]

## geom_path joins points in ROW order and the merges above re-sort by the join key: sort by depth
setorder(prof, lakeid, era, depth)

g <- ggplot(prof, aes(chl, depth, color=era)) +
  geom_path(linewidth=0.8) + geom_point(size=1.4) +
  facet_wrap(~lab, scales="free", nrow=1) +
  scale_y_reverse() +
  scale_color_manual(values=c("1995–2005"=early, "2016–2025"=recent), name=NULL) +
  labs(title="August chlorophyll profile shift, northern NTL lakes: 1995–2005 vs 2016–2025",
       subtitle=paste0(
"Profiles are era MEDIANS; the DCM shift quoted in each panel header is the change in the depth and magnitude of the profile peak. Panels ordered by maximum depth.\n",
"Depths are selected within each era (>", MIN_N, " August measurements at that level); the two eras are drawn on their own depth sets and need not match. Crystal 2012–13 excluded."),
       x="Chlorophyll (µg/L)", y="Depth (m)") +
  theme_minimal(base_size=10) +
  theme(panel.grid.minor=element_blank(), legend.position="top",
        plot.title=element_text(face="bold", size=11.5),
        plot.subtitle=element_text(color="grey40", size=7.6, lineheight=1.2),
        strip.text=element_text(size=7.6, lineheight=1.05))
ggsave("figures/fig15_chl_profiles_aug.png", g, width=11, height=3.9, dpi=500, bg="white")
cat("wrote figures/fig15_chl_profiles_aug.png\n")

print(merge(w, meta[, .(lakeid,name)], by="lakeid")[order(-pz_e), .(
  name, dcm_z_early=pz_e, dcm_z_recent=pz_r, shift_m=pz_r-pz_e,
  dcm_chl_early=round(pc_e,2), dcm_chl_recent=round(pc_r,2), change=round(pc_r-pc_e,2))])
