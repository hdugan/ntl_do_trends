## Fig 4: August chlorophyll profile shift, northern NTL lakes, 1995-2005 vs 2016-2025. Data: EDI 35.
## Chlorophyll counterpart to Figure3_o2_profiles.R. Where that figure tracks the DO profile shape,
## this tracks the DEEP CHLOROPHYLL MAXIMUM: its depth and its magnitude.
##
## DEPTH SELECTION: split into the two eras first, then within each era keep only depths sampled
## more than MIN_N times -- the routinely-occupied levels, not opportunistic one-offs. Each era is
## then drawn on its own surviving depth set; the two need not match.
suppressMessages({library(data.table); library(ggplot2); library(ggtext)})

## Four deep oligotrophic lakes only. Allequash and the two bogs are excluded: Allequash's deepest
## level carries benthic material rather than plankton, and neither bog retains a usable recent-era
## profile (Trout Bog's recent August casts are surface-only), so their era comparison would report
## a sampling change rather than a chlorophyll change.
meta <- data.table(
  lakeid=c("TR","BM","CR","SP"),
  name  =c("Trout","Big Muskellunge","Crystal","Sparkling"),
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
ch[, era := fifelse(year>=1995 & year<=2005, "1995-2005",
             fifelse(year>=2016, "2016-2025", NA_character_))]
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
w <- merge(dcm[era=="1995-2005", .(lakeid, pz_e=pz, pc_e=pchl)],
           dcm[era=="2016-2025", .(lakeid, pz_r=pz, pc_r=pchl)], by="lakeid", all=TRUE)
w[, dcmlab := mapply(function(ze,zr,ce,cr){
  if(is.na(ze)||is.na(zr)) return("")
  dz <- zr-ze; dc <- cr-ce
  mv  <- if(dz < -0.5) sprintf("up %.0fm", -dz) else if(dz > 0.5) sprintf("down %.0fm", dz) else "depth held"
  mag <- if(dc >  0.3) sprintf("+%.1f", dc) else
         if(dc < -0.3) sprintf("%.1f", dc) else "~same"
  sprintf("DCM %s, %sµg/L", mv, mag)
}, pz_e, pz_r, pc_e, pc_r)]

## panel titles: lake name (bold) + DCM shift, same convention as Figure2/Figure3
m2 <- merge(meta, w[, .(lakeid, dcmlab)], by="lakeid")
m2[, striplab := sprintf("<span style='font-size:7pt;font-weight:bold'>%s</span><br>%s", name, dcmlab)]
ord <- m2[order(-zmax), striplab]
prof <- merge(prof, m2[, .(lakeid,striplab)], by="lakeid")[, strip := factor(striplab, levels=ord)]
pk   <- merge(w,    m2[, .(lakeid,striplab)], by="lakeid")[, strip := factor(striplab, levels=ord)]

## geom_path joins points in ROW order and the merges above re-sort by the join key: sort by depth
setorder(prof, lakeid, era, depth)

th <- theme_minimal(base_size=6.5) + theme(
  panel.grid.minor=element_blank(), panel.spacing=unit(0.4,"lines"),
  legend.position="bottom", legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.22,"cm"),
  legend.title=element_text(size=6.5),
  axis.text=element_text(size=5), axis.text.x=element_text(size=6),
  axis.title=element_text(size=7.5),
  strip.text=element_markdown(size=5, lineheight=1.2))

g <- ggplot(prof, aes(chl, depth, color=era)) +
  geom_path(linewidth=0.6) + geom_point(size=0.5) +
  facet_wrap(~strip, scales="free", nrow=1) +
  scale_y_reverse() +
  scale_color_manual(values=c("1995-2005"=early, "2016-2025"=recent), name=NULL) +
  labs(x="Chlorophyll (µg/L)", y="Depth (m)") +
  th
ggsave("figures/fig04_chl_profiles.png", g, width=6.5, height=2.7, dpi=500, bg="white")
cat("wrote figures/fig04_chl_profiles.png\n")

## Title + long explanatory paragraph are NOT drawn on the PNG -- written to the shared
## figures/captions.csv instead (see Figure1_rate_of_change.R for the same convention).
write_captions <- function(new_caps){
  path <- "figures/captions.csv"
  old <- if(file.exists(path)) fread(path) else data.table(file=character(),title=character(),caption=character())
  fwrite(rbind(old[!file %in% new_caps$file], new_caps), path)
  cat("wrote", path, "\n")
}
caption <- paste0(
  "August chlorophyll profile shift, 4 northern deep oligotrophic NTL lakes: 1995-2005 vs 2016-2025. ",
  "Profiles are era MEDIANS; the DCM (deep chlorophyll maximum) shift quoted in each panel is the change ",
  "in the depth and magnitude of the profile peak. Panels ordered by maximum depth.\n",
  "Depths are selected within each era (>", MIN_N, " August measurements at that level); the two eras ",
  "are drawn on their own depth sets and need not match. Allequash and the two bog lakes excluded: ",
  "Allequash's deepest level carries benthic material rather than plankton, and neither bog retains a ",
  "usable recent-era profile. Crystal 2012-2013 excluded (whole-lake mixing experiment).")
write_captions(data.table(file="figures/fig04_chl_profiles.png",
  title="August chlorophyll profile shift", caption=caption))

print(merge(w, meta[, .(lakeid,name)], by="lakeid")[order(-pz_e), .(
  name, dcm_z_early=pz_e, dcm_z_recent=pz_r, shift_m=pz_r-pz_e,
  dcm_chl_early=round(pc_e,2), dcm_chl_recent=round(pc_r,2), change=round(pc_r-pc_e,2))])
