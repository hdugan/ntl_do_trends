## Fig 3: August DO profile by decade, all 11 NTL lakes. Data: EDI 29.
## Four August MEDIAN profiles per lake -- 1985-1995, 1995-2005, 2005-2015, 2015-2025 -- overlaid
## on one panel per lake, coloured blue (earliest) -> green -> orange -> red (most recent) so a
## profile sliding right/left/up/down across the four lines reads as a multi-decade progression,
## not just a two-era before/after snapshot (the earlier version of this figure compared only
## 1995-2005 vs 2016-2025). No anoxic-boundary or metalimnetic-O2-max annotations here -- just the
## raw shape of the four profiles; that annotation layer can come back once the 4-era shape is
## worth quantifying that way.
suppressMessages({library(data.table); library(ggplot2); library(ggtext)})

meta <- data.table(
  lakeid=c("TR","BM","CR","SP","AL","TB","CB","ME","MO","FI","WI"),
  name  =c("Trout","Big Muskellunge","Crystal","Sparkling","Allequash","Trout Bog","Crystal Bog",
           "Mendota","Monona","Fish","Wingra"),
  region=c(rep("Northern",7), rep("Southern",4)),
  zmax  =c(35.7,21.3,20.4,20.0,8.0,7.9,2.5,25.3,22.5,18.9,4.0))

## Crystal 2012-13 dropped everywhere: whole-lake mixing experiment, not natural conditions
p <- fread("data/profiles_clean.csv")[lakeid %in% meta$lakeid & !is.na(o2) &
                                !(lakeid=="CR" & year4 %in% c(2012,2013))]
p[, `:=`(month=month(as.Date(sampledate)), year=year4)]
p <- p[month==8]                                            # August only

## four non-overlapping decade bins; last one includes 2025. findInterval gives 0 for years
## before 1985 and 5 for years >=2026 -- clip to 1..4 first so lv[idx] doesn't silently drop rows
## (lv[0] is character(0), which would shrink the column instead of producing NA)
breaks <- c(1985,1995,2005,2015,2026)
lv <- c("1985-1995","1995-2005","2005-2015","2015-2025")
p[, era := { idx <- findInterval(year, breaks)
             fifelse(idx>=1 & idx<=4, lv[pmin(pmax(idx,1L),4L)], NA_character_) }]
p <- p[!is.na(era)]
p[, era := factor(era, levels=lv)]

## MEDIAN across casts (not mean): robust to the odd bloom-day or post-storm cast, same treatment
## used throughout this project (Figure1/Figure2). Depth rounded to the nearest metre first so
## casts at slightly different depths pool into one profile.
do <- p[, .(o2=median(o2)), by=.(lakeid, era, depth=round(depth))]
do <- merge(do, meta, by="lakeid")
do <- do[depth<=zmax]
setorder(do, lakeid, era, depth)   # geom_path joins in row order; guard against merge re-sorting

## strip labels: lake name (bold) + region tag, same convention as Figure2_clarity_trends.R
meta[, striplab := sprintf("<span style='font-size:7pt;font-weight:bold'>%s</span> (%s)",
                            name, tolower(region))]
## blank placeholder panel right after the 7 northern lakes, so with ncol=4 the northern lakes
## fill their own two rows (4 + 3-plus-a-gap) and the 4 southern lakes start fresh on row 3,
## rather than Mendota filling that trailing gap
ord <- meta[order(region,-zmax), striplab]
ord <- append(ord, "", after=uniqueN(meta[region=="Northern", lakeid]))
do <- merge(do, meta[,.(lakeid,striplab)], by="lakeid")
do[, strip := factor(striplab, levels=ord)]

pal <- setNames(c("#2166ac","#1a9850","#e08214","#b2182b"), lv)

th <- theme_minimal(base_size=6.5) + theme(
  panel.grid.minor=element_blank(), panel.spacing=unit(0.4,"lines"),
  legend.position="bottom", legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.22,"cm"),
  legend.title=element_text(size=6.5), plot.title=element_text(face="bold",size=9.5),
  plot.subtitle=element_text(face="bold", size=6.5, color="grey30"),
  axis.text=element_text(size=5), axis.text.x=element_text(size=6),
  axis.title.y=element_text(size=7.5),
  strip.text=element_markdown(size=5, lineheight=1.2))

g <- ggplot(do, aes(o2, depth, color=era)) +
  geom_path(linewidth=0.6) + geom_point(size=0.5) +
  facet_wrap(~strip, scales="free_y", ncol=4, drop=FALSE) +
  scale_y_reverse() +
  scale_color_manual(values=pal, name=NULL) +
  labs(x="Dissolved oxygen (mg/L)", y="Depth (m)") +
  th
ggsave("figures/fig03_o2_profiles.png", g, width=6.5, height=6, dpi=500, bg="white")
cat("wrote figures/fig03_o2_profiles.png\n")

## Title + long explanatory paragraph are NOT drawn on the PNG -- written to the shared
## figures/captions.csv instead (see Figure1_rate_of_change.R for the same convention).
write_captions <- function(new_caps){
  path <- "figures/captions.csv"
  old <- if(file.exists(path)) fread(path) else data.table(file=character(),title=character(),caption=character())
  fwrite(rbind(old[!file %in% new_caps$file], new_caps), path)
  cat("wrote", path, "\n")
}
caption <- paste0(
  "August dissolved-oxygen depth profile, all 11 lakes, four overlaid decade MEDIANS: 1985-1995, ",
  "1995-2005, 2005-2015, 2015-2025 (blue -> green -> orange -> red, earliest to most recent). Each ",
  "profile is the per-depth median O2 (mg/L) across every August cast at that lake within that decade ",
  "(depth rounded to the nearest metre before pooling). Panels ordered north then south, by max depth ",
  "(zmax); region shown after each lake name. Crystal 2012-13 excluded (whole-lake mixing experiment).\n",
  "No anoxic-boundary or metalimnetic-O2-maximum annotations on this version -- just the raw profile ",
  "shape across the four decades, so a lake's trajectory (deepening/shoaling hypoxia, a strengthening ",
  "or weakening subsurface O2 peak) can be read directly from how the four lines separate.")
write_captions(data.table(file="figures/fig03_o2_profiles.png",
  title="August dissolved-oxygen profiles by decade", caption=caption))
