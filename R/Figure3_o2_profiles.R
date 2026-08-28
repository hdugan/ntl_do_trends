## Fig 3: August DO (% saturation) profile by decade, all 11 NTL lakes. Data: EDI 29.
## Four August MEDIAN profiles per lake -- 1985-1995, 1995-2005, 2005-2015, 2015-2025 -- overlaid
## on one panel per lake, coloured blue (earliest) -> green -> orange -> red (most recent) so a
## profile sliding right/left/up/down across the four lines reads as a multi-decade progression,
## not just a two-era before/after snapshot (the earlier version of this figure compared only
## 1995-2005 vs 2016-2025). No anoxic-boundary or metalimnetic-O2-max annotations here -- just the
## raw shape of the four profiles; that annotation layer can come back once the 4-era shape is
## worth quantifying that way.
suppressMessages({library(data.table); library(ggplot2); library(ggtext); library(patchwork)})

meta <- data.table(
  lakeid=c("TR","BM","CR","SP","AL","TB","CB","ME","MO","FI","WI"),
  name  =c("Trout","Big Muskellunge","Crystal","Sparkling","Allequash","Trout Bog","Crystal Bog",
           "Mendota","Monona","Fish","Wingra"),
  region=c(rep("Northern",7), rep("Southern",4)),
  zmax  =c(35.7,21.3,20.4,20.0,8.0,7.9,2.5,25.3,22.5,18.9,4.0))

## Crystal 2012-13 dropped everywhere: whole-lake mixing experiment, not natural conditions
p <- fread("data/profiles_clean.csv")[lakeid %in% meta$lakeid & !is.na(o2sat) &
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
do <- p[, .(o2sat=median(o2sat)), by=.(lakeid, era, depth=round(depth))]
do <- merge(do, meta, by="lakeid")
do <- do[depth<=zmax]
setorder(do, lakeid, era, depth)   # geom_path joins in row order; guard against merge re-sorting

## panel titles: lake name (bold) + region tag, same convention as Figure2_clarity_trends.R
meta[, striplab := sprintf("<span style='font-size:7pt;font-weight:bold'>%s</span> (%s)",
                            name, tolower(region))]
striplabs <- setNames(meta$striplab, meta$lakeid)
lake_ord <- meta[order(region,-zmax), lakeid]
north <- lake_ord[1:7]; south <- lake_ord[8:11]

pal <- setNames(c("#2166ac","#1a9850","#e08214","#b2182b"), lv)
xr <- range(do$o2sat, na.rm=TRUE)

## Built as individual panels + patchwork (not facet_wrap) so the gap left by the 7th northern
## lake (northern doesn't fill a full row of 4) can be a TRUE blank space via plot_spacer() --
## an empty facet_wrap panel still draws its axes/gridlines, which read as a broken plot rather
## than an intentional gap.
## legend.position="none": patchwork's automatic guides="collect" kept producing two duplicate
## copies of the (structurally identical) legend across these 11 independent ggplots rather than
## merging them into one, for reasons that resisted a few different fixes -- extracted and placed
## manually below instead (see extract_legend()).
th <- theme_minimal(base_size=6.5) + theme(
  panel.grid.minor=element_blank(), panel.spacing=unit(0.4,"lines"),
  legend.position="none",
  plot.title=element_markdown(size=5, lineheight=1.2, hjust=0.5),
  axis.text=element_text(size=5), axis.text.x=element_text(size=6),
  axis.title=element_text(size=7.5))

## x-axis TITLE is a single shared caption below the whole grid (added via plot_annotation), not
## repeated per panel -- but every panel shows its own 0/50/100 x tick numbers and its own depth
## tick numbers, so each profile is readable without hunting for the nearest axis
mkpanel <- function(lk, show_ytitle){
  d <- do[lakeid==lk]
  ggplot(d, aes(o2sat, depth, color=era)) +
    geom_path(linewidth=0.6) + geom_point(size=0.5) +
    expand_limits(x=xr) + scale_x_continuous(breaks=c(0,50,100)) + scale_y_reverse() +
    scale_color_manual(values=pal, name=NULL, drop=FALSE, limits=lv) +
    labs(title=striplabs[[lk]], x=NULL, y=if(show_ytitle) "Depth (m)" else NULL) +
    th
}

## leftmost column (north-index 1 & 5, south-index 1) keeps the "Depth (m)" axis title; every
## panel still gets its own depth tick numbers regardless
panels <- lapply(seq_along(north), function(i) mkpanel(north[i], show_ytitle = i %in% c(1,5)))
## true blank gap panel: no data, no axes, no legend
panels[[8]] <- ggplot() + theme_void()
panels <- c(panels, lapply(seq_along(south), function(i) mkpanel(south[i], show_ytitle = i==1)))

## manually extract the (one, shared) legend from a reference panel and place it as an extra
## row below the 3x4 grid, rather than relying on patchwork's guides="collect"
extract_legend <- function(p){
  pdf(NULL)  # ggplotGrob() opens a graphics device if none is open, leaving an Rplots.pdf behind
  on.exit(dev.off())
  g <- ggplotGrob(p + theme(legend.position="bottom", legend.title=element_text(size=6.5)))
  g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
}
legend_grob <- extract_legend(mkpanel(north[1], show_ytitle=TRUE))

grid <- wrap_plots(panels, ncol=4, nrow=3)
design <- grid / wrap_elements(full=legend_grob) + plot_layout(heights=c(1,0.08))
design <- design + plot_annotation(caption="Dissolved oxygen (% sat)",
  theme=theme(plot.caption=element_text(hjust=0.5, size=7.5, margin=margin(t=2))))
ggsave("figures/fig03_o2_profiles.png", design, width=6.5, height=6, dpi=500, bg="white")
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
  "August dissolved-oxygen (% saturation) depth profile, all 11 lakes, four overlaid decade MEDIANS: ",
  "1985-1995, 1995-2005, 2005-2015, 2015-2025 (blue -> green -> orange -> red, earliest to most recent). ",
  "Each profile is the per-depth median O2 (% sat) across every August cast at that lake within that decade ",
  "(depth rounded to the nearest metre before pooling). Panels ordered north then south, by max depth ",
  "(zmax); region shown after each lake name; northern lakes fill their own two rows, southern lakes start ",
  "fresh on row 3. Crystal 2012-13 excluded (whole-lake mixing experiment).\n",
  "No anoxic-boundary or metalimnetic-O2-maximum annotations on this version -- just the raw profile ",
  "shape across the four decades, so a lake's trajectory (deepening/shoaling hypoxia, a strengthening ",
  "or weakening subsurface O2 peak) can be read directly from how the four lines separate.")
write_captions(data.table(file="figures/fig03_o2_profiles.png",
  title="August dissolved-oxygen profiles by decade", caption=caption))
