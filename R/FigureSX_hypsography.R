## FigureSX_hypsography.R: hypsographic curves for the 9 NTL lakes, from data/Bathymetry/lakehypsometry.csv
## (the same file that weights depth bins by volume in Figure6_volatility.R panel b).
## lakehypsometry.csv gives hp_factor = fraction of lake volume per layer (each lake sums to 1); see the
## layer-convention note below. Area at a layer's mid-depth is proportional to hp_factor / layer thickness;
## it is scaled so the shallowest layer equals the lake's surface area (lake_characteristics.csv).
##   (a) one panel per lake, absolute area (ha) vs depth
##   (b) all lakes on relative axes (area / surface area vs depth / max depth), with a cone for reference
## (data/Bathymetry/lakemorphometry.csv is not used: it covers only the sampled basin of Allequash and Trout
##  and has no areas for Mendota or Monona.)
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(ggrepel)})

meta <- data.table(
  lakeid=c("TR","BM","CR","SP","AL","TB","CB","ME","MO"),
  name  =c("Trout","Big Muskellunge","Crystal","Sparkling","Allequash","Trout Bog","Crystal Bog","Mendota","Monona"),
  lake  =c("Trout Lake","Big Muskellunge Lake","Crystal Lake","Sparkling Lake","Allequash Lake","Trout Bog",
           "Crystal Bog","Lake Mendota","Lake Monona"),
  region=c(rep("Northern",7), rep("Southern",2)))
lc <- fread("data/lake_characteristics.csv")[, .(lake=waterbody_name, A0_ha=waterbody_area_ha, zmean=avg_depth_m, zmax=max_depth_m)]
meta <- merge(meta, lc, by="lake")

h <- fread("data/Bathymetry/lakehypsometry.csv")[lakeid %in% meta$lakeid]
setorder(h, lakeid, depth)
## layer bounds. Northern lakes: `depth` = layer BOTTOM. Mendota, Monona: `depth` = layer CENTRE
## (0.25 = 0-0.5 m, then d-0.5 to d+0.5).
h[, bot := if(.BY$lakeid %in% c("ME","MO")) fifelse(depth == min(depth), 2*depth, depth + 0.5) else depth, by=lakeid]
h[, top := shift(bot, fill=0), by=lakeid]
h <- h[bot > top]                                          # drop the 0-thickness row at the surface
## the 0.5-m layers of the deep northern lakes are 1-m values split in half (pairs repeat), which draws as
## stair-steps; merge them back into 1-m layers for plotting
h[, grp := if(median(bot - top) == 0.5 && !.BY$lakeid %in% c("ME","MO")) ceiling(bot) else bot, by=lakeid]
h <- h[, .(top=min(top), bot=max(bot), hp_factor=sum(hp_factor)), by=.(lakeid, grp)]
h[, depth := bot]
h[, `:=`(zmid=(top+bot)/2, dens=hp_factor/(bot-top))]     # proportional to area
h[, rel_area := dens/max(dens), by=lakeid]
## curve: surface (relative area 1) -> layer mid-depths -> 0 at the deepest layer's bottom
curve <- rbind(h[, .(depth_m=0, rel_area=1), by=lakeid], h[, .(lakeid, depth_m=zmid, rel_area)],
               h[, .(depth_m=max(depth), rel_area=0), by=lakeid])
setorder(curve, lakeid, depth_m)
curve <- merge(curve, meta, by="lakeid")
curve[, `:=`(area_ha=rel_area*A0_ha, rel_depth=depth_m/max(depth_m)), by=lakeid]

ord <- meta[order(region, -zmax), name]                    # north first, then by zmax (repo convention)
meta[, lab := sprintf("%s\n%s ha · mean %.1f m\nmax %.1f m", name, formatC(A0_ha, format="fg", digits=3, big.mark=","), zmean, zmax)]
curve <- merge(curve, meta[, .(lakeid, lab)], by="lakeid")
curve[, lab := factor(lab, levels=meta[match(ord, name), lab])]
curve[, name := factor(name, levels=ord)]

regcol <- c(Northern="#1b7837", Southern="#b2182b")
th <- theme_minimal(base_size=7) + theme(panel.grid.minor=element_blank(),
  strip.text=element_text(size=5.5, hjust=0, lineheight=1.1), plot.tag=element_text(face="bold", size=9),
  legend.position="bottom", legend.key.width=unit(0.6,"cm"))

pa <- ggplot(curve, aes(area_ha, depth_m, color=region)) +
  geom_path(linewidth=0.5) +
  facet_wrap(~lab, ncol=3, scales="free") +
  scale_y_reverse(expand=expansion(mult=c(0.02, 0.04))) +
  scale_x_continuous(expand=expansion(mult=c(0.02, 0.06)), limits=c(0, NA)) +
  scale_color_manual(values=regcol, name=NULL) +
  labs(x="Area (ha)", y="Depth (m)", tag="a") + th

pb <- ggplot(curve, aes(rel_area, rel_depth, group=name, color=region)) +
  geom_line(data=data.table(rel_depth=seq(0,1,0.01))[, rel_area := (1-rel_depth)^2], aes(rel_area, rel_depth),
            inherit.aes=FALSE, color="grey60", linewidth=0.35, linetype="22") +   # cone: A/A0 = (1 - z/zmax)^2
  annotate("text", x=0.12, y=0.55, label="cone", size=1.9, color="grey50", angle=-62) +
  geom_path(linewidth=0.5) +
  geom_text_repel(data=curve[, .SD[which.min(abs(rel_area - 0.5))], by=name],
            aes(rel_area, rel_depth, label=name), size=1.9, show.legend=FALSE, seed=2,
            nudge_x=0.25, direction="y", segment.size=0.2, segment.color="grey60", box.padding=0.12,
            xlim=c(0.75, 1.3), min.segment.length=0) +
  scale_y_reverse(limits=c(1, 0), breaks=seq(0,1,0.25)) +
  scale_x_continuous(limits=c(0, 1.3), breaks=seq(0,1,0.25)) +
  scale_color_manual(values=regcol, guide="none") +
  coord_cartesian(clip="off") +
  labs(x="Area / surface area", y="Depth / maximum depth", tag="b") + th

p <- pa + pb + plot_layout(widths=c(2.2, 1), guides="collect") & theme(legend.position="bottom")
ggsave("figures/figSX_hypsography.png", p, width=6.5, height=4.8, dpi=500, bg="white")
cat("wrote figures/figSX_hypsography.png\n")

write_captions <- function(new_caps){
  path <- "figures/captions.csv"
  old <- if(file.exists(path)) fread(path) else data.table(file=character(),title=character(),caption=character())
  fwrite(rbind(old[!file %in% new_caps$file], new_caps), path); cat("wrote", path, "\n")
}
write_captions(data.table(
  file="figures/figSX_hypsography.png",
  title="Hypsographic curves of the study lakes",
  caption=paste0(
    "Lake area at depth, from layer volume fractions (0.25-1 m layers; data/Bathymetry/lakehypsometry.csv), with area ",
    "proportional to each layer's volume fraction divided by its thickness and scaled to the lake's surface area. ",
    "(a) Absolute area (ha) for each lake, with surface area, mean and maximum depth. Northern lakes in green, ",
    "southern in red. (b) All lakes on relative axes (area / surface area against depth / maximum depth); the dashed ",
    "grey curve is a perfect cone, A/A0 = (1 - z/zmax)^2. Curves further right have relatively more deep area ",
    "(steep-sided basins); curves nearer the cone have most area in shallow water. The same volume fractions weight ",
    "depth bins in Figure 6b.")))
