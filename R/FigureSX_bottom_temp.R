## FigureSX_bottom_temp.R: spread of bottom-water temperature by month and lake -- the bottom-water
## counterpart of FigureSX_surface_temp_spread.R.
## Bottom temperature = mean wtemp over the deepest 1 m of each cast. A cast only counts if it reaches
## at least 90% of that lake's median cast depth, so a profile that stopped short of the bottom is
## not mistaken for a bottom reading. April-November, Crystal Lake 2012-2013 excluded.
## Output: figSX_bottom_temp_spread.png, the 5th-95th percentile range by lake x month (spread only).
suppressMessages({library(data.table); library(ggplot2)})

meta <- data.table(
  lakeid=c("TR","BM","CR","SP","AL","TB","CB","ME","MO"),
  name  =c("Trout Lake","Big Muskellunge Lake","Crystal Lake","Sparkling Lake","Allequash Lake",
           "Trout Bog","Crystal Bog","Lake Mendota","Lake Monona"),
  region=c(rep("Northern",7), rep("Southern",2)),
  zmax  =c(35.7,21.3,20.4,20.0,8.0,7.9,2.5,25.3,22.5))

prof <- fread("data/profiles_clean.csv")[lakeid %in% meta$lakeid & !is.na(wtemp)]
prof[, `:=`(date=as.Date(sampledate), year=year4)]
prof <- prof[!(lakeid=="CR" & year %in% c(2012,2013))]
cast <- prof[, .(zbot=max(depth), btm=mean(wtemp[depth >= max(depth) - 1])), by=.(lakeid, date, year)]
cast[, zref := median(zbot), by=lakeid]
cast <- cast[zbot >= 0.9 * zref]
cast[, month := month(date)]
cast <- cast[month >= 4 & month <= 11]
cast[, mlab := factor(month.abb[month], levels=month.abb[4:11])]
cast <- merge(cast, meta, by="lakeid")
ord <- meta[order(region, -zmax), name]
cast[, strip := factor(name, levels=ord)]

## ---- spread heatmap -------------------------------------------------------------------------
sp <- cast[, .(spread=unname(diff(quantile(btm, c(0.05, 0.95)))), n=.N), by=.(name, mlab)]
sp[, lake := factor(name, levels=rev(ord))]
thr <- mean(range(sp$spread))
g2 <- ggplot(sp, aes(mlab, lake, fill=spread)) +
  geom_tile(color="white", linewidth=0.8) +
  geom_text(aes(label=sprintf("%.1f", spread), color=spread > thr), size=3) +
  scale_color_manual(values=c(`TRUE`="white", `FALSE`="grey15"), guide="none") +
  scale_fill_viridis_c(option="mako", direction=-1, name="Spread (°C, 5th–95th percentile)") +
  scale_x_discrete(position="top", expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  labs(x=NULL, y=NULL) +
  theme_minimal(base_size=10.5) +
  theme(panel.grid=element_blank(), axis.text.x=element_text(size=9), legend.position="bottom",
        legend.key.width=unit(1,"cm"), legend.key.height=unit(0.3,"cm"),
        legend.title.position="top", legend.title=element_text(hjust=0.5, size=9.5))
ggsave("figures/figSX_bottom_temp_spread.png", g2, width=6.5, height=4.2, dpi=500, bg="white")
cat("wrote figures/figSX_bottom_temp_spread.png\n")
print(sp[, .(mean_spread=round(mean(spread),1), min_n=min(n)), by=mlab][order(mlab)])
print(dcast(sp, name ~ mlab, value.var="spread")[match(ord, name)])

write_captions <- function(new_caps) {
  path <- "figures/captions.csv"
  old <- if (file.exists(path)) fread(path) else data.table(file=character(), title=character(), caption=character())
  fwrite(rbind(old[!file %in% new_caps$file], new_caps), path)
  cat("wrote", path, "\n")
}
defn <- paste0("Bottom-water temperature is the mean of readings in the deepest 1 m of each profile; profiles ",
  "that did not reach at least 90% of the lake's median profile depth were excluded. April-November, full ",
  "period of record. Crystal Lake 2012-2013 excluded (whole-lake mixing experiment). ")
write_captions(data.table(
  file="figures/figSX_bottom_temp_spread.png",
  title="Spread of bottom-water temperature by month and lake",
  caption=paste0(defn, "Spread is the difference between the 95th and 5th percentiles (degC) across all profiles in each lake and month, independent of the absolute temperature level; darker cells vary more. Lakes ordered north then south by max depth.")
))
