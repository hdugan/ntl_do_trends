## FigureSX_surface_temp_spread.R: how much surface temperature VARIES by month and lake,
## ignoring its absolute level. Spread = 5th-95th percentile range (degC) of surface temperature
## (mean of 0-1 m readings per profile) across all profiles in that lake x month, full record.
## Same cast definition, window and exclusions as FigureSX_surface_temp.R.
suppressMessages({library(data.table); library(ggplot2)})

meta <- data.table(
  lakeid=c("TR","BM","CR","SP","AL","TB","CB","ME","MO"),
  name  =c("Trout Lake","Big Muskellunge Lake","Crystal Lake","Sparkling Lake","Allequash Lake",
           "Trout Bog","Crystal Bog","Lake Mendota","Lake Monona"),
  region=c(rep("Northern",7), rep("Southern",2)),
  zmax  =c(35.7,21.3,20.4,20.0,8.0,7.9,2.5,25.3,22.5))

prof <- fread("data/profiles_clean.csv")[lakeid %in% meta$lakeid & !is.na(wtemp) & depth <= 1]
prof[, `:=`(date=as.Date(sampledate), year=year4)]
prof <- prof[!(lakeid=="CR" & year %in% c(2012,2013))]
cast <- prof[, .(sst=mean(wtemp)), by=.(lakeid, date, year)]
cast[, month := month(date)]
cast <- cast[month >= 4 & month <= 11]

sp <- cast[, .(spread=unname(diff(quantile(sst, c(0.05, 0.95)))), n=.N), by=.(lakeid, month)]
sp <- merge(sp, meta, by="lakeid")
sp[, mlab := factor(month.abb[month], levels=month.abb[4:11])]
sp[, lake := factor(name, levels=rev(meta[order(region, -zmax), name]))]

g <- ggplot(sp, aes(mlab, lake, fill=spread)) +
  geom_tile(color="white", linewidth=0.8) +
  geom_text(aes(label=sprintf("%.1f", spread), color=spread > 9.5), size=3) +
  scale_color_manual(values=c(`TRUE`="white", `FALSE`="grey15"), guide="none") +
  scale_fill_viridis_c(option="mako", direction=-1, name="Spread (°C, 5th–95th percentile)") +
  scale_x_discrete(position="top", expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
  labs(x=NULL, y=NULL) +
  theme_minimal(base_size=10.5) +
  theme(panel.grid=element_blank(), axis.text.x=element_text(size=9),
        legend.position="bottom", legend.key.width=unit(1,"cm"), legend.title.position="top", legend.title=element_text(hjust=0.5, size=9.5), legend.key.height=unit(0.3,"cm"))

ggsave("figures/figSX_surface_temp_spread.png", g, width=6.5, height=4.2, dpi=500, bg="white")
cat("wrote figures/figSX_surface_temp_spread.png\n")
print(sp[, .(mean_spread=round(mean(spread),1)), by=mlab][order(mlab)])

write_captions <- function(new_caps) {
  path <- "figures/captions.csv"
  old <- if (file.exists(path)) fread(path) else data.table(file=character(), title=character(), caption=character())
  fwrite(rbind(old[!file %in% new_caps$file], new_caps), path)
  cat("wrote", path, "\n")
}
write_captions(data.table(
  file="figures/figSX_surface_temp_spread.png",
  title="Spread of surface water temperature by month and lake",
  caption=paste0(
    "Range of surface water temperature (mean of 0-1 m readings per profile) across all profiles in ",
    "each lake and month, April-November, full period of record. Spread is the difference between the ",
    "95th and 5th percentiles (degC) and is independent of the absolute temperature level; darker ",
    "cells vary more. Lakes ordered north then south by max depth. Crystal Lake 2012-2013 excluded ",
    "(whole-lake mixing experiment)."
  )
))
