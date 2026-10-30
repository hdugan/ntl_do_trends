## FigureSX_surface_temp.R: range of surface water temperature by month, by lake (boxplots).
## Surface temperature = mean wtemp over depths <= 1 m for each cast (one value per profile, so
## casts with several near-surface readings aren't over-weighted). April-November, same window and
## Crystal Lake 2012-2013 exclusion as the other figures; Fish Lake and Wingra are not shown.
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
cast[, mlab := factor(month.abb[month], levels=month.abb[4:11])]
cast <- merge(cast, meta, by="lakeid")
cast[, strip := factor(name, levels=meta[order(region, -zmax), name])]

## Ridgelines: one density curve per month, filled by temperature along the x-axis, so each lake's
## seasonal march (April cold -> July/August peak -> November cold) and how wide each month's
## range is are both visible at a glance. Median marked with a tick; Apr at the top, Nov at the bottom.
suppressMessages(library(ggridges))
cast[, mrev := factor(mlab, levels=rev(levels(mlab)))]
g <- ggplot(cast, aes(x=sst, y=mlab, fill=after_stat(x))) +
  geom_density_ridges_gradient(scale=1.6, rel_min_height=0.01, bandwidth=1.1, linewidth=0.25, color="grey20",
                               quantile_lines=TRUE, quantiles=2, vline_size=0.5, vline_color="white") +
  scale_fill_viridis_c(option="plasma", name="Surface\ntemperature (\u00b0C)", limits=c(0,30), oob=scales::squish) +
  scale_y_discrete(limits=rev(levels(cast$mlab)), expand=expansion(add=c(0.2, 1.4))) +
  scale_x_continuous(limits=c(-2, 32), breaks=seq(0,30,10), expand=c(0,0)) +
  facet_wrap(~strip, ncol=3) +
  labs(x="Surface water temperature (\u00b0C)", y=NULL) +
  theme_minimal(base_size=10.5) +
  theme(panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(),
        strip.text=element_text(face="bold", size=9.5), axis.text.y=element_text(size=7.5),
        legend.position="bottom", legend.key.width=unit(1.6,"cm"), legend.key.height=unit(0.3,"cm"),
        panel.spacing=unit(0.6,"lines"))

ggsave("figures/figSX_surface_temp.png", g, width=9.5, height=7.5, dpi=500, bg="white")
cat("wrote figures/figSX_surface_temp.png\n")

write_captions <- function(new_caps) {
  path <- "figures/captions.csv"
  old <- if (file.exists(path)) fread(path) else data.table(file=character(), title=character(), caption=character())
  fwrite(rbind(old[!file %in% new_caps$file], new_caps), path)
  cat("wrote", path, "\n")
}
write_captions(data.table(
  file="figures/figSX_surface_temp.png",
  title="Surface water temperature by month and lake",
  caption=paste0(
    "Distribution of surface water temperature (mean of readings at 0-1 m depth for each profile) by ",
    "month, April-November, full period of record, shown as ridgeline density curves coloured by ",
    "temperature; the white line in each curve is the monthly median. Each observation is one profile. Panels ordered north then south by max depth (zmax). Crystal Lake 2012-2013 excluded ",
    "(whole-lake mixing experiment)."
  )
))
