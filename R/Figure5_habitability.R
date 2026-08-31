## Fig 5: habitability status by depth window x month, all 11 NTL lakes. Data: EDI 29.
## For each lake x depth bin x month (Apr-Nov), compares median DO (mg/L) in two eras -- baseline
## (pre-2016) vs recent (2016-2025) -- against the 3 mg/L cold/coolwater-fish habitat threshold,
## and classifies the transition into one of four states:
##   remained  < 3 mg/L : chronically hypoxic, unhabitable in both eras
##   became    < 3 mg/L : newly hypoxic -- lost habitat
##   became   >= 3 mg/L : newly oxygenated -- gained habitat
##   remained >= 3 mg/L : stayed habitable in both eras
## Drawn as a depth x month heatmap per lake (facet_wrap), so a lake's whole seasonal cycle of
## habitat change is visible at once, not just its August snapshot.
suppressMessages({library(data.table); library(ggplot2); library(scales)})

meta <- data.table(
  lakeid=c("TR","BM","CR","SP","AL","TB","CB","ME","MO","FI","WI"),
  name  =c("Trout","Big Muskie","Crystal","Sparkling","Allequash","Trout Bog","Crystal Bog",
           "Mendota","Monona","Fish","Wingra"),
  region=c(rep("Northern",7), rep("Southern",4)),
  zmax  =c(35.7,21.3,20.4,20.0,8.0,7.9,2.5,25.3,22.5,18.9,4.0))
zmax_lookup <- setNames(meta$zmax, meta$lakeid)

HYPOXIA_THRESHOLD <- 3   # mg/L, cold/coolwater-fish habitat cutoff
MIN_YEARS <- 4           # a depth x month cell needs data from at least this many distinct years,
                          # in BOTH eras, to be classified; otherwise left blank (insufficient evidence)
ERA_BREAK <- 2016        # same baseline/recent split used elsewhere in this repo

prof <- fread("data/profiles_clean.csv")[lakeid %in% meta$lakeid & !is.na(o2)]
prof[, `:=`(date=as.Date(sampledate), month=month(as.Date(sampledate)), year=year4)]
prof <- prof[month>=4 & month<=11 & !(lakeid=="CR" & year %in% c(2012,2013))]  # Apr-Nov, same season as Figure1

## interpolate each cast to a depth grid: 0.5 m for lakes <10 m zmax, 1 m for deeper lakes --
## same rule used in Figure1_rate_of_change.R, so a "depth window" here means the same thing it
## does there
g <- prof[, { o<-order(depth); dd<-depth[o]; yy<-o2[o]
  dstep <- if(zmax_lookup[[lakeid]] < 10) 0.5 else 1
  if(length(dd)>=3 && diff(range(dd))>=1){ zg<-seq(0,max(dd),dstep); zg<-zg[zg>=min(dd)]
    .(dbin=zg, dstep=dstep, val=approx(dd,yy,zg,rule=2)$y) } else .(dbin=numeric(0),dstep=numeric(0),val=numeric(0)) },
  by=.(lakeid, date, year, month)]

g[, era := fifelse(year < ERA_BREAK, "baseline", "recent")]

## per lake x depth bin x month x era: median DO and how many distinct years back it
agg <- g[, .(med=median(val), n_yr=uniqueN(year)), by=.(lakeid, dbin, dstep, month, era)]
agg[, below := med < HYPOXIA_THRESHOLD]

wide <- dcast(agg, lakeid+dbin+dstep+month ~ era, value.var=c("below","n_yr"))
wide <- wide[!is.na(below_baseline) & !is.na(below_recent) &
             n_yr_baseline>=MIN_YEARS & n_yr_recent>=MIN_YEARS]

state_lab <- c("remained_below"="Remained < 3 mg/L (chronically hypoxic)",
               "became_below"  ="Became < 3 mg/L (newly hypoxic)",
               "became_above"  ="Became ≥ 3 mg/L (newly oxygenated)",
               "remained_above"="Remained ≥ 3 mg/L (oxygenated)")
wide[, state := fifelse(below_baseline & below_recent, "remained_below",
              fifelse(!below_baseline & below_recent, "became_below",
              fifelse(below_baseline & !below_recent, "became_above", "remained_above")))]
wide[, state := factor(state_lab[state], levels=state_lab)]

wide <- merge(wide, meta, by="lakeid")
ord <- meta[order(region, -zmax), name]
wide[, strip := factor(name, levels=ord)]

mlabs <- c("Apr","May","Jun","Jul","Aug","Sep","Oct","Nov")
wide[, mlab := factor(mlabs[month-3], levels=mlabs)]

pal <- c("Remained < 3 mg/L (chronically hypoxic)"="#7f0000",
         "Became < 3 mg/L (newly hypoxic)"        ="#e66101",
         "Became ≥ 3 mg/L (newly oxygenated)" ="#1a9850",
         "Remained ≥ 3 mg/L (oxygenated)"     ="#2166ac")

g_plot <- ggplot(wide, aes(x=mlab, y=dbin, fill=state)) +
  geom_tile(aes(height=dstep), width=0.98) +
  scale_y_reverse(expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  scale_fill_manual(values=pal, name=NULL, drop=FALSE, na.translate=FALSE) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
  facet_wrap(~strip, ncol=4, scales="free_y") +
  labs(x=NULL, y="Depth (m)") +
  theme_minimal(base_size=10.5) +
  theme(panel.grid=element_blank(), legend.position="bottom", legend.key.size=unit(0.4,"cm"),
        legend.text=element_text(size=8), axis.text.x=element_text(size=7, angle=45, hjust=1),
        strip.text=element_text(face="bold", size=9.5), panel.spacing=unit(0.6,"lines"))

ggsave("figures/fig05_habitability.png", g_plot, width=9.5, height=7.5, dpi=500, bg="white")
cat("wrote figures/fig05_habitability.png\n")

## Title + caption off the PNG, into the shared captions.csv (same convention as the other fig*.R)
write_captions <- function(new_caps) {
  path <- "figures/captions.csv"
  old <- if (file.exists(path)) fread(path) else data.table(file=character(), title=character(), caption=character())
  fwrite(rbind(old[!file %in% new_caps$file], new_caps), path)
  cat("wrote", path, "\n")
}
write_captions(data.table(
  file="figures/fig05_habitability.png",
  title="Habitability by depth and month: dissolved oxygen, baseline vs. recent",
  caption=paste0(
    "For each lake, depth window (0.5 m bins for lakes under 10 m zmax, 1 m otherwise -- same ",
    "rule as Figure1_rate_of_change.R), and month (April-November), median dissolved oxygen (mg/L) ",
    "is compared between a baseline era (pre-2016) and a recent era (2016-2025), against a 3 mg/L ",
    "cold/coolwater-fish habitat threshold. Dark red = remained below 3 mg/L in both eras ",
    "(chronically hypoxic); orange = crossed from at-or-above to below 3 mg/L (newly hypoxic, lost ",
    "habitat); green = crossed from below to at-or-above 3 mg/L (newly oxygenated, gained habitat); ",
    "blue = remained at or above 3 mg/L in both eras. A cell is left blank if either era has fewer ",
    "than 4 distinct years of data there. Panels ordered north then south by max depth (zmax); each ",
    "lake's own depth axis runs only to its sampled depth. Crystal Lake 2012-2013 excluded ",
    "(whole-lake mixing experiment)."
  )
))
