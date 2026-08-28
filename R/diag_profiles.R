## Diagnostics: every individual temp & DO(%sat) profile, per lake x half-month bin, split
## baseline (pre-2016) vs recent (2016-2025) and colored by cast date. Data: EDI 29.
## Same half-month bins as fig12 (Apr 1 - Nov 16), so these diagnostics can be read side-by-side
## with the fig12 heatmap cells they correspond to -- extends the old Sep-only diagnostics to
## the full season.
suppressMessages({library(data.table); library(ggplot2); library(scales); library(patchwork)})

meta <- data.table(
  lakeid=c("TR","BM","CR","SP","AL","TB","CB","ME","MO","FI","WI"),
  name  =c("Trout","Big Muskellunge","Crystal","Sparkling","Allequash","Trout Bog","Crystal Bog",
           "Mendota","Monona","Fish","Wingra"))

prof <- fread("data/profiles.csv")[lakeid %in% meta$lakeid]
prof[, `:=`(date=as.Date(sampledate), doy=yday(as.Date(sampledate)), year=year4)]
prof <- prof[doy>=91 & doy<=319 & !(lakeid=="CR" & year %in% c(2012,2013))]  # Apr 1 - Nov 15

## half-month bins, identical edges to fig12
hedges <- c(91,106,121,136,152,167,182,197,213,228,244,259,274,289,305,320)
hlabs  <- c("apr1-15","apr16-30","may1-15","may16-31","jun1-15","jun16-30","jul1-15","jul16-31",
            "aug1-15","aug16-31","sep1-15","sep16-30","oct1-15","oct16-31","nov1-15")
prof[, bi := findInterval(doy, hedges)]
prof <- prof[bi>=1 & bi<=length(hlabs)]
prof[, binlab := factor(hlabs[bi], levels=hlabs)]

th <- theme_minimal(base_size=11) + theme(
  panel.grid.minor=element_blank(), legend.key.width=unit(0.35,"cm"),
  plot.title=element_text(face="bold",size=9.5), strip.text=element_text(face="bold"))

## one date-gradient panel (baseline OR recent) for one variable
panel <- function(d, var, xlab, era_lab, pal, hline100=FALSE){
  dd <- d[is.finite(get(var))]
  ggplot(dd, aes(x=get(var), y=depth, group=date, color=date)) +
    { if(hline100) geom_vline(xintercept=100, linetype="dotted", color="grey40") } +
    geom_path(linewidth=0.5) + geom_point(size=0.7) +
    scale_y_reverse() +
    scale_color_gradientn(colours=pal, trans="date",
      guide=guide_colorbar(barheight=unit(2.6,"cm"), barwidth=unit(0.3,"cm"))) +
    labs(title=era_lab, x=xlab, y="Depth (m)", color=NULL) + th
}

## baseline = dark purple -> teal-green (viridis), recent = yellow -> dark red (inferno-ish)
pal_base <- c("#440154","#3b528b","#21918c","#5ec962")
pal_rec  <- c("#fde725","#f98e09","#bc3754","#57106e")

plot_lake_bin <- function(lk, bl){
  d <- prof[lakeid==lk & binlab==bl]
  casts <- unique(d[,.(date,year)])
  if(nrow(casts) < 2) return(invisible(NULL))
  d[, era := fifelse(year<2016,"baseline","recent")]
  n_base <- sum(casts$year<2016); n_rec <- sum(casts$year>=2016)
  lakename <- meta[lakeid==lk, name]

  db <- d[era=="baseline"]; dr <- d[era=="recent"]
  tb <- if(nrow(db)>0) panel(db,"wtemp","Temperature (°C)","baseline",pal_base) else plot_spacer()
  tr <- if(nrow(dr)>0) panel(dr,"wtemp","Temperature (°C)","recent", pal_rec)  else plot_spacer()
  ob <- if(nrow(db)>0) panel(db,"o2sat","Dissolved oxygen (% sat)","baseline",pal_base,hline100=TRUE) else plot_spacer()
  orr<- if(nrow(dr)>0) panel(dr,"o2sat","Dissolved oxygen (% sat)","recent", pal_rec, hline100=TRUE) else plot_spacer()

  final <- (tb|tr)/(ob|orr) +
    plot_annotation(
      title=sprintf("%s: every individual profile sampled %s, by cast date", lakename, bl),
      subtitle=sprintf("%d casts in baseline (pre-2016) vs %d in recent (2016-2025). Dotted line = 100%% (equilibrium saturation).", n_base, n_rec),
      theme=theme(plot.title=element_text(face="bold",size=13), plot.subtitle=element_text(size=9.5,color="grey35")))

  outfile <- sprintf("figures/diagnostics/diag_%s_%s_profiles.png", lk, bl)
  ggsave(outfile, final, width=11, height=11, dpi=300, bg="white")
  cat("wrote", outfile, "\n")
}

for(lk in meta$lakeid) for(bl in hlabs) plot_lake_bin(lk, bl)
