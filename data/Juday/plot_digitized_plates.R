## Re-draw the digitized 1906 Lake Mendota plates (Birge & Juday 1911, Bulletin XXII,
## Plates I and II) in the same layout as the originals, for checking the digitization
## by eye against Birge_Juday_1911_tempME.png / Birge_Juday_1911_DOME.png.
## Open squares = points flagged "low" confidence; open circles = everything else.
## Run from the repo root: Rscript data/Juday/plot_digitized_plates.R
suppressMessages({library(data.table); library(ggplot2); library(ggrepel)})

dir <- "data/Juday"
month_starts <- c(Apr=91, May=121, Jun=152, Jul=182, Aug=213, Sep=244, Oct=274, Nov=305, Dec=335)

plate <- function(d, yvar, ylim, ybreaks, xmax, label_day, ylab, title, file) {
  d <- copy(d)
  d[, y := get(yvar)]
  d[, low := confidence == "low"]
  d[, dlab := factor(depth_m)]
  setorder(d, depth_m, yday)
  ms <- month_starts[month_starts <= ceiling(xmax)]
  mids <- (head(ms, -1) + tail(ms, -1)) / 2
  # one label per curve, at the cast nearest label_day
  labs <- d[, .SD[which.min(abs(yday - label_day))], by=depth_m]

  g <- ggplot(d, aes(yday, y, group=dlab)) +
    geom_hline(yintercept=ybreaks, color="grey25", linewidth=0.5) +
    geom_vline(xintercept=ms, linetype="22", color="grey25", linewidth=0.4) +
    geom_path(aes(linetype=dlab, linewidth=dlab)) +
    geom_point(data=d[low==FALSE], shape=21, fill="white", size=1.4, stroke=0.5) +
    geom_point(data=d[low==TRUE], shape=22, fill="white", color="red3", size=1.6, stroke=0.6) +
    geom_text_repel(data=labs, aes(label=depth_m), size=3.2, fontface="bold",
                    min.segment.length=0, segment.color="grey50", seed=1, box.padding=0.4) +
    scale_x_continuous(limits=c(91, xmax), breaks=mids, labels=toupper(names(mids)),
                       position="top", expand=c(0.01, 0)) +
    scale_y_continuous(limits=ylim, breaks=ybreaks, expand=c(0, 0)) +
    scale_linetype_manual(values=rep(c("solid","42","12","4212","22","1242"), 3), guide="none") +
    scale_linewidth_manual(values=rep(c(1.1,0.45,0.45,1.1,0.45,0.45,1.1,0.45,0.45,1.1,0.45), 2), guide="none") +
    labs(x=NULL, y=ylab, title=title,
         caption="Digitized re-draw for checking. Open squares (red) = low-confidence points.") +
    theme_minimal(base_size=11) +
    theme(panel.grid=element_blank(), plot.background=element_rect(fill="#fbf3dc", color=NA),
          panel.border=element_rect(fill=NA, color="grey25", linewidth=0.6),
          plot.title=element_text(size=10.5, hjust=0.5), plot.caption=element_text(color="red3", size=8),
          axis.text.x=element_text(face="bold", size=9))
  ggsave(file.path(dir, file), g, width=10, height=6.2, dpi=200)
  cat("wrote", file.path(dir, file), "\n")
}

temp <- fread(file.path(dir, "Mendota_1906_temp.csv"))
plate(temp, "temp_C", c(7, 26), seq(8, 24, 2), 306, 205,
      "Temperature (°C)",
      "Plate I (re-drawn) — Temperature of the water at different depths in Lake Mendota in 1906",
      "check_Plate_I_temp_redrawn.png")

do <- fread(file.path(dir, "Mendota_1906_DO_plate.csv"))
plate(do, "do_mLL", c(0, 8.05), 0:8, 335, 177,
      "Dissolved oxygen (cc per liter)",
      "Plate II (re-drawn) — Dissolved oxygen at different depths in Lake Mendota in 1906",
      "check_Plate_II_DO_redrawn.png")
