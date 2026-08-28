## Fig 11: stratified-season (June-September) water-clarity trends, 7 northern NTL lakes.
## FOUR independently-measured views of the same property, so they cross-validate each other:
##   Kd        -- light-extinction coefficient (m^-1) from PAR profiles;  higher = darker
##   CDOM a440 -- absorption coefficient at 440 nm, optical signature of coloured DOM; higher = darker
##   DOC       -- dissolved organic carbon (mg/L), the lab-measured driver;  higher = darker
##   Secchi    -- Secchi disk depth (m), read by eye;                        lower  = darker
## The Kd axis is REVERSED and the other three are mapped to match, so on this figure DOWN = darker
## for all four lines: a genuinely darkening lake shows all four tracking downward together.
## Kd/a440/DOC are three different instruments (PAR sensor, spectrophotometer, carbon analyser)
## plus a disk read by eye -- agreement across them is strong evidence, disagreement is a flag.
##
## Secchi is recorded twice per visit -- secview (with viewing scope, PRIMARY here) and secnview
## (no scope, reads ~0.5 m shallower, r=0.97). Trusting just one is fragile: at Big Muskellunge they
## disagree in SIGN (secview +0.007 p=0.98 vs secnview -0.116 p=0.059), so a single-column analysis
## can invent a darkening trend that isn't there. A Secchi trend is only drawn solid when BOTH
## columns are significant and agree in sign; otherwise dashed, and not to be read as a trend.
## Southern lakes are excluded: their PAR profiles only start in 2019, too short to trend.
suppressMessages({library(data.table); library(ggplot2); library(ggh4x)})

SEASON  <- c(152, 273)  # Jun 1 - Sep 30, stratified season
MIN_OBS <- 3            # a lake-year is plotted only if >=3 sampling DATES fall in the window

## Annual value = MEDIAN across sampling dates (not mean): robust to the odd bloom-day or
## post-storm cast, consistent with the median-based treatment used in fig09/fig12.
## Counting DATES not rows matters -- colour has 11 wavelength rows per sample and DOC several
## depths per sample, so a raw row count would let one visit satisfy the 3-observation threshold.
##
## EXCEPTION: the >=3-date rule applies to Kd, DOC and Secchi, which are sampled ~4-9 times each
## summer (so the filter drops only the genuinely thin years: 100% of Kd/Secchi and ~95% of DOC
## lake-years pass). CDOM colour is collected as ONE integrated sample per lake per year by
## design -- median 1 date/lake-year, 0% would pass -- so applying the rule there would not
## remove unreliable years, it would remove the variable. Colour is therefore exempt (min_obs=1)
## and its annual point is a single sample; treat it as the noisier of the four series.
annual <- function(d, min_obs=MIN_OBS)
  d[, .(v=median(v), n=uniqueN(date)), by=.(lakeid,year)][n>=min_obs, .(lakeid,year,v)]

meta <- data.table(
  lakeid =c("TR","BM","CR","SP","AL","TB","CB"),
  name   =c("Trout","Big Muskellunge","Crystal","Sparkling","Allequash","Trout Bog","Crystal Bog"),
  trophic=c("oligotrophic","oligotrophic","oligotrophic","oligotrophic","mesotrophic","dystrophic (bog)","dystrophic (bog)"))

## Theil-Sen slope + median intercept + Mann-Kendall p. Returns the intercept too so the drawn
## line IS the reported fit (an earlier version drew an lm smooth but labelled it with a Sen slope).
senfit <- function(y, x){
  ok <- is.finite(y) & is.finite(x); y <- y[ok]; x <- x[ok]; n <- length(y)
  if(n < 8) return(list(slope=NA_real_, int=NA_real_, p=NA_real_, n=n))
  sl <- median(outer(y,y,"-")[lower.tri(diag(n))] / outer(x,x,"-")[lower.tri(diag(n))], na.rm=TRUE)
  list(slope=sl, int=median(y - sl*x),
       p=suppressWarnings(cor.test(x, y, method="kendall")$p.value), n=n)
}
## no p<0.1 tier: a "." marker on a p=0.059 null is what produced a spurious "weak darkening" read
star <- function(p) ifelse(is.na(p),"", ifelse(p<0.001,"***", ifelse(p<0.01,"**", ifelse(p<0.05,"*",""))))

## ---------------- metric 1: Kd (m^-1), slope of ln(fraction surface PAR) vs depth, thru origin ----
p <- fread("data/profiles.csv")[lakeid %in% meta$lakeid & !is.na(frlight) & frlight>0 & depth>0.1 &
                                !(lakeid=="CR" & year4 %in% c(2012,2013))]
p[, `:=`(year=year4, doy=yday(as.Date(sampledate)), date=as.Date(sampledate))]
p <- p[doy>=SEASON[1] & doy<=SEASON[2]]
kd <- p[, { v <- NA_real_
  if(.N>=3){ f <- tryCatch(lm(log(frlight) ~ depth + 0), error=function(e) NULL)
    if(!is.null(f)) v <- as.numeric(-coef(f)[1]) }
  list(v=v) }, by=.(lakeid, date, year)]           # one Kd per cast
kdyr <- annual(kd[is.finite(v) & v>0])

## ---------------- metric 2: CDOM colour, absorption coefficient a440 (m^-1), 1990- ----------------
## CRITICAL: color.csv `value` is raw absorbance, which scales with the spectrophotometer PATH
## LENGTH (`cuvette`, in cm) via Beer-Lambert. The lab switched cuvettes mid-record -- 10 cm through
## the 2000s, then only 5 cm and 1 cm -- so raw `value` falls by ~2x for reasons that have nothing
## to do with the water. Using it raw makes ALL SEVEN lakes look like they are clearing (every
## slope negative, most p<0.001), including the bogs that are unambiguously browning.
## Converting to the Napierian absorption coefficient a440 = 2.303 * A / L(m) removes the artifact
## and flips the bogs to browning, as DOC and Kd independently say. It also puts colour in m^-1,
## the same units as Kd. Flagged values dropped, consistent with the flag handling in fig03/04.
co <- fread("data/color.csv")[lakeid %in% meta$lakeid & is.finite(value) &
                              wavelength>=435 & wavelength<=445 &
                              is.finite(cuvette) & cuvette>0 & (is.na(color_flag) | color_flag=="") &
                              !(lakeid=="CR" & year4 %in% c(2012,2013))]
co[, `:=`(year=year4, date=as.Date(sampledate), doy=yday(as.Date(sampledate)),
          a440 = 2.303 * value * 100 / cuvette)]
## one value per sample date first (median over the 435-445 nm rows), then the annual median
codate <- co[doy>=SEASON[1] & doy<=SEASON[2], .(v=median(a440)), by=.(lakeid,date,year)]
colyr  <- annual(codate, min_obs=1)   # exempt: one integrated sample per year by design

## ---------------- metric 3: DOC (mg/L), surface (0 m exactly) -------------------------------------
## depth==0, NOT depth<=2: the discrete depths NTL samples differ by lake, so a <=2 m filter is not
## a like-for-like comparison -- BM/CR/SP/TR have only 0 m within 2 m (pure surface grab) while
## AL/TB/CB also have a 2 m sample that would get averaged in. That matters most for Crystal Bog,
## only 2.5 m deep, where 2 m is nearly the sediment: including it inflated its DOC trend by ~24%
## (+0.894 -> +0.683 /decade). Pinning to 0 m makes "surface DOC" mean one thing across all lakes
## and costs almost no data (35-38 lake-years either way). Same change applied in fig08.
doc <- fread("data/chem_north.csv")[lakeid %in% meta$lakeid & doc>0 & doc<80 & depth==0 &
                                    (is.na(flagdoc) | flagdoc=="") &
                                    !(lakeid=="CR" & year4 %in% c(2012,2013))]
doc[, `:=`(year=year4, date=as.Date(sampledate), doy=yday(as.Date(sampledate)))]
docdate <- doc[doy>=SEASON[1] & doy<=SEASON[2], .(v=median(doc)), by=.(lakeid,date,year)]
docyr   <- annual(docdate)

## ---------------- metric 4: Secchi (m), BOTH columns, for the robustness test --------------------
sec <- fread("data/secchi.csv")[lakeid %in% meta$lakeid &
                                !(lakeid=="CR" & year4 %in% c(2012,2013))]
sec[, `:=`(year=year4, date=as.Date(sampledate), doy=yday(as.Date(sampledate)))]
sec <- sec[doy>=SEASON[1] & doy<=SEASON[2]]
secyr  <- annual(sec[!is.na(secview)  & secview>0,  .(lakeid,date,year,v=secview)])   # primary
secyr2 <- annual(sec[!is.na(secnview) & secnview>0, .(lakeid,date,year,v=secnview)])  # cross-check

## ---------------- trends -------------------------------------------------------------------------
tr <- function(d, nm) d[, { f<-senfit(v, year); .(metric=nm, slope=f$slope*10, p=f$p, nyr=f$n) }, by=lakeid]
fits <- rbind(tr(kdyr,"Kd"), tr(colyr,"a440"), tr(docyr,"DOC"), tr(secyr,"Secchi"), tr(secyr2,"Secchi_alt"))
fw <- dcast(fits, lakeid ~ metric, value.var=c("slope","p"))
## Secchi counts as robust only if BOTH read methods are significant AND agree in sign
fw[, sec_robust := !is.na(p_Secchi) & !is.na(p_Secchi_alt) &
     p_Secchi<0.05 & p_Secchi_alt<0.05 & sign(slope_Secchi)==sign(slope_Secchi_alt)]
fw <- merge(fw, meta, by="lakeid")
fw[, strip := sprintf(
      "%s (%s)\nKd %+.3f%s   ·   a440 %+.3f%s\nDOC %+.3f%s   ·   Secchi %+.2f%s%s",
      name, trophic, slope_Kd, star(p_Kd), slope_a440, star(p_a440),
      slope_DOC, star(p_DOC), slope_Secchi, star(p_Secchi),
      fifelse(sec_robust, "", " (fragile)"))]
setorder(fw, -slope_Kd)
fw[, strip := factor(strip, levels=strip)]

## ---------------- map all three onto the Kd axis --------------------------------------------------
## Kd and a440 rise as the lake darkens; Secchi falls. Mapping Secchi inversely + reversing the
## primary axis puts "darker" at the bottom for all three, so agreement is read as parallel lines.
rng <- Reduce(function(a,b) merge(a,b,by="lakeid"), list(
  kdyr[,  .(k0=min(v), k1=max(v)), by=lakeid],
  colyr[, .(c0=min(v), c1=max(v)), by=lakeid],
  docyr[, .(d0=min(v), d1=max(v)), by=lakeid],
  secyr[, .(s0=min(v), s1=max(v)), by=lakeid]))
rng[, `:=`(kspan=k1-k0, cspan=c1-c0, dspan=d1-d0, sspan=s1-s0)]

mk <- function(d, nm){ x <- merge(d, rng, by="lakeid")
  x[, yplot := switch(nm,
      Kd     = v,
      a440   = k0 + (v - c0)/cspan * kspan,          # same direction as Kd
      DOC    = k0 + (v - d0)/dspan * kspan,          # same direction as Kd
      Secchi = k0 + (s1 - v)/sspan * kspan)]         # inverted: clear water -> low Kd
  x[, .(lakeid, year, v, yplot, metric=nm)] }
dat <- rbind(mk(kdyr,"Kd"), mk(colyr,"a440"), mk(docyr,"DOC"), mk(secyr,"Secchi"))
dat <- merge(dat, fw[,.(lakeid,strip)], by="lakeid")

## trend segments, fitted in the plotted coordinate so the line matches what is drawn
segs <- dat[, { f<-senfit(yplot, year)
  .(x=min(year), xend=max(year),
    y=f$int+f$slope*min(year), yend=f$int+f$slope*max(year), p=f$p) }, by=.(lakeid,strip,metric)]
segs <- merge(segs, fw[,.(lakeid,sec_robust)], by="lakeid")
## solid = trustworthy: significant, and for Secchi also agreeing across both read methods
segs[, solid := !is.na(p) & p<0.05 & (metric!="Secchi" | sec_robust)]

lv <- c("Kd  (light extinction)","a440  (CDOM colour)","DOC  (mg/L)","Secchi  (disk depth)")
relab <- c(Kd=lv[1], a440=lv[2], DOC=lv[3], Secchi=lv[4])
dat[,  metric := factor(relab[metric], levels=lv)]
segs[, metric := factor(relab[metric], levels=lv)]
pal <- setNames(c("#2b6a3d","#a1622f","#6a4c93","#1f6f9e"), lv)

## per-panel secondary (Secchi) axis; scale_y_reverse => Kd increases downward
rng <- rng[match(fw$lakeid, rng$lakeid)]
yscales <- lapply(seq_len(nrow(rng)), function(i){ r <- rng[i]
  scale_y_reverse(sec.axis = sec_axis(
    transform = ~ r$s1 - (. - r$k0)/r$kspan * r$sspan, name="Secchi depth (m)")) })

g <- ggplot(dat, aes(year, yplot, color=metric)) +
  geom_line(alpha=0.3, linewidth=0.35) +
  geom_point(size=0.8, alpha=0.75) +
  geom_segment(data=segs, aes(x=x, xend=xend, y=y, yend=yend, color=metric, linetype=solid),
               linewidth=0.9, inherit.aes=FALSE) +
  facet_wrap(~strip, scales="free_y", ncol=3) +
  facetted_pos_scales(y = yscales) +
  scale_color_manual(values=pal, name=NULL) +
  scale_linetype_manual(values=c(`TRUE`="solid",`FALSE`="22"), guide="none") +
  labs(title="Open-water water clarity across the 7 northern NTL lakes",
       subtitle=paste0(
         "June–September. Four independently measured views of clarity, all oriented so DOWN = darker water: a truly darkening lake shows all four lines falling together.\n",
         "Each point is the annual MEDIAN across sampling dates. Lake-years with <3 sampling dates are dropped for Kd, DOC and Secchi (sampled 4–9×/summer); CDOM a440 is one integrated sample per year by design and is exempt — the noisiest of the four series.\n",
         "Kd on the left axis (reversed); Secchi (secview) on the right axis; CDOM a440 and DOC rescaled onto the same panel for direction comparison (a440 from 1990, DOC from 1986, others 1981). Crystal 2012–13 excluded (whole-lake mixing experiment).\n",
         "Solid trend = significant (Theil-Sen / Mann-Kendall, p<0.05); dashed = not significant, or — for Secchi — not reproducible across both disk-read methods. Labels: slope/decade (*<0.05 **<0.01 ***<0.001)."),
       x=NULL, y=expression("K"[d]*" (m"^-1*")  \u2014 reversed")) +
  theme_minimal(base_size=11) +
  theme(panel.grid.minor=element_blank(), legend.position="top",
        plot.title=element_text(face="bold",size=13), plot.subtitle=element_text(color="grey40",size=8.3),
        strip.text=element_text(size=8, lineheight=1.05),
        axis.title.y.left=element_text(color=pal[1]), axis.text.y.left=element_text(color=pal[1]),
        axis.title.y.right=element_text(color=pal[4]), axis.text.y.right=element_text(color=pal[4]))
ggsave("figures/fig11_light_attenuation.png", g, width=12.5, height=8.5, dpi=400, bg="white")
cat("wrote figures/fig11_light_attenuation.png\n")

out <- fw[order(-slope_Kd), .(name,
  Kd=round(slope_Kd,4), Kd_p=round(p_Kd,4),
  a440=round(slope_a440,3), a440_p=round(p_a440,4),
  DOC=round(slope_DOC,3), DOC_p=round(p_DOC,4),
  Secchi=round(slope_Secchi,3), Secchi_p=round(p_Secchi,4),
  Secchi_alt=round(slope_Secchi_alt,3), Secchi_alt_p=round(p_Secchi_alt,4), sec_robust)]
print(out)
