## Clean unrealistic bottom-point jumps out of data/profiles.csv -> data/profiles_clean.csv.
## Raw data/profiles.csv (as pulled by 00_pull_data.R) is left untouched; the fig*.R scripts
## read the cleaned file instead.
## Spot-checked in the diagnostic plots (01_plot_profiles.R): a handful of casts have a single
## deepest reading that jumps sharply UP from the next-shallowest reading -- e.g. ME 2022-07-13
## (o2sat 8.3 -> 102.4 %sat at the bottom) and MO 2001-10-22 (wtemp 12.7 -> 21.7 degC at the
## bottom). Physically implausible: temperature doesn't jump up at the very bottom of a
## stratified lake, and DO doesn't suddenly reaerate at the deepest point. Reads as a sensor
## glitch on the last reading of the cast, so only that single value is nulled -- the rest of
## the cast (and the other variable at that same depth) is left untouched.
## Restricted to lakes with zmax > 10 m: in shallow lakes/bogs (Allequash, Trout Bog, Crystal
## Bog, Wingra) a jump like this over 1-2 m of depth is plausible real biogeochemistry, not a
## sensor glitch -- only in a deep, stratified lake is a jump at the very bottom implausible.
## Also manually removes one whole bad-sensor DO profile (Sparkling 2004-08-16), and all
## temp/DO for Crystal Lake in 2012-2013 (whole-lake mixing experiment) -- see the "manual
## removal" blocks below.
suppressMessages(library(data.table))

TEMP_JUMP_MAX <- 5   # degC: max allowed increase at the single deepest reading vs the next-shallowest
DO_JUMP_MAX   <- 20  # %sat: max allowed increase at the single deepest reading vs the next-shallowest
ZMAX_MIN      <- 10  # m: only apply the check to lakes deeper than this

## max depth (m) of each of the 11 NTL primary study lakes, same values used in fig12
zmax <- data.table(
  lakeid=c("TR","BM","CR","SP","AL","TB","CB","ME","MO","FI","WI"),
  zmax  =c(35.7,21.3,20.4,20.0,8.0,7.9,2.5,25.3,22.5,18.9,4.0))
deep_lakes <- zmax[zmax > ZMAX_MIN, lakeid]

prof <- fread("data/profiles.csv")

## for one variable: within each cast, order by depth and flag the deepest valid reading if it
## jumps up by more than `thresh` from the next-shallowest valid reading
flag_bottom_jump <- function(dt, var, thresh){
  d <- dt[is.finite(get(var)), .(lakeid, sampledate, depth, val=get(var))]
  setorder(d, lakeid, sampledate, depth)
  d[, `:=`(prev_val = shift(val), prev_depth = shift(depth)), by=.(lakeid, sampledate)]
  d[, is_deepest := seq_len(.N) == .N, by=.(lakeid, sampledate)]
  bad <- d[is_deepest == TRUE & is.finite(prev_val) & (val - prev_val) > thresh]
  bad[, .(lakeid, sampledate, variable = var, depth, value = val,
          prev_depth, prev_value = prev_val, jump = val - prev_val)]
}

removed <- rbind(
  flag_bottom_jump(prof[lakeid %in% deep_lakes], "wtemp", TEMP_JUMP_MAX),
  flag_bottom_jump(prof[lakeid %in% deep_lakes], "o2sat", DO_JUMP_MAX)
)
setorder(removed, variable, lakeid, sampledate)
removed[, reason := "auto_bottom_jump"]

## o2 (mg/L) and o2sat are the same underlying DO measurement in different units. Checked every
## flagged o2sat row: o2 jumps anomalously there too (same corrupted reading), so it's nulled
## alongside o2sat to keep the two units consistent; logged here for the record.
removed[, o2_also_nulled := NA_real_]
setkey(prof, lakeid, sampledate, depth)
for(i in which(removed$variable == "o2sat")){
  removed$o2_also_nulled[i] <- prof[.(removed$lakeid[i], removed$sampledate[i], removed$depth[i]), o2]
}

## --- manual removal: the whole Sparkling 2004-08-16 DO profile -----------------------------
## Visually odd in 01_plot_profiles.R (diag_SP_aug16-31_profiles.png): every other August cast
## at Sparkling stays >100% o2sat out past 4 m then crashes sharply into the metalimnion, but
## this one starts declining by 4 m (102.2 -> 78.4 %sat, vs a next-largest surface-to-4m drop of
## just 3.0 pts across the whole record) and instead of a sharp crash sits unrealistically low
## and flat (43-67 %sat) through 5-13 m -- a bad DO sensor for the whole cast, not one point.
## Temperature that day is unremarkable, so only o2 / o2sat are removed; wtemp is kept.
manual_lake <- "SP"; manual_date <- as.IDate("2004-08-16")
manual <- prof[lakeid == manual_lake & sampledate == manual_date & is.finite(o2sat),
  .(lakeid, sampledate, variable = "o2sat", depth, value = o2sat,
    prev_depth = NA_real_, prev_value = NA_real_, jump = NA_real_,
    reason = "manual_full_profile", o2_also_nulled = o2)]
removed <- rbind(removed, manual)

## --- manual removal: Crystal Lake 2012-2013, whole-lake mixing experiment ------------------
## Crystal was experimentally whole-lake mixed in 2012-2013 (the source project's fig08/09/11/12
## scripts already exclude these two years for the same reason). Temperature and DO during this
## period reflect the manipulation, not natural dynamics, so all wtemp and o2/o2sat for CR in
## 2012-2013 are removed here too -- an experimental artifact, not a sensor error, but treated
## the same way so nothing downstream has to special-case it.
cr_bad <- prof[lakeid == "CR" & year4 %in% c(2012, 2013)]
manual2 <- rbind(
  cr_bad[is.finite(wtemp), .(lakeid, sampledate, variable = "wtemp", depth, value = wtemp,
    prev_depth = NA_real_, prev_value = NA_real_, jump = NA_real_,
    reason = "manual_crystal_mixing_experiment", o2_also_nulled = NA_real_)],
  cr_bad[is.finite(o2sat), .(lakeid, sampledate, variable = "o2sat", depth, value = o2sat,
    prev_depth = NA_real_, prev_value = NA_real_, jump = NA_real_,
    reason = "manual_crystal_mixing_experiment", o2_also_nulled = o2)]
)
removed <- rbind(removed, manual2)

cat("=== Points removed, by reason x lake x variable ===\n")
print(removed[, .N, by=.(reason, lakeid, variable)][order(reason, lakeid, variable)], nrows=100)
cat(sprintf("\n(%d total; full detail in data/profiles_removed_points.csv)\n", nrow(removed)))

## null out the flagged value (row and other variables at that depth kept), via a keyed join
## rather than a per-row scan -- the Crystal 2012-13 removal alone is ~1800 rows
setkey(prof, lakeid, sampledate, depth)
for(v in unique(removed$variable)){
  hits <- removed[variable == v, .(lakeid, sampledate, depth)]
  prof[hits, (v) := NA_real_]
  if(v == "o2sat") prof[hits, o2 := NA_real_]
}

fwrite(prof, "data/profiles_clean.csv")
fwrite(removed, "data/profiles_removed_points.csv")
cat(sprintf("\nRemoved %d point(s). Wrote data/profiles_clean.csv; logged to data/profiles_removed_points.csv\n",
            nrow(removed)))
