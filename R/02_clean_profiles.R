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

## o2 (mg/L) and o2sat are the same underlying DO measurement in different units. Checked every
## flagged o2sat row: o2 jumps anomalously there too (same corrupted reading), so it's nulled
## alongside o2sat to keep the two units consistent; logged here for the record.
removed[, o2_also_nulled := NA_real_]
setkey(prof, lakeid, sampledate, depth)
for(i in which(removed$variable == "o2sat")){
  removed$o2_also_nulled[i] <- prof[.(removed$lakeid[i], removed$sampledate[i], removed$depth[i]), o2]
}

cat("=== Unrealistic bottom-point jumps removed ===\n")
print(removed, nrows = nrow(removed))

## null out the flagged value (row and other variables at that depth kept)
for(i in seq_len(nrow(removed))){
  v <- removed$variable[i]
  cond <- prof$lakeid == removed$lakeid[i] & prof$sampledate == removed$sampledate[i] &
          prof$depth == removed$depth[i]
  prof[cond, (v) := NA_real_]
  if(v == "o2sat") prof[cond, o2 := NA_real_]
}

fwrite(prof, "data/profiles_clean.csv")
fwrite(removed, "data/profiles_removed_points.csv")
cat(sprintf("\nRemoved %d point(s). Wrote data/profiles_clean.csv; logged to data/profiles_removed_points.csv\n",
            nrow(removed)))
