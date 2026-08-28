## Clean unrealistic bottom-point jumps out of data/profiles.csv (in place).
## Spot-checked in the diagnostic plots (01_plot_profiles.R): a handful of casts have a single
## deepest reading that jumps sharply UP from the next-shallowest reading -- e.g. ME 2022-07-13
## (o2sat 8.3 -> 102.4 %sat at the bottom) and MO 2001-10-22 (wtemp 12.7 -> 21.7 degC at the
## bottom). Physically implausible: temperature doesn't jump up at the very bottom of a
## stratified lake, and DO doesn't suddenly reaerate at the deepest point. Reads as a sensor
## glitch on the last reading of the cast, so only that single value is nulled -- the rest of
## the cast (and the other variable at that same depth) is left untouched.
suppressMessages(library(data.table))

TEMP_JUMP_MAX <- 5   # degC: max allowed increase at the single deepest reading vs the next-shallowest
DO_JUMP_MAX   <- 20  # %sat: max allowed increase at the single deepest reading vs the next-shallowest

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
  flag_bottom_jump(prof, "wtemp", TEMP_JUMP_MAX),
  flag_bottom_jump(prof, "o2sat", DO_JUMP_MAX)
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

fwrite(prof, "data/profiles.csv")
fwrite(removed, "data/profiles_removed_points.csv")
cat(sprintf("\nRemoved %d point(s). Overwrote data/profiles.csv; logged to data/profiles_removed_points.csv\n",
            nrow(removed)))
