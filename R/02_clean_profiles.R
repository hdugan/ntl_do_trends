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
prof[, source := "NTL"]

## --- add Lake Mendota MEMO profiles (EDI knb-lter-ntl.415.3; data/MEMO/) --------------------
## Extra Mendota temperature/DO profiles, 2006-2019, from three files: a YSI sonde on the
## Poseidon (2014-2019; temp + DO) and earlier hand-recorded casts (temp file + DO file,
## YSI MPS 556 / Hydrolab / YSI Pro Plus). Checked against NTL casts on shared days:
## temperature agrees closely (median diff ~0 degC, r = 0.99); DO follows the NTL profile
## shape but calibration offsets vary by day (about +/-1 mg/L), so it is screened:
##   - DO dates whose field notes report meter problems are dropped (DO only; temp kept),
##     and casts noted as "Not at deep hole" are dropped entirely
##   - DO < 0 or > 20 mg/L (incl. -99999 missing codes) and temp <= 0 are set to NA
##   - one cast per day (the earliest); duplicate readings at a depth are averaged
##   - any day with an NTL Mendota cast keeps the NTL cast only (NTL wins)
## MEMO has no o2sat; it is computed here with NTL's own convention, which is exactly
## 1.0613 x 100 x o2 / (standard fresh-water O2 saturation at 1 atm; Benson & Krause form)
## -- checked on every NTL Mendota reading (ratio 1.0613 at all temperatures and decades).
## Set USE_MEMO <- FALSE to rebuild profiles_clean.csv from NTL data alone.
USE_MEMO <- TRUE
if(USE_MEMO){
  mdate <- function(d) d[, date := as.IDate(sprintf("%d-%02d-%02d", Year, Month, Day))]
  first_cast <- function(d){ d[is.na(Hour), Hour := 0L]; d[d[, .I[Hour == min(Hour)], by=date]$V1] }

  ysi <- mdate(fread("data/MEMO/robin_MEMO_YSI-Poseidon_2014-2019.csv"))
  ysi <- first_cast(ysi)[, .(date, depth=as.numeric(Depth.m), wtemp=Temperature.C, o2=Dissolved.Oxygen.mg.L)]

  mt  <- mdate(fread("data/MEMO/robin_MEMO_non-Poseidon_water_temps.csv"))
  mo  <- mdate(fread("data/MEMO/robin_MEMO_non-Poseidon_dissolved_oxygen.csv"))
  offsite  <- mo[grepl("Not at deep hole", Notes.Dissolved.Oxygen, ignore.case=TRUE), unique(date)]
  bad_do   <- mo[grepl(paste("suspicious|unstable|not working|broken|estimated|never stable|out of range",
                             "not stabiliz|sporadic|not properly warmed|wrong pressure|incorrect pressure", sep="|"),
                       Notes.Dissolved.Oxygen, ignore.case=TRUE), unique(date)]
  mt <- first_cast(mt)[, .(wtemp=mean(Water.Temp.C, na.rm=TRUE)), by=.(date, depth=Depth.m)]
  mo <- first_cast(mo[!date %in% bad_do])[, .(o2=mean(Dissolved.Oxygen.mg.L, na.rm=TRUE)), by=.(date, depth=Depth.m)]
  hand <- merge(mt, mo, by=c("date","depth"), all=TRUE)

  memo <- rbind(ysi, hand[!date %in% ysi$date])          # YSI file wins if a day is in both
  memo <- memo[!date %in% offsite]
  memo[!is.finite(o2) | o2 < 0 | o2 > 20, o2 := NA_real_]
  memo[!is.finite(wtemp) | wtemp <= 0, wtemp := NA_real_]
  memo <- memo[is.finite(wtemp) | is.finite(o2)]
  n_shared <- memo[date %in% prof[lakeid == "ME", sampledate], uniqueN(date)]
  memo <- memo[!date %in% prof[lakeid == "ME", sampledate]]   # NTL wins on shared days

  o2_sat_ref <- function(T){ Ts <- log((298.15 - T)/(273.15 + T))
    1.42905 * exp(2.00907 + 3.22014*Ts + 4.0501*Ts^2 + 4.94457*Ts^3 - 0.256847*Ts^4 + 3.88767*Ts^5) }
  memo[, o2sat := round(1.0613 * 100 * o2 / o2_sat_ref(wtemp), 1)]

  memo <- memo[, .(lakeid="ME", year4=year(date), sampledate=date, depth, rep=1L, sta=NA_integer_,
                   event="MEMO", wtemp, o2, o2sat, source="MEMO")]
  cat(sprintf("MEMO: added %d casts (%d with DO); dropped %d off-site casts, DO on %d dates with meter-problem notes, %d days already sampled by NTL\n",
              uniqueN(memo$sampledate), memo[is.finite(o2), uniqueN(sampledate)], length(offsite), length(bad_do), n_shared))
  prof <- rbind(prof, memo, fill=TRUE)
  setorder(prof, lakeid, sampledate, depth)
}

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
