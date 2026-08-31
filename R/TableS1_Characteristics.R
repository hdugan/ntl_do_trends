## TableS1_Characteristics.R: pull static lake-characteristic data (location, morphometry,
## trophic state) for the 11 NTL primary study lakes and write a formatted LaTeX table.
## Data: EDI 434 (https://portal.edirepository.org/nis/mapbrowse?scope=knb-lter-ntl&identifier=434).
## Independent of the other R scripts here -- these are static lake descriptors, not the
## time-series data (data/profiles*.csv etc.) the fig*.R scripts analyze.
suppressMessages({library(EDIutils); library(readr); library(data.table)})
dir.create("data", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

## --- Core Research Lakes Information -> data/lake_characteristics.csv ---
## Lottig, N.R. and H.A. Dugan. 2024. North Temperate Lakes-LTER Core Research Lakes
## Information ver 1. Environmental Data Initiative.
## https://doi.org/10.6073/pasta/b9080c962f552029ee2b43aec1410328. Accessed 2026-08-31.
revision <- list_data_package_revisions(
  scope = "knb-lter-ntl",
  identifier = "434",
  filter = "newest"
)
package_id <- paste("knb-lter-ntl", "434", revision, sep = ".")
res <- read_data_entity_names(package_id)
raw <- read_data_entity(package_id, entityId = res$entityId[1])
chars <- readr::read_csv(file = raw, show_col_types = FALSE)
readr::write_csv(chars, "data/lake_characteristics.csv")
cat("wrote data/lake_characteristics.csv (", package_id, ")\n")

## --- Restrict/order to the 11 primary study lakes, same meta + ordering convention (region,
## then descending max depth) used in Figure1_rate_of_change.R and Figure2_clarity_trends.R ---
meta <- data.table(
  lakeid = c("TR", "BM", "CR", "SP", "AL", "TB", "CB", "ME", "MO", "FI", "WI"),
  waterbody_name = c(
    "Trout Lake", "Big Muskellunge Lake", "Crystal Lake", "Sparkling Lake", "Allequash Lake",
    "Trout Bog", "Crystal Bog", "Lake Mendota", "Lake Monona", "Fish Lake", "Lake Wingra"
  ),
  region = c(rep("Northern", 7), rep("Southern", 4))
)

setDT(chars)
lakes <- merge(meta, chars, by = "waterbody_name", sort = FALSE)
if (nrow(lakes) != 11) stop("Expected all 11 study lakes to match EDI 434 by name; got ", nrow(lakes))
setorder(lakes, region, -max_depth_m)

## --- Build the LaTeX table ---
fmt1 <- function(x) ifelse(is.na(x), "--", formatC(x, format = "f", digits = 1))
fmt2 <- function(x) ifelse(is.na(x), "--", formatC(x, format = "f", digits = 2))
fmt_int <- function(x) ifelse(is.na(x), "--", formatC(round(x), format = "d", big.mark = ","))

lines <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\footnotesize",
  "\\caption{Physical and limnological characteristics of the 11 NTL-LTER primary study lakes \\cite{ntllter_lakes434}.}",
  "\\label{tab:lake_characteristics}",
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{lccrrll}",
  "\\toprule",
  "Lake & Lat (\\textdegree N) & Lon (\\textdegree W) & Max depth (m) & Area (ha) & Drainage type & Trophic state \\\\",
  "\\midrule"
)
region_prev <- NULL
for (i in seq_len(nrow(lakes))) {
  r <- lakes[i]
  if (!is.null(region_prev) && r$region != region_prev) lines <- c(lines, "\\hline")
  lines <- c(lines, sprintf(
    "%s & %s & %s & %s & %s & %s & %s \\\\",
    r$waterbody_name,
    fmt2(r$waterbody_lat_decdeg), fmt2(abs(r$waterbody_lon_decdeg)),
    fmt1(r$max_depth_m),
    fmt_int(r$waterbody_area_ha),
    tools::toTitleCase(r$drainage_type), tools::toTitleCase(r$trophic_state)
  ))
  region_prev <- r$region
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}%", "}", "\\end{table}")

writeLines(lines, "tables/TableS1_Characteristics.tex")
cat("wrote tables/TableS1_Characteristics.tex\n")
