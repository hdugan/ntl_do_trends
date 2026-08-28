## Pull raw data from EDI for the figures in this repo. Run once before the fig*.R / 01_plot_profiles.R
## scripts (they all read data/<file>.csv directly). Re-run any time to refresh to the newest
## published revision of each package.
suppressMessages({library(EDIutils); library(readr)})
dir.create("data", showWarnings = FALSE)

## --- Physical Limnology of Primary Study Lakes (temp, DO, DO %sat, light) -> data/profiles.csv ---
## Magnuson, J.J., S.R. Carpenter, and E.H. Stanley. 2026. North Temperate Lakes LTER: Physical
## Limnology of Primary Study Lakes 1981 - current ver 38. Environmental Data Initiative.
## https://doi.org/10.6073/pasta/eeaa2a029eda1683f004d06da5d8d808. Accessed 2026-08-28.
revision <- list_data_package_revisions(
  scope = "knb-lter-ntl",
  identifier = "29",
  filter = "newest"
)
package_id <- paste("knb-lter-ntl", "29", revision, sep = ".")
res <- read_data_entity_names(package_id)
raw <- read_data_entity(package_id, entityId = res$entityId[1])
profiles <- readr::read_csv(file = raw, show_col_types = FALSE)
readr::write_csv(profiles, "data/profiles.csv")
cat("wrote data/profiles.csv (", package_id, ")\n")

## --- Chemical Limnology of Primary Study Lakes: Nutrients, pH and Carbon (DOC) -> data/chem_north.csv ---
## Magnuson, J.J., S.R. Carpenter, and E.H. Stanley. 2026. North Temperate Lakes LTER: Chemical
## Limnology of Primary Study Lakes: Nutrients, pH and Carbon 1981 - current ver 65. Environmental
## Data Initiative. https://doi.org/10.6073/pasta/ffa652bb3cb1a2c82c674c3e5354a5b2. Accessed 2026-08-28.
revision <- list_data_package_revisions(
  scope = "knb-lter-ntl",
  identifier = "1",
  filter = "newest"
)
package_id <- paste("knb-lter-ntl", "1", revision, sep = ".")
res <- read_data_entity_names(package_id)
raw <- read_data_entity(package_id, entityId = res$entityId[1])
chem_north <- readr::read_csv(file = raw, show_col_types = FALSE)
readr::write_csv(chem_north, "data/chem_north.csv")
cat("wrote data/chem_north.csv (", package_id, ")\n")

## --- Secchi Disk Depth; Other Auxiliary Base Crew Sample Data -> data/secchi.csv ---
## Magnuson, J.J., S.R. Carpenter, and E.H. Stanley. 2026. North Temperate Lakes LTER: Secchi Disk
## Depth; Other Auxiliary Base Crew Sample Data 1981 - current ver 34. Environmental Data Initiative.
## https://doi.org/10.6073/pasta/c85ded1d123a76125690a3d14f773d7a. Accessed 2026-08-28.
revision <- list_data_package_revisions(
  scope = "knb-lter-ntl",
  identifier = "31",
  filter = "newest"
)
package_id <- paste("knb-lter-ntl", "31", revision, sep = ".")
res <- read_data_entity_names(package_id)
raw <- read_data_entity(package_id, entityId = res$entityId[1])
secchi <- readr::read_csv(file = raw, show_col_types = FALSE)
readr::write_csv(secchi, "data/secchi.csv")
cat("wrote data/secchi.csv (", package_id, ")\n")

## --- Color - Trout Lake Area -> data/color.csv ---
## Magnuson, J.J., S.R. Carpenter, and E.H. Stanley. 2026. North Temperate Lakes LTER: Color -
## Trout Lake Area 1989 - current ver 36. Environmental Data Initiative.
## https://doi.org/10.6073/pasta/55b3e50beec0a49dbcb147be607c1f55. Accessed 2026-08-28.
revision <- list_data_package_revisions(
  scope = "knb-lter-ntl",
  identifier = "87",
  filter = "newest"
)
package_id <- paste("knb-lter-ntl", "87", revision, sep = ".")
res <- read_data_entity_names(package_id)
raw <- read_data_entity(package_id, entityId = res$entityId[1])
color <- readr::read_csv(file = raw, show_col_types = FALSE)
readr::write_csv(color, "data/color.csv")
cat("wrote data/color.csv (", package_id, ")\n")

## --- Chlorophyll - Trout Lake Area -> data/chl.csv ---
## Magnuson, J.J., S.R. Carpenter, and E.H. Stanley. 2025. North Temperate Lakes LTER: Chlorophyll -
## Trout Lake Area 1981 - current ver 33. Environmental Data Initiative.
## https://doi.org/10.6073/pasta/659e43f4796f71e3a37673a0f2d7a77a. Accessed 2026-08-28.
revision <- list_data_package_revisions(
  scope = "knb-lter-ntl",
  identifier = "35",
  filter = "newest"
)
package_id <- paste("knb-lter-ntl", "35", revision, sep = ".")
res <- read_data_entity_names(package_id)
raw <- read_data_entity(package_id, entityId = res$entityId[1])
chl <- readr::read_csv(file = raw, show_col_types = FALSE)
readr::write_csv(chl, "data/chl.csv")
cat("wrote data/chl.csv (", package_id, ")\n")

cat("\nAll data pulled to data/. Run the fig*.R and 01_plot_profiles.R scripts next.\n")
