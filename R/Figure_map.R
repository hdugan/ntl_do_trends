# ==============================================================================
# map.R
# NTL-LTER Study Lakes — Southern (Yahara) & Northern (Vilas County) regions
#
# Produces a single, publication-quality graphic with two side-by-side maps:
#   LEFT  — Southern region, Yahara lakes chain (Madison, WI)
#   RIGHT — Northern region, Northwoods lakes (Vilas County, WI)
#
# Basemaps are real map tiles (so actual lake shapes render automatically —
# no shapefiles required), with clean typography, a matching accent-color
# palette per region, scale bars, and north arrows.
#
# NOTE: This script downloads basemap tiles at run time, so it needs an
# internet connection when you run it.
# ==============================================================================

# ---- 0. Packages ------------------------------------------------------------
required_pkgs <- c(
    "sf", # spatial vector data
    "maptiles", # download basemap raster tiles (CartoDB / OSM / Esri, etc.)
    "tidyterra", # ggplot2 geoms for terra SpatRaster objects (the tile output)
    "ggplot2",
    "ggrepel", # non-overlapping point labels
    "ggspatial", # scale bar + north arrow
    "patchwork", # combine the two maps side by side
    "data.table", # captions.csv merge-write
    "maps" # world outline for the globe inset
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
    install.packages(missing_pkgs)
}

library(sf)
library(maptiles)
library(tidyterra)
library(ggplot2)
library(ggrepel)
library(ggspatial)
library(patchwork)
library(data.table)

# ---- 1. Options you may want to tweak ---------------------------------------

# Fish Lake is an official NTL "southern" study lake but sits ~30 miles from
# the other three (Mendota/Monona/Wingra), which widens the southern map a
# lot. Set to FALSE to keep the southern panel tight on the Yahara chain.
include_fish_lake <- TRUE

# Basemap style. "Esri.WorldGrayCanvas" = clean/light, no API key required.
# NOTE: as of 2024 CartoDB's free tile endpoints (CartoDB.Positron etc.) require
# an API key and silently return tiles stamped "API KEY REQUIRED" without one --
# do not switch back to a CartoDB.* provider without adding maptiles::set_key().
# Other nice options to try: "Esri.WorldTopoMap" (more color/labels),
# "Esri.WorldImagery" (satellite), "OpenStreetMap" (classic OSM).
basemap_provider <- "Esri.WorldGrayCanvas"

# Accent colors, one per region.
color_south <- "#0E7C86" # teal — southern (Yahara) lakes
color_north <- "#1F4B3F" # forest green — northern (Northwoods) lakes

output_file <- "figures/figMap_regions.png"

# ---- 2. Lake coordinates -----------------------------------------------------
# Coordinates are the official NTL-LTER sampling-station locations.

southern_lakes <- data.frame(
    name = c("Lake Mendota", "Lake Monona", "Lake Wingra", "Fish Lake"),
    lat = c(43.09885, 43.06337, 43.05258, 43.28733),
    lon = c(-89.40545, -89.36086, -89.42499, -89.65173),
    nudge_x = 0, # metres, applied to the ggrepel label starting position
    stringsAsFactors = FALSE
)

if (!include_fish_lake) {
    southern_lakes <- subset(southern_lakes, name != "Fish Lake")
}

northern_lakes <- data.frame(
    name = c(
        "Allequash Lake", "Big Muskellunge Lake", "Crystal Bog",
        "Crystal Lake", "Sparkling Lake", "Trout Bog", "Trout Lake"
    ),
    lat = c(46.038317, 46.021067, 46.007583, 46.00275, 46.007733, 46.04125, 46.029267),
    lon = c(-89.620617, -89.611783, -89.606183, -89.612233, -89.701183, -89.686283, -89.665017),
    # Big Muskellunge Lake's label otherwise gets placed to the west, straight over
    # the Trout Lake point -- nudge its ggrepel starting position east to clear it.
    nudge_x = c(0, 900, 0, 0, 0, 0, 0),
    stringsAsFactors = FALSE
)

# ---- 3. Convert to sf points, project to Web Mercator (matches tile CRS) ----

to_points <- function(df) {
    sf::st_as_sf(df, coords = c("lon", "lat"), crs = 4326) |>
        sf::st_transform(3857)
}

pts_south <- to_points(southern_lakes)
pts_north <- to_points(northern_lakes)

# ---- 3b. Lake outlines (NHD polygons, data/GIS/NHLDshp) so each study lake's ----
# actual shape can be filled in light blue rather than relying on the basemap's
# own (often faint/absent) water rendering. No shapefile is bundled for Fish
# Lake -- it just keeps its plain point marker.
lake_shp_dir <- "data/GIS/NHLDshp"
lake_shapefiles <- c(
    "Lake Mendota" = "mendota", "Lake Monona" = "monona", "Lake Wingra" = "wingra",
    "Allequash Lake" = "allequash", "Big Muskellunge Lake" = "bigmusky",
    "Crystal Bog" = "CrystalBog", "Crystal Lake" = "crystal", "Sparkling Lake" = "sparkling",
    "Trout Bog" = "TroutBog", "Trout Lake" = "Trout"
)

read_lake_polys <- function(lake_names) {
    shp <- lake_shapefiles[intersect(names(lake_shapefiles), lake_names)]
    polys <- lapply(names(shp), function(nm) {
        g <- sf::st_read(file.path(lake_shp_dir, paste0(shp[[nm]], ".shp")), quiet = TRUE)
        g$study_name <- nm
        g[, "study_name"] # keep just the name + geometry columns
    })
    do.call(rbind, polys) |> sf::st_transform(3857)
}

polys_south <- read_lake_polys(southern_lakes$name)
polys_north <- read_lake_polys(northern_lakes$name)

lake_fill <- "#CFEAF5"
lake_border <- "#7FB8D4"

# ---- 4. Helper: build a padded bounding box around a set of points ---------

padded_bbox <- function(pts, pad_frac = 0.35, target_aspect = NULL) {
    bb <- sf::st_bbox(pts)
    dx <- (bb["xmax"] - bb["xmin"])
    dy <- (bb["ymax"] - bb["ymin"])
    # guard against a zero-width bbox (e.g. a single point)
    dx <- max(dx, 1000)
    dy <- max(dy, 1000)
    bb["xmin"] <- bb["xmin"] - dx * pad_frac
    bb["xmax"] <- bb["xmax"] + dx * pad_frac
    bb["ymin"] <- bb["ymin"] - dy * pad_frac
    bb["ymax"] <- bb["ymax"] + dy * pad_frac

    # The two regions have very different natural footprints (north is a wide
    # east-west scatter, south is a tall north-south chain). Left alone, that
    # forces two side-by-side panels of very different heights into one
    # patchwork row -- the row height stretches to the taller one and the
    # shorter panel centers with dead space above/below it. Pad the shorter
    # dimension further so every panel shares the same target aspect ratio
    # and fully fills its row.
    if (!is.null(target_aspect)) {
        dx <- (bb["xmax"] - bb["xmin"])
        dy <- (bb["ymax"] - bb["ymin"])
        current_aspect <- dy / dx # height / width
        if (current_aspect < target_aspect) {
            new_dy <- dx * target_aspect
            grow <- (new_dy - dy) / 2
            bb["ymin"] <- bb["ymin"] - grow
            bb["ymax"] <- bb["ymax"] + grow
        } else {
            new_dx <- dy / target_aspect
            grow <- (new_dx - dx) / 2
            bb["xmin"] <- bb["xmin"] - grow
            bb["xmax"] <- bb["xmax"] + grow
        }
    }

    sf::st_as_sfc(bb, crs = 3857)
}

# height / width each panel should render at, so the south (tall) and north
# (wide) point clusters end up filling the same row height side by side
panel_target_aspect <- 0.85

bbox_south <- padded_bbox(pts_south, pad_frac = 0.45, target_aspect = panel_target_aspect)
bbox_north <- padded_bbox(pts_north, pad_frac = 0.55, target_aspect = panel_target_aspect)

# The padding above is based on the sampling-point bbox, but Trout Lake's actual
# shoreline (data/GIS/NHLDshp) extends much farther north than its point, so the
# cluster renders crowded against the panel's north edge with a lot of empty
# basemap left over to the south. Shift the whole window north to rebalance.
bbox_north <- sf::st_set_crs(bbox_north + c(0, 1400), 3857)

# ---- 5. Download basemap tiles -----------------------------------------------

tiles_south <- maptiles::get_tiles(bbox_south, provider = basemap_provider, crop = TRUE)
tiles_north <- maptiles::get_tiles(bbox_north, provider = basemap_provider, crop = TRUE)

# ---- 6. Shared map-panel builder --------------------------------------------

build_panel <- function(tiles, pts, polys = NULL, accent, title, subtitle, caption = NULL) {
    pts_df <- pts

    p <- ggplot() +
        tidyterra::geom_spatraster_rgb(data = tiles, maxcell = 5e6) +
        {
            if (!is.null(polys)) geom_sf(data = polys, fill = lake_fill, color = lake_border, linewidth = 0.35, inherit.aes = FALSE)
        } +
        geom_sf(
            data = pts_df,
            shape = 21,
            fill = accent, color = "white", stroke = 1.1, size = 3.36 # 20% smaller than the original 4.2
        ) +
        ggrepel::geom_label_repel(
            data = pts_df,
            aes(label = name, geometry = geometry),
            stat = "sf_coordinates",
            seed = 42,
            nudge_x = pts_df$nudge_x,
            size = 2.6, # ~2pt smaller than the original 3.3mm (size is in mm; 2pt = 2 / .pt mm)
            fontface = "bold",
            color = "grey15",
            fill = alpha("white", 0.7),
            label.size = 0,
            label.r = unit(0.15, "lines"),
            segment.color = accent,
            segment.size = 0.5,
            box.padding = 0.55,
            min.segment.length = 0
        ) +
        ggspatial::annotation_scale(
            location = "br", width_hint = 0.28,
            line_width = 0.8, text_cex = 0.7,
            bar_cols = c(accent, "white")
        ) +
        ggspatial::annotation_north_arrow(
            location = "tr", which_north = "true",
            style = ggspatial::north_arrow_minimal(text_size = 8, line_width = 1.2),
            height = unit(0.9, "cm"), width = unit(0.9, "cm")
        ) +
        # Titles are long, so wrap to 2 lines rather than letting bold text
        # run past the ~2.9in panel width.
        labs(title = paste(strwrap(title, width = 32), collapse = "\n"), subtitle = subtitle, caption = caption) +
        coord_sf(expand = FALSE) +
        theme_void(base_family = "sans") +
        theme(
            plot.title = element_text(face = "bold", size = 10, color = accent, hjust = 0, lineheight = 1.05, margin = margin(b = 3)),
            plot.subtitle = element_text(size = 8.3, color = "grey30", hjust = 0, margin = margin(t = 1, b = 6)),
            plot.caption = element_text(size = 7, color = "grey45", hjust = 0, margin = margin(t = 5)),
            plot.margin = margin(6, 8, 6, 8),
            panel.border = element_rect(color = accent, fill = NA, linewidth = 1.1)
        )

    p
}

p_south <- build_panel(
    tiles_south, pts_south, polys_south,
    accent = color_south,
    title = "Southern Lakes",
    subtitle = "Dane County, Wisconsin"
)

p_north <- build_panel(
    tiles_north, pts_north, polys_north,
    accent = color_north,
    title = "Northern Lakes",
    subtitle = "Vilas County, Wisconsin"
)

# ---- 6b. Globe inset centered on Madison, WI, dropped into the bottom-left ----
# of the southern panel (the scale bar moved to bottom-right above to make room).

center_lon <- -89.4012
center_lat <- 43.0731

world <- sf::st_as_sf(maps::map("world", plot = FALSE, fill = TRUE)) |>
    sf::st_set_crs(4326) |>
    sf::st_make_valid()

center_pt <- sf::st_sfc(sf::st_point(c(center_lon, center_lat)), crs = 4326)

# visible hemisphere = everything within ~90 degrees (great-circle) of the center point
# radius in meters, just under a quarter of Earth's circumference
hemisphere <- sf::st_buffer(center_pt, dist = 9800000)

world_visible <- suppressWarnings(sf::st_intersection(world, hemisphere))

ortho_crs <- sprintf("+proj=ortho +lat_0=%f +lon_0=%f", center_lat, center_lon)

world_globe <- sf::st_transform(world_visible, crs = ortho_crs)
globe_outline <- sf::st_transform(hemisphere, crs = ortho_crs) # the circular "ocean" backdrop
point_globe <- sf::st_transform(center_pt, crs = ortho_crs)

grat <- sf::st_graticule(lon = seq(-180, 180, 20), lat = seq(-80, 80, 20)) |>
    sf::st_intersection(hemisphere) |>
    sf::st_transform(ortho_crs)

us_inset <- ggplot() +
    geom_sf(data = globe_outline, fill = "#eaf3fa", color = NA) + # ocean/sky
    geom_sf(data = grat, color = "grey85", linewidth = 0.15) + # graticule
    geom_sf(data = world_globe, fill = "grey80", color = "white", linewidth = 0.15) +
    geom_sf(data = point_globe, color = "#D7191C", size = 2.2) +
    theme_void() +
    theme(
        panel.background = element_rect(fill = color_south, color = NA),
        plot.background = element_rect(fill = color_south, color = NA),
        # White frame so the inset pops against the basemap it sits on top of.
        panel.border = element_rect(color = "white", fill = NA, linewidth = 1.1),
        plot.margin = margin(3, 3, 3, 3)
    )

# plot.background is set on each base panel here (rather than with a trailing
# `&` on the combined patchwork below) because `&` broadcasts to every nested
# plot in the tree, including the inset -- that would silently paint the
# inset's background white again despite the element_blank() above.
p_south <- p_south + theme(plot.background = element_rect(fill = "white", color = NA))
p_north <- p_north + theme(plot.background = element_rect(fill = "white", color = NA))

p_south <- p_south + patchwork::inset_element(
    us_inset,
    left = 0.00, bottom = 0.02, right = 0.30, top = 0.32,
    align_to = "panel"
)

# ---- 7. Combine side by side --------------------------------------------------

# No overall title/explainer drawn here -- same convention as the other fig*.R
# scripts (see Figure4_chl_profiles.R): that text goes to figures/captions.csv
# only. The per-panel region names stay on the PNG since they're panel labels,
# not the figure's headline/explainer.
combined <- p_south | p_north

# ---- 8. Save ------------------------------------------------------------------
# Height is sized to the actual content (panel titles + map area at
# panel_target_aspect + caption), not an arbitrary tall canvas -- otherwise
# ggsave pads the extra height as dead white space top/bottom of the figure.
panel_width_in <- 6.5 / 2 - (8 + 6) / 72 # half the page minus left/right plot.margin
map_height_in <- panel_width_in * panel_target_aspect
# Titles are now a single short line each ("Southern Lakes" / "Northern Lakes"
# + one-line subtitle), so the title block needs much less headroom than the
# 2-line estimate this budget once assumed -- was measured leaving ~0.28in of
# true dead space above AND below the whole figure.
output_height <- 0.42 + map_height_in + 0.05 # panel title (1 line)/subtitle + map + bottom margin

ggsave(output_file, combined, width = 6.5, height = output_height, dpi = 500, bg = "white")

message("Saved: ", normalizePath(output_file))

## Title + caption are NOT drawn on the PNG -- written to the shared
## figures/captions.csv instead (see Figure1_rate_of_change.R for the same convention).
write_captions <- function(new_caps) {
    path <- "figures/captions.csv"
    old <- if (file.exists(path)) data.table::fread(path) else data.table::data.table(file = character(), title = character(), caption = character())
    data.table::fwrite(rbind(old[!file %in% new_caps$file], new_caps), path)
    cat("wrote", path, "\n")
}
write_captions(data.table::data.table(
    file = output_file,
    title = "NTL-LTER study lake regions",
    caption = paste0(
        "Locations of the 11 NTL-LTER primary study lakes across two contrasting Wisconsin landscapes: ",
        "Dane County (southern, 4 lakes) and Vilas County (northern, 7 lakes). Each lake's outline is filled ",
        "light blue from NHD polygons (data/GIS/NHLDshp); no shapefile is bundled for Fish Lake, which keeps ",
        "a plain point marker. Southern panel includes a globe inset locating Madison, WI globally. ",
        "Basemap: Esri World Light Gray Canvas."
    )
))
