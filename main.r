#######################################################
#                 Zonal statistics of classes with R
#                 Milos Popovic
#                 2023/06/13
########################################################

setwd()
main <- getwd()

install.packages("remotes")
remotes::install_github(
    "ropensci/MODIStsp"
)

remotes::install_github(
    "dickoa/rgeoboundaries"
)

remotes::install_github(
    "isciences/exactextractr"
)

# libraries we need

libs <- c(
    "MODIStsp", "rgeoboundaries",
    "tidyverse", "sf", "terra",
    "exactextractr", "httr"
)

# install missing libraries

installed_libs <- libs %in% rownames(
    installed.packages()
)

if(any(installed_libs == F)){
    install.packages(libs[!installed_libs])
}

# load libraries

invisible(lapply(libs, library, character.only = T))

### 1. SUB-STATE POLYGON
### --------------------

ukraine_poly <- rgeoboundaries::gb_adm1(
    "UKR"
)

plot(
    sf::st_geometry(
        ukraine_poly
    )
)

file_path <- paste0(
    main, "/", "ukraine_lvl1.shp"
)

sf::st_write(
    ukraine_poly,
    file_path
)

### 2. MODIS LAND COVER DATA
### ------------------------
MODIStsp::MODIStsp_get_prodlayers(
    "MCD12Q1"
)

MODIStsp::MODIStsp(
    gui = F,
    out_folder = main,
    out_folder_mod = main,
    selprod = "LandCover_type_Yearly_500m (MCD12Q1)",
    bandsel = "LC1",
    user = "", #please set your user name
    password = "", # please set your password
    start_date = "2020.01.01",
    end_date = "2020.12.31",
    spatmeth = "file",
    spafile = file_path,
    out_format = "GTiff"
)

list.files(path = paste0(
    main, "/", "ukraine_lvl1/LandCover_Type_Yearly_500m_v6/LC1"
))

raster_path <- paste0(
    main, "/",
    "ukraine_lvl1/LandCover_Type_Yearly_500m_v6/LC1",
    "/", "MCD12Q1_LC1_2020_001.tif"
)

ukraine_lc <- terra::rast(
    raster_path
)

plot(ukraine_lc)

### 3. REPROJECT LAND COVER DATA
### ----------------------------

crs_longlat <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
ukraine_lc <- terra::rast(
    raster_path
) |>
    terra::project(
        crs_longlat
    )

### 4. BASIC ZONAL STATISTICS
### -------------------------

ukraine_poly$id <- 1:max(nrow(ukraine_poly))

ukraine_lc_mode <- exactextractr::exact_extract(
    ukraine_lc, ukraine_poly, "mode",
    append_cols = "id"
)|>
dplyr::inner_join(ukraine_poly, by = "id") |>
sf::st_as_sf()

map1 <- ggplot() +
    geom_sf(data = ukraine_lc_mode,
        aes(fill = factor(mode)),
        color = "white",
        size = .125
        ) +
    guides(
        fill = guide_legend(
            direction = "vertical",
            keywidth = unit(5, "mm"),
            keyheight = unit(5, "mm"),
            title.position = "top",
            label.position = "right",
            title.hjust = .5,
            label.hjust = 0,
            ncol = 1,
            byrow = F
        )
    ) +
    labs(
        title = "Land cover - Ukraine",
        caption = "Data: MODIS Land Cover Type Annual L3 Global 500m",
        x = "",
        y = ""
    ) +
    theme_void() +
    theme(
        legend.position = "left",
        plot.title = element_text(
            size = 22, color = "grey10",
            hjust = .5, vjust = 1 
        ),
        plot.caption = element_text(
            size = 10, color = "grey60",
            hjust = .5, vjust = -3 
        ),
        legend.title = element_text(
            size = 12, color = "grey10"
        ),
        legend.text = element_text(
            size = 11, color = "grey10"
        ),
        plot.margin = unit(c(t = 0, b = 0, r = 0, l = 5),
        "lines"
        )
    )

print(map1)

### 5. ZONAL STATISTICS - SINGLE CLASS
### ----------------------------------

ukraine_cropland <- exactextractr::exact_extract(
    ukraine_lc, ukraine_poly, function(value, fraction) {
        100 * sum(fraction[value == 12]) / sum(fraction)
    }
)

# convert to data.frame and transpose
ukraine_cropland_df <- as.data.frame(
    t(ukraine_cropland)
)

# wide to long format
ukraine_cropland_df_long <- ukraine_cropland_df |>
    tidyr::pivot_longer(
        cols = 1:27,
        names_to = "id",
        values_to = "value"
    )

# merge data and polygons

ukraine_cropland_sf <- cbind(
    ukraine_poly,
    ukraine_cropland_df_long
)

# region centroids

ukraine_region_cents <- ukraine_cropland_sf |>
    sf::st_centroid() |>
    dplyr::mutate(
        long = unlist(map(geometry, 1)),
        lat = unlist(map(geometry, 2))
    ) |>
    dplyr::select(
        shapeName, long, lat, value) |>
    sf::st_drop_geometry()


map2 <- ggplot() +
    geom_sf(data = ukraine_cropland_sf,
        aes(fill = value),
        color = "white",
        size = .125
        ) +
    geom_text(
        data = ukraine_region_cents,
        aes(x = long, y = lat,
        label = paste0(
            shapeName, ":", "", 
            round(value, 0), "%"
        ))
    ) +
    scale_fill_gradientn(
        name = "% of fractional cover area",
        colours = rev(terrain.colors(6)),
        limits = c(0, 100)
    ) +
    guides(
        fill = guide_legend(
            direction = "horizontal",
            keywidth = unit(10, "mm"),
            keyheight = unit(2.5, "mm"),
            title.position = "top",
            label.position = "bottom",
            title.hjust = .5,
            label.hjust = .5,
            nrow = 1,
            byrow = T
        )
    ) +
    labs(
        title = "Croplands - Ukraine",
        caption = "Data: MODIS Land Cover Type Annual L3 Global 500m",
        x = "",
        y = ""
    ) +
    theme_void() +
    theme(
        legend.position = "top",
        plot.title = element_text(
            size = 22, color = "grey10",
            hjust = .5, vjust = 1 
        ),
        plot.caption = element_text(
            size = 10, color = "grey60",
            hjust = .5, vjust = -3 
        ),
        legend.title = element_text(
            size = 12, color = "grey10"
        ),
        legend.text = element_text(
            size = 11, color = "grey10"
        ),
        plot.margin = unit(c(t = 1, b = 0, r = 0, l = 0),
        "lines"
        )
    )

ggsave(
    filename = "ukraine_croplands.png",
    width = 8, height = 7, dpi = 600,
    device = "png", bg = "white", map2
)


print(map2)

### 6. ZONAL STATISTICS - MULTICLASS
### --------------------------------

lc_values <- terra::minmax(
    ukraine_lc
)

classes <- lc_values[1]:lc_values[2]

zonal_stats_ukr <- exactextractr::exact_extract(
    ukraine_lc, ukraine_poly
)

ukraine_multiclass <- lapply(
    zonal_stats_ukr, function(x){
        as.data.frame(
            prop.table(
                table(factor(x[,1], 
                levels = classes)
            ))
        )
    }
)

ukraine_multiclass

ukraine_multiclass_df <- dplyr::bind_rows(
    ukraine_multiclass, .id = "id"
)

ukraine_multiclass_df

class(ukraine_multiclass_df$lc_class_id)

names(ukraine_multiclass_df) <- c(
    "id", "lc_class_id", "value"
)

# scrape Land cover codebook

res <- httr::GET(
    "https://ladsweb.modaps.eosdis.nasa.gov/filespec/MODIS/6/MCD12Q1"
)

txt <- httr::content(
    res, "text"
)

# turn text into data.frame
df <- read.csv(text = txt)

fix(df)

d <- df[c(449:475),] |>
    as.data.frame()

names(d) <- "cols"

d

classes_df <- d |>
    tidyr::separate(
        cols, c("keep", "discard"),
        sep = "\\s(?=unit8|uint8)"
    ) |>
    dplyr::select("keep")

unique_classes <- unique(
    classes_df$keep
)

empty_lines <- grepl(
    "^\\s*$", unique_classes
)

lc_classes_list <- unique_classes[!empty_lines]

lc_classes_df <- data.frame(
    lc_class_id = as.factor(
        
    ),
)

ukraine_multiclass_df_labs <- ukraine_multiclass_df |>
    dplyr::left_join(
        lc_classes_df,
        by = "lc_class_id"
    )

fix(ukraine_multiclass_df_labs)

ukraine_multiclass_sf <- ukraine_poly |>
    dplyr::mutate(
        id = as.character(id)
    ) |>
    dplyr::left_join(
        ukraine_multiclass_df_labs,
        by = "id"
    )

### 6. FINAL MAP
### ------------

map3 <- ggplot() +
    geom_sf(
        data = subset(
            ukraine_multiclass_sf, lc_class_id != 17),
        aes(
            fill = value * 100
        ),
        color = "white",
        size = .125
    ) +
    facet_wrap(~class) +
    # geom_text(
    #     data = ukraine_region_cents,
    #     aes(x = long, y = lat, label = paste0(
    #         shapeName, ":", "", round(value, 0), "%")
    #         ),
    #     color = "grey10",
    #     size = 3
    # ) +
    scale_fill_gradientn(
        name = "% of fractional cover area",
        colours = rev(terrain.colors(5)),
        limits = c(0, 100)
    ) +
    guides(
        fill = guide_legend(
            direction = "horizontal",
            keyheight = unit(2.5, "mm"),
            keywidth = unit(10, "mm"),
            title.position = "top",
            label.position = "bottom",
            title.hjust = .5,
            label.hjust = .5,
            nrow = 1,
            byrow = T
        )
    ) +
    labs(
        title = "Land cover class - Ukraine",
        caption = "Data source: MODIS Land Cover Type Yearly L3 Global 500m",
        x = "",
        y = ""
    ) +
    theme_void() +
    theme(
        legend.position = "top",
        plot.title = element_text(
            size = 22, color = "grey10",
            hjust = .5, vjust = 1
        ),
        plot.caption = element_text(
            size = 10, color = "grey60",
            hjust = .5, vjust = -3
        ),
        legend.title = element_text(
            size = 12, color = "grey10",
            hjust = .5
        ),
        legend.text = element_text(
            size = 11, color = "grey10",
            hjust = .5
        ),
        plot.margin = unit(
            c(t = 0, b = 0, r = 0, l = 0),
            "lines"
        )
    )

print(map3)
