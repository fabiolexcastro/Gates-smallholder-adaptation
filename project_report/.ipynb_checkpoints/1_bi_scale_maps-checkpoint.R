# bicolor plot
# remotes::install_github("slu-openGIS/biscale")
dir <- "~/work/hazard_layers/report_maps"

# load dependencies
library(biscale)
library(sf)
library(raster)
library(ggplot2)
library(cowplot)

# load layers
rr <- file.path(dir, c("hazard_count.tif", "hazard_score.tif"))
rr <- stack(rr)

# convert to sf as raster is not supported by biscale
rsp <- rasterToPolygons(rr) # takes some time
rsf <- st_as_sf(rsp)

# saving this for futur use
saveRDS(rsf, file.path(dir, "hazard_score_count_sf_grid.rds"))

# read sf-grid
rsf <- readRDS(file.path(dir, "hazard_score_count_sf_grid.rds"))
st_crs(rsf) <- 4326

# test in smaller area
# v <- getData("GADM", country = "UGA", level = 0, path = dir)
# v <- st_as_sf(v)
# rsf <- st_crop(rsf, v)

# create classes
data <- bi_class(rsf, x = hazard_count, y = hazard_score, style = "equal", dim = 3)

map <- ggplot() +
  geom_sf(data , mapping=aes(fill=bi_class), colour=NA, show.legend=FALSE) +
  bi_scale_fill(pal="GrPink", dim=3) +
  bi_theme()

legend <- bi_legend(pal = "GrPink",
                    dim = 3,
                    xlab = "hazard count ",
                    ylab = "hazard score ",
                    size = 8)

allplot <- ggdraw() +
  draw_plot(map, 0, 0, 1, 1) +
  draw_plot(legend, 0.05, 0.05, 0.3, 0.3)

save_plot(file.path(dir, "bivariate_plot_hazards.png"), allplot,
          base_asp = 1.2)
