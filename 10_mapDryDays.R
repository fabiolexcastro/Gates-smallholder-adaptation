
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, velox, sf, tidyverse, gtools)
rm(list = ls())

lbl <- data.frame(value = factor(1:12), mnth = month.abb) %>% 
  mutate(mnth = factor(month.abb, levels = month.abb))
fls <- list.files('../raster/indicators/dryDays/meanMonthly', full.names = TRUE, pattern = '.tif$') %>% 
  mixedsort()
shp <- shapefile('../data/shp/base/continents.shp')

unique(shp@data$CONTINENT)
afr <- shp[shp@data$CONTINENT %in% 'Africa',]


stk <- stack(fls)
tbl <- rasterToPoints(stk, spatial = TRUE) %>% 
  as_tibble()
tbl <- tbl %>% 
  mutate(id = 1:nrow(.)) %>% 
  gather(var, value, -x, -y)
tbl <- tbl %>% 
  mutate(month = parse_number(var) %>% as.factor()) %>% 
  inner_join(., lbl, by = c('month' = 'value'))

tb2 <- tbl %>% 
  dplyr::select(x, y, mnth, value)

gg <- ggplot(tb2)  +
  geom_tile(aes(x = x, y =  y, fill = value)) +
  facet_wrap(~ mnth) +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 8, name = "GnBu")), 
                       na.value = 'white', limits = c(1, 31)) +
  theme_bw() +
  # scale_x_continuous(breaks = c(-78, -77, -76, -75)) +
  coord_equal() +
  # coord_equal(xlim = extent(vll)[1:2], ylim = extent(vll)[3:4]) +
  labs(title = 'Dry days - Africa', fill = 'days',  x = 'Longitude', y = 'Latitude') +
  theme(legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.width = unit(5, 'line')) +
  guides(shape = guide_legend(override.aes = list(size = 10)))

ggsave(plot = gg, filename = '../dryDays_Africa.png', width = 11, height = 9, units = 'in', dpi = 300)






