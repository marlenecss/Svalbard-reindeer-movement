## Scientific Graphs
# Assignment 08.08.2026
# by Marlene Sehrbrock

library(move2)
library(sf)
library(terra)
library(ggplot2)
library(rnaturalearth)
library(mapview)
library(lubridate)
library(units)
library(magrittr)
library(dplyr)
library(legendry)
library(tidyterra)
library(maptiles)
library(patchwork)
library(showtext)
library(sysfonts)

#setwd("C:/WORKING DIRECTORY")

# font
font_add_google("Source Sans 3", "Source Sans 3")
showtext_auto()
showtext_opts(dpi = 96)  # screen preview; switch to 600 before ggsave

######## 1 get data ############################################################
#movebank_store_credentials("USERNAME", getPass::getPass())
keyring::key_list()

# get study metadata
studies <- movebank_download_study_info()
study_info <- studies[studies$id == 2608802883, ]

# study overview
study_info$sensor_type_ids
study_info$number_of_individuals
study_info$taxon_ids
study_info$principal_investigator_name
study_info$i_can_see_data
study_info$i_have_download_access
study_info$timestamp_first_deployed_location
study_info$timestamp_last_deployed_location

#mapview(study_info$main_location)

study_id <- 2608802883

# deployment metadata
study_indis <- movebank_download_deployment(
  study_id,
  'license-md5' = 'cb787f617ff25cc66359ddf0f2495428'
)

# download all individuals
study_data <- movebank_download_study(study_id)
nrow(study_data)
head(study_data)
#mapview(study_data)

# download one individual for inspection
study_data_indi <- movebank_download_study(
  study_id = study_id,
  individual_local_identifier = "B123"
)
#mapview(study_data_indi)

##### 2 explore move2 object ###################################################
st_coordinates(study_data_indi)
st_crs(study_data_indi)
class(study_data_indi)

study_bbox <- st_bbox(study_data_indi)
mapview(study_bbox) + mapview(study_info$main_location)

# time lags
mt_time(study_data)

# sort by track and time stamp
study_data <- study_data[order(mt_track_id(study_data), mt_time(study_data)), ]

# now filter for even hours
study_data <- study_data %>%
mutate(hour = hour(mt_time(.))) %>%
  filter(hour %% 2 == 0)

study_data_lags <- mt_time_lags(study_data, units = "min")
range(study_data_lags,  na.rm = TRUE)
median(study_data_lags, na.rm = TRUE)
sd(study_data_lags,     na.rm = TRUE)

# track info
mt_track_data(study_data)
mt_track_id(study_data)
mt_n_tracks(study_data)
as.data.frame(table(mt_track_id(study_data)))

# splitting / stacking
study_data_split <- split(study_data, mt_track_id(study_data))
mt_stack(study_data_split)

# coercion
study_data_move <- to_move(study_data)
study_data_sf   <- study_data
class(study_data_sf) <- class(study_data) %>% setdiff("move2")

##### 3 choose individuals and filter to even hours ############################
# inspect individual coverage
as.data.frame(table(mt_track_id(study_data))) %>% arrange(desc(Freq))

# pick 3 with long, continuous records
chosen <- c("B244", "B245", "B251")

data_sub <- study_data %>%
  filter(mt_track_id(.) %in% chosen) %>%
  mutate(
    individual = as.character(mt_track_id(.)),
    timestamp  = mt_time(.),
    date       = as_date(timestamp),
    hour       = hour(timestamp),
    month_num  = month(timestamp),
    season     = case_when(
      month_num %in% c(6, 7, 8)   ~ "Arctic Summer",
      month_num %in% c(12, 1, 2)  ~ "Arctic Winter",
      TRUE                         ~ NA_character_
    )
  )

##### 4 calculate daily total distance #########################################
# reproject to UTM zone 33 for distances in metres
utm_crs  <- st_crs("+proj=utm +zone=33 +datum=WGS84 +units=m")
data_utm <- st_transform(data_sub, utm_crs)

# step lengths between consecutive fixes per individual
data_utm <- data_utm %>%
  group_by(individual) %>%
  arrange(timestamp, .by_group = TRUE) %>%
  mutate(
    step_m = as.numeric(
      st_distance(geometry, lag(geometry), by_element = TRUE)
    )
  ) %>%
  ungroup()

# sum to daily total distance per individual
daily_dist <- data_utm %>%
  st_drop_geometry() %>%
  filter(!is.na(step_m)) %>%
  group_by(individual, date) %>%
  summarise(
    daily_km  = sum(step_m, na.rm = TRUE) / 1000,
    month_num = first(month(date)),
    season    = first(season),
    .groups   = "drop"
  )

# save
st_write(study_data, "svalbard_study_data.gpkg", delete_dsn = TRUE)
write.csv(daily_dist, "svalbard_daily_dist.csv", row.names = FALSE)

# load (next session — skip everything above step 3)
library(move2)
library(sf)

study_data <- st_read("svalbard_study_data.gpkg")

# st_read drops the move2 class, so restore it
study_data <- mt_as_move2(
  study_data,
  time_column      = "timestamp",
  track_id_column  = "individual_local_identifier"
)

daily_dist <- read.csv("svalbard_daily_dist.csv") %>%
  mutate(date = as.Date(date))

##### 5 satellite basemap ######################################################
bbox_4326 <- st_bbox(st_transform(data_sub, 4326))

tiles <- maptiles::get_tiles(
  st_as_sfc(bbox_4326) %>% st_set_crs(4326),
  provider = "Esri.WorldImagery",
  zoom     = 10     # increase to 11-12 for more detail (slower)
)

##### 6 build track lines ######################################################
tracks <- data_sub %>%
  group_by(individual) %>%
  arrange(timestamp, .by_group = TRUE) %>%
  summarise(
    geometry = st_cast(st_combine(geometry), "LINESTRING"),
    .groups  = "drop"
  )

##### 7 colours ################################################################
ind_colors <- setNames(
  c("#440154", "#31688E", "#35B779"),
  chosen
)

# convert date to numeric for continuous colour scale
data_sub_steps <- data_sub_steps %>%
  mutate(date_num = as.numeric(date))

# updated shared range
date_range <- range(data_sub_steps$date_num, na.rm = TRUE)

# nice date labels for the legend
date_breaks <- pretty(as.Date(c(min(data_sub_steps$date), 
                                max(data_sub_steps$date))), n = 4)

##### 8 map plots — one per individual, shared extent & viridis scale #########

txt <- 10

make_map <- function(ind) {
  points_i <- data_sub_steps %>% filter(individual == ind)
  
  ggplot() +
    tidyterra::geom_spatraster_rgb(data = tiles) +
    geom_sf(
      data = points_i,
      aes(color = date_num),
      size = 0.8, alpha = 0.7
    ) +
    scale_color_viridis_c(
      option  = "D",
      name    = "Date",
      limits  = date_range,
      breaks  = as.numeric(date_breaks),
      labels  = format(date_breaks, "%b %Y"),
      guide   = "none"
    ) +
    coord_sf(
      xlim   = c(bbox_shared["xmin"], bbox_shared["xmax"]),
      ylim   = c(bbox_shared["ymin"], bbox_shared["ymax"]),
      expand = FALSE
    ) +
    labs(title = ind) +
    theme_void(base_size = txt) +
    theme(
      text            = element_text(family = "Source Sans 3"),
      plot.title      = element_text(color = "white", size = txt + 1,
                                     face = "bold", hjust = 0.5,
                                     family = "Source Sans 3"),
      plot.background = element_rect(fill = "#1a1a1a", color = NA)
    )
}

map1 <- make_map(chosen[1])
map2 <- make_map(chosen[2])
map3 <- make_map(chosen[3])

maps_row <- map1 + map2 + map3 +
  plot_annotation(
    title    = "Svalbard Reindeer Movement",
    subtitle = "GPS tracks — even-hour fixes, coloured by step length",
    theme    = theme(
      plot.background = element_rect(fill = "#1a1a1a", color = NA),
      plot.title      = element_text(color = "white", size = txt + 3,
                                     face = "bold", family = "Source Sans 3"),
      plot.subtitle   = element_text(color = "white", size = txt + 1,
                                     family = "Source Sans 3")
    )
  )

##### 9 histogram with legendry — same viridis scale as map ###################
# join date_num onto daily_dist
daily_dist <- daily_dist %>%
  mutate(date_num = as.numeric(date))

hist_plot <- ggplot(
  daily_dist,
  aes(x = daily_km, color = individual)
) +
  geom_histogram(
    aes(fill = after_stat(x)),
    bins     = 30,
    alpha    = 0.9,
    position = "identity"
  ) +
  scale_fill_viridis_c(
    option = "D",
    name   = "Daily km",
    guide  = compose_stack(
      gizmo_histogram(
        hist.args = list(breaks = 20, col = NA),
        just      = 0.5,
        alpha     = 0.8
      ),
      "colorbar",
      theme = theme(
        legend.text       = element_text(color = "white", size = txt,
                                         family = "Source Sans 3"),
        legend.title      = element_text(color = "white", size = txt,
                                         family = "Source Sans 3"),
        legend.frame      = element_rect(color = NA),
        legend.key.height = unit(1.5, "cm"),
        legend.key.width  = unit(0.3, "cm")
      )
    )
  ) +
  scale_color_manual(
    values = ind_colors,
    name   = "Individual",
    guide  = guide_legend(
      theme = theme(
        legend.text  = element_text(color = "white", size = txt,
                                    family = "Source Sans 3"),
        legend.title = element_text(color = "white", size = txt,
                                    family = "Source Sans 3")
      )
    )
  ) +
  scale_x_continuous(name = "Daily distance travelled (km)") +
  facet_wrap(~ individual, ncol = 3,
             labeller = as_labeller(setNames(chosen, chosen))) +
  labs(y = "Days") +
  theme_minimal(base_size = txt) +
  theme(
    text             = element_text(family = "Source Sans 3"),
    plot.background  = element_rect(fill = "#1a1a1a", color = NA),
    panel.background = element_rect(fill = "#1a1a1a", color = NA),
    panel.grid       = element_line(color = "#333333"),
    strip.text       = element_text(color = "white", face = "bold",
                                    size = txt, family = "Source Sans 3"),
    axis.text        = element_text(color = "white", size = txt - 1,
                                    family = "Source Sans 3"),
    axis.title       = element_text(color = "white", size = txt,
                                    family = "Source Sans 3"),
    legend.position  = "right"
  )

##### 10 combine and save ######################################################

final_plot <- maps_row / hist_plot +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(
    caption = "Data: Svalbard Reindeer Project – Nordenskiöld | Graphic: Marlene",
    theme   = theme(
      plot.background = element_rect(fill = "#1a1a1a", color = NA),
      plot.caption    = element_text(color = "grey60", size = txt - 2,
                                     hjust = 0, family = "Source Sans 3")
    )
  )

final_plot

showtext_opts(dpi = 600)
ggsave("svalbard_reindeer.png", final_plot,
       width = 20, height = 18, units = "cm", dpi = 600)
showtext_opts(dpi = 96)
