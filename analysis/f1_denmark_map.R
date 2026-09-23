# directories
project_dir = path.expand("~/Desktop/pfd_github")
data_dir = file.path(project_dir, "data")
figure_dir = file.path(project_dir, "figures")

# packages
library(dplyr)
library(data.table)
library(readxl)
library(ggplot2)
library(sf)
library(grid)

# samples
ont = fread(file.path(data_dir, "raw_ont_prediction.tsv")) %>% 
  mutate(technology = "ONT") %>%
  select(sample, technology) %>%
  distinct(sample, .keep_all = T)

illumina = fread(file.path(data_dir, "raw_illumina_plasmid_prediction.tsv")) %>%
  select(sample, technology) %>%
  distinct(sample, .keep_all = T) 

mfd_meta = read_excel(file.path(data_dir, "2025-05-28_mfd_db.xlsx")) %>%
  rename(sample = fieldsample_barcode) %>%
  select(sample, latitude, longitude, mfd_sampletype)

raw_pfd = bind_rows(ont, illumina) %>%
  arrange(sample, desc(technology == "ONT")) %>%
  distinct(sample, .keep_all = TRUE) %>%
  left_join(mfd_meta, by = "sample") %>%
  filter(
    !is.na(longitude),
    !is.na(latitude),
    !is.na(mfd_sampletype),
    !is.na(technology)
  ) %>%
  arrange(technology == "ONT")


shp_path = file.path(data_dir, "CNTR_RG_01M_2024_4326.shp")

dk = st_read(shp_path) %>%
  st_make_valid() %>%
  st_transform(4326) %>%
  filter(CNTR_ID %in% c("DK", "DE", "SE", "NO"))

sampletype.palette = c(
  "Sediment" = "#222222",
  "Soil"     = "#9A7A55",
  "Water"    = "#8EC3DA"
)

col_grid = "grey90"
bg_col = "white"
land_col = "grey90"

# bornholm
map_bh = ggplot() +
  geom_sf(data = dk, fill = "grey90", colour = "black", linewidth = 0.6) +
  geom_point(
    data = raw_pfd,
    aes(x = longitude, y = latitude, colour = mfd_sampletype),
    size = 3.5,
    alpha = 1,
    shape = 16
  ) +
  geom_point(
    data = raw_pfd %>% filter(technology == "ONT"),
    aes(x = longitude, y = latitude),
    shape = 21,
    fill = NA,
    colour = "darkred",
    size = 3.5,
    stroke = 1.1
  ) +
  scale_colour_manual(values = sampletype.palette) +
  coord_sf(
    xlim = c(14, 16),
    ylim = c(54.5, 55.75),
    expand = FALSE,
    clip = "on"
  ) +
  scale_y_continuous(breaks = c(54.5, 55, 55.5)) +
  scale_x_continuous(breaks = c(14, 15, 16)) +
  labs(x = "", y = "") +
  annotate(
    geom = "text",
    x = 14.75,
    y = 54.6,
    label = "Bornholm",
    size = 12 / ggplot2::.pt
  ) +
  theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid = element_line(color = col_grid, linewidth = 0.8),
    axis.title = element_blank(),
    axis.text = element_text(size = 10, colour = "black"),
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0)
  )

# denmark
map_dk = ggplot() +
  geom_sf(data = dk, fill = "grey90", colour = "black", linewidth = 0.6) +
  
  geom_point(
    data = raw_pfd,
    aes(x = longitude, y = latitude, colour = mfd_sampletype),
    size = 3.5,
    alpha = 1,
    shape = 16
  ) +
  
  geom_point(
    data = raw_pfd %>% filter(technology == "ONT"),
    aes(x = longitude, y = latitude, shape = "ONT"),
    fill = NA,
    colour = "darkred",
    size = 3.5,
    stroke = 1.1
  ) +
  
  scale_colour_manual(
    values = sampletype.palette,
    name = NULL,
    guide = guide_legend(
      order = 1,
      override.aes = list(
        shape = 16,
        size = 3.5,
        fill = NA,
        stroke = 0
      )
    )
  ) +
  
  scale_shape_manual(
    values = c("ONT" = 21),
    name = NULL,
    guide = guide_legend(
      order = 2,
      override.aes = list(
        colour = "darkred",
        fill = NA,
        size = 3.5,
        stroke = 1.3
      )
    )
  ) +
  
  coord_sf(
    xlim = c(7.5, 13.1),
    ylim = c(54.5, 58.1),
    expand = FALSE,
    clip = "on"
  ) +
  labs(x = "", y = "") +
  
  annotation_custom(
    ggplotGrob(map_bh),
    xmin = 10.7,
    ymin = 56.9,
    xmax = 13.0,
    ymax = 58.0
  ) +
  
  theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    panel.background = element_rect(fill = bg_col, colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid = element_line(color = col_grid, linewidth = 0.8),
    
    axis.title = element_blank(),
    axis.text = element_text(size = 10, colour = "black"),
    
    legend.position = c(0.23, 0.90),
    legend.box = "vertical",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 10, colour = "black"),
    
    legend.box.spacing = unit(-0.25, "cm"),
    legend.spacing.y = unit(0.2, "cm"),
    legend.key.size = unit(0.35, "cm"),
    legend.margin = margin(0, 0, 0, 0, unit = "cm")
  )


map_dk

# save figure
figure_width = 21
figure_height = 24
figure_margin = margin(5, 5, 5, 5, unit = "mm")

dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = file.path(figure_dir, "f1_denmark_map.png"),
  plot = map_dk + theme(plot.margin = figure_margin),
  width = figure_width,
  height = figure_height,
  units = "cm",
  dpi = 600,
  bg = "white"
)
