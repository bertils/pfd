# directories
project_dir = normalizePath(if (dir.exists("data") && dir.exists("analysis")) {
  "."
} else if (dir.exists("../data") && dir.exists("../analysis")) {
  ".."
} else {
  stop("Cannot locate the data/ and analysis/ directories.")
}, winslash = "/", mustWork = TRUE)
data_dir = file.path(project_dir, "data")
figure_dir = file.path(project_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# packages
library(dplyr)
library(data.table)
library(purrr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(scales)

# theme
figure_theme = theme_bw(base_size = 12, base_family = "sans") +
  theme(
    axis.title = element_text(size = 12),
    plot.margin = margin(0, 5.5, 0, 5.5),
    axis.text = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.66),
    legend.position = "none"
  )

# data
resolutions = c(0.001, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 0.995, 0.999, 0.9999, 1)

# read each resolution-file and combine the files
read_resolution = function(resolution) {
  filename = paste0("ptus_", resolution, ".tsv")
  fread(file.path(data_dir, filename)) %>%
    rename(representative = id) %>%
    mutate(resolution = resolution)
}
ptus = map(resolutions, read_resolution) %>%
  list_rbind()

# metadata df
pfd = fread(file.path(data_dir, "pfd_metadata.tsv")) %>%
  distinct(representative, .keep_all = TRUE) %>%
  select(representative, source)

# total ptus per resolution
ptu_counts = ptus %>%
  group_by(resolution) %>%
  summarise(ptus_per_resolution = n_distinct(ptu), .groups = "drop")

# members per ptu
ptu_sizes = ptus %>%
  group_by(resolution, ptu) %>%
  summarise(members = n_distinct(representative), .groups = "drop")

# multi-member ptus and singletons per resolution
member_counts = ptu_sizes %>%
  group_by(resolution) %>%
  summarise(
    multi_members = sum(members >= 2),
    singletons = sum(members == 1),
    .groups = "drop"
  )

# assign one category to each ptu at each resolution
ptu_categories = ptus %>%
  left_join(pfd, by = "representative") %>%
  group_by(resolution, ptu) %>%
  summarise(
    ptu_category = case_when(
      all(source == "DB") ~ "DB-only",
      all(source == "MFD") ~ "MFD-only",
      any(source == "DB") & any(source == "MFD") ~ "Shared"
    ),
    .groups = "drop"
  )

# count the distinct ptus in each category
category_counts = ptu_categories %>%
  count(resolution, ptu_category, name = "ptus_per_category") %>%
  mutate(ptu_category = factor(ptu_category, levels = c("DB-only", "MFD-only", "Shared")))

# mfd-only ptus relative to db-only and shared ptus
source_percentages = category_counts %>%
  group_by(resolution) %>%
  summarise(
    mfd_fraction = sum(ptus_per_category[ptu_category == "MFD-only"]) /
      sum(ptus_per_category[ptu_category %in% c("DB-only", "Shared")]),
    .groups = "drop"
  ) %>%
  left_join(ptu_counts, by = "resolution")

# summary table
resolution_stats = ptu_counts %>%
  left_join(member_counts, by = "resolution")
print(resolution_stats)

# total ptus
plot_ptus = ggplot(ptu_counts, aes(resolution, ptus_per_resolution)) +
  geom_point(size = 2) +
  scale_y_continuous(
    labels = label_comma(), breaks = seq(0, 160000, 40000),
    limits = c(-8000, 168000), expand = expansion(mult = 0)
  ) +
  labs(x = "Resolution (Leiden -r)", y = "No. of PTUs") +
  figure_theme

# clusters and singletons
plot_clusters = ggplot(member_counts, aes(singletons, multi_members)) +
  geom_point(size = 2) +
  geom_label_repel(
    aes(label = resolution),
    size = 10 / .pt, seed = 13,
    nudge_x = rep(c(-15000, 20000), length.out = nrow(member_counts)),
    nudge_y = rep(c(-1400, 1400), length.out = nrow(member_counts)),
    box.padding = 0.6, point.padding = 0.8,
    point.size = 2, min.segment.length = 0,
    xlim = c(12000, 173000), ylim = c(1000, 19000),
    force = 10, force_pull = 0.1, max.time = 10, max.iter = 100000,
    max.overlaps = Inf
  ) +
  scale_x_continuous(
    labels = label_comma(), breaks = c(0, 50000, 100000, 150000),
    limits = c(0, 180000), expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    labels = label_comma(), breaks = seq(0, 20000, 5000),
    limits = c(-1000, 21000), expand = expansion(mult = 0)
  ) +
  labs(x = "No. of singletons", y = "No. of multi-member PTUs") +
  figure_theme

# ptu sources
# percentages show mfd-only / (db-only + shared)
plot_sources = ggplot(category_counts, aes(factor(resolution, levels = resolutions), ptus_per_category, fill = ptu_category)) +
  geom_col(width = 0.85) +
  geom_text(
    data = source_percentages,
    aes(x = factor(resolution, levels = resolutions), y = ptus_per_resolution,
        label = percent(mfd_fraction, accuracy = 0.1)),
    inherit.aes = FALSE, vjust = -0.4, size = 10 / .pt
  ) +
  scale_fill_manual(values = c("DB-only" = "gray", "MFD-only" = "blue", Shared = "#FFD700")) +
  scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Resolution (Leiden -r)", y = "No. of PTUs", fill = "PTU categories") +
  figure_theme +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "inside",
    legend.position.inside = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = NA)
  )

# combine plots
figure_plot = (plot_ptus + plot_clusters) / plot_sources +
  plot_layout(heights = c(1.2, 1))

# page size and margins
figure_width = 21
figure_height = 18
figure_margin = margin(0, 3, 0, 3, unit = "mm")
figure_plot = figure_plot +
  plot_annotation(tag_levels = "A", theme = theme(plot.margin = figure_margin)) &
  theme(plot.tag = element_text(size = 20, face = "bold"))

# save figure
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(
  filename = file.path(figure_dir, "si_f12.png"),
  plot = figure_plot,
  width = figure_width,
  height = figure_height,
  units = "cm",
  dpi = 300,
  bg = "white"
)


# Trim blank rows at the top and bottom; preserve horizontal whitespace.
source(file.path(dirname(figure_dir), "analysis", "trim_export_height.R"))
trim_export_height(file.path(figure_dir, "si_f12.png"))
