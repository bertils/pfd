data_dir = if (dir.exists("data")) "data" else "../data"
figure_dir = file.path(dirname(data_dir), "figures")

# packages
library(dplyr)
library(data.table)
library(ggplot2)

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
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    legend.position = "none"
  )

gtdb_mp = fread(file.path(data_dir, "meth-vs-halo.csv"), sep = ",")

plot_df = gtdb_mp %>%
  mutate(
    c = factor(c, levels = c("Halobacteria", "Methanosarcinia")),
    source_db = if_else(db %in% c("IMG/PR", "PLSDB"), "IMG/PR & PLSDB", db),
    source_db = factor(source_db, levels = c("GTDB", "IMG/PR & PLSDB", "MFD"))
  ) %>%
  # draw MFD points last so they remain visible
  arrange(source_db == "MFD")

plot_gc = ggplot(plot_df, aes(x = c, y = gc)) +
  geom_violin(
    trim = FALSE,
    fill = "grey90",
    colour = "grey50",
    linewidth = 0.35,
    show.legend = FALSE
  ) +
  geom_point(
    aes(colour = source_db),
    position = position_jitter(width = 0.055, height = 0),
    shape = 16,
    size = 1.8,
    alpha = 1
  ) +
  facet_wrap(~ source, nrow = 1) +
  scale_colour_manual(
    name = NULL,
    values = c("GTDB" = "black", "IMG/PR & PLSDB" = "gray", "MFD" = "blue")
  ) +
  labs(x = NULL, y = "GC [%]") +
  figure_theme +
  theme(
    legend.position = "right",
    legend.box.spacing = unit(3, "mm"),
    legend.title.position = "top",
    legend.title = element_text(hjust = 0)
  )

# page size and margins
figure_width = 21
figure_height = 12
figure_margin = margin(0, 10, 0, 10, unit = "mm")
plot_gc = plot_gc + theme(plot.margin = figure_margin)

# save figure
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(
  filename = file.path(figure_dir, "si_f8.png"),
  plot = plot_gc,
  width = figure_width,
  height = figure_height,
  units = "cm",
  dpi = 300,
  bg = "white"
)


# Trim blank rows at the top and bottom; preserve horizontal whitespace.
source(file.path(dirname(figure_dir), "analysis", "trim_export_height.R"))
trim_export_height(file.path(figure_dir, "si_f8.png"))
