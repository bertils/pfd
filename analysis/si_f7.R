# directories
project_dir = normalizePath(if (dir.exists("data") && dir.exists("analysis")) {
  "."
} else if (dir.exists("../data") && dir.exists("../analysis")) {
  ".."
} else {
  stop("Run from the repository root or its analysis/ directory.")
}, winslash = "/", mustWork = TRUE)
data_dir = file.path(project_dir, "data")
figure_dir = file.path(project_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# packages
library(dplyr)
library(ggplot2)
library(data.table)

pfd = fread(file.path(data_dir, "pfd_metadata.tsv"))

upset = pfd %>%
  select(representative, rps) %>%
  filter(!is.na(rps))

library(tidyr)
library(patchwork)
library(scales)

# Split the comma-separated annotations; count each RP once per plasmid.
rp_long = upset %>%
  separate_longer_delim(rps, delim = ", ") %>%
  distinct(representative, rps)

set_sizes = rp_long %>% count(rps, sort = TRUE)
rp_order = set_sizes$rps

# Canonical combinations ensure that RP order within a string does not matter.
combinations = rp_long %>%
  arrange(rps) %>%
  group_by(representative) %>%
  summarise(rps = paste(rps, collapse = ", "), .groups = "drop") %>%
  count(rps, sort = TRUE) %>%
  slice_head(n = 30) %>%
  mutate(Group = row_number())

members = combinations %>%
  separate_longer_delim(rps, delim = ", ") %>%
  mutate(row = match(rps, rev(rp_order)))
connections = members %>%
  group_by(Group) %>%
  summarise(bottom = min(row), top = max(row), .groups = "drop")
background = expand_grid(Group = 1:30, row = seq_along(rp_order))
stripes = data.frame(row = seq(1, length(rp_order), by = 2))
set_sizes = set_sizes %>% mutate(row = match(rps, rev(rp_order)))

# Same typography and borders as the other supplementary figures.
figure_theme = theme_bw(base_size = 12, base_family = "sans") +
  theme(axis.title = element_text(size = 12),
    plot.margin = margin(0, 5.5, 0, 5.5),
        axis.text = element_text(size = 10, colour = "black"),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.66))

intersection_plot = ggplot(combinations, aes(Group, n)) +
  geom_col(fill = "grey35", width = 0.85) +
  geom_text(data = filter(combinations, Group != 1),
            aes(label = comma(n)), angle = 90, hjust = -0.1, size = 10 / .pt) +
  geom_text(data = filter(combinations, Group == 1),
            aes(y = n - 800, label = comma(n)), angle = 90, hjust = 1,
            colour = "white", size = 10 / .pt) +
  scale_x_continuous(limits = c(0.5, 30.5), breaks = 1:30,
                     minor_breaks = NULL, expand = expansion(mult = 0)) +
  scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.04))) +
  labs(x = NULL, y = NULL) + figure_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()
        #panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.25),
        #panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.5)
        )

matrix_plot = ggplot(background, aes(Group, row)) +
  geom_rect(data = stripes, aes(ymin = row - 0.5, ymax = row + 0.5),
            xmin = -Inf, xmax = Inf, fill = "grey95", colour = NA,
            inherit.aes = FALSE) +
  geom_point(colour = "grey85", size = 3.2) +
  geom_segment(data = connections,
               aes(x = Group, xend = Group, y = bottom, yend = top),
               inherit.aes = FALSE, linewidth = 0.4) +
  geom_point(data = members, colour = "black", size = 3.2) +
  scale_x_continuous(limits = c(0.5, 30.5), breaks = NULL, expand = expansion(mult = 0)) +
  scale_y_continuous(breaks = seq_along(rp_order), labels = rev(rp_order),
                     limits = c(0.5, length(rp_order) + 0.5), expand = expansion(mult = 0)) +
  labs(x = "Group", y = NULL) + figure_theme +
  theme(axis.ticks = element_blank())

set_plot = ggplot(set_sizes, aes(n, row)) +
  geom_rect(data = stripes, aes(ymin = row - 0.5, ymax = row + 0.5),
            xmin = -Inf, xmax = Inf, fill = "grey95", colour = NA,
            inherit.aes = FALSE) +
  geom_col(orientation = "y", width = 0.65, fill = "grey35") +
  scale_x_continuous(labels = label_comma(), breaks = c(0, 25000, 50000, 75000),
                  expand = expansion(mult = c(0, 0.04))) +
  scale_y_continuous(limits = c(0.5, length(rp_order) + 0.5), expand = expansion(mult = 0)) +
  labs(x = "Set size", y = NULL) + figure_theme +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 10, hjust = 1))

# Patchwork layout: bars above the RP matrix; set sizes on the right.
layout = "
A#
BC
"

upset_plot = wrap_plots(
  A = intersection_plot,
  B = matrix_plot,
  C = set_plot,
  design = layout,
  widths = c(3, 1),       # matrix column : set-size column
  heights = c(14, 17)     # top row : bottom row
) +
  plot_annotation(theme = theme(
    plot.margin = margin(0, 0, 0, 0, unit = "mm")
  ))
print(upset_plot)

# Export size in centimetres; panel proportions are controlled above.
figure_width = 24
figure_height = 15
ggsave(file.path(figure_dir, "si_f7.png"), upset_plot,
       width = figure_width, height = figure_height,
       units = "cm", dpi = 300, bg = "white")


# Trim blank rows at the top and bottom; preserve horizontal whitespace.
source(file.path(dirname(figure_dir), "analysis", "trim_export_height.R"))
trim_export_height(file.path(figure_dir, "si_f7.png"))
