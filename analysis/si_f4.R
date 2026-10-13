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
library(data.table)
library(ggplot2)
library(tidyr)
library(patchwork)
library(scales)

# df 
pfd_md = fread(file.path(data_dir, "pfd_metadata.tsv")) 
pfd_md_rd = pfd_md %>%
  mutate(source = if_else(grepl("MFD", member_id), "MFD", "DB")) %>%
  select(representative, member_id, db, source, ptu, dptu, length, gc) %>%
  distinct(representative, .keep_all = T) 

# dptu groups
dptu_df = pfd_md_rd %>%
  filter(!is.na(dptu)) %>%
  group_by(dptu) %>% 
  mutate(dptu_size = n()) %>%
  mutate(
    dptu_category = case_when(
      all(source == "MFD") ~ "MFD",
      all(source == "DB") ~ "DB",
      TRUE ~ "Shared"
    )
  ) %>%
  ungroup()

# colours
col_db = "gray"
col_shared = "#FFD700"
col_mfd = "blue"

# theme
figure_theme = theme_bw(base_size = 12, base_family = "sans") +
  theme(
    strip.text = element_text(size = 12, face = "plain"),
    axis.title = element_text(size = 12),
    plot.margin = margin(0, 5.5, 0, 5.5),
    axis.text = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.66),
    legend.position = "none"
  )

# dptu
# dptu counts
dptu_sum = dptu_df %>%
  group_by(dptu, dptu_category) %>%
  summarise(
    n_db = sum(source == "DB"),
    n_mfd = sum(source == "MFD"),
    total = n(),
    .groups = "drop"
  )

# db only
db_only = dptu_sum %>%
  filter(dptu_category == "DB") %>%
  arrange(desc(total)) %>%
  transmute(
    dptu,
    gain = total,
    group = "DB-only dPTUs"
  ) %>%
  mutate(
    x = row_number(),
    y0 = lag(cumsum(gain), default = 0),
    y1 = cumsum(gain)
  )

db_last_x = max(db_only$x)
db_last_y = max(db_only$y1)

# shared
shared = dptu_sum %>%
  filter(dptu_category == "Shared") %>%
  arrange(desc(total)) %>%
  transmute(
    dptu,
    gain = total,
    group = "Shared dPTUs"
  ) %>%
  mutate(
    x = db_last_x + row_number(),
    y0 = db_last_y + lag(cumsum(gain), default = 0),
    y1 = db_last_y + cumsum(gain)
  )

shared_last_x = max(shared$x)
shared_last_y = max(shared$y1)

# mfd only
mfd_only = dptu_sum %>%
  filter(dptu_category == "MFD") %>%
  arrange(desc(total)) %>%
  transmute(
    dptu,
    gain = total,
    group = "MFD-only dPTUs"
  ) %>%
  mutate(
    x = shared_last_x + row_number(),
    y0 = shared_last_y + lag(cumsum(gain), default = 0),
    y1 = shared_last_y + cumsum(gain)
  )

# join curves
db_only_line = db_only

shared_line = bind_rows(
  tibble(
    dptu = NA_character_,
    gain = 0,
    group = "Shared dPTUs",
    x = db_last_x,
    y0 = db_last_y,
    y1 = db_last_y
  ),
  shared
)

mfd_only_line = bind_rows(
  tibble(
    dptu = NA_character_,
    gain = 0,
    group = "MFD-only dPTUs",
    x = shared_last_x,
    y0 = shared_last_y,
    y1 = shared_last_y
  ),
  mfd_only
)

plot_df = bind_rows(db_only, shared, mfd_only)
line_df = bind_rows(db_only_line, shared_line, mfd_only_line)

dptu_cumu = ggplot() +
  geom_ribbon(
    data = plot_df,
    aes(x = x, ymin = y0, ymax = y1, fill = group),
    alpha = 0.35
  ) +
  geom_line(
    data = line_df,
    aes(x = x, y = y1, color = group),
    linewidth = 2.5
  ) +
  scale_fill_manual(
    name = "dPTU categories",
    guide = "none",
    breaks = c("DB-only dPTUs", "Shared dPTUs", "MFD-only dPTUs"),
    labels = c("DB", "Shared", "MFD"),
    values = c(
      "DB-only dPTUs" = col_db,
      "Shared dPTUs" = col_shared,
      "MFD-only dPTUs" = col_mfd
    )
  ) +
  scale_color_manual(
    name = "dPTU categories",
    breaks = c("DB-only dPTUs", "Shared dPTUs", "MFD-only dPTUs"),
    labels = c("DB-only", "Shared", "MFD-only"),
    values = c(
      "DB-only dPTUs" = col_db,
      "Shared dPTUs" = col_shared,
      "MFD-only dPTUs" = col_mfd
    ),
    guide = guide_legend(override.aes = list(linewidth = 1.5))
  ) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(
    x = "dPTUs",
    y = "Accumulative plasmid count"
  ) +
  figure_theme +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.05),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", color = NA)
  )
dptu_cumu

# shared members
shared_order = shared$dptu

shared_dptu_wide = dptu_df %>%
  filter(dptu_category == "Shared") %>%
  group_by(dptu, source) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = source, values_from = n, values_fill = 0) %>%
  mutate(
    total = DB + MFD,
    frac_mfd = MFD / total,
    dptu = factor(dptu, levels = shared_order)
  ) %>%
  arrange(dptu)

shared_dptu_long = shared_dptu_wide %>%
  pivot_longer(cols = c(DB, MFD), names_to = "source", values_to = "n")

# total size
p_top = ggplot(shared_dptu_wide, aes(x = dptu, y = total)) +
  geom_col(fill = col_shared, width = 1) +
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 0.1),
    breaks = c(2, 10, 100, 500, 2000, 20000),
    labels = comma
  ) +
  labs(
    x = NULL,
    y = "Total count (Shared)"
  ) +
  figure_theme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

# mfd count
p_mid = ggplot(shared_dptu_wide, aes(x = dptu, y = MFD)) +
  geom_col(fill = col_mfd, width = 1) +
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 0.1),
    breaks = c(1, 2, 10, 50, 150),
    labels = comma
  ) +
  labs(
    x = NULL,
    y = "MFD count (Shared)"
  ) +
  figure_theme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

# composition
p_bottom = ggplot(shared_dptu_long, aes(x = dptu, y = n, fill = source)) +
  geom_col(position = "fill", width = 1) +
  scale_fill_manual(
    name = "Source in shared dPTUs",
    values = c("DB" = col_db, "MFD" = col_mfd)
  ) +
  scale_y_continuous(labels = percent) +
  labs(
    x = "Shared dPTUs",
    y = "Composition [%]"
  ) +
  figure_theme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  )

dptu_shared = p_top / p_mid / p_bottom + plot_layout(heights = c(1, 1, 1))
dptu_shared

# size groups
dptu_size_plot = dptu_df %>%
  distinct(dptu, dptu_category, dptu_size) %>%
  filter(dptu_category %in% c("DB", "Shared", "MFD")) %>%
  mutate(
    size_bin = case_when(
      dptu_size == 1 ~ "1",
      dptu_size >= 2 & dptu_size <= 9 ~ "2-9",
      dptu_size >= 10 & dptu_size <= 99 ~ "10-99",
      dptu_size >= 100 ~ "100+"
    ),
    size_bin = factor(size_bin, levels = c("1", "2-9", "10-99", "100+")),
    dptu_category = factor(dptu_category, levels = c("DB", "Shared", "MFD"))
  ) %>%
  count(dptu_category, size_bin)

dptu_bar = ggplot(dptu_size_plot, aes(x = size_bin, y = n, fill = dptu_category)) +
  geom_col(width = 0.75) +
  geom_text(
    aes(label = comma(n)),
    vjust = -0.3,
    hjust = 0.5,
    size = 8 / .pt
  ) +
  facet_wrap(~ dptu_category, nrow = 1,
    labeller = as_labeller(c(DB = "DB-only", Shared = "Shared", MFD = "MFD-only"))) +
  scale_x_discrete(expand = expansion(add = 0.6), guide = guide_axis(n.dodge = 1)) +
  scale_fill_manual(values = c(
    "DB" = col_db,
    "Shared" = col_shared,
    "MFD" = col_mfd
  )) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.25))
  ) +
  labs(
    x = "dPTU size distribution per dPTU category",
    y = "No. of dPTUs"
  ) +
  figure_theme +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )
dptu_bar

# combine plots
# Remove vertical margins from the nested layouts as well as each panel.
block_theme = theme(plot.margin = margin(0, 5.5, 0, 5.5))
left_plots = (dptu_cumu / dptu_bar) +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(theme = block_theme)
right_plots = ((p_top / p_mid / p_bottom) +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(theme = block_theme)) & theme(legend.position = "none")

# a small gap between the two blocks
final_dptu = (
  wrap_plots(
      left_plots, plot_spacer(), right_plots,
      nrow = 1,
      widths = grid::unit(c(4, 0.15, 1.5), c("null", "cm", "null"))
    )
)
final_dptu

# save figure
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = file.path(figure_dir, "si_f4.png"),
  plot = final_dptu + plot_annotation(
    theme = theme(plot.margin = margin(0, 0, 0, 0, unit = "mm"))
  ),
  width = 21,
  height = 14,
  units = "cm",
  dpi = 300,
  bg = "white"
)


# Trim blank rows at the top and bottom; preserve horizontal whitespace.
source(file.path(dirname(figure_dir), "analysis", "trim_export_height.R"))
trim_export_height(file.path(figure_dir, "si_f4.png"))
