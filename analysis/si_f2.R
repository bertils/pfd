# directories
project_dir = path.expand("~/Desktop/pfd_github")
data_dir = file.path(project_dir, "data")
figure_dir = file.path(project_dir, "figures")


# packages
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)
library(patchwork)
library(scales)

# df 
pfd_md = fread(file.path(data_dir, "pfd_md_26_09_09.tsv")) 
pfd_md_rd = pfd_md %>%
  mutate(source = if_else(grepl("MFD", member_id), "MFD", "DB")) %>%
  select(representative, member_id, db, source, ptu, dptu, length, gc) %>%
  distinct(representative, .keep_all = T) 

# ptu groups
ptu_df = pfd_md_rd %>%
  filter(!is.na(ptu)) %>%
  group_by(ptu) %>% 
  mutate(ptu_size = n()) %>%
  mutate(
    ptu_category = case_when(
      all(source == "MFD") ~ "MFD",
      all(source == "DB") ~ "DB",
      TRUE ~ "Shared"
    )
  ) %>%
  ungroup()

# colours
col_db = "#ef6351"
col_shared = "grey50"
col_mfd = "#1f4e8c"

# theme
figure_theme = theme_bw(base_size = 12, base_family = "sans") +
  theme(
    strip.text = element_text(size = 12, face = "plain"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.66),
    legend.position = "none"
  )

# ptu
# ptu counts
ptu_sum = ptu_df %>%
  group_by(ptu, ptu_category) %>%
  summarise(
    n_db = sum(source == "DB"),
    n_mfd = sum(source == "MFD"),
    total = n(),
    .groups = "drop"
  )

# db only
db_only = ptu_sum %>%
  filter(ptu_category == "DB") %>%
  arrange(desc(total)) %>%
  transmute(
    ptu,
    gain = total,
    group = "DB-only PTUs"
  ) %>%
  mutate(
    x = row_number(),
    y0 = lag(cumsum(gain), default = 0),
    y1 = cumsum(gain)
  )

db_last_x = max(db_only$x)
db_last_y = max(db_only$y1)

# shared
shared = ptu_sum %>%
  filter(ptu_category == "Shared") %>%
  arrange(desc(total)) %>%
  transmute(
    ptu,
    gain = total,
    group = "Shared PTUs"
  ) %>%
  mutate(
    x = db_last_x + row_number(),
    y0 = db_last_y + lag(cumsum(gain), default = 0),
    y1 = db_last_y + cumsum(gain)
  )

shared_last_x = max(shared$x)
shared_last_y = max(shared$y1)

# mfd only
mfd_only = ptu_sum %>%
  filter(ptu_category == "MFD") %>%
  arrange(desc(total)) %>%
  transmute(
    ptu,
    gain = total,
    group = "MFD-only PTUs"
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
    ptu = NA_character_,
    gain = 0,
    group = "Shared PTUs",
    x = db_last_x,
    y0 = db_last_y,
    y1 = db_last_y
  ),
  shared
)

mfd_only_line = bind_rows(
  tibble(
    ptu = NA_character_,
    gain = 0,
    group = "MFD-only PTUs",
    x = shared_last_x,
    y0 = shared_last_y,
    y1 = shared_last_y
  ),
  mfd_only
)

plot_df = bind_rows(db_only, shared, mfd_only)
line_df = bind_rows(db_only_line, shared_line, mfd_only_line)

ptu_cumu = ggplot() +
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
    name = "PTU category",
    breaks = c("DB-only PTUs", "Shared PTUs", "MFD-only PTUs"),
    labels = c("DB", "Shared", "MFD"),
    values = c(
      "DB-only PTUs" = col_db,
      "Shared PTUs" = col_shared,
      "MFD-only PTUs" = col_mfd
    )
  ) +
  scale_color_manual(
    values = c(
      "DB-only PTUs" = col_db,
      "Shared PTUs" = col_shared,
      "MFD-only PTUs" = col_mfd
    ),
    guide = "none"
  ) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(
    x = "Number of PTUs added",
    y = "Cumulative plasmid count"
  ) +
  figure_theme +
  theme(legend.position = "right")
ptu_cumu

# shared members
shared_order = shared$ptu

shared_ptu_wide = ptu_df %>%
  filter(ptu_category == "Shared") %>%
  group_by(ptu, source) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = source, values_from = n, values_fill = 0) %>%
  mutate(
    total = DB + MFD,
    frac_mfd = MFD / total,
    ptu = factor(ptu, levels = shared_order)
  ) %>%
  arrange(ptu)

shared_ptu_long = shared_ptu_wide %>%
  pivot_longer(cols = c(DB, MFD), names_to = "source", values_to = "n")

# total size
p_top = ggplot(shared_ptu_wide, aes(x = ptu, y = total)) +
  geom_col(fill = col_shared, width = 1) +
  scale_y_sqrt(breaks = c(0, 500, 1500, 2500), labels = comma) +
  labs(
    x = NULL,
    y = "Plasmids\nper PTU"
  ) +
  figure_theme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

# mfd count
p_mid = ggplot(shared_ptu_wide, aes(x = ptu, y = MFD)) +
  geom_col(fill = col_mfd, width = 1) +
  scale_y_sqrt(breaks = c(0, 1, 3, 5), labels = comma) +
  labs(
    x = NULL,
    y = "MFD plasmids\nper PTU"
  ) +
  figure_theme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

# composition
p_bottom = ggplot(shared_ptu_long, aes(x = ptu, y = n, fill = source)) +
  geom_col(position = "fill", width = 1) +
  scale_fill_manual(
    name = "Source in shared PTUs",
    values = c("DB" = col_db, "MFD" = col_mfd)
  ) +
  scale_y_continuous(labels = percent) +
  labs(
    x = "Shared PTUs\n(largest to smallest)",
    y = "Plasmid source (%)"
  ) +
  figure_theme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  )

ptu_shared = p_top / p_mid / p_bottom + plot_layout(heights = c(1, 1, 3))
ptu_shared

# size groups
ptu_size_plot = ptu_df %>%
  distinct(ptu, ptu_category, ptu_size) %>%
  filter(ptu_category %in% c("DB", "Shared", "MFD")) %>%
  mutate(
    size_bin = case_when(
      ptu_size == 1 ~ "1",
      ptu_size >= 2 & ptu_size <= 9 ~ "2-9",
      ptu_size >= 10 & ptu_size <= 99 ~ "10-99",
      ptu_size >= 100 ~ "100+"
    ),
    size_bin = factor(size_bin, levels = c("1", "2-9", "10-99", "100+")),
    ptu_category = factor(ptu_category, levels = c("DB", "Shared", "MFD"))
  ) %>%
  count(ptu_category, size_bin)

ptu_bar = ggplot(ptu_size_plot, aes(x = size_bin, y = n, fill = ptu_category)) +
  geom_col(width = 0.75) +
  geom_text(
    aes(label = comma(n), vjust = if_else(size_bin == "100+", -1.7, -0.3)),
    hjust = 0.5,
    size = 10 / .pt
  ) +
  facet_wrap(~ ptu_category, nrow = 1) +
  scale_x_discrete(expand = expansion(add = 1.1)) +
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
    x = "Plasmids per PTU",
    y = "Number of PTUs"
  ) +
  figure_theme +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ptu_bar

# combine plots
final_ptu = (
  ptu_cumu + p_top + p_mid + ptu_bar + p_bottom +
    plot_layout(
      design = "
      AB
      AC
      DE
      ",
      widths = c(2.3, 1.2),
      heights = c(1, 1, 1)
    )
) & theme(legend.position = "none")
final_ptu

# save figure
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = file.path(figure_dir, "si_f2.png"),
  plot = final_ptu + plot_annotation(
    theme = theme(plot.margin = margin(0, 10, 0, 10, unit = "mm"))
  ),
  width = 21,
  height = 14,
  units = "cm",
  dpi = 600,
  bg = "white"
)
