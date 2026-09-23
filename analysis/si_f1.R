# directories
project_dir = path.expand("~/Desktop/pfd_github")
data_dir = file.path(project_dir, "data")
figure_dir = file.path(project_dir, "figures")

# packages
library(dplyr)
library(data.table)
library(ggplot2)
library(patchwork)

# theme
figure_theme = theme_bw(base_size = 12, base_family = "sans") +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.66),
    legend.position = "none"
  )

# df
pfd_md = fread(file.path(data_dir, "pfd_md_26_09_09.tsv")) %>%
  distinct(representative, .keep_all = T)

# length
general_stats_db = pfd_md %>%
  mutate(source = ifelse(db %in% c("IMG/PR", "PLSDB"), "DB", "MFD")) %>%
  group_by(db) %>%
  summarise(
    mean_length = mean(length, na.rm = TRUE),
    median_length = median(length, na.rm = TRUE),
    mean_gc = mean(gc, na.rm = TRUE),
    median_gc = median(gc, na.rm = TRUE),
    count = n()
  )

general_stats_source = pfd_md %>%
  mutate(source = ifelse(db %in% c("IMG/PR", "PLSDB"), "DB", "MFD")) %>%
  group_by(source) %>%
  summarise(
    mean_length = mean(length, na.rm = TRUE),
    median_length = median(length, na.rm = TRUE),
    mean_gc = mean(gc, na.rm = TRUE),
    median_gc = median(gc, na.rm = TRUE),
    count = n()
  )

# technology
mfd_technology = pfd_md %>%
  filter(db == "MFD") %>%
  filter(technology %in% c("ONT", "Illumina"))

stats_technology = pfd_md %>%
  filter(db == "MFD") %>%
  group_by(technology) %>%
  summarise(
    mean_length = mean(length, na.rm = TRUE),
    median_length = median(length, na.rm = TRUE),
    count = n()
  )
  

# # wilcoxon test
# length_test = wilcox.test(
#   length ~ technology,
#   data = mfd_technology,
#   alternative = "two.sided"
# )
# 
# length_test$p.value
# 
# format.pval(length_test$p.value, digits = 3)

# boxplot gc
gc = pfd_md %>%
  mutate(source = ifelse(db %in% c("IMG/PR", "PLSDB"), "DB", "MFD")) %>%
  mutate(source = factor(source, levels = c("MFD", "DB")))

mean_gc = gc %>%
  group_by(source) %>%
  summarise(mean_gc = mean(gc, na.rm = TRUE), .groups = "drop")

plot_gc = ggplot(gc, aes(x = source, y = gc, fill = source)) +
  geom_violin(trim = FALSE, alpha = 0.35, color = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.8) +
  geom_crossbar(
    data = mean_gc,
    aes(x = source, y = mean_gc, ymin = mean_gc, ymax = mean_gc),
    inherit.aes = FALSE,
    width = 0.45,
    linewidth = 0.8,
    color = "black"
  ) +
  scale_fill_manual(values = c("MFD" = "navy", "DB" = "tomato")) +
  figure_theme +
  labs(x = NULL, y = "GC [%]", fill = NULL)
plot_gc

# length boxplot source
length_df_source = pfd_md %>%
  mutate(source = ifelse(db %in% c("IMG/PR", "PLSDB"), "DB", "MFD")) %>%
    mutate(source = factor(source, levels = c("MFD", "DB")))

mean_length_source = length_df_source %>%
  group_by(source) %>%
  summarise(mean_length = mean(length, na.rm = TRUE), .groups = "drop")

plot_length_source = ggplot(length_df_source, aes(x = source, y = length, fill = source)) +
  geom_violin(trim = FALSE, alpha = 0.35, color = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.8) +
  geom_crossbar(
    data = mean_length_source,
    aes(x = source, y = mean_length, ymin = mean_length, ymax = mean_length),
    inherit.aes = FALSE,
    width = 0.45,
    linewidth = 0.8,
    color = "black"
  ) +
  scale_y_log10(
    breaks = c(1e3, 1e4, 1e5, 1e6, 1e7),
    labels = scales::label_log()
  ) +
  scale_fill_manual(values = c("MFD" = "navy", "DB" = "tomato")) +
  figure_theme +
  labs(x = NULL, y = expression("Length (log"[10]*") [bp]"), fill = NULL)
plot_length_source

# length boxplot db
length_df_db = pfd_md %>%
  mutate(db = ifelse(is.na(db), "MFD", db)) %>%
  mutate(db = factor(db, levels = c("MFD", "IMG/PR", "PLSDB")))

mean_length_db = length_df_db %>%
  group_by(db) %>%
  summarise(mean_length = mean(length, na.rm = TRUE), .groups = "drop")

plot_length_db = ggplot(length_df_db, aes(x = db, y = length, fill = db)) +
  geom_violin(trim = FALSE, alpha = 0.35, color = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.8) +
  geom_crossbar(
    data = mean_length_db,
    aes(x = db, y = mean_length, ymin = mean_length, ymax = mean_length),
    inherit.aes = FALSE,
    width = 0.45,
    linewidth = 0.8,
    color = "black"
  ) +
  scale_y_log10(
    breaks = c(1e3, 1e4, 1e5, 1e6, 1e7),
    labels = scales::label_log()
  ) +
  scale_fill_manual(values = c(
    "MFD" = "navy",
    "IMG/PR" = "grey55",
    "PLSDB" = "grey80"
  )) +
  figure_theme +
labs(x = NULL, y = expression("Length (log"[10]*") [bp]"), fill = NULL)
plot_length_db

# save figure
supplementary_stats = plot_gc + plot_length_source + plot_length_db +
  plot_layout(nrow = 1, widths = c(1, 1, 1.2))
supplementary_stats

dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = file.path(figure_dir, "si_f1.png"),
  plot = supplementary_stats +
    plot_annotation(
      theme = theme(plot.margin = margin(0, 20, 0, 20, unit = "mm"))
    ),
  width = 21,
  height = 7.5,
  units = "cm",
  dpi = 600,
  bg = "white"
)
