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
library(tidyr)
library(ggplot2)
library(ggpubr)

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

# df
c5 = read.delim(file.path(data_dir, "skani_dist_c5.tsv")) %>%
  mutate(parameter_c = "c_5")
c15 = read.delim(file.path(data_dir, "skani_dist_c15.tsv")) %>%
  mutate(parameter_c = "c_15")
c30 = read.delim(file.path(data_dir, "skani_dist_c30.tsv")) %>%
  mutate(parameter_c = "c_30")

parameter_c = bind_rows(c5, c15, c30) %>%
  filter(Ref_name == Query_name) %>%
  transmute(
    plasmid = trimws(Ref_name),
    alignment_fraction = Align_fraction_ref,
    parameter_c = factor(parameter_c, levels = c("c_5", "c_15", "c_30"))
  )

# mean af
parameter_c_stats = parameter_c %>%
  group_by(parameter_c) %>%
  summarise(n = n(), mean_af = mean(alignment_fraction), .groups = "drop")
print(parameter_c_stats)

# match plasmids
stopifnot(!anyDuplicated(parameter_c[c("plasmid", "parameter_c")]))
paired_af = parameter_c %>%
  pivot_wider(names_from = parameter_c, values_from = alignment_fraction)
stopifnot(all(complete.cases(paired_af)))

# paired t-tests
c5_c30 = t.test(paired_af$c_5, paired_af$c_30, paired = TRUE)
c5_c15 = t.test(paired_af$c_5, paired_af$c_15, paired = TRUE)
c15_c30 = t.test(paired_af$c_15, paired_af$c_30, paired = TRUE)

# comparisons
comparison_stats = data.frame(
  group1 = c("c_5", "c_5", "c_15"),
  group2 = c("c_30", "c_15", "c_30"),
  mean_difference = c(c5_c30$estimate, c5_c15$estimate, c15_c30$estimate),
  ci_lower = c(c5_c30$conf.int[1], c5_c15$conf.int[1], c15_c30$conf.int[1]),
  ci_upper = c(c5_c30$conf.int[2], c5_c15$conf.int[2], c15_c30$conf.int[2]),
  p_value = c(c5_c30$p.value, c5_c15$p.value, c15_c30$p.value),
  y.position = c(101, 103, 105)
) %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "bonferroni"),
    p_label = formatC(p_adjusted, format = "e", digits = 1)
  )
print(comparison_stats)

# plot
set.seed(1)

plot_parameter_c = ggplot(
  parameter_c,
  aes(x = parameter_c, y = alignment_fraction)
) +
  geom_boxplot(
    width = 0.18,
    linewidth = 0.66,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.18,
    height = 0,
    size = 1.4,
    alpha = 0.45
  ) +
  stat_pvalue_manual(
    comparison_stats,
    label = "p_label",
    tip.length = 0.015,
    size = 10 / ggplot2::.pt
  ) +
  scale_x_discrete(
    labels = c("c_5" = "-c 5", "c_15" = "-c 15", "c_30" = "-c 30")
  ) +
  scale_y_continuous(
    breaks = c(70, 80, 90, 100),
    expand = expansion(mult = c(0.02, 0.14))
  ) +
  labs(
    x = "skani -c parameter",
    y = "Alignment fraction [%]"
  ) +
  figure_theme

plot_parameter_c

# save figure
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = file.path(figure_dir, "si_f11.png"),
  plot = plot_parameter_c +
    theme(plot.margin = margin(0, 55, 0, 55, unit = "mm")),
  width = 21,
  height = 9,
  units = "cm",
  dpi = 300,
  bg = "white"
)


# Trim blank rows at the top and bottom; preserve horizontal whitespace.
source(file.path(dirname(figure_dir), "analysis", "trim_export_height.R"))
trim_export_height(file.path(figure_dir, "si_f11.png"))
