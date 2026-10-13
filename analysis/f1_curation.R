# directories
# Run from the repository root, or knit from analysis/.
source(if (file.exists("analysis/project_setup.R")) {
  "analysis/project_setup.R"
} else {
  "project_setup.R"
})
data_dir = file.path(project_dir, "data")
figure_dir = file.path(project_dir, "figures")

# packages
library(dplyr)
library(data.table)
library(ggplot2)

# data
df_drep = fread(file.path(data_dir, "drep_scores.tsv"))

# score distribution of dereplicated representatives
figure_theme = theme_bw(base_size = 6, base_family = "sans") +
  theme(
    text = element_text(size = 6),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6, colour = "black"),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.ticks = element_line(colour = "black", linewidth = 0.3),
    plot.margin = margin(1, 1, 1, 1, unit = "mm")
  )

score_counts = df_drep %>%
  filter(!is.na(total_score_no_mob_double_counts)) %>%
  mutate(source = factor(source, levels = c("IMG/PR & PLSDB", "MFD"))) %>%
  count(total_score_no_mob_double_counts, source, name = "n_plasmids")

score_plot = ggplot(score_counts,
                    aes(x = total_score_no_mob_double_counts,
                        y = n_plasmids, fill = source)) +
  geom_col(width = 0.9, colour = "black", linewidth = 0.25) +
  scale_fill_manual(values = c("IMG/PR & PLSDB" = "gray", "MFD" = "blue"),
                    drop = FALSE) +
  scale_x_continuous(breaks = 0:6, expand = expansion(add = 0.55)) +
  scale_y_continuous(breaks = seq(0, 80000, 20000), limits = c(0, 80000),
                     labels = function(x) format(x, big.mark = ",", scientific = FALSE,
                                                 trim = TRUE),
                     expand = expansion(mult = c(0, 0.06))) +
  labs(x = "Score (points)", y = NULL, fill = "Databases") +
  figure_theme +
  theme(legend.position = c(0.025, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.key.size = grid::unit(4, "mm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 7, hjust = 0))

score_plot
print(score_plot)
ggsave(file.path(figure_dir, "f1_curation.png"), score_plot,
       width = 5.75, height = 3, units = "cm", dpi = 600, bg = "white")
ggsave(file.path(figure_dir, "f1_curation.svg"), score_plot,
       device = svglite::svglite, fix_text_size = FALSE,
       width = 5.75, height = 3, units = "cm", bg = "white")
