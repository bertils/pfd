data_dir = "data"
figure_dir = "figures"

# packages
library(dplyr)
library(data.table)
library(ggplot2)
library(svglite)
library(ggVennDiagram)

# theme
figure_theme = theme_bw(base_size = 5, base_family = "sans") +
  theme(
    axis.title = element_text(size = 5),
    axis.text = element_text(size = 5, color = "black"),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5),
    plot.title = element_text(size = 5),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.66),
    legend.position = "none"
  )

# data
prediction_columns = c("seq_name", "sample", "geNomad", "PLASMe", "Plasmer", "viralVerify")
ont = fread(file.path(data_dir, "raw_ont_prediction.tsv"), select = prediction_columns)
illumina = fread(file.path(data_dir, "raw_illumina_plasmid_prediction.tsv"), select = prediction_columns)

# count each contig and sample once
venn_df = bind_rows(ont, illumina) %>%
  mutate(member_id = paste(seq_name, sample, sep = "__")) %>%
  distinct(member_id, .keep_all = TRUE)

# prediction sets
venn_sets = list(
  geNomad = venn_df %>% filter(geNomad == 1) %>% pull(member_id),
  Plasmer = venn_df %>% filter(Plasmer == 1) %>% pull(member_id),
  viralVerify = venn_df %>% filter(viralVerify == 1) %>% pull(member_id),
  PLASMe = venn_df %>% filter(PLASMe == 1) %>% pull(member_id)
)

# venn diagram
# convert points to the millimetres used for the labels
venn_plot = ggVennDiagram(
  venn_sets,
  label = "both",
  label_alpha = 0,
  label_size = 5 / .pt,
  label_percent_digit = 2,
  set_size = 5 / .pt,
  label_font = "sans",
  edge_size = 0.35
) +
  scale_fill_gradient(low = "#F7FBFF", high = "#73BFE2") +
  figure_theme +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  )

# use the same heading positions as si_f3
venn_plot$layers[[3]]$data = venn_plot$layers[[3]]$data %>%
  mutate(Y = Y + if_else(name %in% c("Plasmer", "viralVerify"), 0.015, 0) - 0.014)

# center each count and percentage at the same position
region_labels = venn_plot$layers[[4]]$data %>%
  mutate(
    count_label = format(count, big.mark = ",", trim = TRUE),
    percent_label = paste0("(", trimws(percent), ")")
  )

# replace the combined labels with two centered lines
venn_plot$layers[[4]] = NULL
venn_plot = venn_plot +
  geom_text(
    data = region_labels,
    aes(x = X, y = Y, label = count_label),
    inherit.aes = FALSE,
    hjust = 0.5, vjust = -0.15, size = 5 / .pt, family = "sans"
  ) +
  geom_text(
    data = region_labels,
    aes(x = X, y = Y, label = percent_label),
    inherit.aes = FALSE,
    hjust = 0.5, vjust = 1.15, size = 5 / .pt, family = "sans"
  )

# leave room for the outer labels without clipping them
venn_plot = venn_plot +
  scale_x_continuous(expand = expansion(mult = 0.07)) +
  scale_y_continuous(expand = expansion(mult = 0.025)) +
  coord_cartesian(clip = "off")

# page size and margins
figure_width = 6.95
figure_height = 5.15
figure_margin = margin(1, 1, 1, 1, unit = "mm")
venn_plot = venn_plot + theme(plot.margin = figure_margin, aspect.ratio = NULL)

# save figure
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(
  filename = file.path(figure_dir, "f1_venn_diagram.png"),
  plot = venn_plot,
  width = figure_width,
  height = figure_height,
  units = "cm",
  dpi = 600,
  bg = "white"
)

# save svg
ggsave(file.path(figure_dir, "f1_venn_diagram.svg"), venn_plot,
       device = svglite, fix_text_size = FALSE,
       width = figure_width, height = figure_height,
       units = "cm", bg = "white")
