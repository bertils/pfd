data_dir = if (dir.exists("data")) "data" else "../data"
figure_dir = file.path(dirname(data_dir), "figures")

# packages
library(dplyr)
library(data.table)
library(ggplot2)
library(ape)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(grid)
library(stringr)
library(tidyverse)

# theme
figure_theme = theme_bw(base_size = 12, base_family = "sans") +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.66),
    legend.position = "none"
  )

# data
pfd_md = fread(file.path(data_dir, "pfd_metadata.tsv")) %>%
  separate_rows(member_id, sep = ",\\s*")
tree = read.tree(file.path(data_dir, "NZ_CP106969.1_91_cluster.faa.id80_c50.msa.contree"))

orc_clusters = fread(file.path(data_dir, "NZ_CP106969.1_91_cluster.faa.id80_c50_cluster.tsv"), header = F) %>%
  rename(orf_representative = V1,
         orf_member_id = V2) %>%
  mutate(member_id = sub("_[^_]+$", "", orf_member_id)) %>%
  mutate(representative = sub("_[^_]+$", "", orf_representative))


mobileog_filter = orc_clusters %>%
  filter(representative == "mobileOG") %>%
  filter(member_id != representative)

clusters_combined = orc_clusters %>% 
  filter(representative != "mobileOG") %>%
  bind_rows(mobileog_filter)

orc_plasmids = pfd_md %>%
  inner_join(clusters_combined %>% select(orf_member_id, member_id), by = "member_id") %>%
  select(orf_member_id, member_id, gc, taxa_gtdb, source) %>%
  rename(label = orf_member_id) %>%
  mutate(source = ifelse(grepl("MFD", member_id), "MFD", "DB"))


orc_plasmids = orc_plasmids %>%
  mutate(phylum_class = str_extract(taxa_gtdb, "p_[^;]+;c_[^;]+"))

# tree
replacement_vec = setNames(mobileog_filter$orf_member_id, mobileog_filter$orf_representative)

tree2 = tree
tree2$tip.label = ifelse(
  tree2$tip.label %in% names(replacement_vec),
  replacement_vec[tree2$tip.label],
  tree2$tip.label
)
tree2$edge.length = sqrt(tree2$edge.length)

c_p_tree = ggtree(tree2, layout = "circular")

c_p_tree$data = c_p_tree$data %>%
  left_join(orc_plasmids, by = "label") %>%
  mutate(
    bootstrap = suppressWarnings(as.numeric(label)),
    bootstrap_cat = case_when(
      !isTip & !is.na(bootstrap) & bootstrap > 90 ~ ">90",
      !isTip & !is.na(bootstrap) & bootstrap > 70 ~ ">70",
      !isTip & !is.na(bootstrap) ~ "<=70",
      TRUE ~ NA_character_
    )
  )

tree_data = c_p_tree$data

internal_nodes = tree_data %>%
  filter(!isTip & !is.na(bootstrap_cat))

phylum_class_cols = c(
  "p_Halobacteriota;c_Halobacteria"      = "#C58C2B",
  "p_Halobacteriota;c_Methanosarcinia"   = "#6FA8DC",
  "p_Halobacteriota;c_Bog-38"            = "#4C9A8A",
  "p_Methanobacteriota B;c_Thermococci"  = "#B07AA1",
  "p_Methanobacteriota;c_Methanobacteria" = "#D66B4D"
)

c_p4 = c_p_tree +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(x = 0, fill = source, color = source),
    width = 1, offset = 0
  ) +
  scale_fill_manual(
    name = "Source",
    values = c(MFD = "blue", DB = "gray"),
    na.value = "white",
    guide = guide_legend(order = 1)
  ) +
  scale_color_manual(
    name = "Source",
    values = c(MFD = "blue", DB = "gray"),
    na.value = "white",
    guide = guide_legend(order = 1)
  ) +
  new_scale_fill() +
  new_scale_color() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, x = 0, fill = gc, color = gc),
    width = 1, offset = 0.04
  ) +
  scale_fill_gradient(
    name = "GC (%)", low = "aliceblue", high = "navy",
    guide = guide_colorbar(order = 2)
  ) +
  scale_color_gradient(
    name = "GC (%)", low = "aliceblue", high = "navy",
    guide = guide_colorbar(order = 2)
  ) +
  new_scale_fill() +
  new_scale_color() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, x = 0, fill = phylum_class, color = phylum_class),
    width = 1, offset = 0.04
  ) +
  scale_fill_manual(
    name = "Taxonomic assignment (class-level)",
    values = phylum_class_cols,
    breaks = names(phylum_class_cols),
    labels = sub("^p_", "p__", sub(";c_", ";c__", names(phylum_class_cols))),
    na.value = "white",
    guide = guide_legend(order = 3)
  ) +
  scale_color_manual(
    name = "Taxonomic assignment (class-level)",
    values = phylum_class_cols,
    breaks = names(phylum_class_cols),
    labels = sub("^p_", "p__", sub(";c_", ";c__", names(phylum_class_cols))),
    na.value = "white",
    guide = guide_legend(order = 3)
  ) +
  new_scale_color() +
  geom_point2(
    data = internal_nodes,
    mapping = aes(color = bootstrap_cat),
    size = 2
  ) +
  scale_color_manual(
    name = "Bootstrap support",
    values = c(">90" = "black", ">70" = "grey55", "<=70" = "darkred"),
    breaks = c(">90", ">70", "<=70"),
    labels = c(">90", "70 ≤ x ≤ 90", "<70"),
    guide = guide_legend(order = 4, override.aes = list(size = 4))
  ) +
  figure_theme +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    legend.position = "right"
  )

# legend
c_p4 = c_p4 + theme(
  legend.key.size = unit(0.35, "cm"),
  legend.spacing.y = unit(0.75, "cm"),
  legend.box.spacing = unit(0.6, "cm"),
  legend.margin = margin(0, 0, 0, 4, unit = "mm"),
  plot.margin = margin(0, 0, 0, 0)
)

# branch scale
# one unit here refers to the square-root transformed branch lengths
scale_length = 1

# get the range used to draw the tree
built_tree = ggplot_build(c_p4)
tree_layout = built_tree$layout
tree_panel = tree_layout$panel_params[[1]]
radial_range = tree_panel$r.range
radial_span = diff(radial_range)

# enlarging the panel below makes the radius half the panel width
scale_width = 0.5 * scale_length / radial_span

# use ggtree to draw a horizontal scale and its label
scale_data = data.frame(x = 0, y = 0)
scale_plot = ggplot(scale_data, aes(x = x, y = y)) +
  geom_treescale(
    x = 0, y = 0.3, width = scale_length,
    offset = 0.5, fontsize = 10 / .pt,
    linesize = 0.5
  ) +
  # add a short vertical line at each end of the scale
  annotate(
    "segment",
    x = c(0, scale_length), xend = c(0, scale_length),
    y = 0.15, yend = 0.45, linewidth = 0.5
  ) +
  scale_x_continuous(limits = c(0, scale_length), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))
scale_bar = ggplotGrob(scale_plot)

# page size and margins
figure_width = 21
figure_margin = margin(0, 0, 0, 0, unit = "mm")
c_p4 = c_p4 + theme(plot.margin = figure_margin)
c_p4
# get the tree part of the figure
figure_plot = ggplotGrob(c_p4)
part_names = figure_plot$layout$name
panel_number = which(part_names == "panel")
tree_drawing = figure_plot$grobs[[panel_number]]

# enlarge only the tree to remove the empty space around the circle
# ggplot uses a radius of 0.4, so 1.25 makes it 0.5
tree_space = viewport(width = 1.25, height = 1.25)
tree_drawing$vp = tree_space

# place the scale at the bottom-left of the tree
bar_space = viewport(
  x = 0.005,
  y = 0.005,
  width = scale_width,
  height = 0.06,
  just = c("left", "bottom")
)
scale_bar$vp = bar_space

# put the tree and scale back into the figure
tree_with_scale = grobTree(tree_drawing, scale_bar)
figure_plot$grobs[[panel_number]] = tree_with_scale
figure_plot$layout$clip[panel_number] = "off"

# measure the space used by the legend, gaps, and margins
fixed_width = sum(figure_plot$widths)
fixed_width_cm = convertWidth(fixed_width, "cm", valueOnly = TRUE)
fixed_height = sum(figure_plot$heights)
fixed_height_cm = convertHeight(fixed_height, "cm", valueOnly = TRUE)

# use the remaining width for a square tree
# this adjusts automatically when the legend text size changes
tree_width = figure_width - fixed_width_cm
figure_height = tree_width + fixed_height_cm

# save figure
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(
  filename = file.path(figure_dir, "f4.png"),
  plot = figure_plot,
  width = figure_width,
  height = figure_height,
  units = "cm",
  dpi = 600,
  bg = "white"
)
