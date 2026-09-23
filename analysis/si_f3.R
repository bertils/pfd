# directories
project_dir = path.expand("~/Desktop/pfd_github")
data_dir = file.path(project_dir, "data")
figure_dir = file.path(project_dir, "figures")

# packages
library(dplyr)
library(data.table)
library(ggplot2)
library(ggVennDiagram)

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
  distinct(representative, .keep_all = TRUE) %>%
  select(member_id, representative, db, ptu) %>%
  filter(!is.na(ptu))

gspr = fread(file.path(data_dir, "ptus_all_gspr_mfd_db.tsv")) %>%
  rename(representative = id) %>%
  left_join(pfd_md %>% select(representative, db), by = "representative") %>%
  mutate(db = if_else(is.na(db), "GSPR", db)) %>%
  group_by(ptu) %>%
  filter(any(db == "GSPR")) %>%
  ungroup()

# venn diagram
venn_list = gspr %>%
  distinct(ptu, db) %>%
  group_by(db) %>%
  summarise(ptus = list(unique(ptu)), .groups = "drop") %>%
  tibble::deframe()

venn_list = venn_list[c("MFD", "GSPR", "IMG/PR", "PLSDB")]

# Venn labels use mm; convert the theme sizes from points.
plot_ptus = ggVennDiagram(
  venn_list,
  label = "count",
  set_size = 12 / ggplot2::.pt,
  label_size = 10 / ggplot2::.pt,
  label_font = "sans"
) +
  labs(fill = "PTU count") +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  figure_theme +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = "right"
  )
# lift the top labels
plot_ptus$layers[[3]]$data = plot_ptus$layers[[3]]$data %>%
  mutate(Y = Y + if_else(name %in% c("GSPR", "IMG/PR"), 0.015, 0))

plot_ptus

# save figure
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = file.path(figure_dir, "si_f3.png"),
  plot = plot_ptus +
    theme(plot.margin = margin(0, 30, 0, 30, unit = "mm")),
  width = 21,
  height = 9,
  units = "cm",
  dpi = 600,
  bg = "white"
)
