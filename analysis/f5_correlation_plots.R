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
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(grid)
library(svglite)
library(systemfonts)
library(ragg)

# Use Inter when installed; otherwise use the platform's sans-serif font.
figure_font = if ("Inter" %in% systemfonts::system_fonts()$family) "Inter" else "sans"


# theme
figure_theme = theme_bw(base_size = 6.5, base_family = figure_font) +
  theme(
    axis.title = element_text(size = 6.5),
    axis.text = element_text(size = 6.5, color = "black"),
    legend.title = element_text(size = 6.5),
    legend.text = element_text(size = 6.5),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.66)
  )

# selected mag and plasmids, in plotting order
mag_id = "MFD01716.bin.1.47"
plasmid_ids = c("contig_90854__MFD01716", "contig_33344__MFD04752",
                "contig_57306__MFD09567")
plasmid_colors = c("#009E73", "#0072B2", "#D55E00")
min_samples = 5

# load dfs
sylph_mp = fread(file.path(data_dir, "sylph.tsv"))

# code for supplementary table
pc = fread(file.path(data_dir, "260515_new_eggnog_df.tsv"),
            select = c("member_id", "representative", "translation_table"))
pc_archaea = pc %>%
   distinct(member_id, .keep_all = TRUE) %>%
   filter(translation_table == "TT34") %>%
   distinct(representative)
df_mp2 = sylph_mp %>%
   filter((origin == "plasmid" & representative %in% pc_archaea$representative) |
            (origin == "mag" & !is.na(bin))) %>%
   mutate(representative = if_else(origin == "mag", bin, representative))

dat = df_mp2 %>%
   filter(is.finite(rel_abund)) %>%
   group_by(sample, representative, origin) %>%
   summarise(value = sum(rel_abund), .groups = "drop") %>%
   filter(value > 0)
tax = df_mp2 %>% distinct(representative, origin, taxa_gtdb)
 stopifnot(!anyDuplicated(select(tax, representative, origin)))

wide = dat %>%
   select(sample, representative, value) %>%
   pivot_wider(names_from = representative, values_from = value)
P = as.matrix(select(wide, all_of(unique(dat$representative[dat$origin == "plasmid"]))))
M = as.matrix(select(wide, all_of(unique(dat$representative[dat$origin == "mag"]))))

n = crossprod(1 * !is.na(P), 1 * !is.na(M))
r = cor(P, M, use = "pairwise.complete.obs", method = "pearson")
res = expand_grid(mag = colnames(M), plasmid = colnames(P)) %>%
   mutate(pearson = as.vector(r), n_samples = as.vector(n)) %>%
   filter(n_samples >= min_samples, is.finite(pearson)) %>%
   mutate(
     pearson = pmax(-1, pmin(1, pearson)),
     t = pearson * sqrt((n_samples - 2) / (1 - pearson^2)),
     p_value = 2 * pt(-abs(t), df = n_samples - 2),
     p_adj = p.adjust(p_value, method = "BH")
   ) %>%
   select(-t) %>%
   left_join(tax %>% filter(origin == "plasmid") %>%
               select(plasmid = representative, plasmid_taxa_gtdb = taxa_gtdb), by = "plasmid") %>%
   left_join(tax %>% filter(origin == "mag") %>%
               select(mag = representative, mag_taxa_gtdb = taxa_gtdb), by = "mag") %>%
   arrange(p_adj, desc(abs(pearson)))

 res2 = res %>% filter(pearson >= 0.5)
#fwrite(res2, file.path(data_dir, "co_corr_all_arc_supplementary.tsv"), sep = "\t")

abundance = sylph_mp %>%
  filter((origin == "plasmid" & representative %in% plasmid_ids) |
           (origin == "mag" & bin == mag_id)) %>%
  mutate(representative = if_else(origin == "mag", bin, representative)) %>%
  filter(is.finite(rel_abund)) %>%
  group_by(sample, representative) %>%
  summarise(abundance = sum(rel_abund), .groups = "drop")

# match the plasmid and mag abundances by sample
mag_data = abundance %>%
  filter(representative == mag_id) %>%
  select(sample, MAG = abundance)

plot_data = lapply(plasmid_ids, function(id) {
  plasmid_data = abundance %>%
    filter(representative == id) %>%
    select(sample, Plasmid = abundance)

  full_join(plasmid_data, mag_data, by = "sample") %>%
    mutate(
      representative = id,
      both = !is.na(Plasmid) & !is.na(MAG) & Plasmid > 0 & MAG > 0,
      Plasmid = replace_na(Plasmid, 0),
      MAG = replace_na(MAG, 0)
    ) %>%
    filter(Plasmid > 0 | MAG > 0)
}) %>% bind_rows()

# test only samples where both are detected, as in the original analysis, zeros added for plotting are not included in the test or regression
comparison_stats = plot_data %>%
  group_by(representative) %>%
  group_modify(function(data, key) {
    paired = filter(data, both)
    stopifnot(nrow(paired) >= min_samples,
              sd(paired$Plasmid) > 0, sd(paired$MAG) > 0)
    test = cor.test(paired$Plasmid, paired$MAG, method = "pearson",
                    alternative = "two.sided")
    data.frame(
      r = unname(test$estimate), p_value = test$p.value,
      n_both = nrow(paired),
      n_plasmid_only = sum(data$Plasmid > 0 & data$MAG == 0),
      n_mag_only = sum(data$MAG > 0 & data$Plasmid == 0),
      ratio_average = mean(paired$Plasmid / paired$MAG),
      ratio_sd = sd(paired$Plasmid / paired$MAG)
    )
  }) %>%
  ungroup() %>%
  arrange(match(representative, plasmid_ids))
print(comparison_stats)

x_max = max(0.8, max(plot_data$Plasmid) * 1.02)
y_max = max(0.4, max(plot_data$MAG) * 1.02)
plots = lapply(seq_along(plasmid_ids), function(i) {
  data = filter(plot_data, representative == plasmid_ids[i])
  ggplot(data, aes(Plasmid, MAG)) +
    geom_smooth(data = filter(data, both), method = "lm", formula = y ~ x,
                colour = plasmid_colors[i], fill = plasmid_colors[i],
                linewidth = 0.6, alpha = 0.15) +
    geom_point(data = filter(data, !both), colour = plasmid_colors[i],
               size = 1.1, alpha = 0.5) +
    geom_point(data = filter(data, both), colour = plasmid_colors[i],
               size = 1.6, alpha = 0.85) +
    coord_cartesian(xlim = c(-0.015, x_max), ylim = c(-0.008, y_max), expand = FALSE) +
    labs(title = bquote(.(plasmid_ids[i]) ~ bold(.(paste0("(P", i, ")")))),
         x = "Plasmid relative abundance",
         y = NULL) +
    figure_theme +
    theme(plot.title = element_text(size = 6.5, hjust = 0.5),
          axis.text.x = element_text(hjust = 1))
})

# p is the unadjusted two-sided pearson p-value
summary_data = comparison_stats %>%
  mutate(
    column = row_number(),
    r = sprintf("%.3f", r),
    p = format.pval(p_value, digits = 3, eps = 1e-99),
    both = as.character(n_both),
    plasmid_only = as.character(n_plasmid_only),
    mag_only = as.character(n_mag_only),
    average = sprintf("%.2f", ratio_average),
    sd = sprintf("± %.2f", ratio_sd)
  ) %>%
  select(column, r, p, both, plasmid_only, mag_only, average, sd) %>%
  pivot_longer(-column, names_to = "statistic", values_to = "value") %>%
  mutate(row = match(statistic, c("sd", "average", "mag_only", "plasmid_only", "both", "p", "r")),
         row = if_else(statistic == "sd", 1.25, as.numeric(row)),
         x = c(0.45, 0.74, 1.03)[column])

row_labels = data.frame(
  row = c(1.25, 2:7),
  label = c("'('*mu*' ± '*sigma*')'", "'Copies'", "italic(n)[m]*' ='",
            "italic(n)[p]*' ='", "italic(n)[p+m]*' ='", "italic(P)*' ='", "italic(r)*' ='")
)
headers = data.frame(x = c(0.45, 0.74, 1.03), label = paste0("P", 1:3))

p_table = ggplot(summary_data, aes(x, row)) +
  geom_text(aes(label = value), family = figure_font, size = 6.5 / .pt, hjust = 1) +
  geom_text(data = row_labels, aes(x = 0.025, y = row, label = label),
            family = figure_font, parse = TRUE, hjust = 0, size = 6.5 / .pt) +
  geom_point(data = headers, aes(x = x - 0.12, y = 8),
             colour = plasmid_colors, size = 2.2) +
  geom_text(data = headers, aes(x = x, y = 8, label = label),
            family = figure_font, fontface = "bold", size = 6.5 / .pt, hjust = 1) +
  coord_cartesian(xlim = c(0, 1.03), ylim = c(-0.5, 9), expand = FALSE, clip = "off") +
  theme_void(base_size = 6.5, base_family = figure_font) +
  theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5))

# combine
# one y-axis title for the three scatter plots
y_title = textGrob(paste0("MAG (", mag_id, ") relative abundance"),
                   x = unit(0.75, "npc"),
                   rot = 90, gp = gpar(fontsize = 6.5, fontfamily = figure_font))
plot_grid = wrap_plots(plots[[1]], plots[[2]], plots[[3]], wrap_elements(full = p_table), ncol = 2)
plot_grid = plot_grid & theme(plot.background = element_rect(fill = NA, colour = NA))
figure_plot = wrap_elements(full = y_title, clip = FALSE) + plot_grid +
  plot_layout(widths = c(0.035, 1))

# page size and margins
figure_width = 10
figure_height = 9.5
figure_margin = margin(2, 2, 2, 2, unit = "mm")
figure_plot = figure_plot +
  plot_annotation(theme = theme(plot.margin = figure_margin))

# save figure
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(figure_dir, "f5_correlation_plots.png"), figure_plot, device = agg_png,
       width = figure_width, height = figure_height,
       units = "cm", dpi = 600, bg = "white")

# save svg
ggsave(file.path(figure_dir, "f5_correlation_plots.svg"), figure_plot,
       device = svglite, fix_text_size = FALSE,
       width = figure_width, height = figure_height,
       units = "cm", bg = "white")

# keep each svg title together so editors preserve the spacing
svg_file = file.path(figure_dir, "f5_correlation_plots.svg")
svg = readLines(svg_file)
for (i in seq_along(plasmid_ids)) {
  title_line = grep(paste0(">", plasmid_ids[i], "</text>"), svg, fixed = TRUE)
  if (length(title_line) != 1L) next
  svg[title_line] = sub("</text>",
    paste0(" <tspan font-weight='bold'>(P", i, ")</tspan></text>"),
    svg[title_line], fixed = TRUE)
  svg = svg[-c(title_line + 1, title_line + 2)]
}
writeLines(svg, svg_file)
