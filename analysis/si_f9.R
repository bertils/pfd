data_dir = if (dir.exists("data")) "data" else "../data"
figure_dir = file.path(dirname(data_dir), "figures")

library(data.table)
library(svglite)
library(ggplot2)
library(dplyr)
library(tidyr)


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
    strip.text = element_text(size = 10),
    strip.background = element_blank(),
    legend.position = "none"
  )

# plasmids and mag
mag_id = "MFD01716.bin.1.47"     # all plasmids are correlated against THIS MAG

plasmid_ids = c(
  "contig_90854__MFD01716",
  "contig_33344__MFD04752",
  "contig_57306__MFD09567",
  "contig_42799__MFD05580",
  "contig_91650__MFD01716",
  "contig_51694__MFD09567",
  "contig_103805__MFD04752"
)

# single neutral colour for all facets (no per-plasmid palette)
point_col = "black"

# load dfs
sylph_mp = fread(file.path(data_dir, "sylph.tsv"),
                 select = c("sample", "representative", "rel_abund", "origin", "bin"))
df_mp = sylph_mp %>%
  filter(is.finite(rel_abund),
         (origin == "mag" & bin == mag_id) |
         (origin == "plasmid" & representative %in% plasmid_ids))

# match abundances by sample
mag_df = df_mp %>%
  filter(origin == "mag", bin == mag_id) %>%
  group_by(sample) %>%
  summarise(MAG = sum(rel_abund), .groups = "drop")

plasmid_df = df_mp %>%
  filter(origin == "plasmid", representative %in% plasmid_ids) %>%
  group_by(sample, representative) %>%
  summarise(Plasmid = sum(rel_abund), .groups = "drop")

all_samples = union(mag_df$sample, plasmid_df$sample)

dat = crossing(sample = all_samples, representative = plasmid_ids) %>%
  left_join(plasmid_df, by = c("sample", "representative")) %>%
  left_join(mag_df,     by = "sample") %>%
  mutate(
    Plasmid  = replace_na(Plasmid, 0),
    MAG      = replace_na(MAG, 0),
    category = case_when(
      Plasmid > 0 & MAG > 0  ~ "Both",
      Plasmid > 0 & MAG == 0 ~ "Plasmid only",
      MAG > 0 & Plasmid == 0 ~ "MAG only",
      TRUE                   ~ "Neither"
    )
  ) %>%
  filter(category != "Neither")

dat$category       = factor(dat$category, levels = c("Both", "Plasmid only", "MAG only"))
dat$representative = factor(dat$representative, levels = plasmid_ids)   # facet order

both_all = dat %>%
  filter(category == "Both") %>%
  mutate(copy_number = Plasmid / MAG)

# pearson tests for samples where both are detected
stat_one = function(p) {
  d = both_all[both_all$representative == p, ]
  n = nrow(d)
  if (n >= 3 && sd(d$Plasmid) > 0 && sd(d$MAG) > 0) {
    ct = cor.test(d$Plasmid, d$MAG, method = "pearson")
    r  = unname(ct$estimate); p_val = ct$p.value
  } else { r = NA_real_; p_val = NA_real_ }
  data.frame(
    representative = p, r = r, p = p_val, n_both = n,
    n_plasmid_only = sum(dat$representative == p & dat$category == "Plasmid only"),
    n_mag_only     = sum(dat$representative == p & dat$category == "MAG only"),
    mean_copy = mean(d$copy_number), sd_copy = sd(d$copy_number),
    stringsAsFactors = FALSE
  )
}
stats = do.call(rbind, lapply(plasmid_ids, stat_one))

# statistics inside each panel
yr = range(dat$MAG)
ystep = 0.12 * diff(yr)
ytop = max(dat$MAG) + 6 * ystep
rng = dat %>% group_by(representative) %>%
  summarise(xa = min(Plasmid) + 0.03 * (max(Plasmid) - min(Plasmid)), .groups = "drop")
xa_map = setNames(rng$xa, as.character(rng$representative))

fmtp = function(x) if (is.na(x)) "NA" else formatC(x, format = "g", digits = 3)

lines_df = do.call(rbind, lapply(seq_len(nrow(stats)), function(i) {
  s = stats[i, ]
  cp = if (is.na(s$mean_copy)) "'Copies = NA'"
  else if (is.na(s$sd_copy))
    paste0("'Copies ('*mu*') = '*'", sprintf("%.2f", s$mean_copy), "'")
  else
    paste0("'Copies ('*mu %+-% sigma*') = '*'", sprintf("%.2f", s$mean_copy),
           "' %+-% '", sprintf("%.2f", s$sd_copy), "'")
  labs = c(
    paste0("italic(r) == ", if (is.na(s$r)) "'NA'" else sprintf("%.3f", s$r)),
    paste0("italic(P) == '", fmtp(s$p), "'"),
    paste0("italic(n)[p+m] == ", s$n_both),
    paste0("italic(n)[p] == ",   s$n_plasmid_only),
    paste0("italic(n)[m] == ",   s$n_mag_only),
    cp
  )
  data.frame(representative = s$representative,
             y = ytop - (0:5) * ystep, lab = labs, stringsAsFactors = FALSE)
}))
lines_df$x = xa_map[as.character(lines_df$representative)]
lines_df$representative = factor(lines_df$representative, levels = plasmid_ids)

# plot
ncol = 2
p_facet = ggplot(dat, aes(x = Plasmid, y = MAG)) +
  geom_smooth(data = both_all,
              method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.6, show.legend = FALSE) +
  geom_point(data = filter(dat, category != "Both"),
             colour = point_col,
             size = 1.2, alpha = 0.5) +
  geom_point(data = filter(dat, category == "Both"),
             colour = point_col,
             size = 1.7, alpha = 0.85) +
  geom_text(data = lines_df, aes(x = x, y = y, label = lab),
            parse = TRUE, hjust = 0, vjust = 1, size = 10 / .pt,
            colour = "black", inherit.aes = FALSE) +
  facet_wrap(~ representative, scales = "free_x", ncol = ncol, axes = "all_x") +
  scale_x_continuous(position = "bottom") +
  labs(x = "Plasmid relative abundance",
       y = paste0("MAG (", mag_id, ") relative abundance")) +
  figure_theme

# page size and margins
figure_width = 21
figure_height = 29.7
figure_margin = margin(0, 3, 0, 3, unit = "mm")
p_facet = p_facet + theme(plot.margin = figure_margin)

# save figures
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(figure_dir, "si_f9.png"), p_facet,
       width = figure_width, height = figure_height,
       units = "cm", dpi = 300, bg = "white")


# Trim blank rows at the top and bottom; preserve horizontal whitespace.
source(file.path(dirname(figure_dir), "analysis", "trim_export_height.R"))
trim_export_height(file.path(figure_dir, "si_f9.png"))
