# Nation-scale recovery of complete environmental plasmids pushes the boundaries of plasmid diversity through long-read sequencing

This repository contains input data tables, R analysis code, and figure files for the paper.

## Repository structure

- **`data/`** — input tables and supporting files. Data-processing methods are described in the manuscript; upstream processing pipelines are not included here.
- **`analysis/`** — R scripts, R Markdown notebooks, and shared helpers for the analyses and figure exports. See [analysis instructions](analysis/README.md) for dependencies and commands.
- **`figures/`** — figure exports, using the paper's numbering. All analysis scripts save their figures here.

## Finding analyses and figures

Filenames follow the order and numbering of the paper. Files beginning with `si_` correspond to supplementary figures.

| Paper figure | Analysis | Figure files |
| --- | --- | --- |
| Figure 1 | `analysis/f1_curation.R`, `analysis/f1_denmark_map.R`, `analysis/f1_venn_diagram.R` | `figures/f1.png` and the `f1_*` component exports |
| Figure 2 | `analysis/f2.Rmd` | `figures/f2.png` |
| Figure 3 | `analysis/f3.Rmd` | `figures/f3.png` |
| Figure 4 | `analysis/f4.R` | `figures/f4.png` |
| Figure 5 correlation plots | `analysis/f5_correlation_plots.R` | `figures/f5_correlation_plots.png` |
| Supplementary Figures 1–12 | Matching `analysis/si_f1` through `analysis/si_f12` scripts or notebooks | Matching `figures/si_f1.png` through `figures/si_f12.png` |

**Figure 5B was made manually, so no code or figure file for that panel is included in this repository.**

The trait summaries used in Figure 2D are computed in `analysis/si_f5.Rmd` and supplied as `data/traits_all.tsv`. The supplementary correlation-table calculation is included in `analysis/f5_correlation_plots.R`. Figure 3's supporting input files are stored directly in `data/`.

## Large input tables omitted from GitHub

The following files exceed GitHub's 100 MiB file-size limit and are not included in this repository. Sizes are in MiB (1,048,576 bytes).

| File | Size (MiB) |
| --- | ---: |
| `data/251201_sylph.tsv` | 3,604.5 |
| `data/260515_new_eggnog_df.tsv` | 3,879.5 |
| `data/IMGPR_plasmid_data.tsv` | 164.0 |
| `data/merged.faa-ripdb.ultra.sensitive-with-real-names.tsv` | 237.1 |
| `data/raw_illumina_plasmid_prediction.tsv` | 272.9 |
| `data/sylph.tsv` | 3,512.2 |

These tables can be found at [Zenodo link]. Once deposited, download the required tables into `data/` before running analyses that use them.

The Zenodo link above is a placeholder until the deposit is available. `251201_sylph.tsv` is a legacy table and is not referenced by the current scripts.
