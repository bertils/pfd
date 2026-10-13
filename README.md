# Nation-scale recovery of complete environmental plasmids pushes the boundaries of plasmid diversity through long-read sequencing

This repository contains input data tables, R analysis code, and figure files for the paper.

## Repository structure

- **`data/`** — input tables and supporting files. Data-processing methods are described in the manuscript; upstream processing pipelines are not included here.
- **`analysis/`** — R scripts for the analyses and figure exports.
- **`figures/`** — figure exports, using the paper's numbering. All analysis scripts save their figures here.

## Finding analyses and figures

Filenames follow the order and numbering of the paper. Files beginning with `si_` correspond to supplementary figures.

| /analysis | /figures |
| --- | --- |
| `analysis/f1_curation.R`, `analysis/f1_denmark_map.R`, `analysis/f1_venn_diagram.R` | `figures/f1.png` and the `f1_*` component exports |
| `analysis/f2.Rmd` | `figures/f2.png` |
| `analysis/f3.Rmd` | `figures/f3.png` |
| `analysis/f4.R` | `figures/f4.png` |
| `analysis/f5_correlation_plots.R` | `figures/f5_correlation_plots.png` |
| Matching `analysis/si_f1` through `analysis/si_f12` R scripts | Matching `figures/si_f1.png` through `figures/si_f12.png` |

**Figure 5B was made manually, so no code or figure file for that panel is included in this repository.**

The trait summaries used in Figure 2D are computed in `analysis/si_f5.Rmd` and supplied as `data/traits_all.tsv`. The supplementary correlation-table calculation is included in `analysis/f5_correlation_plots.R`. Figure 3's supporting input files are stored directly in `data/`.

## Large input tables omitted from GitHub

The following files exceed GitHub's 100 MiB file-size limit and are not included in this repository. Sizes are in GB (1,000,000,000 bytes).

| File | Size (GB) |
| --- | ---: |
| `data/orfs_mfd_and_reference.tsv` | 2.9 |
| `data/IMGPR_plasmid_data.tsv` | 0.2 |
| `data/merged.faa-ripdb.ultra.sensitive-with-real-names.tsv` | 0.2 |
| `data/raw_illumina_plasmid_prediction.tsv` | 0.3 |
| `data/sylph.tsv` | 3.7 |

These files can be found at [10.5281/zenodo.22944851](https://doi.org/10.5281/zenodo.22944851). Download the required tables into `data/` before running analyses that use them.

If you use circure, please cite **doi:** [10.64898/2026.10.06.756670](https://doi.org/10.64898/2026.10.06.756670).
