# PFD manuscript: data, analysis, and figures

This repository contains the input data tables, R analysis code, and figure files used for the figures and statistics in the PFD manuscript.

## Repository structure

- **`data/`** — input tables used by the analyses. Data-processing steps are described in the manuscript Methods; those processing pipelines are not included here.
- **`analysis/`** — R scripts and R Markdown documents for figures and statistical analyses.
- **`figures/`** — resulting figure files.

## Finding analyses and figures

Filenames follow the order and numbering of the paper. Match the numbered filename prefixes in `analysis/` and `figures/` to locate the code and corresponding figure. For example, `analysis/f2.Rmd` corresponds to `figures/f2.png`; files beginning with `si_` contain supplementary analyses and outputs. Original filenames and numbering are preserved.

Two specific cross-references:

- The analysis for the Figure 2C traits is in `analysis/si_f5.Rmd`.
- The code for the supplementary correlation table is in `analysis/f5_correlation_plots.R`.

## Large input tables omitted from GitHub

The following files exceed GitHub's 100 MiB limit for ordinary Git files and are not included in this repository. Sizes below are in MiB (1,048,576 bytes).

| File | Size (MiB) |
| --- | ---: |
| `data/251201_sylph.tsv` | 3,604.5 |
| `data/260515_new_eggnog_df.tsv` | 3,879.5 |
| `data/260921_sylph.tsv` | 181.7 |
| `data/IMGPR_plasmid_data.tsv` | 164.0 |
| `data/raw_illumina_plasmid_prediction.tsv` | 272.9 |
| `data/sylph_260921.tsv` | 3,662.2 |

The original README identifies these tables for a separate data deposit, such as Zenodo. A dataset link is not yet available. Once deposited, download the required tables into `data/` before running analyses that use them.

## Using the analyses

Open the relevant R or R Markdown file in `analysis/` for the figure or statistic of interest. Review its package requirements and input/output paths for your local setup. Some analyses require the large tables listed above.
