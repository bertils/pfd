PFD manuscript: data, analysis, and figures

This repository contains the data tables, R code, and figure files used for the figures and statistics in the PFD manuscript.

- data/ contains input tables used by the analyses. The data-processing steps are described in the manuscript Methods; those pipelines are not included here.
- analysis/ contains the R and R Markdown code for figures and statistical analyses. The code for the Figure 2C traits is in si_f5.Rmd. The code for the supplementary correlation table is in f5_correlation_plots.R.
- figures/ contains the resulting figure files.

The following large input tables are not stored in this GitHub repository. They are still in the local data/ folder and are intended for a separate data deposit, such as Zenodo:

- data/251201_sylph.tsv
- data/260515_new_eggnog_df.tsv
- data/260921_sylph.tsv
- data/IMGPR_plasmid_data.tsv
- data/raw_illumina_plasmid_prediction.tsv
- data/sylph_260921.tsv

Zenodo dataset link: to be added after deposit. Download these tables into data/ before running analyses that use them.
