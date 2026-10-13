# Analysis paths

Run the R scripts from the repository root, for example `source("analysis/f1_curation.R")`, or source them from `analysis/`. R Markdown notebooks can be opened and knitted from `analysis/`.

The shared `project_setup.R` file locates the repository and sets `data_dir` to `data/` and `figure_dir` to `figures/`. Figure 3 inputs are stored directly in `data/`. Keep the shapefile and its companion files together in `data/CNTR_RG_01M_2024_4326.shp/`.

Install the packages listed at the beginning of the relevant script or notebook. Download any required oversized input tables from the data deposit before running the analysis; these are listed in the main README.

Figure 5 correlation plots use Inter when installed and otherwise use the system sans-serif font.
