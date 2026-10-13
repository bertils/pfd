# Analysis paths

Run the R scripts from the repository root, for example `source("analysis/f1_curation.R")`, or source them from `analysis/`. R Markdown notebooks can be opened and knitted from `analysis/`.

Each script and notebook explicitly defines `project_dir`, `data_dir`, and `figure_dir` at the top. Paths are resolved relative to the checkout, with inputs in `data/` and exports in `figures/`. Figure 3 inputs are stored directly in `data/`. Keep the shapefile and its companion files together in `data/CNTR_RG_01M_2024_4326.shp/`.

Install the packages listed at the beginning of the relevant script or notebook. Download any required oversized input tables from the data deposit before running the analysis; these are listed in the main README.

Figure 5 correlation plots use Inter when installed and otherwise use the system sans-serif font.
