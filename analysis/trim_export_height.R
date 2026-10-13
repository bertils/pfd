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

# Remove only blank rows above and below an export; keep its full width and DPI.
trim_export_height = function(filename) {
  image = png::readPNG(filename, native = TRUE, info = TRUE)
  dpi = attr(image, "info")$dpi
  # Native raster pixels are stored row by row.
  pixels = matrix(as.integer(image), nrow = nrow(image), ncol = ncol(image),
                  byrow = TRUE)
  rows = which(rowSums(pixels != -1L) > 0)
  if (!length(rows)) return(invisible(NULL))
  top = min(rows)
  bottom = max(rows)
  original_height = nrow(image)
  cropped = as.integer(t(pixels[top:bottom, , drop = FALSE]))
  dim(cropped) = c(bottom - top + 1L, ncol(image))
  class(cropped) = class(image)
  png::writePNG(cropped, filename, dpi = dpi)

  # Apply the same vertical crop to the SVG, when one was exported.
  svg_file = sub("\\.png$", ".svg", filename)
  if (file.exists(svg_file)) {
    svg = readLines(svg_file, warn = FALSE)
    root = which(grepl("^<svg ", svg))[1]
    match = regexec("viewBox=['\"]([^'\"]+)['\"]", svg[root])
    viewbox = as.numeric(strsplit(regmatches(svg[root], match)[[1]][2], "\\s+")[[1]])
    viewbox[2] = viewbox[2] + (top - 1) / original_height * viewbox[4]
    viewbox[4] = nrow(cropped) / original_height * viewbox[4]
    svg[root] = sub("viewBox=['\"][^'\"]+['\"]",
                    paste0("viewBox='", paste(viewbox, collapse = " "), "'"), svg[root])
    svg[root] = sub("height=['\"][^'\"]+['\"]",
                    paste0("height='", viewbox[4], "pt'"), svg[root])
    writeLines(svg, svg_file)
  }
  message(basename(filename), ": removed ", top - 1, " top and ",
          original_height - bottom, " bottom rows; width and DPI preserved.")
}
