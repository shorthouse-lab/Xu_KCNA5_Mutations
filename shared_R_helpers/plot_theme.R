## ---------------------------------------------------------------------------
## Shared ggplot theme + save-both-formats helper for publication figures.
## ---------------------------------------------------------------------------
library(ggplot2)

theme_kcna5 <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90", colour = NA),
      plot.title = element_text(face = "bold", size = base_size + 1),
      legend.position = "right"
    )
}

## Saves a ggplot (or a base plotting function) as both PDF and 300 DPI PNG
## to results/figures/<analysis>/<name>.{pdf,png}
save_fig <- function(plot, name, analysis, width = FIG_WIDTH_IN, height = FIG_HEIGHT_IN) {
  out_dir <- file.path(RESULTS_DIR, "figures", analysis)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  pdf_path <- file.path(out_dir, paste0(name, ".pdf"))
  png_path <- file.path(out_dir, paste0(name, ".png"))

  if (inherits(plot, "ggplot") || inherits(plot, "gg")) {
    ggsave(pdf_path, plot, width = width, height = height, device = "pdf")
    ggsave(png_path, plot, width = width, height = height, dpi = FIG_DPI, device = "png")
  } else if (is.function(plot)) {
    pdf(pdf_path, width = width, height = height); plot(); dev.off()
    png(png_path, width = width, height = height, units = "in", res = FIG_DPI); plot(); dev.off()
  } else {
    stop("save_fig: `plot` must be a ggplot object or a zero-arg plotting function")
  }
  invisible(list(pdf = pdf_path, png = png_path))
}

## Writes a results table as both .csv and .rds to results/tables/<analysis>/<name>.*
save_table <- function(df, name, analysis) {
  out_dir <- file.path(RESULTS_DIR, "tables", analysis)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(df, file.path(out_dir, paste0(name, ".csv")))
  saveRDS(df, file.path(out_dir, paste0(name, ".rds")))
  invisible(df)
}

## Writes a short plain-English markdown summary to results/summaries/
save_summary <- function(text, name, analysis) {
  out_dir <- file.path(RESULTS_DIR, "summaries")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  writeLines(text, file.path(out_dir, paste0(analysis, "_", name, ".md")))
}
