#!/usr/bin/env Rscript
# Render the official manuscript report from the supplied PSAMs.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 1L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
  setwd(dirname(dirname(script_path)))
}
if (!file.exists("Figures.Rmd")) stop("Run this script from a checkout containing Figures.Rmd.")
suppressPackageStartupMessages(library(FamilyCode))
if (nzchar(Sys.getenv("RSTUDIO_PANDOC"))) {
  pandoc_dir <- Sys.getenv("RSTUDIO_PANDOC")
  pandoc_binary <- file.path(pandoc_dir,
    if (.Platform$OS.type == "windows") "pandoc.exe" else "pandoc")
  if (!file.exists(pandoc_binary)) stop("RSTUDIO_PANDOC does not contain a Pandoc executable.")
  pandoc_version <- system2(pandoc_binary, "--version", stdout = TRUE)[1L]
  rmarkdown:::set_pandoc_info(pandoc_dir,
    numeric_version(sub("^pandoc ", "", pandoc_version)))
}
if (!rmarkdown::pandoc_available()) {
  stop("Pandoc is required. Install Pandoc or set RSTUDIO_PANDOC to RStudio's Pandoc directory.")
}
render_env <- new.env(parent = globalenv())
rmarkdown::render("Figures.Rmd", clean = TRUE, envir = render_env, encoding = "UTF-8")
panel_map <- read.csv("Rfunctions/manuscript_panels.csv", stringsAsFactors = FALSE)
stopifnot(nrow(panel_map) == 60L, sum(panel_map$Kind == "scatter") == 28L,
  length(render_env$scatter_panel_results) == 28L,
  all(file.exists(c(panel_map$PNG, panel_map$PDF))),
  all(render_env$projection_audit$validation$Passed))
writeLines(capture.output(sessionInfo()), "outputs/manuscript_panel_validation/render_session_info.txt")
message("Official report rendered: Figures.html, Figures.md, 60 panel exports and numerical source data.")
