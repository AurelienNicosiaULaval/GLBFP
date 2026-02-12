args <- commandArgs(trailingOnly = TRUE)

`%||%` <- function(x, y) if (is.null(x)) y else x

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- sub("^--file=", "", file_arg[1] %||% "paper/reproduce.R")
script_path <- normalizePath(script_path, mustWork = FALSE)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(repo_root)

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(repo_root, quiet = TRUE)
} else {
  library(GLBFP)
}

source(file.path("inst", "bench", "sim_scenarios.R"))
source(file.path("inst", "bench", "run_validation.R"))

res <- run_validation_suite(reps = 3, grid_size_1d = 80, grid_size_2d = 20)
summary_tbl <- summarise_validation(res)

summary_df <- summary_tbl

if (!dir.exists(file.path("paper", "results"))) {
  dir.create(file.path("paper", "results"), recursive = TRUE)
}

write.csv(summary_df, file.path("paper", "results", "validation_summary.csv"), row.names = FALSE)
write.csv(res$sensitivity, file.path("paper", "results", "sensitivity.csv"), row.names = FALSE)

cat("Wrote:\n")
cat("- paper/results/validation_summary.csv\n")
cat("- paper/results/sensitivity.csv\n")
