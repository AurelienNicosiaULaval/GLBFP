library(GLBFP)

message("No leave-one-out API is currently exported by GLBFP.")
message("This placeholder records that LOO benchmarks are not applicable until such an API exists.")

available_functions <- getNamespaceExports("GLBFP")
loo_functions <- grep("loo|leave", available_functions, ignore.case = TRUE, value = TRUE)

if (length(loo_functions) == 0L) {
  message("Detected LOO-related exports: none.")
} else {
  message("Detected LOO-related exports: ", paste(loo_functions, collapse = ", "))
}
