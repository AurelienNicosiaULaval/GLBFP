# Recreate data/ashua.rda from the original hydrometric source.
#
# The raw file and exact extraction query are not currently included in this
# repository. This script intentionally stops until those inputs are supplied.

stop(
  paste(
    "Raw data provenance for ashua is incomplete.",
    "Add the ECCC raw file, station identifier, download query, and cleaning steps before regenerating data/ashua.rda."
  ),
  call. = FALSE
)
