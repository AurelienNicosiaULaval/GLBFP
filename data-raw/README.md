# Data preparation notes

The package currently includes `data/ashua.rda`, documented in `R/ashua.R`.

The data are described as daily flow and level observations from Environment and
Climate Change Canada Historical Hydrometric Data. The raw source file used to
create `data/ashua.rda` is not included in the repository.

Before CRAN submission or journal publication, the raw-data provenance should be
completed with:

- station identifier;
- exact download URL or reproducible query;
- download date;
- raw file checksum;
- transformation script from raw data to `data/ashua.rda`.
