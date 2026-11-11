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

## CDEM window near Quebec City

`extract_cdem.R` acquires a small 16 m projected window from Natural Resources
Canada's public CDEM COG using GDAL range requests. Run it from the package
source root with terra installed. It verifies the CRS before transforming the
requested longitude-latitude center and cropping. The shipped CSV and provenance
notice are in `inst/extdata/`. No stationarity or mixing is asserted.

Contains information licensed under the Open Government Licence - Canada.
Source catalogue: https://open.canada.ca/data/en/dataset/7f245e4d-76c2-4caa-951a-45d1d2051333
Licence: https://open.canada.ca/en/open-government-licence-canada
