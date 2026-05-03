# River Ashuapmushuan daily flow and level data

Daily observations of river flow and level for the Ashuapmushuan river.

## Usage

``` r
ashua
```

## Format

A data frame with 4,389 rows and 3 variables:

- flow:

  Flow rate in cubic meters per second.

- level:

  Water level in meters.

- day:

  Day code as integer in `YYYYDDD` format.

## Source

Environment and Climate Change Canada, Historical Hydrometric Data. The
exact extraction query still needs to be documented in `data-raw/`.

## Details

Data cover 22 March 1992 to 30 September 2007 with a small number of
missing calendar days.

## Examples

``` r
data(ashua)
summary(ashua)
#>       flow            level            day         
#>  Min.   :  39.5   Min.   :29.33   Min.   :1992081  
#>  1st Qu.:  95.1   1st Qu.:30.00   1st Qu.:1995289  
#>  Median : 192.0   Median :30.31   Median :2001004  
#>  Mean   : 249.0   Mean   :30.42   Mean   :1999735  
#>  3rd Qu.: 318.0   3rd Qu.:30.67   3rd Qu.:2004063  
#>  Max.   :1530.0   Max.   :33.29   Max.   :2007273  
```
