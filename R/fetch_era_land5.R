# ERA5-Land hourly summer temperature for the collection sites
#
# Queries the Open-Meteo historical archive API (models = era5_land) for 2 m air
# temperature at each collection site, June-August, 1950-2024. Writes the table
# read by R/fig1bcd_climate.R.
#
# Input:  data/collection_regions.csv  (site, family, where, lat, long)
# Output: data/era5_land_temperature_2m.csv
#
# The full output is ~500 MB and is not deposited. Requests are rate-limited, so
# a full run takes several hours.

library(dplyr)
library(httr)
library(jsonlite)

coords_csv <- "data/collection_regions.csv"
out_csv    <- "data/era5_land_temperature_2m.csv"

START_DATE <- "1950-06-01"
END_DATE   <- "2024-08-31"
VARS       <- "temperature_2m"

fetch_site <- function(site, family, where, lat, long) {
  url <- paste0(
    "https://archive-api.open-meteo.com/v1/archive?",
    "latitude=", lat,
    "&longitude=", long,
    "&start_date=", START_DATE,
    "&end_date=", END_DATE,
    "&hourly=", paste(VARS, collapse = ","),
    "&models=era5_land",
    "&timezone=auto"
  )
  
  res <- GET(url)
  if (status_code(res) != 200) {
    warning("API error for site ", site, ": status ", status_code(res))
    return(NULL)
  }
  
  hourly <- fromJSON(content(res, as = "text"))$hourly
  if (length(hourly$time) == 0) {
    warning("No data returned for site ", site)
    return(NULL)
  }
  
  df <- data.frame(
    time   = as.POSIXct(hourly$time, format = "%Y-%m-%dT%H:%M"),
    site   = site,
    family = family,
    where  = where,
    lat    = lat,
    long   = long
  )
  for (v in VARS) df[[v]] <- hourly[[v]]
  
  df %>%
    filter(format(time, "%m") %in% c("06", "07", "08")) %>%
    mutate(time = format(time, "%Y-%m-%d %H:%M:%S"))
}

# --- run ---------------------------------------------------------------------
coords <- na.omit(read.csv(coords_csv))
stopifnot(all(c("site", "family", "where", "lat", "long") %in% names(coords)))

weather <- lapply(seq_len(nrow(coords)), function(i) {
  message("Fetching site ", coords$site[i], " (", i, " of ", nrow(coords), ")")
  
  df <- NULL
  for (attempt in 1:5) {
    Sys.sleep(runif(1, 15, 25))          # Open-Meteo rate limit
    df <- fetch_site(coords$site[i], coords$family[i], coords$where[i],
                     coords$lat[i], coords$long[i])
    if (!is.null(df)) break
    Sys.sleep(runif(1, 30, 60) + 30 * attempt)
  }
  df
}) %>% bind_rows()

write.csv(weather, out_csv, row.names = FALSE)
