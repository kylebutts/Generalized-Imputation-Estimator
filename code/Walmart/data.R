#' ---
#' title: "Create Walmart dataset"
#' ---
#' The sample criterion here is based on Emek Basker's paper

# %%
#| message: false
#| warning: false
library(tidyverse)
library(data.table)
library(glue)
library(here)
library(sf)

# 1964 CBP ---------------------------------------------------------------------
# From https://www.fpeckert.me/elmmss/
cbp1964 <- fread(here("raw-data/County_Business_Patterns/CBP1964.csv")) |>
  _[sic == "----", .(emp1964 = emp), by = .(fipstate, fipscty)] |>
  setnames(
    old = c("fipstate", "fipscty"),
    new = c("state_fips", "county_fips")
  ) |>
  _[county_fips != "999", ]

# CBP Imputed employment counts ------------------------------------------------
# Info: http://fpeckert.me/cbp/Imputed%20Files/efsy_readme_panel.txt
cbp <- fread(here("raw-data/County_Business_Patterns/efsy_panel_naics.csv"))

cbp[, let(
  naics2 = stringr::str_sub(naics12, 1, 2),
  naics4 = stringr::str_sub(naics12, 1, 4),
  naics6 = stringr::str_sub(naics12, 1, 6),
  v1 = NULL
)]
setnames(cbp, old = c("fipstate", "fipscty"), new = c("state_fips", "county_fips"))

# NAICS 2-digit codes 44 and 45 are "Retail Trade"
# Note that Eckert et. al. change the data so that all the naics industries (2-digit, 4-digit, etc. are mutually exclusive), so summing across does not double count.
cbp <- cbp |>
  _[
    ,
    let(
      retail = naics2 %in% c("44", "45"),
      # Following Basker
      manufacturing = naics2 %in% c("31"),
      wholesale = naics2 %in% c("42"),
      construction = naics2 %in% c("23"),
      agriculture = naics2 %in% c("11"),
      healthcare = naics2 %in% c("62"),
      accomodations = naics2 %in% c("72")
    )
  ] |>
  _[county_fips != 999, ]

## Retail -------------------------------------------------------––-------------
retail <- cbp |>
  _[,
    .(emp = sum(emp, na.rm = T)),
    by = .(state_fips, county_fips, year, retail)
  ]

# Make balanced panel
retail <- retail |>
  expand(nesting(state_fips, county_fips), year) |>
  left_join(retail, by = c("state_fips", "county_fips", "year")) |>
  as.data.table()

retail[is.na(emp), emp := 0]
retail <- retail[!is.na(retail), ]

# Pivot wider
retail |>
  _[, retail := ifelse(retail, "retail", "nonretail")]

retail <- pivot_wider(
  retail,
  id_cols = c("state_fips", "county_fips", "year"),
  names_from = "retail",
  names_glue = "{retail}_{.value}",
  values_from = c("emp")
)

retail |>
  as.data.table() |>
  _[is.na(retail_emp), retail_emp := 0] |>
  _[is.na(nonretail_emp), nonretail_emp := 0]

## Wholesale ----------------------------------------------------––-------------
wholesale <- cbp |>
  _[,
    .(emp = sum(emp, na.rm = T)),
    by = .(state_fips, county_fips, year, wholesale)
  ]

# Make balanced panel
wholesale <- wholesale |>
  expand(nesting(state_fips, county_fips), year) |>
  left_join(wholesale, by = c("state_fips", "county_fips", "year")) |>
  as.data.table()

wholesale[is.na(emp), emp := 0]
wholesale <- wholesale[!is.na(wholesale), ]

# Pivot wider
wholesale |>
  _[, wholesale := ifelse(wholesale, "wholesale", "nonwholesale")]

wholesale <- pivot_wider(
  wholesale,
  id_cols = c("state_fips", "county_fips", "year"),
  names_from = "wholesale",
  names_glue = "{wholesale}_{.value}",
  values_from = c("emp")
)

wholesale |>
  as.data.table() |>
  _[is.na(wholesale_emp), wholesale_emp := 0] |>
  _[is.na(nonwholesale_emp), nonwholesale_emp := 0]

## Manufacturing ------------------------------------------------––-------------
manufacturing <- cbp |>
  _[,
    .(emp = sum(emp, na.rm = T)),
    by = .(state_fips, county_fips, year, manufacturing)
  ]

# Make balanced panel
manufacturing <- manufacturing |>
  expand(nesting(state_fips, county_fips), year) |>
  left_join(manufacturing, by = c("state_fips", "county_fips", "year")) |>
  as.data.table()

manufacturing[is.na(emp), emp := 0]
manufacturing <- manufacturing[!is.na(manufacturing), ]

# Pivot wider
manufacturing |>
  _[, manufacturing := ifelse(manufacturing, "manufacturing", "nonmanufacturing")]

manufacturing <- pivot_wider(
  manufacturing,
  id_cols = c("state_fips", "county_fips", "year"),
  names_from = "manufacturing",
  names_glue = "{manufacturing}_{.value}",
  values_from = c("emp")
)

manufacturing |>
  as.data.table() |>
  _[is.na(manufacturing_emp), manufacturing_emp := 0] |>
  _[is.na(nonmanufacturing_emp), nonmanufacturing_emp := 0]

## Construction -------------------------------------------------––-------------
construction <- cbp |>
  _[,
    .(emp = sum(emp, na.rm = T)),
    by = .(state_fips, county_fips, year, construction)
  ]

# Make balanced panel
construction <- construction |>
  expand(nesting(state_fips, county_fips), year) |>
  left_join(construction, by = c("state_fips", "county_fips", "year")) |>
  as.data.table()

construction[is.na(emp), emp := 0]
construction <- construction[!is.na(construction), ]

# Pivot wider
construction |>
  _[, construction := ifelse(construction, "construction", "nonconstruction")]

construction <- pivot_wider(
  construction,
  id_cols = c("state_fips", "county_fips", "year"),
  names_from = "construction",
  names_glue = "{construction}_{.value}",
  values_from = c("emp")
)

construction |>
  as.data.table() |>
  _[is.na(construction_emp), construction_emp := 0] |>
  _[is.na(nonconstruction_emp), nonconstruction_emp := 0]

## Agriculture --------------------------------------------------––-------------
agriculture <- cbp |>
  _[,
    .(emp = sum(emp, na.rm = T)),
    by = .(state_fips, county_fips, year, agriculture)
  ]

# Make balanced panel
agriculture <- agriculture |>
  expand(nesting(state_fips, county_fips), year) |>
  left_join(agriculture, by = c("state_fips", "county_fips", "year")) |>
  as.data.table()

agriculture[is.na(emp), emp := 0]
agriculture <- agriculture[!is.na(agriculture), ]

# Pivot wider
agriculture |>
  _[, agriculture := ifelse(agriculture, "agriculture", "nonagriculture")]

agriculture <- pivot_wider(
  agriculture,
  id_cols = c("state_fips", "county_fips", "year"),
  names_from = "agriculture",
  names_glue = "{agriculture}_{.value}",
  values_from = c("emp")
)

agriculture |>
  as.data.table() |>
  _[is.na(agriculture_emp), agriculture_emp := 0] |>
  _[is.na(nonagriculture_emp), nonagriculture_emp := 0]

## Health Care --------------------------------------------------––-------------
healthcare <- cbp |>
  _[,
    .(emp = sum(emp, na.rm = T)),
    by = .(state_fips, county_fips, year, healthcare)
  ]

# Make balanced panel
healthcare <- healthcare |>
  expand(nesting(state_fips, county_fips), year) |>
  left_join(healthcare, by = c("state_fips", "county_fips", "year")) |>
  as.data.table()

healthcare[is.na(emp), emp := 0]
healthcare <- healthcare[!is.na(healthcare), ]

# Pivot wider
healthcare |>
  _[, healthcare := ifelse(healthcare, "healthcare", "nonhealthcare")]

healthcare <- pivot_wider(
  healthcare,
  id_cols = c("state_fips", "county_fips", "year"),
  names_from = "healthcare",
  names_glue = "{healthcare}_{.value}",
  values_from = c("emp")
)

healthcare |>
  as.data.table() |>
  _[is.na(healthcare_emp), healthcare_emp := 0] |>
  _[is.na(nonhealthcare_emp), nonhealthcare_emp := 0]

## Accomodations ------------------------------------------------––-------------
accomodations <- cbp |>
  _[,
    .(emp = sum(emp, na.rm = T)),
    by = .(state_fips, county_fips, year, accomodations)
  ]

# Make balanced panel
accomodations <- accomodations |>
  expand(nesting(state_fips, county_fips), year) |>
  left_join(accomodations, by = c("state_fips", "county_fips", "year")) |>
  as.data.table()

accomodations[is.na(emp), emp := 0]
accomodations <- accomodations[!is.na(accomodations), ]

# Pivot wider
accomodations |>
  _[, accomodations := ifelse(accomodations, "accomodations", "nonaccomodations")]

accomodations <- pivot_wider(
  accomodations,
  id_cols = c("state_fips", "county_fips", "year"),
  names_from = "accomodations",
  names_glue = "{accomodations}_{.value}",
  values_from = c("emp")
)

accomodations |>
  as.data.table() |>
  _[is.na(accomodations_emp), accomodations_emp := 0] |>
  _[is.na(nonaccomodations_emp), nonaccomodations_emp := 0]

## Merge together datasets -----------------------------------------------------
retail <- retail |>
  merge(
    wholesale |> select(-nonwholesale_emp),
    by = c("state_fips", "county_fips", "year")
  ) |>
  merge(
    manufacturing |> select(-nonmanufacturing_emp),
    by = c("state_fips", "county_fips", "year")
  ) |>
  merge(
    construction |> select(-nonconstruction_emp),
    by = c("state_fips", "county_fips", "year")
  ) |>
  merge(
    agriculture |> select(-nonagriculture_emp),
    by = c("state_fips", "county_fips", "year")
  ) |>
  merge(
    healthcare |> select(-nonhealthcare_emp),
    by = c("state_fips", "county_fips", "year")
  ) |>
  merge(
    accomodations |> select(-nonaccomodations_emp),
    by = c("state_fips", "county_fips", "year")
  ) |>
  as.data.table()

## Merge in 1964 cbp -----------------------------------------------------------
retail <- merge(retail, cbp1964, by = c("state_fips", "county_fips"))
setorder(retail, "state_fips", "county_fips", "year")

retail <- retail |>
  as.data.table() |>
  _[,
    emp1977 := retail_emp[year == 1977] + nonretail_emp[year == 1977],
    by = c("state_fips", "county_fips")
  ]

# Load Walmart openings --------------------------------------------------------
openings <- here("raw-data/Walmart_openings/opening_dates.csv") |>
  fread() |>
  subset(!is.na(Latitude)) |>
  subset(!(BuildingStateProv %in% c("AK", "PR", "HI"))) |>
  st_as_sf(coords = c("Longitude", "Latitude"), crs = st_crs(4326)) |>
  st_transform(st_crs(2163)) |>
  as.data.table()

## Place openings within counties ----------------------------------------------

counties <- tigris::counties(cb = TRUE) |>
  st_transform(st_crs(2163)) |>
  as.data.table() |>
  _[, .(
    state_fips = STATEFP, county_fips = COUNTYFP,
    GEOID, geometry
  )]

# Sanity check: plot the points
# openings |>
#   st_geometry() |>
#   plot()

# Place walmarts into counties
openings$which_county <- lapply(
  st_within(openings$geometry, counties$geometry),
  \(x) if (length(x) == 0) NA else x[1]
) |> unlist()

# Manually place the one NA into Dare County, NC
openings[is.na(which_county), ]$which_county <- 1038

openings$fips <- counties[openings$which_county, .(fips = GEOID)]

# Collapse to county x year
openings <- openings[, .(
  store_number = StoreNbr,
  store_type = StoreType,
  store_zip = BuildingPostalCode,
  store_sqft = SizeSqFt,
  store_open_date = lubridate::mdy(OpenDate),
  store_open_year = lubridate::year(lubridate::mdy(OpenDate)),
  fips,
  state_fips = stringr::str_sub(fips, 1, 2),
  county_fips = stringr::str_sub(fips, 3, 5)
)]

counties_panel <- counties |>
  expand(
    nesting(state_fips, county_fips),
    year = 1975:2018
  ) |>
  left_join(counties, by = c("state_fips", "county_fips")) |>
  as.data.table()

# Many-to-many is okay! Multiple openings in a county
openings_panel <- left_join(
  counties_panel, openings,
  by = c("state_fips", "county_fips"),
  relationship = "many-to-many"
)

openings_panel <- openings_panel[,
  .(
    any_open = any(year >= store_open_year, na.rm = T),
    n_open = sum(year >= store_open_year, na.rm = T)
  ),
  by = .(state_fips, county_fips, year)
]

# Remove Alaska, Hawaii, and US Terrirotries
openings_panel <- openings_panel |>
  _[, let(
    state_fips = as.numeric(state_fips),
    county_fips = as.numeric(county_fips)
  )] |>
  _[!(state_fips %in% c(2, 15, 78, 72, 69, 66, 60)), ]

## Merge with CBP data ---------------------------------------------------------
retail <- inner_join(
  retail, openings_panel,
  by = c("state_fips", "county_fips", "year")
)

retail <- retail |>
  as.data.table() |>
  _[, fips := paste0(
    str_pad(state_fips, 2, side = "left", pad = "0"),
    str_pad(county_fips, 3, side = "left", pad = "0")
  )] |>
  _[,
    g := ifelse(
      any(any_open),
      year[min(which(any_open))] |> as.numeric(),
      Inf
    ),
    by = fips
  ] |>
  _[, g_anticipation_2 := g - 2] |>
  _[
    ,
    rel_year := ifelse(g == Inf, -Inf, year - g)
  ] |>
  _[
    ,
    rel_year_anticipation_2 :=
      ifelse(g_anticipation_2 == Inf, -Inf, year - g_anticipation_2)
  ]

retail |>
  _[, log_retail_emp := log(retail_emp + 1)] |>
  _[, log_nonretail_emp := log(nonretail_emp + 1)] |>
  _[, log_construction_emp := log(construction_emp + 1)] |>
  _[, log_wholesale_emp := log(wholesale_emp + 1)] |>
  _[, log_manufacturing_emp := log(manufacturing_emp + 1)] |>
  _[, log_agriculture_emp := log(agriculture_emp + 1)] |>
  _[, log_healthcare_emp := log(healthcare_emp + 1)] |>
  _[, log_accomodations_emp := log(accomodations_emp + 1)]

# Generate instrument: 1975 baseline retail and non-retail share
retail |>
  _[,
    retail_emp_share :=
      retail_emp[year == 1975] /
        (retail_emp[year == 1975] + nonretail_emp[year == 1975]),
    by = "fips"
  ] |>
  _[
    is.nan(retail_emp_share), retail_emp_share := 0
  ] |>
  _[,
    nonretail_emp_share :=
      nonretail_emp[year == 1975] /
        (retail_emp[year == 1975] + nonretail_emp[year == 1975]),
    by = "fips"
  ] |>
  _[
    is.nan(nonretail_emp_share), nonretail_emp_share := 0
  ]

# Merge with 1980 census data --------------------------------------------------
# Merge 3 datasets
{
  census1 <- fread(here("raw-data/1980_census/nhgis_ds104_1980_county.csv"))
  census2 <- fread(here("raw-data/1980_census/nhgis_ds107_1980_county.csv"))
  census3 <- fread(here("raw-data/1980_census/nhgis_ds110_1980_county.csv"))

  census1[
    ,
    fips := paste0(str_pad(STATEA, 2, "left", "0"), str_pad(COUNTYA, 3, "left", "0"))
  ]
  census2[
    ,
    fips := paste0(str_pad(STATEA, 2, "left", "0"), str_pad(COUNTYA, 3, "left", "0"))
  ]
  census3[
    ,
    fips := paste0(str_pad(STATEA, 2, "left", "0"), str_pad(COUNTYA, 3, "left", "0"))
  ]
  vars <- c("GISJOIN", "YEAR", "FSTATUS", "REGIONA", "DIVISIONA", "STATE", "STATEA", "SMSAA", "COUNTY", "COUNTYA", "REGION", "DIVISION", "CTY_SUBA", "PLACEA", "TRACTA", "BLCK_GRPA", "BLOCKA", "EDINDA", "ENUMDISTA", "SCSAA", "URB_AREAA", "CDA", "AIANHHA", "MCDSEQNO", "SEA", "UATYPE", "PLACDESC", "CBD", "INDSUBR", "LONGITUD", "LATITUDE", "LANDAREA", "AREANAME", "SUPFLG01", "SUPFLG02", "SUPFLG03", "SUPFLG04", "SUPFLG05", "SUPFLG06", "SUPFLG07", "SUPFLG08", "SUPFLG09", "SUPFLG10", "SUPFLG11", "SUPFLG12", "SUPFLG13", "SUPFLG14", "SUPFLG15", "SUPFLG16", "SUPFLG17", "SUPFLG18", "SUPFLG19", "SUPFLG20", "SUPFLG21", "SUPFLG22", "SUPFLG23", "SUPFLG24", "SUPFLG25", "SUPFLG26", "SUPFLG27")
  census1[, vars := NULL, env = list(vars = I(vars))]
  census2[, vars := NULL, env = list(vars = I(vars))]
  census3[, vars := NULL, env = list(vars = I(vars))]

  census <- census1 |>
    merge(census2, by = "fips") |>
    merge(census3, by = "fips")

  rm(census1, census2, census3)
}

setcolorder(census, "fips")

# Renames only the columns that show up in the data
rename_subset <- function(df, cols = col_convert) {
  cols <- cols[cols %in% colnames(df)]
  data.table::setnames(df, old = cols, new = names(cols))
  return(df)
}

# Clean variable names
col_convert <- c(
  "fips" = "fips",
  "total_pop" = "C7L001",
  "urban_pop" = "C7M001",
  "suburban_pop" = "C7M002",
  "rural_pop" = "C7M003",
  "total_housing_units" = "C8Y001",
  "male_pop" = "C9C001",
  "female_pop" = "C9C002",
  "white_pop" = "C9D001",
  "black_pop" = "C9D002",
  "median_rent" = "C8O001",
  "n_school_no_hs" = "DHM001",
  "n_school_some_hs" = "DHM002",
  "n_school_hs" = "DHM003",
  "n_school_some_col" = "DHM004",
  "n_school_col" = "DHM005",
  "n_pop_ind_agriculture" = "DIA001",
  "n_pop_ind_construction" = "DIA002",
  "n_pop_ind_manuf_nondurable" = "DIA003",
  "n_pop_ind_manuf_durable" = "DIA004",
  "n_pop_ind_transportation" = "DIA005",
  "n_pop_ind_communication" = "DIA006",
  "n_pop_ind_wholesale" = "DIA007",
  "n_pop_ind_retail" = "DIA008",
  "n_pop_ind_finance" = "DIA009",
  "n_pop_ind_biz_services" = "DIA010",
  "n_pop_ind_entertainment" = "DIA011",
  "n_pop_ind_health_services" = "DIA012",
  "n_pop_ind_ed_services" = "DIA013",
  "n_pop_ind_other_services" = "DIA014",
  "n_pop_ind_public" = "DIA015",
  "n_pop_income_79_bin_1" = "DID001",
  "n_pop_income_79_bin_2" = "DID002",
  "n_pop_income_79_bin_3" = "DID003",
  "n_pop_income_79_bin_4" = "DID004",
  "n_pop_income_79_bin_5" = "DID005",
  "n_pop_income_79_bin_6" = "DID006",
  "n_pop_income_79_bin_7" = "DID007",
  "n_pop_income_79_bin_8" = "DID008",
  "n_pop_income_79_bin_9" = "DID009",
  "n_pop_income_79_bin_10" = "DID010",
  "n_pop_income_79_bin_11" = "DID011",
  "n_pop_income_79_bin_12" = "DID012",
  "n_pop_income_79_bin_13" = "DID013",
  "n_pop_income_79_bin_14" = "DID014",
  "n_pop_income_79_bin_15" = "DID015",
  "n_pop_income_79_bin_16" = "DID016",
  "n_pop_income_79_bin_17" = "DID017",
  "n_pop_poverty_78_above" = "DI8001",
  "n_pop_poverty_78_below" = "DI8002",
  "n_pop_emp_private" = "DR1001",
  "n_pop_emp_government" = "DR1002",
  "n_pop_emp_self_employed" = "DR1003",
  "n_pop_emp_unpaid_family" = "DR1004"
)

census <- census |>
  rename_subset(col_convert) |>
  dplyr::select(tidyselect::all_of(names(col_convert)))

census |>
  _[, let(
    share_school_no_hs = n_school_no_hs / total_pop,
    share_school_some_hs = n_school_some_hs / total_pop,
    share_school_hs = n_school_hs / total_pop,
    share_school_some_col = n_school_some_col / total_pop,
    share_school_col = n_school_col / total_pop,
    share_pop_ind_agriculture = n_pop_ind_agriculture / total_pop,
    share_pop_ind_construction = n_pop_ind_construction / total_pop,
    share_pop_ind_manuf_nondurable = n_pop_ind_manuf_nondurable / total_pop,
    share_pop_ind_manuf_durable = n_pop_ind_manuf_durable / total_pop,
    share_pop_ind_transportation = n_pop_ind_transportation / total_pop,
    share_pop_ind_communication = n_pop_ind_communication / total_pop,
    share_pop_ind_wholesale = n_pop_ind_wholesale / total_pop,
    share_pop_ind_retail = n_pop_ind_retail / total_pop,
    share_pop_ind_finance = n_pop_ind_finance / total_pop,
    share_pop_ind_biz_services = n_pop_ind_biz_services / total_pop,
    share_pop_ind_entertainment = n_pop_ind_entertainment / total_pop,
    share_pop_ind_health_services = n_pop_ind_health_services / total_pop,
    share_pop_ind_ed_services = n_pop_ind_ed_services / total_pop,
    share_pop_ind_other_services = n_pop_ind_other_services / total_pop,
    share_pop_ind_public = n_pop_ind_public / total_pop,
    share_pop_income_79_bin_1 = n_pop_income_79_bin_1 / total_pop,
    share_pop_income_79_bin_2 = n_pop_income_79_bin_2 / total_pop,
    share_pop_income_79_bin_3 = n_pop_income_79_bin_3 / total_pop,
    share_pop_income_79_bin_4 = n_pop_income_79_bin_4 / total_pop,
    share_pop_income_79_bin_5 = n_pop_income_79_bin_5 / total_pop,
    share_pop_income_79_bin_6 = n_pop_income_79_bin_6 / total_pop,
    share_pop_income_79_bin_7 = n_pop_income_79_bin_7 / total_pop,
    share_pop_income_79_bin_8 = n_pop_income_79_bin_8 / total_pop,
    share_pop_income_79_bin_9 = n_pop_income_79_bin_9 / total_pop,
    share_pop_income_79_bin_10 = n_pop_income_79_bin_10 / total_pop,
    share_pop_income_79_bin_11 = n_pop_income_79_bin_11 / total_pop,
    share_pop_income_79_bin_12 = n_pop_income_79_bin_12 / total_pop,
    share_pop_income_79_bin_13 = n_pop_income_79_bin_13 / total_pop,
    share_pop_income_79_bin_14 = n_pop_income_79_bin_14 / total_pop,
    share_pop_income_79_bin_15 = n_pop_income_79_bin_15 / total_pop,
    share_pop_income_79_bin_16 = n_pop_income_79_bin_16 / total_pop,
    share_pop_income_79_bin_17 = n_pop_income_79_bin_17 / total_pop,
    share_pop_poverty_78_above = n_pop_poverty_78_above / total_pop,
    share_pop_poverty_78_below = n_pop_poverty_78_below / total_pop,
    share_pop_emp_private = n_pop_emp_private / total_pop,
    share_pop_emp_government = n_pop_emp_government / total_pop,
    share_pop_emp_self_employed = n_pop_emp_self_employed / total_pop,
    share_pop_emp_unpaid_family = n_pop_emp_unpaid_family / total_pop
  )] |>
  _[, share_pop_ind_manuf := share_pop_ind_manuf_durable + share_pop_ind_manuf_nondurable]

setcolorder(retail, c("fips"))
retail <- merge(retail, census, by = "fips")

# Estimation sample ------------------------------------------------------------

# Matches Emek Basker
YEARS <- 1977:1999

# Only using units treated between 1986-1995
T0 <- min(YEARS) + 8

sample <- retail |>
  _[year %in% YEARS, ] |>
  _[,
    balanced := (length(year) == length(YEARS)),
    by = .(fips)
  ] |>
  _[balanced == TRUE, ] |>
  _[,
    balanced := (sort(year) == YEARS),
    by = .(fips)
  ] |>
  _[balanced == TRUE, ] |>
  # Fix g for the panel
  _[g > max(YEARS), g := Inf] |>
  _[g > T0] |>
  _[g == Inf, rel_year := -Inf]

sample_basker <- sample |>
  _[emp1964 >= 1500, ] |>
  _[emp1977 >= emp1964, ]

# N counties
sample_basker[, fips |> unique() |> length()]

# N counties with <= 1 walmart
# sample_basker[, temp := (max(n_open) <= 2), by = fips][temp == TRUE, fips |> unique() |> length()]

# Export -----------------------------------------------------------------------
fwrite(
  sample,
  "data/County_Business_Patterns/sample_YEARS_{min(YEARS)}_{max(YEARS)}_T0_{T0}.csv" |>
    glue() |> here()
)

fwrite(
  sample_basker,
  "data/County_Business_Patterns/sample_basker_YEARS_{min(YEARS)}_{max(YEARS)}_T0_{T0}.csv" |>
    glue() |> here()
)
