# Download and unzip the EFSY panel data
library(here)

url <- "http://fpeckert.me/cbp/Imputed%20Files/efsy_panel_naics.csv.zip"
zipfile <- here("raw-data/County_Business_Patterns/efsy_panel_naics.csv.zip")

options(timeout = 60 * 20)
download.file(url, destfile = zipfile)
unzip(
  zipfile, 
  exdir = here("raw-data/County_Business_Patterns/")
)
