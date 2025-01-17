library(tidyverse)
library(did2s)
library(glue)
library(here)

outcome <- "retail"
sample <- read_csv(here("data/County_Business_Patterns/sample_basker_YEARS_1977_1999_T0_1985.csv"))

long_diffs <- sample |> 
  reframe(
    .by = fips,
    total_pop = total_pop[1],
    share_school_col = share_school_col[1],
    Delta_log_retail_emp = log_retail_emp[year == 1985] - log_retail_emp[year == 1977]
  )

ggplot(long_diffs, aes(x = share_school_col, y = Delta_log_retail_emp)) + 
  geom_point(
    aes(size = total_pop), alpha = 0.2
  ) + 
  geom_smooth(
    aes(weight = sqrt(total_pop)),
    method = "lm"
  ) +
  geom_smooth(
    aes(weight = total_pop),
    method = "lm"
  ) +
  kfbmisc::theme_kyle(base_size = 14)


