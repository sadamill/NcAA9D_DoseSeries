library(tidyverse)

source("global_functions.R")

doses <- calculate_dose("fwd")

readr::write_csv(doses, "1-fwd_calculation/output/r_output/fwds.csv")
