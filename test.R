library(dplyr)

akg_00 <- rbind(
  # c(5,0.025, NA,1,1,1,1,0,NA),
  # c(4,0.025,NA,0,1,1,1,1,NA),
  c(5,1,NA,1,1,1,1,1,NA)) |>
  data.frame( ) |>
  tibble()
prop_inj <- 0.9
N_cycle <- 10

isotopomers <- cycle_rec(akg_init = akg_00, akg_cyc = akg_00,p_inj = prop_inj)
isotopomers
