# labelled via 13C5-glutamine
akg_00 <- rbind(
  # c(5,0.025, NA,1,1,1,1,0,NA),
  # c(4,0.025,NA,0,1,1,1,1,NA),
  c(0,0.4,NA,0,0,0,0,0,NA),
  c(5,0.6,NA,1,1,1,1,1,NA)
) |>
  data.frame( ) |>
  tibble()
AcCoA <- c(0,0)
prop_inj <- 0.9
N_cycle <- 10

# labelled AcCoA via glucose
akg_00 <- rbind(c(0,1,NA,0,0,0,0,0,NA) ) |>
  data.frame() |>
  tibble()
AcCoA <- c(1,1)
prop_inj <- 0.1
N_cycle <- 20