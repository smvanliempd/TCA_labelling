library(dplyr)
library(ggplot2)

akg_00 <- rbind(
  # c(5,0.025, NA,1,1,1,1,0,NA),
  # c(4,0.025,NA,0,1,1,1,1,NA),
  c(5,1,NA,1,1,1,1,1,NA)) |>
  data.frame( ) |>
  tibble()
prop_inj <- 0.5
N_cycle <- 10

isotopomers <- cycle_rec(akg_init = akg_00, akg_cyc = akg_00,p_inj = prop_inj)
mets <- c("suc","mal","icit","akg")
isos_expand <- sapply(seq_along(isotopomers), function(i) { 
  ll <- sapply(mets, function(m) {
    l <- isotopomers[[i]][[m]] 
    l$cycle <- i
    l$metab <- m
    return(l)
  },simplify = FALSE)
  ll <- do.call(rbind,ll)
  return(ll)
},simplify = FALSE )
isos_expand <- do.call(rbind, isos_expand)

isos_short <- isos_expand |>
  group_by(X1, cycle, metab) |>
  select(X1, X2, cycle, metab) |>
  mutate(X2 = sum(X2)) |>
  distinct() |>
  ungroup()

ggplot(isos_short,
       aes(
         x = cycle,
         y = X2,
         col = factor(X1)
       )) +
  geom_line()+
  facet_wrap(metab~.) +
  theme_bw()