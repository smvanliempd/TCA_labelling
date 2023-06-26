library(dplyr)
library(ggplot2)

# get isotopomer list
isotopomers <- cycle_rec(akg_init = akg_00,
                         akg_cyc = akg_00,
                         p_inj = prop_inj, 
                         ac = AcCoA) |>
  iso_table()

# plot mass changes per cycle
iso_plot(isotopomers)



