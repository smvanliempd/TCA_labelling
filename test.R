icit_init <- list(list(ent = c(0,0,1,1,1,1),
                       reg = c(4,0.5,1)),
                  list(ent = c(0,0,0,0,0,0),
                       reg = c(0,0.5,1)))
akg_init <- icit_to_akg(icit_init)
suc_init <- akg_to_suc(akg_init)
mal_init <- suc_to_mal(suc_init)
icit_fin <- mal_to_icit(mal_init)
