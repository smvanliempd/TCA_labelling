icit_init <- list(list(ent = c(0,0,1,1,1,1),
                       reg = c(4,0.5,1)),
                  list(ent = c(0,0,0,0,0,0),
                       reg = c(0,0.5,1)))
akg_init <- icit_to_akg(icit_init)
suc_init <- akg_to_suc(akg_init)
mal_init <- suc_to_mal(suc_init)
icit_fin <- mal_to_icit(mal_init)

# 2nd cycle
akg_2nd  <- icit_to_akg(icit_fin)
suc_2nd <- akg_to_suc(akg_2nd)
mal_2nd <- suc_to_mal(suc_2nd)

# cycle
N_cycle <- 10
L <- list()
akg_init <- list(list(ent = c(1,1,1,1,1),
                      reg = c(5,1,1)))
cycle(n =1 , akg_init = akg_init )
