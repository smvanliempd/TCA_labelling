clean_isos <- function(l) {
  
  # Add props for identical isotopomers
  #  and delete repeted row
  x <- do.call(rbind, l) |>
    data.frame() |>
    group_by(X3) |>
    mutate(X2 = sum(X2)) |>
    distinct() |>
    ungroup()
  return(x)
}
akg_to_suc <- function(akg) {
  
  # From akg to suc by deleting C1
  suc <- akg
  suc[,4] <- NA  #del
  suc <- apply(suc, 1, function(s) {
    s[1] <- sum(s[5:8])
    s[3] <- strtoi(paste0(s[5:8],collapse = ""), base = 2)
    return(s)
  }, simplify = FALSE)
  
  suc <- clean_isos(suc)
  return(suc)
}
suc_to_mal <- function(suc) {
  
  # from suc to mal by return original and 
  #  mirrored isotopomer. Props divided by 2
  mal <- apply(suc, 1, function(s) {
    
    # forward seq
    mal_fwd <- s
    mal_fwd[2] <- mal_fwd[2]/2
    names(mal_fwd) <- paste0("X",1:9)
    
    # reverse seq
    mal_rev <- s[c(1:3,9:4)] #mir
    mal_rev[2] <- mal_rev[2]/2
    mal_rev[3] <- strtoi(paste0(mal_rev[5:8],collapse = ""), base = 2)
    names(mal_rev) <- paste0("X",1:9)
    
    mal <- rbind(mal_fwd,mal_rev)
    return(mal)
  }, simplify = FALSE )
  
  mal <- clean_isos(mal)
  return(mal)
}
mal_to_icit <- function(mal, ac ) {
  
  # from mal to icit via cit by 1) moving Cs 2) adding Ac
  #  3) mirroring isotopomers and divide prop by 2 
  icit <- apply(mal, 1, function(m) {
    
    cit <- m[c(1:4,9,6:8,5)] # mov
    cit[4:5] <- ac           # add
    
    # forward seq
    icit_fwd <- cit 
    icit_fwd[1] <- sum(icit_fwd[4:9])
    icit_fwd[2] <- icit_fwd[2]/2
    icit_fwd[3] <- strtoi(paste0(icit_fwd[4:9],collapse = ""), base = 2)
    names(icit_fwd) <- paste0("X",1:9)
    
    # reverse seq
    icit_rev <- cit[c(1:3,8:4,9)] # mir
    icit_rev[1] <- sum(icit_rev[4:9])
    icit_rev[2] <- icit_rev[2]/2
    icit_rev[3] <- strtoi(paste0(icit_rev[4:9],collapse = ""), base = 2)
    names(icit_fwd) <- paste0("X",1:9)
    
    icit <- rbind(icit_fwd,icit_rev)
    return(icit)
  }, simplify = FALSE )
  
  icit <- clean_isos(icit)
  return(icit)
}
icit_to_akg <- function(icit) {
  
  # from icit to akg by deleting C6
  akg <- apply(icit, 1, function(i) {
    akg <- i
    akg[9] <- NA    # del
    akg[1] <- sum(akg[4:8])
    akg[3] <- strtoi(paste0(akg[4:8],collapse = ""), base = 2)
    return(akg)
  },simplify = FALSE)
  
  akg <- clean_isos(akg)
  return(akg)
}
inject_akg <- function(akg_cyc, akg_inj, p_inj) {
  
  # inject Gln in sequence via akg
  
  # get akg isotopomers formed from icit
  akg_cyc <- akg_cyc|>
    mutate(X2 = X2 * (1-p_inj) )
  
  # get injected isotopomers
  akg_inj <- apply(akg_inj, 1, function(a) {
    a[2] <- a[2] * p_inj
    a[3] <- strtoi(paste0(a[4:8],collapse = ""), base = 2)
    return(a)
  },simplify = FALSE)
  akg_inj <- clean_isos(akg_inj)
  
  akg <- rbind(akg_cyc,akg_inj) 
  return(akg)
}
cycle_rec <- function(n = 1, akg_init, akg_cyc, L = list(), p_inj, ac ) {
  
  # recursive cycling function
  if(n == N_cycle + 1) {
    cat("done\n")
    return(L)
  } else {
    suc <- akg_to_suc(akg_cyc)
    mal <- suc_to_mal(suc)
    icit <- mal_to_icit(mal, ac = ac)
    akg <- icit_to_akg(icit)
    akg_fin <- inject_akg(akg_cyc = akg,
                          akg_inj = akg_init, 
                          p_inj   = p_inj)
    L[[n]] <- list(suc  = suc,
                   mal  = mal,
                   icit = icit,
                   akg  = akg_fin)
    cycle_rec(n = n + 1,
              akg_init = akg_init,
              akg_cyc = akg_fin,
              L = L,
              p_inj = p_inj,
              ac = ac )
  }
}
iso_table <- function(iso_list) {
  
  # extract isotopomers and put them
  #  in a table for further analysis
  mets <- c("suc","mal","icit","akg")
  isos <- sapply(seq_along(iso_list), function(i) { 
    ll <- sapply(mets, function(m) {
      l <- iso_list[[i]][[m]] 
      l$cycle <- i
      l$metab <- m
      return(l)
    },simplify = FALSE)
    ll <- do.call(rbind,ll)
    return(ll)
  },simplify = FALSE )
  isos <- do.call(rbind, isos)
  names(isos) <- c("delta_mass", "prop", "id", paste0("C",1:6),"cycle","metab")
  isos$metab <- factor(isos$metab, mets)
  return(isos)
}

iso_plot <- function(iso_table, N) {
  
  # plot porportions of equal-mass isotopomers
  #  note that different isotopomers can have
  #  same mass 
  isos_short <- iso_table |>
    filter(metab != "suc") |>
    group_by(delta_mass, cycle, metab) |>
    select(delta_mass, prop, cycle, metab) |>
    mutate(prop = sum(prop)) |>
    distinct() |>
    ungroup()
  brks <- unlist(sapply(seq(0,N, by = 4), function(i) c(i,rep("",3)),simplify = F ))
  p <- ggplot(isos_short,
              aes(
                x = cycle,
                y = prop,
                col = factor(delta_mass)
              )) +
    geom_line()+
    scale_color_brewer(name = "ΔM",palette = "Dark2")+
    scale_x_continuous(breaks = 0:(length(brks)-1), 
                       minor_breaks = NULL, 
                       labels =  brks)+
    scale_y_continuous(breaks =  seq(0,1,0.2) ) +
    facet_wrap(metab~.) +
    labs(x = "Cycle", y = "Proportion")+
    theme_bw() +
    theme()
  ggsave("TCA_label_disp.png",dpi = 300,
         height = 3, width = 6)
  return(p)
}
