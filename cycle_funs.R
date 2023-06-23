clean_isos <- function(l) {
  x <- do.call(rbind, l) |>
    data.frame() |>
    group_by(X3) |>
    mutate(X2 = sum(X2)) |>
    distinct() |>
    ungroup()
  return(x)
}
akg_to_suc <- function(akg) {
  suc <- akg
  suc[,4] <- NA
  suc <- apply(suc, 1, function(s) {
    s[1] <- sum(s[5:8])
    s[3] <- strtoi(paste0(s[5:8],collapse = ""), base = 2)
    return(s)
  }, simplify = FALSE)
  suc <- clean_isos(suc)
  return(suc)
}
suc_to_mal <- function(suc) {
  mal <- apply(suc, 1, function(s) {

    # forward seq
    mal_fwd <- s
    mal_fwd[2] <- mal_fwd[2]/2
    names(mal_fwd) <- paste0("X",1:9)
    
    # reverse seq
    mal_rev <- s[c(1:3,9:4)]
    mal_rev[2] <- mal_rev[2]/2
    mal_rev[3] <- strtoi(paste0(mal_rev[5:8],collapse = ""), base = 2)
    names(mal_rev) <- paste0("X",1:9)
    
    # return
    rbind(mal_fwd,mal_rev)
    
  }, simplify = FALSE )
  mal <- clean_isos(mal)
  return(mal)
}
mal_to_icit <- function(mal, ac = c(0,0)) {
  icit <- apply(mal, 1, function(m) {

    # forward seq
    icit_fwd <- m[c(1:4,9,6:8,5)] # mov
    icit_fwd[4:5] <- ac           # add
    icit_fwd[2] <- icit_fwd[2]/2
    icit_fwd[3] <- strtoi(paste0(icit_fwd[4:9],collapse = ""), base = 2)
    names(icit_fwd) <- paste0("X",1:9)
    
    # reverse seq
    icit_rev <- m[c(1:4,9,7:5,8)] # mov
    icit_rev[4:5] <- ac           # add
    icit_rev[2] <- icit_rev[2]/2
    icit_rev[3] <- strtoi(paste0(icit_rev[4:9],collapse = ""), base = 2)
    names(icit_fwd) <- paste0("X",1:9)
    
    # return
    rbind(icit_fwd,icit_rev)
    
  }, simplify = FALSE )
  icit <- clean_isos(icit)
  return(icit)
}
icit_to_akg <- function(icit) {
  akg <- apply(icit, 1, function(i) {
    akg <- i
    akg[9] <- NA
    akg[1] <- sum(akg[4:8])
    akg[3] <- strtoi(paste0(akg[4:8],collapse = ""), base = 2)
    return(akg)
  },simplify = FALSE)
  akg <- clean_isos(akg)
  return(akg)
}
inject_akg <- function(akg_cyc, akg_inj, p_inj) {
  akg_cyc <- akg_cyc|>
    mutate(X2 = X2 * (1-p_inj) )
  
  akg_inj <- apply(akg_inj, 1, function(a) {
    a[2] <- a[2] * p_inj
    a[3] <- strtoi(paste0(a[4:8],collapse = ""), base = 2)
    return(a)
  },simplify = FALSE)
  akg_inj <- clean_isos(akg_inj)
  akg <- rbind(akg_cyc,akg_inj) 
  return(akg)
}
cycle_rec <- function(n = 1, akg_init, akg_cyc, L = list(), p_inj ) {
  if(n == N_cycle + 1) {
    cat("done")
    return(L)
  } else {
    suc <- akg_to_suc(akg_cyc)
    mal <- suc_to_mal(suc)
    icit <- mal_to_icit(mal)
    akg <- icit_to_akg(icit)
    akg_fin <- inject_akg(akg_cyc = akg,
                          akg_inj = akg_init, 
                          p_inj   = p_inj)
    L[[n]] <- list(suc  = suc,
                   mal  = mal,
                   icit = icit,
                   akg  = akg_fin)
    cycle_rec(n = n + 1, akg_init, akg_fin, L = L, p_inj = p_inj)
  }
}


