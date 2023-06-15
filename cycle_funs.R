# cycle functions
filter_dups <- function(ent_fwd, ent_rev, e) {
  if (all(ent_fwd == ent_rev)) {
    ent <- ent_fwd
    ent_reg <- c(sum(ent), e$reg[2], strtoi(paste0(ent_fwd,collapse = ""), base = 2))
    return(list(list(ent = ent,
                     reg = ent_reg)))
  } else {
    ent_fwd_reg <- c(sum(ent_fwd), e$reg[2]/2, strtoi(paste0(ent_fwd,collapse = ""), base = 2) )
    ent_rev_reg <- c(sum(ent_rev), e$reg[2]/2, strtoi(paste0(ent_rev,collapse = ""), base = 2) )
    return(list(list(ent = ent_fwd,
                     reg = ent_fwd_reg),
                list(ent = ent_rev,
                     reg = ent_rev_reg ) ))
  }
} 
icit_to_akg <- function(icit) {
  akg <- sapply(icit, function(i) {
    akg_fwd <- i$ent[1:5]
    akg_rev <- i$ent[5:1]
    filter_dups(akg_fwd,akg_rev,i)
  }, simplify = FALSE )
  akg <- unlist(akg,recursive = FALSE)
}
akg_to_suc <- function(akg) {
  suc <- sapply(akg, function(a) {
    suc <- a$ent[2:5]
    suc_reg <- c(sum(suc), a$reg[2],strtoi(paste0(suc,collapse = ""), base = 2)  )
    return(list(ent = suc,
                reg = suc_reg))
  },simplify = FALSE )
} 
suc_to_mal <- function(suc) {
  mal <- sapply(suc, function(s) {
    mal_fwd <- s$ent[1:4]
    mal_rev <- s$ent[4:1]
    filter_dups(mal_fwd,mal_rev,s)
  }, simplify = FALSE )
  mal <- unlist(mal,recursive = FALSE)
}
mal_to_icit <- function(mal, accoa = c(0,0) ) {
  icit <- sapply(mal, function(m) {
    icit <- c(accoa, m$ent[2:4], m$ent[1])
    icit_reg <- c(sum(icit),m$reg[2],strtoi(paste0(icit,collapse = ""), base = 2) )
    return(list(ent = icit,
                reg = icit_reg))
  },simplify = FALSE)
}

# recursive cycle function
cycle <- function(n, akg_init) {
  
  if (n == N_cycle) {
    cat("stop")
  } else {
    suc  <- akg_to_suc(akg = akg_init)
    mal  <- suc_to_mal(suc = suc)
    icit <- mal_to_icit(mal = mal)
    akg  <- icit_to_akg(icit = icit)
    
    # TODO clean up dups for akg
    
    L[[n]] <<- list(suc = suc,mal = mal,icit =  icit,akg =akg)
    cycle(n + 1, akg)
  }
  
}


