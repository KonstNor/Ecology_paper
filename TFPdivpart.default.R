###############################################################################################################################################
# TFPdivpart: R function to perform additive partitionning of taxonomic (Simpson index), functional and phylogenetic diversity (Rao index)
#             and null model testing
# The script integrates two functions: "Rao" by de Bello et al.(J Veg Sci, 2010), and "adipart" in the package vegan.
# Author: Konstantinia Dede, konstantina.dede@gmail.com
# 
#
# INPUTS:                                                                                 
#   - "y": matrix of abundances (c x s) of the s species for the c local communities
#   - "x": matrix with c rows as in matrix "y", columns coding the levels of sampling hierarchy. The number of groups within 
#          the hierarchy must decrease from left to right  
#   - "index": character, a diversity index to be calculated
#   - "weights": argument "unif" for uniform weights, "prop" for weighting proportional to sample abundances to use in 
#                weighted averaging of individual alpha values within strata of a given level of the sampling hierarchy 
#   - "relative": logical, if TRUE then alpha and beta diversity values are given relative to the value of gamma 
#   - "nsimul": number of permutations 
#   - "method": null model method: either a name (character string) of a method defined in make.commsim or a commsim function (vegan package)
#	- "fundiv": matrix (s x s) or dist object with pairwise functional trait distances between the s species
#   - "phyldiv": matrix (s x s) or dist object for phylogenetic distances
#   - "VMcor": logical, defining if the correction by Villeger & Mouillot (J Ecol, 2008) is applied or not
#   - "Jost": logical, defining if the correction by Jost (Ecology, 2007) is applied 
#   ...	Other arguments passed to functions, e.g. method, thin or burnin arguments for oecosimu function (vegan package)
#
# DETAILS: 
# Additive partitioning is performed on taxonomic diversity (Simpson's index), functional and phylogenetic diversity (Rao's index).
# The results are compared to expected values from nsimul permutations by individual based randomisation of the community data matrix.
# This is done by the "r2dtable" method in oecosimu by default. 
#
#


TFPdivpart <- function(y, x, index=c("richness", "Simpson", "RaoFD", "RaoPD"),
           weights=c("unif", "prop"), relative = FALSE, nsimul=99,
           method = "r2dtable", fundiv = NULL, phyldiv = NULL, 
           VMcor = FALSE, Jost = FALSE, ...)
  { 
    ## evaluate formula
    lhs <- as.matrix(y)
    if (missing(x))
      x <- cbind(level_1=seq_len(nrow(lhs)),
                 leve_2=rep(1, nrow(lhs)))
    rhs <- data.frame(x)
    rhs[] <- lapply(rhs, as.factor)
    rhs[] <- lapply(rhs, droplevels, exclude = NA)
    nlevs <- ncol(rhs)
    if (nlevs < 2)
      stop("provide at least two level hierarchy")
    if (any(rowSums(lhs) == 0))
      stop("data matrix contains empty rows")
    if (any(lhs < 0))
      stop("data matrix contains negative entries")
    if (is.null(colnames(rhs)))
      colnames(rhs) <- paste("level", 1:nlevs, sep="_")
    tlab <- colnames(rhs)
    
    ## check proper design of the model frame
    l1 <- sapply(rhs, function(z) length(unique(z)))
    if (!any(sapply(2:nlevs, function(z) l1[z] <= l1[z-1])))
      stop("number of levels are inappropriate, check sequence")
    rval <- list()
    rval[[1]] <- rhs[,nlevs]
    nCol <- nlevs - 1
    for (i in 2:nlevs) {
      rval[[i]] <- interaction(rhs[,nCol], rval[[(i-1)]], drop=TRUE)
      nCol <- nCol - 1
    }
    rval <- as.data.frame(rval[rev(seq_along(rval))])
    l2 <- sapply(rval, function(z) length(unique(z)))
    if (any(l1 != l2))
      stop("levels are not perfectly nested")
    
    ## aggregate response matrix
    fullgamma <-if (nlevels(rhs[,nlevs]) == 1)
      TRUE else FALSE
    ftmp <- vector("list", nlevs)
    for (i in seq_len(nlevs)) {
      ftmp[[i]] <- as.formula(paste("~", tlab[i], "- 1"))
    }
    
    ## is there burnin/thin in ... ?
    burnin <- if (is.null(list(...)$burnin))
      0 else list(...)$burnin
    thin <- if (is.null(list(...)$thin))
      1 else list(...)$thin
    base <- if (is.null(list(...)$base))
      exp(1) else list(...)$base
    
   
       
    ## evaluate other arguments
    index <- match.arg(index)  
    weights <- match.arg(weights)
    
    ## checking absence in 'fundiv' for RaoFD index
    if (is.null(fundiv) && index == "RaoFD")  stop("'fundiv': dist object is missing")
    ## checking absence in 'phyldiv' for RaoPD index
    if (is.null(phyldiv) && index == "RaoPD")  stop("'phyldiv': dist object is missing")
    ####function Rao by de Bello et al.(J Veg Sci, 2010)#####
    if(Jost){
      
      switch(index,
             "richness" = {divfun <- function(x) rowSums(x > 0)},
             
             "Simpson" = {divfun <- function(x) Rao(t(x), dfunc = NULL, dphyl = NULL, weight=VMcor, Jost=TRUE, structure=NULL)$TD$Alpha},             
             
             "RaoFD" = {divfun <- function(x) Rao(t(x), dfunc=fundiv, dphyl=NULL, weight=VMcor, Jost=TRUE, structure=NULL)$FD$Alpha},
             
             "RaoPD" = {divfun <- function(x) Rao(t(x), dfunc=NULL, dphyl=phyldiv, weight=VMcor, Jost=TRUE, structure=NULL)$PD$Alpha})
    } else {
      
      switch(index,
             "richness" = {divfun <- function(x) rowSums(x > 0)},
             
             "Simpson" = {divfun <- function(x) Rao(t(x), dfunc = NULL, dphyl = NULL, weight=VMcor, Jost=FALSE, structure=NULL)$TD$Alpha},
             
             "RaoFD" = {divfun <- function(x) Rao(t(x), dfunc=fundiv, dphyl=NULL, weight=VMcor, Jost=FALSE, structure=NULL)$FD$Alpha},
             
             "RaoPD" = {divfun <- function(x) Rao(t(x), dfunc=NULL, dphyl=phyldiv, weight=VMcor, Jost=FALSE, structure=NULL)$PD$Alpha})
      
    }
    
    ## this is the function passed to oecosimu
    wdivfun <- function(x) {
      ## matrix sum *can* change in oecosimu (but default is constant sumMatr)
      sumMatr <- sum(x)
      if (fullgamma) {
        tmp <- lapply(seq_len(nlevs-1), function(i) t(model.matrix(ftmp[[i]], rhs)) %*% x)
        tmp[[nlevs]] <- matrix(colSums(x), nrow = 1, ncol = ncol(x))
      } else {
        tmp <- lapply(seq_len(nlevs), function(i) t(model.matrix(ftmp[[i]], rhs)) %*% x)
      }
      ## weights will change in oecosimu thus need to be recalculated
      if (weights == "prop")
        wt <- lapply(seq_len(nlevs), function(i) apply(tmp[[i]], 1, function(z) sum(z) / sumMatr))
      else wt <- lapply(seq_len(nlevs), function(i) rep(1 / NROW(tmp[[i]]), NROW(tmp[[i]])))
      a <- sapply(seq_len(nlevs), function(i) sum(divfun(tmp[[i]]) * wt[[i]]))
      if (relative)
        a <- a / a[length(a)]
      b <- sapply(2:nlevs, function(i) a[i] - a[(i-1)])
      c(a, b)
    }
    if (nsimul > 0) {
      sim <- oecosimu(lhs, wdivfun, method = method, nsimul=nsimul,
                      burnin=burnin, thin=thin)
    } else {
      sim <- wdivfun(lhs)
      tmp <- rep(NA, length(sim))
      sim <- list(statistic = sim,
                  oecosimu = list(z = tmp, pval = tmp, method = NA, statistic = sim))
    }
    nam <- c(paste("alpha", seq_len(nlevs-1), sep="."), "gamma",
             paste("beta", seq_len(nlevs-1), sep="."))
    names(sim$statistic) <- attr(sim$oecosimu$statistic, "names") <- nam
    call <- match.call()
    call[[1]] <- as.name("TFPdivpart")
    attr(sim, "call") <- call
    attr(sim$oecosimu$simulated, "index") <- index
    attr(sim$oecosimu$simulated, "weights") <- weights
    attr(sim, "n.levels") <- nlevs
    attr(sim, "terms") <- tlab
    attr(sim, "model") <- rhs
    class(sim) <- c("TFPdivpart", class(sim))
    sim
  }
  
  
###################################################################################################################################
#--------------------------------------------------EXAMPLES-----------------------------------------------------------------------#
###################################################################################################################################  
#generate a random trait matrix
trait_mat = matrix(c(0.60, 0.05, 0.20, 0, 0.25, 0.75, 0,0, 0.5, 0, 0.25, 0.5, 0.25, 0.25, 0.45, 0.25) , nrow = 4, ncol = 4, T)  
row.names(trait_mat) = c("sp1","sp2","sp3","sp4")
#find the dissimilarity or dist objest in functional traits
diss_trait_mat = dist(trait_mat)  
#generate a random table of abundances
abund = data.frame("sp1" = c(0,3,0,4), "sp2" = c(1,2,0,0), "sp3" = c(15,0,2,9), "sp4" = c(1,5,2,10))
row.names(abund) = c("S1","S2","S3","S4")
#calculate additive partitioning of functional diversity
TFPdivpart(abund, index="RaoFD", nsimul=19, fundiv = diss_trait_mat, Jost = F)


#with levels of sampling hierarchy
hier = with(abund, data.frame("l1" = 1:nrow(abund),
                "l2" = c(1,2,3,1),
                "l3" = c(1,2,2,1),
                "l4" = c(1,1,1,1)))               
TFPdivpart(abund, hier, index="Simpson", nsimul=19)                
TFPdivpart(abund, hier, index="RaoFD", nsimul=19, fundiv = diss_trait_mat, Jost = F)
