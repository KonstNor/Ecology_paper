TFPdivpart.formula <-
    function(formula, data, index=c("richness", "simpson", "Rao"),
             weights=c("unif", "prop"), relative = FALSE, nsimul=99,
             method = "r2dtable", fundiv = NULL, Jost = FALSE, ...)
{
    ## evaluate formula
    if (missing(data))
        data <- parent.frame()
    tmp <- hierParseFormula(formula, data)
    ## run simulations
    sim <- TFPdivpart.default(tmp$lhs, tmp$rhs, index=index, weights = weights,
                           relative = relative, nsimul = nsimul,
                           method = method, fundiv = NULL, Jost = FALSE, ...)
    
    call <- match.call()
    call[[1]] <- as.name("TFPdivpart")
    attr(sim, "call") <- call
    sim
}