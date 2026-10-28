function (datExpr, weights = NULL, corType = "pearson", networkType = "unsigned", 
    power = 6, TOMType = "signed", TOMDenom = "min", maxPOutliers = 1, 
    quickCor = 0, pearsonFallback = "individual", cosineCorrelation = FALSE, 
    replaceMissingAdjacencies = FALSE, suppressTOMForZeroAdjacencies = FALSE, 
    suppressNegativeTOM = FALSE, useInternalMatrixAlgebra = FALSE, 
    nThreads = 0, verbose = 1, indent = 0) 
{
    corTypeC = as.integer(pmatch(corType, .corTypes) - 1)
    if (is.na(corTypeC)) 
        stop(paste("Invalid 'corType'. Recognized values are", 
            paste(.corTypes, collapse = ", ")))
    TOMTypeC = as.integer(pmatch(TOMType, .TOMTypes) - 1)
    if (is.na(TOMTypeC)) 
        stop(paste("Invalid 'TOMType'. Recognized values are", 
            paste(.TOMTypes, collapse = ", ")))
    TOMDenomC = as.integer(pmatch(TOMDenom, .TOMDenoms) - 1)
    if (is.na(TOMDenomC)) 
        stop(paste("Invalid 'TOMDenom'. Recognized values are", 
            paste(.TOMDenoms, collapse = ", ")))
    if ((maxPOutliers < 0) | (maxPOutliers > 1)) 
        stop("maxPOutliers must be between 0 and 1.")
    if (quickCor < 0) 
        stop("quickCor must be positive.")
    if ((maxPOutliers < 0) | (maxPOutliers > 1)) 
        stop("maxPOutliers must be between 0 and 1.")
    fallback = as.integer(pmatch(pearsonFallback, .pearsonFallbacks))
    if (is.na(fallback)) 
        stop(paste("Unrecognized 'pearsonFallback'. Recognized values are (unique abbreviations of)\n", 
            paste(.pearsonFallbacks, collapse = ", ")))
    if (nThreads < 0) 
        stop("nThreads must be positive.")
    if (is.null(nThreads) || (nThreads == 0)) 
        nThreads = .useNThreads()
    if ((power < 1) | (power > 30)) 
        stop("power must be between 1 and 30.")
    networkTypeC = as.integer(charmatch(networkType, .networkTypes) - 
        1)
    if (is.na(networkTypeC)) 
        stop(paste("Unrecognized networkType argument.", "Recognized values are (unique abbreviations of)", 
            paste(.networkTypes, collapse = ", ")))
    dimEx = dim(datExpr)
    if (length(dimEx) != 2) 
        stop("datExpr has incorrect dimensions.")
    if (length(weights) > 0) {
        if (!isTRUE(all.equal(dimEx, dim(weights)))) 
            stop("When 'weights' are given, they must have the same dimensions as 'datExpr'.")
        if (any(!is.finite(weights))) {
            if (verbose > 0) 
                warning("Found non-finite weights. The corresponding data points will be removed.")
            weights[!is.finite(weights)] = NA
        }
    }
    nGenes = dimEx[2]
    nSamples = dimEx[1]
    warn = as.integer(0)
    datExpr = as.matrix(datExpr)
    if (length(weights) > 0) 
        weights = as.matrix(weights)
    else weights = NULL
    tom = .Call("tomSimilarity_call", datExpr, weights, as.integer(corTypeC), 
        as.integer(networkTypeC), as.double(power), as.integer(TOMTypeC), 
        as.integer(TOMDenomC), as.double(maxPOutliers), as.double(quickCor), 
        as.integer(fallback), as.integer(cosineCorrelation), 
        as.integer(replaceMissingAdjacencies), as.integer(suppressTOMForZeroAdjacencies), 
        as.integer(suppressNegativeTOM), as.integer(useInternalMatrixAlgebra), 
        warn, as.integer(nThreads), as.integer(verbose), as.integer(indent), 
        PACKAGE = "WGCNA")
    diag(tom) = 1
    return(tom)
}
