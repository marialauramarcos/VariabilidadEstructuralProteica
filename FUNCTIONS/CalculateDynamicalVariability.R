CalculateDynamicalVariability <- function(r.p.1,
                                          r.p.2,
                                          aligned.p.1.index,
                                          aligned.p.2.index,
                                          not.aligned.p.1.index,
                                          not.aligned.p.2.index,
                                          R0, 
                                          tolerance, 
                                          K.analysis) {

  # analyse p.1 and p.2
  ## calculate ENMK of p.1
  ENMK.p.1 = CalculateENMKeff(r.p.1, 
                              aligned.p.1.index, 
                              not.aligned.p.1.index, 
                              R0, 
                              tolerance, 
                              K.analysis)

  ## calculate ENMK of p.2
  ENMK.p.2 = CalculateENMKeff(r.p.2, 
                              aligned.p.2.index, 
                              not.aligned.p.2.index, 
                              R0, 
                              tolerance, 
                              K.analysis)

  ## get cov matrices
  cov.p.1 = ENMK.p.1$Cov
  cov.p.2 = ENMK.p.2$Cov

  ## get eigen vectors
  ve.p.1 = ENMK.p.1$ve
  ve.p.2 = ENMK.p.2$ve

  # calculate Mean Square Fluctuation
  ## get the diagonal of the cov matrix
  diag.p.1 = diag(cov.p.1)
  diag.p.2 = diag(cov.p.2)

  ## calculate the factor to divide the diagonal
  factor = sort(rep(seq(1:n.aligned), 3))

  ## separate the diagonal
  s.diag.p.1 = split(diag.p.1, factor)
  s.diag.p.2 = split(diag.p.2, factor)

  ## create matrices to save data
  MSF.p.1 = matrix(rnow = 1, ncol = n.aligned)
  MSF.p.2 = matrix(rnow = 1, ncol = n.aligned)

  ## start a loop to calculate sums for each site
  for (i in (1:n.aligned)) {
    MSF.p.1[, i] = sum(unlist(s.diag.p.1[i]), use.names = F)
    MSF.p.2[, i] = sum(unlist(s.diag.p.2[i]), use.names = F)
  }
 
  ## calculate the difference
  square.dif.MSF = (MSF.p.1 - MSF.p.2) ^ 2
  
  # calculate nH and nR
  overlap = t(ve.p.1) %*% ve.p.2
  nH = exp( - rowSums(overlap ^ 2 * log(overlap ^ 2 + tolerance)))
  nR = 1 / rowSums(overlap ^ 4)
  
  # create a list for the output
  out = list("MSF.p.1" = MSF.p.1,
             "MSF.p.2" = MSF.p.2,
      "square.dif.MSF" = square.dif.MSF, 
                  "nH" = nH,
                  "nR" = nR) 
  out
}
