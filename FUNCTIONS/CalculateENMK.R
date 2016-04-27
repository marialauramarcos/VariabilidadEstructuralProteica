# This function calculates K of r, its Covariance, eigenvalues and eigenvectors.
#
#  Args:
#    r: 3 x nsites matrix containing the equilibrium coordinates of each site.
#    CalculateKij: function that calculates kij between 2 sites i and j.
#    R0: cut-off for ANM.
#    TOLERANCE: 0 tolerance.
#
#  Returns;
#    K: 3*nsites x 3*nsites K matrix of r.
#    kij: nsites x nsites kij matrix of r.
#    cov: 3*nsites x 3*nsites covariance matrix of K.
#    va: eigenvalues > TOLERANCE of K.
#    ve: eigenvectors > TOLERANCE of K.
CalculateENMK <- function(r, CalculateKij, R0, TOLERANCE) {
  # Calculate K.
  nsites <- ncol(r)
  kij <- matrix(0, ncol = nsites, nrow = nsites)
  K <- matrix(0, ncol = 3*nsites, nrow = 3*nsites)
  for (i in (1:(nsites - 1))) {
    ai <- (i + (2 * (i - 1)))
    bi <- ai + 2
    for (j in ((i + 1) : nsites)) {
      aj <- (j + (2 * (j - 1)))				
      bj <- aj + 2	
      rij <- r[1:3, j] - r[1:3, i]
      dij <- drop(sqrt(rij %*% rij))
      eij <- rij / dij
      kij[i, j] <- CalculateKij(dij, R0)
      kij[j, i] <- kij[i, j]					
      K[ai:bi, aj:bj] <- -kij[i, j] * (eij %*% t(eij))
      K[aj:bj, ai:bi] <- K[ai:bi, aj:bj]
      K[ai:bi, ai:bi] <- K[ai:bi, ai:bi] - K[ai:bi, aj:bj]
      K[aj:bj, aj:bj] <- K[aj:bj, aj:bj] - K[aj:bj, ai:bi]
    }	
  }
  # Calculate and order eigenvalues and eigenvectors of K.
  eig <- eigen(K, symmetric = TRUE)
  ord <- sort.list(Mod(eig$values), decreasing = FALSE)
  va <- eig$values[ord]
  ve <- eig$vectors[, ord]
  modes <- va > TOLERANCE
  va <- va[modes]
  ve  <- ve[, modes]
  nmodes <- length(va)
  # Calculate covariance matrix of K.
  cov <-  ve %*% ((1 / va) * t(ve))

  output <- list("K" = K, "kij" = kij, "cov" = cov, "va" = va, "ve" = ve, "nmodes" = nmodes)
  output
}
