# This function calculates K effective of r, its eigenvalues and eigenvectors.
#
#  Args:
#    r: 3 x nsites matrix containing the equilibrium coordinates of each site.
#    aligned.index: aligned sites
#    not.aligned.index: not aligned sites
#    kijFunction: name of the function that calculates kij between 2 sites i and j.
#    R0: cut-off for ANM.
#    TOLERANCE: 0 tolerance.
#
#  Returns:
#    KEFF: 3*naligned x 3*naligned KEFF matrix.
#    va: eigenvalues > TOLERANCE of KEFF.
#    ve: eigenvectors > TOLERANCE of KEFF.
CalculateKeff <- function(r, aligned.index, not.aligned.index, CalculateKij, R0, TOLERANCE) {
  naligned = length(aligned.index)
  n.not.aligned = length(not.aligned.index)
  nsites = ncol(r)
  # Calculate K.
  K <- CalculateK(r, CalculateKij, R0, TOLERANCE)$K 
  # Order K in order to get KPP, KQQ, KPQ and KQP and KEFF.
  sortm <- matrix(0, ncol = 3 * nsites, nrow = 1)
  count = 1
  for (i in (1:nsites)) {
    for (j in (1:naligned)) {
      if (aligned.index[j] == i) {
        sortm[1, (3 * i - 2)] <- 3 * count - 2
        sortm[1, (3 * i - 1)] <- 3 * count - 1
        sortm[1, (3 * i)] <- 3 * count
        count <- count + 1
      }
    }
  }
  if (n.not.aligned > 0) {
    for (i in (1:nsites)) {
      for (j in (1:n.not.aligned)) {
        if (not.aligned.index[j] == i) {
          sortm[1, (3 * i - 2)] <- 3 * count - 2
          sortm[1, (3 * i - 1)] <- 3 * count - 1
          sortm[1, (3 * i)] <- 3 * count
          count <- count + 1
        }
      }
    }
  }
  Kord <- K[sortm, sortm]
  KPP <- Kord[(1:(3 * naligned)), (1:(3 * naligned))]
  if(n.not.aligned > 0) {
    KQQ <- Kord[((3 * naligned + 1):(3 * nsites)), ((3 * naligned + 1):(3 * nsites))]
    KPQ <- Kord[(1:(3 * naligned)), ((3 * naligned + 1):(3 * nsites))]
    KQP <- t(KPQ)
    # Calculate KQQ^-1.
    eigQQ <- eigen(KQQ, symmetric = T)
    veQQ <- eigQQ$vectors
    vaQQ <- eigQQ$values
    modes <- vaQQ > TOLERANCE #there aren´t negative evalues#
    veQQ  <- veQQ[, modes]
    vaQQ <- vaQQ[modes]
    covQQ <-  veQQ %*% ((1 / vaQQ) * t(veQQ))
    # Calculate KEFF.
    KEFF <- KPP - (KPQ %*% covQQ %*% KQP)
  }
  if(n.not.aligned == 0) {
    KEFF <- KPP
  }
  # Calculate and order eigenvalues and eigenvectors of KEFF.
  eigKEFF <- eigen(KEFF, symmetric = T)
  veKEFF <- eigKEFF$vectors
  vaKEFF <- eigKEFF$values
  ord <- sort.list(Mod(vaKEFF), decreasing = FALSE)
  vaKEFF <- vaKEFF[ord]
  veKEFF <- veKEFF[, ord]
  modes <- vaKEFF > TOLERANCE  # there aren´t negative evalues.
  veKEFF <- veKEFF[, modes]
  vaKEFF <- vaKEFF[modes]
  output = list("KEFF" = KEFF, "ve" = veKEFF, "va" = vaKEFF)
  output
}


