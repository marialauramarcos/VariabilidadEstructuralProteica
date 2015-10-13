# This function calculates the force that models a mutation using LF-ENM.
#
#  Args:
#    l: Mutated site.
#    r: 3 x nsites matrix containing the equilibrium coordinates of each site.
#    kij: nsites x nsites matrix containing kij between all sites i and j.
#    fmax: maximum value for f.
#
#  Returns:
#    fvec: 3 * nsites vector containing the force that models the mutation of site l.
CalculateForce = function(l, r, kij, fmax) {
  nsites = ncol(r)
  fvec = matrix(0, 3, nsites)
  for (j in seq(nsites)){
    if(kij[l, j] > 0.9) {
      rij = r[1:3, j] - r[1:3, l]
      eij = rij/sqrt(sum(rij ^ 2))
      fij = fmax * sqrt(3) * runif(1, -1, 1)
      fvec[, j] = fij * eij
      fvec[, l] = fvec[, l] - fvec[, j]
    }
  }
  dim(fvec) = c(3 * nsites)
  fvec
}



