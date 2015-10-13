#This function calculates measures of structural and dynamical variability between 2 K/KEFF objects.
#
#  Args:
#    dr: Diference between equilibrium coordinates of 2 proteins.
#    K.1: an object of P1 obtained using CalculateK() or CalculateKEFF().
#    K.2: an object of P2 obtained using CalculateK() or CalculateKEFF().
#
#  Returns:
#    dr.squarei
#    nH
#    nR
#    d.evalues
#    Pn
#    Qn
CalculateVariability <- function(dr, K.1, K.2) {
  dr.squarei = (colSums(dr ^ 2))
  ov = t(K.1$ve) %*% K.2$ve
  nH = exp( - rowSums(ov ^ 2 * log(ov ^ 2 + TOLERANCE)))
  nR = 1 / rowSums(ov ^ 4)
  d.evalues = K.2$va - K.1$va
  Pn = (t(K.1$ve) %*% matrix(dr[1:3, ], ncol = 1)) ^ 2 / sum((t(K.1$ve) %*% matrix(dr[1:3, ] , ncol = 1)) ^ 2)
  Qn = cumsum(Pn)
  output = list(dr.squarei = dr.squarei, nH = nH, nR = nR, d.evalues = d.evalues, Pn = Pn, Qn = Qn)
  output
}