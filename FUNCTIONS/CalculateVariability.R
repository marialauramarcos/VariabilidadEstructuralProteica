# This function calculates measures of variability between two structures. 
#
#  Args:
#    r.p.1: vector of coordinates of p.1.
#    r.p.2: vector of coordinates of p.2.
#    aligned.p.1.index: aligned sites of p.1.
#    aligned.p.2.index: aligned sites of p.2.
#    not.aligned.p.1.index: not aligned sites of p.1.
#    not.aligned.p.2.index: not aligned sites of p.2.
#    R0: cut-off for the ANM.
#    rotate: it can be "TRUE" or "FALSE". If it is "TRUE", r.p.2 is rotated in order to minimize RMSD with r.p.1.
#    K.analysis: It can be "K" or "Keff". For "K" or "Keff", the analysis is based on normal modes of "K" or "Keff"
#    respectibly. 
#    TOLERANCE: 0 tolerance.
# 
#  Required functions:
#    CalculateENMKeff
#    CalculateENMK
#    CalculateKij
#
#  Returns:
#    Pn
#    dr.squarei
#    va

CalculateVariability <- function(r.p.1, 
                                 r.p.2, 
                                 aligned.p.1.index, 
                                 aligned.p.2.index, 
                                 not.aligned.p.1.index,
                                 not.aligned.p.2.index,
                                 R0,
                                 rotate,
                                 K.analysis,
                                 TOLERANCE) {
  
  # Calculate 3N indexes.
  aligned.p.1.index.3N = sort(c(aligned.p.1.index * 3, aligned.p.1.index * 3 - 2, aligned.p.1.index * 3 - 1))
  aligned.p.2.index.3N = sort(c(aligned.p.2.index * 3, aligned.p.2.index * 3 - 2, aligned.p.2.index * 3 - 1))
  
  # Rotate r.p.2 in order to minimize RMSD with r.p.1.
  if (rotate == TRUE) {
    r.p.2 = fit.xyz(fixed = r.p.1,
                   mobile = r.p.2,
               fixed.inds = aligned.p.1.index.3N, 
              mobile.inds = aligned.p.2.index.3N) 
  }

  # Calculate dr.
  dr = r.p.2[aligned.p.2.index.3N] - r.p.1[aligned.p.1.index.3N]

  # Calculate ENM Keff of p.1.
  ENMKeff.p.1 = CalculateENMKeff(matrix(r.p.1, nrow = 3),
                                 aligned.p.1.index, 
                                 not.aligned.p.1.index,
                                 R0,
                                 TOLERANCE,
                                 K.analysis)
  
  # Get eigenvalues and eigenvectors of p.1.
  va = ENMKeff.p.1$va
  ve = ENMKeff.p.1$ve
  
  # correct ve for K.analysis == "K".
  if (K.analysis == "K") {
    ve = ve[aligned.p.1.index.3N,]
  }
  
  # Calculate measures of variability between p.1 and p.2.
  Pn = ((t(ve) %*% dr) ^ 2) / (sum((t(ve) %*% dr) ^ 2))
  dr.squarei = colSums(matrix(dr, nrow = 3) ^ 2) 

  # Create a list for the output.
  output = list( "va" = va,
                 "ve" = ve,
                 "Pn" = Pn, 
         "dr.squarei" = dr.squarei)
  
  output
}

