# This function calculates measures of variability between two structures. 
#
#  Args:
#    r.p.1: vector with coordinates of p.1.
#    r.p.2: vector with coordinates of p.2.
#    aligned.p.1.index: aligned sites of p.1.
#    aligned.p.2.index: aligned sites of p.2.
#    not.aligned.p.1.index: not.aligned sites of p.1.
#    not.aligned.p.2.index: not.aligned sites of p.2.
#    R0: cut-off for the ANM.
#    rotate: it can be "TRUE" or "FALSE". If it is "TRUE" r.p.2 is rotaded in order to minimize RMSD with r.p.1.
#    TOLERANCE: 0 tolerance.
# 
#  Requires functions:
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
                                 TOLERANCE) {
  
  # Calculate 3N indexes.
  aligned.p.1.index.3N = sort(c(aligned.p.1.index * 3, aligned.p.1.index * 3 - 2, aligned.p.1.index * 3 - 1))
  aligned.p.2.index.3N = sort(c(aligned.p.2.index * 3, aligned.p.2.index * 3 - 2, aligned.p.2.index * 3 - 1))
  
  # Rotate r.p.2 in order to minimize MSD with r.p.1.
  if (rotate == TRUE) {
    r.p.2 = fit.xyz(fixed = r.p.1,
                   mobile = r.p.2,
               fixed.inds = aligned.p.1.index.3N, 
              mobile.inds = aligned.p.2.index.3N) 
  }

  # Calculate dr.
  dr = r.p.2[aligned.p.2.index.3N] - r.p.1[aligned.p.1.index.3N]

  # Calculate ENM Keff of p.1. 
  ENMKeff.p.1 <- CalculateENMKeff(matrix(r.p.1, nrow = 3),
                                  aligned.p.1.index, 
                                  not.aligned.p.1.index,
                                  R0,
                                  TOLERANCE)
  
  # Get evalues and eigenvectors of p.1.
  va = ENMKeff.p.1$va
  ve = ENMKeff.p.1$ve
  
  # Calculate measures of variability between p.1 and p.2.
  Pn = ((t(ve) %*% dr) ^ 2)/ (sum(dr ^ 2))
  dr.squarei = colSums(matrix(dr, nrow = 3) ^ 2) 

  # Create a list for the output.
  output = list("va" = va,
                "Pn" = Pn, 
        "dr.squarei" = dr.squarei)
  
  output
}

