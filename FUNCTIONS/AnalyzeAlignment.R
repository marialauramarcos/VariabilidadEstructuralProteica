# This function takes the alignment of 2 proteins and returs aligned and not aligned sites of P1.
#
# Args:
#   alignment.P1: 1 x lalignment matrix containing the alignment of P1 extracted from the multiple
#   sequence alignment.
#   alignment.P2: 1 x lalignment matrix containing the alignment of P2 extracted from the multiple
#   sequence alignment.
#   nsites.P1: Number of aminoacids of P1.
#
# Returns:
#   aligned.index: Sites of P1 that aligne with P2.
#   not.aligned.index: Sites of P1 that do not aligne with P2.
#   naligned: Number of aminoacides of P1 that aligne with P2.
AnalyzeAlignment <- function(alignment.P1, alignment.P2, nsites.P1) {
  lalignment = length(alignment.P1)
  index.P1 = seq(1, nsites.P1)
  alignment.P1.index = matrix(ncol = lalignment, nrow = 2)
  alignment.P1.index[1, ] = alignment.P1
  count = 1
  for (i in (1:lalignment)) {
    if (alignment.P1[i] != "-") {
      alignment.P1.index[2, i] <- index.P1[count]
      count = count + 1
    }
  }	
  aligned.P1.index = matrix(ncol = 1, nrow = 1)
  not.aligned.P1.index = matrix(ncol = 1, nrow = 1)
  for (i in (1:lalignment)) {
    if (alignment.P2[i] != "-" & alignment.P1[i] != "-") {
      aligned.P1.index <- cbind(aligned.P1.index, alignment.P1.index[2, i])
    }
    if (alignment.P2[i] == "-" & alignment.P1[i] != "-") {
      not.aligned.P1.index <- cbind(not.aligned.P1.index, alignment.P1.index[2, i])
    }
  }	
  aligned.P1.index = matrix(as.numeric(aligned.P1.index[-1]), nrow = 1)
  not.aligned.P1.index = matrix(as.numeric(not.aligned.P1.index[-1]), nrow = 1)	
  naligned <- ncol(aligned.P1.index)
  output = list("aligned.index" = aligned.P1.index, "not.aligned.index" = not.aligned.P1.index, "naligned" = naligned)
  output
}


