# This function calculates % sequence identity between proteins P1 and P2.
#
# Args:
#   alignment.P1: 1 x lalignment matrix containing the alignment of P1 extracted from the multiple
#   sequence alignment.
#   alignment.P2: 1 x lalignment matrix containing the alignment of P2 extracted from the multiple
#   sequence alignment.
#   
# Returns:
#   ID: % sequence identity between P1 and P2.
CalculateID <- function(alignment.P1, alignment.P2) {
  lalignment = length(alignment.P1)
  aligned.P1 = matrix(0)
  aligned.P2 = matrix(0)
  for (i in (1:lalignment)) {
    if (alignment.P1[i] != "-" & alignment.P2[i] != "-") {
      aligned.P1 = cbind(aligned.P1,alignment.P1[i])
      aligned.P2 = cbind(aligned.P2,alignment.P2[i])
    }
  }
  aligned.P2 = aligned.P2[, -1]	
  aligned.P1 = aligned.P1[, -1]	
  naligned = length(aligned.P1)	
  count = 0
  for (i in (1:naligned)) {
    if (aligned.P1[i] == aligned.P2[i]) {
      count = count +1
    }
  }
  ID = count * 100 / naligned
  ID
}

