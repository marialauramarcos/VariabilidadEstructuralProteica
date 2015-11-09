# This function calculates % sequence identity between proteins P1 and P2.
#
# Args:
#   df.alignment: data frame with data of the alignment. 
#   
# Returns:
#   ID: % sequence identity between P1 and P2.
CalculateID <- function(df.alignment) {
  # Prepare data.
  alignment <- df.alignment$alignment
  p.1 <- df.alignment$p.1
  p.2 <- df.alignment$p.2
  
  # Extract data form the alignment.
  lalignment = ncol(alignment)
  alignment.p.1 = alignment[id == p.1,]
  alignment.p.2 = alignment[id == p.2,]
  
  aligned.p.1 = c()
  aligned.p.2 = c()
  for (i in (1:lalignment)) {
    if (alignment.p.1[i] != "-" & alignment.p.2[i] != "-") {
      aligned.p.1 = cbind(aligned.p.1, alignment.p.1[i])
      aligned.p.2 = cbind(aligned.p.2, alignment.p.2[i])
    }
  }
  naligned = length(aligned.p.1)	
  count = 0
  for (i in (1:naligned)) {
    if (aligned.p.1[i] == aligned.p.2[i]) {
      count = count +1
    }
  }
  ID = count * 100 / naligned
  ID
}

