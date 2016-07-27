# Description:
#
# This function calculates aligned and not aligned sites of two proteins p.ref and p.2 in a 
# multiple alignment.
#
#  Args:
#   - alignment: multiple alignment.
#   - pdbid: pdbidS of proteins in the rows of the alignment.
#   - p.ref: pdbid of the reference protein.
#   - p.2: pdbid of the other protein to analyze.
#
#  Returns:
#    n.sites.p.ref: n sites of p.ref.
#    n.sites.p.2: n sites of p.2.
#    n.aligned: n aligned sites of p.ref.
#    aligned.p.ref.index: aligned sites of p.ref.
#    aligned.p.2.index: aligned sites of p.2.
#    not.aligned.p.ref.index: not aligned sites of p.ref.
#    not.aligned.p.2.index: not aligned sites of p.2.
#    n.aligned.mut.p.ref: n aligned but mutated sites of p.ref.
#    aligned.mut.p.ref.index: aligned but mutated sites of p.ref.

AnalyzeAlignment <- function(alignment, pdbid, p.ref, p.2) {
  
  # Extract data form the alignment
  l.alignment = ncol(alignment)
  n.prot = nrow(alignment)
  alignment.p.ref = alignment[which(grepl(p.ref, pdbid)), ]
  alignment.p.2 = alignment[which(grepl(p.2, pdbid)), ]

  # Get aa indexes in the alignment
  aa.p.ref = 0
  aa.p.2 = 0
  p.ref.index = c()
  p.2.index = c()
  
  for (i in (1:l.alignment)) {
    if (alignment.p.ref[i] != "-") {
      aa.p.ref = aa.p.ref + 1
      p.ref.index[i] = aa.p.ref
    }
    if (alignment.p.2[i] != "-") {
      aa.p.2 = aa.p.2 + 1
      p.2.index[i] = aa.p.2
    }
  }	
  
  # Get aligned and not.aligned indexes and aa of p.ref and p.2
  aligned.p.ref.index = c()
  aligned.p.2.index = c()
  not.aligned.p.ref.index = c()
  not.aligned.p.2.index = c()
  
  aligned.p.ref = c()
  aligned.p.2 = c()
  
  for (i in (1:l.alignment)) {
    if (alignment.p.ref[i] != "-" & alignment.p.2[i] != "-") {
      aligned.p.ref.index = cbind(aligned.p.ref.index, p.ref.index[i])
      aligned.p.2.index = cbind(aligned.p.2.index, p.2.index[i])
      
      aligned.p.ref = cbind(aligned.p.ref, alignment.p.ref[i])
      aligned.p.2 = cbind(aligned.p.2, alignment.p.2[i])
    }
    if (alignment.p.ref[i] != "-" & alignment.p.2[i] == "-") {
      not.aligned.p.ref.index = cbind(not.aligned.p.ref.index, p.ref.index[i])
    }
    if (alignment.p.2[i] != "-" & alignment.p.ref[i] == "-") {
      not.aligned.p.2.index = cbind(not.aligned.p.2.index, p.2.index[i])
    }
  }	
  n.aligned = length(aligned.p.ref.index)
  
  # Get aligned but mutated indexes of p.ref
  aligned.mut.p.ref.index = c()

  for (i in (1:n.aligned)) {
    if (aligned.p.ref[i] != aligned.p.2[i]) {
      aligned.mut.p.ref.index = cbind(aligned.mut.p.ref.index, aligned.p.ref.index[i])
    }
  }
  
  n.aligned.mut.p.ref = length(aligned.mut.p.ref.index)
  
  # Get the % sequence identity
  identity = 100 - (n.aligned.mut.p.ref * 100 / n.aligned)
  
  # Create a list for the output.
  output = list(    "n.sites.p.ref" = aa.p.ref,  
                      "n.sites.p.2" = aa.p.2, 
                        "n.aligned" = n.aligned,
              "aligned.p.ref.index" = aligned.p.ref.index,
                "aligned.p.2.index" = aligned.p.2.index,
          "not.aligned.p.ref.index" = not.aligned.p.ref.index,
            "not.aligned.p.2.index" = not.aligned.p.2.index,
              "n.aligned.mut.p.ref" = n.aligned.mut.p.ref,
          "aligned.mut.p.ref.index" = aligned.mut.p.ref.index,
                         "identity" = identity)
  output
}
