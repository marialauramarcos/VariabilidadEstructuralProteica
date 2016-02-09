# This function calculates aligned and not aligned sites of two proteins p.1 and p.2 in a 
# multiple alignment. It has to possibilitys:
#  - core: aligned and not aligned sites of p.1 and p.2 that are in the 
# conserved core of the alignment: positions with no gaps in the alignemnt.
#  - noCore: all aligned and not aligned sites of p.1 and p.2.
#
#  Args:
#      alignment: multiple alignment.
#      id: pdbid of proteins in the rows of the alignment.
#      p.1: pdbid of one of the proteins to analyze.
#      p.2: pdbid of the other protein to analyze.
#
#  Returns:
#    n.sites.p.1: n sites of p.1.
#    n.sites.p.2: n sites of p.2.
#    n.aligned: n aligned sites of p.1.
#    aligned.p.1.index: aligned sites of p.1.
#    aligned.p.2.index: aligned sites of p.2.
#    not.aligned.p.1.index: not.aligned sites of p.1.
#    not.aligned.p.2.index: not.aligned sites of p.2.
#    n.aligned.mut.p.1: n aligned but mutated sites of p.1.
#    aligned.mut.p.1.index: aligned but mutated sites of p.1.
#    n.core: n core sites of p.1.
#    core.p.1.index: core sites of p.1.
#    core.p.2.index: core sites of p.2.
#    no.core.p.1.index: no.core sites of p.1.
#    no.core.p.2.index: no.core sites of p.2.
#    identity: % sequence identity.

AnalyzeAlignment <- function(alignment, pdbid, p.1, p.2) {
  
  # Extract data form the alignment.
  l.alignment = ncol(alignment)
  n.prot = nrow(alignment)
  alignment.p.1 = alignment[which(grepl(p.1, pdbid)), ]
  alignment.p.2 = alignment[which(grepl(p.2, pdbid)), ]

  # Get aa indexes in the alignment.
  p.1.index = c()
  p.2.index = c()
  aa.p.1 = 0
  aa.p.2 = 0
  
  for (i in (1:l.alignment)) {
    if (alignment.p.1[i] != "-") {
      aa.p.1 = aa.p.1 + 1
      p.1.index[i] = aa.p.1
    }
    if (alignment.p.2[i] != "-") {
      aa.p.2 = aa.p.2 + 1
      p.2.index[i] = aa.p.2
    }
  }	
  
  ### No core ###
  
  # Get aligned and not.aligned indexes and aa of p.1 and p.2.
  aligned.p.1.index = c()
  aligned.p.2.index = c()
  not.aligned.p.1.index = c()
  not.aligned.p.2.index = c()
  
  aligned.p.1 = c()
  aligned.p.2 = c()
  
  for (i in (1:l.alignment)) {
    if (alignment.p.1[i] != "-" & alignment.p.2[i] != "-") {
      aligned.p.1.index = cbind(aligned.p.1.index, p.1.index[i])
      aligned.p.2.index = cbind(aligned.p.2.index, p.2.index[i])
      
      aligned.p.1 = cbind(aligned.p.1, alignment.p.1[i])
      aligned.p.2 = cbind(aligned.p.2, alignment.p.2[i])
    }
    if (alignment.p.1[i] != "-" & alignment.p.2[i] == "-") {
      not.aligned.p.1.index = cbind(not.aligned.p.1.index, p.1.index[i])
    }
    if (alignment.p.2[i] != "-" & alignment.p.1[i] == "-") {
      not.aligned.p.2.index = cbind(not.aligned.p.2.index, p.2.index[i])
    }
  }	
  n.aligned = length(aligned.p.1.index)
  
  # Get % sequence identity and aligned but mutated sites of p.1.
  aligned.mut.p.1.index = c()

  for (i in (1:n.aligned)) {
    if (aligned.p.1[i] != aligned.p.2[i]) {
      aligned.mut.p.1.index = cbind(aligned.mut.p.1.index, aligned.p.1.index[i])
    }
  }
  n.aligned.mut.p.1 = length(aligned.mut.p.1.index)
  
  identity = 100 - (n.aligned.mut.p.1 * 100 / n.aligned)
  
  ### Core ### 
  
  # Calculate core indexes.
  gaps = matrix(nrow = l.alignment, ncol = 2)
  gaps[, 2] = seq(1:l.alignment)
  for (i in (1:l.alignment)) {
    ngaps = 0
    for (j in (1:n.prot)) {
      is.gap = is.gap(alignment[j, i])
      if (is.gap == "TRUE") {
        ngaps = ngaps + 1
      }
    }
    gaps[i, 1] = ngaps
  }
  core.index = gaps[gaps[, 1] == 0, 2]
  n.core = length(core.index)
  
  # Get aligned and not.aligned indexes.
  core.p.1.index = p.1.index[core.index]
  core.p.2.index = p.2.index[core.index]

  no.core.p.1.index = p.1.index[-c(core.index)]
  no.core.p.2.index = p.2.index[-c(core.index)]
  no.core.p.1.index = no.core.p.1.index[!is.na(no.core.p.1.index)]
  no.core.p.2.index = no.core.p.2.index[!is.na(no.core.p.2.index)]
  
  # Create a list for the output.
  output = list(      "n.sites.p.1" = aa.p.1,  
                      "n.sites.p.2" = aa.p.2, 
                        "n.aligned" = n.aligned,
                "aligned.p.1.index" = aligned.p.1.index,
                "aligned.p.2.index" = aligned.p.2.index,
            "not.aligned.p.1.index" = not.aligned.p.1.index,
            "not.aligned.p.2.index" = not.aligned.p.2.index,
                "n.aligned.mut.p.1" = n.aligned.mut.p.1,
            "aligned.mut.p.1.index" = aligned.mut.p.1.index,
                           "n.core" = n.core,
                   "core.p.1.index" = core.p.1.index, 
                   "core.p.2.index" = core.p.2.index, 
                "no.core.p.1.index" = no.core.p.1.index,
                "no.core.p.2.index" = no.core.p.2.index,
                         "identity" = identity)
  output
}