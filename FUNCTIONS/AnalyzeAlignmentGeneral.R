AnalyzeAlignmentGeneral <- function(df.alignment)UseMethod("AnalyzeAlignmentGeneral") 

AnalyzeAlignmentGeneral.NoCore <- function(df.alignment) {
  
  # Prepare data.
  alignment <- df.alignment$alignment
  id <- df.alignment$id
  p.1 <- df.alignment$p.1
  p.2 <- df.alignment$p.2
  
  # Extract data form the alignment.
  lalignment = ncol(alignment)
  alignment.p.1 = alignment[id == p.1,]
  alignment.p.2 = alignment[id == p.2,]

  # Get indexes.
  p.1.index = c()
  p.2.index = c()
  aa.p.1 = 1
  aa.p.2 = 1
  
  for (i in (1:lalignment)) {
    if (alignment.p.1[i] != "-") {
      p.1.index[i] <- aa.p.1
      aa.p.1 = aa.p.1 + 1
    }
    if (alignment.p.2[i] != "-") {
      p.2.index[i] <- aa.p.2
      aa.p.2 = aa.p.2 + 1
    }
  }	
  
  # Get aligned and not.aligned indexes.
  aligned.p.1.index = c()
  aligned.p.2.index = c()
  not.aligned.p.1.index = c()
  not.aligned.p.2.index = c()
  
  for (i in (1:lalignment)) {
    if (alignment.p.1[i] != "-" & alignment.p.2[i] != "-") {
      aligned.p.1.index <- cbind(aligned.p.1.index, p.1.index[i])
      aligned.p.2.index <- cbind(aligned.p.2.index, p.2.index[i])
    }
    if (alignment.p.1[i] != "-" & alignment.p.2[i] == "-") {
      not.aligned.p.1.index <- cbind(not.aligned.p.1.index, p.1.index[i])
    }
    if (alignment.p.2[i] != "-" & alignment.p.1[i] == "-") {
      not.aligned.p.2.index <- cbind(not.aligned.p.2.index, p.2.index[i])
    }
  }	
  naligned <- length(aligned.p.1.index)
  output = list("aligned.p.1.index" = aligned.p.1.index,
                "aligned.p.2.index" = aligned.p.2.index,
                "not.aligned.p.1.index" = not.aligned.p.1.index,
                "not.aligned.p.2.index" = not.aligned.p.2.index,
                "naligned" = naligned)
  output
}

AnalyzeAlignmentGeneral.Core <- function(df.alignment) {

  # Prepare data.
  alignment <- df.alignment$alignment
  id <- df.alignment$id
  p.1 <- df.alignment$p.1
  p.2 <- df.alignment$p.2
  
  # Extract data form the alignment.
  nprot = nrow(alignment)
  lalignment = ncol(alignment)
  alignment.p.1 = alignment[id == p.1,]
  alignment.p.2 = alignment[id == p.2,]
  
  # Calculate core.index.
  gaps = matrix(nrow = lalignment, ncol = 2)
  gaps[, 2] = seq(1:lalignment)
  for (i in (1:lalignment)) {
    ngaps = 0
    for (j in (1:nprot)) {
      is.gap = is.gap(alignment[j, i])
      if (is.gap == "TRUE") {
        ngaps = ngaps + 1
      }
    }
    gaps[i, 1] = ngaps
  }
  core.index = gaps[gaps[, 1] == 0, 2]
  
  # Get indexes.
  p.1.index = c()
  p.2.index = c()
  aa.p.1 = 1
  aa.p.2 = 1
  
  for (i in (1:lalignment)) {
    if (alignment.p.1[i] != "-") {
      p.1.index[i] <- aa.p.1
      aa.p.1 = aa.p.1 + 1
    }
    if (alignment.p.2[i] != "-") {
      p.2.index[i] <- aa.p.2
      aa.p.2 = aa.p.2 + 1
    }
  }	
  
  # Get aligned and not.aligned indexes.
  core.index.p.1 = p.1.index[core.index]
  core.index.p.2 = p.2.index[core.index]
  naligned <- length(core.index.p.1)
  
  no.core.index.p.1 = p.1.index[-c(core.index)]
  no.core.index.p.2 = p.2.index[-c(core.index)]
  no.core.index.p.1 = no.core.index.p.1[!is.na(no.core.index.p.1)]
  no.core.index.p.2 = no.core.index.p.2[!is.na(no.core.index.p.2)]
  
  output = list("aligned.p.1.index" = core.index.p.1, 
              "aligned.p.2.index" = core.index.p.2,
              "not.aligned.p.1.index" = no.core.index.p.1,
              "not.aligned.p.2.index" = no.core.index.p.2,
              "naligned" = naligned)
  output
}