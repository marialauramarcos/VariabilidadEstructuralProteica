# This function analyzes a multiple alignment of a family of proteins. It has to possibilitys:
#  - core: aligned and not aligned sites of p.1 and p.2 that are in the 
# conserved core of the alignment: positions with no gaps in whole the alignemnt.
#  - noCore: all aligned and not aligned sites of p.1 and p.2.
#
#  Args:
#   family: the family of p.ref.
#   p.ref: reference protein.
#   data.dir: directory of the data. It must contain the alignment and the dataset. 
#   out.dir: output directory.
#
#  Requires:
#   AnalyzeAlignment()
#
#  Returns:
#   Files with the result of the analyisis in out.dir .

AnalyzeFamily <- function(family, p.ref, data.dir, out.dir) {
  
  # Data filenames.
  dataset.fname <- file.path(data.dir, paste(family, "_dataset.csv", sep = ""))  
  alignment.fname <- file.path(data.dir, paste(family, "_alignment.txt", sep = ""))  
  
  # Read dataset.
  dataset <- read.csv(dataset.fname)
  pdbid.dataset <- as.character(dataset$pdbid) 
  n.prot = length(pdbid.dataset)
  
  # Read alignment.
  alignment.id <- ReadFasta(alignment.fname)
  alignment <- alignment.id$ali[, -ncol(alignment.id$ali)]  # The last column is "*".
  pdbid.alignment <- alignment.id$id
  
  # Create matrices to save data.
  m.n.sites.p.1 = matrix(nrow = n.prot, ncol = 1)
  m.n.sites.p.2 = matrix(nrow = n.prot, ncol = 1)
  m.n.aligned = matrix(nrow = n.prot, ncol = 1)
  m.aligned.p.1.index = matrix(nrow = n.prot, ncol = 500)
  m.aligned.p.2.index = matrix(nrow = n.prot, ncol = 500)
  m.not.aligned.p.1.index = matrix(nrow = n.prot, ncol = 500)
  m.not.aligned.p.2.index = matrix(nrow = n.prot, ncol = 500)
  m.n.aligned.mut.p.1 = matrix(nrow = n.prot, ncol = 1)
  m.aligned.mut.p.1.index = matrix(nrow = n.prot, ncol = 500)
  m.n.core = matrix(nrow = n.prot, ncol = 1)
  m.core.p.1.index = matrix(nrow = n.prot, ncol = 500) 
  m.core.p.2.index = matrix(nrow = n.prot, ncol = 500)
  m.no.core.p.1.index = matrix(nrow = n.prot, ncol = 500)
  m.no.core.p.2.index = matrix(nrow = n.prot, ncol = 500)
  m.identity = matrix(nrow = n.prot, ncol = 1)
  
  # Start a loop to evaluate each protein of the family.
  for (P in (1:n.prot)) {
    print(P)
    
    # Get pdbid of p.2.
    p.2 <- pdbid.dataset[P]
    
    # Anylize the alignment.
    analysis.alignment = AnalyzeAlignment(alignment, 
                                          pdbid.alignment,   
                                          p.ref, 
                                          p.2)
    
    # Save analysis in the matrices.
    m.n.sites.p.1[P, ] = analysis.alignment$n.sites.p.1
    m.n.sites.p.2[P, ] = analysis.alignment$n.sites.p.2
    m.n.aligned[P, ] = analysis.alignment$n.aligned
    m.aligned.p.1.index[P, 1:length(analysis.alignment$aligned.p.1.index)] = analysis.alignment$aligned.p.1.index
    m.aligned.p.2.index[P, 1:length(analysis.alignment$aligned.p.2.index)] = analysis.alignment$aligned.p.2.index
    if (length(analysis.alignment$not.aligned.p.1.index) > 0) {
      m.not.aligned.p.1.index[P, 1:length(analysis.alignment$not.aligned.p.1.index)] = analysis.alignment$not.aligned.p.1.index
    }
    if (length(analysis.alignment$not.aligned.p.2.index) > 0) {
      m.not.aligned.p.2.index[P, 1:length(analysis.alignment$not.aligned.p.2.index)] = analysis.alignment$not.aligned.p.2.index
    }
    m.n.aligned.mut.p.1[P, ] = analysis.alignment$n.aligned.mut.p.1
    if (analysis.alignment$n.aligned.mut.p.1 > 0) {
      m.aligned.mut.p.1.index[P, 1:length(analysis.alignment$aligned.mut.p.1.index)] = analysis.alignment$aligned.mut.p.1.index
    }
    m.n.core[P, ] = analysis.alignment$n.core
    m.core.p.1.index[P, 1:length(analysis.alignment$core.p.1.index)] = analysis.alignment$core.p.1.index 
    m.core.p.2.index[P, 1:length(analysis.alignment$core.p.2.index)] = analysis.alignment$core.p.2.index
    m.no.core.p.1.index[P, 1:length(analysis.alignment$no.core.p.1.index)] = analysis.alignment$no.core.p.1.index
    m.no.core.p.2.index[P, 1:length(analysis.alignment$no.core.p.2.index)] = analysis.alignment$no.core.p.2.index
    m.identity[P, ] = analysis.alignment$identity
  }
  
  
  # Create files to save data.
  write.csv(m.n.sites.p.1, file = file.path(out.dir, paste(family, "_out_m.n.sites.p.1.csv", sep = "")), row.names = FALSE)
  write.csv(m.n.sites.p.2, file = file.path(out.dir, paste(family, "_out_m.n.sites.p.2.csv", sep = "")), row.names = FALSE)
  write.csv(m.n.aligned, file = file.path(out.dir, paste(family, "_out_m.n.aligned.csv", sep = "")), row.names = FALSE)
  write.csv(m.aligned.p.1.index, file = file.path(out.dir, paste(family, "_out_m.aligned.p.1.index.csv", sep = "")), row.names = FALSE)
  write.csv(m.aligned.p.2.index, file = file.path(out.dir, paste(family, "_out_m.aligned.p.2.index.csv", sep = "")), row.names = FALSE)
  write.csv(m.not.aligned.p.1.index, file = file.path(out.dir, paste(family, "_out_m.not.aligned.p.1.index.csv", sep = "")), row.names = FALSE)
  write.csv(m.not.aligned.p.2.index, file = file.path(out.dir, paste(family, "_out_m.not.aligned.p.2.index.csv", sep = "")), row.names = FALSE)
  write.csv(m.n.aligned.mut.p.1, file = file.path(out.dir, paste(family, "_out_m.n.aligned.mut.p.1.csv", sep = "")), row.names = FALSE)
  write.csv(m.aligned.mut.p.1.index, file = file.path(out.dir, paste(family, "_out_m.aligned.mut.p.1.index.csv", sep = "")), row.names = FALSE)
  write.csv(m.n.core, file = file.path(out.dir, paste(family, "_out_m.n.core.csv", sep = "")), row.names = FALSE)
  write.csv(m.core.p.1.index, file = file.path(out.dir, paste(family, "_out_m.core.p.1.index.csv", sep = "")), row.names = FALSE)
  write.csv(m.core.p.2.index, file = file.path(out.dir, paste(family, "_out_m.core.p.2.index.csv", sep = "")), row.names = FALSE)
  write.csv(m.no.core.p.1.index, file = file.path(out.dir, paste(family, "_out_m.no.core.p.1.index.csv", sep = "")), row.names = FALSE)
  write.csv(m.no.core.p.2.index, file = file.path(out.dir, paste(family, "_out_m.no.core.p.2.index.csv", sep = "")), row.names = FALSE)
  write.csv(m.identity, file = file.path(out.dir, paste(family, "_out_m.identity.csv", sep = "")), row.names = FALSE)
}