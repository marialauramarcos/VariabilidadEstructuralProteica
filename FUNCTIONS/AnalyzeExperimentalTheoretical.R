# Description:
#
# This function analyzes experimental and theoretical proteins calculating measures of structural variability in cartesian 
# coordinates and projected on the normal modes of a reference protein. 
# 
# Usage:
#
# AnalyzeExperimentalTheoretical <- function(family, chain.p.ref, n.mut.p, R0, rotate = TRUE/FALSE,
# heme = TRUE/FALSE, K.analysis, data.dir, out.dir, mut.fname.id, analysis.fname.id, tolerance)
#
#  Args:
#    - family: the family of the protein to mutate. It can be "globins", "serinProteases", 
#    "snakesToxin", "sh3", "fabp", "rrm", "phoslip" or "cys".
#    - chain.p.ref: the chain of p.ref in the pdb file obtained from Homstrad.
#    - n.mut.p: the number of mutants to generate for each member of the family. For example, if the family has 20 
#    members, the program generates n.mut.p x 20 mutants.
#    - R0: the cut-off for the ANM ("Anisotropic Network Model") that represents the proteins.
#    - rotate: it can be "TRUE" or "FALSE". If it is "TRUE", r.p.2 is rotaded in order to minimize RMSD with r.p.ref.
#    - heme: argument for globins. It can be "TRUE" or "FALSE". If it is "TRUE", the program considers the heme group. 
#    - K.analysis: It can be "K" or "Keff". For "K" or "Keff", the analysis is based on normal modes of "K" or "Keff"
#    respectibly.
#    - data.dir: directory of the data. It must contain the dataset ("data.dir/family_dataset.csv") and the pdb file 
#    obtained from Homstrad ("data.dir/family_coordinates.csv").
#    - out.dir: directory of the output. It must contain output files generated with AnalyzeFamily() and GenerateMutants().
#    The output of this function is also saved in out.dir.
#    - mut.fname.id: ID of filnames of mutant proteins.
#    - analysis.fname.id: ID of output filenames.
#    - tolerance: 0 tolerance.
#
#  Required libraries
#    {Bio3d}
#
#  Required functions
#    ReadCA
#    ReadHeme
#    CalculateVariability
#    CalculateENMKeff
#    CalculateENMK

AnalyzeExperimentalTheoretical <- function(family,
                                           chain.p.ref,
                                           n.mut.p,
                                           R0, 
                                           rotate,
                                           heme,
                                           K.analysis,
                                           data.dir,
                                           out.dir,
                                           mut.fname.id, 
                                           analysis.fname.id,
                                           tolerance) {
  ### GENERAL ###
  
  # Filenames
  dataset.fname <- file.path(data.dir, paste(family, "_dataset.csv", sep = ""))
  pdbs.fname <- file.path(data.dir, paste(family, "_coordinates.pdb", sep = "")) 
  
  m.n.aligned.fname <- file.path(out.dir, paste(family, "_out_m.n.aligned.csv", sep = ""))
  m.aligned.p.ref.index.fname <- file.path(out.dir, paste(family, "_out_m.aligned.p.ref.index.csv", sep = ""))
  m.aligned.p.2.index.fname <- file.path(out.dir, paste(family, "_out_m.aligned.p.2.index.csv", sep = ""))
  m.not.aligned.p.ref.index.fname <- file.path(out.dir, paste(family, "_out_m.not.aligned.p.ref.index.csv", sep = ""))
  m.not.aligned.p.2.index.fname <- file.path(out.dir, paste(family, "_out_m.not.aligned.p.2.index.csv", sep = ""))
  
  # Read the dataset
  dataset <- read.csv(dataset.fname)
  pdbid.dataset <- as.character(dataset$pdbid)
  chain <- as.character(dataset$chain)
  n.prot = length(pdbid.dataset)
  
  ### THEORETICAL ###
  
  # Filenames
  theo.r.p.ref.fname <- file.path(out.dir, paste(mut.fname.id, "_out_r.p.ref.csv", sep = ""))
  m.r.mut.fname <- file.path(out.dir, paste(mut.fname.id, "_out_m.r.mut.csv", sep = ""))
  
  # Read coordinates
  theo.r.p.ref = read.csv(theo.r.p.ref.fname)$x
  m.r.mut = read.csv(m.r.mut.fname)
  n.sites.p.ref = length(theo.r.p.ref)/3
  n.mut = ncol(m.r.mut)
  
  # Create matrices to save measures of variability of each mutant
  m.theo.Pn = matrix(nrow = n.mut, ncol = 3 * n.sites.p.ref)
  m.theo.va = matrix(nrow = n.mut, ncol = 3 * n.sites.p.ref)
  m.theo.dr.squarei = matrix(nrow = n.mut, ncol = n.sites.p.ref)
  m.theo.norm.dr.squarei = matrix(nrow = n.mut, ncol = n.sites.p.ref)
  m.theo.smooth.dr.squarei = matrix(nrow = n.mut, ncol = n.sites.p.ref)
  m.theo.smooth.norm.dr.squarei = matrix(nrow = n.mut, ncol = n.sites.p.ref)
  
  # Start a loop to analyze each mutant
  for (mut in (1:(n.mut))) {
    print(c(mut))
      
    # Get theo.r.p.2
    theo.r.p.2 = m.r.mut[, mut]
    
    # Calculate measures of variavility
    theo.variability = CalculateVariability(as.vector(theo.r.p.ref), 
                                            as.vector(theo.r.p.2), 
                                            seq(1, n.sites.p.ref),                                               
                                            seq(1, n.sites.p.ref), 
                                            c(),
                                            c(),
                                            R0,
                                            rotate,
                                            K.analysis,
                                            tolerance)
      
    m.theo.va[mut, 1:length(theo.variability$va)] = theo.variability$va
    m.theo.Pn[mut, 1:length(theo.variability$Pn)] = theo.variability$Pn
    m.theo.dr.squarei[mut, 1:length(theo.variability$dr.squarei)] = theo.variability$dr.squarei
      
    # Calculate norm.dr.squarei
    m.theo.norm.dr.squarei[mut, ] = m.theo.dr.squarei[mut, ]/ mean(m.theo.dr.squarei[mut, ], na.rm = T)
  
    # Calculate smooth.norm.dr.squarei
    kij = CalculateENMK(matrix(theo.r.p.ref, nrow = 3), CalculateKij, R0, tolerance)$kij
    m.theo.smooth.dr.squarei[mut, ] = (m.theo.dr.squarei[mut, ] + (kij %*%  m.theo.dr.squarei[mut, ])) / (rowSums(kij) + 1)     
    m.theo.smooth.norm.dr.squarei[mut, ] = m.theo.smooth.dr.squarei[mut, ]/ mean(m.theo.smooth.dr.squarei[mut, ], na.rm = T)
  }
  
  # Create files to save the data
  write.csv(m.theo.va, file = file.path(out.dir, paste(analysis.fname.id, "_out_m.theo.va.csv", sep = "")), row.names = FALSE)
  write.csv(m.theo.Pn, file = file.path(out.dir, paste(analysis.fname.id, "_out_m.theo.Pn.csv", sep = "")), row.names = FALSE)
  write.csv(m.theo.norm.dr.squarei, file = file.path(out.dir, paste(analysis.fname.id, "_out_m.theo.norm.dr.squarei.csv", sep = "")), row.names = FALSE)
  write.csv(m.theo.smooth.norm.dr.squarei, file = file.path(out.dir, paste(analysis.fname.id, "_out_m.theo.smooth.norm.dr.squarei.csv", sep = "")), row.names = FALSE)
  
  ### EXPERIMENTAL ###
  
  # Read pdb of exp.p.ref
  exp.r.p.ref = theo.r.p.ref
  
  if(heme == "TRUE") {
    n.aa.p.ref = n.sites.p.ref - 5 
  } else {
    n.aa.p.ref = n.sites.p.ref
  }
  
  # Read indexes files
  m.n.aligned = read.csv(m.n.aligned.fname)
  m.aligned.p.ref.index = read.csv(m.aligned.p.ref.index.fname)
  m.aligned.p.2.index = read.csv(m.aligned.p.2.index.fname)
  m.not.aligned.p.ref.index = read.csv(m.not.aligned.p.ref.index.fname)
  m.not.aligned.p.2.index = read.csv(m.not.aligned.p.2.index.fname)
  
  # Create matrices to save measures of variability of each mutant
  m.exp.Pn = matrix(nrow = n.prot, ncol = 3 * (n.sites.p.ref)) 
  m.exp.va = matrix(nrow = n.prot, ncol = 3 * (n.sites.p.ref))
  m.exp.dr.squarei = matrix(nrow = n.prot, ncol = (n.sites.p.ref))
  m.exp.norm.dr.squarei = matrix(nrow = n.prot, ncol = (n.sites.p.ref))
  m.exp.smooth.dr.squarei = matrix(nrow = n.prot, ncol = (n.sites.p.ref))
  m.exp.smooth.norm.dr.squarei = matrix(nrow = n.prot, ncol = (n.sites.p.ref))
  
  # Start a loop to evaluate each protein of the family
  for (P in (1:n.prot)) {
    print(P)
    
    # Get aligned and not aligned indexes
    n.aligned = as.numeric(m.n.aligned[P, ])
    aligned.p.ref.index = as.numeric(m.aligned.p.ref.index[P, !is.na(m.aligned.p.ref.index[P, ])])
    aligned.p.2.index = as.numeric(m.aligned.p.2.index[P, !is.na(m.aligned.p.2.index[P, ])])
    not.aligned.p.ref.index = as.numeric(m.not.aligned.p.ref.index[P, !is.na(m.not.aligned.p.ref.index[P, ])])
    not.aligned.p.2.index = as.numeric(m.not.aligned.p.2.index[P, !is.na(m.not.aligned.p.2.index[P, ])])
    
    # Read PDB of exp.p.2
    chain.p.2 <- chain[[P]]
    exp.pdb.p.2 = ReadCA(pdbs.fname, chain.p.2)
    exp.r.p.2 = exp.pdb.p.2$xyz.calpha
    exp.n.aa.p.2 = exp.pdb.p.2$n.sites
    
    # Calculate heme coordinates, add them to CA´s coordinates and calculate the number of sites and not aligned indexes
    if (heme == "TRUE") {
      exp.r.heme.p.2 = ReadHeme(pdbs.fname, chain.p.2)
      exp.r.p.2 = cbind(exp.r.p.2, exp.r.heme.p.2)
      exp.n.sites.p.2 = ncol(exp.r.p.2)
      
      aligned.p.ref.index <- c(aligned.p.ref.index, t(seq((n.aa.p.ref + 1), n.sites.p.ref)))
      aligned.p.2.index <- c(aligned.p.2.index, t(seq((exp.n.aa.p.2 + 1), exp.n.sites.p.2)))
    }
    
    # Calculate measures of variability
    exp.variability = CalculateVariability(as.vector(exp.r.p.ref), 
                                           as.vector(exp.r.p.2), 
                                           aligned.p.ref.index, 
                                           aligned.p.2.index, 
                                           not.aligned.p.ref.index,
                                           not.aligned.p.2.index,
                                           R0,
                                           rotate,
                                           K.analysis,
                                           tolerance)
    
    m.exp.va[P, 1:length(exp.variability$va)] = exp.variability$va
    m.exp.Pn[P, 1:length(exp.variability$Pn)] = exp.variability$Pn
    exp.dr.squarei = exp.variability$dr.squarei
    for (i in (1:n.sites.p.ref)) {
      m.exp.dr.squarei[P, i] = matrix(exp.dr.squarei[aligned.p.ref.index == i], nrow = 1, ncol = 1)
    }
    
    # Calculate norm.dr.squarei
    m.exp.norm.dr.squarei[P, ] = m.exp.dr.squarei[P, ]/ mean(m.exp.dr.squarei[P, ], na.rm = T)
    
    # Calculate smooth.norm.dr.squarei
    m.exp.dr.squarei.0 = m.exp.dr.squarei[P, ]
    m.exp.dr.squarei.0[is.na(m.exp.dr.squarei.0)] = tolerance
    m.exp.smooth.dr.squarei[P, ] = (m.exp.dr.squarei[P, ] + (kij %*%  m.exp.dr.squarei.0)) / (rowSums(kij[, !is.na(m.exp.dr.squarei[P, ])]) + 1)
    m.exp.smooth.norm.dr.squarei[P, ] = m.exp.smooth.dr.squarei[P, ]/ mean(m.exp.smooth.dr.squarei[P, ], na.rm = T)
  }  
  
  # Create files to save the data
  write.csv(m.exp.va, file = file.path(out.dir, paste(analysis.fname.id, "_out_m.exp.va.csv", sep = "")), row.names = FALSE)
  write.csv(m.exp.Pn, file = file.path(out.dir, paste(analysis.fname.id, "_out_m.exp.Pn.csv", sep = "")), row.names = FALSE)
  write.csv(m.exp.norm.dr.squarei, file = file.path(out.dir, paste(analysis.fname.id, "_out_m.exp.norm.dr.squarei.csv", sep = "")), row.names = FALSE)
  write.csv(m.exp.smooth.norm.dr.squarei, file = file.path(out.dir, paste(analysis.fname.id, "_out_m.exp.smooth.norm.dr.squarei.csv", sep = "")), row.names = FALSE)
}






