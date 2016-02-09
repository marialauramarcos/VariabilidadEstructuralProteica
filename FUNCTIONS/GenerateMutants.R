# This function generates  mutants of a given protein using two possible models:
# - "LFENM" (Linearly Forced - Elastic Network Model) considering additive single mutations.
# - "MND" (Multivariate Normal Distribution) considering mu = dr expected value = 0 and Sigma = Cov.
#
#  Args:
#    p.ref: the pdb code (pdbid) of the protein to mutate.
#    mut.model: mutational model. It can be "LFENM" or "MND".
#    theo.chain.p.ref: the chain of p.ref to mutate.
#    n.mut.p: the number of mutants to generate for each protein of the family.
#    fmax: the maximun value for the forces that model the mutations.
#    R0: the cut-off for the ANM.
#    heme: argument for globins. It can be "TRUE" or "FALSE". If it is "TRUE" the program considers the hemo group.
#    data.dir: directory of the data. It must contain the file "_out_m.aligned.mut.p.1.index.csv" with indexes to mutate.
#    This file is generated using AnalyzeFamily().
#    out.dir: output directory.
#    out.name: ID for output filenames.
#    TOLERANCE: 0 tolerance.
#
#  Requires:
#    ReadCA
#    ReadHeme
#    CalculateENMK
#    CalculateForce
#
#  Returns:
#    File with a vector of coordinates of p.ref in out.dir.
#    File with a matrix of coordinates of the mutants in each column in out.dir.

GenerateMutants <- function(family,
                            p.ref, 
                            chain,
                            mut.model,
                            n.mut.p,
                            fmax, 
                            R0,
                            heme, 
                            data.dir,
                            out.dir,
                            out.name,
                            TOLERANCE) {

  # Filenames.
  pdb.fname <- file.path(data.dir, paste(p.ref, ".pdb", sep = ""))
  m.aligned.mut.p.1.index.fname <- file.path(out.dir, paste(family, "_out_m.aligned.mut.p.1.index.csv", sep = ""))
  
  # Get and read PDB of p.ref.
  get.pdb(as.character(p.ref), data.dir)
  
  pdb = ReadCA(pdb.fname, chain)
  r.p.ref = pdb$xyz.calpha
  n.aa = pdb$n.sites

  # Calculate heme coordinates, add them to CA´s coordinates and calculate the number of sites.
  if (heme == "TRUE") {
    r.heme = ReadHeme(pdb.fname, chain)
    r.p.ref = cbind(r.p.ref, r.heme)
    n.sites = ncol(r.p.ref)
  } else {
    n.sites = n.aa
  }

  # Calculate K of p.ref.
  ENMK.p.ref = CalculateENMK(r.p.ref, CalculateKij, R0, TOLERANCE)
  
  # Get mutated sites of p.ref in the alignment.
  m.aligned.mut.p.1.index = read.csv(m.aligned.mut.p.1.index.fname)
  n.prot = nrow(m.aligned.mut.p.1.index)

  # Calculate mutants using "LF-ENM".
  if (mut.model == "LFENM") {

    # Create a matrix to save coordinates of each generated mutant.
    m.r.mut = matrix(nrow = 3 * n.sites, ncol = n.prot * n.mut.p)

    # Start a loop to calculate coordinates of the mutants.
    for (P in (1:n.prot)) {
      
      # Get mutated sites of p.ref.
      aligned.mut.p.1.index = m.aligned.mut.p.1.index[P, ]
      aligned.mut.p.1.index = as.numeric(aligned.mut.p.1.index[, !is.na(aligned.mut.p.1.index)])
      
      for(mut in seq(n.mut.p)) {
        print(c(P,mut))
        f <- rep(0, 3 * n.sites)
        for (l in as.numeric(aligned.mut.p.1.index)) {
          fl = CalculateForce(l, r.p.ref, ENMK.p.ref$kij, fmax)
          f = f + fl
       }
       dr.mut = ENMK.p.ref$cov %*% f
       r.mut = as.vector(r.p.ref) + dr.mut
       m.r.mut[, (n.mut.p * P) - n.mut.p - mut] = r.mut
      }
    }
  }
  
  # Calculate mutants using "MND".
  if (mut.model == "MND") {
    
    # Calculate the expected dr value for each site.
    dr = matrix(0, nrow = 1, ncol = 3 * n.sites)
  
    # Create a matrix to save coordinates of each mutant.
    m.r.mut = matrix(0, nrow = 3 * n.sites, ncol = n.prot * n.mut.p)
  
    # Start a loop to calculate coordinates of the mutants.
    for (mut in seq(n.prot * n.mut.p)) {
      print(mut)
      dr.mut =  mvrnorm(n = 1, dr, ENMK.p.ref$cov)
      r.mut =  as.vector(r.p.ref) + dr.mut
      m.r.mut[, mut] = r.mut
    }
  }
  
  # Create files to save the data.
  write.csv(as.vector(r.p.ref), file = file.path(out.dir, paste(out.name, "_out_r.p.ref.csv", sep = "")), row.names = FALSE)
  write.csv(m.r.mut, file = file.path(out.dir, paste(out.name, "_out_m.r.mut.csv", sep = "")), row.names = FALSE)
}
  
  
