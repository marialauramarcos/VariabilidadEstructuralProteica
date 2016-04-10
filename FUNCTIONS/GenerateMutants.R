# This function generates mutants of a given protein using two possible models:
# - "LFENM" (Linearly Forced - Elastic Network Model), considering additive single mutations.
# - "MND" (Multivariate Normal Distribution), considering mu = dr expected value = 0 and Sigma = Cov.
#
#  Args:
#    family: the family of the protein to mutate. It can be "globins", "serinProteases", 
#    "snakesToxin", "sh3", "fabp", "rrm", "phoslip" or "cys".
#    exp.chain.p.ref: the chain of p.ref in the pdb file obtained from Homstrad.
#    mut.model: mutational model. It can be "LFENM" (Linearly Forced - Elastic Network Model) or "MND" 
#    (Multivariate Normal Distribution).
#    n.mut.p: the number of mutants to generate for each member of the family. For example, if the family has 20 
#    members, the program generates n.mut.p x 20 mutants.
#    fmax: argument for "LFENM". It is the maximun value for the forces that model the mutations.
#    R0: the cut-off for the "ANM" (Anisotropic Network Model) that represents the proteins.
#    heme: argument for "globins". It can be "TRUE" or "FALSE". If it is "TRUE", the program considers the heme group.
#    natural.selection: It can be "TRUE" or "FALSE". If it is "TRUE", the mutants are calculated considering natural 
#    selection. If it is "FALSE", the mutants are calculated in a random manner.
#    data.dir: directory of the data. It must contain the file the pdb file 
#    obtained from Homstrad ("data.dir/family_coordinates.csv").
#    out.dir: output directory. It must contain the output of the function AnalyzeFamily().
#    mut.fname.id: ID for output filenames.
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
                            exp.chain.p.ref,
                            mut.model,
                            n.mut.p,
                            fmax, 
                            R0,
                            heme,
                            natural.selection,
                            data.dir,
                            out.dir,
                            mut.fname.id,
                            TOLERANCE) {
  
  # Filenames.
  pdbs.fname <- file.path(data.dir, paste(family, "_coordinates.pdb", sep = "")) 
  m.identity.fname <- file.path(out.dir, paste(family, "_out_m.identity.csv", sep = ""))
  
  # Read PDB of p.ref.
  pdb = ReadCA(pdbs.fname, exp.chain.p.ref)
  r.p.ref = pdb$xyz.calpha
  n.aa = pdb$n.sites
  
  # Calculate heme coordinates, add them to CA´s coordinates and calculate the number of sites.
  if (heme == "TRUE") {
    r.heme = ReadHeme(pdbs.fname, exp.chain.p.ref)
    r.p.ref = cbind(r.p.ref, r.heme)
    n.sites = ncol(r.p.ref)
  } else {
    n.sites = n.aa
  }
  
  # Enumerate sites of p.ref.
  sites = seq(1:n.sites)
  
  # Calculate K of p.ref.
  ENMK.p.ref = CalculateENMK(r.p.ref, CalculateKij, R0, TOLERANCE)
  
  # Get the % sequence identity between p.ref and the other proteins and get the number of proteins of the family.
  m.identity = read.csv(m.identity.fname)$V1
  n.prot = length(m.identity)
  
  # Create a matrix to save coordinates of each mutant.
  m.r.mut = matrix(0, nrow = 3 * n.sites, ncol = n.prot * n.mut.p)
  
  # Calculate mutants using "LF-ENM".
  if (mut.model == "LFENM") {

    # Start a loop for each P.
    for (P in (1:n.prot)) {
      
      # Get de sequence identity and the number of mutated sites of p.ref for P.
      identity = m.identity[P]
      n.sites.mut = (100 - (identity)) * n.aa / 100
      
      # Decide which sites to mutate for natural.selection == "TRUE".
      if (natural.selection == "TRUE") {
      
        # Calculate the number of contacts of each site.
        CN = ENMK.p.ref$kij
        CN.i = colSums(CN)

        # Calculate the probability of mutation of each site following the Stress Model (with beta = 1).
        prob.i = 1 - CN.i
      
        # Get sites with more probability of mutation.
        decreasing.prob = order(prob.i, decreasing = T)
        mutated.index = decreasing.prob[1:n.sites.mut]
      }
      
      # Start a loop to generate n.mut.p mutants for each P.
      for(mut in seq(n.mut.p)) {
        print(c(P, mut))
        
        # Get mutated sites of p.ref for each P and for natural.selection == "FALSE".
        if (natural.selection == "FALSE") {
          mutated.index = sample(1:n.aa, replace = F)[1:n.sites.mut]
        }
        
        # Calculate forces.
        f = rep(0, 3 * n.sites)
        for (l in as.numeric(mutated.index)) {
          fl = CalculateForce(l, r.p.ref, ENMK.p.ref$kij, fmax)
          f = f + fl
        }
        
        # Calculate dr and r.mut.
        dr.mut = ENMK.p.ref$cov %*% f
        r.mut = as.vector(r.p.ref) + dr.mut
        
        # Keep r.mut.
        m.r.mut[, n.mut.p * P - (n.mut.p - mut)] = r.mut
      }
    }
  }
  
  # Calculate mutants using "MND".
  if (mut.model == "MND") {
    
    # Calculate the expected dr value for each site.
    dr = matrix(0, nrow = 1, ncol = 3 * n.sites)
    
    # Start a loop to calculate coordinates of the mutants.
    for (mut in seq(n.prot * n.mut.p)) {
      print(mut)
      
      # Calculate dr and r.mut.
      dr.mut =  mvrnorm(n = 1, dr, ENMK.p.ref$cov)
      r.mut =  as.vector(r.p.ref) + dr.mut
      
      # Keep r.mut.
      m.r.mut[, mut] = r.mut
    }
  }
  
  # Create files and save the data.
  write.csv(as.vector(r.p.ref), file = file.path(out.dir, paste(mut.fname.id, "_out_r.p.ref.csv", sep = "")), row.names = FALSE)
  write.csv(m.r.mut, file = file.path(out.dir, paste(mut.fname.id, "_out_m.r.mut.csv", sep = "")), row.names = FALSE)
}