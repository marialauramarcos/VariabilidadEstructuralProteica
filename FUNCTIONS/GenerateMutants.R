# Description:
#
# This function generates multiple mutants of a given protein using the LFENM ("Linearly Forced - Elastic Network Model")
# and considering additive single mutations. It calculates the probability of acceptance of each multuple mutant
# according to the "Stress Model" using 4 different values of the parameter "beta" for no selection, weak selection,
# meduim selection and stong selection.
#
# Usage:
#
# GenerateMutants(family, chain.p.ref, n.mut.p, fmax, R0, heme = TRUE/FALSE, data.dir, out.dir,
# mut.fname.id, TOLERANCE)
#
#  Args:
#    family: the family of the protein to mutate. It can be "globins", "serinProteases", 
#    "snakesToxin", "sh3", "fabp", "rrm", "phoslip" or "cys".
#    chain.p.ref: the chain of p.ref in the pdb file obtained from Homstrad.
#    n.mut.p: the number of mutants to generate for each member of the family. For example, if the family has 20 
#    members, the program generates n.mut.p x 20 mutants.
#    fmax: It is the maximun value for the forces that model the mutations.
#    R0: the cut-off for the ANM ("Anisotropic Network Model") that represents the proteins.
#    heme: argument for "globins". It can be "TRUE" or "FALSE". If it is "TRUE", the program considers the heme group.
#    data.dir: directory of the data. It must contain the pdb file obtained from Homstrad ("data.dir/family_coordinates.csv").
#    out.dir: output directory. It must contain the output of the function AnalyzeFamily().
#    mut.fname.id: ID for output filenames.
#    TOLERANCE: 0 tolerance.
#
#  Required libraries:
#    {Bio3d}
#
#  Required function:
#    ReadCA()
#    ReadHeme()
#    CalculateENMK()
#    CalculateForce()

GenerateMutants <- function(family,
                            chain.p.ref,
                            n.mut.p,
                            fmax, 
                            R0,
                            heme,
                            data.dir,
                            out.dir,
                            mut.fname.id,
                            TOLERANCE) {
  
  # Filenames.
  pdbs.fname <- file.path(data.dir, paste(family, "_coordinates.pdb", sep = "")) 
  m.identity.fname <- file.path(out.dir, paste(family, "_out_m.identity.csv", sep = ""))
  
  # Read PDB of p.ref.
  pdb = ReadCA(pdbs.fname, chain.p.ref)
  r.p.ref = pdb$xyz.calpha
  n.aa = pdb$n.sites
  
  # Calculate heme coordinates, add them to CA´s coordinates and calculate the number of sites.
  if (heme == "TRUE") {
    r.heme = ReadHeme(pdbs.fname, chain.p.ref)
    r.p.ref = cbind(r.p.ref, r.heme)
    n.sites = ncol(r.p.ref)
  } else {
    n.sites = n.aa
  }
  
  # Calculate K of p.ref.
  ENMK.p.ref = CalculateENMK(r.p.ref, CalculateKij, R0, TOLERANCE)
  
  # Read data of the family: % sequence identity and number of proteins.
  m.identity = read.csv(m.identity.fname)$V1
  n.prot = length(m.identity)
  
  # Calculate mean identity and the number of sites to mutate.
  identity = mean(m.identity)
  n.sites.mut = (100 - (identity)) * n.aa / 100
  
  # Create a matrices to save calculated data.
  m.r.mut = matrix(0, nrow = 3 * n.sites, ncol = n.prot * n.mut.p)
  m.fl.sum.sum.squares = matrix(0, nrow = 1, ncol = n.prot * n.mut.p)
  
  # Start a loop to generate n.mut.p mutants for each protein of teh family.
  for(mut in (1:(n.prot * n.mut.p))) {
    print(c(mut))
      
    # Mutate random sites. 
    mutated.index = sample(1:n.aa, replace = F)[1:n.sites.mut]
    f = rep(0, 3 * n.sites)
    fl.sum.squares = c()
      
    # Start a loop to calculate each mutation at site l.
    for (l in as.numeric(mutated.index)) {
      fl = CalculateForce(l, r.p.ref, ENMK.p.ref$kij, fmax)
      f = f + fl
      fl.sum.squares = c(fl.sum.squares, sum(fl ^ 2)) 
    }
      
    # Calculate and Keep fl.sum.sum.squares.
    m.fl.sum.sum.squares[, mut] = sum(fl.sum.squares)
      
    # Calculate dr and r.mut.
    dr.mut = ENMK.p.ref$cov %*% f
    r.mut = as.vector(r.p.ref) + dr.mut
      
    # Keep r.mut.
    m.r.mut[, mut] = r.mut
  }
  
  # Calculate betas for the "Stress Model".
  beta.0.1 = - log(0.1) * 2 / mean(m.fl.sum.sum.squares)
  beta.0.5 = - log(0.5) * 2 / mean(m.fl.sum.sum.squares)
  beta.0.9 = - log(0.9) * 2 / mean(m.fl.sum.sum.squares)
  betas = c(beta.0.1, beta.0.5, beta.0.9)
  
  # Calculate the acceptance probability of each mutant.
  prob.accept.mut.all.betas = rbind(exp(- beta.0.1 / 2 * (m.fl.sum.sum.squares)),
                                    exp(- beta.0.5 / 2 * (m.fl.sum.sum.squares)),
                                    exp(- beta.0.9 / 2 * (m.fl.sum.sum.squares)))
  
  # Create files to save the data.
  write.csv(as.vector(r.p.ref), file = file.path(out.dir, paste(mut.fname.id, "_out_r.p.ref.csv", sep = "")), row.names = FALSE)
  write.csv(m.r.mut, file = file.path(out.dir, paste(mut.fname.id, "_out_m.r.mut.csv", sep = "")), row.names = FALSE)
  write.csv(m.fl.sum.sum.squares, file = file.path(out.dir, paste(mut.fname.id, "_out_m.fl.sum.sum.squares.csv", sep = "")), row.names = FALSE)
  write.csv(betas, file = file.path(out.dir, paste(mut.fname.id, "_out_betas.csv", sep = "")), row.names = FALSE)
  write.csv(prob.accept.mut.all.betas, file = file.path(out.dir, paste(mut.fname.id, "_out_prob.accept.mut.all.betas.csv", sep = "")), row.names = FALSE)
}




