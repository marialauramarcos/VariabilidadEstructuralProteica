# This function analises experimental and theoretical data calculating measures of variability 
# of experimental proteins or simulated mutants respectibly. 
# 
#  Args:
#    family: the family of p.ref.
#    exp.chain.p.ref: chain of p.ref in the pdb file containing rotated coordinates of experimental proteins.
#    n.mut.p: the number of simulated mutants to analyze for each protein of the family.
#    R0: cut-off for the ANM.
#    rotate: it can be "TRUE" or "FALSE". If it is "TRUE" r.p.2 is rotaded in order to minimize RMSD with r.p.ref.
#    core: it can be "TRUE" or "FALSE". If it is "TRUE" the program only considers the conserved core of 
#    the alignment. If it is "FALSE" the program analyzes the whole alignment.
#    heme: argument for globins. It can be "TRUE" or "FALSE". If it is "TRUE" the program considers the hemo. 
#    data.dir: directory of the data. It must contain the dataset and pdb file with rotated coordinates of the experimental proteins.
#    out.dir: directory of the output. It must contain files generated with AnalyzeFamily() and GenerateMutants(). The output of this program 
#    is also saved in out.dir.
#    mut.name: ID of filnames of mutant proteins.
#    out.name: ID of output filenames.
#    TOLERANCE: 0 tolerance.
#
#  Requires:
#    ReadCA
#    ReadHeme
#    CalculateVariability
#    CalculateENMKeff
#    CalculateENMK
#
#  Returns:
#    File with m.exp.va in out.dir.
#    File with m.exp.Pn in out.dir.
#    File with m.exp.dr.squarei in out.dir.
#    File with m.theo.va in out.dir.
#    File with m.theo.Pn in out.dir.
#    File with m.theo.dr.squarei in out.dir.

AnalyzeExperimentalTheoretical <- function(family,
                                           exp.chain.p.ref,
                                           n.mut.p,
                                           R0, 
                                           core,
                                           rotate,
                                           heme,
                                           data.dir,
                                           out.dir,
                                           mut.name, 
                                           out.name,
                                           TOLERANCE) {
  # Filenames.
  dataset.fname <- file.path(data.dir, paste(family, "_dataset.csv", sep = ""))
  pdbs.fname <- file.path(data.dir, paste(family, "_coordinates.pdb", sep = "")) 
  
  m.n.aligned.fname <- file.path(out.dir, paste(family, "_out_m.n.aligned.csv", sep = ""))
  m.aligned.p.1.index.fname <- file.path(out.dir, paste(family, "_out_m.aligned.p.1.index.csv", sep = ""))
  m.aligned.p.2.index.fname <- file.path(out.dir, paste(family, "_out_m.aligned.p.2.index.csv", sep = ""))
  m.not.aligned.p.1.index.fname <- file.path(out.dir, paste(family, "_out_m.not.aligned.p.1.index.csv", sep = ""))
  m.not.aligned.p.2.index.fname <- file.path(out.dir, paste(family, "_out_m.not.aligned.p.2.index.csv", sep = ""))
  
  m.n.core.fname <- file.path(out.dir, paste(family, "_out_m.n.core.csv", sep = ""))
  m.core.p.1.index.fname <- file.path(out.dir, paste(family, "_out_m.core.p.1.index.csv", sep = ""))
  m.core.p.2.index.fname <- file.path(out.dir, paste(family, "_out_m.core.p.2.index.csv", sep = ""))
  m.no.core.p.1.index.fname <- file.path(out.dir, paste(family, "_out_m.no.core.p.1.index.csv", sep = ""))
  m.no.core.p.2.index.fname <- file.path(out.dir, paste(family, "_out_m.no.core.p.2.index.csv", sep = ""))
  
  m.r.mut.fname <- file.path(out.dir, paste(mut.name, "_out_m.r.mut.csv", sep = ""))
  theo.r.p.ref.fname <- file.path(out.dir, paste(mut.name, "_out_r.p.ref.csv", sep = ""))
  
  # Read dataset.
  dataset <- read.csv(dataset.fname)
  pdbid.dataset <- as.character(dataset$pdbid) 
  exp.chain <- as.character(dataset$chain)
  n.prot = length(pdbid.dataset)
  
  # Read pdb of exp.p.ref.
  exp.pdb.p.ref = ReadCA(pdbs.fname, exp.chain.p.ref)
  exp.r.p.ref = exp.pdb.p.ref$xyz.calpha
  n.aa.p.ref = exp.pdb.p.ref$n.sites
  
  # Read indexes files.
  m.n.aligned = read.csv(m.n.aligned.fname)
  m.aligned.p.1.index = read.csv(m.aligned.p.1.index.fname)
  m.aligned.p.2.index = read.csv(m.aligned.p.2.index.fname)
  m.not.aligned.p.1.index = read.csv(m.not.aligned.p.1.index.fname)
  m.not.aligned.p.2.index = read.csv(m.not.aligned.p.2.index.fname)
  
  m.n.core = read.csv(m.n.core.fname)
  m.core.p.1.index = read.csv(m.core.p.1.index.fname)
  m.core.p.2.index = read.csv(m.core.p.2.index.fname)
  m.no.core.p.1.index = read.csv(m.no.core.p.1.index.fname)
  m.no.core.p.2.index = read.csv(m.no.core.p.2.index.fname)
  
  # Read coordinates of theo.p.ref and the simulated mutants.
  m.r.mut = read.csv(m.r.mut.fname)
  theo.r.p.ref = read.csv(theo.r.p.ref.fname)$x
  n.sites.p.ref = length(theo.r.p.ref)/3

  # Create matrices to save measures of variability of each mutant.
  m.exp.Pn = matrix(nrow = n.prot, ncol = 3 * n.sites.p.ref)
  m.exp.va = matrix(nrow = n.prot, ncol = 3 * n.sites.p.ref)
  m.exp.dr.squarei = matrix(nrow = n.prot, ncol = n.sites.p.ref)
  
  m.theo.Pn = matrix(nrow = n.prot * n.mut.p, ncol = 3 * n.sites.p.ref)
  m.theo.va = matrix(nrow = n.prot * n.mut.p, ncol = 3 * n.sites.p.ref)
  m.theo.dr.squarei = matrix(nrow = n.prot * n.mut.p, ncol = n.sites.p.ref)
  
  # Start a loop to evaluate each protein of the family.
  for (P in (1:n.prot)) {
    print(P)
    
    # Get aligned and not aligned indexes.
    if (core == "FALSE") {
      n.aligned = as.numeric(m.n.aligned[P, ])
      aligned.p.1.index = as.numeric(m.aligned.p.1.index[P, !is.na(m.aligned.p.1.index[P, ])])
      aligned.p.2.index = as.numeric(m.aligned.p.2.index[P, !is.na(m.aligned.p.2.index[P, ])])
      not.aligned.p.1.index = as.numeric(m.not.aligned.p.1.index[P, !is.na(m.not.aligned.p.1.index[P, ])])
      not.aligned.p.2.index = as.numeric(m.not.aligned.p.2.index[P, !is.na(m.not.aligned.p.2.index[P, ])])
    } else {
      n.aligned = as.numeric(m.n.core[P, ]) 
      aligned.p.1.index = as.numeric(m.core.p.1.index[P, !is.na(m.core.p.1.index[P, ])])
      aligned.p.2.index = as.numeric(m.core.p.2.index[P, !is.na(m.core.p.2.index[P, ])])
      not.aligned.p.1.index = as.numeric(m.no.core.p.1.index[P, !is.na(m.no.core.p.1.index[P, ])])
      not.aligned.p.2.index = as.numeric(m.no.core.p.2.index[P, !is.na(m.no.core.p.2.index[P, ])])
    }
    
    # Read PDB of exp.p.2.  
    exp.chain.p.2 <- exp.chain[[P]]
    exp.pdb.p.2 = ReadCA(pdbs.fname, exp.chain.p.2)
    exp.r.p.2 = exp.pdb.p.2$xyz.calpha
    exp.n.aa.p.2 = exp.pdb.p.2$n.sites
    
    # Calculate heme coordinates, add them to CA´s coordinates and calculate the number of sites and not aligned indexes.
    if (heme == "TRUE") {
      if (P == 1) {
        exp.r.heme.p.ref = ReadHeme(pdbs.fname, exp.chain.p.ref)
        exp.r.p.ref = cbind(exp.r.p.ref, exp.r.heme.p.ref)
        n.sites.p.ref = ncol(exp.r.p.ref)
      }
      exp.r.heme.p.2 = ReadHeme(pdbs.fname, exp.chain.p.2)
      exp.r.p.2 = cbind(exp.r.p.2, exp.r.heme.p.2)
      exp.n.sites.p.2 = ncol(exp.r.p.2)
      
      not.aligned.p.1.index <- c(not.aligned.p.1.index, t(seq((n.aa.p.ref + 1), n.sites.p.ref)))
      not.aligned.p.2.index <- c(not.aligned.p.2.index, t(seq((exp.n.aa.p.2 + 1), exp.n.sites.p.2)))
    }
    
    # Calculate measures of variavility for exp.
    exp.variability = CalculateVariability(as.vector(exp.r.p.ref), 
                                           as.vector(exp.r.p.2), 
                                           aligned.p.1.index, 
                                           aligned.p.2.index, 
                                           not.aligned.p.1.index,
                                           not.aligned.p.2.index,
                                           R0,
                                           rotate,
                                           TOLERANCE)
    
    m.exp.va[P, 1:length(exp.variability$va)]  = exp.variability$va
    m.exp.Pn[P, 1:length(exp.variability$Pn)] = exp.variability$Pn
    exp.dr.squarei = exp.variability$dr.squarei
    for (i in (1:n.sites.p.ref)) {
      m.exp.dr.squarei[P, i] = matrix(exp.dr.squarei[ aligned.p.1.index == i], nrow = 1, ncol = 1)
    }

    for (mut in (1:n.mut.p)) {
    print(c(P,mut))
      
      # Get theo.r.p.2.
      theo.r.p.2 = m.r.mut[, n.mut.p * P - (n.mut.p - mut)]
    
      # Calculate measures of variavility for theo.
      theo.variability = CalculateVariability(as.vector(theo.r.p.ref), 
                                              as.vector(theo.r.p.2), 
                                              aligned.p.1.index, 
                                              aligned.p.1.index, 
                                              not.aligned.p.1.index,
                                              not.aligned.p.1.index,
                                              R0,
                                              rotate,
                                              TOLERANCE)
      
      m.theo.va[n.mut.p * P - (n.mut.p - mut), 1:length(theo.variability$va)]  = theo.variability$va
      m.theo.Pn[n.mut.p * P - (n.mut.p - mut), 1:length(theo.variability$Pn)] = theo.variability$Pn
      theo.dr.squarei = theo.variability$dr.squarei
      for (i in (1:n.sites.p.ref)) {
        m.theo.dr.squarei[n.mut.p * P - (n.mut.p - mut), i] = matrix(theo.dr.squarei[ aligned.p.1.index == i], nrow = 1, ncol = 1)
      }
    }
  }
  
  # Create files to save the data.
  write.csv(m.exp.va, file = file.path(out.dir, paste(out.name, "_out_m.exp.va.csv", sep = "")), row.names = FALSE)
  write.csv(m.exp.Pn, file = file.path(out.dir, paste(out.name, "_out_m.exp.Pn.csv", sep = "")), row.names = FALSE)
  write.csv(m.exp.dr.squarei, file = file.path(out.dir, paste(out.name, "_out_m.exp.dr.squarei.csv", sep = "")), row.names = FALSE)
  write.csv(m.theo.va, file = file.path(out.dir, paste(out.name, "_out_m.theo.va.csv", sep = "")), row.names = FALSE)
  write.csv(m.theo.Pn, file = file.path(out.dir, paste(out.name, "_out_m.theo.Pn.csv", sep = "")), row.names = FALSE)
  write.csv(m.theo.dr.squarei, file = file.path(out.dir, paste(out.name, "_out_m.theo.dr.squarei.csv", sep = "")), row.names = FALSE)
}
