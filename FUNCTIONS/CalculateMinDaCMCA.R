# Description:
#
# This function calculates distances between functional sites of proteins and the other sites of the protein.
# It calculates distances between CAs and CMs.
#
# Usage:
#
# CalculateMinDaCMCA(family,
#                    p.ref,
#                    chain.p.ref,
#                    heme,
#                    data.dir)
#
#  Args:
#    - family: The family of the protein to mutate. It can be "globins", "serinProteases", 
#    "snakesToxin", "sh3", "fabp", "rrm", "phoslip" or "cys".
#    - p.ref: The pdb code (pdbid) of the protein to mutate (example: "1a6m"). The protein must be a member of
#    the selected family. This pdbid must not be included in the dataset ("DATA/family_dataset.csv").
#    - chain.p.ref: The chain of p.ref in the pdb file obtained from Homstrad.
#    - heme: argument for "globins". It can be "TRUE" or "FALSE". If it is "TRUE", the program considers the heme group.
#    - data.dir: Directory of the data. It must contain the pdb file obtained from Homstrad ("data.dir/family_coordinates.csv").
#
#  Required libraries:
#    {Bio3d}
#
#  Required functions:
#    ReadCA()
#    ReadHeme()
#    CalculateSideChainCM()

  CalculateMinDaCMCA <- function(family,
                                 p.ref,
                                 chain.p.ref,
                                 heme,
                                 data.dir) {
  
  # read functional sites
  active = read.csv(file.path(data.dir, paste(p.ref, "_functionalSites.csv", sep = "")), sep = ";")$functional.sites
  active = active[!is.na(active)]
  
  # pdb filename
  pdb.fname <- file.path(data.dir, paste(p.ref, ".pdb", sep = ""))

  # get coordinates of p.ref
  r.ca.p.ref = ReadCA(pdb.fname, chain.p.ref)$xyz
  r.CM.p.ref = CalculateSideChainCM(pdb.fname, chain.p.ref)
  n.aa = ncol(r.ca.p.ref)
  
  if (heme == "TRUE") {
    r.heme = ReadHeme(pdb.fname, chain.p.ref)
    r.ca.p.ref = cbind(r.ca.p.ref, r.heme)
    r.CM.p.ref = cbind(r.CM.p.ref, r.heme)
  }
  
  # calculate the number of sites
  n.sites = ncol(r.ca.p.ref)

  # calculate the distances between each site and each functional site (CA or CM)
  m.da.ca = matrix(nrow = n.sites, ncol = length(active))
  m.da.CM = matrix(nrow = n.sites, ncol = length(active))
  
  for (i in (active)) {
    print(i)
    for (j in (1:n.sites)) {
      m.da.ca[j, which(active == i)] = sqrt(sum((r.ca.p.ref[, j] - r.ca.p.ref[, i]) ^ 2))
      m.da.CM[j, which(active == i)] = sqrt(sum((r.CM.p.ref[, j] - r.CM.p.ref[, i]) ^ 2))
    }
  }
  
  m.da.ca = cbind(seq(1:n.sites), m.da.ca)
  colnames(m.da.ca) = c("site", paste("da", seq(1:length(active)), sep = ""))
  
  m.da.CM = cbind(seq(1:n.sites), m.da.CM)
  colnames(m.da.CM) = c("site", paste("da", seq(1:length(active)), sep = ""))
  
  if (family == "globins") {
    m.da.ca = m.da.ca[1:n.aa, ]
    m.da.CM = m.da.CM[1:n.aa, ]
  }
  
  # calculate the minimum distance of each site to the active sites
  m.min.da.ca = matrix(nrow = n.aa, ncol = 1)
  m.min.da.CM = matrix(nrow = n.aa, ncol = 1)
  
  for (i in (1:n.aa)) {
    m.min.da.ca[i, ] = min(m.da.ca[i, ])
    m.min.da.CM[i, ] = min(m.da.CM[i, ])
  }

  # build a dataframe with minimum distances
  data = data.frame("site" = seq(1:n.aa), 
               "min.da.CM" = as.vector(m.min.da.CM), 
               "min.da.ca" = as.vector(m.min.da.ca))

  # write output files
  out.dir = "OUT/out_subset_CA_ANM"
  write.csv(data, file.path(out.dir, paste(p.ref, "_min.da.CM.ca.csv", sep = "")))
  write.csv(m.da.ca, file.path(out.dir, paste(p.ref, "_m.da.ca.csv", sep = "")))
  write.csv(m.da.CM, file.path(out.dir, paste(p.ref, "_m.da.CM.csv", sep = "")))

  out.dir = "OUT/out_subset_CM_ANM" 
  write.csv(data, file.path(out.dir, paste(p.ref, "_min.da.CM.ca.csv", sep = "")))
  write.csv(m.da.ca, file.path(out.dir, paste(p.ref, "_m.da.ca.csv", sep = "")))
  write.csv(m.da.CM, file.path(out.dir, paste(p.ref, "_m.da.CM.csv", sep = "")))
}
  