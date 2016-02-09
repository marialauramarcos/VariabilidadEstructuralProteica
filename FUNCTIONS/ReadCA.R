# This function reads a pdb file and returns CA equilibrium coordinates, sites and nsites.
#
#  Args:
#    pdb.fname: pdb file name.
#    chain: chain to read.
#
#  Returns:
#    xyz.calpha: Equilibrium coordinates of CAs.
#    site: A vector with the sites.
#    nsites: length(site).
ReadCA <- function(pdb.fname, chain) {
  pdb <- read.pdb(file = pdb.fname)     
  sel <- atom.select(pdb, chain = chain, elety = "CA")    
  site <- as.numeric(pdb$atom[sel$atom, c("resno")])    
  n.sites <- length(site)    
  xyz.calpha <- matrix(pdb$xyz[sel$xyz], ncol = n.sites, nrow = 3, byrow = F)    
  output <- list("xyz.calpha" = xyz.calpha, "site" = site, "n.sites" = n.sites)
  output
}
