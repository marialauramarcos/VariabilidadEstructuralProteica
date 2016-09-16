# Description:
#
# This function reads a pdb file and returns the coordinates of the center of mass of side chains of each aminoacid.
#
#  Args:
#   - pdb.fname: pdb filename.
#   - chain: chain to read.
#
#  Returns:
#    xyz.side.chain.CM: coordinates of the centers of mass of side chains of each aminoacid.

CalculateSideChainCM <- function(pdb.fname, chain) {
  
  # Read pdb
  pdb <- read.pdb(file = pdb.fname)
  
  # Get n.aa of the chain
  sel.ca <- atom.select(pdb, chain = chain, elety = "CA")    
  site.ca = as.numeric(pdb$atom[sel.ca$atom, c("resno")])    
  n.aa = length(site.ca)    

  # Create a matrix to the save the centers of mass
  xyz.side.chain.CM = matrix(nrow = 3, ncol = n.aa)
  
  # Star a loop to analyze each residue
  for (i in (site.ca)) { 
    
    # Extract information of the chain
    sel <- atom.select(pdb, chain = chain, strain = "protein", resno = i)    
    site.elety <- pdb$atom[sel$atom, c("elety")]     
    n.atoms = length(site.elety) 
    xyz.site <- matrix(pdb$xyz[sel$xyz], nrow = 3, byrow = F)     
  
    # Get weights
    weights = matrix(nrow = n.atoms)

    index.C = which(grepl("C", site.elety))
    index.O = which(grepl("O", site.elety))
    index.N = which(grepl("N", site.elety))
    index.S = which(grepl("S", site.elety))
    
    weights[index.C, ] = 12
    weights[index.O, ] = 16
    weights[index.N, ] = 14
    weights[index.S, ] = 32
    
    # Eliminate the weights of CAs
    for(j in (1:n.atoms)) {
      if (site.elety[j] == "CA") weights[j, ] = 0
    }
    
    # Calculate the center of mass
    xyz.side.chain.CM[1, which(site.ca == i)] = sum(xyz.site[1, ] * weights) / sum(weights)
    xyz.side.chain.CM[2, which(site.ca == i)] = sum(xyz.site[2, ] * weights) / sum(weights)
    xyz.side.chain.CM[3, which(site.ca == i)] = sum(xyz.site[3, ] * weights) / sum(weights)
  }
  xyz.side.chain.CM
}

  

