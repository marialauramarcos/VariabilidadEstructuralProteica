#Function that reads a pdb file and returns CA coordinates, sites and nsites
#A pdb file name and the chain must be specified

readCA <- function(pdb.fname,chain) {
	pdb <- read.pdb(file=pdb.fname)     
	sel <- atom.select(pdb,chain=chain,elety="CA")    
	site <- as.numeric(pdb$atom[sel$atom,c("resno")])    
	nsites <- length(site)    
	xyz.calpha <- matrix(pdb$xyz[sel$xyz],ncol=nsites,nrow=3,byrow=F)    
	output <- list("xyz.calpha"=xyz.calpha,"site"=site,"nsites"=nsites)
  output
}
