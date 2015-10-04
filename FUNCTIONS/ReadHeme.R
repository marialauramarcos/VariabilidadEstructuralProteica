#Function that reads a pdb file and returns coordinates of NA, NB, NC, ND and Fe of the Heme group
#A pdb file name and the chain must be specified

readHeme <- function(pdb.fname,chain){
	pdb <- read.pdb(file=pdb.fname,het2atom = TRUE)     
	selNA <- atom.select(pdb,chain=chain,elety="NA") 
	selNB <- atom.select(pdb,chain=chain,elety="NB")    
	selNC <- atom.select(pdb,chain=chain,elety="NC")    
	selND <- atom.select(pdb,chain=chain,elety="ND")    
	selFE <- atom.select(pdb,chain=chain,elety="FE")    

	xyz.heme <- matrix(c(pdb$xyz[selNA$xyz],pdb$xyz[selNB$xyz],pdb$xyz[selNC$xyz],pdb$xyz[selND$xyz],pdb$xyz[selFE$xyz]),nrow=3)
	xyz.heme
}
