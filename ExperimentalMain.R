# INFORMACIÓN:
#
# El programa calcula medidas de variabilidad estructural y dinámica de alineamientos múltiples de 
# familias de proteínas provenientes de la base de datos de alineamientos múltiples "Homstrad".
#
# Para utilizar el programa se debe completar un input ("DATA/Experimental/experimental_input.csv") especificando: 
# -family:la familia de proteínas del alineamiento, que puede ser "serinProteases", "globins" o "plastocyanins"
# -p.ref: la proteína de referncia que se desee utilizar, teniendo en cuenta que la misma debe ser escrita de 
# la misma forma en la que se encuentra en el dataset de la familia correspondiente 
# (ej:"DATA/Experimental/serinProteases_dataset.csv").
# -heme: solo se utiliza para la familia de las globinas, puede ser "TRUE" o "FALSE", dependiendo de si 
# se desea o no considerar al grupo HEMO. 
# -model: el modelo de red elástica a usar, que solo puede ser "ANM" por el momento.
# -R0: cut-off del ANM.
# -core: puede ser "TRUE" o "FALSE" dependiendo de si se quiere analizar solo los sectores
# del alineamiento donde no hay gaps. 

# Remove objects.
rm(list = ls())

# Data.dir.
data.dir <- "DATA/Experimental"

# Read input.
input.fname <- file.path(data.dir, "experimental_input.csv")
input <- read.csv(input.fname)
family <- as.character(input$family) 
p.ref <- as.character(input$p.ref) 
heme <- as.character(input$heme)  
model <- input$model  
R0 = input$R0  
core = input$core  

# Data fnames.
alignment.fname <- file.path(data.dir, paste(family, "_alignment.txt", sep = ""))  
pdbs.fname <- file.path(data.dir, paste(family, "_coordinates.pdb", sep = "")) 
dataset.fname <- file.path(data.dir, paste(family, "_dataset.csv", sep = ""))  

# Functions fnames.
ReadCA.fname <- "FUNCTIONS/ReadCA.R" 
ReadHeme.fname <- "FUNCTIONS/ReadHeme.R"  
AnalyzeAlignmentGeneral.fname <- "FUNCTIONS/AnalyzeAlignmentGeneral.R" 
CalculateID.fname <- "FUNCTIONS/CalculateID.R" 
CalculateKeff.fname <- "FUNCTIONS/CalculateKeff.R"  
CalculateK.fname <- "FUNCTIONS/CalculateK.R"  
CalculateKij.fname <- file.path("FUNCTIONS", paste("CalculateKij", model, ".R", sep = "")) 
CalculateVariability.fname <- "FUNCTIONS/CalculateVariability.R" 

# Output dir.
out.dir <- "OUT/Experimental"

# General parameters.
TOLERANCE = 1e-10 

# Load Librarys.
library(seqinr) 
library(bio3d) 

# Source Functions.
source(ReadCA.fname) 
source(ReadHeme.fname) 
source(AnalyzeAlignmentGeneral.fname)
source(CalculateID.fname)
source(CalculateKij.fname)
source(CalculateK.fname)
source(CalculateKeff.fname)
source(CalculateVariability.fname)

# Read dataset.
dataset <- read.csv(dataset.fname)
pdbid <- as.character(dataset$pdbid) 
chain <- as.character(dataset$chain)
nprot = length(pdbid)

# Read PDB of p.ref.
chain.p.ref <- chain[pdbid == p.ref]
pdb.p.ref <- ReadCA(pdbs.fname, chain.p.ref)
r.p.ref = pdb.p.ref$xyz.calpha
naa.p.ref = pdb.p.ref$nsites

# Read alignment.
alignment.id <- read.fasta(alignment.fname)
alignment <- alignment.id$ali[, -ncol(alignment.id$ali)]  # The last column is "*".
id <- alignment.id$id

# Measures of variability.
m.laligned = matrix(ncol = 1, nrow = nprot)
m.ID = matrix(ncol = 1, nrow = nprot)
m.nH = matrix(ncol = 3 * naa.p.ref, nrow = nprot)
m.Pn = matrix(ncol = 3 * naa.p.ref, nrow = nprot)
m.evalues = matrix(ncol = 3 * naa.p.ref, nrow = nprot)
m.d.evalues = matrix(ncol = 3 * naa.p.ref, nrow = nprot)
m.dr.squarei = matrix(ncol = naa.p.ref, nrow = nprot)

for (P in (1:nprot)) {
  
  # Create a data frame to analyze the alignment.
	df.alignment <- list("alignment" = alignment,
	                            "id" = id,   
	                           "p.1" = p.ref,
	                           "p.2" = pdbid[P])
	if (core == "TRUE"){
	  class(df.alignment) <- "Core"
	} else {
	  class(df.alignment) <- "NoCore"
	}
	
	# Anylize alignment.
	a.alignment <- AnalyzeAlignmentGeneral(df.alignment)
	aligned.p.ref.index <- a.alignment$aligned.p.1.index
	aligned.p.2.index <- a.alignment$aligned.p.2.index
	not.aligned.p.ref.index <- a.alignment$not.aligned.p.1.index
	not.aligned.p.2.index <- a.alignment$not.aligned.p.2.index
	naligned <- a.alignment$naligned
  
  # Read PDB of p.2.  
  chain.p.2 <- chain[[P]]
  pdb.p.2 <- ReadCA(pdbs.fname, chain.p.2)
  r.p.2 <- pdb.p.2$xyz.calpha
  naa.p.2 <- pdb.p.2$nsites
  
  # Calculate r.heme and add to r and to not.aligned.index.
	if (family == "globins" & heme == "TRUE") {
	  if (P == 1) {
	    r.heme.p.ref = ReadHeme(pdbs.fname, chain.p.ref)
	    r.p.ref = cbind(r.p.ref, r.heme.p.ref)
	    nsites.p.ref = ncol(r.p.ref)
	  }
	  not.aligned.p.ref.index <- c(not.aligned.p.ref.index, t(seq((naa.p.ref+1), nsites.p.ref)))
	  
	  r.heme.p.2 = ReadHeme(pdbs.fname, chain.p.2)
	  r.p.2 = cbind(r.p.2, r.heme.p.2)
	  nsites.p.2 = ncol(r.p.2)
	  not.aligned.p.2.index <- c(not.aligned.p.2.index, t(seq((naa.p.2+1), nsites.p.2)))
	} else {
	  nsites.p.ref = naa.p.ref
	  nsites.p.2 = naa.p.2
	}
	
	# Rotate r.p.2 minimizing RMSD with P.ref.
	aligned.p.ref.index3N = c(aligned.p.ref.index * 3,
	                          aligned.p.ref.index * 3 - 2,
	                          aligned.p.ref.index * 3 - 1)
	aligned.p.2.index3N = c(aligned.p.2.index * 3,
	                        aligned.p.2.index * 3 - 2,
	                        aligned.p.2.index * 3 - 1)
	
	r.p.2 <- matrix(fit.xyz(fixed = as.vector(r.p.ref),
	                       mobile = as.vector(r.p.2),
	                   fixed.inds = aligned.p.ref.index3N,
                    mobile.inds = aligned.p.2.index3N), nrow = 3)
	
	# Calculate dr.
	dr = r.p.2[, aligned.p.2.index] - r.p.ref[, aligned.p.ref.index]
	
	# Cakculate K of p.ref and p.2.
	K.p.ref <- CalculateKeff(r.p.ref, 
	                         aligned.p.ref.index, 
	                         not.aligned.p.ref.index,
	                         CalculateK,
	                         CalculateKij, 
	                         R0, 
	                         TOLERANCE)	
	K.p.2 <- CalculateKeff(r.p.2, 
	                       aligned.p.2.index,
	                       not.aligned.p.2.index, 
	                       CalculateK,
	                       CalculateKij, 
	                       R0, 
	                       TOLERANCE)	
	
	# Calculate nmodes.
	nmodes <- length(K.p.ref$va)
	
	# Calculate % sequence identity between p.ref and p.2.
	ID.p.2 <- CalculateID(df.alignment)
	
	# Calculate measures of variability.
	VA <- CalculateVariability(dr, K.p.ref, K.p.2)
	m.ID[P] = ID.p.2
  m.laligned[P] = naligned 
  m.nH[P, 1:nmodes] = t(VA$nH)
  m.Pn[P, 1:nmodes] = t(VA$Pn)
  m.d.evalues[P, 1:nmodes]  = t(VA$d.evalues[1:nmodes])
	m.evalues[P, 1:nmodes]  = t(K.p.ref$va[1:nmodes])
	dr.squarei = rbind(VA$dr.squarei, aligned.p.ref.index)
	for (i in (1:naa.p.ref)){
		m.dr.squarei[P, i] = matrix(dr.squarei[1, dr.squarei[2, ] == i], ncol = 1, nrow = 1)
	}
}

P.index = (pdbid != p.ref)
m.ID <- m.ID[P.index, ]
m.laligned <- m.laligned[P.index, ]
m.nH <- m.nH[P.index, ]
m.Pn <- m.Pn[P.index, ]
m.d.evalues <- m.d.evalues[P.index, ]
m.evalues <- m.evalues[P.index, ]
m.dr.squarei <- m.dr.squarei[P.index, ]

# Calculate means.
mean.nH = colMeans(m.nH, na.rm = T)
mean.Pn = colMeans(m.Pn, na.rm = T)
mean.d.evalues= colMeans(m.d.evalues, na.rm = T)
mean.evalues= colMeans(m.evalues, na.rm = T)
MSDi = colMeans(m.dr.squarei, na.rm = T)
MSD = rowMeans(m.dr.squarei, na.rm = T)

# Save information.
if (family == "globins") {
  family <- paste(family, "_heme_", heme, sep = "")
}
write.csv(m.laligned, file = file.path(out.dir, paste(family, "_core_", core, "_out_m.laligned.csv", sep = "")), row.names = FALSE)
write.csv(m.ID, file = file.path(out.dir, paste(family, "_core_", core, "_out_m.ID.csv", sep = "")), row.names = FALSE)
write.csv(m.nH, file = file.path(out.dir, paste(family, "_core_", core, "_out_m.nH.csv", sep = "")), row.names = FALSE)
write.csv(m.Pn, file = file.path(out.dir, paste(family, "_core_", core, "_out_m.Pn.csv", sep = "")), row.names = FALSE)
write.csv(m.d.evalues, file = file.path(out.dir, paste(family, "_core_", core, "_out_m.d.evalues.csv", sep = "")), row.names = FALSE)
write.csv(m.evalues, file = file.path(out.dir, paste(family, "_core_", core, "_out_m.evalues.csv", sep = "")), row.names = FALSE)
write.csv(m.dr.squarei, file = file.path(out.dir, paste(family, "_core_", core, "_out_m.dr.squarei.csv", sep = "")), row.names = FALSE)

write.csv(mean.nH, file = file.path(out.dir, paste(family, "_core_", core, "_out_nH.mean.csv", sep = "")), row.names = FALSE)
write.csv(mean.Pn, file = file.path(out.dir, paste(family, "_core_", core, "_out_Pn.mean.csv", sep = "")), row.names = FALSE)
write.csv(mean.d.evalues, file = file.path(out.dir, paste(family, "_core_", core, "_out_d.evalues.mean.csv", sep = "")), row.names = FALSE)
write.csv(mean.evalues, file = file.path(out.dir, paste(family, "_core_", core, "_out_evalues.mean.csv", sep = "")), row.names = FALSE)
write.csv(MSDi, file = file.path(out.dir, paste(family, "_core_", core, "_out_MSDi.csv", sep = "")), row.names = FALSE)
write.csv(MSD, file = file.path(out.dir, paste(family , "_core_", core, "_out_MSD.csv", sep = "")), row.names = FALSE)
write.csv(input, file = file.path(out.dir, paste(family , "_input", sep = "")), row.names = FALSE)
