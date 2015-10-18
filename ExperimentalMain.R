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
# Remove objects.
rm(list = ls())
# Read input.
data.dir <- "DATA/Experimental"
input.fname <- file.path(data.dir, "experimental_input.csv")
input <- read.csv(input.fname)
family <- input$family
p.ref <- input$p.ref 
heme <- input$heme
model <- input$model
R0 = input$R0
# Files fnames.
alignment.fname <- file.path(data.dir, paste(family, "_alignment.txt", sep = ""))  # File with the multiple sequence alignmet.
pdbs.fname <- file.path(data.dir, paste(family, "_coordinates.pdb", sep = ""))  # Coordinates of the proteins.
dataset.fname <- file.path(data.dir, paste(family, "_dataset.csv", sep = ""))  # Dataset with pdbids and chains.
# Functions fnames.
ReadCA.fname <- "FUNCTIONS/ReadCA.R" 
ReadHeme.fname <- "FUNCTIONS/ReadHeme.R"  
AnalyzeAlignment.fname <- "FUNCTIONS/AnalyzeAlignment.R" 
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
source(AnalyzeAlignment.fname)
source(CalculateID.fname)
source(CalculateKeff.fname)
source(CalculateK.fname)
source(CalculateKij.fname)
source(CalculateVariability.fname)
# Read dataset.
dataset <- read.csv(dataset.fname)
pdbids <- dataset$pdbid
chains <- dataset$chain
nprot = length(pdbids)
# Read multiple alignment.
alignments.ids <- read.fasta(alignment.fname)
alignments <- alignments.ids$ali[, -ncol(alignments.ids$ali)]
# Read PDB & alignment of p.ref.
chain.p.ref <- chains[pdbids == as.character(p.ref)]
pdb.p.ref <- ReadCA(pdbs.fname, chain.p.ref)
r.p.ref = pdb.p.ref$xyz.calpha
nsites.p.ref = pdb.p.ref$nsites 
naa.p.ref = nsites.p.ref 
if (family == "globins" & heme == "TRUE") {
  r.heme.p.ref = ReadHeme(pdbs.fname, chain.p.ref)
  r.p.ref = cbind(r.p.ref, r.heme.p.ref)
  nsites.p.ref = ncol(r.p.ref)
}
alignment.p.ref <- alignments[alignments.ids$id == as.character(p.ref), ]
# Measures of variability.
m.laligned = matrix(ncol = 1, nrow = nprot)
m.ID = matrix(ncol = 1, nrow = nprot)
m.nH = matrix(ncol = 3 * naa.p.ref, nrow = nprot)
m.Pn = matrix(ncol = 3 * naa.p.ref, nrow = nprot)
m.evalues = matrix(ncol = 3 * naa.p.ref, nrow = nprot)
m.d.evalues = matrix(ncol = 3 * naa.p.ref, nrow = nprot)
m.dr.squarei = matrix(ncol = naa.p.ref, nrow = nprot)
for (P in (1:nprot)) {
	# Read PDB & alignment of p.2.
	chain.p.2 <- chains[[P]]
	pdb.p.2 <- ReadCA(pdbs.fname, chain.p.2)
	r.p.2 <- pdb.p.2$xyz.calpha
	nsites.p.2 <- pdb.p.2$nsites
	naa.p.2 = nsites.p.2
	if (family == "globins" & heme == "TRUE") {
	  r.heme.p.2 = ReadHeme(pdbs.fname, chain.p.2)
	  r.p.2 = cbind(r.p.2, r.heme.p.2)
	  nsites.p.2 = ncol(r.p.2)
	}
	alignment.p.2 <- alignments[alignments.ids$id == as.character(pdbids[P]), ]	
	# Calculate %ID between P.ref y p.2.
	ID.p.2 <- CalculateID(alignment.p.ref, alignment.p.2)
	# Analyze alignment of p.ref.
	a.alignment.p.ref <- AnalyzeAlignment(alignment.p.ref, alignment.p.2, naa.p.ref)
	aligned.p.ref.index <- a.alignment.p.ref$aligned.index
	naligned <- a.alignment.p.ref$naligned
	not.aligned.p.ref.index <- a.alignment.p.ref$not.aligned.index
	# Analyze alignment of p.2.
	a.alignment.p.2 <- AnalyzeAlignment(alignment.p.2, alignment.p.ref, naa.p.2)
	aligned.p.2.index <- a.alignment.p.2$aligned.index
	not.aligned.p.2.index <- a.alignment.p.2$not.aligned.index
	# Add heme to not.aligned.
	if (family == "globins" & heme == "TRUE"){
	  not.aligned.p.ref.index <- cbind(not.aligned.p.ref.index, t(seq((naa.p.ref+1), nsites.p.ref)))
	  not.aligned.p.2.index <- cbind(not.aligned.p.2.index, t(seq((naa.p.2+1), nsites.p.2)))
	}
	# Rotate r.p.2 minimizing RMSD with P.ref.
	aligned.p.ref.index3N = matrix(0, ncol = 3 * naligned, nrow = 1)
	aligned.p.2.index3N = matrix(0, ncol = 3 * naligned, nrow = 1)
  for (i in (1:naligned)) {
    aligned.p.ref.index3N[1, ((3 * i - 2):(3 * i))] = c((3 * aligned.p.ref.index[i] - 2), (3 * aligned.p.ref.index[i] - 1),(3 * aligned.p.ref.index[i]))
    aligned.p.2.index3N[1, ((3 * i - 2):(3 * i))] = c((3 * aligned.p.2.index[i] - 2), (3 * aligned.p.2.index[i] - 1),(3 * aligned.p.2.index[i]))
  }
	r.p.2 <- matrix(fit.xyz(fixed = as.vector(r.p.ref),
	                       mobile = as.vector(r.p.2),
	                   fixed.inds = aligned.p.ref.index3N,
                    mobile.inds = aligned.p.2.index3N), nrow = 3)
	# Get aligned coordinates.
	r.aligned.p.ref <- r.p.ref[, aligned.p.ref.index]
	r.aligned.p.2 <- r.p.2[, aligned.p.2.index]
	# Calculate dr.
	dr = r.aligned.p.2 - r.aligned.p.ref
	# Cakculate KEFF p.ref.
	Keff.p.ref <- CalculateKeff(r.p.ref, aligned.p.ref.index, not.aligned.p.ref.index, CalculateKij, R0, TOLERANCE)	
	nmodes <- length(Keff.p.ref$va)
	# Calculate KEFF p.2.
	Keff.p.2 <- CalculateKeff(r.p.2, aligned.p.2.index, not.aligned.p.2.index, CalculateKij, R0, TOLERANCE)	
	# Calculate measures of variability.
	VA <- CalculateVariability(dr, Keff.p.ref, Keff.p.2)
	m.ID[P] = ID.p.2
  m.laligned[P] = naligned 
  m.nH[P, 1:nmodes] = t(VA$nH)
  m.Pn[P, 1:nmodes] = t(VA$Pn)
  m.d.evalues[P, 1:nmodes]  = t(VA$d.evalues[1:nmodes])
	m.evalues[P, 1:nmodes]  = t(Keff.p.ref$va[1:nmodes])
	dr.squarei = rbind(VA$dr.squarei, aligned.p.ref.index)
	for (i in (1:naa.p.ref)){
		m.dr.squarei[P, i] = matrix(dr.squarei[1, dr.squarei[2, ] == i], ncol = 1, nrow = 1)
	}
}
P.index = (pdbids != as.character(p.ref))
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
write.csv(m.laligned, file = file.path(out.dir, paste(family, "_out_m.laligned.csv", sep = "")), row.names = FALSE)
write.csv(m.ID, file = file.path(out.dir, paste(family, "_out_m.ID.csv", sep = "")), row.names = FALSE)
write.csv(m.nH, file = file.path(out.dir, paste(family, "_out_m.nH.csv", sep = "")), row.names = FALSE)
write.csv(m.Pn, file = file.path(out.dir, paste(family, "_out_m.Pn.csv", sep = "")), row.names = FALSE)
write.csv(m.d.evalues, file = file.path(out.dir, paste(family, "_out_m.d.evalues.csv", sep = "")), row.names = FALSE)
write.csv(m.evalues, file = file.path(out.dir, paste(family, "_out_m.evalues.csv", sep = "")), row.names = FALSE)
write.csv(m.dr.squarei, file = file.path(out.dir, paste(family, "_out_m.dr.squarei.csv", sep = "")), row.names = FALSE)
write.csv(mean.nH, file = file.path(out.dir, paste(family, "_out_nH.mean.csv", sep = "")), row.names = FALSE)
write.csv(mean.Pn, file = file.path(out.dir, paste(family, "_out_Pn.mean.csv", sep = "")), row.names = FALSE)
write.csv(mean.d.evalues, file = file.path(out.dir, paste(family, "_out_d.evalues.mean.csv", sep = "")), row.names = FALSE)
write.csv(mean.evalues, file = file.path(out.dir, paste(family, "_out_evalues.mean.csv", sep = "")), row.names = FALSE)
write.csv(MSDi, file = file.path(out.dir, paste(family, "_out_MSDi.csv", sep = "")), row.names = FALSE)
write.csv(MSD, file = file.path(out.dir, paste(family , "_out_MSD.csv", sep = "")), row.names = FALSE)
write.csv(input, file = file.path(out.dir, paste(family , "_input", sep = "")), row.names = FALSE)
