# INFORMACIÓN:
#
# El programa genera mutantes múltiples de una proteína determinada usando el modelo mutacional Linearly 
# Forced - Elastic Network Model (LF-ENM) y considerando a las fuerzas que modelan cada mutación aditivas 
# entre sí. Además, calcula medidas de variabilidad estructural y dinámica entre las proteínas mutantes 
# generadas y la proteína de referencia.
#
# Para utilizar el programa se debe completar un input ("DATA/Theoretical/theoretical_input.csv") especificando: 
# -family: la familia de proteínas a la que pertenece la proteína a mutar. 
# Tener en cuenta que las familias que pueden ser analizadas experimentalmente son
# "serinProteases", "globins" o "plastocyanins".
# -p.ref: el código de pdb (pdbid) de la proteína a mutar (ej.:"1a6m")
# -chain: la cadena de p.ref que se desea mutar.
# -heme: solo se utiliza para la familia de las globinas, puede ser "TRUE" o "FALSE" dependiendo de 
# si se quiere considerar o no al grupo HEMO. 
# -nmut: la cantidad de mutantes a generar.
# -nsitesmut: la cantidad de sitios mutados por proteína mutante.
# -Fmax: valor máximo de las fuerzas usadas para la simulación de las mutaciones.
# -model: el modelo de red elástica a usar, que solo puede ser "ANM" por el momento.
# -R0: cut-off del ANM.
# Remove objects.
rm(list = ls()) 
# Read input.
data.dir <- "DATA/Theoretical"
input.fname <- file.path(data.dir, "theoretical_input.csv")
input <- read.csv(input.fname)
family <- input$family
p.ref <- input$p.ref 
chain <- input$chain
heme <- input$heme 
nmut = input$nmut 
nsitesmut = input$nsitesmut 
fmax = input$fmax 
model <- input$model
R0 = input$R0
# PDB fname.
get.pdb(as.character(p.ref), data.dir) 
pdb.fname <- file.path(data.dir, paste(p.ref, ".pdb", sep = ""))
# Functions fnames.
ReadCA.fname <- "FUNCTIONS/ReadCA.R" 
ReadHeme.fname <- "FUNCTIONS/ReadHeme.R"  
CalculateKeff.fname <- "FUNCTIONS/CalculateKeff.R"  
CalculateK.fname <- "FUNCTIONS/CalculateK.R" 
CalculateKij.fname <- file.path("FUNCTIONS", paste("CalculateKij", model, ".R", sep = "")) 
CalculateForce.fname <- "FUNCTIONS/CalculateForce.R"  
CalculateVariability.fname <- "FUNCTIONS/CalculateVariability.R"
# Output dir.
out.dir <- "OUT/Theoretical"
# General parameters.
TOLERANCE = 1e-10
# Load librarys.
library(bio3d)
# Source functions.
source(ReadCA.fname)
source(ReadHeme.fname)
source(CalculateKeff.fname)
source(CalculateK.fname)
source(CalculateKij.fname)
source(CalculateForce.fname)
source(CalculateVariability.fname)
# Read PDB of p.ref.
pdb <- ReadCA(pdb.fname, chain)
r.p.ref = pdb$xyz.calpha
nsites = pdb$nsites
naa = nsites  # naa = n.aminoaids, nsites = n.aminoácids + n.atoms we want to consider, ej.: Ns and Fe from heme group in Globins# 
if (family == "globins" & heme == "TRUE") {
  r.heme = ReadHeme(pdb.fname, chain)
  r.p.ref = cbind(r.p.ref, r.heme)
  nsites = ncol(r.p.ref)
  aligned.index = seq(1, naa)
  not.aligned.index = seq((naa+1), nsites)
}
# Calculate K of p.ref.
K.p.ref <- CalculateK(r.p.ref, CalculateKij, R0, TOLERANCE)
# Measures of variability.
m.nH = matrix(ncol = 3 * naa, nrow = nmut)
m.Pn = matrix(ncol = 3 * naa, nrow = nmut)
m.d.evalues = matrix(ncol = 3 * naa, nrow = nmut)
m.dr.squarei = matrix(ncol = naa, nrow = nmut)
# Count the number of mutants to discard.
count = 0
# Calculate mutants.
for(mut in seq(nmut)) {
  f <- rep(0, 3 * nsites)
  for (l in (sample(1:naa, replace = F)[1:nsitesmut])) {
    print(c(mut, l))
    fl = CalculateForce(l, r.p.ref, K.p.ref$kij, fmax)
    f = f + fl
  }
  dr = K.p.ref$cov %*% f
  dim(dr) <- c(3, nsites)
  r.mut = r.p.ref + dr
  # Rotate mutants minimizing RMSD with P.ref.
  r.mut <- matrix(fit.xyz(fixed = as.vector(r.p.ref), 
                         mobile = as.vector(r.mut),
                     fixed.inds = seq(1, 3 * naa),
                    mobile.inds = seq(1, 3 * naa)), nrow = 3) 
  # Calculate K of mutants and p.ref.
  if (naa != nsites) {
    K.mut <- CalculateKeff(r.mut, aligned.index, not.aligned.index, CalculateKij, R0, TOLERANCE)
    K.p.ref.2 <- CalculateKeff(r.p.ref, aligned.index, not.aligned.index, CalculateKij, R0, TOLERANCE)
  }
  if (naa == nsites) {
    K.mut <- CalculateK(r.mut, CalculateKij, R0, TOLERANCE)
    K.p.ref.2 <- K.p.ref
  }
  nmodes.p.ref.2 <- length(K.p.ref.2$va)  
  nmodes.mut <- length(K.mut$va)  
  # Calculate dr.
  dr = (r.mut - r.p.ref)[, 1:naa]   
  # Calculate measures o variability.
  if (nmodes.p.ref.2 == nmodes.mut) { 
    VA <- CalculateVariability(dr, K.p.ref.2, K.mut)
    m.nH[mut, 1:nmodes.p.ref.2] = t(VA$nH)
    m.Pn[mut, 1:nmodes.p.ref.2] = t(VA$Pn)
    m.d.evalues[mut, 1:nmodes.p.ref.2] = t(VA$d.evalues[1:nmodes.p.ref.2])
    dr.squarei = rbind(VA$dr.squarei, seq(1, naa))
    for (i in (1:naa)) {
      m.dr.squarei[mut, i] = matrix(dr.squarei[1, dr.squarei[2, ] == i], ncol = 1, nrow = 1)
    }
  }
  if (nmodes.p.ref.2 != nmodes.mut) { 
    count <- count +1
  }
}
# Calculate Means.
mean.nH = colMeans(m.nH, na.rm = T)
mean.Pn = colMeans(m.Pn, na.rm = T)
mean.d.evalues = colMeans(m.d.evalues, na.rm = T)
MSDi = colMeans(m.dr.squarei, na.rm = T)
MSD = rowMeans(m.dr.squarei, na.rm = T)
# Save information.
if (family == "globins") {
  family <- paste(family, "_heme_", heme, sep = "")
}
write.csv(m.nH, file = file.path(out.dir, paste(family, "_out_m.nH.csv", sep = "")), row.names = FALSE)
write.csv(m.Pn, file = file.path(out.dir, paste(family, "_out_m.Pn.csv", sep = "")), row.names = FALSE)
write.csv(m.d.evalues, file = file.path(out.dir, paste(family, "_out_m.d.evalues.csv", sep = "")), row.names = FALSE)
write.csv(m.dr.squarei, file = file.path(out.dir, paste(family, "_out_m.dr.squarei.csv", sep = "")), row.names = FALSE)
write.csv(K.p.ref.2$va, file = file.path(out.dir, paste(family, "_out_evalues.csv", sep = "")), row.names = FALSE)
write.csv(mean.nH, file = file.path(out.dir, paste(family, "_out_nH.mean.csv", sep = "")), row.names = FALSE)
write.csv(mean.Pn, file = file.path(out.dir, paste(family, "_out_Pn.mean.csv", sep = "")), row.names = FALSE)
write.csv(mean.d.evalues, file = file.path(out.dir, paste(family, "_out_d.evalues.mean.csv", sep = "")), row.names = FALSE)
write.csv(MSDi, file = file.path(out.dir, paste(family, "_out_MSDi.csv", sep = "")), row.names = FALSE)
write.csv(MSD, file = file.path(out.dir, paste(family, "_out_MSD.csv", sep = "")), row.names = FALSE)
write.csv(input, file = file.path(out.dir, paste(family, "_input", sep = "")), row.names = FALSE)
