# INFORMACIÓN:
#
# El programa genera mutantes múltiples de una proteína determinada usando el modelo mutacional Linearly 
# Forced - Elastic Network Model (LF-ENM) y considerando a las fuerzas que modelan cada mutación aditivas 
# entre sí. Además, calcula medidas de variabilidad estructural y dinámica entre las proteínas mutantes 
# generadas y la proteína de referencia.
#
# Para utilizar el programa se debe completar un input ("DATA/Theoretical/inputT.csv") especificando: 
# -family: la familia de proteínas a la que pertenece la proteína a mutar escrita en inglés y plural 
# (ej.:"Globins").
# -p.ref: el código de pdb (pdbid) de la proteína a mutar (ej.:"1a6m")
# -chain: la cadena de p.ref que se desea mutar.
# -heme: solo se utiliza para la familia de las globinas, puede ser "TRUE" o "FALSE" dependiendo de 
# si se quiere considerar o no al grupo HEMO. 
# -nmut: la cantidad de mutantes a generar.
# -nsitesmut: la cantidad de sitios mutados por proteína mutante.
# -Fmax: valor máximo de las fuerzas usadas para la simulación de las mutaciones.
# -model: el modelo de red elástica a usar, que solo puede ser "ANM" por el momento.
# -R0: cut-off del ANM.

#Remove objects#
rm(list = ls()) 

#Read input#
data.dir <- "DATA/Theoretical"

input.dir <- file.path(data.dir,"InputT.csv")
input <- read.csv(input.dir)

#Reference family and protein#
family <- input$family
p.ref <- input$p.ref 
chain <- input$chain
heme <- input$heme 

#Parameters for mutations#
nmut = input$nmut 
nsitesmut = input$nsitesmut 
fmax = input$fmax 

#ENM#
model <- input$model
R0 = input$R0

#PDB dir#
get.pdb(as.character(p.ref),data.dir) 
pdb.fname <- file.path(data.dir,paste(p.ref,".pdb",sep=""))

#Functions dir#
f.readCA <-"FUNCTIONS/ReadCA.R" 
f.readHeme <- "FUNCTIONS/ReadHeme.R"  
f.Keff <- "FUNCTIONS/Keff.R"  
f.K <- "FUNCTIONS/K.R" 
f.kij.function <- file.path("FUNCTIONS",paste("Kij.",model,".R",sep = "")) 
f.force <- "FUNCTIONS/Force.R"  
f.variability <- "FUNCTIONS/Variability.R"  

#Output dir#
OUT <- "OUT/Theoretical"

#General parameters#
TOLERANCE = 1e-10

#Load librarys#
library(bio3d)

#Load Functions#
source(f.readCA) 
source(f.readHeme) 
source(f.Keff)
source(f.K)
source(f.kij.function)
source(f.force)
source(f.variability)

#Read PDB of p.ref#
pdb <- readCA(pdb.fname,chain)
r.p.ref = pdb$xyz.calpha
nsites = pdb$nsites
naa = nsites 

if ((family == "Globins") & (heme == "TRUE")){
  r.heme = readHeme(pdb.fname,chain)
  r.p.ref = cbind(r.p.ref,r.heme)
  nsites = ncol(r.p.ref)
  align.index = seq(1,naa)
  no.align.index = seq((naa+1),nsites)
}

#Calculate K of p.ref#
K.p.ref <- K(r.p.ref,kij.function,R0,TOLERANCE)

#Variability#
m.nH = matrix(ncol = 3*naa , nrow = nmut)
m.Pn = matrix(ncol = 3*naa , nrow = nmut)
m.d.evalues = matrix(ncol = 3*naa , nrow = nmut)
m.dr.squarei = matrix(ncol = naa , nrow = nmut)

#Count the number of mutants to discard#
count = 0

#Calculate mutants#
for(mut in seq(nmut)){
 	f <- rep(0,3*nsites)
	for (l in (sample(1:naa,replace = F)[1:nsitesmut])){
    print(c(mut,l))
    fl = force(l,r.p.ref,K.p.ref$kij,fmax)
		f = f + fl
	}

  dr = K.p.ref$Cov %*% f
  dim(dr) <- c(3,nsites)
  r.mut = r.p.ref + dr

  #Rotate mutants minimizing RMSD with P.ref#
	r.mut <- matrix(fit.xyz(fixed = as.vector(r.p.ref), 
	mobile = as.vector(r.mut),
	fixed.inds = seq(1,3*naa),
	mobile.inds = seq(1,3*naa)),nrow=3) 
	
	#Calculate K of mutants and p.ref#
  if (naa != nsites) {
    K.mut <- keff(r.mut, align.index,no.align.index, kij.function,R0,TOLERANCE)
    K.p.ref.2 <- keff(r.p.ref,align.index,no.align.index,kij.function,R0,TOLERANCE)
  }
  if (naa == nsites) {
    K.mut <- K(r.mut,kij.function,R0,TOLERANCE)
    K.p.ref.2 <- K.p.ref
  }
	
	nmodes <- length(K.p.ref.2$va)  
	
  #Calculate dr#
  dr = (r.mut - r.p.ref)[,1:naa]   
  
	if (ncol(K.mut$ve) == ncol(K.p.ref.2$ve)){ 
	  
	  #Calculate measures o variability#
		VA <- variability(dr,K.p.ref.2,K.mut)
    m.nH[mut,1:nmodes] = t(VA$nH)
  	m.Pn[mut,1:nmodes] = t(VA$Pn)
    m.d.evalues[mut,1:nmodes]  = t(VA$d.evalues[1:nmodes])
		dr.squarei = rbind(VA$dr.squarei,seq(1,naa))
		for (i in (1:naa)){
			m.dr.squarei[mut,i]= matrix(dr.squarei[1,dr.squarei[2,]==i],ncol=1,nrow=1)
		}
	}
	if (ncol(K.mut$ve) != ncol(K.p.ref.2$ve)){ 
	  count <- count +1
	}
}

#Calculate Means#
mean.nH = colMeans(as.matrix(m.nH),na.rm=T)
mean.Pn = colMeans(as.matrix(m.Pn),na.rm=T)
mean.d.evalues= colMeans(as.matrix(m.d.evalues),na.rm=T)
MSDi = colMeans(as.matrix(m.dr.squarei),na.rm=T)
MSD = rowMeans(m.dr.squarei,na.rm = T)

#SAVE INFORMATION#
if (family == "Globins") family <- paste(family,"Heme",heme,sep ="")

#Dataframes#
write.csv(m.nH,file = file.path(OUT , paste(family , ".Out.m.nH.csv" , sep = "")) , row.names = FALSE)
write.csv(m.Pn,file = file.path(OUT , paste(family , ".Out.m.Pn.csv" , sep="")) , row.names = FALSE)
write.csv(m.d.evalues,file = file.path(OUT , paste(family , ".Out.m.d.evalues.csv" , sep = "")) , row.names = FALSE)
write.csv(m.dr.squarei,file = file.path(OUT , paste(family , ".Out.m.dr.squarei.csv" , sep = "")) , row.names = FALSE)
write.csv(K.p.ref.2$va,file = file.path(OUT , paste(family , ".Out.evalues.csv" , sep = "")) , row.names = FALSE)

#Means#
write.csv(mean.nH,file = file.path(OUT , paste(family , ".Out.nH.mean.csv" , sep = "")) , row.names = FALSE)
write.csv(mean.Pn,file = file.path(OUT , paste(family , ".Out.Pn.mean.csv" , sep = "")) , row.names = FALSE)
write.csv(mean.d.evalues,file = file.path(OUT , paste(family , ".Out.d.evalues.mean.csv" , sep = "")) , row.names = FALSE)
write.csv(MSDi,file = file.path(OUT , paste(family , ".Out.dr.squarei.mean.csv" , sep = "")) , row.names = FALSE)
write.csv(MSD,file = file.path(OUT , paste(family , ".Out.dr.square.mean.csv", sep = "")) , row.names = FALSE)

#input#
write.csv(input,file = file.path(OUT , paste(family , ".Input" , sep = "")) , row.names = FALSE)


