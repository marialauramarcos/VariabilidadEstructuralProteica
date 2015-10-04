# INFORMACIÓN:
#
# El programa calcula medidas de variabilidad estructural y dinámica de mutantes múltiples de
# una proteína determinada utilizando el modelo mutacional Linearly Forced - Elastic Network Model 
# (LF-ENM) y considerando a las fuerzas que modelan cada mutación aditivas entre sí. 
#
# Para utilizar el programa se debe completar un input ("DATA/Theoretical/inputT.csv") especificando: 
# -family: la familia de proteínas a la que pertenece la proteína a mutar (en inglés y plural).
# -p.ref: el código de pdb (pdbid) de la proteína a mutar.
# -chain: la cadena de p.ref que se desea mutar.
# -heme: solo se utiliza para la familia de las globinas, puede ser "TRUE" o "FALSE" dependiendo de 
# si se quiere considerar o no al grupo HEMO. 
# -nmut: la cantidad de mutantes a generar.
# -nsitiosmut: la cantidad de sitios mutados por proteína mutante.
# -Fmax: valor máximo de las fuerzas usadas para la simulación de las mutaciones.
# -model: el modelo de red elástica, que solo puede ser "ANM" por el momento.
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

#Files dir#
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
r.wt = pdb$xyz.calpha
nsites = pdb$nsites
naa = nsites 

if ((family == "Globins") & (heme == "TRUE")){
  r.heme = readHeme(pdb.fname,chain)
  r.wt = cbind(r.wt,r.heme)
  nsites = ncol(r.wt)
}

#Calculate K of p.ref#
K.wt <- K(r.wt,kij.function,R0,TOLERANCE)

#Variability#
df.nH = matrix(ncol = 3*naa , nrow = nmut)
df.Pn = matrix(ncol = 3*naa , nrow = nmut)
df.d.evalues = matrix(ncol = 3*naa , nrow = nmut)
df.dr.squarei = matrix(ncol = naa , nrow = nmut)

#Count the number of mutants to discard#
count = 0

#Calculate mutants#
for(mut in seq(nmut)){
 	f <- rep(0,3*nsites)
	for (l in (sample(1:naa,replace = F)[1:nsitesmut])){
    print(c(mut,l))
    fl = force(l,r.wt,K.wt$kij,fmax)
		f = f + fl
	}

  dr = K.wt$Cov %*% f
  dim(dr) <- c(3,nsites)
  r.mut = r.wt + dr

  #Rotate mutants minimizing RMSD with P.ref#
	r.mut <- matrix(fit.xyz(fixed = as.vector(r.wt), 
	mobile = as.vector(r.mut),
	fixed.inds = seq(1,3*naa),
	mobile.inds = seq(1,3*naa)),nrow=3) 
	
	#Calculate K of mutants and p.ref#
	if (nsites != naa){ #Family = "Globins", heme = "TRUE"
    K.mut <- keff(r.mut,seq(1,naa),seq((naa+1),nsites), kij.function,R0,TOLERANCE)
	  K.wt.2 <- keff(r.wt,seq(1,naa),seq((naa+1),nsites),kij.function,R0,TOLERANCE)
	}
	if (nsites == naa){
	  K.mut <- K(r.mut,kij.function,R0,TOLERANCE)
	  K.wt.2 <- K.wt
	}
	nmodes <- length(K.wt.2$va)  
	
  #Calculate dr#
  dr = (r.mut - r.wt)[,1:naa]   
  
	if (ncol(K.mut$ve) == ncol(K.wt.2$ve)){ 
	  
	  #Calculate measures o variability#
		VA <- variability(dr,K.wt.2,K.mut)
    df.nH[mut,1:nmodes] = t(VA$nH)
  	df.Pn[mut,1:nmodes] = t(VA$Pn)
    df.d.evalues[mut,1:nmodes]  = t(VA$d.evalues[1:nmodes])
		dr.squarei = rbind(VA$dr.squarei,seq(1,naa))
		for (i in (1:naa)){
			df.dr.squarei[mut,i]= matrix(dr.squarei[1,dr.squarei[2,]==i],ncol=1,nrow=1)
		}
	}
	if (ncol(K.mut$ve) != ncol(K.wt.2$ve)){ 
	  count <- count +1
	}
}

#Calculate Means#
mean.nH = colMeans(as.matrix(df.nH),na.rm=T)
mean.Pn = colMeans(as.matrix(df.Pn),na.rm=T)
mean.d.evalues= colMeans(as.matrix(df.d.evalues),na.rm=T)
MSDi = colMeans(as.matrix(df.dr.squarei),na.rm=T)
MSD = rowMeans(df.dr.squarei,na.rm = T)

#SAVE INFORMATION#
if (family == "Globins") p.ref <- paste(p.ref,"Heme",heme,sep ="")

#Dataframes#
write.csv(df.nH,file = file.path(OUT , paste(p.ref , ".Out.df.nH.csv" , sep = "")) , row.names = FALSE)
write.csv(df.Pn,file = file.path(OUT , paste(p.ref , ".Out.df.Pn.csv" , sep="")) , row.names = FALSE)
write.csv(df.d.evalues,file = file.path(OUT , paste(p.ref , ".Out.df.d.evalues.csv" , sep = "")) , row.names = FALSE)
write.csv(df.dr.squarei,file = file.path(OUT , paste(p.ref , ".Out.df.dr.squarei.csv" , sep = "")) , row.names = FALSE)
write.csv(K.wt.2$va,file = file.path(OUT , paste(p.ref , ".Out.evalues.csv" , sep = "")) , row.names = FALSE)

#Means#
write.csv(mean.nH,file = file.path(OUT , paste(p.ref , ".Out.nH.mean.csv" , sep = "")) , row.names = FALSE)
write.csv(mean.Pn,file = file.path(OUT , paste(p.ref , ".Out.Pn.mean.csv" , sep = "")) , row.names = FALSE)
write.csv(mean.d.evalues,file = file.path(OUT , paste(p.ref , ".Out.d.evalues.mean.csv" , sep = "")) , row.names = FALSE)
write.csv(MSDi,file = file.path(OUT , paste(p.ref , ".Out.dr.squarei.mean.csv" , sep = "")) , row.names = FALSE)
write.csv(MSD,file = file.path(OUT , paste(p.ref , ".Out.dr.square.mean.csv", sep = "")) , row.names = FALSE)

#input#
write.csv(input,file = file.path(OUT , paste(p.ref , ".Input" , sep = "")) , row.names = FALSE)


