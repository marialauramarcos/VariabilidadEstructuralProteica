# INFORMACIÓN:
#
# El programa calcula medidas de variabilidad estructural y dinámica de alineamientos múltiples de 
# familias de proteínas provenientes de la base de datos de alineamientos múltiples "Homstrad".
#
# Para utilizar el programa se debe completar un input ("DATA/Experimental/inputE.csv") especificando: 
# -family:la familia de proteínas del alineamiento, que puede ser "SerinProteases", "Globins" o "Plastocyanins"
# -p.ref: la proteína de referncia que se desee utilizar, teniendo en cuenta que la misma debe ser escrita de 
# la misma forma en la que se encuentra en el dataset.
# -heme: solo se utiliza para la familia de las globinas, puede ser "TRUE" o "FALSE", dependiendo de si 
# se desea o no considerar al grupo HEMO. 
# -model: el modelo de red elástica a usar, que solo puede ser "ANM" por el momento.
# -R0: cut-off del ANM.

#Remove objects#
rm(list = ls())

#Read input#
data.dir <- "DATA/Experimental"

input.dir <- file.path(data.dir,"InputE.csv")
input <- read.csv(input.dir)

#Reference family and protein#
family <- input$family
p.ref <- input$p.ref 
heme = input$heme

#ENM#
model <- input$model
R0 = input$R0

#Files dir#
align.dir <- file.path(data.dir,paste(family,"Alignments.txt",sep="")) #File with a multiple sequence alignmet#
pdbs.dir <- file.path(data.dir,paste(family,"Coordinates.pdb",sep="")) #Coordinates of the proteins# 
dataset.dir <- file.path(data.dir,paste(family,"Dataset.csv",sep="")) #Dataset with pdbids and chains# 

#Functions dir#
f.readCA <-"FUNCTIONS/ReadCA.R" 
f.readHeme <- "FUNCTIONS/ReadHeme.R"  
f.alignment <- "FUNCTIONS/Alignment.R" 
f.ID <- "FUNCTIONS/ID.R" 
f.Keff <- "FUNCTIONS/Keff.R"  
f.K <- "FUNCTIONS/K.R"  
f.kij.function <- file.path("FUNCTIONS",paste("Kij.",model,".R",sep = "")) 
f.variability <- "FUNCTIONS/Variability.R" 

#Output dir#
OUT <- "OUT/Experimental"

#General parameters#
TOLERANCE = 1e-10 

#Load Librarys#
library(seqinr) 
library(bio3d) 

#Load Functions#
source(f.readCA) 
source(f.readHeme) 
source(f.alignment)
source(f.ID)
source(f.Keff)
source(f.K)
source(f.kij.function)
source(f.variability)

#Read Dataset#
Dataset <- read.csv(dataset.dir)
pdbids <- Dataset$pdbid
chains <- Dataset$chain
nprot = length(pdbids)

#Read multiple alignment#
alignments.id <- read.fasta(align.dir)
alignments <- alignments.id$ali[,-ncol(alignments.id$ali)]

#Read PDB & alignment of p.ref#
chain.p.ref <- chains[pdbids == as.character(p.ref)]
pdb.p.ref <- readCA(pdbs.dir,chain.p.ref)
r.p.ref = pdb.p.ref$xyz.calpha
nsites.p.ref = pdb.p.ref$nsites 
naa.p.ref = nsites.p.ref #naa = n.aminoaids, nsites = n.aminoácids + n.atoms we want to consider, ej.: Ns and Fe from heme group in Globins# 

if ((family == "Globins") & (heme == "TRUE")){
  r.heme.p.ref = readHeme(pdbs.dir,chain.p.ref)
  r.p.ref = cbind(r.p.ref,r.heme.p.ref)
  nsites.p.ref = ncol(r.p.ref)
}

alignment.p.ref <- alignments[alignments.id$id == as.character(p.ref),]

#Variability#
m.lalign = matrix(ncol = 1 , nrow = nprot)
m.ID = matrix(ncol = 1 , nrow = nprot)
m.nH = matrix(ncol = 3*naa.p.ref , nrow = nprot)
m.Pn = matrix(ncol = 3*naa.p.ref , nrow = nprot)
m.evalues = matrix(ncol = 3*naa.p.ref , nrow = nprot)
m.d.evalues = matrix(ncol = 3*naa.p.ref , nrow = nprot)
m.dr.squarei = matrix(ncol = naa.p.ref , nrow = nprot)

for (P in (1:nprot)){
  
	#Read PDB & alignment of p.2#
	chain.p.2 <- chains[[P]]
	pdb.p.2 <- readCA(pdbs.dir,chain.p.2)
	r.p.2 <- pdb.p.2$xyz.calpha
	nsites.p.2 <- pdb.p.2$nsites
	naa.p.2 = nsites.p.2
	
	if ((family == "Globins") & (heme == "TRUE")){
	  r.heme.p.2 = readHeme(pdbs.dir,chain.p.2)
	  r.p.2 = cbind(r.p.2,r.heme.p.2)
	  nsites.p.2 = ncol(r.p.2)
	}

	alignment.p.2 <- alignments[alignments.id$id == as.character(pdbids[P]),]	
	
	#Analyze alignment of p.ref#
	a.align.p.ref <- alignment(alignment.p.ref,alignment.p.2,naa.p.ref)
	align.p.ref.index <- a.align.p.ref$align.index
	nalign <- a.align.p.ref$nalign
	r.align.p.ref <- r.p.ref[,align.p.ref.index]
	no.align.p.ref.index <- a.align.p.ref$no.align.index
	
	if ((family == "Globins") & (heme == "TRUE")){
	  no.align.p.ref.index <- cbind(no.align.p.ref.index,t(seq((naa.p.ref+1),nsites.p.ref)))
	}
	
	#Cakculate KEFF p.ref#
	Keff.p.ref <- keff(r.p.ref,align.p.ref.index,no.align.p.ref.index,kij.function,R0,TOLERANCE)	
	nmodes <- length(Keff.p.ref$va)
	
	#Calculate %ID between P.ref y p.2#
	ID.p.2 <- ID(alignment.p.ref,alignment.p.2)
	
	#Analyze alignment of p.2#
	a.align.p.2 <- alignment(alignment.p.2,alignment.p.ref,naa.p.2)
	align.p.2.index <- a.align.p.2$align.index
	no.align.p.2.index <- a.align.p.2$no.align.index
	
	if ((family == "Globins") & (heme == "TRUE")){
	  no.align.p.2.index <- cbind(no.align.p.2.index,t(seq((naa.p.2+1),nsites.p.2)))
	}
	
	#Rotate p.2 minimizing RMSD with P.ref#
	align.p.ref.index3N = matrix(0 , ncol = 3*nalign, nrow = 1)
	align.p.2.index3N = matrix(0 , ncol = 3*nalign, nrow = 1)
  for (i in (1: nalign)){
    align.p.ref.index3N[1,((3*i-2):(3*i))] = c((3*align.p.ref.index[i]-2),(3*align.p.ref.index[i]-1),(3*align.p.ref.index[i]))
    align.p.2.index3N[1,((3*i-2):(3*i))] = c((3*align.p.2.index[i]-2),(3*align.p.2.index[i]-1),(3*align.p.2.index[i]))
  }

	r.p.2 <- matrix(fit.xyz(fixed = as.vector(r.p.ref),
	mobile = as.vector(r.p.2),
	fixed.inds = align.p.ref.index3N,
  mobile.inds = align.p.2.index3N
	),nrow=3)
	
	r.align.p.2 <- r.p.2[,align.p.2.index]
	
	#Cakculate KEFF p.2#	
	Keff.p.2 <- keff(r.p.2,align.p.2.index,no.align.p.2.index,kij.function,R0,TOLERANCE)	
	
	#dr#	
	dr = r.align.p.2 - r.align.p.ref

	#Variability#
	VA <- variability(dr,Keff.p.ref,Keff.p.2)
	m.ID[P] = ID.p.2
  m.lalign[P] = nalign 
  m.nH[P,1:nmodes] = t(VA$nH)
  m.Pn[P,1:nmodes] = t(VA$Pn)
  m.d.evalues[P,1:nmodes]  = t(VA$d.evalues[1:nmodes])
	m.evalues[P,1:nmodes]  = t(Keff.p.ref$va[1:nmodes])
	dr.squarei = rbind(VA$dr.squarei,align.p.ref.index)
	for (i in (1:naa.p.ref)){
		m.dr.squarei[P,i]= matrix(dr.squarei[1,dr.squarei[2,]==i],ncol=1,nrow=1)
	}
}

P.index = (pdbids != as.character(p.ref))

m.ID <- m.ID[P.index,]
m.lalign <- m.lalign[P.index,]
m.nH <- m.nH[P.index,]
m.Pn <- m.Pn[P.index,]
m.d.evalues <- m.d.evalues[P.index,]
m.evalues <- m.evalues[P.index,]
m.dr.squarei <- m.dr.squarei[P.index,]

#Calculate means#
mean.nH = colMeans(m.nH,na.rm=T)
mean.Pn = colMeans(m.Pn,na.rm=T)
mean.d.evalues= colMeans(m.d.evalues,na.rm=T)
mean.evalues= colMeans(m.evalues,na.rm=T)
MSDi = colMeans(m.dr.squarei,na.rm = T)
MSD = rowMeans(m.dr.squarei,na.rm = T)

#SAVE INFORMATION#
if (family == "Globins") family <- paste(family,"Heme",heme,sep ="")

#Dataframes#
write.csv(m.lalign,file = file.path(OUT , paste(family , ".Out.m.lalign.csv" , sep = "")) , row.names = FALSE)
write.csv(m.ID,file = file.path(OUT, paste(family , ".Out.m.ID.csv" , sep = "")) , row.names = FALSE)
write.csv(m.nH,file = file.path(OUT , paste(family , ".Out.m.nH.csv" , sep = "")) , row.names = FALSE)
write.csv(m.Pn,file = file.path(OUT , paste(family , ".Out.m.Pn.csv" , sep="")) , row.names = FALSE)
write.csv(m.d.evalues,file = file.path(OUT , paste(family , ".Out.m.d.evalues.csv" , sep = "")) , row.names = FALSE)
write.csv(m.evalues,file = file.path(OUT , paste(family , ".Out.m.evalues.csv" , sep = "")), row.names = FALSE)
write.csv(m.dr.squarei,file = file.path(OUT , paste(family , ".Out.m.dr.squarei.csv" , sep = "")) , row.names = FALSE)

#Means#
write.csv(mean.nH,file = file.path(OUT , paste(family , ".Out.nH.mean.csv" , sep = "")) , row.names = FALSE)
write.csv(mean.Pn,file = file.path(OUT , paste(family , ".Out.Pn.mean.csv" , sep = "")) , row.names = FALSE)
write.csv(mean.d.evalues,file = file.path(OUT , paste(family , ".Out.d.evalues.mean.csv" , sep = "")) , row.names = FALSE)
write.csv(mean.evalues,file = file.path(OUT , paste(family , ".Out.evalues.mean.csv" , sep = "")) , row.names = FALSE)
write.csv(MSDi,file = file.path(OUT , paste(family , ".Out.MSDi.csv" , sep = "")) , row.names = FALSE)
write.csv(MSD,file = file.path(OUT , paste(family , ".Out.MSD.csv", sep = "")) , row.names = FALSE)

#input#
write.csv(input,file = file.path(OUT , paste(family , ".Input" , sep = "")) , row.names = FALSE)
