# INFORMACIÓN:
#
# El programa calcula medidas de variabilidad estructural y dinámica de alineamientos múltiples de 
# familias de proteínas provenientes de la base de datos de alineamientos múltiples "Homstrad".
#
# Para utilizar el programa se debe completar un input ("DATA/Experimental/inputE.csv") especificando: 
# -family:la familia de proteínas del alineamiento, que puede ser "SerinProteases", "Globins" o "Plastocyanins"
# -p.ref: la proteína de referncia que se desee utilizar, teniendo en cuenta que la misma debe ser escrita de 
# la msima forma en la que se encuentra en el dataset.
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
align.dir <- file.path(data.dir,paste(family,"Alignments.txt",sep="")) #File with a multiple sequence alignment from Homstrad#
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

#Read PDB & alignment of P.ref#
chain.Pref <- chains[pdbids == as.character(p.ref)]
pdb.Pref <- readCA(pdbs.dir,chain.Pref)
r.Pref = pdb.Pref$xyz.calpha
nsites.Pref = pdb.Pref$nsites 
naa.Pref = nsites.Pref #naa = n.aminoaids, nsites = n.aminoácids + n.atoms we want to consider, ej.: Ns and Fe from heme group in Globins# 

if ((family == "Globins") & (heme == "TRUE")){
  r.heme.Pref = readHeme(pdbs.dir,chain.Pref)
  r.Pref = cbind(r.Pref,r.heme.Pref)
  nsites.Pref = ncol(r.Pref)
}

alignment.Pref <- alignments[alignments.id$id == as.character(p.ref),]

#Variability#
df.lalign = matrix(ncol = 1 , nrow = nprot)
df.ID = matrix(ncol = 1 , nrow = nprot)
df.nH = matrix(ncol = 3*naa.Pref , nrow = nprot)
df.Pn = matrix(ncol = 3*naa.Pref , nrow = nprot)
df.evalues = matrix(ncol = 3*naa.Pref , nrow = nprot)
df.d.evalues = matrix(ncol = 3*naa.Pref , nrow = nprot)
df.dr.squarei = matrix(ncol = naa.Pref , nrow = nprot)

for (P in (1:nprot)){
  
	#Read PDB & alignment of P2#
	chain.P2 <- chains[[P]]
	pdb.P2 <- readCA(pdbs.dir,chain.P2)
	r.P2 <- pdb.P2$xyz.calpha
	nsites.P2 <- pdb.P2$nsites
	naa.P2 = nsites.P2
	
	if ((family == "Globins") & (heme == "TRUE")){
	  r.heme.P2 = readHeme(pdbs.dir,chain.P2)
	  r.P2 = cbind(r.P2,r.heme.P2)
	  nsites.P2 = ncol(r.P2)
	}

	alignment.P2 <- alignments[alignments.id$id == as.character(pdbids[P]),]	
	
	#Analyze alignment of Pref#
	a.align.Pref <- alignment(alignment.Pref,alignment.P2,naa.Pref)
	align.Pref.index <- a.align.Pref$align.index
	nalign <- a.align.Pref$nalign
	r.align.Pref <- r.Pref[,align.Pref.index]
	no.align.Pref.index <- a.align.Pref$no.align.index
	
	if ((family == "Globins") & (heme == "TRUE")){
	  no.align.Pref.index <- cbind(no.align.Pref.index,t(seq((naa.Pref+1),nsites.Pref)))
	}
	
	#Cakculate KEFF Pref#
	Keff.Pref <- keff(r.Pref,align.Pref.index,no.align.Pref.index,kij.function,R0,TOLERANCE)	
	nmodes <- length(Keff.Pref$va)
	
	#Calculate %ID between P.ref y P2#
	ID.P <- ID(alignment.Pref,alignment.P2)
	
	#Analyze alignment of P2#
	a.align.P2 <- alignment(alignment.P2,alignment.Pref,naa.P2)
	align.P2.index <- a.align.P2$align.index
	no.align.P2.index <- a.align.P2$no.align.index
	
	if ((family == "Globins") & (heme == "TRUE")){
	  no.align.P2.index <- cbind(no.align.P2.index,t(seq((naa.P2+1),nsites.P2)))
	}
	
	#Rotate P2 minimizing RMSD with P.ref#
	align.Pref.index3N = matrix(0 , ncol = 3*nalign, nrow = 1)
	align.P2.index3N = matrix(0 , ncol = 3*nalign, nrow = 1)
  for (i in (1: nalign)){
    align.Pref.index3N[1,((3*i-2):(3*i))] = c((3*align.Pref.index[i]-2),(3*align.Pref.index[i]-1),(3*align.Pref.index[i]))
    align.P2.index3N[1,((3*i-2):(3*i))] = c((3*align.P2.index[i]-2),(3*align.P2.index[i]-1),(3*align.P2.index[i]))
  }

	r.P2 <- matrix(fit.xyz(fixed = as.vector(r.Pref),
	mobile = as.vector(r.P2),
	fixed.inds = align.Pref.index3N,
  mobile.inds = align.P2.index3N
	),nrow=3)
	
	r.align.P2 <- r.P2[,align.P2.index]
	
	#Cakculate KEFF P2#	
	Keff.P2 <- keff(r.P2,align.P2.index,no.align.P2.index,kij.function,R0,TOLERANCE)	
	
	#dr#	
	dr = r.align.P2 - r.align.Pref

	#Variability#
	VA <- variability(dr,Keff.Pref,Keff.P2)
	df.ID[P] = ID.P
  df.lalign[P] = nalign 
  df.nH[P,1:nmodes] = t(VA$nH)
  df.Pn[P,1:nmodes] = t(VA$Pn)
  df.d.evalues[P,1:nmodes]  = t(VA$d.evalues[1:nmodes])
	df.evalues[P,1:nmodes]  = t(Keff.Pref$va[1:nmodes])
	dr.squarei = rbind(VA$dr.squarei,align.Pref.index)
	for (i in (1:naa.Pref)){
		df.dr.squarei[P,i]= matrix(dr.squarei[1,dr.squarei[2,]==i],ncol=1,nrow=1)
	}
}

P.index = (pdbids != as.character(p.ref))

df.ID <- df.ID[P.index,]
df.lalign <- df.lalign[P.index,]
df.nH <- df.nH[P.index,]
df.Pn <- df.Pn[P.index,]
df.d.evalues <- df.d.evalues[P.index,]
df.evalues <- df.evalues[P.index,]
df.dr.squarei <- df.dr.squarei[P.index,]

#Calculate means#
mean.nH = colMeans(df.nH,na.rm=T)
mean.Pn = colMeans(df.Pn,na.rm=T)
mean.d.evalues= colMeans(df.d.evalues,na.rm=T)
mean.evalues= colMeans(df.evalues,na.rm=T)
MSDi = colMeans(df.dr.squarei,na.rm = T)
MSD = rowMeans(df.dr.squarei,na.rm = T)

#SAVE INFORMATION#
if (family == "Globins") family <- paste(family,"Heme",heme,sep ="")

#Dataframes#
write.csv(df.lalign,file = file.path(OUT , paste(family , ".Out.df.lalign.csv" , sep = "")) , row.names = FALSE)
write.csv(df.ID,file = file.path(OUT, paste(family , ".Out.df.ID.csv" , sep = "")) , row.names = FALSE)
write.csv(df.nH,file = file.path(OUT , paste(family , ".Out.df.nH.csv" , sep = "")) , row.names = FALSE)
write.csv(df.Pn,file = file.path(OUT , paste(family , ".Out.df.Pn.csv" , sep="")) , row.names = FALSE)
write.csv(df.d.evalues,file = file.path(OUT , paste(family , ".Out.df.d.evalues.csv" , sep = "")) , row.names = FALSE)
write.csv(df.evalues,file = file.path(OUT , paste(family , ".Out.df.evalues.csv" , sep = "")), row.names = FALSE)
write.csv(df.dr.squarei,file = file.path(OUT , paste(family , ".Out.df.dr.squarei.csv" , sep = "")) , row.names = FALSE)

#Means#
write.csv(mean.nH,file = file.path(OUT , paste(family , ".Out.nH.mean.csv" , sep = "")) , row.names = FALSE)
write.csv(mean.Pn,file = file.path(OUT , paste(family , ".Out.Pn.mean.csv" , sep = "")) , row.names = FALSE)
write.csv(mean.d.evalues,file = file.path(OUT , paste(family , ".Out.d.evalues.mean.csv" , sep = "")) , row.names = FALSE)
write.csv(mean.evalues,file = file.path(OUT , paste(family , ".Out.evalues.mean.csv" , sep = "")) , row.names = FALSE)
write.csv(MSDi,file = file.path(OUT , paste(family , ".Out.MSDi.csv" , sep = "")) , row.names = FALSE)
write.csv(MSD,file = file.path(OUT , paste(family , ".Out.MSD.csv", sep = "")) , row.names = FALSE)

#input#
write.csv(input,file = file.path(OUT , paste(family , ".Input" , sep = "")) , row.names = FALSE)
