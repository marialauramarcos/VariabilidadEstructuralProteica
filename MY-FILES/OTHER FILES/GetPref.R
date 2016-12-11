# remove objects
rm(list = ls())

# librarys
library(seqinr)
library(bio3d)

# function fnames
f.readCA.fname <- "FUNCTIONS/readCA.R" 
f.readFasta.fname <- "FUNCTIONS/readFasta.R" 

# load functions
source(f.readCA.fname)
source(f.readFasta.fname)

# data dir
data.dir = "DATA"

# general parameters
TOLERANCE = 1e-10 

# input
input.fname <- "MY-FILES/OTHER FILES/input_GetPref.csv"
input <- read.csv(input.fname)

for (a in (1:nrow(input))) {
  
  family <- as.character(input$family)[a]
    
  # filenames
  pdbs.dir <- file.path(data.dir, paste(family, "_coordinates.pdb", sep = "")) 
  alignment.dir <- file.path(data.dir, paste(family, "_alignment.txt", sep = "")) 
  dataset.fname <- file.path(data.dir, paste(family, "_dataset.csv", sep = ""))  
  
  # read multiple alignment
  alignment.id <- ReadFasta(alignment.dir)
  alignment <- alignment.id$ali[, - ncol(alignment.id$ali)]
  pdbid.alignment <- alignment.id$id
  l.align = ncol(alignment)

  # read the dataset
  dataset = read.csv(dataset.fname)
  pdbid.dataset = as.character(dataset$pdbid)
  chain = as.character(dataset$chain)
  n.prot = length(pdbid.dataset)

  # put xyz in the site
  alignment.xyz = matrix(0, nrow = n.prot, ncol = 3 * l.align)
  n.sites = matrix(0, nrow = n.prot, ncol = 1)

  for (P in (1:n.prot)) {
    print(P)
    pdbid.P = pdbid.dataset[P]
    chain.P = chain[P]
    pdb.P = ReadCA(pdbs.dir, chain.P)
    xyz.P = pdb.P$xyz.calpha
    count = 0
    alignment.P = alignment[which(grepl(pdbid.P, pdbid.alignment)), ]
    for (i in (1:l.align)) {
      if (alignment.P[i] != "-") {
        count = count + 1
        alignment.xyz[P, ((3 * i) - 2)] = xyz.P[1, count]
        alignment.xyz[P, ((3 * i) - 1)] = xyz.P[2, count]
        alignment.xyz[P, (3 * i)] = xyz.P[3, count]
      }
    }
  }

  # calculate means
  mean.xyz <- matrix(0, nrow = 1, ncol = (3 * l.align))

  for (i in (1:(3 * l.align))){
    mean.xyz[, i] = mean(alignment.xyz[abs(alignment.xyz[, i]) > TOLERANCE, i])
  }

  # calculate MSD with respect to the mean xyz
  MSD = matrix(0, nrow = n.prot, ncol = 3)
  for (P in (1:n.prot)){
    print(P)
    count = 0
    for (i in (1:(3 * l.align))){
      if (abs(alignment.xyz[P, i]) > TOLERANCE) {
        count = count + 1
        dr2 = (alignment.xyz[P, i] - mean.xyz[, i]) ^ 2
        MSD[P, 1] = MSD[P, 1] + dr2
        MSD[P, 2] = count / 3
        MSD[P, 3] = MSD[P, 1] / (count / 3) 
      }
    }
  }

  # get p.ref and its chain
  p.ref <- pdbid.dataset[MSD[, 3] == min(MSD[, 3])]
  chain.p.ref <- chain[MSD[, 3] == min(MSD[, 3])]
  
  # write an output file
  out <- data.frame(p.ref, chain.p.ref, MSD[, 3]) 
  write.csv(file = paste("OUT/", family, "_MSD.xyz.mean.csv", sep = ""), out)
}

