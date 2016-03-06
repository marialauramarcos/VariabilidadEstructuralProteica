# remove objects.
rm(list = ls())

# Data dir.
data.dir = "DATA"

# Input.
input.fname <- "input_GetPref.csv"
input <- read.csv(input.fname)
for (a in (1:nrow(input))) {
    family <- as.character(input$family)[a]
    
  # Filenames.
  pdbs.dir <- file.path(data.dir, paste(family, "_coordinates.pdb", sep = "")) 
  alignment.dir <- file.path(data.dir, paste(family, "_alignment.txt", sep = "")) 
  dataset.fname <- file.path(data.dir, paste(family, "_dataset.csv", sep = ""))  

  # Function fnames.
  f.readCA.fname <- "FUNCTIONS/readCA.R" 
  f.readFasta.fname <- "FUNCTIONS/readFasta.R" 

  # General parameters.
  TOLERANCE = 1e-10 
  
  # Librarys.
  library(seqinr)
  library(bio3d)

  # Load functions.
  source(f.readCA.fname)
  source(f.readFasta.fname)

  # Read multiple alignment.
  alignment.id <- read.fasta(alignment.dir)
  alignment <- alignment.id$ali[, - ncol(alignment.id$ali)]
  pdbid.alignment <- alignment.id$id
  l.align = ncol(alignment)

  # Read the dataset.
  dataset = read.csv(dataset.fname)
  pdbid.dataset = as.character(dataset$pdbid)
  chain = as.character(dataset$chain)
  n.prot = length(pdbid.dataset)

  # Put xyz in the site.
  alignment.xyz = matrix(0, nrow = n.prot, ncol = 3 * l.align)
  n.sites = matrix(0, nrow = n.prot, ncol = 1)

  for (P in (1:n.prot)) {
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

  # Calculate means.
  mean.xyz <- matrix(0, nrow = 1, ncol = (3 * l.align))

  for (i in (1:(3 * l.align))){
    mean.xyz[, i] = mean(alignment.xyz[abs(alignment.xyz[, i]) > TOLERANCE, i])
  }

  MSD = matrix(0, nrow = n.prot, ncol = 3)

  for (P in (1:n.prot)){
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

  p.ref <- pdbid.dataset[MSD[, 3] == min(MSD[, 3])]
  chain.p.ref <- chain[MSD[, 3] == min(MSD[, 3])]
  out <- data.frame(p.ref, chain.p.ref, MSD[, 3]) 
  write.csv(file = paste("OUT/", family, "_MSD.xyz.mean.csv", sep = ""), out)
}

