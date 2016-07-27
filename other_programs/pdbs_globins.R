# Directories
out.dir <- "pdbs_globins"

# Read the dataset
dataset <- read.csv(dataset.fname)
pdbid.dataset <- as.character(dataset$pdbid)
chain <- as.character(dataset$chain)
n.prot = length(pdbid.dataset)

for (P in (1:n.prot)) {
  print(P)
  
  chain.p.2 <- chain[[P]]
  pdbid.p.2 <- pdbid.dataset[[P]]
  exp.pdb.p.2 = ReadCA(pdbs.fname, chain.p.2)
  exp.r.p.2 = exp.pdb.p.2$xyz.calpha
  
  write.pdb(xyz = c(exp.r.p.2), file = file.path(out.dir, paste(family, "_", pdbid.p.2, ".pdb", sep = "")))
}