# Download "muscle"
source("https://bioconductor.org/biocLite.R")
biocLite("muscle")

# Library
library(muscle)

# Get sequences
alignment.p.ref = alignment.p.ref[alignment.p.ref != "-"]
alignment.p.2 = alignment.p.2[alignment.p.2 != "-"]

# Create a string set
AAstring <- AAStringSet(c("VLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPE
              TLEKFDRFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHH
              EAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQG
              AMNKALELFRKDIAAKYKELGY"
              , "GLSDGEWHLVLNVWGKVETDLAGHGQEVLIRLFKSHPE
TLEKFDKFK-HLKSEDDMRRSEDLRKHGNTVLTALGGILKKKGHH
              EAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSKHPAEFGADAQA
              AMKKALELFRNDIAAKYKELGFHG
              "))

# Align
align = muscle(AAstring)
align.p.ref <- as.vector(align@unmasked[1])
align.p.2 <- as.vector(align@unmasked[2])

alignment.p.ref <- strsplit(align.p.ref, split = "")[[1]]
alignment.p.2 <- strsplit(align.p.2, split = "")[[1]]

# Run AnalyzeFamily.

# Calculate 3N indexes
aligned.p.ref.index.3N = sort(c(aligned.p.ref.index * 3, aligned.p.ref.index * 3 - 2, aligned.p.ref.index * 3 - 1))
aligned.p.2.index.3N = sort(c(aligned.p.2.index * 3, aligned.p.2.index * 3 - 2, aligned.p.2.index * 3 - 1))

# Structures
pdb.p.ref <- read.pdb(get.pdb("1a6m"))
pdb.p.2 <- read.pdb(get.pdb("1mbs"))

p.ref.inds <- atom.select(pdb.p.ref, chain = "A", elety = "CA")
p.2.inds <- atom.select(pdb.p.2, chain = "A", elety = "CA")

r.p.ref <- matrix(pdb.p.ref$xyz[p.ref.inds$xyz], nrow = 3, byrow = F)    
r.p.2 <- matrix(pdb.p.2$xyz[p.2.inds$xyz], nrow = 3, byrow = F)    

r.p.2 = matrix(fit.xyz(fixed = as.vector(r.p.ref),
        mobile = as.vector(r.p.2),
        fixed.inds = aligned.p.1.index.3N, 
        mobile.inds = aligned.p.2.index.3N), nrow = 3)

r.p.1 = as.vector(r.p.ref)
r.p.2 = as.vector(r.p.2)