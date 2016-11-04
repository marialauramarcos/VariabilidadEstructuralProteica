# input
# n.seq = 10
# p.ref = "1a6m"
# chain = "A"

# libraries
# bio3d

# get the pdb file of coordinates of p.ref
get.pdb(p.ref)

# read the file
pdb <- read.pdb(paste(p.ref, ".pdb", sep = ""))    
sel <- atom.select(pdb, chain = chain, elety = "CA")    
seq <- pdb$atom[sel$atom, c("resid")]

# make a blast
blast.info <- blast.pdb(seq, database = "pdb")
score = blast.info$mlog.evalue
pdbid = blast.info$pdb.id
identity = blast.info$hit.tbl[, 3]

# filter
filter = (pdbid != p.ref)
pdbid = pdbid[filter][1:n.seq]
score = score[filter][1:n.seq]

# align with muscle

# get indeces

# write output