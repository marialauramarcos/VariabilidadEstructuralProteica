strong.pdb <- vec2resno(as.vector(strong), sort(rep(1:ncol(xyz.calpha),3)))
write.pdb(xyz = c(xyz.calpha), b = strong, file = file.path(out.dir, paste(family, "_", p.ref, "seq_strongSel", ".pdb", sep = "")))

strong.pdb <- vec2resno(as.vector(r4s), sort(rep(1:ncol(xyz.calpha),3)))
write.pdb(xyz = c(xyz.calpha), b = r4s, file = file.path(out.dir, paste(family, "_", p.ref, "seq_r4s", ".pdb", sep = "")))

# Read file.
sequence = read.csv("DATA/msf_divergence.csv", header = F)
strong = sequence$V5[-1]
r4s = sequence$V7[-1]

# Write pdb files for MSDi.
pdbs.fname <- file.path(data.dir, paste("1a6m.pdb", sep = ""))

pdb <- read.pdb(pdbs.fname, ATOM.only = "TRUE")
sel <- atom.select(pdb)
b.factor <- pdb$atom$b
site = as.numeric(pdb$atom[sel$atom, c("resno")])[1:1203]
strong.pdb <- vec2resno(as.vector(strong), site)
r4s.pdb <- vec2resno(as.vector(r4s), site)
b.factor.strong = b.factor
b.factor.r4s = b.factor

b.factor.strong[1:1203] = strong.pdb * 100
b.factor.r4s[1:1203] = r4s.pdb * 100

write.pdb(pdb = pdb, b = b.factor.strong, file = file.path(out.dir, paste(family, "_", p.ref, "_seq_strong", ".pdb", sep = "")))
write.pdb(pdb = pdb, b = b.factor.r4s, file = file.path(out.dir, paste(family, "_", p.ref, "_seq_r4s", ".pdb", sep = "")))


# Read file.
dinamic = read.csv("DATA/msf_divergence.csv", header = F)
strong = dinamic$V7[-1][1:151]
no = dinamic$V6[-1][1:151]
exp = dinamic$V5[-1][1:151]

# Make 2 the maximum value.
for (i in (1:length(strong))) {
  if (strong[i] >= 1){
    strong[i] = 1
  }
  if (exp[i] >= 1){
    exp[i] = 1
  }
  if (no[i] >= 1){
    no[i] = 1
  }
}

# Write pdb files for MSDi.
pdbs.fname <- file.path(data.dir, paste("1a6m.pdb", sep = ""))

pdb <- read.pdb(pdbs.fname, ATOM.only = "TRUE")
sel <- atom.select(pdb)
b.factor <- pdb$atom$b
site = as.numeric(pdb$atom[sel$atom, c("resno")])[1:1203]
strong.pdb <- as.numeric(vec2resno(as.vector(strong), site))
exp.pdb <- as.numeric(vec2resno(as.vector(exp), site))
no.pdb <- as.numeric(vec2resno(as.vector(no), site))

b.factor.strong = b.factor
b.factor.exp = b.factor
b.factor.no = b.factor



b.factor.strong[1:1203] = strong.pdb * 1000
b.factor.exp[1:1203] = exp.pdb * 1000
b.factor.no[1:1203] = no.pdb * 1000

write.pdb(pdb = pdb, b = b.factor.strong, file = file.path(out.dir, paste(family, "_", p.ref, "_msf_strong", ".pdb", sep = "")))
write.pdb(pdb = pdb, b = b.factor.exp, file = file.path(out.dir, paste(family, "_", p.ref, "_msf_exp", ".pdb", sep = "")))
write.pdb(pdb = pdb, b = b.factor.no, file = file.path(out.dir, paste(family, "_", p.ref, "_msf_no_sel", ".pdb", sep = "")))

