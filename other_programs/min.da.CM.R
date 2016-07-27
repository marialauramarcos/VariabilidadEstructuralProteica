# Filenames
pdb.fname <- file.path(data.dir, paste(p.ref, ".pdb", sep = "")) 

r.CM.p.ref = CalculateSideChainCM(pdb.fname, chain.p.ref)
r.ca.p.ref = ReadCA(pdb.fname, chain.p.ref)$xyz

# Get heme coordinates, add them to CM´s coordinates and calculate the new number of sites
r.heme = ReadHeme(pdb.fname, chain.p.ref)
r.CM.p.ref = cbind(r.CM.p.ref, r.heme)
r.ca.p.ref = cbind(r.ca.p.ref, r.heme)
n.sites = ncol(r.CM.p.ref)

# Active sites
active = c(28, 29, 32, 64, 68, 93, 107, 152, 153, 154, 155, 156)

# Calculate distance of each center of mass 
m.da.CM = matrix(nrow = n.sites, ncol = length(active))
m.da.ca = matrix(nrow = n.sites, ncol = length(active))

for (i in (active)) {
  print(i)
  for (j in (1:n.sites)) {
    m.da.CM[j, which(active == i)] = sqrt(sum((r.CM.p.ref[, j] - r.CM.p.ref[, i]) ^ 2))
    m.da.ca[j, which(active == i)] = sqrt(sum((r.ca.p.ref[, j] - r.ca.p.ref[, i]) ^ 2))
  }
}

# Calculate the minimum distance
m.min.da.CM = matrix(nrow = n.sites, ncol = 1)
m.min.da.ca = matrix(nrow = n.sites, ncol = 1)

for (i in (1:n.sites)) {
  m.min.da.CM[i, ] = min(m.da.CM[i, ])
  m.min.da.ca[i, ] = min(m.da.ca[i, ])
}

data = data.frame("site" = seq(1:156), "min.da.CM" = as.vector(m.min.da.CM), "min.da.ca" = as.vector(m.min.da.ca))
setwd("C:/Users/Usuario/Desktop/VariabilidadEstructuralProteica")
out.dir = "OUT/out_subset_CA"
write.csv(data, file.path(out.dir, "min.da.CM.ca.csv", sep = ""))
out.dir = "OUT/out_subset_CM"
write.csv(data, file.path(out.dir, "min.da.CM.ca.csv", sep = ""))

