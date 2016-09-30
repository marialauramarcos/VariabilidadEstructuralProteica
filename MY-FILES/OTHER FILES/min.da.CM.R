# Read input
input.fname <- file.path("input_MainProgram.csv")
input <- read.csv(input.fname)

# Data dir
data.dir <- "DATA"

# Start a loop to analyze each family
for (f in (1:nrow(input))) { 
  p.ref <- as.character(input$p.ref)[f]
  heme <- input$heme[f]
  
  # read functional sites
  active <- read.csv(file.path(data.dir, paste(p.ref, "_functionalSites.csv"), sep = ""))$functional.sites

  # pdb filename
  pdb.fname <- file.path(data.dir, paste(p.ref, ".pdb", sep = "")) 

  # get coordinates of p.ref
  r.CM.p.ref = CalculateSideChainCM(pdb.fname, chain.p.ref)
  r.ca.p.ref = ReadCA(pdb.fname, chain.p.ref)$xyz

  if (heme = T) {
    r.heme = ReadHeme(pdb.fname, chain.p.ref)
    r.CM.p.ref = cbind(r.CM.p.ref, r.heme)
    r.ca.p.ref = cbind(r.ca.p.ref, r.heme)
  }
  
  # calculate the number of sites
  n.sites = ncol(r.CM.p.ref)

  # calculate the distances between each site and each functional site (ca or CM)
  m.da.CM = matrix(nrow = n.sites, ncol = length(active))
  m.da.ca = matrix(nrow = n.sites, ncol = length(active))

  for (i in (active)) {
    print(i)
    for (j in (1:n.sites)) {
      m.da.CM[j, which(active == i)] = sqrt(sum((r.CM.p.ref[, j] - r.CM.p.ref[, i]) ^ 2))
      m.da.ca[j, which(active == i)] = sqrt(sum((r.ca.p.ref[, j] - r.ca.p.ref[, i]) ^ 2))
    }
  }

  # calculate the minimum distance of each site to the active sites
  m.min.da.CM = matrix(nrow = n.sites, ncol = 1)
  m.min.da.ca = matrix(nrow = n.sites, ncol = 1)

  for (i in (1:n.sites)) {
    m.min.da.CM[i, ] = min(m.da.CM[i, ])
    m.min.da.ca[i, ] = min(m.da.ca[i, ])
  }

  # build a dataframe
  data = data.frame("site" = seq(1:n.sites), 
               "min.da.CM" = as.vector(m.min.da.CM), 
               "min.da.ca" = as.vector(m.min.da.ca))

  # write output files
  out.dir = "OUT/out_subset_CA"
  write.csv(data, file.path(out.dir, paste(p.ref, "min.da.CM.ca.csv"), sep = ""))
  write.csv(m.da.ca, file.path(out.dir, paste(p.ref, "m.da.ca.csv"), sep = ""))
  write.csv(m.da.CM, file.path(out.dir, paste(p.ref, "m.da.CM.csv"), sep = ""))

  out.dir = "OUT/out_subset_CM" 
  write.csv(data, file.path(out.dir, paste(p.ref, "min.da.CM.ca.csv"), sep = ""))
  write.csv(m.da.ca, file.path(out.dir, paste(p.ref, "m.da.ca.csv"), sep = ""))
  write.csv(m.da.CM, file.path(out.dir, paste(p.ref, "m.da.CM.csv"), sep = ""))
}