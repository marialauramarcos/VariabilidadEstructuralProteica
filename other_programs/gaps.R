gaps = colSums(matrix(as.numeric(is.na(m.exp.smooth.dr.squarei)), nrow=n.prot))[1:151]
m.exp.smooth.dr.squarei = m.exp.smooth.dr.squarei[,1:151]

layout(matrix(1:16,4,4))
for (i in (1:n.prot)) {
  plot(x=gaps, y = m.exp.smooth.dr.squarei[i,])
}

MSDi.gaps = rbind(m.exp.smooth.dr.squarei, gaps)
write.csv(MSDi.gaps, file = file.path(out.dir, paste(analysis.fname.id, "_out_m.MSDi.gaps.csv", sep = "")), row.names = FALSE)
MSDi.gaps = rbind(m.exp.smooth.dr.squarei, gaps)
active.sites = list(type=rep("Non active", 151))
active.sites$type[28] = "Active"
active.sites$type[29] = "Active"
active.sites$type[32] = "Active"
active.sites$type[64] = "Active"
active.sites$type[68] = "Active"
active.sites$type[93] = "Active"
active.sites$type[107] = "Active"
active.sites$type[18] = "Gap"
active.sites$type[19] = "Gap"
active.sites$type[47] = "Gap"
active.sites$type[48] = "Gap"
active.sites$type[49] = "Gap"
active.sites$type[50] = "Gap"
active.sites$type[51] = "Gap"
active.sites$type[52] = "Gap"
active.sites$type[53] = "Gap"
active.sites$type[54] = "Gap"
active.sites$type[55] = "Gap"
active.sites$type[148] = "Gap"
active.sites$type[149] = "Gap"
active.sites$type[150] = "Gap"
active.sites$type[151] = "Gap"

write.csv(active.sites, file = file.path("DATA", paste("globins", "active_sites.csv", sep = "")), row.names = FALSE)
