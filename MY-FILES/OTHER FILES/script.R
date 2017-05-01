for (P in (20:23)){
  pdbid.P = pdbid.dataset[P]
       chain.P = chain[P]
       print(pdbid.P)
       print(chain.P)
       
       pdb.P = ReadCA(pdbs.dir, chain.P)
       xyz.P = pdb.P$xyz.calpha
       count = 0
       alignment.P = alignment[which(grepl(pdbid.P, pdbid.alignment)), ]
   print(sum(alignment.P!="-"))
   print(ncol(xyz.P))}

# write pdb files with mean.z.RSDi as b-factor for building images

## read pdb CAs
pdb.xyz.ca = ReadCA(reference.pdb.fname, chain.p.ref)$xyz

## experimental
exp.mean.z.RSD.vec <- vec2resno(as.vector(d.exp$mean.z.RSD), sort(rep(1:ncol(pdb.xyz.ca),3)))
write.pdb(xyz = c(pdb.xyz.ca), b = exp.mean.z.RSD.vec, file =  file.path(data.dir, paste(p.ref, "_exp.pdb", sep = "")))

## mut
mut.mean.z.RSD.vec <- vec2resno(as.vector(d.mut$mean.z.RSD), sort(rep(1:ncol(pdb.xyz.ca),3)))
write.pdb(xyz = c(pdb.xyz.ca), b = mut.mean.z.RSD.vec, file = file.path(data.dir, paste(p.ref, "_mut.pdb", sep = "")))

## mut + sel / medium
medium.mean.z.RSD.vec <- vec2resno(as.vector(d.medium$mean.z.RSD), sort(rep(1:ncol(pdb.xyz.ca),3)))
write.pdb(xyz = c(pdb.xyz.ca), b = medium.mean.z.RSD.vec, file = file.path(data.dir, paste(p.ref, "_medium.pdb", sep = "")))