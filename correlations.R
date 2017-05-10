# input
R0S = c(7.5, 10, 12.5)
data.dir = paste("OUT/out_subset_CA_ANM", sep = "")

# read input
input.fname <- "input_MainReport.csv"
input <- read.csv(input.fname)

cc.z.RSD.exp.mut = matrix(ncol = 3, nrow = length(input$family))
cc.z.RSD.exp.medium = matrix(ncol = 3, nrow = length(input$family))
cc.z.RSD.mut.medium = matrix(ncol = 3, nrow = length(input$family))

cc.Pn.exp.mut = matrix(ncol = 3, nrow = length(input$family))
cc.Pn.exp.medium = matrix(ncol = 3, nrow = length(input$family))
cc.Pn.mut.medium = matrix(ncol = 3, nrow = length(input$family))

cc.Pn.100.exp.mut = matrix(ncol = 3, nrow = length(input$family))
cc.Pn.100.exp.medium = matrix(ncol = 3, nrow = length(input$family))
cc.Pn.100.mut.medium = matrix(ncol = 3, nrow = length(input$family))

cc.MSF.exp.mut = matrix(ncol = 3, nrow = length(input$family))
cc.MSF.exp.medium = matrix(ncol = 3, nrow = length(input$family))
cc.MSF.mut.medium = matrix(ncol = 3, nrow = length(input$family))

cc.nH.exp.mut = matrix(ncol = 3, nrow = length(input$family))
cc.nH.exp.medium = matrix(ncol = 3, nrow = length(input$family))
cc.nH.mut.medium = matrix(ncol = 3, nrow = length(input$family))

cc.nH.100.exp.mut = matrix(ncol = 3, nrow = length(input$family))
cc.nH.100.exp.medium = matrix(ncol = 3, nrow = length(input$family))
cc.nH.100.mut.medium = matrix(ncol = 3, nrow = length(input$family))

# satart a loop for each family
for (f in (1:nrow(input))) { 
  print(f)
  family <- as.character(input$family)[f]
  type <- as.character(input$type)[f]
  p.ref <- as.character(input$p.ref)[f]
  enm <- as.character(input$enm)[f]
  n.mut.p <- input$n.mut.p[f]
  R0.CA = input$R0.CA[f]
  R0.CM = input$R0.CM[f]
  chain.p.ref <- as.character(input$chain.p.ref)[f]
  print(family)
  
  for (R0 in (R0S)) {

    ### z.RSD ###
    
    # set input filenames
    dri2.exp.fname = file.path(data.dir, paste(family, "_R0_", R0, "_beta_no.sel_K.analysis_Keff_out_m.exp.norm.dr.squarei.csv", sep = ""))
    dri2.medium.fname = file.path(data.dir, paste(family, "_R0_", R0, "_beta_medium.sel_K.analysis_Keff_out_m.theo.norm.dr.squarei.csv", sep = ""))
    dri2.mut.fname = file.path(data.dir,paste(family, "_R0_", R0, "_beta_no.sel_K.analysis_Keff_out_m.theo.norm.dr.squarei.csv", sep = ""))
    
    # read data
    dri2.exp = read.csv(dri2.exp.fname, header = T)
    dri2.mut = read.csv(dri2.mut.fname, header = T)
    dri2.medium = read.csv(dri2.medium.fname, header = T)
    n.sites = length(dri2.exp)
    
    # prepare data
    rownames(dri2.exp) = paste("mut", rownames(dri2.exp), sep = "")
    rownames(dri2.mut) = paste("mut", rownames(dri2.mut), sep = "")
    rownames(dri2.medium) = paste("mut", rownames(dri2.medium), sep = "")
    
    tdri2.exp = as.data.frame(t(as.matrix(dri2.exp)))
    tdri2.mut = as.data.frame(t(as.matrix(dri2.mut)))
    tdri2.medium = as.data.frame(t(as.matrix(dri2.medium)))
    
    site.info = seq(1:length(dri2.mut))
    
    tdri2.exp = cbind(site.info, tdri2.exp)
    tdri2.mut = cbind(site.info, tdri2.mut)
    tdri2.medium = cbind(site.info, tdri2.medium)
    
    siteComparison.exp = melt(tdri2.exp, id.vars = c("site.info"), variable.name = "protein", value.name = "dri2")
    siteComparison.mut = melt(tdri2.mut, id.vars = c("site.info"), variable.name = "protein", value.name = "dri2")
    siteComparison.medium = melt(tdri2.medium, id.vars = c("site.info"), variable.name = "protein", value.name = "dri2")
    
    # bluild the data.frame
    siteComparison.all.datasets = rbind(data.frame("dataset" = "exp", siteComparison.exp),
                                        data.frame("dataset" = "mut", siteComparison.mut),
                                        data.frame("dataset" = "medium", siteComparison.medium))
    
    siteComparison.all.datasets = ddply(siteComparison.all.datasets, c("dataset", "protein"), mutate,
                                        "n.SD" = dri2/mean(dri2, na.rm = T),
                                        "z.SD" = vscale(dri2),
                                        "RSD" = sqrt(dri2),
                                        "n.RSD" = RSD/mean(RSD, na.rm = T),
                                        "z.RSD" = vscale(RSD))
    
    dat = siteComparison.all.datasets
    dat$dataset = revalue(dat$dataset, c("medium" = "selection"))
    dat$dataset = revalue(dat$dataset, c("mut" = "no-selection"))
    
    ## dat with means and standard errors
    d = ddply(dat,c("dataset", "site.info"), function(x) {
        mean.z.RSD = mean(x$z.RSD, na.rm = T)
        ndata = sum(!is.na(x$z.RSD))
        se.z.RSD = sd(x$z.RSD, na.rm = T) / sqrt(ndata)
        data.frame(mean.z.RSD, ndata, se.z.RSD)
    })
    
    # calculate correlations berween profiles
    d.exp = subset(d, dataset == "exp")
    d.mut = subset(d, dataset == "no-selection")
    d.medium = subset(d, dataset == "selection")
    
    r.exp.mut = cor(d.exp$mean.z.RSD, d.mut$mean.z.RSD)
    r.exp.medium = cor(d.exp$mean.z.RSD, d.medium$mean.z.RSD)
    r.mut.medium = cor(d.mut$mean.z.RSD, d.medium$mean.z.RSD)
    
    cc.z.RSD.exp.mut[which(input$family == family), which(R0S == R0)] = r.exp.mut
    cc.z.RSD.exp.medium[which(input$family == family), which(R0S == R0)] = r.exp.medium
    cc.z.RSD.mut.medium[which(input$family == family), which(R0S == R0)] = r.mut.medium
    
    ### Pn ###
    
    # set input filenames
    Pn.exp.fname = file.path(data.dir, paste(family, "_R0_", R0, "_beta_no.sel_K.analysis_Keff_out_m.exp.Pn.csv", sep = ""))
    Pn.medium.fname = file.path(data.dir, paste(family, "_R0_", R0, "_beta_medium.sel_K.analysis_Keff_out_m.theo.Pn.csv", sep = ""))
    Pn.mut.fname = file.path(data.dir,paste(family, "_R0_", R0, "_beta_no.sel_K.analysis_Keff_out_m.theo.Pn.csv", sep = ""))
    
    # prepare the data
    Pn.exp = read.csv(Pn.exp.fname, header = T)
    rownames(Pn.exp) = paste("mut", rownames(Pn.exp), sep = "")
    Pn.mut = read.csv(Pn.mut.fname, header = T)
    rownames(Pn.mut) = paste("mut", rownames(Pn.mut), sep = "")
    Pn.medium = read.csv(Pn.medium.fname, header = T)
    rownames(Pn.medium) = paste("mut", rownames(Pn.medium), sep = "")
    
    mode = seq(ncol(Pn.medium))
    
    # Prepare modeComparison
    tPn.medium = as.data.frame(t(as.matrix(Pn.medium)))
    tPn.medium = cbind(data.frame("mode" = mode), tPn.medium)
    modeComparison.medium = melt(tPn.medium, "mode", variable.name = "protein", value.name = "Pn")
    modeComparison.medium = ddply(modeComparison.medium, "protein", mutate, "Pn" = Pn / mean(Pn, na.rm = T))
    
    tPn.mut = as.data.frame(t(as.matrix(Pn.mut)))
    tPn.mut = cbind(data.frame("mode" = mode), tPn.mut)
    modeComparison.mut = melt(tPn.mut, "mode", variable.name = "protein", value.name = "Pn")
    modeComparison.mut = ddply(modeComparison.mut, "protein", mutate, "Pn" = Pn / mean(Pn, na.rm = T))
    
    tPn.exp = as.data.frame(t(as.matrix(Pn.exp)))
    tPn.exp = cbind(data.frame("mode" = mode), tPn.exp)
    modeComparison.exp = melt(tPn.exp, "mode", variable.name = "protein", value.name = "Pn")
    modeComparison.exp = ddply(modeComparison.exp, "protein", mutate, "Pn" = Pn / mean(Pn, na.rm = T))
    
    # modeComparison.exp has NA's because for proteins with gaps when using Keff there're less modes
    #print("Warning: modeComparison.exp <- na.omit(modeComparison.exp)")
    modeComparison.exp <- na.omit(modeComparison.exp)
    modeComparison.medium <- na.omit(modeComparison.medium)
    modeComparison.mut <- na.omit(modeComparison.mut)
    
    modeComparison.all.datasets = rbind(data.frame("dataset" = "exp", modeComparison.exp),
                                        data.frame("dataset" = "mut", modeComparison.mut),
                                        data.frame("dataset" = "medium", modeComparison.medium))
    
    dat = modeComparison.all.datasets
    dat$dataset = revalue(dat$dataset, c("medium" = "selection"))
    dat$dataset = revalue(dat$dataset, c("mut" = "no-selection"))
    
    # calculate dataset of medians
    d = ddply(dat,c("dataset", "mode"), function(x) {
        median.Pn = median(x$Pn)
        data.frame(median.Pn)
    })
    
    # calculate correlations berween profiles
    d.exp = subset(d, dataset == "exp")
    d.mut = subset(d, dataset == "no-selection")
    d.medium = subset(d, dataset == "selection")
    
    r.exp.mut = cor(d.exp$median.Pn, d.mut$median.Pn)
    r.exp.medium = cor(d.exp$median.Pn, d.medium$median.Pn)
    r.mut.medium = cor(d.mut$median.Pn, d.medium$median.Pn)
    
    cc.Pn.exp.mut[which(input$family == family), which(R0S == R0)] = r.exp.mut
    cc.Pn.exp.medium[which(input$family == family), which(R0S == R0)] = r.exp.medium
    cc.Pn.mut.medium[which(input$family == family), which(R0S == R0)] = r.mut.medium
    
    # 100 modes
    r.exp.mut = cor(d.exp$median.Pn[1:100], d.mut$median.Pn[1:100])
    r.exp.medium = cor(d.exp$median.Pn[1:100], d.medium$median.Pn[1:100])
    r.mut.medium = cor(d.mut$median.Pn[1:100], d.medium$median.Pn[1:100])
    
    cc.Pn.100.exp.mut[which(input$family == family), which(R0S == R0)] = r.exp.mut
    cc.Pn.100.exp.medium[which(input$family == family), which(R0S == R0)] = r.exp.medium
    cc.Pn.100.mut.medium[which(input$family == family), which(R0S == R0)] = r.mut.medium
    
    ### MSF ###
    
    # set input filenames
    square.dif.MSF.exp.fname = file.path(data.dir, paste(family, "_R0_", R0, "_beta_no.sel_K.analysis_Keff_out_m.exp.square.dif.MSF.csv", sep = ""))
    square.dif.MSF.medium.fname = file.path(data.dir, paste(family, "_R0_", R0, "_beta_medium.sel_K.analysis_Keff_out_m.theo.square.dif.MSF.csv", sep = ""))
    square.dif.MSF.mut.fname = file.path(data.dir,paste(family, "_R0_", R0, "_beta_no.sel_K.analysis_Keff_out_m.theo.square.dif.MSF.csv", sep = ""))
    
    # prepare the data
    square.dif.MSF.exp = read.csv(square.dif.MSF.exp.fname, header = T)
    rownames(square.dif.MSF.exp) = paste("mut", rownames(square.dif.MSF.exp), sep = "")
    square.dif.MSF.mut = read.csv(square.dif.MSF.mut.fname, header = T)
    rownames(square.dif.MSF.mut) = paste("mut", rownames(square.dif.MSF.mut), sep = "")
    square.dif.MSF.medium = read.csv(square.dif.MSF.medium.fname, header = T)
    rownames(square.dif.MSF.medium) = paste("mut", rownames(square.dif.MSF.medium), sep = "")
    
    site = seq(ncol(square.dif.MSF.medium))
    
    # Prepare siteComparison
    tsquare.dif.MSF.medium = as.data.frame(t(as.matrix(square.dif.MSF.medium)))
    tsquare.dif.MSF.medium = cbind(data.frame("site" = site), tsquare.dif.MSF.medium)
    siteComparison.medium = melt(tsquare.dif.MSF.medium, "site", variable.name = "protein", value.name = "square.dif.MSF")
    
    tsquare.dif.MSF.mut = as.data.frame(t(as.matrix(square.dif.MSF.mut)))
    tsquare.dif.MSF.mut = cbind(data.frame("site" = site), tsquare.dif.MSF.mut)
    siteComparison.mut = melt(tsquare.dif.MSF.mut, "site", variable.name = "protein", value.name = "square.dif.MSF")
    
    tsquare.dif.MSF.exp = as.data.frame(t(as.matrix(square.dif.MSF.exp)))
    tsquare.dif.MSF.exp = cbind(data.frame("site" = site), tsquare.dif.MSF.exp)
    siteComparison.exp = melt(tsquare.dif.MSF.exp, "site", variable.name = "protein", value.name = "square.dif.MSF")
    
    siteComparison.all.datasets = rbind(data.frame("dataset" = "exp", siteComparison.exp),
                                        data.frame("dataset" = "mut", siteComparison.mut),
                                        data.frame("dataset" = "medium", siteComparison.medium))
    
    siteComparison.all.datasets = ddply(siteComparison.all.datasets, c("dataset", "protein"), mutate,
                                             "r.SDF" = rank(square.dif.MSF, na.last = "keep")) 
    
    dat = siteComparison.all.datasets
    dat$dataset = revalue(dat$dataset, c("medium" = "selection"))
    dat$dataset = revalue(dat$dataset, c("mut" = "no-selection"))
    
    # calculate dataset of means
    d = ddply(dat, c("dataset", "site"), function(x) {
        mean.r.SDF = mean(x$r.SDF, na.rm = T)
        data.frame(mean.r.SDF)
    })
    
    # calculate correlations berween profiles
    d.exp = subset(d, dataset == "exp")
    d.mut = subset(d, dataset == "no-selection")
    d.medium = subset(d, dataset == "selection")
    
    r.exp.mut = cor(d.exp$mean.r.SDF, d.mut$mean.r.SDF)
    r.exp.medium = cor(d.exp$mean.r.SDF, d.medium$mean.r.SDF)
    r.mut.medium = cor(d.mut$mean.r.SDF, d.medium$mean.r.SDF)
    
    cc.MSF.exp.mut[which(input$family == family), which(R0S == R0)] = r.exp.mut
    cc.MSF.exp.medium[which(input$family == family), which(R0S == R0)] = r.exp.medium
    cc.MSF.mut.medium[which(input$family == family), which(R0S == R0)] = r.mut.medium
    
    ### nH ###
    # set input filenames
    nH.exp.fname = file.path(data.dir, paste(family, "_R0_", R0, "_beta_no.sel_K.analysis_Keff_out_m.exp.nH.csv", sep = ""))
    nH.medium.fname = file.path(data.dir, paste(family, "_R0_", R0, "_beta_medium.sel_K.analysis_Keff_out_m.theo.nH.csv", sep = ""))
    nH.mut.fname = file.path(data.dir,paste(family, "_R0_", R0, "_beta_no.sel_K.analysis_Keff_out_m.theo.nH.csv", sep = ""))
    
    # prepare the data
    nH.exp = read.csv(nH.exp.fname, header = T)
    rownames(nH.exp) = paste("mut", rownames(nH.exp), sep = "")
    nH.mut = read.csv(nH.mut.fname, header = T)
    rownames(nH.mut) = paste("mut", rownames(nH.mut), sep = "")
    nH.medium = read.csv(nH.medium.fname, header = T)
    rownames(nH.medium) = paste("mut", rownames(nH.medium), sep = "")
    
    mode = seq(ncol(nH.medium))
    
    # Prepare modeComparison
    tnH.medium = as.data.frame(t(as.matrix(nH.medium)))
    tnH.medium = cbind(data.frame("mode" = mode), tnH.medium)
    modeComparison.medium = melt(tnH.medium, "mode", variable.name = "protein", value.name = "nH")
    modeComparison.medium = ddply(modeComparison.medium, "protein", mutate, "nH" = nH / mean(nH, na.rm = T))
    
    tnH.mut = as.data.frame(t(as.matrix(nH.mut)))
    tnH.mut = cbind(data.frame("mode" = mode), tnH.mut)
    modeComparison.mut = melt(tnH.mut, "mode", variable.name = "protein", value.name = "nH")
    modeComparison.mut = ddply(modeComparison.mut, "protein", mutate, "nH" = nH / mean(nH, na.rm = T))
    
    tnH.exp = as.data.frame(t(as.matrix(nH.exp)))
    tnH.exp = cbind(data.frame("mode" = mode), tnH.exp)
    modeComparison.exp = melt(tnH.exp, "mode", variable.name = "protein", value.name = "nH")
    modeComparison.exp = ddply(modeComparison.exp, "protein", mutate, "nH" = nH / mean(nH, na.rm = T))
    
    # modeComparison.exp has NA's because for proteins with gaps when using Keff there're less modes
    #print("Warning: modeComparison.exp <- na.omit(modeComparison.exp)")
    modeComparison.exp <- na.omit(modeComparison.exp)
    modeComparison.medium <- na.omit(modeComparison.medium)
    modeComparison.mut <- na.omit(modeComparison.mut)
    
    modeComparison.all.datasets = rbind(data.frame("dataset" = "exp", modeComparison.exp),
                                        data.frame("dataset" = "mut", modeComparison.mut),
                                        data.frame("dataset" = "medium", modeComparison.medium))
    
    dat = modeComparison.all.datasets
    dat$dataset = revalue(dat$dataset, c("medium" = "selection"))
    dat$dataset = revalue(dat$dataset, c("mut" = "no-selection"))
    
    # calculate dataset of means
    d = ddply(dat,c("dataset", "mode"), function(x) {
        mean.nH = mean(x$nH)
        data.frame(mean.nH)
    })
    
    # calculate correlations berween profiles
    d.exp = subset(d, dataset == "exp")
    d.mut = subset(d, dataset == "no-selection")
    d.medium = subset(d, dataset == "selection")
    
    r.exp.mut = cor(d.exp$mean.nH, d.mut$mean.nH)
    r.exp.medium = cor(d.exp$mean.nH, d.medium$mean.nH)
    r.mut.medium = cor(d.mut$mean.nH, d.medium$mean.nH)
    
    cc.nH.exp.mut[which(input$family == family), which(R0S == R0)] = r.exp.mut
    cc.nH.exp.medium[which(input$family == family), which(R0S == R0)] = r.exp.medium
    cc.nH.mut.medium[which(input$family == family), which(R0S == R0)] = r.mut.medium
    
    # 100 modes
    r.exp.mut = cor(d.exp$mean.nH[1:100], d.mut$mean.nH[1:100])
    r.exp.medium = cor(d.exp$mean.nH[1:100], d.medium$mean.nH[1:100])
    r.mut.medium = cor(d.mut$mean.nH[1:100], d.medium$mean.nH[1:100])
    
    cc.nH.100.exp.mut[which(input$family == family), which(R0S == R0)] = r.exp.mut
    cc.nH.100.exp.medium[which(input$family == family), which(R0S == R0)] = r.exp.medium
    cc.nH.100.mut.medium[which(input$family == family), which(R0S == R0)] = r.mut.medium
  }
}

cc.z.RSD = data.frame(cc.z.RSD.exp.mut, cc.z.RSD.exp.medium, cc.z.RSD.mut.medium)
cc.z.RSD = cc.z.RSD[1:8, ]
cc.z.RSD = rbind(cc.z.RSD, colMeans(cc.z.RSD))
family = as.character(input$family[1:8])
cc.z.RSD_CA = data.frame(c(family, "mean"), cc.z.RSD)
colnames(cc.z.RSD_CA) = c("family", "exp-mut-R0.1", "exp-mut-R0.2", "exp-mut-R0.3", "exp-sel-R0.1", "exp-sel-R0.2", "exp-sel-R0.3", "mut-sel-R0.1", "mut-sel-R0.2", "mut-sel-R0.3")

write.csv(cc.z.RSD_CA, file = "cc_z.RSD_CA.csv")

cc.Pn = data.frame(cc.Pn.exp.mut, cc.Pn.exp.medium, cc.Pn.mut.medium)
cc.Pn = cc.Pn[1:8, ]
cc.Pn = rbind(cc.Pn, colMeans(cc.Pn))
family = as.character(input$family[1:8])
cc.Pn_CA = data.frame(c(family, "mean"), cc.Pn)
colnames(cc.Pn_CA) = c("family", "exp-mut-R0.1", "exp-mut-R0.2", "exp-mut-R0.3", "exp-sel-R0.1", "exp-sel-R0.2", "exp-sel-R0.3", "mut-sel-R0.1", "mut-sel-R0.2", "mut-sel-R0.3")

write.csv(cc.Pn_CA, file = "cc_Pn_CA.csv")

cc.Pn.100 = data.frame(cc.Pn.100.exp.mut, cc.Pn.100.exp.medium, cc.Pn.100.mut.medium)
cc.Pn.100 = cc.Pn.100[1:8, ]
cc.Pn.100 = rbind(cc.Pn.100, colMeans(cc.Pn.100))
family = as.character(input$family[1:8])
cc.Pn.100_CA = data.frame(c(family, "mean"), cc.Pn.100)
colnames(cc.Pn.100_CA) = c("family", "exp-mut-R0.1", "exp-mut-R0.2", "exp-mut-R0.3", "exp-sel-R0.1", "exp-sel-R0.2", "exp-sel-R0.3", "mut-sel-R0.1", "mut-sel-R0.2", "mut-sel-R0.3")

write.csv(cc.Pn.100_CA, file = "cc_Pn_100_CA.csv")

cc.MSF = data.frame(cc.MSF.exp.mut, cc.MSF.exp.medium, cc.MSF.mut.medium)
cc.MSF = cc.MSF[1:8, ]
cc.MSF = rbind(cc.MSF, colMeans(cc.MSF))
family = as.character(input$family[1:8])
cc.MSF_CA = data.frame(c(family, "mean"), cc.MSF)
colnames(cc.MSF_CA) = c("family", "exp-mut-R0.1", "exp-mut-R0.2", "exp-mut-R0.3", "exp-sel-R0.1", "exp-sel-R0.2", "exp-sel-R0.3", "mut-sel-R0.1", "mut-sel-R0.2", "mut-sel-R0.3")

write.csv(cc.MSF_CA, file = "cc_MSF_CA.csv")

cc.nH = data.frame(cc.nH.exp.mut, cc.nH.exp.medium, cc.nH.mut.medium)
cc.nH = cc.nH[1:8, ]
cc.nH = rbind(cc.nH, colMeans(cc.nH))
family = as.character(input$family[1:8])
cc.nH_CA = data.frame(c(family, "mean"), cc.nH)
colnames(cc.nH_CA) = c("family", "exp-mut-R0.1", "exp-mut-R0.2", "exp-mut-R0.3", "exp-sel-R0.1", "exp-sel-R0.2", "exp-sel-R0.3", "mut-sel-R0.1", "mut-sel-R0.2", "mut-sel-R0.3")

write.csv(cc.nH_CA, file = "cc_nH_CA.csv")

cc.nH.100 = data.frame(cc.nH.100.exp.mut, cc.nH.100.exp.medium, cc.nH.100.mut.medium)
cc.nH.100 = cc.nH.100[1:8, ]
cc.nH.100 = rbind(cc.nH.100, colMeans(cc.nH.100))
family = as.character(input$family[1:8])
cc.nH.100_CA = data.frame(c(family, "mean"), cc.nH.100)
colnames(cc.nH.100_CA) = c("family", "exp-mut-R0.1", "exp-mut-R0.2", "exp-mut-R0.3", "exp-sel-R0.1", "exp-sel-R0.2", "exp-sel-R0.3", "mut-sel-R0.1", "mut-sel-R0.2", "mut-sel-R0.3")

write.csv(cc.nH.100_CA, file = "cc_nH_100_CA.csv")