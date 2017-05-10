# CM
# prepare the data
family = cc.z.RSD_CM$family
n.families = length(family)
n.sets = ncol(cc.z.RSD_CM) - 2

cc.z.RSD_CM = cc.z.RSD_CM[1:8, 3:ncol(cc.z.RSD_CM)]

c.R0 = rep(rep(paste("R0.", seq(n.sets/3), sep = ""), 3), 8)
dataset = rep(c("exp vs mut","exp vs mut","exp vs mut", "exp vs selection", "exp vs selection", "exp vs selection", "selection vs no-selection", "selection vs no-selection", "selection vs no-selection"), 8)

df.cc.z.RSD_CM = cbind(data.frame("dataset" = dataset, "c.R0" = c.R0), cc.z.RSD = c(t(cc.z.RSD_CM)))

dat = df.cc.z.RSD_CM

d = ddply(dat,c("dataset", "c.R0"), function(x) {
  cc.mean = mean(x$cc.z.RSD)
  cc.sd = sd(x$cc.z.RSD)
  data.frame(cc.mean, cc.sd)
})

## prepare plots with mean profiles and deviations
ggplot(d, aes(x = c.R0, y = cc.mean, col = dataset, fill = dataset)) +
  geom_bar() +
  geom_errorbar(aes(ymin = cc.mean - cc.sd, ymax = cc.mean + cc.sd), alpha = .4, show.legend = F) +
  labs(x = "R0", y = "mean cc z.RSD CM")


```{r read-files}
## read files CA
cc.z.RSD_CA = read.csv("cc_z.RSD_CA.csv")
cc.Pn_CA = read.csv("cc_Pn_CA.csv")
cc.Pn.100_CA = read.csv("cc_Pn_100_CA.csv")
cc.MSF_CA = read.csv("cc_MSF_CA.csv")
cc.nH_CA = read.csv("cc_nH_CA.csv")
cc.nH.100_CA = read.csv("cc_nH_100_CA.csv")

## read files CM
cc.z.RSD_CM = read.csv("cc_z.RSD_CM.csv")
cc.Pn_CM = read.csv("cc_Pn_CM.csv")
cc.Pn.100_CM = read.csv("cc_Pn_100_CM.csv")
cc.MSF_CM = read.csv("cc_MSF_CM.csv")
cc.nH_CM = read.csv("cc_nH_CM.csv")
cc.nH.100_CM = read.csv("cc_nH_100_CM.csv")
```


FAKE,FAKE,1mcta,LFENM,ANM,10,10,7.5,Keff,G

sh3,small beta,1lcka,LFENM,ANM,10,12.5,10,Keff,F
snakesToxin,small,1ntx,LFENM,ANM,10,12.5,10,Keff,S
globins,alpha,1a4fb,LFENM,ANM,10,12.5,10,Keff,d
az,beta,1bxva,LFENM,ANM,10,12.5,10,Keff,G
phoslip,alpha,1jiaa,LFENM,ANM,10,12.5,10,Keff,L
fabp,beta,1hmt,LFENM,ANM,10,12.5,10,Keff,D
rrm,alpha and beta,1fxla2,LFENM,ANM,10,12.5,10,Keff,G
serinProteases,beta,1mcta,LFENM,ANM,10,12.5,10,Keff,G

FAKE,FAKE,1mcta,LFENM,ANM,10,10,7.5,Keff,G

sh3,small beta,1lcka,LFENM,ANM,10,7.5,5,Keff,F
snakesToxin,small,1ntx,LFENM,ANM,10,7.5,5,Keff,S
globins,alpha,1a4fb,LFENM,ANM,10,7.5,5,Keff,d
az,beta,1bxva,LFENM,ANM,10,7.5,5,Keff,G
phoslip,alpha,1jiaa,LFENM,ANM,10,7.5,5,Keff,L
fabp,beta,1hmt,LFENM,ANM,10,7.5,5,Keff,D
rrm,alpha and beta,1fxla2,LFENM,ANM,10,7.5,5,Keff,G
serinProteases,beta,1mcta,LFENM,ANM,10,7.5,5,Keff,G

sh3,small beta,1lcka,LFENM,ANM,10,10,7.5,Keff,F
snakesToxin,small,1ntx,LFENM,ANM,10,10,7.5,Keff,S
globins,alpha,1a4fb,LFENM,ANM,10,10,7.5,Keff,d
az,beta,1bxva,LFENM,ANM,10,10,7.5,Keff,G
phoslip,alpha,1jiaa,LFENM,ANM,10,10,7.5,Keff,L
fabp,beta,1hmt,LFENM,ANM,10,10,7.5,Keff,D
rrm,alpha and beta,1fxla2,LFENM,ANM,10,10,7.5,Keff,G
serinProteases,beta,1mcta,LFENM,ANM,10,10,7.5,Keff,G

fer4,small,1dura,LFENM,ANM,10,10,7.5,Keff,J
ldh,alpha-beta,1ldm,LFENM,ANM,10,10,7.5,Keff,G
gluts,multi-domain,1hna,LFENM,ANM,10,10,7.5,Keff,A
cys,alpha and beta,1cqda,LFENM,ANM,10,10,7.5,Keff,H

FAKE,FAKE,d,10,ANM,2,1,7.5,TRUE,FALSE,TRUE,TRUE,TRUE,TRUE,Keff
igvar_l,all beta,1reilv,LFENM,ANM,10,10,7.5,Keff,J
igvar_h,all beta,1fdlhv,LFENM,ANM,10,10,7.5,Keff,S
alpha_amylase_N,multi domain,7taa,LFENM,ANM,10,10,7.5,Keff,L
alpha_amylase,alpga beta barrel,1g94a,LFENM,ANM,10,10,7.5,Keff,D
alpha_amylase_C,all beta,1ciu,LFENM,ANM,10,10,7.5,Keff,J

sh3,small beta,1lcka,LFENM,ANM,10,10,7.5,Keff,F
snakesToxin,small,1ntx,LFENM,ANM,10,10,7.5,Keff,S
globins,alpha,1a4fb,LFENM,ANM,10,10,7.5,Keff,d
az,beta,1bxva,LFENM,ANM,10,10,7.5,Keff,G
phoslip,alpha,1jiaa,LFENM,ANM,10,10,7.5,Keff,L
fabp,beta,1hmt,LFENM,ANM,10,10,7.5,Keff,D
rrm,alpha and beta,1fxla2,LFENM,ANM,10,10,7.5,Keff,G
serinProteases,beta,1mcta,LFENM,ANM,10,10,7.5,Keff,G