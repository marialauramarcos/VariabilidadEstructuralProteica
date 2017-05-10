# input
number.R0 = 3

## read files CA
data.dir = paste("OUT/out_subset_CA_ANM", sep = "")

cc.z.RSD_CA = read.csv("cc_z.RSD_CA.csv")
cc.Pn_CA = read.csv("cc_Pn_CA.csv")
cc.Pn.100_CA = read.csv("cc_Pn_100_CA.csv")
cc.MSF_CA = read.csv("cc_MSF_CA.csv")
cc.nH_CA = read.csv("cc_nH_CA.csv")
cc.nH.100_CA = read.csv("cc_nH_100_CA.csv")

## read files CM
data.dir = paste("OUT/out_subset_CM_ANM", sep = "")

cc.z.RSD_CM = read.csv("cc_z.RSD_CM.csv")
cc.Pn_CM = read.csv("cc_Pn_CM.csv")
cc.Pn.100_CM = read.csv("cc_Pn_100_CM.csv")
cc.MSF_CM = read.csv("cc_MSF_CM.csv")
cc.nH_CM = read.csv("cc_nH_CM.csv")
cc.nH.100_CM = read.csv("cc_nH_100_CM.csv")

# prepare data

family = cc.z.RSD_CA$family
n.families = length(family)
n.sets = ncol(cc.z.RSD_CA) - 2

cc.z.RSD_CA = cc.z.RSD_CA[, 3:ncol(cc.z.RSD_CA)]


c.R0 = rep(rep(paste("R0.", seq(n.sets/3), sep = ""), 3), 9)
dataset = rep(c("exp vs mut","exp vs mut","exp vs mut", "exp vs selection", "exp vs selection", "exp vs selection", "selection vs no-selection", "selection vs no-selection", "selection vs no-selection"), 9)

df.cc.z.RSD_CA = cbind(data.frame("dataset" = dataset, "c.R0" = c.R0), cc.z.RSD = c(t(cc.z.RSD_CA)))

dat = df.cc.z.RSD_CA

d = ddply(dat,c("dataset", "c.R0"), function(x) {
    cc.mean = mean(x$cc.z.RSD)
    cc.sd = sd(x$cc.z.RSD)
    data.frame(cc.mean, cc.sd)
})

## prepare plots with mean profiles and deviations
p = ggplot(d, aes(x = c.R0, y = cc.mean, col = dataset, fill = dataset)) +
  geom_point() +
  geom_errorbar(aes(ymin = cc.mean - cc.sd, ymax = cc.mean + cc.sd), alpha = .4, show.legend = F) +
  labs(x = "R0", y = "mean cc")
