
score.fname = file.path(out.dir, paste(p.ref, "_consurf", sep = ""))
sim.data.no.sel.fname = file.path(out.dir, paste(family, "_R0_", R0, "_beta_no.sel_out_df.data.csv", sep = ""))
sim.data.strong.sel.fname = file.path(out.dir, paste(family, "_R0_", R0, "_beta_strong.sel_out_df.data.csv", sep = ""))

score = read.csv(score.fname)
score = consurf$score

sim.data.no.sel = read.csv(sim.data.no.sel.fname)
mut.site.no.sel = sim.data.no.sel$site
p.accept.no.sel = sim.data.no.sel$p.accept

sim.data.strong.sel = read.csv(sim.data.strong.sel.fname)
mut.site.strong.sel = sim.data.strong.sel$site
p.accept.strong.sel = sim.data.strong.sel$p.accept

mean.p.accept.no.sel = c()
mean.p.accept.strong.sel = c()

for (i in (1:n.sites)) {
  p.accept.no.sel.i = p.accept.no.sel[which(mut.site.no.sel == i)]
  mean.p.accept.no.sel.i = mean(p.accept.no.sel.i)
  mean.p.accept.no.sel[i] = mean.p.accept.no.sel.i
  
  p.accept.strong.sel.i = p.accept.strong.sel[which(mut.site.strong.sel == i)]
  mean.p.accept.strong.sel.i = mean(p.accept.strong.sel.i)
  mean.p.accept.strong.sel[i] = mean.p.accept.strong.sel.i
}

plot(x = mean.p.accept.no.sel, y = CN)
plot(x = mean.p.accept.no.sel, y = score)
plot(x = score, y = CN)

r.mean.p.accept.no.sel.CN = cor(mean.p.accept.no.sel, CN)
r.mean.p.accept.no.sel.score = cor(mean.p.accept.no.sel, score)
r.score.CN = cor(score, CN)

cor.df = data.frame("comparison" = c("p.accept.no.sel vs CN", "p.accept.no.sel vs score", "score vs CN"),
                    "R" = c(r.mean.p.accept.no.sel.CN , r.mean.p.accept.no.sel.score, r.score.CN))

knitr::kable(cor.df, digits = 2, caption = "Correlation between sequence profiles")


plot(x = mean.p.accept.strong.sel, y = CN)
plot(x = mean.p.accept.strong.sel, y = score)
plot(x = score, y = CN)

r.mean.p.accept.strong.sel.CN = cor(mean.p.accept.strong.sel, CN)
r.mean.p.accept.strong.sel.score = cor(mean.p.accept.strong.sel, score)
r.score.CN = cor(score, CN)

cor.df = data.frame("comparison" = c("p.accept.strong.sel vs CN", "p.accept.strong.sel vs score", "score vs CN"),
                    "R" = c(r.mean.p.accept.strong.sel.CN , r.mean.p.accept.strong.sel.score, r.score.CN))

knitr::kable(cor.df, digits = 2, caption = "Correlation between sequence profiles")
