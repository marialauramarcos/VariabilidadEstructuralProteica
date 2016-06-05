library("ggplot2")

data = data.frame("x.c" = c(1, 2, 3),
                  "y.c" = c(10, 20, 30),
                  "s.d" = c(0.5, 0.8, 1.4),
                  "names" = c(1,2,3))

ggplot(data, aes(x = names)) +
geom_point(aes(y = x.c, colour = "x.c")) + 
geom_point(aes(y = y.c, colour = "y.c")) + 
geom_errorbar(aes(y = y.c, ymin = y.c - s.d, ymax = y.c + s.d), width = 0, col = "blue") +
ylim(c(0,40))

library("ggplot2")

data = data.frame("x.c" = c(1, 2, 3),
                  "y.c" = c(10, 20, 30),
                  "s.d" = c(0.5, 0.8, 1.4),
                  "names" = c(1, 2, 3))

MSDi.exp = matrix(runif(151, 1, 151), 151, 1)
MSDi.theo = matrix(runif(151, 1, 151), 151, 1)
error = runif(151, 0, 1)
sites = seq(1:151)

data.MSDi = data.frame("MSDi.exp" = MSDi.exp, "MSDi.theo" = MSDi.theo, "sites" = sites, "error" = error)
ggplot(data = data.MSDi, aes(x = sites)) + 
  geom_point(aes(y = MSDi.exp, colour = "MSDi.exp")) + 
  geom_line(aes(y = MSDi.theo, colour = "MSDi.theo"), col = "black") +
  ylab("MSDi.theo") + 
  geom_errorbar(aes(y = MSDi.exp, ymin = MSDi.exp - error, ymax = MSDi.exp + error), width = 0, col = "red") 