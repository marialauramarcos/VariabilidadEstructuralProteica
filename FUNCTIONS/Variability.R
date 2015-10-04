#Function that calculates measures of structural and dynamical variability between 2 proteins 
#dr (diference between P1 and P2 equilibrium coordinates of atoms), enm.1 and enm.2 (enm containing K, va and ve) must be specified

variability <- function(dr,enm.1,enm.2){
			
	dr.squarei <- (colSums(dr^2))
	ov = t(enm.1$ve) %*% enm.2$ve
  nH = exp(-rowSums(ov^2*log(ov^2+1e-10)))
  nR = 1/rowSums(ov^4)
  d.evalues = enm.2$va - enm.1$va
	Pn = (t(enm.1$ve)%*%matrix(dr[1:3,],ncol=1))^2/sum((t(enm.1$ve)%*%matrix(dr[1:3,],ncol=1))^2)
	Qn = cumsum(Pn)

	output = list(dr.squarei = dr.squarei , nH = nH , nR = nR , d.evalues = d.evalues , Pn = Pn , Qn = Qn)
	output
}