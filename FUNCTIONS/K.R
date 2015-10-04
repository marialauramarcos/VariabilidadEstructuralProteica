#Function that calculates K, its Cov (Covariance), va (eigenvalues) and ve (eigenvecors)
#r (equilibrium coordinates of atoms), kij.function (see for example Kij.ANM), R0 (cut- off for ANM) and TOLERANCE (0) must be specified

K <- function(r,kij.function,R0,TOLERANCE){

	nsites <- ncol(r)

	kij <- matrix( 0 , ncol = nsites , nrow = nsites )
	K <- matrix( 0 , ncol = 3*nsites , nrow = 3*nsites )

	for (i in (1:(nsites-1))) {
		ai <- (i+(2*(i-1)))
		bi <- ai + 2
		for (j in ((i+1):nsites)) {
			aj <- (j+(2*(j-1)))				
			bj <- aj +2	
			rij <- r[1:3,j] - r[1:3,i]
			dij <- drop(sqrt(rij %*% rij))
			eij <- rij/dij
			kij[i,j] <- kij.function(dij,R0)
			kij[j,i] <- kij[i,j]					
			K[ai:bi,aj:bj] <- -kij[i,j]*(eij %*% t(eij))
			K[aj:bj,ai:bi] <- K[ai:bi,aj:bj]
      K[ai:bi,ai:bi] <- K[ai:bi,ai:bi] - K[ai:bi,aj:bj]
      K[aj:bj,aj:bj] <- K[aj:bj,aj:bj] - K[aj:bj,ai:bi]
		}	
	}

  eig <- eigen(K , symmetric = TRUE)
  ord <- sort.list(Mod(eig$values) , decreasing = FALSE)
  va <- eig$values[ord]
  ve <- eig$vectors[,ord]
  modes <- va > TOLERANCE
  va <- va[modes]
  ve  <- ve[,modes]
  Cov <-  ve %*% ((1/va)*t(ve))
  output <- list("K"= K,"kij"= kij,"Cov"= Cov,"va"= va,"ve"= ve)
  output
}
