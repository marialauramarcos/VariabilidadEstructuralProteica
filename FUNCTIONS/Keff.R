#Function that calculates KEFF (K effective), its va (eigenvalues) and ve (eigenvecors)
#r (equilibrium coordinates of atoms), aligned and not aligned indexes, kij.function (see for example Kij.ANM), R0 (cut-off for ANM) and TOLERANCE (0) must be specified

keff <- function(r,align.index,no.align.index,kij.function,R0,TOLERANCE){

	nalign = length(align.index)
	n.no.align = length(no.align.index)
	nsites = ncol(r)
	
	#Calculate K#
	K <- K(r,kij.function,R0,TOLERANCE)$K 

	#Order K in order to get KPP, KQQ, KPQ and KQP and KEFF = KPP - (KPQ%*%CovQQ%*%KQP)#
	sortm <- matrix(0 , ncol = 3*nsites , nrow = 1)
	count = 1
	for (i in (1: nsites)){
		for (j in (1:nalign)){
			if (align.index[j] == i){
				sortm[1,(3*i-2)]<- 3*count-2
				sortm[1,(3*i-1)]<- 3*count-1
				sortm[1,(3*i)]<- 3*count
				count <- count + 1
			}
		}
	}

	if (n.no.align > 0){
		for (i in (1:nsites)){
			for (j in (1:n.no.align)){
				if (no.align.index[j] == i){
					sortm[1,(3*i-2)]<- 3*count-2
					sortm[1,(3*i-1)]<- 3*count-1
					sortm[1,(3*i)]<- 3*count
					count <- count + 1
				}
			}
		}
	}


	Kord <- K[sortm,sortm]

	KPP <- Kord[(1:(3*nalign)),(1:(3*nalign))]

	if(n.no.align > 0){

		KQQ <- Kord[((3*nalign+1):(3*nsites)),((3*nalign+1):(3*nsites))]
		KPQ <- Kord[(1:(3*nalign)),((3*nalign+1):(3*nsites))]
		KQP <- t(KPQ)

		#KQQ^-1#

		eigQQ <- eigen(KQQ,symmetric = T)
		veQQ <- eigQQ$vectors
		vaQQ <- eigQQ$values
		modes <- vaQQ> TOLERANCE #there aren´t negative evalues#
		veQQ  <- veQQ[,modes]
		vaQQ <- vaQQ[modes]
		CovQQ <-  veQQ %*% ((1/vaQQ)*t(veQQ))

		#KEFF#
		KEFF <- KPP - (KPQ%*%CovQQ%*%KQP)
	}

	if(n.no.align == 0){
		KEFF <- KPP
	}

	eigKEFF <- eigen(KEFF,symmetric = T)
	veKEFF <- eigKEFF$vectors
	vaKEFF <- eigKEFF$values
	ord <- sort.list(Mod(vaKEFF),decreasing = FALSE)
	vaKEFF <- vaKEFF[ord]
	veKEFF <- veKEFF[,ord]
	modes <- vaKEFF > TOLERANCE #there aren´t negative evalues#
	veKEFF <- veKEFF[,modes]
	vaKEFF <- vaKEFF[modes]
		
	output = list("KEFF" = KEFF, "ve" = veKEFF, "va" = vaKEFF)
	output
}


