#Function that analyzes the alignment of 2 proteins in terms of aligned and not aligned sites
#The alignment of P1 and P2 (character string) and nsites of P1 must be specified

alignment <- function(alignment.P1,alignment.P2,nsites.P1){
	lalign = length(alignment.P1)
	index.P1 = seq(1,nsites.P1)
	
	align.P1.index.tot = matrix(ncol = lalign , nrow = 2)
	align.P1.index.tot[1,] = alignment.P1

	count = 1
	for (i in (1:lalign)) {
		if (alignment.P1[i] != "-"){
			align.P1.index.tot[2,i] <- index.P1[count]
			count <- count +1
		}
	}	

	align.P1.index = matrix(ncol=1,nrow=1)
	no.align.P1.index = matrix(ncol=1,nrow=1)
	
	for (i in (1:lalign)) {
		if ((alignment.P2[i] != "-") & (alignment.P1[i] != "-")){
			align.P1.index <- cbind(align.P1.index , align.P1.index.tot[2,i])
		}
		if ((alignment.P2[i] == "-") & (alignment.P1[i] != "-")){
		  no.align.P1.index <- cbind(no.align.P1.index , align.P1.index.tot[2,i])
		}
	}	
	
	align.P1.index = matrix(as.numeric(align.P1.index[-1]) , nrow=1)
	no.align.P1.index = matrix(as.numeric(no.align.P1.index[-1]) , nrow=1)	
	
	nalign <- ncol(align.P1.index)

output = list("align.index" = align.P1.index,"no.align.index" = no.align.P1.index,"nalign" = nalign)
output
}


