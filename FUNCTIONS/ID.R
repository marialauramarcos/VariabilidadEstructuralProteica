#Function that calculates %ID of 2 proteins
#The alignment of P1 and P2 (a character string) must be specified

ID <- function(alignment.P1, alignment.P2){
  lalign <- length(alignment.P1)
	align.P1 <- matrix(0)
	align.P2 <- matrix(0)

	for (i in (1:lalign)) {
		if ((alignment.P1[i] != "-") & (alignment.P2[i] != "-")){
		  align.P1 <- cbind(align.P1,alignment.P1[i])
			align.P2 <- cbind(align.P2,alignment.P2[i])
		}
	}	

	align.P2 <- align.P2[,-1]	
	align.P1 <- align.P1[,-1]	
	nalign <- length(align.P1)
	
	count = 0
	for (i in (1:nalign)) {
		if (align.P1[i] == align.P2[i]){
			count <- count +1
		}
	}
	ID = count*100/nalign
	ID
}

