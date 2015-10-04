mean.groups <- function(ngroups,nexp,df.t){
  
  mean.group <- matrix(0,ncol = ncol(df.t),nrow = ngroups)
  quantile <- matrix(0,ncol = ncol(df.t),nrow = 2)
  
  for (i in (1:ngroups)){
    mean.group[i,] <- colMeans(df.t[((i*nexp)-(nexp-1)):(i*nexp),])
  }
  
  mean.total <- colMeans(mean.group)
  
  for (i in (1:ncol(df.t))){
    quantiles <- as.vector(quantile((mean.group[,i]),probs=c(0.05,0.95)))
    quantile[1,i] <- quantiles[1]
    quantile[2,i] <- quantiles[2]
  }
  
  out <- list("mean.group" = mean.group, "mean.total" = mean.total, "quantile" = quantile)
  out 
}