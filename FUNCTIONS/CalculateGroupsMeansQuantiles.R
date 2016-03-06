# This function partitionates data into "ngroups" of "length" members and calculates
# means of means and quantiles of means
#
#  Args:
#    data: The matrix of data to partitionate. nrow(data) must be ngroups*length.
#    ngroups: Number of groups to partitionate the data.
#    length: Number of members for each group. 
#
#  Returns:
#    mean.group: ngroups x ncol(data) matrix with the mean of each group in each row.
#    mean.total: 1 x ncol(data) matrix with colMeans(mean.groups).
#    quantile: 2 x ncol(data) matrix with quantile 0.05 of each column of mean.groups in the first row 
#    and quantile 0.95 of each column of mean.groups in the second row.
CalculateGroupsMeansQuantiles <- function(data, ngroups, length) {
  mean.group = matrix(0, ncol = ncol(data), nrow = ngroups)
  quantile1.group = matrix(0, ncol = ncol(data), nrow = ngroups)
  quantile2.group = matrix(0, ncol = ncol(data), nrow = ngroups)
  
  for (i in (1:ngroups)) {
    mean.group[i, ] = colMeans(data[((i * length) - (length - 1)):(i * length), ])
    for (j in (1:ncol(data))){
      quantiles = as.vector(quantile((data[((i * length) - (length - 1)):(i * length), j]), probs = c(0.05, 0.95)))
      quantile1.group[i, j] = quantiles[1]
      quantile2.group[i, j] = quantiles[2]
    }
  }
  mean.mean.group = colMeans(mean.group)
  mean.quantile1.group = colMeans(quantile1.group)
  mean.quantile2.group = colMeans(quantile2.group)
  
  
  output <- list("mean.group" = mean.group, 
                 "mean.mean.group" = mean.mean.group, 
                 "mean.quantile1.group" = mean.quantile1.group, 
                 "mean.quantile2.group" = mean.quantile2.group)
  output 
}