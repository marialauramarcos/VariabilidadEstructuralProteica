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
  mean.group <- matrix(0, ncol = ncol(data), nrow = ngroups)
  quantile <- matrix(0, ncol = ncol(data), nrow = 2)
  for (i in (1:ngroups)) {
    mean.group[i, ] <- colMeans(data[((i * length) - (length - 1)):(i * length), ])
  }
  mean.total <- colMeans(mean.group)
  for (i in (1:ncol(data))){
    quantiles <- as.vector(quantile((mean.group[, i]), probs = c(0.05, 0.95)))
    quantile[1, i] <- quantiles[1]
    quantile[2, i] <- quantiles[2]
  }
  output <- list("mean.group" = mean.group, "mean.total" = mean.total, "quantile" = quantile)
  output 
}