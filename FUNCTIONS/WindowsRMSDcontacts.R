WindowsRMSDcontacts <- function(kij,
                                r.p.1,
                                r.p.2,
                                aligned.p.1.index,
                                aligned.p.2.index) {
  # prepare r.p.1 and r.p.2
  r.p.1 = matrix(r.p.1, nrow = 3)[, aligned.p.1.index]
  r.p.2 = matrix(r.p.2, nrow = 3)[, aligned.p.2.index]
  
  # change kij
  diag(kij) = 1
  
  # calculate the number of windows
  n.windows = length(aligned.p.1.index)

  # build a vector to save the data
  score.windows = c()
  
  # sart a loop to aligne each window and calculate the score
  for (i in (1:n.windows)) {
    
    index.windows.i = which(kij[i, ] == 1)
    index.windows.i.3N = sort(c(index.windows.i * 3, index.windows.i * 3 - 2,  index.windows.i * 3 - 1))

    r.p.2.window = matrix(fit.xyz(fixed = as.vector(r.p.1),
                                 mobile = as.vector(r.p.2),
                             fixed.inds = index.windows.i.3N,
                            mobile.inds = index.windows.i.3N)[index.windows.i.3N], nrow = 3)
    
    score.windows.i = mean(colSums((r.p.2.window - r.p.1[, index.windows.i]) ^ 2)) 
    score.windows[i] = score.windows.i
  }
  as.vector(score.windows)
}
