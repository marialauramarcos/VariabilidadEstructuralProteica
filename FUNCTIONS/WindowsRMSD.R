WindowsRMSD <- function(length.windows,
                        r.p.1,
                        r.p.2,
                        aligned.p.1.index.3N,
                        aligned.p.2.index.3N) {
        
  # select aligned sites          
  r.p.1 = r.p.1[aligned.p.1.index.3N]
  r.p.2 = r.p.2[aligned.p.2.index.3N] 

  # calculate the number of aligned sites
  n.aligned = length(r.p.1)/3
  
  # calculate the number of windows and the legth of the last window)
  n.windows = as.integer(n.aligned/length.windows)
  resid.last.window = n.aligned - (n.windows * length.windows)
  
  # sart a loop to aligne each window 
  for (i in (1:n.windows)) {
    
    # set resid value
    if (i == n.windows) {
      resid = resid.last.window
    } else {
      resid = 0
    }
    
    index.windows.i = seq(((length.windows * i) - (length.windows - 1)), ((length.windows * i) + resid))
    index.windows.i.3N = sort(c(index.windows.i * 3, index.windows.i * 3 - 2,  index.windows.i * 3 - 1))

    r.p.2[index.windows.i.3N] = as.vector(fit.xyz(fixed = r.p.1[index.windows.i.3N],
                                          mobile = r.p.2[index.windows.i.3N],
                                          fixed.inds = (1:length(index.windows.i.3N)),
                                          mobile.inds = (1:length(index.windows.i.3N))))
                                        
  }
  r.p.2
}
