#Function that calculates the force that models a mutantion using LF-ENM
#i (mutated site), r(equilibrium coordinates of atoms), kij (kij between sites i and j) and fmax (maximum value for f) must be specified

force = function(i,r,kij,fmax){
    nsites = ncol(r)
    fvec = matrix(0,3,nsites)
    for (j in seq(nsites)){
        if(kij[i,j] > 0.9) {
            rij = r[1:3,j] - r[1:3,i]
            eij = rij/sqrt(sum(rij^2))
            fij = fmax*sqrt(3)*runif(1,-1,1)
            fvec[,j] = fij*eij
            fvec[,i] = fvec[,i] - fvec[,j]
        }
    }
    dim(fvec)=c(3*nsites)
    fvec
}



