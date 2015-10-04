#Function that calculates ANM kij
#dij (distance between atoms i and j at their equilibrium position) and R0 (cut-off) must be specified

kij.function = function(dij,R0) { 
    if (dij <= R0) {
        kij = 1
    } else {
       kij = 0
    }
    kij
}


