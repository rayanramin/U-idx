## Uracilation index 
## IMPORTING THE NUCLEOTIDE COUNTS FILES ( USING bam-readcount , SEE README.md ) ##

BRC <- read.table("out.txt", sep="\t", col.names=c("Chr.","POS","REF","COV","A","C","G","T"))

## UX function for calculating Uracilation index at TC context ##
UX_TC <- function(Y) {    
    X <- Y
# defining conditions    
    cond1 <- logical(nrow(X))
    cond2 <- logical(nrow(X))
    for (i in 2:nrow(X)) { cond1[i] <- (( X$REF[i] == "C") & (X$REF[i-1] == "T")) }
    for (i in 1:(nrow(X)-1)) { cond2[i] <- (( X$REF[i] == "G") & (X$REF[i+1] == "A")) }
# TC to TT changes only
    t <- ifelse(cond1, as.integer(X$T), 0) / X$COV
    nc <- sum(cond1)
    AllC <- sum(t) / nc * 1e3
# GA to AA changes only
    a <- ifelse(cond2, as.integer(X$A), 0) / X$COV
    ng <- sum(cond2)
    AllG <- sum(a) / ng * 1e3
# All changes
    u <- (a+t) 
    z <- sum(cond2 | cond1)
    AllU <- sum(u) / z * 1e3
 
    UDF <-  data.frame("POS" = X$POS , "CtoT" = t ,  "GtoA" = a , "CtoT_GtoA" = u)

    return(list(CtoT_indx = AllC , GtoA_indx = AllG ,  U_indx = AllU , Number_of_TCs = nc, Number_of_GAs = ng , Number_of_TC_GA = z,  Out = UDF  ))
}


## UX function for calculating Uracilation index at all cytosines ##
UX_NC <- function(Y) {
    X <- Y
# defining conditions    
    cond1 <- logical(nrow(X))
    cond2 <- logical(nrow(X))
    for (i in 1:nrow(X)) { cond1[i] <- ( X$REF[i] == "C")  } 
    for (i in 1:nrow(X)) { cond2[i] <- ( X$REF[i] == "G")  }
# C to T changes only
    t <- ifelse(cond1, as.integer(X$T), 0) / X$COV
    nc <- sum(cond1)
    AllC <- sum(t) / nc * 1e3
# G to A changes only
    a <- ifelse(cond2, as.integer(X$A), 0) / X$COV
    ng <- sum(cond2)
    AllG <- sum(a) / ng * 1e3
# All changes
    u <- (a+t) 
    z <- sum(cond2 | cond1)
    AllU <- sum(u) / z * 1e3
 
    UDF <-  data.frame("POS" = X$POS , "CtoT" = t ,  "GtoA" = a , "CtoT_GtoA" = u)

    return(list(CtoT_indx = AllC , GtoA_indx = AllG ,  U_indx = AllU , Number_of_Cs = nc, Number_of_Gs = ng , Number_of_C_G = z,  Out = UDF  ))
}

UX_TC(BRC) # for Uracilation index calculations at TC context
UX_NC(BRC) # for Uracilation index calculations for all cytosines
