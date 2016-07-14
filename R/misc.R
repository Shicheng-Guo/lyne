
summarizeDat <- function(dat){
#-----------------------------

        #- matrix-type alignemnt data is summarized to sufficient statistics
        #  in dat rows correspond to "species" and columns to genome positions.

        #- initialze result
        ssDat = vector(3,mode="list")
        names(ssDat) = c("cols","mult","freq")
        ssDat$freq = table(dat)/sum(table(dat))

        #- count unique columns
        dat[is.na(dat)] = -Inf
        tmp = as.data.frame((dat))
        tmp = aggregate(rep(1,nrow(tmp)),tmp,sum)

        #- what goes where
        ssDat$cols = t(as.matrix(tmp[,1:(ncol(tmp)-1)]))
        ssDat$cols[is.infinite(ssDat$cols)] = NA
        ssDat$mult = tmp[,ncol(tmp)]

        return(ssDat)
}



