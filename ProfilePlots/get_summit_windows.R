get_summit_windows<-function(bi,W=10,numBins=1) {

    X2=c()
    X3=c()

    for(ni in seq(numBins)) {
        X2=c(X2,-(ni+1/2)*W)
        X3=c(X3,-(ni-1/2)*W)
    }

    X2=c(X2,-W*(1/2))
    X3=c(X3,W*(1/2)+1)

    for(ni in seq(numBins)) {
        X2=c(X2,(ni-1/2)*W+1)
        X3=c(X3,(ni+1/2)*W+1)
    }

    tibble(X1=bi$X1,X2=bi$X2+X2,X3=bi$X2+X3,X4=bi$X4,X5=bi$X5) %>%
        arrange(X2) %>%
        mutate(II=row_number())

}

