library("MASS")
library("Gmedian")
library("ggplot2")

sigma <- 100.0
D <-2   # The number of dimensions
Nsamples <- 1000

generateSample <- function( N )
{
    r <- mvrnorm( n=N, mu=runif(D), Sigma=diag(runif(D))*runif(1, min=0.5, max=3.0)*sigma )
    r[1:10,1] <- sigma*-20 # create an outlier
    r[1:10,2] <- sigma*20
    return(r)
}

MAD <- function( data )
{
    stopifnot(dim(data)[2]==D)
    #print(dim(data))
    geomedian <- Gmedian( data, init=apply(data, 2, mean), nstart=100 )
    print(geomedian)
    print(median( data ) )
                                        #print(dim(geomedian))
    geomedian.mtx <- matrix(rep(geomedian, each=nrow(data)), nrow=nrow(data))

    x3 <- apply( abs(data - geomedian.mtx), 2, median )
    stopifnot(length(x3)==D)
    
    sqrt(sum(x3*x3))
}

meanAD <- function( data )
{
    stopifnot(dim(data)[2]==D)

    mean.point = apply(data, 2, mean)
    mean.point.mtx <- matrix(rep(mean.point, each=nrow(data)), nrow=nrow(data))
    x3 <- apply( abs(data - mean.point.mtx), 2, mean)

    sqrt( sum(x3*x3) )
}

euclideanDists <- function( data )
{
    stopifnot(dim(data)[2]==D)

    mean.point = apply(data, 2, mean)
    mean.point.mtx <- matrix(rep(mean.point, each=nrow(data)), nrow=nrow(data))
    abs.devs <- abs(data - mean.point.mtx)
    #print(dim(abs.devs))

    x3 <- matrix( apply( abs.devs*abs.devs, 1, sum ), ncol=1)
    #print(x3)
    #print(dim(x3))
    x4 <- apply( x3, 1, sqrt)
    #print(length(x4))

    x4
}


MAD1D <- function( vec )
{
    dists <- abs( vec - median(vec) )
    median( dists )
}

make_circle <- function(cx, cy, rad)
{
    angles = seq(0, 2*pi, 0.1)
    data.frame( x=cx+cos(angles)*rad, y=cy+sin(angles)*rad )
}



testUsingN <- function( N )
{
    td <- generateSample( N )
    #print( median(td))

    mad <- MAD( td )
    print( mad )

    mnad <- meanAD( td )
    print( mnad )
    
    if( D==1 ) { print( MAD1D( td )) }
    else { print( apply( td, 2, var ) ) }

    td.df <- data.frame( td )
    colnames(td.df) <- letters[1:D]

    center <- apply( td, 2, mean )
    geomedian <- Gmedian( td, init=apply(td, 2, mean), nstart=100 )


    p <- ggplot( td.df, aes(x=a,y=b) ) +
        geom_point(alpha=0.4) +
        geom_polygon( data=make_circle(geomedian[1], geomedian[2], mad/2),  aes(x=x,y=y), alpha=0.8, color="green", fill=NA ) +
        geom_polygon( data=make_circle(center[1], center[2], mnad/2), aes(x=x,y=y), alpha=0.8, color="red",   fill=NA ) +
        scale_x_continuous( limits=c(-200,50) ) +
        scale_y_continuous( limits=c(-50,200) ) 
    print(p)


    euc.dists = euclideanDists(td)
    p <- ggplot( data.frame(x=euc.dists), aes(x) ) +
        geom_histogram() +
        geom_vline(xintercept=mad,  color="green") +
        geom_vline(xintercept=mnad, color="red")
    print(p)
    
}

#print( generateSample( 10 ) )
testUsingN(Nsamples)
testUsingN(Nsamples)
testUsingN(Nsamples)
testUsingN(Nsamples)
testUsingN(Nsamples)

#ds <- data.frame( n=integer(), 

dev.off()
quit()
