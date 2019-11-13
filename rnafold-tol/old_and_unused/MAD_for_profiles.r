library("MASS")
library("Gmedian")
library("ggplot2")

#
inputVectors <- read.table("MAD_for_profiles.in.csv", sep=",", header=TRUE ) #, row.names="tax_id")


GmeanAD <- function( data )
{
    #stopifnot(dim(data)[2]==D)

    mean.point = apply(data, 2, mean) # find the mean point
    mean.point.mtx <- matrix(rep(mean.point, each=nrow(data)), nrow=nrow(data)) # duplicate rows
    x3 <- apply( abs(data - mean.point.mtx), 2, mean) # calculate mean absolute deviation from the mean 

    sqrt( sum(x3*x3) ) # return the cartesian length of the means vectors
}



# Calculate the Geometric Median Absolute Deviation of N D-dimensional vectors
# My implementation (tested in test_MAD.r)
# based on: https://en.wikipedia.org/wiki/Median_absolute_deviation#Geometric_median_absolute_deviation
GMAD <- function( data )
{
    #stopifnot(dim(data)[2]==D)
    geomedian <- Gmedian( data, init=apply(data, 2, mean), nstart=100 ) # This seems to ensure finding the right median on test-data
    
    geomedian.mtx <- matrix(rep(geomedian, each=nrow(data)), nrow=nrow(data))

    x3 <- apply( abs(data - geomedian.mtx), 2, median )

    ret <- sqrt(sum(x3*x3))
    
    p <- ggplot( data.frame( y=x3 ), aes(y) ) +
        geom_histogram() +
        geom_vline( xintercept=ret,           color="green" ) +
        geom_vline( xintercept=GmeanAD(data), color="red" )
    ggsave("MAD_for_profiles.out.pdf", plot=p)

    return(ret)
}


# silly late-night python interface (tm)
pythonIOinterface <- function(infile, outfile)
{
    inputVectors <- as.matrix( read.table(infile, sep=",", header=FALSE ) )
    gmad    <- GMAD( inputVectors )
    gmeanad <- GmeanAD( inputVectors )
    write.table( data.frame( N=c(nrow(inputVectors)), dim=c(ncol(inputVectors)), geometric.median.abs.dev=c(gmad), geometric.mean.abs.dev=c(gmeanad) ), file=outfile, sep=",", row.names=FALSE )
}

pythonIOinterface( "MAD_for_profiles.in.csv", "MAD_for_profiles.out.csv" )

