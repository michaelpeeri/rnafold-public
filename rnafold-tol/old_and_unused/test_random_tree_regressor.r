library("randomForest")
library("ggplot2")

# Initialize random seed
seed <- as.integer(difftime( as.POSIXlt(Sys.time()), as.POSIXlt("2010-01-01 12:00:00"), units="secs" ))
set.seed( seed )
print(sprintf("Using random seed %d", seed))

pdf("test_random_tree_regressor.out.pdf")

poly <- function(x, coeffs)
{
    return( sum( c(1,x) * coeffs ) + 10 )  # 1st degree polynomial
}

makeDummyVar <- function(N, K)
{
    return( sample( vapply(1:K, function (n) paste("L", n, sep=""), character(1)), size=N, replace=TRUE ) )
}

Nx <- 1000
xs <- seq(-5, 5, len=Nx)
coeff1 <- runif(2, min=-1, max=1)
ys1 <- vapply( xs, function (x) poly(x, coeff1), numeric(1) )

coeff2 <- runif(2, min=-1, max=1)
ys2 <- vapply( xs, function (x) poly(x, coeff2), numeric(1) )

coeff3 <- runif(2, min=-1, max=1)
ys3 <- vapply( xs, function (x) poly(x, coeff3), numeric(1) )

data <- data.frame( x=rep(xs,3), y=c(ys1, ys2, ys3), cat=factor(c(rep("A", Nx), rep("B", Nx), rep("C", Nx))), dummy1=makeDummyVar(Nx*3,5), dummy2=makeDummyVar(Nx*3,20), dummy3=makeDummyVar(Nx*3,5) )
perm <- sample(Nx*3, replace=FALSE)
data <- data[perm,]  # permute the elements

print(data[data$cat=="A","x"])
print(data[data$cat=="A","y"])

p <- ggplot( data=data, aes(x=x, y=y, color=cat) ) + geom_point()
print(p)

#rf <- randomForest( y ~ x + cat + dummy1 + dummy2 + dummy3 + 1, data=data )
rf <- randomForest( y ~ x + cat + dummy1 + dummy2 + dummy3, data=data, nodesize=0.2, mtry=5, ntree=100 )
varImpPlot( rf )
print( importance( rf ) )


plotPrediction <- function( cat )
{
    out <- predict( rf, data[data$cat==cat,] )

    comp <- data.frame( x=data[data$cat==cat, "x"], actual=data[data$cat==cat,"y"], predicted=out )
    p <- ggplot( data=comp, aes(x=x) ) +
        geom_point( aes(y=actual),    color="black" ) +
        geom_point( aes(y=predicted), color="green" )
    print(p)

}


plotPrediction( "A" )
plotPrediction( "B" )
plotPrediction( "C" )

dev.off()



