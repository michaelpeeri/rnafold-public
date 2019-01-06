library("tseries")
library("StepwiseTest")

#matrix(rnorm(m*n),m,n)

# linear interpolation between 1d points
# x - input vector; n - number of points to interpolate each element of x into
interpolate <- function(x, n)
{
    xi <- rep(x, each=n)
    interp <- filter(xi, rep(1/n, n), sides=2)
    return( interp[!is.na(interp)] )
}

center.points <- runif(3, min=-10, max=10)   # create random points to interpolate
centers <- interpolate( center.points, n=5 )   # interpolate the points (to create dependency)

# draw samples for each point
num.groups <- length(centers)
num.samples <- 500  # how many samples in each group
samples <- matrix(rnorm(num.groups*num.samples, mean=rep(centers, num.samples)), nrow=num.groups)  # arrange samples points as a matrix
tstat.samples <- apply( samples, 1, function (x) { t.test(x)$statistic } )
samples.mean = apply(samples, 1, mean)
samples.sd = apply(samples, 1, sd)

# create bootstrap sets
B <- 100
bs <- tsbootstrap(1:num.groups, B, b=2, type="stationary") # indices for bootstrap sets
stopifnot( all( dim(bs)==c(num.groups, B ) ) )

bs.stat <- c()
for( i in 1:B )
{
    samples.bs <- samples[, bs[,i]]
    #print(samples.bs)
    samples.bs.mean <- apply( samples.bs, 1, mean )
    samples.bs.stat <- sqrt(num.samples) * (samples.bs.mean - samples.mean) / samples.sd
    bs.stat <- cbind(bs.stat, samples.bs.stat)
}

print(centers)
FWERkControl( tstat.samples, bs.stat, 1, 0.05 )





