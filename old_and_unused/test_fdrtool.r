library('fdrtool')

cairo_pdf("test_fdrtool.r.out.pdf")

interpolate <- function(x, n)
{
    xi <- rep(x, each=n)
    interp <- filter(xi, rep(1/n, n), sides=2)
    return( interp[!is.na(interp)] )
}

center.points <- c(0.0, runif(6, min=-10, max=10))   # create random points to interpolate
centers <- interpolate( center.points, n=5 )   # interpolate the points (to create dependency)

# draw samples for each point
num.groups <- length(centers)
num.samples <- 100  # how many samples in each group
samples <- matrix(rnorm(num.groups*num.samples, mean=rep(centers, num.samples)), nrow=num.groups)  # arrange samples points as a matrix
pval.samples  <- apply( samples, 1, function (x) { t.test(x)$p.value } )
tstat.samples <- apply( samples, 1, function (x) { t.test(x)$statistic } )
#samples.mean = apply(samples, 1, mean)
#samples.sd = apply(samples, 1, sd)

#fdrtool( pval.samples, statistic="pvalue", plot=TRUE )
fdrtool( tstat.samples, statistic="normal", plot=TRUE )

fdr.pval <- p.adjust( pval.samples, method="fdr" )
print( cbind(centers, pval.samples, tstat.samples, fdr.pval) )


dev.off()
quit()

