library("ggplot2")
library("nlme")

pdf("test_gibbs_and_mahalanobis.out.pdf")

# Source: http://blog.revolutionanalytics.com/2016/08/simulating-form-the-bivariate-normal-distribution-in-r-1.html
gibbs <- function(n, mu1, s1, mu2, s2, rho)
{
    mat <- matrix( ncol=2, nrow=n )
    x <- 0
    y <- 0
    mat[1, ] <- c(x,y)
    for( i in 2:n)
    {
        x <- rnorm( 1, mu1 +
            (s1/s2) * rho * (y - mu2), sqrt((1 - rho^2)*s1^2))
        y <- rnorm( 1, mu2 +
            (s2/s1) * rho * (y - mu1), sqrt((1 - rho^2)*s2^2))
        mat[i, ] <- c(x,y)
    }
    return( mat)
}

linearrand <- function(n, intercept, slope, mu)
{
    mat <- matrix( ncol=2, nrow=n )
    mat[,1] <- runif(n, min=-5, max=5)

    mat[,2] <- mat[,1] * slope + intercept + rnorm(n, 0, mu)
    
    return( mat )
}


#vcv <- matrix( data=c(1.0, 0.52, 0.52, 1.5), nrow=2, ncol=2)
vcv <- matrix( data=c(1.0, 0.35, 0.35, 2.5), nrow=2, ncol=2)
mu <- c( 0, 0 )

bvn4 <- gibbs( 200000, mu[1], vcv[1,1], mu[2], vcv[2,2], vcv[1,2] )
#bvn4 <- linearrand( 200000, 3, 2, 2 )
colnames(bvn4) <- c("bvn4_X1", "bvn4_X2")


N <- 500
bvn4_s <- data.frame(bvn4[sample(200000, N),])

df1 <- data.frame(bvn4_s)
lm1 <- lm( bvn4_X2 ~ bvn4_X1, df1 )
co1 <- coef(lm1)
summary(lm1)

#gls1 <- gls( as.formula("bvn4_X2 ~ bvn4_X1"), df1, correlation=vcv); summary(gls1)
gls1 <- gls( bvn4_X2 ~ bvn4_X1, df1 ); summary(gls1)
#cs1 <- corAR1( 0.2, form=~bvn4_X1 )
#gls1 <- gls( bvn4_X2 ~ bvn4_X1, df1, correlation=corSymm( ~bvn4_X1 )); summary(gls1)
co2 <- coef(gls1)
print(co2)


cov1 <- cov(bvn4_s)
mu1 <- colMeans(bvn4_s)


design <- cbind( rep(1, N), bvn4_s$bvn4_X1 )
print(dim(design))
print(dim(design %*% c(co2["(Intercept)"], co2["bvn4_X1"]) ) )
print(dim(bvn4_s$bvn_X2))
resid <- bvn4_s$bvn4_X2 - design %*% c(co2["(Intercept)"], co2["bvn4_X1"])
print(dim(resid))

#print(mahalanobis( c(0, 0), mu1, cov1 ) )
#print(mahalanobis( c(co2["(Intercept)"], co2["bvn4_X1"]), mu1, cov1 ) )
print(mahalanobis( resid, mu1, cov1 ) )



p <- ggplot( data=data.frame(bvn4), aes(x=bvn4_X1, y=bvn4_X2) ) +
    stat_density_2d(geom="raster", aes(fill=..density..), contour=FALSE) +
    geom_point( data=bvn4_s, aes(x=bvn4_X1, y=bvn4_X2), colour="yellow", size=0.5, alpha=0.5 ) +
    geom_abline( slope=co1["bvn4_X1"],  intercept=co1["(Intercept)"], colour="red"  ) +
    geom_abline( slope=co2["bvn4_X1"],  intercept=co2["(Intercept)"], colour="green"  ) +
    coord_fixed(ratio=1) +
    theme( aspect.ratio=1 )
print(p)

dev.off()
