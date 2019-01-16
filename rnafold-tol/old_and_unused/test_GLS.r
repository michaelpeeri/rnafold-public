library("nlme")
library("ggplot2")

N <- 30
separation <- 100.0
sigma.int <- 0.05

x.1 <- c( rnorm(N, mean=10.0), rnorm(N, mean=20.0) )
y.1 <- c( rnorm(N, mean= 0.0), rnorm(N, mean=10.0) ) + x.1*0.1
t.1 <- as.factor(c( rep(1, N), rep(2,N) ))
pos.1 <- c( rnorm(N, mean=-separation, sd=sigma.int), rnorm(N, mean=separation, sd=sigma.int) )

d.1 <- data.frame(x=x.1, y=y.1, pos=pos.1, t=t.1)



x.2 <- c( rnorm(N, mean=10.0), rnorm(N, mean=20.0) )
y.2 <- c( rnorm(N, mean= 0.0), rnorm(N, mean=10.0) ) + x.2*0.1
t.2 <- sample( c(-separation,separation), size=N*2, replace=TRUE )
pos.2 <- t.2 + rnorm(N*2, mean=0.0, sd=sigma.int)
t.2 <- as.factor(t.2)

d.2 <- data.frame(x=x.2, y=y.2, pos=pos.2, t=t.2)


make.data.d3 <- function()
{
    
    x.3 <- c( rnorm(N, mean=10.0), rnorm(N, mean=10.0) )
    y.3 <- c( rnorm(N, mean=10.0), rnorm(N, mean=10.0) ) + x.3*0.3
    t.3 <- as.factor(c( rep(1, N), rep(2,N) ))
    pos.3 <- c( rnorm(N, mean=-separation, sd=sigma.int), rnorm(N, mean=separation, sd=sigma.int) )
    
    d.3 <- data.frame(x=x.3, y=y.3, pos=pos.3, t=t.3)
    return(d.3)
}
d.3 <- make.data.d3()

## m <- matrix(0, 16, 16)
## m[1:8, 1:8]  <- 1
## m[9:16,9:16] <- 1
## diag(m) <- 2

#corStruct <- pdSymm(m)

plotData <- function(data)
{
    
    gls1 <- gls( y~x, data, cor=corSpatial(form=~pos) )
    print(summary(gls1))
    co2 <- coef( gls1 )

    lm1 <- lm( y~x, data )
    co <- coef( lm1 )

    p <- ggplot( data, aes(x=x, y=y) ) +
        geom_point( aes(color=t) ) +
        geom_abline( aes( slope=co["x"],  intercept=co["(Intercept)"]),  colour="red"  ) +
        geom_abline( aes( slope=co2["x"], intercept=co2["(Intercept)"]), colour="green"   )
    print(p)

}

plotData(d.1)
plotData(d.2)
plotData(d.3)

dfc <- data.frame(n=integer(), gls.coef=double(), lm.coef=double() )
for(i in 1:5000)
{
    print(i)
    d <- make.data.d3()

    gls.ok <- TRUE
    
    tryCatch( gls1 <- gls( y~x, d, cor=corSpatial(form=~pos) ), error=function(e) gls.ok <- FALSE )
    if( gls.ok )
    {
        print(summary(gls1))
        co2 <- coef( gls1 )

        lm1 <- lm( y~x, d )
        co <- coef( lm1 )
        
        dfc <- rbind( dfc, data.frame(n=c(i), gls.coef=c( co2["x"] ), lm.coef=c( co["x"] ) ) )
    }
}

p <- ggplot( data=dfc, aes(x=gls.coef, y=lm.coef) ) +
    geom_abline( slope=1, intercept=0, color="blue", alpha=0.4, size=0.5 ) +
    geom_point(size=0.5, alpha=0.4)
print(p)


dev.off()
quit()

            
