library("lmtest")
library("ggplot2")

#N <- 2000

test.N <- function(N)
{

    # Create heteroskedistic sample
    a1 <- data.frame( x=runif(N/2, min= 0, max= 50), y=rnorm(N/2, mean=0.0, sd=1.0))
    a2 <- data.frame( x=runif(N/2, min=50, max=100), y=rnorm(N/2, mean=2.0, sd=5.0))
    sample.pos <- rbind(a1, a2)
    lm.pos <- lm(data=sample.pos, y~x)

    # Create homoskedistic sample
    b1 <- data.frame( x=runif(N/2, min= 0, max= 50), y=rnorm(N/2, mean=0.0, sd=1.0))
    b2 <- data.frame( x=runif(N/2, min=50, max=100), y=rnorm(N/2, mean=2.0, sd=1.0))
    sample.neg <- rbind(b1, b2)
    lm.neg <- lm(data=sample.neg, y~x)

    bptest.pos <- bptest( lm.pos )
    bptest.neg <- bptest( lm.neg )

    return( c(test=bptest.pos$p.value, control=bptest.neg$p.value))
}

test.data <- data.frame(N=integer(), test=double(), control=double())

for( i in 1:10000 )
{
    N <- round( runif(1, min=10, max=250) )
    
    r <- test.N(N)

    test.data <- rbind( test.data, data.frame( N=c(N), test=c(r[1]), control=c(r[2]) ) )
}

    
#print(test.data)
p <- ggplot( test.data, aes(x=N) ) +
    geom_point( aes(y=test), color="red", size=0.5 ) +
    geom_point( aes(y=control), color="green", size=0.5 )
print(p)

dev.off()
quit()

