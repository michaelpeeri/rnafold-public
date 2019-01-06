library("ggplot2")
library("reshape")

pdf("brownian_motion.pdf")

sim_bm <- function(size, sigma=1)
{
    y <- rnorm(size-1, 0.0, sqrt(sigma))
    #y[0] = 0.0
    #x = y

    #for (i in 1:size)
    #{
    #    #x[i] <- 1/sqrt(size) * sum(y[1:i]*sqrt(i))
    #    x[i] <- sum(y[1:i])
    #}
    x <- c(0.0, cumsum(y))
    
    return(x)
}

N <- 2000
T <- 150
#f = data.frame(bm1=double(), bm2=double(), bm3=double(), bm4=double(), bm5=double())
df = data.frame(time=1:T, row.names=1:T)

for (i in 1:N)
{
    x <- data.frame( value=sim_bm(T, 0.05), row.names=1:T)
    #plot(x, type="l", xlab="Time", ylab="Value")
    df <- cbind( df, x )
}

#print(dim(df))
#print(dimnames(df))
#dimnames(df)[[2]] <- list("time", "bm1", "bm2", "bm3", "bm4", "bm5")

bmlist <- vapply(1:N, function(n) paste("bm", as.character(n), sep=""), character(1))

dimnames(df)[[2]] <- c("time", bmlist)

#print(dimnames(df)[[2]])
#print(bmlist)
#print(df[bmlist])

df$var = apply( df[bmlist], 1, var )



#print(df)

theme_set(theme_gray(base_size=18))

df2 <- data.frame( melt( df, id.vars=c("time"), measure.vars=bmlist ) )
#print(df2)
p <- ggplot( data=df2, aes(x=time, y=value, colour=variable) ) +
    scale_colour_grey(start=0.0, end=0.0) +
    geom_line(alpha=0.02) +
    annotate( "text", x=Inf,  y=Inf, label=sprintf("n = %d", N), hjust=1, vjust=1, size=6 ) +
    guides(colour=FALSE); p

p <- ggplot( data=df, aes(x=time, y=var) ) +
    labs(y="Variance") +
    geom_line() +
    annotate( "text", x=Inf,  y=Inf, label=sprintf("n = %d", N), hjust=1, vjust=1, size=6 ) +
    guides(colour=FALSE); p


dev.off()
