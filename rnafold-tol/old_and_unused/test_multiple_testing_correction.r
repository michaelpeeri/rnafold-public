library('ggplot2')


# ------------------- Configuration -------------------
alpha <- 0.05
# --- Random data-set ---
N <- 1000
span <- 5
test.mean <- 10.0
test.sd <- 1.0

# Create normally-distributed auto-correlated series with the given mean and sd
# This implementation has several obvious problems:
# 1) The actual sd of the generated series is somewhat smaller than specified
# 2) In addition to the main autocorrelation, there are additional autocorrelation peaks (higher harmonics)
create.autocorrelated.series <- function(N=10000, autocorr.span=5, mean=0.0, sd=1.0)
{
    data <- rnorm( N/autocorr.span, mean=mean, sd=sd )
    data <- rep( data, 1, each=autocorr.span )  # repeat each element 'autocorr.span' times

    data <- filter( data, rep(1/autocorr.span, autocorr.span) ) # moving average
    data <- data[!is.na(data)]   # remove elements near the end (where the moving average is undefined...)

    #ac <- acf( data, lag.max=100 )  # calculate autocorrelation
    #print(ac$acf)  # strong autocorrelation should exit, and drop off gradually after autocorr.span positions
    #print(mean(data))
    #print(sd(data))
                                        #print(shapiro.test(data))
    #qqnorm(data)
    #print(ks.test(data, "pnorm", mean=mean(data), sd=sd(data)))

    return(data)
}


# Generate random data-set
test.data <- create.autocorrelated.series(N=N, autocorr.span=span, mean=test.mean, sd=test.sd )

Pi_uncorr = abs( qnorm( alpha, mean=mean(test.data), sd=sd(test.data) ) )

# Benjamini–Hochberg procedure
Pi <- sort( pnorm( -abs( (test.data - mean(test.data) )/ sd(test.data) ), mean=0.0, sd=1.0 ) )   # convert to z-scores and get p-values
stopifnot(all(Pi >= 0.0 & Pi <= 1.0))

Neff <- length(Pi) / span
cutoff <- 1:Neff / Neff * alpha
k <- sum(Pi < cutoff)
print(k)
print(k/Neff)
print(Pi[max(1,k)])
Pi_corr <- abs( qnorm( Pi[max(1,k)], mean=mean(test.data), sd=sd(test.data) ) )
print(Pi_corr)


# Diagnostic plot
g <- ggplot( data=data.frame(y=test.data), aes(x=1:length(test.data), y=y) ) +
    geom_line() +
    geom_hline( yintercept = mean(test.data) + Pi_uncorr, color="red",   alpha=0.6 ) +
    geom_hline( yintercept = mean(test.data) - Pi_uncorr, color="red",   alpha=0.6 ) +
    geom_hline( yintercept = mean(test.data) + Pi_corr,   color="green", alpha=0.6 ) +
    geom_hline( yintercept = mean(test.data) - Pi_corr,   color="green", alpha=0.6 ) +
    annotate( "text", x=-Inf, y=Inf, label=sprintf("%4g", sum(test.data > Pi_uncorr) + sum(test.data < -Pi_uncorr) ), hjust=0, vjust=1 )
print(g)

