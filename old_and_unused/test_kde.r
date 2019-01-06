library("kde1d")
library("ggplot2")
library("reshape")


cairo_pdf("test_kde.r.out.%d.pdf")


getStatsForProfiles <- function( values, positions )
{
    stats <- data.frame( pos=integer(), iqr0=double(), p25=double(), p50=double(), p75=double(), iqr1=double() )
    outliers <- data.frame( x=double(), y=double() )

    for( rowId in 1:ncol(values))
    {
        xs <- values[,rowId]
        fit <- kde1d( xs )
        quarts <- qkde1d(c(0.25,0.5,0.75), fit)
        iqr <- quarts[3]-quarts[1]
        iqr.vals <- c(quarts[1]-1.5*iqr, quarts[3]+1.5*iqr)

        stats <- rbind( stats, data.frame( pos=positions[rowId], iqr0=iqr.vals[1], p25=quarts[1], p50=quarts[2], p75=quarts[3], iqr1=iqr.vals[2] ) )

        outliers <- c( xs[xs > iqr.vals[2]], xs[xs < iqr.vals[1]] )

        for( i in outliers )
        {
            outliers <- rbind( outliers, data.frame( x=rowId, y=i ) )
        }
    }
    return( list( stats, outliers ) )
}

kdePlot <- function( values, positions )
{
    x <- getStatsForProfiles( values, positions )
    stats    <- x[[1]]
    outliers <- x[[2]]

    print(stats)

    # do plot...
    p <- ggplot( data=stats, aes(x=pos) ) +
        geom_ribbon( aes(ymin=iqr0, ymax=iqr1), color=NA, fill="blue", alpha=0.1 ) +
        geom_ribbon( aes(ymin=p25, ymax=p75),   color=NA, fill="red", alpha=0.5 ) +
        geom_line( aes(y=p50) );
    print(p)
}

num.profiles <- 50
profile.length <- 30

# make random "profiles"
profiles <- array( data=rnorm( num.profiles * profile.length ), dim=c( num.profiles, profile.length ) )
positions <- seq(0, 290, 10)

kdePlot( profiles, positions )

dev.off()
