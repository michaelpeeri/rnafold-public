library("ggplot2")
library("scales")

pdf("test_mcfadden_pseudo_r_squared.out.pdf")


McFadden <- function(Pmodel, Pnull)
{
    stopifnot( Pmodel >= 0.0 )
    stopifnot( Pnull  >= 0.0 )
    
    return( 1 - log(Pmodel)/log(Pnull) )
}

CoxSnell <- function(Pmodel, Pnull, N)
{
    stopifnot( Pmodel >= 0.0 )
    stopifnot( Pnull  >= 0.0 )
    
    return(   (1 - (log(Pnull)/log(Pmodel))**(2.0/N) ) / (1 - log(Pnull)**(2.0/N) ) )
}


plotMcFaddenParams <- function( N=80, Lmodel.range=c(0.001, 5.001), Lnull.range=c(0.001, 5.001) )
{
    #N <- 80
    d <- matrix(rep(0, N*N*3), ncol=3)
    ix <- 0
    for( x in seq(Lmodel.range[1], Lmodel.range[2], length=N) )
    {
        iy <- 1
        for( y in seq(Lnull.range[1], Lnull.range[2], length=N) )
        {
            pos <- ix*N  + iy
            d[pos,1] <- x
            d[pos,2] <- y
            #d[pos,3] <- McFadden(x,y)
            d[pos,3] <- CoxSnell(x,y,200)

            iy <- iy+1
        }
        ix <- ix+1
    }

    plotDf <- data.frame(d)
    colnames(plotDf) <- c("Lmodel", "Lnull", "McFadden")
    #print(plotDf)

    v.min <- min(trunc(plotDf$McFadden)-1, 0.0)
    v.max <- max(trunc(plotDf$McFadden)+1, 1.0)
    valScale <- scale_fill_gradientn(
        breaks=c(v.min, 0.0, 0.5, 1.0, v.max),
        limits=c(v.min, v.max),
        colors=        c("#000000", "#333333", "#7799ff",     "#aab088", "#ffff55"),
        values=rescale(c(    v.min,       0.0,       1.0,   1.0 + (v.max+1)/5,     v.max) ) )


    p <- ggplot( data = plotDf, aes(x = Lmodel, y = Lnull, fill = McFadden) ) +
         geom_raster() +
         labs( title = "McFadden pseudo-R^2" ) +
         coord_fixed() +
         valScale +
         theme( aspect.ratio = 1,
                plot.background  = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
         )
    print(p)
    #grid.newpage()
}

plotMcFaddenParams( N=80, Lmodel.range=c(0.001, 1.101), Lnull.range=c(0.001, 1.101) )
plotMcFaddenParams( N=80, Lmodel.range=c(0.001, 5.001), Lnull.range=c(0.001, 5.001) )
plotMcFaddenParams( N=80, Lmodel.range=c(0.001, 150.001), Lnull.range=c(0.001, 150.001) )
plotMcFaddenParams( N=80, Lmodel.range=c(0.001, 2.001), Lnull.range=c(0.8, 1.2) )
plotMcFaddenParams( N=80, Lmodel.range=c(0.001, 0.101), Lnull.range=c(0.9, 1.1) )



dev.off()