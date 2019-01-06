library("ape")
library("phylobase")
library("nlme")
library("ggplot2")


seed <- as.integer(difftime( as.POSIXlt(Sys.time()), as.POSIXlt("2010-01-01 12:00:00"), units="secs" ))
set.seed( seed )
print(sprintf("Using random seed %d", seed))

pdf("test_tree_partial_correlation.out.pdf")

# Use birth-death model to generate tree
Ntarget = 200
N <- 0
while (N < Ntarget-3 || N >= Ntarget+3)
{
        #tree <- rbdtree(0.1, 0.025)
        tree <- rbdtree(0.1, 0.03)
        N <- tree$Nnode+1
}

N


getR2 <- function(model, yvals)
{
    meany <- mean(yvals)
    tss <- sum( (yvals - meany)**2 )
    res <- residuals(model)
    stopifnot( length(res) == length(yvals) )
    rss <- sum( res**2 )

    print("var")
    print(tss)
    print(tss/length(yvals))
    print(var(yvals))
    print("rss")
    print(rss)

    Rsquared <- 1-rss/tss
    print(Rsquared)
    #stopifnot( Rsquared >= 0.0 && Rsquared <= 1.0 )
    return( Rsquared )
}

performTest <- function(effectTrait, modelTrait, plotResult=TRUE, predictiveSD=0.1)
{
    # Create a GC%-like trait
    bmTrait <- rTraitCont(tree,
                          model="BM", sigma=0.07,
                          ancestor=FALSE, root.value=0.5)
    bmTrait[bmTrait > 1.0] <- 1.0
    bmTrait[bmTrait < 0.0] <- 0.0

    # Create one slightly correlated trait
    mixedTrait <- rTraitCont(tree,
                             model=function (x,l) rnorm(1, mean=x, sd=l*0.05) + rnorm(1, mean=0.0, sd=0.05),   # BM + uncorrelated Gaussian "noise"
                             ancestor=FALSE, root.value=0.5)

    # Create an uncorrelated trait
    uncorTrait <- rTraitCont(tree,
                             model=function (x,l) rnorm(1, mean=0.5, sd=0.2),   # random Gaussian values
                             ancestor=FALSE, root.value=0.5)


    traits <- rbind( bmTrait, mixedTrait, uncorTrait )

    randomFactors <- runif( 2, min=-1, max=1)
    predictiveTrait <- traits[effectTrait,] * randomFactors[1] + randomFactors[2] + rnorm( length(traits[effectTrait,]), mean=0, sd=predictiveSD )

    nonpredictiveTrait <- rTraitCont(tree,
                          model="BM", sigma=0.07,
                          ancestor=FALSE, root.value=0.5)
    nonpredictiveTrait[nonpredictiveTrait > 1.0] <- 1.0
    nonpredictiveTrait[nonpredictiveTrait < 0.0] <- 0.0

    traits <- rbind( traits, predictiveTrait, nonpredictiveTrait )
    traits <- data.frame(t(traits))
    #print(traits)

    predictiveFormula = as.formula( paste(effectTrait, " ~ ", modelTrait) )
    interceptFormula = as.formula( paste(effectTrait, " ~ ", "1") )

    m1 <- lm(  predictiveFormula, traits); summary(m1)
    colm <- coef(m1)
    print(colm)

    bmcorr <- corBrownian( phy=tree )

    print("===== gls1 =====")
    gls1 <- gls( predictiveFormula, traits, correlation=bmcorr, na.action=na.omit, method="ML"); print(summary(gls1))
    co1 <- coef(gls1)
    #print(co1)
    #r2 <- getR2(gls1, traits[,effectTrait] )
    #print(r2)
    #print(sqrt(r2))

    print("===== gls0 =====")
    gls0 <- gls( interceptFormula, traits, correlation=bmcorr, na.action=na.omit, method="ML"); print(summary(gls0))
    co0 <- coef(gls0)
    #print(co0)
    #r2 <- getR2(gls0, traits[,effectTrait] )
    #print(r2)
    #print(sqrt(r2))

    #c <- cor.test( residuals(gls1), residuals(gls0), method="pearson" )
    #print(summary(c))
    #print(c)

    pseudoR2 <- 1 - as.numeric( logLik(gls1) / logLik(gls0) )  # McFadden's pseudo-R^2
    #print(pseudoR2)


    if( plotResult )
    {
        clr2 <- "#55CCB0"

        p <- ggplot( data=traits, aes_string(x=modelTrait, y=effectTrait) ) + geom_point() +
                  geom_abline( aes( slope=colm[modelTrait],  intercept=colm["(Intercept)"]), colour="red"  ) +
                  geom_abline( aes( slope=co1 [modelTrait],  intercept=co1 ["(Intercept)"]), colour=clr2   ) +
                  geom_abline( aes( slope=0.0             ,  intercept=co0 ["(Intercept)"]), colour=clr2   )
        print(p)

}
    
    return( c( pseudoR2, summary(m1)$r.squared ) )
}

print("======================== predictive ========================")
performTest(effectTrait="mixedTrait", modelTrait="predictiveTrait")
print("======================== nonpredictive ========================")
performTest(effectTrait="mixedTrait", modelTrait="nonpredictiveTrait")


df <- data.frame( sd=numeric(), OLS.R2=numeric(), McFadden.R2=numeric() )

for( sd in seq(0.001, 0.5, length=2000) )
{
    rs <- performTest(effectTrait="mixedTrait", modelTrait="predictiveTrait", plotResult=FALSE, predictiveSD=sd)
    df <- rbind( df, data.frame( sd=c(sd), OLS.R2=c(rs[2]), McFadden.R2=c(rs[1]) ) )
}
print(df)
p <- ggplot( data=df, aes(x=OLS.R2, y=McFadden.R2) ) +
     geom_point(alpha=0.1)
print(p)


dev.off()
