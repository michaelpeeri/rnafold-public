library("ape")
library("phylobase")
library("nlme")
library("ggplot2")
library("Matrix")


seed <- as.integer(difftime( as.POSIXlt(Sys.time()), as.POSIXlt("2010-01-01 12:00:00"), units="secs" ))
set.seed( seed )
print(sprintf("Using random seed %d", seed))

pdf("test_buse_pseudo_r_squared.out.pdf")

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

performTest <- function(effectTrait, modelTrait, plotResult=TRUE, predictiveSD=0.1 )
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


    #    A   B   C   D
    # A 0.0 0.9 0.1 0.0
    # B 0.3 0.0 0.7 0.0
    # C 0.2 0.4 0.0 0.4
    # D 0.0 0.1 0.9 0.0
    #
    discreteUncorrelatedTrait <- rTraitDisc(tree,
                                            model=matrix( c(c(0.0,0.9,0.1,0.0),c(0.3,0.0,0.7,0.0),c(0.2,0.4,0.0,0.4),c(0.0,0.1,0.9,0.0)), ncol=4 ),
                                            states=factor(c("TraitVal.A", "TraitVal.B", "TraitVal.C", "TraitVal.D")),
                                            k=4,
                                            root.value=sample(4) )


    traits <- rbind( bmTrait, mixedTrait, uncorTrait, discreteUncorrelatedTrait )

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
    print("===== lm =====")
    print(summary(m1))

    if( TRUE )
    {
        bmcorr <- corBrownian( phy=tree )  # perform phylogentic GLS
    }
    else
    {
        bmcorr <- corSymm( fixed=TRUE )   # No correlations - use GLS to calculate OLS regression (Test mode - regression and R^2 should coincide with OLS model)
    }
        
    print("===== gls1 =====")
    gls1 <- gls( predictiveFormula, traits, correlation=bmcorr, na.action=na.omit, method="REML"); print(summary(gls1))
    co1 <- coef(gls1)

    print("===== gls0 =====")
    gls0 <- gls( interceptFormula, traits, correlation=bmcorr, na.action=na.omit, method="REML"); print(summary(gls0))
    co0 <- coef(gls0)

    #pseudoR2 <- 1 - as.numeric( logLik(gls1) / logLik(gls0) )  # McFadden's pseudo-R^2
    #print(pseudoR2)

    #
    # Buse's R^2
    # Reference: A. Buse, "Goodness of Fit in Generalized Least Squares Estimation", The American Statistician, June 1973, Vol. 27, No. 3.
    # URL: http://www.jstor.org/stable/2683631
    #
    # Notation: u.hat := epsilon (residuals)
    #           V := Omega (Covariance matrix)
    #           Y.bar := Intercept of equivalent intercept-only model (see: http://r.789695.n4.nabble.com/Pseudo-R-squared-in-gls-model-td4641148.html)
    #           e := first column in design matrix X (i.e., rep(1, N) if intercepts are used). See above eq. (11) p. 107.
    #
    u.hat <- residuals(gls1)
    #V <- getVarCov( gls1 )   # this doesn't work (bug?)
    V <- corMatrix( gls1$modelStruct$corStruct )
    inv.V <- solve(V)   # invert V (=Omega). Every positive-definite matrix is invertible (https://en.wikipedia.org/wiki/Positive-definite_matrix#cite_note-5)
    e <- rep(1, 2)
    Y <- traits[,effectTrait]
    Y.bar <- coef(gls0)["(Intercept)"]
    yye <- Y - Y.bar * e
    
    pseudoR2 <- 1 - ( t(u.hat) %*% inv.V %*% u.hat )/( t(yye) %*% inv.V %*% yye )   # See eq. (15) p. 107.
    stopifnot( pseudoR2 >= 0 && pseudoR2 <= 1.0 )


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

print("======================== nonpredictive ========================")
performTest(effectTrait="mixedTrait", modelTrait="discreteUncorrelatedTrait")

dev.off() 
quit()

df <- data.frame( sd=numeric(), OLS.R2=numeric(), Buse.R2=numeric() )

for( sd in seq(0.001, 0.5, length=200) )
{    
    rs <- performTest(effectTrait="mixedTrait", modelTrait="predictiveTrait", plotResult=FALSE, predictiveSD=sd)
    df <- rbind( df, data.frame( sd=c(sd), OLS.R2=c(rs[2]), Buse.R2=c(rs[1]) ) )
}
for( sd in seq(0.3, 1.5, length=100) )
{    
    rs <- performTest(effectTrait="mixedTrait", modelTrait="predictiveTrait", plotResult=FALSE, predictiveSD=sd)
    df <- rbind( df, data.frame( sd=c(sd), OLS.R2=c(rs[2]), Buse.R2=c(rs[1]) ) )
}
print(df)
p <- ggplot( data=df, aes(x=OLS.R2, y=Buse.R2) ) +
     geom_point(alpha=0.1)
print(p)


dev.off()
