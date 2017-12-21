# ==========================================================================================
# Test program for my implementation of Buse[1] R^2 for GLS
# * Simulate trees and random traits and test how GLS and OLS regressions work on them.
# * Plot OLS R^2 values against GLS R^2 values
#
# Reference: [1] A. Buse, "Goodness of Fit in Generalized Least Squares Estimation", The American Statistician, June 1973, Vol. 27, No. 3.
# URL: http://www.jstor.org/stable/2683631

library("ape")
library("phylobase")
library("nlme")
library("ggplot2")
library("Matrix")

# ==========================================================================================
# Configuration

# FALSE - perform gls-based analysis; TRUE - use gls code to simulate OLS (testing mode - should give the same p-value and R^2 as OLS)
testingModeSimulateOLS <- TRUE

# How many points to plot in the validation plot
validationTestPoints = 10

# Tree target (approximate) size
Ntarget = 200

# PDF output
pdf("test_buse_pseudo_r_squared.out.pdf")
# ==========================================================================================

# ==========================================================================================
# Set random seed
seed <- as.integer(difftime( as.POSIXlt(Sys.time()), as.POSIXlt("2010-01-01 12:00:00"), units="secs" ))
set.seed( seed )
print(sprintf("Using random seed %d", seed))


# ==========================================================================================
# Use birth-death model to generate tree
# (since the model returns trees of widely varying size and I would like to continually observe models having roughly the same power, I repeat the process until the tree size falls within an acceptable range)
# Note: this same tree will be used throughout the run of the program
N <- 0
while (N < Ntarget-3 || N >= Ntarget+3)
{
        tree <- rbdtree(0.1, 0.03)
        N <- tree$Nnode+1
}

N


# ==========================================================================================
# Do GLS a single regression between randomly-generated traits
# Note: each test creates a new set of random traits, but the tree is fixed throughout the run of the program
#
# effectTrait  - explained variable (dependent variable)
# modelTrait   - explanatory variable (=independent variable)
# plotResult   - draw regression plot of the data
# predictiveSD - SD of the 'noise' distribution used for the random (used to simulate random traits with varying degrees of predictive quality)
# ==========================================================================================
performTest <- function(effectTrait, modelTrait, plotResult=TRUE, predictiveSD=0.1 )
{
    # ==========================================================================================
    # ============================ Part 1 - Generate Random Traits  ============================
    # ==========================================================================================
    
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

    traits <- NA
    
    if( !any(modelTrait==c('predictiveTrait', 'nonpredictiveTrait') ) )
    {
        # ============================ Discrete regressor (explanatory) trait  ============================
        #
        # Create a discrete trait based on the tree (using a Markovian process)
        #
        # Test transition matrix
        #    A   B   C   D
        # A 0.0 0.9 0.1 0.0
        # B 0.3 0.0 0.7 0.0
        # C 0.2 0.4 0.0 0.4
        # D 0.0 0.1 0.9 0.0
        #
        discreteUncorrelatedTrait <- as.factor( rTraitDisc( tree,
                                                            model=matrix( c(c(0.0,0.9,0.1,0.0),c(0.3,0.0,0.7,0.0),c(0.2,0.4,0.0,0.4),c(0.0,0.1,0.9,0.0)), ncol=4 ),
                                                            states=factor(c("TraitVal.A", "TraitVal.B", "TraitVal.C", "TraitVal.D")),
                                                            k=4,
                                                            root.value=sample(4) ) )

        selectedTraitVec <- as.numeric(discreteUncorrelatedTrait)
        stopifnot(min(selectedTraitVec)==1)
        stopifnot(max(selectedTraitVec)==4)

        # Create a continuous trait that is correlated with the discrete trait
        centers <- runif( 4, min=-1, max=1 )  # Create center values for all trait levels
        correlatedTraitForDiscreteTrait <- rnorm( N, mean=0.0, sd=predictiveSD ) + centers[selectedTraitVec]  # Draw continuous trait values - for each node based on the center for its discrete trait level

        # Store all generated traits in a data-frame
        traits <- data.frame( bmTrait=bmTrait, mixedTrait=mixedTrait, uncorTrait=uncorTrait, discreteUncorrelatedTrait=discreteUncorrelatedTrait, correlatedTraitForDiscreteTrait=correlatedTraitForDiscreteTrait )
        
    }
    else
    {
        # ============================ Continuous regressor (explanatory) trait  ============================
        # Create a predictive trait, based on linear function of the effect trait (=x), y = a*x + b + epsilon
        # The noise component (epsilon) has mean 0 and sd=predictiveSD
        randomFactors <- runif( 2, min=-1, max=1)   # draw 'a' and 'b' from a uniform distribution
        explanatoryTraits <- data.frame( bmTrait=bmTrait, mixedTrait=mixedTrait, uncorTrait=uncorTrait ) # store all explanatory traits in a temporary data-frame so they would be easily accessible
        predictiveTrait <- explanatoryTraits[,effectTrait] * randomFactors[1] + randomFactors[2] + rnorm( length(explanatoryTraits[,effectTrait]), mean=0, sd=predictiveSD )
        #print(predictiveTrait)

        nonpredictiveTrait <- rTraitCont(tree,
                              model="BM", sigma=0.07,
                              ancestor=FALSE, root.value=0.5)
        nonpredictiveTrait[nonpredictiveTrait > 1.0] <- 1.0
        nonpredictiveTrait[nonpredictiveTrait < 0.0] <- 0.0

        # Store all generated traits in a data-frame
        traits <- data.frame( bmTrait=bmTrait, mixedTrait=mixedTrait, uncorTrait=uncorTrait, predictiveTrait=predictiveTrait, nonpredictiveTrait=nonpredictiveTrait )
    }
        
    stopifnot(any(effectTrait==colnames(traits)))   # requested effectTrait not found
    stopifnot(any(modelTrait==colnames(traits)))    # requested modelTrait not found

    # ==========================================================================================
    # ============================ Part 2 - Perform Regressions ================================
    # ==========================================================================================
    
    # Set the predictive (regression) and intecept-only formulas
    # 'predictiveFormula' is the actual regression. It uses an intercept when a continuous explanatory var is used (and none with a discrete var, as in one-way ANOVA)
    # 'inteceptFormula' is an intecept-only model, used as the baseline for calculation of R^2
    predictiveFormula <- NA    
    if( any(class(traits[,modelTrait])==c("factor")) )
    {
        predictiveFormula <- as.formula( paste(effectTrait, " ~ ", modelTrait, " + 0") )
    }
    else
    {
        predictiveFormula <- as.formula( paste(effectTrait, " ~ ", modelTrait) )
    }
    interceptFormula = as.formula( paste(effectTrait, " ~ ", "1") )
    print(predictiveFormula)

    # Perform OLS regression (used as a reference for comparison)
    m1 <- lm(  predictiveFormula, traits); summary(m1)
    colm <- coef(m1)
    print("===== lm =====")
    print(summary(m1))

    # Perform GLS and calculate R^2
    
    # Determine the variance-covariance matrix to be used for GLS
    if( !testingModeSimulateOLS )
    {
        bmcorr <- corBrownian( phy=tree )  # perform phylogentic GLS
    }
    else
    {
        bmcorr <- corSymm( fixed=TRUE )   # No correlations - use GLS to calculate OLS regression (Test mode - regression and R^2 should coincide with OLS model)
    }

    # Perform GLS regression using the predictive and intercept-only models
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
    #                (intercepts are used for continuous explanatory variables (independent variables), but not for discrete ones).
    #
    u.hat <- residuals(gls1)
    #V <- getVarCov( gls1 )   # this doesn't work (bug?)
    V <- corMatrix( gls1$modelStruct$corStruct )
    inv.V <- solve(V)   # invert V (=Omega). Every positive-definite matrix is invertible (https://en.wikipedia.org/wiki/Positive-definite_matrix#cite_note-5)
    if( any(class(traits[,modelTrait])==c("factor")) )
    {
        e <- rep(0, N)
    }
    else
    {
        e <- rep(1, N)
    }
    Y <- traits[,effectTrait]
    Y.bar <- coef(gls0)["(Intercept)"]
    yye <- Y - Y.bar * e
    
    R2 <- 1 - ( t(u.hat) %*% inv.V %*% u.hat )/( t(yye) %*% inv.V %*% yye )   # See eq. (15) p. 107.
    stopifnot( R2 >= 0 && R2 <= 1.0 )

    print("////------")
    print(anova(m1))
    print("////------")
    print(anova(gls1))

    pvalue <- NA
    pvalue.OLS <- NA
    if( any(class(traits[,modelTrait])==c("factor")) )
    {
        pvalue     <- anova(gls1)[1,'p-value']
        pvalue.OLS <- anova(m1)[1,'Pr(>F)']
    }
    else
    {
        pvalue     <- summary(gls1)$tTable[2,4]
        pvalue.OLS <- summary(m1)$coefficients[2,4]
    }

    # ==========================================================================================
    # ========================= Part 3 (optional) - Plot Regression ============================
    # ==========================================================================================
    if( plotResult )
    {
        clr2 <- "#55CCB0"
        max.y <- max(traits[,effectTrait])
        min.y <- min(traits[,effectTrait])
        line.y <- (max.y-min.y)/32

        if( any(class(traits[,modelTrait])==c("factor")) )
        {
            p <- ggplot( data=traits, aes_string(x=modelTrait, y=effectTrait) ) +
                geom_boxplot(outlier.size=0, outlier.alpha=0) +  # don't show outliers on boxplot, since we're plotting all data anyway...
                geom_jitter(colour="grey") +
                annotate( "text", x=Inf, y=max.y - c(0,1,2,3) * line.y, label=c(  "lm",   sprintf("italic(R) ^ 2 == %.3g", summary(m1)$r.squared), sprintf('italic(p)*"-val" ==  %.3g', pvalue.OLS ), sprintf("italic(N) == %d", nrow(traits))), hjust=1, vjust=0, colour="red", parse=TRUE ) +
                annotate( "text", x=-Inf, y=max.y - c(0,1,2,3) * line.y, label=c(  "gls", sprintf("italic(R) ^ 2 == %.3g", R2)                   , sprintf('italic(p)*"-val" ==  %.3g', pvalue     ), sprintf("italic(N) == %d", nrow(traits))), hjust=0, vjust=0, colour=clr2, parse=TRUE )
            print(p)
        }
        else
        {
            p <- ggplot( data=traits, aes_string(x=modelTrait, y=effectTrait) ) + geom_point() +
                      geom_abline( aes( slope=colm[modelTrait],  intercept=colm["(Intercept)"]), colour="red"  ) +
                      geom_abline( aes( slope=co1 [modelTrait],  intercept=co1 ["(Intercept)"]), colour=clr2   ) +
                      geom_abline( aes( slope=0.0             ,  intercept=co0 ["(Intercept)"]), colour=clr2   ) +
                      annotate( "text", x=Inf,  y=max.y+c(0,1,2,3) * line.y, label=c(  "lm", sprintf("italic(R) ^ 2 == %.3g", summary(m1)$r.squared), sprintf('italic(p)*"-val" ==  %.3g', pvalue.OLS ), sprintf("italic(N) == %d", nrow(traits))), hjust=1, vjust=0, colour="red", parse=TRUE ) +
                      annotate( "text", x=-Inf, y=max.y+c(0,1,2,3) * line.y, label=c(  "gls", sprintf("italic(R) ^ 2 == %.3g", R2),                   sprintf('italic(p)*"-val" ==  %.3g', pvalue     ),     sprintf("italic(N) == %d", nrow(traits))), hjust=0, vjust=0, colour=clr2, parse=TRUE )
            print(p)
        }
    }
    
    return( c( R2, summary(m1)$r.squared, pvalue, pvalue.OLS ) )
}


# ==========================================================================================
# Perform regression with multiple sets of simulated traits; plot GLS R^2 values
# against OLS R^2 values
validationRun <- function( useDiscreteExplanatoryVar=FALSE )
{
    # Validation run - create many random traits of various predictive quality, and plot the calculated R^2 for GLS and OLS models

    df <- data.frame( sd=numeric(), OLS.R2=numeric(), Buse.R2=numeric(), OLS.pval=numeric(), GLS.pval=numeric(), predictive=logical() )
    rs.predictive    <- NA
    rs.nonpredictive <- NA

    # Simulate using low values of predictiveSD
    for( sd in seq(0.001, 0.5, length=round(validationTestPoints*0.5)) )
    {
        if( !useDiscreteExplanatoryVar )
        {
            # Test for continuous traits
            rs.predictive    <- performTest(effectTrait="mixedTrait", modelTrait="predictiveTrait"   , plotResult=FALSE, predictiveSD=sd)
            rs.nonpredictive <- performTest(effectTrait="mixedTrait", modelTrait="nonpredictiveTrait", plotResult=FALSE, predictiveSD=sd)
        }
        else
        {
            # Test for discrete traits
            rs.predictive    <- performTest(effectTrait="correlatedTraitForDiscreteTrait", modelTrait="discreteUncorrelatedTrait", plotResult=FALSE, predictiveSD=sd)
            rs.nonpredictive <- performTest(effectTrait="correlatedTraitForDiscreteTrait", modelTrait="mixedTrait",                plotResult=FALSE, predictiveSD=sd)
        }
        df <- rbind( df, data.frame( sd=c(sd), OLS.R2=c(rs.predictive[2]),    Buse.R2=c(rs.predictive[1]),    OLS.pval=rs.predictive[4],    GLS.pval=rs.predictive[3], predictive=TRUE  ) )
        df <- rbind( df, data.frame( sd=c(sd), OLS.R2=c(rs.nonpredictive[2]), Buse.R2=c(rs.nonpredictive[1]), OLS.pval=rs.nonpredictive[4], GLS.pval=rs.nonpredictive[3], predictive=FALSE ) )
    }

    # Simulate using higher values of predictiveSD
    for( sd in seq(0.3, 1.5, length=round(validationTestPoints*0.5)) )
    {
        if( !useDiscreteExplanatoryVar )
        {
            # Test for continuous traits
            rs.predictive    <- performTest(effectTrait="mixedTrait", modelTrait="predictiveTrait"   , plotResult=FALSE, predictiveSD=sd)
            rs.nonpredictive <- performTest(effectTrait="mixedTrait", modelTrait="nonpredictiveTrait", plotResult=FALSE, predictiveSD=sd)
        }
        else
        {
            # Test for discrete traits
            rs.predictive    <- performTest(effectTrait="correlatedTraitForDiscreteTrait", modelTrait="discreteUncorrelatedTrait", plotResult=FALSE, predictiveSD=sd)
            rs.nonpredictive <- performTest(effectTrait="correlatedTraitForDiscreteTrait", modelTrait="mixedTrait",                plotResult=FALSE, predictiveSD=sd)
        }
        df <- rbind( df, data.frame( sd=c(sd), OLS.R2=c(rs.predictive[2]),    Buse.R2=c(rs.predictive[1]),    OLS.pval=rs.predictive[4],    GLS.pval=rs.predictive[3], predictive=TRUE  ) )
        df <- rbind( df, data.frame( sd=c(sd), OLS.R2=c(rs.nonpredictive[2]), Buse.R2=c(rs.nonpredictive[1]), OLS.pval=rs.nonpredictive[4], GLS.pval=rs.nonpredictive[3], predictive=FALSE ) )
    }

    # Plot R^2 values - OLS vs. GLS
    # If testingModeSimulateOLS==TRUE, all values should be equal (x==y)
    alpha.adaptive <- 1/(1+exp((validationTestPoints-400)/100))*0.9+0.1  # use alpha when there are many points
    p <- ggplot( data=df, aes(x=OLS.R2, y=Buse.R2, colour=predictive) ) +
         geom_point(alpha=alpha.adaptive)
    print(p)

    # Plot p-values for predictive and non-predictive 
    #clr1 <- "red"
    #clr2 <- "#55CCB0"
    nbins.adaptive = min( 100, max(20, round(validationTestPoints/4) ) )
    p <- ggplot( data=df, aes( colour=predictive )) +
        geom_freqpoly( aes(OLS.pval), bins=nbins.adaptive ) +
        geom_freqpoly( aes(GLS.pval), bins=nbins.adaptive )
    print(p)
    
}



# ==========================================================================================
# Tests runs with continuous explanatory traits
print("======================== predictive ========================")
performTest(effectTrait="mixedTrait", modelTrait="predictiveTrait")
print("======================== nonpredictive ========================")
performTest(effectTrait="mixedTrait", modelTrait="nonpredictiveTrait")

validationRun( useDiscreteExplanatoryVar=FALSE )

# ==========================================================================================
# Tests runs with discerete explanatory traits
print("======================== predictive ========================")
performTest(effectTrait="correlatedTraitForDiscreteTrait", modelTrait="discreteUncorrelatedTrait", predictiveSD=1.0)
print("======================== nonpredictive ========================")
performTest(effectTrait="mixedTrait",                      modelTrait="discreteUncorrelatedTrait", predictiveSD=1.0)


validationRun( useDiscreteExplanatoryVar=TRUE )


dev.off()
