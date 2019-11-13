library("phylobase")
library("ape")
library("ggplot2")
library("reshape")
library("nlme")


pdf("test_lm_gls_pvalues.out.pdf")

#------------------------------------------
# Create tree
#------------------------------------------
#tree <- rtree(40)
# Create random tree using a birth-death process
# Since the same params create trees with a wide range of sizes, filter by target tree size (e.g., for testing statistical power)
Ntarget = 200
Ntree <- 0
while (Ntree < Ntarget-5 || Ntree >= Ntarget+5)
{
        #tree <- rbdtree(0.1, 0.025)
        tree <- rbdtree(0.1, 0.03)
        Ntree <- tree$Nnode+1
}
Ntree
#write.tree(tree, file="test_lm_gls_pvalues.out.nw")

bmcorr <- corBrownian( phy=tree )


test_regression_methods <- function()
{
    correlatedPhenotypes <- replicate(2, rTraitCont(tree))
    #correlatedPhenData <- phylo4d(tree, correlatedPhenotypes)

    df <- data.frame( c1=correlatedPhenotypes[,1], c2=correlatedPhenotypes[,2] )


    lm1  <- lm(  c1 ~ c2, df); #summary(lm1)
    gls1 <- gls( c1 ~ c2, df, correlation=bmcorr); #summary(gls1)

    return( c( summary(lm1)$coefficients[2,4], summary(gls1)$tTable[2,4]) )
}

repeats <- 10000

pvals <- data.frame( lm.pval = rep(0.0, repeats), gls.pval = rep(0.0, repeats) )

for (i in 1:repeats)
{
    a = test_regression_methods()
    pvals[i, "lm.pval"]  <- a[1]
    pvals[i, "gls.pval"] <- a[2]
}

#print(pvals)
#

plotdf <- data.frame( melt( pvals, measure.vars=c(1,2) ) )
#print(dimnames(plotdf)[[2]])
#print(plotdf)

theme_set(theme_gray(base_size=18))

p <- ggplot( data=plotdf, aes(value, colour=variable ) ) +
    labs( x="P-value" ) +
    geom_freqpoly( bins=200, size=1.8 ) +
    annotate( "text", x=Inf,  y=Inf, label=sprintf("N = %d\n|tree tips| = %d", repeats, Ntree ), hjust=1, vjust=1, size=6 ); p


dev.off()