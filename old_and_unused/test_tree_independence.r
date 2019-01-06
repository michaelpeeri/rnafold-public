library("phylobase")
library("ape")
library("ggplot2")
library("nlme")

pdf("test_tree_independence.out.pdf")

#------------------------------------------
# Create tree
#------------------------------------------
#tree <- rtree(40)
# Create random tree using a birth-death process
# Since the same params create trees with a wide range of sizes, filter by target tree size (e.g., for testing statistical power)
Ntarget = 330
N <- 0
while (N < Ntarget-5 || N >= Ntarget+5)
{
        #tree <- rbdtree(0.1, 0.025)
        tree <- rbdtree(0.135, 0.03)
        N <- tree$Nnode+1
}

N
#write.tree(tree, file="test_phylo4d.out.nw")

#------------------------------------------
# Create random phenotypes
#------------------------------------------
P <- 5   # Number of random phenotypes of each type - must be 5 (or change below)...

# Create phenotypes that are uncorrelated with phylogenetic distance
randomPhenotypes <- matrix(rnorm(N*P), N)
randomPhenData <- phylo4d(tree, randomPhenotypes)

# Create phenotypes that are correlated with phylogenetic distance, using BM
correlatedPhenotypes <- replicate(P, rTraitCont(tree))
correlatedPhenData <- phylo4d(tree, correlatedPhenotypes)

df <- data.frame( r1=randomPhenotypes[,1], r2=randomPhenotypes[,2], r3=randomPhenotypes[,3], r4=randomPhenotypes[,4], r5=randomPhenotypes[,5], c1=correlatedPhenotypes[,1], c2=correlatedPhenotypes[,2], c3=correlatedPhenotypes[,3], c4=correlatedPhenotypes[,4], c5=correlatedPhenotypes[,5])


#-------------------------------------------------
# Plot binary relations between random variables
#-------------------------------------------------
p <- ggplot(df, aes(r1, r2)) + geom_point(); p
p <- ggplot(df, aes(r1, r3)) + geom_point(); p
p <- ggplot(df, aes(r1, r4)) + geom_point(); p
p <- ggplot(df, aes(r1, r5)) + geom_point(); p
p <- ggplot(df, aes(r2, r3)) + geom_point(); p
p <- ggplot(df, aes(r2, r4)) + geom_point(); p
p <- ggplot(df, aes(r2, r5)) + geom_point(); p
p <- ggplot(df, aes(r3, r4)) + geom_point(); p
p <- ggplot(df, aes(r3, r5)) + geom_point(); p
p <- ggplot(df, aes(r4, r5)) + geom_point(); p

p <- ggplot(df, aes(c1, c2)) + geom_point(); p
p <- ggplot(df, aes(c1, c3)) + geom_point(); p
p <- ggplot(df, aes(c1, c4)) + geom_point(); p
p <- ggplot(df, aes(c1, c5)) + geom_point(); p
p <- ggplot(df, aes(c2, c3)) + geom_point(); p
p <- ggplot(df, aes(c2, c4)) + geom_point(); p
p <- ggplot(df, aes(c2, c5)) + geom_point(); p
p <- ggplot(df, aes(c3, c4)) + geom_point(); p
p <- ggplot(df, aes(c3, c5)) + geom_point(); p
p <- ggplot(df, aes(c4, c5)) + geom_point(); p

p <- ggplot(df, aes(r1, c1)) + geom_point(); p
p <- ggplot(df, aes(r2, c2)) + geom_point(); p
p <- ggplot(df, aes(r3, c3)) + geom_point(); p
p <- ggplot(df, aes(r4, c4)) + geom_point(); p
p <- ggplot(df, aes(r5, c5)) + geom_point(); p

#---------------------------------------------------------------------------------
# Test linear regression models :
# 1. Ordinary least-squares models, which assume the residuals are uncorrelated.
#    This assumption holds for r1..r5, but c1..c5 where created using the same
#    tree, so may appear to be correlated.
# 2. Generalized least-square models, which model correlation between samples
#    using a variance-covariance matrix (in this case derived from the tree by
#    assuming a BM model). These models should not find any of the pairs to
#    be significantly correlated (except by chance).
#---------------------------------------------------------------------------------
bmcorr <- corBrownian( phy=tree )

print("------------", quote=FALSE); m1 <- lm(  r1 ~ r2, df); summary(m1)
m1 <- gls( r1 ~ r2, df, correlation=bmcorr); summary(m1)

print("------------", quote=FALSE); m1 <- lm(  r1 ~ r3, df); summary(m1)
m1 <- gls( r1 ~ r3, df, correlation=bmcorr); summary(m1)

print("------------", quote=FALSE); m1 <- lm(  r1 ~ r4, df); summary(m1)
m1 <- gls( r1 ~ r4, df, correlation=bmcorr); summary(m1)

print("------------", quote=FALSE); m1 <- lm(  r1 ~ r5, df); summary(m1)
m1 <- gls( r1 ~ r5, df, correlation=bmcorr); summary(m1)

print("------------", quote=FALSE); m1 <- lm(  r2 ~ r3, df); summary(m1)
m1 <- gls( r2 ~ r3, df, correlation=bmcorr); summary(m1)

print("------------", quote=FALSE); m1 <- lm(  r2 ~ r4, df); summary(m1)
m1 <- gls( r2 ~ r4, df, correlation=bmcorr); summary(m1)

print("------------", quote=FALSE); m1 <- lm(  r2 ~ r5, df); summary(m1)
m1 <- gls( r2 ~ r5, df, correlation=bmcorr); summary(m1)

print("------------", quote=FALSE); m1 <- lm(  r3 ~ r4, df); summary(m1)
m1 <- gls( r3 ~ r4, df, correlation=bmcorr); summary(m1)

print("------------", quote=FALSE); m1 <- lm(  r3 ~ r5, df); summary(m1)
m1 <- gls( r3 ~ r5, df, correlation=bmcorr); summary(m1)

print("------------", quote=FALSE); m1 <- lm(  r4 ~ r5, df); summary(m1)
m1 <- gls( r4 ~ r5, df, correlation=bmcorr); summary(m1)


print("------------", quote=FALSE); m1 <- lm(  c1 ~ c2, df); summary(m1)
m1 <- gls( c1 ~ c2, df, correlation=bmcorr); summary(m1)

print("------------", quote=FALSE); m1 <- lm(  c1 ~ c3, df); summary(m1)
m1 <- gls( c1 ~ c3, df, correlation=bmcorr); summary(m1)

print("------------", quote=FALSE); m1 <- lm(  c1 ~ c4, df); summary(m1)
m1 <- gls( c1 ~ c4, df, correlation=bmcorr); summary(m1)

print("------------", quote=FALSE); m1 <- lm(  c1 ~ c5, df); summary(m1)
m1 <- gls( c1 ~ c5, df, correlation=bmcorr); summary(m1)

print("------------", quote=FALSE); m1 <- lm(  c2 ~ c3, df); summary(m1)
m1 <- gls( c2 ~ c3, df, correlation=bmcorr); summary(m1)

print("------------", quote=FALSE); m1 <- lm(  c2 ~ c4, df); summary(m1)
m1 <- gls( c2 ~ c4, df, correlation=bmcorr); summary(m1)

print("------------", quote=FALSE); m1 <- lm(  c2 ~ c5, df); summary(m1)
m1 <- gls( c2 ~ c5, df, correlation=bmcorr); summary(m1)

print("------------", quote=FALSE); m1 <- lm(  c3 ~ c4, df); summary(m1)
m1 <- gls( c3 ~ c4, df, correlation=bmcorr); summary(m1)

print("------------", quote=FALSE); m1 <- lm(  c3 ~ c5, df); summary(m1)
m1 <- gls( c3 ~ c5, df, correlation=bmcorr); summary(m1)

print("------------", quote=FALSE); m1 <- lm(  c4 ~ c5, df); summary(m1)
m1 <- gls( c4 ~ c5, df, correlation=bmcorr); summary(m1)




dev.off()

