pdf("tree_regressions.out.pdf")

# Source: http://www2.math.su.se/PATHd8/
# Created using: ~/src/PATHd8/PATHd8 ./data/nmicrobiol201648-s6.txt ./data/nmicrobiol201648-s6.txt.nw.PATHd8.nw
tree <- read.tree("./data/nmicrobiol201648-s6.txt.nw.PATHd8.out.d8.nw")

bmcorr <- corBrownian( phy=tree )


# Load phenotypes

P <- 5

randomPhenotypes <- matrix(rnorm(N*P), N)
randomPhenData <- phylo4d(tree, randomPhenotypes)

# Create phenotypes that are correlated with phylogenetic distance, using BM
correlatedPhenotypes <- replicate(P, rTraitCont(tree))
correlatedPhenData <- phylo4d(tree, correlatedPhenotypes)

df <- data.frame( r1=randomPhenotypes[,1], r2=randomPhenotypes[,2], r3=randomPhenotypes[,3], r4=randomPhenotypes[,4], r5=randomPhenotypes[,5], c1=correlatedPhenotypes[,1], c2=correlatedPhenotypes[,2], c3=correlatedPhenotypes[,3], c4=correlatedPhenotypes[,4], c5=correlatedPhenotypes[,5])


#p <- ggplot(df, aes(r1, r2)) + geom_point(); p


# Try the regressions
#m1 <- lm(  r1 ~ r2, df); summary(m1)
#m1 <- gls( r1 ~ r2, df, correlation=bmcorr); summary(m1)


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
