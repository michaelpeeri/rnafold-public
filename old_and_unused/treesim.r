library("ape")
#library("apTreeshape")

pdf('rlineage.pdf')

# Use birth-death model to generate tree
t <- rbdtree(0.1, 0.025)

# Use splitting probabilities proprtional to the size of the clade
# Source: Emmanuel Paradis / Analysis of Phylogenetics and Evolution with R (2nd Edition), Section 7.1, p. 319.
#Q <- function(n,i) if (i > 0 && i < n) n else 0
#t <- rtreeshape(n=1, tip.number=120, FUN = Q)
#plot.treeshape(t[[1]], show.tip.label=FALSE)

x <- rTraitCont(t, model="BM", sigma=0.01, alpha=1, theta=0, ancestor=FALSE, root.value=0.5)
x

#plot(t, show.tip.label=FALSE); axisPhylo()
plot(t, type="u", show.tip.label=FALSE, no.margin=TRUE)
add.scale.bar()

dev.off()
