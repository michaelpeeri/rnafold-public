library("ape")
library("phylobase")

pdf('test_bm_tree_simulation.out.pdf')

# Use birth-death model to generate tree
Ntarget = 100
N <- 0
while (N < Ntarget-3 || N >= Ntarget+3)
{
        #tree <- rbdtree(0.1, 0.025)
        tree <- rbdtree(0.1, 0.03)
        N <- tree$Nnode+1
}

N

x <- rTraitCont(tree, model="BM", sigma=0.07, alpha=1, theta=0, ancestor=TRUE, root.value=0.5)
x[x>1.0] <- 1.0
x[x<0.0] <- 0.0

print(nTips(tree))
print(class(x))
print(length(x))
print("-------------")
print(range(x))
print(mean(x))

traitToString <- function(v)
{
    stopifnot(v >= 0.0 && v <= 1.0)
    intval <- round(v*255)
    return(sprintf("#%02X%02X%02X", intval, intval, intval))
}

y <- vapply(x, traitToString, character(1))
#print(y)
print(class(y))
print(length(y))

plot(tree, type="phylo", show.tip.label=FALSE, y.lim=100.0 )
nodelabels(pch=21, bg=y)
tiplabels(pch=21, bg=y)
add.scale.bar()

dev.off()

warnings()
