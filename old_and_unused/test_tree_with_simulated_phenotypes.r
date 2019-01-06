library("ape")
library("phylobase")
library("ade4")
library("adephylo")
library("picante")
source("read.newick.R")   # Source: http://www.phytools.org/read.newick/v0.4/read.newick.R
                          # See: https://www.mail-archive.com/r-sig-phylo@r-project.org/msg02637.html

pdf("test_phylo4d.out.pdf")

#------------------------------------------
# Create tree
#------------------------------------------
#tree <- rtree(40)
# Create random tree using a birth-death process
# Since the same params create trees with a wide range of sizes, filter by target tree size (e.g., for testing statistical power)
Ntarget = 420
N <- 0
while (N < Ntarget-5 || N >= Ntarget+5)
{
        #tree <- rbdtree(0.1, 0.025)
        tree <- rbdtree(0.135, 0.03)
        N <- tree$Nnode+1
}
#tree <- read.newick(file="nmicro_s6_pruned.nw")
#tree <- collapse.singles(tree)
#N <- tree$Nnode+1

N
write.tree(tree, file="test_phylo4d.out.nw")
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

# Plot the trees along phenotypic variables
table.phylo4d(randomPhenData, box=FALSE)  
table.phylo4d(correlatedPhenData, box=FALSE)

#------------------------------------------------------------------
# Use autocorrelation to infer evolutionary effect on phenotypes
#------------------------------------------------------------------
# Calculate tree distances matrix
w <- 1/cophenetic(tree)
diag(w) <- 0

# Calculate Moran's autocorrelation index I (Gittleman and Kot) for each variable, using the tree distances)
# Ref: Paradis 2e @6.1.2, pp. 209-210.
Moran.I( randomPhenotypes[,1], w)
Moran.I( randomPhenotypes[,2], w)
Moran.I( randomPhenotypes[,3], w)
Moran.I( randomPhenotypes[,4], w)
Moran.I( randomPhenotypes[,5], w)
Moran.I( correlatedPhenotypes[,1], w)
Moran.I( correlatedPhenotypes[,2], w)
Moran.I( correlatedPhenotypes[,3], w)
Moran.I( correlatedPhenotypes[,4], w)
Moran.I( correlatedPhenotypes[,5], w)
N

# Use randomization procedure to test the significance of Moran's I
# Ref: Paradis 2e pp. 210-211
repeats = 1e5 # real
#repeats = 2e3 # draft
gearymoran(w, data.frame( randomPhenotypes ),     nrepet=repeats )
gearymoran(w, data.frame( correlatedPhenotypes ), nrepet=repeats )

# Test phylogenetic signal statistic K, using randomization
# Ref: Paradis 2e pp. 236-237
# Note: This test seems to have more power than the other tests used here.
repeats = 1e5  # real
#repeats = 5e3  # draft
x <- list( v1=randomPhenotypes[,1], v2=randomPhenotypes[,2], v3=randomPhenotypes[,3], v4=randomPhenotypes[,4], v5=randomPhenotypes[,5] )
sapply( x, phylosignal, tree, reps=repeats )
x <- list( v1=correlatedPhenotypes[,1], v2=correlatedPhenotypes[,2], v3=correlatedPhenotypes[,3], v4=correlatedPhenotypes[,4], v5=correlatedPhenotypes[,5] )
sapply( x, phylosignal, tree, reps=repeats )

dev.off()

