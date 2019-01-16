library("ggplot2")
library("reshape")
library("rhdf5")
library("rredis")
library("phylobase")
library("ape")
library("ade4")
library("adephylo")
library("picante")
library("phylosignal")

pdf("tree_trait_correlograms.out.pdf")



#---------------------------------------------------------------
# Configuration

redisConnect(host="power5", password="rnafold")

profileStart <- 0
profileStop <- 1000
profileStep <- 10
profileReference <- "begin"
profileLen <- (profileStop-profileStart)/profileStep
#---------------------------------------------------------------



getH5Filename <- function( taxId )
{
    sprintf("gcdata_v2_taxid_%d_profile_%d_%d_%s_%d.h5", taxId, profileStop, profileStep, profileReference, profileStart)
}

readDeltaLFEProfile <- function( taxId, h5filename )
{
    y <- h5ls(h5filename)
    key <- y[which(startsWith(y[,1], "/df_"))[1], 1]   # find the first element that starts with '/df_'

    x <- h5read(h5filename, key )

    stopifnot(x$block0_items[[2]]=="native")
    native <- x$block0_values[2,]

    stopifnot(x$block0_items[[3]]=="shuffled")
    shuffled <- x$block0_values[3,]

    return( native - shuffled )
}

getGenomicGCContent <- function(taxId)
{
    val <- redisGet(sprintf("species:taxid:%s:properties:gc-content", taxId))
    
    if( !is.null(val)) {
        val <- as.double(val)
        stopifnot( val > 15.0 && val < 85.0 )
        return( val )
    }
    else
    {
        return( NA )  # convert NULL to NA
    }
}

getPairedFraction <- function(taxId)
{
    val <- redisGet(sprintf("species:taxid:%s:properties:paired-mRNA-fraction", taxId))
    
    if( !is.null(val)) {
        val <- as.double(val)
        stopifnot( val > 0.1 && val < 0.9 )
        return( val )
    }
    else
    {
        return( NA )  # convert NULL to NA
    }
}

readAllProfiles <- function( taxIds )
{

    combined <- data.frame( GenomicGC=double(), PairedFraction=double(), Profile=matrix( rep(0.0, profileLen), nrow=0, ncol=profileLen), row.names=integer() )

    for (taxId in taxIds)
    {
        profile <- readDeltaLFEProfile( taxIds, getH5Filename( taxId ) )

        gcContent <- getGenomicGCContent( taxId )

        #optimumTemperature <- getTemperature( taxId )

        #temperatureRange <- getTemperatureCat( taxId )

        pairedFraction <- getPairedFraction( taxId )
        
        #salinity <- getSalinity( taxId )
        
        #habitat <- getHabitat( taxId )

        #oxygenReq <- getOxygenReq( taxId )
                
        newrow <- data.frame( GenomicGC=c(gcContent), PairedFraction=c(pairedFraction), Profile=t(profile), row.names=c(taxId) )
        
        combined <- rbind(combined, newrow)
    }
    combined
}



getAllTaxIds <- function()
{
    f <- function(s) { as.integer( strsplit(s, "_")[[1]][[4]] ) }
    glob1 <- sprintf("gcdata_v2_taxid_*_profile_%d_%d_%s_%d.h5", profileStop, profileStep, profileReference, profileStart )
    lapply( list.files(pattern=glob2rx(glob1)), f )
}


allTaxIds = getAllTaxIds()

traits <- readAllProfiles( allTaxIds )

stopifnot(nrow(traits) > 300 )
stopifnot(ncol(traits) == profileLen+2 )

#----------------------------------


# Source: http://www2.math.su.se/PATHd8/
# Created using: ~/src/PATHd8/PATHd8 ./data/nmicrobiol201648-s6.txt ./data/nmicrobiol201648-s6.txt.nw.PATHd8.nw
#tree <- read.tree("nmicro_s6_pruned_with_taxids.nw")
tree <- read.tree("test_tree.nw")  # TODO - verify this tree

tree <- drop.tip( tree, c("470", "1280", "4932", "2850") )   # Filter species that will prevent analysis from being performed from the tree

N <- nTips(tree)
print(N)
bmcorr <- corBrownian( phy=tree )

summary(tree)

print("------------------------------")
print(nTips(tree))
print(nrow(traits))

#treeTips <- tipLabels(traitsTree)
treeSpecies <- tree$tip.label

traits <- traits[treeSpecies,]   # Discard traits not found in the tree


# Discard tree tips for species missing data
speciesWithMissingData <- row.names(traits[(is.na(traits["GenomicGC"]) | is.na(traits["Profile.1"])),])
tree <- drop.tip( tree, speciesWithMissingData )   # Filter species that will prevent analysis from being performed from the tree

# Discard trait data for species missing from the tree
treeSpecies <- tree$tip.label
traits <- traits[treeSpecies,]   # Discard traits not found in the tree

treeWithTraits <- phylo4d(tree, tip.data=traits, rownamesAsLabels=TRUE)
stopifnot(hasTipData( treeWithTraits ))


getCorrelogram <- function( tree, traitName, n.points=10, ci.bs=20 )
{
    pcCorr <- phyloCorrelogram( treeWithTraits, trait=traitName, n.points=n.points, ci.bs=ci.bs )
    print("-----------------")
    print(class(pcCorr))
    print(class(pcCorr$res))
    print(dim(pcCorr$res))
    return( pcCorr )
}

getCorrelogramsForProfiles <- function( tree, profiles, n.points=10, ci.bs=5 )
{
    df <- data.frame( ProfileStart=ordered(c()), Correlogram=matrix( rep(0.0, n.points), nrow=0, ncol=n.points) )

    for( profile in profiles )
    {
        pos <- as.character((as.integer(strsplit(profile, ".", fixed=TRUE)[[1]][2]) - 1)*10)
        correlogram <- getCorrelogram( tree, profile, n.points, ci.bs )
        newrow <- data.frame( ProfileStart=ordered(c(pos)), Correlogram=t(correlogram$res[,4]) )
        colnames(newrow) <- c("ProfileStart", vapply(1:n.points, function(x) as.character(x), character(1)))
        print(newrow)
        df <- rbind( df, newrow)
    }

    correlogram <- getCorrelogram( tree, "GenomicGC", n.points, ci.bs )
    newrow <- data.frame( ProfileStart=ordered(c("GenomicGC")), Correlogram=t(correlogram$res[,4]) )
    colnames(newrow) <- c("ProfileStart", vapply(1:n.points, function(x) as.character(x), character(1)))
    print(newrow)
    df <- rbind( df, newrow )

    correlogram <- getCorrelogram( tree, "PairedFraction", n.points, ci.bs )
    newrow <- data.frame( ProfileStart=ordered(c("PairedFraction")), Correlogram=t(correlogram$res[,4]) )
    colnames(newrow) <- c("ProfileStart", vapply(1:n.points, function(x) as.character(x), character(1)))
    print(newrow)
    df <- rbind( df, newrow )

    print(df$ProfileStart)
    
    return(df)
}

traitToString <- function(v, minval=15, maxval=85)
{
    #print(class(v))
    #print(v)
    stopifnot(v >= minval && v <= maxval)
    intval <- round((v-minval)/(maxval-minval)*255)
    return(sprintf("#%02X%02X%02X", intval, intval, intval))
}


plotTraitOnTree <- function(tree, trait, minval=15, maxval=85)
{

    x <- tipData(tree)[,trait]
    #print(class(x))
    y <- vapply(x, traitToString, character(1))
    #print(y)
    #print(class(y))
    #print(length(y))
    tipData(tree)$y <- y
    print(class(tipData(tree)))

    plot(tree, type="phylo", show.tip.label=FALSE, y.lim=100.0 )
    #nodelabels(pch=21, bg=y)
    tiplabels(pch=21, bg=y)
    add.scale.bar()

}

plotTraitOnTree(treeWithTraits, "GenomicGC")



#profiles <- vapply(1:31, function (x) paste("Profile.", as.character(x), sep=""), character(1))
#correl <- getCorrelogramsForProfiles( treeWithTraits, profiles, n.points=50 )
#print(correl)

#print( melt( correl, id.vars=c("ProfileStart") ) )

#p <- ggplot( data=data.frame( melt( correl, id.vars=c("ProfileStart") ) ), aes(y=variable, x=ProfileStart, fill=value) ) +
#    geom_raster() +
#    scale_x_discrete() +
#    theme(axis.text.x=element_text(angle=90, hjust=1)) +
#    scale_fill_gradient2(low="#7777ff", mid="#000000", high="#ff7777", midpoint=0)

#print(p)




dev.off()
