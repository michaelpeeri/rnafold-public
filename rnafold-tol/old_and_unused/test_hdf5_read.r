library("rhdf5")
library("rredis")
library("phylobase")
library("ape")
library("ggplot2")
library("nlme")


redisConnect(host="power5", password="rnafold")

profileStart <- 0
profileStop <- 150
profileStep <- 10
profileReference <- "begin"
profileLen <- (profileStop-profileStart)/profileStep


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

    native - shuffled
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


readAllProfiles <- function( taxIds )
{
    #combined <- data.frame( GenomicGC=double(), profile.1=double(), profile.2=double(), profile.3=double(), profile.4=double(), profile.5=double(), profile.6=double(), profile.7=double(), profile.8=double(), profile.9=double(), profile.10=double(), profile.11=double(), profile.12=double(), profile.13=double(), profile.14=double(), profile.15=double(), row.names=factor() )

    combined <- data.frame( GenomicGC=double(), profile=matrix( rep(0.0, profileLen), nrow=0, ncol=profileLen), row.names=factor() )
    
    for (taxId in taxIds)
    {
        profile <- readDeltaLFEProfile( taxIds, getH5Filename( taxId ) )

        gcContent <- getGenomicGCContent( taxId )
                
        newrow <- data.frame( Profile=t(profile), GenomicGC=t(gcContent), row.names=factor(c(taxId)) )
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

x <- readAllProfiles( allTaxIds )

stopifnot(nrow(x) > 300 )
stopifnot(ncol(x) == profileLen+1 )


#x["511145",]
#x[as.character(511145),]


#----------------------------------


# Source: http://www2.math.su.se/PATHd8/
# Created using: ~/src/PATHd8/PATHd8 ./data/nmicrobiol201648-s6.txt ./data/nmicrobiol201648-s6.txt.nw.PATHd8.nw
#tree <- read.tree("nmicro_s6_pruned_with_taxids.nw")
tree <- read.tree("test_tree.nw")
N <- nTips(tree)
print(N)
bmcorr <- corBrownian( phy=tree )

summary(tree)


