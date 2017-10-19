library("nlme")
library("ggplot2")
library("reshape")
library("rhdf5")
library("rredis")
library("phylobase")
library("ape")
library("adephylo")
library("grid")

#---------------------------------------------------------------
# Configuration

pdf("tree_phenotypes_regression.out.pdf")

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

    stopifnot(x$block0_items[[1]]=="gc")
    gc <- x$block0_values[1,]

    return( native - shuffled )
    #return( gc )
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

getTemperature <- function(taxId)
{
    val <- redisGet(sprintf("species:taxid:%s:properties:optimum-temperature", taxId))
    
    if( !is.null(val)) {
        val <- as.double(val)
        stopifnot( val > -20.0 && val < 120.0 )
        return( val )
    }
    else
    {
        return( NA )  # convert NULL to NA
    }
}

temperatureLevels = c(NA, "Psychrophilic", "Mesophilic", "Thermophilic", "Hyperthermophilic")
getTemperatureCat <- function(taxId)
{
    val <- redisGet(sprintf("species:taxid:%s:properties:temperature-range", taxId))
    
    if( !is.null(val)) {
        #val <- as.double(val)
        #stopifnot( val > -20.0 && val < 120.0 )

        if( val=="Unknown") {
            return( ordered(c(NA), temperatureLevels) )
        }
        return( ordered(c(val), temperatureLevels) )
    }
    else
    {
        return( ordered(c(NA), temperatureLevels ) )  # convert NULL to NA
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

salinityLevels = c(NA, "Mesophilic", "ModerateHalophilic", "ExtremeHalophilic")
getSalinity <- function(taxId)
{
    val <- redisGet(sprintf("species:taxid:%s:properties:salinity", taxId))
    
    if( !is.null(val)) {
        if( val=="NonHalophilic" )
        {
            val <- "Mesophilic"
        }           
        val <- ordered(c(val), salinityLevels)
        return( val )
    }
    else
    {
        return( ordered(c(NA), salinityLevels) )  # convert NULL to NA
    }
}

habitatLevels = c(NA, "Aquatic", "Terrestrial", "HostAssociated", "Specialized", "Multiple")
getHabitat <- function(taxId)
{
    val <- redisGet(sprintf("species:taxid:%s:properties:habitat", taxId))
    
    if( !is.null(val)) {
        #if( val=="NonHalophilic" )
        #{
        #    val <- "Mesophilic"
        #}           
        val <- factor(c(val), habitatLevels)
        return( val )
    }
    else
    {
        return( factor(c(NA), habitatLevels) )  # convert NULL to NA
    }
}

oxygenLevels = c(NA, "Aerobic", "Facultative", "Anaerobic")
getOxygenReq <- function(taxId)
{
    val <- redisGet(sprintf("species:taxid:%s:properties:oxygen-req", taxId))
    
    if( !is.null(val)) {
        if( val=="Microaerophilic" )
        {
            val <- NA # TODO
        }           
        val <- factor(c(val), oxygenLevels)
        return( val )
    }
    else
    {
        return( factor(c(NA), oxygenLevels) )  # convert NULL to NA
    }
}


readAllProfiles <- function( taxIds )
{

    combined <- data.frame( GenomicGC=double(), OptimumTemp=double(), TemperatureRange=ordered(c(), temperatureLevels), PairedFraction=double(), Profile=matrix( rep(0.0, profileLen), nrow=0, ncol=profileLen), Salinity=ordered(c(), salinityLevels), Habitat=factor(c(), habitatLevels), OxygenReq=ordered(c(), oxygenLevels), row.names=integer() )

    for (taxId in taxIds)
    {
        profile <- readDeltaLFEProfile( taxIds, getH5Filename( taxId ) )

        gcContent <- getGenomicGCContent( taxId )

        optimumTemperature <- getTemperature( taxId )

        temperatureRange <- getTemperatureCat( taxId )

        pairedFraction <- getPairedFraction( taxId )
        
        salinity <- getSalinity( taxId )
        
        habitat <- getHabitat( taxId )

        oxygenReq <- getOxygenReq( taxId )
                
        newrow <- data.frame( GenomicGC=c(gcContent), OptimumTemp=c(optimumTemperature), TemperatureRange=temperatureRange, PairedFraction=c(pairedFraction), Profile=t(profile), Salinity=salinity, Habitat=habitat, OxygenReq=oxygenReq, row.names=c(taxId) )
        
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
stopifnot(ncol(traits) == profileLen+7 )


#x["511145",]
#x[as.character(511145),]


#traits["580340",]
#x[c(580340),]
#print(x)

print("###################")
print(sum(is.na(traits$GenomicGC)))
print(sum(is.na(traits$Profile.15)))


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
#print(nrow(traits))


#-------------------------------------------------------------------------------------
# Draw VCV matrix

vcv1 <- vcv( phy=tree, model="Brownian")
dimnames(vcv1) <- list(1:ncol(vcv1))
bmcorrmtx <- data.frame( melt( vcv1 ) )
#print(range(bmcorrmtx$X1))
#print(range(bmcorrmtx$X2))
#print(range(bmcorrmtx$value))
#dimnames(bmcorrmtx)
#p <- ggplot( data=bmcorrmtx, aes(x=X1, y=X2, fill=value) ) + geom_raster() + scale_y_reverse(); p



performGLSregression <- function( traits, tree, Xtrait, Ytrait, plotRegression=TRUE )
{
    # Discard tree tips for species missing data
    speciesWithMissingData <- row.names(traits[(is.na(traits[Xtrait]) | is.na(traits[Ytrait])),])
    tree <- drop.tip( tree, speciesWithMissingData )   # Filter species that will prevent analysis from being performed from the tree

    # Discard trait data for species missing from the tree
    treeSpecies <- tree$tip.label
    traits <- traits[treeSpecies,]   # Discard traits not found in the tree

    # Check tree <-> data-frame correlation appears valid
    stopifnot(hasTipData( phylo4d(tree, tip.data=traits, rownamesAsLabels=TRUE) ) )

    N <- nTips(tree)
    print(N)

    # Calculate the correlation matrix (for a BM process)
    bmcorr <- corBrownian( phy=tree )

    YvsX = as.formula( paste(Ytrait, " ~ ", Xtrait) )
    #print(YvsX)

    m1 <- lm(  YvsX, traits); summary(m1)
    co <- coef(m1)
    print(co)

    #print("//////////////")
    #print( nrow(traits) )
    #print( nTips(traitsTree) )

    clr2 <- "#55CCB0"

    gls1 <- gls( YvsX, traits, correlation=bmcorr, na.action=na.omit); summary(gls1)
    #gls1 <- pGLS( YvsX, traits, bmcorr, na.action=na.omit); summary(gls1)
    co2 <- coef(gls1)

    if( plotRegression )
    {
        if( any(class(traits[,Xtrait])==c("factor")) )
        {
            p <- ggplot(traits, aes(get(Xtrait), get(Ytrait))) +
              labs(y=Ytrait, x=Xtrait) +
              geom_boxplot() +
              geom_jitter() +
              geom_hline( yintercept = 0 ) +
              annotate( "text", x=Inf,  y=Inf, label=sprintf("lm\nR^2 = %4g\nP-val = %4g\nn = %d", summary(m1)$r.squared,      summary(m1)$coefficients[2,4], nrow(traits) ), hjust=1, vjust=1, colour="red" ) +
              annotate( "text", x=-Inf, y=Inf, label=sprintf("gls\n\nP-val = %4g\nn = %d", summary(gls1)$tTable[2,4],     nrow(traits) ), hjust=0, vjust=1, colour=clr2 );
              print(p)
        }
        else
        {
            p <- ggplot(traits, aes(get(Xtrait), get(Ytrait))) +
              labs(y=Ytrait, x=Xtrait) +
              geom_point() +
              geom_hline( yintercept = 0 ) +
              geom_abline( aes( slope=co[Xtrait],  intercept=co["(Intercept)"]), colour="red"  ) +
              geom_abline( aes( slope=co2[Xtrait], intercept=co2["(Intercept)"]), colour=clr2  ) +
              annotate( "text", x=Inf,  y=Inf, label=sprintf("lm\nR^2 = %4g\nP-val = %4g\nn = %d", summary(m1)$r.squared,      summary(m1)$coefficients[2,4], nrow(traits) ), hjust=1, vjust=1, colour="red" ) +
              annotate( "text", x=-Inf, y=Inf, label=sprintf("gls\n\nP-val = %4g\nn = %d", summary(gls1)$tTable[2,4],     nrow(traits) ), hjust=0, vjust=1, colour=clr2 );
              print(p)
        }
    }

    #return( c( summary(gls1)$tTable[2,4], co2[Xtrait] ) )
    return( c( summary(gls1)$tTable[2,4], summary(m1)$r.squared ) )
}

performGLSregression_profileRangeMean <- function( traits, tree, Xtrait, profileRange )
{
    # Make list of the selected profile columns
    variables <- vapply(profileRange[1]:profileRange[2], function(x) paste('Profile.', as.character(x), sep=""), character(1))

    ###########################################################################################
    ## DEBUG ONLY ### DEBUG ONLY ### DEBUG ONLY ### DEBUG ONLY ### DEBUG ONLY ### DEBUG ONLY ##
    ###########################################################################################
    #if( profileRange[1]==1 && profileRange[2]==5 )
    #{
    #    return( c(0.05, 0.0) )
    #}
    ###########################################################################################
    ## DEBUG ONLY ### DEBUG ONLY ### DEBUG ONLY ### DEBUG ONLY ### DEBUG ONLY ### DEBUG ONLY ##
    ###########################################################################################

    # Extract the select columns and compute the mean values
    # TODO - learn dplyr...
    traits$RangeMean <- rowMeans( traits[variables] )

    # Perform
    return( performGLSregression( traits, tree, Xtrait, "RangeMean", plotRegression=FALSE ) )
}

glsRegressionRangeAnalysis <- function( traits, tree, Xtrait, profileRange )
{
    # Generate indices for all valid subranges in the given range
    stopifnot(profileRange[2]>=profileRange[1])
    range <- profileRange[1]:profileRange[2]
    subranges <- expand.grid( range, range )
    subranges <- subranges[ subranges$Var2 >= subranges$Var1, ]
    stopifnot(all(subranges$Var2 - subranges$Var1 >= 0))
    stopifnot(nrow(subranges) == length(range)*(length(range)+1)/2 ) # arithmetic series sum

    # Perform regression on each subrange
    regressionResults <- apply( subranges, 1, function(v) performGLSregression_profileRangeMean(traits, tree, Xtrait, v) )
    subranges$Pvalue <- regressionResults[1,]
    subranges$Coef   <- regressionResults[2,]
    subranges$Logpval <- log10( subranges$Pvalue )
    subranges$Var1 <- (subranges$Var1 - 1) * 10
    subranges$Var2 <- (subranges$Var2 - 1) * 10

    minLogpval <- min( subranges$Logpval, -3 )

    #print(subranges)

    plotTitle = sprintf("%s effect on dLFE", Xtrait) 

    grid.newpage()
    plotdf <- melt( subranges, id.vars=c("Var1", "Var2"), measure.vars=c("Logpval") )
    p <- ggplot( data=plotdf, aes(x=Var1, y=Var2, fill=value) ) +
        geom_raster() +
        labs(y="Profile End", x="Profile Start", value="log(p-value)", title=plotTitle) +
        scale_fill_gradient2(low="#cc7777", mid="#ffff88", high="#000000", limits=c(min(subranges$Logpval), 0), midpoint=minLogpval ) + 
        coord_fixed() +
        geom_point(data=plotdf[plotdf$value<=-2.0,], aes(x=Var1,y=Var2), position=position_nudge(y=4, x=-4), color="white", size=0.2 ) +
        theme( plot.margin=unit(c(2,2,2,2), "cm"), aspect.ratio=1 )
    print( p, vp=viewport(angle=-45) )

    grid.newpage()
    p <- ggplot( data=data.frame( melt( subranges, id.vars=c("Var1", "Var2"), measure.vars=c("Coef"   ) ) ), aes(x=Var1, y=Var2, fill=value) ) +
        geom_raster() +
        labs(y="Profile End", x="Profile Start", value="R^2", title=plotTitle) +
        scale_fill_gradient2(low="#7777ff", mid="#000000", high="#ff7777", limits=c(-1, 1)) + 
        coord_fixed() +
        geom_point(data=plotdf[plotdf$value<=-2.0,], aes(x=Var1,y=Var2), position=position_nudge(y=4, x=-4), color="white", size=0.2 ) +
        theme( plot.margin=unit(c(2,2,2,2), "cm"), aspect.ratio=1 )
    print( p, vp=viewport(angle=-45) )

}

#performGLSregression( traits, tree, "GenomicGC",   "Profile.15" )
#print(performGLSregression_profileRangeMean( traits, tree, "GenomicGC", c(15,15) ))
#print(performGLSregression_profileRangeMean( traits, tree, "GenomicGC", c(15,19) ))

pyramidLength <- 31

#--------------------
#glsRegressionRangeAnalysis( traits, tree, "OptimumTemp", c(1,pyramidLength))
#performGLSregression( traits, tree, "OptimumTemp", "Profile.11" )
#performGLSregression( traits, tree, "OptimumTemp", "GenomicGC" )

#glsRegressionRangeAnalysis( traits, tree, "TemperatureRange", c(1,pyramidLength))
#performGLSregression( traits, tree, "TemperatureRange", "Profile.11" )
#performGLSregression( traits, tree, "TemperatureRange", "GenomicGC" )


#dev.off()
#quit()
#--------------------


glsRegressionRangeAnalysis( traits, tree, "GenomicGC", c(1,pyramidLength))

glsRegressionRangeAnalysis( traits, tree, "OptimumTemp", c(1,pyramidLength))

glsRegressionRangeAnalysis( traits, tree, "OxygenReq", c(1,pyramidLength))

glsRegressionRangeAnalysis( traits, tree, "Habitat", c(1,pyramidLength))

glsRegressionRangeAnalysis( traits, tree, "Salinity", c(1,pyramidLength))

glsRegressionRangeAnalysis( traits, tree, "PairedFraction", c(1,pyramidLength))

performGLSregression( traits, tree, "OptimumTemp", "Profile.15" )
performGLSregression( traits, tree, "GenomicGC",   "Profile.1" )
performGLSregression( traits, tree, "OptimumTemp", "Profile.15" )

performGLSregression( traits, tree, "PairedFraction", "Profile.15" )
performGLSregression( traits, tree, "PairedFraction", "Profile.1" )

performGLSregression( traits, tree, "TemperatureRange", "Profile.15" )
performGLSregression( traits, tree, "Salinity",         "Profile.15" )
performGLSregression( traits, tree, "OxygenReq",        "Profile.14" )
performGLSregression( traits, tree, "Habitat",          "Profile.15" )
performGLSregression( traits, tree, "Habitat",          "Profile.1" )

performGLSregression( traits, tree, "TemperatureRange", "GenomicGC" )
performGLSregression( traits, tree, "Salinity",         "GenomicGC" )
performGLSregression( traits, tree, "OxygenReq",        "GenomicGC" )
performGLSregression( traits, tree, "Habitat",          "GenomicGC" )
performGLSregression( traits, tree, "PairedFraction",   "GenomicGC" )


dev.off()