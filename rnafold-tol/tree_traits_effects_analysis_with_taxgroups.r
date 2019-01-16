# Perform GLS analysis for different traits and different taxonomic groups

library("Matrix")
library("nlme")
library("ggplot2")
library("reshape")
library("scales")
library("rhdf5")
library("rredis")
#library("redux")
library("phylobase")
library("ape")
#library("adephylo")
library("phytools")  # for read.newick
library("phylosignal")
library("grid")
library("gridExtra")
library("dplyr")
library("rlang") # for duplicate
library("kde1d")
library("minerva")
library("lmtest") # for Breush-Pagan test (homoskedasticity)
library("Metrics")

#---------------------------------------------------------------
# Configuration

# Source: http://www2.math.su.se/PATHd8/
# Created using: ~/src/PATHd8/PATHd8 ./data/nmicrobiol201648-s6.txt ./data/nmicrobiol201648-s6.txt.nw.PATHd8.nw
inputTree <- "nmicro_s6_pruned_with_taxids.nw"

#redis <- redux::hiredis(host="power5", password="rnafold")
redisConnect(host="compute-0-224", password="rnafold")

#profileStart <- 0
#profileStop <- 1000
profileStep <- 10
#profileReference <- "begin"
#profileLen <- (profileStop-profileStart)/profileStep
seriesNumber     <- 102

#-------------------------------
# Start profile only
#profileSpecifications <- list( list( start=0, stop=1000, step=profileStep, reference="begin", series=102, len=100 ) )
#-------------------------------
# Start- and End-profiles
profileSpecifications <- list( list( start=0, stop=310, step=profileStep, reference="begin", series=102, len=31 ),
                               list( start=0, stop=310, step=profileStep, reference="end",   series=102, len=32 ) )
#-------------------------------
# Start profile with native-only and shuffled-only components
#profileSpecifications <- list( list( start=0, stop=310, step=profileStep, reference="begin",          series=102, len=31 ),
#                               list( start=0, stop=310, step=profileStep, reference="begin:native",   series=102, len=31 ),
#                               list( start=0, stop=310, step=profileStep, reference="begin:shuffled", series=102, len=31 ) )
#-------------------------------
#profileLen <- 30
                              
taxidToKingdomFilename <- "TaxidToKingdom.csv"   # create file using: python2 create_taxid_kingdom_table.py


#pyramidLength <- 30+1  # unit is single windows (not real nt when profileStep>1)
#pyramidLength <- 30+1  # unit is single windows (not real nt when profileStep>1)

#pyramidLength <- 99+1  # unit is single windows (not real nt when profileStep>1)
pyramidLength <- 30+1  # unit is single windows (not real nt when profileStep>1)

maxPyramidHeight <- 1-1     # calculate single-windows only
#maxPyramidHeight <- 32-1   # calculate full pyramid

profileMode <- "dLFE"
#profileMode <- "nativeLFE"
#profileMode <- "shuffledLFE"
#profileMode <- "abs(dLFE)"
#profileMode <- "gc"

loadDeltas <- FALSE

#                                          -r1-  -r2-
groupsTableOutputFile.limitRangeFromNt <- c( 0,   70)   # range is inclusive at either end
groupsTableOutputFile.limitRangeToNt   <- c(20, 1000)   # range is inclusive at either end
#groupsTableOutputFile.limitRangeFromNt <- c( 0)   # range is inclusive at either end
#groupsTableOutputFile.limitRangeToNt   <- c(10)   # range is inclusive at either end
stopifnot( length(groupsTableOutputFile.limitRangeFromNt) == length(groupsTableOutputFile.limitRangeToNt))
stopifnot(groupsTableOutputFile.limitRangeFromNt <= groupsTableOutputFile.limitRangeToNt )
stopifnot( c( groupsTableOutputFile.limitRangeFromNt, 1e8 ) > c( -1, groupsTableOutputFile.limitRangeToNt ) ) # ranges shouldn't overlap
stopifnot( profileStep*(pyramidLength-1) >= groupsTableOutputFile.limitRangeFromNt[length(groupsTableOutputFile.limitRangeFromNt)] ) # pyramidLength must overlap with limits (otherwise peaks csv output will be empty)
groupsTableOutputFile <- sprintf("tree_traits_effects_analysis_with_taxgroups.out.%s.length.%d.%%s.csv", profileMode, profileStep*(pyramidLength-1) )  # peaks CSV output file
# Add range 0 (full range)
groupsTableOutputFile.limitRangeFromNt <- c(0  , groupsTableOutputFile.limitRangeFromNt)
groupsTableOutputFile.limitRangeToNt   <- c(1e6, groupsTableOutputFile.limitRangeToNt  )

lipaMoranReps <- 2000
phylosignalLipaOutputFile <- "phylosignalLipaOutputFile.csv"
lipaUseRankedValues <- TRUE


significanceLevel <- 1e-2

MIC.pval.num.iterations <- 10  # the default number of permutations used for estimating p-value for MIC
#MIC.pval.num.iterations <- 50000  # the default number of permutations used for estimating p-value for MIC

rankTransformResponse <- FALSE

scale.dLFE <- 2.85

max.num.species.to.label <- 80  # if the number of samples is less than this, label each one with the species' nickname


getProfilePositions <- function( profileId=1 )
{
    spec <- profileSpecifications[[profileId]]
    if( spec$reference %in% c("begin", "begin:native", "begin:shuffled" ) )
    {
        return( seq( spec$start, spec$stop-1, spec$step ) )
    }
    else if ( spec$reference == "end" )
    {
        #return( seq( spec$stop-spec$step, spec$start, -(spec$step) ) )
        return(  seq( -spec$stop         , spec$start, spec$step ) )
    }
    else
    {
        # unsupported profile reference
        stopifnot(FALSE);
    }
}

getProfileReference <- function( profileId )
{
    return( profileSpecifications[[profileId]]$reference )
}



getProfileVariables <- function( range, profileId=1 )
{
    stopifnot( length(range)==2 )
    stopifnot( !is.na(profileId) )
    return( vapply(range[1]:range[2], function(x) sprintf("Profile_%d.%d", profileId, x), character(1)) )
}


# Big plot scale (tiny font)
#fontScale <- 8
#cdsScaleResolution <- 35

# Standard scale
fontScale <- 12
cdsScaleResolution <- 50

## Tiny scale (large fonts)
#fontScale <- 16
#cdsScaleResolution <- 150




minimalTaxonSize <- 9 # Note - should match the value in create_taxid_kingdom_table.py
#minimalTaxonSize <- 100 # Testing only


#cairo_pdf(sprintf("tree_phenotypes_regression.out.%s.by.taxgroups.ranked.%%d.pdf", profileMode))
cairo_pdf(sprintf("tree_phenotypes_regression.out.%s.by.taxgroups.%%d.pdf", profileMode))


speciesBlacklist <- c("470", "1280", "4932", "2850", "508771", "195065", "641309")  # Filter species that will prevent analysis from being performed from the tree because of errors or missing data
# TODO - try to fix these...

#---------------------------------------------------------------

profileFilenameMapping <- c("begin"="begin", "end"="end", "begin:native"="begin", "begin:shuffled"="begin")

getH5Filename <- function( taxId, profileSpec )
{
    # Note: this heuristic will have to be changed when the files are recreated...
    if( profileSpec$stop < 350 )
    {
        seriesId <- "_t11"
    }
    else
    {
        seriesId <- ""
    }
    ##seriesId <- "_t11"
    
    return( sprintf("gcdata_v2_taxid_%d_profile_%d_%d_%s_%d%s.h5", taxId, profileSpec$stop, profileSpec$step, profileFilenameMapping[profileSpec$reference], profileSpec$start, seriesId) )
}

# Read native and shuffled LFE profiles from the hdf files
# returns numeric vector.
readDeltaLFEProfile <- function( taxId, h5filename, profileSpec )
{
    if( !file.exists(h5filename) )
    {
        print(sprintf("Warning: file not found: %s", h5filename))
        return( NA );
    }
    
    y <- h5ls(h5filename)
    # Read the profile data
    key <- y[which(startsWith(y[,1], "/df_"))[1], 1]   # find the first element that starts with '/df_'

    x <- h5read(h5filename, key )

    stopifnot(x$block0_items[[2]]=="native")
    native <- x$block0_values[2,]

    stopifnot(x$block0_items[[3]]=="shuffled")
    shuffled <- x$block0_values[3,]

    gc <- NA
    if( profileMode == "gc" )
    {
        stopifnot(x$block0_items[[1]]=="gc")
        gc <- x$block0_values[1,]
    }

    diag.d <- NA
    
    if( loadDeltas )
    {
        # Read the raw differences (for Wilcoxon test)
        key <- y[which(startsWith(y[,1], "/deltas_"))[1], 1]   # find the first element that starts with '/deltas_'
        d <- h5read(h5filename, key )

        #print(d$block0_items)
        #print(d$block1_items)
        stopifnot(d$block0_items[[1]]=="delta")
        deltas <- d$block0_values[1,]
        #print(length(deltas))
        stopifnot(d$block1_items[[1]]=="pos")
        poss <- d$block1_values[1,]
        #print(length(poss))

        if( length(deltas)==length(poss))
        {
            data <- data.frame(delta=deltas, pos=poss)
            print("--")
            p0 <- min(poss)
            p1 <- max(poss)
            N <- sort(unique(poss))
            diag.d <- data.frame( pos=N, count=integer(length(N)), p.val=double(length(N)), mean.delta=double(length(N)) )

            for( i in N )
            {
                sample <- data[data$pos==i,]
                sample <- sample[!is.na(sample$delta),]
                #print( sprintf( "%d] %f", i, mean(data[data$pos==i,]$delta) ) )
                #print( length(sample$delta) )
                #print( mean(sample$delta) )
                diag.d[diag.d$pos==i,]$count <- length(sample$delta)
                diag.d[diag.d$pos==i,]$p.val = wilcox.test(sample$delta, alternative="two.sided")$p.value
                diag.d[diag.d$pos==i,]$mean.delta = mean(sample$delta)

                diag.d$realval  <- native - shuffled
                diag.d$native   <- native
                diag.d$shuffled <- shuffled
            }

            if( length(deltas) > 0 )
            {
                #d <- data.frame(melt( diag.p, id.vars=c("pos") ))
                #print(colnames(d))

                p <- ggplot( data=diag.d, aes(x=pos) ) +
                    geom_line( aes(y=log10(p.val)), color="red" ) +
                    geom_hline( yintercept=-2 )
                print(p)

                p <- ggplot( data=diag.d, aes(x=pos) ) +
                    geom_line( aes(y=count), color="grey" )
                print(p)

                p <- ggplot( data=diag.d, aes(x=pos) ) +
                    geom_line( aes(y=mean.delta), color="grey",  size=2,   alpha=0.3 ) +
                    geom_line( aes(y=realval),    color="black", size=0.6 ) +
                    geom_line( aes(y=native),     color="green" ) +
                    geom_line( aes(y=shuffled),   color="red" ) +
                    geom_hline( yintercept=0 ) +
                    scale_y_continuous( limits=c(-2.84,2.84) ) +
                    labs(title=taxId)
                print(p)
            }
        }
    }

    # Return values as configured...
    if(   profileMode == "dLFE" && profileSpec$reference %in% c("begin", "end") )
    {
        vals <- native - shuffled
    }
    else if( profileMode == "nativeLFE" && profileSpec$reference %in% c("begin", "end") )
    {
        vals <- native 
    }
    else if( profileMode == "shuffledLFE" && profileSpec$reference %in% c("begin", "end") )
    {
        vals <-  shuffled
    }
    else if( profileMode == "abs(dLFE)" && profileSpec$reference %in% c("begin", "end") )
    {
        vals <- abs(native - shuffled)
    }
    else if( profileMode == "gc" && profileSpec$reference %in% c("begin", "end") )
    {
        vals <- gc
    }
    else if( profileSpec$reference=="begin:native" )
    {
        vals <- native
    }
    else if( profileSpec$reference=="begin:shuffled" )
    {
        vals <- shuffled
    }
    else
    {
        stopifnot(FALSE)
        vals <- NA 
    }
    return( list(vals, diag.d) )
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

temperatureLevels <- c(NA, "Psychrophilic", "Mesophilic", "Thermophilic", "Hyperthermophilic")
getTemperatureCat <- function(taxId)
{
    val <- redisGet(sprintf("species:taxid:%s:properties:temperature-range", taxId))
    
    if( !is.null(val)) {
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

salinityLevels <- c(NA, "Mesophilic", "ModerateHalophilic", "ExtremeHalophilic")
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

habitatLevels <- c(NA, "Aquatic", "Terrestrial", "HostAssociated", "Specialized", "Multiple")
getHabitat <- function(taxId)
{
    val <- redisGet(sprintf("species:taxid:%s:properties:habitat", taxId))
    
    if( !is.null(val)) {
        val <- factor(c(val), habitatLevels)
        return( val )
    }
    else
    {
        return( factor(c(NA), habitatLevels) )  # convert NULL to NA
    }
}

oxygenLevels <- c(NA, "Aerobic", "Facultative", "Microaerophilic", "Anaerobic")
getOxygenReq <- function(taxId)
{
    val <- redisGet(sprintf("species:taxid:%s:properties:oxygen-req", taxId))
    
    if( !is.null(val)) {
        #if( val=="Microaerophilic" )
        #{
        #    val <- NA # TODO
        #}           
        val <- factor(c(val), oxygenLevels)
        return( val )
    }
    else
    {
        return( factor(c(NA), oxygenLevels) )  # convert NULL to NA
    }
}

getGenomeSizeMb <- function(taxId)
{
    val <- redisGet(sprintf("species:taxid:%s:properties:genome-size-mb", taxId))
    
    if( !is.null(val)) {
        val <- as.double(val)
        stopifnot( val > 0.1 && val < 1000.0 )
        return( val )
    }
    else
    {
        return( NA )  # convert NULL to NA
    }
}

getProteinCount <- function(taxId)
{
    val <- redisGet(sprintf("species:taxid:%s:properties:protein-count", taxId))
    
    if( !is.null(val)) {
        val <- as.integer(val)
        stopifnot( val > 100 && val < 100000 )
        return( val )
    }
    else
    {
        return( NA )  # convert NULL to NA
    }
}

getGrowthTimeHours <- function(taxId)
{
    val <- redisGet(sprintf("species:taxid:%s:properties:growth-time-hours-v2", taxId))
    
    if( !is.null(val)) {
        val <- as.double(val)
        stopifnot( val > 0.05 && val < 1000.0 )
        return( val )
    }
    else
    {
        return( NA )  # convert NULL to NA
    }
}

getGenomicENc <- function(taxId)
{
    val <- redisGet(sprintf("species:taxid:%s:properties:ENc", taxId))
    
    if( !is.null(val)) {
        val <- as.double(val)
        stopifnot( val > 10.0 && val < 75.0 )  # The extreme values for ENc are not clear to me...
        return( val )
    }
    else
    {
        return( NA )  # convert NULL to NA
    }
}

getGenomicENc.prime <- function(taxId)
{
    val <- redisGet(sprintf("species:taxid:%s:properties:ENc-prime", taxId))
    
    if( !is.null(val)) {
        val <- as.double(val)
        stopifnot( val > 10.0 && val < 75.0 )  # The extreme values for ENc are not clear to me...
        return( val )
    }
    else
    {
        return( NA )  # convert NULL to NA
    }
}

getGenomicCAI <- function(taxId)
{
    val <- redisGet(sprintf("species:taxid:%s:properties:genomic-CAI", taxId))
    
    if( !is.null(val)) {
        val <- as.double(val)
        stopifnot( val >= 0.0 && val <= 1.0 )
        return( val )
    }
    else
    {
        return( NA )  # convert NULL to NA
    }
}

getGenomicCBI <- function(taxId)
{
    val <- redisGet(sprintf("species:taxid:%s:properties:genomic-CBI", taxId))
    
    if( !is.null(val)) {
        val <- as.double(val)
        stopifnot( val >= -0.5 && val <= 1.0 )
        return( val )
    }
    else
    {
        return( NA )  # convert NULL to NA
    }
}

getGenomicFop <- function(taxId)
{
    val <- redisGet(sprintf("species:taxid:%s:properties:genomic-Fop", taxId))
    
    if( !is.null(val)) {
        val <- as.double(val)
        stopifnot( val >= 0.0 && val <= 1.0 )
        return( val )
    }
    else
    {
        return( NA )  # convert NULL to NA
    }
}

getGenomicNc <- function(taxId)
{
    val <- redisGet(sprintf("species:taxid:%s:properties:genomic-Nc-codonw", taxId))
    
    if( !is.null(val)) {
        val <- as.double(val)
        stopifnot( val >= 10.0 && val <= 75.0 )
        return( val )
    }
    else
    {
        return( NA )  # convert NULL to NA
    }
}


getGenomicDCBS <- function(taxId)
{
    val <- redisGet(sprintf("species:taxid:%s:properties:DCBS-geomean", taxId))
    
    if( !is.null(val)) {
        val <- as.double(val)
        stopifnot( val > -1.0 && val < 5.0 )  # I don't know what are the extreme values for DCBS
        return( val )
    }
    else
    {
        return( NA )  # convert NULL to NA
    }
}


readAllProfiles <- function( taxIds, profileSpecifications )
{

    args <- list( GenomicGC=double(), OptimumTemp=double(), TemperatureRange=ordered(c(), temperatureLevels), PairedFraction=double(), Salinity=ordered(c(), salinityLevels), Habitat=factor(c(), habitatLevels), OxygenReq=factor(c(), oxygenLevels), GenomeSizeMb=double(), LogGenomeSize=double(), ProteinCount=integer(), GrowthTimeHours=double(), LogGrowthTime=double(), GenomicENc=double(), GenomicENc.prime=double(), GenomicCAI=double(), GenomicCBI=double(), GenomicFop=double(), GenomicNc=double(), GenomicDCBS=double(), row.names=integer() )
    for (i in 1:length(profileSpecifications) )
    {
        newargs <- lapply(1:profileSpecifications[[i]]$len, function (i) numeric(0) )
        names(newargs) <- sapply(1:length(newargs), function (j) { sprintf("Profile_%d.%d", i, j) } )
        #print("/////\\\\\\")
        #print(newargs)
        
        args <- c( args, newargs )
    }
    #print(args)
    #print("///1///")
    combined <- do.call( data.frame, args )
    #newrow <- data.frame( GenomicGC=c(gcContent), OptimumTemp=c(optimumTemperature), TemperatureRange=temperatureRange, PairedFraction=c(pairedFraction), Profile=t(profile), Salinity=salinity, Habitat=habitat, OxygenReq=oxygenReq, row.names=c(taxId) )
        
        

    #combined <- data.frame( GenomicGC=double(), OptimumTemp=double(), TemperatureRange=ordered(c(), temperatureLevels), PairedFraction=double(), Profile=matrix( rep(0.0, profileLen), nrow=0, ncol=profileLen), Salinity=ordered(c(), salinityLevels), Habitat=factor(c(), habitatLevels), OxygenReq=factor(c(), oxygenLevels), row.names=integer() )

    for (taxId in taxIds)
    {
        profiles <- lapply( profileSpecifications, function(profileSpec) { readDeltaLFEProfile( taxIds, getH5Filename( taxId, profileSpec ), profileSpec )[[1]]; } )
        gcContent <- getGenomicGCContent( taxId )

        optimumTemperature <- getTemperature( taxId )

        temperatureRange <- getTemperatureCat( taxId )

        pairedFraction <- getPairedFraction( taxId )
        
        salinity <- getSalinity( taxId )
        
        habitat <- getHabitat( taxId )

        oxygenReq <- getOxygenReq( taxId )

        genomeSizeMb <- getGenomeSizeMb( taxId )

        proteinCount <- getProteinCount( taxId )

        growthTimeHours <- getGrowthTimeHours( taxId )

        genomicENc <- getGenomicENc( taxId )
        
        genomicENc.prime <- getGenomicENc.prime( taxId )

        genomicCAI <- getGenomicCAI( taxId )
        
        genomicCBI <- getGenomicCBI( taxId )
        
        genomicFop <- getGenomicFop( taxId )

        genomicNc <- getGenomicNc( taxId )

        genomicDCBS <- getGenomicDCBS( taxId )
        

        #args <- list( GenomicGC=c(gcContent), OptimumTemp=c(optimumTemperature), TemperatureRange=temperatureRange, PairedFraction=c(pairedFraction), Profile=t(profile), Salinity=salinity, Habitat=habitat, OxygenReq=oxygenReq, row.names=c(taxId) )
        args <- list( GenomicGC=c(gcContent), OptimumTemp=c(optimumTemperature), TemperatureRange=temperatureRange, PairedFraction=c(pairedFraction), Salinity=salinity, Habitat=habitat, OxygenReq=oxygenReq, GenomeSizeMb=genomeSizeMb, LogGenomeSize=log(genomeSizeMb), ProteinCount=proteinCount, GrowthTimeHours=c(growthTimeHours), LogGrowthTime=log(growthTimeHours), GenomicENc=c(genomicENc), GenomicENc.prime=c(genomicENc.prime), GenomicCAI=c(genomicCAI), GenomicCBI=c(genomicCBI), GenomicFop=c(genomicFop), GenomicNc=c(genomicNc), GenomicDCBS=c(genomicDCBS), row.names=c(taxId) )
        for (i in 1:length(profileSpecifications) )
        {
            #print("--//--//--")
            #print(i)
            #print(length(profiles[[i]]))
            
            if( !isTRUE(is.na(profiles[[i]][1])) )
            {
                #print("((a))")
                newargs <- as.list(profiles[[i]])
            }
            else
            {
                #print("((b))")
                newargs <- as.list(rep(NA, profileSpecifications[[i]]$len))
            }
                
            names(newargs) <- sapply(1:length(newargs), function (j) { sprintf("Profile_%d.%d", i, j) } )

            args <- c( args, newargs )
        }
        #print(args)

        #print("///2///")
        newrow <- do.call( data.frame, args )
                                        #newrow <- data.frame( GenomicGC=c(gcContent), OptimumTemp=c(optimumTemperature), TemperatureRange=temperatureRange, PairedFraction=c(pairedFraction), Profile=t(profile), Salinity=salinity, Habitat=habitat, OxygenReq=oxygenReq, row.names=c(taxId) )
        #print(colnames(newrow))
        #print(colnames(combined))
        #print(nrow(newrow))
        #print(ncol(newrow))
        #print(ncol(combined))
        #print(taxId)
        
        

        #print(colnames(combined))
                                        #print(colnames(newrow))
        combined <- rbind(combined, newrow)
    }
    #print(names(combined))
    return(combined)
}



getAllTaxIds <- function()
{
    f <- function(s) { as.integer( strsplit(s, "_")[[1]][[4]] ) }
    glob1 <- sprintf("gcdata_v2_taxid_*_profile_%d_%d_%s_%d_t11.h5", profileSpecifications[[1]]$stop, profileSpecifications[[1]]$step, profileSpecifications[[1]]$reference, profileSpecifications[[1]]$start )
    #glob1 <- sprintf("gcdata_v2_taxid_*_profile_%d_%d_%s_%d.h5", profileSpecifications[[1]]$stop, profileSpecifications[[1]]$step, profileSpecifications[[1]]$reference, profileSpecifications[[1]]$start )   # work-around for old 1000nt profiles...
    lapply( list.files(pattern=glob2rx(glob1)), f )
}


#=========================================================================================
#================= Add taxon membership indicators to the traits dataset =================
#=========================================================================================

# Read the taxon memberships table
taxidToKingdom <- read.table(taxidToKingdomFilename, sep=",", header=TRUE )    # create file using: python2 create_taxid_kingdom_table.py
taxidToKingdom$Member_all_1 <- 1   # Add a special group which includes all species

allTaxIds <- getAllTaxIds() # get list of all taxids in the current dataset

# Read profiles for all species
traits <- readAllProfiles( allTaxIds, profileSpecifications )
traits$tax_id <- rownames(traits)  # create an explicit tax_id column (rather than an index), by which we'll be able to merge (since left_join can't merge using indices)

taxidToKingdom$tax_id <- as.character(taxidToKingdom$tax_id)  # Convert the taxids column in the taxon memberships table to string, to match traits table

traits <- traits %>% left_join(taxidToKingdom, by="tax_id")  # Perform the merge
rownames(traits) <- traits$tax_id  # Restore the tax_ids index

traits$HighGC <- (traits$GenomicGC >  50) + 0
traits$LowGC  <- (traits$GenomicGC <= 50) + 0

traits$GC.45     <- as.factor(traits$GenomicGC > 45)
traits$Is_low_GC <- as.factor((traits$GenomicGC < 38) + 0)


traits$GC.0.40   <- (traits$GenomicGC <= 40) + 0
traits$GC.40.60  <- (traits$GenomicGC >= 40) & (traits$GenomicGC <= 60) + 0
traits$GC.60.100 <- (traits$GenomicGC >= 60) + 0

traits$TempHighLow75   <- as.factor((traits$OptimumTemp > 75) + 0)
traits$Is_high_temp    <- as.factor((traits$OptimumTemp > 58) + 0)
traits$TempOver75      <- (traits$OptimumTemp >  75)
traits$TempUnder75     <- (traits$OptimumTemp <= 75)

stopifnot(sum(traits[!is.na(traits$TempOver75),]$TempOver75   + 0) > 10)
stopifnot(sum(traits[!is.na(traits$TempUnder75),]$TempUnder75 + 0) > 100)

# Create a random trait
traits$Noise <- rnorm(nrow(traits))
# Create a random trait that only has values for approx. 30 of the species
make.partial.noise.trait <- function(nrows)
{
    out <- rnorm(nrows)
    out[rbinom(nrows,1,0.75)==1] <- NA
    return(out)
}

#traits$Noise30 <- rnorm(nrow(traits))
#traits[rbinom(nrow(traits),1,415/515)==1,"Noise30"] <- NA
#print(traits$Noise30)
traits$Noise30 <- make.partial.noise.trait( nrow(traits) )

stopifnot(nrow(traits) > 300 )  # sanity test
#stopifnot(ncol(traits) == profileLen+7 )


knownEndosymbionts = c(107806,203907,322098,331104,228908,1236703,1116230,218497,115713,353152,347515,5693,36329,99287,400667,83332,203267,1227812,272631,508771,224914,262768,272633,169963,227377,266834,264462,862908,283166,1321371,1208920,138677,331113,1165094)
# https://en.wikipedia.org/wiki/Intracellular_parasite#Obligate
# Buchnera, Blochmania, Aster yellows, Blattabacteria, Nanoarchaeum, Photodesmum, Piscirickettsia, Chlamydia, Cryptosporidium, Leishmania, Trypanosoma
traits$Is_endosymbiont <- factor(0, c(0,1))
for( taxId in knownEndosymbionts)
{
    traits[as.character(taxId),"Is_endosymbiont"] <- 1
}

protectFromNAs <- function(x) { replace(x, is.na(x), 0) }

#Is_endosymbiont GenomicGC GenomicENc.prime OptimumTemp Response precision    recall
#              1     37.99            56.51       58.01     0.14  0.659292 0.8186813
traits$Is_high_ENc_prime   <- as.factor( (traits$GenomicENc.prime > 56.5) + 0 )
print(traits$Is_low_GC)
print(sum(na.omit(traits$Is_low_GC == 1)))
print(sum(na.omit(traits$Is_high_ENc_prime == 1)))
print(sum(na.omit(traits$Is_endosymbiont==1)))
traits$Compound   <- as.factor((protectFromNAs(traits$Is_low_GC == 1) | protectFromNAs(traits$Is_high_ENc_prime==1) | protectFromNAs(traits$Is_endosymbiont==1) | protectFromNAs(traits$Is_high_temp==1)) + 0)
#print(traits$Compound)
#print(traits[order(traits$Profile_1.26),c("Profile_1.26", "Compound", "Is_endosymbiont", "OptimumTemp", "Is_high_temp", "GenomicENc.prime",  "Is_high_ENc_prime", "GenomicGC", "Is_low_GC")])


# ----------------------------------------------------------------------------------------------------------------------
# Scale-shape decomposition

profile.vars.1 <- getProfileVariables( c(1,31), profileId=1 )
profile.vars.2 <- getProfileVariables( c(2,32), profileId=2 )
#profile.vars.2 <- sapply(seq(32,2,-1), function (j) { sprintf("Profile_2.%d", j) } )
#dLFE.sd.1 <- apply( traits[,profile.vars.1], MARGIN=1, FUN=sd )
dLFE.sd.12 <- apply( traits[,c(profile.vars.1,profile.vars.2)], MARGIN=1, FUN=sd )
###dLFE.sd.12[20] <- dLFE.sd.12[20]*1.01 # TEST ONLY
###dLFE.sd.12[20] <- dLFE.sd.12[30]*0.99 # TEST ONLY

traits.normalized <- traits  # create another traits data-set, with profile values of each species normalized by sd 

#traits.normalized[,getProfileVariables( c(1,31), profileId=1 )] <- traits.normalized[,getProfileVariables( c(1,31), profileId=1 )] * (1/dLFE.sd.1)
#traits.normalized[,profile.vars] <- traits.normalized[,profile.vars] * (1/dLFE.sd)
traits.normalized[,profile.vars.1] <- traits.normalized[,profile.vars.1] * (1/dLFE.sd.12)
traits.normalized[,profile.vars.2] <- traits.normalized[,profile.vars.2] * (1/dLFE.sd.12)

traits.normalized$dLFE.sd.12 <- dLFE.sd.12  # save the sd values for each species as an extra trait (dLFE scale)
#traits.normalized
#dLFE.sd.12.check <- apply( traits.normalized[,c(profile.vars.1,profile.vars.2)], MARGIN=1, FUN=sd )

#########dLFE.normalized.sd <- apply( traits.normalized[,sapply(seq(32,2,-1), function (j) { sprintf("Profile_2.%d", j) } )], MARGIN=1, FUN=sd )
#traits.normalized[,getProfileVariables( c(1,31), profileId=1 )] <- traits.normalized[,getProfileVariables( c(1,31), profileId=1 )] * (1/dLFE.sd.1)
#dLFE.normalized.sd <- apply( traits.normalized[,profile.vars], MARGIN=1, FUN=sd )

dLFE.sd.12.check <- apply( traits.normalized[,c(profile.vars.1,profile.vars.2)], MARGIN=1, FUN=sd )
print(dLFE.sd.12.check)
#print(dLFE.normalized.sd)

print(all.equal( dLFE.sd.12.check, rep(1,0, length(dLFE.sd.12.check)), check.names=FALSE ))
stopifnot(all.equal( dLFE.sd.12.check, rep(1,0, length(dLFE.sd.12.check)), check.names=FALSE ))
#print(all.equal(dLFE.normalized.sd, rep(1.0, length(dLFE.normalized.sd)), check.names=FALSE  ))
#print(isTRUE(all.equal(dLFE.normalized.sd, rep(1.0, length(dLFE.normalized.sd)), check.names=FALSE )))
#stopifnot(isTRUE(all.equal(dLFE.normalized.sd, rep(1.0, length(dLFE.normalized.sd)), check.names=FALSE )))
#traits.normalized$Profile_1.sd <- dLFE.sd
################traits.normalized$Profile_2.sd <- dLFE.sd

remove(profile.vars.1)
remove(profile.vars.2)
# ----------------------------------------------------------------------------------------------------------------------


#----------------------------------


#print(traits["580340",])
#print(traits["511145",])


print("###################")
print(sum(is.na(traits$GenomicGC)))  # Print the number of missing items for these columns
print(sum(is.na(traits$Profile.15)))


#=========================================================================================
#===================== Read the tree and create the covariance matrix ====================
#=========================================================================================

tree <- read.tree(inputTree)

tree <- drop.tip( tree, speciesBlacklist )   # Filter species that will prevent analysis from being performed from the tree because of errors or missing data

N <- nTips(tree)
print(N)

#bmcorr <- corBrownian( phy=tree )

#speciesWithMissingData <- row.names(traits[is.na(traits[,"Profile.15"]),])
#print(speciesWithMissingData)
#tree.x <- drop.tip( tree, speciesWithMissingData )   # Filter species that will prevent analysis from being performed from the tree
#print(nTips(tree.x))
#print(nrow(traits[,"Profile.15"]))
#estimates.ou <- compar.ou( traits$Profile.15, tree.x, alpha=1.0 )
#print( estimates.ou$deviance)
#print( estimates.ou$para )
#
#quit()

summary(tree)

#print("------------------------------")
#print(nTips(tree))
#print(nrow(traits))

#treeTips <- tipLabels(traitsTree)
#treeSpecies <- tree$tip.label

#traits <- traits[treeSpecies,]   # Discard traits not found in the tree
#print(nrow(traits))

#-------------------------------------------------------------------------------------
# Draw VCV matrix

#vcv1 <- vcv( phy=tree, model="Brownian")
#dimnames(vcv1) <- list(1:ncol(vcv1))
#bmcorrmtx <- data.frame( melt( vcv1 ) )
#p <- ggplot( data=bmcorrmtx, aes(x=X1, y=X2, fill=value) ) + geom_raster() + scale_y_reverse(); p


#----------------------------------


performGLSregression <- function( traits, tree, Xtrait, Ytrait, plotRegression=TRUE, caption="", traits.full=NA, bm.gamma=1.0, extras="", plot.yrange=NA )
{
    # ==========================================================================================
    # ============================== Part 1 - Prapare Tree ======================================
    # ==========================================================================================
    # Discard tree tips for species missing data

    extraVars = c('Member_all_1')
    if( nchar(extras) )
    {
        extraVars <- all.vars(as.formula(paste(" ~ ", extras, sep="")))
    }
    
    speciesWithMissingData <- row.names(traits[(is.na(traits[Xtrait]) | is.na(traits[Ytrait]) | apply(is.na(traits[extraVars]), MARGIN=1, FUN=any) ),])
    tree <- drop.tip( tree, speciesWithMissingData )   # Filter species that will prevent analysis from being performed from the tree
    treeSpecies <- tree$tip.label

    # Create separate data-frame containing species excluded only because they are not found in the tree (this is for display only and is not related to the regression)
    if( isTRUE(is.na( traits.full )) )
    {
        traits.full <- duplicate( traits )
    }
    traits.total.count <- nrow(traits.full)
    traits.notintree <- traits.full[setdiff(row.names(traits.full), treeSpecies),]  # create a traits matrix of species not included in the tree (in order to plot them separately so they provide additional visual confirmation of the results)
    traits.notintree <- traits.notintree[(!is.na(traits.notintree[Xtrait]) | !is.na(traits.notintree[Ytrait])),]  # the "not-in-tree" dataset should only inlcude species that can be plotted...
    
    # Discard trait data for species missing from the tree
    traits <- traits[treeSpecies,]   # Discard traits not found in the tree

    print(sprintf("%d included  (%d excluded, %d total)", nrow(traits), nrow(traits.notintree), traits.total.count ) )
    #stopifnot( nrow(traits) + nrow(traits.notintree) == traits.total.count )  # incorrect
    
    if( is.null(tree) )  # no species remaining in the tree
    {
        return( c(NA,NA,NA))
    }
    if( all( traits[,Xtrait]==traits[1,Xtrait] ) ) # the regressor has only a single value left (i.e., a factor with one level represented)
    {
        return( c(NA,NA,NA))
    }
    
    # Check tree <-> data-frame correlation appears valid
    x <- phylo4d(tree, tip.data=traits, rownamesAsLabels=TRUE)
    tipsOk <- hasTipData( x )
    #print(tipsOk)
    stopifnot(tipsOk)

    N <- nTips(tree)

    print(sprintf("%s -> %s (N=%d)", Xtrait, Ytrait, N))

    # Calculate the correlation matrix (for a BM process)
    ###print(bm.gamma)
    #bmcorr <- corBrownian( value=bm.gamma, phy=tree )
    bmcorr <- corBrownian( phy=tree )
    
    #estimates.ou <- compar.ou( traits[,Ytrait], tree, alpha=0.5 )
    #print( estimates.ou$deviance)
    #print( estimates.ou$para )
    oucorr <- corMartins( 1.0, tree )
    

    # ==========================================================================================
    # ============================ Part 2 - Perform Regressions ================================
    # ==========================================================================================

    # For discrete variables, we are interested in whether the effect trait is significantly different in any of the groups. Consequently,
    # we set the mean to be 0.
    if( any(class(traits[,Xtrait])==c("factor")) )
    {
        ## print(Ytrait)
        ## # TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY #
        ## #if( Ytrait=="Profile.14" )
        ## #{
        ##     levelToAlter <- levels(traits[,Xtrait])[2]
        ##     print(levelToAlter)
        ##     print(paste("Warning: altering", levelToAlter))
        ##     print(traits[(traits[,Xtrait]==levelToAlter), Ytrait])
        ##     traits[(traits[,Xtrait]==levelToAlter), Ytrait] <- traits[(traits[,Xtrait]==levelToAlter), Ytrait] - 0.2
        ##     print(traits[(traits[,Xtrait]==levelToAlter), Ytrait])
        ## #}
        ## # TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY #

        normYtrait <- paste0(Ytrait, ".norm")
        traits[,normYtrait] <- traits[,Ytrait] - mean(traits[,Ytrait])
        Ytrait <- normYtrait   # proceed in analysis using the normalized trait

    }

    if( rankTransformResponse )
    {
        ###print(Ytrait)
        rankedYtrait <- paste(Ytrait, ".ranks", sep="" )
        traits[,rankedYtrait] <- rank(traits[,Ytrait], ties.method="average")
        #print(traits[,rankedYtrait])
        Ytrait <- rankedYtrait
        ####print(traits[,Ytrait])
    }
    
    
    # Set the predictive (regression) and intercept-only formulas
    # 'predictiveFormula' is the actual regression. It uses an intercept when a continuous explanatory var is used (and none with a discrete var, as in one-way ANOVA)
    # 'inteceptFormula' is an intecept-only model, used as the baseline for calculation of R^2
    predictiveFormula <- NA    
    if( any(class(traits[,Xtrait])==c("factor")) )
    {
        predictiveFormula <- as.formula( paste(Ytrait, " ~ ", Xtrait, " + 0") )  # model without intercept term
    }
    else
    {
        #predictiveFormula <- as.formula( paste(Ytrait, " ~ ", Xtrait) )          # model with intercept term
        if( nchar(extras) )
        {
            predictiveFormula <- as.formula( paste(Ytrait, " ~ ", Xtrait, " + ", extras ) )          # model with intercept term
        }
        else
        {
            predictiveFormula <- as.formula( paste(Ytrait, " ~ ", Xtrait ) )                         # model with intercept term
        }
        
    }
    #interceptFormula <- as.formula( paste(Ytrait, " ~ ", "1") )
    interceptFormula <- as.formula( paste(Ytrait, " ~ ", " 1 ") )

    ###print( all.vars(predictiveFormula) )
    
    
    # Perform OLS regression (used as a reference for comparison)
    m1 <- lm(  predictiveFormula, traits); summary(m1)
    co <- coef(m1)

    # Perform GLS regression using the predictive and intercept-only models
    gls1 <- gls( predictiveFormula,
                traits,
                correlation=bmcorr,
                na.action=na.omit,
                method="REML"); #print( summary(gls1))
    co2 <- coef(gls1)
    #print(summary(gls1))
    #print(anova(gls1))

    gls0 <- gls( interceptFormula,
                 traits,
                 correlation=bmcorr,
                 na.action=na.omit,
                 method="REML")
    #print(bmcorr)

    gls.nocorr <- gls( predictiveFormula,
                      traits,
                      na.action=na.omit,
                      method="REML");

    gls1.ou <- NA
    try(gls1.ou <- gls( predictiveFormula,
                   traits,
                   correlation=oucorr,
                   na.action=na.omit,
                   method="REML") );
    ###print(summary(gls1.ou))
    
    #
    # Buse's R^2
    # Reference: A. Buse, "Goodness of Fit in Generalized Least Squares Estimation", The American Statistician, June 1973, Vol. 27, No. 3.
    # URL: http://www.jstor.org/stable/2683631
    #
    # Notation: u.hat := epsilon (residuals)
    #           V := Omega (Covariance matrix)
    #           Y.bar := Intercept of equivalent intercept-only model (see: http://r.789695.n4.nabble.com/Pseudo-R-squared-in-gls-model-td4641148.html)
    #           e := first column in design matrix X (i.e., rep(1, N) if intercepts are used). See above eq. (11) p. 107.
    #
    u.hat <- residuals(gls1)
    #print(u.hat)
    #V <- getVarCov( gls1 )   # this doesn't work (bug?)
    V <- corMatrix( gls1$modelStruct$corStruct )
    inv.V <- solve(V)   # invert V (=Omega). Every positive-definite matrix is invertible (https://en.wikipedia.org/wiki/Positive-definite_matrix#cite_note-5)
    if( any(class(traits[,Xtrait])==c("factor")) )
    {
        e <- rep(0, N)
    }
    else
    {
        e <- rep(1, N)
    }
    Y <- traits[,Ytrait]
    Y.bar <- coef(gls0)["(Intercept)"]
    yye <- Y - Y.bar * e
    
    R2 <- 1 - ( t(u.hat) %*% inv.V %*% u.hat )/( t(yye) %*% inv.V %*% yye )   # See eq. (15) p. 107.
    stopifnot( R2 >= 0 && R2 <= 1.0 )

    # Slope direction
    slopeDirection <- NA
    if( any(class(traits[,Xtrait])==c("factor")) )
    {
        slopeDirection <- 1  # No slope if the explanatory variable is discrete
    }
    else
    {
        slopeDirection <- sign(coef(gls1)[Xtrait])
    }
    stopifnot( slopeDirection==1 || slopeDirection==-1 )

    # Extract p-values to represent the GLS and OLS models
    pvalue <- NA
    pvalue.OLS <- NA
    if( any(class(traits[,Xtrait])==c("factor")) )
    {
        # For discrete explanatory variables, we perform an "omnibus" F-test, i.e., we determine how much of the variance is explained by the division into the specified levels (the null hypothesis being that the values for all levels are drawn from identical populations)
        pvalue     <- anova(gls1)[1,'p-value']
        pvalue.OLS <- anova(m1)[1,'Pr(>F)']
    }
    else
    {
        # For continuous explanatory variables, we perform a t-test to determine if the coefficient (slope) is different than 0.
        pvalue     <- summary(gls1)$tTable[2,4]
        pvalue.OLS <- summary(m1)$coefficients[2,4]
    }

    ## print("////////////////////////////")
    ## print(sprintf("AIC(m1):         %.3g", AIC(m1)))
    ## print(sprintf("AIC(gls1):       %.3g", AIC(gls1)))
    ## print(sprintf("AIC(gls0):       %.3g", AIC(gls0)))
    ## print(sprintf("AIC(gls.nocorr): %.3g", AIC(gls.nocorr)))
    ## if( !is.na(gls1.ou) )
    ## {
    ##     print(sprintf("AIC(gls1.ou):    %.3g", AIC(gls1.ou)))
    ## }
    ## else
    ## {
    ##     print(sprintf("AIC(gls1.ou):    N/A"))
    ## }
    ## print("////////////////////////////")

    # direction indicator (this is used by glsRegressionRangeAnalysis() to create the third plot in each series, which shows the effect of the different levels)
    directionsIndicator <- NA
    if( any(class(traits[,Xtrait])==c("factor")) )
    {
        # for discrete explanatory variables, the indicator is a string of '-', '0' or '+' characters
        signs <- sign(summary(gls1)$tTable[,1])  # get the sign of the coefficients for all levels
        significant <- (summary(gls1)$tTable[,4] < significanceLevel) + 0
        
        directionsIndicator <- paste( c("-", "0", "+")[((signs * significant)+2)], collapse="")
    }
    else
    {
        # for continous explanatory vars, the indicator is a single character
        if( pvalue > significanceLevel )
        {
            directionsIndicator <- "0"
        }
        else
        {
            signs <- sign(coef(gls1)[Xtrait])
            directionsIndicator <- c("-", "0", "+")[(signs+2)]
        }
    }

    if( !any(class(traits[,Xtrait])==c("factor")) )
    {
        # Calculate MINE statistics (including MIC)
        mine.results <- mine( traits[,Xtrait], traits[,Ytrait])
        # For comparison, calculate Pearson-r
        mine.results$Pearson.r <- cor( traits[,Xtrait], traits[,Ytrait], use="complete.obs", method="pearson" )

        # Estimate p-value for MIC using permutations
        # Ref (not definitive): https://stats.stackexchange.com/a/199681
        M <- length( traits[,Xtrait] )
        num.perms <- MIC.pval.num.iterations
        
        count.more.extreme <- 0
        for( i in 1:num.perms )
        {
            y.perm <- sample( traits[,Ytrait], size=M, replace=TRUE)
            mine.perm <- mine( traits[,Xtrait], y.perm )

            if( mine.perm$MIC > mine.results$MIC )
            {
                count.more.extreme <- count.more.extreme+1
            }
        }
        mine.results$MIC.pval <- count.more.extreme / num.perms
        print(mine.results$MIC.pval)

        # Perform Breush-Pagan homoskedasticity test
        bp.result <- bptest( m1 )
    }
    else
    {
        # TODO implement MIC for categorical vars?
        mine.results <- data.frame( MIC=c(0.0), MAS=c(0.0), MIC.pval=c(1.0), Pearson.r=c(0.0) )
        bp.result <- data.frame( p.value=c(1.0) )
    }
    
    # ==========================================================================================
    # ========================= Part 3 (optional) - Plot Regression ============================
    # ==========================================================================================
    if( plotRegression )
    {
        clr2 <- "#55CCB0"

        max.y  <- max(0, max(traits[,Ytrait])) # y-axis is always included
        min.y  <- min(0, min(traits[,Ytrait])) # y-axis is always included
        line.y <- (max.y-min.y)/32  # line height for diagnostic text

        # Note: this uses plotmath syntax; See: https://www.rdocumentation.org/packages/grDevices/versions/3.4.3/topics/plotmath
    
        if( any(class(traits[,Xtrait])==c("factor")) )
        {
            factor.pvals <- summary(gls1)$tTable[,4]
            vars <- names(factor.pvals)
        
            significanceMarkers <- data.frame(x=1:length(vars),
                                              y=rep(max.y+(max.y-min.y)*0.02, length(vars)),
                                              var=vars,
                                              pval=factor.pvals,
                                              pval.char=vapply(factor.pvals, function(x) sprintf("%.3g", x), character(1)),
                                              label=vapply(factor.pvals, function(pval) { if ( pval < significanceLevel) "*" else ""}, character(1)) )
            
            p <- ggplot(traits, aes(get(Xtrait), get(Ytrait))) +
              labs(y=Ytrait, x=Xtrait, title=caption) +
              geom_boxplot(outlier.size=0, outlier.alpha=0) +  # don't show outliers on boxplot, since we're plotting all data anyway...
              #geom_point( data=traits.notintree, aes_string(x=Xtrait, y=Ytrait), colour="#30f050", alpha=0.4 ) +
              geom_jitter(colour="black") +
              geom_text( data=significanceMarkers, aes(x=x, y=y, label=var),       colour="gray",   size=2, alpha=0.7 ) +
              geom_text( data=significanceMarkers, aes(x=x, y=y, label=pval.char), colour="black",  size=2  ) +
              geom_text( data=significanceMarkers, aes(x=x, y=y, label=label),     colour="black",  size=10 ) +
              geom_hline( yintercept = 0 ) +
              annotate( "text", x=Inf,  y=max.y+c(0, 1, 2, 3) * line.y, label=c(  "lm",  sprintf("italic(R) ^ 2 == %.3g", summary(m1)$r.squared), sprintf('italic(p)*"-val" ==  %.3g', pvalue.OLS), sprintf("italic(N) == %d", nrow(traits))), hjust=1, vjust=0, colour="red", parse=TRUE ) +
                annotate( "text", x=-Inf, y=max.y+c(0, 1, 2, 3) * line.y, label=c(  "gls", sprintf("italic(R) ^ 2 == %.3g", R2),                    sprintf('italic(p)*"-val" ==  %.3g', pvalue), sprintf("italic(N) == %d", nrow(traits))),     hjust=0, vjust=0, colour=clr2,  parse=TRUE )

            if( nrow(traits) < max.num.species.to.label )  # add species names (if there aren't too many points that the graph becomes cluttered)
            {
                p <- p + geom_text( data=traits, aes( x=get(Xtrait), y=get(Ytrait), label=short.name), colour="#404040", size=4, hjust=0, nudge_x=0.05, alpha=0.7 )
            }
            
            print(p)
        }
        else
        {
            p <- ggplot(traits, aes(get(Xtrait), get(Ytrait))) +
                labs(y=Ytrait, x=Xtrait, title=caption) +
                geom_point( data=traits.notintree, aes_string(x=Xtrait, y=Ytrait),  colour="#be7964", alpha=0.8 ) +  #30f050
                geom_abline( aes( slope=co[Xtrait],  intercept=co["(Intercept)"]),  colour="#be7964", size=1.4, alpha=0.8   ) +
                geom_abline( aes( slope=co2[Xtrait], intercept=co2["(Intercept)"]), colour=clr2,      size=1.4, alpha=0.8   ) +
                geom_point() +
                geom_hline( yintercept = 0 ) +
                annotate( "text", x=Inf,  y=max.y+c(3, 2, 1, 0) * line.y, label=c(  "lm",  sprintf("italic(R) ^ 2 == %.3g", summary(m1)$r.squared), sprintf('italic(p)*"-val" == %.3g', pvalue.OLS), sprintf("italic(N) == %d", nrow(traits))), hjust=1, vjust=0, colour="#be7964", parse=TRUE ) +
                annotate( "text", x=-Inf, y=max.y+c(3, 2, 1, 0) * line.y, label=c(  "gls", sprintf("italic(R) ^ 2 == %.3g", R2),                    sprintf('italic(p)*"-val" == %.3g', pvalue),     sprintf("italic(N) == %d", nrow(traits))), hjust=0, vjust=0, colour=clr2,  parse=TRUE ) +
                annotate( "text", x=Inf, y=max.y+c(7, 6, 5, 4) * line.y, label=c(  "MINE", sprintf("MIC == %.3g", mine.results$MIC),                    sprintf('italic(p)*"-val" == %.3g', mine.results$MIC.pval),    sprintf("Pearson-r == %.3g", mine.results$Pearson.r )), hjust=1, vjust=0, colour='#4444cc',  parse=TRUE ) +
                theme( plot.background = element_blank(),   # Hide unnecessary theme elements (background panels, etc.)
                      panel.grid.major.y = element_line(color="grey", size=0.70),
                      panel.grid.major.x = element_line(color="grey", size=0.70),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank()
                      ) # +
                    
            
            if( nrow(traits) < max.num.species.to.label )  # add species names (if there aren't too many points that the graph becomes cluttered)
            {
                p <- p + geom_text( data=traits,           aes( x=get(Xtrait), y=get(Ytrait), label=short.name), colour="#404040", size=3, hjust=0, nudge_x=0.05, alpha=0.7 ) +
                         geom_text( data=traits.notintree, aes( x=get(Xtrait), y=get(Ytrait), label=short.name), colour="#30f050", size=3, hjust=0, nudge_x=0.05, alpha=0.7 )
            }

            if( !isTRUE(is.na(plot.yrange)))
            {
                p <- p + scale_y_continuous( limits=plot.yrange )
            }
            
            
            print(p)
        }
    }

    ## print("???????")
    ## print(mine.results)
    ## print(bp.result)
    ## print(R2)
    ## print(slopeDirection)
    ## print(N)
    ## print(directionsIndicator)
    return( data.frame( pvalue=pvalue, R2=R2 * slopeDirection, N=N, directionsIndicator=directionsIndicator, MIC=mine.results$MIC, MAS=mine.results$MAS, pearson.r=mine.results$Pearson.r, MIC.pvalue=mine.results$MIC.pval, BP.pvalue=bp.result$p.value ) )
}

performOLSregression <- function( traits, Xtrait, Ytrait, plotRegression=TRUE, extras="", caption="", traits.all=NA, colorTrait="Member_all_1" )
{
    #print("******** performOLSregression *********")
    # For discrete variables, we are interested in whether the effect trait is significantly different in any of the groups. Consequently,
    # we set the mean to be 0.
    if( any(class(traits[,Xtrait])==c("factor")) )
    {
        normYtrait <- paste0(Ytrait, ".norm")
        traits[,normYtrait] <- traits[,Ytrait] - mean(traits[,Ytrait], na.rm=TRUE)
        Ytrait <- normYtrait   # proceed in analysis using the normalized trait
    }
    #print(Ytrait)
    
    
    predictiveFormula <- NA    
    if( any(class(traits[,Xtrait])==c("factor")) )
    {
        predictiveFormula <- as.formula( paste(Ytrait, " ~ ", Xtrait, " + 0") )  # model without intercept term
    }
    else
    {
        #predictiveFormula <- as.formula( paste(Ytrait, " ~ ", Xtrait) )          # model with intercept term
        if( nchar(extras) )
        {
            predictiveFormula <- as.formula( paste(Ytrait, " ~ ", Xtrait, " + ", extras ) )          # model with intercept term
        }
        else
        {
            predictiveFormula <- as.formula( paste(Ytrait, " ~ ", Xtrait ) )                         # model with intercept term
        }
        
    }
    #interceptFormula <- as.formula( paste(Ytrait, " ~ ", "1") )
    interceptFormula <- as.formula( paste(Ytrait, " ~ ", " 1 ") )

    extraVars = c('Member_all_1')
    if( nchar(extras) )
    {
        print("##################################################################")
        extraVars <- all.vars(as.formula(paste(" ~ ", extras, sep="")))
    }

    speciesWithMissingData <- row.names(traits[(!is.na(traits[Xtrait]) & !is.na(traits[Ytrait]) & apply(!is.na(traits[extraVars]), MARGIN=1, FUN=all) ),])
    #print(speciesWithMissingData)
    traits <- traits[speciesWithMissingData,]   # Discard species with missing values for one of the traits
    N <- nrow( traits )

    if( isTRUE( is.na(traits.all) ) ) { traits.all <- traits[0,]; stopifnot(nrow(traits.all)==0) }
    
    if( N < 3 )
    {
        print("******** performOLSregression - insufficient data *********")
        print(N)
        print(traits)
        return( NA );
    }
    
    #print("******** performOLSregression - regression *********")
    
    # Perform OLS regression
    m1 <- lm(  predictiveFormula, traits); summary(m1)
    co <- coef(m1)

    # Slope direction
    slopeDirection <- NA
    if( any(class(traits[,Xtrait])==c("factor")) )
    {
        slopeDirection <- 1  # No slope if the explanatory variable is discrete
    }
    else
    {
        slopeDirection <- sign(coef(m1)[Xtrait])
    }
    stopifnot( slopeDirection==1 || slopeDirection==-1 )

    # Extract p-values to represent the GLS and OLS models
    pvalue.OLS <- NA
    if( any(class(traits[,Xtrait])==c("factor")) )
    {
        # For discrete explanatory variables, we perform an "omnibus" F-test, i.e., we determine how much of the variance is explained by the division into the specified levels (the null hypothesis being that the values for all levels are drawn from identical populations)
        pvalue.OLS <- anova(m1)[1,'Pr(>F)']
    }
    else
    {
        # For continuous explanatory variables, we perform a t-test to determine if the coefficient (slope) is different than 0.
        pvalue.OLS <- summary(m1)$coefficients[2,4]
    }

    mine.results <- list( MIC=NA, MAS=NA, MIC.pval=NA, Pearson.r=NA )

    if( !any(class(traits[,Xtrait])==c("factor")) )
    {
        print("Calculating MIC...")
        # Calculate MINE statistics (including MIC)
        mine.results <- mine( traits[,Xtrait], traits[,Ytrait])
        # For comparison, calculate Pearson-r
        mine.results$Pearson.r <- cor( traits[,Xtrait], traits[,Ytrait], use="complete.obs", method="pearson" )

        # Estimate p-value for MIC using permutations
        # Ref (not definitive): https://stats.stackexchange.com/a/199681
        M <- length( traits[,Xtrait] )
        num.perms <- MIC.pval.num.iterations
        
        count.more.extreme <- 0
        for( i in 1:num.perms )
        {
            y.perm <- sample( traits[,Ytrait], size=M, replace=TRUE)
            mine.perm <- mine( traits[,Xtrait], y.perm )

            if( mine.perm$MIC > mine.results$MIC )
            {
                count.more.extreme <- count.more.extreme+1
            }
        }
        mine.results$MIC.pval <- count.more.extreme / num.perms
        print(mine.results$MIC.pval)

        # Perform Breush-Pagan homoskedasticity test
        bp.result <- bptest( m1 )
    }
    else
    {
        # TODO IMPL MIC for factors
        bp.result <- list("p.value"=NA)
    }
    
    #print(mine.results)
    #print("******** performOLSregression - before plotting *********")

    # ==========================================================================================
    # ========================= Part 3 (optional) - Plot Regression ============================
    # ==========================================================================================
    if( plotRegression )
    {
        clr2 <- "#55CCB0"

        max.y  <- max(0, max(traits[,Ytrait])) # y-axis is always included
        min.y  <- min(0, min(traits[,Ytrait])) # y-axis is always included
        line.y <- (max.y-min.y)/32  # line height for diagnostic text

        # Note: this uses plotmath syntax; See: https://www.rdocumentation.org/packages/grDevices/versions/3.4.3/topics/plotmath
    
        if( any(class(traits[,Xtrait])==c("factor")) )
        {
            ## factor.pvals <- summary(gls1)$tTable[,4]
            ## vars <- names(factor.pvals)
        
            ## significanceMarkers <- data.frame(x=1:length(vars),
            ##                                   y=rep(max.y+(max.y-min.y)*0.02, length(vars)),
            ##                                   var=vars,
            ##                                   pval=factor.pvals,
            ##                                   pval.char=vapply(factor.pvals, function(x) sprintf("%.3g", x), character(1)),
            ##                                   label=vapply(factor.pvals, function(pval) { if ( pval < significanceLevel) "*" else ""}, character(1))
            ##                                  )
            
            p <- ggplot(traits, aes(get(Xtrait), get(Ytrait))) +
              labs(y=Ytrait, x=Xtrait, title=caption) +
              geom_boxplot(outlier.size=0, outlier.alpha=0) +  # don't show outliers on boxplot, since we're plotting all data anyway...
              #geom_point( data=traits.notintree, aes_string(x=Xtrait, y=Ytrait), colour="#30f050", alpha=0.4 ) +
              geom_jitter(colour="black") +
              #geom_text( data=significanceMarkers, aes(x=x, y=y, label=var),       colour="gray",   size=2, alpha=0.7 ) +
              #geom_text( data=significanceMarkers, aes(x=x, y=y, label=pval.char), colour="black",  size=2  ) +
              #geom_text( data=significanceMarkers, aes(x=x, y=y, label=label),     colour="black",  size=10 ) +
              geom_hline( yintercept = 0 ) +
                annotate( "text", x=Inf,  y=max.y+c(0, 1, 2, 3) * line.y, label=c(  "lm",  sprintf("italic(R) ^ 2 == %.3g", summary(m1)$r.squared), sprintf('italic(p)*"-val" ==  %.3g', pvalue.OLS), sprintf("italic(N) == %d", nrow(traits))), hjust=1, vjust=0, colour="red", parse=TRUE ) +

            if( nrow(traits) < max.num.species.to.label )  # add species names (if there aren't too many points that the graph becomes cluttered)
            {
                p <- p + geom_text( data=traits, aes( x=get(Xtrait), y=get(Ytrait), label=short.name), colour="#404040", size=4, hjust=0, nudge_x=0.05, alpha=0.7 )
            }
            
            print(p)
        }
        else
        {
            p <- ggplot(traits, aes(get(Xtrait), get(Ytrait))) +
                labs(y=Ytrait, x=Xtrait, title=caption) +
                geom_point( data=traits.all, aes_string(x=Xtrait, y=Ytrait), colour="#be7964", alpha=0.7 ) +
                geom_point( aes(color=get(colorTrait) ) ) +
                geom_hline( yintercept = 0 ) +
                geom_abline( aes( slope=co[Xtrait],  intercept=co["(Intercept)"]),  colour="#be7964"  ) +
                geom_abline( aes( slope=co2[Xtrait], intercept=co2["(Intercept)"]), colour=clr2   ) +
                annotate( "text", x=Inf,  y=max.y+c(3, 2, 1, 0) * line.y, label=c(  "lm",  sprintf("italic(R) ^ 2 == %.3g", summary(m1)$r.squared), sprintf('italic(p)*"-val" == %.3g', pvalue.OLS), sprintf("italic(N) == %d", nrow(traits))), hjust=1, vjust=0, colour="#be7964", parse=TRUE ) +
                annotate( "text", x=Inf, y=max.y+c(7, 6, 5, 4) * line.y, label=c(  "MINE", sprintf("MIC == %.3g", mine.results$MIC),                    sprintf('italic(p)*"-val" == %.3g', mine.results$MIC.pval),    sprintf("Pearson-r == %.3g", mine.results$Pearson.r )), hjust=1, vjust=0, colour='#4444cc',  parse=TRUE ) +
                #scale_colour_manual( values=c("0"="black", "1"="#5658eb") ) +
#                scale_colour_hue( ) +
                guides(color=FALSE)

            
            if( nrow(traits) < max.num.species.to.label )  # add species names (if there aren't too many points that the graph becomes cluttered)
            {
                p <- p +
                    geom_text( data=traits.all, aes( x=get(Xtrait), y=get(Ytrait), label=short.name), colour="#30f050", size=3, hjust=0, nudge_x=0.05, alpha=0.7 ) +
                    geom_text( data=traits,     aes( x=get(Xtrait), y=get(Ytrait), label=short.name), colour="#404040", size=3, hjust=0, nudge_x=0.05, alpha=0.7 ) 
            }
            
            print(p)
        }
    }

    return( data.frame( pvalue=pvalue.OLS, R2=summary(m1)$r.squared * slopeDirection, N=N, directionsIndicator=NA, MIC=mine.results$MIC, MAS=mine.results$MAS, pearson.r=mine.results$Pearson.r, MIC.pvalue=mine.results$MIC.pval, BP.pvalue=bp.result$p.value ) )
}



performGLSregression_profileRangeMean <- function( traits, tree, Xtrait, profileRange, profileId=1, plotRegression=FALSE, caption="", extras="", traits.full=NA, plot.yrange=NA )
{
    #print("performGLSregression_profileRangeMean(): -->")
    #print(profileRange)
    #print(class(profileRange))
    #print("profileRange:")
    #print(profileRange)
    #print(class(profileRange))
    stopifnot( class(profileRange)=="integer" || class(profileRange)=="numeric" )
    stopifnot(length(profileRange)==2)
    # Make list of the selected profile columns
    variables <- getProfileVariables( profileRange, profileId )
    #print(variables)

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

    # Extract the selected columns and compute the mean values
    # TODO - learn dplyr...
    traits$RangeMean <- rowMeans( traits[variables] )

    if( !isTRUE( is.na(traits.full) ) )
    {
        traits.full$RangeMean <- rowMeans( traits.full[variables] )
    }
    

    # Perform the regression (for the mean values)
    return( performGLSregression( traits, tree, Xtrait, "RangeMean", plotRegression=plotRegression, caption=caption, extras=extras, traits.full=traits.full, plot.yrange=plot.yrange ) )
}

performOLSregression_profileRangeMean <- function( traits, Xtrait, profileRange, profileId=1, plotRegression=FALSE, caption="", extras="", colorTrait="Member_all_1" )
{
    #print("******** performOLSregression_profileRangeMean *********")
    #print("performGLSregression_profileRangeMean(): -->")
    #print(profileRange)
    #print(class(profileRange))
    #print("profileRange:")
    #print(profileRange)
    #print(class(profileRange))
    stopifnot( class(profileRange)=="integer" || class(profileRange)=="numeric" )
    stopifnot(length(profileRange)==2)
    # Make list of the selected profile columns
    variables <- getProfileVariables( profileRange, profileId )
    #print(variables)

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

    # Extract the selected columns and compute the mean values
    # TODO - learn dplyr...
    traits$RangeMean <- rowMeans( traits[variables] )

    # Perform the regression (for the mean values)
    return( performOLSregression( traits, Xtrait, "RangeMean", plotRegression=plotRegression, caption=caption, extras=extras, colorTrait=colorTrait ) )
}

# Extract legend from ggplot figure
# Source: https://stackoverflow.com/a/12041779
# Tested with R 3.4.1
getLegend <- function(a.gplot)
{ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
}

# Tested with R 3.4.1
plotRotatedPyramid <- function( plotdf, dotsDF, valLabel, valScale, Xtrait, plotTitle, profileVisualization=NA )
{
    # Create the main pyramid plot
    if( class(plotdf$value)=="factor")
    {
        pyramidPlot <- ggplot( data=plotdf, aes(x=Var1, y=Var2, fill=factor(value) ) )

    }
    else
    {
        pyramidPlot <- ggplot( data=plotdf, aes(x=Var1, y=Var2, fill=       value  ) )
    }
    
    plot.elements <- list( 
        geom_tile( linetype=0, colour=NA ) )   # Use geom_tile (because geom_raster has poor pdf/svg output), with minimal borders

    if( class(plotdf$value)=="factor")
    {
        plot.elements <- c(plot.elements, 
                   geom_text( data=plotdf, aes( x=Var1, y=Var2, label=factor(value) ), size=1, colour="white", angle=45 )
                   )
    }

    plot.elements <- c(plot.elements, list(
         labs( y="Profile End", x="Profile Start", title=plotTitle ) ,
         valScale ,
         coord_fixed() ,
         guides( fill=FALSE ) ,     # Hide the legend (will be plotted separately)
         geom_point( data=dotsDF, aes(x=Var1,y=Var2), position=position_nudge(y=4, x=-4), color="white", size=((fontScale/12)**2)*0.6 ),  # Significance markers (white dots)
         theme( plot.margin = unit(c(0.1,0.1,0.1,0.1), "npc"), # Set the margins (should match the calculation used to create the "diagonal" X axis
               plot.background = element_blank(),   # Hide unnecessary theme elements (background panels, etc.)
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               title       = element_text( size=round(fontScale) ),
               axis.text.y = element_text( size=round(fontScale) ),
               axis.text.x = element_blank(),
               axis.line = element_line(color="black"),
               aspect.ratio=1 )  # Fix 1:1 aspect ratio
         ) )
    
    #----------------------------------------------------------------------
    # Draw heat-map pyramid on canvas
    
    grid.newpage()
    pushViewport( viewport( x=unit(0.5, "npc"), y=unit(0.35, "npc"), angle=-45 ) )

    pyramidPlotGrob <- ggplotGrob( pyramidPlot + plot.elements )
    grid.draw( pyramidPlotGrob )

    
    uu1a <- convertX( unit(1, "native"), "mm" )
    #print(uu1a)
    uu1b <- convertY( unit(1, "native"), "mm" )
    #print(uu1b)
        
    uw <- 0.12  # triangle offset
    grid.polygon( x=c(0+uw, 1.0, 1.0, 0+uw), y=c(0.0, 1-uw, 0.0, 0.0), default.units="native", gp=gpar(fill="white", col=NA ) )

    
    #----------------------------------------------------------------------
    # Draw X-axis

    ww2 <- convertUnit(uu1b * 0.83 * sqrt(2.0), "npc")  # Todo - this returns more than 1 npcs (screen width). Why?
    #print(ww2)
    
    #u2 <- 0.05
    u2 <- 0.05
    #vp2 <- viewport( x=unit(0.5+u2, "native"), y=unit(0.5-u2, "native"), h=unit(0.01, "native"), w=ww2, angle=45 )
    vp2 <- viewport( x=unit(0.5+u2, "native"), y=unit(0.5-u2, "native"), h=unit(0.01, "native"), w=unit(0.9, "native"), angle=45 )
    pushViewport( vp2 )
    
    #grid.lines( c(0.0, 1.0), c(0.5, 0.5) )

    xaxis <- xaxisGrob( at=seq(0,1,(cdsScaleResolution/profileStep)/(pyramidLength-1)) )

    xaxis <- editGrob(xaxis,
                      gPath("major"),
                      x=unit(c(0,1), "npc")
                      )             
    xaxis <- editGrob(xaxis,        
                      gPath("labels"),
                      label=seq( 0, (pyramidLength-1) * profileStep, cdsScaleResolution)
                      )
    xaxis <- editGrob(xaxis,        
                      gPath("labels"),
                      gp=gpar( fontsize=round(fontScale*2) )
                      )
    grid.draw(xaxis)

    #----------------------------------------------------------------------
    # Draw legend
    upViewport()
    
    vp3 <- viewport( x=unit(0.0, "native"), y=unit(0.5, "native"), h=unit(0.5, "native"), w=unit(0.5, "native"), angle=45 )
    pushViewport( vp3 )

    if( class(plotdf$value)=="factor")
    {
        p4 <- ggplot( data=plotdf, aes(x=Var1, y=Var2, fill=factor(value) ) )
    }
    else
    {
        p4 <- ggplot( data=plotdf, aes(x=Var1, y=Var2, fill=       value  ) )
    }

    p4 <- p4 + 
         valScale +
    #         scale_fill_gradient2(low="#cc7777", mid="#ffff88", high="#000000", limits=c(min(subranges$Logpval), 0), midpoint=minLogpval ) +
         labs( value=valLabel ) +
         geom_tile()
    p4leg <- getLegend( p4 )
    grid.draw( p4leg )

    #----------------------------------------------------------------------
    # Annotate stats
    upViewport(2)
    #vp4 <- viewport( x=unit(0.5, "native"), y=unit(1.0, "native"), h=unit(0.5, "native"), w=unit(0.5, "native"), angle=45 )
    vp4 <- viewport( x=unit(0.95, "npc"), y=unit(0.5, "npc"), h=unit(0.5, "npc"), w=unit(0.5, "npc"), angle=0 )
    pushViewport( vp4 )

    if( class(plotdf$value)=="factor" )
    {
        grid.text(sprintf("levels==%d", nlevels(plotdf$value)), x=unit(0.5, "npc"), y=unit(0.5, "npc"), just="right" )
    }
    else
    {
        grid.text(sprintf("min==%.3g", min(plotdf[,"value"])), x=unit(0.5, "npc"),                     y=unit(0.7, "npc"),                      just="right", gp=gpar( fontsize=round(fontScale) ) )
        grid.text(sprintf("max==%.3g", max(plotdf[,"value"])), x=unit(0.5, "npc") - unit(0, "points"), y=unit(0.7, "npc") + unit(15, "points"), just="right", gp=gpar( fontsize=round(fontScale) ) )
    }

    #----------------------------------------------------------------------
    # Draw profile visualization
    if( !is.na(profileVisualization) )
    {
        upViewport()

        vp5 <- viewport( x=unit(0.45, "npc"), y=unit(0.1, "npc"), h=unit(0.2, "npc"), w=unit(1.1, "npc"), angle=0 )
        pushViewport( vp5 )
        #grid.rect()
        grid.draw( ggplotGrob( profileVisualization ) )
    }
    
    
    #fn <- sprintf("tree_phenotypes_regression.ranges.%s_vs_dLFE.out.pdf", Xtrait)
    #ggsave(fn, plot=p )
    #unlink(fn)
    ##fn <- sprintf("tree_phenotypes_regression.ranges.%s_vs_dLFE.out.svg", Xtrait)
    ##ggsave(fn, plot=p )
    ##unlink(fn)
    ##fn <- sprintf("tree_phenotypes_regression.ranges.%s_vs_dLFE.out.wmf", Xtrait)
    ##ggsave(fn, plot=p )
    ##unlink(fn)
}



getStatsForProfiles <- function( values, positions )
{
    stats <- data.frame( pos=integer(), iqr0=double(), p25=double(), p50=double(), p75=double(), iqr1=double() )
    outliers <- data.frame( x=double(), y=double() )

    for( rowId in 1:ncol(values))
    {
        xs <- values[,rowId]
        fit <- kde1d( xs )
        quarts <- qkde1d(c(0.25,0.5,0.75), fit)
        iqr <- quarts[3]-quarts[1]
        iqr.vals <- c(quarts[1]-1.5*iqr, quarts[3]+1.5*iqr)

        stats <- rbind( stats, data.frame( pos=positions[rowId], iqr0=iqr.vals[1], p25=quarts[1], p50=quarts[2], p75=quarts[3], iqr1=iqr.vals[2] ) )

        outliers <- c( xs[xs > iqr.vals[2]], xs[xs < iqr.vals[1]] )

        for( i in outliers )
        {
            outliers <- rbind( outliers, data.frame( x=rowId, y=i ) )
        }
    }
    return( list( stats, outliers ) )
}

profileColorScale <- function( vals )
{
    return( vapply( vals, function  (x) (x+scale.dLFE)/(2*scale.dLFE), FUN.VALUE=double(1) ) )
}


expand_scale <- function(mult = 0, add = 0) {
  stopifnot(is.numeric(mult) && is.numeric(add))
  stopifnot((length(mult) %in% 1:2) && (length(add) %in% 1:2))

  mult <- rep(mult, length.out = 2)
  add <- rep(add, length.out = 2)
  c(mult[1], add[1], mult[2], add[2])
}

kdePlot <- function( values, positions, profileId=1, yrange=NA )
{
    x <- getStatsForProfiles( values, positions )
    stats    <- x[[1]]
    outliers <- x[[2]]

    print(positions)

    interpolated <- approx( stats$pos, stats$p50, method="linear", n=1000 )
    interpolated$color = profileColorScale( interpolated$y )

    if( isTRUE( is.na(yrange) ) )
    {
        min.y <- min( stats$iqr0 )
        max.y <- max( stats$iqr1 )
    }
    else
    {
        stopifnot( yrange[1] < yrange[2] )
        min.y <- yrange[1]
        max.y <- yrange[2]
        print( sprintf("(%f, %f)", min( stats$iqr0 ), max( stats$iqr1 ) ) )
    }
    lines.y <- seq(ceiling(min.y), floor(max.y), 0.5)
                                        #lines.y = lines.y[lines.y != 0]

    if( profileId==1 )
    {
        xtrans <- "identity"
    }
    else
    {
        xtrans <- "reverse"
    }
    

    # create the plot...
    p <- ggplot( data=stats, aes(x=pos) ) +
        #geom_hline( yintercept=lines.y, color="#808080", size=0.3 ) +
        geom_ribbon( aes(ymin=iqr0, ymax=iqr1), color=NA, fill="#a0a0a0" ) +
        geom_ribbon( aes(ymin=p25, ymax=p75),   color=NA, fill="#707070" ) +
        geom_point( data=data.frame(interpolated), aes( x=x, y=y, color=color ), size=0.7 ) +
            guides( color=FALSE ) +  # Hide the legend (will be plotted separately)
            labs(y="", x="") +
        scale_color_gradientn(
            breaks=c(0.0, 0.5, 1.0),
            limits=c(0.0, 1.0),
            colors=        c(                                   "#ff9090", "#a00000", "#ffffff", "#0000a0", "#9090ff"),
            values=rescale(c(                                         0.0,      0.35,       0.5,       0.65,      1.0) ) ) + 
            scale_y_continuous( limits=c( min(min.y, -0.5), max(max.y, 0.5)), breaks=lines.y ) +
            scale_x_continuous( expand=expand_scale( mult=c(0,0) ), trans=xtrans ) +
            geom_hline( yintercept=0, color="#404040", alpha=0.8 ) +
            theme( plot.background = element_blank(),   # Hide unnecessary theme elements (background panels, etc.)
               panel.grid.major.y = element_line(color="grey", size=0.50),
               panel.grid.major.x = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank()
               ) # +
            #scale_x_discrete( labels=stats$pos )
            
        #geom_line( aes(y=p50) );
    #print(p)
    return( p );
}

# Perform regression for pyramid (range robustness) plot
glsRegressionRangeAnalysis <- function( traits, tree, Xtrait, profileRange, plotCaption, profileId=1, addKdePlot=TRUE, extras="" )
{
    stopifnot(class(profileRange)=="numeric")
    stopifnot(length(profileRange)==2)
    stopifnot(profileRange[2] >= profileRange[1])
    # Generate indices for all valid subranges in the given range
    range <- profileRange[1]:profileRange[2]
    subranges <- expand.grid( range, range )
    subranges <- subranges[ subranges$Var2 >= subranges$Var1, ]
    stopifnot(all(subranges$Var2 - subranges$Var1 >= 0))
    stopifnot(nrow(subranges) == length(range)*(length(range)+1)/2 ) # arithmetic series sum

    out <- data.frame( Var1=integer(), Var2=integer(), Pvalue=double(), Buse.R2=double(), Logpval=double(), NumSpecies=integer(), DirectionsIndicator=character() )
    profilePositions <- getProfilePositions( profileId )


    handleSubrange <- function(v)   # Perform regression on each subrange
    {
        out <- get('out', parent.frame())
        #regressionResults <- apply( subranges, 1, function(v) performGLSregression_profileRangeMean(traits, tree, Xtrait, v) )
        #print("-->>")
        stopifnot(class(v)=="data.frame")
        stopifnot(v$Var2 >= v$Var1)

        if( v$Var2 - v$Var1 <= maxPyramidHeight )
        {
            regressionResults <- performGLSregression_profileRangeMean(traits, tree, Xtrait, as.integer(v), profileId=profileId, plotRegression=FALSE, extras=extras )
            v$Pvalue    <- regressionResults$pvalue
            v$Buse.R2   <- regressionResults$R2
            v$Logpval <- log10( v$Pvalue )
            #v$Var1 <- (v$Var1 - 1) * profileStep
            v$Var1 <- profilePositions[v$Var1]
            #v$Var2 <- (v$Var2 - 1) * profileStep
            v$Var2 <- profilePositions[v$Var2]
            v$NumSpecies <- regressionResults$N
            v$DirectionsIndicator <- regressionResults$directionsIndicator
            v$MIC <- regressionResults$MIC
            v$MIC.pvalue <- regressionResults$MIC.pvalue
            v$MAS <- regressionResults$MAS
            v$pearson.r <- regressionResults$pearson.r
            v$BP.pvalue <- regressionResults$BP.pvalue

            #print("regressionResults")
            #print(regressionResults)
            #print("out")
            #print(out)
            #print("v")
            #print(v)
            out <- rbind( out, v )
            assign('out', out, parent.frame() )
        }
        ## else
        ## {
        ##     regressionResults <- data.frame( pvalue=c(NA), R2=c(NA), N=c(NA), directionsIndicator=c(NA) )
        ## }
        
        
        #print("->>->")
        #print(out)
        #print(search())
        
    }
    #by( subranges, 1:nrow(subranges), function(v) handleSubrange(v) )
    for( i in 1:nrow(subranges) ) handleSubrange( subranges[i,] )
    
    subranges <- out
    #print(subranges)

    #minLogpval <- min( subranges$Logpval, -3 )


    dotsDF <- data.frame(melt( subranges, id.vars=c("Var1", "Var2"), measure.vars=c("Logpval") ))
    dotsDF <- dotsDF[dotsDF$value<=-2.0,]


    plotTitle <- sprintf("%s effect on %s (%s)", Xtrait, profileMode, plotCaption)
    
    #grid.newpage()
    #plotdf <- melt( subranges, id.vars=c("Var1", "Var2"), measure.vars=c("Logpval") )
    #p <- ggplot( data=plotdf, aes(x=Var1, y=Var2, fill=value) ) +
    #    geom_raster() +
    #    labs(y="Profile End", x="Profile Start", value="log(p-value)", title=plotTitle) +
    #    scale_fill_gradient2(low="#cc7777", mid="#ffff88", high="#000000", limits=c(min(subranges$Logpval), 0), midpoint=minLogpval ) + 
    #    coord_fixed() +
    #    geom_point(data=plotdf[plotdf$value<=-2.0,], aes(x=Var1,y=Var2), position=position_nudge(y=4, x=-4), color="white", size=0.2 ) +
    #    theme( plot.margin=unit(c(2,2,2,2), "cm"), aspect.ratio=1 )
    #print( p, vp=viewport(angle=-45) )

    y.min <- min( min(subranges$Logpval), -1 )

    valScale <- scale_fill_gradientn(
        breaks=seq( floor( y.min ), 0 ),
        limits=c( floor(y.min), 0),
        colors=c("#7799ff", "#355080", "#000000" ),
        values=rescale(c(y.min, -1.0, 0.0)) )

    #plotdf <- melt( subranges, id.vars=c("Var1", "Var2"), measure.vars=c("Logpval") )
    
    ## plotRotatedPyramid( plotdf, dotsDF, "log(p-val)", valScale, Xtrait, paste(plotTitle, " (logPval)") )

    #----------------------------------------------
    
    plotdf <- data.frame( melt( subranges, id.vars=c("Var1", "Var2"), measure.vars=c("Buse.R2"   ) ) )

    valScale <- NA
    if( any(class(traits[,Xtrait])==c("factor")) )
    {
        valScale <- scale_fill_gradientn(
            breaks=c(0.0, 0.5, 1.0),
            limits=c(0.0, 1.0),
            colors=        c(                                   "#000000", "#602011", "#cc4022", "#ffc075"),
            values=rescale(c(                                         0.0,      0.15,       0.5,       1.0) ) )
    }
    else
    {
        valScale <- scale_fill_gradientn(
            breaks=c(-1.0, -0.5, 0.0, 0.5, 1.0),
            limits=c(-1.0, 1.0),
            colors=        c("#75c0ff", "#2240cc", "#112060",   "#000000", "#602011", "#cc4022", "#ffc075"),
            values=rescale(c(     -1.0,      -0.5,     -0.15,         0.0,      0.15,       0.5,       1.0) ) )
    }


    if( addKdePlot )
    {
        print( getProfileVariables( profileRange, profileId ) )
        kdePlot <- kdePlot( traits[, getProfileVariables( profileRange, profileId ) ], getProfilePositions( profileId )  )
    }
    else
    {
        kdePlot <- NA
    }
       
    plotRotatedPyramid( plotdf, dotsDF, "Buse R^2", valScale, Xtrait, plotTitle=paste(plotTitle, " (Buse R^2)"), kdePlot )


    #----------------------------------------------
    
    #plotdf <- data.frame( melt( subranges, id.vars=c("Var1", "Var2"), measure.vars=c("DirectionsIndicator"   ) ) )
    #print(levels(plotdf$value))
    #levels(plotdf$value) <- as.factor(1:nlevels(plotdf$value))
    #print(levels(plotdf$value))
    #print(plotdf)

    #plotdf$value <- as.factor(plotdf$value)

    #valScale <- scale_colour_manual( breaks=c("1","2"), labels=c("1","2"), values=c("#7799f0", "#f06666") )
    #valScale <- scale_fill_hue( ) # plotdf, aes(fill=factor(value) ) )

    ## plotRotatedPyramid( plotdf, dotsDF, "Directions", valScale, Xtrait, plotTitle=paste(plotTitle, " (order)") )


    # Plot MIC pyramid
    if( FALSE )
    {
        plotdf <- data.frame( melt( subranges, id.vars=c("Var1", "Var2"), measure.vars=c("MIC"  ) ) )
        print(plotdf)
        dotsDF <- data.frame( melt( subranges, id.vars=c("Var1", "Var2"), measure.vars=c("MIC.pvalue") ))
        dotsDF <- dotsDF[dotsDF$value < significanceLevel,]
        plotRotatedPyramid( plotdf, dotsDF, "MIC", valScale, Xtrait, plotTitle=paste(plotTitle, " (MIC)") )
    }
    
    
    return(subranges)
}


olsRegressionRangeAnalysis <- function( traits, Xtrait, profileRange, plotCaption, profileId=1, addKdePlot=TRUE, extras="" )
{
    stopifnot(class(profileRange)=="numeric")
    stopifnot(length(profileRange)==2)
    stopifnot(profileRange[2] >= profileRange[1])
    # Generate indices for all valid subranges in the given range
    range <- profileRange[1]:profileRange[2]
    subranges <- expand.grid( range, range )
    subranges <- subranges[ subranges$Var2 >= subranges$Var1, ]
    stopifnot(all(subranges$Var2 - subranges$Var1 >= 0))
    stopifnot(nrow(subranges) == length(range)*(length(range)+1)/2 ) # arithmetic series sum

    out <- data.frame( Var1=integer(), Var2=integer(), Pvalue=double(), Buse.R2=double(), Logpval=double(), NumSpecies=integer(), DirectionsIndicator=character() )

    handleSubrange <- function(v)   # Perform regression on each subrange
    {
        out <- get('out', parent.frame())
        #regressionResults <- apply( subranges, 1, function(v) performGLSregression_profileRangeMean(traits, tree, Xtrait, v) )
        #print("-->>")
        stopifnot(class(v)=="data.frame")
        stopifnot(v$Var2 >= v$Var1)

        if( v$Var2 - v$Var1 <= maxPyramidHeight )
        {
            regressionResults <- performOLSregression_profileRangeMean(traits, Xtrait, as.integer(v), profileId=profileId, plotRegression=FALSE, extras=extras )
            v$Pvalue    <- regressionResults$pvalue
            v$Buse.R2   <- regressionResults$R2
            v$Logpval <- log10( v$Pvalue )
            v$Var1 <- (v$Var1 - 1) * profileStep
            v$Var2 <- (v$Var2 - 1) * profileStep
            v$NumSpecies <- regressionResults$N
            v$DirectionsIndicator <- regressionResults$directionsIndicator
            v$MIC <- regressionResults$MIC
            v$MIC.pvalue <- regressionResults$MIC.pvalue
            v$MAS <- regressionResults$MAS
            v$pearson.r <- regressionResults$pearson.r
            v$BP.pvalue <- regressionResults$BP.pvalue

            #print("regressionResults")
            #print(regressionResults)
            #print("out")
            #print(out)
            #print("v")
            #print(v)
            out <- rbind( out, v )
            assign('out', out, parent.frame() )
        }
        ## else
        ## {
        ##     regressionResults <- data.frame( pvalue=c(NA), R2=c(NA), N=c(NA), directionsIndicator=c(NA) )
        ## }
        
        
        #print("->>->")
        #print(out)
        #print(search())
        
    }
    #by( subranges, 1:nrow(subranges), function(v) handleSubrange(v) )
    for( i in 1:nrow(subranges) ) handleSubrange( subranges[i,] )
    
    subranges <- out
    #print(subranges)

    minLogpval <- min( subranges$Logpval, -3 )


    dotsDF <- data.frame(melt( subranges, id.vars=c("Var1", "Var2"), measure.vars=c("Logpval") ))
    dotsDF <- dotsDF[dotsDF$value<=-2.0,]


    plotTitle <- sprintf("%s effect on %s (%s)", Xtrait, profileMode, plotCaption)
    
    #grid.newpage()
    #plotdf <- melt( subranges, id.vars=c("Var1", "Var2"), measure.vars=c("Logpval") )
    #p <- ggplot( data=plotdf, aes(x=Var1, y=Var2, fill=value) ) +
    #    geom_raster() +
    #    labs(y="Profile End", x="Profile Start", value="log(p-value)", title=plotTitle) +
    #    scale_fill_gradient2(low="#cc7777", mid="#ffff88", high="#000000", limits=c(min(subranges$Logpval), 0), midpoint=minLogpval ) + 
    #    coord_fixed() +
    #    geom_point(data=plotdf[plotdf$value<=-2.0,], aes(x=Var1,y=Var2), position=position_nudge(y=4, x=-4), color="white", size=0.2 ) +
    #    theme( plot.margin=unit(c(2,2,2,2), "cm"), aspect.ratio=1 )
    #print( p, vp=viewport(angle=-45) )

    y.min <- min( min(subranges$Logpval), -1 )

    valScale <- scale_fill_gradientn(
        breaks=seq( floor( y.min ), 0 ),
        limits=c( floor(y.min), 0),
        colors=c("#7799ff", "#355080", "#000000" ),
        values=rescale(c(y.min, -1.0, 0.0)) )

    #plotdf <- melt( subranges, id.vars=c("Var1", "Var2"), measure.vars=c("Logpval") )
    
    ## plotRotatedPyramid( plotdf, dotsDF, "log(p-val)", valScale, Xtrait, paste(plotTitle, " (logPval)") )

    #----------------------------------------------
    
    plotdf <- data.frame( melt( subranges, id.vars=c("Var1", "Var2"), measure.vars=c("Buse.R2"   ) ) )

    valScale <- NA
    if( any(class(traits[,Xtrait])==c("factor")) )
    {
        valScale <- scale_fill_gradientn(
            breaks=c(0.0, 0.5, 1.0),
            limits=c(0.0, 1.0),
            colors=        c(                                   "#000000", "#602011", "#cc4022", "#ffc075"),
            values=rescale(c(                                         0.0,      0.15,       0.5,       1.0) ) )
    }
    else
    {
        valScale <- scale_fill_gradientn(
            breaks=c(-1.0, -0.5, 0.0, 0.5, 1.0),
            limits=c(-1.0, 1.0),
            colors=        c("#75c0ff", "#2240cc", "#112060",   "#000000", "#602011", "#cc4022", "#ffc075"),
            values=rescale(c(     -1.0,      -0.5,     -0.15,         0.0,      0.15,       0.5,       1.0) ) )
    }


    if( addKdePlot )
    {
        print( getProfileVariables( profileRange, profileId ) )
        kdePlot <- kdePlot( traits[, getProfileVariables( profileRange, profileId ) ], getProfilePositions( profileId )  )
    }
    else
    {
        kdePlot <- NA
    }
       
    plotRotatedPyramid( plotdf, dotsDF, "Buse R^2", valScale, Xtrait, plotTitle=paste(plotTitle, " (Buse R^2)"), kdePlot )


    #----------------------------------------------
    
    #plotdf <- data.frame( melt( subranges, id.vars=c("Var1", "Var2"), measure.vars=c("DirectionsIndicator"   ) ) )
    #print(levels(plotdf$value))
    #levels(plotdf$value) <- as.factor(1:nlevels(plotdf$value))
    #print(levels(plotdf$value))
    #print(plotdf)

    #plotdf$value <- as.factor(plotdf$value)

    #valScale <- scale_colour_manual( breaks=c("1","2"), labels=c("1","2"), values=c("#7799f0", "#f06666") )
    #valScale <- scale_fill_hue( ) # plotdf, aes(fill=factor(value) ) )

    ## plotRotatedPyramid( plotdf, dotsDF, "Directions", valScale, Xtrait, plotTitle=paste(plotTitle, " (order)") )


    # Plot MIC pyramid
    plotdf <- data.frame( melt( subranges, id.vars=c("Var1", "Var2"), measure.vars=c("MIC"  ) ) )
    print(plotdf)
    dotsDF <- data.frame( melt( subranges, id.vars=c("Var1", "Var2"), measure.vars=c("MIC.pvalue") ))
    dotsDF <- dotsDF[dotsDF$value < significanceLevel,]
    plotRotatedPyramid( plotdf, dotsDF, "MIC", valScale, Xtrait, plotTitle=paste(plotTitle, " (MIC)") )
    
    
    
    return(subranges)
}

findUncorrectedVariableCorrelation <- function( traits, filterTrait, Xtrait, Ytrait )
{
    traits <- traits[ (traits[,filterTrait] & (!is.na(traits[,Xtrait])) & (!is.na(traits[,Ytrait]))  ), ]

    d1 <- traits[,Xtrait]
    d2 <- traits[,Ytrait]
    N <- sum( !(is.na(d1) | is.na(d2)) )
    ret <- cor.test( d1, d2, method="spearman", rm.na=TRUE, exact=TRUE )
    #print(ret$estimate)
    return( data.frame( pvalue=ret$p.value, corr=ret$estimate, N=N, directionsIndicator=NA ) )
}



#print(performGLSregression( traits, tree, "GenomicGC",   "Profile_1.15", extras="", plotRegression=FALSE ))
#print(performGLSregression( traits, tree, "GenomicGC",   "Profile_1.15", extras=" + GenomeSizeMb", plotRegression=FALSE ) )
#print(performGLSregression( traits, tree, "GenomicGC",   "Profile_1.15", extras=" + GenomeSizeMb + GenomicCAI", plotRegression=FALSE ) )
#print(performGLSregression( traits, tree, "GenomicGC",   "Profile_1.15", extras=" + GenomeSizeMb + GenomicCAI + GenomicCBI + GenomicFop", plotRegression=FALSE ) )
#print(performGLSregression( traits, tree, "GenomicGC",   "Profile_1.15", extras=" + GenomeSizeMb + GenomicCAI + GenomicCBI + GenomicFop + Noise", plotRegression=FALSE ) )
#print(performGLSregression( traits, tree, "GenomicGC",   "Profile_1.15", extras=" + GenomeSizeMb + GenomicCAI + GenomicCBI + GenomicFop + Profile_1.14", plotRegression=FALSE ) )

TestCorrelationBetweenScalarTraits <- function()
{
    
    print(cor( traits$GenomicGC, traits$GenomicCAI, use="complete.obs"  ) )
    print(cor( traits$GenomicGC, traits$GenomicCBI, use="complete.obs"  ) )
    print(cor( traits$GenomicGC, traits$GenomicFop, use="complete.obs"  ) )
    print(cor( traits$GenomicGC, traits$GenomicNc,  use="complete.obs"  ) )

    print(cor( traits$GenomicCAI, traits$GenomicCBI, use="complete.obs"  ) )
    print(cor( traits$GenomicCAI, traits$GenomicFop, use="complete.obs"  ) )
    print(cor( traits$GenomicCBI, traits$GenomicFop, use="complete.obs"  ) )
    print(cor( traits$GenomicCBI, traits$GenomicNc,  use="complete.obs"  ) )

    print(cor( traits$GenomicENc, traits$GenomicNc,   use="complete.obs"  ) )
    print(cor( traits$GenomicGC,  traits$GenomicNc,   use="complete.obs"  ) )
    print(cor( traits$GenomicGC,  traits$GenomicENc,  use="complete.obs"  ) )
    print(cor( traits$GenomicGC,  traits$GenomicENc.prime,  use="complete.obs"  ) )
    #print(cor( traits.normalized$GenomicENc.prime, traits.normalized$Profile_1.sd, use="complete.obs"  ) )
    #print(cor( traits.normalized$GenomicENc.prime, traits.normalized$Profile_2.sd, use="complete.obs"  ) )

}

#------------------------------------------------------------------------------------------------------------------------
# Compare R^2 of GC and ENc' at different points in the native and "sigma-normalized" profiles

## print(performGLSregression( traits, tree, "GenomicENc.prime",   "Profile_1.26", extras="", plotRegression=FALSE ) )
## print(performGLSregression( traits, tree, "GenomicENc.prime",   "Profile_1.26", extras=" GenomicGC ", plotRegression=FALSE ) )

## print(performGLSregression( traits, tree, "GenomicGC",   "Profile_1.26", extras="", plotRegression=FALSE ) )
## print(performGLSregression( traits, tree, "GenomicGC",   "Profile_1.26", extras=" GenomicENc.prime ", plotRegression=FALSE ) )


## print(performGLSregression( traits.normalized, tree, "GenomicGC",         "Profile_1.sd",  extras="", plotRegression=FALSE ) )
## print(performGLSregression( traits.normalized, tree, "GenomicENc.prime",  "Profile_1.sd",  extras="", plotRegression=FALSE ) )
## print(performGLSregression( traits.normalized, tree, "GenomicENc.prime",  "Profile_1.sd",  extras="GenomicGC", plotRegression=FALSE ) )
## print(performGLSregression( traits.normalized, tree, "GenomicCAI",        "Profile_1.sd",  extras="", plotRegression=FALSE ) )
## print(performGLSregression( traits.normalized, tree, "GenomicDCBS",       "Profile_1.sd",  extras="", plotRegression=FALSE ) )


## print(performGLSregression( traits, tree, "GenomicGC",                "Profile_1.26",  extras="", plotRegression=TRUE ) )
## print(performGLSregression( traits, tree, "GenomicENc.prime",         "Profile_1.26",  extras="", plotRegression=TRUE ) )



## print(performGLSregression( traits.normalized, tree, "GenomicENc.prime",  "Profile_2.sd",  extras="", plotRegression=FALSE ) )
## print(performGLSregression( traits.normalized, tree, "GenomicENc.prime",  "Profile_2.sd",  extras="GenomicGC", plotRegression=FALSE ) )
## print(performGLSregression( traits.normalized, tree, "GenomicCAI",        "Profile_2.sd",  extras="", plotRegression=FALSE ) )
## print(performGLSregression( traits.normalized, tree, "GenomicDCBS",       "Profile_2.sd",  extras="", plotRegression=FALSE ) )

## print(performGLSregression( traits           , tree, "GenomicGC",          "Profile_1.1",  extras="", plotRegression=FALSE ) )
## print(performGLSregression( traits.normalized, tree, "GenomicGC",          "Profile_1.1",  extras="", plotRegression=FALSE ) )
## print(performGLSregression( traits           , tree, "GenomicENc.prime",   "Profile_1.1",  extras="", plotRegression=FALSE ) )
## print(performGLSregression( traits.normalized, tree, "GenomicENc.prime",   "Profile_1.1",  extras="", plotRegression=FALSE ) )
## print(performGLSregression( traits           , tree, "GenomicENc.prime",   "Profile_1.1",  extras="GenomicGC", plotRegression=FALSE ) )

## print(performGLSregression( traits           , tree, "GenomicGC",          "Profile_1.15", extras="", plotRegression=FALSE ) )
## print(performGLSregression( traits.normalized, tree, "GenomicGC",          "Profile_1.15", extras="", plotRegression=FALSE ) )
## print(performGLSregression( traits           , tree, "GenomicENc.prime",   "Profile_1.15", extras="", plotRegression=FALSE ) )
## print(performGLSregression( traits.normalized, tree, "GenomicENc.prime",   "Profile_1.15", extras="", plotRegression=FALSE ) )
## print(performGLSregression( traits           , tree, "GenomicENc.prime",   "Profile_1.15", extras="GenomicGC", plotRegression=FALSE ) )

## print(performGLSregression( traits           , tree, "GenomicGC",          "Profile_1.22", extras="", plotRegression=FALSE ) )
## print(performGLSregression( traits.normalized, tree, "GenomicGC",          "Profile_1.22", extras="", plotRegression=FALSE ) )
## print(performGLSregression( traits           , tree, "GenomicENc.prime",   "Profile_1.22", extras="", plotRegression=FALSE ) )
## print(performGLSregression( traits.normalized, tree, "GenomicENc.prime",   "Profile_1.22", extras="", plotRegression=FALSE ) )
## print(performGLSregression( traits           , tree, "GenomicENc.prime",   "Profile_1.22", extras="GenomicGC", plotRegression=FALSE ) )

## print(performGLSregression( traits           , tree, "GenomicGC",          "Profile_1.26", extras="", plotRegression=FALSE ) )
## print(performGLSregression( traits.normalized, tree, "GenomicGC",          "Profile_1.26", extras="", plotRegression=FALSE ) )
## print(performGLSregression( traits           , tree, "GenomicENc.prime",   "Profile_1.26", extras="", plotRegression=FALSE ) )
## print(performGLSregression( traits.normalized, tree, "GenomicENc.prime",   "Profile_1.26", extras="", plotRegression=FALSE ) )
## print(performGLSregression( traits           , tree, "GenomicENc.prime",   "Profile_1.26", extras="GenomicGC", plotRegression=FALSE ) )

## print(performGLSregression( traits           , tree, "GenomicGC",          "Profile_1.29", extras="", plotRegression=FALSE ) )
## print(performGLSregression( traits.normalized, tree, "GenomicGC",          "Profile_1.29", extras="", plotRegression=FALSE ) )
## print(performGLSregression( traits           , tree, "GenomicENc.prime",   "Profile_1.29", extras="", plotRegression=FALSE ) )
## print(performGLSregression( traits.normalized, tree, "GenomicENc.prime",   "Profile_1.29", extras="", plotRegression=FALSE ) )
## print(performGLSregression( traits           , tree, "GenomicENc.prime",   "Profile_1.29", extras="GenomicGC", plotRegression=FALSE ) )
##quit()
#------------------------------------------------------------------------------------------------------------------------

## print(performGLSregression( traits, tree, "GenomicGC",   "Profile_1.15", extras=" + GenomeSizeMb", plotRegression=FALSE ) )
## print(performGLSregression( traits, tree, "GenomicGC",   "Profile_1.15", extras=" + GenomeSizeMb + GenomicCAI + GenomicCBI + GenomicFop", plotRegression=FALSE ) )
## print(performGLSregression( traits, tree, "GenomicGC",   "Profile_1.15", extras=" + GenomeSizeMb + Noise", plotRegression=FALSE ) )
## print(performGLSregression( traits, tree, "GenomicGC",   "Profile_1.15", extras=" + GenomeSizeMb + Profile_1.14", plotRegression=FALSE ) )

## print(performGLSregression( traits, tree, "GenomeSizeMb",   "Profile_1.15", extras=" + GenomicCAI + GenomicCBI + GenomicFop", plotRegression=FALSE ) )
## print(performGLSregression( traits, tree, "GenomeSizeMb",   "Profile_1.15", extras=" + GenomicCAI + GenomicCBI + GenomicFop + GenomicGC", plotRegression=FALSE ) )
## print(performGLSregression( traits, tree, "GenomeSizeMb",   "Profile_1.15", extras=" + GenomicCAI + GenomicCBI + GenomicFop + Noise", plotRegression=FALSE ) )
## print(performGLSregression( traits, tree, "GenomeSizeMb",   "Profile_1.15", extras=" + GenomicCAI + GenomicCBI + GenomicFop + Profile_1.14", plotRegression=FALSE ) )

## print(performGLSregression( traits, tree, "GenomicGC",   "Profile_1.15", extras=" + GenomicCAI + GenomicCBI + GenomicFop", plotRegression=FALSE ) )
## print(performGLSregression( traits, tree, "GenomicGC",   "Profile_1.15", extras=" + GenomicCAI + GenomicCBI + GenomicFop + GenomeSizeMb", plotRegression=FALSE ) )
## print(performGLSregression( traits, tree, "GenomicGC",   "Profile_1.15", extras=" + GenomicCAI + GenomicCBI + GenomicFop + Noise", plotRegression=FALSE ) )
## print(performGLSregression( traits, tree, "GenomicGC",   "Profile_1.15", extras=" + GenomicCAI + GenomicCBI + GenomicFop + Profile_1.14", plotRegression=FALSE ) )


## print(performGLSregression( traits, tree, "GenomicCAI",   "GenomicGC", extras=" + GenomicCBI + GenomicFop", plotRegression=FALSE ) )
## print(performGLSregression( traits, tree, "GenomicCBI",   "GenomicGC", extras=" + GenomicCAI + GenomicFop", plotRegression=FALSE ) )
## print(performGLSregression( traits, tree, "GenomicFop",   "GenomicGC", extras=" + GenomicCAI + GenomicCBI", plotRegression=FALSE ) )

## print(performGLSregression( traits, tree, "GenomicGC",   "Profile_1.15", extras="", plotRegression=FALSE ) )




#print(performGLSregression( traits, tree, "GenomicGC",   "Profile_1.15", extras=" + GenomeSizeMb + GenomicCAI + GenomicCBI + GenomicFop + GenomicENc.prime", plotRegression=FALSE ) ))
#print(performGLSregression_profileRangeMean( traits, tree, "GenomicGC", c(15,15) ))
#print(performGLSregression_profileRangeMean( traits, tree, "GenomicGC", c(15,19) ))


traitComparisonAnalysis <- function( allTraits, tree, traitsToAnalyze, profileRange, profileId, plotRegressions=TRUE )
{
    stopifnot(length(profileRange)==2)
    stopifnot(profileRange[2]>profileRange[1])
    data <- data.frame()

    pos <- 1

    for( t1 in traitsToAnalyze )
    {
        
        reg.v1 <- performGLSregression_profileRangeMean( allTraits, tree, t1, profileRange, profileId=profileId, plotRegression=plotRegressions, extras="" )
        
        for( t2 in traitsToAnalyze )
        {
            if( t2==t1 ) next;
            
            print(sprintf("%s %s", t1, t2))
            
            #reg.v1.v2 <- performGLSregression_profileRangeMean( allTraits, tree, t1, profileRange, profileId=profileId, plotRegression=plotRegressions, extras=sprintf(" %s ", t2) )
            reg.v1.v2 <- performOLSregression_profileRangeMean( allTraits,        t1, profileRange, profileId=profileId, plotRegression=plotRegressions, extras=sprintf(" %s ", t2) )
            if( !all(is.na(reg.v1.v2) ) )
            {
                pos <- pos+1
                reg.v1.v2$Trait1 <- t1
                reg.v1.v2$Trait2 <- t2
                reg.v1.v2$R2.1 <- reg.v1$R2
                reg.v1.v2$R2.2 <- reg.v1.v2$R2
                reg.v1.v2$pos <- pos
                print(reg.v1.v2)
                data <- rbind(data, reg.v1.v2)
                print("~")
            }
        }
    }

    return( data)
}


getBinarizedTrait <- function( vals, threshold, threshold.logic )
{
    ret <- NA
    if( threshold.logic==-1 )
    {
        ret <- vals <= threshold
    }
    else if ( threshold.logic==0 )
    {
        ret <- vals == threshold
    }
    else if ( threshold.logic==1 )
    {
        ret <- vals >= threshold
    }
    ret[is.na(ret)] <- FALSE # treat NA values as FALSE
    ret
}


testCompoundClassification <- function( traits, tree, values, values.logic, profileRange, profileId )
{
    stopifnot( length(values)==length(values.logic) )
    
    # Prepare composite predictor
    compound <- rep(FALSE, nrow(traits))
    for( i in 1:length(values) )
    {
        threshold <- values[i]
        var.name  <- names(values[i])
        var.logic <- values.logic[i]
        stopifnot(var.logic %in% c(-1,0,1))
        if( var.name == "Response" ) next

        stopifnot( var.name != "Response" )
        vals <- getBinarizedTrait( traits[,var.name], threshold, var.logic )
        
        compound <- compound | vals  # combining logical op is always "or'
    }
    traits$BinaryClass <- compound

    # Prepare response variable
    variables <- getProfileVariables( profileRange, profileId )
    #traits$AbsRangeMean <- abs( rowMeans( traits[variables] ) )
    traits$SD <- apply(traits[,variables], MARGIN=1, FUN=sd)

    responsePos <- which(names(values)=="Response")[1]
    stopifnot( names(values)[responsePos] == "Response" )
    traits$BinaryResponse <- getBinarizedTrait( traits$SD, values["Response"], values.logic[responsePos] )
    #print( traits$BinaryResponse )

    prec   <- precision( traits$BinaryResponse, traits$BinaryClass )
    recall <- recall(    traits$BinaryResponse, traits$BinaryClass )
    #print( f1(        traits$BinaryResponse, traits$BinaryClass ) )

    #ret <- performGLSregression_profileRangeMean( traits, tree, "BinaryClass", c(16,31), profileId=1, plotRegression=TRUE )
    #ret$params <- paste(sprintf("%s", values), collapse=",")
                                        #ret
    ret <- rbind(data.frame(), values)
    colnames(ret) <- names(values)
    ret$precision <- prec
    ret$recall    <- recall
    ret
}

gridSearch <- function( traits, tree, values.from, values.to, values.logic, grid.steps )
{
    stopifnot(length(values.from)==length(values.to))
    stopifnot(length(values.from)==length(values.logic))
    stopifnot(all(names(values.from)==names(values.to)))
    stopifnot(all(names(values.from)==names(values.logic)))
    # create list with all parameter values
    grid.vals <- list()
    for( i in 1:length(values.from) )
    {
        if( values.from[i] == values.to[i] )
        {
            vals <- c(values.from[i])
        }
        else
        {
            vals <- seq( values.from[i], values.to[i], length.out=grid.steps )
        }
        
        newlist <- list( a=vals )
        names(newlist)[1] <- names(values.from)[i]
        grid.vals <- c(grid.vals, newlist)
    }

    # create a data-frame with all combinations
    grid.combinations <- do.call("expand.grid", grid.vals)

    eval.combination <- function(params)
    {
        testCompoundClassification( traits, tree, values=params, values.logic=values.logic, profileRange=c(1,31), profileId=1 )
    }

    results <- data.frame()
    for( i in 1:nrow(grid.combinations) )
    {
        params = as.numeric( grid.combinations[i,] )
        names(params) <- colnames( grid.combinations )
        res <- eval.combination(params)
        results <- rbind( results, res )
    }
    results
}

report_testCompoundClassification <- function()
{
    

    # This is the "official" classification (chosen for use in python plotting code)
    #37.50	56.21	58.01	0.14	0.6550	0.8242
    #38.40	56.61	58.01	0.14	0.6550	0.8242
    print("\"Official\" classification results:")
    print(testCompoundClassification( traits, tree, values=c(Is_endosymbiont=1, GenomicGC=37.99, GenomicENc.prime=56.51, OptimumTemp=58.01, Response=0.14 ), values.logic=c(0, -1, 1, 1, -1), profileRange=c(1,31), profileId=1 ) )

    #print(testCompoundClassification( traits, tree, values=c(Is_endosymbiont=0, GenomicGC=34.99, GenomicENc.prime=55.01, OptimumTemp=65.01, Response=0.14 ), values.logic=c(0, -1, 1, 1, -1), profileRange=c(16,31), profileId=1 ) )
    #print(testCompoundClassification( traits, tree, values=c(Is_endosymbiont=1, GenomicGC=32.10, GenomicENc.prime=56.60, OptimumTemp=74.01, Response=0.14 ), values.logic=c(0, -1, 1, 1, -1), profileRange=c(16,31), profileId=1 ) )
    #print(testCompoundClassification( traits, tree, values=c(Is_endosymbiont=1, GenomicGC=32.10, GenomicENc.prime=56.60, OptimumTemp=74.01, Response=0.12 ), values.logic=c(0, -1, 1, 1, -1), profileRange=c(16,31), profileId=1 ) )
}


writeWeakDLFEBinaryModelGridSearchResults <- function()
{
    
    ## write.csv(
    ##     gridSearch( traits, tree,
    ##                values.from= c(Is_endosymbiont=1, GenomicGC=30.99, GenomicENc.prime=52.01, OptimumTemp=50.01, Response=0.09 ),
    ##                values.to=   c(Is_endosymbiont=1, GenomicGC=38.99, GenomicENc.prime=58.01, OptimumTemp=80.01, Response=0.09 ),
    ##                values.logic=c(Is_endosymbiont=0, GenomicGC=-1,    GenomicENc.prime=1,     OptimumTemp=1,     Response=-1 ),
    ##               grid.steps=11 ),
    ##     file="weak_dlfe_binary_classifier.csv" )

    #37.50	56.21	58.01	0.14	0.6550	0.8242
    #38.40	56.61	58.01	0.14	0.6550	0.8242
    write.csv(
        gridSearch( traits, tree,
                   values.from= c(Is_endosymbiont=1, GenomicGC=36.00, GenomicENc.prime=55.01, OptimumTemp=58.01, Response=0.14 ),
                   values.to=   c(Is_endosymbiont=1, GenomicGC=39.00, GenomicENc.prime=59.01, OptimumTemp=80.01, Response=0.14 ),
                   values.logic=c(Is_endosymbiont=0, GenomicGC=-1,    GenomicENc.prime=1,     OptimumTemp=1,     Response=-1 ),
                   grid.steps=11 ),
        file="weak_dlfe_binary_classifier.csv" )

}
##dev.off()
##quit()

#performGLSregression_profileRangeMean( traits, tree, "Is_endosymbiont", c(16,31), profileId=1, extras=" Is_high_temp ")

#tt2t <- traitComparisonAnalysis( traits, tree, c("Is_endosymbiont", "Is_low_GC", "Is_high_ENc_prime", "Is_high_temp"), c(16,31), profileId=1 )
#print(tt2t)

#print(melt(tt2t, measure.vars=c("R2.2")))


#dev.off()
#quit()


figure_GLS_GC_vs_endosymbionts <- function()
{
    performGLSregression( traits, tree, "Is_endosymbiont", "GenomicGC",  plotRegression=TRUE )
    performGLSregression( traits, tree, "GenomicENc.prime", "GenomicGC", plotRegression=TRUE )
}


report_testRegressionForWeakModelComponents <- function()
{
##print(performGLSregression( traits, tree, "GenomicGC",   "Profile_1.15", extras="", plotRegression=FALSE ))

    print(performGLSregression_profileRangeMean( traits, tree, "Is_endosymbiont", c(1,2), profileId=1, plotRegression=TRUE ))
    print(performOLSregression_profileRangeMean( traits,       "Is_endosymbiont", c(1,2), profileId=1, plotRegression=TRUE ))

    print(performGLSregression_profileRangeMean( traits, tree, "Is_endosymbiont", c(16,31), profileId=1, plotRegression=TRUE ))
    print(performOLSregression_profileRangeMean( traits,       "Is_endosymbiont", c(16,31), profileId=1, plotRegression=TRUE ))

    print(performGLSregression_profileRangeMean( traits, tree, "Is_endosymbiont", c(2,16), profileId=2, plotRegression=TRUE ))
    print(performOLSregression_profileRangeMean( traits,       "Is_endosymbiont", c(2,16), profileId=2, plotRegression=TRUE ))

    print(performGLSregression_profileRangeMean( traits, tree, "Compound", c(16,31), profileId=1, plotRegression=TRUE ))
    print(performOLSregression_profileRangeMean( traits,       "Compound", c(16,31), profileId=1, plotRegression=TRUE ))

    print(performGLSregression_profileRangeMean( traits, tree, "Compound", c(2,16), profileId=2, plotRegression=TRUE ))
    print(performOLSregression_profileRangeMean( traits,       "Compound", c(2,16), profileId=2, plotRegression=TRUE ))



    #traits.fb <- na.exclude(traits[traits$Compound==0,"Profile_1.26"])
    #print(traits.fb)
    #print(nrow(traits.fb))
                                            #print(traits.fb[order(traits.fb$Profile_1.26),"Profile_1.26"])
    #print(traits[order(traits$Profile_1.26),c("Profile_1.26", "Compound", "Is_endosymbiont", "OptimumTemp", "TempHighLow75", "GenomicENc.prime", "GenomicGC")])
    ##dev.off()
    ##quit()


    performGLSregression( traits, tree, "Is_endosymbiont", "GenomicENc.prime", plotRegression=TRUE )
    performOLSregression( traits,       "Is_endosymbiont", "GenomicENc.prime", plotRegression=TRUE )

    performGLSregression( traits, tree, "Is_endosymbiont", "LogGenomeSize",    plotRegression=TRUE )
    performOLSregression( traits,       "Is_endosymbiont", "LogGenomeSize",    plotRegression=TRUE )


    print(performGLSregression_profileRangeMean( traits, tree, "Is_high_ENc_prime",   c(16,31), profileId=1, plotRegression=TRUE ))
    print(performOLSregression_profileRangeMean( traits,       "Is_high_ENc_prime",   c(16,31), profileId=1, plotRegression=TRUE ))
    print(performGLSregression_profileRangeMean( traits, tree, "Is_high_ENc_prime",   c(2,16),  profileId=2, plotRegression=TRUE ))
    print(performOLSregression_profileRangeMean( traits,       "Is_high_ENc_prime",   c(2,16),  profileId=2, plotRegression=TRUE ))

    print(performGLSregression_profileRangeMean( traits, tree, "GenomicENc.prime",   c(16,31), profileId=1, plotRegression=TRUE ))
    print(performOLSregression_profileRangeMean( traits,       "GenomicENc.prime",   c(16,31), profileId=1, plotRegression=TRUE ))
    print(performGLSregression_profileRangeMean( traits, tree, "GenomicENc.prime",   c(2,16),  profileId=2, plotRegression=TRUE ))
    print(performOLSregression_profileRangeMean( traits,       "GenomicENc.prime",   c(2,16),  profileId=2, plotRegression=TRUE ))

    print(performGLSregression_profileRangeMean( traits, tree, "Is_high_temp",       c(16,31), profileId=1, plotRegression=TRUE ))
    print(performOLSregression_profileRangeMean( traits,       "Is_high_temp",       c(16,31), profileId=1, plotRegression=TRUE ))
    print(performGLSregression_profileRangeMean( traits, tree, "Is_high_temp",       c(2,16),  profileId=2, plotRegression=TRUE ))
    print(performOLSregression_profileRangeMean( traits,       "Is_high_temp",       c(2,16),  profileId=2, plotRegression=TRUE ))
}


#dev.off()
#quit()


#-----------------------------------------------------------------------------------------------------------------------------------
# Supp. figure - correlation between traits and separate components of the profile
#-----------------------------------------------------------------------------------------------------------------------------------

#   Var1 Var2       Pvalue   Buse.R2    Logpval NumSpecies DirectionsIndicator          MIC MIC.pvalue        MAS pearson.r    BP.pvalue
#      0    0 2.028521e-63 0.5714344 -62.692821        336                   +    0.6179513          0 0.03328743 0.8063670 5.543088e-07

## rr1 <- data.frame()
## tt11 <- glsRegressionRangeAnalysis( traits, tree, "GenomicGC", c(1,31), "DeltaLFE", profileId=1 )
## tt11$type <- "deltaLFE"
## rr1 <- rbind( rr1, tt11 )
## tt11 <- glsRegressionRangeAnalysis( traits, tree, "GenomicGC", c(1,31), "Native",   profileId=2 )
## tt11$type <- "native"
## rr1 <- rbind( rr1, tt11 )
## tt11 <- glsRegressionRangeAnalysis( traits, tree, "GenomicGC", c(1,31), "Shuffled", profileId=3 )
## tt11$type <- "shuffled"
## rr1 <- rbind( rr1, tt11 )

## p <- ggplot( rr1, aes(x=Var1, y=Buse.R2) ) +
##     geom_hline( yintercept=0, color="black" ) +
##     geom_line( aes( color=type, group=type ), size=3, alpha=0.9 ) +
##     labs( x="CDS position (nt)", y="Buse-R^2 * sign", title="GenomicGC" )
## print(p)

## rr1 <- data.frame()
## tt11 <- glsRegressionRangeAnalysis( traits, tree, "GenomicENc.prime", c(1,31), "DeltaLFE", profileId=1 )
## tt11$type <- "deltaLFE"
## rr1 <- rbind( rr1, tt11 )
## tt11 <- glsRegressionRangeAnalysis( traits, tree, "GenomicENc.prime", c(1,31), "Native",   profileId=2 )
## tt11$type <- "native"
## rr1 <- rbind( rr1, tt11 )
## tt11 <- glsRegressionRangeAnalysis( traits, tree, "GenomicENc.prime", c(1,31), "Shuffled", profileId=3 )
## tt11$type <- "shuffled"
## rr1 <- rbind( rr1, tt11 )

## p <- ggplot( rr1, aes(x=Var1, y=Buse.R2) ) +
##     geom_hline( yintercept=0, color="black" ) +
##     geom_line( aes( color=type, group=type ), size=3, alpha=0.9 ) +
##     labs( x="CDS position (nt)", y="Buse-R^2 * sign", title="GenomicENc.prime" )
## print(p)


## rr1 <- data.frame()
## tt11 <- glsRegressionRangeAnalysis( traits, tree, "GenomicCAI", c(1,31), "DeltaLFE", profileId=1 )
## tt11$type <- "deltaLFE"
## rr1 <- rbind( rr1, tt11 )
## tt11 <- glsRegressionRangeAnalysis( traits, tree, "GenomicCAI", c(1,31), "Native",   profileId=2 )
## tt11$type <- "native"
## rr1 <- rbind( rr1, tt11 )
## tt11 <- glsRegressionRangeAnalysis( traits, tree, "GenomicCAI", c(1,31), "Shuffled", profileId=3 )
## tt11$type <- "shuffled"
## rr1 <- rbind( rr1, tt11 )

## p <- ggplot( rr1, aes(x=Var1, y=Buse.R2) ) +
##     geom_hline( yintercept=0, color="black" ) +
##     geom_line( aes( color=type, group=type ), size=3, alpha=0.9 ) +
##     labs( x="CDS position (nt)", y="Buse-R^2 * sign", title="GenomicCAI" )
## print(p)

## rr1 <- data.frame()
## tt11 <- glsRegressionRangeAnalysis( traits, tree, "LogGenomeSize", c(1,31), "DeltaLFE", profileId=1 )
## tt11$type <- "deltaLFE"
## rr1 <- rbind( rr1, tt11 )
## tt11 <- glsRegressionRangeAnalysis( traits, tree, "LogGenomeSize", c(1,31), "Native",   profileId=2 )
## tt11$type <- "native"
## rr1 <- rbind( rr1, tt11 )
## tt11 <- glsRegressionRangeAnalysis( traits, tree, "LogGenomeSize", c(1,31), "Shuffled", profileId=3 )
## tt11$type <- "shuffled"
## rr1 <- rbind( rr1, tt11 )

## p <- ggplot( rr1, aes(x=Var1, y=Buse.R2) ) +
##     geom_hline( yintercept=0, color="black" ) +
##     geom_line( aes( color=type, group=type ), size=3, alpha=0.9 ) +
##     labs( x="CDS position (nt)", y="Buse-R^2 * sign", title="LogGenomeSize" )
## print(p)

## rr1 <- data.frame()
## tt11 <- glsRegressionRangeAnalysis( traits, tree, "LogGrowthTime", c(1,31), "DeltaLFE", profileId=1 )
## tt11$type <- "deltaLFE"
## rr1 <- rbind( rr1, tt11 )
## tt11 <- glsRegressionRangeAnalysis( traits, tree, "LogGrowthTime", c(1,31), "Native",   profileId=2 )
## tt11$type <- "native"
## rr1 <- rbind( rr1, tt11 )
## tt11 <- glsRegressionRangeAnalysis( traits, tree, "LogGrowthTime", c(1,31), "Shuffled", profileId=3 )
## tt11$type <- "shuffled"
## rr1 <- rbind( rr1, tt11 )

## p <- ggplot( rr1, aes(x=Var1, y=Buse.R2) ) +
##     geom_hline( yintercept=0, color="black" ) +
##     geom_line( aes( color=type, group=type ), size=3, alpha=0.9 ) +
##     labs( x="CDS position (nt)", y="Buse-R^2 * sign", title="LogGrowthTime" )
## print(p)


## dev.off()
## quit()

#-----------------------------------------------------------------------------------------------------------------------------------


#--------------------
#glsRegressionRangeAnalysis( traits, tree, "OptimumTemp", c(1,pyramidLength))
#performGLSregression( traits, tree, "OptimumTemp", "Profile.11" )
#performGLSregression( traits, tree, "OptimumTemp", "GenomicGC" )

#glsRegressionRangeAnalysis( traits, tree, "TemperatureRange", c(1,pyramidLength))
#performGLSregression( traits, tree, "TemperatureRange", "Profile.11" )
#performGLSregression( traits, tree, "TemperatureRange", "GenomicGC" )


prepareFilteredTreeForGLS <- function( traitsx, tree, filterTrait, studyTrait )
{
    stopifnot( any( filterTrait == colnames(traitsx) ) )  # filterTrait not found
    stopifnot( any( studyTrait  == colnames(traitsx) ) )  # studyTrait not found

    # Filter by the requested trait, and also discard tree tips for species missing data
    speciesWithMissingData <- row.names(traitsx[ (!traitsx[,filterTrait] | is.na(traitsx[,studyTrait])), ])
    tree <- drop.tip( tree, speciesWithMissingData )   # Filter species that will prevent analysis from being performed from the tree
    
    if( is.null(tree) )
    {
        return( list(NA, NA) )
    }
    
    # Discard trait data for species missing from the tree
    treeSpecies <- tree$tip.label
    traitsx <- traitsx[treeSpecies,]   # Discard traits not found in the tree

    if( is.null(traitsx) || nrow(traitsx) < 3 )  # no sense in analyzing a tiny tree...
    {
        return( list(NA, NA) )
    }

    # Check tree <-> data-frame correlation appears valid
    x <- phylo4d(tree, tip.data=traitsx, rownamesAsLabels=TRUE)
    tipsOk <- hasTipData( x )
    stopifnot(tipsOk)
    
    return( list( traitsx, tree ) )
}




performGLSregressionWithFilter <- function( traits, tree, filterTrait, Xtrait, Ytrait, plotRegression=TRUE )
{
    stopifnot( any( filterTrait == colnames(traits)  ) )  # filterTrait not found
    stopifnot( any( Xtrait      == colnames(traits)  ) )  # studyTrait  not found

    print("performGLSregressionWithFilter")

    #print(sum(!is.na(traits$GenomicGC)))
    # Filter tree by trait
    orig.traits <- duplicate( traits)
    ret <- prepareFilteredTreeForGLS( traits, tree, filterTrait, Xtrait )
    traits.filtered <- ret[[1]]
    tree.filtered   <- ret[[2]]
    #print(sum(!is.na(traits.filtered$GenomicGC)))
    #print(sum(!is.na(orig.traits$GenomicGC)))

    if( is.null(tree.filtered) || any( class(tree.filtered) == "logical" ) )
    {
        return( NA )
    }

    #N <- nTips(tree.filtered)
    
    return( performGLSregression( traits.filtered, tree.filtered, Xtrait, Ytrait, plotRegression=plotRegression, caption=sprintf("%s effect on %s (%s)", Xtrait, profileMode, filterTrait, traits), orig.traits, extras="" ) )
}

glsRangeAnalysisWithFilter <- function( traits, tree, filterTrait, studyTrait, pyramidSpec, profileId=1, plotCaption="", extras="" )
{
    stopifnot( any( filterTrait==colnames(traits) ) )  # filterTrait not found
    stopifnot( any( studyTrait==colnames(traits)  ) )  # studyTrait not found

    # Filter tree by trait
    ret <- prepareFilteredTreeForGLS( traits, tree, filterTrait, studyTrait )
    if( isTRUE( is.na(ret[[1]]) ) ) { return( NA ) }
    traits.filtered <- ret[[1]]
    tree.filtered   <- ret[[2]]
    stopifnot( nTips(tree.filtered) == nrow(traits.filtered) ) # after filtering, the tree must still match the traits array
        
    if( is.null(tree.filtered) || nTips(tree.filtered) < minimalTaxonSize )
    {
        #return( NA )
        return( data.frame( Var1=integer(), Var2=integer(), Pvalue=double(), Buse.R2=double(), Logpval=double(), NumSpecies=integer(), DirectionsIndicator=character(), MIC=double(), MAS=double(), Pearson.r=double(), MIC.pvalue=double(), BP.pvalue=double() ) )
    }

    if( any(class(traits.filtered[,studyTrait])==c("factor")) && nlevels(as.factor(as.numeric(traits.filtered[,studyTrait]))) < 2 )
    {
        print("Only one level remaining...")
        #return( NA )
        return( data.frame( Var1=integer(), Var2=integer(), Pvalue=double(), Buse.R2=double(), Logpval=double(), NumSpecies=integer(), DirectionsIndicator=character(), MIC=double(), MAS=double(), Pearson.r=double(), MIC.pvalue=double(), BP.pvalue=double() ) )
    }

    return( glsRegressionRangeAnalysis( traits.filtered, tree.filtered, studyTrait, pyramidSpec, sprintf("%s\n(N=%d)", filterTrait, nrow(traits.filtered) ), profileId=profileId, addKdePlot=FALSE, extras=extras ) )
}

olsRangeAnalysisWithFilter <- function( traits, filterTrait, studyTrait, pyramidSpec, profileId=1, plotCaption="", extras="" )
{
    stopifnot( any( filterTrait==colnames(traits) ) )  # filterTrait not found
    stopifnot( any( studyTrait==colnames(traits)  ) )  # studyTrait not found

    # Filter species by trait
    traits.filtered <- traits[as.logical(traits[,filterTrait]),]
        
    if( nrow(traits.filtered) < minimalTaxonSize )
    {
        return( data.frame( Var1=integer(), Var2=integer(), Pvalue=double(), Buse.R2=double(), Logpval=double(), NumSpecies=integer(), DirectionsIndicator=character(), MIC=double(), MAS=double(), Pearson.r=double(), MIC.pvalue=double(), BP.pvalue=double() ) )
    }

    if( any(class(traits.filtered[,studyTrait])==c("factor")) && nlevels(as.factor(as.numeric(traits.filtered[,studyTrait]))) < 2 )
    {
        print("Only one level remaining...")
        #return( NA )
        return( data.frame( Var1=integer(), Var2=integer(), Pvalue=double(), Buse.R2=double(), Logpval=double(), NumSpecies=integer(), DirectionsIndicator=character(), MIC=double(), MAS=double(), Pearson.r=double(), MIC.pvalue=double(), BP.pvalue=double() ) )
    }
    
    return( olsRegressionRangeAnalysis( traits.filtered, studyTrait, pyramidSpec, sprintf("%s\n(N=%d)", filterTrait, nrow(traits.filtered) ), profileId=profileId, addKdePlot=FALSE, extras=extras ) )
}

testPhylosignal <- function( traits, tree, studyTrait, caption="" )
{
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
    #-----------------------
    traits$bm      <- bmTrait
    traits$bm.rand <- mixedTrait
    traits$random  <- uncorTrait

    studyTraitRanks <- paste(studyTrait, ".ranks", sep="")

    traits[,studyTraitRanks] <- rank( traits[,studyTrait], ties.method="average" )
    print("corr:")
    print(cor( traits[,studyTrait], traits[,studyTraitRanks], method="pearson"  ))
    print(cor( traits[,studyTrait], traits[,studyTraitRanks], method="spearman" ))

    
    speciesWithMissingData <- row.names(traits[is.na(traits[studyTrait]),])
    #print("Bad species:")
    #print(speciesWithMissingData)
    #print(length(speciesWithMissingData))
    #print("Before tree pruning:")
    #print(nTips(tree))
    tree <- drop.tip( tree, speciesWithMissingData )   # Filter species that will prevent analysis from being performed from the tree
    #print("After tree pruning:")
    #print(nTips(tree))

    treeSpecies <- tree$tip.label
    #print("Remaining species:")
    #print(treeSpecies)
    traits <- traits[treeSpecies,]   # Discard traits not found in the tree
    #print(nrow(traits))

    #print("Remaining bad species:")
    #print(row.names(traits[is.na(traits["GenomicGC"]),]))

    treeWithTraits <- phylo4d(tree, tip.data=traits, rownamesAsLabels=TRUE)
    stopifnot(hasTipData( treeWithTraits ))
                                        #checkPhylo4d( treeWithTraits )


    #lipa.GenomicGC  <- lipaMoran( treeWithTraits, trait=c("GenomicGC"),  alternative="two-sided", reps=lipaMoranReps )
    #lipa.Profile.1  <- lipaMoran( treeWithTraits, trait=c("Profile.1"),  alternative="two-sided", reps=lipaMoranReps )
    #lipa.Profile.5  <- lipaMoran( treeWithTraits, trait=c("Profile.5"),  alternative="two-sided", reps=lipaMoranReps )

    lipa.studyTrait      <- lipaMoran( treeWithTraits, trait=c(studyTrait),        alternative="two-sided", reps=lipaMoranReps )
    lipa.studyTraitRanks <- lipaMoran( treeWithTraits, trait=c(studyTraitRanks),   alternative="two-sided", reps=lipaMoranReps )

    lipa.bm              <- lipaMoran( treeWithTraits, trait=c("bm"),              alternative="two-sided", reps=lipaMoranReps )
    lipa.bm.rand         <- lipaMoran( treeWithTraits, trait=c("bm.rand"),         alternative="two-sided", reps=lipaMoranReps )
    lipa.random          <- lipaMoran( treeWithTraits, trait=c("random"),          alternative="two-sided", reps=lipaMoranReps )
    
    #traits$lipa_GenomicGC  <- lipa.GenomicGC$lipa
    #traits$lipa_Profile.1  <- lipa.Profile.1$lipa
    #traits$lipa_Profile.5  <- lipa.Profile.5$lipa
    #traits[,paste("lipa.", studyTrait)] <- lipa.studyTrait$lipa

    traits[,"lipa_studyTrait"]         <- lipa.studyTrait$lipa
    traits[,"lipa_studyTraitRanks"]    <- lipa.studyTraitRanks$lipa
    traits[,"lipa_bm"]                 <- lipa.bm$lipa
    traits[,"lipa_bm_rand"]            <- lipa.bm.rand$lipa
    traits[,"lipa_rand"]               <- lipa.random$lipa

     
    #p.vals = matrix( c( lipa.GenomicGC$p.value, lipa.Profile.15$p.value, lipa.Profile.15$p.value, lipa.Profile.15$p.value ), ncol=4 )
    #mat.col <- ifelse( p.vals < 0.001, "red", "grey35")


    #print(names(lipa.Profile.15))

    #range.min <- min( 0, lipa.GenomicGC$lipa, lipa.Profile.1$lipa, lipa.Profile.5$lipa, lipa.Profile.15$lipa )
    #range.max <- max( 0, lipa.GenomicGC$lipa, lipa.Profile.1$lipa, lipa.Profile.5$lipa, lipa.Profile.15$lipa )
    range.min <- min( 0, lipa.studyTrait$lipa, lipa.bm$lipa, lipa.bm.rand$lipa, lipa.random$lipa )
    range.max <- max( 0, lipa.studyTrait$lipa, lipa.bm$lipa, lipa.bm.rand$lipa, lipa.random$lipa )
    
    treeForPlotting <- phylo4d(tree, tip.data=traits, rownamesAsLabels=TRUE)
    stopifnot(hasTipData( treeForPlotting ))



    #barplot( treeForPlotting, trait=c("lipa_studyTrait", "lipa_bm", "lipa_bm_rand", "lipa_rand"), bar.lwd=c(1), tip.cex=0.15, data.xlim=c(range.min, range.max), center=FALSE, scale=FALSE, bar.col=mat.col, main=caption )
    barplot( treeForPlotting, trait=c("lipa_studyTrait", "lipa_studyTraitRanks", "lipa_bm", "lipa_bm_rand", "lipa_rand"), bar.lwd=c(1), tip.cex=0.15, data.xlim=c(range.min, range.max), center=FALSE, scale=FALSE, main=caption )
    
    #return( data.frame( studyTrait.lipa = lipa.studyTrait$lipa, studyTraitRanked.lipa = lipa.studyTraitRanks$lipa ) )
    corrPearson  <- cor( lipa.studyTrait$lipa, lipa.studyTraitRanks$lipa, method="pearson")
    corrSpearman <- cor( lipa.studyTrait$lipa, lipa.studyTraitRanks$lipa, method="spearman")
    print(sprintf("LIPA methods correlation: Pearson: %.3g Spearman: %.3g", corrPearson, corrSpearman ) )
    if( lipaUseRankedValues )
    {
        return( lipa.studyTraitRanks$lipa )
    }
    else
    {
        return( lipa.studyTrait$lipa )
    }
}


testPhylosignalWithFilter <- function( traits, tree, filterTrait='Member_all_1', studyTrait="GenomicGC" )
{
    stopifnot( any( filterTrait==colnames(traits) ) )  # filterTrait not found
    
    # Filter tree by trait
    ret <- prepareFilteredTreeForGLS( traits, tree, filterTrait, studyTrait )
    if( isTRUE( is.na(ret[[1]])) ) { return( NA ) }
    traits.filtered <- ret[[1]]
    tree.filtered   <- ret[[2]]
    stopifnot( nTips(tree.filtered) == nrow(traits.filtered) ) # after filtering, the tree must still match the traits array
        
    if( is.null(tree.filtered) || nTips(tree.filtered) < minimalTaxonSize )
    {
        return( NA )
    }

    return( testPhylosignal( traits.filtered, tree.filtered, studyTrait, sprintf("%s\n(N=%d)", filterTrait, nrow(traits.filtered) ) ) )
}

# These are the names the of taxon membership variables
#taxGroups <- c('Member_all_1', 'Member_Bacteria_2', 'Member_Terrabacteria_group_1783272', 'Member_Proteobacteria_1224', 'Member_Eukaryota_2759', 'Member_Archaea_2157', 'Member_Gammaproteobacteria_1236', 'Member_Firmicutes_1239', 'Member_FCB_group_1783270', 'Member_Euryarchaeota_28890', 'Member_Bacteroidetes_Chlorobi_group_68336', 'Member_Opisthokonta_33154', 'Member_Bacteroidetes_976', 'Member_Actinobacteria_201174', 'Member_Bacilli_91061', 'Member_Actinobacteria_1760', 'Member_Flavobacteriia_117743', 'Member_Flavobacteriales_200644', 'Member_Fungi_4751', 'Member_unclassified_Bacteria_2323', 'Member_Alphaproteobacteria_28211', 'Member_Dikarya_451864', 'Member_Bacteria_candidate_phyla_1783234', 'Member_Flavobacteriaceae_49546', 'Member_Patescibacteria_group_1783273', 'Member_Bacillales_1385', 'Member_TACK_group_1783275', 'Member_Parcubacteria_group_1794811', 'Member_Ascomycota_4890', 'Member_PVC_group_1783257', 'Member_Cyanobacteria_1117', 'Member_Cyanobacteria_Melainabacteria_group_1798711', 'Member_Chloroflexi_200795', 'Member_saccharomyceta_716545', 'Member_Rhizobiales_356', 'Member_Deinococci_188787', 'Member_Enterobacterales_91347', 'Member_Deinococcus_Thermus_1297', 'Member_delta_epsilon_subdivisions_68525', 'Member_Thermoprotei_183924', 'Member_Tenericutes_544448', 'Member_Crenarchaeota_28889', 'Member_Mollicutes_31969', 'Member_Thermotogae_188708', 'Member_Stramenopiles_33634', 'Member_Aquificae_200783', 'Member_Thermotogae_200918', 'Member_Aquificae_187857', 'Member_Enterobacteriaceae_543', 'Member_Thermococcales_2258', 'Member_Thermococcaceae_2259', 'Member_Micrococcales_85006', 'Member_Thermococci_183968', 'Member_Basidiomycota_5204', 'Member_Clostridia_186801', 'Member_Bacillaceae_186817', 'Member_Lactobacillales_186826' )

taxGroups <- colnames(taxidToKingdom)
taxGroups <- taxGroups[startsWith(taxGroups, "Member_")]

taxGroupToTaxId <- vapply(taxGroups, function (x) as.integer(tail(strsplit(x[[1]], '_')[[1]], n=1)), integer(1) )
# Example: taxGroupToTaxId['Member_Cyanobacteria_1117'] == 1117


#majorGroupsTree <- read.newick(taxonTreeFilename)

#performGLSregression <- function( traits, tree, Xtrait, Ytrait, plotRegression=TRUE, caption="", traits.full=NA, bm.gamma=1.0 )

                                        #testPhylosignal( traits, tree )



## for( gr in taxGroups )
## {
##     print(gr)
##     results <- findValueOutliersWithFilter( traits, gr, "Profile.1", "<0" )
##     regressionResultsByTaxGroup <- rbind( regressionResultsByTaxGroup, data.frame( ExplanatoryVar=c("GenomicGC"),  Range=c(i-1), MaxRangeStart=c(maxResult$Var1), MaxRangeEnd=c(maxResult$Var2), TaxGroup=c(taxGroupToTaxId[gr]), TaxGroupName=c(gr), EffectSize=c(maxResult$Buse.R2), Pvalue=c(maxResult$Pvalue), NumSpecies=c(maxResult$NumSpecies) ) )
    
##     results <- findValueOutliersWithFilter( traits, gr, "Profile.16", ">0" )
##     regressionResultsByTaxGroup <- rbind( regressionResultsByTaxGroup, data.frame( ExplanatoryVar=c("GenomicGC"),  Range=c(i-1), MaxRangeStart=c(maxResult$Var1), MaxRangeEnd=c(maxResult$Var2), TaxGroup=c(taxGroupToTaxId[gr]), TaxGroupName=c(gr), EffectSize=c(maxResult$Buse.R2), Pvalue=c(maxResult$Pvalue), NumSpecies=c(maxResult$NumSpecies) ) )

##     # Store peaks (for each range) for this group

## }


#dev.off()
#quit()

## out <- NA

## print(nTips(tree))

## for( i in 1:15 )
## {
##     ret <- testPhylosignalWithFilter( traits, tree, "Member_all_1", sprintf("Profile.%d", i) )
##     if( i > 1 )
##     {
##         out <- cbind( out, ret )
##     }
##     else
##     {
##         out <- ret
##     }
    
##     print(dim(out))
## }
## #print(out)

## write.csv(out, file=phylosignalLipaOutputFile )

#performGLSregression( traits, tree,   "GenomicGC",   "Profile.16", plotRegression=TRUE, caption="profile_150" )
##performGLSregression( traits, tree,   "GenomicGC",   "Profile_1.1" , plotRegression=TRUE, caption="profile_0"   )
#performGLSregression( traits, tree,   "GenomicGC",   "Profile_2.1" , plotRegression=TRUE, caption="profile_0"   )

## performGLSregression( traits, tree,   "Profile_1.1",   "Profile_2.1" , plotRegression=TRUE, caption="profile_0"   )

## performGLSregression( traits, tree,   "Profile_1.1",   "Profile_2.16" , plotRegression=TRUE, caption="profile_1.1x2.16"   )

#performGLSregression( traits, tree,   "Profile_1.14",   "Profile_1.15" , plotRegression=TRUE, caption="profile_1.14x1.15"   )

#performGLSregression( traits, tree,   "Profile_2.1",   "Profile_2.2" , plotRegression=TRUE, caption="profile_2.1x2.2"   )


## performGLSregressionWithFilter( traits, tree, "Member_all_1",                  "Profile_1.14",   "Profile_1.15" )

## performGLSregressionWithFilter( traits, tree, "Member_all_1",                  "Profile_2.1",   "Profile_2.2" )

## performGLSregressionWithFilter( traits, tree, "Member_all_1",                  "Profile_1.1",   "Profile_2.16" )

## performGLSregressionWithFilter( traits, tree, "Member_all_1",                  "Profile_1.14",   "Profile_2.1" )

## performGLSregressionWithFilter( traits, tree, "Member_all_1",                  "Profile_1.1",   "Profile_1.15" )
## performGLSregressionWithFilter( traits, tree, "Member_all_1",                  "Profile_2.1",   "Profile_2.16" )


#------------------------------------------------------------------------------------------------------------------------

calcMedianProfile <- function(traits, profileRange, profileId)
{
    #, getProfilePositions( profileId )
    vars = traits[, getProfileVariables( profileRange, profileId ) ]
    meds <- apply( vars, MARGIN=2, FUN=median )
    print(meds)
    data.frame(pos=getProfilePositions( profileId ), value=meds )
}

plotMedianProfile <- function( profileData, profileId )
{
    stopifnot( profileMode == "dLFE" )
    
    #ampl <- max(abs(min(profileData$value)), abs(max(profileData$value)))
    ampl <- 2.8436249
    profile.norm <- (profileData$value+ampl)/(2*ampl)
    print(profile.norm)
    stopifnot(all(profile.norm >= 0.0))
    stopifnot(all(profile.norm <= 1.0))
    print(profile.norm)
    steepness <- 15.0
    profile.norm.logistic <- 1/(1+exp(-steepness*(profile.norm-0.5)))
    profileData$profile.norm.logistic <- profile.norm.logistic

    print(profileData)
    
    p <- ggplot( profileData, aes(x=pos, y=y) ) +
#        geom_tile( aes(fill=profile.norm.logistic), linetype=0, colour=NA ) +
        geom_raster( aes(fill=profile.norm.logistic), linetype=0, colour=NA ) +
        scale_fill_gradientn(
            breaks=c(0.0, 0.5, 1.0),
            limits=c(0.0, 1.0),
            colors=        c(                                   "#ff0000", "#ff8080", "#ffffff", "#8080ff", "#0000ff"),
            values=rescale(c(                                         0.0,      0.35,       0.5,       0.65,      1.0) ) )
    if( profileId==2 )
    {
        p <- p + scale_x_reverse()
    }
    
    print(p)
}


## x1aa <- calcMedianProfile(traits, c(1,31), profileId=1)
## x1aa$y <- 1
## x2aa <- calcMedianProfile(traits, c(2,32), profileId=2)
## x2aa$y <- 1

## plotMedianProfile(x1aa, profileId=1)
## plotMedianProfile(x2aa, profileId=2)

#dev.off()
#quit()

#------------------------------------------------------------------------------------------------------------------------
# Measure the effects of GC vs. ENc' as predictors of dLFE
#------------------------------------------------------------------------------------------------------------------------

# get fast data - compatible with glsRegressionRangeAnalysis
testDataForRegressionAnalysis <- function( traits, tree, Xtrait, filterTrait, profileRange, plotCaption, profileId=1, addKdePlot=TRUE, extras=extras )
{
    # returns df:
    #    Var1 Var2       Pvalue      Buse.R2     Logpval NumSpecies DirectionsIndicator
    out <- data.frame( Var1=integer(), Var2=integer(), Pvalue=double(), Buse.R2=double(), Logpval=double(), NumSpecies=integer(), DirectionsIndicator=character() )

    phase <- runif(n=1)*15
    testVals <- sin( (seq(profileRange[1], profileRange[2], 1) + phase) * 0.40 )
    pvals <- (1-abs(testVals))*0.05
    positions <- getProfilePositions( profileId=profileId )
    for( i in seq(profileRange[1], profileRange[2], 1))
    {
        out <- rbind( out, data.frame( Var1=c(positions[i]), Var2=c(positions[i]), Pvalue=pvals[i], Buse.R2=testVals[i], Logpval=log(pvals[i]), NumSpecies=c(99), DirectionsIndicator=c("0") ) )
    }
    return( out );
}


make.group <- function(v)
{
    group <- 0;
    in.group <- FALSE;

    if( length(v)==0)
        return( c());

    out <- rep(0, length(v))

    for( i in 1:length(v) )
    {
        if( v[i] )
        {
            if( !in.group )
            {
                in.group <- TRUE
                group <- group + 1
            }
            out[i] <- group
        }
        else
        {
            in.group <- FALSE
            out[i] <- 0
        }
    }
    return(out)
}

## glsRangeAnalysisWithFilter( traits, tree, "Member_Eukaryota_2759", "GenomicGC", c(1, pyramidLength), profileId=1, plotCaption="", extras="")

## dev.off()
## quit()

## glsRegressionRangeAnalysis( traits, tree, "GenomicGC",  c(1,pyramidLength), plotCaption="GenomicGC", extras="")
## glsRegressionRangeAnalysis( traits, tree, "GenomicENc", c(1,pyramidLength), plotCaption="GenomicENc", extras="")
## glsRegressionRangeAnalysis( traits, tree, "GenomicENc.prime", c(1,pyramidLength), plotCaption="GenomicENc.prime", extras="")
## glsRegressionRangeAnalysis( traits, tree, "GenomicGC",  c(1,pyramidLength), plotCaption="GC+ENc'", extras=" GenomicENc.prime ")
#------------------------------------------------------------------------------------------------------------------------

partialDeterminationAnalysis <- function( traits, tree, trait, pyramidSpec, profileId=1, extras="" )
{
    #ps <- list(NA,NA)
    
    results.t1       <- glsRangeAnalysisWithFilter( traits, tree, "Member_all_1",          trait, pyramidSpec, plotCaption="", profileId=profileId, extras=extras)
    results.t2       <- glsRangeAnalysisWithFilter( traits, tree, "Member_Bacteria_2",     trait, pyramidSpec, plotCaption="", profileId=profileId, extras=extras)
    results.t3       <- glsRangeAnalysisWithFilter( traits, tree, "Member_Eukaryota_2759", trait, pyramidSpec, plotCaption="", profileId=profileId, extras=extras)
    results.t4       <- glsRangeAnalysisWithFilter( traits, tree, "Member_Fungi_4751",     trait, pyramidSpec, plotCaption="", profileId=profileId, extras="")
    results.t5       <- glsRangeAnalysisWithFilter( traits, tree, "Member_Archaea_2157",   trait, pyramidSpec, plotCaption="", profileId=profileId, extras=extras)
    results.all <- data.frame()
                                        #results.t1       <- glsRegressionRangeAnalysis( traits, tree, trait1, pyramidSpec, plotCaption="", profileId=profileId, extras="")
    #results.t1       <- testDataForRegressionAnalysis( traits, tree, trait1, pyramidSpec, plotCaption="", profileId=profileId, extras="")
    #results.t2       <- glsRegressionRangeAnalysis( traits, tree, trait2, pyramidSpec, plotCaption="", profileId=profileId, extras="")
    #results.combined <- glsRegressionRangeAnalysis( traits, tree, trait1, pyramidSpec, plotCaption=trait2, profileId=profileId, extras=trait2)
    # returns df:
    #    Var1 Var2       Pvalue      Buse.R2     Logpval NumSpecies DirectionsIndicator

    #print("pret1")
    if( !isTRUE(is.na(results.t1)) ) { if( nrow(results.t1) > 0 ) {
        #print("t1")
        ##with(results.t1, {
        ##    T <- "All"
        ##    significant <- Pvalue<significanceLevel
        ##    significant.group <- make.group(significant)
        ##})
        results.t1$T <- "All"
        
        results.t1$significant <- results.t1$Pvalue<significanceLevel
        results.t1$significant.group <- make.group( results.t1$significant )

        results.t1$MIC.significant <- results.t1$MIC.pvalue<significanceLevel
        results.t1$MIC.significant.group <- make.group( results.t1$MIC.significant )
        results.all <- rbind( results.all, results.t1)
    }}

    if( !isTRUE(is.na(results.t2)) ) { if( nrow(results.t2) > 0 ) {
        #print("t2")
        ##with(results.t2, {
        ##    T <- "Bacteria"
        ##    significant <- Pvalue<significanceLevel
        ##    significant.group <- make.group(significant)
        ##})
        results.t2$T <- "Bacteria"
        results.t2$significant <- results.t2$Pvalue<significanceLevel
        results.t2$significant.group <- make.group( results.t2$significant )

        results.t2$MIC.significant <- results.t2$MIC.pvalue<significanceLevel
        results.t2$MIC.significant.group <- make.group( results.t2$MIC.significant )
        results.all <- rbind( results.all, results.t2)
    }}

    if( !isTRUE(is.na(results.t3)) ) { if( nrow(results.t3) > 0 ) {
        results.t3$T <- "Eukaryota"
        results.t3$significant <- results.t3$Pvalue<significanceLevel
        results.t3$significant.group <- make.group( results.t3$significant )

        results.t3$MIC.significant <- results.t3$MIC.pvalue<significanceLevel
        results.t3$MIC.significant.group <- make.group( results.t3$MIC.significant )
        results.all <- rbind( results.all, results.t3)
    }}

    if( !isTRUE(is.na(results.t4)) ) { if( nrow(results.t4) > 0 ) {
        results.t4$T <- "Fungi"
        results.t4$significant <- results.t4$Pvalue<significanceLevel
        results.t4$significant.group <- make.group( results.t4$significant )

        results.t4$MIC.significant <- results.t4$MIC.pvalue<significanceLevel
        results.t4$MIC.significant.group <- make.group( results.t4$MIC.significant )        
        results.all <- rbind( results.all, results.t4)
    }}
    
    if( !isTRUE(is.na(results.t5)) ) { if( nrow(results.t5) > 0 ) {
        results.t5$T <- "Archaea"
        results.t5$significant <- results.t5$Pvalue<significanceLevel
        results.t5$significant.group <- make.group( results.t5$significant )

        results.t5$MIC.significant <- results.t5$MIC.pvalue<significanceLevel
        results.t5$MIC.significant.group <- make.group( results.t5$MIC.significant )        
        results.all <- rbind( results.all, results.t5)
    }}


    #results.all <- rbind( results.t1, results.t2, results.t3, results.t5 )
    #results.all <- rbind( results.t1, results.t2, results.t3, results.t4, results.t5 )

    return(results.all)
}    
        

plotPartialDetermination <-  function(results.all, traitName, profileId=1, yrange=NA, referenceLine=NA)
{
    if( isTRUE(is.na(yrange)) )
    {
        yrange <- c(-1.0, 1.0)
    }
    stopifnot(yrange[2] > yrange[1])
    stopifnot(length(yrange)==2)
    
    ## xlabels <- NA
    ## if( profileId==1 )
    ## {
    ##     xlabels <- c(0,100,200,300)
    ## }
    ## else
    ## {
    ##     xlabels <- c(-300,-200,-100,0)
    ## }



    ## print("---------1---------")
    ## p <- ggplot(results.t1, aes(x=Var1) ) + 
    ##     geom_area( data=results.t1[(results.t1$Buse.R2>0) & (results.t1$significant),], aes(x=Var1, y=Buse.R2, group=significant.group), position="identity", fill="#cc5533", color=NA,  alpha=0.5 ) +
    ##     geom_area( data=results.t1[(results.t1$Buse.R2<0) & (results.t1$significant),], aes(x=Var1, y=Buse.R2, group=significant.group), position="identity", fill="#3355cc", color=NA, alpha=0.5 ) +
    ##     geom_hline( yintercept=0 ) +
    ##     geom_line( aes(y=Buse.R2), alpha=0.8, size=1 ) +
    ##     geom_line( aes(y=Buse.R2, alpha=factor(significant), group=significant.group), size=0.3, color="white", linetype="dotted" ) +
    ##     scale_alpha_manual( values=c(0.0, 1.0) ) +
    ##     scale_y_continuous( limits=yrange ) +
    ##     scale_x_continuous( breaks=c(0,100,200,300),
    ##                         labels=xlabels, expand=expand_scale( mult=c(0,0) ) ) +
    ##     labs(title="R^2 comparison", x=sprintf("CDS position (relative to %s) (nt)", getProfileReference(profileId)), y="R^2") +
    ##     guides( alpha=FALSE )     # Hide the legend (will be plotted separately)
    ## ggsave(sprintf("R2comparison_%s.pdf", trait), plot=p, width=20, height=8, units="cm")

    p <- ggplot(results.all, aes(x=Var1) ) + 
        geom_area( data=results.all[(results.all$T=="All") & (results.all$Buse.R2>0) & (results.all$significant),], aes(x=Var1, y=Buse.R2, group=significant.group), position="identity", fill="#cc5533", color=NA,  alpha=0.5 ) +
        geom_area( data=results.all[(results.all$T=="All") & (results.all$Buse.R2<0) & (results.all$significant),], aes(x=Var1, y=Buse.R2, group=significant.group), position="identity", fill="#3355cc", color=NA, alpha=0.5 ) +
    geom_hline( yintercept=0 ) +
        geom_line( aes(y=Buse.R2, color=T, size=T), alpha=0.8 ) +
        geom_line( data=results.all[results.all$T=="All",], aes(x=Var1, y=Buse.R2, alpha=factor(significant), group=significant.group), size=0.35, color="white", linetype="dotted" ) +
        geom_line( data=results.all[results.all$T=="Bacteria",], aes(x=Var1, y=Buse.R2, alpha=factor(significant), group=significant.group), size=0.35, color="white", linetype="dotted" ) +
        geom_line( data=results.all[results.all$T=="Eukaryota",], aes(x=Var1, y=Buse.R2, alpha=factor(significant), group=significant.group), size=0.35, color="white", linetype="dotted" ) +
        geom_line( data=results.all[results.all$T=="Fungi",], aes(x=Var1, y=Buse.R2, alpha=factor(significant), group=significant.group), size=0.35, color="white", linetype="dotted" ) +
        geom_line( data=results.all[results.all$T=="Archaea",], aes(x=Var1, y=Buse.R2, alpha=factor(significant), group=significant.group), size=0.35, color="white", linetype="dotted" ) +
        scale_alpha_manual( values=c(0.0, 1.0) ) +
        scale_colour_manual( values=c("All"="black", "Bacteria"="blue", "Eukaryota"="darkgreen", "Archaea"="orange", "Fungi"="darkgreen") ) +
        scale_size_manual(   values=c("All"=1.8, "Bacteria"=0.8, "Eukaryota"=0.8, "Archaea"=0.8, "Fungi"=0.4) ) +
        scale_y_continuous( limits=yrange, breaks=c(-0.25, 0.0, 0.25, 0.5, 0.75, 1.00), labels=c("0.25", "0.00", "0.25", "0.50", "0.75", "1.00" ) ) +
        labs(title=traitName, x=sprintf("CDS position (relative to %s) (nt)", getProfileReference(profileId)), y="R^2") +
        theme( plot.background = element_blank(),   # Hide unnecessary theme elements (background panels, etc.)
              panel.grid.major.y = element_line(color="grey", size=0.50),
              panel.grid.major.x = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank()
              )
        guides( alpha=FALSE )     # Hide the legend (will be plotted separately)

    if( profileId==1 )
    {
        p <- p +
            scale_x_continuous( breaks=c(0,100,200,300),
                                labels=c(0,100,200,300), expand=expand_scale( mult=c(0,0) ) )
    }
    else
    {
        p <- p +
            scale_x_continuous( breaks=c(-300,-200,-100,0),
                                labels=c(-300,-200,-100,0), expand=expand_scale( mult=c(0,0) ) )
    }
    
    if( !is.na(referenceLine))
    {
        if( profileId==1 )
        {
            xpos <- 200
        }
        else
        {
            xpos <- -200
        }
        
        p <- p +
            geom_hline( yintercept=referenceLine, linetype="dashed", color="red", size=1.5, alpha=0.6  ) +
            annotate( "text", x=xpos, y=referenceLine, label=sprintf("Scale R^2 = %.2f", referenceLine), hjust=0, vjust=-1, color="red", alpha=0.6)
    }
    
    ## p1 <- ggplot(results.all, aes(x=Var1) ) + 
    ##     geom_hline( yintercept=0 ) +
    ##     geom_line( aes(y=MIC , color=T, size=T), alpha=0.8 ) +
    ##     geom_line( aes(y=BP.pvalue , color=T), size=0.4, alpha=1, linetype="dashed" ) +
    ##     geom_line( data=results.all[results.all$T=="All",],       aes(x=Var1, y=MIC, alpha=factor(MIC.significant), group=MIC.significant.group), size=0.35, color="white", linetype="dotted" ) +
    ##     geom_line( data=results.all[results.all$T=="Bacteria",],  aes(x=Var1, y=MIC, alpha=factor(MIC.significant), group=MIC.significant.group), size=0.35, color="white", linetype="dotted" ) +
    ##     geom_line( data=results.all[results.all$T=="Eukaryota",], aes(x=Var1, y=MIC, alpha=factor(MIC.significant), group=MIC.significant.group), size=0.35, color="white", linetype="dotted" ) +
    ##     geom_line( data=results.all[results.all$T=="Archaea",],   aes(x=Var1, y=MIC, alpha=factor(MIC.significant), group=MIC.significant.group), size=0.35, color="white", linetype="dotted" ) +
    ##     scale_alpha_manual( values=c(0.0, 1.0) ) +
    ##     scale_colour_manual( values=c("All"="black", "Bacteria"="blue", "Eukaryota"="darkgreen", "Archaea"="orange") ) +
    ##     scale_size_manual(   values=c("All"=1.6, "Bacteria"=1.0, "Eukaryota"=1.0, "Archaea"=1.0) ) +
    ##     scale_y_continuous( limits=c(0.0, 1.0) ) +
    ##     scale_x_continuous( breaks=c(0,100,200,300),
    ##                        labels=xlabels, expand=expand_scale( mult=c(0,0) ) ) +
    ##     labs(title=traitName, x=sprintf("CDS position (relative to %s) (nt)", getProfileReference(profileId)), y="MIC") +
    ##     theme( plot.background = element_blank(),   # Hide unnecessary theme elements (background panels, etc.)
    ##           panel.grid.major.y = element_line(color="grey", size=0.50),
    ##           panel.grid.major.x = element_blank(),
    ##           panel.grid.minor = element_blank(),
    ##           panel.background = element_blank()
    ##           ) +
    ##     guides( alpha=FALSE )     # Hide the legend (will be plotted separately)
     

    ## p2 <- ggplot(results.all, aes(x=Var1) ) + 
    ##     geom_hline( yintercept=0 ) +
    ##     geom_line( aes(y=MAS, color=T, size=T), alpha=1 ) +
    ##     scale_colour_manual( values=c("All"="black", "Bacteria"="blue", "Eukaryota"="darkgreen", "Archaea"="orange") ) +
    ##     scale_size_manual(   values=c("All"=1.6, "Bacteria"=1.0, "Eukaryota"=1.0, "Archaea"=1.0) ) +
    ##     scale_y_continuous( limits=c(0.0, 1.0) ) +
    ##     scale_x_continuous( breaks=c(0,100,200,300),
    ##                        labels=xlabels, expand=expand_scale( mult=c(0,0) ) ) +
    ##     labs(title=traitName, x=sprintf("CDS position (relative to %s) (nt)", getProfileReference(profileId)), y="MAS") +
    ##     guides( alpha=FALSE )     # Hide the legend (will be plotted separately)

    ## p3 <- ggplot(results.all, aes(x=Var1) ) + 
    ##     geom_hline( yintercept=0 ) +
    ##     geom_line( aes(y=pearson.r, color=T, size=T), alpha=1 ) +
    ##     scale_colour_manual( values=c("All"="black", "Bacteria"="blue", "Eukaryota"="darkgreen", "Archaea"="orange") ) +
    ##     scale_size_manual(   values=c("All"=1.6, "Bacteria"=1.0, "Eukaryota"=1.0, "Archaea"=1.0) ) +
    ##     scale_y_continuous( limits=c(-1.0, 1.0) ) +
    ##     scale_x_continuous( breaks=c(0,100,200,300),
    ##                        labels=xlabels, expand=expand_scale( mult=c(0,0) ) ) +
    ##     labs(title=traitName, x=sprintf("CDS position (relative to %s) (nt)", getProfileReference(profileId)), y="pearson.r") +
    ##     guides( alpha=FALSE )     # Hide the legend (will be plotted separately)
    
    if( profileId==1 ) # show dLFE window scale (=40nt)
    {
        p <- p +
            geom_line( data=data.frame(x=c(250,250+40),y=c(0.75,0.75)), aes(x=x,y=y), color="black", size=3) + # window scale
            annotate( "text", x=250+5, y=0.75, label="window", hjust=0, vjust=-1, color="black")
    }

    ## if( trait==trait1)
    ## {
    ##     ps[[1]] <- p
    ## }
    ## else if ( trait==trait2)
    ## {
    ##     ps[[2]] <- p
    ## }
    ## else
    ## {
    ##     stopifnot(FALSE)
    ## }
        
    # }
    
                                        #grid.draw( ggplotGrob( kdePlot ) )
    #grid.newpage()
                                        #pushViewport( viewport( x=unit(0.0, "npc"), y=unit(0.0, "npc"), width=unit(1.0, "npc"), height=unit(1.0, "npc") ) )

    #######kdePlot1 <- kdePlot( traits[, getProfileVariables( pyramidSpec, profileId ) ], getProfilePositions( profileId )  )

    #grid.arrange( ggplotGrob( kdePlot1 ), ps[[1]], ps[[2]], ncol=1, heights=c(unit(0.32, "npc"), unit(0.34, "npc"), unit(0.34, "npc") ) )

    grid.newpage()
    #grid.arrange( p, p1, p2, p3, ncol=1, heights=c(unit(0.34, "npc"), unit(0.22, "npc"), unit(0.22, "npc"), unit(0.22, "npc") ) )
    grid.arrange( p, p, p, p, ncol=1, heights=c(unit(0.34, "npc"), unit(0.22, "npc"), unit(0.22, "npc"), unit(0.22, "npc") ) )
    
    #pushViewport( viewport( x=unit(0.5, "npc"), y=unit(0.25, "npc") ) )
    #grid.draw(ggplotGrob( kdePlot ) )
    #print("--1--")
    #upViewport()
    #print("--2--")

    #pushViewport( viewport( x=unit(0.5, "npc"), y=unit(0.75, "npc") ) )
    #print("--3--")
    #grid.draw(p)
    #upViewport()
    
    ggsave(sprintf("R2comparison_%s_profile_%d_all.pdf", traitName, profileId), plot=p, width=20, height=8, units="cm")


    ## print("---------2---------")
    ## p <- ggplot(results.t2, aes(x=Var1) ) +
    ##     geom_hline( yintercept=0 ) +
    ##     geom_area( aes(y=Buse.R2), position="identity", color="grey", alpha=0.2 ) +
    ##     geom_line( aes(y=Buse.R2), alpha=0.8, size=1 ) +
    ##     #geom_line( aes(y=Buse.R2, alpha=factor(significant)), size=0.3, color="white", linetype="solid" ) +
    ##     #scale_alpha_manual( values=c(0.0, 1.0) ) +
    ##     labs(title="R^2 comparison", x="CDS position (relative to start) (nt)", y="R^2") 
    ## ggsave(sprintf("R2comparison_%s.pdf", trait2), plot=p, width=20, height=5, units="cm")
    
    ## print("---------3---------")
    ## p <- ggplot(results.all, aes(x=Var1) ) +
    ##     geom_hline( yintercept=0 ) +
    ##     geom_area( data=results.all[results.all$T=="both",], aes(x=Var1,y=abs(Buse.R2)), position="identity", color="grey", alpha=0.2 ) +
    ##     geom_line( aes(y=abs(Buse.R2), color=T), alpha=0.8, size=1 ) +
    ##     #geom_line( aes(y=abs(Buse.R2), alpha=factor(significant)), size=0.3, color="white", linetype="solid" ) +
    ##     #scale_alpha_manual( values=c(0.0, 1.0) ) +
    ##     labs(title="R^2 comparison", x="CDS position (relative to start) (nt)", y="R^2")
    ## ggsave(sprintf("R2comparison_%s_%s.pdf", trait1, trait2), plot=p, width=20, height=5, units="cm")
    ## print(p)
    
}


figure_PartialDeterminationAnalysis <- function(profileId=1, pyramidSpec=c(1,31) )
{
    d.all <- data.frame()
    yrange <- c(0.00,0.70)
    
    d.t1   <- glsRangeAnalysisWithFilter( traits, tree, "Member_all_1",          "GenomicGC",       pyramidSpec, plotCaption="GenomicGC", profileId=profileId, extras="" )
    d.t1$Buse.R2 <- abs(d.t1$Buse.R2)
    d.t1$T <- "GenomicGC"
    d.t1$significant <- d.t1$Pvalue<significanceLevel
    d.t1$significant.group <- make.group( d.t1$significant )
    d.t1$MIC.significant <- d.t1$MIC.pvalue<significanceLevel
    d.t1$MIC.significant.group <- make.group( d.t1$MIC.significant )
    d.all <- rbind( d.all, d.t1)
    
    d.t2   <- glsRangeAnalysisWithFilter( traits, tree, "Member_all_1",          "GenomicENc.prime", pyramidSpec, plotCaption="GenomicENc.prime", profileId=profileId, extras="" )
    d.t2$Buse.R2 <- abs(d.t2$Buse.R2)
    d.t2$T <- "GenomicENc.prime"
    d.t2$significant <- d.t2$Pvalue<significanceLevel
    d.t2$significant.group <- make.group( d.t2$significant )
    d.t2$MIC.significant <- d.t2$MIC.pvalue < significanceLevel
    d.t2$MIC.significant.group <- make.group( d.t2$MIC.significant )
    d.all <- rbind( d.all, d.t2)

    d.t3   <- glsRangeAnalysisWithFilter( traits, tree, "Member_all_1",          "GenomicGC",        pyramidSpec, plotCaption="GenomicGC -> GenomicENc.prime", profileId=profileId, extras=" GenomicENc.prime " )
    d.t3$Buse.R2 <- abs(d.t3$Buse.R2)
    d.t3$T <- "GenomicGC+GenomicENc.prime"
    d.t3$significant <- d.t3$Pvalue<significanceLevel
    d.t3$significant.group <- make.group( d.t3$significant )
    d.t3$MIC.significant <- d.t3$MIC.pvalue < significanceLevel
    d.t3$MIC.significant.group <- make.group( d.t3$MIC.significant )
    d.all <- rbind( d.all, d.t3)

    ## d.t4   <- glsRangeAnalysisWithFilter( traits, tree, "Member_all_1",          "GenomicENc.prime", pyramidSpec, plotCaption="GenomicENc.prime -> GenomicGC", profileId=profileId, extras="  GenomicGC  " )
    ## d.t4$Buse.R2 <- abs(d.t4$Buse.R2)
    ## d.t4$T <- "Test"
    ## d.t4$significant <- d.t4$Pvalue<significanceLevel
    ## d.t4$significant.group <- make.group( d.t4$significant )
    ## d.t4$MIC.significant <- d.t4$MIC.pvalue < significanceLevel
    ## d.t4$MIC.significant.group <- make.group( d.t4$MIC.significant )
    ## d.all <- rbind( d.all, d.t4)
    
    #d1 <- partialDeterminationAnalysis( traits,            tree, "GenomicGC", c(1,pyramidLength), profileId=1, extras="" )

    #plotPartialDetermination( d1, "GenomicGC", profileId=1, yrange=c(-0.40,1.00) )


    #d2 <- partialDeterminationAnalysis( traits,            tree, "GenomicGC", c(2,pyramidLength+1), profileId=2, extras="" )

    #plotPartialDetermination( d2, "GenomicGC", profileId=2, yrange=c(-0.40,1.00) )

    #dev.off()
    #quit()

    ##partialDeterminationAnalysis( traits,            tree, "GenomicENc.prime", "GenomicGC", c(1,pyramidLength), profileId=2, extras="")

    ##partialDeterminationAnalysis( traits,            tree, "GenomicGC", "GenomicENc.prime", c(1,pyramidLength), profileId=1, extras="", yrange=c(-0.3,0.6) )
    ##partialDeterminationAnalysis( traits,            tree, "GenomicGC", "GenomicENc.prime", c(1,pyramidLength), profileId=1, extras="", yrange=c(-0.40,1.00) )
    ##partialDeterminationAnalysis( traits,            tree, "GenomicGC", "GenomicENc.prime", c(1,pyramidLength), profileId=2, extras="", yrange=c(-0.40,1.00) )


    #partialDeterminationAnalysis( traits,            tree, "OptimumTemp", "GenomicNc", c(1,pyramidLength), profileId=2, extras="")


    ##partialDeterminationAnalysis( traits,            tree, "GenomicGC", "LogGrowthTime", c(1,pyramidLength), profileId=1, extras="", yrange=c(-0.40,1.00) )


    ## d1 <- partialDeterminationAnalysis( traits,            tree, "LogGrowthTime", c(1,pyramidLength), profileId=1, extras="" )
    ## print(d1)
    ## plotPartialDetermination( d1, "LogGrowthTime", profileId=1, yrange=c(-1.0,1.0) )

    ## #dev.off()
    ## #quit()

    ##d1 <- partialDeterminationAnalysis( traits,            tree, "GenomicENc.prime", c(1,pyramidLength), profileId=1, extras="" )
    ##print(d1)
    ##plotPartialDetermination( d1, "GenomicENc.prime", profileId=1, yrange=c(-0.5,0.5) )

    ##d1 <- partialDeterminationAnalysis( traits,            tree, "GenomicENc.prime", c(2,pyramidLength+1), profileId=2, extras="" )
    ##print(d1)
    ##plotPartialDetermination( d1, "GenomicENc.prime", profileId=2, yrange=c(-0.5,0.5) )
    
    ## d1 <- partialDeterminationAnalysis( traits,            tree, "LogGrowthTime", c(1,pyramidLength), profileId=1, extras="GenomicGC" )
    ## print(d1)
    ## plotPartialDetermination( d1, "LogGrowthTime", profileId=1, yrange=c(-1.0,1.0) )

    p <- ggplot(d.all, aes(x=Var1) ) + 
        geom_line( aes(y=Buse.R2, color=T, size=T), alpha=0.7 ) +
        geom_line( data=d.all[d.all$T=="GenomicGC",],                  aes(x=Var1, y=Buse.R2, alpha=factor(significant), group=significant.group), size=0.35, color="white", linetype="dotted" ) +
        geom_line( data=d.all[d.all$T=="GenomicENc.prime",],           aes(x=Var1, y=Buse.R2, alpha=factor(significant), group=significant.group), size=0.35, color="white", linetype="dotted" ) +
        geom_line( data=d.all[d.all$T=="GenomicGC+GenomicENc.prime",], aes(x=Var1, y=Buse.R2, alpha=factor(significant), group=significant.group), size=0.35, color="white", linetype="dotted" ) +
        geom_line( data=d.all[d.all$T=="Test",],                       aes(x=Var1, y=Buse.R2, alpha=factor(significant), group=significant.group), size=0.35, color="white", linetype="dotted" ) +
        scale_alpha_manual( values=c(0.0, 1.0) ) +
        scale_colour_manual( values=c("GenomicGC"="blue", "GenomicENc.prime"="#208020", "GenomicGC+GenomicENc.prime"="red", "Test"="grey" ) ) +
        scale_size_manual(   values=c("GenomicGC"=1.0,    "GenomicENc.prime"=1.0,   "GenomicGC+GenomicENc.prime"=1.6,     "Test"=0.6    ) ) +
        scale_y_continuous( limits=yrange, breaks=c( 0.0, 0.20, 0.40, 0.60 ), labels=c( "0.0", "0.2", "0.4", "0.6" ) ) +
        labs(x=sprintf("CDS position (relative to %s) (nt)", getProfileReference(profileId)), y="R^2") +
        theme( plot.background = element_blank(),   # Hide unnecessary theme elements (background panels, etc.)
              panel.grid.major.y = element_line(color="grey", size=0.50),
              panel.grid.major.x = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              aspect.ratio = 0.34
              )
        guides( alpha=FALSE )     # Hide the legend (will be plotted separately)

    if( profileId==1 )
    {
        p <- p +
            scale_x_continuous( breaks=c(0,100,200,300),
                                labels=c(0,100,200,300), expand=expand_scale( mult=c(0,0) ) )
    }
    else
    {
        p <- p +
            scale_x_continuous( breaks=c(-300,-200,-100,0),
                                labels=c(-300,-200,-100,0), expand=expand_scale( mult=c(0,0) ) )
    }
    print(p)
    
}

figure_PartialDeterminationAnalysis2 <- function()
{
    
    ## gr1 <- "Member_all_1"
    ## #gr1 <- "Member_Bacteria_2"
    ## #gr1 <- "Member_Proteobacteria_1224"
    ## #gr1 <- "Member_Terrabacteria_group_1783272"

    ## rt1       <- glsRangeAnalysisWithFilter( traits, tree, gr1, "GenomicGC", c(1,31), profileId=1, extras="")
    ## rt1$T <- "Archaea"
    ## rt1$significant <- rt1$Pvalue<significanceLevel
    ## rt1$significant.group <- make.group( rt1$significant )

    ## rt1$MIC.significant <- rt1$MIC.pvalue<significanceLevel
    ## rt1$MIC.significant.group <- make.group( rt1$MIC.significant )
    ## rt1$Buse.R2 <- abs(rt1$Buse.R2)

    ## rt2       <- glsRangeAnalysisWithFilter( traits, tree, gr1, "LogGrowthTime", c(1,31), profileId=1, extras="")
    ## rt2$T <- "Bacteria"
    ## rt2$significant <- rt2$Pvalue<significanceLevel
    ## rt2$significant.group <- make.group( rt2$significant )

    ## rt2$MIC.significant <- rt2$MIC.pvalue<significanceLevel
    ## rt2$MIC.significant.group <- make.group( rt2$MIC.significant )
    ## rt2$Buse.R2 <- abs(rt2$Buse.R2)

    ## rt3       <- glsRangeAnalysisWithFilter( traits, tree, gr1, "LogGrowthTime", c(1,31), profileId=1, extras="GenomicGC")
    ## rt3$T <- "All"
    ## rt3$significant <- rt3$Pvalue<significanceLevel
    ## rt3$significant.group <- make.group( rt3$significant )

    ## rt3$MIC.significant <- rt3$MIC.pvalue<significanceLevel
    ## rt3$MIC.significant.group <- make.group( rt3$MIC.significant )
    ## rt3$Buse.R2 <- abs(rt3$Buse.R2)

    ## results.all <- rbind( rt1, rt2, rt3 )
    ## plotPartialDetermination( results.all, "LogGrowthTime", profileId=1, yrange=c(-1.0,1.0) )

    ## #-----
    ## rt4       <- glsRangeAnalysisWithFilter( traits, tree, gr1, "GenomicGC", c(1,31), profileId=1, extras="LogGrowthTime")
    ## rt4$T <- "Bacteria"
    ## rt4$significant <- rt4$Pvalue<significanceLevel
    ## rt4$significant.group <- make.group( rt4$significant )

    ## rt4$MIC.significant <- rt4$MIC.pvalue<significanceLevel
    ## rt4$MIC.significant.group <- make.group( rt4$MIC.significant )
    ## rt4$Buse.R2 <- abs(rt4$Buse.R2)

    ## results.all <- rbind( rt1, rt3, rt4 )
    ## plotPartialDetermination( results.all, "LogGrowthTime", profileId=1, yrange=c(-1.0,1.0) )


    ## #-----
    ## rt2       <- glsRangeAnalysisWithFilter( traits, tree, gr1, "LogGrowthTime", c(1,31), profileId=1, extras="")
    ## rt2$T <- "Bacteria"
    ## rt2$significant <- rt2$Pvalue<significanceLevel
    ## rt2$significant.group <- make.group( rt2$significant )

    ## rt2$MIC.significant <- rt2$MIC.pvalue<significanceLevel
    ## rt2$MIC.significant.group <- make.group( rt2$MIC.significant )
    ## #rt2$Buse.R2 <- abs(rt2$Buse.R2)

    ## results.all <- rbind( rt1, rt2, rt3 )
    ## plotPartialDetermination( results.all, "LogGrowthTime", profileId=1, yrange=c(-1.0,1.0) )
}


figure_PartialDeterminationAnalysis_GC_and_ENc.prime <- function()
{
    stopifnot( profileMode == "dLFE" )
    
    yrange=c(-0.45, 1.00)
    
    d1 <- partialDeterminationAnalysis( traits,            tree, "GenomicGC", c(1,pyramidLength), profileId=1, extras="" )
    #print(d1)
    plotPartialDetermination( d1, "GenomicGC", profileId=1, yrange=yrange )
    
    d1 <- partialDeterminationAnalysis( traits,            tree, "GenomicGC", c(2,pyramidLength+1), profileId=2, extras="" )
    #print(d1)
    #d1$Var1 <- d1$Var1-10
    plotPartialDetermination( d1, "GenomicGC", profileId=2, yrange=yrange )

    #print(traits$Member_Bacteria_2)
    #print(nrow(traits))
    
    yrange=c(-0.50, 0.50)
    
    d1 <- partialDeterminationAnalysis( traits,            tree, "GenomicENc.prime", c(1,pyramidLength), profileId=1, extras="" )
    #print(d1)
    plotPartialDetermination( d1, "GenomicENc.prime", profileId=1, yrange=yrange )

    d1 <- partialDeterminationAnalysis( traits,            tree, "GenomicENc.prime", c(2,pyramidLength+1), profileId=2, extras="" )
    #print(d1)
    plotPartialDetermination( d1, "GenomicENc.prime", profileId=2, yrange=yrange )

    
}


figure_PartialDeterminationAnalysis_NormalizedProfiles <- function()
{
    stopifnot( profileMode == "dLFE" )
    
    ## dres1 <- partialDeterminationAnalysis( traits.normalized,            tree, "GenomicGC", "GenomicENc.prime", c(1,pyramidLength), profileId=1, extras="", yrange=c(-0.40,1.00) )
    ## partialDeterminationAnalysis( dres1 )
    ##partialDeterminationAnalysis( traits.normalized,            tree, "GenomicGC", "GenomicENc.prime", c(1,pyramidLength), profileId=2, extras="", yrange=c(-0.40,1.00) )

    GLS.dLFE.sd.vs.GenomicGC        <- performGLSregression( traits.normalized, tree, "GenomicGC",                "dLFE.sd.12",  extras="", plotRegression=FALSE )
    GLS.dLFE.sd.vs.GenomicENc.prime <- performGLSregression( traits.normalized, tree, "GenomicENc.prime",         "dLFE.sd.12",  extras="", plotRegression=FALSE )

    yrange=c(-0.45, 0.70)
    
    d1 <- partialDeterminationAnalysis( traits.normalized,            tree, "GenomicGC", c(1,pyramidLength), profileId=1, extras="" )
    #print(d1)
    plotPartialDetermination( d1, "GenomicGC", profileId=1, yrange=yrange, referenceLine=GLS.dLFE.sd.vs.GenomicGC$R2 )
    
    d1 <- partialDeterminationAnalysis( traits.normalized,            tree, "GenomicGC", c(2,pyramidLength+1), profileId=2, extras="" )
    #print(d1)
    #d1$Var1 <- d1$Var1-10
    plotPartialDetermination( d1, "GenomicGC", profileId=2, yrange=yrange, referenceLine=GLS.dLFE.sd.vs.GenomicGC$R2 )

    #print(traits$Member_Bacteria_2)
    #print(nrow(traits))
    
    d1 <- partialDeterminationAnalysis( traits.normalized,            tree, "GenomicENc.prime", c(1,pyramidLength), profileId=1, extras="" )
    #print(d1)
    plotPartialDetermination( d1, "GenomicENc.prime", profileId=1, yrange=yrange, referenceLine=GLS.dLFE.sd.vs.GenomicENc.prime$R2 )

    d1 <- partialDeterminationAnalysis( traits.normalized,            tree, "GenomicENc.prime", c(2,pyramidLength+1), profileId=2, extras="" )
    #print(d1)
    plotPartialDetermination( d1, "GenomicENc.prime", profileId=2, yrange=yrange, referenceLine=GLS.dLFE.sd.vs.GenomicENc.prime$R2 )

    
}


## dev.off()
## quit()


## #--------------


## rt1       <- glsRangeAnalysisWithFilter( traits, tree, gr1, "GenomicGC", c(1,31), profileId=1, extras="")
## rt1$T <- "Archaea"
## rt1$significant <- rt1$Pvalue<significanceLevel
## rt1$significant.group <- make.group( rt1$significant )

## rt1$MIC.significant <- rt1$MIC.pvalue<significanceLevel
## rt1$MIC.significant.group <- make.group( rt1$MIC.significant )
## rt1$Buse.R2 <- abs(rt1$Buse.R2)

## for( i in 1:20) # noisy
## {
##     traits$Noise30 <- make.partial.noise.trait( nrow(traits) )
    
##     rt2       <- glsRangeAnalysisWithFilter( traits, tree, gr1, "Noise30", c(1,31), profileId=1, extras="")
##     rt2$T <- "Bacteria"
##     rt2$significant <- rt2$Pvalue<significanceLevel
##     rt2$significant.group <- make.group( rt2$significant )

##     rt2$MIC.significant <- rt2$MIC.pvalue<significanceLevel
##     rt2$MIC.significant.group <- make.group( rt2$MIC.significant )
##     rt2$Buse.R2 <- abs(rt2$Buse.R2)

##     rt3       <- glsRangeAnalysisWithFilter( traits, tree, gr1, "Noise30", c(1,31), profileId=1, extras="GenomicGC")
##     rt3$T <- "All"
##     rt3$significant <- rt3$Pvalue<significanceLevel
##     rt3$significant.group <- make.group( rt3$significant )

##     rt3$MIC.significant <- rt3$MIC.pvalue<significanceLevel
##     rt3$MIC.significant.group <- make.group( rt3$MIC.significant )
##     rt3$Buse.R2 <- abs(rt3$Buse.R2)

##     results.all <- rbind( rt1, rt2, rt3 )
##     plotPartialDetermination( results.all, "LogGrowthTime", profileId=1, yrange=c(-1.0,1.0) )
## }


## dev.off()
## quit()



## #--------------
## rt1       <- olsRangeAnalysisWithFilter( traits, gr1, "LogGrowthTime", c(1,31), profileId=1, extras="")
## rt1$T <- "Archaea"
## rt1$significant <- rt1$Pvalue<significanceLevel
## rt1$significant.group <- make.group( rt1$significant )

## rt1$MIC.significant <- rt1$MIC.pvalue<significanceLevel
## rt1$MIC.significant.group <- make.group( rt1$MIC.significant )

## rt2       <- glsRangeAnalysisWithFilter( traits, tree, gr1, "LogGrowthTime", c(1,31), profileId=1, extras="")
## rt2$T <- "Bacteria"
## rt2$significant <- rt2$Pvalue<significanceLevel
## rt2$significant.group <- make.group( rt2$significant )

## rt2$MIC.significant <- rt2$MIC.pvalue<significanceLevel
## rt2$MIC.significant.group <- make.group( rt2$MIC.significant )

## results.all <- rbind( rt1, rt2 )
## plotPartialDetermination( results.all, "LogGrowthTime", profileId=1, yrange=c(-1.0,1.0) )



## dev.off()
## quit()

#------------------------------------------------------------------------------------------------------------------------



correlationBetweenRanges <- function( traits, tree, filterTrait, xRange, xCoor, yRange, yCoor, useUncorrectedCorrelation=FALSE, marks=NA, halfOnly=0 )
{
    # 1=top 0=both -1=bottom
    stopifnot( length(xRange)==length(xCoor) )
    stopifnot( length(yRange)==length(yCoor) )
    stopifnot( (halfOnly <= 1) && (halfOnly >= -1) )
    corValues <- data.frame( Var1=character(), Var2=character(), Cor=double(), DirectionsIndicator=character(), Pval=double() )
    
    for( xi in 1:length(xRange) )
    {
        xTrait <- xRange[xi]
        xPos   <- xCoor[xi]
        
        for( yi in 1:length(yRange) )
        {
            yTrait <- yRange[yi]
            yPos   <- yCoor[yi]


            skip <- FALSE
            if( xPos>yPos )
            {
                skip <- halfOnly > 0
            }
            else if( xPos<yPos )
            {
                skip <- halfOnly < 0
            }
            
            
            if( xTrait != yTrait && (!skip) )
            {
                if( !useUncorrectedCorrelation )
                {
                    result <- performGLSregressionWithFilter( traits, tree, filterTrait, xTrait, yTrait, plotRegression=FALSE )
                    method <- "Method: Corrected (GLS)"
                }
                else
                {
                    result <- findUncorrectedVariableCorrelation( traits, filterTrait, xTrait, yTrait )
                    result$R2 <- result$corr  # compatibility hack...
                    #result$pvalue = result$pvalue
                    method <- "Method: Uncorrected (Spearman)"
                }
            }
            else
            {
                result <- list(R2=1.0, directionsIndicator='+', pvalue=0.0)
            }
            print(result)
            
            corValues <- rbind( corValues, data.frame( Var1=xPos, Var2=yPos, Cor=result$R2, DirectionsIndicator=result$directionsIndicator, Pval=result$pvalue ) )
        }
    }
    print(corValues)

    plotTitle <- sprintf("%s (N=%d)", filterTrait, result$N )

    p <- ggplot( corValues, aes(x=Var1, y=Var2 ) ) + #, fill=Cor) ) +
        geom_tile( aes(fill=Cor), linetype=0, colour=NA ) +   # Use geom_tile (because geom_raster has poor pdf/svg output), with minimal borders
        labs( title=plotTitle ) +
        scale_fill_gradientn(
            breaks=c(-1.0, -0.5, 0.0, 0.5, 1.0),
            limits=c(-1.0, 1.0),
            colors=        c("#75c0ff", "#2240cc", "#112060",   "#000000", "#602011", "#cc4022", "#ffc075"),
            values=rescale(c(     -1.0,      -0.5,     -0.15,         0.0,      0.15,       0.5,       1.0) ) ) +
        theme( plot.background = element_blank(),   # Hide unnecessary theme elements (background panels, etc.)
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank()
              )

    dotsDF <- data.frame(melt( corValues, id.vars=c("Var1", "Var2"), measure.vars=c("Pval") ))
    dotsDF <- dotsDF[dotsDF$value<=0.01,]
    
    p <- p +
        geom_point( data=dotsDF, aes(x=Var1,y=Var2), position=position_nudge(y=4, x=-4), color="white", size=((fontScale/12)**2)*0.6 )  # Significance markers (white dots)
    
    
    if( !all(is.na(marks)) )
    {
        marksDf <- data.frame(t(matrix(data=marks, nrow=2, ncol=length(marks)/2)))
        marksDf$label <- letters[1:(length(marks)/2)];
        print(marksDf)
       
        # Overlay the "points of interest" on the plot
        p <- p + 
            geom_text(  data=marksDf, aes(x=X1, y=X2, label=label), color="black", size=6, position=position_nudge(y=0.5, x=-0.5) ) + 
            geom_text(  data=marksDf, aes(x=X1, y=X2, label=label), color="white", size=6 )
        #    geom_point( data=marksDf, aes(x=X1, y=X2), position=position_nudge(y=4, x=-4), color="yellow" ) +

        marksDf <- rename(marksDf, Var1=X1, Var2=X2)
        tblDf <- marksDf %>% left_join(corValues) %>% select(-DirectionsIndicator)
        tblDf$Cor.fmt  <- format(round(tblDf$Cor,  digit=3))
        tblDf$Pval.fmt <- format(round(tblDf$Pval, digit=6))
        print(tblDf)
        
        grid.newpage()
        tbl <- tableGrob( tblDf %>% select(-Cor, -Pval), cols=c("Pos1 (nt)", "Pos2 (nt)", "Label", "Correlation", "p-val") )  # , cols=c("label","X1","X2","Cor") )
        grid.draw(tbl)
        
    }
    
    print(p)

}




##performGLSregressionWithFilter( traits, tree, "Member_all_1",      "GenomicGC",   "Profile_1.1" )
##performGLSregressionWithFilter( traits, tree, "Member_all_1",      "GenomicGC",   "Profile_1.26" )
##performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2", "GenomicGC",   "Profile_1.26" )

##dev.off()
##quit()

figure_CorrelationBetweenRanges <- function()
{

    ## Autocorrelation plots, showing the correlation between different points along the profile

    ## correlationBetweenRanges( traits, tree, "Member_all_1",
    ##                          sapply(seq(1,31,1), function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          seq(0,300,10),
    ##                          sapply(seq(1,31,1), function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          seq(0,300,10),
    ##                          useUncorrectedCorrelation=TRUE,
    ##                          marks=c(10,  0, 200, 0, 300, 0, 300, 100, 0, 10, 0, 200, 0, 300, 100, 300),
    ##                          halfOnly=1 )

    ## correlationBetweenRanges( traits.normalized, tree, "Member_all_1",
    ##                          sapply(seq(1,31,1), function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          seq(0,300,10),
    ##                          sapply(seq(1,31,1), function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          seq(0,300,10),
    ##                          useUncorrectedCorrelation=TRUE,
    ##                          marks=c(10,  0, 200, 0, 300, 0, 300, 100, 0, 10, 0, 200, 0, 300, 100, 300),
    ##                          halfOnly=1 )


    ## correlationBetweenRanges( traits, tree, "Member_all_1",
    ##                          sapply(seq(1,31,1), function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          seq(0,300,10),
    ##                          sapply(seq(1,31,1), function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          seq(0,300,10),
    ##                          useUncorrectedCorrelation=TRUE,
    ##                          marks=c(10,  0, 200, 0, 300, 0, 300, 100, 0, 10, 0, 200, 0, 300, 100, 300),
    ##                          halfOnly=1 )


    ## correlationBetweenRanges( traits, tree, "Member_Bacteria_2", 
    ##                          sapply(seq(1,31,1), function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          seq(0,300,10),
    ##                          sapply(seq(1,31,1), function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          seq(0,300,10),
    ##                          useUncorrectedCorrelation=TRUE,
    ##                          marks=c(10,  0, 200, 0, 300, 0, 300, 100, 0, 10, 0, 200, 0, 300, 100, 300),
    ##                          halfOnly=1 )

    ## correlationBetweenRanges( traits, tree, "Member_Eukaryota_2759", 
    ##                          sapply(seq(1,31,1), function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          seq(0,300,10),
    ##                          sapply(seq(1,31,1), function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          seq(0,300,10),
    ##                          useUncorrectedCorrelation=TRUE,
    ##                          marks=c(10,  0, 200, 0, 300, 0, 300, 100, 0, 10, 0, 200, 0, 300, 100, 300),
    ##                          halfOnly=1 )


    ## correlationBetweenRanges( traits, tree, "Member_Archaea_2157", 
    ##                          sapply(seq(1,31,1), function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          seq(0,300,10),
    ##                          sapply(seq(1,31,1), function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          seq(0,300,10),
    ##                          useUncorrectedCorrelation=TRUE,
    ##                          marks=c(10,  0, 200, 0, 300, 0, 300, 100, 0, 10, 0, 200, 0, 300, 100, 300),
    ##                          halfOnly=1 )


    ## correlationBetweenRanges( traits, tree, "Member_Proteobacteria_1224", 
    ##                          sapply(seq(1,31,1), function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          seq(0,300,10),
    ##                          sapply(seq(1,31,1), function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          seq(0,300,10),
    ##                          useUncorrectedCorrelation=TRUE,
    ##                          marks=c(10,  0, 200, 0, 300, 0, 300, 100, 0, 10, 0, 200, 0, 300, 100, 300),
    ##                          halfOnly=1 )

    ## correlationBetweenRanges( traits, tree, "Member_Proteobacteria_1224", 
    ##                          sapply(seq(1,31,1), function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          seq(0,300,10),
    ##                          sapply(seq(1,31,1), function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          seq(0,300,10),
    ##                          useUncorrectedCorrelation=FALSE,
    ##                          marks=c(10,  0, 200, 0, 300, 0, 300, 100, 0, 10, 0, 200, 0, 300, 100, 300),
    ##                          halfOnly=-1 )



    ## correlationBetweenRanges( traits, tree, "Member_Opisthokonta_33154",
    ##                          sapply(seq(1,31,1), function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          seq(0,300,10),
    ##                          sapply(seq(1,31,1), function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          seq(0,300,10),
    ##                          useUncorrectedCorrelation=TRUE,
    ##                          marks=c(10,  0, 200, 0, 300, 0, 300, 100, 0, 10, 0, 200, 0, 300, 100, 300),
    ##                          halfOnly=1 )


    ## correlationBetweenRanges( traits, tree, "Member_Opisthokonta_33154",
    ##                          sapply(seq(1,31,1), function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          seq(0,300,10),
    ##                          sapply(seq(1,31,1), function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          seq(0,300,10),
    ##                          useUncorrectedCorrelation=FALSE,
    ##                          marks=c(10,  0, 200, 0, 300, 0, 300, 100, 0, 10, 0, 200, 0, 300, 100, 300),
    ##                          halfOnly=-1 )

    ## dev.off()
    ## quit()




    ## correlationBetweenRanges( traits, tree, "Member_all_1",
    ##                          sapply(seq(1,31,1), function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          seq(0,300,10),
    ##                          sapply(seq(1,31,1), function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          seq(0,300,10),
    ##                          useUncorrectedCorrelation=FALSE,
    ##                          marks=c(10,  0, 200, 0, 300, 0, 300, 100),
    ##                          halfOnly=-1)

    ## dev.off()
    ## quit()


    ## correlationBetweenRanges( traits, tree, "Member_all_1",
    ##                          sapply(seq(32,2,-1), function (j) { sprintf("Profile_2.%d", j) } ),
    ##                          seq(300,0,-10),
    ##                          sapply(seq(32,2,-1), function (j) { sprintf("Profile_2.%d", j) } ),
    ##                          seq(300,0,-10),
    ##                          useUncorrectedCorrelation=TRUE )


    ## correlationBetweenRanges( traits, tree, "Member_all_1",
    ##                          sapply(1:31, function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          seq(0,300,10),
    ##                          sapply(2:32, function (j) { sprintf("Profile_2.%d", j) } ),
    ##                          seq(300,0,-10),
    ##                          useUncorrectedCorrelation=TRUE )

    ## correlationBetweenRanges( traits, tree, "Member_Eukaryota_2759",
    ##                          sapply(1:31, function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          sapply(1:31, function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          useUncorrectedCorrelation=TRUE )


    ## correlationBetweenRanges( traits, tree, "Member_Eukaryota_2759",
    ##                          sapply(1:32, function (j) { sprintf("Profile_2.%d", j) } ),
    ##                          sapply(1:32, function (j) { sprintf("Profile_2.%d", j) } ),
    ##                          useUncorrectedCorrelation=TRUE )


    ## correlationBetweenRanges( traits, tree, "Member_Eukaryota_2759",
    ##                          sapply(1:31, function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          sapply(1:32, function (j) { sprintf("Profile_2.%d", j) } ),
    ##                          useUncorrectedCorrelation=TRUE )

    ## correlationBetweenRanges( traits, tree, "Member_all_1",
    ##                          sapply(1:31, function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          seq(0,300,10),
    ##                          sapply(1:31, function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          seq(0,300,10),
    ##                          useUncorrectedCorrelation=FALSE )


    ## correlationBetweenRanges( traits, tree, "Member_Eukaryota_2759",
    ##                          sapply(seq(1,31,1), function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          sapply(seq(1,31,1), function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          useUncorrectedCorrelation=FALSE )

    ## correlationBetweenRanges( traits, tree, "Member_Bacteria_2",
    ##                          sapply(seq(1,31,1), function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          sapply(seq(1,31,1), function (j) { sprintf("Profile_1.%d", j) } ),
    ##                          useUncorrectedCorrelation=FALSE )

    ## correlationBetweenRanges( traits, tree, "Member_all_1",
    ##                          sapply(2:32, function (j) { sprintf("Profile_2.%d", j) } ),c
    ##                          seq(300,0,-10),
    ##                          sapply(2:32, function (j) { sprintf("Profile_2.%d", j) } ),
    ##                          seq(300,0,-10),
    ##                          useUncorrectedCorrelation=FALSE )

    ## correlationBetweenRanges( traits, tree, "Member_Bacteria_2",
    ##                          sapply(seq(1,32,1), function (j) { sprintf("Profile_2.%d", j) } ),
    ##                          sapply(seq(1,32,1), function (j) { sprintf("Profile_2.%d", j) } ),
    ##                          useUncorrectedCorrelation=FALSE )

    ## dev.off()
    ## quit()


    ## for( gr in taxGroups )
    ## {
    ##     print(gr)
    ##     correlationBetweenRanges( traits, tree, gr,
    ##                              sapply(1:31, function (j) { sprintf("Profile_1.%d", j) } ),
    ##                              sapply(1:31, function (j) { sprintf("Profile_1.%d", j) } ),
    ##                              useUncorrectedCorrelation=TRUE )
    ## }
}

##quit()





#-----------------------------------------------------------------------------------------------------------------------------------
# Figure 3E - correlations between model regions
#-----------------------------------------------------------------------------------------------------------------------------------

correlationBetweenModelRegions <- function( traits, tree, filterTrait, ranges, useUncorrectedCorrelation=FALSE, halfOnly=0 )
{
    # 1=top 0=both -1=bottom
    #stopifnot( length(xRange)==length(xCoor) )
    #stopifnot( length(yRange)==length(yCoor) )
    stopifnot( (halfOnly <= 1) && (halfOnly >= -1) )
    corValues <- data.frame( Var1=character(), Var2=character(), Cor=double(), DirectionsIndicator=character(), Pval=double() )
    ##                          sapply(seq(32,2,-1), function (j) { sprintf("Profile_2.%d", j) } ),
##                          seq(300,0,-10),

    for( xi in 1:nrow(ranges) )
    {
        xFrom <- ranges[xi,"from"]
        xTo <- ranges[xi,"to"]
        xProfile <- ranges[xi,"profileId"]
        print(c(xFrom, xTo))
        print(xProfile)
        print(class(xProfile))

        xVars <- getProfileVariables( c(xFrom, xTo), profileId=xProfile )
        print(xVars)
        traits$RangeMean.x <- rowMeans( traits[xVars] )
        
        #xPos   <- xCoor[xi]
        
        for( yi in 1:nrow(ranges) )
        {

            skip <- FALSE
            if( xi>yi )
            {
                skip <- halfOnly > 0
            }
            else if( xi<yi )
            {
                skip <- halfOnly < 0
            }
            
            if( xi != yi && (!skip) )
            {

                yFrom <- ranges[yi,"from"]
                yTo <- ranges[yi,"to"]
                yProfile <- ranges[yi,"profileId"]
                print(c(yFrom, yTo))



                yVars <- getProfileVariables( c(yFrom, yTo), profileId=yProfile )
                print(yVars)
                traits$RangeMean.y <- rowMeans( traits[yVars] )

                print(sprintf("%d-%d   %d-%d", xFrom, xTo, yFrom, yTo))


                if( !useUncorrectedCorrelation )
                {
                    result <- performGLSregressionWithFilter( traits, tree, filterTrait, "RangeMean.x", "RangeMean.y", plotRegression=FALSE )
                    method <- "Method: Corrected (GLS)"
                }
                else
                {
                    result <- findUncorrectedVariableCorrelation( traits, filterTrait, "RangeMean.x", "RangeMean.y" )
                    result$R2 <- result$corr  # compatibility hack...
                    #result$pvalue = result$pvalue
                    method <- "Method: Uncorrected (Spearman)"
                }
            }
            else
            {
                result <- list(R2=1.0, directionsIndicator='+', pvalue=1.0)
            }

            if (xi <= yi )
            {
                corValues <- rbind( corValues, data.frame( Var1=xi, Var2=yi, Cor=result$R2, DirectionsIndicator=result$directionsIndicator, Pval=result$pvalue ) )
            }
            
        }
    }

    print(corValues)

    plotTitle <- sprintf("%s (N=%d)", filterTrait, result$N )


    corValues$Cor.fmt <- format(round(corValues$Cor,  digit=3))

    p <- ggplot( corValues, aes(x=Var1, y=Var2 ) ) + #, fill=Cor) ) +
        geom_tile( aes(fill=Cor), linetype=0, colour=NA ) +   # Use geom_tile (because geom_raster has poor pdf/svg output), with minimal borders
        geom_text( aes(label=Cor.fmt), color="white" ) +
        labs( title=plotTitle ) +
        scale_x_discrete( labels=c("1"="Start", "2"="Mid1", "3"="Mid2", "4"="End") ) +
        scale_y_discrete( labels=c("1"="Start", "2"="Mid1", "3"="Mid2", "4"="End") ) +
        scale_fill_gradientn(
            breaks=c(-1.0, -0.5, 0.0, 0.5, 1.0),
            limits=c(-1.0, 1.0),
            colors=        c("#75c0ff", "#2240cc", "#112060",   "#000000", "#602011", "#cc4022", "#ffc075"),
            values=rescale(c(     -1.0,      -0.5,     -0.15,         0.0,      0.15,       0.5,       1.0) ) ) +
        theme( plot.background = element_blank(),   # Hide unnecessary theme elements (background panels, etc.)
#              panel.grid.major = element_blank(),
#              panel.grid.minor = element_blank(),
              panel.background = element_blank()
              )

    dotsDF <- data.frame(melt( corValues, id.vars=c("Var1", "Var2"), measure.vars=c("Pval") ))
    dotsDF <- dotsDF[dotsDF$value<=0.01,]
    
    p <- p +
        geom_point( data=dotsDF, aes(x=Var1,y=Var2), position=position_nudge(y=0.4, x=-0.4), color="white", size=((fontScale/12)**2)*1.8 )  # Significance markers (white dots)
    
    print(p)

    corValues
}

figure_CorrelationBetweenModelRegions <- function()
{

    #-----------------------------------------------------------------------------------------------------------------------------------
    # Plot autocorrelation graph (by regions)
    #-----------------------------------------------------------------------------------------------------------------------------------
    range.id      <- 1:4
    range.from    <- c(       1,    11,       2,     31 )
    range.to      <- c(       2,    31,      22,     32 )
    range.profile <- c(       1,     1,       2,      2 )
    range.label   <- c( "Start", "Mid1", "Mid2",  "End" )

    ranges <- data.frame( id=range.id, from=range.from, to=range.to, profileId=range.profile, label=range.label )


    corValues <- correlationBetweenModelRegions( traits, tree, "Member_all_1", ranges, useUncorrectedCorrelation=TRUE, halfOnly=0 )

    # print p-values in scientific notation (to see if the unrounded estimated values are all really 0.00000000...)
    x <- corValues[corValues$Var1==1 & corValues$Var2==2, "Pval"]
    print(sprintf("%e", x))
    x <- corValues[corValues$Var1==1 & corValues$Var2==3, "Pval"]
    print(sprintf("%e", x))
    x <- corValues[corValues$Var1==1 & corValues$Var2==4, "Pval"]
    print(sprintf("%e", x))
    x <- corValues[corValues$Var1==2 & corValues$Var2==3, "Pval"]
    print(sprintf("%e", x))
    x <- corValues[corValues$Var1==2 & corValues$Var2==4, "Pval"]
    print(sprintf("%e", x))
    x <- corValues[corValues$Var1==3 & corValues$Var2==4, "Pval"]
    print(sprintf("%e", x))

}

#dev.off()
#quit()
#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------

## performGLSregression_profileRangeMean( traits, tree, "Profile_1.1",  as.integer(c(8,15)), profileId=1, plotRegression=TRUE, caption="Start-referenced, 0nt vs mean(70-140nt)" )

## performGLSregression_profileRangeMean( traits, tree, "Profile_2.31", as.integer(c(17,24)),  profileId=2, plotRegression=TRUE, caption="End-referenced, 0 nt vs mean(70-140nt)" )



#-----------------------------------------------------------------------------------------------------------------------------------
# Plot correlogram - correlation vs. phylogenetic distance for different traits (supp. figure)
#-----------------------------------------------------------------------------------------------------------------------------------
getCorrelogram <- function( tree, traits, traitName, n.points=10, ci.bs=99 )
{
    speciesWithMissingData <- row.names(traits[is.na(traits[traitName]),])
    tree <- drop.tip( tree, speciesWithMissingData )   # Filter species that will prevent analysis from being performed from the tree

    treeWithTraits <- phylo4d(tree, tip.data=traits, rownamesAsLabels=TRUE)

    phyloCorrelogram( treeWithTraits, trait=traitName, n.points=n.points, ci.bs=ci.bs )
}

makeCorrelogramComparisonFigure <- function(tree, traits, n.points=20, ci.bs=99 )
{
    profileVars <- getProfileVariables( c(1, 31), profileId=1 )
    vars <- c( profileVars[ seq(1, 31, 15) ], "GenomicGC", "LogGenomeSize", "Noise" )
    #vars <- c( profileVars[ seq(1, 31, 15) ] )
    
    df <- data.frame()
    
    for( var in vars )
    {
        correlogram <- getCorrelogram( tree, traits, var, n.points=n.points, ci.bs=ci.bs )

        newrow1 <- data.frame( ProfileStart=ordered(c(var)), Correlogram=t(correlogram$res[,4]), type=c("val"))
        colnames(newrow1) <- c("Trait", vapply(1:n.points, function(x) as.character(x), character(1)), "type")
        newrow2 <- data.frame( ProfileStart=ordered(c(var)), Correlogram=t(correlogram$res[,2]), type=c("min"))
        colnames(newrow2) <- c("Trait", vapply(1:n.points, function(x) as.character(x), character(1)), "type")
        newrow3 <- data.frame( ProfileStart=ordered(c(var)), Correlogram=t(correlogram$res[,3]), type=c("max"))
        colnames(newrow3) <- c("Trait", vapply(1:n.points, function(x) as.character(x), character(1)), "type")
        
        print(newrow1)
        print(newrow2)
        print(newrow3)
        df <- rbind( df, newrow1, newrow2, newrow3)
    }

    df$Distance <- df$variable

    return(df)
}

# Calculate the correlgoram data
correlogram.steps <- 20
correlogramData <- NA
## correlogramData <- makeCorrelogramComparisonFigure( tree, traits, n.points=correlogram.steps, ci.bs=500 )
## # Combine the confidence interval with the estimated values
## dfCorrelogram1 <- melt( correlogramData[correlogramData$type=="val",], id.vars=c("Trait", "type") )
## dfCorrelogram2 <- melt( correlogramData[correlogramData$type=="max",], id.vars=c("Trait", "type") )
## dfCorrelogram3 <- melt( correlogramData[correlogramData$type=="min",], id.vars=c("Trait", "type") )
## dfCorrelogram <- merge( merge( dfCorrelogram1, dfCorrelogram2, by=c("Trait", "variable") ),  dfCorrelogram3, by=c("Trait", "variable") )
## dfCorrelogram$pos <- as.numeric(dfCorrelogram$variable)/correlogram.steps*2.0 # convert correlogram "variable names" (indicating the different positions) to positions
## print(dfCorrelogram)

## if( isTRUE(is.na(correlogramData)) )
## {
    
##     p <- ggplot( dfCorrelogram, aes(x=Trait, y=variable, fill=value.x) ) +
##         geom_raster()
##     print(p)


##     #var2treedist_trans <- function () scales::trans_new( "var2treedist", function (x) as.numeric(x)*0.5, function (x) as.numeric(x)/0.5 )
##     #    coord_trans( x="var2treedist" ) +

##     cs1 <- c("Profile_1.1"="#ff4070", "Profile_1.16"="#c04020", "Profile_1.31"="#b06020", "GenomicGC"="#2040d0", "LogGenomeSize"="#20d060", "Noise"="#809080")
##     cl1 <- c( "DeltaLFE (CDS start+0nt)", "DeltaLFE (CDS start+150nt)", "DeltaLFE (CDS start+300nt)", "GenomicGC", "log(GenomeSize)", "Random")
##     p <- ggplot( dfCorrelogram, aes( x=pos ) ) +
##         geom_hline( yintercept=0, color="black" ) +
##         geom_ribbon( aes( ymin=value, ymax=value.y, fill=Trait, group=Trait), color=NA, alpha=0.3 ) +
##         geom_line( aes( y=value,   color=Trait, group=Trait), size=0.5, alpha=0.4 ) +
##         geom_line( aes( y=value.y, color=Trait, group=Trait), size=0.5, alpha=0.4 ) +
##         geom_line( aes( y=value.x, color=Trait, group=Trait,  size=Trait), alpha=0.9 ) +
##         scale_size_manual( values=c("Profile_1.1"=1.3, "Profile_1.6"=1.3, "Profile_1.11"=1.3, "Profile_1.16"=1.3, "Profile_1.21"=1.3, "Profile_1.26"=1.3, "Profile_1.31"=1.3, "GenomicGC"=1.0, "LogGenomeSize"=1.0, "Noise"=1.0 ) ) +
##         scale_colour_manual( values=cs1, labels=cl1 ) +
##         scale_fill_manual(   values=cs1, labels=cl1 ) +
##         guides( size="none" ) +
##         labs( title="Correlogram for traits (correlation vs. phylogenetic distance)",
##              x="Relative phylogenetic distance (arbitrary units)",
##              y="Correlation in trait values (between species)" ) +
##         theme( plot.background = element_blank(),   # Hide unnecessary theme elements (background panels, etc.)
##               panel.grid.major.y = element_blank(), # element_line(color="grey", size=0.50),
##               panel.grid.major.x = element_blank(), # element_line(color="grey", size=0.50),
##               panel.grid.minor = element_blank(),
##               panel.background = element_blank()
##               )
##     print(p)
## }

#dev.off()
#quit()

#-----------------------------------------------------------------------------------------------------------------------------------





performGLSregression_profileRangeMean_withFilter <- function( traits, tree, filterTrait, studyTrait,  profileRange, profileId=1, plotRegression=TRUE , caption="", plot.yrange=NA)
{
    # Filter tree by trait
    ret <- prepareFilteredTreeForGLS( traits, tree, filterTrait, studyTrait )
    if( isTRUE(is.na(ret[[1]])) ) { return( NA ) }
    traits.filtered <- ret[[1]]
    tree.filtered   <- ret[[2]]
    stopifnot( nTips(tree.filtered) == nrow(traits.filtered) ) # after filtering, the tree must still match the traits array
        
    if( is.null(tree.filtered) || nTips(tree.filtered) < minimalTaxonSize )
    {
        return( NA )
    }

    if( any(class(traits.filtered[,studyTrait])==c("factor")) && nlevels(as.factor(as.numeric(traits.filtered[,studyTrait]))) < 2 )
    {
        print("Only one level remaining...")
        return( NA )
    }
    
    ###return( glsRegressionRangeAnalysis( traits.filtered, tree.filtered, studyTrait, pyramidSpec, sprintf("%s\n(N=%d)", filterTrait, nrow(traits) ), profileId=profileId ) )
    return( performGLSregression_profileRangeMean( traits.filtered, tree.filtered, studyTrait,  profileRange, profileId=profileId, plotRegression=plotRegression, caption=sprintf("%s - %s; %s (N=%d)", filterTrait, caption, studyTrait, nrow(traits.filtered) ), traits.full=traits, plot.yrange=plot.yrange ) )
}


performOLSregression_profileRangeMean_withFilter <- function( traits, filterTrait, studyTrait,  profileRange, profileId=1, plotRegression=TRUE , caption="", colorTrait="Member_all_1" )
{
    # Filter species by trait
    traits.filtered <- traits[as.logical(traits[,filterTrait]),]
    
    return( performOLSregression_profileRangeMean( traits.filtered, studyTrait,  profileRange, profileId=profileId, plotRegression=plotRegression, caption=sprintf("%s - %s; %s (N=%d)", filterTrait, caption, studyTrait, nrow(traits.filtered) ), colorTrait=colorTrait ) )
}

regressionResultsByTaxGroup <- data.frame( ExplanatoryVar=character(), Range=integer(), MaxRangeStart=integer(), MaxRangeEnd=integer(), TaxGroup=integer(), TaxGroupName=character(), EffectSize=double(), Pvalue=double(), NumSpecies=integer() )

## for( gr in taxGroups )
## {
##     print(gr)
    

##     performGLSregression_profileRangeMean_withFilter( traits, tree, gr, "Profile_1.1",  as.integer(c(12,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced" )

##     performGLSregression_profileRangeMean_withFilter( traits, tree, gr, "Profile_2.32", as.integer(c(1,15)),  profileId=2, plotRegression=TRUE, caption="End-referenced" )
        
##     #performGLSregressionWithFilter( traits, tree, gr,                  "Profile_1.1",   "Profile_1.15" )

##     ## # Store peaks (for each range) for this group
##     ## if( any( class(results) == "data.frame" ) && nrow(results) )
##     ## {
##     ##     results <- results[results$Var1==results$Var2,]    # Only include single-window results

##     ##     # Iterate over each range to find the relevant peak
##     ##     for( i in 1:length(groupsTableOutputFile.limitRangeFromNt) )
##     ##     {
##     ##         matching <- (results$Var1 >= groupsTableOutputFile.limitRangeFromNt[i]) &
##     ##                     (results$Var2 <= groupsTableOutputFile.limitRangeToNt[i]  )  # Only include ranges within the configured limits

##     ##         maxResult <- results[matching,][which.max( abs(results[matching, "Buse.R2"]) ),]  # Choose the result with the highest R^2 

##     ##         regressionResultsByTaxGroup <- rbind( regressionResultsByTaxGroup, data.frame( ExplanatoryVar=c("GenomicGC"),  Range=c(i-1), MaxRangeStart=c(maxResult$Var1), MaxRangeEnd=c(maxResult$Var2), TaxGroup=c(taxGroupToTaxId[gr]), TaxGroupName=c(gr), EffectSize=c(maxResult$Buse.R2), Pvalue=c(maxResult$Pvalue), NumSpecies=c(maxResult$NumSpecies) ) )
##     ##     }
##     ## }
## }


# Side-by-side comparison of begin- and end-referenced profiles

## for( gr in taxGroups )
## {
##     print(gr)
##     results <- glsRangeAnalysisWithFilter( traits, tree, gr, "GenomicGC", c(1,31),   profileId=1)

##     results <- glsRangeAnalysisWithFilter( traits, tree, gr, "GenomicGC", c(2,32),   profileId=2)
    
##     #performGLSregressionWithFilter( traits, tree, gr,                  "Profile_1.1",   "Profile_1.15" )

##     # Store peaks (for each range) for this group
##     if( any( class(results) == "data.frame" ) && nrow(results) )
##     {
##         results <- results[results$Var1==results$Var2,]    # Only include single-window results

##         # Iterate over each range to find the relevant peak
##         for( i in 1:length(groupsTableOutputFile.limitRangeFromNt) )
##         {
##             matching <- (results$Var1 >= groupsTableOutputFile.limitRangeFromNt[i]) &
##                         (results$Var2 <= groupsTableOutputFile.limitRangeToNt[i]  )  # Only include ranges within the configured limits

##             maxResult <- results[matching,][which.max( abs(results[matching, "Buse.R2"]) ),]  # Choose the result with the highest R^2 

##             regressionResultsByTaxGroup <- rbind( regressionResultsByTaxGroup, data.frame( ExplanatoryVar=c("GenomicGC"),  Range=c(i-1), MaxRangeStart=c(maxResult$Var1), MaxRangeEnd=c(maxResult$Var2), TaxGroup=c(taxGroupToTaxId[gr]), TaxGroupName=c(gr), EffectSize=c(maxResult$Buse.R2), Pvalue=c(maxResult$Pvalue), NumSpecies=c(maxResult$NumSpecies) ) )
##         }
##     },
## }


## dev.off()
## quit()

plotContrastingKDEWithFilters <- function( traits, filterTraits, negation, profileRange, profileId, yrange=NA )
{
    #stopifnot( any( filterTrait==colnames(traits) ) )  # filterTrait not found
                                        #stopifnot( any( studyTrait==colnames(traits)  ) )  # studyTrait not found
    print(filterTraits)
    print(negation)
    stopifnot(length(negation) > 0)
    stopifnot(length(filterTraits) == length(negation) )

    # Filter species by trait
    negationValues <- matrix(rep( negation, each=nrow(traits) ), ncol=length(filterTraits) )  # values for negation calculation
    #print(apply(negationValues, MARGIN=2, FUN=sum))
    #print(dim(negationValues))
    #print(dim(traits))
    #print(length(negation))
    #print(traits[ filterTraits ])  # not specifying rows (',') ensures we get a matrix back regardless of the number of columns specified
    #print(dim(traits[ filterTraits ]))
    #print(dim(negationValues))
    booleanVal <- xor(traits[ filterTraits ], negationValues)
    #print(booleanVal)
    #print(class(booleanVal))
    print(apply( booleanVal, MARGIN=2, FUN=sum))
    traits.filtered <- traits[apply( booleanVal, MARGIN=1, FUN=all), ]
        
    if( nrow(traits.filtered) < minimalTaxonSize ) { return( NA ) }

    #print(traits.filtered)
    
    kdePlot <- kdePlot( traits.filtered[, getProfileVariables( profileRange, profileId ) ], getProfilePositions( profileId ), profileId=profileId, yrange=yrange  )

    return( kdePlot )

}

plotProfileBoxplotWithFilters <- function( traits, filterTraits, negation, profileRange, profileId, yrange=NA )
{
    #stopifnot( any( filterTrait==colnames(traits) ) )  # filterTrait not found
                                        #stopifnot( any( studyTrait==colnames(traits)  ) )  # studyTrait not found
    print(filterTraits)
    print(negation)
    stopifnot(length(negation) > 0)
    stopifnot(length(filterTraits) == length(negation) )

    # Filter species by trait
    negationValues <- matrix(rep( negation, each=nrow(traits) ), ncol=length(filterTraits) )  # values for negation calculation
    #print(apply(negationValues, MARGIN=2, FUN=sum))
    #print(dim(negationValues))
    #print(dim(traits))
    #print(length(negation))
    #print(traits[ filterTraits ])  # not specifying rows (',') ensures we get a matrix back regardless of the number of columns specified
    #print(dim(traits[ filterTraits ]))
    #print(dim(negationValues))
    booleanVal <- xor(traits[ filterTraits ], negationValues)
    #print(booleanVal)
    #print(class(booleanVal))
    print(apply( booleanVal, MARGIN=2, FUN=sum))
    traits.filtered <- traits[apply( booleanVal, MARGIN=1, FUN=all), ]
        
    if( nrow(traits.filtered) < minimalTaxonSize ) { return( NA ) }

    print(traits.filtered$Profile_1.1)
    print(nrow(traits.filtered))
          

    dd <- melt( traits.filtered[, getProfileVariables( profileRange, profileId ) ] )
    #print(dd)
    p <- ggplot( dd, aes(x=variable, y=value) )+
        geom_boxplot()
    print(p)
}


getFilteredProfileValues <- function( traits, filterTraits, negation, profileRange, profileId, yrange=NA )
{
    stopifnot(length(negation) > 0)
    stopifnot(length(filterTraits) == length(negation) )
    #print("--------------------------------------[0]")

    # Filter species by trait
    negationValues <- matrix(rep( negation, each=nrow(traits) ), ncol=length(filterTraits) )  # values for negation calculation
    rawVals        <- traits[ filterTraits ]   # Take all values for filter variables
    rawValsBoolean <- (rawVals != 0) & (rawVals != FALSE)  # convert to logicals (rawVals can contain logicals or factors)
    #print(rawValsBoolean)
    #print(rawValsBoolean)
    #print("~~~~~~~~~~~~~~~~~~~~~")
    booleanResults <- xor(rawValsBoolean, negationValues)
    #print(booleanResults)
    traits.filtered <- traits[apply( booleanResults, MARGIN=1, FUN=all), ]
        
    #if( nrow(traits.filtered) < minimalTaxonSize ) { return( NA ) }

    ret <-traits.filtered[, getProfileVariables( profileRange, profileId ) ]
    #print(getProfileVariables( profileRange, profileId ))
    #print(getProfilePositions( profileId=profileId ))
    #print(ncol(ret))
    #print("--------------------------------------[1]")
    colnames(ret) <- getProfilePositions( profileId=profileId ) # rename columns to numeric positions
    melt( ret )
}

figure_contrastingKDEsForHighLowGC <- function()
{
    #-----------------------------------------------------------------------------------------------------------------------------------
    # Plot KDEs for subgroups
    #-----------------------------------------------------------------------------------------------------------------------------------

    yrangeForKDEs = c(-1.5, 2.6)

    kdePlot1 <- plotContrastingKDEWithFilters( traits, c( "Member_Bacteria_2", "GC.45" ), c(FALSE, FALSE), as.integer(c(1,31)), profileId=1, yrange=yrangeForKDEs )

    
## kdePlot.all.b <- plotContrastingKDEWithFilters( traits, c( "Member_all_1" ), c(FALSE), as.integer(c(1,31)), profileId=1, yrange=yrangeForKDEs )
## kdePlot.all.e <- plotContrastingKDEWithFilters( traits, c( "Member_all_1" ), c(FALSE), as.integer(c(2,32)), profileId=2, yrange=yrangeForKDEs )
##                                         #print(kdePlot.all)

## grid.arrange( ggplotGrob( kdePlot.all.b ), ggplotGrob( kdePlot.all.e ),
##               ggplotGrob( kdePlot.all.b ), ggplotGrob( kdePlot.all.e ),
##               ggplotGrob( kdePlot.all.b ), ggplotGrob( kdePlot.all.e ),
##              ncol=2, widths=c(unit(0.5, "npc"), unit(0.5, "npc")), heights=c(unit(0.333, "npc"), unit(0.334, "npc"), unit(0.333, "npc")) )

    ## kdePlot2 <- plotContrastingKDEWithFilters( traits, c( "Member_Bacteria_2", "GC.45" ), c(FALSE, TRUE),  as.integer(c(1,31)), profileId=1, yrange=yrangeForKDEs )

    ## kdePlot3 <- plotContrastingKDEWithFilters( traits, c( "Member_Archaea_2157", "GC.45" ), c(FALSE, FALSE), as.integer(c(1,31)), profileId=1, yrange=yrangeForKDEs )

    ## kdePlot4 <- plotContrastingKDEWithFilters( traits, c( "Member_Archaea_2157", "GC.45" ), c(FALSE, TRUE),  as.integer(c(1,31)), profileId=1, yrange=yrangeForKDEs )


    ## kdePlot5 <- plotContrastingKDEWithFilters( traits, c( "Member_Eukaryota_2759", "GC.45" ), c(FALSE, FALSE), as.integer(c(1,31)), profileId=1, yrange=yrangeForKDEs )

    ## kdePlot6 <- plotContrastingKDEWithFilters( traits, c( "Member_Eukaryota_2759", "GC.45" ), c(FALSE, TRUE), as.integer(c(1,31)), profileId=1, yrange=yrangeForKDEs )

    ## print(kdePlot1)

    ## ## #kdePlot4 <- plotContrastingKDEWithFilters( traits, c( "Member_Bacteria_2" ), c(FALSE), as.integer(c(1,31)), profileId=1 )
    ## ## #print(kdePlot1)

    ## kdePlot1e <- plotContrastingKDEWithFilters( traits, c( "Member_Bacteria_2", "GC.45" ), c(FALSE, FALSE), as.integer(c(2,32)), profileId=2, yrange=yrangeForKDEs )

    ## kdePlot2e <- plotContrastingKDEWithFilters( traits, c( "Member_Bacteria_2", "GC.45" ), c(FALSE, TRUE),  as.integer(c(2,32)), profileId=2, yrange=yrangeForKDEs )

    ## kdePlot3e <- plotContrastingKDEWithFilters( traits, c( "Member_Archaea_2157", "GC.45" ), c(FALSE, FALSE), as.integer(c(2,32)), profileId=2, yrange=yrangeForKDEs )

    ## kdePlot4e <- plotContrastingKDEWithFilters( traits, c( "Member_Archaea_2157", "GC.45" ), c(FALSE, TRUE),  as.integer(c(2,32)), profileId=2, yrange=yrangeForKDEs )


    ## kdePlot5e <- plotContrastingKDEWithFilters( traits, c( "Member_Eukaryota_2759", "GC.45" ), c(FALSE, FALSE), as.integer(c(2,32)), profileId=2, yrange=yrangeForKDEs )

    ## kdePlot6e <- plotContrastingKDEWithFilters( traits, c( "Member_Eukaryota_2759", "GC.45" ), c(FALSE, TRUE), as.integer(c(2,32)), profileId=2, yrange=yrangeForKDEs )


    ## grid.arrange( ggplotGrob( kdePlot1 ), ggplotGrob( kdePlot1e ),
    ##              ggplotGrob( kdePlot2 ), ggplotGrob( kdePlot2e ),
    ##              ggplotGrob( kdePlot3 ), ggplotGrob( kdePlot3e ),
    ##              ggplotGrob( kdePlot4 ), ggplotGrob( kdePlot4e ),
    ##              ggplotGrob( kdePlot5 ), ggplotGrob( kdePlot5e ),
    ##              ggplotGrob( kdePlot6 ), ggplotGrob( kdePlot6e ),
    ##              ncol=2, widths=c(unit(0.5, "npc"), unit(0.5, "npc")), heights=c(unit(0.166, "npc"), unit(0.166, "npc"), unit(0.166, "npc"), unit(0.166, "npc"), unit(0.166, "npc"), unit(0.166, "npc") ) )

    
}


plotConstrastingProfilesAsBoxplot <- function( d1, d2, labels, ylimits=c(-1,2), profileId=1 )
{
    print(labels)
    print(profileId)
    stopifnot(!isTRUE(is.na(d1)))
    stopifnot(!isTRUE(is.na(d2)))
    # Add grouping variable (will be used to separate the two datasets and plot by color)
    d1$group <- as.factor(1)
    d2$group <- as.factor(2)
    d.combined <- rbind(d1, d2)  # combine
    d.combined <- d.combined[!is.na(d.combined$value),]  # remove NAs (i.e., species without temperature data)

    # Use histogram to find out how many species exists in each group
    d.combined$pos <- as.numeric(d.combined$variable) # add position as numeric (instead of factor) - needed by hist()
    allPositions <- unique( d.combined$pos )
    countsByPosition.gr1 <- hist( d.combined[d.combined$group==1,]$pos, plot=FALSE, breaks=c(-1000,allPositions) )$counts
    countsByPosition.gr2 <- hist( d.combined[d.combined$group==2,]$pos, plot=FALSE, breaks=c(-1000,allPositions) )$counts
    count.gr1 <- min(countsByPosition.gr1)
    count.gr2 <- min(countsByPosition.gr2)

    p <- ggplot( d.combined, aes(x=variable, y=value, color=group, fill=group) )+
        geom_hline( yintercept=0, color="black") +
        geom_boxplot() +
        scale_colour_manual( values=c("1"="#334e73", "2"="#9c180a") ) +
        scale_fill_manual(  values=c("1"="#6f90bb", "2"="#e2520e"), labels=c("1"=sprintf("%s (N=%d)", labels[1], count.gr1), "2"=sprintf("%s (N=%d)", labels[2], count.gr2) ) ) +
        labs( x="CDS position (nt)", y="Delta-LFE (kcal/mol/window)" ) +
        guides( colour=FALSE ) +
        scale_y_continuous( limits=ylimits, breaks=seq(-1.0,2.0,0.5) ) +
        theme( plot.background = element_blank(),   # Hide unnecessary theme elements (background panels, etc.)
              panel.grid.major.y = element_line(color="grey", size=0.50),
              panel.grid.major.x = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              aspect.ratio=0.75
              ) # +

    if( profileId==1 )
    {
        p <- p +
            scale_x_discrete( breaks=c("0","100","200","300") )
    }
    else
    {
        p <- p +
            scale_x_discrete( breaks=c("-300","-200","-100","0") )
    }
    
    print(p)
    p
}


figure_contrastingBoxplotsHighLowTemp <- function()
{
    
    yrangeForKDEs = c(-1.5, 2.6)

    #print(traits$Profile_1.5)
    #print(traits$GC.45)
    #print(class(traits$GC.45))
    #print(traits$GC.45 != 0)

    p1 <- plotConstrastingProfilesAsBoxplot(
        getFilteredProfileValues( traits, c( "Member_Bacteria_2", "GC.45" ), c(FALSE, TRUE ), as.integer(c(1,31)), profileId=1, yrange=yrangeForKDEs ),
        getFilteredProfileValues( traits, c( "Member_Bacteria_2", "GC.45" ), c(FALSE, FALSE ), as.integer(c(1,31)), profileId=1, yrange=yrangeForKDEs ), c("GC>=45", "GC<45"), ylimits=c(-1,2), profileId=1 )

    performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_all_1", "TempHighLow75",  as.integer(c(11,31)), profileId=1, plotRegression=TRUE )
    #performOLSregression_profileRangeMean_withFilter( traits, tree, "Member_all_1", "TempHighLow75",  as.integer(c(11,31)), profileId=1, plotRegression=TRUE )
    performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_all_1", "TempHighLow75",  as.integer(c( 2,22)), profileId=2, plotRegression=TRUE )
    #performOLSregression_profileRangeMean_withFilter( traits, tree, "Member_all_1", "TempHighLow75",  as.integer(c( 2,22)), profileId=2, plotRegression=TRUE )

    performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_all_1", "TempHighLow75",  as.integer(c( 1, 2)), profileId=1, plotRegression=TRUE )
    performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_all_1", "TempHighLow75",  as.integer(c(31,32)), profileId=2, plotRegression=TRUE )


    performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_all_1", "Is_endosymbiont",  as.integer(c( 1, 2)), profileId=1, plotRegression=TRUE, caption="Begin+0-10nt" )
    performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_all_1", "Is_endosymbiont",  as.integer(c(11,31)), profileId=1, plotRegression=TRUE, caption="Begin+100-300nt" )
    performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_all_1", "Is_endosymbiont",  as.integer(c(31,32)), profileId=2, plotRegression=TRUE, caption="End-0-10nt" )
    performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_all_1", "Is_endosymbiont",  as.integer(c( 2,22)), profileId=2, plotRegression=TRUE, caption="End-100-300nt" )

    performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_all_1", "OptimumTemp",  as.integer(c( 11, 31)), profileId=1, plotRegression=TRUE, caption="Begin+100-300nt", plot.yrange=c(-0.75,0.07) )
    performOLSregression_profileRangeMean_withFilter( traits,       "Member_all_1", "OptimumTemp",  as.integer(c( 11, 31)), profileId=1, plotRegression=TRUE, caption="Begin+100-300nt" )
    performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_all_1", "OptimumTemp",  as.integer(c( 2, 22)), profileId=2, plotRegression=TRUE, caption="End-100-300nt", plot.yrange=c(-0.75,0.07) )
    performOLSregression_profileRangeMean_withFilter( traits,       "Member_all_1", "OptimumTemp",  as.integer(c( 2, 22)), profileId=2, plotRegression=TRUE, caption="End-100-300nt" )
    
    p1 <- plotConstrastingProfilesAsBoxplot(
        getFilteredProfileValues( traits, c( "Member_all_1", "TempHighLow75" ), c(FALSE, TRUE ),  as.integer(c(1,31)), profileId=1, yrange=yrangeForKDEs ),
        getFilteredProfileValues( traits, c( "Member_all_1", "TempHighLow75" ), c(FALSE, FALSE ), as.integer(c(1,31)), profileId=1, yrange=yrangeForKDEs ), c("Optimum-temp <= 75", "Optimum-temp  > 75"), ylimits=c(-1,2), profileId=1 )

    p1 <- plotConstrastingProfilesAsBoxplot(
        getFilteredProfileValues( traits, c( "Member_all_1", "TempHighLow75" ), c(FALSE, TRUE ),  as.integer(c(1,32)), profileId=2, yrange=yrangeForKDEs ),
        getFilteredProfileValues( traits, c( "Member_all_1", "TempHighLow75" ), c(FALSE, FALSE ), as.integer(c(1,32)), profileId=2, yrange=yrangeForKDEs ), c("Optimum-temp <= 75", "Optimum-temp  > 75"), ylimits=c(-1,2), profileId=2 )
    
    p1 <- plotConstrastingProfilesAsBoxplot(
        getFilteredProfileValues( traits, c( "Member_all_1", "Is_high_temp" ), c(FALSE, TRUE ),  as.integer(c(1,31)), profileId=1, yrange=yrangeForKDEs ),
        getFilteredProfileValues( traits, c( "Member_all_1", "Is_high_temp" ), c(FALSE, FALSE ), as.integer(c(1,31)), profileId=1, yrange=yrangeForKDEs ), c("Temp<58", "Temp>=58" ), ylimits=c(-1,2), profileId=1 )
    

    p1 <- plotConstrastingProfilesAsBoxplot(
        getFilteredProfileValues( traits, c( "Member_all_1", "Is_endosymbiont" ), c(FALSE, TRUE ),  as.integer(c(1,31)), profileId=1, yrange=yrangeForKDEs ),
        getFilteredProfileValues( traits, c( "Member_all_1", "Is_endosymbiont" ), c(FALSE, FALSE ), as.integer(c(1,31)), profileId=1, yrange=yrangeForKDEs ), c("Rest", "Endosymbiont"), ylimits=c(-1,2), profileId=1 )

    p1 <- plotConstrastingProfilesAsBoxplot(
        getFilteredProfileValues( traits, c( "Member_all_1", "Is_endosymbiont" ), c(FALSE, TRUE ),  as.integer(c(1,32)), profileId=2, yrange=yrangeForKDEs ),
        getFilteredProfileValues( traits, c( "Member_all_1", "Is_endosymbiont" ), c(FALSE, FALSE ), as.integer(c(1,32)), profileId=2, yrange=yrangeForKDEs ), c("Rest", "Endosymbiont"), ylimits=c(-1,2), profileId=2 )

    #print( getFilteredProfileValues( traits, c( "Member_Gammaproteobacteria_1236", "Is_endosymbiont" ), c(FALSE, FALSE ), as.integer(c(1,31)), profileId=1, yrange=yrangeForKDEs ) )
    
    p1 <- plotConstrastingProfilesAsBoxplot(
        getFilteredProfileValues( traits, c( "Member_Gammaproteobacteria_1236", "Is_endosymbiont" ), c(FALSE, TRUE ),  as.integer(c(1,31)), profileId=1, yrange=yrangeForKDEs ),
        getFilteredProfileValues( traits, c( "Member_Gammaproteobacteria_1236", "Is_endosymbiont" ), c(FALSE, FALSE ), as.integer(c(1,31)), profileId=1, yrange=yrangeForKDEs ), c("Rest", "Endosymbiont"), ylimits=c(-1,2), profileId=1 )

    p1 <- plotConstrastingProfilesAsBoxplot(
        getFilteredProfileValues( traits, c( "Member_Gammaproteobacteria_1236", "Is_endosymbiont" ), c(FALSE, TRUE ),  as.integer(c(1,32)), profileId=2, yrange=yrangeForKDEs ),
        getFilteredProfileValues( traits, c( "Member_Gammaproteobacteria_1236", "Is_endosymbiont" ), c(FALSE, FALSE ), as.integer(c(1,32)), profileId=2, yrange=yrangeForKDEs ), c("Rest", "Endosymbiont"), ylimits=c(-1,2), profileId=2 )

    p1 <- plotConstrastingProfilesAsBoxplot(
        getFilteredProfileValues( traits, c( "Member_Proteobacteria_1224", "Is_endosymbiont" ), c(FALSE, TRUE ),  as.integer(c(1,31)), profileId=1, yrange=yrangeForKDEs ),
        getFilteredProfileValues( traits, c( "Member_Proteobacteria_1224", "Is_endosymbiont" ), c(FALSE, FALSE ), as.integer(c(1,31)), profileId=1, yrange=yrangeForKDEs ), c("Rest", "Endosymbiont"), ylimits=c(-1,2), profileId=1 )

    p1 <- plotConstrastingProfilesAsBoxplot(
        getFilteredProfileValues( traits, c( "Member_Proteobacteria_1224", "Is_endosymbiont" ), c(FALSE, TRUE ),  as.integer(c(1,32)), profileId=2, yrange=yrangeForKDEs ),
        getFilteredProfileValues( traits, c( "Member_Proteobacteria_1224", "Is_endosymbiont" ), c(FALSE, FALSE ), as.integer(c(1,32)), profileId=2, yrange=yrangeForKDEs ), c("Rest", "Endosymbiont"), ylimits=c(-1,2), profileId=2 )


}

plotContrastingProfilesAsSideBySideBoxplots <- function( f1, f2, labels, ylimits, show.guides=FALSE )
{
    print("                           (1)")
    f1$Group <- factor(labels[1])
    print("                           (1.0)")
    print(f1)
    if( length(labels) > 1 )
    {
        print("                           (1a)")
        f2$Group <- factor(labels[2])
        ff <- rbind( f1, f2 )
    }
    else
    {
        print("                           (1b)")
        ff <- f1
    }
    print("                           (1c)")
    print(ff)
    print("                           (1d)")
    print(names(ff))
    print("                           (2)")

    if( length(labels) > 1 )
    {
        print("                           (3a)")
        # Create contrasting color scheme
        cc.line <- c("#bca200", "#00139e") 
        names(cc.line) <- labels
        cc.fill <- c("#e2c200", "#0010ea")
        names(cc.fill) <- labels
    }
    else
    {
        print("                           (3b)")
        # Create color scheme for single population
        cc.line <- c("#1e2d46") 
        names(cc.line) <- labels
        cc.fill <- c("#4c72b0")  # 76 114 176
        names(cc.fill) <- labels
    }
    
    
    print("                           (4)")

    p <- ggplot( ff, aes(x=variable, y=value) )+
        geom_hline( yintercept=0, color="black") +
        geom_boxplot( aes(fill=Group, color=Group), outlier.size=0.4 ) +
                                        #       geom_jitter( aes(color=C), alpha=0.3, size=0.4 ) +
        scale_x_discrete( breaks=c("0","100","200","300") ) +
        scale_y_continuous( limits=ylimits ) +
        scale_color_manual( values=cc.line ) +
        scale_fill_manual(  values=cc.fill ) +
        theme( plot.background = element_blank(),   # Hide unnecessary theme elements (background panels, etc.)
              panel.grid.major.y = element_line(color="grey", size=0.50),
              panel.grid.major.x = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank()
              ) # +


    print("                           (5)")
    if( show.guides )
    {
        p <- p + labs( x="CDS position (nt)", y="Delta-LFE (kcal/mol/window)" )
    }
    else
    {
        p <- p + guides( fill=FALSE, color=FALSE ) +     # Hide the legend (will be plotted separately)
            labs( x="", y="" )
    }
    
    print(p)
    return(p)
}

figure_ContrasingProfileBoxplotsForHighLowGC <- function()
{
    
    dd1a <- plotContrastingProfilesAsSideBySideBoxplots(
        getFilteredProfileValues( traits, c( "Member_Bacteria_2", "GC.45" ), c(FALSE, TRUE ), as.integer(c(1,31)), profileId=1, yrange=yrangeForKDEs ),
        getFilteredProfileValues( traits, c( "Member_Bacteria_2", "GC.45" ), c(FALSE, FALSE), as.integer(c(1,31)), profileId=1, yrange=yrangeForKDEs ),
        c("GC<=45", "GC>45"), ylimits=c(-1.3,3) )

    dd2a <- plotContrastingProfilesAsSideBySideBoxplots(
        getFilteredProfileValues( traits, c( "Member_Archaea_2157", "GC.45" ), c(FALSE, TRUE ), as.integer(c(1,31)), profileId=1, yrange=yrangeForKDEs ),
        getFilteredProfileValues( traits, c( "Member_Archaea_2157", "GC.45" ), c(FALSE, FALSE), as.integer(c(1,31)), profileId=1, yrange=yrangeForKDEs ),   c("GC<=45", "GC>45"), ylimits=c(-1.3,2) )

    dd3a <- plotContrastingProfilesAsSideBySideBoxplots(
        getFilteredProfileValues( traits, c( "Member_Eukaryota_2759", "GC.45" ), c(FALSE, TRUE ), as.integer(c(1,31)), profileId=1, yrange=yrangeForKDEs ),
        getFilteredProfileValues( traits, c( "Member_Eukaryota_2759", "GC.45" ), c(FALSE, FALSE), as.integer(c(1,31)), profileId=1, yrange=yrangeForKDEs ),    c("GC<=45", "GC>45"), ylimits=c(-1.3,2) )


    dd1e <- plotContrastingProfilesAsSideBySideBoxplots(
        getFilteredProfileValues( traits, c( "Member_Bacteria_2", "GC.45" ), c(FALSE, TRUE ), as.integer(c(1,32)), profileId=2, yrange=yrangeForKDEs ),
        getFilteredProfileValues( traits, c( "Member_Bacteria_2", "GC.45" ), c(FALSE, FALSE), as.integer(c(1,32)), profileId=2, yrange=yrangeForKDEs ),
       c("GC<=45", "GC>45"), ylimits=c(-1.3,3) )

    dd2e <- plotContrastingProfilesAsSideBySideBoxplots(
        getFilteredProfileValues( traits, c( "Member_Archaea_2157", "GC.45" ), c(FALSE, TRUE ), as.integer(c(1,32)), profileId=2, yrange=yrangeForKDEs ),
        getFilteredProfileValues( traits, c( "Member_Archaea_2157", "GC.45" ), c(FALSE, FALSE), as.integer(c(1,32)), profileId=2, yrange=yrangeForKDEs ),
       c("GC<=45", "GC>45"), ylimits=c(-1.3,2)  )

    dd3e <- plotContrastingProfilesAsSideBySideBoxplots(
        getFilteredProfileValues( traits, c( "Member_Eukaryota_2759", "GC.45" ), c(FALSE, TRUE ), as.integer(c(1,32)), profileId=2, yrange=yrangeForKDEs ),
        getFilteredProfileValues( traits, c( "Member_Eukaryota_2759", "GC.45" ), c(FALSE, FALSE), as.integer(c(1,32)), profileId=2, yrange=yrangeForKDEs ),    c("GC<=45", "GC>45"), ylimits=c(-1.3,2) )


    plotContrastingProfilesAsSideBySideBoxplots(
        getFilteredProfileValues( traits, c( "Member_Eukaryota_2759", "GC.45" ), c(FALSE, TRUE ), as.integer(c(1,32)), profileId=2, yrange=yrangeForKDEs ),
        getFilteredProfileValues( traits, c( "Member_Eukaryota_2759", "GC.45" ), c(FALSE, FALSE), as.integer(c(1,32)), profileId=2, yrange=yrangeForKDEs ),    c("GC<=45", "GC>45"), ylimits=c(-1.3,2), show.guides=TRUE )



    grid.arrange( ggplotGrob( dd1a ), ggplotGrob( dd1e ),
                  ggplotGrob( dd2a ), ggplotGrob( dd2e ),
                  ggplotGrob( dd3a ), ggplotGrob( dd2e ),
                  ggplotGrob( dd3a ), ggplotGrob( dd2e ),
                  ncol=2, widths=c(unit(0.5, "npc"), unit(0.5, "npc")), heights=c(unit(0.2625, "npc"), unit(0.21875, "npc"), unit(0.21875, "npc"), unit(0.3, "npc") ) )


    dd0a <- plotContrastingProfilesAsSideBySideBoxplots(
        getFilteredProfileValues( traits, c( "Member_all_1" ), c(FALSE ), as.integer(c(1,31)), profileId=1, yrange=yrangeForKDEs ),
        data.frame(),
        c("All"), ylimits=c(-1.3,3) )

    dd0e <- plotContrastingProfilesAsSideBySideBoxplots(
        getFilteredProfileValues( traits, c( "Member_all_1" ), c(FALSE ), as.integer(c(1,32)), profileId=2, yrange=yrangeForKDEs ),
        data.frame(),
        c("All"), ylimits=c(-1.3,3) )

    grid.arrange( ggplotGrob( dd0a ), ggplotGrob( dd0e ),
                  ggplotGrob( dd0a ), ggplotGrob( dd0e ),
                  ggplotGrob( dd0a ), ggplotGrob( dd0e ),
                  ggplotGrob( dd0a ), ggplotGrob( dd0e ),
                  ncol=2, widths=c(unit(0.5, "npc"), unit(0.5, "npc")), heights=c(unit(0.32, "npc"), unit(0.227, "npc"), unit(0.226, "npc"), unit(0.227, "npc") ) )
}


##dev.off()
##quit()
 
#-----------------------------------------------------------------------------------------------------------------------------------


##performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_all_1", "GenomicGC",  as.integer(c(12,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced" )
##performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_Eukaryota_2759", "GenomicGC",  as.integer(c(12,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced" )
##performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_Fungi_4751", "GenomicGC",  as.integer(c(12,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced" )
##performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_Bacteria_2", "GenomicGC",  as.integer(c(12,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced" )


## performOLSregression( traits, "GenomicGC", "Profile_1.16" )

## performOLSregression_profileRangeMean( traits, "GenomicGC", as.integer(c(12,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced" )

## performOLSregression_profileRangeMean_withFilter( traits, "Member_all_1", "GenomicGC",  as.integer(c(12,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced", colorTrait="Member_Eukaryota_2759" )
## performOLSregression_profileRangeMean_withFilter( traits, "Member_Bacteria_2", "GenomicGC",  as.integer(c(12,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced" )
##performOLSregression_profileRangeMean_withFilter( traits, "Member_Eukaryota_2759", "GenomicGC",  as.integer(c(12,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced", colorTrait="Member_Fungi_4751" )
## performOLSregression_profileRangeMean_withFilter( traits, "Member_Opisthokonta_33154", "GenomicGC",  as.integer(c(12,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced", colorTrait="Member_Fungi_4751" )
##performOLSregression_profileRangeMean_withFilter( traits, "Member_Fungi_4751", "GenomicGC",  as.integer(c(12,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced", colorTrait="Member_Ascomycota_4890" )


##performOLSregression_profileRangeMean_withFilter( traits, "Member_Proteobacteria_1224", "GenomicGC",  as.integer(c(12,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced" )

## #olsRegressionRangeAnalysis( traits, "GenomicGC", c(1,pyramidLength), "GenomicGC")

## #olsRangeAnalysisWithFilter ( traits, "Member_Fungi_4751", "GenomicGC", c(1,31),   profileId=1)
## olsRangeAnalysisWithFilter ( traits, "Member_Fungi_4751", "GenomicGC", c(2,32),   profileId=2)

## dev.off()
## quit()

## performOLSregression_profileRangeMean_withFilter( traits, "Member_all_1", "GenomicGC",  as.integer(c(1,20)), profileId=2, plotRegression=TRUE, caption="Start-referenced", colorTrait="Member_Eukaryota_2759" )
## performOLSregression_profileRangeMean_withFilter( traits, "Member_Bacteria_2", "GenomicGC",  as.integer(c(1,20)), profileId=2, plotRegression=TRUE, caption="Start-referenced" )
## performOLSregression_profileRangeMean_withFilter( traits, "Member_Eukaryota_2759", "GenomicGC",  as.integer(c(1,20)), profileId=2, plotRegression=TRUE, caption="Start-referenced", colorTrait="Member_Fungi_4751" )
## performOLSregression_profileRangeMean_withFilter( traits, "Member_Opisthokonta_33154", "GenomicGC",  as.integer(c(1,20)), profileId=2, plotRegression=TRUE, caption="Start-referenced", colorTrait="Member_Fungi_4751" )
## performOLSregression_profileRangeMean_withFilter( traits, "Member_Fungi_4751", "GenomicGC",  as.integer(c(1,20)), profileId=2, plotRegression=TRUE, caption="Start-referenced", colorTrait="Member_Ascomycota_4890" )

## performOLSregression_profileRangeMean_withFilter( traits, "Member_Fungi_4751", "GenomicGC",  as.integer(c(6,6)), profileId=2, plotRegression=TRUE, caption="Start-referenced", colorTrait="Member_Ascomycota_4890" )
## performOLSregression_profileRangeMean_withFilter( traits, "Member_Fungi_4751", "GenomicGC",  as.integer(c(13,13)), profileId=2, plotRegression=TRUE, caption="Start-referenced", colorTrait="Member_Ascomycota_4890" )


## dev.off()
## quit()


## ## performGLSregressionWithFilter( traits, tree, "Member_all_1",                  "Profile_2.1",   "Profile_2.16" )
## performGLSregressionWithFilter( traits, tree, "Member_all_1", "GenomicGC", "Profile_1.26" )
## performGLSregressionWithFilter( traits, tree, "Member_Eukaryota_2759", "GenomicGC", "Profile_1.11" )
## performGLSregressionWithFilter( traits, tree, "Member_Eukaryota_2759", "GenomicGC", "Profile_1.16" )
## performGLSregressionWithFilter( traits, tree, "Member_Eukaryota_2759", "GenomicGC", "Profile_1.21" )
## performGLSregressionWithFilter( traits, tree, "Member_Eukaryota_2759", "GenomicGC", "Profile_1.26" )
## performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2", "GenomicGC", "Profile_1.26" )

## performGLSregressionWithFilter( traits, tree, "Member_all_1", "GenomicGC", "Profile_2.6" )
## performGLSregressionWithFilter( traits, tree, "Member_Eukaryota_2759", "GenomicGC", "Profile_2.6" )
## performGLSregressionWithFilter( traits, tree, "Member_Eukaryota_2759", "GenomicGC", "Profile_2.11" )
## performGLSregressionWithFilter( traits, tree, "Member_Eukaryota_2759", "GenomicGC", "Profile_2.16" )
## performGLSregressionWithFilter( traits, tree, "Member_Eukaryota_2759", "GenomicGC", "Profile_2.21" )
## performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2", "GenomicGC", "Profile_2.6" )

## ##performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2", "GenomicGC",   "Profile_1.26" )

## performGLSregressionWithFilter( traits, tree, "Member_all_1", "GenomicENc.prime", "Profile_1.26" )
## performGLSregressionWithFilter( traits, tree, "Member_Eukaryota_2759", "GenomicENc.prime", "Profile_1.26" )
## performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2", "GenomicENc.prime", "Profile_1.26" )

## performGLSregressionWithFilter( traits, tree, "Member_all_1", "GenomicENc.prime", "Profile_2.6" )
## performGLSregressionWithFilter( traits, tree, "Member_Eukaryota_2759", "GenomicENc.prime", "Profile_2.6" )
## performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2", "GenomicENc.prime", "Profile_2.6" )



## dev.off()
## quit()


findValueOutliersWithFilter <- function( traits, filterTrait, criteria, pvalThreshold=0.01, logical="all", epsilon=0.1 )
{
    stopifnot( any( filterTrait == colnames(traits)  ) )  # filterTrait not found

    # Filter by the filter trait
    traits <- traits[ (traits[,filterTrait]>0), ]

    numSpecies <- nrow(traits)

    
    #print(nrow(traits))
    test.val <- sapply(1:length(criteria), function(i) sprintf("test.val.%d", i) )
    test.res <- sapply(1:length(criteria), function(i) sprintf("test.res.%d", i) )
    descr <- sapply(1:length(criteria), function(i) sprintf("profile%d (%d-%d) %s ", criteria[[i]][[1]], criteria[[i]][[2]], criteria[[i]][[3]], criteria[[i]][[4]]) )
    #print(test.val)
    #print(test.res)

    for( i in 1:length(criteria) )
    {
        criterion <- criteria[[i]]


        variables <- getProfileVariables( c(criterion[[2]], criterion[[3]]), profileId=criterion[[1]] )
        #print(variables)

        traits[,test.val[i]] <- rowMeans( traits[variables] )
                                        #traits$RangeMean <- rowMeans( traits[variables] )
        stopifnot( any( test.val[i] == colnames(traits)  ) )  # studyTrait  not found
        condition  <-  criterion[[4]]
        stopifnot(condition==">0" || condition=="<0")
        #ccol <- crits[i]

        total=nrow(traits)
        
        if( condition==">0" )
        {
            traits[,test.res[i]] <- traits[,test.val[i]] >  epsilon
        }
        else
        {
            traits[,test.res[i]] <- traits[,test.val[i]] <  -epsilon
        }
    }

    #print(dimnames(traits))
    #print(dim(traits[test.res]))
    #print(dim(traits[crits]))
    

    if( logical == "all")
    {
        matching <- apply( traits[test.res], MARGIN=1, FUN=all)
    }
    else
    {
        matching <- apply( traits[test.res], MARGIN=1, FUN=any)
    }

    #print(class(matching))
    #print(matching)
    #print(rownames(matching))

    posBDQA <- paste0( sample( names(matching[ matching]), size=min(10, length(matching[ matching])), replace=FALSE ), collapse=",")
    negBDQA <- paste0( sample( names(matching[!matching]), size=min(10, length(matching[!matching])), replace=FALSE ), collapse=",")
    

    return( data.frame( Filter=c(filterTrait), Descr=c(descr), Matches=c(sum(matching>0)), Percent=c(sum(matching>0)/total*100), Total=c(total), BDQA.pos=c(posBDQA), BDQA.neg=c(negBDQA) ) )
}


## print(sum(traits[!is.na(traits$TempOver75),]$TempOver75   + 0))
## print(sum(traits[!is.na(traits$TempUnder75),]$TempUnder75 + 0))

## #quit()


## #results <- glsRangeAnalysisWithFilter( traits, tree, "TempOver75",  "OptimumTemp", c(1,31),   profileId=1)
## #results <- glsRangeAnalysisWithFilter( traits, tree, "TempUnder75", "OptimumTemp", c(1,31),   profileId=1)

## #results <- glsRangeAnalysisWithFilter( traits, tree, "TempUnder75", "GenomicGC", c(2,32),   profileId=2)

#-----------------------------------------------------------------------------------------------------------------------------------
# Create outlier tallies (I eventually re-implemented this in python)
#-----------------------------------------------------------------------------------------------------------------------------------

## outliers <- data.frame( Filter=character(), Descr=character(), Matches=integer(), Percent=double(), Total=integer(), BDQA.pos=character(), BDQA.net=character() )


## #--1--
## outliers <- rbind( outliers, findValueOutliersWithFilter( traits, "Member_all_1",  list(list(1, 1, 3,  ">0")) ) )
## outliers <- rbind( outliers, findValueOutliersWithFilter( traits, "Member_all_1",  list(list(1, 1, 3,  "<0")) ) )

## outliers <- rbind( outliers, findValueOutliersWithFilter( traits, "Member_all_1",  list(list(1, 1, 1,  ">0")) ) )
## outliers <- rbind( outliers, findValueOutliersWithFilter( traits, "Member_all_1",  list(list(1, 2, 2,  ">0")) ) )
## outliers <- rbind( outliers, findValueOutliersWithFilter( traits, "Member_all_1",  list(list(1, 3, 3,  ">0")) ) )


## #--2--
## outliers <- rbind( outliers, findValueOutliersWithFilter( traits, "Member_all_1",  list(list(2, 32, 32, ">0")) ) )
## outliers <- rbind( outliers, findValueOutliersWithFilter( traits, "Member_all_1",  list(list(2, 32, 32, "<0")) ) )

## #--3--
## outliers <- rbind( outliers, findValueOutliersWithFilter( traits, "TempOver75",  list(list(1,  1,  3,  ">0")) ) )
## outliers <- rbind( outliers, findValueOutliersWithFilter( traits, "TempUnder75", list(list(1, 31, 31, ">0")) ) )

## #--4--
## outliers <- rbind( outliers, findValueOutliersWithFilter( traits, "Member_Bacteria_2",  list(list(1,  1,  3,  "<0")) ) )
## outliers <- rbind( outliers, findValueOutliersWithFilter( traits, "Member_Bacteria_2",  list(list(2, 32, 32, "<0")) ) )

## #--5--
## outliers <- rbind( outliers, findValueOutliersWithFilter( traits, "Member_all_1",  list( list(1, 1, 1, ">0") ), logical="all" ) )
## outliers <- rbind( outliers, findValueOutliersWithFilter( traits, "Member_all_1",  list( list(1, 15, 31, "<0") ), logical="all" ) )
## outliers <- rbind( outliers, findValueOutliersWithFilter( traits, "Member_all_1",  list(list(1, 1, 1, ">0"), list(1, 15, 31, "<0") ), logical="all" ) )
## outliers <- rbind( outliers, findValueOutliersWithFilter( traits, "Member_all_1",  list(list(1, 1, 1, "<0"), list(1, 15, 31, ">0") ), logical="any" ) )



## #--6--
## outliers <- rbind( outliers, findValueOutliersWithFilter( traits, "Member_Bacteria_2",  list(list(1, 1, 1, ">0"), list(1, 15, 31, "<0") ), logical="all" ) )
## outliers <- rbind( outliers, findValueOutliersWithFilter( traits, "Member_Bacteria_2",  list(list(1, 1, 1, "<0"), list(1, 15, 31, ">0") ), logical="any" ) )

## #--7--
## outliers <- rbind( outliers, findValueOutliersWithFilter( traits, "Member_Eukaryota_2759",  list(list(1, 1, 1, ">0"), list(1, 15, 31, "<0") ), logical="all" ) )
## outliers <- rbind( outliers, findValueOutliersWithFilter( traits, "Member_Eukaryota_2759",  list(list(1, 1, 1, "<0"), list(1, 15, 31, ">0") ), logical="any" ) )

## #--8--
## outliers <- rbind( outliers, findValueOutliersWithFilter( traits, "Member_Archaea_2157",  list(list(1, 1, 1, ">0"), list(1, 15, 31, "<0") ), logical="all" ) )
## outliers <- rbind( outliers, findValueOutliersWithFilter( traits, "Member_Archaea_2157",  list(list(1, 1, 1, "<0"), list(1, 15, 31, ">0") ), logical="any" ) )


## #--9--
## outliers <- rbind( outliers, findValueOutliersWithFilter( traits, "Member_all_1",          list(list(2, 2, 16, "<0"), list(2, 32, 32, ">0") ), logical="all" ) )
## outliers <- rbind( outliers, findValueOutliersWithFilter( traits, "Member_Bacteria_2",     list(list(2, 2, 16, "<0"), list(2, 32, 32, ">0") ), logical="all" ) )
## outliers <- rbind( outliers, findValueOutliersWithFilter( traits, "Member_Eukaryota_2759", list(list(2, 2, 16, "<0"), list(2, 32, 32, ">0") ), logical="all" ) )
## outliers <- rbind( outliers, findValueOutliersWithFilter( traits, "Member_Archaea_2157",   list(list(2, 2, 16, "<0"), list(2, 32, 32, ">0") ), logical="all" ) )


## print(outliers)
#-----------------------------------------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------------------------------------
# Find the common position of the "0-crossing" (the point at which the positive profile becomes negative
#-----------------------------------------------------------------------------------------------------------------------------------

## get.mode <- function(v)
## {
##     uniqv <- unique(v)
##     uniqv[which.max(tabulate(match(v, uniqv)))]
## }


plotPositiveStretchLengths <- function(traits, filterTrait="Member_all_1")
{
    # Filter species by trait
    traits.filtered <- traits[as.logical(traits[,filterTrait]),]

    stopifnot(profileMode=="dLFE")
    
    variables <- getProfileVariables( c(1, 31), profileId=1 )
    #print(variables)
    #print(traits[,variables])

    #print(".1.1.1.1.1.1.")
    values <- traits.filtered[,variables]
    #print(length(values))
    #print(values)
    #print(vals[5,])
    #x <- values[5:10]
    #print(x)
    #print(x<0)
    
    
    vals <- data.frame( v0=apply( values, MARGIN=1, FUN=function (x) min(c(1:31)[x<0]) ) )  # find first 0-crossing
    vals$v <- (vals$v0-1)*10
    vals$group <- filterTrait

    p <- ggplot( vals, aes(v) ) +
        geom_histogram( ) +
        geom_histogram( aes(y=cumsum(..count..)), alpha=0.3 ) +
        labs(title="0-crossing position (cummulative)", x="Position relative to CDS start (nt)", y="Number of species")
    print(p)

    p <- ggplot( vals, aes(v) ) +
        geom_histogram() +
        labs(title="0-crossing position", x="Position relative to CDS start (nt)", y="Number of species")
    print(p)
    
    
    return( vals)
}

figure_PositiveStretchLengths <- function()
{
    
    groupsForHistogram <- c("Member_all_1", "Member_Bacteria_2", "Member_Eukaryota_2759", "Member_Archaea_2157")
    vals <- data.frame()
    for( filter in groupsForHistogram )
    {
        v <- plotPositiveStretchLengths(traits, filterTrait=filter)
        print(v)
        vals <- rbind( vals, v)
    }

    vals[ vals$v==Inf, ]$v <- 350

    p <- ggplot( vals, aes(v, color=group) ) +
    #    geom_histogram( ) +
        geom_histogram( aes(y=cumsum(..count..)), alpha=0.3, position="identity" ) +
        labs(title="0-crossing position (cummulative)", x="Position relative to CDS start (nt)", y="Number of species")
    print(p)

    p <- ggplot( vals, aes(v, y=..density.., color=group) ) +
        geom_freqpoly(binwidth=10, position="identity", size=1.5, alpha=0.7) +
        scale_x_continuous( limits=c(0,300) ) +
        labs(title="0-crossing position frequency", x="Position relative to CDS start (nt)", y="Number of species") +
        theme( plot.background = element_blank(),   # Hide unnecessary theme elements (background panels, etc.)
               panel.grid.major.y = element_line(color="grey", size=0.50),
               panel.grid.major.x = element_line(color="grey", size=0.50),
               panel.grid.minor = element_blank(),
               panel.background = element_blank()
               )
    print(p)

    hists <- data.frame()
    histBreaks <- seq(0,360,10)
    L <- length(histBreaks)
    for( filter in groupsForHistogram )
    {
        print(filter)
        print(nrow(vals[vals$group==filter,]))
        h <- hist(vals[vals$group==filter,]$v, breaks=histBreaks, plot=FALSE)

        #print(h)
        #h$group <- filter
        hists <- rbind( hists, data.frame(x0=h$breaks[1:(L-1)], density=h$density*10, group=filter) )
    }

    #print(hists)

    p <- ggplot( hists, aes(x=x0, y=density, color=group) ) +
        geom_step(alpha=0.7, size=1.5) +
        scale_x_continuous( limits=c(0,350), breaks=c(0,100,200,300,340), labels=c("0", "100","200","300","No\ncrossing") ) +
        scale_color_brewer(type="seq", palette="PuBuGn") +
        labs(title="", x="Position relative to CDS start (nt)", y="Fraction of species") +
        theme( plot.background = element_blank(),   # Hide unnecessary theme elements (background panels, etc.)
               panel.grid.major.y = element_line(color="grey", size=0.50),
               panel.grid.major.x = element_line(color="grey", size=0.50),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.text.y = element_text( size=round(fontScale*1.3) ),
               axis.text.x = element_text( size=round(fontScale*1.5) ),
               axis.title  = element_text( size=round(fontScale*1.3) ),
               legend.text = element_text( size=round(fontScale*1.0) ),
               aspect.ratio=1.1
               )
    print(p)

    p <- ggplot( hists, aes(x=x0, y=density, color=group) ) +
        geom_step(alpha=0.7, size=1.5) +
        scale_x_continuous( limits=c(0,350), breaks=c(0,100,200,300,340), labels=c("0", "100","200","300","No\ncrossing") ) +
        scale_y_log10() +
        scale_color_brewer(type="seq", palette="PuBuGn") +
        labs(title="", x="Position relative to CDS start (nt)", y="Fraction of species") +
        theme( plot.background = element_blank(),   # Hide unnecessary theme elements (background panels, etc.)
               panel.grid.major.y = element_line(color="grey", size=0.50),
               panel.grid.major.x = element_line(color="grey", size=0.50),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.text.y = element_text( size=round(fontScale*1.3) ),
               axis.text.x = element_text( size=round(fontScale*1.5) ),
               axis.title  = element_text( size=round(fontScale*1.3) ),
               legend.text = element_text( size=round(fontScale*1.0) ),
               aspect.ratio=1.1
               )
    print(p)

    ##for( p in c("Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd"))
    ##{

    p <- ggplot( hists, aes(x=x0, y=density, color=group, size=group) ) +
        geom_step(alpha=0.4) +           # color "fringe"
        geom_step(alpha=0.8, size=1.1) + # color "core"
        scale_x_continuous( limits=c(0,350), breaks=c(0,100,200,300,340), labels=c("0", "100","200","300","No\ncrossing") ) +
        scale_y_sqrt() +
                                            #    scale_color_brewer(type="seq", palette="YlGnBu", direction=-1) +
        scale_size_manual( values=c("Member_all_1"=3.6, "Member_Bacteria_2"=2.4, "Member_Eukaryota_2759"=2.4, "Member_Archaea_2157"=2.4) ) +
        scale_color_manual(values=c("#ac3208", "#15be6c", "#6aa0de", "#591153")) +
        labs(title="", x="Position relative to CDS start (nt)", y="Fraction of species") +
        theme( plot.background = element_blank(),   # Hide unnecessary theme elements (background panels, etc.)
               panel.grid.major.y = element_line(color="grey", size=0.50),
               panel.grid.major.x = element_line(color="grey", size=0.50),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.text.y = element_text( size=round(fontScale*1.3) ),
               axis.text.x = element_text( size=round(fontScale*1.5) ),
               axis.title  = element_text( size=round(fontScale*1.3) ),
               legend.text = element_text( size=round(fontScale*1.0) ),
               aspect.ratio=1.1
               )
    print(p)
    ##}

    p <- ggplot( hists, aes(x=x0, y=density, fill=group) ) +
        geom_bar( stat="identity", position="dodge" ) +
        scale_x_continuous( limits=c(0,300) ) +
        scale_color_brewer(type="seq", palette="PuBuGn") +
        labs(title="", x="Position relative to CDS start (nt)", y="Fraction of species") +
        theme( plot.background = element_blank(),   # Hide unnecessary theme elements (background panels, etc.)
               panel.grid.major.y = element_line(color="grey", size=0.50),
               panel.grid.major.x = element_line(color="grey", size=0.50),
               panel.grid.minor = element_blank(),
               panel.background = element_blank()
               )
    print(p)




    p <- ggplot( vals, aes(v, color=group) ) +
        stat_ecdf( geom="step", size=1.2, alpha=0.7, pad=TRUE ) +
        scale_x_continuous( limits=c(0,310) ) +
        labs(title="0-crossing position (cummulative)", x="Position relative to CDS start (nt)", y="Cummulative frequency") +
        theme( plot.background = element_blank(),   # Hide unnecessary theme elements (background panels, etc.)
              panel.grid.major.y = element_line(color="grey", size=0.50),
              panel.grid.major.x = element_line(color="grey", size=0.50),
              panel.grid.minor = element_blank(),
              panel.background = element_blank()
              )
    print(p)
}


#-----------------------------------------------------------------------------------------------------------------------------------

############################################################################################################################################3

## for( gr in taxGroups )
## {
##     print(gr)
##     ret <- testPhylosignalWithFilter( traits, tree, gr, "Profile_1.15")
##     print( ret )
##     print(class(ret))
## }


#dev.off()
#quit()

# Why doesn't this have any effect?
#performGLSregression( traits, tree,   "GenomicGC",   "Profile_1.15", plotRegression=TRUE, caption="gamma=1.0",  traits.full=NA, bm.gamma=1.0 )
#performGLSregression( traits, tree,   "GenomicGC",   "Profile_1.15", plotRegression=TRUE, caption="gamma=0.8",  traits.full=NA, bm.gamma=0.8 )
#performGLSregression( traits, tree,   "GenomicGC",   "Profile_1.15", plotRegression=TRUE, caption="gamma=0.6",  traits.full=NA, bm.gamma=0.6 )
#performGLSregression( traits, tree,   "GenomicGC",   "Profile_1.15", plotRegression=TRUE, caption="gamma=0.2",  traits.full=NA, bm.gamma=0.2 )
#performGLSregression( traits, tree,   "GenomicGC",   "Profile_1.15", plotRegression=TRUE, caption="gamma=0.1",  traits.full=NA, bm.gamma=0.1 )
#performGLSregression( traits, tree,   "GenomicGC",   "Profile_1.15", plotRegression=TRUE, caption="gamma=0.05", traits.full=NA, bm.gamma=0.05 )

#quit()


###glsRangeAnalysisWithFilter( traits, tree, "Member_all_1", "GenomicGC", c(1,pyramidLength), profileId=1)
#glsRangeAnalysisWithFilter( traits, tree, "Member_all_1", "GenomicGC", c(2,32), profileId=2)

#dev.off()
#quit()


## performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",                  "OxygenReq",   "Profile_1.7" )
## performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",                  "OxygenReq",   "Profile_1.14" )
## performGLSregressionWithFilter( traits, tree, "Member_Terrabacteria_group_1783272", "OxygenReq",   "Profile_1.8" )

## #performGLSregressionWithFilter( traits, tree, "Member_Proteobacteria_1224",         "OxygenReq",   "Profile.16" )

## #performGLSregressionWithFilter( traits, tree, "Member_Gammaproteobacteria_1236",    "OxygenReq",   "Profile.19" )
## #performGLSregressionWithFilter( traits, tree, "Member_Gammaproteobacteria_1236",    "OxygenReq",   "Profile.20" )

## performGLSregressionWithFilter( traits, tree, "Member_Actinobacteria_1760", "OxygenReq",   "Profile_1.9" )
## performGLSregressionWithFilter( traits, tree, "Member_Actinobacteria_1760", "OxygenReq",   "Profile_1.14" )
#performGLSregressionWithFilter( traits, tree, "Member_Actinobacteria_1760", "OxygenReq",   "Profile.20" )
#performGLSregressionWithFilter( traits, tree, "Member_Terrabacteria_group_1783272", "OxygenReq", "Profile.16" )
#performGLSregressionWithFilter( traits, tree, "Member_Terrabacteria_group_1783272", "OxygenReq", "Profile.19" )

#performGLSregressionWithFilter( traits, tree, "Member_Terrabacteria_group_1783272", "GenomicGC", "Profile.15" )
#performGLSregressionWithFilter( traits, tree, "Member_Terrabacteria_group_1783272", "Habitat",   "Profile.15" )

#print("////-////-//-/")
#print(sum(!is.na(traits$GenomicGC)))
#print(sum(!is.na(traits$HighGC)))
#print(sum(!is.na(traits$Profile.1)))

## performGLSregressionWithFilter( traits, tree, "Member_all_1", "GenomicGC", "Profile_1.15" )


## #performGLSregressionWithFilter( traits, tree, "Member_Bacteroidetes_Chlorobi_group_68336", "GenomicGC", "Profile.22" )
## #performGLSregressionWithFilter( traits, tree, "Member_Bacteroidetes_Chlorobi_group_68336", "GenomicGC", "Profile.25" )

## #performGLSregressionWithFilter( traits, tree, "Member_Chloroflexi_200795", "GenomicGC", "Profile.18" )
## #performGLSregressionWithFilter( traits, tree, "Member_Chloroflexi_200795", "GenomicGC", "Profile.26" )
## #performGLSregressionWithFilter( traits, tree, "Member_Chloroflexi_200795", "GenomicGC", "Profile.27" )

## performGLSregressionWithFilter( traits, tree, "Member_Fungi_4751", "GenomicGC", "Profile_1.8" )
## #performGLSregressionWithFilter( traits, tree, "Member_Fungi_4751", "GenomicGC", "Profile.20" )

## performGLSregressionWithFilter( traits, tree, "Member_Proteobacteria_1224",         "OptimumTemp",   "Profile_1.10" )
## #performGLSregressionWithFilter( traits, tree, "Member_Proteobacteria_1224",         "OptimumTemp",   "Profile.17" )

## performGLSregressionWithFilter( traits, tree, "Member_Archaea_2157",         "OptimumTemp",   "Profile_1.11" )
## #performGLSregressionWithFilter( traits, tree, "Member_Archaea_2157",         "OptimumTemp",   "Profile.20" )

## performGLSregressionWithFilter( traits, tree, "Member_all_1",                "OptimumTemp",   "Profile_1.11" )
## #performGLSregressionWithFilter( traits, tree, "Member_all_1",                "OptimumTemp",   "Profile.20" )


## performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",                  "Habitat",   "Profile_1.1" )
## performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",                  "Habitat",   "Profile_1.2" )
## performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",                  "Habitat",   "Profile_1.3" )
## performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",                  "Habitat",   "Profile_1.4" )
## performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",                  "Habitat",   "Profile_1.5" )
## performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",                  "Habitat",   "Profile_1.6" )
#performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",                  "Habitat",   "Profile_1.27" )


####glsRangeAnalysisWithFilter( traits, tree, "Member_All_1", "GenomicGC", c(1,pyramidLength), xxxxxxxxxxx)


# TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY #
#dev.off()
#quit()
# TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY #


figure_DLFEInteractingTraits_RegressionRangeAnalysisByTaxGroup <- function()
{
    
    regressionResultsByTaxGroup <- data.frame( ExplanatoryVar=character(), Range=integer(), MaxRangeStart=integer(), MaxRangeEnd=integer(), TaxGroup=integer(), TaxGroupName=character(), EffectSize=double(), Pvalue=double(), NumSpecies=integer() )


    ## for( gr in taxGroups )
    ## {
    ##     print(gr)
    ##     results <- glsRangeAnalysisWithFilter( traits, tree, gr, "GenomicGC", c(1,pyramidLength))

    ##     # Store peaks (for each range) for this group
    ##     if( any( class(results) == "data.frame" ) && nrow(results) )
    ##     {
    ##         results <- results[results$Var1==results$Var2,]    # Only include single-window results
    ##         results$Pvalue.adj <- p.adjust( results$Pvalue, method="fdr" )

    ##         print(results)

    ##         # Iterate over each range to find the relevant peak
    ##         for( i in 1:length(groupsTableOutputFile.limitRangeFromNt) )
    ##         {
    ##             matching <- (results$Var1 >= groupsTableOutputFile.limitRangeFromNt[i]) &
    ##                         (results$Var2 <= groupsTableOutputFile.limitRangeToNt[i]  )  # Only include ranges within the configured limits

    ##             maxResult <- results[matching,][which.max( abs(results[matching, "Buse.R2"]) ),]  # Choose the result with the highest R^2 
    ##             medianR2  <- median( results[matching, "Buse.R2"] )

    ##             regressionResultsByTaxGroup <- rbind( regressionResultsByTaxGroup, data.frame( ExplanatoryVar=c("GenomicGC"),  Range=c(i-1), MaxRangeStart=c(maxResult$Var1), MaxRangeEnd=c(maxResult$Var2), TaxGroup=c(taxGroupToTaxId[gr]), TaxGroupName=c(gr), EffectSize=c(medianR2), Pvalue=c(maxResult$Pvalue.adj), NumSpecies=c(maxResult$NumSpecies) ) )
    ##         }
    ##     }
    ## }

    ## for( gr in taxGroups )
    ## {
    ##     print(gr)
    ##     results <- glsRangeAnalysisWithFilter( traits, tree, gr, "GenomicENc.prime", c(1,pyramidLength))

    ##     # Store peaks (for each range) for this group
    ##     if( any( class(results) == "data.frame" ) && nrow(results) )
    ##     {
    ##         results <- results[results$Var1==results$Var2,]    # Only include single-window results
    ##         results$Pvalue.adj <- p.adjust( results$Pvalue, method="fdr" )

    ##         print(results)

    ##         # Iterate over each range to find the relevant peak
    ##         for( i in 1:length(groupsTableOutputFile.limitRangeFromNt) )
    ##         {
    ##             matching <- (results$Var1 >= groupsTableOutputFile.limitRangeFromNt[i]) &
    ##                         (results$Var2 <= groupsTableOutputFile.limitRangeToNt[i]  )  # Only include ranges within the configured limits

    ##             maxResult <- results[matching,][which.max( abs(results[matching, "Buse.R2"]) ),]  # Choose the result with the highest R^2 
    ##             medianR2  <- median( results[matching, "Buse.R2"] )

    ##             regressionResultsByTaxGroup <- rbind( regressionResultsByTaxGroup, data.frame( ExplanatoryVar=c("GenomicENc.prime"),  Range=c(i-1), MaxRangeStart=c(maxResult$Var1), MaxRangeEnd=c(maxResult$Var2), TaxGroup=c(taxGroupToTaxId[gr]), TaxGroupName=c(gr), EffectSize=c(medianR2), Pvalue=c(maxResult$Pvalue.adj), NumSpecies=c(maxResult$NumSpecies) ) )
    ##         }
    ##     }
    ## }


    # Plots dependence of ENc and ENc' on GC%
    #p <- ggplot(
    #    data=data.frame(melt(traits, measure.vars=c("GenomicENc", "GenomicENc.prime"))),
    #    aes(x=GenomicGC, y=value, color=variable) ) + geom_point();
                                            #print(p)


    #dev.off()
    #quit()

    ## for( gr in taxGroups )
    ## {
    ##     print(gr)
    ##     results <- glsRangeAnalysisWithFilter( traits, tree, gr, "GenomicENc", c(1,pyramidLength))

    ##     # Store peaks (for each range) for this group
    ##     if( any( class(results) == "data.frame" ) && nrow(results) )
    ##     {
    ##         results <- results[results$Var1==results$Var2,]    # Only include single-window results
    ##         results$Pvalue.adj <- p.adjust( results$Pvalue, method="fdr" )

    ##         print(results)

    ##         # Iterate over each range to find the relevant peak
    ##         for( i in 1:length(groupsTableOutputFile.limitRangeFromNt
    ##         {
    ##             matching <- (results$Var1 >= groupsTableOutputFile.limitRangeFromNt[i]) &
    ##                         (results$Var2 <= groupsTableOutputFile.limitRangeToNt[i]  )  # Only include ranges within the configured limits

    ##             maxResult <- results[matching,][which.max( abs(results[matching, "Buse.R2"]) ),]  # Choose the result with the highest R^2 
    ##             medianR2  <- median( results[matching, "Buse.R2"] )

    ##             regressionResultsByTaxGroup <- rbind( regressionResultsByTaxGroup, data.frame( ExplanatoryVar=c("GenomicENc"),  Range=c(i-1), MaxRangeStart=c(maxResult$Var1), MaxRangeEnd=c(maxResult$Var2), TaxGroup=c(taxGroupToTaxId[gr]), TaxGroupName=c(gr), EffectSize=c(medianR2), Pvalue=c(maxResult$Pvalue.adj), NumSpecies=c(maxResult$NumSpecies) ) )
    ##         }
    ##     }
    ## }

    # Plots dependence of DCBS on ENc, etc.
    ## p <- ggplot(
    ##     data=data.frame(melt(traits, measure.vars=c("GenomicENc", "GenomicENc.prime", "GenomicDCBS") ) ),
    ##     aes(x=GenomicGC, y=value, color=variable) ) + geom_point();
    ## print(p)

    ## p <- ggplot(
    ##     data=data.frame(melt(traits, measure.vars=c("GenomicENc", "GenomicENc.prime") ) ),
    ##     aes(x=GenomicDCBS, y=value, color=variable) ) + geom_point();
    ## print(p)


    ## for( gr in taxGroups )
    ## {
    ##     print(gr)
    ##     results <- glsRangeAnalysisWithFilter( traits, tree, gr, "GenomicDCBS", c(1,pyramidLength))

    ##     # Store peaks (for each range) for this group
    ##     if( any( class(results) == "data.frame" ) && nrow(results) )
    ##     {
    ##         results <- results[results$Var1==results$Var2,]    # Only include single-window results
    ##         results$Pvalue.adj <- p.adjust( results$Pvalue, method="fdr" )

    ##         print(results)

    ##         # Iterate over each range to find the relevant peak
    ##         for( i in 1:length(groupsTableOutputFile.limitRangeFromNt) )
    ##         {
    ##             matching <- (results$Var1 >= groupsTableOutputFile.limitRangeFromNt[i]) &
    ##                         (results$Var2 <= groupsTableOutputFile.limitRangeToNt[i]  )  # Only include ranges within the configured limits

    ##             maxResult <- results[matching,][which.max( abs(results[matching, "Buse.R2"]) ),]  # Choose the result with the highest R^2 
    ##             medianR2  <- median( results[matching, "Buse.R2"] )

    ##             regressionResultsByTaxGroup <- rbind( regressionResultsByTaxGroup, data.frame( ExplanatoryVar=c("GenomicDCBS"),  Range=c(i-1), MaxRangeStart=c(maxResult$Var1), MaxRangeEnd=c(maxResult$Var2), TaxGroup=c(taxGroupToTaxId[gr]), TaxGroupName=c(gr), EffectSize=c(medianR2), Pvalue=c(maxResult$Pvalue.adj), NumSpecies=c(maxResult$NumSpecies) ) )
    ##         }
    ##     }
    ## }


    # TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY #
    ##write.csv(regressionResultsByTaxGroup, file=groupsTableOutputFile )
    ##dev.off()
    ##quit()
    # TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY #

    ## for( gr in taxGroups )
    ## {
    ##     print(gr)
    ##     results <- glsRangeAnalysisWithFilter( traits, tree, gr, "OptimumTemp", c(1,pyramidLength))

    ##     # Store peaks (for each range) for this group
    ##     if( any( class(results) == "data.frame" ) && nrow(results) )
    ##     {
    ##         results <- results[results$Var1==results$Var2,]    # Only include single-window results

    ##         # Iterate over each range to find the relevant peak
    ##         for( i in 1:length(groupsTableOutputFile.limitRangeFromNt) )
    ##         {
    ##             matching <- (results$Var1 >= groupsTableOutputFile.limitRangeFromNt[i]) &
    ##                         (results$Var2 <= groupsTableOutputFile.limitRangeToNt[i]  )  # Only include ranges within the configured limits

    ##             maxResult <- results[matching,][which.max( abs(results[matching, "Buse.R2"]) ),]  # Choose the result with the highest R^2 

    ##             regressionResultsByTaxGroup <- rbind( regressionResultsByTaxGroup, data.frame( ExplanatoryVar=c("OptimumTemp"), Range=c(i-1), MaxRangeStart=c(maxResult$Var1), MaxRangeEnd=c(maxResult$Var2), TaxGroup=c(taxGroupToTaxId[gr]), TaxGroupName=c(gr), EffectSize=c(maxResult$Buse.R2), Pvalue=c(maxResult$Pvalue), NumSpecies=c(maxResult$NumSpecies) ) )
    ##         }
    ##     }
    ## }


    ## for( gr in taxGroups )
    ## {
    ##     print(gr)
    ##     results <- glsRangeAnalysisWithFilter( traits, tree, gr, "TempHighLow75", c(1,pyramidLength))

    ##     # Store peaks (for each range) for this group
    ##     if( any( class(results) == "data.frame" ) && nrow(results) )
    ##     {
    ##         results <- results[results$Var1==results$Var2,]    # Only include single-window results

    ##         # Iterate over each range to find the relevant peak
    ##         for( i in 1:length(groupsTableOutputFile.limitRangeFromNt) )
    ##         {
    ##             matching <- (results$Var1 >= groupsTableOutputFile.limitRangeFromNt[i]) &
    ##                         (results$Var2 <= groupsTableOutputFile.limitRangeToNt[i]  )  # Only include ranges within the configured limits

    ##             maxResult <- results[matching,][which.max( abs(results[matching, "Buse.R2"]) ),]  # Choose the result with the highest R^2 

    ##             regressionResultsByTaxGroup <- rbind( regressionResultsByTaxGroup, data.frame( ExplanatoryVar=c("TempHighLow75"), Range=c(i-1), MaxRangeStart=c(maxResult$Var1), MaxRangeEnd=c(maxResult$Var2), TaxGroup=c(taxGroupToTaxId[gr]), TaxGroupName=c(gr), EffectSize=c(maxResult$Buse.R2), Pvalue=c(maxResult$Pvalue), NumSpecies=c(maxResult$NumSpecies) ) )
    ##         }
    ##     }
    ## }




    # TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY #
    #dev.off()  # TESTING ONLY
    #write.csv(regressionResultsByTaxGroup, file=groupsTableOutputFile )    # TESTING ONLY
    #quit()    # TESTING ONLY
    # TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY #

    ## for( gr in taxGroups )
    ## {
    ##     print(gr)
    ##     results <- glsRangeAnalysisWithFilter( traits, tree, gr, "Habitat", c(1,pyramidLength))

    ##     # Store peaks (for each range) for this group
    ##     if( any( class(results) == "data.frame" ) && nrow(results) )
    ##     {
    ##         results <- results[results$Var1==results$Var2,]    # Only include single-window results

    ##         # Iterate over each range to find the relevant peak
    ##         for( i in 1:length(groupsTableOutputFile.limitRangeFromNt) )
    ##         {
    ##             matching <- (results$Var1 >= groupsTableOutputFile.limitRangeFromNt[i]) &
    ##                         (results$Var2 <= groupsTableOutputFile.limitRangeToNt[i]  )  # Only include ranges within the configured limits

    ##             maxResult <- results[matching,][which.max( abs(results[matching, "Buse.R2"]) ),]  # Choose the result with the highest R^2 

    ##             regressionResultsByTaxGroup <- rbind( regressionResultsByTaxGroup, data.frame( ExplanatoryVar=c("Habitat"),    Range=c(i-1), MaxRangeStart=c(maxResult$Var1), MaxRangeEnd=c(maxResult$Var2), TaxGroup=c(taxGroupToTaxId[gr]), TaxGroupName=c(gr), EffectSize=c(maxResult$Buse.R2), Pvalue=c(maxResult$Pvalue), NumSpecies=c(maxResult$NumSpecies) ) )
    ##         }
    ##     }
    ## }

    ## for( gr in taxGroups )
    ## {
    ##     print(gr)
    ##     results <- glsRangeAnalysisWithFilter( traits, tree, gr, "Salinity", c(1,pyramidLength) )

    ##     # Store peaks (for each range) for this group
    ##     if( any( class(results) == "data.frame" ) && nrow(results) )
    ##     {
    ##         results <- results[results$Var1==results$Var2,]    # Only include single-window results

    ##         # Iterate over each range to find the relevant peak
    ##         for( i in 1:length(groupsTableOutputFile.limitRangeFromNt) )
    ##         {
    ##             matching <- (results$Var1 >= groupsTableOutputFile.limitRangeFromNt[i]) &
    ##                         (results$Var2 <= groupsTableOutputFile.limitRangeToNt[i]  )  # Only include ranges within the configured limits

    ##             maxResult <- results[matching,][which.max( abs(results[matching, "Buse.R2"]) ),]  # Choose the result with the highest R-squared 

    ##             regressionResultsByTaxGroup <- rbind( regressionResultsByTaxGroup, data.frame( ExplanatoryVar=c("Salinity"),   Range=c(i-1), MaxRangeStart=c(maxResult$Var1), MaxRangeEnd=c(maxResult$Var2), TaxGroup=c(taxGroupToTaxId[gr]), TaxGroupName=c(gr), EffectSize=c(maxResult$Buse.R2), Pvalue=c(maxResult$Pvalue), NumSpecies=c(maxResult$NumSpecies) ) )
    ##         }
    ##     }
    ## }

    ## for( gr in taxGroups )
    ## {
    ##     print(gr)
    ##     results <- glsRangeAnalysisWithFilter( traits, tree, gr, "OxygenReq", c(1,pyramidLength))

    ##     # Store peaks (for each range) for this group
    ##     if( any( class(results) == "data.frame" ) && nrow(results) )
    ##     {
    ##         results <- results[results$Var1==results$Var2,]    # Only include single-window results

    ##         # Iterate over each range to find the relevant peak
    ##         for( i in 1:length(groupsTableOutputFile.limitRangeFromNt) )
    ##         {
    ##             matching <- (results$Var1 >= groupsTableOutputFile.limitRangeFromNt[i]) &
    ##                         (results$Var2 <= groupsTableOutputFile.limitRangeToNt[i]  )  # Only include ranges within the configured limits

    ##             maxResult <- results[matching,][which.max( abs(results[matching, "Buse.R2"]) ),]  # Choose the result with the highest R^2 

    ##             regressionResultsByTaxGroup <- rbind( regressionResultsByTaxGroup, data.frame( ExplanatoryVar=c("OxygenReq"),  Range=c(i-1), MaxRangeStart=c(maxResult$Var1), MaxRangeEnd=c(maxResult$Var2), TaxGroup=c(taxGroupToTaxId[gr]), TaxGroupName=c(gr), EffectSize=c(maxResult$Buse.R2), Pvalue=c(maxResult$Pvalue), NumSpecies=c(maxResult$NumSpecies) ) )
    ##         }
    ##     }
    ## }



    ## results <- glsRangeAnalysisWithFilter( traits, tree, "Member_Bacteroidetes_Chlorobi_group_68336", "LogGenomeSize", c(1,pyramidLength), profileId=1)

    ## results <- glsRangeAnalysisWithFilter( traits, tree, "Member_Bacteroidetes_Chlorobi_group_68336", "LogGenomeSize", c(1,pyramidLength), profileId=2)

    ##dev.off()
    ## quit()

    ## results <- glsRangeAnalysisWithFilter( traits, tree, "Member_Bacteria_2", "LogGenomeSize", c(1,pyramidLength), profileId=1)
    ## results <- glsRangeAnalysisWithFilter( traits, tree, "Member_Proteobacteria_1224", "LogGenomeSize", c(1,pyramidLength), profileId=1)
    ## results <- glsRangeAnalysisWithFilter( traits, tree, "Member_Gammaproteobacteria_1236", "LogGenomeSize", c(1,pyramidLength), profileId=1)
    ## results <- glsRangeAnalysisWithFilter( traits, tree, "Member_Bacteria_2", "LogGenomeSize", c(1,pyramidLength), profileId=1)

    ## performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_Bacteria_2", "LogGenomeSize",  as.integer(c(12,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced" )
    ## performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_Bacteria_2", "LogGenomeSize",  as.integer(c(2,21)), profileId=2, plotRegression=TRUE, caption="End-referenced" )

    ## performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_Proteobacteria_1224", "LogGenomeSize",  as.integer(c(12,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced" )
    ## performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_Proteobacteria_1224", "LogGenomeSize",  as.integer(c(2,21)), profileId=2, plotRegression=TRUE, caption="End-referenced" )

    ## performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_Terrabacteria_group_1783272", "LogGenomeSize",  as.integer(c(12,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced" )
    ## performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_Terrabacteria_group_1783272", "LogGenomeSize",  as.integer(c(2,21)), profileId=2, plotRegression=TRUE, caption="End-referenced" )

    ## #--
    ## performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_Bacteria_2", "LogGrowthTime",  as.integer(c(12,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced" )
    ## performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_Bacteria_2", "LogGrowthTime",  as.integer(c( 2,21)), profileId=2, plotRegression=TRUE, caption="End-referenced" )

    ## performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_Proteobacteria_1224", "LogGrowthTime",  as.integer(c(12,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced" )
    ## performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_Proteobacteria_1224", "LogGrowthTime",  as.integer(c( 2,21)), profileId=2, plotRegression=TRUE, caption="End-referenced" )

    ## performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_Terrabacteria_group_1783272", "LogGrowthTime",  as.integer(c(12,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced" )
    ## performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_Terrabacteria_group_1783272", "LogGrowthTime",  as.integer(c(2,21)), profileId=2, plotRegression=TRUE, caption="End-referenced" )

    ##dev.off()
    ##quit()

    ## performGLSregressionWithFilter( traits, tree, "Member_all_1",                  "GenomeSizeMb",   "Profile_1.1" )
    ## performGLSregressionWithFilter( traits, tree, "Member_all_1",                  "LogGenomeSize",  "Profile_1.1" )
    ## performGLSregressionWithFilter( traits, tree, "Member_all_1",                  "ProteinCount",   "Profile_1.1" )

    ## performGLSregressionWithFilter( traits, tree, "Member_all_1",                  "GenomeSizeMb",   "Profile_1.16" )
    ## performGLSregressionWithFilter( traits, tree, "Member_all_1",                  "LogGenomeSize",  "Profile_1.16" )
    ## performGLSregressionWithFilter( traits, tree, "Member_all_1",                  "ProteinCount",   "Profile_1.16" )

    ## performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",                  "GenomeSizeMb",   "Profile_1.1" )
    ## performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",                  "LogGenomeSize",  "Profile_1.1" )
    ## performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",                  "ProteinCount",   "Profile_1.1" )
    ## performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",                  "GenomeSizeMb",   "Profile_1.6" )
    ## performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",                  "LogGenomeSize",  "Profile_1.6" )
    ## performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",                  "ProteinCount",   "Profile_1.6" )

    ## performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",                  "GenomeSizeMb",   "Profile_1.16" )
    ## performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",                  "LogGenomeSize",  "Profile_1.16" )
    ## performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",                  "ProteinCount",   "Profile_1.16" )

    ## performGLSregressionWithFilter( traits, tree, "Member_Terrabacteria_group_1783272",   "GenomeSizeMb",   "Profile_1.4" )
    ## performGLSregressionWithFilter( traits, tree, "Member_Terrabacteria_group_1783272",   "LogGenomeSize",  "Profile_1.4" )
    ## performGLSregressionWithFilter( traits, tree, "Member_Terrabacteria_group_1783272",   "ProteinCount",   "Profile_1.4" )

    ## performGLSregressionWithFilter( traits, tree, "Member_Proteobacteria_1224",   "GenomeSizeMb",   "Profile_1.13" )
    ## performGLSregressionWithFilter( traits, tree, "Member_Proteobacteria_1224",   "LogGenomeSize",  "Profile_1.13" )
    ## performGLSregressionWithFilter( traits, tree, "Member_Proteobacteria_1224",   "ProteinCount",   "Profile_1.13" )

    ## performGLSregressionWithFilter( traits, tree, "Member_Gammaproteobacteria_1236",   "GenomeSizeMb",   "Profile_1.13" )
    ## performGLSregressionWithFilter( traits, tree, "Member_Gammaproteobacteria_1236",   "LogGenomeSize",  "Profile_1.13" )
    ## performGLSregressionWithFilter( traits, tree, "Member_Gammaproteobacteria_1236",   "ProteinCount",   "Profile_1.13" )

    ## performGLSregressionWithFilter( traits, tree, "Member_Bacteroidetes_Chlorobi_group_68336",   "GenomeSizeMb",   "Profile_1.17" )
    ## performGLSregressionWithFilter( traits, tree, "Member_Bacteroidetes_Chlorobi_group_68336",   "LogGenomeSize",  "Profile_1.17" )
    ## performGLSregressionWithFilter( traits, tree, "Member_Bacteroidetes_Chlorobi_group_68336",   "ProteinCount",   "Profile_1.17" )

    ##results <- glsRangeAnalysisWithFilter( traits, tree, "Member_Bacteria_2", "LogGenomeSize", c(1,pyramidLength), profileId=1)
    ##results <- glsRangeAnalysisWithFilter( traits, tree, "Member_Proteobacteria_1224", "LogGenomeSize", c(1,pyramidLength), profileId=1)
    ##results <- glsRangeAnalysisWithFilter( traits, tree, "Member_Eukaryota_2759", "LogGenomeSize", c(1,pyramidLength), profileId=1)



    ##dev.off()
    ##quit()

    ## for( gr in taxGroups )
    ## {
    ##     print(gr)
    ##     results <- glsRangeAnalysisWithFilter( traits, tree, gr, "GenomeSizeMb", c(1,pyramidLength), profileId=1)

    ##     ##performGLSregressionWithFilter( traits, tree, gr,   "GenomeSizeMb",   "Profile_1.21" )

    ##     # Store peaks (for each range) for this group
    ##     if( any( class(results) == "data.frame" ) && nrow(results) )
    ##     {
    ##         results <- results[results$Var1==results$Var2,]    # Only include single-window results

    ##         # Iterate over each range to find the relevant peak
    ##         for( i in 1:length(groupsTableOutputFile.limitRangeFromNt) )
    ##         {
    ##             matching <- (results$Var1 >= groupsTableOutputFile.limitRangeFromNt[i]) &
    ##                         (results$Var2 <= groupsTableOutputFile.limitRangeToNt[i]  )  # Only include ranges within the configured limits

    ##             maxResult <- results[matching,][which.max( abs(results[matching, "Buse.R2"]) ),]  # Choose the result with the highest R^2 

    ##             regressionResultsByTaxGroup <- rbind( regressionResultsByTaxGroup, data.frame( ExplanatoryVar=c("GenomeSizeMb"),  Range=c(i-1), MaxRangeStart=c(maxResult$Var1), MaxRangeEnd=c(maxResult$Var2), TaxGroup=c(taxGroupToTaxId[gr]), TaxGroupName=c(gr), EffectSize=c(maxResult$Buse.R2), Pvalue=c(maxResult$Pvalue), NumSpecies=c(maxResult$NumSpecies) ) )
    ##         }
    ##     }
    ## }


    performOLSregression( traits, "LogGrowthTime", "GenomicENc.prime" )
    performOLSregression( traits, "LogGrowthTime", "GenomicGC" )
    performOLSregression( traits, "LogGrowthTime", "GenomeSizeMb" )
    performOLSregression( traits, "LogGrowthTime", "LogGenomeSize" )
    performOLSregression( traits, "LogGrowthTime", "OptimumTemp" )


    ## performOLSregression_profileRangeMean_withFilter( traits, "Member_all_1", "GenomicGC",  as.integer(c(12,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced", colorTrait="Member_Eukaryota_2759" )

    ## performOLSregression_profileRangeMean_withFilter( traits, "Member_Bacteria_2", "LogGrowthTime",  as.integer(c(12,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced", colorTrait="Member_Bacteria_2" )
    ## performOLSregression_profileRangeMean_withFilter( traits, "Member_Bacteria_2", "LogGrowthTime",  as.integer(c(1,2)), profileId=1, plotRegression=TRUE, caption="Start-referenced", colorTrait="Member_Bacteria_2" )
    ## performOLSregression_profileRangeMean_withFilter( traits, "Member_Proteobacteria_1224", "LogGrowthTime", as.integer(c(12,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced", colorTrait="Member_Proteobacteria_1224" )
    ## performOLSregression_profileRangeMean_withFilter( traits, "Member_Proteobacteria_1224", "LogGrowthTime", as.integer(c(1,2)), profileId=1, plotRegression=TRUE, caption="Start-referenced", colorTrait="Member_Proteobacteria_1224" )
    ## performOLSregression_profileRangeMean_withFilter( traits, "Member_Proteobacteria_1224", "LogGrowthTime", as.integer(c(2,21)), profileId=2, plotRegression=TRUE, caption="End-referenced", colorTrait="Member_Proteobacteria_1224" )

    ## performOLSregression_profileRangeMean_withFilter( traits, "Member_Gammaproteobacteria_1236", "LogGrowthTime", as.integer(c(12,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced", colorTrait="Member_Gammaproteobacteria_1236" )
    ## performOLSregression_profileRangeMean_withFilter( traits, "Member_Gammaproteobacteria_1236", "LogGrowthTime", as.integer(c(1,2)), profileId=1, plotRegression=TRUE, caption="Start-referenced", colorTrait="Member_Gammaproteobacteria_1236" )
    ## performOLSregression_profileRangeMean_withFilter( traits, "Member_Terrabacteria_group_1783272", "LogGrowthTime",  as.integer(c(12,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced", colorTrait="Member_Terrabacteria_group_1783272" )

    ## performOLSregression_profileRangeMean_withFilter( traits, "GC.45", "LogGrowthTime",  as.integer(c(12,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced", colorTrait="GC.45" )

    ## performOLSregression_profileRangeMean_withFilter( traits, "GC.45", "LogGenomeSize",  as.integer(c(12,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced", colorTrait="GC.45" )


    ## performOLSregression_profileRangeMean_withFilter( traits, "Member_all_1", "LogGenomeSize",  as.integer(c(11,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced", colorTrait="Member_Bacteria_2" )
    ## performOLSregression_profileRangeMean_withFilter( traits, "Member_all_1", "LogGenomeSize",  as.integer(c(1,2)), profileId=1, plotRegression=TRUE, caption="Start-referenced", colorTrait="Member_Bacteria_2" )
    ## performOLSregression_profileRangeMean_withFilter( traits, "Member_Bacteria_2", "LogGenomeSize",  as.integer(c(11,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced", colorTrait="Member_Bacteria_2" )
    ## performOLSregression_profileRangeMean_withFilter( traits, "Member_Bacteria_2", "LogGenomeSize",  as.integer(c(1,2)), profileId=1, plotRegression=TRUE, caption="Start-referenced", colorTrait="Member_Bacteria_2" )
    ## performOLSregression_profileRangeMean_withFilter( traits, "Member_Proteobacteria_1224", "LogGenomeSize", as.integer(c(11,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced", colorTrait="Member_Proteobacteria_1224" )
    ## performOLSregression_profileRangeMean_withFilter( traits, "Member_Proteobacteria_1224", "LogGenomeSize", as.integer(c(1,2)), profileId=1, plotRegression=TRUE, caption="Start-referenced", colorTrait="Member_Proteobacteria_1224" )
    ## performOLSregression_profileRangeMean_withFilter( traits, "Member_Proteobacteria_1224", "LogGenomeSize", as.integer(c(2,22)), profileId=2, plotRegression=TRUE, caption="End-referenced", colorTrait="Member_Proteobacteria_1224" )


    for( gr in taxGroups )
    {
        print(gr)
        results <- glsRangeAnalysisWithFilter( traits, tree, gr, "LogGenomeSize", c(1,pyramidLength), profileId=1)

        # Store peaks (for each range) for this group
        if( any( class(results) == "data.frame" ) && nrow(results) )
        {
            results <- results[results$Var1==results$Var2,]    # Only include single-window results

            # Iterate over each range to find the relevant peak
            for( i in 1:length(groupsTableOutputFile.limitRangeFromNt) )
            {
                matching <- (results$Var1 >= groupsTableOutputFile.limitRangeFromNt[i]) &
                            (results$Var2 <= groupsTableOutputFile.limitRangeToNt[i]  )  # Only include ranges within the configured limits

                maxResult <- results[matching,][which.max( abs(results[matching, "Buse.R2"]) ),]  # Choose the result with the highest R^2 

                regressionResultsByTaxGroup <- rbind( regressionResultsByTaxGroup, data.frame( ExplanatoryVar=c("GenomeSizeMb"),  Range=c(i-1), MaxRangeStart=c(maxResult$Var1), MaxRangeEnd=c(maxResult$Var2), TaxGroup=c(taxGroupToTaxId[gr]), TaxGroupName=c(gr), EffectSize=c(maxResult$Buse.R2), Pvalue=c(maxResult$Pvalue), NumSpecies=c(maxResult$NumSpecies) ) )
            }
        }
    }

    #return (NA);

    ## for( gr in taxGroups )
    ## {
    ##     print(gr)
    ##     results <- glsRangeAnalysisWithFilter( traits, tree, gr, "ProteinCount", c(1,pyramidLength), profileId=1)

    ##     # Store peaks (for each range) for this group
    ##     if( any( class(results) == "data.frame" ) && nrow(results) )
    ##     {
    ##         results <- results[results$Var1==results$Var2,]    # Only include single-window results

    ##         # Iterate over each range to find the relevant peak
    ##         for( i in 1:length(groupsTableOutputFile.limitRangeFromNt) )
    ##         {
    ##             matching <- (results$Var1 >= groupsTableOutputFile.limitRangeFromNt[i]) &
    ##                         (results$Var2 <= groupsTableOutputFile.limitRangeToNt[i]  )  # Only include ranges within the configured limits

    ##             maxResult <- results[matching,][which.max( abs(results[matching, "Buse.R2"]) ),]  # Choose the result with the highest R^2 

    ##             regressionResultsByTaxGroup <- rbind( regressionResultsByTaxGroup, data.frame( ExplanatoryVar=c("ProteinCount"),  Range=c(i-1), MaxRangeStart=c(maxResult$Var1), MaxRangeEnd=c(maxResult$Var2), TaxGroup=c(taxGroupToTaxId[gr]), TaxGroupName=c(gr), EffectSize=c(maxResult$Buse.R2), Pvalue=c(maxResult$Pvalue), NumSpecies=c(maxResult$NumSpecies) ) )
    ##         }
    ##     }
    ## }

    #return (NA);

    #glsRangeAnalysisWithFilter( traits, tree, "HighGC",    "OptimumTemp", c(1,pyramidLength) )
    #glsRangeAnalysisWithFilter( traits, tree, "LowGC",     "OptimumTemp", c(1,pyramidLength) )
    #glsRangeAnalysisWithFilter( traits, tree, "GC.0.40",   "OptimumTemp", c(1,pyramidLength) )
    #glsRangeAnalysisWithFilter( traits, tree, "GC.40.60",  "OptimumTemp", c(1,pyramidLength) )
    #glsRangeAnalysisWithFilter( traits, tree, "GC.60.100", "OptimumTemp", c(1,pyramidLength) )



    for( gr in taxGroups )
    {
        print(gr)
        results <- glsRangeAnalysisWithFilter( traits, tree, gr, "LogGrowthTime", c(1,pyramidLength), profileId=1, extras="")

        results <- glsRangeAnalysisWithFilter( traits, tree, gr, "LogGrowthTime", c(2,32), profileId=2, extras="")

        performGLSregressionWithFilter( traits, tree, gr,   "LogGrowthTime",   "Profile_1.21", plotRegression=TRUE )

        # Store peaks (for each range) for this group
        if( any( class(results) == "data.frame" ) && nrow(results) )
        {
            results <- results[results$Var1==results$Var2,]    # Only include single-window results

            # Iterate over each range to find the relevant peak
            for( i in 1:length(groupsTableOutputFile.limitRangeFromNt) )
            {
                matching <- (results$Var1 >= groupsTableOutputFile.limitRangeFromNt[i]) &
                            (results$Var2 <= groupsTableOutputFile.limitRangeToNt[i]  )  # Only include ranges within the configured limits

                maxResult <- results[matching,][which.max( abs(results[matching, "Buse.R2"]) ),]  # Choose the result with the highest R^2 

                regressionResultsByTaxGroup <- rbind( regressionResultsByTaxGroup, data.frame( ExplanatoryVar=c("GenomeSizeMb"),  Range=c(i-1), MaxRangeStart=c(maxResult$Var1), MaxRangeEnd=c(maxResult$Var2), TaxGroup=c(taxGroupToTaxId[gr]), TaxGroupName=c(gr), EffectSize=c(maxResult$Buse.R2), Pvalue=c(maxResult$Pvalue), NumSpecies=c(maxResult$NumSpecies) ) )
            }
        }
    }

    
    performGLSregressionWithFilter( traits, tree, "Member_Flavobacteriales_200644", "GenomicGC", "Profile_1.18" )


    performGLSregressionWithFilter( traits, tree, "Member_FCB_group_1783270",                  "Habitat",   "Profile_1.3" )
    performGLSregressionWithFilter( traits, tree, "Member_FCB_group_1783270",                  "Habitat",   "Profile_1.8" )
    performGLSregressionWithFilter( traits, tree, "Member_FCB_group_1783270",                  "Habitat",   "Profile_1.15" )

    performGLSregressionWithFilter( traits, tree, "Member_TACK_group_1783275",                 "Habitat",   "Profile.8" )
    performGLSregressionWithFilter( traits, tree, "Member_TACK_group_1783275",                 "Habitat",   "Profile.11" )
    #performGLSregressionWithFilter( traits, tree, "Member_TACK_group_1783275",                 "Habitat",   "Profile.27" )

    performGLSregressionWithFilter( traits, tree, "Member_Proteobacteria_1224",                "Habitat",   "Profile.13" )
    #performGLSregressionWithFilter( traits, tree, "Member_Proteobacteria_1224",                "Habitat",   "Profile.25" )

    performGLSregressionWithFilter( traits, tree, "Member_Firmicutes_1239",                    "Habitat",   "Profile.3" )
    performGLSregressionWithFilter( traits, tree, "Member_Firmicutes_1239",                    "Habitat",   "Profile.6" )


    performGLSregressionWithFilter( traits, tree, "Member_Terrabacteria_group_1783272",        "OxygenReq", "Profile.8" )
    performGLSregressionWithFilter( traits, tree, "Member_Terrabacteria_group_1783272",        "OxygenReq", "Profile.10" )
    #performGLSregressionWithFilter( traits, tree, "Member_Terrabacteria_group_1783272",        "OxygenReq", "Profile.19" )

    performGLSregressionWithFilter( traits, tree, "Member_Actinobacteria_201174",              "OxygenReq", "Profile.14" )
    #performGLSregressionWithFilter( traits, tree, "Member_Actinobacteria_201174",              "OxygenReq", "Profile.20" )
    #performGLSregressionWithFilter( traits, tree, "Member_Actinobacteria_201174",              "OxygenReq", "Profile.30" )

    performGLSregressionWithFilter( traits, tree, "Member_FCB_group_1783270",                  "OxygenReq", "Profile.8" )
    performGLSregressionWithFilter( traits, tree, "Member_FCB_group_1783270",                  "OxygenReq", "Profile.15" )

    #performGLSregressionWithFilter( traits, tree, "Member_Proteobacteria_1224",                "OxygenReq", "Profile.17" )
    #performGLSregressionWithFilter( traits, tree, "Member_Proteobacteria_1224",                "OxygenReq", "Profile.28" )

    performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2", "GenomicGC", "Profile.15" )
    #performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2", "GenomicGC", "Profile.18" )
    #performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2", "GenomicGC", "Profile.21" )


    performGLSregressionWithFilter( traits, tree, "Member_all_1", "OptimumTemp", "Profile.1" )
    performGLSregressionWithFilter( traits, tree, "HighGC",       "OptimumTemp", "Profile.1" )
    performGLSregressionWithFilter( traits, tree, "LowGC",        "OptimumTemp", "Profile.1" )

    performGLSregressionWithFilter( traits, tree, "Member_all_1", "OptimumTemp", "Profile.5" )
    performGLSregressionWithFilter( traits, tree, "HighGC",       "OptimumTemp", "Profile.5" )
    performGLSregressionWithFilter( traits, tree, "LowGC",        "OptimumTemp", "Profile.5" )

    performGLSregressionWithFilter( traits, tree, "Member_all_1", "OptimumTemp", "Profile.10" )
    performGLSregressionWithFilter( traits, tree, "HighGC",       "OptimumTemp", "Profile.10" )
    performGLSregressionWithFilter( traits, tree, "LowGC",        "OptimumTemp", "Profile.10" )

    performGLSregressionWithFilter( traits, tree, "Member_all_1", "OptimumTemp", "Profile.15" )
    performGLSregressionWithFilter( traits, tree, "HighGC",       "OptimumTemp", "Profile.15" )
    performGLSregressionWithFilter( traits, tree, "LowGC",        "OptimumTemp", "Profile.15" )

    #performGLSregressionWithFilter( traits, tree, "Member_all_1", "OptimumTemp", "Profile.20" )
    #performGLSregressionWithFilter( traits, tree, "HighGC",       "OptimumTemp", "Profile.20" )
    #performGLSregressionWithFilter( traits, tree, "LowGC",        "OptimumTemp", "Profile.20" )


    performGLSregressionWithFilter( traits, tree, "GC.0.40",   "OptimumTemp", "Profile.1" )
    performGLSregressionWithFilter( traits, tree, "GC.40.60" , "OptimumTemp", "Profile.1" )
    performGLSregressionWithFilter( traits, tree, "GC.60.100", "OptimumTemp", "Profile.1" )

    performGLSregressionWithFilter( traits, tree, "GC.0.40",   "OptimumTemp", "Profile.10" )
    performGLSregressionWithFilter( traits, tree, "GC.40.60" , "OptimumTemp", "Profile.10" )
    performGLSregressionWithFilter( traits, tree, "GC.60.100", "OptimumTemp", "Profile.10" )

    performGLSregressionWithFilter( traits, tree, "GC.0.40",   "OptimumTemp", "Profile.15" )
    performGLSregressionWithFilter( traits, tree, "GC.40.60" , "OptimumTemp", "Profile.15" )
    performGLSregressionWithFilter( traits, tree, "GC.60.100", "OptimumTemp", "Profile.15" )

    #performGLSregressionWithFilter( traits, tree, "GC.0.40",   "OptimumTemp", "Profile.20" )
    #performGLSregressionWithFilter( traits, tree, "GC.40.60" , "OptimumTemp", "Profile.20" )
    #performGLSregressionWithFilter( traits, tree, "GC.60.100", "OptimumTemp", "Profile.20" )

    performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",  "TempHighLow75", "Profile.1"   )
    performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",  "TempHighLow75", "Profile.5"   )
    performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",  "TempHighLow75", "Profile.11"  )

    performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",  "TempHighLow75", "Profile.15"  )
    #performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",  "TempHighLow75", "Profile.19"  )
    #performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",  "TempHighLow75", "Profile.24"  )

    performGLSregressionWithFilter( traits, tree, "Member_Archaea_2157", "TempHighLow75", "Profile.1"  )
    performGLSregressionWithFilter( traits, tree, "Member_Archaea_2157", "TempHighLow75", "Profile.5"  )
    performGLSregressionWithFilter( traits, tree, "Member_Archaea_2157", "TempHighLow75", "Profile.11" )

    performGLSregressionWithFilter( traits, tree, "Member_Archaea_2157", "TempHighLow75", "Profile.15" )
    performGLSregressionWithFilter( traits, tree, "Member_Archaea_2157", "TempHighLow75", "Profile.19" )
    performGLSregressionWithFilter( traits, tree, "Member_Archaea_2157", "TempHighLow75", "Profile.24" )


    performGLSregressionWithFilter( traits, tree, "Member_Opisthokonta_33154", "GenomicGC", "Profile.1" )
    performGLSregressionWithFilter( traits, tree, "Member_Opisthokonta_33154", "GenomicGC", "Profile.8" )
    performGLSregressionWithFilter( traits, tree, "Member_Opisthokonta_33154", "GenomicGC", "Profile.14" )
    #performGLSregressionWithFilter( traits, tree, "Member_Opisthokonta_33154", "GenomicGC", "Profile.20" )
    #performGLSregressionWithFilter( traits, tree, "Member_Opisthokonta_33154", "GenomicGC", "Profile.26" )

    #---------------------------------------------------------------------------------------------------------------------------
    # Save the summarized results, by group, to csv file (plot using: python2 plot_tree_effects_anaylsis_results_on_tree.py)
    write.csv(regressionResultsByTaxGroup, file=groupsTableOutputFile )
}


report_taxonRobustnessForTrait <- function( trait, rangeSpec, profileId=1 )
{
    regressionResultsByTaxGroup <- data.frame( ExplanatoryVar=character(), Range=integer(), MaxRangeStart=integer(), MaxRangeEnd=integer(), TaxGroup=integer(), TaxGroupName=character(), EffectSize=double(), Pvalue=double(), NumSpecies=integer() )

    for( gr in taxGroups )
    {
        print(gr)
        results <- glsRangeAnalysisWithFilter( traits, tree, gr, trait, rangeSpec, profileId=profileId )

        # Store peaks (for each range) for this group
        if( any( class(results) == "data.frame" ) && nrow(results) )
        {
            results <- results[results$Var1==results$Var2,]    # Only include single-window results
            results$Pvalue.adj <- p.adjust( results$Pvalue, method="fdr" )

            print(results)

            # Iterate over each range to find the relevant peak
            for( i in 1:length(groupsTableOutputFile.limitRangeFromNt) )
            {
                matching <- (results$Var1 >= groupsTableOutputFile.limitRangeFromNt[i]) &
                            (results$Var2 <= groupsTableOutputFile.limitRangeToNt[i]  )  # Only include ranges within the configured limits

                maxResult <- results[matching,][which.max( abs(results[matching, "Buse.R2"]) ),]  # Choose the result with the highest R^2 
                medianR2  <- median( results[matching, "Buse.R2"] )

                regressionResultsByTaxGroup <- rbind( regressionResultsByTaxGroup, data.frame( ExplanatoryVar=c(trait),  Range=c(i-1), MaxRangeStart=c(maxResult$Var1), MaxRangeEnd=c(maxResult$Var2), TaxGroup=c(taxGroupToTaxId[gr]), TaxGroupName=c(gr), EffectSize=c(medianR2), Pvalue=c(maxResult$Pvalue.adj), NumSpecies=c(maxResult$NumSpecies) ) )
            }
        }
    }

    #---------------------------------------------------------------------------------------------------------------------------
    # Save the summarized results, by group, to csv file (plot using: python2 plot_tree_effects_anaylsis_results_on_tree.py)
    write.csv(regressionResultsByTaxGroup, file=sprintf(groupsTableOutputFile, trait) )
    
    regressionResultsByTaxGroup
}


collectTraitInfluenceResultsForMidCDS <- function( allTraits, tree, groups, testTraits )
{
    results <- data.frame()
    for( gr in groups )
    {
        for( testTrait in testTraits )
        {
            r.gls.mid.start <- performGLSregression_profileRangeMean_withFilter( allTraits, tree, gr, testTrait, as.integer(c(11,31)), profileId=1, plotRegression=FALSE )
            if( !isTRUE(is.na(r.gls.mid.start)) )  # R can be weird... (is there a better way to do this?)
            {
                r.gls.mid.start$method <- "gls"
                r.gls.mid.start$trait <- testTrait
                r.gls.mid.start$group <- gr
                results <- rbind( results, r.gls.mid.start )
            }

            r.ols.mid.start <- performOLSregression_profileRangeMean_withFilter( allTraits,       gr, testTrait, as.integer(c(11,31)), profileId=1, plotRegression=FALSE )
            if( !isTRUE(is.na(r.ols.mid.start)) )  # R can be weird... (is there a better way to do this?)
            {
                r.ols.mid.start$method <- "ols"
                r.ols.mid.start$trait <- testTrait
                r.ols.mid.start$group <- gr
                results <- rbind( results, r.ols.mid.start )
            }

            #print( r.gls.mid.start )
            #print( r.ols.mid.start )

        }
    }

    return( results )
}


plotInflueceComparison <- function( results, title="" )
{
    if( nrow(results)==0 ) { return( ggplot() ) }

    results$significant1 <- (results$pvalue < 0.05 & results$pvalue >= 0.001)
    results$significant2 <- results$pvalue < 0.001
    print(results)

    s <- 0.4

    p <- ggplot( results, aes(y=trait, x=R2) ) +
        geom_vline( xintercept=0 )
    if( nrow(results[results$method=="gls",]) > 0 )
    {
        p <- p + geom_segment( data=results[results$method=="gls",], aes(y=trait, yend=trait, x=0, xend=R2, color=trait), size=2*s, alpha=1.0 )
    }
    print(s)
    print(4*s)

    p <- p +
        geom_point( aes(color=trait, shape=method), size=4*s ) +
        geom_text( aes(alpha=factor(significant1) ), label="*",  color="black", size=round(fontScale*s*0.8), nudge_y=0.05  ) +
        geom_text( aes(alpha=factor(significant2) ), label="**", color="black", size=round(fontScale*s), nudge_y=0.05  ) +
        scale_alpha_manual( values=c( 0.0, 1.0) ) +
        scale_x_continuous( limits=c(-0.8, 0.8) ) +
        scale_colour_manual( values=c("GenomicGC"="Blue", "GenomicENc.prime"="Purple", "LogGrowthTime"="Red", "LogGenomeSize"="Orange", "OptimumTemp"="Yellow") ) +
        labs(title=title, y="", x="") +
        guides( shape=FALSE, color=FALSE, alpha=FALSE ) +
        theme( plot.background = element_blank(),   # Hide unnecessary theme elements (background panels, etc.)
#              panel.grid.major.y = element_line(color="grey", size=0.50),
              panel.grid.major.y = element_blank(),
              panel.grid.major.x = element_line(color="grey", size=0.90*s),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.text.y = element_text( size=0 ),
              axis.text.x = element_text( size=round(fontScale*s*1.2) ),
              aspect.ratio = 0.12
              )
#        margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")


    print(p)
    return(p)
}


figure_DLFETraitsInfluencesComparisonForMidCDS <- function()
{
    grs <- c("Member_all_1", "Member_Bacteria_2", "Member_Proteobacteria_1224", "Member_Terrabacteria_group_1783272", "Member_FCB_group_1783270", "Member_Archaea_2157", "Member_Eukaryota_2759", "Member_Fungi_4751")
    trs <- c("GenomicGC", "GenomicENc.prime", "LogGrowthTime", "LogGenomeSize", "OptimumTemp")

    pdata <- collectTraitInfluenceResultsForMidCDS( traits, tree, grs, trs)

    pn <- list()
    i <- 1
    for( tr in trs )
    {
        for( gr in grs )
        {
            p <- plotInflueceComparison( pdata[ (pdata$gr==gr)&(pdata$trait==tr), ], title="" )
            pn[[i]] <- ggplotGrob( p ) # return list of grobs
            i <- i+1
        }
    }

    #print(pn)
    #print(length(pn))
    marrangeGrob( pn, nrow=length(grs), ncol=length(trs) ) #, heights=rep( unit(1/length(pn), "npc"), length(pn) ) )
}


figure_GC_vs_dLFE_in_Eukaryotes_GLS_MIC <- function()
{
    performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_Eukaryota_2759", "GenomicGC",  as.integer(c(11,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced (100-300nt)" )

    performOLSregression_profileRangeMean_withFilter( traits,       "Member_Eukaryota_2759", "GenomicGC",  as.integer(c(11,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced (100-300nt)", colorTrait="Member_Fungi_4751" )

    performOLSregression_profileRangeMean_withFilter( traits,       "Member_Eukaryota_2759", "GenomicGC",  as.integer(c(11,21)), profileId=1, plotRegression=TRUE, caption="Start-referenced (100-200nt)", colorTrait="Member_Fungi_4751" )

    performOLSregression_profileRangeMean_withFilter( traits,       "Member_Eukaryota_2759", "GenomicGC",  as.integer(c(21,31)), profileId=1, plotRegression=TRUE, caption="Start-referenced (200-300nt)", colorTrait="Member_Fungi_4751" )
    
    performGLSregression_profileRangeMean_withFilter( traits, tree, "Member_Eukaryota_2759", "GenomicGC",  as.integer(c(1,21)), profileId=2, plotRegression=TRUE, caption="End-referenced (-300 - -100nt)" )

    performOLSregression_profileRangeMean_withFilter( traits,       "Member_Eukaryota_2759", "GenomicGC",  as.integer(c(1,21)), profileId=2, plotRegression=TRUE, caption="End-referenced (-300 - -100nt)",    colorTrait="Member_Fungi_4751" )
    
    print("--------------------------------------------------------------------")
    print(sprintf("# iters for MIC p-val: %d", MIC.pval.num.iterations))
    print("--------------------------------------------------------------------")
}

############################################################

report_taxonRobustnessForTrait( "GenomicGC",          rangeSpec=c(1,31), profileId=1 )
report_taxonRobustnessForTrait( "GenomicENc.prime",   rangeSpec=c(1,31), profileId=1 )
report_taxonRobustnessForTrait( "OptimumTemp",        rangeSpec=c(1,31), profileId=1 )


#figure_PartialDeterminationAnalysis_GC_and_ENc.prime()
#figure_PartialDeterminationAnalysis(profileId=1, pyramidSpec=c(1,31) )
#igure_PartialDeterminationAnalysis(profileId=2, pyramidSpec=c(2,32) )
#figure_PartialDeterminationAnalysis2()
#figure_PartialDeterminationAnalysis_NormalizedProfiles()
#figure_PositiveStretchLengths()
#figure_DLFEInteractingTraits_RegressionRangeAnalysisByTaxGroup()
#figure_DLFETraitsInfluencesComparisonForMidCDS()
#figure_GC_vs_dLFE_in_Eukaryotes_GLS_MIC()
#figure_CorrelationBetweenModelRegions()
#figure_contrastingKDEsForHighLowGC()
#figure_contrastingBoxplotsHighLowTemp()
#figure_ContrasingProfileBoxplotsForHighLowGC()
#figure_CorrelationBetweenRanges()
#writeWeakDLFEBinaryModelGridSearchResults()
#figure_GLS_GC_vs_endosymbionts()
#report_testRegressionForWeakModelComponents()
#report_testCompoundClassification()
#TestCorrelationBetweenScalarTraits()


############################################################

print("----------------------------------------------------------------")
print(paste0("Note: profileMode == ", profileMode))
print("----------------------------------------------------------------")

warnings()
dev.off()
