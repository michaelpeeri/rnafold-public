# Count the number of species with unusual profiles in each taxonomic group. Create a csv file compatible with tree_traits_effects_analysis_with_taxgroups.out.*.length.*.csv (tree_traits_effects_anaylsis_with_taxgroups.r) and plot_tree_effects_anaylsis_result_on_tree.py.
library("rhdf5")
library("rredis")
library("dplyr")

#---------------------------------------------------------------
# Configuration

redisConnect(host="power5", password="rnafold")

#profileStart <- 0
#profileStop <- 1000
#profileStep <- 10
#profileReference <- "begin"
#profileLen <- (profileStop-profileStart)/profileStep

profileStart <- 0
profileStop <- 150
profileStep <- 10
profileReference <- "end"
profileLen <- (profileStop-profileStart)/profileStep

taxidToKingdomFilename <- "TaxidToKingdom.csv"   # create file using: python2 create_taxid_kingdom_table.py


#profileMode <- "nativeLFE"
#profileMode <- "shuffledLFE"
profileMode <- "dLFE"


groupsTableOutputFile <- sprintf("find_trait_values_outliers.out.%s.csv", profileMode )


getH5Filename <- function( taxId )
{
    #                gcdata_v2_taxid_1454006_profile_150_10_end_0_t11.h5
    return( sprintf("gcdata_v2_taxid_%d_profile_%d_%d_%s_%d_t11.h5", taxId, profileStop, profileStep, profileReference, profileStart) )
}

# Read native and shuffled LFE profiles from the hdf files
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

    if(      profileMode == "dLFE" )
    {
        return( native - shuffled )
    }
    else if( profileMode == "nativeLFE" )
    {
        return( native )
    }
    else if( profileMode == "shuffledLFE" )
    {
        return( shuffled )
    }
    else
    {
        stopifnot(FALSE)
        return( NA )
    }
    #return( gc )
}

getAllTaxIds <- function()
{
    f <- function(s) { as.integer( strsplit(s, "_")[[1]][[4]] ) }
    glob1 <- sprintf("gcdata_v2_taxid_*_profile_%d_%d_%s_%d_t11.h5", profileStop, profileStep, profileReference, profileStart )
    lapply( list.files(pattern=glob2rx(glob1)), f )
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

    combined <- data.frame( GenomicGC=double(), Profile=matrix( rep(0.0, profileLen), nrow=0, ncol=profileLen), row.names=integer() )

    for (taxId in taxIds)
    {
        profile <- readDeltaLFEProfile( taxIds, getH5Filename( taxId ) )

        gcContent <- getGenomicGCContent( taxId )

        #optimumTemperature <- getTemperature( taxId )

        #temperatureRange <- getTemperatureCat( taxId )

        #pairedFraction <- getPairedFraction( taxId )
        
        #salinity <- getSalinity( taxId )
        
        #habitat <- getHabitat( taxId )

        #oxygenReq <- getOxygenReq( taxId )
                
        newrow <- data.frame( GenomicGC=c(gcContent), Profile=t(profile), row.names=c(taxId) )
        
        combined <- rbind(combined, newrow)
    }
    return( combined)
}


taxidToKingdom <- read.table(taxidToKingdomFilename, sep=",", header=TRUE )    # create file using: python2 create_taxid_kingdom_table.py
taxidToKingdom$Member_all_1 <- 1   # Add a special group which includes all species

allTaxIds <- getAllTaxIds() # get list of all taxids in the current dataset

# Read profiles for all species
traits <- readAllProfiles( allTaxIds )
traits$tax_id <- rownames(traits)  # create an explicit tax_id column (rather than an index), by which we'll be able to merge (since left_join can't merge using indices)


taxGroups <- colnames(taxidToKingdom)
taxGroups <- taxGroups[startsWith(taxGroups, "Member_")]

taxGroupToTaxId <- vapply(taxGroups, function (x) as.integer(tail(strsplit(x[[1]], '_')[[1]], n=1)), integer(1) )

taxidToKingdom$tax_id <- as.character(taxidToKingdom$tax_id)  # Convert the taxids column in the taxon memberships table to string, to match traits table

traits <- traits %>% left_join(taxidToKingdom, by="tax_id")  # Perform the merge
rownames(traits) <- traits$tax_id  # Restore the tax_ids index


findValueOutliersWithFilter <- function( traits, filterTrait, valuesTrait, criterion, pvalThreshold=0.01 )
{
    #print(filterTrait)
    #print(valuesTrait)
    #print(colnames(traits))
    stopifnot( any( filterTrait == colnames(traits)  ) )  # filterTrait not found
    stopifnot( any( valuesTrait == colnames(traits)  ) )  # studyTrait  not found


    # Filter by the filter trait
    traits <- traits[ (traits[,filterTrait]>0), ]

    numSpecies <- nrow(traits)

    
    #print(nrow(traits))

    if( criterion==">0" )
    {
        matchingRows <- sum(traits[,valuesTrait] >  0.05)
    }
    else
    {
        matchingRows <- sum(traits[,valuesTrait] < -0.05)
    }
    print(matchingRows)

    return( data.frame( matches=c(matchingRows), criterion=c(criterion), numSpecies=c(numSpecies) ) )
}


traits[,"test"] <- traits[,"Profile.1"] - traits[,"Profile.16"]

print(traits["469371",])
print("mean.1=")
print(mean(traits[,"Profile.1"]))
print("mean.15=")
print(mean(traits[,"Profile.15"]))
print("mean.16=")
print(mean(traits[,"Profile.16"]))

print(dimnames(traits))

allResults <- data.frame( ExplanatoryVar=character(), Range=integer(), MaxRangeStart=integer(), MaxRangeEnd=integer(), TaxGroup=integer(), TaxGroupName=character(), EffectSize=double(), Pvalue=double(), NumSpecies=integer() )

gr <- "Member_all_1";
print("Profile.1 < 0 (red)")
print( findValueOutliersWithFilter( traits, gr, "Profile.1", "<0" ) )
print("Profile.1 > 0 (blue)")
print( findValueOutliersWithFilter( traits, gr, "Profile.1", ">0" ) )

print("Profile.15 < 0 (red)")
print( findValueOutliersWithFilter( traits, gr, "Profile.15", "<0" ) )
print("Profile.15 > 0 (blue)")
print( findValueOutliersWithFilter( traits, gr, "Profile.15", ">0" ) )

print("Profile.16 < 0 (red)")
print( findValueOutliersWithFilter( traits, gr, "Profile.16", "<0" ) )
print("Profile.16 > 0 (blue)")
print( findValueOutliersWithFilter( traits, gr, "Profile.16", ">0" ) )


print("test < 0")
print( findValueOutliersWithFilter( traits, gr, "test", "<0" ) )
print("test > 0")
print( findValueOutliersWithFilter( traits, gr, "test", ">0" ) )

quit()


for( gr in taxGroups )
{
    print(gr)
    results <- findValueOutliersWithFilter( traits, gr, "Profile.1", "<0" )

    if( results$matches > 0 )
    {
        # plot_tree_effects_anaylsis_result_on_tree.py expects to have data for two ranges at least (0 and 1). Here we don't actually do the analysis by ranges.
        allResults <- rbind( allResults, data.frame( ExplanatoryVar=c("Profile.1.below.0"),   Range=c(0), MaxRangeStart=c(0), MaxRangeEnd=c(0), TaxGroup=c(taxGroupToTaxId[gr]), TaxGroupName=c(gr), EffectSize=c(results$matches), Pvalue=c(1.0), NumSpecies=c(results$numSpecies) ) )
        allResults <- rbind( allResults, data.frame( ExplanatoryVar=c("Profile.1.below.0"),   Range=c(1), MaxRangeStart=c(0), MaxRangeEnd=c(0), TaxGroup=c(taxGroupToTaxId[gr]), TaxGroupName=c(gr), EffectSize=c(results$matches), Pvalue=c(1.0), NumSpecies=c(results$numSpecies) ) )
    }
            
    results <- findValueOutliersWithFilter( traits, gr, "Profile.16", ">0" )

    if( results$matches > 0 )
    {
        allResults <- rbind( allResults, data.frame( ExplanatoryVar=c("Profile.16.above.0"),  Range=c(0), MaxRangeStart=c(0), MaxRangeEnd=c(0), TaxGroup=c(taxGroupToTaxId[gr]), TaxGroupName=c(gr), EffectSize=c(results$matches), Pvalue=c(1.0), NumSpecies=c(results$numSpecies) ) )
        allResults <- rbind( allResults, data.frame( ExplanatoryVar=c("Profile.16.above.0"),  Range=c(1), MaxRangeStart=c(0), MaxRangeEnd=c(0), TaxGroup=c(taxGroupToTaxId[gr]), TaxGroupName=c(gr), EffectSize=c(results$matches), Pvalue=c(1.0), NumSpecies=c(results$numSpecies) ) )
    }

    # Store peaks (for each range) for this group

}

#print(allResults)

write.csv( allResults, file=groupsTableOutputFile )


