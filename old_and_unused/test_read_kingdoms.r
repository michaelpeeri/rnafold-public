library("rhdf5")
library("rredis")
library("dplyr")

#---------------------------------------------------------------
# Configuration


redisConnect(host="power5", password="rnafold")

profileStart <- 0
profileStop <- 1000
profileStep <- 10
profileReference <- "begin"
profileLen <- (profileStop-profileStart)/profileStep



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

    combined <- data.frame( GenomicGC=double(), OptimumTemp=double(), TemperatureRange=ordered(c(), c(NA)), PairedFraction=double(), Profile=matrix( rep(0.0, 10), nrow=0, ncol=10), Salinity=ordered(c(), c(NA)), Habitat=factor(c(), c(NA)), OxygenReq=ordered(c(), c(NA)), ProteinCount=integer(), GenomeSizeMb=double(), GenomeDensity=double(), tax_id=integer(), row.names=integer() )

    for (taxId in taxIds)
    {
        #profile <- readDeltaLFEProfile( taxIds, getH5Filename( taxId ) )

        gcContent <- getGenomicGCContent( taxId )

        #optimumTemperature <- getTemperature( taxId )

        #temperatureRange <- getTemperatureCat( taxId )

        #pairedFraction <- getPairedFraction( taxId )
        
        #salinity <- getSalinity( taxId )
        
        #habitat <- getHabitat( taxId )

        #oxygenReq <- getOxygenReq( taxId )

        #proteinCount <- getProteinCount( taxId )

        #genomeSizeMb <- getGenomeSizeMb( taxId )

        #genomeDensity <- proteinCount / genomeSizeMb
                
        newrow <- data.frame( GenomicGC=c(gcContent), OptimumTemp=c(21), TemperatureRange=NA, PairedFraction=c(NA), Profile=t(1:10), Salinity=NA, Habitat=NA, OxygenReq=NA, ProteinCount=NA, GenomeSizeMb=NA, GenomeDensity=NA, tax_id=taxId, row.names=c(taxId) )
        
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


#---------------------------------------------------------------
taxidToKingdom <- read.table("TaxidToKingdom.csv", sep=",", header=TRUE)
print(colnames(taxidToKingdom))

#print(class(taxidToKingdom$kingdom))
#print(class(taxidToKingdom$tax_id))
#print(class(traits$tax_id))

#traits$kingdom <- taxidToKingdom

traits$tax_id <- rownames(traits)
print(colnames(traits))
print(traits$tax_id)

traits <- traits %>% left_join(taxidToKingdom, by="tax_id")
rownames(traits) <- traits$tax_id

print(traits$kingdom)

print(rownames(traits))
print(traits["511145",])

print(sum(traits$kingdom=="Bacteria"))
print(sum(traits$kingdom=="Archaea"))
print(sum(traits$kingdom=="Eukaryota"))


warnings()

