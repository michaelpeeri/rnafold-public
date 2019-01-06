library("nlme")
library("ggplot2")
library("reshape")
library("rhdf5")
library("rredis")
library("phylobase")
library("ape")
library("adephylo")


#Ytrait="Profile.15"
#Xtrait="GenomicGC"

#Ytrait="Profile.15"
#Xtrait="OptimumTemp"

#Ytrait="PairedFraction"
#Xtrait="Profile.15"

#Ytrait="Profile.15"
#Xtrait="Salinity"

#Ytrait="Profile.15"
#Xtrait="Habitat"

#Ytrait="Profile.15"
#Xtrait="OxygenReq"

Ytrait="Profile.15"
Xtrait="TemperatureRange"


pdf("tree_phenotypes_anova.out.pdf")

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
    #combined <- data.frame( GenomicGC=double(), profile.1=double(), profile.2=double(), profile.3=double(), profile.4=double(), profile.5=double(), profile.6=double(), profile.7=double(), profile.8=double(), profile.9=double(), profile.10=double(), profile.11=double(), profile.12=double(), profile.13=double(), profile.14=double(), profile.15=double(), row.names=factor() )
    #combined <- data.frame( GenomicGC=double(), profile.1=double(), profile.2=double(), profile.3=double(), profile.4=double(), profile.5=double(), profile.6=double(), profile.7=double(), profile.8=double(), profile.9=double(), profile.10=double(), profile.11=double(), profile.12=double(), profile.13=double(), profile.14=double(), profile.15=double(), row.names=integer() )

    combined <- data.frame( GenomicGC=double(), OptimumTemp=double(), TemperatureRange=ordered(c(), temperatureLevels), PairedFraction=double(), Profile=matrix( rep(0.0, profileLen), nrow=0, ncol=profileLen), Salinity=ordered(c(), salinityLevels), Habitat=factor(c(), habitatLevels), OxygenReq=ordered(c(), oxygenLevels), row.names=integer() )
    #
    #print("Prototype:")
    #print(class(combined$Salinity))
    #print(dimnames(combined))
    #print(combined)
    
    for (taxId in taxIds)
    {
        #print("----------------------------------------------")
        profile <- readDeltaLFEProfile( taxIds, getH5Filename( taxId ) )

        gcContent <- getGenomicGCContent( taxId )

        optimumTemperature <- getTemperature( taxId )

        temperatureRange <- getTemperatureCat( taxId )

        pairedFraction <- getPairedFraction( taxId )

        salinity <- getSalinity( taxId )
        #print("Value:")
        #print(salinity)
        #print(class(salinity))
        #print(dim(salinity))

        habitat <- getHabitat( taxId )

        oxygenReq <- getOxygenReq( taxId )
                
        newrow <- data.frame( GenomicGC=c(gcContent), OptimumTemp=c(optimumTemperature), TemperatureRange=temperatureRange, PairedFraction=c(pairedFraction), Profile=t(profile), Salinity=salinity, Habitat=habitat, OxygenReq=oxygenReq, row.names=c(taxId) )
        #print("New row:")
        #print(newrow)
        
        combined <- rbind(combined, newrow)
        #print("Combined:")
        #print(class(combined$Salinity))
    }
    #print(class(combined$Salinity))

    return(combined)
}



getAllTaxIds <- function()
{
    f <- function(s) { as.integer( strsplit(s, "_")[[1]][[4]] ) }
    glob1 <- sprintf("gcdata_v2_taxid_*_profile_%d_%d_%s_%d.h5", profileStop, profileStep, profileReference, profileStart )
    lapply( list.files(pattern=glob2rx(glob1)), f )
}


allTaxIds = getAllTaxIds()

traits <- readAllProfiles( allTaxIds )
#print(class(traits$Salinity))

stopifnot(nrow(traits) > 300 )
stopifnot(ncol(traits) == profileLen+7 )


#x["511145",]
#x[as.character(511145),]


traits["580340",]
#x[c(580340),]
#print(x)

print("###################")
print(sum(is.na(traits$GenomicGC)))
print(sum(is.na(traits$Profile.15)))


speciesWithMissingData <- row.names(traits[(is.na(traits[Xtrait]) | is.na(traits$GenomicGC)),])




#----------------------------------


# Source: http://www2.math.su.se/PATHd8/
# Created using: ~/src/PATHd8/PATHd8 ./data/nmicrobiol201648-s6.txt ./data/nmicrobiol201648-s6.txt.nw.PATHd8.nw
#tree <- read.tree("nmicro_s6_pruned_with_taxids.nw")
tree <- read.tree("test_tree.nw")

tree <- drop.tip( tree, c(c("470", "1280", "4932", "2850"), speciesWithMissingData) )

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

#print(traits)


traitsTree <- phylo4d(tree, tip.data=traits, rownamesAsLabels=TRUE)


hasTipData(traitsTree)

#table.phylo4d(traitsTree, box=FALSE)

#-------------------------------------------------------------------------------------

vcv1 <- vcv( phy=tree, model="Brownian")
dimnames(vcv1) <- list(1:ncol(vcv1))
bmcorrmtx <- data.frame( melt( vcv1 ) )
#print(range(bmcorrmtx$X1))
#print(range(bmcorrmtx$X2))
#print(range(bmcorrmtx$value))
#dimnames(bmcorrmtx)
p <- ggplot( data=bmcorrmtx, aes(x=X1, y=X2, fill=value) ) + geom_raster() + scale_y_reverse(); p



#-------------------------------------------------------------------------------------
#x <- 1:5; coef(lm(c(1:3, 7, 6) ~ x))


YvsX = as.formula( paste(Ytrait, " ~ ", Xtrait) )
#print(YvsX)

m1 <- lm(  YvsX, traits); summary(m1)
co <- coef(m1)
print(co)

print("//////////////")
print( nrow(traits) )
print( nTips(traitsTree) )

clr2 <- "#55CCB0"

gls1 <- gls( YvsX, traits, correlation=bmcorr, na.action=na.omit); summary(gls1)
#gls1 <- pGLS( YvsX, traits, bmcorr, na.action=na.omit); summary(gls1)
co2 <- coef(gls1)

plot(gls1) # Diagnostic plots



p <- ggplot(traits, aes(get(Xtrait), get(Ytrait))) +
  labs(y=Ytrait, x=Xtrait) +
  geom_point() +
  geom_hline( yintercept = 0 ) +
  geom_abline( aes( slope=co[Xtrait],  intercept=co["(Intercept)"]), colour="red"  ) +
  geom_abline( aes( slope=co2[Xtrait], intercept=co2["(Intercept)"]), colour=clr2  ) +
  annotate( "text", x=Inf,  y=Inf, label=sprintf("lm\nR^2 = %4g\nP-val = %4g\nn = %d", summary(m1)$r.squared,      summary(m1)$coefficients[2,4], nrow(traits) ), hjust=1, vjust=1, colour="red" ) +
  annotate( "text", x=-Inf, y=Inf, label=sprintf("gls\n\n\nP-val = %4g\nn = %d", summary(gls1)$tTable[2,4],     nrow(traits) ), hjust=0, vjust=1, colour=clr2 ); p


#-----------

P <- 5
bmTraits <- replicate(P, rTraitCont(tree))
bmTraits <- data.frame( bm1=bmTraits[,1], bm2=bmTraits[,2], bm3=bmTraits[,3], bm4=bmTraits[,4], bm5=bmTraits[,5] )
bmTraits <- cbind( traits, bmTraits )
#bmTree <- phylo4d(tree, correlatedPhenotypes)
#print(bmTraits)


m1 <- lm(  bm1 ~ GenomicGC, bmTraits); summary(m1)
co <- coef(m1)
print(co)
#print(co)

gls1 <- gls( bm1 ~ GenomicGC, bmTraits, correlation=bmcorr, na.action=na.omit); #summary(gls1)
co2 <- coef(gls1)

p <- ggplot(bmTraits, aes(GenomicGC, bm1)) +
  geom_point() +
  geom_hline( yintercept = 0 ) +
  geom_abline( aes( slope=co["GenomicGC"],  intercept=co["(Intercept)"]), colour="red"  ) +
  geom_abline( aes( slope=co2["GenomicGC"], intercept=co2["(Intercept)"]), colour=clr2  ) +
  annotate( "text", x=Inf,  y=Inf, label=sprintf("lm\nR^2 = %4g\nP-val = %4g\nn = %d", summary(m1)$r.squared,      summary(m1)$coefficients[2,4], nrow(traits) ), hjust=1, vjust=1, colour="red" ) +
  annotate( "text", x=-Inf, y=Inf, label=sprintf("gls\nR^2 = %4g\nP-val = %4g\nn = %d", (summary(gls1)$corBeta[1,2])^2, summary(gls1)$tTable[2,4],     nrow(traits) ), hjust=0, vjust=1, colour=clr2 ); p

#-

m1 <- lm(  bm2 ~ GenomicGC, bmTraits); summary(m1)
co <- coef(m1)
#print(co)

gls1 <- gls( bm2 ~ GenomicGC, bmTraits, correlation=bmcorr, na.action=na.omit); #summary(gls1)
co2 <- coef(gls1)

p <- ggplot(bmTraits, aes(GenomicGC, bm2)) +
  geom_point() +
  geom_hline( yintercept = 0 ) +
  geom_abline( aes( slope=co["GenomicGC"],  intercept=co["(Intercept)"]), colour="red"  ) +
  geom_abline( aes( slope=co2["GenomicGC"], intercept=co2["(Intercept)"]), colour=clr2  ) +
  annotate( "text", x=Inf,  y=Inf, label=sprintf("lm\nR^2 = %4g\nP-val = %4g\nn = %d", summary(m1)$r.squared,      summary(m1)$coefficients[2,4], nrow(traits) ), hjust=1, vjust=1, colour="red" ) +
  annotate( "text", x=-Inf, y=Inf, label=sprintf("gls\nR^2 = %4g\nP-val = %4g\nn = %d", (summary(gls1)$corBeta[1,2])^2, summary(gls1)$tTable[2,4],     nrow(traits) ), hjust=0, vjust=1, colour=clr2 ); p

m1 <- lm(  bm3 ~ GenomicGC, bmTraits); summary(m1)
co <- coef(m1)
#print(co)

gls1 <- gls( bm3 ~ GenomicGC, bmTraits, correlation=bmcorr, na.action=na.omit); #summary(gls1)
co2 <- coef(gls1)

p <- ggplot(bmTraits, aes(GenomicGC, bm3)) +
  geom_point() +
  geom_hline( yintercept = 0 ) +
  geom_abline( aes( slope=co["GenomicGC"],  intercept=co["(Intercept)"]), colour="red"  ) +
  geom_abline( aes( slope=co2["GenomicGC"], intercept=co2["(Intercept)"]), colour=clr2  ) +
  annotate( "text", x=Inf,  y=Inf, label=sprintf("lm\nR^2 = %4g\nP-val = %4g\nn = %d", summary(m1)$r.squared,      summary(m1)$coefficients[2,4], nrow(traits) ), hjust=1, vjust=1, colour="red" ) +
  annotate( "text", x=-Inf, y=Inf, label=sprintf("gls\nR^2 = %4g\nP-val = %4g\nn = %d", (summary(gls1)$corBeta[1,2])^2, summary(gls1)$tTable[2,4],     nrow(traits) ), hjust=0, vjust=1, colour=clr2 ); p

m1 <- lm(  bm4 ~ GenomicGC, bmTraits); summary(m1)
co <- coef(m1)
#print(co)

gls1 <- gls( bm4 ~ GenomicGC, bmTraits, correlation=bmcorr, na.action=na.omit); #summary(gls1)
co2 <- coef(gls1)

p <- ggplot(bmTraits, aes(GenomicGC, bm4)) +
  geom_point() +
  geom_hline( yintercept = 0 ) +
  geom_abline( aes( slope=co["GenomicGC"],  intercept=co["(Intercept)"], colour="red" ) ) +
  geom_abline( aes( slope=co2["GenomicGC"], intercept=co2["(Intercept)"], colour=clr2 ) ) +
  annotate( "text", x=Inf,  y=Inf, label=sprintf("lm\nR^2 = %4g\nP-val = %4g\nn = %d", summary(m1)$r.squared,      summary(m1)$coefficients[2,4], nrow(traits) ), hjust=1, vjust=1, colour="red" ) +
  annotate( "text", x=-Inf, y=Inf, label=sprintf("gls\nR^2 = %4g\nP-val = %4g\nn = %d", (summary(gls1)$corBeta[1,2])^2, summary(gls1)$tTable[2,4],     nrow(traits) ), hjust=0, vjust=1, colour=clr2 ); p

m1 <- lm(  bm5 ~ GenomicGC, bmTraits); summary(m1)
co <- coef(m1)
#print(co)

gls1 <- gls( bm5 ~ GenomicGC, bmTraits, correlation=bmcorr, na.action=na.omit); #summary(gls1)
co2 <- coef(gls1)

p <- ggplot(bmTraits, aes(GenomicGC, bm5)) +
  geom_point() +
  geom_hline( yintercept = 0 ) +
  geom_abline( aes( slope=co["GenomicGC"],  intercept=co["(Intercept)"]), colour="red"  ) +
  geom_abline( aes( slope=co2["GenomicGC"], intercept=co2["(Intercept)"]), colour=clr2  ) +
  annotate( "text", x=Inf,  y=Inf, label=sprintf("lm\nR^2 = %4g\nP-val = %4g\nn = %d", summary(m1)$r.squared,      summary(m1)$coefficients[2,4], nrow(traits) ), hjust=1, vjust=1, colour="red" ) +
  annotate( "text", x=-Inf, y=Inf, label=sprintf("gls\nR^2 = %4g\nP-val = %4g\nn = %d", (summary(gls1)$corBeta[1,2])^2, summary(gls1)$tTable[2,4],     nrow(traits) ), hjust=0, vjust=1, colour=clr2 ); p

#-----------

m1 <- lm(  Profile.15 ~ Profile.5, bmTraits); summary(m1)
co <- coef(m1)
#print(co)

gls1 <- gls( Profile.15 ~ Profile.5, bmTraits, correlation=bmcorr, na.action=na.omit); #summary(gls1)
co2 <- coef(gls1)

p <- ggplot(bmTraits, aes(Profile.5, Profile.15)) +
  geom_point() +
  geom_hline( yintercept = 0 ) +
  geom_abline( aes( slope=co["Profile.5"],  intercept=co["(Intercept)"]), colour="red"  ) +
  geom_abline( aes( slope=co2["Profile.5"], intercept=co2["(Intercept)"]), colour=clr2  ) +
  annotate( "text", x=Inf,  y=Inf, label=sprintf("lm\nR^2 = %4g\nP-val = %4g\nn = %d", summary(m1)$r.squared,      summary(m1)$coefficients[2,4], nrow(traits) ), hjust=1, vjust=1, colour="red" ) +
  annotate( "text", x=-Inf, y=Inf, label=sprintf("gls\nR^2 = %4g\nP-val = %4g\nn = %d", (summary(gls1)$corBeta[1,2])^2, summary(gls1)$tTable[2,4],     nrow(traits) ), hjust=0, vjust=1, colour=clr2 ); p




dev.off()