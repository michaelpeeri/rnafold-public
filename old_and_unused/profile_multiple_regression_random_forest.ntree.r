library("rhdf5")
library("rredis")
library("randomForest")
library("ggplot2")
library("scales") # for pretty_breaks

#---------------------------------------------------------------
# Configuration

pdf("profile_multiple_regression_random_forest.ntree.out.pdf")

redisConnect(host="power5", password="rnafold")

profileStart <- 0
profileStop <- 1000
profileStep <- 10
profileReference <- "begin"
profileLen <- (profileStop-profileStart)/profileStep



#---------------------------------------------------------------

# Initialize random seed
seed <- as.integer(difftime( as.POSIXlt(Sys.time()), as.POSIXlt("2010-01-01 12:00:00"), units="secs" ))
set.seed( seed )
print(sprintf("Using random seed %d", seed))



#---------------------------------------------------------------
taxidToKingdom <- read.table("TaxidToKingdom.csv", sep=",", header=TRUE, row.names="tax_id")

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

getGenomeSizeMb <- function(taxId)
{
    val <- redisGet(sprintf("species:taxid:%s:properties:genome-size-mb", taxId))
    
    if( !is.null(val)) {
        val <- as.double(val)
        stopifnot( val > 0.03 && val < 5000 )
        return( val )
    }
    else
    {
        return( NA )  # convert NULL to NA
    }
}


readAllProfiles <- function( taxIds )
{

    combined <- data.frame( GenomicGC=double(), OptimumTemp=double(), TemperatureRange=ordered(c(), temperatureLevels), PairedFraction=double(), Profile=matrix( rep(0.0, profileLen), nrow=0, ncol=profileLen), Salinity=ordered(c(), salinityLevels), Habitat=factor(c(), habitatLevels), OxygenReq=ordered(c(), oxygenLevels), ProteinCount=integer(), GenomeSizeMb=double(), GenomeDensity=double(), row.names=integer() )

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

        proteinCount <- getProteinCount( taxId )

        genomeSizeMb <- getGenomeSizeMb( taxId )

        genomeDensity <- proteinCount / genomeSizeMb
                
        newrow <- data.frame( GenomicGC=c(gcContent), OptimumTemp=c(optimumTemperature), TemperatureRange=temperatureRange, PairedFraction=c(pairedFraction), Profile=t(profile), Salinity=salinity, Habitat=habitat, OxygenReq=oxygenReq, ProteinCount=proteinCount, GenomeSizeMb=genomeSizeMb, GenomeDensity=genomeDensity, row.names=c(taxId) )
        
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
stopifnot(ncol(traits) == profileLen+10 )


#x["511145",]
#x[as.character(511145),]


#traits["580340",]
#x[c(580340),]
#print(x)

print("###################")
print(sum(is.na(traits$GenomicGC)))
print(sum(is.na(traits$Profile.15)))


print(nrow(traits))


#-----------------------------------------


getRandomForestRegression <- function( formula, data, ... )
{
    rf <- randomForest( as.formula(formula), data=data, ... )
    return( rf )
}

performRegressionWithCrossVal <- function( formula, data, N=10, ... )
{
    Nrows <- nrow(data)
    data <- cbind( data, CrossValGroup = sample( rep(1:N, len=Nrows), Nrows, replace=FALSE ) )  # Divide the data into N approximately equal groups

    s <- 0

    effectVar <- all.vars(as.formula(formula))[1]  # extract the name of the first var (i.e., the one on the left side) from the forumla

    groupStats <- data.frame( Group=1:N, Pearson=rep(0,N), MSE=rep(0,N), MAE=rep(0,N) )
    
    for( i in 1:N )
    {
        trainData <- data[data$CrossValGroup != i, ]
        testData  <- data[data$CrossValGroup == i, ]
        stopifnot( nrow(testData) > 3 )     # min. num or rows in test-set
        stopifnot( nrow(trainData) + nrow(testData) == nrow(data) )    # all data must be used

        rf <- getRandomForestRegression( formula, trainData, ... )

        out <- predict( rf, testData )

        groupStats[i, "Pearson"] = cor (     testData[, effectVar],   out,      method="pearson" )
        groupStats[i, "MSE"]     = mean(    (testData[, effectVar] -  out)**2 )
        groupStats[i, "MAE"]     = mean( abs(testData[, effectVar] -  out) )

        print(sprintf("(%d) %d", i, nrow(testData) ))
        s <- s + nrow(testData)
    }

    print("-------------")
    stopifnot(s == Nrows)
    print(s)

    return(groupStats)
}

performRegressionWithRepeatedCrossVal <- function( formula, data, N=10, M=10, ... )
{
    allStats <- data.frame( Iter=c(), Group=c(), Pearson=c(), MSE=c(), MAE=c() )
    
    for( i in 1:M )
    {
        groupStats <- performRegressionWithCrossVal( formula, data, N=N, ...)
        groupStats$Iter <- rep(i, N)

        allStats <- rbind( allStats, groupStats )
    }
    return( allStats )
}

testRandomTreeParams <- function( formula, data, ... )
{
    out <- data.frame( mtry=integer(), ntree=double(), pearson=double() )
    
    #for( mtry in seq(2, 9, 1) )
    for( mtry in c(7) )
    {
        #for( nodesize in seq(0.5, 21.5, 3) )
        for( ntree in seq(50, 2500, 50) )
        {
            print(sprintf("mtry=%d, ntree=%d", mtry, ntree))
            stats <- performRegressionWithRepeatedCrossVal( formula, data, mtry=mtry, nodesize=10.5, ntree=ntree, ... )
            pearsonByIter <- aggregate( stats$Pearson, by=list(stats$Iter), mean )
            colnames(pearsonByIter) <- c("Iter", "MeanPearson")

            meanPearson <- mean( pearsonByIter$MeanPearson )

            #out <- rbind( out, data.frame( mtry=c(mtry), nodesize=c(nodesize), pearson=c(mtry+nodesize) ) )  # dummy values
            out <- rbind( out, data.frame( mtry=c(mtry), ntree=c(ntree), pearson=c(meanPearson) ) )
        }
    }
    
    p <- ggplot( data=out, aes( x=mtry, y=ntree, fill=pearson ) ) +
        geom_raster() +
        scale_x_continuous(breaks = pretty_breaks())
    print(p)

    print("Optimum:")
    print( out[ which.max( out$pearson ), ])
    
    return( out )
}


#----------------------------------------
N <- 5
M <- 200
modelFormula <- "Profile.26 ~ GenomicGC + OptimumTemp + Habitat + OxygenReq + ProteinCount + GenomeSizeMb + GenomeDensity"


block <- testRandomTreeParams( modelFormula, data=traits[ ( !is.na(traits$Profile.26) & !is.na(traits$GenomicGC)  & !is.na(traits$Habitat) & !is.na(traits$OptimumTemp) & !is.na(traits$OxygenReq) & !is.na(traits$ProteinCount) & !is.na(traits$GenomeSizeMb) & !is.na(traits$GenomeDensity) ), ], N=N, M=M )

dev.off()
quit()


stats <- performRegressionWithRepeatedCrossVal( modelFormula, data=traits[ ( !is.na(traits$Profile.26) & !is.na(traits$GenomicGC)  & !is.na(traits$Habitat) & !is.na(traits$OptimumTemp) & !is.na(traits$OxygenReq) & !is.na(traits$ProteinCount) & !is.na(traits$GenomeSizeMb) & !is.na(traits$GenomeDensity) ), ], N=N, M=M, nodesize=10.5, mtry=7, ntree=2000 )

#--------------
pearsonByIter <- aggregate( stats$Pearson, by=list(stats$Iter), mean )
colnames(pearsonByIter) <- c("Iter", "MeanPearson")

stopifnot(nrow(stats)==M*N)
stopifnot(nrow(pearsonByIter)==M)


#p <- ggplot( data=pearsonByIter, aes(y=MeanPearson, x=1) ) +
#    geom_boxplot(outlier.size=0, outlier.alpha=0) +
#    geom_jitter(color="#777777", alpha=0.5, size=0.8)
#print(p)

p <- ggplot( data=pearsonByIter, aes(x=MeanPearson) ) +
    geom_histogram()
print(p)

#--------------
mseByIter <- aggregate( stats$MSE, by=list(stats$Iter), mean )
colnames(mseByIter) <- c("Iter", "MSE")

stopifnot(nrow(stats)==M*N)
stopifnot(nrow(mseByIter)==M)


p <- ggplot( data=mseByIter, aes(y=MSE, x=1) ) +
    geom_boxplot(outlier.size=0, outlier.alpha=0) +
    geom_jitter(color="#777777", alpha=0.5, size=0.8)
print(p)

p <- ggplot( data=mseByIter, aes(x=MSE) ) +
    geom_histogram()
print(p)

#--------------
maeByIter <- aggregate( stats$MAE, by=list(stats$Iter), mean )
colnames(maeByIter) <- c("Iter", "MAE")

stopifnot(nrow(stats)==M*N)
stopifnot(nrow(maeByIter)==M)


#p <- ggplot( data=maeByIter, aes(y=MAE, x=1) ) +
#    geom_boxplot(outlier.size=0, outlier.alpha=0) +
#    geom_jitter(color="#777777", alpha=0.5, size=0.8)
#print(p)

p <- ggplot( data=maeByIter, aes(x=MAE) ) +
    geom_histogram()
print(p)


#rf <- randomForest( Profile.15 ~ GenomicGC + OptimumTemp + Salinity + Habitat + OxygenReq , data=traits, nodesize=0.02, mtry=6, ntree=100, na.action="na.omit" )
#rf <- randomForest( Profile.15 ~ GenomicGC , data=traits, nodesize=250, mtry=6, ntree=500 )
rf <- getRandomForestRegression( modelFormula, data=traits, nodesize=10.5, mtry=7, ntree=2000, na.action="na.omit" )
varImpPlot( rf )
print( importance( rf ) )

plot( rf )


out <- predict( rf, traits )
print(length(out))
print(sum(!is.na(out)))

comp <- data.frame( actual=traits[,"Profile.26"], predicted=out, GenomicGC=traits[,"GenomicGC"], OptimumTemp=traits[,"OptimumTemp"], GenomeDensity=traits[,"GenomeDensity"], OxygenReq=traits[,"OxygenReq"] )
p <- ggplot( data=comp ) +
    geom_point( aes(x=actual, y=predicted) ) +
    coord_fixed(ratio=1) +
    theme( aspect.ratio=1 )
print(p)

p <- ggplot( data=comp ) +
    geom_point( aes(x=GenomicGC, y=actual) )
print(p)

p <- ggplot( data=comp ) +
    geom_point( aes(x=GenomicGC, y=predicted) )
print(p)

p <- ggplot( data=comp ) +
    geom_point( aes(x=OptimumTemp, y=predicted) )
print(p)

p <- ggplot( data=comp ) +
    geom_point( aes(x=GenomeDensity, y=predicted) )
print(p)

p <- ggplot( data=comp ) +
    geom_point( aes(x=OxygenReq, y=predicted) )
print(p)



print(cor.test(traits[,"Profile.26"], out, method="pearson", na.action="na.omit" ) )

warnings()

dev.off()