# Perform GLS analysis for different traits and different taxonomic groups

library("Matrix")
library("nlme")
library("ggplot2")
library("reshape")
library("scales")
library("rhdf5")
library("rredis")
library("phylobase")
library("ape")
library("adephylo")
library("phytools")  # for read.newick
library("grid")
library("dplyr")

#---------------------------------------------------------------
# Configuration

# Source: http://www2.math.su.se/PATHd8/
# Created using: ~/src/PATHd8/PATHd8 ./data/nmicrobiol201648-s6.txt ./data/nmicrobiol201648-s6.txt.nw.PATHd8.nw
#inputTree <- "nmicro_s6_pruned_with_taxids.nw"
inputTree <- "test_tree.nw"   # TODO - verify this tree (the file is identical to nmicro_s6_pruned_with_taxids.nw)

redisConnect(host="power5", password="rnafold")

profileStart <- 0
profileStop <- 1000
profileStep <- 10
profileReference <- "begin"
profileLen <- (profileStop-profileStart)/profileStep


taxidToKingdomFilename <- "TaxidToKingdom.csv"   # create file using: python2 create_taxid_kingdom_table.py
#taxonTreeFilename <- "TaxidToKingdom.nw"   # create file using: python2 create_taxid_kingdom_table.py

groupsTableOutputFile <- "tree_traits_effects_analysis_with_taxgroups.out.csv"

pyramidLength <- 30+1

significanceLevel = 1e-2

#minimalTaxonSize <- 9 # Note - should match the value in create_taxid_kingdom_table.py
minimalTaxonSize <- 9 # Testing only

#profileMode <- "nativeLFE"
#profileMode <- "shuffledLFE"
profileMode <- "dLFE"


cairo_pdf(sprintf("tree_phenotypes_regression.out.%s.by.taxgroups.%%d.pdf", profileMode))


speciesBlacklist <- c("470", "1280", "4932", "2850")  # Filter species that will prevent analysis from being performed from the tree because of errors or missing data
# TODO - try to fix these...

#---------------------------------------------------------------

getH5Filename <- function( taxId )
{
    sprintf("gcdata_v2_taxid_%d_profile_%d_%d_%s_%d.h5", taxId, profileStop, profileStep, profileReference, profileStart)
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
        val <- factor(c(val), habitatLevels)
        return( val )
    }
    else
    {
        return( factor(c(NA), habitatLevels) )  # convert NULL to NA
    }
}

oxygenLevels = c(NA, "Aerobic", "Facultative", "Microaerophilic", "Anaerobic")
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


readAllProfiles <- function( taxIds )
{

    combined <- data.frame( GenomicGC=double(), OptimumTemp=double(), TemperatureRange=ordered(c(), temperatureLevels), PairedFraction=double(), Profile=matrix( rep(0.0, profileLen), nrow=0, ncol=profileLen), Salinity=ordered(c(), salinityLevels), Habitat=factor(c(), habitatLevels), OxygenReq=factor(c(), oxygenLevels), row.names=integer() )

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


#=========================================================================================
#================= Add taxon membership indicators to the traits dataset =================
#=========================================================================================

# Read the taxon memberships table
taxidToKingdom <- read.table(taxidToKingdomFilename, sep=",", header=TRUE )    # create file using: python2 create_taxid_kingdom_table.py
taxidToKingdom$Member_all_1 <- 1   # Add a special group which includes all species

allTaxIds = getAllTaxIds() # get list of all taxids in the current dataset

# Read profiles for all species
traits <- readAllProfiles( allTaxIds )
traits$tax_id <- rownames(traits)  # create an explicit tax_id column (rather than an index), by which we'll be able to merge (since left_join can't merge using indices)

taxidToKingdom$tax_id <- as.character(taxidToKingdom$tax_id)  # Convert the taxids column in the taxon memberships table to string, to match traits table

traits <- traits %>% left_join(taxidToKingdom, by="tax_id")  # Perform the merge
rownames(traits) <- traits$tax_id  # Restore the tax_ids index


stopifnot(nrow(traits) > 300 )  # sanity test
#stopifnot(ncol(traits) == profileLen+7 )


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

#vcv1 <- vcv( phy=tree, model="Brownian")
#dimnames(vcv1) <- list(1:ncol(vcv1))
#bmcorrmtx <- data.frame( melt( vcv1 ) )
#p <- ggplot( data=bmcorrmtx, aes(x=X1, y=X2, fill=value) ) + geom_raster() + scale_y_reverse(); p


#----------------------------------


performGLSregression <- function( traits, tree, Xtrait, Ytrait, plotRegression=TRUE, caption="" )
{
    # ==========================================================================================
    # ============================== Part 1 - Prapre Tree ======================================
    # ==========================================================================================
    # Discard tree tips for species missing data
    speciesWithMissingData <- row.names(traits[(is.na(traits[Xtrait]) | is.na(traits[Ytrait])),])
    tree <- drop.tip( tree, speciesWithMissingData )   # Filter species that will prevent analysis from being performed from the tree

    # Discard trait data for species missing from the tree
    treeSpecies <- tree$tip.label
    traits <- traits[treeSpecies,]   # Discard traits not found in the tree

    if( is.null(tree) )
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
    bmcorr <- corBrownian( phy=tree )

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

        normYtrait <- paste(Ytrait, ".norm", sep="" )
        traits[,normYtrait] <- traits[,Ytrait] - mean(traits[,Ytrait])
        Ytrait <- normYtrait   # proceed in analysis using the normalized trait

    }
    
    # Set the predictive (regression) and intecept-only formulas
    # 'predictiveFormula' is the actual regression. It uses an intercept when a continuous explanatory var is used (and none with a discrete var, as in one-way ANOVA)
    # 'inteceptFormula' is an intecept-only model, used as the baseline for calculation of R^2
    predictiveFormula <- NA    
    if( any(class(traits[,Xtrait])==c("factor")) )
    {
        predictiveFormula <- as.formula( paste(Ytrait, " ~ ", Xtrait, " + 0") )
    }
    else
    {
        predictiveFormula <- as.formula( paste(Ytrait, " ~ ", Xtrait) )
    }
    interceptFormula = as.formula( paste(Ytrait, " ~ ", "1") )
    
    
    # Perform OLS regression (used as a reference for comparison)
    m1 <- lm(  predictiveFormula, traits); summary(m1)
    co <- coef(m1)

    # Perform GLS regression using the predictive and intercept-only models
    gls1 <- gls( predictiveFormula,
                traits,
                correlation=bmcorr,
                na.action=na.omit,
                method="REML"); summary(gls1)
    co2 <- coef(gls1)

    gls0 <- gls( interceptFormula,
                 traits,
                 correlation=bmcorr,
                 na.action=na.omit,
                 method="REML")
    

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

    # Slop direction
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
              geom_jitter(colour="grey") +
              geom_text( data=significanceMarkers, aes(x=x, y=y, label=var),       colour="gray",   size=2, alpha=0.7 ) +
              geom_text( data=significanceMarkers, aes(x=x, y=y, label=pval.char), colour="black",  size=2  ) +
              geom_text( data=significanceMarkers, aes(x=x, y=y, label=label),     colour="black",  size=10 ) +
              geom_hline( yintercept = 0 ) +
              annotate( "text", x=Inf,  y=max.y+c(0, 1, 2, 3) * line.y, label=c(  "lm",  sprintf("italic(R) ^ 2 == %.3g", summary(m1)$r.squared), sprintf('italic(p)*"-val" ==  %.3g', pvalue.OLS), sprintf("italic(N) == %d", nrow(traits))), hjust=1, vjust=0, colour="red", parse=TRUE ) +
              annotate( "text", x=-Inf, y=max.y+c(0, 1, 2, 3) * line.y, label=c(  "gls", sprintf("italic(R) ^ 2 == %.3g", R2),                    sprintf('italic(p)*"-val" ==  %.3g', pvalue), sprintf("italic(N) == %d", nrow(traits))),     hjust=0, vjust=0, colour=clr2,  parse=TRUE )
              print(p)
        }
        else
        {
            p <- ggplot(traits, aes(get(Xtrait), get(Ytrait))) +
              labs(y=Ytrait, x=Xtrait, title=caption) +
              geom_point() +
              geom_hline( yintercept = 0 ) +
              geom_abline( aes( slope=co[Xtrait],  intercept=co["(Intercept)"]),  colour="red"  ) +
              geom_abline( aes( slope=co2[Xtrait], intercept=co2["(Intercept)"]), colour=clr2   ) +
              annotate( "text", x=Inf,  y=max.y+c(0, 1, 2, 3) * line.y, label=c(  "lm",  sprintf("italic(R) ^ 2 == %.3g", summary(m1)$r.squared), sprintf('italic(p)*"-val" == %.3g', pvalue.OLS), sprintf("italic(N) == %d", nrow(traits))), hjust=1, vjust=0, colour="red", parse=TRUE ) +
              annotate( "text", x=-Inf, y=max.y+c(0, 1, 2, 3) * line.y, label=c(  "gls", sprintf("italic(R) ^ 2 == %.3g", R2),                    sprintf('italic(p)*"-val" == %.3g', pvalue),     sprintf("italic(N) == %d", nrow(traits))), hjust=0, vjust=0, colour=clr2,  parse=TRUE )
              print(p)
        }
    }

    return( data.frame( pvalue=pvalue, R2=R2 * slopeDirection, N=N, directionsIndicator=directionsIndicator ) )
}

performGLSregression_profileRangeMean <- function( traits, tree, Xtrait, profileRange )
{
    #print("performGLSregression_profileRangeMean(): -->")
    #print(profileRange)
    #print(class(profileRange))
    stopifnot(class(profileRange)=="integer")
    stopifnot(length(profileRange)==2)
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

    # Perform the regression (for the mean values)
    return( performGLSregression( traits, tree, Xtrait, "RangeMean", plotRegression=FALSE ) )
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
plotRotatedPyramid <- function( plotdf, dotsDF, valLabel, valScale, Xtrait, plotTitle )
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
         geom_point( data=dotsDF, aes(x=Var1,y=Var2), position=position_nudge(y=4, x=-4), color="white", size=0.6 ) ,  # Significance markers (white dots)
         theme( plot.margin = unit(c(0.1,0.1,0.1,0.1), "npc"), # Set the margins (should match the calculation used to create the "diagonal" X axis
               plot.background = element_blank(),   # Hide unnecessary theme elements (background panels, etc.)
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(color="black"),
               aspect.ratio=1 )  # Fix 1:1 aspect ratio
         ) )
    
    #----------------------------------------------------------------------
    # Draw heat-map pyramid on canvas
    
    grid.newpage()
    pushViewport( viewport( x=unit(0.5, "npc"), y=unit(0.15, "npc"), angle=-45 ) )

    pyramidPlotGrob <- ggplotGrob( pyramidPlot + plot.elements )
    grid.draw( pyramidPlotGrob )

    
    uu1a <- convertX( unit(1, "native"), "mm" )
    #print(uu1a)
    uu1b <- convertY( unit(1, "native"), "mm" )
    #print(uu1b)
        
    uw = 0.12  # triangle offset
    grid.polygon( x=c(0+uw, 1.0, 1.0, 0+uw), y=c(0.0, 1-uw, 0.0, 0.0), default.units="native", gp=gpar(fill="white", col=NA ) )

    
    #----------------------------------------------------------------------
    # Draw X-axis

    ww2 <- convertUnit(uu1b * 0.83 * sqrt(2.0), "npc")  # Todo - this returns more than 1 npcs (screen width). Why?
    #print(ww2)
    
    #u2 = 0.05
    u2 = 0.05
    #vp2 <- viewport( x=unit(0.5+u2, "native"), y=unit(0.5-u2, "native"), h=unit(0.01, "native"), w=ww2, angle=45 )
    vp2 <- viewport( x=unit(0.5+u2, "native"), y=unit(0.5-u2, "native"), h=unit(0.01, "native"), w=unit(0.9, "native"), angle=45 )
    pushViewport( vp2 )
    
    #grid.lines( c(0.0, 1.0), c(0.5, 0.5) )

    xaxis <- xaxisGrob( at=seq(0,1,5/(pyramidLength-1)) )

    xaxis <- editGrob(xaxis,
                      gPath("major"),
                      x=unit(c(0,1), "npc")
                      )             
    xaxis <- editGrob(xaxis,        
                      gPath("labels"),
                      label=seq(0,(pyramidLength-1)*10,50)
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
        grid.text(sprintf("min==%.3g", min(plotdf[,"value"])), x=unit(0.5, "npc"), y=unit(0.5, "npc"), just="right" )
        grid.text(sprintf("max==%.3g", max(plotdf[,"value"])), x=unit(0.5, "npc") - unit(0, "points"), y=unit(0.5, "npc") + unit(15, "points"), just="right" )
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

glsRegressionRangeAnalysis <- function( traits, tree, Xtrait, profileRange, plotCaption )
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
        regressionResults <- performGLSregression_profileRangeMean(traits, tree, Xtrait, as.vector(v, mode="integer") )
        v$Pvalue    <- regressionResults$pvalue
        v$Buse.R2   <- regressionResults$R2
        v$Logpval <- log10( v$Pvalue )
        v$Var1 <- (v$Var1 - 1) * 10
        v$Var2 <- (v$Var2 - 1) * 10
        v$NumSpecies <- regressionResults$N
        v$DirectionsIndicator <- regressionResults$directionsIndicator

        #print("regressionResults")
        #print(regressionResults)
        #print("out")
        #print(out)
        #print("v")
        #print(v)
        out <- rbind( out, v )
        assign('out', out, parent.frame() )
        
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


    plotTitle = sprintf("%s effect on %s (%s)", Xtrait, profileMode, plotCaption)
    
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

    plotdf <- melt( subranges, id.vars=c("Var1", "Var2"), measure.vars=c("Logpval") )
    
    plotRotatedPyramid( plotdf, dotsDF, "log(p-val)", valScale, Xtrait, paste(plotTitle, " (logPval)") )

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
    
    plotRotatedPyramid( plotdf, dotsDF, "Buse R^2", valScale, Xtrait, plotTitle=paste(plotTitle, " (Buse R^2)") )


    #----------------------------------------------
    
    plotdf <- data.frame( melt( subranges, id.vars=c("Var1", "Var2"), measure.vars=c("DirectionsIndicator"   ) ) )
    #print(levels(plotdf$value))
    #levels(plotdf$value) <- as.factor(1:nlevels(plotdf$value))
    #print(levels(plotdf$value))
    #print(plotdf)

    #plotdf$value <- as.factor(plotdf$value)

    #valScale <- scale_colour_manual( breaks=c("1","2"), labels=c("1","2"), values=c("#7799f0", "#f06666") )
    valScale <- scale_fill_hue( ) # plotdf, aes(fill=factor(value) ) )

    #print("-------3-------")

    #print(class(plotdf$value))
    plotRotatedPyramid( plotdf, dotsDF, "Directions", valScale, Xtrait, plotTitle=paste(plotTitle, " (order)") )
    
    
    return(subranges)
}

#performGLSregression( traits, tree, "GenomicGC",   "Profile.15" )
#print(performGLSregression_profileRangeMean( traits, tree, "GenomicGC", c(15,15) ))
#print(performGLSregression_profileRangeMean( traits, tree, "GenomicGC", c(15,19) ))


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


prepareFilteredTreeForGLS <- function( traits, tree, filterTrait, studyTrait )
{
    stopifnot( any( filterTrait == colnames(traits) ) )  # filterTrait not found
    stopifnot( any( studyTrait  == colnames(traits) ) )  # filterTrait not found
        
    # Filter by the requested trait, and also discard tree tips for species missing data
    speciesWithMissingData <- row.names(traits[ (!traits[,filterTrait] | is.na(traits[,studyTrait])), ])
    tree <- drop.tip( tree, speciesWithMissingData )   # Filter species that will prevent analysis from being performed from the tree

    if( is.null(tree) )
    {
        return( c(NA, NA ) )
    }
    
    # Discard trait data for species missing from the tree
    treeSpecies <- tree$tip.label
    traits <- traits[treeSpecies,]   # Discard traits not found in the tree

    if( is.null(traits) || nrow(traits) < 3 )  # no sense in analyzing a tiny tree...
    {
        return( c(NA, NA) )
    }

    # Check tree <-> data-frame correlation appears valid
    x <- phylo4d(tree, tip.data=traits, rownamesAsLabels=TRUE)
    tipsOk <- hasTipData( x )
    stopifnot(tipsOk)
    
    return( list( traits, tree ) )
}




performGLSregressionWithFilter <- function( traits, tree, filterTrait, Xtrait, Ytrait, plotRegression=TRUE )
{
    stopifnot( any( filterTrait == colnames(traits) ) )  # filterTrait not found
    stopifnot( any( Xtrait      == colnames(traits)  ) )  # studyTrait not found

    # Filter tree by trait
    ret <- prepareFilteredTreeForGLS( traits, tree, filterTrait, Xtrait )
    traits <- ret[[1]]
    tree   <- ret[[2]]
        
    if( is.null(tree) )
    {
        return( NA )
    }

    N <- nTips(tree)
    
    return( performGLSregression( traits, tree, Xtrait, Ytrait, plotRegression=TRUE, caption=sprintf("%s effect on %s (%s)", Xtrait, profileMode, filterTrait) ) )
}

glsRangeAnalysisWithFilter <- function( traits, tree, filterTrait, studyTrait, pyramidSpec )
{
    stopifnot( any( filterTrait==colnames(traits) ) )  # filterTrait not found
    stopifnot( any( studyTrait==colnames(traits)  ) )  # studyTrait not found

    # Filter tree by trait
    ret <- prepareFilteredTreeForGLS( traits, tree, filterTrait, studyTrait )
    if( is.na(ret[[1]]) ) { return( NA ) }
    traits <- ret[[1]]
    tree   <- ret[[2]]
        
    if( is.null(tree) || nTips(tree) < minimalTaxonSize )
    {
        return( NA )
    }

    N <- nTips(tree)
    
    return( glsRegressionRangeAnalysis( traits, tree, studyTrait, pyramidSpec, sprintf("%s\n(N=%d)", filterTrait, N ) ) )
}

# These are the names the of taxon membership variables
#taxGroups <- c('Member_all_1', 'Member_Bacteria_2', 'Member_Terrabacteria_group_1783272', 'Member_Proteobacteria_1224', 'Member_Eukaryota_2759', 'Member_Archaea_2157', 'Member_Gammaproteobacteria_1236', 'Member_Firmicutes_1239', 'Member_FCB_group_1783270', 'Member_Euryarchaeota_28890', 'Member_Bacteroidetes_Chlorobi_group_68336', 'Member_Opisthokonta_33154', 'Member_Bacteroidetes_976', 'Member_Actinobacteria_201174', 'Member_Bacilli_91061', 'Member_Actinobacteria_1760', 'Member_Flavobacteriia_117743', 'Member_Flavobacteriales_200644', 'Member_Fungi_4751', 'Member_unclassified_Bacteria_2323', 'Member_Alphaproteobacteria_28211', 'Member_Dikarya_451864', 'Member_Bacteria_candidate_phyla_1783234', 'Member_Flavobacteriaceae_49546', 'Member_Patescibacteria_group_1783273', 'Member_Bacillales_1385', 'Member_TACK_group_1783275', 'Member_Parcubacteria_group_1794811', 'Member_Ascomycota_4890', 'Member_PVC_group_1783257', 'Member_Cyanobacteria_1117', 'Member_Cyanobacteria_Melainabacteria_group_1798711', 'Member_Chloroflexi_200795', 'Member_saccharomyceta_716545', 'Member_Rhizobiales_356', 'Member_Deinococci_188787', 'Member_Enterobacterales_91347', 'Member_Deinococcus_Thermus_1297', 'Member_delta_epsilon_subdivisions_68525', 'Member_Thermoprotei_183924', 'Member_Tenericutes_544448', 'Member_Crenarchaeota_28889', 'Member_Mollicutes_31969', 'Member_Thermotogae_188708', 'Member_Stramenopiles_33634', 'Member_Aquificae_200783', 'Member_Thermotogae_200918', 'Member_Aquificae_187857', 'Member_Enterobacteriaceae_543', 'Member_Thermococcales_2258', 'Member_Thermococcaceae_2259', 'Member_Micrococcales_85006', 'Member_Thermococci_183968', 'Member_Basidiomycota_5204', 'Member_Clostridia_186801', 'Member_Bacillaceae_186817', 'Member_Lactobacillales_186826' )

taxGroups <- colnames(taxidToKingdom)
taxGroups <- taxGroups[startsWith(taxGroups, "Member_")]

taxGroupToTaxId <- vapply(taxGroups, function (x) as.integer(tail(strsplit(x[[1]], '_')[[1]], n=1)), integer(1) )
# Example: taxGroupToTaxId['Member_Cyanobacteria_1117'] == 1117


#majorGroupsTree <- read.newick(taxonTreeFilename)


performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",                  "OxygenReq",   "Profile.7" )
performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",                  "OxygenReq",   "Profile.14" )
performGLSregressionWithFilter( traits, tree, "Member_Terrabacteria_group_1783272", "OxygenReq",   "Profile.8" )

performGLSregressionWithFilter( traits, tree, "Member_Proteobacteria_1224",         "OxygenReq",   "Profile.16" )

performGLSregressionWithFilter( traits, tree, "Member_Gammaproteobacteria_1236",    "OxygenReq",   "Profile.19" )
performGLSregressionWithFilter( traits, tree, "Member_Gammaproteobacteria_1236",    "OxygenReq",   "Profile.20" )

performGLSregressionWithFilter( traits, tree, "Member_Actinobacteria_1760", "OxygenReq",   "Profile.9" )
performGLSregressionWithFilter( traits, tree, "Member_Actinobacteria_1760", "OxygenReq",   "Profile.14" )
performGLSregressionWithFilter( traits, tree, "Member_Actinobacteria_1760", "OxygenReq",   "Profile.20" )
performGLSregressionWithFilter( traits, tree, "Member_Terrabacteria_group_1783272", "OxygenReq", "Profile.16" )
performGLSregressionWithFilter( traits, tree, "Member_Terrabacteria_group_1783272", "OxygenReq", "Profile.19" )

#performGLSregressionWithFilter( traits, tree, "Member_Terrabacteria_group_1783272", "GenomicGC", "Profile.15" )
#performGLSregressionWithFilter( traits, tree, "Member_Terrabacteria_group_1783272", "Habitat",   "Profile.15" )

performGLSregressionWithFilter( traits, tree, "Member_all_1", "GenomicGC", "Profile.15" )

performGLSregressionWithFilter( traits, tree, "Member_Bacteroidetes_Chlorobi_group_68336", "GenomicGC", "Profile.22" )
performGLSregressionWithFilter( traits, tree, "Member_Bacteroidetes_Chlorobi_group_68336", "GenomicGC", "Profile.25" )

performGLSregressionWithFilter( traits, tree, "Member_Chloroflexi_200795", "GenomicGC", "Profile.18" )
performGLSregressionWithFilter( traits, tree, "Member_Chloroflexi_200795", "GenomicGC", "Profile.26" )
performGLSregressionWithFilter( traits, tree, "Member_Chloroflexi_200795", "GenomicGC", "Profile.27" )

performGLSregressionWithFilter( traits, tree, "Member_Fungi_4751", "GenomicGC", "Profile.8" )
performGLSregressionWithFilter( traits, tree, "Member_Fungi_4751", "GenomicGC", "Profile.20" )

performGLSregressionWithFilter( traits, tree, "Member_Proteobacteria_1224",         "OptimumTemp",   "Profile.10" )
performGLSregressionWithFilter( traits, tree, "Member_Proteobacteria_1224",         "OptimumTemp",   "Profile.17" )

performGLSregressionWithFilter( traits, tree, "Member_Archaea_2157",         "OptimumTemp",   "Profile.11" )
performGLSregressionWithFilter( traits, tree, "Member_Archaea_2157",         "OptimumTemp",   "Profile.20" )


performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",                  "Habitat",   "Profile.1" )
performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",                  "Habitat",   "Profile.2" )
performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",                  "Habitat",   "Profile.3" )
performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",                  "Habitat",   "Profile.4" )
performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",                  "Habitat",   "Profile.5" )
performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",                  "Habitat",   "Profile.6" )
performGLSregressionWithFilter( traits, tree, "Member_Bacteria_2",                  "Habitat",   "Profile.27" )


# TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY #
#dev.off()
#quit()
# TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY #


regressionResultsByTaxGroup <- data.frame( ExplanatoryVar=character(), MaxRangeStart=integer(), MaxRangeEnd=integer(), TaxGroup=integer(), TaxGroupName=character(), EffectSize=double(), Pvalue=double(), NumSpecies=integer() )


for( gr in taxGroups )
{
    print(gr)
    results <- glsRangeAnalysisWithFilter( traits, tree, gr, "GenomicGC", c(1,pyramidLength))

    # Store the summarized result for this group
    if( !is.na(results) )
    {
        results <- results[results$Var1==results$Var2,]
        maxResult <- results [which.max( abs(results$Buse.R2) ),]

        regressionResultsByTaxGroup <- rbind( regressionResultsByTaxGroup, data.frame( ExplanatoryVar=c("GenomicGC"), MaxRangeStart=c(maxResult$Var1), MaxRangeEnd=c(maxResult$Var2), TaxGroup=c(taxGroupToTaxId[gr]), TaxGroupName=c(gr), EffectSize=c(maxResult$Buse.R2), Pvalue=c(maxResult$Pvalue), NumSpecies=c(maxResult$NumSpecies) ) )
    }
}

print(regressionResultsByTaxGroup)

for( gr in taxGroups )
{
    print(gr)
    results <- glsRangeAnalysisWithFilter( traits, tree, gr, "OptimumTemp", c(1,pyramidLength))

    # Store the summarized result for this group
    if( !is.na(results) )
    {
        results <- results[results$Var1==results$Var2,]
        maxResult <- results[which.max( abs(results$Buse.R2) ),]

        regressionResultsByTaxGroup <- rbind( regressionResultsByTaxGroup, data.frame( ExplanatoryVar=c("OptimumTemp"), MaxRangeStart=c(maxResult$Var1), MaxRangeEnd=c(maxResult$Var2), TaxGroup=c(taxGroupToTaxId[gr]), TaxGroupName=c(gr), EffectSize=c(maxResult$Buse.R2), Pvalue=c(maxResult$Pvalue), NumSpecies=c(maxResult$NumSpecies) ) )
    }
}

# TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY #
#dev.off()  # TESTING ONLY
#write.csv(regressionResultsByTaxGroup, file=groupsTableOutputFile )    # TESTING ONLY
#quit()    # TESTING ONLY
# TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY ####  TESTING ONLY #


print(regressionResultsByTaxGroup)

#pdf("tree_phenotypes_regression.Habitat.out.pdf")
for( gr in taxGroups )
{
    print(gr)
    results <- glsRangeAnalysisWithFilter( traits, tree, gr, "Habitat", c(1,pyramidLength))
    
    # Store the summarized result for this group
    if( !is.na(results) )
    {
        results <- results[results$Var1==results$Var2,]
        maxResult <- results[which.max( abs(results$Buse.R2) ),]

        regressionResultsByTaxGroup <- rbind( regressionResultsByTaxGroup, data.frame( ExplanatoryVar=c("Habitat"), MaxRangeStart=c(maxResult$Var1), MaxRangeEnd=c(maxResult$Var2), TaxGroup=c(taxGroupToTaxId[gr]), TaxGroupName=c(gr), EffectSize=c(maxResult$Buse.R2), Pvalue=c(maxResult$Pvalue), NumSpecies=c(maxResult$NumSpecies) ) )
    }
}
#dev.off()

print(regressionResultsByTaxGroup)

#pdf("tree_phenotypes_regression.Salinity.out.pdf")
for( gr in taxGroups )
{
    print(gr)
    results <- glsRangeAnalysisWithFilter( traits, tree, gr, "Salinity", c(1,pyramidLength))
    
    # Store the summarized result for this group
    if( !is.na(results) )
    {
        results <- results[results$Var1==results$Var2,]
        maxResult <- results[which.max( abs(results$Buse.R2) ),]

        regressionResultsByTaxGroup <- rbind( regressionResultsByTaxGroup, data.frame( ExplanatoryVar=c("Salinity"), MaxRangeStart=c(maxResult$Var1), MaxRangeEnd=c(maxResult$Var2), TaxGroup=c(taxGroupToTaxId[gr]), TaxGroupName=c(gr), EffectSize=c(maxResult$Buse.R2), Pvalue=c(maxResult$Pvalue), NumSpecies=c(maxResult$NumSpecies) ) )
    }
}
#dev.off()

#pdf("tree_phenotypes_regression.Oxygen.out.pdf")
for( gr in taxGroups )
{
    print(gr)
    results <- glsRangeAnalysisWithFilter( traits, tree, gr, "OxygenReq", c(1,pyramidLength))

    # Store the summarized result for this group
    if( !is.na(results) )
    {
        results <- results[results$Var1==results$Var2,]
        maxResult <- results[which.max( abs(results$Buse.R2) ),]

        regressionResultsByTaxGroup <- rbind( regressionResultsByTaxGroup, data.frame( ExplanatoryVar=c("OxygenReq"), MaxRangeStart=c(maxResult$Var1), MaxRangeEnd=c(maxResult$Var2), TaxGroup=c(taxGroupToTaxId[gr]), TaxGroupName=c(gr), EffectSize=c(maxResult$Buse.R2), Pvalue=c(maxResult$Pvalue), NumSpecies=c(maxResult$NumSpecies) ) )
    }
}


dev.off()

# Save the summarized results, by group, to csv file (plot using: python2 plot_tree_effects_anaylsis_results_on_tree.py)
write.csv(regressionResultsByTaxGroup, file=groupsTableOutputFile )

warnings()
