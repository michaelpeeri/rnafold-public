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
library("grid")

#---------------------------------------------------------------
# Configuration

cairo_pdf("tree_phenotypes_regression.out.%d.pdf")
#svg("tree_phenotypes_regression.out.%d.svg")

redisConnect(host="power5", password="rnafold")

profileStart <- 0
profileStop <- 1000
profileStep <- 10
profileReference <- "begin"
profileLen <- (profileStop-profileStart)/profileStep

pyramidLength <- 30+1

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
    #print(co)

    #print("//////////////")
    #print( nrow(traits) )
    #print( nTips(traitsTree) )

    gls1 <- gls( YvsX,
                 traits,
                 correlation=bmcorr,
                 na.action=na.omit,
                 method="REML"); summary(gls1)
    #gls1 <- pGLS( YvsX, traits, bmcorr, na.action=na.omit);
    print("--")
    print(summary(gls1))
    print(summary(gls1)$tTable[2,4])
    co2 <- coef(gls1)

    gls0 <- gls( as.formula( paste(Ytrait, " ~ ", "1")),
                 traits,
                 correlation=bmcorr,
                 na.action=na.omit,
                 method="REML")


    #Ydata <- traits[,Ytrait]
    #Xdata <- traits[,Xtrait]
    #tss <- sum( (Ydata - mean(Ydata))**2 )
    #var.from.tss <- tss / length(Ydata)
    ##print(var(Ydata))
    ##print(var.from.tss)
    #rss <- sum( residuals(gls1)**2 )
    #Rsquared <- 1 - rss/tss

    #prediction.cor <- cor( Ydata, predict(gls1, Xdata), method="pearson")
    prediction.cor <- cor( traits[,Ytrait], fitted(gls1), method="pearson")
    Rsquared.2 <- prediction.cor ** 2
    #print(Rsquared)
    #print(Rsquared.2)

    #logP.null <- logLik(gls0, REML=FALSE)
    #print( logP.null )
    #stopifnot( logP.null <= 0.0 )
    #logP.fitted <- logLik(gls1, REML=FALSE)
    #stopifnot( logP.fitted <= 0.0 )

    #D.null.minus.D.fitted = -2 * ( logP.null - logP.fitted )
    #D.null =                -2 * ( logP.null - 0.0         )
    
    #pseudoR2 <- 1 - as.numeric( logP.fitted / logP.null )
    #print("psi")
    #print(pseudoR2)
    #print(summary(m1)$r.squared)

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
    e <- rep(1, 2)
    Y <- traits[,Ytrait]
    Y.bar <- coef(gls0)["(Intercept)"]
    yye <- Y - Y.bar * e
    
    R2 <- 1 - ( t(u.hat) %*% inv.V %*% u.hat )/( t(yye) %*% inv.V %*% yye )   # See eq. (15) p. 107.
    stopifnot( R2 >= 0 && R2 <= 1.0 )

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

    if( plotRegression )
    {
        clr2 <- "#55CCB0"

        if( any(class(traits[,Xtrait])==c("factor")) )
        {
            p <- ggplot(traits, aes(get(Xtrait), get(Ytrait))) +
              labs(y=Ytrait, x=Xtrait) +
              geom_boxplot(outlier.size=0, outlier.alpha=0) +  # don't show outliers on boxplot, since we're plotting all data anyway...
              geom_jitter(colour="grey") +
              geom_hline( yintercept = 0 ) +
              annotate( "text", x=Inf,  y=Inf, label=sprintf("lm\nR^2 = %4g\nP-val = %4g\nn = %d",  summary(m1)$r.squared,      summary(m1)$coefficients[2,4], nrow(traits) ), hjust=1, vjust=1, colour="red" ) +
              annotate( "text", x=-Inf, y=Inf, label=sprintf("gls\nR^2 = %4g\nP-val = %4g\nn = %d", R2,                         summary(gls1)$tTable[2,4],     nrow(traits) ), hjust=0, vjust=1, colour=clr2 );
              print(p)
            fn <- sprintf("tree_phenotypes_regression.%s_vs_%s.out.pdf", Ytrait, Xtrait)
            ggsave(fn, plot=p )
            unlink(fn)
            #fn <- sprintf("tree_phenotypes_regression.%s_vs_%s.out.svg", Ytrait, Xtrait)
            #ggsave(fn, plot=p )
            #unlink(fn)
            #fn <- sprintf("tree_phenotypes_regression.%s_vs_%s.out.wmf", Ytrait, Xtrait)
            #ggsave(fn, plot=p )
            #unlink(fn)
        }
        else
        {
            p <- ggplot(traits, aes(get(Xtrait), get(Ytrait))) +
              labs(y=Ytrait, x=Xtrait) +
              geom_point() +
              geom_hline( yintercept = 0 ) +
              geom_abline( aes( slope=co[Xtrait],  intercept=co["(Intercept)"]), colour="red"  ) +
              geom_abline( aes( slope=co2[Xtrait], intercept=co2["(Intercept)"]), colour=clr2  ) +
              #geom_linerange( aes(min=get(Ytrait), max=fitted(gls1)), color="#aaaaaa80" ) +   # show residuals
              annotate( "text", x=Inf,  y=Inf, label=sprintf("lm\nR^2 = %4g\nP-val = %4g\nn = %d", summary(m1)$r.squared,      summary(m1)$coefficients[2,4], nrow(traits) ), hjust=1, vjust=1, colour="red" ) +
              annotate( "text", x=-Inf, y=Inf, label=sprintf("gls\nBuse R^2 = %4g\nP-val = %4g\nn = %d", R2, summary(gls1)$tTable[2,4],     nrow(traits) ), hjust=0, vjust=1, colour=clr2 );
              print(p)
        }
    }

    #return( c( summary(gls1)$tTable[2,4], co2[Xtrait] ) )
    #return( c( summary(gls1)$tTable[2,4],  summary(m1)$r.squared ) )
    #if( Rsquared >= 0.0 && Rsquared <= 1.0 )
    #{
        return( c( summary(gls1)$tTable[2,4],  R2 * slopeDirection, Rsquared.2, summary(m1)$r.squared ) )
    #}
    #else
    #{
    #    return( c( NA, NA ) )
    #}
        
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
plotRotatedPyramid <- function( plotdf, dotsDF, valScale, Xtrait, plotTitle )
{
    pyramidPlot <- ggplot( data=plotdf, aes(x=Var1, y=Var2, fill=value ) ) +
         geom_tile( linetype=0, colour=NA ) +
         labs(y="Profile End", x="Profile Start", value="log(p-value)", title=plotTitle) +
         valScale +
         coord_fixed() +
         guides( fill=FALSE ) +
         geom_point(data=dotsDF, aes(x=Var1,y=Var2), position=position_nudge(y=4, x=-4), color="white", size=0.6 ) +
         theme( plot.margin = unit(c(0.1,0.1,0.1,0.1), "npc"),
               plot.background = element_blank(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(color="black"),
               #legend.position="left",
               aspect.ratio=1 )


    #----------------------------------------------------------------------
    # Draw heat-map pyramid on canvas
    
    grid.newpage()
    pushViewport( viewport( x=unit(0.5, "npc"), y=unit(0.15, "npc"), angle=-45 ) )
    pyramidPlotGrob <- ggplotGrob( pyramidPlot )
    grid.draw( pyramidPlotGrob )
    #print(str(pyramidPlotGrob))

    #print( str( p ) )
    #print( str( ggplotGrob( p ) ) )
    #print( ggplot_build(p)$layout$panel_params[[1]]$x.range )
    #print( ggplot_build(p)$layout$panel_params[[1]]$y.range )



    #print("--------------------------------------")
    #gt <- ggplot_gtable(ggplot_build(p))
    #for( gr in gt$grobs )
    #{
    #    print(gr$name)
    #}
    
    #dev.off()
    #pdf("grobs.pdf")
    #for( i in 1:length(gt$grobs) ) {
    #    grid.draw(gt$grobs[[i]])
    #    grid.text(i, x=unit(0.1, "npc"), y=unit(0.1, "npc"))
    #    grid.newpage()
    #}
    #dev.off()
    #quit()

    #print("-----------------4---------------------")
    #print(attributes(gt$grobs[[3]]$children[1]$axis))
    #print(attributes(gt$grobs[[3]]$children[1]$axis$x))
    #print(attributes(gt$grobs[[3]]$children[1]$axis$y))

    #print("-----------------7---------------------")
    #print(attributes(gt$grobs[[7]]$children[1]$axis))
    #print(attributes(gt$grobs[[7]]$children[1]$axis$x))
    #print(attributes(gt$grobs[[7]]$children[1]$axis$y))

    #grid.ls()
    #grid.grep("axis")
    #current.vpTree()

    #print(getNames())
    #grid.ls(pyramidPlotGrob)

    #downViewport("axis-b.7-4-7-4")
    #grid.rect(gp=gpar(col=NA, fill=rgb(0,1,0,0.3)))
    #upViewport()

    #downViewport("axis-l.6-3-6-3")
    #grid.rect(gp=gpar(col=NA, fill=rgb(1,0,0.5,0.3)))

    #uu1a <- convertX( unit(1, "native"), "npc" )
    #print(uu1a)
    #uu1b <- convertY( unit(1, "native"), "npc" )
    #print(uu1b)
    uu1a <- convertX( unit(1, "native"), "mm" )
    print(uu1a)
    uu1b <- convertY( unit(1, "native"), "mm" )
    print(uu1b)
    
    #upViewport()

    #uu1 <- grobWidth( grid.get("axis-l.6-3-6-3" ) )
    #print(uu1)
    #uu2 <- grobHeight( grid.get("axis-l.6-3-6-3" ) )
    #print(uu2)


    #stop()
    #dev.off()
    #quit()

    print("--------------------------------------")
    #print( str( ggplot_build(p)$layout ) )
    #print("--------------------------------------")
    #print( str( ggplot_build(p)$layout$panel_ranges ) )
    #print("--------------------------------------")
    #print( str( ggplot_build(p)$layout$panel_scales ) )

    #----------------------------------------------------------------------
    # Draw white polygon to hide "unused" triangle
    
    uw = 0.12  # triangle offset
    grid.polygon( x=c(0+uw, 1.0, 1.0, 0+uw), y=c(0.0, 1-uw, 0.0, 0.0), default.units="native", gp=gpar(fill="white", col=NA ) )


    #print(convertUnit( x=unit(1.0, "native"), unitTo=
    #vptt <- viewport( x=unit(0.5, "native"), y=unit(0.5, "native"), h=unit(1, "native"), w=unit(1, "native") )
    #pushViewport( vptt )
    #grid.rect( gp=gpar(fill="yellow") )
    #upViewport()


    #----------------------------------------------------------------------
    # Draw X-axis

    ww2 <- convertUnit(uu1b * 0.83 * sqrt(2.0), "npc")  # Todo - this returns more than 1 npcs (screen width). Why?
    print(ww2)
    
    #u2 = 0.05
    u2 = 0.0
    #vp2 <- viewport( x=unit(0.5+u2, "native"), y=unit(0.5-u2, "native"), h=unit(0.01, "native"), w=ww2, angle=45 )
    vp2 <- viewport( x=unit(0.5+u2, "native"), y=unit(0.5-u2, "native"), h=unit(0.01, "native"), w=unit(0.9, "native"), angle=45 )
    pushViewport( vp2 )
    
    grid.lines( c(0.0, 1.0), c(0.5, 0.5) )

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

    p4 <- ggplot( data=plotdf, aes(x=Var1, y=Var2, fill=value ) ) +
         valScale +
#         scale_fill_gradient2(low="#cc7777", mid="#ffff88", high="#000000", limits=c(min(subranges$Logpval), 0), midpoint=minLogpval ) + 
         geom_tile()
    p4leg <- getLegend( p4 )
    grid.draw( p4leg )
    


    


    #print( p, vp=viewport(angle=-45) )
    fn <- sprintf("tree_phenotypes_regression.ranges.%s_vs_dLFE.out.pdf", Xtrait)
    ggsave(fn, plot=p )
    unlink(fn)
    #fn <- sprintf("tree_phenotypes_regression.ranges.%s_vs_dLFE.out.svg", Xtrait)
    #ggsave(fn, plot=p )
    #unlink(fn)
    #fn <- sprintf("tree_phenotypes_regression.ranges.%s_vs_dLFE.out.wmf", Xtrait)
    #ggsave(fn, plot=p )
    #unlink(fn)
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
    subranges$Pvalue     <- regressionResults[1,]
    subranges$Buse.R2    <- regressionResults[2,]
    subranges$CorSquared <- regressionResults[3,]
    subranges$OLS.R2     <- regressionResults[4,]
    
    subranges$Logpval <- log10( subranges$Pvalue )
    subranges$Var1 <- (subranges$Var1 - 1) * 10
    subranges$Var2 <- (subranges$Var2 - 1) * 10

    print(subranges)

    minLogpval <- min( subranges$Logpval, -3 )

    #----------------------------------------------------------------------
    # Create main pyramid heat-map

    plotTitle = sprintf("%s effect on dLFE", Xtrait)

    dotsDF <- data.frame(melt( subranges, id.vars=c("Var1", "Var2"), measure.vars=c("Logpval") ))
    dotsDF <- dotsDF[dotsDF$value<=-2.0,]

    plotdf <- melt( subranges, id.vars=c("Var1", "Var2"), measure.vars=c("Logpval") )
    #p <- ggplot( data=plotdf, aes(x=Var1, y=Var2, fill=value ) ) +
    #    geom_tile( linetype=0, colour=NA ) +
    #    labs(y="Profile End", x="Profile Start", value="log(p-value)", title=plotTitle) +
    #    coord_fixed() +
#   #     guides( fill=FALSE ) +
    #    geom_point(data=plotdf[plotdf$value<=-2.0,], aes(x=Var1,y=Var2), position=position_nudge(y=4, x=-4), color="white", size=0.6 ) +
    #    theme( plot.margin = unit(c(0.1,0.1,0.1,0.1), "npc"),
    #           plot.background = element_blank(),
    #           panel.grid.major = element_blank(),
    #           panel.grid.minor = element_blank(),
    #           panel.background = element_blank(),
    #           axis.line = element_line(color="black"),
    #           #legend.position="left",
    #           aspect.ratio=1 )
    valScale <- scale_fill_gradient(low="#7799ff", high="#000000", limits=c(min(subranges$Logpval), 0) )


    plotRotatedPyramid( plotdf, dotsDF, valScale, Xtrait, plotTitle )

    print("//////////////////////////////////")
    print( data.frame( melt( subranges, id.vars=c("Var1", "Var2"), measure.vars=c("Buse.R2"   ) ) ) )
    print("//////////////////////////////////")
    
    #grid.newpage()
    plotdf <- data.frame( melt( subranges, id.vars=c("Var1", "Var2"), measure.vars=c("Buse.R2"   ) ) )
    #p <- ggplot( data=plotdf, aes(x=Var1, y=Var2, fill=value) ) +
    #    geom_tile( linetype=0, colour=NA ) +
    #    labs(y="Profile End", x="Profile Start", value="Buse R^2", title=paste(plotTitle, " (Buse R^2)")) +
    #    coord_fixed() +
    #    geom_point(data=plotdf[plotdf$value<=-2.0,], aes(x=Var1,y=Var2), position=position_nudge(y=4, x=-4), color="white", size=0.6 ) +
    #    theme( plot.margin=unit(c(0.1,0.1,0.1,0.1), "npc"), aspect.ratio=1,
    #           plot.background = element_blank(),
    #           panel.grid.major = element_blank(),
    #           panel.grid.minor = element_blank(),
    #           panel.background = element_blank(),
    #           axis.line = element_line(color="black")
    #    )
    #v.min <- min(trunc(subranges$Buse.R2)-1, 0.0)
    #v.max <- max(trunc(subranges$Buse.R2)+1, 1.0)
    valScale <- scale_fill_gradientn(
        breaks=c(-1.0, -0.5, 0.0, 0.5, 1.0),
        limits=c(-1.0, 1.0),
        colors=        c("#75c0ff", "#2240cc", "#112060",   "#000000", "#602011", "#cc4022", "#ffc075"),
        values=rescale(c(     -1.0,      -0.5,     -0.15,         0.0,      0.15,       0.5,       1.0) ) )
    plotRotatedPyramid( plotdf, dotsDF, valScale, Xtrait, plotTitle=paste(plotTitle, " (Buse R^2)") )

    #grid.newpage()
    plotdf <- data.frame( melt( subranges, id.vars=c("Var1", "Var2"), measure.vars=c("CorSquared"   ) ) )
    #p <- ggplot( data=plotdf, , aes(x=Var1, y=Var2, fill=value) ) +
    #    geom_tile( linetype=0, colour=NA ) +
    #    labs(y="Profile End", x="Profile Start", value="CorSquared", title=paste(plotTitle, " (Cor^2)")) +
    #    coord_fixed() +
    #    geom_point(data=plotdf[plotdf$value<=-2.0,], aes(x=Var1,y=Var2), position=position_nudge(y=4, x=-4), color="white", size=0.6 ) +
    #    theme( plot.margin=unit(c(0.1,0.1,0.1,0.1), "npc"), aspect.ratio=1,
    #           plot.background = element_blank(),
    #           panel.grid.major = element_blank(),
    #           panel.grid.minor = element_blank(),
    #           panel.background = element_blank(),
    #           axis.line = element_line(color="black")
    #    )
    valScale <- scale_fill_gradient(low="#000000", high="#7799ff", limits=c(0, max(subranges$CorSquared, 1.0) ) )
    plotRotatedPyramid( plotdf, dotsDF, valScale, Xtrait, plotTitle=paste(plotTitle, " (Cor^2)") )

    #grid.newpage()
    plotdf <- data.frame( melt( subranges, id.vars=c("Var1", "Var2"), measure.vars=c("OLS.R2"   ) ) )
    #p <- ggplot( data=plotdf, , aes(x=Var1, y=Var2, fill=value) ) +
    #    geom_tile( linetype=0, colour=NA ) +
    #    labs(y="Profile End", x="Profile Start", value="OLS-R^2", title=paste(plotTitle, " (OLS.R^2)")) +
    #    coord_fixed() +
    #    geom_point(data=plotdf[plotdf$value<=-2.0,], aes(x=Var1,y=Var2), position=position_nudge(y=4, x=-4), color="white", size=0.6 ) +
    #    theme( plot.margin=unit(c(0.1,0.1,0.1,0.1), "npc"), aspect.ratio=1,
    #                   plot.background = element_blank(),
    #           panel.grid.major = element_blank(),
    #           panel.grid.minor = element_blank(),
    #           panel.background = element_blank(),
    #           axis.line = element_line(color="black")
    #    )
    #print( p, vp=viewport(angle=-45) )
    valScale <- scale_fill_gradient(low="#000000", high="#7799ff", limits=c(0, max(subranges$OLS.R2, 1.0) ) )
    plotRotatedPyramid( plotdf, dotsDF, valScale, Xtrait, plotTitle=paste(plotTitle, " (OLS.R2)") )

}

performGLSregression( traits, tree, "GenomicGC",   "Profile.15" )
performGLSregression( traits, tree, "GenomicGC",   "Profile.1" )
#print(performGLSregression_profileRangeMean( traits, tree, "GenomicGC", c(15,15) ))
#print(performGLSregression_profileRangeMean( traits, tree, "GenomicGC", c(15,19) ))

#--------------------
#glsRegressionRangeAnalysis( traits, tree, "OptimumTemp", c(1,pyramidLength))
performGLSregression( traits, tree, "OptimumTemp", "Profile.11" )
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

performGLSregression( traits, tree, "GenomicGC",   "Profile.1" )
performGLSregression( traits, tree, "GenomicGC",   "Profile.5" )
performGLSregression( traits, tree, "GenomicGC",   "Profile.8" )
performGLSregression( traits, tree, "GenomicGC",   "Profile.15" )
performGLSregression( traits, tree, "GenomicGC",   "Profile.25" )


performGLSregression( traits, tree, "OptimumTemp", "Profile.15" )
performGLSregression( traits, tree, "OptimumTemp", "Profile.15" )

performGLSregression( traits, tree, "PairedFraction", "Profile.15" )
performGLSregression( traits, tree, "PairedFraction", "Profile.5" )
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

warnings()
