library("ggplot2")
library("reshape")


pdf("test_raster.out.pdf")


df <- data.frame( a=c(10,15,20,15), b=c(15,20,25,20), c=c(20,25,35,25), d=c(15,20,25,20) )
colnames(df) <- c(1,2,3,4)
rownames(df) <- c(1,2,3,4)
df$x = 1:4
dfm <- melt(df, id.vars=c("x"))

dg <- data.frame( a=c(0,0,1,0), b=c(0,0,1,1), c=c(0,1,1,1), d=c(0,1,0,1) )
colnames(dg) <- c(1,2,3,4)
rownames(dg) <- c(1,2,3,4)
dg$x = 1:4
dgm <- melt(dg, id.vars=c("x"))

p <- ggplot( data=dfm, aes(x=x, y=variable, fill=value) ) +
    geom_raster() +
    geom_point(data=dgm, aes(x=x,y=variable, alpha=value ), position=position_nudge(y=0.45, x=0.45), color="white" )
print(p)

dev.off()
