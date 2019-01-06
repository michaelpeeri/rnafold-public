library("ggplot2")
library("reshape")


data <- data.frame(matrix(runif(24*4), nrow=24, ncol=4))
data$group <- rownames(data)

data2 <- melt(data)
data2$pos <- rownames(data2)

p <- ggplot( melt(data2) ) +
    geom_col( aes( x=pos, y=value, fill=variable ) )
print(p)

dev.off()
quit()
    
