library("ggplot2")


N <- 20

pdf("test_color_range.out.pdf")
df <- data.frame( x=runif(N), y=rnorm(N), c=as.factor(sample(5, size=N, replace=TRUE)) )
print(df)

p <- ggplot( df, aes(x=x, y=y) ) +
    geom_point( aes(colour=c) )
print(p)

warnings()
dev.off()
