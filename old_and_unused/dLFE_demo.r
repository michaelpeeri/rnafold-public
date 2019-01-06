library('ggplot2')
library('reshape')

N <- 200
nativeLFE <- cumsum(rnorm(n=N, mean=0.1))+1

K <- 50

randLFE <- matrix(rnorm(n=N*K, mean=matrix(rep.int(cumsum(rnorm(n=N, mean=0.1)), K), ncol=K)), ncol=K)
rownames(randLFE) <- 1:N

randLFE.mean <- rowMeans(randLFE)

dLFE <- data.frame( nativeLFE=nativeLFE, randLFE.mean=randLFE.mean, X1=1:N )

p <- ggplot(melt(randLFE), aes(x=X1) ) +
    geom_ribbon( data=dLFE, aes(ymin=nativeLFE, ymax=randLFE.mean), alpha=0.5, color=NA, fill="#efcf60" ) +
    geom_point( aes(y=value), alpha=0.2, size=0.5, color="purple") +
    geom_line( data=melt(matrix(nativeLFE)), aes(y=value), color="#209f40", size=1.5 ) +
    guides( color=guide_legend(), fill=guide_legend() )
    
    
print(p)

dev.off()
