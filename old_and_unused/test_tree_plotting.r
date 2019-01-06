library("ape")

pdf('rlineage.pdf')

tr <- read.tree("./data/nmicrobiol201648-s6.txt")

plot.phylo(tr, show.tip.label=FALSE)
add.scale.bar()

chrono <- chronoMPL(tr);

plot.phylo(chrono, show.tip.label=FALSE)
add.scale.bar()


# Source: http://www2.math.su.se/PATHd8/
# Created using: ~/src/PATHd8/PATHd8 ./data/nmicrobiol201648-s6.txt ./data/nmicrobiol201648-s6.txt.nw.PATHd8.nw
tr2 <- read.tree("./data/nmicrobiol201648-s6.txt.nw.PATHd8.out.d8.nw")

plot.phylo(tr2, show.tip.label=FALSE)
add.scale.bar()


dev.off()
