library("grid")
library("ggplot2")

cairo_pdf("test_grid.out.pdf")


df <- data.frame( a=runif(100), b= runif(100) )

p1 <- ggplot( data=df, aes(x=a, y=b) ) + geom_point();


pushViewport( viewport(x = unit(0.5, "npc"), y=unit(0.1, "npc"), angle=45) )

grid.draw( ggplotGrob(p1) )

#upViewport()

vp <- viewport( x=unit(0.5, "npc"), y=unit(0.5, "npc"), h=unit(1, "npc"), w=unit(1, "npc"), angle=-45 )
pushViewport(vp)

xs <- seq(-5, 5, 0.1)
df2 <- data.frame( a=xs, b=xs*0.2-2 )

p2 <- ggplot( data=df2, aes(x=a, y=b) ) +
   geom_point() +
   theme( plot.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_line(color="black"));


grid.draw( ggplotGrob( p2 ) )

upViewport()

vp2 <- viewport( x=unit(0.5, "native"), y=unit(0.5, "native"), h=unit(0.1, "native"), w=unit(1, "native"), angle=-45 )
pushViewport( vp2 )
#grid.rect()
#grid.xaxis()

xaxis <- xaxisGrob( at=seq(0,1,1/6) )

xaxis <- editGrob(xaxis,
                  gPath("major"),
                  x=unit(c(0,1), "npc")
                  )
xaxis <- editGrob(xaxis,
                  gPath("labels"),
                  label=seq(0,300,50)
                  )
grid.draw(xaxis)

#grid.newpage()
#grid.draw( rbind( ggplotGrob(p1), ggplotGrob(p2), size="last") )
dev.off()



