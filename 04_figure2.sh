library(data.table)
library(plotrix)

# library(Hmisc)
# library(fmsb)


# library(lmtest)

petrolblue <- rgb(62,81,192,maxColorValue=255 )

d1 <- fread('significant_traits_afted_adjustment.csv',data.table=F)

d1$CI_lower <- d1$Beta - 1.96*d1$SE
d1$CI_upper <- d1$Beta + 1.96*d1$SE



d2 <- d1
#d2 <- d1[dim(:1,]


d2$color <- "#429CF6"
d2[d2$adjusted=="Yes",]$color <- "#F69C42"
d2 <- d2[dim(d2)[1]:1,]


pdf('sumner_mr_effectiszes_may25_2022.pdf',7,3.5)

par(mai=c(1,2.5,1,1))


barCenters <- barplot( height= d2$Beta , names.arg= d2$Abbreviation, beside=TRUE ,border=F,horiz=TRUE,las=1,xlim=c(-0.6,0.6),xaxt='n',col="white",space=rep(c(12,2),10),cex.names=.75)


segments(d2$CI_lower, y0= barCenters,  
         d2$CI_upper,y1= barCenters,lwd = 2.5,col=d2$color)
         
points(d2$Beta,barCenters,pch=15,col=d2$color)
       

axis(1,at=c(-0.6,-0.3,0,0.3,0.6),cex=1.25)
abline(v=0)

 dev.off()
 
