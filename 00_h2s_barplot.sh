R
library(data.table)
library(plotrix)
library(scales)

#this needs to be updated for what is in hte mixer plots, i think bipoalr may be old!

res2 <- fread('heritabilities_for_plot.csv',data.table=F)


res2$LCI <- res2$h2- 1.96*res2$se
res2$UCI <- res2$h2 + 1.96*res2$se

res2$color <- NA
caro_orange <- rgb(249,164,24,maxColorValue=255)
beauty_red <- rgb(237,33,36,maxColorValue=255)
adam_blue <- rgb(41,153,208,maxColorValue=255)


res2$color <- "black"

res2[which(res2$Category == "Biomarker"),]$color <- rgb(255,1,0,max=255)
res2[which(res2$Category == "Connective Tissues"),]$color <- rgb(255,116,0,max=255)
res2[which(res2$Category == "Digestive System"),]$color <- rgb(255,171,0,max=255)
res2[which(res2$Category == "Endocrine"),]$color <- rgb(255,213,2,max=255)
res2[which(res2$Category == "Hematological"),]$color <-rgb(255,255,2,max=255)
res2[which(res2$Category == "Inflammatory Arthritis"),]$color <- rgb(160,239,3,max=255)
res2[which(res2$Category == "Nervous System"),]$color <- rgb(3,206,3,max=255)
res2[which(res2$Category == "PTSD"),]$color <- rgb(0,153,154,max=255)
res2[which(res2$Category == "Skin"),]$color <- rgb(18,63,172,max=255)
res2[which(res2$Category == "Vasculitis"),]$color <- rgb(57,20,175,max=255)





#this will be fore the second data


cis2 <- rbind(res2)

pdf('ads_h2s.pdf',7,8)
par(mar=c(15, 5, 5, 5) + 0.5)
plots <- barplot(height=cis2$h2,names.arg=cis2$Phenotype,col=cis2$color,las=2,srt=75,horiz=FALSE,ylab="h2SNP",xaxt='n',cex.axis=1.2,ylim=c(0,0.3))
#abline(h=-log10(0.05/40),col='black',lwd=2,lty=2)

arrows(plots, cis2$h2+2*cis2$se,plots, cis2$h2-2*cis2$se, angle=90, code=3,lwd=2, length=0.05, col=rgb(0,0,0,0.7))
text(plots,rep(-0.01,nrow(cis2)),labels=cis2$Phenotype,srt=70,las=2,xpd=TRUE,adj=1)

dev.off()

polygen <- fread('polygenicities_for_plot.csv',data.table=F)

polygen$color <- "black"
polygen[which(polygen$Category == "Biomarkers"),]$color <- rgb(255,1,0,max=255)
polygen[which(polygen$Category == "Connective Tissues"),]$color <- rgb(255,116,0,max=255)
polygen[which(polygen$Category == "Digestive System"),]$color <- rgb(255,171,0,max=255)
polygen[which(polygen$Category == "Endocrine"),]$color <- rgb(255,213,2,max=255)
polygen[which(polygen$Category == "Hematological"),]$color <-rgb(255,255,2,max=255)
polygen[which(polygen$Category == "Inflammatory Arthritis"),]$color <- rgb(160,239,3,max=255)
polygen[which(polygen$Category == "Nervous System"),]$color <- rgb(3,206,3,max=255)
polygen[which(polygen$Category == "PTSD"),]$color <- rgb(0,153,154,max=255)
polygen[which(polygen$Category == "Skin"),]$color <- rgb(18,63,172,max=255)
polygen[which(polygen$Category == "Vasculitis"),]$color <- rgb(57,20,175,max=255)



pdf('ads_genarch.pdf',7,7)
par(mar=c(5, 5, 5, 5) + 0.5)
plot(log10(polygen$ncausal+1),polygen$discoverability,col="white",xlab="log number of influential variants", ylab="Effect size variance (discoverability)" ,cex.lab=1.45,cex.axis=1.24)
points(log10(polygen$ncausal+1),polygen$discoverability,col=polygen$color,pch=19,cex=polygen$size*20)
#text(log10(polygen$ncausal+1),polygen$discoverability,polygen$Phenotype,cex=1.3)

dev.off()
