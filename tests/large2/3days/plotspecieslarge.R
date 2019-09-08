color1='deepskyblue'
color2='blue'
color3='darkorchid'

gttimes=read.table('~/speciestrees/tests/large2/3days/results')
mediantimes=tapply(gttimes$V7,gttimes[,c(1,2)],median)
meantimes=tapply(gttimes$V7,gttimes[,c(1,2)],mean)
#par(mfrow=c(3,1)) 
plot(c(0,1000),c(0.1,max(gttimes[,7])),xlab='number of genes per species',ylab='time(seconds)',log="y",col='white',cex.lab=1.3,cex.axis=1.3)
points(gttimes[which(gttimes[,1]==2),2],gttimes[which(gttimes[,1]==2),7],col=color1)
points(gttimes[which(gttimes[,1]==3),2],gttimes[which(gttimes[,1]==3),7],col=color2)
points(gttimes[which(gttimes[,1]==4),2],gttimes[which(gttimes[,1]==4),7],col=color3)

lines(c(20,50,100,200,500,1000),meantimes[1,4:9],col=color1)
lines(c(20,50,100,200,500,1000),meantimes[2,4:9],col=color2)
lines(c(20,50,100,200,500,1000),meantimes[3,4:9],col=color3)

legend(800,1,legend=c('2 species','3 species','4 species'),col=c(color1,color2,color3),lty=c(1,1),cex=1.3)

