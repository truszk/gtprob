color1='darkorange'
color2='firebrick1'
color3='brown'

gttimes=read.table('~/speciestrees/tests/large2/3days/results')
mediantimes=tapply(gttimes$V7,gttimes[,c(1,2)],median)
meantimes=tapply(gttimes$V7,gttimes[,c(1,2)],mean)
#par(mfrow=c(3,1)) 
plot(c(0,1024),c(0.1,max(gttimes[,7])),xlab='number of species',ylab='time(seconds)',log="y",col='white',cex.lab=1.3,cex.axis=1.3)
points(gttimes[which(gttimes[,2]==1),1],gttimes[which(gttimes[,2]==1),7],col=color1)
points(gttimes[which(gttimes[,2]==3),1],gttimes[which(gttimes[,2]==3),7],col=color2)
points(gttimes[which(gttimes[,2]==5),1],gttimes[which(gttimes[,2]==5),7],col=color3)

lines(c(32,64,128,256,512,1024),meantimes[4:9,1],col=color1)
lines(c(32,64,128,256,512,1024),meantimes[4:9,2],col=color2)
lines(c(32,64,128,256,512,1024),meantimes[4:9,3],col=color3)

legend(800,1,legend=c('1 gene/sp.','3 genes/sp.','5 genes/sp.'),col=c(color1,color2,color3),lty=c(1,1),cex=1.3)


