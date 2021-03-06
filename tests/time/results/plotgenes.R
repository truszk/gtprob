gttimes=read.table('~/speciestrees/tests/time/results/results')
mediantimes=tapply(gttimes$V7,gttimes[,c(1,2,4)],median)
meantimes=tapply(gttimes$V7,gttimes[,c(1,2,4)],mean)
par(mfrow=c(3,1)) 
plot(1,type='n',xlim=c(0,50),ylim=c(0.001,max(gttimes[,7])),xlab='samples per species, 2 leaf species',ylab='time(seconds)',log="y")
points(gttimes[which(gttimes[,1]==2&gttimes[,4]==1),2],gttimes[which(gttimes[,1]==2&gttimes[,4]==1),7],col='blue')
points(gttimes[which(gttimes[,1]==2&gttimes[,4]==0),2],gttimes[which(gttimes[,1]==2&gttimes[,4]==0),7],col='red')
points(gttimes[which(gttimes[,1]==2&gttimes[,4]==2),2],gttimes[which(gttimes[,1]==2&gttimes[,4]==2),7],col='green')
lines(c(5,10,15,20,50),meantimes[1,3:7,2],col='blue')
lines(c(5,10,15,20,50),meantimes[1,3:7,1],col='red')
lines(c(5,10,15,20,50),meantimes[1,3:7,3],col='green')

plot(1,type='n',xlim=c(0,50),ylim=c(0.001,max(gttimes[,7])),xlab='samples per species, 3 leaf species',ylab='time(seconds)',log="y")
points(gttimes[which(gttimes[,1]==3&gttimes[,4]==1),2],gttimes[which(gttimes[,1]==3&gttimes[,4]==1),7],col='blue')
points(gttimes[which(gttimes[,1]==3&gttimes[,4]==0),2],gttimes[which(gttimes[,1]==3&gttimes[,4]==0),7],col='red')
points(gttimes[which(gttimes[,1]==3&gttimes[,4]==2),2],gttimes[which(gttimes[,1]==3&gttimes[,4]==2),7],col='green')
lines(c(5,10,15,20,50),meantimes[2,3:7,2],col='blue')
lines(c(5,10,15,20,50),meantimes[2,3:7,1],col='red')
lines(c(5,10,15,20,50),meantimes[2,3:7,3],col='green')

plot(1,type='n',xlim=c(0,50),ylim=c(0.001,max(gttimes[,7])),xlab='samples per species, 4 leaf species',ylab='time(seconds)',log="y")
points(gttimes[which(gttimes[,1]==4&gttimes[,4]==1),2],gttimes[which(gttimes[,1]==4&gttimes[,4]==1),7],col='blue')
points(gttimes[which(gttimes[,1]==4&gttimes[,4]==0),2],gttimes[which(gttimes[,1]==4&gttimes[,4]==0),7],col='red')
points(gttimes[which(gttimes[,1]==4&gttimes[,4]==2),2],gttimes[which(gttimes[,1]==4&gttimes[,4]==2),7],col='green')
lines(c(2,3,5,10,15,20,50),meantimes[3,1:7,2],col='blue')
lines(c(2,3,5,10,15,20,50),meantimes[3,1:7,1],col='red')
lines(c(2,3,5,10,15,20,50),meantimes[3,1:7,3],col='green')

