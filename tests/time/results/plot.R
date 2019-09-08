gttimes=read.table('~/speciestrees/tests/time/results/results')
par(mfrow=c(3,1)) 
plot(1,type='n',xlim=c(0,42),ylim=c(0.001,max(gttimes[,7])),xlab='Number of species, 1 gene/species',ylab='time(seconds)',log="y")
points(gttimes[which(gttimes[,2]==1&gttimes[,4]==1),1],gttimes[which(gttimes[,2]==1&gttimes[,4]==1),7],col='blue')
points(gttimes[which(gttimes[,2]==1&gttimes[,4]==0),1],gttimes[which(gttimes[,2]==1&gttimes[,4]==0),7],col='red')

plot(1,type='n',xlim=c(0,42),ylim=c(0.001,max(gttimes[,7])),xlab='Number of species, 3 genes/species',ylab='time(seconds)',log="y")
points(gttimes[which(gttimes[,2]==3&gttimes[,4]==1),1],gttimes[which(gttimes[,2]==3&gttimes[,4]==1),7],col='blue')
points(gttimes[which(gttimes[,2]==3&gttimes[,4]==0),1],gttimes[which(gttimes[,2]==3&gttimes[,4]==0),7],col='red')

plot(1,type='n',xlim=c(0,42),ylim=c(0.001,max(gttimes[,7])),xlab='Number of species, 5 gene/species',ylab='time(seconds)',log="y")
points(gttimes[which(gttimes[,2]==5&gttimes[,4]==1),1],gttimes[which(gttimes[,2]==5&gttimes[,4]==1),7],col='blue')
points(gttimes[which(gttimes[,2]==5&gttimes[,4]==0),1],gttimes[which(gttimes[,2]==5&gttimes[,4]==0),7],col='red')

