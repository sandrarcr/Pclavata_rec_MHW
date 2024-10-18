####################################################################
#### Code accompanying:
#
#   Ramirez-Calero, S. et al. 2024. Recurrent extreme climatic events are driving gorgonian populations 
#   to local extinction: low adaptation and low adaptability to marine heat waves.
#
# Script written by: Aldo Barreiro
#
#This script contains information on null model testing. Associated data can be found in the repository
#########################################################################################################

#load data

dat<-read.csv("data_sensitivity_analysis.csv",stringsAsFactors=T)
head(dat)
str(dat)

plot(dat$Pop,dat$PCScores)

###Significancia estadística stándar
sp.dat<-split(dat,dat$Year)
pvals<-NULL
slope<-NULL
for (i in 1:76){
		mod<-lm(c(sp.dat[[1]]$PCScores[i],sp.dat[[2]]$PCScores[i],sp.dat[[3]]$PCScores[i])~c(sp.dat[[1]]$Pop[i],sp.dat[[2]]$Pop[i],sp.dat[[3]]$Pop[i]))
        pvals[i]<-summary(mod)$coefficients[2,4]
        slope[i]<-summary(mod)$coefficients[2,1]
}
length(which(pvals<=0.05))/length(pvals) #40.7% of p-values are significant
hist(slope)

###Comparación con el modelo nulo
coef.n<-list()
set.seed(123)
for (i in 1:10000){
	    mod<-lm(sample(dat$PCScores,3,replace=T)~sample(dat$Pop,3,replace=T))
        coef.n[[i]]<-summary(mod)$coefficients
}

slope.n<-NULL
for (i in 1:length(coef.n)){
	slope.n[i]<-ifelse(nrow(coef.n[[i]])==2,coef.n[[i]][2,1],0)
}

#Comparar las pendientes reales com las del modelo nulo
d.n<-density(slope.n)
d<-density(slope)
mean(slope.n) #-0.06
mean(slope) #1
sd(slope.n) #5.1
sd(slope) #0.37

length(which(slope>=mean(slope.n)))/length(slope) #98.7% of the slopes are greater than the average null model slope

pdf("Fig.null_slopes.pdf")
plot(d.n$x,d.n$y,type="l",lwd=2,col="blue",ylim=c(0,1.2),xlab="Slope",ylab="Density")
lines(d$x,d$y,lwd=2,col="red")
legend("topright",c("null model","real data"),lty=c(1,1),lwd=c(2,2),col=c("blue","red"),bty="n")
dev.off()


