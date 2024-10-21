####################################################################
#### Code accompanying:
#
#   Ramirez-Calero, S. et al. 2024. Recurrent extreme climatic events are driving gorgonian populations 
#   to local extinction: low adaptation and low adaptability to marine heat waves.
#
# Script written by: Jean-Baptiste Ledoux & Aldo Barreiro
#
#This script contains information of how heterosis was explored in the data
#This script needs to be run after 2_PCAs_GLM.R and 3_microsat_analysis.R 4_heterosis.R
#Cleaned data is provided in the repository. cleaning process was conducted to match individuals that
#had a microsatellite present and that obtained an intercept after making sure necrosis data was present for the three years only
#individuals with data for two years only were removed.
#
#########################################################################################################

#load packages:
library(lme4)
library("inbreedR")
library(tidyr)
library("MASS")
library("fitdistrplus")
library("ggplot2")

#set working directory
setwd(~)

#load data
#data matched between individuals that have microsatellites and intercepts
dat2<-read.csv("Pcla_interc_het.csv")
head(dat2)

reg1<-lm(Intercepts+1~MLH,data=dat2)
summary(reg1)
plot(dat2$Intercepts+1,dat2$MLH,ylab="nec-int",xlab="sMLH")
abline(reg1,col="red")
par(mfrow=c(2,2))
plot(reg1)
shapiro.test(residuals(reg1))

require(MASS)
Box<-MASS::boxcox(reg1,lambda=seq(-6,6,0.0001))
Cox<-data.frame(Box$x,Box$y)
Cox2<-Cox[with(Cox, order(-Cox$Box.y)),]
Cox2[1,]
lambda<-Cox2[1,"Box.x"] #lambda=0.396
int_trans<-((dat2$Intercepts+1)^lambda-1)/lambda
reg2<-lm(int_trans~MLH,data=dat2)

summary(reg2)
par(mfrow=c(2,2))
plot(reg2)
shapiro.test(residuals(reg2))
plot(dat2$MLH,int_trans,ylab="Intercepts transformadas",
     xlab="Heterocigosidad")
plot(dat2$MLH,int_trans,ylab="nec-int",xlab="sMLH")
plot(int_trans,dat2$MLH,ylab="sMLH",xlab="nec-int")
abline(reg2,col="red")

#Regresión por Monte-Carlo
slope_null<-c()
for (i in 1:10000) {
	x<-sample(dat2$MLH,71,replace=T)
	y<-sample(int_trans,71,replace=T)
	mod<-lm(y~x)
	slope_null[i]<-mod$coefficients[2]
}
lower<-which(slope_null<=reg2$coefficients[2])
length(lower)/10000 # p = 0.3147
ordered_null<-order(slope_null)
ordered_null2<-slope_null[ordered_null]
ordered_null2[250] #-0.8670218
hist(slope_null)
abline(v=reg2$coefficients[2],lty=2,lwd=2,col="red")
abline(v=-0.8196676,lwd=2,col="red")
legend("topright",c("alpha = 0.025","observed slope"),
       lty=c(1,2),lwd=c(2,2),col=c("red","red"),bty="n")



################# End of the code.
