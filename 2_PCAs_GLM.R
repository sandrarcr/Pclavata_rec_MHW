####################################################################
#### Code accompanying:
#
#   Ramirez-Calero, S. et al. 2024. Recurrent extreme climatic events are driving gorgonian populations 
#   to local extinction: low adaptation and low adaptability to marine heat waves.
#
# Script written by: Aldo Barreiro
#
#This script contains information of how exploratory PCAs were done as well a the linear model tested
#in the study. Associated data can be found in the repository
#########################################################################################################

#load data
dat<-read.csv2("datos_pclavata.csv",stringsAsFactors=T)
head(dat)
names(dat)[1]<-"Colony_ID"
names(dat)[7]<-"Necrosis"
str(dat)

#data cleaning - individuals with missing data
ind_pl<-which(dat$Colony_ID=="P_Llop_18"|dat$Colony_ID=="P_Llop_23"|dat$Colony_ID=="P_Llop_NI_1_2016"|dat$Colony_ID=="P_Llop_NI_1_2017"|dat$Colony_ID=="P_Llop_NI_2_2016"|dat$Colony_ID=="P_Llop_NI_2_2017"|dat$Colony_ID=="P_Llop_NI_3_2016"|dat$Colony_ID=="P_Llop_NI_3_2017"|dat$Colony_ID=="P_Llop_NI_4_2016"|dat$Colony_ID=="P_Llop_NI_4_2017")
ind_t<-which(dat$Colony_ID=="Tascons_8"|dat$Colony_ID=="Tascons_31")
ind_v<-which(dat$Colony_ID=="Vaca_14"|dat$Colony_ID=="Vaca_19"|dat$Colony_ID=="Vaca_19"|dat$Colony_ID=="Vaca_20"|dat$Colony_ID=="Vaca_27"|dat$Colony_ID=="Vaca_28"|dat$Colony_ID=="Vaca_29"|dat$Colony_ID=="Vaca_4"|dat$Colony_ID=="Vaca_30"|dat$Colony_ID=="Vaca_31"|dat$Colony_ID=="Vaca_32"|dat$Colony_ID=="Vaca_33"|dat$Colony_ID=="Vaca_34"|dat$Colony_ID=="Vaca_35"|dat$Colony_ID=="Vaca_36"|dat$Colony_ID=="Vaca_37"|dat$Colony_ID=="Vaca_38"|dat$Colony_ID=="Vaca_39"|dat$Colony_ID=="Vaca_40"|dat$Colony_ID=="Vaca_41"|dat$Colony_ID=="Vaca_42"|dat$Colony_ID=="Vaca_6"|dat$Colony_ID=="Vaca_7"|dat$Colony_ID=="Vaca_9"|dat$Colony_ID=="Vaca_NI_3_2016"|dat$Colony_ID=="Vaca_NI_3_2017"|dat$Colony_ID=="Vaca_NI_4_2016"|dat$Colony_ID=="Vaca_NI_4_2017"|dat$Colony_ID=="Vaca_NI_5_2017")
inds<-c(ind_pl,ind_t,ind_v)
dat2<-dat[-inds,]

#adapt population names
v<-which(dat3$Population=="La Vaca")
vc<-rep("Vaca",times=length(v))
dat3$Population[v]<-vc

#remove unexistent factors
dat3$Colony_ID<-droplevels(dat3$Colony_ID)
dat3$Population<-droplevels(dat3$Population)

#Add 'Temperature'
newtemp1<-aggregate(Temp~Colony_ID+Year+Conditions,data=dat3,mean,na.rm=T)
combfact<-as.factor(paste(dat3$Colony_ID,dat3$Year,dat3$Conditions))
list1<-split(combfact,f=combfact)
newtemp2<-list(length=length(list1))
for (i in 1:length(list1)) {
  newtemp2[[i]]<-rep(newtemp1$Temp[i],length=length(list1[[i]]))
}
newtemp3<-unlist(newtemp2)
dat4<-data.frame(dat3[,-5],mean.temp=newtemp3)
head(dat4)
str(dat4)

#Logistic model applied to individual colonies
dat4$newlab<-as.factor(paste(dat4$Colony_ID,dat4$Year,dat4$Conditions))
#pars1<-data.frame(k=NA,d0=NA)
#for (i in 4:length(levels(dat4$newlab))){
#	dat<-dat4[dat4$newlab==levels(dat4$newlab)[i],]
#	mod<-nls(Necrosis~(100/(1+exp(-k*(Experiment.Day-#d0)))),start=list(k=0.4,d0=15),data=dat,nls.control(maxiter=5000,minFactor=1/100000))
#	pars1[i,]<-summary(mod)$parameters[,1]
#}
#pars2<-data.frame(k=NA,d0=NA)
#for (i in 6:length(levels(dat4$newlab))){
#	dat<-dat4[dat4$newlab==levels(dat4$newlab)[i],]
#	mod<-nls(Necrosis~(100/(1+exp(-k*(Experiment.Day-#d0)))),start=list(k=0.4,d0=15),data=dat,nls.control(maxiter=5000,minFactor=1/100000))
#	pars2[i,]<-summary(mod)$parameters[,1]
#}

#Individual cases in which the model cannot be used as they have final values of 0% o very low values
#of 30%, 20%
dat4[dat4$newlab==levels(dat4$newlab)[6],]
#Logistic curve cannot be applied. Too many individuals in which no parameters can be estimated
#so all data from 'Control' is ignored

#sort data per days within each category
dat5<-dat4[with(dat4,order(newlab,Experiment.Day)),]
head(dat5)
str(dat5)

#Select first day when necrosis for individuals appear. This was not used but is left as reference.
#dat5.l<-split(dat5,dat5$newlab)
#day.n<-list()
#for (i in 1:length(dat5.l)){
#	day.n[[i]]<-dat5.l[[i]]$Experiment.Day[which(dat5.l[[i]]$Necrosis>0)]
#}
#first.day<-NULL
#for (i in 1:length(day.n)){
#	first.day[i]<-min(day.n[[i]])
#}
#first.day2<-NULL
#for (i in 1:length(first.day)){
#	first.day2[i]<-ifelse(first.day[i]==Inf,29,first.day[i])
#}
#dat1st<-data.frame(ind=unique(dat5$newlab),day=first.day2)

#Percentages of necrosis in columns per day and matching among days.
dat5.l<-split(dat5,dat5$newlab)
nec<-list()
for (i in 1:length(dat5.l)){
  nec[[i]]<-dat5.l[[i]]$Necrosis
}
days<-list()
for (i in 1:length(dat5.l)){
  days[[i]]<-dat5.l[[i]]$Experiment.Day
}
test<-0:28
days.new<-list()
nec.new<-list()
for ( i in 1:length(days)){
  days.test<-days[[i]]
  nec.test<-nec[[i]]
  days.new[[i]]<-days.test[match(test,days.test)]
  nec.new[[i]]<-nec.test[match(test,days.test)]
}

#remove NAs
nec.new2<-nec.new
for (i in 1:length(nec.new)){
  for (j in 1:29) {
    if (is.na(nec.new[[i]][j])){
      nec.new2[[i]][j]<-nec.new[[i]][j-1]
    } else (nec.new2[[i]][j]<-nec.new[[i]][j])
  }
}
nec.new3<-nec.new2
for (i in 1:length(nec.new2)){
  for (j in 1:29) {
    if (is.na(nec.new2[[i]][j])){
      nec.new3[[i]][j]<-nec.new2[[i]][j-1]
    } else (nec.new3[[i]][j]<-nec.new2[[i]][j])
  }
}
nec.new4<-nec.new3
for (i in 1:length(nec.new3)){
  for (j in 1:29) {
    if (is.na(nec.new3[[i]][j])){
      nec.new4[[i]][j]<-nec.new3[[i]][j-1]
    } else (nec.new4[[i]][j]<-nec.new3[[i]][j])
  }
}

which(is.na(nec.new4)) #ok

#Remove NAs
days.new2<-days.new
for (i in 1:length(days.new)){
  for (j in 1:29) {
    if (is.na(days.new[[i]][j])){
      days.new2[[i]][j]<-test[j]
    } else (days.new2[[i]][j]<-days.new[[i]][j])
  }
}

which(is.na(days.new2))

#New data frame
condition<-unique(dat5$newlab)
nec<-t(data.frame(nec.new4))
rownames(nec)<-condition
colnames(nec)<-test
head(nec)
str(nec)

#NAs
which(is.na(nec))
# [1]  5295  5733  6171  6609  7047  7485  7923  8361  8799  9237  9675 10113 10551 10989 11427 11865 12303
#Remove outlier individual P_Llop_2 2015 Control
nec2<-as.data.frame(nec)
nec2$inds<-c(rownames(nec))
rownames(nec2)<-1:438
head(nec2)
nec2[nec2$inds=="P_Llop_2 2015 Control",] #row 39

nec<-nec[-39,]
which(is.na(nec))

#correcting necrosis values
for (j in 2:ncol(nec)){
  for (i in 1:nrow(nec)){
    if (nec[i,j]<nec[i,j-1]){nec[i,j]<-nec[i,j-1]}
    else{nec[i,j]<-nec[i,j]}
  }
}

#write new data frame
write.csv(nec,"necrosis_cor.csv") #file can be found in the repository

#load data
nec<-read.csv("necrosis_cor.csv")
head(nec)
str(nec)

require(FactoMineR)

#divide text to obtain factors: Individual, year y condition
factors<-matrix(NA,ncol=3,nrow=437)
for (i in 1:437) {
	factors[i,]<-strsplit(as.character(nec$X[i])," +")[[1]]
}
#Population
pop<-matrix(NA,ncol=1,nrow=437)
for (i in 1:437) {
	pop[i]<-substr(factors[i,1],1,4)
}

#Data frame removing controls
nec2<-data.frame(nec,pop=as.factor(pop),year=as.factor(factors[,2]),cond=as.factor(factors[,3]),rep=nec$X)
head(nec2)
nec3<-nec2[nec2$cond=="Tractament",]
head(nec3)

nrow(nec3[nec3$year==2015,]) #67
nrow(nec3[nec3$year==2016,]) #78
nrow(nec3[nec3$year==2017,]) #76

#PCA 1
res.pca.tot<-PCA(nec3[,-c(1,31:34)])
sc<-res.pca.tot$var$cor[,1:2]

pdf("pca1.pdf")
plot(sc,type="n",xlab="PCA1 (50.29% variance)",ylab="PCA2 (27.36% variance)")
arrows(x0=0,y0=0,x1=sc[,1],y1=sc[,2],length=0.1)
text(x=sc[,1],y=sc[,2],labels=1:28,pos=4)
dev.off()

#PCA 2 removing days 1-9
res.pca<-PCA(nec3[,-c(1:11,31:34)]) #1st component = 74.81% of variance
sc2<-res.pca$var$cor
pca1<-res.pca$ind$coord[,1]
data.pca<-data.frame(nec3,pca1=pca1)
head(data.pca)

pdf("pca2.pdf")
plot(sc2,type="n",xlab="PCA1 (74.81% variance)",ylab="PCA2 (14.93% variance)",xlim=c(0,1),ylim=c(-0.4,max(sc2[,2])))
arrows(x0=0,y0=0,x1=sc2[,1],y1=sc2[,2],length=0.1)
text(x=sc2[,1],y=sc2[,2],labels=10:28,pos=4)
dev.off()

#add factor individual
inds<-matrix(NA,ncol=3,nrow=221)
for (i in 1:221) {
	inds[i,]<-strsplit(as.character(nec3$rep[i])," +")[[1]]
}
data.pca$ind<-as.factor(inds[,1])
#remove individuals with only one data
data.pca2<-data.pca[-(219:221),]
data.pca2$ind<-droplevels(data.pca2$ind)

#data exploration
boxplot(data.pca2$pca1~data.pca2$pop+data.pca2$year)
boxplot(data.pca2$pca1~data.pca2$pop)
boxplot(data.pca2$pca1~data.pca2$year)
boxplot(data.pca2$pca1~data.pca2$ind)

#Model testing
min(data.pca2$pca1)
mod<-lm(pca1+4~pop*year,data=data.pca2)
summary(mod)
require(car)
Anova(mod) 
mod2<-lm(pca1+4~pop+year,data=data.pca2)
par(mfrow=c(2,2))
plot(mod2)
ks.test(residuals(mod2),"pnorm",mean=mean(residuals(mod2)),sd=sd(residuals(mod2)))
shapiro.test(residuals(mod))
require(MASS)
Box<-MASS::boxcox(mod2,lambda=seq(-6,6,0.0001))
Cox<-data.frame(Box$x,Box$y)
Cox2<-Cox[with(Cox, order(-Cox$Box.y)),]
Cox2[1,]
lambda<-Cox2[1,"Box.x"] #lambda=0.449
vartrans<-((data.pca2$pca1+4)^lambda-1)/lambda
mod3<-lm(vartrans~pop+year,data=data.pca2)
summary(mod3)
Anova(mod3)
par(mfrow=c(2,2))
plot(mod3)
outlierTest(mod3) #No outliers
ks.test(residuals(mod3),"pnorm",mean=mean(residuals(mod3)),sd=sd(residuals(mod3)))

#model using factor individual as random factor
require(lme4)
mod3.1<-lmer(vartrans~pop*year+(1|ind),data=data.pca2) 
AIC(mod3,mod3.1) #mod3.1 best AIC and significant
anova(mod3.1,mod3) 

require(multcomp)
sink("mod_3.1.txt")
Anova(mod3.1)
summary(glht(mod3.1,linfct=mcp(year="Tukey")))
sink()

#model using individuals with at least three replicates
data.pca3<-data.pca2[-(1:20),]
data.pca3$ind<-droplevels(data.pca3$ind)

mod4<-lm(pca1+4~pop*year,data=data.pca3)
summary(mod4)
Anova(mod4)
mod5<-lm(pca1+4~pop+year,data=data.pca3)
par(mfrow=c(2,2))
plot(mod5)
ks.test(residuals(mod5),"pnorm",mean=mean(residuals(mod5)),sd=sd(residuals(mod5)))
shapiro.test(residuals(mod5))
Box<-MASS::boxcox(mod5,lambda=seq(-6,6,0.0001))
Cox<-data.frame(Box$x,Box$y)
Cox2<-Cox[with(Cox, order(-Cox$Box.y)),]
Cox2[1,]
lambda<-Cox2[1,"Box.x"] #lambda=0.4641
vartrans2<-((data.pca3$pca1+4)^lambda-1)/lambda
mod6<-lm(vartrans2~pop+year,data=data.pca3)
summary(mod6)
Anova(mod6)
par(mfrow=c(2,2))
plot(mod6)
outlierTest(mod6) #No outliers
ks.test(residuals(mod6),"pnorm",mean=mean(residuals(mod6)),sd=sd(residuals(mod6)))
mod6.1<-lmer(vartrans2~pop+year+(1|ind),data=data.pca3)  
AIC(mod6,mod6.1) #mod6 best AIC

#Visualization model for fixed factors
require(visreg)
par(mfrow=c(1,2))
visreg(mod3.1)
par(mfrow=c(1,2))
visreg(mod6.1)

#Visualization model for random factors
require(sjPlot)
plot_model(mod3.1, "re")
plot_model(mod6.1, "re")

#extract intercepts for each individual taken as their phenotypic differential reponses. 
intercepts1<-coef(mod3.1)$ind[,1]
ind_intercepts1<-data.frame(ind=levels(data.pca$ind),intercepts=intercepts1)
intercepts2<-coef(mod6.1)$ind[,1]
ind_intercepts2<-data.frame(ind=levels(data.pca2$ind),intercepts=intercepts2)

#save resultant file
write.csv(ind_intercepts2,"Intercepts_ind_3years.csv") #can be found on repository

################# End of the code.