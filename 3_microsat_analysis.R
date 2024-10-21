####################################################################
#### Code accompanying:
#
#   Ramirez-Calero, S. et al. 2024. Recurrent extreme climatic events are driving gorgonian populations 
#   to local extinction: low adaptation and low adaptability to marine heat waves.
#
# Script was written by: Sandra Ramirez
#
# This script includes:
# Data curation and cleaning 
# Basic stats (Allelic richness, Ho, He, Hst
# HWE
# LD (pairwise)
# Inbreeding coefficients: Fst, Fis
# PCA 
# DAPC
# AMOVA
#########################################################################################################

#setwd
setwd("~")

#load libraries
library("adegenet") 
library("pegas") 
library("mmod")
library("reshape2")
library("ggplot2")
library("poppr")
library("devtools")
library("hierfstat")
library("genepop")
library("graph4lg")
library("lattice")
library("magrittr")
library("cowplot")

#load data
P.clav <- read.genalex("P.clavata_microsat.csv") #can be found on repository and corresponds to raw microsatellite data
Pcl <- as.genclone(P.clav) 
Pclavata <- genclone2genind(Pcl)

Pcla_genepop<-genind_to_genepop(P.clav, output = "data.frame") # convert to genepop format
#write.table(Pcla_genepop,"Pcla_genepop.txt",sep="/t",row.names=FALSE) 

#data exploration
#population size
table(pop(Pclavata)) 

#MPL MTA MVA 
#29  33  25

#genind2df(Pclavata, sep = "|")
sum <- summary (Pclavata)
names(sum)

par(mfrow=c(2,2))
plot(sum$n.by.pop, sum$pop.n.all, xlab="Colonies sample size",
     ylab="Number of alleles",main="Alleles numbers and sample sizes",
     type="n")
text(sum$n.by.pop,sum$pop.n.all,lab=names(sum$n.by.pop))

barplot(sum$loc.n.all, ylab="Number of alleles",
        main="Number of alleles per locus")
barplot(sum$Hexp-sum$Hobs, main="Heterozygosity: expected-observed",
        ylab="Hexp - Hobs")
barplot(sum$n.by.pop, main="Sample sizes per population",
        ylab="Number of genotypes",las=3)
#extra
barplot(Pclavata$loc.n.all,col=funky(17), las=3, ylim = c(0,25), cex.names = 0.80,
        xlab="Locus", ylab="Number of allels", besides=TRUE) #nuber of alleles per locus

#check and eliminate missing data:

info_table(Pclavata, plot = TRUE)
pinflt <- locus_table(Pclavata, information = T)
MPota <- locus_table(Pclavata, pop = "MPL")
MVaca <- locus_table(Pclavata, pop = "MVA")
MTas <- locus_table(Pclavata, pop = "MTA")

Pclavata_nomiss <- missingno(Pclavata, type = "loci", cutoff = 0.1, quiet = FALSE, freq = FALSE)
Pclavata_nomiss
shufflepop(Pclavata_nomiss, method = 1)
#write.csv(Pclavata_nomiss, file = "P.clavata_clean.csv")
names(Pclavata_nomiss)

######### Basic population genetics analyses #########

####### Allele richness #######
all.rich <- allelic.richness(Pclavata_nomiss,min.n=174, min.all, diploid=TRUE)
#write.csv(all.rich, file = "allelic-richness.csv")
all.rich.pegas <- allelicrichness(as.loci(Pclavata_nomiss))
#write.csv(all.rich.pegas, file = "allelic-rich per pop.csv")

####### #Basic stats #####

bastats<- basic.stats(Pclavata_nomiss,diploid=TRUE)

pclavatadf <- genind2hierfstat(Pclavata_nomiss)#changed format to use hierfstat package

# Expected and observed heterozygosity

#plot Hexp & Ho
temp <- summary(Pclavata_nomiss)

plot(temp$Hexp, temp$Hobs, pch=20, cex=3, xlim=c(.4,1), ylim=c(.4,1), xlab="Hexp", 
     ylab="Hob")
abline(0,1,lty=2)

bartlett.test(list(temp$Hexp,temp$Hobs))

t.test(temp$Hexp,temp$Hobs, pair=T,var.equal = TRUE,alter="greater")

#Deficit in heterozygosity can be indicative of population structure. In the following, we try
#to assess this possibility using classical population genetics tools

# Hardy-Weinberg equilibrium 

HW_pclavata <- hw.test(Pclavata_nomiss, B=1000) #HWE for each locus
Pclavatahwe.pop <- seppop(Pclavata_nomiss) %>% lapply(hw.test, B = 0) #HWE per population

#Linkage dissequilibrium 

tascons <- popsub(Pclavata_nomiss, "MTA")
tas_bas <- basic.stats(tascons)
ia(tascons, sample = 999)

pota <- popsub(Pclavata_nomiss, "MPL")
ia(pota, sample = 999)

vaca <- popsub(Pclavata_nomiss, "MVA")
ia(vaca, sample = 999)

#Pairwise LD over all loci

tascons_pair <- pair.ia(Pclavata_nomiss[pop="MTA"])
plot_range_tas <- range(c(tascons_pair), na.rm = TRUE)
tas_LD <- plot(tascons_pair, limits = plot_range_tas, 
     low = "lightblue", high = "red", index = "rbarD", 
     index.col = "transparent", pch = NA)
write.csv(tascons_pair, file="tascons_LD_pairwise_values.csv")

pota_pair <- pair.ia(Pclavata_nomiss[pop="MPL"], plot=TRUE)
plot_range_pota <- range(c(pota_pair), na.rm = TRUE)
pota_LD <- plot(pota_pair, limits = plot_range_pota, 
     low = "lightblue", high = "red", index = "rbarD", 
     index.col = "transparent", pch = NA)
write.csv(pota_pair, file="pota_LD_pairwise_values.csv")

vaca_pair <- pair.ia(Pclavata_nomiss[pop="MVA"], plot=TRUE)
plot_range_vaca <- range(c(vaca_pair), na.rm = TRUE)
vaca_LD <- plot(vaca_pair, limits = plot_range_vaca, 
     low = "lightblue", high = "red", index = "rbarD", 
     index.col = "transparent", pch = NA)
write.csv(vaca_pair, file="vaca_LD_pairwise_values.csv")

plot_grid(tas_LD, pota_LD, vaca_LD)

#check final values:

head(tascons_pair)
head(pota_pair)
head(vaca_pair)

#Overall F statistics

#global fst
wc_fst <- wc(Pclavata_nomiss, diploid = TRUE,) #Weir and Cockerham estimates of Fstatistics
#write.csv(wc_fst[["per.loc"]], file = "wc_fst_per_locus.csv")

#pairwise fst
pairwise_pc <- pairwise.WCfst(pclavatadf, diploid = T)
#write.csv(pairwise_pc, file="pairwise_mat_fst.csv")

#Are these values significant?

pclavatadf #hierfstat df

attach(pclavatadf)
loci <- data.frame(Pcla.21,Pcla.22,Pcla.23,Pcla.24,Pcla.a,Pcla12,Pcla14,Pcla17,Pcla25,
                   Pcla26,Pcla27,Pcla28,Pcla81,Pcla9)

levels <- data.frame(pop)

ind <- data.frame("MPL1","MPL10","MPL11","MPL12","MPL13","MPL14", "MPL15", "MPL16", "MPL17", "MPL18", "MPL19", "MPL2", "MPL20","MPL21","MPL22","MPL23","MPL24","MPL25",
                    "MPL26", "MPL27","MPL28","MPL29","MPL3","MPL30","ML4","MPL5","MPL6","MPL7", "MPL9", "MTA1","MTA10","MTA11","MTA12","MTA13","MTA14","MTA15", "MTA16","MTA17",
                    "MTA18", "MTA19","MTA2", "MTA20", "MTA21", "MTA22", "MTA23", "MTA24", "MTA25", "MTA26",
                    "MTA27", "MTA28", "MTA29", "MTA3", "MTA30", "MA32", "MTA33", "MTA34", "MTA4",
                    "MTA5", "MTA6", "MTA7", "MTA8", "MTA9", "MVA10", "MVA11", "MVA12", "MVA13", "MVA14", "MVA15",
                    "MVA16", "MVA17", "MVA18", "MVA19", "MVA20", "MVA21", "MVA22", "MVA23", "MVA24", "MVA25",
                    "MVA26", "MVA27", "MVA28", "MVA29", "MVA30", "MVA6", "MVA7", "MVA8", "MVA9") 

pop_pc <- data.frame("MPL","MVA","MTA")

#significance tests:

#tests the significance of the effect of level on genetic differentiation
global <-test.g(Pclavata_nomiss,level = pop, nperm=1000) #Fst global

#test individuals between subpops but only considering individual level
betweensubpop <-test.within(loci, test.lev = pop_pc, within = ind, nperm=1000)

#test for ind grouped together in subpops if they're different in the whole pop 
betweenpop <- test.between(loci, test = pop, rand.unit = pop_pc, nperm = 1000)

######## Discriminant Analysis of Principal Components (DAPC) #######
#aims to provide an efcient description of genetic clusters using a few 
#synthetic variables. These are constructed as linear combinations of the 
#original variables (alleles) which have the largest between-group variance 
#and the smallest within-group variance.

#provides membership probabilities of each individual for the diferent groups 
#based on the retained discriminant functions. While these are diferent from 
#the admixture coefcients of software like STRUCTURE, they can still be 
#interpreted as proximities of individuals to the diferent clusters. Membership
#probabilities also provide indications of how clear-cut genetic clusters
#are.

# PCA 
#Before running k-means to detect clusters, we transform the data using PCA:
#it reduces the number of variables so as to speed up the clustering algorithm.
#this does not imply a necessary loss of information since all the principal
#components (PCs) can be retained, and therefore all the variation in the
#original data.

Pclavata_nomiss
x <- tab(Pclavata_nomiss, freq=T, NA.method="mean")
pca.x <- dudi.pca(x, center=T, scale = F) #we select 3
pca.x$eig
pca.x$li
s.label(pca.x$li)
s.class(pca.x$li, fac = pop(Pclavata_nomiss), col=transp(funky(3),1), axesell = F,
        cstar=0, cpoint = 3)
add.scatter.eig(pca.x$eig[1:50],3,1,2, ratio=.2)

s.class(pca.x$li, fac = pop(Pclavata_nomiss), xax = 2, yax = 3, 
        col=transp(funky(15),.6), axesell = F,
        cstar=0, cpoint = 3)
add.scatter.eig(pca.x$eig[1:50],3,1,2, ratio=.2)

#PC1
pca.x$eig[1]
pc1 <- pca.x$li[,1]
var(pc1)
var(pc1)*86/87
mean(pc1^2)
n <- length(pc1)
0.5*mean(dist(pc1)^2)*((n-1)/n)
#0.4112776

#PC2
pca.x$eig[2]
pc2 <- pca.x$li[,2]
var(pc2)
var(pc2)*86/87
mean(pc2^2)
n <- length(pc2)
0.5*mean(dist(pc2)^2)*((n-1)/n)
#0.341314

#PC3
pca.x$eig[3]
pc3 <- pca.x$li[,3]
var(pc3)
var(pc3)*86/87
mean(pc3^2)
n <- length(pc3)
0.5*mean(dist(pc3)^2)*((n-1)/n)
#0.328363

eig.perc <- 100*pca.x$eig/sum(pca.x$eig)
head(eig.perc)

#allele contribution
loadingplot(pca.x$c1^1) #contribution of alleles in the variance component 1
loadingplot(pca.x$c1^2) #component 2

#DAPC

Pclavata_nomiss
pop(Pclavata_nomiss)

clus <- find.clusters(Pclavata_nomiss) #
#choosing 100 PCAs to retain
#BIC = 3

names(clus)
head(clus$Kstat, 10)
head(clus$stat, 10)

table(pop(Pclavata_nomiss), clus$grp)#to see how many individuals from each
#pop belong in each cluster: 1, 2 or 3
clus$size

table.value(table(pop(Pclavata_nomiss), clus$grp), col.lab=paste("inf", 1:3),
            row.lab=paste("ori", 1:3)) 

dapc1 <- dapc(Pclavata_nomiss, clus$grp)
#choosing 40 PCAs, where the curve flattens and little info is gained by adding more
#2 eigen value 

myCol <- funky(3)

scatter(dapc1, ratio.pca=0.1, bg="white", pch=20, cell=0,
        cstar=0, col=myCol, solid=0.7, cex=3, clab=0,
        mstree=TRUE, scree.da=FALSE, posi.leg="bottomright",
        leg=TRUE, txt.leg=paste("Cluster",1:3))

par(xpd=TRUE)
points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4,
       cex=1, lwd=5, col="black")
points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4,
       cex=1, lwd=2, col=myCol)

myInset <- function(){
temp <- dapc1$pca.eig
temp <- 100* cumsum(temp)/sum(temp)
plot(temp, col=rep(c("black","lightgrey"),
                   c(dapc1$n.pca,1000)), ylim=c(0,100),
     xlab="PCA axis", ylab="Cumulated variance (%)",
     cex=1, pch=20, type="h", lwd=2)
}
add.scatter(myInset(), posi="topright",
            inset=c(-0.03,-0.030), ratio=.10,
            bg=transp("white"))


#Calculate membership probabilities
class(dapc1$posterior)
dim(dapc1$posterior)
round(head(dapc1$posterior),3)
summary(dapc1)

#structure plot and membership probabilities
dapc2 <- dapc(Pclavata_nomiss, n.da=100, n.pca=100)
temp <- a.score(dapc2)
names(temp)

temp$tab[1:3,1:3]
temp$pop.score
temp$mean

dapc2 <- dapc(Pclavata_nomiss, n.da=100, n.pca=50)
temp <- optim.a.score(dapc2)
dapc2

dapc3 <- dapc(Pclavata_nomiss, n.da=100, n.pca=35)
myCol <- rainbow(15)

par(mar=c(4.1,4.1,2,2), xpd=TRUE)
compoplot(dapc3, lab=NULL, border = NA, posi=list(x=57,y=-.01), cleg=0.75)

#AMOVA

dist <- dist(Pclavata_nomiss)
strata <- strata(Pclavata_nomiss)
amova <- pegas::amova(dist ~ Pop, data = strata, nperm = 1000)
amova #Popultarion differences are sig pvalue = 0

#########

Amova <- poppr.amova(Pclavata_nomiss, ~Pop)
set.seed(1999)
Amovasignif   <- randtest(Amova, nrepet = 1000)
plot(Amovasignif)
Amovasignif #all significant

################# End of the code.
