####################################################################
#### Code accompanying:
#
#   Ramirez-Calero, S. et al. 2024. Recurrent extreme climatic events are driving gorgonian populations 
#   to local extinction: low adaptation and low adaptability to marine heat waves.
#
# Script written by: Sandra Ramirez
#
#This script contains information on how heterosis was obtained from microsatellite and necrosis data (genotype & trait) to estimate fitness. 
#This script needs to be run after 2_PCAs_GLM.R and 3_microsat_analysis
#Cleaned data is provided in the repository. 
#
#########################################################################################################

#load packages:
library("lme4")
library("inbreedR")
library("tidyr")

#set working directory
setwd(dir = "~") 

#data
pcla_data <- (read.csv(file = "P.clavata_microsats.csv", row.names = "Ind"))
ind_nec_year <- (read.csv(file = "Avnec_het.csv", row.names = "Ind")) #average necrosis at the end of the experiment per individual

#define vectors to be used as traits:
nec_2015 <- ind_nec_year$Nec_2015
nec_2016 <- ind_nec_year$Nec_2016
nec_2017 <- ind_nec_year$Nec_2017

#convert in inbreed format

pcla <- convert_raw(pcla_data)
write.csv(pcla,file = "pcla_format_inbreedR.csv")

#correlation of extent to which heterozygosities are correlated across pairs of loci. 
#as a proxy to characterize the distribution of f (inbreeding level) in populations
#calculates confidence intervals by bootstrapping over individuals. 
#It also permutes the genetic data to generate a P-value for the null hypothesis of 
#no variance in inbreeding in the sample (i.e. g2 = 0).

pcla_g2 <- g2_microsats(pcla, nperm = 1000, nboot = 1000, CI = 0.95)
pcla_g2 #To display a summary of the results just print the output of an inbreedR function.

#plot shows the distribution of bootstrap results including the confidence interval.

plot(pcla_g2, main = "Microsatellites",
     col = "cornflowerblue", cex.axis=0.85)

HHC_pcla <- HHC(pcla, reps = 1000) # heterozygosity-heterozygosity correlation coefficients (HHCs)
HHC_pcla

plot(HHC_pcla, main = "Microsatellites",
     col = "cornflowerblue", cex.axis=0.85) #Distribution of heterozygosity-heterozygosity correlations

### Assuming that HFC (heterozygosity-fitness correlations) we can 
#calculate both the expected correlation between heterozygosity and inbreeding level

# r^2 between inbreeding and heterozygosity 
hf <- r2_hf(genotypes = pcla, type = "msats", nboot = 1000, parallel = F)
hf
plot(hf)

# r^2 between inbreeding and fitness
#check for values of necrosis for each marker, not for individual. Dimension error.

Wf1 <- r2_Wf(genotypes = pcla, trait = nec_2015, 
            family = gaussian, type = "msats", nboot = 1000, parallel = F, CI = 0.95)

Wf2 <- r2_Wf(genotypes = pcla, trait = nec_2016,
            family = gaussian, type = "msats", nboot = 1000, parallel = F, CI = 0.95)

Wf3 <- r2_Wf(genotypes = pcla, trait = nec_2017, 
            family = gaussian, type = "msats", nboot = 1000, parallel = F, CI = 0.95)

r2_Wf(genotypes = pcla, trait = nec_2015, 
      family = gaussian, type = "msats", nboot = 1000, parallel = F, CI = 0.95)

#####  estimating the impact of inbreeding on fitness using HFC
# g2
g2 <- g2_microsats(pcla)
# calculate sMLH
het <- sMLH(pcla)
# variance in sMLH
het_var <- var(het)

# Linear model of fitness trait on heterozygosity
pcla_nec <- cbind(pcla, ind_nec_year[, 1]) #add necrosis to genotypes 2015
names(pcla_nec)[18] <- "nec_2015" #rename column
mod_2015 <- lm(pcla_nec$nec_2015 ~ het)

pcla_nec <- cbind(pcla, ind_nec_year[, 2]) #add necrosis to genotypes 2016
names(pcla_nec)[18] <- "nec_2016" #rename column
mod_2016 <- lm(pcla_nec$nec_2016 ~ het)

pcla_nec <- cbind(pcla, ind_nec_year[, 3]) #add necrosis to genotypes 2017
names(pcla_nec)[18] <- "nec_2017" #rename column
mod_2017 <- lm(pcla_nec$nec_2017 ~ het)

# regression slope
beta_15 <- coef(mod_2015)[2]
beta_16 <- coef(mod_2016)[2]
beta_17 <- coef(mod_2017)[2]

# r2 between fitness and heterozygosity
predict <- data.frame(predict(mod_2015))
pcla <- na.omit(pcla)
Wh <- cor(pcla, predict, use="complete.obs")^2 

predict1 <- data.frame(predict(mod_2016))
pcla <- na.omit(pcla)
Wh <- cor(pcla, predict, use="complete.obs")^2 

predict2 <- data.frame(predict(mod_2017))
pcla <- na.omit(pcla)
Wh <- cor(pcla, predict, use="complete.obs")^2 
