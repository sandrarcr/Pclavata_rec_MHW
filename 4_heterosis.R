####################################################################
#### Code accompanying:
#
#   Ramirez-Calero, S. et al. 2024. Recurrent extreme climatic events are driving gorgonian populations 
#   to local extinction: low adaptation and low adaptability to marine heat waves.
#
#This script contains information of how heterosis was explored in the data
#This script needs to be run after 2_PCAs_GLM.R and 2_microsat_analysis
#Cleaned data is provided in the repository. cleaning process was conducted to match individuals that
#had a microsatellite present and that obtained an intercept after making sure necrosis data was present for the three years only
#individuales with data for two years only was removed.
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
pcla_data <- (read.csv(file = "P.clavata_heterosis_clean.csv", row.names = "Ind"))
intercepts <- (read.csv(file = "Intercepts_ind_3years_clean.csv", row.names = "ind"))



#obtain Multilocus heterozigosity
het <- sMLH(pcla) 
#conver vector in dataframe:
het_df <- data.frame(het = het, stringsAsFactors = FALSE)
# Convert het_df dataframe to a format that ggplot2 can use
het_df <- data.frame(colony = rownames(het_df), MLH = het_df[,1])

# Create a scatter plot with gene names as labels
ggplot(het_df, aes(x = colony, y = MLH, label = colony)) +
  geom_point() +
  geom_text(hjust = 1, vjust = 1, size = 2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Individuals", y = "sMLH") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

#set factor to output data
het_vect <- het_df$MLH
#build dataframe
data <- data.frame(intercepts, het_vect)
#write.csv(data,file = "pcla_intercept_het_linearR.csv")

#run linear model to test for heterosis and individual phenotypic response
fit <- lm(het_vect ~ intercepts, data = data)
summary(fit)

# create scatter plot
plot(data$intercepts, data$het,
     xlab = "nec-int",
     ylab = "sMLH",
     main = "Multilocus Heterozygosity vs. Phenotypic response",
     xlim = c(-1.2, 1.2), ylim = c(0, 1.6))

# add trend line
abline(fit, col = "red")

# add text with p-value and R-squared
text(x = 0.6, y = 1.5, 
     labels = paste("R-squared =", round(summary(fit)$r.squared, 1),
                    "p-value =", round(summary(fit)$coefficients[2, 4], 2)),
     pos = 3)

################# End of the code.
