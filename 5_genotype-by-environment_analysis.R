####################################################################
#### Code accompanying:
#
#   Ramirez-Calero, S. et al. 2024. Recurrent extreme climatic events are driving gorgonian populations 
#   to local extinction: low adaptation and low adaptability to marine heat waves.
#
# Script written by: Jean Baptiste Ledoux and Sandra Ramirez-Calero
#
#This script contains information on the estimation of genotype-by-environment sensitivity of
#P. clavata individuals. This script needs to be run after 2_PCAs_GLM.R
#########################################################################################################

#setwd
setwd("~")

#load libraries
library(dplyr)
library(ggplot2)

#load data
d <- read.table("Pcla_interc_het") #can be found on repository
str(d)

#plot

pop_colors <- c("POTA" = "#a6cee3", "VACA" = "#b15928", "TASCONS" = "#510066")

d %>%
  ggplot(aes(x = Pop, y = PCScores, group = Ind, col = Ind)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", fill = NA, size = 0.71) +
  ylim(c(-4.5, 9)) +
  xlim(c(-2.5, 5))+
  labs(x = "Environmental value (mean PCA scores)", y = "Individual phenotypic response") +
  theme_minimal() +  
  theme(axis.text = element_text(size = 12),
        legend.position = "none",  # Remove the legend
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_text(size = 14),  
        axis.title.y = element_text(size = 14))
 scale_color_manual(values = pop_colors)

################# End of the code.
