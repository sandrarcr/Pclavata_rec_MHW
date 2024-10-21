####################################################################
#### Code accompanying:
#
#   Ramirez-Calero, S. et al. 2024. Recurrent extreme climatic events are driving gorgonian populations 
#   to local extinction: low adaptation and low adaptability to marine heat waves.
#
# Script written by: Sandra Ramirez
#
#This script was run to obtain the average tissue necrosis of each population studied. Data sets can be 
#found in repository
#########################################################################################################

#set working directory
setwd("~")

#libraries
library(tidyr)
library(dplyr)                             
library(ggplot2)
library(stringr)

#### This script is run separately for each population ####

# The data for all populations contains the values of average tissue necrosis of all P. clavata
#colonies collected for each day, from 1 to 28. Since all population datasets are the same, the script applies to all

#load data
data <- read.csv("data_ave_nec_[popname_years].csv", header = T, check.names = F)#data files are included in the repository

#data manupulation
spxx <- select(data,-contains("Error"))
error.xx <- select(data,Day,contains("Error"))

colnames(error.xx) <- colnames(spxx)

error.xx2 <- error.xx %>%
  tidyr::gather(key="sp",value="error",-Day)

sp.xx <- spxx %>%
  tidyr::gather(key="sp",value="nec",-Day)

sp.wholex <- cbind(sp.xx,select(error.xx2,error))

#plot

p= sp.wholex %>%
  ggplot(aes(Day,nec)) +                                  
  geom_point(aes(color=sp,shape=sp), size = 3) +
  geom_line(aes(group=sp,color=sp), size = 2) +
    geom_ribbon(aes(ymin = nec - error,
                   ymax = nec + error,
                  group=sp, fill=sp), alpha = 0.2) +
  scale_color_manual(values = c("forestgreen", "darkblue", "orange", "purple"))+
  scale_fill_manual(values=c("forestgreen", "darkblue", "orange", "purple"))+
  scale_shape_manual(values = c(4,4,4,25,4,15,25,25,15,4,25,15))+
  scale_y_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100),limits=c(0,100))+
  scale_x_continuous(breaks=c(0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30), limits = c(0,30))+
  ggtitle("La Vaca")+
  ylab("Mean extent of injury") +
  xlab("Time (days)")+
  theme(plot.title = element_text(hjust = 0.5,size=30))+
  theme(axis.text.y=element_text(colour="black", size=24), 
        axis.title = element_text(size = 28))+
  theme(axis.text.x=element_text(colour="black", size=24),
        strip.text.x=element_text(size=28,face="bold"),
        strip.background=element_rect(color="black",fill="gainsboro",
                                      size=1.1))+
  theme(legend.key=element_blank())+
  theme(legend.title=element_blank())+
  theme(legend.text = element_text(size=18))+
  theme(axis.line.y = element_line(),
        axis.line.x=element_line(),
        panel.grid.major=element_blank(),
        panel.border=element_rect(fill=NA),
        panel.background=element_blank())

p #display plot

dev.off()

#Merge all files for each populatio in a single pdf

pdf("Average tissue necrosis.pdf", width = 15, height = 15, paper = "USr")

print(p1)
print(p2)
print(p3)

dev.off()

################# End of the code.
