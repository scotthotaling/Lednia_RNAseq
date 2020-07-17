###################################################################################
# R code for 2018 work on Lednia tetonica etc., created by A. A. Shah. 
# Title: "Mountain stoneflies may tolerate warming streams: evidence from organismal physiology and gene expression" 
# Journal: Global Change Biology, 2020
# Authors: Scott Hotaling, Alisha A. Shah, et. al.
###################################################################################

#Set working directory, read in data
setwd("PATH TO DATA")
data <- read.csv("Teton_CTmax.csv")

#Check data
head(data)
tail(data)
str(data)

####### TETONICA ONLY / MAIN ANALYSIS #######

#Subset data to include only Lednia tetonica
tet<-subset(data, data$species=="L. Tetonica")
#Check data
head(tet)
str(tet)

#Is there an effect of body size on CTmax in L. tetonica?
size.test1 <- lm(ctmax_C ~ length_mm, data = tet)
summary(size.test1) #size effect is non-significant
#Visualize
plot(tet$ctmax_C~tet$length_mm)

#Because there is no effect of size, final model will only include CTmax and Tmax. Here we ask, do insects experiencing higher maximum stream temperatures (Tmax) have higher heat tolerance (CTmax)?

m1 = lm(ctmax_C~t_max, data=tet)
summary(m1) #Tmax is a significant predictor of CTmax, relationship is positive. 
anova(m1)
#Visualize
require(ggplot2)
require(dplyr)
#Get means
gd <- tet %>% 
  group_by(t_max) %>% 
  summarise(ctmax_C = mean(ctmax_C))
gd
#Plot data and overlay means
ggplot(tet, aes(x = t_max, y = ctmax_C)) +
  geom_point(alpha = .4) +
  geom_point(data = gd, size = 4)+
  geom_smooth(method='lm', se=FALSE)

####### INCLUDE ZAPADA AND L. TUMANA / SUPPLEMENTARY ANALYSIS #######

#Use data = data

#Is there an effect of body size on CTmax in the complete dataset?
size.test2 <- lm(ctmax_C ~ length_mm, data = data)
summary(size.test2) #size effect is non-significant
#Visualize
plot(data$ctmax_C~data$length_mm)

#As with Tetnica only, there is no effect on size so the final model will only include CTmax and Tmax. Here we ask, do insects experiencing higher maximum stream temperatures (Tmax) have higher heat tolerance (CTmax)?

m2 = lm(ctmax_C~t_max, data=data)
summary(m2) #Tmax is a significant predictor of CTmax 
anova(m2)
#Visualize
require(ggplot2)
require(dplyr)
#Get means
gd <- data %>% 
  group_by(t_max) %>% 
  summarise(ctmax_C = mean(ctmax_C))
gd
#Plot data and overlay means
ggplot(data, aes(x = t_max, y = ctmax_C)) +
  geom_point(alpha = .4) +
  geom_point(data = gd, size = 4)+
  geom_smooth(method='lm', se=FALSE)
