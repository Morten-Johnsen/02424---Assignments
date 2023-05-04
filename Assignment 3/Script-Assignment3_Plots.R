rm(list=ls())

library(reshape)
library(tidyverse)
library(magrittr)
library(dplyr)
library(ggplot2)
library(nlme)
library(gridExtra)
library(grid)
library(lattice)
library(ggpubr)

if (Sys.getenv('USER') == "mortenjohnsen"){
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/10. Semester/02424 - Advanced Dataanalysis and Statistical Modellling/02424---Assignments/Assignment 3/")
  figpath <- "/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/10. Semester/02424 - Advanced Dataanalysis and Statistical Modellling/02424---Assignments/Assignment 3/figs/"
}else if (Sys.getenv('USER') == "freja"){
  setwd("~/Documents/Uni/TiendeSemester/Adv. data analysis and stat. modelling/02424---Assignments/Assignment 3")
  figpath <- "~/Documents/Uni/TiendeSemester/Adv. data analysis and stat. modelling/02424---Assignments/Assignment 3/figs/"
}else{
  setwd("C:/Users/catdu/OneDrive/DTU/10. semester/Advanced Dataanalysis and Statistical Modelling/02424---Assignments/Assignment 3/")
  figpath <- "C:/Users/catdu/OneDrive/DTU/10. semester/Advanced Dataanalysis and Statistical Modelling/02424---Assignments/Assignment 3/figs"
}

clothing <- read.csv(file = "clothingFullAss03.csv", header = T)
head(clothing)

unique(clothing$day)
unique(clothing$time)
c.data <- dplyr::select(clothing, -time2, -X)
head(c.data)
     
# Not very pretty plot
melt(c.data,id = c('subjId','day') )%>%
  ggplot()+
  geom_point(aes(x = day, y=value,colour = subjId))+
  facet_wrap(~variable, scales = "free")+
  theme_bw()

mean.data <- c.data %>%
  group_by(subjId) %>%
  summarise(meanClo = mean(clo))
# Okay plot :)))) 
melt(c.data,id = c('subjId','clo'))%>%
  ggplot()+
  geom_point(aes(x = value, y=clo,colour = subjId), show.legend = FALSE)+
  facet_wrap(~variable, scales = "free")+
  theme_bw()

#Another plot looking at the 10 first subjects, clearly showing a difference between subjects depending on temperature
p1 <- melt(c.data,id = c('subjId','clo','day'))%>%
  filter(day == 1 & subjId %in% unique(c.data$subjId)[1:10]) %>%
  ggplot()+
  geom_point(aes(x = value, y=clo, colour = factor(subjId)))+
  geom_line(aes(x = value, y=clo, colour = factor(subjId), linetype = factor(day)), show.legend = FALSE)+
  facet_wrap(~variable, scales = "free")+
  theme_bw()+
  labs(colour = "SubjectID:", y = "Clothing insulation level", x = "Value")+
  theme(legend.position = "top")+
  ggtitle("Between subject differences in Clothing insulation level")
# p1
ggsave(filename = file.path(figpath, "subjIdDifferences.png"), plot = p1)


#   -----------------------------------------------------------------------
# Subset data to contain oberservations for the 10 first sujects
tiny.data <- c.data %>% filter(day == 1 & subjId %in% unique(c.data$subjId)[1:10])

# Calculate mean tOut, tInop and time for each subject for day 1
tiny.mean.data <- tiny.data %>%
  group_by(subjId) %>%
  summarise(meanClo = mean(clo),
            meantOut = mean(tOut),
            meantime = mean(time),
            meantInOp = mean(tInOp))

# print(tiny.mean.data)

p1 <- 
  ggplot(tiny.data, aes(x = tOut, y = clo, group = factor(subjId), color = factor(subjId)), size = 1) + 
  geom_line() +
  geom_point() +
  labs(colour = "SubjectID:", y = "Clothing insulation level", x = "Outdoor temperature")  +
  theme_minimal()
p1
p2 <- 
  ggplot(tiny.data, aes(x = tInOp, y = clo, group = factor(subjId), color = factor(subjId)), size = 1, show.legend = FALSE) + 
  geom_line() +
  geom_point() +
  labs(colour = "SubjectID:", x = "Indoor operating temperature", y = "") +
  theme_minimal() 
p2
p3 <- 
  ggplot(tiny.data, aes(x = time, y = clo, group = factor(subjId), color = factor(subjId)), size = 1, show.legend = FALSE) + 
  geom_line() +
  geom_point() +
  labs(colour = "SubjectID:", x = "Within day time measurement", y = "") +
  theme_minimal()
p3

pall <- grid.arrange(ggarrange(p1, p2, p3, ncol = 3, nrow = 1, common.legend = TRUE, legend="bottom"),
                     top = textGrob("     Between subject differences in Clothing insulation level",
                                    gp = gpar(fontsize = 12, font = 1),
                                    x = 0, hjust = 0))
ggsave(filename = file.path(figpath, "subjIdDifferences.new.png"), plot = pall, height = 5, width = 8)




############ histograms
c.data.temp <- data.frame(clo = c.data$clo[c.data$sex=="female"])
hist1 <- ggplot(c.data.temp, aes(x = clo)) + 
  geom_histogram(aes(y =..density..),
                 breaks = seq(min(c.data.temp), max(c.data.temp), by = 0.01), 
                 colour = "red", 
                 fill = "red",
                 alpha = 0.2) +
  stat_function(fun = dnorm, args = list(mean = mean(c.data.temp$clo), sd = sd(c.data.temp$clo))) +
  theme_minimal()+
  labs(x = "Clothing insulation level", subtitle = "Female")

c.data.temp <- data.frame(clo = c.data$clo[c.data$sex=="male"])
hist2 <- ggplot(c.data.temp, aes(x = clo)) + 
  geom_histogram(aes(y =..density..),
                 breaks = seq(min(c.data.temp), max(c.data.temp), by = 0.01), 
                 colour = "blue", 
                 fill = "blue",
                 alpha = 0.2) +
  stat_function(fun = dnorm, args = list(mean = mean(c.data.temp$clo), sd = sd(c.data.temp$clo))) +
  theme_minimal()+
  labs(x = "Clothing insulation level", subtitle = "Male")

library(ggpubr)
library(gridExtra)
pallhist <- grid.arrange(ggarrange(hist1, hist2, ncol = 2, nrow = 1, common.legend = TRUE, legend="bottom"),
                     top = textGrob("      Distribution of clothing insulation level",
                                    gp = gpar(fontsize = 12, font = 1),
                                    x = 0, hjust = 0))
pallhist
ggsave(filename = file.path(figpath, "sexDifferenceHist.png"), plot = pallhist, height = 5, width = 8)




#### Plot slopes and intercepts for random effects - color by subject ID:
intercept <- random.effects(final.model.ML1)$`(Intercept)`
slope <- random.effects(final.model.ML1)$tOut

alldf <- data.frame(cbind(slope, intercept, unique(c.data$subjId)))
names(alldf)[3] <- "subjId"

ggplot(data = alldf, aes(x, y)) +
  geom_point() +
  theme_minimal()+
  scale_y_continuous(limits = c(min(intercept), 0.4)) +
  geom_abline(data = alldf, aes(slope = slope, intercept = intercept, color = factor(subjId)))
