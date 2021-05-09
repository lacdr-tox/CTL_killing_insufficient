library(reshape2)
library(ggplot2)
library(deSolve)
library(plyr)






times <- c(7, 10,  12,  14)
mouse_one <- c(1.072,0.134,0.0565,0.01)
mouse_two <- c(1.072,0.0167552,0.0565488,0.01)
mouse_three <-c(2.788, 0.1340416, 0.0565488,0.01)
mouse_four <- c(3.619,0.2618,0.01,0.01)

mouse_five<-c(2.7876464, 28.9529856, 56.5488, 93.7010998)
mouse_six<-c(0.1340416, 14.3654896, 32.725, 166.5194608)
mouse_seven<-c(0.1340416, 12.2145408, 56.5488, 277.8242544)

mousedata <- data.frame(times,mouse_one, mouse_two, mouse_three, mouse_four, mouse_five, mouse_six, mouse_seven)

mousedata <- melt(mousedata, id.vars = "times")
group <- c(rep("regression", 16), rep("growth",12))

mousedata <- cbind(mousedata, group)

means <- mousedata
means$value<-log(means$value)
means <- aggregate(means[, 3], list(means$times,means$group), mean)
names(means)<-c("time","group","volume")
means$volume<-exp(means$volume)

write.csv(mousedata,'/home/richard/tumourReg/writing/mousedata/mousedata.csv',row.names = F)
write.csv(means,'/home/richard/tumourReg/writing/mousedata/means.csv',row.names = F)