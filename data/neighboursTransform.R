sl <- function(i) subset(read.csv(paste(dir,"sum_",toString(i),".csv",sep="")),xor(z==0,xor(z==-1,z==1)))

i <- (60*24)*10

slice <-sl(i)

#d <- function(i,j) sqrt(sum((slice[i,4:6]-slice[j,4:6])^2))

maxrad <- max(slice[,4:5])
minrad <- min(slice[,4:5])
rad <- maxrad-minrad

slice[,4:5]<-slice[,4:5]-minrad+1
slice[,6]<-slice[,6]+2


maxrad <- max(slice[,4:5])
minrad <- min(slice[,4:5])

dummy <- array(0,c(maxrad+2,maxrad+2,5))
dummy2<- array(0,c(maxrad+2,maxrad+2,5))


for(i in 1:length(slice$time)){
  xyz<-as.numeric(slice[i,4:6])+1
  
  dummy[xyz[1],xyz[2],xyz[3]]<-1
}

image(dummy[,,3])

for(i in c(-1,0,1)){
  for(j in c(-1,0,1)){
    for(k in c(-1,0,1)){
      if(i!=0||j!=0||k!=0){
        
        dummy2[2:(maxrad+1),2:(maxrad+1),2:4] <- dummy2[2:(maxrad+1),2:(maxrad+1),2:4] +
          dummy[(i +2:(maxrad+1)),(j +2:(maxrad+1)),(k+2:4)]
        
        
      }
      
    }
  }
}

image(dummy2[,,3])

centre <- subset(slice,z==2)
neighbours <- centre$time*0

for(i in 1:length(centre$time)){
  xyz<-as.numeric(centre[i,4:6])+1
  neighbours[i] <- dummy2[xyz[1],xyz[2],3]
}

centre$neighbours <- 26-neighbours



p<-ggplot(data=centre,aes(x=x,y=y,fill=neighbours))+
  geom_raster()+
  scale_fill_continuous("cell age",low="red",high="green")+
  theme(text = element_text(size=20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.title = element_blank())
p
