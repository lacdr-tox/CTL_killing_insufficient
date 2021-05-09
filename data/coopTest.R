setwd("C:/Users/richard/Documents/PhD/Submission/Rev1/data/coocupancyTest")
files <- list.files()

#ctl <- read.table(files[1],sep=",")
## can easily check for multiple occupancy, since duplicated returns 
## ONLY THE SECOND, THIRD... Nth DUPLICATE. Therefore we can iteratively remove unique cells
#ctlstrings <- sapply(1:nrow(ctl),function(i)  paste(ctl[i,1:4],collapse=""))
dupstrings<-ctlstrings

nctl <- length(ctlstrings)
for( i in 1:10){
  dupstrings <- dupstrings[duplicated(dupstrings)]
  ndups <- length(dupstrings)/nctl
  print(ndups)
}

