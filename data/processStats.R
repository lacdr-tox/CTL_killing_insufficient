
## rscript to summarise statistics from multiple runs for plotting
## we plot volume, conjugates, realised killing rate, number of ctls. 
library(zoo)
## set target

inputDir <- "/home/richard/Documents/000_phd_archive/can_res_19/data/stats/"
outputDir <- "/home/richard/Documents/000_phd_archive/can_res_19/data/processedStats/"

files <- list.files(inputDir)
#files <- tail(list.files(inputDir),1)
#files<-c("mixed")

minMaxMed <- function(myVec){
  
  myVec[is.na(myVec)]<-0
  med<-mean(myVec)
  max<-max(myVec)
  min<-min(myVec)  
  
  return(c(min,max,med))
  
}

ratio <- function(num,denom,scale){
  
  myRat <- scale*num/denom
  stats <- minMaxMed(myRat)
  
  return(stats)
  
}

getStats <- function(i,mylist){
  
  
  vol <- sapply(mylist, function(x) x$N[i]+x$dead[i]) 
  tcells <- sapply(mylist, function(x) x$CTL[i]) 
  con <- sapply(mylist, function(x) x$con[i])
  deaths <-  sapply(mylist, function(x) abs(x$deaths[i]))
  rad <-  sapply(mylist, function(x) abs(x$tumrad[i]))
  
  confrac <- con/tcells
  realised <- deaths/tcells
  measured <- deaths/con
  
  
  volS <- minMaxMed(vol)
  tcellsS <- minMaxMed(tcells)
  conS <- minMaxMed(con)
  deathsS <- minMaxMed(deaths)
  radS <- minMaxMed(rad)
  confracS <- minMaxMed(confrac)
  realisedS <- minMaxMed(realised)
  measuredS <- minMaxMed(measured)
  
  
  time<-i/1440
  
  
  out <- c(time,volS,tcellsS,conS,deathsS,radS,confracS,realisedS,measuredS)
  names(out)<-c("time",
                "volmin","volmax","volmed",
                "tcellmin","tcellmax","tcellmed",
                "conmin","conmax","conmed",
                "deathmin","deathmax","deathmed",
                "radmin","radmax","radmed",
                "confracmin","confracmax","confracmed",
                "realisedmin","realisedmax","realisedmed",
                "measuredmin","measuredmax","measuredmed")
  
  return(out)
  
}

processr <- function(name){
  outdir <- paste(outputDir,name,".csv",sep="")
  indir <- paste(inputDir,name,sep="")
  setwd(indir)
  filelist <- list.files()
  files <- lapply(filelist,read.csv)
  
  ## check all data frames are same length:
  nrows <- sapply(1:length(files), function(i) length(files[[i]]$time))
  
  ## handling of different size data frames (early terminating sims or unfinished ones)
  nrows <- min(nrows) ## cut all data to the shortest (may prefer to exclude)
  #files <-lapply(files,head,n=nrows)
  
  vdata <- lapply(1:nrows,getStats,mylist=files)
  vdata <- data.frame(do.call(rbind,vdata))
  vdata$id<-rep(name,length(vdata$time))
  
  write.csv(vdata,outdir,row.names = FALSE)
}

require(parallel)
mclapply(files,processr,mc.cores = length(files))

