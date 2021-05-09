## generate a build folder for the paper figs. 
## option to setwd manually:
setwd('C:/Users/richard/Documents/PhD/Submission/Rev1/templateBuild')

## or set via command line argument
usrSetArgs=FALSE
if(!usrSetArgs) "setting args manually: useage Rscript thisScript /path/to/build"
args <- commandArgs(usrSetArgs)
if(usrSetArgs) setwd(args[1])

defaultParms <- read.table('parmsPaperDefault.csv',stringsAsFactors = F)

#k16 and k4 n1
#k16 and k4 n5
#k16 and k4 n1 edge
#k16 and k4 n5 edge
#k16 and k4 n1      ifn
#k16 and k4 n5      ifn
#k16 and k4 n1 edge ifn
#k16 and k4 n5 edge ifn

## since we are only changing 3 variables this an be done in loops

## careful because different configurations have different impacts on other settings

## if ctl==false the grid should be 1000
# if ctl== true then 800
# if ifn is true the grid can be 600 and nthreads should increase
# if n5 then pdet =0.1 and pcon=0.5




varyPars <- c('k','hits','edge','ifn')
inRange <- list(k=c('4','16'),
                hits=c(1,5),
                edge=c(1,0),
                ifn=c(1,0))

allFolders <- ""

for(i in 1:varyPars){
  subset<-varyPars[1:i]
  for(par in subset){
    for(range in inRange){
      newParmSheet <- defaultParms
      
    }
  }
}


n1parms <- defaultParms
n1parms$V2[n1parms$V1=="gsize"]<-801


