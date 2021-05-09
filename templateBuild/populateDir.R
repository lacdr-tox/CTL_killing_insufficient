setwd('C:/Users/richard/Documents/PhD/Submission/Rev1/templateBuild')
templatedir <- 'configTemplates'
filelist <- list.files(templatedir)
templates <- c(growth=FALSE,
               n1=TRUE,
               n5=TRUE,
               IFNn1=TRUE,
               IFNn5=TRUE)
filelist<-filelist[templates]

allFolderNames <- ""

for (file in filelist){
  for(k in c(4,16)){
    for(edge in c(1,0)){
      infile <- file
      inParms <- paste(templatedir,infile,sep="/")
      inParms <- read.table(inParms,stringsAsFactors = F)
      inParms$V2[inParms$V1=='k']<-k
      inParms$V2[inParms$V1=='edge']<-edge
      infile <- strsplit(infile,split="Parms.csv")[[1]]
      infile<-paste(infile,"_k",toString(k),"_edge_",toString(edge),sep="")
      
      allFolderNames <- paste(allFolderNames,infile,sep=" ")
      
      outputdir <- paste(infile,"output",sep="/")
      outfile <- paste(infile,"parms.csv",sep="/")
      dir.create(infile)
      dir.create(outputdir)
      write.table(inParms,outfile,sep=" ",quote=FALSE, row.names=FALSE, col.names = FALSE)
    }
  }
}

infile <- 'growthParms.csv'
inParms <- paste(templatedir,infile,sep="/")
inParms <- read.table(inParms,stringsAsFactors = F)
infile <- strsplit(infile,split="Parms.csv")[[1]]
allFolderNames <- paste(allFolderNames,infile,sep=" ")
outputdir <- paste(infile,"output",sep="/")
outfile <- paste(infile,"parms.csv",sep="/")
dir.create(infile)
dir.create(outputdir)
write.table(inParms,outfile,sep=" ",quote=FALSE, row.names=FALSE, col.names = FALSE)

allFolderNames
  



