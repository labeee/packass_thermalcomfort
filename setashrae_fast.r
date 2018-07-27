#Script for SET thermal comfort index calculation using any number of cases or csv files

#---Set the working directory to this file path, to ensure the correct linking
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
code_path <- getwd()

source("f(x) setashrae.R")
input_path <- paste0(code_path,"/datasets/setashrae_fast_input")
output_path <- paste0(code_path,"/datasets/setashrae_fast_output")

#---Get file names in full or short length
name.full=list.files(path=input_path,pattern="*.csv",full.names=TRUE)
print(name.full)
name.short=list.files(path=input_path,pattern="*.csv")
print(name.short)

#---Get the colunm length of the first file
num.files=length(x=name.full)
print(num.files)
get.col=length(read.table(file=name.full[1],header=TRUE,sep=",",
                          as.is=TRUE,na.strings=TRUE))

#---Set the class of all colunms to numeric
classes.n=NULL
classes.n[1:get.col]="numeric"
i=1

#---Loop in the input folder files
for (i in 1:num.files) #---change to num.files
{
  
  file.i=read.table(file=name.full[i],header=TRUE,sep=",",
                    as.is=TRUE,na.strings=TRUE,colClasses=classes.n)
  
  num.col=ncol(file.i)
  num.row=nrow(file.i)
  
  #---input(CLO,MET,WME,TA,TR,VEL,RH,PA,APA)
  CLO.j=file.i[,"CLO"]
  MET.j=file.i[,"MET"]
  WME.j=file.i[,"WME"]
  TA.j=file.i[,"TA"]
  TR.j=file.i[,"TR"]
  VEL.j=file.i[,"VEL"]
  RH.j=file.i[,"RH"]
  PA.j=file.i[,"PA"]
  APA.j=101.325 #---Use default
  
  SET.ashrae=mapply(setashrae,CLO.j,MET.j,WME.j,TA.j,TR.j,VEL.j,RH.j,APA.j)

  #---Output files
  output.file=data.frame(SET.ashrae)
  name.output=as.character(paste0(output_path,"/output_",name.short[i]), sep ="")
  write.table(output.file,file=name.output,row.names=FALSE,col.names=TRUE,sep = ",")
  
  print(paste("file ",i," finished",sep=""))
}
