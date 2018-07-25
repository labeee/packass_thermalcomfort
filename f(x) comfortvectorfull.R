###--- Script for performing vectorized operations with PMV, PPD and SET with
###--- some input variables data. The input files must be in csv format with the
###--- respective order of colunms:
#---CLO clothing insulation [clo]
#---MET metabolic rate [met]
#---WME mechanical work [met]
#---TA air temperature [oC]
#---TR mean radiant temperature [oC]
#---VEL air speed [m/s]
#---RH relative humidity [%]
#---PA partial water vapor pressure [kPa] - Default to 0, to use RH
#---APA atmospheric pressure [kPa] - Default to 101.325 kPa (1 atm)

#---Examples of input and output arguments for testing
#input=c("C:/Users/ruhtr/OneDrive/Trabalho_2018-07_Comfort algorithms/validation calculation packass/input")
#output=c("C:/Users/ruhtr/OneDrive/Trabalho_2018-07_Comfort algorithms/validation calculation packass/output")

#comfortvectorfull(input,output)
comfortvectorfull=function (input,output){
      
      time.initial=proc.time()
      
      #---Set the working directory to this file path, to ensure the correct linking
      setwd(dirname(rstudioapi::getSourceEditorContext()$path))
      #---Open thermal comfort functions. Link to the respective folder.
      source("f(x) pmvboth.R")
      source("f(x) setashrae.R")
      source("f(x) optemplimits.R")
      
      #---Get file names in full or short length
      name.full=list.files(path=input,pattern="*.csv",full.names=TRUE)
      name.short=list.files(path=input,pattern="*.csv")
      
      #---Get the colunm length of the first file
      num.files=length(x=name.full)
      get.col=length(read.table(file=name.full[1],header=TRUE,sep=",",
                               as.is=TRUE,na.strings=TRUE))
      
      #---Set the class of all colunms to numeric
      classes.n=NULL
      classes.n[1:get.col]="numeric"
      i=1
      
      #--------Predicted Percentage Dissatisfied - PPD-----------#
      getPPD=function(insertPMV)
      {
            PPD.out=100-95*exp(-0.03353*insertPMV^4-0.2179*insertPMV^2)
            return(PPD.out)
      }
      
      #--------Operative temperatura formula - ASHRAE Standard 55 (2017) Appendix A-------#
      getTO=function(airtemp,radtemp,airvelocity)
      {
        if(airvelocity<0.2) Acoeff=0.5
        if(airvelocity>=0.2 && airvelocity<0.6) Acoeff=0.6
        if(airvelocity>=0.6) Acoeff=0.7 #---the correct limit would be 1.0m/s...
        
        totempcalc=round(Acoeff*airtemp+(1-Acoeff)*radtemp,digits=2)
        return(totempcalc)
      }
      
      #---Column names for the output files (different for each case)
      col.1=c("PMV")
      col.2=c("PPD [%]")
      col.3=c("SET [oC]")
      col.4=c("TO [oC]")
      col.5=c("TO inferior [oC]")
      col.6=c("TO superior [oC]")
      column.names=c(col.1,col.2,col.3,col.4,col.5,col.6)
      
      #---Name of the output folder
      folder.output=output
      
      #---Loop in the input folder files
      for (i in 1:num.files) #---change to num.files
      {
            
            file.i=read.table(file=name.full[i],header=TRUE,sep=",",
                              as.is=TRUE,na.strings=TRUE,colClasses=classes.n)
            
            num.col=ncol(file.i)
            num.row=nrow(file.i)
            
            #---input(CLO,MET,WME,TA,TR,VEL,RH,PA,APA)
            CLO.j=file.i[,1]
            MET.j=file.i[,2]
            WME.j=file.i[,3]
            TA.j=file.i[,4]
            TR.j=file.i[,5]
            VEL.j=file.i[,6]
            RH.j=file.i[,7]
            PA.j=file.i[,8] #---Use default
            APA.j=file.i[,9] #---Use default
            
            each.output=NULL
            
            PMV.PPD.out=mapply(pmvboth,CLO.j,MET.j,WME.j,TA.j,TR.j,VEL.j,RH.j,PA.j,APA.j)
            each.output[[1]]=PMV.PPD.out[1,]  #---PMV
            each.output[[2]]=PMV.PPD.out[2,]  #---PPD
            print("PMV and PPD calculation complete")
            
            SET.ashrae=mapply(setashrae,CLO.j,MET.j,WME.j,TA.j,TR.j,VEL.j,RH.j,APA.j)
            each.output[[3]]=SET.ashrae
            print("SET calculation complete")
            
            TO.calc=mapply(getTO,TA.j,TR.j,VEL.j)
            each.output[[4]]=TO.calc
            print("TO calculation complete")
            
            TO.limit=mapply(optemplimits,CLO.j,MET.j,WME.j,VEL.j,RH.j,PA.j,APA.j)
            each.output[[5]]=TO.limit[1,]
            each.output[[6]]=TO.limit[2,]
            print("Operative temperature limits calculation complete")
            
            #---Output files
            output.file=data.frame(each.output)
            names(output.file)=column.names
            
            name.output=as.character(paste(folder.output,"/Output ",name.short[i],sep=""))
            write.table(output.file,file=name.output,row.names=FALSE,col.names=TRUE,sep=",")
            
            print(paste("file ",i," finished",sep=""))
      }
      
      time.final=proc.time(); time.elapsed=time.final-time.initial
      print(paste("Elapsed time= ",round(time.elapsed[3],4)," seconds"),sep="")
      print(paste("Elapsed time= ",round(time.elapsed[3]/60,4)," minutes"),sep="")
      print(paste("Elapsed time= ",round(time.elapsed[3]/3600,4)," hours"),sep="")
      
      return(print("Done, end script"))
    
}