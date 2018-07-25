###--- Script that calculates 37 long-term thermal comfort indices
###--- The input arguments are: 'seasons', which contains a vector of the hours
###--- of the year that belongs to each season (summer and winter). A value of 0 (zero)
###--- represents the winter season, and value 1 (one) represents the summer season.
###--- the length of the vector should be equal to the length of the input data.

#--- The above variable is an example of 'seasons' argument, with length of 3120 occupied hours.
#seasons=vector(mode="numeric",length=3120)
#seasons[1:1296]=1; seasons[2341:3120]=1

#input=c("Long-term indices/input")

#---Function of the long-term comfort indices
comfortindices=function (input,output,seasons){
      
      #---Set the working directory to this file path, to ensure the correct linking
      setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  
      #---Transformation in summer and winter
      id.summer=seasons; id.winter=seasons
      id.winter[id.winter==0]=2; id.winter[id.winter==1]=0; id.winter[id.winter==2]=1
      
      #---Link to the folder containing the data set
      name.full=list.files(path=input,pattern="*.csv",full.names=TRUE)
      name.short=list.files(path=input,pattern="*.csv")
      
      #---Definition of numeric variables
      num.files=length(x=name.short)
      classes.n=NULL; classes.n[1:7]="numeric"
      
      #---Knowing the lenght of input data
      file.1=read.table(file=name.full[1],header=TRUE,sep=",",
                        as.is=TRUE,na.strings=TRUE,colClasses=classes.n)
      length.f1=nrow(x=file.1)
      
      #---creates a vector with the sequential occupied hours
      id.hour=c(1:length.f1)
      
      #---names of the long-term indices
      names.columns=c("i1","i2","i3","i4","i5","i6","i7","i8","i9","i10","i11","i12","i13","i14",
                      "i15","i16","i17","i18","i19","i20","i21","i22","i23","i24","i25","i26","i27","i28",
                      "i29","i30","i31","i32","i33","i34","i35","i36","i37","PeopleLoad")
      
      #---initial variables
      res=data.frame(row.names=name.short)
      res.norm=data.frame(row.names=name.short)
      
      i=1
      for (i in 1:num.files)
      {
            file.i=read.table(file=name.full[i],header=TRUE,sep=",",
                              as.is=TRUE,na.strings=TRUE,colClasses=classes.n)
            
            n.col=ncol(file.i)
            n.row=nrow(file.i)
            
            #index 1 - [POR(PMV) WY]
            index.1=round(100*(sum(rowSums(file.i[1]>0.5,na.rm=TRUE))+
                                     sum(rowSums(file.i[1]<(-0.5),na.rm=TRUE)))/n.row,digits=1)
            
            #index 2 - [POR(PMV) SS warm]
            index.2=round(100*(sum(rowSums(file.i[1]*id.summer>0.5,na.rm=TRUE)))/n.row,digits=1)
            
            #index 3 - [POR(PMV) SS cold]
            index.3=round(100*(sum(rowSums(file.i[1]*id.winter<(-0.5),na.rm=TRUE)))/n.row,digits=1)
            
            #index 4 - [POR(PMV) SS Int]
            index.4=index.2+index.3
            
            #index.5 - [POR(OP) WY]
            index.5=round(100*(sum(rowSums(file.i[7]>file.i[5],na.rm=TRUE))+
                                     sum(rowSums(file.i[7]<file.i[4],na.rm=TRUE)))/n.row,digits=1)
            
            #index.6 - [POR(OP) SS warm]
            index.6=round(100*(sum(rowSums(file.i[7]*id.summer>file.i[5],na.rm=TRUE)))/n.row,digits=1)
            
            #index.7 - [POR(OP) SS cold]
            index.7.temp=file.i[7]; index.7.temp[[1]][id.winter==0]=NA
            index.7=round(100*(sum(rowSums(index.7.temp<file.i[4],na.rm=TRUE)))/n.row,digits=1)
            
            #index.8 - [POR(OP) SS Int]
            index.8=index.6+index.7
            
            #index 9 - [Exceed(PMV) WY]
            index.9=round(sum(rowSums(file.i[1]>0.5,na.rm=TRUE))+
                                sum(rowSums(file.i[1]<(-0.5),na.rm=TRUE)),digits=0)
            
            #index 10 - [Exceed(PMV) SS warm]
            index.10=round(sum(rowSums(file.i[1]*id.summer>0.5,na.rm=TRUE)),digits=0)
            
            #index 11 - [Exceed(PMV) SS cold]
            index.11=round(sum(rowSums(file.i[1]*id.winter<(-0.5),na.rm=TRUE)),digits=0)
            
            #index 12 - [Exceed(PMV) SS Int]
            index.12=index.10+index.11
            
            #index 13 - [Exceed(OP) WY]
            index.13=round(sum(rowSums(file.i[7]>file.i[5],na.rm=TRUE))+
                                 sum(rowSums(file.i[7]<file.i[4],na.rm=TRUE)),digits=0)
            
            #index 14 - [Exceed(OP) SS warm]
            index.14=round(sum(rowSums(file.i[7]*id.summer>file.i[5],na.rm=TRUE)),digits=0)
            
            #index 15 - [Exceed(OP) SS cold]
            index.15.temp=file.i[7]; index.15.temp[[1]][id.winter==0]=NA
            index.15=round(sum(rowSums(index.15.temp<file.i[4],na.rm=TRUE)),digits=0)
            
            #index 16 - [Exceed(OP) SS Int]
            index.16=index.14+index.15
            
            #index 17 - [AveragePPD WY]
            index.17=round(as.numeric(apply(file.i[2],FUN=mean,MARGIN=2)),1)
            
            #index 18 - [AveragePPD SS warm]
            index.18.temp=file.i[2]
            index.18.temp[[1]][id.summer==0]=NA
            index.18=round(as.numeric(apply(index.18.temp,FUN=mean,MARGIN=2,na.rm=TRUE)),1)
            
            #index 19 - [AveragePPD SS cold]
            index.19.temp=file.i[2]
            index.19.temp[[1]][id.winter==0]=NA
            index.19=round(as.numeric(apply(index.19.temp,FUN=mean,MARGIN=2,na.rm=TRUE)),1)
            
            #index 20 - [SumPPD WY]
            index.20=round(as.numeric(apply(file.i[2],FUN=sum,MARGIN=2,na.rm=TRUE)),1)
            
            #index 21 - [SumPPD SS warm]
            index.21.temp=file.i[2]
            index.21.temp[[1]][id.summer==0]=NA
            index.21=round(as.numeric(apply(index.21.temp,FUN=sum,MARGIN=2,na.rm=TRUE)),1)
            
            #index 22 - [SumPPD SS cold]
            index.22.temp=file.i[2]
            index.22.temp[[1]][id.winter==0]=NA
            index.22=round(as.numeric(apply(index.22.temp,FUN=sum,MARGIN=2,na.rm=TRUE)),1)
            
            #index 23 - [DhC WY Cooling] with annual average
            TOsup.annual=as.numeric(apply(file.i[5],FUN=mean,MARGIN=2,na.rm=TRUE))
            index.23.count=sum(rowSums(file.i[7]>TOsup.annual,na.rm=TRUE))
            index.23.sum=sum(file.i[file.i[7]>TOsup.annual,7])
            index.23=round(index.23.sum-index.23.count*TOsup.annual,0)
            
            #index 24 - [DhC WY Heating] with annual average
            TOinf.annual=as.numeric(apply(file.i[4],FUN=mean,MARGIN=2,na.rm=TRUE))
            index.24.count=sum(rowSums(file.i[7]<TOinf.annual,na.rm=TRUE))
            index.24.sum=sum(file.i[file.i[7]<TOinf.annual,7])
            index.24=round(-index.24.sum+index.24.count*TOinf.annual,0)
            
            #index 25 - [DhC WY Int] with annual average
            index.25=index.23+index.24
            
            #index 26 - [DhC SS warm Cooling] with annual average of warm period
            TOsup.summer.temp=file.i[5]
            TOsup.summer.temp[[1]][id.summer==0]=NA
            TOsup.summer=as.numeric(apply(TOsup.summer.temp,FUN=mean,MARGIN=2,na.rm=TRUE))
            index.26.temp=file.i[7]
            index.26.temp[[1]][id.summer==0]=NA
            index.26.count=sum(rowSums(index.26.temp>TOsup.summer,na.rm=TRUE))
            index.26.sum=sum(index.26.temp[index.26.temp>TOsup.summer],na.rm=TRUE)
            index.26=round(index.26.sum-index.26.count*TOsup.summer,0)
            
            #index 27 - [DhC SS cold Heating] with annual average of cold period
            TOinf.winter.temp=file.i[4]
            TOinf.winter.temp[[1]][id.winter==0]=NA
            TOinf.winter=as.numeric(apply(TOinf.winter.temp,FUN=mean,MARGIN=2,na.rm=TRUE))
            index.27.temp=file.i[7]
            index.27.temp[[1]][id.winter==0]=NA
            index.27.count=sum(rowSums(index.27.temp<TOinf.winter,na.rm=TRUE))
            index.27.sum=sum(index.27.temp[index.27.temp<TOinf.winter],na.rm=TRUE)
            index.27=round(-index.27.sum+index.27.count*TOinf.winter,0)
            
            #index 28 - [DhC SS Int]
            index.28=index.26+index.27
            
            #index 29 - [PPDwC SS warm]
            index.29.temp=file.i[1]
            index.29.temp[[1]][id.summer==0]=NA
            index.29=round(sum(rowSums(index.29.temp>0.5,na.rm=TRUE)*file.i[2]/10),1)
            
            #index 30 - [PPDwC SS cold]
            index.30.temp=file.i[1]
            index.30.temp[[1]][id.winter==0]=NA
            index.30=round(sum(rowSums(index.30.temp<(-0.5),na.rm=TRUE)*file.i[2]/10),1)
            
            #index 31 - [PPDwC SS Int]
            index.31=index.29+index.30
            
            #index 32 - [OPwC SS warm]
            TOinf.summer.temp=file.i[4]
            TOinf.summer.temp[[1]][id.summer==0]=NA
            TOinf.summer=as.numeric(apply(TOinf.summer.temp,FUN=mean,MARGIN=2,na.rm=TRUE))
            index.32.temp=file.i[7]
            index.32.temp[[1]][id.summer==0]=NA
            index.32.count=rowSums(index.32.temp>TOsup.summer,na.rm=FALSE)*
                  (1+(abs(index.32.temp-TOsup.summer))/abs((TOsup.summer+TOinf.summer)/2-TOsup.summer))
            index.32=as.numeric(apply(index.32.count,FUN=sum,MARGIN=2,na.rm=TRUE))
            
            #index 33 - [OPwC SS cold]
            TOsup.winter.temp=file.i[5]
            TOsup.winter.temp[[1]][id.winter==0]=NA
            TOsup.winter=as.numeric(apply(TOsup.winter.temp,FUN=mean,MARGIN=2,na.rm=TRUE))
            index.33.temp=file.i[7]
            index.33.temp[[1]][id.winter==0]=NA
            index.33.count=rowSums(index.33.temp<TOinf.winter,na.rm=FALSE)*
                  (1+(abs(index.33.temp-TOinf.winter))/abs((TOsup.winter+TOinf.winter)/2-TOinf.winter))
            index.33=as.numeric(apply(index.33.count,FUN=sum,MARGIN=2,na.rm=TRUE))
            
            #index 34 - [OPwC SS Int]
            index.34=index.32+index.33
            
            #index 35 - [DhCwPPD SS warm]
            index.35.temp=file.i[7]
            index.35.temp[[1]][id.summer==0]=NA
            index.35.count=rowSums(index.35.temp>file.i[5],na.rm=TRUE)*
                  (index.35.temp-file.i[5])*file.i[2]/10
            index.35=round(as.numeric(apply(index.35.count,FUN=sum,MARGIN=2,na.rm=TRUE)),0)
            
            #index 36 - [DhCwPPD SS cold]
            index.36.temp=file.i[7]
            index.36.temp[[1]][id.winter==0]=NA
            index.36.count=rowSums(index.36.temp<file.i[4],na.rm=TRUE)*
                  (-index.36.temp+file.i[4])*file.i[2]/10
            index.36=round(as.numeric(apply(index.36.count,FUN=sum,MARGIN=2,na.rm=TRUE)),0)
            
            #index 37 - [DhCwPPD SS Int]
            index.37=index.35+index.36
            
            #column 38 - People internal loads (average)
            column.38=apply(file.i[6],FUN=mean,MARGIN=2,na.rm=TRUE)
            
            index=data.frame(index.1,index.2,index.3,index.4,index.5,index.6,index.7,index.8,
                             index.9,index.10,index.11,index.12,index.13,index.14,index.15,index.16,
                             index.17,index.18,index.19,index.20,index.21,index.22,index.23,index.24,
                             index.25,index.26,index.27,index.28,index.29,index.30,index.31,index.32,
                             index.33,index.34,index.35,index.36,index.37,column.38)
            
            row.names(index)=c(name.short[i])
            names(index)=c(names.columns)
            
            #---results of long-term indices
            newline=data.frame(index)
            res=rbind(res,newline)
            
      }
      
      name.output=c(paste(output,"/indices.csv",sep=""))
      write.table(res, file=name.output,row.names=FALSE,col.names=names.columns,sep=",")
      
      #---Vectorial normalization
      j=1
      for(j in 1:37)
      {
            squares=res[j]^2
            index.norm=(res[j]/sqrt(apply(squares,FUN=sum,MARGIN=2,na.rm=TRUE)))
            newcolumn=data.frame(index.norm)
            res.norm=cbind(res.norm,newcolumn)
      }
      
      name.output.norm=c(paste(output,"/indices-norm.csv",sep=""))
      write.table(res.norm,file=name.output.norm,row.names=FALSE,col.names=TRUE,sep=",")
      
      return(list(res,res.norm))
}
