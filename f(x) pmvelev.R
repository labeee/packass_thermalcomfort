###---Function that calculates PMV and PPD according to elevated air speed model from ASHRAE Standard 55 (2017)
###---Input data:
#---CLO clothing insulation [clo]
#---MET metabolic rate [met]
#---WME mechanical work [met]
#---TA air temperature [oC]
#---TR mean radiant temperature [oC]
#---VEL air speed [m/s]
#---RH relative humidity [%]
#---PA partial water vapor pressure [kPa]
#---APA atmospheric pressure [kPa]

#---Use this for testing the FUNCTION pmvelev:
#CLO=0.35; MET=1.1; WME=0; TA=28; TR=28; VEL=0.35; RH=50; PA=0; APA=101.325
#pmvelev(CLO,MET,WME,TA,TR,VEL,RH)

pmvelev=function(CLO,MET,WME,TA,TR,VEL,RH,PA=0,APA=101.325)
{
  #---Set the working directory to this file path, to ensure the correct linking
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  #---Link to "f(x) setashrae.R" script path - SET calculation
  source("f(x) setashrae.R")
  #---Link to "f(x) pmvlow.R" script path - low air speed model
  source("f(x) pmvlow.R")
  
  #---SET calculated with the actual input
  SET=setashrae(CLO,MET,WME,TA,TR,VEL,RH,APA)
  
  ###---Secant method of iteration
  
  #---Parameters of convergence. Depending of the values, the convergence cannot be achieved properly...
  a.1=0
  b.1=40 #---Java CBE script is 40, I set to 100
  count.1=1
  F1=setashrae(CLO,MET,WME,TA-a.1,TR-a.1,0.10,RH,APA)-SET
  if(abs(F1)<0.001) ce=a.1
  F2=setashrae(CLO,MET,WME,TA-b.1,TR-b.1,0.10,RH,APA)-SET
  if(abs(F2)<0.001) ce=b.1
  
  while (count.1<=101)
  {
    slope=(F2-F1)/(b.1-a.1)
    c.1=b.1-F2/slope
    
    if(c.1>999999 | c.1==Inf | c.1==-Inf)
    {
      ce=NaN
      break;
    }
    
    F3=setashrae(CLO,MET,WME,TA-c.1,TR-c.1,0.10,RH,APA)-SET
    
    #print(paste("F3= ",F3,sep=""))
    
    if(abs(F3)<0.001)
    {
      ce=c.1
      break;
    }
    a.1=b.1
    b.1=c.1
    F1=F2
    F2=F3
    
    count.1=count.1+1
        
    if(count.1==101)
    {
      ce=NaN
      #print("Maximum number of iteration exceeded in secant method of PMV elevated air speed model...")
      break;
    }
  }
  
  #---Bisection method of iteration
  a.2=0 #---Java CBE script is 0, I set to -84
  b.2=80 #---Java CBE script is 40, I set to 94
  DIFF=0.5
  count.2=1
  if(ce=='NaN')
  {    
    while (abs(DIFF)>0.002)
    { 
      c.2=(a.2+b.2)/2
      Eq.1=setashrae(CLO,MET,WME,TA-a.2,TR-a.2,0.10,RH,APA)-SET
      Eq.2=setashrae(CLO,MET,WME,TA-b.2,TR-b.2,0.10,RH,APA)-SET
      Eq.3=setashrae(CLO,MET,WME,TA-c.2,TR-c.2,0.10,RH,APA)-SET
      
      if (Eq.1*Eq.3<0) b.2=c.2
      if (Eq.2*Eq.3<0) a.2=c.2 #---changed... it was Eq.1*Eq.3>0
      if (abs(Eq.3)<0.001) break;
      
      #print(paste("Eq3= ",Eq3,sep=""))
      
      DIFF=abs(b.2-a.2)
      ce=c.2
      count.2=count.2+1
            
      if(count.2==101)
      {
        #print("Maximum number of iteration exceeded in bisection method of PMV elevated air speed model...")
        break;
      }
      
    }
  } else { 
    ce=c.1
  }
  
  #---Predicted Mean Vote - PMV - call function pmvlow()
  PMV=pmvlow(CLO,MET,WME,TA-ce,TR-ce,0.10,RH,PA)[[1]]
  
  #---Predicted Percentage Dissatisfied - PPD - call function pmvlow()
  PPD=100-95*exp(-0.03353*PMV^4-0.2179*PMV^2)
  
  result=c(round(PMV,2),round(PPD,0),round(ce,1))
  names(result)=c("PMV","PPD","CE")
  
  return(result)
}