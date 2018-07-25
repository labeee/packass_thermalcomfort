#---Function that calculate PMV and PPD according to elevated air speed model from ASHRAE Standard 55 (2013)
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

#---Use this for testing the FUNCTION pmvboth:
#CLO=0.5; MET=1.0; WME=0; TA=23; TR=24; VEL=0.40; RH=50; PA=0; APA=101.325
#(pmvboth(CLO,MET,WME,TA,TR,VEL,RH,PA,APA))

pmvboth=function(CLO,MET,WME,TA,TR,VEL,RH,PA=0,APA=101.325)
{
  #---Set the working directory to this file path, to ensure the correct linking
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  #---Link to "f(x) pmvlow.R" script path - low air speed model
  source("f(x) pmvlow.R")
  #---Link to "f(x) pmvelevated.R" script path - elevated air speed model
  source("f(x) pmvelev.R")
  
  if (VEL<=0.10)
  {
    PMV=round(pmvlow(CLO,MET,WME,TA,TR,VEL,RH,PA)[[1]],2)
    CE=0
  }
  
  if (VEL>0.10)
  {
    PMV.elev=pmvelev(CLO,MET,WME,TA,TR,VEL,RH,PA,APA)
    PMV=round(PMV.elev[[1]],2)
    CE=round(PMV.elev[[3]],1)
  }

  resultado=PMV
  
  #---Predicted Percentage Dissatisfied - PPD
  PPD=round(100-95*exp(-0.03353*PMV^4-0.2179*PMV^2),0)
  
  result=c(PMV,PPD,CE)
  names(result)=c("PMV","PPD","CE")
  
  return(result)
}