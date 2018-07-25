#---Function that calculates upper and lower operative temperatures according to ASHRAE Standard 55 (2013)
#---low and elevated air speed models.
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

#---Use to teste the function:
#CLO=1.3; MET=1.3; WME=0; VEL=1.2; RH=80; PA=0; APA=101.325
#optemplimits(CLO,MET,WME,VEL,RH)

optemplimits=function (CLO,MET,WME,VEL,RH,PA=0,APA=101.325)
{
  #---Set the working directory to this file path, to ensure the correct linking
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  
  #---Link to "f(x) PMV_low air speed model.R" script path
  source("f(x) pmvboth.R")
  
  ###---Inferior/lower operative temperature limit calculation
 
  #---Secant method
  PMV.inf=-0.5
  a.1=-10 #---mudei de -100 para -50 em 19/06/2017
  b.1=+60 #---It was 114
  count.1=1
  
  F1=pmvboth(CLO,MET,WME,a.1,a.1,VEL,RH,PA,APA)[[1]]-PMV.inf
  if(abs(F1)<0.01) TO.inf=a.1
  F2=pmvboth(CLO,MET,WME,b.1,b.1,VEL,RH,PA,APA)[[1]]-PMV.inf
  if(abs(F2)<0.01) TO.inf=b.1
  
  while (count.1<=151)
  {
    slope=(F2-F1)/(b.1-a.1)
    c.1=b.1-F2/slope
    
    if(c.1>999999 | c.1==Inf | c.1==-Inf)
    {
      TO.inf=NaN
      #print(paste("Exceeded number of iteration in secant method > TO==Inf",count.1,sep=""))
      break;
    }
    
    F3=pmvboth(CLO,MET,WME,c.1,c.1,VEL,RH,PA,APA)[[1]]-PMV.inf
    #print(paste("F3= ",F3,sep=""))
    
    if(abs(F3)<0.01)
    {
      TO.inf=c.1
      break;
    }
    a.1=b.1
    b.1=c.1
    F1=F2
    F2=F3
    count.1=count.1+1
    
    if(count.1==101)
    {
      TO.inf=NaN
      #print("Exceeded maximum number of iteration in secant method")
      break;
    }
  }
  
  #---Bisection method
  a.2=-10
  b.2=60
  DIFF=0.5
  count.2=1
  if(TO.inf=='NaN')
  {    
    while (abs(DIFF)>0.02)
    { 
      c.2=(a.2+b.2)/2
      Eq.1=pmvboth(CLO,MET,WME,a.2,a.2,VEL,RH,PA,APA)[[1]]-PMV.inf
      Eq.2=pmvboth(CLO,MET,WME,b.2,b.2,VEL,RH,PA,APA)[[1]]-PMV.inf
      Eq.3=pmvboth(CLO,MET,WME,c.2,c.2,VEL,RH,PA,APA)[[1]]-PMV.inf
      #print(paste("Eq.3= ",Eq.3,sep=""))
      
      if (Eq.1*Eq.3<0) b.2=c.2
      if (Eq.2*Eq.3<0) a.2=c.2 #---Changed. It was Eq.1*Eq.3>0
      if (abs(Eq.3)<0.0001) break;
      
      DIFF=abs(b.2-a.2)
      TO.inf=c.2
      count.2=count.2+1
      
      if(count.2==101)
      {
        #print("Exceeded maximum number of iteration in bisection method")
        break;
      }
    }
  } else { 
    TO.inf=c.1
  }
  
  ###---Inferior limit calculation
  
  #---Secant method
  PMV.sup=0.5
  a.1=-10 #---mudei de -100 para -50 em 19/06/2017
  b.1=60
  count.1=1
  F1=pmvboth(CLO,MET,WME,a.1,a.1,VEL,RH,PA,APA)[[1]]-PMV.sup
  if(abs(F1)<0.01) TO.sup=a.1
  F2=pmvboth(CLO,MET,WME,b.1,b.1,VEL,RH,PA,APA)[[1]]-PMV.sup
  if(abs(F2)<0.01) TO.sup=b.1
  
  while (count.1<=151)
  {
    slope=(F2-F1)/(b.1-a.1)
    c.1=b.1-F2/slope
    
    if(c.1>999999 | c.1==Inf | c.1==-Inf)
    {
      TO.sup=NaN
      #print(paste("Exceeded number of iteration in secant method > TO==Inf",count.1,sep=""))
      break;
    }
    
    F3=pmvboth(CLO,MET,WME,c.1,c.1,VEL,RH,PA,APA)[[1]]-PMV.sup
    #print(paste("F3= ",F3,sep=""))
    
    if(abs(F3)<0.01)
    {
      TO.sup=c.1
      break;
    }
    a.1=b.1
    b.1=c.1
    F1=F2
    F2=F3
    count.1=count.1+1
    
    if(count.1==101)
    {
      TO.sup=NaN
      #print("Exceeded maximum number of iteration in secant method")
      break;
    }
  }
  
  #---Bisection method
  a.2=-10 #---It was -90
  b.2=60 #---It was 95
  DIFF=0.5
  count.2=1
  if(TO.sup=='NaN')
  {    
    while (abs(DIFF)>0.02)
    { 
      c.2=(a.2+b.2)/2
      Eq.1=pmvboth(CLO,MET,WME,a.2,a.2,VEL,RH,PA,APA)[[1]]-PMV.sup
      Eq.2=pmvboth(CLO,MET,WME,b.2,b.2,VEL,RH,PA,APA)[[1]]-PMV.sup
      Eq.3=pmvboth(CLO,MET,WME,c.2,c.2,VEL,RH,PA,APA)[[1]]-PMV.sup
      #print(paste("Eq.3= ",Eq.3,sep=""))
      
      if (Eq.1*Eq.3<0) b.2=c.2
      if (Eq.1*Eq.3>0) a.2=c.2
      if (abs(Eq.3)<0.0001) break; #It was 0.0001
      
      DIFF=abs(b.2-a.2)
      TO.sup=c.2
      count.2=count.2+1
      
      if(count.2==101)
      {
        #print("Exceeded maximum number of iteration in bisection method")
        break;
      }
    }
  } else { 
    TO.sup=c.1
  }
  
  result=c(round(TO.inf,1),round(TO.sup,1))
  names(result)=c("TO inferior","TO superior")
  
  return(result)
  
}