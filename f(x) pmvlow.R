###---Function that calculates PMV and PPD according to low air speed model from ISO 7730 (2005)
###---Input data:
#---CLO clothing insulation [clo]
#---MET metabolic rate [met]
#---WME mechanical work [met]
#---TA air temperature [oC]
#---TR mean radiant temperature [oC]
#---VEL air speed [m/s], should be lower or equal to 0.10 m/s
#---RH relative humidity [%]
#---PA partial water vapor pressure [kPa]

#---Use this for testing the FUNCTION pmvlow:
#CLO=0.5; MET=1.0; WME=0; TA=23; TR=24; VEL=0.05; RH=50
#pmvlow(CLO,MET,WME,TA,TR,VEL,RH)

pmvlow=function(CLO,MET,WME,TA,TR,VEL,RH,PA=0)
{
  #---Enter only the relative humidity, not the water vapor pressure PA
  PA=0
  
  #---Initial Calculation
  FNPS.T = exp(16.6536-4030.183/(TA+235))    #saturated vapor pressure KPa
  if(PA==0) PA=RH*10*FNPS.T                  #water vapor pressure, Pa
  ICL=0.155*CLO                              #thermal insulation of the clothing in m2K/W
  M=MET*58.15                                #metabolic rate in W/m2
  W=WME*58.15                                #external work in W/m2
  MW=(M-W)                                   #internal heat production in the human body
  
  FCL=NULL
  FCL=ifelse(ICL<=0.078,1+1.29*ICL,1.05+0.645*ICL)     #clothing area factor
  
  #---heat transfer coefficient by forced convection
  hcf=12.1*sqrt(VEL)
  taa=TA+273
  tra=TR+273
  tcla=taa+(35.5-TA)/(3.5*ICL+0.1)
  
  p1=ICL*FCL
  p2=p1*3.96
  p3=p1*100
  p4=p1*taa
  p5=(308.7-0.028*MW)+(p2*(tra/100)^4)
  xn=tcla/100
  xf=tcla/50
  eps=0.00015
  
  n=0
  while (abs(xn-xf)>eps)
  {
    xf=(xf+xn)/2
    hcn=2.38*(abs(100*xf-taa)^0.25)
    if (hcf>hcn) hc=hcf else hc=hcn
    
    xn=(p5+p4*hc-p2*(xf^4))/(100+p3*hc)
    n=n+1
    
    #---If number of iterations exceeds 150
    if (n>150)
    {
      print("Maximum number of iterations exceeded in PMV low air speed model, reaching iterations n==150...")
      break;
    }
  }
      
  tcl=100*xn-273
  
  #---heat loss difference through skin 
  hl1=3.05*0.001*(5733-(6.99*MW)-PA)
  #---heat loss by sweating
  if (MW>58.15) hl2=0.42*(MW-58.15) else hl2=0
  #---latent respiration heat loss 
  hl3=1.7*0.00001*M*(5867-PA)
  #---dry respiration heat loss
  hl4=0.0014*M*(34-TA)
  #---heat loss by radiation  
  hl5=3.96*FCL*(xn^4-(tra/100)^4)
  #---heat loss by convection
  hl6=FCL*hc*(tcl-TA)
  
  ts=0.303*exp(-0.036*M)+0.028
  
  #---Predicted Mean Vote - PMV
  PMV=ts*(MW-hl1-hl2-hl3-hl4-hl5-hl6)
  
  #------Predicted Percentage Dissatisfied - PPD
  PPD=100-95*exp(-0.03353*PMV^4-0.2179*PMV^2)
  
  result=c(round(PMV,2),round(PPD,0))
  names(result)=c("PMV","PPD")
  
  return(result)
}
