###---Function that calculates SET according to ASHRAE Standard 55 (2013)
###---Input data:
#---CLO clothing insulation [clo]
#---MET metabolic rate [met]
#---WME mechanical work [met]
#---TA air temperature [oC]
#---TR mean radiant temperature [oC]
#---VEL air speed [m/s]
#---RH relative humidity [%]
#---APA atmospheric pressure [kPa]

#---Use this for testing the FUNCTION setashrae:
#CLO=0.5; MET=1.0; WME=0; TA=23; TR=24; VEL=0.40; RH=50
#setashrae(CLO,MET,WME,TA,TR,VEL,RH)

setashrae=function(CLO,MET,WME,TA,TR,VEL,RH,APA=101.325)
{

  #APA=101.325   #---Default of 1 atm equal to 101.325 kPa
  
  #---vapor pressure function in [Torr]
  FindSVP=function(Temperature)
  {
  resultado.SVP=exp(18.6686-4030.183/(Temperature+235))
  return (resultado.SVP)
  }

  #---Constants
  BodyWeight=69.9
  BodySurfaceArea=1.8258
  MetFactor=58.2
  StefanB=5.6697*10^(-8) #---constante de Stefan Boltzmann [W/m2K4]
  TTCR=36.8             #---neutral core body temperature in Python was 36.49, in Java is 36.8
  TTSK=33.7              #---neutral skin body temperature
  TTBM=36.49             #---neutral body temperature
  SKBFN=6.3              #---skin blood flow neutral
  PressureInAtmospheres=APA*0.009869 #---presure from kPa to atm

  #---The actual SET code from comfort tool uses WME in W/m²...so it wil be converted
  WME=WME*MetFactor
  
  #---Initial values
  VaporPressure=RH*FindSVP(TA)/100   
  AirVelocity=max(VEL,0.1) #---0.1 m/s at least
  M=MET*MetFactor      
  RM=MET*MetFactor
  RCL=0.155*CLO            #---m2K/W
  FACL=1+0.15*CLO      
  
  LR=2.2/PressureInAtmospheres
  alpha=0.1                #---actual skin to total body mass ratio
  MSHIV=0
  ESK=0.1*MET

  #---Simulation initial values
  TSK=TTSK   #---Temperature Skin
  TCR=TTCR   #---Temperature Core
  SKBF=SKBFN #---Skin Blood Flow
  
  #---Physiological temperatura regulation controls
  CDIL=120   #----liters/(m2.h.K)
  CSTR=0.5   #----dimensionless
  CSW=170    #----g/(m2.h)
  SKBFL=90   #----liter/(m2h) max SKBF
  REGSWL=500
  
  if (CLO<=0) {
    WCRIT=0.38*AirVelocity^(-0.29)     
    ICL=1.0
  } else {
    WCRIT=0.59*AirVelocity^(-0.08)
    ICL=0.45
  }
  
  CHC.t=3*(PressureInAtmospheres)^0.53
  CHCV=8.600001*(AirVelocity*PressureInAtmospheres)^0.53
  CHC=max(CHC.t,CHCV)
  
  #---initial estimate of TCL
  CHR=4.7              
  CTC=CHR+CHC          
  RA=1/(FACL*CTC)   
  TOP=(CHR*TR+CHC*TA)/CTC   
  TCL=TOP+(TSK-TOP)/(CTC*(RA+RCL))
  
  ##########################---Begin 60 minutes simulation---##########################
  LTIME=60   #---Timestep equal to 60
  TIM=1
  TCL.old=TCL
  flag=TRUE
  while (TIM<=LTIME)
  {
    #---Dry heat balance - Solve TCL and CHR
    DIFF.1=0.5
    if(flag==TRUE)
    {
      count.1=0
      while (DIFF.1>0.01)
      {
        TCL.old=TCL
        CHR=4*0.72*StefanB*((TCL+TR)/2+273.15)^3
        CTC=CHR+CHC
        RA=1/(FACL*CTC)
        TOP=(CHR*TR+CHC*TA)/CTC
        DIFF.1=abs(TCL.old-TCL)
        TCL=(RA*TSK+RCL*TOP)/(RA+RCL)
        count.1=count.1+1
        #print(paste("count.1 inside TCL iteration= ",count.1,sep=""))
        
        #---Iteration from bisection method
        if(DIFF.1>99999) #---an alternative control from non-convergente (rarely used)
        {
           print(paste("entered in the bissection method TIME= ",count.1,sep=""))
           DIFF.5=0.5
           o.1=-300
           o.2=+350
           count.2=0
           while(abs(DIFF.5)>0.01)
           {
             o.3=(o.1+o.2)/2
          
             TCL=o.1
             CHR=4*0.72*StefanB*((TCL+TR)/2+273.15)^3
             CTC=CHR+CHC
             RA=1/(FACL*CTC)
             TOP=(CHR*TR+CHC*TA)/CTC
             fun.o1=(RA*TSK+RCL*TOP)/(RA+RCL)-TCL
          
             TCL=o.2
             CHR=4*0.72*StefanB*((TCL+TR)/2+273.15)^3
             CTC=CHR+CHC
             RA=1/(FACL*CTC)
             TOP=(CHR*TR+CHC*TA)/CTC
             fun.o2=(RA*TSK+RCL*TOP)/(RA+RCL)-TCL
          
             TCL=o.3
             CHR=4*0.72*StefanB*((TCL+TR)/2+273.15)^3
             CTC=CHR+CHC
             RA=1/(FACL*CTC)
             TOP=(CHR*TR+CHC*TA)/CTC
             fun.o3=(RA*TSK+RCL*TOP)/(RA+RCL)-TCL
          
             if(fun.o1*fun.o3<0) o.2=o.3
             if(fun.o1*fun.o3>0) o.1=o.3
             if(abs(fun.o3)<0.0001) break;
          
             DIFF.5=abs(o.2-o.1)
             TCL=o.3
             count.2=count.2+1
             
             if(count.2==101)
             {
              print("Maximum number of iteration exceeded in TCL inside SET from bisection method...")
              break;
             } 
           }
         }#---End of control of non-cnovergence
      }
    }
    flag=FALSE
    #---Heat flow from clothing surface to environment (FACL=1 if CLOE used)
    DRY=(TSK-TOP)/(RA+RCL)                
    HFCS=(TCR-TSK)*(5.28+1.163*SKBF)      

    #---Dry and latent respiratory heat losses
    ERES=0.0023*M*(44-VaporPressure)     
    CRES=0.0014*M*(34-TA)                

    SCR=M-HFCS-ERES-CRES-WME
    SSK=HFCS-DRY-ESK
  
    #---Thermal capacities
    TCSK=0.97*alpha*BodyWeight
    TCCR=0.97*(1-alpha)*BodyWeight
      
    #---Temperature changes in 1 minute
    DTSK=(SSK*BodySurfaceArea)/(TCSK*60)   
    DTCR=(SCR*BodySurfaceArea)/(TCCR*60)   
    TSK=TSK+DTSK
    TCR=TCR+DTCR
    TB=alpha*TSK+(1-alpha)*TCR

    #---Definition of vascular control signals (Old implementation from Gagge et al...)
    #if (TSK>TTSK)
    #{
    #  WARMS=TSK-TTSK
    #  COLDS=0
    #} else {
    #  COLDS=TTSK-TSK
    #  WARMS=0
    #}
    #if (TCR>TTCR)
    #{
    #  WARMC=TCR-TTCR
    #  COLDC=0
    #} else {
    #  COLDC=TTCR-TCR
    #  WARMC=0
    #}
    #if (TB>TTBM)
    #{
    #  WARMB=TB-TTBM
    #  COLDB=0
    #} else {
    #  COLDB=TTBM-TB
    #  WARMB=0
    #}
   
    #---Definition of vascular control signals (New implementation from Comfort Tool Python code...)
    SKSIG=TSK-TTSK
    WARMS=(SKSIG>0)*SKSIG
    COLDS=((-1.0*SKSIG)>0)*(-1.0*SKSIG)
    CRSIG=TCR-TTCR
    WARMC=(CRSIG>0)*CRSIG
    COLDC=((-1.0*CRSIG)>0)*(-1.0*CRSIG)
    BDSIG=TB-TTBM
    WARMB=(BDSIG>0)*BDSIG
    COLDB=((-1.0*BDSIG)>0)*(-1.0*BDSIG)
    
    #---Control skin blood flow
    SKBF=((SKBFN+CDIL*WARMC)/(1+CSTR*COLDS))
    #Old code, separating some terms
    #DILAT=CDIL*WARMC              
    #STRIC=CSTR*COLDS             
    #SKBF=(SKBFN+DILAT)/(1+STRIC)
    
    #Limits from Comfort Tool... Gagge et al. has different limits.
    if (SKBF>90) SKBF=90
    if (SKBF<0.5) SKBF=0.5
    #Old code, with older limits
    #SKBF=max(0.5,min(SKBFL,SKBF)) #SKBF is never below 0.5 liter/m2h not above SKBFL
    
    #---control of regulatory sweting
    REGSW=CSW*WARMB*exp(WARMS/10.7)
    REGSW=min(REGSW,REGSWL)
    ERSW=0.68*REGSW
      
    #---Mass transfer equation between skin and environment
    REA=1/(LR*FACL*CHC) #evaporative resistance of air layer
    RECL=RCL/(LR*ICL)
    EMAX=(FindSVP(TSK)-VaporPressure)/(REA+RECL)
    PRSW=ERSW/EMAX

    #---PDIF for nonsweating skin
    PWET=0.06+0.94*PRSW
    EDIF=PWET*EMAX-ERSW
    ESK=ERSW+EDIF

    #---Beginning of dripping (swet not evaporated on skin surface)
    if (PWET>WCRIT)
    {
      PWET=WCRIT
      PRSW=WCRIT/0.94          
      ERSW=PRSW*EMAX
      EDIF=0.06*(1-PRSW)*EMAX
      ESK=ERSW+EDIF
    }

    #---When EMAX<0 condensation on skin ocurs
    if (EMAX<0)
    {
      EDIF=0
      ERSW=0
      PWET=WCRIT
      PRSW=WCRIT
      ESK=EMAX
    }
  
    #---Adjustment of metabolic heat due to shivering
    ESK=ERSW+EDIF
    MSHIV=19.4*COLDS*COLDC
    M=RM+MSHIV
    
    #---Ratio of skin-core masses change with SKBF
    alpha=0.0417737+0.7451833/(SKBF+0.585417)
    
    TIM=TIM+1 #---Next minute
    
    #---To verify the function...
    #print(paste("iteration inside 60min sim TIM= ",TIM,sep=""))
    #print(paste("ESK= ",ESK,sep=""))
  }
  ##########################---End 60 minutes simulation---##########################
  
  STORE=M-WME-CRES-ERES-DRY-ESK #This is not used...
  HSK=DRY+ESK
  RN=M-WME
  ECOMF=0.42*(RN-(1*MetFactor))
  if (ECOMF<0) ECOMF=0 #from Fanger?
  EMAX=EMAX*WCRIT
  W=PWET
  PSSK=FindSVP(TSK)
  
  #---Standard environment
  CHRS=CHR
  
  #---CHCS=standard convective heat transfer coefficient (level walking/still air)
  if (MET<0.85) {
    CHCS=3
  } else {
    CHCS.t=5.66*(MET-0.85)^0.39
    CHCS=max(3,CHCS.t)
  }

  KCLO=0.25
  
  CTCS=CHCS+CHRS
  RCLOS=1.52/((MET-WME/MetFactor)+0.6944)-0.1835
  RCLS=0.155*RCLOS
  FACLS=1+KCLO*RCLOS
  FCLS=1/(1+0.155*FACLS*CTCS*RCLOS)
  IMS=0.45
  ICLS=IMS*CHCS/CTCS*(1-FCLS)/(CHCS/CTCS-FCLS*IMS)
  RAS=1/(FACLS*CTCS)
  REAS=1/(LR*FACLS*CHCS)
  RECLS=RCLS/(LR*ICLS)
  HD.S=1/(RAS+RCLS)
  HE.S=1/(REAS+RECLS)
  
  ##########################---Standard Effective Temperature SET---##########################
  
  #---Newton method of iteration
  a.3=TSK-HSK/HD.S #---Initial estimation of SET
  DIFF.3=100
  count.3=0
  delta=0.0001
  SET.OLD=a.3
  while (abs(DIFF.3)>0.001)
  { 
    ER.1=HSK-HD.S*(TSK-SET.OLD)-W*HE.S*(PSSK-0.5*FindSVP(SET.OLD))
    ER.2=HSK-HD.S*(TSK-SET.OLD-delta)-W*HE.S*(PSSK-0.5*FindSVP(SET.OLD+delta))
    SET=SET.OLD-delta*ER.1/(ER.2-ER.1)
    DIFF.3=SET-SET.OLD
    SET.OLD=SET
    count.3=count.3+1
    
    if(count.3==101)
    {
      print("Maximum number of iteration exceeded in SET iteration for Newton method")
      break;
    }
    
    #---Bissection method implementation
    if(abs(DIFF.3)>999999)
    {
      SET=NULL
      SET.OLD=NULL
      DIFF.2=0.5
      count.4=0
      a2=-301
      b2=486
      while (DIFF.2>0.01)
      {
        c2=(b2+a2)/2
        SET.OLD=a2
        ERR1=HSK-HD.S*(TSK-SET.OLD)-W*HE.S*(PSSK-0.5*FindSVP(SET.OLD))
        SET.OLD=b2
        ERR2=HSK-HD.S*(TSK-SET.OLD)-W*HE.S*(PSSK-0.5*FindSVP(SET.OLD))
        SET.OLD=c2
        ERR3=HSK-HD.S*(TSK-SET.OLD)-W*HE.S*(PSSK-0.5*FindSVP(SET.OLD))
        if (ERR1*ERR3<0) b2=c2
        if (ERR1*ERR3>0) a2=c2
        if (abs(ERR3)<0.0001) break;
        DIFF.2=abs(b2-a2)
        SET.OLD=c2
        count.4=count.4+1
        
        if(count.4==101)
        {
          #print("Maximum number of iteration exceeded in SET iteration for bisection method")
          break;
        }
        
        SET=SET.OLD #---iterative calculation
      }
    }
  }
    
  result=round(SET,1)
  return(result)

}
