#
#
## Code for multievent model applied to Cook Inlet beluga whale population
## As described in Himes Boor, GK, TL McGuire, AJ Warlick, RL Taylor, SJ Converse, JR McClung, AD Stephens: 
## Estimating a reproductive rate when offspring ages are uncertain: a novel multievent mark-recapture model applied to an 
##   endangered beluga whale population
##
## This script loads data and runs the model described in the manuscript. 
## More details about the model structure and data can be found in the manuscript and supplemental information.
##
##
## 2021-11-03
##
## Although these data have been processed successfully on a computer system at the U.S. Geological Survey (USGS), no warranty expressed or implied is made regarding the display or utility of the data for other purposes, nor on all computer systems, nor shall the act of distribution constitute any such warranty. The USGS or the U.S. Government shall not be held liable for improper or incorrect use of the data described and/or contained herein.


##Required libraries
library(readr)
library(jagsUI)  #for easier running of parallel chains
library(R2OpenBUGS)
library(rjags)



##########################################
## LOAD AND PREP DATA
##########################################


## Choose which data to be used
data.file<-"ms_SH_data"   ##left-side CIBW data


# Read in the data: 
mydata<-read_csv(file=paste(data.file,".csv",sep=""))
# Read in starting latent state matrix
start.mat.csv<-read_csv(file=paste("start_mat-",data.file,".csv",sep=""))

# Get data dimensions (N=# individuals, K=# of observation periods (years))
N <- dim(mydata)[1]
K <- dim(mydata)[2]

# Compute the year of first capture for each individual:
get.first.cap <- function(x) min(which(x != "1"))
e <- apply(mydata, 1, get.first.cap)

alive.function <- function(){
  alive <- start.mat.csv
  for (i in 1:N){
       for (j in 1:K) {
            if (j < e[i]) {alive[i,j] = NA}
       }
  }
  alive
}


##########################################
## DEFINE & RUN MODEL USING JAGS
##########################################

# Define the model

#Notation:
# OBSERVATIONS (i.e. Events; n=72):
# 1	=	0	not seen     
# 2	=	P	seen alone     
# 3	=	J0	seen with YOY     
# 4	=	J1-	seen with a calf of indeterminate age (YOY or 1yo)     
# 5	=	J1	seen with 1yo calf (i.e, 1st birthday to 2nd birthday)     
# 6	=	J1+	seen with a calf of indeterminate age (1yo or older)     
# 7	=	J2-	seen with a calf of indeterminate age (YOY, 1yo, or 2yo)     
# 8	=	J2	seen with 2yo calf     
# 9	=	J2+	seen with a calf of indeterminate age (2yo or older)     
# 10	=	J3	seen with 3yo calf     
# 11	=	J3+	seen with a calf of indeterminate age (3yo or older)     
# 12	=	J4	seen with a 4yo calf     
# 13	=	J4+	seen with a calf of indeterminate age (4yo or older)     
# 14	=	Unk	seen with a calf of indeterminate age        
# 15	=	J1+ J0	seen with 2 calves, one a YOY (J0) and one 1yo or older (except must be 2yo or older)   
# 16	=	J1+ J1-	seen with 2 calves, one 1yo or older and one 1yo or younger
# 17	=	J1+ J1	et
# 18	=	J2- J0	
# 19	=	J2- J1-	
# 20	=	J2 J0	
# 21	=	J2 J1-	
# 22	=	J2 J2-	
# 23	=	J2+ J0	
# 24	=	J2+ J1-	
# 25	=	J2+ J1	
# 26	=	J2+ J1+	
# 27	=	J2+ J2-	
# 28	=	J2+ J2	
# 29	=	J3 J0	
# 30	=	J3 J1-	
# 31	=	J3 J1	
# 32	=	J3 J1+	
# 33	=	J3 J2-	
# 34	=	J3+ J0	
# 35	=	J3+ J1-	
# 36	=	J3+ J1	
# 37	=	J3+ J1+	
# 38	=	J3+ J2-	
# 39	=	J3+ J2	
# 40	=	J3+ J2+	
# 41	=	J4 J0	
# 42	=	J4 J1-	
# 43	=	J4 J1	
# 44	=	J4 J1+	
# 45	=	J4 J2-	
# 46	=	J4 J2	
# 47	=	J4 J2+	
# 48	=	J4+ J0	
# 49	=	J4+ J1-     	
# 50	=	J4+ J1     	
# 51	=	J4+ J1+	     
# 52	=	J4+ J2-	     
# 53	=	J4+ J2	     
# 54	=	J4+ J2+	     
# 55	=	J0 Unk	     
# 56	=	J1- Unk	     
# 57	=	J1 Unk	     
# 58	=	J1+ Unk	     
# 59	=	J2- Unk	     
# 60	=	J2 Unk	     
# 61	=	J2+ Unk	     
# 62	=	J3 Unk	     
# 63	=	J3+ Unk	     
# 64	=	J4 Unk	     
# 65	=	J4+ Unk	     
# 66	=	Unk Unk
# 67 =    J1+ J1+
# 68 =    J1+ J2-
# 69 =    J1+ J2
# 70 =    J1+ J2+
# 71 =    J2- J2-
# 72 =    J2+ J2+  
     

# TRUE STATES
# 1  = NB = alive non-breeder     
# 2  = B = alive previous breeder w/no calves
# 3  = Byoy = breeder w/YOY               
# 4  = Bc1 = breeder w/1yo calf           
# 5  = Bc2 = breeder w/2yo calf              
# 6  = Bc2yoy = breeder w/YOY & 2yo calf  
# 7  = Bc3 = breeder w/3yo calf              
# 8  = Bc3yoy = breeder w/YOY & 3yo calf    
# 9  = Bc3c1 = breeder w/3yo & 1yo calf   
# 10 = Bc4 = breeder w/4yo calf               
# 11 = Bc4yoy = breeder w/YOY & 4yo calf    
# 12 = Bc4c1 = breeder w/4yo & 1yo calf     
# 13 = Bc4c2 = breeder w/4yo & 2yo calf  
# 14 = D = dead     


# PARAMETERS
# phiN  survival prob. of non-breeding class of adults/subadults
# phiB  survival prob. of breeding class of adults
# phiYc survival prob. of dependent calves (YOYs & 1yos)
# phiC  survival prob. of calves (2yos and older)
# psiN  transition prob. from non-breeder to breeder
# psiB  transition prob. from breeder-no calf to breeder-with calf
# pN    detection prob. of non-breeding class adults/subadults
# pBn   detection prob. of breeding female adults with no calf
# pBc   detection prob. of breeding female adults with a calf of any age
# deltaYc detection prob of YOY & 1yo calves
# deltaC  detectin prob of <=2yo calves
# gamma   prob of calf getting put in unknown age category
# alphaTy prob of a YOY getting categorized as YOY without uncertainty
# alphaTc prob of a <=1yo calf getting categorized as their true age without uncertainty
# kappaY  prob a YOY is assigned to J1- category
# kappaC  prob of a 3 or 4yo being assigned to the J2+ or J3+ categories respectively
# omegaA  prob of a 1 or 2yo being assigned to one of several categories that match true age with uncertainty
# omegaB  prob of a 3 or 4yo being assigned to category equal to true age or older (J3+ or J4+)
# eta     prob of 1 or 2yo being assigned to category equal to true age or younger (J1-,J2-)
# pi[1:13] prob. of being in initial states 1 to 13 (NB to Bc4c2)
# 

#MATRICES:
#State-Transition Matrices
#There are 4 state-transition matrices. I list them below with their definition and dimensions (nrow,ncol).
# Matrix 1: Adult Survival [14,14]
# Matrix 2: Young survival [14,37]
# Matrix 3: Deterministic calf aging [37,14]
# Matrix 4: Adult breeding transition [14,14]

#Observation-Event Matrices
#There are 3 observation matrices. I list them below with their dimensions (nrow,ncol).
# Matrix 1: Adult detection [14,14]
# Matrix 2: Calf detection [14,25]
# Matrix 3: Calf-age determination [25,72]

# THE MODEL:
M <- function() {
     
  # PRIORS 
  phiN ~ dunif(0,1)    #non-breeding adult/subadult survival
  phiB ~ dunif(0,1)    #breeding adult survival
  phiYc ~ dunif(0,1)   #YOY & 1yo survival
  phiC ~ dunif(0,1)    #older calf survival
  psiN ~ dunif(0,1)   #breeding transition NB->Byoy
  psiB ~ dunif(0,1)   #breeding transition B->Byoy (or other Bc's)
  pN ~ dunif(0,1)      #non-breeder class adult detection
  pBn ~ dunif(0,1)      #detection of female adult with no calf
  pBc ~ dunif(0,1)      #detection of female breeder with a calf of any age
  deltaYc ~ dunif(0,1) #dependent calf detection
  deltaC ~ dunif(0,1)  #older calf detection

  # calf-age-assignment matrix priors
  gamma ~ dunif(0,1)   # P(calf can be put in an age category vs unknown)
  alphaTy ~ dunif(0,1) # P(YOY calf is identified as such without uncertainty)
  alphaTc ~ dunif(0,1) # P(1,2,3, or 4yo calf identified as such without uncertainty)
  kappaY ~ dunif(0,1)  # P(YOY assigned to J1- category)
  kappaC ~ dunif(0,1)  # P(3 or 4yo assigned to J2+ or J3+ categories respectively)
  omegaA ~ dunif(0,1)  # P(1 or 2yo assigned to uncertain category matching true age: J1+/J1- or J2+/J2-)
  omegaB ~ dunif(0,1)  # P(3 or 4yo assigned to J3+ or J4+ category respectively)
  eta ~ dunif(0,1)     # P(1 or 2yo assigned to J1- or J2- category respectively)

  # Dirichlet priors for initial states 1:13 (Byoy to Bc4c2; excludes dead state which = 0)
  for (i in 1:13) {
       beta[i] ~ dgamma(1,1) #induce Dirichlet prior
       pi[i]<-beta[i]/sum(beta[])
  }
  
  # DEFINE PARAMETERS	
  # probabilities for each INITIAL STATE
  px0[1] <- pi[1]   # prob. of initial state NB
  px0[2] <- pi[2]   # prob. of being in initial state B
  px0[3] <- pi[3]   # prob. of being in initial state Byoy  
  px0[4] <- pi[4]   # prob. of being in initial state Bc1   
  px0[5] <- pi[5]   # prob. of being in initial state Bc2   
  px0[6] <- pi[6]   # prob. of being in initial state Bc2yoy
  px0[7] <- pi[7]   # prob. of being in initial state Bc3   
  px0[8] <- pi[8]   # prob. of being in initial state Bc3yoy
  px0[9] <- pi[9]   # prob. of being in initial state Bc3c1 
  px0[10] <- pi[10] # prob. of being in initial state Bc4   
  px0[11] <- pi[11] # prob. of being in initial state Bc4yoy
  px0[12] <- pi[12] # prob. of being in initial state Bc4c1 
  px0[13] <- pi[13] # prob. of being in initial state Bc4c2
  px0[14] <- 0      # prob. of being in initial state dead  
  
 # OBSERVATION PROCESS: probabilities of observations (columns) at a given occasion given 
 #                      states (rows) at this occasion
  
  # Matrix 1: adult detection [14,14]
  po1[1,1:14]<-c(1-pN,pN,0,0,0,0,0,0,0,0,0,0,0,0)
  po1[2,1:14]<-c(1-pBn,0,pBn,0,0,0,0,0,0,0,0,0,0,0)
  po1[3,1:14]<-c(1-pBc,0,0,pBc,0,0,0,0,0,0,0,0,0,0)
  po1[4,1:14]<-c(1-pBc,0,0,0,pBc,0,0,0,0,0,0,0,0,0)
  po1[5,1:14]<-c(1-pBc,0,0,0,0,pBc,0,0,0,0,0,0,0,0)
  po1[6,1:14]<-c(1-pBc,0,0,0,0,0,pBc,0,0,0,0,0,0,0)
  po1[7,1:14]<-c(1-pBc,0,0,0,0,0,0,pBc,0,0,0,0,0,0)
  po1[8,1:14]<-c(1-pBc,0,0,0,0,0,0,0,pBc,0,0,0,0,0)
  po1[9,1:14]<-c(1-pBc,0,0,0,0,0,0,0,0,pBc,0,0,0,0)
  po1[10,1:14]<-c(1-pBc,0,0,0,0,0,0,0,0,0,pBc,0,0,0)
  po1[11,1:14]<-c(1-pBc,0,0,0,0,0,0,0,0,0,0,pBc,0,0)
  po1[12,1:14]<-c(1-pBc,0,0,0,0,0,0,0,0,0,0,0,pBc,0)
  po1[13,1:14]<-c(1-pBc,0,0,0,0,0,0,0,0,0,0,0,0,pBc)
  po1[14,1:14]<-c(1,0,0,0,0,0,0,0,0,0,0,0,0,0)

  # Matrix 2: Calf Detection [14,25]
  po2[1,1:25]<-c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  po2[2,1:25]<-c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  po2[3,1:25]<-c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  po2[4,1:25]<-c(0,1-deltaYc,deltaYc,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  po2[5,1:25]<-c(0,1-deltaYc,0,deltaYc,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  po2[6,1:25]<-c(0,1-deltaC,0,0,deltaC,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  po2[7,1:25]<-c(0,(1-deltaYc)*(1-deltaC),0,0,0,deltaYc*(1-deltaC),deltaC*(1-deltaYc),
                 0,0,0,0,0,0,0,0,0,0,0,0,deltaYc*deltaC,0,0,0,0,0)
  po2[8,1:25]<-c(0,1-deltaC,0,0,0,0,0,deltaC,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  po2[9,1:25]<-c(0,(1-deltaYc)*(1-deltaC),0,0,0,0,0,0,deltaYc*(1-deltaC),deltaC*(1-deltaYc),
                 0,0,0,0,0,0,0,0,0,0,deltaYc*deltaC,0,0,0,0)
  po2[10,1:25]<-c(0,(1-deltaYc)*(1-deltaC),0,0,0,0,0,0,0,0,deltaYc*(1-deltaC),deltaC*(1-deltaYc),
                  0,0,0,0,0,0,0,0,0,deltaYc*deltaC,0,0,0)
  po2[11,1:25]<-c(0,1-deltaC,0,0,0,0,0,0,0,0,0,0,deltaC,0,0,0,0,0,0,0,0,0,0,0,0)
  po2[12,1:25]<-c(0,(1-deltaYc)*(1-deltaC),0,0,0,0,0,0,0,0,0,0,0,deltaYc*(1-deltaC),
                  deltaC*(1-deltaYc),0,0,0,0,0,0,0,deltaYc*deltaC,0,0)
  po2[13,1:25]<-c(0,(1-deltaYc)*(1-deltaC),0,0,0,0,0,0,0,0,0,0,0,0,0,deltaYc*(1-deltaC),
                  deltaC*(1-deltaYc),0,0,0,0,0,0,deltaYc*deltaC,0)
  po2[14,1:25]<-c(0,(1-deltaC)*(1-deltaC),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,deltaC-deltaC^2,
                  deltaC-deltaC^2,0,0,0,0,0,deltaC^2)

  # Matrix 3: Calf-age-assignment [25,72]
   po3[1,1:72]<-c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
   po3[2,1:72]<-c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
   po3[3,1:72]<-c(0,0,gamma*alphaTy,gamma*(1-alphaTy)*kappaY,0,0,gamma*(1-alphaTy)*(1-kappaY),0,
                  0,0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
   po3[4,1:72]<-c(0,0,0,gamma*(1-alphaTc)*omegaA*eta,gamma*alphaTc,gamma*(1-alphaTc)*omegaA*(1-eta),
                  gamma*(1-alphaTc)*(1-omegaA),0,0,0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
   po3[5,1:72]<-c(0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaA),gamma*(1-alphaTc)*omegaA*eta,gamma*alphaTc,
                  gamma*(1-alphaTc)*omegaA*(1-eta),0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
   po3[6,1:72]<-c(0,0,gamma*alphaTy,gamma*(1-alphaTy)*kappaY,0,0,gamma*(1-alphaTy)*(1-kappaY),0,
                  0,0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
   po3[7,1:72]<-c(0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaA),gamma*(1-alphaTc)*omegaA*eta,
                  gamma*alphaTc,gamma*(1-alphaTc)*omegaA*(1-eta),0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0)
   po3[8,1:72]<-c(0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC),0,0,
                  gamma*(1-alphaTc)*(1-omegaB)*kappaC,gamma*alphaTc,gamma*(1-alphaTc)*omegaB,0,0,
                  1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
   po3[9,1:72]<-c(0,0,gamma*alphaTy,gamma*(1-alphaTy)*kappaY,0,0,gamma*(1-alphaTy)*(1-kappaY),0,
                  0,0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
   po3[10,1:72]<-c(0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC),0,0,
                   gamma*(1-alphaTc)*(1-omegaB)*kappaC,gamma*alphaTc,gamma*(1-alphaTc)*omegaB,0,0,
                   1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
   po3[11,1:72]<-c(0,0,0,gamma*(1-alphaTc)*omegaA*eta,gamma*alphaTc,
                   gamma*(1-alphaTc)*omegaA*(1-eta),gamma*(1-alphaTc)*(1-omegaA),0,0,0,0,0,0,
                   1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
   po3[12,1:72]<-c(0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC),0,0,
                   gamma*(1-alphaTc)*(1-omegaB)*kappaC,gamma*alphaTc,gamma*(1-alphaTc)*omegaB,
                   0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
   po3[13,1:72]<-c(0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,0,
                   gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,gamma*(1-alphaTc)*(1-omegaB)*kappaC,
                   gamma*alphaTc,gamma*(1-alphaTc)*omegaB,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0)
   po3[14,1:72]<-c(0,0,gamma*alphaTy,gamma*(1-alphaTy)*kappaY,0,0,gamma*(1-alphaTy)*(1-kappaY),
                   0,0,0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
   po3[15,1:72]<-c(0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,0,
                   gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,
                   gamma*(1-alphaTc)*(1-omegaB)*kappaC,gamma*alphaTc,gamma*(1-alphaTc)*omegaB,
                   1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
   po3[16,1:72]<-c(0,0,0,gamma*(1-alphaTc)*omegaA*eta,gamma*alphaTc,
                   gamma*(1-alphaTc)*omegaA*(1-eta),gamma*(1-alphaTc)*(1-omegaA),0,0,0,0,0,0,
                   1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
   po3[17,1:72]<-c(0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,0,
                   gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,
                   gamma*(1-alphaTc)*(1-omegaB)*kappaC,gamma*alphaTc,gamma*(1-alphaTc)*omegaB,
                   1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
   po3[18,1:72]<-c(0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaA),gamma*(1-alphaTc)*omegaA*eta,
                   gamma*alphaTc,gamma*(1-alphaTc)*omegaA*(1-eta),0,0,0,0,1-gamma,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0)
   po3[19,1:72]<-c(0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,0,
                   gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,
                   gamma*(1-alphaTc)*(1-omegaB)*kappaC,gamma*alphaTc,gamma*(1-alphaTc)*omegaB,
                   1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
   po3[20,1:72]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*alphaTy)*(gamma*(1-alphaTc)*(1-omegaA)),
                   (gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*(1-omegaA)),0,
                   (gamma*alphaTy)*(gamma*(1-alphaTc)*omegaA*eta),
                   (gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*omegaA*eta),
                   (gamma*alphaTy)*(gamma*alphaTc),(gamma*(1-alphaTy)*kappaY)*(gamma*alphaTc),
                   (gamma*(1-alphaTy)*(1-kappaY))*(gamma*alphaTc),
                   (gamma*alphaTy)*(gamma*(1-alphaTc)*omegaA*(1-eta)),
                   (gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*omegaA*(1-eta)),0,0,
                   (gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*omegaA*(1-eta)),0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*alphaTy)*(1-gamma),
                   (gamma*(1-alphaTy)*kappaY)*(1-gamma),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaA)),
                   ((1-gamma)*(gamma*(1-alphaTc)*omegaA*eta))+((gamma*(1-alphaTy)*(1-kappaY))*(1-gamma)),
                   (1-gamma)*(gamma*alphaTc),(1-gamma)*(gamma*(1-alphaTc)*omegaA*(1-eta)),0,
                   0,0,0,(1-gamma)^2,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*(1-omegaA)),
                   0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*omegaA*eta),0)
   po3[21,1:72]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*alphaTy)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),
                   (gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),0,0,0,
                   0,0,0,(gamma*alphaTy)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),
                   (gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),0,0,
                   (gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),0,
                   (gamma*alphaTy)*(gamma*alphaTc),(gamma*(1-alphaTy)*kappaY)*(gamma*alphaTc),0,0,
                   (gamma*(1-alphaTy)*(1-kappaY))*(gamma*alphaTc),
                   (gamma*alphaTy)*(gamma*(1-alphaTc)*omegaB),
                   (gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*omegaB),0,0,
                   (gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*omegaB),0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,(gamma*alphaTy)*(1-gamma),(gamma*(1-alphaTy)*kappaY)*(1-gamma),0,
                   (1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),
                   (gamma*(1-alphaTy)*(1-kappaY))*(1-gamma),0,
                   (1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(1-gamma)*(gamma*alphaTc),
                   (1-gamma)*(gamma*(1-alphaTc)*omegaB),0,0,(1-gamma)^2,0,
                   (gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),0,0,0,0)
   po3[22,1:72]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   (gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),
                   (gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),0,0,0,0,0,0,
                   (gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),
                   (gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),
                   (gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),
                   (gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),0,0,
                   (gamma*(1-alphaTc)*omegaA*eta)*(gamma*alphaTc),(gamma*alphaTc)*(gamma*alphaTc),
                   (gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*alphaTc),
                   (gamma*(1-alphaTc)*(1-omegaA))*(gamma*alphaTc),0,
                   (gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*omegaB),
                   (gamma*alphaTc)*(gamma*(1-alphaTc)*omegaB),
                   (gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*omegaB),
                   (gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*omegaB),0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(1-gamma),
                   (gamma*alphaTc)*(1-gamma),
                   ((1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)))+
                        ((gamma*(1-alphaTc)*omegaA*(1-eta))*(1-gamma)),
                   (gamma*(1-alphaTc)*(1-omegaA))*(1-gamma),0,
                   (1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(1-gamma)*(gamma*alphaTc),
                   (1-gamma)*(gamma*(1-alphaTc)*omegaB),0,0,(1-gamma)^2,
                   (gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),
                   (gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),0,0,0,0)
   po3[23,1:72]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   (gamma*alphaTy)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),
                   (gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),
                   0,0,0,0,0,0,(gamma*alphaTy)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),
                   (gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,
                   0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),
                   0,0,0,0,0,0,(gamma*alphaTy)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),
                   (gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),0,0,
                   (gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),0,0,
                   (gamma*alphaTy)*(gamma*alphaTc),(gamma*(1-alphaTy)*kappaY)*(gamma*alphaTc),
                   0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*alphaTc),0,0,
                   (gamma*alphaTy)*(gamma*(1-alphaTc)*omegaB),
                   (gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*omegaB),0,0,
                   (gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*omegaB),0,0,
                   (gamma*alphaTy)*(1-gamma),(gamma*(1-alphaTy)*kappaY)*(1-gamma),0,
                   (1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),
                   (gamma*(1-alphaTy)*(1-kappaY))*(1-gamma),0,
                   (1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,
                   (1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(1-gamma)*(gamma*alphaTc),
                   (1-gamma)*(gamma*(1-alphaTc)*omegaB),(1-gamma)^2,0,
                   (gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),
                   0,0,0,0)
   po3[24,1:72]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   (gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),
                   (gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,0,0,0,0,0,
                   (gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),
                   (gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),
                   (gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),
                   (gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),
                   0,0,0,0,0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),
                   (gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),
                   (gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),
                   (gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),0,0,0,
                   (gamma*(1-alphaTc)*omegaA*eta)*(gamma*alphaTc),(gamma*alphaTc)*(gamma*alphaTc),
                   (gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*alphaTc),
                   (gamma*(1-alphaTc)*(1-omegaA))*(gamma*alphaTc),0,0,0,
                   (gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*omegaB),
                   (gamma*alphaTc)*(gamma*(1-alphaTc)*omegaB),
                   (gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*omegaB),
                   (gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*omegaB),0,0,0,
                   (gamma*(1-alphaTc)*omegaA*eta)*(1-gamma),(gamma*alphaTc)*(1-gamma),
                   ((1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5))+
                        ((gamma*(1-alphaTc)*omegaA*(1-eta))*(1-gamma)),
                   (gamma*(1-alphaTc)*(1-omegaA))*(1-gamma),0,
                   (1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,
                   (1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(1-gamma)*(gamma*alphaTc),
                   (1-gamma)*(gamma*(1-alphaTc)*omegaB),
                   (1-gamma)^2,(gamma*(1-alphaTc)*omegaA*(1-eta))*
                        (gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),
                   (gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),
                   0,0,0,0)
   po3[25,1:72]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   (gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),
                   (gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),
                   (gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,0,0,0,0,0,0,0,
                   (gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),
                   (gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),
                   (gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),
                   (gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),
                   0,0,0,(gamma*(1-alphaTc)*(1-omegaA))*(gamma*alphaTc),
                   (gamma*(1-alphaTc)*omegaA*eta)*(gamma*alphaTc),(gamma*alphaTc)*(gamma*alphaTc),
                   (gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*alphaTc),0,0,0,
                   (gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*omegaB),
                   (gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*omegaB),
                   (gamma*alphaTc)*(gamma*(1-alphaTc)*omegaB),
                   (gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*omegaB),0,0,0,
                   ((1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5))+
                        ((gamma*(1-alphaTc)*(1-omegaA))*(1-gamma)),
                   (gamma*(1-alphaTc)*omegaA*eta)*(1-gamma),(gamma*alphaTc)*(1-gamma),
                   ((1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5))+
                        ((gamma*(1-alphaTc)*omegaA*(1-eta))*(1-gamma)),0,
                   (1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(1-gamma)*(gamma*alphaTc),
                   (1-gamma)*(gamma*(1-alphaTc)*omegaB),(1-gamma)^2,
                   (gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),
                   (gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),
                   (gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),
                   (gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),
                   0,(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5))

  
  # Initial State Matrix: [14,14]
  # Populate Initial State Matrix: deterministic initial state
  po1.init[1,1:14]<-c(0,1,0,0,0,0,0,0,0,0,0,0,0,0)
  po1.init[2,1:14]<-c(0,0,1,0,0,0,0,0,0,0,0,0,0,0)
  po1.init[3,1:14]<-c(0,0,0,1,0,0,0,0,0,0,0,0,0,0)
  po1.init[4,1:14]<-c(0,0,0,0,1,0,0,0,0,0,0,0,0,0)
  po1.init[5,1:14]<-c(0,0,0,0,0,1,0,0,0,0,0,0,0,0)
  po1.init[6,1:14]<-c(0,0,0,0,0,0,1,0,0,0,0,0,0,0)
  po1.init[7,1:14]<-c(0,0,0,0,0,0,0,1,0,0,0,0,0,0)
  po1.init[8,1:14]<-c(0,0,0,0,0,0,0,0,1,0,0,0,0,0)
  po1.init[9,1:14]<-c(0,0,0,0,0,0,0,0,0,1,0,0,0,0)
  po1.init[10,1:14]<-c(0,0,0,0,0,0,0,0,0,0,1,0,0,0)
  po1.init[11,1:14]<-c(0,0,0,0,0,0,0,0,0,0,0,1,0,0)
  po1.init[12,1:14]<-c(0,0,0,0,0,0,0,0,0,0,0,0,1,0)
  po1.init[13,1:14]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,1)
  po1.init[14,1:14]<-c(1,0,0,0,0,0,0,0,0,0,0,0,0,0)
  
  
  # Form the matrix product (full observation matrix & initial observation matrix)
  po <- po1 %*% po2 %*% po3 
  po.init <- po1.init %*% po2 %*% po3
  
  

  ## STATE PROCESS: probabilities of states at t+1 (columns) given states at t (rows)
 
  # State Matrix 1: adult survival [14,14]
  px1[1,1:14]<-c(phiN,0,0,0,0,0,0,0,0,0,0,0,0,1-phiN)
  px1[2,1:14]<-c(0,phiB,0,0,0,0,0,0,0,0,0,0,0,1-phiB)
  px1[3,1:14]<-c(0,0,phiB,0,0,0,0,0,0,0,0,0,0,1-phiB)
  px1[4,1:14]<-c(0,0,0,phiB,0,0,0,0,0,0,0,0,0,1-phiB)
  px1[5,1:14]<-c(0,0,0,0,phiB,0,0,0,0,0,0,0,0,1-phiB)
  px1[6,1:14]<-c(0,0,0,0,0,phiB,0,0,0,0,0,0,0,1-phiB)
  px1[7,1:14]<-c(0,0,0,0,0,0,phiB,0,0,0,0,0,0,1-phiB)
  px1[8,1:14]<-c(0,0,0,0,0,0,0,phiB,0,0,0,0,0,1-phiB)
  px1[9,1:14]<-c(0,0,0,0,0,0,0,0,phiB,0,0,0,0,1-phiB)
  px1[10,1:14]<-c(0,0,0,0,0,0,0,0,0,phiB,0,0,0,1-phiB)
  px1[11,1:14]<-c(0,0,0,0,0,0,0,0,0,0,phiB,0,0,1-phiB)
  px1[12,1:14]<-c(0,0,0,0,0,0,0,0,0,0,0,phiB,0,1-phiB)
  px1[13,1:14]<-c(0,0,0,0,0,0,0,0,0,0,0,0,phiB,1-phiB)
  px1[14,1:14]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,1)

 
  # State Matrix 2: calf survival [14,37]
  # This matrix describes calf survival, wherein "survival" means the calf survives AND
  #    stays with its mother to the next year.
  # Departure states = original 14 states (rows)
  # Arrival states (columns):
  # 1	NB        = non-breeder
  # 2	B         = previous breeder with no calf
  # 3	Byoy	     = breeder w/YOY survived
  # 4	Byoy-D	= breeder, YOY died
  # 5	Bc1	     = breeder w/1yo calf alive
  # 6	Bc1-D	= breeder w/1yo calf died
  # 7	Bc2	     = breeder w/2yo calf alive
  # 8	Bc2-D	= breeder w/2yo calf died/left
  # 9	Bc2yoy	= breeder w/2yo calf & YOY, both survive
  # 10    Bc2yoy-D	= breeder w/2yo calf, YOY died
  # 11	Byoyc2-D	= breeder w/YOY, 2yo calf died/left
  # 12	Bc2yoy-DD	= breeder w/2yo calf & YOY, both died/left
  # 13	Bc3	     = breeder w/3yo calf alive
  # 14	Bc3-D	= breeder w/3yo calf died/left
  # 15	Bc3yoy	= breeder w/3yo calf & YOY, both survive
  # 16	Bc3yoy-D	= breeder w/3yo calf, YOY died
  # 17	Byoyc3-D	= breeder w/YOY, 3yo calf died/left
  # 18	Bc3yoy-DD	= breeder w/3yo calf & YOY, both died/left
  # 19	Bc3c1	= breeder w/3yo calf & 1yo calf, both alive
  # 20	Bc3c1-D	= breeder w/3yo calf, 1yo calf died
  # 21	Bc1c3-D	= breeder w/1yo calf, 3yo calf died/left
  # 22	Bc3c1-DD	= breeder w/3yo calf & 1yo calf, both died/left
  # 23	Bc4*	     = breeder w/4yo calf alive
  # 24	Bc4*-D	= breeder w/4yo calf that died/left 
  # 25	Bc4*yoy	= breeder w/4yo & YOY
  # 26	Bc4*yoy-D	= breeder w/4yo calf, YOY died
  # 27	Byoyc4*-D	= breeder w/YOY, 4yo died/left
  # 28	Bc4*yoy-DD= breeder w/4yo & YOY, both died/left
  # 29	Bc4*c1	= breeder w/4yo calf & 1yo calf, both alive
  # 30	Bc4*c1-D	= breeder w/4yo calf, 1yo calf died
  # 31	Bc1c4*-D	= breeder w/1yo calf, 4yo died/left
  # 32	Bc4*c1-DD	= breeder w/4yo calf & 1yo calf, both died/left
  # 33	Bc4*c2	= breeder w/4yo calf & 2yo calf, both alive
  # 34	Bc4*c2-D	= breeder w/4yo calf, 2yo calf died/left
  # 35	Bc2c4*-D	= breeder w/2yo calf, 4yo died/left
  # 36	Bc4*c2-DD	= breeder w/2yo dead & 4yo calf died/left
  # 37	D	     = dead

  px2[1,1:37]<-c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  px2[2,1:37]<-c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  px2[3,1:37]<-c(0,0,phiYc,1-phiYc,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  px2[4,1:37]<-c(0,0,0,0,phiYc,1-phiYc,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  px2[5,1:37]<-c(0,0,0,0,0,0,phiC,1-phiC,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  px2[6,1:37]<-c(0,0,0,0,0,0,0,0,phiC*phiYc,phiC*(1-phiYc),(1-phiC)*phiYc,(1-phiC)*(1-phiYc),0,
                 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  px2[7,1:37]<-c(0,0,0,0,0,0,0,0,0,0,0,0,phiC,1-phiC,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  px2[8,1:37]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,phiC*phiYc,phiC*(1-phiYc),(1-phiC)*phiYc,
                 (1-phiC)*(1-phiYc),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  px2[9,1:37]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,phiC*phiYc,phiC*(1-phiYc),(1-phiC)*phiYc,
                 (1-phiYc)*(1-phiC),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  px2[10,1:37]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,phiC,1-phiC,0,0,0,0,0,0,0,0,0,0,0,0,0)
  px2[11,1:37]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,phiC*phiYc,phiC*(1-phiYc),
                  (1-phiC)*phiYc,(1-phiC)*(1-phiYc),0,0,0,0,0,0,0,0,0)
  px2[12,1:37]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,phiC*phiYc,
                  phiC*(1-phiYc),phiYc*(1-phiC),(1-phiC)*(1-phiYc),0,0,0,0,0)
  px2[13,1:37]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  0,phiC,1-phiC,0)
  px2[14,1:37]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)
  
  # State Matrix 3: young aging [37,14]
  # This is a deterministic matrix describing a calf moving to the next age
  #   given that it survived and stayed with its mother.
  
  # Departure states = 37 states from young survival matrix (see above)
  # Arrival states:
  # 1	NB = non-breeder
  # 2	B = previous breeder with no calf
  # 3	Bc1 = breeder w/calf who was a YOY and now a 1yo                              
  # 4	Bc2 = breeder w/a now 2yo calf                                                
  # 5	Bc3 = breeder w/a now 3yo calf                                                
  # 6	Bc3c1 = breeder w/a now 3yo calf & now 1yo calf                               
  # 7	Bc4 = breeder w/a now 4yo (or older) calf                                     
  # 8	Bc4c1 = breeder w/a now 4yo (or older) calf & now 1yo calf                    
  # 9	Bc4c2 = breeder w/a now 4yo (or older) calf & now 2yo calf                    
  # 10    Byoy-D = breeder whose YOY died (no calves remaining)                         
  # 11	Bc-D-L = breeder whose calf (or calves) died (or left) (no calves remaining)
  # 12	Bc3yoy-D = breeder with a now 3yo (and YOY who died)                        
  # 13	Bc4yoy-D = breeder with a now 4yo (or older) and YOY who died               
  # 14 	D = dead adult   
      	                                                           
  px3[1,1:14]<-c(1,0,0,0,0,0,0,0,0,0,0,0,0,0)
  px3[2,1:14]<-c(0,1,0,0,0,0,0,0,0,0,0,0,0,0)
  px3[3,1:14]<-c(0,0,1,0,0,0,0,0,0,0,0,0,0,0)
  px3[4,1:14]<-c(0,0,0,0,0,0,0,0,0,1,0,0,0,0)
  px3[5,1:14]<-c(0,0,0,1,0,0,0,0,0,0,0,0,0,0)
  px3[6,1:14]<-c(0,0,0,0,0,0,0,0,0,0,1,0,0,0)
  px3[7,1:14]<-c(0,0,0,0,1,0,0,0,0,0,0,0,0,0)
  px3[8,1:14]<-c(0,0,0,0,0,0,0,0,0,0,1,0,0,0)
  px3[9,1:14]<-c(0,0,0,0,0,1,0,0,0,0,0,0,0,0)
  px3[10,1:14]<-c(0,0,0,0,0,0,0,0,0,0,0,1,0,0)
  px3[11,1:14]<-c(0,0,1,0,0,0,0,0,0,0,0,0,0,0)
  px3[12,1:14]<-c(0,0,0,0,0,0,0,0,0,1,0,0,0,0)
  px3[13,1:14]<-c(0,0,0,0,0,0,1,0,0,0,0,0,0,0)
  px3[14,1:14]<-c(0,0,0,0,0,0,0,0,0,0,1,0,0,0)
  px3[15,1:14]<-c(0,0,0,0,0,0,0,1,0,0,0,0,0,0)
  px3[16,1:14]<-c(0,0,0,0,0,0,0,0,0,0,0,0,1,0)
  px3[17,1:14]<-c(0,0,1,0,0,0,0,0,0,0,0,0,0,0)
  px3[18,1:14]<-c(0,0,0,0,0,0,0,0,0,1,0,0,0,0)
  px3[19,1:14]<-c(0,0,0,0,0,0,0,0,1,0,0,0,0,0)
  px3[20,1:14]<-c(0,0,0,0,0,0,1,0,0,0,0,0,0,0)
  px3[21,1:14]<-c(0,0,0,1,0,0,0,0,0,0,0,0,0,0)
  px3[22,1:14]<-c(0,0,0,0,0,0,0,0,0,0,1,0,0,0)
  px3[23,1:14]<-c(0,0,0,0,0,0,1,0,0,0,0,0,0,0)
  px3[24,1:14]<-c(0,0,0,0,0,0,0,0,0,0,1,0,0,0)
  px3[25,1:14]<-c(0,0,0,0,0,0,0,1,0,0,0,0,0,0)
  px3[26,1:14]<-c(0,0,0,0,0,0,0,0,0,0,0,0,1,0)
  px3[27,1:14]<-c(0,0,1,0,0,0,0,0,0,0,0,0,0,0)
  px3[28,1:14]<-c(0,0,0,0,0,0,0,0,0,1,0,0,0,0)
  px3[29,1:14]<-c(0,0,0,0,0,0,0,0,1,0,0,0,0,0)
  px3[30,1:14]<-c(0,0,0,0,0,0,1,0,0,0,0,0,0,0)
  px3[31,1:14]<-c(0,0,0,1,0,0,0,0,0,0,0,0,0,0)
  px3[32,1:14]<-c(0,0,0,0,0,0,0,0,0,0,1,0,0,0)
  px3[33,1:14]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  px3[34,1:14]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  px3[35,1:14]<-c(0,0,0,0,1,0,0,0,0,0,0,0,0,0)
  px3[36,1:14]<-c(0,0,0,0,0,0,0,0,0,0,1,0,0,0)
  px3[37,1:14]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,1)

  
  # State Matrix 4: breeding transition [14,14]
  # This matrix defines the probability of an adult giving birth to a YOY.
  #   Because we assume that births cannot happen in consecutive years, an adult
  #   can only transition from the NB state or from one of the "breeder states"
  #   in which they have a calf 2 years old or older. In the latter case, the 
  #   adult will be transitioning to a 2-calf state.

  # Departure states = 14 states from young aging matrix (see above)
  # Arrival states = original 14 states (NB, Byoy, Bc1,...,D)
  px4[1,1:14]<-c(1-psiN,0,psiN,0,0,0,0,0,0,0,0,0,0,0)
  px4[2,1:14]<-c(0,1-psiB,psiB,0,0,0,0,0,0,0,0,0,0,0)
  px4[3,1:14]<-c(0,0,0,1,0,0,0,0,0,0,0,0,0,0)
  px4[4,1:14]<-c(0,0,0,0,1-psiB,psiB,0,0,0,0,0,0,0,0)
  px4[5,1:14]<-c(0,0,0,0,0,0,1-psiB,psiB,0,0,0,0,0,0)
  px4[6,1:14]<-c(0,0,0,0,0,0,0,0,1,0,0,0,0,0)
  px4[7,1:14]<-c(0,0,0,0,0,0,0,0,0,1-psiB,psiB,0,0,0)
  px4[8,1:14]<-c(0,0,0,0,0,0,0,0,0,0,0,1,0,0)
  px4[9,1:14]<-c(0,0,0,0,0,0,0,0,0,0,0,0,1,0)
  px4[10,1:14]<-c(0,1,0,0,0,0,0,0,0,0,0,0,0,0)
  px4[11,1:14]<-c(0,1-psiB,psiB,0,0,0,0,0,0,0,0,0,0,0)
  px4[12,1:14]<-c(0,0,0,0,0,0,1,0,0,0,0,0,0,0)
  px4[13,1:14]<-c(0,0,0,0,0,0,0,0,0,1,0,0,0,0)
  px4[14,1:14]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,1)
  
  # form the matrix product (full state matrix)
  px <- px1 %*% px2 %*% px3 %*% px4
  
  #Likelihoods
  for (i in 1:N)  # loop over each individual
  {
    # estimated probabilities of initial states are the proportions in each state at 
    #    first capture occasion
    alive[i,First[i]] ~ dcat(px0[1:14])
    mydata[i,First[i]] ~ dcat(po.init[alive[i,First[i]],1:72])
    
    for (j in (First[i]+1):Years)  # loop over time
    {
      ## STATE EQUATIONS ##
      # draw states at j given states at j-1
      alive[i,j] ~ dcat(px[alive[i,j-1],1:14])
      
      ## OBSERVATION EQUATIONS ##
      # draw observations at j given states at j
      mydata[i,j] ~ dcat(po[alive[i,j],1:72])
    }
  }

  #Calculate derived calf-survival parameter
  ##  (accounts for calves dying when mother dies)
  phiYc.d <- phiB * phiYc
  phiC.d <- phiB * phiC
  
}

#write model code to file
mod.file.name<-paste("jags_CIBW_ME_Model_",data.file,".txt",sep="")
write.model(M, mod.file.name)
model.file <- paste(getwd(),mod.file.name, sep="/")


## Bundle the data
mydatax = list(N=N,Years=K,mydata=as.matrix(mydata),First=e)


## Initial values (to aid in convergence)
initfunction <- function() list(pN=runif(1, 0.45, 0.5),phiB=runif(1,0.98,0.99),
                                alphaTy=runif(1,0.4,0.6), pBn=runif(1,0.75,0.85),
                                phiYc=runif(1,0.92,0.97), deltaC=runif(1,0.4,0.5),
                                deltaYc=runif(1,0.48,0.6), kappaC=runif(1,0.5,0.7),
                                omegaB=runif(1,0.5,0.7), phiC=runif(1,0.45,0.6),
                                alive=as.matrix(alive.function()))


# Specify the parameters to be monitored
parameters <- c("phiN","phiB","phiC","phiC.d","phiYc","phiYc.d","psiN","psiB","pN","pBn","pBc","deltaYc",
                "deltaC","alphaTy","alphaTc","gamma","kappaY","kappaC","omegaA",
                "omegaB","eta","pi")

#MCMC settings
nc<-4       #number of chains
ni<-55000   #number of iterations per chain
nbi<-10000  #number of burn-in iterations
nt<-10      #thinning rate


######
## Fit model with JAGS:

# for testing (short run-time)
# out <- jags(mydatax,inits=initfunction,parameters,model.file,n.chains=3,n.iter=500,
#             n.burnin=100,n.cores=3)

# Full model run (takes many hours to run)
out <- jags(mydatax,inits=initfunction,parameters,model.file,n.chains=nc,n.iter=ni,
            n.thin=nt,n.burnin=nbi,parallel=T,n.cores=nc)



##########################################
## MODEL OUTPUT 
##########################################

# Check convergence:
plot(out, ask=F)
densityplot(out,ask=F)

samples<-coda::as.mcmc.list(out$samples)
gelman.plot(samples, ask=F)

# Print results
print(out)

