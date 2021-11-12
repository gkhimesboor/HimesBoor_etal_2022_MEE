
library(dplyr)
library(tidyverse)
library(R2OpenBUGS)
library(knitr)
library(jagsUI)
library(coda)
library(here)

source(here::here('PlotTheme.R'))

figs <- captioner(prefix = "Figure")
tbls <- captioner(prefix = "Table")


#function to generate data
# Define function to simulate ME capture-recapture data
ME.sim <- function(STATE, OBS, STATE.INIT, OBS.INIT, marked){
     n.occasions <- dim(STATE)[4] + 1
     CH <- CH.TRUE <- matrix(NA, ncol = n.occasions, nrow = sum(marked))
     # Define a vector with the occasion of marking
     mark.occ <- rep(1:length(marked), marked)
     
     for (i in 1:sum(marked)){
          # Initial state
          CH.TRUE[i,mark.occ[i]] <- which(rmultinom(1, 1, STATE.INIT[,i,mark.occ[i]])==1)
          # Event at first detection 
          CH[i,mark.occ[i]] <- which(rmultinom(1, 1, OBS.INIT[CH.TRUE[i,mark.occ[i]],,i,mark.occ[i]])==1)
          for (t in (mark.occ[i]+1):n.occasions){
               # Multinomial trials for state transitions
               if (mark.occ[i]==n.occasions) next
               state <- which(rmultinom(1, 1, STATE[CH.TRUE[i,t-1],,i,t-1])==1)
               CH.TRUE[i,t] <- state
               # Multinomial trials for observation process
               event <- which(rmultinom(1, 1, OBS[CH.TRUE[i,t],,i,t-1])==1)
               CH[i,t] <- event
          } #t
     } #i
     
     return(list(CH=CH, CH.TRUE=CH.TRUE))  # CH: to be used as data; CH.TRUE: CH with perfect observations
     
} # me.sim fun

# # function to set initial values for z state.... help
# alive.function <- function(){
#      # z <- start.mat.csv #help - how do we generalize or omit this function? 
#      z <- matrix(1, ncol = n.occasions, nrow = totrel) #help - how do we generalize or omit this function?
#      for (i in 1:n.occasions){
#           for (j in 1:totrel) {
#                if (j < fc[i]) {z[i,j] = NA} #
#           }
#      }
#      z
# }


# THE MODEL:
M <- function() {
     
     # PRIORS 
     phiN ~ dunif(0,1)    #non-breeding adult/subadult survival
     phiB ~ dunif(0,1)    #breeding adult survival
     phiYc ~ dunif(0,1)   #YOY & 1yo survival
     phiC ~ dunif(0,1)    #older calf survival
     psiN ~ dunif(0,psiB) #breeding transition NB->Byoy
     psiB ~ dunif(0,0.5)  #breeding transition B->Byoy (or other Bc's)
     pN ~ dunif(0,1)      #non-breeder class adult detection
     pBn ~ dunif(0,1)      #detection of female adult with no calf
     pBc ~ dunif(0,1)      #detection of female breeder with a calf of any age
     deltaYc ~ dunif(0,1) #dependent calf detection
     deltaC ~ dunif(0,1)  #older calf detection
     
     # calf aging parameters 
     gamma ~ dunif(0,1)   # P(calf can be put in an age category vs unknown)
     alphaTy ~ dunif(0,1) # P(a YOY calf is identified as such without uncertainty)
     alphaTc ~ dunif(0,1) # P(a 1,2,3, or 4yo calf is identified as such without uncertainty)
     kappaY ~ dunif(0,1)  # P(a YOY is assigned to J1- category)
     kappaC ~ dunif(0,1)  # P(a 3 or 4yo is assigned to the J2+ or J3+ categories respectively)
     omegaA ~ dunif(0,1)  # P()
     omegaB ~ dunif(0,1)  # P(a 3 or 4yo is assigned to the J3+ or J4+ categories respectively)
     eta ~ dunif(0,1)     # P()
     
     #set Dirichlet priors for initial states 1:12 (Byoy to Bc4c2; excludes dead state which = 0)
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
     
     # Matrix 3: calf-age assignment [25,72]
     po3[1,1:72]<-c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
     po3[2,1:72]<-c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
     po3[3,1:72]<-c(0,0,gamma*alphaTy,gamma*(1-alphaTy)*kappaY,0,0,gamma*(1-alphaTy)*(1-kappaY),0,0,0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
     po3[4,1:72]<-c(0,0,0,gamma*(1-alphaTc)*omegaA*eta,gamma*alphaTc,gamma*(1-alphaTc)*omegaA*(1-eta),gamma*(1-alphaTc)*(1-omegaA),0,0,0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
     po3[5,1:72]<-c(0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaA),gamma*(1-alphaTc)*omegaA*eta,gamma*alphaTc,gamma*(1-alphaTc)*omegaA*(1-eta),0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
     po3[6,1:72]<-c(0,0,gamma*alphaTy,gamma*(1-alphaTy)*kappaY,0,0,gamma*(1-alphaTy)*(1-kappaY),0,0,0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
     po3[7,1:72]<-c(0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaA),gamma*(1-alphaTc)*omegaA*eta,gamma*alphaTc,gamma*(1-alphaTc)*omegaA*(1-eta),0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
     po3[8,1:72]<-c(0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC),0,0,gamma*(1-alphaTc)*(1-omegaB)*kappaC,gamma*alphaTc,gamma*(1-alphaTc)*omegaB,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
     po3[9,1:72]<-c(0,0,gamma*alphaTy,gamma*(1-alphaTy)*kappaY,0,0,gamma*(1-alphaTy)*(1-kappaY),0,0,0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
     po3[10,1:72]<-c(0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC),0,0,gamma*(1-alphaTc)*(1-omegaB)*kappaC,gamma*alphaTc,gamma*(1-alphaTc)*omegaB,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
     po3[11,1:72]<-c(0,0,0,gamma*(1-alphaTc)*omegaA*eta,gamma*alphaTc,gamma*(1-alphaTc)*omegaA*(1-eta),gamma*(1-alphaTc)*(1-omegaA),0,0,0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
     po3[12,1:72]<-c(0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC),0,0,gamma*(1-alphaTc)*(1-omegaB)*kappaC,gamma*alphaTc,gamma*(1-alphaTc)*omegaB,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
     po3[13,1:72]<-c(0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,gamma*(1-alphaTc)*(1-omegaB)*kappaC,gamma*alphaTc,gamma*(1-alphaTc)*omegaB,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
     po3[14,1:72]<-c(0,0,gamma*alphaTy,gamma*(1-alphaTy)*kappaY,0,0,gamma*(1-alphaTy)*(1-kappaY),0,0,0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
     po3[15,1:72]<-c(0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,gamma*(1-alphaTc)*(1-omegaB)*kappaC,gamma*alphaTc,gamma*(1-alphaTc)*omegaB,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
     po3[16,1:72]<-c(0,0,0,gamma*(1-alphaTc)*omegaA*eta,gamma*alphaTc,gamma*(1-alphaTc)*omegaA*(1-eta),gamma*(1-alphaTc)*(1-omegaA),0,0,0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
     po3[17,1:72]<-c(0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,gamma*(1-alphaTc)*(1-omegaB)*kappaC,gamma*alphaTc,gamma*(1-alphaTc)*omegaB,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
     po3[18,1:72]<-c(0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaA),gamma*(1-alphaTc)*omegaA*eta,gamma*alphaTc,gamma*(1-alphaTc)*omegaA*(1-eta),0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
     po3[19,1:72]<-c(0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,gamma*(1-alphaTc)*(1-omegaB)*kappaC,gamma*alphaTc,gamma*(1-alphaTc)*omegaB,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
     po3[20,1:72]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*alphaTy)*(gamma*(1-alphaTc)*(1-omegaA)),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*(1-omegaA)),0,(gamma*alphaTy)*(gamma*(1-alphaTc)*omegaA*eta),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*omegaA*eta),(gamma*alphaTy)*(gamma*alphaTc),(gamma*(1-alphaTy)*kappaY)*(gamma*alphaTc),(gamma*(1-alphaTy)*(1-kappaY))*(gamma*alphaTc),(gamma*alphaTy)*(gamma*(1-alphaTc)*omegaA*(1-eta)),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*omegaA*(1-eta)),0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*omegaA*(1-eta)),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*alphaTy)*(1-gamma),(gamma*(1-alphaTy)*kappaY)*(1-gamma),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaA)),((1-gamma)*(gamma*(1-alphaTc)*omegaA*eta))+((gamma*(1-alphaTy)*(1-kappaY))*(1-gamma)),(1-gamma)*(gamma*alphaTc),(1-gamma)*(gamma*(1-alphaTc)*omegaA*(1-eta)),0,0,0,0,(1-gamma)^2,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*(1-omegaA)),0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*omegaA*eta),0)
     po3[21,1:72]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*alphaTy)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),0,0,0,0,0,0,(gamma*alphaTy)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),0,(gamma*alphaTy)*(gamma*alphaTc),(gamma*(1-alphaTy)*kappaY)*(gamma*alphaTc),0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*alphaTc),(gamma*alphaTy)*(gamma*(1-alphaTc)*omegaB),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*omegaB),0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*omegaB),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*alphaTy)*(1-gamma),(gamma*(1-alphaTy)*kappaY)*(1-gamma),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),(gamma*(1-alphaTy)*(1-kappaY))*(1-gamma),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(1-gamma)*(gamma*alphaTc),(1-gamma)*(gamma*(1-alphaTc)*omegaB),0,0,(1-gamma)^2,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),0,0,0,0)
     po3[22,1:72]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),(gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),0,0,0,0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),0,0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*alphaTc),(gamma*alphaTc)*(gamma*alphaTc),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*alphaTc),(gamma*(1-alphaTc)*(1-omegaA))*(gamma*alphaTc),0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*omegaB),(gamma*alphaTc)*(gamma*(1-alphaTc)*omegaB),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*omegaB),(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*omegaB),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(1-gamma),(gamma*alphaTc)*(1-gamma),((1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)))+((gamma*(1-alphaTc)*omegaA*(1-eta))*(1-gamma)),(gamma*(1-alphaTc)*(1-omegaA))*(1-gamma),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(1-gamma)*(gamma*alphaTc),(1-gamma)*(gamma*(1-alphaTc)*omegaB),0,0,(1-gamma)^2,(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),0,0,0,0)
     po3[23,1:72]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*alphaTy)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,0,0,0,0,0,(gamma*alphaTy)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,0,0,0,0,0,(gamma*alphaTy)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),0,0,(gamma*alphaTy)*(gamma*alphaTc),(gamma*(1-alphaTy)*kappaY)*(gamma*alphaTc),0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*alphaTc),0,0,(gamma*alphaTy)*(gamma*(1-alphaTc)*omegaB),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*omegaB),0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*omegaB),0,0,(gamma*alphaTy)*(1-gamma),(gamma*(1-alphaTy)*kappaY)*(1-gamma),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*(1-alphaTy)*(1-kappaY))*(1-gamma),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(1-gamma)*(gamma*alphaTc),(1-gamma)*(gamma*(1-alphaTc)*omegaB),(1-gamma)^2,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,0,0,0)
     po3[24,1:72]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,0,0,0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,0,0,0,0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*alphaTc),(gamma*alphaTc)*(gamma*alphaTc),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*alphaTc),(gamma*(1-alphaTc)*(1-omegaA))*(gamma*alphaTc),0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*omegaB),(gamma*alphaTc)*(gamma*(1-alphaTc)*omegaB),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*omegaB),(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*omegaB),0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(1-gamma),(gamma*alphaTc)*(1-gamma),((1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5))+((gamma*(1-alphaTc)*omegaA*(1-eta))*(1-gamma)),(gamma*(1-alphaTc)*(1-omegaA))*(1-gamma),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(1-gamma)*(gamma*alphaTc),(1-gamma)*(gamma*(1-alphaTc)*omegaB),(1-gamma)^2,(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,0,0,0)
     po3[25,1:72]<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,0,0,0,0,0,0,0,(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),0,0,0,(gamma*(1-alphaTc)*(1-omegaA))*(gamma*alphaTc),(gamma*(1-alphaTc)*omegaA*eta)*(gamma*alphaTc),(gamma*alphaTc)*(gamma*alphaTc),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*alphaTc),0,0,0,(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*omegaB),(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*omegaB),(gamma*alphaTc)*(gamma*(1-alphaTc)*omegaB),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*omegaB),0,0,0,((1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5))+((gamma*(1-alphaTc)*(1-omegaA))*(1-gamma)),(gamma*(1-alphaTc)*omegaA*eta)*(1-gamma),(gamma*alphaTc)*(1-gamma),((1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5))+((gamma*(1-alphaTc)*omegaA*(1-eta))*(1-gamma)),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(1-gamma)*(gamma*alphaTc),(1-gamma)*(gamma*(1-alphaTc)*omegaB),(1-gamma)^2,(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5))
     
     
     # Initial Matrix: [14,14]
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
     
     
     # form the matrix product
     po <- po1 %*% po2 %*% po3
     po.init <- po1.init %*% po2 %*% po3
     
     # STATE PROCESS: probabilities of states at t+1 (columns) given states at t (rows)
     
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
     # Departure states = original 14 states
     # Arrival states:
     # 1	NB = non-breeder
     # 2	B = previous breeder with no calf
     # 3	Byoy	= breeder w/YOY survived
     # 4	Byoy-D	= breeder, YOY dead
     # 5	Bc1	= breeder w/1yo calf alive
     # 6	Bc1-D	= breeder w/1yo calf dead
     # 7	Bc2	= breeder w/2yo calf alive
     # 8	Bc2-D	= breeder w/2yo calf dead
     # 9	Bc2yoy	= breeder w/2yo calf & YOY, both survive
     # 10    Bc2yoy-D	= breeder w/2yo calf, YOY dead
     # 11	Byoyc2-D	= breeder w/YOY, 2yo calf dead
     # 12	Bc2yoy-DD	= breeder w/2yo calf & YOY, both dead
     # 13	Bc3	= breeder w/3yo calf alive
     # 14	Bc3-D	= breeder w/3yo calf dead
     # 15	Bc3yoy	= breeder w/3yo calf & YOY, both survive
     # 16	Bc3yoy-D	= breeder w/3yo calf, YOY dead
     # 17	Byoyc3-D	= breeder w/YOY, 3yo calf dead
     # 18	Bc3yoy-DD	= breeder w/3yo calf & YOY, both dead
     # 19	Bc3c1	= breeder w/3yo calf & 1yo calf, both alive
     # 20	Bc3c1-D	= breeder w/3yo calf, 1yo calf dead
     # 21	Bc1c3-D	= breeder w/1yo calf, 3yo calf dead
     # 22	Bc3c1-DD	= breeder w/3yo calf & 1yo calf, both dead
     # 23	Bc4*	=breeder w/4yo calf alive
     # 24	Bc4*-D	= breeder w/4yo calf that died/left mother
     # 25	Bc4*yoy	= breeder w/4yo & YOY
     # 26	Bc4*yoy-D	= breeder w/4yo calf, YOY dead
     # 27	Byoyc4*-D	= breeder w/YOY, 4yo died/left
     # 28	Bc4*yoy-DD	= breeder w/4yo & YOY, both dead/left
     # 29	Bc4*c1	= breeder w/4yo calf & 1yo calf, both alive
     # 30	Bc4*c1-D	= breeder w/4yo calf, 1yo calf dead
     # 31	Bc1c4*-D	= breeder w/1yo calf, 4yo died/left
     # 32	Bc4*c1-DD	= breeder w/4yo calf & 1yo calf, both dead/left
     # 33	Bc4*c2	= breeder w/4yo calf & 2yo calf, both alive
     # 34	Bc4*c2-D	= breeder w/4yo calf, 2yo calf dead
     # 35	Bc2c4*-D	= breeder w/2yo calf, 4yo died/left
     # 36	Bc4*c2-DD	= breeder w/2yo dead & 4yo calf dead/left
     # 37	D	= dead
     
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
                     0,phiC,(1-phiC),0)
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
     
     # form the matrix product
     px <- px1 %*% px2 %*% px3 %*% px4
     
     #Likelihoods
     for (i in 1:N)  # loop over each individual
     {
          # estimated probabilities of initial states are the proportions in each state at 
          #    first capture occasion
          z[i,First[i]] ~ dcat(px0[1:14])
          y[i,First[i]] ~ dcat(po.init[z[i,First[i]],1:72])
          
          for (j in (First[i]+1):Years)  # loop over time
          {
               ## STATE EQUATIONS ##
               # draw states at j given states at j-1
               z[i,j] ~ dcat(px[z[i,j-1],1:14])
               
               ## OBSERVATION EQUATIONS ##
               # draw observations at j given states at j
               y[i,j] ~ dcat(po[z[i,j],1:72])
          }
     }
}

#write model code to file
mod.file.name<-paste("jags_models/v4.6-ME-sim",".txt",sep="")
write.model(M, mod.file.name)
model.file <- paste(getwd(),mod.file.name, sep="/")


nsim <- 50

scenario <- '.4.6.1'

pi[1:12] <- c(0.842, 0.003, 0.061, 0.041, 0.008, 0.005, 0.006, 0.005, 0.009, 0.004, 0.004, 0.006)
pi[13] <- 1-sum(pi[1:12])

#### scenarios #### 
#parameter order below
#phiN, phiB, phiC, phiY, psiN, psiB, pN, pBn, pBc, deltaY, deltaC,
#N, gamma, alphaTy, alphaTc, kappaY, kappaC, omegaA, omegaB, eta, pi

#vital rates are held constant based on empirical estimates
#detection rates and sample population size vary

pars.4.6.1 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                0.75, 0.75, 0.75, 0.75, 0.75, #detection: L/M/H
                475, #N: L/M/H
                #gamma, alphaTy, alphaTc, kappaY, kappaC, omegaA, omegaB, eta
                0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)
pars.4.6.2 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                0.5, 0.5, 0.5, 0.5, 0.5, #detection: L/M/H
                475, #N: L/M/H
                0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)
pars.4.6.3 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                0.25, 0.25, 0.25, 0.25, 0.25, #detection: L/M/H
                475, #N: L/M/H
                0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)
pars.4.6.4 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                0.75, 0.75, 0.75, 0.75, 0.75, #detection: L/M/H
                1000, #N: L/M/H
                0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)
pars.4.6.5 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                0.5, 0.5, 0.5, 0.5, 0.5, #detection: L/M/H
                1000, #N: L/M/H
                0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)
pars.4.6.6 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                0.25, 0.25, 0.25, 0.25, 0.25, #detection: L/M/H
                1000, #N: L/M/H
                0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)
pars.4.6.7 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                0.75, 0.75, 0.75, 0.75, 0.75, #detection: L/M/H
                250, #N: L/M/H
                0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)
pars.4.6.8 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                0.5, 0.5, 0.5, 0.5, 0.5, #detection: L/M/H
                250, #N: L/M/H
                0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)
pars.4.6.9 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                0.25, 0.25, 0.25, 0.25, 0.25, #detection: L/M/H
                250, #N: L/M/H
                0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)

pars.4.6.10 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                 0.75, 0.75, 0.75, 0.25, 0.25, #detection: L/M/H
                 475, #N: L/M/H
                 0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)

pars.4.6.11 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                 0.5, 0.5, 0.5, 0.25, 0.25, #detection: L/M/H
                 475, #N: L/M/H
                 0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)

pars.4.6.12 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                 0.25, 0.25, 0.25, 0.5, 0.5, #detection: L/M/H
                 475, #N: L/M/H
                 0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)
pars.4.6.13 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                 0.70, 0.8, 0.7, 0.54, 0.48, #detection: L/M/H
                 475, #N: L/M/H
                 0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)
pars.4.6.14 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                 0.75, 0.75, 0.75, 0.5, 0.5, #detection: L/M/H
                 475, #N: L/M/H
                 0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)
pars.4.6.15 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                 0.5, 0.5, 0.5, 0.75, 0.75, #detection: L/M/H
                 475, #N: L/M/H
                 0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)
pars.4.6.16 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                 0.25, 0.25, 0.25, 0.75, 0.75, #detection: L/M/H
                 475, #N: L/M/H
                 0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)

pars.4.6.17 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                 0.75, 0.75, 0.75, 0.25, 0.25, #detection: L/M/H
                 250, #N: L/M/H
                 0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)
pars.4.6.18 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                 0.75, 0.75, 0.75, 0.5, 0.5, #detection: L/M/H
                 250, #N: L/M/H
                 0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)
pars.4.6.19 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                 0.5, 0.5, 0.5, 0.25, 0.25, #detection: L/M/H
                 250, #N: L/M/H
                 0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)
pars.4.6.20 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                 0.5, 0.5, 0.5, 0.75, 0.75, #detection: L/M/H
                 250, #N: L/M/H
                 0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)
pars.4.6.21 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                 0.25, 0.25, 0.25, 0.5, 0.5, #detection: L/M/H
                 250, #N: L/M/H
                 0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)
pars.4.6.22 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                 0.25, 0.25, 0.25, 0.75, 0.75, #detection: L/M/H
                 250, #N: L/M/H
                 0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)
pars.4.6.23 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                 0.75, 0.75, 0.75, 0.25, 0.25, #detection: L/M/H
                 1000, #N: L/M/H
                 0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)
pars.4.6.24 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                 0.75, 0.75, 0.75, 0.5, 0.5, #detection: L/M/H
                 1000, #N: L/M/H
                 0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)
pars.4.6.25 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                 0.5, 0.5, 0.5, 0.25, 0.25, #detection: L/M/H
                 1000, #N: L/M/H
                 0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)
pars.4.6.26 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                 0.5, 0.5, 0.5, 0.75, 0.75, #detection: L/M/H
                 1000, #N: L/M/H
                 0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)
pars.4.6.27 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                 0.25, 0.25, 0.25, 0.5, 0.5, #detection: L/M/H
                 1000, #N: L/M/H
                 0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)
pars.4.6.28 <- c(0.93, 0.96, 0.51, 0.96, 0.07, 0.28, #phi and psi - empirical
                 0.25, 0.25, 0.25, 0.75, 0.75, #detection: L/M/H
                 1000, #N: L/M/H
                 0.93, 0.5, 0.05, 0.98, 0.71, 0.96, 0.67, 0.03, pi)

#storage
parameters <- c("phiN","phiB", "phiC","phiYc","psiN","psiB", "pN","pBc", "pBn",
                "deltaC", "deltaYc", 
                'gamma', 'alphaTy', 'alphaTc', 'kappaY', 'kappaC', 'omegaA', 'omegaB', 'eta',
                paste0("pi", "[", 1:13, "]"))

n_pars <- length(parameters)

out_median <- matrix(NA, nsim, n_pars-13+1) #not looking at bias for pi params, add space for rhat
out_true <- matrix(NA, nsim, n_pars-13) 
out_cri <- matrix(NA, nsim, n_pars-13)
phiN.bias <- phiB.bias <- phiC.bias <- phiY.bias <- psiN.bias <- psiB.bias <- numeric()
pN.bias <- pBc.bias <- pBn.bias <- deltaY.bias <- deltaC.bias <- numeric()
gamma.bias <- alphaTy.bias <- alphaTc.bias <- kappaY.bias <- kappaC.bias <- eta.bias <- numeric()
omegaA.bias <- omegaB.bias <- numeric()
bias <- rmse <- matrix(NA, nsim, n_pars-13)

#model features
n.states <- 14
n.calfstates <- 37
n.obs <- 72

n.occasions <- 13
Ntot <- eval(parse(text = paste0('pars', scenario)))[12] 
marked <- rep(round(Ntot/n.occasions), n.occasions)
totrel <- sum(marked)


#simulate data
for (s in 1:nsim) {
     
     #Define params and store true values; #phiN, phiB, phiC, phiY, psiN, psiB, pN, pBn, pBc, deltaY, deltaC, N, gamma, alphaT, alphaC, kappaY, kappaC, omegaA, omegaB, eta, pi
     phiN <- phiNTRUE <- eval(parse(text = paste0('pars', scenario)))[1] #survival of NB adults
     phiB <- phiBTRUE <- eval(parse(text = paste0('pars', scenario)))[2] #survival of B adults
     phiC <- phiCTRUE <- eval(parse(text = paste0('pars', scenario)))[3] #survival of older calves (2+)
     phiYc <- phiYTRUE <- eval(parse(text = paste0('pars', scenario)))[4] #survival of younger young
     psiN <- psiNTRUE <- eval(parse(text = paste0('pars', scenario)))[5] #probability of Byoy first time
     psiB <- psiBTRUE <- eval(parse(text = paste0('pars', scenario)))[6] #probability of Byoy repeat
     
     pN <- pNTRUE <- eval(parse(text = paste0('pars', scenario)))[7] #probability of detecting adult
     pBc <- pBcTRUE <- eval(parse(text = paste0('pars', scenario)))[9] #probability of detecting adult w/ calf
     pBn <- pBnTRUE <- eval(parse(text = paste0('pars', scenario)))[8] #probability of detecting adult w/o calf
     deltaYc <- deltaYTRUE <- eval(parse(text = paste0('pars', scenario)))[10] #probability of detecting calf | p(A)
     deltaC <- deltaCTRUE <- eval(parse(text = paste0('pars', scenario)))[11] #probability of detecting calf | p(A)
     gamma <- gammaTRUE <- eval(parse(text = paste0('pars', scenario)))[13]
     alphaTy <- alphaTyTRUE <- eval(parse(text = paste0('pars', scenario)))[14]
     alphaTc <- alphaTcTRUE <- eval(parse(text = paste0('pars', scenario)))[15]
     kappaY <- kappaYTRUE <- eval(parse(text = paste0('pars', scenario)))[16]
     kappaC <- kappaCTRUE <- eval(parse(text = paste0('pars', scenario)))[17]
     omegaA <- omegaATRUE <- eval(parse(text = paste0('pars', scenario)))[18]
     omegaB <- omegaBTRUE <- eval(parse(text = paste0('pars', scenario)))[19]
     eta <- etaTRUE <- eval(parse(text = paste0('pars', scenario)))[20]
     
     pi[1:13] <- eval(parse(text = paste0('pars', scenario)))[21:33]
     
     true_vals <- c(phiNTRUE, phiBTRUE, phiCTRUE, phiYTRUE, psiNTRUE, psiBTRUE, pNTRUE, pBcTRUE, pBnTRUE,
                    deltaCTRUE, deltaYTRUE, gammaTRUE, alphaTyTRUE, alphaTcTRUE, kappaYTRUE,
                    kappaCTRUE, omegaATRUE, omegaBTRUE, etaTRUE)
     
     # Define matrices with survival, transition and recapture probabilities
     # These are 4-dimensional matrices
     # Dimension 1: state of departure
     # Dimension 2: state of arrival
     # Dimension 3: individual
     # Dimension 4: time
     
     # 1. State process 
     
     STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1))
     for (i in 1:totrel){
          for (t in 1:(n.occasions-1)){
               
               STATE[,,i,t] <- matrix(c( #state matrix 1: adult survival [14,14]
                    phiN,0,0,0,0,0,0,0,0,0,0,0,0,1-phiN, #non-breeder
                    0,phiB,0,0,0,0,0,0,0,0,0,0,0,1-phiB,
                    0,0,phiB,0,0,0,0,0,0,0,0,0,0,1-phiB,
                    0,0,0,phiB,0,0,0,0,0,0,0,0,0,1-phiB,
                    0,0,0,0,phiB,0,0,0,0,0,0,0,0,1-phiB,
                    0,0,0,0,0,phiB,0,0,0,0,0,0,0,1-phiB,
                    0,0,0,0,0,0,phiB,0,0,0,0,0,0,1-phiB,
                    0,0,0,0,0,0,0,phiB,0,0,0,0,0,1-phiB,
                    0,0,0,0,0,0,0,0,phiB,0,0,0,0,1-phiB,
                    0,0,0,0,0,0,0,0,0,phiB,0,0,0,1-phiB,
                    0,0,0,0,0,0,0,0,0,0,phiB,0,0,1-phiB,
                    0,0,0,0,0,0,0,0,0,0,0,phiB,0,1-phiB,
                    0,0,0,0,0,0,0,0,0,0,0,0,phiB,1-phiB,
                    0,0,0,0,0,0,0,0,0,0,0,0,0,1), nrow = n.states, byrow = T) %*%
                    
                    matrix(c( #state matrix 2: calf survival [14,37]
                         1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #
                         0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #
                         0,0,phiYc,1-phiYc,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #
                         0,0,0,0,phiYc,1-phiYc,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #
                         0,0,0,0,0,0,phiC,1-phiC,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #
                         0,0,0,0,0,0,0,0,phiC*phiYc,phiC*(1-phiYc),(1-phiC)*phiYc,(1-phiC)*(1-phiYc),0,
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #
                         0,0,0,0,0,0,0,0,0,0,0,0,phiC,1-phiC,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0,phiC*phiYc,phiC*(1-phiYc),(1-phiC)*phiYc,
                         (1-phiC)*(1-phiYc),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,phiC*phiYc,phiC*(1-phiYc),(1-phiC)*phiYc,
                         (1-phiYc)*(1-phiC),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,phiC,1-phiC,0,0,0,0,0,0,0,0,0,0,0,0,0, #
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,phiC*phiYc,phiC*(1-phiYc),
                         (1-phiC)*phiYc,(1-phiC)*(1-phiYc),0,0,0,0,0,0,0,0,0, #
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,phiC*phiYc,
                         phiC*(1-phiYc),phiYc*(1-phiC),(1-phiC)*(1-phiYc),0,0,0,0,0, #
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,phiC,(1-phiC),0, #
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1), #
                         nrow = n.states, byrow = T) %*%
                    
                    matrix(c( #state matrix 3: calf aging [37,14]
                         1,0,0,0,0,0,0,0,0,0,0,0,0,0, #
                         0,1,0,0,0,0,0,0,0,0,0,0,0,0, #
                         0,0,1,0,0,0,0,0,0,0,0,0,0,0, #
                         0,0,0,0,0,0,0,0,0,1,0,0,0,0, #
                         0,0,0,1,0,0,0,0,0,0,0,0,0,0, #
                         0,0,0,0,0,0,0,0,0,0,1,0,0,0, #
                         0,0,0,0,1,0,0,0,0,0,0,0,0,0, #
                         0,0,0,0,0,0,0,0,0,0,1,0,0,0, #
                         0,0,0,0,0,1,0,0,0,0,0,0,0,0, #
                         0,0,0,0,0,0,0,0,0,0,0,1,0,0, #
                         0,0,1,0,0,0,0,0,0,0,0,0,0,0, #
                         0,0,0,0,0,0,0,0,0,1,0,0,0,0, #
                         0,0,0,0,0,0,1,0,0,0,0,0,0,0, #
                         0,0,0,0,0,0,0,0,0,0,1,0,0,0, #
                         0,0,0,0,0,0,0,1,0,0,0,0,0,0, #
                         0,0,0,0,0,0,0,0,0,0,0,0,1,0, #
                         0,0,1,0,0,0,0,0,0,0,0,0,0,0, #
                         0,0,0,0,0,0,0,0,0,1,0,0,0,0, #
                         0,0,0,0,0,0,0,0,1,0,0,0,0,0, #
                         0,0,0,0,0,0,1,0,0,0,0,0,0,0, #
                         0,0,0,1,0,0,0,0,0,0,0,0,0,0, #
                         0,0,0,0,0,0,0,0,0,0,1,0,0,0, #
                         0,0,0,0,0,0,1,0,0,0,0,0,0,0, #
                         0,0,0,0,0,0,0,0,0,0,1,0,0,0, #
                         0,0,0,0,0,0,0,1,0,0,0,0,0,0, #
                         0,0,0,0,0,0,0,0,0,0,0,0,1,0, #
                         0,0,1,0,0,0,0,0,0,0,0,0,0,0, #
                         0,0,0,0,0,0,0,0,0,1,0,0,0,0, #
                         0,0,0,0,0,0,0,0,1,0,0,0,0,0, #
                         0,0,0,0,0,0,1,0,0,0,0,0,0,0, #
                         0,0,0,1,0,0,0,0,0,0,0,0,0,0, #
                         0,0,0,0,0,0,0,0,0,0,1,0,0,0, #
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0, #
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0, #
                         0,0,0,0,1,0,0,0,0,0,0,0,0,0, #
                         0,0,0,0,0,0,0,0,0,0,1,0,0,0, #
                         0,0,0,0,0,0,0,0,0,0,0,0,0,1), nrow = n.calfstates, byrow = T) %*%
                    
                    matrix(c( #state matrix 4: breeding probability [14,14]
                         1-psiN,0,psiN,0,0,0,0,0,0,0,0,0,0,0, #
                         0,1-psiB,psiB,0,0,0,0,0,0,0,0,0,0,0, #
                         0,0,0,1,0,0,0,0,0,0,0,0,0,0, #
                         0,0,0,0,1-psiB,psiB,0,0,0,0,0,0,0,0, #
                         0,0,0,0,0,0,1-psiB,psiB,0,0,0,0,0,0, #
                         0,0,0,0,0,0,0,0,1,0,0,0,0,0, #
                         0,0,0,0,0,0,0,0,0,1-psiB,psiB,0,0,0, #7
                         0,0,0,0,0,0,0,0,0,0,0,1,0,0, #
                         0,0,0,0,0,0,0,0,0,0,0,0,1,0, #
                         0,1,0,0,0,0,0,0,0,0,0,0,0,0, #
                         0,1-psiB,psiB,0,0,0,0,0,0,0,0,0,0,0, #11
                         0,0,0,0,0,0,1,0,0,0,0,0,0,0, #12
                         0,0,0,0,0,0,0,0,0,1,0,0,0,0, #13
                         0,0,0,0,0,0,0,0,0,0,0,0,0,1), #14
                         nrow = n.states, byrow = T)
          } #t
     } #i
     
     # 2.Observation process 
     OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
     for (i in 1:totrel){
          for (t in 1:(n.occasions-1)){
               
               OBS[,,i,t] <- matrix(c( #obs matrix 1: adult detection [14,14]
                    1-pN,pN,0,0,0,0,0,0,0,0,0,0,0,0,
                    1-pBn,0,pBn,0,0,0,0,0,0,0,0,0,0,0,
                    1-pBc,0,0,pBc,0,0,0,0,0,0,0,0,0,0,
                    1-pBc,0,0,0,pBc,0,0,0,0,0,0,0,0,0,
                    1-pBc,0,0,0,0,pBc,0,0,0,0,0,0,0,0,
                    1-pBc,0,0,0,0,0,pBc,0,0,0,0,0,0,0,
                    1-pBc,0,0,0,0,0,0,pBc,0,0,0,0,0,0,
                    1-pBc,0,0,0,0,0,0,0,pBc,0,0,0,0,0,
                    1-pBc,0,0,0,0,0,0,0,0,pBc,0,0,0,0,
                    1-pBc,0,0,0,0,0,0,0,0,0,pBc,0,0,0,
                    1-pBc,0,0,0,0,0,0,0,0,0,0,pBc,0,0,
                    1-pBc,0,0,0,0,0,0,0,0,0,0,0,pBc,0,
                    1-pBc,0,0,0,0,0,0,0,0,0,0,0,0,pBc,
                    1,0,0,0,0,0,0,0,0,0,0,0,0,0), nrow = n.states, byrow = T) %*%
                    matrix(c( #obs matrix 2: calf detection [14, 25]
                         1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
                         0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,1-deltaYc,deltaYc,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,1-deltaYc,0,deltaYc,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,1-deltaC,0,0,deltaC,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,(1-deltaYc)*(1-deltaC),0,0,0,deltaYc*(1-deltaC),deltaC*(1-deltaYc),
                         0,0,0,0,0,0,0,0,0,0,0,0,deltaYc*deltaC,0,0,0,0,0,
                         0,1-deltaC,0,0,0,0,0,deltaC,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,(1-deltaYc)*(1-deltaC),0,0,0,0,0,0,deltaYc*(1-deltaC),deltaC*(1-deltaYc),
                         0,0,0,0,0,0,0,0,0,0,deltaYc*deltaC,0,0,0,0,
                         0,(1-deltaYc)*(1-deltaC),0,0,0,0,0,0,0,0,deltaYc*(1-deltaC),deltaC*(1-deltaYc),
                         0,0,0,0,0,0,0,0,0,deltaYc*deltaC,0,0,0,
                         0,1-deltaC,0,0,0,0,0,0,0,0,0,0,deltaC,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,(1-deltaYc)*(1-deltaC),0,0,0,0,0,0,0,0,0,0,0,deltaYc*(1-deltaC),
                         deltaC*(1-deltaYc),0,0,0,0,0,0,0,deltaYc*deltaC,0,0,
                         0,(1-deltaYc)*(1-deltaC),0,0,0,0,0,0,0,0,0,0,0,0,0,deltaYc*(1-deltaC),
                         deltaC*(1-deltaYc),0,0,0,0,0,0,deltaYc*deltaC,0,
                         0,(1-deltaC)*(1-deltaC),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,deltaC-deltaC^2,
                         deltaC-deltaC^2,0,0,0,0,0,deltaC^2), nrow = n.states, byrow = T) %*%
                    matrix(c( #obs matrix 3: calf age determination [25, 72]
                         1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #
                         0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #
                         0,0,gamma*alphaTy,gamma*(1-alphaTy)*kappaY,0,0,gamma*(1-alphaTy)*(1-kappaY),0,0,0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #3
                         0,0,0,gamma*(1-alphaTc)*omegaA*eta,gamma*alphaTc,gamma*(1-alphaTc)*omegaA*(1-eta),gamma*(1-alphaTc)*(1-omegaA),0,0,0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #4
                         0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaA),gamma*(1-alphaTc)*omegaA*eta,gamma*alphaTc,gamma*(1-alphaTc)*omegaA*(1-eta),0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #5
                         0,0,gamma*alphaTy,gamma*(1-alphaTy)*kappaY,0,0,gamma*(1-alphaTy)*(1-kappaY),0,0,0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #6
                         0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaA),gamma*(1-alphaTc)*omegaA*eta,gamma*alphaTc,gamma*(1-alphaTc)*omegaA*(1-eta),0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #7
                         0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC),0,0,gamma*(1-alphaTc)*(1-omegaB)*kappaC,gamma*alphaTc,gamma*(1-alphaTc)*omegaB,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,0,0,0,0, #8
                         0,0,gamma*alphaTy,gamma*(1-alphaTy)*kappaY,0,0,gamma*(1-alphaTy)*(1-kappaY),0,0,0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #9
                         0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC),0,0,gamma*(1-alphaTc)*(1-omegaB)*kappaC,gamma*alphaTc,gamma*(1-alphaTc)*omegaB,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,0,0,0,0,0,0, #10
                         0,0,0,gamma*(1-alphaTc)*omegaA*eta,gamma*alphaTc,gamma*(1-alphaTc)*omegaA*(1-eta),gamma*(1-alphaTc)*(1-omegaA),0,0,0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #11
                         0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC),0,0,gamma*(1-alphaTc)*(1-omegaB)*kappaC,gamma*alphaTc,gamma*(1-alphaTc)*omegaB,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,0,0,0,0, #12
                         0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,gamma*(1-alphaTc)*(1-omegaB)*kappaC,gamma*alphaTc,gamma*(1-alphaTc)*omegaB,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #13
                         0,0,gamma*alphaTy,gamma*(1-alphaTy)*kappaY,0,0,gamma*(1-alphaTy)*(1-kappaY),0,0,0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #14
                         0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,gamma*(1-alphaTc)*(1-omegaB)*kappaC,gamma*alphaTc,gamma*(1-alphaTc)*omegaB,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #15
                         0,0,0,gamma*(1-alphaTc)*omegaA*eta,gamma*alphaTc,gamma*(1-alphaTc)*omegaA*(1-eta),gamma*(1-alphaTc)*(1-omegaA),0,0,0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #16
                         0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,gamma*(1-alphaTc)*(1-omegaB)*kappaC,gamma*alphaTc,gamma*(1-alphaTc)*omegaB,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #17
                         0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaA),gamma*(1-alphaTc)*omegaA*eta,gamma*alphaTc,gamma*(1-alphaTc)*omegaA*(1-eta),0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #18
                         0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,gamma*(1-alphaTc)*(1-omegaB)*kappaC,gamma*alphaTc,gamma*(1-alphaTc)*omegaB,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #19
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*alphaTy)*(gamma*(1-alphaTc)*(1-omegaA)),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*(1-omegaA)),0,(gamma*alphaTy)*(gamma*(1-alphaTc)*omegaA*eta),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*omegaA*eta),(gamma*alphaTy)*(gamma*alphaTc),(gamma*(1-alphaTy)*kappaY)*(gamma*alphaTc),(gamma*(1-alphaTy)*(1-kappaY))*(gamma*alphaTc),(gamma*alphaTy)*(gamma*(1-alphaTc)*omegaA*(1-eta)),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*omegaA*(1-eta)),0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*omegaA*(1-eta)),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*alphaTy)*(1-gamma),(gamma*(1-alphaTy)*kappaY)*(1-gamma),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaA)),((1-gamma)*(gamma*(1-alphaTc)*omegaA*eta))+((gamma*(1-alphaTy)*(1-kappaY))*(1-gamma)),(1-gamma)*(gamma*alphaTc),(1-gamma)*(gamma*(1-alphaTc)*omegaA*(1-eta)),0,0,0,0,(1-gamma)^2,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*(1-omegaA)),0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*omegaA*eta),0, #20
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*alphaTy)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),0,0,0,0,0,0,(gamma*alphaTy)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),0,(gamma*alphaTy)*(gamma*alphaTc),(gamma*(1-alphaTy)*kappaY)*(gamma*alphaTc),0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*alphaTc),(gamma*alphaTy)*(gamma*(1-alphaTc)*omegaB),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*omegaB),0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*omegaB),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*alphaTy)*(1-gamma),(gamma*(1-alphaTy)*kappaY)*(1-gamma),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),(gamma*(1-alphaTy)*(1-kappaY))*(1-gamma),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(1-gamma)*(gamma*alphaTc),(1-gamma)*(gamma*(1-alphaTc)*omegaB),0,0,(1-gamma)^2,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),0,0,0,0, #21
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),(gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),0,0,0,0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),0,0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*alphaTc),(gamma*alphaTc)*(gamma*alphaTc),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*alphaTc),(gamma*(1-alphaTc)*(1-omegaA))*(gamma*alphaTc),0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*omegaB),(gamma*alphaTc)*(gamma*(1-alphaTc)*omegaB),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*omegaB),(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*omegaB),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(1-gamma),(gamma*alphaTc)*(1-gamma),((1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)))+((gamma*(1-alphaTc)*omegaA*(1-eta))*(1-gamma)),(gamma*(1-alphaTc)*(1-omegaA))*(1-gamma),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(1-gamma)*(gamma*alphaTc),(1-gamma)*(gamma*(1-alphaTc)*omegaB),0,0,(1-gamma)^2,(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),0,0,0,0, #22
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*alphaTy)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,0,0,0,0,0,(gamma*alphaTy)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,0,0,0,0,0,(gamma*alphaTy)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),0,0,(gamma*alphaTy)*(gamma*alphaTc),(gamma*(1-alphaTy)*kappaY)*(gamma*alphaTc),0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*alphaTc),0,0,(gamma*alphaTy)*(gamma*(1-alphaTc)*omegaB),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*omegaB),0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*omegaB),0,0,(gamma*alphaTy)*(1-gamma),(gamma*(1-alphaTy)*kappaY)*(1-gamma),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*(1-alphaTy)*(1-kappaY))*(1-gamma),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(1-gamma)*(gamma*alphaTc),(1-gamma)*(gamma*(1-alphaTc)*omegaB),(1-gamma)^2,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,0,0,0, #23
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,0,0,0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,0,0,0,0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*alphaTc),(gamma*alphaTc)*(gamma*alphaTc),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*alphaTc),(gamma*(1-alphaTc)*(1-omegaA))*(gamma*alphaTc),0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*omegaB),(gamma*alphaTc)*(gamma*(1-alphaTc)*omegaB),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*omegaB),(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*omegaB),0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(1-gamma),(gamma*alphaTc)*(1-gamma),((1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5))+((gamma*(1-alphaTc)*omegaA*(1-eta))*(1-gamma)),(gamma*(1-alphaTc)*(1-omegaA))*(1-gamma),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(1-gamma)*(gamma*alphaTc),(1-gamma)*(gamma*(1-alphaTc)*omegaB),(1-gamma)^2,(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,0,0,0, #24
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,0,0,0,0,0,0,0,(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),0,0,0,(gamma*(1-alphaTc)*(1-omegaA))*(gamma*alphaTc),(gamma*(1-alphaTc)*omegaA*eta)*(gamma*alphaTc),(gamma*alphaTc)*(gamma*alphaTc),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*alphaTc),0,0,0,(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*omegaB),(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*omegaB),(gamma*alphaTc)*(gamma*(1-alphaTc)*omegaB),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*omegaB),0,0,0,((1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5))+((gamma*(1-alphaTc)*(1-omegaA))*(1-gamma)),(gamma*(1-alphaTc)*omegaA*eta)*(1-gamma),(gamma*alphaTc)*(1-gamma),((1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5))+((gamma*(1-alphaTc)*omegaA*(1-eta))*(1-gamma)),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(1-gamma)*(gamma*alphaTc),(1-gamma)*(gamma*(1-alphaTc)*omegaB),(1-gamma)^2,(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,
                         (gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5)),
                         nrow = 25, byrow = T)
          } #t
     } #i
     
     # 3. Initial state probabilities
     # probability of being in 1 of 13 states at first detection, can't be dead
     STATE.INIT <- array(NA, dim=c(n.states-1, totrel, n.occasions))
     for (i in 1:totrel){
          for (t in 1:n.occasions){
               STATE.INIT[1:13,i,t] <- pi # prob of being in any state for all individuals and occasions 
          }
     }
     
     # 4. Probability of initial observation at first capture | initial state
     # equivalent to po.init <- po1.init %*% po2 %*% po3 in model code
     OBS.INIT <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions)) # [14, 66, inds, t]
     
     for (i in 1:totrel){
          for (t in 1:n.occasions){
               OBS.INIT[,,i,t] <- matrix(c( #deterministic initial observation
                    0,1,0,0,0,0,0,0,0,0,0,0,0,0, 
                    0,0,1,0,0,0,0,0,0,0,0,0,0,0,
                    0,0,0,1,0,0,0,0,0,0,0,0,0,0,
                    0,0,0,0,1,0,0,0,0,0,0,0,0,0,
                    0,0,0,0,0,1,0,0,0,0,0,0,0,0,
                    0,0,0,0,0,0,1,0,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,1,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,1,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,0,1,0,0,0,0,
                    0,0,0,0,0,0,0,0,0,0,1,0,0,0,
                    0,0,0,0,0,0,0,0,0,0,0,1,0,0,
                    0,0,0,0,0,0,0,0,0,0,0,0,1,0,
                    0,0,0,0,0,0,0,0,0,0,0,0,0,1,
                    1,0,0,0,0,0,0,0,0,0,0,0,0,0),
                    nrow = n.states, byrow = TRUE) %*%
                    matrix(c( #obs matrix 2: calf detection [14, 25]
                         1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
                         0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,1-deltaYc,deltaYc,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,1-deltaYc,0,deltaYc,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,1-deltaC,0,0,deltaC,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,(1-deltaYc)*(1-deltaC),0,0,0,deltaYc*(1-deltaC),deltaC*(1-deltaYc),
                         0,0,0,0,0,0,0,0,0,0,0,0,deltaYc*deltaC,0,0,0,0,0,
                         0,1-deltaC,0,0,0,0,0,deltaC,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,(1-deltaYc)*(1-deltaC),0,0,0,0,0,0,deltaYc*(1-deltaC),deltaC*(1-deltaYc),
                         0,0,0,0,0,0,0,0,0,0,deltaYc*deltaC,0,0,0,0,
                         0,(1-deltaYc)*(1-deltaC),0,0,0,0,0,0,0,0,deltaYc*(1-deltaC),deltaC*(1-deltaYc),
                         0,0,0,0,0,0,0,0,0,deltaYc*deltaC,0,0,0,
                         0,1-deltaC,0,0,0,0,0,0,0,0,0,0,deltaC,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,(1-deltaYc)*(1-deltaC),0,0,0,0,0,0,0,0,0,0,0,deltaYc*(1-deltaC),
                         deltaC*(1-deltaYc),0,0,0,0,0,0,0,deltaYc*deltaC,0,0,
                         0,(1-deltaYc)*(1-deltaC),0,0,0,0,0,0,0,0,0,0,0,0,0,deltaYc*(1-deltaC),
                         deltaC*(1-deltaYc),0,0,0,0,0,0,deltaYc*deltaC,0,
                         0,(1-deltaC)*(1-deltaC),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,deltaC-deltaC^2,
                         deltaC-deltaC^2,0,0,0,0,0,deltaC^2),
                         nrow = n.states, byrow = T) %*%
                    matrix(c( #obs matrix 3: calf age determination [25, 72]
                         1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #
                         0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #
                         0,0,gamma*alphaTy,gamma*(1-alphaTy)*kappaY,0,0,gamma*(1-alphaTy)*(1-kappaY),0,0,0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #3
                         0,0,0,gamma*(1-alphaTc)*omegaA*eta,gamma*alphaTc,gamma*(1-alphaTc)*omegaA*(1-eta),gamma*(1-alphaTc)*(1-omegaA),0,0,0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #4
                         0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaA),gamma*(1-alphaTc)*omegaA*eta,gamma*alphaTc,gamma*(1-alphaTc)*omegaA*(1-eta),0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #5
                         0,0,gamma*alphaTy,gamma*(1-alphaTy)*kappaY,0,0,gamma*(1-alphaTy)*(1-kappaY),0,0,0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #6
                         0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaA),gamma*(1-alphaTc)*omegaA*eta,gamma*alphaTc,gamma*(1-alphaTc)*omegaA*(1-eta),0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #7
                         0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC),0,0,gamma*(1-alphaTc)*(1-omegaB)*kappaC,gamma*alphaTc,gamma*(1-alphaTc)*omegaB,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,0,0,0,0, #8
                         0,0,gamma*alphaTy,gamma*(1-alphaTy)*kappaY,0,0,gamma*(1-alphaTy)*(1-kappaY),0,0,0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #9
                         0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC),0,0,gamma*(1-alphaTc)*(1-omegaB)*kappaC,gamma*alphaTc,gamma*(1-alphaTc)*omegaB,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,0,0,0,0,0,0, #10
                         0,0,0,gamma*(1-alphaTc)*omegaA*eta,gamma*alphaTc,gamma*(1-alphaTc)*omegaA*(1-eta),gamma*(1-alphaTc)*(1-omegaA),0,0,0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #11
                         0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC),0,0,gamma*(1-alphaTc)*(1-omegaB)*kappaC,gamma*alphaTc,gamma*(1-alphaTc)*omegaB,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         0,0,0,0,0, #12
                         0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,gamma*(1-alphaTc)*(1-omegaB)*kappaC,gamma*alphaTc,gamma*(1-alphaTc)*omegaB,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #13
                         0,0,gamma*alphaTy,gamma*(1-alphaTy)*kappaY,0,0,gamma*(1-alphaTy)*(1-kappaY),0,0,0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #14
                         0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,gamma*(1-alphaTc)*(1-omegaB)*kappaC,gamma*alphaTc,gamma*(1-alphaTc)*omegaB,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #15
                         0,0,0,gamma*(1-alphaTc)*omegaA*eta,gamma*alphaTc,gamma*(1-alphaTc)*omegaA*(1-eta),gamma*(1-alphaTc)*(1-omegaA),0,0,0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #16
                         0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,gamma*(1-alphaTc)*(1-omegaB)*kappaC,gamma*alphaTc,gamma*(1-alphaTc)*omegaB,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #17
                         0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaA),gamma*(1-alphaTc)*omegaA*eta,gamma*alphaTc,gamma*(1-alphaTc)*omegaA*(1-eta),0,0,0,0,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #18
                         0,0,0,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,0,gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5,0,gamma*(1-alphaTc)*(1-omegaB)*kappaC,gamma*alphaTc,gamma*(1-alphaTc)*omegaB,1-gamma,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, #19
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*alphaTy)*(gamma*(1-alphaTc)*(1-omegaA)),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*(1-omegaA)),0,(gamma*alphaTy)*(gamma*(1-alphaTc)*omegaA*eta),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*omegaA*eta),(gamma*alphaTy)*(gamma*alphaTc),(gamma*(1-alphaTy)*kappaY)*(gamma*alphaTc),(gamma*(1-alphaTy)*(1-kappaY))*(gamma*alphaTc),(gamma*alphaTy)*(gamma*(1-alphaTc)*omegaA*(1-eta)),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*omegaA*(1-eta)),0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*omegaA*(1-eta)),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*alphaTy)*(1-gamma),(gamma*(1-alphaTy)*kappaY)*(1-gamma),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaA)),((1-gamma)*(gamma*(1-alphaTc)*omegaA*eta))+((gamma*(1-alphaTy)*(1-kappaY))*(1-gamma)),(1-gamma)*(gamma*alphaTc),(1-gamma)*(gamma*(1-alphaTc)*omegaA*(1-eta)),0,0,0,0,(1-gamma)^2,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*(1-omegaA)),0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*omegaA*eta),0, #20
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*alphaTy)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),0,0,0,0,0,0,(gamma*alphaTy)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),0,(gamma*alphaTy)*(gamma*alphaTc),(gamma*(1-alphaTy)*kappaY)*(gamma*alphaTc),0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*alphaTc),(gamma*alphaTy)*(gamma*(1-alphaTc)*omegaB),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*omegaB),0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*omegaB),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*alphaTy)*(1-gamma),(gamma*(1-alphaTy)*kappaY)*(1-gamma),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),(gamma*(1-alphaTy)*(1-kappaY))*(1-gamma),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(1-gamma)*(gamma*alphaTc),(1-gamma)*(gamma*(1-alphaTc)*omegaB),0,0,(1-gamma)^2,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),0,0,0,0, #21
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),(gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),0,0,0,0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),0,0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*alphaTc),(gamma*alphaTc)*(gamma*alphaTc),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*alphaTc),(gamma*(1-alphaTc)*(1-omegaA))*(gamma*alphaTc),0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*omegaB),(gamma*alphaTc)*(gamma*(1-alphaTc)*omegaB),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*omegaB),(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*omegaB),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(1-gamma),(gamma*alphaTc)*(1-gamma),((1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)))+((gamma*(1-alphaTc)*omegaA*(1-eta))*(1-gamma)),(gamma*(1-alphaTc)*(1-omegaA))*(1-gamma),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(1-gamma)*(gamma*alphaTc),(1-gamma)*(gamma*(1-alphaTc)*omegaB),0,0,(1-gamma)^2,(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)),0,0,0,0, #22
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*alphaTy)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,0,0,0,0,0,(gamma*alphaTy)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,0,0,0,0,0,(gamma*alphaTy)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),0,0,(gamma*alphaTy)*(gamma*alphaTc),(gamma*(1-alphaTy)*kappaY)*(gamma*alphaTc),0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*alphaTc),0,0,(gamma*alphaTy)*(gamma*(1-alphaTc)*omegaB),(gamma*(1-alphaTy)*kappaY)*(gamma*(1-alphaTc)*omegaB),0,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*omegaB),0,0,(gamma*alphaTy)*(1-gamma),(gamma*(1-alphaTy)*kappaY)*(1-gamma),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*(1-alphaTy)*(1-kappaY))*(1-gamma),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(1-gamma)*(gamma*alphaTc),(1-gamma)*(gamma*(1-alphaTc)*omegaB),(1-gamma)^2,0,(gamma*(1-alphaTy)*(1-kappaY))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,0,0,0, #23
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,0,0,0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,0,0,0,0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*alphaTc),(gamma*alphaTc)*(gamma*alphaTc),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*alphaTc),(gamma*(1-alphaTc)*(1-omegaA))*(gamma*alphaTc),0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*omegaB),(gamma*alphaTc)*(gamma*(1-alphaTc)*omegaB),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*omegaB),(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*omegaB),0,0,0,(gamma*(1-alphaTc)*omegaA*eta)*(1-gamma),(gamma*alphaTc)*(1-gamma),((1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5))+((gamma*(1-alphaTc)*omegaA*(1-eta))*(1-gamma)),(gamma*(1-alphaTc)*(1-omegaA))*(1-gamma),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(1-gamma)*(gamma*alphaTc),(1-gamma)*(gamma*(1-alphaTc)*omegaB),(1-gamma)^2,(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,0,0,0, #24
                         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,0,0,0,0,0,0,0,(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),0,0,0,(gamma*(1-alphaTc)*(1-omegaA))*(gamma*alphaTc),(gamma*(1-alphaTc)*omegaA*eta)*(gamma*alphaTc),(gamma*alphaTc)*(gamma*alphaTc),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*alphaTc),0,0,0,(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*omegaB),(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*omegaB),(gamma*alphaTc)*(gamma*(1-alphaTc)*omegaB),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*omegaB),0,0,0,((1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5))+((gamma*(1-alphaTc)*(1-omegaA))*(1-gamma)),(gamma*(1-alphaTc)*omegaA*eta)*(1-gamma),(gamma*alphaTc)*(1-gamma),((1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5))+((gamma*(1-alphaTc)*omegaA*(1-eta))*(1-gamma)),0,(1-gamma)*(gamma*(1-alphaTc)*(1-omegaB)*kappaC),(1-gamma)*(gamma*alphaTc),(1-gamma)*(gamma*(1-alphaTc)*omegaB),(1-gamma)^2,(gamma*(1-alphaTc)*(1-omegaA))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*(1-alphaTc)*omegaA*eta)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*alphaTc)*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),(gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5),0,
                         (gamma*(1-alphaTc)*omegaA*(1-eta))*(gamma*(1-alphaTc)*(1-omegaB)*(1-kappaC)*0.5)),
                         nrow = 25, byrow = T)
          } #t
     } #i
     
     # Execute function
     sim <- ME.sim(STATE, OBS, STATE.INIT, OBS.INIT, marked)
     CH <- sim$CH
     CH.TRUE <- sim$CH.TRUE
     
     # Compute vector with occasion of first capture
     get.first <- function(x) min(which(x!=0))
     fc <- apply(CH, 1, get.first)
     
     # Recode CH matrix after get fc
     #1 = not seen
     rCH <- CH          # Recoded CH
     rCH[is.na(rCH)] <- 1
     CH <- rCH
     
     ##### run model ######
     
     ## Bundle data
     jags.data = list(N=totrel, Years=n.occasions, y=as.matrix(CH), First=fc)
     
     ## Initial values
     inits <- function() list(z = CH.TRUE,
                              gamma = gammaTRUE, alphaTy = alphaTyTRUE,
                              deltaC=runif(1,0.2,0.8))
     
     ## Specify the parameters to be monitored
     params <- c('phiN', 'phiB', 'phiC', 'phiYc', 'pN', 'pBn', 'pBc',
                 'deltaC', 'deltaYc', 'psiN', 'psiB', 'gamma',  'alphaTy', 'alphaTc', 
                 'kappaY', 'kappaC', 'omegaA', 'omegaB', 'eta', 'pi')
     
     #Set up model run parameters
     nc <- 3; ni <- 40000; nbi <- 20000; nt <- 4
     
     # Fit the model
     out <- jags(jags.data, inits=inits, params, model.file, 
                 n.chains=nc, n.burnin=nbi, n.thin=nt, parallel=TRUE, n.cores=3,
                 n.iter=ni)
     
     max_rhat <- max(unlist(out$Rhat))
     out_median[s,] <- c(unlist(out$mean[1:19]), max_rhat)
     colnames(out_median) <- c(names(out$mean[1:19]), 'rhat')
     
     out_true[s,] <- true_vals
     colnames(out_true) <- c(parameters[1:19])
     
     #rel bias: (est-T)/T
     phiN.bias <- (out$mean$phiN-phiNTRUE)/phiNTRUE
     phiB.bias <- (out$mean$phiB-phiBTRUE)/phiBTRUE
     phiC.bias <- (out$mean$phiC-phiCTRUE)/phiCTRUE
     phiY.bias <- (out$mean$phiYc-phiYTRUE)/phiYTRUE
     pN.bias <- (out$mean$pN-pNTRUE)/pNTRUE
     pBn.bias <- (out$mean$pBn-pBnTRUE)/pBnTRUE
     pBc.bias <- (out$mean$pBc-pBcTRUE)/pBcTRUE
     deltaC.bias <- (out$mean$deltaC-deltaCTRUE)/deltaCTRUE
     deltaY.bias <- (out$mean$deltaYc-deltaYTRUE)/deltaYTRUE
     psiN.bias <- (out$mean$psiN-psiNTRUE)/psiNTRUE
     psiB.bias <- (out$mean$psiB-psiBTRUE)/psiBTRUE
     gamma.bias <- (out$mean$gamma-gammaTRUE)/gammaTRUE
     alphaTy.bias <- (out$mean$alphaTy-alphaTyTRUE)/alphaTyTRUE
     alphaTc.bias <- (out$mean$alphaTc-alphaTcTRUE)/alphaTcTRUE
     kappaY.bias <- (out$mean$kappaY-kappaYTRUE)/kappaYTRUE
     kappaC.bias <- (out$mean$kappaC-kappaCTRUE)/kappaCTRUE
     omegaA.bias <- (out$mean$omegaA-omegaATRUE)/omegaATRUE
     omegaB.bias <- (out$mean$omegaB-omegaBTRUE)/omegaBTRUE
     eta.bias <- (out$mean$eta-etaTRUE)/etaTRUE
     
     bias[s,] <- c(phiN.bias, phiB.bias, phiC.bias, phiY.bias, pN.bias, pBn.bias, pBc.bias,
                   deltaC.bias, deltaY.bias, psiN.bias, psiB.bias, gamma.bias,
                   alphaTy.bias, alphaTc.bias, kappaY.bias, kappaC.bias, omegaA.bias, omegaB.bias, eta.bias)
     colnames(bias) <- c('phiN', 'phiB', 'phiC', 'phiY', 'pN', 'pBn', 'pBc',
                         'deltaC', 'deltaY', 'psiN', 'psiB', 'gamma',  'alphaTy', 'alphaTc', 
                         'kappaY', 'kappaC', 'omegaA', 'omegaB', 'eta')
     
     #rmse: sqrt(mean(theta.iter-T)^2))
     phiB.rmse <- sqrt(mean((out$sims.list$phiB-phiBTRUE)^2))
     phiN.rmse <- sqrt(mean((out$sims.list$phiN-phiNTRUE)^2))
     phiC.rmse <- sqrt(mean((out$sims.list$phiC-phiCTRUE)^2))
     phiY.rmse <- sqrt(mean((out$sims.list$phiYc-phiYTRUE)^2))
     pN.rmse <- sqrt(mean((out$sims.list$pN-pNTRUE)^2))
     pBn.rmse <- sqrt(mean((out$sims.list$pBn-pBnTRUE)^2))
     pBc.rmse <- sqrt(mean((out$sims.list$pBc-pBcTRUE)^2))
     deltaC.rmse <- sqrt(mean((out$sims.list$deltaC-deltaCTRUE)^2))
     deltaY.rmse <- sqrt(mean((out$sims.list$deltaYc-deltaYTRUE)^2))
     psiN.rmse <- sqrt(mean((out$sims.list$psiN-psiNTRUE)^2))
     psiB.rmse <- sqrt(mean((out$sims.list$psiB-psiBTRUE)^2))
     gamma.rmse <- sqrt(mean((out$sims.list$gamma-gammaTRUE)^2))
     alphaTy.rmse <- sqrt(mean((out$sims.list$alphaTy-alphaTyTRUE)^2))
     alphaTc.rmse <- sqrt(mean((out$sims.list$alphaTc-alphaTcTRUE)^2))
     kappaY.rmse <- sqrt(mean((out$sims.list$kappaY-kappaYTRUE)^2))
     kappaC.rmse <- sqrt(mean((out$sims.list$kappaC-kappaCTRUE)^2))
     omegaA.rmse <- sqrt(mean((out$sims.list$omegaA-omegaATRUE)^2))
     omegaB.rmse <- sqrt(mean((out$sims.list$omegaB-omegaBTRUE)^2))
     eta.rmse <- sqrt(mean((out$sims.list$eta-etaTRUE)^2))
     
     rmse[s,] <- c(phiN.rmse, phiB.rmse, phiC.rmse, phiY.rmse, pN.rmse, pBn.rmse, pBc.rmse,
                   deltaC.rmse, deltaY.rmse, psiN.rmse, psiB.rmse, gamma.rmse, 
                   alphaTy.rmse, alphaTc.rmse, 
                   kappaY.rmse, kappaC.rmse, omegaA.rmse, omegaB.rmse, eta.rmse)
     colnames(rmse) <- c('phiN', 'phiB', 'phiC', 'phiY', 'pN', 'pBn', 'pBc',
                         'deltaC', 'deltaY', 'psiN', 'psiB', 'gamma',  'alphaTy', 'alphaTc', 
                         'kappaY', 'kappaC', 'omegaA', 'omegaB', 'eta')
     
     #cri coverage
     out_cri[s,] <- out_true[s,parameters[1:19]] >= unlist(out$q2.5[parameters[1:19]]) &
          out_true[s,parameters[1:19]] <= unlist(out$q97.5[parameters[1:19]])
     
     colnames(out_cri) <- c('phiN', 'phiB', 'phiC', 'phiY', 'pN', 'pBn', 'pBc',
                            'deltaC', 'deltaY', 'psiN', 'psiB', 'gamma', 'alphaTy', 'alphaTc', 
                            'kappaY', 'kappaC', 'omegaA', 'omegaB', 'eta')
     
} #sim

#save output
# saveRDS(bias, file = here::here('output', 'Simulations', 'v4.8',
#                                 paste0("S", scenario, '_', 'bias.rds')))
# saveRDS(rmse, file = here::here('output', 'Simulations', 'v4.8',
#                                 paste0("S", scenario,'_', 'rmse.rds')))
# saveRDS(out_cri, file = here::here('output', 'Simulations', 'v4.8',
#                                    paste0("S", scenario, '_', 'cri.rds')))
# saveRDS(out_median, file = here::here('output', 'Simulations', 'v4.8',
#                                       paste0("S", scenario, '_', 'median.rds')))
# saveRDS(out, file = here::here('output', 'Simulations', 'v4.8', paste0("S", scenario,'_', 'chains.rds')))


#r diagnostics

#check convergence with max_rhat in out_median
# out_median <- readRDS(file = here::here('output', 'Simulations', 
#                                         'v4.8', 'S.4.6.25median.rds'))
# max(out_median[,'rhat'],na.rm = T)


