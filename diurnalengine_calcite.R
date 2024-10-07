#Last Modified on Dec 21 2023

#Authors: Drs. Lizzy Trower and Cedric Hagen, 
#         University of Colorado, Boulder   

setwd('~/Desktop/DiurnalEngine')

# gas exchange calculations ------------------------------------------------------

tempC <- 25 #{°C}

tempK <- tempC + 273 #{K}

K0 <- exp(-58.0931 + 90.5697*(100/tempK) + 22.294*log(tempK/100) + 35*(0.027766 - 0.025888*(tempK/100) + 0.005078*(tempK/100)^2)) #{mol/m^3/atm}
K0 = K0 *10^3

Sc = 2116.8 - 136.25*tempC + 4.7353*(tempC^2) - 0.092307*(tempC^3) + 0.0007555*(tempC^4) #{dimensionless}

# windspeed ------------------------------------------------------

u <- 4 #{m/s}

# ambient atmospheric pCO2 ------------------------------------------------------
pCO2_atm <- 1000 #{uatm} [I SET THIS AT 1000 AFTER A VERY QUICK GOOGLE SCHOLAR SEARCH; ADJUST AS NEEDED]

kCO2 <- 0.251*(u^2)*(Sc/660)^0.5 #{cm/hr}
kCO2 <- kCO2/100 #{m/hr}

waterdepth <- 1 #{m}
waterdensity <- 1000 #{kg/m^3} [AGAIN, ADJUST AS NEEDED]

kappa_p <- 150 #{umol/kg}

# define shape of forcing as a sine function ------------------------------------------------------

t_hr <- seq(from = 0, to = 24, by = 0.1)
period <- 24 #{hr}
offset <- 5 #{hr} when the sine curve will cross 0
kappa_p_factor <- kappa_p/(period/pi) #{umol/kg/hr}
photo <- kappa_p_factor*sin((2*pi/period)*(t_hr-offset)) #{umol/kg/hr}

# calciteinterp ------------------------------------------------------

calciteinterp <- function(tempC) {
  library(pracma)
  
  caln <- c(0.6, 1.9, 2.3)
  calk <- c(14.0, 3.9, 3.7)
  calT <- c(5, 25, 37)
  
  if (tempC <= 37){
    if (tempC >= 5){
      k_BR <- interp1(calT, calk, tempC, method = "linear")
      n_BR <- interp1(calT, caln, tempC, method = "linear")
    } else if (tempC < 5){
      k_BR <- 14.0
      n_BR <- 0.6
    }
  } else if (tempC > 37){
    k_BR <- 3.7
    n_BR <- 2.3
  }
  
  output <- list(k_BR,n_BR)
  
  return(output)
}

# precipitation kinetics ------------------------------------------------------

#k_rate <- 9*(10^-9) #{mol/m^2/s} [k_rate replaced with K_BR in this iteration]
#k_rate <- k_rate*60*60 #{mol/m^2/hr}

precipkin <- calciteinterp(tempC)

k_BR <- unlist(precipkin[1])
n_BR <- unlist(precipkin[2])

# dissolution kinetics - Walter and Morse 1985 (there are more recent data that we could also use / compare with) ------------------------------------------------------

n_diss <- 2.86 # [Mean Low-Mg Calcite n from WM'85]
k_diss <- 2239 #{umol/m^2/hr) [Mean Low-Mg Calcite k from WM'85, T8]

# carbon isotope fractionations ------------------------------------------------------

eps_DIC_cal <- 1.0
eps_g_DIC <- -2
eps_DIC_g <- -10.3
eps_DIC_org <- -19

# set up major ion concentrations ------------------------------------------------------

cCa <- 9.8 #{mmol/kg}
cMg <- 2.5 #{mmol/kg}
cNa <- 1.2 #{mmol/kg}
cK <- 0.1 #{mmol/kg}
cCl <- 0.6 #{mmol/kg}
cSO4 <- 9.8 #{mmol/kg}

#[PARAMETERS SET AT GREEN LAKE VALUES FROM HANNA/MIQUELA]

# set up carbonate chemistry variables ------------------------------------------------------

DIC_model <- zeros(1,length(t_hr))
#DIC_model[1] <- 3700 #{umol/kg}
Alk_model <- zeros(1,length(t_hr))
#Alk_model[1] <- 4600 #{umol/kg}

pH_model <- zeros(1,length(t_hr))
pH_model[1] <- 7.3
pCO2_model <- zeros(1,length(t_hr))
pCO2_model[1] <- 1000 #{uatm}
Omega_cal_model <- zeros(1,length(t_hr))
Fcarb_model <- zeros(1,length(t_hr))
Fgas_model <- zeros(1,length(t_hr))
d13C_DIC_model <- zeros(1,length(t_hr))
d13C_DIC_model[1] <- -0.7
d13C_org_model <- -8*ones(1,length(t_hr)) 

# [mean summer pH from Green Lake is about 7.3, so I used that. Adjust as needed.]
# [pco2 at 1000 assuming equilibrium with the value I found. Adjust as needed.]

# run phreeqc to get initial values ------------------------------------------------------

library(phreeqc)                       # Load the PHREEQC library
phrLoadDatabaseString(phreeqc.dat)     # Use the phreeqc database
phrSetOutputStringsOn(TRUE)            # Format the output as a character string 

co2_val <- log10(pCO2_model[1]/10^6)

input <- c(                            # Text string defining composition, units, and temperature.
  '  SOLUTION              '                                            ,   
  '  units         mmol/kgw'                                            ,
  paste('  temp              ',as.character(tempC)),
  paste('  pH                ',as.character(pH_model[1])),
  paste('  C(4)  1   CO2(g)  ',as.character(co2_val)),              
  paste('  Ca                ',as.character(cCa)),
  paste('  Mg                ',as.character(cMg)),
  paste('  Na                ',as.character(cNa)),
  paste('  K                 ',as.character(cK)),
  paste('  Cl                ',as.character(cCl)),
  paste('  S(6)              ',as.character(cSO4)),'
     SELECTED_OUTPUT       ',
  '  -high precision   TRUE',
  '  -pH               TRUE',
  '  -si               calcite',
  '  -Alkalinity       TRUE',
  ' -totals            C(4)')

phrRunString(input)                    # Run the input string
#output <- phrGetOutputStrings()        # Save the results in 'output'
output <- phrGetSelectedOutput()

Omega_cal_model[1] <- 10^output$n1$si_calcite
#pH_model[1] <- output$n1$pH
#pCO2_model[1] <- (10^output$n1$si_CO2.g.)*(10^6)

DIC_model[1] <- output$n1$C.4..mol.kgw.*10^6 #{umol/kg}
Alk_model[1] <- output$n1$Alk.eq.kgw.*10^6 #{ueq/kg}


Fcarb_model[1] <- k_BR*(Omega_cal_model[1] - 1)^n_BR/waterdensity/waterdepth #{umol/kg/hr}
Fgas_model[1] <- kCO2*K0*(pCO2_model[1] - pCO2_atm)/waterdensity/waterdepth #{umol/kg/hr}

#set mean d13C_org that's used for remineralized OM
d13C_org_mean <- -20 #{permil}

#set mean d13C_carb that's used for dissolving carbonate
d13C_carb_mean <- 4 #{permil}

# run diurnal engine ------------------------------------------------------
for (nn in 2:length(t_hr)){
  delta_t_model = t_hr[nn] - t_hr[nn-1] #{hr}
  if (Omega_cal_model[nn-1] == 1){
    Fcarb_model[nn] <- 0
  } else if (Omega_cal_model[nn-1] > 1){
    Fcarb_model[nn] <- k_BR*(Omega_cal_model[nn-1]-1)^n_BR/waterdensity/waterdepth #{umol/kg/hr}
  } else if (Omega_cal_model[nn-1] < 1){
    Fcarb_model[nn] <- -k_diss*(1 - Omega_cal_model[nn-1])^n_diss/waterdensity/waterdepth #{umol/kg/hr}
  }
  Fgas_model[nn] <- kCO2*K0*(pCO2_model[nn-1] - pCO2_atm)/waterdensity/waterdepth #{umol/kg/hr}
  DIC_model[nn] <- DIC_model[nn-1] - (photo[nn] + (Fcarb_model[nn] + Fgas_model[nn]))*delta_t_model
  Alk_model[nn] <- Alk_model[nn-1] - 2*Fcarb_model[nn]*delta_t_model
  
  input <- c(                            
    '  SOLUTION              '                                            ,   
    '  units         mmol/kgw'                                            ,
    paste('  temp              ',as.character(tempC)),
    paste('  Alkalinity        ',as.character(Alk_model[nn]/1000)),
    paste('  C(4)              ',as.character(DIC_model[nn]/1000)),              
    paste('  Ca                ',as.character(cCa)),
    paste('  Mg                ',as.character(cMg)),
    paste('  Na                ',as.character(cNa)),
    paste('  K                 ',as.character(cK)),
    paste('  Cl                ',as.character(cCl)),
    paste('  S(6)              ',as.character(cSO4)),'
     SELECTED_OUTPUT       ',
    '  -high precision   TRUE',
    '  -pH               TRUE',
    '  -si               calcite',
    '  -si               CO2(g)')
  
  phrRunString(input)                    # Run the input string
  #output <- phrGetOutputStrings();
  output <- phrGetSelectedOutput()
  
  Omega_cal_model[nn] <- 10^output$n1$si_calcite
  pH_model[nn] <- output$n1$pH
  pCO2_model[nn] <- (10^output$n1$si_CO2.g.)*(10^6)
  
  if (Fcarb_model[nn] >= 0){
    if (Fgas_model[nn] >= 0 & photo[nn] >= 0){
      d13C_DIC_model[nn] <- (d13C_DIC_model[nn-1]*DIC_model[nn-1] - (d13C_DIC_model[nn-1] + eps_DIC_cal)*Fcarb_model[nn]*delta_t_model - (d13C_DIC_model[nn-1] + eps_DIC_g)*Fgas_model[nn]*delta_t_model - (d13C_DIC_model[nn-1] + eps_DIC_org)*photo[nn]*delta_t_model)/DIC_model[nn]
      d13C_org_model[nn] <- (d13C_DIC_model[nn-1] + eps_DIC_org)*photo[nn]*delta_t_model
    } else if (Fgas_model[nn] >= 0 & photo[nn] < 0){
      d13C_DIC_model[nn] <- (d13C_DIC_model[nn-1]*DIC_model[nn-1] - (d13C_DIC_model[nn-1] + eps_DIC_cal)*Fcarb_model[nn]*delta_t_model - (d13C_DIC_model[nn-1] + eps_DIC_g)*Fgas_model[nn]*delta_t_model - d13C_org_mean*photo[nn]*delta_t_model)/DIC_model[nn]
    } else if (Fgas_model[nn] < 0 & photo[nn] >= 0){
      d13C_DIC_model[nn] <- (d13C_DIC_model[nn-1]*DIC_model[nn-1] - (d13C_DIC_model[nn-1] + eps_DIC_cal)*Fcarb_model[nn]*delta_t_model - (d13C_DIC_model[nn-1] + eps_g_DIC)*Fgas_model[nn]*delta_t_model - (d13C_DIC_model[nn-1] + eps_DIC_org)*photo[nn]*delta_t_model)/DIC_model[nn]
      d13C_org_model[nn] <- (d13C_DIC_model[nn-1] + eps_DIC_org)*photo[nn]*delta_t_model
    } else if (Fgas_model[nn] < 0 & photo[nn] < 0){
      d13C_DIC_model[nn] <- (d13C_DIC_model[nn-1]*DIC_model[nn-1] - (d13C_DIC_model[nn-1] + eps_DIC_cal)*Fcarb_model[nn]*delta_t_model - (d13C_DIC_model[nn-1] + eps_g_DIC)*Fgas_model[nn]*delta_t_model - d13C_org_mean*photo[nn]*delta_t_model)/DIC_model[nn]
    }
  } else if (Fcarb_model[nn] < 0){
    if (Fgas_model[nn] >= 0 & photo[nn] >= 0){
      d13C_DIC_model[nn] <- (d13C_DIC_model[nn-1]*DIC_model[nn-1] - d13C_carb_mean*Fcarb_model[nn]*delta_t_model - (d13C_DIC_model[nn-1] + eps_DIC_g)*Fgas_model[nn]*delta_t_model - (d13C_DIC_model[nn-1] + eps_DIC_org)*photo[nn]*delta_t_model)/DIC_model[nn]
      d13C_org_model[nn] <- (d13C_DIC_model[nn-1] + eps_DIC_org)*photo[nn]*delta_t_model
    } else if (Fgas_model[nn] >= 0 & photo[nn] < 0){
      d13C_DIC_model[nn] <- (d13C_DIC_model[nn-1]*DIC_model[nn-1] - d13C_carb_mean*Fcarb_model[nn]*delta_t_model - (d13C_DIC_model[nn-1] + eps_DIC_g)*Fgas_model[nn]*delta_t_model - d13C_org_mean*photo[nn]*delta_t_model)/DIC_model[nn]
    } else if (Fgas_model[nn] < 0 & photo[nn] >= 0){
      d13C_DIC_model[nn] <- (d13C_DIC_model[nn-1]*DIC_model[nn-1] - d13C_carb_mean*Fcarb_model[nn]*delta_t_model - (d13C_DIC_model[nn-1] + eps_g_DIC)*Fgas_model[nn]*delta_t_model - (d13C_DIC_model[nn-1] + eps_DIC_org)*photo[nn]*delta_t_model)/DIC_model[nn]
      d13C_org_model[nn] <- (d13C_DIC_model[nn-1] + eps_DIC_org)*photo[nn]*delta_t_model
    } else if (Fgas_model[nn] < 0 & photo[nn] < 0){
      d13C_DIC_model[nn] <- (d13C_DIC_model[nn-1]*DIC_model[nn-1] - d13C_carb_mean*Fcarb_model[nn]*delta_t_model - (d13C_DIC_model[nn-1] + eps_g_DIC)*Fgas_model[nn]*delta_t_model - d13C_org_mean*photo[nn]*delta_t_model)/DIC_model[nn]
    }
  }
}

# make figure with data and model fits ------------------------------------------------------

pdf("DE_plots_calcite_model.pdf") 
out_plots <- par(mfrow=c(3,2)) 
out_plots <- plot(t_hr, DIC_model/1000,type='l',lty=2,lwd=2,xlab = "Hour of Day",ylab = "[DIC] (mmol/kg)",xlim=c(0,24))
out_plots <- plot(t_hr, Alk_model/1000,type='l',lty=2,lwd=2,xlab = "Hour of Day",ylab = "Alk (mmol/kg)",xlim=c(0,24))
out_plots <- plot(t_hr, pH_model,type='l',lty=2,lwd=2,xlab = "Hour of Day",ylab = "pH",xlim=c(0,24))
out_plots <- plot(t_hr, pCO2_model,type='l',lty=2,lwd=2,xlab = "Hour of Day",ylab = expression(paste('pCO'[2]*' (',mu,'atm)',sep='')),xlim=c(0,24))
out_plots <- plot(t_hr, Omega_cal_model,type='l',lty=2,lwd=2,xlab = "Hour of Day",ylab = expression(Omega[calc]),xlim=c(0,24))
out_plots <- plot(t_hr, d13C_DIC_model,type='l',lty=2,lwd=2,xlab = "Hour of Day",ylab = expression(delta^13*"C"[DIC]),xlim=c(0,24))
dev.off()

# plot model drivers ------------------------------------------------------

pdf("DE_plots_calcite_drivers.pdf") 
driver_plots <- par(mfrow=c(3,1)) 
driver_plots <- plot(t_hr, photo,type='l',lty=1,lwd=2,xlab = "Hour of Day",ylab = expression("F"[photo]*" ("*mu*"mol/kg/hr"),xlim=c(0,24))
abline(h=0, lty=2,lwd=1, col="grey")
driver_plots <- plot(t_hr, Fcarb_model,type='l',lty=1,lwd=2,xlab = "Hour of Day",ylab = expression("F"[carb]*" ("*mu*"mol/kg/hr"),xlim=c(0,24))
abline(h=0, lty=2,lwd=1, col="grey")
driver_plots <- plot(t_hr, Fgas_model,type='l',lty=1,lwd=2,xlab = "Hour of Day",ylab = expression("F"[gas]*" ("*mu*"mol/kg/hr"),xlim=c(0,24))
abline(h=0, lty=2,lwd=1, col="grey")
dev.off()

# Burton and Walter (1987) ------------------------------------------------------

R_BR <- zeros(1,length(Omega_cal_model))
for (n in 1:length(R_BR)){
  if (Omega_cal_model[n] >= 1){
    R_BR[n] <- k_BR*(Omega_cal_model[n] - 1)^n_BR #{umol/m^2/hr}
  } else if (Omega_cal_model[n] < 1){
    R_BR[n] <- -k_diss*(1 - Omega_cal_model[n])^n_diss #{umol/m^2/hr}
  }
}

# calculate cumulative precipitation for both rate laws ------------------------------------------------------

CR_BR <- zeros(1,length(R_BR))

if (R_BR[1]>=0){
  CR_BR[1]<-R_BR[1]
} else if (R_BR[1]<0){
  CR_BR[1]<-0
}

for (n in 2:length(R_BR)){
  if (R_BR[n]<0){
    CR_BR[n] <- CR_BR[n-1]
  } else if (R_BR[n]>=0){
  delta_t_precip <- t_hr[n]-t_hr[n-1]
  CR_BR[n] <- CR_BR[n-1] + R_BR[n]*delta_t_precip
  }
}

# calculate instantaneous d13C_carb and cumulative d13C_carb ------------------------------------------------------

d13C_carb_inst <- d13C_DIC_model + eps_DIC_cal

d13C_cumulative_BR <- zeros(1,length(CR_BR))
d13C_cumulative_BR[1] <- d13C_carb_inst[1]

for (m in 2:length(CR_BR)){
  delta_t_d13C <- t_hr[m]-t_hr[m-1]
  if (R_BR[m] >= 0){
    d13C_cumulative_BR[m] <- (d13C_cumulative_BR[m-1]*CR_BR[m-1] + d13C_carb_inst[m]*R_BR[m]*delta_t_d13C)/CR_BR[m]
  } else if (R_BR[m] < 0){
    d13C_cumulative_BR[m] <- d13C_cumulative_BR[m-1]
  }
}

# final plots ------------------------------------------------------

pdf("DE_plots_calcite_isotopes.pdf") 
final_plots <- par(mfrow=c(3,1)) 
final_plots <- plot(t_hr, R_BR,type='l',lty=1,lwd=2,xlab = "",ylab = "",xlim=c(0,24),axes=FALSE,col='black')
axis(2, ylim=c(150,250),col="black",las=1)  
mtext(expression("R_a_r ("*mu*"mol/m"^2*"/hr)"),side=2,line=2.5)
box()
par(new=TRUE)
plot(t_hr, CR_BR,type='l',lty=1,lwd=2,xlab = "",ylab = "",xlim=c(0,24),axes=FALSE,col="red")
axis(4, ylim=c(0,5000), col="red",col.axis="red",las=1)
mtext(expression("Cumulative precip. ("*mu*"mol/m"^2*")"),side=4,col="red",line=4) 
axis(1,pretty(range(t_hr),5))
mtext("Hour of Day",side=1,col="black",line=2.5) 

final_plots <- plot(t_hr, d13C_DIC_model,type='l',lty=1,lwd=2,xlab = "Hour of Day",ylab = expression(delta^13*"C"[DIC]),xlim=c(0,24),col='black')
final_plots <- plot(t_hr, d13C_cumulative_BR,type='l',lty=1,lwd=2,xlab = "Hour of Day",ylab = expression("Cumulative "*delta^13*"C"[carb]),xlim=c(0,24),col='black')
dev.off()



