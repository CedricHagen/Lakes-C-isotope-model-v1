#This script contains a wrapper for the diurnal engine code
#here modified as a function. This allows for easier 
#parameter exploration and model simulations. This version
#of the code is for calcite precipitation.

#Authors: Drs. Lizzy Trower and Cedric Hagen, 
#         University of Colorado, Boulder  

#Last Modified on Jan 14 2025

diurnalengine <- function(params,photo_curve,plot_flag) {
  
  tempC <- params[1]
  
  tempK <- tempC + 273 #{K}
  
  K0 <- exp(-58.0931 + 90.5697*(100/tempK) + 22.294*log(tempK/100) + 35*(0.027766 - 0.025888*(tempK/100) + 0.005078*(tempK/100)^2)) #{mol/m^3/atm}
  K0 = K0 *10^3
  
  Sc = 2116.8 - 136.25*tempC + 4.7353*(tempC^2) - 0.092307*(tempC^3) + 0.0007555*(tempC^4) #{dimensionless}
  
  u <- 4 #{m/s}
  
  pCO2_atm <- params[2] #{uatm}
  
  kCO2 <- 0.251*(u^2)*(Sc/660)^0.5 #{cm/hr}
  kCO2 <- kCO2/100 #{m/hr}
  
  waterdepth <- params[3] #{m}
  waterdensity <- params[4] #{kg/m^3}
  
  t_hr <- seq(from = 0, to = 24, by = 0.1)
  
  #kappa_p <- 150 #{umol/kg}
  #period <- 24 #{hr}
  #offset <- 5 #{hr} when the sine curve will cross 0
  #kappa_p_factor <- kappa_p/(period/pi) #{umol/kg/hr}
  #photo <- kappa_p_factor*sin((2*pi/period)*(t_hr-offset)) #{umol/kg/hr}
  photo <- photo_curve
  
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
  
  #k_rate <- 9*(10^-9) #{mol/m^2/s}
  #k_rate <- k_rate*60*60 #{mol/m^2/hr}
  
  precipkin <- calciteinterp(tempC)
  
  k_BR <- unlist(precipkin[1])
  n_BR <- unlist(precipkin[2])
  
  n_diss <- 2.86 # [Mean Low-Mg Calcite n from WM'85]
  k_diss <- 2239 #{umol/m^2/hr) [Mean Low-Mg Calcite k from WM'85, T8]
  
  eps_DIC_cal <- params[5]
  eps_g_DIC <- params[6]
  eps_DIC_g <- params[7]
  eps_DIC_org <- params[8]
  
  cCa <- params[9] #{mmol/kg}
  cMg <- params[10] #{mmol/kg}
  cNa <- params[11] #{mmol/kg}
  cK <- params[12] #{mmol/kg}
  cCl <- params[13] #{mmol/kg}
  cSO4 <- params[14] #{mmol/kg}
  
  DIC_model <- zeros(1,length(t_hr))
  #DIC_model[1] <- params[16] #{umol/kg}
  Alk_model <- zeros(1,length(t_hr))
  Alk_model[1] <- params[15] #{umol/kg}
  
  pH_model <- zeros(1,length(t_hr))
  #pH_model[1] <- params[15] #{umol/kg}
  pCO2_model <- zeros(1,length(t_hr))
  pCO2_model[1] <- params[16] #{umol/kg}
  Omega_cal_model <- zeros(1,length(t_hr))
  Fcarb_model <- zeros(1,length(t_hr))
  Fgas_model <- zeros(1,length(t_hr))
  d13C_DIC_model <- zeros(1,length(t_hr))
  d13C_DIC_model[1] <- params[17]
  d13C_org_model <- params[18]*ones(1,length(t_hr)) 
  
  library(phreeqc)                       # Load the PHREEQC library
  phrLoadDatabaseString(phreeqc.dat)     # Use the phreeqc database
  phrSetOutputStringsOn(TRUE)            # Format the output as a character string 
  
  co2_val <- log10(pCO2_model[1]/10^6)
  
  input <- c(                            # Text string defining composition, units, and temperature.
    '  SOLUTION              '                                            ,   
    '  units         mmol/kgw'                                            ,
    paste('  temp              ',as.character(tempC)),
    paste('  Alkalinity        ',as.character(Alk_model[1]/1000)),
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
    '  -si               aragonite',
    '  -ionic_strength   TRUE',
    '  -activities       Ca+2 Mg+2',
    '  -si               calcite',
    ' -totals            C(4)')
  
  phrRunString(input)                    # Run the input string
  #output <- phrGetOutputStrings()        # Save the results in 'output'
  output <- phrGetSelectedOutput()
  
  pH_model[1] <- output$n1$pH
  
  Omega_cal_model[1] <- 10^output$n1$si_calcite
  #pCO2_model[1] <- (10^output$n1$si_CO2.g.)*(10^6)
  
  DIC_model[1] <- output$n1$C.4..mol.kgw.*10^6 #{umol/kg}
  #Alk_model[1] <- output$n1$Alk.eq.kgw.*10^6 #{ueq/kg}
  
  Fcarb_model[1] <- k_BR*(Omega_cal_model[1] - 1)^n_BR/waterdensity/waterdepth #{umol/kg/hr}
  Fgas_model[1] <- kCO2*K0*(pCO2_model[1] - pCO2_atm)/waterdensity/waterdepth #{umol/kg/hr}
  
  #set mean d13C_org that's used for remineralized OM
  d13C_org_mean <- params[19] #{permil}
  
  #set mean d13C_carb that's used for dissolving carbonate
  d13C_carb_mean <- params[20] #{permil}
  
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
  
  if (plot_flag==1){
    pdf("DE_model.pdf") 
    out_plots <- par(mfrow=c(3,2)) 
    out_plots <- plot(t_hr, DIC_model/1000,type='l',lty=2,lwd=2,xlab = "Hour of Day",ylab = "[DIC] (mmol/kg)",xlim=c(0,24))
    out_plots <- plot(t_hr, Alk_model/1000,type='l',lty=2,lwd=2,xlab = "Hour of Day",ylab = "Alk (mmol/kg)",xlim=c(0,24))
    out_plots <- plot(t_hr, pH_model,type='l',lty=2,lwd=2,xlab = "Hour of Day",ylab = "pH",xlim=c(0,24))
    out_plots <- plot(t_hr, pCO2_model,type='l',lty=2,lwd=2,xlab = "Hour of Day",ylab = expression(paste('pCO_2 (',mu,'atm)',sep='')),xlim=c(0,24))
    out_plots <- plot(t_hr, Omega_cal_model,type='l',lty=2,lwd=2,xlab = "Hour of Day",ylab = expression(Omega[calc]),xlim=c(0,24))
    out_plots <- plot(t_hr, d13C_DIC_model,type='l',lty=2,lwd=2,xlab = "Hour of Day",ylab = expression(delta^13*"C"[DIC]),xlim=c(0,24))
    dev.off()
    
    pdf("DE_drivers.pdf") 
    driver_plots <- par(mfrow=c(3,1)) 
    driver_plots <- plot(t_hr, photo,type='l',lty=1,lwd=2,xlab = "Hour of Day",ylab = expression("F"[photo]*" ("*mu*"mol/kg/hr"),xlim=c(0,24))
    abline(h=0, lty=2,lwd=1, col="grey")
    driver_plots <- plot(t_hr, Fcarb_model,type='l',lty=1,lwd=2,xlab = "Hour of Day",ylab = expression("F"[carb]*" ("*mu*"mol/kg/hr"),xlim=c(0,24))
    abline(h=0, lty=2,lwd=1, col="grey")
    driver_plots <- plot(t_hr, Fgas_model,type='l',lty=1,lwd=2,xlab = "Hour of Day",ylab = expression("F"[gas]*" ("*mu*"mol/kg/hr"),xlim=c(0,24))
    abline(h=0, lty=2,lwd=1, col="grey")
    dev.off()
  }
  
  R_BR <- zeros(1,length(Omega_cal_model))
  for (n in 1:length(R_BR)){
    if (Omega_cal_model[n] >= 1){
      R_BR[n] <- k_BR*(Omega_cal_model[n] - 1)^n_BR #{umol/m^2/hr}
    } else if (Omega_cal_model[n] < 1){
      R_BR[n] <- -k_diss*(1 - Omega_cal_model[n])^n_diss #{umol/m^2/hr}
    }
  }
  
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
  
  if (plot_flag==1){
    pdf("DE_isotopes.pdf") 
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
  }
  
  all_out <- c(DIC_model/1000,Alk_model/1000,pH_model,pCO2_model,Omega_cal_model,d13C_DIC_model,photo,Fcarb_model,Fgas_model,R_BR,CR_BR,d13C_cumulative_BR);
  
  return(all_out);
  
}
# set input parameters ------------------------------------------------------

temp1 <- 20 #deg C
temp2 <- 25 #deg C
temp3 <- 30 #deg C
pCO2_atmo <- 400 #ppm
water_depth <- 1 #m
water_density_gl <- 1000 #kg/m3
water_density_pl <- 1000 #kg/m3
water_density <- 1080 #kg/m3
water_density_cl <- 1047 #kg/m3
eps_DIC_cal <- 1.0 #permil
eps_g_DIC <- -2 #permil
eps_DIC_g <- -10.3 #permil
eps_DIC_org <- -27.23 #permil

Ca <- 4.5 #mmol/kg
Ca_gl <- 9.8 #mmol/kg
Ca_pl <- 1.1 #mmol/kg
Ca_cl <- 24.5 #mmol/l

Mg <- 156 #mmol/kg
Mg_gl <- 2.5 #mmol/kg
Mg_pl <- 0.8 #mmol/kg
Mg_cl <- 102 #mmol/l

Na <- 1433 #mmol/kg
Na_gl <- 1.2 #mmol/kg
Na_pl <- 0.3 #mmol/kg
Na_cl <- 864 #mmol/l

K <- 58 #mmol/kg
K_gl <- 0.1 #mmol/kg
K_pl <- 0.1 #mmol/kg
K_cl <- 13.4 #mmol/l

Cl <- 1659 #mmol/kg
Cl_gl <- 0.6 #mmol/kg
Cl_pl <- 0.1 #mmol/kg
Cl_cl <- 952 #mmol/l

SO4 <- 78 #mmol/kg
SO4_gl <- 9.8 #mmol/kg
SO4_pl <- 0.1 #mmol/kg
SO4_cl <- 42.5 #mmol/l

#DIC <- 3700
#DIC_gl <- 2470
#DIC_pl <- 3097

pCO2_val <- 400 #ppm
pCO2_val_gl <- 400 #ppm
pCO2_val_pl <- 400 #ppm
pCO2_val_cl <- 400 #ppm

Alk <- 4600
Alk_gl <- 2420
Alk_pl <- 1510
#Alk_pl <- 750 #stand-in as 'synthetic lake' test
Alk_cl <- 3322

d13C_DIC <- -0.7 #permil
d13C_DIC_gl <- -7 #permil
d13C_DIC_pl <- 0 #permil
d13C_DIC_cl <- -1.55 #permil

d13C_org <- -7.3 #permil
d13C_org_gl <- -14.3 #permil
d13C_org_pl <- -7.3 #permil
d13C_org_cl <- -8.7 #permil

d13C_org_mean <- -20 #permil
d13C_org_mean_gl <- -27.7 #permil
d13C_org_mean_pl <- -23.85 #permil
d13C_org_mean_cl <- -20 #permil

d13C_carb_mean <- 4 #permil
d13C_carb_mean_gl <- -2.4 #permil
d13C_carb_mean_pl <- 0.8 #permil
d13C_carb_mean_cl <- 5.2 #permil

params_gsl_temp1 <- c(temp1,pCO2_atmo,water_depth,water_density,eps_DIC_cal,eps_g_DIC,eps_DIC_g,eps_DIC_org,Ca,Mg,Na,K,Cl,SO4,Alk,pCO2_val,d13C_DIC,d13C_org,d13C_org_mean,d13C_carb_mean);
params_gsl_temp2 <- c(temp2,pCO2_atmo,water_depth,water_density,eps_DIC_cal,eps_g_DIC,eps_DIC_g,eps_DIC_org,Ca,Mg,Na,K,Cl,SO4,Alk,pCO2_val,d13C_DIC,d13C_org,d13C_org_mean,d13C_carb_mean);
params_gsl_temp3 <- c(temp3,pCO2_atmo,water_depth,water_density,eps_DIC_cal,eps_g_DIC,eps_DIC_g,eps_DIC_org,Ca,Mg,Na,K,Cl,SO4,Alk,pCO2_val,d13C_DIC,d13C_org,d13C_org_mean,d13C_carb_mean);

params_gl_temp1 <- c(temp1,pCO2_atmo,water_depth,water_density_gl,eps_DIC_cal,eps_g_DIC,eps_DIC_g,eps_DIC_org,Ca_gl,Mg_gl,Na_gl,K_gl,Cl_gl,SO4_gl,Alk_gl,pCO2_val_gl,d13C_DIC_gl,d13C_org_gl,d13C_org_mean_gl,d13C_carb_mean_gl);
params_gl_temp2 <- c(temp2,pCO2_atmo,water_depth,water_density_gl,eps_DIC_cal,eps_g_DIC,eps_DIC_g,eps_DIC_org,Ca_gl,Mg_gl,Na_gl,K_gl,Cl_gl,SO4_gl,Alk_gl,pCO2_val_gl,d13C_DIC_gl,d13C_org_gl,d13C_org_mean_gl,d13C_carb_mean_gl);
params_gl_temp3 <- c(temp3,pCO2_atmo,water_depth,water_density_gl,eps_DIC_cal,eps_g_DIC,eps_DIC_g,eps_DIC_org,Ca_gl,Mg_gl,Na_gl,K_gl,Cl_gl,SO4_gl,Alk_gl,pCO2_val_gl,d13C_DIC_gl,d13C_org_gl,d13C_org_mean_gl,d13C_carb_mean_gl);

params_pl_temp1 <- c(temp1,pCO2_atmo,water_depth,water_density_pl,eps_DIC_cal,eps_g_DIC,eps_DIC_g,eps_DIC_org,Ca_pl,Mg_pl,Na_pl,K_pl,Cl_pl,SO4_pl,Alk_pl,pCO2_val_pl,d13C_DIC_pl,d13C_org_pl,d13C_org_mean_pl,d13C_carb_mean_pl);
params_pl_temp2 <- c(temp2,pCO2_atmo,water_depth,water_density_pl,eps_DIC_cal,eps_g_DIC,eps_DIC_g,eps_DIC_org,Ca_pl,Mg_pl,Na_pl,K_pl,Cl_pl,SO4_pl,Alk_pl,pCO2_val_pl,d13C_DIC_pl,d13C_org_pl,d13C_org_mean_pl,d13C_carb_mean_pl);
params_pl_temp3 <- c(temp3,pCO2_atmo,water_depth,water_density_pl,eps_DIC_cal,eps_g_DIC,eps_DIC_g,eps_DIC_org,Ca_pl,Mg_pl,Na_pl,K_pl,Cl_pl,SO4_pl,Alk_pl,pCO2_val_pl,d13C_DIC_pl,d13C_org_pl,d13C_org_mean_pl,d13C_carb_mean_pl);

params_cl_temp1 <- c(temp1,pCO2_atmo,water_depth,water_density_cl,eps_DIC_cal,eps_g_DIC,eps_DIC_g,eps_DIC_org,Ca_cl,Mg_cl,Na_cl,K_cl,Cl_cl,SO4_cl,Alk_cl,pCO2_val_cl,d13C_DIC_cl,d13C_org_cl,d13C_org_mean_cl,d13C_carb_mean_cl);
params_cl_temp2 <- c(temp2,pCO2_atmo,water_depth,water_density_cl,eps_DIC_cal,eps_g_DIC,eps_DIC_g,eps_DIC_org,Ca_cl,Mg_cl,Na_cl,K_cl,Cl_cl,SO4_cl,Alk_cl,pCO2_val_cl,d13C_DIC_cl,d13C_org_cl,d13C_org_mean_cl,d13C_carb_mean_cl);
params_cl_temp3 <- c(temp3,pCO2_atmo,water_depth,water_density_cl,eps_DIC_cal,eps_g_DIC,eps_DIC_g,eps_DIC_org,Ca_cl,Mg_cl,Na_cl,K_cl,Cl_cl,SO4_cl,Alk_cl,pCO2_val_cl,d13C_DIC_cl,d13C_org_cl,d13C_org_mean_cl,d13C_carb_mean_cl);

#define photo curve
t_hr <- seq(from = 0, to = 24, by = 0.1)
kappa_p <- 150 #{umol/kg}
period <- 24 #{hr}
offset <- 5 #{hr}
kappa_p_factor <- kappa_p/(period/pi) #{umol/kg/hr}
photo_curve <- kappa_p_factor*sin((2*pi/period)*(t_hr-offset)) #{umol/kg/hr}

# run single iteration ------------------------------------------------------

all_out <- diurnalengine(params_pl_temp2,photo_curve,1) #Pavilion Lake at T=25C
DIC_model <- all_out[1:length(photo_curve)]
Alk_model <- all_out[(1+length(photo_curve)*1):(length(photo_curve)+length(photo_curve)*1)]
pH_model <- all_out[(1+length(photo_curve)*2):(length(photo_curve)+length(photo_curve)*2)]
pCO2_model <- all_out[(1+length(photo_curve)*3):(length(photo_curve)+length(photo_curve)*3)]
Omega_cal_model <- all_out[(1+length(photo_curve)*4):(length(photo_curve)+length(photo_curve)*4)]
d13C_DIC_model <- all_out[(1+length(photo_curve)*5):(length(photo_curve)+length(photo_curve)*5)]
photo <- all_out[(1+length(photo_curve)*6):(length(photo_curve)+length(photo_curve)*6)]
Fcarb_model <- all_out[(1+length(photo_curve)*7):(length(photo_curve)+length(photo_curve)*7)]
Fgas_model <- all_out[(1+length(photo_curve)*8):(length(photo_curve)+length(photo_curve)*8)]
R_BR <- all_out[(1+length(photo_curve)*9):(length(photo_curve)+length(photo_curve)*9)]
CR_BR <- all_out[(1+length(photo_curve)*10):(length(photo_curve)+length(photo_curve)*10)]
d13C_cumulative_BR<- all_out[(1+length(photo_curve)*11):(length(photo_curve)+length(photo_curve)*11)]


pdf("DE_d13Ccarb.pdf") 
#out_plots <- par(mfrow=c(2,2)) 
out_plots <- plot(t_hr, DIC_model,type='l',lwd=2,xlab = "Hour of Day",ylab = "[DIC] (umol/kg)",xlim=c(0,24))
out_plots <- plot(t_hr, Omega_cal_model,type='l',lwd=2,xlab = "Hour of Day",ylab = expression(Omega[calc]),xlim=c(0,24))
out_plots <- plot(t_hr, d13C_DIC_model,type='l',lwd=2,xlab = "Hour of Day",ylab = expression(delta^13*"C"[DIC]),xlim=c(0,24))
out_plots <- plot(t_hr, d13C_cumulative_BR,type='l',lwd=2,xlab = "Hour of Day",ylab = expression(delta^13*"C"[carb]),xlim=c(0,24))
abline(h = (mean(d13C_DIC_model)+1), col = "red", lwd = 2, lty = 2)
dev.off()

pdf("DE_plots_model.pdf") 
out_plots <- par(mfrow=c(3,2)) 
out_plots <- plot(t_hr, DIC_model/1000,type='l',lty=2,lwd=2,xlab = "Hour of Day",ylab = "[DIC] (mmol/kg)",xlim=c(0,24))
out_plots <- plot(t_hr, Alk_model/1000,type='l',lty=2,lwd=2,xlab = "Hour of Day",ylab = "Alk (mmol/kg)",xlim=c(0,24),ylim=c(0.0032,0.0034))
out_plots <- plot(t_hr, pH_model,type='l',lty=2,lwd=2,xlab = "Hour of Day",ylab = "pH",xlim=c(0,24))
out_plots <- plot(t_hr, pCO2_model,type='l',lty=2,lwd=2,xlab = "Hour of Day",ylab = expression(paste('pCO'[2]*' (',mu,'atm)',sep='')),xlim=c(0,24))
out_plots <- plot(t_hr, Omega_cal_model,type='l',lty=2,lwd=2,xlab = "Hour of Day",ylab = expression(Omega[calc]),xlim=c(0,24))
out_plots <- plot(t_hr, d13C_DIC_model,type='l',lty=2,lwd=2,xlab = "Hour of Day",ylab = expression(delta^13*"C"[DIC]),xlim=c(0,24))
dev.off()

pdf("DE_plots_drivers.pdf") 
driver_plots <- par(mfrow=c(3,1)) 
driver_plots <- plot(t_hr, photo,type='l',lty=1,lwd=2,xlab = "Hour of Day",ylab = expression("F"[photo]*" ("*mu*"mol/kg/hr"),xlim=c(0,24))
abline(h=0, lty=2,lwd=1, col="grey")
driver_plots <- plot(t_hr, Fcarb_model,type='l',lty=1,lwd=2,xlab = "Hour of Day",ylab = expression("F"[carb]*" ("*mu*"mol/kg/hr"),xlim=c(0,24))
abline(h=0, lty=2,lwd=1, col="grey")
driver_plots <- plot(t_hr, Fgas_model,type='l',lty=1,lwd=2,xlab = "Hour of Day",ylab = expression("F"[gas]*" ("*mu*"mol/kg/hr"),xlim=c(0,24))
abline(h=0, lty=2,lwd=1, col="grey")
dev.off()

# iterative run ------------------------------------------------------

pCO2_vals <- seq(from=200,to=1200,by=100)
#pCO2_vals <- seq(from=400,to=400,by=100)
kappa_vals <- seq(from=0, to=800, by=25)
#alk_vals <- seq(from=750, to=3000, by=50)
#pCO2_vals <- alk_vals

all_pco2_gsl1<-numeric()
all_kappa_gsl1<-numeric()
all_offset_gsl1<-numeric()

all_pco2_gsl2<-numeric()
all_kappa_gsl2<-numeric()
all_offset_gsl2<-numeric()

all_pco2_gsl3<-numeric()
all_kappa_gsl3<-numeric()
all_offset_gsl3<-numeric()

all_pco2_gl1<-numeric()
all_kappa_gl1<-numeric()
all_offset_gl1<-numeric()

all_pco2_gl2<-numeric()
all_kappa_gl2<-numeric()
all_offset_gl2<-numeric()

all_pco2_gl3<-numeric()
all_kappa_gl3<-numeric()
all_offset_gl3<-numeric()

all_pco2_pl1<-numeric()
all_kappa_pl1<-numeric()
all_offset_pl1<-numeric()

all_pco2_pl2<-numeric()
all_kappa_pl2<-numeric()
all_offset_pl2<-numeric()

all_pco2_pl3<-numeric()
all_kappa_pl3<-numeric()
all_offset_pl3<-numeric()

all_pco2_cl1<-numeric()
all_kappa_cl1<-numeric()
all_offset_cl1<-numeric()

all_pco2_cl2<-numeric()
all_kappa_cl2<-numeric()
all_offset_cl2<-numeric()

all_pco2_cl3<-numeric()
all_kappa_cl3<-numeric()
all_offset_cl3<-numeric()

all_pco2_syn<-numeric()
all_kappa_syn<-numeric()
all_offset_syn<-numeric()

all_alk_pl2<-numeric()

t_hr <- seq(from = 0, to = 24, by = 0.1)
period <- 24 #{hr}
offset <- 5 #{hr}

library(svMisc)

which_flag <- 9 #Dictates which lake chemistry-temperature scenario is run

  for (j in 1:length(pCO2_vals)){
    for (k in 1:length(kappa_vals)){
      
      extra <- nchar('||100%')
      width <- options()$width
      
      if (k==length(kappa_vals)){
        progress(j*k, (length(pCO2_vals)*length(kappa_vals)))
        Sys.sleep(0.02)
        if (j*k == (length(pCO2_vals)*length(kappa_vals))) message("Done!")
      }
      
      
      kappa_p_factor <- kappa_vals[k]/(period/pi) #{umol/kg/hr}
      photo_curve <- kappa_p_factor*sin((2*pi/period)*(t_hr-offset)) #{umol/kg/hr}
      
      if (which_flag==1){
      params_gsl_temp1[2] <- pCO2_vals[j];
      params_gsl_temp1[16] <- pCO2_vals[j];
      
      all_out_gsl1 <- diurnalengine(params_gsl_temp1,photo_curve,0)
      DIC_model_gsl1 <- all_out_gsl1[1:length(photo_curve)]
      Alk_model_gsl1 <- all_out_gsl1[(1+length(photo_curve)*1):(length(photo_curve)+length(photo_curve)*1)]
      pH_model_gsl1 <- all_out_gsl1[(1+length(photo_curve)*2):(length(photo_curve)+length(photo_curve)*2)]
      pCO2_model_gsl1 <- all_out_gsl1[(1+length(photo_curve)*3):(length(photo_curve)+length(photo_curve)*3)]
      Omega_cal_model_gsl1 <- all_out_gsl1[(1+length(photo_curve)*4):(length(photo_curve)+length(photo_curve)*4)]
      d13C_DIC_model_gsl1 <- all_out_gsl1[(1+length(photo_curve)*5):(length(photo_curve)+length(photo_curve)*5)]
      photo_gsl1 <- all_out_gsl1[(1+length(photo_curve)*6):(length(photo_curve)+length(photo_curve)*6)]
      Fcarb_model_gsl1 <- all_out_gsl1[(1+length(photo_curve)*7):(length(photo_curve)+length(photo_curve)*7)]
      Fgas_model_gsl1 <- all_out_gsl1[(1+length(photo_curve)*8):(length(photo_curve)+length(photo_curve)*8)]
      R_BR_gsl1 <- all_out_gsl1[(1+length(photo_curve)*9):(length(photo_curve)+length(photo_curve)*9)]
      CR_BR_gsl1 <- all_out_gsl1[(1+length(photo_curve)*10):(length(photo_curve)+length(photo_curve)*10)]
      d13C_cumulative_BR_gsl1<- all_out_gsl1[(1+length(photo_curve)*11):(length(photo_curve)+length(photo_curve)*11)]
      
      all_offset_gsl1<-c(all_offset_gsl1,(d13C_cumulative_BR_gsl1[241]-(mean(d13C_DIC_model_gsl1)+1)))
      all_pco2_gsl1<-c(all_pco2_gsl1,pCO2_vals[j])
      all_kappa_gsl1<-c(all_kappa_gsl1,kappa_vals[k])
      }
      
      if (which_flag==2){
      params_gsl_temp2[2] <- pCO2_vals[j];
      params_gsl_temp2[16] <- pCO2_vals[j];
      
      all_out_gsl2 <- diurnalengine(params_gsl_temp2,photo_curve,0)
      DIC_model_gsl2 <- all_out_gsl2[1:length(photo_curve)]
      Alk_model_gsl2 <- all_out_gsl2[(1+length(photo_curve)*1):(length(photo_curve)+length(photo_curve)*1)]
      pH_model_gsl2 <- all_out_gsl2[(1+length(photo_curve)*2):(length(photo_curve)+length(photo_curve)*2)]
      pCO2_model_gsl2 <- all_out_gsl2[(1+length(photo_curve)*3):(length(photo_curve)+length(photo_curve)*3)]
      Omega_cal_model_gsl2 <- all_out_gsl2[(1+length(photo_curve)*4):(length(photo_curve)+length(photo_curve)*4)]
      d13C_DIC_model_gsl2 <- all_out_gsl2[(1+length(photo_curve)*5):(length(photo_curve)+length(photo_curve)*5)]
      photo_gsl2 <- all_out_gsl2[(1+length(photo_curve)*6):(length(photo_curve)+length(photo_curve)*6)]
      Fcarb_model_gsl2 <- all_out_gsl2[(1+length(photo_curve)*7):(length(photo_curve)+length(photo_curve)*7)]
      Fgas_model_gsl2 <- all_out_gsl2[(1+length(photo_curve)*8):(length(photo_curve)+length(photo_curve)*8)]
      R_BR_gsl2 <- all_out_gsl2[(1+length(photo_curve)*9):(length(photo_curve)+length(photo_curve)*9)]
      CR_BR_gsl2 <- all_out_gsl2[(1+length(photo_curve)*10):(length(photo_curve)+length(photo_curve)*10)]
      d13C_cumulative_BR_gsl2<- all_out_gsl2[(1+length(photo_curve)*11):(length(photo_curve)+length(photo_curve)*11)]
      
      all_offset_gsl2<-c(all_offset_gsl2,(d13C_cumulative_BR_gsl2[241]-(mean(d13C_DIC_model_gsl2)+1)))
      all_pco2_gsl2<-c(all_pco2_gsl2,pCO2_vals[j])
      all_kappa_gsl2<-c(all_kappa_gsl2,kappa_vals[k])
      }
      
      if (which_flag==3){
      params_gsl_temp3[2] <- pCO2_vals[j];
      params_gsl_temp3[16] <- pCO2_vals[j];
        
      all_out_gsl3 <- diurnalengine(params_gsl_temp3,photo_curve,0)
      DIC_model_gsl3 <- all_out_gsl3[1:length(photo_curve)]
      Alk_model_gsl3 <- all_out_gsl3[(1+length(photo_curve)*1):(length(photo_curve)+length(photo_curve)*1)]
      pH_model_gsl3 <- all_out_gsl3[(1+length(photo_curve)*2):(length(photo_curve)+length(photo_curve)*2)]
      pCO2_model_gsl3 <- all_out_gsl3[(1+length(photo_curve)*3):(length(photo_curve)+length(photo_curve)*3)]
      Omega_cal_model_gsl3 <- all_out_gsl3[(1+length(photo_curve)*4):(length(photo_curve)+length(photo_curve)*4)]
      d13C_DIC_model_gsl3 <- all_out_gsl3[(1+length(photo_curve)*5):(length(photo_curve)+length(photo_curve)*5)]
      photo_gsl3 <- all_out_gsl3[(1+length(photo_curve)*6):(length(photo_curve)+length(photo_curve)*6)]
      Fcarb_model_gsl3 <- all_out_gsl3[(1+length(photo_curve)*7):(length(photo_curve)+length(photo_curve)*7)]
      Fgas_model_gsl3 <- all_out_gsl3[(1+length(photo_curve)*8):(length(photo_curve)+length(photo_curve)*8)]
      R_BR_gsl3 <- all_out_gsl3[(1+length(photo_curve)*9):(length(photo_curve)+length(photo_curve)*9)]
      CR_BR_gsl3 <- all_out_gsl3[(1+length(photo_curve)*10):(length(photo_curve)+length(photo_curve)*10)]
      d13C_cumulative_BR_gsl3<- all_out_gsl3[(1+length(photo_curve)*11):(length(photo_curve)+length(photo_curve)*11)]
      
      all_offset_gsl3<-c(all_offset_gsl3,(d13C_cumulative_BR_gsl3[241]-(mean(d13C_DIC_model_gsl3)+1)))
      all_pco2_gsl3<-c(all_pco2_gsl3,pCO2_vals[j])
      all_kappa_gsl3<-c(all_kappa_gsl3,kappa_vals[k])
      }
      
      if (which_flag==4){
      params_gl_temp1[2] <- pCO2_vals[j];
      params_gl_temp1[16] <- pCO2_vals[j];
      
      all_out_gl1 <- diurnalengine(params_gl_temp1,photo_curve,0)
      DIC_model_gl1 <- all_out_gl1[1:length(photo_curve)]
      Alk_model_gl1 <- all_out_gl1[(1+length(photo_curve)*1):(length(photo_curve)+length(photo_curve)*1)]
      pH_model_gl1 <- all_out_gl1[(1+length(photo_curve)*2):(length(photo_curve)+length(photo_curve)*2)]
      pCO2_model_gl1 <- all_out_gl1[(1+length(photo_curve)*3):(length(photo_curve)+length(photo_curve)*3)]
      Omega_cal_model_gl1 <- all_out_gl1[(1+length(photo_curve)*4):(length(photo_curve)+length(photo_curve)*4)]
      d13C_DIC_model_gl1 <- all_out_gl1[(1+length(photo_curve)*5):(length(photo_curve)+length(photo_curve)*5)]
      photo_gl1 <- all_out_gl1[(1+length(photo_curve)*6):(length(photo_curve)+length(photo_curve)*6)]
      Fcarb_model_gl1 <- all_out_gl1[(1+length(photo_curve)*7):(length(photo_curve)+length(photo_curve)*7)]
      Fgas_model_gl1 <- all_out_gl1[(1+length(photo_curve)*8):(length(photo_curve)+length(photo_curve)*8)]
      R_BR_gl1 <- all_out_gl1[(1+length(photo_curve)*9):(length(photo_curve)+length(photo_curve)*9)]
      CR_BR_gl1 <- all_out_gl1[(1+length(photo_curve)*10):(length(photo_curve)+length(photo_curve)*10)]
      d13C_cumulative_BR_gl1<- all_out_gl1[(1+length(photo_curve)*11):(length(photo_curve)+length(photo_curve)*11)]
      
      all_offset_gl1<-c(all_offset_gl1,(d13C_cumulative_BR_gl1[241]-(mean(d13C_DIC_model_gl1)+1)))
      all_pco2_gl1<-c(all_pco2_gl1,pCO2_vals[j])
      all_kappa_gl1<-c(all_kappa_gl1,kappa_vals[k])
      }
      
      if (which_flag==5){
      params_gl_temp2[2] <- pCO2_vals[j];
      params_gl_temp2[16] <- pCO2_vals[j];
      
      all_out_gl2 <- diurnalengine(params_gl_temp2,photo_curve,0)
      DIC_model_gl2 <- all_out_gl2[1:length(photo_curve)]
      Alk_model_gl2 <- all_out_gl2[(1+length(photo_curve)*1):(length(photo_curve)+length(photo_curve)*1)]
      pH_model_gl2 <- all_out_gl2[(1+length(photo_curve)*2):(length(photo_curve)+length(photo_curve)*2)]
      pCO2_model_gl2 <- all_out_gl2[(1+length(photo_curve)*3):(length(photo_curve)+length(photo_curve)*3)]
      Omega_cal_model_gl2 <- all_out_gl2[(1+length(photo_curve)*4):(length(photo_curve)+length(photo_curve)*4)]
      d13C_DIC_model_gl2 <- all_out_gl2[(1+length(photo_curve)*5):(length(photo_curve)+length(photo_curve)*5)]
      photo_gl2 <- all_out_gl2[(1+length(photo_curve)*6):(length(photo_curve)+length(photo_curve)*6)]
      Fcarb_model_gl2 <- all_out_gl2[(1+length(photo_curve)*7):(length(photo_curve)+length(photo_curve)*7)]
      Fgas_model_gl2 <- all_out_gl2[(1+length(photo_curve)*8):(length(photo_curve)+length(photo_curve)*8)]
      R_BR_gl2 <- all_out_gl2[(1+length(photo_curve)*9):(length(photo_curve)+length(photo_curve)*9)]
      CR_BR_gl2 <- all_out_gl2[(1+length(photo_curve)*10):(length(photo_curve)+length(photo_curve)*10)]
      d13C_cumulative_BR_gl2<- all_out_gl2[(1+length(photo_curve)*11):(length(photo_curve)+length(photo_curve)*11)]
      
      all_offset_gl2<-c(all_offset_gl2,(d13C_cumulative_BR_gl2[241]-(mean(d13C_DIC_model_gl2)+1)))
      all_pco2_gl2<-c(all_pco2_gl2,pCO2_vals[j])
      all_kappa_gl2<-c(all_kappa_gl2,kappa_vals[k])
      }
      
      if (which_flag==6){
      params_gl_temp3[2] <- pCO2_vals[j];
      params_gl_temp3[16] <- pCO2_vals[j];
      
      all_out_gl3 <- diurnalengine(params_gl_temp3,photo_curve,0)
      DIC_model_gl3 <- all_out_gl3[1:length(photo_curve)]
      Alk_model_gl3 <- all_out_gl3[(1+length(photo_curve)*1):(length(photo_curve)+length(photo_curve)*1)]
      pH_model_gl3 <- all_out_gl3[(1+length(photo_curve)*2):(length(photo_curve)+length(photo_curve)*2)]
      pCO2_model_gl3 <- all_out_gl3[(1+length(photo_curve)*3):(length(photo_curve)+length(photo_curve)*3)]
      Omega_cal_model_gl3 <- all_out_gl3[(1+length(photo_curve)*4):(length(photo_curve)+length(photo_curve)*4)]
      d13C_DIC_model_gl3 <- all_out_gl3[(1+length(photo_curve)*5):(length(photo_curve)+length(photo_curve)*5)]
      photo_gl3 <- all_out_gl3[(1+length(photo_curve)*6):(length(photo_curve)+length(photo_curve)*6)]
      Fcarb_model_gl3 <- all_out_gl3[(1+length(photo_curve)*7):(length(photo_curve)+length(photo_curve)*7)]
      Fgas_model_gl3 <- all_out_gl3[(1+length(photo_curve)*8):(length(photo_curve)+length(photo_curve)*8)]
      R_BR_gl3 <- all_out_gl3[(1+length(photo_curve)*9):(length(photo_curve)+length(photo_curve)*9)]
      CR_BR_gl3 <- all_out_gl3[(1+length(photo_curve)*10):(length(photo_curve)+length(photo_curve)*10)]
      d13C_cumulative_BR_gl3<- all_out_gl3[(1+length(photo_curve)*11):(length(photo_curve)+length(photo_curve)*11)]
      
      all_offset_gl3<-c(all_offset_gl3,(d13C_cumulative_BR_gl3[241]-(mean(d13C_DIC_model_gl3)+1)))
      all_pco2_gl3<-c(all_pco2_gl3,pCO2_vals[j])
      all_kappa_gl3<-c(all_kappa_gl3,kappa_vals[k])
      }
      
      if (which_flag==7){
      params_pl_temp1[2] <- pCO2_vals[j];
      params_pl_temp1[16] <- pCO2_vals[j];
      
      all_out_pl1 <- diurnalengine(params_pl_temp1,photo_curve,0)
      DIC_model_pl1 <- all_out_pl1[1:length(photo_curve)]
      Alk_model_pl1 <- all_out_pl1[(1+length(photo_curve)*1):(length(photo_curve)+length(photo_curve)*1)]
      pH_model_pl1 <- all_out_pl1[(1+length(photo_curve)*2):(length(photo_curve)+length(photo_curve)*2)]
      pCO2_model_pl1 <- all_out_pl1[(1+length(photo_curve)*3):(length(photo_curve)+length(photo_curve)*3)]
      Omega_cal_model_pl1 <- all_out_pl1[(1+length(photo_curve)*4):(length(photo_curve)+length(photo_curve)*4)]
      d13C_DIC_model_pl1 <- all_out_pl1[(1+length(photo_curve)*5):(length(photo_curve)+length(photo_curve)*5)]
      photo_pl1 <- all_out_pl1[(1+length(photo_curve)*6):(length(photo_curve)+length(photo_curve)*6)]
      Fcarb_model_pl1 <- all_out_pl1[(1+length(photo_curve)*7):(length(photo_curve)+length(photo_curve)*7)]
      Fgas_model_pl1 <- all_out_pl1[(1+length(photo_curve)*8):(length(photo_curve)+length(photo_curve)*8)]
      R_BR_pl1 <- all_out_pl1[(1+length(photo_curve)*9):(length(photo_curve)+length(photo_curve)*9)]
      CR_BR_pl1 <- all_out_pl1[(1+length(photo_curve)*10):(length(photo_curve)+length(photo_curve)*10)]
      d13C_cumulative_BR_pl1<- all_out_pl1[(1+length(photo_curve)*11):(length(photo_curve)+length(photo_curve)*11)]
      
      all_offset_pl1<-c(all_offset_pl1,(d13C_cumulative_BR_pl1[241]-(mean(d13C_DIC_model_pl1)+1)))
      all_pco2_pl1<-c(all_pco2_pl1,pCO2_vals[j])
      all_kappa_pl1<-c(all_kappa_pl1,kappa_vals[k])
      }
      
      if (which_flag==8){
      params_pl_temp2[2] <- pCO2_vals[j];
      params_pl_temp2[16] <- pCO2_vals[j];
      
      all_out_pl2 <- diurnalengine(params_pl_temp2,photo_curve,0)
      DIC_model_pl2 <- all_out_pl2[1:length(photo_curve)]
      Alk_model_pl2 <- all_out_pl2[(1+length(photo_curve)*1):(length(photo_curve)+length(photo_curve)*1)]
      pH_model_pl2 <- all_out_pl2[(1+length(photo_curve)*2):(length(photo_curve)+length(photo_curve)*2)]
      pCO2_model_pl2 <- all_out_pl2[(1+length(photo_curve)*3):(length(photo_curve)+length(photo_curve)*3)]
      Omega_cal_model_pl2 <- all_out_pl2[(1+length(photo_curve)*4):(length(photo_curve)+length(photo_curve)*4)]
      d13C_DIC_model_pl2 <- all_out_pl2[(1+length(photo_curve)*5):(length(photo_curve)+length(photo_curve)*5)]
      photo_pl2 <- all_out_pl2[(1+length(photo_curve)*6):(length(photo_curve)+length(photo_curve)*6)]
      Fcarb_model_pl2 <- all_out_pl2[(1+length(photo_curve)*7):(length(photo_curve)+length(photo_curve)*7)]
      Fgas_model_pl2 <- all_out_pl2[(1+length(photo_curve)*8):(length(photo_curve)+length(photo_curve)*8)]
      R_BR_pl2 <- all_out_pl2[(1+length(photo_curve)*9):(length(photo_curve)+length(photo_curve)*9)]
      CR_BR_pl2 <- all_out_pl2[(1+length(photo_curve)*10):(length(photo_curve)+length(photo_curve)*10)]
      d13C_cumulative_BR_pl2<- all_out_pl2[(1+length(photo_curve)*11):(length(photo_curve)+length(photo_curve)*11)]
      
      all_offset_pl2<-c(all_offset_pl2,(d13C_cumulative_BR_pl2[241]-(mean(d13C_DIC_model_pl2)+1)))
      all_pco2_pl2<-c(all_pco2_pl2,pCO2_vals[j])
      all_kappa_pl2<-c(all_kappa_pl2,kappa_vals[k])
      }
      
      if (which_flag==9){
      params_pl_temp3[2] <- pCO2_vals[j];
      params_pl_temp3[16] <- pCO2_vals[j];
      
      all_out_pl3 <- diurnalengine(params_pl_temp3,photo_curve,0)
      DIC_model_pl3 <- all_out_pl3[1:length(photo_curve)]
      Alk_model_pl3 <- all_out_pl3[(1+length(photo_curve)*1):(length(photo_curve)+length(photo_curve)*1)]
      pH_model_pl3 <- all_out_pl3[(1+length(photo_curve)*2):(length(photo_curve)+length(photo_curve)*2)]
      pCO2_model_pl3 <- all_out_pl3[(1+length(photo_curve)*3):(length(photo_curve)+length(photo_curve)*3)]
      Omega_cal_model_pl3 <- all_out_pl3[(1+length(photo_curve)*4):(length(photo_curve)+length(photo_curve)*4)]
      d13C_DIC_model_pl3 <- all_out_pl3[(1+length(photo_curve)*5):(length(photo_curve)+length(photo_curve)*5)]
      photo_pl3 <- all_out_pl3[(1+length(photo_curve)*6):(length(photo_curve)+length(photo_curve)*6)]
      Fcarb_model_pl3 <- all_out_pl3[(1+length(photo_curve)*7):(length(photo_curve)+length(photo_curve)*7)]
      Fgas_model_pl3 <- all_out_pl3[(1+length(photo_curve)*8):(length(photo_curve)+length(photo_curve)*8)]
      R_BR_pl3 <- all_out_pl3[(1+length(photo_curve)*9):(length(photo_curve)+length(photo_curve)*9)]
      CR_BR_pl3 <- all_out_pl3[(1+length(photo_curve)*10):(length(photo_curve)+length(photo_curve)*10)]
      d13C_cumulative_BR_pl3<- all_out_pl3[(1+length(photo_curve)*11):(length(photo_curve)+length(photo_curve)*11)]
      
      all_offset_pl3<-c(all_offset_pl3,(d13C_cumulative_BR_pl3[241]-(mean(d13C_DIC_model_pl3)+1)))
      all_pco2_pl3<-c(all_pco2_pl3,pCO2_vals[j])
      all_kappa_pl3<-c(all_kappa_pl3,kappa_vals[k])
      }
      
      if (which_flag==10){
        params_cl_temp1[2] <- pCO2_vals[j];
        params_cl_temp1[16] <- pCO2_vals[j];
        
        all_out_cl1 <- diurnalengine(params_cl_temp1,photo_curve,0)
        DIC_model_cl1 <- all_out_cl1[1:length(photo_curve)]
        Alk_model_cl1 <- all_out_cl1[(1+length(photo_curve)*1):(length(photo_curve)+length(photo_curve)*1)]
        pH_model_cl1 <- all_out_cl1[(1+length(photo_curve)*2):(length(photo_curve)+length(photo_curve)*2)]
        pCO2_model_cl1 <- all_out_cl1[(1+length(photo_curve)*3):(length(photo_curve)+length(photo_curve)*3)]
        Omega_cal_model_cl1 <- all_out_cl1[(1+length(photo_curve)*4):(length(photo_curve)+length(photo_curve)*4)]
        d13C_DIC_model_cl1 <- all_out_cl1[(1+length(photo_curve)*5):(length(photo_curve)+length(photo_curve)*5)]
        photo_cl1 <- all_out_cl1[(1+length(photo_curve)*6):(length(photo_curve)+length(photo_curve)*6)]
        Fcarb_model_cl1 <- all_out_cl1[(1+length(photo_curve)*7):(length(photo_curve)+length(photo_curve)*7)]
        Fgas_model_cl1 <- all_out_cl1[(1+length(photo_curve)*8):(length(photo_curve)+length(photo_curve)*8)]
        R_BR_cl1 <- all_out_cl1[(1+length(photo_curve)*9):(length(photo_curve)+length(photo_curve)*9)]
        CR_BR_cl1 <- all_out_cl1[(1+length(photo_curve)*10):(length(photo_curve)+length(photo_curve)*10)]
        d13C_cumulative_BR_cl1<- all_out_cl1[(1+length(photo_curve)*11):(length(photo_curve)+length(photo_curve)*11)]
        
        all_offset_cl1<-c(all_offset_cl1,(d13C_cumulative_BR_cl1[241]-(mean(d13C_DIC_model_cl1)+1)))
        all_pco2_cl1<-c(all_pco2_cl1,pCO2_vals[j])
        all_kappa_cl1<-c(all_kappa_cl1,kappa_vals[k])
      }
      
      if (which_flag==11){
        params_cl_temp2[2] <- pCO2_vals[j];
        params_cl_temp2[16] <- pCO2_vals[j];
        
        all_out_cl2 <- diurnalengine(params_cl_temp2,photo_curve,0)
        DIC_model_cl2 <- all_out_cl2[1:length(photo_curve)]
        Alk_model_cl2 <- all_out_cl2[(1+length(photo_curve)*1):(length(photo_curve)+length(photo_curve)*1)]
        pH_model_cl2 <- all_out_cl2[(1+length(photo_curve)*2):(length(photo_curve)+length(photo_curve)*2)]
        pCO2_model_cl2 <- all_out_cl2[(1+length(photo_curve)*3):(length(photo_curve)+length(photo_curve)*3)]
        Omega_cal_model_cl2 <- all_out_cl2[(1+length(photo_curve)*4):(length(photo_curve)+length(photo_curve)*4)]
        d13C_DIC_model_cl2 <- all_out_cl2[(1+length(photo_curve)*5):(length(photo_curve)+length(photo_curve)*5)]
        photo_cl2 <- all_out_cl2[(1+length(photo_curve)*6):(length(photo_curve)+length(photo_curve)*6)]
        Fcarb_model_cl2 <- all_out_cl2[(1+length(photo_curve)*7):(length(photo_curve)+length(photo_curve)*7)]
        Fgas_model_cl2 <- all_out_cl2[(1+length(photo_curve)*8):(length(photo_curve)+length(photo_curve)*8)]
        R_BR_cl2 <- all_out_cl2[(1+length(photo_curve)*9):(length(photo_curve)+length(photo_curve)*9)]
        CR_BR_cl2 <- all_out_cl2[(1+length(photo_curve)*10):(length(photo_curve)+length(photo_curve)*10)]
        d13C_cumulative_BR_cl2<- all_out_cl2[(1+length(photo_curve)*11):(length(photo_curve)+length(photo_curve)*11)]
        
        all_offset_cl2<-c(all_offset_cl2,(d13C_cumulative_BR_cl2[241]-(mean(d13C_DIC_model_cl2)+1)))
        all_pco2_cl2<-c(all_pco2_cl2,pCO2_vals[j])
        all_kappa_cl2<-c(all_kappa_cl2,kappa_vals[k])
      }
      
      if (which_flag==12){
        params_cl_temp3[2] <- pCO2_vals[j];
        params_cl_temp3[16] <- pCO2_vals[j];
        
        all_out_cl3 <- diurnalengine(params_cl_temp3,photo_curve,0)
        DIC_model_cl3 <- all_out_cl3[1:length(photo_curve)]
        Alk_model_cl3 <- all_out_cl3[(1+length(photo_curve)*1):(length(photo_curve)+length(photo_curve)*1)]
        pH_model_cl3 <- all_out_cl3[(1+length(photo_curve)*2):(length(photo_curve)+length(photo_curve)*2)]
        pCO2_model_cl3 <- all_out_cl3[(1+length(photo_curve)*3):(length(photo_curve)+length(photo_curve)*3)]
        Omega_cal_model_cl3 <- all_out_cl3[(1+length(photo_curve)*4):(length(photo_curve)+length(photo_curve)*4)]
        d13C_DIC_model_cl3 <- all_out_cl3[(1+length(photo_curve)*5):(length(photo_curve)+length(photo_curve)*5)]
        photo_cl3 <- all_out_cl3[(1+length(photo_curve)*6):(length(photo_curve)+length(photo_curve)*6)]
        Fcarb_model_cl3 <- all_out_cl3[(1+length(photo_curve)*7):(length(photo_curve)+length(photo_curve)*7)]
        Fgas_model_cl3 <- all_out_cl3[(1+length(photo_curve)*8):(length(photo_curve)+length(photo_curve)*8)]
        R_BR_cl3 <- all_out_cl3[(1+length(photo_curve)*9):(length(photo_curve)+length(photo_curve)*9)]
        CR_BR_cl3 <- all_out_cl3[(1+length(photo_curve)*10):(length(photo_curve)+length(photo_curve)*10)]
        d13C_cumulative_BR_cl3<- all_out_cl3[(1+length(photo_curve)*11):(length(photo_curve)+length(photo_curve)*11)]
        
        all_offset_cl3<-c(all_offset_cl3,(d13C_cumulative_BR_cl3[241]-(mean(d13C_DIC_model_cl3)+1)))
        all_pco2_cl3<-c(all_pco2_cl3,pCO2_vals[j])
        all_kappa_cl3<-c(all_kappa_cl3,kappa_vals[k])
      }
      
      if (which_flag==18){
        params_pl_temp2[2] <- pCO2_vals[j];
        params_pl_temp2[16] <- pCO2_vals[j];
        
        all_out_pl2 <- diurnalengine(params_pl_temp2,photo_curve,0)
        DIC_model_pl2 <- all_out_pl2[1:length(photo_curve)]
        Alk_model_pl2 <- all_out_pl2[(1+length(photo_curve)*1):(length(photo_curve)+length(photo_curve)*1)]
        pH_model_pl2 <- all_out_pl2[(1+length(photo_curve)*2):(length(photo_curve)+length(photo_curve)*2)]
        pCO2_model_pl2 <- all_out_pl2[(1+length(photo_curve)*3):(length(photo_curve)+length(photo_curve)*3)]
        Omega_cal_model_pl2 <- all_out_pl2[(1+length(photo_curve)*4):(length(photo_curve)+length(photo_curve)*4)]
        d13C_DIC_model_pl2 <- all_out_pl2[(1+length(photo_curve)*5):(length(photo_curve)+length(photo_curve)*5)]
        photo_pl2 <- all_out_pl2[(1+length(photo_curve)*6):(length(photo_curve)+length(photo_curve)*6)]
        Fcarb_model_pl2 <- all_out_pl2[(1+length(photo_curve)*7):(length(photo_curve)+length(photo_curve)*7)]
        Fgas_model_pl2 <- all_out_pl2[(1+length(photo_curve)*8):(length(photo_curve)+length(photo_curve)*8)]
        R_BR_pl2 <- all_out_pl2[(1+length(photo_curve)*9):(length(photo_curve)+length(photo_curve)*9)]
        CR_BR_pl2 <- all_out_pl2[(1+length(photo_curve)*10):(length(photo_curve)+length(photo_curve)*10)]
        d13C_cumulative_BR_pl2<- all_out_pl2[(1+length(photo_curve)*11):(length(photo_curve)+length(photo_curve)*11)]
        
        all_offset_pl2<-c(all_offset_pl2,(d13C_cumulative_BR_pl2[241]-(mean(d13C_DIC_model_pl2)+1)))
        all_pco2_pl2<-c(all_pco2_pl2,pCO2_vals[j])
        all_kappa_pl2<-c(all_kappa_pl2,kappa_vals[k])
      }
      
      if (which_flag==13){
        params_pl_temp2[2] <- pCO2_vals[j];
        params_pl_temp2[16] <- pCO2_vals[j];
        
        all_out_syn <- diurnalengine(params_pl_temp2,photo_curve,0)
        DIC_model_syn <- all_out_syn[1:length(photo_curve)]
        Alk_model_syn <- all_out_syn[(1+length(photo_curve)*1):(length(photo_curve)+length(photo_curve)*1)]
        pH_model_syn <- all_out_syn[(1+length(photo_curve)*2):(length(photo_curve)+length(photo_curve)*2)]
        pCO2_model_syn <- all_out_syn[(1+length(photo_curve)*3):(length(photo_curve)+length(photo_curve)*3)]
        Omega_cal_model_syn <- all_out_syn[(1+length(photo_curve)*4):(length(photo_curve)+length(photo_curve)*4)]
        d13C_DIC_model_syn <- all_out_syn[(1+length(photo_curve)*5):(length(photo_curve)+length(photo_curve)*5)]
        photo_syn <- all_out_syn[(1+length(photo_curve)*6):(length(photo_curve)+length(photo_curve)*6)]
        Fcarb_model_syn <- all_out_syn[(1+length(photo_curve)*7):(length(photo_curve)+length(photo_curve)*7)]
        Fgas_model_syn <- all_out_syn[(1+length(photo_curve)*8):(length(photo_curve)+length(photo_curve)*8)]
        R_BR_syn <- all_out_syn[(1+length(photo_curve)*9):(length(photo_curve)+length(photo_curve)*9)]
        CR_BR_syn <- all_out_syn[(1+length(photo_curve)*10):(length(photo_curve)+length(photo_curve)*10)]
        d13C_cumulative_BR_syn<- all_out_syn[(1+length(photo_curve)*11):(length(photo_curve)+length(photo_curve)*11)]
        
        all_offset_syn<-c(all_offset_syn,(d13C_cumulative_BR_syn[241]-(mean(d13C_DIC_model_syn)+1)))
        all_pco2_syn<-c(all_pco2_syn,pCO2_vals[j])
        all_kappa_syn<-c(all_kappa_syn,kappa_vals[k])
      }
    }
  }

library(dplyr)
library(plotly)

#Make contour plots

fig <- plot_ly(x = ~all_kappa_pl3, y = ~all_pco2_pl3, z = ~all_offset_pl3,
               width = 600, height = 500,  type = "contour") %>% 
  colorbar(title = "Isotopic offset (‰)",limits = c(0, 10),breaks = c(2,4,6,8,10)) %>%
  layout(title = 'Pavilion Lake,  T = 30C',
         xaxis = list(title = 'κ'), yaxis = list(title = 'pCO2 (ppm)'))

fig

fig <- fig %>% add_trace(
  x = ~all_kappa_pl2, 
  y = ~all_co2_pl2/1000, 
  z = ~all_offset_pl2, 
  type = 'contour',
  contours = list(
    start = 3.67,
    end = 3.67,
    width = 0.1,
    coloring = 'none'),
  line = list(width = 1,color = 'white'),
  showlegend = FALSE,
  showscale = FALSE)

fig <- fig %>% add_trace(
  x = ~all_kappa_pl2, 
  y = ~all_pco2_pl2/1000, 
  z = ~all_offset_pl2, 
  type = 'contour',
  contours = list(
    start = 1.18,
    end = 1.18,
    width = 0.1,
    coloring = 'none'),
  line = list(width = 1,color = 'white',dash = 'dash'),
  showlegend = FALSE,
  showscale = FALSE)

fig <- fig %>% add_trace(
  x = ~all_kappa_pl2, 
  y = ~all_pco2_pl2/1000, 
  z = ~all_offset_pl2, 
  type = 'contour',
  contours = list(
    start = 7.64,
    end = 7.64,
    width = 0.1,
    coloring = 'none'),
  line = list(width = 1,color = 'white',dash = 'dash'),
  showlegend = FALSE,
  showscale = FALSE)
fig


#Kappa P vs Offset plot
pdf("Kappa_Offset_AllLakes.pdf") 
plot(all_kappa_syn,all_offset_syn,type="l",col="black",xlab='K',ylab="Isotopic offset (‰)",lwd=2,main="T = 25C, pCO2 = 400 ppm",ylim=c(0,12))
lines(all_kappa_pl2,all_offset_pl2,col="orange",lwd=2)
lines(all_kappa_gl2,all_offset_gl2,col="green",lwd=2)
lines(all_kappa_gsl2,all_offset_gsl2,col="blue",lwd=2)
legend("topleft", legend = c("Synthetic","Pavilion Lake","Green Lake","Great Salt Lake"), col = c("black","orange","green","blue"), lwd = 2)
dev.off()
