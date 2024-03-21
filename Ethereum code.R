
rm(list=ls()) 
library(dplyr)
library(MSGARCH)
library(rugarch)
library(copula)
library(GAS)
library(RobGARCHBoot)
library(esback)
library(forecast)

setwd("C:/Users/Tapan Kar/Dropbox/Crypto new/Journal of forecasting Revision 1/Re-submission")
data <- read.csv("Etherum_Binance_data.csv")
close= data$Close
open= data$Open
high= data$High
low= data$Low
l= length(close)

# parkinson volatilityestimator

pk_proxy= {(log(high)-log(low))^2}/(4*log(2))
pk_proxy= pk_proxy[-1]

# computing log returns
log_ret= c(rep(0,l))

for(i in 2:l){
  log_ret[i]= log(close[i])-log(close[(i-1)])
}

log_ret= log_ret[-1]

log_range= c(rep(0,l))

for(i in 1:l){
  log_range[i]= log(high[i])-log(low[(i)])
}

log_range= log_range[-1]
range= sqrt(log_range)

# Rollingwindow for Garch model 904 windows
m= 904
pre_garch_fast1=c(rep(0,m))
fit_mean_garch= c(rep(0,m))

for(q in 1:m){
  start= q
  end= 999+q
  gar_fast1= log_ret[start:end]
  fit <- auto.arima(gar_fast1,d=0,max.p = 2,max.q = 0)
  pre_mean_garch= (predict(fit, 1))
  res_garch= fit$residuals
  fit_mean_garch[q]= as.numeric(pre_mean_garch$pred)
  spec.gjrGARCH = ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0), include.mean=FALSE), distribution.model="ged")
  e_garch <- ugarchfit(res_garch, spec=spec.gjrGARCH)
  e_forecast= ugarchforecast(e_garch,n.ahead=1,data=res_garch)
  fore= e_forecast@forecast
  pre_garch_fast1[q]= as.numeric(fore$sigmaFor)
  
}

etherum_garch_mean= fit_mean_garch
etherum_garch= pre_garch_fast1

write.csv(etherum_garch_mean,"etherum_garch_mean.csv")
etherum_garch_mean= read.csv("etherum_garch_mean.csv")
etherum_garch_mean= etherum_garch_mean[,-1]
write.csv(etherum_garch,"etherum_garch.csv")
etherum_garch_vol= read.csv("etherum_garch.csv")
etherum_garch_vol= etherum_garch_vol[,-1]

# NAGARCH Model

pre_nagarch_binance= c(rep(0,m))

for(q in 1:m){
  start= q
  end= 999+q
  gar_fast1= log_ret[start:end]
  fit <- auto.arima(gar_fast1,d=0,max.p = 2,max.q = 0)
  pre_mean_garch= (predict(fit, 1))
  res_garch= fit$residuals
  spec.gjrGARCH = ugarchspec(variance.model=list(model="fGARCH", garchOrder=c(1,1),submodel= "NAGARCH"), mean.model=list(armaOrder=c(0,0), include.mean=FALSE), distribution.model="ged")
  e_garch <- ugarchfit(res_garch, spec=spec.gjrGARCH,solver = "hybrid")
  e_forecast= ugarchforecast(e_garch,n.ahead=1,data=res_garch)
  fore= e_forecast@forecast
  pre_nagarch_binance[q]= as.numeric(fore$sigmaFor)
}

bid_nagarch_binance_eth= pre_nagarch_binance
write.csv(bid_nagarch_binance_eth,"bid_nagarch_binance_eth.csv")
nagarch_binance_eth= read.csv("bid_nagarch_binance_eth.csv")
nagarch_binance_eth= nagarch_binance_eth[,-1]

#  Rolling window forecast for MSGARCH Model
msgarch_var_01_binance= c(rep(0,m))
msgarch_es_01_binance= c(rep(0,m))
msgarch_var_05_binance= c(rep(0,m))
msgarch_es_05_binance= c(rep(0,m))
msgarch_var_025_binance= c(rep(0,m))
msgarch_es_025_binance= c(rep(0,m))
eth_msgarch_vol_binance= c(rep(0,m))

for(q in 1:m){
  start= q
  end= 999+q
  gar_fast1= log_ret[start:end]
  fit <- auto.arima(gar_fast1,d=0,max.p = 2,max.q = 0)
  pre_mean_garch= (predict(fit, 1))
  res_garch= fit$residuals
  ms2.garch.n= CreateSpec(variance.spec = list(model = "sGARCH"),
                          distribution.spec = list(distribution = "ged"),
                          switch.spec = list(K = 2))
  fit.ml <- FitMCMC(spec = ms2.garch.n, data = res_garch)
  ms_garch= predict(fit.ml,nahead = 1)
  eth_msgarch_vol_binance[q]= ms_garch$vol
  
  e_forecast_01= Risk(fit.ml, alpha = 0.01, nahead = 1)
  e_forecast_05= Risk(fit.ml, alpha = 0.05, nahead = 1)
  e_forecast_025= Risk(fit.ml, alpha = 0.025, nahead = 1)
  msgarch_var_01_binance[q]= e_forecast_01$VaR
  msgarch_es_01_binance[q]= e_forecast_01$ES
  msgarch_var_05_binance[q]= e_forecast_05$VaR
  msgarch_es_05_binance[q]= e_forecast_05$ES
  msgarch_var_025_binance[q]= e_forecast_025$VaR
  msgarch_es_025_binance[q]= e_forecast_025$ES
  
}

eth_msgarch_var_01_binance=  msgarch_var_01_binance
eth_msgarch_es_01_binance=  msgarch_es_01_binance
eth_msgarch_var_05_binance=  msgarch_var_05_binance
eth_msgarch_es_05_binance=  msgarch_es_05_binance
eth_msgarch_vol_binance= eth_msgarch_vol_binance
eth_msgarch_var_025_binance=  msgarch_var_025_binance
eth_msgarch_es_025_binance=  msgarch_es_025_binance

write.csv(eth_msgarch_vol_binance,"eth_msgarch_vol_binance.csv")
write.csv(eth_msgarch_var_01_binance,"eth_msgarch_var_01_binance.csv")
write.csv(eth_msgarch_es_01_binance,"eth_msgarch_es_01_binance.csv")
write.csv(eth_msgarch_var_05_binance,"eth_msgarch_var_05_binance.csv")
write.csv(eth_msgarch_es_05_binance,"eth_msgarch_es_05_binance.csv")
write.csv(eth_msgarch_var_025_binance,"eth_msgarch_var_025_binance.csv")
write.csv(eth_msgarch_es_025_binance,"eth_msgarch_es_025_binance.csv")

eth_msgarch_var_01_binance= read.csv("eth_msgarch_var_01_binance.csv")
eth_msgarch_var_01_binance= eth_msgarch_var_01_binance[,-1]
eth_msgarch_vol_binance= read.csv("eth_msgarch_vol_binance.csv")
eth_msgarch_vol_binance= eth_msgarch_vol_binance[,-1]
eth_msgarch_es_01_binance= read.csv("eth_msgarch_es_01_binance.csv")
eth_msgarch_es_01_binance= eth_msgarch_es_01_binance[,-1]
eth_msgarch_var_05_binance= read.csv("eth_msgarch_var_05_binance.csv")
eth_msgarch_var_05_binance= eth_msgarch_var_05_binance[,-1]
eth_msgarch_es_05_binance= read.csv("eth_msgarch_es_05_binance.csv")
eth_msgarch_es_05_binance= eth_msgarch_es_05_binance[,-1]
eth_msgarch_var_025_binance= read.csv("eth_msgarch_var_025_binance.csv")
eth_msgarch_var_025_binance= eth_msgarch_var_025_binance[,-1]
eth_msgarch_es_025_binance= read.csv("eth_msgarch_es_025_binance.csv")
eth_msgarch_es_025_binance= eth_msgarch_es_025_binance[,-1]


# C-RV model with normal copula, scaling factor estimation
m1= 1000
insam_est_binance_eth= c(rep(0,m1))

train_fast_sam= pk_proxy[1:1000]
ecdf= pobs(train_fast_sam)
N= length(train_fast_sam)
#ecdf= (N/(N+1))*ecd
u1= ecdf[1:(N-1)] 
u2= ecdf[2:N]
u= cbind(u1,u2)
fit.ml= fitCopula(normalCopula(), u, method="mpl")

n <- 1000
for (q in 1:1000) {
  #q = 1
  u_1 <- ecdf[q]
  U <- cCopula(cbind(u_1, runif(n)), copula = fit.ml@copula, inverse = TRUE)
  rea= U[,2]
  my_ecdf_inv_1 = c(rep(0,n))
  for(i in 1:n){
    my_ecdf_inv_1[i] <- quantile(train_fast_sam, probs = rea[i])
  }
  insam_est_binance_eth[q]= mean(my_ecdf_inv_1)
}

insam_est_binance_eth= insam_est_binance_eth
write.csv(insam_est_binance_eth,"insam_est_binance_eth.csv")
insam_est_binance_eth= read.csv("insam_est_binance_eth.csv")
insam_est_binance_eth= insam_est_binance_eth[,-1]
scaling_factor_binance_eth= var(log_ret[1:1000])/mean(insam_est_binance_eth)

# C-RV model with normal copula rolling window forecasting
nor_mean_1= c(rep(0,m))

for(q in 1:m){
  start= q
  end= 999+q
  train_fast1= pk_proxy[start:end]
  ecdf= pobs(train_fast1)
  N= length(train_fast1)
  #ecdf= (N/(N+1))*ecd
  u1= ecdf[1:(N-1)] 
  u2= ecdf[2:N]
  u= cbind(u1,u2)
  fit.ml= fitCopula(normalCopula(), u, method="mpl")
  n <- 1000
  u_1 <- ecdf[1000]
  U <- cCopula(cbind(u_1, runif(n)), copula = fit.ml@copula, inverse = TRUE)
  rea= U[,2]
  my_ecdf_inv_1 = c(rep(0,n))
  
  for(i in 1:n){
    my_ecdf_inv_1[i] <- quantile(train_fast1, probs = rea[i])
  }
  nor_mean_1[q]= mean(my_ecdf_inv_1)
}

etherum_vol_nor= (nor_mean_1)
write.csv(etherum_vol_nor,"etherum_vol_nor_binance.csv")
etherum_vol_nor_binance= read.csv("etherum_vol_nor_binance.csv")
etherum_vol_nor_binance= etherum_vol_nor_binance[,-1]
adj_eth_vol_nor_binance= scaling_factor_binance_eth*(etherum_vol_nor_binance)
adj_eth_vol_nor_binance= sqrt(adj_eth_vol_nor_binance)


# scaling factor for C-RV model with t-copula

m1= 1000
eth_insam_est_t_binance= c(rep(0,m1))
train_fast_sam= pk_proxy[1:1000]
ecdf= pobs(train_fast_sam)
N= length(train_fast_sam)
#ecdf= (N/(N+1))*ecd
u1_t= ecdf[1:(N-1)] 
u2_t= ecdf[2:N]
u_t= cbind(u1_t,u2_t)
fit.ml= fitCopula(tCopula(), u_t, method="mpl",estimate.variance=FALSE)

n <- 1000
for (q in 1:1000) {
  
  u_t <- ecdf[q]
  U_t <- cCopula(cbind(u_t, runif(n)), copula = fit.ml@copula, inverse = TRUE)
  rea= U_t[,2]
  my_ecdf_inv_t = c(rep(0,n))
  for(i in 1:n){
    my_ecdf_inv_t[i] <- quantile(train_fast_sam, probs = rea[i])
  }
  eth_insam_est_t_binance[q]= mean(my_ecdf_inv_t)
}

eth_insam_est_t_binance= eth_insam_est_t_binance
write.csv(eth_insam_est_t_binance,"eth_insam_est_t_binance.csv")
eth_insam_est_t_binance= read.csv("eth_insam_est_t_binance.csv")
eth_insam_est_t_binance= eth_insam_est_t_binance[,-1]
eth_scaling_factor_t_binance= var(log_ret[1:1000])/mean(eth_insam_est_t_binance)

# C-RV model with t copula rolling window forecasting
t_mean_1= c(rep(0,m))

for(q in 1:m){
  start= q
  end= 999+q
  train_fast1= pk_proxy[start:end]
  ecdf= pobs(train_fast1)
  N= length(train_fast1)
  u1= ecdf[1:(N-1)]
  u2= ecdf[2:N]
  u= cbind(u1,u2)
  fit.ml= fitCopula(tCopula(), u, method="mpl",estimate.variance=FALSE)
  n <- 1000
  u_1 <- ecdf[1000]
  U <- cCopula(cbind(u_1, runif(n)), copula = fit.ml@copula, inverse = TRUE)
  rea= U[,2]
  my_ecdf_inv_1 = c(rep(0,n))
  
  for(i in 1:n){
    my_ecdf_inv_1[i] <- quantile(train_fast1, probs = rea[i])
  }
  t_mean_1[q]= mean(my_ecdf_inv_1)
}

etherum_t_vol= (t_mean_1)
write.csv(etherum_t_vol,"etherum_t_vol_binance.csv")
etherum_t_vol= read.csv("etherum_t_vol_binance.csv")
etherum_t_vol= etherum_t_vol[,-1]
etherum_t_vol= sqrt(eth_scaling_factor_t_binance*(etherum_t_vol))

# RobGARCHBoot model forecasting

library(RobGARCHBoot)

VaR_rob_01_binance= c(rep(0,m))
VaR_rob_05_binance= c(rep(0,m))
VaR_rob_025_binance= c(rep(0,m))

ES_rob_01_binance= c(rep(0,m))
ES_rob_05_binance= c(rep(0,m))
ES_rob_025_binance= c(rep(0,m))
sq_vol_rob_binance= c(rep(0,m))

for(q in 1:m){
  start=q
  end= 999+q
  data_test= log_ret[start:end]
  fit <- auto.arima(data_test,d=0,max.p = 2,max.q = 0)
  pre_mean_garch= (predict(fit, 1))
  res_garch= fit$residuals
  boot = RobGARCHBoot(res_garch, n.boot = 1000, n.ahead = 1)
  return= boot[[1]]
  sq_vol_rob_binance[q]= mean(boot[[2]])
  VaR_rob_01_binance[q] = quantile(boot[[1]], prob = 0.01)
  VaR_rob_05_binance[q] = quantile(boot[[1]], prob = 0.05)
  VaR_rob_025_binance[q] = quantile(boot[[1]], prob = 0.025)
  ES_rob_01_binance[q]= mean(return[return<VaR_rob_01_binance[q]])
  ES_rob_05_binance[q]= mean(return[return<VaR_rob_05_binance[q]])
  ES_rob_025_binance[q]= mean(return[return<VaR_rob_025_binance[q]])
}

eth_robgarch_var_01_binance= (VaR_rob_01_binance)
eth_robgarch_var_05_binance= (VaR_rob_05_binance)
eth_robgarch_var_025_binance= (VaR_rob_025_binance)
eth_robgarch_ES_01_binance= (ES_rob_01_binance)
eth_robgarch_ES_05_binance= (ES_rob_05_binance)
eth_robgarch_ES_025_binance= (ES_rob_025_binance)
eth_sq_vol_rob_binance= sq_vol_rob_binance

write.csv(eth_robgarch_ES_01_binance,"eth_robgarch_ES_01_binance.csv")
write.csv(eth_robgarch_ES_05_binance,"eth_robgarch_ES_05_binance.csv")
write.csv(eth_robgarch_ES_025_binance,"eth_robgarch_ES_025_binance.csv")
write.csv(eth_sq_vol_rob_binance,"eth_sq_vol_rob_binance.csv")
write.csv(eth_robgarch_var_01_binance,"eth_robgarch_var_01_binance.csv")
write.csv(eth_robgarch_var_05_binance,"eth_robgarch_var_05_binance.csv")
write.csv(eth_robgarch_var_025_binance,"eth_robgarch_var_025_binance.csv")

eth_robgarch_var_01_binance= read.csv("eth_robgarch_var_01_binance.csv")
eth_robgarch_var_01_binance= eth_robgarch_var_01_binance[,-1]
eth_robgarch_var_05_binance= read.csv("eth_robgarch_var_05_binance.csv")
eth_robgarch_var_05_binance= eth_robgarch_var_05_binance[,-1]
eth_sq_vol_rob_binance= read.csv("eth_sq_vol_rob_binance.csv")
eth_sq_vol_rob_binance= eth_sq_vol_rob_binance[,-1]
eth_sq_vol_rob_binance= sqrt(eth_sq_vol_rob_binance)
eth_robgarch_var_025_binance= read.csv("eth_robgarch_var_025_binance.csv")
eth_robgarch_var_025_binance= eth_robgarch_var_025_binance[,-1]

eth_robgarch_ES_01_binance= read.csv("eth_robgarch_ES_01_binance.csv")
eth_robgarch_ES_01_binance= eth_robgarch_ES_01_binance[,-1]
eth_robgarch_ES_05_binance= read.csv("eth_robgarch_ES_05_binance.csv")
eth_robgarch_ES_05_binance= eth_robgarch_ES_05_binance[,-1]
eth_robgarch_ES_025_binance= read.csv("eth_robgarch_ES_025_binance.csv")
eth_robgarch_ES_025_binance= eth_robgarch_ES_025_binance[,-1]

#CARR model scaling factor
in_sam_ran= range[1:1000]
spec.gjrGARCH_range = ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0), include.mean=FALSE),distribution.model="norm")
e_garch_range <- ugarchfit(in_sam_ran, spec=spec.gjrGARCH_range)
est_range= (e_garch_range@fit$sigma)^2
scaling_carr_binance= sd(log_ret[1:1000])/mean(est_range)

# CARR model forecasting
carr_vol_binance= c(rep(0,m))

for(q in 1:m){
  start= q
  end= 999+q
  train_fast1= range[start:end]
  spec.gjrGARCH_range = ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0), include.mean=FALSE),distribution.model="norm")
  e_garch_range <- ugarchfit(train_fast1, spec=spec.gjrGARCH_range)
  predi_garch_range= ugarchforecast(e_garch_range,n.ahead=1,data=train_fast1)
  fore= (predi_garch_range@forecast)
  carr_vol_binance[q]= as.numeric(fore$sigma)
}

etherum_carr_vol_binance= (carr_vol_binance)
write.csv(etherum_carr_vol_binance,"etherum_carr_vol_binance.csv")
etherum_carr_vol_binance= read.csv("etherum_carr_vol_binance.csv")
etherum_carr_vol_binance= etherum_carr_vol_binance[,-1]
etherum_carr_vol_binance= (etherum_carr_vol_binance)^2
adj_etherum_carr_vol_binance= scaling_carr_binance*(etherum_carr_vol_binance)

# parameter estimation of the innovations with ged distribution

fit_ar= auto.arima(log_ret[1:1000],d=0,max.p = 2,max.q = 0)
spec.gjrGARCH = ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0), include.mean=FALSE),distribution.model="ged")
e_garch <- ugarchfit(fit_ar$residuals, spec=spec.gjrGARCH)


# VaR for GARCH model

VaR_garch =  (etherum_garch_vol)*qdist("ged", p=0.05, mu = 0, sigma = 1, 
                                       shape=coef(e_garch)["shape"])

# VaR for C-RVN model
VaR_nor_cop =  (adj_eth_vol_nor_binance)*qdist("ged", p=0.05, mu = 0, sigma = 1, 
                                               shape=coef(e_garch)["shape"])

# VaR for C-RVt model
VaR_t_cop =  (etherum_t_vol)*qdist("ged", p=0.05, mu = 0, sigma = 1, 
                                   shape=coef(e_garch)["shape"])

# VaR for CARR model

VaR_carr =   (adj_etherum_carr_vol_binance)*qdist("ged", p=0.05, mu = 0, sigma = 1, 
                                                  shape=coef(e_garch)["shape"])

# VaR for NAGRCH model
VaR_nagarch =  (nagarch_binance_eth)*qdist("ged", p=0.05, mu = 0, sigma = 1, 
                                           shape=coef(e_garch)["shape"])

test= log_ret[1001:1904]
actual= test-etherum_garch_mean

# Backtesting, UC,CC, an DC tests
a= BacktestVaR(data = actual, VaR = VaR_nor_cop, alpha = 0.05)
a1= BacktestVaR(data = actual, VaR = VaR_garch, alpha = 0.05)
a5= BacktestVaR(data = actual, VaR = VaR_nagarch, alpha = 0.05)
a2= BacktestVaR(data = actual, VaR = VaR_carr, alpha = 0.05)
a3= BacktestVaR(data = actual, VaR = VaR_t_cop, alpha = 0.05)
a4= BacktestVaR(data = actual, VaR = eth_robgarch_var_05_binance, alpha = 0.05)
a7= BacktestVaR(data = actual, VaR = eth_msgarch_var_05_binance, alpha = 0.05)

# Expected shortfall McF test GARH model 


f1 = function(x) qdist("ged", p=x, mu = 0, sigma = 1,
                       shape=coef(e_garch)["shape"])
ES_garch = (etherum_garch_vol)*integrate(f1, 0, 0.05)$value/0.05
print(ESTest(0.05, actual, ES_garch, VaR_garch, boot = TRUE))

#Expected shortfall McF test NAGARCH model 
f5 = function(x) qdist("ged", p=x, mu = 0, sigma = 1,
                       shape=coef(e_garch)["shape"])
ES_nagarch = (nagarch_binance_eth)*integrate(f5, 0, 0.05)$value/0.05
print(ESTest(0.05, actual, ES_nagarch, VaR_nagarch, boot = TRUE))

#Expected shortfall McF test C-RVN model 
f2 = function(x) qdist("ged", p=x, mu = 0, sigma = 1,
                       shape=coef(e_garch)["shape"])
ES_nor_vol = (adj_eth_vol_nor_binance)*integrate(f2, 0, 0.05)$value/0.05
print(ESTest(0.05, actual, ES_nor_vol, VaR_nor_cop, boot = TRUE))

#Expected shortfall McF test C-RVt model 
f4 = function(x) qdist("ged", p=x, mu = 0, sigma = 1,
                       shape=coef(e_garch)["shape"])
ES_t_vol = (etherum_t_vol)*integrate(f4, 0, 0.05)$value/0.05
print(ESTest(0.05, actual, ES_t_vol,VaR_t_cop, boot = TRUE))

#Expected shortfall McF test RobGARCHBoot model 
print(ESTest(0.05, actual, eth_robgarch_ES_05_binance, eth_robgarch_var_05_binance, boot = TRUE))

#Expected shortfall McF test MSGARCH model 
print(ESTest(0.05, actual, eth_msgarch_es_05_binance, eth_msgarch_var_05_binance, boot = TRUE))

#Expected shortfall McF test CARR model 
f3 = function(x) qdist("ged", p=x, mu = 0, sigma = 1,
                       shape=coef(e_garch)["shape"])
ES_carr =  (adj_etherum_carr_vol_binance)*integrate(f3, 0, 0.05)$value/0.05
print(ESTest(0.05, actual, ES_carr, VaR_carr, boot = TRUE))


library(esback)

#NF test
cc_backtest(actual, VaR_nor_cop, ES_nor_vol, s = adj_eth_vol_nor_binance, alpha=0.05)
cc_backtest(actual, VaR_t_cop, ES_t_vol, s = etherum_t_vol, alpha=0.05)
cc_backtest(actual, VaR_garch, ES_garch, s = etherum_garch_vol,alpha = 0.05)
cc_backtest(actual, VaR_nagarch, ES_nagarch, s = nagarch_binance_eth,alpha = 0.05)
cc_backtest(actual, VaR_carr, ES_carr, s = adj_etherum_carr_vol_binance, alpha=0.05)
cc_backtest(actual, eth_robgarch_var_05_binance, eth_robgarch_ES_05_binance, s =eth_sq_vol_rob_binance, alpha=0.05)
cc_backtest(actual, eth_msgarch_var_05_binance, eth_msgarch_es_05_binance, s =eth_msgarch_vol_binance, alpha=0.05)

# BD test
esr_backtest( actual,VaR_nor_cop,ES_nor_vol,alpha = 0.05, version = 1)
esr_backtest( actual,VaR_t_cop,ES_t_vol,alpha = 0.05, version = 1)
esr_backtest( actual,eth_msgarch_var_05_binance, eth_msgarch_es_05_binance,alpha = 0.05, version = 1)
esr_backtest( actual,VaR_nagarch, ES_nagarch,alpha = 0.05, version = 1)
esr_backtest( actual,VaR_garch, ES_garch,alpha = 0.05, version = 1)
esr_backtest( actual,VaR_carr, ES_carr,alpha = 0.05, version = 1)
esr_backtest( actual,eth_robgarch_var_05_binance, eth_robgarch_ES_05_binance,alpha = 0.05, version = 1)
