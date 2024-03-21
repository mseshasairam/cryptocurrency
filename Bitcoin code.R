rm(list=ls()) 
library(dplyr)
library(MSGARCH)
library(rugarch)
library(copula)
library(GAS)
library(RobGARCHBoot)
library(esback)

setwd("C:/Users/Tapan Kar/Dropbox/Crypto new/Journal of forecasting Revision 1/Re-submission")

data <- read.csv("Bitcoin_Binance_data.csv")
l= length(data$Close)

# computing log returns
log_ret= c(rep(0,l))

for(i in 2:l){
  log_ret[i]= log(data$Close[i]/data$Close[i-1])
}
log_ret= log_ret[-1]
log_ret= log_ret-mean(log_ret)
high= data$High
low= data$Low

# Parkinson volatility estimator

pk_proxy= {(log(high)-log(low))^2}/(4*log(2))
pk_proxy= pk_proxy[-1]

# Range estimator for CARR model
log_range= c(rep(0,l))

for(i in 1:l){
  log_range[i]= log(high[i])-log(low[(i)])
}

log_range=log_range[-1]
range= sqrt(log_range)

# Rolling window forecast for GARCH Model for 904 windows

m= 904
pre_garch_binance=c(rep(0,m))

for(q in 1:m){
  start= q
  end= 999+q
  gar_fast1= log_ret[start:end]
  spec.gjrGARCH = ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0), include.mean=FALSE), distribution.model="ged")
  e_garch <- ugarchfit(gar_fast1, spec=spec.gjrGARCH)
  e_forecast= ugarchforecast(e_garch,n.ahead=1,data=gar_fast1)
  fore= e_forecast@forecast
  pre_garch_binance[q]= as.numeric(fore$sigmaFor)
}

bid_garch_binance= pre_garch_binance
write.csv(bid_garch_binance,"bid_garch_binance.csv")
bid_garch_binance= read.csv("bid_garch_binance.csv")
bid_garch_binance= bid_garch_binance[,-1]

# Rolling window forecast for NAGARCH Model
pre_nagarch_binance= c(rep(0,m))

for(q in 1:m){
  start= q
  end= 999+q
  gar_fast1= log_ret[start:end]
  spec.gjrGARCH = ugarchspec(variance.model=list(model="fGARCH", garchOrder=c(1,1),submodel= "NAGARCH"), mean.model=list(armaOrder=c(0,0), include.mean=FALSE), distribution.model="ged")
  e_garch <- ugarchfit(gar_fast1, spec=spec.gjrGARCH,solver = "hybrid")
  e_forecast= ugarchforecast(e_garch,n.ahead=1,data=gar_fast1)
  fore= e_forecast@forecast
  pre_nagarch_binance[q]= as.numeric(fore$sigmaFor)
}

bid_nagarch_binance= pre_nagarch_binance
write.csv(bid_nagarch_binance,"bid_nagarch_binance.csv")
bid_nagarch_binance= read.csv("bid_nagarch_binance.csv")
bid_nagarch_binance= bid_nagarch_binance[,-1]

#  Rolling window forecast for MSGARCH Model
msgarch_var_01_binance= c(rep(0,m))
msgarch_es_01_binance= c(rep(0,m))
msgarch_var_05_binance= c(rep(0,m))
msgarch_es_05_binance= c(rep(0,m))
msgarch_var_025_binance= c(rep(0,m))
msgarch_es_025_binance= c(rep(0,m))
msgarch_vol_binance= c(rep(0,m))

for(q in 1:m){
  start= q
  end= 999+q
  gar_fast1= log_ret[start:end]
  ms2.garch.n= CreateSpec(variance.spec = list(model = "sGARCH"),
                          distribution.spec = list(distribution = "ged"),
                          switch.spec = list(K = 2))
  fit.ml <- FitMCMC(spec = ms2.garch.n, data = gar_fast1)
  ms_garch= predict(fit.ml,nahead = 1)
  
  msgarch_vol_binance[q]= ms_garch$vol
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

bid_msgarch_var_01_binance=  msgarch_var_01_binance
bid_msgarch_es_01_binance=  msgarch_es_01_binance
bid_msgarch_var_05_binance=  msgarch_var_05_binance
bid_msgarch_es_05_binance=  msgarch_es_05_binance
msgarch_vol_binance= msgarch_vol_binance
bid_msgarch_var_025_binance=  msgarch_var_025_binance
bid_msgarch_es_025_binance=  msgarch_es_025_binance

write.csv(msgarch_vol_binance,"msgarch_vol_binance.csv")
write.csv(bid_msgarch_var_01_binance,"bid_msgarch_var_01_binance.csv")
write.csv(bid_msgarch_es_01_binance,"bid_msgarch_es_01_binance.csv")
write.csv(bid_msgarch_var_05_binance,"bid_msgarch_var_05_binance.csv")
write.csv(bid_msgarch_es_05_binance,"bid_msgarch_es_05_binance.csv")
write.csv(bid_msgarch_var_025_binance,"bid_msgarch_var_025_binance.csv")
write.csv(bid_msgarch_es_025_binance,"bid_msgarch_es_025_binance.csv")

bid_msgarch_var_01_binance= read.csv("bid_msgarch_var_01_binance.csv")
bid_msgarch_var_01_binance= bid_msgarch_var_01_binance[,-1]
msgarch_vol_binance= read.csv("msgarch_vol_binance.csv")
msgarch_vol_binance= msgarch_vol_binance[,-1]
bid_msgarch_es_01_binance= read.csv("bid_msgarch_es_01_binance.csv")
bid_msgarch_es_01_binance= bid_msgarch_es_01_binance[,-1]
bid_msgarch_var_05_binance= read.csv("bid_msgarch_var_05_binance.csv")
bid_msgarch_var_05_binance= bid_msgarch_var_05_binance[,-1]
bid_msgarch_es_05_binance= read.csv("bid_msgarch_es_05_binance.csv")
bid_msgarch_es_05_binance= bid_msgarch_es_05_binance[,-1]
bid_msgarch_var_025_binance= read.csv("bid_msgarch_var_025_binance.csv")
bid_msgarch_var_025_binance= bid_msgarch_var_025_binance[,-1]
bid_msgarch_es_025_binance= read.csv("bid_msgarch_es_025_binance.csv")
bid_msgarch_es_025_binance= bid_msgarch_es_025_binance[,-1]

# In sample estimate For the C-RV with normal copula and the scaling factor estimation

m1= 1000
insam_est_binance= c(rep(0,m1))
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
  insam_est_binance[q]= mean(my_ecdf_inv_1)
}

insam_est_binance= insam_est_binance
write.csv(insam_est_binance,"insam_est_binance.csv")
insam_est_binance= read.csv("insam_est_binance.csv")
insam_est_binance= insam_est_binance[,-1]
scaling_factor_binance= var(log_ret[1:1000])/mean(insam_est_binance)

# Rolling window forecast for C-RV Model with normal copula

nor_mean_binance= c(rep(0,m))

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
  nor_mean_binance[q]= mean(my_ecdf_inv_1)
}

bid_vol_nor_binance= (nor_mean_binance)
write.csv(bid_vol_nor_binance,"bid_vol_nor_binance.csv")
bid_vol_nor_binance= read.csv("bid_vol_nor_binance.csv")
bid_vol_nor_binance= bid_vol_nor_binance[,-1]
adj_bid_vol_nor_binance= scaling_factor_binance*(bid_vol_nor_binance)
adj_bid_vol_nor_binance= sqrt(adj_bid_vol_nor_binance)

# In sample estimate For C-RV with t copula and the scaling factor estimation

m1= 1000
insam_est_t_binance= c(rep(0,m1))
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
  insam_est_t_binance[q]= mean(my_ecdf_inv_t)
}

insam_est_t_binance= insam_est_t_binance
write.csv(insam_est_t_binance,"insam_est_t_binance.csv")
insam_est_t_binance= read.csv("insam_est_t_binance.csv")
insam_est_t_binance= insam_est_t_binance[,-1]
scaling_factor_t_binance= var(log_ret[1:1000])/mean(insam_est_t_binance)

# Rolling window forecast for C-RV Model with t copula
t_mean_1_binance= c(rep(0,m))

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
  t_mean_1_binance[q]= mean(my_ecdf_inv_1)
}

bid_t_vol_binance= (t_mean_1_binance)
write.csv(bid_t_vol_binance,"bid_t_vol_binance.csv")
bid_t_vol_binance= read.csv("bid_t_vol_binance.csv")
bid_t_vol_binance= bid_t_vol_binance[,-1]
bid_t_vol_binance= (scaling_factor_t_binance*bid_t_vol_binance)
adj_bid_vol_t_binance= sqrt(bid_t_vol_binance)

# in sample estimate CARR model and scaling factor

in_sam_ran= range[1:1000]
spec.gjrGARCH_range = ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0), include.mean=FALSE),distribution.model="norm")
e_garch_range <- ugarchfit(in_sam_ran, spec=spec.gjrGARCH_range)
est_range= (e_garch_range@fit$sigma)^2
scaling_carr= sd(log_ret[1:1000])/mean(est_range)

# Rolling window forecast by CARR model
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

bid_carr_vol_binance= (carr_vol_binance)
write.csv(bid_carr_vol_binance,"bid_carr_vol_binance.csv")
bid_carr_vol_binance= read.csv("bid_carr_vol_binance.csv")
bid_carr_vol_binance= bid_carr_vol_binance[,-1]
bid_carr_vol_binance= (bid_carr_vol_binance)^2
adj_bid_carr_vol_binance= scaling_carr*(bid_carr_vol_binance)

#  Rolling window forecast for RobGARCHBoot Model

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
  boot = RobGARCHBoot(data_test, n.boot = 1000, n.ahead = 1)
  return= boot[[1]]
  sq_vol_rob_binance[q]= mean(boot[[2]])
  VaR_rob_01_binance[q] = quantile(boot[[1]], prob = 0.01)
  VaR_rob_05_binance[q] = quantile(boot[[1]], prob = 0.05)
  VaR_rob_025_binance[q] = quantile(boot[[1]], prob = 0.025)
  
  ES_rob_01_binance[q]= mean(return[return<VaR_rob_01_binance[q]])
  ES_rob_05_binance[q]= mean(return[return<VaR_rob_05_binance[q]])
  ES_rob_025_binance[q]= mean(return[return<VaR_rob_025_binance[q]])
}

bid_robgarch_var_01_binance= (VaR_rob_01_binance)
bid_robgarch_var_05_binance= (VaR_rob_05_binance)
bid_robgarch_var_025_binance= (VaR_rob_025_binance)
bid_robgarch_ES_01_binance= (ES_rob_01_binance)
bid_robgarch_ES_05_binance= (ES_rob_05_binance)
bid_robgarch_ES_025_binance= (ES_rob_025_binance)
bid_sq_vol_rob_binance= sq_vol_rob_binance

write.csv(bid_robgarch_ES_01_binance,"bid_robgarch_ES_01_binance.csv")
write.csv(bid_robgarch_ES_05_binance,"bid_robgarch_ES_05_binance.csv")
write.csv(bid_robgarch_ES_025_binance,"bid_robgarch_ES_025_binance.csv")
write.csv(bid_sq_vol_rob_binance,"bid_sq_vol_rob_binance.csv")

write.csv(bid_robgarch_var_01_binance,"bid_robgarch_var_01_binance.csv")
write.csv(bid_robgarch_var_05_binance,"bid_robgarch_var_05_binance.csv")
write.csv(bid_robgarch_var_025_binance,"bid_robgarch_var_025_binance.csv")

bid_sq_vol_rob_binance= read.csv("bid_sq_vol_rob_binance.csv")
bid_sq_vol_rob_binance= bid_sq_vol_rob_binance[,-1]
bid_sq_vol_rob_binance= sqrt(bid_sq_vol_rob_binance)

bid_robgarch_ES_01_binance= read.csv("bid_robgarch_ES_01_binance.csv")
bid_robgarch_ES_01_binance= bid_robgarch_ES_01_binance[,-1]
bid_robgarch_ES_05_binance= read.csv("bid_robgarch_ES_05_binance.csv")
bid_robgarch_ES_05_binance= bid_robgarch_ES_05_binance[,-1]
bid_robgarch_ES_025_binance= read.csv("bid_robgarch_ES_025_binance.csv")
bid_robgarch_ES_025_binance= bid_robgarch_ES_025_binance[,-1]

bid_robgarch_var_025_binance= read.csv("bid_robgarch_var_025_binance.csv")
bid_robgarch_var_025_binance= bid_robgarch_var_025_binance[,-1]
bid_robgarch_var_01_binance= read.csv("bid_robgarch_var_01_binance.csv")
bid_robgarch_var_01_binance= bid_robgarch_var_01_binance[,-1]
bid_robgarch_var_05_binance= read.csv("bid_robgarch_var_05_binance.csv")
bid_robgarch_var_05_binance= bid_robgarch_var_05_binance[,-1]

# the parameterestimation of the innovations with ged distribution

spec.gjrGARCH = ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)), mean.model=list(armaOrder=c(0,0), include.mean=FALSE),distribution.model="ged")
e_garch <- ugarchfit(log_ret[1:1000], spec=spec.gjrGARCH)


# VaR for GARCH model
VaR_garch_binance = (bid_garch_binance)*qdist("ged", p=0.05, mu = 0, sigma = 1, 
                                               shape=coef(e_garch)["shape"])
# VaR for NAGARCH model
VaR_nagarch_binance = (bid_nagarch_binance)*qdist("ged", p=0.05, mu = 0, sigma = 1,
                                                   shape=coef(e_garch)["shape"])
# VaR for C-RV with normal copula model
VaR_nor_cop_binance = (adj_bid_vol_nor_binance)*qdist("ged", p=0.05, mu = 0, sigma = 1, 
                                                       shape=coef(e_garch)["shape"])
# VaR for C-RV with t copula model
VaR_t_cop_binance = (adj_bid_vol_t_binance)*qdist("ged", p=0.05, mu = 0, sigma = 1, 
                                                   shape=coef(e_garch)["shape"])
# VaR for CARR model
VaR_carr_binance = (adj_bid_carr_vol_binance)*qdist("ged", p=0.05, mu = 0, sigma = 1, 
                                                     shape=coef(e_garch)["shape"])
actual_binance= log_ret[1001:1904]


# Backtesting UC,CC and DQ tests
a= BacktestVaR(data = actual_binance, VaR = VaR_nor_cop_binance, alpha = 0.05)
a1= BacktestVaR(data = actual_binance, VaR = VaR_garch_binance, alpha = 0.05)
a2= BacktestVaR(data = actual_binance, VaR = VaR_carr_binance, alpha = 0.05)
a3= BacktestVaR(data = actual_binance, VaR = VaR_t_cop_binance, alpha = 0.05)
a4= BacktestVaR(data = actual_binance, VaR = bid_robgarch_var_05_binance, alpha = 0.05)
a5= BacktestVaR(data = actual_binance, VaR = VaR_nagarch_binance, alpha = 0.05)
a6= BacktestVaR(data = actual_binance, VaR = bid_msgarch_var_05_binance, alpha = 0.05)

# Expected shortfall McF test GARCH model 
f1 = function(x) qdist("ged", p=x, mu = 0, sigma = 1,
                       shape=coef(e_garch)["shape"])

ES_garch_binance = (bid_garch_binance)*integrate(f1, 0, 0.05)$value/0.05
print(ESTest(0.05, actual_binance, ES_garch_binance, VaR_garch_binance, boot = TRUE))

# Expected shortfall McF test C-RVN model 
f2 = function(x) qdist("ged", p=x, mu = 0, sigma = 1,
                       shape=coef(e_garch)["shape"])
ES_nor_vol_binance = (adj_bid_vol_nor_binance)*integrate(f2, 0, 0.05)$value/0.05
print(ESTest(0.05, actual_binance, ES_nor_vol_binance, VaR_nor_cop_binance, boot = TRUE))

# Expected shortfall McF test C-RVt model 
f4 = function(x) qdist("ged", p=x, mu = 0, sigma = 1,
                       shape=coef(e_garch)["shape"])
ES_t_vol_binance = (adj_bid_vol_t_binance)*integrate(f4, 0, 0.05)$value/0.05
print(ESTest(0.05, actual_binance, ES_t_vol_binance,VaR_t_cop_binance, boot = TRUE))

# Expected shortfall McF test RobGARCHBoot model 
print(ESTest(0.05, actual_binance, bid_robgarch_ES_05_binance, bid_robgarch_var_05_binance, boot = TRUE))

# Expected shortfall McF test CARR model 
f3 = function(x) qdist("ged", p=x, mu = 0, sigma = 1,
                       shape=coef(e_garch)["shape"])
ES_carr_binance = (adj_bid_carr_vol_binance)*integrate(f3, 0, 0.05)$value/0.05
print(ESTest(0.05, actual_binance, ES_carr_binance, VaR_carr_binance, boot = TRUE))

# Expected shortfall McF test NAGARCH model 
f5 = function(x) qdist("ged", p=x, mu = 0, sigma = 1,
                       shape=coef(e_garch)["shape"])
ES_nagarch_binance = (bid_nagarch_binance)*integrate(f5, 0, 0.05)$value/0.05
print(ESTest(0.05, actual_binance, ES_nagarch_binance, VaR_nagarch_binance, boot = TRUE))

# Expected shortfall McF test MSGARCH model 
print(ESTest(0.05, actual_binance, bid_msgarch_es_05_binance, bid_msgarch_var_05_binance, boot = TRUE))


#NF test
cc_backtest(actual_binance, VaR_nor_cop_binance, ES_nor_vol_binance,s= adj_bid_vol_nor_binance, alpha=0.05)
cc_backtest(actual_binance, VaR_t_cop_binance, ES_t_vol_binance, s = bid_t_vol_binance, alpha=0.05)
cc_backtest(actual_binance, VaR_garch_binance, ES_garch_binance,s=bid_garch_binance,alpha = 0.05)
cc_backtest(actual_binance, VaR_nagarch_binance, ES_nagarch_binance, s = bid_nagarch_binance,alpha = 0.05)
cc_backtest(actual_binance, VaR_carr_binance, ES_carr_binance, s = bid_carr_vol_binance, alpha=0.05)
cc_backtest(actual_binance, bid_robgarch_var_05_binance, bid_robgarch_ES_05_binance, s =bid_sq_vol_rob_binance, alpha=0.05)
cc_backtest(actual_binance, bid_msgarch_var_05_binance, bid_msgarch_es_05_binance, s =msgarch_vol_binance, alpha=0.05)

# BD test
esr_backtest( actual_binance,VaR_nor_cop_binance,ES_nor_vol_binance,alpha = 0.05, version = 1)
esr_backtest( actual_binance,VaR_t_cop_binance,ES_t_vol_binance,alpha = 0.05, version = 1)
esr_backtest( actual_binance,bid_msgarch_var_05_binance, bid_msgarch_es_05_binance,alpha = 0.05, version = 1)
esr_backtest( actual_binance,VaR_nagarch_binance, ES_nagarch_binance,alpha = 0.05, version = 1)
esr_backtest( actual_binance,VaR_garch_binance, ES_garch_binance,alpha = 0.05, version = 1)
esr_backtest( actual_binance,VaR_carr_binance, ES_carr_binance,alpha = 0.05, version = 1)
esr_backtest( actual_binance,bid_robgarch_var_05_binance, bid_robgarch_ES_05_binance,alpha = 0.05, version = 1)

