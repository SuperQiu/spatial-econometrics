############lecture des données########################""""
library(spdep)

example(baltimore)
###############################
# spatial object and neighborhood matrix
coord=cbind(baltimore$X,baltimore$Y)
balt.sp=  SpatialPoints(coord)
balt.spdf= SpatialPointsDataFrame(balt.sp, baltimore)
W.knn=knearneigh(coord,k=5,longlat = NULL)
W.nb=knn2nb(W.knn)
#W.mat=nb2mat(W.nb)
W.listw=nb2listw(W.nb)
library(GeoXp)
histnbmap(balt.spdf,W.nb,coord)
###############################
# OLS  model fitting
###############################

ols_mod=lm(PRICE~NROOM + DWELL + NBATH + PATIO + FIREPL + AC + BMENT + NSTOR+GAR+AGE+CITCOU+LOTSZ+SQFT,data=baltimore)
#keep only significant variables
ols_mod_2=lm(PRICE~ NBATH + PATIO + FIREPL + AC + BMENT + GAR+CITCOU,data=baltimore)
lm.morantest(ols_mod_2,W.listw, zero.policy=FALSE, alternative = "greater")
###############################
# spatially lagged variables construction
###############################
baltimore$NROOM_lag=lag.listw(W.listw,baltimore$NROOM)
baltimore$NBATH_lag=lag.listw(W.listw,baltimore$NBATH)
baltimore$PATIO_lag=lag.listw(W.listw,baltimore$PATIO)
baltimore$FIREPL_lag=lag.listw(W.listw,baltimore$FIREPL)
baltimore$AC_lag=lag.listw(W.listw,baltimore$AC)
baltimore$BMENT_lag=lag.listw(W.listw,baltimore$BMENT)
baltimore$GAR_lag=lag.listw(W.listw,baltimore$GAR)
baltimore$CITCOU_lag=lag.listw(W.listw,baltimore$CITCOU)

 ###############################
# OLS  model fitting    with lagged X
###############################

ols_croise=lm(PRICE~ NBATH + PATIO + FIREPL + AC + BMENT + GAR+ CITCOU+ FIREPL_lag + AC_lag + NBATH_lag + PATIO_lag +BMENT_lag + GAR_lag + CITCOU_lag,data=baltimore)
ols_croise_2=lm(PRICE~ NBATH + PATIO + FIREPL + AC + BMENT + GAR+ CITCOU+ FIREPL_lag + AC_lag ,data=baltimore)

lm.morantest(ols_croise_2,W.listw, zero.policy=FALSE, alternative = "greater")


 ###############################
# Spatial  models fitting
###############################

lagmodel <- lagsarlm(PRICE~ NBATH + PATIO + FIREPL + AC + BMENT + GAR+CITCOU,data=baltimore,listw=W.listw)
lagmodel2 <- lagsarlm(PRICE~ NBATH + PATIO + FIREPL + AC + BMENT + GAR+CITCOU+FIREPL_lag + AC_lag,data=baltimore,listw=W.listw)
durbin <- lagsarlm(PRICE~ NBATH + PATIO + FIREPL + AC + BMENT + GAR+CITCOU,data=baltimore,listw=W.listw,type="mixed")
semmod <- spautolm(PRICE~ NBATH + PATIO + FIREPL + AC + BMENT + GAR+CITCOU,data=baltimore,listw=W.listw)
semmod2 <- spautolm(PRICE~ NBATH + PATIO + FIREPL + AC + BMENT + GAR+CITCOU+FIREPL_lag+AC_lag,data=baltimore,listw=W.listw)

  ###############################
# Impacts
###############################

impacts(lagmodel, listw=W.listw)

  ###############################
# Model choice
###############################

lm.LMtests(ols_mod_2,listw=W.listw, test="all")
AIC(ols_mod_2,ols_croise_2,lagmodel,durbin,semmod,lagmodel2,semmod2)
#LMlag andLMerr are both significant hence we look at the robust tests
RLMlag is significant but RLMerr is not hence we select the LAG model
lm.LMtests(ols_croise_2,listw=W.listw, test="all")