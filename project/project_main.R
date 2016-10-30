
#load packages 
library(spdep)
library(sp)
library(maptools)
library(GeoXp)
library(glmnet)

***************** READ DATA **********************

#map data

district<- readShapeSpatial("C:/Users/USER/Desktop/TSE M2 Eco STat/spatial/map_district_region.shp", 
                               ID = "CODGEO") 

class(district)
dim(district)
attach(district@data)
head(district@data)
str(district@data) #this would give the structure of the data
plot(district,axes=TRUE)
coord=coordinates(district) #coordinates of center obtained for each district
class(coord)

#explainatory data

Xfile<-read.csv("C:/Users/USER/Desktop/TSE M2 Eco STat/spatial/Xfile.csv", header=TRUE, sep=",", dec=".")
class(Xfile)
class(Xfile$labour_force)
dim(Xfile)
attach(Xfile)

#variable of interest = flow

flux<-read.table("C:/Users/USER/Desktop/TSE M2 Eco STat/spatial/flux_district.txt")
class(flux)
dim(flux)
attach(flux)
View(flux)
y<-flux[,3] # rename the flux i.e V3 as y

#Here we observe that the dim of flux is 3059*3. This means that a few flows are not included as there are no flows in those o-d's predent. These are then taken to be 0.
#corrected for 3600

y<-rep(0,3600)
for(i in 1:3059){ y[(V1[i]-1)*60+V2[i]]<-flux[i,3]}
#View(y)
class(y)

d<-rep(1:60,60)
#View(d)
o<-rep(1:60, each=60)
#View(o)

#newflux

newflux<-cbind(o,d,y)
View(newflux)
class(newflux[,2])
newflux1<-as.data.frame(newflux)

*********EXPLORATORY PHASE***********

##1

# calculate the outgoing flow
out_f<-rep(0,60)
for(i in 1:60){out_f[i]<-sum(flux[V1==i,3])}
#View(out_f)
out_f[1]

# calculate the incoming flow
in_f<-rep(0,60)
for(i in 1:60){in_f[i]<-sum(flux[V2==i,3])}

summary(out_f)
table(in_f)

# flows converted to categorical variables
out=rep(0,60)
for (i in 1:60)
{ if (out_f[i]<1000)
{  out[i]=1 }
if (out_f[i]>1000 & out_f[i]<2000)
{out[i]=2}
if (out_f[i]>2000 & out_f[i]<3000)
{out[i]=3}
if (out_f[i]>3000)
{out[i]=4}
}

inf=rep(0,60)
for (i in 1:60)
{ if (in_f[i]<1000)
{  inf[i]=1 }
if (in_f[i]>1000 & out_f[i]<2000)
{inf[i]=2}
if (in_f[i]>2000 & out_f[i]<3000)
{inf[i]=3}
if (in_f[i]>3000)
{inf[i]=4}
}

# Compute and map the flows
par(mfrow=c(1:2))

col.map<-c("bisque","bisque2","bisque3","bisque4")
plot(district,col=col.map[out])
title("Outgoing flows")

col.map<-c("bisque","bisque2","bisque3","bisque4")
plot(district,col=col.map[inf])
title("Incoming flows")

legend("bottomright", legend=c("<1000","1000-2000","2000-3000",">3000"), cex=.9,fill=col.map[1:4])
flo<-cbind(Xfile,inf,out)
View(flo)

##2

#nbd matrix computation and analyses 

par(mfrow=c(1,1))
class(coord)

#wt matrix based on common border
wc.nb=poly2nb(district) ####### weight matrix in nb ######
class(wc.nb) # class nb
is.symmetric.nb(wc.nb)
str(wc.nb)

wmat<-nb2mat(wc.nb) ####### weight matrix in matrix format ######
View(wmat)
listw=nb2listw(wc.nb)

#here we get the map and the coordinates of the neighbours joined
plot(district, border='grey')
plot(wc.nb,coord,add=TRUE) 
summary(wc.nb)

#analysis by neighbourmap and barnmap
head(district@data)

#here selection can be make by point or polygon 
#neighbourmap graph has latitudes (here..see code below) on both the axis 
#when "selection is made by point", you can select a pt on the map and the nbs of this point are selected on the graph

neighbourmap(district,"LATITUDE",wc.nb) # graph and map obtained

#bar map has a bar with the no of nbs on X axis and with no of districts with those many nbs on the y axis
#when "selection is made by point", you can select a pt on the map, you can see the nbs of this pt. This is also seen on the bar.

barnbmap(district,wc.nb) # bar and map obtained


##3
#spatial exploratory analysis of the incoming and outgoing flows

district@data<-cbind(district@data,out_f,in_f)
district@data<-cbind(district@data,outc,infc)

par(mfrow=c(1,1))
histomap(district,"out_f",col="red",pch=5)
densitymap(district,"out_f", col="green",pch=5,xlab="outflows")
boxplotmap(district,"out_f")
scattermap(district,c("out_f","in_f"),xlab="outflow",ylab="inflow")

histomap(district,"in_f",col="red",pch=5)
densitymap(district,"in_f", col="green",pch=5,xlab="outflows")
boxplotmap(district,"in_f")
scattermap(district,c("in_f","in_f"),xlab="outflow",ylab="inflow")

#center out_f and in_f. 
#H0:no spa autocorr
outc<-as.vector(scale(out_f,center=TRUE))
infc<-as.vector(scale(in_f,center=TRUE))

#moranplotmap for flows
moranplotmap(district,"outc",nb2listw(wc.nb),lablong="Longitude",lablat="Latitude")
moranplotmap(district,"infc",nb2listw(wc.nb),lablong="Longitude",lablat="Latitude")


#################################################################################################

### MODELLING PHASE ###

#Steps
# DISTANCES ARE CALCULATED
# Xo and Xd are found
# regression is done

#lat long origin
X0dist<-cbind(Xfile,district)

View(X0dist)
names(X0dist)
dim(X0dist)
X0dist <- X0dist[c(-3,-4,-5,-6,-7,-8,-9,-12)]

Xodis<-matrix(0,3600,4)
for (j in 1:4){for(i in 1:3059){Xodis[(V1[i]-1)*60+V2[i],j]<-X0dist[which(X0dist$ID_district==V1[i]),j]}}
View(Xodis)
colnames(Xodis)<-c("ori","id_o","lato","longo")

# LAT LONG dest
Xddist<-cbind(Xfile,district)
View(Xddist)
names(Xddist)
Xddist <- Xddist[c(-3,-4,-5,-6,-7,-8,-9,-12)]
View(Xddist)
Xddis<-matrix(0,3600,4)
for(j in 1:4){for(i in 1:3059){Xddis[(V1[i]-1)*60+V2[i],j]<-Xddist[which(Xddist$ID_district==V2[i]),j]}}
View(Xddis)
colnames(Xddis)<-c("dest","id_d","latd","longd")


#calculate the distances

fordist<-cbind(Xodis,Xddis)
View(fordist)
class(fordist)
dim(fordist)

lato<-as.numeric(fordist[,3])
#class(lato)
longo<-as.numeric(fordist[,4])
latd<-as.numeric(fordist[,7])
longd<-as.numeric(fordist[,8])

dist<-function(X1,X2,X3,X4) sqrt(((X1 - X3) ^ 2)+((X2 - X4) ^ 2) )

for (i in 1:3600) dist[i]
{ euc.dist<-dist(lato,longo,latd,longd)
  euc.dist
}

View(euc.dist)### all the distances####
Euc.dist<-cbind(fordist,euc.dist)
View(Euc.dist)
class(Euc.dist)
Euc.dist<-Euc.dist[,-1]
distance<-Euc.dist[,8]
View(distance)

 
#get the matrix Xo and Xd

Xo<-matrix(0,3600,7)
for (j in 1:7){for(i in 1:3059){Xo[(V1[i]-1)*60+V2[i],j]<-Xfile[which(Xfile$ID_district==V1[i]),j]}}
class(Xo)
View(Xo)                               
Xo<-cbind(Xo,newflux)
Xo<-Xo[,-1]
Xori<-as.data.frame(Xo)
names(Xori)<-c("id_o","labo","acro","empo","unempo","huo","origin","destination","flow")
View(Xori)
class(Xori$acro)
Xo<-as.data.frame(Xori)
View(Xo)


Xd<-matrix(0,3600,7)
for (j in 1:7){for(i in 1:3059){Xd[(V1[i]-1)*60+V2[i],j]<-Xfile[which(Xfile$ID_district==V2[i]),j]}}
#View(Xd) 
Xd<-cbind(Xd,newflux)
Xd<-Xd[,-1]
Xdes<-as.data.frame(Xd)
names(Xdes)<-c("id_d","labd","acrd","empd","unempd","hud","origin","destination","flow")
#View(Xdes)
class(Xdes$acrd)
Xd<-as.data.frame(Xdes)                      
#View(Xd)


#intra and inter

Xointer<-Xo
for (i in 1:3600) {for (j in 1:10) {if( Xo$origin[i]==  Xo$destination[i]){ Xointer[i, 2:6]=rep(0,5)}}}
View(Xointer)

Xointra<-Xo
for (i in 1:3600) {for (j in 1:10) {if( Xo$origin[i]!=  Xo$destination[i]){ Xointra[i, 2:6]=rep(0,5)}}}
View(Xointra)


Xdinter<-Xd
for (i in 1:3600) {for (j in 1:10) {if( Xd$origin[i]==  Xd$destination[i]){ Xdinter[i, 2:6]=rep(0,5)}}}
View(Xdinter)

Xdintra<-Xd
for (i in 1:3600) {for (j in 1:10) {if( Xd$origin[i]!=  Xd$destination[i]){ Xdintra[i, 2:6]=rep(0,5)}}}
View(Xdintra)

###Xdintra and X0intra are obviously the same

#final data is obtained name : data 

data<-cbind(y,distance,Xo[,2:6],Xd[,2:6],Xointra[,2:6],Xdintra[,2:6])
View(data)

colnames(data)<-c("y","dist","lfo","acro","empo","unempo","ho","lfd","acrd","empd","unempd","hd",
                  "l_o_intra","acr_o_intra","em_o_intra","unem_o_intra","hou_o_intra",
                  "l_d_intra","acr_d_intra","em_d_intra","unem_d_intra","hou_d_intra")
View(data)               

#1
#Modeling phase

#STEP 1 : OLS model without adjustment. Variables selected by backward


reg1<-lm(y~dist+acro+ho+lfd+acrd+unempd+hd +empd            
                  , data=data)
summary(reg1)
AIC(reg1)

cov(data)
cor(data$lfo,data$acro)

#STEP 2 :OLS with adjustments   
##traditional gravity model

reg1_adj<-lm(y~dist+acro+lfd+acrd+empd+unempd+hd+
l_o_intra+acr_o_intra+em_o_intra+unem_o_intra+hou_o_intra            
                  , data=data)
summary(reg1_adj)
AIC(reg1_adj)
res<-residuals(reg1_adj)

#2

#Wd and Wo created

library(Matrix)
I<-diag(60)
Wd<-kronecker(I,wmat)
Wd[61:120,61:120]

Wo<-kronecker(wmat,I)
Wo[1:60,61:120]

Wt<-Wd+Wo
dim(Wt)
wt<-mat2listw(Wt)
class(wt)

w<-Wt%*%y
w

#analysis of residuals 

moran.test(res,mat2listw(Wt),randomisation=FALSE) # pval <.05 so spatial autocorr hencedo a lag model


#3
#STEP 3 : Model with lag of dep var

lagreg<-lagsarlm(y~dist+acro+lfd+acrd+empd+unempd+hd+
l_o_intra+acr_o_intra+em_o_intra+unem_o_intra+hou_o_intra ,           
             listw= wt    , data=data,tol.solve=1e-20)
summary(lagreg)
AIC(lagreg)

#4
#STEP 4 : Durbin model

l_o_intra_lag<-lag.listw(wt,data$l_o_intra)
acr_o_intra_lag<-lag.listw(wt,data$acr_o_intra)
em_o_intra_lag<-lag.listw(wt,data$em_o_intra)
hou_o_intra_lag<-lag.listw(wt,data$hou_o_intra)
unem_o_intra_lag<-lag.listw(wt,data$unem_o_intra)


acro_lag<-lag.listw(wt,data$acro)
empd_lag<-lag.listw(wt,data$empd)
unempd_lag<-lag.listw(wt,data$unempd)
hd_lag<-lag.listw(wt,data$hd)
acrd_lag<-lag.listw(wt,data$acrd)
lfd_lag<-lag.listw(wt,data$lfd)


reg_durbin_adj<-lagsarlm(y~dist+acro+lfd+acrd+empd+unempd+hd+
l_o_intra+acr_o_intra+em_o_intra+unem_o_intra+hou_o_intra
+em_o_intra_lag+l_o_intra_lag+acr_o_intra_lag+unem_o_intra_lag+hou_o_intra_lag+
acro_lag+unempd_lag+hd_lag+acrd_lag
+empd_lag +lfd_lag ,listw= wt  , data=data,tol.solve=1e-20)
summary( reg_durbin_adj)
AIC(reg_durbin_adj)


#PREDICTION


#1
#First we eliminate 10 % of the flows and work with 3240 observations. We fit the lag model for 3240 observations. 
#We then use the coeffs we obt here to find the predictions for 360 eliminated obs.
#In the formula, Yod hat would give the fitted values for 3240 obs. Y1 are the 3240 obs and the Xod and Wod are the explanatory variables and weight matrix for the corresponding 3240 obs.

library(forecast)

v<-sample(1:3600,360) # v has 3240 observations

#training data (3240 obs)
data_train<-data[-v,] 
Wt_train<-Wt[-v,-v]
wt_train<-mat2listw(Wt_train)

#testing data (360 obs)
data_test<-data[v,]
Wt_test<-Wt[v,v]

#2
##trad gravity prediction for non eliminated flows

#model fitted on training data
reg1_adj_tr<-lm(y~dist+acro+lfd+acrd+empd+unempd+hd+
               l_o_intra+acr_o_intra+em_o_intra+unem_o_intra+hou_o_intra            
             , data=data_train)
summary(reg1_adj_tr)
#coefficients obtained
coef<-as.vector(coefficients(reg1_adj_tr)) # coefs for 3240 obs

#now that the coeffs are obtained, we make the predictions on testing data(360 obs)
#explanatory variables for testing data

X_trad<-cbind(rep(1,360),data_test$dist,data_test$acro, data_test$lfd, data_test$acrd, data_test$empd, data_test$unempd, data_test$hd,
             data_test$l_o_intra, data_test$acr_o_intra, data_test$em_o_intra, data_test$unem_o_intra, data_test$hou_o_intra)

# predictions
pred_test<-X_trad%*%coef

#quadratic mean error   
error_trad<-mean((pred_test-y[v])^2)

#relative quadratic mean error
Y<-y[v]
y_new<-Y[which(Y!=0)]
View(y_new)
pred_test<-pred_test[which(Y!=0)]
rerror_test<-sum(((pred_test-y_new)/y_new)^2)


#3
##lag model prediction for non eliminated flows

# model fitted on training data

lagreg_p<-lagsarlm(y~dist+acro+lfd+acrd+empd+unempd+hd+
l_o_intra+acr_o_intra+em_o_intra+unem_o_intra+hou_o_intra ,           
             listw= wt_train   , data=data_train,tol.solve=1e-20)
summary(lagreg_p)
AIC(lagreg_p)
coef_lag<-as.vector(coefficients(lagreg_p))# coefs for training data

# weights and explanatory variables for test data

wy<-Wt_test%*%y[v]#360 obs

X_lag<-cbind(rep(1,360),data_test$dist,data_test$acro, data_test$lfd, data_test$acrd, data_test$empd, data_test$unempd, data_test$hd,
             data_test$l_o_intra, data_test$acr_o_intra, data_test$em_o_intra, data_test$unem_o_intra, data_test$hou_o_intra)

rhop<-diag(rep(1,360))-coef_lag[1]*Wt_test
rhopinv<-solve(rhop)*rhop

length(coef_lag)
dim(X_lag)
#predictions for test data
pred_lag<-rhopinv%*%(X_lag%*%coef_lag[-1])

#quadratic mean error 
error_lag<-mean((pred_lag-y[v])^2)

#relative quadratic mean error
Y<-y[v]
y_new<-Y[which(Y!=0)]
View(y_new)
pred_l<-pred_lag[which(Y!=0)]
rerror_lag<-sum(((pred_l-y_new)/y_new)^2)

#Relative error in the lag model is much less than in trad gravity model.

# Cross validation

CV(reg1_adj_tr)







































































