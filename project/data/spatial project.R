##########import data
##############################
##### map data
district_poly<- readShapePoly("/Users/luozhaoqiu/Desktop/M2S2/spatial /project/data files for project/map_district_region.shp")
class(district)
dim(district)
attach(district@data)
head(district@data)
str(district@data)
plot(district,axes=TRUE)
coord=coordinates(district)

district_sp<- readShapeSpatial("/Users/luozhaoqiu/Desktop/M2S2/spatial /project/data files for project/map_district_region.shp", 
                                ID = "CODGEO") 
##### explainatory data
Xfile<-read.csv(file=file.choose(), header=TRUE, sep=";", dec=".")
class(Xfile)
dim(Xfile)
attach(Xfile)

##### variable of interest = flow
read.table("/Users/luozhaoqiu/Desktop/M2S2/spatial /project/data files/flux_district.txt")
flux<-read.table(file=file.choose())
class(flux)
dim(flux)
attach(flux)
y<-flux[,3]
#########################################

##########part one
###################################

##### calculate the outgoing flow
out_f<-rep(0,60)
for(i in 1:60){out_f[i]<-sum(flux[V1==i,3])}

##### calculate the incoming flow
in_f<-rep(0,60)
for(i in 1:60){in_f[i]<-sum(flux[V2==i,3])}

##### mapping the flows

out<-rep(0,60)
for (i in 1:60)
{
  if (out_f[i]<1000 )
  {
    out[i]=1
  }
  if (out_f[i]>1000 & out_f[i]<2000)
  {
    out[i]=2
  }
  if (out_f[i]>2000 &  out_f[i]<3000)
  {
    out[i]=3
  }else
  {
    out[i]=4
  }
}

col.map<-c("aquamarine","aquamarine2","aquamarine3","aquamarine4")
plot(district_sp,col=col.map[out])
##########################################


##### Xo and Xd
#name_district<-as.character(name_district)

X<-matrix(0,3059,7)
for(j in 1:7){for(i in 1:3059){X[i,j]<-Xfile[which(Xfile$ID_district==V1[i]),j]}}
Xo<-cbind(X,flux)
colnames(Xo)<-c(names(Xfile),"origin","destination","flow")


X<-matrix(0,3059,7)
for(j in 1:7){for(i in 1:3059){X[i,j]<-Xfile[which(Xfile$ID_district==V2[i]),j]}}
Xd<-cbind(X,flux)
names(Xd)<-c(names(Xfile),"origin","destination","flow")


##############

### MODELLING PHASE ###
# DISTANCES ARE CALCULATED
# Xo and Xd are found
# regression is done

#lat long origin
X0dist<-cbind(Xfile,district)

View(X0dist)
names(X0dist)
X0dist <- X0dist[c(-3,-4,-5,-6,-7,-8,-9,-12,-13,-14,-15,-16)]

name_district<-as.character(name_district)

Xodis<-matrix(0,3059,4)
for(j in 1:4){for(i in 1:3059){Xodis[i,j]<-X0dist[which(X0dist$ID_district==V1[i]),j]}}

View(Xodis)
#colnames(Xodis)<-c("origin","id_o","lato","longo")

######################################################################################
# LAT LONG dest

Xddist<-cbind(Xfile,district)
View(Xddist)
names(Xddist)
name_district<-as.character(name_district)
Xddist <- Xddist[c(-3,-4,-5,-6,-7,-8,-9,-12,-13,-14,-15,-16)]
View(Xddist)
Xddis<-matrix(0,3059,4)
for(j in 2:4){for(i in 1:3059){Xddis[i,j]<-Xddist[which(Xddist$ID_district==V2[i]),j]
                               Xddis[i,1]<-as.character(Xddist[which(Xddist$ID_district==V2[i]),1])
                               cbind(Xddis[i,1],Xddis[i,j])}}

View(Xddis)
#colnames(Xddis)<-c("dest","id_d","latd","longd")

###################################################################################
##calculate the distances

fordist<-cbind(Xodis,Xddis)
View(fordist)
class(fordist)
fordist<-data.matrix(fordist)
lato<-as.numeric(fordist[,3])
class(lato)
longo<-as.numeric(fordist[,4])
latd<-as.numeric(fordist[,7])
longd<-as.numeric(fordist[,8])

dist<-function(X1,X2,X3,X4) sqrt(((X1 - X3) ^ 2)+((X2 - X4) ^ 2) )

for (i in 1:3059) dist[i]
{ euc.dist<-dist(lato,longo,latd,longd)
  euc.dist
}

View(euc.dist)### all the distances####
Euc.dist<-cbind(fordist,euc.dist)
View(Euc.dist)

colnames(Euc.dist)<-c("origin","id_o","lato","longo","dest","id_d","latd","longd","distance")



##################################
##### nb of employment working inside their region
inside<-employment-out_f
inside

##### nb of jobs
nb_job<-inside+in_f

X<-cbind(Xfile,in_f,out_f,inside,nb_job)
X[,c(1,3,5,8,9,10)]
X[,c(1,3,11)]

##### calculate the coverage rate
coverage_rate<-nb_job/labour_force
coverage_rate

#################################
Xointer<-Xo
for(j in 1:10){for(i in 1:3059){if (Xo$origin[i]==Xo$destination[i]){Xointer[i,3:7]=rep(0,5)}}}
Xointer<-Xointer[,3:7]
Xointra<-Xo
for(j in 1:10){for(i in 1:3059){if (Xo$origin[i]!=Xo$destination[i]){Xointra[i,3:7]=rep(0,5)}}}
Xointra[,3:7]
Xdinter<-Xd
for(j in 1:10){for(i in 1:3059){if (Xd$origin[i]==Xo$destination[i]){Xdinter[i,3:7]=rep(0,5)}}}
Xdinter[,3:7]
Xdintra<-Xd
for(j in 1:10){for(i in 1:3059){if (Xd$origin[i]!=Xd$destination[i]){Xdintra[i,3:7]=rep(0,5)}}}
Xdintra[,3:7]



data<-cbind(y,Xointer,Xointra,Xdinter,Xdintra,distance)



v<-sample(1:3600,360)
data[-v,]
Wt<-Wt[-v,-v]



