
#{Lecture de bases de donn?es spatiales- slide 12 ? 16}
#{Les objets de type SpatialPolygonsDataFrame}

library(spdep)
columbus <- readShapePoly(system.file("etc/shapes/columbus.shp",package="spdep")[1])
col.gal.nb <- read.gal(system.file("etc/weights/columbus.gal",package="spdep")[1])
#reads a GAL lattice file into a neighbours list for spatial analysis.
class(columbus)
dim(columbus)
# 49 obs and 20 var
head(columbus@data)
str(columbus@data)
plot(columbus,axes=TRUE)
title("Neighbourhoods in Columbus")
CP<-as.numeric(as.factor(columbus@data$CP))
col.map<-c("royalblue2","red3")
plot(columbus,col=col.map[CP])
legend("topleft", legend = c("0","1"), cex = 0.8,title = "Core-periphery dummy ",fill=col.map[1:2])


#Construction d'un objet de type SpatialPointsDataFrame
library(GeoXp)
data(immob)
immob.sp = SpatialPoints(cbind(immob$longitude, immob$latitude))
class(immob.sp)
immob.spdf = SpatialPointsDataFrame(immob.sp, immob)
class(immob.spdf)
plot(immob.spdf)

#Construction d'un objet de type SpatialPixelsDataFrame
data(meuse.grid)
m = SpatialPixelsDataFrame(points = meuse.grid[c("x", "y")], data = meuse.grid)
class(m)
plot(m)


#Donn?es semis de points Pompiers de Toulouse
library(mgcv)
library(spatstat)
load("Pompiers_janvier+region.Rdata")
load("C:\\Users\\STAT-ECO\\Documents\\EUSTAT-VITORIA\\EUSTAT-SLIDES\\Pompiers_janvier+region.Rdata")
PP=ppp(sinistres_janvier$X,sinistres_janvier$Y,window=Region)
marks(PP)<-sinistres_janvier$M
plot(PP,main="Sinistres  avec dur?e ")
PPu=unmark(PP)
plot(PPu,main="Sinistres dans Region de Toulouse",cex=0.4)
print(PP)
summary(PP)
s=summary(PP)
str(s)



#donn?es columbus

library(classInt)
q5 <- classIntervals(columbus@data$INC , n=4, style="equal")
plot(columbus, col=findColours(q5, c("lightgreen", "darkgreen")))
q5b <- classIntervals(columbus@data$CRIME , n=4, style="equal")
plot(columbus, col=findColours(q5b, c("lightblue", "darkblue")))
q5c <- classIntervals(columbus@data$HOVAL , n=4, style="equal")
plot(columbus, col=findColours(q5c, c("tomato1", "tomato4")))

class(columbus)
columbus <- readShapePoly(system.file("etc/shapes/columbus.shp",
 package="spdep")[1])
class(columbus)

plot(columbus)
coord=coordinates(columbus)
plot(col.gal.nb,coord,add=TRUE)

# donn?es SIDS      (slide 37)
library(spdep)
nc_file <- system.file("etc/shapes/sids.shp", package = "spdep")
llCRS <- CRS("+proj=longlat +datum=NAD27")
nc <- readShapeSpatial(nc_file, ID = "FIPSNO", proj4string = llCRS)


#generation donn?es spatialement corr?l?es
X=rnorm(49)
W=nb2listw(col.gal.nb)
X1=lag.listw(W,X)
plot(X,X1)
Y=invIrW(W,0.8)%*%X
Y1=lag.listw(W,Y)
plot(Y,Y1)
moran.plot(X,W)

#Exercice: Donn?es Baltimore

library(spdep)
data(baltimore)
balt.sp = SpatialPoints(cbind(baltimore$X,baltimore$Y))
balt.spdf = SpatialPointsDataFrame(balt.sp, baltimore)
library(GeoXp)
driftmap(balt.spdf,, interpol=TRUE, nuage=TRUE)
variocloudmap(balt.spdf,"PRICE", quantiles=0.95)
library(geoR)
balt.geo=as.geodata(baltimore, coords.col = 16:17, data.col = 2)
vario.c <- variog(balt.geo,max.dist=90,op="smooth", direction=0.25, band=10)
plot(vario.c, main="variogram cloud")

