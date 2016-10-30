#generation donn�es spatialement corr�l�es
#charger la structure de voisinage des donn�es columbus
X=rnorm(49)
W=nb2listw(col.gal.nb)
X1=lag.listw(W,X)
plot(X,X1)
Y=invIrW(W,0.8)%*%X
Y1=lag.listw(W,Y)
plot(Y,Y1)
moran.plot(X,W)


#matrices de voisinage
#avec les donn�es SIDS

coordnc=coordinates(nc)

#matrice bas�e sur un seuil de distance
wd.nb=dnearneigh(coordnc,0,75,longlat=TRUE)
class(wd.nb)
plot(nc, border='grey',xlim=c(-84.5,-82),ylim=c(35,36),
axes=TRUE)
coord=coordinates(nc)
plot(wd.nb,coord,add=TRUE)
wd.listw=nb2listw(wd.nb)

#matrice bas�e sur un nombre de plus proches voisins
wv.knn=knearneigh(coordnc, k=4, longlat = TRUE)
class(wv.knn)
str(wv.knn)
plot(nc, border='grey',xlim=c(-84.5,-82),ylim=c(35,36),
axes=TRUE)
plot(knn2nb(wv.knn), coord, add=TRUE)


#matrice bas�e sur la contiguit�
wc.nb=poly2nb(nc)
class(wc.nb)
is.symmetric.nb(wc.nb)
str(wc.nb)
plot(nc, border='grey',xlim=c(-84.5,-82),ylim=c(35,36),
axes=TRUE)
coord=coordinates(nc)
plot(wc.nb,coord,add=TRUE)


#matrice de Delaunay
w.tri=tri2nb(coordnc)
class(w.tri)
plot(nc, border='grey',axes=TRUE)
plot(w.tri, coord, add=TRUE)

#variable spatialement d�cal�e 
nc$SID74_lag.B=lag.listw(nb2listw(knn2nb(wv.knn), style="B"),
nc$SID74)
nc$SID74_lag.W=lag.listw(nb2listw(knn2nb(wv.knn), style="W"),
nc$SID74)

#analyse exploratoire d'une matrice de voisinage
summary(w.tri)
library(GeoXp)
neighbourmap(nc, "east", wd.nb)
barnbmap(nc,wd.nb)
histnbmap(nc,knn2nb(wv.knn))




#la classe listw, petit exemple
t=c(1,2,3,4)
u=c(3,2,5,1)
plot(t,u)
co=cbind(t,u)
W.knn=knearneigh(co,k=2,longlat=TRUE)
plot(knn2nb(W.knn), co, add=TRUE)
W.nb=knn2nb(W.knn)
W.listw1=nb2listw(W.nb,style="B")
str(W.listw1)
W.listw1$neighbours[]
W.listw1$weights[]
W.listw2=nb2listw(W.nb,style="W")
W.listw2$neighbours[]
W.listw2$weights[]