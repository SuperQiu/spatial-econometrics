#spatial autocorrelation
#areal data - continuous variable
# load columbus data

moran.test(columbus$HOVAL, nb2listw(col.gal.nb), randomisation=FALSE)
moran.test(columbus$HOVAL, nb2listw(col.gal.nb), randomisation=TRUE)
S=sample(columbus$HOVAL, length(columbus$HOVAL))
moran.test(S, nb2listw(col.gal.nb), randomisation=FALSE)

geary.test(columbus$HOVAL, nb2listw(col.gal.nb), randomisation=FALSE)
geary.test(columbus$HOVAL, nb2listw(col.gal.nb), randomisation=TRUE)

#moranplot 
moranplotmap(columbus,"HOVAL",nb2listw(col.gal.nb),lablong="Longitude",lablat="Latitude")


#areal data - binary variable

med=median(columbus$CRIME)
HICRIME=as.factor(columbus$CRIME>med)
levels(HICRIME)=c("faible","fort")

library(classInt)
q2=classIntervals(as.numeric(HICRIME) ,n=2)
plot(columbus, col=findColours(q2, c("lightgreen", "darkgreen")))
joincount.test(HICRIME,nb2listw(col.gal.nb))
joincount.mc(HICRIME,nb2listw(col.gal.nb),nsim=100)

#local Moran
 
resI <- localmoran(columbus$HOVAL, nb2listw(col.gal.nb))
# on stocke les résultats dans l'objet columbus
columbus@data<-data.frame(columbus@data,data.frame(resI))

require("GeoXp")
# on sélectionne sur  le deuxième graphique (graphique des p-values)
# les observations dont la p-value est significative (<0.05) et on appuie
# sur le bouton "save results" pour stocker les identifiants sélectionnés
highI<-dbledensitymap(columbus, c("Ii","Pr.z...0."))
# on affiche les valeurs du I de moran local sur ces observations
dev.set(2)
text(coordinates(columbus[last.select,]),as.character(round(columbus@data$Ii[last.select],3)))



localM <- localmoran(nc$SID74,nb2listw(knn2nb(wv.knn)),alternative="greater")
Ilocal=localM[,1]
pvalue=localM[,5]
moranplotmap(nc,"nc@data$SIDS74",nb2listw(knn2nb(wv.knn)))


# spatial autocorrelation of residuals
lmmod=lm(HOVAL~CRIME+INC,data=columbus)
lm.morantest(lmmod,nb2listw(col.gal.nb))

#spatial autocorrelation of geostatistical data
#load baltimore data
library(GeoXp)
variocloudmap(balt.spdf,"PRICE", quantiles=0.95)
library(geoR)
balt.geo=as.geodata(baltimore, coords.col = 16:17, data.col = 2)
vario.c <- variog(balt.geo,max.dist=90,op="smooth", direction=0.25, band=10)
plot(vario.c, main="variogram cloud")


#tests about PP
window=owin(c(0,10),c(0,10))
poisson=rpoispp(10,win=window)
quadrat.test(poisson)
tt=quadrat.test(poisson)
plot(tt)
tt=quadrat.test(poisson,4)
plot(tt)
poisson_inhom=rpoispp(function(x,y){100*exp(-3*x)+100*exp(-3*y)},20,win=window)
quadrat.test(poisson_inhom)
plot(envelope(cells,Fest))
plot(envelope(cells,Gest))
plot(envelope(PP,Kest))
plot(envelope(PP,Kinhom))
fitPP=ppm(PPu,~z, covariates=list(z=Z))
quadrat.test(fitPP)
plot(envelope(fitPP,Kest))