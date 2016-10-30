#présentation GeoXp  avec données immob

histomap(immob.spdf,"prix.vente", col="green",pch=5)
histomap(immob.spdf,"prix.vente", col="green",pch=5,identify="TRUE")
histomap(immob.spdf,"prix.vente", type="percent",xlab="Prix de vente")
histomap(immob.spdf,"prix.vente",criteria=(immob.spdf@data$prix.location>12))
densitymap(immob.spdf,"prix.vente", col="green",pch=5,xlab="Prix de vente")
boxplotmap(immob.spdf,"prix.vente")
scattermap(immob.spdf,c("rentabilite","prix.vente"),xlab="rentabilite",ylab="prix de vente")
barmap(columbus,"CP",col=c("orange","violet"),pch=c(2,4))

driftmap(columbus,"HOVAL",interpol=TRUE,nuage=TRUE)
angleplotmap(columbus,"Y",lablong="Longitude",lablat="Latitude",quantiles=0.95)

# avec structure de voisinage, données  columbus et SIDS

moranplotmap(columbus,"HOVAL",nb2listw(col.gal.nb),lablong="Longitude",lablat="Latitude")

# detection hot spots
localM <- localmoran(nc$SID74,nb2listw(knn2nb(wv.knn)),alternative="greater")
Ilocal=localM[,1]
pvalue=localM[,5]
moranplotmap(nc,"nc@data$SIDS74",nb2listw(knn2nb(wv.knn)))


resI <- localmoran(columbus$HOVAL, nb2listw(col.gal.nb))

