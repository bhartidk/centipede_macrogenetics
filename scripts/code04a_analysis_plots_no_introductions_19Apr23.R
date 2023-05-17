#### load required libraries and set working directory ####

#load libraries
library(dplyr)
library(glmmTMB)
library(bbmle)
library(betareg)
library(corrplot)
library(geodist)
library(ape)
library(spdep)
library(ncf)
library(adespatial)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(rgdal)
library(adegraphics)
library(memisc)
library(ggrepel)
library(boot)
library(effects)
library(phytools)
library(broom)
library(gtools)
library(generics)
library(gt)
library(gtsummary)
library(CoordinateCleaner)

#setting the working directory
setwd("D/PostDoc data/centipede_popgen/github_11Apr23_version3_final")

#### exploratory graphs ####

#input sequence summary
seq.sum<-read.csv(file="data/input_betareg_no_introductions_11Apr23.csv", header=TRUE, na.strings=c("","NA"))

#number of species with complete information
nrow(seq.sum[complete.cases(seq.sum),]) #128

#create a column with latitudinal range
seq.sum$sp.lat.range<-seq.sum$int.max.lat-seq.sum$int.min.lat

#taking a subset of seq.sum which have non-NA data
seq.sum.sub<-seq.sum[,c("species", "family", "order", "avg.pair.diff", "lat.mean", "lon.mean", "avg.geo.dist", "sp.lat.range", "body.size.new", "blindness", "maternal.care", "n", "n.unique.sites", "seq.len", "hapl.div")]

#look at the distribution of genetic diversity values by order
pdf(file=paste0("figures/histogram_genetic_diversity_11Apr23.pdf"), height=11, width=8, useDingbats=FALSE)
ggplot(seq.sum.sub, aes(x=avg.pair.diff, fill=order, color=order)) + geom_histogram(alpha=0.5, position="identity") +  facet_grid(order ~ .) + xlab("Average pairwise difference") + ylab("Frequency") + theme_bw()
dev.off()

#look at the relationship between body size and genetic diversity
pdf(file=paste0("figures/genetic_diversity_body_size_11Apr23.pdf"), height=8, width=11, useDingbats=FALSE)
ggplot(seq.sum.sub, aes(y=avg.pair.diff, x=body.size.new)) + geom_point() + geom_text_repel(label=rownames(seq.sum.sub)) + ylab("Average pairwise difference") + xlab("Body size (mm)") + theme_bw()
dev.off()
#look at the data for species which have large body size but low genetic diversity and may be driving the relationship
seq.sum.sub[c(108,109,110,111,113,116,118),] #these are all Scolopendra species (cingulata, dawydoffi, dehaani, hainanum, mojiangica, mutilans, paradoxa) with low genetic diversity values

#look at the relationship between maternal care and genetic diversity
pdf(file=paste0("figures/genetic_diversity_maternal_care_11Apr23.pdf"), height=8, width=11, useDingbats=FALSE)
ggplot(seq.sum.sub, aes(x=maternal.care, y=avg.pair.diff)) + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_bw() + xlab("Maternal care") + ylab("Average pairwise difference")
dev.off()

#look at the relationship between blindness and genetic diversity
pdf(file=paste0("figures/genetic_diversity_blindness_11Apr23.pdf"), height=8, width=11, useDingbats=FALSE)
ggplot(seq.sum.sub, aes(x=blindness, y=avg.pair.diff)) + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_bw() + xlab("Blindness") + ylab("Average pairwise difference")
dev.off()

#look at the relationship between mean latitude and genetic diversity
pdf(file=paste0("figures/genetic_diversity_latitude_11Apr23.pdf"), height=8, width=11, useDingbats=FALSE)
ggplot(seq.sum.sub, aes(y=avg.pair.diff, x=lat.mean)) + geom_point() + geom_text_repel(label=rownames(seq.sum.sub)) + ylab("Average pairwise difference") + xlab("Mean latitude") + theme_bw()
dev.off()
#look at the data for species which are found in the southern hemisphere and have high values of genetic diversity 
seq.sum.sub[c(1,2,26,80,125),] #Craterostigmus species, Henicops maculatus, Cryptops pictus (suspected to be a species complex), Allothereua serrulata
#look at the data for species which are found in equatorial regions and have high values of genetic diversity 
seq.sum.sub[c(77,93,99,114,119),] #Cryptops japonicus, Otostigmus astenus, O sulcipes, Scolopendra morsitans, S pinguis
#look at the data for species which are found in the northern hemisphere and have high values of genetic diversity 
seq.sum.sub[c(16,17,19,29,34,40,48,56,58,107),] #Stenotaenia sorrentina, Strigamia acuminata, S crassipes, Lithobius aeruginosus, L castaneus, L forficatus, L muticus, L tricuspis, L variegatus, Scolopendra canidens

#look at the relationship between average geographic distance and genetic diversity
pdf(file=paste0("figures/genetic_diversity_geographic_distance_11Apr23.pdf"), height=8, width=11, useDingbats=FALSE)
ggplot(seq.sum.sub, aes(y=avg.pair.diff, x=avg.geo.dist)) + geom_point() + geom_text_repel(label=rownames(seq.sum.sub)) + ylab("Average pairwise difference") + xlab("Average geographic distance") + theme_bw()
dev.off()
#look at the data for species where sequences are separated by >1000 km distance
seq.sum.sub[c(18,26,45,93,94,97,100,102,107,112,114,121,123),] #Strigamia chionophila, Henicops maculatus, Lithobius melanops, Otostigmus astenus, Otostigmus multidens, O scaber, Rhysida immarginata, R longipes, Scolopendra canidens, S japonica, S morsitans, Scolopocryptops rubiginosus, Theatops erythrocephalus

#look at the relationship between species latitudinal range and genetic diversity
pdf(file=paste0("figures/genetic_diversity_species_range_11Apr23.pdf"), height=8, width=11, useDingbats=FALSE)
ggplot(seq.sum.sub, aes(y=avg.pair.diff, x=sp.lat.range)) + geom_point() + geom_text_repel(label=rownames(seq.sum.sub)) + ylab("Average pairwise difference") + xlab("Species latitudinal range") + theme_bw()
dev.off()

#### prepare input files for analysis ####

#check number of species where average pairwise difference equals 0
seq.sum.sub %>%
filter(avg.pair.diff==0) %>%
pull(species)

#we have 9 species where the average pairwise difference value is equal to 0, which cannot be modeled by beta regression - using a modification function to deal with this - code from Douma et al., 2019
transform01<-function(x){
((x*(length(x)-1))+0.5)/(length(x))
}

#creating a new column with the transformed genetic diversity values
seq.sum.sub$avg.pair.diff.tf<-transform01(seq.sum.sub$avg.pair.diff)

#standardize the predictor variables so that they have a mean of 0 and SD of 1
#I would like to compare the relative influence of each of the predictors on the response variable. So I need to scale the predictors in such a way that they range between 0 and 1. I can then interpret the change in the response in terms of standard deviations. For example, 1 SD change in my predictor gives rise to X SD change in the response variable

#I can scale all the continuous predictor variables so that they have a mean of 0 and a standard deviation of 1
seq.sum.sub$body.size.new<-as.numeric(as.character(seq.sum.sub$body.size.new))
new.cols<-lapply(seq.sum.sub[,c("lat.mean", "lon.mean", "avg.geo.dist", "sp.lat.range", "body.size.new", "n")], scale)
seq.sum.sub2<-cbind(seq.sum.sub[,c("species", "family", "avg.pair.diff.tf", "blindness", "maternal.care")], new.cols)
colnames(seq.sum.sub2)

#checking the correlation between various predictor variables
#looking at the correlation coefficients - code from http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
colnames(seq.sum.sub)

res.pear<-cor(seq.sum.sub[,c("lat.mean", "avg.geo.dist", "sp.lat.range", "body.size.new", "n")], method="pearson")
pdf(file=paste0("figures/pearson_correlation_predictors_11Apr23.pdf"), height=8, width=11, useDingbats=FALSE)
corrplot(res.pear, type="upper", order="hclust", tl.col="black", tl.srt=45, addCoef.col='black')
dev.off()

#look at association between blindness and maternal care
table(seq.sum.sub$blindness, seq.sum.sub$maternal.care) #all species that are blind also show maternal care

#look at the association between body size and maternal care
pdf(file=paste0("figures/body_size_maternal_care_11Apr23.pdf"), height=8, width=11, useDingbats=FALSE)
ggplot(seq.sum.sub, aes(x=maternal.care, y=body.size.new)) + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_bw() + xlab("Maternal care") + ylab("Body size (mm)") #species with maternal care are generally of larger body size
dev.off()

#look at the association between body size and blindness
pdf(file=paste0("figures/body_size_blindness_11Apr23.pdf"), height=8, width=11, useDingbats=FALSE)
ggplot(seq.sum.sub, aes(x=blindness, y=body.size.new)) + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_bw() + xlab("Blindness") + ylab("Body size (mm)") #species that have vision (are not blind) are generally of larger body size
dev.off()

#### run beta regression models ####

#running different beta regression models using the glmmTMB function
#all the variables included
mod1<-glmmTMB(avg.pair.diff.tf ~ lat.mean + avg.geo.dist + sp.lat.range + body.size.new + blindness + maternal.care + (1|family), dispformula=~n, data=seq.sum.sub2, list(family="beta", link="logit"))
summary(mod1)

#removing the random effect
mod2<-glmmTMB(avg.pair.diff.tf ~ lat.mean + avg.geo.dist + sp.lat.range + body.size.new + blindness + maternal.care, dispformula=~n, data=seq.sum.sub2, list(family="beta", link="logit"))
summary(mod2)

#removing the variable dispersion parameter
mod3<-glmmTMB(avg.pair.diff.tf ~ lat.mean + avg.geo.dist + sp.lat.range + body.size.new + blindness + maternal.care + (1|family), data=seq.sum.sub2, list(family="beta", link="logit"))
summary(mod3)

#removing both the random effect and variable dispersion factor
mod4<-glmmTMB(avg.pair.diff.tf ~ lat.mean + avg.geo.dist + sp.lat.range + body.size.new + blindness + maternal.care, data=seq.sum.sub2, list(family="beta", link="logit"))
summary(mod4)

#compare all the different kinds of models
AICtab(mod1, mod2, mod3, mod4)
#model2 has the lowest AIC value

#since I am not getting r2 values from glmmTMB and since anyway the model with random effect does not turn out to be important - I am using the betareg package for model evaluation and analyzing coefficients. I checked and confirmed that the coefficients across these packages are the same.
mod2.logit<-betareg(avg.pair.diff.tf ~ lat.mean + avg.geo.dist + sp.lat.range + body.size.new + blindness + maternal.care | n, data=seq.sum.sub2)
summary(mod2.logit)

#checking for multi-collinearity - "Any variable with a high VIF value (above 5 or 10) should be removed from the model. This leads to a simpler model without compromising the model accuracy, which is good." http://www.sthda.com/english/articles/39-regression-model-diagnostics/160-multicollinearity-essentials-and-vif-in-r/#:~:text=For%20a%20given%20predictor%20(p,one%20(absence%20of%20multicollinearity).
#explanation of VIF here as well - https://online.stat.psu.edu/stat462/node/180/
car::vif(mod2.logit) #all VIF values are less than 5

#### checking for spatial autocorrelation in model residuals ####

#testing for the presence of spatial auto-correlation in model residuals using inverse of distance matrix
#calculate the distance between centroids of sets of sequences used for each species

#code from here - https://stats.oarc.ucla.edu/r/faq/how-can-i-calculate-morans-i-in-r/
#calculate the distance between centroids of sets of sequences used for each species
#converting the lat lon data into a matrix
coords<-seq.sum.sub[,c("lon.mean", "lat.mean")]
colnames(coords)<-c("long", "lat")
rownames(coords)<-seq.sum.sub$species
sp.dist<-geodist(x=coords, measure="haversine")/1000 #convert to km

#taking the inverse of distance to calculate weights
sp.dist.inv<-1/sp.dist
diag(sp.dist.inv)<-0

#calculate Moran's I
#note that model residuals are of different types, we are using the raw residuals here
Moran.I(residuals(mod2.logit, type="response"), sp.dist.inv, alternative="two.sided")
#Moran's I = 0.0726, p value = 0.0083

#make a correlogram to visualize Moran's I at different distance classes for the coordinates
set.seed(1287)
moran.cor<-correlog(x=coords$long, y=coords$lat, z=seq.sum.sub2$avg.pair.diff.tf, increment=500, resamp=999, latlon=TRUE)
plot(moran.cor, xlim=c(0,5000), ylim=c(-0.8, 0.2), cex=moran.cor$n/100, pch=ifelse(moran.cor$p<0.05, 19, 1))

#make a correlogram to visualize Moran's I at different distance classes for the model residuals
moran.cor.res<-correlog(x=coords$long, y=coords$lat, z=residuals(mod2.logit, type="response"), increment=500, resamp=999, latlon=TRUE)
points(y=moran.cor.res$correlation, x=moran.cor.res$mean.of.class, cex=moran.cor.res$n/100, xlim=c(0,5000), ylim=c(-0.8, 0.2), col="red", pch=ifelse(moran.cor.res$p<0.05, 19, 1))
lines(y=moran.cor.res$correlation, x=moran.cor.res$mean.of.class, col="red", type="l")

#using neighbourhood graph and distance-based weighting matrix 
#https://mgimond.github.io/simple_moransI_example/ was a useful reference for Moran's I test
#Appendix 8 of Benone et al., 2020 was a useful reference to build the neighbourhood matrix and calculate weights
#converting the lat lon data into a matrix
xy<-coords %>%
as.matrix

#correcting for the significant Moran's I
#computing graph-based connectivity schemes, not running Delaunay triangulation as it can give rise to long edges and works poorly with randomly sampled sites over fine-scale SAC (Bauman et al., 2018). Most of the code here is modified from Benone et al., 2020

#create a Gabriel graph using the mean of latitude and longitude for each species
nbgab<-graph2nb(gabrielneigh(xy), sym=TRUE)

#create a relative neighbourhood
nbrel<-graph2nb(relativeneigh(xy), sym=TRUE)

#create a minimum spanning tree - since these are geographic coordinates, using the haversine distance between coordinates instead of straight line distance
nbmst<-mst.nb(as.dist(geodist(x=coords, measure="haversine")))

#visualizing the above connectivity schemes
#convert nb object to listw object - performs row standardizations, all nodes are given the same weight
nbgab_listw<-nb2listw(nbgab)
nbrel_listw<-nb2listw(nbrel)
nbmst_listw<-nb2listw(nbmst)

#plot the above networks, creating a function to visualize the networks
nb2ggplot<-function(nb, coord){
#'coord' must be a matrix/dataframe with two columns (called "long" and "lat")
#'nb' is an object of class nb
#take out the connections from the nb object and assign them the lat and long in a dataframe
n<-length(attributes(nb$neighbours)$region.id)
DA<-data.frame(
from=rep(1:n, sapply(nb$neighbours, length)),
to=unlist(nb$neighbours)
)
DA<-DA[!(DA$to==0),]
DA<-cbind(DA, coord[DA$from, 1:2], coord[DA$to, 1:2])
colnames(DA)[3:6]=c("long","lat","long_to","lat_to")
return(DA)
}

#Gabriel graph
DA_gab<-nb2ggplot(nbgab_listw, xy)
gab_g<-xy %>%
as.data.frame() %>%
ggplot(aes(long, lat)) + geom_point() + coord_fixed() + geom_segment(data=DA_gab, aes(xend=long_to, yend=lat_to), size=0.3, alpha=0.5, colour= "darkred") + labs(title="Gabriel graph") + theme_bw()

#Relative neighbourhood
DA_rel<-nb2ggplot(nbrel_listw, xy)
rel_g<-xy %>%
  as.data.frame() %>%
  ggplot(aes(long, lat)) + geom_point() + coord_fixed() + geom_segment(data = DA_rel, aes(xend=long_to, yend=lat_to), size=0.3, alpha=0.5,  colour="chocolate") + labs(title="Relative neigh.") + theme_bw()

#minimum spanning tree
DA_mst<-nb2ggplot(nbmst_listw, xy)
mst_g<-xy %>%
  as.data.frame() %>%
  ggplot(aes(long, lat)) + geom_point() + coord_fixed() + geom_segment(data=DA_mst, aes(xend=long_to, yend=lat_to), size=0.3, alpha=0.5, colour="goldenrod") + labs(title="MST") + theme_bw()

#look at all the graphs at one go
pdf(file=paste0("figures/neighbourhood_graphs_11Apr23.pdf"), height=8, width=11, useDingbats=FALSE)
grid.arrange(gab_g, rel_g, mst_g, ncol=2)
dev.off()

#calculate the distance between connected points within each of the edited graphs, when longlat=TRUE - distance is measured in km
distgab<-nbdists(nbgab, xy, longlat=TRUE)
distrel<-nbdists(nbrel, xy, longlat=TRUE)
distmst<-nbdists(nbmst, xy, longlat=TRUE)

#calculating the weights for the connectivity matrices - linear decreasing function of distance
distgab<-lapply(distgab, function(x) 1-(x/max(unlist(distgab))))
distrel<-lapply(distrel, function(x) 1-(x/max(unlist(distrel))))
distmst<-lapply(distmst, function(x) 1-(x/max(unlist(distmst))))

#creating a listw object where the neighbourhood matrix is associated with weights
nbgab_w<-nb2listw(nbgab, glist=distgab, zero.policy=TRUE)
nbrel_w<-nb2listw(nbrel, glist=distrel, zero.policy=TRUE)
nbmst_w<-nb2listw(nbmst, glist=distmst, zero.policy=TRUE)

#binary forms of the above - not using the distance-based weights
nbgab_b<-nb2listw(nbgab, zero.policy=TRUE)
nbrel_b<-nb2listw(nbrel, zero.policy=TRUE)
nbmst_b<-nb2listw(nbmst, zero.policy=TRUE)

#defining the candidates created above
candidates<-list(gab_w=nbgab_w, gab_b=nbgab_b, rel_w=nbrel_w, rel_b=nbrel_b, mst_w=nbmst_w, mst_b=nbmst_b)

#select Morgan eigenvectors from the above candidates
select<-listw.select(residuals(mod2.logit, type="response"), candidates, MEM.autocor="all", method="MIR", p.adjust=TRUE)

#looking at the best models
select$best.id #rel_b is the best candidate model
select$candidates
select$best$summary #MEM13, 40 are selected as the eigenvectors that reduce spatial autocorrelation

#plot the spatial eigenvectors that are identified as being important
coastline<-readOGR('data_raw/ne_10m_coastline', layer='ne_10m_coastline')

#plot the spatial eigenvectors
pdf(file=paste0("figures/spatial_eigenvectors_11Apr23.pdf"), height=8, width=11, useDingbats=FALSE)
s.value(xy, cbind(select$best$MEM.select$MEM13, select$best$MEM.select$MEM40), symbol="circle", ppoint.cex=0.6, Sp=coastline)
dev.off()

#appending the spatial eigenvectors to the dataframe
seq.sum.sub2$MEM13<-select$best$MEM.select$MEM13
seq.sum.sub2$MEM40<-select$best$MEM.select$MEM40

#adding the selected MEMs to the model
mod2.logit.sac<-betareg(avg.pair.diff.tf ~ lat.mean + avg.geo.dist + sp.lat.range + body.size.new + blindness + maternal.care + MEM13 + MEM40 | n, data=seq.sum.sub2)
summary(mod2.logit.sac)

#looking for spatial autocorrelation in model residuals using inverse of distance as weights
Moran.I(residuals(mod2.logit.sac, type="response"), sp.dist.inv, alternative="two.sided") 
#two-sided because the local Moran's I plot revealed significant negative Moran's I values as well between 3000-4000 km
#Moran's I = 0.069, p value = 0.0111

#calculating local Moran's I for discrete distance classes
moran.cor.res2<-correlog(x=coords$long, y=coords$lat, z=residuals(mod2.logit.sac, type="response"), increment=500, resamp=999, latlon=TRUE)

pdf(file=paste0("figures/correlogram_distance_11Apr23.pdf"), height=8, width=11, useDingbats=FALSE)
#plot raw genetic diversity values
plot(moran.cor, xlim=c(0,5000), ylim=c(-0.8, 0.2), cex=moran.cor$n/100, pch=ifelse(moran.cor$p<0.05, 19, 1))

#plot the model residuals from beta regression using variables of interest
points(y=moran.cor.res$correlation, x=moran.cor.res$mean.of.class, cex=moran.cor.res$n/100, xlim=c(0,5000), ylim=c(-0.8, 0.2), col="red", pch=ifelse(moran.cor.res$p<0.05, 19, 1))
lines(y=moran.cor.res$correlation, x=moran.cor.res$mean.of.class, col="red", type="l")

#plotting it with previous correlograms
points(y=moran.cor.res2$correlation, x=moran.cor.res2$mean.of.class, xlim=c(0, 5000), ylim=c(-0.8, 0.2), col="green", pch=ifelse(moran.cor.res2$p<0.05, 19, 1), cex=moran.cor.res2$n/100)
lines(y=moran.cor.res2$correlation, x=moran.cor.res2$mean.of.class, col="green", type="l")

#add a horzintal line at y=0
abline(h=0)
dev.off()

#plot raw genetic diversity values and residuals from different models on the same map
#version 1 - this maybe better - uses point size to represent value instead of point colour
pdf(file=paste0("figures/residuals_map_11Apr23.pdf"), height=8, width=11, useDingbats=FALSE)
s.value(xy, cbind(seq.sum.sub2$avg.pair.diff.tf, residuals(mod2.logit, type="response"), residuals(mod2.logit.sac, type="response")), symbol="circle", ppoint.cex=0.6, Sp=coastline)
dev.off()

#### model diagnostics ####

#checking for multi-collinearity - "Any variable with a high VIF value (above 5 or 10) should be removed from the model. This leads to a simpler model without compromising the model accuracy, which is good." http://www.sthda.com/english/articles/39-regression-model-diagnostics/160-multicollinearity-essentials-and-vif-in-r/#:~:text=For%20a%20given%20predictor%20(p,one%20(absence%20of%20multicollinearity).
#can read the explanation of VIF here as well - https://online.stat.psu.edu/stat462/node/180/
car::vif(mod2.logit.sac) #all vifs are less than 5

#obtaining a table of coefficient estimates
mtable(mod2.logit.sac) #the value in brackets is the standard error associated with the parameter estimate

#looking at model diagnostics plots for model evaluation - need to look into the influence of some of the points which have high Cook's distance
#this resource is helpful in understanding some of the plots - https://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/R/R5_Correlation-Regression/R5_Correlation-Regression7.html
#"A leverage point is defined as an observation that has a value of x that is far away from the mean of x."
#"An influential observation is defined as an observation that changes the slope of the line. Thus, influential points have a large influence on the fit of the model. One method to find influential points is to compare the fit of the model with and without each observation"

pdf(file=paste0("figures/model_diagnostics_11Apr23.pdf"), height=11, width=8, useDingbats=FALSE)
par(mfrow=c(3,2))
plot(mod2.logit.sac, which=1:6)
dev.off()
#the half normal plot of residuals seem okay - the residuals largely fall within the bounds
#there seem to be some influential points based on the leverage vs predicted values plot and the Cook's distance plot

#plot the leverage plot separately with row ids
#leverage vs predicted values
diag<-data.frame(fitted=fitted(mod2.logit.sac), residuals=residuals(mod2.logit.sac), std.residuals=sqrt(abs(residuals(mod2.logit.sac))), leverage=gleverage(mod2.logit.sac), hat.diag=hatvalues(mod2.logit.sac))

pdf(file=paste0("figures/model_diagnostics_leverage_predicted_11Apr23.pdf"), height=8, width=11, useDingbats=FALSE)
ggplot(diag, aes(x=fitted, y=leverage)) + geom_point() + theme_bw() + geom_text_repel(label=rownames(diag)) + xlab("Predicted values") + ylab("Leverage")
dev.off()

#save row ids of species with high leverage
high.lev<-c(13,26,43,93,96,113,114)
seq.sum.sub2[high.lev,]
#species with high leverage are Schendyla nemorensis, Henicops maculatus, Lithobius lucifugus, Otostigmus astenus, Otostigmus rugulosus, Scolopendra mojianica, Scolopendra morsitans

#looking at cook's distance more closely. From - https://towardsdatascience.com/identifying-outliers-in-linear-regression-cooks-distance-9e212e9136a
cooksD<-cooks.distance(mod2.logit.sac)
influential<-cooksD[(cooksD>(3*mean(cooksD, na.rm=TRUE)))]
inf.id<-as.numeric(names(influential))
identical(influential, cooksD[names(cooksD) %in% inf.id])
seq.sum.sub2[rownames(seq.sum.sub2) %in% inf.id,] #the influential points come from Lithobius muticus, Otostigmus astenus, Scolopendra dehaani

#Cryptops parisi and Stenotaenia linearis are suspected to be a species complex - additionally removing these two species as well
sp.comp<-which(seq.sum.sub2$species %in% c("Cryptops parisi", "Stenotaenia linearis"))
seq.sum.sub2[sp.comp,]

#removing these influential points and running the regression again to see if anything changes
mod2.logit.sac.inf<-betareg(avg.pair.diff.tf ~ lat.mean + avg.geo.dist + sp.lat.range + body.size.new + blindness + maternal.care + MEM13 + MEM40 | n, data=seq.sum.sub2[-c(inf.id, high.lev, sp.comp),])
summary(mod2.logit.sac.inf)
#the relationships are much stronger and the R-squared value jumps when influential points are removed - the inferences remain the same, blindness is marginally significant

#looking at the model evaluation plots again
pdf(file=paste0("figures/model_diagnostics_infl_11Apr23.pdf"), height=11, width=8, useDingbats=FALSE)
par(mfrow=c(3,2))
plot(mod2.logit.sac.inf, which=1:6)
dev.off()

#since some of the points seem to be influential, instead of dropping the influential points, obtaining bootstrapped coefficients

#### bootstrapping coefficients to reduce the effect of influential points on coefficient estimates ####

#using the Boot package to obtain bootstrapped estimates from here - https://socialsciences.mcmaster.ca/jfox/Books/Companion-2E/appendix/Appendix-Bootstrapping.pdf. Some notes for the type of the confidence intervals is here - https://stats.stackexchange.com/questions/355781/is-it-true-that-the-percentile-bootstrap-should-never-be-used#:~:text=Yes%2C%20if%20we%20can%2C%20we,and%20think%20we%20discovered%20America.
set.seed(3453)
b<-car::Boot(mod2.logit.sac, R=1000)

#creating a function to convert the results into a dataframe
extract.coef<-function(boot.res){
#create an empty dataframe to store results
boot.res.ci<-data.frame(matrix(nrow=ncol(boot.res$t), ncol=4))
colnames(boot.res.ci)<-c("predictor", "beta", "ci.lower", "ci.upper")
#extract confidence intervals
for(i in 1:ncol(boot.res$t)){
boot.int<-boot.ci(boot.res, type="perc", index=i)
boot.res.ci[i,"predictor"]<-attributes(boot.int$t0)$names
boot.res.ci[i,"beta"]<-as.numeric(boot.int$t0)
boot.res.ci[i,"ci.lower"]<-boot.int$percent[4]
boot.res.ci[i,"ci.upper"]<-boot.int$percent[5]
rm(boot.int)
}
boot.res.ci
}

#create a dataframe with the results 
es<-extract.coef(b)
es$sig<-c(1,1,1,0,1,0,1,1,1,1,1)

#to interpret the betas we need to exponentiate them - since the regression is presented as a logit transformation
es$exp.beta<-exp(es$beta)
es$exp.ci.lower<-exp(es$ci.lower)
es$exp.ci.upper<-exp(es$ci.upper)

#plotting standardized coefficients
es$predictor<-factor(es$predictor, levels=c("(Intercept)", "lat.mean", "avg.geo.dist",  "sp.lat.range", "body.size.new", "blindnessnot blind", "maternal.careyes", "MEM13", "MEM40", "(phi)_(Intercept)", "(phi)_n"))

levels(es$predictor)<-c("Intercept", "Latitude", "Geographic distance", "Latitudinal range", "Body size", "Vision:Yes", "Maternal care:Yes", "MEM13", "MEM40", "Dispersion:Intercept", "Dispersion:Slope:N")

es$sig<-factor(es$sig)

#not plotting the dispersion factor related parameters
es.plot<-es[-c(1,8:12),]
es.plot<-droplevels(es.plot)

my.col<-scale_color_manual(name="Significance", labels=c("Not significant", "Significant"), values=c("#bdbdbd", "black"))

#plot effects of standardized coefficients without exponentiating
ggplot(es.plot, aes(x=beta, y=predictor, xmin=ci.lower, xmax=ci.upper, col=sig)) + geom_vline(xintercept=0, color="red") + geom_point() + geom_errorbarh(height=.2) + my.col + theme_bw() + xlab("Standardized coefficients") + ylab("Predictor variables") + theme(text=element_text(size=20)) 

#A 1SD change in body size leads to a relative change of exp(-3.36)=0.035 units change in E(Proportion)/(1−E(Proportion)). Need to interpret this in some sensible way - what we are interested in is the Proportion, but what we are having to interpret is some transformation of this.

#### plot effect size ####

#the predictorEffects function from the effects package only taking in glmmTMB models as input - therefore running the model again using glmmTMB to plot individual predictor effects

mod2.sac<-glmmTMB(avg.pair.diff.tf ~ lat.mean + avg.geo.dist + sp.lat.range + body.size.new + blindness + maternal.care + MEM13 + MEM40, dispformula=~n, data=seq.sum.sub2, list(family="beta", link="logit"))
summary(mod2.sac)

#plot the effect of variation in body size
plot(predictorEffects(mod2.sac, "body.size.new", residuals=FALSE), axes=list(y=list(type="response", cex=1.5), x=list(cex=1.5)), main="", xlab="Body size", ylab="Genetic diversity", cex.lab=4, cex.axis=4, cex.main=4, cex.sub=4)

#plot the effect of maternal care
plot(predictorEffects(mod2.sac, "maternal.care", residuals=FALSE), axes=list(y=list(type="response", cex=1.5), x=list(cex=1.5)), main="", xlab="Genetic diversity", ylab="Maternal care", cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3)

#plot the effect of average geographic distance
plot(predictorEffects(mod2.sac, "avg.geo.dist", residuals=FALSE), axes=list(y=list(type="response", cex=1.5), x=list(cex=1.5)), main="", xlab="Average geographic distance", ylab="Genetic diversity", cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3)

#plot the effect of variation in mean latitude
trellis.par.set(list(label=list(cex=1.5)))
plot(predictorEffects(mod2.sac, "lat.mean", residuals=FALSE), axes=list(y=list(type="response", cex=1.5), x=list(cex=1.5)), main="", xlab="Mean latitude (degrees)", ylab="Genetic diversity")

#### test for phylogenetic signal in model residuals ####

#input the family level phylogenetic tree
cen.tree<-ape::read.nexus("data_raw/Node_calibration.nexus.tre")
plot(cen.tree, cex=0.5)

#drop branches that do not correspond to centipedes
#look at the label names of the tips to decide which ones need to be removed
sort(cen.tree$tip.label)

#remove the tips which are not centipedes
drop.sp<-c("Abacion", "Cleidogona", "Cambala", "Cylindroiulus", "Narceus", "Prostemmiulus", "Pseudopolydesmus", "Brachycybe", "Petaserpes", "Cyliosoma", "Glomeridesmus", "Glomeris", "Polyxenida", "Scutigerella", "Symphyla_NZ", "Calanus", "Daphnia", "Drosophila", "Anoplodactylus", "Centruroides", "Damon", "Mastigoproctus", "Liphistius", "Metasiro", "Limulus", "Peripatopsis_long_iso")
cen.tree2<-drop.tip(cen.tree, drop.sp)
plot(cen.tree2)

#replace the genus names with the family which it belongs to 
cen.tree2$tip.label
family.label<-c("Pselliodidae", "Scutigerinidae", "Scutigeridae", "Craterostigmidae", "Craterostigmidae", "Lithobiidae", "Lithobiidae", "Henicopidae", "Henicopidae", "Mecistocephalidae", "Mecistocephalidae", "Schendylidae", "Himantariidae", "Oryidae", "Geophilidae", "Geophilidae", "Geophilidae", "Plutoniumidae", "Scolopocryptopidae", "Scolopocryptopidae", "Cryptopidae", "Scolopendridae", "Scolopendridae", "Scolopendridae", "Scolopendridae")
length(unique(family.label)) #15 unique families

#check if there are any differences in the family names present in the phylogenetic tree vs. those present in our data
setdiff(unique(seq.sum.sub2$family), family.label) #Mimopidae is missing in the tree
setdiff(family.label, unique(seq.sum.sub2$family)) #Pselliodidae, Himantariidae and Oryidae are present in the tree but missing from our data

#replace tip names
cen.tree3<-cen.tree2
cen.tree3$tip.label<-family.label

#compare the original and renamed tree next to each other, to check if the assignments were made correctly
par(mfrow=c(1,2))
plot(cen.tree2)
plot(cen.tree3)
par(mfrow=c(1,1))

#removing multiple instances of the same family and families not present in our data from the tree by removing specific genera
drop.sp<-c("Akymnopellis", "Scolopendropsis", "Alipes", "Newportia", "Henia", "Stenotaenia", "Notiphilides", "Himantarium", "Tygarrup", "Anopsobius", "Eupolybothrus", "C_crabilli", "Sphendononema")
cen.tree4<-drop.tip(cen.tree2, drop.sp)
plot(cen.tree4)

#replace the genus names with the family which it belongs to 
cen.tree4$tip.label
cen.tree4$tip.label<-c("Scutigerinidae", "Scutigeridae", "Craterostigmidae", "Lithobiidae", "Henicopidae", "Mecistocephalidae", "Schendylidae", "Geophilidae", "Plutoniumidae", "Scolopocryptopidae", "Cryptopidae", "Scolopendridae")

#plot the different versions of the tree together to cross-check if everything looks okay
par(mfrow=c(1,3))
plot(cen.tree2)
plot(cen.tree3)
plot(cen.tree4)
par(mfrow=c(1,1))

#obtain model residuals from the beta regression model that was run
res<-residuals(mod2.logit.sac, type="response")
names(res)<-as.character(seq.sum.sub2$family)
res.df<-data.frame(res=residuals(mod2.logit.sac, type="response"), family=as.character(seq.sum.sub2$family))
head(res.df)
res.df.sum<-aggregate(res.df$res, by=list(res.df$family), mean)
res.mean<-res.df.sum$x
names(res.mean)<-res.df.sum$Group.1

#only subsetting values that correspond to tips on the tree
res.mean2<-res.mean[names(res.mean) %in% cen.tree4$tip.label]

#ordering the data based on the order of tip labels on the tree
res.mean2<-res.mean2[cen.tree4$tip.label]

#test for phylogenetic signal
#the data must have the same taxon names as the phylogenetic tree - the tree needs to have the same number of tips
phylosig(cen.tree4, res.mean2, method="lambda", test=TRUE)
#lambda = 6.6107*10^-5, p-value=1

#### figures and tables ####

## table 1 - summary of beta regression ##
#check if the parameter names are in the same order
cbind(tidy(mod2.logit.sac)$term, as.character(es$predictor))

betareg_res<-cbind(tidy(mod2.logit.sac), es[,c("ci.lower", "ci.upper")])
betareg_tbl<-betareg_res %>%
mutate(term=case_match(term, '(Intercept)'~'Intercept', 'lat.mean'~'Mean latitude', 'avg.geo.dist'~'Average geographic distance', 'sp.lat.range'~'Species latitudinal range', 'body.size.new'~'Body size', 'blindnessnot blind'~'Vision: Yes', 'maternal.careyes'~'Maternal care: Yes', 'MEM13'~'MEM 13', 'MEM40'~'MEM 40', 'n'~'Number of sequences')) %>%
mutate(component=case_match(component, 'mean'~'Mean', 'precision'~'Precision')) %>%
mutate(across(c(3:5, 7:8), round, 3)) %>%
mutate(across(c(3:5, 7:8), format, nsmall=3)) %>%
mutate(boot.ci=paste0(ci.lower, " - ", ci.upper)) %>%
mutate(estimate=paste0(estimate, stars.pval(p.value))) 

#obtain model summary stats
footer<-glance(mod2.logit.sac)

#create gt table for manuscript
betareg_gt<-betareg_tbl %>%
select(component, term, estimate, std.error, statistic, boot.ci) %>%
gt(groupname_col="component") %>%
row_group_order(groups=c("Mean", "Precision")) %>%
cols_align(align="center", columns=c(std.error, statistic, boot.ci)) %>%
cols_align(align="left", columns=estimate) %>%
cols_label(term="Parameters", estimate="Estimate", std.error="SE", statistic="z-value", boot.ci="Bootstrap 95% CI") %>%
cols_move(boot.ci, after=std.error) %>%
tab_footnote(paste0("Pseudo R-squared = ", round(footer$pseudo.r.squared, 4))) %>%
tab_footnote(paste0("Log-likelihood = ", round(footer$logLik, 4))) %>%
tab_footnote(paste0("N = ", footer$nobs)) %>%
tab_footnote("*** = p < 0.001, ** = p < 0.01, * = p < 0.05, . = p < 0.1", locations=cells_column_labels(columns=estimate)) %>%
tab_footnote("Standard error", locations=cells_column_labels(columns=std.error)) %>%
tab_options(footnotes.font.size="16px")
  
betareg_gt %>%
gtsave("results/table1_model_summary_11Apr23.rtf")

## fig 2. plot of sequence locations ##
#loading information of locations
acc.sp<-read.csv("data/COI_final_locations_11Apr23.csv", header=TRUE)

#keep only the native occurrences
acc.sp<-acc.sp %>%
filter(Introduced==0) %>%
as.data.frame

#loading the database of only coordinate associated occurrences from GBIF
dat.gbif<-read.table("data_raw/GBIF_0420414-210914110416597/0420414-210914110416597.csv", sep='\t', quote="", header=TRUE, fill=TRUE)
colnames(dat.gbif)

#remove problematic coordinates
flags<-clean_coordinates(x=dat.gbif[!is.na(dat.gbif$decimalLatitude),], lon="decimalLongitude", lat="decimalLatitude", species="species", tests=c("capitals", "centroids", "equal", "gbif", "institutions", "zeros"))

#removing the flagged entries and those with coordinate uncertainty greater than 50 km
clean.gbifID<-flags$gbifID[flags$.summary]
dat.gbif.clean<-dat.gbif %>%
filter(!is.na(decimalLatitude)) %>%
filter(gbifID %in% clean.gbifID) %>%
filter(!(coordinateUncertaintyInMeters > 50000)) %>%
as.data.frame()
nrow(dat.gbif.clean) #78007

#write gbif.dat to disk
write.csv(dat.gbif.clean, "data/gbif_chilopoda_species_records_11Apr23.csv")

#input the coastline file
coastline<-readOGR('data_raw/ne_10m_coastline', layer='ne_10m_coastline')

#converting the spatial lines object into a spatial lines dataframe for plotting
world.sdf<-SpatialLinesDataFrame(coastline, coastline@data)

#colouring points based on order
my.theme<-theme(axis.title.x=element_text(size=16), axis.text.x=element_text(size=14), axis.title.y=element_text(size=16), axis.text.y=element_text(size=14))

acc.sp$Order<-as.factor(acc.sp$Order)

my.fill<-scale_fill_colorblind(name="Order", labels=levels(acc.sp$Order))

pdf(file=paste0("figures/fig2_COI_locations_colorblind_11Apr23.pdf"),  height=8, width=11, useDingbats=FALSE)
plot.dat<-ggplot() + geom_path(data=world.sdf, aes(x=long, y=lat, group=group), size=0.2) + geom_point(data=dat.gbif.clean, aes(x=decimalLongitude, y=decimalLatitude), colour="gray", alpha=0.5, size=0.1)
plot.dat + geom_point(data=acc.sp, aes(x=x, y=y, fill=Order), col="black", pch=21, size=2) + xlab("Longitude") + ylab("Latitude") + theme_bw() + geom_jitter() + scale_x_continuous(breaks=seq(-180, 180, by=20)) + scale_y_continuous(breaks= seq(-90, 90, by=20)) + coord_fixed() + my.theme + my.fill
dev.off()

## fig 3. genetic diversity distribution of arthropods ##
nd.comp<-read.csv("data_raw/genetic_diversity_arthropods_no_introductions_11Apr23.csv", header=TRUE)
str(nd.comp)
colnames(nd.comp)
nd.comp$avg.pair.diff<-as.numeric(as.character(nd.comp$avg.pair.diff))

#some of the genetic diversity values from Pelletier & Carstens, 2018 have NA values - removing those
nd.comp<-nd.comp[!is.na(nd.comp$avg.pair.diff), 1:12]

#replacing "_" with a space for species names from Dapporto et al., 2019
nd.comp<-nd.comp %>%
mutate(species.1=gsub("_", " ", species.1))

#checking if there are any duplicated species
nd.comp %>%
group_by(species.1) %>% 
filter(n()>1)

#Dapporto et al., 2019 and Pelletier & Carstens, 2018 share many species - retaining data from Pelletier & Carstens, 2018 as the other study is limited to Europe
lepi.dups<-nd.comp %>% 
filter(order=="Lepidoptera") %>%
group_by(species.1) %>%
filter(n()>1) %>%
distinct(species.1) %>%
pull()

#removing centipede species names from Pelletier & Carstens, 2018 from our data
chilo.dups<-nd.comp %>% 
filter(class=="Chilopoda") %>%
group_by(species.1) %>%
filter(n()>1) %>%
distinct(species.1) %>%
pull()

#remove rows with the reference Pelletier & Carstens, 2018 if they have Chilopoda species names already present in our data
nd.comp2<-nd.comp %>%
filter(!(species.1 %in% chilo.dups & reference=="Pelletier & Carstens, 2018")) %>%
filter(!(species.1 %in% lepi.dups & reference=="Dapporto et al., 2019")) 

#checking if there are any other duplicates
nd.comp2 %>%
group_by(species.1) %>% 
filter(n()>1) %>%
as.data.frame()
#in the Pelletier & Carstens, 2018 study, duplicate species names represent different populations - keeping these as is

#find mean and range for the different arthropod classes
nd.comp2 %>%
group_by(class) %>%
summarise(mean.gd=mean(avg.pair.diff), min.gd=min(avg.pair.diff), max.gd=max(avg.pair.diff), n=n()) %>%
as.data.frame()

#convert class to a factor
nd.comp2$class<-as.factor(nd.comp2$class)

my.fill2<-scale_fill_manual(name="Class", labels=levels(nd.comp2$class), values=c("#E69F00", "#F0E442", "#009E73", "#56B4E9", "#CC79A7"))

#plotting the genetic diversity values across arthropod classes
pdf(file=paste0("figures/fig3_genetic_diversity_arthropods_no_introductions_colorblind_11Apr23.pdf"), height=6, width=11, useDingbats=FALSE)
ggplot(data=nd.comp2, aes(x=avg.pair.diff, y=..scaled.., fill=class)) + geom_density(alpha=0.5) + theme(axis.text=element_text(size=20), axis.title=element_text(size=20), panel.background=element_rect(fill="white", colour="black"), panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + xlab("Average pairwise difference") + ylab("Scaled density") + my.fill2
dev.off()

## fig 4. effect size of model coefficients  ##
my.col<-scale_color_manual(name="Significance", labels=c("Not significant", "Significant"), values=c("#bdbdbd", "black"))

#plot effects of standardized coefficients without exponentiating
pdf(file=paste0("figures/fig4_effect_size_unexp_no_int_sac_11Apr23.pdf"), height=8, width=11, useDingbats=FALSE)
ggplot(es.plot, aes(x=beta, y=predictor, xmin=ci.lower, xmax=ci.upper, col=sig)) + geom_vline(xintercept=0, color="red") + geom_point() + geom_errorbarh(height=.2) + my.col + theme_bw() + xlab("Standardized coefficients") + ylab("Predictor variables") + theme(text=element_text(size=20)) 
dev.off()

## fig 5. individual effect size plots ##
pdf(file=paste0("figures/fig5a_effect_body_size_no_int_sac_11Apr23.pdf"), height=8, width=11, useDingbats=FALSE)
plot(predictorEffects(mod2.sac, "body.size.new", residuals=FALSE), axes=list(y=list(type="response", cex=1.5), x=list(cex=1.5)), main="", xlab="Body size", ylab="Genetic diversity", cex.lab=4, cex.axis=4, cex.main=4, cex.sub=4)
dev.off()

#plot the effect of maternal care
pdf(file=paste0("figures/fig5b_effect_maternal_care_no_int_sac_11Apr23.pdf"), height=8, width=11, useDingbats=FALSE)
plot(predictorEffects(mod2.sac, "maternal.care", residuals=FALSE), axes=list(y=list(type="response", cex=1.5), x=list(cex=1.5)), main="", ylab="Genetic diversity", xlab="Maternal care", cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3)
dev.off()

#plot the effect of average geographic distance
pdf(file=paste0("figures/fig5c_effect_distance_no_int_sac_11Apr23.pdf"), height=8, width=11, useDingbats=FALSE)
plot(predictorEffects(mod2.sac, "avg.geo.dist", residuals=FALSE), axes=list(y=list(type="response", cex=1.5), x=list(cex=1.5)), main="", xlab="Average geographic distance", ylab="Genetic diversity", cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3)
dev.off()

#plot the effect of variation in mean latitude
trellis.par.set(list(label=list(cex=1.5)))
pdf(file=paste0("figures/fig5d_effect_latitude_no_int_sac_11Apr23.pdf"), height=8, width=11, useDingbats=FALSE)
plot(predictorEffects(mod2.sac, "lat.mean", residuals=FALSE), axes=list(y=list(type="response", cex=1.5), x=list(cex=1.5)), main="", xlab="Mean latitude (degrees)", ylab="Genetic diversity")
dev.off()

#### statistics in main text - methods and results ####
#total number of sequences used on analysis
sum(seq.sum.sub$n) #1245 sequences

#number of species
nrow(seq.sum.sub) #128 species

#range of genetic diversity
range(seq.sum.sub$avg.pair.diff) # 0 to 0.1713

#total number of published studies and sequence datasets
#input raw data
dat<-read.csv(file="data_raw/sheet4_popgen_database_analysis_clean_edit_11Apr23.csv", header=TRUE)

sort(unique(dat$Reference))

#look at the published studies - filter out references where 'unpublished' is mentioned
toMatch<-c("unpublished", "in prep")

dat %>%
filter(!grepl(paste(toMatch, collapse="|"), Reference)) %>%
distinct(Reference) %>%
pull(Reference) %>%
sort() #50 published studies

#unpublished datasets
dat %>%
filter(grepl(paste(toMatch, collapse="|"), Reference)) %>%
distinct(Reference) %>%
pull(Reference) %>%
sort() #11 unpublished datasets - two of the entries consist of a combination of the existing references

#number of families
length(unique(seq.sum.sub$family)) #13 families

#reproductive strategy - number of species showing maternal care
seq.sum.sub %>%
group_by(maternal.care) %>%
summarise(n()) #83 species with maternal care of 128 total species

#vision - number of species with vision
seq.sum.sub %>%
group_by(blindness) %>%
summarise(n()) #98 species with vision of 128 total species

#mean and range of body size, number of sequences
seq.sum.sub %>%
summarise(across(c(body.size.new, n, seq.len, n.unique.sites, avg.geo.dist, avg.pair.diff), list(mean=mean, min=min, max=max), .names="{col}_{fn}"))
#body size: mean=48 mm, min=8.5 mm, max=250 mm
#n sequences: mean=9.73, min=3, max=68
#alignment length: mean=648, min=465, max=840
#unique sites: mean=7.5, min=1, max=53
#average geographic distance: mean=439.95 km, min=0, max=5065.76

#statistics from all sequence locations
acc.final<-acc.sp %>%
filter(Introduced=="0") %>%
filter(Species %in% seq.sum.sub$species)
nrow(acc.final) #1245 locations
sum(seq.sum.sub$n) #1245 sequences

#total number of unique sequence locations
nrow(unique(acc.final[,c("x", "y")])) #774 unique geographic locations

#degree span of sequence locations
min(acc.final$y) #min latitude=-46.89
max(acc.final$y) #max latitude=60.46
max(acc.final$y)-min(acc.final$y) #107.36 degrees of latitudinal span

#mean latitude by order
acc.final %>%
group_by(Order) %>%
summarise(mean.lat=mean(y)) %>%
as.data.frame()
#Scolopendromorpha=18.73 and Lithobiomorpha=45.62

#number of sequences - northern vs southern hemisphere
acc.final %>%
filter(y>0) %>%
nrow() #1072 sequences from northern hemisphere

acc.final %>%
filter(y<0) %>%
nrow() #173 sequences from the southern hemisphere

#mean and range genetic diversity of centipedes
seq.sum.sub %>%
summarise(mean(avg.pair.diff), min(avg.pair.diff), max(avg.pair.diff))
#genetic diversity: mean=0.0721, min=0, max=0.1713

#find mean and range for the different arthropod classes
mean.range<-nd.comp2 %>%
group_by(class) %>%
summarise(mean.gd=mean(avg.pair.diff), min.gd=min(avg.pair.diff), max.gd=max(avg.pair.diff), n=n()) %>%
as.data.frame()
#insects: mean gd=0.0098, millipedes: mean.gd=0.0445

#Moran's I value
Moran.I(residuals(mod2.logit, type="response"), sp.dist.inv, alternative="two.sided")
#Moran's I = 0.0726, p value = 0.0083

#SAC corrected model % variance explained
summary(mod2.logit.sac)$pseudo.r.squared #Pseudo r-squared=0.2757 

#lambda value of phylogenetic signal and its p value
phylosig(cen.tree4, res.mean2, method="lambda", test=TRUE)
#lambda = 6.6107*10^-5, p-value=1

#4 sequence cut-off - total number of species
seq.sum.sub %>%
filter(n>3) %>%
nrow() #91 species

#number of Geophilomorpha species sampled
seq.sum.sub %>%
filter(order=="Geophilomorpha") %>%
summarise(n()) #17 species

#### supplementary tables and figures ####
## supplementary table S3.1 ##
seq.sum.tbl<-seq.sum.sub %>%
select(order, family, n, seq.len, n.unique.sites, avg.pair.diff, body.size.new, maternal.care, blindness, sp.lat.range, lat.mean, avg.geo.dist) %>%
rename(blindness=vision) %>%
mutate(family=as.numeric(as.factor(family))) %>%
mutate(vision=case_match(vision, "not blind"~"yes", "blind"~"no")) %>%
tbl_summary(by=order, statistic=list(c(avg.pair.diff, lat.mean, avg.geo.dist, sp.lat.range, body.size.new, n, n.unique.sites, seq.len)~paste0("{mean}", "\n", "({min} - {max})"), c(vision, maternal.care)~"{n} / {N}", family~"{n_distinct}"), digits=list(c(avg.pair.diff, lat.mean, avg.geo.dist, sp.lat.range, body.size.new)~2, family~0), label=list(avg.pair.diff~"Average pairwise difference per species", lat.mean~"Mean latitude of sequence data per species (degrees)", avg.geo.dist~"Mean geographic distance between sequences per species (km)", sp.lat.range~"Mean latitudinal range per species (degrees)", body.size.new~"Mean body size per species (mm)", n~"Mean number of sequences per species", n.unique.sites~"Mean number of unique locations per species", seq.len~"Mean alignment length per species (bp)", vision~"Vision: Present", maternal.care~"Maternal care: Present", family~"No. families")) %>%
modify_header(label="") %>%
add_overall() %>%
modify_footnote(update=everything()~NA) %>%
as_gt() %>%
gt::fmt_markdown(columns = everything())

seq.sum.tbl %>%
gtsave("results/tables3.1_data_summary_11Apr23.docx")

## supplementary figure S3.1 ##

#loading information of locations
acc.sp<-read.csv("data/COI_final_locations_11Apr23.csv", header=TRUE)

#set plot theme
my.theme<-theme(axis.title.x=element_text(size=16), axis.text.x=element_text(size=14), axis.title.y=element_text(size=16), axis.text.y=element_text(size=14))

#convert the Introduced column to a factor so that it can be plotted by symbol
acc.sp$Introduced<-as.factor(as.character(acc.sp$Introduced))
levels(acc.sp$Introduced)<-c("No", "Yes")
acc.sp$Introduced<-factor(acc.sp$Introduced, c("Yes", "No"))

acc.sp$Order<-as.factor(acc.sp$Order)

#using a colorblind colour palette
my.fill<-scale_fill_manual(name="Order", values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442"), labels=levels(acc.sp$Order))
my.shape<-scale_shape_manual(name="Introduction", values=c(24, 21), labels=levels(acc.sp$Introduced))

pdf(file=paste0("figures/COI_locations_int_colorblind_11Apr23.pdf"),  height=8, width=11, useDingbats=FALSE)
plot.dat<-ggplot() + geom_path(data=world.sdf, aes(x=long, y=lat, group=group)) + geom_point(data=dat.gbif.clean, aes(x=decimalLongitude, y=decimalLatitude), colour="gray", alpha=0.5, size=0.1)
plot.dat + geom_point(data=acc.sp, aes(x=x, y=y, fill=Order, shape=Introduced), col="black", size=2) + xlab("Longitude") + ylab("Latitude") + theme_bw() + geom_jitter() + scale_x_continuous(breaks=seq(-180, 180, by=20)) + scale_y_continuous(breaks= seq(-90, 90, by=20)) + coord_fixed() + my.theme + my.fill + my.shape + guides(fill=guide_legend(override.aes=list(shape=21), order=1), shape=guide_legend(order=2))
dev.off()

####







