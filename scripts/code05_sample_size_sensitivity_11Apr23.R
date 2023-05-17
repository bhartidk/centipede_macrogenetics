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
library(ggrepel)
library(gridExtra)
library(rgdal)
library(adegraphics)
library(memisc)
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

#### resampling sequence alignments based on sample size ####
#resampling procedure based on Barrow et al., 2021

#load input data
#input the master sequence database
dat<-read.csv(file="data_raw/sheet4_popgen_database_analysis_clean_edit_11Apr23.csv", header=TRUE)
nrow(dat) #1429 rows

#removing rows that do not have to be included in analysis
dat<-dat %>%
filter(!is.na(Latitude)) %>%
filter(!is.na(COI)) %>%
filter(COI!="") %>%
filter(Introduced==0) %>% #adding this condition here itself
filter(!(Species %in% c("Cormocephalus sp.", "Scolopendra calcarata", "Scolopendra subspinipes"))) %>%
group_by(Species) %>%
filter(n()>2)
nrow(dat) #1244 rows

#creating a unique list of species
sp<-as.character(unique(dat$Species))
length(sp) #122 species

#load phylogatr database for missing species
#load the dataframe for missing species from phylogatr
dat.pg.missing<-read.csv("data/phylogatr_missing_species_11Apr23.csv", header=TRUE)

#save unique species names; remove Anopsobius neozelanicus because it may not represent a single species based on existing literature
missing.sp<-dat.pg.missing %>%
filter(species!="Anopsobius neozelanicus") %>%
distinct(species) %>%
pull()

#combining species names
all.sp<-c(sp, missing.sp)

#create a dataframe to store results
n<-100
sample.size<-10
pi.sampled<-data.frame()
file.aln<-"data/seq_aln_final/"

#looping across species alignments to save sequence summary statistics
for(i in 1:length(all.sp)){
#finding the name of the species
sp.i<-all.sp[i]
#assigning the file name based on whether this species is in the previous lot of species or from the current lot from phylogatr, where alignments are stored in different
if(sp.i %in% missing.sp){
#creating the file name for the alignment pointer from the species name
file.sp<-paste0(gsub(" ", "-", all.sp[i]), "-COI.afa")

#reading in the aligned sequences for sp[h]
aln.sub<-read.dna(paste0(file.aln, "phylogatR_species_only/", file.sp), format="fasta")
}else
{
#creating the file name for the alignment pointer from the species name
file.sp<-paste0("aln_", gsub(" ", "_", all.sp[i]), ".fasta")

#find the accession numbers of the sequences for which location data is available
acc<-dat %>%
filter(Introduced != 1) %>%
filter(Species %in% sp[i]) %>%
pull(COI)

#reading in the aligned sequences for sp[h]
aln.sub<-read.dna(paste0(file.aln, "/", file.sp), format="fasta")

#take a subset of aln.sub based on what accession numbers have location data
ids<-which(labels(aln.sub) %in% acc)
aln.sub<-aln.sub[ids,]
}

#calculating nucleotide diversity of the original dataset
for(j in 2:sample.size){
for(k in 1:n){
#sample fasta
sampled_fasta<-aln.sub[sample(1:length(labels(aln.sub)), j, replace=TRUE),]
#calc nucleotide diversity of the sampled dataset
pi.new<-pegas::nuc.div(sampled_fasta, variance=FALSE, pairwise.deletion=FALSE)
pi.row<-c(sp.i, j, pi.new)
pi.sampled<-rbind(pi.sampled, pi.row)
}
}
}

colnames(pi.sampled)<-c("species", "sample.size", "pi.sample")

#checking if the number of rows match what we expect
nrow(pi.sampled)==(length(all.sp)*length(c(2:10))*100) #TRUE

#writing pi.sampled to disk
write.csv(pi.sampled, file="data/pi_sampled_11Apr23.csv",row.names=FALSE)

#since the summary comes from random sub-samples of an alignment, reruns of the code would give slightly different results. So writing the results to disk and reloading it from there again.
pi.sampled<-read.csv("data/pi_sampled_11Apr23.csv", header=TRUE)
str(pi.sampled)

#getting statistics out of the dataframe
pi.sampled.sum<-pi.sampled %>%
group_by(species, sample.size) %>%
summarise(pi.variance=var(pi.sample)) %>%
mutate(sample.size=as.numeric(sample.size)) %>%  
as.data.frame()
head(pi.sampled.sum)

#looping through species names and saving the plots of sample.size vs variance(pi) to disk
pdf(file=paste0("figures/pi_variance_species_wise_11Apr23.pdf"), height=8, width=11, useDingbats=FALSE)
for(i in 1:length(all.sp)){
#subset the species of interest
pi.sp<-pi.sampled.sum %>%
filter(species==all.sp[i])

#plot variance of pi against sample size
plot(pi.sp$pi.variance ~ pi.sp$sample.size, type='b', xlab="sample size", ylab="pi variance", main=substitute(italic(x), list(x=all.sp[i])))
}
dev.off()

#plot all graphs in the same plot
pdf(file=paste0("figures/pi_variance_all_species_11Apr23.pdf"), height=8, width=11, useDingbats=FALSE)
#make the plot for the first species
pi.sp.1<-pi.sampled.sum %>%
filter(species==all.sp[1])
plot(pi.sp.1$pi.variance ~ pi.sp.1$sample.size, type='l', xlab="sample size", ylab="pi variance", col=alpha("brown", 0.2), ylim=c(0,0.01))

for(i in 2:length(all.sp)){
#subset the species of interest
pi.sp<-pi.sampled.sum %>%
filter(species==all.sp[i])

#plot variance of pi against sample size
lines(pi.sp$pi.variance ~ pi.sp$sample.size, type='l', xlab="sample size", ylab="pi variance", col=alpha("brown", 0.2), ylim=c(0,0.01))
}
dev.off()

#from the graphs it is apparent that the variance in nucleotide diversity estimates comes down considerably at a sample size of 4 - therfore, taking the median nucleotide diversity from 100 random 4 sequence subsamples from the original dataset to rerun analysis
pi.median.4n<-pi.sampled %>%
filter(sample.size==4) %>%
group_by(species) %>%
summarise(pi.median=median(pi.sample)) %>%
as.data.frame()
which(!complete.cases(pi.median.4n)) #complete cases

#input sequence summary
seq.sum<-read.csv(file="data/input_betareg_no_introductions_11Apr23.csv", header=TRUE, na.strings=c("","NA"))

#number of species with complete information
nrow(seq.sum[complete.cases(seq.sum),]) #128

#create a column with latitudinal range
seq.sum$sp.lat.range<-seq.sum$int.max.lat-seq.sum$int.min.lat

#taking a subset of seq.sum which have non-NA data
seq.sum.sub<-seq.sum[,c("species", "family", "order", "avg.pair.diff", "lat.mean", "lon.mean", "avg.geo.dist", "sp.lat.range", "body.size.new", "blindness", "maternal.care", "n", "n.unique.sites", "seq.len", "hapl.div")]

#adding a new column to seq.sum.sub with median pi values
seq.sum.sub.4nmed<-left_join(seq.sum.sub, pi.median.4n, by="species")
seq.sum.sub.4nmed<-seq.sum.sub.4nmed[seq.sum.sub.4nmed$n>3,]
nrow(seq.sum.sub.4nmed) #91 rows

#response values with 0 cannot be modeled by beta regression - using a modification function to deal with this - code from Douma et al., 2019
transform01<-function(x){
((x*(length(x)-1))+0.5)/(length(x))
}

#creating a new column with the transformed genetic diversity values
seq.sum.sub.4nmed$pi.median.tf<-transform01(seq.sum.sub.4nmed$pi.median)

#correlation between median and estimated value of full dataset
cor.test(seq.sum.sub.4nmed$pi.median, seq.sum.sub.4nmed$avg.pair.diff) #R2 value is 0.9627

#running analysis with this dataset
new.cols.4nmed<-lapply(seq.sum.sub.4nmed[,c("lat.mean", "lon.mean", "avg.geo.dist", "sp.lat.range", "body.size.new", "n")], scale)
seq.sum.sub.4nmed.sc<-cbind(seq.sum.sub.4nmed[,c("species", "family", "pi.median.tf", "blindness", "maternal.care")], new.cols.4nmed)
nrow(seq.sum.sub.4nmed) #91 datapoints

#running different beta regression models using the glmmTMB function
#all the variables included, have added longitudinal mean as well
mod1.4nmed<-glmmTMB(pi.median.tf ~ lat.mean + avg.geo.dist + sp.lat.range + body.size.new + blindness + maternal.care + (1|family), dispformula=~n, data=seq.sum.sub.4nmed.sc, list(family="beta", link="logit"))
summary(mod1.4nmed)

#removing the random effect
mod2.4nmed<-glmmTMB(pi.median.tf ~ lat.mean  + avg.geo.dist + sp.lat.range + body.size.new + blindness + maternal.care, dispformula=~n, data=seq.sum.sub.4nmed.sc, list(family="beta", link="logit"))
summary(mod2.4nmed) #no NANs

#removing the variable dispersion parameter
mod3.4nmed<-glmmTMB(pi.median.tf ~ lat.mean + avg.geo.dist + sp.lat.range + body.size.new + blindness + maternal.care + (1|family), data=seq.sum.sub.4nmed.sc, list(family="beta", link="logit"))
summary(mod3.4nmed)

#removing both the random effect and variable dispersion factor
mod4.4nmed<-glmmTMB(pi.median.tf ~ lat.mean + avg.geo.dist + sp.lat.range + body.size.new + blindness + maternal.care, data=seq.sum.sub.4nmed.sc, list(family="beta", link="logit"))
summary(mod4.4nmed)

#compare all the different kinds of models
bbmle::AICtab(mod1.4nmed, mod2.4nmed, mod3.4nmed, mod4.4nmed)
#model2 has the lowest AIC value

#since I am not getting r2 values from glmmTMB, and anyway the model with random effect does not turn out to be important - I am using the betareg package for model evaluation and analyzing coefficients
mod2.logit.4nmed<-betareg(pi.median.tf ~ lat.mean + avg.geo.dist + sp.lat.range + body.size.new + blindness + maternal.care | n, link="logit", data=seq.sum.sub.4nmed.sc)
summary(mod2.logit.4nmed)
mtable(mod2.logit.4nmed)

#### checking for spatial autocorrelation in model residuals ####

#converting the lat lon data into a matrix
coords<-seq.sum.sub.4nmed[,c("lon.mean", "lat.mean")]
colnames(coords)<-c("long", "lat")
rownames(coords)<-seq.sum.sub.4nmed$species
sp.dist<-geodist(x=coords, measure="haversine")/1000 #convert to km

#taking the inverse of distance to calculate weights
sp.dist.inv<-1/sp.dist
diag(sp.dist.inv)<-0

#calculate Moran's I
#note that model residuals are of different types, we are using the raw residuals here, however beta regression guidelines indicate the use of other suitable residuals - need to read about it more
Moran.I(residuals(mod2.logit.4nmed, type="response"), sp.dist.inv, alternative="two.sided")
#Moran's I = 0.06575, p value = 0.0466

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
select<-listw.select(residuals(mod2.logit.4nmed, type="response"), candidates, MEM.autocor="all", method="MIR", p.adjust=TRUE)

#looking at the best models
select$best.id
select$candidates
select$best$summary

#I think the selection procedure uses a critical value of 0.01 to check if the spatial autocorrelation is significant - so the Moran's I value is already not significant to begin with - and none of the candidate MEMs do better than the existing model, hence no model is selected

#### model diagnostics ####

#checking for multi-collinearity - "Any variable with a high VIF value (above 5 or 10) should be removed from the model. This leads to a simpler model without compromising the model accuracy, which is good." http://www.sthda.com/english/articles/39-regression-model-diagnostics/160-multicollinearity-essentials-and-vif-in-r/#:~:text=For%20a%20given%20predictor%20(p,one%20(absence%20of%20multicollinearity).
#can read the explanation of VIF here as well - https://online.stat.psu.edu/stat462/node/180/
car::vif(mod2.logit.4nmed) #all vifs are less than 5

#looking at model diagnostics plots for model evaluation - need to look into the influence of some of the points which have high Cook's distance
#this resource is helpful in understanding some of the plots - https://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/R/R5_Correlation-Regression/R5_Correlation-Regression7.html
#"A leverage point is defined as an observation that has a value of x that is far away from the mean of x."
#"An influential observation is defined as an observation that changes the slope of the line. Thus, influential points have a large influence on the fit of the model. One method to find influential points is to compare the fit of the model with and without each observation"

pdf(file=paste0("figures/model_diagnostics_4n_11Apr23.pdf"), height=11, width=8, useDingbats=FALSE)
par(mfrow=c(3,2))
plot(mod2.logit.4nmed, which=1:6)
dev.off()
#the half normal plot of residuals seem okay - the residuals largely fall within the bounds
#there seem to be some influential points based on the leverage vs predicted values plot and the Cook's distance plot

#plot the leverage plot separately with row ids
#leverage vs predicted values
diag<-data.frame(fitted=fitted(mod2.logit.4nmed), residuals=residuals(mod2.logit.4nmed), std.residuals=sqrt(abs(residuals(mod2.logit.4nmed))), leverage=gleverage(mod2.logit.4nmed), hat.diag=hatvalues(mod2.logit.4nmed))

pdf(file=paste0("figures/model_diagnostics_leverage_predicted_4n_11Apr23.pdf"), height=8, width=11, useDingbats=FALSE)
ggplot(diag, aes(x=fitted, y=leverage)) + geom_point() + theme_bw() + geom_text_repel(label=rownames(diag)) + xlab("Predicted values") + ylab("Leverage")
dev.off()

#save row names of species with high leverage
high.lev<-c(2,93,108,114,116,124)
seq.sum.sub.4nmed.sc[rownames(seq.sum.sub.4nmed.sc) %in% high.lev,]
#species with high leverage are Craterostigmus tasmanianus, Otostigmus astenus, Scolopendra cingulata, Scolopendra morsitans, Scolopendra mutilans, Allothereua maculata

#looking at cook's distance more closely. From - https://towardsdatascience.com/identifying-outliers-in-linear-regression-cooks-distance-9e212e9136a
cooksD<-cooks.distance(mod2.logit.4nmed)
influential<-cooksD[(cooksD>(3*mean(cooksD, na.rm=TRUE)))]
inf.id<-as.numeric(names(influential))
identical(influential, cooksD[names(cooksD) %in% inf.id])
seq.sum.sub.4nmed[rownames(seq.sum.sub.4nmed) %in% inf.id,] #the influential points come from Otostigmus astenus, Scolopendra dehaani, Allothereua maculata

#Cryptops parisi and Stenotaenia linearis are suspected to be a species complex - additionally removing these two species as well - saving their rownames
sp.comp<-row.names(seq.sum.sub.4nmed)[which(seq.sum.sub.4nmed$species %in% c("Cryptops parisi", "Stenotaenia linearis"))]
seq.sum.sub.4nmed[rownames(seq.sum.sub.4nmed) %in% sp.comp,]

#removing these influential points and running the regression again to see if anything changes
mod2.logit.4nmed.inf<-betareg(pi.median.tf ~ lat.mean + avg.geo.dist + sp.lat.range + body.size.new + blindness + maternal.care | n, data=seq.sum.sub.4nmed.sc[!(rownames(seq.sum.sub.4nmed.sc) %in% c(high.lev, inf.id, sp.comp)),])
summary(mod2.logit.4nmed.inf)
#latitudinal range is no longer significant, mean latitude and body size are significant, maternal care is marginally significant

#looking at the model evaluation plots again
pdf(file=paste0("figures/model_diagnostics_4n_infl_11Apr23.pdf"), height=11, width=8, useDingbats=FALSE)
par(mfrow=c(3,2))
plot(mod2.logit.4nmed.inf, which=1:6)
dev.off() #the half normal plots often show slight deviations from the expected range based on the distribution

#since some of the points seem to be influential, instead of dropping the influential points, obtaining bootstrapped coefficients

#### bootstrapping coefficients to reduce the effect of influential points on coefficient estimates ####

#calculating bootstrapped coefficients
set.seed(3453)
b.4nmed<-car::Boot(mod2.logit.4nmed, R=1000)

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
es.4nmed<-extract.coef(b.4nmed)
es.4nmed$sig<-c(1,1,0,0,1,0,0,1,1)

#plotting standardized coefficients
es.4nmed$predictor<-factor(es.4nmed$predictor, levels=c("(Intercept)", "lat.mean", "avg.geo.dist",  "sp.lat.range", "body.size.new", "blindnessnot blind", "maternal.careyes", "(phi)_(Intercept)", "(phi)_n"))

levels(es.4nmed$predictor)<-c("Intercept", "Latitude", "Geographic distance", "Latitudinal range", "Body size", "Vision:Yes", "Maternal care:Yes", "Dispersion:Intercept", "Dispersion:Slope:N")

es.4nmed$sig<-factor(es.4nmed$sig)

#not plotting the dispersion factor related parameters
es.plot<-es.4nmed[c(2:7),]
es.plot<-droplevels(es.plot)

my.col<-scale_color_manual(name="Significance", labels=c("Not significant", "Significant"), values=c("#bdbdbd", "black"))

#plot effects of standardized coefficients without exponentiating
pdf(file=paste0("figures/effect_size_unexp_no_int_4n_11Apr23.pdf"), height=8, width=11, useDingbats=FALSE)
ggplot(es.plot, aes(x=beta, y=predictor, xmin=ci.lower, xmax=ci.upper, col=sig)) + geom_vline(xintercept=0, color="red") + geom_point() + geom_errorbarh(height=.2) + my.col + theme_bw() + xlab("Standardized coefficients") + ylab("Predictor variables") + theme(text=element_text(size=20)) 
dev.off()

#### tables and figures ####

## table 1 - summary of beta regression ##
#check if the parameter names are in the same order
cbind(tidy(mod2.logit.4nmed)$term, as.character(es.4nmed$predictor))

betareg_res<-cbind(tidy(mod2.logit.4nmed), es.4nmed[,c("ci.lower", "ci.upper")])
betareg_tbl<-betareg_res %>%
mutate(term=case_match(term, '(Intercept)'~'Intercept', 'lat.mean'~'Mean latitude', 'avg.geo.dist'~'Average geographic distance', 'sp.lat.range'~'Species latitudinal range', 'body.size.new'~'Body size', 'blindnessnot blind'~'Vision: Yes', 'maternal.careyes'~'Maternal care: Yes', 'n'~'Number of sequences')) %>%
mutate(component=case_match(component, 'mean'~'Mean', 'precision'~'Precision')) %>%
mutate(across(c(3:5, 7:8), round, 3)) %>%
mutate(across(c(3:5, 7:8), format, nsmall=3)) %>%
mutate(boot.ci=paste0(ci.lower, " - ", ci.upper)) %>%
mutate(estimate=paste0(estimate, stars.pval(p.value))) 

#obtain model summary stats
footer<-glance(mod2.logit.4nmed)

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
gtsave("results/tableS5.4_model_summary_4n_11Apr23.rtf")

####