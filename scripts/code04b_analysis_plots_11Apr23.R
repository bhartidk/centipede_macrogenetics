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

#### prepare input file for analysis ####

#input sequence summary
seq.sum<-read.csv(file="data/input_betareg_11Apr23.csv", header=TRUE, na.strings=c("","NA"))
colnames(seq.sum)

#number of species with complete information
nrow(seq.sum[complete.cases(seq.sum),]) #134

#create a column with latitudinal range
seq.sum$sp.lat.range<-seq.sum$max.lat-seq.sum$min.lat

#taking a subset of seq.sum which have non-NA data
seq.sum.sub<-seq.sum[,c("species", "family", "order", "avg.pair.diff", "lat.mean", "lon.mean", "avg.geo.dist", "sp.lat.range", "body.size.new", "blindness", "maternal.care", "n", "n.unique.sites", "seq.len", "hapl.div")]

#check number of species where average pairwise difference equals 0
seq.sum.sub %>%
filter(avg.pair.diff==0) %>%
pull(species)

#we have 10 species where the average pairwise difference value is equal to 0, which cannot be modeled by beta regression - using a modification function to deal with this - code from Douma et al., 2019
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
colnames(seq.sum.sub2)

res.pear<-cor(seq.sum.sub2[,c(6,8,9,10,11)], method="pearson")
pdf(file=paste0("figures/pearson_correlation_predictors_int_11Apr23.pdf"), height=8, width=11, useDingbats=FALSE)
corrplot(res.pear, type="upper", order="hclust", tl.col="black", tl.srt=45, addCoef.col='black')
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

#### model diagnostics ####

#checking for multi-collinearity - "Any variable with a high VIF value (above 5 or 10) should be removed from the model. This leads to a simpler model without compromising the model accuracy, which is good." http://www.sthda.com/english/articles/39-regression-model-diagnostics/160-multicollinearity-essentials-and-vif-in-r/#:~:text=For%20a%20given%20predictor%20(p,one%20(absence%20of%20multicollinearity).
#can read the explanation of VIF here as well - https://online.stat.psu.edu/stat462/node/180/
car::vif(mod2.logit) #all VIF values are less than 5

#obtaining a table of coefficient estimates
mtable(mod2.logit) #the value in brackets is the standard error associated with the parameter estimate

#looking at model diagnostics plots for model evaluation - need to look into the influence of some of the points which have high Cook's distance
#this resource is helpful in understanding some of the plots - https://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/R/R5_Correlation-Regression/R5_Correlation-Regression7.html
#"A leverage point is defined as an observation that has a value of x that is far away from the mean of x."
#"An influential observation is defined as an observation that changes the slope of the line. Thus, influential points have a large influence on the fit of the model. One method to find influential points is to compare the fit of the model with and without each observation"

pdf(file=paste0("figures/model_diagnostics_int_11Apr23.pdf"), height=11, width=8, useDingbats=FALSE)
par(mfrow=c(3,2))
plot(mod2.logit, which=1:6)
dev.off()
#the half normal plot of residuals show some deviations - but they fall within or close to the boundss
#there seem to be some influential points based on the leverage vs predicted values plot and the Cook's distance plot

#plot the leverage plot separately with row ids
#leverage vs predicted values
diag<-data.frame(fitted=fitted(mod2.logit), residuals=residuals(mod2.logit), std.residuals=sqrt(abs(residuals(mod2.logit))), leverage=gleverage(mod2.logit), hat.diag=hatvalues(mod2.logit))

pdf(file=paste0("figures/model_diagnostics_leverage_predicted_int_11Apr23.pdf"), height=8, width=11, useDingbats=FALSE)
ggplot(diag, aes(x=fitted, y=leverage)) + geom_point() + theme_bw() + geom_text_repel(label=rownames(diag)) + xlab("Predicted values") + ylab("Leverage")
dev.off()

#save row ids of species with high leverage
high.lev<-c(1,2,15,32,80,98)
seq.sum.sub2[high.lev,]
#species with high leverage are Craterostigmus crabilli, C tasmanianus, Schendyla nemorensis, Lamyctes coeculus, Cryptops hortensis, Otostigmus astenus 

#looking at cook's distance more closely. From - https://towardsdatascience.com/identifying-outliers-in-linear-regression-cooks-distance-9e212e9136a
cooksD<-cooks.distance(mod2.logit)
influential<-cooksD[(cooksD>(3*mean(cooksD, na.rm=TRUE)))]
inf.id<-as.numeric(names(influential))
identical(influential, cooksD[names(cooksD) %in% inf.id])
seq.sum.sub2[rownames(seq.sum.sub2) %in% inf.id,] #the influential points come from Lamyctes coeculus, Cryptops hortensis, Cryptops japonicus, Scolopendra dehaani, Scolopendra morsitans

#Cryptops parisi and Stenotaenia linearis are suspected to be a species complex - additionally removing these two species as well
sp.comp<-which(seq.sum.sub2$species %in% c("Cryptops parisi", "Stenotaenia linearis"))
seq.sum.sub2[sp.comp,]

#removing these influential points and running the regression again to see if anything changes
mod2.logit.inf<-betareg(avg.pair.diff.tf ~ lat.mean + avg.geo.dist + sp.lat.range + body.size.new + blindness + maternal.care | n, data=seq.sum.sub2[-unique(c(inf.id, sp.comp)),])
summary(mod2.logit.inf)
#the relationships are much stronger and the R-squared value jumps by 5% when influential points are removed - the inferences remain the same.

#looking at the model evaluation plots again
pdf(file=paste0("figures/model_diagnostics_infl_int_11Apr23.pdf"), height=11, width=8, useDingbats=FALSE)
par(mfrow=c(3,2))
plot(mod2.logit.inf, which=1:6)
dev.off()

#since some of the points seem to be influential, instead of dropping the influential points, obtaining bootstrapped coefficients

#### bootstrapping coefficients to reduce the effect of influential points on coefficient estimates ####

#using the Boot package to obtain bootstrapped estimates from here - https://socialsciences.mcmaster.ca/jfox/Books/Companion-2E/appendix/Appendix-Bootstrapping.pdf. Some notes for the type of the confidence intervals is here - https://stats.stackexchange.com/questions/355781/is-it-true-that-the-percentile-bootstrap-should-never-be-used#:~:text=Yes%2C%20if%20we%20can%2C%20we,and%20think%20we%20discovered%20America.
set.seed(3453)
b<-car::Boot(mod2.logit, R=1000)

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
es$sig<-c(1,1,0,0,1,1,1,1,1)

#to interpret the betas we need to exponentiate them - since the regression is presented as a logit transformation
es$exp.beta<-exp(es$beta)
es$exp.ci.lower<-exp(es$ci.lower)
es$exp.ci.upper<-exp(es$ci.upper)

#plotting standardized coefficients
es$predictor<-factor(es$predictor, levels=c("(Intercept)", "lat.mean", "avg.geo.dist",  "sp.lat.range", "body.size.new", "blindnessnot blind", "maternal.careyes", "(phi)_(Intercept)", "(phi)_n"))

levels(es$predictor)<-c("Intercept", "Latitude", "Geographic distance", "Latitudinal range", "Body size", "Vision:Yes", "Maternal care:Yes", "Dispersion:Intercept", "Dispersion:Slope:N")

es$sig<-factor(es$sig)

#not plotting the dispersion factor related parameters
es.plot<-es[-c(1,8:9),]
es.plot<-droplevels(es.plot)

my.col<-scale_color_manual(name="Significance", labels=c("Not significant", "Significant"), values=c("#bdbdbd", "black"))

#plot effects of standardized coefficients without exponentiating
pdf(file=paste0("figures/effect_size_unexp_int_11Apr23.pdf"), height=8, width=11, useDingbats=FALSE)
ggplot(es.plot, aes(x=beta, y=predictor, xmin=ci.lower, xmax=ci.upper, col=sig)) + geom_vline(xintercept=0, color="red") + geom_point() + geom_errorbarh(height=.2) + my.col + theme_bw() + xlab("Standardized coefficients") + ylab("Predictor variables") + theme(text=element_text(size=20))
dev.off()

#A 1SD change in body size leads to a relative change of exp(-3.36)=0.035 units change in E(Proportion)/(1−E(Proportion)). Need to interpret this in some sensible way - what we are interested in is the Proportion, but what we are having to interpret is some transformation of this.

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

#add the species and family information to model residuals
#obtain model residuals from the beta regression model that was run
res<-residuals(mod2.logit, type="response")
names(res)<-as.character(seq.sum.sub2$family)
res.df<-data.frame(res=residuals(mod2.logit, type="response"), family=as.character(seq.sum.sub2$family))
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

#summary statistics
nrow(seq.sum.sub) #134
max(seq.sum.sub$avg.geo.dist) #9585.598
max(seq.sum.sub$sp.lat.range) #103.2196

#reproductive strategy - number of species showing maternal care
seq.sum.sub %>%
group_by(maternal.care) %>%
summarise(n()) #85 species with maternal care of 128 total species

#vision - number of species with vision
seq.sum.sub %>%
group_by(blindness) %>%
summarise(n()) #101 species with vision of 134 total species

#loading information of locations
acc.sp<-read.csv("data/COI_final_locations_11Apr23.csv", header=TRUE)
#statistics from all sequence locations
acc.final<-acc.sp %>%
filter(Species %in% seq.sum.sub$species)
nrow(acc.final) #1245 locations
sum(seq.sum.sub$n) #1245 sequences
#total number of unique sequence locations
nrow(unique(acc.final[,c("x", "y")])) #834 unique geographic locations

#mean and range genetic diversity of centipedes
seq.sum.sub %>%
summarise(mean(avg.pair.diff), min(avg.pair.diff), max(avg.pair.diff))

## table 1 - summary of beta regression ##
#check if the parameter names are in the same order
cbind(tidy(mod2.logit)$term, as.character(es$predictor))

betareg_res<-cbind(tidy(mod2.logit), es[,c("ci.lower", "ci.upper")])
betareg_tbl<-betareg_res %>%
mutate(term=case_match(term, '(Intercept)'~'Intercept', 'lat.mean'~'Mean latitude', 'avg.geo.dist'~'Average geographic distance', 'sp.lat.range'~'Species latitudinal range', 'body.size.new'~'Body size', 'blindnessnot blind'~'Vision: Yes', 'maternal.careyes'~'Maternal care: Yes', 'n'~'Number of sequences')) %>%
mutate(component=case_match(component, 'mean'~'Mean', 'precision'~'Precision')) %>%
mutate(across(c(3:5, 7:8), round, 3)) %>%
mutate(across(c(3:5, 7:8), format, nsmall=3)) %>%
mutate(boot.ci=paste0(ci.lower, " - ", ci.upper)) %>%
mutate(estimate=paste0(estimate, stars.pval(p.value))) 

#obtain model summary stats
footer<-glance(mod2.logit)

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
gtsave("results/tableS5.2_model_summary_int_11Apr23.rtf")

####