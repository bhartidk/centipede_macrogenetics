#### load required libraries and set working directory ####

#load libraries
library(dplyr)
library(ape)
library(ips)
library(pegas)
library(geodist)

#setting the working directory
setwd("D/PostDoc data/centipede_popgen/github_11Apr23_version3_final")

#### calculate sequence summaries for sequences from our database ####

#input the master sequence database
dat<-read.csv(file="data_raw/sheet4_popgen_database_analysis_clean_edit_11Apr23.csv", header=TRUE)
nrow(dat) #1429 rows

#removing rows that do not have to be included in analysis
dat<-dat %>%
filter(!is.na(Latitude)) %>%
filter(!is.na(COI)) %>%
filter(COI!="") %>%
group_by(Species) %>%
filter(n()>2)
nrow(dat) #1405 rows

#find rows which have duplicated Species, COI, Order, Longitude and Latitude
dat %>%
group_by(Species, COI, Order, Longitude, Latitude) %>% filter(n()>1)
#no duplicated rows

#creating a unique list of species
sp<-as.character(unique(dat$Species))
length(sp) #131 species

#create an empty dataframe to store the coordinate information
aln.coords<-data.frame() #sequences from native range
aln.coords.int<-data.frame() #including introductions

#create a dataframe that saves the summary statistics of the sequence alignment
seq.sum<-as.data.frame(matrix(nrow=length(sp), ncol=16))
colnames(seq.sum)<-c("species", "order", "n", "n.sites", "n.unique.sites", "lat.mean", "lon.mean", "lat.median", "seq.lat.range", "avg.geo.dist", "seq.len", "n.seg", "n.pars.inform", "avg.pair.diff", "n.hap", "hapl.div")

#create another dataframe that saves summary statistics including introduced instances of species
seq.sum.int<-seq.sum

#converting order into class 'character'
dat$Order<-as.character(dat$Order)

for(i in 1:length(sp)){
#find the accession numbers of the sequences for which location data is available 
acc<-dat %>%
filter(Introduced != 1) %>%
filter(Species %in% sp[i]) %>%
pull(COI) #only sequences from native range

acc.int<-dat %>%
filter(Species %in% sp[i]) %>%
pull(COI) #all available sequences, including instances of introductions

#find the order
sp.order<-dat %>%
filter(Species %in% sp[i]) %>%
pull(Order) %>%
unique()

#reading in the aligned sequences for sp[h]
file.name<-gsub(" ", "_", sp[i])
aln<-read.dna(paste0("data/seq_aln_final/aln_", file.name, ".fasta"), format="fasta")

#take a subset of aln.sub based on what accession numbers have location data
ids<-which(labels(aln) %in% acc)
aln.sub<-aln[ids,] #native range

ids.int<-which(labels(aln) %in% acc.int)
aln.sub.int<-aln[ids.int,] #including introductions

#get the accession numbers from the alignment - there were some accession numbers I removed because the sequence was either short or not aligning properly
acc.aln<-labels(aln.sub)
acc.aln.int<-labels(aln.sub.int)

#print species names where the accession numbers between acc.aln are different from acc.aln.int
if(!identical(acc.aln, acc.aln.int)){
print(paste0(sp[i], " has ", length(acc.aln), " sequences from its native range and ", length(setdiff(acc.aln.int, acc.aln)), " total number of  sequences, including introductions."))
}

##get the sequence statistics from the alignment
#number of sequences
n<-nrow(aln.sub)
n.int<-nrow(aln.sub.int)

#sequence length
seq.len<-ncol(aln.sub)
seq.len.int<-ncol(aln.sub.int)

n.seg<-length(seg.sites(aln.sub, strict=FALSE))
n.seg.int<-length(seg.sites(aln.sub.int, strict=FALSE))
#in the code above, "if strict = FALSE (the default), the following rule is used to determine if a site is polymorphic or not in the presence of ambiguous bases: `A' and `R' are not interpreted as different, `A' and `Y' are interpreted as different, and `N' and any other base (ambiguous or not) are interpreted as not different. If strict = TRUE, all letters are considered different."

n.pars.inform<-pis(aln.sub, what="absolute", use.ambiguities=FALSE)
n.pars.inform.int<-pis(aln.sub.int, what="absolute", use.ambiguities=FALSE)
#the use.ambiguities argument is not active in the function yet
#from the MEGA manual - "a site is parsimony-informative if it contains at least two types of nucleotides (or amino acids), and at least two of them occur with a minimum frequency of two."

avg.pair.diff<-pegas::nuc.div(aln.sub, variance=FALSE, pairwise.deletion=FALSE)
avg.pair.diff.int<-pegas::nuc.div(aln.sub.int, variance=FALSE, pairwise.deletion=FALSE)
#"The first method uses the sum of the number of differences between pairs of sequences divided by the number of comparisons (i.e. n(n - 1)/2, where n is the number of sequences)." and the reference cited is Nei, 1987
#usually nucleotide diversity (Pi) is defined as the average number of nucleotide differences per site between two sequences
#pairwise.deletion - "a logical indicating whether to delete the sites with missing data in a pairwise way. The  default  is  to  delete  the  sites  with  at  least  one  missing  data  for all sequences."

#detect haplotypes from the alignment subset
hap<-haplotype(aln.sub, labels=c(1:nrow(haplotype(aln.sub))))
hap.int<-haplotype(aln.sub.int, labels=c(1:nrow(haplotype(aln.sub.int))))
#look at the warnings a little more closely later
#strict=FALSE - "a logical value; if TRUE, ambiguities and gaps in the sequences are ignored"
#How ambiguities are handled is explained here - https://www.mail-archive.com/r-sig-phylo@r-project.org/msg05541.html

#number of haplotypes
n.hap<-nrow(hap)
n.hap.int<-nrow(hap.int)

if(n>0){
hapl.div<-hap.div(aln.sub)
}else{
hapl.div<-NA
}
hapl.div.int<-hap.div(aln.sub.int)
#"currently, only Nei and Tajima’s (1981) method is available". the method argument within the function is therefore not used.
#it can be thought of as the probability that the same haplotype is not picked twice and is given by (N/N-1)*(1-summation((pi)^2), where pi is the frequency of haplotype i in the population - ranges from 0-1, larger values indicate higher haplotype diversity

#finding the locations corresponding to the accession numbers in the alignment
#native range
coords<-dat %>%
filter(Introduced != 1) %>%
filter(Species %in% sp[i]) %>%
filter(COI %in% acc.aln) %>%
rename("x"="Longitude", "y"="Latitude") %>%
select(COI, Species, Voucher_no, Order, Family, x, y, Introduced) %>%
mutate(x=as.numeric(x), y=as.numeric(y)) %>%
as.data.frame()

aln.coords<-rbind(aln.coords, coords)
  
coords.sum<-coords %>%  
summarise(lat.mean=mean(y), lat.median=median(y), lon.mean=mean(x), seq.lat.range=max(y)-min(y), all.coords=nrow(coords), unq.coords=n_distinct(coords[,c("x", "y")]))

#including introductions
coords.int<-dat %>%
filter(Species %in% sp[i]) %>%
filter(COI %in% acc.aln.int) %>%
rename("x"="Longitude", "y"="Latitude") %>%
select(COI, Species, Voucher_no, Order, Family, x, y, Introduced) %>%
mutate(x=as.numeric(x), y=as.numeric(y)) %>%
as.data.frame()

aln.coords.int<-rbind(aln.coords.int, coords.int)

if(n!=nrow(coords)){
print(paste0("The number of sequences is not equal to the number of native sites for ", sp[i], "."))
}

if(n.int!=nrow(coords.int)){
print(paste0("The number of sequences is not equal to the number of total sites for ", sp[i], " ."))
}

coords.sum.int<-coords.int %>%  
summarise(lat.mean=mean(y), lat.median=median(y), lon.mean=mean(x), seq.lat.range=max(y)-min(y), all.coords=nrow(coords.int), unq.coords=n_distinct(coords.int[,c("x","y")]))

if(n>0){
#finding the pairwise distances between coordinates
#native range
geo.pdist<-geodist(x=coords[,c("x", "y")], measure="haversine") #this gives a matrix
avg.geo.pdist<-mean(geo.pdist[lower.tri(geo.pdist)])/1000
}else{
avg.geo.dist<-NA
}

#including introductions
geo.pdist.int<-geodist(x=coords.int[,c("x", "y")], measure="haversine") #this gives a matrix
avg.geo.pdist.int<-mean(geo.pdist.int[lower.tri(geo.pdist.int)])/1000

#adding the sequence statistics to the dataframe
seq.sum[i,]<-c(sp[i], sp.order, n, coords.sum$all.coords, coords.sum$unq.coords, coords.sum$lat.mean, coords.sum$lon.mean, coords.sum$lat.median, round(coords.sum$seq.lat.range, 4), round(avg.geo.pdist, 4), seq.len, n.seg, n.pars.inform, round(avg.pair.diff,4), n.hap, round(hapl.div,4))

#adding the sequence statistics to the dataframe
seq.sum.int[i,]<-c(sp[i], sp.order, n.int, coords.sum.int$all.coords, coords.sum.int$unq.coords, coords.sum.int$lat.mean, coords.sum.int$lon.mean, coords.sum.int$lat.median, round(coords.sum.int$seq.lat.range, 4), round(avg.geo.pdist.int, 4), seq.len.int, n.seg.int, n.pars.inform.int, round(avg.pair.diff.int,4), n.hap.int, round(hapl.div.int,4))
}
rm(sp, sp.order, n, n.int, coords.sum, coords.sum.int, avg.geo.pdist, avg.geo.pdist.int, seq.len, seq.len.int, n.seg, n.seg.int, n.pars.inform, n.pars.inform.int, avg.pair.diff, avg.pair.diff.int, n.hap, n.hap.int, hapl.div, hapl.div.int)

#change sequence summary columns to numeric
str(seq.sum)
seq.sum[,3:16]<-lapply(seq.sum[,3:16], as.numeric)
seq.sum.int[,3:16]<-lapply(seq.sum.int[,3:16], as.numeric)

#remove NA values from seq.sum
seq.sum<-seq.sum %>%
na.omit %>%
as.data.frame
nrow(seq.sum) #125 species

seq.sum.int<-seq.sum.int %>%
na.omit %>%
as.data.frame
nrow(seq.sum.int) #130 species

#remove NAs from aln.coords
aln.coords<-aln.coords %>%
na.omit(Longitude, Latitude) %>%
group_by(Species)
nrow(aln.coords) #1232 #Lamyctes africanus has a single native occurrence which is captured in the coordinates, but takes NA value for sequence summary since it is a single sequence, same for Scolopendra calcarata and therefore the difference in the number of sequences and coordinates
sum(seq.sum$n) #1230

aln.coords.int<-aln.coords.int %>%
na.omit(Longitude, Latitude) %>%
group_by(Species) %>%
as.data.frame
nrow(aln.coords.int) #1367 - difference due to Scolopendra calcarata having just one sequence
sum(seq.sum.int$n) #1366

#### calculate sequence summaries for the additional species sequences from phylogatr database ####

#load the dataframe for missing species from phylogatr
dat.pg.missing<-read.csv("data/phylogatr_missing_species_11Apr23.csv", header=TRUE)

#save unique species names; remove Anopsobius neozelanicus because it may not represent a single species based on existing literature
missing.sp<-dat.pg.missing %>%
filter(species!="Anopsobius neozelanicus") %>%
distinct(species) %>%
pull()

#we've stored the names of the 11 missing species in the variable missing.sp and the new aligments are stored in the following address
file.aln<-"data/seq_aln_final/phylogatR_species_only/"

#create an empty dataframe to store the coordinate information
aln.coords.pg<-data.frame()

#create a dataframe that saves the summary statistics of the sequence alignment
seq.sum.pg<-as.data.frame(matrix(nrow=length(missing.sp), ncol=16))
colnames(seq.sum.pg)<-c("species", "order", "n", "n.sites", "n.unique.sites", "lat.mean", "lon.mean", "lat.median", "seq.lat.range", "avg.geo.dist", "seq.len", "n.seg", "n.pars.inform", "avg.pair.diff", "n.hap", "hapl.div")

#looping across missing.sp alignments to save sequence summary statistics
for(i in 1:length(missing.sp)){
#creating the file name for the alignment pointer from the species name
file.sp<-paste0(gsub(" ", "-", missing.sp[i]), "-COI.afa")

#find the order
sp.order<-dat.pg.missing %>%
filter(species %in% missing.sp[i]) %>%
pull(order) %>%
unique()

#reading in the aligned sequences for sp[h]
aln.sub<-read.dna(paste0(file.aln, file.sp), format="fasta")

#get the labels of the sequences in the alignment, this would correspond to the phylogatr id
phylogatr.id<-labels(aln.sub)

#get the sequence statistics from the alignment
n<-length(phylogatr.id)
seq.len<-ncol(aln.sub)
n.seg<-length(ape::seg.sites(aln.sub, strict=FALSE))
n.pars.inform<-ips::pis(aln.sub, what="absolute", use.ambiguities=FALSE)
avg.pair.diff<-pegas::nuc.div(aln.sub, variance=FALSE, pairwise.deletion=FALSE)
hap<-pegas::haplotype(aln.sub, labels=c(1:nrow(haplotype(aln.sub))))
n.hap<-nrow(hap)
hapl.div<-hap.div(aln.sub)

#finding the locations corresponding to the accession numbers in the alignment
coords<-dat.pg.missing %>%
filter(species %in% missing.sp[i]) %>%
filter(phylogatr_id %in% phylogatr.id) %>%
rename("COI"="accession", "Species"="species", "Voucher_no"="phylogatr_id", "Order"="order", "Family"="family", "x"="longitude", "y"="latitude") %>%
select(COI, Species, Voucher_no , Order, Family, x, y) %>%
mutate(x=as.numeric(x), y=as.numeric(y)) %>%
as.data.frame()
coords$Introduced<-rep(0, times=nrow(coords))

aln.coords.pg<-rbind(aln.coords.pg, coords)

coords.sum<-coords %>%  
summarise(lat.mean=mean(y), lat.median=median(y), lon.mean=mean(x), seq.lat.range=max(y)-min(y), all.coords=nrow(coords), unq.coords=n_distinct(coords[,c("x", "y")]))

if(n!=nrow(coords)){
print(paste0("The number of sequences is not equal to the number of native sites for ", missing.sp[i], " ."))
}

#finding the pairwise distances between coordinates
geo.pdist<-geodist(x=coords, measure="haversine") #this gives a matrix
avg.geo.pdist<-mean(geo.pdist[lower.tri(geo.pdist)]/1000)

#adding the sequence statistics to the dataframe
seq.sum.pg[i,]<-c(missing.sp[i], sp.order, n, coords.sum$all.coords, coords.sum$unq.coords, coords.sum$lat.mean, coords.sum$lon.mean, coords.sum$lat.median, round(coords.sum$seq.lat.range, 4), round(avg.geo.pdist, 4), seq.len, n.seg, n.pars.inform, round(avg.pair.diff,4), n.hap, round(hapl.div,4))
}

rm(sp.order, n, coords.sum, avg.geo.pdist, seq.len, n.seg, n.pars.inform, avg.pair.diff, n.hap, hapl.div)

#change sequence summary columns to numeric
seq.sum.pg[,3:16]<-lapply(seq.sum.pg[,3:16], as.numeric)

#reorder species names
seq.sum.pg<-seq.sum.pg[order(seq.sum.pg$order, seq.sum.pg$species),]

#write the sequence summary to disk
write.csv(seq.sum.pg, file=paste0("data/sequence_summary_missing_species_phylogatr_11Apr23.csv"), row.names=FALSE)

#### combine sequence summaries from the two databases and add life history information ####

#combining sequence summaries from the two databases
seq.sum.final<-rbind(seq.sum, seq.sum.pg)
seq.sum.final<-seq.sum.final[order(seq.sum.final$order, seq.sum.final$species),]
nrow(seq.sum.final) #135 species

seq.sum.int.final<-rbind(seq.sum.int, seq.sum.pg)
seq.sum.int.final<-seq.sum.int.final[order(seq.sum.int.final$order, seq.sum.int.final$species),] 
nrow(seq.sum.int.final) #140 species

#input the the life history data
lh<-read.csv("data_raw/centipedes_life_history_11Apr23.csv", header=TRUE)
str(lh)

#add species traits to the sequence summary information
seq.sum.lh<-left_join(seq.sum.final, lh[,-c(2,7,8)], by=c("species"))
nrow(seq.sum.lh) #135 species

seq.sum.int.lh<-left_join(seq.sum.int.final, lh[,-c(2,9,10)], by=c("species"))
nrow(seq.sum.int.lh) #140 species

#reorder columns
colnames(seq.sum.int.lh)
seq.sum.lh<-seq.sum.lh[,c(1:2, 17, 3:16, 18:22)]
seq.sum.int.lh<-seq.sum.int.lh[,c(1:2, 17, 3:16, 18:22)]

#remove information for species where monophyly is in question
seq.sum.lh<-seq.sum.lh %>%
filter(species!="Scolopendra subspinipes") %>%
filter(species!="Cormocephalus sp.")  %>%
as.data.frame

#input for genetic diversity arthropods
write.csv(seq.sum.lh, "data/genetic_diversity_arthropods_centipedes_9May23.csv", row.names=FALSE)

#clean up seq.sum.lh for statistical analysis by filtering samples with less than three sequence representatives and missing life history information
seq.sum.lh.clean<-seq.sum.lh %>%
filter(n>2) %>%
filter(!is.na(body.size.new)) %>%
as.data.frame
nrow(seq.sum.lh.clean) #128 species
str(seq.sum.lh.clean)

#clean up seq.sum.int.lh for statistical analysis - remove Scolopendra subspinipes and Cormocephalus sp., filter samples with less than three sequence representatives
seq.sum.int.lh.clean<-seq.sum.int.lh %>%
filter(species!="Scolopendra subspinipes") %>%
filter(species!="Cormocephalus sp.") %>%
filter(n>2) %>%
filter(!is.na(body.size.new)) %>%
as.data.frame
nrow(seq.sum.int.lh.clean) #134 species

#find the species absent in the native range dataset
setdiff(seq.sum.int.lh.clean$species, seq.sum.lh.clean$species) # six species: Geophilus proximus, Pachymerium ferrugineum, Lamyctes africanus, Lamyctes coeculus, Lamyctes emarginatus, Scutigera coleoptrata

#write the input file to disk
write.csv(seq.sum.lh.clean, "data/input_betareg_no_introductions_11Apr23.csv", row.names=FALSE)
write.csv(seq.sum.int.lh.clean, "data/input_betareg_11Apr23.csv", row.names=FALSE)

#save the sequence coordinates to disk
aln.coords.final<-rbind(aln.coords, aln.coords.pg)
aln.coords.int.final<-rbind(aln.coords.int, aln.coords.pg)

#write these coordinates to disk
write.csv(aln.coords.final, file="data/COI_final_locations_no_introductions_11Apr23.csv", row.names=FALSE)
write.csv(aln.coords.int.final, file="data/COI_final_locations_11Apr23.csv", row.names=FALSE)

####
