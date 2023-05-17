#### load required libraries and set working directory ####

#load libraries
library(dplyr)
library(ape)
library(rgdal)

#setting the working directory
setwd("D/PostDoc data/centipede_popgen/github_11Apr23_version3_final")

#use path to the phylogatr data
file.add<-"data_raw/phylogatr_results_job_19833727_Chilopoda/phylogatr-results/"

#### combining phylogatr species summary information with sequence information ####

#reading in the phylogatr datafile that has species information
dat.pg<-read.table(file=paste0(file.add, "genes.txt"), sep='\t', header=TRUE, na.strings=c("","NA"))
str(dat.pg)

#only keep rows corresponding to COI data - COI is used as the suffix in the column named 'gene'
dat.pg<-dat.pg[grep("COI", dat.pg$gene),]
head(dat.pg)
length(unique(dat.pg$species)) #73
nrow(dat.pg) #73

#create a new empty dataframe and assign relevant column names corresponding to the occurrence.txt file in the phylogatr output for each species
occ.eg<-read.table(file=paste0(file.add, "Chilopoda/Lithobiomorpha/Lithobiidae/Lithobius-latro/occurrences.txt"), sep='\t', header=TRUE, na.strings=c("","NA"))
dat.occ<-cbind(dat.pg[0,], occ.eg[0,])
str(dat.occ)

#running a loop to aggregate the occurrence locations for all species together
for(i in 1:nrow(dat.pg)){
#read the occurrence data for the first species in dat.pg
occ<-read.table(file=paste0(file.add, dat.pg$dir[i],"/occurrences.txt"), sep='\t', header=TRUE, na.strings=c("","NA"))
occ$source_id<-as.character(occ$source_id)

#read the fasta file of COI sequences
seq<-read.dna(file=paste0(file.add, dat.pg$dir[i],"/", dat.pg$gene[i], ".fa"), format="fasta")

#read COI sequence labels
seq.lab<-labels(seq)

#only keep the accession numbers in occ which correspond to the COI labels above
occ<-occ %>%
filter(phylogatr_id %in% seq.lab)

#check if the labels match
if(!identical(seq.lab, occ$phylogatr_id)){
print(paste0("Incorrect phylogatr ids subset for "), dat.pg$species[i])
}

#combine data
dat.merge<-cbind(dat.pg[i,], occ)

#copy the occurrence information
dat.occ<-bind_rows(dat.occ, dat.merge)

#clear the variables that will be reused in the loop
rm(seq, seq.lab, occ, dat.merge)
}

#check if all the species are covered
length(unique(dat.occ$species)) #73

#look at the summary of this newly populated dataframe
summary(dat.occ)

#write the data to disk
write.csv(dat.occ, file=paste0("data/phylogatr_genes_occurrence_11Apr23.csv"), row.names=FALSE)

#### comparing our database with the phylogatr database ####

#reading in our database
dat<-read.csv(file="data_raw/sheet4_popgen_database_analysis_clean_edit_11Apr23.csv", header=TRUE)

#keep taxa identified to the species level, remove entries where coordinates are missing
dat<-dat %>%
filter(Species!="Geophilomorpha sp.") %>%
filter(!is.na(Latitude))  %>%
group_by(Species) %>%
filter(n()>2)
nrow(dat) #1420 matching the species in the final input file

length(unique(dat$Species)) #130 species

#what are the total number of distinct accession numbers in each database (phylogatR vs our compilation)
length(unique(dat$COI)) #1401
length(unique(dat.occ$accession)) #620

#what are the total number of distinct species in each database
length(unique(dat$Species)) #130
length(unique(dat.occ$species)) #73

#filtering accession numbers that are common and exclusive between the two datasets
pg.unq<-setdiff(dat.occ$accession, dat$COI) #unique to phylogatr
length(pg.unq) #71 accessions - many of the accessions in phylogatr are blanks as the sequences have not been submitted to GenBank and hence the number of exclusive sequences in phylogatr is a larger number
dat.unq<-setdiff(dat$COI, dat.occ$accession) #unique to our database
length(dat.unq) #852 accessions
common<-intersect(dat.occ$accession, dat$COI) #common to phylogatr and our database
length(common) #549 accessions
#common and dat.unq add up to 1401 sequences

#check if the species assignments are correct for the common accession numbers
dat.occ.common<-dat.occ %>%
filter(accession %in% common) %>%
select(accession, species) %>%
as.data.frame()

dat.common<-dat %>%
filter(COI %in% common) %>%
select(COI, Species) %>%
as.data.frame()

dat.common.merged<-left_join(dat.occ.common, dat.common, by=c('accession'='COI'))

ids<-which(dat.common.merged$species != dat.common.merged$Species)
dat.common.merged[ids,]

write.csv(dat.common.merged, "data/common_accession_species_comparison_11Apr23.csv", row.names=FALSE)
#AB612267, AB612268, AB612889 are assigned as Scolopendra mutilans in phylogatr but are called Scolopendra subspinipes by us - since there is some confusion in the assignment of species names, removing Scolopendra subspinipes from analysis
#AB614400, AB672740, AB672741 are assigned as Scolopendra subspinipes in phylogatr but are called Scolopendra japonica by us - these come from Kang et al., 2017 and one of the sequences here is very similar to a sequence from Siriwut et al., 2015b, which gives support to the identification, and hence retaining our assignment of Scolopendra japonica

#saving the database of unique accession numbers to disk
dat.unique<-dat %>%
filter(COI %in% dat.unq) %>%
as.data.frame()
write.csv(dat.unique, "data/our_data_unique_accessions_11Apr23.csv", row.names=FALSE)

#pulling out names of species corresponding to unique accession numbers in phylogatr
pg.sp<-dat.occ %>%
filter(accession %in% pg.unq) %>%
distinct(species) %>%
pull()
length(pg.sp) #34 species

#number of species that are unique to phylogatr from the above list
missing.sp<-setdiff(pg.sp, unique(dat$Species))
length(missing.sp) #11 species are only found in phylogatr
#These are Anopsobius neozelanicus, Clinopodes carinthiacus,  Cryptops hortensis_rucneri, Dicellophilus carniolensis, Hemiscolopendra marginata, Lithobius subtilis, Lithobius variegatus, Schendyla armata, Schendyla tyrolensis, Scutigerina weberi and Strigamia chionophila

#save data corresponding to the missing species to disk
dat.occ.missing<-dat.occ %>%
filter(species %in% missing.sp) %>%
as.data.frame

#write data for missing species to disk
write.csv(dat.occ.missing, "data/phylogatr_missing_species_11Apr23.csv", row.names=FALSE)

#### look at the common species, where there are additional accession numbers in phylogatr that are missing from our data ####

#subset accession numbers from phylogatr, which are missing in the common species 
missing.acc.sp<-intersect(unique(dat$Species), pg.sp)
length(missing.acc.sp) #23 species have unique accessions in phylogatr that are missing from the common species between these two datasets 
dat.missing<-dat.occ %>%
filter(accession %in% pg.unq) %>%
filter(species %in% missing.acc.sp)

#write the above dataframe to disk
write.csv(dat.missing, file=paste0("data/phylogatr_common_species_missing_accessions_11Apr23.csv"), row.names=FALSE) #missing accessions from our data, which are present in phylogatR for the common species

#in the sheet above, entries marked as being native or introduced based on distribution records for each species and saved as phylogatr_common_species_missing_accessions_edit_11Apr23.csv

#input this edited dataframe as dat.missing and do the distribution plots again to make sure we've marked the entries correctly
dat.missing<-read.csv("data/phylogatr_common_species_missing_accessions_edit_11Apr23.csv", header=TRUE)

#filter data to remove entries which are flagged, remove Scolopendra subspinipes based on difference species name in phylogatr and GenBank record, remove Pachymerium ferrugineum because it may be introduced in many places and hence is hard to determine accurate latitudinal range limits
rm.sp<-c("Scolopendra subspinipes", "Pachymerium ferrugineum")
dat.missing<-dat.missing %>% 
filter(introduced==0) %>%
filter(is.na(flag)) %>%
filter(!(species %in% rm.sp)) %>%
filter(!is.na(accession))

#write the above dataframe to disk
write.csv(dat.missing, file=paste0("data/phylogatr_common_species_missing_accessions_final_11Apr23.csv"), row.names=FALSE) #missing accessions from our data, which are present in phylogatR for the common species

#find the summary of this dataset
dat.missing.sum<-dat.missing %>%
group_by(species) %>%
summarise(n=n_distinct(phylogatr_id), min.lat=min(latitude), max.lat=max(latitude), dir=unique(dir), gene=unique(gene)) %>%
as.data.frame()
nrow(dat.missing.sum) #9 species
sum(dat.missing.sum$n) #24 sequences, most number for Cryptops parisi (7)

#### I aligned the new phylogatr sequences for the 10 common species with the existing alignments for these species from our database. Of the 10 species, 9 species come from a single barcoding study from Slovenia. Based on the realignment and BLAST searches of these sequences, looks like around 5-6 species from phylogatr may not have been correctly identified in the study - details below. Since the number of additional sequences from phylogatr is small and the species identification may be uncertain, not adding these to the existing alignments. However, since for additional species from phylogatr, the sequences are the only ones available, retaining only the additional species in the analysis, even when they are from BOLD and haven't yet been assigned GenBank ids. ####

####