#### load required libraries and set working directory ####

#load libraries
library(dplyr)
library(rentrez)
library(ape)
library(Biostrings)
library(muscle)

#setting the working directory
setwd("D/PostDoc data/centipede_popgen/github_11Apr23_version3_final")

####load and clean data, query accession numbers and save to disk ####

#loading the clean sequence database - it was dated 23Apr21 in code5
dat<-read.csv("data_raw/sheet4_popgen_database_analysis_clean_11Apr23.csv", header=TRUE)
nrow(dat) #1438 rows

#only retain rows that have non-NA GPS coordinates
dat<-dat %>%
na.omit(Latitude)
nrow(dat) #1430 rows

#check if there are any duplicated COI accession numbers
dat %>% 
filter(!COI=="") %>% #some of the sequences from Jahnavi's work have not been submitted to GenBank and therefore do not have accession numbers
group_by(COI) %>%
filter(n()>1) %>%
as.data.frame()
#AB612892 has been assigned to both Scolopendra mutilans and S subspinipes. Since mutilans was considered a subspecies of subspinipes, keeping the sub-species identity as the species name, and removing the assignment to S subspinipes

#removing the row where AB612892 is assigned to Scolopendra subspinipes
dat<-dat %>%
filter(!(COI=="AB612892" & Species=="Scolopendra subspinipes"))
nrow(dat) #1429 rows

#save the edited file to disk
write.csv(dat, "data_raw/sheet4_popgen_database_analysis_clean_edit_11Apr23.csv", row.names=FALSE)

#check if any of the rows with the columns Species, COI, Order, Longitude and Latitude are duplicated
dat %>%
group_by(Species, Voucher_no, COI, Order, Longitude, Latitude) %>%
filter(n()>1) %>%
as.data.frame()
#no duplicates

dat %>%
group_by(Voucher_no, COI) %>%
filter(n()>1) %>%
as.data.frame()
#no duplicates

#getting a vector of unique COI accession numbers, which can be queried against GenBank - also removing new sequences from our work where CCMB ids have been given against the accession number, or Jahnavi's previous work where CES ids have been used
seq.acc<-dat %>%
filter(COI!="") %>%
filter(!grepl("CCMB|CES", COI)) %>%
distinct(COI) %>%
pull()

length(seq.acc) #1376 sequences

#get a vector of unique species names for which GenBank accession numbers are available
sp<-dat %>%
filter(COI!="") %>%
filter(!grepl("CCMB|CES", COI)) %>%
distinct(Species) %>%
pull()
length(sp) #132 species

#splitting the input accession numbers into chunks of 300 sequences so that these can be downloaded as the server doesn't handle large searches. From https://github.com/ropensci/rentrez/issues/97
chunk_size<-300
seq.acc<-split(seq.acc, ceiling(seq_along(seq.acc)/chunk_size))
length(seq.acc) #5 chunks of accession numbers

#retrieve the sequences from NCBI
#create an object where the sequences will be stored
seq<-character()
for(i in 1:length(seq.acc)){
seq.chunk<-entrez_fetch(db="nuccore", id=seq.acc[[i]], web_history=NULL, rettype="fasta")
seq<-c(seq, seq.chunk)
}

#create a directory to store processed data
dir.create("data")

#saving the sequences to disk
write(seq, file="data/COI_11Apr23.fasta")

#re-input the sequence file using read.dna so that individual sequences are identified
seq<-read.dna("data/COI_11Apr23.fasta", format="fasta")
length(seq) #1376 sequences

#retain only the accession number in the sequence name
names(seq)<-substr(names(seq), 1, 8)

#export this alignment to disk
write.FASTA(seq, file="data/COI_new_labels_11Apr23.fasta", append=FALSE)

### perform multiple sequence alignment for each species ####

#the for loop subsets sequences based on accession numbers corresponding to each species, performs sequence alignment and exports the unaligned and aligned sequences to disk

#create folders where sequence subsets of each species can be exported
dir.create("data/seq_subset")
dir.create("data/seq_aln")

for(i in 1:length(sp)){
#get a vector that lists the COI accession numbers corresponding to a species
labels.sub<-unique(dat$COI[dat$Species==sp[i]])

#find the ids of these accession numbers in the master fasta file
ids<-which(labels(seq) %in% labels.sub)
seq.sub<-seq[ids] #since the sequences are not aligned, there are no rows, each sequence is a separate object and therefore they are subset as a vector

#check if we've got the correct sequences
if(!identical(labels(seq.sub), labels.sub)){
print(paste0("The accession numbers in the database do not match the accession numbers selected from the master fasta file downloaded from GenBank for ", sp[i])) #it is not identical for species exclusively added from our database (not listed in GenBank) and for Scolopendra morsitans, Scutigera coleoptrata, because our database has an empty cell in place of the COI accession number in one of the entries
}

#export this unaligned subset of sequences as a fasta file
file.name<-paste0(gsub(" ", "_", sp[i]), ".fasta")
write.FASTA(seq.sub, file=paste0("data/seq_subset/", file.name), append=FALSE)

#import this file in as a StringSet format to perform multiple sequence alignment
seq.sub<-readDNAStringSet(paste0("data/seq_subset/", file.name), format="fasta")

#running multiple sequence alignment using muscle
aln.sub<-muscle::muscle(seq.sub) #I am using default parameters here
class(aln.sub)

#convert aln.sub which is in a DNAMultipleAlignment format into the DNAbin format
aln.sub<-as.DNAbin(aln.sub)

#export alignment in fasta format
write.FASTA(aln.sub, file=paste0("data/seq_aln/aln_", file.name), append=FALSE)
}

#### I viewed the alignments in Aliview, cropped ends with many Ns for a species alignment and removed sequences with a lot of missing data. In some of the alignments had to reverse, complement or reverse-complement some of the sequences, which were not in the same reading frame as the rest. For species where we had unpublished sequences from our work in CCMB or from Jahnavi's earlier work from CES, we realigned sequences including them. These edited alignments are stored in the folder 'final' and are used to calculate sequence summaries in code03. ####

####