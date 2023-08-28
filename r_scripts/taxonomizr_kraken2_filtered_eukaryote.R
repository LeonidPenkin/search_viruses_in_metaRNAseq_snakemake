#install.packages("taxonomizr")
#install.packages("writexl")
#install.packages("dplyr")

library(taxonomizr)
library(writexl)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)[1]
#args <- "/mnt/disk1/PROJECTS/SURPRISE/old_bats/results/kraken2/contigs/bat2-2015.contigs.taxo.txt"
if (file.size(args) == 0L) {
  print("File is empty")
  args2 = gsub(".taxo.txt","_cellular_contigs.txt",args)
  cat(NULL,file=args2)
  quit() }
file <- read.csv(args, sep='\t',header = F)
#file <- read.csv("/mnt/disk1/PROJECTS/SURPRISE/old_bats/results/kraken2/contigs/22hed15.contigs.taxo.txt",sep='\t',header=F)
file$taxaid <- as.numeric(gsub("\\D", "", file[,2]))
colnames(file) <- c('contig_name','taxa', 'taxaid')
#file<-file[file$taxaid > 0,]
file[,4:10] <- as.data.frame(getTaxonomy(file[[3]],"/mnt/disk1/DATABASES/taxonomizr_data/taxa.sql"))
file <- file %>% filter(superkingdom != "Viruses" | taxaid != 1)
output_name = gsub(".taxo.txt","_cellular_contigs.txt",args)
cat(file$contig_name, file=output_name, sep="\n")



