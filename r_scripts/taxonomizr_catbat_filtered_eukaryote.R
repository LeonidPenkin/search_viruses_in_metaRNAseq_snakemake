#install.packages("taxonomizr")
library(taxonomizr)
library(writexl)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)[1]
if (file.size(args) == 0L) {
  print("File is empty")
  output_name = gsub(".taxo.txt","_catbat_anno.xlsx",args)
  cat(NULL,file=args2)
  quit() }
file <- read.csv(args, sep='\t',header = F)
#file <- read.csv("/mnt/disk1/PROJECTS/SURPRISE/old_bats/results/catbat/catbat_res_22hed15.contigs.taxo.txt",sep='\t',header=F)
colnames(file) <- c('contig_name','taxaid')
file[,3:9] <- as.data.frame(getTaxonomy(file[[2]],"/mnt/disk1/DATABASES/taxonomizr_data/taxa.sql"))
file$contig_length <- as.numeric(lapply(strsplit(as.character(file$contig_name),'_'),"[",4))
file$coverage <- as.numeric(lapply(strsplit(as.character(file$contig_name),'_'),"[",6))
file <- file %>% filter(superkingdom != "Viruses" & taxaid != 1)
output_name = gsub(".taxo.txt","_catbat_cellular_contigs.txt",args)
cat(file$contig_name, file=output_name, sep="\n")

#args="/mnt/disk1/PROJECTS/SURPRISE/old_bats/results/catbat/catbat_res_22hed15.contigs.taxo.txt"
