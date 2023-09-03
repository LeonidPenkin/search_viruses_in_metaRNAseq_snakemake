library(GenomicRanges)
library(plyr)
library(dplyr)
library(writexl)


data_frames <-  list()
ICTV_host <-read.csv2("/mnt/disk1/PROJECTS/SURPRISE/old_bats/ICTV_host.txt", sep = "\t", head = T, na.strings=c("","NA"))
NCBI_host <-read.csv2("/mnt/disk1/PROJECTS/SURPRISE/old_bats/NA_refseq_NCBI_viruses_11072023_hosts.txt", sep = "|", head = T, na.strings=c("","NA"))

base2_location_output_file = commandArgs(trailingOnly=TRUE)[1]
sample_names <- commandArgs(trailingOnly=TRUE)[2:length(commandArgs(trailingOnly=TRUE))]
filtered_len <- 500

bases <-c(paste0(base2_location_output_file, "results/blast_ICTV/contigs_ICTV_na_vir_blastn/"), 
          paste0(base2_location_output_file, "results/blast_RVDB/contigs_RVDB_na_vir_blastn/"), 
          paste0(base2_location_output_file, "results/blast_NCBI/contigs_NCBI_na_vir_blastn/"))
metods <- c("ICTV","RVDB","NCBI")


for(location in 1:length(bases)){
  base <-bases[location]
  metod <-metods[location]
  cat(metod, sep="\n")
  file_names <- sample_names
  for(file_name in file_names){
    filename=paste0(base, file_name, "_l500_blastres.out")
    out_blast <- read.csv2(filename, header = F, sep = "\t")
    colnames(out_blast) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "slen", "qlen", "stitle")
    out_blast$start <- ifelse(out_blast$sstart < out_blast$send, out_blast$sstart, out_blast$send)
    out_blast$end <- ifelse(out_blast$sstart > out_blast$send, out_blast$sstart, out_blast$send)
    out_blast$strand <- ifelse(out_blast$sstart < out_blast$send, "+", "-")
    ifelse(metod == "NCBI", sseqid_index <- 1, sseqid_index <- 3)
    out_blast$Virus.GENBANK.accession.blastres <-data.frame(do.call('rbind', strsplit(out_blast$sseqid, split="[|]")))[,sseqid_index]
    out_blast$Virus.GENBANK.accession.blastres <- sub('..$*', '', out_blast$Virus.GENBANK.accession.blastres)
    out_blast$Virus.GENBANK.accession.blastres <- sub('[.]', '', out_blast$Virus.GENBANK.accession.blastres)
    gr <- GRanges(
      seqnames = Rle(c(out_blast$Virus.GENBANK.accession.blastres)),
      ranges = IRanges(out_blast$start, end = out_blast$end),
      strand = Rle(out_blast$strand))
    viruses <- as.data.frame(reduce(gr))
    colnames(viruses) <- c("Virus.GENBANK.accession.blastres", "sstart", "send", "length", "strand")
    out_blast2 <- out_blast[,c(1:2,13,15,19)]
    out_blast3 <- out_blast2[!duplicated(out_blast2$Virus.GENBANK.accession.blastres), ]
    viruses <- ddply(viruses, .(Virus.GENBANK.accession.blastres), summarise, sum_nucleotides=sum(length))
    viruses <- merge(x = viruses, y = out_blast3, by = "Virus.GENBANK.accession.blastres",  all.x = TRUE, sort = TRUE)
    ifelse(grepl("len", viruses$qseqid[1]), 
           viruses$contig_len <- as.numeric(data.frame(do.call('rbind', strsplit(as.character(viruses$qseqid),'_len=',fixed=TRUE)))[,2]), 
           viruses$contig_len <- as.numeric(data.frame(do.call('rbind', strsplit(as.character(viruses$qseqid),'_',fixed=TRUE)))[,4]))
    
    viruses$coverege <- viruses$sum_nucleotides / viruses$slen
    
    viruses$locathion_file <- filename
    viruses$filtered_len <- filtered_len
    viruses$metod <- metod
    cat(file_name, " ", metod, sep="\n")
    for(i in 1:nrow(viruses)){
      cat(nrow(viruses)-i, " ")
      host_info <- ICTV_host[(grepl(viruses$Virus.GENBANK.accession.blastres[i], ICTV_host$Virus.GENBANK.accession) | grepl(viruses$Virus.GENBANK.accession.blastres[i], ICTV_host$Virus.REFSEQ.accession)), 26]
      ifelse(length(host_info) == 0, viruses[i, 12:35] <- NA , viruses[i, 12:35] <- ICTV_host[(grepl(viruses$Virus.GENBANK.accession.blastres[i], ICTV_host$Virus.GENBANK.accession) | grepl(viruses$Virus.GENBANK.accession.blastres[i], ICTV_host$Virus.REFSEQ.accession)), 3:26])
      
      NCBI_host_info <- NCBI_host[grepl(viruses$Virus.GENBANK.accession.blastres[i], NCBI_host$Virus.GENBANK.accession), 3]
      ifelse(length(NCBI_host_info) == 0, viruses[i, 36:37] <- NA , viruses[i, 36:37] <- NCBI_host[grepl(viruses$Virus.GENBANK.accession.blastres[i], NCBI_host$Virus.GENBANK.accession), c(1,3)])
      
      
    }
    colnames(viruses) <- c("Virus.GENBANK.accession.blastres", "sum_nucleotides", "qseqid", "sseqid", "slen", "stitle", "contig_len", "coverege", "locathion_file", "filtered_len", "metod", "Realm", "Subrealm", "Kingdom", "Subkingdom", "Phylum", "Subphylum", "Class", "Subclass", "Order", "Suborder", "Family", "Subfamily", "Genus", "Subgenus", "Species", "Exemplar.or.additional.isolate", "Virus.name.s.", "Virus.name.abbreviation.s.", "Virus.isolate.designation", "Virus.GENBANK.accession", "Virus.REFSEQ.accession", "Genome.coverage", "Genome.composition", "Host.source", "NCBI_virus_name", "NCBI_host")
    viruses$sample_name <- sub('_l500_blastres.out', '', file_name)
    
    ifelse(metod == "RVDB", stitle_index <- 4, stitle_index <- 3)
    ifelse(metod == "NCBI", stitle_index <- 2, stitle_index <- 3)    
    viruses$name_viruses_in_database <- data.frame(do.call('rbind', strsplit(as.character(viruses$stitle),'|',fixed=TRUE)))[,stitle_index]
 
    viruses$coverege_round <- round(viruses$coverege, digits = 1)
    df <-  data.frame(file_name, viruses)
    file_name_df=paste0(file_name, "_", filtered_len, "_", metod)
    data_frames[[file_name_df]] <- df
  }
}

super_mega_ultra_tabl <- do.call("rbind", data_frames)
write_xlsx(super_mega_ultra_tabl, paste0(base2_location_output_file, "results/contigs_blastres_ICTV_RVDB_NCBI_with_host_ICTV_and_NCBI.xlsx"))
