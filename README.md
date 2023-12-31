# search_viruses_in_metaRNAseq_snakemake

Snakemake file and R scripts for search viruses of metagenomic RNA reads made in Laboratory of Bioinformatics at Research Institute for Systems Biology and Medicine (RISBM).
This algorithm is a continuation of the idea of the taxo_class_snakemake algorithm (https://github.com/ISonets/taxo_class_snakemake/tree/main). This algorithm is adapted to search for viruses in meta RNAseq samples. The algorithm uses a script to search three virus sequence databases (NCBI, RVDB, ICTV). Sequences not found in this data are additionally checked using the algorithm viralverify.

## Prerequisites:

  - Snakemake (https://snakemake.github.io/);  
  
  - Kraken2 >= 2.1.2 (https://github.com/DerrickWood/kraken2);
  
    **Note**: kraken2 needs to be installed into another conda env **or** installed from source due to dependecy conflicts between snakemake and kraken2.
    
    **Note**: use pluspf database for kraken2 (https://benlangmead.github.io/aws-indexes/k2).
    
  - CAT/BAT and its local database (https://github.com/dutilh/CAT);
  - Hisat2 (http://daehwankimlab.github.io/hisat2/);
  - samtools (https://github.com/samtools/samtools);
  - SPAdes (https://github.com/ablab/spades);
  - BLAST and its local nt database (https://blast.ncbi.nlm.nih.gov/Blast.cgi);
  - International Committee on Taxonomy of Viruses: ICTV (https://ictv.global/);
  - Reference Viral DataBase (RVDB) (https://rvdb.dbi.udel.edu/);
  - viralVerify: viral contig verification tool (https://github.com/ablab/viralVerify);
  - seqkit (https://github.com/shenwei356/seqkit)
  - R and its libraries:
    - taxonomizr (https://github.com/sherrillmix/taxonomizr),
    - dplyr (https://dplyr.tidyverse.org/), 
    - readxl (https://readxl.tidyverse.org/),
    - plyr (http://had.co.nz/plyr/ ),
    - reshape (https://cran.r-project.org/web/packages/reshape/index.html);
  
    **Note**: for taxonomizr library you need local database, but for our purposes accession ids do not needed, to achieve this, use `getAccessions=FALSE`
 
 ## Usage:
 
 1. Install all tools and verify if it's working;
 2. Open Snakemake file and provide paths to FASTQ files;
 3. Launch it and make some tea while it's doing science;
 4. Analyse results provided for you in human-readable form:)
 
 **Note**: you might also need to modify paths in R scripts to accomodate your data and folder structure, because they are hardcoded now.

## Authors

Leonid Penkin

Ignat Sonets

Alexander Manolov, PhD
