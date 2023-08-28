import glob



base = "/mnt/disk1/PROJECTS/SURPRISE/2022_bats/unit_filtered/"
suffix_1 = "_R1.fastq"
suffix_2 = "_R2.fastq"
base2 = "/mnt/disk1/PROJECTS/SURPRISE/test/"




SAMPLES = [s.split('/')[-1].replace(suffix_1, '') for s in glob.glob(base+'*'+suffix_1)]


print(SAMPLES)
print(base)
print(base2)
Output_blasts = ["results/blast_ICTV/contigs_ICTV_na_vir_blastn/", "results/blast_ICTV/contigs_ICTV_na_vir_blastn/", "results/blast_ICTV/contigs_ICTV_na_vir_blastn/"]
want_all = []
want_all.append(expand(base2+"results/fastp/{sample}_R1.fastq.gz", sample=SAMPLES))
want_all.append(expand(base2+"results/contigs/{sample}.fasta", sample=SAMPLES))
want_all.append(expand(base2+"results/contigs/filt_contigs/{sample}_l500.fasta", sample=SAMPLES))
want_all.append(expand(base2+"results/kraken2/contigs/{sample}.out", sample=SAMPLES))
want_all.append(expand(base2+"results/kraken2/contigs/{sample}.contigs.taxo.txt", sample=SAMPLES))
want_all.append(expand(base2+"results/catbat/catbat_res_{sample}_c2c.txt", sample=SAMPLES))
want_all.append(expand(base2+"results/catbat/catbat_res_{sample}_true.names", sample=SAMPLES))
want_all.append(expand(base2+"results/catbat/catbat_res_{sample}.contigs.taxo.txt", sample=SAMPLES))
want_all.append(expand(base2+"results/viruses_filtered_contigs/{sample}_l500_viruses_filtered_contigs.fasta", sample=SAMPLES))
want_all.append(expand(base2+"results/blast_ICTV/contigs_ICTV_na_vir_blastn/{sample}_l500_blastres.out", sample=SAMPLES))
want_all.append(expand(base2+"results/blast_RVDB/contigs_RVDB_na_vir_blastn/{sample}_l500_blastres.out", sample=SAMPLES))
want_all.append(expand(base2+"results/blast_NCBI/contigs_NCBI_na_vir_blastn/{sample}_l500_blastres.out", sample=SAMPLES))
want_all.append(expand(base2+"results/contigs_viralverify/{sample}/Prediction_results_fasta/{sample}_l500_viruses_filtered_contigs_all_viruses.txt", sample=SAMPLES))
want_all.append(expand(base2+"results/contigs_blastres_ICTV_RVDB_NCBI_with_host_ICTV_and_NCBI.xlsx"))
want_all.append(expand(base2+"results/unknown_contigs/{sample}_l500_blastres_viruses_contig.txt", sample=SAMPLES))



rule all:
    input: want_all

rule fastp:
    input: read1=base+"{sample}"+suffix_1,
           read2=base+"{sample}"+suffix_2
    output: read1_trim=base2+"results/fastp/{sample}_R1.fastq.gz",
            read2_trim=base2+"results/fastp/{sample}_R2.fastq.gz"
    threads: 30
    conda: "envs/fastp.yaml"
    shell: """
            fastp -i {input.read1} -I {input.read2} -o {output.read1_trim} -O {output.read2_trim} -q 20 -l 30 -D
    """

rule contigs:
    input: read1_trim=base2+"results/fastp/{sample}_R1.fastq.gz",
           read2_trim=base2+"results/fastp/{sample}_R2.fastq.gz"
    output: res_file=base2+"results/contigs/{sample}.fasta"
    params: outfolder=base2+"results/contigs/{sample}"
    threads: 30
    conda: "envs/contigs_env.yaml"
    shell: """
            megahit -t {threads} -m 400 -1 {input.read1_trim} -2 {input.read2_trim} -o {params.outfolder}
            #spades.py --meta -t {threads} -m 400 --only-assembler -1 {input.read1_trim} -2 {input.read2_trim} -o {params.outfolder}
            #spades.py --metaviral -t {threads}  -m 400 --only-assembler -1 {input.read1_trim} -2 {input.read2_trim} -o {params.outfolder}
            #spades.py --rnaviral -t {threads} -m 400 --only-assembler -1 {input.read1_trim} -2 {input.read2_trim} -o {params.outfolder}
            cp {params.outfolder}/final.contigs.fa {output.res_file}
    """

rule contigs_filtering:
    input: contigs_raw=base2+"results/contigs/{sample}.fasta"
    output: contigs_filt=base2+"results/contigs/filt_contigs/{sample}_l500.fasta"
    threads: 30
    conda: "envs/seqkit_env.yaml"    
    shell: """
           seqkit seq -m 500 -j {threads} -g {input.contigs_raw} > {output.contigs_filt}           
    """       

rule kraken2:
    input: contigs_filt=base2+"results/contigs/filt_contigs/{sample}_l500.fasta"
    output: kraken2_report=base2+"results/kraken2/contigs/{sample}.report",
            kraken2_out=base2+"results/kraken2/contigs/{sample}.out"
    params: kraken2_db="/mnt/disk1/DATABASES/kraken2/pluspf"
    conda: "envs/kraken2_env.yaml"
    threads: 30
    shell: """     
        /usr/bin/perl $(which kraken2) --threads {threads} --confidence 0.8 --db {params.kraken2_db} {input.contigs_filt} --minimum-hit-groups 3 --report-minimizer-data --use-names --report {output.kraken2_report} --output {output.kraken2_out}
    """
    
rule kraken_contigs_2nd_taxo:
    input:  kraken2_contigs_out=base2+"results/kraken2/contigs/{sample}.out"
    output: kraken2_contigs_taxo=base2+"results/kraken2/contigs/{sample}.contigs.taxo.txt"
    shell: """
          grep 'C' {input.kraken2_contigs_out} | cut -f2,3 > {output.kraken2_contigs_taxo}
    """

rule kraken2_taxonomy_contigs:
    input: kraken2_spcs=base2+"results/kraken2/contigs/{sample}.contigs.taxo.txt"
    output: kraken2_contigs_taxo=base2+"results/kraken2/contigs/{sample}.contigs_cellular_contigs.txt"
    conda: "envs/R-packages_env.yaml"
    shell: """
          Rscript r_scripts/taxonomizr_kraken2_filtered_eukaryote.R {input.kraken2_spcs}
    """

rule catbat_contigs:
    input: contigs_filt=base2+"results/contigs/filt_contigs/{sample}_l500.fasta"
    output: c2c= base2+"results/catbat/catbat_res_{sample}_c2c.txt"
    params: ref_db="/mnt/disk1/DATABASES/CAT/CAT_prepare_20210107/2021-01-07_CAT_database",
            taxo_ref="/mnt/disk1/DATABASES/CAT/CAT_prepare_20210107/2021-01-07_taxonomy",
            prefix= base2+"results/catbat/catbat_res_{sample}"
    conda: "envs/catbat_env.yaml"
    threads: 30
    shell: """
          CAT contigs -c {input.contigs_filt} --index_chunks 1 --out_prefix {params.prefix} -f 0.9 -r 5 -n {threads} -d {params.ref_db} -t {params.taxo_ref} --force
          mv {params.prefix}.contig2classification.txt {output.c2c}
    """

rule catbat_human_names:
    input:c2c=base2+"results/catbat/catbat_res_{sample}_c2c.txt"
    output:true_names=base2+"results/catbat/catbat_res_{sample}_true.names"
    params: taxo_ref="/mnt/disk1/DATABASES/CAT/CAT_prepare_20210107/2021-01-07_taxonomy"
    conda: "envs/catbat_env.yaml"
    shell: """
          CAT add_names -i {input.c2c} -o {output.true_names} -t {params.taxo_ref} --only_official
   """

rule catbat_true_taxo:
    input: catbat_truenames=base2+"results/catbat/catbat_res_{sample}_true.names"
    output: catbat_taxo=base2+"results/catbat/catbat_res_{sample}.contigs.taxo.txt"
    shell: """
         cut -f1,4 {input.catbat_truenames} | awk '{{(n=split ($2,a,";"));print $1,"\t",a[n]}}' | sed 's/ \t /\t/g' | sed '1d' > {output.catbat_taxo}
    """
    
rule catbat_taxonomy_contigs:
    input: catbat_taxo=base2+"results/catbat/catbat_res_{sample}.contigs.taxo.txt"
    output: catbat_contigs_anno=base2+"results/catbat/catbat_res_{sample}.contigs_catbat_cellular_contigs.txt"
    conda: "envs/R-packages_env.yaml"
    shell: """
          Rscript r_scripts/taxonomizr_catbat_filtered_eukaryote.R {input.catbat_taxo}

    """
    
rule contigs_negativ_selection:
    input:  kraken2_contigs_taxo=base2+"results/kraken2/contigs/{sample}.contigs_cellular_contigs.txt",
            catbat_contigs_anno=base2+"results/catbat/catbat_res_{sample}.contigs_catbat_cellular_contigs.txt",
            seqkit_out=base2+"results/contigs/filt_contigs/{sample}_l500.fasta"
    output: seqkit_out=base2+"results/viruses_filtered_contigs/{sample}_l500_viruses_filtered_contigs.fasta"
    params: cellular=base2+"results/viruses_filtered_contigs/{sample}_cellular.contigs.txt"
    threads: 30
    conda: "envs/seqkit_env.yaml"
    shell: """
           cat {input.kraken2_contigs_taxo} {input.catbat_contigs_anno} | sort | uniq > {params.cellular}
           seqkit grep -v -f {params.cellular} {input.seqkit_out} -j {threads} -o {output.seqkit_out}
           sed -i -e 's/ /_/g' {output.seqkit_out}
    """

rule contigs_blasts_ICTV:
    input: contigs_filt=base2+"results/viruses_filtered_contigs/{sample}_l500_viruses_filtered_contigs.fasta"
    output: blast_out_na=base2+"results/blast_ICTV/contigs_ICTV_na_vir_blastn/{sample}_l500_blastres.out"
    params: ref="/mnt/disk1/DATABASES/ICTV/na"
    threads: 5
    conda: "envs/blast_env.yaml"
    shell: """
           blastn -query {input.contigs_filt} -db {params.ref} -num_threads {threads} -evalue 1e-5 -task 'blastn' -outfmt '6 std slen qlen stitle' -out {output.blast_out_na}
    """
    
rule contigs_blasts_RVDB:
    input: contigs_filt=base2+"results/viruses_filtered_contigs/{sample}_l500_viruses_filtered_contigs.fasta"
    output: blast_out_na=base2+"results/blast_RVDB/contigs_RVDB_na_vir_blastn/{sample}_l500_blastres.out"
    params: ref="/mnt/disk1/DATABASES/RVDBase/C-RVDBv25.0.filtered_HIV_NA"
    threads: 5
    conda: "envs/blast_env.yaml"
    shell: """
           blastn -query {input.contigs_filt} -db {params.ref} -num_threads {threads} -evalue 1e-5 -task 'blastn' -outfmt '6 std slen qlen stitle' -out {output.blast_out_na}
    """
    
rule contigs_blasts_NCBI_viruses:
    input: contigs_filt=base2+"results/viruses_filtered_contigs/{sample}_l500_viruses_filtered_contigs.fasta"
    output: blast_out_na=base2+"results/blast_NCBI/contigs_NCBI_na_vir_blastn/{sample}_l500_blastres.out"
    params: ref="/mnt/disk1/DATABASES/NA_refseq_NCBI_viruses_11072023/NA_refseq_NCBI_viruses_11072023"
    threads: 5
    conda: "envs/blast_env.yaml"
    shell: """
           blastn -query {input.contigs_filt} -db {params.ref} -num_threads {threads} -evalue 1e-5 -task 'blastn' -outfmt '6 std slen qlen stitle' -out {output.blast_out_na}
    """
    
rule contigs_viralverify:
    input: contigs_filt=base2+"results/viruses_filtered_contigs/{sample}_l500_viruses_filtered_contigs.fasta"
    output: contigs_viralverify1=base2+"results/contigs_viralverify/{sample}/Prediction_results_fasta/{sample}_l500_viruses_filtered_contigs_virus.fasta",
            contigs_viralverify2=base2+"results/contigs_viralverify/{sample}/Prediction_results_fasta/{sample}_l500_viruses_filtered_contigs_virus_uncertain.fasta",
            contigs_viralverify3=base2+"results/contigs_viralverify/{sample}/Prediction_results_fasta/{sample}_l500_viruses_filtered_contigs_all_viruses.fasta",
            contigs_viralverify4=base2+"results/contigs_viralverify/{sample}/Prediction_results_fasta/{sample}_l500_viruses_filtered_contigs_all_viruses.txt"
    params: db="/mnt/disk1/DATABASES/nbc_hmms/nbc_hmms.hmm",
            outfolder=base2+"results/contigs_viralverify/{sample}"
    threads: 30
    conda: "envs/viralverify_env.yaml"
    shell: """
           viralverify -f {input.contigs_filt} -t {threads} --hmm {params.db} -o {params.outfolder}
           cat {output.contigs_viralverify1} {output.contigs_viralverify2} > {output.contigs_viralverify3} || true
           cat {output.contigs_viralverify3} | grep ">" | sed 's/>//' > {output.contigs_viralverify4} || true
    """

rule make_xlsx:
    input: expand(base2+"{output_blasts}{sample}_l500_blastres.out", sample=SAMPLES, output_blasts = Output_blasts)
    params: location=base2,
            cutoff=SAMPLES
    output: final_table=base2+"results/contigs_blastres_ICTV_RVDB_NCBI_with_host_ICTV_and_NCBI.xlsx"
    shell: """
           Rscript r_scripts/final_table.R {params.location} {params.cutoff}
    """
    
    
rule contigs_blast_uniq:
    input:  blast_out_na_ICTV=base2+"results/blast_ICTV/contigs_ICTV_na_vir_blastn/{sample}_l500_blastres.out",
            blast_out_na_RVDB=base2+"results/blast_RVDB/contigs_RVDB_na_vir_blastn/{sample}_l500_blastres.out",
            blast_out_na_NCBI=base2+"results/blast_NCBI/contigs_NCBI_na_vir_blastn/{sample}_l500_blastres.out",
            contigs_filt=base2+"results/contigs_viralverify/{sample}/Prediction_results_fasta/{sample}_l500_viruses_filtered_contigs_all_viruses.fasta"
    output: blast_out=base2+"results/unknown_contigs/{sample}_l500_blastres_viruses_contig.txt",
            seqkit_out=base2+"results/unknown_contigs/{sample}_l500_unknown_viruses.fasta"
    params: awk="'{print $1}'"
    threads: 30
    conda: "envs/seqkit_env.yaml"
    shell: """
           cat {input.blast_out_na_ICTV} {input.blast_out_na_RVDB}  {input.blast_out_na_NCBI} | awk {params.awk} | sort | uniq  > {output.blast_out}   
           seqkit grep -f {output.blast_out} {input.contigs_filt} -j {threads} -o {output.seqkit_out} -v
    """
    