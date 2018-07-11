# Genomes_project
## Statistical analysis of genomes
### Goals and objectives
Find correlation between genomic features (like SNPs, methylation, TFBS) and functional genomic regions in different genomes

1. Plot sequence features such as TFBS, SNPs, methylation, RNA-seq coverage
2. Map it on functional genomic regions
3. Find correlation and check reproducibility for different genomes
4. Consider annotation quality and outcomes for functional features (like promoters)prediction for not annotated genomes


### Data:
#### *Arabidopsis thaliana*
1. Graphs for *Oryza sativa* [[1](https://www.ncbi.nlm.nih.gov/pubmed/27774999)]
2. gff3 file (**TAIR10_GFF3_genes.gff**) was taken [here](https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FTAIR10_genome_release%2FTAIR10_gff3)
2. [vcf file 1001 genome TAIR]()
3. *A. thaliana* reference genome (TAIR10_toplevel) was taken from ftp://ftp.ensemblgenomes.org/pub/plants/release-39/fasta/arabidopsis_thaliana/dna/
4. methylation data was taken from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE105000

#### *Medicago truncata*
1. annotation (.gff) and assemly (.fasta) from http://www.medicagogenome.org/downloads
2. SNP files also from http://www.medicagogenome.org/downloads

#### *Homo sapiens*
1. annotation [Release 28 (GRCh38.p12)](https://www.gencodegenes.org/releases/current.html) (CHR) in .gff3 format
2. [.fasta of primary assembly (PRI)](https://www.gencodegenes.org/releases/current.html)


### Scripts for data preprocessing:
1. get_ATGs.py
2. get_4tss.py
4. get_promoters.py
3. get_fin_anno.py 

### Data preprocessing:
1. to create file with ATGs: ``` python3 get_ATGs.py annotation.gff ``` 
2. to create file with tss: ``` python3 get_4tss.py annotation.gff ``` 
3. to create files with promoter regions (.bed + .txt): ``` python3 get_promoters.py 4tss.txt ``` 
4. to obtain promoter regions sequences: ``` sed 's/^>1.*$/>Chr1/' Arabidopsis_thaliana.TAIR10.dna.toplevel.fa | sed 's/^>2.*$/>Chr2/' | sed 's/^>3.*$/>Chr3/'| sed 's/^>4.*$/>Chr4/'| sed 's/^>5.*$/>Chr5/'| sed 's/^>Mt.*$/>ChrM/'| sed 's/^>Pt.*$/>ChrC/' > new_ref.fa ``` in order to get names of chromosomes in fasta consistent with names in bed file, then ``` bedtools getfasta -fi corrected_reference.fasta -bed promoters.bed -name -s -fo promoters_sequences.fasta ``` 
5. to create fin_anno: ``` python3 get_fin_anno.py annotation.gff ``` 

### Plots visualization:
1. first (and the most important) file is **snp_custom_annotation.r**, which contains a function that create custom annotation of snps - all other scripts use these function
2. **ATG_plot.r** is used for visualization SNP distribution around start codon (required packages are dplyr, scales)
3. **intron_exon_junctions.r** is used for visualization of SNP distribution around exon-intron boundary
4. **promoter-terminator.r** is used for visualization of SNP distribution around terminator
5. **transcr_stop_plot.r** is used for visualization of SNP distribution around transcription stop codon
6. **transfac.r** is used for visualization distribution of TFBSs in promoter region (+-500 nucleotides around TSS)

### Several results:
1. *Arabidopsis thaliana*
  ![](https://raw.githubusercontent.com/danchurova/Genomes_project/master/Arabidopsis/pictures/OKintron_exon_junction_arab.jpeg?raw=true)
  
![](https://github.com/danchurova/Genomes_project/blob/master/Arabidopsis/pictures/OK-ATG.jpeg?raw=true)

![](https://github.com/danchurova/Genomes_project/blob/master/Arabidopsis/pictures/OKpromoter-terminator-arab.jpeg?raw=true)

![](https://github.com/danchurova/Genomes_project/blob/master/Arabidopsis/pictures/sp_cells-total.jpeg?raw=true)

![](https://github.com/danchurova/Genomes_project/blob/master/Arabidopsis/pictures/veg-cell-total.jpeg?raw=true)

2. *Medicago truncata*
![](https://github.com/danchurova/Genomes_project/blob/master/Medicago/pictures/ATG.jpeg?raw=true)
![](https://github.com/danchurova/Genomes_project/blob/master/Medicago/pictures/intron_exon_med.jpeg?raw=true)
![](https://github.com/danchurova/Genomes_project/blob/master/Medicago/pictures/promoter-terminator-med.png?raw=true)
