# Genomes_project
## Statistical analysis of annotated genomes
### Goals and objectives
Find correlation between genomic features (like SNPs, methylation, TFBS) and functional genomic regions in different genomes

1. Plot sequence features such as TFBS, SNPs, methylation, RNA-seq coverage
2. Map it on functional genomic regions
3. Find correlation and check reproducibility for different genomes
4. Consider annotation quality and outcomes for functional features (like promoters)prediction for not annotated genomes


### Data:
Graphs for *Oryza sativa* [[1](https://www.ncbi.nlm.nih.gov/pubmed/27774999)]
#### *Arabidopsis thaliana*
1. reference genome TAIR10_toplevel (ftp://ftp.ensemblgenomes.org/pub/plants/release-39/fasta/arabidopsis_thaliana/dna/)
2. annotation [**TAIR10_GFF3_genes.gff3**](https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FTAIR10_genome_release%2FTAIR10_gff3)
2. variation [vcf file 1001 genome TAIR]()
4. [methylation data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE105000)

#### *Medicago truncata*
1. annotation (.gff) and assemly (.fasta) from http://www.medicagogenome.org/downloads
2. SNP files also from http://www.medicagogenome.org/downloads

#### *Homo sapiens*
1. annotation [Release 28 (GRCh38.p12)](https://www.gencodegenes.org/releases/current.html) (CHR) in .gff3 format
2. [.fasta of primary assembly (PRI)](https://www.gencodegenes.org/releases/current.html)

#### *Mus musculus*
1. annotation [Release M17 (GRCm38.p6)](https://www.gencodegenes.org/mouse_releases/current.html) (CHR) in .gff3 format
2. [.fasta of primary assembly (PRI)](https://www.gencodegenes.org/mouse_releases/current.html)

#### *Felis catus*
1. annotation [assembly Felis_catus_9.0](https://www.ncbi.nlm.nih.gov/genome?term=felis%20catus) in .gff format (ID 78)
2. [.fasta of assembly 9.0](https://www.ncbi.nlm.nih.gov/genome?term=felis%20catus) (ID 78)

#### *Drosophila melanogaster*
1. reference [assembly dmel_r5.57_FB2014_03 from FlyBase, **dmel-all-chromosome-r5.57.fasta.gz**](ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.57_FB2014_03/fasta/)
2. annotation [dmel_r5.57_FB2014_03 **dmel-all-filtered-r5.57.gff.gz**](ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.57_FB2014_03/gff/)
3. variation downloaded for each chromosome for all populations in one file in .vcf format[PopFly Browser](https://popfly.uab.cat) [Hervas S, Sanz E, Casillas S, Pool JE, and Barbadilla A (2017) PopFly: the Drosophila population genomics browser.
Bioinformatics, 33, 2779-2780;](https://doi.org/10.1093/bioinformatics/btx301)


#### *Danio rerio*
1. 





### Scripts for data preprocessing:
1. get_ATGs.py
2. get_4tss.py
3. get_4tts.py
4. get_promoters.py
5. get_fin_anno.py 

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
  ![](https://raw.githubusercontent.com/danchurova/Genomes_project/master/Arabidopsis_thaliana/pictures/OKintron_exon_junction_arab.jpeg?raw=true)
  
![](https://github.com/danchurova/Genomes_project/blob/master/Arabidopsis_thaliana/pictures/OK-ATG.jpeg?raw=true)

![](https://github.com/danchurova/Genomes_project/blob/master/Arabidopsis_thaliana/pictures/OKpromoter-terminator-arab.jpeg?raw=true)

![](https://github.com/danchurova/Genomes_project/blob/master/Arabidopsis_thaliana/pictures/sp_cells-total.jpeg?raw=true)

![](https://github.com/danchurova/Genomes_project/blob/master/Arabidopsis_thaliana/pictures/veg-cell-total.jpeg?raw=true)

2. *Medicago truncatula*
![](https://github.com/danchurova/Genomes_project/blob/master/Medicago_truncatula/pictures/ATG.jpeg?raw=true)
![](https://github.com/danchurova/Genomes_project/blob/master/Medicago_truncatula/pictures/intron_exon_med.jpeg?raw=true)
![](https://github.com/danchurova/Genomes_project/blob/master/Medicago_truncatula/pictures/promoter-terminator-med.png?raw=true)
