# Genomes_project
## Statistical analysis of genomes
### Goals and objectives
Find correlation between sequence features and functional regions in different genomes

1. Plot sequence features such as TFBS, SNPs, methylation, RNA-seq coverage
2. Map it on promoter regions
3. Find correlation
4. Consider outcomes for promoter prediction for complex and not annotated genomes


Data:
* *Arabidopsis thaliana*
1. Graphs for *Arabidopsis thaliana* as done for *Oryza sativa* [[1](https://www.ncbi.nlm.nih.gov/pubmed/27774999)]
* gff3 file (**TAIR10_GFF3_genes.gff**) was taken [here](https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FTAIR10_genome_release%2FTAIR10_gff3)
2. [vcf file 1001 genome TAIR]()
3. *A. thaliana* reference genome (TAIR10_toplevel) was taken from ftp://ftp.ensemblgenomes.org/pub/plants/release-39/fasta/arabidopsis_thaliana/dna/

* *Medicago truncata*
1. annotation (.gff) and assemly (.fasta) from http://www.medicagogenome.org/downloads


Scripts:
1. get_ATGs.py
2. get_4tss.py
4. get_promoters.py
3. get_anno_fin.py
 Usage:
 1. to create file with ATGs: ``` python3 get_ATGs.py annotation.gff ```
 2. to create file with tss: ``` python3 get_4tss.py annotation.gff ```
 3. to create files with promoter regions (.bed + .txt): ``` python3 get_promoters.py 4tss.txt ```
