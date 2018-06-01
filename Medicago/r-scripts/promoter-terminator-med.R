snps=read.table('SNPS-qual100.txt',header=T, sep='\t', stringsAsFactors = F)
snps = dplyr::select(snps, -global)
colnames(snps)[c(2,3)]=c('snp_position','snp_chrom')
unique(snps$snp_chrom)
snps$snp_chrom[snps$snp_chrom == "chr1"] = 1
snps$snp_chrom[snps$snp_chrom == "chr2"] = 2
snps$snp_chrom[snps$snp_chrom == "chr3"] = 3
snps$snp_chrom[snps$snp_chrom == "chr4"] = 4
snps$snp_chrom[snps$snp_chrom == "chr5"] = 5
chrom1_2 = dplyr::filter(snps, snp_chrom == 1)


anno = read.table('files/mt_anno_fin.txt', sep='\t',head=T,stringsAsFactors = F)
anno = dplyr::select(anno, c(Chr,Frag,Start,Start.1,Strand))
a = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5')
anno_mod = dplyr::filter(anno, Chr %in% a)
anno_mod$Chr[anno_mod$Chr == 'chr1'] = 'Chr1'
anno_mod$Chr[anno_mod$Chr == 'chr2'] = 'Chr2'
anno_mod$Chr[anno_mod$Chr == 'chr3'] = 'Chr3'
anno_mod$Chr[anno_mod$Chr == 'chr4'] = 'Chr4'
anno_mod$Chr[anno_mod$Chr == 'chr5'] = 'Chr5'
anno = anno_mod
###
anno$End = anno$Start.1
###
anno$len=anno$End-anno$Start


promoters=read.table("files/Promoters_final_mt.txt",sep="\t",header=T, stringsAsFactors = F)
a = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5')
prom_mod = dplyr::filter(promoters, CHROM %in% a)
prom_mod$Chr = prom_mod$CHROM
prom_mod$Chr[prom_mod$Chr == 'chr1'] = 'Chr1'
prom_mod$Chr[prom_mod$Chr == 'chr2'] = 'Chr2'
prom_mod$Chr[prom_mod$Chr == 'chr3'] = 'Chr3'
prom_mod$Chr[prom_mod$Chr == 'chr4'] = 'Chr4'
prom_mod$Chr[prom_mod$Chr == 'chr5'] = 'Chr5'

prom_mod = dplyr::select(prom_mod, c(Chr,MODEL,TSS,STRAND))

promoters = prom_mod



terminator=read.table("files/4tss_mt_mrnas.txt",sep="\t",header=T, stringsAsFactors = F)
term_mod = dplyr::filter(terminator, CHROM %in% a)
term_mod$Chr = term_mod$CHROM
term_mod$Chr[term_mod$Chr == 'chr1'] = 'Chr1'
term_mod$Chr[term_mod$Chr == 'chr2'] = 'Chr2'
term_mod$Chr[term_mod$Chr == 'chr3'] = 'Chr3'
term_mod$Chr[term_mod$Chr == 'chr4'] = 'Chr4'
term_mod$Chr[term_mod$Chr == 'chr5'] = 'Chr5'
term_mod = dplyr::select(term_mod, c(Chr,FRAG_START,FRAG_STOP,STRAND))
terminator = term_mod

#######################
##calculate promoters
#######################
promoters$Start=promoters$TSS-1000
promoters$End=promoters$TSS+1000
prom=promoters
#colnames(prom)[1]='Chr'


snp_custom_annotation(chrom1_2,prom)

minus=subset(snp_annotated_list,STRAND=="-")
plus=subset(snp_annotated_list,STRAND=="+")

#обрабатываем плюс-цепь
temp=plus
temp=cbind(temp$Start,temp$snp_position)
temp=as.data.frame(temp)
temp$V1=temp$V2-temp$V1
temp$V2=1
temp=na.omit(temp)
temp=temp[order(temp),]
temp_graph_df_plus=data.frame()
temp_df_plus=data.frame()
for (i in 0:2000){
  temp_df=subset(temp,V1==i)
  temp_df[1,2]=sum(temp_df$V2)
  temp_df=temp_df[1,]
  temp_graph_df_plus=rbind(temp_graph_df_plus,temp_df)}
plot(temp_graph_df_plus)
temp_graph_df_plus$V1=temp_graph_df_plus$V1-1000
colnames(temp_graph_df_plus)=c('snp_position','plus')


#обрабатываем минус-цепь
temp=minus
temp=cbind(temp$Start,temp$snp_position)
temp=as.data.frame(temp)
temp$V1=temp$V2-temp$V1
temp$V2=1
temp=na.omit(temp)
temp=temp[order(temp),]
#x=data.frame()
temp_graph_df_minus=data.frame()
temp_df=data.frame()
for (i in 0:2000){
  temp_df=subset(temp,V1==i)
  temp_df[1,2]=sum(temp_df$V2)
  temp_df=temp_df[1,]
  temp_graph_df_minus=rbind(temp_graph_df_minus,temp_df)}
plot(temp_graph_df_minus)
temp_graph_df_minus$V2=rev(temp_graph_df_minus$V2) #reverse minus-strand
temp_graph_df_minus$V1=temp_graph_df_minus$V1-1000
colnames(temp_graph_df_minus)=c("V1","V2")


graph_df=temp_graph_df_plus
graph_df[,2]=graph_df[,2]+temp_graph_df_minus[,2]
graph_df[,2]=graph_df[,2]/20367*1000

promoter_plot=graph_df #create plot for transcription promoter



#######################
##calculate terminators
#######################
plus=subset(terminator,STRAND=="+")
minus=subset(terminator,STRAND=="-")
minus$Start=minus$FRAG_START-1000
minus$End=minus$FRAG_START+1000
plus$Start=plus$FRAG_STOP-1000
plus$End=plus$FRAG_STOP+1000
ter=rbind(plus,minus)
colnames(ter)[1]='Chr'


snp_custom_annotation(chrom1_2,ter)

minus=subset(snp_annotated_list,STRAND=="-")
plus=subset(snp_annotated_list,STRAND=="+")

#обрабатываем плюс-цепь
temp=plus
temp=cbind(temp$Start,temp$snp_position)
temp=as.data.frame(temp)
temp$V1=temp$V2-temp$V1
temp$V2=1
temp=na.omit(temp)
temp=temp[order(temp),]
temp_graph_df_plus=data.frame()
temp_df_plus=data.frame()
for (i in 0:2000){
  temp_df=subset(temp,V1==i)
  temp_df[1,2]=sum(temp_df$V2)
  temp_df=temp_df[1,]
  temp_graph_df_plus=rbind(temp_graph_df_plus,temp_df)}
plot(temp_graph_df_plus)
temp_graph_df_plus$V1=temp_graph_df_plus$V1-1000
colnames(temp_graph_df_plus)=c('snp_position','plus')


#обрабатываем минус-цепь
temp=minus
temp=cbind(temp$Start,temp$snp_position)
temp=as.data.frame(temp)
temp$V1=temp$V2-temp$V1
temp$V2=1
temp=na.omit(temp)
temp=temp[order(temp),]
#x=data.frame()
temp_graph_df_minus=data.frame()
temp_df=data.frame()
for (i in 0:2000){
  temp_df=subset(temp,V1==i)
  temp_df[1,2]=sum(temp_df$V2)
  temp_df=temp_df[1,]
  temp_graph_df_minus=rbind(temp_graph_df_minus,temp_df)}
plot(temp_graph_df_minus)
temp_graph_df_minus$V2=rev(temp_graph_df_minus$V2) #reverse minus-strand
temp_graph_df_minus$V1=temp_graph_df_minus$V1-1000
colnames(temp_graph_df_minus)=c("V1","V2")


graph_df=temp_graph_df_plus
graph_df[,2]=graph_df[,2]+temp_graph_df_minus[,2]
graph_df[,2]=graph_df[,2]/20367*1000

terminator_plot=graph_df #create plot for transcription terminator

#############
##Build a graph
#############


#tiff("prom_trancr_ter_coverage_filtered1.tiff",height = 18, width = 12, units = 'in',res=600)
dev.off()
par(mfrow=c(1,2))
par(mar = c(5,4,4,4) + 0.1)
plot(promoter_plot,ylab="Number of SNPs per 1000 genes",
     xlab="Distance from TSS")
y=(density(subset(subset(anno,len<1000),Frag=='five_prime_UTR')$len))$y*2000+10
x=(density(subset(subset(anno,len<1000),Frag=='five_prime_UTR')$len))$x
x=x[19:500]
y=y[19:500]
lines(x,y,col='blue',lwd=3)
axis(4,at=seq(0,0.005,by=0.001)*1300+3,label=seq(0,0.005,by=0.001))
mtext("Distribution of 5'UTR lengths", side=4, line=2)
lines(predict(smooth.spline(promoter_plot,df=10)),col='red',lwd=3)
#legend("topright", c("Average number of SNPs","Distribution of 5'UTR lengths"),fill=c("red","blue"))
#mtext("A)", adj=0, line=1, font=2,cex=1.5)

plot(terminator_plot,ylab="Number of SNPs per 1000 genes"
     ,xlab="Distance from TTS")
y=(density(subset(subset(anno,len<1000),Frag=='three_prime_UTR')$len))$y*2500+10
x=(density(subset(subset(anno,len<1000),Frag=='three_prime_UTR')$len))$x
x=x[28:500]
y=y[28:500]
lines(predict(smooth.spline(terminator_plot,df=10)),col='red',lwd=3)
lines(predict(smooth.spline(x*(-1),y,df=7)),col='blue',lwd=3)
axis(4,at=seq(0,0.005,by=0.001)*1300+3,label=seq(0,0.005,by=0.001))
mtext("Distribution of 3'UTR lengths", side=4, line=2)
#legend("topright", c("Average number of SNPs","Distribution of 3'UTR lengths"),fill=c("red","blue"))
#mtext("B)", adj=0, line=1, font=2,cex=1.5)
