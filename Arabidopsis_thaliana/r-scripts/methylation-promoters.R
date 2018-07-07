meth_arab_CH = read.table('WT sperm BS-seq-arab-CHH.gff',sep='\t',head=F,stringsAsFactors = F)

meth_arab_CH$V1[meth_arab_CH$V1 == "MITOCHONDRIA"] ="M"
meth_arab_CH$V1[meth_arab_CH$V1 == "CHR1"] = 1
meth_arab_CH$V1[meth_arab_CH$V1 == "CHR2"] = 2
meth_arab_CH$V1[meth_arab_CH$V1 == "CHR3"] = 3
meth_arab_CH$V1[meth_arab_CH$V1 == "CHR4"] = 4
meth_arab_CH$V1[meth_arab_CH$V1 == "CHR5"] = 5
meth_arab_CH = dplyr::filter(meth_arab_CH, V1!='CHLOROPLAST')
meth_arab_CH$V6 = as.numeric(meth_arab_CH$V6)
meth_arab_CH$V7 = meth_arab_CH$V4
colnames(meth_arab_CH)[c(1,7)] = c('snp_chrom','snp_position')
meth_arab_CH$snp_position[is.na(meth_arab_CH$V6)] = NA
meth_arab_CH2 = dplyr::filter(meth_arab_CH, !is.na(snp_position))





promoters = read.table('E:/PROJECT-SPRING/arab/arab_files_N/Promoters_final_arab.txt',sep='\t',head=T)
promoters$Start = promoters$TSS-1000
promoters$End = promoters$TSS+1000
colnames(promoters)[1] = 'Chr'

mean(meth_arab_CH$V6, na.rm=T)

snp_custom_annotation(meth_arab_CH2, promoters)

minus=subset(snp_annotated_list,STRAND=="-")
plus=subset(snp_annotated_list,STRAND=="+")

#обрабатываем плюс-цепь
temp=plus
temp=cbind(temp$Start,temp$snp_position)
temp=as.data.frame(temp)
temp$V1=temp$V2-temp$V1
temp$V2=plus$V6
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
temp$V2=minus$V6
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

#promoter_plot_CHG=graph_df
#promoter_plot_CG=graph_df
promoter_plot_CHH=graph_df

promoters_all_meth = cbind(pos = promoter_plot_CG$snp_position, CHG = promoter_plot_CHG$plus,
                           CG = promoter_plot_CG$plus, CHH = promoter_plot_CHH$plus)
promoters_all_meth = as.data.frame(promoters_all_meth)
dev.off()
#par(mfrow=c(1,2))
#par(mar = c(5,4,4,4) + 0.1)
plot(promoter_plot_CHG,ylab="Methylation profile",
     xlab="Distance from TSS", col = 'red', ylim = c(0,5), main = 'Methylation profile for sperm cell promoters')
points(promoter_plot_CG, col='green')
points(promoter_plot_CHH, col='blue')
legend("topleft", c("CHG","CG","CHH"),fill=c("red","green","blue"))
