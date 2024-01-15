setwd("/mnt/NEOGENE1/projects/donkey_2020/modified_MKT")

#MKT=read.delim("counts_new2.tsv")
MKT=read.delim("counts_reviewed.tsv")

#data manipulation

MKT$CDS=as.character(MKT$CDS)

genes=sapply(MKT$CDS,function(x){unlist(strsplit(x,"_"))[2]},USE.NAMES = F)
MKT$genes=genes

# If there are multiple transcripts of a gene, take the longest
#If there are multiple longest ones, thake the first one.

uniqMKT2=data.frame()
for (i in 1:length(unique(genes))){
  a=MKT[MKT$genes==unique(genes)[i] ,]
  b=a[a$len==max(a$len),]
  if (nrow(b) != 1){b=b[1,]}
  uniqMKT2=rbind(uniqMKT2,b)
}
uniqMKT=uniqMKT2
rownames(uniqMKT)=1:nrow(uniqMKT)

#quality control

mcount=apply(uniqMKT,1,function(x){mean(as.numeric(x[1:10]))})

plot(mcount[c(-6963,-12624)]~uniqMKT$len[c(-6963,-12624)])
abline(lm(mcount[c(-6963,-12624)]~uniqMKT$len[c(-6963,-12624)]),col="red")

plot(mcount~uniqMKT$len)
abline(lm(mcount~uniqMKT$len),col="red")

mean(uniqMKT$Af_As_ns<=uniqMKT$Af_As_syn)
mean(uniqMKT$Af_An_ns<=uniqMKT$Af_An_syn)
mean(uniqMKT$As_An_ns<=uniqMKT$As_An_syn)

#add 1 to whole dataset to avoid zero divisions
MKT_vals=uniqMKT[1:11] +1

#Calculate the pairwise MKT statistic
MKT_vals$stat=((MKT_vals$As_An_ns/MKT_vals$As_An_syn)/(MKT_vals$As_ns_pol/MKT_vals$As_syn_pol)) +
  ((MKT_vals$Af_An_ns/MKT_vals$Af_An_syn)/(MKT_vals$Afr_ns_pol/MKT_vals$Afr_syn_pol)) - 
  ((MKT_vals$Af_As_ns/MKT_vals$Af_As_syn)/(MKT_vals$Afr_ns_pol/MKT_vals$Afr_syn_pol))


uniqMKT$stat=MKT_vals$stat
write.table(uniqMKT[,c(1:11,14,12)],file="stat_reviewed.tsv",quote = F,row.names = F,col.names = T,sep = "\t")

plot(uniqMKT$stat~uniqMKT$len)
#order the table
oMKT=uniqMKT[order(uniqMKT$stat),]
#add a z-score
oMKT$Zscore=(oMKT$stat-mean(oMKT$stat))/sd(oMKT$stat)

plot(oMKT$Zscore)
abline(h=3,col="red",lwd=2)
hist(oMKT$Zscore)
abline(v=3,col="red",lwd=2)
#take the genes with z-score bigger than 3
z3=which(oMKT$Zscore>3)[1]

oMKT$pval=2*pnorm(-abs(oMKT$Zscore))
oMKT$padj=p.adjust(oMKT$pval,method = "BH")
sum(oMKT[,"padj"]<0.01)
sum(oMKT[,"pval"]<0.01)
length(z3:length(oMKT$CDS))
sum(oMKT$Zscore>3&oMKT$padj<0.01)

condition=oMKT$Zscore>3&oMKT$padj<0.01
write.table(oMKT[z3:length(oMKT$CDS),"CDS"],file="high_stat_genes_reviewed.txt",quote = F,row.names = F,col.names = F,sep = "\t")
write.table(oMKT[condition,"CDS"],file="multiple_testing_passed_reviewed.txt",quote = F,row.names = F,col.names = F,sep = "\t")



common=matrix(c(7,174,107,15696),nrow = 2)
row.names(common)=c("PMKsig","PMKnon-sig")
colnames(common)=c("PAMLsig","PAMLnon-sig")

common

fisher.test(common)
chisq.test(common)$expected

par(mfcol=c(1,2))
hist(oMKT$stat,col="darkred",main="PMK")
hist(oMKT$Zscore,col="indianred4",xlab="PMK Z-score (n=17327)",main = "",breaks = 70)
abline(v=3,lwd=2, col="cornflowerblue",lty=2)

