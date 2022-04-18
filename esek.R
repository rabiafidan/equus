setwd("/mnt/NEOGENE1/projects/donkey_2020/modified_MKT")

MKT=read.delim("counts_new2.tsv")

MKT$CDS=as.character(MKT$CDS)

genes=sapply(MKT$CDS,function(x){unlist(strsplit(x,"_"))[2]},USE.NAMES = F)
MKT$genes=genes

uniqMKT2=data.frame()
for (i in 1:length(unique(genes))){
  a=MKT[MKT$genes==unique(genes)[i] ,]
  b=a[a$len==max(a$len),]
  if (nrow(b) != 1){b=b[1,]}
  uniqMKT2=rbind(uniqMKT2,b)
}
old_uniqMKT=uniqMKT
uniqMKT=uniqMKT2

rownames(uniqMKT)=1:nrow(uniqMKT)

mcount=apply(uniqMKT,1,function(x){mean(as.numeric(x[1:10]))})

plot(mcount[c(-6963,-12624)]~uniqMKT$len[c(-6963,-12624)])
abline(lm(mcount[c(-6963,-12624)]~uniqMKT$len[c(-6963,-12624)]),col="red")

plot(mcount~uniqMKT$len)
abline(lm(mcount~uniqMKT$len),col="red")

mean(uniqMKT$Af_As_ns<uniqMKT$Af_As_syn)
mean(uniqMKT$Af_An_ns<uniqMKT$Af_An_syn)
mean(uniqMKT$As_An_ns<uniqMKT$As_An_syn)


MKT_vals=uniqMKT[1:11] +1
MKT_vals$stat=((MKT_vals$As_An_ns/MKT_vals$As_An_syn)/(MKT_vals$As_ns_pol/MKT_vals$As_syn_pol)) +
  ((MKT_vals$Af_An_ns/MKT_vals$Af_An_syn)/(MKT_vals$Afr_ns_pol/MKT_vals$Afr_syn_pol)) - 
  ((MKT_vals$Af_As_ns/MKT_vals$Af_As_syn)/(MKT_vals$As_ns_pol/MKT_vals$As_syn_pol))


uniqMKT$stat=MKT_vals$stat
write.table(uniqMKT[,c(1:11,14,12)],file="stat2.tsv",quote = F,row.names = F,col.names = T,sep = "\t")

plot(uniqMKT$stat~uniqMKT$len)
plot(uniqMKT$stat)
