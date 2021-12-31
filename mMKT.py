from Bio import SeqIO
import os
import pandas as pd
import glob
os.chdir("/mnt/NEOGENE1/projects/donkey_2020/modified_MKT")

Aw2_records = list(SeqIO.parse("data/cufflinks_fastas/Aw2_ref.fa", "fasta"))
Aw3_records = list(SeqIO.parse("data/cufflinks_fastas/Aw3_ref.fa", "fasta"))
BayanNur_records = list(SeqIO.parse("data/cufflinks_fastas/BayanNur_ref.fa", "fasta"))
Kiang_records = list(SeqIO.parse("data/cufflinks_fastas/kiang_ref.fa", "fasta"))
Kia_records = list(SeqIO.parse("data/cufflinks_fastas/Kia2_ref.fa", "fasta"))
hem_records = list(SeqIO.parse("data/cufflinks_fastas/hemionus_ref.fa", "fasta"))


for i,_ in enumerate(Aw2_records):
    with open(f"data/Asia_homalt_fasta/{Aw2_records[i].description.replace(' ','_')}","w") as out:
        out.write(f">Aw2 {Aw2_records[i].description}\n")
        out.write(str(Aw2_records[i].seq)+"\n")
        out.write(f">Aw3 {Aw2_records[i].description}\n")
        out.write(str(Aw3_records[i].seq)+"\n")
        out.write(f">BayanNur {Aw2_records[i].description}\n")
        out.write(str(BayanNur_records[i].seq)+"\n")
        out.write(f">Kiang {Aw2_records[i].description}\n")
        out.write(str(Kiang_records[i].seq)+"\n")
        out.write(f">Kia {Aw2_records[i].description}\n")
        out.write(str(Kia_records[i].seq)+"\n")
        out.write(f">hemionus {Aw2_records[i].description}\n")
        out.write(str(hem_records[i].seq))



Chby1_records = list(SeqIO.parse("data/cufflinks_fastas/Ch-by1_ref.fa", "fasta"))
Eg1_records = list(SeqIO.parse("data/cufflinks_fastas/Eg-1_ref.fa", "fasta"))
Ir5_records = list(SeqIO.parse("data/cufflinks_fastas/Ir-5_ref.fa", "fasta"))
Ke14_records = list(SeqIO.parse("data/cufflinks_fastas/Ke-14_ref.fa", "fasta"))
Ky5_records = list(SeqIO.parse("data/cufflinks_fastas/Ky-5_ref.fa", "fasta"))
Sp5_records = list(SeqIO.parse("data/cufflinks_fastas/Sp-5_ref.fa", "fasta"))
somalicus_records = list(SeqIO.parse("data/cufflinks_fastas/somalicus_ref.fa", "fasta"))
    
for i,_ in enumerate(Chby1_records):
    with open(f"data/Africa_homalt_fasta/{Chby1_records[i].description.replace(' ','_')}","w") as out:
        out.write(f">Chby1 {Chby1_records[i].description}\n")
        out.write(str(Chby1_records[i].seq)+"\n")
        out.write(f">Eg1 {Chby1_records[i].description}\n")
        out.write(str(Eg1_records[i].seq)+"\n")
        out.write(f">Ir5 {Chby1_records[i].description}\n")
        out.write(str(Ir5_records[i].seq)+"\n")
        out.write(f">Ke14 {Chby1_records[i].description}\n")
        out.write(str(Ke14_records[i].seq)+"\n")
        out.write(f">Ky5 {Chby1_records[i].description}\n")
        out.write(str(Ky5_records[i].seq)+"\n")
        out.write(f">Sp5 {Chby1_records[i].description}\n")
        out.write(str(Sp5_records[i].seq)+"\n")
        out.write(f">somalicus {Chby1_records[i].description}\n")
        out.write(str(somalicus_records[i].seq))

Anatolia_records = list(SeqIO.parse("data/cufflinks_fastas/cdh008_ref.fa", "fasta"))
for i,_ in enumerate(Anatolia_records):
    with open(f"data/Anatolia_homalt_fasta/{Anatolia_records[i].description.replace(' ','_')}","w") as out:
        out.write(f">cdh008 {Anatolia_records[i].description}\n")
        out.write(str(Anatolia_records[i].seq))
    
filenames=[]
for i in Aw2_records:
    filenames.append(i.description.replace(' ','_'))



codon_dict=dict(TTT="Phe",TTC="Phe",TTA="Leu",TTG="Leu",
                CTT="Leu",CTC="Leu",CTA="Leu",CTG="Leu",
                ATT="Ile",ATC="Ile",ATA="Ile",ATG="Met",
                GTT="VAL",GTC="VAL",GTA="VAL",GTG="VAL",
                TCT="Ser",TCC="Ser",TCA="Ser",TCG="Ser",
                CCT="Pro",CCC="Pro",CCA="Pro",CCG="Pro",
                ACT="Thr",ACC="Thr",ACA="Thr",ACG="Thr",
                GCT="Ala",GCC="Ala",GCA="Ala",GCG="Ala",
                TAT="Tyr",TAC="Tyr",TAA="Stop",TAG="Stop",
                CAT="His",CAC="His",CAA="Gln",CAG="Gln",
                AAT="Asn",AAC="Asn",AAA="Lys",AAG="Lys",
                GAT="Asp",GAC="Asp",GAA="Glu",GAG="Glu",
                TGT="Cys",TGC="Cys",TGA="Stop",TGG="Trp",
                CGT="Arg",CGC="Arg",CGA="Arg",CGG="Arg",
                AGT="Ser",AGC="Ser",AGA="Arg",AGG="Arg",
                GGT="Gly",GGC="Gly",GGA="Gly",GGG="Gly"
                )

def syn(*argv):
    aa_list=[]
    for codon in argv:
        aa_list.append(codon_dict[codon])
    return len(set(aa_list))==1

#subs
Af_As_syn_list=[]
Af_As_ns_list=[]
Af_An_syn_list=[]
Af_An_ns_list=[]
As_An_syn_list=[]
As_An_ns_list=[]
valid_transcript_list=[]

for filename in filenames:
    Af_As_syn=0
    Af_As_ns=0
    Af_An_syn=0
    Af_An_ns=0
    As_An_syn=0
    As_An_ns=0
    asia_records=list(SeqIO.parse(f"data/Asia_homalt_fasta/{filename}", "fasta"))
    africa_records=list(SeqIO.parse(f"data/Asia_homalt_fasta/{filename}", "fasta"))
    anatolia_records=list(SeqIO.parse(f"data/Anatolia_homalt_fasta/{filename}", "fasta"))
    as_seq=[str(x.seq) for x in asia_records]
    af_seq=[str(x.seq) for x in africa_records] 
    an_seq=[str(x.seq) for x in anatolia_records]    
    #if it is multiple of three and the CDS are not identical
    if len(as_seq[0]) %3==0 and len(set(as_seq+af_seq+an_seq))!=1:
        for codon_idx in range(0,len(as_seq[0])//3,3):
            as_codons=[a[codon_idx:codon_idx+3] for a in as_seq]
            af_codons=[a[codon_idx:codon_idx+3] for a in af_seq]
            an_codons=[a[codon_idx:codon_idx+3] for a in an_seq] 
            #if there are polymorphism within populations, pass
            if len(set(as_codons))!=1 or len(set(af_codons))!=1 or len(set(an_codons))!=1:
                pass
            else:
                #Asia-Africa
                if as_codons[0]!=af_codons[0]:
                    if syn(as_codons[0],af_codons[0]):
                        Af_As_syn+=1
                    else:
                        Af_As_ns+=1
                #Africa-Anatolia
                if an_codons[0]!=af_codons[0]:
                    if syn(an_codons[0],af_codons[0]):
                        Af_An_syn+=1
                    else:
                        Af_An_ns+=1
                #Anatolia-Asia
                if an_codons[0]!=as_codons[0]:
                    if syn(an_codons[0],as_codons[0]):
                        As_An_syn+=1
                    else:
                        As_An_ns+=1
 
        Af_As_syn_list.append(Af_As_syn)
        Af_As_ns_list.append(Af_As_ns) 
        Af_An_syn_list.append(Af_An_syn)
        Af_An_ns_list.append(Af_An_ns)
        As_An_syn_list.append(As_An_syn)
        As_An_ns_list.append(As_An_ns)
        valid_transcript_list.append(filename)     


counts_df=pd.DataFrame(dict(CDS=valid_transcript_list,Af_As_syn=Af_As_syn_list,Af_As_ns=Af_As_ns_list,
                                                Af_An_syn=Af_An_syn_list,Af_An_ns=Af_An_ns_list,
                                                As_An_syn=As_An_syn_list,As_An_ns=As_An_ns_list))



Aw2_records_alt = list(SeqIO.parse("data/cufflinks_fastas/Aw2_alt.fa", "fasta"))
Aw3_records_alt = list(SeqIO.parse("data/cufflinks_fastas/Aw3_alt.fa", "fasta"))
BayanNur_records_alt = list(SeqIO.parse("data/cufflinks_fastas/BayanNur_alt.fa", "fasta"))
Kiang_records_alt = list(SeqIO.parse("data/cufflinks_fastas/kiang_alt.fa", "fasta"))
Kia_records_alt = list(SeqIO.parse("data/cufflinks_fastas/Kia2_alt.fa", "fasta"))
hem_records_alt = list(SeqIO.parse("data/cufflinks_fastas/hemionus_alt.fa", "fasta"))


for i,r in enumerate(Aw2_records):
    with open(f"data/Asia_polymorphic_fasta/{r.description.replace(' ','_')}","w") as out:
        out.write(f">Aw2 {Aw2_records[i].description}\n")
        out.write(str(Aw2_records[i].seq)+"\n")
        out.write(f">Aw2_alt {Aw2_records[i].description}\n")
        out.write(str(Aw2_records_alt[i].seq)+"\n")
        out.write(f">Aw3 {Aw2_records[i].description}\n")
        out.write(str(Aw3_records[i].seq)+"\n")
        out.write(f">Aw3_alt {Aw2_records[i].description}\n")
        out.write(str(Aw3_records_alt[i].seq)+"\n")
        out.write(f">BayanNur {Aw2_records[i].description}\n")
        out.write(str(BayanNur_records[i].seq)+"\n")
        out.write(f">BayanNur_alt {Aw2_records[i].description}\n")
        out.write(str(BayanNur_records_alt[i].seq)+"\n")
        out.write(f">Kiang {Aw2_records[i].description}\n")
        out.write(str(Kiang_records[i].seq)+"\n")
        out.write(f">Kiang_alt {Aw2_records[i].description}\n")
        out.write(str(Kiang_records_alt[i].seq)+"\n")
        out.write(f">Kia {Aw2_records[i].description}\n")
        out.write(str(Kia_records[i].seq)+"\n")
        out.write(f">Kia_alt {Aw2_records[i].description}\n")
        out.write(str(Kia_records_alt[i].seq)+"\n")
        out.write(f">hemionus {Aw2_records[i].description}\n")
        out.write(str(hem_records[i].seq))
        out.write(f">hemionus_alt {Aw2_records[i].description}\n")
        out.write(str(hem_records_alt[i].seq))


Chby1_records_alt = list(SeqIO.parse("data/cufflinks_fastas/Ch-by1_alt.fa", "fasta"))
Eg1_records_alt = list(SeqIO.parse("data/cufflinks_fastas/Eg-1_alt.fa", "fasta"))
Ir5_records_alt = list(SeqIO.parse("data/cufflinks_fastas/Ir-5_alt.fa", "fasta"))
Ke14_records_alt = list(SeqIO.parse("data/cufflinks_fastas/Ke-14_alt.fa", "fasta"))
Ky5_records_alt = list(SeqIO.parse("data/cufflinks_fastas/Ky-5_alt.fa", "fasta"))
Sp5_records_alt = list(SeqIO.parse("data/cufflinks_fastas/Sp-5_alt.fa", "fasta"))
somalicus_records_alt = list(SeqIO.parse("data/cufflinks_fastas/somalicus_alt.fa", "fasta"))
    
for i,_ in enumerate(Chby1_records):
    with open(f"data/Africa_polymorphic_fasta/{Chby1_records[i].description.replace(' ','_')}","w") as out:
        out.write(f">Chby1 {Chby1_records[i].description}\n")
        out.write(str(Chby1_records[i].seq)+"\n")
        out.write(f">Chby1_alt {Chby1_records[i].description}\n")
        out.write(str(Chby1_records_alt[i].seq)+"\n")
        out.write(f">Eg1 {Chby1_records[i].description}\n")
        out.write(str(Eg1_records[i].seq)+"\n")
        out.write(f">Eg1_alt {Chby1_records[i].description}\n")
        out.write(str(Eg1_records_alt[i].seq)+"\n")
        out.write(f">Ir5 {Chby1_records[i].description}\n")
        out.write(str(Ir5_records[i].seq)+"\n")
        out.write(f">Ir5_alt {Chby1_records[i].description}\n")
        out.write(str(Ir5_records_alt[i].seq)+"\n")
        out.write(f">Ke14 {Chby1_records[i].description}\n")
        out.write(str(Ke14_records[i].seq)+"\n")
        out.write(f">Ke14_alt {Chby1_records[i].description}\n")
        out.write(str(Ke14_records_alt[i].seq)+"\n")
        out.write(f">Ky5 {Chby1_records[i].description}\n")
        out.write(str(Ky5_records[i].seq)+"\n")
        out.write(f">Ky5_alt {Chby1_records[i].description}\n")
        out.write(str(Ky5_records_alt[i].seq)+"\n")
        out.write(f">Sp5 {Chby1_records[i].description}\n")
        out.write(str(Sp5_records[i].seq)+"\n")
        out.write(f">Sp5_alt {Chby1_records[i].description}\n")
        out.write(str(Sp5_records_alt[i].seq)+"\n")
        out.write(f">somalicus {Chby1_records[i].description}\n")
        out.write(str(somalicus_records[i].seq))
        out.write(f">somalicus_alt {Chby1_records[i].description}\n")
        out.write(str(somalicus_records_alt[i].seq))

def syn_list(codon_list):
    aa_list=[]
    for codon in codon_list:
        aa_list.append(codon_dict[codon])
    return len(set(aa_list))==1

def pols(pop_name,transcript_ls):
    syn_ls=[]
    ns_ls=[]
    for gene in transcript_ls:
        syn=0
        ns=0
        records=list(SeqIO.parse(f"data/{pop_name}_polymorphic_fasta/{gene}", "fasta"))
        seq_ls=[str(a.seq) for a in records]
        if len(set(seq_ls))!=1:
            for codon_idx in range(0,len(seq_ls[0])//3,3):
                codon_ls=[a[codon_idx:codon_idx+3] for a in seq_ls]
                if len(set(codon_ls))>2:
                    print("more than 2 codons!")
                if len(set(codon_ls))>1:
                    if syn_list(list(set(codon_ls))):
                        syn+=1
                    else:
                        ns+=1
        syn_ls.append(syn)
        ns_ls.append(ns)
    return [syn_ls,ns_ls]

As_syn_pol,As_ns_pol=pols("Asia",valid_transcript_list)
Afr_syn_pol,Afr_ns_pol=pols("Africa",valid_transcript_list)

counts_df['As_syn_pol']=As_syn_pol
counts_df['As_ns_pol']=As_syn_pol
counts_df['Afr_syn_pol']=As_syn_pol
counts_df['Afr_ns_pol']=As_syn_pol
counts2_df=counts_df[["Af_As_syn","Af_As_ns","Af_An_syn","Af_An_ns","As_An_syn","As_An_ns",'As_syn_pol','As_ns_pol','Afr_syn_pol','Afr_ns_pol',"CDS"]]

counts2_df.to_csv("counts.tsv",sep="\t",index=False)

    
