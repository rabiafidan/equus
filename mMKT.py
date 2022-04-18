from Bio import SeqIO
import os
import pandas as pd
os.chdir("/mnt/NEOGENE1/projects/donkey_2020/modified_MKT")

#######################################   DIVERGENCE   ############################################################

#read individual pop based fasta files and write back CDS-based fasta files
#Asia

Asia_records = list(SeqIO.parse("/mnt/NEOGENE1/projects/donkey_2020/paml/withgff/asian_alts.vcf.fa.out", "fasta"))
for i,_ in enumerate(Asia_records):
    with open(f"data/Asia_single_fasta/{Asia_records[i].description.replace(' ','_')}","w") as out:
        out.write(f">Asia {Asia_records[i].description}\n")
        out.write(str(Asia_records[i].seq)+"\n")

##Africa

Africa_records = list(SeqIO.parse("/mnt/NEOGENE1/projects/donkey_2020/paml/withgff/african_alts.vcf.fa.out", "fasta"))
for i,_ in enumerate(Africa_records):
    with open(f"data/Africa_single_fasta/{Africa_records[i].description.replace(' ','_')}","w") as out:
        out.write(f">Africa {Africa_records[i].description}\n")
        out.write(str(Africa_records[i].seq)+"\n")


#Anatolia
    
Anatolia_records = list(SeqIO.parse("/mnt/NEOGENE1/projects/donkey_2020/paml/withgff/cdh008_alts.vcf.fa.out", "fasta"))
for i,_ in enumerate(Anatolia_records):
    with open(f"data/Anatolia_single_fasta/{Anatolia_records[i].description.replace(' ','_')}","w") as out:
        out.write(f">Anatolia {Anatolia_records[i].description}\n")
        out.write(str(Anatolia_records[i].seq)+"\n")



filenames=[]
for i in Asia_records:
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
    """
    takes codons as strings and returns true if they are synoymous and false if not
    """
    aa_list=[]
    for codon in argv:
        aa_list.append(codon_dict[codon])
    return len(set(aa_list))==1

#count all substitution events between Africa-Anatolia-Asia
Af_As_syn_list=[]
Af_As_ns_list=[]
Af_An_syn_list=[]
Af_An_ns_list=[]
As_An_syn_list=[]
As_An_ns_list=[]
valid_transcript_list=[]
lengths=[]

for filename in filenames:
    Af_As_syn=0
    Af_As_ns=0
    Af_An_syn=0
    Af_An_ns=0
    As_An_syn=0
    As_An_ns=0
    asia_records=list(SeqIO.parse(f"data/Asia_single_fasta/{filename}", "fasta"))
    africa_records=list(SeqIO.parse(f"data/Africa_single_fasta/{filename}", "fasta"))
    anatolia_records=list(SeqIO.parse(f"data/Anatolia_single_fasta/{filename}", "fasta"))
    as_seq=[str(x.seq) for x in asia_records]
    af_seq=[str(x.seq) for x in africa_records] 
    an_seq=[str(x.seq) for x in anatolia_records]    
    #if it is multiple of three and the CDS are not identical
    if len(as_seq[0]) %3==0 and len(set(as_seq+af_seq+an_seq))!=1:
        for codon_idx in range(0,len(as_seq[0]),3):
            as_codons=[a[codon_idx:codon_idx+3] for a in as_seq]
            af_codons=[a[codon_idx:codon_idx+3] for a in af_seq]
            an_codons=[a[codon_idx:codon_idx+3] for a in an_seq] 
            #if there are polymorphism within populations, pass
            if len(set(as_codons))!=1 or len(set(af_codons))!=1 or len(set(an_codons))!=1 or ("N" in as_codons[0]) or ("N" in af_codons[0]) or ("N" in an_codons[0]) :
                pass
            else:
                #Asia-Africa
                if as_codons[0]!=af_codons[0]: #since all within pop codons are the same we can use the first one
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
        lengths.append(len(as_seq[0]))

#store them in a data frame
counts_df=pd.DataFrame(dict(CDS=valid_transcript_list,Af_As_syn=Af_As_syn_list,Af_As_ns=Af_As_ns_list,
                                                Af_An_syn=Af_An_syn_list,Af_An_ns=Af_An_ns_list,
                                                As_An_syn=As_An_syn_list,As_An_ns=As_An_ns_list,len=lengths))



################################## POLYMORPHISM ##############################################################

#for polymorphism info, we add the second set of consensus fastas which contain alternative alleles for heterozgous positions
#we write them back as CDS based, together with the ref.

Asian_poly = list(SeqIO.parse("data/cufflinks_pop_fastas/Asian_biallelic_polymorphic.fa", "fasta"))
African_poly = list(SeqIO.parse("data/cufflinks_pop_fastas/African_biallelic_polymorphic.fa", "fasta"))
Anatolian_poly = list(SeqIO.parse("data/cufflinks_pop_fastas/Anatolian_biallelic_polymorphic.fa", "fasta"))
ref = list(SeqIO.parse("/mnt/NEOGENE1/projects/donkey_2020/paml/withgff/horse.fa.out", "fasta"))

#Asia
for i,_ in enumerate(Asian_poly):
    with open(f"data/Asia_single_fasta/poly_{Asian_poly[i].description.replace(' ','_')}","w") as out:
        out.write(f">Asia_poly {Asian_poly[i].description}\n")
        out.write(str(Asian_poly[i].seq)+"\n")
        out.write(f">ref {Asian_poly[i].description}\n")
        out.write(str(ref[i].seq))


#Africa

for i,_ in enumerate(African_poly):
    with open(f"data/Africa_single_fasta/poly_{African_poly[i].description.replace(' ','_')}","w") as out:
        out.write(f">Africa_poly {African_poly[i].description}\n")
        out.write(str(African_poly[i].seq)+"\n")
        out.write(f">ref {African_poly[i].description}\n")
        out.write(str(ref[i].seq))



def syn_list(codon_list):
    """
    takes a list of codons and returns true if they are synonymous; false if not
    """
    aa_list=[]
    for codon in codon_list:
        aa_list.append(codon_dict[codon])
    return len(set(aa_list))==1

def pols(pop_name,transcript_ls):
    """
    takes a population name and a transcript name list
    returns 2 lists. One has synonymous polymorphism counts; one has non-synonymous polymorphism counts
    for the given list of CDSs
    """
    syn_ls=[]
    ns_ls=[]
    for gene in transcript_ls:
        syn=0
        ns=0
        records=list(SeqIO.parse(f"data/{pop_name}_single_fasta/poly_{gene}", "fasta"))
        seq_ls=[str(a.seq) for a in records]
        if len(set(seq_ls))!=1: #if all sequences are identical, counts will remain 0
            for codon_idx in range(0,len(seq_ls[0]),3):
                codon_ls=[a[codon_idx:codon_idx+3] for a in seq_ls]
                if len(set(codon_ls))>2: #there were none in our dataset
                    print("more than 2 codons!")
                if "N" in codon_ls[0]:
                    pass
                else:
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
counts_df['As_ns_pol']=As_ns_pol
counts_df['Afr_syn_pol']=Afr_syn_pol
counts_df['Afr_ns_pol']=Afr_ns_pol
counts2_df=counts_df[["Af_As_syn","Af_As_ns","Af_An_syn","Af_An_ns","As_An_syn","As_An_ns",'As_syn_pol','As_ns_pol','Afr_syn_pol','Afr_ns_pol',"len","CDS"]]

counts2_df.to_csv("counts_new2.tsv",sep="\t",index=False)

