configfile: "config.yaml"

localrules: all
"""
(base) rabia@host3:/mnt/NEOGENE1/projects/donkey_2020/dataset/denovoSnpCall$ bcftools +fill-tags merge-vcf-all.vcf.gz -Oz -o merge-vcf-all_with_NS_tags.vcf.gz -- -t NS -S ../../modified_MKT/sample_list.txt
(base) rabia@host3:/mnt/NEOGENE1/projects/donkey_2020/dataset/denovoSnpCall$ bcftools view -m2 -M2 -i 'INFO/NS_Asian>0 & INFO/NS_African>0 & INFO/NS_Anatolian>0' merge-vcf-all_with_NS_tags.vcf.gz -Oz -o merge-vcf-all_with_NS_tags_filtered.vcf.gz 


"""

rule all:
	input:
		expand("data/cufflinks_fastas/{sample}_{H}.fa",sample=config["All"],H=["ref","alt"]),
		#"data/cufflinks_fastas/EquCab2.fa",
		#"data/pop/African_biallelic.vcf.gz",
		#"data/pop/Asian_biallelic.vcf.gz",
		#"data/pop/Anatolian_biallelic.vcf.gz",
		"data/pop/African_biallelic_polymorphic.vcf.gz",
		"data/pop/Asian_biallelic_polymorphic.vcf.gz",
		"data/pop/Anatolian_biallelic_polymorphic.vcf.gz",
		"data/pop/African_biallelic_polymorphic_alt.fa",
		"data/pop/Asian_biallelic_polymorphic_alt.fa",
		"data/pop/Anatolian_biallelic_polymorphic_alt.fa",
		expand("data/cufflinks_pop_fastas/{pops}_biallelic_polymorphic.fa",pops=["Asian","Anatolian","African"])
		#expand("data/consensus_fastas/{sample}_ref.fa",sample=config["All"]),
		#expand("data/consensus_fastas/{sample}_alt.fa",sample=config["All"])

rule split_pop:
	input:
		"/mnt/NEOGENE1/projects/donkey_2020/dataset/denovoSnpCall/merge-vcf-all_with_NS_tags_filtered.vcf.gz"
	output:
		"data/pop/{pops}_biallelic.vcf.gz"
	threads:
		1
	resources:
		mem_mb=2000
	params:
		err=lambda wildcards: "logs/pop_split/err." +wildcards.pops,
		out=lambda wildcards: "logs/pop_split/out." +wildcards.pops,
		s= lambda wildcards: ",".join(config[wildcards.pops])
	shell:
		"bcftools view -s {params.s} -v snps -Oz -o {output} {input}"

rule split_pop2:
	input:
		"/mnt/NEOGENE1/projects/donkey_2020/dataset/denovoSnpCall/merge-vcf-all_with_AN_AC_tags.vcf.gz"
	output:
		v="data/pop/{pops}_biallelic_polymorphic.vcf.gz",
		idx="data/pop/{pops}_biallelic_polymorphic.vcf.gz.csi"
	threads:
		1
	resources:
		mem_mb=2000
	params:
		err=lambda wildcards: "logs/pop_split2/err." +wildcards.pops,
		out=lambda wildcards: "logs/pop_split2/out." +wildcards.pops,
		s= lambda wildcards: ",".join(config[wildcards.pops])
	shell:
		"""
		bcftools view -s {params.s} {input} | bcftools view -v snps -m2 -M2 -i 'INFO/AC_{wildcards.pops}>0 & INFO/AN_{wildcards.pops}>INFO/AC_{wildcards.pops}' -Oz -o {output.v} 
		bcftools index {output.v}
		"""

rule consensus2:
	input:
		"data/pop/{pops}_biallelic_polymorphic.vcf.gz.csi",
		v="data/pop/{pops}_biallelic_polymorphic.vcf.gz",
		f="/mnt/NEOGENE3/share/ref/genomes/eca/Equus_caballus.EquCab2.dna_rm.toplevel.fa"
	output:
		f="data/pop/{pops}_biallelic_polymorphic_alt.fa",
		fidx="data/pop/{pops}_biallelic_polymorphic_alt.fa.fai"
	threads:
		1
	resources:
		mem_mb=1500
	params:
		err=lambda wildcards: "logs/consensus2/err." +wildcards.pops,
		out=lambda wildcards: "logs/consensus2/out." +wildcards.pops
	shell:
		"""
		bcftools consensus -f {input.f} {input.v} -o {output.f}
		samtools faidx {output.f}
		"""


rule cufflinks2:
	input:
		"data/pop/{pops}_biallelic_polymorphic_alt.fa.fai",
		f="data/pop/{pops}_biallelic_polymorphic_alt.fa",
		gff="/mnt/NEOGENE1/projects/donkey_2020/donkey_selection/ref_EquCab2.0_top_level_chr.gff3"
	output:
		"data/cufflinks_pop_fastas/{pops}_biallelic_polymorphic.fa"
	threads:
		1
	resources:
		mem_mb=3000
	params:
		err=lambda wildcards: "logs/cufflinks2/err." +wildcards.pops,
		out=lambda wildcards: "logs/cufflinks2/out." +wildcards.pops
	conda:
		"cufflinks_env.yaml"
	shell:
		"gffread -C {input.gff} -g {input.f} -x {output}"


rule split_vcf:
	input:
		"/mnt/NEOGENE1/projects/donkey_2020/dataset/denovoSnpCall/merge-vcf-all.vcf.gz"

	output:
		"data/ind/{sample}_biallelic.vcf.gz"
	threads:
		4
	resources:
		mem_mb=3000
	params:
		err=lambda wildcards: "logs/vcf_split/err." +wildcards.sample,
		out=lambda wildcards: "logs/vcf_split/out." +wildcards.sample
	shell:
		"bcftools view -m2 -M2 -s {wildcards.sample} -v snps -e 'GT=\"mis\"' -Oz -o {output} {input}"

rule vcf_index:
	input:
		"data/{sample}_biallelic.vcf.gz"
	output:
		"data/{sample}_biallelic.vcf.gz.csi"
	threads:
		1
	resources:
		mem_mb=1500
	params:
		err=lambda wildcards: "logs/index/err." +wildcards.sample,
		out=lambda wildcards: "logs/index/out." +wildcards.sample
	shell:
		"bcftools index {input}"

rule consensus:
	input:
		"data/{sample}_biallelic.vcf.gz.csi",
		v="data/{sample}_biallelic.vcf.gz",
		f="/mnt/NEOGENE3/share/ref/genomes/eca/Equus_caballus.EquCab2.dna_rm.toplevel.fa"
	output:
		r="data/consensus_fastas/{sample}_ref.fa",
		a="data/consensus_fastas/{sample}_alt.fa"
	threads:
		1
	resources:
		mem_mb=1500
	params:
		err=lambda wildcards: "logs/consensus/err." +wildcards.sample,
		out=lambda wildcards: "logs/consensus/out." +wildcards.sample
	shell:
		"""
		bcftools consensus -f {input.f} -s {wildcards.sample} -H R {input.v} -o {output.r}
		bcftools consensus -f {input.f} -s {wildcards.sample} -H A {input.v} -o {output.a}
		"""


rule samtools_index:
	input:
		"data/consensus_fastas/{sample}_{H}.fa"
	output:
		"data/consensus_fastas/{sample}_{H}.fa.fai"
	threads:
		1
	resources:
		mem_mb=1500
	params:
		err=lambda wildcards: "logs/samtools_index/err." +wildcards.sample+wildcards.H,
		out=lambda wildcards: "logs/samtools_index/out." +wildcards.sample+wildcards.H
	shell:
		"samtools faidx {input}"


rule cufflinks:
	input:
		"data/consensus_fastas/{sample}_{H}.fa.fai",
		f="data/consensus_fastas/{sample}_{H}.fa",
		gff="Equus_caballus.EquCab2.94.gff"
	output:
		"data/cufflinks_fastas/{sample}_{H}.fa"
	threads:
		1
	resources:
		mem_mb=3000
	params:
		err=lambda wildcards: "logs/cufflinks/err." +wildcards.sample+wildcards.H,
		out=lambda wildcards: "logs/cufflinks/out." +wildcards.sample+wildcards.H
	conda:
		"cufflinks_env.yaml"
	shell:
		"gffread -C {input.gff} -g {input.f} -x {output}"

use rule cufflinks as ref_cuf with:
	input:
		"/mnt/NEOGENE3/share/ref/genomes/eca/Equus_caballus.EquCab2.dna_rm.toplevel.fa.fai",
		f="/mnt/NEOGENE3/share/ref/genomes/eca/Equus_caballus.EquCab2.dna_rm.toplevel.fa",
		gff="Equus_caballus.EquCab2.94.gff"
	output:
		"data/cufflinks_fastas/EquCab2.fa"
	params:
		err=lambda wildcards: "logs/cufflinks/err.EquCab2",
		out=lambda wildcards: "logs/cufflinks/out.EquCab2"
