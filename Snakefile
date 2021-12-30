configfile: "config.yaml"

localrules: all


rule all:
	input:
		expand("data/cufflinks_fastas/{sample}_{H}.fa",sample=config["All"],H=["ref","alt"]),
		"data/cufflinks_fastas/EquCab2.fa"
		#expand("data/consensus_fastas/{sample}_ref.fa",sample=config["All"]),
		#expand("data/consensus_fastas/{sample}_alt.fa",sample=config["All"])


rule split_vcf:
	input:
		"/mnt/NEOGENE1/projects/donkey_2020/dataset/denovoSnpCall/merge-vcf-all.vcf.gz"

	output:
		"data/{sample}_biallelic.vcf.gz"
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
