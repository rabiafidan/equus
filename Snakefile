configfile: "config.yaml"

localrules: all


rule all:
	input:
		expand("data/consensus_fastas/{sample}_ref.fa",sample=config["Asian"]),
		expand("data/consensus_fastas/{sample}_ref.fa",sample=config["African"]),
		expand("data/consensus_fastas/{sample}_ref.fa",sample=config["Anatolian"]),
		expand("data/consensus_fastas/{sample}_alt.fa",sample=config["Asian"]),
		expand("data/consensus_fastas/{sample}_alt.fa",sample=config["African"]),
		expand("data/consensus_fastas/{sample}_alt.fa",sample=config["Anatolian"])


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
		err=lambda wildcards: "logs/samtools_index/err." +wildcards.sample,
		out=lambda wildcards: "logs/samtools_index/out." +wildcards.sample
	shell:
		"samtools faidx {input}"


rule cufflinks:
	input:
		"data/consensus_fastas/{sample}_{H}.fa.fai",
		f="data/consensus_fastas/{sample}_{H}.fa",
		gff="/mnt/NEOGENE1/projects/donkey_2020/donkey_selection/ref_EquCab2.0_top_level_chr.gff3"
	output:
		"data/cufflinks_fastas/{sample}_{H}.fa"
	threads:
		1
	resources:
		mem_mb=3000
	params:
		err=lambda wildcards: "logs/cufflinks/err." +wildcards.sample,
		out=lambda wildcards: "logs/cufflinks/out." +wildcards.sample
	conda:
		"cufflinks.yaml"
	shell:
		"gffread -C {input.gff} -g {input.f} -x {output}"