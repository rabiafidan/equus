configfile: "config.yaml"

localrules: all

rule all:
	input:
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



rule ref_cuf:
	input:
		"/mnt/NEOGENE3/share/ref/genomes/eca/Equus_caballus.EquCab2.dna_rm.toplevel.fa.fai",
		f="/mnt/NEOGENE3/share/ref/genomes/eca/Equus_caballus.EquCab2.dna_rm.toplevel.fa",
		gff="Equus_caballus.EquCab2.94.gff"
	output:
		"data/cufflinks_fastas/EquCab2.fa"
	threads:
		1
	resources:
		mem_mb=3000
	params:
		err=lambda wildcards: "logs/cufflinks/err.EquCab2",
		out=lambda wildcards: "logs/cufflinks/out.EquCab2"
	conda:
		"cufflinks_env.yaml"
	shell:
		"gffread -C {input.gff} -g {input.f} -x {output}"