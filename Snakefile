import os

configfile: "Diss_Assembly.json"

# Shared paths
adapter_path = "misc/adapter_sequence.fa"
temp_directory = "temp/"
fastq_directory = "/scratch/general/lustre/u6035429/Diss-fastq"

# Runtime values
very_short = "6:00:00"
medium = "12:00:00"
day = "24:00:00"
long = "48:00:00"
very_long = "72:00:00"

# Parameters for splitting reference
num_chunks = 25
chunk_range = [x for x in range(1, num_chunks + 1)]
chunk_range2 = [0] + chunk_range

# Tool paths
bbduk_path = "bbduk.sh"
bcftools_path = "bcftools"
bedtools_path = "bedtools"
bgzip_path = "bgzip"
bwa_path = "bwa"
fastq_dump_path = "fastq-dump"
fastqc_path = "fastqc"
gatk_path = "gatk"
gff2bed_path = "gff2bed"
mosdepth_path = "mosdepth"
multiqc_path = "multiqc"
picard_path = "picard"
prefetch_path = "prefetch"
rename_sh_path = "rename.sh"
samtools_path = "samtools"
tabix_path = "tabix"
blast_path = "blastn"
vcflib_path = "vcflib"
vcftools_path = "vcftools"

# Samples and genomes
mapping_genomes = ["pcoq-HiC"]

sra_ids = [
	"SRR1657028", "SRR1657029", "SRR1575526", "SRR1575545", "SRR1575527",
	"SRR1575528", "SRR1575543", "SRR1575544", "SRR1575541", "SRR1575542", "SRR1575539",
	"SRR1575540", "SRR1575532", "SRR1575531", "SRR1575534", "SRR1575533", "SRR1575538",
	"SRR1575537", "SRR1575536", "SRR1575535", "SRR1575530", "SRR1575529"]

new_samples = [
	"RANO330-OMH_Pedw1", "RANO332-OMH_Pedw2", "RANO54-OMH_Pedw3", "TSINJ32-OMH_Pdia1",
	"TSINJ38-OMH_Pdia2", "TSINJ47-OMH_Pdia3", "JEJ01-OMH_Pcan1", "JEJ3-11-OMH_Pcan2",
	"MERY3-OMH_Pcan3", "ANAL10-OMH_Pper1", "TOBI5-1-OMH_Pper2", "TOBI5-3-OMH_Pper3",
	"DAR4-11-OMH_Ptat1", "DAR4-39-OMH_Ptat2", "DAR4-5-OMH_Ptat3", "JAM4-16-OMH_Pcor1",
	"JAM4-20-OMH_Pcor2", "JAM4-7-OMH_Pcor3", "KIBO15-OMH_Pdec1", "KIBO36-OMH_Pdec2",
	"KIBO44-OMH_Pdec3", "KMTEA7-10-OMH_Pver1", "KMTEA7-2-OMH_Pver2", "KMTEA7-4-OMH_Pver3",
	"DASI5-08-OMH_Alang1", "DASI5-16-OMH_Alang2", "DASI5-21-OMH_Alang3"]

initial_sample_list = new_samples + sra_ids

processed_sample_list = new_samples + [
	"Pcoq1-Marcella", "Pcoq5-Cornelia", "Pcoq6-Octavia", "Pcoq2-Trajan", "Pcoq3-Hadrian",
	"Pcoq4-Flavia", "Pdia1-Romeo", "Pdia2-Titania", "Ptat1-Agrippa", "Ptat2-Cicero",
	"Pver1-Smoke"]

QCfailed_samples = "Pcoq1-Marcella"
QCpassed_sample_list = [i for i in processed_sample_list if i != QCfailed_samples]

read1_dict = {}
read2_dict = {}

for i in new_samples:
	read1_dict[i] = os.path.join(fastq_directory, "{}_read1.fastq.gz".format(i))
	read2_dict[i] = os.path.join(fastq_directory, "{}_read2.fastq.gz".format(i))

for i in sra_ids:
	read1_dict[i] = "renamed_fastqs/{}_fixed_1.fastq.gz".format(i)
	read2_dict[i] = "renamed_fastqs/{}_fixed_2.fastq.gz".format(i)


# Filter parameters
filter_depths = ["10"]
filter_mapqs = ["40"]

rule all:
	input:
		"multiqc/multiqc_report.html",
		"multiqc_trimmed/multiqc_report.html",
		expand(
			"stats/{sample}.{genome}.sorted.mkdup.bam.stats",
			sample=processed_sample_list, genome=mapping_genomes),
		expand(
			"stats/{genome}.gatk.called.filtered_mq{mq}_dp{dp}.vcf.stats",
			genome=mapping_genomes, mq=filter_mapqs, dp=filter_depths),
		expand(
			"mosdepth_results/{sample}.{genome}.total.per-base.dp{dp}.merged.bed",
			sample=QCpassed_sample_list, genome=mapping_genomes, dp=filter_depths),

rule get_annotation:
	output:
		"reference_genomes/{genome}.gff"
	params:
		web_address = lambda wildcards: config["annotation_address"][wildcards.genome],
		initial_output = "reference_genomes/{genome}.gff.gz",
		threads = 1,
		mem = 4,
		t = day
	run:
		shell("wget {params.web_address} -O {params.initial_output}")
		shell("gunzip {params.initial_output}")

rule extract_cds_from_gff:
	input:
		"reference_genomes/{genome}.gff"
	output:
		"regions/{genome}.cds.gff"
	params:
		threads = 1,
		mem = 4,
		t = very_short
	shell:
		"""awk '($3 == "CDS")' {input} > {output}"""

rule gff2bed:
	input:
		"regions/{genome}.cds.gff"
	output:
		"regions/{genome}.cds.merged.bed"
	params:
		gff2bed = gff2bed_path,
		bedtools = bedtools_path,
		threads = 2,
		mem = 8,
		t = medium
	shell:
		"cat {input} | {params.gff2bed} | {params.bedtools} merge > {output}"

rule get_fasta:
	output:
		"reference_genomes/{genome}.fa"
	params:
		web_address = lambda wildcards: config["fasta_address"][wildcards.genome],
		initial_output = "reference_genomes/{genome}.fa.gz",
		threads = 1,
		mem = 4,
		t = very_short
	run:
		shell("wget {params.web_address} -O {params.initial_output}")
		shell("gunzip {params.initial_output}")

rule prepare_reference_fai:
	input:
		ref = "reference_genomes/{genome}.fa"
	output:
		fai = "reference_genomes/{genome}.fa.fai",
	params:
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = medium
	shell:
		"{params.samtools} faidx {input.ref}"

rule prepare_reference_dict:
	input:
		ref = "reference_genomes/{genome}.fa"
	output:
		dict = "reference_genomes/{genome}.dict"
	params:
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = medium
	shell:
		"{params.samtools} dict -o {output.dict} {input.ref}"

rule bwa_indexing:
	input:
		ref = "reference_genomes/{genome}.fa"
	output:
		amb = "reference_genomes/{genome}.fa.amb"
	params:
		bwa = bwa_path,
		threads = 4,
		mem = 16,
		t = medium
	shell:
		"{params.bwa} index {input.ref}"

rule chunk_reference:
	input:
		fai = "reference_genomes/{genome}.fa.fai"
	output:
		expand(
			"reference_genomes/{{genome}}_split_chunk{num}.bed",
			num=chunk_range)
	params:
		chunks = num_chunks,
		out_prefix = "reference_genomes/{genome}_split",
		threads = 2,
		mem = 8,
		t = medium
	shell:
		"python scripts/Chunk_fai.py --fai {input.fai} "
		"--out_prefix {params.out_prefix} --chunks {params.chunks}"

rule prefetch_sra:
	output:
		os.path.join(temp_directory, "{id}/{id}.sra")
	params:
		tool = prefetch_path,
		tmp_dir = temp_directory,
		use_id = "{id}",
		threads = 1,
		mem = 4,
		t = day
	shell:
		"{params.tool} {params.use_id} -O {params.tmp_dir} --max-size 100GB"

rule fastq_dump_paired:
	input:
		sra = os.path.join(temp_directory, "{sample}/{sample}.sra")
	output:
		fq1 = "paired_fastqs/{sample}_1.fastq.gz",
		fq2 = "paired_fastqs/{sample}_2.fastq.gz"
	params:
		output_dir = "paired_fastqs",
		fastq_dump = fastq_dump_path,
		threads = 1,
		mem = 4,
		t = day
	shell:
		"{params.fastq_dump} --outdir {params.output_dir} --gzip --readids --split-files {input.sra}"

rule fix_read_IDs_for_paired_fastqs_from_SRA_paired:
	# The fastq-dump created issues with read ID names for paired files so that
	# they give bwa issues.  This rule will go through and rename them so that
	# they're compatible with bwa
	input:
		fq1 = "paired_fastqs/{sample}_1.fastq.gz",
		fq2 = "paired_fastqs/{sample}_2.fastq.gz"
	output:
		out1 = "renamed_fastqs/{sample}_fixed_1.fastq.gz",
		out2 = "renamed_fastqs/{sample}_fixed_2.fastq.gz"
	params:
		rename_sh = rename_sh_path,
		read_name = "{sample}",
		threads = 1,
		mem = 4,
		t = day
	shell:
		"{params.rename_sh} in={input.fq1} in2={input.fq2} out={output.out1} out2={output.out2} prefix={params.read_name}"

#rule consolidate_fastqs:
#	input:
#		sra1 = expand(
#			"renamed_fastqs/{sample}_fixed_1.fastq.gz",
#			sample=sra_ids),
#		sra2 = expand(
#			"renamed_fastqs/{sample}_fixed_2.fastq.gz",
#			sample=sra_ids),
#		new1 = expand(
#			os.path.join(fastq_directory, "{sample}_read1.fastq.gz"), sample=new_samples),
#		new2 = expand(
#			os.path.join(fastq_directory, "{sample}_read2.fastq.gz"), sample=new_samples)
#	output:
#		expand(
#			"fastqs_consolidated/{sample}_{read}.fastq.gz",
#			sample=initial_sample_list,
#			read=["read1", "read2"])
#	params:
#		threads = 1,
#		mem = 4,
#		t = very_short
#	run:
#		for i in input.sra1:
#			original = i
#			basename = i.split("/")[-1].split("_")[0]
#			new_name = "fastqs_consolidated/{}_read1.fastq.gz".format(basename)
#			shell(
#				"ln -srf {original} {new_name} && touch -h {new_name}")
#		for i in input.sra2:
#			original = i
#			basename = i.split("/")[-1].split("_")[0]
#			new_name = "fastqs_consolidated/{}_read2.fastq.gz".format(basename)
#			shell(
#				"ln -srf {original} {new_name} && touch -h {new_name}")
#		for i in input.new1:
#			original = i
#			basename = i.split("/")[-1].split("_")[0]
#			new_name = "fastqs_consolidated/{}_read1.fastq.gz".format(basename)
#			shell(
#				"ln -sf {original} {new_name} && touch -h {new_name}")
#		for i in input.new2:
#			original = i
#			basename = i.split("/")[-1].split("_")[0]
#			new_name = "fastqs_consolidated/{}_read2.fastq.gz".format(basename)
#			shell(
#				"ln -sf {original} {new_name} && touch -h {new_name}")

rule fastqc_analysis_sra:
	input:
		fq1 = lambda wildcards: read1_dict[wildcards.sample],
		fq2 = lambda wildcards: read2_dict[wildcards.sample]
	output:
		fq1 = "fastqc/{sample}_fixed_1_fastqc.html",
		fq2 = "fastqc/{sample}_fixed_2_fastqc.html"
	params:
		fastqc = fastqc_path,
		threads = 1,
		mem = 4,
		t = very_short
	shell:
		"{params.fastqc} -o fastqc {input.fq1} {input.fq2}"

rule fastqc_analysis_new_samples:
	input:
		fq1 = lambda wildcards: read1_dict[wildcards.sample],
		fq2 = lambda wildcards: read2_dict[wildcards.sample]
	output:
		fq1 = "fastqc/{sample}_read1_fastqc.html",
		fq2 = "fastqc/{sample}_read2_fastqc.html"
	params:
		fastqc = fastqc_path,
		threads = 1,
		mem = 4,
		t = very_short
	shell:
		"{params.fastqc} -o fastqc {input.fq1} {input.fq2}"

rule multiqc_analysis_dna:
	input:
		sra = expand(
			"fastqc/{sample}_fixed_{read}_fastqc.html",
			sample=sra_ids,
			read=["1", "2"]),
		new = expand(
			"fastqc/{sample}_{read}_fastqc.html",
			sample=new_samples,
			read=["read1", "read2"])
	output:
		"multiqc/multiqc_report.html"
	params:
		multiqc = multiqc_path,
		threads = 1,
		mem = 4,
		t = very_short
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} --interactive -f -o multiqc fastqc"

rule trim_adapters_paired_bbduk_dna:
	input:
		fq1 = lambda wildcards: read1_dict[wildcards.sample],
		fq2 = lambda wildcards: read2_dict[wildcards.sample]
	output:
		out_fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		out_fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
	params:
		adapter = adapter_path,
		bbduk = bbduk_path,
		threads = 2,
		mem = 8,
		t = very_short
	shell:
		"{params.bbduk} -Xmx3g in1={input.fq1} in2={input.fq2} out1={output.out_fq1} out2={output.out_fq2} ref={params.adapter} ktrim=r k=21 mink=11 hdist=2 tbo tpe qtrim=rl trimq=10 trimpolyg=10"

rule fastqc_analysis_trimmed:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
	output:
		html1 = "fastqc_trimmed/{sample}_trimmed_read1_fastqc.html",
		html2 = "fastqc_trimmed/{sample}_trimmed_read2_fastqc.html"
	params:
		fastqc = fastqc_path,
		threads = 1,
		mem = 4,
		t = very_short
	shell:
		"{params.fastqc} -o fastqc_trimmed {input.fq1} {input.fq2}"

rule multiqc_analysis_trimmed_dna:
	input:
		expand(
			"fastqc_trimmed/{sample}_trimmed_{reads}_fastqc.html",
			sample=initial_sample_list, reads=["read1", "read2"])
	output:
		"multiqc_trimmed/multiqc_report.html"
	params:
		threads = 1,
		mem = 4,
		t = very_short
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"multiqc --interactive -f -o multiqc_trimmed fastqc_trimmed"

# Process BAM files

rule map_and_process_trimmed_reads:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz",
		ref = "reference_genomes/{genome}.fa",
		fai = "reference_genomes/{genome}.fa.fai",
		amb = "reference_genomes/{genome}.fa.amb",
		dict = "reference_genomes/{genome}.dict"
	output:
		"processed_bams/{sample}.{genome}.sorted.bam"
	params:
		id = lambda wildcards: config[wildcards.sample]["ID"],
		sm = lambda wildcards: config[wildcards.sample]["SM"],
		lb = lambda wildcards: config[wildcards.sample]["LB"],
		pu = lambda wildcards: config[wildcards.sample]["PU"],
		pl = lambda wildcards: config[wildcards.sample]["PL"],
		bwa = bwa_path,
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = very_long
	threads: 4
	shell:
		"{params.bwa} mem -t {params.threads} -R "
		"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
		"{input.ref} {input.fq1} {input.fq2}"
		"| {params.samtools} fixmate -O bam - - | {params.samtools} sort "
		"-O bam -o {output}"

rule index_bam:
	input:
		"processed_bams/{sample}.{genome}.sorted.bam"
	output:
		"processed_bams/{sample}.{genome}.sorted.bam.bai"
	params:
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = very_short
	shell:
		"{params.samtools} index {input}"

rule merge_bams:
	input:
		bams = lambda wildcards: expand(
			"processed_bams/{sample}.{genome}.sorted.bam",
			sample=config["to_merge"][wildcards.sample], genome=wildcards.genome),
		bais = lambda wildcards: expand(
			"processed_bams/{sample}.{genome}.sorted.bam.bai",
			sample=config["to_merge"][wildcards.sample], genome=wildcards.genome)
	output:
		"processed_bams/{sample}.{genome}.sorted.merged.bam"
	params:
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = long
	shell:
		"""
		input_bams=({input.bams})
		merged_output={output}
		if [ "${{#input_bams[@]}}" -gt "1" ]; then
			{params.samtools} merge {output} {input.bams}
		else
			ln -s ../${{input_bams[0]}} $merged_output && touch -h $merged_output
		fi
		"""

rule index_merged_bam:
	input:
		"processed_bams/{sample}.{genome}.sorted.merged.bam"
	output:
		"processed_bams/{sample}.{genome}.sorted.merged.bam.bai"
	params:
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = very_short
	shell:
		"{params.samtools} index {input}"

rule picard_mkdups:
	input:
		bam = "processed_bams/{sample}.{genome}.sorted.merged.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.merged.bam.bai"
	output:
		bam = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam",
		metrics = "metrics/{sample}.{genome}.picard_mkdup_metrics.txt"
	params:
		picard = picard_path,
		temp_dir = temp_directory,
		threads = 8,
		mem = 24,
		t = long,
	shell:
		"{params.picard} -Xmx14g -Djava.io.tmpdir={params.temp_dir} MarkDuplicates I={input.bam} O={output.bam} "
		"M={output.metrics}"

rule index_mkdup_bam:
	input:
		"processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam"
	output:
		"processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam.bai"
	params:
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = very_short
	shell:
		"{params.samtools} index {input}"

rule bam_stats:
	input:
		bam = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam.bai"
	output:
		"stats/{sample}.{genome}.sorted.mkdup.bam.stats"
	params:
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = very_short
	shell:
		"{params.samtools} stats {input.bam} | grep ^SN | cut -f 2- > {output}"

# Call and filter variants
rule gatk_gvcf_per_chunk:
	input:
		ref = "reference_genomes/{genome}.fa",
		dict = "reference_genomes/{genome}.dict",
		bam = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam",
		bai = "processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam.bai",
		chunkfile = "reference_genomes/{genome}_split_chunk{chunk}.bed"
	output:
		"gvcf/{sample}.{genome}.{chunk}.g.vcf.gz"
	benchmark:
		"benchmark/{sample}.{genome}.{chunk}.HaplotypeCaller.benchmark.txt"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path,
		threads = 8,
		mem = 32,
		t = very_long
	shell:
		"""{params.gatk} --java-options "-Xmx30g -Djava.io.tmpdir={params.temp_dir}" """
		"""HaplotypeCaller -R {input.ref} -I {input.bam} -L {input.chunkfile} """
		"""-ERC GVCF --do-not-run-physical-phasing -O {output} """

rule genomicsdbimport_combine_gvcfs_per_chunk:
	input:
		ref = "reference_genomes/{genome}.fa",
		gvcfs = lambda wildcards: expand(
			"gvcf/{sample}.{genome}.{chunk}.g.vcf.gz",
			sample=QCpassed_sample_list,
			genome=[wildcards.genome],
			chunk=[wildcards.chunk]),
		chunkfile = "reference_genomes/{genome}_split_chunk{chunk}.bed"
	output:
		directory("gvcf_databases/{genome}-{chunk}")
	benchmark:
		"benchmark/{genome}.{chunk}.GenomicsDBImport.benchmark.txt"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path,
		threads = 8,
		mem = 86,
		t = very_long
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx84g -Djava.io.tmpdir={params.temp_dir}" """
			"""GenomicsDBImport -R {input.ref} {variant_files} """
			"""--genomicsdb-workspace-path {output} -L {input.chunkfile} --genomicsdb-shared-posixfs-optimizations true""")

rule gatk_genotypegvcf_genomicsdb:
	input:
		gvcf = "gvcf_databases/{genome}-{chunk}",
		ref = "reference_genomes/{genome}.fa"
	output:
		"vcf_genotyped/{genome}.{chunk}.gatk.called.raw.vcf.gz"
	benchmark:
		"benchmark/{genome}.{chunk}.GenotypeGVCFs.benchmark.txt"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path,
		threads = 8,
		mem = 84,
		t = very_long
	shell:
		"""{params.gatk} --java-options "-Xmx82g -Djava.io.tmpdir={params.temp_dir}" """
		"""GenotypeGVCFs -R {input.ref} -V gendb://{input.gvcf} -O {output}"""

rule concatenate_split_vcfs:
	input:
		vcf = lambda wildcards: expand(
			"vcf_genotyped/{genome}.{chunk}.gatk.called.raw.vcf.gz",
			genome=wildcards.genome,
			chunk=chunk_range2)
	output:
		"combined_vcfs/combined.{genome}.raw.vcf.gz"
	benchmark:
		"benchmark/{genome}.concatVCF.benchmark.txt"
	params:
		bcftools = bcftools_path,
		threads = 2,
		mem = 8,
		t = medium
	shell:
		"{params.bcftools} concat -O z -o {output} {input.vcf}"

rule index_concatenated_vcf:
	input:
		"combined_vcfs/combined.{genome}.raw.vcf.gz"
	output:
		"combined_vcfs/combined.{genome}.raw.vcf.gz.tbi"
	benchmark:
		"benchmark/{genome}.IndexConcatVCF.benchmark.txt"
	params:
		tabix = tabix_path,
		threads = 4,
		mem = 16,
		t = very_short
	shell:
		"{params.tabix} -p vcf {input}"


rule filter_vcfs:
	input:
		vcf = "combined_vcfs/combined.{genome}.raw.vcf.gz",
		idx = "combined_vcfs/combined.{genome}.raw.vcf.gz.tbi"
	output:
		"vcf_filtered/{genome}.gatk.called.filtered_mq{mq}_dp{dp}.vcf.gz"
	benchmark:
		"benchmark/{genome}.filterVCF_mq{mq}_dp{dp}.benchmark.txt"
	params:
		bgzip = bgzip_path,
		bcftools = bcftools_path,
		mapq = "{mq}",
		dp = "{dp}",
		threads = 4,
		mem = 16,
		t = long
	shell:
		"{params.bcftools} filter -i "
		"'QUAL >= 30 && MQ >= {params.mapq} && QD > 2' {input} | "
		"{params.bcftools} filter -i 'FMT/DP >= {params.dp} & FMT/GQ >= 30' -S . - | "
		"{params.bgzip} > {output}"

rule index_filtered_vcf:
	input:
		"vcf_filtered/{genome}.gatk.called.filtered_mq{mq}_dp{dp}.vcf.gz"
	output:
		"vcf_filtered/{genome}.gatk.called.filtered_mq{mq}_dp{dp}.vcf.gz.tbi"
	benchmark:
		"benchmark/{genome}.IndexFilterVCF_mq{mq}_dp{dp}.benchmark.txt"
	params:
		tabix = tabix_path,
		threads = 4,
		mem = 16,
		t = very_short
	shell:
		"{params.tabix} -p vcf {input}"

rule vcf_stats:
	input:
		vcf = "vcf_filtered/{genome}.gatk.called.filtered_mq{mq}_dp{dp}.vcf.gz",
		tbi = "vcf_filtered/{genome}.gatk.called.filtered_mq{mq}_dp{dp}.vcf.gz.tbi"
	output:
		"stats/{genome}.gatk.called.filtered_mq{mq}_dp{dp}.vcf.stats"
	benchmark:
		"benchmark/{genome}.StatsFilteredVCF_mq{mq}_dp{dp}.benchmark.txt"
	params:
		bcftools = bcftools_path,
		threads = 4,
		mem = 16,
		t = very_short
	shell:
		"{params.bcftools} stats {input.vcf} | grep ^SN > {output}"

rule mosdepth_depth:
	input:
		"processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam"
	output:
		dist = "mosdepth_results/{sample}.{genome}.total.mosdepth.global.dist.txt",
		per_base = "mosdepth_results/{sample}.{genome}.total.per-base.bed.gz"
	params:
		mosdepth = mosdepth_path,
		prefix = "mosdepth_results/{sample}.{genome}.total",
		threads = 4,
		mem = 16,
		t = long
	shell:
		"{params.mosdepth} --fast-mode -F 1024 {params.prefix} {input}"

rule filter_mosdepth:
	input:
		"mosdepth_results/{sample}.{genome}.total.per-base.bed.gz"
	output:
		"mosdepth_results/{sample}.{genome}.total.per-base.dp{dp}.merged.bed"
	params:
		dp = "{dp}",
		bedtools = bedtools_path,
		threads = 4,
		mem = 16,
		t = very_short
	shell:
		"zcat {input} | awk '$4 >= {params.dp}' | {params.bedtools} merge -i - > {output}"
