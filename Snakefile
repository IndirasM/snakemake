configfile: "config.yaml"

rule all:
    input:
        expand(config["out_fastqc"] + "{sample}_R{read}_001_fastqc.html", sample=config["samples"], read=config["reads"]),
        config["out_multiqc"] + "multiqc_report.html",
        expand(config["data_trimmed"] + "{sample}_R1_001.fastq", sample=config['samples']),
        expand(config["data_trimmed"] + "{sample}_R2_001.fastq", sample=config['samples']),
        expand(config["out_fastqc_trimmed"] + "{sample}_R{read}_001_fastqc.html", sample=config["samples"], read=config["reads"]),
        config["out_multiqc_trimmed"] + "multiqc_report.html"

rule run_fastqc:
    input:
        expand(config["data"] + "{sample}_R{read}_001.fastq.gz", sample=config["samples"], read=config["reads"])
    output:
        expand(config["out_fastqc"] + "{sample}_R{read}_001_fastqc.html", sample=config["samples"], read=config["reads"])
    shell:
        "fastqc {input} --outdir=" + config['out_fastqc']

rule run_multiqc:
    input:
        expand(config["out_fastqc"] + "{sample}_R{read}_001_fastqc.html", sample=config['samples'], read=config["reads"])
    output:
        config["out_multiqc"] + "multiqc_report.html"
    shell:
        "multiqc " + config["out_fastqc"] + " -o " + config["out_multiqc"]

rule run_bbduk_trim:
    input:
        read1 = config["data"] + "{sample}_R1_001.fastq.gz",
        read2 = config["data"] + "{sample}_R2_001.fastq.gz"
    output:
        output1 = config["data_trimmed"] + "{sample}_R1_001.fastq",
        output2 = config["data_trimmed"] + "{sample}_R2_001.fastq"
    shell:
        "scripts/bbmap/bbduk.sh in1={input.read1} out1={output.output1} in2={input.read2} out2={output.output2} ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10"

rule run_fastqc_trimmed:
    input:
        expand(config["data_trimmed"] + "{sample}_R{read}_001.fastq", sample=config["samples"], read=config["reads"])
    output:
        expand(config["out_fastqc_trimmed"] + "{sample}_R{read}_001_fastqc.html", sample=config["samples"], read=config["reads"])
    shell:
        "fastqc {input} --outdir=" + config['out_fastqc_trimmed']

rule run_multiqc_trimmed:
    input:
        expand(config["out_fastqc_trimmed"] + "{sample}_R{read}_001_fastqc.html", sample=config['samples'], read=config["reads"])
    output:
        config["out_multiqc_trimmed"] + "multiqc_report.html"
    shell:
        "multiqc " + config["out_fastqc_trimmed"] + " -o " + config["out_multiqc_trimmed"]
# rule start:
#     shell:
#         "curl -L https://github.com/snakemake/snakemake-tutorial-data/archive/v5.4.5.tar.gz -o v5.4.5.tar.gz && tar -xf v5.4.5.tar.gz --strip 1 '*/data'"
#
# rule bwa_map:
#     input:
#         "data/genome.fa",
#         "data/samples/{sample}.fastq"
#     output:
#         "mapped_reads/{sample}.bam"
#     shell:
#         "bwa mem {input} | samtools view -Sb - > {output}"
#
# rule samtools_sort:
#     input:
#         "mapped_reads/{sample}.bam"
#     output:
#         "sorted_reads/{sample}.bam"
#     shell:
#         "samtools sort -T sorted_reads/{wildcards.sample} "
#         "-O bam {input} > {output}"
#
# rule samtools_index:
#     input:
#         "sorted_reads/{sample}.bam"
#     output:
#         "sorted_reads/{sample}.bam.bai"
#     shell:
#         "samtools index {input}"
#
# rule bcftools_call:
#     input:
#         fa="data/genome.fa",
#         bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
#         bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
#     output:
#         "calls/all.vcf"
#     shell:
#         "samtools mpileup -g -f {input.fa} {input.bam} | "
#         "bcftools call -mv - > {output}"
#
# rule plot_quals:
#     input:
#         "calls/all.vcf"
#     output:
#         "plots/quals.svg"
#     script:
#         "scripts/plot-quals.py"
