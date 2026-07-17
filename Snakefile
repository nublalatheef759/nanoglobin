import glob
import os

configfile: "config.yml"

if config["mode"] == "amplicon":
    REFERENCE = config["ref_amplicon"]
else:
    REFERENCE = config["ref_wgs"]

SAMPLES = config["samples"]

rule all:
    input:
        "results/clinical_reports/sv_annotated.csv"

rule minimap_align:
    input:
        ref=REFERENCE,
        fastq="fastq/{sample}.fastq"
    output:
        temp("aligned/{sample}.bam")
    log:
        "logs/minimap/{sample}.log"
    threads: 8
    shell:
        "(minimap2 -a -x map-ont -t {threads} {input.ref} {input.fastq} | samtools view -Sb - > {output}) 2> {log}"

rule samtools_sort:
    input:
        "aligned/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} -O bam {input} > {output}"

rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"

rule copy_number:
    input:
        bam="sorted_reads/{sample}.bam",
        bai="sorted_reads/{sample}.bam.bai"
    output:
        "variants/{sample}/{sample}.coverage.tsv"
    shell:
        "samtools depth -r chr16:170000-178000 {input.bam} | "
        "awk '{{sum+=$3; n++}} END {{print \"HBA_mean_depth\\t\"sum/n}}' > {output} && "
        "samtools depth -r chr11:5225000-5228000 {input.bam} | "
        "awk '{{sum+=$3; n++}} END {{print \"HBB_mean_depth\\t\"sum/n}}' >> {output}"

rule sniffles_sv:
    input:
        bam="sorted_reads/{sample}.bam",
        bai="sorted_reads/{sample}.bam.bai",
        ref="reference/hg38_globin.fa"
    output:
        "variants/{sample}/{sample}.sniffles.vcf"
    threads: 4
    log:
        "logs/sniffles/{sample}.log"
    shell:
        "{config[sniffles_path]} --input {input.bam} --reference {input.ref} "
        "--vcf {output} --threads {threads} 2> {log}"

rule cutesv_sv:
    input:
        bam="sorted_reads/{sample}.bam",
        bai="sorted_reads/{sample}.bam.bai",
        ref="reference/hg38_globin.fa"
    output:
        "variants/{sample}/{sample}.cutesv.vcf"
    threads: 4
    log:
        "logs/cutesv/{sample}.log"
    shell:
        "rm -rf variants/{wildcards.sample}/cutesv_tmp && "
        "mkdir -p variants/{wildcards.sample}/cutesv_tmp && "
        "{config[cutesv_path]} {input.bam} {input.ref} {output} "
        "variants/{wildcards.sample}/cutesv_tmp "
        "--threads {threads} --sample {wildcards.sample} "
        "--min_support {config[cutesv_min_support]} --min_size 50 --genotype 2> {log}"
        
rule clair3_call:
    input:
        bam="sorted_reads/{sample}.bam",
        bai="sorted_reads/{sample}.bam.bai",
        ref="reference/hg38_globin.fa"
    output: "variants/{sample}/phased_merge_output.vcf.gz"
    threads: 4
    log:
        "logs/clair3/{sample}.log"
        
    shell:
        "{config[clair3_path]} "
        "--bam_fn={input.bam} --ref_fn={input.ref} "
        "--output=variants/{wildcards.sample} --threads={threads} "
        "--platform=ont --model_path={config[clair3_models]} --enable_phasing 2> {log}"

rule filter_variants:
    input:
        "variants/{sample}/phased_merge_output.vcf.gz"
    output:
        "results/{sample}.filtered.vcf"
    shell:
        "bcftools filter -i 'QUAL>{config[min_quality]} && DP>{config[min_depth]}' {input} > {output}"

rule summary_table:
    input:
        expand("results/{sample}.filtered.vcf", sample=SAMPLES)
    output:
        "results/variant_summary.csv"
    shell:
        """
        echo "Sample,Chromosome,Position,Ref,Alt,Quality,Genotype,Depth,AlleleFreq" > {output}
        for vcf in {input}; do
            sample=$(basename $vcf .filtered.vcf)
            grep -v "^#" $vcf | awk -v s="$sample" '{{
                split($10, gt, ":");
                print s","$1","$2","$4","$5","$6","gt[1]","gt[4]","gt[6]
            }}' >> {output}
        done
        """
rule region_filter:
    input:
        "results/{sample}.filtered.vcf"
    output:
        "results/{sample}.globin_only.vcf"
    shell:
        """
        awk '/^#/ {{print; next}} \
             ($1=="chr11" && $2>=5225000 && $2<=5228000) || \
             ($1=="chr16" && $2>=170000 && $2<=178000)' {input} > {output}
        """
        
rule annotate_variants:
    input:
        "results/{sample}.globin_only.vcf"
    output:
        "results/{sample}.annotated.csv"
    log:
        "logs/annotate/{sample}.log"
    shell:
        "python scripts/vep_annotate.py {input} {output} 2> {log}"

rule merge_annotations:
    input:
        expand("results/{sample}.annotated.csv", sample=SAMPLES)
    output:
        "results/all_variants_annotated.csv"
    shell:
        """
        head -1 {input[0]} > {output}
        for f in {input}; do
            tail -n +2 $f >> {output}
        done
        """
rule comprehensive_summary:
    input:
        annotated="results/all_variants_annotated.csv",
        sniffles=expand("variants/{sample}/{sample}.sniffles.vcf", sample=SAMPLES),
        cutesv=expand("variants/{sample}/{sample}.cutesv.vcf", sample=SAMPLES),
        coverage=expand("variants/{sample}/{sample}.coverage.tsv", sample=SAMPLES)
    output:
        "results/comprehensive_report.csv"
    shell:
        "python scripts/comprehensive_report.py {output}"

rule patient_summary:
    input:
        "results/comprehensive_report.csv"
    output:
        "results/patient_summary.csv"
    shell:
        "python scripts/patient_summary.py {input} {output}"
        
rule clinical_annotation:
    input: 
        annotated="results/all_variants_annotated.csv",
        patients="results/patient_summary.csv", 
    output:
        "results/clinical_reports/sv_annotated.csv"
    shell:
        "python scripts/annotation/annotate_variants.py"
        
    

    